#
# Some functions
#

#
# The following function computes the tal-dependence coefficient between two vectors
#
library(Rcpp)
tdep <- function(x,y,p,se_return=F){
  L = c();
  for (pi in p){
    if (pi==1){
      L = c(L, 0)}
    else{
    qx = quantile(x,pi,na.rm=TRUE);
    qy = quantile(y,pi,na.rm=TRUE);
    L = c(L, sum((x>qx)&(y>qy))/(length(x)*(1-pi)))
    }
  }
  n = length(x);
  if (se_return){
    return(list("L"=L,"se"=sqrt(L*(1-L*(1-p))/(n*(1-p)) )))
  } else{
  return(L)
  }
}

#
# The foloowing function computes the *partial* tail dependence coefficient, i.e.,
#
# P[ F_X(X)>p, F_Y(p)>p, F_{Z_i}(Z_i)<=p, i=1,..,k]/(1-p)
#

tdep_partial <- function(x,y,z,p){
  L = c();
  n = length(unlist(x))
  if(length(z)==0){
    return(tdep(x,y,p))
  }
  if (length(dim(z))==0){
    z = matrix(z,n,1)
  }
  uz = apply(z,2,rank)/(n+1);
  uz = apply(uz,1,max);
  for (pi in p){
    if (pi==1){
      L = c(L, 0)
    } else{
    qx = quantile(x,pi,na.rm=TRUE);
    qy = quantile(y,pi,na.rm=TRUE);
    L = c(L, sum((x>qx)&(y>qy)&(uz<=pi))/(n*(1-pi)))
    }
  }
  return(L)
}

tdep_partial_select <- function(y,x,p){
#
# Performs sequential variable selection by maximizing the partial tail 
# dependence coefficients
#
if (length(dim(x))==0){
  errorCondition("Nothing to select from.")
} else{
  k = dim(x)[2];
  keep =c();
  for (j in c(1:k)){
    Lmax = 0;
    for (i in setdiff(c(1:k),keep)){
      Lnew = tdep_partial(y,x[,i],z=x[,keep],p);
      if(Lmax<Lnew){imax = i; Lmax = Lnew} 
    }
    if (Lmax==0){ return(keep) } else {
      keep = c(keep, imax);
    }
  }
  return(keep)
}}

#
# The following function computes efficiently block-maxima (both disjoint and sliding)
#
bmax <- cppFunction('NumericVector get_b(NumericVector a, long m, LogicalVector sliding_blocks) {
  long n = a.size();
  long nm  = floor(n/m);
  long idx;
  if (sliding_blocks[0]==true) {
    nm = n-m+1;
  } 
  NumericVector b(nm);
  for (long j=0; j<nm; j++){
   if (sliding_blocks[0]==true) 
     b(j) = a(j);
   else 
     b(j) = a(j*m);
   for (int i=0; i<m; i++){
     idx = j*m+i;
     if (sliding_blocks[0]==true){
       idx = j+i;
     }
     if (a(idx)>b(j)) 
      b(j) = a(idx);
     }
   }
   return b;
   }')

bmaxCol = function(M,block=1,to_slide=F){
  return(apply(M,2,bmax,block,to_slide))
}
  
block_tdep <- function(x,y,p,m,sliding_blocks=FALSE, to_plot=FALSE){
  #
  # Each row consists of tail-dependence coefficients computed from the non-overlapping 
  # block-maxima corresponding to the block-sizes in the vector "m".  
  #
  # NOTE: consider adding the option of sliding block-maxima
  #
  L = c()
  for (i in c(1:length(m))){
    L = rbind(L,tdep(bmax(x,m[i],sliding_blocks),bmax(y,m[i],sliding_blocks),p)); 
  }
  rownames(L) = as.character(m)
  colnames(L) = as.character(p)
  
  if (to_plot){
    library(fields)
    image.plot(x=m,y=p,(L),
               ylab="Prob",
               xlab="Block size",
               main=paste0("Tail Dependence: Sliding Blocks=",sliding_blocks))
  }
  return(L)
}

#
# Tail dependence  
#

max_tdep_lag <- function(x,y,p,m=1,sliding_blocks=TRUE, max_not_sums= TRUE, lag_windows = c(-10:10),to_plot=TRUE){
 if (max_not_sums==TRUE){
   x = bmax(x,m,sliding_blocks = sliding_blocks)
   y = bmax(y,m,sliding_blocks = sliding_blocks)
 } else {
   x = cumsum(x);
   y = cumsum(y);
   nx = length(x);
   if (sliding_blocks==TRUE){
     x= x[(m+1):nx] - x[1:(nx-m)];
     y= y[(m+1):nx] - y[1:(nx-m)];
   } else {
     x = diff(x[m*c(1:floor(nx/m))])
     y = diff(y[m*c(1:floor(nx/m))])
   }
 }
 n = length(x);
  L = c(); 
  for (lag in lag_windows){
    idx_x = c(1:n)+lag;
    idx_x = idx_x[(idx_x>0)&(idx_x<=n)];
    if (lag>0){
     idx_y = c(1:length(idx_x));
    } else {
      idx_y =abs(lag)+c(1:length(idx_x));
    }
    L = rbind(L, tdep(x[idx_x],y[idx_y],p))
  }
  rownames(L) = lag_windows;
  colnames(L) = p;
  
  if (to_plot){
    library(fields)
    image.plot(x=lag_windows,y=p,L,
               ylab="Prob",
               xlab="Lag",
               main=paste0("Tail Dependence Lag"))
  }
  return(L)
}

granger_tdep_graph <- function(data, p, block_size = 1,max_not_sums= TRUE,sliding_blocks=TRUE,lag_windows = c(-10:10),to_plot=TRUE){
  n_data = dim(data)[2];
  names = colnames(data);
  if (length(names)<n_data){
    names = c(1:n_data)
    colnames(data) = names
  }
  A = matrix(0,n_data,n_data);
  L = matrix(0,n_data,n_data);
  diag(L) = 1;
  Bdata = c();
 for (i in c(1:n_data)){
   Bdata = cbind( Bdata, bmax(data[,i],m=block_size,sliding_blocks = sliding_blocks))
 }
  colnames(Bdata) = names
 for (i in c(1:n_data)){
   cat("\r",i," out of ",n_data)
   for (j in c(1:n_data)){
     if (i!=j){
     Lij = max_tdep_lag(data[,i],data[,j],p,lag_windows = lag_windows,max_not_sums = max_not_sums, to_plot = F);
     idx = which.max(Lij);
     A[i,j] = lag_windows[idx];
     L[i,j] = Lij[idx];
     }
   }
 } 
  colnames(A) = names;
  rownames(A) = names;
  colnames(L) = names;
  rownames(L) = names;
  if (to_plot==TRUE){
    library(fields)
    diag(L) = NA
    image.plot(A,main="Adjacency Matrix")
    image.plot(L,main="Max Tail-dependence")
  }
  return(list("max_tdep"=L,"max_tdep_lag"=A))
}

plot_the_Granger_tail_dependence_graph <- function(A,L){
  # Install and load the igraph package if not already installed
  library(igraph)
  
   # Convert adjacency matrix to an igraph object
  g <- graph_from_adjacency_matrix(A,
                                   mode = "directed",  # Directed graph (oriented)
                                   weighted = TRUE)    # Use the matrix values as edge weights
  #E(g)$weight <- L
  # Plot the graph
  LO = layout_with_fr(g)
  
  Isolated = which(degree(g)==0)
  #cat(degree(g))
  g2 = delete_vertices(g, Isolated)
  LO2 = LO[-Isolated,]
  plot(g2, layout=LO2)
  
  plot(g2, 
       edge.label = paste0(E(g)$weight,L),  # Display edge weights as labels
       edge.width = E(g)$weight,     # Scale edge width by weight
       vertex.size = 10,                 # Vertex size
       vertex.label.cex = 0.7,           # Label size for vertices
       edge.arrow.size = 0.2,            # Arrow size for edge direction
       edge.label.cex=0.6,
       layout = layout_with_fr,          # Layout for positioning vertices
       main = "Weighted Directed Graph from Adjacency Matrix")
}
