library(arrow)
library(MASS)
library(fields)

da = read_parquet("~/Dropbox/doc/data/harp_xray_flux_data/high-qual_near-center-70_no-nas_flares.parquet")
idx = which(is.na(da$flare_class)==FALSE);
da0 = da[idx,];

idx_M = which(da0$flare_class=="M");
idx_X = which(da0$flare_class=="X");
idx_MX = union(idx_M,idx_X);

peak_int_tps = Tps(x=cbind(da0$LON_FWT[idx_MX],da0$LAT_FWT[idx_MX]),Y=log10(da0$peak_intensity[idx_MX]),lambda=0);
peak_int_surf = predictSurface(peak_int_tps);
image.plot(peak_int_surf,xlab="LON",ylab="LAT")

latlon_to_unit_vector <- function(lat, lon) {
  # Convert degrees to radians
  lat_rad <- lat * pi / 180
  lon_rad <- lon * pi / 180
  
  # Compute unit vector components
  x <- cos(lat_rad) * cos(lon_rad)
  y <- cos(lat_rad) * sin(lon_rad)
  z <- sin(lat_rad)
  
  # Return as data frame or list
  return(data.frame(x = x, y = y, z = z))
}
