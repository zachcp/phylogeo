#' @importFrom(sp, spDists)
#' @importFrom(reshape2, melt)
#' @importFrom(phyloseq, distances)
plot_greatcircle_distance <- function(physeq, distancemethod="jaccard"){
    latlon <- .check_physeq(physeq)
    latcol <- as.character( latlon[1] )
    loncol <- as.character( latlon[2] )
    data   <- data.frame( sample_data(physeq) )
    data   <- .check_NA(data, latcol)
    data   <- .coerce_numeric(data,latcol)
    data   <- .check_NA(data, loncol)
    data   <- .coerce_numeric(data,loncol)
    names  <- names(data)
    
    #get bigcircle distances using spDists
    df2 <- data[ c(loncol, latcol) ]
    names(df2) <- c("lon", "lat")
    df2$lat <- sapply(df2$lat, .degree_to_radian)
    df2$lon <- sapply(df2$lon, .degree_to_radian)
    df2 <- as.matrix(df2)
    geodistances <- spDists(df2, longlat=TRUE)
    colnames(geodistances) <- row.names(df2)
    row.names(geodistances) <- row.names(df2)
    geodistances <- melt(geodistances)
    names(geodistances) <- c("Var1", "Var2", "geodist")

    #get ecologicaldistances
    ecodistance <- distance(physeq, method = distancemethod)
    ecodistance <- as.matrix(ecodistance)
    ecodistance <- melt(ecodistance)
    names(ecodistance) <- c("Var1", "Var2", "ecodist")
    
    #make mergeable names for the two distance functions and merge
    concatvals <- function(x,y){ return(paste(x,"_",y,sep=""))}
    geodistances['pairs'] <- mapply(concatvals, geodistances$Var1, geodistances$Var2)
    ecodistance['pairs'] <- mapply(concatvals, ecodistance$Var1, ecodistance$Var2)
    df <- merge(geodistances, ecodistance, by="pairs")
    
    #make the plot
    p <- ggplot(df, aes(y = ecodist,x=geodist)) + geom_point() 
    return(p)
}
 

.great_circle_distance <- function(row1,row2) {
  
  ### get great circle distance from two lat/lon coordinates
  ### data format is in a row where the first column is lat an dthe second is lon
  lat1 <- as.numeric( row1[1] )
  lon1 <- as.numeric( row1[2] )
  lat2 <- as.numeric( row2[1] )
  lon2 <- as.numeric( row2[2] )
  
  print(as.numeric(lat1))
  #check data is small enough to be radians
  max_radian <- pi * 2 
  message <- "lat and Long must be in radians"
  for( x in c(lat1,lon1,lat2,lon2)){
    if(x >max_radian){
      stop(message)
    }
  }
  
  ## http://www.movable-type.co.uk/scripts/latlong.html
  # Haversine formula:  
  # a = sin²(Δφ/2) + cos φ1 ⋅ cos φ2 ⋅ sin²(Δλ/2)
  # c = 2 ⋅ atan2( √a, √(1−a) )
  # d = R ⋅ c
  # where  φ is latitude, λ is longitude, R is earth’s radius (mean radius = 6,371km);
  # note that angles need to be in radians to pass to trig functions!
  #   JavaScript:  
  #   var R = 6371; // km
  # var φ1 = lat1.toRadians();
  # var φ2 = lat2.toRadians();
  # var Δφ = (lat2-lat1).toRadians();
  # var Δλ = (lon2-lon1).toRadians();
  # 
  # var a = Math.sin(Δφ/2) * Math.sin(Δφ/2) +
  #   Math.cos(φ1) * Math.cos(φ2) *
  #   Math.sin(Δλ/2) * Math.sin(Δλ/2);
  # var c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1-a));
  # 
  # var d = R * c;
  #Note in these scripts, I generally use lat/lon for latitude/longitude in degrees, and φ/λ for latitude/longitude in radians – having found that mixing degrees & radians is often the easiest route to head-scratching bugs...
  
  r= 6378137 #meters
  delta_lon = abs(lon1 - lon2)
  delta_lat = abs(lat1 - lat2)
  print(delta_lon)
  a <- sin(delta_lat/2) * sin(delta_lat/2) + cos(lon1) * cos(lon2) * sin( delta_lon/2) * sin(delta_lon/2)
  c <- 2 * atan2( sqrt(a), sqrt(1-a))
  dist <- r*c
  return(dist)
}
  
.degree_to_radian <- function(degree) {
  ### angle in radians = angle in degrees * Pi / 180
  radian <- degree * pi / 180
  return( radian )
}
.radian_to_degree <- function(radian) {
  ### angle in radians * 180 / Pi = angle in degrees
  degree <- radian * 180 / pi
  return(degree)
}