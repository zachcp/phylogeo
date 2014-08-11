#' Calculate and Plot Sample Distances by Geogrpahy/Ecological Distance
#'
#' Scatterplot generation of samples using geogrpahic and eoclogical distances
#'
#' @usage plot_greatcircle_distance(physeq, distancemethod="jaccard")
#'
#' @param physeq (Required). 
#'  The name of the phyloseq object. This must have sample data with Latitude and Longitude Columns.
#'  
#' @param distancemethod (Optional). Default \code{"jaccard"}.
#'  The name of an ecological distance method. See ?distance for more information
#'  
#' @importFrom sp spDists
#' @importFrom reshape2 melt
#' @import phyloseq
#' @export
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