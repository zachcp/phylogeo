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
#' @import phyloseq
#' @importFrom reshape2 melt
#' @export
#' @example
#' data(batfecal)
#' plot_greatcircle_distance(batfecal)
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
    geodistances <- .dist_to_edge_table(geodistances, dname = "geodist")
    
    #geodistances <- melt(geodistances)
    #names(geodistances) <- c("Var1", "Var2", "geodist")
    
    #get ecologicaldistances
    ecodistance <- distance(physeq, method = distancemethod)
    ecodistance <- .dist_to_edge_table(ecodistance, dname="ecodist" )
    #ecodistance <- as.matrix(ecodistance)
    #ecodistance <- melt(ecodistance)
    #names(ecodistance) <- c("Var1", "Var2", "ecodist")
        
    
    #make mergeable names for the two distance functions and merge
    concatvals <- function(x,y){ return(paste(x,"_",y,sep=""))}
    geodistances['pairs'] <- mapply(concatvals, geodistances$Var1, geodistances$Var2)
    ecodistance['pairs'] <- mapply(concatvals, ecodistance$Var1, ecodistance$Var2)
    df <- merge(geodistances, ecodistance, by="pairs")
    
    #make the plot
    p <- ggplot(df, aes(y = ecodist,x=geodist)) + geom_point() 
    return(p)
}
#' Utility Function for Converting Distance Matrices to 
#' three column distances while removing all of the duplicates
#' lifted/modified from here: https://github.com/joey711/phyloseq/blob/master/R/plot-methods.R
#' 
#' @importFrom reshape2 melt
.dist_to_edge_table = function(Dist, dname = "dist"){
  dmat <- as.matrix(Dist)
  # Set duplicate entries and self-links to Inf
  dmat[upper.tri(dmat, diag = TRUE)] <- Inf
  df_3col = reshape2::melt(dmat, as.is = TRUE)
  # Eliminate Inf Values (melt's third column is "value")
  df_3col <- df_3col[is.finite(df_3col$value), ]
  #change names
  names(df_3col) <- c("Var1","Var2", dname)
  return(df_3col)
}
#' Utility Function for Converting Degrees to Radians
.degree_to_radian <- function(degree) {
  ### angle in radians = angle in degrees * Pi / 180
  radian <- degree * pi / 180
  return( radian )
}
#' Utility Function for Converting Radians to Degrees
.radian_to_degree <- function(radian) {
  ### angle in radians * 180 / Pi = angle in degrees
  degree <- radian * 180 / pi
  return(degree)
}