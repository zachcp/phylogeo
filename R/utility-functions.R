#
# A set of functions used by map and plot functions to test for validity of data
#
################################################################################
#' Check for Latitude and Longitude Columns in a Dataframe and return the column values


################################################################################
#' Data: Projectionlist
projlist <- c("aitoff", "albers", "azequalarea", "azequidist",
              "bicentric", "bonne", "conic", "cylequalarea", "cylindrical",
              "eisenlohr", "elliptic", "fisheye", "gall", "gilbert", "guyou",
              "harrison", "hex", "homing", "lagrange", "lambert", "laue", "lune",
              "mercator", "mollweide", "newyorker", "orthographic", "perspective",
              "polyconic", "rectangular", "simpleconic", "sinusoidal", "tetra",
              "trapezoidal")


################################################################################
#' Helper Functions
.check_physeq <- function(physeq){
  #check phyloseq objects for Lat/Lon
  if (!"sam_data" %in% phyloseq::getslots.phyloseq(physeq)){
    stop("Mapping requires that phylos objects have Sample_Data with Latitude and Longitude")
  } 
  #check that sampledata has latitude and longitude columns
  lat <- c('latitude', 'lat', 'lattitude')
  lon <- c('longitude', 'lon', 'long')
  lat_present = FALSE
  lon_present = FALSE
  sampledata <- sample_data(physeq)
  names <- names(sampledata)
  
  for (name in names){
    if( tolower(name) %in% lat){
      latcol <- name
      lat_present <- TRUE
    }
    if( tolower(name) %in% lon){
      loncol <- name
      lon_present <- TRUE
    }
  }
  if (lat_present == FALSE) { stop("sampledata must have a valid latitude column") }
  if (lon_present == FALSE) { stop("sampledata must have a valid longitude column")  }
  list(latcol, loncol)
}
#' Create a basemap from the maps() worldmap focusing on a region
#' projection defaults to mercator, but others can be selected 
#' http://www.inside-r.org/packages/cran/mapproj/docs/mapproject
.create_basemap <-function(region, df, latcol, loncol, proj, parameter, orientation){
  
  # check that the projection is null or is in the projectionlist
  # print out a warning about projections
  if(!is.null(proj)){
    if(!(proj %in% projlist)){
      stop("The projection is not valid. Please use null or one of the following: aitoff, albers, 
        azequalarea, azequidist, bicentric, bonne, conic, cylequalarea, cylindrical, eisenlohr, 
         elliptic, fisheye, gall, gilbert, guyou, harrison, hex, homing, lagrange, lambert, laue, lune,
         mercator, mollweide, newyorker, orthographic, perspective, polyconic, rectangular,
         simpleconic, sinusoidal, tetra, trapezoidal")
    }else{print("You are using a non-default projection that may require additional parameters. 
                See http://www.inside-r.org/packages/cran/mapproj/docs/mapproject for more information")}
  }
  
  
  if(is.null(region)){
    #default worldmap cuts out Antarctica by filtering everythign below -59 Latitude
    world <- ggplot2::map_data("world")
    #world <- world[world$lat > -59,]
    #ToDO: allow subsetting of samples by region. Is there a point-in-polygon library?
    #this is a quick filter based on latitude and longitude not point-in-polygon
    maxlat  = max(world$lat)
    minlat  = min(world$lat)
    maxlong = max(world$long)
    minlong = min(world$long)
    df <- df[ df[, loncol] < maxlong, ]
    df <- df[ df[, loncol] > minlong, ]
    df <- df[ df[, latcol] < maxlat, ]
    df <- df[ df[, latcol] > minlat, ]
  }else if(region=="world"){
    world <- ggplot2::map_data("world")
  }else {
    world <- ggplot2::map_data("world", region = region)
  }
  
  worldmap <- ggplot(world, aes(x=long, y=lat, group=group)) +
    geom_polygon( fill="grey",alpha=0.6) +
    scale_y_continuous(breaks=(-2:2) * 30) +
    scale_x_continuous(breaks=(-4:4) * 45) +
    theme_classic() +
    theme( axis.text = element_blank(), 
           axis.ticks = element_blank(), 
           axis.line = element_blank(), 
           axis.title=element_blank())
  
  if(is.null(proj)) { 
      return(worldmap)
  } else if(is.numeric(parameter)) {
      return(worldmap + coord_map( projection=proj, parameter=parameter ))
  } else if(is.numeric(orientation)) {
      return(worldmap + coord_map( projection=proj, parameter=parameter, orientation=orientation))
  } else {
      return(worldmap + coord_map(projection=proj))
  }
}
#' utility function to check the validity of arguments
.check_names <- function(member, df, allownumeric=FALSE){
  message <- paste(member, " variable must be a valid column name of a Phyloseq table",sep="")
  names   <- names(df)
  if(!is.null(member)){
    if(!allownumeric){
      if(!member %in% names){
        stop(message)
      }
    } else {
      if(!is.numeric(member)){
        if(!member %in% names){
          stop(message)
        }
      }
    }
  }
}
#' utility function to move the x and y positions of the dataset
.jitter_df <- function(df, xcol, ycol, jitter.x, jitter.y){
  df <- data.frame(df)
  dflength <- length(df[,1])
  distx    <- runif(dflength, min= -jitter.x, max=jitter.x)
  disty    <- runif(dflength, min= -jitter.y, max=jitter.y)
  df[xcol] <- df[,xcol] + distx
  df[ycol] <- df[,ycol] + disty
  df
}

.check_NA <- function(df, col){
  colvals <- df[col]
  colvals[ colvals == "None"] <- NA  #some data has "None" so be sure to replace with NA
  truth    <- lapply(colvals, is.na)
  if( any(as.character(truth) == TRUE)){
    warning(paste("Null Values in ",col, sep=""))
    df <- df[ !is.na(df[col,]), ]
  }
  df
}
.coerce_numeric <- function(df, col){
  df[col] <- lapply( lapply(df[col], as.character), as.numeric)
  df
}
#' Utility Function for Converting Distance Matrices to 
#' three column distances while removing all of the duplicates
#' lifted/modified from here: https://github.com/joey711/phyloseq/blob/master/R/plot-methods.R
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
