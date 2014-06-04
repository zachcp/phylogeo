#' Draw A Map from a Phyloseq Object
#'
#' In this case, edges in the network are created if the distance between
#' nodes is below a potentially arbitrary threshold,
#' and special care should be given to considering the choice of this threshold.
#'
#' @usage map_phylo(physeq, type="samples", 
#'   color=NULL, shape=NULL, point_size=4, alpha=1,
#'   label="value", hjust = 1.35, 
#'   line_weight=0.5, line_color=color, line_alpha=0.4,
#' 	layout.method=layout.fruchterman.reingold, title=NULL)
map_phyloseq <- function(physeq, region=NULL, color=NULL, pointalpha = 0.8, pointsize=4){
  #check that phylobject has sample data
  if (!"sam_data" %in% getslots.phyloseq(physeq)){
    stop("Mapping requires that phyloseq objects have Sample_Data with Latitude and Longitude ")
  } 
  sampledata <- sample_data(physeq)
  
  #check that sampledata has latitude and longitude columns
  lat <- c('latitude', 'lat', 'lattitude')
  lon <- c('longitude', 'lon', 'long')
  lat_present = FALSE
  lon_present = FALSE
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
  
  #basemap
  if (!is.null(region)){
    world <- map_data("world", region = region)
    
    #ToDO: allow subsetting of sampels by region. Is there a point-in-polygon library?
    #this is a quick filter based on latitude and longitude not point-in-polygon
    maxlat  = max(world$lat)
    minlat  = min(world$lat)
    maxlong = max(world$long)
    minlong = min(world$long)
    data2 <- data2[ data2$Longitude < maxlong, ]
    data2 <- data2[ data2$Longitude > minlong, ]
    data2 <- data2[ data2$Latitude < maxlat, ]
    data2 <- data2[ data2$Latitude > minlat, ]
    
  } else {
    world <- map_data("world")
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
  
  data = sample_data(physeq)
  worldmap + geom_point(data=data, aes_string( x=loncol, y=latcol, group=names(data)[1], color=color), size=pointsize, alpha= pointalpha)
  
}
#' Draw A Map from a Phyloseq Object
#'
#' In this case, edges in the network are created if the distance between
#' nodes is below a potentially arbitrary threshold,
#' and special care should be given to considering the choice of this threshold.
#'
#' @usage draw_map(physeq, type="samples", 
#'   color=NULL, shape=NULL, point_size=4, alpha=1,
#'   label="value", hjust = 1.35, 
#' 	line_weight=0.5, line_color=color, line_alpha=0.4,
#' 	layout.method=layout.fruchterman.reingold, title=NULL)
draw_map <- function(network, sampledata, region=NULL, pointsize=4, pointalpha = 0.8){
  #get data and plot
  samples <- as.character(network$data$value) 
  data2 <- sampledata[ row.names(sampledata) %in% samples ]
  
  #basemap
  if (!is.null(region)){
    world <- map_data("world", region = region)
    
    #ToDO: allow subsetting of sampels by region. Is there a point-in-polygon library?
    #this is a quick filter based on latitude and longitude not point-in-polygon
    maxlat  = max(world$lat)
    minlat  = min(world$lat)
    maxlong = max(world$long)
    minlong = min(world$long)
    data2 <- data2[ data2$Longitude < maxlong, ]
    data2 <- data2[ data2$Longitude > minlong, ]
    data2 <- data2[ data2$Latitude < maxlat, ]
    data2 <- data2[ data2$Latitude > minlat, ]
    
  } else {
    world <- map_data("world")
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
  
  worldmap + geom_point(data=data2, aes( x=Longitude, y=Latitude, color=Geotype, group=Geotype), size=pointsize, alpha= pointalpha)
  
}
#' Plot a network using ggplot2 (represent microbiome)
#'
#' There are many useful examples of phyloseq network graphics in the
#' \href{http://joey711.github.io/phyloseq/plot_network-examples}{phyloseq online tutorials}.
#' A custom plotting function for displaying networks
#' using advanced \code{\link[ggplot2]{ggplot}}2 formatting.
#' The network itself should be represented using
#' the \code{igraph} package.
#' For the \code{\link{phyloseq-package}} it is suggested that the network object
#' (argument \code{g})
#' be created using the
#'  \code{\link{make_network}} function, 
#' and based upon sample-wise or taxa-wise microbiome ecological distances 
#' calculated from a phylogenetic sequencing experiment 
#' (\code{\link{phyloseq-class}}).
#' In this case, edges in the network are created if the distance between
#' nodes is below a potentially arbitrary threshold,
#' and special care should be given to considering the choice of this threshold.
#'
#' @usage plot_network(g, physeq=NULL, type="samples", 
#'   color=NULL, shape=NULL, point_size=4, alpha=1,
#'   label="value", hjust = 1.35, 
#' 	line_weight=0.5, line_color=color, line_alpha=0.4,
#' 	layout.method=layout.fruchterman.reingold, title=NULL)
#'
map_richness <- function() {
  
}
#allow facetting of subcategories subsamples - facet_wrap(~Phylum)
#plot connections between samples at different values
# allow coloring by either a sample_data field
# or by a an igraph cluster
# make sure to allow facteting of these structures - like facte by the cluster type
# or environemtnela type
map_network  <- function() {
  
}
#plot a tree and the location of the samples
plot_tree    <- function() {
  
}
# do this?
plot_heatmap <- function() {
  
}
##
points_from_network <- function(network, sampledata) {
  samples <- as.character(network$data$value) 
  data <- sampledata[ row.names(sampledata) %in% samples ]
  points <- points(data$Longitude, data$Latitude, pch=20, col="gray50", cex=1) 
}