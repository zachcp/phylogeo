#' Draw A Map from a Phyloseq Object
#'
#' In this case, edges in the network are created if the distance between
#' nodes is below a potentially arbitrary threshold,
#' and special care should be given to considering the choice of this threshold.
#'
#' @usage map_phylseq(physeq, type="samples",  region=NULL,
#'   color=NULL, shape=NULL, pointsize=4, pointalpha=1,
#'   label="value", hjust = 1.35, 
#'   line_weight=0.5, line_color=color, line_alpha=0.4,
#' 	layout.method=layout.fruchterman.reingold, title=NULL)
map_phyloseq <- function(physeq, region=NULL, color=NULL, pointalpha = 0.8, pointsize=NULL){
  #check basic physeq and lat/lon
  latlon <- .check_physeq(physeq)
  latcol <- latlon[1]
  loncol <- latlon[2]
  
  sampledata <- sample_data(physeq)
  names <- names(sampledata)
  
  #check that color and pointsize are present and that pointsize is numeric 
  if (!is.null(color){
    if(!color %in% names) { stop("color variable must be a sampledata column") }
  }
  if (!is.null(color){
    if (!pointsize %in% names) { stop("color variable must be a sampledata column") }
  }
  if (sample_data[pointsize][1] %in% names) { stop("color variable must be a sampledata column") }
  
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
  #how to hande when pointsize information can be either global (outside of aes), orper-sample (inseide of aes)
  if (!is.null(pointsize) ) {
    worldmap <- worldmap + geom_point(data=data, aes_string( x=loncol, y=latcol, group=names(data)[1], color=color, size = pointsize), alpha= pointalpha)    
  } else {
    worldmap <- worldmap + geom_point(data=data, aes_string( x=loncol, y=latcol, group=names(data)[1], color=color), size=4 ,alpha=pointalpha)    
  }
  worldmap
}
#' Create a Network from the Phyloseq Objects and Draw A Map of the Clusters
#'
#' In this case, edges in the network are created if the distance between
#' nodes is below a potentially arbitrary threshold,
#' and special care should be given to considering the choice of this threshold.
#'
#' @usage dmap_network(physeq,maxdist=0.9, distance="jaccard"
#'   color=NULL, shape=NULL, region=NULL, point_size=4, alpha=1,
#'   label="value", hjust = 1.35, 
#' 	line_weight=0.5, line_color=color, line_alpha=0.4,
#' 	layout.method=layout.fruchterman.reingold, title=NULL)
map_network <- function(physeq, maxdist=0.9, distance="jaccard", color=NULL, region=NULL, pointsize=4, pointalpha = 0.8){
  latlon <- .check_physeq(physeq)
  latcol <- latlon[1]
  loncol <- latlon[2]
  sampledata <- sample_data(physeq)
  
  #helper function to assign names in the df to clusters
  get_clusters <- function(num, graph=ig){
    #num is the index of the cluster
    clusts  <- clusters(graph)
    members <- which(clusts$membership == num) #get membership
    names   <- get.vertex.attribute(graph, 'name', members)
    df = data.frame(names)
    df['cluster'] <- num
    rownames(df) <- df$names
   #return a df with name/cluster columns
    df
  }
  
  #create df from cluster membership
  ig <- make_network(physeq, max.dist = maxdist, distance=distance)
  clusts <- seq(clusters(ig)$no)
  clustdf <- Reduce( rbind, Map(get_clusters, clusts))
    
  #merge the original dataframe with the cluster info
  data <- merge(clustdf, data.frame(sampledata), by="row.names", all.x=T)
  
#   #basemap
#   if (!is.null(region)){
#     world <- map_data("world", region = region)
#     
#     #ToDO: allow subsetting of sampels by region. Is there a point-in-polygon library?
#     #this is a quick filter based on latitude and longitude not point-in-polygon
#     maxlat  = max(world$lat)
#     minlat  = min(world$lat)
#     maxlong = max(world$long)
#     minlong = min(world$long)
#     data2 <- data2[ data2$Longitude < maxlong, ]
#     data2 <- data2[ data2$Longitude > minlong, ]
#     data2 <- data2[ data2$Latitude < maxlat, ]
#     data2 <- data2[ data2$Latitude > minlat, ]
#     
#   } else {
#     world <- map_data("world")
#   }
  
  world <- map_data("world")
  worldmap <- ggplot(world, aes(x=long, y=lat, group=group)) +
    geom_polygon( fill="grey",alpha=0.6) +
    scale_y_continuous(breaks=(-2:2) * 30) +
    scale_x_continuous(breaks=(-4:4) * 45) +
    theme_classic() +
    theme( axis.text = element_blank(), 
           axis.ticks = element_blank(), 
           axis.line = element_blank(), 
           axis.title=element_blank())
worldmap <- worldmap + geom_point(data=data, aes_string( x=loncol, y=latcol, group=names(data)[1], color=color), size=4 ,alpha=pointalpha)    
  worldmap + geom_point(data=data, aes( x=Longitude, y=Latitude, color=cluster, group=Geotype), size=pointsize, alpha= pointalpha)
  
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
#
.check_physeq <- function(physeq){
  #check phyloseq objects for Lat/Lon
  if (!"sam_data" %in% getslots.phyloseq(physeq)){
    stop("Mapping requires that phylos objects have Sample_Data with Latitude and Longitude ")
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