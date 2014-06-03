library(phyloseq)
library(ggplot2)
library(igraph)
library(maps)

#plot richness with points as size of richness estimates
map_richness <- function() {}
#allow facetting of subcategories subsamples - facet_wrap(~Phylum)


#plot connections between samples at different values
map_network  <- function() {}

#plot a tree and the location of the samples
plot_tree    <- function() {}

#
plot heatmap <- function() {}

points_from_network <- function(network, sampledata) {
  samples <- as.character(network$data$value) 
  data <- sampledata[ row.names(sampledata) %in% samples ]
  points <- points(data$Longitude, data$Latitude, pch=20, col="gray50", cex=1) 
}

draw_map <- function(network, sampledata, region=NULL, pointsize=4, pointalpha = 0.8){
  #get data and ploe
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