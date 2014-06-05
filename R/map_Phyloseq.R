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
#'
#' @import ggplot2  
#' @export
#' @examples 
#' data(AD)
#' map_phylo(AD)
#' map_phylo(AD, region="bra") 
#' map_phylo(AD, color="Geotype", pointsize="richness") 
map_phyloseq <- function(physeq, region=NULL, color=NULL, pointsize=NULL, pointalpha = 0.8){
  #check basic physeq and lat/lon
  latlon <- .check_physeq(physeq)
  latcol <- as.character( latlon[1] )
  loncol <- as.character( latlon[2] )
  data <- sample_data(physeq)
  names <- names(data)
  
  #check plot options
  .check_names(color,mdf)
  .check_names(pointsize,mdf, allownumeric=T)
  
  #create map
  ############################################################################################################
  worldmap <- .create_basemap(region=region, df=data, latcol=latcol,loncol=loncol)
  
  #how to hande when pointsize information can be either global (outside of aes), orper-sample (inseide of aes)
  if(is.numeric(pointsize)){
    worldmap <- worldmap + geom_point(data=mdf, aes_string( x=loncol, y=latcol, group=NULL, color=color),size = pointsize, alpha= pointalpha) 
  }else{
    worldmap <- worldmap + geom_point(data=mdf, aes_string( x=loncol, y=latcol, group=NULL, color=color, size = pointsize), alpha= pointalpha) 
  } 
  
  worldmap
}
#' Create a Network from the Phyloseq Objects and Draw A Map of the Clusters
#'
#' In this case, edges in the network are created if the distance between
#' nodes is below a potentially arbitrary threshold,
#' and special care should be given to considering the choice of this threshold.
#'
#' @usage map_network(physeq,maxdist=0.9, distance="jaccard"
#'   color=NULL, shape=NULL, region=NULL, point_size=4, alpha=1,
#'   label="value", hjust = 1.35, 
#' 	line_weight=0.5, line_color=color, line_alpha=0.4,
#' 	layout.method=layout.fruchterman.reingold, title=NULL)
#'   
#' @import ggplot2
#' @import igraph
#' @importFrom igraph get.data.frame
#' @importFrom igraph get.vertex.attribute
#' @importFrom igraph clusters  
#' @export
#' @examples 
#' data(AD)
#' map_phylo(AD)
#' map_phylo(AD, region="bra") 
#' map_phylo(AD, color="Geotype", pointsize="richness") 
map_network <- function(physeq, maxdist=0.9, distance="jaccard", color=NULL, region=NULL, pointsize=4, pointalpha = 0.8, lines=FALSE){

  #helper functions to calculate membership in clusters or lines
  ######################################################################################################
  get_clusters <- function(num, graph=ig){
    #get cluster membership info from igraph object from cluster with clusterid of "num
    clusts  <- clusters(graph)
    members <- which(clusts$membership == num) #get membership
    names   <- get.vertex.attribute(graph, 'name', members)
    df = data.frame(names)
    df['cluster'] <- num
    rownames(df) <- df$names
    df    #return a df with name/cluster columns
  }
  
  get_lines <- function(graph=ig, df=data){
    #get each edge of the network and return a list of dataframes with the node info
 
    getline_df <- function(i, l=links, df1=df){
      #subset data frmae using node info
      temprow   <- l[i,]
      tempnames <-c(temprow$from,temprow$to)
      smalldf <- df1[rownames(df1) %in% tempnames, ]
      smalldf
    }
    
    links <- get.data.frame(ig)
    links_range <- seq( 1:dim(links)[1])
    lines_dfs <- Map(getline_df, links_range)
    lines_dfs
  }
  
  #add lines to ggplot object
  draw_lines <- function(plt, df2){
    df2 <- data.frame(df2) #to ensure list returns a df object
    plt <- plt + geom_line(data=df2,  aes_string( x=loncol, y=latcol, group=names(df2)[1]))
  }
  ######################################################################################################
  
  #check basic physeq and lat/lon
  latlon <- .check_physeq(physeq)
  latcol <- as.character( latlon[1] )
  loncol <- as.character( latlon[2] )
  data <- sample_data(physeq)
  names <- names(data)
  
  #make network, get cluster information, and add that to the  original dataframe. 
  ig <- make_network(physeq, max.dist = maxdist, distance=distance)
  clusts <- seq(clusters(ig)$no)
  clustdf <- Reduce( rbind, Map(get_clusters, clusts))
  mdf <- merge(clustdf, data.frame(data), by="row.names", all.x=T)
  
  #check plot options
  .check_names(color,mdf)
  .check_names(pointsize,mdf, allownumeric=T)
  
  #create map
  ############################################################################################################
  worldmap <- .create_basemap(region=region, df=mdf, latcol=latcol, loncol=loncol)
 
  #how to hande when pointsize information can be either global (outside of aes), orper-sample (inseide of aes)
  if(is.numeric(pointsize)){
    worldmap <- worldmap + geom_point(data=mdf, aes_string( x=loncol, y=latcol, group=NULL, color=color),size = pointsize, alpha= pointalpha) 
  }else{
    worldmap <- worldmap + geom_point(data=mdf, aes_string( x=loncol, y=latcol, group=NULL, color=color, size = pointsize), alpha= pointalpha) 
  } 

  #addlines
  if(lines){
    linelist <- get_lines()
    worldmap <- worldmap  + lapply(linelist, geom_line, mapping = aes_string(x=loncol,y=latcol, group=NULL))
  }
  ###########################################################################################################
  
  worldmap
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
#
#plot a tree and the location of the samples
plot_tree    <- function() {
  
}
# do this?
plot_heatmap <- function() {
  
}
#
#' Check for Latitude and Longitude Columns in a Dataframe and return the column values
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
#' Create a basemap from the maps() worldmap focusing on a region
.create_basemap <-function(region, df, latcol, loncol){
  if (!is.null(region)){
    world <- map_data("world", region = region)
    
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
  
  worldmap
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