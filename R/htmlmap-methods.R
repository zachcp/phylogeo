################################################################################
#' Draw An Interactive Leaflet Map from a Phyloseq Object
#'
#' In this case, edges in the network are created if the distance between
#' nodes is below a potentially arbitrary threshold,
#' and special care should be given to considering the choice of this threshold.
#'
#' @return an HTML/JS Leaflet Map
#' @seealso http://rstudio.github.io/leaflet/
#' 
#' @param physeq (Required). 
#'  The name of the phyloseq object. This must have sample data with 
#'  Latitude and Longitude Columns.
#'  
#' @param size (Optional). Default \code{NULL}. 
#'  The size of the vertex points."Abundance" is a special code that will scale 
#'  points by the number of reads in a sample
#' 
#' @export
#' @examples 
#' 
#' data(batfecal)
#' htmlmap_phyloseq(batfecal, size=3)
#' data(batmicrobiome)
#' htmlmap_phyloseq(batmicrobiome, jitter=TRUE, color="SCIENTIFIC_NAME")
htmlmap_phyloseq <- function(physeq, size=NULL, color=NULL){
  
  #install leaflet using devtools
  #replace when leaflet is on CRAN
  if(!require("leaflet")){
    devtools::install_github('rstudio/leaflet')
    library("leaflet")
  }
  
  #get data
  data = phyloseq::sample_data(physeq)
  
  #customize circle size
  if(size == "Abundance"){
    data$circlesize <- phyloseq::sample_sums(physeq)
  }else if(!is.null(size)) {
    data$circlesize <- size
  }
  
  #basemap
  map = leaflet(data) %>% addTiles()
  
  #add custom circles
  if(!is.null(size)) {
    map <- map %>% addCircleMarkers(radius=~circlesize)
  }else{
    map <- map %>% addCircles()
  }

  return(map)
}

#############################################################################
#' Create a Network from the Phyloseq Objects and Draw An HTML Map of the Clusters
#'
#' In this case, edges in the network are created if the distance between
#' nodes is below a potentially arbitrary threshold,
#' and special care should be given to considering the choice of this threshold.
#'   
#' @return an htmlwidgets plot
#' @param physeq (Required). 
#'  The name of the phyloseq object. This must have sample data with 
#'  Latitude and Longitude Columns.
#'  
#' @param igraph  (Optional). Default \code{NULL}
#'  An optional igraph object. Will reduce plotting time to use 
#'  a precalculated network 
#'   
#' @param color (Optional). Default \code{NULL}.
#'  The name of the sample variable in \code{physeq} to use for color mapping
#'  of points (graph vertices). Note: "cluster" can be used to show igraph 
#'  clusters
#'  
#' @param size (Optional). Default \code{4}. 
#'  The size of the vertex points.
#'  
#' @param lines (Optional). Default \code{FALSE}. 
#'  Boolean value. Determines whether lines are drawn between samples
#'  
#' @param distance (Optional). Default \code{"jaccard"}. 
#'  Distance metric used to calculate between-sample distances.
#'  
#' @param maxdist (Optional). Default \code{0.9}. 
#'  Cutoff of the \code{distance} used to detmine whether a sample is 
#'  included in the network.
#'  
#' @param line_weight (Optional). Default \code{0.3}.
#'  The line thickness to use to label graph edges.
#' 
#' @param line_color (Optional). Default \code{Black}.
#'  The name of the sample variable in \code{physeq} to use for color mapping
#'  of lines (graph edges).
#' 
#' @param line_alpha (Optional). Default \code{0.4}.
#'  The transparency level for graph-edge lines.
#'  
#' @param base_data (Optional). Default \code{FALSE}.
#'  Boolean to determine whether to include dat points that aren't in a network.
#' 
#' @param base_data_color (Optional). Default \code{grey}.
#'  named color to determine base data color
#'  
#' \href{https://joey711.github.io/phyloseq/distance}{phyloseq's distance command}.
#' 
#' @import phyloseq
#' @importFrom igraph get.data.frame
#' @importFrom igraph get.vertex.attribute
#' @importFrom igraph clusters  
#' @export
#' @examples
#' library(phyloseq)
#' data(batfecal)
#' htmlmap_network(batfecal)
#' htmlmap_network(batfecal)
#' htmlmap_network(batfecal, maxdist=0.9)
#' data(batmicrobiome)
#' htmlmap_network(batmicrobiome, lines=TRUE)
#' ig <- make_network(batmicrobiome)
#' htmlmap_network(batmicrobiome, igraph= ig )
#' htmlmap_network(batmicrobiome, igraph= ig, color="SCIENTIFIC_NAME")
#' htmlmap_network(batmicrobiome, igraph= ig, lines=TRUE, color="SCIENTIFIC_NAME")
#' htmlmap_network(batmicrobiome, igraph= ig, color="SCIENTIFIC_NAME")
htmlmap_network <- function(physeq, igraph=NULL, maxdist=0.9, distance="jaccard", 
                        color=NULL,  size=4, alpha = 0.8,line_weight=1, 
                        line_color ="Black",line_alpha=0.4 , base_data=FALSE, 
                        base_data_color="grey"){
  
  #install leaflet using devtools
  #replace when leaflet is on CRAN
  if(!require("leaflet")){
    devtools::install_github('rstudio/leaflet')
    library("leaflet")
  }
  
  #helper functions to calculate membership in clusters or lines
  ##############################################################################
  get_clusters <- function(num, graph=igraph){
    #get cluster membership info from igraph object from cluster with clusterid of 'num'
    clusts  <- igraph::clusters(graph)
    members <- which(clusts$membership == num) #get membership
    names   <- get.vertex.attribute(graph, 'name', members)
    df = data.frame(names)
    df['cluster'] <- as.character(num)
    rownames(df) <- df$names
    df    #return a df with name/cluster columns
  }
  
  get_lines <- function(graph=igraph, df=data){
    #get each edge of the network and return a list of dataframes with the node info
    
    getline_df <- function(i, l=links, df1=df){
      #subset data frame using node info
      temprow   <- l[i,]
      tempnames <-c(temprow$from,temprow$to)
      smalldf <- df1[rownames(df1) %in% tempnames, ]
      smalldf['link'] <- i
      smalldf
    }
    
    links <- get.data.frame(graph)
    links_range <- seq( 1:dim(links)[1])
    lines_dfs <- Map(getline_df, links_range)
    lines_df <- Reduce(rbind, lines_dfs)
    lines_df
  }
  
  addlines <- function(map, df, latcol, loncol){
    #df must have link column dennoting belonging to the same line
    lines = unique(df$link)
    for (line in lines){
      line_df <- df[df$link == line,]
      line_df['LAT'] <- line_df[latcol]
      line_df['LON'] <- line_df[loncol]
      
      map <- map %>% addPolylines(data = line_df,
                                  lng  = ~LON,
                                  lat  = ~LAT,
                                  fill = TRUE)
    }
    return(map)
  }
  
  #####################################

  #make network, get cluster information, and add thamesat to the  original dataframe. 
  if(is.null(igraph)){
    igraph <- make_network(physeq, max.dist = maxdist, distance=distance)
  }else{
    if( !"igraph" %in% class(igraph) ){
      stop("igraph must be an igraph network object")} 
  }
  #check basic physeq and lat/lon and make clusters
  latlon <- phylogeo:::.check_physeq(physeq)
  latcol <- as.character( latlon[1] )
  loncol <- as.character( latlon[2] )
  
  clusts <- seq( igraph::clusters(igraph)$no )
  clustdf <- Reduce( rbind, Map(get_clusters, clusts))
  mdf <- merge(clustdf, data.frame(sample_data(physeq)), by="row.names", all.x=TRUE)
  rownames(mdf) <- mdf$Row.names
  
  #create map
  ############################################
  map = leaflet(mdf) %>% addTiles()
  
  #addlines
  linedf <- get_lines(df=mdf)
  map <- addlines(map, linedf, latcol, loncol)
  return(map)
}