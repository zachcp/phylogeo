################################################################################
#' Draw An Interactive Leaflet Map from a Phyloseq Object
#'
#' In this case, edges in the network are created if the distance between
#' nodes is below a potentially arbitrary threshold,
#' and special care should be given to considering the choice of this threshold.
#'
#' @return  \code{\link[htmlwidgets]{htmlwidgets}} Leaflet Map
#' @seealso http://rstudio.github.io/leaflet/
#' 
#' @param physeq (Required). 
#'  The name of the \code{\link[phyloseq]{phyloseq}} phyloseq object. This must have sample data with 
#'  Latitude and Longitude Columns.
#'  
#' @param size (Optional). Default \code{NULL}. 
#'  The size of the vertex points."Abundance" is a special code that will scale 
#'  points by the number of reads in a sample
#' 
#' @param color (Optional). Default \code{"blue"}. 
#'  Color of points using color names (e.g. red, blue)
#' 
#' @seealso
#'  \code{\link[leaflet]{leaflet}}
#'  \code{\link[phyloseq]{phyloseq}}
#' 
#' @import dplyr
#' @import phyloseq
#' @importFrom dplyr %>%
#' @export
#' @examples 
#' data(batfecal)
#' htmlmap_phyloseq(batfecal, size=3)
#' data(batmicrobiome)
#' htmlmap_phyloseq(batmicrobiome, color="blue")
htmlmap_phyloseq <- function(physeq, size=NULL, color="blue"){
  
  #install leaflet using devtools
  #replace when leaflet is on CRAN
  if(!require("leaflet")){
    devtools::install_github('rstudio/leaflet')
    library("leaflet")
  }
  
  #get data
  data = sample_data(physeq)
  
  #customize circle size
  if(!is.null(size)){
    if(size == "Abundance") data$circlesize <- phyloseq::sample_sums(physeq) else
                            data$circlesize <- size
  }
  #basemap
  map = leaflet(data) %>% addTiles()
  
  #add custom circles
  if(!is.null(size)) map <- map %>% addCircleMarkers(radius=~circlesize, color = color) else
                     map <- map %>% addCircles(color = color)
  
  return(map)
}

#############################################################################
#' Create a Network from the Phyloseq Objects and Draw An HTML Map of the Clusters
#'
#' In this case, edges in the network are created if the distance between
#' nodes is below a potentially arbitrary threshold,
#' and special care should be given to considering the choice of this threshold.
#'   
#' @return an  \code{\link[htmlwidgets]{htmlwidgets}} plot
#' @seealso http://rstudio.github.io/leaflet/
#' 
#' @param physeq (Required). 
#'  The name of the \code{\link[phyloseq]{phyloseq}} phyloseq object. This must have sample data with 
#'  Latitude and Longitude Columns.
#'  
#' @param igraph  (Optional). Default \code{NULL}
#'  An optional \code{\link[igraph]{igraph}} igraph object. Will reduce plotting time to use 
#'  a precalculated network 
#'  
#' @param distance (Optional). Default \code{"jaccard"}. 
#'  Distance metric used to calculate between-sample distances.
#'  
#' @param maxdist (Optional). Default \code{0.9}. 
#'  Cutoff of the \code{distance} used to detmine whether a sample is 
#'  included in the network.
#'   
#' @param line_color (Optional). Default \code{black}.
#'  The name of the sample variable in \code{physeq} to use for color mapping
#'  of lines (graph edges).
#' 
#' @param line_alpha (Optional). Default \code{0.4}.
#'  The transparency level for graph-edge lines.
#'  
#' @param line_weight (Optional). Default \code{1}.
#'  The line thickness to use to label graph edges.
#'  
#' @param color (Optional). Default \code{black}.
#'   Color of the points
#'   
#' @param circle_alpha (Optional). Default \code{1}. 
#'  The opacity of the points.
#'  
#' @param fill (Optional). Default \code{TRUE}.
#'  Boolean. Whether to fill in the points or not.
#'  
#' @param fillOpacity (Optional). Default \code{1}.
#'  opacity of circle fills
#'  
#' @param fillColor (Optional). Default \code{color}.
#'  Color to be used for filling in the circle
#'  
#' @param size (Optional). Default \code{1}. 
#'  The size of the vertex points. If "Abundance" is supplied as the argument
#'  the size will be scaled to the abundance of the OTUs in the sample.
#'  
#'
#'  @seealso
#'    \code{\link[leaflet]{leaflet}}
#'    \code{\link[phyloseq]{phyloseq}}
#'  
#'  @seealso  
#'    \href{https://joey711.github.io/phyloseq/distance}{phyloseq's distance command}.
#' 
#' @import phyloseq
#' @import dplyr
#' @importFrom dplyr %>%
#' @importFrom igraph get.data.frame
#' @importFrom igraph get.vertex.attribute
#' @importFrom igraph clusters  
#' @export
#' @examples
#' htmlmap_network(batfecal)
#' htmlmap_network(batfecal, maxdist=0.9)
#' 
#' htmlmap_network(batmicrobiome, maxdist=0.5)
#' ig <- make_network(batmicrobiome)
#' htmlmap_network(batmicrobiome, igraph= ig)
#' htmlmap_network(epoxamicin_KS, maxdist=0.99, line_color = "red", line_weight = 4, line_alpha=0.5)
htmlmap_network <- function(physeq, 
                            #distance related 
                            igraph=NULL, 
                            maxdist=0.9, 
                            distance="jaccard",
                            #linerelated
                            line_color ="black",
                            line_alpha=0.4 , 
                            line_weight=1,
                            #point related
                            color="blue", 
                            circle_alpha = 0.8,
                            fill = FALSE,
                            fillOpacity = 1,
                            fillColor = color,
                            size=1){
  
  # install leaflet using devtools replace when leaflet is on CRAN
  #############################################################################
  if(!require("leaflet")){
    devtools::install_github('rstudio/leaflet')
    library("leaflet")
  }
  
  #helper functions to calculate membership in clusters or lines
  #############################################################################
  get_clusters <- function(num, graph=igraph){
    #get cluster membership info from igraph object 
    #from cluster with clusterid of 'num'
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
      tempnames <- c(temprow$from,temprow$to)
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
      #currently I add extra columns here so I can quote them below.
      # it would be cleaner to call directly.
      # I also have to make sure the columns are numeric
      line_df <- df[df$link == line,]
      line_df['LAT'] <- line_df[latcol]
      line_df['LON'] <- line_df[loncol]
      line_df$LON <- as.numeric(as.character(line_df$LON))
      line_df$LAT <- as.numeric(as.character(line_df$LAT))
      
      map <- map %>% addPolylines(data = line_df,
                                  lng  = ~LON,
                                  lat  = ~LAT,
                                  color = line_color,
                                  weight = line_weight,
                                  opacity = line_alpha,
                                  fill = fill,
                                  fillOpacity =fillOpacity,
                                  fillColor = fillColor)
    }
    return(map)
  }
  
  ############################################################################

  #make network, get cluster information, and add thamesat to the  original dataframe. 
  if(is.null(igraph)){
    igraph <- make_network(physeq, max.dist = maxdist, distance=distance)
  }else{
    if( !"igraph" %in% class(igraph) ){
      stop("igraph must be an igraph network object")} 
  }
  
  #check basic physeq and lat/lon and make clusters
  latlon <- .check_physeq(physeq)
  latcol <- as.character( latlon[1] )
  loncol <- as.character( latlon[2] )
  #get clusters and make a dataframe from them
  clusts <- seq( igraph::clusters(igraph)$no )
  clustdf <- Reduce( rbind, Map(get_clusters, clusts))
  #get sample data
  data <- data.frame(sample_data(physeq))
  
  #customize circle size prior to using the cluster
  if(size == "Abundance") data$circlesize <- phyloseq::sample_sums(physeq) else
                          data$circlesize <- size
  
  #merge sample data with cluster data
  mdf <- merge(clustdf, data, by="row.names", all.x=TRUE)
  rownames(mdf) <- mdf$Row.names
  
  #create map
  ############################################
  map = leaflet(mdf) %>% addTiles()
  
  #addlines
  linedf <- get_lines(df=mdf)
  map <- addlines(map, linedf, latcol, loncol)
  
  #add points
  map <- map %>% addCircleMarkers(radius=~circlesize, 
                                  color=color, 
                                  opacity = circle_alpha, 
                                  fillOpacity = fillOpacity)
  
  return(map)
}