###########################################################################################################
#' Draw A Map from a Phyloseq Object
#'
#' In this case, edges in the network are created if the distance between
#' nodes is below a potentially arbitrary threshold,
#' and special care should be given to considering the choice of this threshold.
#'
#' @usage map_phyloseq(physeq, region=NULL, color=NULL, shape=NULL,
#'   point_size=4, alpha=0.5 )
#'
#' @param physeq (Required). 
#'  The name of the phyloseq object. This must have sample data with Latitude and Longitude Columns.
#'  
#' @param region (Optional). Default \code{NULL}.
#'  The name of geographic region that can be used to zoom.
#' 
#' @param color (Optional). Default \code{NULL}.
#'  The name of the sample variable in \code{physeq} to use for color mapping
#'  of points (graph vertices).
#'  
#' @param shape (Optional). Default \code{NULL}.
#'  The name of the sample variable in \code{physeq} to use for shape mapping.
#'  of points (graph vertices).
#'  
#' @param point_size (Optional). Default \code{4}. 
#'  The size of the vertex points.
#'  
#' @param alpha (Optional). Default \code{0.8}. 
#'  A value between 0 and 1 for the alpha transparency of the vertex points.
#'  
#' @param jitter (Optional). Default \code{False}. 
#'  Determines whether or not to jitter your points.
#'
#' @param jitter.x (Optional). Default \code{3}. 
#'  Value for X jitter
#'
#' @param jitter.y (Optional). Default \code{3}. 
#'  Value for Y jitter
#'
#'  
#' @import ggplot2  
#' @export
#' @examples 
#' 
#' data(batfecal)
#' map_phyloseq(batfecal, region="china", jitter=T, jitter.x=2,jitter.y=2, color="PH")
#' data(batmicrobiome)
#' map_phyloseq(batmicrobiome, jitter=TRUE, color="SCIENTIFIC_NAME")
map_phyloseq <- function(physeq, region=NULL, color=NULL, shape=NULL, point_size=4, alpha = 0.8, 
                         jitter=FALSE, jitter.x=3, jitter.y=3){
  #check basic physeq and lat/lon
  latlon <- .check_physeq(physeq)
  latcol <- as.character( latlon[1] )
  loncol <- as.character( latlon[2] )
  data   <- data.frame( sample_data(physeq) )
  data   <- .check_NA(data, latcol)
  data   <- .coerce_numeric(data,latcol)
  data   <- .check_NA(data, loncol)
  data   <- .coerce_numeric(data,loncol)
  names  <- names(data)
  
  
  #check plot options
  .check_names(color,data)
  .check_names(point_size,data, allownumeric=T)
  
  #create map
  ############################################################################################################
  worldmap <- .create_basemap(region=region, df=data, latcol=latcol,loncol=loncol)
  
  if(jitter){
    data <- .jitter_df(df=data,xcol=loncol,ycol=latcol,jitter.x=jitter.x,jitter.y=jitter.y)
  }
  
  #how to hande when point_size information can be either global (outside of aes), orper-sample (inseide of aes)
  if(is.numeric(point_size)){
    worldmap <- worldmap + geom_point(data=data, aes_string( x=loncol, y=latcol, group=NULL, color=color), 
                                      size = point_size, alpha= alpha) 
  }else{
    worldmap <- worldmap + geom_point(data=data, aes_string( x=loncol, y=latcol, group=NULL, color=color, size = point_size),
                                      alpha= alpha) 
  } 
  
  worldmap
}
###########################################################################################################
#' Create a Network from the Phyloseq Objects and Draw A Map of the Clusters
#'
#' In this case, edges in the network are created if the distance between
#' nodes is below a potentially arbitrary threshold,
#' and special care should be given to considering the choice of this threshold.
#'
#' @usage map_network(physeq,maxdist=0.9, distance="bray", point_size=6, alpha=0.5, lines=TRUE,
#'  line_color= "Red", line_alpha=0.4, title="Awesome Netwrok Graph")
#'   
#' @param physeq (Required). 
#'  The name of the phyloseq object. This must have sample data with Latitude and Longitude Columns.
#'  
#' @param igraph  (Optional). Default \code{NULL}
#'  An optional igraph object. Will reduce plotting time to use a precalculated network 
#'  
#' @param region (Optional). Default \code{NULL}.
#'  The name of geographic region that can be used to zoom.
#' 
#' @param color (Optional). Default \code{NULL}.
#'  The name of the sample variable in \code{physeq} to use for color mapping
#'  of points (graph vertices). Note: "cluster" can be used to show igraph clusters
#'  
#' @param shape (Optional). Default \code{NULL}.
#'  The name of the sample variable in \code{physeq} to use for shape mapping.
#'  of points (graph vertices).  Note: "cluster" can be used to show igraph clusters
#'  
#' @param point_size (Optional). Default \code{4}. 
#'  The size of the vertex points.
#'  
#' @param alpha (Optional). Default \code{0.8}. 
#'  A value between 0 and 1 for the alpha transparency of the vertex points.
#'
#' @param lines (Optional). Default \code{FALSE}. 
#'  Boolean value. Determines whether lines are drawn between samples
#'  
#' @param distance (Optional). Default \code{"jaccard"}. 
#'  Distance metric used to calculate between-sample distances.
#'  
#' @param maxdist (Optional). Default \code{0.9}. 
#'  Cutoff of the \code{distance} used to detmine whether a sample is included in the network.
#'  
#' @param jitter (Optional). Default \code{FALSE}. 
#'  Boolean value. Determines whether to use geom_jitter() to reposition points.
#'  
#' @param jitter.x (Optional). Default \code{5}. 
#'  Value to jitter in the X direction if using \code{jitter}
#' 
#' @param jitter.y (Optional). Default \code{5}. 
#'  Value to jitter in the Y direction if using \code{jitter}
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
#'  @param base_data_color (Optional). Default \code{grey}.
#'  named color to determine base data coloe
#'
#' @import ggplot2
#' @importFrom igraph get.data.frame
#' @importFrom igraph get.vertex.attribute
#' @importFrom igraph clusters  
#' @export
#' @examples
#' data(batfecal)
#' map_network(batfecal)
#' map_network(batfecal, region="china", jitter=TRUE, lines=TRUE)
#' map_network(batfecal, region="china", jitter=TRUE, lines=TRUE, maxdist=0.9)
#' data(batmicrobiome)
#' map_network(batmicrobiome, lines=TRUE)
#' ig <- make_network(batmicrobiome)
#' map_network(batmicrobiome, igraph= ig, lines=TRUE)
#' map_network(batmicrobiome, igraph= ig, lines=TRUE, color="SCIENTIFIC_NAME")
#' map_network(batmicrobiome, igraph= ig, lines=TRUE, color="SCIENTIFIC_NAME", jitter=TRUE)
map_network <- function(physeq, igraph=NULL, maxdist=0.9, distance="jaccard", color=NULL, region=NULL, 
                        point_size=4, alpha = 0.8, jitter=FALSE, jitter.x=3, jitter.y=3, shape=NULL, 
                        lines=FALSE, line_weight=1, line_color ="Black" ,line_alpha=0.4 , base_data=FALSE, base_data_color="grey"){

  #helper functions to calculate membership in clusters or lines
  ######################################################################################################
  get_clusters <- function(num, graph=igraph){
    #get cluster membership info from igraph object from cluster with clusterid of "num
    clusts  <- clusters(graph)
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
  
  #add lines to ggplot object
  draw_lines <- function(plt, df2){
    df2 <- data.frame(df2) #to ensure list returns a df object
    plt <- plt + geom_line(data=df2,  aes_string( x=loncol, y=latcol, group=names(df2)[1]))
  }
  #####################################
  
  #check basic physeq and lat/lon
  latlon <- .check_physeq(physeq)
  latcol <- as.character( latlon[1] )
  loncol <- as.character( latlon[2] )
  data <- data.frame( sample_data(physeq) )
  data <- .check_NA(data, latcol)
  data <- .coerce_numeric(data,latcol)
  data <- .check_NA(data, loncol)
  data <- .coerce_numeric(data,loncol)
  names <- names(data)
  
  #make network, get cluster information, and add thamesat to the  original dataframe. 
  if(is.null(igraph)){
    igraph <- make_network(physeq, max.dist = maxdist, distance=distance)
  }else{
    if( !"igraph" %in% class(igraph) ){
      stop("igraph must be an igraph network object")} 
  }
  clusts <- seq(clusters(igraph)$no)
  clustdf <- Reduce( rbind, Map(get_clusters, clusts))
  mdf <- merge(clustdf, data.frame(data), by="row.names", all.x=T)
  rownames(mdf) <- mdf$Row.names
  
  #check plot options
  .check_names(color,mdf)
  .check_names(point_size,mdf, allownumeric=T)
  
  #create map
  ############################################
  
  #basemap
  worldmap <- .create_basemap(region=region, df=mdf, latcol=latcol, loncol=loncol)
  
  #modify points if using jitter
  if(jitter){
    mdf <- .jitter_df( df=mdf, xcol=loncol, ycol=latcol, jitter.x=jitter.x, jitter.y=jitter.y)
  }
  
  #add points that aren't part of a network
  if(base_data){
    network_points <- rownames(mdf)
    nonetworkdf <- data[!rownames(data) %in% network_points, ] 
    worldmap <- worldmap + geom_point(data = nonetworkdf, aes_string(x=loncol,y=latcol, group=NULL), color=base_data_color, size=point_size)
  }
  
  #addlines
  if(lines){
    linedf <- get_lines(df=mdf) 
    worldmap <- worldmap + geom_line(data=linedf,aes_string(x=loncol,y=latcol, group="link"), size=line_weight, alpha=line_alpha, color=line_color)
  }
 
  #add points
  #how to hande when point_size information can be either global (outside of aes), orper-sample (inseide of aes)
  if(is.numeric(point_size)){
   points <- geom_point(data=mdf, aes_string( x=loncol, y=latcol, group=NULL, color=color, shape=shape), 
                                      size = point_size, alpha= alpha) 
  }else{
    points <- geom_point(data=mdf, aes_string( x=loncol, y=latcol, group=NULL, color=color, size = point_size, shape=shape), 
                                      alpha= alpha) 
  } 
  worldmap <- worldmap + points
  


  points <- geom_point(data=mdf, aes_string( x=loncol, y=latcol, group=NULL, color=NULL, shape=NULL),  size = point_size, alpha= alpha) 
  ###########################
  
  return(worldmap)
}
###########################################################################################################
#' Map a Phyloseq Object while also drawing a phlogenetic tree of the taxa
#'
#' @usage map_tree(physeq, region=NULL, color = NULL, size= NULL, point_size=4, alpha=0.8,jitter= FALSE,
#'                jitter.x=3, jitter.y=3)
#'   
#' @param physeq (Required). 
#'  The name of the phyloseq object. This must have sample data with Latitude and Longitude Columns.
#'  
#' @param color (Optional). Default \code{NULL}.
#'  The name of the sample variable in \code{physeq} to use for color mapping
#'  of points (graph vertices).
#'  
#' @param size (Optional). Default \code{NULL}.
#'  The name of the sample variable in \code{physeq} to use for size mapping.
#'  of points (graph vertices).
#'  
#' @param point_size (Optional). Default \code{4}. 
#'  The size of the vertex points.
#'  
#' @param alpha (Optional). Default \code{0.8}. 
#'  A value between 0 and 1 for the alpha transparency of the vertex points.
#'  
#' @param jitter (Optional). Default \code{False}. 
#'  Determines whether or not to jitter your points.
#'
#' @param jitter.x (Optional). Default \code{3}. 
#'  Value for X jitter
#'
#' @param jitter.y (Optional). Default \code{3}. 
#'  Value for Y jitter

#'  @param width_ratio (Optional). Default \code{2}.
#'  relative widths of tree and map
#'  
#'  @param map_on_left (Optional). Default \code{TRUE}.
#'  determine whether the map is on the left
#'  
#' @import ggplot2
#' @export
#' @examples
#' data(epoxamicin_KS)
#' map_tree(epoxamicin_KS)
#' map_tree(epoxamicin_KS, color="Geotype", jitter=TRUE)
map_tree <- function(physeq,  region=NULL, color = NULL, size= NULL, point_size=4, alpha=0.8,
                    jitter= FALSE, jitter.x=3, jitter.y=3, method = "sampledodge", 
                    nodelabf =nodeplotblank, treesize = NULL, min.abundance = Inf, 
                    label.tips = NULL, text.size = NULL, sizebase = 5, base.spacing = 0.02, ladderize = TRUE,
                    plot.margin = 0.2, title = NULL, treetheme = NULL, justify = "jagged",
                    width_ratio = 2, map_on_left = TRUE) {
    #trim samples that are not in the tree
    physeq2 <- prune_samples(sample_sums(physeq) > 0, physeq)
    
    mapplot  <- map_phyloseq(physeq2, region=region, color= color, point_size=point_size, alpha = alpha, jitter=jitter, 
                              jitter.x=jitter.x, jitter.y=jitter.y)  + 
                              theme(legend.position="none") 
    treeplot <- phyloseq::plot_tree(physeq2, color=color, label.tips=label.tips, text.size=text.size, 
                          sizebase=sizebase, base.spacing = base.spacing, ladderize = ladderize,
                          plot.margin = plot.margin, title = title, treetheme=treetheme, justify = justify, nodelab =nodelabf ) +
                          theme(legend.key = element_rect(fill = "white"))
    # # trim space by setting xlims
    # xvals <- treeplot$data$x
    # xvals <- xvals[!is.na(xvals)]
    # xmin <- min(xvals)
    # xmax <- max(xvals) * 1.5
    # treeplot <- treeplot + xlim( min(xvals), max(xvals))
    
    if(map_on_left){
        combinedplot <- gridExtra::grid.arrange(mapplot + theme(legend.position="none") ,treeplot, ncol=2, widths=c(width_ratio,1))
    } else{
        combinedplot <- gridExtra::grid.arrange(treeplot + theme(legend.position="none"),mapplot, ncol=2, widths=c(1,width_ratio))    
    }
    return(combinedplot)
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
    world <- ggplot2::map_data("world", region = region)
    
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
    world <- ggplot2::map_data("world")
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
  if( any(as.character(truth) == T)){
    warning(paste("Null Values in ",col, sep=""))
    df <- df[ !is.na(df[col,]), ]
  }
  df
}
.coerce_numeric <- function(df, col){
  df[col] <- lapply( lapply(df[col], as.character), as.numeric)
  df
}