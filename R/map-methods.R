#
# methods for drawing maps of phyloseq objects
#
###############################################################################
#' Draw A Map from a Phyloseq Object
#'
#' In this case, edges in the network are created if the distance between
#' nodes is below a potentially arbitrary threshold,and special care should 
#' be given to considering the choice of this threshold.
#'
#' @return a ggplot object
#' @param physeq (Required). 
#'  The name of the \code{\link[phyloseq]{phyloseq}} phyloseq object. This must have sample data with 
#'  Latitude and Longitude Columns.
#'  
#' @param size (Optional). Default \code{4}. 
#'  The size of the vertex points."Abundance" is a special code that will scale 
#'  points by the number of reads in a sample
#'  
#' @param region (Optional). Default \code{NULL}.
#'  The name of geographic region that can be used to zoom.
#'  The default worldmap cuts out Antartica. To get it back use region="world"
#' 
#' @param color (Optional). Default \code{NULL}.
#'  The name of the sample variable in \code{physeq} to use for color mapping
#'  of points (graph vertices).
#'  
#' @param shape (Optional). Default \code{NULL}.
#'  The name of the sample variable in \code{physeq} to use for shape mapping.
#'  of points (graph vertices).
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
#' @param projection (Optional). Default \code{NULL}. 
#'  Projection. Default of NULL will result in meractor projection. Non-default
#'  projection can be specified here but may require additional arguments 
#'  specified by the `parameter` and `orientation` arguments
#'
#' @param orientation (Optional). Default \code{NULL}. 
#'  Additional arguments for the map projection.
#'  
#' @param lat0 (Optional). Default \code{NULL}. 
#'  Additional arguments for nonstandard map projection.
#'  
#' @param lat1 (Optional). Default \code{NULL}. 
#'  Additional arguments for nonstandard map projection.
#'  
#' @param n (Optional). Default \code{NULL}. 
#'  Additional arguments for nonstandard map projection.
#'  
#' @param r (Optional). Default \code{NULL}. 
#'  Additional arguments for nonstandard map projection.
#'  
#' @param lon0 (Optional). Default \code{NULL}. 
#'  Additional arguments for nonstandard map projection.
#'  
#' @param seed (Optional). Default \code{1234}. 
#'  seed is used for repeatable randomness if you are using the jitter functions
#'  
#' @seealso 
#'  \href{http://zachcp.github.io/phylogeo/phylogeo_basics}{phylogeo basics}.
#'  \href{http://zachcp.github.io/phylogeo/phylogeo_projections}{phylogeo projections}.
#'
#' @import ggplot2  
#' @import maps
#' @import mapproj
#' @export
#' @examples 
#' map_phyloseq(batfecal, region="china", jitter=TRUE, jitter.x=2,jitter.y=2, color="PH")
#' map_phyloseq(batfecal, jitter=TRUE, seed=23454)
#' map_phyloseq(batfecal, projection="mollweide")
#' 
#' map_phyloseq(batfecal, projection="conic", lat0=15)
#' map_phyloseq(batfecal, projection="fisheye", n=0.5)
#' map_phyloseq(batfecal, projection="newyorker", r=0.3)
#' map_phyloseq(batfecal, projection="elliptic", lon0=10) 
#' map_phyloseq(batfecal, projection="albers", lat0=-20 , lat1=50)
#' map_phyloseq(batmicrobiome, jitter=TRUE, color="SCIENTIFIC_NAME")
map_phyloseq <- function(physeq, size=4, region=NULL, color=NULL, 
                         shape=NULL, alpha = 0.8, 
                         jitter=FALSE, jitter.x=3, jitter.y=3,
                         projection=NULL,orientation=NULL,
                         lat0=NULL, lat1=NULL, lon0=NULL,n=NULL, r=NULL,
                         seed=1234){
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
  
  #check plot options. "Abundance" is a special method for plotting by size
  .check_names(color,data)
  if(!size == "Abundance"){
    .check_names(size,data, allownumeric=T)
    print(size)
  }
  
  #create map
  #############################################################################
  worldmap <- .create_basemap(region=region, df=data, 
                              latcol=latcol,loncol=loncol,
                              projection=projection,orientation=orientation,
                              lat0=lat0, lat1=lat1, lon0=lon0,n=n, r=r)
  
  if(jitter){
    set.seed(seed)
    data <- .jitter_df(df=data,xcol=loncol,ycol=latcol,jitter.x=jitter.x,
                       jitter.y=jitter.y, seed=seed)
  }
  
  # how to hande when size information can be either global (outside of aes), 
  # per-sample (inside of aes)
  if(is.numeric(size)){
    worldmap <- worldmap + geom_point(data=data, 
                                      aes_string(x=loncol, y=latcol, 
                                                 group=NULL, color=color), 
                                      size = size, alpha= alpha) 
  }else if(size == "Abundance"){
    reads <- data.frame(sample_sums(physeq)); names(reads)<- "Abundance"
    data2 <- merge(data,reads,by="row.names")
    worldmap <- worldmap + geom_point(data=data2, 
                                      aes_string( x=loncol, y=latcol, group=NULL, 
                                                  color=color, size = "Abundance"),
                                      alpha= alpha) 
  }else{
    worldmap <- worldmap + geom_point(data=data, 
                             aes_string(x=loncol, y=latcol, group=NULL, 
                                        color=color, size = size),
                             alpha= alpha) 
  }
  
  return(worldmap)
}
################################################################################
#' Create a Network from the Phyloseq Objects and Draw A Map of the Clusters
#'
#' In this case, edges in the network are created if the distance between
#' nodes is below a potentially arbitrary threshold,
#' and special care should be given to considering the choice of this threshold.
#'   
#' @return a ggplot object
#' @param physeq (Required). 
#'  The name of the \code{\link[phyloseq]{phyloseq}} phyloseq object. This must have sample data with 
#'  Latitude and Longitude Columns.
#'  
#' @param igraph  (Optional). Default \code{NULL}
#'  An optional igraph object. Will reduce plotting time to use 
#'  a precalculated network 
#'  
#' @param region (Optional). Default \code{NULL}.
#'  The name of geographic region that can be used to zoom.
#' 
#' @param color (Optional). Default \code{NULL}.
#'  The name of the sample variable in \code{physeq} to use for color mapping
#'  of points (graph vertices). Note: "cluster" can be used to show igraph 
#'  clusters
#'  
#' @param shape (Optional). Default \code{NULL}.
#'  The name of the sample variable in \code{physeq} to use for shape mapping.
#'  of points (graph vertices).  Note: "cluster" can be used to show igraph 
#'  clusters
#'  
#' @param size (Optional). Default \code{4}. 
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
#'  Cutoff of the \code{distance} used to detmine whether a sample is 
#'  included in the network.
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
#' @param base_data_color (Optional). Default \code{grey}.
#'  named color to determine base data color
#'  
#' @param projection (Optional). Default \code{NULL}. 
#'  Projection. Default of NULL will result in meractor projection. Non-default
#'  projection can be specified here but may require additional arguments 
#'  specified by the `parameter` and `orientation` arguments
#'
#' @param orientation (Optional). Default \code{NULL}. 
#'  Additional arguments for the map projection.
#'  
#' @param lat0 (Optional). Default \code{NULL}. 
#'  Additional arguments for nonstandard map projection.
#'  
#' @param lat1 (Optional). Default \code{NULL}. 
#'  Additional arguments for nonstandard map projection.
#'  
#' @param n (Optional). Default \code{NULL}. 
#'  Additional arguments for nonstandard map projection.
#'  
#' @param r (Optional). Default \code{NULL}. 
#'  Additional arguments for nonstandard map projection.
#'  
#' @param lon0 (Optional). Default \code{NULL}. 
#'  Additional arguments for nonstandard map projection.
#'  
#' @param seed (Optional). Default \code{1234}. 
#'  seed is used for repeatable randomness if you are using the jitter functions
#'
#'  
#'
#' @seealso
#' \href{https://joey711.github.io/phyloseq/distance}{phyloseq's distance command}.
#' @seealso 
#'  \href{http://zachcp.github.io/phylogeo/phylogeo_basics}{phylogeo basics}.
#'  \href{http://zachcp.github.io/phylogeo/phylogeo_projections}{phylogeo projections}.
#' 
#' @import ggplot2
#' @import phyloseq
#' @import maps
#' @import mapproj
#' @importFrom igraph get.data.frame
#' @importFrom igraph get.vertex.attribute
#' @importFrom igraph clusters  
#' @export
#' @examples
#' map_network(batfecal)
#' map_network(batfecal, region="china", jitter=TRUE, lines=TRUE)
#' map_network(batfecal, region="china", jitter=TRUE, seed=3453)
#' map_network(batfecal, region="china", jitter=TRUE, lines=TRUE, maxdist=0.9)
#' map_network(batmicrobiome, lines=TRUE)
#' 
#' ig <- make_network(batmicrobiome)
#' map_network(batmicrobiome, igraph= ig, lines=TRUE)
#' map_network(batmicrobiome, igraph= ig, lines=TRUE, color="SCIENTIFIC_NAME")
#' map_network(batmicrobiome, igraph= ig, lines=TRUE, color="SCIENTIFIC_NAME", jitter=TRUE)
#' map_network(batmicrobiome, projection="mollweide", igraph= ig, lines=TRUE, color="SCIENTIFIC_NAME", jitter=TRUE)
map_network <- function(physeq, igraph=NULL, maxdist=0.9, distance="jaccard", 
                        color=NULL, region=NULL, size=4, alpha = 0.8, 
                        jitter=FALSE, jitter.x=3, jitter.y=3, shape=NULL, 
                        lines=FALSE, line_weight=1, line_color ="Black",
                        line_alpha=0.4 , base_data=FALSE, 
                        base_data_color="grey",projection=NULL, 
                        orientation=NULL,
                        lat0=NULL, lat1=NULL, lon0=NULL,n=NULL, r=NULL,
                        seed=1234){

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
  
  #add lines to ggplot object
  draw_lines <- function(plt, df2){
    df2 <- data.frame(df2) #to ensure list returns a df object
    plt <- plt + geom_line(data=df2,  
                           aes_string( x=loncol, 
                                       y=latcol, 
                                       group=names(df2)[1]))
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
  clusts <- seq( igraph::clusters(igraph)$no )
  clustdf <- Reduce( rbind, Map(get_clusters, clusts))
  mdf <- merge(clustdf, data.frame(data), by="row.names", all.x=TRUE)
  rownames(mdf) <- mdf$Row.names
  
  #check plot options
  .check_names(color,mdf)
  .check_names(size,mdf, allownumeric=TRUE)
  
  #create map
  ############################################
  worldmap <- .create_basemap(region=region, df=mdf,latcol=latcol, 
                              loncol=loncol, projection=projection, orientation=orientation,
                              lat0=lat0, lat1=lat1, lon0=lon0,n=n, r=r)
  
  #modify points if using jitter
  if(jitter){
    mdf <- .jitter_df(df=mdf, xcol=loncol, ycol=latcol, jitter.x=jitter.x,
                      jitter.y=jitter.y, seed=seed)
  }
  
  #add points that aren't part of a network
  if(base_data){
    network_points <- rownames(mdf)
    nonetworkdf <- data[!rownames(data) %in% network_points, ] 
    worldmap <- worldmap + geom_point(data = nonetworkdf, 
                                      aes_string(x=loncol,
                                                 y=latcol, 
                                                 group=NULL), 
                                      color=base_data_color, 
                                      size=size)
  }
  
  #addlines
  if(lines){
    linedf <- get_lines(df=mdf) 
    worldmap <- worldmap + 
                  geom_line(data=linedf,
                            aes_string(x=loncol,y=latcol, group="link"), 
                            size=line_weight, 
                            alpha=line_alpha, 
                            color=line_color)
  }
 
  # add points
  # how to hande when point_size information can be either global (outside of aes),
  # or per-sample (inside of aes)
  if(is.numeric(size)){
   points <- geom_point(data=mdf, size = size, alpha= alpha,
                        aes_string( x=loncol, y=latcol, group=NULL, 
                                    color=color, shape=shape)) 
  }else{
    points <- geom_point(data=mdf, alpha= alpha,
                         aes_string( x=loncol, y=latcol, group=NULL, 
                                     color=color, size = size, shape=shape)) 
  } 
  worldmap <- worldmap + points
  
  points <- geom_point(data=mdf, size = size, alpha= alpha,
                       aes_string(x=loncol, y=latcol, 
                                  group=NULL, color=NULL, shape=NULL)) 
  ###########################
  
  return(worldmap)
}
################################################################################
#' Map a Phyloseq Object while also drawing a phlogenetic tree of the taxa
#'
#' @return a ggplot object
#' 
#' @param physeq (Required). 
#'  The name of the \code{\link[phyloseq]{phyloseq}} phyloseq object. This must have sample data with 
#'  Latitude and Longitude Columns.
#'  
#'  @param region (Optional). Default \code{NULL}.
#'  The name of geographic region that can be used to zoom.
#'  The default worldmap cuts out Antartica. To get it back use region="world"
#'  
#' @param color (Optional). Default \code{NULL}.
#'  The name of the sample variable in \code{physeq} to use for color mapping
#'  of points (graph vertices).
#'  
#' @param size (Optional). Default \code{4}. 
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
#'  @param method (Optional). Degault \code{"sampledodge"}

#'  @param width_ratio (Optional). Default \code{2}.
#'  relative widths of tree and map
#'  
#'  @param map_on_left (Optional). Default \code{TRUE}.
#'  determine whether the map is on the left
#'  
#' @param nodelabf (Optional) Default \code{nodeplotblank}
#' @param treesize (Optional) Default \code{NULL}
#' @param min.abundance (Optional) Default \code{Inf}
#' @param label.tips (Optional) Default \code{NULL}
#' @param text.size (Optional) Default \code{NULL}
#' @param sizebase (Optional) Default \code{5}
#' @param base.spacing (Optional) Default \code{0.02}
#' @param ladderize (Optional) Default \code{TRUE}
#' @param plot.margin (Optional) Default \code{0.2}
#' @param title (Optional) Default \code{NULL}
#' @param treetheme (Optional) Default \code{NULL}
#' @param justify (Optional) Default \code{"jagged"}
#' whether to place the map or the tree on the left.
#' @param projection (Optional). Default \code{NULL}. 
#'  Projection. Default of NULL will result in meractor projection. Non-default
#'  projection can be specified here but may require additional arguments
#'  specified by the `parameter` and `orientation` arguments
#'
#' @param orientation (Optional). Default \code{NULL}. 
#'  Additional arguments for the map projection.
#'  
#' @param lat0 (Optional). Default \code{NULL}. 
#'  Additional arguments for nonstandard map projection.
#'  
#' @param lat1 (Optional). Default \code{NULL}. 
#'  Additional arguments for nonstandard map projection.
#'  
#' @param n (Optional). Default \code{NULL}. 
#'  Additional arguments for nonstandard map projection.
#'  
#' @param r (Optional). Default \code{NULL}. 
#'  Additional arguments for nonstandard map projection.
#'  
#' @param lon0 (Optional). Default \code{NULL}. 
#'  Additional arguments for nonstandard map projection.
#'  
#' @seealso \code{\link[phyloseq]{plot_tree}}
#' @seealso \code{\link[ggplot2]{map_data}}
#' @seealso 
#'  \href{http://zachcp.github.io/phylogeo/phylogeo_basics}{phylogeo basics}.
#'  \href{http://zachcp.github.io/phylogeo/phylogeo_projections}{phylogeo projections}.
#'
#' @import ggplot2
#' @import maps
#' @import mapproj
#' @export
#' @examples
#' map_tree(epoxamicin_KS)
#' map_tree(epoxamicin_KS, color="Geotype", jitter=TRUE)
#' map_tree(epoxamicin_KS, projection="gilbert", color="Geotype", jitter=TRUE)
map_tree <- function(physeq,  region=NULL, color = NULL,size=4, alpha=0.8,
                    jitter= FALSE, jitter.x=3, jitter.y=3, 
                    method = "sampledodge", nodelabf = nodeplotblank, 
                    treesize = NULL, min.abundance = Inf, label.tips = NULL,
                    text.size = NULL, sizebase = 5, base.spacing = 0.02, 
                    ladderize = TRUE,plot.margin = 0.2, title = NULL, 
                    treetheme = NULL, justify = "jagged",width_ratio = 2, 
                    map_on_left = FALSE, projection=NULL, orientation=NULL,
                    lat0=NULL, lat1=NULL, lon0=NULL,n=NULL, r=NULL) {
    #check for the existence of a tree: lifted from phyloseq's plot_tree
    if(!"phy_tree" %in% phyloseq::getslots.phyloseq(physeq)){
      stop("tree missing or invalid. map-tree requires a phylogenetic tree")
    }
    #trim samples that are not in the tree
    physeq2 <- phyloseq::prune_samples(phyloseq::sample_sums(physeq) > 0, physeq)
    
    mapplot  <- map_phyloseq(physeq2, region=region, color= color, size=size, 
                             alpha = alpha, jitter=jitter, jitter.x=jitter.x, 
                             jitter.y=jitter.y,projection=projection,orientation=orientation,
                             lat0=lat0, lat1=lat1, lon0=lon0,n=n, r=r)  + 
                    theme(legend.position="none") 
    
    treeplot <- phyloseq::plot_tree(physeq2, color=color, label.tips=label.tips,
                                    text.size=text.size, sizebase=sizebase, 
                                    base.spacing = base.spacing, 
                                    ladderize = ladderize, 
                                    plot.margin = plot.margin, 
                                    title = title, 
                                    treetheme=treetheme, 
                                    justify = justify, 
                                    nodelab =nodelabf) +
      theme(legend.key = element_rect(fill = "white")) +
      scale_y_continuous(expand = c(0,0)) + 
      scale_x_continuous(expand = c(0,0))
    # # trim space by setting xlims
    # xvals <- treeplot$data$x
    # xvals <- xvals[!is.na(xvals)]
    # xmin <- min(xvals)
    # xmax <- max(xvals) * 1.5
    # treeplot <- treeplot + xlim( min(xvals), max(xvals))
    
    if(map_on_left){
        combinedplot <- gridExtra::arrangeGrob(mapplot + 
                                               theme(legend.position="none"),
                                               treeplot, 
                                               ncol=2, 
                                               widths=c(width_ratio,1))
    } else{
        combinedplot <- gridExtra::arrangeGrob(treeplot + 
                                                 theme(legend.position="none"),
                                               mapplot, 
                                               ncol=2, 
                                               widths=c(1,width_ratio))    
    }
    return(combinedplot)
}
################################################################################
#' Explore the spatial distribution of subsets of your sequence data 
#'   
#' @param physeq (Required). 
#'  The name of the phyloseq object. This must have sample data with 
#'  Latitude and Longitude Columns.
#'  
#' @param clusternum (Optional). Default \code{3}.
#'  Number of kmeans clusters to divide your phylogenetic tree into
#'  
#' @import ggplot2
#' @import maps
#' @import mapproj
#' @import gridExtra
#' @seealso 
#'  \href{http://zachcp.github.io/phylogeo/phylogeo_basics}{phylogeo basics}.
#' @export
#' @examples
#' map_clusters(epoxamicin_KS, clusternum=6)
#' map_clusters(epoxamicin_KS, clusternum=10)
map_clusters <- function(physeq, clusternum=3){
  # check for the existence of a tree: lifted from phyloseq's plot_tree
  if(!"phy_tree" %in% phyloseq::getslots.phyloseq(physeq)){
    stop("tree missing or invalid. map-tree requires a phylogenetic tree")
  }
  
  
  ################################################
  # Helper Functions
  ################################################
  
  # get a list of otus in a cluster
  otus_in_a_cluster <- function(clusternum, clusterdf=clusters){
    df <- clusterdf[clusterdf$cluster == clusternum, ]
    return(rownames(df))
  }
  
  #add the cluster information to the taxonomy table in a "Cluster" column
  markcluster <- function(x, otus){
    if(is.na(x) ){ 0
    }else{ ifelse(x %in% otus, 1, 0)}
  }
  
  # add cluster information to the taxonomy tree and return the modified physeq
  add_cluster_to_taxtree <- function(physeq, otus){
    tax <- data.frame(tax_table(physeq))
    tax$Cluster <- sapply(row.names(tax),markcluster, otus=otus)
    tax_table(physeq) <- as.matrix(tax)
    return(physeq)
  }
  
  # creates a tree where the deignated cluster is colored red
  make_cluster_tree <- function(physeq, clusternum){
    otulist <- otus_in_a_cluster(clusternum)
    physeq2 <- add_cluster_to_taxtree(physeq, otulist)
    p <-  plot_tree(physeq2, color="Cluster", ladderize=TRUE) +
      scale_color_manual(values = c("black","red")) +
      theme(legend.position="none") +
      scale_y_continuous(expand = c(0,0)) + 
      scale_x_continuous(expand = c(0,0))
    return(p)
  }
  
  make_map_of_cluster <- function(physeq, clusternum){
    otulist <- otus_in_a_cluster(clusternum)
    physeq2 <- add_cluster_to_taxtree(physeq, otulist)
    physeq2 <- subset_taxa(physeq2, Cluster=="1")
    physeq2 <- prune_samples(sample_sums(physeq2)>0, physeq2)
    p <-  map_phyloseq(physeq2, size="Abundance")
    return(p)

  }
  
  # make biplot of the tow together
  make_diplot <- function(clusternum, physeq=physeq){
     p1 <- make_cluster_tree(physeq, clusternum)
     p2 <- make_map_of_cluster(physeq, clusternum)
     combinedplot <- gridExtra::arrangeGrob(p1,p2, ncol=2, widths=c(1,2))
  }

  ################################################
  # Do the work
  ################################################
  
  # get kmeans data from the phylogenetic tree
  distances <- ape::cophenetic.phylo( phy_tree(physeq))
  kmeans_out <- kmeans(distances, centers=clusternum)
  clusters <- data.frame(kmeans_out$cluster)
  names(clusters) <- "cluster"
  clusters$clusterOTU <- row.names(clusters)
  
  # plot everthing out
  plots <- lapply(1:clusternum, make_diplot, physeq=physeq)
  combinedplot <- do.call(gridExtra::arrangeGrob,plots)
  
  return(combinedplot)
}
