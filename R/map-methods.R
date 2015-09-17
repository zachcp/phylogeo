#
# methods for drawing maps of phyloseq datasets
#
###############################################################################
#' Simple Maps for \code{\link[phyloseq]{phyloseq}} Samples
#'
#' map_phyloseq is a plotting function that draws a map of your microbiome
#' dataset. This function acts on \code{\link[phyloseq]{phyloseq}}
#' objects and requires that the \code{\link[phyloseq]{sample_data}}
#' sample_data table contains Latitude and Longitude columns. The resulting
#' map can be customized by a nubmer of parameters including size="Abundance"
#' which will scale the size of hte circle according to the nubmer of reads in
#' the \code{\link[phyloseq]{otu_table}} otu table.
#'
#' This plotting function will draw a map of the samples in your microbiome
#' project, using the \code{\link[phyloseq]{phyloseq}} phyloseq package. Most
#' aspect of the map are customizable using the parameters below.
#'
#'
#' @return a ggplot object
#' @param physeq (Required).
#'  The name of the \code{\link[phyloseq]{phyloseq}} object.
#'  This must have sample data with Latitude and Longitude Columns.
#'
#' @param mapdata (Optional). Default \code{"world"}
#'  The name of the \code{\link[maps]{maps}} or \code{\link[mapdata]{mapdata}}
#'  to use for drawing.
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
#' @param projection (Optional). Default \code{"mercator"}.
#' See the documentation of \code{\link[mapproj]{mapproject}} for the list
#' of available projections and their descriptions.
#'
#' @param ...  (Optional Arguments). Default \code{NULL}.
#' Arguments passed internally to the \code{\link[mapproj]{mapproject}}
#' function from the \pkg{mapproj} package to control the projection including
#' "orientation, lat0, lat1, n, r, lon0" See the documentation of
#' \code{\link[mapproj]{mapproject}} for more information
#' (?\code{\link[mapproj]{mapproject}}).
#'
#' @param seed (Optional). Default \code{1234}.
#'  seed is used for repeatable randomness if you are using the jitter functions
#'
#' @seealso
#'  \href{http://zachcp.github.io/phylogeo/phylogeo_basics}{phylogeo basics}.
#'  \href{http://zachcp.github.io/phylogeo/phylogeo_projections}{phylogeo
#'  projections}.
#'
#'@seealso
#'  \code{\link[mapproj]{mapproject}}
#'  \code{\link[ggplot2]{coord_map}}
#'  \code{\link[phyloseq]{phyloseq}}
#'
#'
#' @import ggplot2
#' @import maps
#' @import mapproj
#' @export
#' @examples
#' map_phyloseq(mountainsoil, region="china", jitter=TRUE, jitter.x=2,jitter.y=2, color="PH")
#' map_phyloseq(mountainsoil, jitter=TRUE, seed=23454)
#' map_phyloseq(mountainsoil, projection="mollweide")
#'
#' map_phyloseq(mountainsoil, projection="conic", lat0=15)
#' map_phyloseq(mountainsoil, projection="fisheye", n=0.5)
#' map_phyloseq(mountainsoil, projection="newyorker", r=0.3)
#' map_phyloseq(mountainsoil, projection="elliptic", lon0=10)
#' map_phyloseq(mountainsoil, projection="albers", lat0=-20 , lat1=50)
#' map_phyloseq(batmicrobiome, jitter=TRUE, color="SCIENTIFIC_NAME")
map_phyloseq <- function(physeq,
                         mapdata="world",
                         size=4,
                         region=NULL,
                         color=NULL,
                         shape=NULL,
                         alpha = 0.8,
                         jitter=FALSE,
                         jitter.x=3,
                         jitter.y=3,
                         projection="mercator",
                         orientation=NULL,
                         seed=1234,
                         ...) {
  #convert to phylogeo
  phygeo <- phylogeo(physeq)

  #check plot options. "Abundance" is a special method for plotting by size
  check_names(color, sample_data(phygeo))

  if (!size == "Abundance") check_names(size,data, allownumeric = TRUE)

  #create map
  #############################################################################
  worldmap <- create_basemap(mapdata = mapdata,
                             region = region,
                             df = sample_data(phygeo),
                             latcol = phygeo@latitude,
                             loncol = phygeo@longitude,
                             projection = projection,
                             orientation = orientation,...)

  if (jitter) {
    set.seed(seed)
    data <- jitter_df(df = sample_data(phygeo),
                      xcol = phygeo@longitude,
                      ycol = phygeo@latitude,
                      jitter.x = jitter.x,
                      jitter.y = jitter.y,
                      seed = seed)
  }

  # how to hande when size information can be either global (outside of aes),
  # per-sample (inside of aes)
  if (is.numeric(size)) {
    worldmap <- worldmap +
        geom_point(data = sample_data(phygeo),
                   aes_string(x = phygeo@longitude,
                       y = phygeo@latitude,
                       group = NULL,
                       color = color),
                   size = size,
                   alpha = alpha)
  } else if (size == "Abundance") {
    reads <- data.frame(sample_sums(physeq))
    names(reads) <- "Abundance"
    data2 <- merge(sample_data(phygeo), reads, by = "row.names")
    worldmap <- worldmap +
        geom_point(data = data2,
                   aes_string(x = phygeo@longitude,
                              y = phygeo@latitude,
                              group = NULL,
                              color = color,
                              size = "Abundance"),
                   alpha = alpha)
  }else{
    worldmap <- worldmap +
        geom_point(data = sample_data(phygeo),
                   aes_string(x = phygeo@longitude,
                              y = phygeo@latitude,
                              group = NULL,
                              color = color,
                              size = size),
                   alpha = alpha)
  }

  return(worldmap)
}
################################################################################
#' Mapping Functions to Visualize Phyloseq Datasets Sample Similarity.
#'
#' map_network is a plotting function that draws a map of your microbiome
#' dataset, highlighting the degree of similarity between sample sites.
#' This funciton acts on \code{\link[phyloseq]{phyloseq}}
#' datasets and requires that the \code{\link[phyloseq]{sample_data}}
#' table contains Latitude and Longitude columns. map_network
#' will calculate the ecological simiarity of each of your sample sites using
#' the distance metrics available in the \code{\link[phyloseq]{make_network}}
#' phyloseq package. These distances are then plotted as lines. Many aspects
#' of the plots can be customized using the parameters below.
#'
#' @return a ggplot object
#' @param physeq (Required).
#'  The name of the \code{\link[phyloseq]{phyloseq}} object.
#'  This must have sample data with Latitude and Longitude Columns.
#'
#' @param mapdata (Optional). Default \code{"world"}
#'  The name of the \code{\link[maps]{maps}} or \code{\link[mapdata]{mapdata}}
#'  to use for drawing.
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
#'  Boolean to determine whether to include data points that aren't in a network.
#'
#' @param base_data_color (Optional). Default \code{grey}.
#'  named color to determine base data color
#'
#' @param projection (Optional). Default \code{"mercator"}.
#' See the documentation of \code{\link[mapproj]{mapproject}} for the list
#' of available projections and their descriptions.
#'
#' @param ...  (Optional Arguments). Default \code{NULL}.
#' Arguments passed internally to the \code{\link[mapproj]{mapproject}}
#' function from the \pkg{mapproj} package to control the projection including
#' "orientation, lat0, lat1, n, r, lon0" See the documentation of
#' \code{\link[mapproj]{mapproject}} for more information
#' (?\code{\link[mapproj]{mapproject}}).
#'
#' @param seed (Optional). Default \code{1234}.
#'  seed is used for repeatable randomness if you are using the jitter functions
#'
#' @seealso
#' \href{https://joey711.github.io/phyloseq/distance}{phyloseq's distance command}.
#' @seealso
#'  \href{http://zachcp.github.io/phylogeo/phylogeo_basics}{phylogeo basics}.
#'  \href{http://zachcp.github.io/phylogeo/phylogeo_projections}{phylogeo projections}.
#'@seealso
#'  \code{\link[mapproj]{mapproject}}
#'  \code{\link[ggplot2]{coord_map}}
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
#' map_network(mountainsoil)
#' map_network(mountainsoil, region="china", jitter=TRUE, lines=TRUE)
#' map_network(mountainsoil, region="china", jitter=TRUE, seed=3453)
#' map_network(mountainsoil, region="china", jitter=TRUE, lines=TRUE, maxdist=0.9)
#' map_network(batmicrobiome, lines=TRUE)
#'
#' ig <- make_network(batmicrobiome)
#' map_network(batmicrobiome, igraph= ig, lines=TRUE)
#' map_network(batmicrobiome, igraph= ig, lines=TRUE, color="SCIENTIFIC_NAME")
#' map_network(batmicrobiome, igraph= ig, lines=TRUE, color="SCIENTIFIC_NAME",
#'            jitter=TRUE)
#' map_network(batmicrobiome, projection="mollweide", igraph= ig,
#'             lines=TRUE, color="SCIENTIFIC_NAME", jitter=TRUE)
map_network <- function(physeq,
                        mapdata="world",
                        igraph=NULL,
                        maxdist=0.9,
                        distance="jaccard",
                        color=NULL,
                        region=NULL,
                        size=4,
                        alpha = 0.8,
                        jitter=FALSE,
                        jitter.x=3,
                        jitter.y=3,
                        shape=NULL,
                        lines=FALSE,
                        line_weight=1,
                        line_color ="Black",
                        line_alpha=0.4,
                        base_data=FALSE,
                        base_data_color="grey",
                        projection="mercator",
                        orientation=NULL,
                        ...,
                        seed=1234) {

  #helper functions to calculate membership in clusters or lines
  ##############################################################################
  get_clusters <- function(num, graph=igraph){
    # get cluster membership info from igraph object from cluster
    # with clusterid of 'num'
    clusts  <- clusters(graph)
    members <- which(clusts$membership == num) #get membership
    names   <- get.vertex.attribute(graph, 'name', members)
    df = data.frame(names)
    df['cluster'] <- as.character(num)
    rownames(df) <- df$names
    df    #return a df with name/cluster columns
  }

  get_lines <- function(graph=igraph, df=data){
    # get each edge of the network and return a list of dataframes with
    # the node info

    getline_df <- function(i, l = links, df1 = df){
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

    #####################################

  #convert to phylogeo
  phygeo <- phylogeo(physeq)

  # make network, get cluster information, and add thamesat to the
  # original dataframe.
  if (is.null(igraph)) {
    igraph <- make_network(physeq, max.dist = maxdist, distance = distance)
  } else {
    if (!"igraph" %in% class(igraph)) {
      stop("igraph must be an igraph network object")}
  }
  clusts <- seq(igraph::clusters(igraph)$no)
  clustdf <- Reduce( rbind, Map(get_clusters, clusts))
  mdf <- merge(clustdf, sample_data(phygeo), by = "row.names", all.x = TRUE)
  rownames(mdf) <- mdf$Row.names

  #check plot options
  check_names(color, mdf)
  check_names(size, mdf, allownumeric = TRUE)

  #create map
  ############################################
  worldmap <- create_basemap(mapdata = mapdata,
                             region = region,
                             df = mdf,
                             latcol = phygeo@latitude,
                             loncol = phygeo@longitude,
                             projection = projection,
                             orientation = orientation,
                             ...)

  #modify points if using jitter
  if (jitter) {
    mdf <- jitter_df(df = mdf,
                     xcol = phygeo@longitude,
                     ycol = phygeo@latitude,
                     jitter.x = jitter.x,
                     jitter.y = jitter.y,
                     seed = seed)
  }

  #add points that aren't part of a network
  if (base_data) {
    network_points <- rownames(mdf)
    nonetworkdf <- data[!rownames(data) %in% network_points, ]
    worldmap <- worldmap +
        geom_point(data = nonetworkdf,
                   aes_string(x = phygeo@longitude,
                              y = phygeo@latitude,
                              group = NULL),
                   color = base_data_color,
                   size = size)
  }

  #addlines
  if (lines) {
    linedf <- get_lines(df = mdf)
    worldmap <- worldmap +
      geom_line(data = linedf,
                aes_string(x = phygeo@longitude,
                           y = phygeo@latitude,
                           group = "link"),
                size = line_weight,
                alpha = line_alpha,
                color = line_color)
  }

  # add points
  # how to hande when point_size information can be either global
  # (outside of aes), or per-sample (inside of aes)
  if (is.numeric(size)) {
    points <- geom_point(data = mdf,
                         size = size,
                         alpha = alpha,
                         aes_string(x = phygeo@longitude,
                                    y = phygeo@latitude,
                                    group = NULL,
                                    color = color,
                                    shape = shape))
  } else {
    points <- geom_point(data = mdf,
                         alpha = alpha,
                         aes_string(x = phygeo@longitude,
                                    y = phygeo@latitude,
                                    group = NULL,
                                    color = color,
                                    size = size,
                                    shape = shape))
  }
  worldmap <- worldmap + points
  return(worldmap)
}
################################################################################
#' Plot a phlogenetic tree and map the location of the tree's end-nodes.
#'
#' map_tree is a plotting function that acts on \code{\link[phyloseq]{phyloseq}}
#' datasets to plot a figure that is a composite of a phylogetic tree
#' and a geographic map. The tips on the tree and the locations of the map
#' correspond to each other, allowing for an easy way to look for the geographic
#' distrivution of clades in your tree.
#'
#' @return a ggplot object
#'
#' @param physeq (Required).
#'  The name of the \code{\link[phyloseq]{phyloseq}} dataset.
#'  This must have sample data with Latitude and Longitude Columns.
#'
#'  @param mapdata (Optional). Default \code{"world"}
#'  The name of the \code{\link[maps]{maps}} or \code{\link[mapdata]{mapdata}}
#'  to use for drawing.
#'
#'  @param region (Optional). Default \code{NULL}.
#'  The name of geographic region that can be used to zoom.The default worldmap
#'  cuts out Antartica. To get it back use region="world"
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
#'
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
#' @param projection (Optional). Default \code{"mercator"}.
#' See the documentation of \code{\link[mapproj]{mapproject}} for the list
#' of available projections and their descriptions.
#' @param ...  (Optional Arguments). Default \code{NULL}.
#' Arguments passed internally to the \code{\link[mapproj]{mapproject}}
#' function from the \pkg{mapproj} package to control the projection including
#' "orientation, lat0, lat1, n, r, lon0" See the documentation of
#' \code{\link[mapproj]{mapproject}} for more information
#' (?\code{\link[mapproj]{mapproject}}).
#'
#' @seealso
#'  \href{http://zachcp.github.io/phylogeo/phylogeo_basics}{phylogeo basics}.
#'  \href{http://zachcp.github.io/phylogeo/phylogeo_projections}{phylogeo projections}.
#' @seealso
#'  \code{\link[mapproj]{mapproject}}
#'  \code{\link[ggplot2]{coord_map}}
#' @seealso
#'   \code{\link[phyloseq]{plot_tree}}
#'
#' @import ggplot2
#' @import phyloseq
#' @import maps
#' @import mapproj
#' @importFrom ape ladderize
#' @importFrom ape cophenetic.phylo
#' @importFrom ape reorder.phylo
#' @importFrom cowplot plot_grid
#' @export
#' @examples
#' map_tree(epoxomicin_KS)
#' map_tree(epoxomicin_KS, color="Geotype", jitter=TRUE)
#' map_tree(epoxomicin_KS, projection="gilbert", color="Geotype", jitter=TRUE)
map_tree <- function(physeq,
                     mapdata="world",
                     region=NULL,
                     color = NULL,
                     size=4,
                     alpha=0.8,
                     jitter= FALSE,
                     jitter.x=3,
                     jitter.y=3,
                     method = "sampledodge",
                     nodelabf = nodeplotblank,
                     treesize = NULL,
                     min.abundance = Inf,
                     label.tips = NULL,
                     text.size = NULL,
                     sizebase = 5,
                     base.spacing = 0.02,
                     ladderize = TRUE,
                     plot.margin = 0.2,
                     title = NULL,
                     treetheme = NULL,
                     justify = "jagged",
                     width_ratio = 2,
                     map_on_left = FALSE,
                     projection="mercator",
                     orientation=NULL, ...) {
  #check for the existence of a tree: lifted from phyloseq's plot_tree
  if (!"phy_tree" %in% getslots.phyloseq(physeq)) {
    stop("tree missing or invalid. map-tree requires a phylogenetic tree")
  }
  #trim samples that are not in the tree
  physeq2 <- prune_samples(sample_sums(physeq) > 0, physeq)

  mapplot <- map_phyloseq(physeq2,
                          mapdata = mapdata,
                          region = region,
                          color = color,
                          size = size,
                          alpha = alpha,
                          jitter = jitter,
                          jitter.x = jitter.x,
                          jitter.y = jitter.y,
                          projection = projection,
                          orientation = orientation, ...)  +
    theme(legend.position = "none") +
    labs(x = NULL, y = NULL)

  treeplot <- plot_tree(physeq2,
                        color = color,
                        label.tips = label.tips,
                        text.size = text.size,
                        sizebase = sizebase,
                        base.spacing = base.spacing,
                        ladderize = ladderize,
                        plot.margin = plot.margin,
                        title = title,
                        treetheme = treetheme,
                        justify = justify,
                        nodelabf = nodelabf) +
    theme(legend.key = element_rect(fill = "white"),
          axis.line = element_blank()) +
    labs(x = NULL, y = NULL)

  # # trim space by setting xlims
  # xvals <- treeplot$data$x
  # xvals <- xvals[!is.na(xvals)]
  # xmin <- min(xvals)
  # xmax <- max(xvals) * 1.5
  # treeplot <- treeplot + xlim( min(xvals), max(xvals))
  if (map_on_left) {
    combinedplot <- plot_grid(mapplot + theme(legend.position = "none"), treeplot)
  } else{
    combinedplot <- plot_grid(treeplot,mapplot + theme(legend.position = "none"))
  }
  return(combinedplot)
}
###############################################################################
#' Explore the spatial distribution of clustered-subsets of your sequence data
#'
#' map_clusters is a plotting function for  \code{\link[phyloseq]{phyloseq}}
#' datasets which requires a phylogentic tree. map_clusters will use
#' the phylogenetic tree (required) of your phyloseq dataset and will cluster
#' these sequence into a specified number of clusters using kmeans clustering.
#' A map is then generated for each of these clusters showing the cluster and
#' locaiton of the sequences in each cluster on a map.
#'
#' @return a ggplot object
#' @param physeq (Required).
#'  The name of the \code{\link[phyloseq]{phyloseq}} object.
#'  This must have sample data with Latitude and Longitude Columns.
#'
#' @param clusternum (Optional). Default \code{3}.
#'  Number of clusters to divide your phylogenetic tree into using
#'  kmeans clustering.
#'
#'  @param relwidth (Optional). Default \code{2}
#'  The relative width of the map compaed withthe tree (per cluster.)
#'
#' @import ggplot2
#' @import phyloseq
#' @import maps
#' @import mapproj
#' @importFrom ape ladderize
#' @importFrom ape cophenetic.phylo
#' @importFrom ape reorder.phylo
#' @importFrom cowplot plot_grid
#' @seealso
#'  \href{http://zachcp.github.io/phylogeo/phylogeo_basics}{phylogeo basics}.
#' @seealso
#'  \code{\link[mapproj]{mapproject}}
#'  \code{\link[ggplot2]{coord_map}}
#'
#' @export
#' @examples
#' map_clusters(epoxomicin_KS, clusternum=2)
#' map_clusters(epoxomicin_KS, clusternum=6)
map_clusters <- function(physeq, clusternum=3, relwidth=2){
  # check for the existence of a tree: lifted from phyloseq's plot_tree
  if (!"phy_tree" %in% getslots.phyloseq(physeq)) {
    stop("tree missing or invalid. map-tree requires a phylogenetic tree")
  }

  ################################################
  # Helper Functions
  ################################################

  # get a list of otus in a cluster
  otus_in_a_cluster <- function(clusternum, clusterdf){
    df <- clusterdf[clusterdf$cluster == clusternum, ]
    return(rownames(df))
  }

  #add the cluster information to the taxonomy table in a "Cluster" column
  markcluster <- function(x, otus){ ifelse(is.na(x),0,ifelse(x %in% otus, 1, 0))}

  # add cluster information to the taxonomy tree and return the modified physeq
  add_cluster_to_taxtree <- function(physeq, otus){
    tax <- data.frame(tax_table(physeq))
    tax$Cluster <- sapply(row.names(tax),markcluster, otus = otus)
    tax_table(physeq) <- as.matrix(tax)
    return(physeq)
  }

  # creates a tree where the designated cluster is colored red
  make_cluster_tree <- function(clusternum, clusterdf, physeq){
    otulist <- otus_in_a_cluster(clusternum, clusterdf = clusterdf)
    physeq2 <- add_cluster_to_taxtree(physeq, otulist)
    p <-  plot_tree(physeq2,
                    color = "Cluster",
                    ladderize = TRUE) +
      scale_color_manual(values = c("black","red")) +
      theme(legend.position = "none",
            axis.line = element_blank()) +
      labs(x = NULL,y = NULL)
    return(p)
  }

  make_map_of_cluster <- function(clusternum,clusterdf,physeq=physeq){
    otulist <- otus_in_a_cluster(clusternum = clusternum, clusterdf = clusterdf)
    physeq2 <- add_cluster_to_taxtree(physeq, otulist)
    physeq2 <- subset_taxa(physeq2, Cluster == "1")
    physeq2 <- prune_samples(sample_sums(physeq2) > 0, physeq2)
    p <-  map_phyloseq(physeq2, size = "Abundance")
    return(p)
  }

  ################################################
  # Do the work
  ################################################

  # get kmeans data from the phylogenetic tree
  distances <- cophenetic.phylo(phy_tree(physeq))
  kmeans_out <- kmeans(distances, centers = clusternum)
  clusters <- data.frame(kmeans_out$cluster)
  names(clusters) <- "cluster"
  clusters$clusterOTU <- row.names(clusters)

  #make both plots fro each cluster
  plots <- list()
  for (i in 1:clusternum) {
      plots <- c(plots,
                 list(make_cluster_tree(clusternum = i,
                                        clusterdf = clusters,
                                        physeq = physeq)),
                 list(make_map_of_cluster(clusternum = i,
                                          clusterdf = clusters,
                                          physeq = physeq)))
  }
  return(plot_grid(plotlist = plots,
                   ncol = 2,
                   rel_widths = c(1,relwidth)))
}
