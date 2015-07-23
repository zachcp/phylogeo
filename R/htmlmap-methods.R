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
#'  The name of the \code{\link[phyloseq]{phyloseq}} phyloseq object.
#'  This must have sample data with Latitude and Longitude Columns.
#'
#' @param size (Optional). Default \code{5}.
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
#' @import phyloseq
#' @import leaflet
#' @import dplyr
#' @importFrom magrittr %<>%
#' @export
#' @examples
#' data(mountainsoil)
#' htmlmap_phyloseq(mountainsoil, size=3)
#' data(batmicrobiome)
#' htmlmap_phyloseq(batmicrobiome, color="blue")
htmlmap_phyloseq <- function(physeq, size = 5, color = "blue"){
  #get data
  data = sample_data(physeq) %>%
      data.frame() %>%
      add_rownames(var = "samplename")

  #customize circle size
  if (is.numeric(size)) {
      data$circlesize <- size
  } else if (size == "Abundance") {
      data$circlesize <- phyloseq::sample_sums(physeq)
  } else {
      stop("Size can be a numeric value or the word 'Abundance' ")
  }

  #basemap
  map = leaflet(data) %>%
      addTiles() %>%
      addCircleMarkers(radius = ~circlesize,
                       color = makecolors(data,color),
                       popup = ~samplename)

  return(map)
}

################################################################################
#' Create a Distance Network from Phyloseq Data and Draw An HTML Map of it.
#'
#' In this case, edges in the network are created if the distance between
#' nodes is below a potentially arbitrary threshold,
#' and special care should be given to considering the choice of this threshold.
#'
#' @return an  \code{\link[htmlwidgets]{htmlwidgets}} plot
#' @seealso http://rstudio.github.io/leaflet/
#'
#' @param physeq (Required).
#'  The name of the \code{\link[phyloseq]{phyloseq}} phyloseq object. This must
#'  have sample data with Latitude and Longitude Columns.
#'
#' @param igraph  (Optional). Default \code{NULL}
#'  An optional \code{\link[igraph]{igraph}} igraph object. Will reduce plotting
#'  time to use a precalculated network
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
#' @param size (Optional). Default \code{5}.
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
#' @import leaflet
#' @import dplyr
#' @import leaflet
#' @importFrom igraph get.data.frame
#' @importFrom igraph get.vertex.attribute
#' @importFrom igraph clusters
#' @export
#' @examples
#' htmlmap_network(mountainsoil)
#' htmlmap_network(mountainsoil, maxdist=0.9)
#'
#' htmlmap_network(batmicrobiome, maxdist=0.5)
#' ig <- make_network(batmicrobiome)
#' htmlmap_network(batmicrobiome, igraph= ig)
#' htmlmap_network(epoxamicin_KS, maxdist=0.99, line_color = "red",
#'                 line_weight = 4, line_alpha=0.5)
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
                            size=5,
                            rescale=TRUE){
    # Calculate Distance
    if (inherits(distance, "dist")) {
        # If distance a distance object, use it rather than re-calculate
        Distance <- distance
        # Check that it at least has (a subset of) the correct labels
        if (!all(attributes(distance)$Labels %in% sample_names(physeq))) {
            stop("Some or all `distance` index labels do not match sample names in `physeq`")
            }
    } else {
        # Coerce to character and attempt distance calculation
        scaled_distance = function(physeq, method, type, rescale){
            Dist = distance(physeq, method, type)
            if (rescale) {
                # rescale the distance matrix to be [0, 1]
                Dist <- Dist / max(Dist, na.rm = TRUE)
                Dist <- Dist - min(Dist, na.rm = TRUE)
            }
            return(Dist)
        }
        distance <- as(distance[1], "character")
        Distance = scaled_distance(physeq, distance, type="samples", rescale=rescale)
    }

    #check basic physeq and lat/lon and make clusters
    physeqdata <- check_phyloseq(physeq)
    data <- physeqdata$sampledata %>% add_rownames(var = "samplename")
    rownames(data) <- data$samplename

    #customize circle size
    if (is.numeric(size)) {
        data$circlesize <- size
    } else if (size == "Abundance") {
            data$circlesize <- phyloseq::sample_sums(physeq)
    } else {
          stop("Size can be a numeric value or the word 'Abundance' ")
    }

    # convert distances to lines
    distdf = dist_to_edge_table(Distance, maxdist) %>%
        edgetable_to_linedf(physeqdata = physeqdata)

    #create base map
    ############################################
    map <- leaflet(data) %>% addTiles()

    # add lines to map
    for (g in unique(distdf$rowname)) {
        sdf <- distdf[distdf$rowname == g, ]
        map <- map %>%
            addPolylines(data = sdf,
                         lng = ~lng,
                         lat = ~lat,
                         weight = ~distance)
      }

    #add points to map
    map <- map %>% addCircleMarkers(radius = ~circlesize,
                                    color = makecolors(data, color),
                                    opacity = circle_alpha,
                                    fillOpacity = fillOpacity,
                                    popup = ~samplename)
    return(map)
}

#' makecolors
#'
#' handles the color values and passes correct values to leaflet
#' @param data
#' @param color
#'
#' @import leaflet
#'
#' @return a color string or a vector of color strings
#' @keywords internal
#' https://github.com/rstudio/leaflet/issues/80
makecolors <- function(data, color){
  columns <- names(data)

  #test if the string is a column value
  if (!color %in% columns) {
    return(color)
  } else {
    testdata <- data[[color]]
    #get colors depending on the columntypes
    if (is.factor(testdata)) {
      return(leaflet::colorFactor("RdYlBu", NULL)(data[[color]]))
    } else if (is.numeric(testdata)) {
      return(leaflet::colorNumeric("Blues", NULL)(data[[color]]))
    } else if (is.character(testdata)) {
       fac = as.factor(testdata)
       return(leaflet::colorFactor("RdYlBu", NULL)(fac))
    } else {
      stop("Error in makecolor function. Could not ascertain column datatype")
    }
  }
}
