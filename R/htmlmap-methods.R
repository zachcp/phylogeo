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
htmlmap_phyloseq <- function(physeq, size = 5, color = "blue", map = NULL) {
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
  if (is.null(map))  map <- leaflet() %>% addTiles()


  map <- map %>%
      addCircleMarkers(data = data,
                       radius = ~circlesize,
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
#' htmlmap_network(batmicrobiome, maxdist=0.5)
#' htmlmap_network(epoxamicin_KS, maxdist=0.99, line_color = "red",
#'                 line_weight = 4, line_alpha=0.5)
htmlmap_network <- function(physeq,
                            map = NULL,
                            #distance related
                            maxdist=0.9,
                            distance="jaccard",
                            #linerelated
                            line_color ="black",
                            line_alpha=0.4,
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
#     scaled_distance = function(physeq, method, rescale=rescale){
#         Dist = phyloseq::distance(physeq, method, type = "samples")
#         if (rescale) {
#             # rescale the distance matrix to be [0, 1]
#             Dist <- Dist / max(Dist, na.rm = TRUE)
#             Dist <- Dist - min(Dist, na.rm = TRUE)
#         }
#         return(Dist)
#     }
    #Distance = scaled_distance(physeq, distance)
    Distance = phyloseq::distance(physeq, distance)

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
    if (is.null(map)) map <- leaflet(data) %>% addTiles()

    map <- leaflet()

    # add lines to map
    for (g in unique(distdf$rowname)) {
        sdf <- distdf[distdf$rowname == g, ]
        map <- map %>%
            addPolylines(data = sdf,
                         lng = ~lng,
                         lat = ~lat,
                         weight = ~distance*5,
                         color = line_color,
                         opacity = line_alpha)
      }

    #add points to map
    map <- map %>% addCircleMarkers(data = data,
                                    radius = ~circlesize,
                                    color = makecolors(data, color),
                                    opacity = circle_alpha,
                                    fillOpacity = fillOpacity,
                                    popup = ~samplename)

    pal <- colorBin(palette = "YlGnBu",
                    domain = distdf$distance,
                    bins = 5)

    map <- map %>% addLegend("bottomright",
                             pal = pal,
                             values = ~distdf$distance,
                             labels = c("Test Labels"),
                             title = "Ecological Distance",
                             opacity = 1)
    return(map)
}

#' add a geojson worldmap layer
#'
#' @import dplyr
#' @import leaflet
#'
#' @param m  a Leaflet Map
#' @param ... paramters to pass to \code{[leaflet] addGeoJSON}
#'
#' @return a Leaflet Map
#'
#' @export
addWorld <- function(m, ...){
    if (!require("phylogeo")) library("phylogeo")
    worldjson <- paste0(path.package("phylogeo"),"/data/world-110m.json")
    world <-  readLines(worldjson) %>% paste(collapse = "\n")
    m %>% addGeoJSON(world, options = ...)
}
#' add a geojson NYC layer
#'
#' @import dplyr
#' @import leaflet
#'
#' @param m  a Leaflet Map
#' @param ... paramters to pass to \code{[leaflet] addGeoJSON} fro path sylting
#'
#' @return a Leaflet Map
#'
#' @export
addNYC <- function(m, ... ){
    if (!require("phylogeo")) library("phylogeo")
    nycjson <- paste0(path.package("phylogeo"),"/data/nycboroughs.json")
    nyc <-  readLines(nycjson) %>% paste(collapse = "\n")
    m <- leaflet()
    m %>% addGeoJSON(nyc, options = ...) %>% setView(lng = -73.97, lat = 40.72, zoom = 10)
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
