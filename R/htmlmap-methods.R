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
    data$circlesize <- sample_sums(physeq)
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

