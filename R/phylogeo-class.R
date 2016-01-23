################################################################################
#' Build the phylogeo class object
#'
#' This is the suggested method for constructing phylogeo objects. 
#' 
#' @usage phylogeo(object, errorIfNULL=TRUE)
#'
#' @param object (Required). A \code{\link{data.frame-class}}, 
#'  or a \code{\link{phyloseq-class}} object.
#'
#' @return A \code{\link{phylogeo-class}} object
#' representing a \code{\link{phyloseq-class}} object with some
#' additonal Latitude and Longitude Information.
#'
#' @import phyloseq
#' @import dplyr
#' 
#' @rdname phylogeo-methods
#' @docType methods
#' @export
#'
setGeneric("phylogeo", function(object) standardGeneric("phylogeo"))
#' @rdname phylogeo-methods
#' @aliases phylogeo,data.frame-method
setMethod("phylogeo", "data.frame", function(object){
  #creates a phylogeo object from a data.frame
  
  #create dummyOTU table
  otutab <- matrix(data=1, nrow=dim(object)[[1]],ncol=2)
  colnames(otutab) <- c('dummy_col1','dummy_col2')
  rownames(otab) <- rownames(object)
  
  #load into phyloseq/phylogeo
  physeq <- phyloseq(sample_data(object), otu_table(otutab)) 
  phygeo <- phylogeo(physeq)
  warning("phylogeo objects created from data.frames are convenient
           for mapping smapel locations but should not be used for 
           network analyses.")
  return(phygeo)
})
#' @rdname phylogeo-methods
#' @aliases phylogeo,phyloseq-method
setMethod("phylogeo", "phyloseq", function(object){
 # check for sample_data slot
  if (!"sam_data" %in% phyloseq::getslots.phyloseq(physeq)) {
    stop("Mapping requires that phyloseq objects have sample_data 
          with Latitude and Longitude")
  }
  
  # get lat/long columns. (based on https://github.com/rstudio/leaflet/blob/4ef0023c9fefa00a64e382ccd77d34d1413c47dc/R/normalize.R)
  sampledata <- data.frame(sample_data(physeq))
  sdfnames   <- names(sampledata)
  lats       <- sdfnames[grep("^(lat|latitude)$", sdfnames, ignore.case = TRUE)]
  lngs       <- sdfnames[grep("^(lon|lng|long|longitude)$", sdfnames, ignore.case = TRUE)]
  if (!(length(lats) == 1 && length(lngs) == 1)) stop("Couldn't infer longitude/latitude columns")
  
  # Update sample_data and add it back to the phyloseq object
  sampledata <- sampledata %>%
    coerce_latlon_columns(lats) %>%
    coerce_latlon_columns(lngs)
  sample_data(physeq) <- sampledata
  
  # get samples without NA values in Lat/Lng columns
  samples_to_keep <- row.names(sampledata[ !is.na(sampledata[[lats]]) & !is.na(sampledata[[lngs]]), ])
  
  # if samples have been dropped prune OTUs belonging to them
  if (length(row.names) > 0 ) {
    #update sample data
    prune_samples(samples_to_keep, physeq)
    
    #update otu/sequence data
    if (.hasSlot(physeq, "otu_table")) {
      prune_taxa(taxa_sums(physeq) > 0, physeq)
    }
  }
  
  # return the class
  new('phylogeo',
      latitude  = lats,
      longitude = lngs,
      otu_table = physeq@otu_table,
      tax_table = physeq@tax_table,
      sam_data  = physeq@sam_data,
      phy_tree  = physeq@phy_tree,
      refseq    = physeq@refseq)
})


#' Address Lat//Long Columns:
#'  "None" -> NA
#'  factor -> vector
#'  vector -> numeric
#'
#' @keywords internal
coerce_latlon_columns <- function(df, col){
  colvals <- df[[col]]
  if (is.factor(colvals)) colvals <- as.vector(colvals) #convert to vector
  colvals[ colvals == "None"] <- NA  #some data has "None" so be sure to replace with NA
  df[[col]] <- as.numeric(as.character(colvals))
  df
}

