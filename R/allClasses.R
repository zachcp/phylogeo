################################################################################
#' The phylogeo class
#'
#' Inherits from the phyloseq class and adds the Latitude and longitude
#'
#' @import phyloseq
#' @name phylogeo-class
#' @rdname phylogeo-class
#' @exportClass phylogeo
setClass("phylogeo",
         representation(
             latitude = "character",
             longitude = "character"),
         contains = "phyloseq",
         prototype = prototype(otu_table = NULL,
                             tax_table = NULL,
                             sam_data = NULL,
                             phy_tree = NULL,
                             refseq = NULL,
                             latitude = NULL,
                             longitude = NULL)
)



#' initialize a phylogeo object from a phyloseq object
#'
#'
#' \code{phylogeo()} is a constructor method used internally to create a phylogeo
#' object and to ensure data integrity.
phylogeo <- function(physeq) {
    #get lat/lon and sampledata_table
    physeqdata = check_phyloseq(physeq)
    otus <- subset_otu_table(physeqdata, physeq)

    tree <- subset_tree(otus, physeq)
    seqs <- subset_seqs()

    new('phylogeo',
        latitude = physeqdata$lats,
        longitude = physeqdata$lngs,
        otu_table = otu_table(physeq, taxa_are_rows = TRUE),
        tax_table = tax_table(physeq),
        sam_data = physeqdata$sampledata,
        phy_tree = tree,
        refseq = seqs)
}

############################
#
# Data Validation
#
############################

#' check phyloseq
#'
#' function to check phyloseq data for lat non columns and NA values
#' return a list contiang the lat and lon columnnames and the sampledata
#'
#' @keywords internal
#' @import dplyr
check_phyloseq <- function(physeq){

    #check phyloseq objects for Lat/Lon
    if (!"sam_data" %in% phyloseq::getslots.phyloseq(physeq)) {
        stop("Mapping requires that phyloseq objects have Sample_Data with Latitude
             and Longitude")
    }

    #get lat/long
    #' From Rstudio/leaflet/R/normalize.R
    #' https://github.com/rstudio/leaflet/blob/4ef0023c9fefa00a64e382ccd77d34d1413c47dc/R/normalize.R
    sampledata <- data.frame(sample_data(physeq))
    sdfnames <- names(sampledata)
    lats = sdfnames[grep("^(lat|latitude)$", sdfnames, ignore.case = TRUE)]
    lngs = sdfnames[grep("^(lon|lng|long|longitude)$", sdfnames, ignore.case = TRUE)]

    if (!(length(lats) == 1 && length(lngs) == 1)) stop("Couldn't infer longitude/latitude columns")

    #get sample data info
    sampledata <- sampledata %>%
        check_NA(lats) %>%
        check_NA(lngs) %>%
        coerce_numeric(lats) %>%
        coerce_numeric(lngs)

    return(list(lat = lats,
                lng = lngs,
                sampledata = sampledata))
    }




subset_otu_table <- function(){
}

subset_tree <- function(){
}

subset_seqs <- function(){
}

