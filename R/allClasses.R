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
                               sam_data  = NULL,
                               phy_tree  = NULL,
                               refseq    = NULL,
                               latitude  = NULL,
                               longitude = NULL)
)