###############################################
#'  Mapping microbiome data.
#'
#' @name phylogeo
#' @author Zach Charlop-Powers \email{zcharlop@@rockefeller.edu}
#' @docType package
#' @import phyloseq
#' @import ggplot2
#' @import gridExtra
#' @keywords package
#' @description a package for mapping microbiome data
NULL

#' \code{\link[phyloseq]{phyloseq}} Phyloseq Object for a microbiome study focusing on bat guano
#'
#' A phyloseq object with data collected from the bat guano experiments
#' 
#' @docType data
#' @name batfecal
#' @description this is the Earth Microbiome Data study 1702 concernign bat guano from Chinese caves. 
#' @format a \code{\link{phyloseq}} object
NULL

#' \code{\link[phyloseq]{phyloseq}} Phyloseq Object for Bat Microbiome Data 
#'
#' A phyloseq object with data collected from the bat guano experiments
#' 
#' @docType data
#' @name batmicrobiome
#' @description this is the Earth Microbiome Data study 1734 about bat guano from the US, Ecuador and Costa Rica
#' @format a \code{\link{phyloseq}} object
NULL

#' \code{\link[phyloseq]{phyloseq}} Phyloseq Object for a soil microbiome study using degenerate primers
#' targeting ketosynthase domains (KS), a conserved domain in the biosynthesis of polyketides. This data is
#' the subset of KS amplicons that map to the epoxyketone natural product Epoxamicin.
#'
#' A phyloseq object with OTU data generated from PCR amplifying the KS domain from environmental DNA and keeping the
#' hits for Epoxamicin
#' 
#' @docType data
#' @name epoxamicin_KS
#' @format a \code{\link{phyloseq}} object
NULL