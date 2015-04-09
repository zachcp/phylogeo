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
#' @description a package for mapping microbiome data built on top of 
#' the \code{\link[phyloseq]{phyloseq}} phylsoeq pacakge.
NULL

#' @title \code{\link[phyloseq]{phyloseq}} Phyloseq Object for a microbiome study focusing on bat guano
#' @docType data
#' @name batfecal
#' 
#' @description This dataset is Study 1702 of the  \href{http://www.earthmicrobiome.org/}{Earth Microbiome Project} consisting
#' of data derived from 22soil samples in different elevations in Changbai Mountain, northeastern China. The 
#' microbial communities in these samples were investigated for differences due to elevation and microbial 
#' succession after post-fire.
#' 
#' @format a \code{\link[phyloseq]{phyloseq}} object containing an \code{\link[phyloseq]{otu_table}}, 
#'  a \code{\link[phyloseq]{tax_table}} and \code{\link[phyloseq]{sample_data}}. is a dataframe with the following variables:
#' \describe{
#' \item{X.SampleID}{Unique Sample identifier}
#' \item{BarcodeSequence}{Sample-Specific Barcode}
#' \item{LinkerPrimerSequence}{The Primer Used for Amplification}
#' \item{NUMBER_SAMPLES_PROMISED}{Number of Samples; constant}
#' \item{NITRO_ORG_CARB_UNIT}{units}
#' \item{ASSIGNED_FROM_GEO}{how is the sample assigned}
#' \item{LAB_PERSON}{Who prepped the sample?}
#' \item{EXPERIMENT_CENTER}{Where was the work done}
#' \item{SAMPLE_PROGRESS}{Denoting that the sample was sequenced}
#' \item{RUN_PREFIX}{Internal Processing Variable}
#' \item{INVESTIGATION_TYPE}{Internal Processing Variable}
#' \item{TAXON_ID}{Internal Taxon Value}
#' \item{DEPTH}{Internal Depth Value}
#' \item{TOT_ORG_CARB}{Total Organic Carbon}
#' \item{ILLUMINA_TECHNOLOGY}{Sequencing Technology Used}
#' \item{COMMON_NAME}{metagenome}
#' \item{INCLUDES_TIMESERIES}{Using a tiem series?}
#' \item{EXTRACTED_DNA_AVAIL_NOW}{DNA Availability}
#' \item{STUDY_ABSTRACT}{Study Abstract}
#' \item{TARGET_SUBFRAGMENT}{16S Fragment}
#' \item{SAMPLE_LOCATION}{Where is the sample now}
#' \item{PROJECT_NAME}{Chu_Changbai_mountain_soils}
#' \item{ELEVATION}{Elevation}
#' \item{RUN_DATE}{When was the data acquired}
#' \item{PH_METH}{How was pH measured}
#' \item{PCR_PRIMERS}{Primer sequences}
#' \item{COLLECTION_DATE}{Collection date}
#' \item{ALTITUDE}{Altitude}
#' \item{TITLE}{Project Title}
#' \item{ENV_BIOME}{Environemental Biome}
#' \item{WATER_CONTENT_SOIL}{Water Content}
#' \item{PLATFORM}{Sequencing Platform}
#' \item{CARB_NITRO_RATIO}{Carbon to Nitroen Ratio}
#' \item{ANNUAL_SEASON_TEMP}{Annual Season Temperature}
#' \item{COUNTRY}{Country}
#' \item{PH}{Soil pH}
#' \item{STUDY_TITLE}{Study Title}
#' \item{STUDY_ALIAS}{Study Alias}
#' \item{ANONYMIZED_NAME}{Anonymized Sample Name}
#' \item{MICROBIAL_N_C_UNIT}{Nitrogen and Carbon Units}
#' \item{SAMPLE_CENTER}{Where the sample was processed}
#' \item{SAMP_SIZE}{Size of the sample}
#' \item{PRINCIPAL_INVESTIGATOR}{Principl Investigator}
#' \item{STUDY_DESCRIPTION}{Study Description}
#' \item{PHYSICAL_SAMP_AVAIL_NOW}{Internal Sample Variable}
#' \item{LONGITUDE}{Longitude}
#' \item{MIENS_COMPLIANT}{MIENS Compliant}
#' \item{STUDY_ID}{Study ID}
#' \item{EXPERIMENT_DESIGN_DESCRIPTION}{Experimental Design}
#' \item{Description_duplicate}{Soil Description (duplicate)}
#' \item{SEQUENCING_METH}{Sequencing Method}
#' \item{HAS_EXTRACTED_DATA}{Has Extracted Data}
#' \item{ENV_MATTER}{Environmental Source Material}
#' \item{TARGET_GENE}{Target Gene}
#' \item{SUBMIT_TO_INSDC}{Submit to INSDC}
#' \item{ENV_FEATURE}{Environmental Features}
#' \item{KEY_SEQ}{Internal Variable}
#' \item{SPATIAL_SERIES}{Spatial Series}
#' \item{MICROBIAL_NITRO}{Microbial Nitrogen}
#' \item{MICROBIAL_CARB}{Microbial Carbon}
#' \item{LIBRARY_CONSTRUCTION_PROTOCOL}{Library Construction Protocol}
#' \item{REGION}{Region}
#' \item{RUN_CENTER}{where was the sequencing done}
#' \item{SAMPLE_TYPE}{sample type}
#' \item{EMP_STATUS}{EMP Status}
#' \item{DEFAULT_EMP_STATUS}{Default EMP Status}
#' \item{EMP_PERSON}{Janet Jansson}
#' \item{TOT_NITRO}{Total Nitrogen}
#' \item{LATITUDE}{Latitude}
#' \item{HAS_PHYSICAL_SPECIMEN}{is hte physical sampel available}
#' \item{NUMBER_SAMPLES_COLLECTED}{Number of Samples Collected}
#' \item{STUDY_CENTER}{Study Center}
#' \item{Description}{Sample Description}
#' }
#'
#' @seealso
#'  \code{\link[phyloseq]{phyloseq}}
#' @seealso
#'  \code{\link[phyloseq]{otu_table}}
#' @seealso
#'  \code{\link[phyloseq]{tax_table}}
#' @seealso
#'  \code{\link[phyloseq]{sample_data}}


NULL

#' @title \code{\link[phyloseq]{phyloseq}} Phyloseq Object for Bat Microbiome Data 
#' @docType data
#' @name batmicrobiome
#' @description this is the \href{http://www.earthmicrobiome.org/}{Earth Microbiome Project} Study 1734 about bat guano from the US, Ecuador and Costa Rica
#' @format a \code{\link{phyloseq}} object
#' @seealso
#'  \code{\link[phyloseq]{phyloseq}}
NULL

#' @title \code{\link[phyloseq]{phyloseq}} Phyloseq object indicating the presence of \code{\href{http://en.wikipedia.org/wiki/Epoxomicin}Epoxamicin}-like
#' biosynthetic clusters.
#' @docType data
#' @name epoxamicin_KS
#' @description \code{\link[phyloseq]{phyloseq}} Phyloseq Object for a soil microbiome study using degenerate primers
#' targeting ketosynthase domains (KS), a conserved domain in the biosynthesis of polyketides. This data is
#' the subset of KS amplicons that map to the epoxyketone natural product Epoxamicin.
#' @format a \code{\link{phyloseq}} object
#' @seealso
#'  \code{\link[phyloseq]{phyloseq}}
#'  \code{\link[phyloseq]{plot_tree}}
#' @seealso 
#'   \href{http://www.ncbi.nlm.nih.gov/pubmed/25831524}{Multiplexed metagenome mining using short DNA sequence tags facilitates targeted discovery of epoxyketone proteasome inhibitors}.
NULL