#' @title A phyloseq dataset from a microbiome study focusing on microbial community succession following forest fieres
#' @docType data
#' @name mountainsoil
#' @description This \code{\link[phyloseq]{phyloseq}} dataset is Study 1702 
#'     of the  \href{http://www.earthmicrobiome.org/}{Earth Microbiome Project}
#'     consisting of data derived from 22 soil samples in different elevations 
#'     in Changbai Mountain, northeastern China. The microbial communities in 
#'     these samples were investigated for differences due to elevation and 
#'     microbial succession after post-fire.
#' 
#' @format a \code{\link[phyloseq]{phyloseq}} object containing the following:
#' \describe{
#' \item{otu_table(mountainsoil)}{The observed taxonomic unit (OTU) table containing
#'       species abundance counts for each sample. For documentation on OTU 
#'       tables see the phyloseq's \code{\link[phyloseq]{otu_table}}}
#' \item{tax_table(mountainsoil)}{The taxonomy table providing taxonomy information
#'       for each OTU in the dataset. There are 5815 taxa and 7 taxonomic ranks
#'       labeled Rank1 - Rank7 and corresponding to Kingdom, Phylum... Species.
#'       For documentation on taxonomy tables see the phyloseq's 
#'       \code{\link[phyloseq]{tax_table}}}
#' \item{sample_data(mountainsoil)}{Data for each of the 22 samples in the dataset.
#'       The \code{\link[phyloseq]{sample_data}} dataframe has 73 columns
#'       containing information about sample collection, sample processing, 
#'       the physical-chemistry of the sample, and the location of the sample.
#'       For documentation about sample data see the phyloseq's 
#'       \code{\link[phyloseq]{sample_data}}}
#' }
#' 
#'                                                         
#' @examples
#' data(mountainsoil)
#' map_phyloseq(mountainsoil)
#' map_network(mountainsoil)
#' plot_distance(mountainsoil)
#' 
#' @usage
#' data(mountainsoil)
#'
#'@source 
#' \url{http://qiita.ucsd.edu/study/description/1702} 
#' 
#' @seealso
#'  \code{\link[phyloseq]{phyloseq}}
#'  \code{\link[phyloseq]{otu_table}}
#'  \code{\link[phyloseq]{tax_table}}
#'  \code{\link[phyloseq]{sample_data}}
NULL

#' @title  a phyloseq dataset from a bat-guano microbiome study
#' @docType data
#' @name batmicrobiome
#' @description This is data from Study 1734 of the 
#' \href{http://www.earthmicrobiome.org/}{Earth Microbiome Project}
#' which investigates the gut microbiota of phyllostomid bats that span a
#' breadth of diets. Bat guano samples were obtained from the US, 
#' Ecuador and Belize.
#' 
#' @format a \code{\link[phyloseq]{phyloseq}} object containing the following:
#' \describe{
#' \item{otu_table(batmicrobiome)}{The observed taxonomic unit (OTU) table containing
#'       species abundance counts for each sample. For documentation on OTU 
#'       tables see the phyloseq's \code{\link[phyloseq]{otu_table}}}
#' \item{tax_table(batmicrobiome)}{The taxonomy table providing taxonomy information
#'       for each OTU in the dataset. There are 6174 taxa and 7 taxonomic ranks
#'       labeled Rank1 - Rank7 and corresponding to Kingdom, Phylum... Species.
#'       For documentation on taxonomy tables see the phyloseq's 
#'       \code{\link[phyloseq]{tax_table}}}
#' \item{sample_data(batmicrobiome)}{Data for each of the 94 samples in the dataset.
#'       The \code{\link[phyloseq]{sample_data}} dataframe has 58 columns
#'       containing information about sample collection, sample processing, 
#'       the physical-chemistry of the sample, and the location of the sample.
#'       For documentation about sample data see the phyloseq's 
#'       \code{\link[phyloseq]{sample_data}}}
#'}
#'       
#' @source \url{http://www.earthmicrobiome.org/} 
#' @examples
#' data(batmicrobiome)
#' map_phyloseq(batmicrobiome)
#' map_network(batmicrobiome)
#' plot_distance(batmicrobiome)
#' 
#' @usage
#' data(batmicrobiome)
#' 
#' @seealso
#'  \code{\link[phyloseq]{phyloseq}}
#'  \code{\link[phyloseq]{otu_table}}
#'  \code{\link[phyloseq]{tax_table}}
#'  \code{\link[phyloseq]{sample_data}}
NULL

#' @title a phyloseq dataset of Epoxomicin-like ketosynthase (KS) domains.
#' @docType data
#' @name epoxomicin_KS
#' @references 
#' Owen JG, Charlop-Powers Z., Smith AG., Ternei MA., Calle PY., Reddy BV., Montiel D. & Brady SF (2015)
#' Multiplexed metagenome mining using short DNA sequence tags facilitates targeted discovery 
#' of epoxyketone proteasome inhibitors. 
#' Proceedings of the National Academy of Sciences of the United States of America 112(14):4221-6
#' \url{http://www.ncbi.nlm.nih.gov/pubmed/25831524}
#' 
#' 
#' @details abstract from research article (quoted):
#' 
#'  In molecular evolutionary analyses, short DNA sequences are used to infer
#'  phylogenetic relationships among species. Here we apply this principle to
#'  the study of bacterial biosynthesis, enabling the targeted isolation of
#'  previously unidentified natural products directly from complex metagenomes.
#'  Our approach uses short natural product sequence tags derived from
#'  conserved biosynthetic motifs to profile biosynthetic diversity in the
#'  environment and then guide the recovery of gene clusters from metagenomic
#'  libraries. The methodology is conceptually simple, requires only a small
#'  investment in sequencing, and is not computationally demanding. To
#'  demonstrate the power of this approach to natural product discovery we
#'  conducted a computational search for epoxyketone proteasome inhibitors
#'  within 185 globally distributed soil metagenomes. This led to the
#'  identification of 99 unique epoxyketone sequence tags, falling into
#'  6 phylogenetically distinct clades. Complete gene clusters associated with
#'  nine unique tags were recovered from four saturating soil metagenomic
#'  libraries. Using heterologous expression methodologies, seven potent
#'  epoxyketone proteasome inhibitors (clarepoxcins A-E and landepoxcins
#'  A and B) were produced from these pathways, including compounds with
#'  different warhead structures and a naturally occurring halohydrin
#'  prodrug. This study provides a template for the targeted expansion of
#'  bacterially derived natural products using the global metagenome.
#'  
#'  (end quote)
#' 
#' @description 
#'  a \code{\link[phyloseq]{phyloseq}}object for a soil microbiome
#'  study using degenerate primers targeting ketosynthase domains (KS), a 
#'  conserved domain in the biosynthesis of 
#'  \href{http://en.wikipedia.org/wiki/Polyketide}{polyketides}. This data is
#'  the subset of KS amplicons that map to the epoxyketone natural product 
#'  \href{https://en.wikipedia.org/wiki/Epoxomicin}{epoxomicin} and were used
#'  to guide the discovery of a nubmer of Epoxomicin-like compounds. Note that
#'  there are 102 sequences due to the addition of reference sequences (BCW, SAL
#'   and Epn) from three. Also, there only 162 samples in the data because these  
#'   are the sampels with one or more Epoxomicin-like sequence.
#' 
#' @format a \code{\link[phyloseq]{phyloseq}} object containing the following:
#' \describe{
#' \item{otu_table(epoxomicin_KS)}{The observed taxonomic unit (OTU) table containing
#'       species abundance counts for each sample. For documentation on OTU 
#'       tables see phyloseq's \code{\link[phyloseq]{otu_table}}} command.
#' \item{tax_table(epoxomicin_KS)}{The taxonomy table providing taxonomy information
#'       for each OTU in the dataset. There are 257 taxa and 7 taxonomic ranks.
#'       For documentation on taxonomy tables see the phyloseq's 
#'       \code{\link[phyloseq]{tax_table}}}
#' \item{sample_data(epoxomicin_KS)}{Data for each of the 162 samples in the dataset.
#'       The \code{\link[phyloseq]{sample_data}} dataframe has 10 columns
#'       containing basic information about the sample's location.
#'       For documentation about sample data see the phyloseq's 
#'       \code{\link[phyloseq]{sample_data}}}
#' \item{phys_tree(epoxomicin_KS)}{a phylogenetic tree of the 257 
#'       OTUs/sequences in the dataset. 
#'       For documentation about tree data see the phyloseq's 
#'       \code{\link[phyloseq]{phy_tree}}}
#'       
#'}
#' @examples
#' data(epoxomicin_KS)
#' map_phyloseq(epoxomicin_KS)
#' map_tree(epoxomicin_KS)
#' map_network(epoxomicin_KS)
#' map_clusters(epoxomicin_KS)
#' plot_distance(epoxomicin_KS)
#' 
#' @usage
#' data(epoxomicin_KS)
#' 
#' @seealso
#'  \code{\link[phyloseq]{phyloseq}}
#'  \code{\link[phyloseq]{otu_table}}
#'  \code{\link[phyloseq]{tax_table}}
#'  \code{\link[phyloseq]{sample_data}}
#'  \code{\link[phyloseq]{phy_tree}}
#'  \code{\link[phyloseq]{plot_tree}}
NULL