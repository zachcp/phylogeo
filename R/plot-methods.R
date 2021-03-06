#
# methods for plotting geographic-distance plots
#
###############################################################################
#' Calculate and Plot Sample Distances by Geography/Ecological Distance
#'
#' plot_distance is a plotting function for visualizing the geographic
#' and ecological distances between all pairs of samples in a microbiome study.
#' This function acts on \code{\link[phyloseq]{phyloseq}} datasets and
#' requires that the \code{\link[phyloseq]{sample_data}} table
#' contains Latitude and Longitude columns. This will calcualte the pairwise
#' distances between each set of samples and plot them as a scatter plot.
#'
#' Any of the ecological distances supported by phyloseq
#' \code{\link[phyloseq]{distance}} will work here and can be specified with
#' "distancemethod="
#'
#' @usage plot_distance(physeq, distancemethod="jaccard")
#' @return a ggplot object
#'
#' @param physeq (Required).
#'  The name of the \code{\link[phyloseq]{phyloseq}} dataset. This must have
#'  sample data with Latitude and Longitude Columns.
#'
#' @param distancemethod (Optional). Default \code{"jaccard"}.
#'  The name of an ecological distance method.
#'  See ?distance for more information
#'
#' @seealso \code{\link[phyloseq]{distance}}
#' @import phyloseq
#' @importFrom reshape2 melt
#' @importFrom sp spDists
#' @export
#' @examples
#' data(mountainsoil)
#' plot_distance(mountainsoil)
plot_distance <- function(physeq, distancemethod="jaccard"){
    #convert to phylogeo
    phygeo <- phylogeo(physeq)

    #get bigcircle distances using spDists
    #spDists expects the first column to be longitude
    df <- data.frame(sample_data(phygeo))
    df2 <- df[ c(phygeo@longitude, phygeo@latitude) ]
    names(df2) <- c("lon", "lat")
    df2 <- as.matrix(df2)
    geodistances <- spDists(df2, longlat = TRUE)
    colnames(geodistances)  <- row.names(df2)
    row.names(geodistances) <- row.names(df2)
    geodistances <- dist_to_edge_table(geodistances)
    names(geodistances) <- c("Var1","Var2","geodist")

    #get ecologicaldistances
    ecodistance <- distance(physeq, method = distancemethod)
    ecodistance <- dist_to_edge_table(ecodistance)
    names(ecodistance) <- c("Var1","Var2","ecodist")

    #make mergeable names for the two distance functions and merge
    concatvals <- function(x,y){ return(paste(x,"_",y,sep = ""))}
    geodistances['pairs'] <- mapply(concatvals, geodistances$Var1, geodistances$Var2)
    ecodistance['pairs'] <- mapply(concatvals, ecodistance$Var1, ecodistance$Var2)
    df <- merge(geodistances, ecodistance, by = "pairs")

    #make the plot
    p <- ggplot(df, aes(y = ecodist, x = geodist)) +
         geom_point() +
         xlab("Km") +
         ylab(distancemethod) +
         ggtitle( paste("Pairwise Sample Distance:",
                        "Km vs.", distancemethod,"Distance.", sep = " "))

    return(p)
}
