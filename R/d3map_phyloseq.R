#' create d3 Map
#'
#' @import phyloseq
#' @importFrom jsonlite toJSON
#' @export
#'
#' @examples
#' #' \dontrun{
#' d3map_phyloseq(batmicrobiome, tax_column="Rank1", mapdata = "World50m")
#' d3map_phyloseq(epoxomicin_KS, tax_column="Reference")
#' }
d3map_phyloseq <- function(physeq, tax_column, mapdata="USA", omitNA=TRUE){
    if (!tax_column %in% colnames(tax_table(physeq))) {
        stop("you must split on a column inthe taxtable")
    }
    possible_maps = c("USA","World50m", "World110m", "NYC", "SouthAfrica")
    if (!mapdata %in% possible_maps) {
        stop(paste0("mapdata must be one of ", paste0(possible_maps, collapse = " ")))
    }

    #create phygeo object
    phygeo <- phylogeo(physeq)

    # Tabulate the counts for each value in the tax_table column for each sample
    genfac = factor(tax_table(phygeo)[, tax_column])
    gentab = apply(otu_table(phygeo), MARGIN = 2, function(x) {
        tapply(x, INDEX = genfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
    })

    #create a dataframe where count data is a list of dataframe elements
    # needed for JSON serializatio
    # see https://cran.r-project.org/web/packages/jsonlite/vignettes/json-mapping.pdf
    gendf <- data.frame(gentab)
    samdata  <- data.frame(sample_data(phygeo))
    # add explicit Latitude and Longitude (needed/used by D3.js)
    colnames(samdata)[colnames(samdata)==phygeo@latitude] <- "Latitude"
    colnames(samdata)[colnames(samdata)==phygeo@longitude] <- "Longitude"
    samdata['samplename'] <- rownames(samdata)

    #add samplesumdata
    sampleSumData = data.frame(totalreads = sample_sums(phygeo))
    sampleSumData['samplename'] <- rownames(sampleSumData)
    samdata <- merge(samdata, sampleSumData, by = "samplename")

    if (omitNA) {samdata <- na.omit(samdata)}
    #samplename and countdata will be used in d3.js. they should not be changed.
    samdata['countdata'] <- lapply(samdata['samplename'],
                                   function(x){ data.frame(t(gendf[x]))})

    #mapdata
    nyc = system.file("static/maps/nyc.js", package="phylogeo");
    SA = system.file("static/maps/southafrica.js", package="phylogeo");
    usa_10m = system.file("static/maps/us-10m.js", package="phylogeo");
    world_50m = system.file("static/maps/world-50m.js", package="phylogeo");
    world_100m = system.file("static/maps/world-110m.js", package="phylogeo");

    #background data
    deps <- system.file("static/js_css_dependencies.html", package="phylogeo");
    basecss <- system.file("static/style.css", package="phylogeo");

    # specific data
    nyc_basicmap_js = system.file("static/nycmap.js", package="phylogeo");
    SA_basicmap_js = system.file("static/southafricamap.js", package="phylogeo");
    usa_basicmap_js = system.file("static/usamap.js", package="phylogeo");
    world_basicmap_js = system.file("static/worldmap.js", package="phylogeo");


    #choose the map
    #"USA","World50m", "World110m, SouthAfrica"
    if (mapdata == "USA") {
        mapfile <- usa_10m
        mapJS <- usa_basicmap_js
    } else if(mapdata == "World50m") {
        mapfile <- world_50m
        mapJS <- world_basicmap_js
    } else if (mapdata == "World100m") {
        mapfile <- world_100m
        mapJS <- world_basicmap_js
    } else if (mapdata == "NYC") {
        mapfile <- nyc
        mapJS <- nyc_basicmap_js
    } else if (mapdata == "SouthAfrica") {
        mapfile <- SA
        mapJS <- SA_basicmap_js
    }

    body = paste(includeHTML(deps),
                 includeScript(mapfile),
                 paste0("<script>var pointdata = ", toJSON(samdata) ,"</script>"),
                 includeScript(mapJS),
                 includeCSS(basecss),
                 paste0("<body onload=drawMap()>"),
                 paste0("<select></select>"),
                 paste0("</body>"),
                 collapse="\n")

    servedata(body=body)
}


#' create d3 Map
#'
#' @import phyloseq
#' @importFrom jsonlite toJSON
#' @importFrom reshape2 melt
#' @export
#'
#' @examples
#' #' \dontrun{
#' d3map_phyloseq(batmicrobiome, tax_column="Rank1", mapdata = "World50m")
#' d3map_phyloseq(epoxomicin_KS, tax_column="Reference")
#' }
d3map_network <- function(physeq,
                          distance="jaccard",
                          mapdata="USA",
                          maxdist=0.8,
                          mindist=0,
                          omitNA=TRUE){

    possible_maps = c("USA","World50m", "World110m", "NYC")
    if (!mapdata %in% possible_maps) {
        stop(paste0("mapdata must be one of ", possible_maps))
    }

    #get distance
    if (inherits(distance,"dist")) {
        three_col_dist <- reshape2::melt(as.matrix(distance))
        names(three_col_dist) <- c("source","target","value")
    } else {
        Distance = phyloseq::distance(physeq, distance)
        three_col_dist <- reshape2::melt(as.matrix(Distance))
        names(three_col_dist) <- c("source","target","value")
    }
    three_col_dist <- three_col_dist[(three_col_dist$value > mindist) & (three_col_dist$value <= maxdist),]

    #create phygeo object
    phygeo <- phylogeo(physeq)
    samdata  <- data.frame(sample_data(phygeo))
    # add explicit Latitude and Longitude (needed/used by D3.js)
    colnames(samdata)[colnames(samdata)==phygeo@latitude] <- "Latitude"
    colnames(samdata)[colnames(samdata)==phygeo@longitude] <- "Longitude"
    samdata['samplename'] <- rownames(samdata)

    if (omitNA) {samdata <- na.omit(samdata)}

    #mapdata
    nyc = system.file("static/maps/nyc.js", package="phylogeo");
    SA = system.file("static/maps/southafrica.js", package="phylogeo");
    usa_10m = system.file("static/maps/us-10m.js", package="phylogeo");
    world_50m = system.file("static/maps/world-50m.js", package="phylogeo");
    world_100m = system.file("static/maps/world-110m.js", package="phylogeo");

    #background data
    deps <- system.file("static/js_css_dependencies.html", package="phylogeo");
    basecss <- system.file("static/style.css", package="phylogeo");

    # specific data
    nyc_basicmap_js = system.file("static/nycmap_network.js", package="phylogeo");
    usa_basicmap_js = system.file("static/usamap_network.js", package="phylogeo");
    world_basicmap_js = system.file("static/worldmap_network.js", package="phylogeo");


    #choose the map
    #"USA","World50m", "World110m"
    if (mapdata == "USA") {
        mapfile <- usa_10m
        mapJS <- usa_basicmap_js
    } else if(mapdata == "World50m") {
        mapfile <- world_50m
        mapJS <- world_basicmap_js
    } else if (mapdata == "World100m") {
        mapfile <- world_100m
        mapJS <- world_basicmap_js
    } else if (mapdata == "NYC") {
        mapfile <- nyc
        mapJS <- nyc_basicmap_js
    }

    body = paste(includeHTML(deps),
                 includeScript(mapfile),
                 paste0("<script>var pointdata = ", toJSON(samdata) ,"</script>"),
                 paste0("<script>var linedata = ", toJSON(three_col_dist) ,"</script>"),
                 includeScript(mapJS),
                 includeCSS(basecss),
                 paste0("<body onload=drawMap()>"),
                 paste0("</body>"),
                 collapse="\n")

    servedata(body=body)
}


#' Wrapper to setup a D3 Map
#'
#' Starts the httpuv web server and hosts a simple form including a file
#' upload to demo the multipart parser.
#' Adapted from https://github.com/jeroenooms/webutils/blob/master/R/demo_httpuv.R
#
#' @keywords internal
#' @importFrom stats runif
#' @importFrom utils browseURL getFromNamespace head str tail
#' @importFrom htmltools includeScript includeHTML includeCSS
#' @param jsondata
#' @param mapfile
#' @param htmlfile
#' @param mapjsfile
#' @param basecss
#' @param basecss
servedata <- function(body, port=8000){
    rook_handler <- function(env){
        # See Rook spec
        content_type <- env[["CONTENT_TYPE"]]
        http_method <- env[["REQUEST_METHOD"]]
        request_body <- env[["rook.input"]]$read()
        path <- env[["PATH_INFO"]]

        # Show HTML page for GET requests.
        if(http_method == "GET"){
            message("Received HTTP GET request: ", path)
            list (status = 200,
                  body = body,
                  headers = c("content-type" = "text/html"))}}

    # Start httpuv
    if(missing(port))
        port <- round(runif(1, 2e4, 5e4));
    url <- paste0("http://localhost:", port, "/")
    message("Opening ", url)
    browseURL(url)
    httpuv::runServer("0.0.0.0", port, list(call = rook_handler))
}
