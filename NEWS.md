#0.99.6.3
add the d3map_phyloseq() funciton which will create and load an interactive d3 map displaying
abundance by a variable in the tax table

#0.99.6.2
refactor. create a phylogeo class that handles data conversion/verification. No change to API.

#0.99.6.1
leaflet is now on CRAN. remove the hacks around importing it.
new documentation on gh-pages based on RStudio's leaflet page.
Change from gridExtra to cowplot for compound plots.

#0.99.6
Documentation fixes for BioC. Add mapdata package so that the worldHires map can be used. Leaflet hits CRAn so remove the import shims.

#0.99.5
More documentation changes and small fixes after another round of bioconductor code review.
Fixes to the projection output, function and dataset documentation, naespace import fixes.

#0.99.4
Moved all of the htmlwidgets/leaflet stuff to another branch in order to build for bioconductor.
when Rstudio/leaflet is available on CRAN I can merge i back to the master branch

#0.99.3
fixes for Bioconductor. Improved documentation, some code refactoring

#0.99.2
no new features. added email to biocdev and bumping so automated builds will pass

#0.99.1
cleanup related to projections, added test for projections, 
a vignette for projections, and an Rmd file documenting proecjtions intended for the web.

#0.99.0
bump version for submission to bioconductor. cleaned up docs, Rmd files, etc.

# 0.0.16
add projections to the mapping funtions by explosing map_proj() code.
Idea an initial coding courtest of @ryneches.

# 0.0.15
Basic Release with 4 major functions
- map_phyloseq
- map_network
- map_tree
- map_clusters
- plot_distance
