################################################################################
# Use testthat to test phyloseq data objects plotting functions
################################################################################
library("phylogeo") 
library("testthat") 
library("igraph")

test_that("map_phyloseq: results in a ggplot", {
  data(batfecal)
	p1 <- map_phyloseq(batfecal)
	# ggplot-class tests
	expect_is(p1, "ggplot")
	})

test_that("map_clusters: returns a ggplot", {
  data(batfecal) #no tree shoud fail 
  data(epoxamicin_KS) #has tree should work
  p1 <- map_clusters(epoxamicin_KS)
  # ggplot-class tests
  expect_error(map_clusters(batfecal))
  expect_is(p1, "ggplot")
})

test_that("map_tree: tree_objects are present when mapping a phylogenetic tree", {
  data(batfecal)
  data(epoxamicin_KS)
	p1 <- map_tree(epoxamicin_KS)
	# ggplot-class tests
  expect_error(map_tree(batfecal))
	expect_is(p1, "ggplot")
	})
    
test_that("map_network: can use a precalculated network", {
  data(batmicrobiome)
  p1 <- map_network(batmicrobiome)
  ig <- make_network(batmicrobiome)
  p2 <- map_network(batmicrobiome, igraph=ig)
	# ggplot-class tests
  expect_is(p1, "ggplot")
  expect_is(p2, "ggplot")
  expect_is(ig, "igraph")
	})

