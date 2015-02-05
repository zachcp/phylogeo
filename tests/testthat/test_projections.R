################################################################################
# Use testthat to test phyloseq use of mapproj projections
################################################################################
library("phylogeo")
library("testthat")
context("Check Projections")

test_that("basic projections (that don't need extra values) work", {
  data(batfecal)
  expect_is(map_phyloseq(batfecal, projection="aitoff") , "ggplot")
  expect_is(map_phyloseq(batfecal, projection="azequalarea") , "ggplot")
  expect_is(map_phyloseq(batfecal, projection="cylindrical") , "ggplot")
  expect_is(map_phyloseq(batfecal, projection="gilbert") , "ggplot")
  expect_is(map_phyloseq(batfecal, projection="guyou") , "ggplot")
  expect_is(map_phyloseq(batfecal, projection="laue") , "ggplot")
  expect_is(map_phyloseq(batfecal, projection="mercator") , "ggplot")
  expect_is(map_phyloseq(batfecal, projection="mollweide") , "ggplot")
  expect_is(map_phyloseq(batfecal, projection="orthographic") , "ggplot")
  expect_is(map_phyloseq(batfecal, projection="polyconic") , "ggplot")
  expect_is(map_phyloseq(batfecal, projection="sinusoidal") , "ggplot")
  expect_is(map_phyloseq(batfecal, projection="tetra") , "ggplot")
})

test_that("Projections with Lat0 work", {
  data(batfecal)
  expect_is(map_phyloseq(batfecal, projection="cylequalarea", lat0=15), "ggplot")
  expect_is(map_phyloseq(batfecal, projection="rectangular", lat0=15), "ggplot")
  expect_is(map_phyloseq(batfecal, projection="conic", lat0=15), "ggplot")
  expect_is(map_phyloseq(batfecal, projection="mecca", lat0=15), "ggplot")
  expect_is(map_phyloseq(batfecal, projection="homing", lat0=15) , "ggplot")
})

test_that("Fisheye projection works", {
  data(batfecal)
  expect_is(map_phyloseq(batfecal, projection="fisheye", n=0.5) , "ggplot")
})

test_that("NewYorker projection works", {
  data(batfecal)
  expect_is(map_phyloseq(batfecal, projection="newyorker", r=0.3) , "ggplot")
})

test_that("Projections with two latitude parameters work: ", {
  data(batfecal)
  expect_is(map_phyloseq(batfecal, projection="simpleconic", lat0=10, lat1=10) , "ggplot")
  expect_is(map_phyloseq(batfecal, projection="lambert", lat0=10, lat1=10) , "ggplot")
  expect_is(map_phyloseq(batfecal, projection="albers", lat0=10, lat1=10) , "ggplot")
  expect_is(map_phyloseq(batfecal, projection="trapezoidal", lat0=10, lat1=10) , "ggplot")
})

test_that("Projections with longitude component work: ", {
  data(batfecal)
  expect_is(map_phyloseq(batfecal, projection="bicentric", lon0=10) , "ggplot")
  expect_is(map_phyloseq(batfecal, projection="elliptic", lon0=10) , "ggplot")
})

# test_that("Test the Harrison Projection: ", {
#   data(batfecal)
#   expect_is(map_phyloseq(batfecal, projection="harrison", dist=10, angle=10) , "ggplot")
# })
# 
# test_that("Test the Lune Projection: ", {
#   data(batfecal)
#   expect_is(map_phyloseq(batfecal, projection="harrison", lat=10, angle=10) , "ggplot")
# })