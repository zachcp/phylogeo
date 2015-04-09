################################################################################
# Use testthat to test phyloseq use of mapproj projections
################################################################################
library("phylogeo")
library("testthat")
context("Check Projections")


#######
#still need the following
#"eisenlohr", "harrison","lune","perspective","stereographic"

test_that("basic projections (that don't need extra values) work", {
  data(mountainsoil)
  expect_is(map_phyloseq(mountainsoil, projection="aitoff") , "ggplot")
  #expect_is(map_phyloseq(mountainsoil, projection="bonne") , "ggplot")
  expect_is(map_phyloseq(mountainsoil, projection="azequalarea") , "ggplot")
  expect_is(map_phyloseq(mountainsoil, projection="azequidistant") , "ggplot")
  #expect_is(map_phyloseq(mountainsoil, projection="cylindrical") , "ggplot")
  expect_is(map_phyloseq(mountainsoil, projection="gilbert") , "ggplot")
  expect_is(map_phyloseq(mountainsoil, projection="globular") , "ggplot")
  expect_is(map_phyloseq(mountainsoil, projection="gnomonic") , "ggplot")  
  expect_is(map_phyloseq(mountainsoil, projection="guyou") , "ggplot")
  expect_is(map_phyloseq(mountainsoil, projection="hex") , "ggplot")
  expect_is(map_phyloseq(mountainsoil, projection="lagrange") , "ggplot")
  expect_is(map_phyloseq(mountainsoil, projection="laue") , "ggplot")
  expect_is(map_phyloseq(mountainsoil, projection="mercator") , "ggplot")
  expect_is(map_phyloseq(mountainsoil, projection="mollweide") , "ggplot")
  expect_is(map_phyloseq(mountainsoil, projection="orthographic") , "ggplot")
  expect_is(map_phyloseq(mountainsoil, projection="polyconic") , "ggplot")
  expect_is(map_phyloseq(mountainsoil, projection="sinusoidal") , "ggplot")
  expect_is(map_phyloseq(mountainsoil, projection="square") , "ggplot")
  expect_is(map_phyloseq(mountainsoil, projection="tetra") , "ggplot")
})

test_that("Projections with Lat0 work", {
  data(mountainsoil)
  expect_is(map_phyloseq(mountainsoil, projection="conic", lat0=15), "ggplot")
  expect_is(map_phyloseq(mountainsoil, projection="cylequalarea", lat0=15), "ggplot")
  #expect_is(map_phyloseq(mountainsoil, projection="gall", lat0=15) , "ggplot")
  expect_is(map_phyloseq(mountainsoil, projection="homing", lat0=15) , "ggplot")
  expect_is(map_phyloseq(mountainsoil, projection="mecca", lat0=15), "ggplot")
  expect_is(map_phyloseq(mountainsoil, projection="rectangular", lat0=15), "ggplot")
})

test_that("Fisheye projection works", {
  data(mountainsoil)
  expect_is(map_phyloseq(mountainsoil, projection="fisheye", n=0.5) , "ggplot")
})

test_that("NewYorker projection works", {
  data(mountainsoil)
  expect_is(map_phyloseq(mountainsoil, projection="newyorker", r=0.3) , "ggplot")
})

test_that("Projections with two latitude parameters work: ", {
  data(mountainsoil)
  expect_is(map_phyloseq(mountainsoil, projection="albers", lat0=10, lat1=10) , "ggplot")
  expect_is(map_phyloseq(mountainsoil, projection="lambert", lat0=10, lat1=10) , "ggplot")
  expect_is(map_phyloseq(mountainsoil, projection="simpleconic", lat0=10, lat1=10) , "ggplot")
  expect_is(map_phyloseq(mountainsoil, projection="trapezoidal", lat0=10, lat1=10) , "ggplot")
})

test_that("Projections with longitude component work: ", {
  data(mountainsoil)
  expect_is(map_phyloseq(mountainsoil, projection="bicentric", lon0=10) , "ggplot")
  expect_is(map_phyloseq(mountainsoil, projection="elliptic", lon0=10) , "ggplot")
})