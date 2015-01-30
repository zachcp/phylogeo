################################################################################
# Use testthat to test phyloseq data objects for the presence of Lat/Long
################################################################################
library("phylogeo")
library("testthat")
context("Check Projections")

test_that("basic projections (that don't need extra values) work", {
  data(batfecal)
  expect_is(map_phyloseq(batfecal, proj="aitoff") , "ggplot")
  expect_is(map_phyloseq(batfecal, proj="azequalarea") , "ggplot")
  expect_is(map_phyloseq(batfecal, proj="cylindrical") , "ggplot")
  expect_is(map_phyloseq(batfecal, proj="gilbert") , "ggplot")
  expect_is(map_phyloseq(batfecal, proj="guyou") , "ggplot")
  expect_is(map_phyloseq(batfecal, proj="laue") , "ggplot")
  expect_is(map_phyloseq(batfecal, proj="mercator") , "ggplot")
  expect_is(map_phyloseq(batfecal, proj="mollweide") , "ggplot")
  expect_is(map_phyloseq(batfecal, proj="orthographic") , "ggplot")
  expect_is(map_phyloseq(batfecal, proj="polyconic") , "ggplot")
  expect_is(map_phyloseq(batfecal, proj="sinusoidal") , "ggplot")
  expect_is(map_phyloseq(batfecal, proj="tetra") , "ggplot")
})

test_that("Albers Projection works", {
  data(batfecal)
  p1 <- map_phyloseq(batfecal, proj="albers", paramaters = c(10,50))
  expect_is(p1, "ggplot")
})

#other projections

#expect_is(map_phyloseq(batfecal, proj="albers") , "ggplot")
#expect_is(map_phyloseq(batfecal, proj="bicentric") , "ggplot")
#expect_is(map_phyloseq(batfecal, proj="bonne") , "ggplot")
#expect_is(map_phyloseq(batfecal, proj="conic") , "ggplot")
#map_phyloseq(batfecal, proj="cylequalarea")
#expect_is(map_phyloseq(batfecal, proj="eisenlohr") , "ggplot")
#expect_is(map_phyloseq(batfecal, proj="elliptic") , "ggplot")
# expect_is(map_phyloseq(batfecal, proj="fisheye") , "ggplot")
# expect_is(map_phyloseq(batfecal, proj="gall") , "ggplot")
#expect_is(map_phyloseq(batfecal, proj="harrison") , "ggplot")
#expect_is(map_phyloseq(batfecal, proj="hex") , "ggplot")
#expect_is(map_phyloseq(batfecal, proj="homing") , "ggplot")
#expect_is(map_phyloseq(batfecal, proj="lambert") , "ggplot")
#expect_is(map_phyloseq(batfecal, proj="lagrange") , "ggplot")
#expect_is(map_phyloseq(batfecal, proj="lune") , "ggplot")
#expect_is(map_phyloseq(batfecal, proj="newyorker") , "ggplot")
#  expect_is(map_phyloseq(batfecal, proj="perspective") , "ggplot")
# expect_is(map_phyloseq(batfecal, proj="rectangular") , "ggplot")
#  expect_is(map_phyloseq(batfecal, proj="simpleconic") , "ggplot")
#  expect_is(map_phyloseq(batfecal, proj="trapezoidal") , "ggplot")
