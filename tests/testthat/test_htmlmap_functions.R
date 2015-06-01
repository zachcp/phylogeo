################################################################################
# Use testthat to test phyloseq data objects for the presence of Lat/Long
################################################################################
library("phylogeo")
library("testthat")
library("dplyr")
context("HTML Plotting Functions")


test_that("makecolors utlity function correctly handles column data", {
    data(epoxomicin_KS)
    #check a string color and nuemric factor
    blue = phylogeo:::makecolors(sample_data(epoxomicin_KS),"blue")
    state = phylogeo:::makecolors(sample_data(epoxomicin_KS),"state")
    
    expect_equal(blue, "blue")
    expect_equal(length(state), 59)
    expect_equal(length(unique(state)), 16)
    
})


