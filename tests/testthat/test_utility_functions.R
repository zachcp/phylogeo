################################################################################
# Use testthat to test phyloseq data objects for the presence of Lat/Long
################################################################################
library("phylogeo")
library("testthat")
library("dplyr")
context("Utility Functions")

# Test the following funcitons
# .check_physeq <- function(physeq){
# .create_basemap <-function(region, df, latcol, loncol){
# .check_NA <- function(df, col){
# .coerce_numeric <- function(df, col)
# .dist_to_edge_table = function(Dist, dname = "dist")
# .degree_to_radian <- function(degree) 
# .radian_to_degree <- function(radian) 


test_that("phyloseq objects have latitude and longitude columns", {
  data(batfecal)
  data(batmicrobiome)
  data(epoxamicin_KS)
  latlist  <- list( c("LATITUDE"),c("LONGITUDE"))
  latlist2 <- list( c("Latitude"),c("Longitude"))
  expect_equal(.check_physeq(batfecal), latlist )
  expect_equal(.check_physeq(batmicrobiome), latlist )
  expect_equal(.check_physeq(epoxamicin_KS), latlist2 )
  
  # remove lat.long column and be sure it trips an error
  df <- data.frame(sample_data(batfecal))
  dflat <- df %>% dplyr::select(-LATITUDE)
  dflon <- df %>% dplyr::select(-LONGITUDE)
  bat_lat <- batfecal
  sample_data(bat_lat) <- dflat
  expect_error(.check_physeq(bat_lat))
  bat_lon <- batfecal
  sample_data(bat_lon) <- dflon
  expect_error(.check_physeq(bat_lon))
})

test_that("jittering works and uses the seed correctly", {
  df1 <- data.frame(x= c(1,2,3,4,5), y=c(3,4,5,6,7))
  df2 <- data.frame(x= c(1,2,3,4,5), y=c(3,4,5,6,7))
  df1_j <- .jitter_df(df1, "x","y",3,4,123)
  df2_j <- .jitter_df(df2, "x","y",3,4,123)
  
  expect_is(df1_j , "data.frame")
  expect_is(df2_j , "data.frame")
  expect_equal(df1_j, df2_j)
})

test_that("checknames trips an error when its supposed to",{
  df <- data.frame(col1=c(1,2,3),
                   col2=c(1,2,3))
  
  expect_output(.check_names("col1", df),"")
  expect_output(.check_names("col2", df),"")
  expect_warning(.check_names("col3", df))
})