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
  data(mountainsoil)
  data(batmicrobiome)
  data(epoxomicin_KS)
  latlist  <- list( c("LATITUDE"),c("LONGITUDE"))
  latlist2 <- list( c("Latitude"),c("Longitude"))
  expect_equal(.check_physeq(mountainsoil), latlist)
  expect_equal(.check_physeq(batmicrobiome), latlist)
  expect_equal(.check_physeq(epoxomicin_KS), latlist2)
  
  # remove lat.long column and be sure it trips an error
  df <- data.frame(sample_data(mountainsoil))
  dflat <- df %>% dplyr::select(-LATITUDE)
  dflon <- df %>% dplyr::select(-LONGITUDE)
  bat_lat <- mountainsoil
  sample_data(bat_lat) <- dflat
  expect_error(.check_physeq(bat_lat))
  bat_lon <- mountainsoil
  sample_data(bat_lon) <- dflon
  expect_error(.check_physeq(bat_lon))
})

test_that("jittering works and uses the seed correctly", {
  df <- data.frame(x= c(1,2,3,4,5), y=c(3,4,5,6,7))
  df_j1 <- .jitter_df(df, "x","y",3,4,123)
  df_j2 <- .jitter_df(df, "x","y",3,4,123)
  expect_is(df_j1 , "data.frame")
  expect_is(df_j2 , "data.frame")
  #jittered dfs with same seed should be equal
  expect_equal(df_j1, df_j2)
})

test_that("checknames trips an error when its supposed to",{
  df <- data.frame(col1=c(1,2,3),
                   col2=c(1,2,3))
  
  expect_output(.check_names("col1", df),"")
  expect_output(.check_names("col2", df),"")
  expect_error(.check_names("col3", df))
})

test_that(".check_NA can handle columns that are factors", {
  #col2 has "NONE"
  df1 <- data.frame(col1=c(1,2,3,4,1,2,3),
                   col2=c(1,2,3,"None",1,2,3))
  #col2 is a factor with "None"
  df2 <- data.frame(col1=c(1,2,3,4,1,2,3),
                   col2=as.factor(c(1,2,3,"None",1,2,3)))
  expect_is(df2$col2, "factor")
  #check for warning
  #http://stackoverflow.com/questions/20885845/testthat-handling-both-warning-and-value
  expect_warning(val <-  .check_NA(df1, "col2"))
  expect_true(is.data.frame(val))
  expect_warning(val <-  .check_NA(df2, "col2"))
  expect_true(is.data.frame(val))
  
  #check that the none is removed  values
  testdf <- .check_NA(df1, "col2")
  expect_equal(testdf$col2, c("1","2","3","1","2","3"))
  testdf <- .check_NA(df2, "col2")
  expect_equal(testdf$col2, c("1","2","3","1","2","3"))
})