################################################################################
# Use testthat to test phyloseq data objects for the presence of Lat/Long
################################################################################
library("phylogeo")
library("phyloseq")
library("testthat")
library("dplyr")
context("Utility Functions")

# Test the following funcitons
# .check_physeq <- function(physeq){
# .create_basemap <-function(region, df, latcol, loncol){
# .check_names <- function(member, df, allownumeric=FALSE){
# .jitter_df <- function(df, xcol, ycol, jitter.x, jitter.y){
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