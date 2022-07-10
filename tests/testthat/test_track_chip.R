testthat::context("track_chip")
library(ssvTracks)

bam_files = dir(system.file(package = "ssvTracks", "extdata"), full.names = TRUE, pattern = "bam$")
bw_files = dir(system.file(package = "seqsetvis", "extdata"), full.names = TRUE, pattern = "MCF.+bw$")
peak_files = dir(system.file(package = "seqsetvis", "extdata"), full.names = TRUE, pattern = "MCF.+Peak$")
peak_grs = seqsetvis::easyLoad_narrowPeak(peak_files)
query_gr = peak_grs[[1]][1]

gtf_file = dir(system.file(package = "ssvTracks", "extdata"), full.names = TRUE, pattern = "at_peaks.gtf$")

ex_gr = rtracklayer::import.gff(gtf_file, feature.type = "exon")


files_list = list(
  bw_files,
  bam_files,
  bam_files
)

fetch_fun_list = list(
  seqsetvis::ssvFetchBigwig,
  seqsetvis::ssvFetchBam,
  seqsetvis::ssvFetchBamPE
)

stopifnot(length(files_list) == length(fetch_fun_list))

for(i in seq_along(files_list)){
  fetch_files = files_list[[i]]
  fetch_fun = fetch_fun_list[[i]]

  #QcConfigFeatures
  test_that("track_chip no spline", {
    p = track_chip(fetch_files, query_gr, fetch_fun = fetch_fun, nspline = 1)
    testthat::expect_is(p, "ggplot")
  })

  test_that("track_chip default spline", {
    p = track_chip(fetch_files, query_gr, fetch_fun = fetch_fun)
    testthat::expect_is(p, "ggplot")
  })

  test_that("track_chip low win", {
    p = track_chip(fetch_files, resize(query_gr, 50*width(query_gr), fix = "center"), fetch_fun = fetch_fun, nwin = 10)
    testthat::expect_is(p, "ggplot")
  })

  test_that("track_chip max win (single bp)", {
    p = track_chip(fetch_files, resize(query_gr, 50*width(query_gr), fix = "center"), fetch_fun = fetch_fun, nwin = Inf, nspline = 1)
    testthat::expect_is(p, "ggplot")
  })

  bw_cfg_dt = data.table(file = fetch_files, sample = sub("_CTCF.+", "", basename(fetch_files)))

  test_that("track_chip by cfg", {
    p = track_chip(bw_cfg_dt, resize(query_gr, 50*width(query_gr), fix = "center"), fetch_fun = fetch_fun, nwin = 10)+
      labs(y = "FE", title = "CTCF")
    testthat::expect_is(p, "ggplot")
  })

  test_that("track_chip by cfg - show_lines", {
    p = track_chip(bw_cfg_dt, resize(query_gr, 50*width(query_gr), fix = "center"), fetch_fun = fetch_fun, nwin = 10, show_fill = FALSE, show_lines = TRUE)+
      labs(y = "FE", title = "CTCF")
    testthat::expect_is(p, "ggplot")
  })

  test_that("track_chip by cfg - fill_outline_color and legend.position", {
    p = track_chip(bw_cfg_dt,
                   resize(query_gr, 50*width(query_gr), fix = "center"),
                   fetch_fun = fetch_fun,
                   nwin = 10,
                   show_fill = TRUE,
                   show_lines = FALSE,
                   fill_outline_color = "black", legend.position = "none") +
      labs(y = "FE", title = "CTCF")
    testthat::expect_is(p, "ggplot")
  })

  test_that("track_chip by cfg - names_on_right", {
    p = track_chip(bw_cfg_dt,
                   resize(query_gr, 50*width(query_gr), fix = "center"),
                   fetch_fun = fetch_fun,
                   nwin = 10,
                   show_fill = TRUE,
                   show_lines = FALSE,
                   fill_outline_color = "black", legend.position = "none", names_on_right = FALSE) +
      labs(y = "FE", title = "CTCF")
    testthat::expect_is(p, "ggplot")
  })

}
