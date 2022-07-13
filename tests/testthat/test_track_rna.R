testthat::context("track_rna")
library(ssvTracks)
library(testthat)

bam_files = dir(system.file(package = "ssvTracks", "extdata"), full.names = TRUE, pattern = "RNA.+bam$")
bed_gr = rtracklayer::import.bed(system.file(package = "ssvTracks", "extdata/ESR1.bed"))

query_gr = bed_gr

gtf_file = dir(system.file(package = "ssvTracks", "extdata"), full.names = TRUE, pattern = "at_peaks.gtf$")
ex_gr = rtracklayer::import.gff(gtf_file, feature.type = "exon")


fetch_fun = seqsetvis::ssvFetchBam

options(mc.cores = 10)

#QcConfigFeatures
test_that("track_rna no spline", {
  p = track_rna(bam_files, query_gr, fetch_fun = fetch_fun, nspline = 1)
  p
  testthat::expect_is(p, "ggplot")
})

test_that("track_rna other strand, default spline", {
  p = track_rna(bam_files, query_gr, fetch_fun = fetch_fun, target_strand = "-")
  p
  testthat::expect_is(p, "ggplot")
})

test_that("track_rna low win", {
  p = track_rna(bam_files, query_gr, fetch_fun = fetch_fun, nwin = 10, nspline = 5)
  p
  testthat::expect_is(p, "ggplot")
})

test_that("track_rna max win (single bp)", {
  p = track_rna(bam_files, resize(query_gr, 500), fetch_fun = fetch_fun, nwin = Inf, nspline = 1, x_scale = "bp")
  p
  testthat::expect_is(p, "ggplot")
})

bw_cfg_dt = data.table(file = bam_files, sample = sub("_CTCF.+", "", basename(bam_files)))
bw_cfg_dt[, c("cell", "rep") := tstrsplit(sample, "[_\\.]", keep = 1:2)]
bw_cfg_dt[, sample := paste(cell, rep, sep = "\n")]

test_that("track_rna by cfg", {
  p = track_rna(bw_cfg_dt, query_gr, fetch_fun = fetch_fun, nwin = 10)+
    labs(y = "FE", title = "CTCF")
  p
  testthat::expect_is(p, "ggplot")
})

test_that("track_rna by cfg - show_lines", {
  p = track_rna(bw_cfg_dt, query_gr, fetch_fun = fetch_fun, nwin = 10, show_fill = FALSE, show_lines = TRUE, fill_VAR = "cell", color_VAR = "rep", color_mapping = c("R1" = "black"))+
    labs(y = "FE", title = "CTCF")
  p
  testthat::expect_is(p, "ggplot")
})

test_that("track_rna by cfg - fill_outline_color and legend.position", {
  p = track_rna(bw_cfg_dt,
                resize(query_gr, 50*width(query_gr), fix = "center"),
                fetch_fun = fetch_fun,
                nwin = 10,
                show_fill = TRUE,
                show_lines = FALSE,
                fill_outline_color = "black", legend.position = "none") +
    labs(y = "FE", title = "CTCF")
  p
  testthat::expect_is(p, "ggplot")
})

test_that("track_rna by cfg - names_on_right", {
  p = track_rna(bw_cfg_dt,
                resize(query_gr, 50*width(query_gr), fix = "center"),
                fetch_fun = fetch_fun,
                nwin = 10,
                show_fill = TRUE,
                show_lines = FALSE,
                fill_outline_color = "black", legend.position = "none", names_on_right = FALSE) +
    labs(y = "FE", title = "CTCF")
  p
  testthat::expect_is(p, "ggplot")
})

