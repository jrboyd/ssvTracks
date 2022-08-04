testthat::context("track_rna.SE")
library(ssvTracks)
library(testthat)

bam_files = dir(system.file(package = "ssvTracks", "extdata"), full.names = TRUE, pattern = "^M.+RNA.+bam$")
bam_files = bam_files[!grepl("_wt_", bam_files) & !grepl("_RUNX1-ko1_", bam_files)]
bed_gr = rtracklayer::import.bed(system.file(package = "ssvTracks", "extdata/ESR1.bed"))

query_gr = bed_gr

gtf_file = dir(system.file(package = "ssvTracks", "extdata"), full.names = TRUE, pattern = "at_peaks.gtf$")
ex_gr = rtracklayer::import.gff(gtf_file, feature.type = "exon")


options(mc.cores = 10)

#QcConfigFeatures
test_that("track_rna.SE no spline", {
  p = track_rna.SE(bam_files, query_gr, nspline = 1)
  p
  testthat::expect_is(p, "ggplot")
})

test_that("track_rna.SE other strand, default spline", {
  p = track_rna.SE(bam_files, query_gr, target_strand = "-")
  p
  testthat::expect_is(p, "ggplot")
})

test_that("track_rna.SE low win", {
  p = track_rna.SE(bam_files, query_gr, nwin = 10, nspline = 5)
  p
  testthat::expect_is(p, "ggplot")
})

test_that("track_rna.SE max win (single bp)", {
  p = track_rna.SE(bam_files, resize(query_gr, 500), nwin = Inf, nspline = 1, x_scale = "bp")
  p
  testthat::expect_is(p, "ggplot")
})

bw_cfg_dt = data.table(file = bam_files, sample = sub("_CTCF.+", "", basename(bam_files)))
bw_cfg_dt[, c("cell", "rep") := tstrsplit(sample, "[_\\.]", keep = 1:2)]
bw_cfg_dt[, sample := paste(cell, rep, sep = "\n")]

test_that("track_rna.SE by cfg", {
  p = track_rna.SE(bw_cfg_dt, query_gr, nwin = 10)+
    labs(y = "FE", title = "CTCF")
  p
  testthat::expect_is(p, "ggplot")
})

test_that("track_rna.SE by cfg - show_lines", {
  p = track_rna.SE(bw_cfg_dt, query_gr, nwin = 10, show_fill = FALSE, show_lines = TRUE, fill_VAR = "cell", color_VAR = "rep", color_mapping = c("R1" = "black"))+
    labs(y = "FE", title = "CTCF")
  p
  testthat::expect_is(p, "ggplot")
})

test_that("track_rna.SE by cfg - fill_outline_color and legend.position", {
  p = track_rna.SE(bw_cfg_dt,
                   query_gr,
                   nwin = 10,
                   show_fill = TRUE,
                   show_lines = FALSE,
                   fill_outline_color = "black", legend.position = "none") +
    labs(y = "FE", title = "CTCF")
  p
  testthat::expect_is(p, "ggplot")
})

test_that("track_rna.SE by cfg - names_on_right and facet_VAR", {
  bw_cfg_dt[, facet := "all"]
  p = track_rna.SE(bw_cfg_dt,
                   query_gr,
                   nwin = 5e3, nspline = 10,
                   show_fill = TRUE,
                   show_lines = FALSE,
                   facet_VAR = "facet",
                   fill_alpha = .5,
                   fill_outline_color = "black", legend.position = "none", names_on_right = FALSE) +
    labs(y = "FE", title = "CTCF")
  p
  testthat::expect_is(p, "ggplot")
})

