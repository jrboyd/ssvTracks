testthat::context("track_rna")
library(ssvTracks)
library(testthat)

bam_files = dir(system.file(package = "ssvTracks", "extdata"), full.names = TRUE, pattern = "RUNX1_RNA.+bam$")

query_gr = GRanges("chr21:34781678-35060256")
strand(query_gr) = "-"

if(FALSE){
  prof_dt = ssvFetchBamPE.RNA(bam_files, query_gr, return_data.table = TRUE)
  prof_dt[, x := (start + end)/2]
  splice_dt = ssvFetchBamPE.RNA_splice(bam_files, query_gr, return_data.table = TRUE)

  ggplot(prof_dt, aes(x = x, y = y, color = strand)) +
    geom_path() +
    facet_grid(sample~strand) +
    ggbio::geom_arch(data = splice_dt, aes(x = start, xend = end, height = N), color = "black") +
    coord_cartesian(xlim = c(start(query_gr), end(query_gr)))

  prof_dt = ssvFetchBamPE.RNA(bam_files, query_gr, return_data.table = TRUE, target_strand = "-")
  prof_dt[, x := (start + end)/2]
  splice_dt = ssvFetchBamPE.RNA_splice(bam_files, query_gr, return_data.table = TRUE, target_strand = "-")

  ggplot(prof_dt, aes(x = x, y = y, color = strand)) +
    geom_path() +
    facet_grid(sample~strand) +
    ggbio::geom_arch(data = splice_dt, aes(x = start, xend = end, height = N), color = "black") +
    coord_cartesian(xlim = c(start(query_gr), end(query_gr)))

  prof_dt = ssvFetchBamPE.RNA(bam_files, query_gr, return_data.table = TRUE, target_strand = "+")
  prof_dt[, x := (start + end)/2]
  splice_dt = ssvFetchBamPE.RNA_splice(bam_files, query_gr, return_data.table = TRUE, target_strand = "+")

  ggplot(prof_dt, aes(x = x, y = y, color = strand)) +
    geom_path() +
    facet_grid(sample~strand) +
    ggbio::geom_arch(data = splice_dt, aes(x = start, xend = end, height = N), color = "black") +
    coord_cartesian(xlim = c(start(query_gr), end(query_gr)))
}

test_that("ssvFetchBamPE.RNA", {
  prof_dt = ssvFetchBamPE.RNA(bam_files, query_gr, return_data.table = TRUE)
  prof_dt[, x := (start + end)/2]
  testthat::expect_is(prof_dt, "data.table")
})

test_that("ssvFetchBamPE.RNA_splice", {
  splice_dt = ssvFetchBamPE.RNA_splice(bam_files, query_gr, return_data.table = TRUE)
  testthat::expect_is(splice_dt, "data.table")
  testthat::expect_true(!is.null(splice_dt$sample))
})

#QcConfigFeatures
test_that("track_rna.PE no spline", {
  p1 = track_rna.PE(bam_files, query_gr)
  p1
  testthat::expect_is(p1, "ggplot")
})

test_that("track_rna.PE other strand, default spline", {
  p2 = track_rna.PE(bam_files, query_gr, target_strand = "+")
  p2
  p2 = track_rna.PE(bam_files, query_gr, target_strand = "-")
  p2
  testthat::expect_is(p2, "ggplot")
})

test_that("track_rna.PE low win", {
  p = track_rna.PE(bam_files, query_gr, nwin = 10, nspline = 5)
  p
  testthat::expect_is(p, "ggplot")
})

test_that("track_rna.PE max win (single bp)", {
  p = track_rna.PE(bam_files, resize(query_gr, 500), nwin = Inf, nspline = 1, x_scale = "bp")
  p
  testthat::expect_is(p, "ggplot")
})

bw_cfg_dt = data.table(file = bam_files, sample = sub("_CTCF.+", "", basename(bam_files)))
bw_cfg_dt[, c("cell", "rep") := tstrsplit(sample, "[_\\.]", keep = 1:2)]
bw_cfg_dt[, sample := paste(cell, rep, sep = "\n")]

test_that("track_rna.PE by cfg", {
  p = track_rna.PE(bw_cfg_dt, query_gr, nwin = 10)+
    labs(y = "FE", title = "CTCF")
  p
  testthat::expect_is(p, "ggplot")
})

test_that("track_rna.PE by cfg - show_lines", {
  p = track_rna.PE(
    bw_cfg_dt,
    query_gr,
    nwin = 10,
    show_fill = FALSE,
    show_lines = TRUE,
    fill_VAR = "cell",
    color_VAR = "rep",
    color_mapping = c("RUNX1-ko1" = "red", "wt" = "black")
  )
  p
  testthat::expect_is(p, "ggplot")
})

test_that("track_rna.PE by cfg - no splice", {
  p = track_rna.PE(
    bw_cfg_dt,
    query_gr,
    nwin = 1000,
    show_fill = FALSE,
    show_lines = TRUE,
    show_splice = FALSE,
    fill_VAR = "cell",
    color_VAR = "rep",
    color_mapping = c("RUNX1-ko1" = "red", "wt" = "black")
  )
  p
  testthat::expect_is(p, "ggplot")
})

test_that("track_rna.PE by cfg - fill_outline_color and legend.position", {
  p = track_rna.PE(bw_cfg_dt,
                   resize(query_gr, 50*width(query_gr), fix = "center"),
                   nwin = 10,
                   show_fill = TRUE,
                   show_lines = FALSE,
                   fill_outline_color = "black", legend.position = "none")
  p
  testthat::expect_is(p, "ggplot")
})

test_that("track_rna.PE by cfg - names_on_right and facet_VAR", {
  bw_cfg_dt[, facet := "all"]
  p = track_rna.PE(bw_cfg_dt,
                   resize(query_gr, 50*width(query_gr), fix = "center"),
                   nwin = 10,
                   facet_VAR = "facet",
                   show_fill = TRUE,
                   show_lines = FALSE,
                   fill_outline_color = "black", legend.position = "none", names_on_right = FALSE) +
    labs(y = "FE", title = "CTCF")
  p
  testthat::expect_is(p, "ggplot")
})

