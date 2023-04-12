testthat::context("track_features")
library(ssvTracks)
library(testthat)

peak_files = dir(system.file(package = "seqsetvis", "extdata"), full.names = TRUE, pattern = "MCF.+Peak$")
peak_grs = seqsetvis::easyLoad_narrowPeak(peak_files)

olaps = seqsetvis::ssvOverlapIntervalSets(peak_grs)

query_gr = olaps[3]
query_gr = resize(query_gr, 10e4, fix = "center")

test_that("track_features defaults", {
  p = track_features(peak_grs, query_gr)
  p
  testthat::expect_is(p, "ggplot")
})

test_that("track_features other args", {
  p = track_features(peak_grs,
                     query_gr,
                     pad = .2,
                     flip_x = FALSE,
                     x_scale = "Mbp",
                     manual_levels = names(peak_grs)[c(2, 1, 3)])
  p
  testthat::expect_is(p, "ggplot")
})

test_that("track_features sample_info", {
  cfg_df = data.frame(sample = names(peak_grs))
  cfg_df$cell = sapply(strsplit(cfg_df$sample, "_"), function(x)x[1])
  cfg_df$mark = sapply(strsplit(cfg_df$sample, "_"), function(x)x[2])
  p = track_features(peak_grs,
                     query_gr,
                     sample_info_df = cfg_df,
                     sample_info_df.color_VAR = "mark",
                     sample_info_df.fill_VAR = "cell",
                     color_mapping = c("CTCF" = "black"))
  p
  testthat::expect_is(p, "ggplot")
})

test_that("track_features.numeric no VAR error", {
  testthat::expect_error(
    track_features.numeric(
      peak_grs,
      resize(olaps[3], 10e4, fix = "center"))
  )
})

test_that("track_features.numeric fill only", {
  p = track_features.numeric(
    peak_grs,
    resize(olaps[3], 10e4, fix = "center"),
    fill_VAR = "pValue")
  p
  testthat::expect_is(p, "ggplot")
})

test_that("track_features.numeric fill only new scale", {
  p = track_features.numeric(
    peak_grs,
    resize(olaps[3], 10e4, fix = "center"),
    fill_VAR = "signalValue") +
    scale_color_viridis_c() +
    scale_fill_viridis_c(option = "magma")
  p
  testthat::expect_is(p, "ggplot")
})

test_that("track_features.numeric fill and color", {
  p = track_features.numeric(
    peak_grs,
    resize(olaps[3], 10e4, fix = "center"),
    color_VAR = "pValue",
    fill_VAR = "pValue") +
    scale_color_viridis_c(option = "magma") +
    scale_fill_viridis_c(option = "magma")
  p
  testthat::expect_is(p, "ggplot")
})

test_that("track_features.numeric fill and color new scales", {
  p = track_features.numeric(
    peak_grs,
    resize(olaps[3], 10e4, fix = "center"),
    color_VAR = "pValue",
    fill_VAR = "pValue") +
    scale_color_gradientn(colors = rep("black", 2)) +
    scale_fill_viridis_c(option = "magma") +
    guides(color = "none")
  p
  testthat::expect_is(p, "ggplot")
})

test_that("track_features.numeric fill and color new scales2", {
  p = track_features.numeric(
    peak_grs,
    resize(olaps[3], 10e4, fix = "center"),
    color_VAR = "pValue",
    fill_VAR = "pValue") +
    scale_fill_gradientn(colors = rep("white", 2)) +
    scale_color_viridis_c(option = "magma") +
    guides(fill = "none")
  p
  testthat::expect_is(p, "ggplot")
})

test_that("track_features.numeric discrete color", {
  clust_dt = seqsetvis::ssvSignalClustering(seqsetvis::CTCF_in_10a_profiles_dt)
  assign_dt = unique(clust_dt[, .(id, cluster_id)])
  olap_gr = seqsetvis::CTCF_in_10a_overlaps_gr
  olap_gr$cluster_id = ""
  olap_gr[assign_dt$id]$cluster_id = assign_dt$cluster_id
  p = track_features.numeric(olap_gr,
                         query_gr = GRanges("chr1", IRanges(40e6, 70e6)),
                         color_VAR = "cluster_id",
                         fill_VAR = "cluster_id",
                         attrib = "peaks")
  p
  testthat::expect_is(p, "ggplot")

  grps = seqsetvis::ssvFactorizeMembTable(seqsetvis::CTCF_in_10a_overlaps_gr)
  plot_grs = split(olap_gr, grps$group)
  p2 = track_features.numeric(plot_grs,
                         query_gr = GRanges("chr1", IRanges(40e6, 70e6)),
                         color_VAR = "cluster_id",
                         fill_VAR = "cluster_id",
                         attrib = "peaks")
  p2
  testthat::expect_is(p2, "ggplot")
})
