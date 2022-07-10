testthat::context("track_ref")
library(ssvTracks)

bam_files = dir(system.file(package = "ssvTracks", "extdata"), full.names = TRUE, pattern = "bam$")
bw_files = dir(system.file(package = "seqsetvis", "extdata"), full.names = TRUE, pattern = "MCF.+bw$")
peak_files = dir(system.file(package = "seqsetvis", "extdata"), full.names = TRUE, pattern = "MCF.+Peak$")
peak_grs = seqsetvis::easyLoad_narrowPeak(peak_files)
query_gr = peak_grs[[1]][1]

gtf_file = dir(system.file(package = "ssvTracks", "extdata"), full.names = TRUE, pattern = "at_peaks.gtf$")


hit_goi = unique(ex_gr$gene_name)

# ref_gr.all = rtracklayer::import.gff("~/gencode.v36.annotation.gtf")
# ref_gr = subset(ref_gr.all, gene_name %in% hit_goi)
#
# rtracklayer::export.gff(ref_gr, gtf_file)

query_gr = subsetByOverlaps(peak_grs[[1]], ex_gr)[7]
query_gr = resize(query_gr, 10e4, fix = "center")
ex_gr = rtracklayer::import.gff(gtf_file, feature.type = "exon")

track_ref(ex_gr, query_gr = query_gr)


cowplot::plot_grid(
  nrow = 1,
  track_ref(ex_gr,
            query_gr = query_gr,
            show_tss = TRUE) +
    labs(title = "no flip"),
  track_ref(ex_gr,
            query_gr = query_gr,
            show_tss = TRUE,
            flip_x = TRUE) +
    labs(title = "flip")
)


track_ref(ex_gr,
          query_gr = query_gr,
          minus_strand_color = "black",
          plus_strand_color = "black",
          legend.position = "none",
          show_tss = TRUE)

#QcConfigFeatures
test_that("track_chip no spline", {
  p = track_chip(fetch_files, query_gr, fetch_fun = fetch_fun, nspline = 1)
  testthat::expect_is(p, "ggplot")
})
