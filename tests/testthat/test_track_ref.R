testthat::context("track_gene_transcripts")
library(ssvTracks)

bam_files = dir(system.file(package = "ssvTracks", "extdata"), full.names = TRUE, pattern = "bam$")
bw_files = dir(system.file(package = "seqsetvis", "extdata"), full.names = TRUE, pattern = "MCF.+bw$")
peak_files = dir(system.file(package = "seqsetvis", "extdata"), full.names = TRUE, pattern = "MCF.+Peak$")
peak_grs = seqsetvis::easyLoad_narrowPeak(peak_files)

gtf_file = dir(system.file(package = "ssvTracks", "extdata"), full.names = TRUE, pattern = "at_peaks.gtf$")

# hit_goi = unique(ex_gr$gene_name)
# ref_gr.all = rtracklayer::import.gff("~/gencode.v36.annotation.gtf")
# ref_gr = subset(ref_gr.all, gene_name %in% hit_goi)
#
# rtracklayer::export.gff(ref_gr, file.path("inst/extdata/", basename(gtf_file)))

ex_gr = rtracklayer::import.gff(gtf_file, feature.type = "exon")

query_gr = subsetByOverlaps(peak_grs[[1]], ex_gr)[7]
query_gr = resize(query_gr, 10e4, fix = "center")

if(FALSE){
  subsetByOverlaps(ex_gr, query_gr)

  cowplot::plot_grid(
    nrow = 1,
    track_gene_reference(ex_gr,
                         query_gr = query_gr,
                         show_tss = TRUE) +
      labs(title = "no flip"),
    track_gene_reference(ex_gr,
                         query_gr = query_gr,
                         show_tss = TRUE,
                         flip_x = TRUE) +
      labs(title = "flip")
  )

  track_gene_transcripts(ex_gr, sel_gene_name = "GALNT12")
  track_gene_transcripts(ex_gr, sel_gene_name = "GALNT12", flip_x = TRUE)

  track_gene_transcripts(ex_gr, sel_gene_name = "DZANK1")
  track_gene_transcripts(ex_gr, sel_gene_name = "DZANK1", flip_x = TRUE)
  track_gene_transcripts(ex_gr, sel_gene_name = "DZANK1", flip_x = FALSE)
}

test_that("track_gene_reference defaults", {
  p = track_gene_reference(ex_gr,
                           query_gr = query_gr)
  testthat::expect_is(p, "ggplot")
})

#QcConfigFeatures
test_that("track_gene_reference all args", {
  p = track_gene_reference(ex_gr,
                           query_gr = query_gr,
                           flip_x = TRUE,
                           exon_height = .3,
                           intron_thickness = .1,
                           tss_size = 2,
                           tss_color = "green",
                           tss_up = .4,
                           tss_over = .02,
                           tss_arrow_size = 2,
                           x_scale = "kbp",
                           minus_strand_color = "black",
                           plus_strand_color = "black",
                           legend.position = "none",
                           show_tss = TRUE)
  testthat::expect_is(p, "ggplot")
})

test_that("track_gene_reference no hit", {
  p = track_gene_reference(ex_gr,
                           query_gr = shift(query_gr, 2e5),
  )
  testthat::expect_is(p, "ggplot")
})

test_that("track_gene_reference return_data", {
  res = track_gene_reference(ex_gr,
                             query_gr = query_gr,
                             return_data = TRUE
  )
  testthat::expect_is(res$ref_dt, "data.table")
  testthat::expect_is(res$tss_dt, "data.table")
  testthat::expect_is(res$label_dt, "data.table")
})
