testthat::context("track_chip")
library(ssvTracks)

bam_files = dir(system.file(package = "ssvTracks", "extdata"), full.names = TRUE, pattern = "bam$")
bw_files = dir(system.file(package = "seqsetvis", "extdata"), full.names = TRUE, pattern = "MCF.+bw$")
peak_files = dir(system.file(package = "seqsetvis", "extdata"), full.names = TRUE, pattern = "MCF.+Peak$")
peak_grs = seqsetvis::easyLoad_narrowPeak(peak_files)
query_gr = peak_grs[[1]][1]

gtf_file = dir(system.file(package = "ssvTracks", "extdata"), full.names = TRUE, pattern = "at_peaks.gtf$")

ex_gr = rtracklayer::import.gff(gtf_file, feature.type = "exon")

#QcConfigFeatures
test_that("QcConfigFeatures.parse", {
  track_chip(bw_files, query_gr, nwin = 100, fetch_fun = seqsetvis::ssvFetchBigwig)
})
