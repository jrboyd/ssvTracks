library(GenomicRanges)
library(ssvTracks)
options(mc.cores = 10)
pkg_dir = system.file(package = "ssvTracks", "extdata", mustWork = TRUE)
bam_files_runx1 = dir(pkg_dir, pattern = "RUNX1_RNA.+bam$", full.names = TRUE)
names(bam_files_runx1) = sub("_rep.+", "", basename(bam_files_runx1))
bed_file_runx1 = dir(pkg_dir, pattern = "RUNX1.bed", full.names = TRUE)
query_gr = rtracklayer::import.bed(bed_file_runx1)
strand(query_gr) = "-"
strand(query_gr) = "+"

cfg_dt = data.table(file = bam_files_runx1, name = names(bam_files_runx1))
track_rna.PE(cfg_dt,
             query_gr,
             show_splice = TRUE,
             flip_strand = TRUE, fill_VAR = NULL, color_VAR = "strand", facet_VAR = "name", target_strand = "both")

track_rna.PE(cfg_dt,
             GRanges("chr21", IRanges(34850e3, 34885e3)),
             show_splice = TRUE,
             flip_strand = TRUE, fill_VAR = NULL, color_VAR = "strand", facet_VAR = "name", target_strand = "both")

track_rna.PE(cfg_dt,
             GRanges("chr21", IRanges(34859.4e3, 34859.6e3)),
             show_splice = TRUE,
             flip_strand = TRUE, fill_VAR = NULL, color_VAR = "strand", facet_VAR = "name", target_strand = "both")

track_rna.PE(cfg_dt,
             query_gr,
             show_splice = TRUE,
             flip_strand = TRUE, fill_VAR = NULL, color_VAR = "name", facet_VAR = "strand", target_strand = "both")


track_rna.PE(cfg_dt,
             query_gr,
             show_splice = TRUE,
             flip_strand = FALSE, fill_VAR = NULL, color_VAR = "strand", facet_VAR = "name", target_strand = "both")

track_rna.PE(cfg_dt,
             query_gr,
             show_splice = TRUE,
             flip_strand = FALSE, fill_VAR = NULL, color_VAR = "strand", facet_VAR = "name")

track_rna.PE(cfg_dt,
             query_gr,
             show_splice = TRUE,
             flip_strand = TRUE, fill_VAR = NULL, color_VAR = "strand", facet_VAR = "name")



pkg_dir = system.file(package = "ssvTracks", "extdata", mustWork = TRUE)
bam_files_esr1 = dir(pkg_dir, pattern = "M.+R1.ESR1_RNA.+bam$", full.names = TRUE)
names(bam_files_esr1) = sub("_R.+", "", basename(bam_files_esr1))
bed_file_esr1 = dir(pkg_dir, pattern = "ESR1.bed", full.names = TRUE)
query_gr = rtracklayer::import.bed(bed_file_esr1)
strand(query_gr) = "-"
strand(query_gr) = "+"

cfg_dt = data.table(file = bam_files_esr1, name = names(bam_files_esr1))
track_rna.SE(cfg_dt,
  query_gr,
  show_splice = TRUE,
  flip_strand = TRUE, fill_VAR = NULL, color_VAR = "strand", facet_VAR = "name", target_strand = "both")

track_rna.SE(cfg_dt,
             query_gr,
             show_splice = TRUE,
             flip_strand = TRUE, fill_VAR = NULL, color_VAR = "name", facet_VAR = "strand", target_strand = "both")


track_rna.SE(cfg_dt,
             GRanges("chr6", IRanges(151800e3, 151850e3)),
             show_splice = TRUE,
             flip_strand = TRUE, fill_VAR = NULL, color_VAR = "name", facet_VAR = "strand", target_strand = "both")

track_rna.SE(cfg_dt,
             GRanges("chr6", IRanges(151807e3, 151809e3)),
             show_splice = TRUE,
             flip_strand = TRUE, fill_VAR = NULL, color_VAR = "name", facet_VAR = "strand", target_strand = "both")

track_rna.SE(cfg_dt,
             GRanges("chr6", IRanges(151807e3, 151809e3)),
             show_splice = TRUE,
             flip_strand = TRUE, fill_VAR = NULL, color_VAR = "name", facet_VAR = "strand", target_strand = "+")

#
track_rna.SE(cfg_dt,
             GRanges("chr6", IRanges(151807e3, 151809e3)),
             show_pileup = FALSE,
             show_splice = TRUE, splice_within_range_only = GRanges("chr6", IRanges(151807e3, 151809e3)),
             flip_strand = TRUE, fill_VAR = NULL, color_VAR = "name", facet_VAR = "strand", target_strand = "+")

track_rna.SE(cfg_dt,
             GRanges("chr6", IRanges(151807e3, 151809e3)),
             show_pileup = FALSE,
             show_splice = TRUE,
             flip_strand = TRUE, fill_VAR = "name", color_VAR = "name", facet_VAR = "name", target_strand = "+")



track_rna.SE(cfg_dt,
             GRanges("chr6", IRanges(151807e3, 151809e3)),
             show_splice = TRUE,
             flip_strand = TRUE, fill_VAR = NULL, color_VAR = "name", facet_VAR = "strand", target_strand = "-")
