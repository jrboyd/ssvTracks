#' track_rna.PE
#'
#'
#' @template ssvTracks_signal_params
#'
#' @return
#' @export
#' @import seqsetvis
#'
#' @examples
#' pkg_dir = system.file(package = "ssvTracks", "extdata", mustWork = TRUE)
#' bam_files_runx1 = dir(pkg_dir, pattern = "RUNX1_RNA.+bam$", full.names = TRUE)
#' names(bam_files_runx1) = sub("_rep.+", "", basename(bam_files_runx1))
#' bed_file_runx1 = dir(pkg_dir, pattern = "RUNX1.bed", full.names = TRUE)
#' query_gr = rtracklayer::import.bed(bed_file_runx1)
#'
#' track_rna.PE(bam_files_runx1, query_gr)
#'
#' track_rna.PE(bam_files_runx1,
#'   query_gr,
#'   show_splice = TRUE,
#'   flip_strand = TRUE)
#'
#' track_rna.PE(bam_files_runx1,
#'   GRanges("chr21", IRanges(34800000, 34900000)),
#'   flip_strand = TRUE,
#'   show_splice = TRUE)
#'
#' track_rna.PE(bam_files_runx1,
#'   GRanges("chr21", IRanges(34800000, 34900000)),
#'   flip_strand = TRUE,
#'   show_splice = TRUE,
#'   color_VAR = "sample")
#'
#' track_rna.PE(bam_files_runx1,
#' query_gr,
#' show_splice = TRUE,
#' flip_strand = TRUE,
#' splice_within_range_only = resize(query_gr, 1e5)
#' )
#'
#' track_rna.PE(bam_files_runx1,
#'              query_gr,
#'              show_splice = TRUE,
#'              flip_strand = TRUE,
#'              splice_within_range_only = list(resize(query_gr, 1e5), shift(resize(query_gr, 1e5), 1e5))
#' )
track_rna.PE = function(
    signal_files,
    query_gr,
    summary_FUN = c("mean", "max")[2],
    flip_x = NULL,
    nwin = 3000,
    nspline = 1,
    fill_outline_color = NA,
    fill_alpha = 1,
    color_alpha = 1,
    y_label = "signal",
    x_scale = c("bp", "kbp", "Mbp")[2],
    floor_value = 0,
    ceiling_value = Inf,
    color_VAR = NULL,
    color_mapping = NULL,
    fill_VAR = "sample",
    fill_mapping = NULL,
    facet_VAR = "sample",
    legend.position = "right",
    names_on_right = TRUE,
    show_splice = TRUE,
    min_splice_count = 0,
    splice_within_range_only = NULL,
    target_strand = NULL,
    flip_strand = FALSE,
    return_data = FALSE,
    ...){
  env = as.list(sys.frame(sys.nframe()))
  args = c(as.list(env), list(...))
  args2 = do.call(.track_rna_common_before_fetch, args)
  for(var_name in names(args2)){
    assign(var_name, args2[[var_name]])
  }

  bw_dt.raw = ssvFetchBamPE.RNA(
    signal_files, query_gr,
    win_method = "summary",
    win_size = nwin,
    return_data.table = TRUE,
    target_strand = target_strand,
    flip_strand = flip_strand
  )
  args2$bw_dt.raw = bw_dt.raw

  if(show_splice){
    splice_dt.raw = ssvFetchBamPE.RNA_splice(
      signal_files,
      query_gr,
      return_data.table = TRUE,
      target_strand = target_strand,
      flip_strand = flip_strand
    )
    setnames(splice_dt.raw, "N", "y")
    args2$splice_dt.raw = splice_dt.raw
  }else{
    args2["splice_dt.raw"] = list(NULL)
  }

  do.call(.track_rna_common_after_fetch, args2)
}

