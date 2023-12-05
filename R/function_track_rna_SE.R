#' track_rna.SE
#'
#' @template ssvTracks_signal_params
#'
#' @return
#' @export
#' @import seqsetvis
#'
#' @examples
#' pkg_dir = system.file(package = "ssvTracks", "extdata", mustWork = TRUE)
#' bam_files_esr1 = dir(pkg_dir, pattern = "M.+R1.ESR1_RNA.+bam$", full.names = TRUE)
#' names(bam_files_esr1) = sub("_R.+", "", basename(bam_files_esr1))
#' bed_file_esr1 = dir(pkg_dir, pattern = "ESR1.bed", full.names = TRUE)
#' query_gr = rtracklayer::import.bed(bed_file_esr1)
#'
#' track_rna.SE(bam_files_esr1, query_gr)
#'
#' track_rna.SE(bam_files_esr1,
#'   query_gr,
#'   show_splice = TRUE,
#'   flip_strand = TRUE)
#'
#' track_rna.SE(bam_files_esr1,
#'   query_gr,
#'   show_splice = TRUE,
#'   flip_strand = FALSE)
#'
#' track_rna.SE(bam_files_esr1,
#'   GRanges("chr6", IRanges(151800000, 151900000)),
#'   flip_strand = TRUE,
#'   show_splice = TRUE)
#'
#' track_rna.SE(bam_files_esr1,
#'   GRanges("chr6", IRanges(151800000, 151900000)),
#'   flip_strand = TRUE,
#'   show_splice = TRUE,
#'   color_VAR = "sample")
track_rna.SE = function(
    signal_files,
    query_gr,
    fetch_fun = seqsetvis::ssvFetchBam,
    fetch_fun_splice = fetch_fun,
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
    color_over_fill = TRUE,
    facet_VAR = "sample",
    legend.position = "right",
    names_on_right = TRUE,
    show_splice = TRUE,
    show_pileup = TRUE,
    min_splice_count = 0,
    splice_within_range_only = NULL,
    target_strand = NULL,
    flip_strand = TRUE,
    flip_strand_splice = flip_strand,
    return_data = FALSE,
    ...){
  env = as.list(sys.frame(sys.nframe()))
  args = c(as.list(env), list(...))
  args2 = do.call(.track_rna_common_before_fetch, args)
  for(var_name in names(args2)){
    assign(var_name, args2[[var_name]])
  }
  strand(query_gr) = "*" #seqsetvis handles query strand in a relative fashion (appropriate for ChIP), we need absolute strand
  bw_dt.raw = fetch_fun(signal_files, query_gr,
                        win_method = "summary",  win_size = nwin,
                        summary_FUN = summary_FUN,
                        return_data.table = TRUE,
                        anchor = "left",
                        target_strand = target_strand,
                        flip_strand = flip_strand,
                        fragLens = NA,
                        splice_strategy = "ignore")

  if(show_splice){
    splice_dt.raw = fetch_fun_splice(signal_files, query_gr,
                              win_method = "summary",  win_size = nwin,
                              summary_FUN = summary_FUN,
                              return_data.table = TRUE,
                              anchor = "left",
                              target_strand = target_strand,
                              flip_strand = flip_strand_splice,
                              fragLens = NA,
                              splice_strategy = "splice_count")
    setnames(splice_dt.raw, "N", "y")
  }else{
    splice_dt.raw = NULL
  }

  args2["bw_dt.raw"] = list(bw_dt.raw)
  args2["splice_dt.raw"] = list(splice_dt.raw)
  do.call(.track_rna_common_after_fetch, args2)
}
