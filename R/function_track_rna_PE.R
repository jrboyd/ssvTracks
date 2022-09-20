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
    min_splice_count = 10,
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
    target_strand = target_strand,
    flip_strand = flip_strand,
    return_data.table = TRUE
  )
  args2$bw_dt.raw = bw_dt.raw

  if(show_splice){
    splice_dt.raw = ssvFetchBamPE.RNA_splice(
      signal_files,
      query_gr,
      target_strand = target_strand,
      flip_strand = flip_strand,
      return_data.table = TRUE
    )
    # if(target_strand %in% c("+", "-")){
    #   splice_dt.raw = splice_dt.raw[strand == target_strand]
    # }
    setnames(splice_dt.raw, "N", "y")
    args2$splice_dt.raw = splice_dt.raw
  }else{
    args2["splice_dt.raw"] = list(NULL)
  }

  do.call(.track_rna_common_after_fetch, args2)
}

