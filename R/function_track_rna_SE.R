#' track_rna.SE
#'
#' @template ssvTracks_signal_params
#'
#' @return
#' @export
#' @import seqsetvis
#'
#' @examples
track_rna.SE = function(
    signal_files,
    query_gr,
    fetch_fun = seqsetvis::ssvFetchBam,
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
    flip_strand = TRUE,
    return_data = FALSE,
    ...){
  env = as.list(sys.frame(sys.nframe()))
  args = c(as.list(env), list(...))
  args2 = do.call(.track_rna_common_before_fetch, args)
  for(var_name in names(args2)){
    assign(var_name, args2[[var_name]])
  }
  bw_dt.raw = fetch_fun(signal_files, query_gr,
                        win_method = "summary",  win_size = nwin,
                        summary_FUN = summary_FUN,
                        return_data.table = TRUE,
                        anchor = "left",
                        target_strand = target_strand,
                        flip_strand = flip_strand,
                        fragLens = NA)

  if(show_splice){
    splice_dt.raw = fetch_fun(signal_files, query_gr,
                              win_method = "summary",  win_size = nwin,
                              summary_FUN = summary_FUN,
                              return_data.table = TRUE,
                              anchor = "left",
                              fragLens = NA,
                              splice_strategy = "splice_count")
    if(flip_strand){
      splice_dt.raw[strand == "-", strand := "tmp"]
      splice_dt.raw[strand == "+", strand := "-"]
      splice_dt.raw[strand == "tmp", strand := "+"]
    }
    if(target_strand %in% c("+", "-")){
      splice_dt.raw = splice_dt.raw[strand == target_strand]
    }
    setnames(splice_dt.raw, "N", "y")
  }else{
    splice_dt.raw = NULL
  }

  args2$bw_dt.raw = bw_dt.raw
  args2$splice_dt.raw = splice_dt.raw
  do.call(.track_rna_common_after_fetch, args2)
}
