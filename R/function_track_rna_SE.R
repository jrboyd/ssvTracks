#' track_rna.SE
#'
#' @param signal_files
#' @param query_gr
#' @param fetch_fun
#' @param win_FUN
#' @param sum_FUN
#' @param flip_x
#' @param nwin
#' @param nspline
#' @param fill_outline_color
#' @param fill_alpha
#' @param color_alpha
#' @param y_label
#' @param x_scale
#' @param floor_value
#' @param ceiling_value
#' @param color_VAR
#' @param color_mapping
#' @param fill_VAR
#' @param fill_mapping
#' @param facet_VAR
#' @param legend.position
#' @param names_on_right
#' @param show_splice
#' @param min_splice_count
#' @param target_strand
#' @param flip_strand logical. Should aligned strand be flipped relative to query_gr.  The vast majority of SE libraries are actually flip_strand, so default flip_strand is TRUE.
#' @param return_data
#' @param ...
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
    win_FUN = c("mean", "max")[2],
    sum_FUN = NULL,
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
                        summary_FUN = sum_FUN,
                        return_data.table = TRUE,
                        anchor = "left",
                        target_strand = target_strand,
                        flip_strand = flip_strand,
                        fragLens = NA)

  if(show_splice){
    splice_dt.raw = fetch_fun(signal_files, query_gr,
                              win_method = "summary",  win_size = nwin,
                              summary_FUN = sum_FUN,
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
