
#' track_chip
#'
#' @param signal_files Character or data.frame. Paths to files. Either supplied as a simple character vector (names are use for plot labels if set) or as a data.frame
#' @param query_gr GRanges. Defines regions to be fetched and plotted.
#' @param fetch_fun An ssvFetc* function from seqsetvis. ssvFetchBam, ssvFetchBamPE or ssvFetchBigwig are likely choices.
#' @param summary_FUN Either a function or character "mean" or "max".  If a custom function it must follow the form of weighted.mean() and accept 2 arguments, values and weights.
#' @param flip_x If TRUE, x-axis is flipped to be decreasing.  Default of NULL will flip_x automatically when query_gr strand is negative.
#' @param nwin Numeric. Higher numbers increase resolution but increase plotting time and size of .pdf files.
#' @param nspline Numeric. If higher than 1, splines will be used to interpolate and smooth between windows.
#' @param fill_outline_color Character. Color applied to outline for geom_ribbon used for filled tracks.
#' @param y_label Character. Label for y-axis.
#' @param x_scale One of "bp", "kbp", or "Mbp".  Scales x-axis labels accordingly.
#' @param floor_value Numeric.  Values below floor will be increased to floor. Default is 0.
#' @param ceiling_value Numeric. Values above ceiling will be decreased to ceiling. Default is Inf.
#' @param color_VAR Character. Color variable if supplying signal_files as data.frame. Default of NULL disables colored line.
#' @param color_mapping Named character. Maps values of color_VAR to valid colors.
#' @param fill_VAR Character. Fill variable if supplying signal_files as data.frame. Default of "sample" results in 1 fill color per file.
#' @param fill_mapping Named character. Maps values of fill_VAR to valid colors.
#' @param facet_VAR Character.  Files that share a facet value will appear in the same track row. Default of "sample" results in 1 file per row.
#' @param legend.position Charactter. Position for legend, see ggplot2::theme.  Most likely "right", "bottom", or "none".
#' @param names_on_right logical. If TRUE (default) facet/row names appear on the right. If FALSE facet/row names appear on the left.
#' @param ... Currently not used.
#'
#' @return
#' @export
#' @import seqsetvis
#'
#' @examples
#'
#' bam_files = dir(system.file(package = "ssvTracks", "extdata"), full.names = TRUE, pattern = "100peaks.bam$")
#' bw_files = dir(system.file(package = "seqsetvis", "extdata"), full.names = TRUE, pattern = "MCF.+bw$")
#' peak_files = dir(system.file(package = "seqsetvis", "extdata"), full.names = TRUE, pattern = "MCF.+Peak$")
#' peak_grs = seqsetvis::easyLoad_narrowPeak(peak_files)
#' query_gr = peak_grs[[1]][1]
#'
#' track_chip(bam_files, query_gr)
#'
#' track_chip(
#'   bam_files,
#'   query_gr,
#'   floor_value = 0,
#'   ceiling_value = 5
#' )
#'
#' bam_cfg = data.frame(bam_files)
#' track_chip(bam_cfg, query_gr)
#'
#' bam_cfg = data.frame(not_file = bam_files)
#' track_chip(bam_cfg, query_gr)
#'
#' bam_cfg = data.frame(bam_files, name = basename(bam_files))
#' track_chip(bam_cfg, query_gr, fill_VAR = "name", facet_VAR = "name")
#'
#' bam_cfg = data.frame(bam_files, name = basename(bam_files))
#' #use separate to parse out metadata
#' bam_cfg = bam_cfg %>%
#'   separate(.,
#'   col = "name",
#'   sep = "_",
#'   into = c("cell", "mark"),
#'   remove = FALSE,
#'   extra = "drop")
#'
#' track_chip(
#'   bam_cfg,
#'   query_gr,
#'   fill_VAR = "mark",
#'   color_VAR = "mark",
#'   facet_VAR = "cell",
#'   nspline = 1,
#'   x_scale = "bp",
#'   fill_alpha = .2
#' )
#'
#' track_chip(bw_files, query_gr, target_strand = "*", fetch_fun = seqsetvis::ssvFetchBigwig, nspline = 1)
#'
#' bam_files.PE = dir(system.file(package = "ssvTracks", "extdata"), full.names = TRUE, pattern = "PE.+bam$")
#' bed_gr = rtracklayer::import.bed(system.file(package = "ssvTracks", "extdata/ESR1.bed"))
#' track_chip(bam_files.PE, bed_gr, target_strand = "*", fetch_fun = seqsetvis::ssvFetchBamPE)
#'
track_chip = function(signal_files,
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
                      ...){
  env = as.list(sys.frame(sys.nframe()))
  args = c(as.list(env), list(...))
  args2 = do.call(.track_all_common_before_fetch, args)
  for(var_name in names(args2)){
    assign(var_name, args2[[var_name]])
  }

  bw_dt.raw = fetch_fun(signal_files,
                        query_gr,
                        win_method = "summary",
                        win_size = nwin,
                        summary_FUN = summary_FUN,
                        return_data.table = TRUE,
                        anchor = "left",
                        fragLens = NA
  )

  args2$bw_dt.raw = bw_dt.raw

  do.call(.track_all_common_after_fetch, args2)
}

