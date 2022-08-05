
#' @param signal_files Character or data.frame. Paths to files. Either supplied as a simple character vector (names are use for plot labels if set) or as a data.frame
#' @template ssvTracks_common_params
#' @param fetch_fun An ssvFetc* function from seqsetvis. ssvFetchBam, ssvFetchBamPE or ssvFetchBigwig are likely choices.
#' @param summary_FUN Either a function or character "mean" or "max".  If a custom function it must follow the form of weighted.mean() and accept 2 arguments, values and weights.
#' @param nwin Numeric. Higher numbers increase resolution but increase plotting time and size of .pdf files.
#' @param nspline Numeric. If higher than 1, splines will be used to interpolate and smooth between windows.
#' @param fill_outline_color Character. Color applied to outline for geom_ribbon used for filled tracks.
#' @param y_label Character. Label for y-axis.
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
