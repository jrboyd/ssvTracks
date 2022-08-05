
#' track_chip
#'
#' @template ssvTracks_signal_params
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

