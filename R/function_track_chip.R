
#' track_chip
#'
#' @param signal_files
#' @param query_gr
#' @param fetch_fun
#' @param win_FUN
#' @param sum_FUN
#' @param flip_x
#' @param nwin
#' @param nspline
#' @param show_lines
#' @param show_fill
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
track_chip = function(signal_files,
                      query_gr,
                      fetch_fun = ssvFetchBam,
                      win_FUN = c("mean", "max")[2],
                      sum_FUN = NULL,
                      flip_x = NULL,
                      nwin = 3000,
                      nspline = 10,
                      show_lines = FALSE,
                      show_fill = TRUE,
                      y_label = "signal",
                      x_scale = "kbp",
                      ...){
  stopifnot(x_scale %in% c("bp", "kbp", "Mbp"))

  rng = c(start(query_gr), end(query_gr))
  if(is.null(flip_x)){
    flip_x = as.character(strand(query_gr) == "-")
  }
  if(is.null(sum_FUN)){
    sum_FUN = switch (win_FUN,
                      max = function(x, w)max(x),
                      mean = weighted.mean
    )
  }
  bw_dt.raw = fetch_fun(signal_files, resize(query_gr, 1.2*width(query_gr), fix = "center"),
                        win_method = "summary",  win_size = nwin,
                        summary_FUN = sum_FUN,
                        return_data.table = TRUE,
                        anchor = "left",
                        fragLens = NA)
  if(!is.null(bw_dt.raw$mapped_reads)){
    bw_dt.raw[, y_raw := y]
    bw_dt.raw[, y := y_raw / mapped_reads * 1e6]
    if(y_label == "signal") y_label = "RPM"
  }

  if(is.null(bw_dt.raw$cell)){
    bw_dt.raw$cell = ""
  }
  if(is.null(bw_dt.raw$mark)){
    bw_dt.raw$mark = bw_dt.raw$sample
  }

  bw_dt = bw_dt.raw[, .(y = mean(y)), .(cell, mark, x, id, start, end)]
  bw_dt[, sample := paste(cell, mark)]
  bw_dt = bw_dt[order(cell)][order(mark)]

  bw_dt$sample = factor(bw_dt$sample, levels = unique(bw_dt$sample))

  if(nspline > 1){
    bw_dt = seqsetvis::applySpline(bw_dt, n = nspline, by_ = c("sample", "id"))
    if(flip_x){
      bw_dt[, x := max(end) - (max(end) - min(start))*x]
    }else{
      bw_dt[, x := min(start) + (max(end) - min(start))*x]
    }
  }

  mark_colors = safeBrew(bw_dt$mark)
  mark_colors["input"] = "gray"

  p_chip = ggplot(bw_dt)
  if(show_lines){
    p_chip = p_chip + geom_path(aes_string(x = "x", y = "y", color = "mark"))
  }
  if(show_fill){
    p_chip = p_chip + geom_ribbon(aes_string(x = "x", ymin = 0, ymax = "y", fill = "mark")) +
      scale_fill_manual(values = mark_colors)
  }

  if(flip_x){
    rng = rev(rng)

  }

  x_label_FUN = switch (
    x_scale,
    bp = {
      function(x)x
    },
    kbp = {
      function(x)x/1e3
    },
    Mbp = {
      function(x)x/1e6
    }
  )

  p_chip = p_chip +
    guides(color = "none", fill = "none") +
    labs(x = x_scale, y = y_label) +
    facet_grid(formula("sample~.")) +
    theme_classic() +
    theme(strip.background = element_blank(), strip.text.y = element_text(angle = 0)) +
    scale_x_continuous(labels = x_label_FUN) +
    coord_cartesian(xlim = rng, expand = TRUE)
  p_chip
}
