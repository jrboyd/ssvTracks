
#' my_plot_goi
#'
#' @param goi
#' @param manual_query_gr
#' @param rel_heights.plots
#' @param rel_heights.title
#' @param peak_grs
#' @param nwin
#' @param nspline
#' @param my_query_dt
#' @param show_all_x_axis
#'
#' @return
#' @export
#'
#' @examples
my_plot_goi = function(goi,
                       manual_query_gr = NULL,
                       rel_heights.plots = c(3, 1, 3),
                       rel_heights.title = c(1, 10),
                       peak_grs = NULL,
                       nwin = 3000,
                       nspline = 10,
                       my_query_dt = query_dt,
                       show_all_x_axis = FALSE){
  if(is.null(manual_query_gr)){
    qgr = range(subset(ex_gr, gene_name == goi))
    qgr = resize(qgr, 1.3*width(qgr), fix = "center")
    if(width(qgr) < 20e3){
      qgr = resize(qgr, 20e3, fix = "center")
    }
  }else{
    qgr = manual_query_gr
  }


  show_lines = FALSE
  show_fill = TRUE
  flip_x = as.character(strand(qgr) == "-")
  rng = c(start(qgr), end(qgr))

  p_ref = track_ref(ex_gr, qgr, show_tss = TRUE, flip_x = flip_x)

  ####

  p_chip = track_chip(my_query_dt,
                      qgr,
                      fetch_fun = ssvFetchBam,
                      win_FUN = "max",
                      nwin = nwin,
                      nspline = nspline) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3))
  if(!show_all_x_axis){
    p_chip = p_chip +
      theme(axis.line.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()) +
      labs(x = "")
  }

  if(!is.null(peak_grs)){
    p_gr = track_gr(peak_grs, qgr)
    pg_ass = track_assembly(list(p_chip, p_ref, p_gr), qgr, rel_heights = rel_heights.plots)
  }else{
    pg_ass = track_assembly(list(p_chip, p_ref), qgr, rel_heights = rel_heights.plots[1:2])
  }

  cowplot::plot_grid(ncol = 1, rel_heights = rel_heights.title,
                     ggplot() +
                       labs(title = goi) +
                       theme_void() +
                       theme(plot.title = element_text(hjust = .5, vjust = .5)),
                     pg_ass
  )
}

#### Track Helper Functions Assembly ####

