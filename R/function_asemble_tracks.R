#' assemble_tracks
#'
#' @param plot_list list of ggplots and/or grobs
#' @param query_gr
#' @param rel_heights
#'
#' @return
#' @export
#'
#' @examples
assemble_tracks = function(plot_list, query_gr, rel_heights = rep(1, length(plot_list))){
  stopifnot(length(query_gr) == 1)
  # set x-limits on all plots
  xlim = c(start(query_gr), end(query_gr))
  if(as.character(strand(query_gr)) == "-") xlim = rev(xlim)
  for(i in seq(1, length(plot_list))){
    plot_list[[i]] = plot_list[[i]] +
      labs(x = "") +
      coord_cartesian(xlim = xlim)
  }
  # remove x-axis labels from all plots but final
  if(length(plot_list) > 1){
    for(i in seq(1, length(plot_list) - 1)){
      plot_list[[i]] = plot_list[[i]] + labs(x = "")
    }
  }
  # add axis labels to final plot
  plot_list[[length(plot_list)]] = plot_list[[length(plot_list)]] + labs(x = paste(as.character(seqnames(query_gr)), "kbp")) + scale_x_continuous(labels = function(x)x/1e3)
  plot_list = sync_width(plot_list)
  pg = cowplot::plot_grid(plotlist = plot_list, ncol = 1,
                          rel_heights = rel_heights, scale = 1)
  pg
}



#' sync_width
#'
#' @param grob_list list containing ggplot and/or grobs
#'
#' @return
#' @export
#'
#' @examples
sync_width = function(grob_list){
  stopifnot(class(grob_list) == "list")
  is_ok = sapply(grob_list, function(x){
    "ggplot" %in% class(x) | "grob" %in% class(x)
  })
  stopifnot(all(is_ok))
  my_grobs = lapply(grob_list, function(x){
    if(grid::is.grob(x)){
      x
    }else{
      ggplotGrob(x)
    }
  })

  my_widths = lapply(my_grobs, function(gt){
    gt$widths
  })
  maxWidth = my_widths[[1]]
  if(length(my_widths) > 1){
    for(i in 2:length(my_widths)){
      maxWidth = grid::unit.pmax(maxWidth, my_widths[[i]])
    }
  }
  for(j in 1:length(my_grobs)){
    my_grobs[[j]]$widths = maxWidth
  }
  my_grobs
}

#' sync_height
#'
#' @param grob_list list containing ggplot and/or grobs
#'
#' @return
#' @export
#'
#' @examples
sync_height = function(grob_list){
  stopifnot(class(grob_list) == "list")
  is_ok = sapply(grob_list, function(x){
    "ggplot" %in% class(x) | "grob" %in% class(x)
  })
  stopifnot(all(is_ok))
  my_grobs = lapply(grob_list, function(x){
    if(grid::is.grob(x)){
      x
    }else{
      ggplotGrob(x)
    }
  })

  my_widths = lapply(my_grobs, function(gt){
    gt$heights
  })
  maxWidth = my_widths[[1]]
  if(length(my_widths) > 1){
    for(i in 2:length(my_widths)){
      maxWidth = grid::unit.pmax(maxWidth, my_widths[[i]])
    }
  }
  for(j in 1:length(my_grobs)){
    my_grobs[[j]]$heights = maxWidth
  }
  my_grobs
}
