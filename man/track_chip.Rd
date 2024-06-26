% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function_track_chip.R
\name{track_chip}
\alias{track_chip}
\title{track_chip}
\usage{
track_chip(
  signal_files,
  query_gr,
  fetch_fun = seqsetvis::ssvFetchBam,
  win_FUN = c("mean", "max")[2],
  sum_FUN = NULL,
  flip_x = NULL,
  nwin = 3000,
  nspline = 10,
  show_lines = FALSE,
  show_fill = TRUE,
  fill_outline_color = NA,
  y_label = "signal",
  x_scale = c("bp", "kbp", "Mbp")[2],
  floor_value = 0,
  ceiling_value = Inf,
  mark_colors = NULL,
  legend.position = "right",
  names_on_right = TRUE,
  ...
)
}
\arguments{
\item{...}{}
}
\description{
track_chip
}
