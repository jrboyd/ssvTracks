% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function_track_rna.R
\name{track_rna.SE}
\alias{track_rna.SE}
\title{track_rna.SE}
\usage{
track_rna.SE(
  signal_files,
  query_gr,
  fetch_fun = seqsetvis::ssvFetchBam,
  win_FUN = c("mean", "max")[2],
  sum_FUN = NULL,
  flip_x = NULL,
  nwin = 3000,
  nspline = 1,
  fill_outline_color = NA,
  y_label = "signal",
  x_scale = c("bp", "kbp", "Mbp")[2],
  floor_value = 0,
  ceiling_value = Inf,
  color_VAR = NULL,
  color_mapping = NULL,
  fill_VAR = "sample",
  fill_mapping = NULL,
  legend.position = "right",
  names_on_right = TRUE,
  show_splice = TRUE,
  min_splice_count = 10,
  target_strand = NULL,
  flip_strand = TRUE,
  ...
)
}
\arguments{
\item{...}{}
}
\description{
track_rna.SE
}
