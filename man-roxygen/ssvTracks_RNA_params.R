
#' @template ssvTracks_signal_params
#' @param show_splice If TRUE, show splicing as arches in track. Default is
#'   TRUE.
#' @param min_splice_count Splicing events must have a value (possibly RPM
#'   normalized) greater than this value in any sample to be included in every
#'   track. Default is 0.
#' @param splice_within_range_only If single GRanges, only splice events with
#'   both ends inside GRanges will be included. If  list of length 2 GRanges,
#'   one end must intersect GRanges 1 and the other must intersect GRanges 2.
#'   Default of NULL imposes no such restriction.
