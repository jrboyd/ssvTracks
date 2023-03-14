
#' track_features
#'
#'
#' @param feature_grs List. Must contain GRanges.  To include meta data as color or fill, specify sample_info_df.
#' @param query_gr
#' @param sample_info_df
#' @param sample_info_df.name_VAR Character. Default is "sample". Use for mapping sample_info to feature_grs. Values of sample_info_df for this attribute must be setequal to names of feature_grs.
#' @param sample_info_df.color_VAR Character. Default is NULL (no color mapping). Use for mapping sample_info to color. Values of sample_info_df for this attribute must be setequal to names of color_mapping.
#' @param sample_info_df.fill_VAR Character. Default is "sample". Use for mapping sample_info to fill Values of sample_info_df for this attribute must be setequal to names of fill_mapping.
#' @param attrib
#' @param pad
#' @param flip_x
#' @param manual_levels
#' @param x_scale
#' @param color_mapping
#' @param fill_mapping
#' @param legend.position
#'
#' @return
#' @export
#'
#' @examples
#' peak_grs = seqsetvis::easyLoad_narrowPeak(peak_files)
#' olaps = seqsetvis::ssvOverlapIntervalSets(peak_grs)
#' query_gr = olaps[1]
#' query_gr = resize(query_gr, 10e4, fix = "center")
#' feature_grs = peak_grs
#'
#'
#'
#'
#'
#' track_features(peak_grs, query_gr)
#'
#' peak_grs.is_sig = lapply(peak_grs, function(x){
#'   x$is_sig = x$pValue > 20
#'   x
#' })
#'
#' track_features(peak_grs.is_sig, query_gr, attrib = "is_sig")
#'
track_features = function(feature_grs,
                          query_gr,
                          sample_info_df = NULL,
                          sample_info_df.name_VAR = "sample",
                          sample_info_df.color_VAR = NULL,
                          sample_info_df.fill_VAR = "sample",
                          pad = .1,
                          flip_x = NULL,
                          manual_levels = NULL,
                          x_scale = c("bp", "kbp", "Mbp")[2],
                          color_mapping = NULL,
                          fill_mapping = NULL,
                          legend.position = "right"){
  .check_query_gr(query_gr)
  if(is(feature_grs, "GRangesList")) feature_grs = as.list(feature_grs)
  if(is.list(feature_grs)){
    for(i in seq_along(feature_grs)){
      mcols(feature_grs[[i]])[[sample_info_df.name_VAR]] = names(feature_grs)[i]
    }
    if(is.null(manual_levels)) manual_levels = rev(names(feature_grs))
    feature_grs = unlist(GRangesList(feature_grs))
  }
  feature_grs.hit = subsetByOverlaps(feature_grs, query_gr, ignore.strand = TRUE)
  if(length(feature_grs.hit) == 0){
    warning("no features found in query_gr.")
  }
  if(!is.factor(mcols(feature_grs.hit)[[sample_info_df.name_VAR]])){
    mcols(feature_grs.hit)[[sample_info_df.name_VAR]] = factor(mcols(feature_grs.hit)[[sample_info_df.name_VAR]])
  }else{
    mcols(feature_grs.hit)[[sample_info_df.name_VAR]] = factor(mcols(feature_grs.hit)[[sample_info_df.name_VAR]], levels = rev(levels(mcols(feature_grs.hit)[[sample_info_df.name_VAR]])))
  }
  c_id = mcols(feature_grs.hit)[[sample_info_df.name_VAR]]
  mcols(feature_grs.hit) = NULL
  names(feature_grs.hit) = NULL
  mcols(feature_grs.hit)[[sample_info_df.name_VAR]] = c_id
  # feature_grs.hit = unique(feature_grs.hit)
  feature_grs.hit = as.data.table(feature_grs.hit)

  if(!is.null(manual_levels)){
    stopifnot(feature_grs.hit[[sample_info_df.name_VAR]] %in% manual_levels)
    set(feature_grs.hit, j = sample_info_df.name_VAR, value = factor(feature_grs.hit[[sample_info_df.name_VAR]], levels = manual_levels))
  }

  feature_grs.hit[, ymin := as.numeric(get(sample_info_df.name_VAR))+pad-1]
  feature_grs.hit[, ymax := as.numeric(get(sample_info_df.name_VAR))-pad]

  lev = levels(feature_grs.hit[[sample_info_df.name_VAR]])

  if(!is.null(sample_info_df)){
    if(is.null(sample_info_df[[sample_info_df.name_VAR]])){
      sample_info_df.name_VAR = colnames(sample_info_df)[1]
    }
  }else{
    sample_info_df = data.frame(lev)
    colnames(sample_info_df) = sample_info_df.name_VAR
    if(!is.null(sample_info_df.fill_VAR)) sample_info_df[[sample_info_df.fill_VAR]] = sample_info_df[[sample_info_df.name_VAR]]
    if(!is.null(sample_info_df.color_VAR)) sample_info_df[[sample_info_df.color_VAR]] = sample_info_df[[sample_info_df.name_VAR]]
  }
  sample_info_df[[sample_info_df.name_VAR]] = sample_info_df[[sample_info_df.name_VAR]]
  feature_grs.hit = merge(feature_grs.hit, sample_info_df, by = sample_info_df.name_VAR)
  if(nrow(feature_grs.hit) == 0){
    warning("all features lost merging with sample_info_df. check values of sample_info_df.name_VAR ", sample_info_df.name_VAR)
  }

  if(is.null(color_mapping)){
    if(is.null(sample_info_df.color_VAR)){
      color_mapping = NULL
    }else{
      color_mapping = seqsetvis::safeBrew(feature_grs.hit[[sample_info_df.color_VAR]])
    }
  }
  if(is.null(fill_mapping)){
    if(is.null(sample_info_df.fill_VAR)){
      fill_mapping = NULL
    }else{
      fill_mapping = seqsetvis::safeBrew(feature_grs.hit[[sample_info_df.fill_VAR]])
    }
  }

  p_gr = ggplot(feature_grs.hit,
                aes_string(
                  xmin = "start",
                  xmax = "end",
                  ymin = "ymin",
                  ymax = "ymax",
                  fill = sample_info_df.fill_VAR,
                  color = sample_info_df.color_VAR))

  if(is.null(sample_info_df.color_VAR)){
    p_gr = p_gr + geom_rect(color = NA)
  }else if(is.null(sample_info_df.fill_VAR)){
    p_gr = p_gr + geom_rect(fill = NA)
  }else{
    p_gr = p_gr + geom_rect()
  }
  p_gr = p_gr +
    scale_y_continuous(breaks = seq_along(lev)-.5, labels = lev, limits = c(0, length(lev))) +
    theme_classic() +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          legend.position = legend.position) +
    scale_color_manual(values = color_mapping) +
    scale_fill_manual(values = fill_mapping)
  p_gr = .apply_x_scale(p_gr, x_scale = x_scale, as.character(seqnames(query_gr)))
  p_gr = .apply_x_lim(p_gr, query_gr, flip_x = flip_x)
  p_gr
}

#' track_features.numeric
#'
#'
#' @param feature_grs
#' @param query_gr
#' @param color_VAR
#' @param fill_VAR
#' @param attrib
#' @param pad
#' @param flip_x
#' @param manual_levels
#' @param x_scale
#' @param color_mapping
#' @param fill_mapping
#' @param legend.position
#'
#' @return
#' @export
#'
#' @examples
track_features.numeric = function(feature_grs,
                                  query_gr,
                                  color_VAR = NULL,
                                  fill_VAR = NULL,
                                  attrib = "cluster_id",
                                  pad = .1,
                                  flip_x = NULL,
                                  manual_levels = NULL,
                                  x_scale = c("bp", "kbp", "Mbp")[2],
                                  color_mapping = NULL,
                                  fill_mapping = NULL,
                                  legend.position = "right"){
  .check_query_gr(query_gr)
  if(is.null(color_VAR) & is.null(fill_VAR)) stop("One of color_VAR or fill_VAR must be set.")

  if(is(feature_grs, "GRangesList")) feature_grs = as.list(feature_grs)
  if(is.list(feature_grs)){
    for(i in seq_along(feature_grs)){
      mcols(feature_grs[[i]])[[attrib]] = names(feature_grs)[i]
    }
    if(is.null(manual_levels)) manual_levels = rev(names(feature_grs))
    feature_grs = unlist(GRangesList(feature_grs))
  }
  feature_grs.hit = subsetByOverlaps(feature_grs, query_gr, ignore.strand = TRUE)
  if(!is.factor(mcols(feature_grs.hit)[[attrib]])){
    mcols(feature_grs.hit)[[attrib]] = factor(mcols(feature_grs.hit)[[attrib]])
  }else{
    mcols(feature_grs.hit)[[attrib]] = factor(mcols(feature_grs.hit)[[attrib]], levels = rev(levels(mcols(feature_grs.hit)[[attrib]])))
  }
  # c_id = mcols(feature_grs.hit)[[attrib]]
  # mcols(feature_grs.hit) = NULL
  # names(feature_grs.hit) = NULL
  # mcols(feature_grs.hit)[[attrib]] = c_id
  # feature_grs.hit = unique(feature_grs.hit)
  feature_grs.hit = as.data.table(feature_grs.hit)

  if(!is.null(color_VAR)) if(is.null(feature_grs.hit[[color_VAR]])){
    stop("color_VAR ", color_VAR, " must be present in input GRanges")
  }
  if(!is.null(fill_VAR)) if(is.null(feature_grs.hit[[fill_VAR]])){
    stop("fill_VAR ", fill_VAR, " must be present in input GRanges")
  }

  if(!is.null(manual_levels)){
    stopifnot(feature_grs.hit[[attrib]] %in% manual_levels)
    set(feature_grs.hit, j = attrib, value = factor(feature_grs.hit[[attrib]], levels = manual_levels))
  }

  feature_grs.hit[, ymin := as.numeric(get(attrib))+pad-1]
  feature_grs.hit[, ymax := as.numeric(get(attrib))-pad]

  lev = levels(feature_grs.hit[[attrib]])

  p_gr = ggplot(feature_grs.hit,
                aes_string(
                  xmin = "start",
                  xmax = "end",
                  ymin = "ymin",
                  ymax = "ymax",
                  fill = fill_VAR,
                  color = color_VAR))

  p_gr = p_gr + geom_rect()

  p_gr = p_gr +
    scale_y_continuous(breaks = seq_along(lev)-.5, labels = lev, limits = c(0, length(lev))) +
    theme_classic() +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          legend.position = legend.position)

  p_gr = .apply_x_scale(p_gr, x_scale = x_scale, as.character(seqnames(query_gr)))
  p_gr = .apply_x_lim(p_gr, query_gr, flip_x = flip_x)
  p_gr
}

