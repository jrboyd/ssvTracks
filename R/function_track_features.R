
#' track_features
#'
#' @param feature_grs
#' @param query_gr
#' @param attrib
#' @param pad
#' @param flip_x
#' @param manual_levels
#'
#' @return
#' @export
#'
#' @examples
track_features = function(feature_grs, query_gr, attrib = "cluster_id", pad = .1, flip_x = NULL, manual_levels = NULL){
  rng = c(start(query_gr), end(query_gr))
  if(is.null(flip_x)){
    flip_x = as.character(strand(query_gr) == "-")
  }

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
  c_id = mcols(feature_grs.hit)[[attrib]]
  mcols(feature_grs.hit) = NULL
  names(feature_grs.hit) = NULL
  mcols(feature_grs.hit)[[attrib]] = c_id
  feature_grs.hit = unique(feature_grs.hit)
  feature_grs.hit = as.data.table(feature_grs.hit)



  if(!is.null(manual_levels)){
    stopifnot(feature_grs.hit[[attrib]] %in% manual_levels)
    # feature_grs.hit[[attrib]] = factor(feature_grs.hit[[attrib]], levels = manual_levels)
    set(feature_grs.hit, j = attrib, value = factor(feature_grs.hit[[attrib]], levels = manual_levels))
  }

  feature_grs.hit[, ymin := as.numeric(get(attrib))+pad-1]
  feature_grs.hit[, ymax := as.numeric(get(attrib))-pad]

  lev = levels(feature_grs.hit[[attrib]])

  if(flip_x){
    rng = rev(rng)

  }

  p_gr = ggplot(feature_grs.hit, aes_string(xmin = "start", xmax = "end", ymin = "ymin", ymax = "ymax", fill = attrib)) +
    geom_rect(color = NA) +
    scale_y_continuous(breaks = seq_along(lev)-.5, labels = lev, limits = c(0, length(lev))) +
    scale_x_continuous(labels = function(x)x/10^3) +
    coord_cartesian(xlim = rng, expand = TRUE) +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          legend.position = "bottom")
  p_gr
}
