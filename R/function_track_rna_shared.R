#' do most error catching and set dynamic args
#' returns modified copy of args
.track_rna_common_before_fetch = function(
    signal_files,
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
    target_strand = NULL,
    flip_strand = FALSE,
    return_data = FALSE,
    show_splice = TRUE,
    min_splice_count = 0,
    splice_within_range_only = NULL,
    ...){
  env = as.list(sys.frame(sys.nframe()))
  args = c(as.list(env), list(...))
  args2 = do.call(.track_all_common_before_fetch, args = args)
  args2
}

.track_rna_common_after_fetch = function(
    bw_dt.raw,
    splice_dt.raw,
    signal_files,
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
    facet_VAR = "sample",
    fill_mapping = NULL,
    legend.position = "right",
    names_on_right = TRUE,
    target_strand = NULL,
    flip_strand = FALSE,
    return_data = FALSE,
    show_splice = TRUE,
    min_splice_count = 0,
    splice_within_range_only = NULL,
    ...){
  env = as.list(sys.frame(sys.nframe()))
  args = c(as.list(env), list(...))
  pre_out = do.call(.track_all_common_after_fetch, args = args)

  if(show_splice){
    #### color and fill  ####
    #default color and fill attribute variable
    #determine if color or fill get shown
    #enforce proper color and fill variables in bw_dt.raw
    splice_dt.raw = .check_dt_for_attribute(
      target_dt = splice_dt.raw,
      ATTRIB_VAR = color_VAR,
      DEFAULT_VALUE = DEF_COLOR_
    )
    splice_dt.raw = .check_dt_for_attribute(
      target_dt = splice_dt.raw,
      ATTRIB_VAR = fill_VAR,
      DEFAULT_VALUE = DEF_FILL_
    )
  }

  group_vars = .get_group_vars(
    color_VAR = color_VAR,
    fill_VAR = fill_VAR,
    facet_VAR = facet_VAR
  )

  if(show_splice & !is.null(splice_dt.raw)){
    splice_dt = splice_dt.raw[, list(y = mean(y)), c(unique(c(color_VAR, fill_VAR, facet_VAR, "start", "end", "seqnames")))]
    splice_dt$sample = apply(splice_dt[, group_vars, with = FALSE], 1, paste, collapse = " ")
    splice_dt = splice_dt[order(get(color_VAR))][order(get(fill_VAR))][order(get(facet_VAR))]
    splice_dt$sample = factor(splice_dt$sample, levels = unique(splice_dt$sample))
  }else{
    splice_dt = NULL
  }

  if(return_data){
    return(c(pre_out, list(splice = splice_dt)))
  }else{
    p_rna = pre_out
  }

  if(show_splice){
    #ensure that is splice event in 1 sample it's plotted in all
    #avoids plotting 100 but not plotting 99
    valid_start_end = unique(splice_dt[y >= min_splice_count][, .(start, end)])
    valid_splice_dt = merge(splice_dt, valid_start_end, by = c('start', 'end'))
    if(nrow(valid_splice_dt) > 0){
      if(!is.null(splice_within_range_only)){
        gr = GRanges(valid_splice_dt)
        strand(gr) = "+"
        gr.starts = promoters(gr, 1, 1)
        gr.ends = promoters(invertStrand(gr), 1, 1)
        if(is(splice_within_range_only, "GRanges")){
          olaps1 = findOverlaps(gr.starts, splice_within_range_only, ignore.strand = TRUE)
          olaps2 = findOverlaps(gr.ends, splice_within_range_only, ignore.strand = TRUE)
          hit_i = intersect(queryHits(olaps1), queryHits(olaps2))
        }else if(is.list(splice_within_range_only) || is(splice_within_range_only, "GRangesList")){
          stopifnot(length(splice_within_range_only) == 2)
          stopifnot(is(splice_within_range_only[[1]], "GRanges"))
          stopifnot(is(splice_within_range_only[[2]], "GRanges"))
          olaps1a = findOverlaps(gr.starts, splice_within_range_only[[1]], ignore.strand = TRUE)
          olaps2a = findOverlaps(gr.ends, splice_within_range_only[[2]], ignore.strand = TRUE)


          olaps1b = findOverlaps(gr.starts, splice_within_range_only[[2]], ignore.strand = TRUE)
          olaps2b = findOverlaps(gr.ends, splice_within_range_only[[1]], ignore.strand = TRUE)
          hit_i = union(
            intersect(queryHits(olaps1a), queryHits(olaps2a)),
            intersect(queryHits(olaps1b), queryHits(olaps2b))
          )
        }else{
          stop("splice_within_range_only must be either a GRanges or a list of length 2 containing GRanges")
        }
        valid_splice_dt = valid_splice_dt[hit_i,]
      }
    }
    if(color_VAR == DEF_COLOR_){
      p_rna = p_rna + ggbio::geom_arch(
        data = valid_splice_dt, aes(x = start, xend = end, height = y), show.legend = FALSE,
        color = "black"
      )
    }else{
      color_VAR = ensym(color_VAR)
      p_rna = p_rna + ggbio::geom_arch(
        data = valid_splice_dt, aes(x = start, xend = end, height = y, color = !!color_VAR), show.legend = FALSE#,
        # color = "black"
      )
    }

  }
  p_rna
}
