.track_all_common_before_fetch = function(
    signal_files,
    query_gr,
    fetch_fun = seqsetvis::ssvFetchBam,
    win_FUN = c("mean", "max")[2],
    sum_FUN = NULL,
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
    ...
){
  env = as.list(sys.frame(sys.nframe()))
  args = c(as.list(env), list(...))
  .check_query_gr(query_gr)
  if(is.null(color_VAR) & is.null(fill_VAR)){
    stop("At least one of color_VAR or fill_VAR must be set.")
  }
  if(is.null(flip_x)){
    flip_x = as.character(strand(query_gr) == "-")
  }
  if(is.null(sum_FUN)){
    sum_FUN = switch (win_FUN,
                      max = function(x, w)max(x),
                      mean = weighted.mean
    )
  }

  if(nwin > width(query_gr)){
    nwin = width(query_gr)
  }

  if(!is.data.frame(signal_files)){
    if(!is.null(color_VAR)){
      if(!color_VAR == "sample"){
        stop("With file paths as signal_files, color_VAR must be \"sample\"")
      }
    }
    if(!is.null(fill_VAR)){
      if(!fill_VAR == "sample"){
        stop("With file paths as signal_files, fill_VAR must be \"sample\"")
      }
    }
    if(!is.null(facet_VAR)){
      if(!facet_VAR == "sample"){
        stop("With file paths as signal_files, facet_VAR must be \"sample\"")
      }
    }
    if(is.null(names(signal_files))){
      signal_files = data.frame(
        file = signal_files,
        sample = basename(signal_files)
      )
    }else{
      signal_files = data.frame(
        file = signal_files,
        sample = names(signal_files)
      )
    }

  }
  if(is.null(target_strand)){
    target_strand = as.character(strand(query_gr))
  }

  #update args and return
  args$flip_x = flip_x
  args$target_strand = target_strand
  args$signal_files = signal_files
  args$nwin = nwin
  args$sum_FUN = sum_FUN
  args
}

.track_all_common_before_fetch = function(
    signal_files,
    query_gr,
    fetch_fun = seqsetvis::ssvFetchBam,
    win_FUN = c("mean", "max")[2],
    sum_FUN = NULL,
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
    ...){

}

#' do most error catching and set dynamic args
#' returns modified copy of args
.track_rna_common_before_fetch = function(
    signal_files,
    query_gr,
    fetch_fun = seqsetvis::ssvFetchBam,
    win_FUN = c("mean", "max")[2],
    sum_FUN = NULL,
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
    min_splice_count = 10,
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
    win_FUN = c("mean", "max")[2],
    sum_FUN = NULL,
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
    show_splice = TRUE,
    min_splice_count = 10,
    target_strand = NULL,
    flip_strand = FALSE,
    return_data = FALSE,
    ...){

  if(!is.null(bw_dt.raw$mapped_reads)){
    bw_dt.raw[, y_raw := y]
    bw_dt.raw[, y := y_raw / mapped_reads * 1e6]
    if(y_label == "signal") y_label = "RPM"
    if(show_splice){
      splice_dt.raw[, y_raw := y]
      splice_dt.raw[, y := y_raw / mapped_reads * 1e6]
    }
  }

  #### TODO fix bloody mess ####
  DEF_COLOR_ = "default_color__"
  if(is.null(color_VAR)){
    color_VAR = DEF_COLOR_
    show_color = FALSE
  }else{
    show_color = TRUE
  }
  if(is.null(bw_dt.raw[[color_VAR]])){
    if(color_VAR == DEF_COLOR_){
      bw_dt.raw[[color_VAR]] = rep("", nrow(bw_dt.raw))
      if(show_splice){
        splice_dt.raw[[color_VAR]] = rep("", nrow(splice_dt.raw))
      }
    }else{
      bw_dt.raw[[color_VAR]] = bw_dt.raw$sample
      if(show_splice){
        splice_dt.raw[[color_VAR]] = splice_dt.raw$sample
      }
    }
  }

  DEF_FILL_ = "default_fill__"
  if(is.null(fill_VAR)){
    fill_VAR = DEF_FILL_
    show_fill = FALSE
  }else{
    show_fill = TRUE
  }
  if(is.null(bw_dt.raw[[fill_VAR]])){
    if(fill_VAR == DEF_FILL_){
      bw_dt.raw[[fill_VAR]] = rep("", nrow(bw_dt.raw))
      if(show_splice){
        splice_dt.raw[[fill_VAR]] = rep("", nrow(splice_dt.raw))
      }
    }else{
      bw_dt.raw[[fill_VAR]] = bw_dt.raw$sample
      if(show_splice){
        splice_dt.raw[[fill_VAR]] = splice_dt.raw$sample
      }
    }
  }
  ####  ####
  cn = unique(c(color_VAR, fill_VAR, facet_VAR))
  cn = setdiff(cn, c(DEF_COLOR_, DEF_FILL_))

  bw_dt = bw_dt.raw[, list(y = mean(y)), c(unique(c(color_VAR, fill_VAR, facet_VAR, "x", "start", "end")))]
  bw_dt$sample = apply(bw_dt[, cn, with = FALSE], 1, paste, collapse = " ")
  bw_dt = bw_dt[order(get(color_VAR))][order(get(fill_VAR))][order(get(facet_VAR))]
  bw_dt$sample = factor(bw_dt$sample, levels = unique(bw_dt$sample))

  if(show_splice & !is.null(splice_dt.raw)){
    splice_dt = splice_dt.raw[, list(y = mean(y)), c(unique(c(color_VAR, fill_VAR, facet_VAR, "start", "end")))]
    splice_dt$sample = apply(splice_dt[, cn, with = FALSE], 1, paste, collapse = " ")
    splice_dt = splice_dt[order(get(color_VAR))][order(get(fill_VAR))][order(get(facet_VAR))]
    splice_dt$sample = factor(splice_dt$sample, levels = unique(splice_dt$sample))
  }else{
    splice_dt = NULL
  }

  if(nspline > 1){
    bw_dt = seqsetvis::applySpline(bw_dt, n = nspline, by_ = c("sample"))
  }
  #### SOMETHING WRONG WITH FLIP_X and x calc of bw_dt
  # if(flip_x){
  #   bw_dt[, x := max(end) - (max(end) - min(start))*x]
  # }else{
  #   bw_dt[, x := min(start) + (max(end) - min(start))*x]
  # }
  bw_dt[, x := (end + start)/2]

  bw_dt[y > ceiling_value, y := ceiling_value]
  bw_dt[y < floor_value, y := floor_value]

  if(return_data){
    return(list(reads = bw_dt, splices = splice_dt))
  }

  p_rna = ggplot(bw_dt)
  if(show_color){
    if(is.null(color_mapping)){
      color_mapping = seqsetvis::safeBrew(bw_dt[[color_VAR]])
      if(!is.na(color_mapping["input"])){
        color_mapping["input"] = "gray"
      }
    }
    stopifnot(all(bw_dt[[color_VAR]] %in% names(color_mapping)))

    p_rna = p_rna +
      geom_path(aes_string(x = "x", y = "y", color = color_VAR), alpha = color_alpha) +
      scale_color_manual(values = color_mapping)
  }
  if(show_fill){
    if(is.null(fill_mapping)){
      fill_mapping = seqsetvis::safeBrew(bw_dt[[fill_VAR]])
      if(!is.na(fill_mapping["input"])){
        fill_mapping["input"] = "gray"
      }
    }
    stopifnot(all(bw_dt[[fill_VAR]] %in% names(fill_mapping)))

    p_rna = p_rna +
      geom_ribbon(aes_string(x = "x", ymin = 0, ymax = "y", fill = fill_VAR), color = fill_outline_color, alpha = fill_alpha) +
      scale_fill_manual(values = fill_mapping)
  }

  if(show_splice){
    p_rna =   p_rna + ggbio::geom_arch(
      data = splice_dt[y >= min_splice_count], aes(x = start, xend = end, height = y),
      color = "black"
    )
  }

  facet_switch = if(names_on_right){
    NULL
  }else{
    "y"
  }

  p_rna = p_rna +
    labs(y = paste0(y_label, " (", target_strand, ")")) +
    facet_grid(formula(paste0(facet_VAR, "~.")), switch = facet_switch) +
    theme_classic() +
    theme(strip.background = element_blank(), strip.text.y = element_text(angle = 0)) +
    theme(legend.position = legend.position)
  p_rna = .apply_x_scale(p_rna, x_scale, as.character(seqnames(query_gr)))
  p_rna = .apply_x_lim(p_rna, query_gr, flip_x)
  p_rna
}

#' track_rna.SE
#'
#' @param signal_files
#' @param query_gr
#' @param fetch_fun
#' @param win_FUN
#' @param sum_FUN
#' @param flip_x
#' @param nwin
#' @param nspline
#' @param fill_outline_color
#' @param y_label
#' @param x_scale
#' @param floor_value
#' @param ceiling_value
#' @param color_VAR
#' @param color_mapping
#' @param fill_VAR
#' @param fill_mapping
#' @param legend.position
#' @param names_on_right
#' @param show_splice
#' @param min_splice_count
#' @param target_strand
#' @param flip_strand logical. Should aligned strand be flipped relative to query_gr.  The vast majority of SE libraries are actually flip_strand, so default flip_strand is TRUE.
#' @param ...
#'
#' @return
#' @export
#' @import seqsetvis
#'
#' @examples
track_rna.SE = function(
    signal_files,
    query_gr,
    fetch_fun = seqsetvis::ssvFetchBam,
    win_FUN = c("mean", "max")[2],
    sum_FUN = NULL,
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
    show_splice = TRUE,
    min_splice_count = 10,
    target_strand = NULL,
    flip_strand = TRUE,
    return_data = FALSE,
    ...){
  env = as.list(sys.frame(sys.nframe()))
  args = c(as.list(env), list(...))
  args2 = do.call(.track_rna_common_before_fetch, args)
  for(var_name in names(args2)){
    assign(var_name, args2[[var_name]])
  }
  bw_dt.raw = fetch_fun(signal_files, query_gr,
                        win_method = "summary",  win_size = nwin,
                        summary_FUN = sum_FUN,
                        return_data.table = TRUE,
                        anchor = "left",
                        target_strand = target_strand,
                        flip_strand = flip_strand,
                        fragLens = NA)

  if(show_splice){
    splice_dt.raw = fetch_fun(signal_files, query_gr,
                              win_method = "summary",  win_size = nwin,
                              summary_FUN = sum_FUN,
                              return_data.table = TRUE,
                              anchor = "left",
                              fragLens = NA,
                              splice_strategy = "splice_count")
    if(flip_strand){
      splice_dt.raw[strand == "-", strand := "tmp"]
      splice_dt.raw[strand == "+", strand := "-"]
      splice_dt.raw[strand == "tmp", strand := "+"]
    }
    if(target_strand %in% c("+", "-")){
      splice_dt.raw = splice_dt.raw[strand == target_strand]
    }
    setnames(splice_dt.raw, "N", "y")
  }

  args2$bw_dt.raw = bw_dt.raw
  args2$splice_dt.raw = splice_dt.raw
  do.call(.track_rna_common_after_fetch, args2)
}

#' track_rna.PE
#'
#' @param signal_files
#' @param query_gr
#' @param fetch_fun
#' @param win_FUN
#' @param sum_FUN
#' @param flip_x
#' @param nwin
#' @param nspline
#' @param fill_outline_color
#' @param y_label
#' @param x_scale
#' @param floor_value
#' @param ceiling_value
#' @param color_VAR
#' @param color_mapping
#' @param fill_VAR
#' @param fill_mapping
#' @param legend.position
#' @param names_on_right
#' @param show_splice
#' @param min_splice_count
#' @param target_strand
#' @param flip_strand
#' @param ...
#'
#' @return
#' @export
#' @import seqsetvis
#'
#' @examples
track_rna.PE = function(
    signal_files,
    query_gr,
    win_FUN = c("mean", "max")[2],
    # sum_FUN = NULL,
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
    show_splice = TRUE,
    min_splice_count = 10,
    target_strand = NULL,
    flip_strand = FALSE,
    return_data = FALSE,
    ...){
  env = as.list(sys.frame(sys.nframe()))
  args = c(as.list(env), list(...))
  args2 = do.call(.track_rna_common_before_fetch, args)
  for(var_name in names(args2)){
    assign(var_name, args2[[var_name]])
  }

  bw_dt.raw = ssvFetchBamPE.RNA(
    signal_files, query_gr,
    win_method = "summary",
    win_size = nwin,
    return_data.table = TRUE,
    target_strand = target_strand,
    flip_strand = flip_strand
  )
  args2$bw_dt.raw = bw_dt.raw

  if(show_splice){
    splice_dt.raw = ssvFetchBamPE.RNA_splice(
      signal_files,
      query_gr,
      return_data.table = TRUE
    )
    if(target_strand %in% c("+", "-")){
      splice_dt.raw = splice_dt.raw[strand == target_strand]
    }
    setnames(splice_dt.raw, "N", "y")
    args2$splice_dt.raw = splice_dt.raw
  }else{
    args2["splice_dt.raw"] = list(NULL)
  }
  #
  #   browser()


  do.call(.track_rna_common_after_fetch, args2)
}

