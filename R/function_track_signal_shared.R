DEF_COLOR_ = "default_color__"
DEF_FILL_ = "default_fill__"



#' the idea is that all signal type tracks have a basic framework
#' 1) do stuff before fetch
#' 2) fetch
#' 3) do stuff after fetch
#'
#' .track_all_common_before_fetch
#'   .track_rna_common_before_fetch
#'   .track_chip_common_before_fetch
#' .track_all_common_after_fetch
#'   .track_rna_common_after_fetch
#'   .track_chip_common_after_fetch
#'
#' Used in:
#' track_rna.PE
#' track_rna.SE
#' track_chip.SE
#' track_chip.PE
.track_all_common_before_fetch = function(
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
  if(is.character(summary_FUN)){
    summary_FUN = switch (summary_FUN,
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

  #### check color and fill  ####
  color_VAR = .check_attribute(
    ATTRIB_VAR = color_VAR,
    DEFAULT_VALUE = DEF_COLOR_
  )
  signal_files = .check_dt_for_attribute(
    target_dt = signal_files,
    ATTRIB_VAR = color_VAR,
    DEFAULT_VALUE = DEF_COLOR_
  )
  fill_VAR = .check_attribute(
    ATTRIB_VAR = fill_VAR,
    DEFAULT_VALUE = DEF_FILL_
  )
  signal_files = .check_dt_for_attribute(
    target_dt = signal_files,
    ATTRIB_VAR = fill_VAR,
    DEFAULT_VALUE = DEF_FILL_
  )

  if(is.null(target_strand)){
    target_strand = as.character(strand(query_gr))
  }

  show_color = .check_show_aes(
    ATTRIB_VAR = color_VAR,
    DEFAULT_VALUE = DEF_COLOR_
  )
  if(show_color){
    if(is.null(color_mapping)){
      if(color_VAR == "strand"){
        color_mapping = seqsetvis::safeBrew(c("-", "+"))
      }else{
        color_mapping = seqsetvis::safeBrew(signal_files[[color_VAR]])
      }
      if(!is.na(color_mapping["input"])){
        color_mapping["input"] = "gray"
      }
    }
  }

  show_fill = .check_show_aes(
    ATTRIB_VAR = fill_VAR,
    DEFAULT_VALUE = DEF_FILL_
  )
  if(show_fill){
    if(is.null(fill_mapping)){
      if(fill_VAR == "strand"){
        fill_mapping = seqsetvis::safeBrew(c("-", "+"))
      }else{
        fill_mapping = seqsetvis::safeBrew(signal_files[[fill_VAR]])
      }
      if(!is.na(fill_mapping["input"])){
        fill_mapping["input"] = "gray"
      }
    }
  }

  #update args and return
  args$flip_x = flip_x
  args$target_strand = target_strand
  args$signal_files = signal_files
  args$nwin = nwin
  args$summary_FUN = summary_FUN
  args$color_VAR = color_VAR
  args$fill_VAR = fill_VAR
  args$color_mapping = color_mapping
  args$fill_mapping = fill_mapping
  args
}

.track_all_common_after_fetch = function(
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
    fill_mapping = NULL,
    color_over_fill = TRUE,
    facet_VAR = "sample",

    legend.position = "right",
    names_on_right = TRUE,
    target_strand = NULL,
    flip_strand = FALSE,
    return_data = FALSE,
    show_splice = FALSE,
    show_pileup = TRUE,
    ...){
  if(!is.null(bw_dt.raw$y_cap_value)){
    bw_dt.raw[, y_raw := y]
    bw_dt.raw[, y := y_raw / y_cap_value]
    if(y_label == "signal") y_label = "normalized signal"
    if(show_splice){
      splice_dt.raw[, y_raw := y]
      splice_dt.raw[, y := y_raw / y_cap_value]
    }
    if(!is.null(bw_dt.raw$mapped_reads)){
      warning("Using y_cap_value normalization and ignoring mapped_reads")
    }
  }else if(!is.null(bw_dt.raw$mapped_reads)){
    bw_dt.raw[, y_raw := y]
    bw_dt.raw[, y := y_raw / mapped_reads * 1e6]
    if(y_label == "signal") y_label = "RPM"
    if(show_splice){
      splice_dt.raw[, y_raw := y]
      splice_dt.raw[, y := y_raw / mapped_reads * 1e6]
    }
  }

  ####  grouping ####
  group_vars = .get_group_vars(
    color_VAR = color_VAR,
    fill_VAR = fill_VAR,
    facet_VAR = facet_VAR
  )

  bw_dt = bw_dt.raw[, list(y = mean(y)), c(unique(c(color_VAR, fill_VAR, facet_VAR, "x", "start", "end")))]
  # bw_dt$sample = apply(bw_dt[, group_vars, with = FALSE], 1, paste, collapse = " ")
  bw_dt = bw_dt[order(get(color_VAR))][order(get(fill_VAR))][order(get(facet_VAR))]
  # bw_dt$sample = factor(bw_dt$sample, levels = unique(bw_dt$sample))
  bw_dt[[facet_VAR]] = factor(bw_dt[[facet_VAR]], levels = unique(bw_dt[[facet_VAR]]))

  bw_dt[, x := (end + start)/2]
  if(nspline > 1){
    bw_dt = seqsetvis::applySpline(bw_dt, n = nspline, by_ = unique(c(color_VAR, fill_VAR, facet_VAR)))
    #bw_dt[, x := (end + start)/2]
    # set(bw_dt, j = "x", value = start(query_gr) + width(query_gr)*bw_dt$x)
    #does x need to be end(query_gr) - width(query_gr) if strand is negative?
  }

  bw_dt[y > ceiling_value, y := ceiling_value]
  bw_dt[y < floor_value, y := floor_value]

  if(return_data){
    return(list(reads = bw_dt))
  }

  p_rna = ggplot(bw_dt)

  #### show color and fill  ####
  show_color = .check_show_aes(
    ATTRIB_VAR = color_VAR,
    DEFAULT_VALUE = DEF_COLOR_
  )
  show_fill = .check_show_aes(
    ATTRIB_VAR = fill_VAR,
    DEFAULT_VALUE = DEF_FILL_
  )
  todo = character()
  if(show_color){
    todo = c(todo, "show_color")
  }
  if(show_fill){
    todo = c(todo, "show_fill")
  }
  if(color_over_fill){
    todo = rev(todo)
  }
  for(td in todo){
    if(td == "show_color"){
      if(!all(bw_dt[[color_VAR]] %in% names(color_mapping))){
        if(length(color_mapping) == 1){
          color_mapping = rep(color_mapping, length(unique(bw_dt[[color_VAR]])))
          if(is.factor(bw_dt[[color_VAR]])){
            names(color_mapping) = levels(bw_dt[[color_VAR]])
          }else{
            names(color_mapping) = unique(bw_dt[[color_VAR]])
          }
        }else{
          stop("Missing values from color_mapping for color_VAR \"", color_VAR, "\":\n",
               paste(setdiff(unique(bw_dt[[color_VAR]]), names(color_mapping)), collapse = ", "))
        }
      }

      path_show.legend = length(unique(color_mapping)) > 1
      if(show_pileup){
        p_rna = p_rna +
          geom_path(aes_string(x = "x", y = "y", color = color_VAR), alpha = color_alpha, show.legend = path_show.legend) +
          scale_color_manual(values = color_mapping)
      }
    }
    if(td == "show_fill"){
      if(!all(bw_dt[[fill_VAR]] %in% names(fill_mapping))){
        if(length(fill_mapping) == 1){
          fill_mapping = rep(fill_mapping, length(unique(bw_dt[[fill_VAR]])))
          if(is.factor(bw_dt[[fill_VAR]])){
            names(fill_mapping) = levels(bw_dt[[fill_VAR]])
          }else{
            names(fill_mapping) = unique(bw_dt[[fill_VAR]])
          }
        }else{
          stop("Missing values from fill_mapping for fill_VAR \"", fill_VAR, "\":\n",
               paste(setdiff(unique(bw_dt[[fill_VAR]]), names(fill_mapping)), collapse = ", "))
        }
      }

      ribbon_show.legend = length(unique(fill_mapping)) > 1
      if(show_pileup){
        p_rna = p_rna +
          geom_ribbon(aes_string(x = "x", ymin = 0, ymax = "y", fill = fill_VAR), color = fill_outline_color, alpha = fill_alpha, show.legend = ribbon_show.legend) +
          scale_fill_manual(values = fill_mapping)
      }
    }
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


