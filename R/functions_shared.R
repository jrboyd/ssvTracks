
.apply_x_scale = function(p, x_scale = c("bp", "kbp", "Mbp")[2], prefix = NULL){
  stopifnot(x_scale %in% c("bp", "kbp", "Mbp"))
  x_label_FUN = switch (
    x_scale,
    bp = {
      function(x)x
    },
    kbp = {
      function(x)x/1e3
    },
    Mbp = {
      function(x)x/1e6
    }
  )

  p = p +
    labs(x = ifelse(is.null(prefix), x_scale, paste(prefix, x_scale))) +
    scale_x_continuous(labels = x_label_FUN)
  p
}

.apply_x_lim = function(p, query_gr, flip_x = NULL){
  rng = c(start(query_gr), end(query_gr))
  if(is.null(flip_x)){
    flip_x = as.character(strand(query_gr) == "-")
  }
  if(flip_x){
    rng = rev(rng)
  }
  p = p + coord_cartesian(xlim = rng, expand = TRUE)
  p
}

.check_query_gr = function(query_gr){
  if(length(query_gr) > 1){
    stop("query_gr must be a single region")
  }
}

.check_attribute = function(ATTRIB_VAR, DEFAULT_VALUE){
  if(is.null(ATTRIB_VAR)){
    ATTRIB_VAR = DEFAULT_VALUE
  }
  ATTRIB_VAR
}


.check_dt_for_attribute = function(
    target_dt,
    ATTRIB_VAR,
    DEFAULT_VALUE){
  if(is.null(target_dt[[ATTRIB_VAR]])){
    if(ATTRIB_VAR == DEFAULT_VALUE){
      target_dt[[ATTRIB_VAR]] = rep("", nrow(target_dt))
    }else{
      target_dt[[ATTRIB_VAR]] = target_dt$sample
    }
  }
  target_dt
}

.check_show_aes = function(ATTRIB_VAR, DEFAULT_VALUE){
  #determine if color or fill get shown
  if(ATTRIB_VAR == DEFAULT_VALUE){
    show_color = FALSE
  }else{
    show_color = TRUE
  }
  show_color
}

.get_group_vars = function(color_VAR, fill_VAR, facet_VAR){
  cn = unique(c(color_VAR, fill_VAR, facet_VAR))
  cn = setdiff(cn, c(DEF_COLOR_, DEF_FILL_))
  cn
}
