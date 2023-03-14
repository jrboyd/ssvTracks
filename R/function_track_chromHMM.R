load_state_dt = function(return_GRanges = FALSE){
  chrom_dir = "/slipstream/home/conggao/ChromHMM_v1.23/output18model_20220919"
  state_files = dir(chrom_dir, pattern = "dense.bed", full.names = TRUE)
  state_files = state_files[!grepl("sorted_dense.bed", state_files)]
  state_dt.l = lapply(state_files, function(f){
    nam = sub("_seg.+", "", basename(f))
    state_dt = suppressWarnings({fread(f)})
    colnames(state_dt) = c("seqnames", "start", "end", "state", "zero", "dot", 'start2', "end2", "color")
    state_dt$sample = nam
    state_dt
  })
  if(return_GRanges){
    names(state_dt.l) = sub("_seg.+", "", basename(state_files))
    lapply(state_dt.l, GRanges)
  }else{
    rbindlist(state_dt.l)
  }
}

track_chromHMM_states = function(){

  state_grs = load_state_dt(return_GRanges = TRUE)

  chr_state_dt = ssvFetchGRanges(state_grs, query_gr, attrib_var = "state", return_data.table = TRUE, win_size = 100)
  chr_state_dt$score = 1
  chr_state_l = split(GRanges(chr_state_dt), chr_state_dt$y)
  gr = chr_state_l[[1]]
  res.by_state = lapply(chr_state_l, function(gr){
    gr.by_sample = split(gr, sub("_18.+", "", gr$sample))
    res.by_sample = lapply(gr.by_sample, function(gr_i){
      viewGRangesWinSummary_dt(gr_i, query_gr, n_tiles = 200)
    })
    res = rbindlist(res.by_sample, idcol = "sample")
    res
  })
  chr_state_dt = rbindlist(res.by_state, idcol = "state")

  chr_state_dt[, x := (start + end)/2]
  chr_state_dt$state = factor(chr_state_dt$state, levels = sort(as.numeric(unique(chr_state_dt$state))))
  chr_state_dt[, state_name := state_annotation[state]]
  chr_state_dt$state_name = factor(chr_state_dt$state_name, levels = state_annotation)

  cell_colors = safeBrew(chr_state_dt$sample, pal = "Set1")

  chr_state_dt$x %>% table

  p_chromHMM = ggplot(chr_state_dt, aes(x = x, y = y, color = sample)) +
    geom_path() +
    facet_grid(state_name~., scales = "free_y") +
    scale_x_continuous(labels = function(x)x/1e3) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
    scale_color_manual(values = cell_colors) +
    labs(y = "fraction", x = "kbp") +
    theme(panel.spacing.y = unit(.7, "lines")) +
    theme(strip.text.y = element_text(angle= 0, hjust = 0), legend.position = "bottom") +
    guides(color = guide_legend(ncol = 2) )

}
