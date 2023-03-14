track_hic.fetch = function(hic_files.u2os){
  hic_dt = dthic::fetch_hic_query_gr_data(sel_dt = data.table(file = hic_files.u2os),
                                   query_gr = query_gr.chr6,
                                   sel_binsize = 5e4,
                                   return_data = TRUE)
  hmat.l = lapply(hic_files.u2os, function(f){
    dthic::HiC_matrix.from_hic(hic_f = f,
                               query_gr = query_gr.chr6,
                               norm = "NONE",
                               matrix = "oe",
                               bin_size = 5e4)
  })
  shared_1d = unique(rbind(
    hmat.l[[1]]@hic_1d,
    hmat.l[[2]]@hic_1d
  ))
  shared_1d[, seqnames := paste0("chr", seqnames)]
  hmat.l[[1]]@hic_1d = shared_1d
  hmat.l[[2]]@hic_1d = shared_1d

  hmat = hmat.l[[1]]
  for(i in seq(2, length(hmat.l))){
    hmat = hmat + hmat.l[[i]]
  }
  hmat = HiC_matrix_wInsulation(hmat)
  hmat
}
