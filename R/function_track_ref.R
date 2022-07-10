#' track_ref
#'
#' @param ref
#' @param query_gr
#' @param flip_x
#' @param exon_height
#' @param intron_thickness
#' @param show_tss
#' @param tss_size
#' @param tss_color
#' @param return_data
#' @param tss_up
#' @param tss_over
#' @param tss_arrow_size
#'
#' @return
#' @export
#'
#' @examples
track_ref = function(ref = "~/gencode.v28.annotation.gtf.gz", query_gr, flip_x = FALSE,
                     exon_height = .4, intron_thickness = 2,
                     show_tss = FALSE,
                     tss_size = 1.5,
                     tss_color = "red",
                     return_data = FALSE,
                     tss_up = .4,
                     tss_over = .05,
                     tss_arrow_size = .4,
                     plus_strand_color = "black",
                     minus_strand_color = "darkgray",
                     legend.position = "right"){
  if(!is(ref, "GRanges")){
    if(file.exists(ref)){
      ref = rtracklayer::import.gff(ref, format =  "gtf", feature.type = "exon")
    }else{
      stop("ref must be gtf loaded as GRanges or path to gtf")
    }
  }else{
      ref = subset(ref, type == "exon")
  }
  rng = c(IRanges::start(query_gr), IRanges::end(query_gr))
  hit_genes = unique(IRanges::subsetByOverlaps(ref, query_gr, ignore.strand = TRUE)$gene_id)
  ref_dt = as.data.table(subset(ref, gene_id %in% hit_genes))
  yvar = "gene_name"

  ref_dt[[yvar]] = factor(ref_dt[[yvar]])
  ref_dt$ymin = as.numeric(ref_dt[[yvar]]) - .5 * exon_height
  ref_dt$ymax = as.numeric(ref_dt[[yvar]]) + .5 * exon_height

  ref_dt_base = ref_dt[, .(start = min(start), end = max(end), y = mean(c(ymin, ymax)), strand = unique(strand)), by = yvar]
  ref_dt_base = melt(ref_dt_base, id.vars = c("strand", "gene_name", "y"), value.name = "x")

  p_ref = ggplot() +
    geom_line(data = ref_dt_base, aes_string(x = "x", y = "y", color = "strand", group = "gene_name"), size = intron_thickness) +
    geom_rect(data = ref_dt, aes(fill = strand, color = strand, xmin = start, xmax = end, ymin = ymin, ymax = ymax)) +
    scale_y_continuous(breaks = seq_along(levels(ref_dt[[yvar]])),
                       labels = function(x)levels(ref_dt[[yvar]])[round(x)]) +
    theme_classic() +
    theme(legend.position = legend.position) +
    labs(y = "gene\nannotation", x = "kbp")

  if(flip_x){
    p_ref = p_ref +
      scale_x_reverse(labels = function(x)x/10^3) +
      coord_cartesian(xlim = rev(rng), expand = TRUE) +
      scale_fill_manual(values = c("-" = plus_strand_color, "+" = minus_strand_color)) +
      scale_color_manual(values = c("-" = plus_strand_color, "+" = minus_strand_color))
  }else{
    p_ref = p_ref +
      scale_x_continuous(labels = function(x)x/10^3) +
      coord_cartesian(xlim = rng, expand = TRUE ) +
      scale_fill_manual(values = c("+" = plus_strand_color, "-" = minus_strand_color)) +
      scale_color_manual(values = c("+" = plus_strand_color, "-" = minus_strand_color))
  }
  tss_dt = ref_dt[,
                  .(
                    tss = ifelse(strand == "+", min(start), max(end)),
                    y = mean(c(ymin, ymax)),
                    strand = unique(strand),
                    gene_name = unique(gene_name)
                  ),
                  by = "transcript_id"]
  if(show_tss){
    if(FALSE){
      p_ref = p_ref +
        geom_point(data = tss_dt, aes(x = tss, y = y), color = tss_color, size = tss_size) +
        labs(caption = "tss") +
        theme(plot.caption = element_text(size = 14, color = "red"), legend.position = "bottom")
    }else{
      tss_dt[, y_tss_up := y + tss_up]
      tss_dt[, x_tss_over := tss + (ifelse(strand == "-", -1, 1) * tss_over) * diff(rng)]
      p_ref = p_ref +
        geom_segment(data = tss_dt,
                     aes(x = tss,
                         y = y,
                         xend = tss,
                         yend = y_tss_up,
                         color = strand),
                     size = tss_size,
                     lineend = "round",
                     show.legend = FALSE) +
        geom_segment(data = tss_dt,
                     aes(x = tss,
                         y = y_tss_up,
                         xend = x_tss_over,
                         yend = y_tss_up,
                         color = strand),
                     size = tss_size,
                     lineend = "round",
                     arrow = arrow(length = unit(tss_arrow_size, "cm")),
                     show.legend = FALSE)
    }
  }
  if(return_data){return(list(ref_dt = ref_dt, tss_dt = tss_dt, label_dt = ref_dt_base))}
  p_ref
}
