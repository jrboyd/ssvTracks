#' track_gene_reference
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
track_gene_reference = function(ref = "~/../joeboyd/gencode.v36.annotation.gtf",
                                query_gr,
                                flip_x = FALSE,
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
                                legend.position = "right",
                                x_scale = c("bp", "kbp", "Mbp")[2]){
  .check_query_gr(query_gr)
  stopifnot(x_scale %in% c("bp", "kbp", "Mbp"))
  if(!is(ref, "GRanges")){
    if(file.exists(ref)){
      ref = rtracklayer::import.gff(ref, format =  "gtf", feature.type = "exon")
    }else{
      stop("ref must be gtf loaded as GRanges or path to gtf")
    }
  }else{
    ref = subset(ref, type == "exon")
  }
  hit_genes = unique(IRanges::subsetByOverlaps(ref, query_gr, ignore.strand = TRUE)$gene_id)
  if(length(hit_genes) < 1){
    p_ref = ggplot() +
      scale_y_continuous(breaks = 1,
                         labels = "no genes\nin view")
    ref_dt = data.table()
    tss_dt = data.table()
  }else{
    ref_dt = as.data.table(subset(ref, gene_id %in% hit_genes))
    yvar = "gene_name"

    ref_dt[[yvar]] = factor(ref_dt[[yvar]])
    ref_dt$ymin = as.numeric(ref_dt[[yvar]]) - .5 * exon_height
    ref_dt$ymax = as.numeric(ref_dt[[yvar]]) + .5 * exon_height

    ref_dt_base = ref_dt[, list(start = min(start), end = max(end), y = mean(c(ymin, ymax)), strand = unique(strand)), by = yvar]
    ref_dt_base = melt(ref_dt_base, id.vars = c("strand", "gene_name", "y"), value.name = "x")

    p_ref = ggplot() +
      geom_line(data = ref_dt_base, aes_string(x = "x", y = "y", color = "strand", group = "gene_name"), size = intron_thickness) +
      geom_rect(data = ref_dt, aes(fill = strand, color = strand, xmin = start, xmax = end, ymin = ymin, ymax = ymax)) +
      scale_y_continuous(breaks = seq_along(levels(ref_dt[[yvar]])),
                         labels = function(x)levels(ref_dt[[yvar]])[round(x)]) +
      theme_classic() +
      theme(legend.position = legend.position) +
      labs(y = "gene\nannotation", x = "kbp")

    tss_dt = ref_dt[,
                    list(
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
        tss_dt[, x_tss_over := tss + (ifelse(strand == "-", -1, 1) * tss_over) * abs(end(query_gr) - start(query_gr))]
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
  }
  if(flip_x){
    p_ref = p_ref +
      scale_fill_manual(values = c("-" = plus_strand_color, "+" = minus_strand_color)) +
      scale_color_manual(values = c("-" = plus_strand_color, "+" = minus_strand_color))
  }else{
    p_ref = p_ref +
      scale_fill_manual(values = c("+" = plus_strand_color, "-" = minus_strand_color)) +
      scale_color_manual(values = c("+" = plus_strand_color, "-" = minus_strand_color))
  }

  if(return_data){return(list(ref_dt = ref_dt, tss_dt = tss_dt, label_dt = ref_dt_base))}

  p_ref = .apply_x_scale(p_ref, x_scale, as.character(seqnames(query_gr)))
  p_ref = .apply_x_lim(p_ref, query_gr, flip_x)

  p_ref
}

#' track_gene_transcripts
#'
#' @param ref Either a GRanges object containing "exon" features coded in "type" variable and with "gene_name" entries or a path to gtf file including this information.
#' @param sel_gene_name The gene to plot transcripts for.  Must be present as "gene_name" entry in ref.
#' @param query_gr Determines plot x-range and impacts flip_x.  Default of NULL adjusts to 10\% wider than all transcripts plotted.
#' @param transcript_subset Character vector of transcript_id to plot. Also controls plot order of transcripts if provided. Default of NULL plots all transcripts.
#' @param flip_x If TRUE, x-axis is flipped to be decreasing.  Default of NULL will flip_x automatically when query_gr strand is negative.
#' @param x_scale One of "bp", "kbp", or "Mbp".  Scales x-axis labels accordingly.
#' @param intron_thickness Line thickness for introns.  0 disables intron drawing. Default is 0.
#' @param exon_height Relative height of exons. 1 fills entire row. Default is 0.8.
#' @param exon_color color for exon lines and fill. Default is "black"
#' @param intron_color color for intron lines. Default is the same as exon_color.
#' @param highlight_transcripts transcript_id items to highlight. Must be present. Default of NULL is no highlighting.
#' @param highlight_color line color for highlighted transcripts. Default is "red".
#' @param show_tss If TRUE, TSS is indicated by a leading in-line arrow. Default is FALSE.
#'
#' @return A ggplot in which exon structure of transcripts is shown.
#' @export
#'
#' @examples
#' gtf_file = system.file(package = "ssvTracks", "extdata/ESR1_exons.gtf", mustWork = TRUE)
#' ref_gr = rtracklayer::import.gff(gtf_file)
#'
#' track_gene_transcripts(ref_gr,
#'   sel_gene_name = "ESR1",
#'   intron_thickness = 2,
#'   exon_height = .8,
#'   show_tss = TRUE)
#'
#' track_gene_transcripts(ref_gr,
#'   sel_gene_name = "ESR1",
#'   intron_thickness = 2,
#'   highlight_transcripts = c("ENST00000338799.9"),
#'   highlight_color = "darkorange",
#'   show_tss = TRUE)
track_gene_transcripts = function(ref = "~/../joeboyd/gencode.v36.annotation.gtf",
                                  sel_gene_name,
                                  query_gr = NULL,
                                  transcript_subset = NULL,
                                  flip_x = NULL,
                                  x_scale = c("bp", "kbp", "Mbp")[2],
                                  intron_thickness = 0,
                                  exon_height = .8,
                                  exon_color = "black",
                                  intron_color = exon_color,
                                  highlight_transcripts = NULL,
                                  highlight_color = "red",
                                  show_tss = FALSE

){
  if(!is.null(query_gr)) .check_query_gr(query_gr)
  stopifnot(x_scale %in% c("bp", "kbp", "Mbp"))
  stopifnot(exon_height > 0 & exon_height <= 1)
  if(!is(ref, "GRanges")){
    if(file.exists(ref)){
      ref = rtracklayer::import.gff(ref, format =  "gtf", feature.type = "exon")
    }else{
      stop("ref must be gtf loaded as GRanges or path to gtf")
    }
  }else{
    ref = subset(ref, type == "exon")
  }

  ex_dt = data.table::as.data.table(subset(ref, gene_name == sel_gene_name))

  if(nrow(ex_dt[, .N, list(strand, seqnames)]) > 1){
    stop(sel_gene_name, " has mixed strand and/or seqnames in reference.  Specify transcript_subset such that only one strand and seqnames is present.")
  }

  if(is.null(query_gr)){
    query_gr = range(GenomicRanges::GRanges(ex_dt))
    query_gr = GenomicRanges::resize(query_gr, 1.1 * GenomicRanges::width(query_gr), fix = "center")
  }

  stopifnot(length(query_gr) == 1)
  if(is.null(flip_x)) flip_x = as.character(strand(query_gr)) == "-"

  if(!is.null(transcript_subset)){
    ex_dt = subset(ex_dt, transcript_id %in% transcript_subset)
    transcript_lev = rev(transcript_subset)
  }else{
    tx_sp = GenomicRanges::split(GenomicRanges::GRanges(ex_dt), ex_dt$transcript_id)
    transcript_lev = rev(names(sort(sapply(range(tx_sp), function(x)GenomicRanges::start(GenomicRanges::promoters(x, 1, 0))))))
    if(flip_x) transcript_lev = rev(transcript_lev)
  }

  ex_dt = ex_dt[, c("seqnames", "start", "end", "gene_id", "transcript_id", "exon_id", "gene_name", "strand")]
  ex_dt$transcript_id = factor(ex_dt$transcript_id, levels = transcript_lev)

  pad = (1-exon_height)/2
  set(ex_dt, j = "ymin", value = as.numeric(ex_dt$transcript_id)-1+pad)
  set(ex_dt, j = "ymax", value = as.numeric(ex_dt$transcript_id)-pad)

  ex_dt$color = exon_color
  if(!is.null(highlight_transcripts)){
    ex_dt[transcript_id %in% highlight_transcripts, color := highlight_color]
  }

  p_track_ref = ggplot()
  if(intron_thickness > 0){
    in_dt = ex_dt[, .(start = min(start), end = max(end), y = (ymin + ymax) / 2), .(transcript_id)]
    in_dt$color = intron_color
    if(!is.null(highlight_transcripts)){
      in_dt[transcript_id %in% highlight_transcripts, color := highlight_color]
    }
    p_track_ref = p_track_ref +
      geom_segment(data = in_dt, aes(x = start, xend = end, y = y, yend = y, color = color), size = intron_thickness)
    if(show_tss){
      if(ex_dt$strand[1] == "-"){
        p_track_ref = p_track_ref +
          geom_segment(data = in_dt, aes(x = end, xend = end-1, y = y, yend = y, color = color), size = intron_thickness, arrow = arrow(length = unit(1 / length(transcript_lev) * exon_height * .8, "npc")))
      }else{
        p_track_ref = p_track_ref +
          geom_segment(data = in_dt, aes(x = start, xend = start+1, y = y, yend = y, color = color), size = intron_thickness, arrow = arrow(length = unit(1 / length(transcript_lev) * exon_height * .8, "npc")))
      }
      p_track_ref = p_track_ref + labs(caption = "arrow indicates TSS")

    }
  }
  p_track_ref = p_track_ref +
    geom_rect(data = ex_dt, aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, color = color, fill = color)) +
    scale_color_identity() +
    scale_fill_identity() +
    scale_y_continuous(labels = levels(ex_dt$transcript_id), breaks = seq_along(levels(ex_dt$transcript_id))-.5) +
    theme_classic() +
    theme(panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.y = element_line(size = 1, color = "gray")) +
    labs(title = sel_gene_name, substitle = unique(ex_dt$seqnames)) +
    labs(y = "")

  p_track_ref = .apply_x_scale(p_track_ref, x_scale, as.character(seqnames(query_gr)))
  p_track_ref = .apply_x_lim(p_track_ref, query_gr, flip_x)

  p_track_ref
}
