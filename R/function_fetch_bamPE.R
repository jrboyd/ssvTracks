#' ssvFetchBamPE.RNA
#'
#' @param file_paths
#' @param qgr
#' @param win_size
#' @param target_strand
#' @param splice_strategy
#' @param return_data.table
#' @param win_method
#' @param max_dupes
#' @param flip_strand
#' @param sum_reads
#' @param n_cores
#' @param force_skip_centerFix
#' @param n_region_splits
#'
#' @return
#' @export
#'
#' @examples
ssvFetchBamPE.RNA = function(
    file_paths,
    qgr,
    win_size = 50,
    target_strand = "both",
    splice_strategy = "ignore",
    return_data.table = FALSE,
    win_method = "sample",
    max_dupes = Inf,
    flip_strand = FALSE,
    sum_reads = TRUE,
    n_cores = getOption("mc.cores", 1),
    force_skip_centerFix = TRUE,
    n_region_splits = 1){
  y = cn = NULL #reserve bindings
  strand(qgr) = "*"
  bam_r1 = seqsetvis::ssvFetchBam(
    file_paths = file_paths,
    qgr = qgr,
    target_strand = target_strand,
    splice_strategy = splice_strategy,
    return_data.table = TRUE,
    fragLens = NA,
    win_size = win_size,
    flag = Rsamtools::scanBamFlag(isFirstMateRead = TRUE),
    win_method = win_method,
    flip_strand = !flip_strand,
    max_dupes = max_dupes,
    n_cores = n_cores,
    force_skip_centerFix = force_skip_centerFix,
    n_region_splits = n_region_splits
  )
  cn = colnames(bam_r1)
  bam_r1$read = "r1"
  bam_r2 =
    seqsetvis::ssvFetchBam(
      file_paths = file_paths,
      qgr = qgr,
      target_strand = target_strand,
      splice_strategy = splice_strategy,
      return_data.table = TRUE,
      fragLens = NA,
      win_size = win_size,
      flag = Rsamtools::scanBamFlag(isSecondMateRead = TRUE),
      win_method = win_method,
      flip_strand = flip_strand,
      max_dupes = max_dupes,
      n_cores = n_cores,
      force_skip_centerFix = force_skip_centerFix,
      n_region_splits = n_region_splits
    )
  bam_r2$read = "r2"

  if(sum_reads){
    bam_dt = rbind(bam_r1, bam_r2)[, cn, with = FALSE]
    bam_dt = bam_dt[, list(y = sum(y)), by = c(cn[cn != "y"])][, cn, with = FALSE]
  }else{
    bam_dt = rbind(bam_r1, bam_r2)
  }

  if(!return_data.table){
    bam_dt = GRanges(bam_dt)
  }
  bam_dt
}

#' ssvFetchBamPE.RNA_splice
#'
#' @param file_paths
#' @param qgr
#' @param win_size
#' @param target_strand
#' @param return_data.table
#' @param win_method
#' @param max_dupes
#' @param flip_strand
#' @param n_cores
#' @param force_skip_centerFix
#' @param n_region_splits
#'
#' @return
#' @export
#'
#' @examples
ssvFetchBamPE.RNA_splice = function(
    file_paths,
    qgr, win_size = 50,
    target_strand = "both",
    return_data.table = FALSE,
    win_method = "sample",
    max_dupes = Inf,
    flip_strand = FALSE,
    n_cores = getOption("mc.cores", 1),
    force_skip_centerFix = TRUE,
    n_region_splits = 1,
    ...){
  env = as.list(sys.frame(sys.nframe()))
  final_target_strand = target_strand
  args = c(as.list(env), list(...))
  args$return_data.table = TRUE
  args$splice_strategy = "splice_count"
  args$sum_reads = FALSE
  args$target_strand = "both"
  splice_dt = do.call(ssvFetchBamPE.RNA, args)

  splice_dt[read == "r1" & strand == "+", strand := "tmp"]
  splice_dt[read == "r1" & strand == "-", strand := "+"]
  splice_dt[read == "r1" & strand == "tmp", strand := "-"]
  cn = setdiff(colnames(splice_dt), c("read", "N"))
  splice_dt = splice_dt[, .(N = sum(N)), c(cn)]

  if(final_target_strand %in% c("-", "+")){
    splice_dt = splice_dt[strand == final_target_strand]
  }
  # if(!flip_strand){
  #   splice_dt[strand == "+", strand := "tmp"]
  #   splice_dt[strand == "-", strand := "+"]
  #   splice_dt[strand == "tmp", strand := "-"]
  # }

  splice_dt[]
}

