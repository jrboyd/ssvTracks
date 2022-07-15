library(magrittr)
library(GenomicRanges)
library(seqsetvis)
library(data.table)
bam_files = dir("~/R_workspace.t630/", pattern = "RNAseq_KT_Runx1_KO_MCF10A", full.names = TRUE) %>%
  dir(full.names = TRUE) %>%
  dir(pattern = "bam$", full.names = TRUE)

ref_gr = rtracklayer::import.gff("~/gencode.v36.annotation.gtf", feature.type = "gene")
query_gr = range(subset(ref_gr, gene_name == "ESR1"))

options(mc.cores = 1)
all_prof = pbmcapply::pbmclapply(bam_files, mc.cores = 20, function(f){
  ssvRecipes::ssvFetchBamPE.RNA(f, query_gr, win_size = 2000, win_method = "summary", return_data.table = TRUE, sum_reads = FALSE)
})

k = sapply(all_prof, is, class2 = "data.table")
prof_dt = rbindlist(all_prof[k])

sum(k)

todo = levels(prof_dt$sample)
todo_groups = split(levels(prof_dt$sample), ceiling(seq_along(todo)/9))

td = todo_groups[[2]]
pdf("tmp.pdf", width = 6, height = 6)
for(i in seq_along(todo_groups)){
  td = todo_groups[[i]]
  p = ggplot(prof_dt[sample %in% td], aes(x = x, y = y, color = strand)) +
    geom_path() +
    facet_grid(sample~strand+read)
  plot(p)
}
dev.off()
