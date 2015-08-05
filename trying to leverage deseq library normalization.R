library(DESeq2)

??DESeq2
getwd()
setwd("/Users/annamcgeachy/Google Drive/post-trx reg/datafiles_new")

frag_deseq = read.csv("all unique frag with counts per sample.csv", header = TRUE)


colnames(frag_deseq)
frag_play = data.frame(
  row.names = colnames(frag_deseq[,3:6]),
  recombination = c("recomb", "recomb", "no_recomb", "recomb"),
  up_pop = c("up", "no_up", "up", "up"),
  down_pop = c("no_down", "down", "down", "down"),
  up_selected =c("up_selected", "not_up_selection", "not_up_selection", "not_up_selection"),
  down_selected= c("not_down_selected", "down_selected", "not_down_selected", "not_down_selected")
)
frag_play

frag_test = data.frame(
  row.names=colnames(frag_deseq[,3:6]),
  sample=c("up", "down", "no", "post"))

frag_play
library(DESeq)
cds_3way_nopseu_cut = newCountDataSet(frag_deseq[,3:6], frag_play)
cds_3way_nopseu_cut = estimateSizeFactors(cds_3way_nopseu_cut)
sizeFactors(cds_3way_nopseu_cut)

other_test = newCountDataSet(frag_deseq[,3:6], frag_test)
other_test = estimateSizeFactors(other_test)
sizeFactors(other_test) #doesnt change even with changing the desc data frame


ncounts <- t( t( counts(other_test) ) / sizeFactors(other_test) )
plot( 
  ( ncounts[,2] + ncounts[,1] )/2,
  ncounts[,2] / ncounts[,1],
  log="xy", pch="." )
abline( h = 1)

means <- rowMeans( ncounts )
plot( 
  means,
  ncounts[,2] / means,
  log="xy", pch="." )
abline( h = 1)