getwd()
setwd("/Users/annamcgeachy/Google Drive/post trx reg data/datafiles_screen4/")

xref = read.delim("../SGD_features.tab", header=FALSE, quote="")

#make unique identifiers for each fragment
pasting_together = function(x){
  paste(x["chr_read"], as.numeric(x["start_read"]), as.numeric(x["end_read"]), x["strand_read"], sep="_")
}


##set up A
  one_full_down_hi = read.csv("./datafiles_screen4_hiseq/hi_one_full_down")
  
  one_full_down_hi_uniques = apply(one_full_down_hi, 1, pasting_together)
  one_full_down_hi$uniques = one_full_down_hi_uniques
  
  #add in the annotation
  one_full_down_hi$gene_useful = xref[match(one_full_down_hi$gene_name, xref$V4),"V5"]
  one_full_down_hi$gene_desc = xref[match(one_full_down_hi$gene_name, xref$V4),"V16"]
  
  #pull in frame fragments
  one_full_down_hi_01 = filter(one_full_down_hi, joint_frame=="0,1")

##set up B
  one_full_down_mi = read.csv("./datafiles_screen4_miseq/one_full_down.csv")
  
  one_full_down_mi_uniques = apply(one_full_down_mi, 1, pasting_together)
  one_full_down_mi$uniques = one_full_down_mi_uniques
  
  #add in the annotation
  one_full_down_mi$gene_useful = xref[match(one_full_down_mi$gene_name, xref$V4),"V5"]
  one_full_down_mi$gene_desc = xref[match(one_full_down_mi$gene_name, xref$V4),"V16"]
  
  #pull in frame fragments
  one_full_down_mi_01 = filter(one_full_down_mi, joint_frame=="0,1")


head(one_full_down_hi_01)
head(one_full_down_mi_01)

ones = c(as.character(one_full_down_hi_01$uniques), as.character(one_full_down_mi_01$uniques))
length(ones)
head(ones)
ones = unique(ones)
ones = as.data.frame(ones)
nrow(ones)

head(ones)
ones$hi_frags = one_full_down_hi_01[match(ones$ones, one_full_down_hi_01$uniques), "frag_count"]
table(ones$hi_frags)
ones$mi_frags = one_full_down_mi_01[match(ones$ones, one_full_down_mi_01$uniques), "frag_count"]
table(ones$mi_frags)
ones$down_frags = one_full_down_01_gene_counts[match(ones$ones, one_full_down_01_gene_counts$Var1), "Freq"]

png("1-100N hi v mi frag count - logged.png")
plot(log(ones$hi_frags, 2), log(ones$mi_frags, 2), pch=16, col="#808D8D4A", 
     xlab="log2(unique fragment counts in hi_seq library)",
     ylab="log2(unique fragment counts in mi_seq library)",
     main="unique fragments counts 1-100N down")
dev.off()