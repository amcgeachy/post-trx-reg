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


#####
#generalized version

getwd()
setwd("/Users/annamcgeachy/Google Drive/post trx reg data/datafiles_screen4/")

xref = read.delim("../SGD_features.tab", header=FALSE, quote="")

#make unique identifiers for each fragment
pasting_together = function(x){
  paste(x["chr_read"], as.numeric(x["start_read"]), as.numeric(x["end_read"]), x["strand_read"], sep="_")
}

a_v_b_plots = function(location_of_file_a, condition_a, location_of_file_b, condition_b, comparison){
  ##set up A
  general_a = read.csv(location_of_file_a)
  
  general_a_uniques = apply(general_a, 1, pasting_together)
  general_a$uniques = general_a_uniques
  
  #add in the annotation
  general_a$gene_useful = xref[match(general_a$gene_name, xref$V4),"V5"]
  general_a$gene_desc = xref[match(general_a$gene_name, xref$V4),"V16"]
  
  #pull in frame fragments
  general_a_01 = filter(general_a, joint_frame=="0,1")
  
  ##set up B
  general_b = read.csv(location_of_file_b)
  
  general_b_uniques = apply(general_b, 1, pasting_together)
  general_b$uniques = general_b_uniques
  
  #add in the annotation
  general_b$gene_useful = xref[match(general_b$gene_name, xref$V4),"V5"]
  general_b$gene_desc = xref[match(general_b$gene_name, xref$V4),"V16"]
  
  #pull in frame fragments
  general_b_01 = filter(general_b, joint_frame=="0,1")
  
  head(general_a_01)
  head(general_b_01)
  
  #make joint file of all uniques between libraries
  a_and_b = c(as.character(general_a_01$uniques), as.character(general_b_01$uniques))
  length(a_and_b)
  head(a_and_b)
  a_and_b = unique(a_and_b)
  a_and_b = as.data.frame(a_and_b)
  nrow(a_and_b)
  
  #pull in fragment counts in each library
  head(a_and_b)
  a_and_b[,condition_a] = general_a_01[match(a_and_b$a_and_b, general_a_01$uniques), "frag_count"]
  table(a_and_b[,condition_a])
  a_and_b[,condition_b] = general_b_01[match(a_and_b$a_and_b, general_b_01$uniques), "frag_count"]
  table(a_and_b[,condition_b])
  
  #plot a v b
  png(sprintf("%s.png", comparison))
  plot(log(a_and_b[,condition_a], 2), log(a_and_b[,condition_b], 2), pch=16, col="#808D8D4A", 
#        xlim=max(max(log(a_and_b[,condition_a], 2)), max(log(a_and_b[,condition_b], 2))),
#        ylim=max(max(log(a_and_b[,condition_a], 2)), max(log(a_and_b[,condition_b], 2))),
       xlab=sprintf("log2(%s)", condition_a),
       ylab=sprintf("log2(%s)", condition_b),
       main=sprintf("%s", comparison))
  abline(a=0,b=1)
  abline(a=-4.6,b=1)
  abline(a=4.6,b=1)
  abline(a=-2.3,b=1)
  abline(a=2.3,b=1)
  dev.off()
}
max(max(c(1,2)), max(c(3,4)))

#a_v_b_plots = function(location_of_file_a, condition_a, 
#               location_of_file_b, condition_b, 
#               comparison){
a_v_b_plots("./datafiles_screen4_hiseq/hi_one_full_down.csv", "1-100N-down", 
            "./datafiles_screen4_hiseq/hi_one_full_up.csv", "1-100N-up", 
            "1-100N-up-v-down")