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
  
  a_and_b$a_adj = ifelse(is.na(a_and_b[,condition_a]), 0.99, a_and_b[,condition_a])
  a_and_b$b_adj = ifelse(is.na(a_and_b[,condition_b]), 0.99, a_and_b[,condition_b])
  
  #plot a v b
  png(sprintf("%s.png", comparison))
  plot(log(a_and_b$a_adj, 2), log(a_and_b$b_adj, 2), pch=16, col="#808D8D4A", 
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

write.csv(a_and_b, sprintf("%s.csv", comparison))
}

#a_v_b_plots = function(location_of_file_a, condition_a, 
#               location_of_file_b, condition_b, 
#               comparison){

getwd()
setwd("/Users/annamcgeachy/Google Drive/post trx reg data/datafiles_screen2_miseq_repiped//")

#an example
a_v_b_plots("screen2-mi-down.csv", "screen2-mi-down", 
            "screen2-mi-up.csv", "screen2-mi-up", 
            "screen2-mi-down-v-up")

a_v_b_plots("screen2-mi-down.csv", "screen2-mi-down", 
            "screen2-mi-norecomb.csv", "screen2-mi-norecomb", 
            "screen2-mi-down-v-norecomb")

a_v_b_plots("screen2-mi-down.csv", "screen2-mi-down", 
            "screen2-mi-unsort.csv", "screen2-mi-unsort", 
            "screen2-mi-down-v-unsort")

##  

a_v_b_plots("screen2-mi-up.csv", "screen2-mi-up", 
            "screen2-mi-norecomb.csv", "screen2-mi-norecomb", 
            "screen2-mi-up-v-norecomb")

a_v_b_plots("screen2-mi-up.csv", "screen2-mi-up", 
            "screen2-mi-unsort.csv", "screen2-mi-unsort", 
            "screen2-mi-up-v-unsort")
##

a_v_b_plots("screen2-mi-norecomb.csv", "screen2-mi-norecomb", 
            "screen2-mi-unsort.csv", "screen2-mi-unsort", 
            "screen2-mi-norecomb-v-unsort")