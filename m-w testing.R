getwd()
setwd("/Users/annamcgeachy/Google Drive/post trx reg data/screen1_miseq_repiped/")

library("dplyr")

up = read.csv("screen1_up.csv")
down = read.csv("screen1_down.csv")

up_inframe = filter(up, joint_frame=="0,1")
down_inframe = filter(down, joint_frame=="0,1")

head(up_inframe)

#make a list of all the gene names
up_inframe_genes = as.data.frame(table(up_inframe$gene_name), ncol=2)

#notice that there are a lot of 0s (so they occur 0 times)
table(up_inframe_genes$Freq)
  #so where are they coming from
  head(up_inframe_genes)
  grep("Q0045", up_inframe$gene_name) #they in face don't exist in up_inframe
  grep("Q0045", up$gene_name) #looks like they existed in up,
  up[1065,] # but are out of frame

  #so how do we get rid of these 0's for the in frame so they don't complicate things
  up_inframe_genes_actual = filter(up_inframe_genes, Freq!=0)
  nrow(up_inframe_genes_actual)
  table(up_inframe_genes_actual$Freq) #ok. got all of the 0s removed.

head(up_inframe) # first gene is YAL017W so try that
  #make list of all YAL017W and all of NOT YAL017W
  YAL017W_yes = filter(up_inframe, gene_name=="YAL017W")
  YAL017W_no = filter(up_inframe, gene_name!="YAL017W")
  
  #check to see if the beta substitution works as well as the generated lists
  wilcox.test(filter(up_inframe, gene_name=="YAL017W")[,"frag_count"], filter(up_inframe, joint_frame=="0,1", gene_name!="YAL017W")[,"frag_count"],)
  wilcox.test(YAL017W_yes$frag_count, YAL017W_no$frag_count) #it does




table(filter(up, joint_frame=="0,1", gene_name=="YPR204W")[,"frag_count"])
table(filter(down, joint_frame=="0,1", gene_name=="YCR077C")[,"frag_count"])
nrow(as.matrix(table(up$gene_name)))

colnames(up)
table(up$joint_frame)

down_v_up = read.csv("screen1-down-v-up.csv")
head(down_v_up)
