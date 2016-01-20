getwd()
setwd("/Users/annamcgeachy/Google Drive/post trx reg data/screen1_miseq_repiped/")

library("dplyr")

#read in files
up = read.csv("screen1_up.csv")
down = read.csv("screen1_down.csv")
xref = read.delim("../SGD_features.tab", header=FALSE, quote="")

#filter lists to get only in frame fragments
up_inframe = filter(up, joint_frame=="0,1")
down_inframe = filter(down, joint_frame=="0,1")

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

  #try returning this to an object so I can grab the 
  blorb = wilcox.test(filter(up_inframe, gene_name=="YAL017W")[,"frag_count"], filter(up_inframe, joint_frame=="0,1", gene_name!="YAL017W")[,"frag_count"],) #it does
  typeof(blorb)
  unlist(blorb)

  blorb[3][[1]] #how to get the pvalue out of the fucking ridiculous wilcox list

#so now we have a way to use a gene name to get a wilcox p value. 
#but how do we do this for a list of gene names? not just a given one.

#how to get a gene name out of the in frame gene list
as.character(up_inframe_genes_actual$Var1)[1]

#test this using beta substitution
  test_thing = wilcox.test(
    filter(up_inframe, gene_name==as.character(up_inframe_genes_actual$Var1)[1])[,"frag_count"],
    filter(up_inframe, gene_name!=as.character(up_inframe_genes_actual$Var1)[1])[,"frag_count"])
  
  test_thing
  
  blorb2 = wilcox.test(filter(up_inframe, gene_name=="Q0055")[,"frag_count"], filter(up_inframe, gene_name!="Q0055")[,"frag_count"],) #it does
  blorb2 #make sure that works the same as above; it does

#fetch the p value from this mess
  test_thing_p = wilcox.test(
    filter(up_inframe, gene_name==as.character(up_inframe_genes_actual$Var1)[1])[,"frag_count"],
    filter(up_inframe, gene_name!=as.character(up_inframe_genes_actual$Var1)[1])[,"frag_count"])[3][[1]]
  test_thing_p

#can we use apply to get this to work? instead of a for loop
  just_one = up_inframe_genes_actual[1,]
  
just_one_test_thing_p = wilcox.test(
  filter(up_inframe, gene_name==as.character(just_one$Var1))[,"frag_count"],
  filter(up_inframe, gene_name!=as.character(just_one$Var1))[,"frag_count"])[3][[1]]

#function to do pvalues using apply and mw test
mw_test = function(x){
  wilcox.test(
    filter(up_inframe, gene_name==as.character(x["Var1"]))[,"frag_count"],
    filter(up_inframe, gene_name!=as.character(x["Var1"]))[,"frag_count"])[3][[1]]
}

#generate the pvalues, then adjust
up_p_vals = apply(up_inframe_genes_actual, 1, mw_test)
up_p_vals_adj = p.adjust(up_p_vals)

#histogram of pvals pre and post adjustment
hist(up_p_vals, col=2)
hist(up_p_vals_adj, add=TRUE, col=3)

#add the pvalues to the data frame
up_inframe_genes_actual$mw_p_vals = up_p_vals
up_inframe_genes_actual$p_adj = up_p_vals_adj

head(up_inframe_genes_actual)

#reorder by smallest adjusted p value
up_inframe_genes_actual = up_inframe_genes_actual[order(up_inframe_genes_actual$p_adj),]
table(up_inframe_genes_actual$p_adj)

#add in useful gene information 
up_inframe_genes_actual$gene_useful = xref[match(up_inframe_genes_actual$Var1, xref$V4),"V5"]
up_inframe_genes_actual$gene_desc = xref[match(up_inframe_genes_actual$Var1, xref$V4),"V16"]

head(up_inframe_genes_actual)
up_inframe_genes_actual[grep("PAT1", up_inframe_genes_actual$gene_useful),]


#########

#make a list of all the gene names
down_inframe_genes = as.data.frame(table(down_inframe$gene_name), ncol=2)

#notice that there are a lot of 0s (so they occur 0 times)
table(down_inframe_genes$Freq)
#so where are they coming from
head(down_inframe_genes)
grep("Q0050", down_inframe$gene_name) #they in face don't exist in down_inframe
grep("Q0050", down$gene_name) #looks like they existed in down,
down[576,] # but are out of frame

#so how do we get rid of these 0's for the in frame so they don't complicate things
down_inframe_genes_actual = filter(down_inframe_genes, Freq!=0)
nrow(down_inframe_genes_actual)
table(down_inframe_genes_actual$Freq) #ok. got all of the 0s removed.
head(down_inframe_genes_actual)
head(up_in)

#
mw_test = function(x){
  wilcox.test(
    filter(down_inframe, gene_name==as.character(x["Var1"]))[,"frag_count"],
    filter(down_inframe, gene_name!=as.character(x["Var1"]))[,"frag_count"])[3][[1]]
}

#generate the pvalues, then adjust
down_p_vals = apply(down_inframe_genes_actual, 1, mw_test)
down_p_vals_adj = p.adjust(down_p_vals)

#histogram of pvals pre and post adjustment
hist(down_p_vals, col=2)
hist(down_p_vals_adj, add=TRUE, col=3)
table(down_p_vals_adj)

#add the pvalues to the data frame
down_inframe_genes_actual$mw_p_vals = down_p_vals
down_inframe_genes_actual$p_adj = down_p_vals_adj

#reorder by smallest adjusted p value
down_inframe_genes_actual = down_inframe_genes_actual[order(down_inframe_genes_actual$p_adj),]
table(down_inframe_genes_actual$p_adj)

#add in useful gene information 
down_inframe_genes_actual$gene_useful = xref[match(down_inframe_genes_actual$Var1, xref$V4),"V5"]
down_inframe_genes_actual$gene_desc = xref[match(down_inframe_genes_actual$Var1, xref$V4),"V16"]

head(down_inframe_genes_actual)
down_inframe_genes_actual[grep("PAT1", down_inframe_genes_actual$gene_useful),]

