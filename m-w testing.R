getwd()
setwd("/Users/annamcgeachy/Google Drive/post trx reg data/datafiles_screen1_and_2_hiseq_repiped/")

library("dplyr")

#read in files
up = read.table("screen1-hi-up.csv", sep=",", stringsAsFactors = FALSE, header=TRUE)
head(up)
down = read.table("screen1-hi-down.csv", sep=",", stringsAsFactors = FALSE, header=TRUE)
head(down)
xref = read.delim("../SGD_features.tab", header=FALSE, quote="")

#filter lists to get only in frame fragments
up_inframe = filter(up, joint_frame=="0,1")
down_inframe = filter(down, joint_frame=="0,1")

#make a list of all the gene names
up_inframe_genes = as.data.frame(table(up_inframe$gene_name), ncol=2)
#because we used stringsAsFactors=FALSE above, we no longer have that bullshit extra genes as 0s because they were factors

head(up_inframe) # first gene is YAL019W so try that; NOTE in the miseq the first gene was YAL017W

#run the mw test
wilcox.test(filter(up_inframe, gene_name=="YAL019W")[,"frag_count"], filter(up_inframe, gene_name!="YAL019W")[,"frag_count"],)

  #try returning this to an object so I can grab the 
  blorb = wilcox.test(filter(up_inframe, gene_name=="YAL019W")[,"frag_count"], filter(up_inframe, gene_name!="YAL019W")[,"frag_count"],) #it does
  typeof(blorb)
  unlist(blorb)

  blorb[3][[1]] #how to get the pvalue out of the fucking ridiculous wilcox list

#so now we have a way to use a gene name to get a wilcox p value. 
#but how do we do this for a list of gene names? not just a given one.

#how to get a gene name out of the in frame gene list
head(up_inframe_genes)
up_inframe_genes$Var1[1]
#test this using beta substitution
  test_thing = wilcox.test(
    filter(up_inframe, gene_name==as.character(up_inframe_genes$Var1[1]))[,"frag_count"],
    filter(up_inframe, gene_name!=as.character(up_inframe_genes$Var1[1]))[,"frag_count"])
    #still have to use the as.character because for some idiotic reason it's calling them integers

  test_thing
  
  blorb2 = wilcox.test(filter(up_inframe, gene_name=="Q0045")[,"frag_count"], filter(up_inframe, gene_name!="Q0055")[,"frag_count"],) #it does
  blorb2 #make sure that works the same as above; it does

#fetch the p value from this mess
  test_thing_p = wilcox.test(
    filter(up_inframe, gene_name==as.character(up_inframe_genes$Var1)[1])[,"frag_count"],
    filter(up_inframe, gene_name!=as.character(up_inframe_genes$Var1)[1])[,"frag_count"])[3][[1]]
  test_thing_p

#can we use apply to get this to work? instead of a for loop
  just_one = up_inframe_genes[1,]
  
just_one_test_thing_p = wilcox.test(
  filter(up_inframe, gene_name==as.character(just_one$Var1))[,"frag_count"],
  filter(up_inframe, gene_name!=as.character(just_one$Var1))[,"frag_count"])[3][[1]]

  just_one_test_thing_p

#function to do pvalues using apply and mw test
mw_test_up = function(x){
  wilcox.test(
    filter(up_inframe, gene_name==as.character(x["Var1"]))[,"frag_count"],
    filter(up_inframe, gene_name!=as.character(x["Var1"]))[,"frag_count"])[3][[1]]
}

#generate the pvalues, then adjust
up_p_vals = apply(up_inframe_genes, 1, mw_test_up)
up_p_vals_adj = p.adjust(up_p_vals)

#histogram of pvals pre and post adjustment
hist(up_p_vals, col=2, add=TRUE)
hist(up_p_vals_adj, add=TRUE, col=3)

#add the pvalues to the data frame
up_inframe_genes$mw_p_vals = up_p_vals
up_inframe_genes$p_adj = up_p_vals_adj

head(up_inframe_genes)

#reorder by smallest adjusted p value
up_inframe_genes = up_inframe_genes[order(up_inframe_genes$p_adj),]
table(up_inframe_genes$p_adj)

#add in useful gene information 
up_inframe_genes$gene_useful = xref[match(up_inframe_genes$Var1, xref$V4),"V5"]
up_inframe_genes$gene_desc = xref[match(up_inframe_genes$Var1, xref$V4),"V16"]

head(up_inframe_genes)

#########

#make a list of all the gene names
down_inframe_genes = as.data.frame(table(down_inframe$gene_name), ncol=2)

#write the mw function
mw_test_down = function(x){
  wilcox.test(
    filter(down_inframe, gene_name==as.character(x["Var1"]))[,"frag_count"],
    filter(down_inframe, gene_name!=as.character(x["Var1"]))[,"frag_count"])[3][[1]]
}

#generate the pvalues, then adjust

down_p_vals = apply(down_inframe_genes, 1, mw_test_down)
down_p_vals_adj = p.adjust(down_p_vals)

#histogram of pvals pre and post adjustment
hist(down_p_vals_adj, col=3)
hist(down_p_vals, add=TRUE, col=2)

#add the pvalues to the data frame
down_inframe_genes$mw_p_vals = down_p_vals
down_inframe_genes$p_adj = down_p_vals_adj

#reorder by smallest adjusted p value
down_inframe_genes = down_inframe_genes[order(down_inframe_genes$p_adj),]
table(down_inframe_genes$p_adj)

#add in useful gene information 
down_inframe_genes$gene_useful = xref[match(down_inframe_genes$Var1, xref$V4),"V5"]
down_inframe_genes$gene_desc = xref[match(down_inframe_genes$Var1, xref$V4),"V16"]

head(down_inframe_genes)

####
#try to do this for up-v-down

#read in table that already has counts for each of down and up (in frame only)
down_v_up = read.csv("screen1-hi-down-v-up.csv", stringsAsFactors=FALSE)

#need to relate back to gene names
#to do this, need to get them from down_inframe and up_inframe
#but need to give those _inframe the same type of unique identifiers
#build them the same was as we did for down_v_up
pasting_together = function(x){
  paste(x["chr_read"], as.numeric(x["start_read"]), as.numeric(x["end_read"]), x["strand_read"], sep="_")
}

down_inframe$uniques = apply(down_inframe, 1, pasting_together)
up_inframe$uniques = apply(up_inframe, 1, pasting_together)

#then cross reference to pull gene names
down_v_up$down_genes = down_inframe[match(down_v_up$a_and_b, down_inframe$uniques), "gene_name"]
down_v_up$up_genes = up_inframe[match(down_v_up$a_and_b, up_inframe$uniques), "gene_name"]
head(down_v_up)

#make a column that has gene names (filling in the blanks where it's only in one versus the other)
down_v_up$gene_name = ifelse(is.na(down_v_up$down_genes), as.character(down_v_up$up_genes), as.character(down_v_up$down_genes))
nrow(down_v_up)
sum(!is.na(down_v_up$gene_name))

  #make a list of genes for use in the mw later
  down_v_up_gene_names = as.data.frame(table(down_v_up$gene_name))
  head(down_v_up_gene_names)
  nrow(down_v_up_gene_names)

#need to make a single measurement for use in MW 
#that is, can't use two sets of numbers (up frags and down frags)
#need it to be one set of numbers
#do ratio between them using the pseudocount (otherwise you div by 0 and nope)

down_v_up$down_div_up = down_v_up$a_adj / down_v_up$b_adj
down_v_up$up_div_down = down_v_up$b_adj / down_v_up$a_adj

#now make the mw for it
mw_test_dvu_d = function(x){
  wilcox.test(
    filter(down_v_up, gene_name==(x["Var1"]))[,"down_div_up"],
    filter(down_v_up, gene_name!=(x["Var1"]))[,"down_div_up"])[3][[1]]
}

down_v_up_gene_names$down_div_up_pvals = apply(down_v_up_gene_names, 1, mw_test_dvu_d)
down_v_up_gene_names$down_div_up_padj = p.adjust(down_v_up_gene_names$down_div_up_pvals)
hist(down_v_up_gene_names$down_div_up_pvals)
hist(down_v_up_gene_names$down_div_up_padj, col=2, add=TRUE)

table(down_v_up_gene_names$down_div_up_padj)
head(down_v_up)

##now do it for up div down
mw_test_dvu_u = function(x){
  wilcox.test(
    filter(down_v_up, gene_name==(x["Var1"]))[,"up_div_down"],
    filter(down_v_up, gene_name!=(x["Var1"]))[,"up_div_down"])[3][[1]]
}

down_v_up_gene_names$up_div_down_pvals = apply(down_v_up_gene_names, 1, mw_test_dvu_u)
down_v_up_gene_names$up_div_down_padj = p.adjust(down_v_up_gene_names$up_div_down_pvals)

hist(down_v_up_gene_names$up_div_down_pvals)
hist(down_v_up_gene_names$up_div_down_padj, add=TRUE)

table(down_v_up_gene_names$up_div_down_padj) #get the exact same thing as the reverse, 1948 at 1 and 1 at .97795
#fuck

####
#ok. what if we try m-w test per fragment

head(down_v_up)

mw_test_frag = function(x){
  wilcox.test(
  filter(down_v_up, a_and_b==(x["Var1"]))[,"up_div_down"],
  filter(down_v_up, a_and_b!=(x["Var1"]))[,"up_div_down"])[3][[1]]
}
uniques = as.data.frame(table(down_v_up$a_and_b))

blab = apply(uniques, 1, mw_test_frag)
hist(p.adjust(blab))

###

down_inframe_bi = filter(down_inframe, frag_count>64)

down_inframe_bi_genes = as.data.frame(table(down_inframe_bi$gene_name))
head(down_inframe_bi_genes)
down_inframe_bi_genes = filter(down_inframe_bi_genes, Freq!=0)
head(down_inframe_bi_genes)

mw_test_down_bi = function(x){
  wilcox.test(
    filter(down_inframe_bi, gene_name==as.character(x["Var1"]))[,"frag_count"],
    filter(down_inframe_bi, gene_name!=as.character(x["Var1"]))[,"frag_count"])[3][[1]]
}

down_inframe_bi_genes$pvals = apply(down_inframe_bi_genes, 1, mw_test_down_bi)
down_inframe_bi_genes$padj = p.adjust(down_inframe_bi_genes$pvals)

hist(down_inframe_bi_genes$padj)

head(down_inframe)
as.matrix(filter(down_inframe, gene_name=="YER165W")[,c("uniques","frag_count")])
as.matrix(filter(up_inframe, gene_name=="YER165W")[,c("uniques","frag_count")])

filter(up_inframe, gene_name=="YBL113C")

?max

YPR204W = filter(up_inframe, gene_name=="YPR204W")
YPR204W = YPR204W[order(YPR204W$frag_count, decreasing = TRUE),]
head(YPR204W)

?pmax
