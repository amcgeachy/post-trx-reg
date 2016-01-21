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
mw_test_up = function(x){
  wilcox.test(
    filter(up_inframe, gene_name==as.character(x["Var1"]))[,"frag_count"],
    filter(up_inframe, gene_name!=as.character(x["Var1"]))[,"frag_count"])[3][[1]]
}

#generate the pvalues, then adjust
up_p_vals = apply(up_inframe_genes_actual, 1, mw_test_up)
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

#write the mw function
mw_test_down = function(x){
  wilcox.test(
    filter(down_inframe, gene_name==as.character(x["Var1"]))[,"frag_count"],
    filter(down_inframe, gene_name!=as.character(x["Var1"]))[,"frag_count"])[3][[1]]
}

#generate the pvalues, then adjust

down_p_vals = apply(down_inframe_genes_actual, 1, mw_test_down)
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

####
#try to do this for up-v-down

#read in table that already has counts for each of down and up (in frame only)
down_v_up = read.csv("screen1-down-v-up.csv", stringsAsFactors=FALSE)

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
head(down_v_up)
sum(!is.na(down_v_up$down_genes))

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
mw_test_general = function(input_data_frame, gene_list, condition){
  wilcox.test(
    filter(input_data_frame, gene_name==as.character(gene_list["gene_name"]))[,condition],
    filter(input_data_frame, gene_name!=as.character(gene_list["gene_name"]))[,condition])[3][[1]]
}

head(down_v_up)



down_v_up$down_div_up_p_val = apply(down_v_up_gene_names, 1, mw_test_dvu)
down_v_up_gene_names[1]
wilcox.test(filter(down_v_up, gene_name==down_v_up_gene_names[1])[,"down_div_up"],
            filter(down_v_up, gene_name!=down_v_up_gene_names[1])[,"down_div_up"])[3][[1]]

mw_test_down = function(x){
  wilcox.test(
    filter(down_inframe, gene_name==as.character(x["Var1"]))[,"frag_count"],
    filter(down_inframe, gene_name!=as.character(x["Var1"]))[,"frag_count"])[3][[1]]
}

mw_test_dvu = function(x){
  wilcox.test(
    filter(down_v_up, gene_name==(x["Var1"]))[,"down_div_up"],
    filter(down_v_up, gene_name!=(x["Var1"]))[,"down_div_up"])[3][[1]]
}

blorb = apply(down_v_up_gene_names, 1, mw_test_dvu)
head(blorb)

head(down_v_up_gene_names)

wilcox.test(
  filter(down_v_up, gene_name=="YAL002W")[,"down_div_up"],
  filter(down_v_up, gene_name!="YAL002W")[,"down_div_up"])[3][[1]]
