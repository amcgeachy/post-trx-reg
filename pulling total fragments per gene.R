#make unique identifiers for each fragment
pasting_together = function(x){
  paste(x["chr_read"], as.numeric(x["start_read"]), as.numeric(x["end_read"]), x["strand_read"], sep="_")
}

##one full up
one_full_up = read.csv("one_full_up.csv")

one_full_up_uniques = apply(one_full_up, 1, pasting_together)
one_full_up$uniques = one_full_up_uniques

#add in the annotation
one_full_up$gene_useful = xref[match(one_full_up$gene_name, xref$V4),"V5"]
one_full_up$gene_desc = xref[match(one_full_up$gene_name, xref$V4),"V16"]

#pull in frame fragments
one_full_up_01 = filter(one_full_up, joint_frame=="0,1")
nrow(one_full_up_01)
head(one_full_up_01)

#make a table of the genes
one_full_up_01_gene_counts = as.data.frame(table(one_full_up_01$gene_name))
head(one_full_up_01_gene_counts)

#add in easier to understand names and desc 
one_full_up_01_gene_counts$gene_useful = xref[match(one_full_up_01_gene_counts$Var1, xref$V4),"V5"]
one_full_up_01_gene_counts$gene_desc = xref[match(one_full_up_01_gene_counts$Var1, xref$V4),"V16"]
head(one_full_up_01_gene_counts)

#reorder so its from highest to lowest
one_full_up_01_gene_counts = one_full_up_01_gene_counts[order(one_full_up_01_gene_counts$Freq, decreasing = TRUE),]
write.csv(one_full_up_01_gene_counts, "one_full_up_01_gene_counts.csv")

##one full down
one_full_down = read.csv("one_full_down.csv")

one_full_down_uniques = apply(one_full_down, 1, pasting_together)
one_full_down$uniques = one_full_down_uniques

#add in the annotation
one_full_down$gene_useful = xref[match(one_full_down$gene_name, xref$V4),"V5"]
one_full_down$gene_desc = xref[match(one_full_down$gene_name, xref$V4),"V16"]

#pull in frame fragments
one_full_down_01 = filter(one_full_down, joint_frame=="0,1")
nrow(one_full_down_01)
head(one_full_down_01)

#make a table of the genes
one_full_down_01_gene_counts = as.data.frame(table(one_full_down_01$gene_name))
head(one_full_down_01_gene_counts)


#make joint up and down table
ones = c(as.character(one_full_down_01_gene_counts$Var1), as.character(one_full_up_01_gene_counts$Var1))
length(ones)
head(ones)
ones = unique(ones)
ones = as.data.frame(ones)
head(ones)
ones$up_frags = one_full_up_01_gene_counts[match(ones$ones, one_full_up_01_gene_counts$Var1), "Freq"]
ones$down_frags = one_full_down_01_gene_counts[match(ones$ones, one_full_down_01_gene_counts$Var1), "Freq"]
ones$gene_useful = xref[match(ones$ones, xref$V4),"V5"]
ones$gene_desc = xref[match(ones$ones, xref$V4),"V16"]
head(ones)


table(ones$up_frags)
ones_down_only = filter(ones, up_frags==0, down_frags!=0)
head(ones_down_only)
nrow(ones_down_only)

ones_up_only = filter(ones, up_frags!=0, down_frags==0)
head(ones_up_only)
nrow(ones_up_only)
write.csv(ones_up_only, "ones_up_only.csv")
write.csv(ones_down_only, "ones_down_only.csv")

hist(ones$up_frags)
hist(ones$down_frags)

hist(ones_up_only$up_frags)
hist(ones_down_only$down_frags)
hist(ones$up_frags)


#add in easier to understand names and desc 
one_full_down_01_gene_counts$gene_useful = xref[match(one_full_down_01_gene_counts$Var1, xref$V4),"V5"]
one_full_down_01_gene_counts$gene_desc = xref[match(one_full_down_01_gene_counts$Var1, xref$V4),"V16"]
head(one_full_down_01_gene_counts)

#reorder so its from highest to lowest
one_full_down_01_gene_counts = one_full_down_01_gene_counts[order(one_full_down_01_gene_counts$Freq, decreasing = TRUE),]
write.csv(one_full_down_01_gene_counts, "one_full_down_01_gene_counts.csv")

length(one_full_down_01_gene_counts$Var1)
head(one_full_up_01_gene_counts)
head(one_full_down_01_gene_counts)


plot(ones$up_frags, ones$down_frags, pch=16, col="#808D8D4A")
