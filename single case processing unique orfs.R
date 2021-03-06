getwd()
setwd("/Users/annamcgeachy/Google Drive/post-trx reg/datafiles_20140910_seq/")

unique_orfs = read.table("up_inside_orf_unique.bed", header=FALSE, stringsAsFactors = FALSE)
head(unique_orfs, n =10)
colnames(unique_orfs) = c("chr_read", "start_read", "end_read", "frag_count", "arbitrary_value", "strand_read",
                          "chr_gene", "start_cds", "end_cds", "gene_name", "bed_score", "strand_cds", 
                          "thick_start", "thick_end", "RGB", "exon_number", "exon_sizes", "exon_start", "overlap")

head(unique_orfs)
typeof(two_exon$exon_sizes)
as.numeric(strsplit(two_exon$exon_sizes[1], ",")[[1]][1])
sizes = strsplitAsListOfIntegerVectors(two_exon$exon_sizes)
starts = strsplitAsListOfIntegerVectors(two_exon$exon_start)
sizes[[1]]
starts[[1]]


thing = matrix(rep(NA, 236), nrow=236)
thing
for(i in 1:236){
thing[i] = test[i][[1]]
}
thing
?strsplit  
.libPaths()
library(GenomicRanges)
library(rtracklayer)
unique_orfs[,"read_length"] = unique_orfs$end_read - unique_orfs$start_read
head(unique_orfs)
two_exon = unique_orfs[which(unique_orfs$exon_number==2),]
head(two_exon)
??rtracklayer

number_multi_exon_genes = sum((summary(unique_orfs$exon_number)[3:length(summary(unique_orfs$exon_number))]))
length(summary(unique_orfs$exon_number))
number_single_exon_genes = summary(unique_orfs$exon_number)[2]
number_all_exon_genes = sum(summary(unique_orfs$exon_number)[2:length(summary(unique_orfs$exon_number))])
number_multi_exon_genes
number_single_exon_genes
number_all_exon_genes

percent_single_exon_genes = number_single_exon_genes/number_all_exon_genes
percent_single_exon_genes
percent_multi_exon_genes = number_multi_exon_genes/number_all_exon_genes
percent_multi_exon_genes
exon_percents = matrix(c(percent_single_exon_genes, percent_multi_exon_genes))
rownames(exon_percents) = c("single exon genes", "multi exon genes")
exon_percents

  #make list of single intron genes 
  no_introns = unique_orfs[unique_orfs$exon_number==1,]
  head(no_introns)
  
    #subset it by + strand and make metrics
    no_introns_pos = no_introns[no_introns$strand_read=="+",]
    head(no_introns_pos)
    
    no_introns_pos[,"read_start_in_cds"] = no_introns_pos$start_read - no_introns_pos$start_cds
    head(no_introns_pos)

    #subset it by - strand and make metrics
    no_introns_neg = no_introns[no_introns$strand_read=="-",]
    head(no_introns_neg)
    
    no_introns_neg[,"read_start_in_cds"] = no_introns_neg$end_cds - no_introns_neg$end_read 
    head(no_introns_neg)
    
    #combine - and + into one table
    no_introns_both = rbind(no_introns_neg, no_introns_pos)
    dim(no_introns_both)
    dim(no_introns_neg) + dim(no_introns_pos) #checks out

    no_introns_both[,"read_start_aa_in_cds"] = no_introns_both$read_start_in_cds / 3
    no_introns_both[,"start_readframe_aa_in_cds"] = no_introns_both$read_start_in_cds %% 3
    hist(no_introns_both$start_readframe_aa_in_cds)
  
    no_introns_both[,"read_end_in_cds"] = no_introns_both$read_start_in_cds + no_introns_both$read_length
    head(no_introns_both)
    no_introns_both[,"read_end_aa_in_cds"] = no_introns_both$read_end_in_cds / 3
    no_introns_both[,"end_readframe_aa_in_cds"] = no_introns_both$read_end_in_cds %% 3
    hist(no_introns_both$end_readframe_aa_in_cds, add=TRUE, col=NULL, border=2)
    
    no_introns_both[,"readframe_tot"] = no_introns_both$end_readframe_aa_in_cds + no_introns_both$start_readframe_aa_in_cds
    hist(no_introns_both$readframe_tot)

#subset the table by what reading frame the fragments are in
zero_zero = no_introns_both[which(no_introns_both$start_readframe_aa_in_cds==0 & no_introns_both$end_readframe_aa_in_cds==0),]
head(zero_zero)
zero_one = no_introns_both[which(no_introns_both$start_readframe_aa_in_cds==0 & no_introns_both$end_readframe_aa_in_cds==1),]
zero_two = no_introns_both[which(no_introns_both$start_readframe_aa_in_cds==0 & no_introns_both$end_readframe_aa_in_cds==2),]
head(zero_two)
nrow(zero_two)
one_zero = no_introns_both[which(no_introns_both$start_readframe_aa_in_cds==1 & no_introns_both$end_readframe_aa_in_cds==0),]
head(one_zero)
one_one = no_introns_both[which(no_introns_both$start_readframe_aa_in_cds==1 & no_introns_both$end_readframe_aa_in_cds==1),]
one_two = no_introns_both[which(no_introns_both$start_readframe_aa_in_cds==1 & no_introns_both$end_readframe_aa_in_cds==2),]

two_zero = no_introns_both[which(no_introns_both$start_readframe_aa_in_cds==2 & no_introns_both$end_readframe_aa_in_cds==0),]
head(two_zero)
two_one = no_introns_both[which(no_introns_both$start_readframe_aa_in_cds==2 & no_introns_both$end_readframe_aa_in_cds==1),]
two_two = no_introns_both[which(no_introns_both$start_readframe_aa_in_cds==2 & no_introns_both$end_readframe_aa_in_cds==2),]

read_frame_dist = c(nrow(zero_zero), nrow(zero_one), nrow(zero_two),
                    nrow(one_zero), nrow(one_one), nrow(one_two),
                    nrow(two_zero), nrow(two_one), nrow(two_two))

read_frame_dist_table = matrix(read_frame_dist, nrow =3)
rownames(read_frame_dist_table) = c("ends in zero", "ends in one", "ends in two")
colnames(read_frame_dist_table) = c("ends in zero", "ends in one", "ends in two")
read_frame_dist_table

pdf(sprintf("readingframe %s with numbers.pdf", "dataset_name"), useDingbats = FALSE)
barplot(read_frame_dist, xaxt = "n", main=sprintf("Single intron read frame distribution, %s", dataset_name))
labels=c("0,0", "0,1", "0,2",
         "1,0", "1,1", "1,2",
         "2,0", "2,1", "2,2")
axis(1, at=(1:9), labels=labels, cex=.05)

read_frame_dist_weighted = c(sum(zero_zero$frag_count), sum(zero_one$frag_count), sum(zero_two$frag_count),
                             sum(one_zero$frag_count), sum(one_one$frag_count), sum(one_two$frag_count),
                             sum(two_zero$frag_count), sum(two_one$frag_count), sum(two_two$frag_count))  
barplot(read_frame_dist_weighted, xaxt = "n", main=sprintf("Single intron weighted read frame distribution, %s", dataset_name))
labels=c("0,0", "0,1", "0,2",
         "1,0", "1,1", "1,2",
         "2,0", "2,1", "2,2")
axis(1, at=(1:9), labels=labels, cex=.05)

read_frame_dist_weighted_table = matrix(read_frame_dist_weighted, nrow =3)
rownames(read_frame_dist_weighted_table) = c("ends in zero", "ends in one", "ends in two")
colnames(read_frame_dist_weighted_table) = c("ends in zero", "ends in one", "ends in two")
read_frame_dist_weighted_table
dev.off()


#make counts for pie chart in terms of genic or not
intergenic = unique_orfs[which(unique_orfs$exon_number=="."),]
intergenic_weighted_count = sum(intergenic$frag_count)
genic = unique_orfs[which(unique_orfs$exon_number!="."),]
genic_weighted_count = sum(genic$frag_count)
all_genicity_weighted_count = sum(unique_orfs$frag_count)
intergenic_weighted_count + genic_weighted_count
all_genicity_weighted_count # checks out

#make them into tidy percents
percent_genic = round(genic_weighted_count/all_genicity_weighted_count, digits = 3)
percent_intergenic = round(intergenic_weighted_count/all_genicity_weighted_count, digits = 3)

#make the pie chart
pie(c(intergenic_weighted_count, genic_weighted_count), labels=c(sprintf("intergenic, %s", percent_intergenic), sprintf("genic, %s", percent_genic)), main="percent in gene, up")



