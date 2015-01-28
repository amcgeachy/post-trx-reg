getwd()
setwd("/Users/annamcgeachy/Google Drive/post-trx reg/datafiles/")

unique_orfs = read.table("up_orf_unique.bed", header=FALSE)
colnames(unique_orfs) = c("frag_count", "chr_read", "start_read", "end_read", "strand_read",
                          "chr_gene", "start_cds", "end_cds", "gene_name", "bed_score", "strand_cds", 
                          "thick_start", "thick_end", "RGB", "exon_number", "exon_start", "exon_end", "overlap")

head(unique_orfs)

unique_orfs[,"read_length"] = unique_orfs$end_read - unique_orfs$start_read
head(unique_orfs)

number_multi_exon_genes = sum((summary(unique_orfs$exon_number)[3:8]))
number_single_exon_genes = summary(unique_orfs$exon_number)[2]
number_all_exon_genes = sum(summary(unique_orfs$exon_number)[2:8])
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
    dim(no_introns_neg) + dim(no_introns_pos)
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

pdf(sprintf("readingframe %s with numbers.pdf", dataset_name), useDingbats = FALSE)
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

new_list = list(no_introns, no_introns_pos)