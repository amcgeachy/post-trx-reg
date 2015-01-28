getwd()
setwd("/Users/annamcgeachy/Google Drive/post-trx reg/datafiles/")

reading_frame_single_intron = function(inputfile, dataset_name){
  unique_orfs = read.table(inputfile, header=FALSE)
  colnames(unique_orfs) = c("frag_count", "chr_read", "start_read", "end_read", "strand_read",
                            "chr_gene", "start_cds", "end_cds", "gene_name", "bed_score", "strand_cds", 
                            "thick_start", "thick_end", "RGB", "exon_number", "exon_start", "exon_end", "overlap")
  
  unique_orfs[,"read_length"] = unique_orfs$end_read - unique_orfs$start_read
  head(unique_orfs)
  summary(unique_orfs$read_length)
  summary(unique_orfs$exon_number)
  
  no_introns = unique_orfs[unique_orfs$exon_number==1,]
  head(no_introns)
  
  no_introns_pos = no_introns[no_introns$strand_read=="+",]
  head(no_introns_pos)
  
  summary(no_introns_pos$strand_read)
  no_introns_pos[,"read_start_in_cds"] = no_introns_pos$start_read - no_introns_pos$start_cds
  head(no_introns_pos)
  no_introns_pos[,"read_start_aa_in_cds"] = no_introns_pos$read_start_in_cds / 3
  no_introns_pos[,"start_readframe_aa_in_cds"] = no_introns_pos$read_start_in_cds %% 3
  hist(no_introns_pos$start_readframe_aa_in_cds)
  
  no_introns_pos[,"read_end_in_cds"] = no_introns_pos$read_start_in_cds + no_introns_pos$read_length
  head(no_introns_pos)
  no_introns_pos[,"read_end_aa_in_cds"] = no_introns_pos$read_end_in_cds / 3
  no_introns_pos[,"end_readframe_aa_in_cds"] = no_introns_pos$read_end_in_cds %% 3
  no_introns_pos[,"readframe_tot"] = no_introns_pos$end_readframe_aa_in_cds + no_introns_pos$start_readframe_aa_in_cds
  hist(no_introns_pos$end_readframe_aa_in_cds, add=TRUE, col=NULL, border=2)
  
  hist(no_introns_pos$readframe_tot)
  
  head(no_introns_pos)
  
  zero_zero = no_introns_pos[which(no_introns_pos$start_readframe_aa_in_cds==0 & no_introns_pos$end_readframe_aa_in_cds==0),]
  head(zero_zero)
  write.csv(zero_zero, file=sprintf("%s zero_zero pos.csv", dataset_name))
  zero_one = no_introns_pos[which(no_introns_pos$start_readframe_aa_in_cds==0 & no_introns_pos$end_readframe_aa_in_cds==1),]
  zero_two = no_introns_pos[which(no_introns_pos$start_readframe_aa_in_cds==0 & no_introns_pos$end_readframe_aa_in_cds==2),]
  
  one_zero = no_introns_pos[which(no_introns_pos$start_readframe_aa_in_cds==1 & no_introns_pos$end_readframe_aa_in_cds==0),]
  head(one_zero)
  one_one = no_introns_pos[which(no_introns_pos$start_readframe_aa_in_cds==1 & no_introns_pos$end_readframe_aa_in_cds==1),]
  one_two = no_introns_pos[which(no_introns_pos$start_readframe_aa_in_cds==1 & no_introns_pos$end_readframe_aa_in_cds==2),]
  
  two_zero = no_introns_pos[which(no_introns_pos$start_readframe_aa_in_cds==2 & no_introns_pos$end_readframe_aa_in_cds==0),]
  head(two_zero)
  two_one = no_introns_pos[which(no_introns_pos$start_readframe_aa_in_cds==2 & no_introns_pos$end_readframe_aa_in_cds==1),]
  two_two = no_introns_pos[which(no_introns_pos$start_readframe_aa_in_cds==2 & no_introns_pos$end_readframe_aa_in_cds==2),]
  
  read_frame_dist = c(nrow(zero_zero), nrow(zero_one), nrow(zero_two),
           nrow(one_zero), nrow(one_one), nrow(one_two),
           nrow(two_zero), nrow(two_one), nrow(two_two))
    
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
  dev.off()
  
  new_list = list(no_introns, no_introns_pos)
}
a = NULL
b= NULL
test = reading_frame_single_intron("up_orf_unique.bed", "up")
up_no_introns = test[[1]]
up_no_introns_pos = test[[2]]

nrow(up_no_introns_pos) / nrow(up_no_introns)
nrow(up_no_introns_pos)
tail(test)
reading_frame_single_intron("down_orf_unique.bed", "down")
reading_frame_single_intron("no_recomb_orf_unique.bed", "no_recomb")
reading_frame_single_intron("post_recomb_orf_unique.bed", "post_recomb")

?type

typeof(unique_orfs)
head(zero_two)
plot(head(zero_two$frag_count, n=200))
summary(unique_orfs$exon_number)[[1]]

no_recomb_orfs = read.table("no_recomb_orf_unique.bed", header=FALSE)
colnames(no_recomb_orfs) = c("frag_count", "chr_read", "start_read", "end_read", "strand_read",
                          "chr_gene", "start_cds", "end_cds", "gene_name", "bed_score", "strand_cds", 
                          "thick_start", "thick_end", "RGB", "exon_number", "exon_start", "exon_end", "overlap")
head(no_recomb_orfs)
sum_no_recomb = summary(no_recomb_orfs$exon_number)
sum(sum_no_recomb)
sum_no_recomb[[1]] /sum(sum_no_recomb) #nongenic
1-sum_no_recomb[[1]]/sum(sum_no_recomb) #genic

post_recomb_orfs = read.table("post_recomb_orf_unique.bed", header=FALSE)
colnames(post_recomb_orfs) = c("frag_count", "chr_read", "start_read", "end_read", "strand_read",
                             "chr_gene", "start_cds", "end_cds", "gene_name", "bed_score", "strand_cds", 
                             "thick_start", "thick_end", "RGB", "exon_number", "exon_start", "exon_end", "overlap")
head(post_recomb_orfs)
summary(post_recomb_orfs$exon_number)
sum_post_recomb = summary(post_recomb_orfs$exon_number)
sum(sum_post_recomb)
sum_post_recomb[[1]] /sum(sum_post_recomb) #nongenic
1-sum_post_recomb[[1]]/sum(sum_post_recomb) #genic

sum_up = summary(unique_orfs$exon_number)
sum(sum_up)
sum_up[[1]] /sum(sum_up) #nongenic
1-sum_up[[1]]/sum(sum_up) #genic
pie(
  x=c(
    sum_up[[1]] /sum(sum_up),
    1-sum_up[[1]]/sum(sum_up)),
  labels=c("nongenic", "genic"), main = "up")


pie(
  x=c(
    sum_no_recomb[[1]] /sum(sum_no_recomb),
    1-sum_no_recomb[[1]]/sum(sum_no_recomb)),
  labels=c("nongenic", "genic"), main = "no recomb")

pie(
  x=c(
    sum_post_recomb[[1]] /sum(sum_post_recomb),
    1-sum_post_recomb[[1]]/sum(sum_post_recomb)),
  labels=c("nongenic", "genic"), main = "post recomb")

barplot(blah)
?text

2833925/148362
sum(136232 ,  94092, 2833925  ,148362 ,  76146  , 78148   ,99759 , 133444  ,102088)
2833925/sum(136232 ,  94092, 2833925  ,148362 ,  76146  , 78148   ,99759 , 133444  ,102088)

145224/sum(161245, 146153, 145224, 142945, 130254, 127734 ,148782 ,136660, 133554)
