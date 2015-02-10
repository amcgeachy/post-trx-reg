getwd()
setwd("/Users/annamcgeachy/Google Drive/post-trx reg/datafiles_new/")

reading_frame_single_intron = function(inputfile, dataset_name){
  #import data
  unique_orfs = read.table(inputfile, header=FALSE)
  colnames(unique_orfs) = c("chr_read", "start_read", "end_read", "frag_count", "arbitrary_value", "strand_read",
                            "chr_gene", "start_cds", "end_cds", "gene_name", "bed_score", "strand_cds", 
                            "thick_start", "thick_end", "RGB", "exon_number", "exon_start", "exon_end", "overlap")
  
  head(unique_orfs)
  
  unique_orfs[,"read_length"] = unique_orfs$end_read - unique_orfs$start_read
  head(unique_orfs)
  
  #note % with single exon, since we're only doing 1 exon currently
  number_multi_exon_genes = sum((summary(unique_orfs$exon_number)[3:length(summary(unique_orfs$exon_number))]))
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
    pdf(sprintf("pie chart of genecity %s.pdf", dataset_name), useDingbats = FALSE)
    pie(c(intergenic_weighted_count, genic_weighted_count), 
        labels=c(sprintf("intergenic, %s", percent_intergenic), sprintf("genic, %s", percent_genic)),
        main=sprintf("percent in gene, %s", dataset_name))
    dev.off()
  
  #make a subset of data that we'll work with for now (single intron) and add metrics
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
    
    pdf(sprintf("start end reading frame for %s.pdf", dataset_name), useDingbats = FALSE)
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
    dev.off()
    
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
    
    #make a table of unweighted counts for each reading frame combo
      read_frame_dist = c(nrow(zero_zero), nrow(zero_one), nrow(zero_two),
                          nrow(one_zero), nrow(one_one), nrow(one_two),
                          nrow(two_zero), nrow(two_one), nrow(two_two))
      
      read_frame_dist_table = matrix(read_frame_dist, nrow =3)
      rownames(read_frame_dist_table) = c("ends in zero", "ends in one", "ends in two")
      colnames(read_frame_dist_table) = c("ends in zero", "ends in one", "ends in two")
      read_frame_dist_table
      
      #plot it
      pdf(sprintf("readingframe %s with numbers.pdf", dataset_name), useDingbats = FALSE)
      barplot(read_frame_dist, xaxt = "n", main=sprintf("Single intron read frame distribution, %s", dataset_name))
      labels=c("0,0", "0,1", "0,2",
               "1,0", "1,1", "1,2",
               "2,0", "2,1", "2,2")
      axis(1, at=(1:9), labels=labels, cex=.05)
    
    #now do it again, but weighted
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
  
  #return information into useful objects to use outside of the function
  list(exon_percents = exon_percents, no_introns_both = no_introns_both, zero_two = zero_two)
}

#subset into lists of in frame fragments
  up = reading_frame_single_intron("up_inside_orf_unique.bed", "up")
  typeof(up)
  up$exon_percents
  head(up$no_introns_both)
  head(up$zero_two)
  
  down = reading_frame_single_intron("down_inside_orf_unique.bed", "down")
  typeof(down)
  down$exon_percents
  head(down$no_introns_both)
  head(down$zero_two)
  
  no_recomb = reading_frame_single_intron("no_recomb_inside_orf_unique.bed", "no_recomb")
  typeof(no_recomb)
  no_recomb$exon_percents
  head(no_recomb$no_introns_both)
  head(no_recomb$zero_two)
  
  post_recomb = reading_frame_single_intron("post_recomb_inside_orf_unique.bed", "post_recomb")
  typeof(post_recomb)
  post_recomb$exon_percents
  head(post_recomb$no_introns_both)
  head(post_recomb$zero_two)

#####

#Now combine the data into a single data frame. 

#make unique identifier for each fragment
head(up$no_introns_both)
up$no_introns_both[,"unique"] = paste(up$no_introns_both$chr_read, up$no_introns_both$start_read, 
                                      up$no_introns_both$end_read, up$no_introns_both$strand_read, sep="_")
down$no_introns_both[,"unique"] = paste(down$no_introns_both$chr_read, down$no_introns_both$start_read, 
                                        down$no_introns_both$end_read, down$no_introns_both$strand_read, sep="_")
no_recomb$no_introns_both[,"unique"] = paste(no_recomb$no_introns_both$chr_read, no_recomb$no_introns_both$start_read,
                                             no_recomb$no_introns_both$end_read, no_recomb$no_introns_both$strand_read, sep="_")
post_recomb$no_introns_both[,"unique"] = paste(post_recomb$no_introns_both$chr_read, post_recomb$no_introns_both$start_read,
                                               post_recomb$no_introns_both$end_read, post_recomb$no_introns_both$strand_read, sep="_")
head(up$no_introns_both)

up_unique = up$no_introns_both$unique
up_dup = up_unique[duplicated(up_unique)]
head(up_dup)
length(up_dup)

no_recomb_unique = no_recomb$no_introns_both$unique
no_recomb_dup = no_recomb_unique[duplicated(no_recomb_unique)]
head(no_recomb_dup)
length(no_recomb_dup)

post_recomb_unique = post_recomb$no_introns_both$unique
post_recomb_dup = post_recomb_unique[duplicated(post_recomb_unique)]
head(post_recomb_dup)
length(post_recomb_dup)

down_unique = down$no_introns_both$unique
down_dup = down_unique[duplicated(down_unique)]
head(down_dup)
length(down_dup)

up_duplicates =lapply(up_dup, function(x) up$no_introns_both[which(up$no_introns_both$unique==x),])
no_recomb_duplicates = lapply(no_recomb_dup, function(x) no_recomb$no_introns_both[which(no_recomb$no_introns_both$unique==x),])
post_recomb_duplicates = lapply(post_recomb_dup, function(x) post_recomb$no_introns_both[which(post_recomb$no_introns_both$unique==x),])
down_duplicates =lapply(down_dup, function(x) down$no_introns_both[which(down$no_introns_both$unique==x),])

up$no_introns_both[,"most_unique"] = paste(up$no_introns_both$chr_read, up$no_introns_both$start_read, 
                                      up$no_introns_both$end_read, up$no_introns_both$strand_read, 
                                      up$no_introns_both$gene_name, sep="_")
down$no_introns_both[,"most_unique"] = paste(down$no_introns_both$chr_read, down$no_introns_both$start_read, 
                                        down$no_introns_both$end_read, down$no_introns_both$strand_read, 
                                        down$no_introns_both$gene_name, sep="_")
no_recomb$no_introns_both[,"most_unique"] = paste(no_recomb$no_introns_both$chr_read, no_recomb$no_introns_both$start_read,
                                             no_recomb$no_introns_both$end_read, no_recomb$no_introns_both$strand_read, 
                                             no_recomb$no_introns_both$gene_name, sep="_")
post_recomb$no_introns_both[,"most_unique"] = paste(post_recomb$no_introns_both$chr_read, post_recomb$no_introns_both$start_read,
                                               post_recomb$no_introns_both$end_read, post_recomb$no_introns_both$strand_read, 
                                               post_recomb$no_introns_both$gene_name, sep="_")
head(up$no_introns_both)

up_most_unique = up$no_introns_both$most_unique
down_most_unique = down$no_introns_both$most_unique
no_recomb_most_unique = no_recomb$no_introns_both$most_unique
post_recomb_most_unique = post_recomb$no_introns_both$most_unique


all_uniquers = c(up_most_unique, down_most_unique, no_recomb_most_unique, post_recomb_most_unique)

#length(up_most_unique)+length(down_most_unique)+length(no_recomb_most_unique)+length(post_recomb_most_unique)
#length(all_uniquers)

unique_all_uniquers = unique(all_uniquers)
length(unique_all_uniquers)
test = as.data.frame(unique_all_uniquers)
head(test)
test[,"up"] = 0 
test[,"down"] = 0
test[,"no_recomb"] = 0 
test[,"post_recomb"] = 0 

ifelse(up$no_introns_both$most_unique==test$unique_all_uniquers, up$no_introns_both$frag_count, 0)
head(test)
up$no_introns_both[which(up$no_introns_both$most_unique=="chrI_100399_100568_-_YAL025C"),]

?ifelse
