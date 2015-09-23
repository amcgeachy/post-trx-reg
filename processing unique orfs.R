getwd()
setwd("/Users/annamcgeachy/Google Drive/post-trx reg/datafiles_20140910_seq/")

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
  list(exon_percents = exon_percents, no_introns_both = no_introns_both, zero_two = zero_two, read_frame_dist_weighted = read_frame_dist_weighted)
}

#subset into lists of in frame fragments
  up = reading_frame_single_intron("up_inside_orf_unique.bed", "up")
  typeof(up)
  up$exon_percents
  head(up$no_introns_both)
  head(up$zero_two)
  nrow(up$zero_two)
  
  down = reading_frame_single_intron("down_inside_orf_unique.bed", "down")
  typeof(down)
  down$exon_percents
  head(down$no_introns_both)
  head(down$zero_two)
  nrow(down$zero_two)
  
  no_recomb = reading_frame_single_intron("no_recomb_inside_orf_unique.bed", "no_recomb")
  typeof(no_recomb)
  no_recomb$exon_percents
  head(no_recomb$no_introns_both)
  head(no_recomb$zero_two)
  nrow(no_recomb$zero_two)
  
  post_recomb = reading_frame_single_intron("post_recomb_inside_orf_unique.bed", "post_recomb")
  typeof(post_recomb)
  post_recomb$exon_percents
  head(post_recomb$no_introns_both)
  head(post_recomb$zero_two)
  nrow(post_recomb$zero_two)


read_frame_dist_weighted = c(sum(zero_zero$frag_count), sum(zero_one$frag_count), sum(zero_two$frag_count),
                             sum(one_zero$frag_count), sum(one_one$frag_count), sum(one_two$frag_count),
                             sum(two_zero$frag_count), sum(two_one$frag_count), sum(two_two$frag_count))  
barplot(read_frame_dist_weighted, xaxt = "n", main=sprintf("Single intron weighted read frame distribution, %s", dataset_name))
labels=c("0,0", "0,1", "0,2",
         "1,0", "1,1", "1,2",
         "2,0", "2,1", "2,2")
axis(1, at=(1:9), labels=labels, cex=.05)
barplot(up$read_frame_dist_weighted, xaxt = "n")
labels=c("0,0", "0,1", "0,2",
         "1,0", "1,1", "1,2",
         "2,0", "2,1", "2,2")
axis(1, at=(1:9), labels=labels, cex=.05)
no_recomb$read_frame_dist_weighted
barplot(no_recomb$read_frame_dist_weighted, xaxt = "n")
labels=c("0,0", "0,1", "0,2",
         "1,0", "1,1", "1,2",
         "2,0", "2,1", "2,2")
axis(1, at=(1:9), labels=labels, cex=.05)

up$read_frame_dist_weighted
#####

#Now combine the data into a single data frame. 

#make unique identifier for each fragment (not worrying now about intersecting duplicates)
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

  #and find the duplicates
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

  #take the unique fragments and put them into one table
    all_uniques = c(up_unique, down_unique, post_recomb_unique)
    all_uniques = as.data.frame(unique(all_uniques))
    colnames(all_uniques) = c("identifiers")
    head(all_uniques)
    nrow(all_uniques)

  #then add in the fragment counts per unique fragment
    all_uniques[,"up_frag"] = up$no_introns_both[match(all_uniques$identifiers, up$no_introns_both$unique), "frag_count"]
    all_uniques[,"down_frag"] = down$no_introns_both[match(all_uniques$identifiers, down$no_introns_both$unique), "frag_count"]
    all_uniques[,"no_recomb_frag"] = no_recomb$no_introns_both[match(all_uniques$identifiers, no_recomb$no_introns_both$unique), "frag_count"]
    all_uniques[,"post_recomb_frag"] = post_recomb$no_introns_both[match(all_uniques$identifiers, post_recomb$no_introns_both$unique), "frag_count"]
    head(all_uniques)
    typeof(all_uniques$up_frag)
    all_uniques[is.na(all_uniques)] = 0 #set NAs to 0
    write.csv(all_uniques, "all unique frag with counts per sample.csv")

  #find the highest count up and down (ends up being a glycolysis enzyme and PAT1 (decapping) respectively)
    all_uniques[which(all_uniques$up_frag==max(all_uniques$up_frag)),"identifiers"]
    all_uniques[which(all_uniques$down_frag==max(all_uniques$down_frag)),"identifiers"]

  #now that we have this table of condition a and condition b, lets plot to compare a v b
    plot_a_v_b_counts = function(data, condition1, condition2){
      temp = data[data[,condition1] + data[,condition2] >=1, c(condition1, condition2)]
      upper_lim = log(max(c(max(temp[,condition1]), max(temp[,condition2]))))
      
      plot(log(temp[,condition1]), log(temp[,condition2]), ylim=c(0,upper_lim+1), xlim=c(0,upper_lim+1),
           xlab=sprintf("log(%s)", condition1), ylab=sprintf("log(%s)", condition2), main=sprintf("%s v %s", condition1, condition2))
      abline(a=0,b=1)
      abline(a=-4.6,b=1)
      abline(a=4.6,b=1)
      abline(a=-2.3,b=1)
      abline(a=2.3,b=1)
      text(x=c(upper_lim-4.6,upper_lim-2.3,upper_lim), y=upper_lim+.5, labels=(c("100X", "10X", "0X")))
      }
    pdf("a v b logged plots.pdf", useDingbats = FALSE)
      plot_a_v_b_counts(all_uniques, "up_frag", "down_frag")
      plot_a_v_b_counts(all_uniques, "post_recomb_frag", "up_frag")
      plot_a_v_b_counts(all_uniques, "post_recomb_frag", "down_frag")
      plot_a_v_b_counts(all_uniques, "no_recomb_frag", "up_frag")
      plot_a_v_b_counts(all_uniques, "no_recomb_frag", "down_frag")
      plot_a_v_b_counts(all_uniques, "no_recomb_frag", "post_recomb_frag")
    dev.off()

a_v_b_counts = function(data, condition1, condition2){
  temp = data[data[,condition1] + data[,condition2] >=2^6, c("identifiers", condition1, condition2)]
  temp[,condition1] = log2(temp[,condition1])
  temp[,condition2] = log2(temp[,condition2])
  return(temp)
}

up_v_down = a_v_b_counts(all_uniques, "up_frag", "down_frag")
head(up_v_down)
dim(up_v_down)

up_v_down[,"enriched_in_up_10x"] = ifelse(up_v_down$up_frag > 2.3 + up_v_down$down_frag, up_v_down$up_frag, NA)
up_v_down[,"enriched_in_up_100x"] = ifelse(up_v_down$up_frag > 4.6 + up_v_down$down_frag, up_v_down$up_frag, NA)
up_v_down[,"enriched_in_down_10x"] = ifelse(up_v_down$down_frag > 2.3 + up_v_down$up_frag, up_v_down$down_frag, NA)
up_v_down[,"enriched_in_down_100x"] = ifelse(up_v_down$down_frag > 4.6 + up_v_down$up_frag, up_v_down$down_frag, NA)

head(up_v_down)
tail(up_v_down)

enriched_in_up = up_v_down[which(!is.na(up_v_down$enriched_in_up_10x) & !is.na(up_v_down$enriched_in_up_100x)),]
enriched_in_down = up_v_down[which(!is.na(up_v_down$enriched_in_down_10x) & !is.na(up_v_down$enriched_in_down_100x)),]
dim(enriched_in_down)
dim(enriched_in_up)
background = up_v_down[which(is.na(up_v_down$enriched_in_up_10x) & is.na(up_v_down$enriched_in_up_100x) 
                             & is.na(up_v_down$enriched_in_down_10x) & is.na(up_v_down$enriched_in_down_100x)),]

#plot boxplots
pdf("boxplots of frag count by enrichment.pdf", useDingbats = FALSE)
boxplot(background$up_frag, background$down_frag, xaxt='n', ylab='log(frag_count)', main="background")
axis(1, at=1:2, labels=c("up_frag", "down_frag"))
boxplot(enriched_in_up$up_frag, enriched_in_up$down_frag, xaxt='n', ylab='log(frag_count)', main="enriched in up")
axis(1, at=1:2, labels=c("up_frag", "down_frag"))
boxplot(enriched_in_down$up_frag, enriched_in_down$down_frag, xaxt='n', ylab='log(frag_count)', main="enriched in down")
axis(1, at=1:2, labels=c("up_frag", "down_frag"))
dev.off()

#make ratio of fragment counts
enriched_in_up$ratio = enriched_in_up$up_frag/enriched_in_up$down_frag
enriched_in_down$ratio = enriched_in_down$down_frag/enriched_in_down$up_frag

#pull the highest count 
#picks those with the highest ratio, that being inf (because the alt strain has 0 counts)
enriched_in_up[which(enriched_in_up$up_frag==max(enriched_in_up[which(enriched_in_up$ratio==max(enriched_in_up$ratio)),"up_frag"])),]
enriched_in_down[which(enriched_in_down$down_frag==max(enriched_in_down[which(enriched_in_down$ratio==max(enriched_in_down$ratio)),"down_frag"])),]
#picks the highest count, period (alt strain count not = 0 )
enriched_in_down[which(enriched_in_down$down_frag==max(enriched_in_down$down_frag)),]
enriched_in_up[which(enriched_in_up$up_frag==max(enriched_in_up$up_frag)),]

#read in data file with xref potential (has gene names and descriptions)
xref = read.delim("SGD_features.tab", header=FALSE, quote="")

  #up
    #make a table with all of the fragment information (from no_intron_both) with enrichment information (10X, 100X, and ratio)
    enriched_up_more_info = up$no_introns_both[match(enriched_in_up$identifiers, up$no_introns_both$unique),]
    enriched_up_more_info$enriched_in_up_10x = enriched_in_up[match(enriched_up_more_info$unique, enriched_in_up$identifiers),"enriched_in_up_10x"]
    enriched_up_more_info$enriched_in_up_100x = enriched_in_up[match(enriched_up_more_info$unique, enriched_in_up$identifiers),"enriched_in_up_100x"]
    enriched_up_more_info$ratio = enriched_in_up[match(enriched_up_more_info$unique, enriched_in_up$identifiers),"ratio"]
    
head(enriched_up_more_info)
    #now cross reference it with xref to add gene description
    enriched_up_more_info$gene_useful = xref[match(enriched_up_more_info$gene_name, xref$V4),"V5"]
    enriched_up_more_info$gene_desc = xref[match(enriched_up_more_info$gene_name, xref$V4),"V16"]

    #then reorder it based on ratio of enrichment
    enriched_up_more_info = enriched_up_more_info[order(enriched_up_more_info$ratio),]
    head(enriched_up_more_info)
    write.csv(enriched_up_more_info, "up enriched fragments 20150218.csv")

    #now make a frag count per gene table
    frag_per_gene_up = as.data.frame(table(enriched_up_more_info$gene_name))
    frag_per_gene_up = frag_per_gene_up[which(frag_per_gene_up$Freq!=0),]
    frag_per_gene_up = frag_per_gene_up[order(frag_per_gene_up$Freq, decreasing = TRUE),]
    frag_per_gene_up$gene_useful = xref[match(frag_per_gene_up$Var1, xref$V4), "V5"]
    frag_per_gene_up$desc = xref[match(frag_per_gene_up$Var1, xref$V4), "V16"]
    head(frag_per_gene_up)
    write.csv(frag_per_gene_up, "up frag table count.csv")

  #down
    #make a table with all of the fragment information (from no_intron_both) with enrichment information (10X, 100X, and ratio)
    enriched_down_more_info = down$no_introns_both[match(enriched_in_down$identifiers, down$no_introns_both$unique),]
    enriched_down_more_info$enriched_in_down_10x = enriched_in_down[match(enriched_down_more_info$unique, enriched_in_down$identifiers),"enriched_in_down_10x"]
    enriched_down_more_info$enriched_in_down_100x = enriched_in_down[match(enriched_down_more_info$unique, enriched_in_down$identifiers),"enriched_in_down_100x"]
    enriched_down_more_info$ratio = enriched_in_down[match(enriched_down_more_info$unique, enriched_in_down$identifiers),"ratio"]
    
    #now cross reference it with xref to add gene description
    enriched_down_more_info$gene_useful = xref[match(enriched_down_more_info$gene_name, xref$V4),"V5"]
    enriched_down_more_info$gene_desc = xref[match(enriched_down_more_info$gene_name, xref$V4),"V16"]
    
    #then reorder it based on ratio of enrichment
    enriched_down_more_info = enriched_down_more_info[order(enriched_down_more_info$ratio),]
    head(enriched_down_more_info)
    write.csv(enriched_down_more_info, "down enriched fragments 20150218.csv")
    
  #now make a frag count per gene table
    frag_per_gene_down = as.data.frame(table(enriched_down_more_info$gene_name))
    frag_per_gene_down = frag_per_gene_down[which(frag_per_gene_down$Freq!=0),]
    frag_per_gene_down = frag_per_gene_down[order(frag_per_gene_down$Freq, decreasing = TRUE),]
    frag_per_gene_down$gene_useful = xref[match(frag_per_gene_down$Var1, xref$V4), "V5"]
    frag_per_gene_down$desc = xref[match(frag_per_gene_down$Var1, xref$V4), "V16"]
    head(frag_per_gene_down)
    write.csv(frag_per_gene_down, "down frag table count.csv")

#subset the above by those that are in frame (0,2)
up$zero_two[,"unique"] = paste(up$zero_two$chr_read, up$zero_two$start_read, 
                               up$zero_two$end_read, up$zero_two$strand_read, sep="_")
down$zero_two[,"unique"] = paste(down$zero_two$chr_read, down$zero_two$start_read, 
                                 down$zero_two$end_read, down$zero_two$strand_read, sep="_")

up_zero_two_list  = match(up$zero_two$unique, enriched_up_more_info$unique)
up_zero_two_list_noNA = up_zero_two_list[!is.na(up_zero_two_list)]
enriched_up_zero_two = enriched_up_more_info[up_zero_two_list_noNA,]
enriched_up_zero_two = enriched_up_zero_two[order(enriched_up_zero_two$ratio),]
head(enriched_up_zero_two)
dim(enriched_up_more_info)
dim(enriched_up_zero_two)
write.csv(enriched_up_zero_two, "in frame up enriched.csv")

#make fragment table for this as well
up_zero_two_frag = as.data.frame(table(enriched_up_zero_two$gene_name))
up_zero_two_frag = up_zero_two_frag[which(up_zero_two_frag$Freq!=0),]
up_zero_two_frag = up_zero_two_frag[order(up_zero_two_frag$Freq, decreasing = TRUE),]
up_zero_two_frag$gene_useful = xref[match(up_zero_two_frag$Var1, xref$V4), "V5"]
up_zero_two_frag$desc = xref[match(up_zero_two_frag$Var1, xref$V4), "V16"]
write.csv(up_zero_two_frag, "in frame up fragments table.csv")

down_zero_two_list  = match(down$zero_two$unique, enriched_down_more_info$unique)
down_zero_two_list_noNA = down_zero_two_list[!is.na(down_zero_two_list)]
enriched_down_zero_two = enriched_down_more_info[down_zero_two_list_noNA,]
enriched_down_zero_two = enriched_down_zero_two[order(enriched_down_zero_two$ratio),]
head(enriched_down_zero_two)
dim(enriched_down_more_info)
dim(enriched_down_zero_two)
write.csv(enriched_down_zero_two, "in frame down enriched.csv")

#make fragment table for this as well
down_zero_two_frag = as.data.frame(table(enriched_down_zero_two$gene_name))
down_zero_two_frag = down_zero_two_frag[which(down_zero_two_frag$Freq!=0),]
down_zero_two_frag = down_zero_two_frag[order(down_zero_two_frag$Freq, decreasing = TRUE),]
down_zero_two_frag$gene_useful = xref[match(down_zero_two_frag$Var1, xref$V4), "V5"]
down_zero_two_frag$desc = xref[match(down_zero_two_frag$Var1, xref$V4), "V16"]
write.csv(down_zero_two_frag, "in frame down fragments table.csv")

thing = unique(post_recomb$no_introns_both$gene_name)
thing = thing[which(thing)]
write.csv(thing, "all genes.csv")

head(all_uniques)
summary(all_uniques$up_frag)
hist(as.data.frame(table(all_uniques[which(all_uniques$up_frag>3 & all_uniques$up_frag<100),"up_frag"])))

#now look at these (up and down) versus post recombination
head(all_uniques)
up_v_post = a_v_b_counts(all_uniques, "up_frag", "post_recomb_frag")
down_v_post = a_v_b_counts(all_uniques, "down_frag", "post_recomb_frag")

#A v B plot for up v down using opacity to show density
all_uniques_copy = all_uniques
head(all_uniques_copy)
all_uniques_copy$up_frag = ifelse(all_uniques_copy$up_frag==0, .5, all_uniques_copy$up_frag)
all_uniques_copy$down_frag = ifelse(all_uniques_copy$down_frag==0, .5, all_uniques_copy$down_frag)

pdf("up v down opacity with genes of interest.pdf", useDingbats = FALSE)
plot(log2(all_uniques_copy$down_frag), log2(all_uniques_copy$up_frag), pch=16, col="#808D8D4A",
     xlab="log2(unique fragment counts in down library)",
     ylab="log2(unique fragment counts in up library)",
     main="unique fragments counts in up v down")
abline(a=0,b=1)
abline(a=-4.6,b=1)
abline(a=4.6,b=1)
abline(a=-2.3,b=1)
abline(a=2.3,b=1)

#pull unique identifiers
PAT1_unique = enriched_down_more_info[grep("PAT1", enriched_down_more_info$gene_useful),"unique"]
YMR295C_unique = enriched_down_more_info[grep("YMR295C", enriched_down_more_info$gene_name),"unique"]

CDC19_unique = enriched_up_more_info[grep("CDC19", enriched_up_more_info$gene_useful),"unique"]
YPR204W_unique = enriched_up_more_info[grep("YPR204W", enriched_up_more_info$gene_name),"unique"]

#and the coordinations
PAT1_coord = all_uniques_copy[match(PAT1_unique, all_uniques_copy$identifiers),]
YMR295C_coord = all_uniques_copy[match(YMR295C_unique, all_uniques_copy$identifiers),]

CDC19_coord = all_uniques_copy[match(CDC19_unique, all_uniques_copy$identifiers),]
YPR204W_coord = all_uniques_copy[match(YPR204W_unique, all_uniques_copy$identifiers),]

PAT1_coord[order(PAT1_coord$identifiers),]

#then plot
points(log2(PAT1_coord$down_frag), log2(PAT1_coord$up_frag),
       pch=16, col="red")
points(log2(YMR295C_coord$down_frag), log2(YMR295C_coord$up_frag),
       pch=16, col="orange")

points(log2(CDC19_coord$down_frag), log2(CDC19_coord$up_frag),
       pch=16, col="green")
points(log2(YPR204W_coord$down_frag), log2(YPR204W_coord$up_frag),
       pch=16, col="blue")
dev.off()

# up$no_introns_both[,"most_unique"] = paste(up$no_introns_both$chr_read, up$no_introns_both$start_read, 
#                                       up$no_introns_both$end_read, up$no_introns_both$strand_read, 
#                                       up$no_introns_both$gene_name, sep="_")
# down$no_introns_both[,"most_unique"] = paste(down$no_introns_both$chr_read, down$no_introns_both$start_read, 
#                                         down$no_introns_both$end_read, down$no_introns_both$strand_read, 
#                                         down$no_introns_both$gene_name, sep="_")
# no_recomb$no_introns_both[,"most_unique"] = paste(no_recomb$no_introns_both$chr_read, no_recomb$no_introns_both$start_read,
#                                              no_recomb$no_introns_both$end_read, no_recomb$no_introns_both$strand_read, 
#                                              no_recomb$no_introns_both$gene_name, sep="_")
# post_recomb$no_introns_both[,"most_unique"] = paste(post_recomb$no_introns_both$chr_read, post_recomb$no_introns_both$start_read,
#                                                post_recomb$no_introns_both$end_read, post_recomb$no_introns_both$strand_read, 
#                                                post_recomb$no_introns_both$gene_name, sep="_")
# head(up$no_introns_both)

# up_most_unique = up$no_introns_both$most_unique
# down_most_unique = down$no_introns_both$most_unique
# no_recomb_most_unique = no_recomb$no_introns_both$most_unique
# post_recomb_most_unique = post_recomb$no_introns_both$most_unique
# 
# 
# all_uniquers = c(up_most_unique, down_most_unique, no_recomb_most_unique, post_recomb_most_unique)