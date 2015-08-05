getwd()
setwd("/Users/annamcgeachy/Google Drive/post-trx reg/datafiles_20140910_seq/")

unique_orfs = read.table("up_inside_orf_unique.bed", header=FALSE, stringsAsFactors = FALSE)
head(unique_orfs)
colnames(unique_orfs) = c("chr_read", "start_read", "end_read", "frag_count", "arbitrary_value", "strand_read",
                          "chr_gene", "start_cds", "end_cds", "gene_name", "bed_score", "strand_cds", 
                          "thick_start", "thick_end", "RGB", "exon_number", "exon_size", "exon_start", "overlap")

colnames(unique_orfs)

unique_orfs[,"read_length"] = unique_orfs$end_read - unique_orfs$start_read
head(unique_orfs)

#another method of pulling the non-exon numbers, a bit messier than below
# non_period_or_1_exon_numbers = as.numeric(row.names(table(unique_orfs$exon_number)[3:(length(table(unique_orfs$exon_number)))]))
# thing = unlist(sapply(non_period_or_1_exon_numbers, grep, unique_orfs$exon_number))
# multi_exon_unique_orfs = unique_orfs[thing,]
# table(multi_exon_unique_orfs$exon_number)

#another another method, without using dplyr
# zero_exons = grep(".", unique_orfs$exon_number, fixed = TRUE)
# one_exon = grep(1, unique_orfs$exon_number, fixed = TRUE)
# zero_and_one_exons = c(zero_exons, one_exon)
# all_rows = 1:nrow(unique_orfs)
# multi_exon_rows = setdiff(all_rows, zero_and_one_exons)
# 
# multi_exons_unique_orfs = unique_orfs[multi_exon_rows,]
# head(multi_exons_unique_orfs)
# table(multi_exons_unique_orfs$exon_number)

library("dplyr")
# library("nycflights13")
# colnames(table(unique_orfs$exon_number))
# 
# multi_exons_unique_orfs = filter(unique_orfs, exon_number!=".", exon_number!="1")
# table(multi_exons_unique_orfs$exon_number)[[1]]
# 
# baby = head(multi_exons_unique_orfs)
# 
# library("tidyr")
# 
# baby = separate(baby, exon_start, into=c("start1", "start2", "start3"))
# baby = separate(baby, exon_end, into=c("end1", "end2", "end3"))
# baby
# 
# baby2 = head(multi_exons_unique_orfs)
# baby2
# separate(baby2, exon_end)

##^----this works fine for things with just 2 exons, but what about more? 
#try a mixed population of intron numbers

#proving to myself that there are at max 8 exons in a gene
# saccer3 = read.table("/Users/annamcgeachy/saccer3.bed")
# head(saccer3)
# table(saccer3$V10)
# saccer3[which(saccer3$V10==8),]

# mixed_case = multi_exons_unique_orfs[c(13,51,20,49,243),]
# mixed_case
# 
# tmp = strsplit(mixed_case$exon_start, ",")
# 
# empty=matrix(data=NA, nrow=5, ncol=8)
# empty
# colnames(empty) = c(1:8)
# 
# for (i in 1:5){
#   for (j in 1:8){
#     empty[i,j]=as.numeric(tmp[[i]][j])
#     colnames(empty)[j]=c(sprintf("start_%s", j))
#   }
# }
# 
# empty
# 
# tmp_end = strsplit(mixed_case$exon_end, ",")
# 
# empty_end=matrix(data=NA, nrow=5, ncol=8)
# empty_end
# colnames(empty_end) = c(1:8)
# 
# for (i in 1:5){
#   for (j in 1:8){
#     empty_end[i,j]=as.numeric(tmp[[i]][j])
#     colnames(empty_end)[j]=c(sprintf("end_%s", j))
#   }
# }
# 
# empty_end
# 
# cbind(mixed_case, empty, empty_end)

#it works! that's SO EXCITING!

#try it for the big data set
  head(unique_orfs)
  nrow(unique_orfs)

#first filter out only the things that are genic
  genic = unique_orfs[which(unique_orfs$exon_number!="."),]
  table(genic$exon_number)
  
  head(genic)
  
  
  neg = filter(genic, strand_read=="-")
  pos = filter(genic, strand_read=="+")
  
  nrow(neg) + nrow(pos)
  nrow(genic)

#split exons from column with commas into distinct columns
  start_split_genic = strsplit(genic$exon_start, ",")
  
  empty_start_genic=matrix(data=NA, nrow=nrow(genic), ncol=8)
  colnames(empty_start_genic) = c(1:8)
  
  for (i in 1:nrow(genic)){
    for (j in 1:8){
      empty_start_genic[i,j]=as.numeric(start_split_genic[[i]][j])
      colnames(empty_start_genic)[j]=c(sprintf("start_%s", j))
    }
  }
  
  exon_size_split_genic = strsplit(genic$exon_size, ",")
  
  empty_exon_size_genic=matrix(data=NA, nrow=nrow(genic), ncol=8)
  colnames(empty_exon_size_genic) = c(1:8)
  
  for (i in 1:nrow(genic)){
    for (j in 1:8){
      empty_exon_size_genic[i,j]=as.numeric(exon_size_split_genic[[i]][j])
      colnames(empty_exon_size_genic)[j]=c(sprintf("exon_size_%s", j))
    }
  }
  
  proc.time()
  
  genic_exons = cbind(genic, empty_start_genic, empty_exon_size_genic)

#define absolute genomic coordinates of exons
  starts=NULL
  ends=NULL

  for (i in 1:8){
    starts[i] = paste("start_exon_", i, sep="")
    ends[i] = paste("end_exon_", i, sep="")
  }
 

#for positive strand genes
pos_exons = filter(genic_exons, strand_read=="+")

  for (i in 1:8){
    pos_exons[,starts[i]] = pos_exons$start_cds + pos_exons[,sprintf("start_%s", i)]
    pos_exons[,ends[i]] = pos_exons$start_cds + pos_exons[,sprintf("start_%s", i)] + pos_exons[,sprintf("exon_size_%s", i)]
  }

#for negative strand genes
neg_exons = filter(genic_exons, strand_read=="-")

  for (i in 1:8){
    neg_exons[,starts[i]] = neg_exons$end_cds - neg_exons[,sprintf("start_%s", i)]
    neg_exons[,ends[i]] = neg_exons$end_cds - neg_exons[,sprintf("start_%s", i)] - neg_exons[,sprintf("exon_size_%s", i)]
  }

#now we need to determine which (if any) exon the fragment occurs in

  # 
  # if (two_exon_pos$start_read > two_exon_pos$start_exon_1 & two_exon_pos$end_read < two_exon_pos$end_exon_1) { 
  #   two_exon_pos$occurs_in_exon = "start_exon_1"
  #   } else if (!is.na(two_exon_pos$start_exon_2) & (two_exon_pos$start_read > two_exon_pos$start_exon_2 & two_exon_pos$end_read < two_exon_pos$end_exon_2)){
  #     two_exon_pos$occurs_in_exon = "start_exon_2"
  #   } else if (!is.na(two_exon_pos$start_exon_3) & (two_exon_pos$start_read > two_exon_pos$start_exon_3 & two_exon_pos$end_read < two_exon_pos$end_exon_3)){
  #     two_exon_pos$occurs_in_exon = "start_exon_3"
  #   } else if (!is.na(two_exon_pos$start_exon_4) & (two_exon_pos$start_read > two_exon_pos$start_exon_4 & two_exon_pos$end_read < two_exon_pos$end_exon_4)){
  #     two_exon_pos$occurs_in_exon = "start_exon_4"
  #   } else if (!is.na(two_exon_pos$start_exon_5) & (two_exon_pos$start_read > two_exon_pos$start_exon_5 & two_exon_pos$end_read < two_exon_pos$end_exon_5)){
  #     two_exon_pos$occurs_in_exon = "start_exon_5"
  #   } else if (!is.na(two_exon_pos$start_exon_6) & (two_exon_pos$start_read > two_exon_pos$start_exon_6 & two_exon_pos$end_read < two_exon_pos$end_exon_6)){
  #     two_exon_pos$occurs_in_exon = "start_exon_6"
  #   } else if (!is.na(two_exon_pos$start_exon_7) & (two_exon_pos$start_read > two_exon_pos$start_exon_7 & two_exon_pos$end_read < two_exon_pos$end_exon_7)){
  #     two_exon_pos$occurs_in_exon = "start_exon_7"
  #   } else if (!is.na(two_exon_pos$start_exon_8) & (two_exon_pos$start_read > two_exon_pos$start_exon_8 & two_exon_pos$end_read < two_exon_pos$end_exon_8)){
  #     two_exon_pos$occurs_in_exon = "start_exon_8"
  #   } else {
  #     two_exon_pos$occurs_in_exon = "intronic"
  # }
  # 
  # two_exon_pos
  
  #function to check if exonic and if so which exon for positive strand
  positive_exon = function(generic_positive){
    if (generic_positive$start_read > generic_positive$start_exon_1 & generic_positive$end_read < generic_positive$end_exon_1) { 
      generic_positive$occurs_in_exon = "start_exon_1"
    } else if (!is.na(generic_positive$start_exon_2) & (generic_positive$start_read > generic_positive$start_exon_2 & generic_positive$end_read < generic_positive$end_exon_2)){
      generic_positive$occurs_in_exon = "start_exon_2"
    } else if (!is.na(generic_positive$start_exon_3) & (generic_positive$start_read > generic_positive$start_exon_3 & generic_positive$end_read < generic_positive$end_exon_3)){
      generic_positive$occurs_in_exon = "start_exon_3"
    } else if (!is.na(generic_positive$start_exon_4) & (generic_positive$start_read > generic_positive$start_exon_4 & generic_positive$end_read < generic_positive$end_exon_4)){
      generic_positive$occurs_in_exon = "start_exon_4"
    } else if (!is.na(generic_positive$start_exon_5) & (generic_positive$start_read > generic_positive$start_exon_5 & generic_positive$end_read < generic_positive$end_exon_5)){
      generic_positive$occurs_in_exon = "start_exon_5"
    } else if (!is.na(generic_positive$start_exon_6) & (generic_positive$start_read > generic_positive$start_exon_6 & generic_positive$end_read < generic_positive$end_exon_6)){
      generic_positive$occurs_in_exon = "start_exon_6"
    } else if (!is.na(generic_positive$start_exon_7) & (generic_positive$start_read > generic_positive$start_exon_7 & generic_positive$end_read < generic_positive$end_exon_7)){
      generic_positive$occurs_in_exon = "start_exon_7"
    } else if (!is.na(generic_positive$start_exon_8) & (generic_positive$start_read > generic_positive$start_exon_8 & generic_positive$end_read < generic_positive$end_exon_8)){
      generic_positive$occurs_in_exon = "start_exon_8"
    } else {
      generic_positive$occurs_in_exon = "intronic"
    }}
TRUE && FALSE 
#negative strand reads


#testing it on a small population
# two_exon_neg = head(filter(neg_exons, strand_cds=="-", exon_number==2), n=1)
# two_exon_neg
# 
# (end) 1 .... 3 ... 5... 9... 15 (start)
# to be in read, it has to be less than the start and greater than the end 
# eg a read from (start)3-(end)9 would be in the above exon
# so start_read>end_exon, end_read<start_exon
# 
# two_exon_neg$occurs_in_exon = NA
#   if (!is.na(two_exon_neg$start_exon_1) & two_exon_neg$start_read > two_exon_neg$end_exon_1 & two_exon_neg$end_read < two_exon_neg$start_exon_1){
#     two_exon_neg$occurs_in_exon = "start_exon_1"
#   } else if (!is.na(two_exon_neg$start_exon_2) & two_exon_neg$start_read > two_exon_neg$end_exon_2 & two_exon_neg$end_read < two_exon_neg$start_exon_2){
#     two_exon_neg$occurs_in_exon = "start_exon_2"
#   } else if (!is.na(two_exon_neg$start_exon_3) & two_exon_neg$start_read > two_exon_neg$end_exon_3 & two_exon_neg$end_read < two_exon_neg$start_exon_3){
#     two_exon_neg$occurs_in_exon = "start_exon_3"
#   } else if (!is.na(two_exon_neg$start_exon_4) & two_exon_neg$start_read > two_exon_neg$end_exon_4 & two_exon_neg$end_read < two_exon_neg$start_exon_4){
#     two_exon_neg$occurs_in_exon = "start_exon_4"
#   } else if (!is.na(two_exon_neg$start_exon_5) & two_exon_neg$start_read > two_exon_neg$end_exon_5 & two_exon_neg$end_read < two_exon_neg$start_exon_5){
#     two_exon_neg$occurs_in_exon = "start_exon_5"
#   } else if (!is.na(two_exon_neg$start_exon_6) & two_exon_neg$start_read > two_exon_neg$end_exon_6 & two_exon_neg$end_read < two_exon_neg$start_exon_6){
#     two_exon_neg$occurs_in_exon = "start_exon_6"
#   } else if (!is.na(two_exon_neg$start_exon_7) & two_exon_neg$start_read > two_exon_neg$end_exon_7 & two_exon_neg$end_read < two_exon_neg$start_exon_7){
#     two_exon_neg$occurs_in_exon = "start_exon_7"
#   } else if (!is.na(two_exon_neg$start_exon_8) & two_exon_neg$start_read > two_exon_neg$end_exon_8 & two_exon_neg$end_read < two_exon_neg$start_exon_8){
#     two_exon_neg$occurs_in_exon = "start_exon_8"
#   } else {
#     two_exon_neg$occurs_in_exon = "intronic"
#   }

negative_exon = function(generic_negative){
  if (!is.na(generic_negative$start_exon_1) & generic_negative$start_read > generic_negative$end_exon_1 & generic_negative$end_read < generic_negative$start_exon_1){
  generic_negative$occurs_in_exon = "start_exon_1"
  } else if (!is.na(generic_negative$start_exon_2) & generic_negative$start_read > generic_negative$end_exon_2 & generic_negative$end_read < generic_negative$start_exon_2){
  generic_negative$occurs_in_exon = "start_exon_2"
  } else if (!is.na(generic_negative$start_exon_3) & generic_negative$start_read > generic_negative$end_exon_3 & generic_negative$end_read < generic_negative$start_exon_3){
  generic_negative$occurs_in_exon = "start_exon_3"
  } else if (!is.na(generic_negative$start_exon_4) & generic_negative$start_read > generic_negative$end_exon_4 & generic_negative$end_read < generic_negative$start_exon_4){
  generic_negative$occurs_in_exon = "start_exon_4"
  } else if (!is.na(generic_negative$start_exon_5) & generic_negative$start_read > generic_negative$end_exon_5 & generic_negative$end_read < generic_negative$start_exon_5){
  generic_negative$occurs_in_exon = "start_exon_5"
  } else if (!is.na(generic_negative$start_exon_6) & generic_negative$start_read > generic_negative$end_exon_6 & generic_negative$end_read < generic_negative$start_exon_6){
  generic_negative$occurs_in_exon = "start_exon_6"
  } else if (!is.na(generic_negative$start_exon_7) & generic_negative$start_read > generic_negative$end_exon_7 & generic_negative$end_read < generic_negative$start_exon_7){
  generic_negative$occurs_in_exon = "start_exon_7"
  } else if (!is.na(generic_negative$start_exon_8) & generic_negative$start_read > generic_negative$end_exon_8 & generic_negative$end_read < generic_negative$start_exon_8){
  generic_negative$occurs_in_exon = "start_exon_8"
  } else {
  generic_negative$occurs_in_exon = "intronic"
}}

#now apply the above functions to the respective tables
neg_exons$occurs_in_exon = NA
pos_exons$occurs_in_exon = NA
head(pos_exons)
apply(pos_exons, 1, positive_exon)
?apply
positive_exon

head(neg_exons, n=1)
apply(head(neg_exons, n=1), 1, class)
