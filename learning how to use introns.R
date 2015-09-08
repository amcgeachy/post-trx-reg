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

library("dplyr")

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
head(pos)

#split exons from column with commas into distinct columns
#positive
  start_split_pos = strsplit(pos$exon_start, ",")
  
  empty_start_pos=matrix(data=NA, nrow=nrow(pos), ncol=8)
  colnames(empty_start_pos) = c(1:8)
  
  for (i in 1:nrow(pos)){
    for (j in 1:8){
      empty_start_pos[i,j]=as.numeric(start_split_pos[[i]][j])
      colnames(empty_start_pos)[j]=c(sprintf("start_%s", j))
    }
  }
  
  exon_size_split_pos = strsplit(pos$exon_size, ",")
  
  empty_exon_size_pos=matrix(data=NA, nrow=nrow(pos), ncol=8)
  colnames(empty_exon_size_pos) = c(1:8)
  
  for (i in 1:nrow(pos)){
    for (j in 1:8){
      empty_exon_size_pos[i,j]=as.numeric(exon_size_split_pos[[i]][j])
      colnames(empty_exon_size_pos)[j]=c(sprintf("exon_size_%s", j))
    }
  }
head(pos)
empty_start_pos  
  proc.time()

  pos_exons = cbind(pos, empty_start_pos, empty_exon_size_pos)
  head(pos_exons)
#negative

  start_split_neg = lapply(strsplit(neg$exon_start, ","), rev)

  empty_start_neg=matrix(data=NA, nrow=nrow(neg), ncol=8)
  colnames(empty_start_neg) = c(1:8)
  
  for (i in 1:nrow(neg)){
    for (j in 1:8){
      empty_start_neg[i,j]=as.numeric(start_split_neg[[i]][j])
      colnames(empty_start_neg)[j]=c(sprintf("start_%s", j))
    }
  }
  
  exon_size_split_neg = lapply(strsplit(neg$exon_size, ","), rev)
  
  empty_exon_size_neg=matrix(data=NA, nrow=nrow(neg), ncol=8)
  colnames(empty_exon_size_neg) = c(1:8)
  
  for (i in 1:nrow(neg)){
    for (j in 1:8){
      empty_exon_size_neg[i,j]=as.numeric(exon_size_split_neg[[i]][j])
      colnames(empty_exon_size_neg)[j]=c(sprintf("exon_size_%s", j))
    }
  }
  
  proc.time()
  
  neg_exons = cbind(neg, empty_start_neg, empty_exon_size_neg)
  
#define absolute genomic coordinates of exons
  starts=NULL
  ends=NULL

  for (i in 1:8){
    starts[i] = paste("start_exon_", i, sep="")
    ends[i] = paste("end_exon_", i, sep="")
  }
 
#for both positive and negative strand genes
  genic_exons = rbind(pos_exons, neg_exons)
  
    for (i in 1:8){
      genic_exons[,starts[i]] = genic_exons$start_cds + genic_exons[,sprintf("start_%s", i)]
      genic_exons[,ends[i]] = genic_exons$start_cds + genic_exons[,sprintf("start_%s", i)] + genic_exons[,sprintf("exon_size_%s", i)]
    }
  
#separate them again for checking which exon they occur in since math is slightly different
  neg_exons = filter(genic_exons, strand_read=="-")
  pos_exons = filter(genic_exons, strand_read=="+")

#function to check if exonic and if so which exon for positive strand
  positive_exon = function(generic_positive){
    for (i in 1:nrow(generic_positive)){
    if (!is.na(generic_positive[i,"start_exon_1"]) & (generic_positive[i,"start_read"] == generic_positive[i,"start_exon_1"] & generic_positive[i,"end_read"] <= generic_positive[i,"end_exon_1"])) { 
      generic_positive[i,"occurs_in_exon"] = "start_exon_1"
    } else if (!is.na(generic_positive[i,"start_exon_2"]) & (generic_positive[i,"start_read"] == generic_positive[i,"start_exon_2"] & generic_positive[i,"end_read"] <= generic_positive[i,"end_exon_2"])){
      generic_positive[i,"occurs_in_exon"] = "start_exon_2"
    } else if (!is.na(generic_positive[i,"start_exon_3"]) & (generic_positive[i,"start_read"] == generic_positive[i,"start_exon_3"] & generic_positive[i,"end_read"] <= generic_positive[i,"end_exon_3"])){
      generic_positive[i,"occurs_in_exon"] = "start_exon_3"
    } else if (!is.na(generic_positive[i,"start_exon_4"]) & (generic_positive[i,"start_read"] == generic_positive[i,"start_exon_4"] & generic_positive[i,"end_read"] <= generic_positive[i,"end_exon_4"])){
      generic_positive[i,"occurs_in_exon"] = "start_exon_4"
    } else if (!is.na(generic_positive[i,"start_exon_5"]) & (generic_positive[i,"start_read"] == generic_positive[i,"start_exon_5"] & generic_positive[i,"end_read"] <= generic_positive[i,"end_exon_5"])){
      generic_positive[i,"occurs_in_exon"] = "start_exon_5"
    } else if (!is.na(generic_positive[i,"start_exon_6"]) & (generic_positive[i,"start_read"] == generic_positive[i,"start_exon_6"] & generic_positive[i,"end_read"] <= generic_positive[i,"end_exon_6"])){
      generic_positive[i,"occurs_in_exon"] = "start_exon_6"
    } else if (!is.na(generic_positive[i,"start_exon_7"]) & (generic_positive[i,"start_read"] == generic_positive[i,"start_exon_7"] & generic_positive[i,"end_read"] <= generic_positive[i,"end_exon_7"])){
      generic_positive[i,"occurs_in_exon"] = "start_exon_7"
    } else if (!is.na(generic_positive[i,"start_exon_8"]) & (generic_positive[i,"start_read"] == generic_positive[i,"start_exon_8"] & generic_positive[i,"end_read"] <= generic_positive[i,"end_exon_8"])){
      generic_positive[i,"occurs_in_exon"] = "start_exon_8"
    } else {
      generic_positive[i,"occurs_in_exon"] = "intronic"
    }} 
    generic_positive }

itty = head(pos_exons, n=1)

(!is.na(itty[1,"start_exon_1"]) & (itty[1,"start_read"] == itty[1,"start_exon_1"] & itty[i,"end_read"] <= itty[i,"end_exon_1"])) { 
  
1 == 2

itty[1,"start_read"]
itty[1,"start_exon_1"]
pos_exons_copy = positive_exon(pos_exons)
table(pos_exons_copy$occurs_in_exon)
filter(pos_exons_copy, occurs_in_exon=="start_exon_1")
sample_n(pos_exons_copy, size=5)
# # FXN FROM JAKE THAT SHOULD BE FASTER THAN THE ABOVE, BUT I DONT KNOW HOW OR IF IT WORKS YET
# # fooBarFunc: Finds the first range that the input value of each row in inside of. Note this is an exclusive range (the
# # value cannot start or end at the ends of the range it is being compared with.
# #
# # df: the data frame of the input data
# # featureNameVector: a vector of the names of the features, prefixed by startPrefix and endPrefix to create the column names
# # startPrefix: The prefix prepended  to the feature name to create the column that represents the start of the feature
# # endPrefix: The prefix prepended  to the feature name to create the column that represents the end of the feature
# # testValueName: The name of the value that is under test, when combined with the prefix, 
# 
# fooBarFunc <- 
#   function(df, featureNameVector, startPrefix, endPrefix, testValueName) {
#     getStart <- function(valueString) {
#       paste(startPrefix, valueString, sep = "")
#     }
#     getEnd <- function(valueString) {
#       paste(endPrefix, valueString, sep = "")
#     }
#     features <- length(featureNameVector)
#     out <- data.frame(dropMe=rep(0, nrow(df))) # just creating the output data.frame with rows so its ok to append to
#     mask <- rep(T, nrow(df))
#     for (i in 1:features) {
#       temp <- mask & df[[getStart(featureNameVector[i])]] < df[[getStart(testValueName)]] & df[[getEnd(testValueName)]] < df[[getEnd(featureNameVector[i])]]
#       mask <- (mask & !temp)
#       out[[featureNameVector[i]]] <- ifelse(temp, i, 0)
#     }
#     sumOut <- rowSums(out)
#     intronicFeatureNames <- c("intronic", featureNameVector)
#     intronicFeatureNames[sumOut +1]
#   }


#function to check if exonic and if so which exon for negative strand
  negative_exon = function(generic_negative){
    for (i in 1:nrow(generic_negative)){
      if (!is.na(generic_negative[i,"start_exon_1"]) & generic_negative[i,"start_read"] >= generic_negative[i,"start_exon_1"] & generic_negative[i,"end_read"] <= generic_negative[i,"end_exon_1"]){
        generic_negative[i,"occurs_in_exon"] = "start_exon_1"
      } else if (!is.na(generic_negative[i,"start_exon_2"]) & generic_negative[i,"start_read"] >= generic_negative[i,"start_exon_2"] & generic_negative[i,"end_read"] <= generic_negative[i,"end_exon_2"]){
        generic_negative[i,"occurs_in_exon"] = "start_exon_2"
      } else if (!is.na(generic_negative[i,"start_exon_3"]) & generic_negative[i,"start_read"] >= generic_negative[i,"start_exon_3"] & generic_negative[i,"end_read"] <= generic_negative[i,"end_exon_3"]){
        generic_negative[i,"occurs_in_exon"] = "start_exon_3"
      } else if (!is.na(generic_negative[i,"start_exon_4"]) & generic_negative[i,"start_read"] >= generic_negative[i,"start_exon_4"] & generic_negative[i,"end_read"] <= generic_negative[i,"end_exon_4"]){
        generic_negative[i,"occurs_in_exon"] = "start_exon_4"
      } else if (!is.na(generic_negative[i,"start_exon_5"]) & generic_negative[i,"start_read"] >= generic_negative[i,"start_exon_5"] & generic_negative[i,"end_read"] <= generic_negative[i,"end_exon_5"]){
        generic_negative[i,"occurs_in_exon"] = "start_exon_5"
      } else if (!is.na(generic_negative[i,"start_exon_5"]) & generic_negative[i,"start_read"] >= generic_negative[i,"start_exon_6"] & generic_negative[i,"end_read"] <= generic_negative[i,"end_exon_5"]){
        generic_negative[i,"occurs_in_exon"] = "start_exon_6"
      } else if (!is.na(generic_negative[i,"start_exon_7"]) & generic_negative[i,"start_read"] >= generic_negative[i,"start_exon_7"] & generic_negative[i,"end_read"] <= generic_negative[i,"end_exon_7"]){
        generic_negative[i,"occurs_in_exon"] = "start_exon_7"
      } else if (!is.na(generic_negative[i,"start_exon_8"]) & generic_negative[i,"start_read"] >= generic_negative[i,"start_exon_8"] & generic_negative[i,"end_read"] <= generic_negative[i,"end_exon_8"]){
        generic_negative[i,"occurs_in_exon"] = "start_exon_8"
      } else {
        generic_negative[i,"occurs_in_exon"] = "intronic"
      }}
    generic_negative}

neg_exons_copy = negative_exon(neg_exons)
table(neg_exons_copy$occurs_in_exon)

## now we need to see if the exonic fragments are in frame or not
# first, pull all the exonic fragments for positive strand
test_pos = filter(pos_exons_copy, occurs_in_exon!="intronic")
nrow(test_pos)

pos_intronic = filter(pos_exons_copy, occurs_in_exon=="intronic")
nrow(pos_intronic)

pos_tiny = rbind(head(filter(test_pos, occurs_in_exon=="start_exon_1"), n=1),
                 head(filter(test_pos, occurs_in_exon=="start_exon_2"), n=1))

pos_tiny

CATGTTGCTGAGTGAACTCGTAGCAACCGCCTCCTCTCTGCCATACACGG
CCATCAGCATCCACAATAACTGTCGTGTCCCAGCCGCACGCCACATCCAC
CACGGGTGCCGGTACTTCCACGGGCCTCCAGTCATGCACCTGCCGCAGTG
CTTGCGCACTATCCAGTTCTCCCCGTCTGTTATCTCCACATCCTACCAGA
TTCCCGTCATTTGTCAGCATCACGCTGTGGTTCCCACCGCACGCTATCTT
CCTGACTATTGCTCCATCATCTCCTGGCACAGACCTCTGTGGGGTATCCA
TATCCTCATCGTGCCCCAGTCCCAGTTGCCTTTGCCCATTAGACCCAAAC
GCATACACACAACTCATCGATACAAGCCTGTTATAGCCTTTAATGATCAC
ATTCCATCACTTGCGCTTTGGATCTGCCTGCATTATCAAGGCTCAAACGG
CTGCGTTACCCCCGTCGCCGCGAAATTTTTCATAATTTTTCACTTTGTAG
GATTAAAAGAGATCATGAGCCCATCTCGCAATGCAACACGTAACTTAAAT
CAGTACTGGCGTGTGCTATAG
#"start" 114249
A 114250 T 114251 G 114252 

#end 114819
T 114817 A 114818 G 114819


# first, pull all the exonic fragments for negative strand
test = filter(neg_exons_copy, occurs_in_exon!="intronic")
nrow(test)
head(test)

neg_intronic = filter(neg_exons_copy, occurs_in_exon=="intronic")
nrow(neg_intronic)

neg_tiny = rbind(head(filter(test, occurs_in_exon=="start_exon_1"), n=1),
 head(filter(test, occurs_in_exon=="start_exon_2"), n=1))

neg_tiny

#now that we have absolute genomic coordinates for the cds we can determine frame information
#figuring out how to do math for what's in frame for negative strand
#START
#G 101143 T 101144 A 101145
#2        1        0  -- intuitive/biological frame
#1        2        0  -- computation frame (end_read - start_exon)
(101145 - 101145 ) %% 3 #A, in frame --> 0
(101144 - 101145 ) %% 3 #T 
(101143 - 101145 ) %% 3 #G

#END
#A 100224 A 100223 T 100222
#2        1        0  -- intuitive/biological frame
#1        0        2  -- computation frame (start_read - end_exon)
(100222 - 100224 + 1) %% 3 #T
(100223 - 100224 + 1) %% 3 #A
(100224 - 100224 + 1) %% 3 #A, in frame --> 1


#at the start 
neg_tiny[1,"read_start_frame"] = (neg_tiny[1,"end_read"] - neg_tiny[1,neg_tiny[1,"occurs_in_exon"]] )%%3
neg_tiny[2,"read_start_frame"] = (neg_tiny[2,"end_read"] - neg_tiny[2,neg_tiny[2,"occurs_in_exon"]] )%%3

#...and at the end
neg_tiny[1,"read_end_frame"] = (neg_tiny[1,"start_read"] - neg_tiny[1,sprintf("end_exon_%s", unlist(strsplit(neg_tiny[1,"occurs_in_exon"], "_"))[3])])%%3
neg_tiny[2,"read_end_frame"] = (neg_tiny[2,"start_read"] - neg_tiny[2,sprintf("end_exon_%s", unlist(strsplit(neg_tiny[2,"occurs_in_exon"], "_"))[3])])%%3


seq(from=312950, to=313362, by=3)
neg_tiny
(100224 - 100567 )%% 3
seq(from=100224, to=100567, by=3)
seq(from=100224, to=101145, by=3)
?seq

?split
split(neg_tiny[1,"occurs_in_exon"], "_")
#gets the exon number from exon determination (start_exon_#) and writes corresponding end exon column name
sprintf("end_exon_%s", unlist(strsplit(neg_tiny[1,"occurs_in_exon"], "_"))[3])


neg_tiny(neg_tiny[1,"start_exon_1"] - neg_tiny[1,"end_read"])%%3
neg_tiny[2,"occurs_in_exon"]


(neg_tiny[1,"start_exon_1"] - neg_tiny[1,"start_exon_1"]) %% 3
(neg_tiny[1,"start_exon_1"] - neg_tiny[1,"start_exon_1"] + 1) %% 3
(neg_tiny[1,"start_exon_1"] - neg_tiny[1,"start_exon_1"] + 2) %% 3
