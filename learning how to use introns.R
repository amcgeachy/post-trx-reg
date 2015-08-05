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

  #function to check if exonic and if so which exon for positive strand
#   positive_exon = function(generic_positive){
#     if (generic_positive$start_read > generic_positive$start_exon_1 & generic_positive$end_read < generic_positive$end_exon_1) { 
#       generic_positive$occurs_in_exon = "start_exon_1"
#     } else if (!is.na(generic_positive$start_exon_2) & (generic_positive$start_read > generic_positive$start_exon_2 & generic_positive$end_read < generic_positive$end_exon_2)){
#       generic_positive$occurs_in_exon = "start_exon_2"
#     } else if (!is.na(generic_positive$start_exon_3) & (generic_positive$start_read > generic_positive$start_exon_3 & generic_positive$end_read < generic_positive$end_exon_3)){
#       generic_positive$occurs_in_exon = "start_exon_3"
#     } else if (!is.na(generic_positive$start_exon_4) & (generic_positive$start_read > generic_positive$start_exon_4 & generic_positive$end_read < generic_positive$end_exon_4)){
#       generic_positive$occurs_in_exon = "start_exon_4"
#     } else if (!is.na(generic_positive$start_exon_5) & (generic_positive$start_read > generic_positive$start_exon_5 & generic_positive$end_read < generic_positive$end_exon_5)){
#       generic_positive$occurs_in_exon = "start_exon_5"
#     } else if (!is.na(generic_positive$start_exon_6) & (generic_positive$start_read > generic_positive$start_exon_6 & generic_positive$end_read < generic_positive$end_exon_6)){
#       generic_positive$occurs_in_exon = "start_exon_6"
#     } else if (!is.na(generic_positive$start_exon_7) & (generic_positive$start_read > generic_positive$start_exon_7 & generic_positive$end_read < generic_positive$end_exon_7)){
#       generic_positive$occurs_in_exon = "start_exon_7"
#     } else if (!is.na(generic_positive$start_exon_8) & (generic_positive$start_read > generic_positive$start_exon_8 & generic_positive$end_read < generic_positive$end_exon_8)){
#       generic_positive$occurs_in_exon = "start_exon_8"
#     } else {
#       generic_positive$occurs_in_exon = "intronic"
#     }}

positive_exon = function(generic_positive){
  for (i in 1:nrow(generic_positive)){
  if (generic_positive[i,"start_read"] > generic_positive[i,"start_exon_1"] & generic_positive[i,"end_read"] < generic_positive[i,"end_exon_1"]) { 
    generic_positive[i,"occurs_in_exon"] = "start_exon_1"
  } else if (!is.na(generic_positive[i,"start_exon_2"]) & (generic_positive[i,"start_read"] > generic_positive[i,"start_exon_2"] & generic_positive[i,"end_read"] < generic_positive[i,"end_exon_2"])){
    generic_positive[i,"occurs_in_exon"] = "start_exon_2"
  } else if (!is.na(generic_positive[i,"start_exon_3"]) & (generic_positive[i,"start_read"] > generic_positive[i,"start_exon_3"] & generic_positive[i,"end_read"] < generic_positive[i,"end_exon_3"])){
    generic_positive[i,"occurs_in_exon"] = "start_exon_3"
  } else if (!is.na(generic_positive[i,"start_exon_4"]) & (generic_positive[i,"start_read"] > generic_positive[i,"start_exon_4"] & generic_positive[i,"end_read"] < generic_positive[i,"end_exon_4"])){
    generic_positive[i,"occurs_in_exon"] = "start_exon_4"
  } else if (!is.na(generic_positive[i,"start_exon_5"]) & (generic_positive[i,"start_read"] > generic_positive[i,"start_exon_5"] & generic_positive[i,"end_read"] < generic_positive[i,"end_exon_5"])){
    generic_positive[i,"occurs_in_exon"] = "start_exon_5"
  } else if (!is.na(generic_positive[i,"start_exon_6"]) & (generic_positive[i,"start_read"] > generic_positive[i,"start_exon_6"] & generic_positive[i,"end_read"] < generic_positive[i,"end_exon_6"])){
    generic_positive[i,"occurs_in_exon"] = "start_exon_6"
  } else if (!is.na(generic_positive[i,"start_exon_7"]) & (generic_positive[i,"start_read"] > generic_positive[i,"start_exon_7"] & generic_positive[i,"end_read"] < generic_positive[i,"end_exon_7"])){
    generic_positive[i,"occurs_in_exon"] = "start_exon_7"
  } else if (!is.na(generic_positive[i,"start_exon_8"]) & (generic_positive[i,"start_read"] > generic_positive[i,"start_exon_8"] & generic_positive[i,"end_read"] < generic_positive[i,"end_exon_8"])){
    generic_positive[i,"occurs_in_exon"] = "start_exon_8"
  } else {
    generic_positive[i,"occurs_in_exon"] = "intronic"
  }} 
  generic_positive }
tester = sample_n(pos_exons, size = 5)
complex = rbind(head(filter(pos_exons, exon_number==2), n=2),
head(filter(pos_exons, exon_number==3), n=1),
head(filter(pos_exons, exon_number==4), n=1),
head(filter(pos_exons, exon_number==1), n=1))

blah =  positive_exon(tester)
nrow(tester)
tester
blah
positive_exon(complex)


#
# fooBarFunc: Finds the first range that the input value of each row in inside of. Note this is an exclusive range (the
# value cannot start or end at the ends of the range it is being compared with.
#
# df: the data frame of the input data
# featureNameVector: a vector of the names of the features, prefixed by startPrefix and endPrefix to create the column names
# startPrefix: The prefix prepended  to the feature name to create the column that represents the start of the feature
# endPrefix: The prefix prepended  to the feature name to create the column that represents the end of the feature
# testValueName: The name of the value that is under test, when combined with the prefix, 

fooBarFunc <- 
  function(df, featureNameVector, startPrefix, endPrefix, testValueName) {
    getStart <- function(valueString) {
      paste(startPrefix, valueString, sep = "")
    }
    getEnd <- function(valueString) {
      paste(endPrefix, valueString, sep = "")
    }
    features <- length(featureNameVector)
    out <- data.frame(dropMe=rep(0, nrow(df))) # just creating the output data.frame with rows so its ok to append to
    mask <- rep(T, nrow(df))
    for (i in 1:features) {
      temp <- mask & df[[getStart(featureNameVector[i])]] < df[[getStart(testValueName)]] & df[[getEnd(testValueName)]] < df[[getEnd(featureNameVector[i])]]
      mask <- (mask & !temp)
      out[[featureNameVector[i]]] <- ifelse(temp, i, 0)
    }
    sumOut <- rowSums(out)
    intronicFeatureNames <- c("intronic", featureNameVector)
    intronicFeatureNames[sumOut +1]
  }


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
