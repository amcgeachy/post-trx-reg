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
 
#for both positive genes  
    for (i in 1:8){
      pos_exons[,starts[i]] = pos_exons$start_cds + pos_exons[,sprintf("start_%s", i)] + 1 
      #this +1 is important because start_cds = n-1, where n == A of ATG (eg YAL003W ATG has A at 142174 but start cds is 142713; seq shown below in frameness)
      pos_exons[,ends[i]] = pos_exons$start_cds + pos_exons[,sprintf("start_%s", i)] + pos_exons[,sprintf("exon_size_%s", i)]
    }
  
#and for negative strand genes
    for (i in 1:8){
      neg_exons[,ends[i]] = neg_exons$start_cds + neg_exons[,sprintf("start_%s", i)] + 1 
      neg_exons[,starts[i]] = neg_exons$start_cds + neg_exons[,sprintf("start_%s", i)] + neg_exons[,sprintf("exon_size_%s", i)]   
      #this +1 is important because start_cds = n-1, where n == A of ATG (eg YAL003W ATG has A at 142174 but start cds is 142713; seq shown below in frameness)
    }

#function to check if exonic and if so which exon for positive strand
  positive_exon = function(generic_positive){
    for (i in 1:nrow(generic_positive)){
    if (!is.na(generic_positive[i,"start_exon_1"]) & (generic_positive[i,"start_read"] >= (generic_positive[i,"start_exon_1"] -1) & generic_positive[i,"end_read"] <= generic_positive[i,"end_exon_1"])) { 
      generic_positive[i,"occurs_in_exon"] = "start_exon_1"
    } else if (!is.na(generic_positive[i,"start_exon_2"]) & (generic_positive[i,"start_read"] >= (generic_positive[i,"start_exon_2"] -1) & generic_positive[i,"end_read"] <= generic_positive[i,"end_exon_2"])){
      generic_positive[i,"occurs_in_exon"] = "start_exon_2"
    } else if (!is.na(generic_positive[i,"start_exon_3"]) & (generic_positive[i,"start_read"] >= (generic_positive[i,"start_exon_3"] -1) & generic_positive[i,"end_read"] <= generic_positive[i,"end_exon_3"])){
      generic_positive[i,"occurs_in_exon"] = "start_exon_3"
    } else if (!is.na(generic_positive[i,"start_exon_4"]) & (generic_positive[i,"start_read"] >= (generic_positive[i,"start_exon_4"] -1) & generic_positive[i,"end_read"] <= generic_positive[i,"end_exon_4"])){
      generic_positive[i,"occurs_in_exon"] = "start_exon_4"
    } else if (!is.na(generic_positive[i,"start_exon_5"]) & (generic_positive[i,"start_read"] >= (generic_positive[i,"start_exon_5"] -1) & generic_positive[i,"end_read"] <= generic_positive[i,"end_exon_5"])){
      generic_positive[i,"occurs_in_exon"] = "start_exon_5"
    } else if (!is.na(generic_positive[i,"start_exon_6"]) & (generic_positive[i,"start_read"] >= (generic_positive[i,"start_exon_6"] -1) & generic_positive[i,"end_read"] <= generic_positive[i,"end_exon_6"])){
      generic_positive[i,"occurs_in_exon"] = "start_exon_6"
    } else if (!is.na(generic_positive[i,"start_exon_7"]) & (generic_positive[i,"start_read"] >= (generic_positive[i,"start_exon_7"] -1) & generic_positive[i,"end_read"] <= generic_positive[i,"end_exon_7"])){
      generic_positive[i,"occurs_in_exon"] = "start_exon_7"
    } else if (!is.na(generic_positive[i,"start_exon_8"]) & (generic_positive[i,"start_read"] >= (generic_positive[i,"start_exon_8"] -1) & generic_positive[i,"end_read"] <= generic_positive[i,"end_exon_8"])){
      generic_positive[i,"occurs_in_exon"] = "start_exon_8"
    } else {
      generic_positive[i,"occurs_in_exon"] = "intronic"
    }} 
    generic_positive }

pos_exons_copy = positive_exon(pos_exons)
table(pos_exons_copy$occurs_in_exon)

#function to check if exonic and if so which exon for negative strand
negative_exon = function(generic_negative){
  for (i in 1:nrow(generic_negative)){
    if (!is.na(generic_negative[i,"start_exon_1"]) & generic_negative[i,"end_read"] <= (generic_negative[i,"start_exon_1"] +1) & generic_negative[i,"start_read"] >= generic_negative[i,"end_exon_1"]){
      generic_negative[i,"occurs_in_exon"] = "start_exon_1"
    } else if (!is.na(generic_negative[i,"start_exon_2"]) & generic_negative[i,"end_read"] <= (generic_negative[i,"start_exon_2"] +1) & generic_negative[i,"start_read"] >= generic_negative[i,"end_exon_2"]){
      generic_negative[i,"occurs_in_exon"] = "start_exon_2"
    } else if (!is.na(generic_negative[i,"start_exon_3"]) & generic_negative[i,"end_read"] <= (generic_negative[i,"start_exon_3"] +1) & generic_negative[i,"start_read"] >= generic_negative[i,"end_exon_3"]){
      generic_negative[i,"occurs_in_exon"] = "start_exon_3"
    } else if (!is.na(generic_negative[i,"start_exon_4"]) & generic_negative[i,"end_read"] <= (generic_negative[i,"start_exon_4"] +1) & generic_negative[i,"start_read"] >= generic_negative[i,"end_exon_4"]){
      generic_negative[i,"occurs_in_exon"] = "start_exon_4"
    } else if (!is.na(generic_negative[i,"start_exon_5"]) & generic_negative[i,"end_read"] <= (generic_negative[i,"start_exon_5"] +1) & generic_negative[i,"start_read"] >= generic_negative[i,"end_exon_5"]){
      generic_negative[i,"occurs_in_exon"] = "start_exon_5"
    } else if (!is.na(generic_negative[i,"start_exon_5"]) & generic_negative[i,"end_read"] <= (generic_negative[i,"start_exon_6"] +1) & generic_negative[i,"start_read"] >= generic_negative[i,"end_exon_5"]){
      generic_negative[i,"occurs_in_exon"] = "start_exon_6"
    } else if (!is.na(generic_negative[i,"start_exon_7"]) & generic_negative[i,"end_read"] <= (generic_negative[i,"start_exon_7"] +1) & generic_negative[i,"start_read"] >= generic_negative[i,"end_exon_7"]){
      generic_negative[i,"occurs_in_exon"] = "start_exon_7"
    } else if (!is.na(generic_negative[i,"start_exon_8"]) & generic_negative[i,"end_read"] <= (generic_negative[i,"start_exon_8"] +1) & generic_negative[i,"start_read"] >= generic_negative[i,"end_exon_8"]){
      generic_negative[i,"occurs_in_exon"] = "start_exon_8"
    } else {
      generic_negative[i,"occurs_in_exon"] = "intronic"
    }}
  generic_negative}

neg_exons_copy = negative_exon(neg_exons)
table(neg_exons_copy$occurs_in_exon)

## now we need to see if the exonic fragments are in frame or not
#pull exonic fragments
pos_exonic = filter(pos_exons_copy, occurs_in_exon!="intronic")
pos_intronic = filter(pos_exons_copy, occurs_in_exon=="intronic") #just saving these for later
head(pos_exonic)

neg_exonic = filter(neg_exons_copy, occurs_in_exon!="intronic")
neg_intronic = filter(neg_exons_copy, occurs_in_exon=="intronic") #just saving these for later
head(neg_exonic)

#for positive strand genes
#we'll do this by utilizing the idea that a reads CDS coordinates can be defined as:
# start read CDS loc == sum(exon sizes from 1 to exon before the read occurs in) + (start read genomic location - start exon genomic location for the exon the read is in)
# or when s is size, e is exon start site, g is start read genomic location, CDS location is c
# there are 1 to n exons in a gene, and the read is in exon j
# c = Σ(s[sub 1:j-1]) + (g - e[sub j])

  #sum previous exons
  for (i in 1:nrow(pos_exonic)){
    if (pos_exonic[i,"occurs_in_exon"]=="start_exon_1"){
      pos_exonic[i,"prior_exon_sums"] = 0
    } else if (pos_exonic[i,"occurs_in_exon"]=="start_exon_2"){
      pos_exonic[i,"prior_exon_sums"] = pos_exonic[i,"exon_size_1"]
    } else if (pos_exonic[i,"occurs_in_exon"]=="start_exon_3"){
      pos_exonic[i,"prior_exon_sums"] = pos_exonic[i,"exon_size_1"] + pos_exonic[i,"exon_size_2"]
    } else if (pos_exonic[i,"occurs_in_exon"]=="start_exon_4"){
      pos_exonic[i,"prior_exon_sums"] = pos_exonic[i,"exon_size_1"] + pos_exonic[i,"exon_size_2"] + pos_exonic[i,"exon_size_3"] 
    } else if (pos_exonic[i,"occurs_in_exon"]=="start_exon_5"){
      pos_exonic[i,"prior_exon_sums"] = pos_exonic[i,"exon_size_1"] + pos_exonic[i,"exon_size_2"] + pos_exonic[i,"exon_size_3"] + pos_exonic[i,"exon_size_4"] 
    } else if (pos_exonic[i,"occurs_in_exon"]=="start_exon_6"){
      pos_exonic[i,"prior_exon_sums"] = pos_exonic[i,"exon_size_1"] + pos_exonic[i,"exon_size_2"] + pos_exonic[i,"exon_size_3"] + pos_exonic[i,"exon_size_4"] + pos_exonic[i,"exon_size_5"] 
    } else if (pos_exonic[i,"occurs_in_exon"]=="start_exon_7"){
      pos_exonic[i,"prior_exon_sums"] = pos_exonic[i,"exon_size_1"] + pos_exonic[i,"exon_size_2"] + pos_exonic[i,"exon_size_3"] + pos_exonic[i,"exon_size_4"] + pos_exonic[i,"exon_size_5"] + pos_exonic[i,"exon_size_6"]
    } else if (pos_exonic[i,"occurs_in_exon"]=="start_exon_8"){
      pos_exonic[i,"prior_exon_sums"] = pos_exonic[i,"exon_size_1"] + pos_exonic[i,"exon_size_2"] + pos_exonic[i,"exon_size_3"] + pos_exonic[i,"exon_size_4"] + pos_exonic[i,"exon_size_5"] + pos_exonic[i,"exon_size_6"] + pos_exonic[i,"exon_size_7"]
    } else {
      pos_exonic[i,"prior_exon_sums"] = NA
    }}
  
  head(pos_exonic)

  #find start and end location in current exon
  for (i in 1:nrow(pos_exonic)){
    pos_exonic[i,"read_start_in_exon"] = pos_exonic[i,"start_read"] - pos_exonic[i,pos_exonic[i,"occurs_in_exon"]] + 1 #g[start] - e[sub j], +1 is b/c oddities of 0 counting
    pos_exonic[i,"read_end_in_exon"] = pos_exonic[i,"end_read"] - pos_exonic[i,pos_exonic[i,"occurs_in_exon"]] + 1 #g[end] - e[sub j]
  }
  
  head(pos_exonic)
  
  pos_exonic$read_start_in_cds = pos_exonic$prior_exon_sums + pos_exonic$read_start_in_exon 
  pos_exonic$read_end_in_cds = pos_exonic$prior_exon_sums + pos_exonic$read_end_in_exon 
  
  pos_exonic$read_start_cds_frame = (pos_exonic$read_start_in_cds - 1) %% 3 #because counting in 0 space
  pos_exonic$read_end_cds_frame = (pos_exonic$read_end_in_cds - 1) %% 3 #because counting in 0 space
  
  table(pos_exonic$read_start_cds_frame)
  filter(pos_exonic, start_read==351124) # one example that checks out
  filter(pos_exonic, start_read==189856) # another example thats checks out for frame info
  head(pos_exonic)[1,]

#for negative strand genes
#we'll do this by utilizing the idea that a reads CDS coordinates can be defined as:
# start read CDS loc == sum(exon sizes from 1 to exon before the read occurs in) + (start exon genomic location for the exon the read is in - end read genomic location)
# or when s is size, e is exon start site, g is start read genomic location, CDS location is c
# there are 1 to n exons in a gene, and the read is in exon j
# c = Σ(s[sub 1:j-1]) + (e[sub j] - g)

#summing the preceding exon sizes should be the same
  for (i in 1:nrow(neg_exonic)){
    if (neg_exonic[i,"occurs_in_exon"]=="start_exon_1"){
      neg_exonic[i,"prior_exon_sums"] = 0
    } else if (neg_exonic[i,"occurs_in_exon"]=="start_exon_2"){
      neg_exonic[i,"prior_exon_sums"] = neg_exonic[i,"exon_size_1"]
    } else if (neg_exonic[i,"occurs_in_exon"]=="start_exon_3"){
      neg_exonic[i,"prior_exon_sums"] = neg_exonic[i,"exon_size_1"] + neg_exonic[i,"exon_size_2"]
    } else if (neg_exonic[i,"occurs_in_exon"]=="start_exon_4"){
      neg_exonic[i,"prior_exon_sums"] = neg_exonic[i,"exon_size_1"] + neg_exonic[i,"exon_size_2"] + neg_exonic[i,"exon_size_3"] 
    } else if (neg_exonic[i,"occurs_in_exon"]=="start_exon_5"){
      neg_exonic[i,"prior_exon_sums"] = neg_exonic[i,"exon_size_1"] + neg_exonic[i,"exon_size_2"] + neg_exonic[i,"exon_size_3"] + neg_exonic[i,"exon_size_4"] 
    } else if (neg_exonic[i,"occurs_in_exon"]=="start_exon_6"){
      neg_exonic[i,"prior_exon_sums"] = neg_exonic[i,"exon_size_1"] + neg_exonic[i,"exon_size_2"] + neg_exonic[i,"exon_size_3"] + neg_exonic[i,"exon_size_4"] + neg_exonic[i,"exon_size_5"] 
    } else if (neg_exonic[i,"occurs_in_exon"]=="start_exon_7"){
      neg_exonic[i,"prior_exon_sums"] = neg_exonic[i,"exon_size_1"] + neg_exonic[i,"exon_size_2"] + neg_exonic[i,"exon_size_3"] + neg_exonic[i,"exon_size_4"] + neg_exonic[i,"exon_size_5"] + neg_exonic[i,"exon_size_6"]
    } else if (neg_exonic[i,"occurs_in_exon"]=="start_exon_8"){
      neg_exonic[i,"prior_exon_sums"] = neg_exonic[i,"exon_size_1"] + neg_exonic[i,"exon_size_2"] + neg_exonic[i,"exon_size_3"] + neg_exonic[i,"exon_size_4"] + neg_exonic[i,"exon_size_5"] + neg_exonic[i,"exon_size_6"] + neg_exonic[i,"exon_size_7"]
    } else {
      neg_exonic[i,"prior_exon_sums"] = NA
    }}
  
  head(neg_exonic) #check to make sure things look good

  #find location in current exon, order of g and e different than positive strand
  for (i in 1:nrow(neg_exonic)){
    neg_exonic[i,"read_start_in_exon"] = neg_exonic[i,neg_exonic[i,"occurs_in_exon"]] - neg_exonic[i,"end_read"] + 1 # e[sub j] - g[start] , +1 is b/c oddities of 0 counting
    neg_exonic[i,"read_end_in_exon"] = neg_exonic[i,neg_exonic[i,"occurs_in_exon"]] - neg_exonic[i,"start_read"] + 1 #e[sub j] - g[end], +1 is b/c oddities of 0 counting
  }

  head(neg_exonic)
  sample_n(filter(neg_exonic, exon_number!=1), size=1)
  
  #translate position in exon to position in CDS, same as in positive strand?
  neg_exonic$read_start_in_cds = neg_exonic$prior_exon_sums + neg_exonic$read_start_in_exon 
  neg_exonic$read_end_in_cds = neg_exonic$prior_exon_sums + neg_exonic$read_end_in_exon
  
  neg_exonic$read_start_cds_frame = (neg_exonic$read_start_in_cds - 1) %% 3 #because counting in 0 space
  neg_exonic$read_end_cds_frame = (neg_exonic$read_end_in_cds - 1) %% 3 #because counting in 0 space
  
  #head(neg_exonic, n=1) #YAL025C with start at 100399 seems to check out
  #filter(neg_exonic, gene_name=="YHR203C", start_read==504614) #also seems to check out, woo.
  table(neg_exonic$read_start_cds_frame)

#now lets try to see how the framing distribution works out (0,0 vs 0,1, vs 0,2... etc)
head(pos_exonic)
pos_exonic$joint_frame = paste(pos_exonic$read_start_cds_frame, pos_exonic$read_end_cds_frame, sep=",")
table(pos_exonic$joint_frame)
barplot(table(pos_exonic$joint_frame))

head(neg_exonic)
neg_exonic$joint_frame = paste(neg_exonic$read_start_cds_frame, neg_exonic$read_end_cds_frame, sep=",")
table(neg_exonic$joint_frame)
barplot(table(neg_exonic$joint_frame))

filter(pos_exonic, joint_frame=="2,1")
filter(pos_exonic, gene_name=="YAL017W", start_read==121508)
filter(pos_exonic, gene_name=="YFR019W", read_start_in_cds==5355)
sum(filter(pos_exonic, joint_frame=="2,1")[,"frag_count"])
sum(filter(pos_exonic, joint_frame=="0,2")[,"frag_count"])
hist(filter(pos_exonic, joint_frame=="2,1")[,"frag_count"])

sample_n(filter(pos_exonic, joint_frame=="2,1"), size=1)
table(filter(pos_exonic, gene_name=="YPR204W")[,"joint_frame"])
filter()