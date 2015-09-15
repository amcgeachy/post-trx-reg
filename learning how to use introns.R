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
      genic_exons[,starts[i]] = genic_exons$start_cds + genic_exons[,sprintf("start_%s", i)] + 1 
      #this +1 is important because start_cds = n-1, where n == A of ATG (eg YAL003W ATG has A at 142174 but start cds is 142713; seq shown below in frameness)
      genic_exons[,ends[i]] = genic_exons$start_cds + genic_exons[,sprintf("start_%s", i)] + genic_exons[,sprintf("exon_size_%s", i)]
    }
  
#separate them again for checking which exon they occur in since math is slightly different
  neg_exons = filter(genic_exons, strand_read=="-")
  pos_exons = filter(genic_exons, strand_read=="+")

#function to check if exonic and if so which exon for positive strand
  positive_exon = function(generic_positive){
    for (i in 1:nrow(generic_positive)){
    if (!is.na(generic_positive[i,"start_exon_1"]) & (generic_positive[i,"start_read"] >= generic_positive[i,"start_exon_1"] & generic_positive[i,"end_read"] <= generic_positive[i,"end_exon_1"])) { 
      generic_positive[i,"occurs_in_exon"] = "start_exon_1"
    } else if (!is.na(generic_positive[i,"start_exon_2"]) & (generic_positive[i,"start_read"] >= generic_positive[i,"start_exon_2"] & generic_positive[i,"end_read"] <= generic_positive[i,"end_exon_2"])){
      generic_positive[i,"occurs_in_exon"] = "start_exon_2"
    } else if (!is.na(generic_positive[i,"start_exon_3"]) & (generic_positive[i,"start_read"] >= generic_positive[i,"start_exon_3"] & generic_positive[i,"end_read"] <= generic_positive[i,"end_exon_3"])){
      generic_positive[i,"occurs_in_exon"] = "start_exon_3"
    } else if (!is.na(generic_positive[i,"start_exon_4"]) & (generic_positive[i,"start_read"] >= generic_positive[i,"start_exon_4"] & generic_positive[i,"end_read"] <= generic_positive[i,"end_exon_4"])){
      generic_positive[i,"occurs_in_exon"] = "start_exon_4"
    } else if (!is.na(generic_positive[i,"start_exon_5"]) & (generic_positive[i,"start_read"] >= generic_positive[i,"start_exon_5"] & generic_positive[i,"end_read"] <= generic_positive[i,"end_exon_5"])){
      generic_positive[i,"occurs_in_exon"] = "start_exon_5"
    } else if (!is.na(generic_positive[i,"start_exon_6"]) & (generic_positive[i,"start_read"] >= generic_positive[i,"start_exon_6"] & generic_positive[i,"end_read"] <= generic_positive[i,"end_exon_6"])){
      generic_positive[i,"occurs_in_exon"] = "start_exon_6"
    } else if (!is.na(generic_positive[i,"start_exon_7"]) & (generic_positive[i,"start_read"] >= generic_positive[i,"start_exon_7"] & generic_positive[i,"end_read"] <= generic_positive[i,"end_exon_7"])){
      generic_positive[i,"occurs_in_exon"] = "start_exon_7"
    } else if (!is.na(generic_positive[i,"start_exon_8"]) & (generic_positive[i,"start_read"] >= generic_positive[i,"start_exon_8"] & generic_positive[i,"end_read"] <= generic_positive[i,"end_exon_8"])){
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
#pull exonic fragments
pos_exonic = filter(pos_exons_copy, occurs_in_exon!="intronic")
head(pos_exonic)

neg_exonic = filter(neg_exons_copy, occurs_in_exon!="intronic")
head(neg_exonic)


#for positive strand genes
#we'll do this by utilizing the idea that a reads CDS coordinates can be defined as:
# start read CDS loc == sum(exon sizes from 1 to exon before the read occurs in) + (start read genomic location - start exon genomic location for the exon the read is in)
# or when s is size, e is exon start site, g is start read genomic location, CDS location is c
# there are 1 to n exons in a gene, and the read is in exon j
# c = Σ(s[sub 1:j-1]) + (g - e[sub j])

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
    pos_exonic[i,"read_start_in_exon"] = g_minus_e_start = pos_exonic[i,"start_read"] - pos_exonic[i,pos_exonic[i,"occurs_in_exon"]] + 1 #g[start] - e[sub j], +1 is b/c oddities of 0 counting
    pos_exonic[i,"read_end_in_exon"] = pos_exonic[i,"end_read"] - pos_exonic[i,pos_exonic[i,"occurs_in_exon"]] #g[end] - e[sub j]
  }
  
  head(pos_exonic)
  
  pos_exonic$read_start_in_cds = pos_exonic$prior_exon_sums + pos_exonic$read_start_in_exon 
  pos_exonic$read_end_in_cds = pos_exonic$prior_exon_sums + pos_exonic$read_end_in_exon
  
  pos_exonic$read_start_cds_frame = (pos_exonic$read_start_in_cds - 1) %% 3
  pos_exonic$read_end_cds_frame = pos_exonic$read_end_in_cds %% 3
  
  table(pos_exonic$read_start_cds_frame)
  #filter(pos_exonic, start_read==351124) # one example that checks out


#and what we *really* care about is if it's in frame

# ## BUT exons can enter or leave in ANY frame, so we need to figure out what the frame is of each exon (UGH)
# ## start by defining exons in transcript coordinates
# 
# #define columns that will exist
# start_exon_n = NULL
# start_exon_in_trans_n = NULL
# end_exon_n = NULL
# end_exon_in_trans_n = NULL
# for (i in 1:8){
#   start_exon_n[i] = sprintf("start_exon_%s", i)
#   start_exon_in_trans_n[i] = sprintf("start_exon_in_trans_%s", i)
#   end_exon_n[i] = sprintf("end_exon_%s", i)
#   end_exon_in_trans_n[i] = sprintf("end_exon_in_trans_%s", i)
# }
# 
# #for positive exons, it goes like this
# exons_in_trans_coords = function(gen_pos){
#   for (i in 1:8){
#     for (j in 1:nrow(gen_pos)){
#     gen_pos[j,start_exon_in_trans_n[i]] = gen_pos[j,start_exon_n[i]] - gen_pos[j,"start_cds"]
#     gen_pos[j,end_exon_in_trans_n[i]] = gen_pos[j,end_exon_n[i]] - gen_pos[j,"start_cds"]
#   }}
# return(gen_pos)}
# 
# pos_exon_trans_coords = exons_in_trans_coords(pos_exonic)
# head(pos_exon_trans_coords)
# 
# #then we actually calculate frameness information
# 
# #example of start_cds to ATG from YAL003W
# #chr I 142173-142253
# #non coding A is 142173, coding starts at 142174 (a of atg)
# # Aatggcatccaccgatttctccaagattgaaactttgaaacaattaaacgcttctttggct
# # M  A  S  T  D  F  S  K  I  E  T  L  K  Q  L  N  A  S  L  A 
# # gacaagtcatacattgaa gg (so exits in frame 1)
# # D  K  S  Y  I  E 
# 
# pos_exon_frameness = function(gen_pos){
#   gen_pos$exon_enter_frame_1 = (gen_pos$start_exon_in_trans_1 - 1) %% 3
#     #this adjustment is because we're counting by base 0, not base 1
#   gen_pos$exon_leave_frame_1 = (gen_pos$end_exon_in_trans_1 - gen_pos$start_exon_in_trans_1) %% 3
#   
#   gen_pos$exon_enter_frame_2 = ifelse(!is.na(gen_pos$start_exon_in_trans_2), ((gen_pos$exon_leave_frame_1 + 1) %% 3), NA)
#     #the ifelse checks to make sure this next exon even exists, since the frame info here really is determined by the previous exon
#   gen_pos$exon_leave_frame_2 = (gen_pos$end_exon_in_trans_2 - gen_pos$start_exon_in_trans_2 + gen_pos$exon_enter_frame_2 ) %% 3
#   
#   gen_pos$exon_enter_frame_3 = ifelse(!is.na(gen_pos$start_exon_in_trans_3), ((gen_pos$exon_leave_frame_2 + 1) %% 3), NA)
#     #the ifelse checks to make sure this next exon even exists, since the frame info here really is determined by the previous exon
#   gen_pos$exon_leave_frame_3 = (gen_pos$end_exon_in_trans_3 - gen_pos$start_exon_in_trans_3 + gen_pos$exon_enter_frame_3 ) %% 3
#   
#   gen_pos$exon_enter_frame_4 = ifelse(!is.na(gen_pos$start_exon_in_trans_4), ((gen_pos$exon_leave_frame_3 + 1) %% 3), NA)
#   #the ifelse checks to make sure this next exon even exists, since the frame info here really is determined by the previous exon
#   gen_pos$exon_leave_frame_4 = (gen_pos$end_exon_in_trans_4 - gen_pos$start_exon_in_trans_4 + gen_pos$exon_enter_frame_4 ) %% 3
#   
#   gen_pos$exon_enter_frame_5 = ifelse(!is.na(gen_pos$start_exon_in_trans_5), ((gen_pos$exon_leave_frame_4 + 1) %% 3), NA)
#     #the ifelse checks to make sure this next exon even exists, since the frame info here really is determined by the previous exon
#   gen_pos$exon_leave_frame_5 = (gen_pos$end_exon_in_trans_5 - gen_pos$start_exon_in_trans_5 + gen_pos$exon_enter_frame_5 ) %% 3
#   
#   gen_pos$exon_enter_frame_6 = ifelse(!is.na(gen_pos$start_exon_in_trans_6), ((gen_pos$exon_leave_frame_5 + 1) %% 3), NA)
#     #the ifelse checks to make sure this next exon even exists, since the frame info here really is determined by the previous exon
#   gen_pos$exon_leave_frame_6 = (gen_pos$end_exon_in_trans_6 - gen_pos$start_exon_in_trans_6 + gen_pos$exon_enter_frame_6 ) %% 3
#   
#   gen_pos$exon_enter_frame_7 = ifelse(!is.na(gen_pos$start_exon_in_trans_7), ((gen_pos$exon_leave_frame_6 + 1) %% 3), NA)
#     #the ifelse checks to make sure this next exon even exists, since the frame info here really is determined by the previous exon
#   gen_pos$exon_leave_frame_7 = (gen_pos$end_exon_in_trans_7 - gen_pos$start_exon_in_trans_7 + gen_pos$exon_enter_frame_7 ) %% 3
#   
#   gen_pos$exon_enter_frame_8 = ifelse(!is.na(gen_pos$start_exon_in_trans_8), ((gen_pos$exon_leave_frame_7 + 1) %% 3), NA)
#     #the ifelse checks to make sure this next exon even exists, since the frame info here really is determined by the previous exon
#   gen_pos$exon_leave_frame_8 = (gen_pos$end_exon_in_trans_8 - gen_pos$start_exon_in_trans_8 + gen_pos$exon_enter_frame_8 ) %% 3
#   
#   return(gen_pos)
# }
# 
# pos_exon_frames = pos_exon_frameness(pos_exon_trans_coords)
# head(pos_exon_frames)
# 
# # #checking a handful of examples to make sure things work
# # sample_n(filter(pos_exon_frames, exon_number!=1), size = 1)
# # filter(pos_exon_frames, gene_name==c("YGR214W")) # checks out
# # filter(pos_exon_frames, gene_name==c("YML034W")) # checks out
#  #filter(pos_exon_frames, gene_name==c("YBR186W")) # checks out
# # sample_n(filter(pos_exon_frames, exon_number!=1, exon_leave_frame_1!=2), size = 1)
#  #filter(pos_exon_frames, gene_name==c("YPL081W")) # checks out
# 
# #so we have the transcript coordinates and frameness information. 
# #how do we then relate the reads to this information?
# 
# #start by putting the reads into transcript coordinates
# pos_exon_frames$start_read_in_cds = pos_exon_frames$start_read - pos_exon_frames$start_cds
# pos_exon_frames$end_read_in_cds = pos_exon_frames$end_read - pos_exon_frames$start_cds
# 
# sample_n(pos_exon_frames, size = 1)
# ?cut
# strtrim("happy", width = -1)
# table(pos_exon_frames$occurs_in_exon)
# strtoi(strsplit("start_exon_1", split = "_")[[1]][3])
