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

## BUT exons can enter or leave in ANY frame, so we need to figure out what the frame is of each exon (UGH)
## start by defining exons in transcript coordinates

#define columns that will exist
start_exon_n = NULL
start_exon_in_cds_n = NULL
end_exon_n = NULL
end_exon_in_cds_n = NULL
for (i in 1:8){
  start_exon_n[i] = sprintf("start_exon_%s", i)
  start_exon_in_cds_n[i] = sprintf("start_exon_in_cds_%s", i)
  end_exon_n[i] = sprintf("end_exon_%s", i)
  end_exon_in_cds_n[i] = sprintf("end_exon_in_cds_%s", i)
}

#for positive exons, it goes like this
exons_in_cds_coords = function(gen_pos){
  for (i in 1:8){
    for (j in 1:nrow(gen_pos)){
    gen_pos[j,start_exon_in_cds_n[i]] = gen_pos[j,start_exon_n[i]] - gen_pos[j,"start_cds"]
    gen_pos[j,end_exon_in_cds_n[i]] = gen_pos[j,end_exon_n[i]] - gen_pos[j,"start_cds"]
  }}
return(gen_pos)}

pos_exon_cds_coords = exons_in_cds_coords(pos_exonic)
head(pos_exon_cds_coords)


lil_pos$exon_enter_frame_1 = (lil_pos$start_exon_in_cds_1 - 1) %% 3
#this adjustment is because we're counting by base 0, not base 1
lil_pos$exon_leave_frame_1 = (lil_pos$end_exon_in_cds_1 - lil_pos$start_exon_in_cds_1) %% 3
lil_pos

(570) %% 3
(570 - 1 ) %% 3

#example of start_cds to ATG from YAL003W
#chr I 142173-142253
#non coding A is 142173, coding starts at 142174 (a of atg)
# Aatggcatccaccgatttctccaagattgaaactttgaaacaattaaacgcttctttggct
# M  A  S  T  D  F  S  K  I  E  T  L  K  Q  L  N  A  S  L  A 
# gacaagtcatacattgaa gg (so exits in frame 1)
# D  K  S  Y  I  E 

lil_pos$exon_enter_frame_2 = ifelse(!is.na(lil_pos$start_exon_in_cds_2), ((lil_pos$exon_enter_frame_1 + 1) %% 3), NA)
#the ifelse checks to make sure this next exon even exists, since the frame info here really is determined by the previous exon
lil_pos$exon_leave_frame_2 = (lil_pos$end_exon_in_cds_2 - lil_pos$start_exon_in_cds_2 + lil_pos$exon_enter_frame_2 + 1) %% 3
987-447


lil_pos
>sacCer3_dna range=chrI:142620-143160 5'pad=0 3'pad=0 strand=+ repeatMasking=none
(143160 - 142620 + 1 + lil_pos[2,"exon_enter_frame_2"]) %% 3 #where + 1 is it coming out of exon 1 in frame 1 
TACTGCTGTTTCTCAAGCTGACGTCACTGTCTTCAAGGCTTTCCAATCT
GCTTACCCAGAATTCTCCAGATGGTTCAACCACATCGCTTCCAAGGCCGA
TGAATTCGACTCTTTCCCAGCTGCCTCTGCTGCCGCTGCCGAAGAAGAAG
AAGATGACGATGTCGATTTATTCGGTTCCGACGATGAAGAAGCTGACGCT
GAAGCTGAAAAGTTGAAGGCTGAAAGAATTGCCGCATACAACGCTAAGAA
GGCTGCTAAGCCAGCTAAGCCAGCTGCTAAGTCCATTGTCACTCTAGATG
TCAAGCCATGGGATGATGAAACCAATTTGGAAGAAATGGTTGCTAACGTC
AAGGCCATCGAAATGGAAGGTTTGACCTGGGGTGCTCACCAATTTATCCC
AATTGGTTTCGGTATCAAGAAGTTGCAAATTAACTGTGTTGTCGAAGATG
ACAAGGTTTCCTTGGATGACTTGCAACAAAGCATTGAAGAAGACGAAGAC
CACGTCCAATCTACCGATATTGCTGCTATGCAAAAATTATAA


lil_pos[2,]

(987 - 447) %%3
#lets start by putting everything back from genomic coordinates to transcript coordinates
lil_pos$start_read_in_cds = lil_pos$start_read - lil_pos$start_cds
lil_pos$end_read_in_cds = lil_pos$end_read - lil_pos$start_cds


lil_pos
0 %% 3
1 %% 3 
ATGGCATCCACCGATTTCTCCAAGATTGAAACTTTGAAACAATTAAACGC
TTCTTTGGCTGACAAGTCATACATTGAAGGgtatgttccgatttagttta
ctttatagatcgttgtttttctttcttttttttttttcctatggttacat
gtaaagggaagttaactaataatgattactttttttcgcttatgtgaatg
atgaatttaattctttggtccgtgtttatgatgggaagtaagacccccga
tatgagtgacaaaagagatgtggttgactatcacagtatctgacgatagc
acagagcagagtatcattattagttatctgttatttttttttcctttttt
gttcaaaaaaagaaagacagagtctaaagattgcattacaagaaaaaagt
tctcattactaacaagcaaaatgttttgtttctccttttaaaatagTACT
GCTGTTTCTCAAGCTGACGTCACTGTCTTCAAGGCTTTCCAATCTGCTTA
CCCAGAATTCTCCAGATGGTTCAACCACATCGCTTCCAAGGCCGATGAAT
TCGACTCTTTCCCAGCTGCCTCTGCTGCCGCTGCCGAAGAAGAAGAAGAT
GACGATGTCGATTTATTCGGTTCCGACGATGAAGAAGCTGACGCTGAAGC
TGAAAAGTTGAAGGCTGAAAGAATTGCCGCATACAACGCTAAGAAGGCTG
CTAAGCCAGCTAAGCCAGCTGCTAAGTCCATTGTCACTCTAGATGTCAAG
CCATGGGATGATGAAACCAATTTGGAAGAAATGGTTGCTAACGTCAAGGC
CATCGAAATGGAAGGTTTGACCTGGGGTGCTCACCAATTTATCCCAATTG
GTTTCGGTATCAAGAAGTTGCAAATTAACTGTGTTGTCGAAGATGACAAG
GTTTCCTTGGATGACTTGCAACAAAGCATTGAAGAAGACGAAGACCACGT
CCAATCTACCGATATTGCTGCTATGCAAAAATTATAA



baby = lil_pos[2,]
baby$start_1 %% 3
baby$end_
baby$start_2 %% 3
#########
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
