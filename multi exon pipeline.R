getwd()

library("dplyr")

reading_frame_multi_intron = function(inputfile, dataset_name){
  #import data
  unique_orfs = read.table(inputfile, header=FALSE, stringsAsFactors = FALSE)
  colnames(unique_orfs) = c("chr_read", "start_read", "end_read", "frag_count", "arbitrary_value", "strand_read",
                            "chr_gene", "start_cds", "end_cds", "gene_name", "bed_score", "strand_cds", 
                            "thick_start", "thick_end", "RGB", "exon_number", "exon_size", "exon_start", "overlap")
  head(unique_orfs)
  
  unique_orfs[,"read_length"] = unique_orfs$end_read - unique_orfs$start_read
  pdf(sprintf("read length %s.pdf", dataset_name), useDingbats = FALSE)
  hist(unique_orfs$read_length, main=sprintf("read length %s", dataset_name))
  dev.off()
  head(unique_orfs)
  
  
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
  
  #processing exon genes. 
    #separate by + and -
    neg = filter(genic, strand_read=="-")
    pos = filter(genic, strand_read=="+")
  
  #split exons from column with commas into distinct columns
  #positive
  
  start_split_pos = strsplit(pos$exon_start, ",")
  empty_start_pos=matrix(data=NA, nrow=nrow(pos), ncol=8)
  colnames(empty_start_pos) = c(1:8)
  
  for (j in 1:8){
    empty_start_pos[,j]=as.numeric(sapply(start_split_pos, function(x){x[j]}))
    colnames(empty_start_pos)[j]=c(sprintf("start_%s", j))
  }
  
  exon_size_split_pos = strsplit(pos$exon_size, ",")
  empty_exon_size_pos=matrix(data=NA, nrow=nrow(pos), ncol=8)
  colnames(empty_exon_size_pos) = c(1:8)
  
  for (j in 1:8){
    empty_exon_size_pos[,j]=as.numeric(sapply(exon_size_split_pos, function(x){x[j]}))
    colnames(empty_exon_size_pos)[j]=c(sprintf("exon_size_%s", j))
  }
  
  pos_exons = cbind(pos, empty_start_pos, empty_exon_size_pos)
  
  #negative
  
  start_split_neg = lapply(strsplit(neg$exon_start, ","), rev)
  
  empty_start_neg=matrix(data=NA, nrow=nrow(neg), ncol=8)
  colnames(empty_start_neg) = c(1:8)
  
  for (j in 1:8){
    empty_start_neg[,j]=as.numeric(sapply(start_split_neg, function(x){x[j]}))
    colnames(empty_start_neg)[j]=c(sprintf("start_%s", j))
  }
  
  exon_size_split_neg = lapply(strsplit(neg$exon_size, ","), rev)
  
  empty_exon_size_neg=matrix(data=NA, nrow=nrow(neg), ncol=8)
  colnames(empty_exon_size_neg) = c(1:8)
  
  for (j in 1:8){
    empty_exon_size_neg[,j]=as.numeric(sapply(exon_size_split_neg, function(x){x[j]}))
    colnames(empty_exon_size_neg)[j]=c(sprintf("exon_size_%s", j))
  }
  
  neg_exons = cbind(neg, empty_start_neg, empty_exon_size_neg)
  head(neg_exons)
  
    #define absolute genomic coordinates of exons
    starts=NULL
    ends=NULL
    
    for (i in 1:8){
      starts[i] = paste("start_exon_", i, sep="")
      ends[i] = paste("end_exon_", i, sep="")
    }
    
    #for both positive genes  
    for (i in 1:8){
      pos_exons[,starts[i]] = pos_exons$start_cds + pos_exons[,sprintf("start_%s", i)]  
      #DON'T DO +1 b/c it's 0 base
      pos_exons[,ends[i]] = pos_exons$start_cds + pos_exons[,sprintf("start_%s", i)] + pos_exons[,sprintf("exon_size_%s", i)] -1
      #need -1 because it's 0 base
    }
    
    #and for negative strand genes
    for (i in 1:8){
      neg_exons[,ends[i]] = neg_exons$start_cds + neg_exons[,sprintf("start_%s", i)]  
      #don't do +1 b/c its 0 based
      neg_exons[,starts[i]] = neg_exons$start_cds + neg_exons[,sprintf("start_%s", i)] + neg_exons[,sprintf("exon_size_%s", i)] -1
      #need -1 because it's 0 base
    }
    
    head(neg_exons)[1,]
  
    #function to check if exonic and if so which exon for positive strand
    positive_exon_rows = function(row_entry_positive){
      where_positive = NULL
      if (!is.na(as.numeric(row_entry_positive["start_exon_1"])) & as.numeric(row_entry_positive["start_read"]) >= as.numeric(row_entry_positive["start_exon_1"]) & as.numeric(row_entry_positive["end_read"]) <= as.numeric(row_entry_positive["end_exon_1"])) { 
        where_positive = "start_exon_1"
      } else if (!is.na(as.numeric(row_entry_positive["start_exon_2"])) & as.numeric(row_entry_positive["start_read"]) >= as.numeric(row_entry_positive["start_exon_2"]) & as.numeric(row_entry_positive["end_read"]) <= as.numeric(row_entry_positive["end_exon_2"])){
        where_positive = "start_exon_2"
      } else if (!is.na(as.numeric(row_entry_positive["start_exon_3"])) & as.numeric(row_entry_positive["start_read"]) >= as.numeric(row_entry_positive["start_exon_3"]) & as.numeric(row_entry_positive["end_read"]) <= as.numeric(row_entry_positive["end_exon_3"])){
        where_positive = "start_exon_3"
      } else if (!is.na(as.numeric(row_entry_positive["start_exon_4"])) & as.numeric(row_entry_positive["start_read"]) >= as.numeric(row_entry_positive["start_exon_4"]) & as.numeric(row_entry_positive["end_read"]) <= as.numeric(row_entry_positive["end_exon_4"])){
        where_positive = "start_exon_4"
      } else if (!is.na(as.numeric(row_entry_positive["start_exon_5"])) & as.numeric(row_entry_positive["start_read"]) >= as.numeric(row_entry_positive["start_exon_5"]) & as.numeric(row_entry_positive["end_read"]) <= as.numeric(row_entry_positive["end_exon_5"])){
        where_positive = "start_exon_5"
      } else if (!is.na(as.numeric(row_entry_positive["start_exon_6"])) & as.numeric(row_entry_positive["start_read"]) >= as.numeric(row_entry_positive["start_exon_6"]) & as.numeric(row_entry_positive["end_read"]) <= as.numeric(row_entry_positive["end_exon_6"])){
        where_positive = "start_exon_6"
      } else if (!is.na(as.numeric(row_entry_positive["start_exon_7"])) & as.numeric(row_entry_positive["start_read"]) >= as.numeric(row_entry_positive["start_exon_7"]) & as.numeric(row_entry_positive["end_read"]) <= as.numeric(row_entry_positive["end_exon_7"])){
        where_positive = "start_exon_7"
      } else if (!is.na(as.numeric(row_entry_positive["start_exon_8"])) & as.numeric(row_entry_positive["start_read"]) >= as.numeric(row_entry_positive["start_exon_8"]) & as.numeric(row_entry_positive["end_read"]) <= as.numeric(row_entry_positive["end_exon_8"])){
        where_positive = "start_exon_8"
      } else {
        where_positive = "intronic"
      } 
      where_positive }
    
    pos_exons_occurs = apply(pos_exons, 1, positive_exon_rows)
    
    pos_exons_copy = pos_exons
    pos_exons_copy$occurs_in_exon = pos_exons_occurs
  
    #function to check if exonic and if so which exon for negative strand
    negative_exon_rows = function(row_entry_negative){
      where_negative = NULL
      if (!is.na(as.numeric(row_entry_negative["start_exon_1"])) & as.numeric(row_entry_negative["end_read"]) <= as.numeric(row_entry_negative["start_exon_1"]) & as.numeric(row_entry_negative["start_read"]) >= as.numeric(row_entry_negative["end_exon_1"])){
        where_negative = "start_exon_1"
      } else if (!is.na(as.numeric(row_entry_negative["start_exon_2"])) & as.numeric(row_entry_negative["end_read"]) <= as.numeric(row_entry_negative["start_exon_2"]) & as.numeric(row_entry_negative["start_read"]) >= as.numeric(row_entry_negative["end_exon_2"])){
        where_negative = "start_exon_2"
      } else if (!is.na(as.numeric(row_entry_negative["start_exon_3"])) & as.numeric(row_entry_negative["end_read"]) <= as.numeric(row_entry_negative["start_exon_3"]) & as.numeric(row_entry_negative["start_read"]) >= as.numeric(row_entry_negative["end_exon_3"])){
        where_negative = "start_exon_3"
      } else if (!is.na(as.numeric(row_entry_negative["start_exon_4"])) & as.numeric(row_entry_negative["end_read"]) <= as.numeric(row_entry_negative["start_exon_4"]) & as.numeric(row_entry_negative["start_read"]) >= as.numeric(row_entry_negative["end_exon_4"])){
        where_negative = "start_exon_4"
      } else if (!is.na(as.numeric(row_entry_negative["start_exon_5"])) & as.numeric(row_entry_negative["end_read"]) <= as.numeric(row_entry_negative["start_exon_5"]) & as.numeric(row_entry_negative["start_read"]) >= as.numeric(row_entry_negative["end_exon_5"])){
        where_negative = "start_exon_5"
      } else if (!is.na(as.numeric(row_entry_negative["start_exon_5"])) & as.numeric(row_entry_negative["end_read"]) <= as.numeric(row_entry_negative["start_exon_6"]) & as.numeric(row_entry_negative["start_read"]) >= as.numeric(row_entry_negative["end_exon_5"])){
        where_negative = "start_exon_6"
      } else if (!is.na(as.numeric(row_entry_negative["start_exon_7"])) & as.numeric(row_entry_negative["end_read"]) <= as.numeric(row_entry_negative["start_exon_7"]) & as.numeric(row_entry_negative["start_read"]) >= as.numeric(row_entry_negative["end_exon_7"])){
        where_negative = "start_exon_7"
      } else if (!is.na(as.numeric(row_entry_negative["start_exon_8"])) & as.numeric(row_entry_negative["end_read"]) <= as.numeric(row_entry_negative["start_exon_8"]) & as.numeric(row_entry_negative["start_read"]) >= as.numeric(row_entry_negative["end_exon_8"])){
        where_negative = "start_exon_8"
      } else {
        where_negative = "intronic"
      }
      where_negative}
    
    neg_exons_occurs = apply(neg_exons, 1, negative_exon_rows)
    
    neg_exons_copy = neg_exons
    neg_exons_copy$occurs_in_exon = neg_exons_occurs
  
    ## now we need to see if the exonic fragments are in frame or not
    #pull exonic fragments
    pos_exonic = filter(pos_exons_copy, occurs_in_exon!="intronic")
    pos_intronic = filter(pos_exons_copy, occurs_in_exon=="intronic") #just saving these for later
    write.csv(pos_intronic, sprintf("intronic reads pos %s.csv", dataset_name))
    head(pos_exonic)
    
    neg_exonic = filter(neg_exons_copy, occurs_in_exon!="intronic")
    neg_intronic = filter(neg_exons_copy, occurs_in_exon=="intronic") #just saving these for later
    write.csv(neg_intronic, sprintf("intronic reads neg %s.csv", dataset_name))
    head(neg_exonic)
    
  #for positive strand genes
  #we'll do this by utilizing the idea that a reads CDS coordinates can be defined as:
  # start read CDS loc == sum(exon sizes from 1 to exon before the read occurs in) + (start read genomic location - start exon genomic location for the exon the read is in)
  # or when s is size, e is exon start site, g is start read genomic location, CDS location is c
  # there are 1 to n exons in a gene, and the read is in exon j
  # c = Σ(s[sub 1:j-1]) + (g - e[sub j])
  
  #sum previous exons
  pos_exon_summer = function(pos_summing){
    pos_prior_sums = NULL
    if (pos_summing["occurs_in_exon"]=="start_exon_1"){
      pos_prior_sums = 0
    } else if (pos_summing["occurs_in_exon"]=="start_exon_2"){
      pos_prior_sums = as.numeric(pos_summing["exon_size_1"])
    } else if (pos_summing["occurs_in_exon"]=="start_exon_3"){
      pos_prior_sums = as.numeric(pos_summing["exon_size_1"]) + as.numeric(pos_summing["exon_size_2"])
    } else if (pos_summing["occurs_in_exon"]=="start_exon_4"){
      pos_prior_sums = as.numeric(pos_summing["exon_size_1"]) + as.numeric(pos_summing["exon_size_2"]) + as.numeric(pos_summing["exon_size_3"]) 
    } else if (pos_summing["occurs_in_exon"]=="start_exon_5"){
      pos_prior_sums = as.numeric(pos_summing["exon_size_1"]) + as.numeric(pos_summing["exon_size_2"]) + as.numeric(pos_summing["exon_size_3"]) + as.numeric(pos_summing["exon_size_4"]) 
    } else if (pos_summing["occurs_in_exon"]=="start_exon_6"){
      pos_prior_sums = as.numeric(pos_summing["exon_size_1"]) + as.numeric(pos_summing["exon_size_2"]) + as.numeric(pos_summing["exon_size_3"]) + as.numeric(pos_summing["exon_size_4"]) + as.numeric(pos_summing["exon_size_5"]) 
    } else if (pos_summing["occurs_in_exon"]=="start_exon_7"){
      pos_prior_sums = as.numeric(pos_summing["exon_size_1"]) + as.numeric(pos_summing["exon_size_2"]) + as.numeric(pos_summing["exon_size_3"]) + as.numeric(pos_summing["exon_size_4"]) + as.numeric(pos_summing["exon_size_5"]) + as.numeric(pos_summing["exon_size_6"])
    } else if (pos_summing["occurs_in_exon"]=="start_exon_8"){
      pos_prior_sums = as.numeric(pos_summing["exon_size_1"]) + as.numeric(pos_summing["exon_size_2"]) + as.numeric(pos_summing["exon_size_3"]) + as.numeric(pos_summing["exon_size_4"]) + as.numeric(pos_summing["exon_size_5"]) + as.numeric(pos_summing["exon_size_6"]) + as.numeric(pos_summing["exon_size_7"])
    } else {
      pos_prior_sums = NA
    }}
  
  pos_priors = apply(pos_exonic, 1, pos_exon_summer)
  pos_exonic$prior_exon_sums = pos_priors
  
  #check if it worked
  head(pos_exonic)[1,] #checked
  filter(pos_exonic, exon_number!=1)[1,] #checked
  filter(pos_exonic, occurs_in_exon=="start_exon_3")[1,] #checked
  
  #find start and end location in current exon
  pos_math_start = function(pos_maths){
    read_start_in_exon = as.numeric(pos_maths["start_read"]) - as.numeric(pos_maths[pos_maths["occurs_in_exon"]]) #g[start] - e[sub j]
  }
  
  pos_math_end = function(pos_maths){
    read_end_in_exon = (as.numeric(pos_maths["end_read"]) -1) - as.numeric(pos_maths[pos_maths["occurs_in_exon"]]) #g[end] - e[sub j], b/c its first nt not in read
  }
  
  pos_exonic$read_start_in_exon = apply(pos_exonic, 1, pos_math_start)
  pos_exonic$read_end_in_exon = apply(pos_exonic, 1, pos_math_end)  
  
  #check to make sure they worked
  head(pos_exonic)[1,] #checked
  filter(pos_exonic, exon_number!=1)[1,] #checked
  
  #define the location in the cds
  pos_exonic$read_start_in_cds = pos_exonic$prior_exon_sums + pos_exonic$read_start_in_exon 
  pos_exonic$read_end_in_cds = pos_exonic$prior_exon_sums + pos_exonic$read_end_in_exon 
  
  #check to make sure it worked
  filter(pos_exonic, exon_number!=1)[1,] #checked
  
  pos_exonic$read_start_cds_frame = (pos_exonic$read_start_in_cds) %% 3 #because counting in 0 space
  pos_exonic$read_end_cds_frame = (pos_exonic$read_end_in_cds) %% 3 #because counting in 0 space
  
  #check to make sure it worked
  head(pos_exonic)[1,] #checked
  filter(pos_exonic, exon_number!=1)[1,] #checked
  filter(pos_exonic, read_start_cds_frame==2, read_end_cds_frame==1)[1,]
  filter(pos_exonic, read_start_cds_frame==2, read_end_cds_frame==1)[2,]
  
  table(pos_exonic$read_start_cds_frame)
  
  #for negative strand genes
  #we'll do this by utilizing the idea that a reads CDS coordinates can be defined as:
  # start read CDS loc == sum(exon sizes from 1 to exon before the read occurs in) + (start exon genomic location for the exon the read is in - end read genomic location)
  # or when s is size, e is exon start site, g is start read genomic location, CDS location is c
  # there are 1 to n exons in a gene, and the read is in exon j
  # c = Σ(s[sub 1:j-1]) + (e[sub j] - g)
  
  #summing the preceding exon sizes should be the same
  neg_exon_summer = function(neg_summing){
    neg_prior_sums = NULL
    if (neg_summing["occurs_in_exon"]=="start_exon_1"){
      neg_prior_sums = 0
    } else if (neg_summing["occurs_in_exon"]=="start_exon_2"){
      neg_prior_sums = as.numeric(neg_summing["exon_size_1"])
    } else if (neg_summing["occurs_in_exon"]=="start_exon_3"){
      neg_prior_sums = as.numeric(neg_summing["exon_size_1"]) + as.numeric(neg_summing["exon_size_2"])
    } else if (neg_summing["occurs_in_exon"]=="start_exon_4"){
      neg_prior_sums = as.numeric(neg_summing["exon_size_1"]) + as.numeric(neg_summing["exon_size_2"]) + as.numeric(neg_summing["exon_size_3"]) 
    } else if (neg_summing["occurs_in_exon"]=="start_exon_5"){
      neg_prior_sums = as.numeric(neg_summing["exon_size_1"]) + as.numeric(neg_summing["exon_size_2"]) + as.numeric(neg_summing["exon_size_3"]) + as.numeric(neg_summing["exon_size_4"]) 
    } else if (neg_summing["occurs_in_exon"]=="start_exon_6"){
      neg_prior_sums = as.numeric(neg_summing["exon_size_1"]) + as.numeric(neg_summing["exon_size_2"]) + as.numeric(neg_summing["exon_size_3"]) + as.numeric(neg_summing["exon_size_4"]) + as.numeric(neg_summing["exon_size_5"]) 
    } else if (neg_summing["occurs_in_exon"]=="start_exon_7"){
      neg_prior_sums = as.numeric(neg_summing["exon_size_1"]) + as.numeric(neg_summing["exon_size_2"]) + as.numeric(neg_summing["exon_size_3"]) + as.numeric(neg_summing["exon_size_4"]) + as.numeric(neg_summing["exon_size_5"]) + as.numeric(neg_summing["exon_size_6"])
    } else if (neg_summing["occurs_in_exon"]=="start_exon_8"){
      neg_prior_sums = as.numeric(neg_summing["exon_size_1"]) + as.numeric(neg_summing["exon_size_2"]) + as.numeric(neg_summing["exon_size_3"]) + as.numeric(neg_summing["exon_size_4"]) + as.numeric(neg_summing["exon_size_5"]) + as.numeric(neg_summing["exon_size_6"]) + as.numeric(neg_summing["exon_size_7"])
    } else {
      neg_prior_sums = NA
    }}
  
  neg_priors = apply(neg_exonic, 1, neg_exon_summer)
  neg_exonic$prior_exon_sums = neg_priors
  
  #find location in current exon, order of g and e different than positive strand
  neg_math_start = function(neg_math){
    read_start_in_exon = as.numeric(neg_math[neg_math["occurs_in_exon"]]) - (as.numeric(neg_math["end_read"]) -1 )# e[sub j] - g[start], -1 b/c 1st nt outside of read
  }
  
  neg_math_end = function(neg_math){
    neg_math["read_end_in_exon"] = as.numeric(neg_math[neg_math["occurs_in_exon"]]) - as.numeric(neg_math["start_read"]) #e[sub j] - g[end]
  }
  
  neg_exonic$read_start_in_exon = apply(neg_exonic, 1, neg_math_start)
  neg_exonic$read_end_in_exon = apply(neg_exonic, 1, neg_math_end)
  
  #define the location in the cds
  neg_exonic$read_start_in_cds = neg_exonic$prior_exon_sums + neg_exonic$read_start_in_exon 
  neg_exonic$read_end_in_cds = neg_exonic$prior_exon_sums + neg_exonic$read_end_in_exon 
  
  #translate position in exon to position in CDS
  neg_exonic$read_start_cds_frame = (neg_exonic$read_start_in_cds) %% 3 #because counting in 0 space
  neg_exonic$read_end_cds_frame = (neg_exonic$read_end_in_cds) %% 3 #because counting in 0 space
    
    #now lets try to see how the framing distribution works out (0,0 vs 0,1, vs 0,2... etc)
    head(pos_exonic)
    pos_exonic$joint_frame = paste(pos_exonic$read_start_cds_frame, pos_exonic$read_end_cds_frame, sep=",")
    table(pos_exonic$joint_frame)
    pdf(sprintf("reading frame distribution %s.pdf", dataset_name), useDingbats = FALSE)
    barplot(table(pos_exonic$joint_frame), main=sprintf("reading frame dist, +, %s", dataset_name))
      #do one weighted by frag count
      weighted_pos_read_frame = matrix(c(table(pos_exonic$joint_frame),
                                         c(sum(pos_exonic[which(pos_exonic$joint_frame=="0,0"),"frag_count"]),
                                           sum(pos_exonic[which(pos_exonic$joint_frame=="0,1"),"frag_count"]),
                                           sum(pos_exonic[which(pos_exonic$joint_frame=="0,2"),"frag_count"]),
                                           sum(pos_exonic[which(pos_exonic$joint_frame=="1,0"),"frag_count"]),
                                           sum(pos_exonic[which(pos_exonic$joint_frame=="1,1"),"frag_count"]),
                                           sum(pos_exonic[which(pos_exonic$joint_frame=="1,2"),"frag_count"]),
                                           sum(pos_exonic[which(pos_exonic$joint_frame=="2,0"),"frag_count"]),
                                           sum(pos_exonic[which(pos_exonic$joint_frame=="2,1"),"frag_count"]),
                                           sum(pos_exonic[which(pos_exonic$joint_frame=="2,2"),"frag_count"]))), ncol=2)
      barplot(weighted_pos_read_frame[,2], main=sprintf("weighted reading frame dist, +, %s", dataset_name))
              
    neg_exonic$joint_frame = paste(neg_exonic$read_start_cds_frame, neg_exonic$read_end_cds_frame, sep=",")
    table(neg_exonic$joint_frame)
    barplot(table(neg_exonic$joint_frame), main=sprintf("reading frame dist, -, %s", dataset_name))
      #do one weighted by frag count
      weighted_neg_read_frame = matrix(c(table(neg_exonic$joint_frame),
                                         c(sum(neg_exonic[which(neg_exonic$joint_frame=="0,0"),"frag_count"]),
                                           sum(neg_exonic[which(neg_exonic$joint_frame=="0,1"),"frag_count"]),
                                           sum(neg_exonic[which(neg_exonic$joint_frame=="0,2"),"frag_count"]),
                                           sum(neg_exonic[which(neg_exonic$joint_frame=="1,0"),"frag_count"]),
                                           sum(neg_exonic[which(neg_exonic$joint_frame=="1,1"),"frag_count"]),
                                           sum(neg_exonic[which(neg_exonic$joint_frame=="1,2"),"frag_count"]),
                                           sum(neg_exonic[which(neg_exonic$joint_frame=="2,0"),"frag_count"]),
                                           sum(neg_exonic[which(neg_exonic$joint_frame=="2,1"),"frag_count"]),
                                           sum(neg_exonic[which(neg_exonic$joint_frame=="2,2"),"frag_count"]))), ncol=2)
      barplot(weighted_neg_read_frame[,2], main=sprintf("weighted reading frame dist, -, %s", dataset_name))
  
    dev.off()
  
    #because of the untemplated T at the start, we get pushed into nucleotide frame 0
    #BUT there's also an untempalted T at the end, so things that end in frame 1 are actually going to 
    #be the things that translate downstream. hence 0,1 being our frame here.
    
    #reorder the tables so we can join them
    matrix(c(colnames(neg_exonic[,colnames(pos_exonic)]), colnames(pos_exonic)), ncol=2)
    
    neg_exonic_copy = neg_exonic[,colnames(pos_exonic)]
    matrix(c(colnames(neg_exonic_copy), colnames(pos_exonic)), ncol=2)
    
    tot_exonic_with_cds = rbind(pos_exonic, neg_exonic_copy)
    
    #put in the aa coordinates
    tot_exonic_with_cds$aa_start = trunc(tot_exonic_with_cds$read_start_in_cds/3)
    tot_exonic_with_cds$aa_end = trunc(tot_exonic_with_cds$read_end_in_cds/3)
    
    #see how the starts and stops look
    pdf(sprintf("aa_start_and_end %s.pdf", dataset_name), useDingbats = FALSE)
    hist(tot_exonic_with_cds$aa_start, breaks=20, main=sprintf("aa_start %s", dataset_name)) #want smaller hist bins
    hist(tot_exonic_with_cds$aa_end, breaks=20, main=sprintf("aa_end %s", dataset_name)) #want smaller hist bins
    dev.off()
    
    #but more relevant to see starts and stops in context of whole protein
    cds_length = function(exonic_with_cds){
    total_cds_length = sum(ifelse(!is.na(as.numeric(exonic_with_cds["exon_size_1"])), as.numeric(exonic_with_cds["exon_size_1"]), 0),
                           ifelse(!is.na(as.numeric(exonic_with_cds["exon_size_2"])), as.numeric(exonic_with_cds["exon_size_2"]), 0),
                           ifelse(!is.na(as.numeric(exonic_with_cds["exon_size_3"])), as.numeric(exonic_with_cds["exon_size_3"]), 0),
                           ifelse(!is.na(as.numeric(exonic_with_cds["exon_size_4"])), as.numeric(exonic_with_cds["exon_size_4"]), 0),
                           ifelse(!is.na(as.numeric(exonic_with_cds["exon_size_5"])), as.numeric(exonic_with_cds["exon_size_5"]), 0),
                           ifelse(!is.na(as.numeric(exonic_with_cds["exon_size_6"])), as.numeric(exonic_with_cds["exon_size_6"]), 0),
                           ifelse(!is.na(as.numeric(exonic_with_cds["exon_size_7"])), as.numeric(exonic_with_cds["exon_size_7"]), 0),
                           ifelse(!is.na(as.numeric(exonic_with_cds["exon_size_8"])), as.numeric(exonic_with_cds["exon_size_8"]), 0))  
    }
  
    tot_exonic_with_cds$total_cds_length = apply(tot_exonic_with_cds, 1, cds_length)
    tot_exonic_with_cds$total_cds_length_in_aa = tot_exonic_with_cds$total_cds_length/3
    tot_exonic_with_cds$dist_aa_start = tot_exonic_with_cds$aa_start/tot_exonic_with_cds$total_cds_length_in_aa
    tot_exonic_with_cds$dist_aa_end = tot_exonic_with_cds$aa_end/tot_exonic_with_cds$total_cds_length_in_aa
      
    pdf(sprintf("aa_start_and_end fract of cds %s.pdf", dataset_name), useDingbats = FALSE)
    hist(tot_exonic_with_cds$dist_aa_start, main=sprintf("aa_start fract of orf %s", dataset_name))
    hist(tot_exonic_with_cds$dist_aa_end, main=sprintf("aa_end fract of orf %s", dataset_name))
    dev.off()
    
    #out of curiousity, see length of read vs count
    pdf(sprintf("frag count v readlength %s.pdf", dataset_name), useDingbats = FALSE)
    plot(log(tot_exonic_with_cds$frag_count, base = 10), tot_exonic_with_cds$read_length,
         main=sprintf("frag count v length %s", dataset_name))
    dev.off()
    
    return(tot_exonic_with_cds)
}

#An example
# getwd()
# setwd("/Users/annamcgeachy/Google Drive/post trx reg data/datafiles_screen1_miseq/")    

#up = reading_frame_multi_intron("up_inside_orf_unique.bed", "up")

getwd()
setwd("/Users/annamcgeachy/Google Drive/post trx reg data/screen5_sorted/")

two_t_up = reading_frame_multi_intron("2tup_inside_orf_unique.bed", "two_t_up")
two_t_mid = reading_frame_multi_intron("2tmid_inside_orf_unique.bed", "two_t_mid")
two_t_down = reading_frame_multi_intron("2tdown_inside_orf_unique.bed", "two_t_down")

three_t_up = reading_frame_multi_intron("3tup_inside_orf_unique.bed", "three_t_up")
three_t_mid = reading_frame_multi_intron("3tmid_inside_orf_unique.bed", "three_t_mid")
three_t_down = reading_frame_multi_intron("3tdown_inside_orf_unique.bed", "three_t_down")

four_t_one = reading_frame_multi_intron("4t1_inside_orf_unique.bed", "four_t_one")
four_t_mid = reading_frame_multi_intron("4tmid_inside_orf_unique.bed", "four_t_mid")
four_t_two = reading_frame_multi_intron("4t2_inside_orf_unique.bed", "four_t_two")

hh = reading_frame_multi_intron("hh_inside_orf_unique.bed", "hh")
undet = reading_frame_multi_intron("undetermined_inside_orf_unique.bed", "undet")

write.csv(two_t_up, "two_t_up.csv")
write.csv(two_t_mid, "two_t_mid.csv")
write.csv(two_t_down, "two_t_down.csv")

write.csv(three_t_up, "three_t_up.csv")
write.csv(three_t_mid, "three_t_mid.csv")
write.csv(three_t_down, "three_t_down.csv")

write.csv(four_t_one, "four_t_one.csv")
write.csv(four_t_mid, "four_t_mid.csv")
write.csv(four_t_two, "four_t_two.csv")

write.csv(hh, "hh.csv")
write.csv(undet, "undet.csv")

screen5_genes = c(two_t_up$gene_name, two_t_mid$gene_name, two_t_down$gene_name,
                  three_t_up$gene_name, three_t_mid$gene_name, three_t_down$gene_name,
                  four_t_one$gene_name, four_t_mid, four_t)

#post processing

all_genes_pre = c(norecomb_bigger1$gene_name, norecomb_bigger2$gene_name, norecomb_bigger3$gene_name,
  norecomb_bigger4$gene_name, norecomb_bigger5$gene_name, norecomb_bigger6$gene_name,
  norecomb_bigger7$gene_name, norecomb_bigger8$gene_name, norecomb_bigger9$gene_name,
  norecomb_bigger10$gene_name, norecomb_bigger11$gene_name, norecomb_bigger12$gene_name,
  norecomb_bigger13$gene_name, norecomb_bigger14$gene_name)

all_genes_pre_uniqe = unique(all_genes_pre)
all_genes_pre_uniqe_data = data.frame(all_genes_pre_uniqe)

turbo2_all_uniques$time2 = turbo2_time2[match(turbo2_all_uniques$unique.turbo2_all_uniques, turbo2_time2$unique), "frag_count"]
all_genes_pre_uniqe_data$table1 = table1[match(all_genes_pre_uniqe_data$all_genes_pre_uniqe, table1$Var1), "Freq"]
all_genes_pre_uniqe_data$table2 = table2[match(all_genes_pre_uniqe_data$all_genes_pre_uniqe, table2$Var1), "Freq"]
all_genes_pre_uniqe_data$table3 = 0
all_genes_pre_uniqe_data$table3 = table3[match(all_genes_pre_uniqe_data$all_genes_pre_uniqe, table3$Var1), "Freq"]

ifelse(head(all_genes_pre_uniqe_data)[,2]  
all_genes_pre_uniqe_data[1,1] + all_genes_pre_uniqe_data[1,2]
is.na(all_genes_pre_uniqe_data$table2) = 0 

head(all_genes_pre_uniqe_data)

megea_norecomb = rbind(norecomb_bigger1, norecomb_bigger2, norecomb_bigger3, norecomb_bigger4,
      norecomb_bigger5, norecomb_bigger6, norecomb_bigger7, norecomb_bigger8,
      norecomb_bigger9, norecomb_bigger10, norecomb_bigger11, norecomb_bigger12,
      norecomb_bigger13, norecomb_bigger14)

mega_norecomb_genecount = data.frame(table(megea_norecomb$gene_name))
head(mega_norecomb_genecount)
mega_norecomb_genecount_test = mega_norecomb_genecount[order(mega_norecomb_genecount$Freq),]

head(mega_norecomb_genecount_test)
write.csv(mega_norecomb_genecount_test, "mega_norecomb_genecount_test.csv")

write.csv(mega_norecomb_genecount_test[which(mega_norecomb_genecount_test$Freq==1),], "one_counts_pre.txt")
max(mega_norecomb_genecount_test$Freq)
tail(mega_norecomb_genecount_test)

plot(mega_norecomb_genecount_test$Freq)
mega_cdf = ecdf(mega_norecomb_genecount_test$Freq)

?ecdf
pdf("pre-recomb_cdf_stretch.pdf", height = 3, width = 5, useDingbats = FALSE)
plot(mega_cdf, xlab="Fragments per Gene", main="CDF of pre-recombination coverage by gene", do.points = FALSE, lwd=5)
abline(v=median(mega_norecomb_genecount_test$Freq))
# axis(1, at=median(mega_norecomb_genecount_test$Freq), 
#      labels = as.character(median(mega_norecomb_genecount_test$Freq)), 
#      col.axis="#62866a")
abline(h=.5)
dev.off()
lines(mega_cdf)
nrow(mega_norecomb_genecount_test)
getwd()

post_recomb = read.csv("screen1-hi-unsort.csv")
head(post_recomb)
post_gene_counts = data.frame(table(post_recomb$gene_name))
post_gene_counts = post_gene_counts[order(post_gene_counts$Freq),]

head(post_gene_counts)
unsort_cdf = ecdf(post_gene_counts$Freq)

median(post_gene_counts$Freq)
pdf("post-recomb_cdf_stretch.pdf", height = 3, width = 5, useDingbats = FALSE)
plot(unsort_cdf, xlab="Fragments per Gene", main="CDF of post-recombination coverage by gene", do.points = FALSE, lwd = 5)
abline(v=median(post_gene_counts$Freq))
# axis(1, at=median(post_gene_counts$Freq), labels = as.character(median(post_gene_counts$Freq)))
abline(h=.5)
dev.off()
lines(unsort_cdf)


head(post_gene_counts)
tail(post_gene_counts)

# 
# write.csv(unsorted_2, "screen2-hi-unsort.csv")
# write.csv(up_2, "screen2-hi-up.csv")
# write.csv(mid_2, "screen2-hi-mid.csv")
# write.csv(down_2, "screen2-hi-down.csv")
# write.csv(norecomb_2, "screen2-hi-norecomb.csv")
# 
# write.csv(unsorted_1, "screen1-hi-unsort.csv")
# write.csv(up_1, "screen1-hi-up.csv")
# write.csv(mid_1, "screen1-hi-mid.csv")
# write.csv(down_1, "screen1-hi-down.csv")
# write.csv(norecomb_1, "screen1-hi-norecomb.csv")

setwd("Google Drive/post trx reg data/miseq_screen5_20160504/")
turbo2_time2 = reading_frame_multi_intron("2t_2_inside_orf_unique.bed", "2t_2")
turbo2_time3 = reading_frame_multi_intron("2t_3_inside_orf_unique.bed", "2t_3")
turbo2_time4 = reading_frame_multi_intron("2t_4_inside_orf_unique.bed", "2t_4")
turbo2_time5 = reading_frame_multi_intron("2t_5_inside_orf_unique.bed", "2t_5")
turbo2_time6 = reading_frame_multi_intron("2t_6_inside_orf_unique.bed", "2t_6")

turbo3_time2 = reading_frame_multi_intron("3t_2_inside_orf_unique.bed", "3t_2")
turbo3_time3 = reading_frame_multi_intron("3t_3_inside_orf_unique.bed", "3t_3")
turbo3_time4 = reading_frame_multi_intron("3t_4_inside_orf_unique.bed", "3t_4")
turbo3_time5 = reading_frame_multi_intron("3t_5_inside_orf_unique.bed", "3t_5")
turbo3_time6 = reading_frame_multi_intron("3t_6_inside_orf_unique.bed", "3t_6")

turbo4_time2 = reading_frame_multi_intron("4t_2_inside_orf_unique.bed", "4t_2")
turbo4_time3 = reading_frame_multi_intron("4t_3_inside_orf_unique.bed", "4t_3")
turbo4_time4 = reading_frame_multi_intron("4t_4_inside_orf_unique.bed", "4t_4")
turbo4_time5 = reading_frame_multi_intron("4t_5_inside_orf_unique.bed", "4t_5")
turbo4_time6 = reading_frame_multi_intron("4t_6_inside_orf_unique.bed", "4t_6")

turbo5_time2 = reading_frame_multi_intron("5t_2_inside_orf_unique.bed", "5t_2")
turbo5_time3 = reading_frame_multi_intron("5t_3_inside_orf_unique.bed", "5t_3")
turbo5_time4 = reading_frame_multi_intron("5t_4_inside_orf_unique.bed", "5t_4")
turbo5_time5 = reading_frame_multi_intron("5t_5_inside_orf_unique.bed", "5t_5")
turbo5_time6 = reading_frame_multi_intron("5t_6_inside_orf_unique.bed", "5t_6")

screen5_prerecomb = reading_frame_multi_intron("pre_recomb_inside_orf_unique.bed", "screen5_prerecomb")
screen5_undetermined = reading_frame_multi_intron("undertermined_inside_orf_unique.bed", "screen5_undetermined")

write.csv(turbo2_time2, "screen5_2t_2.csv")
write.csv(turbo2_time3, "screen5_2t_3.csv")
write.csv(turbo2_time4, "screen5_2t_4.csv")
write.csv(turbo2_time5, "screen5_2t_5.csv")
write.csv(turbo2_time6, "screen5_2t_6.csv")

write.csv(turbo3_time2, "screen5_3t_2.csv")
write.csv(turbo3_time3, "screen5_3t_3.csv")
write.csv(turbo3_time4, "screen5_3t_4.csv")
write.csv(turbo3_time5, "screen5_3t_5.csv")
write.csv(turbo3_time6, "screen5_3t_6.csv")

write.csv(turbo4_time2, "screen5_4t_2.csv")
write.csv(turbo4_time3, "screen5_4t_3.csv")
write.csv(turbo4_time4, "screen5_4t_4.csv")
write.csv(turbo4_time5, "screen5_4t_5.csv")
write.csv(turbo4_time6, "screen5_4t_6.csv")

write.csv(turbo5_time2, "screen5_5t_2.csv")
write.csv(turbo5_time3, "screen5_5t_3.csv")
write.csv(turbo5_time4, "screen5_5t_4.csv")
write.csv(turbo5_time5, "screen5_5t_5.csv")
write.csv(turbo5_time6, "screen5_5t_6.csv")

write.csv(screen5_prerecomb, "screen5_prerecomb.csv")
write.csv(screen5_undetermined, "screen5_undetermined.csv")

