getwd()
setwd("~/Google Drive/post trx reg data/datafiles_screen1_and_2_hiseq_redone_in_R/")
read.csv("")
screen1_up = read.csv("screen1-hi-up.csv", stringsAsFactors=FALSE)
screen1_mid = read.csv("screen1-hi-mid.csv", stringsAsFactors=FALSE)
screen1_down = read.csv("screen1-hi-down.csv", stringsAsFactors=FALSE)

screen2_up = read.csv("screen2-hi-up.csv", stringsAsFactors=FALSE)
screen2_mid = read.csv("screen2-hi-mid.csv", stringsAsFactors=FALSE)
screen2_down = read.csv("screen2-hi-down.csv", stringsAsFactors=FALSE)

#make unique identifiers for each fragment
pasting_together = function(x){
  paste(x["chr_read"], as.numeric(x["start_read"]), as.numeric(x["end_read"]), x["strand_read"], sep="_")
}

#example of usage
#general_a_uniques = apply(general_a, 1, pasting_together)

screen1_up$uniques = apply(screen1_up, 1, pasting_together)
screen1_mid$uniques = apply(screen1_mid, 1, pasting_together)
screen1_down$uniques = apply(screen1_down, 1, pasting_together)

screen2_up$uniques = apply(screen2_up, 1, pasting_together)
screen2_mid$uniques = apply(screen2_mid, 1, pasting_together)
screen2_down$uniques = apply(screen2_down, 1, pasting_together)

#make a v b plots

library("dplyr")

xref = read.delim("../SGD_features.tab", header=FALSE, quote="")
setwd("~/Google Drive/post trx reg data/datafiles_screen1_and_2_hiseq_redone_in_R/analysis_20160926/")

#make unique identifiers for each fragment
pasting_together = function(x){
  paste(x["chr_read"], as.numeric(x["start_read"]), as.numeric(x["end_read"]), x["strand_read"], sep="_")
}

a_v_b_plots = function(location_of_file_a, condition_a, location_of_file_b, condition_b, comparison){
  ##set up A
  general_a = read.csv(location_of_file_a)
  
  general_a_uniques = apply(general_a, 1, pasting_together)
  general_a$uniques = general_a_uniques
  
  #add in the annotation
  general_a$gene_useful = xref[match(general_a$gene_name, xref$V4),"V5"]
  general_a$gene_desc = xref[match(general_a$gene_name, xref$V4),"V16"]
  
  #pull in frame fragments
  general_a_01 = filter(general_a, joint_frame=="0,1")
  
  ##set up B
  general_b = read.csv(location_of_file_b)
  
  general_b_uniques = apply(general_b, 1, pasting_together)
  general_b$uniques = general_b_uniques
  
  #add in the annotation
  general_b$gene_useful = xref[match(general_b$gene_name, xref$V4),"V5"]
  general_b$gene_desc = xref[match(general_b$gene_name, xref$V4),"V16"]
  
  #pull in frame fragments
  general_b_01 = filter(general_b, joint_frame=="0,1")
  
  head(general_a_01)
  head(general_b_01)
  
  #make joint file of all uniques between libraries
  a_and_b = c(as.character(general_a_01$uniques), as.character(general_b_01$uniques))
  length(a_and_b)
  head(a_and_b)
  a_and_b = unique(a_and_b)
  a_and_b = as.data.frame(a_and_b)
  nrow(a_and_b)
  
  #pull in fragment counts in each library
  head(a_and_b)
  a_and_b[,condition_a] = general_a_01[match(a_and_b$a_and_b, general_a_01$uniques), "frag_count"]
  table(a_and_b[,condition_a])
  a_and_b[,condition_b] = general_b_01[match(a_and_b$a_and_b, general_b_01$uniques), "frag_count"]
  table(a_and_b[,condition_b])
  
  a_and_b$a_adj = ifelse(is.na(a_and_b[,condition_a]), 0.99, a_and_b[,condition_a])
  a_and_b$b_adj = ifelse(is.na(a_and_b[,condition_b]), 0.99, a_and_b[,condition_b])
  
  #plot a v b
  png(sprintf("%s.png", comparison))
  plot(log(a_and_b$a_adj, 2), log(a_and_b$b_adj, 2), pch=16, col="#808D8D4A", 
       #        xlim=max(max(log(a_and_b[,condition_a], 2)), max(log(a_and_b[,condition_b], 2))),
       #        ylim=max(max(log(a_and_b[,condition_a], 2)), max(log(a_and_b[,condition_b], 2))),
       xlab=sprintf("log2(%s)", condition_a),
       ylab=sprintf("log2(%s)", condition_b),
       main=sprintf("%s", comparison))
  abline(a=0,b=1)
  abline(a=-4.6,b=1)
  abline(a=4.6,b=1)
  abline(a=-2.3,b=1)
  abline(a=2.3,b=1)
  dev.off()
  
  write.csv(a_and_b, sprintf("%s.csv", comparison))
}

#a_v_b_plots = function(location_of_file_a, condition_a, 
#               location_of_file_b, condition_b, 
#               comparison)

# a_v_b_plots("screen3-4-down.csv", "screen3-4-down", 
#             "screen3-4-up.csv", "screen3-4-up", 
#             "screen3-4-down-v-up")

#screen 1
a_v_b_plots("screen1-hi-up.csv", "screen1_up",
            "screen1-hi-down.csv", "screen1_down",
            "screen1_up_V_down")

a_v_b_plots("screen1-hi-up.csv", "screen1_up",
            "screen1-hi-mid.csv", "screen1_mid",
            "screen1_up_V_mid")

a_v_b_plots("screen1-hi-down.csv", "screen1_down",
            "screen1-hi-mid.csv", "screen1_mid",
            "screen1_down_V_mid")

#screen 2
a_v_b_plots("screen2-hi-up.csv", "screen2_up",
            "screen2-hi-down.csv", "screen2_down",
            "screen2_up_V_down")

a_v_b_plots("screen2-hi-up.csv", "screen2_up",
            "screen2-hi-mid.csv", "screen2_mid",
            "screen2_up_V_mid")

a_v_b_plots("screen2-hi-down.csv", "screen2_down",
            "screen2-hi-mid.csv", "screen2_mid",
            "screen2_down_V_mid")

##Look at what genes have most enrichment

#screen 1
screen1_up_V_down = read.csv("screen1_up_V_down.csv", stringsAsFactors=FALSE)

screen1_up_V_down$gene_name_up = screen1_up[match(screen1_up_V_down$a_and_b, screen1_up$uniques),"gene_name"]
screen1_up_V_down$gene_name_down = screen1_down[match(screen1_up_V_down$a_and_b, screen1_down$uniques),"gene_name"]
screen1_up_V_down$gene_name = 
  ifelse(is.na(screen1_up_V_down$gene_name_up), screen1_up_V_down$gene_name_down, 
         screen1_up_V_down$gene_name_up)


head(screen1_up_V_down)
blurgh = as.data.frame(table(screen1_up_V_down$gene_name))
head(screen1_up_V_down)

screen1_up_V_down$up_10x = ifelse((screen1_up_V_down$a_adj/screen1_up_V_down$b_adj) > 10, 
                                screen1_up_V_down$gene_name, NA)
screen1_up_V_down$up_100x = ifelse((screen1_up_V_down$a_adj/screen1_up_V_down$b_adj) > 100, 
                                 screen1_up_V_down$gene_name, NA)

screen1_up_V_down$down_10x = ifelse((screen1_up_V_down$b_adj/screen1_up_V_down$a_adj) > 10, 
                                  screen1_up_V_down$gene_name, NA)
screen1_up_V_down$down_100x = ifelse((screen1_up_V_down$b_adj/screen1_up_V_down$a_adj) > 100, 
                                   screen1_up_V_down$gene_name, NA)

pdf("screen1 gene counts hists.pdf", useDingbats = FALSE)
hist(table(screen1_up_V_down$gene_name))
hist(table(screen1_up_V_down$down_10x))
hist(table(screen1_up_V_down$down_100x))
hist(table(screen1_up_V_down$up_10x))
hist(table(screen1_up_V_down$up_100x))
dev.off()

screen1_down_10x = as.data.frame(table(screen1_up_V_down$down_10x))
screen1_down_10x$useful = xref[match(screen1_down_10x$Var1, xref$V4), "V5"]
screen1_down_10x$desc = xref[match(screen1_down_10x$Var1, xref$V4), "V16"]
head(screen1_down_10x)
screen1_down_10x = screen1_down_10x[order(screen1_down_10x$Freq, decreasing = TRUE),]

write.csv(screen1_down_10x, "screen1_10x_down_genes.csv")

screen1_up_10x = as.data.frame(table(screen1_up_V_down$up_10x))
screen1_up_10x$useful = xref[match(screen1_up_10x$Var1, xref$V4), "V5"]
screen1_up_10x$desc = xref[match(screen1_up_10x$Var1, xref$V4), "V16"]
head(screen1_up_10x)
screen1_up_10x = screen1_up_10x[order(screen1_up_10x$Freq, decreasing = TRUE),]

write.csv(screen1_up_10x, "screen1_10x_up_genes.csv")

#screen 2
screen2_up_V_down = read.csv("screen2_up_V_down.csv", stringsAsFactors=FALSE)

screen2_up_V_down$gene_name_up = screen2_up[match(screen2_up_V_down$a_and_b, screen2_up$uniques),"gene_name"]
screen2_up_V_down$gene_name_down = screen2_down[match(screen2_up_V_down$a_and_b, screen2_down$uniques),"gene_name"]
screen2_up_V_down$gene_name = 
  ifelse(is.na(screen2_up_V_down$gene_name_up), screen2_up_V_down$gene_name_down, 
         screen2_up_V_down$gene_name_up)


head(screen2_up_V_down)
blurgh = as.data.frame(table(screen2_up_V_down$gene_name))
head(screen2_up_V_down)

screen2_up_V_down$up_10x = ifelse((screen2_up_V_down$a_adj/screen2_up_V_down$b_adj) > 10, 
                                  screen2_up_V_down$gene_name, NA)
screen2_up_V_down$up_100x = ifelse((screen2_up_V_down$a_adj/screen2_up_V_down$b_adj) > 100, 
                                   screen2_up_V_down$gene_name, NA)

screen2_up_V_down$down_10x = ifelse((screen2_up_V_down$b_adj/screen2_up_V_down$a_adj) > 10, 
                                    screen2_up_V_down$gene_name, NA)
screen2_up_V_down$down_100x = ifelse((screen2_up_V_down$b_adj/screen2_up_V_down$a_adj) > 100, 
                                     screen2_up_V_down$gene_name, NA)

pdf("screen2 gene counts hists.pdf", useDingbats = FALSE)
hist(table(screen2_up_V_down$gene_name))
hist(table(screen2_up_V_down$down_10x))
hist(table(screen2_up_V_down$down_100x))
hist(table(screen2_up_V_down$up_10x))
hist(table(screen2_up_V_down$up_100x))
dev.off()

screen2_down_10x = as.data.frame(table(screen2_up_V_down$down_10x))
screen2_down_10x$useful = xref[match(screen2_down_10x$Var1, xref$V4), "V5"]
screen2_down_10x$desc = xref[match(screen2_down_10x$Var1, xref$V4), "V16"]
head(screen2_down_10x)
screen2_down_10x = screen2_down_10x[order(screen2_down_10x$Freq, decreasing = TRUE),]

write.csv(screen2_down_10x, "screen2_10x_down_genes.csv")

screen2_up_10x = as.data.frame(table(screen2_up_V_down$up_10x))
screen2_up_10x$useful = xref[match(screen2_up_10x$Var1, xref$V4), "V5"]
screen2_up_10x$desc = xref[match(screen2_up_10x$Var1, xref$V4), "V16"]
head(screen2_up_10x)
screen2_up_10x = screen2_up_10x[order(screen2_up_10x$Freq, decreasing = TRUE),]

write.csv(screen2_up_10x, "screen2_10x_up_genes.csv")

##analysis talked about with Nick (sort out high abundance fragments, MA plots, etc)

#combine all the sorts into one turbido mega table
screen1_up$source_set = "up"
screen1_mid$source_set = "mid"
screen1_down$source_set = "down"

screen1_total = rbind(screen1_up, screen1_mid, screen1_down)
nrow(screen1_total)
head(screen1_total)

#make a histogram of all counts in frame and out of frame
pdf("hist_frag_counts_two_t.pdf", useDingbats = FALSE)
hist(log(filter(screen1_total, joint_frame=="0,1")[,"frag_count"], base=2), col="#808D8D4A",
     main="screen5_2t", xlab="log(frag_count)") 
hist(log(filter(screen1_total, joint_frame!="0,1")[,"frag_count"], base=2), col="#90ffff66", add=TRUE)
#hist(log(screen1_total[,"frag_count"], base=2), col="#90ff8866", add=TRUE)
text(x=14,y=1800,"in frame", col="#998D8D4A")
text(x=14,y=1700,"out of frame", col="#99ffff66")
text(x=14,y=1700,"out of frame", col="#99ffff66")
dev.off()

head(screen1_total)

#make a list of the enriched fragments in frame (only unique names)
screen1_high_frag = as.data.frame(unique(filter(screen1_total, joint_frame=="0,1" & frag_count>(2^6))[,"uniques"]))
colnames(screen1_high_frag) = "uniques"
head(screen1_high_frag)

#add in counts
screen1_high_frag$up = screen1_up[match(screen1_high_frag$uniques, screen1_up$uniques),"frag_count"]
screen1_high_frag$mid = screen1_mid[match(screen1_high_frag$uniques, screen1_mid$uniques),"frag_count"]
screen1_high_frag$down = screen1_down[match(screen1_high_frag$uniques, screen1_down$uniques),"frag_count"]

table(is.na(screen1_high_frag$up))
table(is.na(screen1_high_frag$mid))
table(is.na(screen1_high_frag$down))

filter(screen1_high_frag, is.na(up))
filter(screen1_high_frag, is.na(mid))
filter(screen1_high_frag, is.na(down))

screen1_high_frag$up = ifelse(is.na(screen1_high_frag$up), 0.1, screen1_high_frag$up)
screen1_high_frag$mid = ifelse(is.na(screen1_high_frag$mid), 0.1, screen1_high_frag$mid)
screen1_high_frag$down = ifelse(is.na(screen1_high_frag$down), 0.1, screen1_high_frag$down)


#MA plots

pdf("screen1_MA_plots.pdf", useDingbats = FALSE)
plot(y=log(screen1_high_frag$up/screen1_high_frag$mid,2), 
     x=log(((screen1_high_frag$up +  screen1_high_frag$mid)/2),2), log="xy",
     main="screen1_MA_up_v_mid", ylab="log2(up/mid)", xlab="avg(up,mid)")
abline(h=1)

plot(y=log(screen1_high_frag$down/screen1_high_frag$mid,2), 
     x=log(((screen1_high_frag$down +  screen1_high_frag$mid)/2),2), log="xy",
     main="screen1_MA_down_v_mid", ylab="log2(down/mid)", xlab="avg(down,mid)")
abline(h=1)
dev.off()

#up/mid ratio v down/mid ratio
pdf("screen1_mid_ratios_plot.pdf", useDingbats = FALSE)
plot(y=log(screen1_high_frag$up/screen1_high_frag$mid,2),
     x=log(screen1_high_frag$down/screen1_high_frag$mid,2),
     pch=16, col="#808D8D4A",
     ylab="log2(up_ratio/down_ratio)", xlab="log2(down_ratio/up_ratio)",
     main="screen1 up/mid v down/mid") 
abline(a=0,b=1)
abline(a=-3.3,b=1)
abline(a=3.3,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
dev.off()

#pull the things that are >=2X enriched

screen1_high_frag$up_div_mid = screen1_high_frag$up/screen1_high_frag$mid
screen1_high_frag$down_div_mid = screen1_high_frag$down/screen1_high_frag$mid

screen1_high_frag$up_ratio_div_down_ratio = screen1_high_frag$up_div_mid/screen1_high_frag$down_div_mid
screen1_high_frag$down_ratio_div_up_ratio = screen1_high_frag$down_div_mid/screen1_high_frag$up_div_mid

screen1_high_frag$gene_up = screen1_up[match(screen1_high_frag$uniques, screen1_up$uniques),"gene_name"]
screen1_high_frag$gene_down = screen1_down[match(screen1_high_frag$uniques, screen1_down$uniques),"gene_name"]
screen1_high_frag$gene_mid = screen1_mid[match(screen1_high_frag$uniques, screen1_mid$uniques),"gene_name"]

screen1_high_frag$gene_name =
  ifelse(is.na(screen1_high_frag$gene_up), 
         ifelse(is.na(screen1_high_frag$gene_down), screen1_high_frag$gene_mid, screen1_high_frag$gene_down), 
         screen1_high_frag$gene_up)
sum(is.na(screen1_high_frag$gene_name))

screen1_high_frag$gene_useful = xref[match(screen1_high_frag$gene_name, xref$V4),"V5"]
screen1_high_frag$gene_desc = xref[match(screen1_high_frag$gene_name, xref$V4),"V16"]

head(screen1_high_frag)

write.csv(filter(screen1_high_frag, down_ratio_div_up_ratio>=2), "screen1_down_ratio_enriched.csv")
write.csv(filter(screen1_high_frag, up_ratio_div_down_ratio>=2), "screen1_up_ratio_enriched.csv")

#same analysis for screen 2
#combine all the sorts into one turbido mega table
screen2_up$source_set = "up"
screen2_mid$source_set = "mid"
screen2_down$source_set = "down"

screen2_total = rbind(screen2_up, screen2_mid, screen2_down)
nrow(screen2_total)
head(screen2_total)

#make a histogram of all counts in frame and out of frame
pdf("hist_frag_counts_two_t.pdf", useDingbats = FALSE)
hist(log(filter(screen2_total, joint_frame=="0,1")[,"frag_count"], base=2), col="#808D8D4A",
     main="screen5_2t", xlab="log(frag_count)") 
hist(log(filter(screen2_total, joint_frame!="0,1")[,"frag_count"], base=2), col="#90ffff66", add=TRUE)
#hist(log(screen2_total[,"frag_count"], base=2), col="#90ff8866", add=TRUE)
text(x=14,y=1800,"in frame", col="#998D8D4A")
text(x=14,y=1700,"out of frame", col="#99ffff66")
text(x=14,y=1700,"out of frame", col="#99ffff66")
dev.off()

head(screen2_total)

#make a list of the enriched fragments in frame (only unique names)
screen2_high_frag = as.data.frame(unique(filter(screen2_total, joint_frame=="0,1" & frag_count>(2^6))[,"uniques"]))
colnames(screen2_high_frag) = "uniques"
head(screen2_high_frag)

#add in counts
screen2_high_frag$up = screen2_up[match(screen2_high_frag$uniques, screen2_up$uniques),"frag_count"]
screen2_high_frag$mid = screen2_mid[match(screen2_high_frag$uniques, screen2_mid$uniques),"frag_count"]
screen2_high_frag$down = screen2_down[match(screen2_high_frag$uniques, screen2_down$uniques),"frag_count"]

table(is.na(screen2_high_frag$up))
table(is.na(screen2_high_frag$mid))
table(is.na(screen2_high_frag$down))

filter(screen2_high_frag, is.na(up))
filter(screen2_high_frag, is.na(mid))
filter(screen2_high_frag, is.na(down))

screen2_high_frag$up = ifelse(is.na(screen2_high_frag$up), 0.1, screen2_high_frag$up)
screen2_high_frag$mid = ifelse(is.na(screen2_high_frag$mid), 0.1, screen2_high_frag$mid)
screen2_high_frag$down = ifelse(is.na(screen2_high_frag$down), 0.1, screen2_high_frag$down)


#MA plots

pdf("screen2_MA_plots.pdf", useDingbats = FALSE)
plot(y=log(screen2_high_frag$up/screen2_high_frag$mid,2), 
     x=log(((screen2_high_frag$up +  screen2_high_frag$mid)/2),2), log="xy",
     main="screen2_MA_up_v_mid", ylab="log2(up/mid)", xlab="avg(up,mid)")
abline(h=1)

plot(y=log(screen2_high_frag$down/screen2_high_frag$mid,2), 
     x=log(((screen2_high_frag$down +  screen2_high_frag$mid)/2),2), log="xy",
     main="screen2_MA_down_v_mid", ylab="log2(down/mid)", xlab="avg(down,mid)")
abline(h=1)
dev.off()

#up/mid ratio v down/mid ratio
pdf("screen2_mid_ratios_plot.pdf", useDingbats = FALSE)
plot(y=log(screen2_high_frag$up/screen2_high_frag$mid,2),
     x=log(screen2_high_frag$down/screen2_high_frag$mid,2),
     pch=16, col="#808D8D4A",
     ylab="log2(up_ratio/down_ratio)", xlab="log2(down_ratio/up_ratio)",
     main="screen2 up/mid v down/mid") 
abline(a=0,b=1)
abline(a=-3.3,b=1)
abline(a=3.3,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
dev.off()

#pull the things that are >=2X enriched

screen2_high_frag$up_div_mid = screen2_high_frag$up/screen2_high_frag$mid
screen2_high_frag$down_div_mid = screen2_high_frag$down/screen2_high_frag$mid

screen2_high_frag$up_ratio_div_down_ratio = screen2_high_frag$up_div_mid/screen2_high_frag$down_div_mid
screen2_high_frag$down_ratio_div_up_ratio = screen2_high_frag$down_div_mid/screen2_high_frag$up_div_mid

screen2_high_frag$gene_up = screen2_up[match(screen2_high_frag$uniques, screen2_up$uniques),"gene_name"]
screen2_high_frag$gene_down = screen2_down[match(screen2_high_frag$uniques, screen2_down$uniques),"gene_name"]
screen2_high_frag$gene_mid = screen2_mid[match(screen2_high_frag$uniques, screen2_mid$uniques),"gene_name"]

screen2_high_frag$gene_name =
  ifelse(is.na(screen2_high_frag$gene_up), 
         ifelse(is.na(screen2_high_frag$gene_down), screen2_high_frag$gene_mid, screen2_high_frag$gene_down), 
         screen2_high_frag$gene_up)
sum(is.na(screen2_high_frag$gene_name))

screen2_high_frag$gene_useful = xref[match(screen2_high_frag$gene_name, xref$V4),"V5"]
screen2_high_frag$gene_desc = xref[match(screen2_high_frag$gene_name, xref$V4),"V16"]

head(screen2_high_frag)

write.csv(filter(screen2_high_frag, down_ratio_div_up_ratio>=2), "screen2_down_ratio_enriched.csv")
write.csv(filter(screen2_high_frag, up_ratio_div_down_ratio>=2), "screen2_up_ratio_enriched.csv")



#CDFs on gene count

#screen 1
screen1_up_gene_table = as.data.frame(table(screen1_up$gene_name))
nrow(screen1_up_gene_table)
screen1_up_gene_cdf = ecdf(screen1_up_gene_table$Freq)
plot(screen1_up_gene_cdf)

screen1_up_gene_inframe_cdf = ecdf(as.data.frame(table(filter(screen1_up, joint_frame=="0,1")[,"gene_name"]))[,"Freq"])
plot(screen1_up_gene_inframe_cdf)

#screen 2
screen2_up_gene_table = as.data.frame(table(screen2_up$gene_name))
nrow(screen2_up_gene_table)
screen2_up_gene_cdf = ecdf(screen2_up_gene_table$Freq)
plot(screen2_up_gene_cdf)

screen2_up_gene_inframe_cdf = ecdf(as.data.frame(table(filter(screen2_up, joint_frame=="0,1")[,"gene_name"]))[,"Freq"])
plot(screen2_up_gene_inframe_cdf)

#look for repeats in up/mid v down/mid ratios
#repeats in down
filter(as.data.frame(table(filter(screen1_high_frag, down_ratio_div_up_ratio>=2)[,"gene_name"])), Freq!=1)
filter(as.data.frame(table(filter(three_high_frag, down_ratio_div_up_ratio>=2)[,"gene_name"])), Freq!=1)
filter(as.data.frame(table(filter(four_high_frag, two_ratio_div_one_ratio>=2)[,"gene_name"])), Freq!=1)

filter(as.data.frame(table(rbind(as.data.frame(table(filter(screen1_high_frag, down_ratio_div_up_ratio>=2)[,"gene_name"])),
                                 as.data.frame(table(filter(three_high_frag, down_ratio_div_up_ratio>=2)[,"gene_name"])),
                                 as.data.frame(table(filter(four_high_frag, two_ratio_div_one_ratio>=2)[,"gene_name"])))[,"Var1"])), Freq!=1)

filter(screen1_high_frag, gene_name=="YIL030C")
filter(screen1_high_frag, gene_name=="YMR295C")

#repeats in up
filter(as.data.frame(table(filter(screen1_high_frag, up_ratio_div_down_ratio>=2)[,"gene_name"])), Freq!=1)
filter(as.data.frame(table(filter(three_high_frag, up_ratio_div_down_ratio>=2)[,"gene_name"])), Freq!=1)
filter(as.data.frame(table(filter(four_high_frag, one_ratio_div_two_ratio>=2)[,"gene_name"])), Freq!=1)

filter(as.data.frame(table(rbind(as.data.frame(table(filter(screen1_high_frag, up_ratio_div_down_ratio>=2)[,"gene_name"])),
                                 as.data.frame(table(filter(three_high_frag, up_ratio_div_down_ratio>=2)[,"gene_name"])),
                                 as.data.frame(table(filter(four_high_frag, one_ratio_div_two_ratio>=2)[,"gene_name"])))[,"Var1"])), Freq!=1)

#read post made repeat file
repeats = read.csv("repeat_genes_count.csv", stringsAsFactors=FALSE, header=TRUE)
repeats  

repeats$gene_useful = xref[match(repeats$X, xref$V4),"V5"]
repeats$gene_desc = xref[match(repeats$X, xref$V4),"V16"]

write.csv(repeats, "repeats_genes_deets.csv")

screen1_high_frag[grep(pattern = "YMR295C", screen1_high_frag$gene_name),1:4]
three_high_frag[grep(pattern = "YMR295C", three_high_frag$gene_name),1:4]
four_high_frag[grep(pattern = "YMR295C", four_high_frag$gene_name),1:4]
head(screen1_high_frag)
# ##try using deseq
# library("DESeq2")
# 
#   #make the table with counts
#   #pull just the fragment counts (fragment read counts)
#   screen1_high_frag  
#   two_deseq = (screen1_high_frag)[,2:4]
#   #make the rownames the unique fragment identifiers
#   rownames(two_deseq) = c(as.character(screen1_high_frag[,"uniques"]))
#   two_deseq
#   #have to adjust the pseudocounts back from 0.1 to 0 (DESeq doesn't like non integers b/c expecting counts)
#   two_deseq$up = ifelse(screen1_high_frag$up==0.1, 0, screen1_high_frag$up)
#   two_deseq$mid = ifelse(screen1_high_frag$mid==0.1, 0, screen1_high_frag$mid)
#   two_deseq$down = ifelse(screen1_high_frag$down==0.1, 0, screen1_high_frag$down)
# 
# 
#   #make the table with descriptors
#   two_deseq_desc = data.frame(colnames(two_deseq))
#   colnames(two_deseq_desc) = "lib"
#   rownames(two_deseq_desc) = two_deseq_desc$lib
#   two_deseq_desc
# 
# 
#   #run DESeq  
#   ?DESeqDataSetFromMatrix
#   two_deseq_dataset =  DESeqDataSetFromMatrix(countData=two_deseq,
#                                               colData=two_deseq_desc,
#                                               design = ~ lib)
#   
#   two_deseq_dataset
#   design(two_deseq_dataset) = formula(~lib)
#   
#   two_deseqed = DESeq(two_deseq_dataset, test="LRT", reduced=~1)
#   two_deseqed_res = results(two_deseqed)
#   two_deseqed_res
#   
#   two_deseqed_res[order(two_deseqed_res$stat, decreasing = TRUE),]
#   
#   two_deseq["chrIV_892910_893488_+",]
#   
#   sum(two_deseqed_res$padj<.1, na.rm=TRUE)
#   
#   plotDispEsts(two_deseqed)
#  ----->> DIDN'T WORK
