one_full_down = read.csv("one_full_down.csv")

one_full_down_uniques = apply(one_full_down, 1, pasting_together)
one_full_down$uniques = one_full_down_uniques

#add in the annotation
one_full_down$gene_useful = xref[match(one_full_down$gene_name, xref$V4),"V5"]
one_full_down$gene_desc = xref[match(one_full_down$gene_name, xref$V4),"V16"]

#pull in frame fragments
one_full_down_01 = filter(one_full_down, joint_frame=="0,1")

#
one_up_down_uniques = c(one_full_down_01$uniques, one_full_up_01$uniques)
one_up_down_uniques = as.data.frame(one_up_down_uniques)
head(one_up_down_uniques)
one_up_down_uniques$up_count = one_full_up_01[match(one_full_up_01$gene_name, one_up_down_uniques),"V16"]

one_full_down_01[match(one_full_up_01$uniques, one_up_down_uniques$one_up_down_uniques),"frag_count"][1]

one_up_down_uniques[,"up_frag"] = one_full_up_01[match(one_up_down_uniques$one_up_down_uniques, one_full_up_01$uniques), "frag_count"]
one_up_down_uniques[,"down_frag"] = one_full_down_01[match(one_up_down_uniques$one_up_down_uniques, one_full_down_01$uniques), "frag_count"]

head(one_up_down_uniques)

plot(log(one_up_down_uniques$adj_down, 2), log(one_up_down_uniques$adj_up, 2), pch=16, col="#808D8D4A", 
     xlab="log2(unique fragment counts in down library)",
     ylab="log2(unique fragment counts in up library)",
     main="unique fragments counts in up v down")
abline(a=0,b=1)
abline(a=-4.6,b=1)
abline(a=4.6,b=1)
abline(a=-2.3,b=1)
abline(a=2.3,b=1)

filter(one_up_down_uniques, is.na(up_frag), down_frag=max(one_up_down_uniques$down_frag))

max(one_up_down_uniques$down_frag, na.rm=TRUE)
grep(15800, one_up_down_uniques$down_frag)
one_up_down_uniques[6716,]
one_up_down_uniques[23210,]

max(one_up_down_uniques$up_frag, na.rm=TRUE)
grep(6827, one_up_down_uniques$up_frag)
one_up_down_uniques[13473,]
one_up_down_uniques[29207,]

one_full_down_01[grep("chrXIV_337015_337122_+", one_full_down_01$uniques),]

one_full_up_01[grep("chrX_728919_729173_-", one_full_up_01$uniques),]

log(.1,2)
one_up_down_uniques$adj_up = 
  
switch_NAs_up = function(x){ 
  ifelse(is.na(x["up_frag"]), .1, as.numeric(x["up_frag"]))
}

switch_NAs_down = function(x){ 
  ifelse(is.na(x["down_frag"]), .1, as.numeric(x["down_frag"]))
}


one_up_down_uniques$adj_up = apply(one_up_down_uniques, 1, switch_NAs_up)


plot(log(one_up_down_uniques$adj_down, 2), log(one_up_down_uniques$adj_up, 2), pch=16, col="#808D8D4A", 
     xlab="log2(unique fragment counts in down library)",
     ylab="log2(unique fragment counts in up library)",
     main="unique fragments counts in up v down")
abline(a=0,b=1)
abline(a=-4.6,b=1)
abline(a=4.6,b=1)
abline(a=-2.3,b=1)
abline(a=2.3,b=1)


?ifelse
log(NA, 2)
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
