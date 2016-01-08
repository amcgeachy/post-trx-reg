getwd()
setwd("/Users/annamcgeachy/Google Drive/post trx reg data/datafiles_screen4/")

##set up A
  one_full_down_hi = read.csv("one_full_down_hi.csv")
  
  one_full_down_hi_uniques = apply(one_full_down_hi, 1, pasting_together)
  one_full_down_hi$uniques = one_full_down_hi_uniques
  
  #add in the annotation
  one_full_down_hi$gene_useful = xref[match(one_full_down_hi$gene_name, xref$V4),"V5"]
  one_full_down_hi$gene_desc = xref[match(one_full_down_hi$gene_name, xref$V4),"V16"]
  
  #pull in frame fragments
  one_full_down_hi_01 = filter(one_full_down_hi, joint_frame=="0,1")

##set up B
  one_full_down_mi = read.csv("one_full_down_mi.csv")
  
  one_full_down_mi_uniques = apply(one_full_down_mi, 1, pasting_together)
  one_full_down_mi$uniques = one_full_down_mi_uniques
  
  #add in the annotation
  one_full_down_mi$gene_useful = xref[match(one_full_down_mi$gene_name, xref$V4),"V5"]
  one_full_down_mi$gene_desc = xref[match(one_full_down_mi$gene_name, xref$V4),"V16"]
  
  #pull in frame fragments
  one_full_down_mi_01 = filter(one_full_down_mi, joint_frame=="0,1")