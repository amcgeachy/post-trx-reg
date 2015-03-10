setwd("/Users/annamcgeachy/")

recomb_frag_unique = read.delim(file = "nonrecomb_frag_count_sort.txt", header=F, sep = "")
?read.delim
head(recomb_frag_unique)
recomb_frag_unique = recomb_frag_unique[2:dim(recomb_frag_unique)[1],]

recomb_frag_unique[,"length_nt"] = (recomb_frag_unique$V7 - recomb_frag_unique$V3 + 1)
recomb_frag_unique[,"length_aa"] = (recomb_frag_unique$V7 - recomb_frag_unique$V3 + 1)/3
recomb_frag_unique[,"length_codon_nt"] = (recomb_frag_unique$V7 - recomb_frag_unique$V3 + 1) %% 3

hist(recomb_frag_unique[,"length_nt"], xlim=c(0,400))
table(recomb_frag_unique$length_nt)

summary(recomb_frag_unique$length_nt)
?boxplot
boxplot(recomb_frag_unique$length_nt)
blahhh = recomb_frag_unique[recomb_frag_unique$length_nt>0, ]
head(blahhh)
blahhh2 = recomb_frag_unique[recomb_frag_unique$length_nt<400 & recomb_frag_unique$length_nt>0,]
hist(blahhh2$length_nt,
     main = "Distribution of recombined insert sizes",
     xlab="insert size", col="grey")
head(blahblah$length_nt)
length(blahblah)
hist(blahblah$length_nt)
summary(recomb_frag_unique$length_nt)
hist_breaks = seq(from=0, to=52000, by=2000)
hist_breaks
hist(recomb_frag_unique[recomb_frag_unique$length_codon_nt==0, "V1"], ylim=c(0,600), xlim=c(0,50000),
     xlab="number of sequencing reads per unique fragment", ylab="frequency of sequencing read number",
     main="Histogram of sequencing read number per unique fragments", col="blue", breaks=hist_breaks)
hist(recomb_frag_unique[recomb_frag_unique$length_codon_nt==1, "V1"], ylim=c(0,600), xlim=c(0,50000),
     xlab="number of sequencing reads per unique fragment", ylab="frequency of sequencing read number",
     main="Histogram of sequencing read number for +1 frame fragments", col="red", breaks=hist_breaks, add=T)
hist(recomb_frag_unique[recomb_frag_unique$length_codon_nt==2, "V1"], ylim=c(0,600), xlim=c(0,50000),
     xlab="number of sequencing reads per unique fragment", ylab="frequency of sequencing read number",
     main="Histogram of sequencing read number for +1 frame fragments", col="grey", breaks=hist_breaks, add=T)
legend("topright", c("In frame", "+1 frame", "+2 frame"), fill=c("blue", "red", "grey"), cex=.8)
?legend
?hist


hist(recomb_frag_unique$length_codon_nt)
?hist
hist(recomb_frag_unique[1,])
re
head(recomb_frag_unique[,1])
colnames(recomb_frag_unique)
summary(recomb_frag_unique$length)
head(recomb_frag_unique[,"V2"])
