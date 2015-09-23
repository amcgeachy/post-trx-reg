getwd()
setwd("/Users/annamcgeachy/Google Drive/post-trx reg/")

gen_table = read.csv("generations to make majority be selected simple.csv")
head(gen_table)
gen_table = gen_table[0:31,]
gen_table[,"log_nontransformed"] = log10(gen_table$non.transformed)
gen_table[,"log_transformed"] = log10(gen_table$transform)

#uncapped growth
pdf("transform v untransformed uncapped.pdf", useDingbats = FALSE, width=10, height=4)
plot(gen_table$log_nontransformed, type='n', yaxt = 'n',
     xlab='generations', ylab='# of cells', ylim=c(3,14))
lines(gen_table$log_nontransformed, col="black", cex=2)
lines(gen_table$log_transformed, col="blue")
axis(2, at=seq(from=(0), to=14, by=2), labels=c("10e0","10e2","10e4","10e6","10e8","10e10","10e12", "10e14"), las=2)
dev.off()


#capped growth
gen_table
?pdf
pdf("transform v untransformed capped.pdf", useDingbats = FALSE, width=10, height=4)
plot(gen_table$fract.non.trans, type='n',
         xlab='generations', ylab='fraction of population', ylim=c(0,1), xlim=c(0,30))
lines(gen_table$fract.non.trans, col="black", cex=2)
lines(gen_table$fract.trans, col="blue")
dev.off()