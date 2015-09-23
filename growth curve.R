getwd()
setwd("/Users/annamcgeachy/Google Drive/post-trx reg/")

tecan = read.csv("growth_curve_20150331.csv")

#make column for times even across samples
times = apply(tecan[,c(grep("Time", colnames(tecan)))],1,mean) 
minutes = round(times/1000/60, digits = 0)
minutes

colnames(tecan)

#pull just the columns with samples
  #generate row/column combos
  col_of_interest = c(paste("A",c(1:9,12),sep=""), 
    paste("B",c(1:9,12),sep=""), 
    paste("C",c(1:9),sep=""),
    paste("D",c(1:9),sep=""))

  #subset data
  samples = tecan[,col_of_interest]

  #then add in time
  samples$time = minutes
  head(samples)

#plot the background change over time to figure out how to handle changes
plot(samples$time, samples$A9, pch=16, ylim=c(.11, .14),
     main="background change over time")
points(samples$time, samples$B9, pch=16, col="blue")
points(samples$time, samples$C9, pch=16, col="red")
points(samples$time, samples$D9, pch=16, col="green")

  #--> there is some fluctuation and differences between rows, so probably should do it row by row

#subtract background

  #make matrix of background
  bg = matrix(c(rep(samples$A9, 10),
           rep(samples$B9, 10),
           rep(samples$C9, 9),
           rep(samples$D9, 9),
           rep(0, 96)),
           nrow=96)       
    head(bg)
  #subtract it out
  adjusted = samples - bg
  head(adjusted)


#plot unlogged
  #unlogged for 116
  pdf("unlogged 116 growth curves, 20150331.pdf", useDingbats = FALSE)
  plot(adjusted$time, adjusted$A1, pch=16, main="116 strains, unlogged", ylim=c(0,1))
  points(adjusted$time, adjusted$C1, pch=16)
  points(adjusted$time, adjusted$A2, pch=16, col=2)
  points(adjusted$time, adjusted$C2, pch=16, col=2)
  points(adjusted$time, adjusted$A3, pch=16, col=3)
  points(adjusted$time, adjusted$C3, pch=16, col=3)
  points(adjusted$time, adjusted$A4, pch=16, col=4)
  points(adjusted$time, adjusted$C4, pch=16, col=4)
  points(adjusted$time, adjusted$A5, pch=16, col=5)
  points(adjusted$time, adjusted$C5, pch=16, col=5)
  points(adjusted$time, adjusted$A6, pch=16, col=6)
  points(adjusted$time, adjusted$C6, pch=16, col=6)
  points(adjusted$time, adjusted$A7, pch=16, col=7)
  points(adjusted$time, adjusted$C7, pch=16, col=7)
  points(adjusted$time, adjusted$A8, pch=16, col=8)
  points(adjusted$time, adjusted$C8, pch=16, col=8)
  
  legend(1000, .4, "116,BFP-LN", pch=16, col=1, bty='n', cex=.75)
  legend(1000, .375, "116,DHH1-LN", pch=16, col=2, bty='n', cex=.75)
  legend(1000, .35, "116,POP2-LN", pch=16, col=3, bty='n', cex=.75)
  legend(1000, .325, "116,PAB1-LN", pch=16, col=4, bty='n', cex=.75)
  legend(1000, .3, "116,CDC19-LN", pch=16, col=5, bty='n', cex=.75)
  legend(1000, .275, "116,YPR204W-LN", pch=16, col=6, bty='n', cex=.75)
  legend(1000, .25, "116,PAT1-LN", pch=16, col=7, bty='n', cex=.75)
  legend(1000, .225, "116,YMR295C-LN", pch=16, col=8, bty='n', cex=.75)
  dev.off()
  
  #unlogged for 117
  pdf("unlogged 117 growth curves, 20150331.pdf", useDingbats = FALSE)
  plot(adjusted$time, adjusted$B1, pch=16, main="117 strains, unlogged", ylim=c(0,1))
  points(adjusted$time, adjusted$D1, pch=16)
  points(adjusted$time, adjusted$B2, pch=16, col=2)
  points(adjusted$time, adjusted$D2, pch=16, col=2)
  points(adjusted$time, adjusted$B3, pch=16, col=3)
  points(adjusted$time, adjusted$D3, pch=16, col=3)
  points(adjusted$time, adjusted$B4, pch=16, col=4)
  points(adjusted$time, adjusted$D4, pch=16, col=4)
  points(adjusted$time, adjusted$B5, pch=16, col=5)
  points(adjusted$time, adjusted$D5, pch=16, col=5)
  points(adjusted$time, adjusted$B6, pch=16, col=6)
  points(adjusted$time, adjusted$D6, pch=16, col=6)
  points(adjusted$time, adjusted$B7, pch=16, col=7)
  points(adjusted$time, adjusted$D7, pch=16, col=7)
  points(adjusted$time, adjusted$B8, pch=16, col=8)
  points(adjusted$time, adjusted$D8, pch=16, col=8)
  
  legend(1000, .4, "117,BFP-LN", pch=16, col=1, bty='n', cex=.75)
  legend(1000, .375, "117,DHH1-LN", pch=16, col=2, bty='n', cex=.75)
  legend(1000, .35, "117,POP2-LN", pch=16, col=3, bty='n', cex=.75)
  legend(1000, .325, "117,PAB1-LN", pch=16, col=4, bty='n', cex=.75)
  legend(1000, .3, "117,CDC19-LN", pch=16, col=5, bty='n', cex=.75)
  legend(1000, .275, "117,YPR204W-LN", pch=16, col=6, bty='n', cex=.75)
  legend(1000, .25, "117,PAT1-LN", pch=16, col=7, bty='n', cex=.75)
  legend(1000, .225, "117,YMR295C-LN", pch=16, col=8, bty='n', cex=.75)
  dev.off()

#plot logged

  #log things
  adjusted_log = apply(adjusted[,1:38], 2, log2)
  adjusted_log = as.data.frame(adjusted_log)
  adjusted_log[,"time"] = adjusted$time
  head(adjusted_log)

  #logged for 116
  pdf("logged 116 growth curves, 20150331.pdf", useDingbats = FALSE)
  plot(adjusted_log$time, adjusted_log$A1, pch=16, main="116 strains, unlogged", ylim=c(-9,0))
  points(adjusted_log$time, adjusted_log$C1, pch=16)
  points(adjusted_log$time, adjusted_log$A2, pch=16, col=2)
  points(adjusted_log$time, adjusted_log$C2, pch=16, col=2)
  points(adjusted_log$time, adjusted_log$A3, pch=16, col=3)
  points(adjusted_log$time, adjusted_log$C3, pch=16, col=3)
  points(adjusted_log$time, adjusted_log$A4, pch=16, col=4)
  points(adjusted_log$time, adjusted_log$C4, pch=16, col=4)
  points(adjusted_log$time, adjusted_log$A5, pch=16, col=5)
  points(adjusted_log$time, adjusted_log$C5, pch=16, col=5)
  points(adjusted_log$time, adjusted_log$A6, pch=16, col=6)
  points(adjusted_log$time, adjusted_log$C6, pch=16, col=6)
  points(adjusted_log$time, adjusted_log$A7, pch=16, col=7)
  points(adjusted_log$time, adjusted_log$C7, pch=16, col=7)
  points(adjusted_log$time, adjusted_log$A8, pch=16, col=8)
  points(adjusted_log$time, adjusted_log$C8, pch=16, col=8)
  
  legend(1000, -5, "116,BFP-LN", pch=16, col=1, bty='n', cex=.75)
  legend(1000, -5.25, "116,DHH1-LN", pch=16, col=2, bty='n', cex=.75)
  legend(1000, -5.5, "116,POP2-LN", pch=16, col=3, bty='n', cex=.75)
  legend(1000, -5.75, "116,PAB1-LN", pch=16, col=4, bty='n', cex=.75)
  legend(1000, -6, "116,CDC19-LN", pch=16, col=5, bty='n', cex=.75)
  legend(1000, -6.25, "116,YPR204W-LN", pch=16, col=6, bty='n', cex=.75)
  legend(1000, -6.5, "116,PAT1-LN", pch=16, col=7, bty='n', cex=.75)
  legend(1000, -6.75, "116,YMR295C-LN", pch=16, col=8, bty='n', cex=.75)
  dev.off()
  
  #unlogged for 117
  pdf("logged 117 growth curves, 20150331.pdf", useDingbats = FALSE)
  plot(adjusted_log$time, adjusted_log$B1, pch=16, main="117 strains, unlogged", ylim=c(-9,0))
  points(adjusted_log$time, adjusted_log$D1, pch=16)
  points(adjusted_log$time, adjusted_log$B2, pch=16, col=2)
  points(adjusted_log$time, adjusted_log$D2, pch=16, col=2)
  points(adjusted_log$time, adjusted_log$B3, pch=16, col=3)
  points(adjusted_log$time, adjusted_log$D3, pch=16, col=3)
  points(adjusted_log$time, adjusted_log$B4, pch=16, col=4)
  points(adjusted_log$time, adjusted_log$D4, pch=16, col=4)
  points(adjusted_log$time, adjusted_log$B5, pch=16, col=5)
  points(adjusted_log$time, adjusted_log$D5, pch=16, col=5)
  points(adjusted_log$time, adjusted_log$B6, pch=16, col=6)
  points(adjusted_log$time, adjusted_log$D6, pch=16, col=6)
  points(adjusted_log$time, adjusted_log$B7, pch=16, col=7)
  points(adjusted_log$time, adjusted_log$D7, pch=16, col=7)
  points(adjusted_log$time, adjusted_log$B8, pch=16, col=8)
  points(adjusted_log$time, adjusted_log$D8, pch=16, col=8)
  
  legend(1000, -5, "117,BFP-LN", pch=16, col=1, bty='n', cex=.75)
  legend(1000, -5.25, "117,DHH1-LN", pch=16, col=2, bty='n', cex=.75)
  legend(1000, -5.5, "117,POP2-LN", pch=16, col=3, bty='n', cex=.75)
  legend(1000, -5.75, "117,PAB1-LN", pch=16, col=4, bty='n', cex=.75)
  legend(1000, -6, "117,CDC19-LN", pch=16, col=5, bty='n', cex=.75)
  legend(1000, -6.25, "117,YPR204W-LN", pch=16, col=6, bty='n', cex=.75)
  legend(1000, -6.5, "117,PAT1-LN", pch=16, col=7, bty='n', cex=.75)
  legend(1000, -6.75, "117,YMR295C-LN", pch=16, col=8, bty='n', cex=.75)
  dev.off()

#calculate growth rates

  #pull the region where things are growing exponentially
  temp = adjusted_log[(adjusted_log>(-5) & adjusted_log<(-3)),c("B1","time")]
  (temp)
  
  #change do regression to get slope
  blah = lm(temp$B1~temp$time)
  blah
  blah$coefficients[2] #the slope
  gen_time = 1/blah$coefficients[2] #the generation time from slope
  gen_time
  
  plot(adjusted_log$time, adjusted_log$B1)
  abline(a=blah$coefficients[1], b=blah$coefficients[2])
  
  plot(temp$time, temp$B1)
  abline(a=blah$coefficients[1], b=blah$coefficients[2])

gen_time_fxn = function(well){
  #pull the region where things are growing exponentially
  temp = adjusted_log[(adjusted_log>(-5) & adjusted_log<(-3)),c(well,"time")]
  temp
  
  #change do regression to get slope
  blah = lm(temp[,well]~temp$time)
  blah
  blah$coefficients[2] #the slope
  plot(temp)
  
  pdf(sprintf("%s regression.pdf", well), useDingbats = FALSE)
  plot(adjusted_log$time, adjusted_log[,well], main=sprintf("%s", well))
  legend(x=250, y=-4.5, round(1/blah$coefficients[2], digits = 1), bty='n')
  abline(a=blah$coefficients[1], b=blah$coefficients[2])
  
  plot(temp$time, temp[,well], main=sprintf("%s", well))
  abline(a=blah$coefficients[1], b=blah$coefficients[2])
  legend(x=250, y=-4.5, round(1/blah$coefficients[2], digits = 1), bty='n')
  dev.off()
  
  gen_time = 1/blah$coefficients[2] #the generation time from slope
  gen_time
  
}


gen_time_fxn("B1")
colnames(adjusted_log)
#need to pull all wells that aren't currently -inf (the blanks, well 9)
growth_results = as.matrix(sapply(colnames(adjusted_log)[grep("9", colnames(adjusted_log), invert=TRUE)], gen_time_fxn))

key = c("116,BFP-LN", "116,DHH1-LN", "116,POP2-LN", "116,PAB1-LN", 
"116,CDC19-LN", "116,YPR204W-LN", "116,PAT1-LN" ,"116,YMR295C-LN",
"117-BFP-LN*",
"116,BFP-LN", "116,DHH1-LN", "116,POP2-LN", "116,PAB1-LN", 
"116,CDC19-LN", "116,YPR204W-LN", "116,PAT1-LN" ,"116,YMR295C-LN",
"117-BFP-LN*",
"117,BFP-LN", "117,DHH1-LN", "117,POP2-LN", "117,PAB1-LN", 
"117,CDC19-LN", "117,YPR204W-LN", "117,PAT1-LN" ,"117,YMR295C-LN",
"117,BFP-LN", "117,DHH1-LN", "117,POP2-LN", "117,PAB1-LN", 
"117,CDC19-LN", "117,YPR204W-LN", "117,PAT1-LN" ,"117,YMR295C-LN",
"time")

fits = c("ok", "not great", "ok", "ok", "messy", "ok", "ok", "terrible", "messy", 
         "ok", "ok", "messy", "ok", "messy", "ok", "messy", "ok", "really messy", 
         "ok", "ok", "ok", "ok", "a little messy", "a little messy", "ok", "ok", 
         "pretty messy", "really messy", "pretty messy", "messy", "messy", "really messy", "really messy", "ok", "whatever") 

growth_results = as.data.frame(cbind(growth_results, key, fits))
growth_results
colnames(growth_results)[1] = "doubling_time"
growth_results
growth_results = growth_results[order(growth_results$key),]

write.csv(growth_results, "20150331 growth results.csv")
?write.csv
