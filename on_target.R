dir<-Sys.getenv("outDir")
sample<-Sys.getenv("sample")

setwd(paste0(dir,"/depth"))

cover<-read.table(paste0(sample,".all.cov"))
cov_cumul<-1-cumsum(cover[,5])
pdf(file=paste0("on_target_depth_",sample,".pdf"))
plot(cover[1:200,2], cov_cumul[1:200], type ='l', xlab="Depth", ylab="Fraction of capture target bases", ylim=c(0,1.0), main=paste0("Target Region coverage sample ",sample))
abline(v=25, col="gray60")
abline(v=50, col="gray60")
abline(v=80, col="gray60")
abline(h=0.9, col="gray60")
dev.off()
