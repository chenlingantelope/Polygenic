args = commandArgs(trailingOnly = TRUE);
filename=args[1];


library("minpack.lm")
look_file=paste(filename,".txt",sep="")
out_file=paste(filename,"_inter.txt",sep="")
log_file=paste(filename,"_inter.log",sep="")
pdf_file=paste(filename,"_inter.pdf",sep="")

look=read.table(file=look_file, sep="\t", stringsAsFac=F)[,1:5]
#colnames(look)=c("nloci", "selcoeff", "currfreq", "mean", "sd")

nchroms=sort(unique(look[,1]))
x=sort(unique(look[,2])) # sel coeffs
cfreq=sort(unique(look[,3]))

cat("", file=out_file)
cat("", file=log_file)

for (i in nchroms) {
	for (j in cfreq) {

		ind = which(look[,1]==i & look[,3]==j)
		y1 = look[ind,4] # nSL mean values
		y2 = look[ind,5] # nSL sd values

		data=data.frame(y1,y2,x)
		fit_mean=nlsLM
		fit_mean  = nlsLM( y1 ~ p*exp(x*lambda1) + (1-p)*exp(x*lambda2) + cons, data=data, start= list(p=0.9,lambda1=-0.1, lambda2=-0.1,cons=-4))
        fit_sd  = nlsLM( y2 ~ p*exp(x*lambda1) + (1-p)*exp(x*lambda2) + cons, data=data, start= list(p=0.9,lambda1=-0.1, lambda2=-0.1,cons=-3))
		y1_pred=as.numeric(predict(fit_mean, data.frame(x=x)))
        y2_pred=as.numeric(predict(fit_sd, data.frame(x=x)))

		for (b in 1:length(y1_pred)) cat(i, x[b], j, y1_pred[b], y2_pred[b], sep="\t", "\n", file=out_file, append=T)
	}
}

int=read.table(out_file)
cols=rainbow(length(cfreq))
titles=c("Mean","Var")
pdf(pdf_file)

for(kind in c(4,5)){
	plot(x=-100,y=-100,pch=16,xlab="selection coefficient",ylab="nSL",main=titles[kind-3],xlim=range(look[,2]),ylim=range(look[,kind]))
	for(i in c(1:length(cfreq))){
		freq=cfreq[i]
		points(x=look[look[,3]==freq,2],y=look[look[,3]==freq,kind],pch=16,col=cols[i])
		points(x=int[int[,3]==freq,2],y=int[int[,3]==freq,kind],type="l",col=cols[i])
	}
	legend(1000, 0,c(0.1,0.3,0.5,0.7,0.9), cex=0.8, pch=16,col=cols)
}

dev.off()