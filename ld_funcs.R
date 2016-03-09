
trim<-function (sequence, inner = FALSE) {
## Author: Giorgia Menozzi
    sequence <- sub("^ +", "", sequence)
    sequence <- sub(" +$", "", sequence)
    if (inner) {
        sequence <- gsub(" +", "", sequence)
    }
    sequence
}

readMs<-function(filein, nsam, len=c()) {

 lfile<-readLines(con=filein)
 inds<-which(lfile=="//")+3
 inde<-which(lfile=="//")+3+nsam-1
 #cat("Sample:",length(inds),"\n")
 indpos<-inds-1

 res<-pos<-list()
 for (i in 1:length(inds)) {
  res[[i]]<-trim(lfile[inds[i]:inde[i]])
  ss<-strsplit(lfile[indpos[i]],split=" ")[[1]]
  pos[[i]]<-as.numeric(ss[2:length(ss)])
  if (length(len)>0) pos[[i]]<-as.integer(pos[[i]]*len)
 }

 readMs<-list(hap=res, pos=pos)

}


linkage<-function(haplos, pos) {

 if (!is.list(haplos)) haplos<-list(haplos)
 npop<-length(haplos)

 res<-list()

 for (p in 1:npop) {

  polyss<-NA

  len<-nchar(haplos[[p]][1])
  nsam<-length(haplos[[p]])
 # cat("Sample:",nsam,"Sites:",len,"\n")
 # cat("\n",rep("",(floor(len/10)-1)),"|\n")
  Dprime<-rsquare<-rep(NA, len)

  i=pos;
  subi<-substring(haplos[[p]],i,i)
  ali<-unique(subi)
  if (length(ali)==1) ali<-c(ali, "N")

  for (j in 1:len) {

    # alleli
    subj<-substring(haplos[[p]],j,j)
    # forme alleliche
    alj<-unique(subj)
    if (length(alj)==1) alj<-c(alj, "N")

    # freq alleliche
    pp<-rep(NA,4)
    pp[1]<-length(which(subi==ali[1]))/nsam
    pp[2]<-length(which(subi==ali[2]))/nsam
    pp[3]<-length(which(subj==alj[1]))/nsam
    pp[4]<-length(which(subj==alj[2]))/nsam

    # aplotipi
    haps<-paste(subi,subj,sep="")

    # poss aplotipi
    posshap<-c()
    for (ai in 1:2) {
     for (aj in 1:2) {
      posshap<-c(posshap, paste(ali[ai],alj[aj],sep="",collapse=""))
     }
    }

    # freq aplotipiche
    P<-rep(NA,4) # prob aplotipi
    for (c in 1:4) {
     P[c]<-length(which(haps==posshap[c]))/nsam
    }

    D<-P[1]*P[4]-P[2]*P[3]
    #Dmax<-min(c(pp[1]*pp[4]),c(pp[2]*pp[3]))
    #Dmin<-max(c(-pp[1]*pp[3]),c(-pp[2]*pp[4]))
    #if (D<0) {
    # Dprime[j]<-abs(D/Dmin)
    #} else {
    # Dprime[j]<-abs(D/Dmax)
    #}

    # rsquare
    rsquare[j]<-(D^2)/prod(pp)

    # test
    #chiquadro<-rsquare[i,j]*nsam
    #pv[i,j]<-pchisq(q=chiquadro, df=1, lower.tail=F)

  } # fine for in j
 
  #res[[p]]<-list(Dprime=Dprime, rsquare=rsquare, pv=pv, snp=polyss)

  res[[p]]=rsquare

 } # fine for in p

 linkage<-res

} # fine function




