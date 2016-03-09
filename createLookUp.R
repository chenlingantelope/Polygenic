
# this is to create a look up table with columns: current alle freq, selection coefficient, nSL mean, nSL std

# use msms to simulate data under selection
# use nSL to compute statistic

# Rscript params_file.txt output_file.txt

args = commandArgs(trailingOnly = TRUE);
params=args[1] 
# file with all parameter list:
# tab separated
##
#NCHROMS	100
#NE	10000
#NITER	1000
#SELCOEFFS	0,100,10
#SELCOEFFS	150,1200,50
#CURRFREQS	0.1,0.9,0.2
#MSMS	~/Documents/Software/msms/bin/msms
#NSL	~/Documents/Software/nSL/nSL_matteo
#TEMPFILE	out
#THETA	88
#RECOMB	101
#NSITES	200000
#MAXLEN	200
#VERBOSE	0
outfile=args[2]
rm(args)

# READ AND ASSIGN PARAMETERS

params=read.table(params, head=F, stringsAsFactors=F, sep="\t")
# assign:
# population effect size, we also assume std neutral const-size pop model
Ne=as.numeric(params[which(params[,1]=="NE"),2])[1]
# how many chromosomes
nchroms=as.numeric(unlist(strsplit(params[which(params[,1]=="NCHROMS"),2],split=",")))
# how many simulations per scenario
niter=as.numeric(params[which(params[,1]=="NITER"),2])[1]
# selection coefficients in 2NeS format!
tmp=params[which(params[,1]=="SELCOEFFS"),2]
selcoeffs=c()
for (i in 1:length(tmp)) {
	ints=as.numeric(unlist(strsplit(tmp[i],split=",")))
	selcoeffs=c(selcoeffs, seq(ints[1],ints[2],ints[3]) )
}
selcoeffs=sort(unique(selcoeffs))
rm(tmp)
rm(ints)
# current allele frequencies of selected site
tmp=params[which(params[,1]=="CURRFREQS"),2]
currentfreq=c()
for (i in 1:length(tmp)) {
        ints=as.numeric(unlist(strsplit(tmp[i],split=",")))
        currentfreq=c(currentfreq, seq(ints[1],ints[2],ints[3]) )
}
currentfreq=sort(unique(currentfreq))
rm(tmp)
rm(ints)
# msms and nSL directory
msms_dir=as.character(params[which(params[,1]=="MSMS"),2])[1]
nsl_dir=as.character(params[which(params[,1]=="NSL"),2])[1]
# temp file for msms and nsl
ftemp=as.character(params[which(params[,1]=="TEMPFILE"),2])[1]
# verbose
verbose=as.numeric(params[which(params[,1]=="VERBOSE"),2])[1]
# theta
theta=as.numeric(params[which(params[,1]=="THETA"),2])[1]
# rho
rho=as.numeric(params[which(params[,1]=="RHO"),2])[1]
# nsites
nsites=as.numeric(params[which(params[,1]=="NSITES"),2])[1]
# maxlen for nSL
maxlen=as.numeric(params[which(params[,1]=="MAXLEN"),2])[1]

# keep all floats
options(scipen=999)

if (verbose) print(params)

# this will be the output for each scenario
lookup=c()

# initialise the output file
cat("", file=outfile);

# cycle over all ncrhoms, selcoeffs, current freqs
for (c in nchroms) {

	for (s in selcoeffs) {

		for (f in currentfreq) {

			nsl=c()
			n=0;

			if (verbose) cat("\n",c,s,f)

			while (n<niter) {

				# simulate
				cmd=paste(msms_dir," -N ", Ne, " -ms ", c, " 1 -t ", theta, " -r ", rho, " ", nsites, " -SAA ", s, " -SAa ", s/2, " -Sp 0.50 -SF 0 ", f, " -Smark > ", ftemp, ".msms", sep="", collapse="");
				system(cmd);
				system(paste("sed -i '1s/.*/-ms ", c, " 1 -t ", theta, " -r ", rho, " ", nsites, "/' ", ftemp, ".msms", sep="", collapse=""));

				# compute stat
				system(paste(nsl_dir, " -msfile ", ftemp, ".msms -maxLen ", maxlen, " -msLen ", nsites, " 2> ", ftemp, ".log.nsl > ", ftemp, ".nsl", sep="", collapse=""))
				values=read.table(paste(ftemp, ".nsl", sep="", collapse=""), header=T, stringsAsFactors=F)
				ind=which(values[,2]==nsites/2);
				if (length(ind)==1) {
					nsl=c(nsl, values[ind,3])
					n=n+1;
				} 
				if (verbose & length(ind)>1) {
					cat("\nProblems here:")
					print(values[(ind[1]-3):(ind[length(ind)]+3),])
					cat("\n")
				}
				if (verbose & length(ind)==0) cat("\nNone") 
			}

			lookup=c(c,s,f,mean(nsl),sd(nsl))
			if (verbose) cat(":",mean(nsl),sd(nsl))
			cat(lookup, sep="\t", "\n", append=T, file=outfile);
		}
	}
}

system(paste("rm ", ftemp,".*", sep="", collapse=""), ignore.stderr=T)


