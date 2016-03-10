
# This script is to simulate genotypes and phenotypes

options(scipen=999)

# input 

args <- commandArgs(trailingOnly = TRUE);
K = 		as.numeric(args[1]); # number of loci 
nchroms = 	as.numeric(args[2]); # number of chroms
selcoeff = 	as.numeric(args[3]); # selection coefficient in 2Ns units (assuming Ne=1e4)
currentfreq = 	as.numeric(args[4]); # current allele frequency (use only 0.6 0.7 0.8 0.9), we assume it is known (it is indeed observable)
h2 =		as.numeric(args[5]); # hereditability (the same across all sites)
mu=c()
mu[1] =		as.numeric(args[6]); # mean phenotype for genotype 0
mu[2] =        	as.numeric(args[7]); # mean phenotype for genotype 1
mu[3] =        	as.numeric(args[8]); # mean phenotype for genotype 2
af = 		as.numeric(args[9]); # allelic effect
model =	0 	#as.numeric(args[10]); # model additive (0), dominant (1), complete (2)
lookup_file = 	args[11]; # lookup table file
output_file = 	args[12]; # output file text
niter = args[13]; # how many iterations per scenario
params = args[14]; # param file
cat(args[c(1:10,13)], sep="\t");
cat("\t")
rm(args);

source("ld_funcs.R")
K=20;nchroms=100;selcoeff=100;currentfreq=0.8;h2=0.75;mu=c(0,5,10);af=5;model=0;lookup_file="lookup_chenling_inter.txt"
output_file="03092016/sc_es.txt";niter=50;params="params_lookup.txt"
# read param file
params=read.table(params, head=F, stringsAsFactors=F, sep="\t")
Ne=as.numeric(params[which(params[,1]=="NE"),2])[1]
nchroms=as.numeric(unlist(strsplit(params[which(params[,1]=="NCHROMS"),2],split=",")))
msms_dir=as.character(params[which(params[,1]=="MSMS"),2])[1]
nsl_dir=as.character(params[which(params[,1]=="NSL"),2])[1]
theta=as.numeric(params[which(params[,1]=="THETA"),2])[1]
rho=as.numeric(params[which(params[,1]=="RHO"),2])[1]
nsites=as.numeric(params[which(params[,1]=="NSITES"),2])[1]
maxlen=as.numeric(params[which(params[,1]=="MAXLEN"),2])[1]

# read lookup table 
lookup=read.table(lookup_file, sep="\t", stringsAsFactors=F);
print("lookup ok")
# take only the subset 
lookup=lookup[which((lookup[,1]==nchroms) & (lookup[,3]==currentfreq)),] 

# interval in which calculate likelihood sel coeff
# THIS MUST BE THE SAME AS THE ONE USED FOR THE LOOKUP TABLE!
#selcoeff_bins=c(seq(0,100,10),seq(150,1300,50))
# here let's check it
selcoeff_bins=sort(unique(as.numeric(lookup[,2])))

# prepare header output file
cat("iter", "h2", "mu0", "mu1", "mu2", "af", "selcoeff", "selcoeff_simul", "beta_inter", "beta_es", "beta_esHOvsWT", selcoeff_bins, "\n", sep="\t", file=output_file)

cat("iter", "h2", "mu0", "mu1", "mu2", "af", "selcoeff", "selcoeff_simul", "beta_inter", "beta_es", "beta_esHOvsWT", selcoeff_bins, "\n", sep="\t", file=paste(output_file, "_ld", sep="", collapse=""))

# for each iteration
for (n in 1:niter) {

	k=0;

	while(k<K) {

		if (selcoeff>0) {
			# pick a selcoeff from a very skewed normal distribution
			nselcoeff=rnorm(1, selcoeff, selcoeff/100);
			#nselcoeff=selcoeff;
		} else nselcoeff=0;

		## SIMULATE genotypes according to a certain selection coefficient and current allele frequency

		# standard constant-size model
		# e.g. Ne 1e4
		# len of 500k
		# mut rate 1.1e-8 per site per generation (Roach et al Science 2010) -> theta 4NuL = 220
		# recomb rate 1.26 cM/Mb (Jensen-Seaman et al Genome research 2004) -> rho 0.000504 

		cmd=paste(msms_dir, " -N ", Ne, " -ms ", nchroms, " 1 -t ", theta, " -r ", rho, " ", format(nsites,digits=), " -SAA ", nselcoeff, " -SAa ", nselcoeff/2, " -Sp 0.50 -SF 0 ", currentfreq, " -Smark > out2.msms", sep="", collapse="");
		system(cmd);
#		system(paste("sed -i \"\" '1s/.*/-ms ", nchroms, " 1 -t ", theta, " -r ", rho, " ", nsites, "/' out2.msms", sep="", collapse=""));

		## ESTIMATE selection coefficient likelihood surface from nSL and a lookup table

		#system(paste(nsl_dir, " -msfile out2.msms -maxLen ", maxlen, " -msLen ", nsites, " 2> out2.log.nsl > out2.nsl", sep="", collapse=""))
		system(paste("python ms2nSL.py ", nsites, sep=""))
		cmd=paste(nsl_dir, 'nSL -samfile sample.txt -hapfile Haplotype.txt -adfile Ans_Der.txt -maxLen ', nsites,' > out2.nsl',sep="")
		system(cmd)
		system('rm sample.txt Haplotype.txt Ans_Der.txt')

		values=read.table("out2.nsl", header=T, stringsAsFactors=F)[,2:3]
		
		# get value for selected site (in the middle)
		midpos=nsites/2
		ind=which(values[,1]==midpos); 
		if (length(ind)==1) {
			nsl=values[ind,2];
			sc_like=dnorm(nsl, lookup[,4], lookup[,5])
			#the probability of the nSL being in each selection coefficient category is proportinal to the probability density 
			#added divide by the sum: Chenling03032016
			sc_like=sc_like/sum(sc_like)
			k=k+1;

			# also get mean nSL values for SNPs in LD with target (r2>0.8), here 0.7
			LD_r2_th=0.7
			haps=readMs(filein="out2.msms", nsam=nchroms, len=nsites)
			lin=linkage(haps$hap[[1]], which(haps$pos[[1]]==midpos))[[1]]
			nsl_lds=abs(c(values[match(haps$pos[[1]][which(lin>LD_r2_th)], values[,1]),2]))
			if (nsl<0) nsl_lds=-nsl_lds
			nsl_ld=mean(nsl_lds, na.rm=T)
			sc_like_ld=dnorm(nsl_ld, lookup[,4], (lookup[,5]))
			#same, divide by the sum: Chenling03032016
			sc_like_ld=sc_like_ld/(sum(sc_like_ld))
			
			## RETRIEVE  allele frequencies and genotypes from simulated data

			ms=readLines("out2.msms"); 
			pos=strsplit(ms[6], split=" ")[[1]];
			pos=as.numeric(pos[2:length(pos)]);
			ms=ms[7:(length(ms)-1)]
			genos=matrix(as.numeric(unlist(strsplit(ms, split=""))),nrow=length(ms), byrow=T)[,which(pos==0.5)]
			genos=genos[seq(1, length(genos), 2)]+genos[seq(2, length(genos), 2)];
			p=sum(genos)/length(ms) # -> allele frequency
			#rm(pos); rm(ms);
			q=1-p;

			## ASSIGN PHENOTYPES

			# additive model

			# INPUT h^2, hereditability
			# INPUT mean phenotypes and allelic effect

			#af should be calculated 
			# from genotype data for each site compute Va=2pqa^2 where a is the allelic effect (e.g. 10)
			#chenling03032016
			Emu=sum(c(p^2, 2*p*q, q^2)*mu)
			af=(mu[1]-mu[2])*p + mu[2]-Emu
			Va=2*p*q*(af^2)

			# h^2 = Va/Vp so compute Vp=Va/h^2
			Vp=Va/h2
			#Vp is additive phenotypic variation: chenling03032016
			#commented out the parts calculating the actual mean phenotype 
			# Vp = sig^2 + mu^2_bar  -mu_bar^2 , solve for sig (this is the variance in the gassians); mu^2_bar is the mean of squared means; mu_bar^2 is the square of the mean phenotypic value
			
			geno.props <- c((1-p)^2,(2*p*(1-p)),p^2); # or use observed proportions?
			mu.bar.sq <- sum(mu*geno.props)^2;
			mu.sq.bar <- sum(geno.props*(mu^2));
			ph.var <- Vp-mu.sq.bar+mu.bar.sq;
			#ph.var.all <- mu.sq.bar-mu.bar.sq
			# then for each individual (depending ot its genotype) pick up a phenotype accoridng to the corresponding normal distribution

			phenos <- NULL;
			for(d in 1:length(genos)) {
				if(genos[d]==0) {
					#phenos[d] <- rnorm(1,mean=mu[1],sd=sqrt(ph.var));
					phenos[d] <- rnorm(1,mean=mu[1],sd=sqrt(Vp));
				}else if(genos[d]==1){
					#phenos[d] <- rnorm(1,mean=mu[2],sd=sqrt(ph.var));
					phenos[d] <- rnorm(1,mean=mu[2],sd=sqrt(Vp));
        		}else if(genos[d]==2){
                	#phenos[d] <- rnorm(1,mean=mu[3],sd=sqrt(ph.var));
                	phenos[d] <- rnorm(1,mean=mu[3],sd=sqrt(Vp));
        		}
			}

			# COMPUTE effect size (code given by Anders Albrechtsen)
			# calculating the ML effect size assuming normal distributed errors (least squares)

			if (model==0) {
					# additive model
					X<-cbind(1,genos)
					beta<-solve(t(X)%*%X)%*%t(X)%*%phenos # beta[1,1] is intercept and beta[2,1] is the effect size per allele
				}

			if (model==1) {
				cat("\n not implemented. stop"); stop;
                        	# dominant model
                        	X<-cbind(1,ifelse(genos==0,0,1))
                        	beta<-solve(t(X)%*%X)%*%t(X)%*%phenos # beta[1,1] is intercept and beta[2,1] is the effect size
                	}

	     		if (model==2) {
				cat("\n not implemented. stop"); stop;
                        	# full model
				X<-cbind(1,ifelse(genos==1,1,0),ifelse(genos==2,1,0))
                        	beta<-solve(t(X)%*%X)%*%t(X)%*%phenos # beta[1,1] is intercept and beta[2,1] is the effect size HE vs WT and beta[3,1] is the effect size HO vs WT
				beta3=beta[3,1]
                	} else beta3=NA;

			# PRINT results
			#c(n, h2, mu, af, selcoeff, nselcoeff, beta[1,1], beta[2,1], beta3, sc_like)
			#c(n, h2, mu, af, selcoeff, nselcoeff, beta[1,1], beta[2,1], beta3, sc_like_ld)

			cat(n, h2, mu, af, selcoeff, nselcoeff, beta[1,1], beta[2,1], beta3, sc_like, "\n", sep="\t", file=output_file, append=T)
			cat(n, h2, mu, af, selcoeff, nselcoeff, beta[1,1], beta[2,1], beta3, sc_like_ld, "\n", sep="\t", file=paste(output_file,"_ld",sep="",collapse=""), append=T)

		} # end if nsl value

	} # for K loci

# do this for K loci and record disribution of selcoeff and single effect sizes -> estimate MLE of I -> use this I to re-estimate selcoeff from e*I -> measure accuracy

} # for n iterations

system("rm out*")




