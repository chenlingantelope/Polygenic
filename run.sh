
## This is the main file with the entire pipeline.

## 1) The first step is to create a look-up table, where we are tabulating mean and sd of nSL (our primary summary statistics for selection) for different scenarios. 

# This script takes two files: the first one contains the parameters for the simulations while the second one is the output file name.
# I provide you with a sample parameter file with all values that need to be set up. Open it (params_lookup.txt):
#NCHROMS 100,200 # how many chromosomes to simulate, if multiple choices then they should be comma-separated
#NE      10000 # effective pop size (Ne), for now leave it at 10k
#NITER   1000 # how many iterations per scenario (per combination of sel coeff and current allele frequency, see below)
#SELCOEFFS       20,1200,20 # intervals (start,end,set) of selection coefficients to simulate in 2*Ne units
#CURRFREQS       0.1,0.9,0.2 # intervals (start,end,set) of current allele frequency of the selected variant to simulate
#MSMS    ~/Documents/Software/msms/bin/msms # path to msms software
#NSL     ~/Documents/Software/nSL/nSL_matteo # path to nSL software (please use the one provided, as I changed something)
#TEMPFILE        out # name for temporary file, leave it like that
#THETA   88 # population-scaled mutation rate = 4*Ne*mu*L with mu mutation rate per site per generation and L length of simulated region (in bp)
#RHO     101 # recombination rate in 4*Ne*r*(L-1)
#NSITES  200000 # length of simulated region (in bp)
#VERBOSE 1 # prints out messages while running or not
#MAXLEN  200 # option for nSL, maximum length of the window, you may want to keep it like that at the beginning

# Once you create and edit this file you can run this R script.

Rscript createLookUp.R params_lookup.txt lookup.txt

# Results are stored in the output file as: number of chromosomes, selection coefficient, current allele frequency, mean nSL, sd nSL

## 2) The second step is to interpolate the distribution of mean and sd of nSL against the selection coefficient.
# This can be achieved by the following script. However, first you need to define how to interpolate.

DEGREE_POLY=1
# if positive (e.g. 1,2,3) then it uses a polynomial with DEGREE_POLY degrees; if 0 it uses an exponential; if -1 it uses a mixture of exponential
# the exponential fitting hasn't been fully tested.

Rscript interpolate.R lookup.txt $DEGREE_POLY lookup_inter.txt inter_log.txt

# This script produces 2 files: the first one (here called lookup_inter.txt) is equal to the look-up table file but values are not taken from the fitted curve; the log file (the second one) reports fitting coefficients.

## 3) Now we can run the main part of the pipeline. This involves simulating genotypes-phenotypes and estimate S, the selection differential.

# First, we need to simulated how many iterations per scenario (allelic effect, selection coefficients and so on, see below) we want to simulate and test.
ITER=50

> results.txt # initialise the output file

# Iterate across all possible parameters

for NC in 100 200 # nr of chromosomes
do
	for K in 20 50 # nr of loci
	do
		for CF in 0.7 0.8 0.9 # current alle freqs
		do
			for H in 0.6 0.75 0.9 # hereditability (h2)
			do
				for AF in 5 10 20 # effect size
				do
					for M in 0 # model, 0 is additive, only available now
					do
						for SC in 30 50 70 100 200 300 400 600 # sel coeffs to test
						do

							# initialize intermediate files
							> sc_es.txt
							> sc_es.txt_ld

							# this script will simulate genotypes and phenotypes
							Rscript simulate.R $K $NC $SC $CF $H 0 $((0+$AF)) $((0+$AF+$AF)) $AF $M lookup_inter.txt sc_es.txt $ITER params_lookup.txt >> results.txt
							# it works like that (for each iteration):
							# - pick a sel coeff from a normal distribution N(s,s/100) (you can change this)
							# - simulate genotypes according to the sampled sel coeff and current allele frequency
							# - estimate sel coeff likelihood surface from nSL and a lookup table
							# - in addition, take the mean nSL for SNPs in LD with target one (these results will be in the _ld file), this is done to assume the case you don't type the true casual variant
							# - assign phenotypes based on genotype and model, h2, allelic effect
							# - again phenotypes are taken from a norma distribution
							# - compute actual effect size

							# the output files (sc_es.*) have the following columns:
							# n, h2, mu, af, selcoeff, nselcoeff, beta[1,1], beta[2,1], beta3, sc_like (likelihood surface)

							# ML estimate of S
							# this script takes as input the file previously generated, while 0/1 is boolean for verbosity
							Rscript mleS.R sc_es.txt 0 >> results.txt
							Rscript mleS.R sc_es.txt_ld 0 >> results.txt
							# it will print: mean hat{S}, std hat{S}, std bias, rmsd, std bias using median, rmsd using median, + the same assuming known sel coeffs
							# so basically we test the method by re-estimating sel coeff assuming the estimated S and compute bias

							# all results are appended to this file
							echo >> results.txt # add new line

							# rm sc_es.* # uncomment if you want to clean up
						done
					done
				done
			done
		done
	done
done

## 4) Finally, you can plot some results.
# You may need to edit this file to save the plots in the desired directory
Rscript plotRes.R lookup.txt results.txt




