
# this script accepts a lookup table file, and interpolate points with a curve

args = commandArgs(trailingOnly = TRUE);
look_file = args[1]; # lookup table file
degree_poly = as.numeric(args[2]); # degree of polynomial
out_file = args[3]; # output file
log_file = args[4]; # log interpolation

# if degree is 0 then use exponential

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

		if (degree_poly>0) { # polynomial

			fit_mean = lm( y1 ~ poly(x, degree_poly) )
			fit_sd = lm( y2 ~ poly(x, degree_poly) )
			cat(c(i, j, as.numeric(fit_mean$coefficients), as.numeric(fit_sd$coefficients ) ), sep="\t", file=log_file, append=T )
			cat("\n", file=log_file, append=T)

			y1_pred=as.numeric(predict(fit_mean, data.frame(x=x)))
			y2_pred=as.numeric(predict(fit_sd, data.frame(x=x)))

		} else { # mixture of exponential

			data=data.frame(y1,y2,x)

			if (degree_poly==0) { # 1 exponential

				fit_mean  = nls( (y1+2) ~ exp(a+b*x), data=data, start= list(a = 0, b = 0))
				fit_sd  = nls( y2 ~ exp(a+b*x), data=data, start= list(a = 0, b = 0))

			}

			if (degree_poly==(-1)) { # mix of 2 exponentials
                                fit_mean  = nls( (y1+2) ~ exp(a+b*x)*w + exp(c+d*x)*(1-w), data=data, start= list(a = 0, b = -0.5, c=0, d=-0.5, w=0.5))
                                fit_sd  = nls( y2 ~ exp(a+b*x)*w + exp(c+d*x)*(1-w), data=data, start= list(a = 0, b = 0, c=0, d=0, w=0.5))

                        }

			y1_pred=as.numeric(predict(fit_mean, data.frame(x=x)))-2
                        y2_pred=as.numeric(predict(fit_sd, data.frame(x=x)))

			#nls(y1 ~ (a*(exp(1)^(-a*x))*w + b*(exp(1)^(-b*x))*(1-w)), data=data, start=list(a=1,b=1,w=0.5), trace=T)

			# old try
			# mle=optim( par=c(1,1,0.5), fn=likeF, method="L-BFGS-B", lower=c(0,0,0), upper=c(10,10,1), v=y1 )


		}

		# write to file
		for (b in 1:length(y1_pred)) cat(i, x[b], j, y1_pred[b], y2_pred[b], sep="\t", "\n", file=out_file, append=T)

	}
}

# alternatively I can use b-splines with R function 'ns'

