
ml <- function(S, es, sc, h2, denom=10, cats) {
        nloci=length(es);
        if (nloci!=nrow(sc)) stop("Dimensions do not match.");
        ncat=ncol(sc);

        value=0;

        for (i in 1:nloci) {

                subvalue=0;
                esI_mean=es[i]*(S/h2);
		esI_sd= (es[i]/denom) * (S/h2);
                #esI=dnorm(cats, esI_mean, esI_sd);
                width=(cats[2]-cats[1])/2;
                esI=pnorm((cats+width),esI_mean, esI_sd)-pnorm((cats-width),esI_mean, esI_sd);
                for (j in 1:ncat) {
                        subvalue=subvalue+(sc[i,j]*esI[j]);
                }
                value = value + log10(subvalue);
        }
        -value;
}

# unknown: S/h2
mlr <- function(rat, es, sc, denom=10, cats) {

        nloci=length(es);
        if (nloci!=nrow(sc)) stop("Dimensions do not match.");
        ncat=ncol(sc);

        value=0;

        for (i in 1:nloci) {

                subvalue=0;
                esI_mean=es[i]*(rat);
                esI_sd= (es[i]/denom) * (rat);
                esI=dnorm(cats, esI_mean, esI_sd);

                for (j in 1:ncat) {
                        subvalue=subvalue+(sc[i,j]*esI[j]);
                }
                value = value + log10(subvalue);
        }
        -value;
}

# unknown: S and h2
ml2 <- function(par, es, sc, denom=10, cats) {

	S=par[1]
	h2=par[2]

        nloci=length(es);
        if (nloci!=nrow(sc)) stop("Dimensions do not match.");
        ncat=ncol(sc);

        value=0;

        for (i in 1:nloci) {

                subvalue=0;
                esI_mean=es[i]*(S/h2);
                esI_sd= (es[i]/denom) * (S/h2);
                esI=dnorm(cats, esI_mean, esI_sd);

                for (j in 1:ncat) {
                        subvalue=subvalue+(sc[i,j]*esI[j]);
                }
                value = value + log10(subvalue);
        }
	cat("\n",S,h2,-value)
        -value;
}




ml_sd <- function(S, es, sc, h2, es_sd, cats) {
        nloci=length(es);
        if (nloci!=nrow(sc)) stop("Dimensions do not match.");
        ncat=ncol(sc);
        value=0;
        for (i in 1:nloci) {
                subvalue=0;
                esI_mean=es[i]*(S/h2);
                esI_sd= es_sd[i] * (S/h2);
                esI=dnorm(cats, esI_mean, esI_sd);
                for (j in 1:ncat) {
                        subvalue=subvalue+(sc[i,j]*esI[j]);
                }
                value = value + log10(subvalue);
        }
        -value;
}



