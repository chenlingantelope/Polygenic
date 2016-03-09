
multiplot <- function(..., plotlist=NULL, cols) {
    require(grid)

    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)
    # Make the panel
    plotCols = cols                       # Number of columns of plots
    plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols

    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
    vplayout <- function(x, y)
        viewport(layout.pos.row = x, layout.pos.col = y)

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
        curRow = ceiling(i/plotCols)
        curCol = (i-1) %% plotCols + 1
        print(plots[[i]], vp = vplayout(curRow, curCol ))
    }

}


options(stringsAsFactors = FALSE)

library(ggplot2)

args=commandArgs(T)
lookfile=args[1]
infile=args[2]
rm(args)

# lookup

look=read.table(lookfile, fill=TRUE)
colnames(look)=c("nchroms","selcoeff","currentfreq","nSL_mean","nSL_sd")

qq=c()
cf=seq(0.1,0.9,0.1)
for (i in 1:length(cf)) {
	look.subset <- subset(look, currentfreq==cf[i])
	qq[[i]]=qplot(data=look.subset, x=selcoeff, y=nSL_mean) #+ facet_grid(nchroms ~ currentfreq, labeller=label_both)
}
multiplot(qq, cols=3)
ggsave("Figures/lookup_mean.eps")

qq=list()
cf=seq(0.1,0.9,0.1)
for (i in 1:length(cf)) {
        look.subset <- subset(look, currentfreq==cf[i])
        qq[[i]]=qplot(data=look.subset, x=selcoeff, y=nSL_sd) #+ facet_grid(nchroms ~ currentfreq, labeller=label_both)
}
multiplot(qq, cols=3)
ggsave("Figures/lookup_sd.eps")

# res

res=read.table(infile, header=F, fill=TRUE)

colnames(res)=c("nloci","nchroms","selcoeff","currentfreq","h2","muAA","muAa","muaa","allelic_effect","model","niter","S","S_std","bias","rmsd", "bias_median", "rmsd_median", "S_known","S_std_known", "bias_known","rmsd_known", "bias_median_known", "rmsd_median_known", "S_LD","S_std_LD","bias_LD","rmsd_LD", "bias_median_LD", "rmsd_median_LD", "S_known_LD","S_std_known_LD", "bias_known_LD","rmsd_known_LD", "bias_median_known_LD", "rmsd_median_known_LD")

res=as.data.frame(res)

res$h2=as.character(res$h2)
res$allelic_effect=as.character(res$allelic_effect)
res$nloci=as.character(res$nloci)
res$nchroms=as.character(res$nchroms)
res$currentfreq=as.character(res$currentfreq)

# plot

qplot(data=res[which(res$nchroms==100),], x=selcoeff, y=S, shape=nloci, colour=currentfreq) + facet_grid(allelic_effect ~ h2, labeller=label_both)
ggsave("Figures/S.eps")

qplot(data=res[which(res$nchroms==100),], x=selcoeff, y=S_known, shape=nloci, colour=currentfreq) + facet_grid(allelic_effect ~ h2, labeller=label_both)
ggsave("Figures/S_known.eps")

qplot(data=res[which(res$nchroms==100),], x=selcoeff, y=S_std, shape=nloci, colour=currentfreq) + facet_grid(allelic_effect ~ h2, labeller=label_both)
ggsave("Figures/S_std.eps")

qplot(data=res[which(res$nchroms==100),], x=selcoeff, y=S_std_known, shape=nloci, colour=currentfreq) + facet_grid(allelic_effect ~ h2, labeller=label_both)
ggsave("Figures/S_std_known.eps")

qplot(data=res[which(res$nchroms==100),], x=selcoeff, y=bias, shape=nloci, colour=currentfreq) + facet_grid(allelic_effect ~ h2, labeller=label_both)
ggsave("Figures/bias.eps")

qplot(data=res[which(res$nchroms==100),], x=selcoeff, y=bias_known, shape=nloci, colour=currentfreq) + facet_grid(allelic_effect ~ h2, labeller=label_both)
ggsave("Figures/bias_known.eps")

qplot(data=res[which(res$nchroms==100),], x=selcoeff, y=rmsd, shape=nloci, colour=currentfreq) + facet_grid(allelic_effect ~ h2, labeller=label_both)
ggsave("Figures/rmsd.eps")

qplot(data=res[which(res$nchroms==100),], x=selcoeff, y=rmsd_known, shape=nloci, colour=currentfreq) + facet_grid(allelic_effect ~ h2, labeller=label_both)
ggsave("Figures/rmsd_known.eps")

qplot(data=res[which(res$nchroms==100),], x=selcoeff, y=bias_median, shape=nloci, colour=currentfreq) + facet_grid(allelic_effect ~ h2, labeller=label_both)
ggsave("Figures/bias_median.eps")

qplot(data=res[which(res$nchroms==100),], x=selcoeff, y=bias_median_known, shape=nloci, colour=currentfreq) + facet_grid(allelic_effect ~ h2, labeller=label_both)
ggsave("Figures/bias_median_known.eps")

qplot(data=res[which(res$nchroms==100),], x=selcoeff, y=rmsd_median, shape=nloci, colour=currentfreq) + facet_grid(allelic_effect ~ h2, labeller=label_both)
ggsave("Figures/rmsd_median.eps")

qplot(data=res[which(res$nchroms==100),], x=selcoeff, y=rmsd_median_known, shape=nloci, colour=currentfreq) + facet_grid(allelic_effect ~ h2, labeller=label_both)
ggsave("Figures/rmsd_median_known.eps")

# LD

qplot(data=res[which(res$nchroms==100),], x=selcoeff, y=S_LD, shape=nloci, colour=currentfreq) + facet_grid(allelic_effect ~ h2, labeller=label_both)
ggsave("Figures/S_LD.eps")

qplot(data=res[which(res$nchroms==100),], x=selcoeff, y=S_known_LD, shape=nloci, colour=currentfreq) + facet_grid(allelic_effect ~ h2, labeller=label_both)
ggsave("Figures/S_known_LD.eps")

qplot(data=res[which(res$nchroms==100),], x=selcoeff, y=S_std_LD, shape=nloci, colour=currentfreq) + facet_grid(allelic_effect ~ h2, labeller=label_both)
ggsave("Figures/S_std_LD.eps")

qplot(data=res[which(res$nchroms==100),], x=selcoeff, y=S_std_known_LD, shape=nloci, colour=currentfreq) + facet_grid(allelic_effect ~ h2, labeller=label_both)
ggsave("Figures/S_std_known_LD.eps")

qplot(data=res[which(res$nchroms==100),], x=selcoeff, y=bias_LD, shape=nloci, colour=currentfreq) + facet_grid(allelic_effect ~ h2, labeller=label_both)
ggsave("Figures/bias_LD.eps")

qplot(data=res[which(res$nchroms==100),], x=selcoeff, y=bias_known_LD, shape=nloci, colour=currentfreq) + facet_grid(allelic_effect ~ h2, labeller=label_both)
ggsave("Figures/bias_known_LD.eps")

qplot(data=res[which(res$nchroms==100),], x=selcoeff, y=rmsd_LD, shape=nloci, colour=currentfreq) + facet_grid(allelic_effect ~ h2, labeller=label_both)
ggsave("Figures/rmsd_LD.eps")

qplot(data=res[which(res$nchroms==100),], x=selcoeff, y=rmsd_known_LD, shape=nloci, colour=currentfreq) + facet_grid(allelic_effect ~ h2, labeller=label_both)
ggsave("Figures/rmsd_known_LD.eps")

qplot(data=res[which(res$nchroms==100),], x=selcoeff, y=bias_median_LD, shape=nloci, colour=currentfreq) + facet_grid(allelic_effect ~ h2, labeller=label_both)
ggsave("Figures/bias_median_LD.eps")

qplot(data=res[which(res$nchroms==100),], x=selcoeff, y=bias_median_known_LD, shape=nloci, colour=currentfreq) + facet_grid(allelic_effect ~ h2, labeller=label_both)
ggsave("Figures/bias_median_known_LD.eps")

qplot(data=res[which(res$nchroms==100),], x=selcoeff, y=rmsd_median_LD, shape=nloci, colour=currentfreq) + facet_grid(allelic_effect ~ h2, labeller=label_both)
ggsave("Figures/rmsd_median_LD.eps")

qplot(data=res[which(res$nchroms==100),], x=selcoeff, y=rmsd_median_known_LD, shape=nloci, colour=currentfreq) + facet_grid(allelic_effect ~ h2, labeller=label_both)
ggsave("Figures/rmsd_median_known_LD.eps")






