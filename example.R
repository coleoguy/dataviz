# simple seq def example code
library(ape)
library(viridis)
source("function.R")
######### Simple Example Presented ##############
# tree <- read.tree(text="(taxa1:100, (taxa2:90,(taxa3:80,(taxa4:70, (taxa5:60, (taxa6:50,(taxa7:40,(taxa8:30, (taxa9:10, taxa10:10):20):10):10):10):10):10):10):10):10;")
# plot(tree, cex=.8)
# df <- read.csv("seqdef.csv", header=F)
# result <- SeqDef(tree, df, 2, invert=T)
# The SeqDef Function has 5 arguments
## tree: a phylogentic tree in phylo format
## df: a dataframe that has a column that matches the species names on the 
##     phylogeny and at least one column of data that you wish to analyze
## data.col: a numeric variable specifying the column that contains 
##           data to be analyzed
## invert: TRUE or FALSE to indicate whether tips with the least data should be
##         1 (FALSE) or 0 (TRUE)
## scale: Rescale values to be on a 0 to 1 scale defaults to TRUE. Use of 
##        unscaled values have not been explored.

# Results from this simple example shows us that we see the pattern that we expect where 
# taxa have a higher sequencing deficiency statistic as we move further
# from the species that has sequencing data.
######## End of simple example presented #########

# Now we also want to show performance on a more 
# realistically complex phylogeny and 
# data distribution

# simulate a phylogeny under a Yule model
set.seed(3)
tree <- rcoal(50)
# pick two species to have sequenced genomes
# while the rest have no full genome sequences
dat <- as.data.frame(matrix(data=0,nrow=50,ncol=2))
colnames(dat) <- c("species","genome.assembly")
dat$species <- tree$tip.label
dat$genome.assembly[c(10,20)] <- 1
# make pest status with most having low values
x <- rexp(10000, rate = 3)
x <- sample(x[x<=1], 50)
dat$pest <- x
dat$pest[c(21,49)] <- 1

# Plot the base tree
par(mar = c(1, 1, 2, 1))
plot(tree,show.tip.label=F, x.lim = c(0, 3))
tiplabels(tip=which(dat[,2]==1),pch=16,adj=.52,cex=.75)
points(x=.01,y=48,pch=16)
text(x=.01,y=48,"genome assemblies",pos=4,cex=.6)

# calculate basic seqdef
res <- SeqDef(tree=tree, df=dat, data.col=2)
cols <- heat.colors(101)[round(res$seqdef*100+1)]
tiplabels(pch=21,adj=.62, bg=cols,cex=.75)
text(x=2.33,y=51.3,"SeqDef",cex=.6)
points(x=seq(from=1,to=0,length.out=101),y=rep(44.5,101),
       pch=15,col=heat.colors(101))
text(x=1,y=44.5,"no data", cex=.6,pos=1)
text(x=0.1,y=44.5,"complete data", cex=.6,pos=1)

# incorporate traits
tiplabels(pch=16,adj=.90,cex=.75, col=viridis(101)[round(dat$pest*100)+1])
text(x=2.65,y=51.3,"pest",cex=.6)
points(x=seq(from=0,to=1,length.out=101),y=rep(39,101),
       pch=15,col=viridis(101))
text(x=0.01,y=39,"nonpest", cex=.6,pos=1)
text(x=1,y=39,"pest", cex=.6,pos=1)

# calculate final priority
importance <- (1-res$seqdef) * dat$pest
priority <- round(100*importance)+1
tiplabels(pch=21,adj=1.2, 
          bg=viridis(101, option="C")[priority],cex=.75)
text(x=2.93,y=51.3,"priority",cex=.6)
points(x=seq(from=0,to=1,length.out=101),y=rep(34,101),
       pch=15,col=viridis(101,option="C"))
text(x=0.05,y=34,"low priority", cex=.6,pos=1)
text(x=1,y=34,"high priority", cex=.6,pos=1)


