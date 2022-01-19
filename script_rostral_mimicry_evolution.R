
###### Evolution of rostral scale and mimicry in the genus Xenodon Boie, 1826 (Serpentes, Dipsadidae, Xenodontinae) #####
################################# Hugo Cabral, Pier Cacciali & Diego J. Santana #########################################
# Library
library(caper)
library(phytools)
library(geiger)
library(ade4)
library(picante)
library(phytools)
library(caper)
library(ape)
rm(list=ls()); gc()
setwd("")

# Load table with information of characters
X<-read.csv("",row.names=1, sep=";")
rostral<-setNames(X[,1],rownames(X))
mimic<-setNames(X[,2], rownames(X))
habitat <-setNames(X[,4], rownames(X))

# Load tree
Tree<-read.tree("")
Tree
head(xenodon_out)
plot(Tree, lwd=4)
axisPhylo()

####ANCESTRAL CHARACTER RECONSTRUCTION####
# Perform the ancestral reconstruction of rostral and mimicry with fitDiscrate. Here you can changes between rostral or mimicry. 
# In $data the next numbers is the column in or traits table, would depend on the location of or trait##
# When data are binary, just use ER and ARD model, as the case of rostral, in the case of mimicry use the three models
xen_geiger <- treedata(Tree, X)
xen_geiger$phy
xen_geiger$data

xengeigerer <- fitDiscrete(xenodon, xen_geiger$data [,1], type = "discrete", model = "ER")
xengeigerer

xengeigersym <- fitDiscrete(xenodon, xen_geiger$data [,1], type = "discrete", model = "SYM")
xengeigersym

xengeigerard <- fitDiscrete(xenodon, xen_geiger$data [,1], type = "discrete", model = "ARD")
xengeigerard

# Model Selection
modSel.geiger(xengeigerer, xengeigersym, xengeigerard, type="AICc")

ER_SYM <- abs(xengeigerer$opt$aicc - xengeigersym$opt$aicc)
ER_SYM

ER_ARD <- abs(xengeigerer$opt$aicc - xengeigerard$opt$aicc)
ER_ARD

SYM_ARD <- abs(xengeigersym$opt$aicc - xengeigerard$opt$aicc)
SYM_ARD


dER_ard <- abs(2*(xengeigerer$opt$lnL-xengeigerard$opt$lnL))
dER_ard

pvalueER_ard <- pchisq(dER_ard, 2-1, lower.tail = FALSE)
pvalueER_ard

# Phylogenetic signal  for discrate characters####


# First, a lambda tree with no signal###
tree_lambda_0 <- rescale(Tree, model = "lambda", 0)

# Fitting lambda transformation to determinate if there is phylogenetic signal####

xen_lambda_geiger <- fitDiscrete(xenodon, xen_geiger$data [,1], type = "discrete", model = "ER", transform = "lambda", niter = 500)
xen_lambda_geiger

xen_lambda_geiger_0 <- fitDiscrete(xendon_tree_lambda_0, xen_geiger$data [,1], type="discrete", model = "ER", transform = "lambda", niter = 500)
xen_lambda_geiger_0

par(mfrow=c(1,2))
plot(Tree, edge.width = 1, show.tip.label = TRUE)
add.scale.bar(cex = 0.7, font = 2, col = "red")
plot(Tree_lambda_0, edge.width = 1, show.tip.label = FALSE)
add.scale.bar(cex = 0.7, font = 2, col= "red")
plot(xen_lambda_geiger_0, edge.width = 1, show.tip.label = FALSE)

# Model comparision for lambda, lowest aic best model##

xen_lambda_geiger_0_aicc <- abs(xen_lambda_geiger_0$opt$aicc - xen_lambda_geiger$opt$aicc)
xen_lambda_geiger_0_aicc

d_xen_lambda_geiger_0 <- abs(2*(xen_lambda_geiger_0$opt$lnL - xen_lambda_geiger$opt$lnL))
d_xen_lambda_geiger_0

p_d_xen_lambda_geiger_0 <- pchisq(d_xen_lambda_geiger_0, 1, lower.tail = FALSE)
p_d_xen_lambda_geiger_0


### MCMC stochastic character mapping. To perform the mimicry mapping, changes rostral to mimicry, following the original table

cols<- c("black", "#DF536B")
mtrees<-make.simmap(Tree,rostral,model="ER")
mtrees
plot(mtrees)
plot(mtrees,type="phylogram",fsize=0.8,ftype="i")

add.simmap.legend(colors=cols2,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(xenodon)),fsize=0.8)

mtrees2<-make.simmap(xenodon,rostral,model="ER",nsim=100)

par(mfrow=c(10,10))
null<-sapply(mtrees2,plot,lwd=1,ftype="off")
pd<-summary(mtrees2,plot=FALSE)
pd
plot(pd,fsize=0.6,ftype="i")

plot(mtrees2[[1]],type="phylogram",fsize=0.8,ftype="i", lwd=4)

nodelabels(pie=pd$ace,piecol=cols,cex=0.5)
add.simmap.legend(colors=cols2,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(xenodon)),fsize=0.8)
axisPhylo()

# Other options of representation

plot(fitER$lik.anc,pd$ace,xlab="marginal ancestral states",
     ylab="posterior probabilities from stochastic mapping")
lines(c(0,1),c(0,1),lty="dashed",col="red",lwd=2)

# Threshold model

# Testing correlation between rostral and habitat
matrix <- read.csv2("binary_matrix.csv", row.names = 1) #load Binary Matrix
hab.ros<-matrix[,c(1,2)]
hab.ros<-as.matrix(hab.ros)
da<- match.phylo.data(xenodon, hab.ros)
da$data

# set some parameters for the MCMC
sample <- 100
ngen <- 1000000
burnin <- 0.2 * ngen


# Correlation between rostral and habitat, This procedure would take at least 8 hours, depending on your computer specifications,
# and the number of generations ngen
hab.ros <- threshBayes(xenodon, da$data, types = c("disc", "disc"),
                       ngen = 1000000,control=list(sample=sample))

plot(hab_ros2)
d <- density(hab_ros2)
dev.off()
print(hab_ros2)
plot(d)
title(main=paste("posterior density of correlation",
                 "coefficient, r,\nfrom threshold model"),
      font.main=3)
# what is our estimate of the correlation?
mean(hab_ros2$par[(burnin/sample + 1):nrow(hab_ros2$par), "r"])

# plot our likelihood profile
plot(hab_ros2$par[, "logL"], type = "l", xlab = "generation", 
     ylab = "logL")

# plot our posterior sample for the correlation
h<-hist(hab_ros2$par[(burnin/sample + 1):nrow(hab_ros2$par), "r"],
        xlab = "r", ylab = "frequency",  main = NULL)
lines(c(0.03,0.03),c(0,max(h$density)),lty="dashed")
plot(density(hab_ros2$par[(burnin/sample+1):nrow(hab_ros2$par),
                          "r"],bw=0.1),xlab="r",main="posterior density for r")
lines(c(r,r),c(0,1000),lty="dashed")
