
rm(list=ls())
library("limma")

# load data
fname_X <- 'normalized_intensity.txt'
X <- as.matrix(read.table(fname_X, header=F, row.names=NULL, sep="\t", check.names=F))
X2 <- X
X2[X2==0] <- NA
X2 <- log2(X2)

fname_Y <- 'sampleInfo.txt'
Y <- as.matrix(read.table(fname_Y, header=T, row.names=NULL, sep="\t", check.names=F))


# select samples
drug="Clozapine"
# drug="Imipramine"
# drug="Flouxetine"
ind_cont <- which(Y[,2] %in% "Control")
ind_drug <- which(Y[,2] %in% drug)
ind <- c(ind_cont,ind_drug)
X2 <- X2[,ind]
Y <- Y[ind,]


# build the design matrix
group <- factor(Y[,2])
batch <- factor(Y[,1])
design <- model.matrix(~0+group+batch)


# fit the model
fit <- lmFit(X2, design)
cont.matrix <- makeContrasts(diseasevscontrol=paste("group",drug,"-groupControl",sep=''), levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
result <- topTable(fit2, number="all", adjust.method = "BH", sort.by="none")
write.table(result, paste0("adjbatchlimma_",drug,".txt"), sep="\t", col.names=NA)

