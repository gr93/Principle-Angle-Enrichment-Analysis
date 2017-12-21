### R code from vignette source 'GSAR.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: GSAR.Rnw:116-117
###################################################
library(GSAR)


###################################################
### code chunk number 2: GSAR.Rnw:164-165
###################################################
options(width=80, digits=3)


###################################################
### code chunk number 3: GSAR.Rnw:168-179
###################################################
library(MASS)
set.seed(123)
nf <- 20
nobs <- 60
zero_vector <- array(0,c(1,nf))
cov_mtrx <- diag(nf)
dataset <- mvrnorm(nobs, zero_vector, cov_mtrx)
Wmat <- as.matrix(dist(dataset, method="euclidean", diag=TRUE, 
upper=TRUE, p=2))
gr <- graph_from_adjacency_matrix(Wmat, weighted=TRUE, mode="undirected")
MST <- mst(gr)


###################################################
### code chunk number 4: GSAR.Rnw:209-214
###################################################
## The input of findMST2 must be a matrix with rows and columns 
## respectively corresponding to genes and columns.
## Therefore, dataset must be transposed first.
dataset <- aperm(dataset, c(2,1))
MST2 <- findMST2(dataset)


###################################################
### code chunk number 5: GSAR.Rnw:314-315
###################################################
HDP.ranking(MST)


###################################################
### code chunk number 6: GSAR.Rnw:318-319
###################################################
radial.ranking(MST)


###################################################
### code chunk number 7: mst
###################################################
V(MST)$color <- c(rep("green",20), rep("yellow",20))
par(mfrow=c(1,1), mar=c(1,1,1,1), oma=c(1,1,1,1))
plot(MST, vertex.label.cex=0.8, vertex.size=9)


###################################################
### code chunk number 8: GSAR.Rnw:522-525
###################################################
library(GSVAdata)
data(p53DataSet)
data(c2BroadSets)


###################################################
### code chunk number 9: GSAR.Rnw:539-576
###################################################
library(org.Hs.eg.db)
library(GSEABase)
C2 <- as.list(geneIds(c2BroadSets))
len <- length(C2)
genes.entrez <- unique(unlist(C2))
genes.symbol <- array("",c(length(genes.entrez),1))
x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])
for (ind in 1:length(genes.entrez)){
    if (length(xx[[genes.entrez[ind]]])!=0)
        genes.symbol[ind] <- xx[[genes.entrez[ind]]]
                                   }
## discard genes with no mapping to gene symbol identifiers
genes.no.mapping <- which(genes.symbol == "")
if(length(genes.no.mapping) > 0){
    genes.entrez <- genes.entrez[-genes.no.mapping]
    genes.symbol <- genes.symbol[-genes.no.mapping]
                                }
names(genes.symbol) <- genes.entrez
## discard genes in C2 pathways which do not exist in p53 dataset
p53genes <- rownames(p53DataSet)
remained <- array(0,c(1,len))
for (k in seq(1, len, by=1)) {
    remained[k] <- sum((genes.symbol[C2[[k]]] %in% p53genes) & 
    (C2[[k]] %in% genes.entrez))
                             }
## discard C2 pathways which have less than 10 or more than 500 genes
C2 <- C2[(remained>=10)&&(remained<=500)]
pathway.names <- names(C2)
c2.pathways <- list()
for (k in seq(1, length(C2), by=1)) {
    selected.genes <- which(p53genes %in% genes.symbol[C2[[k]]])
    c2.pathways[[length(c2.pathways)+1]] <- p53genes[selected.genes]
                                    }
names(c2.pathways) <- pathway.names
path.index <- which(names(c2.pathways) == "LU_TUMOR_VASCULATURE_UP")


###################################################
### code chunk number 10: GSAR.Rnw:585-601
###################################################
target.pathway <- p53DataSet[c2.pathways[["LU_TUMOR_VASCULATURE_UP"]],]
group.label <- c(rep(1,17), rep(2,33))
WW_pvalue <- WWtest(target.pathway, group.label)
KS_pvalue <- KStest(target.pathway, group.label)
MD_pvalue <- MDtest(target.pathway, group.label)
RKS_pvalue <- RKStest(target.pathway, group.label)
RMD_pvalue <- RMDtest(target.pathway, group.label)
F_pvalue <- AggrFtest(target.pathway, group.label)
GSNCA_pvalue <- GSNCAtest(target.pathway, group.label)
WW_pvalue
KS_pvalue
MD_pvalue
RKS_pvalue
RMD_pvalue
F_pvalue
GSNCA_pvalue


###################################################
### code chunk number 11: mst2plot
###################################################
plotMST2.pathway(object=p53DataSet[c2.pathways[[path.index]],],
group=c(rep(1,17), rep(2,33)), name="LU_TUMOR_VASCULATURE_UP", 
legend.size=1.2, leg.x=-1.2, leg.y=2, 
label.size=1, label.dist=0.8, cor.method="pearson")


###################################################
### code chunk number 12: GSAR.Rnw:647-650
###################################################
results <- TestGeneSets(object=p53DataSet, group=group.label, 
geneSets=c2.pathways[1:3], min.size=10, max.size=100, test="GSNCAtest")
results


###################################################
### code chunk number 13: GSAR.Rnw:677-712
###################################################
library(Biobase)
library(genefilter)
library(annotate)
library(hgu95av2.db)
library(ALL)
data(ALL)
bcell = grep("^B", as.character(ALL$BT))
types = c("NEG", "BCR/ABL")
moltyp = which(as.character(ALL$mol.biol) %in% types)
ALL_bcrneg = ALL[, intersect(bcell, moltyp)]
ALL_bcrneg$mol.biol = factor(ALL_bcrneg$mol.biol)
ALL_bcrneg$BT = factor(ALL_bcrneg$BT)
nBCR <- sum(ALL_bcrneg$mol.biol == "BCR/ABL")
nNEG <- sum(ALL_bcrneg$mol.biol == "NEG")
BCRsamples <- which(ALL_bcrneg$mol.biol == "BCR/ABL")
NEGsamples <- which(ALL_bcrneg$mol.biol == "NEG")
ALL_bcrneg <- ALL_bcrneg[,c(BCRsamples,NEGsamples)]
platform <- annotation(ALL_bcrneg)
annType <- c("db", "env")
entrezMap <- getAnnMap("ENTREZID", annotation(ALL_bcrneg), 
type=annType, load=TRUE)
symbolMap <- getAnnMap("SYMBOL", annotation(ALL_bcrneg), 
type=annType, load=TRUE)
filtered <- nsFilter(ALL_bcrneg, require.entrez=TRUE, 
remove.dupEntrez=FALSE, require.symbol=TRUE, require.GOBP=FALSE, 
var.func=IQR, var.filter=FALSE, var.cutof=0.5)
filtered.set <- filtered$eset
probe.names <- featureNames(filtered.set)
rr <- rowttests(filtered.set, as.factor(ALL_bcrneg$mol.biol), tstatOnly=TRUE)
fL <- findLargest(probe.names, abs(rr$statistic), platform)
filtset2 <- filtered.set[fL,]
affymetrix.probe.names <- featureNames(filtset2)
gene.symbols <- lookUp(affymetrix.probe.names, platform, "SYMBOL")
featureNames(filtset2) <- gene.symbols
ALLdataset <- exprs(filtset2)


###################################################
### code chunk number 14: GSAR.Rnw:726-761
###################################################
C2 <- as.list(geneIds(c2BroadSets))
len <- length(C2)
genes.entrez <- unique(unlist(C2))
genes.symbol <- array("",c(length(genes.entrez),1))
x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])
for (ind in 1:length(genes.entrez)){
    if (length(xx[[genes.entrez[ind]]])!=0)
        genes.symbol[ind] <- xx[[genes.entrez[ind]]]
                                   }
## discard genes with no mapping to gene symbol identifiers
genes.no.mapping <- which(genes.symbol == "")
if(length(genes.no.mapping) > 0){
    genes.entrez <- genes.entrez[-genes.no.mapping]
    genes.symbol <- genes.symbol[-genes.no.mapping]
                                }
names(genes.symbol) <- genes.entrez
## discard genes in C2 pathways which do not exist in ALL dataset
ALLgenes <- rownames(ALLdataset)
remained <- array(0,c(1,len))
for (k in seq(1, len, by=1)) {
    remained[k] <- sum((genes.symbol[C2[[k]]] %in% ALLgenes) & 
    (C2[[k]] %in% genes.entrez))
                             }
## discard C2 pathways which have less than 10 or more than 500 genes
C2 <- C2[(remained>=10)&&(remained<=500)]
pathway.names <- names(C2)
c2.pathways <- list()
for (k in seq(1, length(C2), by=1)) {
    selected.genes <- which(ALLgenes %in% genes.symbol[C2[[k]]])
    c2.pathways[[length(c2.pathways)+1]] <- ALLgenes[selected.genes]
                                    }
names(c2.pathways) <- pathway.names
path.index <- which(names(c2.pathways) == "KEGG_CHRONIC_MYELOID_LEUKEMIA")


###################################################
### code chunk number 15: GSAR.Rnw:768-784
###################################################
KCMLpathway <- ALLdataset[c2.pathways[["KEGG_CHRONIC_MYELOID_LEUKEMIA"]],]
group.label <- c(rep(1,37), rep(2,42))
WW_pvalue <- WWtest(KCMLpathway, group.label)
KS_pvalue <- KStest(KCMLpathway, group.label)
MD_pvalue <- MDtest(KCMLpathway, group.label)
RKS_pvalue <- RKStest(KCMLpathway, group.label)
RMD_pvalue <- RMDtest(KCMLpathway, group.label)
F_pvalue <- AggrFtest(KCMLpathway, group.label)
GSNCA_pvalue <- GSNCAtest(KCMLpathway, group.label)
WW_pvalue
KS_pvalue
MD_pvalue
RKS_pvalue
RMD_pvalue
F_pvalue
GSNCA_pvalue


###################################################
### code chunk number 16: KCMLplot
###################################################
plotMST2.pathway(object=KCMLpathway, group=group.label, 
name="KEGG_CHRONIC_MYELOID_LEUKEMIA", legend.size=1.2, leg.x=-1, 
leg.y=2, label.size=0.8, cor.method="pearson")


###################################################
### code chunk number 17: GSAR.Rnw:837-848
###################################################
library(tweeDEseqCountData)
data(pickrell)
data(annotEnsembl63)
data(genderGenes)
gender <- pickrell.eset$gender
pickrell.eset
sampleNames(pickrell.eset)[gender == "male"]
sampleNames(pickrell.eset)[gender == "female"]
head(annotEnsembl63)
length(msYgenes)
length(XiEgenes)


###################################################
### code chunk number 18: GSAR.Rnw:855-858
###################################################
allXgenes <- rownames(annotEnsembl63)[annotEnsembl63$Chr == "X"]
Xigenes <- allXgenes[!(allXgenes %in% XiEgenes)]
length(Xigenes)


###################################################
### code chunk number 19: GSAR.Rnw:874-887
###################################################
library(edgeR)
gene.indices <- which(!(is.na(annotEnsembl63$EntrezID) | 
is.na(annotEnsembl63$Length)))
PickrellDataSet <- exprs(pickrell.eset)
PickrellDataSet <- PickrellDataSet[gene.indices,]
genes.length <- annotEnsembl63$Length[gene.indices]
cpm.matrix <- cpm(PickrellDataSet)
cpm.means <- rowMeans(cpm.matrix)
cpm.filter <- which(cpm.means > 0.1)
PickrellDataSet <- PickrellDataSet[cpm.filter,]
genes.length <- genes.length[cpm.filter]
rpkm.set <- rpkm(PickrellDataSet, genes.length)
rpkm.set <- log2(1 + rpkm.set)


###################################################
### code chunk number 20: GSAR.Rnw:898-905
###################################################
gene.space <- rownames(rpkm.set)
msYgenes <- msYgenes[msYgenes %in% gene.space]
XiEgenes <- XiEgenes[XiEgenes %in% gene.space]
Xigenes <- Xigenes[Xigenes %in% gene.space]
XYgenes <- c(msYgenes, XiEgenes)
length(XYgenes)
length(Xigenes)


###################################################
### code chunk number 21: GSAR.Rnw:911-917
###################################################
XYpathway <- rpkm.set[XYgenes,]
group.label.pickrell <- (gender == "male") + 1
WW_pvalue <- WWtest(XYpathway, group.label.pickrell)
KS_pvalue <- KStest(XYpathway, group.label.pickrell)
WW_pvalue
KS_pvalue


###################################################
### code chunk number 22: GSAR.Rnw:922-927
###################################################
Xipathway <- rpkm.set[Xigenes,]
WW_pvalue <- WWtest(Xipathway, group.label.pickrell)
KS_pvalue <- KStest(Xipathway, group.label.pickrell)
WW_pvalue
KS_pvalue


###################################################
### code chunk number 23: GSAR.Rnw:933-950
###################################################
nrow(XYpathway)
nrow(Xipathway)
tiny.sd.XY.female <- which(apply(XYpathway[, group.label.pickrell == 1], 1, "sd") < 1e-3)
tiny.sd.XY.male <- which(apply(XYpathway[, group.label.pickrell == 2], 1, "sd") < 1e-3)
tiny.sd.Xi.female <- which(apply(Xipathway[, group.label.pickrell == 1], 1, "sd") < 1e-3)
tiny.sd.Xi.male <- which(apply(Xipathway[, group.label.pickrell == 2], 1, "sd") < 1e-3)
length(tiny.sd.XY.female)
length(tiny.sd.XY.male)
length(tiny.sd.Xi.female)
length(tiny.sd.Xi.male)
apply(XYpathway[, group.label.pickrell == 1], 1, "sd")
if(length(tiny.sd.XY.male) > 0) XYpathway <- XYpathway[-tiny.sd.XY.male,]
if(length(tiny.sd.XY.female) > 0) XYpathway <- XYpathway[-tiny.sd.XY.female,]
if(length(tiny.sd.Xi.male) > 0) Xipathway <- Xipathway[-tiny.sd.Xi.male,]
if(length(tiny.sd.Xi.female) > 0) Xipathway <- Xipathway[-tiny.sd.Xi.female,]
nrow(XYpathway)
nrow(Xipathway)


###################################################
### code chunk number 24: GSAR.Rnw:973-974
###################################################
sessionInfo()


