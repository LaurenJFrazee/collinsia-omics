#### Prep DESeq2 package/data and preparation ####

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

#load
library(DESeq2)

#read in and format count data for analysis
CollinsiaRNAreadcountsonly <- read.delim("CollinsiaRNAreadcountsonlyERGO.txt", header = TRUE)
CollinsiaRNAreadcountsonlyDES <- CollinsiaRNAreadcountsonly[,-1]
row.names(CollinsiaRNAreadcountsonlyDES) <- CollinsiaRNAreadcountsonly$CollinsiaRNAreadcounts.annotation...c..GeneID...

#re-order & re-name sample columns for ease of interpretation
CollinsiaRNAreadcountsonlyDES <- CollinsiaRNAreadcountsonlyDES[,c(10,1:9,20,11:19)]
colnames(CollinsiaRNAreadcountsonlyDES) <- c("C. linearis all", "C. linearis Stage B1 plant 1",
                                              "C. linearis Stage B1 plant 2", "C. linearis Stage B1 plant 3",
                                              "C. linearis Stage B2 plant 1", "C. linearis Stage B2 plant 2",
                                              "C. linearis Stage B2 plant 3", "C. linearis Stage B3 plant 1",
                                              "C. linearis Stage B3 plant 2", "C. linearis Stage B3 plant 3",
                                              "C. rattanii all", "C. rattanii Stage B1 plant 1",
                                              "C. rattanii Stage B1 plant 2", "C. rattanii Stage B1 plant 3",
                                              "C. rattanii Stage B2 plant 1", "C. rattanii Stage B2 plant 2",
                                              "C. rattanii Stage B2 plant 3", "C. rattanii Stage B3 plant 1",
                                              "C. rattanii Stage B3 plant 2", "C. rattanii Stage B3 plant 3")

#remove bulk-processed samples (i.e., not analyzed for this project)
CollinsiaRNAreadcountsonlyDES <- CollinsiaRNAreadcountsonlyDES[,-c(1,11)]

#create a dataframe that has the experimental design data
species <- as.factor(c(rep("linearis", 9), rep("rattanii", 9)))
stage <- as.factor(c("B1","B1","B1","B2","B2","B2","B3","B3","B3",
                     "B1","B1","B1","B2","B2","B2","B3","B3","B3"))
CollinsiaRNAexpdesignDES <- data.frame(species, stage)
row.names(CollinsiaRNAexpdesignDES) <- colnames(CollinsiaRNAreadcountsonlyDES)

#create a "DESeqDataSet" to use in DESeq2 modeling
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = CollinsiaRNAreadcountsonlyDES,
  colData = CollinsiaRNAexpdesignDES,
  design = ~ species + stage + species:stage)

#conduct and save rlog-transformed data to use for visualizations
ddsRlogged <- rlog(ddsFullCountTable)

#### Genome-wide principle components analysis (PCA) & visualization ####

#ref: http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/
#load packages
install.packages("ggplot2") #if not done before
library(ggplot2)
install.packages("ggfortify") #if not done before
library(ggfortify)
install.packages("ggpubr") #if not done before
library(ggpubr)
install.packages("factoextra") #if not done before
library(factoextra)
install.packages("FactoMineR") #if not done before
library(FactoMineR)
install.packages("dplyr") #if not done before
library(dplyr)

#check out raw read count data first to test variability differences between species
#remove columns with only 0's
rawcounts <- CollinsiaRNAreadcountsonlyDES[rowSums(CollinsiaRNAreadcountsonlyDES[]) > 0,] 

#calculate standard deviations, means, and coefficients of variation
testraw <- rawcounts %>%
  rowwise() %>%
  mutate(CratSD = sd(c(`C. rattanii Stage B1 plant 1`, `C. rattanii Stage B1 plant 2`, `C. rattanii Stage B1 plant 3`,
                       `C. rattanii Stage B2 plant 1`, `C. rattanii Stage B2 plant 2`, `C. rattanii Stage B2 plant 3`,
                       `C. rattanii Stage B3 plant 1`, `C. rattanii Stage B3 plant 2`, `C. rattanii Stage B3 plant 3`), na.rm = TRUE),
         ClinSD = sd(c(`C. linearis Stage B1 plant 1`, `C. linearis Stage B1 plant 2`, `C. linearis Stage B1 plant 3`,
                       `C. linearis Stage B2 plant 1`, `C. linearis Stage B2 plant 2`, `C. linearis Stage B2 plant 3`,
                       `C. linearis Stage B3 plant 1`, `C. linearis Stage B3 plant 2`, `C. linearis Stage B3 plant 3`), na.rm = TRUE),
         CratMean = mean(c(`C. rattanii Stage B1 plant 1`, `C. rattanii Stage B1 plant 2`, `C. rattanii Stage B1 plant 3`,
                           `C. rattanii Stage B2 plant 1`, `C. rattanii Stage B2 plant 2`, `C. rattanii Stage B2 plant 3`,
                           `C. rattanii Stage B3 plant 1`, `C. rattanii Stage B3 plant 2`, `C. rattanii Stage B3 plant 3`), na.rm = TRUE),
         ClinMean = mean(c(`C. linearis Stage B1 plant 1`, `C. linearis Stage B1 plant 2`, `C. linearis Stage B1 plant 3`,
                           `C. linearis Stage B2 plant 1`, `C. linearis Stage B2 plant 2`, `C. linearis Stage B2 plant 3`,
                           `C. linearis Stage B3 plant 1`, `C. linearis Stage B3 plant 2`, `C. linearis Stage B3 plant 3`), na.rm = TRUE),
         CratCV = CratSD/CratMean,
         ClinCV = ClinSD/ClinMean)

#find C. linearis mean coefficient of variation
summary(testraw$ClinCV, na.rm = TRUE) #1.0196
sd(testraw$ClinCV, na.rm = TRUE)

#find C. rattanii mean coefficient of variation
summary(testraw$CratCV, na.rm = TRUE) #0.830
sd(testraw$CratCV, na.rm = TRUE)

#conduct t-test on coefficients of variation
t.test(testraw$ClinCV, testraw$CratCV)

#get data in shape for PCA
data <- assay(ddsRlogged)
data <- data[1:46127, 1:18]
datat <- t(data)
datatdf <- as.data.frame(datat)

#remove columns with only 0's, leaving 38,589 expressed genes' columns
datatdf <- datatdf[, colSums(datatdf) != 0]

#compute PCA with prcomp
res.pca <- prcomp(datatdf, scale = FALSE)

#visualize eigenvalues (i.e, scree plot)
fviz_eig(res.pca) 

#show exact percentages of variances explained by each PC
eig.val <- get_eigenvalue(res.pca)
eig.val #variance.percent = ~50, ~14, ~7

#create plot of PCs 1&2 then 1&3 and combine
PCA12figure <- ggplot(res.pca, aes(x = PC1, y = PC2, color = factor(stage), shape = factor(species))) +
  geom_jitter(size = 4,width = 7, height = 7) +
  scale_shape_manual(values=c(17,19)) + 
  scale_color_manual(values=c("gold", "orange", "tomato")) + 
  xlab("PC1: 50% variance") +
  ylab("PC2: 14% variance") +
  labs(color = "stage", shape = "species") +
  theme_classic() +
  theme(legend.position = c(0.8, 0.2), legend.box = "horizontal",
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black"))

PCA13figure <- ggplot(res.pca, aes(x = PC1, y = PC3, color = factor(stage), shape = factor(species))) + 
  geom_point(size = 4) +
  scale_shape_manual(values=c(17,19)) + 
  scale_color_manual(values=c("gold", "orange", "tomato")) + 
  xlab("PC1: 50% variance") +
  ylab("PC3: 7% variance") +
  labs(color = "stage", shape = "species") +
  theme_classic()

PCAdualfigure <- ggarrange(PCA12figure, PCA13figure, labels = c("a", "b"), ncol = 2, nrow = 1,
                           common.legend = TRUE, legend = "top")
PCAdualfigure

#re-running PCA to extract more results info with FactoMineR
#refs: http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/ &
#http://www.sthda.com/english/wiki/get-pca-extract-the-results-for-individuals-variables-in-principal-component-analysis-r-software-and-data-mining
res.pca <- PCA(datatdf, graph = FALSE, scale = FALSE)
res.pca$eig #same results as prcomp

#find genes with expression levels that are correlated with PC coordinates
res.desc <- dimdesc(res.pca, axes = c(1:3), proba = 0.01)
res.desc$Dim.1 #genes with expression levels that are correlated with PC1 coordinates

#find genes with expression levels that are correlated with PC2 coordinates
corsDim2 <- res.desc$Dim.2
corsDim2 <- corsDim2[["quanti"]]
corsDim2 <- corsDim2[(corsDim2[,1] > 0.7 | corsDim2[,1] < -0.7),]
PC2correlatedgenes <- row.names(corsDim2)
PC2correlatedgenes #1,420 genes with > |0.7| correlation

#is there a significant correlation between sample bud stages and PC2?
ind <- get_pca_ind(res.pca)
coords <- ind$coord #coordinates by sample
cor.test(coords[,2], c(1,1,1,2,2,2,3,3,3,1,1,1,2,2,2,3,3,3),
         method = "pearson", conf.level = 0.95) #yes, marginally for both species' data, p-value = 0.01498 & r = 0.563

#what about between linearis sample bud stages and PC2?
cor.test(coords[1:9,2], c(1,1,1,2,2,2,3,3,3),
         method = "pearson", conf.level = 0.95) #yes, for linearis alone: sig p-value = 0.00748 & r=0.815

#graph it
x <- c(1,1,1,2,2,2,3,3,3)
y <- coords[1:9,2]
plot(x, y, xaxt="n", xlab = "bud stage", ylab = "PC2 coordinates", main = "C. linearis")
axis(side=1,at=x)
abline(lm(y~x), col="red")

#what about between rattanii sample bud stages and PC2?
cor.test(coords[10:18,2], c(1,1,1,2,2,2,3,3,3),
         method = "pearson", conf.level = 0.95) #no for rattanii alone: p-value = 0.8245 & r=.087

#graph it
x <- c(1,1,1,2,2,2,3,3,3)
y <- coords[10:18,2]
plot(x, y, xaxt="n", xlab = "bud stage", ylab = "PC2 coordinates", main = "C. rattanii")
axis(side=1,at=x)
abline(lm(y~x), col="blue")

#find genes with expression levels that are correlated with PC3 coordinates
corsDim3 <- res.desc$Dim.3
corsDim3 <- corsDim3[["quanti"]]
corsDim3 <- corsDim3[(corsDim3[,1] > 0.7 | corsDim3[,1] < -0.7),]
PC3correlatedgenes <- row.names(corsDim3)
PC3correlatedgenes #529 genes with > |0.7| correlation

#is there a significant correlation between sample bud stages and PC3?
cor.test(coords[,3], c(1,1,1,2,2,2,3,3,3,1,1,1,2,2,2,3,3,3),
         method = "pearson", conf.level = 0.95) #yes, negatively for both species' data, sig p-value = 0.009556 & r = -0.593

#what about between linearis sample bud stages and PC3?
cor.test(coords[1:9,3], c(1,1,1,2,2,2,3,3,3),
         method = "pearson", conf.level = 0.95) #no for linearis alone: p-value = 0.1362 & r=-0.5367686

#graph it
x <- c(1,1,1,2,2,2,3,3,3)
y <- coords[1:9,3]
plot(x, y, xaxt="n", xlab = "bud stage", ylab = "PC3 coordinates", main = "C. linearis")
axis(side=1,at=x)
abline(lm(y~x), col="blue")

#what about between rattanii sample bud stages and PC3?
cor.test(coords[10:18,3], c(1,1,1,2,2,2,3,3,3),
         method = "pearson", conf.level = 0.95) #yes, for rattanii alone: sig p-value = 0.005571 & r=-0.8304857

#graph it
x <- c(1,1,1,2,2,2,3,3,3)
y <- coords[10:18,3]
plot(x, y, xaxt="n", xlab = "bud stage", ylab = "PC3 coordinates", main = "C. rattanii")
axis(side=1,at=x)
abline(lm(y~x), col="red")


#### DESeq2 full model & likelihood ratio test (LRT) for interaction term ####

dds <- DESeq(ddsFullCountTable) #actually running DESeq2 modeling

#"In this case, using the likelihood ratio test with a reduced model which does not contain
#the interaction term will test whether the condition induces a change in gene expression
#at any time point after the base-level time point." --Love et al. 2014
ddsTimeSeries <- DESeq(dds, test="LRT", reduced = ~ species + stage)
resultsTimeSeries <- results(ddsTimeSeries, alpha = 0.01)

#get stats about independent filtering, outlier removal, etc.
summary(resultsTimeSeries) #38,589 genes with non-zero total read counts

#subset out 855 significant DEGs, i.e."interaction effect genes" here, & create gene list
resultsTimeSeriesSig <- resultsTimeSeries[ which(resultsTimeSeries$padj < 0.01 ), ]
resultsTimeSeriesSiggenes <- row.names(resultsTimeSeriesSig)


#### Contrasts to find early developmental (ED) differential expression ####

#"This term is a test for if the RAT vs LIN fold change is different at [STAGE B2] than
#at [STAGE B1]." (adapted from M.Love via Bioconductor)
resultsTimeSeriesX12 <- results(dds, name="speciesrattanii.stageB2", alpha = 0.01, test="Wald")

#get stats about results, independent filtering, outlier removal, etc.
summary(resultsTimeSeriesX12)

#subset out 164 significant DEGs, i.e."early dev. contrast effect genes" here, & create gene list
resultsTimeSeriesX12Sig <- resultsTimeSeriesX12[ which(resultsTimeSeriesX12$padj < 0.01 ), ] 
resultsTimeSeriesX12Siggenes <- row.names(resultsTimeSeriesX12Sig)

#get stats on average significant differential expression in ED
summary(resultsTimeSeriesX12Sig$log2FoldChange)
sd(resultsTimeSeriesX12Sig$log2FoldChange)

#116 genes have an positive interaction term (i.e., rattanii's lfc is + relative to linearis's)
genesupED <- resultsTimeSeriesX12Sig[resultsTimeSeriesX12Sig$log2FoldChange > 0,] 
summary(genesupED$log2FoldChange) #7.85 mean lfc
sd(genesupED$log2FoldChange) #w/4.35 sd

#48 genes have a negative interaction term (i.e., rattanii's lfc is - relative to linearis's)
genesdownED <- resultsTimeSeriesX12Sig[resultsTimeSeriesX12Sig$log2FoldChange < 0,]
summary(genesdownED$log2FoldChange) #-9.44 mean lfc
sd(genesdownED$log2FoldChange) #w/5.98 sd


#### Contrasts to find late developmental (LD) differential expression ####

#this tests if the species effect is different in III compared to II
resultsTimeSeriesX23 <- results(dds, contrast=list("speciesrattanii.stageB3", "speciesrattanii.stageB2"),
       alpha = 0.01, test="Wald")

#get stats about results, independent filtering, outlier removal, etc.
summary(resultsTimeSeriesX23)

#subset out 418 significant DEGs, i.e."late dev. contrast effect genes" here, & create gene list
resultsTimeSeriesX23Sig <- resultsTimeSeriesX23[ which(resultsTimeSeriesX23$padj < 0.01 ), ]
resultsTimeSeriesX23Siggenes <- row.names(resultsTimeSeriesX23Sig)

#get stats on average significant differential expression in LD
summary(resultsTimeSeriesX23Sig$log2FoldChange)
sd(resultsTimeSeriesX23Sig$log2FoldChange)

#165 genes have an positive interaction term (i.e., rattanii's lfc is + relative to linearis's)
genesupLD <- resultsTimeSeriesX23Sig[resultsTimeSeriesX23Sig$log2FoldChange > 0,]
summary(genesupLD$log2FoldChange) #6.37 mean lfc
sd(genesupLD$log2FoldChange) #w/7.28 sd

#253 genes have a negative interaction term (i.e., rattanii's lfc is - relative to linearis's)
genesdownLD <- resultsTimeSeriesX23Sig[resultsTimeSeriesX23Sig$log2FoldChange < 0,]
summary(genesdownLD$log2FoldChange) #-3.34 mean lfc
sd(genesdownLD$log2FoldChange) #w/2.77 sd

##create list of 569 total contrast effect genes (would be 418 + 164 = 582,
#but 13 genes are in both ED & LD results)
resultsNPDEpairwisegenes <- union(resultsTimeSeriesX12Siggenes, resultsTimeSeriesX23Siggenes)
intersect(resultsTimeSeriesX12Siggenes, resultsTimeSeriesX23Siggenes)


#### DESeq2 modeling & LRT to find parallel gene expression ####

##first find all genes with a significant ‘stage’ effect

#set full model to account for species & stage, reduced = ~species & run LRT
design(dds) <- ~ species + stage
ddsLRTstage <- DESeq(dds, test="LRT", reduced = ~ species)
resultsLRTstage <- results(ddsLRTstage, alpha = 0.01)

#get stats about independent filtering, outlier removal, etc.
summary(resultsLRTstage)

#subset out 1101 significant "stage effect" genes here, & create gene list
resultsLRTstageSig <- resultsLRTstage[ which(resultsLRTstage$padj < 0.01 ), ]
resultsLRTstageSiggenes <- row.names(resultsLRTstageSig)

##then exclude the subset of the ‘stage effect’ genes that are also "interaction effect"
#genes, based on the full model used to analyze the 'species x stage' interaction term, above

#this leaves 1022 'shared effect' genes
resultsPDELRTgenes <- base::setdiff(resultsLRTstageSiggenes, resultsTimeSeriesSiggenes)









#### Prep topGO package and GO data ####

BiocManager::install("topGO") #if not done before
library(topGO)

#ref: https://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html
#read in GO annotations for all genes in genome
geneID2GO <- readMappings(file = "bp_cora.tab")

#create list of all genes in genome
geneUniverse <- names(geneID2GO)

#### Find (855) interaction effect genes' enriched GO terms ####

genesOfInterest <- resultsTimeSeriesSiggenes

#"tell TopGO where these interesting genes appear in the 'geneUniverse' vector" --avrilomics 2015
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

#create object of type 'topGOdata'
myGOdata <- new("topGOdata", description="interaction effect genes", ontology="BP",
                allGenes=geneList, nodeSize = 3, annot = annFUN.gene2GO, gene2GO = geneID2GO)
myGOdata
sigGenes(myGOdata) #368 out of 855 genes have GO annotations to assess for enrichment

#run enrichment test using default 'weight01' topGO algorithm and Fisher statistic
resultinteractionGO <- runTest(myGOdata, statistic="fisher")
resultinteractionGO #11 terms with p < 0.01

#create output in table format for 11 significant terms
ResinteractionGO <- GenTable(myGOdata, fisher = resultinteractionGO,
                             orderBy = resultinteractionGO,
                             ranksOf = resultinteractionGO, topNodes =11)
ResinteractionGO

#list implicated/"observed" genes that are annotated (directly or indirectly
#through the GO hierarchy) to these significant GO terms (as adapted from avrilomics)
myterms = ResinteractionGO$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
for (i in 1:length(myterms)) {
  myterm <- myterms[i]
  mygenesforterm <- mygenes[myterm][[1]]
  mygenesforterm <- mygenesforterm[mygenesforterm %in% resultsTimeSeriesSiggenes]
  mygenesforterms <- paste(mygenesforterm, collapse=',') 
  print(paste("Term",myterm,"genes:",mygenesforterms))
}

#and/or list only genes that are directly annotated to these significant GO terms
myterms = ResinteractionGO$GO.ID
GOtermsbp <- read.delim("bp_cora.tab", header = FALSE)
GOtermsbp <- GOtermsbp[GOtermsbp$V1 %in% resultsTimeSeriesSiggenes,]
GOtermsbp <- droplevels(GOtermsbp)
for (i in 1:length(myterms)) {
  directlyannotatedgenes <- as.character(GOtermsbp[grep((myterms[i]), GOtermsbp$V2), 1])
  print(noquote(paste((myterms[i]), "is directly annotated to", sep = " ")))
  print(noquote(directlyannotatedgenes))
}

#### Find (1022) shared effect genes' enriched GO terms ####

genesOfInterest <- resultsPDELRTgenes

#"tell TopGO where these interesting genes appear in the 'geneUniverse' vector" --avrilomics 2015
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

#create object of type 'topGOdata'
myGOdata <- new("topGOdata", description="shared effect genes", ontology="BP",
                allGenes=geneList, nodeSize = 3, annot = annFUN.gene2GO, gene2GO = geneID2GO)
myGOdata
sigGenes(myGOdata) #415 out of 1022 genes have GO annotations to assess for enrichment

#run enrichment test using default 'weight01' topGO algorithm and Fisher statistic
resultsharedGO <- runTest(myGOdata, statistic="fisher")
resultsharedGO #15 terms with p < 0.01

#create output in table format for 15 significant terms
RessharedGO <- GenTable(myGOdata, fisher = resultsharedGO, orderBy = resultsharedGO,
                      ranksOf = resultsharedGO, topNodes =15)
RessharedGO

#list implicated/"observed" genes that are annotated (directly or indirectly
#through the GO hierarchy) to these significant GO terms (as adapted from avrilomics)
myterms = RessharedGO$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
for (i in 1:length(myterms)) {
  myterm <- myterms[i]
  mygenesforterm <- mygenes[myterm][[1]]
  mygenesforterm <- mygenesforterm[mygenesforterm %in% resultsPDELRTgenes]
  mygenesforterms <- paste(mygenesforterm, collapse=',') 
  print(paste("Term",myterm,"genes:",mygenesforterms))
}

#and/or list only genes that are directly annotated to these significant GO terms
myterms = RessharedGO$GO.ID
GOtermsbp <- read.delim("bp_cora.tab", header = FALSE)
GOtermsbp <- GOtermsbp[GOtermsbp$V1 %in% resultsPDELRTgenes,]
GOtermsbp <- droplevels(GOtermsbp)
for (i in 1:length(myterms)) {
  directlyannotatedgenes <- as.character(GOtermsbp[grep((myterms[i]), GOtermsbp$V2), 1])
  print(noquote(paste((myterms[i]), "is directly annotated to", sep = " ")))
  print(noquote(directlyannotatedgenes))
}

#### Find (569) ED + LD contrast effect genes' enriched GO terms ####

genesOfInterest <- resultsNPDEpairwisegenes

#"tell TopGO where these interesting genes appear in the 'geneUniverse' vector" --avrilomics 2015
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

#create object of type 'topGOdata'
myGOdata <- new("topGOdata", description="resultsNPDEpairwisegenes", ontology="BP",
                allGenes=geneList, nodeSize = 3, annot = annFUN.gene2GO, gene2GO = geneID2GO)
myGOdata
sigGenes(myGOdata) #231 out of 569 genes have GO annotations to assess for enrichment

#run enrichment test using default 'weight01' topGO algorithm and Fisher statistic
resultcontrastsGO <- runTest(myGOdata, statistic="fisher")
resultcontrastsGO #6 terms with p < 0.01

#create output in table format for 6 significant terms
RescontrastsGO <- GenTable(myGOdata, fisher = resultcontrastsGO, orderBy = resultcontrastsGO,
                            ranksOf = resultcontrastsGO, topNodes = 6)
RescontrastsGO

#list implicated/"observed" genes that are annotated (directly or indirectly
#through the GO hierarchy) to these significant GO terms (as adapted from avrilomics)
myterms = RescontrastsGO$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
for (i in 1:length(myterms)) {
  myterm <- myterms[i]
  mygenesforterm <- mygenes[myterm][[1]]
  mygenesforterm <- mygenesforterm[mygenesforterm %in% resultsNPDEpairwisegenes]
  mygenesforterms <- paste(mygenesforterm, collapse=',') 
  print(paste("Term",myterm,"genes:",mygenesforterms))
}

#and/or list only genes that are directly annotated to these significant GO terms
myterms = RescontrastsGO$GO.ID
GOtermsbp <- read.delim("bp_cora.tab", header = FALSE)
GOtermsbp <- GOtermsbp[GOtermsbp$V1 %in% resultsNPDEpairwisegenes,]
GOtermsbp <- droplevels(GOtermsbp)
for (i in 1:length(myterms)) {
  directlyannotatedgenes <- as.character(GOtermsbp[grep((myterms[i]), GOtermsbp$V2), 1])
  print(noquote(paste((myterms[i]), "is directly annotated to", sep = " ")))
  print(noquote(directlyannotatedgenes))
}

#### Find (164) ED contrast effect genes' enriched GO terms ####

genesOfInterest <- resultsTimeSeriesX12Siggenes

#"tell TopGO where these interesting genes appear in the 'geneUniverse' vector" --avrilomics 2015
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

#create object of type 'topGOdata'
myGOdata <- new("topGOdata", description="resultsTimeSeriesX12Siggenes", ontology="BP",
                allGenes=geneList, nodeSize = 3, annot = annFUN.gene2GO, gene2GO = geneID2GO)
myGOdata
sigGenes(myGOdata) #59 out of 164 genes have GO annotations to assess for enrichment

#run enrichment test using default 'weight01' topGO algorithm and Fisher statistic
resultEDcontrastGO <- runTest(myGOdata, statistic="fisher")
resultEDcontrastGO #4 terms with p < 0.01

#create output in table format for 4 significant terms
ResEDcontrastGO <- GenTable(myGOdata, fisher = resultEDcontrastGO,
                             orderBy = resultEDcontrastGO,
                             ranksOf = resultEDcontrastGO, topNodes = 4)
ResEDcontrastGO

#### Find (418) LD contrast effect genes' enriched GO terms ####

genesOfInterest <- resultsTimeSeriesX23Siggenes

#"tell TopGO where these interesting genes appear in the 'geneUniverse' vector" --avrilomics 2015
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

#create object of type 'topGOdata'
myGOdata <- new("topGOdata", description="resultsTimeSeriesX23Siggenes", ontology="BP",
                allGenes=geneList, nodeSize = 3, annot = annFUN.gene2GO, gene2GO = geneID2GO)
myGOdata
sigGenes(myGOdata) #175 out of 418 genes have GO annotations to assess for enrichment

#run enrichment test using default 'weight01' topGO algorithm and Fisher statistic
resultLDcontrastGO <- runTest(myGOdata, statistic="fisher")
resultLDcontrastGO #2 terms with p < 0.01

#create output in table format for 2 significant terms
ResLDcontrastGO <- GenTable(myGOdata, fisher = resultLDcontrastGO,
                             orderBy = resultLDcontrastGO,
                             ranksOf = resultLDcontrastGO, topNodes = 2)
ResLDcontrastGO


#### Prep karyoploteR package & annotation files ####

BiocManager::install("karyoploteR") #if not done before
library(karyoploteR)
install("dplyr") #if not done before
library(dplyr)

#set C. rattanii major scaffold genome characteristics
custom.genome <- toGRanges(data.frame(chr=c("1","2","3","4","5","6","7"),
                                      start=c(1, 1, 1, 1, 1, 1, 1),
                                      end=c(79579559,69253743,69988318,56907126,56249891,95162470,60324063)))

#read in and format locations of all C. rattanii genes across genome
allfeats <- read.delim("CORA.gff", sep = "\t", header = FALSE)

#### Make GRanges object out of all genes' data ####
#organizing this data for GRangres object creation
colnames(allfeats)[7] <- "byStrand"
allfeats <- allfeats[,c(1,4,5,7,9,3,2)]
allfeats <- allfeats[allfeats$V2 == "ERGO",] #getting rid of "#" commented rows
allgenes <- allfeats[allfeats$V3 == "gene",] #47210 genes across whole genome are in here now!
allgenes <- droplevels(allgenes)

#reformatting column with ID so that I know which row is which gene
install.packages("stringr") #if not done before
library(stringr)
allgenes$V9 <- str_split_fixed(allgenes$V9, "=", 2)
allgenes$V9 <- allgenes$V9[,2]
allgenes$V9 <- as.factor(allgenes$V9)

#subset out major chromosome genes only & rename conveniently
majorscafallgenes <- allgenes[allgenes$V1 %in% c("PGA_scaffold_1__167_contigs__length_79579559",
                                                 "PGA_scaffold_2__194_contigs__length_69253743",
                                                 "PGA_scaffold_3__178_contigs__length_69988318",
                                                 "PGA_scaffold_4__128_contigs__length_56907126",
                                                 "PGA_scaffold_5__124_contigs__length_56249891",
                                                 "PGA_scaffold_6__239_contigs__length_95162470",
                                                 "PGA_scaffold_7__160_contigs__length_60324063"),]
majorscafallgenes$V1 <- dplyr::recode(majorscafallgenes$V1, PGA_scaffold_1__167_contigs__length_79579559 = "1",
                               PGA_scaffold_2__194_contigs__length_69253743 = "2",
                               PGA_scaffold_3__178_contigs__length_69988318 = "3",
                               PGA_scaffold_4__128_contigs__length_56907126 = "4",
                               PGA_scaffold_5__124_contigs__length_56249891 = "5",
                               PGA_scaffold_6__239_contigs__length_95162470 = "6",
                               PGA_scaffold_7__160_contigs__length_60324063 = "7")

#create GRanges object
majorscafallgenes <- droplevels(majorscafallgenes)
majorscafallgenes <- majorscafallgenes[,c("V1", "V4", "V5", "byStrand", "V9", "V3", "V2")]
majorscafallgenesGRanges <- toGRanges(majorscafallgenes) 
head(majorscafallgenesGRanges)


#### Physically map all gene density along 7 chromosomes ####

kp <- plotKaryotype(genome=custom.genome, plot.type=4, ideogram.plotter = NULL, labels.plotter = NULL)

#place representative chromosome lines & their numbers
kpAddCytobandsAsLine(kp)
kpAddChromosomeNames(kp, srt=0)

#plot density and add y-axis counts
kp <- kpPlotDensity(kp, majorscafallgenesGRanges)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, cex=0.8)

#### Make GRanges object out of all genes' data with GO annotations ####

#load annotations
CORAGOBP <- read.delim("bp_cora.tab", sep = "\t", header = FALSE)

#rename gene name & annotation columns to merge data based on gene names
allgenes <- allgenes %>% dplyr::rename(gene = V9)
CORAGOBP <- CORAGOBP %>% dplyr::rename(gene = V1, GOterms = V2)
allgenesannot <- merge(allgenes, CORAGOBP, by = "gene", all = TRUE)

#remove genes (rows) with 0 GO terms
allgenesannot <- allgenesannot[!(allgenesannot$GOterms == ""),]

#subset to major scaffolds & rename conveniently
majorscafallgenesannot <- allgenesannot[allgenesannot$V1 %in% c("PGA_scaffold_1__167_contigs__length_79579559",
                                                                "PGA_scaffold_2__194_contigs__length_69253743",
                                                                "PGA_scaffold_3__178_contigs__length_69988318",
                                                                "PGA_scaffold_4__128_contigs__length_56907126",
                                                                "PGA_scaffold_5__124_contigs__length_56249891",
                                                                "PGA_scaffold_6__239_contigs__length_95162470",
                                                                "PGA_scaffold_7__160_contigs__length_60324063"),]
majorscafallgenesannot$V1 <- dplyr::recode(majorscafallgenesannot$V1, PGA_scaffold_1__167_contigs__length_79579559 = "1",
                                    PGA_scaffold_2__194_contigs__length_69253743 = "2",
                                    PGA_scaffold_3__178_contigs__length_69988318 = "3",
                                    PGA_scaffold_4__128_contigs__length_56907126 = "4",
                                    PGA_scaffold_5__124_contigs__length_56249891 = "5",
                                    PGA_scaffold_6__239_contigs__length_95162470 = "6",
                                    PGA_scaffold_7__160_contigs__length_60324063 = "7")

#create GRanges object
majorscafallgenesannot <- droplevels(majorscafallgenesannot)
majorscafallgenesannot$GOterms <- as.character(majorscafallgenesannot$GOterms)
majorscafallgenesannot <- majorscafallgenesannot[,c("V1", "V4", "V5", "byStrand", "V3", "V2", "gene", "GOterms")]
majorscafallgenesannotGRanges <- toGRanges(majorscafallgenesannot)
head(majorscafallgenesannotGRanges)


#### Physically map only genes with GO annotations along 7 chromosomes ####

#place representative chromosome numbers & lines
kp <- plotKaryotype(genome = custom.genome, ideogram.plotter = NULL)
kpAddCytobandsAsLine(kp)

#plot genes as vertical lines
kpPlotRegions(kp, data=majorscafallgenesannotGRanges)

##plot pollination, reproduction, development, and growth-associated genes

#make a pollination GO term GRanges data set (=15 genes)
majorscafallPollGRanges <- toGRanges(majorscafallgenesannot[grep(("GO:0009856"),
                                                                 majorscafallgenesannot$GOterms),])

#make a reproduction GO term GRanges data set (=57 genes)
majorscafallReproGRanges <- toGRanges(majorscafallgenesannot[grep(("GO:0000003"),
                                                                  majorscafallgenesannot$GOterms),])

#make a developmental process GO term GRanges data set (=31 genes)
majorscafallDevelGRanges <- toGRanges(majorscafallgenesannot[grep(("GO:0032502"),
                                                                  majorscafallgenesannot$GOterms),])

#make a growth GO term GRanges data set (=20 genes)
majorscafallGrowGRanges <- toGRanges(majorscafallgenesannot[grep(("GO:0040007"),
                                                                 majorscafallgenesannot$GOterms),])

#create plot & lines
kp <- plotKaryotype(genome = custom.genome, ideogram.plotter = NULL)
kpAddCytobandsAsLine(kp)

#add all color-coded gene lines
kpPlotRegions(kp, data=majorscafallPollGRanges, avoid.overlapping = FALSE, col = "goldenrod3")
kpPlotRegions(kp, data=majorscafallReproGRanges, avoid.overlapping = FALSE, col = "orchid")
kpPlotRegions(kp, data=majorscafallDevelGRanges, avoid.overlapping = FALSE, col = "green3")
kpPlotRegions(kp, data=majorscafallGrowGRanges, avoid.overlapping = FALSE, col = "black")





#### Code for select heatmaps ####

#ref: https://www.bioconductor.org/help/course-materials/2015/CSAMA2015/lab/rnaseqCSAMA.html
#install and load packages
install.packages("pheatmap") #if not done before
library("pheatmap")
install.packages("RColorBrewer") #if not done before
library("RColorBrewer")


## HEATMAP OF 855 INTERACTION EFFECT GENES
mat <- assay(ddsRlogged)[ resultsTimeSeriesSiggenes, ]
mat <- t(scale(t(mat))) 
col.order <- c("C. linearis Stage B1 plant 1", "C. linearis Stage B1 plant 2", "C. linearis Stage B1 plant 3",
               "C. rattanii Stage B1 plant 1", "C. rattanii Stage B1 plant 2", "C. rattanii Stage B1 plant 3",
               "C. linearis Stage B2 plant 1", "C. linearis Stage B2 plant 2", "C. linearis Stage B2 plant 3",
               "C. rattanii Stage B2 plant 1", "C. rattanii Stage B2 plant 2", "C. rattanii Stage B2 plant 3",
               "C. linearis Stage B3 plant 1", "C. linearis Stage B3 plant 2", "C. linearis Stage B3 plant 3",
               "C. rattanii Stage B3 plant 1", "C. rattanii Stage B3 plant 2", "C. rattanii Stage B3 plant 3")
matstagespaired <- mat[,col.order]
df <- as.data.frame(colData(ddsRlogged)[,c("species","stage")])
my_colour = list(stage = c(B1 = "gold", B2 = "orange", B3 = "tomato"),
  species = c(linearis = "purple4", rattanii = "mediumpurple1"))
pheatmap(matstagespaired, cluster_cols = FALSE, annotation_col=df, annotation_color = my_colour,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100), show_rownames = FALSE)


## HEATMAP OF 569 CONTRAST EFFECT GENES
mat <- assay(ddsRlogged)[ resultsNPDEpairwisegenes, ]
col.order <- c("C. linearis Stage B1 plant 1", "C. linearis Stage B1 plant 2", "C. linearis Stage B1 plant 3",
               "C. rattanii Stage B1 plant 1", "C. rattanii Stage B1 plant 2", "C. rattanii Stage B1 plant 3",
               "C. linearis Stage B2 plant 1", "C. linearis Stage B2 plant 2", "C. linearis Stage B2 plant 3",
               "C. rattanii Stage B2 plant 1", "C. rattanii Stage B2 plant 2", "C. rattanii Stage B2 plant 3",
               "C. linearis Stage B3 plant 1", "C. linearis Stage B3 plant 2", "C. linearis Stage B3 plant 3",
               "C. rattanii Stage B3 plant 1", "C. rattanii Stage B3 plant 2", "C. rattanii Stage B3 plant 3")
mat <- mat[,col.order]
mat <- as.data.frame(mat)
mat2 <- mat %>% transmute(linearis_B1 = rowMeans(.[, c("C. linearis Stage B1 plant 1",
                                                       "C. linearis Stage B1 plant 2",
                                                       "C. linearis Stage B1 plant 3")]),
                          rattanii_B1 = rowMeans(.[, c("C. rattanii Stage B1 plant 1",
                                                       "C. rattanii Stage B1 plant 2",
                                                       "C. rattanii Stage B1 plant 3")]),
                          linearis_B2 = rowMeans(.[, c("C. linearis Stage B2 plant 1",
                                                       "C. linearis Stage B2 plant 2",
                                                       "C. linearis Stage B2 plant 3")]),
                          rattanii_B2 = rowMeans(.[, c("C. rattanii Stage B2 plant 1",
                                                       "C. rattanii Stage B2 plant 2",
                                                       "C. rattanii Stage B2 plant 3")]),
                          linearis_B3 = rowMeans(.[, c("C. linearis Stage B3 plant 1",
                                                       "C. linearis Stage B3 plant 2",
                                                       "C. linearis Stage B3 plant 3")]),
                          rattanii_B3 = rowMeans(.[, c("C. rattanii Stage B3 plant 1",
                                                       "C. rattanii Stage B3 plant 2",
                                                       "C. rattanii Stage B3 plant 3")]))

mat3 <- mat2 %>% transmute(linearis_B2B1 = linearis_B2 - linearis_B1,
                           rattanii_B2B1 = rattanii_B2 - rattanii_B1,
                           linearis_B3B2 = linearis_B3 - linearis_B2,
                           rattanii_B3B2 = rattanii_B3 - rattanii_B2)

colnames(mat3) <- c("linearis stage B1-B2", "rattanii stage B1-B2",
                    "linearis stage B2-B3", "rattanii stage B2-B3")

#workin on ED side
rownames(mat3) <- rownames(mat)
mat3ED <- as.matrix(subset(mat3, rownames(mat3) %in% resultsTimeSeriesX12Siggenes))
mat3ED <- mat3ED[,1:2]
dfED <- data.frame(species = c("linearis", "rattanii"))
rownames(dfED) <- c("linearis stage B1-B2", "rattanii stage B1-B2")

#working on LD side
rownames(mat3) <- rownames(mat)
mat3LD <- as.matrix(subset(mat3, rownames(mat3) %in% resultsTimeSeriesX23Siggenes))
mat3LD <- mat3LD[,3:4]
dfLD <- data.frame(species = c("linearis", "rattanii"))
rownames(dfLD) <- c("linearis stage B2-B3", "rattanii stage B2-B3")

#settings
my_colour = list(
  #stage = c(B1 = "gold", B2 = "orange", B3 = "tomato"),
  species = c(linearis = "purple4", rattanii = "mediumpurple1"))
paletteLength = 100
myBreaks <- c(seq(min(mat3), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(mat3)/paletteLength, max(mat3), length.out=floor(paletteLength/2)))

EDheatmap <- pheatmap(mat3ED, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(paletteLength),
         cluster_cols = FALSE, annotation_col=dfED, show_rownames = FALSE, show_colnames = FALSE,
         annotation_colors = my_colour, #annotation_legend = FALSE, legend = FALSE, 
         breaks = myBreaks, annotation_names_col = FALSE)

LDheatmap <- pheatmap(mat3LD, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(paletteLength),
         cluster_cols = FALSE, annotation_col=dfLD, show_rownames = FALSE, show_colnames = FALSE,
         breaks = myBreaks, annotation_colors = my_colour, annotation_names_col = FALSE)

#combination ED+LD finale
ggarrange(EDheatmap[[4]], LDheatmap[[4]], labels = c("EARLY", "LATE"), hjust = c(0,-0.25))





#### my sessionInfo() ####
sessionInfo()
#> sessionInfo()
#R version 3.6.2 (2019-12-12)
#Platform: x86_64-apple-darwin15.6.0 (64-bit)
#Running under: macOS Mojave 10.14.6

#Matrix products: default
#BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
#LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

#locale:
#  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#attached base packages:
#  [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
#[8] methods   base     

#other attached packages:
#  [1] RColorBrewer_1.1-2          pheatmap_1.0.12            
#[3] karyoploteR_1.12.4          regioneR_1.18.1            
#[5] topGO_2.38.1                SparseM_1.78               
#[7] GO.db_3.10.0                AnnotationDbi_1.48.0       
#[9] graph_1.64.0                FactoMineR_2.3             
#[11] ggfortify_0.4.10            factoextra_1.0.7.999       
#[13] DESeq2_1.26.0               SummarizedExperiment_1.16.1
#[15] DelayedArray_0.12.2         BiocParallel_1.20.1        
#[17] matrixStats_0.56.0          Biobase_2.46.0             
#[19] GenomicRanges_1.38.0        GenomeInfoDb_1.22.0        
#[21] IRanges_2.20.2              S4Vectors_0.24.3           
#[23] BiocGenerics_0.32.0         rstatix_0.6.0              
#[25] ggpubr_0.4.0.999            forcats_0.5.0              
#[27] stringr_1.4.0               purrr_0.3.4                
#[29] readr_1.3.1                 tidyr_1.1.1                
#[31] tibble_3.0.3                ggplot2_3.3.2              
#[33] tidyverse_1.3.0             car_3.0-8                  
#[35] carData_3.0-4               logspline_2.1.16           
#[37] fitdistrplus_1.1-1          survival_3.2-3             
#[39] MASS_7.3-51.6               nortest_1.0-4              
#[41] dplyr_1.0.1               

#loaded via a namespace (and not attached):
#  [1] readxl_1.3.1             backports_1.1.8          Hmisc_4.4-0             
#[4] BiocFileCache_1.10.2     lazyeval_0.2.2           splines_3.6.2           
#[7] digest_0.6.25            ensembldb_2.10.2         htmltools_0.5.0         
#[10] fansi_0.4.1              magrittr_1.5             checkmate_2.0.0         
#[13] memoise_1.1.0            BSgenome_1.54.0          cluster_2.1.0           
#[16] openxlsx_4.1.5           Biostrings_2.54.0        annotate_1.64.0         
#[19] modelr_0.1.8             askpass_1.1              prettyunits_1.1.1       
#[22] jpeg_0.1-8.1             colorspace_1.4-1         rappdirs_0.3.1          
#[25] blob_1.2.1               rvest_0.3.6              ggrepel_0.8.2           
#[28] haven_2.3.1              xfun_0.16                crayon_1.3.4            
#[31] RCurl_1.98-1.2           jsonlite_1.7.0           genefilter_1.68.0       
#[34] VariantAnnotation_1.32.0 glue_1.4.1               gtable_0.3.0            
#[37] zlibbioc_1.32.0          XVector_0.26.0           abind_1.4-5             
#[40] scales_1.1.1             bezier_1.1.2             DBI_1.1.0               
#[43] Rcpp_1.0.5               progress_1.2.2           xtable_1.8-4            
#[46] htmlTable_2.0.1          flashClust_1.01-2        foreign_0.8-75          
#[49] bit_4.0.3                Formula_1.2-3            htmlwidgets_1.5.1       
#[52] httr_1.4.2               acepack_1.4.1            ellipsis_0.3.1          
#[55] pkgconfig_2.0.3          XML_3.99-0.3             farver_2.0.3            
#[58] nnet_7.3-14              dbplyr_1.4.4             locfit_1.5-9.4          
#[61] utf8_1.1.4               tidyselect_1.1.0         labeling_0.3            
#[64] rlang_0.4.7              munsell_0.5.0            cellranger_1.1.0        
#[67] tools_3.6.2              cli_2.0.2                generics_0.0.2          
#[70] RSQLite_2.2.0            broom_0.7.0              knitr_1.29              
#[73] bit64_4.0.2              fs_1.5.0                 zip_2.0.4               
#[76] AnnotationFilter_1.10.0  leaps_3.1                xml2_1.3.2              
#[79] biomaRt_2.42.1           compiler_3.6.2           rstudioapi_0.11         
#[82] curl_4.3                 png_0.1-7                ggsignif_0.6.0          
#[85] reprex_0.3.0             geneplotter_1.64.0       stringi_1.4.6           
#[88] GenomicFeatures_1.38.2   lattice_0.20-41          ProtGenerics_1.18.0     
#[91] Matrix_1.2-18            vctrs_0.3.2              pillar_1.4.6            
#[94] lifecycle_0.2.0          BiocManager_1.30.10      data.table_1.12.8       
#[97] cowplot_1.0.0            bitops_1.0-6             rtracklayer_1.46.0      
#[100] R6_2.4.1                 latticeExtra_0.6-29      gridExtra_2.3           
#[103] rio_0.5.16               dichromat_2.0-0          assertthat_0.2.1        
#[106] openssl_1.4.2            withr_2.2.0              GenomicAlignments_1.22.1
#[109] Rsamtools_2.2.3          GenomeInfoDbData_1.2.2   hms_0.5.3               
#[112] grid_3.6.2               rpart_4.1-15             bamsignals_1.18.0       
#[115] biovizBase_1.34.1        scatterplot3d_0.3-41     lubridate_1.7.9         
#[118] base64enc_0.1-3