library(data.table)
library(Rsamtools)
library(tictoc)
library(stringr)
library(testthat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(class)
library(mlbench)
library(caret)

hg38=TRUE
cv <- clinvarCaller(newPileups=TRUE,newFalsePileups=TRUE,newLociAllFalse=TRUE,runsAtImid=FALSE,hg38=TRUE,returnMergedDb=FALSE,
                       returnStats=FALSE,returnPileups=TRUE,returnLociAll=TRUE,returnLociAllFalse=TRUE,returnHistogramData=TRUE,
                       verbose=TRUE,minBaseQuality=0,minMapq=0,fullPileupStats=TRUE,
                       varDbFormat="tsv",
                       chromosomeLengthsPath="./input/chromosomeLengths.csv",
                       bamPath="./input/corriell_S7-ready.bam",
                       bedPath="./input/corriell_S7-ready.per-base.bed.gz",
                       inputVariantsDbPath="./input/clinvar_alleles.single.b38.tsv.gz",
                       vcfHcGeneralPath="./input/corriell_S7-HCScan.vcf.gz",
                       vcfHcFcPath="./input/corriell_S7-HCScanFC.vcf.gz",
                       outputLociAllPath="./output/lociAll.positions",
                       outputLociAllFalsePath="./output/lociAllFalse.positions",
                       outputPileupsPath="./output/pileups.mpileup",
                       outputPileupsFalsePath="./output/pileupsFalse.mpileup",
                       outputStatsPath="./output/stats.csv",
                       outputVarDbMergedPath="./output/varDbMerged-reduced.csv")


#rozne progi


pileups2 <- copy(cv$pileups)
pileupsRandom2 <- copy(cv$falsePileups)

pileups2[,"VAR":="variant"]
pileupsRandom2[,"VAR":="non-variant"]
pileupsAll <- rbind(pileups2,pileupsRandom2)
pileupsAll <- addRefFromReferenceGenome(pileupsAll)
pileupsAll <- getVRF(pileupsAll,getDP=TRUE,getRD=TRUE)
#pileupsAll[,"bias":=runif(nrow(pileupsAll),0.0,0.00001)]

pileupsAllShuffled <- pileupsAll[sample(nrow(pileupsAll)),]

rm(pileups2,pileupsRandom2,pileupsAll)

#knn_result_5_3 <- knn(pileupsAllShuffled[(1+nrow(pileupsAllShuffled)/2):(nrow(pileupsAllShuffled)),c("DP","RD","VRF","bias")], pileupsAllShuffled[1:(nrow(pileupsAllShuffled)/2),c("DP","RD","VRF","bias")], pileupsAllShuffled[(1+nrow(pileupsAllShuffled)/2):(nrow(pileupsAllShuffled))]$VAR, k=5, l=3)
#print("k=5, l=3")
#print(paste0("TRUE: ",(round(sum((knn_result_5_3 == pileupsAllShuffled[1:(nrow(pileupsAllShuffled)/2)]$VAR)==TRUE,na.rm=TRUE)/(nrow(pileupsAllShuffled)/2),digits=4)),"%, (",(sum((knn_result_5_3 == pileupsAllShuffled[1:(nrow(pileupsAllShuffled)/2)]$VAR)==TRUE,na.rm=TRUE)),"/",(nrow(pileupsAllShuffled)/2),")."))
#print(paste0("FALSE: ",(round(sum((knn_result_5_3 == pileupsAllShuffled[1:(nrow(pileupsAllShuffled)/2)]$VAR)==FALSE,na.rm=TRUE)/(nrow(pileupsAllShuffled)/2),digits=4)),"%, (",(sum((knn_result_5_3 == pileupsAllShuffled[1:(nrow(pileupsAllShuffled)/2)]$VAR)==FALSE,na.rm=TRUE)),"/",(nrow(pileupsAllShuffled)/2),")."))
#confusionMatrix(knn_result_5_3,as.factor(unlist(pileupsAllShuffled[1:(nrow(pileupsAllShuffled)/2),"VAR"])))

control <- trainControl(method="repeatedcv", number=5, repeats=1)
seed <- 7
metric <- "Accuracy"
preProcessp =c("center", "scale")

set.seed(seed)
fit.lda <- train(VAR~., data=pileupsAllShuffled[,!c("CHROM","POS")], method="lda", metric=metric, preProc=c("center", "scale"), trControl=control)
set.seed(seed)
fit.knn <- train(VAR~., data=pileupsAllShuffled[,!c("CHROM","POS")], method="knn", metric=metric, preProc=c("center", "scale"), trControl=control)
results <- resamples(list(lda=fit.lda, knn=fit.knn))
# Table comparison
summary(results)

# boxplot comparison
bwplot(results)
# Dot-plot comparison
dotplot(results)

