library(data.table)
library(Rsamtools)
library(tictoc)
library(stringr)
library(testthat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(class)
library(caret)

#pileups2 <- copy(pileups)
#pileupsRandom2 <- copy(pileupsRandom)
cv00 <- clinvarCaller(newPileups=TRUE,newFalsePileups=TRUE,newNotVariantPositions=TRUE,runsAtImid=FALSE,hg38=TRUE,returnMergedDb=TRUE,
                    returnStats=FALSE,returnPileups=TRUE,returnTmpPileups=TRUE,returnNotVariantPositions=TRUE,returnHistogramData=TRUE,
                    verbose=TRUE,minBaseQuality=0,minMapq=0,varDbFormat="tsv",chromosomeLengthsPath="./input/chromosomeLengths.csv",
                    bamPath="./input/corriell_S7-ready.bam",bedPath="./input/corriell_S7-ready.per-base.bed.gz",
                    inputVariantsDbPath="./input/clinvar_alleles.single.b38.tsv.gz",vcfHcGeneralPath="./input/corriell_S7-HCScan.vcf.gz",
                    vcfHcFcPath="./input/corriell_S7-HCScanFC.vcf.gz",outputPileupsPath="./output/pileups.csv",
                    outputFalsePileupsPath="./output/falsePileups.csv",outputStatsPath="./output/stats.csv",
                    outputVarDbMergedPath="./output/varDbMerged-reduced.csv")

cv0bqmq <- clinvarCaller(newPileups=TRUE,newFalsePileups=TRUE,newNotVariantPositions=FALSE,runsAtImid=FALSE,hg38=TRUE,returnMergedDb=FALSE,
                    returnStats=FALSE,returnPileups=TRUE,returnTmpPileups=FALSE,returnNotVariantPositions=FALSE,returnHistogramData=TRUE,
                    verbose=TRUE,minBaseQuality=0,minMapq=0,varDbFormat="tsv",chromosomeLengthsPath="./input/chromosomeLengths.csv",
                    bamPath="./input/corriell_S7-ready.bam",bedPath="./input/corriell_S7-ready.per-base.bed.gz",
                    inputVariantsDbPath="./input/clinvar_alleles.single.b38.tsv.gz",vcfHcGeneralPath="./input/corriell_S7-HCScan.vcf.gz",
                    vcfHcFcPath="./input/corriell_S7-HCScanFC.vcf.gz",outputPileupsPath="./output/pileups_cv0bq.csv",
                    outputFalsePileupsPath="./output/falsePileups_cv0bq.csv",outputStatsPath="./output/stats.csv",
                    outputVarDbMergedPath="./output/varDbMerged-reduced.csv",notVariantPositionsTmp=cv00$notVariantPositions)

cv1bq <- clinvarCaller(newPileups=TRUE,newFalsePileups=TRUE,newNotVariantPositions=FALSE,runsAtImid=FALSE,hg38=TRUE,returnMergedDb=FALSE,
                    returnStats=FALSE,returnPileups=TRUE,returnTmpPileups=FALSE,returnNotVariantPositions=FALSE,returnHistogramData=TRUE,
                    verbose=TRUE,minBaseQuality=10,minMapq=0,varDbFormat="tsv",chromosomeLengthsPath="./input/chromosomeLengths.csv",
                    bamPath="./input/corriell_S7-ready.bam",bedPath="./input/corriell_S7-ready.per-base.bed.gz",
                    inputVariantsDbPath="./input/clinvar_alleles.single.b38.tsv.gz",vcfHcGeneralPath="./input/corriell_S7-HCScan.vcf.gz",
                    vcfHcFcPath="./input/corriell_S7-HCScanFC.vcf.gz",outputPileupsPath="./output/pileups_cv1bq.csv",
                    outputFalsePileupsPath="./output/falsePileups_cv1bq.csv",outputStatsPath="./output/stats.csv",
                    outputVarDbMergedPath="./output/varDbMerged-reduced.csv",notVariantPositionsTmp=cv00$notVariantPositions)

cv2bq <- clinvarCaller(newPileups=TRUE,newFalsePileups=TRUE,newNotVariantPositions=FALSE,runsAtImid=FALSE,hg38=TRUE,returnMergedDb=FALSE,
                    returnStats=FALSE,returnPileups=TRUE,returnTmpPileups=FALSE,returnNotVariantPositions=FALSE,returnHistogramData=TRUE,
                    verbose=TRUE,minBaseQuality=20,minMapq=0,varDbFormat="tsv",chromosomeLengthsPath="./input/chromosomeLengths.csv",
                    bamPath="./input/corriell_S7-ready.bam",bedPath="./input/corriell_S7-ready.per-base.bed.gz",
                    inputVariantsDbPath="./input/clinvar_alleles.single.b38.tsv.gz",vcfHcGeneralPath="./input/corriell_S7-HCScan.vcf.gz",
                    vcfHcFcPath="./input/corriell_S7-HCScanFC.vcf.gz",outputPileupsPath="./output/pileups_cv2bq.csv",
                    outputFalsePileupsPath="./output/falsePileups_cv2bq.csv",outputStatsPath="./output/stats.csv",
                    outputVarDbMergedPath="./output/varDbMerged-reduced.csv",notVariantPositionsTmp=cv00$notVariantPositions)

cv3bq <- clinvarCaller(newPileups=TRUE,newFalsePileups=TRUE,newNotVariantPositions=FALSE,runsAtImid=FALSE,hg38=TRUE,returnMergedDb=FALSE,
                    returnStats=FALSE,returnPileups=TRUE,returnTmpPileups=FALSE,returnNotVariantPositions=FALSE,returnHistogramData=TRUE,
                    verbose=TRUE,minBaseQuality=30,minMapq=0,varDbFormat="tsv",chromosomeLengthsPath="./input/chromosomeLengths.csv",
                    bamPath="./input/corriell_S7-ready.bam",bedPath="./input/corriell_S7-ready.per-base.bed.gz",
                    inputVariantsDbPath="./input/clinvar_alleles.single.b38.tsv.gz",vcfHcGeneralPath="./input/corriell_S7-HCScan.vcf.gz",
                    vcfHcFcPath="./input/corriell_S7-HCScanFC.vcf.gz",outputPileupsPath="./output/pileups_cv3bq.csv",
                    outputFalsePileupsPath="./output/falsePileups_cv3bq.csv",outputStatsPath="./output/stats.csv",
                    outputVarDbMergedPath="./output/varDbMerged-reduced.csv",notVariantPositionsTmp=cv00$notVariantPositions)

cv4bq <- clinvarCaller(newPileups=TRUE,newFalsePileups=TRUE,newNotVariantPositions=FALSE,runsAtImid=FALSE,hg38=TRUE,returnMergedDb=FALSE,
                    returnStats=FALSE,returnPileups=TRUE,returnTmpPileups=FALSE,returnNotVariantPositions=FALSE,returnHistogramData=TRUE,
                    verbose=TRUE,minBaseQuality=40,minMapq=0,varDbFormat="tsv",chromosomeLengthsPath="./input/chromosomeLengths.csv",
                    bamPath="./input/corriell_S7-ready.bam",bedPath="./input/corriell_S7-ready.per-base.bed.gz",
                    inputVariantsDbPath="./input/clinvar_alleles.single.b38.tsv.gz",vcfHcGeneralPath="./input/corriell_S7-HCScan.vcf.gz",
                    vcfHcFcPath="./input/corriell_S7-HCScanFC.vcf.gz",outputPileupsPath="./output/pileups_cv4bq.csv",
                    outputFalsePileupsPath="./output/falsePileups_cv4bq.csv",outputStatsPath="./output/stats.csv",
                    outputVarDbMergedPath="./output/varDbMerged-reduced.csv",notVariantPositionsTmp=cv00$notVariantPositions)

cv5bq <- clinvarCaller(newPileups=TRUE,newFalsePileups=TRUE,newNotVariantPositions=FALSE,runsAtImid=FALSE,hg38=TRUE,returnMergedDb=FALSE,
                    returnStats=FALSE,returnPileups=TRUE,returnTmpPileups=FALSE,returnNotVariantPositions=FALSE,returnHistogramData=TRUE,
                    verbose=TRUE,minBaseQuality=50,minMapq=0,varDbFormat="tsv",chromosomeLengthsPath="./input/chromosomeLengths.csv",
                    bamPath="./input/corriell_S7-ready.bam",bedPath="./input/corriell_S7-ready.per-base.bed.gz",
                    inputVariantsDbPath="./input/clinvar_alleles.single.b38.tsv.gz",vcfHcGeneralPath="./input/corriell_S7-HCScan.vcf.gz",
                    vcfHcFcPath="./input/corriell_S7-HCScanFC.vcf.gz",outputPileupsPath="./output/pileups_cv5bq.csv",
                    outputFalsePileupsPath="./output/falsePileups_cv5bq.csv",outputStatsPath="./output/stats.csv",
                    outputVarDbMergedPath="./output/varDbMerged-reduced.csv",notVariantPositionsTmp=cv00$notVariantPositions)

cv6bq <- clinvarCaller(newPileups=TRUE,newFalsePileups=TRUE,newNotVariantPositions=FALSE,runsAtImid=FALSE,hg38=TRUE,returnMergedDb=FALSE,
                    returnStats=FALSE,returnPileups=TRUE,returnTmpPileups=FALSE,returnNotVariantPositions=FALSE,returnHistogramData=TRUE,
                    verbose=TRUE,minBaseQuality=60,minMapq=0,varDbFormat="tsv",chromosomeLengthsPath="./input/chromosomeLengths.csv",
                    bamPath="./input/corriell_S7-ready.bam",bedPath="./input/corriell_S7-ready.per-base.bed.gz",
                    inputVariantsDbPath="./input/clinvar_alleles.single.b38.tsv.gz",vcfHcGeneralPath="./input/corriell_S7-HCScan.vcf.gz",
                    vcfHcFcPath="./input/corriell_S7-HCScanFC.vcf.gz",outputPileupsPath="./output/pileups_cv6bq.csv",
                    outputFalsePileupsPath="./output/falsePileups_cv6bq.csv",outputStatsPath="./output/stats.csv",
                    outputVarDbMergedPath="./output/varDbMerged-reduced.csv",notVariantPositionsTmp=cv00$notVariantPositions)


cv1mq <- clinvarCaller(newPileups=TRUE,newFalsePileups=TRUE,newNotVariantPositions=FALSE,runsAtImid=FALSE,hg38=TRUE,returnMergedDb=FALSE,
                       returnStats=FALSE,returnPileups=TRUE,returnTmpPileups=FALSE,returnNotVariantPositions=FALSE,returnHistogramData=TRUE,
                       verbose=TRUE,minBaseQuality=0,minMapq=10,varDbFormat="tsv",chromosomeLengthsPath="./input/chromosomeLengths.csv",
                       bamPath="./input/corriell_S7-ready.bam",bedPath="./input/corriell_S7-ready.per-base.bed.gz",
                       inputVariantsDbPath="./input/clinvar_alleles.single.b38.tsv.gz",vcfHcGeneralPath="./input/corriell_S7-HCScan.vcf.gz",
                       vcfHcFcPath="./input/corriell_S7-HCScanFC.vcf.gz",outputPileupsPath="./output/pileups_cv1mq.csv",
                       outputFalsePileupsPath="./output/falsePileups_cv1mq.csv",outputStatsPath="./output/stats.csv",
                       outputVarDbMergedPath="./output/varDbMerged-reduced.csv",notVariantPositionsTmp=cv00$notVariantPositions)

cv2mq <- clinvarCaller(newPileups=TRUE,newFalsePileups=TRUE,newNotVariantPositions=FALSE,runsAtImid=FALSE,hg38=TRUE,returnMergedDb=FALSE,
                       returnStats=FALSE,returnPileups=TRUE,returnTmpPileups=FALSE,returnNotVariantPositions=FALSE,returnHistogramData=TRUE,
                       verbose=TRUE,minBaseQuality=0,minMapq=20,varDbFormat="tsv",chromosomeLengthsPath="./input/chromosomeLengths.csv",
                       bamPath="./input/corriell_S7-ready.bam",bedPath="./input/corriell_S7-ready.per-base.bed.gz",
                       inputVariantsDbPath="./input/clinvar_alleles.single.b38.tsv.gz",vcfHcGeneralPath="./input/corriell_S7-HCScan.vcf.gz",
                       vcfHcFcPath="./input/corriell_S7-HCScanFC.vcf.gz",outputPileupsPath="./output/pileups_cv2mq.csv",
                       outputFalsePileupsPath="./output/falsePileups_cv2mq.csv",outputStatsPath="./output/stats.csv",
                       outputVarDbMergedPath="./output/varDbMerged-reduced.csv",notVariantPositionsTmp=cv00$notVariantPositions)

cv3mq <- clinvarCaller(newPileups=TRUE,newFalsePileups=TRUE,newNotVariantPositions=FALSE,runsAtImid=FALSE,hg38=TRUE,returnMergedDb=FALSE,
                       returnStats=FALSE,returnPileups=TRUE,returnTmpPileups=FALSE,returnNotVariantPositions=FALSE,returnHistogramData=TRUE,
                       verbose=TRUE,minBaseQuality=0,minMapq=30,varDbFormat="tsv",chromosomeLengthsPath="./input/chromosomeLengths.csv",
                       bamPath="./input/corriell_S7-ready.bam",bedPath="./input/corriell_S7-ready.per-base.bed.gz",
                       inputVariantsDbPath="./input/clinvar_alleles.single.b38.tsv.gz",vcfHcGeneralPath="./input/corriell_S7-HCScan.vcf.gz",
                       vcfHcFcPath="./input/corriell_S7-HCScanFC.vcf.gz",outputPileupsPath="./output/pileups_cv3mq.csv",
                       outputFalsePileupsPath="./output/falsePileups_cv3mq.csv",outputStatsPath="./output/stats.csv",
                       outputVarDbMergedPath="./output/varDbMerged-reduced.csv",notVariantPositionsTmp=cv00$notVariantPositions)

cv4mq <- clinvarCaller(newPileups=TRUE,newFalsePileups=TRUE,newNotVariantPositions=FALSE,runsAtImid=FALSE,hg38=TRUE,returnMergedDb=FALSE,
                       returnStats=FALSE,returnPileups=TRUE,returnTmpPileups=FALSE,returnNotVariantPositions=FALSE,returnHistogramData=TRUE,
                       verbose=TRUE,minBaseQuality=0,minMapq=40,varDbFormat="tsv",chromosomeLengthsPath="./input/chromosomeLengths.csv",
                       bamPath="./input/corriell_S7-ready.bam",bedPath="./input/corriell_S7-ready.per-base.bed.gz",
                       inputVariantsDbPath="./input/clinvar_alleles.single.b38.tsv.gz",vcfHcGeneralPath="./input/corriell_S7-HCScan.vcf.gz",
                       vcfHcFcPath="./input/corriell_S7-HCScanFC.vcf.gz",outputPileupsPath="./output/pileups_cv4mq.csv",
                       outputFalsePileupsPath="./output/falsePileups_cv4mq.csv",outputStatsPath="./output/stats.csv",
                       outputVarDbMergedPath="./output/varDbMerged-reduced.csv",notVariantPositionsTmp=cv00$notVariantPositions)

cv5mq <- clinvarCaller(newPileups=TRUE,newFalsePileups=TRUE,newNotVariantPositions=FALSE,runsAtImid=FALSE,hg38=TRUE,returnMergedDb=FALSE,
                       returnStats=FALSE,returnPileups=TRUE,returnTmpPileups=FALSE,returnNotVariantPositions=FALSE,returnHistogramData=TRUE,
                       verbose=TRUE,minBaseQuality=0,minMapq=50,varDbFormat="tsv",chromosomeLengthsPath="./input/chromosomeLengths.csv",
                       bamPath="./input/corriell_S7-ready.bam",bedPath="./input/corriell_S7-ready.per-base.bed.gz",
                       inputVariantsDbPath="./input/clinvar_alleles.single.b38.tsv.gz",vcfHcGeneralPath="./input/corriell_S7-HCScan.vcf.gz",
                       vcfHcFcPath="./input/corriell_S7-HCScanFC.vcf.gz",outputPileupsPath="./output/pileups_cv5mq.csv",
                       outputFalsePileupsPath="./output/falsePileups_cv5mq.csv",outputStatsPath="./output/stats.csv",
                       outputVarDbMergedPath="./output/varDbMerged-reduced.csv",notVariantPositionsTmp=cv00$notVariantPositions)

cv6mq <- clinvarCaller(newPileups=TRUE,newFalsePileups=TRUE,newNotVariantPositions=FALSE,runsAtImid=FALSE,hg38=TRUE,returnMergedDb=FALSE,
                       returnStats=FALSE,returnPileups=TRUE,returnTmpPileups=FALSE,returnNotVariantPositions=FALSE,returnHistogramData=TRUE,
                       verbose=TRUE,minBaseQuality=0,minMapq=60,varDbFormat="tsv",chromosomeLengthsPath="./input/chromosomeLengths.csv",
                       bamPath="./input/corriell_S7-ready.bam",bedPath="./input/corriell_S7-ready.per-base.bed.gz",
                       inputVariantsDbPath="./input/clinvar_alleles.single.b38.tsv.gz",vcfHcGeneralPath="./input/corriell_S7-HCScan.vcf.gz",
                       vcfHcFcPath="./input/corriell_S7-HCScanFC.vcf.gz",outputPileupsPath="./output/pileups_cv6mq.csv",
                       outputFalsePileupsPath="./output/falsePileups_cv6mq.csv",outputStatsPath="./output/stats.csv",
                       outputVarDbMergedPath="./output/varDbMerged-reduced.csv",notVariantPositionsTmp=cv00$notVariantPositions)



#rozne progi


pileups2 <- copy(cv00$pileupsTmp)
pileupsRandom2 <- copy(cv00$falsePileupsTmp)

pileups2 <- unifyColNamesInPileups(pileups2)
pileups2 <- addNucColumnsToPileups(pileups2)

pileupsRandom2 <- unifyColNamesInPileups(pileupsRandom2)
pileupsRandom2 <- addNucColumnsToPileups(pileupsRandom2)


pileups2[,"VAR":="variant"]
pileupsRandom2[,"VAR":="non-variant"]
pileupsAll <- rbind(pileups2,pileupsRandom2)
hg38=TRUE
pileupsAll <- addRefFromReferenceGenome(pileupsAll)



pileupsAllShuffled <- pileupsAll[sample(nrow(pileupsAll)),]
pileupsAllShuffled <- getVRF(pileupsAllShuffled,getDP=TRUE,getRD=TRUE)
pileupsAllShuffled[,"bias":=runif(nrow(pileupsAllShuffled),0.0,0.00001)]


knn_result_5_3 <- knn(pileupsAllShuffled[(1+nrow(pileupsAllShuffled)/2):(nrow(pileupsAllShuffled)),c("DP","RD","VRF","bias")], pileupsAllShuffled[1:(nrow(pileupsAllShuffled)/2),c("DP","RD","VRF","bias")], pileupsAllShuffled[(1+nrow(pileupsAllShuffled)/2):(nrow(pileupsAllShuffled))]$VAR, k=5, l=3)



print("k=5, l=3")
print(paste0("TRUE: ",(round(sum((knn_result_5_3 == pileupsAllShuffled[1:(nrow(pileupsAllShuffled)/2)]$VAR)==TRUE,na.rm=TRUE)/(nrow(pileupsAllShuffled)/2),digits=4)),"%, (",(sum((knn_result_5_3 == pileupsAllShuffled[1:(nrow(pileupsAllShuffled)/2)]$VAR)==TRUE,na.rm=TRUE)),"/",(nrow(pileupsAllShuffled)/2),")."))
print(paste0("FALSE: ",(round(sum((knn_result_5_3 == pileupsAllShuffled[1:(nrow(pileupsAllShuffled)/2)]$VAR)==FALSE,na.rm=TRUE)/(nrow(pileupsAllShuffled)/2),digits=4)),"%, (",(sum((knn_result_5_3 == pileupsAllShuffled[1:(nrow(pileupsAllShuffled)/2)]$VAR)==FALSE,na.rm=TRUE)),"/",(nrow(pileupsAllShuffled)/2),")."))

confusionMatrix(knn_result_5_3,as.factor(unlist(pileupsAllShuffled[1:(nrow(pileupsAllShuffled)/2),"VAR"])))

set.seed(seed)
fit.lda <- train(VAR~., data=pileupsAllShuffled[,c("DP","RD","VRF","bias","VAR")], method="lda", metric=metric, preProc=c("center", "scale"), trControl=control)
set.seed(seed)
fit.knn <- train(VAR~., data=pileupsAllShuffled[,c("DP","RD","VRF","bias","VAR")], method="knn", metric=metric, preProc=c("center", "scale"), trControl=control)
results <- resamples(list(lda=fit.lda, knn=fit.knn))
# Table comparison
summary(results)

# boxplot comparison
bwplot(results)
# Dot-plot comparison
dotplot(results)

