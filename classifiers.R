library(data.table)
library(Rsamtools)
library(tictoc)
library(stringr)
library(testthat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(class)
library(mlbench)
library(caret)
library(plotROC)
library(doMC)




#registerDoMC(cores=0)
hg38=TRUE
cvTest <- clinvarCaller(newCoverage=FALSE,newPileups=TRUE,newFalsePileups=TRUE,newPositionsAllFalse=TRUE,
                        newPileupsWithStats=TRUE,newFalsePileupsWithStats=TRUE,
                        newVcfHcGeneral=FALSE, newVcfHcFc=FALSE,runsAtImid=FALSE,hg38=TRUE,
                        returnMergedDb=FALSE,returnStats=FALSE,returnPileups=TRUE,returnPositionsAll=TRUE,returnPositionsAllFalse=TRUE,
                        returnHistogramData=TRUE,verbose=TRUE,minBaseQuality=0,minMapq=0,fullPileupStats=TRUE,
                        chromosomeLengthsPath="./input/chromosomeLengths.csv",
                        bamPath="./input/corriell_S7-ready.bam",
                        #inputVariantsDbPath="./input/clinvar_alleles.single.b38.tsv.gz",
                        inputVariantsDbPath="./input/cvSample.tsv",
                        #inputVariantsDbPath="./input/goldenStandard38.vcf.gz",
                        vcfHcGeneralPath="./input/corriell_S7-HCScan.vcf.gz",
                        vcfHcFcPath="./input/corriell_S7-HCScanFC.vcf.gz",
                        outputCoveragePath="./output/corriell_S7-ready",
                        outputPositionsAllPath="./output/positionsAll1.positions",
                        outputPositionsAllFalsePath="./output/positionsAllFalse1.positions",
                        outputPileupsPath="./output/pileups1.mpileup",
                        outputPileupsFalsePath="./output/pileupsFalse1.mpileup",
                        outputPileupsWithStatsPath="./output/pileupsWithStats1.csv",
                        outputFalsePileupsWithStatsPath="./output/falsePileupsWithStats1.csv",
                        outputStatsPath="./output/stats1.csv",
                        outputVarDbMergedPath="./output/varDbMerged-reduced1.csv")

#rozne progi

#knn_result_5_3 <- knn(pileupsAllShuffled[(1+nrow(pileupsAllShuffled)/2):(nrow(pileupsAllShuffled)),c("DP","RD","VRF","bias")], pileupsAllShuffled[1:(nrow(pileupsAllShuffled)/2),c("DP","RD","VRF","bias")], pileupsAllShuffled[(1+nrow(pileupsAllShuffled)/2):(nrow(pileupsAllShuffled))]$VAR, k=5, l=3)
#print("k=5, l=3")
#print(paste0("TRUE: ",(round(sum((knn_result_5_3 == pileupsAllShuffled[1:(nrow(pileupsAllShuffled)/2)]$VAR)==TRUE,na.rm=TRUE)/(nrow(pileupsAllShuffled)/2),digits=4)),"%, (",(sum((knn_result_5_3 == pileupsAllShuffled[1:(nrow(pileupsAllShuffled)/2)]$VAR)==TRUE,na.rm=TRUE)),"/",(nrow(pileupsAllShuffled)/2),")."))
#print(paste0("FALSE: ",(round(sum((knn_result_5_3 == pileupsAllShuffled[1:(nrow(pileupsAllShuffled)/2)]$VAR)==FALSE,na.rm=TRUE)/(nrow(pileupsAllShuffled)/2),digits=4)),"%, (",(sum((knn_result_5_3 == pileupsAllShuffled[1:(nrow(pileupsAllShuffled)/2)]$VAR)==FALSE,na.rm=TRUE)),"/",(nrow(pileupsAllShuffled)/2),")."))
#confusionMatrix(knn_result_5_3,as.factor(unlist(pileupsAllShuffled[1:(nrow(pileupsAllShuffled)/2),"VAR"])))





#dataSubset <- subset(pileupsAllShuffled, (DP>=20 & DP<=50))

performTraining <- function(formula=NULL, data=NULL, metric=NULL, control=NULL, seed=7, mc=0){
  registerDoMC(cores=mc)
  set.seed(seed)
  preProc=c("center","zv","scale")
  fit.lda <- train(formula,data=data, method="lda", metric=metric, preProc=preProc, trControl=control)
  set.seed(seed)
  preProc=c("center","zv","pca")
  fit.ldaPca <- train(formula,data=data, method="lda", metric=metric, preProc=preProc, trControl=control)
  set.seed(seed)
  preProc=c("center","zv","scale")
  fit.rpart <- train(formula,data=data, method="rpart", metric=metric, preProc=preProc, trControl=control)
  set.seed(seed)
  preProc=c("center","scale")
  fit.rpartV <- train(formula,data=data, method="rpart", metric=metric, preProc=preProc, trControl=control)
  set.seed(seed)
  preProc=c("center","zv","scale")
  fit.knn <- train(formula,data=data, method="knn", metric=metric, preProc=preProc, trControl=control)
  set.seed(seed)
  preProc=c("center","zv","pca")
  fit.knnPca <- train(formula,data=data, method="knn", metric=metric, preProc=preProc, trControl=control)
  set.seed(seed)
  preProc=c("center","zv","scale")
  fit.nb <- train(formula,data=data, method="nb", metric=metric, preProc=preProc, trControl=control)
  set.seed(seed)
  preProc=c("center","zv","ica")
  fit.nbIca <- train(formula,data=data, method="nb", metric=metric, preProc=preProc, trControl=control)
  set.seed(seed)
  preProc=c("center","zv","pca")
  fit.nbPca <- train(formula,data=data, method="nb", metric=metric, preProc=preProc, trControl=control)
  set.seed(seed)
  preProc=c("center","zv","scale")
  fit.rf <- train(formula,data=data, method="rf", metric=metric, preProc=preProc, trControl=control)
  set.seed(seed)
  preProc=c("center","zv","scale")
  fit.svmRadial <- train(formula,data=data, method="svmRadial", metric=metric, preProc=preProc, trControl=control)
  set.seed(seed)
  preProc=c("center","zv","pca")
  fit.svmRadialPca <- train(formula,data=data, method="svmRadial", metric=metric, preProc=preProc, trControl=control)

  
  models <- list(lda=fit.lda, ldaPca=fit.ldaPca, rpart=fit.rpart, rpartPca=fit.rpartPca, rpartV=fit.rpartV, 
                 knn=fit.knn,knn=fit.knn,nb=fit.nb,nbIca=fit.nbIca,nbPca=fit.nbPca,rf=fit.rf,
                 fit.svmRadial=fit.svmRadial,fit.svmRadialPca=fit.svmRadialPca)
  
  results <- resamples(models)
  
  res <- list()
  registerDoMC(cores=0)
  return(append(res,list(models=models, results=resamples(models), summary=summary(resamples(models)))))
}

createROC <- function(trainObject=NULL, title=""){
  roc <- ggplot(trainObject$pred, aes(m=variant, d=factor(obs, levels = c("nonVariant", "variant")))) + 
    geom_roc(n.cuts=0) + 
    coord_equal() +
    style_roc() +
    ggtitle(title)
    # + annotate("text", x=0.75, y=0.25, label=paste("AUC =", round((calc_auc(g))$AUC, 4)))
  return(roc)
}

createCombinedROC <- function(listOfTrainObjects=NULL, title="", method=""){
  listOfPreds <- lapply(seq_along(listOfTrainObjects), function(n){
    #pred <- as.data.table(listOfTrainObjects[[n]]$models$lda$pred[,c("nonVariant", "variant", "obs")])
    #View(res1[["models"]][["lda"]][["pred"]])
    pred <- as.data.table(listOfTrainObjects[[n]][["models"]][[method]][["pred"]][,c("nonVariant", "variant", "obs")])
    pred[,"series":=paste(names(listOfTrainObjects)[[n]])]
  })
  boundData <- rbindlist(listOfPreds)
  #View(boundData)
  roc <- ggplot(boundData, aes(m=variant, d=obs, color=series)) + 
    geom_roc(n.cuts=0) + 
    coord_equal() +
    style_roc() + 
    ggtitle(title)
  return(roc)
}
  

pileups2 <- copy(cv$pileups)
pileupsRandom2 <- copy(cv$falsePileups)
pileups2[,"VAR":=1]
pileupsRandom2[,"VAR":=0]
pileupsAll <- rbind(pileups2,pileupsRandom2)
pileupsAll <- addRefFromReferenceGenome(pileupsAll)
pileupsAllShuffled <- pileupsAll[sample(nrow(pileupsAll)),]
rm(pileups2,pileupsRandom2,pileupsAll)


metric <- "ROC"
preProc=c("center","zv","pca")
control <- trainControl(method="repeatedcv", summaryFunction=twoClassSummary, classProbs=TRUE, savePredictions = TRUE, number=5, repeats=3, allowParallel = TRUE) #or list(pcaComp = 7))
seed <- 7
formula <- as.formula(c("VAR~", paste(names(dataSubset[,!c("CHROM","POS","NUC","BQ","MQ","PR","VAR","REF")]), collapse = "+")))

#all
res_all <- performTraining(formula=formula,data=pileupsAllShuffled,metric=metric,control=control,seed=7)
#1
res1 <- performTraining(formula=formula,data=subset(pileupsAllShuffled,DP==1),metric=metric,control=control,seed=7)
roc1 <- createROC(trainObject = res1$models$lda, title="lda, DP: [1]")
#2-5
res2_5 <- performTraining(formula=formula,data=subset(pileupsAllShuffled,(DP>=2 & DP<=5)),metric=metric,control=control,seed=7)
roc2_5 <- createROC(trainObject = res2_5$models$lda, title="lda, DP: [2:5]")
#6-10
res6_10 <- performTraining(formula=formula,data=subset(pileupsAllShuffled,(DP>=6 & DP<=10)),metric=metric,control=control,seed=7)
roc6_10 <- createROC(trainObject = res6_10$models$lda, title="lda, DP: [6:10]")
#11-20
res11_20 <- performTraining(formula=formula,data=subset(pileupsAllShuffled,(DP>=11 & DP<=20)),metric=metric,control=control,seed=7)
roc11_20 <- createROC(trainObject = res11_20$models$lda, title="lda, DP: [11:20]")
#20-50
res21_50 <- performTraining(formula=formula,data=subset(pileupsAllShuffled,(DP>=21 & DP<=50)),metric=metric,control=control,seed=7)
roc21_50 <- createROC(trainObject = res21_50$models$lda, title="lda, DP: [21:50]")
#51-100
res51_100 <- performTraining(formula=formula,data=subset(pileupsAllShuffled,(DP>=51 & DP<=100)),metric=metric,control=control,seed=7)
roc51_100 <- createROC(trainObject = res51_100$models$lda, title="lda, DP: [51:100]")
#101-200
res101_200 <- performTraining(formula=formula,data=subset(pileupsAllShuffled,(DP>=101 & DP<=200)),metric=metric,control=control,seed=7)
roc101_200 <- createROC(trainObject = res101_200$models$lda, title="lda, DP: [101:200]")
#201+
res201 <- performTraining(formula=formula,data=subset(pileupsAllShuffled,DP>=201),metric=metric,control=control,seed=7)
roc201 <- createROC(trainObject = res201$models$lda, title="lda, DP: [201+]")

#rm(list = ls(all.names = TRUE, pattern = "res.*"))



metric <- "ROC"
preProc=c("center","zv","pca")
control <- trainControl(method="repeatedcv", summaryFunction=twoClassSummary, classProbs=TRUE, savePredictions = TRUE, number=5, repeats=3) #or list(pcaComp = 7))
seed <- 7
formula <- as.formula(c("VAR~", paste(names(dataSubset[,!c("CHROM","POS","DP","NUC","BQ","MQ","PR","VAR")]), collapse = "+")))

tic()
set.seed(seed)
fit.lda <- train(formula,data=pileupsAllShuffled, method="lda", 
                    metric=metric, preProc=c("center","zv","scale"), trControl=
                      trainControl(method="repeatedcv", summaryFunction=twoClassSummary, classProbs=TRUE, savePredictions = TRUE, number=5, repeats=3))
toc()
#rocLda <- createROC(trainObject = fit.lda, title="lda")
tic()
set.seed(seed)
fit.ldaPca <- train(formula,data=pileupsAllShuffled, method="lda", 
                 metric=metric, preProc=c("center","zv","pca"), trControl=
                   trainControl(method="repeatedcv", summaryFunction=twoClassSummary, classProbs=TRUE, savePredictions = TRUE, number=5, repeats=3))
toc()
#rocLdaPca <- createROC(trainObject = fit.ldaPca, title="ldaPca")


methodsList <- list(lda=fit.lda, ldaPca=fit.ldaPca)

results <- resamples(methodsList)
# Table comparison
summary(results)


tic()
set.seed(seed)
fit.rpart2_5 <- train(formula,data=subset(pileupsAllShuffled,(DP>=2 & DP<=5)), method="rpart", 
                      metric=metric, preProc=c("center","zv","scale"), trControl=
                        trainControl(method="repeatedcv", summaryFunction=twoClassSummary, classProbs=TRUE, savePredictions = TRUE, number=5, repeats=3, allowParallel = TRUE))

toc()
tic()
set.seed(seed)
fit.knn2_5 <- train(formula,data=subset(pileupsAllShuffled,(DP>=2 & DP<=5)), method="knn", 
                    metric=metric, preProc=c("center","zv","scale"), trControl=
                      trainControl(method="repeatedcv", summaryFunction=twoClassSummary, classProbs=TRUE, savePredictions = TRUE, number=5, repeats=3, allowParallel = TRUE))
#roc2_5knn <- createROC(trainObject = fit.knn2_5, title="knn, DP: [2:5]")
toc()
tic()
set.seed(seed)
fit.knn2_5pca <- train(formula,data=subset(pileupsAllShuffled,(DP>=2 & DP<=5)), method="knn", 
                       metric=metric, preProc=c("center","zv","pca"), trControl=
                         trainControl(method="repeatedcv", summaryFunction=twoClassSummary, classProbs=TRUE, savePredictions = TRUE, number=5, repeats=3, allowParallel = TRUE))
#roc2_5knnPca <- createROC(trainObject = fit.knn2_5pca, title="knnPCA, DP: [2:5]")
toc()
tic()
set.seed(seed)
fit.nb2_5 <- train(formula,data=subset(pileupsAllShuffled,(DP>=2 & DP<=5)), method="nb", 
                    metric=metric, preProc=c("center","zv", "ica"), trControl=
                     trainControl(method="repeatedcv", summaryFunction=twoClassSummary, classProbs=TRUE, savePredictions = TRUE, number=5, repeats=3, allowParallel = TRUE))
#roc2_5nb <- createROC(trainObject = fit.nb2_5, title="nb, DP: [2:5]")
toc()
tic()
set.seed(seed)
fit.nb2_5pca <- train(formula,data=subset(pileupsAllShuffled,(DP>=2 & DP<=5)), method="nb", 
                       metric=metric, preProc=c("center","zv","pca"), trControl=
                        trainControl(method="repeatedcv", summaryFunction=twoClassSummary, classProbs=TRUE, savePredictions = TRUE, number=5, repeats=3, allowParallel = TRUE))
toc()
tic()
set.seed(seed)
fit.rf2_5 <- train(formula,data=subset(pileupsAllShuffled,(DP>=2 & DP<=5)), method="rf", 
                   metric=metric, preProc=c("center", "scale", "zv"), trControl=
                     trainControl(method="repeatedcv", summaryFunction=twoClassSummary, classProbs=TRUE, savePredictions = TRUE, number=5, repeats=3, allowParallel = TRUE))
#roc2_5rf <- createROC(trainObject = fit.rf2_5, title="rf, DP: [2:5]")
toc()
tic()
set.seed(seed)
fit.svmRadial2_5 <- train(formula,data=subset(pileupsAllShuffled,(DP>=2 & DP<=5)), method="svmRadial", 
                   metric=metric, preProc=c("center","scale","zv"), trControl=
                     trainControl(method="repeatedcv", summaryFunction=twoClassSummary, classProbs=TRUE, savePredictions = TRUE, number=5, repeats=3, allowParallel = TRUE))
#roc2_5svmRadial <- createROC(trainObject = fit.svmRadial2_5, title="SVMRadial, DP: [2:5]")
toc()
tic()
set.seed(seed)
fit.svmRadial2_5pca <- train(formula,data=subset(pileupsAllShuffled,(DP>=2 & DP<=5)), method="svmRadial", 
                      metric=metric, preProc=c("center","pca","zv"), trControl=
                        trainControl(method="repeatedcv", summaryFunction=twoClassSummary, classProbs=TRUE, savePredictions = TRUE, number=5, repeats=3, allowParallel = TRUE))
#roc2_5svmRadialPca <- createROC(trainObject = fit.svmRadial2_5pca, title="SVMRadialPCA, DP: [2:5]")
toc()





methodsList <- list(lda2_5=fit.lda2_5, lda2_5pca=fit.lda2_5pca, rpart2_5=fit.rpart2_5, rpart2_5pca=fit.rpart2_5pca,
                    knn2_5=fit.knn2_5, knn2_5pca=fit.knn2_5pca, nb2_5=fit.nb2_5, nb2_5pca=fit.nb2_5pca,
                    rf2_5=fit.rf2_5, svmRadial2_5=fit.svmRadial2_5, svmRadial2_5pca=fit.svmRadial2_5pca)

results <- resamples(methodsList)
# Table comparison
summary(results)


#g <- subset(pileupsAllShuffled,(DP==1))
#g[,"predict":=predict(fit.lda2_5,newdata = subset(pileupsAllShuffled,(DP==1)))]



                      
#set.seed(seed)
#fit.knn <- train(VAR~., data=dataSubset[,!c("CHROM","POS","DP","NUC","BQ","MQ","PR")], method="knn", metric=metric, preProc=preProc, trControl=control)
#set.seed(seed)
#fit.nb <- train(VAR~., data=dataSubset[,!c("CHROM","POS","DP","NUC","BQ","MQ","PR")], method="nb", metric=metric, preProc=preProc, trControl=control)
#set.seed(seed)
#fit.rf <- train(VAR~., data=dataSubset[,!c("CHROM","POS","DP","NUC","BQ","MQ","PR")], method="rf", metric=metric, preProc=preProc, trControl=control)
#set.seed(seed)
#fit.svmRadial <- train(VAR~., data=dataSubset[,!c("CHROM","POS","DP","NUC","BQ","MQ","PR")], method="svmRadial", metric=metric, preProc=preProc, trControl=control)

#results <- resamples(list(lda=fit.lda, knn=fit.knn, nb=fit.nb, svmRadial=fit.svmRadial))
# Table comparison
#summary(results)

## boxplot comparison
#bwplot(results)
## Dot-plot comparison
#dotplot(results)


l <- createCombinedROC(listOfTrainObjects = list(res001=res1,res002_005=res2_5,res006_010=res6_10,res011_020=res11_20,
                                                 res021_050=res21_50,res051_100=res51_100,res101_200=res101_200,
                                                 res201=res201), title="lda", method="lda")

k <- createCombinedROC(listOfTrainObjects = list(res001=res1,res002_005=res2_5,res006_010=res6_10,res011_020=res11_20,
                                                 res021_050=res21_50,res051_100=res51_100,res101_200=res101_200,
                                                 res201=res201), title="knn", method="knn")




# naive approach:
# pas$AllNonRefCount > 0 -> variant observed
# pas$AllNonRefCount == 0 -> no variant observed
# pas$VAR == "variant" -> variant is present
# pas$VAR == "nonVariant" -> variant is not present

# TP = pas$AllNonRefCount  > 0 & pas$VAR == "variant"
# FN = pas$AllNonRefCount == 0 & pas$VAR == "variant"
# TN = pas$AllNonRefCount == 0 & pas$VAR == "nonVariant"
# FP = pas$AllNonRefCount  > 0 & pas$VAR == "nonVariant"

# TP + FN = pas$VAR == "variant"
# TN + FP = pas$VAR == "nonVariant"

#TODO: try-catch
getNaiveBaseline <- function(subsetBounds=NULL, data=NULL){
  pas <- data
  res <- lapply(1:nrow(subsetBounds), function(i){
    lower <- subsetBounds[i,lower]
    upper <- subsetBounds[i,upper]
    params <- list(
      lower <- lower,
      upper <- upper,
      # TP = pas$AllNonRefCount  > 0 & pas$VAR == "variant"
      TP = length(which(pas$AllNonRefCount  > 0 & pas$VAR == "variant" & pas$DP >= lower & pas$DP <= upper)),
      
      # FN = pas$AllNonRefCount == 0 & pas$VAR == "variant"
      FN = length(which(pas$AllNonRefCount == 0 & pas$VAR == "variant" & pas$DP >= lower & pas$DP <= upper)),
      
      # TN = pas$AllNonRefCount == 0 & pas$VAR == "nonVariant"
      TN = length(which(pas$AllNonRefCount == 0 & pas$VAR == "nonVariant" & pas$DP >= lower & pas$DP <= upper)),
      
      # FP = pas$AllNonRefCount  > 0 & pas$VAR == "nonVariant"
      FP = length(which(pas$AllNonRefCount  > 0 & pas$VAR == "nonVariant" & pas$DP >= lower & pas$DP <= upper)),
      
      # P = TP + FN = pas$VAR == "variant"
      P = length(which(pas$VAR == "variant" & pas$DP >= lower & pas$DP <= upper)),
      
      # N = TN + FP = pas$VAR == "nonVariant"
      N = length(which(pas$VAR == "nonVariant" & pas$DP >= lower & pas$DP <= upper)),
      
      # TPR = sensitivity = TP/P = TP/(TP+FN)
      TPR = length(which(pas$AllNonRefCount  > 0 & pas$VAR == "variant" & pas$DP >= lower & pas$DP <= upper)) / length(which(pas$VAR == "variant" & pas$DP >= lower & pas$DP <= upper)),
      # TNR = specificity = TN/N = TN/(TN+FP) = 1-FPR = 1-FP/(TN+FP)
      # FPR = fall-out = FP/N = FP/(TN+FP)
      FPR = length(which(pas$AllNonRefCount  > 0 & pas$VAR == "nonVariant" & pas$DP >= lower & pas$DP <= upper)) / length(which(pas$VAR == "nonVariant" & pas$DP >= lower & pas$DP <= upper))
    )
    return(params)
  })
  res <- rbindlist(res)
  return(res) 
}
