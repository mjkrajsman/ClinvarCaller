loadModel <- function(folder, file, pileups){
	model <- readRDS(paste0(folder,file))
	model$pred$obs <- as.character(model$pred$obs)
	model$pred$pred <- as.character(model$pred$pred)
	model_final <- data.frame(cbind(model$pred, pileups[model$pred$rowIndex,]))
	model_final
}

initializeSubsetDpBins <- function(lowerBounds = c(1,1,2,6,11,21,51,101,201), upperBounds = c(9999,1,5,10,20,50,100,200,9999)){
  out <- tryCatch(
    {
      if(length(lowerBounds)!=length(upperBounds)){
        stop("lowerBounds and upperBounds have different sizes!")
      }
      if(!all(lowerBounds<=upperBounds)){
        stop("All values in lowerBounds have to be lower than corresponding values in upperBounds!")
      }
      subsetDpDataTable <- data.table("lower"=lowerBounds,
                                      "upper"=upperBounds)		
    },
    error=function(cond) {
      message("Something went wrong. Function: initializeSubsetDpBins")
      message("This function creates a subsetDpBins data.table, which contains lower and upper values of DP ranges used in calculations.")
      message("Returns: data.table with DP values.")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
      return(subsetDpDataTable)
    })
  return(out)
}



getScores <- function(res, dpMin, dpMax, modelDesc, paramName =NULL)
{
  tpr <- length(which(res$pred == "variant"  & res$VAR == "variant" & res$DP <= dpMax  & res$DP >= dpMin))/ length(which(res$VAR == "variant" & res$DP <= dpMax & res$DP >= dpMin))
  fpr <- length(which(res$pred == "variant"  & res$VAR == "nonVariant" & res$DP <= dpMax  & res$DP >= dpMin))/ length(which(res$VAR == "nonVariant" & res$DP <= dpMax & res$DP >= dpMin))
  tpr_baseline <- length(which(res$AllNonRefCount >0  & res$VAR == "variant" & res$DP <= dpMax  & res$DP >= dpMin))/ length(which(res$VAR == "variant" & res$DP <= dpMax & res$DP >= dpMin))
  fpr_baseline <- length(which(res$AllNonRefCount >0 & res$VAR == "nonVariant" & res$DP <= dpMax  & res$DP >= dpMin))/ length(which(res$VAR == "nonVariant" & res$DP <= dpMax & res$DP >= dpMin))
  
  tpr_imp <- 100*((tpr - tpr_baseline) / tpr_baseline)
  fpr_imp <- 100*(-(fpr - fpr_baseline) / fpr_baseline)
  
  tprTab <- tpr
  fprTab <- fpr
  
  if (!is.na(paramName)){   
    tprTab <- lapply(unique(res[,paramName]), function(x){
      length(which(res[,paramName] == x & res$pred == "variant"  & res$VAR == "variant" & res$DP <= dpMax  & res$DP >= dpMin))/ length(which(res$VAR == "variant" & res$DP <= dpMax & res$DP >= dpMin & res[,paramName] == x))
    })

    fprTab <- lapply(unique(res[,paramName]), function(x){
      length(which(res[,paramName] == x & res$pred == "variant"  & res$VAR == "nonVariant" & res$DP <= dpMax  & res$DP >= dpMin))/ length(which(res$VAR == "nonVariant" & res$DP <= dpMax & res$DP >= dpMin  & res[,paramName] == x))
    })
   
  }
  
  tpr_imp_min <- 100*((min(unlist(tprTab)) - tpr_baseline) / tpr_baseline)
  tpr_imp_max <- 100*((max(unlist(tprTab)) - tpr_baseline) / tpr_baseline)
  fpr_imp_min <- 100*(-(min(unlist(fprTab)) - fpr_baseline) / fpr_baseline)
  fpr_imp_max <- 100*(-(max(unlist(fprTab)) - fpr_baseline) / fpr_baseline)

  
  nsample <-  length(which( res$DP <= dpMax  & res$DP >= dpMin))
  c(modelDesc=modelDesc,bin = paste0(dpMin, "-", dpMax),
    nsample=nsample,tpr=tpr, fpr=fpr, tpr_baseline=tpr_baseline,
    fpr_baseline=fpr_baseline , tpr_imp= tpr_imp, fpr_imp=fpr_imp,
   tpr_imp_min=tpr_imp_min, tpr_imp_max=tpr_imp_max,
   fpr_imp_min=fpr_imp_min,fpr_imp_max=fpr_imp_max)
}



getFinalScores <- function(model, bins,  modelDesc, paramName=NULL){
  do.call(rbind, mclapply(1:nrow(bins), function(i){
   getScores(model, bins$lower[i],bins$upper[i], as.character(modelDesc), as.character(paramName))
  }, mc.cores=4))
}



library(parallel)
pileups <- readRDS("pileupsAllShuffled.rds")
folder <- "./model_benchmark/"
files <- dir(folder)
#ff <- files[-grep("pileups", files )]
ff <- files
models <- lapply(files, function(x){loadModel(folder, x, pileups)})
bins <- initializeSubsetDpBins(lowerBounds = c(1,1,2,6,11,21,51,101,201), upperBounds = c(9999,1,5,10,20,50,100,200,9999))

names(models) <- ff
names(models) <- files
modelSig <- data.frame(do.call(rbind,strsplit(ff, "[_\\.]")))
colnames(modelSig) <- c("modelName", "center", "zv", "pca", "pca_thresh", "cv", "rds")
modelName2paramName_modelDesc <- data.frame(rbind(c("knn", "KNN", "k"), 
					   c("lda", "LDA", NA),
					   c("nb", "NB", "usekernel"),
					   c("rf", "RF", "mtry"),
					   c("rpart", "RPART", "cp")))

colnames(modelName2paramName_modelDesc) <- c("modelName", "modelDesc", "paramName")

finalScores <- do.call(rbind, lapply(1:length(models),  function(i){
	model <- models[[i]]
	paramName <- modelName2paramName_modelDesc$paramName[which(modelName2paramName_modelDesc$modelName ==  modelSig$modelName[i])]
	modelDesc <- modelName2paramName_modelDesc$modelDesc[which(modelName2paramName_modelDesc$modelName ==  modelSig$modelName[i])]
	print(modelDesc)
	mfs <- getFinalScores(model, bins, modelDesc, paramName)
	mfs
}))
       
       
#finalScores <- data.frame(rbind(,  
#                                getFinalScores(ld_final, bins, "LDA")), stringsAsFactors=F)

library(reshape2)
evalMelt <- melt(data.frame(finalScores, stringsAsFactors=F), id.vars=c("bin", "modelDesc", "nsample", "tpr_imp_min", "tpr_imp_max","fpr_imp_min", "fpr_imp_max", "tpr", "fpr", "tpr_baseline", "fpr_baseline"),variable.name = "metric")
evalMelt$imp_min <- evalMelt$fpr_imp_min 
evalMelt$imp_max <- evalMelt$fpr_imp_max 
evalMelt$imp_min  [which(evalMelt$metric == "tpr_imp")]<- evalMelt$tpr_imp_max [which(evalMelt$metric == "tpr_imp")]
evalMelt$imp_max  [which(evalMelt$metric == "tpr_imp")]<- evalMelt$tpr_imp_max [which(evalMelt$metric == "tpr_imp")]
evalMelt$imp_min <- as.numeric(evalMelt$imp_min)
evalMelt$imp_max <- as.numeric(evalMelt$imp_max)

library(ggplot2)
library(plyr)
evalMelt$bin <- factor(evalMelt$bin, levels=c("1-9999","1-1","2-5","6-10","11-20","21-50","51-100","101-200","201-9999"))
evalMelt$metric <- factor(evalMelt$metric, levels=c("fpr_imp","tpr_imp"))
evalMelt$modelDesc <- factor(evalMelt$modelDesc, levels=c("LDA", "KNN", "RF", "RPART", "NB"))

evalMelt$value <- as.numeric(evalMelt$value)
plotTitle <- "Comparison of classification methods"
evalPlot <- ggplot(data = evalMelt, aes(bin, value, fill=metric))  +
  geom_hline(yintercept=100, linetype="dashed", color="#808080") +
  geom_hline(yintercept=-100, linetype="dashed", color="#808080") +
  geom_hline(yintercept=0) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes( ymax = imp_max, ymin=imp_min), position="dodge", colour="#808080")+
  facet_grid(modelDesc~., scales="free") +
  theme(axis.text.x = element_text(face="bold", color="#993333", angle=45),
  axis.text.y = element_text(face="bold", color="#993333", angle=45)) + 
  geom_text(aes(label=round(value,digits = 2),fontface="bold"), vjust=-0.3, position=position_dodge(1.0), size=2.5) +
  scale_fill_manual(values = c("#0000FF","#FF0000")) +
  ylab("value [%]") +
  xlab("DP subset") +
  ggtitle(str_replace_all(plotTitle," ","_"))


evalPlotScaled <- evalPlot +  coord_cartesian(ylim=c(-10,110)) + scale_y_continuous(breaks = c(0,25,50,75,100))
  


evalPlot
evalPlotScaled
rm(models)


plotFolderPath <- "./output/plots/"
if(!dir.exists(plotFolderPath)) dir.create(plotFolderPath)
pdf(paste0(plotFolderPath,str_replace_all(plotTitle," ","_"),".pdf"))
plot(evalPlot)
dev.off()
if(verbose==TRUE){print(paste0("Plot ",plotFolderPath,str_replace_all(plotTitle," ","_"),".pdf saved."))}  
if(!dir.exists(plotFolderPath)) dir.create(plotFolderPath)
pdf(paste0(plotFolderPath,str_replace_all(plotTitle," ","_"),"Scaled",".pdf"))
plot(evalPlotScaled)
dev.off()
if(verbose==TRUE){print(paste0("Plot ",plotFolderPath,str_replace_all(plotTitle," ","_"),"Scaled",".pdf saved."))}  

