library(data.table)
library(Rsamtools)
library(tictoc)
library(stringr)
library(testthat)
library(tools)
library(BSgenome.Hsapiens.UCSC.hg38)
#library(BSgenome.Hsapiens.UCSC.hg19)
library(class)
library(mlbench)
library(ggplot2)
library(caret)
library(plotROC)
library(doMC)
#source("./tests.R")

#TODO: optimize new(...)=TRUE/FALSE usage
#TODO: optional: multi-threading?
clinvarCaller <- function(newCoverage=TRUE,newPileups=TRUE,newFalsePileups=TRUE,newPostitionsAll=FALSE,
                          newPositionsAllFalse=TRUE,newPileupsWithStats=TRUE,newFalsePileupsWithStats=TRUE,
                          newVcfHcGeneral=FALSE, newVcfHcFc=FALSE,runsAtImid=FALSE,hg38=TRUE,
                          returnMergedDb=FALSE,returnStats=TRUE,returnPileups=TRUE,returnPositionsAll=TRUE,returnPositionsAllFalse=TRUE,
                          returnHistogramData=TRUE,verbose=TRUE,minBaseQuality=0,minMapq=0,fullPileupStats=TRUE,
                          train=FALSE,HCTest=TRUE, mc=0,
						  lowerSubsetBounds = c(1,2,6,11,21,51,101,201), upperSubsetBounds = c(1,5,10,20,50,100,200,9999),
                          chromosomeLengthsPath="./input/chromosomeLengths.csv",
                          bamPath="./input/corriell_S7-ready.bam",
                          inputVariantsDbPath="./input/clinvar_alleles.single.b38.tsv.gz",
                          vcfHcGeneralPath="./input/corriell_S7-HCScan.vcf.gz",
                          vcfHcFcPath="./input/corriell_S7-HCScanFC.vcf.gz",
                          modelPath="./input/model.rds",
                          listOfTrainArguments=list(method="lda",
                                                       preProc=c("center","zv","scale"),
                                                       control=trainControl(method="repeatedcv", summaryFunction=twoClassSummary, classProbs=TRUE, 
                                                                            savePredictions = TRUE, number=5, repeats=3, allowParallel = TRUE), 
                                                       formula=as.formula(c("VAR~", paste(names(dataSubset[,!c("CHROM","POS","NUC","BQ","MQ","PR","VAR","REF")]), 
                                                                                          collapse = "+")))     
                          ),
                          outputCoveragePath="./output/corriell_S7-ready",
                          outputPositionsAllPath="./output/positionsAll.positions",
                          outputPositionsAllFalsePath="./output/positionsAllFalse.positions",
                          outputPileupsPath="./output/pileups.mpileup",
                          outputPileupsFalsePath="./output/pileupsFalse.mpileup",
                          outputPileupsWithStatsPath="./output/pileupsWithStats.csv",
                          outputFalsePileupsWithStatsPath="./output/falsePileupsWithStats.csv",
                          outputStatsPath="./output/stats.csv",
                          outputVarDbMergedPath="./output/varDbMerged-reduced.csv"){
  #TODO: make returnMergedDb true in release

if(hg38==TRUE){
  library(BSgenome.Hsapiens.UCSC.hg38)
}else{
  library(BSgenome.Hsapiens.UCSC.hg19)
}

#TODO: try-catch
gunzipPath <- function(string){
  
  if(file_ext(string)=="gz"){
    res <- paste0("gunzip -c ",string)
  }else{
    res <- string
  }
  return(res)
}


#TODO: try-catch
#gunzipPath
loadBamFile <- function(bamPath=""){
	# BAM
	test_that("Does BAM file exist?", { expect_that(file.exists(bamPath), is_true()) })
	bamPath <- gunzipPath(bamPath)
	bamFile <- BamFile(bamPath)
	if(verbose==TRUE){print("BAM file loaded.")}
	#seqinfo(bamFile)
	return(bamFile)
}

#TODO: try-catch, gunzip-path?
getIndexFile <- function(bamPath=""){
	if(exists(paste0(bamPath, ".bai"))){
	  indexFile <- paste0(bamPath, ".bai")
	  if(verbose==TRUE){print("Index (BAI) file loaded.")}
	}else{
	  indexFile <- indexBam(bamFile)
	  if(verbose==TRUE){print("Index (BAI) file created and loaded.")}
	}
	return(indexFile)
}

#TODO: try-catch
makeBed <- function(coveragePath=""){
  cmd <- paste0("mosdepth '",coveragePath,"' ",bamPath)
  system(cmd)
  return(paste0(coveragePath,".per-base.bed.gz"))
}

#TODO: try-catch
filterChroms <- function(input){
  res <- input[which(input$CHROM %in% c(1:22, "X", "Y", "M")),]
  return(res)
}

#TODO: try-catch
#filterChroms
getCoverageFromBed <- function(bedPath=""){
  coverage <- fread(bedPath)
  colnames(coverage) <- c("CHROM","start","stop","readsCount")
  coverage <- subset(coverage, readsCount > 0)
  coverage$CHROM <- gsub("chr", "", coverage$CHROM)
  coverage <- filterChroms(coverage)
  coverage[,"start":=start+1L]
  #coverage[,"length":=stop-start]
  return(coverage)
}

#TODO: try-catch
#gunzipPath, makeBed, getCoverageFromBed, filterChroms
getCoverage <- function(outputCoveragePath=outputCoveragePath, new=TRUE){
	# BED
	if(new==TRUE){
		bedPath <- makeBed(outputCoveragePath)
		if(verbose==TRUE){print("Coverage file (BED) created.")}
	}else{
	  bedPath <- paste0(outputCoveragePath,".per-base.bed.gz")
	}
  
	test_that("Does BED file exist?", { expect_that(file.exists(bedPath), is_true())})
	bedPath <- gunzipPath(bedPath)
	coverage <- getCoverageFromBed(bedPath)
	if(verbose==TRUE){print("BED file loaded.")}
	return(coverage)
}

#TODO: try-catch
getFileExtension <- function(string, removeLastExt = FALSE){
  if(removeLastExt==TRUE & str_count(str_remove(string,"[.]{1,}"), "\\.")>1){
    string <- file_path_sans_ext(string)
  }
  res <- file_ext(string)
  return(res)
}

#TODO: try-catch, merge with loadTsvData?
loadVcfData <- function(vcfPath=""){
  res <- fread(vcfPath, skip = "CHROM")
  colnames(res)[colnames(res)=="#CHROM"] <- "CHROM"
  res$CHROM <- gsub("chr", "", res$CHROM)
  res <- res[!like(ALT,",")] # remove multiallelic sites
  res <- subset(unique(res, by=c("CHROM","POS")))
  res[CHROM == "MT" ,"CHROM":="M"]
  return(res)
}

#TODO: try-catch
loadTsvData <- function(tsvPath=""){
  res <- fread(tsvPath)
  colnames(res)[colnames(res)=="chrom"] <- "CHROM"
  colnames(res)[colnames(res)=="pos"] <- "POS"
  colnames(res)[colnames(res)=="ref"] <- "REF"
  colnames(res)[colnames(res)=="alt"] <- "ALT"
  res[CHROM == "MT" ,"CHROM":="M"]
  #normalizeDeletions(res)
  return(res)
}

#TODO: try-catch
normalizeDeletions <- function(vcfData=""){
  res <- copy(vcfData)
  res[nchar(as.character(ALT))==1 & nchar(as.character(REF))>1 & substr(as.character(REF),1,1)==as.character(ALT),`:=`("POS"=POS+1L, "ALT"="-", "REF"=substring(as.character(REF),2))]
  return(res)
}

#TODO: try-catch
#gunzipPath, getFileExtension, loadVcfData, loadTsvData, normalizeDeletions
loadVarDb <- function(inputVariantsDbPath=""){
	inputVariantsDbPath <- gunzipPath(inputVariantsDbPath)
	varDbFormat <- getFileExtension(inputVariantsDbPath, removeLastExt = TRUE)
	if(varDbFormat=="vcf"){
		varDb <- loadVcfData(inputVariantsDbPath)
		if(verbose==TRUE){print("varDb VCF data loaded.")}
	}
	if(varDbFormat=="tsv"){
		varDb <- loadTsvData(inputVariantsDbPath) 
		if(verbose==TRUE){print("varDb TSV data loaded.")}
	}
	if(varDbFormat!="tsv" & varDbFormat!="vcf"){print("varDb file must be in tsv or vcf format!")}
	test_that("Is varDb is loaded as data.table?", {
	  expect_that(is.data.table(varDb), is_true())
	})
	test_that("Are CHROM,POS,REF,ALT present in varDb colnames?", {
	  expect_that(colnames(varDb[,c("CHROM","POS","REF","ALT")]), equals(c("CHROM","POS","REF","ALT")))
	})
	varDb <- normalizeDeletions(varDb)
	if(verbose==TRUE){print("Deletions normalized.")}
	#varDb <- setTenthColumnNameInVcfToGenotype(vcfTable = varDb)
	#test_that("varDb 10th colname", {
	#  expect_that(colnames(varDb)[10], equals(c("GEONTYPE")))
	#})

	##TG
	if(runsAtImid==TRUE){
		varDb <- varDb [which(varDb$pathogenic > 0 | varDb$likely_pathogenic > 0),]
	}
	return(varDb)
}

#TODO: try-catch
getUniqueLoci <- function(vcfData, variationType="", saveNucleotideData=FALSE){
  res <- vcfData[,c("CHROM","POS","REF","ALT")]
  #res[,"lengthDiff":=nchar(as.character(REF))-nchar(as.character(ALT))] # snp: 0, in: <0, del: >0
  if(variationType=="in"){res <- res[nchar(as.character(ALT))>1 & nchar(as.character(REF))==1 & as.character(REF)==substr(as.character(ALT),1,1)]}
  if(variationType=="del"){res <- res[(nchar(as.character(ALT))==1 & nchar(as.character(REF))>1 & substr(as.character(REF),1,1)==as.character(ALT)) | ALT=="-"]}
  if(variationType=="snp"){res <- res[nchar(as.character(ALT))==1 & nchar(as.character(REF))==1]}
  if(variationType=="notMatching"){res <- res[nchar(as.character(ALT))>1 & nchar(as.character(REF))>1]}
  if(saveNucleotideData==TRUE){
    res <- subset(unique(res[,c("CHROM","POS","REF","ALT")], by=c("CHROM","POS")))
  }else{
    res <- subset(unique(res[,c("CHROM","POS")], by=c("CHROM","POS")))
  }
  return(res)
}

#TODO: try-catch
#getUniqueLoci
getPositions <- function(varDb=NULL, new=FALSE){
  if(new==TRUE){
  	positionsAll <- getUniqueLoci(varDb)
  	positionsAll2 <- copy(positionsAll)
  	positionsAll2[,"CHROM":=paste0("chr",CHROM)]
  	fwrite(positionsAll2, outputPositionsAllPath, sep= " ")
  	rm(positionsAll2)
  	if(verbose==TRUE){print("positionsAll exported to external file.")}
  }else{
    if(file.exists(outputPositionsAllPath)){
      positionsAll <- fread(outputPositionsAllPath)
      positionsAll$CHROM <- gsub("chr", "", positionsAll$CHROM)
      if(verbose==TRUE){print("positionsAll restored from outputPositionsAllPath file.")}
    }else{
      stop("outputPositionsAllPath file is not present! Please provide a valid outputPositionsAllPath or set newPositionsAll=TRUE")
    }
  }
	return(positionsAll)
}

#TODO: try-catch
addRefFromReferenceGenome <- function(pileups){
  ir <- IRanges(as.matrix(pileups[,"POS"])[,1],width=1)
  gr <- GRanges(paste0("chr",c(as.matrix(pileups[,"CHROM"])[,1])), ir)
  if(hg38==TRUE){
    grseq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)
  }else{
    grseq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, gr)
  }
  grseq <- data.table(as.character(grseq))
  colnames(grseq)[colnames(grseq)=="V1"] <- "REF"
  res <- pileups[, "REF" := grseq[["REF"]]]
  return(res)
}

#TODO: try-catch
phredQualityDecode <- function(n){
  return(utf8ToInt(as.character(n))-33)
}

#TODO: try-catch
phredQualityEncode <- function(n){
  return(intToUtf8(as.integer(n)+33))
}

#TODO: try-catch, fix name, merge with doPileups?
#addRefFromReferenceGenome, phredQualityDecode, phredQualityEncode
pileupClinVar <- function(lociPath="",pileupPath="",bamPath="",minBaseQuality=0,minMapq=0,newPileups=TRUE,normalizeData=TRUE){
  if(newPileups==TRUE){
    cmd <- paste0("samtools mpileup -A -d 0 -l ",lociPath," -q ",minMapq," -Q ",minBaseQuality," -O -s ",bamPath," > ",pileupPath)
    system(cmd)
  }
  pileups <- fread(pileupPath, quote="")
  colnames(pileups) <- c("CHROM","POS","REF","DP","NUC","BQ","MQ","PR")
  pileups[,"REF":=NULL]
  pileups$CHROM <- gsub("chr", "", pileups$CHROM)
  pileups <- addRefFromReferenceGenome(pileups)
  #pileups$NUC <- gsub("$", "", pileups$NUC)
  if(normalizeData==TRUE){  
    pileups[,"InDelLengths":=lapply(str_extract_all(NUC,"[0-9]+"), sort, decreasing=TRUE)]
    pileups[,"InDelMaxLength":=as.numeric(lapply(InDelLengths,head, n=1))]
    pileups[is.na(InDelMaxLength)==TRUE,"InDelMaxLength":=0]
    inDelMaxLength <- pileups[,max(InDelMaxLength, na.rm = TRUE)]   
    pileups[,InDelLengths:=NULL]

    pileups$NUC <- str_replace_all(pileups$NUC,"\\^.","")
    for(n in 1:inDelMaxLength){
      pileups$NUC <- str_replace_all(pileups$NUC,paste0("\\-",n,"[ACGTNacgtn]{",n,"}"),"") # deletions: scan next position and look for "*" (remove "-" from this read)
      pileups$NUC <- str_replace_all(pileups$NUC,paste0("[ACGTNacgtn]{1}\\+",n,"[ACGTNacgtn]{",n,"}"),"+") # insertions: change base to "+"
    }
    pileups$NUC <- str_replace_all(pileups$NUC,"[^[ACGTNacgtn*-0123456789+]]","")
    pileups$NUC <- str_replace_all(pileups$NUC,"\\*","-")
    pileups$NUC <- toupper(pileups$NUC)
    #pileups[,"NCHAR":=(nchar(NUC)==DP)]
    
    # encode row
    #pileups[,"PR2":=phredQualityEncode(as.list(strsplit(PR,",")[[1]])), by=1:nrow(pileups)]
    # decode row
    #pileups[,"PR3":=paste(phredQualityDecode(PR2),sep=",",collapse = ","), by=1:nrow(pileups)]
    
    # decode BQ, MQ
    pileups[,"BQ":=paste(phredQualityDecode(BQ),sep=",",collapse = ","), by=1:nrow(pileups)]
    pileups[,"MQ":=paste(phredQualityDecode(MQ),sep=",",collapse = ","), by=1:nrow(pileups)]
    pileups[,`:=`(NUC=strsplit(NUC,""),BQ=strsplit(BQ,","),MQ=strsplit(MQ,","),PR=strsplit(PR,","))]
  }
  return(pileups)
}

#TODO: try-catch
calculateNucStats <- function(chrom="", pos="",ref="", nuc="", countsOnly=FALSE, bq="", mq="", pr=""){
  if(countsOnly==TRUE){
    dt <- data.table(
      CHROM=unlist(chrom)
      ,POS=unlist(pos)
      ,NUC=strsplit(unlist(nuc),"")[[1]]
    )
    
    nucCounts <- data.table(CHROM=dt[1]$CHROM,POS=dt[1]$POS,countA=nrow(dt[NUC=="A",]),countC=nrow(dt[NUC=="C",])
                            ,countG=nrow(dt[NUC=="G",]),countT=nrow(dt[NUC=="T",]),countIn=nrow(dt[NUC=="+",])
                            ,countDel=nrow(dt[NUC=="-",],countN=nrow(dt[NUC=="N",]))
    )
    nucCounts[is.na(nucCounts)] <- 0 
    return(nucCounts)
  }else{
    #x <- 5805
    dt <- data.table(
      CHROM=unlist(rep(chrom,length(nuc)))
      ,POS=unlist(rep(pos,length(nuc)))
      ,REF=unlist(rep(ref,length(nuc)))
      ,NUC=unlist(nuc)
      ,BQ=as.integer(unlist(bq)) # BQ decoded!
      ,MQ=as.integer(unlist(mq)) # MQ decoded!
      ,PR=as.integer(unlist(pr))
    )
    #View(dt)
    
    
    subsetRef <- subset(dt, NUC==REF)
    subsetNonRef <- subset(dt, NUC!=REF) 
    
    if(nrow(subsetRef)==0){
      minMQ_REF <- 0
      maxMQ_REF <- 0
      medianMQ_REF <- 0
      minBQ_REF <- 0
      maxBQ_REF <- 0
      medianBQ_REF <- 0
      minPR_REF <- 0
      maxPR_REF <- 0
      medianPR_REF <- 0
    }else{
      minMQ_REF <- as.double(min(subsetRef[,MQ]))
      maxMQ_REF <- as.double(max(subsetRef[,MQ]))
      medianMQ_REF <- as.double(median(subsetRef[,MQ]))
      minBQ_REF <- as.double(min(subsetRef[,BQ]))
      maxBQ_REF <- as.double(max(subsetRef[,BQ]))
      medianBQ_REF <- as.double(median(subsetRef[,BQ]))
      minPR_REF <- as.double(min(subsetRef[,PR]))
      maxPR_REF <- as.double(max(subsetRef[,PR]))
      medianPR_REF <- as.double(median(subsetRef[,PR]))
    }
    
    
    if(nrow(subsetNonRef)==0){
      minMQ_NREF <- 0
      maxMQ_NREF <- 0
      medianMQ_NREF <- 0
      minBQ_NREF <- 0
      maxBQ_NREF <- 0
      medianBQ_NREF <- 0
      minPR_NREF <- 0
      maxPR_NREF <- 0
      medianPR_NREF <- 0
      
      maxNonRefCount <- 0
      minMQ_MNREF <- 0
      maxMQ_MNREF <- 0
      medianMQ_MNREF <- 0
      minBQ_MNREF <- 0
      maxBQ_MNREF <- 0
      medianBQ_MNREF <- 0
      minPR_MNREF <- 0
      maxPR_MNREF <- 0
      medianPR_MNREF <- 0
    }else{
      minMQ_NREF <- as.double(min(subsetNonRef[,MQ]))
      maxMQ_NREF <- as.double(max(subsetNonRef[,MQ]))
      medianMQ_NREF <- as.double(median(subsetNonRef[,MQ]))
      minBQ_NREF <- as.double(min(subsetNonRef[,BQ]))
      maxBQ_NREF <- as.double(max(subsetNonRef[,BQ]))
      medianBQ_NREF <- as.double(median(subsetNonRef[,BQ]))
      minPR_NREF <- as.double(min(subsetNonRef[,PR]))
      maxPR_NREF <- as.double(max(subsetNonRef[,PR]))
      medianPR_NREF <- as.double(median(subsetNonRef[,PR]))
      
      
      maxNonRefChrom <- names(sort(table(subsetNonRef$NUC),decreasing=TRUE)[1])
      subsetMaxNonRef <- subset(subsetNonRef,NUC==maxNonRefChrom)
      
      maxNonRefCount <- nrow(subsetMaxNonRef)#max(subsetNonRef[, .(count = .N), by = NUC]$count)
      minMQ_MNREF <- as.double(min(subsetMaxNonRef[,MQ]))
      maxMQ_MNREF <- as.double(max(subsetMaxNonRef[,MQ]))
      medianMQ_MNREF <- as.double(median(subsetMaxNonRef[,MQ]))
      minBQ_MNREF <- as.double(min(subsetMaxNonRef[,BQ]))
      maxBQ_MNREF <- as.double(max(subsetMaxNonRef[,BQ]))
      medianBQ_MNREF <- as.double(median(subsetMaxNonRef[,BQ]))
      minPR_MNREF <- as.double(min(subsetMaxNonRef[,PR]))
      maxPR_MNREF <- as.double(max(subsetMaxNonRef[,PR]))
      medianPR_MNREF <- as.double(median(subsetMaxNonRef[,PR]))
    }
    
    summaryDt <- data.table(#CHROM=dt[1]$CHROM,POS=dt[1]$POS,#REF=dt[1]$REF,
      #NUC=nuc,BQ=bq,MQ=mq,PR=pr,
      countA=nrow(dt[NUC=="A",]),countC=nrow(dt[NUC=="C",]),countG=nrow(dt[NUC=="G",]),
      countT=nrow(dt[NUC=="T",]),countIn=nrow(dt[NUC=="+",]),countDel=nrow(dt[NUC=="-",]),countN=nrow(dt[NUC=="N",]),
      RefCount=as.integer(nrow(subsetRef)),#dt[NUC==REF, .(count = .N), by = NUC], #policz w NUC wystąpienia REF
      AllNonRefCount=as.integer(nrow(subsetNonRef)),#dt[NUC!=REF, .(count = .N), by = NUC],
      MaxNonRefCount=as.integer(maxNonRefCount),
      
      #DP2=nrow(dt),
      #MaxNonRefRate=maxNonRefCount/nrow(dt[NUC==REF]),
      #VAF=MaxNonRefCount/DP2,
      #VRF=RefCount/DP2,
      
      MinMQ=as.double(min(dt[,MQ])),#all ref maxnonref allnonref
      MaxMQ=as.double(max(dt[,MQ])),
      MedianMQ=as.double(median(dt[,MQ])),
      MinBQ=as.double(min(dt[,BQ])),
      MaxBQ=as.double(max(dt[,BQ])),
      MedianBQ=as.double(median(dt[,BQ])),
      MinPR=as.double(min(dt[,PR])),
      MaxPR=as.double(max(dt[,PR])),
      MedianPR=as.double(median(dt[,PR])),
      
      
      MinMQ_REF=minMQ_REF,#all ref maxnonref allnonref
      MaxMQ_REF=maxMQ_REF,
      MedianMQ_REF=medianMQ_REF,
      MinBQ_REF=minBQ_REF,
      MaxBQ_REF=maxBQ_REF,
      MedianBQ_REF=medianBQ_REF,
      MinPR_REF=minPR_REF,
      MaxPR_REF=maxPR_REF,
      MedianPR_REF=medianPR_REF,
      
      
      MinMQ_NREF=minMQ_NREF,#all ref maxnonref allnonref
      MaxMQ_NREF=maxMQ_NREF,
      MedianMQ_NREF=medianMQ_NREF,
      MinBQ_NREF=minBQ_NREF,
      MaxBQ_NREF=maxBQ_NREF,
      MedianBQ_NREF=medianBQ_NREF,
      MinPR_NREF=minPR_NREF,
      MaxPR_NREF=maxPR_NREF,
      MedianPR_NREF=medianPR_NREF,
      
      MinMQ_MNREF=minMQ_MNREF,#all ref maxnonref allnonref
      MaxMQ_MNREF=maxMQ_MNREF,
      MedianMQ_MNREF=medianMQ_MNREF,
      MinBQ_MNREF=minBQ_MNREF,
      MaxBQ_MNREF=maxBQ_MNREF,
      MedianBQ_MNREF=medianBQ_MNREF,
      MinPR_MNREF=minPR_MNREF,
      MaxPR_MNREF=maxPR_MNREF,
      MedianPR_MNREF=medianPR_MNREF
      #class
    )
    
    summaryDt[is.na(summaryDt)] <- 0 
    #summaryDt[is.infinite(summaryDt)] <- 0 
    #for (j in 1:ncol(summaryDt)){set(summaryDt, which(!is.finite(summaryDt[[j]])), j, 0)}
    
    return(summaryDt)
  }
}

#TODO: try-catch
#calculateNucStats
addStatsColumnsToPileup <- function(pileups=NULL, outputPileupsWithStatsPath="", countsOnly=FALSE, rmCols=FALSE, newPileupsWithStats=TRUE){
  #TODO: READ pileups with stats
  if(newPileupsWithStats==TRUE){
    res <- copy(pileups)
    tic(paste0("Stats calculated. Pileups, nrow: ", nrow(res), " time"))
    if(countsOnly==TRUE){
      res[,c("CHROM","POS","countA","countC","countG","countT","countIn","countN","countDel"):=
            calculateNucStats(chrom=CHROM, pos=POS, nuc=NUC, countsOnly=TRUE, bq=BQ, mq=MQ, pr=PR), by=1:nrow(res)]
    }else{
      #View(res)
      res[,c(#"CHROM","POS","REF",
        #"NUC","BQ","MQ","PR",
        "countA","countC","countG","countT","countIn","countDel","countN",
        "RefCount",
        "AllNonRefCount",
        "MaxNonRefCount",
        
        #"DP2",
        #"MaxNonRefRate",
        
        "MinMQ",
        "MaxMQ",
        "MedianMQ",
        "MinBQ",
        "MaxBQ",
        "MedianBQ",
        "MinPR",
        "MaxPR",
        "MedianPR",
        "MinMQ_REF",
        "MaxMQ_REF",
        "MedianMQ_REF",
        "MinBQ_REF",
        "MaxBQ_REF",
        "MedianBQ_REF",
        "MinPR_REF",
        "MaxPR_REF",
        "MedianPR_REF",
        "MinMQ_NREF",
        "MaxMQ_NREF",
        "MedianMQ_NREF",
        "MinBQ_NREF",
        "MaxBQ_NREF",
        "MedianBQ_NREF",
        "MinPR_NREF",
        "MaxPR_NREF",
        "MedianPR_NREF",
        "MinMQ_MNREF",
        "MaxMQ_MNREF",
        "MedianMQ_MNREF",
        "MinBQ_MNREF",
        "MaxBQ_MNREF",
        "MedianBQ_MNREF",
        "MinPR_MNREF",
        "MaxPR_MNREF",
        "MedianPR_MNREF"
      ):=calculateNucStats(chrom=CHROM, pos=POS, ref=REF, nuc=NUC, 
                           countsOnly=FALSE, bq=BQ, mq=MQ, pr=PR), by=1:nrow(res)]
    }  
    toc()
    if(rmCols==TRUE){
      res[, `:=`(NUC = NULL,BQ = NULL,MQ = NULL,PR = NULL)]
    }
    fwrite(res, outputPileupsWithStatsPath, sep= " ")
  }else{
    if(file.exists(outputPileupsWithStatsPath)){
      res <- fread(outputPileupsWithStatsPath)
      if(verbose==TRUE){print("pileupsWithStatsPath restored from outputPileupsWithStatsPath file.")}
    }else{
      stop("outputPileupsWithStatsPath file is not present! Please provide a valid outputPileupsWithStatsPath or set newPileupsWithStatsPath=TRUE")
    }
    
  }
  
  return(res)
}

#TODO: try-catch
#pileupClinVar, calculateNucStats, addStatsColumnsToPileup
doPileups <- function(outputPositionsPath="",outputPileupsPath="",outputPileupsWithStatsPath="",bamPath="",minBaseQuality=0,minMapq=0,newPileups=TRUE,newPileupsWithStats=TRUE){
	if(newPileups==TRUE & newPileupsWithStats!=TRUE){
	  pileups <- pileupClinVar(outputPositionsPath,outputPileupsPath,bamPath,minBaseQuality,minMapq,newPileups=TRUE)
	  if(verbose==TRUE){print("pileups done.")}
	}else{
	  if(file.exists(outputPileupsPath)){
		pileups <- pileupClinVar(outputPositionsPath,outputPileupsPath,bamPath,minBaseQuality,minMapq,newPileups=FALSE)
		if(verbose==TRUE){print("pileups restored from outputPileupsPath file.")}
	  }else{
		stop("outputPileupsPath file is not present! Please provide a valid outputPileupsPath or set newPileups=TRUE")
	  }
	}
  #pileups <- pileupsBackup[1:100]
  #colnames(pileups)[colnames(pileups)=="posRead"] <- "PR"
  #pileups[,`:=`(NUC=strsplit(NUC,""),BQ=strsplit(BQ,","),MQ=strsplit(MQ,","),PR=strsplit(PR,","))]
	pileups <- addStatsColumnsToPileup(pileups=pileups, outputPileupsWithStatsPath=outputPileupsWithStatsPath, countsOnly = !(fullPileupStats), rmCols = FALSE, newPileupsWithStats=newPileupsWithStats) #TODO: make it true in release
	if(verbose==TRUE){print("Pileups: New columns created.")}
	#pileupsTests(pileups)
	return(pileups)
}


#======================================================train

#TODO: try-catch, make CC entry argument from this
initializeHistogramData <- function(lowerBounds = c(1,2,6,11,21,51,101,201), upperBounds = c(1,5,10,20,50,100,200,9999)){
	if(length(lowerBounds)!=length(upperBounds)){
		stop("lowerBounds and upperBounds have different sizes!")
    }
	if(!all(lowerBounds<=upperBounds)){
		stop("All values in lowerBounds have to be lower than corresponding values in upperBounds!")
    }
	histogramData <- data.table("lower"=lowerBounds,
                           "upper"=upperBounds,
                           "DP"=c(matrix(1,1,length(lowerBounds))),
                           "DP_rand"=c(matrix(1,1,length(lowerBounds)))
						   )
	return(histogramData)
}


#TODO: try-catch, merge with something else?
getDP <- function(pileupsWithNucCols=NULL, removeNucColumns=FALSE){
  res <- copy(pileupsWithNucCols)
  res[,"DP":=countA+countC+countG+countT+countIn+countDel+countN]
  if(removeNucColumns==TRUE){
    res[, `:=`(countA = NULL, countC = NULL, countG = NULL, countT = NULL, countIn = NULL, countDel = NULL, countN= NULL)]
  }
  return(res)
}

#TODO: try-catch
#getDP
getSubsetDepths <- function(subsetBounds=NULL, pileupsWithNucCols=NULL){
  res <- copy(subsetBounds)
  totalDepths <- getDP(pileupsWithNucCols=pileupsWithNucCols, removeNucColumns=TRUE)
  res <- sapply(1:nrow(subsetBounds), function(i){
    return(nrow(subset(totalDepths, (DP >= subsetBounds[i,lower] & DP <= subsetBounds[i,upper]))))
  })
  return(res)
}

#TODO: try-catch
#initializeHistogramData, getSubsetDepths, getDP
updateDepthsColumn <- function(subsetBounds=NULL, updateDp=TRUE, updateDpRand=TRUE, pileupsWithNucCols=NULL, falsePileupsWithNucCols=NULL){
  if(!exists("subsetBounds")){
    subsetBounds <- initializeHistogramData()
  }
  res <- copy(subsetBounds)
  if(updateDp==TRUE){
    subsetDepths <- getSubsetDepths(subsetBounds=subsetBounds, pileupsWithNucCols=pileupsWithNucCols)
    res[,"DP":=subsetDepths]
    if(verbose==TRUE){print("subsetBounds DP updated.")}
  }
  if(updateDpRand==TRUE){
    subsetDepthsF <- getSubsetDepths(subsetBounds=subsetBounds, pileupsWithNucCols=falsePileupsWithNucCols)
    res[,"DP_rand":=subsetDepthsF]
    if(verbose==TRUE){print("subsetBounds DP_rand updated.")}
  }
  return(res)
}

#TODO: try-catch
getCoverageSubsetByReadsCount <- function(coverage, min=1, max=1){
  if(min<=max){
    return(subset(coverage, (readsCount >= min & readsCount <= max)))
  }else{
    print("max must be greater or equal to min!")
    return(NULL)
  }  
}

#TODO: try-catch
getPositionsFromCoverageSubset <- function(coverageSubset,subsetMaxLength=200000){
  if(nrow(coverageSubset)>subsetMaxLength){coverageSubset <- coverageSubset[sample(.N, subsetMaxLength, replace=F)]}
  #tic()
  #res2 <- rbindlist(lapply(1:nrow(f), function(i){data.table("CHROM"=f[i,CHROM],"POS"=f[i,start]:f[i,stop],"readsCount"=f[i,readsCount])}))
  #toc()
  
  #tic()
  res <- rbindlist(
      lapply(c(1:22, "X", "Y", "M"), function(i){
      chromSubset <- subset(coverageSubset, CHROM==i)
      #print(paste0("chr",i))
      if(nrow(chromSubset)){
        chromPositions <- rbindlist(
          lapply(min(chromSubset[,readsCount]):max(chromSubset[,readsCount]), function(j){
          countSubset <- subset(chromSubset, readsCount==j)
          #print(paste0("readsCount",j))
          if(nrow(countSubset)){
            countPositions <- rbindlist(
              lapply(1:nrow(countSubset), function(k){
              data.table("POS"=countSubset[k,start]:countSubset[k,stop])
              }))
            countPositions[,"readsCount":=j]
            return(countPositions)
          }  
        }))
        chromPositions[,"CHROM":=i] 
        return(chromPositions)
      }
    }))
  #toc()  
  #setequal(res[order("CHROM","POS","readsCount")],res2[order("CHROM","POS","readsCount")])
  return(res)
}

#TODO: try-catch
getRandomPositions <- function(listOfPositions, randomSubsetSize=10000){
  if(nrow(listOfPositions)>=randomSubsetSize){
    res <- listOfPositions[sample(.N, randomSubsetSize, replace=F)]
  }else{
    res <- listOfPositions  
  }
  return(res)
}

#TODO: try-catch
#getCoverageSubsetByReadsCount, getPositionsFromCoverageSubset, getRandomPositions
getRandomPositionsFromCoverage <- function(coverage,subsetBounds,subsetMaxLength=200000){
  res <- lapply(1:nrow(subsetBounds), function(i){
    coverageSubset <- getCoverageSubsetByReadsCount(coverage,min=subsetBounds[i,lower],max=subsetBounds[i,upper])
    listOfPositions <- getPositionsFromCoverageSubset(coverageSubset,subsetMaxLength)
    randomPositions <- getRandomPositions(listOfPositions, subsetBounds[i,DP])
    return(randomPositions)
  })
  res <- rbindlist(res)
  return(res) #res size: nrow(subsetBounds) * randomSubsetSize
}

#TODO: try-catch, clear comments, fix name
#getRandomPositionsFromCoverage, getCoverageSubsetByReadsCount, getPositionsFromCoverageSubset, getRandomPositions
getRandomPositionsExcludeSubset2 <- function(coverage,variantsSetToExclude,subsetBounds,returnReadsCount=FALSE){ 
  #count - samples per chromosome
  #variantsSetToExclude - set of variants to exclude from selection
  posCount <- sum(subsetBounds[,"DP"])
  #print(posCount)
  #===draw, check, draw again until success===
  #calculate initial randomSubsetSize
  randomSubsetSize <- ceiling(posCount/nrow(subsetBounds))
  subsetMaxLength <- randomSubsetSize
  #print(paste0(randomSubsetSize,", ",subsetMaxLength))
  #set initial tmpCount == count
  tmpPosCount <- posCount
  #get first set of random positions
  res <- getRandomPositionsFromCoverage(coverage,subsetBounds,subsetMaxLength) #res -> nrow(subsetBounds) * randomSubsetSize
  repeat{
    #remove duplicates
    res <- subset(unique(res, by=c("CHROM","POS")))
    # remove loci present in variantsSetToExclude
    #print("=======================================================================================================")
    #dts_disjoint(res[,c("CHROM","POS")], variantsSetToExclude[,c("CHROM","POS")])
    variantsSetToExclude[,"CHROM":=as.character(CHROM)]
    res <- res[!variantsSetToExclude, on=.(CHROM,POS)]
    #dts_disjoint(res[,c("CHROM","POS")], variantsSetToExclude[,c("CHROM","POS")])
    #end if number of drawn loci is equal or greater than count
    if(nrow(res)>=posCount){
      break
    }
    #update tmpCount
    tmpPosCount <- posCount - nrow(res)
    #print(tmpPosCount)
    #update randomSubsetSize
    randomSubsetSize <- ceiling(tmpPosCount/nrow(subsetBounds))
    subsetMaxLength <- randomSubsetSize*3
    #draw next set of random positions
    res <- rbind(res,getRandomPositionsFromCoverage(coverage,subsetBounds,subsetMaxLength))
  }
  #ensure that the resulting loci set will not be larger than specified in count
  if(returnReadsCount==TRUE){
    res <- res[sample(.N, posCount, replace=F),c("CHROM","POS","readsCount")]
  }else{
    res <- res[sample(.N, posCount, replace=F),c("CHROM","POS")]
  }
  return(res)
}

dts_disjoint <- function(dt1,dt2){
  if(nrow(fsetdiff(dt1,dt2)) == nrow(dt1)){
    print('true')
  }else{
    print('false')
  }
}

#TODO: try-catch
#getRandomPositionsExcludeSubset2, getRandomPositionsFromCoverage, getCoverageSubsetByReadsCount, getPositionsFromCoverageSubset, getRandomPositions, dts_disjoint
getFalsePositions <- function(outputPositionsAllFalsePath="", new=TRUE, coverage=NULL, variantsSetToExclude=NULL, histogramData=NULL){
	if(new==TRUE){
	  positionsAllFalse <- getRandomPositionsExcludeSubset2(coverage, variantsSetToExclude=variantsSetToExclude, histogramData,returnReadsCount=FALSE)
	  if(verbose==TRUE){print("Random positions generated.")}
	  
	  positionsAllFalse2 <- copy(positionsAllFalse)
	  positionsAllFalse2[,"CHROM":=paste0("chr",CHROM)]
	  fwrite(positionsAllFalse2, outputPositionsAllFalsePath, sep= " ")
	  rm(positionsAllFalse2)
	  if(verbose==TRUE){print("positionsAllFalse exported to external file.")}
	  
	  test_that('Are positionsAllFalse and positionsAll disjoint?', {
		expect_that(dts_disjoint(positionsAllFalse[,c("CHROM","POS")],positionsAll[,c("CHROM","POS")]), prints_text('true'))
	  })
	  #falsePileups <- pileupClinVar(positionsAllFalse, 25000) # loci w/o variants - FPR
	}else{
	  if(file.exists(outputPositionsAllFalsePath)){
		positionsAllFalse <- fread(outputPositionsAllFalsePath)
		positionsAllFalse$CHROM <- gsub("chr", "", positionsAllFalse$CHROM)
		if(verbose==TRUE){print("positionsAllFalse restored from outputPositionsAllFalsePath file.")}
	  }else{
		stop("outputPositionsAllFalsePath file is not present! Please provide a valid outputPositionsAllFalsePath or set newPositionsAllFalse=TRUE")
	  }
	}
	return(positionsAllFalse)
}

#TODO: try-catch
#addRefFromReferenceGenome
prepareTrainingData <- function(pileups=pileups, falsePileups=falsePileups){
  pileups2 <- copy(pileups)
  pileupsRandom2 <- copy(falsePileups)
  pileups2[,"VAR":="variant"]
  pileupsRandom2[,"VAR":="nonVariant"]
  pileupsAll <- rbind(pileups2,pileupsRandom2)
  pileupsAll <- pileupsAll[,"VAR":=as.factor(VAR)]
  pileupsAll <- addRefFromReferenceGenome(pileupsAll)
  pileupsAllShuffled <- pileupsAll[sample(nrow(pileupsAll)),]
  rm(pileups2,pileupsRandom2,pileupsAll)
  return(pileupsAllShuffled)
}

#TODO: try-catch
trainModel <- function(modelPath="", data=NULL, listOfTrainArguments=NULL, mc=0){
  registerDoMC(cores=mc)
  seed <- 7
  
  set.seed(seed)
  model <- train(listOfTrainArguments$formula, data=data, method=listOfTrainArguments$method, 
                 metric=listOfTrainArguments$metric, preProc=listOfTrainArguments$preProc, 
                 trControl=listOfTrainArguments$control)
  
  saveRDS(model, modelPath)
  registerDoMC(cores=0)
  return(model)
}


classifyNaive <- function(newData=NULL){
  out <- tryCatch(
    {
      data <- copy(newData)
      data[AllNonRefCount>0,result:="variant"]
      data[AllNonRefCount==0,result:="nonVariant"]
      data[,result:=as.factor(result)]
    },
    error=function(cond) {
      message("Something went wrong with classification")
      message("Here's the original error message:")
      message(cond)
      return(NA)
    },
    finally={
      return(data)
    }
  )  
  if(verbose==TRUE){print("Naive variant classification complete.")}
  return(out)
}


classify <- function(model=NULL, newData=NULL){
  out <- tryCatch(
    {
      if(typeof(model)=="character"){
        if(file.exists(model)){
          externalModel <- readRDS(model)
          result <- predict(externalModel, newdata = newData)
        }else{
          stop("modelPath file is not present! Please provide a valid classification model!")
        }
      }else{
        result <- predict(model, newdata = newData)
      }
    },
    error=function(cond) {
      message("Something went wrong with classification. Original error message:")
      message(cond)
      return(NA)
    },
    finally={
      result <- cbind(result,newData)
      return(result)
    }
  )  
  if(verbose==TRUE){print("Variant classification complete.")}
  return(out)
}

#TODO: try-catch, clear comments
mergePileupsWithVarDb <- function(pileups=NULL, varDb=NULL, allx=TRUE, ally=TRUE){
  tmp <- copy(pileups)
  if("REF" %in% colnames(tmp)){tmp[, `:=`(REF = NULL)]}
  if("NUC" %in% colnames(tmp)){tmp[, `:=`(NUC = NULL)]}
  if("BQ" %in% colnames(tmp)){tmp[, `:=`(BQ = NULL)]}
  if("MQ" %in% colnames(tmp)){tmp[, `:=`(MQ = NULL)]}
  if("PR" %in% colnames(tmp)){tmp[, `:=`(PR = NULL)]}
  #print(colnames(pileups))
  #print(colnames(varDb))

  #varDb[,"CHROMtmp":=CHROM]
  #varDb[,"CHROM":=NULL]
  #varDb[,"CHROM":=as.character(CHROMtmp)]
  #varDb[,"CHROMtmp":=NULL]
  varDb[,"CHROM":=as.character(CHROM)]
  #View(varDb)
  res <- merge(varDb, tmp[,c("CHROM","POS","RefCount","AllNonRefCount","MaxNonRefCount","InDelMaxLength",
                             "countA","countC","countG","countT","countIn","countDel","countN","result")], by=c("CHROM","POS"), all.x=allx, all.y=ally)
  res$result[is.na(res$result)] <- "nonVariant"  #WARN: potentially breaks confusion matrix
  rm(tmp)
  #TODO: maybe shorter? 
  res[,][is.na(res[,])] <- 0 
  
  #res[,c("countA","countC","countG","countT","countIn","countDel","countN")][is.na(res[,c("countA","countC","countG","countT","countIn","countDel","countN")])] <- 0
  #res[,c("countA","countC","countG","countT","countIn","countDel","countN")][is.na(res[,c("countA","countC","countG","countT","countIn","countDel","countN")])] <- 0 
  return(res)
}

#TODO: try-catch, clear comments
calculateDepths <- function(variantsTable=""){
  
  # nucleotides matching REF
  #if(!("RD" %in% colnames(variantsTable))){
  #  variantsTable[, `:=`(RD = integer(0))]
  #  }
  if("RD" %in% colnames(variantsTable)){variantsTable[, `:=`(RD = NULL)]}
  variantsTable[, `:=`(RD = integer(0))]
  variantsTable[,c("RD")][is.na(variantsTable[,c("RD")])] <- 0
  variantsTable[substr(REF,1,1) == "A" ,"RD":=RD+countA]
  variantsTable[substr(REF,1,1) == "C" ,"RD":=RD+countC]
  variantsTable[substr(REF,1,1) == "G" ,"RD":=RD+countG]
  variantsTable[substr(REF,1,1) == "T" ,"RD":=RD+countT]
  ## variantsTable[nchar(as.character(alt))<nchar(as.character(ref)) & nchar(as.character(alt))>1, "RD":=RD+countIn]
  ## variantsTable[nchar(as.character(alt))>nchar(as.character(ref)) & nchar(as.character(ref))>1, "RD":=RD+countDel]
  #variantsTable[substr(REF,1,1) == "N" ,"RD":=RD+countN] # REF is never N 
  
  # nucleotides matching ALT
  #if(!("AD" %in% colnames(variantsTable))){
  #  variantsTable[, `:=`(AD = integer(0))]
  #  }
  if("AD" %in% colnames(variantsTable)){variantsTable[, `:=`(AD = NULL)]}
  variantsTable[, `:=`(AD = integer(0))]
  variantsTable[,c("AD")][is.na(variantsTable[,c("AD")])] <- 0
  variantsTable[ALT == "A" ,"AD":=AD+countA]
  variantsTable[ALT == "C" ,"AD":=AD+countC]
  variantsTable[ALT == "G" ,"AD":=AD+countG]
  variantsTable[ALT == "T" ,"AD":=AD+countT]
  #variantsTable[nchar(as.character(ALT))>nchar(as.character(REF)) & substr(as.character(REF),1,1)==substr(as.character(ALT),1,1),"AD":=AD+countIn]
  #variantsTable[(nchar(as.character(ALT))<nchar(as.character(REF)) & substr(as.character(REF),1,1)==substr(as.character(ALT),1,1)) | ALT:="-","AD":=AD+countDel]
  #variantsTable[substr(ALT,1,1) == "N" ,"AD":=AD+countN] # ALT is never N 
  
  # nucleotides not matching REF or ALT
  #if(!("nullDepth" %in% colnames(variantsTable))){
  #  variantsTable[, `:=`(nullDepth = integer(0))]
  #  }
  if("nullDepth" %in% colnames(variantsTable)){variantsTable[, `:=`(nullDepth = NULL)]}
  variantsTable[, `:=`(nullDepth = integer(0))]
  variantsTable[,c("nullDepth")][is.na(variantsTable[,c("nullDepth")])] <- 0
  variantsTable[substr(REF,1,1) != "A" & ALT != "A" ,"nullDepth":=nullDepth+countA]
  variantsTable[substr(REF,1,1) != "C" & ALT != "C" ,"nullDepth":=nullDepth+countC]
  variantsTable[substr(REF,1,1) != "G" & ALT != "G" ,"nullDepth":=nullDepth+countG]
  variantsTable[substr(REF,1,1) != "T" & ALT != "T" ,"nullDepth":=nullDepth+countT]
  #variantsTable[nchar(as.character(ALT))<=nchar(as.character(REF)),"nullDepth":=nullDepth+countIn]
  #variantsTable[nchar(as.character(ALT))>=nchar(as.character(REF)),"nullDepth":=nullDepth+countDel]
  variantsTable[,"nullDepth":=nullDepth+countN] # REF and ALT are never N 
  
  # coverage depth
  variantsTable[,"DP":=countA+countC+countG+countT+countN]
  #variantsTable[,"depthCtrl":=AD+RD+nullDepth]
  
  # error flag
  #variantsTable[DP!=depthCtrl, "errFlag":=1]
  variantsTable[nchar(as.character(REF))>1 & nchar(as.character(ALT))>1, "errFlag":=1]
  variantsTable[,c("errFlag")][is.na(variantsTable[,c("errFlag")])] <- 0 
  
  return(variantsTable)
}



#TODO: try-catch, merge with getVcfHcFc
#gunzipPath, loadVcfData, filterChroms, normalizeDeletions
getVcfHcGeneral <- function(vcfHcGeneralPath=""){
  # GATK HAPLOTYPECALLER - general scan
  test_that("Does vcfHcGeneral file exist?", {
    expect_that(file.exists(vcfHcGeneralPath), is_true())
  })
  
  vcfHcGeneralPath <- gunzipPath(vcfHcGeneralPath)
  
  vcfHcGeneral <- loadVcfData(vcfHcGeneralPath)
  if(verbose==TRUE){print("vcfHcGeneral file loaded.")}
  vcfHcGeneral <- filterChroms(vcfHcGeneral)
  vcfHcGeneral <- normalizeDeletions(vcfHcGeneral)
  
  test_that("Are only valid chromosomes (1:22, X, Y, M) present in vcfHcGeneral?", {
    expect_that(unique(unique(vcfHcGeneral$CHROM) %in% c(1:22, "X", "Y", "M")), is_true())
  })
  return(vcfHcGeneral)
}

#TODO: try-catch
setTenthColumnNameInVcfToGenotype <- function(vcfTable=""){
  colnames(vcfTable)[10] <- "GENOTYPE"
  return(vcfTable)
}

#TODO: try-catch, merge with getVcfHcGeneral
#setTenthColumnNameInVcfToGenotype
getVcfHcFc <- function(vcfHcFcPath=""){
  # GATK HAPLOTYPECALLER - forcecalling
  test_that("Does vcfHcFc file exist?", {
    expect_that(file.exists(vcfHcFcPath), is_true())
  })
  
  vcfHcFcPath <- gunzipPath(vcfHcFcPath)

  vcfHcFc <- loadVcfData(vcfHcFcPath)
  if(verbose==TRUE){print("vcfHcFc file loaded.")}
  vcfHcFc <- setTenthColumnNameInVcfToGenotype(vcfHcFc)
  vcfHcFc01 <- vcfHcFc[like(GENOTYPE,"0/1")]
  vcfHcFc11 <- vcfHcFc[like(GENOTYPE,"1/1")]
  vcfHcFc <- rbind(vcfHcFc01, vcfHcFc11)
  rm(vcfHcFc01, vcfHcFc11)
  
  vcfHcFc <- filterChroms(vcfHcFc)
  vcfHcFc <- normalizeDeletions(vcfHcFc)
  
  
  test_that("Are only 0/1 and 1/1 present in GENOTYPE?", {
    expect_that(unique(substr(vcfHcFc$GENOTYPE,1,3)), equals(c("0/1","1/1")))
  })
  return(vcfHcFc)
}

#TODO: try-catch
getClassificationStats <- function(subsetBounds=NULL, classificationData=NULL){
  cd <- classificationData
  bounds <- copy(subsetBounds)
  rowAll <- bounds[1]
  rowAll[,c("lower","upper"):=data.table(min(bounds$lower),max(bounds$upper))]
  bounds <- rbind(rowAll, bounds)
  res <- lapply(1:nrow(bounds), function(i){
    lower <- bounds[i,lower]
    upper <- bounds[i,upper]
    params <- list(
      lower = lower,
      upper = upper,
      TP = length(which(cd$result == "variant" & cd$VAR == "variant" & cd$DP >= lower & cd$DP <= upper)),
      FN = length(which(cd$result == "nonVariant" & cd$VAR == "variant" & cd$DP >= lower & cd$DP <= upper)),
      TN = length(which(cd$result == "nonVariant" & cd$VAR == "nonVariant" & cd$DP >= lower & cd$DP <= upper)),
      FP = length(which(cd$result == "variant" & cd$VAR == "nonVariant" & cd$DP >= lower & cd$DP <= upper)),
      P = length(which(cd$VAR == "variant" & cd$DP >= lower & cd$DP <= upper)),
      N = length(which(cd$VAR == "nonVariant" & cd$DP >= lower & cd$DP <= upper)),
      TPR = length(which(cd$result == "variant" & cd$VAR == "variant" & cd$DP >= lower & cd$DP <= upper)) / length(which(cd$VAR == "variant" & cd$DP >= lower & cd$DP <= upper)),
      FPR = length(which(cd$result == "variant" & cd$VAR == "nonVariant" & cd$DP >= lower & cd$DP <= upper)) / length(which(cd$VAR == "nonVariant" & cd$DP >= lower & cd$DP <= upper)),
      samplesTotal = length(which(cd$DP >= lower & cd$DP <= upper))
    )
    return(params)
  })
  res <- rbindlist(res)
  return(res) 
}

#TODO: try-catch
evaluateClassifier <- function(naiveDataStats=NULL, classifierDataStats=NULL){
  if(nrow(naiveDataStats)!=nrow(classifierDataStats)){
    stop("naiveDataStats and classifierDataStats have different sizes!")
  }
  if(!setequal(naiveDataStats$lower,classifierDataStats$lower) | !setequal(naiveDataStats$upper,classifierDataStats$upper)){
    stop("naiveDataStats and classifierDataStats do not match!")
  }
  res <- lapply(1:nrow(naiveDataStats), function(i){
    if(naiveDataStats[[i,"FPR"]]!=0){
      FPRdiffProc = round(100*((naiveDataStats[[i,"FPR"]] - classifierDataStats[[i,"FPR"]])/naiveDataStats[[i,"FPR"]]),digits = 4)
    }else{
      FPRdiffProc = Inf
    }
    if(naiveDataStats[[i,"TPR"]]!=0){
      TPRdiffProc = -round(100*((naiveDataStats[[i,"TPR"]] - classifierDataStats[[i,"TPR"]])/naiveDataStats[[i,"TPR"]]),digits = 4)
    }else{
      TPRdiffProc = -Inf
    }
    params <- list(
      FPRdiffProc = FPRdiffProc,
      TPRdiffProc = TPRdiffProc
    )
    return(params)
  })
  res <- rbindlist(res)
  res <- cbind(naiveDataStats[,c("lower","upper")],res)
  return(res)
}

#TODO: try-catch, input: plot title etc.
evaluateClassifierPlot <- function(evaluationData=NULL){
  eval <- copy(evaluationData)
  #lower-upper -> make label
  labelHalfLength <- nchar(as.character(max(eval$upper)))
  eval[,subset:=paste(str_pad(lower,labelHalfLength,pad = "0"),str_pad(upper,labelHalfLength,pad = "0"),sep = "-")]
  eval[, `:=`(lower=NULL,upper=NULL)]
  
  eval <- melt(eval, id.vars="subset", measure.vars = c("FPRdiffProc","TPRdiffProc"), variable.name = "metric")
  ggplot(eval, aes(subset, value, fill=metric)) + 
    geom_bar(stat="identity", position="dodge") + 
    geom_text(aes(label=round(value,digits = 2)), vjust=-0.3, position=position_dodge(1.0), size=2.5) +
    scale_fill_manual(values = c("#0000FF","#FF0000")) +
    ggtitle("lda")
}







#main part of the script
#TODO: put HC/stats and everything else that still needs to be done here
bamFile <- loadBamFile(bamPath=bamPath)		#gunzipPath, libs: data.table, Rsamtools, testthat
indexFile <- getIndexFile(bamPath=bamPath)		#libs: data.table, Rsamtools, testthat
coverage <- getCoverage(outputCoveragePath=outputCoveragePath, new=newCoverage)		#gunzipPath, makeBed, getCoverageFromBed libs: data.table, testthat
varDb <- loadVarDb(inputVariantsDbPath=inputVariantsDbPath)		#gunzipPath, getFileExtension, loadVcfData, loadTsvData, normalizeDeletions, libs: data.table, testthat
positionsAll <- getPositions(varDb=varDb, new=newPostitionsAll)		#getUniqueLoci, libs: data.table, testthat

pileups <- doPileups(outputPositionsPath=outputPositionsAllPath,outputPileupsPath=outputPileupsPath,outputPileupsWithStatsPath=outputPileupsWithStatsPath,bamPath=bamPath,minBaseQuality=minBaseQuality,minMapq=minMapq,newPileups=newPileups,newPileupsWithStats=newPileupsWithStats)		#pileupClinVar, addStatsColumnsToPileup libs: data.table, testthat


if(train==TRUE){
  histogramData <- initializeHistogramData(lowerBounds = lowerSubsetBounds, upperBounds = upperSubsetBounds)		#libs: data.table
  histogramData <- updateDepthsColumn(subsetBounds=histogramData, updateDp=TRUE, updateDpRand=FALSE, pileupsWithNucCols=pileups)

  positionsAllFalse <- getFalsePositions(outputPositionsAllFalsePath=outputPositionsAllFalsePath, new=newPositionsAllFalse, coverage=coverage,
  										variantsSetToExclude=positionsAll, histogramData=histogramData)		#getRandomPositionsExcludeSubset2, libs: data.table, testthat
  print("false pileups:")
  falsePileups <- doPileups(outputPositionsPath=outputPositionsAllFalsePath, outputPileupsPath=outputPileupsFalsePath,outputPileupsWithStatsPath=outputFalsePileupsWithStatsPath,bamPath=bamPath,minBaseQuality=minBaseQuality,minMapq=minMapq,new=newFalsePileups,newPileupsWithStats=newFalsePileupsWithStats)		#pileupClinVar, addStatsColumnsToPileup libs: data.table, testthat
  histogramData <- updateDepthsColumn(subsetBounds=histogramData, updateDp=TRUE, updateDpRand=TRUE, pileupsWithNucCols=pileups, falsePileupsWithNucCols = falsePileups)
  #createHistogram(histogramData=histogramData, pileups=pileups, falsePileups=falsePileups)		#initializeHistogramData, getSubsetDepths, libs: data.table

  trainingData <- prepareTrainingData(pileups=pileups, falsePileups=falsePileups)
  model <- trainModel(modelPath=modelPath, data=trainingData,listOfTrainArguments=listOfTrainArguments, mc=0)
  
  #test the model
  baseline <- classifyNaive(newData = trainingData)
  baselineStats <- getClassificationStats(subsetBounds = histogramData, classificationData = baseline)
  newModelClassificationData <- classify(model=model, newData = trainingData)
  newModelStats <- getClassificationStats(subsetBounds = histogramData, classificationData = newModelClassificationData)
  
  
  eval <- evaluateClassifier(naiveDataStats = baselineStats, classifierDataStats = newModelStats)
  #eval <- evaluateClassifier(naiveDataStats = newModelStats, classifierDataStats = baselineStats)
  plot <- evaluateClassifierPlot(evaluationData = eval)
  #get the regular classification result
  classificationResult <- classify(model=model, newData=pileups)
}else{
  classificationResult <- classify(model=modelPath, newData=pileups)
}
varDb.merged <- mergePileupsWithVarDb(pileups=classificationResult, varDb=varDb, allx=TRUE,ally=TRUE)
if(verbose==TRUE){print("Pileup merged with varDb.")}
varDb.merged <- calculateDepths(varDb.merged)
if(verbose==TRUE){print("Depths calculated.")}

#TODO: HC fc/general should be loaded and processed here
if(HCTest==TRUE){
  vcfHcGeneral <- getVcfHcGeneral(vcfHcGeneralPath=vcfHcGeneralPath)
  vcfHcFc <- getVcfHcFc(vcfHcFcPath=vcfHcFcPath)
}






















#TODO: tidy up everything below this line, put into functions, place before the main part of the script etc.
#TODO: correct stats calculation
if(returnStats==TRUE){
  FNPlusTP <- nrow(varDb.merged) # P = TP + FN, total number of positives
  ccTP <- sum(varDb.merged[,result=="variant"], na.rm=TRUE) # pileups on P set only; 
  ccFN <- sum(varDb.merged[,result=="nonVariant"], na.rm=TRUE) #Error in eval(jsub, SDenv, parent.frame()) : object 'result' not found

  # ===== TPR ===== # sensitivity: TPR = TP/(TP+FN)
  ccTPR <- ccTP/FNPlusTP # ClinvarCaller, forcecalling
  
  if(HCTest==TRUE){
    hcGeneralTP <- nrow(vcfHcGeneral[CHROM %in% varDb$CHROM & POS %in% varDb$POS])
    hcFcTP <- nrow(vcfHcFc[CHROM %in% varDb$CHROM & POS %in% varDb$POS]) # pileups on P set only, could be nrow(vcfHcFc)
    hcGeneralTPR <- hcGeneralTP/FNPlusTP # GATK HaplotypeCaller, general scan
    hcFcTPR <- hcFcTP/FNPlusTP # GATK HaplotypeCaller, forcecalling
  }

  #====================================================================== TPR ^ ======================================================================
  #====================================================================== FPR v ======================================================================
  FPPlusTN <- nrow(falsePileups) # N = FP + TN, total number of negatives
  
  #TODO: repair or delete
  #ccFP <- sum(falsePileups[,nucleotide!=REF])
  ccFP <- 0
  #ccTN <- sum(falsePileups[,nucleotide==REF])
  ccTN <- 0
  
  # ===== FPR ===== # miss rate: FPR = FP/(FP+TN)
  #ccFPR <- ccFP/FPPlusTN # ClinvarCaller, forcecalling
  ccFPR <- 0

#====================================================================== FPR ^ ======================================================================
#====================================================================== CC + HC sum v ======================================================================

  if(HCTest==TRUE){
    variantsCcPos<-subset(varDb.merged, result=="variant")
    variantsCcPos<-variantsCcPos[,c("CHROM","POS","REF","ALT")]
    variantsHcPos<-vcfHcFc[,c("CHROM","POS","REF","ALT")]
    ccHcOR<-rbind(variantsCcPos,variantsHcPos)
    setkeyv(ccHcOR, c("CHROM","POS"))
    ccHcOR <- subset(unique(ccHcOR, by=c("CHROM","POS")))
    ccHcAND <- variantsHcPos[CHROM %in% variantsCcPos$CHROM & POS %in% variantsCcPos$POS] # variants in both cc and hc
    variantsHcOnly <- ccHcOR[!variantsCcPos]
    variantsCcOnly <- ccHcOR[!variantsHcPos]
    rm(variantsCcPos,variantsHcPos)
    
    variantsHcOnlyCount <- nrow(variantsHcOnly)
    variantsCcOnlyCount <- nrow(variantsCcOnly)
    variantsHcOnlyInDelCount <- nrow(variantsHcOnly[nchar(as.character(ALT))!=nchar(as.character(REF)) & substr(as.character(REF),1,1)==substr(as.character(ALT),1,1)]) + nrow(variantsHcOnly[variantsHcOnly[,ALT=="-"]])
    variantsCcOnlyInDelCount <- nrow(variantsCcOnly[nchar(as.character(ALT))!=nchar(as.character(REF)) & substr(as.character(REF),1,1)==substr(as.character(ALT),1,1)]) + nrow(variantsCcOnly[variantsCcOnly[,ALT=="-"]])
    variantsHcOnlySnpCount <- nrow(variantsHcOnly[nchar(as.character(ALT))==nchar(as.character(REF)) & variantsHcOnly[,ALT!="-"]])
    #unique(variantsHcOnly[nchar(as.character(ALT))==nchar(as.character(REF))][,REF]) # SNP check
    variantsCcOnlySnpCount <- nrow(variantsCcOnly[nchar(as.character(ALT))==nchar(as.character(REF)) & variantsCcOnly[,ALT!="-"]])
    #unique(variantsCcOnly[nchar(as.character(ALT))==nchar(as.character(REF))][,REF]) # SNP check
  }
#====================================================================== CC + HC sum ^ ======================================================================
  
  stats <- data.table(FNPlusTP=FNPlusTP,ccTP=ccTP,ccFN=ccFN,ccTPR=ccTPR,FPPlusTN=FPPlusTN,ccFP=ccFP,ccTN=ccTN,ccFPR=ccFPR)
  
  if(HCTest==TRUE){
    stats <- cbind(stats,data.table(hcGeneralTP=hcGeneralTP,hcFcTP=hcFcTP,hcGeneralTPR=hcGeneralTPR,hcFcTPR=hcFcTPR,variantsHcOnlyCount=variantsHcOnlyCount,
    variantsCcOnlyCount=variantsCcOnlyCount,variantsHcOnlyInDelCount=variantsHcOnlyInDelCount,
    variantsCcOnlyInDelCount=variantsCcOnlyInDelCount,variantsHcOnlySnpCount=variantsHcOnlySnpCount,
    variantsCcOnlySnpCount=variantsCcOnlySnpCount))
  }
  View(stats)
  fwrite(stats, outputStatsPath)
  if(verbose==TRUE){print("Stats ready.")}
}

fwrite(varDb.merged[,c("CHROM","POS","REF","ALT","countA","countC","countG","countT","countIn","countDel","RD","AD","nullDepth","DP","errFlag")], outputVarDbMergedPath)
if(verbose==TRUE){print("varDb.merged exported to external file.")}




res <- list()
if(returnPileups==TRUE){res <- append(res,list(pileups=pileups, falsePileups=falsePileups))}
if(returnPositionsAll==TRUE){res <- append(res, list(positionsAll=positionsAll))}
if(returnPositionsAllFalse==TRUE){res <- append(res, list(positionsAllFalse=positionsAllFalse))}
if(returnStats==TRUE){res <- append(res,list(stats=stats))}
if(returnMergedDb==TRUE){res <- append(res,list(varDb.merged=varDb.merged))}
if(returnHistogramData==TRUE){res <- append(res,list(histogramData=histogramData))}
  
return(res)

}