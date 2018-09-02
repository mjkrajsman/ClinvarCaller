library(data.table)
library(Rsamtools)
library(tictoc)
library(stringr)
library(testthat)
library(BSgenome.Hsapiens.UCSC.hg38)


pileupClinVar <- function(lociPath="",pileupPath="",bamPath="",newPileups=TRUE,normalizeData=TRUE){
  if(newPileups==TRUE){
    cmd <- paste0("samtools mpileup -A -d 0 -l ",lociPath," -q 0 -Q 0 -O -s ",bamPath," > ",pileupPath)
    system(cmd)
  }
  pileups <- fread(pileupPath, quote="")
  colnames(pileups) <- c("CHROM","POS","REF","DP","NUC","BQ","MQ","PR")
  pileups[,"REF":=NULL]
  pileups$CHROM <- gsub("chr", "", pileups$CHROM)
  pileups <- addRefFromReferenceGenome(pileups)
  #pileups$NUC <- gsub("$", "", pileups$NUC)
  if(normalizeData==TRUE){  
    #maxReadLength <- max(as.numeric(pileups[,unlist(PR)]))
    pileups$NUC <- str_replace_all(pileups$NUC,"\\^.","")
    for(n in 1:99){
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

phredQualityDecode <- function(n){
  return(utf8ToInt(as.character(n))-33)
}

phredQualityEncode <- function(n){
  return(intToUtf8(as.integer(n)+33))
}




calculateNucStats <- function(chrom="", pos="",ref="", nuc="", countsOnly=FALSE, bq="", mq="", pr=""){
  if(countsOnly==TRUE){
    dt <- data.table(
      CHROM=unlist(chrom)
      ,POS=unlist(pos)
      ,NUC=strsplit(unlist(nuc),"")[[1]]
    )
    
    nucCounts <- data.table(CHROM=dt[1]$CHROM,POS=dt[1]$POS,nucCount_A=nrow(dt[NUC=="A",]),nucCount_C=nrow(dt[NUC=="C",])
                            ,nucCount_G=nrow(dt[NUC=="G",]),nucCount_T=nrow(dt[NUC=="T",]),nucCount_In=nrow(dt[NUC=="+",])
                            ,nucCount_Del=nrow(dt[NUC=="-",])
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
                            countT=nrow(dt[NUC=="T",]),countIn=nrow(dt[NUC=="+",]),countDel=nrow(dt[NUC=="-",]),
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

addStatsColumnsToPileup <- function(pileups, countsOnly=FALSE, rmCols=FALSE){
  res <- copy(pileups)
  tic(paste0("Stats calculated. Pileups, nrow: ", nrow(res), " time"))
  if(countsOnly==TRUE){
    res[,c("CHROM","POS","countA","countC","countG","countT","countIn","countDel"):=
          calculateNucStats(chrom=CHROM, pos=POS, nuc=NUC, countsOnly=TRUE, bq=BQ, mq=MQ, pr=PR), by=1:nrow(res)]
  }else{
    #View(res)
    res[,c(#"CHROM","POS","REF",
           #"NUC","BQ","MQ","PR",
           "countA","countC","countG","countT","countIn","countDel",
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
  return(res)
}

pileups <- pileupClinVar(lociPath="./output/positionsAll.positions",pileupPath = "./output/pileups.mpileup",bamPath="./input/corriell_S7-ready.bam",newPileups=TRUE,normalizeData=TRUE)
rowPlusStats <- addStatsColumnsToPileup(pileups = pileups[3990:4000],countsOnly = FALSE,rmCols = FALSE)

# pileupsBackup <- pileupClinVar(lociPath="./output/positionsAll.positions",pileupPath = "./output/pileups.mpileup",bamPath="./input/corriell_S7-ready.bam",newPileups=TRUE,normalizeData=TRUE)
# pileups <- pileupClinVar(lociPath="./output/positionsAll.positions",pileupPath = "./output/pileups.mpi1leup",bamPath="./input/corriell_S7-ready.bam",newPileups=TRUE,normalizeData=FALSE)
# pileups <- pileupsBackup[1:1000]
# pileups[,`:=`(NUC=strsplit(NUC,""),BQ=strsplit(BQ,","),MQ=strsplit(MQ,","),PR=strsplit(PR,","))]
# 
# 
# #pileups <- pileupClinVar("./output/positionsAll.positions","./output/pileups.mpileup","./input/corriell_S7-ready.bam",newPileups=TRUE)
# 
# row <- pileups[100]
# rowMod <- copy(row)
# #rowMod$NUC <- paste0(rowMod$NUC,rowMod$NUC)
# rowMod$NUC <- "GAGAGAGGGGAAGGGATCCCGAGGGGAAGGGA"
# rowMod$DP <- nchar("GAGAGAGGGGAAGGGATCCCGAGGGGAAGGGA")
# rowMod$BQ <- paste0(rowMod$BQ,",",rowMod$BQ)
# rowMod$MQ <- paste0(rowMod$MQ,",",rowMod$MQ)
# rowMod$PR <- paste0(rowMod$PR,",",rowMod$PR)
# 
# rowPlusStats <- addStatsColumnsToPileup(pileups = row,countsOnly = FALSE,rmCols = FALSE)
# rowModPlusStats <- addStatsColumnsToPileup(pileups = rowMod,countsOnly = FALSE,rmCols = FALSE)
# 
# dt <- data.table(CHROM=unlist(row$CHROM)
#                  ,POS=unlist(row$POS)
#                  ,REF=unlist(row$REF)
#                  ,NUC=strsplit(unlist(row$NUC),"")[[1]]
#                  ,BQ=as.numeric(strsplit(unlist(row$BQ),",")[[1]]) # BQ decoded!
#                  ,MQ=as.numeric(strsplit(unlist(row$MQ),",")[[1]]) # MQ decoded!
#                  ,PR=as.numeric(strsplit(unlist(row$PR),",")[[1]])
# )
# dtMod <- data.table(CHROM=unlist(rowMod$CHROM)
#                  ,POS=unlist(rowMod$POS)
#                  ,REF=unlist(rowMod$REF)
#                  ,NUC=strsplit(unlist(rowMod$NUC),"")[[1]]
#                  ,BQ=as.numeric(strsplit(unlist(rowMod$BQ),",")[[1]]) # BQ decoded!
#                  ,MQ=as.numeric(strsplit(unlist(rowMod$MQ),",")[[1]]) # MQ decoded!
#                  ,PR=as.numeric(strsplit(unlist(rowMod$PR),",")[[1]])
# )
# 
# #pileupsWithCols <- addStatsColumnsToPileup(pileups, countsOnly=FALSE, rmCols=FALSE) 
# 
# 
# rowPlusStats <- addStatsColumnsToPileup(pileups = pileups[1:1000],countsOnly = FALSE,rmCols = FALSE)
