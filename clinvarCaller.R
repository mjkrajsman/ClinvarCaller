library(data.table)
library(Rsamtools)
library(tictoc)
library(stringr)
library(testthat)
library(BSgenome.Hsapiens.UCSC.hg38)
#library(BSgenome.Hsapiens.UCSC.hg19)
source("./tests.R")

clinvarCaller <- function(newPileups=TRUE,newFalsePileups=TRUE,newLociAllFalse=TRUE,runsAtImid=FALSE,hg38=TRUE,returnMergedDb=TRUE,
                          returnStats=TRUE,returnPileups=TRUE,returnLociAll=TRUE,returnLociAllFalse=TRUE,returnHistogramData=TRUE,
                          verbose=TRUE,minBaseQuality=0,minMapq=0,fullPileupStats=TRUE,varDbFormat="tsv",chromosomeLengthsPath="./input/chromosomeLengths.csv",
                          bamPath="./input/corriell_S7-ready.bam",bedPath="./input/corriell_S7-ready.per-base.bed.gz",
                          inputVariantsDbPath="./input/clinvar_alleles.single.b38.tsv.gz",vcfHcGeneralPath="./input/corriell_S7-HCScan.vcf.gz",
                          vcfHcFcPath="./input/corriell_S7-HCScanFC.vcf.gz",
                          outputLociAllPath="./output/lociAll.positions",
                          outputLociAllFalsePath="./output/lociAllFalse.positions",
                          outputPileupsPath="./output/pileups.mpileup",
                          outputPileupsFalsePath="./output/pileupsFalse.mpileup",
                          outputStatsPath="./output/stats.csv",
                          outputVarDbMergedPath="./output/varDbMerged-reduced.csv"){

if(hg38==TRUE){
library(BSgenome.Hsapiens.UCSC.hg38)
}else{
library(BSgenome.Hsapiens.UCSC.hg19)
}
# # data path
# if(runsAtImid==TRUE){
# inputPath <- "/ibm-storage/NGS/results/clinvar_caller/input/"
# outputPath <- "/ibm-storage/NGS/results/clinvar_caller/output/"
# }else{
# inputPath <- "./input/"
# outputPath <- "./output/"
# }
  
getChromLength <- function(chrom){
  chromLengths <- fread(chromosomeLengthsPath)
  return(chromLengths[CHROM==chrom,POS_MAX])
}

gunzipPath <- function(string){
  if(str_sub(string,-3,-1)==".gz"){
    res <- paste0("gunzip -c ",string)
  }else{
    res <- string
  }
  return(res)
}


randomPositionsChrom <- function(chrom, count){
  res <- data.table("POS"=sample(1:getChromLength(chrom), count, replace=F))
  res[,"CHROM":=chrom]   
  res <- res[, c("CHROM","POS")]
  res <- res[order(POS)]
  return(res)
}

#count - samples per chromosome
#variantsSetToExclude - set of variants to exclude from selection
getRandomPositionsExcludeSubset <- function(positionsCount, variantsSetToExclude){ 
  chromLengths <- fread(chromosomeLengthsPath)
  chroms <- lapply(chromLengths$CHROM, function(i){
    #print(paste0("chr",c(i)))
    # draw, check, draw again until success
    variantsSubset<-subset(variantsSetToExclude, CHROM==i)
    tmpPositionsCount <- positionsCount
    resRandom <- randomPositionsChrom(i,tmpPositionsCount)
    repeat{
      #print(paste0("tmpCount: ",c(tmpCount)))
      # remove loci present in variantsSetToExclude
      resRandom <- resRandom[!variantsSubset, on=.(POS)]
      tmpPositionsCount <- positionsCount - nrow(resRandom)
      if(nrow(resRandom)==positionsCount){
        break
      }
      resRandom <- rbind(resRandom,randomPositionsChrom(i,tmpPositionsCount))
    }
    return(resRandom)
  })
  res <- rbindlist(chroms)
  return(res)
}

loadVcfData <- function(vcfPath=""){
  res <- fread(vcfPath, skip = "CHROM")
  colnames(res)[colnames(res)=="#CHROM"] <- "CHROM"
  res$CHROM <- gsub("chr", "", res$CHROM)
  res <- res[!like(ALT,",")] # remove multiallelic sites
  res <- subset(unique(res, by=c("CHROM","POS")))
  res[CHROM == "MT" ,"CHROM":="M"]
  return(res)
}

loadClinvarData <- function(vcfPath=""){
  res <- fread(vcfPath)
  colnames(res)[colnames(res)=="chrom"] <- "CHROM"
  colnames(res)[colnames(res)=="pos"] <- "POS"
  colnames(res)[colnames(res)=="ref"] <- "REF"
  colnames(res)[colnames(res)=="alt"] <- "ALT"
  res[CHROM == "MT" ,"CHROM":="M"]
  #normalizeDeletions(res)
  return(res)
}

normalizeDeletions <- function(vcfData=""){
  res <- copy(vcfData)
  res[nchar(as.character(ALT))==1 & nchar(as.character(REF))>1 & substr(as.character(REF),1,1)==as.character(ALT),`:=`("POS"=POS+1L, "ALT"="-", "REF"=substring(as.character(REF),2))]
  return(res)
}

setTenthColumnNameInVcfToGenotype <- function(vcfTable=""){
  colnames(vcfTable)[10] <- "GENOTYPE"
  return(vcfTable)
}



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

mergeLoci <- function(loci1, loci2){
  res <- rbind(loci1, loci2)
  res[order("CHROM","POS")]                            						# sort by chrom, pos
  return(res)
}

#getAllDelLoci <- function(preceedingLoci, saveNucleotideData=FALSE){
#  res <- preceedingLoci[,c("CHROM","POS","lengthDiff","REF","ALT")]
#  #for every loci with lenghtDiff>0 add
#}

# filterChroms <- function(input, removeY=FALSE, removeM=FALSE){
  # res <- lapply(1:22, function(i){
    # x<-subset(input, CHROM==i)
    # return(x)
  # })
  # res <- rbindlist(res)
  
  # XX<-subset(input, CHROM=="X")
  # YY<-subset(input, CHROM=="Y")
  # MM<-subset(input, CHROM=="M")
  # MTMT<-subset(input, CHROM=="MT")
  
  # res <- na.omit(rbind(res, XX))
  # res <- na.omit(rbind(res, YY))
  # res <- na.omit(rbind(res, MM))
  # res <- na.omit(rbind(res, MTMT))
  # rm(XX,YY,MM,MTMT)
  
  # if(removeY==TRUE){
    # res <- res[!res[,CHROM=="Y"]]
  # }
  # if(removeM==TRUE){
    # res <- res[!res[,CHROM=="M"]]
    # res <- res[!res[,CHROM=="MT"]]
  # }
  # return(res)
# }

filterChroms <- function(input){
  res <- input[which(input$CHROM %in% c(1:22, "X", "Y", "M")),]
  return(res)
}

#================================================================== Pileups v ==================================================================



pileupClinVar <- function(lociPath="",pileupPath="",bamPath="",newPileups=TRUE){
  if(newPileups==TRUE){
    cmd <- paste0("samtools mpileup -A -d 0 -l ",lociPath," -q 0 -Q 0 -O -s ",bamPath," > ",pileupPath)
    system(cmd)
  }
  pileups <- fread(pileupPath, quote="")
  colnames(pileups) <- c("CHROM","POS","REF","DP","NUC","BQ","MQ","posRead")
  pileups[,"REF":=NULL]
  pileups$CHROM <- gsub("chr", "", pileups$CHROM)
  hg38=TRUE
  pileups <- addRefFromReferenceGenome(pileups)
  #pileups$NUC <- gsub("$", "", pileups$NUC)
  pileups$NUC <- str_replace_all(pileups$NUC,"\\^.","")
  for(n in 1:9){
    pileups$NUC <- str_replace_all(pileups$NUC,paste0("\\-",n,"[ACGTNacgtn]{",n,"}"),"") # deletions: scan next position and look for "*" (remove "-" from this read)
    pileups$NUC <- str_replace_all(pileups$NUC,paste0("[ACGTNacgtn]{1}\\+",n,"[ACGTNacgtn]{",n,"}"),"+") # insertions: change base to "+"
  }
  pileups$NUC <- str_replace_all(pileups$NUC,"[^[ACGTNacgtn*-0123456789+]]","")
  pileups$NUC <- str_replace_all(pileups$NUC,"\\*","-")
  pileups$NUC <- toupper(pileups$NUC)
  #pileups[,"NCHAR":=(nchar(NUC)==DP)]

  # encode row
  #pileups[,"posRead2":=phredQualityEncode(as.list(strsplit(posRead,",")[[1]])), by=1:nrow(pileups)]
  # decode row
  #pileups[,"posRead3":=paste(phredQualityDecode(posRead2),sep=",",collapse = ","), by=1:nrow(pileups)]
  
  # decode BQ, MQ
  pileups[,"BQ":=paste(phredQualityDecode(BQ),sep=",",collapse = ","), by=1:nrow(pileups)]
  pileups[,"MQ":=paste(phredQualityDecode(MQ),sep=",",collapse = ","), by=1:nrow(pileups)]
  return(pileups)
}
phredQualityDecode <- function(n){
  return(utf8ToInt(as.character(n))-33)
}

phredQualityEncode <- function(n){
  return(intToUtf8(as.integer(n)+33))
}

calculateNucStats <- function(chrom="", pos="", nuc="", countsOnly=FALSE, bq="", mq="", pr=""){
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
      CHROM=unlist(chrom)
      ,POS=unlist(pos)
      ,NUC=strsplit(unlist(nuc),"")[[1]]
      ,BQ=as.numeric(strsplit(unlist(bq),",")[[1]]) # BQ decoded!
      ,MQ=as.numeric(strsplit(unlist(mq),",")[[1]]) # MQ decoded!
      ,posRead=as.numeric(strsplit(unlist(pr),",")[[1]])
    )
    
    summaryDt <- data.table(CHROM=dt[1]$CHROM,POS=dt[1]$POS
                            ,nucCount_A=nrow(dt[NUC=="A",]),meanBQ_A=dt[NUC=="A",mean(BQ)],sdBQ_A=dt[NUC=="A",sd(BQ)],meanMQ_A=dt[NUC=="A",mean(MQ)]
                            ,sdMQ_A=dt[NUC=="A",sd(MQ)],meanPosRead_A=dt[NUC=="A",mean(posRead)],sdPosRead_A=dt[NUC=="A",sd(posRead)]
                            ,nucCount_C=nrow(dt[NUC=="C",]),meanBQ_C=dt[NUC=="C",mean(BQ)],sdBQ_C=dt[NUC=="C",sd(BQ)],meanMQ_C=dt[NUC=="C",mean(MQ)]
                            ,sdMQ_C=dt[NUC=="C",sd(MQ)],meanPosRead_C=dt[NUC=="C",mean(posRead)],sdPosRead_C=dt[NUC=="C",sd(posRead)]  
                            ,nucCount_G=nrow(dt[NUC=="G",]),meanBQ_G=dt[NUC=="G",mean(BQ)],sdBQ_G=dt[NUC=="G",sd(BQ)],meanMQ_G=dt[NUC=="G",mean(MQ)]
                            ,sdMQ_G=dt[NUC=="G",sd(MQ)],meanPosRead_G=dt[NUC=="G",mean(posRead)],sdPosRead_G=dt[NUC=="G",sd(posRead)]  
                            ,nucCount_T=nrow(dt[NUC=="T",]),meanBQ_T=dt[NUC=="T",mean(BQ)],sdBQ_T=dt[NUC=="T",sd(BQ)],meanMQ_T=dt[NUC=="T",mean(MQ)]
                            ,sdMQ_T=dt[NUC=="T",sd(MQ)],meanPosRead_T=dt[NUC=="T",mean(posRead)],sdPosRead_T=dt[NUC=="T",sd(posRead)]
                            ,nucCount_In=nrow(dt[NUC=="+",]),meanBQ_In=dt[NUC=="+",mean(BQ)],sdBQ_In=dt[NUC=="+",sd(BQ)],meanMQ_In=dt[NUC=="+",mean(MQ)]
                            ,sdMQ_In=dt[NUC=="+",sd(MQ)],meanPosRead_In=dt[NUC=="+",mean(posRead)],sdPosRead_In=dt[NUC=="+",sd(posRead)]
                            ,nucCount_Del=nrow(dt[NUC=="-",]),meanBQ_Del=dt[NUC=="-",mean(BQ)],sdBQ_Del=dt[NUC=="-",sd(BQ)],meanMQ_Del=dt[NUC=="-",mean(MQ)]
                            ,sdMQ_Del=dt[NUC=="-",sd(MQ)],meanPosRead_Del=dt[NUC=="-",mean(posRead)],sdPosRead_Del=dt[NUC=="-",sd(posRead)]
    )
    summaryDt[is.na(summaryDt)] <- 0 
    return(summaryDt)
  }
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

addStatsColumnsToPileup <- function(pileups, countsOnly=FALSE, rmCols=TRUE){
  res <- copy(pileups)
  tic(paste0("Stats calculated. Pileups, nrow: ", nrow(res), " time"))
  if(countsOnly==TRUE){
    res[,c("CHROM","POS","nucCount_A","nucCount_C","nucCount_G","nucCount_T","nucCount_In","nucCount_Del"):=calculateNucStats(chrom=CHROM, pos=POS, nuc=NUC, countsOnly=TRUE, bq=BQ, mq=MQ, pr=posRead), by=1:nrow(res)]
  }else{
    res[,c("CHROM","POS","nucCount_A","meanBQ_A","sdBQ_A","meanMQ_A","sdMQ_A","meanPosRead_A","sdPosRead_A","nucCount_C","meanBQ_C","sdBQ_C","meanMQ_C","sdMQ_C","meanPosRead_C","sdPosRead_C","nucCount_G","meanBQ_G","sdBQ_G","meanMQ_G","sdMQ_G","meanPosRead_G","sdPosRead_G","nucCount_T","meanBQ_T","sdBQ_T","meanMQ_T","sdMQ_T","meanPosRead_T","sdPosRead_T","nucCount_In","meanBQ_In","sdBQ_In","meanMQ_In","sdMQ_In","meanPosRead_In","sdPosRead_In","nucCount_Del","meanBQ_Del","sdBQ_Del","meanMQ_Del","sdMQ_Del","meanPosRead_Del","sdPosRead_Del"):=
            calculateNucStats(chrom=CHROM, pos=POS, nuc=NUC, countsOnly=FALSE, bq=BQ, mq=MQ, pr=posRead), by=1:nrow(res)]
  }  
  toc()
  if(rmCols==TRUE){
    res[, `:=`(NUC = NULL,BQ = NULL,MQ = NULL,posRead = NULL)]
  }
  return(res)
}



mergePileupsWithVarDb <- function(pileups="", varDb="",allx=TRUE,ally=TRUE){
  tmp <- copy(pileups)
  if("REF" %in% colnames(tmp)){tmp[, `:=`(REF = NULL)]}
  if("NUC" %in% colnames(tmp)){tmp[, `:=`(NUC = NULL)]}
  if("BQ" %in% colnames(tmp)){tmp[, `:=`(BQ = NULL)]}
  if("MQ" %in% colnames(tmp)){tmp[, `:=`(MQ = NULL)]}
  if("posRead" %in% colnames(tmp)){tmp[, `:=`(posRead = NULL)]}
  #print(colnames(pileups))
  #print(colnames(varDb))
  res <- merge(varDb, tmp, by=c("CHROM","POS"), all.x=allx, all.y=ally)
  rm(tmp)
  res[,c("nucCount_A","nucCount_C","nucCount_G","nucCount_T","nucCount_In","nucCount_Del")][is.na(res[,c("nucCount_A","nucCount_C","nucCount_G","nucCount_T","nucCount_In","nucCount_Del")])] <- 0 
  return(res)
}


#================================================================== Pileups ^ ==================================================================
#================================================================== varDb.merged v ==================================================================

calculateDepths <- function(variantsTable=""){
  
  # nucleotides matching REF
  #if(!("RD" %in% colnames(variantsTable))){
  #  variantsTable[, `:=`(RD = integer(0))]
  #  }
  if("RD" %in% colnames(variantsTable)){variantsTable[, `:=`(RD = NULL)]}
  variantsTable[, `:=`(RD = integer(0))]
  variantsTable[,c("RD")][is.na(variantsTable[,c("RD")])] <- 0
  variantsTable[substr(REF,1,1) == "A" ,"RD":=RD+nucCount_A]
  variantsTable[substr(REF,1,1) == "C" ,"RD":=RD+nucCount_C]
  variantsTable[substr(REF,1,1) == "G" ,"RD":=RD+nucCount_G]
  variantsTable[substr(REF,1,1) == "T" ,"RD":=RD+nucCount_T]
  ## variantsTable[nchar(as.character(alt))<nchar(as.character(ref)) & nchar(as.character(alt))>1, "RD":=RD+nucCount_In]
  ## variantsTable[nchar(as.character(alt))>nchar(as.character(ref)) & nchar(as.character(ref))>1, "RD":=RD+nucCount_Del]
  
  # nucleotides matching ALT
  #if(!("AD" %in% colnames(variantsTable))){
  #  variantsTable[, `:=`(AD = integer(0))]
  #  }
  if("AD" %in% colnames(variantsTable)){variantsTable[, `:=`(AD = NULL)]}
  variantsTable[, `:=`(AD = integer(0))]
  variantsTable[,c("AD")][is.na(variantsTable[,c("AD")])] <- 0
  variantsTable[ALT == "A" ,"AD":=AD+nucCount_A]
  variantsTable[ALT == "C" ,"AD":=AD+nucCount_C]
  variantsTable[ALT == "G" ,"AD":=AD+nucCount_G]
  variantsTable[ALT == "T" ,"AD":=AD+nucCount_T]
  #variantsTable[nchar(as.character(ALT))>nchar(as.character(REF)) & substr(as.character(REF),1,1)==substr(as.character(ALT),1,1),"AD":=AD+nucCount_In]
  #variantsTable[(nchar(as.character(ALT))<nchar(as.character(REF)) & substr(as.character(REF),1,1)==substr(as.character(ALT),1,1)) | ALT:="-","AD":=AD+nucCount_Del]
  
  # nucleotides not matching REF or ALT
  #if(!("nullDepth" %in% colnames(variantsTable))){
  #  variantsTable[, `:=`(nullDepth = integer(0))]
  #  }
  if("nullDepth" %in% colnames(variantsTable)){variantsTable[, `:=`(nullDepth = NULL)]}
  variantsTable[, `:=`(nullDepth = integer(0))]
  variantsTable[,c("nullDepth")][is.na(variantsTable[,c("nullDepth")])] <- 0
  variantsTable[substr(REF,1,1) != "A" & ALT != "A" ,"nullDepth":=nullDepth+nucCount_A]
  variantsTable[substr(REF,1,1) != "C" & ALT != "C" ,"nullDepth":=nullDepth+nucCount_C]
  variantsTable[substr(REF,1,1) != "G" & ALT != "G" ,"nullDepth":=nullDepth+nucCount_G]
  variantsTable[substr(REF,1,1) != "T" & ALT != "T" ,"nullDepth":=nullDepth+nucCount_T]
  #variantsTable[nchar(as.character(ALT))<=nchar(as.character(REF)),"nullDepth":=nullDepth+nucCount_In]
  #variantsTable[nchar(as.character(ALT))>=nchar(as.character(REF)),"nullDepth":=nullDepth+nucCount_Del]
  
  # coverage depth
  variantsTable[,"DP":=nucCount_A+nucCount_C+nucCount_G+nucCount_T]
  #variantsTable[,"depthCtrl":=AD+RD+nullDepth]
  
  # error flag
  #variantsTable[DP!=depthCtrl, "errFlag":=1]
  variantsTable[nchar(as.character(REF))>1 & nchar(as.character(ALT))>1, "errFlag":=1]
  variantsTable[,c("errFlag")][is.na(variantsTable[,c("errFlag")])] <- 0 
  
  return(variantsTable)
}





#================================================================== varDb.merged ^ ==================================================================
#================================================================== tests v ==================================================================
dts_disjoint <- function(dt1,dt2){
  if(nrow(fsetdiff(dt1,dt2)) == nrow(dt1)){
    print('true')
  }else{
    print('false')
  }
}
#================================================================== tests ^ ==================================================================


#================================================================== randomPositions, histogram v ==================================================================
getCoverageFromBed <- function(bedPath="./input/corriell_S7-ready.per-base.bed.gz"){
  coverage <- fread(bedPath)
  colnames(coverage) <- c("CHROM","start","stop","readsCount")
  coverage <- subset(coverage, readsCount > 0)
  coverage$CHROM <- gsub("chr", "", coverage$CHROM)
  coverage <- filterChroms(coverage)
  coverage[,"start":=start+1L]
  #coverage[,"length":=stop-start]
  return(coverage)
}

getCoverageSubsetByReadsCount <- function(coverage, min=1, max=1){
  if(min<=max){
    return(subset(coverage, (readsCount >= min & readsCount <= max)))
  }else{
    print("max must be greater or equal to min!")
    return(NULL)
  }  
}



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

getRandomPositions <- function(listOfPositions, randomSubsetSize=10000){
  #print(randomSubsetSize)
  #print(nrow(listOfPositions))
  if(nrow(listOfPositions)>=randomSubsetSize){
    res <- listOfPositions[sample(.N, randomSubsetSize, replace=F)]
  }else{
    res <- listOfPositions  
  }
  return(res)
}

getRandomPositionsFromCoverage <- function(coverage,subsetBounds,subsetMaxLength=200000){
  res <- lapply(1:nrow(subsetBounds), function(i){
    coverageSubset <- getCoverageSubsetByReadsCount(coverage,min=subsetBounds[i,lower],max=subsetBounds[i,upper])
    listOfPositions <- getPositionsFromCoverageSubset(coverageSubset,subsetMaxLength)
    randomPositions <- getRandomPositions(listOfPositions, subsetBounds[i,DP])
    return(randomPositions)
  })
  res <- rbindlist(res)
  return(res)
}
#res -> nrow(subsetBounds) * randomSubsetSize

#getPositionsFromCoverage <- function(coverage,subsetBounds,subsetMaxLength=200000){
#  res <- lapply(1:nrow(subsetBounds), function(i){
#    coverageSubset <- getCoverageSubsetByReadsCount(coverage,min=subsetBounds[i,lower],max=subsetBounds[i,upper])
#    listOfPositions <- getPositionsFromCoverageSubset(coverageSubset,subsetMaxLength)
#    return(listOfPositions)
#  })
#  return(res)
#}

getVRF <- function(pileupsWithNucAndRefCols, getDP = TRUE, getRD = TRUE){
  res <- copy(pileupsWithNucAndRefCols)
  res[,"DP":=nucCount_A+nucCount_C+nucCount_G+nucCount_T+nucCount_In+nucCount_Del]
  
  if("RD" %in% colnames(res)){res[, `:=`(RD = NULL)]}
  res[, `:=`(RD = integer(0))]
  res[,c("RD")][is.na(res[,c("RD")])] <- 0
  res[substr(REF,1,1) == "A" ,"RD":=RD+nucCount_A]
  res[substr(REF,1,1) == "C" ,"RD":=RD+nucCount_C]
  res[substr(REF,1,1) == "G" ,"RD":=RD+nucCount_G]
  res[substr(REF,1,1) == "T" ,"RD":=RD+nucCount_T]
  
  res["DP"!=0,"VRF":=(RD/DP)]
  res["DP"==0,"VRF":=0]
  if(getRD==FALSE){res[, `:=`(RD = NULL)]}
  if(getDP==FALSE){res[, `:=`(DP = NULL)]}
  return(res)
}

getDP <- function(pileupsWithNucCols, removeNucColumns=FALSE){
  res <- copy(pileupsWithNucCols)
  res[,"DP":=nucCount_A+nucCount_C+nucCount_G+nucCount_T+nucCount_In+nucCount_Del]
  if(removeNucColumns==TRUE){
    res[, `:=`(nucCount_A = NULL, nucCount_C = NULL, nucCount_G = NULL, nucCount_T = NULL, nucCount_In = NULL, nucCount_Del = NULL)]
  }
  return(res)
}

getSubsetDepths <- function(pileupsWithNucCols, subsetBounds){
  res <- subsetBounds
  totalDepths <- getDP(pileupsWithNucCols, removeNucColumns=TRUE)
  res <- sapply(1:nrow(subsetBounds), function(i){
    return(nrow(subset(totalDepths, (DP >= subsetBounds[i,lower] & DP <= subsetBounds[i,upper]))))
  })
  return(res)
}
  
updateDepthsColumn <- function(pileupsWithNucCols, subsetBounds){
  res <- subsetBounds
  subsetDepths <- getSubsetDepths(pileupsWithNucCols, subsetBounds)
  res[,"DP":=subsetDepths]
  return(res)
}

#count - samples per chromosome
#variantsSetToExclude - set of variants to exclude from selection
getRandomPositionsExcludeSubset2 <- function(coverage,variantsSetToExclude,subsetBounds,returnReadsCount=FALSE){ 
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

#================================================================== randomPositions, histogram ^ ==================================================================












test_that("Does BAM file exist?", {
  expect_that(file.exists(bamPath), is_true())
})

bamFile <- BamFile(bamPath)
if(verbose==TRUE){print("BAM file loaded.")}
if(exists(paste0(bamPath, ".bai"))){
indexFile <- paste0(bamPath, ".bai")
  if(verbose==TRUE){print("Index (BAI) file loaded.")}
}else{
indexFile <- indexBam(bamFile)
if(verbose==TRUE){print("Index (BAI) file created and loaded.")}
}

seqinfo(bamFile)

test_that("Does BED file exist?", {
  expect_that(file.exists(bedPath), is_true())
})
bedPath <- gunzipPath(bedPath)
coverage <- getCoverageFromBed(bedPath)
if(verbose==TRUE){print("BED file loaded.")}

histogramData <- data.table("lower"=c(1,2,6,11,21,51,101,201),
                           "upper"=c(1,5,10,20,50,100,200,9999),
                           "DP"=c(1,1,1,1,1,1,1,1),
                           "DP_rand"=c(1,1,1,1,1,1,1,1))


#================================================================== varDB: Choose one v ==================================================================

inputVariantsDbPath <- gunzipPath(inputVariantsDbPath)

if(varDbFormat=="vcf"){
	varDb <- loadVcfData(inputVariantsDbPath)
	if(verbose==TRUE){print("varDb VCF data loaded.")}
}
if(varDbFormat=="tsv"){
	varDb <- loadClinvarData(inputVariantsDbPath) 
	if(verbose==TRUE){print("varDb TSV data loaded.")}
}
if(varDbFormat!="tsv" & varDbFormat!="vcf"){print("varDb file must be in tsv (clinvar) or vcf format!")}


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

#================================================================== varDB: Choose one ^ ==================================================================
# varDb; leave only unique chrom+pos
#snpLoci <- ggetUniqueLoci(varDb, variationType="snp", saveNucleotideData=TRUE)
#inLoci <- getUniqueLoci(varDb, variationType="in", saveNucleotideData=TRUE)
#delLoci <- getUniqueLoci(varDb, variationType="del", saveNucleotideData=TRUE)
#notMatchingLoci <- getUniqueLoci(varDb, variationType="notMatching", saveNucleotideData=TRUE)
lociAll <- getUniqueLoci(varDb)
lociAll2 <- copy(lociAll)
lociAll2[,"CHROM":=paste0("chr",CHROM)]
fwrite(lociAll2, outputLociAllPath, sep= " ")
rm(lociAll2)
if(verbose==TRUE){print("lociAll exported to external file.")}

#================================================================== Test sets v ==================================================================

if(returnStats==TRUE){
  # GATK HAPLOTYPECALLER - general scan
  test_that("Does vcfHcGeneral file exist?", {
    expect_that(file.exists(vcfHcGeneralPath), is_true())
  })
  
  vcfHcGeneralPath <- gunzipPath(vcfHcGeneralPath)
  
  vcfHcGeneral <- loadVcfData(vcfHcGeneralPath)
  if(verbose==TRUE){print("vcfHcGeneral file loaded.")}
  vcfHcGeneral <- filterChroms(vcfHcGeneral)
  
  test_that("Are only valid chromosomes (1:22, X, Y, M) present in vcfHcGeneral?", {
   expect_that(unique(unique(vcfHcGeneral$CHROM) %in% c(1:22, "X", "Y", "M")), is_true())
  })
  
  
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
  
  test_that("Are only 0/1 and 1/1 present in GENOTYPE?", {
    expect_that(unique(substr(vcfHcFc$GENOTYPE,1,3)), equals(c("0/1","1/1")))
  })
}

#================================================================== Test sets ^ ==================================================================

if(newPileups==TRUE){
  pileups <- pileupClinVar(outputLociAllPath,outputPileupsPath,bamPath,newPileups=TRUE)
  if(verbose==TRUE){print("pileups done.")}
  
}else{
  # ================================ B (no new falsePileups) ================================
  if(file.exists(outputPileupsPath)){
    pileups <- pileupClinVar(outputLociAllPath,outputPileupsPath,bamPath,newPileups=FALSE)
    if(verbose==TRUE){print("pileups restored from outputPileupsPath file.")}
  }else{
    stop("outputPileupsPath file is not present! Please provide a valid foutputPileupsPath or set newFalsePileups=TRUE")
  }
  
}


# ================================ pileups calculations
pileups <- addStatsColumnsToPileup(pileups, countsOnly = !(fullPileupStats), rmCols = TRUE)
if(verbose==TRUE){print("New columns created in pileups.")}


#pileupsTests(pileups)

#View(pileups)
varDb.merged <- mergePileupsWithVarDb(pileups, varDb, TRUE, TRUE)
if(verbose==TRUE){print("Pileup merged with varDb.")}
#varDb.merged <- varDb.merged[, c("CHROM","POS","REF","ALT","nucCount_A","nucCount_C","nucCount_G","nucCount_T","nucCount_In","nucCount_Del")]

#varDb.merged[, `:=`(RD = NULL, AD = NULL, nullDepth = NULL, DP = NULL, depthCtrl = NULL, errFlag = NULL)]

#fwrite(pileups, outputPileupsPath)
#if(verbose==TRUE){print("pileups exported to external file.")}

varDb.merged <- calculateDepths(varDb.merged)
#varDb.merged.f <- calculateDepths(varDb.merged.f)
if(verbose==TRUE){print("Depths calculated.")}
#View(varDb.merged)


histogramData <- updateDepthsColumn(pileups, histogramData)
if(verbose==TRUE){print("histogramData DP updated.")}


if(newLociAllFalse==TRUE){
  lociAllFalse <- getRandomPositionsExcludeSubset2(coverage, variantsSetToExclude=lociAll, histogramData,returnReadsCount=FALSE)
  if(verbose==TRUE){print("Random positions generated.")}
  
  lociAllFalse2 <- copy(lociAllFalse)
  lociAllFalse2[,"CHROM":=paste0("chr",CHROM)]
  fwrite(lociAllFalse2, outputLociAllFalsePath, sep= " ")
  rm(lociAllFalse2)
  if(verbose==TRUE){print("lociAllFalse exported to external file.")}
  
  test_that('Are lociAllFalse and lociAll disjoint?', {
    expect_that(dts_disjoint(lociAllFalse[,c("CHROM","POS")],lociAll[,c("CHROM","POS")]), prints_text('true'))
  })
  
  #falsePileups <- pileupClinVar(lociAllFalse, 25000) # loci w/o variants - FPR
}else{
  if(file.exists(outputLociAllFalsePath)){
    lociAllFalse <- fread(outputLociAllFalsePath)
    lociAllFalse$CHROM <- gsub("chr", "", lociAllFalse$CHROM)
    if(verbose==TRUE){print("lociAllFalse restored from outputLociAllFalsePath file.")}
  }else{
    stop("outputLociAllFalsePath file is not present! Please provide a valid outputLociAllFalsePath or set newLociAllFalse=TRUE")
  }
}


if(newFalsePileups==TRUE){
  falsePileups <- pileupClinVar(outputLociAllFalsePath,outputPileupsFalsePath,bamPath,newPileups=TRUE)
	if(verbose==TRUE){print("falsePileups done.")}
  
}else{
  # ================================ B (no new falsePileups) ================================
  if(file.exists(outputPileupsFalsePath)){
    falsePileups <- pileupClinVar(outputLociAllFalsePath,outputPileupsFalsePath,bamPath,newPileups=FALSE)
    if(verbose==TRUE){print("falsePileups restored from outputPileupsFalsePath file.")}
  }else{
    stop("outputPileupsFalsePath file is not present! Please provide a valid outputPileupsFalsePath or set newFalsePileups=TRUE")
  }

}

# ================================ falsePileups calculations ================================ 
#falsePileups <- addRefFromReferenceGenome(falsePileups)
#if(verbose==TRUE){print("REF added to falsePileups.")}
falsePileups <- addStatsColumnsToPileup(falsePileups, countsOnly = !(fullPileupStats), rmCols = TRUE)
if(verbose==TRUE){print("New columns created in falsePileups.")}


#fwrite(falsePileups, outputPileupsFalsePath)
#if(verbose==TRUE){print("falsePileups exported to external file.")}


# ================================ histogram ================================ 
if(returnHistogramData==TRUE){
  histogramP <- getSubsetDepths(pileups, histogramData)
  histogramN <- getSubsetDepths(falsePileups, histogramData)
  if(verbose==TRUE){print("Histogram data ready.")}
  
  print(histogramP)
  print(histogramN)
  histogramData[,"DP_rand":=histogramN]
  if(verbose==TRUE){print("histogramData DP_rand updated.")}
}

# if(runsAtImid==TRUE){
	# #bamPath <- "/ibm-storage/NGS/analysis/Roche_RASGENODERM_V1/miseq reporter/44028_S1.bam"
	# getResFinal <- function(dd, varDb, lociAll){
	# res <- mclapply((1:length(dd)), function(i){ 
	  # bamPath <- dd[i]
	 
	  # print(paste0("sample nr=", i))
	  # bamFile <- BamFile(bamPath)
	  # indexFile <- paste0(bamPath, ".bai")
	  # seqinfo(bamFile)
	  # samples <- bamPath
	  
	  # pp <- PileupParam(include_insertions=T, distinguish_strand=F) #ignore_query=T -> count N

	  # pileups <- pileupClinVar(lociAll, 25000)
    # pileups <- addStatsColumnsToPileup(pileups)
	  # varDb.merged <- mergePileupsWithVarDb(pileups, varDb, TRUE, TRUE)
	  # varDb.merged <- calculateDepths(varDb.merged)
	  # varDb.merged$Identifier <- bamPath
	  # #varDb.merged[which(varDb.merged$pathogenic > 0)]
	  # varDb.merged
	 # }, mc.cores=40)

	# resFinal <- rbindlist(mclapply (res, function(x){x}, mc.cores=40))

	# resFinal
	# }


	# bamPath <- paste0(c(inputPath),"NA12878/",c(bamFilename))
	# varDb <- loadVcfData(paste0("gunzip -c ",c(inputPath),"NA12878/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz"))
	# varDb <- loadVcfData(paste0("gunzip -c ",c(inputPath),"NA12878/hg38.hybrid.vcf.gz"))
	# vcfHcGeneral <- loadVcfDataloadVcfData(paste0("gunzip -c ",c(inputPath),"NA12878/corriell_S7-gatk-haplotype.vcf.gz"))
	# vcfHcGeneral <- filterChroms(vcfHcGeneral)


	# varDb$key <- paste0(varDb$CHROM, ":" , varDb$POS, "_", varDb$REF, ">", varDb$ALT)
	# vcfHcGeneral$key <- paste0(vcfHcGeneral$CHROM, ":" , vcfHcGeneral$POS, "_", vcfHcGeneral$REF, ">", vcfHcGeneral$ALT)
	# lociAll <- getUniqueLoci(varDb)

	# pgCC <- getResFinal (bamPath, varDb,lociAll)

	# resFinal1 <- getResFinal(dd1)


	# dd <- dir("/ibm-storage/NGS/analysis/Roche_RASGENODERM_V1/miseq reporter", "*bam$", full.names=T)
	# dd1 <- dir("/ibm-storage/NGS/analysis/Roche_RASGENODERM_V2/miseq reporter", "*bam$", full.names=T)
	# dd2 <- dir("/ibm-storage/NGS/analysis/Roche_RASGENODERM_V2/miseq reporter/hyperplus", "*bam$", full.names=T)

	# dd <- dir ("/ibm-storage/NGS/analysis/Roche_HLEPI/COVARIS/miseq_reporter", "*bam$", full.names=T)
	# dd1 <- dir ("/ibm-storage/NGS/analysis/Roche_HLEPI/HYPERPLUS/miseq_reporter", "*bam$", full.names=T)
	# dd2 <- dir("/ibm-storage/NGS/analysis/Roche_HLEPI/HYPERPLUS/miseq_reporter/run36", "*bam$", full.names=T)

	# dd <-dir( "/ibm-storage/NGS/analysis/SureSelect_V6_CEGAT/CEGAT/run6/CEGAT-run6-Nov_14_2017/final/46015IMID_S572Nr6", "*bam$", full.names=T) 

	# dd <- dir("/ibm-storage/NGS/secondary_analysis/wiotkie/clinvar_caller/input", "*bam$", , full.names=T)
	# resFinalWiotkie <- getResFinal(dd)

	# resFinal1 <- getResFinal(dd1)
	# resFinal2 <- getResFinal(dd2)

	# fwrite(resFinal, paste0(c(outputPath),"rasgenoderm_v1_pathogenic.csv"))
	# fwrite(resFinal1, paste0(c(outputPath),"rasgenoderm_v2_pathogenic.csv"))

	# rr <- rbindlist(list(resFinal, resFinal1, resFinal2))
	# sel <- rr [which(   rr$AD > 0  & ((nchar(rr$REF) ==1)  | (rr$nucCount_Del > 0))),]
					
	# #sel [ grep("icht", tolower(sel$all_traits) ) ,]
	# #ids <- c("44028", "45688", "47986" , "48747",  "50939" , "51030",  "51649")
	# sel$key <- paste0(sel$CHROM, ":", sel$POS, "_", sel$REF, ">", sel$ALT)
	# keyDf <-  as.data.frame(table(sel$key))
	# sel$Freq <- keyDf$Freq [match(sel$key, keyDf$Var1)        ]
					# write.table(sel,  paste0(c(outputPath),"hlepi_306.xls"), sep="\t", row.names=F, quote=F)

	# sel$IDs <- sapply(strsplit(sel$Identifier, "[/_]"), function(x){x[9]})
	# which( grepl("icht", tolower(sel$all_traits) ))
# }




#====================================================================== TPR v ======================================================================
if(returnStats==TRUE){
  FNPlusTP <- nrow(varDb.merged) # P = TP + FN, total number of positives
  ccTP <- sum(varDb.merged[,AD>0]) # pileups on P set only; AD > 0 -> variant found # AD == 0 -> no variant
  ccFN <- sum(varDb.merged[,AD==0])
  hcGeneralTP <- nrow(vcfHcGeneral[CHROM %in% varDb$CHROM & POS %in% varDb$POS])
  
  hcFcTP <- nrow(vcfHcFc[CHROM %in% varDb$CHROM & POS %in% varDb$POS]) # pileups on P set only, could be nrow(vcfHcFc)
  
  # ===== TPR ===== # sensitivity: TPR = TP/(TP+FN)
  ccTPR <- ccTP/FNPlusTP # ClinvarCaller, forcecalling
  hcGeneralTPR <- hcGeneralTP/FNPlusTP # GATK HaplotypeCaller, general scan
  hcFcTPR <- hcFcTP/FNPlusTP # GATK HaplotypeCaller, forcecalling
  
  #====================================================================== TPR ^ ======================================================================
  #====================================================================== FPR v ======================================================================
  FPPlusTN <- nrow(falsePileups) # N = FP + TN, total number of negatives
  
  #ccFP <- sum(falsePileups[,nucleotide!=REF])
  ccFP <- 0
  #ccTN <- sum(falsePileups[,nucleotide==REF])
  ccTN <- 0
  
  # ===== FPR ===== # miss rate: FPR = FP/(FP+TN)
  #ccFPR <- ccFP/FPPlusTN # ClinvarCaller, forcecalling
  ccFPR <- 0

#====================================================================== FPR ^ ======================================================================
#====================================================================== CC + HC sum v ======================================================================


  variantsCcPos<-subset(varDb.merged, AD>0)
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
  variantsHcOnlyInDelCount <- nrow(variantsHcOnly[nchar(as.character(ALT))!=nchar(as.character(REF)) & substr(as.character(REF),1,1)==substr(as.character(ALT),1,1)])
  variantsCcOnlyInDelCount <- nrow(variantsCcOnly[nchar(as.character(ALT))!=nchar(as.character(REF)) & substr(as.character(REF),1,1)==substr(as.character(ALT),1,1)])
  variantsHcOnlySnpCount <- nrow(variantsHcOnly[nchar(as.character(ALT))==nchar(as.character(REF))])
  #unique(variantsHcOnly[nchar(as.character(ALT))==nchar(as.character(REF))][,REF]) # SNP check
  variantsCcOnlySnpCount <- nrow(variantsCcOnly[nchar(as.character(ALT))==nchar(as.character(REF))])
  #unique(variantsCcOnly[nchar(as.character(ALT))==nchar(as.character(REF))][,REF]) # SNP check

#====================================================================== CC + HC sum ^ ======================================================================
  
  stats <- data.table(FNPlusTP=FNPlusTP,ccTP=ccTP,ccFN=ccFN,hcGeneralTP=hcGeneralTP,hcFcTP=hcFcTP,ccTPR=ccTPR,hcGeneralTPR=hcGeneralTPR,
                hcFcTPR=hcFcTPR,FPPlusTN=FPPlusTN,ccFP=ccFP,ccTN=ccTN,ccFPR=ccFPR,variantsHcOnlyCount=variantsHcOnlyCount,
                variantsCcOnlyCount=variantsCcOnlyCount,variantsHcOnlyInDelCount=variantsHcOnlyInDelCount,
                variantsCcOnlyInDelCount=variantsCcOnlyInDelCount,variantsHcOnlySnpCount=variantsHcOnlySnpCount,
                variantsCcOnlySnpCount=variantsCcOnlySnpCount)
  
  fwrite(stats, outputStatsPath)
  if(verbose==TRUE){print("Stats ready.")}
}

fwrite(varDb.merged[,c("CHROM","POS","REF","ALT","nucCount_A","nucCount_C","nucCount_G","nucCount_T","nucCount_In","nucCount_Del","RD","AD","nullDepth","DP","errFlag")], outputVarDbMergedPath)
if(verbose==TRUE){print("varDb.merged exported to external file.")}




res <- list()
if(returnPileups==TRUE){res <- append(res,list(pileups=pileups, falsePileups=falsePileups))}
if(returnLociAll==TRUE){res <- append(res, list(lociAll=lociAll))}
if(returnLociAllFalse==TRUE){res <- append(res, list(lociAllFalse=lociAllFalse))}
if(returnStats==TRUE){res <- append(res,list(stats=stats))}
if(returnMergedDb==TRUE){res <- append(res,list(varDb.merged=varDb.merged))}
if(returnHistogramData==TRUE){res <- append(res,list(histogramData=histogramData))}
  
return(res)
}