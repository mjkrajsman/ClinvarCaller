library(data.table)
library(Rsamtools)
library(tictoc)
library(stringr)
library(testthat)
library(BSgenome.Hsapiens.UCSC.hg38)
#library(BSgenome.Hsapiens.UCSC.hg19)
source("./tests.R")

gunzipPath <- function(string){
  if(str_sub(string,-3,-1)==".gz"){
    res <- paste0("gunzip -c ",string)
  }else{
    res <- string
  }
  return(res)
}

filterChroms <- function(input){
  res <- input[which(input$CHROM %in% c(1:22, "X", "Y", "M")),]
  return(res)
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
    randomPositions <- getRandomPositions(listOfPositions, subsetBounds[i,TD])
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

getTotalDepths <- function(pileupsWithNucCols){
  res <- copy(pileupsWithNucCols)
  res[,"TD":=nucA+nucC+nucG+nucT]
  res[, `:=`(nucA = NULL, nucC = NULL, nucG = NULL, nucT = NULL, nucIn = NULL, nucDel = NULL)]
  return(res)
}

getSubsetDepths <- function(pileupsWithNucCols, subsetBounds){
  res <- subsetBounds
  totalDepths <- getTotalDepths(pileupsWithNucCols)
  res <- sapply(1:nrow(subsetBounds), function(i){
    return(nrow(subset(totalDepths, (TD >= subsetBounds[i,lower] & TD <= subsetBounds[i,upper]))))
  })
  return(res)
}
  
updateDepthsColumn <- function(pileupsWithNucCols, subsetBounds){
  res <- subsetBounds
  subsetDepths <- getSubsetDepths(pileupsWithNucCols, subsetBounds)
  res[,"TD":=subsetDepths]
  return(res)
}

#count - samples per chromosome
#variantsSetToExclude - set of variants to exclude from selection
getRandomPositionsExcludeSubset2 <- function(coverage,variantsSetToExclude,subsetBounds){ 
  posCount <- sum(subsetBounds[,"TD"])
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
  res <- res[sample(.N, posCount, replace=F),c("CHROM","POS","readsCount")]
  return(res)
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

loadVcfData <- function(vcfPath=""){
  res <- fread(vcfPath, skip = "CHROM")
  colnames(res)[colnames(res)=="#CHROM"] <- "CHROM"
  res$CHROM <- gsub("chr", "", res$CHROM)
  res <- res[!like(ALT,",")] # remove multiallelic sites
  res <- subset(unique(res, by=c("CHROM","POS")))
  res[CHROM == "MT" ,"CHROM":="M"]
  return(res)
}

dts_disjoint <- function(dt1,dt2){
  if(nrow(fsetdiff(dt1,dt2)) == nrow(dt1)){
    print('true')
  }else{
    print('false')
  }
}

#inputVariantsDbPath<-gunzipPath("./input/clinvar_alleles.single.b38.tsv.gz")
inputVariantsDbPath<-gunzipPath("./input/goldenStandard38.vcf.gz")
bamPath="./input/corriell_S7-ready.bam"
bamFile <- BamFile(bamPath)


#varDb <- loadClinvarData(inputVariantsDbPath) 
varDb <- loadVcfData(inputVariantsDbPath) 

bedPath <- gunzipPath("./input/corriell_S7-ready.per-base.bed.gz")
coverage <- getCoverageFromBed(bedPath)

lociAll <- getUniqueLoci(varDb)


subsetBounds <- data.table("lower"=c(1,2,6,11,21,51,101,201),
                           "upper"=c(1,5,10,20,50,100,200,9999),
                           "TD"=c(10,10,10,10,10,10,10,10))


pileups <- pileupClinVar(lociAll, 25000)
pileups <- unifyColNamesInPileups(pileups)
pileups <- addNucColumnsToPileups(pileups)

subsetBounds <- updateDepthsColumn(pileups, subsetBounds)
randomPositions <- getRandomPositionsExcludeSubset2(coverage, variantsSetToExclude=lociAll, subsetBounds)

pileupsRandom <- pileupClinVar(randomPositions, nrow(randomPositions))
pileupsRandom <- unifyColNamesInPileups(pileupsRandom)
pileupsRandom <- addNucColumnsToPileups(pileupsRandom)
dts_disjoint(pileupsRandom[,c("CHROM","POS")], pileups[,c("CHROM","POS")])

#hg38=TRUE
#pileupsRandom <- addRefFromReferenceGenome(pileupsRandom)

hist <- getSubsetDepths(pileups, subsetBounds)
histRand <- getSubsetDepths(pileupsRandom, subsetBounds)

print(hist)
print(histRand)
subsetBounds[,"TD_rand":=histRand]

#TP == lociAll, pileups: +, variant found
#TP <- sum(pileups[,altDepth>0])
#FN == lociAll: +, no variant found
#FN <- sum(pileups[,altDepth==0])


#TN == randomPositions, pileupsRandom: +, no variant found
#TN <- sum(pileupsRandom[,altDepth>0])
#FP == randomPositions: +, variant found
#FP <- sum(pileupsRandom[,altDepth==0])





#================================================================== coverage ^ ==================================================================
