library(data.table)
library(Rsamtools)
library(tictoc)
library(stringr)
library(testthat)
#library(BSgenome.Hsapiens.UCSC.hg38)
#library(BSgenome.Hsapiens.UCSC.hg19)
source("./tests.R")

clinvarCaller <- function(newPileups=TRUE,newFalsePileups=TRUE,runsAtImid=FALSE,hg38=TRUE,returnMergedDb=TRUE,calculateStats=TRUE,varDbFormat="tsv",chromosomeLengthsPath="./input/chromosomeLengths.csv",bamPath="./input/corriell_S7-ready.bam",inputVariantsDbPath="./input/clinvar_alleles.single.b38.tsv.gz",vcfHcGeneralPath="./input/corriell_S7-HCScan.vcf.gz",vcfHcFcPath="./input/corriell_S7-HCScanFC.vcf.gz",outputPileupsPath="./output/pileups.csv",outputStatsPath="./output/stats.csv",outputVarDbMergedPath="./output/varDbMerged-reduced.csv",pileupsTmp,falsePileupsTmp,notVariantPositionsTmp){

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
#pileupChrom <- function(chr,pParam){
#  chrLoci <-subset(lociSNP, chrom==as.character(chr))
#  ir <- IRanges(as.matrix(chrLoci[,"POS"])[,1],width=1)
#  res <- data.table(pileup (bamFile, indexFile, scanBamParam=ScanBamParam(which=GRanges(paste0("CHR",c(chr)), ir)), pileupParam=pParam))
#  return(res)
#}

pileupPosition <- function(chrom,pos,pParam){
  ir <- IRanges(pos,width=1)
  res <- data.table(pileup (bamFile, indexFile, scanBamParam=ScanBamParam(which=GRanges(paste0("chr",c(chrom)), ir)), pileupParam=pParam))
  return(res)
}



pileupClinVarPart <- function(loci,begin,end,pParam){
  ir <- IRanges(as.matrix(loci[begin:end,"POS"])[,1],width=1)
  res <- data.table(pileup (bamFile, indexFile, scanBamParam=ScanBamParam(which=GRanges(paste0("chr",c(as.matrix(loci[begin:end,"CHROM"])[,1])), ir)), pileupParam=pParam))
  return(res)
}


pileupClinVar <- function(loci, batchSize){
  if(batchSize>nrow(loci)) {
    print("batchSize must be less or equal to number of rows in loci")
    return(NULL)
  }
  tic(paste0("BAM, Batch size: ", batchSize, ", variants scanned: ", nrow(loci), " time"))
  n <- floor(nrow(loci)/batchSize)
  lastPileupSize <- nrow(loci)-n*batchSize
  #pp <- PileupParam(include_insertions=T, distinguish_strand=F, max_depth=1) #ignore_query=T -> count N
  pp <- PileupParam(include_insertions=T, include_deletions=T, distinguish_strand=F) #ignore_query=T -> count N
  pileups <- lapply(1:n, function(i){
    #pileups <- mclapply(1:n, function(i){
    #print(as.numeric(c(i, (i-1)*batchSize+1, i*batchSize)))
    pileup <- pileupClinVarPart(loci, (i-1)*batchSize+1, i*batchSize, pp)
    if(nrow(pileup)==0){
      blankRow <- data.table(NA,NA,NA,NA,NA)
      colnames(blankRow) <- colnames(pileup)
      pileup <- rbind(pileup,blankRow)
    }
    return(pileup)
  })
  #}, mc.cores = 4)
  pileups <- rbindlist(pileups) # data table
  if(lastPileupSize!=0){
    #print(as.numeric(c(n+1, n*batchSize+1, nrow(loci))))
    lastPileup <- pileupClinVarPart(loci, n*batchSize+1,nrow(loci),pp)
    res <- na.omit(rbind(pileups, lastPileup))
  }else{
    res <- na.omit(pileups)
  }
  toc()
  return(res)
}

unifyColNamesInPileups <- function(pileups=""){
  colnames(pileups)[colnames(pileups)=="seqnames"] <- "CHROM"
  colnames(pileups)[colnames(pileups)=="pos"] <- "POS"
  pileups$CHROM <- gsub("chr", "", pileups$CHROM)
  pileups[, `:=`(which_label = NULL)]
  return(pileups)
}



addNucColumnsToPileups <- function(pileups=""){
  pileups[, `:=`(nucA = integer(0), nucC = integer(0), nucG = integer(0), nucT = integer(0), nucIn = integer(0), nucDel = integer(0))]
  pileupsA<-subset(pileups, nucleotide=="A")
  pileupsA[,"nucA":=count]
  pileupsC<-subset(pileups, nucleotide=="C")
  pileupsC[,"nucC":=count]
  pileupsG<-subset(pileups, nucleotide=="G")
  pileupsG[,"nucG":=count]
  pileupsT<-subset(pileups, nucleotide=="T")
  pileupsT[,"nucT":=count]
  pileupsIn<-subset(pileups, nucleotide=="+")
  pileupsIn[,"nucIn":=count]
  pileupsDel<-subset(pileups, nucleotide=="-")
  pileupsDel[,"nucDel":=count]
  
  pileups <- rbind(pileupsA, pileupsC, pileupsG, pileupsT, pileupsIn, pileupsDel) # merge subsets, all blank reads are lost at this point
  rm(pileupsA,pileupsC,pileupsG,pileupsT,pileupsIn,pileupsDel)                    # remove subsets
  pileups[order(CHROM,POS)]                            						# sort by chrom, pos
  pileups[, `:=`(count = NULL, nucleotide = NULL)]         
  pileups <- pileups[, list("nucA"=sum(nucA,na.rm = TRUE), "nucC"=sum(nucC,na.rm = TRUE), "nucG"=sum(nucG,na.rm = TRUE), "nucT"=sum(nucT,na.rm = TRUE)
                            , "nucIn"=sum(nucIn,na.rm = TRUE), "nucDel"=sum(nucDel,na.rm = TRUE)), by=c("CHROM","POS")] # merge rows with matching pairs c("CHROM","POS")
  return(pileups)
}

mergePileupsWithVarDb <- function(pileups="", varDb="",allx=TRUE,ally=TRUE){
  res <- merge(varDb, pileups, by=c("CHROM","POS"), all.x=allx, all.y=ally)
  res[,c("nucA","nucC","nucG","nucT","nucIn","nucDel")][is.na(res[,c("nucA","nucC","nucG","nucT","nucIn","nucDel")])] <- 0 
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



#================================================================== Pileups ^ ==================================================================
#================================================================== varDb.merged v ==================================================================

calculateDepths <- function(variantsTable=""){
  
  # nucleotides matching REF
  variantsTable[, `:=`(refDepth = integer(0))]
  variantsTable[,c("refDepth")][is.na(variantsTable[,c("refDepth")])] <- 0
  variantsTable[substr(REF,1,1) == "A" ,"refDepth":=refDepth+nucA]
  variantsTable[substr(REF,1,1) == "C" ,"refDepth":=refDepth+nucC]
  variantsTable[substr(REF,1,1) == "G" ,"refDepth":=refDepth+nucG]
  variantsTable[substr(REF,1,1) == "T" ,"refDepth":=refDepth+nucT]
  ## variantsTable[nchar(as.character(alt))<nchar(as.character(ref)) & nchar(as.character(alt))>1, "refDepth":=refDepth+nucIn]
  ## variantsTable[nchar(as.character(alt))>nchar(as.character(ref)) & nchar(as.character(ref))>1, "refDepth":=refDepth+nucDel]
  
  # nucleotides matching ALT
  variantsTable[, `:=`(altDepth = integer(0))]
  variantsTable[,c("altDepth")][is.na(variantsTable[,c("altDepth")])] <- 0
  variantsTable[ALT == "A" ,"altDepth":=altDepth+nucA]
  variantsTable[ALT == "C" ,"altDepth":=altDepth+nucC]
  variantsTable[ALT == "G" ,"altDepth":=altDepth+nucG]
  variantsTable[ALT == "T" ,"altDepth":=altDepth+nucT]
  #variantsTable[nchar(as.character(ALT))>nchar(as.character(REF)) & substr(as.character(REF),1,1)==substr(as.character(ALT),1,1),"altDepth":=altDepth+nucIn]
  #variantsTable[(nchar(as.character(ALT))<nchar(as.character(REF)) & substr(as.character(REF),1,1)==substr(as.character(ALT),1,1)) | ALT:="-","altDepth":=altDepth+nucDel]
  
  # nucleotides not matching REF or ALT
  variantsTable[, `:=`(nullDepth = integer(0))]
  variantsTable[,c("nullDepth")][is.na(variantsTable[,c("nullDepth")])] <- 0
  variantsTable[substr(REF,1,1) != "A" & ALT != "A" ,"nullDepth":=nullDepth+nucA]
  variantsTable[substr(REF,1,1) != "C" & ALT != "C" ,"nullDepth":=nullDepth+nucC]
  variantsTable[substr(REF,1,1) != "G" & ALT != "G" ,"nullDepth":=nullDepth+nucG]
  variantsTable[substr(REF,1,1) != "T" & ALT != "T" ,"nullDepth":=nullDepth+nucT]
  #variantsTable[nchar(as.character(ALT))<=nchar(as.character(REF)),"nullDepth":=nullDepth+nucIn]
  #variantsTable[nchar(as.character(ALT))>=nchar(as.character(REF)),"nullDepth":=nullDepth+nucDel]
  
  # coverage depth
  variantsTable[,"totalDepth":=nucA+nucC+nucG+nucT]
  #variantsTable[,"depthCtrl":=altDepth+refDepth+nullDepth]
  
  # error flag
  #variantsTable[totalDepth!=depthCtrl, "errFlag":=1]
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

test_that("Does BAM file exist?", {
  expect_that(file.exists(bamPath), is_true())
})

bamFile <- BamFile(bamPath)
if(exists(paste0(bamPath, ".bai"))){
indexFile <- paste0(bamPath, ".bai")
}else{
indexFile <- indexBam(bamFile)
}
seqinfo(bamFile)

#================================================================== varDB: Choose one v ==================================================================

inputVariantsDbPath <- gunzipPath(inputVariantsDbPath)

if(varDbFormat=="vcf"){
	varDb <- loadVcfData(inputVariantsDbPath)
}
if(varDbFormat=="tsv"){
	varDb <- loadClinvarData(inputVariantsDbPath) 
}
if(varDbFormat!="tsv" & varDbFormat!="vcf"){print("varDb file must be in tsv (clinvar) or vcf format!")}


test_that("Is varDb is loaded as data.table?", {
  expect_that(is.data.table(varDb), is_true())
})

test_that("Are CHROM,POS,REF,ALT present in varDb colnames?", {
  expect_that(colnames(varDb[,c("CHROM","POS","REF","ALT")]), equals(c("CHROM","POS","REF","ALT")))
})

varDb <- normalizeDeletions(varDb)

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
lociUniqueAll <- getUniqueLoci(varDb)

#================================================================== Test sets v ==================================================================


# GATK HAPLOTYPECALLER - general scan
test_that("Does vcfHcGeneral file exist?", {
  expect_that(file.exists(vcfHcGeneralPath), is_true())
})

vcfHcGeneralPath <- gunzipPath(vcfHcGeneralPath)

vcfHcGeneral <- loadVcfData(vcfHcGeneralPath)
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
vcfHcFc <- setTenthColumnNameInVcfToGenotype(vcfHcFc)
vcfHcFc01 <- vcfHcFc[like(GENOTYPE,"0/1")]
vcfHcFc11 <- vcfHcFc[like(GENOTYPE,"1/1")]
vcfHcFc <- rbind(vcfHcFc01, vcfHcFc11)
rm(vcfHcFc01, vcfHcFc11)

test_that("Are only 0/1 and 1/1 present in GENOTYPE?", {
  expect_that(unique(substr(vcfHcFc$GENOTYPE,1,3)), equals(c("0/1","1/1")))
})


#================================================================== Test sets ^ ==================================================================


if(newPileups==TRUE){
  # ================================ A (new pileups) ================================
  pileups <- pileupClinVar(lociUniqueAll, 25000)
  pileupsTmp <- pileups
}else{
  # ================================ B (no new pileups) ================================
  if(exists("pileupsTmp")){
    pileups <- pileupsTmp
  }else{
    stop("pileupsTmp is not present! Please provide a valid pileupsTmp or set newPileups=TRUE")
  }
}  
# ================================ pileups calculations
pileups <- unifyColNamesInPileups(pileups)
pileups <- addNucColumnsToPileups(pileups)

pileupsTests(pileups)

#View(pileups)
varDb.merged <- mergePileupsWithVarDb(pileups, varDb, TRUE, TRUE)
#varDb.merged <- varDb.merged[, c("CHROM","POS","REF","ALT","nucA","nucC","nucG","nucT","nucIn","nucDel")]

#varDb.merged[, `:=`(refDepth = NULL, altDepth = NULL, nullDepth = NULL, totalDepth = NULL, depthCtrl = NULL, errFlag = NULL)]

fwrite(pileups, outputPileupsPath)

varDb.merged <- calculateDepths(varDb.merged)
#varDb.merged.f <- calculateDepths(varDb.merged.f)
View(varDb.merged)


if(calculateStats==TRUE){
  if(newFalsePileups==TRUE){
    test_that("Does chromosomeLengths file exist?", {
      expect_that(file.exists(chromosomeLengthsPath), is_true())
    })
    # ================================ A (new falsePileups) ================================
    if(nrow(lociUniqueAll)>1000000){
      notVariantPositions <- (getRandomPositionsExcludeSubset(150000,lociUniqueAll)) #149149 * 24 -> 3,6m
    }else{
      notVariantPositions <- (getRandomPositionsExcludeSubset(12500,lociUniqueAll)) #12500 * 24 -> 300k
    }
    test_that('Are notVariantPositions and lociUniqueAll disjoint?', {
      expect_that(dts_disjoint(notVariantPositions[,c("CHROM","POS")],lociUniqueAll[,c("CHROM","POS")]), prints_text('true'))
    })
    
    notVariantPositionsTmp <- notVariantPositions
    falsePileups <- pileupClinVar(notVariantPositions, 25000) # loci w/o variants - FPR
    falsePileupsTmp <- falsePileups
  }else{
    # ================================ B (no new falsePileups) ================================
    if(exists("falsePileupsTmp") & exists("notVariantPositionsTmp")){
      notVariantPositions <- notVariantPositionsTmp
      falsePileups <- falsePileupsTmp
    }else{
      stop("falsePileupsTmp and/or notVariantPositionsTmp are not present! Please provide a valid falsePileupsTmp and notVariantPositionsTmp or set newFalsePileups=TRUE")
    }
  
  }
  
  # ================================ falsePileups calculations
  falsePileups <- unifyColNamesInPileups(falsePileups)
  falsePileups <- addRefFromReferenceGenome(falsePileups)
}



# if(runsAtImid==TRUE){
	# #bamPath <- "/ibm-storage/NGS/analysis/Roche_RASGENODERM_V1/miseq reporter/44028_S1.bam"
	# getResFinal <- function(dd, varDb, lociUniqueAll){
	# res <- mclapply((1:length(dd)), function(i){ 
	  # bamPath <- dd[i]
	 
	  # print(paste0("sample nr=", i))
	  # bamFile <- BamFile(bamPath)
	  # indexFile <- paste0(bamPath, ".bai")
	  # seqinfo(bamFile)
	  # samples <- bamPath
	  
	  # pp <- PileupParam(include_insertions=T, distinguish_strand=F) #ignore_query=T -> count N

	  # pileups <- pileupClinVar(lociUniqueAll, 25000)
	  # pileups <- unifyColNamesInPileups(pileups)
	  # pileups <- addNucColumnsToPileups(pileups)
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
	# lociUniqueAll <- getUniqueLoci(varDb)

	# pgCC <- getResFinal (bamPath, varDb,lociUniqueAll)

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
	# sel <- rr [which(   rr$altDepth > 0  & ((nchar(rr$REF) ==1)  | (rr$nucDel > 0))),]
					
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
if(calculateStats==TRUE){
  FNPlusTP <- nrow(varDb.merged) # P = TP + FN, total number of positives
  ccTP <- sum(varDb.merged[,altDepth>0]) # pileups on P set only; altDepth > 0 -> variant found # altDepth == 0 -> no variant
  ccFN <- sum(varDb.merged[,altDepth==0])
  hcGeneralTP <- nrow(vcfHcGeneral[CHROM %in% varDb$CHROM & POS %in% varDb$POS])
  
  hcFcTP <- nrow(vcfHcFc[CHROM %in% varDb$CHROM & POS %in% varDb$POS]) # pileups on P set only, could be nrow(vcfHcFc)
  
  # ===== TPR ===== # sensitivity: TPR = TP/(TP+FN)
  ccTPR <- ccTP/FNPlusTP # ClinvarCaller, forcecalling
  hcGeneralTPR <- hcGeneralTP/FNPlusTP # GATK HaplotypeCaller, general scan
  hcFcTPR <- hcFcTP/FNPlusTP # GATK HaplotypeCaller, forcecalling
  
  #====================================================================== TPR ^ ======================================================================
  #====================================================================== FPR v ======================================================================
  FPPlusTN <- nrow(falsePileups) # N = FP + TN, total number of negatives
  ccFP <- sum(falsePileups[,nucleotide!=REF])
  ccTN <- sum(falsePileups[,nucleotide==REF])
  
  # ===== FPR ===== # miss rate: FPR = FP/(FP+TN)
  ccFPR <- ccFP/FPPlusTN # ClinvarCaller, forcecalling

#====================================================================== FPR ^ ======================================================================
#====================================================================== CC + HC sum v ======================================================================


  variantsCcPos<-subset(varDb.merged, altDepth>0)
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
}

fwrite(varDb.merged[,c("CHROM","POS","REF","ALT","nucA","nucC","nucG","nucT","nucIn","nucDel","refDepth","altDepth","nullDepth","totalDepth","errFlag")], outputVarDbMergedPath)



#return(list(pileupsTmp=pileupsTmp, falsePileupsTmp=falsePileupsTmp, notVariantPositionsTmp=notVariantPositionsTmp))

if(calculateStats==TRUE & returnMergedDb==FALSE){return(list(stats=stats, pileupsTmp=pileupsTmp, falsePileupsTmp=falsePileupsTmp, notVariantPositionsTmp=notVariantPositionsTmp))}
if(calculateStats==TRUE & returnMergedDb==TRUE){return(list(varDb.merged=varDb.merged, stats=stats, pileupsTmp=pileupsTmp, falsePileupsTmp=falsePileupsTmp, notVariantPositionsTmp=notVariantPositionsTmp))}

return(list(varDb.merged=varDb.merged))
}