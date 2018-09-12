library(data.table)
library(Rsamtools)
library(tictoc)
library(stringr)
library(tools)
library(BSgenome.Hsapiens.UCSC.hg38)
#library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(caret)

#TODO: optimize new(...)=TRUE/FALSE usage
#TODO: optional: multi-threading?
clinvarCaller <- function(newCoverage=TRUE,newPileups=TRUE,newFalsePileups=TRUE,newPostitionsAll=FALSE,
                          newPositionsAllFalse=TRUE,newPileupsWithStats=TRUE,newFalsePileupsWithStats=TRUE,hg38=TRUE,
                          returnMergedDb=FALSE,returnPileups=TRUE,returnPositionsAll=TRUE,returnPositionsAllFalse=TRUE,
                          returnSubsetDpDataTable=TRUE,verbose=TRUE,minBaseQuality=0,minMapq=0,fullPileupStats=TRUE,
                          fullMergedDb=TRUE,train=FALSE,phredOffset=33,
						              lowerSubsetBounds = c(1,2,6,11,21,51,101,201), 
						              upperSubsetBounds = c(1,5,10,20,50,100,200,9999),
                          chromosomeLengthsPath="./input/chromosomeLengths.csv",
                          bamPath="./input/corriell_S7-ready.bam",
                          inputVariantsDbPath="./input/clinvar_alleles.single.b38.tsv.gz",
                          modelPath="./input/model.rds",
                          listOfTrainArguments=list(method="lda",
                                                       preProc=c("center","zv","scale"),
                                                       control=trainControl(method="repeatedcv", summaryFunction=twoClassSummary, classProbs=TRUE, 
                                                                            savePredictions = TRUE, number=5, repeats=3, allowParallel = TRUE), 
                                                       formula=as.formula(c("VAR~", paste(names(dataSubset[,!c("CHROM","POS","NUC","BQ","MQ","PR","VAR","REF")]), 
                                                                                          collapse = "+")))     
                          ),
						              plotFolderPath="./output/plots/", plotTitle="Classifier evaluation",
                          outputCoveragePath="./output/corriell_S7-ready",
                          outputPositionsAllPath="./output/positionsAll.positions",
                          outputPositionsAllFalsePath="./output/positionsAllFalse.positions",
                          outputPileupsPath="./output/pileups.mpileup",
                          outputPileupsFalsePath="./output/pileupsFalse.mpileup",
                          outputPileupsWithStatsPath="./output/pileupsWithStats.csv",
                          outputFalsePileupsWithStatsPath="./output/falsePileupsWithStats.csv",
                          outputVarDbMergedPath="./output/varDbMerged-reduced.csv"){
  #TODO: make returnMergedDb true in release

if(hg38==TRUE){
  library(BSgenome.Hsapiens.UCSC.hg38)
}else{
  library(BSgenome.Hsapiens.UCSC.hg19)
}

gunzipPath <- function(string){
  out <- tryCatch(
    {
  	  if(file_ext(string)=="gz"){
  	  	res <- paste0("gunzip -c ",string)
  	  }else{
  	  	res <- string
  	  }
    },
    error=function(cond) {
	  message("Something went wrong. Function: gunzipPath")
	  message("This function adds 'gunzip -c ' before paths to gzipped files (string). Gzip has to be installed - otherwise the gunzip command will not execute properly. ")
	  message("Returns: gunzip command + path string")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		return(res)
    })
  return(out)
}

#gunzipPath
loadBamFile <- function(bamPath=""){
  out <- tryCatch(
    {
      if(!file.exists(bamPath)){
        stop("bamPath file is not present! Please provide a valid bamPath!")
      }
		  bamFile <- BamFile(bamPath)
    },
    error=function(cond) {
	  message("Something went wrong. Function: loadBamFile")
	  message("This function loads a data from a BAM file at specified bamPath.")
	  message("Returns: BAM data")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		if(verbose==TRUE){print("BAM file loaded.")}
		return(bamFile)
    })
  return(out)
}

#loadBamFile
loadIndexFile <- function(bamPath=""){
  out <- tryCatch(
    {
		if(exists(paste0(bamPath, ".bai"))){
		  indexFile <- paste0(bamPath, ".bai")
		  if(verbose==TRUE){print("Index (BAI) file loaded.")}
		}else{
		  if(!file.exists(bamPath)){
		    stop("bamPath file is not present! Please provide a valid bamPath!")
		  }
		  bamFile <- BamFile(bamPath)
		  indexFile <- indexBam(bamFile)[[1]]
		  if(verbose==TRUE){print("Index (BAI) file created and loaded.")}
		}
    },
    error=function(cond) {
	  message("Something went wrong. Function: loadIndexFile")
	  message("This function loads a BAM index file (BAI) from specified bamPath. If the file does not exist, this function also creates it.")
	  message("Returns: path to BAM index file (BAI)")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		return(indexFile)
    })
  return(out)
}

makeBed <- function(coveragePath="", bamPath=""){
  out <- tryCatch(
    {
		if(!file.exists(bamPath)){
		  stop("bamPath file is not present! Please provide a valid bamPath!")
		}
	  cmd <- paste0("mosdepth '",coveragePath,"' ",bamPath)
	  system(cmd)
    },
    error=function(cond) {
	  message("Something went wrong. Function: makeBed")
	  message("This function requires mosdepth installed to work properly. It creates a coverage (BED) file at specified coveragePath for a specified BAM file (bamPath). Path to the coverage file is returned.")
	  message("Returns: path to the gzipped coverage (BED) file")
    message("Original error message:\n")
    message(cond)
      return(NA)
    },
    finally={
		  return(paste0(coveragePath,".per-base.bed.gz"))
    })
  return(out)
}

filterChroms <- function(input){
  out <- tryCatch(
    {
	  res <- input[which(input$CHROM %in% c(1:22, "X", "Y", "M")),]
    },
    error=function(cond) {
	  message("Something went wrong. Function: filterChroms")
	  message("This function trims records in the input data.table. Only records with CHROM 1:22, X, Y or M are left.")
	  message("Returns: trimmed data.table")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		  return(res)
    })
  return(out)
}


#filterChroms
getCoverageFromBed <- function(bedPath=""){
    out <- tryCatch(
    {
		if(!file.exists(bedPath)){
			stop("bedPath file is not present! Please provide a valid bedPath!")
		}
		coverage <- fread(gunzipPath(bedPath))
		colnames(coverage) <- c("CHROM","start","stop","readsCount")
		coverage <- subset(coverage, readsCount > 0)
		coverage$CHROM <- gsub("chr", "", coverage$CHROM)
		coverage <- filterChroms(coverage)
		coverage[,"start":=start+1L]
		#coverage[,"length":=stop-start]
    },
    error=function(cond) {
  	  message("Something went wrong. Function: getCoverageFromBed")
  	  message("This function reads the coverage from the specified BED file.")
  	  message("Returns: data.table containing coverage data from specified BED file (with headers)")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		return(coverage)
    })
  return(out)
}

#TODO: try-catch
#gunzipPath, makeBed, getCoverageFromBed, filterChroms
getCoverage <- function(outputCoveragePath="", bamPath="", new=FALSE){
    out <- tryCatch(
		{
			if(new==TRUE | !file.exists(paste0(outputCoveragePath,".per-base.bed.gz"))){
				bedPath <- makeBed(outputCoveragePath=outputCoveragePath, bamPath=bamPath)
				if(verbose==TRUE){print("Coverage file (BED) created.")}
			}else{
			  bedPath <- paste0(outputCoveragePath,".per-base.bed.gz")
			}
			if(!file.exists(bedPath)){
			  stop("BED file is not present!")
			}
			coverage <- getCoverageFromBed(bedPath)
		},
		error=function(cond) {
		  message("Something went wrong. Function: getCoverage")
		  message("This function reads the coverage for the specified BAM file (bamPath) from its coverage file (outputCoveragePath). If the coverage file does not exist or the 'new' atrribute is set to TRUE, the file is created beforehand.")
		  message("Returns: data.table containing coverage data for the specified BAM file")
		  message("Original error message:\n")
		  message(cond)
		  return(NA)
		},
		finally={
			if(verbose==TRUE){print("Coverage file (BED) loaded.")}
			return(coverage)
    })
  return(out)
}

getFileExtension <- function(path="", removeLastExt = FALSE){
    out <- tryCatch(
    {
	  if(removeLastExt==TRUE & str_count(str_remove(path,"[.]{1,}"), "\\.")>1){
		path <- file_path_sans_ext(path)
	  }
	  res <- file_ext(path)
    },
    error=function(cond) {
	  message("Something went wrong. Function: getFileExtension")
	  message("This function checks the extension of the file (path). If removeLastExt argument is set to TRUE, it omits the last segment of the extension (eg. input 'file.vcf.gz' -> output 'vcf')")
	  message("Returns: extension of the specified file (without the leading dot)")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		return(res)
    })
  return(out)
}

#TODO: merge with loadTsvData?
#gunzipPath
loadVcfData <- function(vcfPath=""){
   out <- tryCatch(
    {
  	  res <- fread(gunzipPath(vcfPath), skip = "CHROM")
  	  colnames(res)[colnames(res)=="#CHROM"] <- "CHROM"
  	  res$CHROM <- gsub("chr", "", res$CHROM)
  	  res <- res[!like(ALT,",")] # remove multiallelic sites
  	  res <- subset(unique(res, by=c("CHROM","POS")))
  	  res[CHROM == "MT" ,"CHROM":="M"]
    },
    error=function(cond) {
	  message("Something went wrong. Function: loadVcfFile")
	  message("This function loads a data from a VCF file at specified vcfPath.")
	  message("Returns: VCF data")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		return(res)
    })
  return(out)
}

#gunzipPath
loadTsvData <- function(tsvPath=""){
   out <- tryCatch(
    {
	  res <- fread(gunzipPath(tsvPath))
	  colnames(res)[colnames(res)=="chrom"] <- "CHROM"
	  colnames(res)[colnames(res)=="pos"] <- "POS"
	  colnames(res)[colnames(res)=="ref"] <- "REF"
	  colnames(res)[colnames(res)=="alt"] <- "ALT"
	  res[CHROM == "MT" ,"CHROM":="M"]
	  #normalizeDeletions(res)
    },
    error=function(cond) {
	  message("Something went wrong. Function: loadTsvFile")
	  message("This function loads a data from a TSV file at specified tsvPath.")
	  message("Returns: TSV data")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		return(res)
    })
  return(out)
}

#TODO: check if positions arent already normalized, make insertions=="+"?
normalizeDeletions <- function(vcfData=""){
  out <- tryCatch(
    {
	  res <- copy(vcfData)
	  res[nchar(as.character(ALT))==1 & nchar(as.character(REF))>1 & substr(as.character(REF),1,1)==as.character(ALT),`:=`("POS"=POS+1L, "ALT"="-", "REF"=substring(as.character(REF),2))]
    },
    error=function(cond) {
	  message("Something went wrong. Function: normalizeDeletions")
	  message("This function normalizes the deletion format in data.table. The deletion POS should be incremented. Otherwise pileup would not find it. ALT is then turned into '-'.")
	  message("Returns: input data.table with incremented deletion positions and ALT='-'")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		return(res)
    })
  return(out)
}

#gunzipPath, getFileExtension, loadVcfData, loadTsvData, normalizeDeletions
loadVarDb <- function(inputVariantsDb=NULL){
  out <- tryCatch(
    {
		if(typeof(inputVariantsDb)=="character"){
			if(file.exists(inputVariantsDb)){
					varDbFormat <- getFileExtension(inputVariantsDb, removeLastExt = TRUE)					
					if(varDbFormat=="vcf"){
						varDb <- loadVcfData(inputVariantsDb)
						if(verbose==TRUE){print("varDb VCF data loaded.")}
					}
					if(varDbFormat=="tsv"){
						varDb <- loadTsvData(inputVariantsDb) 
						if(verbose==TRUE){print("varDb TSV data loaded.")}
					}
					if(varDbFormat!="tsv" & varDbFormat!="vcf"){
						stop("varDb file must be in tsv or vcf format!")
					}
			}else{
			  stop("inputVariantsDb file is not present! Please provide a valid path!")
			}
		}else{
			varDb <- inputVariantsDb
		}
		if(is.data.table(varDb)!=TRUE){
		  stop("varDb is not a data.table!")
		}
		varDb <- normalizeDeletions(varDb)
		if(verbose==TRUE){print("Deletions normalized.")}
		#varDb <- setTenthColumnNameInVcfToGenotype(vcfTable = varDb)
	    },
    error=function(cond) {
	  message("Something went wrong. Function: loadVarDb")
	  message("This function reads variant database from external file (inputVariantsDb). Database can be in VCF or TSV format.")
	  message("Returns: data.table with variant database")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		return(varDb)
    })
  return(out)
}

getUniqueLoci <- function(vcfData=NULL, variationType="", saveNucleotideData=FALSE){
  out <- tryCatch(
    {
		  res <- vcfData[,c("CHROM","POS","REF","ALT")]
		  #res[,"lengthDiff":=nchar(as.character(REF))-nchar(as.character(ALT))] # snp: 0, in: <0, del: >0
		  if(variationType=="in"){res <- res[(nchar(as.character(ALT))>1 & nchar(as.character(REF))==1 & as.character(REF)==substr(as.character(ALT),1,1)) | ALT=="+"]}
		  if(variationType=="del"){res <- res[(nchar(as.character(ALT))==1 & nchar(as.character(REF))>1 & substr(as.character(REF),1,1)==as.character(ALT)) | ALT=="-"]}
		  if(variationType=="snp"){res <- res[nchar(as.character(ALT))==1 & nchar(as.character(REF))==1]}
		  if(variationType=="notMatching"){res <- res[nchar(as.character(ALT))>1 & nchar(as.character(REF))>1]}
		  if(saveNucleotideData==TRUE){
			res <- subset(unique(res[,c("CHROM","POS","REF","ALT")], by=c("CHROM","POS")))
		  }else{
			res <- subset(unique(res[,c("CHROM","POS")], by=c("CHROM","POS")))
		  }
  	},
    error=function(cond) {
	  message("Something went wrong. Function: getUniqueLoci")
	  message("This function extracts unique loci from provided data.table (vcfData). If variationType is specified, the result will contain only given variation type (snp, in, del, notMatching). If saveNucleotideData is set to TRUE, REF and ALT columns will be included in the result.")
	  message("Returns: data.table with CHROM + POS of unique loci in the dataset. REF and ALT can be also included.")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		return(res)
    })
  return(out)
}

#TODO: check if works
#getUniqueLoci
getPositions <- function(positionsAll="", varDb=NULL, new=FALSE){
  out <- tryCatch(
	{
		if(typeof(positionsAll)=="character"){
			if(new==TRUE | !file.exists(outputPositionsAllPath)){
				res <- getUniqueLoci(varDb)
				positionsAll2 <- copy(res)
				positionsAll2[,"CHROM":=paste0("chr",CHROM)]
				fwrite(positionsAll2, positionsAll, sep= " ")
				rm(positionsAll2)
				if(verbose==TRUE){print("positionsAll exported to external file.")}
			}else{
				if(file.exists(positionsAll)){
				  res <- fread(positionsAll)
				  res$CHROM <- gsub("chr", "", res$CHROM)
				  if(verbose==TRUE){print("positionsAll restored from outputPositionsAllPath file.")}
				}else{
				  stop("outputPositionsAll file is not present! Please provide a valid path!")
				}
			}
		}else{
			if(new==TRUE){
				stop("outputPositionsAllPath has to be a string if you want to create a new file!")
			}else{
				res <- positionsAll
			}
		}
	},
    error=function(cond) {
	  message("Something went wrong. Function: getPositions")
	  message("This function extracts positions from external file (outputPositionsAllPath) or data.table with CHROM and POS (varDb) columns. In the second case, it also writes them to the external file (outputPositionsAllPath). If argument new is set to TRUE, new file will be created from varDb, regardless of its previous existence.")
	  message("Returns: data.table with CHROM + POS extracted from varDb. Also writes an external file.")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		return(res)
    })
  return(out) 
}

addRefFromReferenceGenome <- function(dt="", hg38=TRUE){
  out <- tryCatch(
	{
			ir <- IRanges(as.matrix(dt[,"POS"])[,1],width=1)
			gr <- GRanges(paste0("chr",c(as.matrix(dt[,"CHROM"])[,1])), ir)
			if(hg38==TRUE){
			grseq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)
			}else{
			grseq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, gr)
			}
			grseq <- data.table(as.character(grseq))
			colnames(grseq)[colnames(grseq)=="V1"] <- "REF"
			res <- dt[, "REF" := grseq[["REF"]]]
  
  	},
    error=function(cond) {
	  message("Something went wrong. Function: addRefFromReferenceGenome")
	  message("This function checks REF nucleotides in the reference genome and adds them to the data.table dt, basing on CHROM+POS pairs in existing columns of dt. Input data.table must contain CHROM and POS columns. For reference genome hg38 - set hg38 to TRUE. Otherwise hg19 will be used.")
	  message("Returns: data.table with additional REF column.")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		return(res)
    })
  return(out) 
}

#TODO: offset parameter
phredQualityDecode <- function(n=0, offset=33){
  out <- tryCatch(
	{
		res <- utf8ToInt(as.character(n))-offset
  	},
    error=function(cond) {
	  message("Something went wrong. Function: phredQualityDecode")
	  message("This function transforms phred-scaled quality utf8 symbol into numeric value.")
	  message("Returns: a number.")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		return(res)
    })
  return(out) 
}

phredQualityEncode <- function(n, offset=33){
  out <- tryCatch(
	{
		res <- intToUtf8(as.integer(n)+offset)
  	},
    error=function(cond) {
	  message("Something went wrong. Function: phredQualityEncode")
	  message("This function transforms numeric value into phred-scaled quality utf8 symbol.")
	  message("Returns: an utf8 symbol.")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		return(res)
    })
  return(out) 
}

#TODO: merge with doPileups?
#addRefFromReferenceGenome, phredQualityDecode, phredQualityEncode
pileupVariantDatabase <- function(lociPath="",pileupPath="",bamPath="",minBaseQuality=0,minMapq=0,phredOffset=33,newPileups=TRUE,normalizeData=TRUE,hg38=TRUE){
  out <- tryCatch(
	{
	  if(newPileups==TRUE | !file.exists(pileupPath)){
		cmd <- paste0("samtools mpileup -A -d 0 -l ",lociPath," -q ",minMapq," -Q ",minBaseQuality," -O -s ",bamPath," > ",pileupPath)
		system(cmd)
	  }
	  pileups <- fread(pileupPath, quote="")
	  colnames(pileups) <- c("CHROM","POS","REF","DP","NUC","BQ","MQ","PR")
	  pileups[,"REF":=NULL]
	  pileups$CHROM <- gsub("chr", "", pileups$CHROM)
	  pileups <- addRefFromReferenceGenome(pileups, hg38=hg38)
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

		# decode BQ, MQ
		pileups[,"BQ":=paste(phredQualityDecode(BQ, offset=phredOffset),sep=",",collapse = ","), by=1:nrow(pileups)]
		pileups[,"MQ":=paste(phredQualityDecode(MQ, offset=phredOffset),sep=",",collapse = ","), by=1:nrow(pileups)]
		pileups[,`:=`(NUC=strsplit(NUC,""),BQ=strsplit(BQ,","),MQ=strsplit(MQ,","),PR=strsplit(PR,","))]
	  }
  
    },
    error=function(cond) {
	  message("Something went wrong. Function: pileupVariantDatabase")
	  message("This function required samtools installed to work properly. It uses samtools to create pileup file (pileupPath) from BAM data (bamPath). Positions are taken from variant database - lociPath. Afterwards the new file is imported into a data.table. If the file is already present, the creation step can be omitted by making newPileups=FALSE. Minimal base and mapping qualities can be specified - only higher quality reads will be used in pileup. Hg38 should be set to TRUE if the sample matches human reference genome hg38. Otherwise, hg19 will be used. In order to make the output processible by this tool, normalizeData should be set to TRUE.
	  ")
	  message("Returns: data.table with (normalized) pileup outcome.")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		return(pileups)
    })
  return(out) 
}

calculateNucStats <- function(chrom="", pos="",ref="", nuc="", countsOnly=FALSE, bq="", mq="", pr=""){
  out <- tryCatch(
	{
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
    },
    error=function(cond) {
	  message("Something went wrong. Function: calculateNucStats")
	  message("This function calculates simple nucleotide statistics for one set of CHROM, POS, RES, NUC, BQ, MQ and PR values. These values should be later used in variant calssification.")
	  message("Returns: 1-row data.table with nucleotide statistics ()")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		return(summaryDt)
    })
  return(out)	  
}

#TODO: try-catch
#calculateNucStats
addStatsColumnsToPileup <- function(pileups=NULL, outputPileupsWithStatsPath="", countsOnly=FALSE, rmCols=FALSE, newPileupsWithStats=TRUE){
  out <- tryCatch(
	{
		  #TODO: READ pileups with stats
		  if(newPileupsWithStats==TRUE | !file.exists(outputPileupsWithStatsPath)){
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
      },
    error=function(cond) {
	  message("Something went wrong. Function: addStatsColumnsToPileup")
	  message("This function creates new columns in pileup data.table. These values are later used to classify variants.")
	  message("Returns: pileup data.table extended with numerous columns with simple pileup statistics")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		return(res)
    })
  return(out)	
}

#TODO: input from variable
#TODO: check if works
#pileupVariantDatabase, calculateNucStats, addStatsColumnsToPileup
doPileups <- function(outputPositionsPath="",outputPileupsPath="",outputPileupsWithStatsPath="",bamPath="",minBaseQuality=0,minMapq=0,phredOffset=33,newPileups=TRUE,newPileupsWithStats=TRUE){
  out <- tryCatch(
	{
		if(newPileups==TRUE & newPileupsWithStats!=TRUE){
		  pileups <- pileupVariantDatabase(outputPositionsPath,outputPileupsPath,bamPath,phredOffset=phredOffset,minBaseQuality,minMapq,newPileups=TRUE)
		  if(verbose==TRUE){print("pileups done.")}
		}else{
		  if(file.exists(outputPileupsPath)){
			pileups <- pileupVariantDatabase(outputPositionsPath,outputPileupsPath,bamPath,phredOffset=phredOffset,minBaseQuality,minMapq,newPileups=FALSE)
			if(verbose==TRUE){print("pileups restored from outputPileupsPath file.")}
		  }else{
			stop("outputPileupsPath file is not present! Please provide a valid outputPileupsPath or set newPileups=TRUE")
		  }
		}
		pileups <- addStatsColumnsToPileup(pileups=pileups, outputPileupsWithStatsPath=outputPileupsWithStatsPath, countsOnly = !(fullPileupStats), rmCols = FALSE, newPileupsWithStats=newPileupsWithStats) #TODO: make rmCols true in release
		if(verbose==TRUE){print("Pileups: New columns created.")}	
	},
    error=function(cond) {
	  message("Something went wrong. Function: doPileups")
	  message("This function wraps together pileupVariantDatabase and addStatsColumnsToPileup.")
	  message("Returns: pileup data.table extended with numerous columns with simple pileup statistics")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		return(pileups)
    })
  return(out)
}


#======================================================train

initializeSubsetDpDataTable <- function(lowerBounds = c(1,2,6,11,21,51,101,201), upperBounds = c(1,5,10,20,50,100,200,9999)){
  out <- tryCatch(
	{
		if(length(lowerBounds)!=length(upperBounds)){
			stop("lowerBounds and upperBounds have different sizes!")
		}
		if(!all(lowerBounds<=upperBounds)){
			stop("All values in lowerBounds have to be lower than corresponding values in upperBounds!")
		}
		subsetDpDataTable <- data.table("lower"=lowerBounds,
							   "upper"=upperBounds,
							   "DP"=c(matrix(0,1,length(lowerBounds))),
							   "DP_rand"=c(matrix(0,1,length(lowerBounds)))
							   )		
	},
    error=function(cond) {
	  message("Something went wrong. Function: initializeSubsetDpDataTable")
	  message("This function creates initial subsetDpDataTable data.table, which contains lower and upper values of DP ranges used in calculations. Created data.table contains 2 additional columns - number of positions with specified DP's will be stored there.")
	  message("Returns: data.table with DP values and 2 blank columns to store position counts")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		return(subsetDpDataTable)
    })
  return(out)
}

getDP <- function(pileupsWithNucCols=NULL, removeNucColumns=FALSE){
  out <- tryCatch(
	{
	  res <- copy(pileupsWithNucCols)
	  res[,"DP":=countA+countC+countG+countT+countIn+countDel+countN]
	  if(removeNucColumns==TRUE){
		res[, `:=`(countA = NULL, countC = NULL, countG = NULL, countT = NULL, countIn = NULL, countDel = NULL, countN= NULL)]
	  }
  	},
    error=function(cond) {
	  message("Something went wrong. Function: getDP")
	  message("This function calculates DP from nucleotide counts in 1 row of data.table and adds it to the input data.table. If removeNucColumns is set to TRUE, columns with nucleotide counts are later removed.")
	  message("Returns: data.table with additional DP column. May remove the nucleotide count columns on demand.")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		return(res)
    })
  return(out)
}

#getDP
getSubsetDepths <- function(subsetBounds=NULL, pileupsWithNucCols=NULL){
  out <- tryCatch(
	{
	  res <- copy(subsetBounds)
	  totalDepths <- getDP(pileupsWithNucCols=pileupsWithNucCols, removeNucColumns=TRUE)
	  res <- sapply(1:nrow(subsetBounds), function(i){
		return(nrow(subset(totalDepths, (DP >= subsetBounds[i,lower] & DP <= subsetBounds[i,upper]))))
	  })
  	},
    error=function(cond) {
	  message("Something went wrong. Function: getSubsetDepths")
	  message("This function calculates number of reads belonging to each subset (with specified DP bounds).")
	  message("Returns: data.table with single column, that contains numbers of reads belonging to each subset.")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		return(res)
    })
  return(out)
}

#initializeSubsetDpDataTable, getSubsetDepths, getDP
updateDepthsColumn <- function(subsetBounds=NULL, updateDp=TRUE, updateDpRand=TRUE, pileupsWithNucCols=NULL, falsePileupsWithNucCols=NULL){
  out <- tryCatch(
	{
	  if(!exists("subsetBounds")){
		subsetBounds <- initializeSubsetDpDataTable()
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
  	},
    error=function(cond) {
	  message("Something went wrong. Function: getSubsetDepths")
	  message("This function uses getSubsetDepths to update Dp and DpRand columns in subset bounds data.table. Data.tables with pileups sholud be provided as arguments. Function can update both Dp and DpRand or just one of them.")
	  message("Returns: updated subset bounds data.table with number of reads belonging to each subset.")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		return(res)
    })
  return(out)
}

getCoverageSubsetByReadsCount <- function(coverage=NULL, min=1, max=1){
  out <- tryCatch(
	{
	  if(min<=max){
		res <- subset(coverage, (readsCount >= min & readsCount <= max))
	  }else{
		stop("max must be greater or equal to min!")
	  }  
    },
    error=function(cond) {
	  message("Something went wrong. Function: getCoverageSubsetByReadsCount")
	  message("This function creates a subset of rows from coverage (BED) data.table, basing on positions (rows) readsCount.")
	  message("Returns: data.table containing subset of input data.table.")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		return(res)
    })
  return(out)
}

getPositionsFromCoverageSubset <- function(coverageSubset=NULL,subsetMaxLength=200000){
  out <- tryCatch(
	{
	  if(nrow(coverageSubset)>subsetMaxLength){coverageSubset <- coverageSubset[sample(.N, subsetMaxLength, replace=F)]}
	  res <- rbindlist(
		  lapply(c(1:22, "X", "Y", "M"), function(i){
		  chromSubset <- subset(coverageSubset, CHROM==i)
		  if(nrow(chromSubset)){
			chromPositions <- rbindlist(
			  lapply(min(chromSubset[,readsCount]):max(chromSubset[,readsCount]), function(j){
			  countSubset <- subset(chromSubset, readsCount==j)
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
    },
    error=function(cond) {
	  message("Something went wrong. Function: getPositionsFromCoverageSubset")
	  message("This function creates a data.table containing all of the positions from coverage (BED) file. subsetMaxLength parameter specifies maximum numbr of positions for function to check. ")
	  message("Returns: data.table with CHROM, POS and readsCount columns.")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		return(res)
    })
  return(out)
}

getRandomPositions <- function(listOfPositions=NULL, randomSubsetSize=10000){
  out <- tryCatch(
	{
	  if(nrow(listOfPositions)>=randomSubsetSize){
		res <- listOfPositions[sample(.N, randomSubsetSize, replace=F)]
	  }else{
		res <- listOfPositions  
	  }
    },
    error=function(cond) {
	  message("Something went wrong. Function: getRandomPositions")
	  message("This function takes a random sample (size: randomSubsetSize) from specified list of positions (listOfPositions).")
	  message("Returns: data.table, sample taken from another data.table with list of positions.")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		return(res)
    })
  return(out)
}

#getCoverageSubsetByReadsCount, getPositionsFromCoverageSubset, getRandomPositions
getRandomPositionsFromCoverage <- function(coverage=NULL,subsetBounds=NULL,subsetMaxLength=200000){
  out <- tryCatch(
	{
	  res <- lapply(1:nrow(subsetBounds), function(i){
		coverageSubset <- getCoverageSubsetByReadsCount(coverage,min=subsetBounds[i,lower],max=subsetBounds[i,upper])
		listOfPositions <- getPositionsFromCoverageSubset(coverageSubset,subsetMaxLength)
		randomPositions <- getRandomPositions(listOfPositions, subsetBounds[i,DP])
		return(randomPositions)
	  })
	  res <- rbindlist(res)
    },
    error=function(cond) {
	  message("Something went wrong. Function: getRandomPositionsFromCoverage")
	  message("This function wraps together coverageSubsetand, listOfPositions and randomPositions. It generates a set of random positions, with coverage similar to the one specified in subsetBounds.")
	  message("Returns: a set of random positions with coverage similar to the one specified in DP column of (subsetBounds), taken from coverage (BED) file.")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		return(res)
    })
  return(out)
}

#getRandomPositionsFromCoverage, getCoverageSubsetByReadsCount, getPositionsFromCoverageSubset, getRandomPositions
getRandomPositionsExcludeSubset <- function(coverage=NULL,variantsSetToExclude=NULL,subsetBounds=NULL,returnReadsCount=FALSE){ 
  out <- tryCatch(
	{
	  posCount <- sum(subsetBounds[,"DP"])
	  if(posCount<=0){
		stop("posCount must be greater than 0!")
	  }
	  randomSubsetSize <- ceiling(posCount/nrow(subsetBounds))
	  subsetMaxLength <- randomSubsetSize
	  tmpPosCount <- posCount
	  res <- getRandomPositionsFromCoverage(coverage,subsetBounds,subsetMaxLength)
	  repeat{
		res <- subset(unique(res, by=c("CHROM","POS")))
		variantsSetToExclude[,"CHROM":=as.character(CHROM)]
		res <- res[!variantsSetToExclude, on=.(CHROM,POS)]
		if(nrow(res)>=posCount){
		  break
		}
		tmpPosCount <- posCount - nrow(res)
		randomSubsetSize <- ceiling(tmpPosCount/nrow(subsetBounds))
		subsetMaxLength <- randomSubsetSize*3
		res <- rbind(res,getRandomPositionsFromCoverage(coverage,subsetBounds,subsetMaxLength))
	  }
	  if(returnReadsCount==TRUE){
		res <- res[sample(.N, posCount, replace=F),c("CHROM","POS","readsCount")]
	  }else{
		res <- res[sample(.N, posCount, replace=F),c("CHROM","POS")]
	  }
    },
    error=function(cond) {
	  message("Something went wrong. Function: getRandomPositionsExcludeSubset")
	  message("This function generates a set of random positions, with coverage similar to the one specified in subsetBounds. Positions from variantsSetToExclude cannot appear in the generated set. If returnReadsCount is TRUE, readsCount column will be included in the outcome.")
	  message("Returns: a set (data.table) of random positions, with coverage similar to the one specified in subsetBounds and not included in variantsSetToExclude ")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		return(res)
    })
  return(out)
}

#TODO: try-catch
dts_disjoint <- function(dt1,dt2){
    res <- (nrow(fsetdiff(dt1,dt2)) == nrow(dt1))
  return(res)
}

#TODO: check if works
#getRandomPositionsExcludeSubset, getRandomPositionsFromCoverage, getCoverageSubsetByReadsCount, getPositionsFromCoverageSubset, getRandomPositions, dts_disjoint
getFalsePositions <- function(positionsAllFalse=NULL, new=TRUE, coverage=NULL, variantsSetToExclude=NULL, subsetDpDataTable=NULL){
  out <- tryCatch(
	{
		if(typeof(positionsAllFalse)=="character"){
			if(new==TRUE | !file.exists(positionsAllFalse)){
			  res <- getRandomPositionsExcludeSubset(coverage, variantsSetToExclude=variantsSetToExclude, subsetDpDataTable,returnReadsCount=FALSE)
			  if(verbose==TRUE){print("Random positions generated.")}
			  positionsAllFalse2 <- copy(res)
			  positionsAllFalse2[,"CHROM":=paste0("chr",CHROM)]
			  fwrite(positionsAllFalse2, positionsAllFalse, sep= " ")
			  rm(positionsAllFalse2)
			  if(verbose==TRUE){print("positionsAllFalse exported to external file.")}
			  
			  #TODO: input arguments
			  if(!dts_disjoint(res[,c("CHROM","POS")],variantsSetToExclude[,c("CHROM","POS")])){
			    stop("positionsAllFalse and variantsSetToExclude should be disjoint!")
			  }
			  #falsePileups <- pileupVariantDatabase(positionsAllFalse, 25000) # loci w/o variants - FPR
			}else{
			  if(file.exists(positionsAllFalse)){
				res <- fread(positionsAllFalse)
				res$CHROM <- gsub("chr", "", res$CHROM)
				if(verbose==TRUE){print("positionsAllFalse restored from positionsAllFalse file.")}
			  }else{
				stop("positionsAllFalse file is not present! Please provide a valid positionsAllFalse path or set newPositionsAllFalse=TRUE")
			  }
			}
		}else{
			if(new==TRUE){
				stop("positionsAllFalse has to be a string if you want to create a new file!")
			}else{
				res <- positionsAllFalse
			}
		}
	},
    error=function(cond) {
	  message("Something went wrong. Function: getFalsePositions")
	  message("This function provides a complete pipeline needed to construct false (without variants) positions dataset. Writes the result to external file if specified or if the file does not exist.")
	  message("Returns: false positions. Writes the result to external file.")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		return(positionsAllFalse)
    })
  return(out)
}

#addRefFromReferenceGenome
prepareTrainingData <- function(pileups=NULL, falsePileups=NULL, hg38=TRUE){
  out <- tryCatch(
	{
	  pileups2 <- copy(pileups)
	  pileupsRandom2 <- copy(falsePileups)
	  pileups2[,"VAR":="variant"]
	  pileupsRandom2[,"VAR":="nonVariant"]
	  pileupsAll <- rbind(pileups2,pileupsRandom2)
	  pileupsAll <- pileupsAll[,"VAR":=as.factor(VAR)]
	  pileupsAll <- addRefFromReferenceGenome(pileupsAll,hg38=hg38)
	  pileupsAllShuffled <- pileupsAll[sample(nrow(pileupsAll)),]
	  rm(pileups2,pileupsRandom2,pileupsAll)
	},
    error=function(cond) {
	  message("Something went wrong. Function: prepareTrainingData")
	  message("This function prepares the training data from pileups and false pileups. hg38 should be set if the reference genome is hg38. Otherwise, hg19 will be used.")
	  message("Returns: data.table with set of balanced training data.")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		return(pileupsAllShuffled)
    })
  return(out)
}

trainModel <- function(modelPath="", data=NULL, listOfTrainArguments=NULL){
  out <- tryCatch(
	{
	  seed <- 7
	  set.seed(seed)
	  model <- train(listOfTrainArguments$formula, data=data, method=listOfTrainArguments$method, 
					 metric=listOfTrainArguments$metric, preProc=listOfTrainArguments$preProc, 
					 trControl=listOfTrainArguments$control)
	  saveRDS(model, modelPath)
	},
    error=function(cond) {
	  message("Something went wrong. Function: trainModel")
	  message("This function trains the model and saves it to the external file, under specified modelPath. Training data should be provided via data argument. Training arguments sholud be provided via listOfTrainArguments. The argument list should match the caret library train arguments format: ' train(listOfTrainArguments$formula, data=data, method=listOfTrainArguments$method, metric=listOfTrainArguments$metric, preProc=listOfTrainArguments$preProc, trControl=listOfTrainArguments$control)' ")
	  message("Returns: trained model. Also saves it to external file.")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
		return(model)
    })
  return(out)
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
	  message("Something went wrong. Function: classifyNaive")
	  message("This function classifies the data with a naive assumption, that if 1 or more nucleotides does not match the REF, there is a variant under given position.")
	  message("Returns: data.table with the result of naive classification.")
      message("Original error message:\n")
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
          result <- predict(externalModel, newdata = newData, type="prob")
        }else{
          stop("modelPath file is not present! Please provide a valid classification model!")
        }
      }else{
        result <- predict(model, newdata = newData, type="prob")
      }
      result <- data.table(result)
      result[variant>=nonVariant,result:="variant"]
      result[variant<nonVariant,result:="nonVariant"]
      if("VAR" %in% colnames(newData)){
        result <- cbind(result,newData[,c("VAR")],newData[,!c("VAR")])
      }else{
        result <- cbind(result,newData)
      } 
    },
    error=function(cond) {
	  message("Something went wrong. Function: classify")
	  message("This function classifies the data (newData). Classification model should be specified via 'model' argument (either path or R variable).")
	  message("Returns: data.table with the result of classification.")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
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

  varDb[,"CHROM":=as.character(CHROM)]

  res <- merge(varDb, tmp[,c("CHROM","POS","RefCount","AllNonRefCount","MaxNonRefCount","InDelMaxLength",
                             "countA","countC","countG","countT","countIn","countDel","countN","result")], by=c("CHROM","POS"), all.x=allx, all.y=ally)
  res$result[is.na(res$result)] <- "nonVariant"  #WARN: potentially breaks confusion matrix
  rm(tmp)
  #TODO: maybe shorter? 
  #res[,c("countA","countC","countG","countT","countIn","countDel","countN")][is.na(res[,c("countA","countC","countG","countT","countIn","countDel","countN")])] <- 0
  res[,][is.na(res[,])] <- 0 
  
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
  
  if("varCount" %in% colnames(variantsTable)){variantsTable[, `:=`(varCount = NULL)]}
  variantsTable[, `:=`(varCount = integer(0))]
  variantsTable[,c("varCount")][is.na(variantsTable[,c("varCount")])] <- 0
  #AD + indele
  variantsTable[AD>0,"varCount":=varCount+AD]
  
  variantsTable[(nchar(as.character(ALT))>nchar(as.character(REF)) & substr(as.character(REF),1,1)==substr(as.character(ALT),1,1)) | ALT=="-", "varCount":=varCount+countIn]
  variantsTable[(nchar(as.character(ALT))<nchar(as.character(REF)) & substr(as.character(REF),1,1)==substr(as.character(ALT),1,1)) | ALT=="-","varCount":=varCount+countDel]
  #variantsTable[substr(ALT,1,1) == "N" ,"varCount":=varCount+countN] # ALT is never N 
  
  
  
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
      FPRdiffProc = -round(100*((classifierDataStats[[i,"FPR"]] - naiveDataStats[[i,"FPR"]])/naiveDataStats[[i,"FPR"]]),digits = 4)
    }else{
      FPRdiffProc = -Inf
    }
    if(naiveDataStats[[i,"TPR"]]!=0){
      TPRdiffProc = round(100*((classifierDataStats[[i,"TPR"]] - naiveDataStats[[i,"TPR"]])/naiveDataStats[[i,"TPR"]]),digits = 4)
    }else{
      TPRdiffProc = Inf
    }
    params <- list(
      FPRdiff = FPRdiffProc,
      TPRdiff = TPRdiffProc
    )
    return(params)
  })
  res <- rbindlist(res)
  res <- cbind(naiveDataStats[,c("lower","upper")],res)
  return(res)
}

#TODO: input: plot title etc.
evaluateClassifierPlot <- function(evaluationData=NULL, plotTitle="", plotFolderPath=""){
  out <- tryCatch(
    {
      eval <- copy(evaluationData)
      series <- colnames(eval)[colnames(eval)!=c("lower","upper")]
      eval[,subset:=as.factor(paste(lower,upper,sep = "-"))]
      eval$subset <- factor(eval$subset, levels=eval[,subset,by=.GRP]$subset)
      eval[, `:=`(lower=NULL,upper=NULL)]
      #eval <- melt(eval, id.vars="subset", measure.vars = c("FPRdiff","TPRdiff"), variable.name = "metric")
      eval <- melt(eval, id.vars="subset", measure.vars = series, variable.name = "metric")
      evalPlot <- ggplot(eval, aes(subset, value, fill=metric)) + 
        geom_bar(stat="identity", position="dodge") + 
        geom_text(aes(label=round(value,digits = 2)), vjust=-0.3, position=position_dodge(1.0), size=2.5) +
        scale_fill_manual(values = c("#0000FF","#FF0000")) +
        theme(axis.text.x = element_text(face="bold", color="#993333", angle=45),
              axis.text.y = element_text(face="bold", color="#993333", angle=45)) + 
        #scale_y_continuous(labels=scales::percent_format()) +
        ylab("value [%]") +
        xlab("DP subset") +
        ggtitle(plotTitle)
      
      if(!dir.exists(plotFolderPath)) dir.create(plotFolderPath)
      pdf(paste0(plotFolderPath,str_replace_all(plotTitle," ","_"),".pdf"))
      plot(evalPlot)
      dev.off()
      if(verbose==TRUE){print(paste0("Plot ",plotFolderPath,str_replace_all(plotTitle," ","_"),".pdf saved."))}  
    },
    error=function(cond) {
      message("Something went wrong. Function: evaluateClassifierPlot")
      message("This function plots the results of classifier evaluation.")
      message("Returns: the plot with results of classifier evaluation. Saves the plot to external file.")
      message("Original error message:\n")
      message(cond)
      return(NA)
    },
    finally={
      return(evalPlot)
    }
  )  
  return(out)
}



#main part of the script
#TODO: put HC/stats and everything else that still needs to be done here
bamFile <- loadBamFile(bamPath=bamPath)		#gunzipPath, libs: data.table, Rsamtools, testthat
indexFile <- loadIndexFile(bamPath=bamPath)		#libs: data.table, Rsamtools, testthat
coverage <- getCoverage(outputCoveragePath=outputCoveragePath, bamPath=bamPath, new=newCoverage)		#gunzipPath, makeBed, getCoverageFromBed libs: data.table, testthat
varDb <- loadVarDb(inputVariantsDb=inputVariantsDbPath)		#gunzipPath, getFileExtension, loadVcfData, loadTsvData, normalizeDeletions, libs: data.table, testthat
positionsAll <- getPositions(positionsAll=outputPositionsAllPath, varDb=varDb, new=newPostitionsAll)		#getUniqueLoci, libs: data.table, testthat

pileups <- doPileups(outputPositionsPath=outputPositionsAllPath,outputPileupsPath=outputPileupsPath,outputPileupsWithStatsPath=outputPileupsWithStatsPath,bamPath=bamPath,phredOffset=phredOffset,minBaseQuality=minBaseQuality,minMapq=minMapq,newPileups=newPileups,newPileupsWithStats=newPileupsWithStats)		#pileupVariantDatabase, addStatsColumnsToPileup libs: data.table, testthat


if(train==TRUE){
  subsetDpDataTable <- initializeSubsetDpDataTable(lowerBounds = lowerSubsetBounds, upperBounds = upperSubsetBounds)		#libs: data.table
  subsetDpDataTable <- updateDepthsColumn(subsetBounds=subsetDpDataTable, updateDp=TRUE, updateDpRand=FALSE, pileupsWithNucCols=pileups)

  positionsAllFalse <- getFalsePositions(positionsAllFalse=outputPositionsAllFalsePath, new=newPositionsAllFalse, coverage=coverage,
  										variantsSetToExclude=positionsAll, subsetDpDataTable=subsetDpDataTable)		#getRandomPositionsExcludeSubset, libs: data.table, testthat
  print("false pileups:")
  falsePileups <- doPileups(outputPositionsPath=outputPositionsAllFalsePath, outputPileupsPath=outputPileupsFalsePath,outputPileupsWithStatsPath=outputFalsePileupsWithStatsPath,bamPath=bamPath,phredOffset=phredOffset,minBaseQuality=minBaseQuality,minMapq=minMapq,new=newFalsePileups,newPileupsWithStats=newFalsePileupsWithStats)		#pileupVariantDatabase, addStatsColumnsToPileup libs: data.table, testthat
  subsetDpDataTable <- updateDepthsColumn(subsetBounds=subsetDpDataTable, updateDp=TRUE, updateDpRand=TRUE, pileupsWithNucCols=pileups, falsePileupsWithNucCols = falsePileups)
  #createHistogram(subsetDpDataTable=subsetDpDataTable, pileups=pileups, falsePileups=falsePileups)		#initializeSubsetDpDataTable, getSubsetDepths, libs: data.table

  trainingData <- prepareTrainingData(pileups=pileups, falsePileups=falsePileups)
  model <- trainModel(modelPath=modelPath, data=trainingData,listOfTrainArguments=listOfTrainArguments)
  
  #test the model
  baseline <- classifyNaive(newData = trainingData)
  baselineStats <- getClassificationStats(subsetBounds = subsetDpDataTable, classificationData = baseline)
  newModelClassificationData <- classify(model=model, newData = trainingData)
  newModelStats <- getClassificationStats(subsetBounds = subsetDpDataTable, classificationData = newModelClassificationData)
  
  eval <- evaluateClassifier(naiveDataStats = baselineStats, classifierDataStats = newModelStats)
  #eval <- evaluateClassifier(naiveDataStats = newModelStats, classifierDataStats = baselineStats)
  plot <- evaluateClassifierPlot(evaluationData = eval, plotTitle=plotTitle, plotFolderPath=plotFolderPath)
  plot
  #get the regular classification result
  classificationResult <- classify(model=model, newData=pileups)
}else{
  classificationResult <- classify(model=modelPath, newData=pileups)
}
varDb.merged <- mergePileupsWithVarDb(pileups=classificationResult, varDb=varDb, allx=TRUE,ally=TRUE)
if(verbose==TRUE){print("Pileup merged with varDb.")}
varDb.merged <- calculateDepths(varDb.merged)
if(verbose==TRUE){print("Depths calculated.")}

#TODO: sprawdzic kolumny
if(fullMergedDb==TRUE){ 
  fwrite(varDb.merged, outputVarDbMergedPath)
}else{
  fwrite(varDb.merged[,c("CHROM","POS","REF","ALT","countA","countC","countG","countT","countIn","countDel","countN","RD","AD","varCount","nullDepth","DP","errFlag")], outputVarDbMergedPath)
}
  
if(verbose==TRUE){print("varDb.merged exported to external file.")}




res <- list()
if(returnPileups==TRUE){res <- append(res,list(pileups=pileups, falsePileups=falsePileups))}
if(returnPositionsAll==TRUE){res <- append(res, list(positionsAll=positionsAll))}
if(returnPositionsAllFalse==TRUE){res <- append(res, list(positionsAllFalse=positionsAllFalse))} #TODO: tylko jak liczone
#if(returnStats==TRUE){res <- append(res,list(stats=stats))}
if(returnMergedDb==TRUE){res <- append(res,list(varDb.merged=varDb.merged))}
if(returnSubsetDpDataTable==TRUE){res <- append(res,list(subsetDpDataTable=subsetDpDataTable))} #TODO: tylko jak liczone
  
return(res)

}