library(data.table)
library(Rsamtools)
library(tools)

#wczytaj wynik CC - varDB.merged, zrob podejscie naiwne

#wez DP z CC, podziel na kubełki
#wez DP z HC, podziel na kubelki


HCCompare <- function(vcfHcGeneralPath="./input/corriell_S7-HCScan.vcf.gz",
                      vcfHcFcPath="./input/corriell_S7-HCScanFC.vcf.gz",
                      varDbMergedPath="./output/varDbMerged-reduced.csv",
                      outputStatsPath="./output/stats.csv",
                      verbose=TRUE){
  
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
  
  setTenthColumnNameInVcfToGenotype <- function(vcfTable=""){
    out <- tryCatch(
      {
        colnames(vcfTable)[10] <- "GENOTYPE"
      },
      error=function(cond) {
        message("Something went wrong. Function: setTenthColumnNameInVcfToGenotype")
        message("This function changes the input data.tables 10th column name to 'GENOTYPE'.")
        message("Returns: data.table with 10th column name changed to 'GENOTYPE'")
        message("Original error message:\n")
        message(cond)
        return(NA)
      },
      finally={
        return(vcfTable)
      }
    )  
    if(verbose==TRUE){print("Variant classification complete.")}
    return(out)
  }
  
  #TODO: try-catch
  #gunzipPath, loadVcfData, filterChroms, normalizeDeletions, setTenthColumnNameInVcfToGenotype
  getVcfHcGeneral <- function(vcfHcGeneralPath=""){
    # GATK HAPLOTYPECALLER - general scan
    if(file.exists(vcfHcGeneralPath)!=TRUE){
      stop("vcfHcGeneral file does not exist!")
    }
    vcfHcGeneral <- loadVcfData(vcfHcGeneralPath)
    if(verbose==TRUE){print("vcfHcGeneral file loaded.")}
    vcfHcGeneral <- setTenthColumnNameInVcfToGenotype(vcfHcGeneral)
    vcfHcGeneral <- filterChroms(vcfHcGeneral)
    vcfHcGeneral <- normalizeDeletions(vcfHcGeneral)
    if(unique(unique(vcfHcGeneral$CHROM) %in% c(1:22, "X", "Y", "M"))!=TRUE){
      stop("Invalid chromosomes (not: 1:22, X, Y, M) present in vcfHcGeneral!")
    }
    return(vcfHcGeneral)
  }
  
  #TODO: try-catch
  #gunzipPath, loadVcfData, filterChroms, normalizeDeletions, setTenthColumnNameInVcfToGenotype
  getVcfHcFc <- function(vcfHcFcPath=""){
    # GATK HAPLOTYPECALLER - forcecalling
    if(file.exists(vcfHcFcPath)!=TRUE){
      stop("vcfHcFc file does not exist!")
    }
    vcfHcFc <- loadVcfData(vcfHcFcPath)
    if(verbose==TRUE){print("vcfHcFc file loaded.")}
    vcfHcFc <- setTenthColumnNameInVcfToGenotype(vcfHcFc)
    vcfHcFc01 <- vcfHcFc[like(GENOTYPE,"0/1")]
    vcfHcFc11 <- vcfHcFc[like(GENOTYPE,"1/1")]
    vcfHcFc <- rbind(vcfHcFc01, vcfHcFc11)
    rm(vcfHcFc01, vcfHcFc11)
    
    vcfHcFc <- filterChroms(vcfHcFc)
    vcfHcFc <- normalizeDeletions(vcfHcFc)
    
    if(!setequal(unique(substr(vcfHcFc$GENOTYPE,1,3)),c("0/1","1/1"))){
      stop("GENOTYPE should be only 0/1 or 1/1!")
    }
    
    return(vcfHcFc)
  }
  
  
  #TODO: try-catch
  getDpFromVcf <- function(vcf=NULL){
    DP <- data.table(DP=as.integer(do.call(rbind,str_split(vcf[,"GENOTYPE"][[1]],":"))[,3]))
    res <- cbind(vcf,DP)
    return(res)
  }
  
  #TODO: try-catch
  #getDpFromVcf
  getHcTpInSubsets <- function(subsetBounds=NULL, haplotypeCallerVcf=NULL){
    cd <- copy(haplotypeCallerVcf)
    cd <- getDpFromVcf(vcf=cd)
    bounds <- copy(subsetBounds[,1:4])
    rowAll <- bounds[1]
    rowAll[,c("lower","upper"):=data.table(min(bounds$lower),max(bounds$upper))]
    bounds <- rbind(rowAll, bounds)
    res <- lapply(1:nrow(bounds), function(i){
      lower <- bounds[i,lower]
      upper <- bounds[i,upper]
      params <- list(
        lower = lower,
        upper = upper,
        TP = length(which(cd$DP >= lower & cd$DP <= upper))
      )
      return(params)
    })
    res <- rbindlist(res)
    return(res) 
  }
  
  #vcfHcGeneralPath="./input/corriell_S7-HCScan.vcf.gz"
  #vcfHcFcPath="./input/corriell_S7-HCScanFC.vcf.gz"
  #varDbMergedPath="./output/varDbMerged-reduced.csv"
  #outputStatsPath="./output/stats.csv"
  


    varDb.merged <- fread(varDbMergedPath)
  #TODO: HC fc/general should be loaded and processed here
    vcfHcGeneral <- getVcfHcGeneral(vcfHcGeneralPath=vcfHcGeneralPath)
    vcfHcFc <- getVcfHcFc(vcfHcFcPath=vcfHcFcPath)
    vcfHcGeneral <- getDpFromVcf(vcf=vcfHcGeneral)
    vcfHcFc <- getDpFromVcf(vcf=vcfHcFc)
    xtp <- getHcTpInSubsets(subsetBounds = histogramData,haplotypeCallerVcf = vcfHcFc)

  #TODO: tidy up everything below this line, put into functions, place before the main part of the script etc.
  #TODO: correct stats calculation
    # ===== TPR ===== # sensitivity: TPR = TP/(TP+FN)
    FNPlusTP <- nrow(varDb.merged) # P = TP + FN, total number of positives
    
    ccTP <- sum(varDb.merged[,varCount>0], na.rm=TRUE) # pileups on P set only; 
    ccFN <- sum(varDb.merged[,varCount==0], na.rm=TRUE)
    ccTPR <- ccTP/FNPlusTP # ClinvarCaller, forcecalling
    
    hcGeneralTP <- nrow(vcfHcGeneral[CHROM %in% varDb.merged$CHROM & POS %in% varDb.merged$POS])
    hcGeneralTPR <- hcGeneralTP/FNPlusTP # GATK HaplotypeCaller, general scan
      
    hcFcTP <- nrow(vcfHcFc[CHROM %in% varDb.merged$CHROM & POS %in% varDb.merged$POS]) # pileups on P set only, could be nrow(vcfHcFc)
    hcFcTPR <- hcFcTP/FNPlusTP # GATK HaplotypeCaller, forcecalling

    
    #====================================================================== TPR ^ ======================================================================
    #====================================================================== CC + HC sum v ======================================================================

      variantsCcPos<-subset(varDb.merged, varCount>0)
      variantsCcPos<-variantsCcPos[,c("CHROM","POS","REF","ALT")]
      
      variantsHcPos<-vcfHcFc[CHROM %in% varDb.merged$CHROM & POS %in% varDb.merged$POS,c("CHROM","POS","REF","ALT")]
      ccHcOR<-rbind(variantsCcPos,variantsHcPos)
      setkeyv(ccHcOR, c("CHROM","POS"))
      ccHcOR <- subset(unique(ccHcOR, by=c("CHROM","POS")))
      ccHcAND <- variantsHcPos[CHROM %in% variantsCcPos$CHROM & POS %in% variantsCcPos$POS] # variants in both cc and hc
      variantsHcOnly <- ccHcOR[!variantsCcPos]
      variantsCcOnly <- ccHcOR[!variantsHcPos]
      rm(variantsCcPos,variantsHcPos)
      
      variantsHcOnlyCount <- nrow(variantsHcOnly)
      variantsCcOnlyCount <- nrow(variantsCcOnly)
      ccHcORCount <- nrow(ccHcOR)
      ccHcANDCount <- nrow(ccHcAND)
      variantsHcOnlyInDelCount <- nrow(variantsHcOnly[nchar(as.character(ALT))!=nchar(as.character(REF)) & substr(as.character(REF),1,1)==substr(as.character(ALT),1,1)]) + nrow(variantsHcOnly[variantsHcOnly[,ALT=="-"]]) + nrow(variantsHcOnly[variantsHcOnly[,ALT=="+"]])
      variantsCcOnlyInDelCount <- nrow(variantsCcOnly[nchar(as.character(ALT))!=nchar(as.character(REF)) & substr(as.character(REF),1,1)==substr(as.character(ALT),1,1)]) + nrow(variantsCcOnly[variantsCcOnly[,ALT=="-"]]) + nrow(variantsCcOnly[variantsCcOnly[,ALT=="+"]])
      variantsHcOnlySnpCount <- nrow(variantsHcOnly[nchar(as.character(ALT))==nchar(as.character(REF)) & variantsHcOnly[,ALT!="-"] & variantsHcOnly[,ALT!="+"]])
      #unique(variantsHcOnly[nchar(as.character(ALT))==nchar(as.character(REF))][,REF]) # SNP check
      variantsCcOnlySnpCount <- nrow(variantsCcOnly[nchar(as.character(ALT))==nchar(as.character(REF)) & variantsCcOnly[,ALT!="-"] & variantsCcOnly[,ALT!="+"]])
      #unique(variantsCcOnly[nchar(as.character(ALT))==nchar(as.character(REF))][,REF]) # SNP check

    #====================================================================== CC + HC sum ^ ======================================================================
    

    stats <- data.table(FNPlusTP=FNPlusTP,ccTP=ccTP,ccFN=ccFN,ccTPR=ccTPR,hcGeneralTP=hcGeneralTP,hcFcTP=hcFcTP,
                        hcGeneralTPR=hcGeneralTPR,hcFcTPR=hcFcTPR,ccHcANDCount=ccHcANDCount, ccHcORCount=ccHcORCount,
                        variantsHcOnlyCount=variantsHcOnlyCount,variantsCcOnlyCount=variantsCcOnlyCount,
                        variantsHcOnlyInDelCount=variantsHcOnlyInDelCount,variantsCcOnlyInDelCount=variantsCcOnlyInDelCount,
                        variantsHcOnlySnpCount=variantsHcOnlySnpCount,variantsCcOnlySnpCount=variantsCcOnlySnpCount)

    fwrite(stats, outputStatsPath)
    if(verbose==TRUE){print("Stats ready.")}

    return(stats)

}