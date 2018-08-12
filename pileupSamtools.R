library(data.table)
library(Rsamtools)
library(tictoc)
library(stringr)
library(testthat)
library(BSgenome.Hsapiens.UCSC.hg38)

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
addStatsColumnsToPileup <- function(pileups, countsOnly=FALSE){
  res <- copy(pileups)
  tic(paste0("Stats calculated. Pileups, nrow: ", nrow(res), " time"))
  if(countsOnly==TRUE){
    res[,c("CHROM","POS","nucCount_A","nucCount_C","nucCount_G","nucCount_T","nucCount_In","nucCount_Del"):=calculateNucStats(chrom=CHROM, pos=POS, nuc=NUC, countsOnly=TRUE, bq=BQ, mq=MQ, pr=posRead), by=1:nrow(res)]
  }else{
    res[,c("CHROM","POS","nucCount_A","meanBQ_A","sdBQ_A","meanMQ_A","sdMQ_A","meanPosRead_A","sdPosRead_A","nucCount_C","meanBQ_C","sdBQ_C","meanMQ_C","sdMQ_C","meanPosRead_C","sdPosRead_C","nucCount_G","meanBQ_G","sdBQ_G","meanMQ_G","sdMQ_G","meanPosRead_G","sdPosRead_G","nucCount_T","meanBQ_T","sdBQ_T","meanMQ_T","sdMQ_T","meanPosRead_T","sdPosRead_T","nucCount_In","meanBQ_In","sdBQ_In","meanMQ_In","sdMQ_In","meanPosRead_In","sdPosRead_In","nucCount_Del","meanBQ_Del","sdBQ_Del","meanMQ_Del","sdMQ_Del","meanPosRead_Del","sdPosRead_Del"):=
          calculateNucStats(chrom=CHROM, pos=POS, nuc=NUC, countsOnly=FALSE, bq=BQ, mq=MQ, pr=posRead), by=1:nrow(res)]
  }  
  toc()
  return(res)
}

cmd <- "samtools mpileup -A -d 0 -l ./output/lociAll.positions -q 0 -Q 0 -O -s ./input/corriell_S7-ready.bam > ./output/samtoolsPileups.mpileup"
system(cmd)

pileups <- fread("./output/samtoolsPileups.mpileup", quote="")
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
#[ACGTNacgtn]{1}
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

View(pileups)

pileups <- addStatsColumnsToPileup(pileups, countsOnly = TRUE)

tic()
pileups[,c("CHROM","POS","nucCount_A","meanBQ_A","sdBQ_A","meanMQ_A","sdMQ_A","meanPosRead_A","sdPosRead_A","nucCount_C","meanBQ_C","sdBQ_C","meanMQ_C","sdMQ_C","meanPosRead_C","sdPosRead_C","nucCount_G","meanBQ_G","sdBQ_G","meanMQ_G","sdMQ_G","meanPosRead_G","sdPosRead_G","nucCount_T","meanBQ_T","sdBQ_T","meanMQ_T","sdMQ_T","meanPosRead_T","sdPosRead_T","nucCount_In","meanBQ_In","sdBQ_In","meanMQ_In","sdMQ_In","meanPosRead_In","sdPosRead_In","nucCount_Del","meanBQ_Del","sdBQ_Del","meanMQ_Del","sdMQ_Del","meanPosRead_Del","sdPosRead_Del"):=calculateNucStats(chrom=CHROM, pos=POS, nuc=NUC, bq=BQ, mq=MQ, pr=posRead), by=1:nrow(pileups)]
toc()


tic()
pilSample[,c("CHROM","POS","nucCount_A","meanBQ_A","sdBQ_A","meanMQ_A","sdMQ_A","meanPosRead_A","sdPosRead_A","nucCount_C","meanBQ_C","sdBQ_C","meanMQ_C","sdMQ_C","meanPosRead_C","sdPosRead_C","nucCount_G","meanBQ_G","sdBQ_G","meanMQ_G","sdMQ_G","meanPosRead_G","sdPosRead_G","nucCount_T","meanBQ_T","sdBQ_T","meanMQ_T","sdMQ_T","meanPosRead_T","sdPosRead_T","nucCount_In","meanBQ_In","sdBQ_In","meanMQ_In","sdMQ_In","meanPosRead_In","sdPosRead_In","nucCount_Del","meanBQ_Del","sdBQ_Del","meanMQ_Del","sdMQ_Del","meanPosRead_Del","sdPosRead_Del"):=calculateNucStats(chrom=CHROM, pos=POS, nuc=NUC, bq=BQ, mq=MQ, pr=posRead), by=1:nrow(pilSample)]
toc()