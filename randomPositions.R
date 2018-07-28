library(data.table)
library(Rsamtools)
library(tictoc)

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

getCoverageFromBed <- function(bedPath="./input/corriell_S7-ready.per-base.bed.gz"){
  coverage <- fread(bedPath)
  colnames(coverage) <- c("CHROM","start","stop","count")
  coverage <- subset(coverage, count > 0)
  coverage$CHROM <- gsub("chr", "", coverage$CHROM)
  coverage <- filterChroms(coverage)
  coverage[,"start":=start+1L]
  #coverage[,"length":=stop-start]
  return(coverage)
}

getCoverageSubsetByCount <- function(coverage, min=1, max=1){
  if(min<=max){
    return(subset(coverage, (count >= min & count <= max)))
  }else{
    print("max must be greater or equal to min!")
    return(NULL)
  }  
}



getPositionsFromCoverageSubset <- function(coverageSubset,subsetMaxLength=200000){
  if(nrow(coverageSubset)>subsetMaxLength){coverageSubset <- coverageSubset[sample(.N, subsetMaxLength, replace=F)]}
  #tic()
  #res2 <- rbindlist(lapply(1:nrow(f), function(i){data.table("CHROM"=f[i,CHROM],"POS"=f[i,start]:f[i,stop],"count"=f[i,count])}))
  #toc()
  
  #tic()
  res <- rbindlist(
      lapply(c(1:22, "X", "Y", "M"), function(i){
      chromSubset <- subset(coverageSubset, CHROM==i)
      print(paste0("chr",i))
      if(nrow(chromSubset)){
        chromPositions <- rbindlist(
          lapply(min(chromSubset[,count]):max(chromSubset[,count]), function(j){
          countSubset <- subset(chromSubset, count==j)
          print(paste0("count",j))
          if(nrow(countSubset)){
            countPositions <- rbindlist(
              lapply(1:nrow(countSubset), function(k){
              data.table("POS"=countSubset[k,start]:countSubset[k,stop])
              }))
            countPositions[,"count":=j]
            return(countPositions)
          }  
        }))
        chromPositions[,"CHROM":=i] 
        return(chromPositions)
      }
    }))
  #toc()  
  #setequal(res[order("CHROM","POS","count")],res2[order("CHROM","POS","count")])
  return(res)
}

getRandomPositions <- function(listOfPositions, randomSubsetSize=10000){
  res <- listOfPositions[sample(.N, randomSubsetSize, replace=F)]
  return(res)
}

getRandomPositionsFromCoverage <- function(coverage,subsetBounds,subsetMaxLength=200000,randomSubsetSize=10000){
  res <- lapply(1:nrow(subsetBounds), function(i){
    coverageSubset <- getCoverageSubsetByCount(coverage,min=subsetBounds[i,lower],max=subsetBounds[i,upper])
    listOfPositions <- getPositionsFromCoverageSubset(coverageSubset,subsetMaxLength)
    randomPositions <- getRandomPositions(listOfPositions, randomSubsetSize)
    return(randomPositions)
  })
  res <- rbindlist(res)
  return(res)
}

getPositionsFromCoverage <- function(coverage,subsetBounds,subsetMaxLength=200000){
  res <- lapply(1:nrow(subsetBounds), function(i){
    coverageSubset <- getCoverageSubsetByCount(coverage,min=subsetBounds[i,lower],max=subsetBounds[i,upper])
    listOfPositions <- getPositionsFromCoverageSubset(coverageSubset,subsetMaxLength)
    return(listOfPositions)
  })
  return(res)
}





bedPath <- gunzipPath("./input/corriell_S7-ready.per-base.bed.gz")
coverage <- getCoverageFromBed(bedPath)
subsetBounds <- data.table("lower"=c(1,2,6,11,21,51,101,201),"upper"=c(1,5,10,20,50,100,200,9999))

#subsetBounds <- data.table("lower"=c(21,26),"upper"=c(25,30))

randomPositions <- getRandomPositionsFromCoverage(coverage,subsetBounds,subsetMaxLength=300000,randomSubsetSize=200000)







#rm(cov)
#cov <- list(
#  "cov1"=subset(coverage, count == 1),
#  "cov2to5"=subset(coverage, (count >= 2 & count <= 5)),
#  "cov6to10"=subset(coverage, (count >= 6 & count <= 10)),
#  "cov11to20"=subset(coverage, (count >= 11 & count <= 20)),
#  "cov21to50"=subset(coverage, (count >= 21 & count <= 50)),
#  "cov51to100"=subset(coverage, (count >= 51 & count <= 100)),
#  "cov101to200"=subset(coverage, (count >= 101 & count <= 200)),
#  "cov200plus"=subset(coverage, count > 200)
#)
#cov$lengths <- lapply(cov, function(i){
#    i[,(sum(stop-start+1))]
#})

#fwrite(randomPositions, "./output/randomPositions1.csv")
#tic()
#allPositionsFromCoverage <- getPositionsFromCoverage(coverage, subsetBounds, subsetMaxLength = 1000000)
#toc()
#fwrite(allPositionsFromCoverage, "./output/allPositionsFromCoverage.csv")


#================================================================== coverage ^ ==================================================================
