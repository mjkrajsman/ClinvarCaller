library(testthat)

pileupsTests <- function(pileups){
  # deletions
  # rbind(subset(pileups,(CHROM=="X" & POS=="45076680")),subset(pileups,(CHROM=="1" & POS=="201372048")),subset(pileups,(CHROM=="17" & POS=="3513850")),subset(pileups,(CHROM=="X" & POS=="70028047")))
  test_that("Is chr1:201372048 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="201372048"))[,"nucA"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="201372048"))[,"nucC"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="201372048"))[,"nucG"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="201372048"))[,"nucT"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="201372048"))[,"nucIn"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="201372048"))[,"nucDel"]), equals(1)) #1
  })
  test_that("Is chr17:3513850 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="17" & POS=="3513850"))[,"nucA"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="17" & POS=="3513850"))[,"nucC"]), equals(60)) #61
    expect_that(as.integer(subset(pileups,(CHROM=="17" & POS=="3513850"))[,"nucG"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="17" & POS=="3513850"))[,"nucT"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="17" & POS=="3513850"))[,"nucIn"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="17" & POS=="3513850"))[,"nucDel"]), equals(48)) #48
  })
  test_that("Is chr22:41178966 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="22" & POS=="41178966"))[,"nucA"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="22" & POS=="41178966"))[,"nucC"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="22" & POS=="41178966"))[,"nucG"]), equals(112)) #114
    expect_that(as.integer(subset(pileups,(CHROM=="22" & POS=="41178966"))[,"nucT"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="22" & POS=="41178966"))[,"nucIn"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="22" & POS=="41178966"))[,"nucDel"]), equals(70)) #71
  })
  test_that("Is chrX:45076680 valid?", {
    #expect_that(as.integer(subset(pileups,(seqnames=="chrX" & pos=="45076680" & nucleotide=="-"))[,"count"]), equals(79))
    expect_that(as.integer(subset(pileups,(CHROM=="X" & POS=="45076680"))[,"nucA"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="X" & POS=="45076680"))[,"nucC"]), equals(5)) #5
    expect_that(as.integer(subset(pileups,(CHROM=="X" & POS=="45076680"))[,"nucG"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="X" & POS=="45076680"))[,"nucT"]), equals(38)) #38
    expect_that(as.integer(subset(pileups,(CHROM=="X" & POS=="45076680"))[,"nucIn"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="X" & POS=="45076680"))[,"nucDel"]), equals(79)) #80
  })
  test_that("Is chrX:70028047 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="X" & POS=="70028047"))[,"nucA"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="X" & POS=="70028047"))[,"nucC"]), equals(59)) #59
    expect_that(as.integer(subset(pileups,(CHROM=="X" & POS=="70028047"))[,"nucG"]), equals(6)) #6
    expect_that(as.integer(subset(pileups,(CHROM=="X" & POS=="70028047"))[,"nucT"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="X" & POS=="70028047"))[,"nucIn"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="X" & POS=="70028047"))[,"nucDel"]), equals(58)) #58
  })
  # insertions
  # rbind(subset(pileups,(CHROM=="1" & POS=="209618670")),subset(pileups,(CHROM=="2" & POS=="38997474")),subset(pileups,(CHROM=="8" & POS=="143931520")),subset(pileups,(CHROM=="19" & POS=="8585059")),subset(pileups,(CHROM=="19" & POS=="45363958")))  
  test_that("Is chr1:209618670 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="209618670"))[,"nucA"]), equals(124)) #124
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="209618670"))[,"nucC"]), equals(1)) #1
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="209618670"))[,"nucG"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="209618670"))[,"nucT"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="209618670"))[,"nucIn"]), equals(119)) #119
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="209618670"))[,"nucDel"]), equals(0)) #0
  })
  test_that("Is chr2:38997474 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="2" & POS=="38997474"))[,"nucA"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="2" & POS=="38997474"))[,"nucC"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="2" & POS=="38997474"))[,"nucG"]), equals(106)) #108
    expect_that(as.integer(subset(pileups,(CHROM=="2" & POS=="38997474"))[,"nucT"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="2" & POS=="38997474"))[,"nucIn"]), equals(84)) #84
    expect_that(as.integer(subset(pileups,(CHROM=="2" & POS=="38997474"))[,"nucDel"]), equals(0)) #0
  })
  test_that("Is chr8:143931520 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="8" & POS=="143931520"))[,"nucA"]), equals(83)) #83
    expect_that(as.integer(subset(pileups,(CHROM=="8" & POS=="143931520"))[,"nucC"]), equals(1)) #1
    expect_that(as.integer(subset(pileups,(CHROM=="8" & POS=="143931520"))[,"nucG"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="8" & POS=="143931520"))[,"nucT"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="8" & POS=="143931520"))[,"nucIn"]), equals(74)) #74
    expect_that(as.integer(subset(pileups,(CHROM=="8" & POS=="143931520"))[,"nucDel"]), equals(5)) #5
  })
  test_that("Is chr19:8585059 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="19" & POS=="8585059"))[,"nucA"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="19" & POS=="8585059"))[,"nucC"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="19" & POS=="8585059"))[,"nucG"]), equals(1)) #1
    expect_that(as.integer(subset(pileups,(CHROM=="19" & POS=="8585059"))[,"nucT"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="19" & POS=="8585059"))[,"nucIn"]), equals(1)) #1
    expect_that(as.integer(subset(pileups,(CHROM=="19" & POS=="8585059"))[,"nucDel"]), equals(0)) #0
  })
  test_that("Is chr19:45363958 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="19" & POS=="45363958"))[,"nucA"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="19" & POS=="45363958"))[,"nucC"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="19" & POS=="45363958"))[,"nucG"]), equals(2)) #3
    expect_that(as.integer(subset(pileups,(CHROM=="19" & POS=="45363958"))[,"nucT"]), equals(59)) #62
    expect_that(as.integer(subset(pileups,(CHROM=="19" & POS=="45363958"))[,"nucIn"]), equals(28)) #28
    expect_that(as.integer(subset(pileups,(CHROM=="19" & POS=="45363958"))[,"nucDel"]), equals(0)) #0
  })
  # SNP A
  # rbind(subset(pileups,(CHROM=="10" & POS=="121479598")),subset(pileups,(CHROM=="12" & POS=="8851916")),subset(pileups,(CHROM=="16" & POS=="30723642")))  
  test_that("Is chr10:121479598 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="10" & POS=="121479598"))[,"nucA"]), equals(203)) #204
    expect_that(as.integer(subset(pileups,(CHROM=="10" & POS=="121479598"))[,"nucC"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="10" & POS=="121479598"))[,"nucG"]), equals(1)) #1
    expect_that(as.integer(subset(pileups,(CHROM=="10" & POS=="121479598"))[,"nucT"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="10" & POS=="121479598"))[,"nucIn"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="10" & POS=="121479598"))[,"nucDel"]), equals(0)) #0
  })
  test_that("Is chr12:8851916 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="12" & POS=="8851916"))[,"nucA"]), equals(251)) #287
    expect_that(as.integer(subset(pileups,(CHROM=="12" & POS=="8851916"))[,"nucC"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="12" & POS=="8851916"))[,"nucG"]), equals(0)) #1
    expect_that(as.integer(subset(pileups,(CHROM=="12" & POS=="8851916"))[,"nucT"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="12" & POS=="8851916"))[,"nucIn"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="12" & POS=="8851916"))[,"nucDel"]), equals(0)) #0
  })
  test_that("Is chr16:30723642 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="16" & POS=="30723642"))[,"nucA"]), equals(205)) #206
    expect_that(as.integer(subset(pileups,(CHROM=="16" & POS=="30723642"))[,"nucC"]), equals(0)) #1
    expect_that(as.integer(subset(pileups,(CHROM=="16" & POS=="30723642"))[,"nucG"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="16" & POS=="30723642"))[,"nucT"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="16" & POS=="30723642"))[,"nucIn"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="16" & POS=="30723642"))[,"nucDel"]), equals(0)) #0
  })
  # SNP C
  # rbind(subset(pileups,(CHROM=="5" & POS=="157472861")),subset(pileups,(CHROM=="6" & POS=="7585734")),subset(pileups,(CHROM=="6" & POS=="36307681")))  
  test_that("Is chr5:157472861 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="5" & POS=="157472861"))[,"nucA"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="5" & POS=="157472861"))[,"nucC"]), equals(252)) #270
    expect_that(as.integer(subset(pileups,(CHROM=="5" & POS=="157472861"))[,"nucG"]), equals(0)) #1
    expect_that(as.integer(subset(pileups,(CHROM=="5" & POS=="157472861"))[,"nucT"]), equals(0)) #2
    expect_that(as.integer(subset(pileups,(CHROM=="5" & POS=="157472861"))[,"nucIn"]), equals(0)) #2
    expect_that(as.integer(subset(pileups,(CHROM=="5" & POS=="157472861"))[,"nucDel"]), equals(0)) #2
  })
  test_that("Is chr6:7585734 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="7585734"))[,"nucA"]), equals(0)) #1
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="7585734"))[,"nucC"]), equals(238)) #239
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="7585734"))[,"nucG"]), equals(1)) #1
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="7585734"))[,"nucT"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="7585734"))[,"nucIn"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="7585734"))[,"nucDel"]), equals(0)) #0
  })  
  test_that("Is chr6:36307681 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="36307681"))[,"nucA"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="36307681"))[,"nucC"]), equals(230)) #235
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="36307681"))[,"nucG"]), equals(0)) #1
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="36307681"))[,"nucT"]), equals(0)) #2
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="36307681"))[,"nucIn"]), equals(0)) #2
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="36307681"))[,"nucDel"]), equals(0)) #2
  }) 
  # SNP G
  # rbind(subset(pileups,(CHROM=="1" & POS=="155910782")),subset(pileups,(CHROM=="3" & POS=="190312891")),subset(pileups,(CHROM=="11" & POS=="88312415")))  
  test_that("Is chr1:155910782 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="155910782"))[,"nucA"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="155910782"))[,"nucC"]), equals(0)) #1
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="155910782"))[,"nucG"]), equals(246)) #258
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="155910782"))[,"nucT"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="155910782"))[,"nucIn"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="155910782"))[,"nucDel"]), equals(0)) #0
  }) 
  test_that("Is chr3:190312891 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="3" & POS=="190312891"))[,"nucA"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="3" & POS=="190312891"))[,"nucC"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="3" & POS=="190312891"))[,"nucG"]), equals(245)) #247
    expect_that(as.integer(subset(pileups,(CHROM=="3" & POS=="190312891"))[,"nucT"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="3" & POS=="190312891"))[,"nucIn"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="3" & POS=="190312891"))[,"nucDel"]), equals(0)) #0
  })
  test_that("Is chr11:88312415 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="11" & POS=="88312415"))[,"nucA"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="11" & POS=="88312415"))[,"nucC"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="11" & POS=="88312415"))[,"nucG"]), equals(204)) #206
    expect_that(as.integer(subset(pileups,(CHROM=="11" & POS=="88312415"))[,"nucT"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="11" & POS=="88312415"))[,"nucIn"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="11" & POS=="88312415"))[,"nucDel"]), equals(0)) #0
  })
  # SNP T
  # rbind(subset(pileups,(CHROM=="2" & POS=="219565481")),subset(pileups,(CHROM=="6" & POS=="7584384")),subset(pileups,(CHROM=="8" & POS=="143919156")))  
  test_that("Is chr2:219565481 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="2" & POS=="219565481"))[,"nucA"]), equals(0)) #1
    expect_that(as.integer(subset(pileups,(CHROM=="2" & POS=="219565481"))[,"nucC"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="2" & POS=="219565481"))[,"nucG"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="2" & POS=="219565481"))[,"nucT"]), equals(229)) #233
    expect_that(as.integer(subset(pileups,(CHROM=="2" & POS=="219565481"))[,"nucIn"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="2" & POS=="219565481"))[,"nucDel"]), equals(0)) #0
  })
  test_that("Is chr6:7584384 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="7584384"))[,"nucA"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="7584384"))[,"nucC"]), equals(1)) #3
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="7584384"))[,"nucG"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="7584384"))[,"nucT"]), equals(247)) #280
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="7584384"))[,"nucIn"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="7584384"))[,"nucDel"]), equals(0)) #0
  })
  test_that("Is chr8:143919156 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="8" & POS=="143919156"))[,"nucA"]), equals(0)) #1
    expect_that(as.integer(subset(pileups,(CHROM=="8" & POS=="143919156"))[,"nucC"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="8" & POS=="143919156"))[,"nucG"]), equals(0)) #1
    expect_that(as.integer(subset(pileups,(CHROM=="8" & POS=="143919156"))[,"nucT"]), equals(248)) #254
    expect_that(as.integer(subset(pileups,(CHROM=="8" & POS=="143919156"))[,"nucIn"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="8" & POS=="143919156"))[,"nucDel"]), equals(0)) #0
  })
  # Regular A
  # rbind(subset(pileups,(CHROM=="6" & POS=="43046304")),subset(pileups,(CHROM=="6" & POS=="7579851")),subset(pileups,(CHROM=="17" & POS=="31327637")))
  test_that("Is chr6:43046304 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046304"))[,"nucA"]), equals(268)) #422
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046304"))[,"nucC"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046304"))[,"nucG"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046304"))[,"nucT"]), equals(0)) #2
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046304"))[,"nucIn"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046304"))[,"nucDel"]), equals(0)) #0
  })
  test_that("Is chr6:7579851 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="7579851"))[,"nucA"]), equals(266)) #319
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="7579851"))[,"nucC"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="7579851"))[,"nucG"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="7579851"))[,"nucT"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="7579851"))[,"nucIn"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="7579851"))[,"nucDel"]), equals(0)) #0
  })
  test_that("Is chr17:31327637 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="17" & POS=="31327637"))[,"nucA"]), equals(267)) #317
    expect_that(as.integer(subset(pileups,(CHROM=="17" & POS=="31327637"))[,"nucC"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="17" & POS=="31327637"))[,"nucG"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="17" & POS=="31327637"))[,"nucT"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="17" & POS=="31327637"))[,"nucIn"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="17" & POS=="31327637"))[,"nucDel"]), equals(0)) #0
  })
  # Regular C
  # rbind(subset(pileups,(CHROM=="6" & POS=="43046251")),subset(pileups,(CHROM=="6" & POS=="43046310")),subset(pileups,(CHROM=="6" & POS=="43046341")))
  test_that("Is chr6:43046251 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046251"))[,"nucA"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046251"))[,"nucC"]), equals(270)) #342
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046251"))[,"nucG"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046251"))[,"nucT"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046251"))[,"nucIn"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046251"))[,"nucDel"]), equals(0)) #0
  })
  test_that("Is chr6:43046310 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046310"))[,"nucA"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046310"))[,"nucC"]), equals(274)) #421
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046310"))[,"nucG"]), equals(0)) #1
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046310"))[,"nucT"]), equals(0)) #1
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046310"))[,"nucIn"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046310"))[,"nucDel"]), equals(0)) #0
  })
    test_that("Is chr6:43046341 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046341"))[,"nucA"]), equals(1)) #1
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046341"))[,"nucC"]), equals(278)) #406
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046341"))[,"nucG"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046341"))[,"nucT"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046341"))[,"nucIn"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046341"))[,"nucDel"]), equals(0)) #0
  })
  # Regular G
  # rbind(subset(pileups,(CHROM=="1" & POS=="152308647")),subset(pileups,(CHROM=="1" & POS=="152308777")),subset(pileups,(CHROM=="1" & POS=="152308956")))
  test_that("Is chr1:152308647 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="152308647"))[,"nucA"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="152308647"))[,"nucC"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="152308647"))[,"nucG"]), equals(281)) #476
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="152308647"))[,"nucT"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="152308647"))[,"nucIn"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="152308647"))[,"nucDel"]), equals(0)) #0
  })
  test_that("Is chr1:152308777 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="152308777"))[,"nucA"]), equals(1)) #1
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="152308777"))[,"nucC"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="152308777"))[,"nucG"]), equals(274)) #452
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="152308777"))[,"nucT"]), equals(0)) #1
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="152308777"))[,"nucIn"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="152308777"))[,"nucDel"]), equals(0)) #0
  })
  test_that("Is chr1:152308956 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="152308956"))[,"nucA"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="152308956"))[,"nucC"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="152308956"))[,"nucG"]), equals(275)) #421
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="152308956"))[,"nucT"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="152308956"))[,"nucIn"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="152308956"))[,"nucDel"]), equals(0)) #0
  })
  # Regular T
  # rbind(subset(pileups,(CHROM=="1" & POS=="152308940")),subset(pileups,(CHROM=="6" & POS=="43046293")),subset(pileups,(CHROM=="17" & POS=="31327615")))
  test_that("Is chr1:152308940 valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="152308940"))[,"nucA"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="152308940"))[,"nucC"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="152308940"))[,"nucG"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="152308940"))[,"nucT"]), equals(275)) #402
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="152308940"))[,"nucIn"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="1" & POS=="152308940"))[,"nucDel"]), equals(0)) #0
  })
  test_that("Is chr6: valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046293"))[,"nucA"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046293"))[,"nucC"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046293"))[,"nucG"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046293"))[,"nucT"]), equals(274)) #414
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046293"))[,"nucIn"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="6" & POS=="43046293"))[,"nucDel"]), equals(0)) #0
  })
  test_that("Is chr17: valid?", {
    expect_that(as.integer(subset(pileups,(CHROM=="17" & POS=="31327615"))[,"nucA"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="17" & POS=="31327615"))[,"nucC"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="17" & POS=="31327615"))[,"nucG"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="17" & POS=="31327615"))[,"nucT"]), equals(268)) #332
    expect_that(as.integer(subset(pileups,(CHROM=="17" & POS=="31327615"))[,"nucIn"]), equals(0)) #0
    expect_that(as.integer(subset(pileups,(CHROM=="17" & POS=="31327615"))[,"nucDel"]), equals(0)) #0
  })
  print("Pileups tests passed.")
}