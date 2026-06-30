test_that("micro2intData (with min-max default) creates an microintData object with the correct dimensions and attributes", {
  MicroDt <- data.frame(X=1:9,Y=9:1)
  agrfct <- factor(c("A","B","A","C","B","C","B","A","A"))
  Idt <- micro2intData(MicroDt,agrfct)

  expect_s4_class(Idt,"intData")
  expect_equal(nrow(Idt),3)
  expect_equal(ncol(Idt),2)
  expect_equal(NbMicroUnits(Idt),c(A=4,B=3,C=2))
  expect_equal(names(Idt),c("X","Y"))
} )

test_that("micro2intData performs a correct agregation into (min-max based) intData objects", {
  MicroDt <- data.frame(X=1:9,Y=9:1)
  agrfct <- factor(c("A","B","A","C","B","C","B","A","A"))
  Idt <- micro2intData(MicroDt,agrfct)
  TrueIDt <- intData(data.frame(list(c(1,2,4),c(9,7,6),c(1,3,4),c(9,8,6))),
                   VarNames=c("X","Y"),ObsNames = c("A","B","C"),
                   NbMicroUnits=c(A=4,B=3,C=2), Seq = "LbUb_VarbyVar")
  
  expect_identical(Idt,TrueIDt)
} )

test_that("quantile based micro2intData creates an intData object with the correct dimensions and attributes", {
  MicroDt <- data.frame(X=rep(0:10,3),Y=rep(seq(0,100,by=10),3))
  agrfct <- factor(c(rep("A",6),rep("B",11),rep("A",5),rep("C",11)))
  Idt <- micro2intData(MicroDt,agrby=agrfct,agrcrt=c(0.1,0.9))

  expect_s4_class(Idt,"intData")
  expect_equal(nrow(Idt),3)
  expect_equal(ncol(Idt),2)
  expect_equal(NbMicroUnits(Idt),c(A=11,B=11,C=11))
  expect_equal(names(Idt),c("X","Y"))
} )

test_that("micro2intData performs a correct agregation into quantile based intData objects", {
  MicroDt <- data.frame(X=rep(0:10,3),Y=rep(seq(0,100,by=10),3))
  agrfct <- factor(c(rep("A",6),rep("B",11),rep("A",5),rep("C",11)))
  Idt <- micro2intData(MicroDt,agrby=agrfct,agrcrt=c(0.1,0.9))
  TrueIDt <- intData(data.frame(list(rep(1,3),rep(9,3),rep(10,3),rep(90,3))),
                   VarNames=c("X","Y"),ObsNames = c("A","B","C"),
                   NbMicroUnits=c(A=11,B=11,C=11), Seq = "LbUb_VarbyVar")
  
  expect_identical(Idt,TrueIDt)
} )

test_that("micro2intData stops if microdata is not a data frame", {
  MicroDt <- matrix(1:6, nrow=3, ncol=2)
  agrfct <- factor(c("A","B","C"))
  expect_error(micro2intData(MicroDt, agrfct), 
               "First argument of AgMicroData must be a data frame")
})

test_that("micro2intData stops if non-numeric columns exist in microdata", {
  MicroDt <- data.frame(X=c(1,2,3), Y=c("a","b","c"))
  agrfct <- factor(c("A","B","C"))
  expect_error(micro2intData(MicroDt, agrfct),
               "Some of the columns of the")
})

test_that("micro2intData stops if agrby is not a factor", {
  MicroDt <- data.frame(X=1:3, Y=3:1)
  agrby <- c("A", "B", "C")  # character vector, not a factor
  expect_error(micro2intData(MicroDt, agrby),
               "Argument agrby is not a factor")
})

test_that("micro2intData stops if length of agrby does not match number of rows", {
  MicroDt <- data.frame(X=1:5, Y=5:1)
  agrby <- factor(c("A", "B", "C"))  # only 3 levels but 5 rows
  expect_error(micro2intData(MicroDt, agrby),
               "Size of the agrby argument does not agree")
})

test_that("micro2intData stops if agrcrt has invalid length", {
  MicroDt <- data.frame(X=1:3, Y=3:1)
  agrby <- factor(c("A", "B", "C"))
  expect_error(micro2intData(MicroDt, agrby, agrcrt=c(0.1, 0.5, 0.9)),
               "Wrong value for the agrcrt argument")
})

test_that("micro2intData stops if agrcrt has invalid quantile bounds", {
  MicroDt <- data.frame(X=1:3, Y=3:1)
  agrby <- factor(c("A", "B", "C"))
  expect_error(micro2intData(MicroDt, agrby, agrcrt=c(0.5, 0.3)),
               "Wrong value for the agrcrt argument")
})

test_that("micro2intData stops if agrcrt quantile is out of range", {
  MicroDt <- data.frame(X=1:3, Y=3:1)
  agrby <- factor(c("A", "B", "C"))
  expect_error(micro2intData(MicroDt, agrby, agrcrt=c(-0.1, 0.9)),
               "Wrong value for the agrcrt argument")
})

test_that("micro2intData removes rows with all non-finite values with warning", {
  MicroDt <- data.frame(X=c(1, NA, 3), Y=c(4, NA, 2))
  agrby <- factor(c("A", "B", "A"))
  
  expect_warning(micro2intData(MicroDt, agrby),
                 "rows of the")
  result <- suppressWarnings(micro2intData(MicroDt, agrby))
  expect_equal(nrow(result), 1)
})

test_that("micro2intData removes groups with all NA in a variable with warning", {
  MicroDt <- data.frame(X=c(1, NA, 3), Y=c(4, 2, 1))
  agrby <- factor(c("A", "B", "A"))
  
  expect_warning(micro2intData(MicroDt, agrby),
                 "Removed .* groups with at least one variable fully NA")
  result <- suppressWarnings(micro2intData(MicroDt, agrby))
  expect_equal(nrow(result), 1)
  expect_equal(result@ObsNames, "A")
})

test_that("micro2intData returns NULL when all units lead to degenerate intervals", {
  # Create data where all intervals are degenerate (min == max)
  MicroDt <- data.frame(X=c(5,5,5,2,2,2), Y=c(10,10,10,20,20,20))
  agrby <- factor(c("A","A","A","B","B","B"))
  
  result <- suppressWarnings(micro2intData(MicroDt, agrby))
  expect_null(result)
})

test_that("micro2intData removes single degenerate interval with appropriate warning", {
  MicroDt <- data.frame(X=c(1, 2, 5, 5, 5), Y=c(1, 2, 10, 10, 10))
  agrby <- factor(c("A", "A", "B", "B", "B"))
  
  expect_warning(micro2intData(MicroDt, agrby),
                 "Data unit.*was eliminated because it lead to some degenerate intervals")
  result <- suppressWarnings(micro2intData(MicroDt, agrby))
  expect_equal(nrow(result), 1)
})

test_that("micro2intData removes multiple degenerate intervals (<10) with appropriate warning", {
  MicroDt <- data.frame(X=c(1, 2, 3, 3, 3, 4, 4, 4), 
                        Y=c(1, 2, 10, 10, 10, 20, 20, 20))
  agrby <- factor(c("A", "A", "B", "B", "B", "C", "C", "C"))
  
  expect_warning(micro2intData(MicroDt, agrby),
                 "Data units.*were eliminated because they lead to some degenerate intervals")
  result <- suppressWarnings(micro2intData(MicroDt, agrby))
  expect_equal(nrow(result), 1)
})

test_that("micro2intData removes many degenerate intervals (>=10) with count warning", {
  # Create 11 degenerate groups and 1 valid
  X <- c(rep(1, times=12), 2:20)
  Y <- c(rep(10, times=12), 20:2)
  MicroDt <- data.frame(X=X, Y=Y)
  agrby <- factor(c(rep(LETTERS[1:12], each=1), rep("M", 19)))
  
  expect_warning(micro2intData(MicroDt, agrby),
                 "12.*were eliminated because they lead to some degenerate intervals")
  result <- suppressWarnings(micro2intData(MicroDt, agrby))
  expect_equal(nrow(result), 1)
})

test_that("micro2intData detects U_id_symmetric case when LatentParam is NULL", {
  MicroDt <- data.frame(X=1:6, Y=6:1)
  agrby <- factor(c("A", "A", "B", "B", "C", "C"))
  
  result <- micro2intData(MicroDt, agrby, LatentParam=NULL)
  expect_s4_class(result, "intData")
  expect_equal(result@LatentCase, "U_id_symmetric")
})

test_that("micro2intData detects U_id_symmetric case when LatentParam has is a list of two numbers", {
  MicroDt <- data.frame(X=1:6, Y=6:1)
  agrby <- factor(c("A", "A", "B", "B", "C", "C"))
  
  result <- micro2intData(MicroDt, agrby, LatentParam=list(0, 1), LatentDist = "Triang")
  expect_s4_class(result, "intData")
  expect_equal(result@LatentCase, "U_id")
})

test_that("micro2intData detects U_id case when LatentParam is numeric", {
  MicroDt <- data.frame(X=1:6, Y=6:1)
  agrby <- factor(c("A", "A", "B", "B", "C", "C"))
  
  result <- micro2intData(MicroDt, agrby, LatentParam=list(0.5), LatentDist = "Unif")
  expect_s4_class(result, "intData")
  expect_equal(result@LatentCase, "U_id_symmetric")
})

test_that("micro2intData explicitly accepts LatentCase parameter", {
  MicroDt <- data.frame(X=1:6, Y=6:1)
  agrby <- factor(c("A", "A", "B", "B", "C", "C"))
  
  result <- micro2intData(MicroDt, agrby, LatentCase="U_id_symmetric")
  expect_s4_class(result, "intData")
  expect_equal(result@LatentCase, "U_id_symmetric")
})

test_that("micro2intData accepts General LatentCase", {
  MicroDt <- data.frame(X=1:9, Y=9:1)
  agrby <- factor(c("A", "A", "A", "B", "B", "B", "C", "C", "C"))
  
  result <- suppressWarnings(micro2intData(MicroDt, agrby, LatentCase="General"))
  expect_s4_class(result, "intData")
  expect_equal(result@LatentCase, "General")
})

test_that("micro2intData defaults LatentDist to KDE when LatentCase is General", {
  MicroDt <- data.frame(X=1:9, Y=9:1)
  agrby <- factor(c("A", "A", "A", "B", "B", "B", "C", "C", "C"))
  
  result <- suppressWarnings(micro2intData(MicroDt, agrby, LatentCase="General"))
  expect_equal(result@LatentDist, "KDE")
})

test_that("micro2intData accepts TriangParam", {
  MicroDt <- data.frame(X=1:6, Y=6:1)
  agrby <- factor(c("A", "A", "B", "B", "C", "C"))
  
  result <- micro2intData(MicroDt, agrby, 
                         LatentDist="Triang", TriangParam=0.3)
  expect_s4_class(result, "intData")
})

test_that("micro2intData accepts BetaParam.a and BetaParam.b", {
  MicroDt <- data.frame(X=1:6, Y=6:1)
  agrby <- factor(c("A", "A", "B", "B", "C", "C"))
  
  result <- micro2intData(MicroDt, agrby, LatentCase = "U_id",
                         LatentDist="Beta", 
                         BetaParam.a=2, BetaParam.b=3)
  expect_s4_class(result, "intData")
})

test_that("micro2intData accepts estimate.DistParam with General case", {
  MicroDt <- data.frame(X=1:9, Y=9:1)
  agrby <- factor(c("A", "A", "A", "B", "B", "B", "C", "C", "C"))
  
  result <- suppressWarnings(micro2intData(MicroDt, agrby,
                                          LatentCase="General",
                                          estimate.DistParam=TRUE))
  expect_s4_class(result, "intData")
})

test_that("micro2intData handles single variable data", {
  MicroDt <- data.frame(X=1:6)
  agrby <- factor(c("A", "A", "B", "B", "C", "C"))
  
  result <- micro2intData(MicroDt, agrby)
  expect_s4_class(result, "intData")
  expect_equal(ncol(result), 1)
  expect_equal(names(result), "X")
})

test_that("micro2intData handles many variables (3+)", {
  MicroDt <- data.frame(X=1:9, Y=9:1, Z=5:13)
  agrby <- factor(c("A", "A", "A", "B", "B", "B", "C", "C", "C"))
  
  result <- micro2intData(MicroDt, agrby)
  expect_s4_class(result, "intData")
  expect_equal(ncol(result), 3)
  expect_equal(names(result), c("X", "Y", "Z"))
})

test_that("micro2intData handles negative values", {
  MicroDt <- data.frame(X=1:9, Y=-9:-1, Z=5:13)
  agrby <- factor(c("A", "A", "A", "B", "B", "B", "C", "C", "C"))
  
  result <- micro2intData(MicroDt, agrby)
  expect_s4_class(result, "intData")
  expect_true(all(is.finite(as.matrix(result@Ranges))))
})

test_that("micro2intData handles mixed positive and negative values", {
  MicroDt <- data.frame(X=c(-5, 0, 5, -3, 2, 8), Y=c(1, -1, 3, -3, 2, 5))
  agrby <- factor(c("A", "A", "A", "B", "B", "B"))
  
  result <- micro2intData(MicroDt, agrby)
  expect_s4_class(result, "intData")
  expect_equal(nrow(result), 2)
})

test_that("micro2intData handles large numeric values", {
  MicroDt <- data.frame(X=c(1e6, 2e6, 3e6, 4e6), Y=c(1e5, 2e5, 3e5, 4e5))
  agrby <- factor(c("A", "A", "B", "B"))
  
  result <- micro2intData(MicroDt, agrby)
  expect_s4_class(result, "intData")
  expect_true(all(is.finite(as.matrix(result@Ranges))))
})

test_that("micro2intData handles very small numeric values", {
  MicroDt <- data.frame(X=c(1e-6, 2e-6, 3e-6, 4e-6), Y=c(1e-5, 2e-5, 3e-5, 4e-5))
  agrby <- factor(c("A", "A", "B", "B"))
  
  result <- micro2intData(MicroDt, agrby)
  expect_s4_class(result, "intData")
  expect_true(all(is.finite(as.matrix(result@Ranges))))
})

test_that("micro2intData preserves row names", {
  MicroDt <- data.frame(X=1:6, Y=6:1, row.names=c("r1","r2","r3","r4","r5","r6"))
  agrby <- factor(c("A", "A", "B", "B", "C", "C"))
  
  result <- micro2intData(MicroDt, agrby)
  expect_equal(result@ObsNames, c("A", "B", "C"))
})

test_that("micro2intData with quantiles correctly trims outliers", {
  MicroDt <- data.frame(X=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
  agrby <- factor(rep("A", 10))
  
  result <- micro2intData(MicroDt, agrby, agrcrt=c(0.2, 0.8))
  expect_s4_class(result, "intData")
  # With 0.2-0.8 quantiles on 1:10, should get [3,8]
  expect_true(result@Centers[1,1] == 5.5)
  expect_true(result@Ranges[1,1] == 5)
})

test_that("micro2intData with quantiles creates correct NbMicroUnits", {
  MicroDt <- data.frame(X=rep(1:10, 2), Y=rep(1:10, 2))
  agrby <- factor(rep(c("A", "B"), each=10))
  
  result <- micro2intData(MicroDt, agrby, agrcrt=c(0.1, 0.9))
  expect_equal(NbMicroUnits(result), c(A=10, B=10))
})
