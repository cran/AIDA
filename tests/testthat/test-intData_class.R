library(testthat)

test_that("intData constructor and accessors", {
	Data <- data.frame(L1 = c(1,2,3), U1 = c(2,3,4), L2 = c(0,1,2), U2 = c(1,2,3))
	obj <- intData(Data, Seq = "LbUb_VarbyVar", VarNames = c("A","B"))
obj[1,1]
	expect_s4_class(obj, "intData")
	expect_equal(nrow(obj), 3)
	expect_equal(ncol(obj), 2)
	expect_equal(colnames(obj), c("A","B"))
	expect_equal(rownames(obj), as.character(1:3))

	expected_centers <- data.frame(A.Centers = c(1.5,2.5,3.5), B.Centers = c(0.5,1.5,2.5), row.names = as.character(1:3), check.names = FALSE)
	expect_equal(Centers(obj), expected_centers)

	expected_ranges <- data.frame(A.Ranges = c(1,1,1), B.Ranges = c(1,1,1), row.names = as.character(1:3), check.names = FALSE)
	expect_equal(Ranges(obj), expected_ranges)

	expected_lb <- data.frame(A.Lbnd = c(1,2,3), B.Lbnd = c(0,1,2), row.names = as.character(1:3), check.names = FALSE)
	expect_equal(LowerBounds(obj), expected_lb)

	expected_ub <- data.frame(A.Ubnd = c(2,3,4), B.Ubnd = c(1,2,3), row.names = as.character(1:3), check.names = FALSE)
	expect_equal(UpperBounds(obj), expected_ub)
})

test_that("intData constructor errors on invalid Data inputs", {
  expect_error(intData(1:3), "must be a data frame or a matrix")

  bad <- data.frame(a = 1:3, b = 2:4, c = 3:5)
  expect_error(intData(bad), "must be an even number")

  bad2 <- data.frame(L1 = c(1,2), U1 = c(2, NA), L2 = c(0,1), U2 = c(1,2))
  expect_error(intData(bad2, Seq = "LbUb_VarbyVar"), "Invalid data")
})

test_that("intData reads correctly a data frame of lower and upper bounds presented variable by variable", {

  Xbnds <- data.frame(XLB=c(1,2,4),XUB=c(9,7,6),row.names=paste("Unit",1:3))
  Ybnds <- data.frame(YLB=c(10,3,5),YUB=c(15,8,6),row.names=paste("Unit",1:3))
  Allbnds <- cbind(Xbnds,Ybnds)

  Idt <- intData(Allbnds,VarNames = c("X","Y"),Seq="LbUb_VarbyVar")
  
  expect_equal(nrow(Idt),3)
  expect_equal(ncol(Idt),2)
  expect_equal(row.names(Idt),row.names(Allbnds))
  expect_equal(rownames(Idt),rownames(Allbnds))
  
  X.Centers <- (Xbnds$XLB + Xbnds$XUB) /2 
  Y.Centers <- (Ybnds$YLB + Ybnds$YUB) /2 
  AllCenters <- as.data.frame(cbind(X.Centers,Y.Centers))
  row.names(AllCenters) <- row.names(Allbnds)

  expect_identical(Centers(Idt),AllCenters)
    
  X.Ranges <- Xbnds$XUB - Xbnds$XLB 
  Y.Ranges <- Ybnds$YUB - Ybnds$YLB
  AllRanges <- as.data.frame(cbind(X.Ranges,Y.Ranges))
  row.names(AllRanges) <- row.names(Allbnds)
  
  expect_identical(Ranges(Idt),AllRanges)
  
} )

test_that("intData reads correctly a data frame of all lower bounds followed by all upper bounds", {
  
  Lbnds <- data.frame(XLB=c(1,2,4),YLB=c(10,3,5),row.names=paste("Unit",1:3))
  Ubnds <- data.frame(XUB=c(9,7,6),YUB=c(15,8,6),row.names=paste("Unit",1:3))
  Allbnds <- cbind(Lbnds,Ubnds)
  
  Idt <- intData(Allbnds,Seq="AllLb_AllUb",VarNames = c("X","Y"))
  
  expect_equal(nrow(Idt),3)
  expect_equal(ncol(Idt),2)
  expect_equal(row.names(Idt),row.names(Allbnds))
  expect_equal(rownames(Idt),rownames(Allbnds))
  
  X.Centers <- (Lbnds$XLB + Ubnds$XUB) /2 
  Y.Centers <- (Lbnds$YLB + Ubnds$YUB) /2 
  AllCenters <- as.data.frame(cbind(X.Centers,Y.Centers))
  row.names(AllCenters) <- row.names(Allbnds)
  
  expect_identical(Centers(Idt),AllCenters)
  
  X.Ranges <- Ubnds$XUB - Lbnds$XLB
  Y.Ranges <- Ubnds$YUB - Lbnds$YLB
  AllRanges <- as.data.frame(cbind(X.Ranges,Y.Ranges))
  row.names(AllRanges) <- row.names(Allbnds)
  
  expect_identical(Ranges(Idt),AllRanges)
  
} )

test_that("intData reads correctly a data frame of Centersoints and Rangesanges presented variable by variable", {
  
  Xbnds <- data.frame(XLB=c(1,2,4),XUB=c(9,7,6),row.names=paste("Unit",1:3))
  Ybnds <- data.frame(YLB=c(10,3,5),YUB=c(15,8,6),row.names=paste("Unit",1:3))
  X.Centers <- (Xbnds$XLB + Xbnds$XUB) /2 
  Y.Centers <- (Ybnds$YLB + Ybnds$YUB) /2 
  X.Ranges <- Xbnds$XUB - Xbnds$XLB
  Y.Ranges <- Ybnds$YUB - Ybnds$YLB
  RangesDF <- data.frame(cbind(X.Centers,X.Ranges,Y.Centers,Y.Ranges),row.names=rownames(Xbnds))
  
  Idt <- intData(RangesDF,Seq="CenRng_VarbyVar",VarNames = c("X","Y"))
  
  expect_equal(nrow(Idt),3)
  expect_equal(ncol(Idt),2)
  expect_equal(row.names(Idt),row.names(Xbnds))
  expect_equal(rownames(Idt),rownames(Xbnds))
  
  AllCenters <- as.data.frame(cbind(X.Centers,Y.Centers))
  row.names(AllCenters) <- row.names(Xbnds)
  expect_identical(Centers(Idt),AllCenters)
  
  AllRanges <- as.data.frame(cbind(X.Ranges,Y.Ranges))
  row.names(AllRanges) <- row.names(Xbnds)
  expect_identical(Ranges(Idt),AllRanges)
  
} )

test_that("intData reads correctly a data frame of all Centers followed by all Ranges", {
  
  Xbnds <- data.frame(XLB=c(1,2,4),XUB=c(9,7,6),row.names=paste("Unit",1:3))
  Ybnds <- data.frame(YLB=c(10,3,5),YUB=c(15,8,6),row.names=paste("Unit",1:3))
  X.Centers <- (Xbnds$XLB + Xbnds$XUB) /2 
  Y.Centers <- (Ybnds$YLB + Ybnds$YUB) /2 
  X.Ranges <- Xbnds$XUB - Xbnds$XLB
  Y.Ranges <- Ybnds$YUB - Ybnds$YLB 

  AllCenters <- as.data.frame(cbind(X.Centers,Y.Centers))
  row.names(AllCenters) <- row.names(Xbnds)
  AllRanges <- as.data.frame(cbind(X.Ranges,Y.Ranges))
  row.names(AllRanges) <- row.names(Xbnds)
  RangesDF <- cbind(AllCenters,AllRanges)
  
  Idt <- intData(RangesDF,Seq="AllCen_AllRng",VarNames = c("X","Y"))
  
  expect_equal(nrow(Idt),3)
  expect_equal(ncol(Idt),2)
  expect_equal(row.names(Idt),row.names(Xbnds))
  expect_equal(rownames(Idt),rownames(Xbnds))
  
  expect_identical(Centers(Idt),AllCenters)
  expect_identical(Ranges(Idt),AllRanges)
  
} )

test_that("rbind works correctly with two arguments", {
  
  Idt1 <- intData(data.frame(list(c(1,2,4),c(9,7,6),c(1,3,4),c(9,8,6))),VarNames=c("X","Y"), ObsNames = c("A","B","C"))
  Idt2 <- intData(data.frame(list(c(3,8,4),c(5,9,7)),c(2,6,3),c(7,8,4)),VarNames=c("X","Y"), ObsNames = c("D","E","F"))
  Idt12 <- intData(data.frame(list(c(1,2,4,3,8,4),c(9,7,6,5,9,7),c(1,3,4,2,6,3),c(9,8,6,7,8,4))),
                 VarNames=c("X","Y"),ObsNames = c("A","B","C","D","E","F"))
  
  expect_identical(rbind(Idt1,Idt2),Idt12)
  
  
  Idt1 <- intData(data.frame(list(c(1,2,4),c(9,7,6),c(1,3,4),c(9,8,6))),VarNames=c("X","Y"))
  Idt2 <- intData(data.frame(list(c(3,8,4),c(5,9,7)),c(2,6,3),c(7,8,4)),VarNames=c("X","Y"))
  Idt12 <- intData(data.frame(list(c(1,2,4,3,8,4),c(9,7,6,5,9,7),c(1,3,4,2,6,3),c(9,8,6,7,8,4))),
                 VarNames=c("X","Y"),ObsNames = as.character(c(1,2,3,11,21,31)))
  
  expect_identical(rbind(Idt1,Idt2),Idt12)
} )

test_that("rbind works correctly with more than two arguments", {
  
  Idt1 <- intData(data.frame(list(c(1,2,4),c(9,7,6),c(1,3,4),c(9,8,6))),VarNames=c("X","Y"), ObsNames = c("A","B","C"))
  Idt2 <- intData(data.frame(list(c(3,8,4),c(5,9,7)),c(2,6,3),c(7,8,4)),VarNames=c("X","Y"), ObsNames = c("D","E","F"))
  Idt3 <- intData(data.frame(list(2,3,6,8)),VarNames=c("X","Y"), ObsNames = "G")
  Idt123 <- intData(data.frame(list(c(1,2,4,3,8,4,2),c(9,7,6,5,9,7,3),c(1,3,4,2,6,3,6),c(9,8,6,7,8,4,8))),
                 VarNames=c("X","Y"),ObsNames = c("A","B","C","D","E","F","G"))
  
  expect_identical(rbind(Idt1,Idt2,Idt3),Idt123)
  
} )

test_that("row-only indexing works correctly for intData objects", {

  Xbnds <- data.frame(XLB=c(1,2,4),XUB=c(9,7,6),row.names=paste("Unit",1:3))
  Ybnds <- data.frame(YLB=c(10,3,5),YUB=c(15,8,6),row.names=paste("Unit",1:3))
  Allbnds <- cbind(Xbnds,Ybnds)
  Idt <- intData(Allbnds,VarNames = c("X","Y"))

  Idt2 <- intData(Allbnds[2,],VarNames = c("X","Y"))

  Idt13 <- intData(Allbnds[c(1,3),],VarNames = c("X","Y"))
  Idt31 <- intData(Allbnds[c(3,1),],VarNames = c("X","Y"))
  Idt23 <- intData(Allbnds[2:3,],VarNames = c("X","Y"))
  Idt32 <- intData(Allbnds[3:2,],VarNames = c("X","Y"))
  
  expect_identical(Idt[2,],Idt2)
  expect_identical(Idt[-c(1,3),],Idt2)
  expect_identical(Idt[c(1,3),],Idt13)
  expect_identical(Idt[-2,],Idt13)
  expect_identical(Idt[c(3,1),],Idt31)
  expect_identical(Idt[2:3,],Idt23)
  expect_identical(Idt[-1,],Idt23)
  expect_identical(Idt[3:2,],Idt32)

} )

test_that("column-only indexing works correctly for intData objects", {
  
  Lbnds <- data.frame(XLB=c(1,2,4),YLB=c(10,3,5),ZLB=c(8,5,4),row.names=paste("Unit",1:3))
  Ubnds <- data.frame(XUB=c(9,7,6),YUB=c(15,8,6),ZUB=c(12,9,7),row.names=paste("Unit",1:3))
  Idt <- intData(cbind(Lbnds,Ubnds),Seq="AllLb_AllUb",VarNames = c("X","Y","Z"))
  
  Idt2 <- intData(cbind(Lbnds[,2,drop=FALSE],Ubnds[,2,drop=FALSE]),Seq="AllLb_AllUb",VarNames = "Y")
  Idt13 <- intData(cbind(Lbnds[,c(1,3)],Ubnds[,c(1,3)]),Seq="AllLb_AllUb",VarNames = c("X","Z"))
  Idt31 <- intData(cbind(Lbnds[,c(3,1)],Ubnds[,c(3,1)]),Seq="AllLb_AllUb",VarNames = c("Z","X"))
  Idt23 <- intData(cbind(Lbnds[,2:3],Ubnds[,2:3]),Seq="AllLb_AllUb",VarNames = c("Y","Z"))
  Idt32 <- intData(cbind(Lbnds[,3:2],Ubnds[,3:2]),Seq="AllLb_AllUb",VarNames = c("Z","Y"))
  
  expect_identical(Idt[,2],Idt2)
  expect_identical(Idt[,-c(1,3)],Idt2)
  expect_identical(Idt[,c(1,3)],Idt13)
  expect_identical(Idt[,-2],Idt13)
  expect_identical(Idt[,c(3,1)],Idt31)
  expect_identical(Idt[,2:3],Idt23)
  expect_identical(Idt[,-1],Idt23)
  expect_identical(Idt[,3:2],Idt32)
  
} )

test_that("indexing by both rows and columns works correctly for intData objects", {
  
  Lbnds <- data.frame(XLB=c(1,2,4),YLB=c(10,3,5),ZLB=c(8,5,4),row.names=paste("Unit",1:3))
  Ubnds <- data.frame(XUB=c(9,7,6),YUB=c(15,8,6),ZUB=c(12,9,7),row.names=paste("Unit",1:3))
  Idt <- intData(cbind(Lbnds,Ubnds),Seq="AllLb_AllUb",VarNames = c("X","Y","Z"))
  
  Idtr2c3 <- intData(cbind(Lbnds[2,3,drop=FALSE],Ubnds[2,3,drop=FALSE]),Seq="AllLb_AllUb",VarNames = "Z")
  Idtr2c13 <- intData(cbind(Lbnds[2,c(1,3),drop=FALSE],Ubnds[2,c(1,3)]),Seq="AllLb_AllUb",VarNames = c("X","Z"))
  Idtr13c2 <- intData(cbind(Lbnds[c(1,3),2,drop=FALSE],Ubnds[c(1,3),2]),Seq="AllLb_AllUb",VarNames = "Y")
  Idtr13c21 <- intData(cbind(Lbnds[c(1,3),2:1],Ubnds[c(1,3),2:1]),Seq="AllLb_AllUb",VarNames = c("Y","X"))
  Idtr32c12 <- intData(cbind(Lbnds[3:2,1:2],Ubnds[3:2,1:2]),Seq="AllLb_AllUb",VarNames = c("X","Y"))
  
  expect_identical(Idt[2,3],Idtr2c3)
  expect_identical(Idt[-c(1,3),3],Idtr2c3)
  expect_identical(Idt[2,-(1:2)],Idtr2c3)
  expect_identical(Idt[-c(1,3),-c(1,2)],Idtr2c3)

  expect_identical(Idt[2,c(1,3)],Idtr2c13)
  expect_identical(Idt[-c(1,3),c(1,3)],Idtr2c13)
  expect_identical(Idt[2,-2],Idtr2c13)
  expect_identical(Idt[-c(1,3),-2],Idtr2c13)
  
  expect_identical(Idt[c(1,3),2],Idtr13c2)
  expect_identical(Idt[-2,2],Idtr13c2)
  expect_identical(Idt[c(1,3),-c(1,3)],Idtr13c2)
  expect_identical(Idt[-2,-c(1,3)],Idtr13c2)
  
  expect_identical(Idt[c(1,3),2:1],Idtr13c21)
  expect_identical(Idt[-2,2:1],Idtr13c21)

  expect_identical(Idt[3:2,1:2],Idtr32c12)
  expect_identical(Idt[3:2,-3],Idtr32c12)

} )

test_that("intData LatentParam validation", {
  Data <- data.frame(L1 = c(1,2), U1 = c(2,3))
  expect_error(intData(Data, LatentParam = 1, LatentDist = "Unif"), "`LatentParam` must be either")
  expect_error(intData(Data, LatentParam = list(1)), "Error: If LatentParam is provided, LatentDist must also be provided.")
  obj <- intData(Data, LatentParam = list(1), LatentDist = "Unif")
  expect_equal(LatentCase(obj), "U_id_symmetric")
})

test_that("NbMicroUnits handling", {
  Data <- data.frame(L1 = 1:3, U1 = 2:4)
  obj <- intData(Data, VarNames = "X", NbMicroUnits = 3)
  expect_identical(NbMicroUnits(obj), as.integer(3))
  obj2 <- intData(Data, VarNames = "X")
  expect_null(NbMicroUnits(obj2))
})

test_that("bounds and log ranges", {
  Data <- data.frame(L1 = c(1,2), U1 = c(3,4))
  obj <- intData(Data, VarNames = "X")
  expect_equal(names(LowerBounds(obj)), "X.Lbnd")
  expect_equal(names(UpperBounds(obj)), "X.Ubnd")
  expect_equal(as.numeric(LogRanges(obj)[,1]), as.numeric(log(Ranges(obj)[,1])))
})

test_that("head and tail behavior", {
  Data <- data.frame(L1 = 1:5, U1 = 2:6)
  obj <- intData(Data, VarNames = "X")
  expect_identical(head(obj, n = 2), obj[1:2, ])
  expect_identical(tail(obj, n = 1), obj[nrow(obj), ])
})

test_that("plot error conditions for too many variables", {
  Data <- data.frame(L1 = 1:3, U1 = 2:4, L2 = 1:3, U2 = 2:4, L3 = 1:3, U3 = 2:4)
  obj3 <- intData(Data, VarNames = c("X","Y","Z"))
  expect_error(plot(obj3), "Currently method plot can only plot at most two interval variables")

  obj2 <- intData(Data[,1:4], VarNames = c("X","Y"))
  expect_error(plot(obj2, obj2), "Currently intData method plot can plot only one integer variable on the horizontal axis")
})

test_that("summary and show output", {
  Data <- data.frame(L1 = 1:2, U1 = 2:3)
  obj <- intData(Data, VarNames = "X", ObsNames = c("A","B"))
  expect_output(summary(obj), "Centers summary:")
  expect_output(show(obj), "A")
})

test_that("equality comparison and errors", {
  Data <- data.frame(L1 = 1:2, U1 = 2:3)
  a <- intData(Data, VarNames = "X")
  b <- intData(Data, VarNames = "X")
  eq <- a == b
  expect_true(all(eq))
  expect_error(a == 1, "comparison")
  c <- intData(data.frame(L1 = 1:3, U1 = 2:4), VarNames = "X")
  expect_error(a == c, "== only defined for equally-sized intData objects")
})

test_that("indexing by character and numeric names", {
  Data <- data.frame(L1 = c(1,2,3), U1 = c(2,3,4))
  obj <- intData(Data, VarNames = "X", ObsNames = c("one","two","three"))
  expect_identical(obj["two", ], obj[2, ])
  expect_identical(obj[, "X"], obj[,1])
})

test_that("General latent case subsetting preserves LatentParam matrices", {
  Data <- data.frame(L1 = 1:3, U1 = 2:4, L2 = 2:4, U2 = 3:5)
  mat1 <- matrix(c(1,0,0,1),2,2)
  mat2 <- matrix(c(2,0,0,2),2,2)
  obj <- intData(Data, VarNames = c("A","B"), LatentParam = list(mat1, mat2), LatentCase = "General", LatentDist = c("Unif","Triang"))
  sub <- obj[, "A"]
  lp <- LatentParam(sub)
  expect_true(is.list(lp) && is.matrix(lp[[1]]) && dim(lp[[1]])[1]==1)
  expect_equal(as.numeric(lp[[1]]), as.numeric(mat1[1,1]))
  expect_equal(LatentDist(sub), "Unif")
})

test_that("plotting branches run without error (two-object and single-object)", {
  Data <- data.frame(L1 = c(1,2,3), U1 = c(2,3,4))
  x <- intData(Data, VarNames = "X")
  y <- intData(Data, VarNames = "Y")
  # two-object plots (NIVar must be 1)
  png(filename = tempfile(fileext = ".png"))
  expect_silent(plot(x, x, type = "crosses"))
  expect_silent(plot(x, x, type = "rectangles"))
  expect_silent(plot(x, x, type = "crosses2"))
  dev.off()

  # single-object plots vertical/horizontal with append
  png(filename = tempfile(fileext = ".png"))
  expect_silent(plot(x, layout = "vertical"))
  expect_silent(plot(x, layout = "horizontal", append = TRUE))
  dev.off()
})

test_that("rbind errors on mismatched properties", {
  Data1 <- data.frame(L1 = 1:2, U1 = 2:3)
  Data2 <- data.frame(L1 = 1:3, U1 = 2:4, L2 = 8:6, U2 = 11:13)
  a <- intData(Data1, VarNames = "X", LatentParam = list(1), LatentDist = "Unif")
  a@ObsNames <- c("a.1","a.2")
  b <- intData(Data2, VarNames = "X", LatentParam = list(1), LatentDist = "Unif")
  expect_error(rbind(a, b), "Arguments x and y have a different number of interval-valued variables|== only defined for equally-sized intData objects")

  # different LatentDist
  c <- intData(Data1, VarNames = "X", LatentParam = list(1), LatentDist = "Triang")
  expect_error(rbind(a, c), "Arguments x and y have different LatentDist")
})

test_that("basic accessors and names/dim methods", {
  Data <- data.frame(L1 = 1:4, U1 = 2:5, L2 = 3:6, U2 = 4:7)
  obj <- intData(Data, VarNames = c("X","Y"), ObsNames = c("a","b","c","d"), Seq = "LbUb_VarbyVar")
  expect_equal(nrow(obj), 4)
  expect_equal(ncol(obj), 2)
  expect_equal(dim(obj), c(4,2))
  expect_equal(rownames(obj), c("a","b","c","d"))
  expect_equal(colnames(obj), c("X","Y"))
  expect_equal(names(obj), c("X","Y"))
  expect_equal(Centers(obj)[1,1], (1+2)/2)
})
