
#New Annotations file
#exons = readRDS("../data/exons.rds")
cdss = readRDS("../data/cdss.rds")
new_annotations = readRDS("../data/treatment_res_annot.rds")

context("Classify differences between two transcripts")

test_that("RASSF1 has only upstream coding changes", {
	rassf1 = classifyDifference("ENST00000327761", "ENST00000359365", new_annotations, cdss)
	expect_that(rassf1$difference, is_equivalent_to( c(756,0,0) )) #756 bp upstream change 
	expect_that(rassf1$coding, is_equivalent_to( c(1,0,0) )) #is coding
	})

test_that("ZC3HAV1 has downstream coding changes", {
  zc3hav1 = classifyDifference("ENST00000242351.2","ENST00000471652", new_annotations, cdss)
  expect_equivalent(zc3hav1$difference, c(0,5385,0)) #Upstream change 
  expect_equivalent(zc3hav1$coding, c(0,1,0)) #is coding
})

test_that("RASSF5 has upstream and downstream changes, but only upstream is coding",{
  rassf5 = classifyDifference("ENST00000355294","ENST00000304534", new_annotations, cdss)
  expect_equivalent(rassf5$difference, c(1165,1555,0))
  expect_equivalent(rassf5$coding, c(1,0,0)) 
})


test_that("OSBPL9 has upstream and coding changes", {
  osbpl9 = classifyDifference("ENST00000428468","ENST00000428468.2", new_annotations, cdss)
  expect_equivalent(osbpl9$difference, c(507,0,0))
  expect_equivalent(osbpl9$coding, c(1,0,0)) 
})


test_that("PNPT1 downstream and contained, but not coding",{
  pnpt1 = classifyDifference("ENST00000447944.2","ENST00000415374", new_annotations, cdss)
  expect_equivalent(pnpt1$difference, c(0,304,712))
  expect_equivalent(pnpt1$coding, c(0,0,0)) 
})

test_that("SUN2 has all coding changes", {
  sun2 = classifyDifference("ENST00000216064","ENST00000469086", new_annotations, cdss)
  expect_equivalent(sun2$difference, c(1104,2029,574))
  expect_equivalent(sun2$coding, c(1,1,1))
})

test_that("RMI2 has no overlapping transcripts so the comparison does not make semse",{
  class = classifyDifference("ENST00000572992","ENST00000312499", new_annotations, cdss)
  expect_equal(class, NULL)
  warn = tryCatch({x = classifyDifference("ENST00000572992","ENST00000312499", new_annotations, cdss)}, warning = function(w){return(w)})
  expect_equal(class(warn)[1], "simpleWarning")
})

test_that("NCOA7 - changes both upstreatm and downstream",{
  ncoa7 = classifyDifference("ENST00000438495","ENST00000392477", new_annotations, cdss)
  expect_equivalent(ncoa7$difference, c(2648,9,0))
  expect_equivalent(ncoa7$coding, c(1,0,0))
})

context("Expand two transcripts into a set with only unit changes")




