#New Annotations file
exons = readRDS("../data/exons_test.rds")
cdss = readRDS("../data/cdss_test.rds")
new_annotations = readRDS("../data/treatment_res_annot.rds")

context("Expand two transcripts into a set with only unit changes")

test_that("Revsing RASSF1 generates no new transcripts", {
  rassf1 = constructAlternativeTranscripts("ENST00000327761", "ENST00000359365",exons)
  expect_equal(length(rassf1),2) #Two transcripts
  expect_true(all(names(rassf1) %in% c("ENST00000327761", "ENST00000359365"))) #Names unchanges
})

test_that("Alternative promoter vs exon (OSBPL9)",{
  txs = c("ENST00000428468","ENST00000361556")
  new_tx = constructAlternativeTranscripts(txs[1], txs[2],exons)
  expect_equal(length(new_tx),4) #New annotations have four transcripts
  
  #Five exon different between the initial two transcripts
  expect_equal(sum(width(setdiff(new_tx$ENST00000428468, new_tx$ENST00000361556))),359)
  expect_equal(length(width(setdiff(new_tx$ENST00000428468, new_tx$ENST00000361556))),5)
  
  #Alternative exon removed in the new transcript
  expect_equal(width(setdiff(new_tx$ENST00000428468, new_tx$ENST00000428468.3)),39)
  
  #Four exons removed from the first alternative transcript
  expect_equal(length(setdiff(new_tx$ENST00000428468, new_tx$ENST00000428468.2)),4)
  expect_equal(sum(width(setdiff(new_tx$ENST00000428468, new_tx$ENST00000428468.2))),320)
})

test_that("RMI2 has two non-overlapping transcripts. Should issue a warning.",{
  rmi2_txs = c("ENST00000572992","ENST00000312499")
  new_tx = constructAlternativeTranscripts(rmi2_txs[1], rmi2_txs[2],exons)
  warn = tryCatch(constructAlternativeTranscripts(rmi2_txs[1], rmi2_txs[2],exons), warning = function(w) w)
  expect_equal(warn$message,"No shared exons between two transcripts.")
  expect_equal(length(new_tx),2)
})

context("Classify differences between two transcripts")

test_that("RASSF1 has only upstream coding changes", {
	rassf1 = classifyDifference("ENST00000327761", "ENST00000359365", new_annotations, cdss)
	expect_that(rassf1$transcribe$diff, is_equivalent_to( c(756,0,0) )) #756 bp upstream change 
	expect_that(rassf1$coding$diff, is_equivalent_to( c(504,0,0) )) #is coding
	})

test_that("ZC3HAV1 has upstream and downstream changes", {
  zc3hav1 = classifyDifference("ENST00000242351","ENST00000471652", new_annotations, cdss)
  expect_equivalent(zc3hav1$transcribed$diff, c(71,5385,0)) #Upstream change 
  expect_equivalent(zc3hav1$coding$diff, c(0,617,0)) #is coding
})

test_that("ZC3HAV1 modified transcript has only downsteam changes", {
  zc3hav1 = classifyDifference("ENST00000242351.2","ENST00000471652", new_annotations)
  expect_equivalent(zc3hav1$transcribed$diff, c(0,5385,0)) #Upstream change 
})

test_that("RASSF5 has upstream and downstream changes, but only upstream is coding",{
  rassf5 = classifyDifference("ENST00000355294","ENST00000304534", new_annotations, cdss)
  expect_equivalent(rassf5$transcribed$diff, c(1165,1555,0))
  expect_equivalent(rassf5$coding$diff, c(699,0,0)) 
})

test_that("OSBPL9 has upstream and coding changes", {
  osbpl9 = classifyDifference("ENST00000428468","ENST00000361556", new_annotations, cdss)
  expect_equivalent(osbpl9$transcribed$diff, c(507,0,39))
  expect_equivalent(osbpl9$coding$diff, c(345,0,39)) 
})

test_that("OSBPL9 has upstream and coding changes", {
  osbpl9 = classifyDifference("ENST00000428468","ENST00000428468.2", new_annotations)
  expect_equivalent(osbpl9$transcribed$diff, c(507,0,0))
})


test_that("PNPT1 downstream and contained, but not coding",{
  pnpt1 = classifyDifference("ENST00000447944.2","ENST00000415374", new_annotations)
  expect_equivalent(pnpt1$transcribed$diff, c(0,304,712))
})

test_that("SUN2 has all coding changes", {
  sun2 = classifyDifference("ENST00000216064","ENST00000469086", new_annotations, cdss)
  expect_equivalent(sun2$transcribed$diff, c(1104,2029,574))
  expect_equivalent(sun2$coding$diff, c(887,492,574))
})

test_that("RMI2 has no overlapping transcripts so the comparison does not make semse",{
  class = classifyDifference("ENST00000572992","ENST00000312499", new_annotations, cdss)
  expect_equal(class, NULL)
  warn = tryCatch({x = classifyDifference("ENST00000572992","ENST00000312499", new_annotations, cdss)}, 
                  warning = function(w){return(w)})
  expect_equal(class(warn)[1], "simpleWarning")
})

test_that("NCOA7 - changes both upstreatm and downstream",{
  ncoa7 = classifyDifference("ENST00000438495","ENST00000392477", new_annotations, cdss)
  expect_equivalent(ncoa7$transcribed$diff, c(2648,9,0))
  expect_equivalent(ncoa7$coding$diff, c(2371,0,0))
})

context("Apply classifyDifference to more than two txs")

test_that("NCOA7 - changes both upstreatm and downstream",{
  #Run with old annotations
  tx_list = list(NCOA7 = c("ENST00000438495","ENST00000392477"), OSBPL9 = c("ENST00000428468","ENST00000361556"))
  class1 = applyClassifyDifference(tx_list, new_annotations, cdss)
  
  #Run with new annotations
  tx_list_new = list(NCOA7 = c("ENST00000438495","ENST00000438495.2"), 
                     OSBPL9 = c("ENST00000428468","ENST00000428468.2"))
  class2 = applyClassifyDifference(tx_list_new, new_annotations)
  
  coding_sum1 = colSums(class1$coding$diff * sign(class1$transcribed$diff)) 
  coding_sum2 = colSums(class1$coding$diff * sign(class2$transcribed$diff)) 
  
  expect_equivalent(coding_sum1, c(2716,0,39))
  expect_equivalent(coding_sum2, c(2716,0,0))
})



