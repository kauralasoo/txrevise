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
	expect_that(dplyr::filter(rassf1, type == "transcribed")$diff, is_equivalent_to(756)) #756 bp upstream change 
	expect_that(dplyr::filter(rassf1, type == "coding")$diff, is_equivalent_to( 504 )) #is coding
	})

test_that("ZC3HAV1 has upstream and downstream changes", {
  zc3hav1 = classifyDifference("ENST00000242351","ENST00000471652", new_annotations, cdss)
  expect_equivalent(dplyr::filter(zc3hav1, type == "transcribed")$diff, c(71,5385)) #Upstream change 
  expect_equivalent(dplyr::filter(zc3hav1, type == "coding")$diff, c(617)) #is coding
})

test_that("ZC3HAV1 modified transcript has only downsteam changes", {
  zc3hav1 = classifyDifference("ENST00000242351.2","ENST00000471652", new_annotations)
  expect_equivalent(zc3hav1$diff, 5385) #Upstream change 
})

test_that("RASSF5 has upstream and downstream changes, but only upstream is coding",{
  rassf5 = classifyDifference("ENST00000355294","ENST00000304534", new_annotations, cdss)
  tx_diff = dplyr::filter(rassf5, type == "transcribed") %>% 
    summarize(diff = sum(diff)) %>% 
    as.numeric()
  expect_equivalent(tx_diff, 2720)
  cds_diff = dplyr::filter(rassf5, type == "coding") %>% 
    summarize(diff = sum(diff)) %>% 
    as.numeric()
  expect_equivalent(cds_diff, 699)
})

test_that("OSBPL9 has upstream and coding changes", {
  osbpl9 = classifyDifference("ENST00000428468","ENST00000361556", new_annotations, cdss)
  tx_test = dplyr::filter(osbpl9, type == "transcribed") %>% select(diff)
  expect_equivalent(tx_test, data.frame(c(507,39)))
  cds_test = dplyr::filter(osbpl9, type == "coding") %>% select(diff)
  expect_equivalent(cds_test, data.frame(c(345,39))) 
})

test_that("OSBPL9 has upstream and coding changes", {
  osbpl9 = classifyDifference("ENST00000428468","ENST00000428468.2", new_annotations)
  expect_equivalent(osbpl9$diff, 507)
})


test_that("PNPT1 downstream and contained, but not coding",{
  pnpt1 = classifyDifference("ENST00000447944.2","ENST00000415374", new_annotations)
  expect_equivalent(pnpt1$diff, c(304,712))
})

test_that("SUN2 has all coding changes", {
  sun2 = classifyDifference("ENST00000216064","ENST00000469086", new_annotations, cdss)
  tx_test = dplyr::filter(sun2, type == "transcribed") %>% select(diff)
  expect_equivalent(tx_test, data.frame(c(1104,2029,574)))
  cds_test = dplyr::filter(sun2, type == "coding") %>% select(diff)
  expect_equivalent(cds_test, data.frame(c(887,492,574)))
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
  
  test1 = dplyr::filter(ncoa7, type == "transcribed", position == "upstream") %>% select(diff) %>% as.numeric()
  expect_equivalent(test1, 2648)
  test2 = dplyr::filter(ncoa7, type == "coding", position == "upstream") %>% select(diff) %>% as.numeric()
  expect_equivalent(test2, 2371)
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
  
  expect_equivalent(nrow(class1), 7)
  expect_equivalent(nrow(class2), 2)
})

test_that("Genes with no overlapping transcripts are not included.",{
  tx_list_3 = list(RMI2 = c("ENST00000572992","ENST00000312499"),
                   NCOA7 = c("ENST00000438495","ENST00000392477"))
  class3 = applyClassifyDifference(tx_list_3, new_annotations)
  expect_equivalent(unique(class3$gene_name), "NCOA7")
})



