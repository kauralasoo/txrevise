exons = readRDS("../data/exons.rds")
cdss = readRDS("../data/cdss.rds")
new_annotations = readRDS("../data/treatment_res_annot.rds")

test_transcripts = c("ENST00000327761","ENST00000359365",
                     "ENST00000428468","ENST00000361556",
                     "ENST00000572992","ENST00000312499",
                     "ENST00000327761","ENST00000359365",
                     "ENST00000242351","ENST00000471652", #ZC3HAV1
                     "ENST00000242351.2","ENST00000471652",#ZC3HAV1
                     "ENST00000355294","ENST00000304534", #RASSF5
                     "ENST00000428468","ENST00000428468.2", #OSBPL9
                     "ENST00000428468","ENST00000361556", #OSBPL9 
                     "ENST00000447944.2","ENST00000415374",
                     "ENST00000216064","ENST00000469086",
                     "ENST00000572992","ENST00000312499",
                     "ENST00000438495","ENST00000392477",
                     "ENST00000438495","ENST00000438495.2")
test_transcripts = unique(test_transcripts)
#Filter annotations
exons = exons[intersect(test_transcripts, names(exons))]
cdss = cdss[intersect(test_transcripts, names(exons))]
print(length(exons))
saveRDS(exons, "../data/exons_test.rds")
saveRDS(cdss, "../data/cdss_test.rds")