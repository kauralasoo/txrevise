library("GenomicFeatures")
library("biomaRt")
library("dplyr")
library("devtools")
load_all("../txrevise")
load_all("../seqUtils")

#Import transcript tags
transcript_tags = readr::read_tsv("data-raw/Homo_sapiens.GRCh38.92.transcript_tags.txt.gz") %>%
  dplyr::rename(ensembl_transcript_id = transcript_id)
transcript_meta = txrevise::importTranscriptMetadataFromGTF("data-raw/Homo_sapiens.GRCh38.92.gtf.gz", transcript_tags)

#Filter the metadata
filtered_metadata = txrevise::filterTranscriptMetadata(transcript_meta)

#Construct TxDb
txdb = GenomicFeatures::makeTxDbFromGFF("data-raw/Homo_sapiens.GRCh38.92.gtf.gz")
exons = GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
cdss = GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)

#Export data for txrevise event construction
txrevise_data = list(transcript_metadata = filtered_metadata, 
                     exons = exons[filtered_metadata$ensembl_transcript_id], 
                     cdss = cdss[intersect(names(cdss),filtered_metadata$ensembl_transcript_id)])
saveRDS(txrevise_data, "data-raw/Homo_sapiens.GRCh38.92.txrevise_data.rds")

#Make events  
events = txrevise::constructAlternativeEventsWrapper("ENSG00000128604", txrevise_data$transcript_metadata, 
                                                  txrevise_data$exons, txrevise_data$cdss)

