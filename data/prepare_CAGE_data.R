library("dplyr")
library("wiggleplotr")
library("devtools")
library("tidyr")

#Import prepared transcript annotations
txrevise_data = readRDS("scripts/processed/Homo_sapiens.GRCh38.96.txrevise_annotations.rds")

#Import CAGE data
cage_transcripts = readRDS("data/new_transcripts.rds")

#Make transcript metadata for cage peaks
cage_metadata = dplyr::tibble(cage_id = names(cage_transcripts)) %>% 
  tidyr::separate(cage_id, c("ensembl_gene_id","dot1", "dot2", "index")) %>%
  dplyr::select(ensembl_gene_id) %>%
  dplyr::group_by(ensembl_gene_id) %>%
  dplyr::mutate(index = c(1:n())) %>%
  dplyr::mutate(ensembl_transcript_id = paste0(ensembl_gene_id, "_CAGE", index)) %>%
  dplyr::select(-index) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(longest_start = 0, longest_end = 0, cds_start_NF = 0, cds_end_NF = 1, cds_start_end_NF = 0)
cage_metadata = dplyr::semi_join(cage_metadata, txrevise_data$transcript_metadata, by = "ensembl_gene_id")

names(cage_transcripts) = cage_metadata$ensembl_transcript_id

#Make list with cage metadata
cage_list = list(exons = cage_transcripts, transcript_metadata = cage_metadata)
saveRDS(cage_list, "data/CAGE_promoter_annotations_240420.rds")
