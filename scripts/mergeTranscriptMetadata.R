library("GenomicFeatures")
library("biomaRt")
library("dplyr")

#Import trancript data from biomart
col_names = c("ensembl_gene_id", "ensembl_transcript_id", "chromosome_name", "gene_end", 
              "strand", "transcript_tsl","transcript_length","transcript_appris",
              "transcript_gencode_basic","percentage_gene_gc_content","gene_biotype","transcript_biotype",
              "transcript_start","transcript_end","gene_start", "external_gene_name", "gene_version", "transcript_version")
biomart_raw = readr::read_tsv("data-raw/Homo_sapiens.GRCh38.92.biomart_export.txt.gz", col_names = col_names, skip = 1)

#Import transcript tags
transcript_tags = readr::read_tsv("data-raw/Homo_sapiens.GRCh38.92.transcript_tags.txt.gz") %>%
  dplyr::rename(ensembl_transcript_id = transcript_id)

#Merge the two datasets
transcript_metadata = dplyr::left_join(biomart_raw, transcript_tags, by = "ensembl_transcript_id")
saveRDS(transcript_metadata, "data-raw/Homo_sapiens.GRCh38.92.transcript_metadata.rds")