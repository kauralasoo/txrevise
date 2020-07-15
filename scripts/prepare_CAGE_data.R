library("dplyr")
library("devtools")
library("tidyr")
library("optparse")

#Read command-line options
option_list <- list(
  make_option(c("--new_transcripts"), type="character", default=NULL,
              help="Path to new annotations made based on CAGE.", metavar = "path"),
  make_option(c("--txrevise_annotations"), type="character", default=NULL,
              help="Path to the annotations file made by prepareAnnotations.R.", metavar = "path"),
  make_option(c("--output"), type="character", default=NULL,
              help="Where to output metadata.", metavar = "path")
)
opt <- parse_args(OptionParser(option_list=option_list))

#Import prepared transcript annotations
txrevise_data = readRDS(opt$txrevise_annotations)

#Import CAGE data
cage_transcripts = readRDS(opt$new_transcripts)

#Make transcript metadata for cage peaks
cage_metadata = dplyr::tibble(cage_id = names(cage_transcripts)) %>%
  tidyr::separate(cage_id, c("ensembl_gene_id","dot1", "dot2", "index")) %>%
  dplyr::select(ensembl_gene_id) %>%
  dplyr::group_by(ensembl_gene_id) %>%
  dplyr::summarise(transcript_number = c(1:n())) %>%
  dplyr::mutate(ensembl_transcript_id = paste0(ensembl_gene_id, "_CAGE", transcript_number)) %>%
  dplyr::select(-transcript_number) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(longest_start = 0, longest_end = 0, cds_start_NF = 0, cds_end_NF = 1, cds_start_end_NF = 0)
cage_metadata = dplyr::semi_join(cage_metadata, txrevise_data$transcript_metadata, by = "ensembl_gene_id")

names(cage_transcripts) = cage_metadata$ensembl_transcript_id

#Make list with cage metadata
cage_list = list(exons = cage_transcripts, transcript_metadata = cage_metadata)
saveRDS(cage_list, opt$output)
