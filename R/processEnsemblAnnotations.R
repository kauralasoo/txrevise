
#' Extract transcript metadata from Ensembl GTF file
#'
#' @param gtf_path Path to the GTF file
#' @param transcript_tags tibble of transcript tags
#'
#' @return A tibble with transcript metadata.
#' @export
importTranscriptMetadataFromGTF <- function(gtf_path, transcript_tags){
  
  #Check that transcript tags tibble has all of the right columns
  assertthat::assert_that(assertthat::has_name(transcript_tags, "ensembl_transcript_id"))
  assertthat::assert_that(assertthat::has_name(transcript_tags, "CCDS"))
  assertthat::assert_that(assertthat::has_name(transcript_tags, "basic"))
  assertthat::assert_that(assertthat::has_name(transcript_tags, "cds_start_NF"))
  assertthat::assert_that(assertthat::has_name(transcript_tags, "mRNA_start_NF"))
  assertthat::assert_that(assertthat::has_name(transcript_tags, "cds_end_NF"))
  assertthat::assert_that(assertthat::has_name(transcript_tags, "mRNA_end_NF"))
  
  #Import the GTF file
  gtf_file = rtracklayer::import(gtf_path, format = "gtf")
  
  #Calculate transcript lengths (sum of exons)
  transcript_lengths = gtf_file[gtf_file$type == "exon",] %>% as.data.frame() %>% as_tibble() %>%
    dplyr::select(transcript_id, width) %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::summarise(transcript_length = sum(width)) %>%
    dplyr::rename(ensembl_transcript_id = transcript_id)
  
  #Extract transcript data and merge with other datasets
  transcript_meta = gtf_file[gtf_file$type == "transcript",] %>% as.data.frame() %>% as_tibble()
  transcript_df = dplyr::transmute(transcript_meta, ensembl_gene_id = gene_id, 
                                   ensembl_transcript_id = transcript_id, transcript_biotype, gene_biotype,
                                   chromosome = seqnames, transcript_start = start, transcript_emd = end, strand) %>%
    dplyr::left_join(transcript_tags, by = "ensembl_transcript_id") %>%
    dplyr::left_join(transcript_lengths, by = "ensembl_transcript_id") %>%
    dplyr::mutate(CCDS = ifelse(is.na(CCDS), 0, CCDS)) %>%
    dplyr::mutate(basic = ifelse(is.na(basic), 0, basic)) %>%
    dplyr::mutate(cds_start_NF = ifelse(is.na(cds_start_NF), 0, cds_start_NF)) %>%
    dplyr::mutate(mRNA_start_NF = ifelse(is.na(mRNA_start_NF), 0, mRNA_start_NF)) %>%
    dplyr::mutate(cds_end_NF = ifelse(is.na(cds_end_NF), 0, cds_end_NF)) %>%
    dplyr::mutate(mRNA_end_NF = ifelse(is.na(mRNA_end_NF), 0, mRNA_end_NF)) %>%
    dplyr::rename(is_gencode_basic = basic)
  return(transcript_df)
}

#' For each gene in a gene metadata data frame, mark the transcript with longes sequence.
#'
#' @param gene_metadata Data frame with original gene metadata (required columns:
#' ensembl_gene_id, ensembl_transcript_id, transcript_start, transcript_end, strand)
#'
#' @return Original metadata df with longest_start and longest_end columns added.
markLongestGencodeTranscript <- function(gene_metadata){
  #Make sure that the data frame has all of the required columns
  assertthat::assert_that(assertthat::has_name(gene_metadata, "ensembl_gene_id"))
  assertthat::assert_that(assertthat::has_name(gene_metadata, "ensembl_transcript_id"))
  assertthat::assert_that(assertthat::has_name(gene_metadata, "transcript_length"))
  assertthat::assert_that(assertthat::has_name(gene_metadata, "is_good_reference"))

  #For each gene mark the transcript that has the longest sequence
  longest_tx = dplyr::filter(gene_metadata, is_good_reference == 1) %>% 
    dplyr::group_by(ensembl_gene_id) %>% 
    dplyr::arrange(-transcript_length) %>% #Find smallest possible transcript_start coordinate
    dplyr::select(ensembl_gene_id, ensembl_transcript_id, strand, transcript_start) %>% 
    dplyr::filter(row_number() == 1) %>% 
    dplyr::ungroup() %>%
    dplyr::transmute(ensembl_transcript_id, longest_start = 1, longest_end = 1)
  
  #Combine the start and end coordinates 
  marked_data = dplyr::left_join(gene_metadata, longest_tx, by = "ensembl_transcript_id") %>% 
    dplyr::mutate(longest_start = ifelse(is.na(longest_start), 0, longest_start),
                  longest_end = ifelse(is.na(longest_end), 0, longest_end))
  return(marked_data)
}


#' Prepare transcript metadata before extending truncated transcripts and constructing alternative events
#' 
#' Retain only genes with >= 2 transcripts and mark CDS start and end 'not found' for incomplete transcripts.
#'
#' @param transcript_metadata Data frame of transcript metadata.
#' @param complete_transcripts Biotypes of complete transcripts (CDS start and end NF flags extracted grom GFF)
#' @param incomplete_transcripts Biotypes of trnascripts whose start and end are likely to be missing.
#'
#' @return Filtered transcript metadata df.
#' @export
filterTranscriptMetadata <- function(transcript_metadata, complete_transcripts = c("protein_coding", "lincRNA"),
                                     incomplete_transcripts = c("nonsense_mediated_decay","processed_transcript","retained_intron")){
  
  #Make sure that all necessary columns are present
  assertthat::assert_that(assertthat::has_name(transcript_metadata, "ensembl_gene_id"))
  assertthat::assert_that(assertthat::has_name(transcript_metadata, "ensembl_transcript_id"))
  assertthat::assert_that(assertthat::has_name(transcript_metadata, "transcript_biotype"))
  assertthat::assert_that(assertthat::has_name(transcript_metadata, "transcript_length"))
  assertthat::assert_that(assertthat::has_name(transcript_metadata, "cds_start_NF"))
  assertthat::assert_that(assertthat::has_name(transcript_metadata, "cds_end_NF"))
  assertthat::assert_that(assertthat::has_name(transcript_metadata, "is_gencode_basic"))
  
  filtered_metadata = dplyr::filter(transcript_metadata, transcript_biotype %in% c(complete_transcripts, incomplete_transcripts)) %>% 
    dplyr::group_by(ensembl_gene_id) %>% 
    dplyr::mutate(n_transcripts = length(ensembl_gene_id)) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(n_transcripts > 1) %>%
    dplyr::mutate(cds_start_NF = ifelse(transcript_biotype %in% incomplete_transcripts, 1, cds_start_NF), 
                  cds_end_NF = ifelse(transcript_biotype %in% incomplete_transcripts, 1, cds_end_NF)) %>% 
    dplyr::mutate(cds_start_end_NF = pmax(cds_start_NF, cds_end_NF)) %>%
    dplyr::mutate(is_good_reference = ifelse(cds_start_end_NF == 0, is_gencode_basic, 0)) %>%
    markLongestGencodeTranscript()
  
  return(filtered_metadata)
}
