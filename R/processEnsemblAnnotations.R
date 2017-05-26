
#' For each gene in a gene metadata data frame, mark the transcripts transcripts
#' that have the longest start and end positions.
#'
#' @param gene_metadata Data frame with original gene metadata (required columns:
#' ensembl_gene_id, ensembl_transcript_id, transcript_start, transcript_end, strand)
#'
#' @return Original metadata df with longest_start and longest_end columns added.
#' @export
markLongestTranscripts <- function(gene_metadata){
  #Make sure that the data frame has all of the required columns
  assertthat::assert_that(assertthat::has_name(gene_metadata, "ensembl_gene_id"))
  assertthat::assert_that(assertthat::has_name(gene_metadata, "ensembl_transcript_id"))
  assertthat::assert_that(assertthat::has_name(gene_metadata, "transcript_start"))
  assertthat::assert_that(assertthat::has_name(gene_metadata, "transcript_end"))
  assertthat::assert_that(assertthat::has_name(gene_metadata, "strand"))
  
  #For each gene mark the transcripts with longest starts and ends
  by_tx_start = dplyr::filter(gene_metadata, is_good_reference == 1) %>% 
    dplyr::group_by(ensembl_gene_id) %>% 
    dplyr::arrange(transcript_start) %>% #Find smallest possible transcript_start coordinate
    dplyr::select(ensembl_gene_id, ensembl_transcript_id, strand, transcript_start) %>% 
    dplyr::filter(row_number() == 1) %>% 
    dplyr::ungroup() %>%
    #Use strand infromation to decide whether it is at the start or the end of the transcript
    dplyr::transmute(ensembl_transcript_id, longest_start = sign(strand +1), longest_end = sign(abs(strand-1)))
  
  by_tx_end = dplyr::filter(gene_metadata, is_good_reference == 1) %>% 
    dplyr::group_by(ensembl_gene_id) %>% 
    dplyr::arrange(desc(transcript_end)) %>% #Find largest possible transcript_end coordinate
    dplyr::select(ensembl_gene_id, ensembl_transcript_id, strand, transcript_end) %>% 
    dplyr::filter(row_number() == 1) %>% 
    dplyr::ungroup() %>%
    #Use strand infromation to decide whether it is at the start or the end of the transcript
    dplyr::transmute(ensembl_transcript_id, longest_start = sign(abs(strand-1)), longest_end = sign(strand +1))
  
  #Combine the start and end coordinates
  both_ends = dplyr::left_join(gene_metadata, by_tx_start, by = "ensembl_transcript_id") %>% 
    dplyr::left_join(by_tx_end, by = "ensembl_transcript_id") %>%
    dplyr::transmute(ensembl_transcript_id, 
                     longest_start = ifelse(is.na(longest_start.x), 0, longest_start.x) + ifelse(is.na(longest_start.y), 0, longest_start.y),
                     longest_end = ifelse(is.na(longest_end.x), 0, longest_end.x) + ifelse(is.na(longest_end.y), 0, longest_end.y))
  
  marked_data = dplyr::left_join(gene_metadata, both_ends, by = "ensembl_transcript_id")
  return(marked_data)
}


#' For each gene in a gene metadata data frame, mark the transcript with longes sequence.
#'
#' @param gene_metadata Data frame with original gene metadata (required columns:
#' ensembl_gene_id, ensembl_transcript_id, transcript_start, transcript_end, strand)
#'
#' @return Original metadata df with longest_start and longest_end columns added.
#' @export
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
    reviseAnnotations::markLongestGencodeTranscript()
  
  return(filtered_metadata)
}
