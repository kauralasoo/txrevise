
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
  by_tx_start = dplyr::filter(filtered_data, is_good_reference == 1) %>% 
    dplyr::group_by(ensembl_gene_id) %>% 
    dplyr::arrange(transcript_start) %>% #Find smallest possible transcript_start coordinate
    dplyr::select(ensembl_gene_id, ensembl_transcript_id, strand, transcript_start) %>% 
    dplyr::filter(row_number() == 1) %>% 
    dplyr::ungroup() %>%
    #Use strand infromation to decide whether it is at the start or the end of the transcript
    dplyr::transmute(ensembl_transcript_id, longest_start = sign(strand +1), longest_end = sign(abs(strand-1)))
  
  by_tx_end = dplyr::filter(filtered_data, is_good_reference == 1) %>% 
    dplyr::group_by(ensembl_gene_id) %>% 
    dplyr::arrange(desc(transcript_end)) %>% #Find largest possible transcript_end coordinate
    dplyr::select(ensembl_gene_id, ensembl_transcript_id, strand, transcript_end) %>% 
    dplyr::filter(row_number() == 1) %>% 
    dplyr::ungroup() %>%
    #Use strand infromation to decide whether it is at the start or the end of the transcript
    dplyr::transmute(ensembl_transcript_id, longest_start = sign(abs(strand-1)), longest_end = sign(strand +1))
  
  #Combine the start and end coordinates
  both_ends = dplyr::left_join(filtered_data, by_tx_start, by = "ensembl_transcript_id") %>% 
    dplyr::left_join(by_tx_end, by = "ensembl_transcript_id") %>%
    dplyr::transmute(ensembl_transcript_id, 
                     longest_start = ifelse(is.na(longest_start.x), 0, longest_start.x) + ifelse(is.na(longest_start.y), 0, longest_start.y),
                     longest_end = ifelse(is.na(longest_end.x), 0, longest_end.x) + ifelse(is.na(longest_end.y), 0, longest_end.y))
  
  marked_data = dplyr::left_join(filtered_data, both_ends, by = "ensembl_transcript_id")
  return(marked_data)
}
