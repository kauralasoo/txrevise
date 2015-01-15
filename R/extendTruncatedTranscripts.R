extendTranscriptsPerGene <- function(metadata, exons, cdss){
  #Extend truncated transcripts for a single gene
  
  #Identify IDs of transcripts with longest starts and ends
  longest_start_id = dplyr::filter(metadata, longest_start == 1) %>% dplyr::select(ensembl_transcript_id) %>% unlist()
  longest_end_id = dplyr::filter(metadata, longest_end == 1) %>% dplyr::select(ensembl_transcript_id) %>% unlist()
  
  #Identify transcripts with missing ends
  if(length(longest_end_id) == 1){
    missing_ends = dplyr::filter(metadata, cds_end_NF == 1)$ensembl_transcript_id
    missing_tx_list = as.list(missing_ends)
    names(missing_tx_list) = missing_ends
    
    #Extend txs with missing ends to the longest transcripts
    missing_ends_new = lapply(missing_tx_list, extendSingleTranscript, longest_end_id, direction = "downstream", exons)
    print(lapply(missing_ends_new, width))
  }
  
  #Identify transcripts with missing starts
  if(length(longest_start_id) == 1){
    missing_starts = dplyr::filter(metadata, cds_start_NF == 1)$ensembl_transcript_id
    missing_starts_list = as.list(missing_starts)
    names(missing_starts_list) = missing_starts
    
    #Add missing starts
    missing_starts_new = lapply(missing_starts_list, extendSingleTranscript, longest_end_id, direction = "upstream", exons)
    print(lapply(missing_starts_new, width))
  }
}

extendSingleTranscript <- function(tx_id, longest_tx_id, direction, exons){
  #Extend single transcript tx_id based on reference transcript (longest_tx_id) in a specifed direction
  print(tx_id)
  tx_exons = exons[[tx_id]]
  longest_tx_exons = exons[[longest_tx_id]]
  diff = reviseAnnotations::indentifyAddedRemovedRegions(tx_id, longest_tx_id, exons[c(tx_id, longest_tx_id)])
  
  #Identify potential unique exon in the truncated end and set a flag
  extend = TRUE
  tx_spec_exons = diff[[tx_id]]
  if(length(tx_spec_exons) > 0){
    tx_direction_exons = tx_spec_exons[elementMetadata(tx_spec_exons)[,direction] == 1]
    if(length(tx_direction_exons) > 0){
      extend = FALSE
    }
  }
  
  #If truncated transcript has unique exon (or a bit of exon) in the truncated end, then do not extend
  if(!extend){
    return(GRanges())
  }
  #Otherwise proceed with the extension
  else{
    #Find the exons that are missing in truncated transcript
    missing_exons = diff[[longest_tx_id]]
    missing_direction_exons = missing_exons[elementMetadata(missing_exons)[,direction] == 1,]
    
    #Add new exons to the original transcript
    new_exons = mergeGRangesIgnoreMeta(tx_exons, missing_direction_exons)
    return(new_exons)
  }
}