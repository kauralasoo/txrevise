extendTranscriptsPerGene <- function(metadata, exons, cdss){
  #Extend truncated transcripts for a single gene
  
  #Identify IDs of transcripts with longest starts and ends
  longest_start_id = dplyr::filter(metadata, longest_start == 1) %>% dplyr::select(ensembl_transcript_id) %>% unlist()
  longest_end_id = dplyr::filter(metadata, longest_end == 1) %>% dplyr::select(ensembl_transcript_id) %>% unlist()
  
  #Go through all mRNA ends first
  truncated_transcripts = dplyr::filter(metadata, cds_start_end_NF > 0)
  new_exons = extendTranscripts(truncated_transcripts, longest_start_id, longest_end_id, exons[metadata$ensembl_transcript_id])
  
  #Modify CDS for the transcripts that have been extended
  extended_transcripts = names(new_exons)
  truncated_cds = dplyr::filter(metadata, ensembl_transcript_id %in% extended_transcripts)
  new_cdss = extendTranscripts(truncated_cds, longest_start_id, longest_end_id, cdss[metadata$ensembl_transcript_id])
  
  return(list(exons = new_exons, cdss = new_cdss))
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
      #Identify original exons of the turncated transcript
      tx_dir_full_exons = tx_exons[subjectHits(findOverlaps(tx_direction_exons, tx_exons))]
      #Find if the orginal exons overlap exons int the longest transcript
      query_hits = queryHits(findOverlaps(tx_dir_full_exons, longest_tx_exons))
      if(length(query_hits) > length(tx_direction_exons)){ #Some query exons do not overlap exons in longest transcript
        extend = FALSE
      } else{
        if (sum(width(tx_direction_exons)) > 15){ #All query exons overlap the longest trancript, but the differences is greater than 10 nucleotides
          extend = FALSE
        }
      }
    }
  }
  
  #If truncated transcript has unique exon (or a bit of exon) in the truncated end, then do not extend
  print(extend)
  if(!extend){
    return(GRanges())
  }
  #Otherwise proceed with the extension
  else{
    #Find the exons that are missing in truncated transcript
    missing_exons = diff[[longest_tx_id]]
    missing_direction_exons = missing_exons[elementMetadata(missing_exons)[,direction] == 1,]
    #print(missing_direction_exons)
    
    #Add new exons to the original transcript
    new_exons = mergeGRangesIgnoreMeta(tx_exons, missing_direction_exons)
    return(new_exons)
  }
}

extendTranscripts <- function(truncated_transcripts_meta, longest_start_id, longest_end_id, feature_list){
  #Extend all transcripts specifed in the truncated_transcripts_meta file
  #truncated_transcripts_meta - metadata for the truncated tranacripts of the genes
  #longest_end_id - id of the transcript with the longest end flag
  #longest_start_id - id of the transcript with the longest start flag
  #feature_list - GRangesList object containg either mRNA or CDS regions
  
  results = list()
  #Identify transcripts with missing ends
  if(length(longest_end_id) == 1){
    missing_ends = dplyr::filter(truncated_transcripts_meta, cds_end_NF == 1)$ensembl_transcript_id
    missing_tx_list = as.list(missing_ends)
    names(missing_tx_list) = missing_ends
    
    #Extend txs with missing ends to the longest transcripts
    missing_ends_new = lapply(missing_tx_list, extendSingleTranscript, longest_end_id, direction = "downstream", feature_list)
    new_ends_exon_counts = lapply(missing_ends_new, length) %>% unlist()
    missing_ends_new = missing_ends_new[new_ends_exon_counts > 0]
    results = missing_ends_new
    modified_features = c(GRangesList(results), removeMetadata(feature_list[setdiff(names(feature_list), names(results))]))
  }
  
  #Identify transcripts with missing starts
  if(length(longest_start_id) == 1){
    missing_starts = dplyr::filter(truncated_transcripts_meta, cds_start_NF == 1)$ensembl_transcript_id
    missing_starts_list = as.list(missing_starts)
    names(missing_starts_list) = missing_starts
    
    #Add missing starts
    missing_starts_new = lapply(missing_starts_list, extendSingleTranscript, longest_end_id, direction = "upstream", modified_features)
    new_starts_exon_counts = lapply(missing_starts_new, length) %>% unlist()
    missing_starts_new = missing_starts_new[new_starts_exon_counts > 0]
    results = c(missing_starts_new, results[setdiff(names(results), names(missing_starts_new))])
    
  }
  return(GRangesList(results))
}

