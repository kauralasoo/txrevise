#' Construct alternative transcription events form a list of transcripts.
#'
#' Also works if not all transcripts of a gene overlap with each other. In this case the transcripts are first 
#' split into overlapping sets (cliques in the transcript overlap graph) and events are calculated within
#' each set independently.
#' 
#' @param granges_list List of GRanges objects containing the exon coordinates of each transcript 
#' (One object per transcript).
#' @param gene_id ID of the the gene to which the transcripts belong to.
#' @param max_internal_diff maximal internal difference between two evets for them to be considered the same.
#' @param max_start_end_diff maximal difference at the start and end of two events for them to be considered 
#' the same.
#' @return Three-level list. First level correspond the each set of overlapping transcripts, second 
#' level contains lists of downstream, upstream and contained alternative transcription events.
#' @author Kaur Alasoo
#' @export 
constructAlternativeEvents <- function(granges_list, gene_id, max_internal_diff = 10, max_start_end_diff = 25){
  
  #Identify groups of overlapping transcripts
  tx_groups = identifyTranscriptGroups(granges_list)
  group_list = list()

  #Construct alternative events for each clique separately
  group_number = 1
  for (tx_group in tx_groups){
    shared_exons = listIntersect(tx_group)
    
    #Only porceed with the analysis in the clique if there are at least some shared exons
    if(length(shared_exons) > 0){
      #Identify all changes between transcripts and shared exons
      changes_list = differenceFromShared(tx_group)
      
      #Identify all types of alternative events
      event_types = c("upstream", "downstream","contained")
      event_list = list()
      for (event_type in event_types){
        event_changes = lapply(changes_list, extractExonsByType, event_type) %>% removeMetadata()
        event_transcripts = lapply(event_changes, function(x){event = c(x, shared_exons) %>% sort() %>% reduce()})
        #Remove duplicates
        event_transcripts = event_transcripts[!duplicated(event_transcripts)]
        #Remove events that are too similar to other events
        event_transcripts = mergeByMaxDifference(event_transcripts, max_internal_diff, max_start_end_diff)
        #Add events to the event list
        event_list[[event_type]] = event_transcripts
      }
      #Add events to clique list
      group_id = paste(gene_id, paste("grp_",group_number,sep =""), sep = ".")
      group_number = group_number + 1
      group_list[[group_id]] = event_list
    }
  }
  
  return(group_list)
}

extractExonsByType <- function(granges, type){
  #From a indentifyAddedRemovedRegions object extract changes by the part of the transcript that is affectes
  
  #Check that type is specified correctly
  correct_event_types = c("upstream","downstream", "contained")
  if (!(type %in% correct_event_types)){
    stop("Type has to be either 'upstream', 'downstream' or 'contained'")
  }
  
  #Deal with granges objects that have zero rows
  if(length(granges) == 0){
    return(granges)
  } else{
    filtered_exons = granges[elementMetadata(granges)[,type] == 1,]
    return(filtered_exons)
  }
}

#' If two transcripts in a GRanges list a are not sufficiently different from each other, then keep only one.
mergeByMaxDifference <- function(granges_list, max_internal_diff = 10, max_start_end_diff = 25){
  
  #Mergeing only makes sense for lists that have more than one element
  if (length(granges_list) < 2){
    return(granges_list)
  }

  #Extract transcript names
  tx_names = names(granges_list)
  new_list = granges_list[1]
  
  #Iterate over all transcripts
  for(i in 2:length(tx_names)){
    current_name = tx_names[i]
    current_tx = granges_list[[current_name]]
    
    #Calculate the difference between the current transcript and previosuly added transcripts
    diff_mat = dplyr::data_frame()
    for(old_tx in new_list){
      diff = basesDifferent(current_tx, old_tx)
      if(nrow(diff_mat) == 0){
        diff_mat = diff
      }else{
        diff_mat = rbind(diff_mat, diff)
      }
    }    
    #Calculate internal diff
    diff_mat = dplyr::mutate(diff_mat, internal_diff = total_diff - start_end_diff)
    #If the current transcript is sufficiently different then added it to the list of transcripts
    if ((min(diff_mat$internal_diff) > max_internal_diff) | (min(diff_mat$start_end_diff) > max_start_end_diff)){
      new_list[[current_name]] = current_tx
    }
  }
  return(new_list)
}

#' Apply constructAlternativeEvents() to an gene_id and metadata objects. Extend truncated transcripts first.
#'
#' @param gene_id Focal gene id.
#' @param gene_metadata Tibble containing gene metadata from the prepareAnnotations script.
#' @param exons List of annotated exons
#' @param cdss List of annotated CDSs
#' @param max_internal_diff maximal internal difference between two evets for them to be considered the same.
#' @param max_start_end_diff maximal difference at the start and end of two events for them to be considered 
#' the same.
#'
#' @return Three-level list. First level correspond the each set of overlapping transcripts, second 
#' level contains lists of downstream, upstream and contained alternative transcription events.
#' @export
constructAlternativeEventsWrapper <- function(gene_id, gene_metadata, exons, cdss, max_internal_diff = 10, max_start_end_diff = 25){
  
  #Print current gene_id
  print(gene_id)
  
  #Extract gene data from annotations
  gene_data = extractGeneData(gene_id, gene_metadata, exons, cdss)
  
  #Extend truncated transcripts until the longest transcript
  gene_extended_tx = extendTranscriptsPerGene(gene_data$metadata, gene_data$exons, gene_data$cdss)
  gene_data_ext = replaceExtendedTranscripts(gene_data, gene_extended_tx)
  
  #Construct alternative events
  alt_events = constructAlternativeEvents(gene_data_ext$exons, gene_id, max_internal_diff, max_start_end_diff)
  return(alt_events)
}

#' Combine all alternative events into a single list with appropriate names
#' @export
flattenAlternativeEvents <- function(alt_events){
  flat_event_list = list()
  #Iterate through all events
  gene_cliques = names(alt_events)
  for (gene_clique in gene_cliques){
    clique_events = alt_events[[gene_clique]]
    event_types = names(clique_events)
    for (event in event_types){
      event_list = clique_events[[event]]
      #Only add events to the final list if there are at least 2 alternative transcripts
      if (length(event_list) > 1){
        new_names = paste(gene_clique, event, names(event_list), sep =".")
        names(event_list) = new_names
        flat_event_list = c(flat_event_list, event_list)
      }
    }
  }
  return(flat_event_list)
}

#' Construct event metadata
#' @export
constructEventMetadata <- function(transcript_ids){
  event_metadata = data.frame(transcript_id = transcript_ids) %>% 
    tidyr::separate(transcript_id, c('ensembl_gene_id', 'grp_id', 'event_type','ensembl_transcript_id'), 
                    sep = "\\.", remove = F) %>%
    dplyr::mutate(gene_id = paste(ensembl_gene_id, grp_id, event_type, sep = ".")) %>%
    dplyr::mutate(transcript_id = as.character(transcript_id)) %>%
    tbl_df()
  return(event_metadata)
}

makeBinary <- function(ids, len){
  result_vector = rep(0, len)
  result_vector[ids] = 1
  return(result_vector)
}

identifyTranscriptGroups <- function(granges_list){
  
  #Construct disjoint exons and count overlaps
  disjoint_exons = purrr::reduce(granges_list, c) %>% 
    GenomicRanges::disjoin() %>% 
    GenomicRanges::sort()
  overlaps = purrr::map(S4Vectors::as.list(granges_list), ~GenomicRanges::findOverlaps(.,disjoint_exons) %>% 
                          S4Vectors::subjectHits() %>%
                          makeBinary(.,length(disjoint_exons))) %>%
    bind_rows()
  
  #Construct unique configurations
  unique_conf = unique(overlaps) %>% as.matrix()
  rownames(unique_conf) = apply(unique_conf, 1, paste, collapse = "")
  
  #Count exons and transcripts in each configuration
  config_df = dplyr::data_frame(sharing_config = apply(overlaps, 1, paste, collapse = ""), exon_in_tx_count = rowSums(overlaps))
  config_counts = dplyr::group_by(config_df, sharing_config) %>%
    dplyr::summarise(n_exons_shared = length(sharing_config), n_transcripts_sharing = exon_in_tx_count[[1]]) %>%
    dplyr::arrange(-n_exons_shared, -n_transcripts_sharing) %>%
    dplyr::filter(n_transcripts_sharing > 1)
  
  #Extract transcripts
  if(nrow(config_counts) == 1){
    grp_1_ids = names(which(unique_conf[config_counts$sharing_config[[1]],] == 1))
    return(list(grp_1 = granges_list[grp_1_ids], grp_2 = NULL))
  } else if(nrow(config_counts) >= 2){
    grp_1_ids = names(which(unique_conf[config_counts$sharing_config[[1]],] == 1))
    grp_2_ids = names(which(unique_conf[config_counts$sharing_config[[2]],] == 1))
    return(list(grp_1 = granges_list[grp_1_ids], grp_2 = granges_list[grp_2_ids]))
  } else {
    return(list(grp_1 = NULL, grp_2 = NULL))
  }
}


#Identify all changes between transcripts and shared exons
differenceFromShared <- function(transcript_list){
  #Identify shared exons
  shared_exons = listIntersect(transcript_list)
  
  #Only proceed if there are at least some shared exons
  if(length(shared_exons) > 0){
    
    #Add shared exons to the exon list
    exon_list = transcript_list
    exon_list[["INTERSECTION"]] = shared_exons
    
    #Identify all changes between transcripts and shared exons
    changes_list = list()
    for(tx_id in names(transcript_list)){
      tx_changes = indentifyAddedRemovedRegions(tx_id, "INTERSECTION", exon_list)[[1]]
      changes_list[[tx_id]] = tx_changes
    }
    return(changes_list)
  } else{
    warning("No exon is shared between all transcripts.")
    return(NULL)
  }
}
