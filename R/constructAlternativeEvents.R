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
  
  #Identify cliques of overlapping transcripts
  tx_cliques = findTranscriptCliques(granges_list)
  clique_list = list()

  #Construct alternative events for each clique separately
  clique_number = 1
  for (tx_clique in tx_cliques){
    clique_exons = granges_list[tx_clique]
    shared_exons = listIntersect(clique_exons)
    
    #Only porceed with the analysis in the clique if there are at least some shared exons
    if(length(shared_exons) > 0){
      #Add shared exons to the exon list
      exon_list = clique_exons
      exon_list[["INTERSECTION"]] = shared_exons
      
      #Identify all changes between transcripts and chared exons
      changes_list = list()
      for(tx_id in names(clique_exons)){
        tx_changes = reviseAnnotations::indentifyAddedRemovedRegions(tx_id, "INTERSECTION", exon_list)[[1]]
        changes_list[[tx_id]] = tx_changes
      }
      
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
      clique_id = paste(gene_id, paste("clique_",clique_number,sep =""), sep = ".")
      clique_number = clique_number + 1
      clique_list[[clique_id]] = event_list
    }
  }
  
  return(clique_list)
}

findTranscriptCliques <- function(granges_list){
  tx_names = names(granges_list)
  
  #Construct a graph
  edges = findOverlaps(granges_list) %>% as.data.frame() %>% 
    dplyr::filter(queryHits != subjectHits) %>% as.matrix() %>% t() %>% as.vector()
  tx_overlap_graph = igraph::graph(edges, directed = FALSE)
  tx_cliques = igraph::max_cliques(tx_overlap_graph)
  clique_list = lapply(tx_cliques, function(x){tx_names[as.vector(x)]}) %>% rev()
  return(clique_list)
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

# Apply constructAlternativeEvents to gene_id and metadata objects
constructAlternativeEventsWrapper <- function(gene_id, gene_metadata, exons, cdss){
  
  #Print current gene_id
  print(gene_id)
  
  #Extract gene data from annotations
  gene_data = reviseAnnotations::extractGeneData(gene_id, gene_metadata, exons, cdss)
  
  #Extend truncated transcripts until the longest transcript
  gene_extended_tx = reviseAnnotations::extendTranscriptsPerGene(gene_data$metadata, gene_data$exons, gene_data$cdss)
  gene_data_ext = reviseAnnotations::replaceExtendedTranscripts(gene_data, gene_extended_tx)
  
  #Construct alternative events
  alt_events = reviseAnnotations::constructAlternativeEvents(gene_data_ext$exons, gene_id)
  return(alt_events)
}

#Combine all alternative events into a single list with appropriate names
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
