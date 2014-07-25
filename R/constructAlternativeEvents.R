constructAlternativeEvents <- function(canonical_name, tx_names, exons){
  #Construct alternative transcription events from 2 or more transcripts
  
  event_list = GRangesList()
  
  #Iterate over all alternative transcripts
  upstream_list = GRangesList()
  downstream_list = GRangesList()
  contained_list = GRangesList()
  
  for (tx_name in tx_names){
    #Iterate over alternative transcripts
    changes = indentifyAddedRemovedRegions(tx_name, canonical_name, exons)
    
    if(!is.null(changes)){ 
      #Extract changes
      canonical_changes = changes[[canonical_name]]
      tx_changes = changes[[tx_name]]
      shared_exons = changes$shared_exons
      
      upstream_events = constructEventsByType(canonical_changes, tx_changes, shared_exons, type = "upstream", id_base = tx_name)
      downstream_events = constructEventsByType(canonical_changes, tx_changes, shared_exons, type = "downstream", id_base = tx_name)
      contained_events = constructEventsByType(canonical_changes, tx_changes, shared_exons, type = "contained", id_base = tx_name)
      
      #Merge into respective lists
      upstream_list = mergeGRangesList(upstream_list, upstream_events)
      downstream_list = mergeGRangesList(downstream_list, downstream_events)
      contained_list = mergeGRangesList(contained_list, contained_events)
    }
  }
  #Merge all events together
  all_events = c(upstream_list, downstream_list, contained_list)
  return(all_events)
}

constructEventsByType <- function(canonical_changes, tx_changes, shared_exons, type, id_base){
  #Construct two new alternative transcription events based on changes and type.
  
  #Check that type is specified correctly
  if (!(type %in% c("upstream","downstream", "contained"))){
    stop("Type has to be either 'upstream', 'downstream' or 'contained'")
  }
  
  #Keep changes of specific type
  canonical_changes = canonical_changes[as.numeric(values(canonical_changes)[,type]) == 1]
  tx_changes = tx_changes[as.numeric(values(tx_changes)[,type]) == 1] 
  total_length = length(tx_changes) + length(canonical_changes)
  if(total_length > 0){
    #Remove annotation columns before concatenation
    values(canonical_changes) = c()
    values(tx_changes) = c()
    
    #Make two new transcripts
    tx1 = c(shared_exons, canonical_changes)
    tx1 = tx1[order(tx1)]
    values(tx1) = data.frame(modification = type, stringsAsFactors = FALSE)
    tx2 = c(shared_exons, tx_changes)
    tx2 = tx2[order(tx2)]
    values(tx2) = data.frame(modification = type, stringsAsFactors = FALSE)
    
    #Return as list
    result = list()
    result[[paste(id_base, ".", type, "1", sep = "")]] = tx1
    result[[paste(id_base, ".", type, "2", sep = "")]] = tx2
    return(result)
  } else {
    return(NULL)
  }
}

mergeGRangesList <- function(list1, list2){
  #Merges entries form list2 into list1 and ignore duplicates
  names2 = names(list2)
  for(name in names2){
    if(length(list1) == 0){
      list1[[name]] = list2[[name]]
    } else if (!txIsInList(list2[[name]], list1)){
      list1[[name]] = list2[[name]]
    }
  }
  return(list1)
}
