constructAlternativeTranscripts <- function(canonical_name, tx_names, exons){
  #Contruct alternative transcripts by only modyfing single part of the annotations
  
  #Add canonical tx to the list
  canonical_tx = unlist(exons[canonical_name])
  tx_list = GRangesList()
  values(canonical_tx) = data.frame(modification = "none", stringsAsFactors = FALSE)
  names(ranges(canonical_tx)) = c()
  canonical_tx = canonical_tx[order(canonical_tx)] #Order exon coordinates by positon
  tx_list[[canonical_name]] = canonical_tx
  
  #Iterate over known alternative txs
  for (tx_name in tx_names){
    tx = unlist(exons[tx_name])
    values(tx) = data.frame(modification = "none", stringsAsFactors = FALSE)
    names(ranges(tx)) = c()
    tx = tx[order(tx)] #Order exon coordinates by positon
    tx_list[[tx_name]] = tx
    
    #Identify all changes between two txs
    changes = indentifyAddedRemovedRegions(tx_name, canonical_name, exons)
    
    #Only proceed if there are any changes
    if(!is.null(changes)){ 
      #Construct alternative txs where only 5', 3' or middle exons have changed
      events = c("upstream", "downstream", "contained")
      for (event in events){
        new_tx = modifyTranscript(canonical_tx, changes[[canonical_name]], changes[[tx_name]], event)
        #Add new tx to the list if it's already not there
        if (!txIsInList(new_tx, tx_list)){
          new_tx_name = paste(canonical_name, ".", length(tx_list), sep = "")
          tx_list[[new_tx_name]] = new_tx
        }
      }
    } 
  }
  return(tx_list)
}

modifyTranscript <- function(canonical_transcript, removedFromTarget, addedInTarget, modification_type){
  #Moify canonical transcript eiter at 5' end, 3'end or inside.
  #print(removedFromTarget)
  #print(addedInTarget)
  
  result = canonical_transcript
  if(length(removedFromTarget) > 0){
    result = setdiff(canonical_transcript, removedFromTarget[values(removedFromTarget)[,modification_type] == 1,])
  }
  if(length(addedInTarget) > 0){
    result = union(result, addedInTarget[values(addedInTarget)[,modification_type] == 1,])   
  }
  
  values(result) = data.frame(modification = rep(modification_type, length(result)), stringsAsFactors = FALSE )
  result = result[order(result)] #Order exon coordinates by positon
  return(result)
}