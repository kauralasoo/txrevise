makeTranscriptList <- function(transcript_ids){
  #Convert vector with comma-separeted transcript ids into a list suitable for makeNewTranscripts
  split_txs = strsplit(transcript_ids,",")
  primary_tx = unlist(lapply(split_txs, function(x){x[1]}))
  alternative_tx = lapply(split_txs, function(x){x[2]})
  names(alternative_tx) = primary_tx
  return(alternative_tx)
}

extractGeneTxMap <- function(annotations){
  #extract GeneTxMap from newly created annotations file
  new_tx_map = values(annotations[annotations$type == "mRNA",])
  new_tx_map = new_tx_map[,c("Parent","ID")]
  colnames(new_tx_map) = c("gene_id", "transcript_id")
  new_tx_map = as.data.frame(new_tx_map)
  return(new_tx_map)
}

extractFeature <- function(vector){
  #Return the one and only unqiue value in a vector, factor, RLE, etc. Complain if there is more than one unique value.
  value = unique(as.vector(vector))
  if(length(value) == 1){
    return(value)
  } else{
    warning("More than one unqiue value in the vector.")
    return(NULL)
  }
}

txIsInList <- function(tx, tx_list){
  #Test if a transcript is already in a GRangesList object
  result = FALSE
  itr = 1
  if(length(tx_list) == 0){
    warning("tx_list object is empty")
    return(FALSE)
  }
  for (itr in 1:length(tx_list)){
    current_tx = tx_list[[itr]]
    if (length(tx) == length(current_tx)){
      if(all(tx == current_tx)){
        result = TRUE
      }
    }
  }
  return(result)
}

mergeGRangesIgnoreMeta <- function(gr1, gr2){
  #Merge two GRanges object and ignore all metadata information
  
  #Remove metadata from the GRanges
  elementMetadata(gr1) <- c()
  elementMetadata(gr2) <- c()
  result = reduce(sort(c(gr1, gr2)))
  return(result)
}

classifySplicingTable <- function(splicing_table, annotations, cdss = NULL){
  #Decided what kind of splicing events are happening based on a splcing table
  split_txs = strsplit(splicing_table$transcript_ids,",")
  names(split_txs) = splicing_table$gene_id
  classification = applyClassifyDifference(split_txs, annotations, cdss, printProgress = TRUE)
  return(classification)
}

constructEventMask <- function(upstream_res, contained_res, downstream_res, diff_matrix){
  #Combine all upstream and downstream into a single mask to filter inital classification results
  #TODO: This code does not work anymore after tidy implementation of classification code
  
  #All significant hits
  significant_hits = unique(c(downstream_res$collapsed$gene_id, upstream_res$collapsed$gene_id, contained_res$collapsed$gene_id))
  
  #Make an empty mask
  mask = diff_matrix
  mask = mask[significant_hits,]
  mask[,] = 0
  
  #Populate with ones where there is a significant event
  mask[contained_res$collapsed$gene_id,"contained"] = 1
  mask[upstream_res$collapsed$gene_id,"upstream"] = 1
  mask[downstream_res$collapsed$gene_id,"downstream"] = 1
  
  return(mask)
}

applyEventMask <- function(classification_table, mask){
  #Apply a mask from constructEventMask to result from applyClassifyDifference
  #TODO: This code does not work anymore after tidy implementation of classification code
  classification_table$transcribed$diff = classification_table$transcribed$diff[rownames(mask),]*sign(mask)
  classification_table$transcribed$abs_diff = classification_table$transcribed$abs_diff[rownames(mask),]*sign(mask)
  classification_table$coding$diff = classification_table$coding$diff[rownames(mask),]*sign(mask)
  classification_table$coding$abs_diff = classification_table$coding$abs_diff[rownames(mask),]*sign(mask)
  
  return(classification_table)
}