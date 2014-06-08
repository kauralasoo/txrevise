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
  return(new_tx_map)
}

classifySplicingTable <- function(splicing_table, annotations, cdss){
  #Decided what kind of splicing events are happening based on a splcing table
  split_txs = strsplit(splicing_table$transcript_ids,",")
  names(split_txs) = splicing_table$gene_id
  classification = applyClassifyDifference(split_txs, annotations, cdss)
  classification = filterClassification(classification)
  return(classification)
}
