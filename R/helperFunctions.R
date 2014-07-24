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

classifySplicingTable <- function(splicing_table, annotations, cdss = NULL){
  #Decided what kind of splicing events are happening based on a splcing table
  split_txs = strsplit(splicing_table$transcript_ids,",")
  names(split_txs) = splicing_table$gene_id
  classification = applyClassifyDifference(split_txs, annotations, cdss, printProgress = TRUE)
  return(classification)
}

filterClassification <- function(classification_list, ratio_cutoff = 0.05){
  #Filter out small upstream and dowstram changes, because they cannot be estimated reliably from RNA-Seq
  classification = classification_list$transcribed$diff
  abs_diff = classification_list$transcribed$abs_diff
  multiple = rep(0, nrow(classification))
  
  ambig_filter = rowSums(sign(classification)) > 1
  upstream_ratio = classification$upstream/apply(classification, 1, max) < ratio_cutoff
  downstream_ratio = classification$downstream/apply(classification,1, max) < ratio_cutoff
  classification$upstream[ambig_filter & upstream_ratio] = 0
  classification$downstream[ambig_filter & downstream_ratio] = 0
  classification$upstream[classification$upstream <= 5] = 0
  classification$downstream[classification$downstream <= 5] = 0
  
  #Mark ambiguous genes
  ambig_filter2 = rowSums(sign(classification)) > 1
  multiple[ambig_filter2] = 1
  
  abs_diff[sign(classification[,1:3]) == 0] = 0
  return(list(transcribed = list(diff = classification, abs_diff = abs_diff), multiple = multiple))
}

prepareDataForPlotting <- function(filtered_classification_list, remove_multiple = FALSE){
  #Extract data
  diff = filtered_classification_list$difference
  diff$gene_id = rownames(diff)
  code = filtered_classification_list$coding
  code$gene_id = rownames(code)
  
  #Filter out genes with multiple changes
  if(remove_multiple == TRUE){
    diff = diff[diff$multiple == 0,]
    code = code[diff$multiple == 0,]
  }
  
  #Melt diff
  diff_melt = melt(diff[,c("upstream","downstream", "contained", "gene_id")])
  
  #Melt code
  if(remove_multiple == TRUE){
    diff_melt = diff_melt[diff_melt$value > 0,]
    code_melt = melt(code[,c("combined", "gene_id")])
    filter = code_melt$value == 1
    code_melt[filter,]$value = "coding"
    code_melt[!filter,]$value = "non-coding" 
    diff_melt$coding = code_melt[match(diff_melt$gene_id, code_melt$gene_id),]$value
  }
  else{
    code_melt = melt(code[,c("upstream","downstream", "contained", "gene_id")])
    filter = code_melt$value == 1
    code_melt[filter,]$value = "coding"
    code_melt[!filter,]$value = "non-coding" 
    
    #Add coding information to the diff data.frame
    diff_melt$coding = code_melt$value
    diff_melt = diff_melt[diff_melt$value > 0,]
  }
  return(diff_melt)
}

makeClassificationFigure <- function(classification_list){
  type = colSums(sign(a$difference[a$difference$multiple == 0,]))
  colSums(sign(a$difference[a$difference$multiple == 1,]))
}

copyCodingChanges <- function(modified_tx_classification, initial_tx_classification){
  #Copy coding changes from initial classification of transcripts
  coding = initial_tx_classification$coding
  
  #Keep common genes
  gene_ids = rownames(modified_tx_classification$transcribed$diff)
  coding$diff = coding$diff[gene_ids,] * sign(modified_tx_classification$transcribed$diff)
  coding$abs_diff = coding$abs_diff[gene_ids,] * sign(modified_tx_classification$transcribed$diff)
  
  modified_tx_classification[["coding"]] = coding
  return(modified_tx_classification)
}

