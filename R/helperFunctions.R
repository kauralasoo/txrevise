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
  return(classification)
}

filterClassification <- function(classification_list, ratio_cutoff = 0.05){
  #Filter out small upstream and dowstram changes, because they cannot be estimated reliably from RNA-Seq
  classification = classification_list$difference
  coding = classification_list$coding
  abs_diff = classification_list$coding
  classification$multiple = 0
  
  ambig_filter = rowSums(sign(classification)) > 1
  upstream_ratio = classification$upstream/apply(classification, 1, max) < ratio_cutoff
  downstream_ratio = classification$downstream/apply(classification,1, max) < ratio_cutoff
  classification$upstream[ambig_filter & upstream_ratio] = 0
  classification$downstream[ambig_filter & downstream_ratio] = 0
  classification$upstream[classification$upstream <= 5] = 0
  classification$downstream[classification$downstream <= 5] = 0
  
  #Mark ambiguous genes
  ambig_filter2 = rowSums(sign(classification)) > 1
  classification$multiple[ambig_filter2] = 1
  
  coding[sign(classification[,1:3]) == 0] = 0
  abs_diff[sign(classification[,1:3]) == 0] = 0
  coding$combined = sign(rowSums(coding))
  abs_diff$combined = sign(rowSums(abs(abs_diff)))
  return(list(difference = classification, coding = coding, abs_diff = abs_diff))
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

