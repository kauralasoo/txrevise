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

removeMetadata <- function(granges_list){
  list = lapply(granges_list, function(x){elementMetadata(x) <- c(); return(x)})
  return(GRangesList(list))
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


#' Extract gene-specifc data from annotations data frame and exons/cdss lists.
#' 
#' @param gene_id Ensembl id of the gene.
#' @param annotations_df data frame containing transcript metadata.
#' @param exons named list of GRanges objects containing the exon coordinates for each transcript.
#' @param cdss named list of GRanges objects containing the CDS coordinates for each transcript.
#' @return List of 3 objects: filtered gene metadata, filtered exon coordinates (gene trancripts only),
#' filtered CDS coordinates (gene transcripts only)
#' @author Kaur Alasoo
#' @export 
extractGeneData <- function(gene_id, annotations_df, exons, cdss){
  gene_data = dplyr::filter(annotations_df, ensembl_gene_id == gene_id)
  tx_ids = dplyr::select(gene_data, ensembl_transcript_id) %>% unlist()
  gene_exons = exons[tx_ids]
  gene_cdss = cdss[intersect(tx_ids, names(cdss))]
  gene_data = list(metadata = gene_data, exons = gene_exons, cdss = gene_cdss)
}

#' Replace extended trancripts in gene-specific data list.
#' 
#' @param gene_data List produced by the extractGeneData function.
#' @param extended_transcripts List produced by the extendTranscriptsPerGene function.
#' @return Update gene_data list where original transcripts have been replaced by extended transcripts.
#' @author Kaur Alasoo
#' @export 
replaceExtendedTranscripts <- function(gene_data, extended_transcripts){
  gene_data$exons = removeMetadata(gene_data$exons)
  gene_data$cdss = removeMetadata(gene_data$cdss)
  gene_data$exons[names(extended_transcripts$exons)] = extended_transcripts$exons
  gene_data$cdss[names(extended_transcripts$cdss)] = extended_transcripts$cdss
  return(gene_data)
}