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

classifySplicingTable <- function(splicing_table, annotations, cdss = NULL){
  #Decided what kind of splicing events are happening based on a splcing table
  split_txs = strsplit(splicing_table$transcript_ids,",")
  names(split_txs) = splicing_table$gene_id
  classification = applyClassifyDifference(split_txs, annotations, cdss, printProgress = TRUE)
  return(classification)
}

prepareDataForPlotting <- function(filtered_classification_list, remove_multiple = FALSE){
  #Extract data
  diff = filtered_classification_list$transcribed$diff
  diff$gene_id = rownames(diff)
  code = filtered_classification_list$coding$diff
  code$gene_id = rownames(code)
  
  #Filter out genes with multiple changes
  if(remove_multiple == TRUE){
    diff$multiple = sign(rowSums(sign(filtered_classification_list$transcribed$diff)) -1)
    code$combined = rowSums(sign(filtered_classification_list$coding$diff))
    diff = diff[diff$multiple == 0,]
    code = code[diff$multiple == 0,]
  }
  
  #Melt diff
  diff_melt = melt(diff[,c("upstream","downstream", "contained", "gene_id")])
  
  #Melt code
  if(remove_multiple == TRUE){
    diff_melt = diff_melt[diff_melt$value > 0,]
    code_melt = melt(code[,c("combined", "gene_id")])
    filter = code_melt$value > 0
    code_melt[filter,]$value = "coding"
    code_melt[!filter,]$value = "non-coding" 
    diff_melt$coding = code_melt[match(diff_melt$gene_id, code_melt$gene_id),]$value
  }
  else{
    code_melt = melt(code[,c("upstream","downstream", "contained", "gene_id")])
    filter = code_melt$value > 0
    code_melt[filter,]$value = "coding"
    code_melt[!filter,]$value = "non-coding" 
    
    #Add coding information to the diff data.frame
    diff_melt$coding = code_melt$value
    diff_melt = diff_melt[diff_melt$value > 0,]
  }
  return(diff_melt)
}

constructEventMask <- function(upstream_res, contained_res, downstream_res, diff_matrix){
  #Combine all upstream and downstream into a single mask to filter inital classification results
  
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
  classification_table$transcribed$diff = classification_table$transcribed$diff[rownames(mask),]*sign(mask)
  classification_table$transcribed$abs_diff = classification_table$transcribed$abs_diff[rownames(mask),]*sign(mask)
  classification_table$coding$diff = classification_table$coding$diff[rownames(mask),]*sign(mask)
  classification_table$coding$abs_diff = classification_table$coding$abs_diff[rownames(mask),]*sign(mask)
  
  return(classification_table)
}

plotCompareProportions <- function(classification_list, facet_by){
  
  #Prepare data for plotting
  data = prepareDataListForPlotting(classification_list)
  
  #Facet by type of splicing event
  if(facet_by == "type"){
    plot = ggplot(data, aes(x = comparison, fill = coding)) + 
      geom_bar(aes(y = (..count..)/tapply(..count..,..x..,sum)[..x..])) + 
      facet_grid(~variable)
  } else if (facet_by == "comparison"){ #Facet by comparison
    plot = ggplot(data, aes(x = variable, fill = coding)) + 
      geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + 
      facet_grid(~comparison) 
  }
  
  #Common code:
  plot = plot + scale_y_continuous(labels = percent_format()) +
    theme(legend.position = c(0,1), legend.justification = c(0, 1), legend.key.size = unit(0.4, "line")) + 
    ylab("Percentage") + 
    scale_fill_discrete(name = "Type") + 
    theme(axis.text.x = element_text(angle = 20), axis.title.x = element_blank())
  return(plot)
}

prepareDataListForPlotting <- function(classification_list, type_labels = c("5' End", "3' End", "Internal")){
  #Apply prepareDataForPlotting to a list of classification objects
  plot_data = lapply(classification_list, prepareDataForPlotting, remove_multiple = FALSE)
  names = names(plot_data)
  complete_data = c()
  for (name in names){
    data_set = plot_data[[name]]
    data_set$comparison = name
    complete_data = rbind(complete_data, data_set)
  }
  complete_data$comparison = factor(complete_data$comparison, levels = names)
  levels(complete_data$variable) = type_labels
  return(complete_data)
}

extractEventLengths <- function(class_object, max_length = 7500){
  #Convert classification object into a data.frame of event lengths
  transcribed = melt(class_object$transcribed$diff)
  coding = melt(class_object$coding$diff)
  coding$type = "coding"
  transcribed$type = "all"
  
  #Combine the two
  data = rbind(transcribed, coding)
  data = data[data$value != 0,]
  levels(data$variable) = c("5'","3'","Internal")
  data$joint = as.factor(paste(data$variable, data$type))
  
  #Filter very long events
  data$value[data$value > max_length] = max_length
  
  return(data)
}

plotEventLengths <- function(event_lengths, jitter_width = .2, jitter_size = 0.8){
  #Plot event lengths for different types of events
  plot = ggplot(event_lengths, aes(x = joint, y = value, fill = type, color = type)) + 
    geom_violin(adjust = 0.8, scale = "width") + 
    scale_y_continuous(limit = c(0,7500)) +
    theme(axis.title.x = element_blank(), legend.position = c(1,1), 
          legend.justification = c(1, 1), 
          legend.key.size = unit(1, "line")) +
    ylab("Length (nt)") + 
    geom_jitter(position = position_jitter(width=jitter_width), color = "black", size = jitter_size)
  return(plot)
}

extractSplicingEvents <- function(splicing_class, unique = TRUE){
  #Extract gene names with different alternative transcription events
  
  #Filter unique events
  single_events = rep(TRUE, nrow(splicing_class$transcribed$diff))
  if(unique){
    single_events = rowSums(sign(splicing_class$transcribed$diff)) == 1
  }
  transcribed_mat = splicing_class$transcribed$diff[single_events,]
  coding_mat = splicing_class$coding$diff[single_events,]
  
  #Extract coding and non-coding events for each subtype
  result_list = list()
  event_types = colnames(transcribed_mat)
  for (event in event_types){
    print(event)
    coding_label = paste(event, "coding", sep = "_")
    noncoding_label = paste(event, "non-coding", sep ="_")
    
    coding_genes = transcribed_mat[transcribed_mat[,event] > 0 & coding_mat[,event] > 0,]
    noncoding_genes = transcribed_mat[transcribed_mat[,event] > 0 & coding_mat[,event] == 0,]
    
    result_list[[coding_label]] = rownames(coding_genes)
    result_list[[noncoding_label]] = rownames(noncoding_genes)
  }
  return(result_list)
}

mergeGRangesIgnoreMeta <- function(gr1, gr2){
  #Merge two GRanges object and ignore all metadata information
  
  #Remove metadata from the GRanges
  elementMetadata(gr1) <- c()
  elementMetadata(gr2) <- c()
  result = reduce(sort(c(gr1, gr2)))
  return(result)
}
