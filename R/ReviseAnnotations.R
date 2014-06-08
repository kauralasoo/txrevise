
compareRanges <- function(target_region, shared_region){
  #Compare two Granges to each other and decide if the target is left, right or contained in the target region
  tx_strand = extractFeature(strand(shared_region))
  results = c(0,0,0)
  names(results) = c("upstream", "downstream", "contained")
  
  target_start = start(target_region)
  target_end = end(target_region)
  shared_start = start(shared_region)
  shared_end = end(shared_region)
    
  if ((target_start > shared_start & target_start < shared_end) | (target_end > shared_start & target_end < shared_end)){
    results["contained"] = 1
  }
  else if (target_end < shared_start){
    if(tx_strand == "+") {results["upstream"] = 1}
    if(tx_strand == "-") {results["downstream"] = 1}
  }
  else if (target_start > shared_end){
    if(tx_strand == "+") {results["downstream"] = 1}
    if(tx_strand == "-") {results["upstream"] = 1}
  }
  return(results)
} 

compareGRanges <- function(target_regions, shared_region){
  #Compare multiple regions in a GRanges target_regions object to a single shared_region object.
  #Decide the type of change happening in each region.
  
  #If no target regions the return empty GRanges
  if(length(target_regions) == 0){
    return(target_regions)
  } 
  
  if(length(shared_region) == 0){
    warning("Shared regio is empty.")
    return(shared_region)
  }
  
  #Otherwise construct data.frame for metadata
  results = c()
  n = length(target_regions)
  for (i in 1:n){
    comparison = compareRanges(target_regions[i], shared_region)
    results = rbind(results, comparison)
  }
  rownames(results) = c()
  values(target_regions) = results
  return(target_regions)
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
    changes = indentifyAddedRemovedRegions(tx, canonical_tx)
    #print(changes)
    
    #Only proceed if there are any changes
    if(!is.null(changes)){ 
      #Construct alternative txs where only 5', 3' or middle exons have changed
      events = c("upstream", "downstream", "contained")
      for (event in events){
        new_tx = modifyTranscript(canonical_tx, changes$removedFromFirst, changes$addedInSecond, event)
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

addTranscriptAnnotations <- function(transcript_list, gene_id){
  transcript_names = names(transcript_list)
  canonical_tx_name = transcript_names[1]
  
  annotations = GRanges()
  for (tx_name in transcript_names){
    tx = transcript_list[[tx_name]]
    #tx_list[[tx_name]] = tx
    
    #Add metadata columns for each transcript
    metadata = values(tx)
    metadata$Parent = tx_name
    metadata$type = "exon"
    metadata$source = "reviseAnnotations"
    metadata$ID = paste("exon:", tx_name, ":", c(1:length(tx)), sep = "")
    values(tx) = metadata
    
    #Make transcript record
    tx_record = makeGeneTranscriptFromExons(tx, gene_id)
    annotations = c(annotations, tx_record)
    annotations = c(annotations, tx)
  }
  
  #Construct gene record
  transcripts = annotations[annotations$type == "mRNA"]
  gene = makeGeneTranscriptFromExons(transcripts, "NA", make_gene = TRUE)
  return(c(gene, annotations))
}

makeGeneTranscriptFromExons <- function(exons, parent_id, make_gene = FALSE){
  tx_id = extractFeature(exons$Parent)
  seqname = extractFeature(seqnames(exons))
  
  #Only use the modification field for transcripts
  if (make_gene){
    modification = "NA"
    type = "gene"
  } else {
    modification = extractFeature(exons$modification) 
    type = "mRNA"
  }
  
  ranges = IRanges(start = min(start(exons)), end = max(end(exons)))
  strand = extractFeature(strand(exons))
  metadata = data.frame(modification = modification, Parent = parent_id, 
                        type = type, source = "reviseAnnotations", ID = tx_id, stringsAsFactors = FALSE)
  
  #Create the GRanges obejct
  tx_object = GRanges(seqnames = seqname, ranges = ranges, 
                      strand = strand, seqlengths = seqlengths(exons))
  values(tx_object) = metadata
  return(tx_object)
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

extractTranscriptSets <- function(gene_set, gene_annot){
  #Identify primary transcript and additional transcripts for each gene in the gene_set
  
  #Filter out genes of interes
  genes = gene_annot[gene_annot$gene_id %in% gene_set,]
  
  #Keep only protein coding genes and transcripts
  genes = genes[genes$gene_biotype == "protein_coding" & genes$transcript_biotype == "protein_coding",]
  
  #Iterate over all genes
  transcript_list = list()
  for(gene in unique(genes$gene_id)){
    gene_data = genes[genes$gene_id == gene,]
    
    #Decide on principal isoform
    appris_pos = gene_data[!is.na(gene_data$appris_status),]
    if(nrow(appris_pos) > 0){
      principal_isoform = rownames(appris_pos)[1]
    } else{
      principal_isoform = rownames(gene_data)[1]
    }
    other_isoforms = setdiff(rownames(gene_data), principal_isoform)
    transcript_list[[principal_isoform]] = other_isoforms
  }
  return(transcript_list)
}

makeNewTranscripts <- function(transcript_sets, exons, gene_annot){
  #Create new transcripts based on transcript sets
  principal_transcripts = names(transcript_sets)
  result = GRanges()
  iteration = 1
  for (tx in principal_transcripts){
    print(iteration)
    iteration = iteration + 1
    other_tx = transcript_sets[[tx]]
    transcript_list = constructAlternativeTranscripts(tx, other_tx, exons)
    gene_id = gene_annot[tx,]$gene_id
    transcripts = addTranscriptAnnotations(transcript_list, gene_id)
    result = c(result, transcripts)
  }
  return(result)
}

classifyDifference <- function(tx1_id, tx2_id, annotations, cdss){
  #Compare two trancripts and decide where they differ.
  
  #Extract exons of the transcripts
  exons_tx1 = annotations[annotations$Parent == tx1_id,]
  exons_tx2 = annotations[annotations$Parent == tx2_id,]
  
  changes = indentifyAddedRemovedRegions(exons_tx1, exons_tx2)
  if(is.null(changes)){ #No shared regions between transcripts
    difference = c(1,1,1)
    names(difference) = c("upstream", "downstream", "contained")
    return(NULL)
  }
  else{
    changes = c(changes$addedInSecond, changes$removedFromFirst)
    difference = calulateChangesLength(changes)
    coding = decideCoding(c(tx1_id,tx2_id), cdss, changes)
  }
  return(list(difference  = difference, coding = coding))
}

applyClassifyDifference <- function(txs_list, new_annotations, cdss){
  #Classify all transcripts in a list
  names = names(txs_list)
  
  difference = c()
  coding = c()
  used_names = c()
  for (name in names){
    print(name)
    txs = txs_list[[name]]
    classification = tryCatch(classifyDifference(txs[1], txs[2], new_annotations, cdss), warning = function(w) w)
    if (class(classification)[1] != "simpleWarning"){ #Avoid txs with warnings (no shared exons)
      difference = rbind(difference, classification$difference)
      coding = rbind(coding, classification$coding)
      used_names = c(used_names, name)
    }
  }
  difference = data.frame(difference)
  coding = data.frame(coding)
  rownames(difference) = used_names
  rownames(coding) = used_names
  
  return(list(difference = difference, coding = coding))
}

indentifyAddedRemovedRegions <- function(exon_set1, exon_set2){
  #Identify regions that have been added or removed betweeb two transcripts. 
  #Return NULL if there is no shared region between two transcripts.
  
  #Indentify shared exons
  shared_exons = intersect(exon_set1, exon_set2)
  if(length(shared_exons) == 0){
    warning("No shared exons between two transcripts.")
    return(NULL)
  }
  
  #Idetinfy added or removed regions
  shared_region = union(gaps(shared_exons, start = min(start(shared_exons))), shared_exons)  
  addedInSecond = setdiff(exon_set1, shared_exons)
  addedInSecond = compareGRanges(addedInSecond, shared_region)
  removedFromFirst = setdiff(exon_set2, shared_exons)
  removedFromFirst = compareGRanges(removedFromFirst, shared_region)
  
  return(list(addedInSecond = addedInSecond, removedFromFirst = removedFromFirst))
}

calulateChangesLength <- function(granges){
  vector = c(0,0,0)
  names(vector) = c("upstream", "downstream", "contained")
  vector["upstream"] = as.numeric(sum(width(granges[granges$upstream == 1,])))
  vector["downstream"] = as.numeric(sum(width(granges[granges$downstream == 1,])))
  vector["contained"] = as.numeric(sum(width(granges[granges$contained == 1,])))
  return(vector)
}

filterClassification <- function(classification_list, ratio_cutoff = 0.05){
  #Filter out small upstream and dowstram changes, because they cannot be estimated reliably from RNA-Seq
  classification = classification_list$difference
  coding = classification_list$coding
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
  coding$combined = sign(rowSums(coding))
  return(list(difference = classification, coding = coding))
}

decideCoding <- function(tx_ids, cdss, changes){
  #Identify the original trancripts
  tx_ids = sapply(strsplit(tx_ids, ".", fixed = TRUE), function(x) x[1])
  
  #Contruct CDS regions
  cds1 = cdss[[tx_ids[1] ]]
  cds1 = union(gaps(cds1, start = min(start(cds1))), cds1) 
  cds2 = cdss[[tx_ids[2] ]]
  cds2 = union(gaps(cds2, start = min(start(cds2))), cds2) 
  cds_regions = c(cds1, cds2)

  #Check if the changes overlap with the cds regions
  cds_changes = changes[unique(queryHits(findOverlaps(changes,cds_regions)))]
  result = sign(colSums(as.data.frame(values(cds_changes))))
  return(result)
}

makeClassificationFigure <- function(classification_list){
  type = colSums(sign(a$difference[a$difference$multiple == 0,]))
  colSums(sign(a$difference[a$difference$multiple == 1,]))
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