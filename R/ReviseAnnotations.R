
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

classifyDifference <- function(tx1_id, tx2_id, annotations, cdss = NULL){
  #Compare two trancripts and decide where they differ.
  
  #Create exons list to idenitify added and removed regions
  exons_list = list()
  exons_list[[tx1_id]] = annotations[annotations$Parent == tx1_id,]
  exons_list[[tx2_id]] = annotations[annotations$Parent == tx2_id,]
  changes = indentifyAddedRemovedRegions(tx1_id, tx2_id, exons_list)
  if(is.null(changes)){ #No shared regions between transcripts
    return(NULL)
  }
  else{
    result = list(transcribed = NULL, coding = NULL)
    diff = calulateChangesLength(changes)
    abs_diff = calulateChangesLength(changes, abs_diff = TRUE)
    result$transcribed = list(diff = diff, abs_diff = abs_diff)
    
    #Calculate coding changes as well
    if(!is.null(cdss)){
      coding_changes = intersectCoding(tx1_id, tx2_id, changes, cdss)
      if(!is.null(coding_changes)){ #Proceed only if CDS information available for both transcripts
        coding_diff = calulateChangesLength(coding_changes)
        coding_abs_diff = calulateChangesLength(coding_changes, abs_diff = TRUE)
        result$coding = list(diff = coding_diff, abs_diff = coding_abs_diff)
      }
    }
    return(result)
  }
}

applyClassifyDifference <- function(txs_list, new_annotations, cdss = NULL, printProgress = FALSE){
  #Classify all transcripts in a list
  names = names(txs_list)
  
  #Prepare results list
  result = list(transcribed = list(diff = c(), abs_diff = c()))
  if(!is.null(cdss)){ result[["coding"]] = list(diff = c(), abs_diff = c()) }
  used_names = c()
  
  for (name in names){
    if(printProgress){ print(name) }
    txs = txs_list[[name]]
    classification = tryCatch(classifyDifference(txs[1], txs[2], new_annotations, cdss), warning = function(w) w)
    if (class(classification)[1] != "simpleWarning"){ #Avoid txs with warnings (no shared exons)
      used_names = c(used_names, name)
      #Bind trancribed changes
      result$transcribed$diff = rbind(result$transcribed$diff, classification$transcribed$diff)
      result$transcribed$abs_diff = rbind(result$transcribed$abs_diff, classification$transcribed$abs_diff)     
      #Bind coding changes if cdss is present
      if(!is.null(cdss)){
        result$coding$diff = rbind(result$coding$diff, classification$coding$diff)
        result$coding$abs_diff = rbind(result$coding$abs_diff, classification$coding$abs_diff)        
      }
    }
  }
  #Convert to data.frame and add rownames
  result$transcribed$diff = as.data.frame(result$transcribed$diff)
  result$transcribed$abs_diff = as.data.frame(result$transcribed$abs_diff)
  rownames(result$transcribed$diff) = used_names
  rownames(result$transcribed$abs_diff) = used_names
  
  #Do the same for coding changes
  if(!is.null(cdss)){ 
    result$coding$diff = as.data.frame(result$coding$diff)
    result$coding$abs_diff = as.data.frame(result$coding$abs_diff)
    rownames(result$coding$diff) = used_names
    rownames(result$coding$abs_diff) = used_names 
  }
  
  return(result)
}

indentifyAddedRemovedRegions <- function(tx1_id, tx2_id, exons_list){
  #Identify regions that have been added or removed betweeb two transcripts. 
  #Return NULL if there is no shared region between two transcripts.
  
  #Extract exons of the transcripts
  exon_set1 = exons_list[[tx1_id]]
  exon_set2 = exons_list[[tx2_id]]
  
  #Indentify shared exons
  shared_exons = intersect(exon_set1, exon_set2)
  if(length(shared_exons) == 0){
    warning("No shared exons between two transcripts.")
    return(NULL)
  }
  
  #Idetinfy added or removed regions
  shared_region = union(gaps(shared_exons, start = min(start(shared_exons))), shared_exons)  
  addedInSecond = setdiff(exon_set2, shared_exons)
  addedInSecond = compareGRanges(addedInSecond, shared_region)
  removedFromFirst = setdiff(exon_set1, shared_exons)
  removedFromFirst = compareGRanges(removedFromFirst, shared_region)
  
  #Build list with tx ids as names
  result = list()
  result[[tx1_id]] = removedFromFirst
  result[[tx2_id]] = addedInSecond
  return(result)
}

calulateChangesLength <- function(changes, abs_diff = FALSE){ 
  vector = c(0,0,0)
  names(vector) = c("upstream", "downstream", "contained")
  addedInSecond = changes[[2]]
  removedFromFirst = changes[[1]]
  
  if(!abs_diff){
    granges = c(addedInSecond, removedFromFirst)
    vector["upstream"] = as.numeric(sum(width(granges[granges$upstream == 1,])))
    vector["downstream"] = as.numeric(sum(width(granges[granges$downstream == 1,])))
    vector["contained"] = as.numeric(sum(width(granges[granges$contained == 1,])))
  }
  else{
    vector["upstream"] = as.numeric(sum(width(addedInSecond[addedInSecond$upstream == 1,]))) - 
      as.numeric(sum(width(removedFromFirst[removedFromFirst$upstream == 1,])))
    vector["downstream"] = as.numeric(sum(width(addedInSecond[addedInSecond$downstream == 1,]))) -
      as.numeric(sum(width(removedFromFirst[removedFromFirst$downstream == 1,])))
    vector["contained"] = as.numeric(sum(width(addedInSecond[addedInSecond$contained == 1,]))) -
      as.numeric(sum(width(removedFromFirst[removedFromFirst$contained == 1,])))
  }
  return(vector)
}

intersectCoding <- function(tx1_id, tx2_id, changes, cdss){
  #Intersect changes with known CDS to find the changes that affect coding sequence
    
  #Contruct CDS regions
  intersection = intersect(c(tx1_id, tx2_id),names(cdss))
  if(length(intersection) < 2){
    warning("Some transcripts not found in the cdss database")
    return(NULL)
  }

  cds1 = cdss[[tx1_id]]
  cds1 = union(gaps(cds1, start = min(start(cds1))), cds1) 
  cds2 = cdss[[tx2_id]]
  cds2 = union(gaps(cds2, start = min(start(cds2))), cds2) 

  #Find changes affecting the coding sequence
  coding_changes = list()
  tx1_coding = intersect(changes[[tx1_id]], cds1)
  tx2_coding = intersect(changes[[tx2_id]], cds2)
  
  #Put back metadata for overlapping exons
  tx1_metadata = mcols(changes[[tx1_id]][subjectHits(findOverlaps(tx1_coding, changes[[tx1_id]])),])
  tx2_metadata = mcols(changes[[tx2_id]][subjectHits(findOverlaps(tx2_coding, changes[[tx2_id]])),])
  mcols(tx1_coding) = tx1_metadata
  mcols(tx2_coding) = tx2_metadata
  
  #Return
  coding_changes[[tx1_id]] = tx1_coding
  coding_changes[[tx2_id]] = tx2_coding
  return(coding_changes)
}