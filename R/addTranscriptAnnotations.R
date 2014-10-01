
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
  table = c(gene, annotations)
  
  #Add gene id to each entry for easy filtering
  val = values(table)
  val$gene_id = gene_id
  values(table) = val
  return(table)
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

transcriptsToAnnotations <- function(transcript_list, gene_tx_map){
  #Convert a list of transcripts into annotations that can be used by 
  map_filtered = dplyr::filter(tx_gene_map, transcript_id %in% names(transcript_list))
  gene_ids = unique(map_filtered$gene_id)
  annotations = GRanges()
  
  for(gene in gene_ids){
    tx_ids = dplyr::filter(map_filtered, gene_id == gene)$transcript_id
    tx_annot = addTranscriptAnnotations(transcript_list[tx_ids], gene)
    annotations = c(annotations, tx_annot)
  }
  return(annotations)
}