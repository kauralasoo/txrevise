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
  result[["shared_exons"]] = shared_exons
  return(result)
}