
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
    #Calculated the length of transcribed difference between two transcripts
    diff = calulateChangesLength(changes)
    abs_diff = calulateChangesLength(changes, abs_diff = TRUE)
    transcribed_diff = data.frame(tx1_id = tx1_id, tx2_id = tx2_id, 
                         position = names(diff), type = "transcribed", 
                         diff = diff, abs_diff = abs_diff)
    diff_df = dplyr::filter(transcribed_diff, diff > 0)
    
    #Calculate coding changes as well
    if(!is.null(cdss)){
      coding_changes = intersectCoding(tx1_id, tx2_id, changes, cdss)
      if(!is.null(coding_changes)){ #Proceed only if CDS information available for both transcripts
        coding_diff = calulateChangesLength(coding_changes)
        coding_abs_diff = calulateChangesLength(coding_changes, abs_diff = TRUE)
        cds_diff = data.frame(tx1_id = tx1_id, tx2_id = tx2_id, 
                              position = names(diff), type = "coding", 
                              diff = coding_diff, abs_diff = coding_abs_diff)
        cds_df = dplyr::filter(cds_diff, diff > 0)
        
        #Bind all results together
        diff_df = rbind(diff_df, cds_df)
      }
    }
    return(diff_df)
  }
}

applyClassifyDifference <- function(txs_list, new_annotations, cdss = NULL, printProgress = FALSE){
  #Classify all transcripts in a list
  
  #Prepare results data.frame
  result = data.frame()
  
  #Iterate over genes
  names = names(txs_list)
  for (name in names){
    if(printProgress){ print(name) }
    txs = txs_list[[name]]
    classification = classifyDifference(txs[1], txs[2], new_annotations, cdss)
    if (!is.null(classification)){ #Avoid genes with no overlaping transcripts
      classification = dplyr::mutate(classification, gene_name = name)
      result = rbind(result, classification)
    }
  }
  return(result)
}



calulateChangesLength <- function(changes, abs_diff = FALSE){ 
  #Calculate the number of bases that differ between two transcripts.
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