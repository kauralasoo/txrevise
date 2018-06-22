#' Fill missing internal exons into alternative transcript start events
#'
#' @param transcripts list of GRanges objects containing alternative transcript starts
#' @param type either 'start' or 'end' to indicate whether alternative internal exons 
#' should be masked in promoter or 3'end events, respectively.
#'
#' @return New list of alternative transcript starts where the missing exons have been filled in.
#' @author Kaur Alasoo
#' @export
fillMissingInternalExons <- function(transcripts, type){
  #Identify all unique exons
  all_diff = differenceFromShared(transcripts)
  
  #Keep unique exons that are not first or last
  if(type == "start"){
    not_first_or_last = listUnion(purrr::map(all_diff, removeFirstExon))
  } else if(type == "end"){
    not_first_or_last = listUnion(purrr::map(all_diff, removeLastExon))
  } else{
    error("Type has to be either 'start' or end 'end'.")
  }
  
  #Add missing internal exons
  new_transcripts = purrr::map(as.list(transcripts), ~addMissingInternalExons(., not_first_or_last))
  
  #Get rid of newly created redundant transcripts
  non_redundant = mergeByMaxDifference(new_transcripts, max_internal_diff = 0, max_start_end_diff = 0)
  return(non_redundant)
}

#Remove first exon from a GRanges object using the correct strand information
removeFirstExon <- function(granges){
  if(length(granges) > 0){
    gene_strand = as.character(strand(granges))[1]
    if(gene_strand == "+"){
      new_granges = sort(granges)
      new_granges = new_granges[-1]
    } else{
      #First sort decereasing, then increasing
      new_granges = sort(granges, decreasing = TRUE)
      new_granges = new_granges[-1]
      new_granges = sort(new_granges)
    }
    return(new_granges)
  } else {
    return(granges)
  }
}

#Remove first exon from a GRanges object using the correct strand information
removeLastExon <- function(granges){
  if(length(granges) > 0){
    gene_strand = as.character(strand(granges))[1]
    if(gene_strand == "+"){
      new_granges = sort(granges, decreasing = TRUE)
      new_granges = new_granges[-1]
    } else{
      #First sort decereasing, then increasing
      new_granges = sort(granges)
      new_granges = new_granges[-1]
      new_granges = sort(new_granges)
    }
    return(new_granges)
  } else {
    return(granges)
  }
}

constructTranscriptBoundaries <- function(granges){
  gps = GenomicRanges::gaps(granges, start = min(start(granges)))
  boundaries = c(granges, gps) %>% GenomicRanges::reduce() %>% GenomicRanges::sort()
  return(boundaries)
}

addMissingInternalExons <- function(granges, missing_exons){
  boundaries = constructTranscriptBoundaries(granges)
  overlapping_exons = GenomicRanges::intersect(boundaries, missing_exons)
  joint_exons = c(granges, overlapping_exons) %>% GenomicRanges::sort() %>% GenomicRanges::reduce()
  return(joint_exons)
}