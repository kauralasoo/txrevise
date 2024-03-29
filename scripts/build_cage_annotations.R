library("data.table") # %like%
library("rtracklayer")
library("GenomicRanges")
library("GenomicFeatures")

library("dplyr")
library("readr")
library("optparse")

option_list <- list(
  make_option(c("--grp1"), type="character", default="",
              help="Txrevise group 1 upstream gff file", metavar = "path"),
  make_option(c("--grp2"), type="character", default="",
              help="Txrevise group 2 upstream gff file", metavar = "path"),
  make_option(c("--promoters"), type="character", default="",
              help="TSV file with annotations of promoters", metavar = "path"),
  make_option(c("--N"), type="integer", default=20,
              help="Integer specifying how far new exon starts need to be
              from exixting and new annotations", metavar = "integer"),
  make_option(c("--output_transcripts"), type="character", default="",
              help="RDS file to where to write the new annotations", metavar = "path"),
  make_option(c("--output_genes"), type="character", default="",
              help="RDS file to where to write the gene ID-s of the new annotations", metavar = "path")
)
opt <- parse_args(OptionParser(option_list=option_list))

PROMOTER_BUFFER_RANGE = opt$N
EXON_START_BUFFER_RANGE = PROMOTER_BUFFER_RANGE
EXON_BACK_RANGE = 1000

promoter_annots = read_tsv(opt$promoters, col_types="ccciiicii")

upstream1 = GenomicFeatures::makeTxDbFromGFF(opt$grp1)
upstream2 = GenomicFeatures::makeTxDbFromGFF(opt$grp2)
exons_list1 = GenomicFeatures::exonsBy(upstream1, by = "tx", use.names = TRUE)
exons_list2 = GenomicFeatures::exonsBy(upstream2, by = "tx", use.names = TRUE)


#for every gene

all_new_transcripts = list()
new_transcript_genes = c()

start_time = Sys.time()
for (gene in unique(promoter_annots$gene_name)) {
  gene_new_transcripts = list()

  #for every promoter belonging to the gene

  promoters = promoter_annots %>%
    filter(gene_name == gene)

  annot1 = exons_list1[names(exons_list1) %like% paste0(gene, ".grp_1.upstream")]
  annot2 = exons_list2[names(exons_list2) %like% paste0(gene, ".grp_2.upstream")]
  exons = c(annot1, annot2)

  new_transcript_id = 0

  utilized_promoters = list() # have been used for creating transcripts
  for (i in 1 : dim(promoters)[1]) {

    peak_start = promoters$peak_start[i]
    peak_end = promoters$peak_end[i]
    strand = promoters$strand[i]
    chr = promoters$chr[i]

    #if promoter +- PROMOTER_BUFFER_RANGE does not overlap with an already utilized promoter

    promoter = GRanges(seqnames=chr, IRanges(start=peak_start, end=peak_end))
    extended_promoter = GRanges(seqnames=chr, IRanges(start=peak_start - PROMOTER_BUFFER_RANGE,
                                                    end=peak_end + PROMOTER_BUFFER_RANGE))

    intersections = lapply(utilized_promoters, pintersect, extended_promoter, drop.nohit.ranges=TRUE)
    if (sum(unlist(lapply(intersections, length))) != 0) {next}

    #for every exon in transcripts

    if (length(exons) < 1) {next}

    a = peak_start - EXON_START_BUFFER_RANGE
    b = peak_end + EXON_START_BUFFER_RANGE

    overlaps_with_exon_start = FALSE
    for (j in 1 : length(exons@unlistData)) {

      #if promoter +- EXON_START_BUFFER_RANGE does not overlap with an exon's start

      exon_start = exons@unlistData@ranges@start[j]
      if (strand == "-") {
        exon_start = exon_start + exons@unlistData@ranges@width[j]
      }

      if (a < exon_start & b > exon_start) {
        overlaps_with_exon_start = TRUE
        break
      }

    }
    if (overlaps_with_exon_start) {next}

    #for every transcript belonging to the gene

    for (j in 1 : length(ranges(exons))) {
      transcript = unlist(ranges(exons)[j])

      #for every exon belonging to the transcript

      for (k in 1 : length(transcript)) {
        exon_start = transcript@start[k]
        exon_end = transcript@start[k] + transcript@width[k]

        #if promoter overlapping exon + EXON_BACK_RANGE bp back

        if (strand == "+") {exon_start = exon_start - EXON_BACK_RANGE}
        if (strand == "-") {exon_end = exon_end + EXON_BACK_RANGE}
        if ((peak_start > exon_end) | (peak_end < exon_start)) {next}

        #create a new exon from start of promoter to end of overlapping exon

        if (strand == "+") {
          new_exon_start = peak_start
          new_exon_end = exon_end
        }
        if (strand == "-") {
          new_exon_start = exon_start
          new_exon_end = peak_end
        }

        new_exon = IRanges(start=new_exon_start, end=new_exon_end)

        #take all exons in the transcript after the overlapping exon and create new transcript

        following = IRanges()

        if (strand == "+") {
          if (k < length(transcript)) {
            following = unname(transcript[(k+1) : length(transcript)])
          }
          new_transcript = c(new_exon, following)
        }

        if (strand == "-") {
          if (k > 1) {
            following = unname(transcript[1 : (k-1)])
          }
          new_transcript = c(following, new_exon)
        }

        ns = 1 : length(new_transcript)
        exon_rank = ns
        if (strand == "-") {exon_rank = rev(exon_rank)}

        new_transcript_id = new_transcript_id + 1
        exon_names = paste0(gene, ".new.upstream.", new_transcript_id, ":", ns)

        new_transcript = GRanges(seqnames=chr, ranges=new_transcript, strand=strand,
                                 exon_id=ns, exon_name=exon_names, exon_rank=exon_rank)
        new_is_duplicate = FALSE
        if (length(gene_new_transcripts) > 0) {
          new_is_duplicate = any(sapply(new_transcript %in% gene_new_transcripts, all))
        }


        if (!new_is_duplicate) {
          gene_new_transcripts = c(gene_new_transcripts, list(new_transcript))
          new_transcript_name = paste(gene, "new", "upstream", length(all_new_transcripts)+1, sep=".")
          all_new_transcripts[new_transcript_name] = new_transcript
          new_transcript_genes = c(new_transcript_genes, gene)
        }

      }
    }

    utilized_promoters = c(utilized_promoters, list(promoter))
  }
}

Sys.time() - start_time
print("Created")
length(all_new_transcripts)
print("transcripts over")
new_transcript_genes = unique(new_transcript_genes)
length(new_transcript_genes)
print("genes using distance")
PROMOTER_BUFFER_RANGE

all_new_transcripts = GRangesList(all_new_transcripts)
saveRDS(all_new_transcripts, opt$output_transcripts)
saveRDS(new_transcript_genes, opt$output_genes)
