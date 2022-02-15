# original: https://github.com/kauralasoo/RNAseq_pipeline/blob/master/import_export/makeTxreviseCageMetadata.R

library("rtracklayer")
library("tidyr")
library("purrr")
library("readr")
library("dplyr")
library("optparse")

option_list <- list(
  make_option(c("--grp1"), type="character", default="",
              help="Txrevise group 1 upstream gff file", metavar = "path"),
  make_option(c("--grp2"), type="character", default="",
              help="Txrevise group 2 upstream gff file", metavar = "path"),
  make_option(c("--gene_metadata"), type="character", default="",
              help="TSV metadata for genes", metavar = "path"),
  make_option(c("--out"), type="character", default="",
              help="Gzip filename to write output to", metavar = "path")
)
opt <- parse_args(OptionParser(option_list=option_list))


# Import gene metadata
gene_metadata = read.table(opt$gene_metadata, stringsAsFactors = FALSE, header = TRUE) %>%
  select(-phenotype_id, -phenotype_gc_content, -group_id, -quant_id, -phenotype_length)

# Specify required phenotype metadata columns
required_phenotype_meta_columns = c("phenotype_id","quant_id","group_id","gene_id","chromosome","gene_start",
                                    "gene_end","strand","gene_name","gene_type","gene_version","phenotype_pos")

gff_a = import.gff3(opt$grp1) %>%
  as.data.frame()
gff_b = import.gff3(opt$grp2) %>%
  as.data.frame()

gff = rbind(gff_a, gff_b) %>%
  filter(type=="mRNA") %>%
  rename(phenotype_id=ID) %>%
  select(phenotype_id) %>%
  unique() %>%
  separate(phenotype_id, c("gene_id", "grp", "pos", "transcript"), sep = "\\.", remove = F) %>%
  mutate(quant_id = paste(gene_id, grp, pos, sep = "."),
         group_id = paste(gene_id,      pos, sep = ".")) %>%
  select(phenotype_id, quant_id, group_id, gene_id) %>%
  left_join(gene_metadata, by = "gene_id") %>%
  select(required_phenotype_meta_columns) %>%
  na.omit()

# Save expression matrix
gz2 = gzfile(opt$out, "w")
write.table(gff, gz2, sep = "\t", quote = FALSE, row.names = FALSE)
close(gz2)
