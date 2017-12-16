## ---- echo=FALSE---------------------------------------------------------
suppressMessages(library("devtools"))
suppressMessages(library("dplyr"))
suppressMessages(library("wiggleplotr"))
suppressMessages(load_all("../../txrevise/"))

## ------------------------------------------------------------------------
IRF5_data = readRDS("../data/IRF5.rds")
plotting_annotations = dplyr::select(IRF5_data$metadata, ensembl_transcript_id, ensembl_gene_id, external_gene_name, strand) %>% 
  dplyr::rename(transcript_id = ensembl_transcript_id, gene_id = ensembl_gene_id, gene_name = external_gene_name)
knitr::kable(plotting_annotations)

## ----fig.width=6, fig.height=6-------------------------------------------
wiggleplotr::plotTranscripts(IRF5_data$exons, IRF5_data$cdss, plotting_annotations, rescale_introns = TRUE)


## ------------------------------------------------------------------------
missing_ends = dplyr::select(IRF5_data$metadata, ensembl_transcript_id, mRNA_start_NF, mRNA_end_NF, cds_start_NF, cds_end_NF)
knitr::kable(missing_ends)

## ---- fig.width=6, fig.height=6------------------------------------------
gene_extended_tx = txrevise::extendTranscriptsPerGene(IRF5_data$metadata, IRF5_data$exons, IRF5_data$cdss)
gene_data_ext = txrevise::replaceExtendedTranscripts(IRF5_data, gene_extended_tx)
wiggleplotr::plotTranscripts(gene_data_ext$exons, gene_data_ext$cdss, plotting_annotations, rescale_introns = TRUE)

## ---- fig.width=6, fig.height=6------------------------------------------
transcript_groups = txrevise::identifyTranscriptGroups(gene_data_ext$exons)
wiggleplotr::plotTranscripts(transcript_groups$grp_1, rescale_introns = TRUE)


## ---- fig.width=6, fig.height=6------------------------------------------
wiggleplotr::plotTranscripts(transcript_groups$grp_2, rescale_introns = TRUE)

## ------------------------------------------------------------------------
alt_events = txrevise::constructAlternativeEvents(gene_data_ext$exons, "ENSG00000128604")

## ---- fig.width=5, fig.height=4------------------------------------------
wiggleplotr::plotTranscripts(alt_events$ENSG00000128604.grp_1$upstream, transcript_annotations = plotting_annotations, rescale_introns = TRUE)

## ---- fig.width=5, fig.height=4------------------------------------------
wiggleplotr::plotTranscripts(alt_events$ENSG00000128604.grp_1$contained, transcript_annotations = plotting_annotations, rescale_introns = TRUE)

## ---- fig.width=5, fig.height=4------------------------------------------
wiggleplotr::plotTranscripts(alt_events$ENSG00000128604.grp_1$downstream, transcript_annotations = plotting_annotations, rescale_introns = TRUE)


## ---- fig.width=5, fig.height=4------------------------------------------
wiggleplotr::plotTranscripts(alt_events$ENSG00000128604.grp_2$upstream, transcript_annotations = plotting_annotations, rescale_introns = TRUE)

## ---- fig.width=5, fig.height=4------------------------------------------
wiggleplotr::plotTranscripts(alt_events$ENSG00000128604.grp_2$contained, transcript_annotations = plotting_annotations, rescale_introns = TRUE)

## ---- fig.width=5, fig.height=4------------------------------------------
wiggleplotr::plotTranscripts(alt_events$ENSG00000128604.grp_2$downstream, transcript_annotations = plotting_annotations, rescale_introns = TRUE)

