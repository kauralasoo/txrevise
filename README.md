# _txrevise_
_txrevise_ R package provides utilites to pre-process Ensembl transcript annotations to quantify differences in transcript strucuture (alternative promoters, alternative splicing, alternative poly-adenylation) either between experimental conditions or genotypes (e.g. for transcript usage quntitative trait loci (tuQTL) mapping). 

## Constructing transcription events
This section contains step-by-step instruction for how to construct transcriptional events based from Ensembl transcrtipt annotations.
### Dependencies
Make sure that you have R 3.5 installed together with the following packages:

 - optparse
 - dplyr
 - purrr
 - rtracklayer
 - GenomicFeatures
 - txrevise

You can install _txrevise_ directly from GitHub using the following command:

	library("devtools")
	install_github("kauralasoo/txrevise")

### Step 1: Download the GTF file
First, you need to download the GTF file from the Ensembl website. For example, if you want to use Enssembl

## Getting started
See the [vignette](http://htmlpreview.github.io/?https://github.com/kauralasoo/txrevise/blob/master/inst/doc/construct_events.html) for examples of the basic functionality.

## Pre-computed transcript annotations
Running _txrevise_ on the latest version of Ensembl can be quite timeconsuming. Thus, to make it easier to get started, we have pre-computed alternatve transcription events in the GFF3 format for both GRCh37 and GRCh38 reference genomes:
* [GRCh38 + Ensembl 87](https://zenodo.org/record/997492#.Wcqa3tMjHOQ)
* [GRCh37(hg19) + Ensembl 90](https://zenodo.org/record/997251#.Wco2Q9MjHUJ)
<!--stackedit_data:
eyJoaXN0b3J5IjpbNDkzNDkzNzUwLC0xODMzMTgyODM0LC0xNT
g1ODk1MDQ4XX0=
-->