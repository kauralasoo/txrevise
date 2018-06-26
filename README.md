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
First, you need to download the GTF file from the Ensembl website. For example, if you want to use Ensembl version 92, you should download the following file:

	wget ftp://ftp.ensembl.org/pub/release-92/gtf/homo_sapiens/Homo_sapiens.GRCh38.92.gtf.gz

### Step 2: Extract tanscript tags from the GTF file
Ensembl GTF file contains a tags field marking protein coding transcript that are truncated at 5' or 3' ends. Txrevise uses this information to extend truncated transcripts. Unfortunately, the `import` function from rtracklayer does not import the tags field correctly. Thus, we need to first extract the transcript tags manually using a custom Python 3 script that comes with txrevise.

	python scripts/extractTranscriptTags.py --gtf Homo_sapiens.GRCh38.92.gtf.gz > Homo_sapiens.GRCh38.92.transcript_tags.txt

### Step 3: Prepare transcript annotations for event construction
Next, we need to convert the transcript annotations (in GTF format) to a binary representation that can be efficiently used by txrevise.

	Rscript scripts/prepareAnnotations.R --gtf Homo_sapiens.GRCh38.92.gtf.gz --tags Homo_sapiens.GRCh38.92.transcript_tags.txt --out Homo_sapiens.GRCh38.92.txrevise_annotations.rds

### Step 4: Construct transcription events
Finally, we can use the constructEvents.R script to construct alternative transcription events. Since 
	
	Rscript scripts/constructEvents.R --annot Homo_sapiens.GRCh38.92.txrevise_annotations.rds --batch '1 2000' --out txrevise_events  --fill TRUE

## Getting started
See the [vignette](http://htmlpreview.github.io/?https://github.com/kauralasoo/txrevise/blob/master/inst/doc/construct_events.html) for examples of the basic functionality.

## Pre-computed transcript annotations
Running _txrevise_ on the latest version of Ensembl can be quite timeconsuming. Thus, to make it easier to get started, we have pre-computed alternatve transcription events in the GFF3 format for both GRCh37 and GRCh38 reference genomes:
* [GRCh38 + Ensembl 87](https://zenodo.org/record/997492#.Wcqa3tMjHOQ)
* [GRCh37(hg19) + Ensembl 90](https://zenodo.org/record/997251#.Wco2Q9MjHUJ)
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTc5MjU3ODQ5NiwyMDMwNDc5MDAzLDE1MD
gxOTU4MzVdfQ==
-->