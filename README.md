# _txrevise_
_txrevise_ R package provides utilites to pre-process Ensembl transcript annotations to quantify differences in transcript strucuture (alternative promoters, alternative splicing, alternative poly-adenylation) either between experimental conditions or genotypes (e.g. for transcript usage quntitative trait loci (tuQTL) mapping). 

## Constructing transcription events
This section contains step-by-step instruction for how to construct transcriptional events based from Ensembl transcrtipt annotations.
### Dependencies
Make sure that you have R 3.5 installed together with the following packages:

 - optparse
 - dplyr
 - purrr
 - [rtracklayer](https://bioconductor.org/packages/rtracklayer/)
 - [GenomicFeatures](https://bioconductor.org/packages/GenomicFeatures/)
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
Finally, we can use the `constructEvents.R` script to construct alternative transcription events. Since the implementation of the various algorithms in txrevise have not been optimised for efficiency, processing the full Ensembl GTF file can take several days. However, different genes can trivially processed in parallel. To simplify parallel execution, `constructEvents.R` has the `--batch` option that takes as an input two integers separated by a space. For example, '1 200' would mean that all genes are split into 200 batches and and only genes in the first are be processed. To process all genes, you simply need to iterate from '1 200' to '200 200'. 
	
	Rscript scripts/constructEvents.R --annot Homo_sapiens.GRCh38.92.txrevise_annotations.rds --batch '1 2000' --out txrevise_events  --fill TRUE

### Step 5: Merge output files
Each run of `constructEvents.R` produces up to six output files: alternative promoter, internal exon and 3' end events (labeled as upstream, contained and downstream) for two possible sets of shared exons (grp_1 and grp_2). See the [vignette](http://htmlpreview.github.io/?https://github.com/kauralasoo/txrevise/blob/master/inst/doc/construct_events.html) for more details on how the events are constructed.

To be able to use these events with transcript quantification software such as [Salmon](http://salmon.readthedocs.io/en/latest/) or [kallisto](https://pachterlab.github.io/kallisto/), we first need to merge all files from different batches:

	cat txrevise_events/txrevise.grp_1.upstream.* | grep -v "^#" > txrevise.grp_1.upstream.gff3
	cat txrevise_events/txrevise.grp_1.contained.* | grep -v "^#" > txrevise.grp_1.contained.gff3
	cat txrevise_events/txrevise.grp_1.downstream.* | grep -v "^#" > txrevise.grp_1.downstream.gff3
	
	cat txrevise_events/txrevise.grp_2.upstream.* | grep -v "^#" > txrevise.grp_2.upstream.gff3
	cat txrevise_events/txrevise.grp_2.contained.* | grep -v "^#" > txrevise.grp_2.contained.gff3
	cat txrevise_events/txrevise.grp_2.downstream.* | grep -v "^#" > txrevise.grp_2.downstream.gff3

### Step 6: Quantify event expression
The GFF3 files constructed by _txrevise_ can be used with any transcript quantification algorithm (e.g. [Salmon](http://salmon.readthedocs.io/en/latest/) or [kallisto](https://pachterlab.github.io/kallisto/)). **Importantly, since the promoter, internal exon and 3' end events constructed by txrevise share some of the same exons, they should always quantified independently**. For example, if you have six GFF3 files from txrevise (txrevise.grp_1.upstream.gff3, txrevise.grp_1.contained.gff3, txrevise.grp_1.downstream.gff3, txrevise.grp_2.upstream.gff3, txrevise.grp_2.contained.gff3, txrevise.grp_2.downstream.gff3), you should also run the quantification algorithm six times with each of the GFF3 file independently.

## Pre-computed transcript annotations
Running _txrevise_ on the latest version of Ensembl can be quite timeconsuming. Thus, to make it easier to get started, we have pre-computed alternatve transcription events in the GFF3 format for both GRCh37 and GRCh38 reference genomes:

### Current transcription events
TODO

### Deprecated transcription events
We previously made pre-computed sets of transcription events available here, but these are now outdated and we strongly recommend to construct new events using the instructions above.
* [GRCh38 + Ensembl 87](https://zenodo.org/record/997492#.Wcqa3tMjHOQ)
* [GRCh37(hg19) + Ensembl 90](https://zenodo.org/record/997251#.Wco2Q9MjHUJ)
<!--stackedit_data:
eyJoaXN0b3J5IjpbMTczNzQwOTQ4OCwyMDAxOTE1NTA1LDE0Mj
k3ODk5NzIsMTU1NDYyOTIyMSwxNjQxOTM2Mzk5LDc1NjI1MDcw
LC0xMzU0MjI0NTAsLTE0MDcxMjc3MTUsMTY1MzMxOTM2NSwtMT
Y1NTA0MDQzOCwtODg0MjM4NjMzLC0yMDAzNDA1NjM5LDE1MDgx
OTU4MzVdfQ==
-->