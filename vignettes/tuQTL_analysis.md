# Transcript usage QTL analysis with txrevise and qtlmap

## Quantification

You can use the [eQTL-Catalogue/rnaseq](https://github.com/eQTL-Catalogue/rnaseq) pipeline to quantify the expression of txrevise events. Here's an example command for the GEUVADIS dataset:

```
NXF_VER=18.10.1 nextflow run main.nf\
 -profile eqtl_catalogue\
 --readPathsFile /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/readPaths/GEUVADIS_readPaths.tsv\
 --unstranded\
 --skip_qc\
 --skip_multiqc\
 --skip_stringtie\
 --saveAlignedIntermediates\
 --run_txrevise\
 --txrevise_gffs '/gpfs/hpc/projects/genomic_references/annotations/txrevise/Homo_sapiens.GRCh38.96_CAGE_10bp/*.gff3'\
 -resume\
 -executor.queueSize 100
 ```
 
## Normalisation

Txrevise event usage data can be normalised using the kerimoff/qcnorm pipeline. Make sure that you have the txrevise phenotype metdata file file. If you are using the eQTL Catalogue [txrevise annotations](https://zenodo.org/record/3366280#.XqnnZJMzaGh), then you can also use our [phenotype metadata files](https://zenodo.org/record/3366011#.XqnnPpMzaGg). 

```
nextflow run normalisation.nf -profile tartu_hpc -resume\
 --study_name GEUVADIS_EUR\
 --quant_results_path /gpfs/hpc/home/a72094/projects/rnaseq/results/\
 --sample_meta_path /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/studies/cleaned/GEUVADIS_EUR.tsv\
 --txrev_pheno_meta_path /gpfs/hpc/projects/genomic_references/annotations/txrevise/Homo_sapiens.GRCh38.96_CAGE_10bp/txrevise_Ensembl_96_CAGE_10bp_phenotype_metadata.tsv.gz\
 --skip_exon_norm\
 --skip_tx_norm\
 --skip_leafcutter_norm\
 --outdir GEUVADIS_EUR
```

## QTL analysis


