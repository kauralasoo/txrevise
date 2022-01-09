### Written for Homo_sapiens.GRCh38.96, but other versions are analogous

# Requires files:
# - txrevise.img - Singularity image from build_container.sh
# - processed/input/Homo_sapiens.GRCh38.96.gtf.gz - genome annotation from Ensembl
# - processed/input/gene_metadata.tsv.gz - gene metadata, such as https://zenodo.org/record/3366011/files/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz .

# Fill internal exons for promoter and 3'end events

snakemake -p processed/Homo_sapiens.GRCh38.96_regular/completed.txt --cluster ./snakemake_submit_UT.py --jobs 200 --config fill=TRUE --use-singularity
# Revised transcripts in processed/Homo_sapiens.GRCh38.96_regular/

# Fill internal exons for promoter and 3'end events and generate additional transcripts based on promoter annotations, using minimum distances between new transcripts of 10, 20, 50 and 100 base pairs
# Requires in addition:
# - processed/input/promoters.tsv - promoter annotations, such as created in https://github.com/andreasvija/cage/blob/master/qtlmap_prep/clean_promoters.R

snakemake -p processed/Homo_sapiens.GRCh38.96_all_completed.txt --cluster ./snakemake_submit_UT.py --jobs 200 --config fill=TRUE Ns="10 20 50 100" --use-singularity
# Revised transcripts in processed/Homo_sapiens.GRCh38.96_CAGE-N/ for N in Ns
