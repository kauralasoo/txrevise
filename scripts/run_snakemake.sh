### Homo_sapiens.GRCh38.96, other versions are analogous

# Fill internal exons for promoter and 3'end events
snakemake -p processed/Homo_sapiens.GRCh38.96_regular/completed.txt --cluster ./snakemake_submit_UT.py --jobs 200 --config fill=TRUE --use-singularity
# Revised transcripts in processed/Homo_sapiens.GRCh38.96_regular/

# Generate additional transcripts based on promoter annotations
snakemake -p processed/Homo_sapiens.GRCh38.96_all_completed.txt --cluster ./snakemake_submit_UT.py --jobs 200 --config fill=TRUE Ns="20" --use-singularity
# Revised transcripts in processed/Homo_sapiens.GRCh38.96_CAGE-20/
