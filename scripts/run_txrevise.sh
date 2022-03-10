#!/bin/bash

#SBATCH -p amd
#SBATCH -J integrate_annots
#SBATCH -t 3-00:00:00
#SBATCH -c 4
#SBATCH --mem=8G

ANNOTATION=Homo_sapiens.GRCh38.105
IFS='.' read -r -a array <<< "$ANNOTATION"

module load any/singularity/3.5.3
module load squashfs/4.4

OUTPUT=processed/${ANNOTATION}_all_completed.txt
rm $OUTPUT
snakemake -p $OUTPUT --cluster ./snakemake_submit_UT.py --jobs 200 --config fill=TRUE Ns="25" --use-singularity
