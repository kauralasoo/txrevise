# Requires files:
# - txrevise.img Singularity image (always)
# - processed/input/{annotation}.gtf.gz (always)
# - processed/input/gene_metadata.tsv.gz (always)
# - processed/input/promoters.tsv (cage)
# - processed/input/genes.rds (optional) (cage)

N_BATCHES = 200

#Extract transcript tags from the GTF file
rule extract_tags:
	input:
		gtf = "processed/input/{annotation}.gtf.gz"
	output:
		tags = "processed/{annotation}.transcript_tags.txt"
	threads: 1
	resources:
		mem = 1000
	singularity:
		"./txrevise.img"
	shell:
		"""
		python3 extractTranscriptTags.py --gtf {input.gtf} > {output.tags}
		"""

#Prepare annotationns
rule prepare_annotations:
	input:
		gtf = "processed/input/{annotation}.gtf.gz",
		tags = "processed/{annotation}.transcript_tags.txt"
	output:
		annot = "processed/{annotation}_regular/{annotation}_regular.txrevise_annotations.rds"
	threads: 1
	resources:
		mem = 6000
	singularity:
		"./txrevise.img"
	shell:
		"""
		Rscript prepareAnnotations.R --gtf {input.gtf} --tags {input.tags} --out {output.annot}
		"""

#Create new transcript annotations from promoter annotations and regular txrevise annotations
rule build_cage_annotations:
	input:
		grp1 = "processed/{annotation}_regular/txrevise_regular.grp_1.upstream.gff3",
		grp2 = "processed/{annotation}_regular/txrevise_regular.grp_2.upstream.gff3",
		promoters = "processed/input/promoters.tsv",
		genes = "processed/input/genes.rds"
	output:
		annot = "processed/{annotation}_CAGE-{N}/new_transcripts_{N}.rds"
		genes = "processed/{annotation}_CAGE-{N}/new_transcript_genes_{N}.rds"
	threads: 4
	resources:
		mem = 8000
	singularity:
		"./txrevise.img"
	shell:
		"""
		Rscript --verbose build_cage_annotations.R --grp1 {input.grp1} --grp2 {input.grp2} --promoters {input.promoters} --genes {input.genes} --N {wildcards.N} --output_transcripts {output.annot} --output_genes {output.genes}
		"""

#Prepare CAGE annotations for integration
rule prepare_cage:
	input:
		transcripts = "processed/{annotation}_CAGE-{N}/new_transcripts_{N}.rds",
		annot = "processed/{annotation}_regular/{annotation}_regular.txrevise_annotations.rds"
	output:
		cage_annots = "processed/{annotation}_CAGE-{N}/{annotation}_CAGE-{N}.txrevise_annotations.rds"
	threads: 1
	resources:
		mem = 6000
	singularity:
		"./txrevise.img"
	shell:
		"""
		Rscript prepare_CAGE_data.R --new_transcripts {input.transcripts} --txrevise_annotations {input.annot} --output {output.cage_annots}
		"""

#Construct events for regular annotations
rule construct_events_regular:
	input:
		annot = "processed/{annotation}_regular/{annotation}_regular.txrevise_annotations.rds"
	output:
		"processed/{annotation}_regular/batch/txrevise.grp_1.upstream.{batch}_{n_batches}.gff3",
		"processed/{annotation}_regular/batch/txrevise.grp_2.upstream.{batch}_{n_batches}.gff3",
		"processed/{annotation}_regular/batch/txrevise.grp_1.contained.{batch}_{n_batches}.gff3",
		"processed/{annotation}_regular/batch/txrevise.grp_2.contained.{batch}_{n_batches}.gff3",
		"processed/{annotation}_regular/batch/txrevise.grp_1.downstream.{batch}_{n_batches}.gff3",
		"processed/{annotation}_regular/batch/txrevise.grp_2.downstream.{batch}_{n_batches}.gff3"
	params:
		batch_str = "'{batch} {n_batches}'",
		outdir = "processed/{annotation}_regular/batch"
	threads: 1
	resources:
		mem = 4000
	singularity:
		"./txrevise.img"
	shell:
		"""
		Rscript constructEvents.R --annot {input.annot} --batch {params.batch_str} --out {params.outdir} --fill {config[fill]}
		"""

#Construct events for CAGE
rule construct_events_cage:
	input:
		annot = "processed/{annotation}_regular/{annotation}_regular.txrevise_annotations.rds",
		cage_annots = "processed/{annotation}_CAGE-{N}/{annotation}_CAGE-{N}.txrevise_annotations.rds"
	output:
		"processed/{annotation}_CAGE-{N}/batch/txrevise.grp_1.upstream.{batch}_{n_batches}.gff3",
		"processed/{annotation}_CAGE-{N}/batch/txrevise.grp_2.upstream.{batch}_{n_batches}.gff3",
		"processed/{annotation}_CAGE-{N}/batch/txrevise.grp_1.contained.{batch}_{n_batches}.gff3",
		"processed/{annotation}_CAGE-{N}/batch/txrevise.grp_2.contained.{batch}_{n_batches}.gff3",
		"processed/{annotation}_CAGE-{N}/batch/txrevise.grp_1.downstream.{batch}_{n_batches}.gff3",
		"processed/{annotation}_CAGE-{N}/batch/txrevise.grp_2.downstream.{batch}_{n_batches}.gff3"
	params:
		batch_str = "'{batch} {n_batches}'",
		outdir = "processed/{annotation}_CAGE-{N}/batch"
	threads: 1
	resources:
		mem = 4000
	singularity:
		"./txrevise.img"
	shell:
		"""
		Rscript constructEvents.R --annot {input.annot} --batch {params.batch_str} --out {params.outdir} --fill {config[fill]} --start_end_diff {wildcards.N} --cage {input.cage_annots}
		"""

#Merge txrevise output files
rule merge_gff_files:
	input:
		gff = expand("processed/{{annotation}}_{{kind}}/batch/txrevise.{{group}}.{{position}}.{batch}_{n_batches}.gff3",
			batch = [i for i in range(1, N_BATCHES + 1)],
			n_batches = N_BATCHES)
	output:
		gff = "processed/{annotation}_{kind}/txrevise_{kind}.{group}.{position}.gff3"
	threads: 1
	resources:
		mem = 1000
	singularity:
		"./txrevise.img"
	shell:
		'cat {input.gff} | grep -v "^#" > {output.gff}'

#Build metadata for annotations
rule build_metadata:
	input:
		grp1 = "processed/{annotation}_{kind}/txrevise_{kind}.grp_1.upstream.gff3",
		grp2 = "processed/{annotation}_{kind}/txrevise_{kind}.grp_2.upstream.gff3",
		gene_metadata = "processed/input/gene_metadata.tsv.gz"
	output:
		out = "processed/{annotation}_{kind}/txrevise_{kind}_phenotype_metadata.tsv.gz"
	threads: 1
	resources:
		mem = 4000
	singularity:
		"./txrevise.img"
	shell:
		"""
		Rscript make_annotation_metadata.R --grp1 {input.grp1} --grp2 {input.grp2} --out {output.out} --gene_metadata {input.gene_metadata}
		"""

#Iterate over groups and positions
rule make_one:
	input:
		gff = expand("processed/{{annotation}}_{{kind}}/txrevise_{{kind}}.{group}.{position}.gff3",
			group = ["grp_1", "grp_2"],
			position = ["upstream", "contained", "downstream"]),
		metadata = "processed/{annotation}_{kind}/txrevise_{kind}_phenotype_metadata.tsv.gz"
	output:
		"processed/{annotation}_{kind}/completed.txt"
	threads: 1
	resources:
		mem = 1000
	shell:
		"echo 'Done!' > {output}"

#Create cage annotations for every N
rule make_multiple_cage:
	input:
		N_out = expand("processed/{{annotation}}_CAGE-{N}/completed.txt", N = str(config['Ns']).split(' ')),
		regular_out = "processed/{annotation}_regular/completed.txt" # to get all regular gff files
	output:
		"processed/{annotation}_all_completed.txt"
	threads: 1
	resources:
		mem = 1000
	shell:
		"echo 'Done!' > {output}"
