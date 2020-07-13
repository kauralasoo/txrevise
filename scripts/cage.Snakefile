N_BATCHES = 200
annotation = config["annotation"]


#Iterate over groups and positions
rule make_all:
	input:
		gff = expand("processed/{{annotation}}/merged/txrevise.{group}.{position}.gff3",
			group = ["grp_1", "grp_2"],
			position = ["upstream", "contained", "downstream"])
	output:
		"processed/{annotation}_log.txt"
	threads: 1
	resources:
		mem = 1000
	shell:
		"echo 'Done!' > {output}"


#Extract trascript tags from the GTF file
rule extract_tags:
	input:
		gtf = "processed/{annotation}.gtf.gz"
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
		gtf = "processed/{annotation}.gtf.gz",
		tags = "processed/{annotation}.transcript_tags.txt"
	output:
		annot = "processed/{annotation}.txrevise_annotations.rds"
	threads: 1
	resources:
		mem = 6000
	singularity:
		"./txrevise.img"
	shell:
		"""
		Rscript prepareAnnotations.R --gtf {input.gtf} --tags {input.tags} --out {output.annot}
		"""

#Prepare CAGE annotations for integration
rule prepare_cage:
	input:
		annot = "processed/{annotation}.txrevise_annotations.rds",
		transcripts = "../data/new_transcripts_25.rds"
	output:
		cage_annots = "processed/CAGE_promoter_annotations_25.rds"
	threads: 1
	resources:
		mem = 6000
	singularity:
		"./txrevise.img"
	shell:
		"""
		Rscript prepare_CAGE_data.R --txrevise_annotations {input.annot}
		"""

#Construct events
rule construct_events:
	input:
		annot = "processed/{annotation}.txrevise_annotations.rds",
		cage_annots = "processed/CAGE_promoter_annotations_25.rds"
	output:
		"processed/{annotation}/txrevise.grp_1.upstream.{batch}_{n_batches}.gff3",
		"processed/{annotation}/txrevise.grp_2.upstream.{batch}_{n_batches}.gff3",
		"processed/{annotation}/txrevise.grp_1.contained.{batch}_{n_batches}.gff3",
		"processed/{annotation}/txrevise.grp_2.contained.{batch}_{n_batches}.gff3",
		"processed/{annotation}/txrevise.grp_1.downstream.{batch}_{n_batches}.gff3",
		"processed/{annotation}/txrevise.grp_2.downstream.{batch}_{n_batches}.gff3"
	params:
		batch_str = "'{batch} {n_batches}'",
		outdir = "processed/{annotation}"
	threads: 1
	resources:
		mem = 4000
	singularity:
		"./txrevise.img"
	shell:
		"""
		Rscript constructEvents.R --annot {input.annot} --batch {params.batch_str} --out {params.outdir} --fill {config[fill]} --start_end_diff {config[start_end_diff]} --cage {input.cage_annots}
		"""

#Merge txrevise output files
rule merge_gff_files:
	input:
		gff = expand("processed/{{annotation}}/txrevise.{{group}}.{{position}}.{batch}_{n_batches}.gff3",
			batch = [i for i in range(1, N_BATCHES + 1)],
			n_batches = N_BATCHES)
	output:
		gff = "processed/{annotation}/merged/txrevise.{group}.{position}.gff3"
	threads: 1
	resources:
		mem = 1000
	singularity:
		"./txrevise.img"
	shell:
		'cat {input.gff} | grep -v "^#" > {output.gff}'
