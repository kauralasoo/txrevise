N_BATCHES = 200

#Extract trascript tags from the GTF file
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
		annot = "processed/{annotation}_{kind}/{annotation}_regular.txrevise_annotations.rds"
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
		transcripts = "processed/input/new_transcripts_{N}.rds",
		annot = "processed/{annotation}_CAGE_{N}/{annotation}_regular.txrevise_annotations.rds"
	output:
		cage_annots = "processed/{annotation}_CAGE_{N}/{annotation}_CAGE_{N}.txrevise_annotations.rds"
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
		annot = "processed/{annotation}.txrevise_annotations.rds",
		cage_annots = "processed/{annotation}.CAGE_promoter_annotations_25.rds"
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
		annot = "processed/{annotation}.txrevise_annotations.rds",
		cage_annots = "processed/{annotation}.CAGE_promoter_annotations_25.rds"
	output:
		"processed/{annotation}_CAGE_{N}/batch/txrevise.grp_1.upstream.{batch}_{n_batches}.gff3",
		"processed/{annotation}_CAGE_{N}/batch/txrevise.grp_2.upstream.{batch}_{n_batches}.gff3",
		"processed/{annotation}_CAGE_{N}/batch/txrevise.grp_1.contained.{batch}_{n_batches}.gff3",
		"processed/{annotation}_CAGE_{N}/batch/txrevise.grp_2.contained.{batch}_{n_batches}.gff3",
		"processed/{annotation}_CAGE_{N}/batch/txrevise.grp_1.downstream.{batch}_{n_batches}.gff3",
		"processed/{annotation}_CAGE_{N}/batch/txrevise.grp_2.downstream.{batch}_{n_batches}.gff3"
	params:
		batch_str = "'{batch} {n_batches}'",
		outdir = "processed/{annotation}_CAGE_{N}/batch"
		n = "{N}"
	threads: 1
	resources:
		mem = 4000
	singularity:
		"./txrevise.img"
	shell:
		"""
		Rscript constructEvents.R --annot {input.annot} --batch {params.batch_str} --out {params.outdir} --fill {config[fill]} --start_end_diff {n} --cage {input.cage_annots}
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
		grp2 = "processed/{annotation}_{kind}/txrevise_{kind}.grp_2.upstream.gff3"
	output:
		out = "processed/{annotation}_{kind}/txrevise_{kind}_phenotype_metadata.tsv.gz"
	threads: 1
	resources:
		mem = 4000
	singularity:
		"./txrevise.img"
	shell:
		"""
		Rscript make_annotation_metadata.R --grp1 {input.grp1} --grp2 {input.grp2} --out {output.out}
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
		N_out = expand("processed/{{annotation}}_CAGE_{N}/completed.txt", N = config['Ns'].split(' '))
	output:
		"processed/{annotation}_all_completed.txt"
	threads: 1
	resources:
		mem = 1000
	shell:
		"echo 'Done!' > {output}"
