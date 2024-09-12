#!/usr/bin/env python

# An Automated Pipeline to Transcriptome Assembly and Quality Assessment
#
# Author: Felipe Vaz Peres
#
# Preparation:
# 1) Setup 'config.yaml' file  with softwares path
# 2) Create 'samples_{genotype}.csv' file and add SRA identifiers (e.g SRR5258954,SRR5258955,SRR5258994,SRR5258995)
# 3) Create 'parts.csv' file and add the value of parts you want the kraken file to be split (e.g 00,01,02 for 3 equal parts)
# 4) Setup 'GENOTYPE' variable with the genotype name (e.g GENOTYPE=QN05-1509)
#
# Usage (remove the -n to dont show dry-run and start the jobs):
# 1) Load modules: BUSCO/3.0; transrate/1.0.3
# 2) Run the following command:
# snakemake -np -s Snakefile \
# --cluster "qsub -q all.q -V -cwd -l h={params.server} -pe smp {threads} -l mem_free={resources.mem_free}G" \
# --jobs 10
#
# Build DAG:
#
# snakemake -s Snakefile --dag | dot -Tsvg > dag.svg
import pandas as pd
import yaml

sample_lists = pd.read_csv("accessions/sugarcane_accessions.txt", sep="\t")
genotypes = sample_lists["genotype"]
samples = sample_lists["sample"]

rule all:
	input:
#		[f"../data/{genotypes[i]}/1_raw_reads_in_fastq_format/{samples[i]}_1.fastq.gz" for i in range(len(genotypes))],
#		[f"../data/{genotypes[i]}/1_raw_reads_in_fastq_format/{samples[i]}_2.fastq.gz" for i in range(len(genotypes))],
#               [f"../data/{genotypes[i]}/2_trimmed/{samples[i]}_trimmed_1.fastq.gz" for i in range(len(genotypes))],
#               [f"../data/{genotypes[i]}/2_trimmed/{samples[i]}_trimmed_2.fastq.gz" for i in range(len(genotypes))],
#		"/home/dmpachon/jorge/comparative_cane/data/scpan_longest_salmon_index",
               "/home/dmpachon/jorge/comparative_cane/data/scpan_all_transcripts_salmon_index",
		[f"../data/{genotypes[i]}/3_salmon_quant_all/{samples[i]}" for i in range(len(genotypes))]


rule sra_fetch:
        output:
                SRA = "../data/{genotype}/1_raw_reads_in_fastq_format/{sample}/{sample}.sra"
        threads: 1
        resources:
                load=60
        params:
                genotype="{genotype}",
                partition="long",
                jobname="fetch",
                mem="20gb"
        conda:
                "conda_envs/sra-toolkit.yml"
        log:
                "logs/{genotype}_{sample}_fetch.log"
        shell:
                """ 
                cd ../data/{wildcards.genotype}/1_raw_reads_in_fastq_format && \
                prefetch {wildcards.sample}
                """

rule sra_split:
        input:
                SRA = "../data/{genotype}/1_raw_reads_in_fastq_format/{sample}/{sample}.sra"
	output:
		R1 = "../data/{genotype}/1_raw_reads_in_fastq_format/{sample}_1.fastq.gz",
		R2 = "../data/{genotype}/1_raw_reads_in_fastq_format/{sample}_2.fastq.gz"
	threads: 1
	resources:
		load=1
	params:
		genotype="{genotype}",
		partition="long",
		jobname="split",
                mem="8gb"
	conda:
		"conda_envs/sra-toolkit.yml"
	log:
		"logs/{genotype}_{sample}_split.log"
	shell:
		"""
		cd ../data/{wildcards.genotype}/1_raw_reads_in_fastq_format && \
                fastq-dump --split-files --gzip {wildcards.sample} && \
                rm -rf {wildcards.sample}
		"""
rule bbduk:
        input:
                R1="../data/{genotype}/1_raw_reads_in_fastq_format/{sample}_1.fastq.gz",
                R2="../data/{genotype}/1_raw_reads_in_fastq_format/{sample}_2.fastq.gz"
        output:
                R1 = "../data/{genotype}/2_trimmed/{sample}_trimmed_1.fastq.gz",
                R2 = "../data/{genotype}/2_trimmed/{sample}_trimmed_2.fastq.gz",
                refstats = "../data/{genotype}/2_trimmed/{sample}_trimmed.refstats",
                stats = "../data/{genotype}/2_trimmed/{sample}_trimmed.stats"
        threads: 10
        resources:
                load=10
        params:
                genotype="{genotype}",
                partition="long",
                jobname="bbduk",
                mem="20gb"
        conda:
                "conda_envs/bbduk.yml"
        log:
                "logs/{genotype}_{sample}_bbmap.log"
        shell:
                "bbduk.sh -Xmx20g threads={threads} in1={input.R1} in2={input.R2} "
                "refstats={output.refstats} stats={output.stats} "
                "out1={output.R1} out2={output.R2} "
                "ref=/home/dmpachon/jorge/comparative_cane/data/dbs/adapters.fa,"
                "/home/dmpachon/jorge/comparative_cane/data/dbs/rfam-5.8s-database-id98.fasta,"
                "/home/dmpachon/jorge/comparative_cane/data/dbs/silva-bac-16s-id90.fasta,"
                "/home/dmpachon/jorge/comparative_cane/data/dbs/rfam-5s-database-id98.fasta,"
                "/home/dmpachon/jorge/comparative_cane/data/dbs/silva-bac-23s-id98.fasta,"
                "/home/dmpachon/jorge/comparative_cane/data/dbs/silva-arc-16s-id95.fasta,"
                "/home/dmpachon/jorge/comparative_cane/data/dbs/silva-euk-18s-id95.fasta,"
                "/home/dmpachon/jorge/comparative_cane/data/dbs/silva-arc-23s-id98.fasta,"
                "/home/dmpachon/jorge/comparative_cane/data/dbs/silva-euk-28s-id98.fasta "
                "minlength=50 qtrim=w trimq=20 tpe tbo 2> {log}"
rule salmon_index:
        input:
                transcriptome="/home/dmpachon/jorge/comparative_cane/data/all_transcripts_checked_with_cds.fasta"
        output:
                salmon_index=directory("/home/dmpachon/jorge/comparative_cane/data/sugarcane_all_transcripts_salmon_index")
        params:
                jobname="index",
                mem="15gb",
                partition="long"
        resources:
                load=10
        conda:
                "conda_envs/salmon.yml"
        threads: 30
        log:
                "logs/index_salmon.log"
        shell:
                "salmon index -t {input.transcriptome} -p {threads} -i {output.salmon_index} 2> {log}"
rule salmon_quant:
        input:
                R1="../data/{genotype}/2_trimmed/{sample}_trimmed_1.fastq.gz",
                R2="../data/{genotype}/2_trimmed/{sample}_trimmed_2.fastq.gz",
                index="/home/dmpachon/jorge/comparative_cane/data/sugarcane_all_transcripts_salmon_index"
        output:
                salmon_out=directory("../data/{genotype}/3_salmon_quant_all/{sample}")
        params:
                jobname="quant",
                mem="100gb",
                partition="long"
        resources:
                load=35
        conda:
                "conda_envs/salmon.yml"
        threads: 25
        log:
                "logs/{genotype}_{sample}_salmon.log"
        shell:
                "salmon quant -i {input.index} -l A -1 {input.R1} -2 {input.R2} -p {threads} -o {output.salmon_out} >> {log} 2>&1"

#rule quant_results:
#        input:
#            salmon_quant="../data/{genotype}/3_salmon_quant/{sample}"
#        output:
#            pretty_table="../results/pretty_table.csv"
#        params:
#                jobname=""
