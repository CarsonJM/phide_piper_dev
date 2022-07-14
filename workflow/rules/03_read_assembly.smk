# -------------------------------------
# Read assembly (Only runs if data_type: "reads")
# -------------------------------------
import pandas as pd


# Load sample information and validate
configfile: "config/config.yaml"
samples_df = pd.read_csv("config/samples.tsv", sep="\t")
samples = samples_df['sample']
coassemblies = list(set(samples_df['assembly']))


# load results path
results = config["results"]


# load resources path
resources = config["resources"]


# load report
report: "report/workflow.rst"


# -------------------------------------
# Read assembly rules
# -------------------------------------
# -----------------------------------------------------
# 03 metaSPAdes
# -----------------------------------------------------
# Assemble reads using metaspades (single)
checkpoint metaspades_single:
    input:
        R1=results
        + "01_READ_PREPROCESSING/04_kneaddata/{assembly}_paired_1.fastq",
        R2=results
        + "01_READ_PREPROCESSING/04_kneaddata/{assembly}_paired_2.fastq",
    output:
        results
        + "03_READ_ASSEMBLY/03_metaspades/{assembly}/"
        + config["read_assembly"]["assembly_output"]
        + ".fasta",
    params:
        output_dir=results + "03_READ_ASSEMBLY/03_metaspades/{assembly}",
        extra_args=config["read_assembly"]["metaspades_arguments"],
    log:
        results + "00_LOGS/03_read_assembly_{assembly}.metaspades_single.log",
    conda:
        "../envs/metaspades.yml"
    threads: config["read_assembly"]["metaspades_threads"]
    shell:
        """
        # assemble reads using metaspades
        spades.py \
        --meta \
        -1 {input.R1} \
        -2 {input.R2} \
        -o {params.output_dir} \
        --threads {threads} \
        {params.extra_args}

        # copy spades.log to log file
        cp {params.output_dir}/spades.log {log}
        """

# Assemble reads using metaspades (coassembly)
checkpoint metaspades_coassembly:
    input:
        R1=results + "03_READ_ASSEMBLY/01_combine_coassembly_reads/{assembly}_coassembly_1.fastq",
        R2=results + "03_READ_ASSEMBLY/01_combine_coassembly_reads/{assembly}_coassembly_2.fastq",
    output:
        results
        + "03_READ_ASSEMBLY/03_metaspades/{assembly}_coassembly/"
        + config["read_assembly"]["assembly_output"]
        + ".fasta",
    params:
        output_dir=results + "03_READ_ASSEMBLY/03_metaspades/{assembly}_coassembly",
        extra_args=config["read_assembly"]["metaspades_arguments"],
    log:
        results + "00_LOGS/03_read_assembly_{assembly}.metaspades_coassembly.log",
    conda:
        "../envs/metaspades.yml"
    threads: config["read_assembly"]["metaspades_threads"]
    shell:
        """
        # assemble reads using metaspades
        spades.py \
        --meta \
        -1 {input.R1} \
        -2 {input.R2} \
        -o {params.output_dir} \
        --threads {threads} \
        {params.extra_args}

        # copy spades.log to log file
        cp {params.output_dir}/spades.log {log}
        """

# -----------------------------------------------------
# 04 QUAST
# -----------------------------------------------------
# run quast to determine the quality of single assemblies
rule quast_single:
    input:
        results
        + "03_READ_ASSEMBLY/03_metaspades/{sample}/"
        + config["read_assembly"]["assembly_output"]
        + ".fasta",
    output:
        results + "03_READ_ASSEMBLY/04_quast/{sample}/transposed_report.tsv",
    params:
        output_dir=results + "03_READ_ASSEMBLY/04_quast/{sample}",
        min_len=config["read_assembly"]["min_contig_length"],
        labels="{sample}",
        extra_args=config["read_assembly"]["quast_arguments"],
    log:
        results + "00_LOGS/03_read_assembly_{sample}.quast_single.log",
    conda:
        "../envs/quast.yml"
    shell:
        """
        # assembly analysis using quast
        metaquast.py \
        {input} \
        -o {params.output_dir} \
        --threads {threads} \
        --min-contig {params.min_len} \
        --contig-thresholds 0,1000,5000,10000,{params.min_len} \
        --labels {params.labels} \
        {params.extra_args}

        # copy spades.log to log file
        cp {params.output_dir}/quast.log {log}
        """
    
# combine quast outputs
rule combine_quast_single_across_samples:
    input:
        expand(results + "03_READ_ASSEMBLY/04_quast/{sample}/transposed_report.tsv", sample=samples),
    output:
        results + "03_READ_ASSEMBLY/read_assembly_report_single.tsv",
    shell:
        """
        # combine quast reports for all assemblies, only keeping the header from one file
        awk 'FNR>1 || NR==1' {input} > {output}
        """

# run quast to determine the quality of coassemblies
rule quast_coassembly:
    input:
        results
        + "03_READ_ASSEMBLY/03_metaspades/{assembly}_coassembly/"
        + config["read_assembly"]["assembly_output"]
        + ".fasta",
    output:
        results + "03_READ_ASSEMBLY/04_quast/{assembly}_coassembly/transposed_report.tsv",
    params:
        output_dir=results + "03_READ_ASSEMBLY/04_quast/{assembly}_coassembly",
        min_len=config["read_assembly"]["min_contig_length"],
        labels="{assembly}_coassembly",
        extra_args=config["read_assembly"]["quast_arguments"],
    log:
        results + "00_LOGS/03_read_assembly_{assembly}.quast_coassembly.log",
    conda:
        "../envs/quast.yml"
    shell:
        """
        # assembly analysis using quast
        metaquast.py \
        {input} \
        -o {params.output_dir} \
        --threads {threads} \
        --min-contig {params.min_len} \
        --contig-thresholds 0,1000,5000,10000,{params.min_len} \
        --labels {params.labels} \
        {params.extra_args}

        # copy spades.log to log file
        cp {params.output_dir}/quast.log {log}
        """

# combine quast outputs
rule combine_quast_coassembly_across_assemblies:
    input:
        expand(results + "03_READ_ASSEMBLY/04_quast/{assembly}_coassembly/transposed_report.tsv", assembly=coassemblies),
    output:
        results + "03_READ_ASSEMBLY/read_assembly_report_coassembly.tsv",
    shell:
        """
        # combine quast reports for all assemblies, only keeping the header from one file
        awk 'FNR>1 || NR==1' {input} > {output}
        """

# # -----------------------------------------------------
# # 05 Contig length filter (single)
# # -----------------------------------------------------
# # filter contigs based on contig length
# rule contig_length_filter_single:
#     input:
#         results
#         + "03_READ_ASSEMBLY/03_metaspades/{sample}/"
#         + config["read_assembly"]["assembly_output"]
#         + ".fasta",
#     output:
#         results
#         + "03_READ_ASSEMBLY/05_contig_length_filter/{sample}_"
#         + config["read_assembly"]["assembly_output"]
#         + ".fasta",
#     params:
#         min_length=config["read_assembly"]["min_contig_length"],
#     conda:
#         "../envs/jupyter.yml"
#     notebook:
#         "../notebooks/03_contig_length_filter.py.ipynb"


# # filter contigs based on contig length
# rule contig_length_filter_coassembly:
#     input:
#         results
#         + "03_READ_ASSEMBLY/03_metaspades/{assembly}_coassembly/"
#         + config["read_assembly"]["assembly_output"]
#         + ".fasta",
#     output:
#         results
#         + "03_READ_ASSEMBLY/05_contig_length_filter/{assembly}_coassembly_"
#         + config["read_assembly"]["assembly_output"]
#         + ".fasta",
#     params:
#         min_length=config["read_assembly"]["min_contig_length"],
#     conda:
#         "../envs/jupyter.yml"
#     notebook:
#         "../notebooks/03_contig_length_filter.py.ipynb"


# # -----------------------------------------------------
# # Analyze assemblies
# # -----------------------------------------------------
# # analyze quast results to visualize assembly quality
# rule read_assembly_analysis:
#     input:
#         single=results + "03_READ_ASSEMBLY/read_assembly_report_single.tsv",
#         coassembly=results + "03_READ_ASSEMBLY/read_assembly_report_coassembly.tsv",
#     output:
#         report(
#             results + "03_READ_ASSEMBLY/read_assembly_figure.png",
#             caption="../report/read_assembly_analysis_contig_count.rst",
#             category="Step 02: Read assembly",
#         ),
#     params:
#         min_len=config["read_assembly"]["min_contig_length"],
#     conda:
#         "../envs/jupyter.yml"
#     notebook:
#         "../notebooks/03_read_assembly_analysis.py.ipynb"
