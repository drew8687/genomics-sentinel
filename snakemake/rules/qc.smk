# ──────────────────────────────────────────────────────────
#  qc.smk — FastQC + Trimmomatic + MultiQC
# ──────────────────────────────────────────────────────────
rule fastqc_raw:
    input:
        r1 = "data/raw/{sample}_R1.fastq.gz",
        r2 = "data/raw/{sample}_R2.fastq.gz"
    output:
        html_r1 = "results/qc/fastqc/{sample}_R1_fastqc.html",
        html_r2 = "results/qc/fastqc/{sample}_R2_fastqc.html",
    threads: config["threads"]["fastqc"]
    shell:
        "fastqc -t {threads} -o results/qc/fastqc/ {input.r1} {input.r2}"

rule trimmomatic:
    input:
        r1 = "data/raw/{sample}_R1.fastq.gz",
        r2 = "data/raw/{sample}_R2.fastq.gz"
    output:
        r1p = "data/trimmed/{sample}_R1_paired.fastq.gz",
        r2p = "data/trimmed/{sample}_R2_paired.fastq.gz",
        r1u = "data/trimmed/{sample}_R1_unpaired.fastq.gz",
        r2u = "data/trimmed/{sample}_R2_unpaired.fastq.gz",
        log = "logs/trimmomatic/{sample}.log"
    threads: config["threads"]["trimmomatic"]
    shell:
        """
        trimmomatic PE -threads {threads} -phred33 \
            {input.r1} {input.r2} \
            {output.r1p} {output.r1u} \
            {output.r2p} {output.r2u} \
            ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 \
            2> {output.log}
        """

rule multiqc:
    input:
        expand("results/qc/fastqc/{sample}_R1_fastqc.html", sample=config["samples"]),
        expand("logs/trimmomatic/{sample}.log", sample=config["samples"])
    output:
        "results/qc/multiqc/multiqc_report.html"
    shell:
        "multiqc results/qc/ logs/ -o results/qc/multiqc/ --force"
