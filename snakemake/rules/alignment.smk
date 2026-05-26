# ──────────────────────────────────────────────────────────
#  alignment.smk — STAR splice-aware alignment
# ──────────────────────────────────────────────────────────
rule star_align:
    input:
        r1  = "data/trimmed/{sample}_R1_paired.fastq.gz",
        r2  = "data/trimmed/{sample}_R2_paired.fastq.gz",
        idx = config["genome"]["star_index"]
    output:
        bam = temp("results/alignment/{sample}.Aligned.sortedByCoord.out.bam"),
        log = "logs/star/{sample}.log"
    threads: config["threads"]["star"]
    params:
        prefix = "results/alignment/{sample}."
    shell:
        """
        STAR --runThreadN {threads} \
            --genomeDir {input.idx} \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI NM MD AS \
            --outFileNamePrefix {params.prefix} \
            --outFilterMultimapNmax 20 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --sjdbScore 1 \
            > {output.log} 2>&1
        """

rule picard_mark_duplicates:
    input:  "results/alignment/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        bam     = "results/alignment/{sample}.sorted.bam",
        metrics = "logs/picard/{sample}.dup_metrics.txt"
    shell:
        """
        picard MarkDuplicates \
            I={input} O={output.bam} \
            M={output.metrics} \
            REMOVE_DUPLICATES=false \
            CREATE_INDEX=true
        """

rule samtools_index:
    input:  "results/alignment/{sample}.sorted.bam"
    output: "results/alignment/{sample}.sorted.bam.bai"
    threads: config["threads"]["samtools"]
    shell:  "samtools index -@ {threads} {input}"

rule samtools_flagstat:
    input:  "results/alignment/{sample}.sorted.bam"
    output: "logs/flagstat/{sample}.flagstat"
    shell:  "samtools flagstat {input} > {output}"


# ──────────────────────────────────────────────────────────
#  variants.smk — GATK HaplotypeCaller → VQSR → SnpEff
# ──────────────────────────────────────────────────────────
rule gatk_haplotype_caller:
    input:
        bam   = "results/alignment/{sample}.sorted.bam",
        ref   = config["genome"]["fasta"],
        dbsnp = config["variants"]["dbsnp"]
    output:
        gvcf = "results/variants/gvcf/{sample}.g.vcf.gz"
    threads: config["threads"]["gatk"]
    shell:
        """
        gatk HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.gvcf} \
            -D {input.dbsnp} \
            --emit-ref-confidence GVCF \
            --native-pair-hmm-threads {threads}
        """

rule gatk_genotype_gvcfs:
    input:
        gvcfs = expand("results/variants/gvcf/{sample}.g.vcf.gz", sample=config["samples"]),
        ref   = config["genome"]["fasta"]
    output:
        vcf = "results/variants/raw.vcf.gz"
    params:
        gvcf_args = lambda w, input: " ".join([f"-V {v}" for v in input.gvcfs])
    shell:
        """
        gatk CombineGVCFs -R {input.ref} {params.gvcf_args} -O combined.g.vcf.gz
        gatk GenotypeGVCFs -R {input.ref} -V combined.g.vcf.gz -O {output.vcf}
        """

rule gatk_variant_filter:
    input:
        vcf = "results/variants/raw.vcf.gz",
        ref = config["genome"]["fasta"]
    output:
        vcf = "results/variants/filtered.vcf.gz"
    shell:
        """
        gatk VariantFiltration \
            -R {input.ref} -V {input.vcf} \
            --filter-expression "QD < 2.0" --filter-name "QD2" \
            --filter-expression "FS > 60.0" --filter-name "FS60" \
            --filter-expression "MQ < 40.0" --filter-name "MQ40" \
            -O {output.vcf}
        """

rule snpeff_annotate:
    input:  "results/variants/filtered.vcf.gz"
    output: "results/variants/annotation/annotated.vcf"
    shell:
        """
        java -jar /opt/snpEff/snpEff.jar \
            GRCh38.105 {input} \
            -stats results/variants/annotation/snpeff_stats.html \
            > {output}
        """
