configfile: "config.yaml"

SAMPLES, = glob_wildcards("data/raw_fastq/{sample}_1.fq.gz")

if not SAMPLES:
	raise ValueError("No FASTQ files found in data/raw_fastq matching *_1.fq.gz")

CHR = "chrX"

# Decide final outputs based on ploidy
if config["ploid"] > 2:
    FINAL_VCFS = expand("results/vcf/{sample}_chrX_diploid.vcf.gz", sample=SAMPLES)
else:
    FINAL_VCFS = expand("results/vcf/{sample}_chrX_with_reference.vcf.gz", sample=SAMPLES)

rule all:
    input:
        expand(
            "results/vcf/{sample}_chrX_diploid.vcf.gz",
            sample=SAMPLES
        ) if config["ploid"] > 2 else
        expand(
            "results/vcf/{sample}_chrX_with_reference.vcf.gz",
            sample=SAMPLES
        )

###############################################################################
# Step 1: fastp
###############################################################################

rule fastp:
    input:
        r1="data/raw_fastq/{sample}_1.fq.gz",
        r2="data/raw_fastq/{sample}_2.fq.gz"
    output:
        r1="results/fastp/{sample}_Trimmed_DeDuped_1.fq.gz",
        r2="results/fastp/{sample}_Trimmed_DeDuped_2.fq.gz",
        html="results/fastp/{sample}.html",
        json="results/fastp/{sample}.json"
    threads: config["threads"]
    conda: "envs/fastp.yaml"
    log:
        "logs/fastp/{sample}.log"
    shell:
        """
        fastp -D \
              -w {threads} \
              -i {input.r1} \
              -I {input.r2} \
              -o {output.r1} \
              -O {output.r2} \
              -h {output.html} \
              -j {output.json} \
              > {log} 2>&1
        """

###############################################################################
# Step 2: BWA MEM
###############################################################################

rule bwa_mem:
    input:
        r1="results/fastp/{sample}_Trimmed_DeDuped_1.fq.gz",
        r2="results/fastp/{sample}_Trimmed_DeDuped_2.fq.gz",
        ref=config["genome"]
    output:
        bam="results/bam/{sample}_BWA.sorted.bam"
    threads: config["threads"]
    conda: "envs/alignment.yaml"
    log:
        "logs/bwa/{sample}.log"
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.r1} {input.r2} | \
        samtools sort -@ {threads} -o {output.bam} \
        > {log} 2>&1
        """

###############################################################################
# Step 3: Extract chrX
###############################################################################

rule extract_chrX:
    input:
        bam="results/bam/{sample}_BWA.sorted.bam"
    output:
        bam="results/bam/{sample}_chrX.bam"
    conda: "envs/alignment.yaml"
    log:
        "logs/chrX/{sample}.log"
    shell:
        """
        samtools view -h {input.bam} {CHR} | \
        samtools sort -o {output.bam} \
        > {log} 2>&1
        """

###############################################################################
# Step 4: Add read groups
###############################################################################

rule add_read_groups:
    input:
        bam="results/bam/{sample}_chrX.bam"
    output:
        bam="results/bam/{sample}_chrX_RG.bam"
    conda: "envs/gatk.yaml"
    log:
        "logs/rg/{sample}.log"
    shell:
        """
        gatk AddOrReplaceReadGroups \
            -I {input.bam} \
            -O {output.bam} \
            -RGID 1 \
            -RGLB lib1 \
            -RGPL illumina \
            -RGPU unit1 \
            -RGSM {wildcards.sample} \
            > {log} 2>&1
        """

###############################################################################
# Step 5: Mark duplicates
###############################################################################

rule mark_duplicates:
    input:
        bam="results/bam/{sample}_chrX_RG.bam"
    output:
        bam="results/bam/{sample}_chrX_RG_dedup.bam",
        metrics="results/bam/{sample}_chrX_RG_dedup.metrics.txt"
    conda: "envs/gatk.yaml"
    log:
        "logs/dedup/{sample}.log"
    shell:
        """
        gatk MarkDuplicates \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} \
            > {log} 2>&1
        """

###############################################################################
# Step 6: HaplotypeCaller (chrX only)
###############################################################################

rule haplotypecaller:
    input:
        bam="results/bam/{sample}_chrX_RG_dedup.bam",
        ref=config["genome"]
    output:
        vcf="results/vcf/{sample}_chrX.vcf.gz"
    conda: "envs/gatk.yaml"
    log:
        "logs/gatk/{sample}.log"
    shell:
        """
        gatk --java-options "-Xmx{config[java_mem]}" HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.vcf} \
            --sample-ploidy {config[ploid]} \
            -L {CHR} \
            > {log} 2>&1
        """

###############################################################################
# Step 7: Add reference header
###############################################################################

rule add_reference:
    input:
        vcf="results/vcf/{sample}_chrX.vcf.gz"
    output:
        vcf="results/vcf/{sample}_chrX_with_reference.vcf.gz"
    conda: "envs/vcf.yaml"
    shell:
        """
        zcat {input.vcf} | \
        sed "/^#CHROM/i ##reference={config[reference_name]}" | \
        bgzip > {output.vcf}
        tabix -p vcf {output.vcf}
        """

###############################################################################
# Step 8: Diploidisation (polyploid chrX only)
###############################################################################

rule diploidise:
    input:
        vcf="results/vcf/{sample}_chrX_with_reference.vcf.gz"
    output:
        vcf="results/vcf/{sample}_chrX_diploid.vcf.gz"
    conda: "envs/vcf.yaml"
    shell:
        r"""
        bcftools view -h {input.vcf} > tmp.vcf
        bcftools view -H {input.vcf} | \
        awk -F'\t' '{{
            split($10,a,":");
            if (a[1] ~ /1/ && a[1] ~ /0/) a[1]="0/1";
            else if (a[1]=="1/1/1") a[1]="1/1";
            else a[1]="0/0";
            $10=a[1];
            print
        }}' OFS='\t' >> tmp.vcf
        bgzip -c tmp.vcf > {output.vcf}
        tabix -p vcf {output.vcf}
        rm tmp.vcf
        """

