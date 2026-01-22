# chrX Polyploid Variant Calling Pipeline

Snakemake pipeline for calling SNPs on chromosome X from WGS data in diploid or polyploid samples.

## Features
- Fastp QC, trimming, deduplication
- BWA-MEM alignment
- chrX-only BAM filtering (hardcoded)
- GATK HaplotypeCaller with configurable ploidy
- Optional diploidisation for polyploid samples

## Requirements
- Conda / Mamba
- Snakemake ≥ 7
- Indexed reference genome (bwa, samtools, GATK)

## Input
Paired-end FASTQs:
sample_1.fq.gz
sammple_2.fq.gz

Place in:
data/raw_fastq


## Configuration
Edit `config.yaml`:
```yaml
genome: "/path/to/genome.fa"
ploid: 3
threads: 8
```
## Run
Dry run:
snakemake -n --use-conda

Run:
snakemake --use-conda --cores 8


## Notes
- chrX is intentionally hardcoded
- Diploidisation only runs if ploid > 2
