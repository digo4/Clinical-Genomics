[![Static Badge](https://img.shields.io/badge/LICENSE-MIT-yellow)](https://opensource.org/license/mit/) [![Static Badge](https://img.shields.io/badge/FASTQC-v.0.12.0-blue)](https://github.com/s-andrews/FastQC) [![Static Badge](https://img.shields.io/badge/fastp-v.0.20.1-blue)](https://github.com/OpenGene/fastp/releases)  [![Static Badge](https://img.shields.io/badge/bwa-v.0.7.17-blue)](https://github.com/lh3/bwa)  [![Static Badge](https://img.shields.io/badge/java-%3E%3Dv.17-pink)
](https://openjdk.org/projects/jdk/17/) [![Static Badge](https://img.shields.io/badge/gatk-%3E%3Dv.4.4.0.0-teal)
](https://hub.docker.com/r/broadinstitute/gatk/) [![Static Badge](https://img.shields.io/badge/snpEff%26SnpSift-%3Dv.5.2-purple)
](https://pcingola.github.io/SnpEff/)    [![Static Badge](https://img.shields.io/badge/InterVar-%3E%3Dv.2.2.2-green)](https://github.com/WGLab/InterVar)


## Introduction
Clinical genomics analysis refers to the application of genomic technologies and computational methods to understand the genetic basis of diseases and inform clinical decision-making. It involves the analysis of genomic data obtained from patients to identify genetic variations, mutations, and other molecular features that may contribute to disease development or progression. The goal is to use this information to guide diagnosis, treatment decisions, and personalized medicine approaches.

## A step-by-step guide to Clinical Genomics analysis:
### 1. **Pre-processing reads:** 
Pre-processing the reads involves checking the quality (Phred-scores) using tools like [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) or [MultiQC](https://github.com/MultiQC/MultiQC). After that, we need to remove adapter contamination and trim quality reads using either one of the following tools : [fastp](https://github.com/OpenGene/fastp?tab=readme-ov-file), [cutadapt](https://cutadapt.readthedocs.io/en/stable/), [TrimGalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) or [Trimmomatic](/http://www.usadellab.org/cms/?page=trimmomatic).

The following code can be used to check the quality of reads obtained from the sequencer.

#### Running FastQC

```
mkdir fastq_bt  # bt stands for before trimming
#Running FASTQC
fastqc -o fastq_bt gatk_demo1.fastq.gz gatk_demo2.fastq.gz

```
#### Running FASTP
```
fastp -i gatk_demo1.fastq.gz -I gatk_demo2.fastq.gz -o gatk_demo1_trimmed.fastq.gz -O gatk_demo2_trimmed.fastq.gz \
       -R gatk_demo -h gatk_demo.html -j gatk_demo.json --detect_adapter_for_pe

In case you want to use cutadapt or trimgalore ( which is essentially a wrapper around cutadapt, 
```
#### Running FastQC again after Adapter and Quality Trimming
```
mkdir fastq_at  # at stands for after trimming
#Running FASTQC
fastqc -o fastq_bt gatk_demo1.fastq.gz gatk_demo2.fastq.gz

```
### 2. Aligning reads to reference genome:
Aligning sequencing reads to a reference genome is a crucial step in many bioinformatics workflows, including variant calling and downstream genomic analyses. The Burrows-Wheeler Aligner (BWA) is a popular tool for aligning short DNA sequences to a large reference genome efficiently.

  
### 3. Alignment Post-Processing:
5. **Base Quality score recalibration (not mandatory, but highly recommended)**
6. **Variant Calling:**
7. **Variant Filtering:**
   - 6.1 **Splitting variants into SNPs and INDELs-**
   - 6.2 **Variant Quality Score Recalibration-**
   - 6.3 **Hard Filtering Variants-**
8. **Variant Annotation:**
9. **Clinical Interpretation:**
   



