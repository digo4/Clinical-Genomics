[![Static Badge](https://img.shields.io/badge/LICENSE-MIT-yellow)](https://opensource.org/license/mit/) [![Static Badge](https://img.shields.io/badge/FASTQC-v.0.12.0-blue)](https://github.com/s-andrews/FastQC) [![Static Badge](https://img.shields.io/badge/fastp-v.0.20.1-blue)](https://github.com/OpenGene/fastp/releases)  [![Static Badge](https://img.shields.io/badge/bwa-v.0.7.17-blue)](https://github.com/lh3/bwa)  [![Static Badge](https://img.shields.io/badge/java-%3E%3Dv.17-pink)
](https://openjdk.org/projects/jdk/17/) [![Static Badge](https://img.shields.io/badge/gatk-%3E%3Dv.4.4.0.0-teal)
](https://hub.docker.com/r/broadinstitute/gatk/) [![Static Badge](https://img.shields.io/badge/snpEff%26SnpSift-%3Dv.5.2-purple)
](https://pcingola.github.io/SnpEff/)    [![Static Badge](https://img.shields.io/badge/InterVar-%3E%3Dv.2.2.2-green)](https://github.com/WGLab/InterVar)


## Introduction
Clinical genomics analysis refers to the application of genomic technologies and computational methods to understand the genetic basis of diseases and inform clinical decision-making. It involves the analysis of genomic data obtained from patients to identify genetic variations, mutations, and other molecular features that may contribute to disease development or progression. The goal is to use this information to guide diagnosis, treatment decisions, and personalized medicine approaches.

## Prerequisites:
Before we begin, please have all the files from the [GATK Resource Bundle](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false) and keep all the files downloaded in one folder (let's name the folder resource). We can create a folder using the `mkdir` command.


## A step-by-step guide to Clinical Genomics analysis:

### 1. Pre-processing reads: 
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
fastqc -o fastq_at gatk_demo1_trimmed.fastq.gz gatk_demo2_trimmed.fastq.gz

```
### 2. Aligning reads to reference genome:
Aligning sequencing reads to a reference genome is a crucial step in many bioinformatics workflows, including variant calling and downstream genomic analyses. The Burrows-Wheeler Aligner [BWA](https://github.com/lh3/bwa) is a popular tool for aligning short DNA sequences to a large reference genome efficiently. Before performing the alignment, we first need to index the reference genome.
#### Indexing the reference genome
```
cd resource
bwa index Homo_sapiens_assembly38.fasta 
```
Note: Since the Homo spaiens reference genome is very huge in size (3GB), the indexing process can take anywhere between 45min to 2hrs, depending on your system configurations. Alternatively, if you have downloaded all the files in the GATK resource bundle folder, then the respective index files from BWA are already downloaded, and you do not need to index your reference genome again.

#### Alignment to the indexed reference genome
```
bwa mem Homo_sapiens_assembly38.fasta gatk_demo1_trimmed.fastq.gz gatk_demo2_trimmed.fastq.gz > demo.sam

```

### 3. Alignment Post-Processing:
After aligning your sequencing reads to a reference genome using BWA, there are several post-processing steps you may want to perform to analyze and manipulate the alignment results. For this purpose, we will be utilising the [Picard](https://broadinstitute.github.io/picard/) tool. Picard is a collection of command-line tools for manipulating high-throughput sequencing (HTS) data files, such as SAM and BAM files. These tools are developed by the Broad Institute and are widely used in genomics and bioinformatics workflows. Picard tools are often employed for quality control, data preprocessing, and various manipulations of sequence data. Here are some common post-processing steps:
- **3.1 Adding one or more read groups to your SAM file:** The [AddOrReplaceReadGroups](https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-) tool in Picard is used to add or replace read group information in a SAM or BAM file. This is important for downstream applications that require proper grouping and identification of reads, especially when dealing with data from multiple sequencing libraries or samples. Even though for this demo data purpose, we have a single sample bam file, still this step is absolutely mandatory for the GATK pipeline! The GATK Pipeline requires at least one ReadGroup id to be specified in the  Here's an example of how to use the AddOrReplaceReadGroups tool:
  ```
  java -jar picard.jar AddOrReplaceReadGroups \
         I=demo.sam \
         O=demo_rg.sam \
         RGID=ID123 \
         RGLB=library1 \
         RGPL=illumina \
         RGPU=unit1 \
         RGSM=sample1

  ```

  Note: For those of you who are running [GATK in Docker](https://gatk.broadinstitute.org/hc/en-us/articles/360035889991--How-to-Run-GATK-in-a-Docker-container), you dont need to install picard separately. All the functionalities of picard are already present within GATK, check out the script (clinical_genomics_demo.sh) for clarity. However, those of you who will have downloaded the GATK https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip. Make sure you have Java (v.>= 17) installed, and then you can execute the Picard command as shown in the example. Additionally, specify the full path to the picard.jar file in the command, depending on your system configuratio
- **3.2 SAM to BAM conversion:** The [SamFormatConverter](https://gatk.broadinstitute.org/hc/en-us/articles/360037058992-SamFormatConverter-Picard-) tool in Picard can be used to convert a SAM file to BAM format. Here's an example of how to use SamFormatConverter for this purpose:
  ```
  java -jar picard.jar SamFormatConverter \
         I=demo_rg.sam\
         O=demo_rg.bam

  ```
- **3.3 Sorting BAM file:** The [SortSam](https://gatk.broadinstitute.org/hc/en-us/articles/360037594291-SortSam-Picard-) tool in Picard is used to sort a SAM or BAM file by coordinate order. Sorting is a necessary step for downstream analyses and visualization tools that require data to be ordered by genomic coordinates. Here's an example of how to use SortSam:
  ```
  java -jar picard.jar SortSam \
         I=demo_rg.bam \
         O=demo_sorted.bam \
         SORT_ORDER=coordinate

  ```
- **3.4 Marking (and optionally deleting) duplicates:** In Picard, the [MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-) tool is used to identify and mark duplicate reads in a SAM or BAM file. Marking duplicates is a common step in the preprocessing of sequencing data, especially when dealing with PCR-amplified libraries. Here's an example of how to use MarkDuplicates:
  ```
  java -jar picard.jar MarkDuplicates \
      I=demo_sorted.bam \
      O=demo_sorted_md.bam\
      M=marked_dup_metrics.txt
  ```
In case we want to remove the duplicated reads the command goes as below:
```
  java -jar picard.jar MarkDuplicates \
      I=demo_sorted.bam \
      O=demo_sorted_dedup.bam\
      M=marked_dup_metrics.txt
       --REMOVE_DUPLICATES
  ```
### 4. Base Quality score recalibration (not mandatory, but highly recommended) :
Base Quality Score Recalibration (BQSR) is primarily used to improve the accuracy of base quality scores assigned to individual nucleotides in the sequencing data.
In NGS, the quality scores associated with each base in a sequence represent the estimated probability of an incorrect base call. However, these scores are not always perfectly accurate, and recalibration is performed to correct for systematic errors and biases in the sequencing data. The general process of Base Quality Score Recalibration involves the following steps:

- **Collecting Metrics:** The first step involves collecting various metrics from the sequencing data, including the actual base calls, base quality scores, and additional information such as machine cycle, base context, and sequence context.

- **Model Training:** A statistical model is then trained using the collected metrics to estimate the true error probability for each base call. The model takes into account various factors that can influence the accuracy of base calls, such as the quality scores of neighboring bases, machine-specific biases, and sequence context.

- **Recalibration:** Based on the trained model, the base quality scores are recalibrated to reflect a more accurate estimate of the true error probability. This recalibration involves adjusting the original quality scores assigned during the sequencing process.

- **Quality Control:** After recalibration, quality control metrics are often generated to assess the success of the recalibration process. These metrics help ensure that the recalibrated data meets certain quality standards.
In the GATK pipeline, we can perform BQSR using two commands : [BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator) and [ApplyBQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360037055712-ApplyBQSR). BaseRecalibrator builds a machine learning model on the basis of known variant sites and generates a recalibration table on the statistical accuracy of sequencer-predicted base quality scores.
```
gatk BaseRecalibrator \
       -I demo_sorted_dedup.bam \
       -R resource/Homo_sapiens_assembly38.fasta \
       --known-sites resource/Homo_sapiens_assembly38.known_indels.vcf.gz \
       --known-sites resource/Homo_sapiens_assembly38.dbsnp.vcf.gz \
       -O recal_data.table
```
ApplyBQSR then applies the recalibration table generated in the previous step to recalibrate the base quality scores of the processed bam file.
```
gatk ApplyBQSR \
       -I demo_sorted_dedup.bam \
       -R resource/Homo_sapiens_assembly38.fasta \
       --bqsr-recal-file recal_data.table \
       -O demo_sorted_dedup_recal.bam

```
### 5. Variant Calling:
Variant calling is a crucial step in clinical genomics analysis, where the goal is to identify genetic variations (variants) such as single nucleotide polymorphisms (SNPs), insertions, deletions, and structural variants in a patient's genome. These variations can be associated with diseases, provide insights into individual susceptibility to certain conditions, and guide treatment decisions.  

The GATK Variant Calling pipeline employs 2 SNV callers : [HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) and [Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2). HaplotypeCaller is used for germline variant calling, while Mutect2 is specialized for somatic variant calling in cancer samples. For this tutorial we will be focussing on HaplotypeCaller for the calling the variants from clinical genomic data.
```
gatk HaplotypeCaller \
       -I demo_sorted_dedup_recal.bam \
       -O unfiltered.vcf \
       -R resource/Homo_sapiens_assembly38.fasta \
       --dbsnp resource/Homo_sapiens_assembly38.dbsnp.vcf.gz \
       --annotation AlleleFraction \
       --annotation VariantType
```
### 6. Variant Filtering:
Variant Filtering is necessary to obtain high confidence variants from the list obtained from the variant calling step. There are two ways to filter variants in the GATK pipeline : 
(i) Variant Quality Score Recalibration and (ii) Hard-Filtering.

**(i) Variant Quality Score Recalibration (VQSR):**

This step involves applying machine learning to the variant call set to recalibrate variant quality scores based on multiple features such as depth of coverage, mapping quality, and strand bias.A model is trained using known variants and is then applied to the dataset to filter variants based on the recalibrated scores. 

**(ii) Hard Filtering:**

GATK provides recommended hard filtering criteria as an alternative to VQSR. This involves setting specific thresholds for certain variant quality metrics.
Common filters include read depth (DP), variant allele fraction (AF), strand bias, and mapping quality.

- **6.1 Splitting variants into SNPs and INDELs (optional): ** The variants called during the HaplotypeCaller step can be split into separate vcf files containing SNPs and InDels each, using the [SelectVariants](https://gatk.broadinstitute.org/hc/en-us/articles/360036362532-SelectVariants) command.
```
gatk SelectVariants \
       -V unfiltered.vcf \
       -R resource/Homo_sapiens_assembly38.fasta \
       -O unfiltered_snp.vcf \
       -select-type SNP

gatk SelectVariants \
       -V unfiltered.vcf \
       -R resource/Homo_sapiens_assembly38.fasta \
       -O unfiltered_indel.vcf \
       -select-type INDEL
```
- **6.2 Variant Quality Score Recalibration:**
- **6.3 Hard Filtering Variants:**
### 7. Variant Annotation:
### 8. Clinical Interpretation:
   



