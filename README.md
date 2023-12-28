[![Static Badge](https://img.shields.io/badge/LICENSE-MIT-yellow)](https://opensource.org/license/mit/) [![Static Badge](https://img.shields.io/badge/FASTQC-v.0.12.0-blue)](https://github.com/s-andrews/FastQC) [![Static Badge](https://img.shields.io/badge/fastp-v.0.20.1-blue)](https://github.com/OpenGene/fastp/releases)  [![Static Badge](https://img.shields.io/badge/bwa-v.0.7.17-blue)](https://github.com/lh3/bwa)  [![Static Badge](https://img.shields.io/badge/java-%3E%3Dv.17-pink)
](https://openjdk.org/projects/jdk/17/) [![Static Badge](https://img.shields.io/badge/gatk-%3E%3Dv.4.4.0.0-teal)
](https://hub.docker.com/r/broadinstitute/gatk/) [![Static Badge](https://img.shields.io/badge/snpEff%26SnpSift-%3Dv.5.2-purple)
](https://pcingola.github.io/SnpEff/)    [![Static Badge](https://img.shields.io/badge/InterVar-%3E%3Dv.2.2.2-green)](https://github.com/WGLab/InterVar)





## Introduction
Clinical genomics analysis refers to the application of genomic technologies and computational methods to understand the genetic basis of diseases and inform clinical decision-making. It involves the analysis of genomic data obtained from patients to identify genetic variations, mutations, and other molecular features that may contribute to disease development or progression. The goal is to use this information to guide diagnosis, treatment decisions, and personalized medicine approaches.

Key steps and considerations for clinical genomics analysis :

- **Sample Collection and DNA Extraction:**

Obtain biological samples (e.g., blood, tissue) from patients. Extract genomic DNA from the samples using standardized procedures.

- **Library Preparation and Sequencing:**

Prepare DNA libraries suitable for sequencing. Perform next-generation sequencing (NGS) to generate genomic data. Common techniques include whole-genome sequencing (WGS) or targeted sequencing of specific genomic regions.

- **Data Quality Control:**

Assess the quality of raw sequencing data. Perform quality control steps, including checking for sequencing errors, read depth, coverage uniformity, and other metrics.

- **Read Alignment and Variant Calling:**

Align sequenced reads to a reference genome. Identify genetic variants (e.g., single nucleotide variants, insertions, deletions) using variant calling algorithms.

- **Variant Annotation:**

Annotate identified variants with information about their functional impact, population frequency, and association with known diseases. Use databases such as dbSNP, ClinVar, and others.

- **Variant Filtering:**

Apply filtering criteria to prioritize variants based on their relevance to the patient's phenotype and the suspected disease. Consider factors like variant type, predicted pathogenicity, and allele frequency.

- **Interpretation and Clinical Reporting:**

Review and interpret the filtered variants in the context of the patient's clinical history and phenotype. Determine the clinical significance of each variant. Generate a comprehensive clinical report.

- **Clinical Validation:**

Validate clinically relevant variants using orthogonal methods, such as Sanger sequencing or other established techniques. Confirm the presence of pathogenic variants.

- **Variant Prioritization:**

Prioritize variants based on their clinical significance, potential impact on disease, and relevance to treatment options. Identify actionable variants that may guide therapeutic interventions.

- **Clinical Decision Support:**

Provide clinical decision support to healthcare providers based on the genomic findings. Offer recommendations for patient management, including potential therapeutic options and surveillance strategies.

- **Genetic Counseling:**

Engage genetic counselors to communicate genomic results to patients and their families. Provide information about the implications of identified variants, recurrence risks, and potential preventive measures.

- **Integration with Electronic Health Records (EHR):**

Integrate genomic data and clinical interpretations into the patient's electronic health record (EHR) to facilitate ongoing healthcare management and collaboration among healthcare providers.

- **Continuous Learning and Research:**

Stay informed about the latest genomic research and clinical guidelines. Continuously update the interpretation pipeline based on evolving scientific knowledge.

- **Ethical and Legal Considerations:**

Address ethical and legal considerations related to genomic data, including patient consent, privacy, and adherence to regulatory requirements such as the Health Insurance Portability and Accountability Act (HIPAA).

- **Follow-Up and Long-Term Monitoring:**

Establish a follow-up plan for long-term monitoring of patients based on genomic findings. Periodically reassess and update clinical recommendations as new information becomes available.

## A step-by-step guide to Clinical Genomics analysis:
1. **Pre-processing reads:** 
Pre-processing the reads involves checking the quality (Phred-scores) using tools like [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) or [MultiQC](https://github.com/MultiQC/MultiQC). After that, we need to remove adapter contamination and trim quality reads using either one of the following tools : [fastp](https://github.com/OpenGene/fastp?tab=readme-ov-file), [cutadapt](https://cutadapt.readthedocs.io/en/stable/), [TrimGalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) or [Trimmomatic](/http://www.usadellab.org/cms/?page=trimmomatic).

```
mkdir fastq_bt
#Running FASTQC
fastqc -o fastq_bt gatk_demo1.fastq.gz gatk_demo2.fastq.gz

```

3. **Aligning reads to reference genome:**
4. Alignment Post-Processing
5. Base Quality score recalibration (optional)
6. Variant Calling
7. Variant Filtering
   - 6.1 Splitting variants into SNPs and INDELs
   - 6.2 Variant Quality Score Recalibration
   - 6.3 Hard Filtering Variants
8. Variant Annotation
9. Clinical Interpretation
   



