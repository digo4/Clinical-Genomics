fastqc -o fastq_bt gatk_demo1.fastq.gz gatk_demo2.fastq.gz
fastp -i gatk_demo1.fastq.gz -I gatk_demo2.fastq.gz -o gatk_demo1_trimmed.fastq.gz -O gatk_demo2_trimmed.fastq.gz \
	-R gatk_demo -h gatk_demo.html -j gatk_demo.json --detect_adapter_for_pe
fastqc -o fastq_at gatk_demo1_trimmed.fastq.gz gatk_demo2_trimmed.fastq.gz
bwa mem ref.fa gatk_demo1_trimmed.fastq.gz gatk_demo2_trimmed.fastq.gz > demo.sam
docker run -v /home/mediomix/Desktop/Clinical_Genomics_Demo:/gatk/data -it broadinstitute/gatk:4.4.0.0
gatk AddOrReplaceReadGroups -I demo.sam -O demo.bam -RGID SRR21388960 -RGLB SRR21388960 -RGPL ILLUMINA -RGPU unit1 -RGSM SRR21388960
gatk SortSam -I demo.bam -O demo_sorted.bam -SO coordinate
gatk CollectAlignmentSummaryMetrics -R resource/Homo_sapiens_assembly38.fasta -I demo_sorted.bam -O alignment_metrics.txt
gatk MarkDuplicates -I demo_sorted.bam -O demo_sorted_dedup.bam -M marked_dup_metrics.txt --REMOVE_DUPLICATES true
gatk IndexFeatureFile -I resource/Homo_sapiens_assembly38.dbsnp.vcf.gz
gatk BaseRecalibrator -I demo_sorted_dedup.bam -R resource/Homo_sapiens_assembly38.fasta --known-sites resource/Homo_sapiens_assembly38.known_indels.vcf.gz --known-sites resource/Homo_sapiens_assembly38.dbsnp.vcf.gz -O recal_data.table
gatk ApplyBQSR -I demo_sorted_dedup.bam -R resource/Homo_sapiens_assembly38.fasta --bqsr-recal-file recal_data.table -O demo_sorted_dedup_recal.bam
screen
gatk HaplotypeCaller -I demo_sorted_dedup_recal.bam -O unfiltered.vcf -R resource/Homo_sapiens_assembly38.fasta --dbsnp resource/Homo_sapiens_assembly38.dbsnp.vcf.gz --annotation AlleleFraction --annotation VariantType
gatk VariantFiltration -V unfiltered.vcf -filter "QD < 2.0" --filter-name "QD2" \                                                                  
	-filter "QUAL < 30.0" --filter-name "QUAL30"     \
	-filter "SOR > 3.0" --filter-name "SOR3"     \
	-filter "FS > 60.0" --filter-name "FS60"     \
	-filter "MQ < 40.0" --filter-name "MQ40" \
	-O hardfiltered.vcf
bcftools filter -i '%FILTER="PASS"' hardfiltered.vcf > hardfilter_PASS.vcf
docker run -v /home/mediomix/Desktop/Clinical_Genomics_Demo:/gatk/data -it broadinstitute/gatk:4.4.0.0
gatk --java-options "-Xmx4g -Xms4g" VariantRecalibrator \
	-V hardfilter_PASS.vcf \
	--trust-all-polymorphic \
	-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
	-tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
	-an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
	-mode INDEL \
	--max-gaussians 4 \
	-resource:mills,known=false,training=true,truth=true,prior=12 resource/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
	-resource:axiomPoly,known=false,training=true,truth=false,prior=10 resource/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
	-resource:dbsnp,known=true,training=false,truth=false,prior=2 resource/Homo_sapiens_assembly38.dbsnp138.vcf \
	-O hardfiltered_indel.recal \
	--tranches-file hardfiltered_indel.tranches
gatk --java-options "-Xmx4g -Xms4g" ApplyVQSR \
	-V hardfilter_PASS.vcf \
	--recal-file hardfiltered_indel.recal \
	--tranches-file hardfiltered_indel.tranches \
	--truth-sensitivity-filter-level 99.7 \
	--create-output-variant-index true \
	-mode INDEL \
	-O hardfiltered_indel_recal.vcf
gatk --java-options "-Xmx3g -Xms3g" VariantRecalibrator \
	-V hardfiltered_indel_recal.vcf \
	--trust-all-polymorphic \
	-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 \
	-tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
	-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
	-mode SNP \
	--max-gaussians 6 \
	-resource:hapmap,known=false,training=true,truth=true,prior=15 resource/hapmap_3.3.hg38.vcf.gz \
	-resource:omni,known=false,training=true,truth=true,prior=12 resource/1000G_omni2.5.hg38.vcf.gz \
	-resource:1000G,known=false,training=true,truth=false,prior=10 resource/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
	-resource:dbsnp,known=true,training=false,truth=false,prior=7 resource/Homo_sapiens_assembly38.dbsnp138.vcf \
	-O hardfiltered_snp.recal \
	--tranches-file hardfiltered_snp.tranches \
	--rscript-file hardfiltered_snp.R
gatk --java-options "-Xmx4g -Xms4g" ApplyVQSR \
	-V hardfiltered_indel_recal.vcf \
	--recal-file hardfiltered_snp.recal \
	--tranches-file hardfiltered_snp.tranches \
	--truth-sensitivity-filter-level 99.7 \
	--create-output-variant-index true \
	-mode SNP \
	-O hardfiltered_recalibrated.vcf
bcftools filter -i '%FILTER="PASS"' hardfiltered_recalibrated.vcf > hardfiltered_recalibrated_PASS.vcf
docker run -v /home/mediomix/Desktop/Clinical_Genomics_Demo:/gatk/data -it broadinstitute/gatk:4.4.0.0
gatk SelectVariants -V hardfiltered_recalibrated_PASS.vcf -R resource/Homo_sapiens_assembly38.fasta -O hardfiltered_recalibrated_snp.vcf -select-type SNP
gatk SelectVariants -V hardfiltered_recalibrated_PASS.vcf -R resource/Homo_sapiens_assembly38.fasta -O hardfiltered_recalibrated_indel.vcf -select-type INDEL
