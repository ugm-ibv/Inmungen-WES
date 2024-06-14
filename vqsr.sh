#!/bin/sh
#SBATCH --job-name="vqsr"
#SBATCH --mem=15G
#SBATCH --cpus-per-task 8 

### Añadir archivos control para saltar funciones ###
FILE=${outDir}/recalibrated_variants.vcf
if [ -f "$FILE" ]; then
    echo "$FILE exists. Going to filtering process!"
    sbatch -p normal ${scrDir}/filter.sh
    exit 0
fi

## Filtering SNPs ##
gatk VariantRecalibrator \
   -R ${refDir}/${assembly}/Homo_sapiens_assembly.fasta \
   -V ${outDir}/raw_variants.vcf \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${refDir}/${assembly}/hapmap_3.3.vcf \
   --resource:omni,known=false,training=true,truth=true,prior=12.0 ${refDir}/${assembly}/1000G_omni2.5.vcf \
   --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${refDir}/${assembly}/1000G_phase1.snps.high_confidence.vcf \
   --resource:dbsnp,known=true,training=false,truth=false,prior=3.0 ${refDir}/${assembly}/dbsnp_138.vcf \
   -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -an MQ\
   -mode SNP \
   -O ${outDir}/output_snps.recal \
   --tranches-file ${outDir}/output_snps.tranches \
   --rscript-file ${outDir}/output_snps.plots.R


### Filtering indels ##
gatk VariantRecalibrator \
   -R ${refDir}/${assembly}/Homo_sapiens_assembly.fasta \
   -V ${outDir}/raw_variants.vcf \
   --resource:mills,known=false,training=true,truth=true,prior=12.0 ${refDir}/${assembly}/Mills_and_1000G_gold_standard.indels.vcf \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${refDir}/${assembly}/dbsnp_138.vcf \
   -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR \
   -mode INDEL \
   --max-gaussians 4 \
   -O ${outDir}/output_indels.recal \
   --tranches-file ${outDir}/output_indels.tranches \
   --rscript-file ${outDir}/output_indels.plots.R

### Applying ###
gatk ApplyVQSR \
   -R ${refDir}/${assembly}/Homo_sapiens_assembly.fasta \
   -V ${outDir}/raw_variants.vcf \
   -O ${outDir}/recalibrated_snps_raw_indels.vcf \
   --truth-sensitivity-filter-level 100.0 \
   --tranches-file ${outDir}/output_snps.tranches \
   --recal-file ${outDir}/output_snps.recal \
   -mode SNP

gatk ApplyVQSR \
   -R ${refDir}/${assembly}/Homo_sapiens_assembly.fasta \
   -V ${outDir}/recalibrated_snps_raw_indels.vcf \
   -O ${outDir}/recalibrated_variants.vcf \
   --truth-sensitivity-filter-level 100.0 \
   --tranches-file ${outDir}/output_indels.tranches \
   --recal-file ${outDir}/output_indels.recal \
   -mode INDEL


## Filtering by depth, biallelic, ALT count ##
sbatch -p normal --dependency=afterok:${SLURM_JOB_ID} ${scrDir}/filter.sh
