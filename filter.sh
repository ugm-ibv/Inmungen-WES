#!/bin/sh
#SBATCH --job-name="filter"

echo $(date)

for i in ALL severe moderate assymp; do
		FILE=${outDir}/recalibrated_variants_DP20_gatk_het_${i}.vcf
		if [ -f "$FILE" ]; then
		    echo "$FILE exists. Going to VEP annotation!"
			sbatch -p normal ${scrDir}/annotate_${i}.sh
		else

			### Filter variants by Depth and Genome Quality##
			gatk VariantFiltration \
		        	-R ${refDir}/${assembly}/Homo_sapiens_assembly.fasta \
			        -V ${outDir}/recalibrated_variants.vcf \
			        -O ${outDir}/recalibrated_variants_DP20_gatk.vcf \
			        -genotype-filter-name "Depth" -genotype-filter-expression "DP < 20" \
				--set-filtered-genotype-to-no-call true \
				-genotype-filter-name "GQ" -genotype-filter-expression "GQ < 20" 
		
			## Getting assymp genotypes (at least 1 called && at least 1 VAR allele & BIALLELIC)
			gatk SelectVariants \
				-R ${refDir}/${assembly}/Homo_sapiens_assembly.fasta \
				-V ${outDir}/recalibrated_variants_DP20_gatk.vcf \
				-L ${refDir}/${assembly}/Exome_V6.bed \
				-O ${outDir}/recalibrated_variants_DP20_gatk_het_assymp.vcf \
				--exclude-sample-name ${outDir}/moderate.args \
				--exclude-sample-name ${outDir}/severe.args \
			 	-select 'vc.getCalledChrCount() != 0 && (vc.getHetCount() >= 1 or vc.getHomVarCount() >= 1)'  \
				--restrict-alleles-to BIALLELIC \
				--exclude-filtered false
			
			## Getting severe genotypes (at least 1 called && at least 1 VAR allele & BIALLELIC)
			gatk SelectVariants \
			        -R ${refDir}/${assembly}/Homo_sapiens_assembly.fasta \
		        	-V ${outDir}/recalibrated_variants_DP20_gatk.vcf \
			        -L ${refDir}/${assembly}/Exome_V6.bed \
			        -O ${outDir}/recalibrated_variants_DP20_gatk_het_severe.vcf \
			        --exclude-sample-name ${outDir}/assymp.args \
		        	--exclude-sample-name ${outDir}/moderate.args \
				-select 'vc.getCalledChrCount() != 0 && (vc.getHetCount() >= 1 or vc.getHomVarCount() >= 1)' \
			        --restrict-alleles-to BIALLELIC \
			        --exclude-filtered false
		
			## Getting moderate genotypes (at least 1 called && at least 1 VAR allele & BIALLELIC)
			gatk SelectVariants \
			        -R ${refDir}/${assembly}/Homo_sapiens_assembly.fasta \
			        -V ${outDir}/recalibrated_variants_DP20_gatk.vcf \
			        -L ${refDir}/${assembly}/Exome_V6.bed \
			        -O ${outDir}/recalibrated_variants_DP20_gatk_het_moderate.vcf \
			        --exclude-sample-name ${outDir}/assymp.args \
			        --exclude-sample-name ${outDir}/severe.args \
				-select 'vc.getCalledChrCount() != 0 && (vc.getHetCount() >= 1 or vc.getHomVarCount() >= 1)' \
		        	--restrict-alleles-to BIALLELIC \
			        --exclude-filtered false

			gatk SelectVariants \
                                -R ${refDir}/${assembly}/Homo_sapiens_assembly.fasta \
                                -V ${outDir}/recalibrated_variants_DP20_gatk.vcf \
                                -L ${refDir}/${assembly}/Exome_V6.bed \
                                -O ${outDir}/recalibrated_variants_DP20_gatk_het_ALL.vcf \
                                -select 'vc.getCalledChrCount() != 0 && (vc.getHetCount() >= 1 or vc.getHomVarCount() >= 1)' \
                                --restrict-alleles-to BIALLELIC \
                                --exclude-filtered false

			sbatch -p normal --dependency=afterok:${SLURM_JOB_ID} ${scrDir}/annotate_${i}.sh
		fi
done

echo $(date)
