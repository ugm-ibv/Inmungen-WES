#!/bin/sh
#SBATCH --job-name="PlinkALL"
#SBATCH -c 4

if [ ! -d ${outDir}/plink ]; then
        mkdir -p ${outDir}/plink 
fi

FILE=${outDir}/plink/recalibrated_ALL_af1_clean_pruned.bim
if [ -f "$FILE" ]; then
    echo "$FILE exists. Going to ethseq!"
    sbatch -p normal ${scrDir}/ethseq.sh
    exit 0
fi

### QC
plink --vcf ${outDir}/recalibrated_variants_DP20_gatk_het_ALL_vep_reheader.vcf --double-id --make-bed --out ${outDir}/plink/recalibrated_ALL
###Check -sex
#plink --bed ${outDir}/plink/recalibrated_ALL.bed --bim ${outDir}/plink/recalibrated_ALL.bim --fam ${refDir}/fam/covid_nou_Diseasecontrol.fam --check-sex --out ${outDir}/plink/recalibrated_ALL_sex
###Duplicate samples
#plink --bed ${outDir}/plink/recalibrated_ALL.bed --bim ${outDir}/plink/recalibrated_ALL.bim --fam ${refDir}/fam/covid_nou_Diseasecontrol.fam --allow-no-sex --genome --out ${outDir}/plink/recalibrated_ALL_genome
###Heterozygosity
#plink --bed ${outDir}/plink/recalibrated_ALL.bed --bim ${outDir}/plink/recalibrated_ALL.bim --fam ${refDir}/fam/covid_nou_Diseasecontrol.fam --allow-no-sex --het --out ${outDir}/plink/recalibrated_ALL_het
###Missing
#plink --bed ${outDir}/plink/recalibrated_ALL.bed --bim ${outDir}/plink/recalibrated_ALL.bim --fam ${refDir}/fam/covid_nou_Diseasecontrol.fam --allow-no-sex --missing --out ${outDir}/plink/recalibrated_ALL_missing
### HW control
#plink --bfile ${outDir}/plink/recalibrated_ALL --keep ${refDir}/fam/controls_AD.fam --hardy --out ${outDir}/plink/recalibrated_ALL_hw_controls
###HW cases
#plink --bfile ${outDir}/plink/recalibrated_ALL --keep ${refDir}/fam/casos_AD.fam --hardy --out ${outDir}/plink/recalibrated_ALL_hw_casos
#


# Plink recalibrated ALL files
for freq in 1 5 15; do
#	python /home/jperez/scripts/retrieveID.py ${outDir}/recalibrated_variants_DP20_gatk_het_ALL_vep_af${freq}.vcf ${outDir}/recalibrated_variants_DP20_gatk_het_ALL_vep_af${freq}_mod.vcf
	plink --vcf ${outDir}/recalibrated_variants_DP20_gatk_het_ALL_vep_af${freq}.vcf --make-bed --double-id --out ${outDir}/plink/recalibrated_ALL_af${freq}

	#plink --vcf ${outDir}/recalibrated_variants_DP20_gatk_het_ALL_vep_notAF.vcf --make-bed --out ${outDir}/plink/recalibrated_ALL_notAF
	#plink --vcf ${outDir}/recalibrated_variants_DP20_gatk_het_ALL_vep_AF15.vcf --make-bed --out ${outDir}/plink/recalibrated_ALL_AF15
done

## Filtering by missing genotypes ###
	plink --bfile ${outDir}/plink/recalibrated_ALL_af5 --geno 0.25 --recode --make-bed --hwe 0.00001 --set-missing-var-ids @:# --min-ac 2 --set-hh-missing  --out ${outDir}/plink/recalibrated_ALL_af5_clean
	plink --bfile ${outDir}/plink/recalibrated_ALL_af1 --geno 0.05 --recode --make-bed --hwe 0.00001 --set-missing-var-ids @:# --min-ac 2 --set-hh-missing  --out ${outDir}/plink/recalibrated_ALL_af1_clean
#	plink --bfile ${outDir}/plink/recalibrated_ALL_notaf.bed --recode --make-bed --hwe 0.00001 --set-missing-var-ids @:# --set-hh-missing --out ${outDir}/plink/recalibrated_ALL_notaf_clean
	plink --bfile ${outDir}/plink/recalibrated_ALL_af15  --geno 0.1 --recode --make-bed --hwe 0.00001 --set-missing-var-ids @:# --min-ac 2 --set-hh-missing --out ${outDir}/plink/recalibrated_ALL_af15_clean
        plink --bfile ${outDir}/plink/recalibrated_ALL --geno 0.1 --recode --make-bed --hwe 0.00001 --set-missing-var-ids @:# --min-ac 2 --set-hh-missing  --out ${outDir}/plink/recalibrated_ALL_clean

## Filtrado por LD ###
plink --bfile ${outDir}/plink/recalibrated_ALL_af5_clean --allow-no-sex --indep-pairwise 200 50 0.2 --make-bed --out ${outDir}/plink/recalibrated_ALL_af5_clean_prune

for freq in 1 15; do
	plink --bfile ${outDir}/plink/recalibrated_ALL_af${freq}_clean --allow-no-sex --indep-pairwise 200 50 0.5 --make-bed --out ${outDir}/plink/recalibrated_ALL_af${freq}_clean_prune
done

plink --bfile ${outDir}/plink/recalibrated_ALL_clean --allow-no-sex --indep-pairwise 200 50 0.5 --make-bed --out ${outDir}/plink/recalibrated_ALL_clean_prune

for freq in 1 5 15; do
plink --bed ${outDir}/plink/recalibrated_ALL_af${freq}_clean.bed --bim ${outDir}/plink/recalibrated_ALL_af${freq}_clean.bim --fam ${refDir}/fam/covid_nou_Diseasecontrol.fam --allow-no-sex --extract ${outDir}/plink/recalibrated_ALL_af${freq}_clean_prune.prune.in --make-bed --out ${outDir}/plink/recalibrated_ALL_af${freq}_clean_pruned
done

plink --bed ${outDir}/plink/recalibrated_ALL_clean.bed --bim ${outDir}/plink/recalibrated_ALL_clean.bim --fam ${refDir}/fam/covid_nou_Diseasecontrol.fam --allow-no-sex --extract ${outDir}/plink/recalibrated_ALL_clean_prune.prune.in --make-bed --out ${outDir}/plink/recalibrated_ALL_clean_pruned

### PCA ###
plink --bfile ${outDir}/plink/recalibrated_ALL_clean_pruned --pca var-wts --allow-no-sex --out ${outDir}/plink/recalibrated_ALL_PCA
plink --bfile ${outDir}/plink/recalibrated_ALL_af5_clean_pruned --pca var-wts --allow-no-sex --out ${outDir}/plink/recalibrated_ALL_PCA_af5
plink --bfile ${outDir}/plink/recalibrated_ALL_af15_clean_pruned --pca var-wts --allow-no-sex --out ${outDir}/plink/recalibrated_ALL_PCA_af15

### Association ####
for class in SevereAssymp AssympModerate Diseasecontrol SevereModerate; do
	## Fisher ##
	plink --bed ${outDir}/plink/recalibrated_ALL_af5_clean_pruned.bed --bim ${outDir}/plink/recalibrated_ALL_af5_clean_pruned.bim --fam ${refDir}/fam/covid_nou_${class}.fam --remove ${outDir}/remove_individuals_rs.txt --model fisher --allow-no-sex --out ${outDir}/plink/${class}/recalibrated_ALL_af5_clean_pruned_fisher
	plink --bed ${outDir}/plink/recalibrated_ALL_clean_pruned.bed --bim ${outDir}/plink/recalibrated_ALL_clean_pruned.bim --fam ${refDir}/fam/covid_nou_${class}.fam --remove ${outDir}/remove_individuals_rs.txt --model fisher --allow-no-sex --out ${outDir}/plink/${class}/recalibrated_ALL_clean_pruned_fisher
	plink --bed ${outDir}/plink/recalibrated_ALL_af15_clean_pruned.bed --bim ${outDir}/plink/recalibrated_ALL_af15_clean_pruned.bim --fam ${refDir}/fam/covid_nou_${class}.fam --remove ${outDir}/remove_individuals_rs.txt --model fisher --allow-no-sex --out ${outDir}/plink/${class}/recalibrated_ALL_af15_clean_pruned_fisher


	## Logistic ##
	plink --bed ${outDir}/plink/recalibrated_ALL_af5_clean_pruned.bed --bim ${outDir}/plink/recalibrated_ALL_af5_clean_pruned.bim --fam ${refDir}/fam/covid_nou_${class}.fam --remove ${outDir}/remove_individuals_rs.txt --logistic beta --covar ${outDir}/plink/recalibrated_ALL_PCA_af5.eigenvec --adjust --allow-no-sex --out ${outDir}/plink/${class}/recalibrated_ALL_af5_clean_pruned_logistic
	plink --bed ${outDir}/plink/recalibrated_ALL_clean_pruned.bed --bim ${outDir}/plink/recalibrated_ALL_clean_pruned.bim --fam ${refDir}/fam/covid_nou_${class}.fam --remove ${outDir}/remove_individuals_rs.txt --logistic beta --covar ${outDir}/plink/recalibrated_ALL_PCA.eigenvec --adjust --allow-no-sex --out ${outDir}/plink/${class}/recalibrated_ALL_clean_pruned_logistic
	plink --bed ${outDir}/plink/recalibrated_ALL_af15_clean_pruned.bed  --bim ${outDir}/plink/recalibrated_ALL_af15_clean_pruned.bim --fam ${refDir}/fam/covid_nou_${class}.fam --remove ${outDir}/remove_individuals_rs.txt --logistic beta --covar ${outDir}/plink/recalibrated_ALL_PCA_af15.eigenvec --adjust --allow-no-sex --out ${outDir}/plink/${class}/recalibrated_ALL_af15_clean_pruned_logistic
done

sbatch -p normal --dependency=afterok:${SLURM_JOB_ID} ${scrDir}/ethseq.sh
echo $(date)
