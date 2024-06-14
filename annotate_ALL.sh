#!/bin/sh
#SBATCH --job-name="AnnCOVALL"
#SBATCH -c 4
#SBATCH --mem=100G


FILE=${outDir}/recalibrated_variants_DP20_gatk_het_ALL_vep_af5.vcf
if [ -f "$FILE" ]; then
    echo "$FILE exists. Going to PLINK!"
    sbatch -p normal ${scrDir}/plink.sh
    exit 0
fi

if [ ${assembly} == "hg19" ]
then
        lofteedir=/home/jperez/.vep/Plugins/loftee_37
        lof_database=phylocsf_gerp.sql
	port=3337
elif [ ${assembly} == "hg38" ]
then
        lofteedir=/home/jperez/.vep/Plugins/loftee_38
        lof_database=loftee.sql
	port=5306
fi

export PERL5LIB=${PERL5LIB}:${lofteedir}

# VEP ##
vep \
    --cache \
    --dir  ${refDir}/${assembly} \
    --port ${port} \
    --sift b \
    --fork 4 \
    --variant_class \
    --polyphen b \
    --humdiv \
    --gene_phenotype \
    --regulatory \
    --numbers \
    --total_length \
    --keep_csq \
     --vcf \
    --vcf_info_field CSQ \
    --symbol \
    --uniprot \
    --check_existing \
    --clin_sig_allele 1 \
    --af_1kg \
    --af_esp \
    --af_gnomad \
    --canonical \
    --distance 300 \
    --pubmed \
     --force_overwrite \
    --plugin REVEL,/home/jperez/.vep/new_tabbed_revel_${assembly}.tsv.gz \
    --plugin GeneSplicer,/home/jperez/Jordi/programes/GeneSplicer/sources/genesplicer,/home/jperez/Jordi/programes/GeneSplicer/human \
    --plugin SpliceAI,snv=/home/jperez/ref/JordiD/spliceai_scores.raw.snv.${assembly}.vcf.gz,indel=/home/jperez/ref/JordiD/spliceai_scores.raw.indel.${assembly}.vcf.gz \
    --plugin CADD,/home/jperez/.vep/whole_genome_SNVs_${assembly}.tsv.gz \
    --plugin Mastermind,${refDir}/${assembly}/mastermind_cited_variants_reference-2022.04.02-${assembly}.vcf.gz \
    -i ${outDir}/recalibrated_variants_DP20_gatk_het_ALL.vcf \
    -o ${outDir}/recalibrated_variants_DP20_gatk_het_ALL_vep.vcf


# Filter VEP ##
filter_vep \
	-i ${outDir}/recalibrated_variants_DP20_gatk_het_ALL_vep.vcf \
	-o ${outDir}/gene_conseq_variants_DP20_gatk_het_ALL_vep.vcf \
	--format vcf \
	--only_matched \
	--filter "Consequence in /home/jperez/ref/consequences.txt" \
	--filter "SYMBOL in /home/jperez/ref/genes.txt" \
	--filter "CANONICAL is YES" \
	--force_overwrite


## Just filtering by AF ###
filter_vep \
	-i ${outDir}/recalibrated_variants_DP20_gatk_het_ALL_vep.vcf \
	-o ${outDir}/recalibrated_variants_DP20_gatk_het_ALL_vep_af5.vcf \
	--format vcf \
	--only_matched \
	--filter "Consequence in /home/jperez/ref/consequences.txt" \
	--filter "gnomAD_NFE_AF >= 0.05 or gnomAD_AMR_AF >= 0.05" \
	--filter "CANONICAL is YES" \
	--force_overwrite

filter_vep \
	-i ${outDir}/recalibrated_variants_DP20_gatk_het_ALL_vep.vcf \
	-o ${outDir}/recalibrated_variants_DP20_gatk_het_ALL_vep_af1.vcf \
	--format vcf \
	--only_matched \
	--filter "Consequence in /home/jperez/ref/consequences.txt" \
	--filter "gnomAD_NFE_AF <= 0.01 or gnomAD_AMR_AF <= 0.01" \
	--filter "CANONICAL is YES" \
	--force_overwrite

filter_vep \
	-i ${outDir}/recalibrated_variants_DP20_gatk_het_ALL_vep.vcf \
	-o ${outDir}/recalibrated_variants_DP20_gatk_het_ALL_vep_af15.vcf \
	--format vcf \
	--only_matched \
	--filter "Consequence in /home/jperez/ref/consequences.txt" \
	--filter "(gnomAD_NFE_AF < 0.05 and gnomAD_NFE_AF > 0.01) or (gnomAD_AMR_AF < 0.05 and gnomAD_AMR_AF > 0.01)"\
	--filter "CANONICAL is YES" \
	--force_overwrite

filter_vep \
        -i ${outDir}/recalibrated_variants_DP20_gatk_het_ALL_vep.vcf \
        -o ${outDir}/recalibrated_variants_DP20_gatk_het_ALL_vep_notaf.vcf \
        --format vcf \
        --only_matched \
        --filter "Consequence in /home/jperez/ref/consequences.txt" \
        --filter "not gnomAD_NFE_AF or not gnomAD_AMR_AF"\
        --filter "CANONICAL is YES" \
        --force_overwrite

sbatch -p normal --dependency=afterok:${SLURM_JOB_ID} ${scrDir}/plink.sh

echo $(date)
