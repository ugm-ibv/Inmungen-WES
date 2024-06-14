#!/bin/sh
#SBATCH --job-name="EthSEQ"
#SBATCH --mem=100G
#SBATCH -n 4

module load R/4.0.2-foss-2016b-X11-20160819
module load BCFtools/1.3-foss-2016b


if [ ! -d "${outDir}/ethseq" ]; then
        mkdir "${outDir}/ethseq"
fi

FILE=${outDir}/ethseq/Report.pdf
if [ -f "$FILE" ]; then
    echo "$FILE exists."
    sbatch -p normal ${scrDir}/tred.sh
    exit 0
fi

bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\tGT\t.[\t%GT]\n' ${outDir}/recalibrated_variants_DP20_gatk_het_ALL_vep_reheader.vcf > ${outDir}/ethseq/recalibrated_variants_DP20_gatk_het_ALL_vep_modified_temp.vcf
sed 's/|/\//g' ${outDir}/ethseq/recalibrated_variants_DP20_gatk_het_ALL_vep_modified_temp.vcf > ${outDir}/ethseq/recalibrated_variants_DP20_gatk_het_ALL_vep_modified.vcf
rm ${outDir}/ethseq/recalibrated_variants_DP20_gatk_het_ALL_vep_modified_temp.vcf

if [ ${assembly} == "hg19" ]
then
	Rscript ${scrDir}/ethseq.R

elif [ ${assembly} == "hg38" ]
then
	Rscript ${scrDir}/ethseq_38.R
fi

sbatch -p normal ${scrDir}/tred.sh
