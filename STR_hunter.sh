#!/bin/sh
#SBATCH --job-name="STR_hunter"
#SBATCH --mem=50G
#SBATCH -n 8


cp ${outDir}/${group}/${sample}/${sample}_recal_reads.bai ${outDir}/${group}/${sample}/${sample}_recal_reads.bam.bai

${scrDir}/ExpansionHunter-v5.0.0-linux_x86_64/bin/ExpansionHunter \
        --reads ${outDir}/${group}/${sample}/${sample}_recal_reads.bam \
        --reference ${refDir}/${assembly}/Homo_sapiens_assembly.fasta \
        --variant-catalog ${scrDir}/ExpansionHunter-v5.0.0-linux_x86_64/variant_catalog/${assembly}/variant_catalog.json \
        --output-prefix ${outDir}/STR/Hunter/${sample}_ExpHunter \
        --threads 8

