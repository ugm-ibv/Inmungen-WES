#!/bin/bash
#SBATCH --job-name="cnv_Post"
#SBATCH --output=out/main_%A_%a.out
#SBATCH --error=err/main_%A_%a.err

source /opt/sci-soft/software/Anaconda3/5.3.0/bin/activate gatk

gatk_local PostprocessGermlineCNVCalls \
                --model-shard-path ${outDir}/cnv/gatk/cohort/${project}-model \
                --calls-shard-path ${outDir}/cnv/gatk/cohort/${project}-calls \
                --contig-ploidy-calls ${outDir}/cnv/gatk/${project}-calls \
                --output-genotyped-intervals ${outDir}/cnv/gatk/cohort/output/sample_${SLURM_ARRAY_TASK_ID}_genotyped-intervals.vcf.gz \
                --output-genotyped-segments ${outDir}/cnv/gatk/cohort/output/sample_${SLURM_ARRAY_TASK_ID}_genotyped-segments.vcf.gz \
                --sequence-dictionary ${refDir}/${assembly}/Homo_sapiens_assembly.dict \
                --sample-index  ${SLURM_ARRAY_TASK_ID} \
                --output-denoised-copy-ratios ${outDir}/cnv/gatk/cohort/output/${project}_denoised_copy_ratios.tsv \
                --allosomal-contig 1 \
                --allosomal-contig 2 \
                --allosomal-contig 3 \
                --allosomal-contig 4 \
                --allosomal-contig 5 \
                --allosomal-contig 6 \
                --allosomal-contig 7 \
                --allosomal-contig 8 \
                --allosomal-contig 9 \
                --allosomal-contig 10 \
                --allosomal-contig 11 \
                --allosomal-contig 12 \
                --allosomal-contig 13 \
                --allosomal-contig 14 \
                --allosomal-contig 15 \
                --allosomal-contig 16 \
                --allosomal-contig 17 \
                --allosomal-contig 18 \
                --allosomal-contig 19 \
                --allosomal-contig 20 \
                --allosomal-contig 21 \
		--allosomal-contig X \
                --allosomal-contig Y

/opt/sci-soft/software/Anaconda3/5.3.0/bin/deactivate
