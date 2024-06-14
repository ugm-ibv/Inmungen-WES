#!/bin/bash
#SBATCH --job-name=ExomeCOVID
#SBATCH --cpus-per-task 8
#SBATCH --mem 100G
#SBATCH --array=1-196%5
#SBATCH --output=out/main_%A_%a.out
#SBATCH --error=err/main_%A_%a.err

module load GATK/4.1.7.0-Java-1.8.0_92
module load BWA/0.7.17-foss-2016b
module load BCFtools/1.10.2-foss-2016b
module load SAMtools/1.10-foss-2016b
module load DBD-mysql/4.039-foss-2016b-Perl-5.24.0
module load VEP/104.3-foss-2016b-Perl-5.24.0
module load fastp/0.20.1-foss-2016b
module load picard/2.9.2-Java-1.8.0_92
module load tabix/0.2.6-foss-2016b
module load BEDTools/2.29.2-foss-2016b
module load SIFT4G/2.4-Java-1.8.0_92
module load R/3.6.2-foss-2016b-X11-20160819
module load PLINK/1.9b_6.24-x86_64
module load FastQC/0.11.8-Java-1.8.0_92
module load MultiQC/1.7-foss-2016b-Python-2.7.12
module load Anaconda3/5.3.0
module load Pysam/0.9.0-goolf-1.7.20-Python-2.7.11

SEEDFILE=/home/jperez/COVID19/Exome_scripts/com_file_prova.txt
SEED=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SEEDFILE)

$SEED
