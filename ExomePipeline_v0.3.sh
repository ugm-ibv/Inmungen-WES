#!/bin/bash

doc1="
#################################################################
#                                                               #
#       Programa para la anotación de exomas:                   #
#                                                               #
#       Pipeline para la anotación de exomas de un directorio   #
#       concreto mediante las herramientas de Funcotator,       #
#       SIFT4G y VEP a partir de los datos crudos obtenidos     #
#       por NGS.                                                #
#                                                               #
#       Versión: 0.3                                            #
#       Autor: Rubén Cañas Cañas                                #
#       Fecha: 28/07/21                                         #
#################################################################

"
doc2="opciones:

-h -> muestra la documentación
-i -> Entrada de muestra 
-d -> Directorio de la muestra a tratar (por defecto es el directorio donde se 
encuentra el programa)
-g -> Grupo de la muestra
-r -> Directorio con las secuencias de referencia (por defecto es el directorio ref
del cluster)
-s -> directorio fuente de los scripts
-o -> Directorio de salida, será donde se almacenen los resultados (por defecto es el directorio
resultados del cluster)
-a -> hg38 o hg19

## Requisitos previos a la pipeline:
# Descarga de los recursos de Funcotator
# Descarga de la base de datos de SIFT4G
# Descarga del repositorio de VEP
# Indexar el genoma de referencia (esto está incluido en la pipeline pero será necesario
# comprobar que la dirección al repositorio sea correcta)"

##### ACTUALIZACIONES DE LA VERSIÓN 0.3 #####
# Ahora existe entrada de opciones para diferentes parámetros (muestra, grupo, 
# directorios de referencia, muestras y salida, además de un comando de ayuda).
# En principio las opciones funcionan tanto con rutas absolutas como relativas, no 
# debería de dar ningún problema.
# Se eliminan todas las acciones que sean iterativas sin necesidad (e.g. indexación 
# del genoma de referencia). Ahora será NECESARIO realizarlo previamente.
# La iteración sobre las muestras se tendrá que hacer con un script que controle al programa (master/slave).
# de manera que los procesos se puedan paralelizar.


### Muchas ideas para esta pipeline han sido sacadas de la siguiente URL:
### https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/


### Entrada de opciones ###

# Tendremos valores por determinado pero se pueden cambiar fácilmente
export sample="default"
export group="default"
export sampDir=$(pwd)
export refDir=""
export outDir=""
export scrDir=""
export assembly=""
export project=""
export PATH=/home/jperez/Jordi/programes/bin:$PATH
export PYTHONPATH=/home/jperez/Jordi/programes/lib64/python3.6/site-packages
export nombre_treballs=$SLURM_ARRAY_TASK_COUNT

### Creación de directorios ###
while [ -n "$1" ]; do
	case "$1" in
		-a)
                        assembly="$2"
                        shift
                        ;;
		-s)
                        scrDir="$2"
                        shift
                        ;;
		-i)
			sample="$2"
			shift
			;;
		-g)
			group="$2"
			shift
			;;
		-d)
			sampDir=$2
			shift
			;;
		-r)
			refDir="$2"
			shift
			;;
		
		-o)
			outDir="$2"
			shift
			;;
		
		-p)
                        project="$2"
                        shift
                        ;;

		-h)
			echo "$doc1$doc2"
			exit
			;;
		--)
			shift
			break
			;;

		*) echo "Option $1 not recognized" ;;

	esac
	shift
done

# Documentación (realmente esto se podría eliminar, pero como es algo que se suele hacer
# si no da ningún problema lo voy a dejar ahí)
echo "$doc1"
echo $(date)

### Tratamiento de la entrada de las opciones ###
if [ "$sample" = "default" ]; then
	echo 'no se ha introducido ninguna muestra'
	echo 'se paraliza el proceso'
	exit
fi

if [ "$group" = 'default' ]; then
	echo 'no se ha introducido ningun grupo'
	echo 'se paraliza el proceso'
	exit
fi

### Creación de directorios de salida si fuese necesario ###
dir(){

	if [ ! -d "$outDir" ]; then
		mkdir "$outDir"
		#echo "se ha creado el directorio $outDir"
	fi

	if [ ! -d "${outDir}/${group}/${sample}" ]; then
		mkdir -p "${outDir}/${group}/${sample}"
	fi

	if [ ! -d "${outDir}/gvcfs" ]; then
		mkdir "${outDir}/gvcfs"
	fi

	if [ ! -d "${outDir}/STR" ]; then
		mkdir "${outDir}/STR"
	fi
	
	if [ ! -d "${outDir}/STR/TRED" ]; then
                mkdir -p "${outDir}/STR/TRED"
                echo "#SampleKey,BAM" >> ${outDir}/STR/TRED/TRED.csv
        fi
	
	if [ ! -d "${outDir}/STR/Hunter" ]; then
                mkdir "${outDir}/STR/Hunter"
        fi

	if [ ! -d "${outDir}/depth" ]; then
		mkdir "${outDir}/depth"
	fi

	if [ ! -d "${outDir}/trash" ]; then
                mkdir "${outDir}/trash"
        fi

	if [ ! -d "${outDir}/bam" ]; then
                mkdir "${outDir}/bam"
        fi

	if [ ! -d "${outDir}/cnv" ]; then
                mkdir -p "${outDir}/cnv/gatk"
		mkdir -p "${outDir}/cnv/exomecopy"
		mkdir -p "${outDir}/cnv/gatk/cohort"
		mkdir -p "${outDir}/cnv/gatk/cohort/output"
        fi
	
### .args file, needed to Select_Variants ###
echo  ${sample} >> ${outDir}/${group}.args

### File needed to TREDPARSE ###
echo ${sample}",${outDir}/${group}/${sample}/${sample}_recal_reads.bam" >> ${outDir}/STR/TRED/TRED.csv

}


### Establecemos un punto de control en caso de existir el archivo raw_variants.vcf ###
control(){
FILE=${outDir}/raw_variants.vcf
if [ -f "$FILE" ]; then
	if [ $SLURM_ARRAY_TASK_ID == 1 ]; then
		echo "$FILE exists. Going to VQSR"
		sbatch -p normal ${scrDir}/vqsr.sh
	else
		scancel ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
	fi
    exit 
fi
} 

## Recogemos el nombre del identificador ###
prefix=$(ls ${sampDir}/${sample}/*.gz | head -n 1 | awk -F "/" '{print $NF}' | sed 's/_R1_001.fastq.gz//g')
echo $prefix


### Ejecutamos fastp para procesar las lecturas ###
trimming(){
	
fastp -i ${sampDir}/${sample}/${prefix}_R1_001.fastq.gz -I ${sampDir}/${sample}/${prefix}_R2_001.fastq.gz -o ${outDir}/${group}/${sample}/${sample}_R1_out.fq -O ${outDir}/${group}/${sample}/${sample}_R2_out.fq --cut_tail --cut_window_size=10 --cut_mean_quality=20 --length_required=50 -c

}

### Ejecutamos bwa mem ###
mapping(){
	bwa mem \
	    -t 8 \
	    -K 100000000 \
	    -R "@RG\tID:L1\tLB:$prefix\tPL:ILLUMINA\tPM:HISEQ\tSM:${sample}" \
	    ${refDir}/${assembly}/Homo_sapiens_assembly.fasta \
	    ${outDir}/${group}/${sample}/${sample}_R1_out.fq \
	    ${outDir}/${group}/${sample}/${sample}_R2_out.fq \
	    > ${outDir}/${group}/${sample}/${sample}.sam

	# Compresss .fq files
	tar czvf ${outDir}/${group}/${sample}/${sample}_R1_out.fq
	tar czvf ${outDir}/${group}/${sample}/${sample}_R2_out.fq

	rm ${outDir}/${group}/${sample}/${sample}_R1_out.fq ${outDir}/${group}/${sample}/${sample}_R2_out.fq

	## Samtools
	samtools view -bS ${outDir}/${group}/${sample}/${sample}.sam | samtools sort -o ${outDir}/${group}/${sample}/${sample}.sorted
	rm ${outDir}/${group}/${sample}/${sample}*.sam

	## Marcamos duplicados
	java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${outDir}/${group}/${sample}/${sample}.sorted M=${outDir}/${group}/${sample}/${sample}_dedup_metrics.txt O=${outDir}/${group}/${sample}/${sample}_sorted_dedup_reads.bam
	rm ${outDir}/${group}/${sample}/*.sorted
}


### Calculem stats ###
stats(){
	module load GCC/5.4.0-2.26
	samtools depth -b ${refDir}/${assembly}/Exome_V6.bed ${outDir}/${group}/${sample}/${sample}_recal_reads.bam > ${outDir}/depth/${sample}.depth
	samtools flagstat ${outDir}/${group}/${sample}/${sample}_recal_reads.bam > ${outDir}/depth/${sample}.mapped
	bedtools coverage -hist -a ${refDir}/${assembly}/Exome_V6.bed -b ${outDir}/${group}/${sample}/${sample}_recal_reads.bam > ${outDir}/depth/${sample}.bed.cov
	grep ^all ${outDir}/depth/${sample}.bed.cov > ${outDir}/depth/${sample}.all.cov
	python /home/jperez/scripts/JordiD/reads_on_target_bedtools.py -i ${outDir}/${group}/${sample}/${sample}_recal_reads.bam -p ${refDir}/${assembly}/Exome_V6.bed -f ${outDir}/trash -o ${outDir}/depth/${sample}.ontarget
	ontarget=`cat ${outDir}/depth/${sample}.ontarget | awk '{print $NF}'`
	mean=`cat ${outDir}/depth/${sample}.depth | awk '{c++;s+=$3}END{print s/c}'`
	total_reads=`cat ${outDir}/depth/${sample}.mapped | awk -F " " 'NR == 1 {print $1}'`
	mapped=`cat ${outDir}/depth/${sample}.mapped | awk -F "[(|%]" 'NR == 5 {print $2}'`
	echo "Sample $sample,$total_reads,$mapped,$ontarget,$mean" >> ${outDir}/depth/stats.csv
	Rscript ${scrDir}/on_target.R
	rm ${outDir}/trash/*.bam ${outDir}/depth/${sample}.depth ${outDir}/depth/${sample}.bed.cov ${outDir}/depth/${sample}.all.cov ${outDir}/depth/${sample}.ontarget
}

### Base Quality Score Recalibration (BQSR) ###
bqsr(){
	gatk BaseRecalibrator \
		-R ${refDir}/${assembly}/Homo_sapiens_assembly.fasta \
		-I ${outDir}/${group}/${sample}/${sample}_sorted_dedup_reads.bam \
		--known-sites ${refDir}/${assembly}/dbsnp_138.vcf \
		--known-sites ${refDir}/${assembly}/Homo_sapiens_assembly.known_indels.vcf \
		-O ${outDir}/${group}/${sample}/${sample}_recal_data.table

	gatk ApplyBQSR \
		-R ${refDir}/${assembly}/Homo_sapiens_assembly.fasta \
		-I ${outDir}/${group}/${sample}/${sample}_sorted_dedup_reads.bam \
		-bqsr ${outDir}/${group}/${sample}/${sample}_recal_data.table \
		-O ${outDir}/${group}/${sample}/${sample}_recal_reads.bam

	rm ${outDir}/${group}/${sample}/${sample}_sorted_dedup_reads.bam

	## Base Quality Score Recalibration (BQSR) ##
	gatk BaseRecalibrator \
		-R ${refDir}/${assembly}/Homo_sapiens_assembly.fasta \
		-I ${outDir}/${group}/${sample}/${sample}_recal_reads.bam \
		--known-sites ${refDir}/${assembly}/dbsnp_138.vcf \
		--known-sites ${refDir}/${assembly}/Homo_sapiens_assembly.known_indels.vcf \
		-O ${outDir}/${group}/${sample}/${sample}_post_recal_data.table

	## Analyze Covariates ##
	gatk AnalyzeCovariates \
		-before ${outDir}/${group}/${sample}/${sample}_recal_data.table \
		-after ${outDir}/${group}/${sample}/${sample}_post_recal_data.table \
		-plots ${outDir}/${group}/${sample}/${sample}_recalibration_plots.pdf
}

### Variant Calling (GATK) ###
variant_calling(){

	gatk HaplotypeCaller -ERC GVCF -R ${refDir}/${assembly}/Homo_sapiens_assembly.fasta -I ${outDir}/${group}/${sample}/${sample}_recal_reads.bam -O ${outDir}/gvcfs/${sample}.g.vcf -L ${refDir}/${assembly}/Exome_V6.bed 

	## Compressing and indexing ##
	bgzip -c ${outDir}/gvcfs/${sample}.g.vcf > ${outDir}/gvcfs/${sample}.g.vcf.gz
	tabix -f -p vcf ${outDir}/gvcfs/${sample}.g.vcf.gz
	rm ${outDir}/gvcfs/${sample}.g.vcf
}


### CNV (GATK) ###
cnv(){
#gatk_interval=${outDir}/cnv/gatk/ExomeV6.preprocessed.interval_list
#if [ $SLURM_ARRAY_TASK_ID == 1 ]; then
#        if [ ! -f "$gatk_interval" ]; then
#           ## Binning Targets Section 1 ## 
#                       gatk PreprocessIntervals \
#                                  -R ${refDir}/${assembly}/Homo_sapiens_assembly.fasta \
#                                  -L ${refDir}/${assembly}/Exome_V6.bed \
#                                  --bin-length 0 \
#                                  --padding 250 \
#                                  -O ${outDir}/cnv/gatk/ExomeV6.preprocessed.interval_list \
#                                  -imr OVERLAPPING_ONLY
#          ## Annotate intervals Section 2 ##
#                        gatk AnnotateIntervals \
#                                -L ${outDir}/cnv/gatk/ExomeV6.preprocessed.interval_list \
#                                -R ${refDir}/${assembly}/Homo_sapiens_assembly.fasta \
#                                -imr OVERLAPPING_ONLY \
#                                -O ${outDir}/cnv/gatk/ExomeV6.annotated.tsv
#      	fi 
#			### Read counts sample 1 Section 1 ##
#			gatk CollectReadCounts \
#				-L ${outDir}/cnv/gatk/ExomeV6.preprocessed.interval_list \
#				-R ${refDir}/${assembly}/Homo_sapiens_assembly.fasta \
#				-imr OVERLAPPING_ONLY \
#				-I ${outDir}/${group}/${sample}/${sample}_recal_reads.bam \
#				--format TSV \
#				-O ${outDir}/cnv/gatk/${sample}.tsv	
#	
#else 
#	if [ ! -f "$gatk_interval" ]; then
#		sleep 60 #Dono temps al primer job per a que crei ExomeV6.preprocessed.interval_list
#		### Read counts rest of first samples Section 1 ##
#                gatk CollectReadCounts \
#                        -L ${outDir}/cnv/gatk/ExomeV6.preprocessed.interval_list \
#                        -R ${refDir}/${assembly}/Homo_sapiens_assembly.fasta \
#                        -imr OVERLAPPING_ONLY \
#                        -I ${outDir}/${group}/${sample}/${sample}_recal_reads.bam \
#                        --format TSV \
#                        -O ${outDir}/cnv/gatk/${sample}.tsv
#	else
#		### Read counts rest of samples with ExomeV6.preprocessed.interval_list already done Section 1 ##
#		gatk CollectReadCounts \
#			-L ${outDir}/cnv/gatk/ExomeV6.preprocessed.interval_list \
#			-R ${refDir}/${assembly}/Homo_sapiens_assembly.fasta \
#			-imr OVERLAPPING_ONLY \
#			-I ${outDir}/${group}/${sample}/${sample}_recal_reads.bam \
#			--format TSV \
#			-O ${outDir}/cnv/gatk/${sample}.tsv
#	fi
#fi

if [ $SLURM_ARRAY_TASK_ID == 1 ]; then
	sbatch -p normal --dependency=afterok:${SLURM_ARRAY_JOB_ID} ${scrDir}/cnv_det_germ_prova.sh
fi 

}

### STR ###
str_hunter(){
# ExpansionHunter
cp ${outDir}/${group}/${sample}/${sample}_recal_reads.bai ${outDir}/${group}/${sample}/${sample}_recal_reads.bam.bai

${scrDir}/ExpansionHunter-v5.0.0-linux_x86_64/bin/ExpansionHunter \
	--reads ${outDir}/${group}/${sample}/${sample}_recal_reads.bam \
	--reference ${refDir}/${assembly}/Homo_sapiens_assembly.fasta \
	--variant-catalog ${scrDir}/ExpansionHunter-v5.0.0-linux_x86_64/variant_catalog/${assembly}/variant_catalog.json \
	--output-prefix ${outDir}/STR/Hunter/${sample}_ExpHunter \
	--threads 8
}

### Stretch ###
stretch(){
ln -s ${outDir}/${group}/${sample}/*_recal_reads.bam ${outDir}/bam
ln -s ${outDir}/${group}/${sample}/*_recal_reads.bam.bai ${outDir}/bam
if [ ${assembly} == "hg19" ]
	then
	bpipe run -p input_regions=${scrDir}/${assembly}/STR/hg19.simpleRepeat_period1-6.dedup.sorted.bed ${scrDir}/STRetch/pipelines/STRetch_exome_bam_pipeline_hg19.groovy ${outDir}/bam/*.bam

elif [ ${assembly} == "hg38" ]
then
	bpipe run -p input_regions=${scrDir}/${assembly}/STR/hg38.simpleRepeat_period1-6.dedup.sorted.bed ${scrDir}/STRetch/pipelines/STRetch_exome_bam_pipeline_h38.groovy ${outDir}/bam/*.bam
fi

}

### Combine GVCFs ###
combine(){
	if [ $SLURM_ARRAY_TASK_ID == 1 ]; then
	      
	      sbatch -p normal --dependency=afterok:${SLURM_ARRAY_JOB_ID} ${scrDir}/combine.sh
	fi
}


# ----------------------------------------------------
# ------- WES pipeline --------------------
# -----------------------------------------------------

dir ${@}
#control ${@}
#trimmming ${@}
#mapping ${@}
#stats ${@}
#bqsr ${@}
#variant_calling ${@}
cnv ${@}
#combine ${@}
#str_hunter ${@}
#stretch ${@}
