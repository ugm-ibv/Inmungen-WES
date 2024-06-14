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


## #Creación de directorios ###

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

## Tratamiento de la entrada de las opciones
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

# Creación de directorios de salida si fuese necesario
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

## Añadir archivos control para saltar funciones ###
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

# Documentación (realmente esto se podría eliminar, pero como es algo que se suele hacer
# si no da ningún problema lo voy a dejar ahí)
echo "$doc1"
echo $(date)

java -jar $EBROOTPICARD/picard.jar ReorderSam I=${outDir}/${group}/${sample}/blood01_${sample}_merged.mdup.bam O=${outDir}/${group}/${sample}/${sample}_merged.mdup.bam R=${refDir}/${assembly}/Homo_sapiens_assembly.fasta S=TRUE

rm -rf /home/jperez/COVID19/bam_files/wes/${sample}

## Base Quality Score Recalibration (BQSR) ##
gatk BaseRecalibrator \
        -R ${refDir}/${assembly}/Homo_sapiens_assembly.fasta \
        -I ${outDir}/${group}/${sample}/${sample}_merged.mdup.bam \
        --known-sites ${refDir}/${assembly}/dbsnp_138.vcf \
        --known-sites ${refDir}/${assembly}/Homo_sapiens_assembly.known_indels.vcf \
        -O ${outDir}/${group}/${sample}/${sample}_recal_data.table

gatk ApplyBQSR \
        -R ${refDir}/${assembly}/Homo_sapiens_assembly.fasta \
        -I ${outDir}/${group}/${sample}/${sample}_merged.mdup.bam \
        -bqsr ${outDir}/${group}/${sample}/${sample}_recal_data.table \
        -O ${outDir}/${group}/${sample}/${sample}_recal_reads.bam

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

## Variant Calling (GATK) ##
gatk HaplotypeCaller -ERC GVCF -R ${refDir}/${assembly}/Homo_sapiens_assembly.fasta -I ${outDir}/${group}/${sample}/${sample}_recal_reads.bam -O ${outDir}/gvcfs/${sample}.g.vcf -L ${refDir}/${assembly}/Exome_V6.bed 

## Compressing and indexing ##
bgzip -c ${outDir}/gvcfs/${sample}.g.vcf > ${outDir}/gvcfs/${sample}.g.vcf.gz
tabix -f -p vcf ${outDir}/gvcfs/${sample}.g.vcf.gz
rm ${outDir}/gvcfs/${sample}.g.vcf
rm -rf /home/jperez/COVID19/bam_files/wes/${sample}

# Combine GVCFs ##
if [ $SLURM_ARRAY_TASK_ID == 1 ]; then
      
      sbatch  -p normal ${scrDir}/combine.sh
else
      scancel ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
fi

# .args file, needed to Select_Variants ##
echo  ${sample} >> ${outDir}/${group}.args
