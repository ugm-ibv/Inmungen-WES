#!/bin/sh
#SBATCH --job-name="TRED"
#SBATCH --mem=50G
#SBATCH -n 4

export PYTHONPATH=/home/jperez/Jordi/programes/lib/python2.7/site-packages

FILE=${outDir}/STR/TRED/TRED.tsv.report.txt
if [ -f "$FILE" ]; 
	then
        	echo "$FILE exists. Going to Hunter"
		exit 0
#                sbatch -p normal ${scrDir}/STR_hunter.sh
else

	tred.py ${outDir}/STR/TRED/TRED.csv --ref ${assembly}_nochr --workdir $outDir/STR/TRED --log INFO --cpus 4 
	tredreport.py $outDir/STR/TRED/*.json --tsv $outDir/STR/TRED/TRED.tsv --ref ${assembly}_nochr --cpus 4
	
fi
