#!/bin/bash

#=============================================================================
# slurm batch script to extract reads from kraken2 analysis by taxid
#
# written by Ugur Cabuk
# ugur.cabuk@awi.de
#
# modified by Lars Harms
#
# slurm options and variables under >set required variables<
# have to be modified by the user
#=============================================================================

#SBATCH --account=envi.envi
#SBATCH --job-name=test
#SBATCH --partition=smp
#SBATCH --time=12:00:00
#SBATCH --qos=12h
#SBATCH --array=1-42%6
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ugur.cabuk@awi.de


# set required variables (adapt according to your own requirements)
#===================================================================

# taxa ID should be changed depending on what you want, General ID
# 3193 --> Embryophyta
# 35493 --> Streptophyta
# 33090 --> Viridiplantae
#3688 Salicaeae

# Therefore, this variable should be changed
TAXA="33090"

# Give path to fastp output
FASTP="output/out.fastp"

# Give path to kraken2 output
OUT_KRAKEN="output/out.kraken2"

# Give the confidence value of the kraken2 analysis
CONFIDENCE="0.6"

# Give version of krakentools
VER="1.2"



# given variables (please do not change)
#===================================================================
WORK=${PWD}

# FASTP output
END_MERGED="_fastp_merged_R2.fq.gz"
END_R1="_fastp_R1.fq.gz"
END_R2="_fastp_R2.fq.gz"

# KRAKEN output
END_KRAKEN_MERGED="_conf${CONFIDENCE}_merged.kraken"
END_KRAKEN_M_REPORT="_conf${CONFIDENCE}_merged.kraken.report"
END_KRAKEN_PAIRED="_conf${CONFIDENCE}_paired.kraken"
END_KRAKEN_P_REPORT="_conf${CONFIDENCE}_paired.kraken.report"

# Read extract output
EXT_KRAKEN="out.read_extract_${TAXA}"

END_MERGED_TAX_EXTRACT="_tax${TAXA}_extract_merged.fq"
END_R1_TAX_EXTRACT="_tax${TAXA}_extract_R1.fq"
END_R2_TAX_EXTRACT="_tax${TAXA}_extract_R2.fq"

# Read length output
LENGTH_KRAKEN="out.length_extract_${TAXA}"

LENGTH_MERGED="_merged_tax${TAXA}_extract_length.tsv"
LENGTH_PAIRED="_paired_tax${TAXA}_extract_length.tsv"

LENGTH_MERGED_TAX="_merged_tax${TAXA}_separated_extract_length.tsv"
LENGTH_R1_TAX="_paired_tax${TAXA}_separated_extract_length.tsv"
LENGTH_R2_TAX="_paired_tax${TAXA}_separated_extract_length.tsv"

# prepare environment
#===================================================================
module load krakentools/${VER}

mkdir -p ${OUT_KRAKEN}/${EXT_KRAKEN}
mkdir -p ${OUT_KRAKEN}/${LENGTH_KRAKEN}


cd ${FASTP}

R1_FILE=$(ls *${END_R1} | sed -n ${SLURM_ARRAY_TASK_ID}p)
R2_FILE=$(ls *${END_R2} | sed -n ${SLURM_ARRAY_TASK_ID}p)
MERGED=$(ls *${END_MERGED} | sed -n ${SLURM_ARRAY_TASK_ID}p)

FILEBASE=${R1_FILE%${END_R1}}

cd ${WORK}

# tasks to be performed
#===================================================================

# Extract reads by taxa
# ---------------------

__next_merged="
next inputs: ${FILEBASE}
-k: ${OUT_KRAKEN}/${FILEBASE}${END_KRAKEN_MERGED}
-s: ${FASTP}/${MERGED}
-r: ${OUT_KRAKEN}/${FILEBASE}${END_KRAKEN_M_REPORT}


command: extract_kraken_reads.py -k ${OUT_KRAKEN}/${FILEBASE}${END_KRAKEN_MERGED} -s ${FASTP}/${MERGED} -t ${TAXA} -r ${OUT_KRAKEN}/${FILEBASE}${END_KRAKEN_M_REPORT} -o ${OUT_KRAKEN}/${EXT_KRAKEN}/${FILEBASE}${END_MERGED_TAX_EXTRACT} --include-children --fastq-output
"

echo "${__next_merged}"

srun --exclusive --ntasks 1 extract_kraken_reads.py -k ${OUT_KRAKEN}/${FILEBASE}${END_KRAKEN_MERGED} -s ${FASTP}/${MERGED} -t ${TAXA} -r ${OUT_KRAKEN}/${FILEBASE}${END_KRAKEN_M_REPORT} -o ${OUT_KRAKEN}/${EXT_KRAKEN}/${FILEBASE}${END_MERGED_TAX_EXTRACT} --include-children --fastq-output &

__next_paired="
next inputs: ${FILEBASE}
-k: ${OUT_KRAKEN}/${FILEBASE}${END_KRAKEN_PAIRED}
-s1: ${FASTP}/${R1_FILE}
-s2: ${FASTP}/${R2_FILE}
-r: ${OUT_KRAKEN}/${FILEBASE}${END_KRAKEN_P_REPORT}

command: extract_kraken_reads.py -k ${OUT_KRAKEN}/${FILEBASE}${END_KRAKEN_PAIRED} -s1 ${FASTP}/${R1_FILE} -s2 ${FASTP}/${R2_FILE} -t ${TAXA} -r ${OUT_KRAKEN}/${FILEBASE}${END_KRAKEN_P_REPORT} -o ${OUT_KRAKEN}/${EXT_KRAKEN}/${FILEBASE}${END_R1_TAX_EXTRACT} -o2 ${OUT_KRAKEN}/${EXT_KRAKEN}/${FILEBASE}${END_R2_TAX_EXTRACT} --include-children --fastq-output
"

echo "${__next_paired}"

srun --exclusive --ntasks 1 extract_kraken_reads.py -k ${OUT_KRAKEN}/${FILEBASE}${END_KRAKEN_PAIRED} -s1 ${FASTP}/${R1_FILE} -s2 ${FASTP}/${R2_FILE} -t ${TAXA} -r ${OUT_KRAKEN}/${FILEBASE}${END_KRAKEN_P_REPORT} -o ${OUT_KRAKEN}/${EXT_KRAKEN}/${FILEBASE}${END_R1_TAX_EXTRACT} -o2 ${OUT_KRAKEN}/${EXT_KRAKEN}/${FILEBASE}${END_R2_TAX_EXTRACT} --include-children --fastq-output &

wait

# Extract read length by taxa/sample
# ----------------------------------

cd ${OUT_KRAKEN}/
cat ${EXT_KRAKEN}/${FILEBASE}${END_MERGED_TAX_EXTRACT} | awk '{if(NR%4==2) print length($1)}' > ${LENGTH_KRAKEN}/__tmp.${FILEBASE}_merged.len
grep @ ${EXT_KRAKEN}/${FILEBASE}${END_MERGED_TAX_EXTRACT} | awk {'print $1}' > ${LENGTH_KRAKEN}/__tmp.${FILEBASE}_merged.name
paste -d"\t" ${LENGTH_KRAKEN}/__tmp.${FILEBASE}_merged.name ${LENGTH_KRAKEN}/__tmp.${FILEBASE}_merged.len > ${LENGTH_KRAKEN}/${FILEBASE}${LENGTH_MERGED}

cat ${EXT_KRAKEN}/${FILEBASE}${END_R1_TAX_EXTRACT} | awk '{if(NR%4==2) print length($1)}' > ${LENGTH_KRAKEN}/__tmp.${FILEBASE}_R1.len
grep @ ${EXT_KRAKEN}/${FILEBASE}${END_R1_TAX_EXTRACT} | awk {'print $1}' > ${LENGTH_KRAKEN}/__tmp.${FILEBASE}_R1.name
cat ${EXT_KRAKEN}/${FILEBASE}${END_R2_TAX_EXTRACT} | awk '{if(NR%4==2) print length($1)}' > ${LENGTH_KRAKEN}/__tmp.${FILEBASE}_R2.len
grep @ ${EXT_KRAKEN}/${FILEBASE}${END_R2_TAX_EXTRACT} | awk {'print $1}' > ${LENGTH_KRAKEN}/__tmp.${FILEBASE}_R2.name
paste -d"\t" ${LENGTH_KRAKEN}/__tmp.${FILEBASE}_R1.name ${LENGTH_KRAKEN}/__tmp.${FILEBASE}_R1.len ${LENGTH_KRAKEN}/__tmp.${FILEBASE}_R2.name ${LENGTH_KRAKEN}/__tmp.${FILEBASE}_R2.len > ${LENGTH_KRAKEN}/${FILEBASE}${LENGTH_PAIRED}



cut -f1 ${LENGTH_KRAKEN}/${FILEBASE}${LENGTH_MERGED} | sed 's/@//' > ${LENGTH_KRAKEN}/__tmp.${FILEBASE}.merged.ids
grep -f ${LENGTH_KRAKEN}/__tmp.${FILEBASE}.merged.ids ${FILEBASE}${END_KRAKEN_MERGED} | cut -f2,3,4 | sort -k2,2 > ${LENGTH_KRAKEN}/${FILEBASE}${LENGTH_MERGED_TAX}
cut -f1 ${LENGTH_KRAKEN}/${FILEBASE}${LENGTH_PAIRED} | sed 's/@//' > ${LENGTH_KRAKEN}/__tmp.${FILEBASE}.R1.ids
grep -f ${LENGTH_KRAKEN}/__tmp.${FILEBASE}.R1.ids ${FILEBASE}${END_KRAKEN_MERGED} | cut -f2,3,4 | sort -k2,2 > ${LENGTH_KRAKEN}/${FILEBASE}${LENGTH_R1_TAX}
cut -f3 ${LENGTH_KRAKEN}/${FILEBASE}${LENGTH_PAIRED} | sed 's/@//' > ${LENGTH_KRAKEN}/__tmp.${FILEBASE}.R2.ids
grep -f ${LENGTH_KRAKEN}/__tmp.${FILEBASE}.R1.ids ${FILEBASE}${END_KRAKEN_MERGED} | cut -f2,3,4 | sort -k2,2 > ${LENGTH_KRAKEN}/${FILEBASE}${LENGTH_R1_TAX}


rm ${LENGTH_KRAKEN}/__tmp.*

cd ${WORK}
