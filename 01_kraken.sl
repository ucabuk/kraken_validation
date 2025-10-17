#!/bin/bash

#===========================================================================
# slurm batch script to run the shotgun pipeline step 2
# on several sample sequentially on a fat node
# Version 0.4
#
# by Lars Harms
#
# contact: lars.harms@awi.de
#
# slurm options and variables under >set required variables<
# have to be modified by the user
#=============================================================================

#SBATCH --account=envi.envi
#SBATCH --job-name=Lama_step2_nt_0.0
#SBATCH --partition=fat
#SBATCH --time=48:00:00
#SBATCH --qos=48h
#SBATCH --cpus-per-task=32
#SBATCH --mem=700G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ugur.cabuk@awi.de

# set required variables (adapt according to your own requirements)
#===================================================================
#KRAKEN
DB="/albedo/work/projects/p_biodiv_dbs/nt_2022_10_db"
CONFIDENCE="0.6"
DEDUPE="FALSE"

# given variables (please do not change)
#===================================================================

WORK=${PWD}
OUTDIR="output"
OUT_FASTP="out.fastp"
OUT_DEDUPE="out.dedupe"
OUT_KRAKEN="out.kraken2"
OUT_KRONA="out.krona"

KRAKEN2="kraken2/2.1.2"
KRONA="krona/2.8.1"
END_R1="_fastp_R1.fq.gz"
END_R2="_fastp_R2.fq.gz"
END_MERGED="_fastp_merged_R2.fq.gz"

END_R1_DD="_fastp_dedupe__R1.fq.gz"
END_R2_DD="_fastp_dedupe__R2.f,q.gz"
END_MERGED_DD="$_fastp_depupe_merged.fq.gz"

CPU=${SLURM_CPUS_PER_TASK}

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}


# prepare environment
#===================================================================
mkdir -p ${OUTDIR}/${OUT_KRAKEN}
mkdir -p ${OUTDIR}/${OUT_KRONA}

# tasks to be performed
#===================================================================

# KRAKEN2
#----------
if [ ${DEDUPE} == "FALSE" ]; then
        module load ${KRAKEN2}
        for fq in ${OUTDIR}/${OUT_FASTP}/*${END_MERGED}
do
                BASE=${fq##*/}
                ID=${BASE%${END_MERGED}}
 srun kraken2 --confidence ${CONFIDENCE} --db ${DB} ${fq} --threads ${CPU} --gzip-compressed --output ${OUTDIR}/${OUT_KRAKEN}/${ID}_conf${CONFIDENCE}_merged.kraken --report ${OUTDIR}/${OUT_KRAKEN}/${ID}_conf${CONFIDENCE}_merged.kraken.report
                srun kraken2 --confidence ${CONFIDENCE} --db ${DB} --paired ${OUTDIR}/${OUT_FASTP}/${ID}${END_R1} ${OUTDIR}/${OUT_FASTP}/${ID}${END_R2} --threads ${CPU} --gzip-compressed --output ${OUTDIR}/${OUT_KRAKEN}/${ID}_conf${CONFIDENCE}_paired.kraken --report ${OUTDIR}/${OUT_KRAKEN}/${ID}_conf${CONFIDENCE}_paired.kraken.report
        done
        module unload ${KRAKEN2}
else
        module load ${KRAKEN2}
        for fq in ${OUTDIR}/${OUT_DEDUPE}/*${END_MERGED_DD}
        do
                BASE=${fq##*/}
                ID=${BASE%${END_MERGED_DD}}
  srun kraken2 --confidence ${CONFIDENCE} --db ${DB} ${fq} --threads ${CPU} --gzip-compressed --output ${OUTDIR}/${OUT_KRAKEN}/${ID}_conf${CONFIDENCE}_merged.kraken --report ${OUTDIR}/${OUT_KRAKEN}/${ID}_conf${CONFIDENCE}_merged.kraken.report
                srun kraken2 --confidence ${CONFIDENCE} --db ${DB} --paired ${OUTDIR}/${OUT_DEDUPE}/${ID}${END_R1_DD} ${OUTDIR}/${OUT_DEDUPE}/${ID}${END_R2_DD} --threads ${CPU} --gzip-compressed --output ${OUTDIR}/${OUT_KRAKEN}/${ID}_conf${CONFIDENCE}_paired.kraken --report ${OUTDIR}/${OUT_KRAKEN}/${ID}_conf${CONFIDENCE}_paired.kraken.report
        done
        module unload ${KRAKEN2}
fi
