#!/bin/bash
#SBATCH --job-name=index
#SBATCH -t 02-00:00:00
#SBATCH --mem=16g
#SBATCH --ntasks=1
#SBATCH --output=/pine/scr/d/e/dennisg/MS-Demix/log/index_%A.out
#SBATCH --error=/pine/scr/d/e/dennisg/MS-Demix/log/index_%A.err

source ./config.sh

cd ${CRUX_PATH}/../

${CRUX_PATH}/crux tide-index --missed-cleavages 1 --mods-spec 1M+15.994914 --overwrite T ${FASTA_PATH}/human_sp_iso_020117_crap.fasta human




