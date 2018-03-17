#!/bin/bash
#SBATCH --job-name=hardklor
#SBATCH -t 01-00:00:00
#SBATCH --mem=4g
#SBATCH --ntasks=1
#SBATCH --output=/pine/scr/d/e/dennisg/MS-Demix/log/hardklor_%A.out
#SBATCH --error=/pine/scr/d/e/dennisg/MS-Demix/log/hardklor_%A.err

source ./config.sh

filename=HELA_2017-10-24_CID_OT

outPath=${ROOT_OUT_DIR}/${filename}/hardklor_v2/

mkdir -p $outPath

${CRUX_PATH}/crux hardklor --hardklor-algorithm version2 --instrument orbitrap --resolution 85000 --overwrite T --output-dir ${outPath} ${DATA_DIR}/${filename}.mzML

