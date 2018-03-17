#!/bin/bash
#SBATCH --job-name=bullseye
#SBATCH -t 00-02:00:00
#SBATCH --mem=4g
#SBATCH --ntasks=1
#SBATCH --output=/pine/scr/d/e/dennisg/MS-Demix/log/bullseye_%A.out
#SBATCH --error=/pine/scr/d/e/dennisg/MS-Demix/log/bullseye_%A.err

source ./config.sh

outPath=${ROOT_OUT_DIR}/${filename}/bullseye/

mkdir -p $outPath

${CRUX_PATH}/crux bullseye --exact-match T --exact-tolerance 20 --hardklor-algorithm version2 --instrument orbitrap --resolution 85000 --overwrite T --output-dir ${outPath} ${DATA_DIR}/${filename}.mzML ${DATA_DIR}/${filename}.mzML

