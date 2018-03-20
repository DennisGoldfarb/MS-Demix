#!/bin/bash
#SBATCH --job-name=bullseye
#SBATCH -t 00-04:00:00
#SBATCH --mem=4g
#SBATCH --ntasks=1
#SBATCH --output=/pine/scr/d/e/dennisg/MS-Demix/log/bullseye_%A.out
#SBATCH --error=/pine/scr/d/e/dennisg/MS-Demix/log/bullseye_%A.err

source ./config.sh

outPath=${ROOT_OUT_DIR}/${filename}/bullseye_0.5_1_3/

mkdir -p $outPath

${CRUX_PATH}/crux bullseye --retention-tolerance 0.5 --gap-tolerance 1 --scan-tolerance 3 --exact-match T --exact-tolerance 20 --hardklor-algorithm version2 --instrument orbitrap --resolution 85000 --overwrite T --output-dir ${outPath} ${DATA_DIR}/${filename}.mzML ${DATA_DIR}/${filename}.mzML

