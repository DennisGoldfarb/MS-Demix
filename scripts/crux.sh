#!/bin/bash
#SBATCH --job-name=crux
#SBATCH -t 02-00:00:00
#SBATCH --mem=64g
#SBATCH --ntasks=8
#SBATCH --output=/pine/scr/d/e/dennisg/MS-Demix/log/crux_%A.out
#SBATCH --error=/pine/scr/d/e/dennisg/MS-Demix/log/crux_%A.err

source ./config.sh

tol=0.02
outPath=${ROOT_OUT_DIR}/${filename}/tol_${tol}_v1/

mkdir -p $outPath

${CRUX_PATH}/crux pipeline --bullseye T --retention-tolerance 0.02 --gap-tolerance 0 --scan-tolerance 1 --hardklor-algorithm version2 --instrument orbitrap --resolution 85000 --precursor-window $tol --overwrite T --num-threads 8 --pin-output T --output-dir ${outPath} ${DATA_DIR}/${filename}.mzML $TIDE_INDEX 

#${CRUX_PATH}/crux pipeline --precursor-window $tol --overwrite T --num-threads 8 --pin-output T --output-dir ${outPath} ${DATA_DIR}/${filename}.mgf $TIDE_INDEX



