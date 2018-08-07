#!/bin/bash
#SBATCH --job-name=crux
#SBATCH -t 02-00:00:00
#SBATCH --mem=32g
#SBATCH --ntasks=8
#SBATCH --output=/pine/scr/d/e/dennisg/MS-Demix/log/crux_%A.out
#SBATCH --error=/pine/scr/d/e/dennisg/MS-Demix/log/crux_%A.err

source ./config.sh

outPath=${ROOT_OUT_DIR}/${filename}/tol_${globalTol}/

mkdir -p $outPath

#${CRUX_PATH}/crux pipeline --bullseye T --retention-tolerance 0.02 --gap-tolerance 0 --scan-tolerance 1 --hardklor-algorithm version2 --instrument orbitrap --resolution 85000 --precursor-window $tol --overwrite T --num-threads 8 --pin-output T --output-dir ${outPath} ${DATA_DIR}/${filename}.mzML $TIDE_INDEX 

#--mz-bin-width 0.02

${CRUX_PATH}/crux pipeline --post-processor none --precursor-window-type mass --precursor-window $globalTol --compute-sp T --overwrite T --num-threads 8 --pin-output T --output-dir ${outPath} ${DATA_DIR}/${filename}.mzML --min-peaks 5 --score-function both --exact-p-value T $TIDE_INDEX

${CRUX_PATH}/crux make-pin --overwrite T --output-dir ${outPath}/ ${outPath}/tide-search.target.txt

python fixPIN.py ${outPath}/make-pin.pin > ${outPath}/fixed.pin

${CRUX_PATH}/crux percolator --overwrite T --protein T --output-dir ${outPath}/ ${outPath}/fixed.pin



