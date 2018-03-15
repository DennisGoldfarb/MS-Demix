#!/bin/bash
#SBATCH --job-name=crux
#SBATCH -t 01-00:00:00
#SBATCH --mem=16g
#SBATCH --ntasks=4
#SBATCH --output=/pine/scr/d/e/dennisg/MS-Demix/log/crux_%A.out
#SBATCH --error=/pine/scr/d/e/dennisg/MS-Demix/log/crux_%A.err

source ./config.sh

filename=HELA_2017-10-24_CID_OT
tol=0.02
outPath=${ROOT_OUT_DIR}/${filename}/tol_${tol}/

mkdir -p $outPath

${CRUX_PATH}/crux tide-search --precursor-window $tol --overwrite T --num-threads 4 --pin-output T --output-dir ${outPath} ${DATA_DIR}/${filename}.mzML $TIDE_INDEX

cat ${outPath}/tide-search.target.pin >> ${outPath}/target_and_decoy.pin
tail -n +2 -q ${outPath}/tide-search.decoy.pin >> ${outPath}/target_and_decoy.pin

cut -f19,20 --complement ${outPath}/target_and_decoy.pin > ${outPath}/target_and_decoy_no_mass.pin

${CRUX_PATH}/crux percolator --overwrite T --output-dir ${outPath} ${outPath}/target_and_decoy_no_mass.pin

