#!/bin/bash
#SBATCH --job-name=perc
#SBATCH -t 01-00:00:00
#SBATCH --mem=16g
#SBATCH --ntasks=16
#SBATCH --output=/pine/scr/d/e/dennisg/MS-Demix/log/perc_%A.out
#SBATCH --error=/pine/scr/d/e/dennisg/MS-Demix/log/perc_%A.err

source ./config.sh

outPath=${ROOT_OUT_DIR}/${filename}/${algName}/${lambda1}_${lambda2}_${alpha}_${deisotope}_${calcPrecursorMass}/

# merge mgf files
cat ${outPath}/mgf/*mgf > ${outPath}/${filename}.mgf
cat ${outPath}/mgf/*tab > ${outPath}/${filename}.tab

# execute crux demixed spectrum                                                                                                                                                     
# --mz-bin-width 0.02
${CRUX_PATH}/crux pipeline --post-processor none --score-function both --exact-p-value T --precursor-window-type mass --precursor-window $globalTol --compute-sp T --overwrite T --num-threads 16 --pin-output T --output-dir ${outPath}/ ${outPath}/${filename}.mgf --remove-precursor-peak F --remove-precursor-tolerance 0.1 --min-peaks 5 $TIDE_INDEX

${CRUX_PATH}/crux make-pin --overwrite T --output-dir ${outPath}/ ${outPath}/tide-search.target.txt

python fixPIN.py ${outPath}/make-pin.pin > ${outPath}/fixed.pin

${CRUX_PATH}/crux percolator --overwrite T --protein T --output-dir ${outPath}/ ${outPath}/fixed.pin

python statsMGF.py ${outPath}/${filename}.mgf ${outPath}/percolator.target.psms.txt ${ROOT_OUT_DIR}/${filename}/tol_${globalTol}/percolator.target.psms.txt > ${outPath}/stats.txt
