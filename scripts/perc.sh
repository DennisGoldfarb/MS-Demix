#!/bin/bash
#SBATCH --job-name=perc
#SBATCH -t 00-02:00:00
#SBATCH --mem=16g
#SBATCH --ntasks=4
#SBATCH --output=/pine/scr/d/e/dennisg/MS-Demix/log/perc_%A.out
#SBATCH --error=/pine/scr/d/e/dennisg/MS-Demix/log/perc_%A.err

source ./config.sh

outPath=${ROOT_OUT_DIR}/${filename}/${algName}/${lambda1}_${lambda2}_${alpha}_${deisotope}_${calcPrecursorMass}/

# merge mgf files
cat ${outPath}/mgf/*mgf > ${outPath}/${filename}.mgf

# execute crux demixed spectrum                                                                                                                                                     
${CRUX_PATH}/crux pipeline --precursor-window $globalTol --overwrite T --num-threads 4 --pin-output T --output-dir ${outPath}/ ${outPath}/${filename}.mgf $TIDE_INDEX


