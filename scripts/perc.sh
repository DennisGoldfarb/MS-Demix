#!/bin/bash
#SBATCH --job-name=perc
#SBATCH -t 01-00:00:00
#SBATCH --mem=16g
#SBATCH --ntasks=4
#SBATCH --output=/pine/scr/d/e/dennisg/MS-Demix/log/perc_%A.out
#SBATCH --error=/pine/scr/d/e/dennisg/MS-Demix/log/perc_%A.err

source ./config.sh

filename=HELA_2017-10-24_CID_OT

algName='NNLS-sparseGroupLasso'
lambda1=0.1
lambda2=0.1
alpha=0.5
deisotope=1
calcPrecursorMass=1

outPath=${ROOT_OUT_DIR}/${filename}/${algName}/${lambda1}_${lambda2}_${alpha}_${deisotope}_${calcPrecursorMass}/

# merge tide-search results
skip=1
for f in ${outPath}/tide-search.*.pin;
do
    tail -n +${skip} -q $f >> ${outPath}/target_and_decoy.pin
    #rm $f
    skip=2
done

cut -f19,20 --complement ${outPath}/target_and_decoy.pin > ${outPath}/target_and_decoy_no_mass.pin

# execute percolator
${CRUX_PATH}/crux percolator --overwrite T --output-dir ${outPath} ${outPath}/target_and_decoy_no_mass.pin

