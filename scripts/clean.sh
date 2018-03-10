#!/bin/bash
#SBATCH --job-name=clean
#SBATCH -t 00-01:00:00
#SBATCH --mem=100m
#SBATCH --ntasks=1
#SBATCH --output=/pine/scr/d/e/dennisg/MS-Demix/log/clean_%A.out
#SBATCH --error=/pine/scr/d/e/dennisg/MS-Demix/log/clean_%A.err

source ./config.sh

filename=HELA_2017-10-24_CID_OT

# evaluate model in matlab
algName='NNLS-sparseGroupLasso'
lambda1=0.1
lambda2=0.1
alpha=0.5
deisotope=1
calcPrecursorMass=1
globalTol=0.02

outPath=${ROOT_OUT_DIR}/${filename}/${algName}/${lambda1}_${lambda2}_${alpha}_${deisotope}_${calcPrecursorMass}/

rm -r $outPath
