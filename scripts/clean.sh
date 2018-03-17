#!/bin/bash
#SBATCH --job-name=clean
#SBATCH -t 00-01:00:00
#SBATCH --mem=100m
#SBATCH --ntasks=1
#SBATCH --output=/pine/scr/d/e/dennisg/MS-Demix/log/clean_%A.out
#SBATCH --error=/pine/scr/d/e/dennisg/MS-Demix/log/clean_%A.err

source ./config.sh

outPath=${ROOT_OUT_DIR}/${filename}/${algName}/${lambda1}_${lambda2}_${alpha}_${deisotope}_${calcPrecursorMass}/

rm -r $outPath
