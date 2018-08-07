#!/bin/bash
#SBATCH --job-name=stats
#SBATCH -t 00-02:00:00
#SBATCH --mem=4g
#SBATCH --ntasks=1
#SBATCH --output=/pine/scr/d/e/dennisg/MS-Demix/log/stats_%A.out
#SBATCH --error=/pine/scr/d/e/dennisg/MS-Demix/log/stats_%A.err

source ./config.sh

outPath=${ROOT_OUT_DIR}/${filename}/${algName}/${lambda1}_${lambda2}_${alpha}_${deisotope}_${calcPrecursorMass}/

python statsMGF.py ${outPath}/${filename}.mgf ${outPath}/percolator.target.psms.txt ${ROOT_OUT_DIR}/${filename}/tol_${globalTol}/percolator.target.psms.txt
