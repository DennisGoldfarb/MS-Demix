#!/bin/bash
#SBATCH --job-name=demix
#SBATCH --array=1-500
#SBATCH -t 00-02:00:00
#SBATCH --mem=4g
#SBATCH --ntasks=1
#SBATCH --output=/pine/scr/d/e/dennisg/MS-Demix/log/Deconvolution_%A_%a.out
#SBATCH --error=/pine/scr/d/e/dennisg/MS-Demix/log/Deconvolution_%A_%a.err

module load matlab

source ./config.sh

NUMSPECTRA=103284

outPath=${ROOT_OUT_DIR}/${filename}/${algName}/${lambda1}_${lambda2}_${alpha}_${deisotope}_${calcPrecursorMass}/
mkdir -p $outPath/mgf/

arrayOutPath=${ROOT_OUT_DIR}/${SLURM_ARRAY_TASK_ID}/
mkdir -p $arrayOutPath

cd ${SOURCE_DIR}/matlab/

# generate model
${BUILD_DIR}/ProcessSingleSpectrum ${BULLSEYE_DIR}/${filename}.ms2 1.6 $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT 1 5 3 $arrayOutPath

# demix spectra
param="demix('"${outPath}/mgf/"','"${arrayOutPath}"',"${SLURM_ARRAY_TASK_ID}","${SLURM_ARRAY_TASK_COUNT}","${NUMSPECTRA}",'"${algName}"',"${lambda1}","${lambda2}","${alpha}","${deisotope}","${calcPrecursorMass}","${globalTol}");quit force"
matlab -nodesktop -nosplash -r $param

# clean up
rm -r ${arrayOutPath}
