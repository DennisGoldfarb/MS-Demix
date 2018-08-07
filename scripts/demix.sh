#!/bin/bash
#SBATCH --job-name=demix
#SBATCH --array=1-500
#SBATCH -t 01-00:00:00
#SBATCH --mem=6g
#SBATCH --ntasks=1
#SBATCH --output=/pine/scr/d/e/dennisg/MS-Demix/log/Deconvolution_%A_%a.out
#SBATCH --error=/pine/scr/d/e/dennisg/MS-Demix/log/Deconvolution_%A_%a.err

#module load matlab
module load matlab/2018a

source ./config.sh

outPath=${ROOT_OUT_DIR}/${filename}/${algName}/${lambda1}_${lambda2}_${alpha}_${deisotope}_${calcPrecursorMass}/
mkdir -p $outPath/mgf/

arrayOutPath=${ROOT_OUT_DIR}/${SLURM_ARRAY_TASK_ID}/
mkdir -p $arrayOutPath

cd ${SOURCE_DIR}/matlab/

# determine number of spectra in file
#NUMSPECTRA=`less $1 | sed -n 's/.*<spectrumList count=\"\([0-9]*\).*/\1/p'`
NUMSPECTRA=103284
#103284

# process each spectrum
start=$(($SLURM_ARRAY_TASK_ID-1))
echo "STARTING SCRIPT AT:" $start

# generate model
#${BUILD_DIR}/ProcessSingleSpectrum ${DATA_DIR}/${filename}.mzML ${ROOT_OUT_DIR}/${filename}/hardklor_v2/hardklor.mono.txt $start $SLURM_ARRAY_TASK_COUNT 2 5 3 $arrayOutPath
${BUILD_DIR}/ProcessSingleSpectrum ${DATA_DIR}/${filename}.mzML ${ROOT_OUT_DIR}/${filename}/hardklor_v2/hardklor.mono.txt $start 500 2 5 2 $arrayOutPath

# demix spectra
#param="demix('"${outPath}/mgf/"','"${arrayOutPath}"',"${SLURM_ARRAY_TASK_ID}","${SLURM_ARRAY_TASK_COUNT}","${NUMSPECTRA}",'"${algName}"',"${lambda1}","${lambda2}","${alpha}","${deisotope}","${calcPrecursorMass}","${globalTol}","${copy}");quit force"
param="demix('"${outPath}/mgf/"','"${arrayOutPath}"',"${SLURM_ARRAY_TASK_ID}",500,"${NUMSPECTRA}",'"${algName}"',"${lambda1}","${lambda2}","${alpha}","${deisotope}","${calcPrecursorMass}","${globalTol}","${copy}");quit force"
matlab -nodesktop -nosplash -r $param

# clean up
rm -r ${arrayOutPath}
