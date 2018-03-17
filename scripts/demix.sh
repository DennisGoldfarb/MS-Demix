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

outPath=${ROOT_OUT_DIR}/${filename}/${algName}/${lambda1}_${lambda2}_${alpha}_${deisotope}_${calcPrecursorMass}/
mkdir -p $outPath/mgf/

arrayOutPath=${ROOT_OUT_DIR}/${SLURM_ARRAY_TASK_ID}/
mkdir -p $arrayOutPath

cd ${SOURCE_DIR}/matlab/

# determine number of spectra in file
#NUMSPECTRA=`less $1 | sed -n 's/.*<spectrumList count=\"\([0-9]*\).*/\1/p'`
NUMSPECTRA=103284
#103284

printHeader=F
# process each spectrum
start=$(($SLURM_ARRAY_TASK_ID-1))

echo "STARTING SCRIPT AT:" $start

# generate model
${BUILD_DIR}/ProcessSingleSpectrum ${DATA_DIR}/${filename}.mzML ${ROOT_OUT_DIR}/${filename}/hardklor_v2/hardklor.mono.txt $start $SLURM_ARRAY_TASK_COUNT 2 3 4 $arrayOutPath

# demix spectra
param="demix('"${outPath}/mgf/"','"${arrayOutPath}"',"${SLURM_ARRAY_TASK_ID}","${SLURM_ARRAY_TASK_COUNT}","${NUMSPECTRA}",'"${algName}"',"${lambda1}","${lambda2}","${alpha}","${deisotope}","${calcPrecursorMass}","${globalTol}");quit force"
matlab -nodesktop -nosplash -r $param

# execute crux on each demixed spectrum
#${CRUX_PATH}/crux tide-search --precursor-window $globalTol --overwrite T --num-threads 1 --pin-output T --output-dir ${outPath}/crux-output/${SLURM_ARRAY_TASK_ID} ${outPath}/${SLURM_ARRAY_TASK_ID}.mgf $TIDE_INDEX 

#python ${SOURCE_DIR}/scripts/parseTide.py ${outPath}/crux-output/${SLURM_ARRAY_TASK_ID}/tide-search.target.pin $printHeader >>  ${outPath}/tide-search.${SLURM_ARRAY_TASK_ID}.pin
#printHeader=F
#python ${SOURCE_DIR}/scripts/parseTide.py ${outPath}/crux-output/${SLURM_ARRAY_TASK_ID}/tide-search.decoy.pin $printHeader >>  ${outPath}/tide-search.${SLURM_ARRAY_TASK_ID}.pin

# clean up
rm -r ${arrayOutPath}
