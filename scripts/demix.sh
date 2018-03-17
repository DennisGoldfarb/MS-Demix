#!/bin/bash
#SBATCH --job-name=demix
#SBATCH --array=1-800
#SBATCH -t 01-00:00:00
#SBATCH --mem=4g
#SBATCH --ntasks=1
#SBATCH --output=/pine/scr/d/e/dennisg/MS-Demix/log/Deconvolution_%A_%a.out
#SBATCH --error=/pine/scr/d/e/dennisg/MS-Demix/log/Deconvolution_%A_%a.err

module load matlab

source ./config.sh

outPath=${ROOT_OUT_DIR}/${filename}/${algName}/${lambda1}_${lambda2}_${alpha}_${deisotope}_${calcPrecursorMass}/
mkdir -p $outPath

cd ${SOURCE_DIR}/matlab/

# determine number of spectra in file
#NUMSPECTRA=`less $1 | sed -n 's/.*<spectrumList count=\"\([0-9]*\).*/\1/p'`
NUMSPECTRA=103284
#103284

printHeader=T
# process each spectrum
numPer=20
step=$(($SLURM_ARRAY_TASK_COUNT * $numPer))
start=$((($SLURM_ARRAY_TASK_ID-1) * $numPer))
for i in `seq $start $step $NUMSPECTRA`; 
do
    echo "STARTING SCRIPT AT:" $i

    # generate model
    ${BUILD_DIR}/ProcessSingleSpectrum ${DATA_DIR}/${filename}.mzML ${ROOT_OUT_DIR}/${filename}/hardklor_v2/hardklor.mono.txt $i $numPer 2 3 4 ${ROOT_OUT_DIR}

    if [[ "$?" -gt 0 ]]; then
        # demix spectra
        param="demix('"${outPath}"','"${ROOT_OUT_DIR}"',"${i}","${numPer}",'"${algName}"',"${lambda1}","${lambda2}","${alpha}","${deisotope}","${calcPrecursorMass}","${globalTol}");quit force"
        matlab -nodesktop -nosplash -r $param
    fi

    for j in `seq 0 $(($numPer-1))`;
    do
        scanID=$(($i+$j))

        echo "CURRENT SCAN ID:" $scanID
        if [ ! -f ${ROOT_OUT_DIR}/A_${scanID}.bin ]; then
            continue
        fi
    
	# execute crux on each demixed spectrum
	for f in ${outPath}/${scanID}_*;
	do
	    name=$(basename "$f" .mgf)
	    tol="${name##*_}"
	    ${CRUX_PATH}/crux tide-search --precursor-window $tol --overwrite T --num-threads 1 --pin-output T --output-dir ${outPath}/crux-output/${name} $f $TIDE_INDEX 

	    python ${SOURCE_DIR}/scripts/parseTide.py ${outPath}/crux-output/${name}/tide-search.target.pin $printHeader >>  ${outPath}/tide-search.${SLURM_ARRAY_TASK_ID}.pin
	    printHeader=F
            python ${SOURCE_DIR}/scripts/parseTide.py ${outPath}/crux-output/${name}/tide-search.decoy.pin $printHeader >>  ${outPath}/tide-search.${SLURM_ARRAY_TASK_ID}.pin

	    # clean up
	    rm -r ${outPath}/crux-output/${name}
	    rm $f
	done


	# clean up
	rm ${ROOT_OUT_DIR}/A_${scanID}.bin
	rm ${ROOT_OUT_DIR}/b_${scanID}.bin
	#rm ${ROOT_OUT_DIR}/groupWeights_${scanID}.bin
	rm ${ROOT_OUT_DIR}/indices_${scanID}.bin
	rm ${ROOT_OUT_DIR}/mz_${scanID}.tab
	rm ${ROOT_OUT_DIR}/precursorOptions_${scanID}.tab
	rm ${ROOT_OUT_DIR}/${scanID}.tab
        
    done
done
