#!/bin/bash
#SBATCH --job-name=demix
#SBATCH --array=1-100
#SBATCH -t 01-00:00:00
#SBATCH --mem=16g
#SBATCH --ntasks=4
#SBATCH --output=/pine/scr/d/e/dennisg/MS-Demix/log/Deconvolution_%A_%a.out
#SBATCH --error=/pine/scr/d/e/dennisg/MS-Demix/log/Deconvolution_%A_%a.err

module load matlab

source ./config.sh

mkdir -p ${ROOT_OUT_DIR}

cd ${SOURCE_DIR}/matlab/

# determine number of spectra in file
#NUMSPECTRA=`less $1 | sed -n 's/.*<spectrumList count=\"\([0-9]*\).*/\1/p'`
NUMSPECTRA=1000
#103284

# process each spectrum
numPer=10
step=$(($SLURM_ARRAY_TASK_COUNT * $numPer))
start=$((($SLURM_ARRAY_TASK_ID-1) * $numPer + 1))
for i in `seq $start $step $NUMSPECTRA`; 
do
    echo "STARTING SCRIPT AT:"
    echo $i

    # generate model
    ${BUILD_DIR}/ProcessSingleSpectrum ${DATA_DIR}/HELA_2017-10-24_CID_OT.mzML $i $numPer 2 3 4 ${ROOT_OUT_DIR}

    for j in `seq 0 $(($numPer-1))`;
    do	
	scanID=$(($i+$j))
	
	echo "CURRENT SCAN ID:"
	echo $scanID
	if [ ! -f ${ROOT_OUT_DIR}/A_${scanID}.bin ]; then
	    continue
	fi

	# use return code to figure out if it was an MS2
	#if [[ $? != 0 ]]; then continue; fi

	# evaluate model in matlab
	algName='NNLS-sparseGroupLasso'
	lambda1=0.1
	lambda2=0.1
	alpha=0.5
	deisotope=0
	calcPrecursorMass=1

	param="demix('"${ROOT_OUT_DIR}"','"${scanID}"','"${algName}"',"${lambda1}","${lambda2}","${alpha}","${deisotope}","${calcPrecursorMass}");quit force"
	matlab -nodesktop -nosplash -r $param
    
	# execute crux

	# clean up
	rm ${ROOT_OUT_DIR}/A_${scanID}.bin
	rm ${ROOT_OUT_DIR}/b_${scanID}.bin
	rm ${ROOT_OUT_DIR}/groupWeights_${scanID}.bin
	rm ${ROOT_OUT_DIR}/indices_${scanID}.bin
	rm ${ROOT_OUT_DIR}/mz_${scanID}.tab
	rm ${ROOT_OUT_DIR}/precursorOptions_${scanID}.tab
    done
done
