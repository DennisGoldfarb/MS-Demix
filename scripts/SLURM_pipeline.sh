#!/bin/bash
#SBATCH --job-name=Demix
#SBATCH --array=1-100
#SBATCH --time-min=60
#SBATCH --mem-per-cpu=4000
#SBATCH --ntasks=4
#SBATCH --output=/pine/scr/d/e/dennisg/MS-Demix/log/Deconvolution_%A_%a.out
#SBATCH --error=/pine/scr/d/e/dennisg/MS-Demix/log/Deconvolution_%A_%a.err

#module load matlab

#source ./config.sh

#mkdir -p ${ROOT_OUT_DIR}

# determine number of spectra in file
#NUMSPECTRA=`less $1 | sed -n 's/.*<spectrumList count=\"\([0-9]*\).*/\1/p'`
#NUMSPECTRA=3000#103284

echo "HELLOW SLURM!"$NUMSPECTRA

# process each spectrum
#i = $SLURM_ARRAY_TASK_ID
#while ($i <= $NUMSPECTRA)
    # generate model
#    ${BUILD_DIR}/ProcessSingleSpectrum ${DATA_DIR}/HELA_2017-10-24_CID_OT.mzML $i 2 3 4 ${ROOT_OUT_DIR}

    # evaluate model in matlab

    # execute crux

    # clean up
#    rm ${ROOT_OUT_DIR}/A_$i.bin
#    rm ${ROOT_OUT_DIR}/b_$i.bin
#    rm ${ROOT_OUT_DIR}/groupWeights_$i.bin
#    rm ${ROOT_OUT_DIR}/indices_$i.bin
#    rm ${ROOT_OUT_DIR}/mz_$i.tab
#    rm ${ROOT_OUT_DIR}/precursorOptions_$i.tab

#    i = $i + $SLURM_ARRAY_TASK_COUNT
#end