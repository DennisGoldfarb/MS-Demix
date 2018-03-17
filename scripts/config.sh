#!/bin/bash
SOURCE_DIR="/nas/longleaf/home/dennisg/MS-Demix/"
BUILD_DIR=${SOURCE_DIR}"/build/"
DATA_DIR="/nas/longleaf/home/dennisg/data/"

ROOT_OUT_DIR="/pine/scr/d/e/dennisg/MS-Demix/"
LOG_DIR=${ROOT_OUT_DIR}"/log/"

CRUX_PATH="/nas/longleaf/home/dennisg/crux/bin/"
TIDE_INDEX="/nas/longleaf/home/dennisg/crux/human"

#algName='NNLS-sparseGroupLasso'
algName='NNLS-L1'
lambda1=0.1
lambda2=0.1
alpha=0.5
deisotope=1
calcPrecursorMass=1
globalTol=0.02

filename=HELA_2017-10-24_CID_OT
