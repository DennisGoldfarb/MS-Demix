#!/bin/bash
module load gcc/6.3.0

SOURCE_DIR="/nas/longleaf/home/dennisg/MS-Demix/"
BUILD_DIR=${SOURCE_DIR}"/build/"
DATA_DIR="/nas/longleaf/home/dennisg/data/"

ROOT_OUT_DIR="/pine/scr/d/e/dennisg/MS-Demix/"
#ROOT_OUT_DIR="/nas/longleaf/home/dennisg/results/" 
LOG_DIR=${ROOT_OUT_DIR}"/log/"

FASTA_PATH=${SOURCE_DIR}"/FASTA/"
CRUX_PATH="/nas/longleaf/home/dennisg/crux/bin/"
TIDE_INDEX="/nas/longleaf/home/dennisg/crux/human"

algName='NNLS-sparseGroupLasso2'
#algName='NNLS-L1'
#algName='NNLS'
lambda1=0.08
lambda2=0.08
alpha=1
deisotope=1
calcPrecursorMass=1
globalTol=3
copy=0


filename=HELA_2017-10-24_CID_OT
#filename=161117_PS_SW480_1
#filename=20160302_HELA400ng_1
#filename=HELA_2018-04-02_OT
#filename=20131106_Q2_SDC_120MIN_HELA1
#filename=20131202_hela1_120min_131202214433
#filename=20131202_hela1_120min_2MZ
#filename=HELA_2018-07-09_66
