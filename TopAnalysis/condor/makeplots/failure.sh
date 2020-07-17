#! /bin/sh
EXIT_CODE_SCRAM=$1
EXIT_CODE_PYTHON=$2
CLUSTERID=$3
PROCID=$4
BEGINNAME=$5
PROJECT=/afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/
echo -e "${EXIT_CODE_SCRAM}\t${EXIT_CODE_PYTHON}\t${PROCID}\t${BEGINNAME}\t\t"`hostname`"\t"`pwd`"\t"`date` >> $PROJECT/condor/makeplots/failed_jobs${CLUSTERID}_${BEGINNAME}.txt
