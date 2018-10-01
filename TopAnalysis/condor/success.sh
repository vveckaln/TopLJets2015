#! /bin/sh
EXIT_CODE_SCRAM=$1
EXIT_CODE_PYTHON=$2
CLUSTERID=$3
PROCID=$4
JOB_TAG=$5
PROJECT=/afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/
echo -e "$1\t$2\t$4\t$5\t\t"`hostname`"\t"`pwd`"\t"`date` >> $PROJECT/condor/successful_jobs$3.txt
