#! /bin/sh
WORKDIR=`pwd`
PROJECT=/afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/
EOS_STORE_PREFIX=root://eoscms//eos/cms//store/cmst3/group/top/ReReco2016/b312177/
EOS=/eos/user/v/vveckaln/
OUTPUTDIR=$1
INPUT_FILE=$2
OUTPUT_FILE=$3
RUN_OPTS=$4
TAG=$5
SYS=$6
cd /afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1
echo $PATH
MY_SCRAM=$(which scram)
echo "SCRAM is " $MY_SCRAM
scram_result=`scram r -sh`
EXIT_CODE_SCRAM=$?
echo "EXIT_CODE_SCRAM" $EXIT_CODE_SCRAM >&2
eval $scram_result
cd ${WORKDIR}
#valgrind --track-origins=yes --leak-check=yes --suppressions=$ROOTSYS/etc/valgrind-root.supp 
python ${PROJECT}/condor/runLocalAnalysis.py -i ${INPUT_FILE} -o ${WORKDIR}/${OUTPUT_FILE}.root ${RUN_OPTS}
EXIT_CODE_PYTHON=$?
python ${PROJECT}/condor/testfilesanity.py ${WORKDIR}/${OUTPUT_FILE}.root
python ${PROJECT}/condor/testfilesanity.py ${WORKDIR}/migration_${OUTPUT_FILE}.root
EXIT_CODE_PYTHON=$(($? + $EXIT_CODE_PYTHON))
echo "EXIT_CODE_PYTHON " $EXIT_CODE_PYTHON >&2

DIR=""
if [[ $SYS == "nominal" ]];
then
    DIR=$TAG
else
    DIR=${TAG}_${SYS}
fi
if [ $EXIT_CODE_PYTHON ]; 
eos root://eosuser.cern.ch// mkdir -p ${OUTPUTDIR}/Chunks/${DIR}
eos root://eosuser.cern.ch// mkdir -p ${OUTPUTDIR}/migration/${DIR}
then
    sh ${PROJECT}/scripts/EOS_file_copy.sh ${WORKDIR}/${OUTPUT_FILE}.root ${OUTPUTDIR}/Chunks/${DIR}/${OUTPUT_FILE}.root
    sh ${PROJECT}/scripts/EOS_file_copy.sh ${WORKDIR}/migration_${OUTPUT_FILE}.root ${OUTPUTDIR}/migration/${DIR}/migration_${OUTPUT_FILE}.root
fi
#exit $EXIT_CODE_PYTHON || $EXIT_CODE_SCRAM
