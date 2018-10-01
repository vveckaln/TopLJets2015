#! /bin/sh
WORKDIR=`pwd`
PROJECT=/afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/
INPUTDIR=$1
OUTPUTDIR=$2
SAMPLENAME=$3
CLUSTERID=$4
mkdir -p $OUTPUTDIR/HADDChunks
mkdir -p $OUTPUTDIR/HADDmigration
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
python ${PROJECT}/condor/mergeOutputs/mergejob.py -i ${INPUTDIR}/Chunks -o ${OUTPUTDIR}/HADDChunks -s ${SAMPLENAME} -T
EXIT_CODE_PYTHON=$?
python ${PROJECT}/condor/testfilesanity.py ${OUTPUTDIR}/HADDChunks/${SAMPLENAME}.root
EXIT_CODE_PYTHON=$(($? | $EXIT_CODE_PYTHON))

if [ -f badfilereport.txt ];
then
    cat badfilereport.txt >> ${PROJECT}/condor/mergeOutputs/badfilereportChunks_${CLUSTERID}.txt
    rm badfilereport.txt
fi

python ${PROJECT}/condor/mergeOutputs/mergejob.py -i ${INPUTDIR}/migration -o ${OUTPUTDIR}/HADDmigration -s ${SAMPLENAME} 
EXIT_CODE_PYTHON=$?
mv ${OUTPUTDIR}/HADDmigration/${SAMPLENAME}.root ${OUTPUTDIR}/HADDmigration/migration_${SAMPLENAME}.root
python ${PROJECT}/condor/testfilesanity.py ${OUTPUTDIR}/HADDmigration/migration_${SAMPLENAME}.root
EXIT_CODE_PYTHON=$(($? | $EXIT_CODE_PYTHON))
if [ -f badfilereport.txt ];
then
    cat badfilereport.txt >> ${PROJECT}/condor/mergeOutputs/badfilereportmigration_${CLUSTERID}.txt
    rm badfilereport.txt
fi

echo "EXIT_CODE_PYTHON " $EXIT_CODE_PYTHON >&2
