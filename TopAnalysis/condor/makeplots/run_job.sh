#! /bin/sh
WORKDIR=`pwd`
PROJECT=/afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/
EOS=/eos/user/v/vveckaln/
INPUTDIR=$1
OUTPUTDIR=$2
BEGINNAME=$3
ENDNAME=$4
METHOD=$5
LISTFILE=$6
cd /afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1
echo $PATH
MY_SCRAM=$(which scram)
echo "SCRAM is " $MY_SCRAM
scram_result=`scram r -sh`
EXIT_CODE_SCRAM=$?
echo "EXIT_CODE_SCRAM" $EXIT_CODE_SCRAM >&2
eval $scram_result
cd ${WORKDIR}
total_lumi=35874.8
python $PROJECT/condor/makeplots/plotter.py -i $INPUTDIR -O $OUTPUTDIR -j $PROJECT/data/era2016/samples.json,$PROJECT/data/era2016/qcd_samples.json  --systTheorJson=$PROJECT/data/era2016/syst_samples.json --systExpJson=$PROJECT/data/era2016/expsyst_samples.json -l $total_lumi -m $METHOD --end=$ENDNAME --begin=$BEGINNAME --outName="plotter_$BEGINNAME.root" --lfile=$LISTFILE
EXIT_CODE_PYTHON=$?
# if [ ! $EXIT_CODE_PYTHON ]; 
# then
#     sh ${PROJECT}/scripts/EOS_file_copy.sh ${WORKDIR}/${OUTPUT_FILE}.root ${OUTPUTDIR}/Chunks/${OUTPUT_FILE}.root
#     sh ${PROJECT}/scripts/EOS_file_copy.sh ${WORKDIR}/migration_${OUTPUT_FILE}.root ${OUTPUTDIR}/migration/migration_${OUTPUT_FILE}.root
# fi
#exit $EXIT_CODE_PYTHON || $EXIT_CODE_SCRAM
