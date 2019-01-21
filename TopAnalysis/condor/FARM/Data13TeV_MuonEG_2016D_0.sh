#!/bin/bash
CLUSTERID=$1
PROCID=$2
echo -e "Data13TeV_MuonEG_2016D_0\t${PROCID}\t"`date` >> /afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/condor/registry_$CLUSTERID.txt
echo "hostname " `hostname`
echo "pwd" `pwd`
echo "klist" >&2
klist >&2
echo "end klist" >&2
echo "ls EOS" >&2
ls /eos/user/v/vveckaln >&2
echo "end ls EOS" >&2

source /afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/condor/run_job.sh /eos/user/v/vveckaln/analysis_MC13TeV_TTJets root://eoscms//eos/cms//store/cmst3/group/top/ReReco2016/be52dbe_03Feb2017/Data13TeV_MuonEG_2016D/MergedMiniEvents_0_ext0.root Data13TeV_MuonEG_2016D_0 "--charge 1 --ch 13 --era /afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/data/era2016 --tag Data13TeV_MuonEG_2016D --flav 0 --method TOPJetPull::RunTopJetPull --systVar nominal --runSysts" Data13TeV_MuonEG_2016D nominal
#EXIT_CODE=$?
if [ $(( $EXIT_CODE_SCRAM || $EXIT_CODE_PYTHON )) = 0 ];
then
sh /afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/condor/success.sh $EXIT_CODE_SCRAM $EXIT_CODE_PYTHON $CLUSTERID $PROCID Data13TeV_MuonEG_2016D_0 
else
sh /afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/condor/failure.sh $EXIT_CODE_SCRAM $EXIT_CODE_PYTHON $CLUSTERID $PROCID Data13TeV_MuonEG_2016D_0
fi
