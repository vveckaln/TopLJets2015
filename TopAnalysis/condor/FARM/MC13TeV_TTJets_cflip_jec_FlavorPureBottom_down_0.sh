#!/bin/bash
CLUSTERID=$1
PROCID=$2
echo -e "MC13TeV_TTJets_cflip_jec_FlavorPureBottom_down_0\t${PROCID}\t"`date` >> /afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/condor/registry_$CLUSTERID.txt
echo "hostname " `hostname`
echo "pwd" `pwd`
echo "klist" >&2
klist >&2
echo "end klist" >&2
echo "ls EOS" >&2
ls /eos/user/v/vveckaln >&2
echo "end ls EOS" >&2

source /afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/condor/run_job.sh /eos/user/v/vveckaln/analysis_MC13TeV_TTJets_cflip root://eoscms//eos/cms//store/cmst3/group/top/ReReco2016/b312177/MC13TeV_TTJets_cflip/MergedMiniEvents_0_ext0.root MC13TeV_TTJets_cflip_jec_FlavorPureBottom_down_0 "--charge 0 --ch 13 --era /afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/data/era2016 --tag MC13TeV_TTJets_cflip --flav 0 --method TOPJetPull::RunTopJetPull --systVar jec_FlavorPureBottom_down" MC13TeV_TTJets_cflip jec_FlavorPureBottom_down
#EXIT_CODE=$?
if [ $(( $EXIT_CODE_SCRAM || $EXIT_CODE_PYTHON )) = 0 ];
then
sh /afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/condor/success.sh $EXIT_CODE_SCRAM $EXIT_CODE_PYTHON $CLUSTERID $PROCID MC13TeV_TTJets_cflip_jec_FlavorPureBottom_down_0 
else
sh /afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/condor/failure.sh $EXIT_CODE_SCRAM $EXIT_CODE_PYTHON $CLUSTERID $PROCID MC13TeV_TTJets_cflip_jec_FlavorPureBottom_down_0
fi
