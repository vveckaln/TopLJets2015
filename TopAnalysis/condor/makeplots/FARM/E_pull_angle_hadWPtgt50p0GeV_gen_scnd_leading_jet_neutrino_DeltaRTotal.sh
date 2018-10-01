#!/bin/bash
CLUSTERID=$1
PROCID=$2
echo -e "E_pull_angle_hadWPtgt50p0GeV_gen_scnd_leading_jet_neutrino_DeltaRTotal\t${PROCID}\t"`date` >> /afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/condor/makeplots/registry_$CLUSTERID.txt
echo "hostname " `hostname`
echo "pwd" `pwd`
echo "klist" >&2
klist >&2
echo "end klist" >&2
echo "ls EOS" >&2
ls /eos/user/v/vveckaln >&2
echo "end ls EOS" >&2

source /afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/condor/makeplots/run_job.sh /eos/user/v/vveckaln/analysis_MC13TeV_TTJets_cflip/HADDChunks/ /eos/user/v/vveckaln/analysis_MC13TeV_TTJets_cflip/plots E_pull_angle_hadWPtgt50p0GeV_gen_scnd_leading_jet_neutrino_DeltaRTotal E_pull_angle_hadWPtle50p0GeV_gen_leading_jet_lept_w_DeltaRTotal cflip /eos/user/v/vveckaln/analysis_MC13TeV_TTJets/HADDChunks/MC13TeV_TTJets.root $CLUSTERID 
#EXIT_CODE=$?
if [ $(( $EXIT_CODE_SCRAM || $EXIT_CODE_PYTHON )) = 0 ];
then
sh /afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/condor/makeplots/success.sh $EXIT_CODE_SCRAM $EXIT_CODE_PYTHON $CLUSTERID $PROCID E_pull_angle_hadWPtgt50p0GeV_gen_scnd_leading_jet_neutrino_DeltaRTotal
else
sh /afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/condor/makeplots/failure.sh $EXIT_CODE_SCRAM $EXIT_CODE_PYTHON $CLUSTERID $PROCID E_pull_angle_hadWPtgt50p0GeV_gen_scnd_leading_jet_neutrino_DeltaRTotal
fi
