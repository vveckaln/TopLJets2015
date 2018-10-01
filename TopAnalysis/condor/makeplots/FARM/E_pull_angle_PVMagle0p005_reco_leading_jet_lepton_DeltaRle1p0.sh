#!/bin/bash
CLUSTERID=$1
PROCID=$2
echo -e "E_pull_angle_PVMagle0p005_reco_leading_jet_lepton_DeltaRle1p0\t${PROCID}\t"`date` >> /afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/condor/makeplots/registry_$CLUSTERID.txt
echo "hostname " `hostname`
echo "pwd" `pwd`
echo "klist" >&2
klist >&2
echo "end klist" >&2
echo "ls EOS" >&2
ls /eos/user/v/vveckaln >&2
echo "end ls EOS" >&2

source /afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/condor/makeplots/run_job.sh /eos/user/v/vveckaln/analysis_MC13TeV_TTJets_cflip/HADDChunks/ /eos/user/v/vveckaln/analysis_MC13TeV_TTJets_cflip/plots E_pull_angle_PVMagle0p005_reco_leading_jet_lepton_DeltaRle1p0 E_pull_angle_allconst_gen_scnd_leading_jet_leading_jet_DeltaRle1p0 cflip /eos/user/v/vveckaln/analysis_MC13TeV_TTJets/HADDChunks/MC13TeV_TTJets.root $CLUSTERID 
#EXIT_CODE=$?
if [ $(( $EXIT_CODE_SCRAM || $EXIT_CODE_PYTHON )) = 0 ];
then
sh /afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/condor/makeplots/success.sh $EXIT_CODE_SCRAM $EXIT_CODE_PYTHON $CLUSTERID $PROCID E_pull_angle_PVMagle0p005_reco_leading_jet_lepton_DeltaRle1p0
else
sh /afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/condor/makeplots/failure.sh $EXIT_CODE_SCRAM $EXIT_CODE_PYTHON $CLUSTERID $PROCID E_pull_angle_PVMagle0p005_reco_leading_jet_lepton_DeltaRle1p0
fi
