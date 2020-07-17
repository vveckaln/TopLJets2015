#!/bin/bash
CLUSTERID=$1
PROCID=$2
echo -e "E_mag_PV_chconst_reco_had_b\t${PROCID}\t"`date` >> /afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/condor/makeplots/registry_nominal_$CLUSTERID.txt
echo "hostname " `hostname`
echo "pwd" `pwd`
echo "klist" >&2
klist >&2
echo "end klist" >&2
echo "ls EOS" >&2
eos root://eosuser.cern.ch/ ls /eos/user/v/vveckaln >&2
echo "end ls EOS" >&2

source /afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/condor/makeplots/run_job.sh /eos/user/v/vveckaln/analysis_MC13TeV_TTJets/HADDChunks/ /eos/user/v/vveckaln/analysis_MC13TeV_TTJets/plots/ E_mag_PV_chconst_reco_had_b E_pull_angle_PFNle20_reco_scnd_leading_jet_lept_t_DeltaRgt1p0 nominal /eos/user/v/vveckaln/analysis_MC13TeV_TTJets/HADDChunks/MC13TeV_TTJets.root $CLUSTERID 
#EXIT_CODE=$?
if [ $(( $EXIT_CODE_SCRAM || $EXIT_CODE_PYTHON )) = 0 ];
then
sh /afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/condor/makeplots/success.sh $EXIT_CODE_SCRAM $EXIT_CODE_PYTHON $CLUSTERID $PROCID E_mag_PV_chconst_reco_had_b
else
sh /afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/condor/makeplots/failure.sh $EXIT_CODE_SCRAM $EXIT_CODE_PYTHON $CLUSTERID $PROCID E_mag_PV_chconst_reco_had_b
fi
