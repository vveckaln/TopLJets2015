#!/bin/bash
CLUSTERID=$1
PROCID=$2
echo -e "L_eta_PV_bckg_chconst_gen_lept_b\t${PROCID}\t"`date` >> /afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/condor/makeplots/registry_cflip_$CLUSTERID.txt
echo "hostname " `hostname`
echo "pwd" `pwd`
echo "klist" >&2
klist >&2
echo "end klist" >&2
echo "ls EOS" >&2
eos root://eosuser.cern.ch/ ls /eos/user/v/vveckaln >&2
echo "end ls EOS" >&2

source /afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/condor/makeplots/run_job.sh /eos/user/v/vveckaln/analysis_MC13TeV_TTJets_cflip/HADDChunks/ /eos/user/v/vveckaln/analysis_MC13TeV_TTJets_cflip/plots/ L_eta_PV_bckg_chconst_gen_lept_b L_phi_PV_hadWPtle50p0GeV_gen_leading_b cflip /eos/user/v/vveckaln/analysis_MC13TeV_TTJets/HADDChunks/MC13TeV_TTJets.root $CLUSTERID 
#EXIT_CODE=$?
if [ $(( $EXIT_CODE_SCRAM || $EXIT_CODE_PYTHON )) = 0 ];
then
sh /afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/condor/makeplots/success.sh $EXIT_CODE_SCRAM $EXIT_CODE_PYTHON $CLUSTERID $PROCID L_eta_PV_bckg_chconst_gen_lept_b
else
sh /afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/condor/makeplots/failure.sh $EXIT_CODE_SCRAM $EXIT_CODE_PYTHON $CLUSTERID $PROCID L_eta_PV_bckg_chconst_gen_lept_b
fi