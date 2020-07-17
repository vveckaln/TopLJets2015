#!/bin/bash
CLUSTERID=$1
PROCID=$2
echo -e "@BEGINNAME\t${PROCID}\t"`date` >> @PROJECT/condor/makeplots/registry_@METHOD_$CLUSTERID.txt
echo "hostname " `hostname`
echo "pwd" `pwd`
echo "klist" >&2
klist >&2
echo "end klist" >&2
echo "ls EOS" >&2
eos root://eosuser.cern.ch/ ls /eos/user/v/vveckaln >&2
echo "end ls EOS" >&2

source @PROJECT/condor/makeplots/run_job.sh @INPUTDIR @OUTPUTDIR @BEGINNAME @ENDNAME @METHOD @LISTFILE $CLUSTERID 
#EXIT_CODE=$?
if [ $(( $EXIT_CODE_SCRAM || $EXIT_CODE_PYTHON )) = 0 ];
then
sh @PROJECT/condor/makeplots/success.sh $EXIT_CODE_SCRAM $EXIT_CODE_PYTHON $CLUSTERID $PROCID @BEGINNAME
else
sh @PROJECT/condor/makeplots/failure.sh $EXIT_CODE_SCRAM $EXIT_CODE_PYTHON $CLUSTERID $PROCID @BEGINNAME
fi
