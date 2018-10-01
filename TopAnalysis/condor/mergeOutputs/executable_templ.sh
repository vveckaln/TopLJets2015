#!/bin/bash
CLUSTERID=$1
PROCID=$2
echo -e "@SAMPLENAME\t${PROCID}\t"`date` >> @PROJECT/condor/mergeOutputs/registry_$CLUSTERID.txt
echo "hostname " `hostname`
echo "pwd" `pwd`
echo "klist" >&2
klist >&2
echo "end klist" >&2
echo "ls EOS" >&2
ls /eos/user/v/vveckaln >&2
echo "end ls EOS" >&2

source @PROJECT/condor/mergeOutputs/run_job.sh @INPUTDIR @OUTPUTDIR @SAMPLENAME $CLUSTERID 
#EXIT_CODE=$?
if [ $(( $EXIT_CODE_SCRAM || $EXIT_CODE_PYTHON )) = 0 ];
then
sh @PROJECT/condor/mergeOutputs/success.sh $EXIT_CODE_SCRAM $EXIT_CODE_PYTHON $CLUSTERID $PROCID @SAMPLENAME
else
sh @PROJECT/condor/mergeOutputs/failure.sh $EXIT_CODE_SCRAM $EXIT_CODE_PYTHON $CLUSTERID $PROCID @SAMPLENAME
fi
