#!/bin/bash


user=`whoami`

if [ `hostname | cut -d "." -f2-4` == "ncg.ingrid.pt" ] ; then  
    buffer=$(qstat 2>/dev/null)
    i=0
    while read line
    do
	if [ $i -ge 2 ] ; then
	    arr=($line)
	    qdel ${arr[0]}
	fi
	((i++))
    done <<< "$buffer"
    
elif [ `hostname | cut -d "." -f2-3` == "cern.ch" ] ; then  
    buffer=$(bjobs 2>/dev/null)
    i=0
    while read line
    do
	if [ $i -ne 0 ] ; then
	    arr=($line)
	    bkill ${arr[0]}
	fi
	((i++))
    done <<< "$buffer"
    
fi