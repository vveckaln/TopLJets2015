#!/bin/bash


# "grep vischia" selects only the lines corresponding to a job. It might as well be "grep job-" or stuff like that.
user=`whoami`
echo "RUNNING        PENDING      TOTAL"
if [ `hostname | cut -d "." -f2-4` == "ncg.ingrid.pt" ] ; then  

    if [ "$1" == "" ]; then
	while true ;
	do
	    echo "running: " `qstat -u $user | grep $user | grep " r " | wc -l` "    total: " `qstat -u $user | grep $user | wc -l` 


	done
    else
	while true;
	do
	    buffer=$(qstat 2> /dev/null)
	    i=0
	    running=0
	    total=0
	    while read line ;
	    do
		if [ $i -ge 2 ]; 
		then
		    arr=($line)
		    job_number=${arr[0]}
		    status=${arr[4]}
		    l=`qstat -j $job_number 2> /dev/null | sed '/job_name:/!d'`
		    if [ "$l" != "" ]; then
			l_arr=($l)
			job_name=${l_arr[1]}
			dtag=`echo $job_name | sed 's/^.\{,9\}//' | sed 's/.....$//g'`
#	       dtag=`echo $dtag1 ' `
#| sed 's/.....$//'`
			if [ "$1" == "$dtag" ]; 
			then
			    ((total++))
			    if [ "$status" == "r" ];
			    then
				((running++))
			    fi
			fi
		    fi
		fi
		((i++))
	    done <<<"$buffer"
	    echo $running"     " $total
	done

    fi
    
elif [ `hostname | cut -d "." -f2-3` == "cern.ch" ] ; then  
    while true ;
    do
        
	#berr=$(bjobs > buffer 2>&1 )
        buffer=$(bjobs 2> /dev/null )
	i=0
	R=0
	P=0
	while read line ; 
	do
	    if [ $i -ne 0 ] ; then
		arr=($line)
		var=${arr[2]}
		#echo "dtag "$dtag
		if [ "$var" == "RUN" ] ; then
			((R++))
		elif [ "$var" == "PEND" ] ; then
			((P++))
		fi
		
	    fi
	    ((i++))
	done <<< "$buffer" #so-called "here variable"
	
	echo $R"       "$P"     "$(($R+$P))
    done
fi


