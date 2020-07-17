#! /bin/bash
#sh condor/makeplots/haddplotter.sh /eos/user/v/vveckaln/analysis_MC13TeV_TTJets/plots/
dir=$1
target=${dir}/plotter.root
sources=""

for file in ${dir}/plotter*.root;
do
    echo $file
    if [[ $file == $target ]];
    then
	echo "probe"
	continue

    fi
    sources=${sources}" "${file}
    echo "sources="$sources
done

hadd -T -f $target $sources
