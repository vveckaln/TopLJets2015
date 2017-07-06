#! /bin/sh
INDIR=$EOS/analysis
./scripts/mergeOutputs.py $INDIR
#./scripts/mergeOutputs.py analysis_muminus
#./scripts/mergeOutputs.py analysis_munoniso/
rm $INDIR/plots/plotter.py
python scripts/plotter.py -i $INDIR   -j data/era2016/samples.json  -l 36500
#python scripts/plotter.py -i analysis_muplus/   -j data/samples_Run2015.json                           -l 2267.84
#python scripts/plotter.py -i analysis_muminus/  -j data/samples_Run2015.json                           -l 2093.6
#python scripts/plotter.py -i analysis_muplus/   -j data/syst_samples_Run2015.json -o syst_plotter.root -l 2093.6
##python scripts/plotter.py -i analysis_muminus/  -j data/syst_samples_Run2015.json -o syst_plotter.root -l 2093.6
#python scripts/plotter.py -i analysis_munoniso/ -j data/samples_Run2015.json 

#python scripts/runQCDEstimation.py --iso analysis_muplus/plots/plotter.root --noniso analysis_munoniso/plots/plotter.root --out analysis_muplus/
#python scripts/runQCDEstimation.py --iso analysis_muminus/plots/plotter.root --noniso analysis_munoniso/plots/plotter.root --out analysis_muminus/