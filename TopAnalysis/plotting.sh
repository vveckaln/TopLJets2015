#! /bin/sh

./scripts/mergeOutputs.py analysis
#./scripts/mergeOutputs.py analysis_muminus
#./scripts/mergeOutputs.py analysis_munoniso/

python scripts/plotter.py -i analysis/   -j data/era2016/samples.json  -l 12870
#python scripts/plotter.py -i analysis_muplus/   -j data/samples_Run2015.json                           -l 2267.84
#python scripts/plotter.py -i analysis_muminus/  -j data/samples_Run2015.json                           -l 2093.6
#python scripts/plotter.py -i analysis_muplus/   -j data/syst_samples_Run2015.json -o syst_plotter.root -l 2093.6
##python scripts/plotter.py -i analysis_muminus/  -j data/syst_samples_Run2015.json -o syst_plotter.root -l 2093.6
#python scripts/plotter.py -i analysis_munoniso/ -j data/samples_Run2015.json 

#python scripts/runQCDEstimation.py --iso analysis_muplus/plots/plotter.root --noniso analysis_munoniso/plots/plotter.root --out analysis_muplus/
#python scripts/runQCDEstimation.py --iso analysis_muminus/plots/plotter.root --noniso analysis_munoniso/plots/plotter.root --out analysis_muminus/