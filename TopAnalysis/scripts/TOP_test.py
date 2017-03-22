import os
import sys
import optparse
import ROOT
import pickle
import json
from TopLJets2015.TopAnalysis.storeTools import *

"""
Wrapper to be used when run in parallel
"""
    #configuration
usage = 'usage: %prog [options]'

    #compile macro
ROOT.AutoLibraryLoader.enable()
ROOT.gSystem.Load('libTopLJets2015TopAnalysis.so')
ROOT.gROOT.LoadMacro('src/TOP-16-006.cc+')
from ROOT import RunTop16006

cache = 'data/genweights.pck'
    #read normalization
cachefile = open(cache, 'r')
genWgts   = pickle.load(cachefile)
cachefile.close()        
print 'Normalization read from cache (%s)' % cache

    #process tasks
#inF="root://eoscms//eos/cms//store/cmst3/user/psilva/LJets2015/8c1e7c9/MC13TeV_SingleT_t/MergedMiniEvents_0.root"
#inF="root://eoscms//eos/cms//store/cmst3/user/psilva/LJets2015/8c1e7c9/MC13TeV_TTJets/MergedMiniEvents_106.root"
inF="root://eoscms//eos/cms//store/cmst3/user/psilva/LJets2015/8c1e7c9/MC13TeV_TTJets/MergedMiniEvents_60.root"

outF="test_file.root"
channel = 13
charge = 1
flav = 0
runSysts = True
tag = "MC13TeV_SingleT_t"
wgtH = genWgts[tag]
print str(inF), str(outF), channel, charge, flav, wgtH, runSysts
ROOT.RunTop16006(str(inF), str(outF), channel, charge, flav, wgtH, runSysts)
        
