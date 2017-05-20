#!/usr/bin/env python
import os,sys
import json
import ROOT
from collections import OrderedDict
jsonPath = 'data/era2016/samples.json'
eos = os.getenv('EOS')
print eos
luminosity = 35922
samplesList = []
jsonFile = open(jsonPath, 'r')
samplesList += json.load(jsonFile, encoding = 'utf-8', object_pairs_hook = OrderedDict).items()
jsonFile.close()
fOut = ROOT.TFile.Open("MC_sum.root", "RECREATE")
histos = dict()
ind = 0
for tag, sample in samplesList:
    xsec   = sample[0]
    isData = sample[1]
    if isData:
        continue
    print "tag %s, xsec %f" % (tag, sample[0])
    fIn = ROOT.TFile.Open(eos + '/analysis/' + tag + '.root')
    if not fIn:
        continue
    print "opened %s" % fIn.GetName()
    check = 0
    for tkey in fIn.GetListOfKeys():
        key = tkey.GetName()
        obj = fIn.Get(key)
        if obj.InheritsFrom("TH2"):
            continue
        check +=1
        #if check > 5:
        #    break
        print "key %s object name %s" % (key, obj.GetName())
        obj.Scale(xsec * luminosity)
        #obj.SetDirectory(0)
        #obj.SetName(obj.GetName() + str(ind))
        if not key in histos:
            histos[key] = obj.Clone()
            print "stored %s"% histos[key].GetName()
            histos[key].SetDirectory(0)
        else:
      #      hist = ROOT.TH1F(histos[key])
            histos[key].Add(obj)
        #obj.Delete()
       #     histos[key] = hist
    ind += 1    
    #fIn.Close()
for key in histos:
    fOut.cd()
    histos[key].Write()
fOut.Close()
