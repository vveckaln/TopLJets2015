#! /usr/bin/env python
import sys, os.path
import ROOT

filename = sys.argv[1]
print >> sys.stderr, "checking sanity of file %s" % filename
if os.path.isfile(filename):
    fIn = ROOT.TFile.Open(filename)

    goodFile = False
    if fIn and not fIn.IsZombie() and not fIn.TestBit(ROOT.TFile.kRecovered):
        goodFile = True
        print >> sys.stderr, "file is good\nSanity test done" 
        fIn.Close()
    else:
        print >> sys.stderr, "problem with file\nSanity test done"
    if not goodFile:
        sys.exit(100)
else:
    print >> sys.stderr, "file  not found\nSanity test done" % filename
    sys.exit(100)

