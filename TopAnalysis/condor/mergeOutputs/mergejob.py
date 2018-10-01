#! /usr/bin/env python

import sys, os, ROOT
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i", dest = "inputdir",   help = "input directory")
parser.add_option("-o", dest = "outputdir",     help = "output directory")
parser.add_option("-s", dest = "samplename", help = "sample name")
parser.add_option("-T", dest = "notrees",    help = "supply if no trees", action = "store_true", default = False)
(options, args) = parser.parse_args()
inputdir   = options.inputdir
outputdir  = options.outputdir
samplename = options.samplename
notrees    = options.notrees
sampledir = os.path.join(inputdir, samplename)
badfiles = []
goodfiles = ""
for item in os.listdir(sampledir):
    fIn = ROOT.TFile.Open(sampledir + '/' +item)
    if fIn and not fIn.IsZombie() and not fIn.TestBit(ROOT.TFile.kRecovered):
        if goodfiles == "":
            goodfiles = os.path.join(sampledir, item)
        else:
            goodfiles += " " + os.path.join(sampledir, item)
    else:
        badfiles.append(item)
haddtarget = os.path.join(outputdir, "%s.root" % samplename)
if notrees:
    cmd = 'hadd -f -T %s %s' % (haddtarget, goodfiles)
else:
    cmd = 'hadd -f %s %s' % (haddtarget, goodfiles)

os.system(cmd)

if len(badfiles) > 0:
    badfilereport = open("badfilereport.txt", "w")
    for badfile in badfiles:
        badfilereport.write(badfile + '\n')
    badfilereport.close()

