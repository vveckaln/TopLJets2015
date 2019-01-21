#! /usr/bin/env python
import sys, os.path
import ROOT
print "importing testfilesanity"
def testsanity(filename):
    print >> sys.stderr, "checking sanity of file %s" % filename
    if os.path.isfile(filename):
        fIn = ROOT.TFile.Open(filename)

        goodFile = False
        if fIn and not fIn.IsZombie() and not fIn.TestBit(ROOT.TFile.kRecovered):
            goodFile = True
            print >> sys.stderr, "file is good\nSanity test done" 
            fIn.Close()
            return 0
        else:
            print >> sys.stderr, "problem with file\nSanity test done"
            if not goodFile:
                return 100
            else:
                print >> sys.stderr, "file  not found\nSanity test done" % filename
                return 100
    else:
        print "no file found"
        return 15

if __name__ == '__main__':
    # execute only if run as the entry point into the program
    print "calling testfilesanity from shell"
    filename = sys.argv[1]
    sys.exit(testsanity(filename))
