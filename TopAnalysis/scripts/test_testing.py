import ROOT
import sys

def main():
        #compile macro
    ROOT.AutoLibraryLoader.enable()
    ROOT.gSystem.Load('libTopLJets2015TopAnalysis.so')
    ROOT.gROOT.LoadMacro('src/myfunction.cc+')
    from ROOT import myfunction

    ROOT.myfunction()

if __name__ == "__main__":
    sys.exit(main())
