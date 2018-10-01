#! /usr/bin/env python

import os, sys, ROOT

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i",               dest = "InputDir",   help = "input directory")
parser.add_option("-o",               dest = "OutputDir",  help = "output directory")
parser.add_option(       "--lfile",   dest = "ListSource", help = "file to generate the list of keys")
parser.add_option("-m",               dest = "Method",     help = "method")

(opt, args) = parser.parse_args()

listfile = ROOT.TFile.Open(opt.ListSource)
if not (listfile and not listfile.IsZombie() and not listfile.TestBit(ROOT.TFile.kRecovered)):
    print "bad listfile"
    sys.exit(-1)

keylist = listfile.GetListOfKeys()

next = ROOT.TIter(keylist)
key = next()
ind = 0
begin = []
end = []
while key:
    name = ROOT.TString(key.GetName()) 
    #if not listfile.Get(name).InheritsFrom("TH1"):
     #   continue
    print "%u %s" % (ind, name.Data())
    if ind % 100 == 0:
        begin.append(name)
    if ind % 100 == 99:
        end.append(name)
    key = next()
    ind += 1

if ind - 1 % 100 != 99:
    end.append(name)

ind = 0
print "*"*50
for key in keylist:
    print "ind %u %s" % (ind, key.GetName())
    ind += 1

queue = "longlunch"
PROJECT = "/afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis"
FarmCfgDirectory = os.path.join(PROJECT, "condor", "makeplots", "FARM")
os.system("mkdir -p %s" % FarmCfgDirectory)
CONDOROUTDirectory = os.path.join(PROJECT, "condor", "makeplots", "CONDOROUT")
os.system("mkdir -p %s" % CONDOROUTDirectory)
njob = 0
with open ('%s/condor.sub' % FarmCfgDirectory, 'w') as condor:
    condor.write('executable  = {0}/$(cfgFile).sh\n'.format(FarmCfgDirectory))
    condor.write('output      = {0}/$(cfgFile).out\n'.format(CONDOROUTDirectory))
    condor.write('error       = {0}/$(cfgFile).err\n'.format(CONDOROUTDirectory))
    condor.write('log         = {0}/*.$(ClusterId).log\n'.format(FarmCfgDirectory))
    condor.write('arguments   = $(ClusterId) $(ProcId)\n')
    condor.write('+JobFlavour = "{0}"\n'.format(queue))
    for k in range(0, len(begin)):
        print "%u %s %s" % (k, begin[k].Data(), end[k].Data())
        njob += 1    
        cfgFile = '%s' % begin[k].Data()
        condor.write('cfgFile=%s\n' % cfgFile)
        condor.write('queue 1\n')
        sedcmd  = 'sed \''
        sedcmd += 's%@OUTPUTDIR%'            + opt.OutputDir       + '%g;'
        sedcmd += 's%@INPUTDIR%'             + opt.InputDir        + '%g;'
        sedcmd += 's%@BEGINNAME%'            + begin[k].Data()     + '%g;'
        sedcmd += 's%@ENDNAME%'              + end[k].Data()       + '%g;'
        sedcmd += 's%@PROJECT%'              + PROJECT             + '%g;'
        sedcmd += 's%@METHOD%'               + opt.Method          + '%g;'
        sedcmd += 's%@LISTFILE%'             + opt.ListSource      + '%g;'
        sedcmd += '\''
        os.system('cat ' + PROJECT + '/condor/makeplots/executable_templ.sh | ' + sedcmd + ' > ' + FarmCfgDirectory + '/' + cfgFile + '.sh')
        os.system('chmod u+x %s/%s.sh' % (FarmCfgDirectory, cfgFile))
        # if k == 0:
        #     break
print 'Submitting %u jobs to condor, flavour "%s"' % (njob, queue)
os.system('condor_config_val -dump | grep CONDOR; condor_submit %s/condor.sub' % FarmCfgDirectory)
