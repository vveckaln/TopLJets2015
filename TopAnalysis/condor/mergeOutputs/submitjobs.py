#! /usr/bin/env python
import os, sys
sys.path.append('/afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/condor/')
import testfilesanity

def getSampleNames(dirname):
    names = set()
    for item in os.listdir(dirname):
        names.add(item)
    return names

try:
    InputDir = sys.argv[1]
    OutputDir = InputDir
    if not os.path.isdir(InputDir):
        print "Input directory not found:", InputDir
        exit(-1)
except IndexError:
    print "Need to provide an input directory."
    exit(-1)
queue = "workday"
samplenames = getSampleNames(os.path.join(InputDir, 'Chunks'))
PROJECT = "/afs/cern.ch/work/v/vveckaln/private/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis"
FarmCfgDirectory = os.path.join(PROJECT, "condor", "mergeOutputs", "FARM")
os.system("mkdir -p %s" % FarmCfgDirectory)
CONDOROUTDirectory = os.path.join(PROJECT, "condor", "mergeOutputs", "CONDOROUT")
os.system("mkdir -p %s" % CONDOROUTDirectory)
njob = 0
with open ('%s/condor.sub' % FarmCfgDirectory, 'w') as condor:
    condor.write('executable  = {0}/$(cfgFile).sh\n'.format(FarmCfgDirectory))
    condor.write('output      = {0}/$(cfgFile).out\n'.format(CONDOROUTDirectory))
    condor.write('error       = {0}/$(cfgFile).err\n'.format(CONDOROUTDirectory))
    condor.write('log         = {0}/*.$(ClusterId).log\n'.format(FarmCfgDirectory))
    condor.write('arguments   = $(ClusterId) $(ProcId)\n')
    condor.write('+JobFlavour = "{0}"\n'.format(queue))
    for samplename in samplenames:
        print samplename
        sanity = testfilesanity.testsanity(os.path.join(InputDir, "HADDChunks", samplename +  '.root'))
        migration_sanity = testfilesanity.testsanity(os.path.join(InputDir, "HADDmigration", 'migration_' + samplename +  '.root'))
        if os.path.isfile(os.path.join(InputDir, "HADDChunks", samplename +  '.root')) and os.path.isfile(os.path.join(InputDir, "HADDmigration", 'migration_' + samplename +  '.root')) and sanity == 0 and migration_sanity == 0:
#            print "skipping %s" % samplename
            continue
        else:
#            raw_input ("penter")
            pass
        njob += 1    
        cfgFile = '%s' % samplename
        condor.write('cfgFile=%s\n' % cfgFile)
        condor.write('queue 1\n')
        sedcmd  = 'sed \''
        sedcmd += 's%@OUTPUTDIR%'            + OutputDir  + '%g;'
        sedcmd += 's%@INPUTDIR%'             + InputDir   + '%g;'
        sedcmd += 's%@SAMPLENAME%'           + samplename + '%g;'
        sedcmd += 's%@PROJECT%'              + PROJECT    + '%g;'
        sedcmd += '\''
        os.system('cat ' + PROJECT + '/condor/mergeOutputs/executable_templ.sh | ' + sedcmd + ' > ' + FarmCfgDirectory + '/' + cfgFile + '.sh')
        os.system('chmod u+x %s/%s.sh' % (FarmCfgDirectory, cfgFile))
print 'Submitting %u jobs to condor, flavour "%s"' % (njob, queue)
os.system('condor_config_val -dump | grep CONDOR; condor_submit %s/condor.sub' % FarmCfgDirectory)
        


