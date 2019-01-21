import os
import sys
import optparse
import ROOT
import pickle
import testfilesanity
from collections import OrderedDict
import json
import re
import commands
from TopLJets2015.TopAnalysis.storeTools import *
from TopLJets2015.TopAnalysis.batchTools import *

"""
Wrapper to be used when run in parallel
"""
def RunMethodPacked(args):

    method, inF, outF, channel, charge, flav, runSysts, systVar, era, tag, debug=args
    print 'Running ',method,' on ',inF
    print 'Output file',outF
    print 'Selection ch=',channel,' charge=',charge,' flavSplit=',flav,' systs=',runSysts
    print 'Normalization applied from tag=',tag
    print 'Corrections will be retrieved for era=',era

    try:
        args = '--era %s --normTag %s --in %s --out %s --method %s --charge %d --channel %d --flav %d --systVar %s' % (era, tag, inF, outF, method, charge, channel, flav, systVar)
        if runSysts: 
            args += ' --runSysts'
        cmd = 'analysisWrapper ' + args
        #cmd = 'gdb -ex \'break src/TOPJetPull.cc:254\' -ex \'break src/TOPJetPull.cc:567\' -ex \'break src/TOPJetPull.cc:767\' -ex \'run ' + args + '\' analysisWrapper' 
        if debug : cmd += ' --debug'
        print(cmd)
        exit_code = os.system(cmd)
        print "exit code in RunMethodPacked %s" % str(exit_code)
        os.system('echo "exit code in RunMethodPacked ' + str(exit_code) + '" 1>&2')
        return exit_code

    except :
        print 50*'<'
        print "  Problem  (%s) with %s continuing without"%(sys.exc_info()[1],inF)
        print 50*'<'
        return False
    return True

"""
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-m', '--method',      dest='method',      help='method to run [%default]',                   default='TOP-16-006::RunTop16006',  type='string')
    parser.add_option('-i', '--in',          dest='input',       help='input directory with files or single file [%default]',  default=None,       type='string')
    parser.add_option('-o', '--out',         dest='output',      help='output directory (or file if single file to process)  [%default]',  default='analysis', type='string')
    parser.add_option(      '--only',        dest='only',        help='csv list of samples to process  [%default]',             default=None,       type='string')
    parser.add_option(      '--skip',        dest='skip',        help='csv list of samples to skip  [%default]',             default=None,       type='string')
    parser.add_option(      '--runSysts',    dest='runSysts',    help='run systematics  [%default]',                            default=False,      action='store_true')
    parser.add_option(      '--systVar',     dest='systVar',     help='single systematic variation  [%default]',   default='nominal',       type='string')
    parser.add_option(      '--debug',       dest='debug',      help='debug mode  [%default]',                            default=False,      action='store_true')
    parser.add_option(      '--flav',        dest='flav',        help='split according to heavy flavour content  [%default]',   default=0,          type=int)
    parser.add_option(      '--ch',          dest='channel',     help='channel  [%default]',                                    default=13,         type=int)
    parser.add_option(      '--charge',      dest='charge',      help='charge  [%default]',                                     default=0,          type=int)
    parser.add_option(      '--era',         dest='era',         help='era to use for corrections/uncertainties  [%default]',   default='era2016',       type='string')
    parser.add_option(      '--tag',         dest='tag',         help='normalize from this tag  [%default]',                    default=None,       type='string')
    parser.add_option('-q', '--queue',       dest='queue',       help='if not local send to batch with condor. queues are now called flavours, see http://batchdocs.web.cern.ch/batchdocs/local/submit.html#job-flavours   [%default]',     default='local',    type='string')    
    parser.add_option('-n', '--njobs',       dest='njobs',       help='# jobs to run in parallel  [%default]',                  default=0,    type='int')
    parser.add_option(      '--skipexisting',dest='skipexisting',help='skip jobs with existing output files  [%default]',       default=False,      action='store_true')
    parser.add_option(      '--exactonly',   dest='exactonly',   help='match only exact sample tags to process  [%default]',    default=False,      action='store_true')
    parser.add_option(      '--outputonly',        dest='outputonly',        help='filter job submission for a csv list of output files  [%default]',             default=None,       type='string')
    parser.add_option(      '--farmappendix',        dest='farmappendix',        help='Appendix to condor FARM directory [%default]',             default=None,       type='string')
    (opt, args) = parser.parse_args()
    PROJECT = os.getenv('PROJECT')

    skipexisting = opt.skipexisting
    #parse selection lists
    onlyList = []
    try:
        for t in opt.only.split(','):
            if '.json' in t:
                jsonFile = open(t,'r')
                samplesList = json.load(jsonFile, encoding='utf-8', object_pairs_hook=OrderedDict).items()
                for s,_ in samplesList: onlyList.append(s)
                jsonFile.close()
            else:
                onlyList.append(t)
    except:
        pass
    skipList=[]
    try:
        skipList=opt.skip.split(',')
    except:
        pass
    outputOnlyList=[]
    try:
        outputOnlyList=opt.outputonly.split(',')
    except:
        pass
    #parse list of systematic variations
    varList=[]
    if opt.systVar == 'all':
        allSystVars = ['jec_CorrelationGroupMPFInSitu', 'jec_RelativeFSR',
                       'jec_CorrelationGroupUncorrelated', 'jec_FlavorPureGluon', 'jec_FlavorPureQuark',
                       'jec_FlavorPureCharm', 'jec_FlavorPureBottom', 'jer',
                       'btag_heavy', 'btag_light', 'csv_heavy', 'csv_light', 'tracking']
        for var in allSystVars:
            varList.append(var+'_up')
            varList.append(var+'_down')
    else:
        try:
            varList=opt.systVar.split(',')
        except:
            pass
    print 'Running following variations: ', varList

    #prepare output if a directory
    if not '.root' in opt.output :
        print opt.output
        if not '/store/' in opt.output:
            os.system('mkdir -p %s/Chunks' % opt.output)
        else:
            os.system('eos mkdir %s' % opt.output)
            os.system('eos mkdir %s/Chunks' % opt.output)
    #correct location of corrections to be used using cmsswBase, if needed
    cmsswBase = os.environ['CMSSW_BASE']
    if not cmsswBase in opt.era : opt.era=cmsswBase+'/src/TopLJets2015/TopAnalysis/data/'+opt.era

    #process tasks
    task_list = []
    processedTags=[]
    if '.root' in opt.input:
        inF=opt.input
        if '/store/' in inF and not 'root:' in inF : 
            inF = 'root://eoscms//eos/cms' + opt.input        
        for systVar in varList:
            outF = opt.output
            if systVar != 'nominal' and not systVar in opt.output: 
                outF=opt.output[:-5]+'_'+systVar+'.root'
            task_list.append( (opt.method, inF, outF, opt.channel, opt.charge, opt.flav, opt.runSysts, systVar, opt.era, opt.tag, opt.debug) )
    else:

        inputTags = getEOSlslist(directory = opt.input, prepend='')
        for baseDir in inputTags:

            tag = os.path.basename(baseDir)
            if tag=='backup' : continue

            #filter tags
            if len(onlyList)>0:
                processThisTag=False
                for itag in onlyList:
                    if ((opt.exactonly and itag == tag) or 
                        (not opt.exactonly and itag in tag)):
                        processThisTag=True
                        break
                if not processThisTag : continue
            if len(skipList)>0:
                processThisTag=True
                for itag in skipList:
                    if itag in tag:
                        processThisTag=False
                        break
                if not processThisTag : continue
                
            input_list = getEOSlslist(directory = '%s/%s' % (opt.input, tag) )
                
            for systVar in varList:
                nexisting = 0
                for ifile in xrange(0, len(input_list)):
                    inF = input_list[ifile]

                    outF           = os.path.join(opt.output, 'Chunks', tag, '%s_%d.root' % (tag, ifile))
                    migration_outF = os.path.join(opt.output, 'migration', tag, 'migration_%s_%d.root' %(tag, ifile))
                    if systVar != 'nominal' and not systVar in tag: 
                        outF           = os.path.join(opt.output,'Chunks', tag + '_' + systVar, '%s_%s_%d.root' % (tag, systVar, ifile))
                        migration_outF = os.path.join(opt.output, 'migration', tag + '_' + systVar, 'migration_%s_%s_%d.root' % (tag, systVar, ifile))
#                    print "%s %s %s %s" % (outF, migration_outF, os.path.isfile(outF), os.path.isfile(migration_outF))
#                    raw_input("press Enter")
                    sanity_outF = testfilesanity.testsanity(outF)
                    sanity_migration_outF = testfilesanity.testsanity(migration_outF)
                    if sanity_outF != 0 or sanity_migration_outF != 0:
                        raw_input("penter")
                        pass
                    if opt.skipexisting and os.path.isfile(outF) and os.path.isfile(migration_outF) and sanity_outF == 0 and sanity_migration_outF == 0:
                        nexisting += 1
                        continue

                    if (len(outputOnlyList) > 1 and not outF in outputOnlyList):
                        continue
                    task_list.append( (opt.method, inF, outF, opt.channel, opt.charge, opt.flav, opt.runSysts, systVar, opt.era, tag, opt.debug) )
                if (opt.skipexisting and nexisting): 
                    print '--skipexisting: %s - skipping %d of %d tasks as files already exist' % (systVar, nexisting, len(input_list))

    #run the analysis jobs
                
    if opt.queue=='local':
        print 'launching %d tasks in %d parallel jobs'%(len(task_list),opt.njobs)
        if opt.njobs == 0:
            for args in task_list: 
                catch_error_code = RunMethodPacked(args)
                print "caught error code %u" % catch_error_code
                os.system('echo "caught error code ' + str(catch_error_code) + '" 1>&2')
                sys.exit(catch_error_code >> 8)
        else:
            from multiprocessing import Pool
            pool = Pool(opt.njobs)
            pool.map(RunMethodPacked, task_list)
    else:
        MigrationDirectory = '%s/migration' % opt.output
        FarmOutputDirectory = '%s/FARM%s' % (opt.output, opt.farmappendix)
        FarmCfgDirectory = '%s/condor/FARM' % PROJECT
        os.system('mkdir -p %s' % FarmOutputDirectory)
        os.system('mkdir -p %s' % FarmCfgDirectory)
        os.system('mkdir -p %s' % MigrationDirectory)

        print 'Preparing %d tasks to submit to the batch'%len(task_list)
        print 'Executables and condor wrapper are stored in %s'%FarmCfgDirectory

        with open ('%s/condor.sub' % FarmCfgDirectory, 'w') as condor:

            condor.write('executable  = {0}/$(cfgFile).sh\n'.format(FarmCfgDirectory))
            condor.write('output      = {0}/condor/CONDOROUT/$(cfgFile).out\n'.format(PROJECT))
            condor.write('error       = {0}/condor/CONDOROUT/$(cfgFile).err\n'.format(PROJECT))
            condor.write('log         = {0}/*.$(ClusterId).log\n'.format(FarmCfgDirectory))
            condor.write('arguments   = $(ClusterId) $(ProcId)\n')
            condor.write('+JobFlavour = "{0}"\n'.format(opt.queue))

            jobNb = 0
            for method, inF, outF, channel, charge, flav, runSysts, systVar, era, tag, debug in task_list:
#                if (jobNb > 1500):
#                    break
#                jobNb += 1
#                if (jobNb <1450):
#                    continue

                job_tag = os.path.splitext(os.path.basename(outF))[0]
                outfile_basename_noext = job_tag
                localOutF = os.path.basename(outF)
                isTTbar = False
                isTTbar_string = "false"
                if "TTJets" in localOutF:
                    isTTbar = True
                    isTTbar_string = "true"
                runOpts='--charge %d --ch %d --era %s --tag %s --flav %d --method %s --systVar %s'\
                        %(charge, channel, era, tag, flav, method, systVar)
                if runSysts : runOpts += ' --runSysts'
                if debug :    runOpts += ' --debug'

                cfgFile = '%s' % job_tag

                condor.write('cfgFile=%s\n' % cfgFile)
                condor.write('queue 1\n')
                sedcmd  = 'sed \''
                sedcmd += 's%@OUTPUTDIR%'            + opt.output             + '%g;'
                sedcmd += 's%@INPUT_FILE%'           + inF                    + '%g;'
                sedcmd += 's%@OUTPUT_FILE_BASENAME%' + outfile_basename_noext + '%g;'
                sedcmd += 's%@JOB_TAG%'              + job_tag                + '%g;'
                sedcmd += 's%@PROJECT%'              + PROJECT                + '%g;'
                sedcmd += 's%@RUN_OPTS%'             + runOpts                + '%g;'
                sedcmd += 's%@TAG%'                  + tag                    + '%g;'
                sedcmd += 's%@SYS%'                  + systVar                + '%g;'
                sedcmd += '\''
                
                os.system('cat ' + PROJECT + '/condor/executable.templ | ' + sedcmd + ' > condor/FARM/' + cfgFile + '.sh')
                
                os.system('chmod u+x %s/%s.sh' % (FarmCfgDirectory, cfgFile))

        print 'Submitting jobs to condor, flavour "%s"'%(opt.queue)
        os.system('condor_config_val -dump | grep CONDOR; condor_submit %s/condor.sub' % FarmCfgDirectory)
        

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
