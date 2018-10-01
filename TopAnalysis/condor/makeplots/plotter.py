import optparse
import os,sys
import json
import pickle
from collections import OrderedDict
from TopLJets2015.TopAnalysis.Plot import *

"""
steer the script
"""
def main():
    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option(      '--mcUnc',            dest='mcUnc'  ,      help='common MC related uncertainty (e.g. lumi)',        default=0,              type=float)
    parser.add_option(      '--com',              dest='com'  ,        help='center of mass energy',                            default='13 TeV',       type='string')
    parser.add_option('-j', '--json',             dest='json'  ,      help='json with list of files',        default=None,              type='string')
    parser.add_option(      '--systTheorJson',    dest = 'systTheorJson',     help = 'json with list of theoretical systematics',       default = None, type='string')
    parser.add_option(      '--systExpJson',      dest = 'systExpJson',         help='json with list of experimental systematics',                    default = None, type='string')
    parser.add_option(      '--signalJson',       dest='signalJson',  help='signal json list',               default=None,              type='string')
    parser.add_option(      '--overlayJson',      dest='overlayJson',  help='overlay json list',               default=None,              type='string')
    parser.add_option('-i', '--inDir',            dest='inDir' ,      help='input directory',                default=None,              type='string')
    parser.add_option('-O', '--outDir',           dest='outDir' ,     help='output directory',                default=None,              type='string')
    parser.add_option('-o', '--outName',          dest='outName' ,    help='name of the output file',        default='plotter.root',    type='string')
    parser.add_option(      '--noStack',          dest='noStack',     help='don\'t stack distributions',     default=False,             action='store_true')
    parser.add_option(      '--saveLog',          dest='saveLog' ,    help='save log versions of the plots', default=False,             action='store_true')
    parser.add_option(      '--silent',           dest='silent' ,     help='only dump to ROOT file',         default=False,             action='store_true')
    parser.add_option(      '--onlyData',         dest ='onlyData' ,   help='only plots containing data',     default=False,             action='store_true')
    parser.add_option(      '--saveTeX',          dest ='saveTeX' ,    help='save as tex file as well',       default=False,             action='store_true')
    parser.add_option(      '--rebin',            dest ='rebin',       help='rebin factor',                   default=1,                 type=int)
    parser.add_option('-l', '--lumi',             dest ='lumi' ,       help='lumi to print out',              default=12900,              type=float)
    parser.add_option(      '--lumiSpecs',        dest ='lumiSpecs',   help='lumi specifications for some channels [tag:lumi,tag2:lumi2,...]', default=None,       type=str)
    parser.add_option(      '--only',             dest ='only',        help='plot only these (csv)',          default='',                type='string')
    parser.add_option(      '--skip',             dest ='skip',        help='skip these samples (csv)',          default='',                type='string')
    parser.add_option(      '--puNormSF',         dest ='puNormSF',    help='Use this histogram to correct pu weight normalization', default=None, type='string')
    parser.add_option(      '--procSF',           dest ='procSF',      help='Use this to scale a given process component e.g. "W":.wjetscalefactors.pck,"DY":dyscalefactors.pck', default=None, type='string')
    parser.add_option('-m', '--method',           dest = 'method',    help = 'method - nominal, cflip, amc@nlo',  default = 'nominal',       type = 'string')
    parser.add_option(      '--begin',            dest = 'begin',     help = 'start name')
    parser.add_option(      '--end',              dest = 'end',       help = 'end name')
    parser.add_option(      '--lfile',            dest = 'listfile',   help = 'file wherefrom to generate the list of keys')
    
    (opt, args) = parser.parse_args()
    print "opt.lumi ",  opt.lumi
    method = opt.method
    samplesList=[]
    jsonList = opt.json.split(',')
    for jsonPath in jsonList:
        jsonFile = open(jsonPath,'r')
        samplesList += json.load(jsonFile, encoding='utf-8', object_pairs_hook=OrderedDict).items()
        jsonFile.close()
    
    #read lists of syst samples
    systTheorSamplesList = []
    if opt.systTheorJson:
        systTheorJsonList = opt.systTheorJson.split(',')
        for jsonPath in systTheorJsonList:
            jsonFile = open(jsonPath, 'r')
            systTheorSamplesList += json.load(jsonFile, encoding = 'utf-8').items()
            jsonFile.close()

    systExpSamplesList = []
    if opt.systExpJson:
        systExpJsonList = opt.systExpJson.split(',')
        for jsonPath in systExpJsonList:
            jsonFile = open(jsonPath, 'r')
            systExpSamplesList += json.load(jsonFile, encoding = 'utf-8').items()
            jsonFile.close()

    #read list of signal samples
    signalSamplesList=None
    try:
        jsonFile = open(opt.signalJson,'r')
        signalSamplesList=json.load(jsonFile, encoding='utf-8', object_pairs_hook=OrderedDict).items()
        jsonFile.close()
    except:
        pass
    overlaySamplesList = []
    ind = 0
    while ind < len(systTheorSamplesList):
        tag = systTheorSamplesList[ind][0]
        dellist = ["MC13TeV_TTJets2l2nu_noSC", "MC13TeV_TTJets_m166v5", "MC13TeV_TTJets_m169v5", "MC13TeV_TTJets_m175v5", "MC13TeV_TTJets_m178v5", "MC13TeV_TTJets_widthx0p2", "MC13TeV_TTJets_widthx0p5", "MC13TeV_TTJets_widthx4", "MC13TeV_TTJets_widthx8", "MC13TeV_TTTT"]
        if method != "amcatnlo":
            dellist.append("MC13TeV_TTJets2l2nu_amcatnlo")
        if tag in dellist:
            del systTheorSamplesList[ind]
        elif "ext" in systTheorSamplesList[ind][0]:
            del systTheorSamplesList[ind]
        elif not "t#bar{t}" in systTheorSamplesList[ind][1][3]:
            del systTheorSamplesList[ind]
        else:
            ind += 1
    samplesind = 0
    while samplesind < len(samplesList):
        tag = samplesList[samplesind][0]
        if "ext" in tag:
            del samplesList[samplesind]
        else:
            samplesind += 1
    signalTitle = 0
    methodtag = 0
    if method != 'nominal':
        if method == "amcatnlo":
            signalTitle   = "t#bar{t} aMC@NLO"
            methodtag     = "2l2nu_amcatnlo"
        if method == "cflip":
            signalTitle   = "t#bar{t} cflip"
            methodtag     = "_cflip"
        if method == "herwig":
            signalTitle   = "t#bar{t} Herwig++"
            methodtag     = "_herwig"
        ind_TTJets = 0
        found = False
        while ind_TTJets < len(samplesList) and not found:
            tag, sample = samplesList[ind_TTJets]
            if tag == "MC13TeV_TTJets":
                found = True
            if not found:
                ind_TTJets += 1
        TTJets = samplesList[ind_TTJets]
        del samplesList[ind_TTJets]
        overlaySamplesList.append((TTJets[0], TTJets[1]))
        ind_TTJets_method = 0
        found = False
        while ind_TTJets_method < len(systTheorSamplesList) and not found:
            tag, sample = systTheorSamplesList[ind_TTJets_method]
            if tag == "MC13TeV_TTJets" + methodtag:
                found = True
            if not found:
                ind_TTJets_method += 1
        TTJets_method = systTheorSamplesList[ind_TTJets_method]
        del systTheorSamplesList[ind_TTJets_method]
        samplesList.append((TTJets_method[0], TTJets_method[1]))
        for ind in range(0, len(systExpSamplesList)):
            systExpSamplesList[ind] = (systExpSamplesList[ind][0].replace("TTJets", "TTJets" + methodtag), systExpSamplesList[ind][1])
    else:
        signalTitle   = "t#bar{t}"
        methodtag     = "nominal"
        deltags = ["MC13TeV_TTJets_cflip"]
        ind = 0
        for tag, sample in systTheorSamplesList:
            if tag in deltags:
                del systTheorSamplesList[ind]
            elif "width" in tag or "ext" in tag:
                 del systTheorSamplesList[ind]
            else:
              ind += 1
    skipList=opt.skip.split(',')

    #lumi specifications per tag
    lumiSpecs={}
    if opt.lumiSpecs:
        for spec in opt.lumiSpecs.split(','):
            tag,lumi=spec.split(':')
            lumiSpecs[tag]=float(lumi)

    #proc SF
    procSF={}
    if opt.procSF:
        procList=opt.procSF.split(',')
        for newProc in procList:
            proc,cacheUrl=newProc.split(':')
            if not os.path.isfile(cacheUrl) : continue
            cache=open(cacheUrl,'r')
            procSF[proc]=pickle.load(cache)
            cache.close()
            print 'Scale factors added for',proc

    onlyList = opt.only.split(',')

    #read plots 
    plots = OrderedDict()

    report = ''
    listfile = ROOT.TFile.Open(opt.listfile)
    if not listfile:
        print >> sys.stderr, "failed to open listfile %s" % opt.listfile
        raise Exception('no listfile')
    keylist = listfile.GetListOfKeys()

    inputfiles = []
    inputfiles.append(listfile)
    
    for slist, isSignal, isTheorSyst, isExpSyst, isOverlay in [ (samplesList, False, False, False, False), (signalSamplesList, True, False, False, False), (systTheorSamplesList, False, True, False, False), (systExpSamplesList, False, False, True, False), (overlaySamplesList, False, False, False, True) ]:
        if slist is None: continue
        isSyst = isTheorSyst or isExpSyst
        for tag, sample in slist: 
            if isSyst and not 't#bar{t}' in sample[3] : continue
            if tag in skipList:
              print("SKIPPED "+tag)
              continue
            xsec                = sample[0]
            isData              = sample[1]
            doFlavourSplitting  = sample[6]
            subProcs            = [(tag, sample[3], sample[4])]
            if doFlavourSplitting:
                subProcs = []
                for flav in [(1, sample[3] + '+l'), (4,sample[3]+'+c'),(5,sample[3]+'+b',sample[4])]:
                    subProcs.append(('%d_%s'%(flav[0],tag),flav[1],sample[4]+3*len(subProcs)))
            for sp in subProcs:
#                raw_input("penter")
                file_name = '%s/%s.root' % ( opt.inDir, sp[0])
                if method != "nominal" and not isExpSyst:
                    file_name = '%s/%s.root' % ( opt.inDir.replace("MC13TeV_TTJets" + methodtag, "MC13TeV_TTJets"), sp[0])
  #                  print "file name %s" % file_name
 #                   raw_input("penter")

                fIn = None
                if file_name != opt.listfile:
                    fIn = ROOT.TFile.Open(file_name )
                    inputfiles.append(fIn)
                else:
                    fIn = listfile
                if not fIn : 
                    print >> sys.stderr, "failed to open %s" % file_name
                    continue

                #fix pileup weighting normalization
                puNormSF=1
                if opt.puNormSF and not isData:
                    puCorrH=fIn.Get(opt.puNormSF)
                    nonWgt=puCorrH.GetBinContent(1)
                    wgt=puCorrH.GetBinContent(2)
                    if wgt>0 :
                        puNormSF = nonWgt/wgt
                        if puNormSF > 1.3 or puNormSF < 0.7 : 
                            puNormSF = 1
                            report += '%s wasn\'t be scaled as too large SF was found (probably low stats)\n' % sp[0]
                        else :
                            report += '%s was scaled by %3.3f for pileup normalization\n' % (sp[0],puNormSF)
                kStart = False
                kEnd = False
                find = 0
                kind = 0
                for tkey in keylist:
                    keyIsSyst = False
                    try:
                        key = tkey.GetName()
                        if key == opt.begin:
                            kStart = True
                        if not kStart:
                            continue
                        if kEnd:
                            break
                        if key == opt.end:
                            kEnd = True
#                        print "kind %u key %s" % (kind, key)
                        kind += 1
#                        if not (key == "L_pull_angle_allconst_reco_leading_jet_scnd_leading_jet_DeltaRTotal" ):
#                            continue
                         
#                        print tag
                        #filter plots using a selection list
                        histos = []
                        obj = fIn.Get(key)
                        if obj.InheritsFrom('TH2'):
                            if key[-5:]=='_syst':
                                if sample[3]=='t#bar{t}':
                                    keyIsSyst=True
                                    key = key[:-5]
                                    for ybin in xrange(1,obj.GetNbinsY()+1):
                                        for xbin in xrange(0,obj.GetNbinsX()+2):
                                            if math.isnan(obj.GetBinContent(xbin, ybin)):
                                                obj.SetBinContent(xbin, ybin, 0)
                                                obj.SetBinError(xbin, ybin, 0)
                                        weighthist = obj.ProjectionX('_px'+str(ybin), ybin, ybin)
                                        weighthist.SetTitle(sp[1]+' weight '+str(ybin))
                                        fixExtremities(weighthist, False, False)
                                        weighthist.Draw()
                                        if (weighthist.Integral() > 0): histos.append(weighthist)
                                else:
                                    continue
                            else:
                                histos.append(obj)
                                histos[-1].SetTitle(sp[1])
                        else:
                            if method != "nominal" and isTheorSyst:
                                objnom = fTTJets.Get(key)
                                obj_method = fTTJets_method.Get(key)
                                for binindex in range(1, obj.GetNbinsX() + 1):
                                    obj.SetBinContent(binindex, obj.GetBinContent(binindex) * obj_method.GetBinContent(binindex)/objnom.GetBinContent(binindex))
                            fixExtremities(obj, False, False)
                            histos.append(obj)
                            histos[-1].SetTitle(sp[1])
                        for hist in histos:
                            if not isData and not '(data)' in sp[1]: 

                                #check if a special scale factor needs to be applied
                                sfVal=1.0                            
                                for procToScale in procSF:
                                    if sp[1]==procToScale:
                                        for pcat in procSF[procToScale]:                                    
                                            if pcat not in key: continue
                                            sfVal=procSF[procToScale][pcat][0]
                                            break

                                #scale by lumi
                                lumi = opt.lumi
                                for tag in lumiSpecs:
                                    if not tag in key: 
                                        continue
                                    lumi = lumiSpecs[tag]
                                    break
                                hist.Scale(xsec * lumi * puNormSF * sfVal)                    
                            #rebin if needed
                            if opt.rebin > 1:
                                hist.Rebin(opt.rebin)

                            #create new plot if needed
                            if not key in plots: 
                                plots[key] = Plot(key, signalTitle = signalTitle, com = opt.com)
                            #add process to plot
                            title = hist.GetTitle()
                            isSystN = isSyst or keyIsSyst
                            plots[key].add(h = hist, title = hist.GetTitle(), color = sp[2], isData = sample[1], spImpose = isSignal, isSyst = (isSyst or keyIsSyst), isOverlay = isOverlay)
                    
                    except:
                        pass
    #show plots
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)
    outDir = opt.outDir
    os.system('mkdir -p %s' % outDir)
    os.system('rm %s/%s' % (outDir, opt.outName))
    for p in plots :
        plots[p].mcUnc = opt.mcUnc
        if opt.saveLog: 
            plots[p].savelog=True
        skipPlot = False
        if opt.onlyData and plots[p].dataH is None: 
            skipPlot=True 
        if opt.silent: 
            skipPlot=True
        if not skipPlot : 
            plots[p].show(outDir = outDir, lumi = opt.lumi, noStack = opt.noStack, saveTeX = opt.saveTeX)
        plots[p].appendTo('%s/%s' % (outDir, opt.outName))
        plots[p].reset()
            #raw_input("Press Enter to continue...")
    for inputfile in inputfiles:
        if inputfile:
            inputfile.Close()
    print '-'*50
    print 'Plots and summary ROOT file can be found in %s' % outDir
    if len(report) : print report
    print '-'*50

        
"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())

