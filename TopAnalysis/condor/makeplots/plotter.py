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
    parser.add_option(      '--mcUnc',            dest='mcUnc'  ,          help='common MC related uncertainty (e.g. lumi)',        default=0,              type=float)
    parser.add_option(      '--com',              dest='com'  ,            help='center of mass energy',                            default='13 TeV',       type='string')
    parser.add_option('-j', '--json',             dest='json'  ,           help='json with list of files',        default=None,              type='string')
    parser.add_option(      '--systTheorJson',    dest = 'systTheorJson',  help = 'json with list of theoretical systematics',       default = None, type='string')
    parser.add_option(      '--systExpJson',      dest = 'systExpJson',    help='json with list of experimental systematics',                    default = None, type='string')
    parser.add_option(      '--signalJson',       dest='signalJson',       help='signal json list',               default=None,              type='string')
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
    method     = opt.method
    signal     = "MC13TeV_TTJets"
    samplesList= {}
    jsonList = opt.json.split(',')
    for jsonPath in jsonList:
        jsonFile = open(jsonPath, 'r')
        samplesList.update(json.load(jsonFile))#, encoding='utf-8', object_pairs_hook=OrderedDict).items()
        jsonFile.close()
    
    #read lists of syst samples
    systTheorSamplesList = {}
    if opt.systTheorJson:
        systTheorJsonList = opt.systTheorJson.split(',')
        for jsonPath in systTheorJsonList:
            jsonFile = open(jsonPath, 'r')
            systTheorSamplesList.update(json.load(jsonFile))#, encoding = 'utf-8').items())
            jsonFile.close()

    systExpSamplesList = {}
    if opt.systExpJson:
        systExpJsonList = opt.systExpJson.split(',')
        for jsonPath in systExpJsonList:
            jsonFile = open(jsonPath, 'r')
            systExpSamplesList.update(json.load(jsonFile))#, encoding = 'utf-8').items()
            jsonFile.close()

    #read list of signal samples
    signalSamplesList = {}
    try:
        jsonFile = open(opt.signalJson, 'r')
        signalSamplesList.update(json.load(jsonFile))#, encoding='utf-8', object_pairs_hook=OrderedDict).items()
        jsonFile.close()
    except:
        pass
    overlaySamplesList = {}

    

    # ind = 0
    # print signal
    # print systTheorSamplesList[signal]
    # print "---------"
    # print systTheorSamplesList[signal]

    for key in list(systTheorSamplesList[signal]):
        tag = systTheorSamplesList[signal][key]
        dellist = ["MC13TeV_TTJets2l2nu_noSC", "MC13TeV_TTJets_m166v5", "MC13TeV_TTJets_m169v5", "MC13TeV_TTJets_m175v5", "MC13TeV_TTJets_m178v5", "MC13TeV_TTJets_widthx0p2", "MC13TeV_TTJets_widthx0p5", "MC13TeV_TTJets_widthx4", "MC13TeV_TTJets_widthx8", "MC13TeV_TTTT"]
        if method != "amcatnlo":
            dellist.append("MC13TeV_TTJets2l2nu_amcatnlo")
        if tag in dellist:
            del systTheorSamplesList[signal][key]
        elif "ext" in systTheorSamplesList[signal][key]:
            del systTheorSamplesList[signal][key]
        # elif not "t#bar{t}" in systTheorSamplesList[signal][ind][1][3]:
        #     del systTheorSamplesList[signal][ind]
        # else:
        #     ind += 1
    for key in list(samplesList):
        if "ext" in key:
            del samplesList[key]
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
        TTJets = samplesList["MC13TeV_TTJets"]
        overlaySamplesList.update({"MC13TeV_TTJets": samplesList["MC13TeV_TTJets"]})
        del samplesList["MC13TeV_TTJets"]
        deltag = "MC13TeV_TTJets" + methodtag
        TTJets_method = systTheorSamplesList[signal][deltag]
        samplesList.update({deltag: systTheorSamplesList[signal][deltag]})
        del systTheorSamplesList[signal][deltag]
        for key in list(systExpSamplesList[signal]):
            newkey = key + "methodtag"
            systExpSamplesList[newkey] = systExpSamplesList[key]
            del systExpSamplesList[key]
    else:
        signalTitle   = "t#bar{t}"
        methodtag     = "nominal"
        deltag        = "MC13TeV_TTJets_cflip"
        TTJets_cflip  = systTheorSamplesList[signal][deltag]
        print TTJets_cflip
        overlaySamplesList.update({deltag: systTheorSamplesList[signal][deltag]})
        for key in list(systTheorSamplesList[signal]):
            if key == deltag:
                del systTheorSamplesList[signal][key]
            elif "width" in key or "ext" in key:
                del systTheorSamplesList[signal][key]
    skipList = opt.skip.split(',')

    #lumi specifications per tag
    lumiSpecs = {}
    if opt.lumiSpecs:
        for spec in opt.lumiSpecs.split(','):
            tag, lumi = spec.split(':')
            lumiSpecs[tag] = float(lumi)

    #proc SF
    procSF = {}
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
    listfile = ROOT.TFile.Open("root://eosuser.cern.ch//" + opt.listfile)
    if not listfile:
        print >> sys.stderr, "failed to open listfile %s" % opt.listfile
        raise Exception('no listfile')
    keylist = listfile.GetListOfKeys()

    inputfiles = []
    inputfiles.append(listfile)
    print overlaySamplesList
    # try:
    #     raw_input("penter")
    # except (EOFError):
    #     pass
    greatlist = [(samplesList, False, False, False, False, ""), (signalSamplesList, True, False, False, False, "")]
    for key in systTheorSamplesList:
        greatlist += [(systTheorSamplesList[key], False, True, False, False, key)]
    for key in systExpSamplesList:
        greatlist += [(systExpSamplesList[key], False, True, False, False, key)]
    # for key in overlaySamplesList:
    #     greatlist += ((overlaySamplesList[key], False, False, False, True, ""))    
    greatlist += [(overlaySamplesList, False, False, False, True, "")]
    for slist, isSignal, isTheorSyst, isExpSyst, isOverlay, sampleref in greatlist:
        if len(slist) == 0: 
            continue
        TESTPLOT = Plot("a", "b", "c", True)
        TESTPLOT.Dot()

        isSyst = isTheorSyst or isExpSyst
        for tag in slist:
            #print "tag ", tag
                # if isSyst and not 't#bar{t}' in sample[3] : continue
            if tag in skipList:
                print("SKIPPED " + tag)
                continue
            sample              = slist[tag]
            xsec                = sample[0]
            isData              = sample[1]
            doFlavourSplitting  = sample[7]
            subProcs            = [(tag, sample[4], sample[5])]
            if doFlavourSplitting:
                subProcs = []
                for flav in [(1, sample[4] + '+l'), (4, sample[4] + '+c'), (5, sample[4] + '+b', sample[5])]:
                    subProcs.append(('%d_%s' % (flav[0], tag), flav[1], sample[5] + 3*len(subProcs)))
            # print "subProces", subProcs 
            for sp in subProcs:
               # raw_input("penter")
                file_name = '%s/%s.root' % ( opt.inDir, sp[0])

                if method != "nominal" and not isExpSyst:
                    file_name = '%s/%s.root' % ( opt.inDir.replace("MC13TeV_TTJets" + methodtag, "MC13TeV_TTJets"), sp[0])
                   # print "file name %s" % file_name
                   # raw_input("penter")
                if method == "nominal" and isOverlay:
                    file_name = '%s/%s.root' % ( opt.inDir.replace("MC13TeV_TTJets" + methodtag, "MC13TeV_TTJets_cflip"), sp[0])

                fIn = None
                if os.path.abspath(file_name) != os.path.abspath(opt.listfile):
                    fIn = ROOT.TFile.Open("root://eosuser.cern.ch/" + file_name )
                    # print "reading %s " % file_name
                    inputfiles.append(fIn)
                else:
                    fIn = listfile
                if not fIn : 
                    print >> sys.stderr, "failed to open %s" % file_name
                    continue

                #fix pileup weighting normalization
                puNormSF=1
                if opt.puNormSF and not isData:
                    puCorrH = fIn.Get(opt.puNormSF)
                    nonWgt  = puCorrH.GetBinContent(1)
                    wgt     = puCorrH.GetBinContent(2)
                    if wgt > 0 :
                        puNormSF = nonWgt/wgt
                        if puNormSF > 1.3 or puNormSF < 0.7 : 
                            puNormSF = 1
                            report += '%s wasn\'t be scaled as too large SF was found (probably low stats)\n' % sp[0]
                        else :
                            report += '%s was scaled by %3.3f for pileup normalization\n' % (sp[0], puNormSF)
                kStart = False
                kEnd   = False
                find   = 0
                kind   = 0
                for tkey in keylist:
                    #keyIsSyst = False
                    try:
                        key = tkey.GetName()
                        # if key == "L_pull_angle_PFNgt20_reco_leading_jet_had_w_DeltaRle1p0":
                        if key == opt.begin:
                            kStart = True
                        if not kStart:
                            continue
                        if kEnd:
                            break
                        if key == opt.end:
                            kEnd = True
                        # print "kind %u key %s" % (kind, key)
                        obj = fIn.Get(key)
                        if obj.InheritsFrom("TH2"):
                            #                            print "inherits from th2"
                            continue
                        kind += 1
                        if not (key == "L_pull_angle_allconst_reco_leading_jet_scnd_leading_jet_DeltaRTotal"):
                            continue
                        print key

                        # print tag 
                        #filter plots using a selection list
                        histos = []
                        if sample[4] == signalTitle:
                            systobj = fIn.Get(key + '_syst')
                            if systobj:
                                #keyIsSyst = True
                                    #key = key[:-5]
                                for ybin in xrange(1, systobj.GetNbinsY() + 1):
                                    # if ybin>0:
                                    #     break
                                    for xbin in xrange(0, systobj.GetNbinsX() + 2):
                                        if math.isnan(systobj.GetBinContent(xbin, ybin)):
                                            systobj.SetBinContent(xbin, ybin, 0)
                                            systobj.SetBinError(xbin, ybin, 0)
                                    weighthist = systobj.ProjectionX('_pxsyst' + str(ybin), ybin, ybin)
                                    weighthist.SetTitle(sp[1] + ' weight ' + str(ybin))
                                    fixExtremities(weighthist, False, False)
                                    weighthist.Draw()
                                    if (weighthist.Integral() > 0): 
                                        histos.append(weighthist)
                        if method != "nominal" and isTheorSyst:
                            objnom = fTTJets.Get(key)
                            obj_method = fTTJets_method.Get(key)
                            for binindex in range(1, obj.GetNbinsX() + 1):
                                obj.SetBinContent(binindex, obj.GetBinContent(binindex) * obj_method.GetBinContent(binindex)/objnom.GetBinContent(binindex))
                        fixExtremities(obj, False, False)
                        histos.append(obj)
                        histos[-1].SetTitle(sp[1])
                        print len(histos)
                        for hist in histos:
                            print hist
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
                            #rebin if n!eeded
                            if opt.rebin > 1:
                                hist.Rebin(opt.rebin)

                            #create new plot if needed
                            print "key in plots", (key in plots)
                            if not key in plots:
                                print "probe"
                                isLogY = False
                                if "chi" in key or "_selection" in key:
                                    isLogY = True
                                #print plots[key], " creating plots[key] ", key
                                plots[key] = Plot(key, signalTitle = signalTitle, com = opt.com, isLogY = isLogY)
                                TEST = Plot(key, signalTitle = signalTitle, com = opt.com, isLogY = isLogY)
                                TEST.Do()
                                print plots[key]
                                plots[key].Do()
                                print "added key"
                            #add process to plot
                            title = hist.GetTitle()
                            keyIsSyst = False
                            if ROOT.TString(hist.GetName()).Contains('_pxsyst'):
                                keyIsSyst = True
                            isSystN = isSyst or keyIsSyst
                            samplesubtitle = ""
                            if not isData and not isSyst:
                                sampelsubtitle = tag
                            if isSyst:
                                samplesubtitle = sampleref
                            print "samplesubtitle [", samplesubtitle, "]"
                            print "adding", key, "  ", title#, " ", plots[key]
                            plots[key].Do()
                            plots[key].addtestt(title = hist.GetTitle())
                            # plots[key].add(h = hist, title = hist.GetTitle(), color = sp[2], isData = sample[1], spImpose = isSignal, isSyst = (isSyst or keyIsSyst), isOverlay = isOverlay, samplesubtititle = samplesubtitle)
                            raw_input("penter")
                            print "probe"
                    except:
                        pass
    #show plots
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)
    outDir = opt.outDir
    #os.system('mkdir -p %s' % outDir)
    #os.system('rm %s/%s' % (outDir, opt.outName))
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

