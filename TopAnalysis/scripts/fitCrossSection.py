#!/usr/bin/env python

import optparse
import os,sys
import commands
from array import array
import ROOT

POItitles={'r':'#mu=#sigma/#sigma_{th}',
           'BtagEff':'kx#sigma_{varepsilon_{b}}',
           'Mtop':'m_{t} [GeV]'}

"""
common CMS label
"""
def drawCMSlabel():        
    cmsLabel=ROOT.TLatex()
    cmsLabel.SetTextFont(42)
    cmsLabel.SetTextSize(0.035)
    cmsLabel.SetNDC()
    cmsLabel.DrawLatex(0.12,0.97,'#bf{CMS} #it{preliminary}')
    cmsLabel.DrawLatex(0.85,0.97,'#scale[0.8]{(13 TeV)}')
    cmsLabel.Draw()

"""
Prepares the fitting script
"""
def prepareFitScript(datacard,POIs,unblind=False):
    baseDir=os.path.dirname(datacard)
    fitScriptUrl='%s/runCombine.sh'%baseDir
    fitScript=open(fitScriptUrl,'w')

    fitScript.write('cd %s\n'%baseDir)
    
    fitScript.write('\n# convert datacard to workspace\n')
    fitScript.write('echo \"Converting datacard to workspace\"\n')
    fitScript.write('text2workspace.py %s -m 0 -o workspace.root\n' % os.path.basename(datacard))

    fitScript.write('\n# likelihood scans\n')
    for parameter in POIs:

        minParVal=0.8 if parameter=='r' else -2.0
        maxParVal=1.2 if parameter=='r' else +2.0
        rangeOpt='--setPhysicsModelParameterRanges %s=%f,%f'%(parameter,minParVal,maxParVal)

        poiOpt='' if parameter=='r' else '--redefineSignalPOIs %s'%parameter

        if parameter=='r':
            fitScript.write('\n## max likelihood fit\n')
            fitScript.write('echo \"Running MaxLikelihoodFit for r\"\n')
            fitScript.write('combine workspace.root -M MaxLikelihoodFit -t -1 --expectSignal=1 -m 0\n')
            fitScript.write('mv mlfit.root mlfit_exp.root\n')
            if unblind:
                fitScript.write('combine workspace.root -M MaxLikelihoodFit -m 0\n')
                fitScript.write('mv mlfit.root mlfit_obs.root\n')
                            
        fitScript.write('\n## function of %s\n'%parameter)
        fitScript.write('echo \"Running likelihood scan for %s\"\n'%parameter)
        fitScript.write('combine workspace.root -M MultiDimFit -P %s -t -1 --expectSignal=1 --algo=grid --points=100 %s %s -m 0\n'%(parameter,rangeOpt,poiOpt))
        fitScript.write('mv higgsCombineTest.MultiDimFit.mH0.root exp_plr_scan_%s.root\n'%parameter)
        fitScript.write('combine workspace.root -M MultiDimFit -P %s -t -1 --expectSignal=1 --algo=grid --points=100 %s %s -m 0 -S 0\n'%(parameter,rangeOpt,poiOpt))
        fitScript.write('mv higgsCombineTest.MultiDimFit.mH0.root exp_plr_scan_stat_%s.root\n'%parameter)
        if unblind:
            fitScript.write('combine workspace.root -M MultiDimFit -P %s --algo=grid --points=200 %s %s -m 0\n'%(parameter,rangeOpt,poiOpt))
            fitScript.write('mv higgsCombineTest.MultiDimFit.mH0.root obs_plr_scan_%s.root\n'%parameter)
            fitScript.write('combine workspace.root -M MultiDimFit -P %s --algo=grid --points=200 %s %s -m 0 -S 0\n'%(parameter,rangeOpt,poiOpt))
            fitScript.write('mv higgsCombineTest.MultiDimFit.mH0.root obs_plr_scan_stat_%s.root\n'%parameter)

    fitScript.write('\n# 2D likelihood scans\n')
    for i in xrange(0,len(POIs)):
        for j in xrange(i+1,len(POIs)):
            minPiVal=0.8 if POIs[i]=='r' else -2.0
            maxPiVal=1.2 if POIs[i]=='r' else +2.0
            minPjVal=0.8 if POIs[j]=='r' else -2.0
            maxPjVal=1.2 if POIs[j]=='r' else +2.0
            rangeOpt='--setPhysicsModelParameterRanges %s=%f,%f:%s=%f,%f'%(POIs[i],minPiVal,maxPiVal,POIs[j],minPjVal,maxPjVal)
            
            fitScript.write('## function of %s,%s\n'%(POIs[i],POIs[j]))
            fitScript.write('echo \"Running 2D likelihood scan for %s vs %s\"\n'%(POIs[i],POIs[j]))
            fitScript.write('combine workspace.root -M MultiDimFit --redefineSignalPOIs %s,%s -P %s -P %s  -t -1 --expectSignal=1 --algo=grid --points=1000 %s -m 0\n'%(POIs[i],POIs[j],POIs[i],POIs[j],rangeOpt)) 
            fitScript.write('mv higgsCombineTest.MultiDimFit.mH0.root  exp_plr_scan_%svs%s.root\n'%(POIs[i],POIs[j]))
            if unblind:
                fitScript.write('combine workspace.root -M MultiDimFit --redefineSignalPOIs %s,%s -P %s -P %s --algo=grid --points=1000  %s -m 0\n'%(POIs[i],POIs[j],POIs[i],POIs[j],rangeOpt)) 
                fitScript.write('mv higgsCombineTest.MultiDimFit.mH0.root obs_plr_scan_%svs%s.root\n'%(POIs[i],POIs[j]))

    fitScript.write('\ncd -\n')
    
    fitScript.close()
    return fitScriptUrl

"""
Opens the output ROOT files from the fits and plots the results for comparison
"""
def show1DLikelihoodScan(resultsSet,parameter='r',output='./'):
   
    #likelihood scans
    nllGrs={}
    colors=[ROOT.kMagenta, ROOT.kAzure+4,  ROOT.kMagenta-9, ROOT.kRed+1, ROOT.kBlue-7]
    ires=0
    for title,datacard in resultsSet:
        ires+=1
        dir=os.path.dirname(datacard)
        files=[('#splitline{expected}{#scale[0.8]{(stat+syst)} }', 'exp_plr_scan',      1, 1),
               ('#splitline{expected}{#scale[0.8]{(stat)}}',       'exp_plr_scan_stat', 3, 1),
               ('#splitline{observed}{#scale[0.8]{(stat+syst)}}',  'obs_plr_scan',      1, 3),
               ('#splitline{observed}{#scale[0.8]{stat only}}',    'obs_plr_scan_stat', 3, 3)]

        for ftitle,f,lstyle,lwidth in files:

            #if file not available continue
            fIn=ROOT.TFile.Open('%s/%s_%s.root' % (dir,f,parameter) )
            if not fIn : continue
        
            #create new graph for the likelihood scan
            if not ftitle in nllGrs: nllGrs[ftitle]=[]
            nllGrs[ftitle].append(ROOT.TGraph())
            nllGrs[ftitle][-1].SetTitle(title)
            nllGrs[ftitle][-1].SetLineStyle(lstyle)
            nllGrs[ftitle][-1].SetMarkerStyle(1)
            nllGrs[ftitle][-1].SetFillStyle(0)
            nllGrs[ftitle][-1].SetLineWidth(lwidth)
            nllGrs[ftitle][-1].SetLineColor(colors[ires-1])
            nllGrs[ftitle][-1].SetMarkerColor(colors[ires-1])

            #fill graph
            tree=fIn.Get('limit')
            for n in xrange(0,tree.GetEntriesFast()):
                tree.GetEntry(n)
                nll=tree.deltaNLL
                if nll>20 : continue
                npoint=nllGrs[ftitle][-1].GetN()
                parVal=getattr(tree,parameter)
                if parameter=='Mtop':
                    parVal=parVal*3+172.5
                nllGrs[ftitle][-1].SetPoint(npoint,parVal,nll)
            nllGrs[ftitle][-1].Sort()
            fIn.Close()

    #show 1D likelihood scan
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.1)
    c.SetRightMargin(0.05)
    minParVal=0.8 if parameter=='r' else -2.0
    maxParVal=1.2 if parameter=='r' else +2.0
    if parameter=='Mtop' : minParVal,maxParVal=166.5,178.5
    frame=ROOT.TH1F('frame',';%s;-2#DeltalnL'%POItitles[parameter],100,minParVal,maxParVal)
    frame.Draw()
    frame.GetYaxis().SetRangeUser(0,10)
    allLegs=[]
    ileg=0
    for ftitle in nllGrs:
        allLegs.append( ROOT.TLegend(0.15+ileg*0.15,0.92,0.3+ileg*0.15,0.85-0.04*len(nllGrs[ftitle]) ) )
        allLegs[-1].SetTextFont(42)
        allLegs[-1].SetTextSize(0.035)
        allLegs[-1].SetBorderSize(0)
        allLegs[-1].SetFillStyle(1001)
        allLegs[-1].SetFillColor(0)
        allLegs[-1].SetHeader(ftitle)
        for gr in nllGrs[ftitle]:
            gr.Draw('c')
            allLegs[-1].AddEntry(gr,gr.GetTitle(),'l')
        ileg+=1
    for leg in allLegs: leg.Draw()

    cl=ROOT.TLine()
    cl.SetLineStyle(2)
    cl.SetLineColor(ROOT.kGray+2)
    txt=ROOT.TLatex()
    txt.SetTextFont(42)
    txt.SetTextSize(0.025)
    txt.SetTextColor(ROOT.kGray+3)
    for delta,title in [(1.0,'68% CL'),(3.84,'95% CL')]:
        txt.DrawLatex(maxParVal-10*frame.GetBinWidth(1),delta+0.25,title)
        cl.DrawLine(minParVal,delta,maxParVal,delta)

    drawCMSlabel()
   
    c.Modified()
    c.Update()
    for ext in ['png','pdf','C']:
        c.SaveAs('%s/nll1dscan_%s.%s'%(output,parameter,ext))
    

"""
2D likelihood scan
"""
def show2DLikelihoodScan(resultsSet,parameters):

    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.1)
    c.SetRightMargin(0.05)
    
    contours = array('d',[1.0,3.84])


    #likelihood scans
    nllGrs={}
    colors=[ROOT.kMagenta, ROOT.kAzure+4,  ROOT.kRed+1,  ROOT.kMagenta-9, ROOT.kBlue-7]
    ires=0
    for title,datacard in resultsSet:
        ires+=1
        dir=os.path.dirname(datacard)
        files=[('#splitline{expected}{#scale[0.8]{(stat+syst)} }', 'exp_plr_scan',      1, 1),
               ('#splitline{expected}{#scale[0.8]{(stat)}}',       'exp_plr_scan_stat', 3, 1),
               ('#splitline{observed}{#scale[0.8]{(stat+syst)}}',  'obs_plr_scan',      1, 3),
               ('#splitline{observed}{#scale[0.8]{stat only}}',    'obs_plr_scan_stat', 3, 3)]

        for ftitle,f,lstyle,lwidth in files:

            if not ftitle in nllGrs: nllGrs[ftitle]=[]

            #if file not available continue
            fIn=ROOT.TFile.Open('%s/%s_%svs%s.root' % (dir,f,parameters[0],parameters[1]) )
            if not fIn : continue
            tree=fIn.Get('limit')

            for ll,ul,tag in [(1-0.99,1.0,'99cl')] : #(1-0.68,1.0,'68cl'),(1-0.95,1.0,'95cl'),(1-0.99,1.0,'99cl')]:
                c.Clear()    
                nllGrs[ftitle].append( ROOT.ll2dContourPlot(tree,parameters[0],parameters[1],ll,ul) )
                nllGrs[ftitle][-1].SetTitle(title)
                nllGrs[ftitle][-1].SetLineStyle(lstyle)
                nllGrs[ftitle][-1].SetMarkerStyle(1)
                nllGrs[ftitle][-1].SetFillStyle(1001)
                nllGrs[ftitle][-1].SetFillColor(colors[ires-1])
                nllGrs[ftitle][-1].SetLineWidth(lwidth)
                nllGrs[ftitle][-1].SetLineColor(colors[ires-1])
                nllGrs[ftitle][-1].SetMarkerColor(colors[ires-1])
                            
    #show 2D likelihood scan
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.1)
    c.SetRightMargin(0.05)
    frame=ROOT.TH2F('frame','frame',10,0.8,1.2,10,-3,3)
    frame.Draw()
    allLegs=[]
    ileg=0
    for ftitle in nllGrs:
        if len(nllGrs[ftitle])==0 : continue
        allLegs.append( ROOT.TLegend(0.15+ileg*0.15,0.92,0.3+ileg*0.15,0.85-0.04*len(nllGrs[ftitle]) ) )
        allLegs[-1].SetTextFont(42)
        allLegs[-1].SetTextSize(0.035)
        allLegs[-1].SetBorderSize(0)
        allLegs[-1].SetFillStyle(1001)
        allLegs[-1].SetFillColor(0)
        allLegs[-1].SetHeader(ftitle)
        for gr in nllGrs[ftitle]:
            gr.Draw('f')
            allLegs[-1].AddEntry(gr,gr.GetTitle(),'f')
        ileg+=1
    for leg in allLegs: leg.Draw()

    drawCMSlabel()
   
    c.Modified()
    c.Update()
    raw_input()

"""
compare prefit and postfit nuisances
"""
def compareNuisances(resultsSet,output):
   
    colors=[ROOT.kMagenta, ROOT.kAzure+4,  ROOT.kRed+1,  ROOT.kMagenta-9, ROOT.kBlue-7]
    ires=0
    frame=None
    gr1s,gr2s=ROOT.TGraph(),ROOT.TGraph()
    postFitNuisGr={}
    dx=1./(len(resultsSet)+2.)
    for title,datacard in resultsSet:
        ires+=1
        dir=os.path.dirname(datacard)
        inF=ROOT.TFile.Open('%s/mlfit_exp.root'%dir)
        fit_s=inF.Get('fit_s')
        npars=fit_s.floatParsFinal().getSize()-1

        #init frames if not yet available
        if frame is None:
            frame=ROOT.TH1F('frame',';Nuisance parameter;N x #sigma_{pre-fit}',npars,0,npars)
            frame.SetDirectory(0)
            
            gr1s.SetMarkerStyle(1)
            gr1s.SetMarkerColor(ROOT.kGreen-8)
            gr1s.SetLineColor(ROOT.kGreen-8)
            gr1s.SetFillStyle(1001)
            gr1s.SetFillColor(ROOT.kGreen-8)
            gr1s.SetPoint(0,0,-1)
            gr1s.SetPoint(1,npars,-1)
            gr1s.SetPoint(2,npars,1)
            gr1s.SetPoint(3,0,1)
            gr1s.SetPoint(4,0,-1)
            
            gr2s.SetMarkerStyle(1)
            gr2s.SetMarkerColor(ROOT.kYellow-10)
            gr2s.SetLineColor(ROOT.kYellow-10)
            gr2s.SetFillStyle(1001)
            gr2s.SetFillColor(ROOT.kYellow-10)
            gr2s.SetPoint(0,0,-2)
            gr2s.SetPoint(1,npars,-2)
            gr2s.SetPoint(2,npars,2)
            gr2s.SetPoint(3,0,2)
            gr2s.SetPoint(4,0,-2)

        #save post fit parameter values
        postFitNuisGr[title]=ROOT.TGraphErrors()
        postFitNuisGr[title].SetTitle(title)
        postFitNuisGr[title].SetMarkerStyle(19+ires)
        postFitNuisGr[title].SetMarkerColor(colors[ires-1])
        postFitNuisGr[title].SetLineColor(colors[ires-1])
        postFitNuisGr[title].SetFillStyle(0)
        for ipar in range(npars):
            var=fit_s.floatParsFinal().at(ipar)
            pname=var.GetName()
            if pname=='r': continue
            np=postFitNuisGr[title].GetN()
            postFitNuisGr[title].SetPoint(np,ipar+0.2+ires*dx,var.getVal())
            postFitNuisGr[title].SetPointError(np,0,var.getError())
            if ires==1:
                frame.GetXaxis().SetBinLabel(ipar+1,pname)
        inF.Close()

    #show 1D likelihood scan
    c=ROOT.TCanvas('c','c',1000,500)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.12)
    c.SetRightMargin(0.05)
    frame.Draw()
    frame.GetYaxis().SetRangeUser(-3,3)
    gr2s.Draw('f')
    gr1s.Draw('f')
    leg=ROOT.TLegend(0.15,0.92,0.6,0.85)
    leg.SetNColumns(len(resultsSet))
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    leg.SetBorderSize(0)
    leg.SetFillStyle(1001)
    leg.SetFillColor(0)
    for ftitle in postFitNuisGr:
        postFitNuisGr[ftitle].Draw('p')
        leg.AddEntry(postFitNuisGr[ftitle],postFitNuisGr[ftitle].GetTitle(),'p')
    leg.Draw()

    txt=ROOT.TLatex()
    txt.SetTextFont(42)
    txt.SetTextSize(0.025)
    txt.SetTextColor(ROOT.kGray+3)
    for delta,title in [(1.0,'1#sigma'),(2,'2#sigma')]:
        txt.DrawLatex(frame.GetXaxis().GetXmax()-0.5,delta+0.2,title)      

    drawCMSlabel()
   
    c.Modified()
    c.Update()
    for ext in ['png','pdf','C']:
        c.SaveAs('%s/nuisances.%s'%(output,ext))

        
"""
main
"""
def main():
  
    #configuration
    usage = 'usage: %prog category1=datacard1 category2=datacard2 ....'
    parser = optparse.OptionParser(usage)
    parser.add_option('-o', '--output',       dest='output',       help='output directory',       default='./',           type='string')
    parser.add_option(      '--noFit',        dest='noFit',        help='don\'t run the fits',    action='store_true')
    parser.add_option(      '--POIs',         dest='POIs',         help='parameters of interest', default='r,BtagEff',       type='string')
    parser.add_option(      '--unblind',      dest='unblind',      help='unblind',                action='store_true')
    (opt, args) = parser.parse_args()

    POIs=opt.POIs.split(',')

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True) #False)
       
    resultsSet=[]
    for newCat in args:
        cat,datacard=newCat.split('=')
        resultsSet.append( (cat,datacard) )
        if not opt.noFit:
            fitScriptUrl=prepareFitScript(datacard=datacard,POIs=POIs,unblind=opt.unblind)
            print 'Fit script for %s available at %s'%(cat,fitScriptUrl)
            os.system('sh %s'%fitScriptUrl)
            
    compareNuisances(resultsSet=resultsSet,output=opt.output)
    
    for parameter in POIs:
        show1DLikelihoodScan(resultsSet=resultsSet,parameter=parameter,output=opt.output)


    #ROOT.gROOT.LoadMacro('src/RootTools.cc+')

    #for i in xrange(0,len(POIs)):
    #    for j in xrange(i+1,len(POIs)):
    #        show2DLikelihoodScan(resultsSet,parameters=[POIs[i],POIs[j]])
    
            

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
