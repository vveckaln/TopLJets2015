#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/CorrectionTools.h"
#include "TopLJets2015/TopAnalysis/interface/TOPJetShape.h"

#include "TopLJets2015/TopAnalysis/interface/TOPJetPull.h"
#include "TopLJets2015/TopAnalysis/interface/LeptonEfficiencyWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/BtagUncertaintyComputer.h"

#include <vector>
#include <set>
#include <iostream>
#include <algorithm>

#include "TMath.h"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

#include "fastjet/tools/Recluster.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
#include "TopLJets2015/TopAnalysis/interface/CFAT_cmssw.hh"

#include "TopLJets2015/TopAnalysis/interface/JetConstituentAnalysisTool.hh"

#include "CFAT_Event.hh"
#include "TopLJets2015/TopAnalysis/interface/CFAT_AssignHistograms.hh"
#include "TopLJets2015/TopAnalysis/interface/CFAT_Core_cmssw.hh"
#include "TopLJets2015/TopAnalysis/interface/Definitions_cmssw.hh"
#include "TopQuarkAnalysis/TopTools/interface/MEzCalculator.h"

using namespace fastjet;
using namespace fastjet::contrib;

#include "Rivet/Math/MatrixN.hh"
#include "Rivet/Math/MatrixDiag.hh"
using Rivet::Matrix;
using Rivet::EigenSystem;

using namespace std;
//using namespace TOPJetPull_ns;

//
void RunTopJetPull(TString filename,
		    TString outname,
		    Int_t channelSelection, 
		    Int_t chargeSelection, 
		    SelectionTool::FlavourSplitting flavourSplitting,
		    TH1F *normH, 
		    Bool_t runSysts,
		    std::string systVar,
		    TString era,
		    Bool_t debug)
{

  /////////////////////
  // INITIALIZATION //
  ///////////////////
  TH1F h("histo", "histo", 1000, 0, 100);
  TH1F * h_selection_reco = new TH1F("selection_reco", "selection RECO; Selection stage; Events", 2, -0.5, 1.5);
  TH1F * h_selection_gen = new TH1F("selection_gen", "selection GEN; Selection stage; Events", 2, -0.5, 1.5);

  TH1F * const h_selection[] = {h_selection_reco, h_selection_gen};
  const unsigned short Nbin_labels = 4;
  const char* bin_labels[Nbin_labels] = {"1l", "1l + #geq4j", "1l + #geq4j(2b)", "1l + #geq4j(2b, 2lj)"}; 

  for (unsigned short h_ind = 0; h_ind < 2; h_ind++)
    {
      h_selection[h_ind] -> SetDirectory(0);
      for (unsigned short bin_ind = 1; bin_ind <= Nbin_labels; bin_ind ++)
	{
	  h_selection[h_ind] -> GetXaxis() -> SetBinLabel(bin_ind, bin_labels[bin_ind -1]); 
	  //	  h_selection[h_ind] -> GetXaxis() -> SetBinLabel(bin_ind, bin_labels[bin_ind]); 

	}

    }
  TRandom* random = new TRandom(0); // random seed for period selection
  std::vector<RunPeriod_t> runPeriods=getRunPeriods(era);

  bool isTTbar( filename.Contains("_TTJets") or (normH and TString(normH->GetTitle()).Contains("_TTJets")));
  bool isData( filename.Contains("Data") );
  
  // explicit systematics
  std::vector<std::string> vSystVar;
  boost::split(vSystVar, systVar, boost::is_any_of("_"));

  //PREPARE OUTPUT
  TopJetShapeEvent_t tjsev;
  const TString baseName=gSystem->BaseName(outname); 
  const TString dirName=gSystem->DirName(outname);
  printf("outname %s, output file name %s\n", outname.Data(), (dirName+"/"+baseName).Data());
  TFile *fOut=TFile::Open(dirName+"/"+baseName, "RECREATE");
  fOut->cd();
  TTree *outT=new TTree("tjsev","tjsev");
  createTopJetShapeEventTree(outT,tjsev);
  outT->SetDirectory(fOut);


  //READ TREE FROM FILE
  MiniEvent_t ev;
  printf("opening file %s\n", filename.Data());

  TFile *f = TFile::Open(filename);
  TH1 *genPU=(TH1 *)f->Get("analysis/putrue");
  TH1 *triggerList=(TH1 *)f->Get("analysis/triggerList");
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true);
  Int_t nentries(t->GetEntriesFast());
  if (debug) nentries = 10000; //restrict number of entries for testing
  t->GetEntry(0);

  cout << "...producing " << outname << " from " << nentries << " events" << endl;

  //lumi
  TH1F *ratevsrunH=0;
  std::map<Int_t,Float_t> lumiMap;
  if( isData )  
    {
      std::pair<std::map<Int_t,Float_t>, TH1F *> result=parseLumiInfo(era);
      lumiMap   = result.first;
      ratevsrunH = result.second;
    }
  
  //PILEUP WEIGHTING
  std::map<TString, std::vector<TGraph *> > puWgtGr;
  if( !isData ) puWgtGr=getPileupWeightsMap(era,genPU);
  
  
  //LEPTON EFFICIENCIES
  LeptonEfficiencyWrapper lepEffH(filename.Contains("Data13TeV"),era);


  //B-TAG CALIBRATION
  BTagSFUtil* myBTagSFUtil = new BTagSFUtil();
  std::map<TString, std::map<BTagEntry::JetFlavor, BTagCalibrationReader *> > btvsfReaders = getBTVcalibrationReadersMap(era, BTagEntry::OP_MEDIUM);

  //dummy calls
  btvsfReaders[runPeriods[0].first][BTagEntry::FLAV_B]->eval_auto_bounds("central", BTagEntry::FLAV_B,   0., 30.);
  btvsfReaders[runPeriods[0].first][BTagEntry::FLAV_UDSG]->eval_auto_bounds("central", BTagEntry::FLAV_UDSG,   0., 30.);

  std::map<BTagEntry::JetFlavor, TGraphAsymmErrors *>    expBtagEffPy8 = readExpectedBtagEff(era);
  TString btagExpPostFix("");
  if(isTTbar) {
    if(filename.Contains("_herwig"))    btagExpPostFix="_herwig";
    if(filename.Contains("_scaleup"))   btagExpPostFix="_scaleup";
    if(filename.Contains("_scaledown")) btagExpPostFix="_scaledown";
  }
  std::map<BTagEntry::JetFlavor, TGraphAsymmErrors *> expBtagEff=readExpectedBtagEff(era,btagExpPostFix);
  
  
  //JET ENERGY UNCERTAINTIES
  std::string jecVar = "Total";
  if (vSystVar[0] == "jec") {
    if (vSystVar.size() != 3) {
      std::cout << "Wrong format of JEC uncertainty, expected jec_Source_up/down. Exiting..." << std::endl;
      return;
    }
    jecVar = vSystVar[1];
  }
  TString jecUncUrl(era+"/Summer16_23Sep2016V4_MC_UncertaintySources_AK4PFchs.txt");
  gSystem->ExpandPathName(jecUncUrl);
  JetCorrectorParameters *jecParam = new JetCorrectorParameters(jecUncUrl.Data(), jecVar);
  JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty( *jecParam );

  MEzCalculator neutrinoPzComputer;
  
  
  //BFRAG UNCERTAINTIES
  std::map<TString, TGraph*> bfrag = getBFragmentationWeights(era);
  std::map<TString, std::map<int, double> > semilepbr = getSemilepBRWeights(era);
  

  //BOOK HISTOGRAMS
  HistTool ht;
  if (isData) ht.setNsyst(0);
  std::map<TString, TH1 *> allPlots;
  std::map<TString, TH2 *> all2dPlots;
  allPlots["puwgtctr"] = new TH1F("puwgtctr","Weight sums",2,0,2);

  std::vector<TString> stageVec = { "0_pre", "1_1l", "2_1l4j", "3_1l4j2b", "4_1l4j2b2w" };
  std::vector<TString> chTags = { "L", "E", "M" };
  for(auto& stage : stageVec) {
    for(auto& channel : chTags) {  
      TString tag(channel+stage+"_");
      
      if(ratevsrunH)
        allPlots[tag+"ratevsrun"] = (TH1 *)ratevsrunH->Clone(tag+"ratevsrun");
      ht.addHist(tag+"nvtx", new TH1F(tag+"nvtx",";Vertex multiplicity;Events",55,0,55));
      ht.addHist(tag+"nleps", new TH1F(tag+"nleps",";Lepton multiplicity;Events",5,-0.5,4.5));
      ht.addHist(tag+"njets", new TH1F(tag+"njets",";Jet multiplicity;Events",15,-0.5,14.5));
      ht.addHist(tag+"nbjets", new TH1F(tag+"nbjets",";b jet multiplicity;Events",5,-0.5,4.5));
      ht.addHist(tag+"nwjets", new TH1F(tag+"nwjets",";W jet multiplicity;Events",10,-0.5,9.5));
      ht.addHist(tag+"wcandm", new TH1F(tag+"wcandm",";W candidate mass;W candidates",60,0,300));
      ht.addHist(tag+"tcandm", new TH1F(tag+"tcandm",";top candidate mass;top candidates",70,50,400));
      ht.addHist(tag+"tcandwcutm", new TH1F(tag+"tcandwcutm",";top candidate mass;top candidates (W cut)",70,50,400));
      for(int i=0; i<2; i++) {
        TString pf(Form("l%d",i));          
        ht.addHist(tag+pf+"pt", new TH1F(tag+pf+"pt",";Lepton p_{t} [GeV];Events",50,0,250));
        ht.addHist(tag+pf+"eta", new TH1F(tag+pf+"eta",";Lepton pseudo-rapidity;Events",50,-2.5,2.5));
      }
      for(int i=0; i<6; i++) {
        TString pf(Form("j%d",i));
        ht.addHist(tag+pf+"pt", new TH1F(tag+pf+"pt",";Jet transverse momentum [GeV];Events",50,0,250));
        ht.addHist(tag+pf+"eta", new TH1F(tag+pf+"eta",";Jet pseudo-rapidity;Events",50,-5,5));
      }
      ht.addHist(tag+"met", new TH1F(tag+"met",";MET [GeV];Events",50,0,250));
    }
  }
  CFAT_cmssw colour_flow_analysis_tool; 
  AssignHistograms(allPlots);
  AssignSpecificHistograms2D(all2dPlots);
  colour_flow_analysis_tool.plots_ptr_ = & allPlots;
  colour_flow_analysis_tool.plots2D_ptr_ = & all2dPlots;
  colour_flow_analysis_tool.SetMigrationOutput(dirName + "/migration_" + baseName);

  for (auto& it : allPlots)   { it.second->Sumw2(); it.second->SetDirectory(0); }
  for (auto& it : all2dPlots) { it.second->Sumw2(); it.second->SetDirectory(0); }

  std::cout << "init done" << std::endl;

  ///////////////////////
  // LOOP OVER EVENTS //
  /////////////////////
  
  //EVENT SELECTION WRAPPER


  SelectionTool selector(filename, false, triggerList);
  unsigned long nGEN_events = 0;  
  unsigned long nRECO_events = 0; 
  unsigned long nBOTH_events = 0;

  for (Int_t iev = 0;iev < nentries%100; iev ++)
    {
      t->GetEntry(iev);
      resetTopJetShapeEvent(tjsev);
      if(iev%10==0) printf ("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries);
      
      //assign randomly a run period
      TString period = assignRunPeriod(runPeriods,random);
      
      //////////////////
      // CORRECTIONS //
      ////////////////

      
      printf("!!!!!! Event number %u !!!!!!\n", iev);
      
      bool GEN_selected = false;
      bool RECO_selected = false;
      /*         for (int k=0; k<ev.nj; k++) 
	{
	  TLorentzVector jp4;
	  jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);

	  int flavor = 0;
	  if (ev.j_btag[k]) {
	    flavor = 5;
	  }

	  
	  }
*/

      
      /*      for (int part_ind = 0; part_ind <ev.npf; part_ind ++)
	{
	  TLorentzVector particle;
	  particle.SetPtEtaPhiM(ev.pf_pt[part_ind], ev.pf_eta[part_ind], ev.pf_phi[part_ind], ev.pf_m[part_ind]);
	  printf("part_ind %u ind %d rapidity %f, phi %f\n", part_ind, ev.pf_j[part_ind], particle.Rapidity(), particle.Phi());
	  }*/
      double csvm = 0.8484;
      if (vSystVar[0] == "csv") {
          if (vSystVar[1] == "heavy") {
              //heavy flavor uncertainty +/-3.5%
              if (vSystVar[2] == "up")   addBTagDecisions(ev, 0.8726, csvm);
              if (vSystVar[2] == "down") addBTagDecisions(ev, 0.8190, csvm);
          }
          if (vSystVar[1] == "light") {
              //light flavor uncertainty +/-10%
              if (vSystVar[2] == "up")   addBTagDecisions(ev, csvm, 0.8557);
              if (vSystVar[2] == "down") addBTagDecisions(ev, csvm, 0.8415);
          }
      }
      else addBTagDecisions(ev, csvm, csvm);
      
      if(!ev.isData) {
        //jec
        if (vSystVar[0] == "jec") {
          applyJetCorrectionUncertainty(ev, jecUnc, jecVar, vSystVar[2]);
        }
        //jer
        if (vSystVar[0] == "jer") {
          smearJetEnergies(ev, vSystVar[1]);
        }
        else smearJetEnergies(ev);
        //b tagging
        if (vSystVar[0] == "btag") {
          if (vSystVar[1] == "heavy") {
            updateBTagDecisions(ev, btvsfReaders[period],expBtagEff,expBtagEffPy8,myBTagSFUtil,vSystVar[2],"central");
          }
          if (vSystVar[1] == "light") {
            updateBTagDecisions(ev, btvsfReaders[period],expBtagEff,expBtagEffPy8,myBTagSFUtil,"central",vSystVar[2]);
          }
        }
        else updateBTagDecisions(ev, btvsfReaders[period],expBtagEff,expBtagEffPy8,myBTagSFUtil);
      }
      
      ///////////////////////////
      // RECO LEVEL SELECTION //
      /////////////////////////
      
      //decide the lepton channel and get selected objects
      TString chTag = selector.flagFinalState(ev);
      if (chTag != "M" and chTag != "E")
	continue;
      std::vector<Particle> &leptons     = selector.getSelLeptons(); 
      printf("leptons size %lu chTag %s\n", leptons.size(), chTag.Data());
      std::vector<Jet>      &jets        = selector.getJets();  
      std::vector<unsigned int> jet_indices = selector.getJetIndices();
      //count b and W candidates
      int sel_nbjets = 0;
      int sel_nwjets = 0;

      vector<TLorentzVector> bJets, lightJets;
 
      //     vector<TLorentzVector> bJets_gen,lightJets_gen;

      vector<unsigned short> bJets_index, lightJets_index;
      unsigned char jet_index = 0;
      //std::vector<Jet> exp_light_jets;
      
      for (auto& jet : jets) 
	{
	  TLorentzVector jp4 = jet.p4();
	  if (jet.flavor() == 5) 
	    {
	      bJets.push_back(jp4);
	      bJets_index.push_back(jet_indices[jet_index]);
	      ++sel_nbjets;
	    }
	  if (jet.flavor() == 1) 
	    {
	      lightJets.push_back(jp4);
	      lightJets_index.push_back(jet_indices[jet_index]);
	      //      exp_light_jets.push_back(jet);
	      ++sel_nwjets;
	    }
	  jet_index ++;
	}

      //printf("bJets.size() %lu, lightJets.size() %lu, jets.size() %lu\n", bJets.size(), lightJets.size(), jets.size());
      //event selected on reco level?
      bool preselected          (true);
      bool singleLepton         ((chTag=="E" or chTag=="M") and
                                 (selector.getVetoLeptons().size() == 0));
      bool singleLepton4Jets    (singleLepton and jets.size()>=4);
      bool singleLepton4Jets2b  (singleLepton4Jets and sel_nbjets==2);
      bool singleLepton4Jets2b2W(singleLepton4Jets2b and sel_nwjets>0);
      
      std::vector<bool> recoPass; recoPass.push_back(preselected); recoPass.push_back(singleLepton); recoPass.push_back(singleLepton4Jets); recoPass.push_back(singleLepton4Jets2b); recoPass.push_back(singleLepton4Jets2b2W); 
      
      if (singleLepton4Jets2b2W) tjsev.reco_sel = 1;
      
      tjsev.nj=jets.size();
      
      ////////////////////
      // EVENT WEIGHTS //
      //////////////////
      
      float wgt(1.0);
      // Pairs for systematic uncertainty weights
      // double: weight value (divided by central weight)
      // bool: use weight for plotting, otherwise just save to tree
      std::vector<std::pair<double, bool> > varweights;
      std::vector<double> plotwgts;
      
      allPlots["puwgtctr"]->Fill(0.,1.0);
      if (!ev.isData) {
        // norm weight
        wgt  = (normH? normH->GetBinContent(1) : 1.0);
        
        // pu weight
        double puWgt(puWgtGr[period][0]->Eval(ev.g_pu));
        allPlots["puwgtctr"]->Fill(1,puWgt);
        wgt *= puWgt;
        varweights.push_back(std::make_pair(puWgtGr[period][1]->Eval(ev.g_pu), true));
        varweights.push_back(std::make_pair(puWgtGr[period][2]->Eval(ev.g_pu), true));
        
        // lepton trigger*selection weights
        if (singleLepton) {
          EffCorrection_t trigSF = lepEffH.getTriggerCorrection(leptons, period);
          varweights.push_back(std::make_pair(1.+trigSF.second, true));
          varweights.push_back(std::make_pair(1.-trigSF.second, true));
          EffCorrection_t selSF = lepEffH.getOfflineCorrection(leptons[0], period);
          varweights.push_back(std::make_pair(1.+selSF.second, true));
          varweights.push_back(std::make_pair(1.-selSF.second, true));
          wgt *= trigSF.first*selSF.first;
        }
        else varweights.insert(varweights.end(), 4, std::make_pair(1., true));
        
        // bfrag weights
        varweights.push_back(std::make_pair(computeBFragmentationWeight(ev, bfrag["upFrag"]), true));
        varweights.push_back(std::make_pair(computeBFragmentationWeight(ev, bfrag["downFrag"]), true));
        varweights.push_back(std::make_pair(computeBFragmentationWeight(ev, bfrag["PetersonFrag"]), true));
        // weights for semilep BR
        // simultaneous variation for all hadrons, average over particle and antiparticle
        // divide by mean weight from 100k events to avoid change in cross section
        varweights.push_back(std::make_pair(computeSemilepBRWeight(ev, semilepbr["semilepbrUp"], 0, true)/1.00480, true));
        varweights.push_back(std::make_pair(computeSemilepBRWeight(ev, semilepbr["semilepbrDown"], 0, true)/0.992632, true));
        
        // lhe weights
        wgt *= (ev.g_nw>0 ? ev.g_w[0] : 1.0);
        
        std::set<std::string> scalesForPlotter = {
          "id1002muR1muF2hdampmt272.7225",
          "id1003muR1muF0.5hdampmt272.7225",
          "id1004muR2muF1hdampmt272.7225",
          "id1005muR2muF2hdampmt272.7225",
          "id1007muR0.5muF1hdampmt272.7225",
          "id1009muR0.5muF0.5hdampmt272.7225"
        };
        for (int i = 1; i < ev.g_nw; ++i) {
          bool forPlotter = (normH and scalesForPlotter.count(normH->GetXaxis()->GetBinLabel(i)) != 0);
          varweights.push_back(std::make_pair(ev.g_w[i]/ev.g_w[0], forPlotter));
        }
        
        tjsev.nw = 1 + varweights.size();
        tjsev.weight[0]=wgt;
        for (unsigned int i = 0; i < varweights.size(); ++i) {
          tjsev.weight[i+1] = varweights[i].first;
        }
        plotwgts.push_back(wgt);
        for (auto vw : varweights)
          if (vw.second)
            plotwgts.push_back(vw.first);
      }
      else {
        tjsev.nw=1;
        tjsev.weight[0]=wgt;
        plotwgts.push_back(wgt);
      }
      
      //////////////////////////
      // RECO LEVEL ANALYSIS //
      ////////////////////////
      
      //W and top masses
      std::vector<TLorentzVector> wCands;
      for (unsigned int i = 0; i < jets.size(); i++) {
        for (unsigned int j = i+1; j < jets.size(); j++) {
          if (jets[i].flavor()==5 or jets[j].flavor()==5) continue;
          TLorentzVector wCand = jets[i].p4() + jets[j].p4();
          wCands.push_back(wCand);
        }
      }
      std::vector<TLorentzVector> tCands;
      for (unsigned int i = 0; i < jets.size(); i++) {
        if (jets[i].flavor()!=5) continue;
        for (auto& wCand : wCands) {
          TLorentzVector tCand = jets[i].p4() + wCand;
          tCands.push_back(tCand);
        }
      }
      std::vector<TLorentzVector> tCandsWcut;
      for (unsigned int i = 0; i < jets.size(); i++) {
        if (jets[i].flavor()!=5) continue;
        for (auto& wCand : wCands) {
          if (abs(wCand.M()-80.4) > 15.) continue;
          TLorentzVector tCand = jets[i].p4() + wCand;
          tCandsWcut.push_back(tCand);
        }
      }
      //control histograms
      for(size_t istage=0; istage<stageVec.size(); istage++) { 
        for(auto& channel : chTags) { 
          if (not recoPass[istage]) continue;
          if (channel == "E" and chTag != "E") continue;
          if (channel == "M" and chTag != "M") continue;
          TString tag(channel+stageVec[istage]+"_");
          
          std::map<Int_t,Float_t>::iterator rIt=lumiMap.find(ev.run);
          if(rIt!=lumiMap.end() && ratevsrunH) allPlots[tag+"ratevsrun"]->Fill(std::distance(lumiMap.begin(),rIt),1./rIt->second);
          
          ht.fill(tag+"nvtx", ev.nvtx, plotwgts);
          ht.fill(tag+"nleps", leptons.size(), plotwgts);
          ht.fill(tag+"njets", jets.size(), plotwgts);
          ht.fill(tag+"nbjets", sel_nbjets, plotwgts);
          ht.fill(tag+"nwjets", sel_nwjets, plotwgts);
          for (auto& wCand : wCands) ht.fill(tag+"wcandm", wCand.M(), plotwgts);
          for (auto& tCand : tCands) ht.fill(tag+"tcandm", tCand.M(), plotwgts);
          for (auto& tCand : tCandsWcut) ht.fill(tag+"tcandwcutm", tCand.M(), plotwgts);
          for(unsigned int i=0; i<leptons.size(); i++) {
            if (i>1) break;
            TString pf(Form("l%d",i));          
            ht.fill(tag+pf+"pt", leptons[i].pt(), plotwgts);
            ht.fill(tag+pf+"eta", leptons[i].eta(), plotwgts);
          }
          for(unsigned int i=0; i<jets.size(); i++) {
            if (i>5) break;
            TString pf(Form("j%d",i));
            ht.fill(tag+pf+"pt", jets[i].pt(), plotwgts);
            ht.fill(tag+pf+"eta", jets[i].eta(), plotwgts);
          }
          ht.fill(tag+"met", ev.met_pt[0], plotwgts);
        }
      }

      //fill leptons
      if (leptons.size() == 1)
	{
	  h_selection_reco -> Fill("1l", wgt);
	}
      tjsev.nl=leptons.size();
      int il = 0;
      for(auto& lepton : leptons) 
	{
	  tjsev.l_pt[il]  = lepton.pt();
	  tjsev.l_eta[il] = lepton.eta();
	  tjsev.l_phi[il] = lepton.phi();
	  tjsev.l_m[il]   = lepton.m();
	  tjsev.l_id[il]  = lepton.id();
	  il++;
	}
      
      TLorentzVector lp4 = leptons[0].p4();
      //fill MET
      TLorentzVector met(0.0, 0.0, 0.0, 0.0);
      met.SetPtEtaPhiM(ev.met_pt[0], 0.0, ev.met_phi[0], 0.0);
      met.SetPz(0.0); 
      met.SetE(met.Pt());
      const float mt(computeMT(lp4, met));

      //compute neutrino kinematics
      neutrinoPzComputer.SetMET(met);
      neutrinoPzComputer.SetLepton(lp4);
      
      const float nupz = neutrinoPzComputer.Calculate();
      const TLorentzVector neutrinoHypP4(met.Px(), met.Py(), nupz, TMath::Sqrt(TMath::Power(met.Pt(), 2) + TMath::Power(nupz, 2)));
      tjsev.met_pt=ev.met_pt[0];
      tjsev.met_phi=ev.met_phi[0];
      
      //fill jets (with jet shapes)
      for(int ij=0; ij<(int)jets.size(); ij++) {
        tjsev.j_pt[ij]      = jets[ij].p4().Pt();
        tjsev.j_eta[ij]     = jets[ij].p4().Eta();
        tjsev.j_phi[ij]     = jets[ij].p4().Phi();
        tjsev.j_m[ij]       = jets[ij].p4().M(); 
        tjsev.j_flavor[ij]  = jets[ij].flavor();
        tjsev.j_overlap[ij] = jets[ij].overlap();
        
        if (tjsev.reco_sel != 1) continue;
      }
      
      
      ///////////////////////
      // GENERATOR LEVEL  //
      /////////////////////
     
      if (isTTbar) {
        //////////////////////////
        // GEN LEVEL SELECTION //
        ////////////////////////
       
        
        //decide the lepton channel at particle level
        std::vector<Particle> genVetoLeptons = selector.getGenLeptons(ev,15.,2.5);
        std::vector<Particle> genLeptons     = selector.getGenLeptons(ev,30.,2.1);
	
	if (genLeptons.size() == 1)
	  {
	    h_selection_gen -> Fill("1l", wgt);
	  }
	//  const TLorentzVector lp4 = leptons[0].p4()
	//	const TLorentzVector gen_lp4 = genLeptons[0].p4();

        TString genChTag = selector.flagGenFinalState(ev, genLeptons);
        printf("genLeptons.size() %lu genChTag %s\n", genLeptons.size(), genChTag.Data());
	std::vector<Jet> genJets = selector.getGenJets();
	std::vector<unsigned int> & genJets_indices = selector.getGenJetIndices();
        //count b and W candidates
        int sel_ngbjets = 0;
        int sel_ngwcand = 0;

	vector<TLorentzVector> gen_bJets, gen_lightJets;
 
	//     vector<TLorentzVector> bJets_gen,lightJets_gen;

	vector<unsigned short> gen_bJets_index, gen_lightJets_index;
	unsigned char gen_jet_index = 0;
	if (genJets.size() >= 4)
	  {
	    h_selection_gen -> Fill( "1l + #geq4j", wgt);
	  }

	for (auto& genjet : genJets) 
	  {
	    TLorentzVector gen_jp4 = genjet.p4();
	    if (genjet.flavor() == 5) 
	      {
		gen_bJets.push_back(gen_jp4);
		gen_bJets_index.push_back(genJets_indices[gen_jet_index]);
		++sel_ngbjets;
	      }
	    if (genjet.flavor() == 1) 
	      {
		gen_lightJets.push_back(gen_jp4);
		gen_lightJets_index.push_back(genJets_indices[gen_jet_index]);
		++sel_ngwcand;
	      }
	    gen_jet_index ++;
	  }

       
	//	printf("gen_bJets.size() %lu, gen_lightJets.size() %lu, genJets.size() %lu\n", gen_bJets.size(), gen_lightJets.size(), genJets.size());
	static const unsigned char gtop_size = 25;
	TLorentzVector gen_nu; 
	bool gen_nu_found = false;
	for (unsigned char gtop_ind = 0; gtop_ind < gtop_size; gtop_ind ++)
	  {
	    if (abs(ev.gtop_id[gtop_ind]) == 12000 or abs(ev.gtop_id[gtop_ind]) == 14000)
	      {
		gen_nu.SetPtEtaPhiM(ev.gtop_pt[gtop_ind], ev.gtop_eta[gtop_ind], ev.gtop_phi[gtop_ind], ev.gtop_m[gtop_ind]);
		if (gen_nu_found)
		  {
		    gen_nu_found = false;
		    break;
		  }
		gen_nu_found = true;
	      }
	    
	  }
	colour_flow_analysis_tool.ResetMigrationValues();
	if (genJets.size() >= 4 and gen_bJets.size() == 2)
	  {
	    h_selection_gen -> Fill( "1l + #geq4j(2b)", wgt);
	  }
	if (gen_lightJets.size() == 2 and gen_bJets.size() == 2 and gen_nu_found)
	  {
	    h_selection_gen -> Fill( "1l + #geq4j(2b, 2lj)", wgt);
	    GEN_selected = true;
	    //printf("Running GEN\n");
	    CFAT_Core_cmssw core_gen;
	    CFAT_Event event_gen;
	    core_gen.SetEvent(ev);
	    event_gen.SetCore(core_gen);
 
	    core_gen.AddLightJets(gen_lightJets, gen_lightJets_index);
	    core_gen.AddVector(Definitions::LEPTON, lp4);
	    core_gen.AddVector(Definitions::NEUTRINO, gen_nu);
	    core_gen.AddBJets(gen_bJets, gen_bJets_index);
	    core_gen.SetEventDisplayMode(0);
	    event_gen.CompleteVectors();
	    event_gen.SetWeight(wgt);
	    event_gen.SetEventNumber(iev);


	    colour_flow_analysis_tool.SetEvent(event_gen);
	    //	    printf("cfat %p, core_gen %p, event_gen %p\n", & colour_flow_analysis_tool, &event_gen, &core_gen);
	    //printf("CHECK cfat %p, core_gen %p, event_gen %p\n", event_gen.GetCFAT(), core_gen.GetCFATEvent(), event_gen.GetCore());

	    colour_flow_analysis_tool.SetWorkMode(Definitions::GEN);
	   
	    colour_flow_analysis_tool.Work();
	    
	  }
    
      //event selected on gen level?
        bool genSingleLepton((genChTag=="E" or genChTag=="M") and
                             (genVetoLeptons.size() == 1)); // only selected lepton in veto collection
        if (sel_ngbjets==2 && sel_ngwcand>0 && genSingleLepton) tjsev.gen_sel = 1;
        
        tjsev.ngj = genJets.size();
            
        /////////////////////////
        // GEN LEVEL ANALYSIS //
        ///////////////////////

        //store jets to tree
        for (int i = 0; i < tjsev.ngj; i++) {
          tjsev.gj_pt     [i] = genJets[i].p4().Pt();
          tjsev.gj_eta    [i] = genJets[i].p4().Eta();
          tjsev.gj_phi    [i] = genJets[i].p4().Phi();
          tjsev.gj_m      [i] = genJets[i].p4().M();
          tjsev.gj_flavor [i] = genJets[i].flavor();
          tjsev.gj_overlap[i] = genJets[i].overlap();
          
          //matching to reco jet
          for(unsigned int ij = 0; ij< jets.size(); ij++) {
            int ig = i;
            if(jets[ij].p4().DeltaR(genJets[ig].p4())>0.4) continue;
            tjsev.j_gj[ij] = ig;
            tjsev.gj_j[ig] = ij;
            break;
          }
          
          if (tjsev.gen_sel != 1) continue;
          
        }
      }
      else tjsev.gen_sel = -1;
      
      //proceed only if event is selected on gen or reco level
      if (tjsev.gen_sel + tjsev.reco_sel == -2) continue;
      if (jets.size() >= 4)
	{
	  h_selection_reco -> Fill("1l + #geq4j", wgt);
	  
	}


      if (jets.size() >= 4 and bJets.size() == 2)
	{
	  h_selection_reco -> Fill( "1l + #geq4j(2b)", wgt);
	}
      if(bJets.size() == 2 and lightJets.size() == 2)
	{
	  h_selection_reco -> Fill( "1l + #geq4j(2b, 2lj)", wgt);
	  RECO_selected = true;
	  
	  CFAT_Core_cmssw core_reco;
	  CFAT_Event event_reco;
	  core_reco.SetEvent(ev);
	  event_reco.SetCore(core_reco);
	  core_reco.AddLightJets(lightJets, lightJets_index);
 
	  core_reco.AddVector(Definitions::LEPTON, lp4);
	  core_reco.AddVector(Definitions::NEUTRINO, neutrinoHypP4);

	  core_reco.AddBJets(bJets, bJets_index);
	  core_reco.SetEventDisplayMode(0);
	  event_reco.CompleteVectors();
	  event_reco.SetWeight(wgt);
	  event_reco.SetEventNumber(iev);

	  
	  colour_flow_analysis_tool.SetEvent(event_reco);

	  colour_flow_analysis_tool.SetWorkMode(Definitions::RECO);

	  //core_reco.ls_particles(Definitions::SCND_LEADING_JET);
	  //      printf("*** event %u ***** \n", iev);
	  colour_flow_analysis_tool.Work();
	}
      if (isTTbar)
	{
	  //  getchar();

	  colour_flow_analysis_tool.PlotMigrationValues();
	}
      //printf("******************* !!! **********************\n");
      /*if (iev == 41)
	getchar();*/
      if (GEN_selected and not RECO_selected)
	{
	  nGEN_events ++;
	  
	}
      if (not GEN_selected and RECO_selected)
	{
	  nRECO_events ++;
	  
	}
      if (GEN_selected and RECO_selected)
	{
	  nBOTH_events ++;
	}
      //getchar();
   }
  
  //close input file
  f->Close();

  //save histos to file  
  fOut->cd();
  //  outT->Write();
  for (auto& it : ht.getPlots())  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : ht.get2dPlots())  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : allPlots)  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : all2dPlots)  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  if (isTTbar)
    {
      colour_flow_analysis_tool.WriteMigrationTree();
    }

  //  h_selection_gen -> SetDirectory(fOut);
  //h_selection_reco -> SetDirectory(fOut);
  fOut -> cd();
  h_selection_gen -> Write();
  h_selection_reco -> Write();
  fOut->Close();
}
