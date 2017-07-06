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
#include "TopLJets2015/TopAnalysis/interface/MttbarAnalyzer.h"
#include "TopLJets2015/TopAnalysis/interface/LeptonEfficiencyWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/TOPJetShape.h"


#include <vector>
#include <set>
#include <iostream>
#include <algorithm>

#include "TMath.h"
#include "TopQuarkAnalysis/TopTools/interface/MEzCalculator.h"


using namespace std;


//
void RunMttbarAnalyzer(TString filename,
                       TString outname,
                       Int_t channelSelection, 
                       Int_t chargeSelection, 
                       TH1F *normH, 
                       TString era,
                       Bool_t debug)
{
  /////////////////////
  // INITIALIZATION //
  ///////////////////
  TRandom* random = new TRandom(0); // random seed for period selection
  std::vector<RunPeriod_t> runPeriods=getRunPeriods(era);
  bool isTTbar( filename.Contains("_TTJets") or (normH and TString(normH->GetTitle()).Contains("_TTJets")));
  bool isData( filename.Contains("Data") );
  
  //PREPARE OUTPUT
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
  fOut->cd();
  TTree *tree=new TTree("data","data");
  MttbarSummary_t summary;
  createMttbarSummaryTree(tree,summary);
  tree->SetDirectory(fOut);

  //READ TREE FROM FILE
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);
  TH1 *genPU=(TH1 *)f->Get("analysis/putrue");
  TH1 *triggerList=(TH1 *)f->Get("analysis/triggerList");
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true);
  Int_t nentries(t->GetEntriesFast());
  if (debug) nentries = 10000; //restrict number of entries for testing
  t->GetEntry(0);

  cout << "...producing " << outname << " from " << nentries << " events" << endl;
  
  //auxiliary to solve neutrino pZ using MET
  MEzCalculator neutrinoPzComputer;

  //PILEUP WEIGHTING
  std::map<TString, std::vector<TGraph *> > puWgtGr;
  if( !isData ) puWgtGr=getPileupWeightsMap(era,genPU);
    
  //LEPTON EFFICIENCIES
  LeptonEfficiencyWrapper lepEffH(filename.Contains("Data13TeV"),era);

  //B-TAG CALIBRATION
  std::map<TString, std::map<BTagEntry::JetFlavor, BTagCalibrationReader *> > btvsfReaders = getBTVcalibrationReadersMap(era, BTagEntry::OP_MEDIUM);
  std::map<BTagEntry::JetFlavor, TGraphAsymmErrors *>    expBtagEffPy8 = readExpectedBtagEff(era);
  
   //BOOK HISTOGRAMS
  HistTool ht;
  ht.setNsyst(0);
  ht.addHist("puwgtctr",     new TH1F("puwgtctr",    ";Weight sums;Events",2,0,2));
  ht.addHist("nvtx",         new TH1F("nvtx",        ";Vertex multiplicity;Events",55,-0.5,49.5));
  ht.addHist("njets",        new TH1F("njets",       ";Jet multiplicity;Events",15,-0.5,14.5));
  ht.addHist("nbjets",       new TH1F("nbjets",      ";b jet multiplicity;Events",10,-0.5,9.5));
  ht.addHist("ht",           new TH1F("ht",          ";H_{T} [GeV];Events",50,0,250));
  ht.addHist("mttbar",   new TH1F("mttbar",  ";M_{ttbar} [GeV];Events",100,300,1000));

  std::cout << "init done" << std::endl;

  ///////////////////////
  // LOOP OVER EVENTS //
  /////////////////////
  
  //EVENT SELECTION WRAPPER
  SelectionTool selector(filename, false, triggerList);
  
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      if(iev%10==0) printf ("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries);
      
      //assign randomly a run period
      TString period = assignRunPeriod(runPeriods,random);
      
      //////////////////
      // CORRECTIONS //
      ////////////////
      double csvm = 0.8484;
      addBTagDecisions(ev, csvm, csvm);
      if(!ev.isData) ev = smearJetEnergies(ev);
           
      ///////////////////////////
      // RECO LEVEL SELECTION //
      /////////////////////////
      TString chTag = selector.flagFinalState(ev);
      if(chTag=="") continue;
      std::vector<Particle> &leptons     = selector.getSelLeptons(); 
      std::vector<Jet>      &jets        = selector.getJets();  

      //count n b-jets
      std::vector<Jet> bJets,lightJets;
      float scalarht(0.);
      for(size_t ij=0; ij<jets.size(); ij++) 
        {
          if(jets[ij].flavor()==5) bJets.push_back(jets[ij]);
          else                     lightJets.push_back(jets[ij]);
          scalarht += jets[ij].pt();
        }

      //require one good lepton
      if(leptons.size()!=1) continue;
      bool passJets(lightJets.size()>=2);
      bool passBJets(bJets.size()>=2);
      if(!passJets || !passBJets) continue;
      
      ////////////////////
      // EVENT WEIGHTS //
      //////////////////
      float wgt(1.0);
      std::vector<double>plotwgts(1,wgt);
      ht.fill("puwgtctr",0,plotwgts);
      if (!ev.isData) {

        // norm weight
        wgt  = (normH? normH->GetBinContent(1) : 1.0);
        
        // pu weight
        double puWgt(puWgtGr[period][0]->Eval(ev.g_pu));
        std::vector<double>puPlotWgts(1,puWgt);
        ht.fill("puwgtctr",1,puPlotWgts);
        
        // lepton trigger*selection weights
        EffCorrection_t trigSF = lepEffH.getTriggerCorrection(leptons, period);
        EffCorrection_t  selSF= lepEffH.getOfflineCorrection(leptons[0], period);

        wgt *= puWgt*trigSF.first*selSF.first;
      
        // generator level weights
        wgt *= (ev.g_nw>0 ? ev.g_w[0] : 1.0);

        //update weight for plotter
        plotwgts[0]=wgt;
      }

      //visible system
      TLorentzVector visSystem(leptons[0].p4()+bJets[0].p4()+bJets[1].p4()+lightJets[0].p4()+lightJets[1].p4());
      
      //determine the neutrino kinematics
      TLorentzVector met(0,0,0,0);
      met.SetPtEtaPhiM(ev.met_pt[0],0,ev.met_phi[0],0.);
      neutrinoPzComputer.SetMET(met);
      neutrinoPzComputer.SetLepton(leptons[0].p4());
      float nupz=neutrinoPzComputer.Calculate();
      TLorentzVector neutrinoP4(met.Px(),met.Py(),nupz ,TMath::Sqrt(TMath::Power(met.Pt(),2)+TMath::Power(nupz,2)));
          
      //ttbar system
      TLorentzVector ttbarSystem(visSystem+neutrinoP4);

      //control histograms
      ht.fill("nvtx",     ev.nvtx,        plotwgts);
      ht.fill("nbjets", bJets.size(), plotwgts);
      ht.fill("njets",  jets.size(),  plotwgts);
      ht.fill("ht",         scalarht, plotwgts);          
      ht.fill("mttbar", ttbarSystem.M(),   plotwgts);

  
      //event weight
      summary.weight=wgt;
      summary.nvtx=ev.nvtx;
      summary.rho=ev.rho;

      summary.mttbar=ttbarSystem.M();

      //save lepton
      summary.l_pt=leptons[0].p4().Pt();
      summary.l_eta=leptons[0].p4().Eta();
      summary.l_phi=leptons[0].p4().Phi();
      summary.l_m=leptons[0].p4().M();

      //save met
      summary.met_pt=ev.met_pt[0];
      summary.met_phi=ev.met_phi[0];

      //save jets
      summary.nj=jets.size();
      summary.nb=bJets.size();
      summary.nl=lightJets.size();
      for(size_t ij=0; ij<jets.size(); ij++)
        {
          summary.j_pt[ij]=jets[ij].p4().Pt();
          summary.j_eta[ij]=jets[ij].p4().Eta();
          summary.j_phi[ij]=jets[ij].p4().Phi();
          summary.j_m[ij]=jets[ij].p4().M();
          summary.j_csv[ij]=jets[ij].getCSV();
          summary.j_nch[ij]=getMult(jets[ij]);
          summary.j_PtD[ij]=getPtD(jets[ij]);
          summary.j_PtDs[ij]=getPtDs(jets[ij]);
          summary.j_width[ij]=getWidth(jets[ij]);
          summary.j_tau21[ij]=getTau(2, 1, jets[ij]);
          summary.j_tau32[ij]=getTau(3, 2, jets[ij]);
          summary.j_tau43[ij]=getTau(4, 3, jets[ij]);
          std::vector<double> zgResult_charged = getZg(jets[ij]);
          summary.j_zg[ij]=zgResult_charged[0];
        }

      //ue summary
      summary.ue_nch=0;
      summary.ue_chsumpt=0;
      summary.ue_chsumpz=0;
      for(int ipf = 0; ipf < ev.npf; ipf++)
        {
          TLorentzVector tkP4(0,0,0,0);
          tkP4.SetPtEtaPhiM(ev.pf_pt[ipf],ev.pf_eta[ipf],ev.pf_phi[ipf],0.);
          if(ev.pf_c[ipf]==0) continue; 
          if( tkP4.Pt()<0.9 || fabs(tkP4.Eta())>2.4) continue;
          summary.ue_nch++;
          summary.ue_chsumpt+=tkP4.Pt();
          summary.ue_chsumpz+=fabs(tkP4.Pz());
        }

      //save ttbar and pseudo-ttbar kinematics
      if(ev.ngtop>0)
        {
          TLorentzVector ttbar(0,0,0,0);
          for(Int_t it=0; it<ev.ngtop; it++)
            {
              int absid(abs(ev.gtop_id[it]));
              if(absid!=6) continue;
              TLorentzVector p4(0,0,0,0);
              p4.SetPtEtaPhiM(ev.gtop_pt[it],ev.gtop_eta[it],ev.gtop_phi[it],ev.gtop_m[it]);
              ttbar += p4;
            }
          summary.gen_mttbar=ttbar.M();
        }

      tree->Fill();
    }
  
  //close input file
  f->Close();
  
  //save histos to file  
  fOut->cd();
  tree->Write();
  for (auto& it : ht.getPlots())  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : ht.get2dPlots())  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  fOut->Close();
}


//
void createMttbarSummaryTree(TTree *t,MttbarSummary_t &summary)
{
  //event category
  t->Branch("weight",      &summary.weight,      "weight/F");
  t->Branch("nj",          &summary.nj,          "nj/I");
  t->Branch("nb",          &summary.nb,          "nb/I");
  t->Branch("nl",          &summary.nl,          "nl/I");
  t->Branch("nvtx",        &summary.nvtx,        "nvtx/I");
  t->Branch("rho",         &summary.rho,         "rho/F");
  t->Branch("gen_mttbar",         &summary.gen_mttbar,      "gen_mttbar/F");
  t->Branch("mttbar",              &summary.mttbar,           "mttbar/F");
  t->Branch("l_pt",              &summary.l_pt,           "l_pt/F");
  t->Branch("l_eta",              &summary.l_eta,           "l_eta/F");
  t->Branch("l_phi",              &summary.l_phi,           "l_phi/F");
  t->Branch("l_m",              &summary.l_m,           "l_m/F");
  t->Branch("met_pt",              &summary.met_pt,           "met_pt/F");
  t->Branch("met_phi",              &summary.met_phi,           "met_phi/F");
  t->Branch("nj",              &summary.nj,           "nj/I");
  t->Branch("j_pt",          summary.j_pt,           "j_pt[nj]/F");
  t->Branch("j_eta",              summary.j_eta,           "j_eta[nj]/F");
  t->Branch("j_phi",              summary.j_phi,           "j_phi[nj]/F");
  t->Branch("j_m",              summary.j_m,           "j_m[nj]/F");
  t->Branch("j_csv",              summary.j_csv,           "j_csv[nj]/F");
  t->Branch("j_PtD",              summary.j_PtD,           "j_PtD[nj]/F");
  t->Branch("j_nch",              summary.j_nch,           "j_nch[nj]/F");
  t->Branch("j_width",              summary.j_width,           "j_width[nj]/F");
  t->Branch("j_tau21",              summary.j_tau21,           "j_tau21[nj]/F");
  t->Branch("j_tau32",              summary.j_tau32,           "j_tau32[nj]/F");
  t->Branch("j_tau43",              summary.j_tau43,           "j_tau43[nj]/F");
  t->Branch("j_zg",              summary.j_zg,           "j_zg[nj]/F");
  t->Branch("ue_nch",              &summary.ue_nch,           "ue_nch/F");
  t->Branch("ue_chsumpt",              &summary.ue_chsumpt,           "ue_chsumpt/F");
  t->Branch("ue_chsumpz",              &summary.ue_chsumpz,           "ue_chsumpz/F");
}
  
