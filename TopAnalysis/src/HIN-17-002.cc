#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "TSystem.h"

#include "TopLJets2015/TopAnalysis/interface/LeptonEfficiencyWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/BtagUncertaintyComputer.h"
#include "TopLJets2015/TopAnalysis/interface/TOP-HIForest.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include <string>
#include <vector>
#include <algorithm>
#include <limits>

using namespace std;

//
void RunHin17002(TString inFileName,
                 TString outFileName,
                 Int_t channelSelection,
                 Int_t chargeSelection,
                 TH1F *normH,
                 Bool_t runSysts,
                 TString era)
{		 
  
  Float_t JETPTTHRESHOLD=25;
  Float_t DRJJTHRESHOLD=2.0;

  if(inFileName=="") 
    {
      std::cout << "No inputs specified. return" << std::endl;
      return;
    }

  bool isMC(false);
  if(inFileName.Contains("/MC") || inFileName.Contains("PYQUEN") || inFileName.Contains("Pyquen")) isMC=true;
  bool isTTJets(false);
  if(inFileName.Contains("/MCTTNominal") || inFileName.Contains("TTBAR")) isTTJets=true;

  float totalEvtNorm(1.0);
  if(isMC && normH) totalEvtNorm=normH->GetBinContent(1);
  if(!isMC)
    {
      if( (channelSelection==11 || channelSelection==1100) && inFileName.Contains("FilteredSingleMuHighPt")) return;
      if( (channelSelection==13 || channelSelection==1300) && inFileName.Contains("HighPtLowerPhotons"))     return;
      runSysts=false;
    }
  std::cout << "Will process " << inFileName << " and save the results in " << outFileName << endl
	    << "Sample will be treated as MC=" << isMC <<  std::endl
	    << "Systematics will be run=" << runSysts << std::endl
	    << "Corrections to be retrieved from era=" << era << std::endl
	    << "Total normalization factor=" << totalEvtNorm << std::endl;

  std::map<TString, TGraphAsymmErrors *> expBtagEff;
  BTagSFUtil myBTagSFUtil;
  if(isMC)
    { 
      TString btagEffExpUrl(era+"/expTageff.root");
      gSystem->ExpandPathName(era+"/expTageff.root");   
      TFile *beffIn=TFile::Open(btagEffExpUrl);
      expBtagEff["b"]=(TGraphAsymmErrors *)beffIn->Get("b");
      expBtagEff["c"]=(TGraphAsymmErrors *)beffIn->Get("c");
      expBtagEff["udsg"]=(TGraphAsymmErrors *)beffIn->Get("udsg");
      beffIn->Close();
      cout << "Read " << expBtagEff.size() << " b-tag efficiency expectations from " << btagEffExpUrl << endl;
    }

  std::vector<TString>* inFileNames_p = new std::vector<TString>;
  inFileNames_p->push_back(inFileName);
  const Int_t nFiles = (Int_t)inFileNames_p->size();

  LeptonEfficiencyWrapper lepEffH(!isMC,era.ReplaceAll("era5TeV","era2015"));

  Int_t nSysts(0);
  TString lselTxt( channelSelection==13 ? "m" : (channelSelection==11 ? "e" : "") );
  TString expSysts[]={"btagup","btagdn","othertagup","othertagdn","jesup","jesdn","jerup","jerdn",lselTxt+"effup",lselTxt+"effdn"};
  
  //prepare output
  TFile* outFile_p = new TFile(outFileName, "RECREATE");

  //book tree
  outFile_p->cd();
  
  LJEvent_t ljev;
  TTree *outT=new TTree("data","data");
  defineTreeBranches(outT,ljev,isMC);
  outT->SetDirectory(outFile_p);

  //book histograms
  std::map<TString,TH1 *> histos;
  histos["wgtcounter"] = new TH1F("wgtcounter",";Weight;Events;",200,0,200);
  histos["fidcounter"] = new TH1F("fidcounter",";Weight;Events;",200,0,200);
  histos["gencounter"] = new TH1F("gencounter",";Step;Events;",4, 1,5);
  TString  genstep[4] = {"Initial", "#req 1Lepton", "#req 1Lepton fiducial", "#req 2jets"};
  for (int i=1; i < histos["gencounter"]->GetNbinsX()+1; i++)
    histos["gencounter"]->GetXaxis()->SetBinLabel(i,genstep[i-1]);
  histos["recocounter"] = new TH1F("recocounter",";Step;Events;",8, 1,9);
  TString  recostep[8] = {"Initial", "Trigger", "#req 1Lepton", "#equiv 1Lepton", "#req 4jets","#req 1 b-tags","#equiv 1 b-tag", "#equiv 2btag"};
  for (int i=1; i < histos["recocounter"]->GetNbinsX()+1; i++)
     histos["recocounter"]->GetXaxis()->SetBinLabel(i,recostep[i-1]);
  histos["trig"] = new TH1F("trig",";Trigger;Events",2,0,2);
  histos["lpt"]  = new TH1F("lpt",";Transverse momentum [GeV];Events",20.,0.,200.);
  histos["leta"] = new TH1F("leta",";Pseudo-rapidity;Events",20.,0.,2.1);
  histos["mt"]   = new TH1F("mt",";Transverse Mass [GeV];Events" ,20,0.,200.);
  histos["metpt"]= new TH1F("metpt",";Missing transverse energy [GeV];Events" ,20,0.,200.);
  histos["drlj"] = new TH1F("drlj",";min#DeltaR(l,j) [GeV];Events" ,63,0.,6.3);
  histos["drlb"] = new TH1F("drlb",";min#DeltaR(l,j) [GeV];Events" ,63,0.,6.3);
  histos["ttevttype"] = new TH1F("ttevttype",";t#bar{t} event type;Events" ,9,0.,9.);

   //electron selection control plots
  for(int ireg=0; ireg<2; ireg++)
    {
      TString pf(ireg==0 ? "_ee" : "_eb");
      histos["sigmaietaieta"+pf]=new TH1F("sigmaietaieta"+pf,";#sigma(i#eta,i#eta);Electrons",25,0,0.04);
      histos["detain"+pf]=new TH1F("detain"+pf,";#Delta#eta(in);Electrons",25,-0.015,0.15);
      histos["dphiin"+pf]=new TH1F("dephiin"+pf,";#Delta#phi(in);Electrons",25,-0.05,0.05);
      histos["hovere"+pf]=new TH1F("hovere"+pf,";H/E;Electrons",25,0,0.05);
      histos["eoverpinv"+pf]=new TH1F("eoverpinv"+pf,";1/E-1/p [GeV^{-1}];Electrons",25,0,0.1);
      histos["d0"+pf]=new TH1F("d0"+pf,";d_{0} [cm];Electrons",25,-0.05,0.05);
      histos["dz"+pf]=new TH1F("dz"+pf,";d_{z} [cm];Electrons",25,-0.5,0.5);
      histos["misshits"+pf]=new TH1F("misshits"+pf,";Missing hits;Electrons",4,0,4);
      histos["convveto"+pf]=new TH1F("convveto"+pf,";Conversion veto flag;Electrons",2,0,2);
      histos["reliso"+pf]=new TH1F("reliso"+pf,";Relative isolation;Electrons",25,0,0.1);
      histos["lpt"+pf]  = new TH1F("lpt"+pf,";Transverse momentum [GeV];Events",20.,0.,200.);
      histos["leta"+pf] = new TH1F("leta"+pf,";Pseudo-rapidity;Events",20.,0.,2.1);
    }

  //per b-tag multiplicity control plots
  for(int ij=0; ij<=2; ij++)
    {
      TString pf(Form("%db",ij));
      histos["lpt_"+pf]    = new TH1F("lpt_"+pf,";Transverse momentum [GeV];Events",5.,20.,120.);
      histos["leta_"+pf]   = new TH1F("leta_"+pf,";Pseudo-rapidity;Events",5.,0.,2.2);
      TString region1("ee_"); TString region2("eb_");
      if( (channelSelection==11 || channelSelection==1100) )
	{
	  histos["lpt_"+region1+pf]    = new TH1F("lpt_"+region1+pf,";Transverse momentum [GeV];Events",5.,20.,120.); 
	  histos["lpt_"+region2+pf]    = new TH1F("lpt_"+region2+pf,";Transverse momentum [GeV];Events",5.,20.,120.);
	  histos["leta_"+region1+pf]   = new TH1F("leta_"+region1+pf,";Pseudo-rapidity;Events",5.,0.,2.2);
	  histos["leta_"+region2+pf]   = new TH1F("leta_"+region2+pf,";Pseudo-rapidity;Events",5.,0.,2.2); 
	}
      histos["jpt_"+pf]    = new TH1F("jpt_"+pf,";Transverse momentum [GeV];Events",5.,0.,250.);
      histos["jeta_"+pf]   = new TH1F("jeta_"+pf,";Pseudo-rapidity;Events",5.,0.,2.5);
      histos["ht_"+pf]     = new TH1F("ht_"+pf,";H_{T} [GeV];Events",10.,0.,800.);
      histos["metpt_"+pf]  = new TH1F("metpt_"+pf,";Missing transverse energy [GeV];Events" ,5,0.,200.);
      histos["metphi_"+pf] = new TH1F("metphi_" + pf,";MET #phi [rad];Events" ,5,-3.2,3.2);
      histos["mt_"+pf]     = new TH1F("mt_"+pf,";Transverse Mass [GeV];Events" ,10,0.,300.);
      histos["mjj_"+pf]    = new TH1F("mjj_"+pf,";Mass(j,j') [GeV];Events" ,20,0.,400.);
      histos["rankedmjj_"+pf]    = new TH1F("rankedmjj_"+pf,";Mass(j,j') [GeV];Events" ,20,0.,400.);
      histos["rankedq70mjj_"+pf]    = new TH1F("rankedq70mjj_"+pf,";Mass(j,j') [GeV];Events" ,20,0.,400.);
      histos["drjj_"+pf]    = new TH1F("drjj_"+pf,";min#DeltaR(j,j') [GeV];Events" ,12,0.,6.3);      
      histos["ptjj_"+pf]    = new TH1F("ptjj_"+pf,";p_{T}(j,j') [GeV];Events" ,20,0.,300.);
      histos["etajj_"+pf]    = new TH1F("etajj_"+pf,";#eta(j,j');Events" ,10,-3.,3.);
      histos["mlb_"+pf]    = new TH1F("mlb_"+pf,";Mass(l,b) [GeV];Events" ,20,0.,300.);
      histos["njets_"+pf]  = new TH1F("njets_"+pf,";Jet multiplicity;Events" ,6,2.,8.);
      histos["mbjj_"+pf]    = new TH1F("mbjj_"+pf,";Mass(bjj') [GeV];Events" ,20,0.,400.);
      histos["rankedmbjj_"+pf]    = new TH1F("rankedmbjj_"+pf,";Mass(bjj') [GeV];Events" ,20,0.,400.);

      if(isMC && runSysts)
	{
	  nSysts=sizeof(expSysts)/sizeof(TString);
	  histos["mjjshapes_"+pf+"_exp"]=new TH2F("mjjshapes_"+pf+"_exp",";Mass(j,j');Systematic uncertainty;Events",20,0,400,nSysts,0,nSysts);
	  histos["drjjshapes_"+pf+"_exp"]=new TH2F("drjjshapes_"+pf+"_exp",";min#DeltaR(j,j');Systematic uncertainty;Events",12,0,6.3,nSysts,0,nSysts);
	  histos["rankedmjjshapes_"+pf+"_exp"]=new TH2F("rankedmjjshapes_"+pf+"_exp",";Mass(j,j');Systematic uncertainty;Events",20,0,400,nSysts,0,nSysts);
	  histos["rankedq70mjjshapes_"+pf+"_exp"]=new TH2F("rankedq70mjjshapes_"+pf+"_exp",";Mass(j,j');Systematic uncertainty;Events",20,0,400,nSysts,0,nSysts);
	  for(int i=0; i<nSysts; i++)
	    {
	      histos["drjjshapes_"+pf+"_exp"]->GetYaxis()->SetBinLabel(i+1,expSysts[i]);
	      histos["mjjshapes_"+pf+"_exp"]->GetYaxis()->SetBinLabel(i+1,expSysts[i]);
	      histos["rankedmjjshapes_"+pf+"_exp"]->GetYaxis()->SetBinLabel(i+1,expSysts[i]);
	      histos["rankedq70mjjshapes_"+pf+"_exp"]->GetYaxis()->SetBinLabel(i+1,expSysts[i]);
	    }
	  
	  histos["mjjshapes_"+pf+"_gen"]=new TH2F("mjjshapes_"+pf+"_gen",";Mass(j,j') [GeV];Systematic uncertainty;Events",20,0,400,200,0,200);
	  histos["drjjshapes_"+pf+"_gen"]=new TH2F("drjjshapes_"+pf+"_gen",";min#DeltaR(j,j');Systematic uncertainty;Events",12,0,6.3,200,0,200);
	  histos["rankedmjjshapes_"+pf+"_gen"]=new TH2F("rankedmjjshapes_"+pf+"_gen",";Mass(j,j') [GeV];Systematic uncertainty;Events",20,0,400,200,0,200);
	  histos["rankedq70mjjshapes_"+pf+"_gen"]=new TH2F("rankedq70mjjshapes_"+pf+"_gen",";Mass(j,j') [GeV];Systematic uncertainty;Events",20,0,400,200,0,200);
	  for(int i=0; i<200;i++)
	    {
	      histos["drjjshapes_"+pf+"_gen"]->GetYaxis()->SetBinLabel(i+1,Form("genUnc%d",i));
	      histos["mjjshapes_"+pf+"_gen"]->GetYaxis()->SetBinLabel(i+1,Form("genUnc%d",i));
	      histos["rankedmjjshapes_"+pf+"_gen"]->GetYaxis()->SetBinLabel(i+1,Form("genUnc%d",i));
	      histos["rankedq70mjjshapes_"+pf+"_gen"]->GetYaxis()->SetBinLabel(i+1,Form("genUnc%d",i));
	    }
	}
    }

  //prepare histograms
  for(std::map<TString,TH1 *>::iterator it=histos.begin();
      it!=histos.end();
      it++)
    {
      it->second->Sumw2();
      it->second->SetDirectory(0);
    }

  for(Int_t fileIter = 0; fileIter < nFiles; fileIter++){

    TString inF(inFileNames_p->at(fileIter));
    if(inF.Contains("/store") && !inF.Contains("root:")) inF="root://eoscms//eos/cms/"+inF;
    TFile* inFile_p = TFile::Open(inF, "READ");
    
    TTree* lepTree_p = (TTree*)inFile_p->Get("ggHiNtuplizer/EventTree");
    TTree* jetTree_p = (TTree*)inFile_p->Get("ak4PFJetAnalyzer/t");
    TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");
    TTree* hltTree_p = (TTree*)inFile_p->Get("hltanalysis/HltTree");
    TTree *pfCand_p  = (TTree *)inFile_p->Get("pfcandAnalyzer/pfTree");
    TTree *hiFJRho_p  = (TTree *)inFile_p->Get("hiFJRhoAnalyzer/t");
    TTree *pseudotop_p = (TTree *)inFile_p->Get("topGenAnalyzer/t");
    TTree *mc_p = (TTree *)inFile_p->Get("HiGenParticleAna/hi");

    //PF candidates
    std::vector<int> *pfId_p=0;
    std::vector<float> *pfPt_p=0,*pfEta_p=0,*pfPhi_p=0,*pfEnergy_p=0;
    pfCand_p->SetBranchStatus("pfId",     1);
    pfCand_p->SetBranchStatus("pfPt",     1);
    pfCand_p->SetBranchStatus("pfEta",    1);
    pfCand_p->SetBranchStatus("pfPhi",    1);
    pfCand_p->SetBranchStatus("pfEnergy", 1);
    pfCand_p->SetBranchAddress("pfId",     &pfId_p);
    pfCand_p->SetBranchAddress("pfPt",     &pfPt_p);
    pfCand_p->SetBranchAddress("pfEta",    &pfEta_p);
    pfCand_p->SetBranchAddress("pfPhi",    &pfPhi_p);
    pfCand_p->SetBranchAddress("pfEnergy", &pfEnergy_p);


    //muon variables
    std::vector<float>* muPt_p = 0;
    std::vector<float>* muPhi_p = 0;
    std::vector<float>* muEta_p = 0;
    std::vector<int>* muChg_p = 0;
    std::vector<float>* muChi2NDF_p = 0;
    std::vector<float>* muInnerD0_p = 0;
    std::vector<float>* muInnerDz_p = 0;
    std::vector<int>* muMuonHits_p = 0;
    std::vector<int>* muStations_p = 0;
    std::vector<int>* muTrkLayers_p = 0;
    std::vector<int>* muPixelHits_p = 0;    
    std::vector<float> *muPFChIso_p=0,*muPFPhoIso_p=0,*muPFNeuIso_p=0,*muPFPUIso_p=0;
    lepTree_p->SetBranchStatus("*", 0);
    lepTree_p->SetBranchStatus("muPt", 1);
    lepTree_p->SetBranchStatus("muPhi", 1);
    lepTree_p->SetBranchStatus("muEta", 1);
    lepTree_p->SetBranchStatus("muCharge", 1);
    lepTree_p->SetBranchStatus("muChi2NDF", 1);
    lepTree_p->SetBranchStatus("muInnerD0", 1);
    lepTree_p->SetBranchStatus("muInnerDz", 1);
    lepTree_p->SetBranchStatus("muMuonHits", 1);
    lepTree_p->SetBranchStatus("muStations", 1);
    lepTree_p->SetBranchStatus("muTrkLayers", 1);
    lepTree_p->SetBranchStatus("muPixelHits", 1);    
    lepTree_p->SetBranchStatus("muPFChIso", 1);
    lepTree_p->SetBranchStatus("muPFPhoIso", 1);
    lepTree_p->SetBranchStatus("muPFNeuIso", 1);
    lepTree_p->SetBranchStatus("muPFPUIso", 1);
    lepTree_p->SetBranchAddress("muPt", &muPt_p);
    lepTree_p->SetBranchAddress("muPhi", &muPhi_p);
    lepTree_p->SetBranchAddress("muEta", &muEta_p);
    lepTree_p->SetBranchAddress("muCharge", &muChg_p);
    lepTree_p->SetBranchAddress("muChi2NDF", &muChi2NDF_p);
    lepTree_p->SetBranchAddress("muInnerD0", &muInnerD0_p);
    lepTree_p->SetBranchAddress("muInnerDz", &muInnerDz_p);
    lepTree_p->SetBranchAddress("muMuonHits", &muMuonHits_p);
    lepTree_p->SetBranchAddress("muStations", &muStations_p);
    lepTree_p->SetBranchAddress("muTrkLayers", &muTrkLayers_p);
    lepTree_p->SetBranchAddress("muPixelHits", &muPixelHits_p);    
    lepTree_p->SetBranchAddress("muPFChIso", &muPFChIso_p);
    lepTree_p->SetBranchAddress("muPFPhoIso", &muPFPhoIso_p);
    lepTree_p->SetBranchAddress("muPFNeuIso", &muPFNeuIso_p);
    lepTree_p->SetBranchAddress("muPFPUIso", &muPFPUIso_p);

    //electron variables
    std::vector<float>* elePt_p = 0;
    std::vector<float>* elePhi_p = 0;
    std::vector<float>* eleEta_p = 0;
    std::vector<float>* eleSigmaIEtaIEta_p = 0;
    std::vector<float>* eledEtaAtVtx_p = 0;
    std::vector<float>* eledPhiAtVtx_p = 0;
    std::vector<float>* eleHoverE_p = 0;
    std::vector<float>* eleEoverP_p = 0;
    std::vector<float>* eleD0_p = 0;
    std::vector<float>* eleDz_p = 0;
    std::vector<float>* eleMissHits_p = 0;
    std::vector<float>* elepassConversionVeto_p = 0;
    std::vector<float>* elePFChIso_p=0, *elePFPhoIso_p=0, *elePFNeuIso_p=0, *elePFPUIso_p=0, *eleEffAreaTimesRho_p=0;
    std::vector<int>*   eleIDVeto_p=0;
    std::vector<int>*   eleIDLoose_p=0;
    std::vector<int>*   eleIDMedium_p=0;
    std::vector<int>*   eleIDTight_p=0;
    std::vector<int>*   eleCharge_p=0;
    lepTree_p->SetBranchStatus("elePt", 1);
    lepTree_p->SetBranchStatus("elePhi", 1);
    lepTree_p->SetBranchStatus("eleEta", 1);
    lepTree_p->SetBranchStatus("eleSigmaIEtaIEta", 1);
    lepTree_p->SetBranchStatus("eledEtaAtVtx", 1);
    lepTree_p->SetBranchStatus("eledPhiAtVtx", 1);
    lepTree_p->SetBranchStatus("eleHoverE", 1);
    lepTree_p->SetBranchStatus("eleEoverPInv", 1);
    lepTree_p->SetBranchStatus("eleD0", 1);
    lepTree_p->SetBranchStatus("eleDz", 1);
    lepTree_p->SetBranchStatus("eleMissHits", 1);
    lepTree_p->SetBranchStatus("elepassConversionVeto", 1);
    lepTree_p->SetBranchStatus("elePFChIso", 1);
    lepTree_p->SetBranchStatus("elePFPhoIso", 1);
    lepTree_p->SetBranchStatus("elePFNeuIso", 1);
    lepTree_p->SetBranchStatus("elePFPUIso", 1);
    lepTree_p->SetBranchStatus("eleEffAreaTimesRho", 1);
    lepTree_p->SetBranchStatus("eleID*", 1);
    lepTree_p->SetBranchStatus("eleCharge", 1);
    lepTree_p->SetBranchAddress("elePt", &elePt_p);
    lepTree_p->SetBranchAddress("elePhi", &elePhi_p);
    lepTree_p->SetBranchAddress("eleEta", &eleEta_p);
    lepTree_p->SetBranchAddress("eleSigmaIEtaIEta", &eleSigmaIEtaIEta_p);
    lepTree_p->SetBranchAddress("eledEtaAtVtx", &eledEtaAtVtx_p);
    lepTree_p->SetBranchAddress("eledPhiAtVtx", &eledPhiAtVtx_p);
    lepTree_p->SetBranchAddress("eleHoverE", &eleHoverE_p);
    lepTree_p->SetBranchAddress("eleEoverPInv", &eleEoverP_p);
    lepTree_p->SetBranchAddress("eleD0", &eleD0_p);
    lepTree_p->SetBranchAddress("eleDz", &eleDz_p);
    lepTree_p->SetBranchAddress("eleMissHits", &eleMissHits_p);
    lepTree_p->SetBranchAddress("elepassConversionVeto", &elepassConversionVeto_p);
    lepTree_p->SetBranchAddress("elePFChIso", &elePFChIso_p);
    lepTree_p->SetBranchAddress("elePFPhoIso", &elePFPhoIso_p);
    lepTree_p->SetBranchAddress("elePFNeuIso", &elePFNeuIso_p);
    lepTree_p->SetBranchAddress("elePFPUIso", &elePFPUIso_p);
    lepTree_p->SetBranchAddress("eleEffAreaTimesRho", &eleEffAreaTimesRho_p);
    lepTree_p->SetBranchAddress("eleIDVeto", &eleIDVeto_p);
    lepTree_p->SetBranchAddress("eleIDLoose", &eleIDLoose_p);
    lepTree_p->SetBranchAddress("eleIDMedium", &eleIDMedium_p);
    lepTree_p->SetBranchAddress("eleIDTight", &eleIDTight_p);
    lepTree_p->SetBranchAddress("eleCharge", &eleCharge_p);

    //gen-level variables
    std::vector<int> *mcPID=0,*mcStatus=0;
    std::vector<float> *mcPt=0,*mcEta=0,*mcPhi=0,*mcMass=0;
    if(mc_p)
      {
        mc_p->SetBranchStatus("pdg", 1);
        mc_p->SetBranchStatus("sta", 1);
        mc_p->SetBranchStatus("pt", 1);
        mc_p->SetBranchStatus("eta", 1);
        mc_p->SetBranchStatus("phi", 1);
        mc_p->SetBranchStatus("mass", 1);
        mc_p->SetBranchAddress("sta", &mcStatus);
        mc_p->SetBranchAddress("pdg", &mcPID);
        mc_p->SetBranchAddress("pt", &mcPt);
        mc_p->SetBranchAddress("eta", &mcEta);
        mc_p->SetBranchAddress("phi", &mcPhi);
        mc_p->SetBranchAddress("mass", &mcMass);
      }

    //jet variables
    const int maxJets = 5000;
    Int_t   nref,ngen;
    Float_t jtpt[maxJets],genpt[maxJets];
    Float_t jteta[maxJets],geneta[maxJets];
    Float_t jtphi[maxJets],genphi[maxJets];
    Float_t jtarea[maxJets]; 
    Float_t jtm[maxJets]; 
    Float_t discr_csvV2[maxJets];
    Float_t refpt[maxJets], refarea[maxJets], refdrjt[maxJets], gendrjt[maxJets];
    Float_t jPfCHF[maxJets],jPfNHF[maxJets],jPfCEF[maxJets],jPfNEF[maxJets], jPfMUF[maxJets];//jet id
    Int_t jPfCHM[maxJets],jPfNHM[maxJets], jPfCEM[maxJets],jPfNEM[maxJets], jPfMUM[maxJets]; //jet id
    Int_t refparton_flavorForB[maxJets], refparton_flavor[maxJets], genmatchindex[maxJets];
    jetTree_p->SetBranchStatus("*", 0);
    jetTree_p->SetBranchStatus("ngen", 1);
    jetTree_p->SetBranchStatus("genpt", 1);
    jetTree_p->SetBranchStatus("genphi", 1);
    jetTree_p->SetBranchStatus("geneta", 1);
    jetTree_p->SetBranchStatus("gendrjt", 1);
    jetTree_p->SetBranchStatus("genmatchindex", 1);
    jetTree_p->SetBranchStatus("nref", 1);
    jetTree_p->SetBranchStatus("jtpt", 1);
    jetTree_p->SetBranchStatus("jtphi", 1);
    jetTree_p->SetBranchStatus("jteta", 1);
    jetTree_p->SetBranchStatus("jtarea", 1);
    jetTree_p->SetBranchStatus("jtm", 1);
    jetTree_p->SetBranchStatus("jtPfCHF", 1);
    jetTree_p->SetBranchStatus("jtPfNHF", 1);
    jetTree_p->SetBranchStatus("jtPfCEF", 1);
    jetTree_p->SetBranchStatus("jtPfNEF", 1);
    jetTree_p->SetBranchStatus("jtPfMUF", 1);
    jetTree_p->SetBranchStatus("jtPfCHM", 1);
    jetTree_p->SetBranchStatus("jtPfNHM", 1);
    jetTree_p->SetBranchStatus("jtPfCEM", 1);
    jetTree_p->SetBranchStatus("jtPfNEM", 1);
    jetTree_p->SetBranchStatus("jtPfMUM", 1);
    jetTree_p->SetBranchStatus("discr_csvV2", 1);
    jetTree_p->SetBranchStatus("refpt", 1);
    jetTree_p->SetBranchStatus("refarea", 1);
    jetTree_p->SetBranchStatus("refdrjt", 1);
    jetTree_p->SetBranchStatus("refparton_flavorForB", 1);
    jetTree_p->SetBranchStatus("refparton_flavor", 1);
    jetTree_p->SetBranchAddress("nref", &nref);
    jetTree_p->SetBranchAddress("jtpt", jtpt);
    jetTree_p->SetBranchAddress("jtphi", jtphi);
    jetTree_p->SetBranchAddress("jteta", jteta);
    jetTree_p->SetBranchAddress("jtarea", jtarea);
    jetTree_p->SetBranchAddress("jtm", jtm);
    jetTree_p->SetBranchAddress("jtPfCHF", jPfCHF);
    jetTree_p->SetBranchAddress("jtPfNHF", jPfNHF);
    jetTree_p->SetBranchAddress("jtPfCEF", jPfCEF);
    jetTree_p->SetBranchAddress("jtPfNEF", jPfNEF);
    jetTree_p->SetBranchAddress("jtPfMUF", jPfMUF);
    jetTree_p->SetBranchAddress("jtPfCHM", jPfCHM);
    jetTree_p->SetBranchAddress("jtPfNHM", jPfNHM);
    jetTree_p->SetBranchAddress("jtPfCEM", jPfCEM);
    jetTree_p->SetBranchAddress("jtPfNEM", jPfNEM);
    jetTree_p->SetBranchAddress("jtPfMUM", jPfMUM);
    jetTree_p->SetBranchAddress("discr_csvV2", discr_csvV2);
    jetTree_p->SetBranchAddress("refpt", refpt);
    jetTree_p->SetBranchAddress("refarea", refarea);
    jetTree_p->SetBranchAddress("refdrjt", refdrjt);
    jetTree_p->SetBranchAddress("refparton_flavorForB", refparton_flavorForB);
    jetTree_p->SetBranchAddress("refparton_flavor", refparton_flavor);
    jetTree_p->SetBranchAddress("ngen", &ngen);
    jetTree_p->SetBranchAddress("genpt", genpt);
    jetTree_p->SetBranchAddress("genphi", genphi);
    jetTree_p->SetBranchAddress("geneta", geneta);
    jetTree_p->SetBranchAddress("gendrjt", gendrjt);
    jetTree_p->SetBranchAddress("genmatchindex", genmatchindex);

    //jet related variables from the hiFJRho tree
    std::vector<float>* etaMin_hiFJRho_p = 0;
    std::vector<float>* etaMax_hiFJRho_p = 0;
    std::vector<float>* rho_hiFJRho_p = 0;
    std::vector<float>* rhom_hiFJRho_p = 0;
    std::vector<float>* rhoCorr_hiFJRho_p = 0;
    std::vector<float>* rhomCorr_hiFJRho_p = 0;
    std::vector<float>* rhoCorr1Bin_hiFJRho_p = 0;
    std::vector<float>* rhomCorr1Bin_hiFJRho_p = 0;
    
    hiFJRho_p->SetBranchStatus("*", 0);
    hiFJRho_p->SetBranchStatus("etaMin", 1);
    hiFJRho_p->SetBranchStatus("etaMax", 1);
    hiFJRho_p->SetBranchStatus("rho", 1);
    hiFJRho_p->SetBranchStatus("rhom", 1);
    hiFJRho_p->SetBranchStatus("rhoCorr", 1);
    hiFJRho_p->SetBranchStatus("rhomCorr", 1);
    hiFJRho_p->SetBranchStatus("rhoCorr1Bin", 1);
    hiFJRho_p->SetBranchStatus("rhomCorr1Bin", 1);
    
    hiFJRho_p->SetBranchAddress("etaMin", &etaMin_hiFJRho_p);
    hiFJRho_p->SetBranchAddress("etaMax", &etaMax_hiFJRho_p);
    hiFJRho_p->SetBranchAddress("rho", &rho_hiFJRho_p);
    hiFJRho_p->SetBranchAddress("rhom", &rhom_hiFJRho_p);
    hiFJRho_p->SetBranchAddress("rhoCorr", &rhoCorr_hiFJRho_p);
    hiFJRho_p->SetBranchAddress("rhomCorr", &rhomCorr_hiFJRho_p);
    hiFJRho_p->SetBranchAddress("rhoCorr1Bin", &rhoCorr1Bin_hiFJRho_p);
    hiFJRho_p->SetBranchAddress("rhomCorr1Bin", &rhomCorr1Bin_hiFJRho_p);



    //event variables
    UInt_t run_, lumi_;
    ULong64_t evt_;
    Int_t hiBin_;
    Float_t hiHFplus_, hiHFminus_, hiHFplusEta4_, hiHFminusEta4_;
    Float_t vz_;
    Float_t weight;
    std::vector<float> *ttbar_w_p=0;
    hiTree_p->SetBranchStatus("*", 0);
    hiTree_p->SetBranchStatus("run", 1);
    hiTree_p->SetBranchStatus("evt", 1);
    hiTree_p->SetBranchStatus("lumi", 1);
    hiTree_p->SetBranchStatus("hiBin", 1);
    hiTree_p->SetBranchAddress("hiHFplus", &hiHFplus_);
    hiTree_p->SetBranchAddress("hiHFminus", &hiHFminus_);
    hiTree_p->SetBranchAddress("hiHFplusEta4", &hiHFplusEta4_);
    hiTree_p->SetBranchAddress("hiHFminusEta4", &hiHFminusEta4_);
    hiTree_p->SetBranchStatus("vz", 1);
    hiTree_p->SetBranchStatus("weight", 1);
    hiTree_p->SetBranchStatus("ttbar_w",1);
    hiTree_p->SetBranchAddress("run", &run_);
    hiTree_p->SetBranchAddress("evt", &evt_);
    hiTree_p->SetBranchAddress("lumi", &lumi_);
    hiTree_p->SetBranchAddress("hiBin", &hiBin_);
    hiTree_p->SetBranchAddress("vz", &vz_);
    hiTree_p->SetBranchAddress("weight", &weight);
    hiTree_p->SetBranchAddress("ttbar_w",&ttbar_w_p);

    //trigger
    int trig = 0;
    std::string triggerName;
    if(channelSelection==13 || channelSelection==1300)
      {
	triggerName = isMC ? "HLT_PAL3Mu15_v1" : "HLT_PAL3Mu15_v1";
      }
    
    if(channelSelection==11 || channelSelection==1100) 
      {
	triggerName = isMC ? "HLT_PASinglePhoton30_Eta3p1_v1" : "HLT_PASinglePhoton30_Eta3p1_v1"; 
      }
    
    hltTree_p->SetBranchStatus(triggerName.data(),1);
    hltTree_p->SetBranchAddress(triggerName.data(),&trig);

    //pseudo-top (if available)
    std::vector<int> *pt_pdgid=0;
    std::vector<double> *pt_pt=0,*pt_eta=0,*pt_phi=0, *pt_m=0;
    std::vector<std::vector<int> > *pt_daughterArr=0;

    if(pseudotop_p)
      {
        pseudotop_p->SetBranchStatus("*", 1);
        pseudotop_p->SetBranchAddress("gtop_id", &pt_pdgid);
        pseudotop_p->SetBranchAddress("gtop_pt", &pt_pt);
        pseudotop_p->SetBranchAddress("gtop_eta", &pt_eta);
        pseudotop_p->SetBranchAddress("gtop_phi", &pt_phi);
        pseudotop_p->SetBranchAddress("gtop_mass", &pt_m);
        pseudotop_p->SetBranchAddress("gtop_daughterArr", &pt_daughterArr);
      }
    
    Int_t nEntries = (Int_t)lepTree_p->GetEntries();
    
    std::cout << "Analysing " << nEntries << " events isMC=" << isMC
	      << " trigger=" << triggerName << endl;

    for(Int_t entry = 0; entry < nEntries; entry++)
      {
	if(entry%1000==0)
	  {
	    printf("\r [%d/%d] done",entry,nEntries);
	    cout << flush;
	  }
	histos["recocounter"]->Fill(1);
	//reset summary tree
	ljev.nj=0; ljev.ngj=0; ljev.ngp=0; ljev.nb=0; ljev.l_id=0; ljev.w=0;	ljev.npt=0;

	//readout this event
	lepTree_p->GetEntry(entry);
	jetTree_p->GetEntry(entry);
	hiTree_p->GetEntry(entry);
	hltTree_p->GetEntry(entry);
	pfCand_p->GetEntry(entry);
	hiFJRho_p->GetEntry(entry);
        if(mc_p) mc_p->GetEntry(entry);

	//pseudo-top
	if(pseudotop_p)
	  {
	    pseudotop_p->GetEntry(entry);
	    ljev.npt=pt_pdgid->size();
	    for(size_t ipt=0; ipt<pt_pdgid->size(); ipt++)
	      {
		ljev.pt_id[ipt]= ipt;
		ljev.pt_pdgid[ipt]= (*pt_pdgid)[ipt];
		ljev.pt_pt[ipt]= (*pt_pt)[ipt];
		ljev.pt_eta[ipt]= (*pt_eta)[ipt];
		ljev.pt_phi[ipt]= (*pt_phi)[ipt];
		ljev.pt_m[ipt]= (*pt_m)[ipt];
	      }
	  }


	//assign an event weight
	float evWeight(1.0);
        float ttevtype(0);
	if(isMC)
	  {
	    histos["gencounter"]->Fill(1);
	    if(ttbar_w_p->size()) evWeight = ttbar_w_p->at(0);
	    evWeight *= totalEvtNorm;

	    //fiducial region analysis
	    std::vector<TLorentzVector> selGenLeptons,otherLeptons;
	    if(channelSelection==13 || channelSelection==11)
	      {
		for(size_t imc=0; imc<mcPID->size(); imc++)
		  {
		    int status=abs(mcStatus->at(imc));
		    int abspid=abs(mcPID->at(imc));		
		    TLorentzVector p4;
		    p4.SetPtEtaPhiM(mcPt->at(imc),mcEta->at(imc),mcPhi->at(imc),mcMass->at(imc));
                    //hardprocess
		    if(status==3 && (abspid<=6 || abspid==24 || abspid==11 || abspid==13))
                      {
                        {
                          ljev.gp_pdgid[ljev.ngp]=mcPID->at(imc);
                          ljev.gp_pt[ljev.ngp]=p4.Pt();
                          ljev.gp_eta[ljev.ngp]=p4.Eta();
                          ljev.gp_phi[ljev.ngp]=p4.Phi();
                          ljev.gp_m[ljev.ngp]=p4.M();
                          ljev.ngp++;
                          if(abspid==11 || abspid==13) ttevtype+=1;
                        }
                      }
                    
		    //final state leptons
		    if(status!=1) continue;
		    if(channelSelection==13 && abspid!=13) continue;
		    if(channelSelection==11 && abspid!=11) continue;
		    if(fabs(p4.Eta())<2.1)
		      {
			if(abspid==13 && p4.Pt()>25) selGenLeptons.push_back(p4);
			if(abspid==11 && p4.Pt()>35) selGenLeptons.push_back(p4);
		      }
		    else if( (abspid==13 || abspid==11) && p4.Pt()>15 && fabs(p4.Eta())<2.5)
		      {
			otherLeptons.push_back(p4);
		      }
		  }
	      }

	    if(selGenLeptons.size()>0 || otherLeptons.size()>0)histos["gencounter"]->Fill(2);
	    if(selGenLeptons.size()>0)
	      {
		histos["gencounter"]->Fill(3);
		ljev.gl_pt=selGenLeptons[0].Pt();
		ljev.gl_eta=selGenLeptons[0].Eta();
		ljev.gl_phi=selGenLeptons[0].Phi();
		ljev.gl_m=selGenLeptons[0].M();
	      }
	    
	    //select gen jets cross-cleaning with leading muon
	    int nGenJets(0);
	    for(int imcj=0; imcj<ngen; imcj++)
	      {
		TLorentzVector p4;
		p4.SetPtEtaPhiM(genpt[imcj],geneta[imcj],genphi[imcj],0.);
		if(selGenLeptons.size() && p4.DeltaR(selGenLeptons[0])<0.4) continue;
		if(p4.Pt()<25 || fabs(p4.Eta())>2.4) continue;
		nGenJets++;
		
		if(ljev.ngj<20)
		  {
		    ljev.gj_pt[ljev.ngj]=p4.Pt();
		    ljev.gj_eta[ljev.ngj]=p4.Eta();
		    ljev.gj_phi[ljev.ngj]=p4.Phi();
		    ljev.gj_m[ljev.ngj]=p4.M();
		    ljev.gj_dr[ljev.ngj]=gendrjt[imcj];
		    ljev.gj_index[ljev.ngj]=genmatchindex[imcj];
		    ljev.ngj++;
		  }

	      }

	    //check if it passes the gen level acceptance
	    bool passFid(nGenJets>=2 && selGenLeptons.size()==1);
	    if(passFid) histos["gencounter"]->Fill(4);
	    for(size_t iw=0; iw<ttbar_w_p->size(); iw++)
	      {
		histos["wgtcounter"]->Fill(iw,ttbar_w_p->at(iw));
		if(passFid) histos["fidcounter"]->Fill(iw,ttbar_w_p->at(iw));
	      }
	  }

	//require trigger for the event
	histos["trig"]->Fill(trig,evWeight);
	if(!isMC && trig==0) continue;

	histos["recocounter"]->Fill(2);

	//select good muons
	//cf. details in https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
	std::vector<TLorentzVector> tightMuons,looseMuons,tightMuonsNonIso;
	std::vector<int> muonCharge;
	const Int_t nMu = (Int_t)muPt_p->size();
	for(Int_t muIter = 0; muIter < nMu; muIter++)
	  {
	    bool passLooseKin( muPt_p->at(muIter) > 15. && TMath::Abs(muEta_p->at(muIter))<2.4);
	    bool passTightKin( muPt_p->at(muIter) > 20. && TMath::Abs(muEta_p->at(muIter))<2.1);
	    bool passLooseId(true);
	    bool passTightId( passLooseId
			      && muChi2NDF_p->at(muIter) < 10
			      && muMuonHits_p->at(muIter) >0 		
			      && muStations_p->at(muIter) >1
			      && TMath::Abs(muInnerD0_p->at(muIter))<0.2
			      && TMath::Abs(muInnerDz_p->at(muIter))<0.5
			      && muPixelHits_p->at(muIter)>0		
			    && muTrkLayers_p->at(muIter)>5);
	    float relIso=(muPFChIso_p->at(muIter)+TMath::Max(muPFPhoIso_p->at(muIter)+muPFNeuIso_p->at(muIter)-0.5*muPFPUIso_p->at(muIter),0.))/muPt_p->at(muIter);
	    bool passTightIso( relIso<0.15);
	    bool passLooseIso( relIso<0.25);
	    
	    //save muon if good
	    TLorentzVector p4(0,0,0,0);
	    p4.SetPtEtaPhiM(muPt_p->at(muIter),muEta_p->at(muIter),muPhi_p->at(muIter), 0.1056583715);
	    if(passTightKin && passTightId && !passTightIso && relIso>0.2) 
	      {
		tightMuonsNonIso.push_back(p4);
		muonCharge.push_back(muChg_p->at(muIter));
	      }
	    if(passTightKin && passTightId && passTightIso)     
	      {
		tightMuons.push_back( p4 );
		muonCharge.push_back(muChg_p->at(muIter));
	      }
	    else if(passLooseKin && passLooseId && passLooseIso) looseMuons.push_back( p4 );
	  }
  

	//select good electronss
	//cf. details in https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
	std::vector<TLorentzVector> mediumElectrons,vetoElectrons,mediumElectronsNonIso,mediumElectronsFailId;
	std::vector<int> elCharge;
	const Int_t nEl = (Int_t)elePt_p->size();
	for(Int_t elIter = 0; elIter < nEl; elIter++)
	  {
	    bool passMediumPt(elePt_p->at(elIter) > 40.0);  
	    bool passMediumEta(fabs(eleEta_p->at(elIter)) < 2.1 && (fabs(eleEta_p->at(elIter)) < 1.4442 || fabs(eleEta_p->at(elIter)) > 1.5660));
	    bool passPt(elePt_p->at(elIter) > 15.0);  
	    bool passEta(fabs(eleEta_p->at(elIter)) < 2.5);
	    if(!passPt || !passEta) continue;
	    bool passMedium2015Id ((fabs(eleEta_p->at(elIter)) <= 1.4479
				&& fabs(eleSigmaIEtaIEta_p->at(elIter)) < 0.0101
				&& fabs(eledEtaAtVtx_p->at(elIter)) < 0.0103
				&& fabs(eledPhiAtVtx_p->at(elIter)) < 0.0336
				&& fabs(eleHoverE_p->at(elIter)) < 0.0876
				&& fabs(eleEoverP_p->at(elIter)) < 0.0174
				&& fabs(eleD0_p->at(elIter)) < 0.0118
				&& fabs(eleDz_p->at(elIter)) < 0.373
				&& fabs(eleMissHits_p->at(elIter)) <= 2 
				&& elepassConversionVeto_p->at(elIter)==1)
			       ||
			       (fabs(eleEta_p->at(elIter)) > 1.4479
				&& fabs(eleSigmaIEtaIEta_p->at(elIter)) < 0.0283
				&& fabs(eledEtaAtVtx_p->at(elIter)) < 0.00733
				&& fabs(eledPhiAtVtx_p->at(elIter)) < 0.114
				&& fabs(eleHoverE_p->at(elIter)) < 0.0678
				&& fabs(eleEoverP_p->at(elIter)) < 0.0898
				&& fabs(eleD0_p->at(elIter)) < 0.0739
				&& fabs(eleDz_p->at(elIter)) < 0.602
				&& fabs(eleMissHits_p->at(elIter)) <= 1 
				&& elepassConversionVeto_p->at(elIter)==1)
			       );
	    bool passVeto2015Id ((fabs(eleEta_p->at(elIter)) <= 1.4479
			      && fabs(eleSigmaIEtaIEta_p->at(elIter)) < 0.0114
			      && fabs(eledEtaAtVtx_p->at(elIter)) < 0.0152
			      && fabs(eledPhiAtVtx_p->at(elIter)) < 0.216
			      && fabs(eleHoverE_p->at(elIter)) < 0.181
			      && fabs(eleEoverP_p->at(elIter)) < 0.207
			      && fabs(eleD0_p->at(elIter)) < 0.0564
			      && fabs(eleDz_p->at(elIter)) < 0.472
			      && fabs(eleMissHits_p->at(elIter)) <= 2 
			      && fabs(elepassConversionVeto_p->at(elIter)))
			     ||
			     (fabs(eleEta_p->at(elIter)) > 1.4479
			      && fabs(eleSigmaIEtaIEta_p->at(elIter)) < 0.0352
			      && fabs(eledEtaAtVtx_p->at(elIter)) < 0.0113
			      && fabs(eledPhiAtVtx_p->at(elIter)) < 0.237
			      && fabs(eleHoverE_p->at(elIter)) < 0.116
			      && fabs(eleEoverP_p->at(elIter)) < 0.174
			      && fabs(eleD0_p->at(elIter)) < 0.222
			      && fabs(eleDz_p->at(elIter)) < 0.921
			      && fabs(eleMissHits_p->at(elIter)) <= 3 
			      && fabs(elepassConversionVeto_p->at(elIter)))
			     );
	    bool passMedium2016Id ((fabs(eleEta_p->at(elIter)) <= 1.4479
				    && fabs(eleSigmaIEtaIEta_p->at(elIter)) < 0.011
				    && fabs(eledEtaAtVtx_p->at(elIter)) < 0.0047
				    && fabs(eledPhiAtVtx_p->at(elIter)) < 0.222 
				    && fabs(eleHoverE_p->at(elIter)) < 0.298
				    && fabs(eleEoverP_p->at(elIter)) < 0.241
				    && fabs(eleD0_p->at(elIter)) < 0.05
				    && fabs(eleDz_p->at(elIter)) < 0.373
				    && fabs(eleMissHits_p->at(elIter)) <= 1
				    && elepassConversionVeto_p->at(elIter)==1)
				   ||
				   (fabs(eleEta_p->at(elIter)) > 1.4479
				    && fabs(eleSigmaIEtaIEta_p->at(elIter)) < 0.0314
				    && fabs(eledEtaAtVtx_p->at(elIter)) < 0.00868
				    && fabs(eledPhiAtVtx_p->at(elIter)) < 0.213 
				    && fabs(eleHoverE_p->at(elIter)) < 0.101
				    && fabs(eleEoverP_p->at(elIter)) < 0.14 
				    && fabs(eleD0_p->at(elIter)) < 0.10
				    && fabs(eleDz_p->at(elIter)) < 0.602
				    && fabs(eleMissHits_p->at(elIter)) <= 1 
				    && elepassConversionVeto_p->at(elIter)==1)
				   );
	    bool passVeto2016Id ((fabs(eleEta_p->at(elIter)) <= 1.4479
				  && fabs(eleSigmaIEtaIEta_p->at(elIter)) < 0.0115
				  && fabs(eledEtaAtVtx_p->at(elIter)) < 0.00749
				  && fabs(eledPhiAtVtx_p->at(elIter)) < 0.228
				  && fabs(eleHoverE_p->at(elIter)) <0.356 
				  && fabs(eleEoverP_p->at(elIter)) < 0.299
				  && fabs(eleD0_p->at(elIter)) < 0.0564
				  && fabs(eleDz_p->at(elIter)) < 0.472
				  && fabs(eleMissHits_p->at(elIter)) <= 2 
				  && fabs(elepassConversionVeto_p->at(elIter)))
				 ||
				 (fabs(eleEta_p->at(elIter)) > 1.4479
				  && fabs(eleSigmaIEtaIEta_p->at(elIter)) < 0.037 
				  && fabs(eledEtaAtVtx_p->at(elIter)) < 0.00895
				  && fabs(eledPhiAtVtx_p->at(elIter)) < 0.213 
				  && fabs(eleHoverE_p->at(elIter)) <0.211
				  && fabs(eleEoverP_p->at(elIter)) < 0.15
				  && fabs(eleD0_p->at(elIter)) < 0.222
				  && fabs(eleDz_p->at(elIter)) < 0.921
				  && fabs(eleMissHits_p->at(elIter)) <= 3 
				  && fabs(elepassConversionVeto_p->at(elIter)))
);
	    
	    double deposit, corrEA_isolation;
	    deposit =  fabs(elePFPhoIso_p->at(elIter)+elePFNeuIso_p->at(elIter)-eleEffAreaTimesRho_p->at(elIter));
	    corrEA_isolation = (elePFChIso_p->at(elIter) + TMath::Max (0.0, deposit )) / elePt_p->at(elIter);

	    bool passMediumIso( (corrEA_isolation < 0.0994 && fabs(eleEta_p->at(elIter)) <= 1.4479) || (corrEA_isolation < 0.107 && fabs(eleEta_p->at(elIter)) > 1.4479) );
	    bool passVetoIso( (corrEA_isolation < 0.175 && fabs(eleEta_p->at(elIter)) <= 1.4479) || (corrEA_isolation < 0.159 && fabs(eleEta_p->at(elIter)) > 1.4479) );

	    TString pf(fabs(eleEta_p->at(elIter)) > 1.4479 ? "_ee" : "_eb");
	    if(passMedium2016Id)
	      {
		histos["reliso"+pf]->Fill( corrEA_isolation,evWeight);
	      }
	    if(passMediumIso)
	      {
		histos["sigmaietaieta"+pf]->Fill(eleSigmaIEtaIEta_p->at(elIter),evWeight);
		histos["detain"+pf]->Fill(eledEtaAtVtx_p->at(elIter),evWeight);
		histos["dphiin"+pf]->Fill(eledPhiAtVtx_p->at(elIter),evWeight);
		histos["hovere"+pf]->Fill(eleHoverE_p->at(elIter),evWeight);
		histos["eoverpinv"+pf]->Fill(eleEoverP_p->at(elIter),evWeight);
		histos["d0"+pf]->Fill(eleD0_p->at(elIter),evWeight);
		histos["dz"+pf]->Fill(eleDz_p->at(elIter),evWeight);
		histos["misshits"+pf]->Fill(eleMissHits_p->at(elIter),evWeight);
		histos["convveto"+pf]->Fill(elepassConversionVeto_p->at(elIter),evWeight);
	      }
	    
	
	    //save electron if good
	    TLorentzVector p4(0,0,0,0);
	    p4.SetPtEtaPhiM(elePt_p->at(elIter),eleEta_p->at(elIter),elePhi_p->at(elIter), 0.0510);
	        
	    if (passMediumPt && passMediumEta && passMedium2016Id && corrEA_isolation > 0.2)
	      {
		mediumElectronsNonIso.push_back( p4 );
		elCharge.push_back(eleCharge_p->at(elIter));
	      }
	    else if(passMediumPt && passMediumEta && passMedium2016Id && passMediumIso)
	      {
		mediumElectrons.push_back( p4 );
		elCharge.push_back(eleCharge_p->at(elIter));
	      }
	    else if (passVeto2016Id && passVetoIso)
	      {
		vetoElectrons.push_back( p4 );
		elCharge.push_back(eleCharge_p->at(elIter));
	      }
	    else if (passMediumPt && passMediumEta && !passVeto2016Id)
              {
                mediumElectronsFailId.push_back( p4 );
		elCharge.push_back(eleCharge_p->at(elIter));
	      }

	  }
	


	//CHANNEL SELECTION
	//muons
	if(channelSelection==1300)
	  {
	    if(tightMuonsNonIso.size()==0) continue;
	    histos["recocounter"]->Fill(3);
	    if(tightMuons.size()+looseMuons.size()+mediumElectrons.size()+vetoElectrons.size()!=0) continue;
	    histos["recocounter"]->Fill(4);
	    tightMuons=tightMuonsNonIso;
	  }
	if(channelSelection==13)
	  {
            histos["ttevttype"]->Fill(ttevtype);
            bool fail2ndLepton(false);
            if(tightMuons.size()) histos["ttevttype"]->Fill(ttevtype+3);
	    if(tightMuons.size()!=1) fail2ndLepton=true;
	    if(!fail2ndLepton) histos["recocounter"]->Fill(3);
	    if(looseMuons.size()+mediumElectrons.size()+vetoElectrons.size()!=0) fail2ndLepton=true;
            if(fail2ndLepton)
              {
                histos["ttevttype"]->Fill(ttevtype+6);
                continue;
              }
	    histos["recocounter"]->Fill(4);
	    if(chargeSelection!=0)
	      {
		if(muonCharge[0]!=chargeSelection) continue;
	      }
	  }
	//electrons
	EffCorrection_t eselSF(1.0,0.0);
	if(channelSelection==1100)
	  {            
	    if(mediumElectronsFailId.size()==0) continue;
	     histos["recocounter"]->Fill(3);
	    if(mediumElectrons.size()+vetoElectrons.size()+tightMuons.size()+looseMuons.size()!=0) continue;
	    histos["recocounter"]->Fill(4);
	    mediumElectrons=mediumElectronsFailId;

	    //from Georgios studies, endcap electron fakes would still benefit from a ~15% extra weight
	    if(!isMC && fabs(mediumElectrons[0].Eta())>1.4479) evWeight *= 1.5;
	  }
	if(channelSelection==11)
	  {
            histos["ttevttype"]->Fill(ttevtype);
            bool fail2ndLepton(false);
            if(mediumElectrons.size()) histos["ttevttype"]->Fill(ttevtype+3);
            if(mediumElectrons.size()!=1) fail2ndLepton=true;
	    if(!fail2ndLepton) histos["recocounter"]->Fill(3);
	    if(vetoElectrons.size()+tightMuons.size()+looseMuons.size()!=0)fail2ndLepton=true;
            if(fail2ndLepton)
              {
                histos["ttevttype"]->Fill(ttevtype+6);
                continue;
              }
	    histos["recocounter"]->Fill(4);
	    if(chargeSelection!=0)
	      {
		if(elCharge[0]!=chargeSelection) continue;
	      }

	    //update event weight
	    eselSF=lepEffH.getOfflineCorrection(11,mediumElectrons[0].Pt(),mediumElectrons[0].Eta());
	    eselSF.second=sqrt(pow(0.03,2)+pow(eselSF.second,2));
	    evWeight*=eselSF.first;
	  }


	// combined leptons
	std::vector<TLorentzVector> goodLeptons;
	goodLeptons = (channelSelection==1300 || channelSelection==13) ?  tightMuons : mediumElectrons;  
	if(goodLeptons.size()==0) continue;

	//raw MET (noHF)
	TLorentzVector rawMET(0,0,0,0);	
	ljev.ntracks=0;
	ljev.ntracks_hp=0;
	for(size_t ipf=0; ipf<pfId_p->size(); ipf++)
	  {
	    Float_t abseta=TMath::Abs(pfEta_p->at(ipf));
	    if(abseta>3.0) continue;
	    rawMET += TLorentzVector(-pfPt_p->at(ipf)*TMath::Cos(pfPhi_p->at(ipf)),
				     -pfPt_p->at(ipf)*TMath::Sin(pfPhi_p->at(ipf)),
				     0,
				     0);	    
	  }
	
	//transverse mass
	float mt(computeMT(goodLeptons[0],rawMET));

	histos["lpt"]->Fill(goodLeptons[0].Pt(),evWeight);
	histos["leta"]->Fill(fabs(goodLeptons[0].Eta()),evWeight);
	if( (channelSelection==11 || channelSelection==1100) )
	  {
	    TString pf(fabs(goodLeptons[0].Eta()) > 1.4479 ? "_ee" : "_eb");
	    histos["lpt"+pf]->Fill(goodLeptons[0].Pt(),evWeight);
	    histos["leta"+pf]->Fill(fabs(goodLeptons[0].Eta()),evWeight);
	  }
	histos["mt"]->Fill(mt,evWeight);
	histos["metpt"]->Fill(rawMET.Pt(),evWeight);

	//jet counting
	typedef std::vector<TLorentzVector> JetColl_t;
	std::vector<JetColl_t> bJets(9),lightJets(9);
	for (Int_t jetIter = 0; jetIter < nref; jetIter++)
	  {
	    //cross clean with trigger muon
	    TLorentzVector jp4(0,0,0,0);
	    jp4.SetPtEtaPhiM(jtpt[jetIter],jteta[jetIter],jtphi[jetIter],jtm[jetIter]);

	    Int_t jflav(abs(refparton_flavor[jetIter]));	    

	    float UE_correction=0;
	    if (jteta[jetIter]>etaMin_hiFJRho_p->at(0) && jteta[jetIter]<etaMax_hiFJRho_p->at(0))
	      UE_correction= rhoCorr1Bin_hiFJRho_p->at(0);
	    if (jteta[jetIter]>etaMin_hiFJRho_p->at(1) && jteta[jetIter]<etaMax_hiFJRho_p->at(1))
	      UE_correction= rhoCorr1Bin_hiFJRho_p->at(1);
	    if (jteta[jetIter]>etaMin_hiFJRho_p->at(2) && jteta[jetIter]<etaMax_hiFJRho_p->at(2))
	      UE_correction= rhoCorr1Bin_hiFJRho_p->at(2);
	    if (jteta[jetIter]>etaMin_hiFJRho_p->at(3) && jteta[jetIter]<etaMax_hiFJRho_p->at(3))
	      UE_correction= rhoCorr1Bin_hiFJRho_p->at(3);
	    if (jteta[jetIter]>etaMin_hiFJRho_p->at(4) && jteta[jetIter]<etaMax_hiFJRho_p->at(4))
	      UE_correction= rhoCorr1Bin_hiFJRho_p->at(4);
	    if (jteta[jetIter]>etaMin_hiFJRho_p->at(5) && jteta[jetIter]<etaMax_hiFJRho_p->at(5))
	      UE_correction= rhoCorr1Bin_hiFJRho_p->at(5);
	    if (jteta[jetIter]>etaMin_hiFJRho_p->at(6) && jteta[jetIter]<etaMax_hiFJRho_p->at(6))
	      UE_correction= rhoCorr1Bin_hiFJRho_p->at(6);
	    if (jteta[jetIter]>etaMin_hiFJRho_p->at(7) && jteta[jetIter]<etaMax_hiFJRho_p->at(7))
	      UE_correction= rhoCorr1Bin_hiFJRho_p->at(7);
	    if (jteta[jetIter]>etaMin_hiFJRho_p->at(8) && jteta[jetIter]<etaMax_hiFJRho_p->at(8))
	      UE_correction= rhoCorr1Bin_hiFJRho_p->at(8);

	    //apply UE event subtraction for jet pT
	    if(!isMC)
	      jp4.SetPtEtaPhiM(jtpt[jetIter]-UE_correction*jtarea[jetIter],jteta[jetIter],jtphi[jetIter],jtm[jetIter]);
	    else
	      jp4.SetPtEtaPhiM(jtpt[jetIter],jteta[jetIter],jtphi[jetIter],jtm[jetIter]);
	    if(jp4.DeltaR(goodLeptons[0])<0.4) continue;
	    
	    //in tracker region
	    if(TMath::Abs(jp4.Eta())>2.4) continue;
            float drlj(jp4.DeltaR(goodLeptons[0]));
	    if(jp4.Pt()>JETPTTHRESHOLD)
              {
                if(jflav==5) histos["drlb"]->Fill(drlj);
                histos["drlj"]->Fill(drlj);
              }
	    if(drlj<0.4) continue;
	    
	    //systematic variations
	    bool passCSVL(discr_csvV2[jetIter]>0.460), passCSVM(discr_csvV2[jetIter]>0.8),passCSVMUp(passCSVM),passCSVMDn(passCSVM);	    
	    std::vector<float> jerSmear(3,1.0),jesScaleUnc(3,1.0);
	    if(isMC)
	      {
		//jet energy resolution smearing	
		//if(refpt[jetIter]>0) jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Eta(),refpt[jetIter]);
		//TLorentzVector rawjp4(jp4);
		//jp4 *= jerSmear[0];

		jesScaleUnc[1]=1.028;
		jesScaleUnc[2]=0.972;

		//b-tagging
		/*
		float jptforBtag(jp4.Pt()>1000. ? 999. : jp4.Pt());
		if(jflav==5)
		  {
		    float expEff    = expBtagEff["b"]->Eval(jptforBtag); 
		    myBTagSFUtil.modifyBTagsWithSF(passCSVMUp,1.1,expEff);	
		    myBTagSFUtil.modifyBTagsWithSF(passCSVMDn,0.9,expEff);	
		  }
		else if(jflav==4)
		  {
		    float expEff    = expBtagEff["c"]->Eval(jptforBtag); 
		    myBTagSFUtil.modifyBTagsWithSF(passCSVMUp,1.1,expEff);	
		    myBTagSFUtil.modifyBTagsWithSF(passCSVMDn,0.9,expEff);	
		  }
		else
		  {
		    float expEff    = expBtagEff["udsg"]->Eval(jptforBtag); 
		    myBTagSFUtil.modifyBTagsWithSF(passCSVMUp,1.3,expEff);	
		    myBTagSFUtil.modifyBTagsWithSF(passCSVMDn,0.7,expEff);	
		  }
		*/
	      }
	    
	    if(jp4.Pt()>JETPTTHRESHOLD)
	      {
		//nominal selection
		if(passCSVM) bJets[0].push_back(jp4);
		else         lightJets[0].push_back(jp4);

		if(ljev.nj<20)
		  {
		    ljev.j_pt[ljev.nj]  = jp4.Pt();
		    ljev.j_eta[ljev.nj] =jp4.Eta();
		    ljev.j_phi[ljev.nj] =jp4.Phi();
		    ljev.j_m[ljev.nj]   =jp4.M();
		    ljev.j_btag[ljev.nj]=passCSVL | (passCSVM<<1);
		    ljev.j_area[ljev.nj]= jtarea[jetIter];
		    ljev.j_refpt[ljev.nj]= refpt[jetIter];
		    ljev.j_refarea[ljev.nj]= refarea[jetIter];
		    ljev.j_refdr[ljev.nj]= refdrjt[jetIter]; 
		    ljev.j_refparton_flavorForB[ljev.nj]= refparton_flavorForB[jetIter];
		    ljev.j_refparton_flavor[ljev.nj]= refparton_flavor[jetIter];
		    ljev.j_PfCHF[ljev.nj]=jPfCHF[jetIter];
		    ljev.j_PfNHF[ljev.nj]=jPfNHF[jetIter];
		    ljev.j_PfCEF[ljev.nj]=jPfCEF[jetIter];
		    ljev.j_PfNEF[ljev.nj]=jPfNEF[jetIter];
		    ljev.j_PfMUF[ljev.nj]=jPfMUF[jetIter];
		    ljev.j_PfCHM[ljev.nj]=jPfCHM[jetIter];
		    ljev.j_PfNHM[ljev.nj]=jPfNHM[jetIter];
		    ljev.j_PfCEM[ljev.nj]=jPfCEM[jetIter];
		    ljev.j_PfNEM[ljev.nj]=jPfNEM[jetIter];
		    ljev.j_PfMUM[ljev.nj]=jPfMUM[jetIter];
		    ljev.nj++;
		  }

		//tag variations affect differently depending on the flavour
		if(jflav==5 || jflav==4)
		  {
		    if(passCSVMUp)
		      {
			bJets[1].push_back(jp4);
		      }
		    else
		      {
			lightJets[1].push_back(jp4);
		      }
		    if(passCSVMDn) 
		      {
			bJets[2].push_back(jp4);
		      }
		    else
		      {
			lightJets[2].push_back(jp4);
		      }
		    if(passCSVM)   
		      {
			bJets[3].push_back(jp4);
			bJets[4].push_back(jp4);
		      }
		    else      
		      {
			lightJets[3].push_back(jp4);
			lightJets[4].push_back(jp4);
		      }
		  }
		else
		  {
		    if(passCSVM)   
		      {
			bJets[1].push_back(jp4);
			bJets[2].push_back(jp4);
		      }
		    else      
		      {
			lightJets[1].push_back(jp4);
			lightJets[2].push_back(jp4);
		      }
		    if(passCSVMUp) 
		      {
			bJets[3].push_back(jp4);
		      }
		    else
		      {
			lightJets[3].push_back(jp4);
		      }
		    if(passCSVMDn) 
		      {
			bJets[4].push_back(jp4);
		      }
		    else
		      {
			lightJets[4].push_back(jp4);
		      }
		  }
	      }
	    
	    for(size_t ivar=0; ivar<2; ivar++)
	      {
		//JES varied selections
		TLorentzVector jesVarP4(jp4); jesVarP4*=jesScaleUnc[ivar+1];
		if(jesVarP4.Pt()>JETPTTHRESHOLD)
		  {
		    if(passCSVM) bJets[5+ivar].push_back(jesVarP4);
		    else         lightJets[5+ivar].push_back(jesVarP4);
		  }

		//JER varied selections
		TLorentzVector jerVarP4(jp4); jerVarP4*=jerSmear[ivar+1]/jerSmear[0];     
		if(jerVarP4.Pt()>JETPTTHRESHOLD)
		  {
		    if(passCSVM) bJets[7+ivar].push_back(jerVarP4);
		    else         lightJets[7+ivar].push_back(jerVarP4);
		  }
	      }
	  }

	//
	for(Int_t ivar=0; ivar<=nSysts; ivar++)
	  {
	    Int_t jetIdx(0);
	    if(ivar>=1 && ivar<=8) jetIdx=ivar;

	    //require at least two light jet acompanying the lepton
	    Int_t nljets(lightJets[jetIdx].size());
	    Int_t nbtags(bJets[jetIdx].size());
	    Int_t njets(nljets+nbtags);
	    
	    if(nljets<2) continue;
	    if (njets>=4)
	      histos["recocounter"]->Fill(5);
	    if (njets>=4&&nbtags>=1)
	      histos["recocounter"]->Fill(6);
	    if (njets>=4&&nbtags==1)
	      histos["recocounter"]->Fill(7);
	    if (njets>=4&&nbtags==2)
	      histos["recocounter"]->Fill(8);
	    TString pf(Form("%db",TMath::Min(nbtags,2)));
	    
	    //jet-related quantities
	    TLorentzVector jjp4((lightJets[jetIdx][0]+lightJets[jetIdx][1]));
	    Float_t mjj( jjp4.M() );
	    std::pair<int,int> jjLegsIdx=getDijetsSystemCandidate(lightJets[jetIdx]);
	    int idx1(jjLegsIdx.first),idx2(jjLegsIdx.second);
	    TLorentzVector rankedjjp4(lightJets[jetIdx][idx1]+lightJets[jetIdx][idx2]);
	    Float_t rankedmjj( rankedjjp4.M() );
	    Float_t drjj( lightJets[jetIdx][idx1].DeltaR( lightJets[jetIdx][idx2]) );
	    Float_t ptjj( rankedjjp4.Pt() );
	    Float_t etajj( rankedjjp4.Eta() );
	    Float_t htsum(0);
	    for(Int_t ij=0; ij<nbtags; ij++) htsum += bJets[jetIdx][ij].Pt();
	    for(Int_t ij=0; ij<nljets; ij++) htsum += lightJets[jetIdx][ij].Pt();


	    Float_t mbjj(-1),rankedmbjj(-1);
	    Float_t mlb( TMath::Min( (goodLeptons[0]+lightJets[jetIdx][idx1]).M(),
				     (goodLeptons[0]+lightJets[jetIdx][idx2]).M()) );
	    if(nbtags>0)
	      {
		mbjj=(jjp4+bJets[jetIdx][0]).M();
		rankedmbjj=(rankedjjp4+bJets[jetIdx][0]).M();
		mlb=(goodLeptons[0]+bJets[jetIdx][0]).M();
		if(nbtags>1)
		  {
		    mlb=TMath::Min( mlb, Float_t((goodLeptons[0]+bJets[jetIdx][1]).M()) );
		    
		    mbjj = jjp4.DeltaR(bJets[jetIdx][0])<jjp4.DeltaR(bJets[jetIdx][1]) ?
		      (rankedjjp4+bJets[jetIdx][0]).M():
		      (rankedjjp4+bJets[jetIdx][1]).M();
		    rankedmbjj = rankedjjp4.DeltaR(bJets[jetIdx][0])<rankedjjp4.DeltaR(bJets[jetIdx][1]) ? 
		      (rankedjjp4+bJets[jetIdx][0]).M() :
		      (rankedjjp4+bJets[jetIdx][1]).M();
		  }
	      }

	    //update event weight if needed
	    Float_t iweight(evWeight);
	    if(channelSelection==13 && ivar==9)  iweight*=1.03;
	    if(channelSelection==13 && ivar==10) iweight*=0.97;
	    if(channelSelection==11 && ivar==9)  iweight*=(1.0+eselSF.second);
	    if(channelSelection==11 && ivar==10) iweight*=(1.0-eselSF.second);

	    //fill histos and tree
	    if(ivar==0)
	      {
		ljev.run=run_;
		ljev.lumi=lumi_;
		ljev.event=evt_;

		//useful for UE dependecy
		ljev.hiHFplus=hiHFplus_;
		ljev.hiHFminus=hiHFminus_;
		ljev.hiHFplusEta4=hiHFplusEta4_;
		ljev.hiHFminusEta4=hiHFminusEta4_;

		//number of tracks and tracks overlapping with hard process objects
		ljev.ntracks=0;
		ljev.ntracks_hp=0;
		for(size_t ipf=0; ipf<pfId_p->size(); ipf++)
		  {

		    if(fabs(pfEta_p->at(ipf))>2.5 || pfPt_p->at(ipf)<0.9) continue;

		    int pfid(abs(pfId_p->at(ipf)));		    
		    if(pfid!=reco::PFCandidate::h && pfid!=reco::PFCandidate::e && pfid!=reco::PFCandidate::mu) continue;
		    
		    ljev.ntracks++;
		    TLorentzVector pfp4(0,0,0,0);
		    pfp4.SetPtEtaPhiE(pfPt_p->at(ipf),pfEta_p->at(ipf),pfPhi_p->at(ipf),pfEnergy_p->at(ipf));
		    bool overlapsWithHPobjects(false);
		    if(goodLeptons[0].DeltaR(pfp4)<0.01)              overlapsWithHPobjects=true;
		    if(lightJets[jetIdx][idx1].DeltaR(pfp4)<0.4)      overlapsWithHPobjects=true;
		    if(lightJets[jetIdx][idx2].DeltaR(pfp4)<0.4)      overlapsWithHPobjects=true;
		    if(nbtags>0 && bJets[jetIdx][0].DeltaR(pfp4)<0.4) overlapsWithHPobjects=true;
		    if(nbtags>1 && bJets[jetIdx][1].DeltaR(pfp4)<0.4) overlapsWithHPobjects=true;
		    if(!overlapsWithHPobjects) continue;
		    ljev.ntracks_hp++;
		  }

		ljev.w=iweight;
		ljev.l_id=channelSelection;
		ljev.l_pt=goodLeptons[0].Pt();
		ljev.l_eta=goodLeptons[0].Eta();
		ljev.l_phi=goodLeptons[0].Phi();
		ljev.l_m=goodLeptons[0].M();
		ljev.nb=nbtags;
		ljev.met_pt=rawMET.Pt();
		ljev.met_phi=rawMET.Phi();
		if(channelSelection==1300 || channelSelection==1100)
		  {
		    if(ljev.nj>2)	outT->Fill();
		  }
		else
		  {
		    outT->Fill();
		  }
		histos["mjj_"+pf]->Fill(mjj,iweight);
		histos["rankedmjj_"+pf]->Fill(rankedmjj,iweight);
		histos["mbjj_"+pf]->Fill(mbjj,iweight);
		histos["rankedmbjj_"+pf]->Fill(rankedmbjj,iweight);
		histos["drjj_"+pf]->Fill(drjj,iweight);
		histos["lpt_"+pf]->Fill(goodLeptons[0].Pt(),iweight);
		histos["leta_"+pf]->Fill(fabs(goodLeptons[0].Eta()),iweight);
		if( (channelSelection==11 || channelSelection==1100) )
		  {
		    TString region(fabs(goodLeptons[0].Eta()) > 1.4479 ? "ee_" : "eb_");
		    histos["lpt_"+region+pf]->Fill(goodLeptons[0].Pt(),evWeight);
		    histos["leta_"+region+pf]->Fill(fabs(goodLeptons[0].Eta()),evWeight);
		  }
		histos["ht_"+pf]->Fill(htsum,iweight);
		histos["mlb_"+pf]->Fill(mlb,iweight);
		histos["metpt_"+pf]->Fill(rawMET.Pt(),iweight);
		histos["metphi_"+pf]->Fill(rawMET.Phi(),iweight);
		histos["mt_"+pf]->Fill(mt,iweight);    
		if(nbtags)
		  {
		    histos["jpt_"+pf]->Fill(bJets[jetIdx][0].Pt(),iweight);
		    histos["jeta_"+pf]->Fill(fabs(bJets[jetIdx][0].Eta()),iweight);
		  }
		else
		  {
		    histos["jpt_"+pf]->Fill(lightJets[jetIdx][idx1].Pt(),iweight);
		    histos["jeta_"+pf]->Fill(fabs(lightJets[jetIdx][idx1].Eta()),iweight);
		    histos["jpt_"+pf]->Fill(lightJets[jetIdx][idx2].Pt(),iweight);
		    histos["jeta_"+pf]->Fill(fabs(lightJets[jetIdx][idx2].Eta()),iweight);
		  }
		histos["njets_"+pf]->Fill(njets,iweight);
		if(drjj<DRJJTHRESHOLD)
		  {
		    histos["rankedq70mjj_"+pf]->Fill(rankedmjj,iweight);
		    histos["ptjj_"+pf]->Fill(ptjj,iweight);
		    histos["etajj_"+pf]->Fill(etajj,iweight);
		  }
				  	       
		if(runSysts)
		  {
		    //theory uncertainties (by matrix-element weighting)
		    for(size_t igs=0; igs<ttbar_w_p->size(); igs++)
		      {
			float newWeight( iweight );
			if(isTTJets && normH && normH->GetBinContent(igs+1))
			  {
			    newWeight *= (ttbar_w_p->at(igs)/ttbar_w_p->at(0)) * ( normH->GetBinContent(1)/normH->GetBinContent(igs+1));
			  }
			else
			  {
			    newWeight *= (ttbar_w_p->at(igs)/ttbar_w_p->at(0));
			  }

			((TH2 *)histos["mjjshapes_"+pf+"_gen"])->Fill(mjj,igs,newWeight);
			((TH2 *)histos["drjjshapes_"+pf+"_gen"])->Fill(drjj,igs,newWeight);			
			((TH2 *)histos["rankedmjjshapes_"+pf+"_gen"])->Fill(rankedmjj,igs,newWeight);
			if(drjj<DRJJTHRESHOLD) ((TH2 *)histos["rankedq70mjjshapes_"+pf+"_gen"])->Fill(rankedmjj,igs,newWeight);
		      }
		  }
	      }
	    else if (runSysts)
	      {
		((TH2 *)histos["mjjshapes_"+pf+"_exp"])->Fill(mjj,ivar-1,iweight);
		((TH2 *)histos["drjjshapes_"+pf+"_exp"])->Fill(drjj,ivar-1,iweight);			
		((TH2 *)histos["rankedmjjshapes_"+pf+"_exp"])->Fill(rankedmjj,ivar-1,iweight);
		if(drjj<DRJJTHRESHOLD) 
		  ((TH2 *)histos["rankedq70mjjshapes_"+pf+"_exp"])->Fill(rankedmjj,ivar-1,iweight);
	      }
	  }
      }
    
    inFile_p->Close();
    delete inFile_p;    
  }
  
  //dump histograms
  outFile_p->cd();
  outT->Write();
  for(std::map<TString, TH1 *>::iterator it=histos.begin();
      it!=histos.end();
      it++)
    {
      it->second->SetDirectory(outFile_p);
      it->second->Write();
    }
  outFile_p->Close();
  delete outFile_p;
  
  return;
}
