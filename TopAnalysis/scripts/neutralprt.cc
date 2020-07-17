const unsigned char nfiles(16);
const float luminosity(35874.8);
const float xsec(832);
const char * datafiles[nfiles] = {"Data13TeV_SingleElectron_2016B.root",
			    "Data13TeV_SingleElectron_2016C.root",
			    "Data13TeV_SingleElectron_2016D.root",
			    "Data13TeV_SingleElectron_2016E.root",
			    "Data13TeV_SingleElectron_2016F.root",
			    "Data13TeV_SingleElectron_2016G.root",
			    "Data13TeV_SingleElectron_2016Hv2.root",
			    "Data13TeV_SingleElectron_2016Hv3.root",
			    "Data13TeV_SingleMuon_2016B.root",
			    "Data13TeV_SingleMuon_2016C.root",
			    "Data13TeV_SingleMuon_2016D.root",
			    "Data13TeV_SingleMuon_2016E.root",
			    "Data13TeV_SingleMuon_2016F.root",
			    "Data13TeV_SingleMuon_2016G.root",
			    "Data13TeV_SingleMuon_2016Hv2.root",
			    "Data13TeV_SingleMuon_2016Hv3.root"};
int neutralprt()
{
  TH1F * hall = nullptr;
  TH1F * hcharged = nullptr;
  for (unsigned char nfile = 0; nfile < nfiles; nfile ++)
    {
      TFile * f = new TFile(TString("$EOS/analysis_MC13TeV_TTJets/HADDChunks/") + datafiles[nfile]);
      if (nfile == 0)
	{
	  hall = (TH1F *) f -> Get("L_pull_angle_allconst_reco_leading_jet_scnd_leading_jet_DeltaRgt1p0") -> Clone();
	  hall -> SetDirectory(nullptr);
	  hcharged = (TH1F *) f -> Get("L_pull_angle_chconst_reco_leading_jet_scnd_leading_jet_DeltaRgt1p0") -> Clone();
	  hcharged -> SetDirectory(nullptr);
	}
      else
	{
	  hall -> Add((TH1F *) f -> Get("L_pull_angle_allconst_reco_leading_jet_scnd_leading_jet_DeltaRgt1p0"), 1.0);
	  hcharged -> Add((TH1F *) f -> Get("L_pull_angle_chconst_reco_leading_jet_scnd_leading_jet_DeltaRgt1p0"), 1.0);
	}
      f -> Close();
    }
  TH1F * hall_mc = nullptr;
  //  TH1F * hcharged_mc = nullptr;
  TFile * f = new TFile(TString("$EOS/analysis_MC13TeV_TTJets/HADDChunks/MC13TeV_TTJets.root"));
  hall_mc = (TH1F *) f -> Get("L_pull_angle_allconst_reco_leading_jet_scnd_leading_jet_DeltaRgt1p0") -> Clone();
  hall_mc -> Add((TH1F *) f -> Get("L_pull_angle_chconst_reco_leading_jet_scnd_leading_jet_DeltaRgt1p0"), -1);
  hall_mc -> SetLineColor(kRed);
  TCanvas * c = new TCanvas();
  hall -> Add(hcharged, -1);
  hall -> Scale(1.0/8.0);
  
  hall -> Rebin(8);
  hall -> SetMaximum(hall -> GetMaximum() * 1.2);
  hall -> SetMinimum(0.0);
  hall -> Draw("HIST");
  hall_mc -> Draw("SAMEHIST");
  hall_mc -> Rebin(8);
  hall_mc -> Scale(luminosity * xsec/8.0);
  gStyle -> SetOptTitle(kFALSE);
  TLegend * legend = new TLegend(0.1, 0.7, 0.3, 0.9);
  legend -> AddEntry(hall_mc, "mc t#bar{t}");
  legend -> AddEntry(hall, "data");
  legend -> Draw("SAME");
  c -> SaveAs("neutralprt.png");
  return 0;
}
