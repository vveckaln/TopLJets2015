int checkerror()
{
  // const TString mcfilename = "$EOS/analysis_MC13TeV_TTJets/HADDChunks/MC13TeV_TTJets.root";
  const TString mcfilename = "/eos/user/v/vveckaln/analysis_MC13TeV_TTJets/Chunks/MC13TeV_TTJets/MC13TeV_TTJets_0.root";
  TFile * fmcttbar = TFile::Open(mcfilename.Data());
  TH1 * hmc = (TH1*) fmcttbar -> Get("L_pull_angle_allconst_reco_leading_jet_scnd_leading_jet_DeltaRTotal");
  hmc -> Scale(35874.8 * 832);
  for (unsigned char bind = 1; bind < 3; bind ++)
    {
      printf("%.9f error %.9f sqrt %.9f\n", hmc -> GetBinContent(bind), hmc -> GetBinError(bind), TMath::Sqrt(hmc -> GetBinContent(bind)));
    }
  printf("*******\n");
  // const TString datafilename = "/eos/user/v/vveckaln/analysis_MC13TeV_TTJets/Chunks/Data13TeV_SingleMuon_2016B/Data13TeV_SingleMuon_2016B_0.root";
  const TString datafilename = "/eos/user/v/vveckaln/analysis_MC13TeV_TTJets/HADDChunks/Data13TeV_SingleMuon_2016B.root";

  TFile * fdata = TFile::Open(datafilename);
  printf("%s\n", fdata -> GetName());
  //"$EOS/analysis_MC13TeV_TTJets/HADDChunks/Data13TeV_SingleMuon_2016B.root");
  TH1 * hdata = (TH1*) fdata -> Get("L_pull_angle_allconst_reco_leading_jet_scnd_leading_jet_DeltaRTotal");
  for (unsigned char bind = 1; bind < 3; bind ++)
    {
      printf("%f %f %f\n", hdata -> GetBinContent(bind), hdata -> GetBinError(bind), TMath::Sqrt(hdata -> GetBinContent(bind)));
    }
  printf("***************\n");
  TH1F h("h", "h", 5, 0.0, 5.0);
  //h.FillRandom("gaus", 5000);
  h.Sumw2();
  h.Fill(1.0, 1.0);
  h.Fill(3.0, 0.5);
  h.Fill(3.0, 0.5);
  for (unsigned char bind = 1; bind < h.GetNbinsX() + 1; bind ++)
    {
      printf("%f error %f sqrt %f\n", h.GetBinContent(bind), h . GetBinError(bind), TMath::Sqrt(h . GetBinContent(bind)));
    }

  printf("*** SCALING ***\n");
  h.Scale(0.5);
  // h.Fill(1.0, 1.0);
  // h.Fill(3.0, 0.5);
  // h.Fill(3.0, 0.5);
  for (unsigned char bind = 1; bind < h.GetNbinsX() + 1; bind ++)
    {
      printf(" %f %f %f\n", h.GetBinContent(bind), h . GetBinError(bind), TMath::Sqrt(h . GetBinContent(bind)));
    }

  return 0;
}
