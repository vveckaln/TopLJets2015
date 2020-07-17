int checkeventselection()
{
  TFile *f = TFile::Open("/eos/user/v/vveckaln/analysis_MC13TeV_TTJets/plots/plotter_L_pull_angle_allconst_reco_scnd_leading_jet_leading_jet_DeltaRTotal.root", "READ");
  TDirectory * dir = f -> GetDirectory("L_reco_selection");
  dir -> ls();
  TH1 * h = (TH1 *) dir -> Get("L_reco_selection");
  h -> SetDirectory(nullptr);
  // h -> Print("all");
  for (unsigned char bind = 1; bind < h -> GetNbinsX() + 1; bind ++)
    {
      printf("%u, %f %f\n", bind, h -> GetBinContent(bind), TMath::Log(h -> GetBinContent(bind))/TMath::Ln10());
    }
  TCanvas * c = new TCanvas;
  h -> Draw("HIST");
  h -> SetMaximum(20 * h -> GetMaximum());
  h -> SetMinimum(1);
  c -> SetLogy();
  f -> Close();
  return 0;
}
