const unsigned char nsamples(3);
const char * samples[nsamples] = {"sherpa", "herwig7", "dire2002"};
const Color_t scolor[nsamples] = {kBlue, kGreen, kMagenta};
int dire()
{
  TFile * fgen[nsamples];
  TH1 * hgen[nsamples];
  THStack * stack = new THStack;
  TLegend * legend = new TLegend(0.6, 0.6, 1.0, 1.0);
  for (unsigned char sind = 0; sind < nsamples; sind ++)
  {
    // if (TString(samples[sind]) == "sherpa")
    //   continue;
    fgen[sind] = TFile::Open(TString("$EOS/analysis_MC13TeV_TTJets") + "_" + samples[sind] + "/HADDChunks/MC13TeV_TTJets_" + samples[sind] + ".root");
    hgen[sind] = (TH1 *) fgen[sind] -> Get("L_pull_angle_allconst_gen_leading_jet_scnd_leading_jet_DeltaRTotal");
    hgen[sind] -> SetDirectory(nullptr);
    fgen[sind] -> Close();
    hgen[sind] -> SetLineColor(scolor[sind]);
    printf("%s entries %f\n", samples[sind], hgen[sind] -> GetEntries());
    stack -> Add(hgen[sind]);
    legend -> AddEntry(hgen[sind], samples[sind], "l");
  }
  TFile * fttbar = TFile::Open("$EOS/analysis_MC13TeV_TTJets/plots/plotter.root");
  TH1 * httbar = (TH1 *) fttbar -> GetDirectory("L_pull_angle_allconst_gen_leading_jet_scnd_leading_jet_DeltaRTotal") 
    -> Get("L_pull_angle_allconst_gen_leading_jet_scnd_leading_jet_DeltaRTotal_t#bar{t}");
  printf("ttbar entries %f \n", httbar -> GetEntries());
  httbar -> SetDirectory(nullptr);
  fttbar -> Close();
  httbar -> SetLineColor(kRed);
  stack -> Add(httbar);
  legend -> AddEntry(httbar, "Pythia8", "l");
  TCanvas * canvas = new TCanvas;
  stack -> Draw("nostack");
  legend -> Draw("SAME");
  //stack -> SetMinimum(0.0);
  stack -> SetMaximum(stack -> GetMaximum("nostack") * 1.2);
  gPad -> SetLogy();
  canvas -> SaveAs("GEN.png");
  return 0;
}
