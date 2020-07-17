int compareeventyields()
{
  TFile * fv1 = TFile::Open("$EOS/eventnumbercheck/MC13TeV_TTJets_0_v1.root");
  TH1 * hv1 = (TH1*) fv1 -> Get("L_reco_selection");
  hv1 -> Scale(0.5);
  TFile * fv2 = TFile::Open("$EOS/eventnumbercheck/MC13TeV_TTJets_0_v2.root");
  TH1 * hv2 = (TH1*) fv2 -> Get("L_reco_selection");
  hv2 -> SetLineColor(kRed);
  THStack * stack = new THStack;
  stack -> Add(hv1);
  stack -> Add(hv2);
  TLegend * l = new TLegend(0.6, 0.6, 0.95, 0.95);
  l -> AddEntry(hv1, "old", "l");
  l -> AddEntry(hv2, "new", "l");
  
  TCanvas * c1 = new TCanvas;
  stack -> Draw("nostack");
  l -> Draw("same");
  c1 -> SaveAs("compareeventyields.png");
  return 0;
}
