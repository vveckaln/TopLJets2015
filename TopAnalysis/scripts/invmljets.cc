int invmljets()
{
  TFile * f = TFile::Open("/eos/user/v/vveckaln/analysis_MC13TeV_TTJets/Chunks/MC13TeV_TTJets/MC13TeV_TTJets_0.root");
  // TH1 * h1 = (TH1 * ) f -> Get("L0_pre_invmljets");
  //h1 -> SetDirectory(nullptr);
  //invmleadingjets
  TH1 * h2 = (TH1 * ) f -> Get("L4_1l4j2b2w_invmleadingjets");
  printf("%p\n", h2);
  h2 -> SetDirectory(nullptr);
  h2 -> SetLineColor(kGreen);
  TH1 * h3 = (TH1 * ) f -> Get("L4_1l4j2b2w_invmwjets");
  printf("h3 %p\n", h3);
  h3 -> SetDirectory(nullptr);
  h3 -> SetLineColor(kRed);
  f -> Close();
  THStack * stack = new THStack;
  //  stack -> Add(h1);
  stack -> Add(h2, "HIST");
  stack -> Add(h3, "HIST");
  stack -> SetMinimum(0.0);
  TLegend * legend = new TLegend(0.6, 0.6, 0.95, 0.95);
  //  legend -> AddEntry(h1, h1 -> GetName(), "l");
  legend -> AddEntry(h2, h2 -> GetName(), "l");
  legend -> AddEntry(h3, "marked as W", "l");
  TCanvas *c = new TCanvas;
  stack -> Draw("nostack");
  legend -> Draw("SAME");
  c -> SaveAs("invmljets.png");
  return 0;
}
