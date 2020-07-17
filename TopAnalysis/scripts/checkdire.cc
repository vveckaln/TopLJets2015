const unsigned char nsamples(3);
const char * samples[nsamples] = {"sherpa", "herwig7", "dire2002"};
const Color_t scolor[nsamples] = {kBlue, kGreen, kMagenta};
int checkdire()
{
  const unsigned char sind(2);
  TFile * f = TFile::Open(TString("$EOS/analysis_MC13TeV_TTJets") + "_" + samples[sind] + "/HADDChunks/MC13TeV_TTJets_" + samples[sind] + ".root");
  TH1 * h = (TH1*) f -> Get("L_pull_angle_allconst_gen_leading_jet_scnd_leading_jet_DeltaRgt1p0");
  h -> SetDirectory(nullptr);
  f -> Close();
  TChain chain("migration");
  chain.Add(TString("$EOS/analysis_MC13TeV_TTJets") + "_" + samples[sind] + "/HADDmigration/migration_MC13TeV_TTJets_" + samples[sind] + ".root/M_allconst_leading_jet_migration");
  chain.Add(TString("$EOS/analysis_MC13TeV_TTJets") + "_" + samples[sind] + "/HADDmigration/migration_MC13TeV_TTJets_" + samples[sind] + ".root/E_allconst_leading_jet_migration");
  Float_t reco;
  Float_t gen;
  vector<double> * weights = nullptr;
  const char * observable = "pull_angle";
  chain.SetBranchAddress(TString(observable) + "_reco", & reco);
  chain.SetBranchAddress(TString(observable) + "_gen",  & gen);
  chain.SetBranchAddress("weight",          & weights);
  const unsigned int dbcut = 1;
  TH2F * hm = new TH2F("m", "m", h -> GetNbinsX(), h -> GetXaxis() -> GetXmin(), h -> GetXaxis() -> GetXmax(), h -> GetNbinsX(), h -> GetXaxis() -> GetXmin(), h -> GetXaxis() -> GetXmax());
  //      h -> FillFromTree(, tag_jet_types_[jetcode], tag_charge_types_[chargecode], observable);                                                                                                      
  for (unsigned long event_ind = 0; event_ind < chain.GetEntries()/dbcut; event_ind ++)
    {
      if (event_ind % 10000 == 0)
	printf("%u\r", event_ind);
      chain.GetEntry(event_ind);
      if (string(observable).compare("pull_angle") == 0)
	hm -> Fill(reco, gen, (*weights)[0]);
    }
  TCanvas * c = new  TCanvas;
  h -> SetMinimum(0.0);
  h -> Draw("HIST");
  c -> Modified();
  c-> Update();
  getchar();
  TH1D * hmproj = hm -> ProjectionY();
  hmproj -> SetLineColor(kRed);
  hmproj -> Draw("HISTSAME");
  return 0;
}
