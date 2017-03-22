#include "CFAT_AssignHistograms_raw.hh"
#include "TopLJets2015/TopAnalysis/interface/CFAT_AssignHistograms.hh"

void AssignHistograms(map_hist_1D & target)
{
  AssignHistograms_raw<map_hist_1D>(target, CreateHistogram1D);
}
void AssignSpecificHistograms2D(map<TString, TH2*> & target)
{
  CreateHistogram2D(target, "Migration", "Migration : rec : gen", 100, -TMath::Pi(), TMath::Pi(), 30, -TMath::Pi(), TMath::Pi());
}
void CreateHistogram1D(map_hist_1D & target, const char * title, const char * axes, unsigned int nbins, double min , double max)
{
  target[title] = new TH1F(title, TString(title) + axes, nbins, min, max);
}
void CreateHistogram2D(map<TString, TH2*> & target, const char * title, const char * axes, unsigned int nbins_x, double min_x , double max_x, unsigned int nbins_y, double min_y , double max_y)
{
  target[title] = new TH2F(title, TString(title) + axes, nbins_x, min_x, max_x, nbins_y, min_y, max_y);
}

