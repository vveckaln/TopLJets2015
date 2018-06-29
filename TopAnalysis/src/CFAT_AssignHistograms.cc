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
  for (unsigned short channel_ind = 0; channel_ind < N_channels_types_; channel_ind ++)
    {
      printf("creatin %s\n", (TString(tag_channels_types_[channel_ind]) + "_" + title).Data());
      target[TString(tag_channels_types_[channel_ind]) + "_" + title] = new TH1F(TString(tag_channels_types_[channel_ind]) + "_" + title, TString(title) + axes, nbins, min, max);
    }
}
void CreateHistogram2D(map<TString, TH2*> & target, const char * title, const char * axes, unsigned int nbins_x, double min_x , double max_x, unsigned int nbins_y, double min_y , double max_y)
{
  for (unsigned short channel_ind = 0; channel_ind < N_channels_types_; channel_ind ++)
    {
      target[TString(tag_channels_types_[channel_ind]) + "_" + title] = new TH2F(TString(tag_channels_types_[channel_ind]) + "_" + title, TString(title) + axes, nbins_x, min_x, max_x, nbins_y, min_y, max_y);
    }
}

