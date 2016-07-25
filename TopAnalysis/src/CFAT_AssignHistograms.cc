#include "CFAT_AssignHistograms_raw.hh"
#include "TopLJets2015/TopAnalysis/interface/CFAT_AssignHistograms.hh"

void AssignHistograms(mapp & target)
{
  //sum<double, int>( a, b);
  AssignHistograms_raw<mapp>(target, CreateHistogram1D);
}

void CreateHistogram1D(mapp & target, const char * title, const char * axes, unsigned int nbins, double min , double max)
{
  target[title] = new TH1F(title, TString(title) + axes, nbins, min, max);
}

/*void sum(double a, int b)
{
  printf("%f\n", a + b);
}
*/
