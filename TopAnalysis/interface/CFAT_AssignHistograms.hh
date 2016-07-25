#ifndef _CFAT_AssignHistograms_hh_
#define _CFAT_AssignHistograms_hh_
#include "TString.h"
#include "TH1F.h"
#include <map>
using namespace std;
typedef /*double */map<TString, TH1*> mapp;
void AssignHistograms(mapp &);
void CreateHistogram1D(mapp &, const char * title, const char * axes, unsigned int nbins, double min , double max);
//void sum(double a, int b);
#endif
