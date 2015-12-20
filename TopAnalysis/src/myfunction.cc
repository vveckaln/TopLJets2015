#include "TopLJets2015/TopAnalysis/interface/myfunction.h"

#include "TLorentzVector.h"
#include "TVector2.h"
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include <sys/time.h>
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
/*
struct timeval_private:public timeval {
};
*/
void myfunction()
{
  printf("%lu\n", sizeof(MiniEvent_t));
  printf("%lu\n", sizeof(TH1F));
  printf("%lu\n", sizeof(TFile));
  printf("%lu\n", sizeof(double));
  printf("%lu\n", sizeof(float));
  TFile * f = TFile::Open("test_root.root", "RECREATE");
  TH1F * h = new TH1F("test", "test", 5, 0, 5);
  TH1F *h1 = new TH1F("test1", "test1", 5, 0, 5);
  h-> SetDirectory(0);
  h1 -> SetDirectory(0);
  h-> SetDirectory(f);
  h1-> SetDirectory(f);
  f -> Close();

}
