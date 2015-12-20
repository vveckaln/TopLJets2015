#ifndef _ColourFlowAnalysisTool_hh_
#define _ColourFlowAnalysisTool_hh_
#include "TLorentzVector.h"
#include "TH1.h"
#include "TString.h"
#include <map>

class MiniEvent_t;

using namespace std;
class ColourFlowAnalysisTool
{
  struct PullVector: public TVector2
  {
    PullVector(float, float);
    Double_t & phi_component = fX;
    Double_t & eta_component = fY;
  };
public:
  
  const MiniEvent_t            * event_ptr_;
  const vector<TLorentzVector> * light_jets_ptr_;
  const vector<unsigned char>  * light_jets_indices_ptr_;
  const vector<TLorentzVector> * b_jets_ptr_;
  const vector<unsigned char>  * b_jets_indices_ptr_;
  map<TString, TH1*>     * plots_ptr_;
  float weight_;
  void AssignHistograms() const;
  void Work() const;
  float PullAngle(const PullVector & pull_vector, const TLorentzVector & jet_difference) const;
  PullVector CalculatePullVector(const TLorentzVector & jet, unsigned char index, bool OnlyChargedConstituents) const; 
  void Do();
};

#endif
