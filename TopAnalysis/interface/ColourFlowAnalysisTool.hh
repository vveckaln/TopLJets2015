#ifndef _ColourFlowAnalysisTool_hh_
#define _ColourFlowAnalysisTool_hh_
#include "TLorentzVector.h"
#include "TH1.h"
#include "TString.h"
class TCanvas;
class TH2F;
#include <map>

extern TH1F * test;
extern TH1F * test2;
extern TH1F * eta;
extern TH1F * phi;

extern TH2F * test2D;


class MiniEvent_t;

using namespace std;
class ColourFlowAnalysisTool
{
  struct PullVector: public TVector2
  {
    PullVector(Double_t, Double_t);
    Double_t & phi_component = fX;
    Double_t & eta_component = fY;
  };
  static const unsigned char N_jet_types_;
  static const char* tag_channel_;
  static const char * tag_charge_types_[2];
  static const char* tag_jet_types_[];
  static const unsigned char N_DeltaR_types_;
  static const char * tag_DeltaR_types_[];
  static const char * tag_levels_[];
  static const TLorentzVector beam_;
  TCanvas * canvas_;
  TH2F    * PtRadProf_;
public:
  unsigned char                work_mode_;
  unsigned char                event_display_mode_;
  unsigned char                PtRadiation_mode_;
  unsigned long                event_number_;
  ColourFlowAnalysisTool();
  const MiniEvent_t            * event_ptr_;
  const vector<TLorentzVector> * light_jets_ptr_;
  const vector<unsigned short>  * light_jets_indices_ptr_;
  unsigned char                  jet_indices_[2] = {255, 255}; 
  const vector<TLorentzVector> * b_jets_ptr_;
  const vector<unsigned short>  * b_jets_indices_ptr_;
  const TLorentzVector         * neutrino_ptr_;
  const TLorentzVector         * lepton_ptr_;
  const TLorentzVector         * leading_light_jet_ptr_;
  vector<const TLorentzVector *> vect_jets_;
  float                          leading_light_jet_index_;
  float                          second_leading_light_jet_index_;
  const TLorentzVector         * second_leading_light_jet_ptr_;
  TLorentzVector had_t_;
  TLorentzVector lept_t_;
  
  map<TString, TH1*>     * plots_ptr_;
  //map<TString, TH1*>     & plots_ = *plots_ptr_;
  float weight_;
  TLorentzVector GetChargedJet(unsigned char jet_index) const;
  vector<const TLorentzVector *> IdentifyJets() ;
  
  void AssignHistograms();
  void Work();
  double PullAngle(const PullVector & pull_vector, const TVector2 & jet_difference) const;
  PullVector CalculatePullVector(const TLorentzVector & jet, unsigned char index, bool OnlyChargedConstituents) const; 
  TVector2 CalculatePullVectorEuclidian(const TLorentzVector & jet, unsigned char index, bool OnlyChargedConstituents) const;
  void PlotPUPPIWeight(unsigned char jet_vector_index, unsigned char jet_index, unsigned char charge_index) const;
  void PlotAngleBetweenJets() const;
  void PlotJetPhiAndEta() const;
  void Do();
  void AnalyseParticleFlow() const;
  void EventDisplay(const PullVector &, bool, const char *) const;
  void PtRadiationProfile() const;
  void ConfigurePtRadiationProfileTool();
  void EndPtRadiationProfile();
};

#endif
