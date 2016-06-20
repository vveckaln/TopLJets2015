#ifndef _ColourFlowAnalysisTool_hh_
#define _ColourFlowAnalysisTool_hh_
#include "TopLJets2015/TopAnalysis/interface/CFAT_Event.hh"

#include "TLorentzVector.h"
#include "TH1.h"
#include "TString.h"
class TCanvas;
#include "TH2F.h"
#include <map>

extern TH1F * test;
extern TH1F * test2;
extern TH1F * eta;
extern TH1F * phi;

extern TH2F * test2D;

using namespace std;

class ColourFlowAnalysisTool
{
  static const char          * tag_channel_;
  static const char          * tag_charge_types_[2];
  static const char          * tag_jet_types_[];
  static const unsigned char   N_DeltaR_types_;
  
  static const char          * tag_DeltaR_types_[];
  static const unsigned char   N_PF_Pt_cuts_ = 3;
  static const char          * PF_Pt_cuts_types_[];// = {"PFPtTotal", "PFPtle0p5GeV", "PFPtgt0p5GeV"};

  static const unsigned char   N_HadW_Pt_cuts_ = 3;
  static const char          * HadW_Pt_cuts_types_[];// = {"hadWPtTotal", "hadWPtle50p0GeV", "hadWPtGt50p0GeV"}; 
  static const unsigned char   N_PF_N_cuts_ = 3;
  static const char          * PF_N_cuts_types_[];// = {"PFN_Total", "PFNle20", "PFNgt20"};
  static const unsigned char   N_PVMag_cuts_ = 3;
  static const char *          PVMag_cuts_types_[];// = { "PVMag_Total", "PVMag_le_0p005", "PVMag_gt_0p005"};


static const char          * tag_levels_[];
  TCanvas * canvas_;
  TH2F    * PtRadProf_;
  WorkCode_t work_mode_;
public:
  CFAT_Event * cfat_event_ptr_;
  
  unsigned char event_display_mode_;
  unsigned char                PtRadiation_mode_;
  ColourFlowAnalysisTool();
  
  map<TString, TH1*>     * plots_ptr_;
  map<TString, TH2*>     * plots2D_ptr_;
  void SetWorkMode(WorkCode_t);
  void AssignHistograms();
  void Work();
  void Fill1D(const TString &, const double) const;
  void CreateHistogram1D(const char *, const char *, unsigned int, double, double) const;
  void PlotAngleBetweenJets() const;
  void PlotJetDimensions() const;
  void Do();
  void AnalyseParticleFlow() const;
  void EventDisplay(const PullVector &, bool, const char *) const;
  void PtRadiationProfile() const;
  void ConfigurePtRadiationProfileTool();
  void EndPtRadiationProfile();
};

#endif
