#ifndef _CFAT_Event_hh_
#define _CFAT_Event_hh_

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"

#include "TVector2.h"
#include "TLorentzVector.h"
using namespace std;

enum VectorEnum_t     {LEADING_JET, SCND_LEADING_JET, HAD_B, HAD_W, HAD_T, LEPTON, NEUTRINO, LEPT_B, LEPT_W, LEPT_T, BEAM, FAKE };
  
enum ChargeEnum_t     {ALLCOMP, CHARGED};
enum DeltaRCutEnum_t  {DELTAR_TOTAL, DELTAR_LE_1p0, DELTAR_GT_1p0};
enum PF_PTCutEnum_t   {PF_PT_TOTAL, PF_PT_LE_0p5_GEV, PF_PT_GT_0p5_GEV};
enum HadW_PtCutEnum_t {HADW_PT_TOTAL, HADW_PT_LE_50p0_GEV, HADW_PT_GT_50p0_GEV}; 
enum PFNEnum_t        {PFN_TOTAL, PFN_LE_20, PFN_GT_20};
enum PVMagEnum_t      {PVMAG_TOTAL, PVMAG_LE_0p005, PVMAG_GT_0p005};
enum WorkEnum_t       {RECO, GEN};
typedef unsigned char VectorCode_t;
typedef unsigned char ChargeCode_t;
typedef unsigned char WorkCode_t;
typedef unsigned char PF_PTCutCode_t;
class CFAT_Event;
struct PullVector: public TVector2
{
  PullVector(Double_t, Double_t);
  Double_t           & phi_component = fX;
  Double_t           & eta_component = fY;
  unsigned short     Ncomponents;
  VectorCode_t       origin_jet;
  const CFAT_Event * origin_event_ptr_;
};

class CFAT_Event
{

  friend class ColourFlowAnalysisTool;
  
  WorkCode_t                         work_mode_;
  struct PF_struct
  {
    unsigned int                     size;
    const int                        * jet_index;
    const int                        * id; 
    const int                        * charge;
    const float                      * pt;
    const float                      * eta; 
    const float                      * phi; 
    const float                      * mass;
                                     PF_struct();
    TLorentzVector                   GetPF(unsigned short index) const;
  } PF;
  
  double                           weight_;
  unsigned long                    event_number_;
  const MiniEvent_t              * event_ptr_;
  const TLorentzVector           * leading_light_jet_ptr_;
  unsigned short                   leading_light_jet_index_;

  const TLorentzVector           * second_leading_light_jet_ptr_;
  unsigned short                   second_leading_light_jet_index_;
  const TLorentzVector           * had_W_ptr_;
  const TLorentzVector           * had_b_ptr_;
  unsigned short                   had_b_index_;
  const TLorentzVector           * had_t_ptr_;
  
  const TLorentzVector           * lepton_ptr_;
  const TLorentzVector           * neutrino_ptr_;
  const TLorentzVector           * lept_W_ptr_;
  const TLorentzVector           * lept_b_ptr_;
  unsigned short                   lept_b_index_;
  const TLorentzVector           * lept_t_ptr_;
  const TLorentzVector           * fake_ptr_;
  static const TLorentzVector    * beam_ptr_;
  const TLorentzVector *&          GetVectorRef(VectorCode_t) ;
  void                             SetWorkMode(WorkCode_t = RECO);
  
public:
static const unsigned char       N_charge_types_ = 2;
  static const unsigned char       N_jet_types_ = 12;  
                                   CFAT_Event();
  void                             AddLightJets(const vector<TLorentzVector> &, const vector<unsigned short> &);
  void                             AddBJets(const vector<TLorentzVector> &, const vector<unsigned short> &);
  //void                             CompleteVectors();
  TLorentzVector                   GetChargedJet(VectorCode_t) const;
  void                             SetEvent(MiniEvent_t * ev);
  void                             SetWeight(double);
  void                             SetEventNumber(unsigned long);
  void                             AddVector(VectorCode_t, const TLorentzVector *);
  const TLorentzVector * const     GetVector(VectorCode_t) const;

  unsigned short                   GetIndex(VectorCode_t) const;
  double                           Angle(VectorCode_t, VectorCode_t) const;
  double                           DeltaR(VectorCode_t, VectorCode_t) const;
  double                           PullAngle(const PullVector &, VectorCode_t) const;
  PullVector                       CalculatePullVector(VectorCode_t, ChargeCode_t = ALLCOMP, PF_PTCutCode_t = PF_PT_TOTAL) const; 
  TVector2                         CalculatePullVectorEuclidian(const TLorentzVector & jet, unsigned char index, bool OnlyChargedConstituents) const;

  };

#endif
