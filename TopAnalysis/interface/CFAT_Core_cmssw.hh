#ifndef _CFAT_Core_cmssw_hh_
#define _CFAT_Core_cmssw_hh_
#include "CFAT_Core.hh"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"

class CFAT_Core_cmssw;
class pf_cmssw: public pf
{
  friend class CFAT_Core_cmssw; 
  friend class CFAT_Core;
  CFAT_Core_cmssw * core_ptr_;
  unsigned short index_;
public: 
  virtual double GetCharge() const;
  virtual TLorentzVector GetLorentzVector() const;
  unsigned short GetJetConstituentIndex() const;
  unsigned short GetJetIndex() const;
  
  void SetCore(CFAT_Core_cmssw &);
  pf_cmssw();
};


class CFAT_Core_cmssw: public CFAT_Core
{
  const TLorentzVector           * leading_light_jet_ptr_;
  unsigned short                   leading_light_jet_index_;

  const TLorentzVector           * second_leading_light_jet_ptr_;
  unsigned short                   second_leading_light_jet_index_;
  const TLorentzVector           * had_b_ptr_;
  unsigned short                   had_b_index_;
  
  const TLorentzVector           * lepton_ptr_;
  const TLorentzVector           * neutrino_ptr_;
  const TLorentzVector           * lept_b_ptr_;
  unsigned short                   lept_b_index_;

  friend class pf_cmssw;
  friend class CFAT_Core;
  struct
  {
    unsigned short                     size;
    const Int_t                      * jet_index;
    const int                        * id; 
    const int                        * charge;
    const float                      * pt;
    const float                      * eta; 
    const float                      * phi; 
    const float                      * mass;
  } PF;
  MiniEvent_t * const GetEvent() const;

  class pf_iter_core: public CFAT_Core::pf_iter_core
  {
    friend class CFAT_Core_cmssw;
    friend class CFAT_Core;

  protected:
    pf_cmssw  particle;
    virtual pf       * operator -> ();
    pf_cmssw * GetPFCMSSW();
    pf_cmssw * const GetPFCMSSW() const;
   
    virtual bool       compare (const CFAT_Core::pf_iter_core *) ;
    pf_iter_core();
    virtual ~pf_iter_core();
  public:
    virtual void WhoAmI();
    virtual VectorCode_t   GetVectorCode() const;
  
  };
  static pf_iter_core * const END ; 

protected:
  virtual const TLorentzVector *& GetVectorRef(VectorCode_t code); 
  virtual        void increment(pf_iter &);

  unsigned short                   GetIndex(VectorCode_t) const;
  
  virtual const TLorentzVector * GetVector(VectorCode_t);
  virtual pf_iter begin(VectorCode_t);
  virtual pf_iter end(VectorCode_t);
  virtual pf_iter begin();
  virtual pf_iter end();
  virtual void SetWorkMode(WorkCode_t);
public:
  void SetEvent(const MiniEvent_t &);
  CFAT_Core_cmssw();
  void                             AddLightJets(const vector<TLorentzVector> &, const vector<unsigned short> &);
  void                             AddBJets(const vector<TLorentzVector> &, const vector<unsigned short> &);
  void                             AddVector(VectorCode_t, const TLorentzVector &);

  virtual void check(const char * = "") const;                                   

};

#endif
