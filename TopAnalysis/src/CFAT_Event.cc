#include "TopLJets2015/TopAnalysis/interface/CFAT_Event.hh"

//const int N_jet_types_ = 12;
//const unsigned char N_charge_types_ = 2;
const TLorentzVector * CFAT_Event::beam_ptr_ = new TLorentzVector();


PullVector::PullVector(Double_t phi, Double_t eta): TVector2(phi, eta)
{
  origin_jet = 255;
  Ncomponents = 0;
  origin_event_ptr_ = NULL;
}


CFAT_Event::PF_struct::PF_struct()
{
  size      = 0;
  jet_index = NULL;
  id        = NULL; 
  charge    = NULL;
  pt        = NULL;
  eta       = NULL; 
  phi       = NULL; 
  mass      = NULL;

}

TLorentzVector CFAT_Event::PF_struct::GetPF(unsigned short index) const
{
  TLorentzVector ret;
  ret.SetPtEtaPhiM(pt[index], eta[index], phi[index], mass[index]);
  return ret;
}

CFAT_Event::CFAT_Event()
{
  event_number_                    = 0;
  weight_                          = 0;
  work_mode_                       = RECO;
  event_ptr_                       = NULL;
  leading_light_jet_ptr_           = NULL;
  leading_light_jet_index_         = 0;

  second_leading_light_jet_ptr_    = NULL;
  second_leading_light_jet_index_  = 0;
  had_W_ptr_                       = NULL;
  had_b_ptr_                       = NULL;
  had_b_index_                     = 0;
  had_t_ptr_                       = NULL;
  
  lepton_ptr_                      = NULL;
  neutrino_ptr_                    = NULL;
  lept_W_ptr_                      = NULL;
  lept_b_ptr_                      = NULL;
  lept_b_index_                    = 0;
  lept_t_ptr_                      = NULL;
  fake_ptr_                        = NULL;

}

void CFAT_Event::SetEvent(MiniEvent_t * ev)
{
  event_ptr_ = ev;
}

void CFAT_Event::SetWeight(double w)
{
  weight_ = w;
}

void CFAT_Event::SetEventNumber(unsigned long numb)
{
  event_number_ = numb;
}

void CFAT_Event::AddLightJets(const vector<TLorentzVector> & light_jets, const vector<unsigned short> & light_jets_indices)
{
  const bool IsFirstLeading = light_jets[0].Pt() >= light_jets[1].Pt(); 
  leading_light_jet_ptr_ = IsFirstLeading ? 
    &light_jets[0] : &light_jets[1];
  leading_light_jet_index_ = IsFirstLeading ? 
    light_jets_indices[0] : light_jets_indices[1];
  second_leading_light_jet_ptr_ = IsFirstLeading ? 
    &light_jets[1] : &light_jets[0];
  second_leading_light_jet_index_ = IsFirstLeading ? 
    light_jets_indices[1] : light_jets_indices[0];
}

void CFAT_Event::AddBJets(const vector<TLorentzVector> & b_jets, const vector<unsigned short> & b_jets_indices)
{
  static TLorentzVector had_W(0.0, 0.0, 0.0, 0.0);
  had_W_ptr_ = &had_W;
  had_W = *leading_light_jet_ptr_ + *second_leading_light_jet_ptr_;
  static TLorentzVector lept_W(0.0, 0.0, 0.0, 0.0);
  lept_W_ptr_ = &lept_W;
  lept_W = *lepton_ptr_ + *neutrino_ptr_;
  static const float t_mass = 173.34;
  double mass_dif[2][2] = 
    {
      {1000, 1000}, 
      {1000, 1000}
    };
  unsigned char min_index[2] = {2, 2}; 
  for (unsigned char index = 0; index < 2; index ++)
    {
      mass_dif[index][0] = (b_jets[index] + had_W).M() - t_mass;
      mass_dif[index][1] = (b_jets[index] + lept_W).M() - t_mass;
    }
  for (unsigned char index = 0; index < 2; index ++)
    {
      min_index[index] = fabs(mass_dif[index][0]) >= fabs(mass_dif[index][1]) ? 1 : 0;
    }
  unsigned char had_b_local_index = 2;
  unsigned char lept_b_local_index = 2;
  if (min_index[0] != min_index[1])
    {
      had_b_local_index  = min_index[0] == 0 ? 0 : 1;
      lept_b_local_index = min_index[0] == 0 ? 1 : 0;
      had_b_ptr_         = &b_jets.at(had_b_local_index);
      lept_b_ptr_        = &b_jets.at(lept_b_local_index);
      had_b_index_       = b_jets_indices.at(had_b_local_index);
      lept_b_index_      = b_jets_indices.at(lept_b_local_index);
    }
  static TLorentzVector had_t(0.0, 0.0, 0.0, 0.0);
  if (had_b_ptr_)
    {
      had_t = *had_b_ptr_ + *had_W_ptr_;
      had_t_ptr_ = &had_t;
    }
  static TLorentzVector lept_t(0.0, 0.0, 0.0, 0.0);
  if (lept_b_ptr_)
    {
      lept_t = *lept_b_ptr_ + *lept_W_ptr_;
      lept_t_ptr_ = &lept_t;
    }
}

void CFAT_Event::AddVector(VectorCode_t code, const TLorentzVector * other)
{
  GetVectorRef(code) = other;

}


const TLorentzVector *& CFAT_Event::GetVectorRef(VectorCode_t code) 
{
  switch(code)
    {
    case LEADING_JET:
      return leading_light_jet_ptr_;
    case SCND_LEADING_JET:
      return second_leading_light_jet_ptr_;
    case HAD_B:
      return had_b_ptr_;
    case HAD_W:
      return had_W_ptr_;
    case HAD_T:
      return had_t_ptr_;
    case LEPTON:
      return lepton_ptr_;
    case NEUTRINO:
      return neutrino_ptr_;
    case LEPT_B:
      return lept_b_ptr_;
    case LEPT_W:
      return lept_W_ptr_;
    case LEPT_T:
      return lept_t_ptr_;
    case BEAM:
      return beam_ptr_;
    case FAKE:
      return fake_ptr_;
    default:
      throw "CFAT_Event::GetVectorRef(VectorCode_t) : check vector code";
    }

}

const TLorentzVector * const CFAT_Event::GetVector(VectorCode_t vector_code) const
{
  return (const TLorentzVector * const) const_cast< CFAT_Event * >(this) -> GetVectorRef(vector_code);
} 

unsigned short CFAT_Event::GetIndex(VectorCode_t vector_code) const
{
  unsigned int index = 32767;
  switch(vector_code)
    {
    case LEADING_JET:
      index = leading_light_jet_index_;
      break;
    case SCND_LEADING_JET:
      index = second_leading_light_jet_index_;
      break;
    case HAD_B:
      index = had_b_index_;
      break;
    case LEPT_B:
      index = lept_b_index_;
      break;
    default:
      char error[128];
      sprintf(error, "CFAT_Event::GetIndex(VectorCode_t) : Invalid vector code %u,", vector_code);
      throw error;
    }
  switch (work_mode_)
    {
    case GEN:
      return event_ptr_ -> j_g[index];
    default:
      return index;
    }
}

void CFAT_Event::SetWorkMode(WorkCode_t work_code)
{
  work_mode_ = work_code;
  switch (work_mode_)
    {
    case GEN:
      
      PF.size        = event_ptr_ -> ngpf;
      PF.jet_index   = event_ptr_ -> gpf_g;
      PF.id          = event_ptr_ -> gpf_id; 
      PF.charge      = event_ptr_ -> gpf_c;
      PF.pt          = event_ptr_ -> gpf_pt;
      PF.eta         = event_ptr_ -> gpf_eta; 
      PF.phi         = event_ptr_ -> gpf_phi; 
      PF.mass        = event_ptr_ -> gpf_m;
    
      break;
    default:
    
      PF.size        = event_ptr_ -> npf;
      PF.jet_index   = event_ptr_ -> pf_j;
      PF.id          = event_ptr_ -> pf_id; 
      PF.charge      = event_ptr_ -> pf_c;
      PF.pt          = event_ptr_ -> pf_pt;
      PF.eta         = event_ptr_ -> pf_eta; 
      PF.phi         = event_ptr_ -> pf_phi; 
      PF.mass        = event_ptr_ -> pf_m;
    
    }
}

double CFAT_Event::Angle(VectorCode_t code1, VectorCode_t code2) const
{
  const TLorentzVector * jet1 = GetVector(code1);
  const TLorentzVector * jet2 = GetVector(code2);
  if (not jet1 or not jet2)
    throw "CFAT_Event::Angle(VectorCode_t, VectorCode_t) : null vectors. Please Check!";
  if (jet1 == beam_ptr_)
    {
      if (jet2 == beam_ptr_)
	return 0.0;
      else
	return TMath::ACos(jet2 -> Pz() / jet2 -> Vect().Mag());
    }
  else
    { 
      if (jet2 == beam_ptr_)
	return TMath::ACos(jet1 -> Pz() / jet1 -> Vect().Mag());
      else
	return jet1 -> Angle(jet2 -> Vect());
    }
}

double CFAT_Event::DeltaR(VectorCode_t code1, VectorCode_t code2) const
{
  const TLorentzVector * jet1 = GetVector(code1);
  const TLorentzVector * jet2 = GetVector(code2);
  if (not jet1 or not jet2)
    throw "CFAT_Event::DeltaR(VectorCode_t, VectorCode_t) : null vectors. Please Check!";
  if (jet1 == beam_ptr_)
    {
      if (jet2 == beam_ptr_)
	return 0.0;
      else
	return 100.0;
    }
  else
    { 
      if (jet2 == beam_ptr_)
	return 100.0;
      else
	return jet1 -> DeltaR(*jet2);
    }
}

double CFAT_Event::PullAngle(const PullVector & pv, VectorCode_t code2) const
{
  const TLorentzVector * jet1 = GetVector(pv.origin_jet);
  const TLorentzVector * jet2 = GetVector(code2);
  if (not jet1 or not jet2)
    throw "CFAT_Event::PullAngle(const PullVector &, VectorCode_t) : null vectors. Please Check!";
 
      if (jet2 == beam_ptr_)
	return TMath::ACos(pv.eta_component / pv.Mod());
      else
	{
	  const TLorentzVector jet_difference = *jet2 - *jet1;
	  const TVector2 jet_difference_PhiEta(jet_difference.Phi(), jet_difference.Eta());
	  return fabs(TVector2::Phi_mpi_pi(jet_difference_PhiEta.Phi() - pv.Phi()));
    }
}

PullVector CFAT_Event::CalculatePullVector(VectorCode_t vector_code, ChargeCode_t charge_code, PF_PTCutCode_t pf_ptcut_code) const
{
  PullVector ret(0.0, 0.0);
  ret.origin_jet = vector_code;
  ret.origin_event_ptr_ = this;
  TLorentzVector charged_jet;
  if (charge_code == CHARGED)
    charged_jet = GetChargedJet(vector_code);
  const TLorentzVector * jet = charge_code == ALLCOMP ? GetVector(vector_code) : & charged_jet;
  const TVector2 jetV2(jet -> Phi(), jet -> Rapidity());
  double Pt_jet_constituents = 0.0;
  for (unsigned int jet_const_index = 0; jet_const_index < PF.size; jet_const_index ++)
    {
      TLorentzVector constituent_4vector ;
      if (PF.jet_index[jet_const_index] != GetIndex(vector_code))
	continue;
      if (charge_code == CHARGED and PF.charge[jet_const_index] == 0)
	continue;
      constituent_4vector = PF.GetPF(jet_const_index);
      if (pf_ptcut_code == PF_PT_LE_0p5_GEV and constituent_4vector.Pt() > 0.5)
	continue;
      if (pf_ptcut_code == PF_PT_GT_0p5_GEV and constituent_4vector.Pt() <= 0.5)
	continue;
      Pt_jet_constituents += constituent_4vector.Pt();
      const TVector2 componentV2 = TVector2(
					  TVector2::Phi_mpi_pi(constituent_4vector.Phi()), 
					  constituent_4vector.Rapidity())
	- jetV2;
      const double mag = componentV2.Mod();
      ret += mag * componentV2 * constituent_4vector.Pt();
      ret.Ncomponents ++;
    }
  if (Pt_jet_constituents < 1E-10)
    throw "PullVector CFAT_Event::CalculatePullVector(VectorCode_t, ChargeCode_t, PF_PTCutCode_t) const: Zero components";
  const double scale = /*OnlyChargedConstituents ?*/ Pt_jet_constituents;// : */jet.Pt();
  ret /= scale;
  return ret;
}

TLorentzVector CFAT_Event::GetChargedJet(VectorCode_t vector_code) const
{
  TLorentzVector charged_jet(0, 0, 0, 0);
  for (unsigned int jet_const_index = 0; jet_const_index < PF.size; jet_const_index ++)
    {
      TLorentzVector constituent_4vector ;
      if (PF.jet_index[jet_const_index] != GetIndex(vector_code))
	continue;
      if (PF.charge[jet_const_index] == 0)
	continue;
      constituent_4vector = PF.GetPF(jet_const_index);
      charged_jet += constituent_4vector;
      
    }
    
  return charged_jet;
}

