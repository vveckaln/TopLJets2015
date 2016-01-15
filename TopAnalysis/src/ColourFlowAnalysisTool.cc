#include "TopLJets2015/TopAnalysis/interface/ColourFlowAnalysisTool.hh"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"

const char* ColourFlowAnalysisTool::tag_channel_              = "4j2t";
const char * ColourFlowAnalysisTool::tag_charge_type_[2]      = {"allconst", "chconst"};
const unsigned char ColourFlowAnalysisTool::Njettypes_        = 5;
const char* ColourFlowAnalysisTool::tag_jet_type_[Njettypes_] = {"q1q2", "q1b", "q1_leptb", "q1t", "q1_leptt"};

  

ColourFlowAnalysisTool::PullVector::PullVector(Double_t phi, Double_t eta): TVector2(phi, eta)
{
}

void ColourFlowAnalysisTool::AssignHistograms() const
{
  static const unsigned char ncategories_PA         = 2;
  static const char * cat_title_PA[ncategories_PA]     = {"pull_angle",     "cos_pull_angle" };
  static const unsigned short nbins_PA[ncategories_PA] = {100,              100              };
  static const double min_PA[ncategories_PA]           = {-1.2*TMath::Pi(), -1.2             };
  static const double max_PA[ncategories_PA]           = {- min_PA[0],         -min_PA[1]          };
  static const char* axes_PA[ncategories_PA]           = {"; rad; Events",  "; cos; Events"  };
  for (unsigned char cat_index = 0; cat_index < ncategories_PA; cat_index++)
    {
      for (unsigned char charge_index = 0; charge_index < 2; charge_index ++)
	{
	  for (unsigned char jet_index = 0; jet_index < Njettypes_; jet_index ++)
	    {
	      const TString hash_key = TString(cat_title_PA[cat_index]) + "_" + tag_charge_type_[charge_index] + "_" + tag_jet_type_[jet_index] + "_" + tag_channel_;
	      plots_ptr_ -> operator[](hash_key) = new TH1F(hash_key, axes_PA[cat_index], nbins_PA[cat_index], min_PA[cat_index], max_PA[cat_index]);
	    }
    
	}
    }
  static const unsigned char ncategories_PV = 3;
  static const char * cat_title_PV[ncategories_PV]     = {"phi_PV",         "eta_PV",         "mag_PV"      };
  static const unsigned short nbins_PV[ncategories_PV] = {100,              100,              100           };
  static const double min_PV[ncategories_PV]           = {-1.2*TMath::Pi(), -10,              7             };
  static const double max_PV[ncategories_PV]           = {- min_PV[0],         -min_PV[1],          0             };
  static const char* axes_PV[ncategories_PV]           = {"; rad; Events",  "; a.u.; Events", "a.u.; Events"};
  for (unsigned char cat_index = 0; cat_index < ncategories_PV; cat_index++)
    {
      for (unsigned char charge_index = 0; charge_index < 2; charge_index ++)
	{
	  const TString hash_key = TString(cat_title_PV[cat_index]) + "_" + tag_charge_type_[charge_index] + "_" + tag_channel_;
	  plots_ptr_ -> operator[](hash_key) = new TH1F(hash_key, axes_PV[cat_index], nbins_PV[cat_index], min_PV[cat_index], max_PV[cat_index]);
    
	}
    }
}

void ColourFlowAnalysisTool::Work()
{
  if( b_jets_ptr_ -> size() != 2 or light_jets_ptr_ -> size() != 2)
    return;
  const TLorentzVector * leading_jet = (*light_jets_ptr_)[0].Pt() >= (*light_jets_ptr_)[1].Pt() ? 
    &(*light_jets_ptr_)[0] : &(*light_jets_ptr_)[1];
  const unsigned char leading_jet_index = (*light_jets_ptr_)[0].Pt() >= (*light_jets_ptr_)[1].Pt() ? 
    (*light_jets_indices_ptr_)[0] : (*light_jets_indices_ptr_)[1];
  static const bool OnlyChargedConstituents[2] = {false, true};
  const TLorentzVector charged_jet = GetChargedJet(leading_jet_index);
  const TLorentzVector * jet[2] = {leading_jet, &charged_jet};
  const vector<const TLorentzVector *> other_jets = IdentifyOtherJets();
  for (unsigned char charge_ind = 0; charge_ind < 2; charge_ind ++)
    {
      try
	{
	  const PullVector pull_vector = CalculatePullVector(*jet[charge_ind], leading_jet_index, OnlyChargedConstituents[charge_ind]);
	  plots_ptr_ -> operator[](TString("phi_PV_") + tag_charge_type_[charge_ind] + "_" + tag_channel_) -> Fill(pull_vector.phi_component, weight_);
	  plots_ptr_ -> operator[](TString("eta_PV_") + tag_charge_type_[charge_ind] + "_" + tag_channel_) -> Fill(pull_vector.eta_component, weight_);
	  plots_ptr_ -> operator[](TString("mag_PV_") + tag_charge_type_[charge_ind] + "_" + tag_channel_) -> Fill(pull_vector.Mod(), weight_);
	  
	  
	  for (unsigned char jet_ind = 0; jet_ind < Njettypes_; jet_ind ++)
	    {
	      const TLorentzVector * other_jet = other_jets[jet_ind];
	      if (not other_jet)
		continue;
	      const TVector2 jet_difference(TVector2::Phi_mpi_pi(other_jet -> Phi() - jet[charge_ind] -> Phi()), other_jet -> Eta() - jet[charge_ind] -> Eta());
	      //	const TString hash_key = TString("mag_pull_vector_") + tag3[ind] + "_q1q2_" + tag;
	      //plots_ptr_ -> operator[](hash_key) -> Fill(pull_vector.Mod(), weight_); 
	      try
		{
		  const float pull_angle_q1q2 = PullAngle(pull_vector, jet_difference);
		  plots_ptr_ -> operator[](TString("pull_angle_") + tag_charge_type_[charge_ind] + "_" + tag_jet_type_[jet_ind] + "_" + tag_channel_) -> Fill(pull_angle_q1q2, weight_);
		  const float cos_pull_angle_q1q2 = TMath::Cos(pull_angle_q1q2);
		  plots_ptr_ -> operator[](TString("cos_pull_angle_") + tag_charge_type_[charge_ind] + "_" + tag_jet_type_[jet_ind] + "_" + tag_channel_)-> Fill(cos_pull_angle_q1q2, weight_);
		}catch (const char * e)
		{
		  printf("q1q2 %s\n", e);
		}
      
	      
	      
	    }
	}catch(const char *e)
	{

	  printf("%s\n", e);
	}
    }
}

vector<const TLorentzVector*> ColourFlowAnalysisTool::IdentifyOtherJets() 
{
  vector<const TLorentzVector*> other_jets;
  //second leading jet
  const TLorentzVector * second_leading_jet =  (*light_jets_ptr_)[0].Pt() >= (*light_jets_ptr_)[1].Pt() ? 
    &(*light_jets_ptr_)[1] : &(*light_jets_ptr_)[0];
  other_jets.push_back(second_leading_jet);
  //b jet
  const TLorentzVector * b_jet = NULL;
  const TLorentzVector * lept_bjet = NULL;
  const TLorentzVector Wboson = (*light_jets_ptr_)[0] + (*light_jets_ptr_)[1];
  const float t_mass = 173.34;
  float mass_dif = 5;
  float mass_dif_lept = 5;
  //unsigned char b_jets_index;
  for (unsigned char index = 0; index < 2; index ++)
    {
      //  b_jets_index = b_jets_indices_ptr_ -> operator[](index);
      float dif = fabs((b_jets_ptr_ -> operator[](index) + Wboson).M() - t_mass);
      if (dif < mass_dif)
	{
	  b_jet = &b_jets_ptr_ -> operator[](index);
	  mass_dif = dif;
	}
      if (dif > mass_dif_lept)
	{
	  lept_bjet = &b_jets_ptr_ -> operator[](index);
	  mass_dif_lept = dif;
	}
      /*const char jet_flavour = event_ptr_ -> j_hadflav[b_jets_index];
      const char lepton_charge = event_ptr_ -> l_charge;
      if (lepton_charge * jet_flavour < 0) 
      continue;
      b_jet = &b_jets_ptr_ -> operator[](index);*/
    }
  other_jets.push_back(b_jet);
  other_jets.push_back(lept_bjet);
  if (b_jet)
    {
      t_ = *b_jet + Wboson;
      other_jets.push_back(&t_);
    }
  else
    other_jets.push_back(NULL);
  if (lept_bjet)
    {
      lept_t_ = *lept_bjet + *lepton_ptr_ + *neutrino_ptr_;
      other_jets.push_back(&lept_t_);
    }
  else
    other_jets.push_back(NULL);
  return other_jets;

}

TLorentzVector ColourFlowAnalysisTool::GetChargedJet(unsigned char jet_index) const
{
  TLorentzVector charged_jet(0, 0, 0, 0);
  for (int jet_const_index = 0; jet_const_index < event_ptr_ -> npf; jet_const_index ++)
    {
      if (event_ptr_ -> pf_j[jet_const_index] != jet_index)
	continue;
      if (event_ptr_ -> pf_charge[jet_const_index] == 0)
	continue;
      const float jet_const_energy = sqrt(
					  event_ptr_ -> pf_px[jet_const_index]*event_ptr_ -> pf_px[jet_const_index]+
					  event_ptr_ -> pf_py[jet_const_index]*event_ptr_ -> pf_py[jet_const_index]+
					  event_ptr_ -> pf_pz[jet_const_index]*event_ptr_ -> pf_pz[jet_const_index]					                                                   );
      const TLorentzVector constituent_4vector(event_ptr_ -> pf_px[jet_const_index], 
					       event_ptr_ -> pf_py[jet_const_index], 
					       event_ptr_ -> pf_pz[jet_const_index], 
					       jet_const_energy);
      charged_jet += constituent_4vector;		
    }
  return charged_jet;
}

float ColourFlowAnalysisTool::PullAngle(const PullVector & pull_vector, const TVector2 & jet_difference) const
{
  const float magnitude_pull = pull_vector.Mod();
  const float phi_dif = jet_difference.Px();
  const float eta_dif = jet_difference.Py();
  const float magnitude_dif = sqrt(phi_dif*phi_dif + eta_dif*eta_dif);
  float pull_angle = 0;
  if (magnitude_pull > 1E-6 and magnitude_dif > 1E-6)
    {
      const float cos_pullangle = (pull_vector.phi_component*phi_dif + pull_vector.eta_component*eta_dif)/
	(magnitude_pull * magnitude_dif);
      pull_angle = TMath::ACos(cos_pullangle);
      if (pull_vector.eta_component - eta_dif < 0) 
	pull_angle *= -1;
    }
  else throw "Null vector";
  return pull_angle;
}

 ColourFlowAnalysisTool::PullVector ColourFlowAnalysisTool::CalculatePullVector(const TLorentzVector & jet, unsigned char index, bool OnlyChargedConstituents) const
{
  	     
  float phi_component = 0;
  float eta_component = 0;
  const float jet_phi = jet.Phi();
  const float jet_eta = jet.Eta();
  float Pt_jet_constituents = 0;
  for (int jet_const_index = 0; jet_const_index < event_ptr_ -> npf; jet_const_index ++)
    {
      if (event_ptr_ -> pf_j[jet_const_index] != index)
	continue;
      if (OnlyChargedConstituents and event_ptr_ -> pf_charge[jet_const_index] == 0)
	continue;
      const float jet_const_energy = sqrt(
					  event_ptr_ -> pf_px[jet_const_index]*event_ptr_ -> pf_px[jet_const_index]+
					  event_ptr_ -> pf_py[jet_const_index]*event_ptr_ -> pf_py[jet_const_index]+
					  event_ptr_ -> pf_pz[jet_const_index]*event_ptr_ -> pf_pz[jet_const_index]					                                                   );
      const TLorentzVector constituent_4vector(event_ptr_ -> pf_px[jet_const_index], 
					       event_ptr_ -> pf_py[jet_const_index], 
					       event_ptr_ -> pf_pz[jet_const_index], 
					       jet_const_energy);
      Pt_jet_constituents += constituent_4vector.Pt();
      const float delta_phi = TVector2::Phi_mpi_pi(constituent_4vector.Phi() - jet_phi);
      const float delta_eta = constituent_4vector.Eta() - jet_eta;
      const float mag = sqrt(delta_phi*delta_phi + delta_eta*delta_eta);
      phi_component += mag * delta_phi * constituent_4vector.Pt();
      eta_component += mag * delta_eta * constituent_4vector.Pt();
		
    }
  if (Pt_jet_constituents < 1E-6)
    throw "Zero components";
  const float scale = /*OnlyChargedConstituents ?*/ Pt_jet_constituents;// : */jet.Pt();
  phi_component /= scale;
  eta_component /= scale;
  return PullVector(phi_component, eta_component);
}

void ColourFlowAnalysisTool::Do()
{

}
