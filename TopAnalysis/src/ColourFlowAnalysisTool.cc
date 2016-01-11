#include "TopLJets2015/TopAnalysis/interface/ColourFlowAnalysisTool.hh"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"

ColourFlowAnalysisTool::PullVector::PullVector(Double_t phi, Double_t eta): TVector2(phi, eta)
{
}

void ColourFlowAnalysisTool::AssignHistograms() const
{
  static const char* tag                         = "4j2t";
  static const char* tag2[2]                     = {"allconst", "chconst"};
  static const char* tag3[2]                     = {"q1q2", "q1b"};
  static const unsigned char ncategories         = 3;
  static const char * cat_title[ncategories]     = {"pull_angle",     "cos_pull_angle", "mag_pull_vector"  };
  static const unsigned short nbins[ncategories] = {100,              100,              100               };
  static const double min[ncategories]           = {-1.2*TMath::Pi(), -1.2,             0                 };
  static const double max[ncategories]           = {- min[0],         -min[1],          3                 };
  static const char* axes[ncategories]           = {"; rad; Events",  "; cos; Events",  "; a.u.; Events"  };
  for (unsigned char cat_index = 0; cat_index < ncategories; cat_index++)
    {
      for (unsigned char tag2_index = 0; tag2_index < 2; tag2_index ++)
	{
	  for (unsigned char tag3_index = 0; tag3_index < 2; tag3_index ++)
	    {
	      const TString hash_key = TString(cat_title[cat_index]) + "_" + tag2[tag2_index] + "_" + tag3[tag3_index] + "_" + tag;
	      plots_ptr_ -> operator[](hash_key) = new TH1F(hash_key, axes[cat_index], nbins[cat_index], min[cat_index], max[cat_index]);
	    }
    
	}
    }
  plots_ptr_ -> operator[]("sum_light_q_masses") = new TH1F("sum_light_q_masses", "sum of light quark masses; GeV; Events", 100, 0, 150);
}

void ColourFlowAnalysisTool::Work() const
{
  if( b_jets_ptr_ -> size() != 2 or light_jets_ptr_ -> size() != 2)
    return;
  plots_ptr_ -> operator[]("sum_light_q_masses") -> Fill((*light_jets_ptr_)[0].M() + (*light_jets_ptr_)[1].M(), weight_);
  const char* tag = "4j2t";
  const TLorentzVector * leading_jet = (*light_jets_ptr_)[0].Pt() >= (*light_jets_ptr_)[1].Pt() ? 
    &(*light_jets_ptr_)[0] : &(*light_jets_ptr_)[1];
  const unsigned char leading_jet_index = (*light_jets_ptr_)[0].Pt() >= (*light_jets_ptr_)[1].Pt() ? 
    (*light_jets_indices_ptr_)[0] : (*light_jets_indices_ptr_)[1];
  const TLorentzVector * second_leading_jet =  (*light_jets_ptr_)[0].Pt() >= (*light_jets_ptr_)[1].Pt() ? 
    &(*light_jets_ptr_)[1] : &(*light_jets_ptr_)[0];
  const TLorentzVector * b_jet = NULL;
  unsigned char b_jets_index;
  for (unsigned char index = 0; index < 2; index ++)
    {
      b_jets_index = b_jets_indices_ptr_ -> operator[](index);
      const char jet_flavour = event_ptr_ -> j_hadflav[b_jets_index];
      const char lepton_charge = event_ptr_ -> l_charge;
      if (lepton_charge * jet_flavour < 0) 
	continue;
      b_jet = &b_jets_ptr_ -> operator[](index);
    }

  static const char * tag3[2]  = {"allconst", "chconst"};
  static const bool OnlyChargedConstituents[2] = {false, true};
  for (unsigned char ind = 0; ind < 2; ind ++)
    {
      try
	{
	  const PullVector pull_vector = CalculatePullVector(*leading_jet, leading_jet_index, OnlyChargedConstituents[ind]);
	  {
	
	    const TVector2 jet_difference(TVector2::Phi_mpi_pi(second_leading_jet -> Phi() - leading_jet -> Phi()), second_leading_jet -> Eta() - leading_jet -> Eta());
	    const TString hash_key = TString("mag_pull_vector_") + tag3[ind] + "_q1q2_" + tag;
	    plots_ptr_ -> operator[](hash_key) -> Fill(pull_vector.Mod(), weight_); 
	    try
	      {
		const float pull_angle_q1q2 = PullAngle(pull_vector, jet_difference);
		plots_ptr_ -> operator[](TString("pull_angle_") + tag3[ind] + "_q1q2_" + tag) -> Fill(pull_angle_q1q2, weight_);
		const float cos_pull_angle_q1q2 = TMath::Cos(pull_angle_q1q2);
		plots_ptr_ -> operator[](TString("cos_pull_angle_") + tag3[ind] + "_q1q2_" + tag) -> Fill(cos_pull_angle_q1q2, weight_);
	      }catch (const char * e)
	      {
		printf("q1q2 %s\n", e);
	      }
      
	  }
	  if(b_jet)
	    {
	      const TVector2 jet_difference(TVector2::Phi_mpi_pi(b_jet -> Phi() - leading_jet -> Phi()), b_jet -> Eta() - leading_jet -> Eta());
	      try
		{
		  const float pull_angle_q1b = PullAngle(pull_vector, jet_difference);
		  plots_ptr_ -> operator[](TString("pull_angle_") + tag3[ind] + "_q1b_" + tag) -> Fill(pull_angle_q1b, weight_);
		  const float cos_pull_angle_q1b = TMath::Cos(pull_angle_q1b);
		  plots_ptr_ -> operator[](TString("cos_pull_angle_") + tag3[ind] + "_q1b_" + tag) -> Fill(cos_pull_angle_q1b, weight_);
		}catch (const char * e)
		{
		  printf("q1b %s\n", e);
		}
	    }
      
	}catch(const char *e)
	{

	  printf("%s\n", e);
	}
      }
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
  //  printf("%f %f phi %f\n", magnitude_pull, sqrt(pull_vector.phi_component*pull_vector.phi_component + pull_vector.eta_component * pull_vector.eta_component), pull_vector.phi_component);
  const float phi_dif = jet_difference.Px();
  //  printf("%f %f\n", phi_dif, jet_difference.Phi());
  //  getchar();
  const float eta_dif = jet_difference.Py();
  const float magnitude_dif = sqrt(phi_dif*phi_dif + eta_dif*eta_dif);
  float pull_angle = 0;
  if (magnitude_pull > 1E-6 and magnitude_dif > 1E-6)
    {
      const float cos_pullangle = (pull_vector.phi_component*phi_dif + pull_vector.eta_component*eta_dif)/
	(magnitude_pull * magnitude_dif);
      printf("cos angle %f\n", cos_pullangle);
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
