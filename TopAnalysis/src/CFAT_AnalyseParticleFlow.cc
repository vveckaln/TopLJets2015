#include "TopLJets2015/TopAnalysis/interface/ColourFlowAnalysisTool.hh"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"


void ColourFlowAnalysisTool::AnalyseParticleFlow() const
{
  const TVector3 leading_jet_vect = leading_light_jet_ptr_ -> Vect();
  const TVector3 second_leading_jet_vect = second_leading_light_jet_ptr_ -> Vect();
  const TVector3 normal = leading_jet_vect.Cross(second_leading_jet_vect).Unit();
  const float angle = TMath::ACos(leading_jet_vect.Dot(second_leading_jet_vect)/(leading_jet_vect.Mag() * second_leading_jet_vect.Mag()));
  const float DeltaR = leading_jet_vect.DeltaR(second_leading_jet_vect);
  const char * DeltaR_tag = DeltaR <= 1.0 ? tag_DeltaR_types_[0] : tag_DeltaR_types_[1];
  
  static const bool OnlyChargedConstituents[2] = {false, true};
  unsigned int size = work_mode_ == 0 ? event_ptr_ -> npf : event_ptr_ -> ngen;
  
  for (unsigned int jet_const_index = 0; jet_const_index < size; jet_const_index ++)
    {
      for (unsigned char charge_index = 0; charge_index < 2; charge_index ++)
	{  
	  TLorentzVector constituent_4vector;
	  TVector3 constituent_vector;
	  bool const_jet_specific = true;
	  if (work_mode_ == 0)
	    {
	    
	      if (OnlyChargedConstituents[charge_index] and event_ptr_ -> pf_charge[jet_const_index] == 0)
		continue;
	       if (event_ptr_ -> pf_j[jet_const_index] != leading_light_jet_index_ and event_ptr_ -> pf_j[jet_const_index] != second_leading_light_jet_index_)
		 const_jet_specific = false;
	      const float jet_const_energy = sqrt(pow(event_ptr_ -> pf_px[jet_const_index], 2) +
						  pow(event_ptr_ -> pf_py[jet_const_index], 2) + 
						  pow(event_ptr_ -> pf_pz[jet_const_index], 2));
	      constituent_4vector = TLorentzVector    (event_ptr_ -> pf_px[jet_const_index], 
						       event_ptr_ -> pf_py[jet_const_index], 
						       event_ptr_ -> pf_pz[jet_const_index], 
						       jet_const_energy);
	    }
	  else
	    {
	      if (OnlyChargedConstituents[charge_index] and event_ptr_ -> g_charge[jet_const_index] == 0)
		continue;
	      if (event_ptr_ -> g_j[jet_const_index] != leading_light_jet_index_ and event_ptr_ -> g_j[jet_const_index] != second_leading_light_jet_index_)
		 const_jet_specific = false;
	      const float jet_const_energy = sqrt(pow(event_ptr_ -> g_px[jet_const_index], 2) +
						  pow(event_ptr_ -> g_py[jet_const_index], 2) + 
						  pow(event_ptr_ -> g_pz[jet_const_index], 2));
	      if (jet_const_energy == 0)
		continue;
	      constituent_4vector = TLorentzVector    (event_ptr_ -> g_px[jet_const_index], 
						       event_ptr_ -> g_py[jet_const_index], 
						       event_ptr_ -> g_pz[jet_const_index], 
						       jet_const_energy);

	    }
	  constituent_vector = constituent_4vector.Vect();

	  const TVector3 projection = normal.Cross(constituent_vector.Cross(normal));
	  try
	    {
	      if (projection.Mag() < 1E-6)
		throw "Zero projection";
	      
	      const float angle1 = TMath::ACos(leading_jet_vect.Dot(projection)/(leading_jet_vect.Mag() * projection.Mag()));
	      const float angle2 = TMath::ACos(second_leading_jet_vect.Dot(projection)/(second_leading_jet_vect.Mag() * projection.Mag()));
	      const TString prefix = TString("chi_") + tag_charge_types_[charge_index] + "_" + tag_levels_[work_mode_] + "_";
	      if (fabs(angle1 + angle2 - angle) < 1E-9)
		{
		  plots_ptr_ -> operator[](prefix + DeltaR_tag + "_all") -> Fill(angle1/angle, weight_);
		  plots_ptr_ -> operator[](prefix + tag_DeltaR_types_[2] + "_all") -> Fill(angle1/angle, weight_);
		  
		  if (const_jet_specific)
		    {
		      plots_ptr_ -> operator[](prefix + DeltaR_tag + "_jet") -> Fill(angle1/angle, weight_);
		      plots_ptr_ -> operator[](prefix + tag_DeltaR_types_[2] + "_jet") -> Fill(angle1/angle, weight_);
		 
		    }
		  
		  //		  printf("filling 1 %f \n", angle1/angle);
		}
	      else
		{
		  plots_ptr_ -> operator[](prefix + DeltaR_tag + "_all") -> Fill(angle1/(2*TMath::Pi() - angle) + 1, weight_);
		  plots_ptr_ -> operator[](prefix + tag_DeltaR_types_[2] + "_all") -> Fill(angle1/(2*TMath::Pi() - angle) + 1, weight_);

		  if (const_jet_specific)
		    {
		      plots_ptr_ -> operator[](prefix + DeltaR_tag + "_jet") -> Fill(angle1/(2*TMath::Pi() - angle) + 1, weight_);
		      plots_ptr_ -> operator[](prefix + tag_DeltaR_types_[2] + "_jet") -> Fill(angle1/(2*TMath::Pi() - angle) + 1, weight_);
		    }
		  //printf("filling 2 %f\n", angle1/(2*TMath::Pi() - angle));
		}
	    }
	  catch (const char *e)
	    {
	      printf("%s\n", e);
	    }
	}
    }
  
}
