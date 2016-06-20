#include "TopLJets2015/TopAnalysis/interface/ColourFlowAnalysisTool.hh"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"

#include "TH2F.h"

TH1F * test = NULL;
TH1F * test2 = NULL;
TH1F * eta = NULL;
TH1F * phi = NULL;
TH2F * test2D = NULL;
/*TLorentzVector ColourFlowAnalysisTool::GetChargedJet(unsigned char jet_index) const
{
  TLorentzVector charged_jet(0, 0, 0, 0);
  if (work_mode_ == 0)
    for (int jet_const_index = 0; jet_const_index < event_ptr_ -> npf; jet_const_index ++)
      {
	if (event_ptr_ -> pf_j[jet_const_index] != jet_index)
	  continue;
	if (event_ptr_ -> pf_c[jet_const_index] == 0)
	  continue;
	TLorentzVector constituent_4vector;
	constituent_4vector.SetPtEtaPhiM(event_ptr_ -> pf_pt[jet_const_index], 
					 event_ptr_ -> pf_eta[jet_const_index], 
					 event_ptr_ -> pf_phi[jet_const_index], 
			 		 event_ptr_ -> pf_m[jet_const_index] 
					 );
	charged_jet += constituent_4vector;		
      }
  if (work_mode_ == 1)
    for (int jet_const_index = 0; jet_const_index < event_ptr_ -> ngpf; jet_const_index ++)
      {
	if (event_ptr_ -> gpf_g[jet_const_index] != jet_index)
	  continue;
	if (event_ptr_ -> gpf_c[jet_const_index] == 0)
	  continue;
	
	TLorentzVector constituent_4vector;
	constituent_4vector.SetPtEtaPhiM(event_ptr_ -> gpf_pt[jet_const_index], 
					 event_ptr_ -> gpf_eta[jet_const_index], 
					 event_ptr_ -> gpf_phi[jet_const_index], 
					 event_ptr_ -> gpf_m[jet_const_index] 
					 );
	charged_jet += constituent_4vector;		
      }
  return charged_jet;
  }*/
 /*
double ColourFlowAnalysisTool::PullAngle(const PullVector & pull_vector, const TVector2 & jet_difference) const
{
  const double magnitude_pull = pull_vector.Mod();
  const double phi_dif = jet_difference.Px();
  const double eta_dif = jet_difference.Py();
  const double magnitude_dif = sqrt(phi_dif*phi_dif + eta_dif*eta_dif);

  double pull_angle = 0;
  if (magnitude_pull > 1E-12 and magnitude_dif > 1E-12)
    {
      const double cos_pullangle = (pull_vector.phi_component*phi_dif + pull_vector.eta_component*eta_dif)/
	(magnitude_pull * magnitude_dif);
      pull_angle = TMath::ACos(cos_pullangle);
      if (pull_vector.eta_component - eta_dif < 0) 
	pull_angle *= -1;
    }
  
  else throw "Null vector";
  return pull_angle;
}*/


  /*
TVector2 ColourFlowAnalysisTool::CalculatePullVectorEuclidian(const TLorentzVector & jet, unsigned char index, bool OnlyChargedConstituents) const
{
  	     
  TVector2 pull_vector(0, 0);
  double x_component = 0;
  double y_component = 0;
  double Pt_jet_constituents = 0;
  unsigned int size = work_mode_ == 0 ? event_ptr_ -> npf : event_ptr_ -> ngpf;
  for (unsigned int jet_const_index = 0; jet_const_index < size; jet_const_index ++)
    {
      TLorentzVector constituent_4vector;
      if (work_mode_ == 0)
	{
	  if (event_ptr_ -> pf_j[jet_const_index] != index)
	    continue;
	  if (OnlyChargedConstituents and event_ptr_ -> pf_c[jet_const_index] == 0)
	    continue;
	  constituent_4vector.SetPtEtaPhiM(event_ptr_ -> pf_pt[jet_const_index], 
					   event_ptr_ -> pf_eta[jet_const_index], 
					   event_ptr_ -> pf_phi[jet_const_index], 
					   event_ptr_ -> pf_m[jet_const_index]);
	}
      else
	{
	  //index = event_ptr_ -> j_g[index];
	  if (event_ptr_ -> gpf_g[jet_const_index] != index)
	    continue;
	  if (OnlyChargedConstituents and event_ptr_ -> gpf_c[jet_const_index] == 0)
	    continue;
	  
	  constituent_4vector.SetPtEtaPhiM(event_ptr_ -> gpf_pt[jet_const_index], 
					   event_ptr_ -> gpf_eta[jet_const_index], 
					   event_ptr_ -> gpf_phi[jet_const_index], 
					   event_ptr_ -> gpf_m[jet_const_index]);

	}
      Pt_jet_constituents += constituent_4vector.Pt();

      const double mag = constituent_4vector.Pt();
      x_component += mag * (constituent_4vector.Px() - jet.Px()) * constituent_4vector.Pt();
      y_component += mag * (constituent_4vector.Py() - jet.Py()) * constituent_4vector.Pt();
		
    }
  if (Pt_jet_constituents < 1E-10)
    throw "Zero components Euclidian calculation";
  const float scale =   OnlyChargedConstituents ? Pt_jet_constituents;// :  et.Pt();
  x_component /= scale;
  y_component /= scale;
  return TVector2(x_component, y_component);
}*/
