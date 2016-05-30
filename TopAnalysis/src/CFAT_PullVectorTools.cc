#include "TopLJets2015/TopAnalysis/interface/ColourFlowAnalysisTool.hh"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"

#include "TH2F.h"

TH1F * test = NULL;
TH1F * test2 = NULL;
TH1F * eta = NULL;
TH1F * phi = NULL;
TH2F * test2D = NULL;
TLorentzVector ColourFlowAnalysisTool::GetChargedJet(unsigned char jet_index) const
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
}

double ColourFlowAnalysisTool::PullAngle(const PullVector & pull_vector, const TVector2 & jet_difference) const
{
  const double magnitude_pull = pull_vector.Mod();
  const double phi_dif = jet_difference.Px();
  const double eta_dif = jet_difference.Py();
  const double magnitude_dif = sqrt(phi_dif*phi_dif + eta_dif*eta_dif);

  double pull_angle = 0;
  /*const double nom = pull_vector.phi_OAcomponent*phi_dif + pull_vector.eta_component*eta_dif;
  const double denom = magnitude_pull * magnitude_dif;
  const double cos_angle = nom/denom;
  const double angle = TMath::ACos(cos_angle);
  const double angle_test = TMath::ACos(pull_vector.eta_component/magnitude_pull);
  */
  /*  if (fabs(pull_vector.phi_component) > 0.02)
    {
      printf("ABNORMAL %lu %f work_mode_ %u\n", event_number_, pull_vector.phi_component, work_mode_);
      EventDisplay(pull_vector, false, "ABNORMAL");
      getchar();
      throw "Abnromal phi component";
      }*/
  if (magnitude_pull > 1E-12 and magnitude_dif > 1E-12)
    {
      const double cos_pullangle = (pull_vector.phi_component*phi_dif + pull_vector.eta_component*eta_dif)/
	(magnitude_pull * magnitude_dif);
      pull_angle = TMath::ACos(cos_pullangle);
      if (pull_vector.eta_component - eta_dif < 0) 
	pull_angle *= -1;
      if (fabs(fabs(TMath::ACos(pull_vector.eta_component/magnitude_pull)) - TMath::Pi()/2.0) < 1E-1 )
	{
	  //printf("eta %f phi %f\n", pull_vector.eta_component, pull_vector.phi_component);
	  //getchar();
	}




      //printf("eta_dif %f nom %f denom %f cos_Angle %f angle %f angle test %f pull_angle %f\n", eta_dif, nom, denom, cos_angle, angle, angle_test, pull_angle);
      //getchar();
      /*test -> Fill(pull_angle, weight_);
      test2 -> Fill(- TMath::ACos(pull_vector.eta_component/magnitude_pull), weight_);
      test2D -> Fill(pull_vector.phi_component, pull_vector.eta_component, weight_);
      eta -> Fill(pull_vector.eta_component, weight_);
      phi -> Fill(pull_vector.phi_component, weight_);
      */
    }
  
  else throw "Null vector";
  return pull_angle;
}

 ColourFlowAnalysisTool::PullVector ColourFlowAnalysisTool::CalculatePullVector(const TLorentzVector & jet, unsigned char index, bool OnlyChargedConstituents) const
{
  //printf("calculate pull vector BEGIN \n");	     
  double phi_component = 0;
  double eta_component = 0;
  const double jet_phi = jet.Phi();
  const double jet_eta = jet.Rapidity();
  double Pt_jet_constituents = 0;
  unsigned int size = work_mode_ == 0 ? event_ptr_ -> npf : event_ptr_ -> ngpf;
  //  printf("index %u work_mode_ %u jet_rapidity %f jet_phi %f\n", index, work_mode_, jet_eta, jet_phi);
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
      const double delta_phi = TVector2::Phi_mpi_pi(constituent_4vector.Phi() - jet_phi);
      const double delta_eta = constituent_4vector.Rapidity() - jet_eta;
      const double mag = sqrt(delta_phi*delta_phi + delta_eta*delta_eta);
      phi_component += mag * delta_phi * constituent_4vector.Pt();
      //  printf("component_4vector.Rapidity() %f, constituent_4vector.Phi() %f delta phi %f \n", constituent_4vector.Rapidity(), constituent_4vector.Phi(), delta_phi);
      eta_component += mag * delta_eta * constituent_4vector.Pt();
		
    }
  if (Pt_jet_constituents < 1E-10)
    throw "Zero components phi eta calculation";
  const float scale = /*OnlyChargedConstituents ?*/ Pt_jet_constituents;// : */jet.Pt();
  phi_component /= scale;
  eta_component /= scale;
  //printf("Calculate pull vector END\n");
  return PullVector(phi_component, eta_component);
}


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
  const float scale = /*OnlyChargedConstituents ?*/ Pt_jet_constituents;// : */jet.Pt();
  x_component /= scale;
  y_component /= scale;
  return TVector2(x_component, y_component);
}
