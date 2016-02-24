#include "TopLJets2015/TopAnalysis/interface/ColourFlowAnalysisTool.hh"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
ColourFlowAnalysisTool::PullVector::PullVector(Double_t phi, Double_t eta): TVector2(phi, eta)
{
}

ColourFlowAnalysisTool::ColourFlowAnalysisTool()
{
  event_display_mode_ = 0;
  canvas_ = 0;
  PtRadProf_ = 0;
  PtRadiation_mode_ = 0;
}

void ColourFlowAnalysisTool::Work()
{
  leading_light_jet_ptr_ = NULL;
  second_leading_light_jet_ptr_ = NULL;
  if( b_jets_ptr_ -> size() != 2 or light_jets_ptr_ -> size() != 2)
    return;
  static const bool OnlyChargedConstituents[2] = {false, true};
  
  vect_jets_ = IdentifyJets();
  PlotAngleBetweenJets();   
  for (unsigned char jet1_index = 0; jet1_index < 2; jet1_index ++)
    {
      const TLorentzVector charged_jet = GetChargedJet(jet_indices_[jet1_index]);
      const TLorentzVector * jet1_array[2] = {vect_jets_[jet1_index], &charged_jet};
      for (unsigned char charge_index = 0; charge_index < 2; charge_index ++)
	{
	  const TLorentzVector * jet1 = jet1_array[charge_index];
	  try
	    {
	      const PullVector pull_vector = CalculatePullVector(*jet1, jet_indices_[jet1_index], OnlyChargedConstituents[charge_index]);
	      const TString suffix =  TString("_") + tag_charge_types_[charge_index] + "_" + 
		tag_levels_[work_mode_] + "_" + 
		tag_jet_types_[jet1_index] + "_" + 
		tag_channel_;
	      plots_ptr_ -> operator[](TString("phi_PV") + suffix) -> Fill(pull_vector.phi_component, weight_);
	      plots_ptr_ -> operator[](TString("eta_PV") + suffix) -> Fill(pull_vector.eta_component, weight_);
	      plots_ptr_ -> operator[](TString("mag_PV") + suffix) -> Fill(pull_vector.Mod(), weight_);
	  
	  
	      for (unsigned char jet2_index = 0; jet2_index < N_jet_types_; jet2_index ++)
		{
		  if (jet2_index == jet1_index)
		    continue;
		  const TLorentzVector * jet2 = vect_jets_[jet2_index];
		  if (not jet2)
		    continue;
		  //		  printf("jet1 %u jet2 %u\n", jet1_index, jet2_index);
		  float DeltaR = 0;
		  if (jet2_index != 6)
		    DeltaR = jet1 -> DeltaR(*jet2);
		  else
		    {
		      const float Theta = jet2 -> Theta();
		      const float Eta = -0.5 * TMath::Log(TMath::Tan(Theta/2));
		      DeltaR = (pow(jet1 -> Phi() - jet2 -> Phi(), 2) + pow(jet1 -> Eta() - Eta, 2));

		    }
		  const char * DeltaR_tag = DeltaR <= 1.0 ? tag_DeltaR_types_[0] : tag_DeltaR_types_[1];
		  const TVector2 jet_difference(/*TVector2::Phi_mpi_pi*/(jet2 -> Phi() - jet1 -> Phi()), jet2 -> Rapidity() - jet1 -> Rapidity());
		  if (jet1_index == 0 and jet2_index == 1 )
		    {
		      const TString suffix =  TString("_") + tag_charge_types_[charge_index] + "_" + 
			tag_levels_[work_mode_] + "_" + 
			tag_jet_types_[jet1_index] + "_" + 
			tag_channel_;
		      plots_ptr_ -> operator[](TString("phi_jet_dif") + suffix) -> Fill(jet_difference.Px(), weight_);
		      plots_ptr_ -> operator[](TString("eta_jet_dif") + suffix) -> Fill(jet_difference.Py(), weight_);
		      plots_ptr_ -> operator[](TString("mag_jet_dif") + suffix) -> Fill(jet_difference.Mod(), weight_);
		    }
		  try
		    {
		      const float pull_angle = PullAngle(pull_vector, jet_difference);
		      static unsigned short count = 0;

		      if (count < 60 )
			if (jet1_index == 0 and jet2_index == 1)
			  {
			    //	    event_display_mode_ = 1;
			    EventDisplay(pull_vector, pull_angle, jet_difference, OnlyChargedConstituents[charge_index]);
			    count ++;
			    //  getchar();
			  }
		      
		      const float cos_pull_angle = TMath::Cos(pull_angle);
		      const TString infix = TString("_") + tag_charge_types_[charge_index] + "_" + 
			tag_levels_[work_mode_] + "_" + 
			tag_jet_types_[jet1_index] + "_:_" + 
			tag_jet_types_[jet2_index] + "_"; 
		      
		      plots_ptr_ -> operator[](TString("pull_angle")     + infix + DeltaR_tag       + "_" + tag_channel_) 
			-> Fill(pull_angle, weight_);
		      plots_ptr_ -> operator[](TString("cos_pull_angle") + infix + DeltaR_tag       + "_" + tag_channel_) 
			-> Fill(cos_pull_angle, weight_);
		      plots_ptr_ -> operator[](TString("pull_angle")     + infix + tag_DeltaR_types_[2] + "_" + tag_channel_) 
			-> Fill(pull_angle, weight_);
		      plots_ptr_ -> operator[](TString("cos_pull_angle") + infix + tag_DeltaR_types_[2] + "_" + tag_channel_) 
			-> Fill(cos_pull_angle, weight_);
		      
		    }
		  catch (const char * e)
		    {
		      printf("%s\n", e);
		    }
		}
	    }
	  catch(const char *e)
	    {

	      printf("%s\n", e);
	    }
	}
    }
  AnalyseParticleFlow();
  PtRadiationProfile();
  
}

vector<const TLorentzVector*> ColourFlowAnalysisTool::IdentifyJets() 
{
  vector<const TLorentzVector*> vect_jets;
  vect_jets.reserve(N_jet_types_);
  const TLorentzVector * leading_jet =  (*light_jets_ptr_)[0].Pt() >= (*light_jets_ptr_)[1].Pt() ? 
    &(*light_jets_ptr_)[0] : &(*light_jets_ptr_)[1];
  const unsigned char leading_jet_index =  (*light_jets_ptr_)[0].Pt() >= (*light_jets_ptr_)[1].Pt() ? 
    (*light_jets_indices_ptr_)[0] : (*light_jets_indices_ptr_)[1];
  leading_light_jet_index_ = leading_jet_index;
  jet_indices_[0] = leading_jet_index;
  if (work_mode_ == 0)
    plots_ptr_ -> operator[]("leading_jet_flavour") -> Fill(event_ptr_ -> j_hadflav[leading_jet_index], weight_);
  vect_jets.push_back(leading_jet);

  //second leading jet
  const TLorentzVector * second_leading_jet =  (*light_jets_ptr_)[0].Pt() >= (*light_jets_ptr_)[1].Pt() ? 
    &(*light_jets_ptr_)[1] : &(*light_jets_ptr_)[0];
  const unsigned char second_leading_jet_index =  (*light_jets_ptr_)[0].Pt() >= (*light_jets_ptr_)[1].Pt() ? 
    (*light_jets_indices_ptr_)[1] : (*light_jets_indices_ptr_)[0];
  jet_indices_[1] = second_leading_jet_index;
  leading_light_jet_ptr_ = leading_jet;
  second_leading_light_jet_ptr_ = second_leading_jet;
  second_leading_light_jet_index_ = second_leading_jet_index;

  if (work_mode_ == 0)
    plots_ptr_ -> operator[]("second_leading_jet_flavour") -> Fill(event_ptr_ -> j_hadflav[second_leading_jet_index], weight_);
  vect_jets.push_back(second_leading_jet);
  //b jet
  const TLorentzVector had_W_boson = (*light_jets_ptr_)[0] + (*light_jets_ptr_)[1];
  const TLorentzVector lept_W_boson = *lepton_ptr_ + *neutrino_ptr_;
  plots_ptr_ -> operator[](TString("had_") + tag_levels_[work_mode_] + "_W_mass") -> Fill(had_W_boson.M(), weight_);
  plots_ptr_ -> operator[](TString("lept_") + tag_levels_[work_mode_] + "_W_mass") -> Fill(lept_W_boson.M(), weight_);
  
  const float t_mass = 173.34;
  float mass_dif[2][2] = {{1000, 1000}, {1000, 1000}};
  unsigned char min_index[2] = {2, 2}; 
  for (unsigned char index = 0; index < 2; index ++)
    {
      mass_dif[index][0] = (b_jets_ptr_ -> operator[](index) + had_W_boson).M() - t_mass;
      mass_dif[index][1] = (b_jets_ptr_ -> operator[](index) + lept_W_boson).M() - t_mass;
    }
  for (unsigned char index = 0; index < 2; index ++)
    {
      min_index[index] = fabs(mass_dif[index][0]) >= fabs(mass_dif[index][1]) ? 1 : 0;
    }
  unsigned char had_b_jet_local_index = 2;
  unsigned char lept_b_jet_local_index = 2;
  const TLorentzVector * had_b_jet = NULL;
  const TLorentzVector * lept_b_jet = NULL; 

  if (min_index[0] != min_index[1])
    {
      had_b_jet_local_index = min_index[0] == 0 ? 0 : 1;
      lept_b_jet_local_index = min_index[0] == 0 ? 1 : 0;
      had_b_jet = &b_jets_ptr_ -> at(had_b_jet_local_index);
      lept_b_jet = &b_jets_ptr_ -> at(lept_b_jet_local_index);
      const unsigned char had_b_jet_index = b_jets_indices_ptr_ -> at(had_b_jet_local_index);
      const unsigned char lept_b_jet_index = b_jets_indices_ptr_ -> at(lept_b_jet_local_index);
      if (work_mode_ == 0)
	{
	  plots_ptr_ -> operator[]("had_b_flavour") -> Fill(event_ptr_ -> j_hadflav[had_b_jet_index], weight_);
	  plots_ptr_ -> operator[]("lept_b_flavour") -> Fill(event_ptr_ -> j_hadflav[lept_b_jet_index], weight_);;
	}
      /*const char jet_flavour = event_ptr_ -> j_hadflav[b_jets_index];
      const char lepton_charge = event_ptr_ -> l_charge;
      if (lepton_charge * jet_flavour < 0) 
      continue;
      b_jet = &b_jets_ptr_ -> operator[](index);*/
    }
  vect_jets.push_back(had_b_jet);
  vect_jets.push_back(lept_b_jet);
  if (had_b_jet)
    {
      had_t_ = *had_b_jet + had_W_boson;
      plots_ptr_ -> operator[](TString("had_") + tag_levels_[work_mode_] + "_t_mass") -> Fill(had_t_.M(), weight_); 
      vect_jets.push_back(&had_t_);
    }
  else
    vect_jets.push_back(NULL);
  if (lept_b_jet)
    {
      lept_t_ = *lept_b_jet + lept_W_boson;
      plots_ptr_ -> operator[](TString("lept_") + tag_levels_[work_mode_] + "_t_mass") -> Fill(lept_t_.M(), weight_); 
     
      vect_jets.push_back(&lept_t_);
    }
  else
    vect_jets.push_back(NULL);
  vect_jets.push_back(&beam_);
  return vect_jets;

}

TLorentzVector ColourFlowAnalysisTool::GetChargedJet(unsigned char jet_index) const
{
  TLorentzVector charged_jet(0, 0, 0, 0);
  if (work_mode_ == 0)
    for (int jet_const_index = 0; jet_const_index < event_ptr_ -> npf; jet_const_index ++)
      {
	if (event_ptr_ -> pf_j[jet_const_index] != jet_index)
	  continue;
	if (event_ptr_ -> pf_charge[jet_const_index] == 0)
	  continue;
	const float jet_const_energy = sqrt(pow(event_ptr_ -> pf_px[jet_const_index], 2) +
					    pow(event_ptr_ -> pf_py[jet_const_index], 2) + 
					    pow(event_ptr_ -> pf_pz[jet_const_index], 2));
	const TLorentzVector constituent_4vector(event_ptr_ -> pf_px[jet_const_index], 
						 event_ptr_ -> pf_py[jet_const_index], 
						 event_ptr_ -> pf_pz[jet_const_index], 
						 jet_const_energy);
	charged_jet += constituent_4vector;		
      }
  if (work_mode_ == 1)
    for (int jet_const_index = 0; jet_const_index < event_ptr_ -> ngen; jet_const_index ++)
      {
	if (event_ptr_ -> g_j[jet_const_index] != jet_index)
	  continue;
	if (event_ptr_ -> g_charge[jet_const_index] == 0)
	  continue;
	const float jet_const_energy = sqrt(pow(event_ptr_ -> g_px[jet_const_index], 2) +
					    pow(event_ptr_ -> g_py[jet_const_index], 2) +
					    pow(event_ptr_ -> g_pz[jet_const_index], 2));
	const TLorentzVector constituent_4vector(event_ptr_ -> g_px[jet_const_index], 
						 event_ptr_ -> g_py[jet_const_index], 
						 event_ptr_ -> g_pz[jet_const_index], 
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
  const float jet_eta = jet.Rapidity();
  float Pt_jet_constituents = 0;
   unsigned int size = work_mode_ == 0 ? event_ptr_ -> npf : event_ptr_ -> ngen;
  for (unsigned int jet_const_index = 0; jet_const_index < size; jet_const_index ++)
    {
      TLorentzVector constituent_4vector;
      if (work_mode_ == 0)
	{
	  if (event_ptr_ -> pf_j[jet_const_index] != index)
	    continue;
	  if (OnlyChargedConstituents and event_ptr_ -> pf_charge[jet_const_index] == 0)
	    continue;
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
	  if (event_ptr_ -> g_j[jet_const_index] != index)
	    continue;
	  if (OnlyChargedConstituents and event_ptr_ -> g_charge[jet_const_index] == 0)
	    continue;
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
      Pt_jet_constituents += constituent_4vector.Pt();
      const float delta_phi = /*TVector2::Phi_mpi_pi*/(constituent_4vector.Phi() - jet_phi);
      const float delta_eta = constituent_4vector.Rapidity() - jet_eta;
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



void ColourFlowAnalysisTool::PlotAngleBetweenJets() const
{
  for (unsigned char jet1_index = 0; jet1_index < vect_jets_.size(); jet1_index ++)
    {
      const TLorentzVector * jet1 = vect_jets_[jet1_index];
      if (not jet1)
	continue;
      for (unsigned char jet2_index = jet1_index + 1; jet2_index < vect_jets_.size(); jet2_index ++)
	{
	  const TLorentzVector * jet2 = vect_jets_[jet2_index];
	  if (not jet2)
	    continue;
	  float DeltaR = 0;
	  if (jet2_index != 6)
	    DeltaR = jet1 -> DeltaR(*jet2);
	  else
	    {
	      const float Theta = jet2 -> Theta();
	      const float Eta = -0.5 * TMath::Log(TMath::Tan(Theta/2));
	      DeltaR = (pow(jet1 -> Phi() - jet2 -> Phi(), 2) + pow(jet1 -> Eta() - Eta, 2));

	    }
	  const unsigned char DeltaR_index = DeltaR < 1 ? 0 : 1;
	  const float angle = jet1 -> Angle(jet2 -> Vect());
	  const TString hash_key1 = TString("angle_") + 
	    tag_jet_types_[jet1_index] + "_:_" + 
	    tag_jet_types_[jet2_index] + "_" + 
	    tag_levels_[work_mode_] + "_" + 
	    tag_DeltaR_types_[DeltaR_index];
	  plots_ptr_ -> operator[](hash_key1) -> Fill(angle, weight_);
	  const TString hash_key2 = TString("angle_") + 
	    tag_jet_types_[jet1_index] + "_:_" + 
	    tag_jet_types_[jet2_index] + "_" + 
	    tag_levels_[work_mode_] + "_" + 
	    tag_DeltaR_types_[N_DeltaR_types_ - 1];
	  plots_ptr_ -> operator[](hash_key2) -> Fill(angle, weight_);
	  
	}
    }
}

void ColourFlowAnalysisTool::ConfigurePtRadiationProfileTool()
{
  if (PtRadiation_mode_ == 0)
    return;
  gROOT -> SetBatch(kTRUE);
  gStyle -> SetOptStat(0);
  canvas_ = new TCanvas("canvas", "canvas");
  PtRadProf_ = new TH2F("Pt_radiation_profile", "Pt_radiation_profile; X [a.u.]; Y [a.u]", 100, -0.2, 0.2, 50, -0.2, 0.2);
  //  PtRadProf_ -> Draw();
  
}

void ColourFlowAnalysisTool::PtRadiationProfile() const
{
  if (PtRadiation_mode_ == 0)
    return;
  unsigned int size = work_mode_ == 0 ? event_ptr_ -> npf : event_ptr_ -> ngen;
  const TLorentzVector * jets[2] = {leading_light_jet_ptr_, second_leading_light_jet_ptr_};
  const float indices[2] = {leading_light_jet_index_, second_leading_light_jet_index_};
  const TVector2 jet_difference(jets[1] -> Phi() - jets[0] -> Phi(), jets[1] -> Rapidity() - jets[0] -> Rapidity());
  //const float mag = jet_difference.Mod();
  //  const TVector2 jet_difference_norm = jet_difference *2/mag;
  const float tan = jet_difference.Py()/jet_difference.Px();
  const TVector2 normal(-tan/sqrt(1 + tan*tan), 1/sqrt(1 + tan*tan));
  const char direction[2] = {1, -1};
  const float coordinates[2][2] = {{-0.05, 0}, {0.05, 0}};
  for (unsigned char jet_index = 0; jet_index < 2; jet_index ++)
    {
      const unsigned char index = indices[jet_index];
      for (unsigned int jet_const_index = 0; jet_const_index < size; jet_const_index ++)
	{
	  TLorentzVector constituent_4vector;
	  if (work_mode_ == 0)
	    {
	      if (event_ptr_ -> pf_j[jet_const_index] != index)
		continue;
	      /* if (OnlyChargedConstituents and event_ptr_ -> pf_charge[jet_const_index] == 0)
		 continue;*/
	      const float jet_const_energy = sqrt(pow(event_ptr_ -> pf_px[jet_const_index], 2) +
						  pow(event_ptr_ -> pf_py[jet_const_index], 2) + 
						  pow(event_ptr_ -> pf_pz[jet_const_index], 2));
	      constituent_4vector = TLorentzVector    (event_ptr_ -> pf_px[jet_const_index], 
						       event_ptr_ -> pf_py[jet_const_index], 
						       event_ptr_ -> pf_pz[jet_const_index], 
						       jet_const_energy);
	      const float Pt = constituent_4vector.Pt();
	      TVector2 jet_const_phi_eta = TVector2(constituent_4vector.Phi(), constituent_4vector.Rapidity());
	      const TVector2 jet_const_phi_eta_from_jet(jet_const_phi_eta.Px() - jets[jet_index] -> Phi(), 
							jet_const_phi_eta.Py() - jets[jet_index] -> Rapidity());
	      const float parallel_projection = direction[jet_index]*
		(jet_const_phi_eta_from_jet.Px()*jet_difference.Px() + jet_const_phi_eta_from_jet.Py()*jet_difference.Py())/(jet_difference.Mod());
	      const float perpendicular_projection = jet_const_phi_eta_from_jet.Px()*normal.Px() + jet_const_phi_eta_from_jet.Py()*normal.Py();
	      PtRadProf_ -> Fill(coordinates[jet_index][0] + parallel_projection, coordinates[jet_index][1] + perpendicular_projection, Pt*weight_);
	    }
	}
    }
  /*  printf("weight %.8f\n", weight_);
  if (weight_ < 0)
    printf("************************\n");
  */
}

void ColourFlowAnalysisTool::EndPtRadiationProfile() 
{
  if (not canvas_)
    return;
  canvas_ -> cd();
  PtRadProf_ -> Draw("COLZ");
  
  canvas_ -> SaveAs("Radiation_profile.png");
  delete canvas_;

}

void ColourFlowAnalysisTool::Do()
{

}
