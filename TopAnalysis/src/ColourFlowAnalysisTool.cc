#include "TopLJets2015/TopAnalysis/interface/ColourFlowAnalysisTool.hh"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
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
  PlotJetPhiAndEta();
  for (unsigned char jet1_index = 0; jet1_index < 2; jet1_index ++)
    {
      const unsigned int jet_index = work_mode_ == 0 ? jet_indices_[jet1_index] : event_ptr_ -> j_g[jet_indices_[jet1_index]];

      const TLorentzVector charged_jet = GetChargedJet(jet_index);
      const TLorentzVector * jet1_array[2] = {vect_jets_[jet1_index], &charged_jet};
      for (unsigned char charge_index = 0; charge_index < 2; charge_index ++)
	{
	  const TLorentzVector * jet1 = jet1_array[charge_index];
	  try
	    {
	      PlotPUPPIWeight(jet1_index, jet_index, charge_index);
	      
	      const PullVector pull_vector = CalculatePullVector(*jet1, jet_index, OnlyChargedConstituents[charge_index]);
	      
	      const TString suffix =  TString("_") + tag_charge_types_[charge_index] + "_" + 
		tag_levels_[work_mode_] + "_" + 
		tag_jet_types_[jet1_index] + "_" + 
		tag_channel_;
	      plots_ptr_ -> operator[](TString("phi_PV") + suffix) -> Fill(pull_vector.phi_component, weight_);
	      plots_ptr_ -> operator[](TString("eta_PV") + suffix) -> Fill(pull_vector.eta_component, weight_);
	      plots_ptr_ -> operator[](TString("mag_PV") + suffix) -> Fill(pull_vector.Mod(), weight_);
	      if (pull_vector.Mod() > 0.015)
		{
		  plots_ptr_ -> operator[](TString("phi_PV_bckg") + suffix) -> Fill(pull_vector.phi_component, weight_);
		  plots_ptr_ -> operator[](TString("eta_PV_bckg") + suffix) -> Fill(pull_vector.eta_component, weight_);
		  plots_ptr_ -> operator[](TString("mag_PV_bckg") + suffix) -> Fill(pull_vector.Mod(), weight_);
	     

		}
	      const TVector2 pull_vector_Euclidian = CalculatePullVectorEuclidian(*jet1, jet_index, OnlyChargedConstituents[charge_index]);
	      {
		const TString phi_PV_Euclidian_hash_key = TString("phi_PV_Euclidian_") +
		  tag_charge_types_[charge_index] + "_" + 
		  tag_levels_[work_mode_] + "_" +
		  tag_jet_types_[jet1_index];
		plots_ptr_ -> operator[](phi_PV_Euclidian_hash_key) -> Fill(pull_vector_Euclidian.Phi(), weight_);
	      }
	      

	      for (unsigned char jet2_index = 0; jet2_index < N_jet_types_; jet2_index ++)
		{
		  const TLorentzVector * jet2 = vect_jets_[jet2_index];
		  
		  if (jet2_index == jet1_index)
		    continue;
		  if (not jet2)
		    continue;
		  //		  printf("jet1 %u jet2 %u\n", jet1_index, jet2_index);
		  double DeltaR = 0;
		  if (jet2 != & beam_)
		    DeltaR = jet1 -> DeltaR(*jet2);
		  else
		    {
		      const double Theta = jet2 -> Theta();
		      const double Eta = -0.5 * TMath::Log(TMath::Tan(Theta/2));
		      DeltaR = (pow(jet1 -> Phi() - jet2 -> Phi(), 2) + pow(jet1 -> Eta() - Eta, 2));

		    }
		  const char * DeltaR_tag = DeltaR <= 1.0 ? tag_DeltaR_types_[0] : tag_DeltaR_types_[1];
		  
		  try
		    {
		      /*
		      if (jet1_index == 0 and jet2_index == 6 )
			{
			  event_display_mode_ = 1;
			  PullAngle(pull_vector, jet_difference);
			}*/
		      double pull_angle = 0;
		      if (jet2 == & beam_)
			pull_angle = TMath::ACos(pull_vector.eta_component/pull_vector.Mod());
		      else
			{
			  const TVector2 jet_difference(/*TVector2::Phi_mpi_pi*/(jet2 -> Phi() - jet1 -> Phi()), jet2 -> Rapidity() - jet1 -> Rapidity());

			  pull_angle = PullAngle(pull_vector, jet_difference);
			}
		      
		      static unsigned short count = 0;

		      if (count < 60 )
			if (jet1_index == 0 and jet2_index == 1)
			  {
			    event_display_mode_ = 1;
			    //			    EventDisplay(pull_vector, pull_angle, OnlyChargedConstituents[charge_index]);
			    count ++;
			    //  getchar();
			  }
		      
		      const double cos_pull_angle = TMath::Cos(pull_angle);
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
  // AnalyseParticleFlow();
  //PtRadiationProfile();
  
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
  double mass_dif[2][2] = {{1000, 1000}, {1000, 1000}};
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
  vect_jets.push_back(lepton_ptr_);
  vect_jets.push_back(neutrino_ptr_);

  const double pt = sqrt(pow(second_leading_jet -> Px(), 2) + pow(second_leading_jet -> Py(), 2));
  const double x = pt * TMath::Cos(second_leading_jet -> Rapidity());
  const double y = sqrt(pow(pt, 2) - pow(x, 2));
  const double A = TMath::Exp(- 2 * second_leading_jet -> Phi());
  const double E = second_leading_jet -> E();
  const double z = E * (A - 1) / (A + 1);
  static TLorentzVector fakelv;
  fakelv.SetPxPyPzE(x, y, z, E);
  
  if (second_leading_jet -> Rapidity() < 0)
    fakelv .SetPz( -fakelv.Pz());
  //vect_jets.push_back(&fakelv);
  const TVector2 v1(fakelv.Phi(), fakelv.Rapidity());
  if (v1.Mod() > 1E-12)
    vect_jets.push_back(&fakelv);
  else
    vect_jets.push_back(NULL);
  //const TVector2 v2(second_leading_jet -> Phi(), second_leading_jet -> Rapidity());
  //const double cos_angle = (v1.Px()*v2.Px() + v1.Py()*v2.Py())/(v1.Mod()*v2.Mod());
  //printf("cos_angle %f\n", cos_angle); getchar();

  return vect_jets;

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
	  double DeltaR = 0;
	  if (jet1 != & beam_ and jet2 != & beam_)
	    DeltaR = jet1 -> DeltaR(*jet2);
	  else
	    {
	      DeltaR = 100.0;
	      /*
	      if (jet2 == &beam_)
		{
		  const float Theta1 = jet1 -> Theta();
		  const float Eta1 = -0.5 * TMath::Log(TMath::Tan(Theta1/2));
		  DeltaR = (pow(jet1 -> Phi() - jet2 -> Phi(), 2) + pow(Eta1 - jet2 -> Eta(), 2));
		}
	    
	      if (jet1 == &beam_)
		{
		  const float Theta2 = jet2 -> Theta();
		  const float Eta2 = -0.5 * TMath::Log(TMath::Tan(Theta2/2));
		  DeltaR = (pow(jet1 -> Phi() - jet2 -> Phi(), 2) + pow(jet1 -> Eta() - Eta2, 2));
		  }*/
	    }
	  const unsigned char DeltaR_index = DeltaR < 1 ? 0 : 1;
          //const TLorentzVector jet_difference = *jet2 - *jet1;
	  
	  const TString infix = TString(tag_jet_types_[jet1_index]) + "_:_" + 
	    tag_jet_types_[jet2_index] + "_" + 
	    tag_levels_[work_mode_] + "_"; 
	    
	  const TString hash_key1_angle = TString("angle_") + 
	    infix +
	    tag_DeltaR_types_[DeltaR_index];


	  double angle = 0;
	  
	  if (jet1 == & beam_)
	    angle = TMath::ACos(jet2 -> Pz() / jet2 -> Mag());
	  else if (jet2 == & beam_)
	    angle = TMath::ACos(jet1 -> Pz() / jet1 -> Mag());
	  else
	    angle = jet1 -> Angle(jet2 -> Vect());
	  plots_ptr_ -> operator[](hash_key1_angle) -> Fill(angle, weight_);
	  const TString hash_key2_angle = TString("angle_") + 
	    infix +
	    tag_DeltaR_types_[N_DeltaR_types_ - 1];
	  plots_ptr_ -> operator[](hash_key2_angle) -> Fill(angle, weight_);
	  const TVector2 jet_difference_phi_eta(TVector2::Phi_mpi_pi(jet2 -> Phi() - jet1 -> Phi()), jet2 -> Rapidity() - jet1 -> Rapidity());
	  {
	    const TString hash_key = TString("jet_dif_phi_") + 
	      infix +
	      tag_DeltaR_types_[DeltaR_index];
	    plots_ptr_ -> operator[](hash_key) -> Fill(jet_difference_phi_eta.Px(), weight_);
	  					     
	  }
	
	  {
	    const TString hash_key = TString("jet_dif_phi_") + 
	      infix +
	      tag_DeltaR_types_[N_DeltaR_types_ - 1];
	    plots_ptr_ -> operator[](hash_key) -> Fill(jet_difference_phi_eta.Px(), weight_);
	  					     
	  }
	
	  {
	    const TString hash_key = TString("jet_dif_eta_") + 
	      infix +
	      tag_DeltaR_types_[DeltaR_index];
	    plots_ptr_ -> operator[](hash_key) -> Fill(jet_difference_phi_eta.Py(), weight_);
	  					     
	  }
	
	  {
	    const TString hash_key = TString("jet_dif_eta_") + 
	      infix +
	      tag_DeltaR_types_[N_DeltaR_types_ - 1];
	    plots_ptr_ -> operator[](hash_key) -> Fill(jet_difference_phi_eta.Py(), weight_);
	  					     
	  }


       	}
    }
}

void ColourFlowAnalysisTool::PlotJetPhiAndEta() const
{
  for (unsigned char jet_index = 0; jet_index < vect_jets_.size(); jet_index ++)
    {
      if (not vect_jets_[jet_index])
	continue;

      //printf("index %u\n", jet_index);
      const TString sufix = TString(tag_levels_[work_mode_]) + "_" +
	tag_jet_types_[jet_index];
      //printf("test %s %s \n", (TString("jet_phi_") +  sufix).Data(), (TString("jet_rapidity_") +  sufix).Data()); 
      const TLorentzVector *jet = vect_jets_[jet_index];
      if (jet == &beam_)
	continue;
      plots_ptr_ -> operator[](TString("jet_phi_") +  sufix) -> Fill(jet -> Phi(), weight_);
      plots_ptr_ -> operator[](TString("jet_rapidity_") +  sufix) -> Fill(jet -> Rapidity(), weight_);
      plots_ptr_ -> operator[](TString("jet_eta_") +  sufix) -> Fill(jet -> Eta(), weight_);
      const double P = sqrt(pow(jet -> Pt(), 2) + pow(jet -> Pz(), 2));
      plots_ptr_ -> operator[](TString("jet_mass_norm_") +  sufix) -> Fill(jet -> M()/P, weight_);
      plots_ptr_ -> operator[](TString("jet_pt_norm_") +  sufix) -> Fill(jet -> Pt()/P, weight_);
      plots_ptr_ -> operator[](TString("jet_pz_norm_") +  sufix) -> Fill(jet -> Pz()/P, weight_);
      plots_ptr_ -> operator[](TString("jet_px_norm_") +  sufix) -> Fill(jet -> Px()/P, weight_);

      plots_ptr_ -> operator[](TString("jet_mass_") +  sufix) -> Fill(jet -> M(), weight_);
      plots_ptr_ -> operator[](TString("jet_pt_") +  sufix) -> Fill(jet -> Pt(), weight_);
      plots_ptr_ -> operator[](TString("jet_pz_") +  sufix) -> Fill(jet -> Pz(), weight_);
      plots_ptr_ -> operator[](TString("jet_px_") +  sufix) -> Fill(jet -> Px(), weight_);

    }
}

void ColourFlowAnalysisTool::PlotPUPPIWeight(unsigned char jet_vector_index, unsigned char jet_index, unsigned char charge_index) const
{
  const TString hash_key = TString("PUPPI_weight_") +
		tag_charge_types_[charge_index] + "_" + 
		tag_levels_[work_mode_] + "_" +
		tag_jet_types_[jet_vector_index]; 
  //printf("test %s\n", hash_key.Data());
  static const bool OnlyChargedConstituents[2] = {false, true};
  unsigned int size = work_mode_ == 0 ? event_ptr_ -> npf : event_ptr_ -> ngpf;
  for (unsigned int jet_const_index = 0; jet_const_index < size; jet_const_index ++)
    {
      if (work_mode_ == 0)
	{
	  if (event_ptr_ -> pf_j[jet_const_index] != jet_index)
	    continue;
	  if (OnlyChargedConstituents[charge_index] and event_ptr_ -> pf_id[jet_const_index] == 0)
	    continue;
	  plots_ptr_ -> operator[](hash_key) -> Fill(event_ptr_ -> pf_puppiWgt[jet_const_index], weight_);

	}
      else
	{
	  return;

	}
		
    }
  
}

void ColourFlowAnalysisTool::Do()
{

}
