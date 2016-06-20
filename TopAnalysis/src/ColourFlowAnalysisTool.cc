#include "TopLJets2015/TopAnalysis/interface/ColourFlowAnalysisTool.hh"
#include "TopLJets2015/TopAnalysis/interface/CFAT_Event.hh"
ColourFlowAnalysisTool::ColourFlowAnalysisTool()
{
  event_display_mode_ = 0;
  canvas_ = 0;
  PtRadProf_ = 0;
  PtRadiation_mode_ = 0;
}

void ColourFlowAnalysisTool::Work()
{
  PlotAngleBetweenJets();
  PlotJetDimensions();
  for (VectorCode_t jet1_code = 0; jet1_code < 2; jet1_code ++)
    {
      for (ChargeCode_t charge_code = 0; charge_code < 2; charge_code ++)
	{
	  try
	    {
	      const PullVector pull_vector = cfat_event_ptr_ -> CalculatePullVector(jet1_code, charge_code);
	      
	      const TString suffix =  TString("_") + tag_charge_types_[charge_code] + "_" + 
		tag_levels_[work_mode_] + "_" + 
		tag_jet_types_[jet1_code] + "_" + 
		tag_channel_;
	      
	      if (fabs(pull_vector.phi_component) > 0.015 or fabs(pull_vector.eta_component) > 0.015)
		{
		  Fill1D(TString("phi_PV_bckg") + suffix, pull_vector.phi_component);
		  Fill1D(TString("eta_PV_bckg") + suffix, pull_vector.eta_component);
		  Fill1D(TString("mag_PV_bckg") + suffix, pull_vector.Mod());
	     	}
	      else
		{
		  Fill1D(TString("phi_PV") + suffix, pull_vector.phi_component);
		  Fill1D(TString("eta_PV") + suffix, pull_vector.eta_component);
		  Fill1D(TString("mag_PV") + suffix, pull_vector.Mod());
		  
		}
	      for (VectorCode_t jet2_code = 0; jet2_code < CFAT_Event::N_jet_types_; jet2_code ++)
		{
		  const TLorentzVector * jet2 = cfat_event_ptr_ -> GetVector(jet2_code);
		  
		  if (jet2_code == jet1_code)
		    continue;
		  if (not jet2)
		    continue;
		  const double DeltaR = cfat_event_ptr_ -> DeltaR(jet1_code, jet2_code);
		  
		  const char * DeltaR_tag = DeltaR <= 1.0 ? tag_DeltaR_types_[DELTAR_LE_1p0] : tag_DeltaR_types_[DELTAR_GT_1p0];
		  
		  try
		    {
		      const double pull_angle = cfat_event_ptr_ -> PullAngle(pull_vector, jet2_code);
		      
		      const double cos_pull_angle = TMath::Cos(pull_angle);
		      const TString infix = TString("_") + tag_charge_types_[charge_code] + "_" + 
			tag_levels_[work_mode_] + "_" + 
			tag_jet_types_[jet1_code] + "_:_" + 
			tag_jet_types_[jet2_code] + "_"; 
		      if (fabs(pull_vector.phi_component) < 0.015 and fabs(pull_vector.eta_component) < 0.015)
			{
			  Fill1D(TString("pull_angle")     + infix + DeltaR_tag       + "_" + tag_channel_, pull_angle);
			  Fill1D(TString("cos_pull_angle") + infix + DeltaR_tag       + "_" + tag_channel_, cos_pull_angle);
			  Fill1D(TString("pull_angle")     + infix + tag_DeltaR_types_[DELTAR_TOTAL] + "_" + tag_channel_, pull_angle);
			  Fill1D(TString("cos_pull_angle") + infix + tag_DeltaR_types_[DELTAR_TOTAL] + "_" + tag_channel_, cos_pull_angle);
			}
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

      for (PF_PTCutCode_t pf_ptcut_code = 1; pf_ptcut_code < 3; pf_ptcut_code ++)
	{
	  try
	    {
	      const PullVector pull_vector = cfat_event_ptr_ -> CalculatePullVector(jet1_code, ALLCOMP, pf_ptcut_code);
	      
	      const TString suffix =  TString("_") + PF_Pt_cuts_types_[pf_ptcut_code] + "_" + 
		tag_levels_[work_mode_] + "_" + 
		tag_jet_types_[jet1_code] + "_" + 
		tag_channel_;
	      
	      if (fabs(pull_vector.phi_component) > 0.015 or fabs(pull_vector.eta_component) > 0.015)
		{
		  Fill1D(TString("phi_PV_bckg") + suffix, pull_vector.phi_component);
		  Fill1D(TString("eta_PV_bckg") + suffix, pull_vector.eta_component);
		  Fill1D(TString("mag_PV_bckg") + suffix, pull_vector.Mod());
	     	}
	      else
		{
		  Fill1D(TString("phi_PV") + suffix, pull_vector.phi_component);
		  Fill1D(TString("eta_PV") + suffix, pull_vector.eta_component);
		  Fill1D(TString("mag_PV") + suffix, pull_vector.Mod());
		  
		}
	      for (VectorCode_t jet2_code = 0; jet2_code < CFAT_Event::N_jet_types_; jet2_code ++)
		{
		  const TLorentzVector * jet2 = cfat_event_ptr_ -> GetVector(jet2_code);
		  
		  if (jet2_code == jet1_code)
		    continue;
		  if (not jet2)
		    continue;
		  const double DeltaR = cfat_event_ptr_ -> DeltaR(jet1_code, jet2_code);
		  const char * DeltaR_tag = DeltaR <= 1.0 ? tag_DeltaR_types_[DELTAR_LE_1p0] : tag_DeltaR_types_[DELTAR_GT_1p0];
		  
		  try
		    {
		      const double pull_angle = cfat_event_ptr_ -> PullAngle(pull_vector, jet2_code);
		      
		      const double cos_pull_angle = TMath::Cos(pull_angle);
		      const TString infix = TString("_") + PF_Pt_cuts_types_[pf_ptcut_code] + "_" + 
			tag_levels_[work_mode_] + "_" + 
			tag_jet_types_[jet1_code] + "_:_" + 
			tag_jet_types_[jet2_code] + "_"; 
		      if (fabs(pull_vector.phi_component) < 0.015 and fabs(pull_vector.eta_component) < 0.015)
			{
			  Fill1D(TString("pull_angle")     + infix + DeltaR_tag       + "_" + tag_channel_, pull_angle);
			  Fill1D(TString("cos_pull_angle") + infix + DeltaR_tag       + "_" + tag_channel_, cos_pull_angle);
			  Fill1D(TString("pull_angle")     + infix + tag_DeltaR_types_[DELTAR_TOTAL] + "_" + tag_channel_, pull_angle);
			  Fill1D(TString("cos_pull_angle") + infix + tag_DeltaR_types_[DELTAR_TOTAL] + "_" + tag_channel_, cos_pull_angle);
			}
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


	  try
	    {
	      const PullVector pull_vector = cfat_event_ptr_ -> CalculatePullVector(jet1_code);
	      const char * PF_N_cuts_tag = pull_vector.Ncomponents <= 20 ? PF_N_cuts_types_[PFN_LE_20] : PF_N_cuts_types_[PFN_GT_20];
	      const char * HadW_Pt_cuts_tag = "TBD";
	      const char * PVMag_cuts_tag = pull_vector.Mod() <= 0.005 ? PVMag_cuts_types_[PVMAG_LE_0p005] :PVMag_cuts_types_[PVMAG_GT_0p005];
		
      //const char * tags[3] = {PF_N_cuts_tag, HadW_Pt_cuts_tag, PVMag_cuts_tag};
	      {
		const TString suffix =  TString("_") + PF_N_cuts_tag + "_" + 
		  tag_levels_[work_mode_] + "_" + 
		  tag_jet_types_[jet1_code] + "_" + 
		  tag_channel_;
	      
		if (fabs(pull_vector.phi_component) > 0.015 or fabs(pull_vector.eta_component) > 0.015)
		  {
		    Fill1D(TString("phi_PV_bckg") + suffix, pull_vector.phi_component);
		    Fill1D(TString("eta_PV_bckg") + suffix, pull_vector.eta_component);
		    Fill1D(TString("mag_PV_bckg") + suffix, pull_vector.Mod());
		  }
		else
		  {
		    Fill1D(TString("phi_PV") + suffix, pull_vector.phi_component);
		    Fill1D(TString("eta_PV") + suffix, pull_vector.eta_component);
		    Fill1D(TString("mag_PV") + suffix, pull_vector.Mod());
		  
		  }
	      }
	      if (cfat_event_ptr_ -> GetVector(HAD_W))
		{
		  HadW_Pt_cuts_tag = cfat_event_ptr_ -> GetVector(HAD_W) -> Pt() <= 50.0 ? HadW_Pt_cuts_types_[HADW_PT_LE_50p0_GEV] : HadW_Pt_cuts_types_[HADW_PT_GT_50p0_GEV];

	      
		  const TString suffix =  TString("_") + HadW_Pt_cuts_tag + "_" + 
		    tag_levels_[work_mode_] + "_" + 
		    tag_jet_types_[jet1_code] + "_" + 
		    tag_channel_;
	      
		  if (fabs(pull_vector.phi_component) > 0.015 or fabs(pull_vector.eta_component) > 0.015)
		    {
		      Fill1D(TString("phi_PV_bckg") + suffix, pull_vector.phi_component);
		      Fill1D(TString("eta_PV_bckg") + suffix, pull_vector.eta_component);
		      Fill1D(TString("mag_PV_bckg") + suffix, pull_vector.Mod());
		    }
		  else
		    {
		      Fill1D(TString("phi_PV") + suffix, pull_vector.phi_component);
		      Fill1D(TString("eta_PV") + suffix, pull_vector.eta_component);
		      Fill1D(TString("mag_PV") + suffix, pull_vector.Mod());
		  
		    }
		}
	      for (VectorCode_t jet2_code = 0; jet2_code < CFAT_Event::N_jet_types_; jet2_code ++)
		{
		  const TLorentzVector * jet2 = cfat_event_ptr_ -> GetVector(jet2_code);
		  
		  if (jet2_code == jet1_code)
		    continue;
		  if (not jet2)
		    continue;
		  const double DeltaR = cfat_event_ptr_ -> DeltaR(jet1_code, jet2_code);
		  const char * DeltaR_tag = DeltaR <= 1.0 ? tag_DeltaR_types_[DELTAR_LE_1p0] : tag_DeltaR_types_[DELTAR_GT_1p0];
		  
		  try
		    {
		      const double pull_angle = cfat_event_ptr_ -> PullAngle(pull_vector, jet2_code);
		      
		      const double cos_pull_angle = TMath::Cos(pull_angle);
		      {
			const TString infix = TString("_") + PF_N_cuts_tag + "_" + 
			  tag_levels_[work_mode_] + "_" + 
			  tag_jet_types_[jet1_code] + "_:_" + 
			  tag_jet_types_[jet2_code] + "_"; 
			if (fabs(pull_vector.phi_component) < 0.015 and fabs(pull_vector.eta_component) < 0.015)
			  {
			    Fill1D(TString("pull_angle")     + infix + DeltaR_tag       + "_" + tag_channel_, pull_angle);
			    Fill1D(TString("cos_pull_angle") + infix + DeltaR_tag       + "_" + tag_channel_, cos_pull_angle);
			    Fill1D(TString("pull_angle")     + infix + tag_DeltaR_types_[DELTAR_TOTAL] + "_" + tag_channel_, pull_angle);
			    Fill1D(TString("cos_pull_angle") + infix + tag_DeltaR_types_[DELTAR_TOTAL] + "_" + tag_channel_, cos_pull_angle);
			  }
		      }
		
		      if (cfat_event_ptr_ -> GetVector(HAD_W))
			{
			  const TString infix = TString("_") + HadW_Pt_cuts_tag + "_" + 
			    tag_levels_[work_mode_] + "_" + 
			    tag_jet_types_[jet1_code] + "_:_" + 
			    tag_jet_types_[jet2_code] + "_"; 
			  if (fabs(pull_vector.phi_component) < 0.015 and fabs(pull_vector.eta_component) < 0.015)
			    {
			      Fill1D(TString("pull_angle")     + infix + DeltaR_tag       + "_" + tag_channel_, pull_angle);
			      Fill1D(TString("cos_pull_angle") + infix + DeltaR_tag       + "_" + tag_channel_, cos_pull_angle);
			      Fill1D(TString("pull_angle")     + infix + tag_DeltaR_types_[DELTAR_TOTAL] + "_" + tag_channel_, pull_angle);
			      Fill1D(TString("cos_pull_angle") + infix + tag_DeltaR_types_[DELTAR_TOTAL] + "_" + tag_channel_, cos_pull_angle);
			    }
			}
		

		      {
			const TString infix = TString("_") + PVMag_cuts_tag + "_" + 
			  tag_levels_[work_mode_] + "_" + 
			  tag_jet_types_[jet1_code] + "_:_" + 
			  tag_jet_types_[jet2_code] + "_"; 
			if (fabs(pull_vector.phi_component) < 0.015 and fabs(pull_vector.eta_component) < 0.015)
			  {
			    Fill1D(TString("pull_angle")     + infix + DeltaR_tag       + "_" + tag_channel_, pull_angle);
			    Fill1D(TString("cos_pull_angle") + infix + DeltaR_tag       + "_" + tag_channel_, cos_pull_angle);
			    Fill1D(TString("pull_angle")     + infix + tag_DeltaR_types_[DELTAR_TOTAL] + "_" + tag_channel_, pull_angle);
			    Fill1D(TString("cos_pull_angle") + infix + tag_DeltaR_types_[DELTAR_TOTAL] + "_" + tag_channel_, cos_pull_angle);
			  }
		      }
		
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
  // AnalyseParticleFlow();
  //PtRadiationProfile();
  
}

void ColourFlowAnalysisTool::Fill1D(const TString & key, double value) const
{
  plots_ptr_ -> operator[](key) -> Fill(value, cfat_event_ptr_ -> weight_);

} 

void ColourFlowAnalysisTool::SetWorkMode(WorkCode_t mode)
{
  work_mode_ = mode;
  cfat_event_ptr_ -> SetWorkMode(mode);
}

/*
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
  const TVector2 v1(fakelv.Phi(), fakelv.Rapidity());
  if (v1.Mod() > 1E-12)
    vect_jets.push_back(&fakelv);
  else
    vect_jets.push_back(NULL);
  
  return vect_jets;

}


*/
void ColourFlowAnalysisTool::PlotAngleBetweenJets() const
{
  for (VectorCode_t jet1_code = 0; jet1_code < CFAT_Event::N_jet_types_; jet1_code ++)
    {
      const TLorentzVector * jet1 = cfat_event_ptr_ -> GetVector(jet1_code);
      if (not jet1)
	continue;
      for (VectorCode_t jet2_code = jet1_code + 1; jet2_code < CFAT_Event::N_jet_types_; jet2_code ++)
	{
	  const TLorentzVector * jet2 = cfat_event_ptr_ -> GetVector(jet2_code);
	  if (not jet2)
	    continue;
	  const double DeltaR =  cfat_event_ptr_ -> DeltaR(jet1_code, jet2_code);
	  const unsigned char DeltaR_index = DeltaR < 1.0 ? DELTAR_LE_1p0 : DELTAR_GT_1p0;
          
	  const TString infix = TString(tag_jet_types_[jet1_code]) + "_:_" + 
	    tag_jet_types_[jet2_code] + "_" + 
	    tag_levels_[work_mode_] + "_"; 
	    
	  const TString hash_key1_angle = TString("angle_") + 
	    infix +
	    tag_DeltaR_types_[DeltaR_index];


	  const double angle = cfat_event_ptr_ -> Angle(jet1_code, jet2_code);
	  Fill1D(hash_key1_angle, angle);
	  const TString hash_key2_angle = TString("angle_") + 
	    infix +
	    tag_DeltaR_types_[DELTAR_TOTAL];
	  Fill1D(hash_key2_angle, angle);
	  const TVector2 jet_difference_phi_eta(TVector2::Phi_mpi_pi(jet2 -> Phi() - jet1 -> Phi()), jet2 -> Rapidity() - jet1 -> Rapidity());
	  {
	    const TString hash_key = TString("jet_dif_phi_") + 
	      infix +
	      tag_DeltaR_types_[DeltaR_index];
	    Fill1D(hash_key, jet_difference_phi_eta.Px());
	  					     
	  }
	
	  {
	    const TString hash_key = TString("jet_dif_phi_") + 
	      infix +
	      tag_DeltaR_types_[DELTAR_TOTAL];
	    Fill1D(hash_key, jet_difference_phi_eta.Px());
	  					     
	  }
	
	  {
	    const TString hash_key = TString("jet_dif_eta_") + 
	      infix +
	      tag_DeltaR_types_[DeltaR_index];
	    Fill1D(hash_key, jet_difference_phi_eta.Py());
	  					     
	  }
	
	  {
	    const TString hash_key = TString("jet_dif_eta_") + 
	      infix +
	      tag_DeltaR_types_[DELTAR_TOTAL];
	    Fill1D(hash_key, jet_difference_phi_eta.Py());
	  					     
	  }


       	}
    }
}

void ColourFlowAnalysisTool::PlotJetDimensions() const
{
  for (VectorCode_t jet_code = 0; jet_code < CFAT_Event::N_jet_types_; jet_code ++)
    {
      const TLorentzVector *jet = cfat_event_ptr_ -> GetVector(jet_code); 
      
      if (not jet)
	continue;
      const TString sufix = TString(tag_levels_[work_mode_]) + "_" +
	tag_jet_types_[jet_code];
      //printf("test %s %s \n", (TString("jet_phi_") +  sufix).Data(), (TString("jet_rapidity_") +  sufix).Data()); 
      if (jet == CFAT_Event::beam_ptr_)
	continue;
      Fill1D(TString("jet_phi_")       + sufix,  jet -> Phi());
      Fill1D(TString("jet_rapidity_")  + sufix,  jet -> Rapidity());
      Fill1D(TString("jet_eta_")       + sufix,  jet -> Eta());
      const double P = sqrt(pow(jet -> Pt(), 2) + pow(jet -> Pz(), 2));
      Fill1D(TString("jet_mass_norm_") + sufix,  jet -> M()/P);
      Fill1D(TString("jet_pt_norm_")   + sufix,  jet -> Pt()/P);
      Fill1D(TString("jet_pz_norm_")   + sufix,  jet -> Pz()/P);
      Fill1D(TString("jet_px_norm_")   + sufix,  jet -> Px()/P);

      Fill1D(TString("jet_mass_")      + sufix,  jet -> M());
      Fill1D(TString("jet_pt_")        + sufix,  jet -> Pt());
      Fill1D(TString("jet_pz_")        + sufix,  jet -> Pz());
      Fill1D(TString("jet_px_")        + sufix,  jet -> Px());

    }
}
/*
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
*/
void ColourFlowAnalysisTool::Do()
{

}
