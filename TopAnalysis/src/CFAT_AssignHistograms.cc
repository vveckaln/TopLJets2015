#include "TopLJets2015/TopAnalysis/interface/ColourFlowAnalysisTool.hh"

const char* ColourFlowAnalysisTool::tag_channel_                    = "4j2t";
const char * ColourFlowAnalysisTool::tag_charge_types_[2]           = {"allconst", "chconst"};
const unsigned char ColourFlowAnalysisTool::N_jet_types_            = 10;
const char* ColourFlowAnalysisTool::tag_jet_types_[N_jet_types_]    = {"leading_jet", "2nd_leading_jet", "had_b", "lept_b", "had_t", "lept_t", "beam", "lepton", "neutrino", "fake"};
const unsigned char ColourFlowAnalysisTool::N_DeltaR_types_         = 3;
const char * ColourFlowAnalysisTool::tag_DeltaR_types_[N_DeltaR_types_] = {"DeltaRle1.0", "DeltaRgt1.0", "DeltaRTotal"}; 
const char * ColourFlowAnalysisTool::tag_levels_[2]                 = {"reco", "gen"};
const TLorentzVector ColourFlowAnalysisTool::beam_                   = TLorentzVector(1E-2, 0, 1E10, sqrt(1E-4 + 1E20)+1E-6);


void ColourFlowAnalysisTool::AssignHistograms() 
{
  static const unsigned char ncategories_PA         = 2;
  static const char * cat_title_PA[ncategories_PA]     = {"pull_angle",     "cos_pull_angle" };
  static const unsigned short nbins_PA[ncategories_PA] = {100,               100              };
  static const double min_PA[ncategories_PA]           = {-1.2*TMath::Pi(), -1.2             };
  static const double max_PA[ncategories_PA]           = {- min_PA[0],         -min_PA[1]          };
  static const char* axes_PA[ncategories_PA]           = {"; pull angle [rad]; Events",  "; cos; Events"  };
  for (unsigned char cat_index = 0; cat_index < ncategories_PA; cat_index++)
    {
      for (unsigned char level_index = 0; level_index < 2; level_index ++)
	{
	  for (unsigned char charge_index = 0; charge_index < 2; charge_index ++)
	    {
	      for (unsigned char jet1_index = 0; jet1_index < 2; jet1_index ++)
		{
		  for (unsigned char jet2_index = 0; jet2_index < N_jet_types_; jet2_index ++)
		    {
		      for (unsigned char DeltaR_index = 0; DeltaR_index < N_DeltaR_types_; DeltaR_index ++)
			{
	
			  if (jet1_index == jet2_index)
			    continue;
			  const TString hash_key = TString(cat_title_PA[cat_index]) + "_" + 
			    tag_charge_types_[charge_index] + "_" + 
			    tag_levels_[level_index] + "_" +
			    tag_jet_types_[jet1_index] + "_:_" + 
			    tag_jet_types_[jet2_index] + "_" + 
			    tag_DeltaR_types_[DeltaR_index] + "_" + 
			    tag_channel_;
			  
			  plots_ptr_ -> operator[](hash_key) = new TH1F(hash_key, 
									hash_key + axes_PA[cat_index], 
									nbins_PA[cat_index], 
									min_PA[cat_index], 
									max_PA[cat_index]);
			}
		    }
		}
    
	    }
	}
    }

  static const unsigned char ncategories_PV = 6;
  static const char * cat_title_PV[ncategories_PV]     = {"phi_PV",                "eta_PV",               "mag_PV",                 "phi_PV_bckg",           "eta_PV_bckg",            "mag_PV_bckg"};
  static const unsigned short nbins_PV[ncategories_PV] = {100,                     100,                     100,                      100,                     100,                      100          };
  static const double min_PV[ncategories_PV]           = {-0.01*TMath::Pi(),       -0.05,                   0.0,                      - TMath::Pi() -0.2,      - 5.0,                    0.0          };
  static const double max_PV[ncategories_PV]           = {- min_PV[0],             -min_PV[1],              0.05,                     TMath::Pi() + 0.2,       5.0,                      6.0          };
  static const char* axes_PV[ncategories_PV]           = {"; #phi [rad]; Events",  "; #eta [a.u.]; Events", "; magn. [a.u.]; Events", "; #phi [rad]; Events",   "; #eta [a.u.]; Events", "; magn. [a.u.]; Events"};
  for (unsigned char cat_index = 0; cat_index < ncategories_PV; cat_index++)
    {
      for (unsigned char level_index = 0; level_index < 2; level_index ++)
	{ 
	  for (unsigned char charge_index = 0; charge_index < 2; charge_index ++)
	    {
	      for (unsigned char jet1_index = 0; jet1_index < 2; jet1_index ++)
		{
		  const TString hash_key = TString(cat_title_PV[cat_index]) + "_" + 
		    tag_charge_types_[charge_index] + "_" + 
		    tag_levels_[level_index] + "_" +
		    tag_jet_types_[jet1_index] + "_" + 
		    tag_channel_;
		  plots_ptr_ -> operator[](hash_key) = new TH1F(hash_key, 
								hash_key + axes_PV[cat_index], 
								nbins_PV[cat_index], 
								min_PV[cat_index], 
								max_PV[cat_index]);
    
		
		}
	    }
	}
    }

  
  static const unsigned char ncategories_CTRL              = 2;
  static const char * branch[2]                            = {"had", "lept"};
  static const char * cat_title_CTRL[ncategories_CTRL]     = {"W_mass",         "t_mass"    };
  static const unsigned short nbins_CTRL[ncategories_CTRL] = {100,              100         };
  static const double min_CTRL[ncategories_CTRL]           = {0,                50          };
  static const double max_CTRL[ncategories_CTRL]           = {150,              400         };
  static const char* axes_CTRL[ncategories_CTRL]           = {"; GeV; Events",  "; GeV; Events"};
  for (unsigned char branch_index = 0; branch_index < 2; branch_index ++)
    {
      for (unsigned char level_index = 0; level_index < 2; level_index ++)
	{
	  for (unsigned char cat_index = 0; cat_index < ncategories_CTRL; cat_index ++)
	    { 
	      const TString hash_key = TString(branch[branch_index]) + "_" + tag_levels_[level_index] + "_" + cat_title_CTRL[cat_index]; 
	      plots_ptr_ -> operator[](hash_key) = new TH1F(hash_key, hash_key + axes_CTRL[cat_index], nbins_CTRL[cat_index], min_CTRL[cat_index], max_CTRL[cat_index]);
	    }
	}
    }

  const char * particles[2] = {"all", "jet"};
  
  for (unsigned char DeltaR_index = 0; DeltaR_index < N_DeltaR_types_; DeltaR_index ++)
    {
      for (unsigned char level_index = 0; level_index < 2; level_index ++)
	{ 
	  for (unsigned char jet1_index = 0; jet1_index < N_jet_types_; jet1_index ++)
	    {
	      for (unsigned char jet2_index = jet1_index + 1; jet2_index < N_jet_types_; jet2_index ++)
		{
		  {
		    const TString hash_key = TString("angle_") + tag_jet_types_[jet1_index] + "_:_" + tag_jet_types_[jet2_index] + "_" + tag_levels_[level_index] + "_" + tag_DeltaR_types_[DeltaR_index];
		    plots_ptr_ -> operator[](hash_key) = new TH1F(hash_key, hash_key + "; angle [rad]; Events", 50, -0.1, TMath::Pi() + 0.1);
		  }
		  {
		    const TString hash_key = TString("jet_dif_phi_") + tag_jet_types_[jet1_index] + "_:_" + tag_jet_types_[jet2_index] + "_" + tag_levels_[level_index] + "_" + tag_DeltaR_types_[DeltaR_index];
		    plots_ptr_ -> operator[](hash_key) = new TH1F(hash_key, hash_key + "; #phi [rad]; Events", 50, -1.2 * TMath::Pi(), 1.2 * TMath::Pi());
		  
		  }
		  {
		    const TString hash_key = TString("jet_dif_eta_") + tag_jet_types_[jet1_index] + "_:_" + tag_jet_types_[jet2_index] + "_" + tag_levels_[level_index] + "_" + tag_DeltaR_types_[DeltaR_index];
		    plots_ptr_ -> operator[](hash_key) = new TH1F(hash_key, hash_key + "; #eta; Events", 50, -5.2, 5.2);
		  
		  }
		}
	    }
	  for (unsigned char charge_index = 0; charge_index < 2; charge_index ++)
	    {
	      for (unsigned char particle_index = 0; particle_index < 2; particle_index ++)
		{
		  const TString hash_key = TString("chi_") + tag_charge_types_[charge_index] + "_" + tag_levels_[level_index] + "_" + tag_DeltaR_types_[DeltaR_index] + "_" + particles[particle_index];
		  plots_ptr_ -> operator[](hash_key) = new TH1F(hash_key, hash_key + "; #chi; Events", 50, -0.1, 2.1);
		}
	    }
	}
    }
  plots_ptr_ -> operator[]("had_b_flavour")              = new TH1F("had_b_flavour", "had_b_flavour; flavour; Events", 21, -10.5, 10.5);
  plots_ptr_ -> operator[]("lept_b_flavour")             = new TH1F("lept_b_flavour", "lept_b_flavour; flavour; Events", 21, -10.5, 10.5);

  plots_ptr_ -> operator[]("leading_jet_flavour")        = new TH1F("leading_jet_flavour", "leading_jet_flavour; flavour; Events", 21, -10.5, 10.5);
  plots_ptr_ -> operator[]("second_leading_jet_flavour") = new TH1F("second_leading_jet_flavour", "second_leading_jet_flavour; flavour; Events", 21, -10.5, 10.5);
  for (unsigned char level_index = 0; level_index < 2; level_index ++)
    { 
      for (unsigned char charge_index = 0; charge_index < 2; charge_index ++)
	{
	  for (unsigned char jet1_index = 0; jet1_index < 2; jet1_index ++)
	    {

	      {
		const TString hash_key = TString("phi_PV_Euclidian_") +
		  tag_charge_types_[charge_index] + "_" + 
		  tag_levels_[level_index] + "_" +
		  tag_jet_types_[jet1_index]; 
		plots_ptr_ -> operator[](hash_key)       = new TH1F(hash_key, hash_key + "; #phi[rad]; Events", 50, -0.2, 2 * TMath::Pi() + 0.2);
	      }
	      
	    }
	}
    }
  
  for (unsigned char level_index = 0; level_index < 1; level_index ++)
    { 
      for (unsigned char charge_index = 0; charge_index < 2; charge_index ++)
	{
	  for (unsigned char jet1_index = 0; jet1_index < 2; jet1_index ++)
	    {
	      const TString hash_key = TString("PUPPI_weight_") +
		tag_charge_types_[charge_index] + "_" + 
		tag_levels_[level_index] + "_" +
		tag_jet_types_[jet1_index]; 
	      //printf("%s\n", hash_key.Data());
	      plots_ptr_ -> operator[](hash_key)       = new TH1F(hash_key, hash_key + "; weight; Events", 25, -0.1, 1.1);
	    }
	}
    }

  for (unsigned char level_index = 0; level_index < 2; level_index ++)
    { 
      for (unsigned char jet1_index = 0; jet1_index < N_jet_types_; jet1_index ++)
	{
	  const TString hash_key_phi = TString("jet_phi_") +
		tag_levels_[level_index] + "_" +
		tag_jet_types_[jet1_index]; 
	  plots_ptr_ -> operator[](hash_key_phi)       = new TH1F(hash_key_phi, hash_key_phi + "; #phi [rad]; Events", 50, -0.2 - TMath::Pi(), TMath::Pi() + 0.2);
	  const TString hash_key_rapidity = TString("jet_rapidity_") +
		tag_levels_[level_index] + "_" +
		tag_jet_types_[jet1_index]; 
	  plots_ptr_ -> operator[](hash_key_rapidity)       = new TH1F(hash_key_rapidity, hash_key_rapidity + "; #eta; Events", 100, -5.2, 5.2);
	  const TString hash_key_eta = TString("jet_eta_") +
		tag_levels_[level_index] + "_" +
		tag_jet_types_[jet1_index]; 
	  plots_ptr_ -> operator[](hash_key_eta)       = new TH1F(hash_key_eta, hash_key_eta + "; #eta; Events", 100, -5.2, 5.2);
	  {
	    const TString hash_key = TString("jet_mass_norm_") +
	      tag_levels_[level_index] + "_" +
	      tag_jet_types_[jet1_index]; 
	    plots_ptr_ -> operator[](hash_key)       = new TH1F(hash_key, hash_key + "; mass_norm; Events", 100, 0, 1.1);
	  }
	  {
	    const TString hash_key = TString("jet_pt_norm_") +
	      tag_levels_[level_index] + "_" +
	      tag_jet_types_[jet1_index]; 
	    plots_ptr_ -> operator[](hash_key)       = new TH1F(hash_key, hash_key + "; pt_norm; Events", 100, 0, 1.1);
	  }
	  {
	    const TString hash_key = TString("jet_pz_norm_") +
	      tag_levels_[level_index] + "_" +
	      tag_jet_types_[jet1_index]; 
	    plots_ptr_ -> operator[](hash_key)       = new TH1F(hash_key, hash_key + "; pz_norm; Events", 100, -1.1, 1.1);
	  }
	  {
	    const TString hash_key = TString("jet_px_norm_") +
	      tag_levels_[level_index] + "_" +
	      tag_jet_types_[jet1_index]; 
	    plots_ptr_ -> operator[](hash_key)       = new TH1F(hash_key, hash_key + "; px_norm; Events", 100, -1.1, 1.1);
	  }
	  {
	    const TString hash_key = TString("jet_mass_") +
	      tag_levels_[level_index] + "_" +
	      tag_jet_types_[jet1_index]; 
	    plots_ptr_ -> operator[](hash_key)       = new TH1F(hash_key, hash_key + "; mass; Events", 100, 0.0, 30.0);
	  }
	  {
	    const TString hash_key = TString("jet_pt_") +
	      tag_levels_[level_index] + "_" +
	      tag_jet_types_[jet1_index]; 
	    plots_ptr_ -> operator[](hash_key)       = new TH1F(hash_key, hash_key + "; pt; Events", 100, -1.0, 200.0);
	  }
	  {
	    const TString hash_key = TString("jet_pz_") +
	      tag_levels_[level_index] + "_" +
	      tag_jet_types_[jet1_index]; 
	    plots_ptr_ -> operator[](hash_key)       = new TH1F(hash_key, hash_key + "; pz; Events", 100, - 250.0, 250.0);
	  }
	  {
	    const TString hash_key = TString("jet_px_") +
	      tag_levels_[level_index] + "_" +
	      tag_jet_types_[jet1_index]; 
	    plots_ptr_ -> operator[](hash_key)       = new TH1F(hash_key, hash_key + "; px; Events", 100, - 200.0, 200.0);
	  }



	  //printf("%s %s\n", hash_key_phi.Data(), hash_key_eta.Data());
	}
    }
    
  //  ConfigurePtRadiationProfileTool();

}
