#include "TopLJets2015/TopAnalysis/interface/CFAT_cmssw.hh"
#include "TopLJets2015/TopAnalysis/interface/Definitions_cmssw.hh"
#include "Definitions.hh"

#include "TTree.h"
#include "TFile.h"
#include <iostream>
CFAT_cmssw::CFAT_cmssw():ColourFlowAnalysisTool()
{
  for (ChannelCode_t channel_code = 1; channel_code < N_channels_types_; channel_code ++)
    {
      for (ChargeCode_t charge_code = 0; charge_code < 2; charge_code ++)
	{
	  for (unsigned short jet1_iter = 0; jet1_iter < 4; jet1_iter ++)
	    {
	      const TString name = 
		TString(tag_channels_types_[channel_code]) + "_" + 
		tag_charge_types_[charge_code] + "_" + 
		tag_jet_types_[jet1_iter] + "_migration";
	      migration_tree_[channel_code - 1][charge_code][jet1_iter] = 
		new TTree(name, name);
	      migration_tree_[channel_code - 1][charge_code][jet1_iter] -> SetDirectory(NULL);
	      migration_tree_[channel_code - 1][charge_code][jet1_iter] -> 
		Branch("pull_angle_reco", &pull_angle_[0][charge_code][jet1_iter],  "pull_angle_reco/F");
	      migration_tree_[channel_code - 1][charge_code][jet1_iter] -> 
		Branch("pull_angle_gen",  &pull_angle_[1][charge_code][jet1_iter],   "pull_angle_gen/F");
	      migration_tree_[channel_code - 1][charge_code][jet1_iter] -> 
		Branch("pvmag_reco",      &pvmag_[0][charge_code][jet1_iter],        "pvmag_reco/F");
	      migration_tree_[channel_code - 1][charge_code][jet1_iter] -> 
		Branch("pvmag_gen",       &pvmag_[1][charge_code][jet1_iter],         "pvmag_gen/F");
	      migration_tree_[channel_code - 1][charge_code][jet1_iter] -> Branch("weight",     &weights_);
	    }
	}
    }
}

void CFAT_cmssw::SetMigrationOutput(const char * migration_output)
{
  file_tag_ = TString(migration_output);
}

void CFAT_cmssw::ResetMigrationValues()
{
      for (ChargeCode_t charge_code = 0; charge_code < 2; charge_code ++)
	{
	  for (unsigned short jet1_iter = 0; jet1_iter < 4; jet1_iter ++)
	    {
	      fill_[charge_code][jet1_iter] = false;
	      for (WorkCode_t level_code = 0; level_code < N_levels_types_; level_code ++)
		{
		  pull_angle_[level_code][charge_code][jet1_iter] = - 10.0;
		  pvmag_[level_code][charge_code][jet1_iter] = - 1.0;
		}
	    }
    }
}

void CFAT_cmssw::StoreMigrationValues(ChargeCode_t chargecode, VectorCode_t jetcode, double pa, double mag)
{
  weights_ = GetEvent() -> weights_;
  printf("storing %lu\n", GetEvent() -> weights_.size());
  fill_[chargecode][jetcode] = true;
  pull_angle_[work_mode_][chargecode][jetcode] = pa;
  pvmag_[work_mode_][chargecode][jetcode] = mag;
}


void CFAT_cmssw::PlotMigrationValues()
{ 
  for (ChargeCode_t charge_code = 0; charge_code < 2; charge_code ++)
    {
      for (unsigned short jet1_iter = 0; jet1_iter < 4; jet1_iter ++)
	{
	  if (fill_[charge_code][jet1_iter])
	    {
	      migration_tree_[channel_code_ - 1][charge_code][jet1_iter] -> Fill();
	      if (charge_code == 1 and jet1_iter == 3)
		{
		  //		  printf("filling values wgt %.9f pa %f pvmag %f\n", weight_, pull_angle_[0][1][3], pvmag_[0][1][3]);
		  //getchar();
		}
	    }
	}
    }
}


void CFAT_cmssw::Fill1D(const TString & key, double value, double weight) const                                                                                                                        
{                                                                                                                                                                                                       
  {
    pair<TH1*, TH2*> * p = ((map<TString, pair<TH1*, TH2*> *>*)plots_ptr_) -> operator[](TString(tag_channels_types_[channel_code_]) + "_" + key);
    if (TString(tag_channels_types_[channel_code_]) + "_" + key = "L_pull_angle_allconst_reco_leading_jet_scnd_leading_jet_DeltaRTotal")
      {
	printf("plotting %lu \n", weights_.size()); 
      }
    p -> first -> Fill(value, GetEvent() -> weights_[0] * weight);
    //printf("%u\n", Definitions::nsyst_);
    if (Definitions::nsyst_ > 0) 
      {
	if (GetEvent() -> weights_.size() > nsyst_)
	  cout << "WARNING: Size of uncertainty weight vector larger than uncertainty histogram size." << endl;
	p -> second ->Fill(value, 0.0, GetEvent() -> weights_[0] * weight);
	for (unsigned int i = 1; i < GetEvent () -> weights_.size(); ++i) 
	  {
	    p -> second -> Fill(value, i, GetEvent() -> weights_[0] * GetEvent() -> weights_[i] * weight);
	  }
      }
  }
  {
    pair<TH1*, TH2*> * p = ((map<TString, pair<TH1*, TH2*>* >*)plots_ptr_) -> operator[](TString(tag_channels_types_[L]) + "_" + key);
    p -> first -> Fill(value, GetEvent() -> weights_[0] * weight);
    if (Definitions::nsyst_ > 0) 
      {
	if (GetEvent() -> weights_.size() > nsyst_)
	  cout << "WARNING: Size of uncertainty weight vector larger than uncertainty histogram size." << endl;
	p -> second ->Fill(value, 0.0, GetEvent() -> weights_[0] * weight);
	for (unsigned int i = 1; i < GetEvent() -> weights_.size(); ++i) 
	  {
	    p -> second -> Fill(value, i, GetEvent() -> weights_[0] * GetEvent() -> weights_[i] * weight);
	  }
      }
  }
}

void CFAT_cmssw::Fill2D(const TString & key, double value_x, double value_y, double weight) const
{                                                                                                                                                                                                         ((map<TString, TH2*>*)plots2D_ptr_) -> operator[](TString(tag_channels_types_[channel_code_]) + "_" + key) -> Fill(value_x, value_y, weight);                                                               ((map<TString, TH2*>*)plots2D_ptr_) -> operator[](TString(tag_channels_types_[L]) + "_" + key) -> Fill(value_x, value_y, weight);                                                                                                             
}


void CFAT_cmssw::WriteMigrationTree()
{
  //migration_tree_ -> SetBranchAddress("pull_angle_gen", &y);
  migration_file_ = TFile::Open(file_tag_, "RECREATE");
  printf("Opened migration file %p %s\n", migration_file_, file_tag_.Data());
  for (ChannelCode_t channel_code = 1; channel_code < N_channels_types_; channel_code ++)
    {
      for (ChargeCode_t charge_code = 0; charge_code < 2; charge_code ++)
	{
	  for (unsigned short jet1_iter = 0; jet1_iter < 4; jet1_iter ++)
	    {
	      TTree * tree = migration_tree_[channel_code - 1][charge_code][jet1_iter];
	      tree -> SetDirectory(migration_file_);
	      printf("%s collected %lld entries\n", tree -> GetName(), tree -> GetEntries());
	      tree -> Write();
	    }
	}
    }
  migration_file_ -> Close();
}

