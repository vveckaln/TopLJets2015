#include "TopLJets2015/TopAnalysis/interface/CFAT_cmssw.hh"
#include "TopLJets2015/TopAnalysis/interface/Definitions_cmssw.hh"
#include "Definitions.hh"

#include "TTree.h"
#include "TFile.h"

CFAT_cmssw::CFAT_cmssw():ColourFlowAnalysisTool()
{
  
  migration_tree_[L] = NULL;
  for (unsigned char migration_tree_ind = 1; migration_tree_ind < N_channels_types_; migration_tree_ind ++)
    {
      migration_tree_[migration_tree_ind] = new TTree(TString(tag_channels_types_[migration_tree_ind]) + "_migration", TString(tag_channels_types_[migration_tree_ind]) + "_migration");
      migration_tree_[migration_tree_ind] -> SetDirectory(NULL);
      migration_tree_[migration_tree_ind] -> Branch("pull_angle_reco", &pull_angle_reco_,  "pull_angle_reco/D");
      migration_tree_[migration_tree_ind] -> Branch("pull_angle_gen",  &pull_angle_gen_,   "pull_angle_gen/D");
      migration_tree_[migration_tree_ind] -> Branch("weight_reco",     &weight_reco_,      "weight_reco/D");
      migration_tree_[migration_tree_ind] -> Branch("weight_gen",      &weight_gen_,       "weight_gen/D");
    }
}

void CFAT_cmssw::SetMigrationOutput(const char * migration_output)
{
  file_tag_ = TString(migration_output);
}

void CFAT_cmssw::ResetMigrationValues()
{

  pull_angle_reco_           = 0;
  pull_angle_reco_plottable_ = false;
  weight_reco_               = 0;
  pull_angle_gen_            = 0;
  pull_angle_gen_plottable_  = false;
  weight_gen_                = 0;
}

void CFAT_cmssw::StoreMigrationValues(double value)
{
  if (work_mode_ == RECO)
    {
      pull_angle_reco_ = value;
      pull_angle_reco_plottable_ = true;
      weight_reco_ = GetEvent() -> weight_;
    }
  if (work_mode_ == GEN)
    {
      pull_angle_gen_ = value;
      pull_angle_gen_plottable_ = true;
      weight_gen_ = GetEvent() -> weight_;
    }

  //printf("mode %u, value %f, weight %.9f\n", work_mode_, value, GetEvent() -> weight_);
  //getchar();
}


void CFAT_cmssw::PlotMigrationValues()
{ 
  const TString histoname = "Migration";
  if (pull_angle_reco_plottable_ and pull_angle_gen_plottable_)
    {
      const double weight = (weight_reco_ + weight_gen_) * 0.5 ;
      Fill2D(histoname, pull_angle_reco_, pull_angle_gen_, weight);
      //      printf("PLOTTED pull_angle_reco_ %f, pull_angle_gen_ %f, weight %.9f\n", pull_angle_reco_, pull_angle_gen_, weight);
      migration_tree_[channel_code_] -> Fill();
    }
  if (not pull_angle_reco_plottable_ and pull_angle_gen_plottable_)
    {
      pull_angle_reco_ = - 10;
      const double weight = (weight_reco_ + weight_gen_) * 0.5 ;
      Fill2D(histoname, pull_angle_reco_, pull_angle_gen_, weight);
      //printf("PLOTTED pull_angle_reco_ %f, pull_angle_gen_ %f, weight %.9f\n", pull_angle_reco_, pull_angle_gen_, weight);
      migration_tree_[channel_code_] -> Fill();

    }
  if (pull_angle_reco_plottable_ and not pull_angle_gen_plottable_)
    {
      pull_angle_gen_ = - 10;
      const double weight = (weight_reco_ + weight_gen_) * 0.5 ;
      Fill2D(histoname, pull_angle_reco_, pull_angle_gen_, weight);
      //printf("PLOTTED pull_angle_reco_ %f, pull_angle_gen_ %f, weight %.9f\n", pull_angle_reco_, pull_angle_gen_, weight);
      migration_tree_[channel_code_] -> Fill();

    }
  
}

void CFAT_cmssw::Fill1D(const TString & key, double value, double weight) const                                                                                                                        
{                                                                                                                                                                                                       
  ((map<TString, TH1*>*)plots_ptr_) -> operator[](TString(tag_channels_types_[channel_code_]) + "_" + key) -> Fill(value, GetEvent() -> weight_*weight);                                                             ((map<TString, TH1*>*)plots_ptr_) -> operator[](TString(tag_channels_types_[L]) + "_" + key) -> Fill(value, GetEvent() -> weight_*weight);                                                 
}

void CFAT_cmssw::Fill2D(const TString & key, double value_x, double value_y, double weight) const
{                                                                                                                                                                                                         ((map<TString, TH2*>*)plots2D_ptr_) -> operator[](TString(tag_channels_types_[channel_code_]) + "_" + key) -> Fill(value_x, value_y, weight);                                                               ((map<TString, TH2*>*)plots2D_ptr_) -> operator[](TString(tag_channels_types_[L]) + "_" + key) -> Fill(value_x, value_y, weight);                                                                                                             
}


void CFAT_cmssw::WriteMigrationTree()
{
  //migration_tree_ -> SetBranchAddress("pull_angle_gen", &y);
  migration_file_ = TFile::Open(file_tag_, "RECREATE");
  printf("Opened migration file %p %s\n", migration_file_, file_tag_.Data());
  for (unsigned char migration_tree_ind = 1; migration_tree_ind < N_channels_types_; migration_tree_ind ++)
    {
      migration_tree_[migration_tree_ind] -> SetDirectory(migration_file_);
      printf("migration tree collected %lld entries\n", migration_tree_[migration_tree_ind] -> GetEntries());
      migration_tree_[migration_tree_ind] -> Write();
    }
  migration_file_ -> Close();
}

