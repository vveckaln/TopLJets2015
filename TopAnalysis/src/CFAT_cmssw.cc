#include "TopLJets2015/TopAnalysis/interface/CFAT_cmssw.hh"
#include "TopLJets2015/TopAnalysis/interface/Definitions_cmssw.hh"
#include "TTree.h"
#include "TFile.h"

CFAT_cmssw::CFAT_cmssw():ColourFlowAnalysisTool()
{
  
  migration_tree_ = new TTree("migration", "migration");

  migration_tree_ -> SetDirectory(NULL);

  migration_tree_ -> Branch("pull_angle_reco", &pull_angle_reco_,  "pull_angle_reco/D");
  migration_tree_ -> Branch("pull_angle_gen",  &pull_angle_gen_,   "pull_angle_gen/D");
  migration_tree_ -> Branch("weight_reco",     &weight_reco_,      "weight_reco/D");
  migration_tree_ -> Branch("weight_gen",      &weight_gen_,       "weight_gen/D");
  
}

void CFAT_cmssw::ResetMigrationValues()
{

  pull_angle_reco_ = 0;
  pull_angle_reco_plottable_ = false;
  weight_reco_ = 0;
  pull_angle_gen_ = 0;
  pull_angle_gen_plottable_ = false;
  weight_gen_ = 0;
  printf("******* Reset Migration Values *********\n");
}

void CFAT_cmssw::StoreMigrationValues(double value)
{
  printf("Entering CFAT_cmssw::StoreMigrationValues\n");
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
  printf("work mode %u pull_angle_reco %f pull_angle_reco_plottable %s weight_reco %.9f pull_angle_gen %f pull_angle_gen_plottable %s weight_gen %.9f\n", work_mode_, pull_angle_reco_,
	 pull_angle_reco_plottable_ ? "true" : "false", weight_reco_, pull_angle_gen_, pull_angle_gen_plottable_ ? "true" : "false", weight_gen_);
  printf("Exiting CFAT_cmssw::StoreMigrationValues\n");

  //printf("mode %u, value %f, weight %.9f\n", work_mode_, value, GetEvent() -> weight_);
  //getchar();
}


void CFAT_cmssw::PlotMigrationValues()
{ 
  if (pull_angle_reco_plottable_ and pull_angle_gen_plottable_)
    {
      const double weight = (weight_reco_ + weight_gen_) * 0.5 ;
      Fill2D("Migration", pull_angle_reco_, pull_angle_gen_, weight);
      printf("PLOTTED pull_angle_reco_ %f, pull_angle_gen_ %f, weight %f", pull_angle_reco_, pull_angle_gen_, weight);
      migration_tree_ -> Fill();
    }
  
}

void CFAT_cmssw::Fill1D(const TString & key, double value) const                                                                                                                            
{                                                                                                                                                                                                         
  ((map<TString, TH1*>*)plots_ptr_) -> operator[](key) -> Fill(value, GetEvent() -> weight_);                                                                                                                                    
}


void CFAT_cmssw::Fill2D(const TString & key, double value_x, double value_y, double weight) const
{                                                                                                                                                                                                        
  ((map<TString, TH2*>*)plots2D_ptr_) -> operator[](key) -> Fill(value_x, value_y, weight);                                                                                                                        //printf("filling\n");
  //getchar();
}


void CFAT_cmssw::WriteMigrationTree()
{
  //migration_tree_ -> SetBranchAddress("pull_angle_gen", &y);
  migration_file_ = TFile::Open(file_tag_, "RECREATE");
  migration_tree_ -> SetDirectory(migration_file_);
  printf("migration tree collected %lld entries\n", migration_tree_ -> GetEntries());
  migration_tree_ -> Write();
  migration_file_ -> Close();
}

