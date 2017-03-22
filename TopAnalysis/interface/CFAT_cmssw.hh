#include "ColourFlowAnalysisTool.hh"
class TFile;
class TTree;
class CFAT_cmssw: public ColourFlowAnalysisTool
{
  
  Double_t pull_angle_reco_;
  bool pull_angle_reco_plottable_;
  Double_t weight_reco_;
  Double_t pull_angle_gen_;
  bool pull_angle_gen_plottable_;
  Double_t weight_gen_;
public:
  CFAT_cmssw();
  TString file_tag_;
  TTree * migration_tree_;
  TFile * migration_file_;
  virtual void WriteMigrationTree();
  virtual void ResetMigrationValues();
  virtual void StoreMigrationValues(double);
  virtual void PlotMigrationValues();
  virtual inline void Fill1D(const TString & key, double value) const;
  virtual inline void Fill2D(const TString & key, double value_x, double value_y, double weight) const;

};
