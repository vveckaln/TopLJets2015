#include "ColourFlowAnalysisTool.hh"
#include "Definitions.hh"
class TFile;
class TTree;
using namespace Definitions;
class CFAT_cmssw: public ColourFlowAnalysisTool
{
  float      weight_;
  float      pull_angle_[2][2][4]; //recolevel/charge/1stjet
  bool       fill_[2][4];
  float      pvmag_[2][2][4];
  TString    file_tag_;

public:
  CFAT_cmssw();
  TTree * migration_tree_[N_channels_types - 1][2][4];
  TFile * migration_file_;
  void SetMigrationOutput(const char *);
  void WriteMigrationTree();
  void ResetMigrationValues();
  void StoreMigrationValues(ChargeCode_t chargecode, VectorCode_t jetcode, double pa, double mag);
  void PlotMigrationValues();
  inline void Fill1D(const TString & key, double value, double = 1.0) const;
  inline void Fill2D(const TString & key, double value_x, double value_y, double weight) const;

};
