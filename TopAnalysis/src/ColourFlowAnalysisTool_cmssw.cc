#include "ColourFlowAnalysisTool.hh"
void ColourFlowAnalysisTool::Fill1D(const TString & key, double value) const                                                                                                                            
{                                                                                                                                                                                                         
  ((map<TString, TH1*>*)plots_ptr_) -> operator[](key) -> Fill(value, GetEvent() -> weight_);                                                                                                                                    
}
