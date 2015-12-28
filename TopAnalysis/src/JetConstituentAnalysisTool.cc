#include "TopLJets2015/TopAnalysis/interface/JetConstituentAnalysisTool.hh"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"

void JetConstituentAnalysisTool::AssignHistograms() const
{
  plots_ptr_ -> operator[]("JetPtMinusJetConstituentPt")     = new TH1F("JetPtMinusJetConstituentPt", "; Pt [GeV]; Events", 100, -15, 15);
  plots_ptr_ -> operator[]("JetChargeContent")               = new TH1F("JetChargeContent",           "; Percentage; Events", 100, 0, 1.2);
  plots_ptr_ -> operator[]("NoAllConstituents")              = new TH1F("NoAllConstituents",          "; No; Events", 301, -0.5, 300.5);
  plots_ptr_ -> operator[]("NoChargedConstituents")          = new TH1F("NoChargedConstituents",      "; No; Events", 301, -0.5, 300.5);
  plots2D_ptr_ -> operator[]("ConstituentFlavours")          = new TH2F("ConstituentFlavours",        "; Flavour; Charge", 41, -20.5, 20.5, 7, -1 - 1/6, 1 + 1/6);
  plots2D_ptr_ -> operator[]("ConstituentFlavoursLightJets") = new TH2F("ConstituentFlavoursLightJets", "; Flavour; Charge", 41, -20.5, 20.5, 7, -1 - 1/6, 1 + 1/6);
  plots2D_ptr_ -> operator[]("ConstituentFlavoursBJets") = new TH2F("ConstituentFlavoursBJets",       "; Flavour; Charge", 41, -20.5, 20.5, 7, -1 - 1/6, 1 + 1/6);
  
}

void JetConstituentAnalysisTool::AnalyseAllJets() const
{
  float jet_const_total_pt = 0;  
  float jet_chconst_total_pt = 0;
  unsigned short NoAllConstituents = 0;
  unsigned short NoChargedConstituents = 0;
  for (int jet_const_index = 0; jet_const_index < event_ptr_ -> npf; jet_const_index ++)
    {
      if (event_ptr_ -> pf_j[jet_const_index] != index_)
	continue;
      const float jet_const_pt = sqrt(
	event_ptr_->pf_px[jet_const_index] * event_ptr_->pf_px[jet_const_index] +
	event_ptr_->pf_py[jet_const_index] * event_ptr_->pf_py[jet_const_index]);
      jet_const_total_pt += jet_const_pt;
      NoAllConstituents ++;
      plots2D_ptr_ -> operator[]("ConstituentFlavours") -> Fill(event_ptr_ -> pf_id[jet_const_index], event_ptr_ -> pf_charge[jet_const_index], weight_);
      
      if (event_ptr_ -> pf_charge[jet_const_index] != 0)
	{
	  jet_chconst_total_pt += jet_const_pt;
	  NoChargedConstituents ++;
	}
    }
  plots_ptr_ -> operator[]("JetPtMinusJetConstituentPt") -> Fill(jet_const_total_pt - jet_ptr_ -> Pt(), weight_);
  const float charge_content = jet_chconst_total_pt/jet_const_total_pt;
  plots_ptr_ -> operator[]("JetChargeContent")           -> Fill(charge_content, weight_);
  plots_ptr_ -> operator[]("NoAllConstituents")          -> Fill(NoAllConstituents, weight_);
  plots_ptr_ -> operator[]("NoChargedConstituents")      -> Fill(NoChargedConstituents, weight_);

}

void JetConstituentAnalysisTool::AnalyseLightJets() const
{
  for (int jet_const_index = 0; jet_const_index < event_ptr_ -> npf; jet_const_index ++)
    {
      if (event_ptr_ -> pf_j[jet_const_index] != index_)
	continue;
      plots2D_ptr_ -> operator[]("ConstituentFlavoursLightJets") -> Fill(event_ptr_ -> pf_id[jet_const_index], event_ptr_ -> pf_charge[jet_const_index], weight_);
    }
}

void JetConstituentAnalysisTool::AnalyseBJets() const
{
  for (int jet_const_index = 0; jet_const_index < event_ptr_ -> npf; jet_const_index ++)
    {
      if (event_ptr_ -> pf_j[jet_const_index] != index_)
	continue;
      plots2D_ptr_ -> operator[]("ConstituentFlavoursBJets") -> Fill(event_ptr_ -> pf_id[jet_const_index], event_ptr_ -> pf_charge[jet_const_index], weight_);
    }
}

void JetConstituentAnalysisTool::Do()
{
}
