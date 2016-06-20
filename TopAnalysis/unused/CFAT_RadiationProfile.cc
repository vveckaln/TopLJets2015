#include "TopLJets2015/TopAnalysis/interface/ColourFlowAnalysisTool.hh"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"

#include "TH2F.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"

void ColourFlowAnalysisTool::ConfigurePtRadiationProfileTool()
{
  if (PtRadiation_mode_ == 0)
    return;
  gROOT -> SetBatch(kTRUE);
  gStyle -> SetOptStat(0);
  canvas_ = new TCanvas("canvas", "canvas");
  PtRadProf_ = new TH2F("Pt_radiation_profile", "Pt_radiation_profile; X [a.u.]; Y [a.u]", 100, -0.2, 0.2, 50, -0.2, 0.2);
  //  PtRadProf_ -> Draw();
  
}

void ColourFlowAnalysisTool::PtRadiationProfile() const
{
  if (PtRadiation_mode_ == 0)
    return;
  unsigned int size = work_mode_ == 0 ? event_ptr_ -> npf : event_ptr_ -> ngpf;
  const TLorentzVector * jets[2] = {leading_light_jet_ptr_, second_leading_light_jet_ptr_};
  const float indices[2] = {leading_light_jet_index_, second_leading_light_jet_index_};
  const TVector2 jet_difference(jets[1] -> Phi() - jets[0] -> Phi(), jets[1] -> Rapidity() - jets[0] -> Rapidity());
  //const float mag = jet_difference.Mod();
  //  const TVector2 jet_difference_norm = jet_difference *2/mag;
  const float tan = jet_difference.Py()/jet_difference.Px();
  const TVector2 normal(-tan/sqrt(1 + tan*tan), 1/sqrt(1 + tan*tan));
  const char direction[2] = {1, -1};
  const float coordinates[2][2] = {{-0.05, 0}, {0.05, 0}};
  for (unsigned char jet_index = 0; jet_index < 2; jet_index ++)
    {
      const unsigned char index = indices[jet_index];
      for (unsigned int jet_const_index = 0; jet_const_index < size; jet_const_index ++)
	{
	  TLorentzVector constituent_4vector;
	  if (work_mode_ == 0)
	    {
	      if (event_ptr_ -> pf_j[jet_const_index] != index)
		continue;
	      /* if (OnlyChargedConstituents and event_ptr_ -> pf_charge[jet_const_index] == 0)
		 continue;*/
	      constituent_4vector.SetPtEtaPhiM(event_ptr_ -> pf_pt[jet_const_index], 
					       event_ptr_ -> pf_eta[jet_const_index], 
					       event_ptr_ -> pf_phi[jet_const_index], 
					       0);
	      const float Pt = constituent_4vector.Pt();
	      TVector2 jet_const_phi_eta = TVector2(constituent_4vector.Phi(), constituent_4vector.Rapidity());
	      const TVector2 jet_const_phi_eta_from_jet(jet_const_phi_eta.Px() - jets[jet_index] -> Phi(), 
							jet_const_phi_eta.Py() - jets[jet_index] -> Rapidity());
	      const float parallel_projection = direction[jet_index]*
		(jet_const_phi_eta_from_jet.Px()*jet_difference.Px() + jet_const_phi_eta_from_jet.Py()*jet_difference.Py())/(jet_difference.Mod());
	      const float perpendicular_projection = jet_const_phi_eta_from_jet.Px()*normal.Px() + jet_const_phi_eta_from_jet.Py()*normal.Py();
	      PtRadProf_ -> Fill(coordinates[jet_index][0] + parallel_projection, coordinates[jet_index][1] + perpendicular_projection, Pt*weight_);
	    }
	}
    }
  /*  printf("weight %.8f\n", weight_);
  if (weight_ < 0)
    printf("************************\n");
  */
}

void ColourFlowAnalysisTool::EndPtRadiationProfile() 
{
  if (not canvas_)
    return;
  canvas_ -> cd();
  PtRadProf_ -> Draw("COLZ");
  
  canvas_ -> SaveAs("Radiation_profile.png");
  delete canvas_;

}
