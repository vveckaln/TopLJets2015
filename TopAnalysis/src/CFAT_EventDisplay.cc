#include "TopLJets2015/TopAnalysis/interface/ColourFlowAnalysisTool.hh"
#include "TEllipse.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TArrow.h"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"

void ColourFlowAnalysisTool::EventDisplay(const PullVector & pv, float pull_angle, const TVector2 & jet_difference, bool charged) const
{
  if (event_display_mode_ == 0)
    return;
  const unsigned char charge_ind = charged ? 1 : 0;
  gStyle -> SetOptStat(0);
  gROOT -> SetBatch(kTRUE); 
  const char* histogram_name[4] = {"leading_jet", "2nd_leading_jet", "other", "total"};
  TLine pull_vector(leading_light_jet_ptr_ -> Rapidity(), 
		    leading_light_jet_ptr_ -> Phi(),
		    leading_light_jet_ptr_ -> Rapidity() + 250 * pv.eta_component,
		    leading_light_jet_ptr_ -> Phi() + 250 * pv.phi_component);
  TLine jet_dif(leading_light_jet_ptr_ -> Rapidity(), 
		leading_light_jet_ptr_ -> Phi(),
		leading_light_jet_ptr_ -> Rapidity() + jet_difference.Py(),
		leading_light_jet_ptr_ -> Phi() + jet_difference.Px());
  jet_dif.SetLineStyle(kDashed);
  TCanvas *canvas[4];
  for (unsigned char ind = 0; ind < 4; ind ++)
    {
      canvas[ind] = new TCanvas(histogram_name[ind], histogram_name[ind]);
    }
  
  float DeltaR = leading_light_jet_ptr_ -> DeltaR(*second_leading_light_jet_ptr_);
  vector<TH2F*> TH2F_ptr_collector;
  char dir_name[128];
  sprintf(dir_name, "event_displays/event_%lu_%s_%s_DeltaR_%f_pull_angle_%f", 
	  event_number_, 
	  tag_levels_[work_mode_], 
	  tag_charge_types_[charge_ind], 
	  DeltaR, 
	  pull_angle);
  system(TString("mkdir ") + dir_name);
  for (unsigned char hist_ind = 0; hist_ind < 4; hist_ind ++)
    {
      char name[128];
      sprintf(name, "event_display_%s_event_%lu_%s_%s_#DeltaR_%f_pull_angle_%f", 
	      histogram_name[hist_ind], 
	      event_number_, 
	      tag_levels_[work_mode_], 
	      tag_charge_types_[charge_ind],
	      DeltaR, 
	      pull_angle);
      TH2F * hist = new TH2F(name, TString(name) + "; Rapidity[a.u.]; #phi[rad]", 10, -5, 5, 10, -1.1*TMath::Pi(), 1.1*TMath::Pi());
      TH2F_ptr_collector.push_back(hist);
      canvas[hist_ind] -> cd();
      hist -> Draw("9");
    }
 
 
  TEllipse * leading_jet = new TEllipse(
					leading_light_jet_ptr_ -> Rapidity(), 
					leading_light_jet_ptr_ -> Phi(), 
					leading_light_jet_ptr_ -> Pt()/150);
  canvas[0] -> cd();
  leading_jet -> Draw();
  canvas[3] -> cd();
  leading_jet -> Draw();


  TEllipse * second_leading_jet = new TEllipse(
					second_leading_light_jet_ptr_ -> Rapidity(), 
					second_leading_light_jet_ptr_ -> Phi(), 
					second_leading_light_jet_ptr_ -> Pt()/150);
  second_leading_jet -> SetLineStyle(kDotted);
  canvas[1] -> cd();
  second_leading_jet -> Draw();
  canvas[3] -> cd();
  second_leading_jet -> Draw();
  TEllipse * b_jets[2];
  for (unsigned char ind = 0; ind < b_jets_ptr_ -> size(); ind ++)
    {
      const TLorentzVector * b_jet = &b_jets_ptr_ -> operator[](ind);
      b_jets[ind] = new TEllipse(b_jet -> Rapidity(), b_jet -> Phi(), b_jet -> Pt()/150);
      b_jets[ind] -> SetLineStyle(kDashed);
      canvas[2] -> cd();
      b_jets[ind] -> Draw();
      canvas[3] -> cd();
      b_jets[ind] -> Draw();
    }
  const Color_t fill_colour[3] = {kBlue, kRed, kGreen};
  const float ljet_indices[2] = {leading_light_jet_index_, second_leading_light_jet_index_};
  vector<TEllipse*> pointer_collector;
  unsigned int size = work_mode_ == 0 ? event_ptr_ -> npf : event_ptr_ -> ngen;
  
  printf(" ********************** EVENT %lu ******************** \n", event_number_);
  for (unsigned char part_type_ind = 0; part_type_ind < 3; part_type_ind ++)
    {
      unsigned char index = 0;
      if (part_type_ind < 2)
	index = ljet_indices[part_type_ind];
      
      for (unsigned short jet_const_index = 0; jet_const_index < size; jet_const_index ++)
	{
	  unsigned short particle_index = work_mode_ == 0 ? event_ptr_ -> pf_j[jet_const_index] : event_ptr_ -> g_j[jet_const_index];
	  if (part_type_ind < 2)
	    {
	      if (particle_index != index)
		{
		  continue;
		}
	    }
	  if (part_type_ind == 2)
	    {
	      if (particle_index == ljet_indices[0] or particle_index == ljet_indices[1])
		{
		  continue;
		}
	    }
	  float charge = work_mode_ == 0 ? event_ptr_ -> pf_charge[jet_const_index] : event_ptr_ -> g_charge[jet_const_index];
	  if (charged and charge == 0)
	    continue;
	  TLorentzVector component;
	  if (work_mode_ == 0)
	    {
	      const float energy = sqrt(
					pow(event_ptr_ -> pf_px[jet_const_index], 2) + 
					pow(event_ptr_ -> pf_py[jet_const_index], 2) + 
					pow(event_ptr_ -> pf_pz[jet_const_index], 2));
	      component = TLorentzVector(event_ptr_ -> pf_px[jet_const_index], 
					     event_ptr_ -> pf_py[jet_const_index], 
					     event_ptr_ -> pf_pz[jet_const_index],
					     energy);
	    }
	  else
	    {
	      const float energy = sqrt(
					pow(event_ptr_ -> g_px[jet_const_index], 2) + 
					pow(event_ptr_ -> g_py[jet_const_index], 2) + 
					pow(event_ptr_ -> g_pz[jet_const_index], 2));
	      component = TLorentzVector(event_ptr_ -> g_px[jet_const_index], 
					 event_ptr_ -> g_py[jet_const_index], 
					 event_ptr_ -> g_pz[jet_const_index],
					 energy);
	    }
	  const float Rapidity = component.Rapidity();
	  const float phi   = component.Phi();
	  const float radius = component.Pt()/150.0;
	  TEllipse * ellipse = new TEllipse(Rapidity, phi, radius, radius);
	  pointer_collector.push_back(ellipse);
	  //printf("Type %u Drawing %f %f %f\n", part_type_ind, Rapidity, phi, radius);
	  ellipse -> SetFillColor(fill_colour[part_type_ind]);
	  canvas[part_type_ind] -> cd();
	  ellipse -> Draw();
	  canvas[3] -> cd();
	  ellipse -> Draw();
	    
	}
    }
  canvas[3] -> cd();
  pull_vector.Draw();
  jet_dif.Draw();
  for (unsigned char ind = 0; ind < 4; ind ++)
    {
      char name[128];
      sprintf(name, "%s/%s_%s_%s_%lu_DeltaR_%f_pull_angle_%f", 
	      dir_name, 
	      tag_levels_[work_mode_], 
	      tag_charge_types_[charge_ind], 
	      histogram_name[ind], 
	      event_number_, 
	      DeltaR, 
	      pull_angle); 
      canvas[ind] -> SaveAs(TString(name) + ".pdf");
      canvas[ind] -> SaveAs(TString(name) + ".png");
    }
  for (unsigned short ind = 0; ind < pointer_collector.size(); ind ++)
    {
      delete pointer_collector[ind];
    }
  delete leading_jet; 
  delete second_leading_jet;
  delete b_jets[0];
  delete b_jets[1];
  for (unsigned short ind = 0; ind < TH2F_ptr_collector.size(); ind ++)
    {
      delete TH2F_ptr_collector[ind];
      delete canvas[ind];
    }
}
