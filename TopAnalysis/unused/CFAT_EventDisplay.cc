#include "TopLJets2015/TopAnalysis/interface/ColourFlowAnalysisTool.hh"
#include "TEllipse.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TArrow.h"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"

void ColourFlowAnalysisTool::EventDisplay(const PullVector & pv, bool charged, const char * tag) const
{
  if (event_display_mode_ == 0)
    return;
  static const float radius_scale = 25.0; 
  const TVector2 jet_difference(second_leading_light_jet_ptr_ -> Phi() - leading_light_jet_ptr_ -> Phi(), second_leading_light_jet_ptr_ -> Rapidity() - leading_light_jet_ptr_ -> Rapidity()); 
  const double pull_angle = TMath::ACos((pv.phi_component * jet_difference.Px() + pv.eta_component * jet_difference.Py())/(pv.Mod() * jet_difference.Mod())); 
  const unsigned char charge_ind = charged ? 1 : 0;
  gStyle -> SetOptStat(0);
  gROOT -> SetBatch(kTRUE); 
  const char* histogram_name[4] = {"leading_jet", "2nd_leading_jet", "other", "total"};
  TLine pull_vector(leading_light_jet_ptr_ -> Rapidity(), 
		    leading_light_jet_ptr_ -> Phi(),
		    leading_light_jet_ptr_ -> Rapidity() + 150 * pv.eta_component,
		    leading_light_jet_ptr_ -> Phi() + 150 * pv.phi_component);
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
  sprintf(dir_name, "event_displays/event_%lu_%s_%s_DeltaR_%f_pull_angle_%f_%s", 
	  event_number_, 
	  tag_levels_[work_mode_], 
	  tag_charge_types_[charge_ind], 
	  DeltaR, 
	  pull_angle, 
	  tag);
  system(TString("mkdir ") + dir_name);
  for (unsigned char hist_ind = 0; hist_ind < 4; hist_ind ++)
    {
      char name[128];
      sprintf(name, "event_display_%s_event_%lu_%s_%s_#DeltaR_%f_pull_angle_%f_%s", 
	      histogram_name[hist_ind], 
	      event_number_, 
	      tag_levels_[work_mode_], 
	      tag_charge_types_[charge_ind],
	      DeltaR, 
	      pull_angle,
	      tag);
      TH2F * hist = new TH2F(name, TString(name) + "; Rapidity[a.u.]; #phi[rad]", 10, -5, 5, 10, -1.1*TMath::Pi(), 1.1*TMath::Pi());
      TH2F_ptr_collector.push_back(hist);
      canvas[hist_ind] -> cd();
      hist -> Draw("9");
    }
 
  const Color_t fill_colour[5] = {kBlue, kRed, kGreen, kMagenta, kCyan};
 
  TEllipse * leading_jet = new TEllipse(
					leading_light_jet_ptr_ -> Rapidity(), 
					leading_light_jet_ptr_ -> Phi(), 
					leading_light_jet_ptr_ -> Pt()/(5 * radius_scale));
  leading_jet -> SetLineColor(fill_colour[0]);
  canvas[0] -> cd();
  leading_jet -> Draw();
  canvas[3] -> cd();
  leading_jet -> Draw();


  TEllipse * second_leading_jet = new TEllipse(
					second_leading_light_jet_ptr_ -> Rapidity(), 
					second_leading_light_jet_ptr_ -> Phi(), 
					second_leading_light_jet_ptr_ -> Pt()/(5 * radius_scale));
  second_leading_jet -> SetLineStyle(kDotted);
  second_leading_jet -> SetLineColor(fill_colour[1]);
  canvas[1] -> cd();
  second_leading_jet -> Draw();
  canvas[3] -> cd();
  second_leading_jet -> Draw();
  TEllipse * b_jets[2];
  for (unsigned char ind = 0; ind < b_jets_ptr_ -> size(); ind ++)
    {
      const TLorentzVector * b_jet = &b_jets_ptr_ -> operator[](ind);
      b_jets[ind] = new TEllipse(b_jet -> Rapidity(), b_jet -> Phi(), b_jet -> Pt()/(5 * radius_scale));
      b_jets[ind] -> SetLineStyle(kDashed);
      b_jets[ind] -> SetLineColor(fill_colour[ind + 2] + 10); 
      canvas[2] -> cd();

      b_jets[ind] -> Draw();
      canvas[3] -> cd();
      b_jets[ind] -> Draw();
    }
  unsigned short jet_indices[4] = {(unsigned short) leading_light_jet_index_, (unsigned short) second_leading_light_jet_index_, b_jets_indices_ptr_ -> operator[](0), b_jets_indices_ptr_ -> operator [](1)};
  if (work_mode_ == 1)
    {
      jet_indices[0] = event_ptr_ -> j_g[jet_indices[0]];
      jet_indices[1] = event_ptr_ -> j_g[jet_indices[1]];
      jet_indices[2] = event_ptr_ -> j_g[jet_indices[2]];
      jet_indices[3] = event_ptr_ -> j_g[jet_indices[3]];
   
    }
  vector<TEllipse*> pointer_collector;
  unsigned int size = work_mode_ == 0 ? event_ptr_ -> npf : event_ptr_ -> ngpf;
  printf("Event display leading jet rap %f phi %f\n", leading_light_jet_ptr_ -> Rapidity(), leading_light_jet_ptr_ -> Phi()); 
  printf(" ********************** EVENT %lu ******************** \n", event_number_);
  for (unsigned char part_type_ind = 0; part_type_ind < 3; part_type_ind ++)
    {
      unsigned char index = 0;
      if (part_type_ind < 2)
	{
	  index = jet_indices[part_type_ind];
	}
      for (unsigned short jet_const_index = 0; jet_const_index < size; jet_const_index ++)
	{
	  unsigned short particle_index = work_mode_ == 0 ? event_ptr_ -> pf_j[jet_const_index] : event_ptr_ -> gpf_g[jet_const_index];
	  
	   if (part_type_ind < 2)
	    {
	      if (particle_index != index)
		{
		  continue;
		}
	    }
	  if (part_type_ind == 2)
	    {
	      if (particle_index == jet_indices[0] or particle_index == jet_indices[1])
		{
		  continue;
		}
	    }
	  float charge = work_mode_ == 0 ? event_ptr_ -> pf_c[jet_const_index] : event_ptr_ -> gpf_c[jet_const_index];
	  if (charged and charge == 0)
	    continue;
	  TLorentzVector component;
	  if (work_mode_ == 0)
	    {
	      component.SetPtEtaPhiM(event_ptr_ -> pf_pt[jet_const_index], 
				  event_ptr_ -> pf_eta[jet_const_index], 
				  event_ptr_ -> pf_phi[jet_const_index],
				  event_ptr_ -> pf_m[jet_const_index]
				  );
	    }
	  else
	    {
	      component.SetPtEtaPhiM(event_ptr_ -> gpf_pt[jet_const_index], 
				  event_ptr_ -> gpf_eta[jet_const_index], 
				  event_ptr_ -> gpf_phi[jet_const_index],
				  event_ptr_ -> gpf_m[jet_const_index]
				  );	    
	    }
	  const float Rapidity = component.Rapidity();
	  const float phi   = component.Phi();
	  if (part_type_ind == 0 and particle_index == index)
	    {

	      printf("EventDisplay index %u component eta %f phi%f\n", index, Rapidity, phi); 
	    }
	 
	  const float radius = component.Pt()/radius_scale;
	  TEllipse * ellipse = new TEllipse(Rapidity, phi, radius, radius);
	  pointer_collector.push_back(ellipse);
	  //printf("Type %u Drawing %f %f %f\n", part_type_ind, Rapidity, phi, radius);
	  Color_t sel_fill_colour = kWhite;
	  if (particle_index == jet_indices[0])
	    {
	      sel_fill_colour = fill_colour[0];
	    }
	  else if (particle_index == jet_indices[1])
	    {
	      sel_fill_colour = fill_colour[1];
	    }
	  else if (particle_index == jet_indices[2])
	    {
	      sel_fill_colour = fill_colour[2];
	    }
	  else if (particle_index == jet_indices[3])
	    {
	      sel_fill_colour = fill_colour[3];
	    } 
	  else 
	    {
	      sel_fill_colour = fill_colour[4];
	    }
	  ellipse -> SetFillColor(sel_fill_colour);
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
      sprintf(name, "%s/%s_%s_%s_%lu_DeltaR_%f_pull_angle_%f_%s", 
	      dir_name, 
	      tag_levels_[work_mode_], 
	      tag_charge_types_[charge_ind], 
	      histogram_name[ind], 
	      event_number_, 
	      DeltaR, 
	      pull_angle,
	      tag); 
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
