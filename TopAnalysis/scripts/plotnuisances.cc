#include "TFile.h"
#include "JsonParser.hh"
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
#include "TStyle.h"
#include "TH2.h"
TH1F * desymetrise(TH1 *);
static const unsigned char Nnondedicated = 21;
static const char * nondedicated_titles[Nnondedicated] = {"pileup up", "pileup down",
							  "trig efficiency correction up", "trig efficiency correction down",
							  "sel efficiency correction up", "sel efficiency correction down",
							  "b fragmentation up", "b fragmentation down",
							  "Peterson Frag",
							  "semilep BR up", "semilep BR down",
							  "QCD scale 1", "QCD scale 2", "QCD scale 3", "QCD scale 4", "QCD scale 5", "QCD scale 6", "QCD scale 7", "QCD scale 8", "QCD scale 9", "QCD scale 10"};
static const char * nondedicated_names[Nnondedicated] = {"pileup_up", "pileup_down",
							 "trig_efficiency_correction_up", "trig_efficiency_correction_down",
							 "sel_efficiency_correction_up", "sel_efficiency_correction_down",
							 "b_fragmentation_up", "b_fragmentation_down",
							 "Peterson_frag",
							 "semilep_BR_up", "semilep_BR_down",
							 "QCD_scale_1", "QCD_scale_2", "QCD_scale_3", "QCD_scale_4", "QCD_scale_5", "QCD_scale_6", "QCD_scale_7", "QCD_scale_8", "QCD_scale_9", "QCD_scale_10"};
int main()
{
  gStyle -> SetOptStat(0);
  gStyle -> SetOptTitle(0);
  static const float xsec = 832.0;
  JsonParser parser;
  parser.Load("data/era2016/syst_samples.json");
  parser.Load("data/era2016/expsyst_samples.json");
  TFile * fnominal = TFile::Open("$EOS/analysis_MC13TeV_TTJets/HADDChunks/MC13TeV_TTJets.root");
  TH1F * hnominal = desymetrise((TH1F*) fnominal -> Get("L_pull_angle_allconst_reco_leading_jet_scnd_leading_jet_DeltaRTotal"));
  hnominal -> SetDirectory(nullptr);
  hnominal -> Scale(luminosity * xsec);
  hnominal -> SetLineColor(kRed);
  for (unsigned int ind = 0; ind < parser.GetSize(); ind ++)
    {
      TCanvas c;
      SampleDescriptor *sd = parser.GetSample(ind);
      if (not TString(sd -> GetTitle()).Contains("t#bar{t}"))
	continue;
      TString tag = sd -> GetTag();
      if (tag.Index("up") == tag.Length() - 2 or tag.Contains("m173v5"))
	continue;
      TString uptag;
      bool isdown = false;
      TString title = tag;
      TString uptitle = sd -> GetTitle();
      if (tag.Contains("dn"))
	{
	  uptag = tag;
	  uptag .ReplaceAll("dn", "up");
	  uptitle.ReplaceAll("dn", "up");
	  title = title.ReplaceAll("dn", "");
	  isdown = true;
	}
      else if (tag.Contains("_down"))
	{
	  uptag = tag;
	  uptag .ReplaceAll("_down", "_up");
	  uptitle.ReplaceAll("down", "up");
	  title = title.ReplaceAll("_down", "");
	  isdown = true;
	}
      else if (tag.Contains("m171v5"))
	{
	  uptag = tag;
	  uptag .ReplaceAll("m171v5", "m173v5");
	  uptitle.ReplaceAll("171.5", "173.5");
	  title = title.ReplaceAll("m171v5", "topmass");
	  isdown = true;
	}

      TFile *f = TFile::Open(TString("$EOS/analysis_MC13TeV_TTJets/HADDChunks/") + sd -> GetTag() + ".root");
      if (not f) continue;
      TH1F *h = desymetrise((TH1F*) f -> Get("L_pull_angle_allconst_reco_leading_jet_scnd_leading_jet_DeltaRTotal"));
      TFile * fup = nullptr;
      TH1F * hup = nullptr;
      THStack stack;
      stack.Add(h, "P");
      stack.Add(hnominal);
      TLegend legend(0.6, 0.6, 0.9, 0.9);
      legend.SetLineWidth(0);
      legend.SetFillStyle(0);
      legend.AddEntry(hnominal, "t#bar{t}");
      legend.AddEntry(h, sd -> GetTitle());
      if (isdown)
	{
	  h -> SetMarkerStyle(kOpenTriangleDown);
	  fup = TFile::Open(TString("$EOS/analysis_MC13TeV_TTJets/HADDChunks/") + uptag + ".root");
	  if (not fup) continue;
	  hup = desymetrise((TH1F*) fup -> Get("L_pull_angle_allconst_reco_leading_jet_scnd_leading_jet_DeltaRTotal"));
	  hup -> Scale(luminosity * sd -> GetXsec());
	  hup -> SetMarkerStyle(kOpenTriangleUp);
	  hup -> SetDirectory(nullptr);
	  legend.AddEntry(hup, uptitle);
	  stack.Add(hup, "P");
	  hup -> SetLineWidth(0);
	  fup -> Close();
	  
	}
      else
	h -> SetMarkerStyle(kPlus);
      h -> SetLineWidth(0);
      h -> Scale(luminosity * sd -> GetXsec());
      TH1 * hframe = (TH1*)h -> Clone(title);
      hframe -> Reset("ICE");
      hframe ->SetDirectory(nullptr);
      printf("tag %s xsec %f\n", sd -> GetTag(), sd -> GetXsec());
      h -> SetDirectory(nullptr);
      f -> Close();
      stack.ls();
      stack.SetMaximum(1.4 * stack.GetMaximum("nostack"));
      hframe -> SetMaximum(1.4 * stack.GetMaximum("nostack"));
      hframe -> Draw();
      stack.Draw("NOSTACKSAME");
      legend.Draw("SAME");
      c.SaveAs(TString("nuisanceplots_SM/") + title + ".png");
      //      delete hframe;
    }
  /*JsonParser parser_cflip;
  parser_cflip.Load("data/era2016/expsyst_samples_cflip.json");
  TFile * fnominal_cflip = TFile::Open("$EOS/analysis/MC13TeV_TTJets_cflip.root");
  TH1F * hnominal_cflip = (TH1F*) fnominal_cflip -> Get("L_pull_angle_allconst_reco_leading_jet_2nd_leading_jet_DeltaRTotal");
  hnominal_cflip -> SetDirectory(nullptr);
  hnominal_cflip -> Scale(luminosity * 832.0);
  hnominal_cflip -> SetLineColor(kRed);
  hnominal_cflip -> SetMinimum(0.0);
  fnominal_cflip -> Close();
  for (unsigned int ind = 0; ind < parser_cflip.GetSize(); ind ++)
    {
      TCanvas c;
      SampleDescriptor *sd = parser_cflip.GetSample(ind);
      if (not TString(sd -> _title).Contains("t#bar{t}"))
	continue;
      TFile *f = TFile::Open(TString("$EOS/analysis/") + sd -> _tag + ".root");
      if (not f) continue;
      TH1F *h = (TH1F*) f -> Get("L_pull_angle_allconst_reco_leading_jet_2nd_leading_jet_DeltaRTotal");
      h -> Scale(luminosity * sd -> _xsec);
      printf("tag %s xsec %f\n", sd -> _tag, sd -> _xsec);
      h -> SetDirectory(nullptr);
      f -> Close();
      THStack stack;
      stack.Add(h);
      stack.Add(hnominal_cflip);
      stack.Draw("NOSTACK");
      TLegend legend(0.6, 0.6, 1.0, 1.0);
      legend.AddEntry(h, "nuisance");
      legend.AddEntry(hnominal, "ttbar cflip");
      legend.Draw("SAME");
      c.SaveAs(TString("nuisanceplots_cflip/") + sd -> _tag + ".png");
    }*/
  TH2 * h2D = (TH2F*) fnominal -> Get("L_pull_angle_allconst_reco_leading_jet_scnd_leading_jet_DeltaRTotal_syst"); 
  h2D -> SetDirectory(nullptr);
  fnominal -> Close();
  h2D -> Scale(xsec * luminosity);
  for (unsigned char bin_ind = 1; bin_ind < h2D -> GetNbinsY() + 1; bin_ind ++)
    {
      TCanvas c;
      const TString title = TString(nondedicated_names[bin_ind - 1]).ReplaceAll("_up", "").ReplaceAll("_1", "");
      TH1F * hup = desymetrise((TH1F*) h2D -> ProjectionX(nondedicated_names[bin_ind - 1], bin_ind, bin_ind));
      hup -> SetLineWidth(0);
      TH1 * hframe = (TH1*) hup -> Clone(TString(nondedicated_names[bin_ind-1]) + "_frame");
      hframe -> Reset("ICE");
      hframe -> SetDirectory(nullptr);
      THStack stack;
      stack.Add(hup, "P");
      stack.Add(hnominal);
      TLegend legend(0.6, 0.6, 0.9, 0.9);
      legend.SetLineWidth(0);
      legend.SetFillStyle(0);
      legend.AddEntry(hnominal, "t#bar{t}");
      legend.AddEntry(hup, nondedicated_titles[bin_ind - 1]);

      if (bin_ind == 1 or bin_ind == 3 or bin_ind == 5 or bin_ind == 7 or bin_ind == 10)
	{
	  hup -> SetMarkerStyle(kOpenTriangleUp);
	  hup -> SetDirectory(nullptr);
	  bin_ind ++;
	  TH1F * hdown = desymetrise((TH1F*) h2D -> ProjectionX(nondedicated_names[bin_ind - 1], bin_ind, bin_ind));
	      hdown -> Print("all");
	  hdown -> SetLineWidth(0);
	  hdown -> SetMarkerStyle(kOpenTriangleDown);
	  hdown -> SetDirectory(nullptr);
	  stack.Add(hdown, "P");
	  legend.AddEntry(hdown, nondedicated_titles[bin_ind - 1]);
	}
      else if (bin_ind == 9)
	{
	  hup -> SetMarkerStyle(kPlus);
	}
      else
	{
	  unsigned char ind = 0;
	  hup -> SetMarkerStyle(kPlus);
	  Style_t styles[] = {kStar, kCircle, kMultiply, 27, 28, 34, 30, 43, 48, 49 };
	  while (bin_ind < h2D -> GetNbinsY()  )  
	    {
	      bin_ind ++;
	      printf("bin_ind %u %u\n", bin_ind, h2D -> GetNbinsY());
	      TH1F * hqs = desymetrise((TH1F*) h2D -> ProjectionX(nondedicated_names[bin_ind - 1], bin_ind, bin_ind));
	      if (hqs -> Integral() == 0)
		continue;
	      hqs -> SetLineWidth(0);
	      hqs -> SetMarkerStyle(styles[ind]);
	      hqs -> SetDirectory(nullptr);
	      stack.Add(hqs, "P");
	      legend.AddEntry(hqs, nondedicated_titles[bin_ind - 1]);
	      ind ++;
	    } 
	}
      stack.SetMaximum(1.4 * stack.GetMaximum("nostack"));
      hframe -> SetMaximum(1.4 * stack.GetMaximum("nostack"));
      hframe -> Draw();
      stack.ls();
      stack.Draw("NOSTACKSAME");
      legend.Draw("SAME");
      c.SaveAs(TString("nuisanceplots_SM/") + title + ".png");

    }
}

TH1F * desymetrise(TH1 *h)
{
  TH1F * hret = new TH1F(TString(h-> GetName()) + "_desym", h -> GetTitle(), h -> GetNbinsX()/2, 0.0, h -> GetXaxis() -> GetXmax());
  for (unsigned char bin_ind = 1; bin_ind < h -> GetNbinsX()/2 + 1; bin_ind ++)
    {
      hret -> SetBinContent(bin_ind, h -> GetBinContent(h -> GetNbinsX()/2 - bin_ind + 1) + h -> GetBinContent(bin_ind + h -> GetNbinsX()/2));
    }
  hret -> Rebin(8);
  return hret;
}
