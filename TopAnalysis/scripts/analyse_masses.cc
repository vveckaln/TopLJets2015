#include "TFile.h"
#include "TH1F.h"
#include "TApplication.h"
#include "Definitions.hh"
#include "interface/Definitions_cmssw.hh"
#include "TMath.h"
#include "TF1.h"
#include <map>
#include "TCanvas.h"
#include <string>
#include "TMarker.h"
#include "TROOT.h"
using namespace std;
using namespace Definitions;

static const unsigned short N_MC_histo_names = 7;
static const char * MC_histo_names[N_MC_histo_names] = {"t#bar{t}", "Single top", "W", "DY", "Multiboson", "t#bar{t}+V", "QCD"};
TH1F * sum_MC(const char *);
enum model_t {SM, FLIP};
unsigned char model;
TFile * plotter_file;// = model == SM ? TFile::Open("$EOS/analysis_24_12_2017/plots/plotter.root") : TFile::Open("$EOS/analysis_24_12_2017/MC13TeV_TTJets_cflip.root");
const char * model_tag[] = {"", "_flip"};
const char * model_title[] = {"", " for the colour octet $W$ model"};
const unsigned char N_interesting_jets = 4;
const unsigned char interesting_jets[N_interesting_jets] = {HAD_W, HAD_T, LEPT_W, LEPT_T};
double CalculateMass(TH1F *);
double CalculateWidth(TH1F *);
double determinant(unsigned int, void *);

struct Result
{
  char title[128];
  double value;
  double error[2][2]     = {{0.0, 0.0}, {0.0, 0.0}};
  bool   error_set[2][2] = {{false, false}, {false, false}};
  Result(const char * title)
  {
    sprintf(this -> title, "%s", title);
    value           = 0.0;
  };
  void ls()
  {
    printf("Listing %s\n", title);
    printf("value %f\n", value);
    printf("stat high %.3e, low %.3e\n", error[0][0], error[0][1]);
    printf("syst high %.3e, low %.3e\n", error[1][0], error[1][1]);

  };
  double GetAverageStats() const
  {
    return 0.5 * TMath::Abs(error[0][0] + error[0][1]);
  };
  double GetAverageSyst() const
  {
    return 0.5 * TMath::Abs(error[1][0] + error[1][1]);
  };
};
void prepare_header(FILE *, const char *, const char *);
void prepare_footer(FILE *);
void add_entry(FILE *, Result, Result);
Result CalculateResult(TH1F *, TH1F *, double (*)(TH1F*));

int main(int argc, char* argv[])
{
  gROOT -> SetBatch(kTRUE);
  if (argc < 2)
    throw "Included must be at least 1 program argument";
  model = stoi(argv[1]);
  plotter_file = model == SM ? TFile::Open("$EOS/analysis_MC13TeV_TTJets/plots/plotter.root") : TFile::Open("$EOS/analysis_MC13TeV_TTJets_cflip/plots/plotter.root");

  map<TString, FILE * > file_map;
  for (unsigned short ch_ind = L; ch_ind < N_channels_types_; ch_ind ++)
    {
      for (unsigned short level_ind = RECO; level_ind < N_levels_types_; level_ind ++)
	{
	  for (unsigned short type_ind = 0; type_ind < 2; type_ind ++)
	    {
	      const TString file_name = TString("mass_") + tag_channels_types_[ch_ind] + "_" + tag_levels_types_[level_ind] + "_" + tag_sources_types_[type_ind] + model_tag[model] + ".txt";
	      printf("creating %s\n", file_name.Data());
	      file_map[file_name] = fopen((TString("masses/tables/") + file_name).Data(), "w");
	      char caption[128];
	      sprintf(caption, "Observed masses of objects for the %s channel at %s level for %s%s.", channel_titles_[ch_ind], level_titles_[level_ind], title_sources_types_[type_ind], model_title[model]); 
	      char label[128];
	      sprintf(label, "mass_%s_%s_%s%s", tag_channels_types_[ch_ind], tag_levels_types_[level_ind], tag_sources_types_[type_ind], model_tag[model]); 
	      prepare_header(file_map[file_name], caption, label);
	    }
	}
    }
  TApplication app("myapp", 0, 0);
  for (unsigned short ch_ind = L; ch_ind < N_channels_types_; ch_ind ++)
    {
      for (unsigned short level_ind = RECO; level_ind < N_levels_types_; level_ind ++)
	{
	  for (unsigned short type_ind = 0; type_ind < (model == SM ? 2 : 1); type_ind ++)
	    {
	      for (unsigned short jet_ind = 0; jet_ind < N_interesting_jets; jet_ind ++)
		{	      //	      printf("model %u type_ind %u tag_sources_types_ %s model_tag %s, \n", model, type_ind, tag_sources_types_[type_ind], model_tag[model]);
		  const unsigned char interesting_jet_ind = interesting_jets[jet_ind];
		  printf("analysing jet %s\n", jet_titles_[interesting_jet_ind]);

		  const TString file_name = TString("mass_") + tag_channels_types_[ch_ind] + "_" + tag_levels_types_[level_ind] + "_" + tag_sources_types_[type_ind] + model_tag[model] + ".txt";
		  FILE * table_file = file_map[file_name];
		  printf("Opening %s result %p\n", file_name.Data(), table_file);
		  fprintf(table_file, "\t\t%s ", jet_titles_[interesting_jet_ind]);
		  if (level_ind == GEN and type_ind == DATA)
		    {
		      fprintf(table_file, "\t& x\\\\\n");
		      continue;
		    }
		  TH1F * h_chi = NULL;
		  TH1F * h_chi_syst = NULL;
		  const TString dir = TString(tag_channels_types_[ch_ind]) + "_jet_mass_" + tag_levels_types_[level_ind] + "_" + tag_jet_types_[interesting_jet_ind];

		  h_chi = type_ind == 0 ? sum_MC(dir) : (TH1F*) plotter_file -> GetDirectory(dir) -> Get(dir);
		  h_chi_syst = type_ind == 0 ? (TH1F *) plotter_file -> Get(dir + "/totalmcunc") : NULL;
		  Result mass = CalculateResult(h_chi, h_chi_syst, CalculateMass);
		  Result width = CalculateResult(h_chi, h_chi_syst, CalculateWidth);
		  add_entry(table_file, mass, width);
		  mass.ls();
		  width.ls();
		}
	    }
	}
    }
  for (unsigned short ch_ind = L; ch_ind < N_channels_types_; ch_ind ++)
    {
      for (unsigned short level_ind = RECO; level_ind < N_levels_types_; level_ind ++)
	{
	  for (unsigned short type_ind = 0; type_ind < 2; type_ind ++)
	    {

	      const TString file_name = TString("mass_") + tag_channels_types_[ch_ind] + "_" + tag_levels_types_[level_ind] + "_" + tag_sources_types_[type_ind] + model_tag[model] + ".txt";
	      prepare_footer(file_map[file_name]);
	      fclose(file_map[file_name]);
	    }
	}
    }

  app.Run(kTRUE);
  app.Terminate();

  plotter_file -> Close();

}

TH1F * sum_MC(const char * dir)
{
  TH1F * h = (TH1F*) (plotter_file -> GetDirectory(dir) -> Get(TString(dir) + "_" +  MC_histo_names[0]) -> Clone(TString(dir) + "_MCsum"));
  printf("hname %s h int %f\n", h-> GetName(), h -> Integral()); 
  for (unsigned short MC_ind = 1; MC_ind < N_MC_histo_names; MC_ind ++)
    {
      h -> Add((TH1F*)  (plotter_file -> GetDirectory(dir) -> Get(TString(dir) + "_" + MC_histo_names[MC_ind])));

    }
  return h;
	       
}

double CalculateMass(TH1F * h)
{
  const double y_max = h -> GetMaximum();
  int bin_max = 0;
  h -> GetBinWithContent(y_max, bin_max);
  const double x_max = h -> GetBinCenter(bin_max);
  //  printf("x_max %f y_max %f underflow %f\n", x_max, y_max, h -> GetBinContent(0));
  const double a = 1.0;
  const double b = -2*x_max;
  const double c = y_max + 2 * x_max * x_max;

  TF1 polynomial("polynomial", "[0]*x*x + [1]*x + [2]");
  polynomial.SetParameters(a, b, c);
  polynomial.SetLineColor(kGreen);
  TCanvas canvas(TString("mass_histo") + h -> GetName(), TString("mass_histo") + h -> GetTitle());
  h -> Draw();
  h -> Fit(& polynomial, "", "", x_max - 8.0, x_max + 8.0);
  double results[3];
  polynomial.GetParameters(results);
  const double mass = - results[1]/(2.0 * results[0]);
  const double peak = results[2] - results[1] * results[1]/(4 * results[0]);
  canvas.SaveAs(TString("masses/plots/") + h -> GetName() + ".png");
  return mass; 
}


double CalculateWidth(TH1F * h)
{
  enum side {LEFT, RIGHT};
  double x_FWHM[2];

  const double y_max = h -> GetMaximum();
  int bin_max = 0;
  h -> GetBinWithContent(y_max, bin_max);
  const double x_max = h -> GetBinCenter(bin_max);
  // printf("%s\n", h -> GetName());
  if (TString(h -> GetName()).Contains(tag_jet_types_[LEPT_W]))
    {
      for (unsigned int side_ind = LEFT; side_ind <= RIGHT; side_ind ++)
	{
	  unsigned short bin_ind = side_ind == LEFT ? 0 : h-> GetNbinsX(); 
	  do
	    {
	      if (side_ind == LEFT)
		{
		  bin_ind ++;
		}
	      else
		{ 
		  bin_ind --;
		}
	    }while (h -> GetBinContent(side_ind == LEFT ? bin_ind + 1 : bin_ind -1) < y_max/2.0);
	  printf("side_ind %u bin_ind %u low edge %f width %f\n", side_ind, bin_ind, h -> GetBinLowEdge(bin_ind), h-> GetBinWidth(bin_ind)); 
	  x_FWHM[side_ind] = side_ind == LEFT ? h -> GetBinLowEdge(bin_ind) + h -> GetBinWidth(bin_ind) : h -> GetBinLowEdge(bin_ind); 
	}
      const double width = x_FWHM[RIGHT] - x_FWHM[LEFT];
      printf("leptonic W width %f\n", width);
      
      return width;
    }
  //getchar();
  // printf("x_max %f y_max %f underflow %f\n", x_max, y_max, h -> GetBinContent(0));
  const double a = 1.0;
  const double b = -2*x_max;
  const double c = y_max + 2 * x_max * x_max;

  TF1 polynomial("polynomial", "[0]*x*x + [1]*x + [2]");
  polynomial.SetParameters(a, b, c);
  //  TCanvas canvas("canvas_width", "canvas_width");
  //h -> Draw();
  polynomial.SetLineColor(kBlue);
  h -> Fit(& polynomial, "", "", x_max - 8.0, x_max + 8.0);

  /*  canvas.Update();
      canvas.Modified();*/
  double results[3];
  polynomial.GetParameters(results);
  
  const double mass = - results[1]/(2.0 * results[0]);
  const double peak = results[2] - results[1] * results[1]/(4.0 * results[0]);
  //  printf("y_max %f peak %f mass %f\n", y_max, peak, mass);
  double y[3] = {peak*1.6/2.0, peak/2.0, peak*0.8/2.0};
  
  const char sign[] = {-1, 1};
  //TCanvas canvas;
  for (unsigned char side_ind = LEFT; side_ind <= RIGHT; side_ind ++)
    {
      double bin[3]; 
      double A[3][3];
      double coeff[3]; 
      double x[3];
      for (unsigned char ind = 0; ind < 3; ind ++)
	{
	  unsigned short bin_ind = side_ind == LEFT ? 0 : h -> GetNbinsX();
	  do
	    {
	      if (side_ind == LEFT)
		{
		  bin_ind ++;
		}
	      else
		{
		  bin_ind --;
		}
	    } while(h -> GetBinContent(side_ind == LEFT ? bin_ind + 1 : bin_ind - 1) < y[ind]);
	  //	  h -> GetBinWithContent(y[ind], bin, side_ind == LEFT ? 0 : bin_max);
	  x[ind] = side_ind == LEFT ?  h -> GetBinLowEdge(bin_ind) + h -> GetBinWidth(bin_ind) : h -> GetBinLowEdge(bin_ind);
	  if (ind == 2)
	    x[ind] += sign[side_ind]*5;
	  //	  printf("side %u x[ind] %f y[ind] %f \n", side_ind, x[ind], y[ind]);
	  A[ind][0] = x[ind] * x[ind];
	  A[ind][1] = x[ind];
	  A[ind][2] = 1.0;
	}
      double x_min = x[0] < x[2] ? x[0] : x[2];
      double x_max = x[0] < x[2] ? x[2] : x[0];
      for (unsigned char coeff_ind = 0; coeff_ind < 3; coeff_ind ++)
	{
	  double M[3][3];
	  for (unsigned char ind_M_row = 0; ind_M_row < 3; ind_M_row ++)
	    {
	      for (unsigned char ind_M_col = 0; ind_M_col < 3; ind_M_col ++)
		{
		  M[ind_M_row][ind_M_col] = (ind_M_col == coeff_ind) ? y[ind_M_row] : A[ind_M_row][ind_M_col];
		}
	    }
	  //printf("detM %f detA %f\n", determinant(3, M), determinant(3, A));
	  coeff[coeff_ind] = determinant(3, M)/determinant(3, A);
	}
      
      
      TF1 FWHM_polynomial("FWHM_polynomial", "[0]*x*x + [1]*x + [2]", x[0], x[2]);
      //  printf("before fit %f %f %f\n", coeff[0], coeff[1], coeff[2]);
      FWHM_polynomial.SetParameters(coeff[0], coeff[1], coeff[2]);
      FWHM_polynomial.SetLineColor(kGreen);
      //printf("check %f\n", coeff[0] * x[1]*x[1] + coeff[1]*x[1] + coeff[2]);

      /*if (side_ind == LEFT) 
       */

      //h -> Draw();
      
      FWHM_polynomial.Draw("SAME");
      h -> Fit(& FWHM_polynomial, "", "", x_min, x_max);
      /* canvas.Update();
	 canvas.Modified();*/
      double results[3];
      FWHM_polynomial.GetParameters(results);
      //printf("after fit %f %f %f\n", results[0], results[1], results[2]);
      //printf("after fit check %f\n", results[0] * x[1]*x[1] + results[1]*x[1] + results[2]);
      //printf("x_min %f x_max %f\n", x_min, x_max);
      const double a = results[0]; const double b = results[1]; const double c = results[2] - y[1];
      const double discriminant = b*b - 4 * a * c;
      const double x1 = (- b + TMath::Sqrt(discriminant)) / (2 * a);
      const double x2 = (- b - TMath::Sqrt(discriminant)) / (2 * a);
      if (x1 > x_min and x1 < x_max and not (x2 > x_min and x2 < x_max))
	x_FWHM[side_ind] = x1; 
      else if (not (x1 > x_min and x1 < x_max) and x2 > x_min and x2 < x_max)
	x_FWHM[side_ind] = x2; 
      else
	x_FWHM[side_ind] = x[1]; 

	//throw "Not in range";
      /*printf("side %u solution %f\n", side_ind, x_FWHM[side_ind]);
      TMarker marker(x_FWHM[side_ind], y[1], 2);
      marker.SetMarkerColor(kRed);
      marker.Draw("SAME");
      canvas.Update();
      canvas.Modified();
      
      getchar();
      */      
    }
  const double width = x_FWHM[RIGHT] - x_FWHM[LEFT];
  //printf("width %f\n", width);
  //getchar();
  return width; 
}

Result CalculateResult(TH1F *h, TH1F * h_syst, double ( *funcptr )(TH1F*))
{
  Result result("");
  result.value = funcptr(h);
  //  static const char * error_directions[2] = {"HIGH", "LOW"};
  static const char   error_signs[2]      = {1, -1};
  //  static const char * error_type[2]       = {"STAT", "SYST"};
  for (unsigned short type_ind = 0; type_ind < (h_syst ? 2 : 1); type_ind ++)
    {
      for (unsigned short direction_ind = 0; direction_ind < 2; direction_ind ++)
	{

	  TH1F * h_clone = (TH1F*) h -> Clone(TString("clone") + h -> GetName());
	  for (unsigned short bin_ind = 0; bin_ind < h -> GetNbinsX(); bin_ind ++)
	    {
	      const double bin_error = type_ind == 0 ? TMath::Sqrt(h -> GetBinContent(bin_ind)) : h_syst -> GetBinError(bin_ind);
	      //	      printf("type_ind %u, GetBinContent(bin_ind) %f, bin_error %f \n", type_ind, h -> GetBinContent(bin_ind), bin_error);
	      h_clone -> SetBinContent(bin_ind, h_clone -> GetBinContent(bin_ind) + error_signs[direction_ind] * bin_error); 
	    }
	  const double error = result.value - funcptr(h_clone);
	  printf("type_ind %u error %f\n", type_ind, error); 
	  unsigned short high_low = error > 0 ? 0 : 1; 
	  if (not result.error_set[type_ind][high_low])
	    {
	      result.error_set[type_ind][high_low] = true;
  	      result.error[type_ind][high_low] = error;
	    }
	  else
	    {
	      result.error[type_ind][high_low] = TMath::Sqrt(TMath::Power(error, 2) + TMath::Power(result.error[type_ind][high_low], 2));
	    }


	  delete h_clone;
	}
			       
    }
  return result;
}

void prepare_header(FILE * file, const char * caption, const char * label)
{
  fprintf(file, "\\begin{table}[htp]\n");
  fprintf(file, "\t\\begin{center}\n");
  fprintf(file, "\t\\caption{%s}\n", caption);
  fprintf(file, "\t\\label{tab:%s}\n", label);
  fprintf(file, "\t\t\\begin{tabular}{l|rr}\n");
  fprintf(file, "\t\t\\hline\n");
  fprintf(file, "\t\t\\makecell[c]{Object}&\\makecell[c]{Mass}&\\makecell[c]{FWHM}\\\\\n");
  fprintf(file, "\t\t\\hline\n");
}

void add_entry(FILE * table_file, Result mass_result, Result width_result)
{
  if (mass_result.error_set[1][0] and mass_result.error_set[1][1])
    {
      fprintf(table_file, "\t& %.3e \t$\\pm$ %.3e \t$\\pm$ %.3e \t& %.3e \t$\\pm$ %.3e \t$\\pm$ %.3e \\\\\n", mass_result.value, mass_result.GetAverageStats(), mass_result.GetAverageSyst(), width_result.value, width_result.GetAverageStats(), width_result.GetAverageSyst());
    }
  else
    {
      fprintf(table_file, "\t& %.3e \t$\\pm$ %.3e \t& %.3e \t$\\pm$ %.3e \\\\\n", mass_result.value, mass_result.GetAverageStats(), width_result.value, width_result.GetAverageStats());

    }
}

void prepare_footer(FILE *file)
{
  fprintf(file, "\t\t\\hline\n");
  fprintf(file, "\t\t\\end{tabular}\n");
  fprintf(file, "\t\\end{center}\n");
  fprintf(file, "\\end{table}\n");

}

double determinant(unsigned int N, void * pointer)
{
  double (* matrix)[N] = (double (*)[N]) pointer;
  double (* test)[N] = matrix;
  /*printf("listing matrix \n");
  for (unsigned char row_ind = 0; row_ind < N; row_ind ++)
    {
      for (unsigned char col_ind = 0; col_ind < N; col_ind ++)
	{
	  printf("%f ", *(*(matrix + row_ind) + col_ind)); 
	}
      printf("\n");
    }*/
  
  if (N == 1)
    return ** matrix;
  double det = 0.0;
  char coeff = 1;

  for (unsigned int col_ind = 0; col_ind < N; col_ind ++)
    {
      double minor[N - 1][N - 1];
      unsigned int minor_col_ind = 0;
      for (unsigned int minor_col_iter = 0; minor_col_iter < N; minor_col_iter ++)
	{
	  if (minor_col_iter == col_ind)
	    continue;

	  for (unsigned int minor_row_iter = 1; minor_row_iter < N; minor_row_iter ++)
	    {
	      const unsigned int minor_row_ind = minor_row_iter - 1;
	      minor[minor_row_ind][minor_col_ind] = *(*(matrix + minor_row_iter) + minor_col_iter); 
	      //printf("adding to minor %f minor_row_ind %u minor_col_ind %u \n", *(*(matrix + minor_row_iter) + minor_col_iter), minor_row_ind, minor_col_ind);
	    }
	  minor_col_ind ++;

	}    
      det += coeff * *(*matrix + col_ind) * determinant(N - 1, minor);
      // printf("coeff %d *(*matrix + ind) %f determinant(N - 1, minor) %f \n", coeff, *(*matrix + col_ind),  determinant(N - 1, minor));
      coeff *= -1;
      
    }
  return det;
}
