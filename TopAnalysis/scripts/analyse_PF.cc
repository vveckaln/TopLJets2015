#include "TFile.h"
#include "TH1F.h"
#include "TApplication.h"
#include "Definitions.hh"
#include "interface/Definitions_cmssw.hh"
#include "TMath.h"
#include <map>
//#include <stdlib.h>
#include "TString.h"
#include <assert.h>
TFile * plotter_file = nullptr;
static const unsigned short N_MC_histo_names = 7;
static const char * MC_histo_names[N_MC_histo_names] = {"t#bar{t}", "Single top", "W", "DY", "Multiboson", "t#bar{t}+V", "QCD"};
TH1F * sum_MC(const char *);
enum model_t {SM, FLIP};
enum flow_t {N, E, PT};
const unsigned char nflows = 3;
const char* tag_flow[nflows] = {"N", "E", "Pt"};
const char * title_flow[nflows] = {"particle", "energy", "$P_{T}$"};
const char * model_tag[] = {"SM", "cflip"};
const char * model_title[] = {"SM", " colour octet $W$"};
const unsigned char npairs = 4;
const char * jetpairs[npairs] = {"blb2l", "hbqf", "qcqf", "hbqc"};
const char * jetpairstitle[npairs] = {"$j_{1}^{b}$,$j_{2}^{b}$", 
				      "$j_{h}^{b}$,$j_{f}^{W}$", 
				      "$j_{h}^{b}$,$j_{c}^{W}$", 
				      "$j_{1}^{W}$,$j_{2}^{W}$"};
const unsigned char nsources = 2;
const char * tag_source_type[nsources] = {"MC", "data"};
double R(TH1F *);
struct Result
{
  double value;
  double error[2][2]     = {{0.0, 0.0}, {0.0, 0.0}};
  bool   error_set[2][2] = {{false, false}, {false, false}};
  Result()
  {
    value           = 0.0;
  };
  void ls()
  {
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
  void normalise(float);
};
void prepare_header(FILE *, const char *, const char *);
void prepare_footer(FILE *);
void add_entry(FILE *, Result, Result);
Result CalculateResult(TH1F *, TH1F *, double (*)(TH1F*));
using namespace std;
using namespace Definitions;
int main(int argc, char *argv[])
{
  assert (argc == 2);
  const unsigned char model = stoi(argv[1]);
  const TString dir(TString("Rvalues_") + model_tag[model]);
  plotter_file = model == SM ? TFile::Open("$EOS/analysis_MC13TeV_TTJets/plots/plotter.root") : 
    TFile::Open("$EOS/analysis_MC13TeV_TTJets_cflip/plots/plotter.root");
  system(TString("mkdir -p ") + dir);
  system(TString("rm ") + dir + "/*");
  map<TString, FILE * > file_map;
  for (unsigned short ch_ind = L; ch_ind < N_channels_types_; ch_ind ++)
    {
      for (unsigned short level_ind = RECO; level_ind < N_levels_types_; level_ind ++)
	{
	  for (unsigned short type_ind = 0; type_ind < nsources; type_ind ++)
	    {
	      for (unsigned short flow_ind = 0; flow_ind < nflows; flow_ind ++)
		{
		  const TString file_name = dir + "/R_" + 
		    tag_channels_types_[ch_ind] + "_" + 
		    tag_levels_types_[level_ind] + "_" + 
		    tag_source_type[type_ind] + "_" + 
		    tag_flow[flow_ind] + "_" + 
		    model_tag[model] + ".txt";
		  file_map[file_name] = fopen(file_name.Data(), "w");
		  char caption[128];
		  sprintf(caption, "Value of $R'$ for the %s channel at %s level for %s, %s flow for the %s model.", channel_titles_[ch_ind], level_titles_[level_ind], tag_source_type[type_ind], title_flow[flow_ind], model_title[model]); 
		  char label[128];
		  sprintf(label, "R_%s_%s_%s_%s_%s", tag_channels_types_[ch_ind], tag_levels_types_[level_ind], tag_flow[type_ind], tag_source_type[type_ind], model_tag[model]); 
		  prepare_header(file_map[file_name], caption, label);
		}
	    }
	}
    }

  for (unsigned short ch_ind = L; ch_ind < N_channels_types_; ch_ind ++)
    {
      for (unsigned short level_ind = RECO; level_ind < N_levels_types_; level_ind ++)
	{
	  for (unsigned short type_ind = 0; type_ind < (model == SM ? 2 : 1); type_ind ++)
	    {
	      for (unsigned short flow_ind = 0; flow_ind < nflows; flow_ind ++)
		{
		  //	      printf("model %u type_ind %u tag_sources_types_ %s model_tag %s, \n", model, type_ind, tag_sources_types_[type_ind], model_tag[model]);
		  const TString file_name = dir + "/R_" + 
		    tag_channels_types_[ch_ind] + "_" + 
		    tag_levels_types_[level_ind] + "_" + 
		    tag_source_type[type_ind] + "_" + 
		    tag_flow[flow_ind] + "_" + 
		    model_tag[model] + ".txt";
		  FILE * table_file = file_map[file_name];
		  printf("Opening %s result %p\n", file_name.Data(), table_file);
		  Result Int[N_charge_types_][npairs];
		  Result Rvalue[N_charge_types_][npairs];
		  for (unsigned short pair_ind = 0; pair_ind < npairs; pair_ind ++)
		    {
		      for (unsigned short charge_ind = 0; charge_ind < N_charge_types_; charge_ind ++)
			{
			  TH1F * h_chi = NULL;
			  TH1F * h_chi_syst = NULL;
			  const TString dir = TString(tag_channels_types_[ch_ind]) + "_chi" + jetpairs[pair_ind] + 
			    "_" + tag_flow[flow_ind] + "_" + 
			    tag_charge_types_[charge_ind] + "_" + 
			    tag_levels_types_[level_ind] + "_jetprt";
			  printf("dir %s\n", dir.Data());
			  if (model == SM)
			    {
			      h_chi = type_ind == 0 ? sum_MC(dir) : (TH1F*) plotter_file -> GetDirectory(dir) -> Get(dir);
			      h_chi_syst = type_ind == 0 ? (TH1F *) plotter_file -> Get(dir + "/totalmcunc") : nullptr;
			    }
			  if (model == FLIP)
			    {
			      h_chi = (TH1F*) plotter_file -> GetDirectory(dir) -> Get(dir + "_t#bar{t} cflip");
			      h_chi_syst = (TH1F *) plotter_file -> Get(dir + "/totalmcunc") ;

			    }
			  if (h_chi_syst)
			    h_chi_syst -> Scale(1.0 / h_chi_syst -> Integral());
			  h_chi -> Scale(1.0 / h_chi -> Integral());
			  Int[charge_ind][pair_ind] = CalculateResult(h_chi, h_chi_syst, R);
			  if (flow_ind == 0 and charge_ind == 0)
			    printf("%f\n", Int[charge_ind][pair_ind].value);
			}

		    }
		  printf("*****");
		  for (unsigned short pair_ind = 0; pair_ind < npairs; pair_ind ++)
		    {
		      for (unsigned short charge_ind = 0; charge_ind < N_charge_types_; charge_ind ++)
			{
			  Rvalue[charge_ind][pair_ind] = Result(Int[charge_ind][pair_ind]);
			  
			  Rvalue[charge_ind][pair_ind].normalise(Int[charge_ind][3].value);
			}
		    }
		  for (unsigned short pair_ind = 0; pair_ind < npairs; pair_ind ++)
		    {
		      fprintf(table_file, "\t\t\\multicolumn{3}{|c|}{%s}\\\\\n", jetpairstitle[pair_ind]);
		      fprintf(table_file, "\t\t\\hline\n");
		      for (unsigned short charge_ind = 0; charge_ind < N_charge_types_; charge_ind ++)
			{
			  fprintf(table_file, "\t\t %s ", charge_titles_[charge_ind]);
			  if (level_ind == GEN and type_ind == DATA)
			    {
			      fprintf(table_file, "\t& x \t & x \\\\\n");
			      continue;
			    }
			  add_entry(table_file, Int[charge_ind][pair_ind], Rvalue[charge_ind][pair_ind]);
			  Rvalue[charge_ind][pair_ind].ls();

			}
		      fprintf(table_file, "\t\t\\hline\n");

		    }

		}

	    }
	}
    }
  for (unsigned short ch_ind = L; ch_ind < N_channels_types_; ch_ind ++)
    {
      for (unsigned short level_ind = RECO; level_ind < N_levels_types_; level_ind ++)
	{
	  for (unsigned short type_ind = 0; type_ind < nsources; type_ind ++)
	    {
	      for (unsigned short flow_ind = 0; flow_ind < nflows; flow_ind ++)
		{
		  const TString file_name = dir + "/R_" + 
		    tag_channels_types_[ch_ind] + "_" + 
		    tag_levels_types_[level_ind] + "_" + 
		    tag_source_type[type_ind] + "_" + 
		    tag_flow[flow_ind] + "_" + 
		    model_tag[model] + ".txt";
		  prepare_footer(file_map[file_name]);
		  fclose(file_map[file_name]);
		}
	    }
	}
    }


  plotter_file -> Close();
}

TH1F * sum_MC(const char * dir)
{
  printf("dir %s\n", dir);
  
  TH1F * h = (TH1F*) (plotter_file -> GetDirectory(dir) -> Get(TString(dir) + "_" +  MC_histo_names[0]) -> Clone(TString(dir) + "_MCsum"));
  printf("hname %s h int %f\n", h-> GetName(), h -> Integral()); 
  for (unsigned short MC_ind = 1; MC_ind < N_MC_histo_names; MC_ind ++)
    {
      h -> Add((TH1F*)  (plotter_file -> GetDirectory(dir) -> Get(TString(dir) + "_" + MC_histo_names[MC_ind])));

    }
  return h;
	       
}

double R(TH1F * h)
{
  TAxis * axis = h -> GetXaxis();
  return h -> Integral(axis -> FindBin(0.3), axis -> FindBin(0.7))/h -> Integral();
}

Result CalculateResult(TH1F *h, TH1F * h_syst, double ( *funcptr )(TH1F*))
{
  Result result;
  result.value = funcptr(h);
  //  static const char * error_directions[2] = {"HIGH", "LOW"};
  static const char   error_signs[2]      = {1, -1};
  //  static const char * error_type[2]       = {"STAT", "SYST"};
  for (unsigned short type_ind = 0; type_ind < (h_syst ? 2 : 1); type_ind ++)
    {
      for (unsigned short direction_ind = 0; direction_ind < 2; direction_ind ++)
	{

	  TH1F * h_clone = (TH1F*) h -> Clone("clone");
	  for (unsigned short bin_ind = 0; bin_ind < h -> GetNbinsX(); bin_ind ++)
	    {
	      const double bin_error = type_ind == 0 ? TMath::Sqrt(h -> GetBinContent(bin_ind)) : h_syst -> GetBinError(bin_ind);
	      
	      h_clone -> SetBinContent(bin_ind, h_clone -> GetBinContent(bin_ind) + error_signs[direction_ind] * bin_error); 
	    }
	  const double error = result.value - funcptr(h_clone);
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

 void Result::normalise(float norm)
 {
  value /= norm;
  for (unsigned char j = 0; j < 2; j ++)
    for (unsigned char k =0; k < 2; k ++)
      error[j][k] /= norm;
 }

void prepare_header(FILE * file, const char * caption, const char * label)
{
  fprintf(file, "\\begin{table}[htp]\n");
  fprintf(file, "\t\\begin{center}\n");
  fprintf(file, "\t\\caption{%s}\n", caption);
  fprintf(file, "\t\\label{tab:%s}\n", label);
  fprintf(file, "\t\t\\begin{tabular}{|l|cc|}\n");
  fprintf(file, "\t\t\\hline\n");
  fprintf(file, "\t\t\\makecell[c]{Jet constituents} & $I$ & $R^{-1}$ \\\\\n");
  fprintf(file, "\t\t\\hline\n");
  

}

void add_entry(FILE * table_file, Result Int, Result Rval)
{
  if (Rval.error_set[1][0] and Rval.error_set[1][1])
    {
      fprintf(table_file, "\t& %.2f \t$\\pm$ %.2f \t$\\pm$ %.2f \t& %.2f \t$\\pm$ %.2f \t$\\pm$ %.2f \\\\\n", 
	      Int.value, Int.GetAverageStats(), Int.GetAverageSyst(),
	      Rval.value, Rval.GetAverageStats(), Rval.GetAverageSyst()
	      );
    }
  else
    {
      fprintf(table_file, "\t& %.2f \t$\\pm$ %.2f \t& %.2f \t$\\pm$ %.2f \\\\\n", Int.value, Int.GetAverageStats(),
	      Rval.value, Rval.GetAverageStats());

    }
}

void prepare_footer(FILE *file)
{
  fprintf(file, "\t\t\\end{tabular}\n");
  fprintf(file, "\t\\end{center}\n");
  fprintf(file, "\\end{table}\n");

}
