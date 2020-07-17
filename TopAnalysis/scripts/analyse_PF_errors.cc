//run as ./analyse_PF_errors 0 lx
#include "TFile.h"
#include "TH1F.h"
#include "TApplication.h"
#include "Definitions.hh"
#include "interface/Definitions_cmssw.hh"
#include "TMath.h"
#include <map>
#include "ERRORS.h"
//#include <stdlib.h>
#include "TString.h"
#include <assert.h>
using namespace std;
using namespace Definitions;

TFile * plotter_file = nullptr;
static const unsigned short N_MC_histo_names = 7;
static const char * MC_histo_names[N_MC_histo_names] = {"t#bar{t}", "Single top", "W", "DY", "Multiboson", "t#bar{t}+V", "QCD"};
TH1F * sum_MC(const char *);
enum model_t {SM, FLIP};
enum flow_t {N, E, PT};
const unsigned char nflows = 3;
const char* tag_flow[nflows] = {"N", "E", "Pt"};
const char * title_flow[nflows] = {"particle", "energy", "$p_{\\rm T}$"};
const char * model_tag[] = {"SM", "cflip"};
const char * model_title[] = {"SM", "colour octet \\PW"};
const unsigned char npairs = 4;
const char * jetpairs[npairs] = {"blb2l", "qfhb", "hbqc", "qlq2l"};
const char * jetpairstitle[npairs] = {"$j_{1}^{b}$,\\ $j_{2}^{b}$", 
				      "$j_{\\text{f}}^{W}$,\\ $j_{\\text{h}}^{b}$", 
				      "$j_{\\text{h}}^{b}$,\\ $j_{\\text{c}}^{W}$", 
				      "$j_{1}^{W}$,\\ $j_{2}^{W}$"};
const unsigned char nsources = 2;
const char * tag_source_type[nsources] = {"MC", "data"};
double R(TH1 *, void *);
void prepare_header(FILE *, const char *, const char *, SourceCode_t);
void prepare_footer(FILE *);
void add_entry(FILE *, Result);
const char * env = nullptr;
TH1 * join (TH1* hintra, TH1 * hinter);
int main(int argc, char *argv[])
{
  assert (argc == 3);
  const unsigned char model = stoi(argv[1]);
  env = argv[2];
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
		  char caption[256];
		  sprintf(caption, "Values of $R$ for the %s channel at %s level for %s flow in %s for the %s model.", channel_titles_[ch_ind], level_titles_[level_ind], title_flow[flow_ind], tag_source_type[type_ind], model_title[model]); 
		  char label[128];
		  sprintf(label, "R_%s_%s_%s_%s_%s", tag_channels_types_[ch_ind], tag_levels_types_[level_ind], tag_flow[flow_ind], tag_source_type[type_ind], model_tag[model]); 
		  prepare_header(file_map[file_name], caption, label, type_ind);
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
		  printf("model %u type_ind %u tag_sources_types_ %s model_tag %s, \n", model, type_ind, tag_sources_types_[type_ind], model_tag[model]);
		  const TString file_name = dir + "/R_" + 
		    tag_channels_types_[ch_ind] + "_" + 
		    tag_levels_types_[level_ind] + "_" + 
		    tag_source_type[type_ind] + "_" + 
		    tag_flow[flow_ind] + "_" + 
		    model_tag[model] + ".txt";
		  FILE * table_file = file_map[file_name];
		  printf("Opening %s result %p\n", file_name.Data(), table_file);
		  Result Rvalue[N_charge_types_][npairs - 1];
		  TH1F * h_chi[npairs][N_charge_types_];
		  TH1F * h_chi_syst[npairs][N_charge_types_];

		  for (unsigned short pair_ind = 0; pair_ind < npairs; pair_ind ++)
		    {
		      for (unsigned short charge_ind = 0; charge_ind < N_charge_types_; charge_ind ++)
			{
			  h_chi[pair_ind][charge_ind] = nullptr;
			  h_chi_syst[pair_ind][charge_ind] = nullptr;
			  const TString dir = TString(tag_channels_types_[ch_ind]) + "_chi" + jetpairs[pair_ind] + 
			    "_" + tag_flow[flow_ind] + "_" + 
			    tag_charge_types_[charge_ind] + "_" + 
			    tag_levels_types_[level_ind] + "_jetprt";
			  //			  printf("dir %s\n", dir.Data());
			  if (model == SM)
			    {
			      h_chi[pair_ind][charge_ind] = type_ind == 0 ? sum_MC(dir) : (TH1F*) plotter_file -> GetDirectory(dir) -> Get(dir);
			      h_chi_syst[pair_ind][charge_ind] = type_ind == 0 ? (TH1F *) plotter_file -> Get(dir + "/totalmcunc") : nullptr;
			    }
			  if (model == FLIP)
			    {
			      h_chi[pair_ind][charge_ind] = (TH1F*) plotter_file -> GetDirectory(dir) -> Get(dir + "_t#bar{t} cflip");
			      h_chi_syst[pair_ind][charge_ind] = (TH1F *) plotter_file -> Get(dir + "/totalmcunc") ;

			    }
			  // if (h_chi_syst)
			  //   {
			  //     printf("%f %f\n", h_chi_syst -> Integral(), h_chi -> Integral());
			  //     h_chi_syst -> Scale(1.0 / h_chi_syst -> Integral());
			  //     getchar();
			  //   }
			  //			  h_chi -> Scale(1.0 / h_chi -> Integral());
			}

		    }
		  printf("*****\n");
		      for (unsigned short charge_ind = 0; charge_ind < N_charge_types_; charge_ind ++)
			{
			  for (unsigned short pair_ind = 0; pair_ind < npairs - 1; pair_ind ++)
			    {
			      if (type_ind == DATA and level_ind == GEN)
				continue;
			      const TString dir = TString(tag_channels_types_[ch_ind]) + "_chi" + jetpairs[pair_ind] + 
				"_" + tag_flow[flow_ind] + "_" + 
				tag_charge_types_[charge_ind] + "_" + 
				tag_levels_types_[level_ind] + "_jetprt";
			      // if (not (type_ind == 0 and pair_ind == 0 and  charge_ind == 0 and ch_ind == 0 and level_ind == GEN and flow_ind == 0))
			      // 	continue;
			      //type 0, pair 2 charge 0 ch 2 level 1 flow 2
			      //type 0, pair 0 charge 0 ch 0 level 0 flow 1
			      
			      // if (dir != "L_chiblb2l_E_allconst_gen_jetprt")
			      // 	continue;
			      TH1 * hjoined = join(h_chi[pair_ind][charge_ind], h_chi[3][charge_ind]);
			      printf("type_ind %u level_ind %u\n", type_ind, level_ind);
			      // printf("h_chi[pair_ind][charge_ind]\n");
			      // h_chi[pair_ind][charge_ind] -> Print("all");
			      // printf("h_chi[3][charge_ind]\n");
			      // h_chi[3][charge_ind] -> Print("all");
			      // printf("hjoined\n");
			      // hjoined -> Print("all");
			      // getchar();
			      
			      TH1 * hjoined_syst = nullptr;
			      
			      if (type_ind == MC)
				{
				  hjoined_syst = join(h_chi_syst[pair_ind][charge_ind], h_chi_syst[3][charge_ind]);				   
				}
			      // for (unsigned char bind = 1; bind <= hjoined -> GetNbinsX(); bind ++)
			      // 	{
			      // 	  printf("bind %u \t %f \t%f \t syst cont %f \t syst error%f\n", bind, hjoined -> GetBinContent(bind), hjoined -> GetBinError(bind), hjoined_syst -> GetBinContent(bind), hjoined_syst -> GetBinError(bind));
			      // 	}
			      // getchar();
			      
			      printf("%s\n", dir.Data());
			      Rvalue[charge_ind][pair_ind] = calculate_result(hjoined, hjoined_syst, R, nullptr, 1.0, 1.0);
			      Rvalue[charge_ind][pair_ind] . SetName("R");
//			      printf("type %u, charge %u pair %u\n", type_ind, charge_ind, pair_ind);
			      Rvalue[charge_ind][pair_ind].ls();
			      //			      getchar();
			      // printf("type %u, pair %u charge %u ch %u level %u flow %u\n", type_ind, pair_ind, charge_ind, ch_ind, level_ind, flow_ind);
			      // getchar();
			      delete hjoined;
			      delete hjoined_syst;
			    }
		    }
		  for (unsigned short pair_ind = 0; pair_ind < npairs - 1; pair_ind ++)
		    {
		      fprintf(table_file, "\t\t\\multicolumn{2}{c}{\\boldmath\\bf %s}\\\\\n", jetpairstitle[pair_ind]);
		      fprintf(table_file, "\t\t\\hline\n");
		      for (unsigned short charge_ind = 0; charge_ind < N_charge_types_; charge_ind ++)
			{
			  fprintf(table_file, "\t\t %s ", charge_titles_[charge_ind]);
			  if (level_ind == GEN and type_ind == DATA)
			    {
			      fprintf(table_file, "\t& x \\\\\n");
			      continue;
			    }
			  add_entry(table_file, Rvalue[charge_ind][pair_ind]);

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
  //  printf("dir %s\n", dir);
  
  TH1F * h = (TH1F*) (plotter_file -> GetDirectory(dir) -> Get(TString(dir) + "_" +  MC_histo_names[0]) -> Clone(TString(dir) + "_MCsum"));
  //printf("hname %s h int %f\n", h-> GetName(), h -> Integral()); 
  for (unsigned short MC_ind = 1; MC_ind < N_MC_histo_names; MC_ind ++)
    {
      h -> Add((TH1F*)  (plotter_file -> GetDirectory(dir) -> Get(TString(dir) + "_" + MC_histo_names[MC_ind])));

    }
  return h;
	       
}

double R(TH1 * h, void * params)
{
  //inter W : intra W
  //h -> Print("all");
  
  TAxis * axis = h -> GetXaxis();
  const double A(h -> Integral(axis -> FindBin(1.2), axis -> FindBin(1.8)));
  const double B(h -> Integral(axis -> FindBin(0.0), axis -> FindBin(1.0) - 1));
  char error[256];
  const double criterion(1E-6);
  if (h -> Integral() == 0.0)
    {
      throw "Integral of input histogram is zero\n";
    }
  if (A < criterion * h -> Integral() and B >= criterion * h -> Integral())
    {
      sprintf(error, "Integral from 1.2 to 1.8 in the intra W region too low");
      throw error;
    }
  else if (A >= criterion * h -> Integral() and B < criterion * h -> Integral())
    {
      sprintf(error, "Integral from 0.2 to 0.8 in the inter W region too low");
      throw error;
    }
  else if (A < criterion * h -> Integral() and B < criterion * h -> Integral())
    {
      sprintf(error, "Integral from 0.2 to 0.8 in the inter W region and integral from 1.0 to 2.0 in the intra W region too low");
      throw error;
    }
 
  double ret = h -> Integral(axis -> FindBin(0.2), axis -> FindBin(0.8))/A * h -> Integral(axis -> FindBin(1.0), axis -> FindBin(2.0))/B;
  if (TMath::IsNaN(ret))
    {
      printf("A %.9f B %.9f integral %.9f\n", A, B, h-> Integral());
      getchar();
    }
  //printf ("ret %f\n", ret);
  return ret;
}

void prepare_header(FILE * file, const char * caption, const char * label, SourceCode_t source)
{
  fprintf(file, "\\begin{table}[h!btp]\n");
  fprintf(file, "\t\\centering\n");
  fprintf(file, "\t\\caption{%s}\n", caption);
  fprintf(file, "\t\\label{tab:%s}\n", label);
  fprintf(file, "\t\\begin{tabular}{l|c}\n");
  fprintf(file, "\t\t\\noalign{\\global\\arrayrulewidth=0.5mm}\\hline\\noalign{\\global\\arrayrulewidth=0.4pt}\n");
  char err[128];
  if (source == MC)
    {
      sprintf(err, "%s", "$\\pm$(stat)$\\pm$(syst)");
    }
  else
    {
      sprintf(err, "%s", "$\\pm$(stat)");
    }
  if (string(env).compare("lx") == 0)
    {
      fprintf(file, "\t\tJet constituents & $R$%s \\\\\n", err, err);
    }
  else
    {
      fprintf(file, "\t\t\\makecell[c]{\\bf Jet constituents} & {\\boldmath\\makecell{\\bf$I$%s [rad] }} & {\\boldmath\\makecell{\\bf$R^{-1}$%s}}\\\\\n", err, err);
    }
    fprintf(file, "\t\t\\hline\n");
  

}

void add_entry(FILE * table_file,  Result Rval)
{
  // printf("adding entry %s\n", Rval.latex(".3f").Data());
  fprintf(table_file, "\t& %s \\\\\n", Rval.latex(".3f").Data());
}

void prepare_footer(FILE *file)
{
  fprintf(file, "\t\\noalign{\\global\\arrayrulewidth=0.5mm}\\hline");
  fprintf(file, "\t\\end{tabular}\n");
  fprintf(file, "\\end{table}\n");

}

TH1 * join (TH1* hintra, TH1 * hinter)
  
{
  if (not hintra or not hinter) return nullptr;
  const unsigned char nbins ( hinter -> GetNbinsX());
  TH1F * hjoined = new TH1F (TString(hinter -> GetName()) + "_joined", 
			     TString(hinter -> GetName()) + "_joined", 
			     2 * hinter -> GetNbinsX(), 
			     hintra -> GetXaxis() -> GetBinLowEdge(1), 
			     2 *hintra -> GetXaxis() -> GetBinUpEdge(nbins));
  for (unsigned char bind = 1; bind <= nbins; bind ++)
    {
      hjoined -> SetBinContent(bind, hinter -> GetBinContent(bind));
      hjoined -> SetBinContent(bind + nbins, hintra -> GetBinContent(bind));
      hjoined -> SetBinError(bind, hinter -> GetBinError(bind));
      hjoined -> SetBinError(bind + nbins, hintra -> GetBinError(bind));
    }
  return hjoined;
}
