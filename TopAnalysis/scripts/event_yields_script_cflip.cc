#include "TFile.h"
#include "TApplication.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TMath.h"
const unsigned short N_bins = 4;
const unsigned short N_MC_histo_names = 6;
const char * MC_histo_table_label = "\\ttbar cflip";


void insert_header(FILE *);
void insert_channel_header(FILE *, const char *);
void insert_MC(FILE *, double[N_bins]);
void insert_MC_unc(FILE *, TH1 *);
void insert_footer(FILE *);
int main()
{
  FILE * file = fopen("event_yields_tables/event_yields_table_cflip.txt", "w");
  insert_header(file);
  TFile * plotter = TFile::Open("$EOS/analysis_MC13TeV_TTJets_cflip/plots/plotter.root");
  const unsigned short N_levels = 2;
  const char * level_names[N_levels] = {"reco", "gen"};
  const unsigned short N_ch = 3;
  const char * ch_names[N_ch] = {"E", "M", "L"};
  const unsigned short N_histo_names = 7;
  const char * histo_names[N_histo_names] = {"", "t#bar{t}", "Single top", "W", "DY", "Multiboson", "t#bar{t}+V"};
  const char * MC_histo_name = "t#bar{t} cflip";
  
  for (unsigned short level_ind  = 0; level_ind < N_levels - 1; level_ind ++)
    {
      for (unsigned short ch_ind = 0; ch_ind < N_ch; ch_ind ++)
	{
	  insert_channel_header(file, ch_names[ch_ind]);
	  double MC_results[N_bins];
	
	  for (unsigned short bin_ind = 0; bin_ind < N_bins; bin_ind ++)
	    {
		  MC_results[bin_ind] = 0.0;
	    }
	  TH1 * h_MC = NULL;
	  TH1 * h_MC_Unc = NULL;
	  const TString dir(TString(ch_names[ch_ind]) + "_" + level_names[level_ind] + "_selection");
	  printf("dir %s\n", dir.Data());
	  plotter -> GetDirectory(dir) -> ls();
     
	  TList * list = plotter -> GetDirectory(dir) -> GetListOfKeys();
	  TObject * it = list -> First();
	  printf("last %s\n", list -> Last() -> GetName());
	  do
	    {
	      TObject * obj = plotter  -> Get(dir + "/" + it -> GetName()); 
	      printf("[%s] %s\n", it -> GetName(), obj -> ClassName());

	      if (obj -> InheritsFrom("TH1"))
		{
		  TH1 * h = (TH1 *) obj;
		   if (h -> GetName() == TString("totalmcunc"))
		     {
		       printf("MC unc histo found\n");
		       h_MC_Unc = h;
		     }
		   if (h -> GetName() == dir + "_" + MC_histo_name)
		     {
		       h_MC = h;
		     }

		}
	      it = list -> After(it);
	    } while (it != list -> After(list -> Last()));
	  for (unsigned short bin_ind = 1; bin_ind <= h_MC -> GetNbinsX(); bin_ind ++)
	    {
	      printf("bin %u content %f check %f\n", bin_ind, h_MC -> GetBinContent(bin_ind), TMath::Sqrt(h_MC -> GetBinContent(bin_ind)));
	      if (h_MC)
		{
		  MC_results[bin_ind - 1] = h_MC -> GetBinContent(bin_ind);
		}
	    }
	  insert_MC(file, MC_results);
	  insert_MC_unc(file, h_MC_Unc);
	}

    }
  /*  TApplication app("myapp", 0, 0);
  app.Run(kTRUE);
  app.Terminate();*/
  plotter -> Close();
  insert_footer(file);
  fclose(file);
}

void insert_header(FILE * file)
{
  fprintf(file, "\\begin{table}[htp]\n");
  fprintf(file, "\\centering\n");
  fprintf(file, "\\caption{Event yields for the colour octet $W$ sample}\n");
  fprintf(file, "\\label{tab:yields}\n");
  fprintf(file, "\\begin{longtable}{lrrrr}\n");
  fprintf(file, "Process & $1 \\ell$ & $1 \\ell + \\geq 4 j$ & $1 \\ell + \\geq 4 j (2 b)$ & $1 \\ell + 4 j (2 b, 2 lj)$ \\\\\n");
  fprintf(file, "\\hline\n");

}

void insert_channel_header(FILE * file, const char * channel)
{
  fprintf(file, "%s + jets channel\\\\\n", channel);
  fprintf(file, "\\hline\n");
}

void insert_MC(FILE * file, double MC_results[N_bins])
{
  fprintf(file, "%*s\t", 25, MC_histo_table_label);
  for (unsigned short bin_ind = 0; bin_ind < N_bins; bin_ind ++)
    {
      fprintf(file, "& %20.1f\t\t", MC_results[bin_ind]);
    }
  fprintf(file, "\\\\\n");
}



void insert_MC_unc(FILE * file, TH1 * h_MC_Unc)
{
  fprintf(file, "%*s", 25, "(\\ttbar cflip\\ uncertainty)\t");
  for (unsigned short bin_ind = 0; bin_ind < N_bins; bin_ind ++)
    {
      fprintf(file, "& $\\pm$ %16.1f\t", h_MC_Unc -> GetBinError(bin_ind + 1));
    }
  fprintf(file, "\\\\\n");
  fprintf(file, "\\hline\n");


}

void insert_footer(FILE * file)
{
  fprintf(file, "\\end{longtable}\n");
  fprintf(file, "\\end{table}\n");
}
