#include "TFile.h"
#include "TApplication.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TMath.h"
const unsigned short N_bins = 4;
const unsigned short N_MC_histo_names = 7;
const char * MC_histo_table_labels[N_MC_histo_names] = {"\\ttbar", "Single top", "W", "DY", "Multiboson", "\\ttbar+V", "QCD"};


void insert_header(FILE *);
void insert_channel_header(FILE *, const char *);
void insert_MC(FILE *, double (*)[N_bins], double *);
void insert_MC_sum(FILE *, double *);
void insert_MC_unc(FILE *, TH1 *);
void insert_data(FILE *, TH1 *);
void insert_footer(FILE *);
int main()
{
  FILE * file = fopen("event_yields_tables/event_yields_table.txt", "w");
  insert_header(file);
  TFile * plotter = TFile::Open("$EOS/analysis_MC13TeV_TTJets/plots/plotter.root");
  TFile * plotter_cflip = TFile::Open("$EOS/analysis_MC13TeV_TTJets_cflip/plots/plotter.root");
  const unsigned short N_levels = 2;
  const char * level_names[N_levels] = {"reco", "gen"};
  const unsigned short N_ch = 3;
  const char * ch_names[N_ch] = {"E", "M", "L"};
  const unsigned short N_histo_names = 8;
  const char * histo_names[N_histo_names] = {"", "t#bar{t}", "Single top", "W", "DY", "Multiboson", "t#bar{t}+V", "QCD"};
  const char * MC_histo_names[N_MC_histo_names] = {"t#bar{t}", "Single top", "W", "DY", "Multiboson", "t#bar{t}+V", "QCD"};
  
  for (unsigned short level_ind  = 0; level_ind < N_levels - 1; level_ind ++)
    {
      for (unsigned short ch_ind = 0; ch_ind < N_ch; ch_ind ++)
	{
	  insert_channel_header(file, ch_names[ch_ind]);
	  double MC_sum[N_bins];
	  double MC_results[N_MC_histo_names][N_bins];
	  double MC_results_cflip[N_bins];
	  for (unsigned short bin_ind = 0; bin_ind < N_bins; bin_ind ++)
	    {
	      MC_sum[bin_ind] = 0.0;
	      MC_results_cflip[bin_ind] = 0.0;
	      for (unsigned short h_ind = 0; h_ind < N_MC_histo_names; h_ind ++)
		{
	      
		  MC_results[h_ind][bin_ind] = 0.0;
		}
	    }
	  TH1 * h_data = NULL;
	  TH1 * h_MC_Unc = NULL;
	  const TString dir(TString(ch_names[ch_ind]) + "_" + level_names[level_ind] + "_selection");
	  printf("dir %s\n", dir.Data());
	  TList * list = plotter -> GetDirectory(dir) -> GetListOfKeys();
	  TObject * it = list -> First();
	  do
	    {
	      TObject * obj = plotter  -> Get(dir + "/" + it -> GetName()); 
	      printf("[%s] %s\n", it -> GetName(), obj -> ClassName());

	      if (obj -> InheritsFrom("TH1"))
		{
		  TH1 * h = (TH1 *) obj;
		  bool isMC = false;
		   if (h -> GetName() == dir)
		     {
		       printf("data histo found \n");
		       h_data = h;
		     }
		   if (h -> GetName() == TString("totalmcunc"))
		     {
		       printf("MC unc histo found\n");
		       h_MC_Unc = h;
		     }
		  unsigned short MC_sample_ind = 0;
		  unsigned short order = 0;
		  do
		    {
		      const TString comp = dir + "_" + MC_histo_names[MC_sample_ind];
		      if (comp == h -> GetName())
			{
			  isMC = true;
			  
			  order = MC_sample_ind;
			  printf("is MC %u\n", order);
			}
		      MC_sample_ind ++;
		    } while (MC_sample_ind < N_MC_histo_names); 
		  printf("\n Listing %s \n", h -> GetName());
		  for (unsigned short bin_ind = 1; bin_ind <= h -> GetNbinsX(); bin_ind ++)
		    {
		      printf("bin %u content %f check %f\n", bin_ind, h -> GetBinContent(bin_ind), TMath::Sqrt(h -> GetBinContent(bin_ind)));
		      if (isMC)
			{
			  MC_sum[bin_ind - 1] += h -> GetBinContent(bin_ind);
			  MC_results[order][bin_ind - 1] = h -> GetBinContent(bin_ind);
			}
		    }
		  printf("\n errors \n");
		  for (unsigned short bin_ind = 0; bin_ind <= h -> GetNbinsX() + 1; bin_ind ++)
		    {
		      printf("bin %u content %f\n", bin_ind, h -> GetBinError(bin_ind));
		    }

		  printf("\nhisto listed \n");
		}
	      it = list -> After(it);
	    } while (it != list -> After(list -> Last()));
	  TH1F * h_cflip = (TH1F*) plotter_cflip -> GetDirectory(dir) -> Get(dir + "_t#bar{t} cflip");
	  for (unsigned short bin_ind = 1; bin_ind < h_cflip -> GetNbinsX() + 1; bin_ind ++)
	    {
	      MC_results_cflip[bin_ind - 1] = h_cflip -> GetBinContent(bin_ind);
	    }

	  insert_MC(file, MC_results, MC_results_cflip);
	  insert_MC_sum(file, MC_sum);
	  insert_MC_unc(file, h_MC_Unc);
	  insert_data(file, h_data);
	}

    }
  /*  TApplication app("myapp", 0, 0);
  app.Run(kTRUE);
  app.Terminate();*/
  plotter -> Close();
  plotter_cflip -> Close();
  insert_footer(file);
  fclose(file);
}

void insert_header(FILE * file)
{
  fprintf(file, "\\begin{table}[htp]\n");
  fprintf(file, "\\centering\n");
  fprintf(file, "\\caption{Event yields}\n");
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

void insert_MC(FILE * file, double MC_results[][N_bins], double * MC_results_cflip)
{
  for (unsigned short MC_samples_ind = 0; MC_samples_ind < N_MC_histo_names; MC_samples_ind ++)
    {
      fprintf(file, "%*s\t", 25, MC_histo_table_labels[MC_samples_ind]);
       for (unsigned short bin_ind = 0; bin_ind < N_bins; bin_ind ++)
	 {
	   fprintf(file, "& %20.1f\t\t", MC_results[MC_samples_ind][bin_ind]);
	 }
       fprintf(file, "\\\\\n");			
       /*       if (MC_samples_ind == 0)
	 {
	   fprintf(file, "\\small %*s\t", 25, "t\\overline{t} cflip");
	   for (unsigned short bin_ind = 0; bin_ind < N_bins; bin_ind ++)
	     {
	       fprintf(file, "& \\small %20.1f\t\t", MC_results_cflip[bin_ind]);
	     }
	   fprintf(file, "\\\\\n");	
	   }*/
    }
  fprintf(file, "\\hline\n");
}


void insert_MC_sum(FILE * file, double * MC_sum)
{
  fprintf(file, "%*s", 25, "Total MC\t");
  for (unsigned short bin_ind  = 0; bin_ind < N_bins; bin_ind ++)
    {
      fprintf(file, "& %20.1f\t\t", MC_sum[bin_ind]);

    }
  fprintf(file, "\\\\\n");
}

void insert_MC_unc(FILE * file, TH1 * h_MC_Unc)
{
  fprintf(file, "%*s", 25, "(\\ttbar\\ uncertainty)\t");
  for (unsigned short bin_ind = 0; bin_ind < N_bins; bin_ind ++)
    {
      fprintf(file, "& $\\pm$ %16.1f\t", h_MC_Unc -> GetBinError(bin_ind + 1));
    }
  fprintf(file, "\\\\\n");
}

void insert_data(FILE * file, TH1 * h_data)
{
  fprintf(file, "%*s", 25, "Data\t");
  for (unsigned short bin_ind = 0; bin_ind < N_bins; bin_ind ++)
    {
      fprintf(file, "& %20.1f\t\t", h_data -> GetBinContent(bin_ind + 1));
    }
  fprintf(file, "\\\\\n\\hline\n");
  
}

void insert_footer(FILE * file)
{
  fprintf(file, "\\end{longtable}\n");
  fprintf(file, "\\end{table}\n");
}
