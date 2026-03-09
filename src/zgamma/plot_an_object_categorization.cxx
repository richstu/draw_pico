#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <bitset>
#include <unistd.h>
#include <getopt.h>
#include "TError.h"
#include "TColor.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/utilities.hpp"
#include "core/named_func_utilities.hpp"
#include "zgamma/zg_utilities.hpp"
#include "zgamma/zg_functions.hpp"
#include "zgamma/categorization_utilities.hpp"
#include "zgamma/controlregion_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace ZgUtilities;
using namespace ZgFunctions;
using namespace CatUtilities;
using namespace CRUtilities;
using namespace NamedFuncUtilities;

std::unordered_map<std::string,double>  year_lumi{{"2016",19.5},{"2016APV",16.8},{"2017",41.48},{"2018",59.83}};

int ZHtype =200400; 
bool ztonunu_b(const Baby &b){  
  int N_momZ = 0;
  int N_Zdecayproducts = 0;
  for(unsigned int idx_mc = 0; idx_mc < b.mc_id() -> size(); idx_mc++){
    if(b.mc_id() -> at(idx_mc) == 23){N_momZ++; }
    if(b.mc_mom() -> at(idx_mc) == 23){N_Zdecayproducts++; }
  }
  
  //bool isZnunu = N_momZ > N_Zdecayproducts-2; 
  //cout << "Event has Z to nunu(" << b.SampleType() << "): " << isZnunu << endl;
  return N_momZ > N_Zdecayproducts-2;//In ZH one Z always has Z-->ll the -2 accounts for that
}

NamedFunc Znunu("Znunu", [](const Baby &b) -> NamedFunc::ScalarType{return (ztonunu_b(b) && b.type()==ZHtype) || (b.type()!=ZHtype);});
NamedFunc only_ztonunu_signal("only_ztonunu_signal",[](const Baby &b) -> NamedFunc::ScalarType{return b.type()==ZHtype ? ztonunu_b(b) : true;});


//Currently the skeleton has the year as an input variable, but this can be removed
//This is only used in defining processes
//int main() {
int main() {
  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");

  //Uses txt/zg_samples.txt to define needed processes
  
  //Run2 + Run3 
  vector<shared_ptr<Process>> procs = ZgUtilities::procs_with_sig_scale(CRUtilities::SamplesFile, "AllMoreTrigsSignalSplit", 10); 
  string year_label = "_all";
  string era_label  = CRUtilities::R23Lumi;
 
  //These are all the selections used as a part of the baseline selection from the zg_functions.cpp file
  NamedFunc tight_baseline = ZgFunctions::tightened_baseline_refit;

  //Defines the plot options for 1D and 2D plots
  vector<PlotOpt> ops    = CRUtilities::an_plotting_options("soverb"); ops[0].Title(PlotOptTypes::TitleType::simulation);
  vector<PlotOpt> ops_1D = CRUtilities::an_plotting_options("1D"); ops_1D[0].Title(PlotOptTypes::TitleType::simulation);

  //This block of code handles all of the plot making
  PlotMaker pm;

  //Functions to shorten the command to add plots to PlotMaker
  string label = year_label;
  NamedFunc selection = tight_baseline && higgs_window_refit;

  //Add plots to plot maker
  pm.Push<Hist1D>(Axis(5, 2, 7, "nlep",  "N_{l}",     {}), selection, procs, ops).Weight(wgt10).Tag("ShortName:an_ob_categorization_nl");
  pm.Push<Hist1D>(Axis(3, 0, 3, "nbdfm", "N_{b,med}", {}), selection, procs, ops).Weight(wgt10).Tag("ShortName:an_ob_categorization_nbdfm");
  pm.Push<Hist1D>(Axis(8, 0, 8, "njet",  "N_{jet}",   {}), selection, procs, ops).Weight(wgt10).Tag("ShortName:an_ob_categorization_njet");
  pm.Push<Hist1D>(Axis(75, 0, 150, "met", "p_{T}^{miss} [GeV]", {90}), selection && only_ztonunu_signal, procs, ops).Weight(wgt10).Tag("ShortName:an_ob_categorization_ptmiss");

  //ggF -- Convener comment
  NamedFunc sel_noptgmllg = Nminusk(ZgFunctions::vector_tightened_baseline,{4}) && higgs_window_refit;
  vector<shared_ptr<Process>> procs_ptg = ZgUtilities::procs_with_sig_scale(CRUtilities::SamplesFile, "AllMCRun2And3", 10); 
  pm.Push<Hist1D>(Axis(50, 0.05, 0.3, "photon_pt[0]/llphoton_m[0]", "p_{T}(#gamma)/m_{ll#gamma}", {15.0/110.0, 15.0/100.0}), sel_noptgmllg, procs_ptg, ops).Weight(wgt10).Tag("ShortName:an_ob_categorization_ptg_mllgamma");
  pm.Push<Hist1D>(Axis(50, 0.05, 0.3, "photon_pt[0]/llphoton_m[0]", "p_{T}(#gamma)/m_{ll#gamma}", {15.0/110.0, 15.0/100.0}), selection, procs_ptg, ops).Weight(wgt10).Tag("ShortName:an_ob_categorization_ptg_mllgamma_wcurrsel");




  //VBF
  pm.Push<Hist1D>(Axis(8, 0, 8, "njet", "N_{jet}",   {2}), selection && "nbdfm==0", procs, ops).Weight(wgt10).Tag("ShortName:an_ob_categorization_vbf_njet");

  //ttH leptonic plots
  pm.Push<Hist1D>(Axis(5, 0, 5, "nbdfm", "N_{b,med}", {1}), selection && "nlep==3",             procs, ops).Weight(wgt10).Tag("ShortName:an_ob_categorization_tthlep_nbdfm_nlep_e3");
  pm.Push<Hist1D>(Axis(8, 0, 8, "njet",  "N_{jet}",   {3}), selection && "nlep==3 && nbdfm>=1", procs, ops).Weight(wgt10).Tag("ShortName:an_ob_categorization_tthlep_njet_nlep_e3_nbdfm_ge1");

  pm.Push<Hist1D>(Axis(5, 0, 5, "nbdfm", "N_{b,med}", {1}), selection && "nlep>=4",             procs, ops).Weight(wgt10).Tag("ShortName:an_ob_categorization_tthlep_nbdfm_nlep_ge4");
  pm.Push<Hist1D>(Axis(8, 0, 8, "njet",  "N_{jet}",   {}),  selection && "nlep>=4 && nbdfm>=1", procs, ops).Weight(wgt10).Tag("ShortName:an_ob_categorization_tthlep_njet_nlep_ge4_nbdfm_ge1");

  //ttH hadronic plots
  pm.Push<Hist1D>(Axis(5, 0, 5, "nbdfm", "N_{b,med}", {1}), selection && "nlep==2", procs, ops).Weight(wgt10).Tag("ShortName:an_ob_categorization_tthhad_nbdfm_nlep_e2");
  pm.Push<Hist1D>(Axis(8, 0, 8, "njet",  "N_{jet}",   {5}), selection && "nlep==2 && nbdfm>=1", procs, ops).Weight(wgt10).Tag("ShortName:an_ob_categorization_tthhad_njet_nlep_e2_nbdfm_ge1");
  pm.Push<Hist1D>(Axis(8, 0, 8, "njet",  "N_{jet}",   {5}), selection && "nlep==2 && nbdfm>=2", procs, ops).Weight(wgt10).Tag("ShortName:an_ob_categorization_tthhad_njet_nlep_e2_nbdfm_ge2");


  //VH Nlep >= 3
  pm.Push<Hist1D>(Axis(5, 2, 7, "nlep",  "N_{l}",     {3}), selection && "nbdfm==0", procs, ops).Weight(wgt10).Tag("ShortName:an_ob_categorization_vh3l_nl_nbdfm");

  //ZH ptmiss
  pm.Push<Hist1D>(Axis(75, 0, 150, "met", "p_{T}^{miss} [GeV]", {90}), selection && "njet<=1" && Znunu, procs, ops).Weight(wgt10).Tag("ShortName:an_ob_categorization_zhptmiss_ptmiss");
 
  //The full set of categories
  vector<NamedFunc> categories_object_selections = CatUtilities::run3_category_vector;
  vector<string>    categories_object_labels     = CatUtilities::run3_category_labels;

  //Add a m_llphoton plot for each category with object selections
  NamedFunc category       = "1";
  string    category_label = "";
  vector<int> mlly_binning = {20,20,20,25,40,80,25};
  selection = tight_baseline; 
  for(unsigned int idx = 0; idx < categories_object_selections.size(); idx++){
    category       = categories_object_selections[idx];
    category_label =  categories_object_labels[idx];
    int nbins = mlly_binning[idx];
    pm.Push<Hist1D>(Axis(nbins, 100, 180, "llphoton_refit_m", "m_{ll#gamma} [GeV]", {}), selection && category, procs, ops_1D).Weight(wgt10).Tag("ShortName:an_ob_categorization_" + category_label);
  }

  //This block controls the verbose or non-verbose plotting option and whether or not to use multiple threads
  pm.min_print_ = true;
  pm.multithreaded_ = true;

  //Luminosity option for the plots
  pm.SetLuminosityTag(era_label).MakePlots(1); 
}


