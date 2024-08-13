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
#include "zgamma/zg_utilities.hpp"
#include "zgamma/zg_functions.hpp"
#include "zgamma/categorization_utilities.hpp"
using namespace std;
using namespace PlotOptTypes;
using namespace ZgUtilities;
using namespace ZgFunctions;
using namespace CatUtilities;

//These NamedFuncs are to define the different photon eta regions used in below
NamedFunc photon_barrel("photon_barrel", [](const Baby &b) -> NamedFunc::ScalarType{
  return fabs(b.photon_eta()-> at(0)) < 1.3;
});

NamedFunc photon_near_crack("photon_near crack", [](const Baby &b) -> NamedFunc::ScalarType{
  return fabs(b.photon_eta()-> at(0)) >  1.3 && fabs(b.photon_eta()-> at(0)) < 1.7;
});

NamedFunc photon_endcap("photon_endcap", [](const Baby &b) -> NamedFunc::ScalarType{
  return fabs(b.photon_eta()-> at(0)) >  1.7;
});


//Currently the skeleton has the year as an input variable, but this can be removed
//This is only used in defining processes
//int main() {
int main() {

  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");


  //Defines the plot options for 1D and 2D plots
  PlotOpt style("txt/plot_styles.txt","Eff2D");
  vector<PlotOpt> bkg_hist = {style().Stack(StackType::data_norm)
                                     .YTitleOffset(1.)
                                     .Overflow(OverflowType::both)
                                     .LabelSize(0.04)
                                     .UseCMYK(true)  
                                     .YAxis(YAxisType::log)
                                     .Title(TitleType::info)};
  vector<PlotOpt> ops_2D = {bkg_hist};

  //Output root files will be used to make overlay plots
  PlotOpt lin_lumi("txt/plot_styles.txt","CMSPaper");
  lin_lumi.YAxis(YAxisType::linear)
          .Stack(StackType::signal_overlay)
          .YTitleOffset(1.75)
          .AutoYAxis(false)
          .UseCMYK(false)
          .LeftMargin(0.17)
          .CanvasWidth(800)
          .CanvasHeight(800)
          .Title(TitleType::info)
          .FileExtensions({"pdf","root"});
  vector<PlotOpt> ops = {lin_lumi()};

  //Uses txt/zg_samples.txt to define needed processes
  vector<shared_ptr<Process>> procs = ZgUtilities::procs_with_sig_scale("txt/samples_zgamma.txt","AllMoreTrigs",10);

  //These are all the selections used as a part of the baseline selection from the zg_functions.cpp file
  NamedFunc tight_baseline = ZgFunctions::tightened_baseline;
  NamedFunc mll_selection  = "ll_m[0] > 80 && ll_m[0] < 100"; 

  //This defines the categorization that will be used for plotting
  vector<NamedFunc> cat_wsel_vec = CatUtilities::run3_catwsel_vector;
  vector<string>    cat_vec_str  = CatUtilities::run3_category_labels;                          
  
  //This defines the splits based on lepton flavor that will be plotted
  vector<NamedFunc> lep     = { "ll_lepid[0]==11","ll_lepid[0]==13","1"};
  vector<string>    lep_lab = {"_ee", "_mumu", "_ll"};

  //This defines the splits based on lepton flavor that will be plotted
  vector<NamedFunc> gamma_eta_vec = {photon_barrel, photon_near_crack, photon_endcap};
  vector<string>    gamma_eta_lab = {"_barrel", "_crack", "_endcap"};

  //This block of code loops through to create the NamedFuncs for each category
  //The purpose is to only require one loop while making plots
  vector<NamedFunc> NamedFunc_loop  = {};
  vector<NamedFunc> NamedFunc_loop_nomll  = {};
  vector<string>    string_loop_label = {};
  for(unsigned int idx_lep = 0; idx_lep < lep.size(); idx_lep++){
    for(unsigned int idx_ph = 0; idx_ph < lep.size(); idx_ph++){
      for(unsigned int idx_cat = 0; idx_cat < cat_wsel_vec.size(); idx_cat++){

        //Vectors used for plots with the tightened baseline selection
        NamedFunc_loop.push_back(tight_baseline && cat_wsel_vec[idx_cat]  && lep[idx_lep]);
        string_loop_label.push_back(cat_vec_str[idx_cat] + lep_lab[idx_lep]);
      }
    }
  }

  //This block of code handles all of the plot making
  PlotMaker pm;
  string    plt_lab = "";
  NamedFunc sel_cat = "1";

  //Makes plots with the tightened baseline selection
  for(unsigned int idx_plt = 0; idx_plt < NamedFunc_loop.size(); idx_plt++){
    //These are not necessary but sometimes make it a bit more concise when making plots
    plt_lab = string_loop_label[idx_plt];
    sel_cat  = NamedFunc_loop[idx_plt];

    //Refit no mll selection
    pm.Push<Hist1D>(Axis(70, 50, 120,  "ll_refit_m[0]",       "m_{ll,refit} [GeV]", {}),       sel_cat, procs, ops).Weight(wgt).Tag("ShortName:Skeleton_" + plt_lab + "_refit_ll_m_nomll");
    pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_refit_m[0]", "m_{ll#gamma,refit} [GeV]", {}), sel_cat, procs, ops).Weight(wgt).Tag("ShortName:Skeleton_" + plt_lab + "_refit_llphoton_m_nomll");

    pm.Push<Hist2D>(
      Axis(70,  50, 120,  "ll_refit_m[0]",       "m_{ll,refit} [GeV]", {}),
      Axis(40, 100, 180,  "llphoton_refit_m[0]", "m_{ll#gamma,refit} [GeV]", {}),
      sel_cat, procs, ops_2D).Weight(wgt).Tag("ShortName:Skeleton_" + plt_lab+"_mll_mlly_refit");

    //Refit w/ mll selection
    sel_cat  = NamedFunc_loop[idx_plt] && mll_selection;
    pm.Push<Hist1D>(Axis(70, 50, 120,  "ll_refit_m[0]",       "m_{ll,refit} [GeV]", {}),       sel_cat, procs, ops).Weight(wgt).Tag("ShortName:Skeleton_" + plt_lab + "_refit_ll_m_nomll");
    pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_refit_m[0]", "m_{ll#gamma,refit} [GeV]", {}), sel_cat, procs, ops).Weight(wgt).Tag("ShortName:Skeleton_" + plt_lab + "_refit_llphoton_m_nomll");

    pm.Push<Hist2D>(
      Axis(70,  50, 120,  "ll_refit_m[0]",       "m_{ll,refit} [GeV]", {}),
      Axis(40, 100, 180,  "llphoton_refit_m[0]", "m_{ll#gamma,refit} [GeV]", {}),
      sel_cat, procs, ops_2D).Weight(wgt).Tag("ShortName:Skeleton_" + plt_lab+"_mll_mlly_refit");

  }

  //This block controls the verbose or non-verbose plotting option and whether or not to use multiple threads
  pm.min_print_ = true;
  pm.multithreaded_ = true;

  //Luminosity option for the plots
  pm.SetLuminosityTag("200").MakePlots(1.0);
}
