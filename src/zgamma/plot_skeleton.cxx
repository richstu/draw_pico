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

//--This is a skeleton file meant to be used to easily make plotting files

  std::unordered_map<std::string,double>  year_lumi{{"2016",19.5},{"2016APV",16.8},{"2017",41.48},{"2018",59.83}};


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
          .FileExtensions({"pdf"});
  vector<PlotOpt> ops = {lin_lumi()};

  //Uses txt/zg_samples.txt to define needed processes
  vector<shared_ptr<Process>> procs = ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_zgamma.txt","All");

  //This defines the splits based on lepton flavor that will be plotted
  vector<NamedFunc> lep = { "ll_lepid[0]==11","ll_lepid[0]==13","1"};
  vector<string> lep_lab= {"_ee", "_mumu", "_ll"};

  //This defines the categorization that will be used for plotting
  vector<NamedFunc> cat_vec = CatUtilities::run3_category_vector;
  vector<string> cat_vec_str= CatUtilities::run3_category_labels;                          

  //These are all the selections used as a part of the baseline selection from the zg_functions.cpp file
  NamedFunc r2_baseline    = ZgFunctions::run2_baseline;
  NamedFunc tight_baseline = ZgFunctions::tightened_baseline;

  //This block of code loops through to create the NamedFuncs for each category
  //The purpose is to only require one loop while making plots
  vector<NamedFunc> NamedFunc_loop_run2  = {};
  vector<string>    string_loop_label_run2 = {};
  vector<NamedFunc> NamedFunc_loop_tight  = {};
  vector<string>    string_loop_label_tight = {};
  for(unsigned int idx_i = 0; idx_i < lep.size(); idx_i++){
    for(unsigned int idx_j = 0; idx_j < cat_vec.size(); idx_j++){

      //Vectors used for plots with Run 2 baseline selection
      NamedFunc_loop_run2.push_back( run2_baseline && cat_vec[idx_j] && lep[idx_i]);
      string_loop_label_run2.push_back(cat_vec_str[idx_j] + lep_lab[idx_i] + "_RUN2_BASELINE");

      //Vectors used for plots with the tightened baseline selection
      NamedFunc_loop_tight.push_back( run2_baseline && cat_vec[idx_j] && lep[idx_i]);
      string_loop_label_tight.push_back(cat_vec_str[idx_j] + lep_lab[idx_i] + "_TIGHT_BASELINE");
    }
  }

  //This block of code handles all of the plot making
  PlotMaker pm;
  string pltNames = "";
  NamedFunc sel_cat = "1";


  //Makes plots with the Run 2 baseline selection
  for(unsigned int idx_plt = 0; idx_plt < NamedFunc_loop_run2.size(); idx_plt++){
    //These are not necessary but sometimes make it a bit more concise when making plots
    pltNames = string_loop_label_run2[idx_plt];
    sel_cat  = NamedFunc_loop_run2[idx_plt];

    pm.Push<Hist1D>(Axis(70, 50, 120,  "ll_m[0]",       "m_{ll} [GeV]", {}),       sel_cat, procs, ops).Weight(wgt_1invfb).Tag("ShortName:Skeleton_" + pltNames+"_ll_m");
    pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), sel_cat, procs, ops).Weight(wgt_1invfb).Tag("ShortName:Skeleton_" + pltNames+"_llphoton_m");

    pm.Push<Hist2D>(
      Axis(70,50,120,  "ll_m[0]", "m_{ll} [GeV]", {}),
      Axis(40, 100, 180,  "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
      sel_cat, procs, ops_2D).Tag("ShortName:Skeleton_" + pltNames+"_mll_mlly");
  }

  //Makes plots with the tightened baseline selection
  for(unsigned int idx_plt = 0; idx_plt < NamedFunc_loop_tight.size(); idx_plt++){
    //These are not necessary but sometimes make it a bit more concise when making plots
    pltNames = string_loop_label_tight[idx_plt];
    sel_cat  = NamedFunc_loop_tight[idx_plt];

    pm.Push<Hist1D>(Axis(70, 50, 120,  "ll_m[0]",       "m_{ll} [GeV]", {}),       sel_cat, procs, ops).Weight(wgt_1invfb).Tag("ShortName:Skeleton_" + pltNames+"_ll_m");
    pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), sel_cat, procs, ops).Weight(wgt_1invfb).Tag("ShortName:Skeleton_" + pltNames+"_llphoton_m");

    pm.Push<Hist2D>(
      Axis(70,50,120,  "ll_m[0]", "m_{ll} [GeV]", {}),
      Axis(40, 100, 180,  "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
      sel_cat, procs, ops_2D).Tag("ShortName:Skeleton_" + pltNames+"_mll_mlly");
  }


  //This block controls the verbose or non-verbose plotting option and whether or not to use multiple threads
  pm.min_print_ = true;
  pm.multithreaded_ = true;

  //Luminosity option for the plots
  pm.MakePlots(137.61);

}
