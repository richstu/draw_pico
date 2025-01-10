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

//Currently the an_kinematicfit has the year as an input variable, but this can be removed
//This is only used in defining processes
//int main() {
int main(int argc, char *argv[]) {
  cout << "Plotting kinematic fit plots" << endl;

  bool plot_data = false; bool era_select = true; 
  string era_select_string  = (argc > 1) ? argv[1] : "";
  string data_select_string = (argc > 2) ? argv[2] : ""; 

  if(era_select_string=="run3"){  era_select = true; }
  if(data_select_string=="data"){ plot_data  = true; } 

  //Uses txt/zg_samples.txt to define needed processes
  vector<shared_ptr<Process>> procs = ZgUtilities::procs_with_sig_scale("txt/samples_zgamma.txt","AllMoreTrigsData",10);
  if(!plot_data && !era_select){
    cout << "Plotting Run 2 MC only." << endl;
    procs = ZgUtilities::procs_with_sig_scale("txt/samples_zgamma.txt","AllMoreTrigs",10);
  } else if(!plot_data && era_select){
    cout << "Plotting Run 3 MC only." << endl;
    procs = ZgUtilities::procs_with_sig_scale("txt/samples_zgamma.txt","AllMoreTrigsRun3",10);
  } else if(plot_data && era_select){
    cout << "Plotting Run 3 data and MC only." << endl;
    procs = ZgUtilities::procs_with_sig_scale("txt/samples_zgamma.txt","AllMoreTrigsDataRun3",10);
  } else {
    cout << "Plotting data and MC." << endl;
  }

  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");

  //Defines the plot options for 1D and 2D plots
  PlotOpt lin_lumi("txt/plot_styles.txt","CMSPaper");
  lin_lumi.Title(TitleType::info)
          .YAxis(YAxisType::linear)
          .Stack(StackType::data_norm)
          .Overflow(OverflowType::both)
          .YTitleOffset(1.)
          .AutoYAxis(true)
          .UseCMYK(false)
          .LeftMargin(0.17)
          .LegendColumns(1)
          .CanvasWidth(1077)
          .FileExtensions({"pdf","root"});

  PlotOpt style("txt/plot_styles.txt","Eff2D");
  vector<PlotOpt> ops_2D = {style().Stack(StackType::data_norm)
                                   .YTitleOffset(1.)
                                   .LogMinimum(0.1)
                                   .Overflow(OverflowType::both)
                                   .CanvasWidth(1077)
                                   .LabelSize(0.04)
                                   .UseCMYK(true)  
                                   .YAxis(YAxisType::log)
                                   .Title(TitleType::info)};

  vector<PlotOpt> ops = {lin_lumi};

  cout << "Adding samples." << endl;


  //This defines the splits based on lepton flavor that will be plotted
  vector<NamedFunc> lep = {"ll_lepid[0]==11","ll_lepid[0]==13","1"};
  vector<string> lep_lab= {"_ee", "_mumu", "_ll"};

  //This defines the categorization that will be used for plotting
  vector<NamedFunc> cat_nomll_vec = CatUtilities::run3_catwsel_nomll_vector;
  vector<NamedFunc> cat_wsel_vec  = CatUtilities::run3_catwsel_vector;
  vector<string> cat_vec_str      = CatUtilities::run3_category_labels;                          

  //These are all the selections used as a part of the baseline selection from the zg_functions.cpp file
  NamedFunc tight_baseline = ZgFunctions::tightened_baseline;
  NamedFunc tight_baseline_minus_mll = Nminusk(ZgFunctions::vector_tightened_baseline,{5});

  //This block of code loops through to create the NamedFuncs for each category
  //The purpose is to only require one loop while making plots
  vector<NamedFunc> NamedFunc_loop  = {};
  vector<NamedFunc> NamedFunc_loop_nomll  = {};
  vector<string>    string_loop_label = {};
  for(unsigned int idx_i = 0; idx_i < lep.size(); idx_i++){
    for(unsigned int idx_j = 0; idx_j < cat_wsel_vec.size(); idx_j++){

      //Vectors used for plots with the tightened baseline selection
      NamedFunc_loop_nomll.push_back( tight_baseline_minus_mll && cat_nomll_vec[idx_j] && lep[idx_i]);
      NamedFunc_loop.push_back(       tight_baseline           && cat_wsel_vec[idx_j]  && lep[idx_i]);
      string_loop_label.push_back(cat_vec_str[idx_j] + lep_lab[idx_i]);
    }
  }

  //This block of code handles all of the plot making
  PlotMaker pm;
  string pltNames = "";
  NamedFunc sel_cat    = "1";
  NamedFunc sel_cat_rf = "1";

  cout << "Adding the plots to the PlotMaker object." << endl;

  //Makes plots with the tightened baseline selection
  for(unsigned int idx_plt = 0; idx_plt < NamedFunc_loop.size(); idx_plt++){
    //These are not necessary but sometimes make it a bit more concise when making plots
    pltNames = string_loop_label[idx_plt];
    if(plot_data){pltNames = pltNames + "_DATA"; }

    NamedFunc not_hwin = "llphoton_m[0] < 120 || llphoton_m[0] > 130";
    sel_cat  = NamedFunc_loop_nomll[idx_plt];
    sel_cat_rf = sel_cat;
    if(plot_data){
      sel_cat_rf = sel_cat && not_hwin;
      sel_cat = sel_cat && not_hwin;
    }

    //No refit no mll selection
    pm.Push<Hist1D>(Axis(70, 50, 120,  "ll_m[0]",       "m_{ll} [GeV]", {}),       sel_cat, procs, ops).Weight(wgt).Tag("ShortName:an_kinematicfit_" + pltNames + "_ll_m_nomll");
    pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), sel_cat, procs, ops).Weight(wgt).Tag("ShortName:an_kinematicfit_" + pltNames + "_llphoton_m_nomll");
   
    pm.Push<Hist2D>(
      Axis(70,50,120,  "ll_m[0]", "m_{ll} [GeV]", {}),
      Axis(40, 100, 180,  "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
      sel_cat, procs, ops_2D).Tag("ShortName:an_kinematicfit_" + pltNames + "_mll_mlly");

    //Refit no mll selection
    pm.Push<Hist1D>(Axis(70, 50, 120,  "ll_refit_m",       "m_{ll,refit} [GeV]", {}),       sel_cat_rf, procs, ops).Weight(wgt).Tag("ShortName:an_kinematicfit_" + pltNames + "_refit_ll_m_nomll");
    pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_refit_m", "m_{ll#gamma,refit} [GeV]", {}), sel_cat_rf, procs, ops).Weight(wgt).Tag("ShortName:an_kinematicfit_" + pltNames + "_refit_llphoton_m_nomll");
    sample_kinrefit_plots(pm, sel_cat_rf, procs, ops, wgt, pltNames + "_nomll");

    pm.Push<Hist2D>(
      Axis(70,50,120,  "ll_refit_m", "m_{ll,refit} [GeV]", {}),
      Axis(40, 100, 180,  "llphoton_refit_m", "m_{ll#gamma,refit} [GeV]", {}),
      sel_cat_rf, procs, ops_2D).Tag("ShortName:an_kinematicfit_" + pltNames + "_mll_mlly_refit");


    //Adding mll selection sel_cat
    sel_cat = NamedFunc_loop[idx_plt];
    sel_cat_rf = sel_cat;
    if(plot_data){
      sel_cat_rf = sel_cat && "llphoton_refit_m < 122 || llphoton_refit_m > 128";
      sel_cat = sel_cat && "llphoton_m[0] < 122 || llphoton_m[0] > 128";
    }

    //No refit with mll selection
    pm.Push<Hist1D>(Axis(70, 50, 120,  "ll_m[0]",       "m_{ll} [GeV]", {}),       sel_cat, procs, ops).Weight(wgt).Tag("ShortName:an_kinematicfit_" + pltNames + "_ll_m");
    pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), sel_cat, procs, ops).Weight(wgt).Tag("ShortName:an_kinematicfit_" + pltNames + "_llphoton_m");
    sample_kinrefit_plots(pm, sel_cat_rf, procs, ops, wgt, pltNames);

  }

  //This block controls the verbose or non-verbose plotting option and whether or not to use multiple threads
  pm.min_print_ = true;
  pm.multithreaded_ = true;

  //Luminosity option for the plots
  pm.SetLuminosityTag("137.61").MakePlots(1);
}
