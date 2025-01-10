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

int main(int argc, char *argv[]) {
  bool era_select = select_era(argc,argv);

  //Uses txt/zg_samples.txt to define needed processes
  vector<shared_ptr<Process>> procs     = ZgUtilities::procs_with_sig_scale("txt/samples_zgamma.txt", "AllMoreTrigs", 10);
  vector<shared_ptr<Process>> procs_sig = ZgUtilities::procs_with_sig_scale("txt/samples_zgamma.txt", "SignalSplit",  10);
  string lumi_label = "137.61";
  if(era_select){
    procs      = ZgUtilities::procs_with_sig_scale("txt/samples_zgamma.txt", "AllMoreTrigsRun3", 10);
    procs_sig  = ZgUtilities::procs_with_sig_scale("txt/samples_zgamma.txt", "SignalSplitRun3",  10);
    lumi_label = "62.32";
  }


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
          //.Bottom(BottomType::sorb)//TEST
          .CanvasWidth(800)
          .CanvasHeight(800)
          .Title(TitleType::info)
          .FileExtensions({"pdf"});
  vector<PlotOpt> ops = {lin_lumi()};


  //This defines the splits based on lepton flavor that will be plotted
  vector<NamedFunc> lep = {"1"};   // "ll_lepid[0]==11","ll_lepid[0]==13","1"};
  vector<string> lep_lab= {"_ll"}; //"_ee", "_mumu", "_ll"};

  //This defines the categorization that will be used for plotting
  vector<NamedFunc> cat_vec      = CatUtilities::run3_category_vector;
  vector<NamedFunc> cat_wsel_vec = CatUtilities::run3_catwsel_vector;
  vector<string>    cat_vec_str  = CatUtilities::run3_category_labels;                          

  //These are all the selections used as a part of the baseline selection from the zg_functions.cpp file
  NamedFunc tight_baseline = ZgFunctions::tightened_baseline;

  //This block of code loops through to create the NamedFuncs for each category
  //The purpose is to only require one loop while making plots
  vector<NamedFunc> NamedFunc_loop      = {};
  vector<NamedFunc> NamedFunc_refit_loop = {};
  vector<NamedFunc> NamedFunc_wsel_loop = {};
  vector<NamedFunc> NamedFunc_wsel_refit_loop = {};
  vector<string>    string_loop_label   = {};
  for(unsigned int idx_i = 0; idx_i < lep.size(); idx_i++){
    for(unsigned int idx_j = 0; idx_j < cat_vec.size(); idx_j++){

      //Vectors used for plots with the tightened baseline selection
      NamedFunc_loop.push_back( tightened_baseline && cat_vec[idx_j] && lep[idx_i]);
      NamedFunc_wsel_loop.push_back( tightened_baseline && cat_wsel_vec[idx_j] && lep[idx_i]);
      
      //Vectors used for plotting
      NamedFunc_refit_loop.push_back( tightened_baseline_refit && cat_vec[idx_j] && lep[idx_i]);
      NamedFunc_wsel_refit_loop.push_back( tightened_baseline_refit && cat_wsel_vec[idx_j] && lep[idx_i]);
      string_loop_label.push_back(cat_vec_str[idx_j] + lep_lab[idx_i]);
    }
  }

  //This block of code handles all of the plot making
  PlotMaker pm;
  string plot_lab = "";
  NamedFunc sel_cat = "1";
  //int Ncategories = 6;
  vector<TableRow> selection_tablerows = {};
  vector<vector<TableRow>> category_tablerows = {{},{},{},{},{},{}};

  vector<int> bins = {20,20,25,25,40,80};
  //Makes plots with the tightened baseline selection
  
  
  for(unsigned int idx_plt = 0; idx_plt < NamedFunc_loop.size(); idx_plt++){
    //These are not necessary but sometimes make it a bit more concise when making plots
    plot_lab = string_loop_label[idx_plt];
    sel_cat  = NamedFunc_loop[idx_plt];

    pm.Push<Hist1D>(Axis(70, 50, 120,  "ll_m[0]",       "m_{ll} [GeV]", {}),       sel_cat, procs, ops).Weight(wgt).Tag("ShortName:an_categorization_" + plot_lab + "_ll_m");
    pm.Push<Hist1D>(Axis(bins[idx_plt], 100, 180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), sel_cat, procs, ops).Weight(wgt).Tag("ShortName:an_categorization_" + plot_lab + "_llphoton_m");

    sel_cat  = NamedFunc_refit_loop[idx_plt];
    pm.Push<Hist1D>(Axis(bins[idx_plt], 100, 180, "llphoton_refit_m", "m_{ll#gamma,refit} [GeV]", {}), sel_cat, procs, ops).Weight(wgt).Tag("ShortName:an_categorization_" + plot_lab + "_llphoton_refit_m");


    pm.Push<Hist2D>(
      Axis(70,50,120,  "ll_m[0]", "m_{ll} [GeV]", {}),
      Axis(bins[idx_plt], 100, 180,  "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
      sel_cat, procs, ops_2D).Tag("ShortName:an_categorization_" + plot_lab + "_mll_mlly");

    //These plots include all selections (which does include mll selections)
    sel_cat = NamedFunc_wsel_loop[idx_plt];
    pm.Push<Hist1D>(Axis(70, 50, 120,  "ll_m[0]",       "m_{ll} [GeV]", {}),       sel_cat, procs, ops).Weight(wgt).Tag("ShortName:an_categorization_" + plot_lab + "_fullsel_ll_m");
    pm.Push<Hist1D>(Axis(bins[idx_plt], 100, 180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), sel_cat, procs, ops).Weight(wgt).Tag("ShortName:an_categorization_" + plot_lab + "_fullsel_llphoton_m");

    sel_cat  = NamedFunc_wsel_refit_loop[idx_plt];
    pm.Push<Hist1D>(Axis(bins[idx_plt], 100, 180, "llphoton_refit_m", "m_{ll#gamma,refit} [GeV]", {}), sel_cat, procs, ops).Weight(wgt).Tag("ShortName:an_categorization_" + plot_lab + "_fullsel_llphoton_refit_m");

    constructCutflowTable(category_tablerows[idx_plt], wgt, 47, true, cat_vec[idx_plt%6] && "ntrulep >= 2" );
    
    selection_tablerows.push_back(TableRow(plot_lab, sel_cat && "llphoton_refit_m > 120 && llphoton_refit_m < 130", 0, 1, wgt));
  }

  int cnt = 0; 
  for(vector<TableRow> tablerow : category_tablerows){
    pm.Push<Table>("an_cutflow_category_" + cat_vec_str[cnt++] + "_" + lep_lab[0], tablerow, procs, false,true,false,false,false,true);
  }
  pm.Push<Table>("r3_signal_yields", selection_tablerows, procs_sig, false,true,false,false,false,true);

  //This block controls the verbose or non-verbose plotting option and whether or not to use multiple threads
  pm.min_print_ = true;
  pm.multithreaded_ = true;

  //Luminosity option for the plots
  pm.SetLuminosityTag(lumi_label).MakePlots(1.0);
}

