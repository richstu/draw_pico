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

//--This is a an_signal_shape file meant to be used to easily make plotting files

  std::unordered_map<std::string,double>  year_lumi{{"2016",19.5},{"2016APV",16.8},{"2017",41.48},{"2018",59.83}};


//Currently the an_signal_shape has the year as an input variable, but this can be removed
//This is only used in defining processes
//int main() {
int main() {

  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");

  PlotOpt lin_lumi("txt/plot_styles.txt","CMSPaper");
  lin_lumi.YAxis(YAxisType::linear)
          .Stack(StackType::shapes)
          .YTitleOffset(1.75)
          .AutoYAxis(false)
          .UseCMYK(false)
          .LeftMargin(0.17)
          .CanvasWidth(800)
          .CanvasHeight(800)
          .Title(TitleType::info)
          .FileExtensions({"pdf"});
  vector<PlotOpt> ops = {lin_lumi(), lin_lumi().Stack(StackType::signal_overlay)};


  //Uses txt/zg_samples.txt to define needed processes
  vector<shared_ptr<Process>> procs_sig_r2     =  ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_test_speed.txt","SigMCSplitRun2");
  vector<shared_ptr<Process>> procs_sig_2022   =  ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_test_speed.txt","SigMCSplit2022"); 
  vector<shared_ptr<Process>> procs_sig_2022EE =  ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_test_speed.txt","SigMCSplit2022EE"); 
  vector<vector<shared_ptr<Process>>> procs_vector = {procs_sig_r2, procs_sig_2022, procs_sig_2022EE};
  vector<string> procs_string = {"_procs_sig_r2", "_procs_sig_2022", "_procs_sig_2022EE"};



  //This defines the splits based on lepton flavor that will be plotted
  vector<NamedFunc> lep = { "ll_lepid[0]==11","ll_lepid[0]==13","1"};
  vector<string> lep_lab= {"_ee", "_mumu", "_ll"};

  //This defines the categorization that will be used for plotting
  vector<NamedFunc> cat_vec = CatUtilities::run3_category_vector;
  vector<NamedFunc> cat_wsel_vec = CatUtilities::run3_category_vector;
  vector<string> cat_vec_str= CatUtilities::run3_category_labels;                          

  //These are all the selections used as a part of the baseline selection from the zg_functions.cpp file
  NamedFunc tight_baseline = ZgFunctions::tightened_baseline;
 

  //This block of code loops through to create the NamedFuncs for each category
  //The purpose is to only require one loop while making plots
  vector<NamedFunc> NamedFunc_loop  = {};
  vector<NamedFunc> NF_catsel_loop  = {};
  vector<string>    string_loop_label_tight = {};
  vector<string>    string_loop_lab_sel_tight = {};

  for(unsigned int idx_i = 0; idx_i < lep.size(); idx_i++){
    for(unsigned int idx_j = 0; idx_j < cat_vec.size(); idx_j++){
      //Vectors used for plots with the tightened baseline selection
      NamedFunc_loop.push_back( tightened_baseline && cat_vec[idx_j] && lep[idx_i]);
      NF_catsel_loop.push_back( tightened_baseline && cat_wsel_vec[idx_j] && lep[idx_i]);
      string_loop_label_tight.push_back(cat_vec_str[idx_j] + lep_lab[idx_i]);
      string_loop_lab_sel_tight.push_back(cat_vec_str[idx_j] + lep_lab[idx_i] + "catsel");

    }
  }

  //This block of code handles all of the plot making
  PlotMaker pm;
  string pltNames = "all";
  NamedFunc sel_cat = tight_baseline;
  string pltNstot = "all";
  NamedFunc sel_tot = tight_baseline;

  int bins_mlly = 160;
  for(unsigned int idx_proc = 0; idx_proc < procs_vector.size(); idx_proc++){
    string procs_lab = procs_string[idx_proc];

    //Makes plots with the tightened baseline selection
    vector<vector<TableRow>> higgs_window_yields = {{},{},{}};
    vector<vector<TableRow>> higgs_window_yields_totsel = {{},{},{}};


    pm.Push<Hist1D>(Axis(bins_mlly, 100, 180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {120,130}), sel_cat, procs_vector[idx_proc], ops).Weight(wgt_nodiff).Tag("ShortName:an_signal_shape_" + pltNames+"_llphoton_m" + procs_lab);

    for(unsigned int idx_plt = 0; idx_plt < NamedFunc_loop.size(); idx_plt++){
      //These are not necessary but sometimes make it a bit more concise when making plots
      pltNames = string_loop_label_tight[idx_plt];
      pltNstot = string_loop_lab_sel_tight[idx_plt];
      sel_cat  = NamedFunc_loop[idx_plt];
      sel_tot  = NF_catsel_loop[idx_plt];

      //plots without category selections
      pm.Push<Hist1D>(Axis(bins_mlly, 100, 180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {120,130}), sel_cat, procs_vector[idx_proc], ops).Weight(wgt_nodiff).Tag("ShortName:an_signal_shape_" + pltNames+"_llphoton_m" + procs_lab);
      pm.Push<Hist1D>(Axis(bins_mlly, 100, 180, "llphoton_refit_m", "m_{ll#gamma,refit} [GeV]", {120,130}), sel_cat, procs_vector[idx_proc], ops).Weight(wgt_nodiff).Tag("ShortName:an_signal_shape_" + pltNames+"_llphoton_m" + procs_lab);

      //plots with category selections
      pm.Push<Hist1D>(Axis(bins_mlly, 100, 180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {120,130}), sel_tot, procs_vector[idx_proc], ops).Weight(wgt_nodiff).Tag("ShortName:an_signal_shape_" + pltNstot+"_llphoton_m" + procs_lab);
      pm.Push<Hist1D>(Axis(bins_mlly, 100, 180, "llphoton_refit_m", "m_{ll#gamma,refit} [GeV]", {120,130}), sel_tot, procs_vector[idx_proc], ops).Weight(wgt_nodiff).Tag("ShortName:an_signal_shape_" + pltNstot+"_llphoton_m" + procs_lab);



      //Adds rows to tables, CHANGE ROW NAME
      higgs_window_yields[static_cast<int>(idx_plt/6)].push_back(TableRow(pltNames, sel_cat, 0,0, wgt_nodiff));
      higgs_window_yields_totsel[static_cast<int>(idx_plt/6)].push_back(TableRow(pltNstot, sel_tot, 0,0, wgt_nodiff));

    }

    //Tables of yields for each final state signal MC
    pm.Push<Table>("an_signal_shape_ee"   + procs_lab, higgs_window_yields[0], procs_vector[idx_proc],0, 1, 0, 1, 0).Precision(3);
    pm.Push<Table>("an_signal_shape_mumu" + procs_lab, higgs_window_yields[1], procs_vector[idx_proc],0, 1, 0, 1, 0).Precision(3);
    pm.Push<Table>("an_signal_shape_ll"   + procs_lab, higgs_window_yields[2], procs_vector[idx_proc],0, 1, 0, 1, 0).Precision(3);

    pm.Push<Table>("an_signal_shape_ee"   + procs_lab, higgs_window_yields_totsel[0], procs_vector[idx_proc],0, 1, 0, 1, 0).Precision(3);
    pm.Push<Table>("an_signal_shape_mumu" + procs_lab, higgs_window_yields_totsel[1], procs_vector[idx_proc],0, 1, 0, 1, 0).Precision(3);
    pm.Push<Table>("an_signal_shape_ll"   + procs_lab, higgs_window_yields_totsel[2], procs_vector[idx_proc],0, 1, 0, 1, 0).Precision(3);

  }

  //This block controls the verbose or non-verbose plotting option and whether or not to use multiple threads
  pm.min_print_ = true;
  pm.multithreaded_ = true;

  //Luminosity option for the plots
  pm.SetLuminosityTag("1").MakePlots(1);

}
