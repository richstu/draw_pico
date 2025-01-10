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

//--This is a an_peak_H_bkg file meant to be used to easily make plotting files

  std::unordered_map<std::string,double>  year_lumi{{"2016",19.5},{"2016APV",16.8},{"2017",41.48},{"2018",59.83}};


//Currently the an_peak_H_bkg has the year as an input variable, but this can be removed
//This is only used in defining processes
//int main() {
int main(int argc, char *argv[]) {

  bool run3_plots = false;
  string era_select = "";
  if(argc>1){ era_select=argv[1];}
  if(era_select == "run3"){run3_plots=true;}

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
  vector<shared_ptr<Process>> procs_allhiggs =  ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_zgamma.txt","AllHiggsRun2"); 
  vector<shared_ptr<Process>> procs_peakinghiggs =  ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_zgamma.txt","HiggsBkgRun2");
  string lumi_tag = "137.61";

  if(run3_plots){
    cout << "add the ability to plot run3 samples" << endl;
    procs_allhiggs     =  ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_zgamma.txt","AllHiggsRun3"); 
    procs_peakinghiggs =  ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_zgamma.txt","HiggsBkgRun3");
    lumi_tag = "62.32";
  }

  vector<vector<shared_ptr<Process>>> procs_vector = {procs_allhiggs,  procs_peakinghiggs};
  vector<string> procs_string = {"_procs_allhiggs", "_procs_peakinghiggs"};

  //This defines the splits based on lepton flavor that will be plotted
  vector<NamedFunc> lep = { "ll_lepid[0]==11","ll_lepid[0]==13","1"};
  vector<string> lep_lab= {"_ee", "_mumu", "_ll"};

  //This defines the categorization that will be used for plotting
  vector<NamedFunc> cat_vec = CatUtilities::run3_category_vector;
  vector<string> cat_vec_str= CatUtilities::run3_category_labels;                          

  //These are all the selections used as a part of the baseline selection from the zg_functions.cpp file
  NamedFunc tight_baseline = ZgFunctions::tightened_baseline;
 

  //This block of code loops through to create the NamedFuncs for each category
  //The purpose is to only require one loop while making plots
  vector<NamedFunc> NamedFunc_loop  = {};
  vector<string>    string_loop_label_tight = {};
  for(unsigned int idx_i = 0; idx_i < lep.size(); idx_i++){
    for(unsigned int idx_j = 0; idx_j < cat_vec.size(); idx_j++){
      //Vectors used for plots with the tightened baseline selection
      NamedFunc_loop.push_back( tightened_baseline && cat_vec[idx_j] && lep[idx_i]);
      string_loop_label_tight.push_back(cat_vec_str[idx_j] + lep_lab[idx_i]);
    }
  }

  //This block of code handles all of the plot making
  PlotMaker pm;
  string pltNames = "all";
  NamedFunc selection = tight_baseline;


  pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {120,130}), selection, procs_vector[0], ops).Weight(wgt_nodiff).Tag("ShortName:an_peak_H_bkg_nocat_llphoton_m_allh");
  pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {120,130}), selection, procs_vector[1], ops).Weight(wgt_nodiff).Tag("ShortName:an_peak_H_bkg_nocat_llphoton_m_bkgh");


  for(unsigned int idx_proc = 0; idx_proc < procs_vector.size(); idx_proc++){
    string procs_lab = procs_string[idx_proc];

    //Makes plots with the tightened baseline selection
    vector<vector<TableRow>> higgs_window_yields = {{},{},{}};
    //pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {120,130}), selection, procs, ops).Weight(wgt_nodiff).Tag("ShortName:an_peak_H_bkg_" + pltNames+"_llphoton_m" + procs_lab);

    for(unsigned int idx_plt = 0; idx_plt < NamedFunc_loop.size(); idx_plt++){
      //These are not necessary but sometimes make it a bit more concise when making plots
      pltNames = string_loop_label_tight[idx_plt];
      selection  = NamedFunc_loop[idx_plt];

      pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {120,130}), selection, procs_vector[idx_proc], ops).Weight(wgt_nodiff).Tag("ShortName:an_peak_H_bkg_" + pltNames+"_llphoton_m" + procs_lab);
      pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_refit_m", "m_{ll#gamma,refit} [GeV]", {120,130}), selection, procs_vector[idx_proc], ops).Weight(wgt_nodiff).Tag("ShortName:an_peak_H_bkg_" + pltNames+"_llphoton_m" + procs_lab);

      //Adds rows to tables, CHANGE ROW NAME
      higgs_window_yields[static_cast<int>(idx_plt/6)].push_back(TableRow(pltNames, selection, 0,0, wgt_nodiff));
    }

    //Tables of yields for each final state signal MC
    pm.Push<Table>("an_peaking_higgs_bkg_ee" + procs_lab,   higgs_window_yields[0], procs_vector[idx_proc],0, 1, 0, 1, 0).Precision(3);
    pm.Push<Table>("an_peaking_higgs_bkg_mumu" + procs_lab, higgs_window_yields[1], procs_vector[idx_proc],0, 1, 0, 1, 0).Precision(3);
    pm.Push<Table>("an_peaking_higgs_bkg_ll" + procs_lab,   higgs_window_yields[2], procs_vector[idx_proc],0, 1, 0, 1, 0).Precision(3);
  }

  //This block controls the verbose or non-verbose plotting option and whether or not to use multiple threads
  pm.min_print_ = true;
  pm.multithreaded_ = true;

  //Luminosity option for the plots
  pm.SetLuminosityTag(lumi_tag).MakePlots(1);

}
