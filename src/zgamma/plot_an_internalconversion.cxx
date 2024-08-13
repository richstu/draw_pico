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

  NamedFunc mom_is_photon("mom_is_photon", [](const Baby &b) -> NamedFunc::ScalarType{
    for (unsigned int idx_i = 0; idx_i < b.mc_id()->size(); idx_i++){
      if(abs(b.mc_id()->at(idx_i)) < 16 && b.mc_mom() -> at(idx_i)==22){
        return true;
      }
    }
    return false;
  });


  NamedFunc ll_mom_is_Z("ll_mom_is_Z", [](const Baby &b) -> NamedFunc::ScalarType{
    for (unsigned int idx_i = 0; idx_i < b.mc_id()->size(); idx_i++){
      if( (abs(b.mc_id()->at(idx_i)) == 11 || abs(b.mc_id()->at(idx_i)) == 13) && b.mc_mom() -> at(idx_i)==23){
        return true;
      }
    }
    return false;
  });

  NamedFunc id_mom_is_photon("id_mom_is_photon", [](const Baby &b) -> NamedFunc::ScalarType{
    for (unsigned int idx_i = 0; idx_i < b.mc_id()->size(); idx_i++){
      if(abs(b.mc_id()->at(idx_i)) < 16 && b.mc_mom() -> at(idx_i)==22){
        return abs(b.mc_id() -> at(idx_i));
      }
    }
    return -1;
  });

  NamedFunc photon_from_hard("photon_from_hard", [](const Baby &b) -> NamedFunc::ScalarType{
    for (unsigned int idx_i = 0; idx_i < b.mc_id()->size(); idx_i++){
      int bitmap = b.mc_statusflag() ->at(idx_i);

      if(abs(b.mc_id()->at(idx_i)) == 22 && (CatUtilities::checkBit(bitmap,8))){
      //b.mc_statusflag() -> at(idx_i))[8]){
        return true;
      }
    }
    return false;
  });

//Currently the an_peak_H_bkg has the year as an input variable, but this can be removed
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
  vector<shared_ptr<Process>> procs =  ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_zgamma.txt","SignalUnskimmed"); 
  //vector<shared_ptr<Process>> procs =  ZgUtilities::procs_with_sig_scale("txt/samples_zgamma.txt","OtherHiggsBkg",1);

  //This defines the splits based on lepton flavor that will be plotted
  vector<NamedFunc> lep = {ll_mom_is_Z};// "ll_lepid[0]==11","ll_lepid[0]==13","1"};
  vector<string> lep_lab= {"_ll"};//ee", "_mumu", "_ll"};

  //These are all the selections used as a part of the baseline selection from the zg_functions.cpp file
  NamedFunc tight_baseline = "1";//ll_mom_is_Z; //nlep >= 2";
 
  //This block of code loops through to create the NamedFuncs for each category
  //The purpose is to only require one loop while making plots
  vector<NamedFunc> NamedFunc_loop  = {};
  vector<string>    string_loop_label_tight = {};
  for(unsigned int idx_i = 0; idx_i < lep.size(); idx_i++){
    //Vectors used for plots with the tightened baseline selection
    NamedFunc_loop.push_back( tightened_baseline  && lep[idx_i]);
    string_loop_label_tight.push_back(lep_lab[idx_i]);
  }

  //This block of code handles all of the plot making
  PlotMaker pm;
  string pltNames = "all";
  NamedFunc sel_cat = tight_baseline;

  //Makes plots with the tightened baseline selection
  vector<TableRow> internal_conversion_table = {};
  for(unsigned int idx_plt = 0; idx_plt < NamedFunc_loop.size(); idx_plt++){
    //These are not necessary but sometimes make it a bit more concise when making plots
    pltNames = string_loop_label_tight[idx_plt];
    sel_cat  =  mom_is_photon;

    cout << pltNames << endl;
    pm.Push<Hist1D>(Axis(17, -1, 16, id_mom_is_photon, "MC PDG ID", {}), sel_cat, procs, ops).Weight(wgt_nodiff).Tag("ShortName:an_intconv_" + pltNames+"_pdg_id__mom_is_ph");

    //Adds rows to tables, CHANGE ROW NAME
    internal_conversion_table.push_back(TableRow(pltNames + "mom is photon", sel_cat, 0,0, wgt_nodiff));
    internal_conversion_table.push_back(TableRow(pltNames + "photon from hard", ll_mom_is_Z && !photon_from_hard, 0,0, wgt_nodiff));

    sel_cat  = NamedFunc_loop[idx_plt] && !mom_is_photon;
    pm.Push<Hist1D>(Axis(17, -1, 16, id_mom_is_photon, "MC PDG ID", {}), sel_cat, procs, ops).Weight(wgt_nodiff).Tag("ShortName:an_intconv_" + pltNames+"_pdg_id__mom_not_ph");

    //Adds rows to tables, CHANGE ROW NAME
    internal_conversion_table.push_back(TableRow(pltNames + "mom is not photon", sel_cat, 0,0, wgt_nodiff));
    internal_conversion_table.push_back(TableRow(pltNames + "photon from hard", ll_mom_is_Z && !photon_from_hard, 0,0, wgt_nodiff));
  }

  //Tables of yields for each final state signal MC
  pm.Push<Table>("an_internal_conversion_ll",   internal_conversion_table, procs, 0, 1, 0, 1, 0).Precision(3);

  //This block controls the verbose or non-verbose plotting option and whether or not to use multiple threads
  pm.min_print_ = true;
  pm.multithreaded_ = true;

  //Luminosity option for the plots
  pm.SetLuminosityTag("137.61").MakePlots(1);

}
