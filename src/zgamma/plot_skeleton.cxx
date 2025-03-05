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
#include "zgamma/controlregion_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace ZgUtilities;
using namespace ZgFunctions;
using namespace CatUtilities;
using namespace CRUtilities;
//--This is a skeleton file meant to be used to easily make plotting files

std::unordered_map<std::string,double>  year_lumi{{"2016",19.5},{"2016APV",16.8},{"2017",41.48},{"2018",59.83}};

//Currently the skeleton has the year as an input variable, but this can be removed
//This is only used in defining processes
//int main() {
int main() {
  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");

  //Defines the plot options for 1D and 2D plots
  vector<PlotOpt> ops    = CRUtilities::an_plotting_options("1D");
  vector<PlotOpt> ops_2D = CRUtilities::an_plotting_options("2D");

  //Uses txt/zg_samples.txt to define needed processes
  
  //Run2
  vector<shared_ptr<Process>> procs = ZgUtilities::procs_with_sig_scale("txt/samples_zgamma_lassen.txt", "AllMoreTrigs", 10); 
  string year_label = "_run2";
  string era_label  = "137.61";

  //This defines the splits based on lepton flavor that will be plotted
  vector<NamedFunc> lep = { "ll_lepid[0]==11","ll_lepid[0]==13","1"};
  vector<string> lep_lab= {"_ee", "_mumu", "_ll"};
 
  //These are all the selections used as a part of the baseline selection from the zg_functions.cpp file
  NamedFunc tight_baseline = ZgFunctions::tightened_baseline;

  //This block of code handles all of the plot making
  PlotMaker pm;

  //Functions to shorten the command to add plots to PlotMaker
  string label = year_label;
  NamedFunc selection = tight_baseline;

  //Add plots to plot maker
  pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), selection, procs, ops).Weight(wgt).Tag("ShortName:Skeleton_" + label + "_llphoton_m");

  //This block controls the verbose or non-verbose plotting option and whether or not to use multiple threads
  pm.min_print_ = true;
  pm.multithreaded_ = true;

  //Luminosity option for the plots
  pm.SetLuminosityTag(era_label).MakePlots(1); 
}


