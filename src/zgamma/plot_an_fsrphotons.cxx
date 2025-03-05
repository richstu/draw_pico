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
#include "core/named_func_utilities.hpp"
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

NamedFunc llfsrg_m("llfsrg_m",[](const Baby &b) -> NamedFunc::ScalarType{
  std::map<unsigned int, TLorentzVector> fsrphotons = fsrphoton_ret(b);
  TLorentzVector ll = AssignZ(b);
  return (ll + fsrphotons[0] + fsrphotons[1]).M();
});

NamedFunc llgfsrg_m("llfsrg_m",[](const Baby &b) -> NamedFunc::ScalarType{
  std::map<unsigned int, TLorentzVector> fsrphotons = fsrphoton_ret(b);
  TLorentzVector ll = AssignZ(b);
  TLorentzVector photon = AssignGamma(b);

  return (ll + photon + fsrphotons[0] + fsrphotons[1]).M();
});



//Currently the skeleton has the year as an input variable, but this can be removed
//This is only used in defining processes
//int main() {
int main() {
  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");

  vector<PlotOpt> ops = CRUtilities::an_plotting_options("basic_1D");

  //Uses txt/zg_samples.txt to define needed processes
  
  //Run2
  vector<shared_ptr<Process>> procs = ZgUtilities::procs_with_sig_scale("txt/samples_zgamma_lassen.txt", "SignalSplit", 1); 

  //These are all the selections used as a part of the baseline selection from the zg_functions.cpp file
  NamedFunc baseline = NamedFuncUtilities::Nminusk(ZgFunctions::vector_tightened_baseline,{5}) && "ll_lepid[0]==13";

  //This block of code handles all of the plot making
  PlotMaker pm;

  NamedFunc selection = baseline && "nfsrphoton>0";
  pm.Push<Hist1D>(Axis(80, 100, 180, "llphoton_m[0]",   "m_{ll#gamma} [GeV]",                 {}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_fsrphoton_llphoton_m");
  pm.Push<Hist1D>(Axis(50,  50, 100, "ll_m[0]",         "m_{ll} [GeV]",                       {}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_fsrphoton_ll_m");
  pm.Push<Hist1D>(Axis(30,   0,  30, "fsrphoton_pt[0]", "p_{T}(#gamma_{FSR}) [GeV]",          {}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_fsrphoton_fsrphoton_pt");
  pm.Push<Hist1D>(Axis(80,  50, 140, llfsrg_m,          "m_{ll#gamma} w/ FSR Recovery [GeV]", {}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_fsrphoton_llfsrphoton_m");
  pm.Push<Hist1D>(Axis(50, 100, 180, llgfsrg_m,         "m_{ll}       w/ FSR Recovery [GeV]", {}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_fsrphoton_llphotonfsrphoton_m");

  //This block controls the verbose or non-verbose plotting option and whether or not to use multiple threads
  pm.min_print_ = true;
  pm.multithreaded_ = true;

  //Luminosity option for the plots
  pm.SetLuminosityTag("137.61").MakePlots(1); 
}


