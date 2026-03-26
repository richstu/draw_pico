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


int main() {
  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");

  vector<PlotOpt> ops = CRUtilities::an_plotting_options("shapes"); ops[0].Bottom(PlotOptTypes::BottomType::off).Title(PlotOptTypes::TitleType::simulation);
  vector<PlotOpt> ops_shape = CRUtilities::an_plotting_options("shapes"); ops_shape[0].Title(PlotOptTypes::TitleType::simulation);
  vector<PlotOpt> ops_2D = an_plotting_options("2D"); ops_2D[0].YAxis(YAxisType::log).LogMinimum(0.001).LogMaximum(10);
  vector<PlotOpt> ops_2D_bkg = an_plotting_options("2D"); ops_2D_bkg[0].YAxis(YAxisType::log).LogMinimum(0.1).LogMaximum(10000);
 
  vector<PlotOpt> ops_1D = CRUtilities::an_plotting_options("1D"); ops_1D[0].Title(PlotOptTypes::TitleType::simulation);
  
  //Run2
  string samples_file = CRUtilities::SamplesFile;
  vector<shared_ptr<Process>> procs        = ZgUtilities::procs_with_sig_scale(samples_file, "SignalAllYearsSplit", 1); 
  vector<shared_ptr<Process>> procs_r2vsr3 = ZgUtilities::procs_with_sig_scale(samples_file, "SignalRun2vsRun3", 1); 
  vector<shared_ptr<Process>> procs_svb    = ZgUtilities::procs_with_sig_scale(samples_file, "SignalBackgroundAllYears", 1); 
  vector<shared_ptr<Process>> procs_unskim = ZgUtilities::procs_with_sig_scale(samples_file, "SignalAllYearsSplitUnskimmed", 1); 
  vector<shared_ptr<Process>> procs_usbkg  = ZgUtilities::procs_with_sig_scale(samples_file, "SignalAllYearsBkg", 1); 

  //These are all the selections used as a part of the baseline selection from the zg_functions.cpp file
  NamedFunc baseline = tightened_baseline_pinnacles;
  
  //This block of code handles all of the plot making
  PlotMaker pm;

  NamedFunc selection = baseline && "w_pu <  10.";
  pm.Push<Hist1D>(Axis(40,   80,  100, "ll_m[0]",   "m_{ll} [GeV]",    {}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_sigkin_ll_m");
  pm.Push<Hist1D>(Axis(60,    0,  180, "ll_pt[0]",  "p_{T}(ll) [GeV]", {}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_sigkin_ll_pt");
  pm.Push<Hist1D>(Axis(50, -5.0,  5.0, "ll_eta[0]", "#eta(ll)",        {}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_sigkin_ll_eta");
  pm.Push<Hist1D>(Axis(63,-3.15, 3.15, "ll_phi[0]", "#phi(ll)",        {}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_sigkin_ll_phi");
 
  pm.Push<Hist1D>(Axis(80,    0,   80, "photon_pt[0]", "p_{T}(#gamma) [GeV]", {}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_sigkin_photon_pt");
  pm.Push<Hist1D>(Axis(50, -2.5,  2.5, "photon_eta[0]","#eta(#gamma)",        {}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_sigkin_photon_eta");
  pm.Push<Hist1D>(Axis(63,-3.15, 3.15, "photon_phi[0]","#phi(#gamma)",        {}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_sigkin_photon_phi");
  pm.Push<Hist1D>(Axis(80,    0,    1, "photon_idmva[0]", "#gamma IDMVA",     {}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_sigkin_photon_idmva");
  
  pm.Push<Hist1D>(Axis(80,  110,  140, "llphoton_refit_m",   "m_{ll#gamma} [GeV]",  {}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_sigkin_llphoton_m"); 
  pm.Push<Hist1D>(Axis(60,   0,   180, "llphoton_refit_pt",  "p_{T}(ll#gamma) [GeV]", {}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_sigkin_llphoton_pt");
  pm.Push<Hist1D>(Axis(50, -5.0,  5.0, "llphoton_refit_eta", "#eta(ll#gamma)",        {}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_sigkin_llphoton_eta");
  pm.Push<Hist1D>(Axis(63,-3.15, 3.15, "llphoton_refit_phi", "#phi(ll#gamma)",        {}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_sigkin_llphoton_phi");
  
  pm.Push<Hist1D>(Axis(50,   -1,    1, "llphoton_refit_cosTheta", "cos#Theta", {}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_sigkin_cosTheta");
  pm.Push<Hist1D>(Axis(50,   -1,    1, "llphoton_refit_costheta", "cos#theta", {}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_sigkin_costheta");
  pm.Push<Hist1D>(Axis(63,-3.15, 3.15, "llphoton_refit_psi",      "#phi",      {}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_sigkin_psi");


  vector<NamedFunc> vector_baseline = {(Zee || Zmumu), "pass && nll>0", pass_trigs_and_pt, "nphoton > 0",   "ll_m[0] > 80 && ll_m[0] < 100", "photon_pt[0]/llphoton_m[0] > 15.0/110",
                                       "llphoton_m[0] + ll_m[0] > 185", "llphoton_refit_m > 100 && llphoton_refit_m < 180"};
  //Baseline selection 
  vector<NamedFunc> baseline_progressive = NamedFuncUtilities::progressive_cuts(vector_baseline);
  pm.Push<Hist1D>(Axis(3,   0,  3, "nll",             "N_{ll}",                      {1}), baseline_progressive[0],                      procs_unskim, ops_1D).Weight(wgt).Tag("ShortName:an_sigkin_baseline_0");
  pm.Push<Hist1D>(Axis(2,   0,  2, pass_trigs_and_pt, "triggers + p_{T} thresholds", {1}), baseline_progressive[1], procs_unskim, ops_1D).Weight(wgt).Tag("ShortName:an_sigkin_baseline_1");
  pm.Push<Hist1D>(Axis(4,   0,  4, "nphoton",         "N_{#gamma}",                  {1}), baseline_progressive[2], procs_unskim, ops_1D).Weight(wgt).Tag("ShortName:an_sigkin_baseline_2");
  pm.Push<Hist1D>(Axis(60, 60,120, "ll_m[0]",         "m_{ll} [GeV]",           {80,100}), baseline_progressive[3], procs_unskim, ops_1D).Weight(wgt).Tag("ShortName:an_sigkin_baseline_3");
  
  pm.Push<Hist2D>(
    Axis(60, 60, 120, "ll_m[0]",       "m_{ll} [GeV]",       {80,100}),
    Axis(100, 70, 170, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {95}),
    baseline_progressive[3], procs_usbkg, ops_2D).Weight(wgt).Tag("ShortName:an_sigkin_2D_baseline3_signal");

  vector<shared_ptr<Process>> procs_bkg  = ZgUtilities::procs_with_sig_scale(samples_file, "BkgAllYears", 1); 
  pm.Push<Hist2D>(
    Axis(60, 60, 120, "ll_m[0]",       "m_{ll} [GeV]",       {80,100}),
    Axis(100, 70, 170, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {95}),
    baseline_progressive[3], procs_bkg, ops_2D_bkg).Weight(wgt).Tag("ShortName:an_sigkin_2D_baseline3_bkg");


  pm.Push<Hist1D>(Axis(50, 0, 0.5, "photon_pt[0]/llphoton_m[0]", "p_{T}(#gamma)/m_{ll#gamma}", {15.0/110}), baseline_progressive[4], procs_unskim, ops_1D).Weight(wgt).Tag("ShortName:an_sigkin_baseline_4");
 
  pm.Push<Hist1D>(Axis(100,  100,  200, "llphoton_refit_m",   "m_{ll#gamma} [GeV]",  {180}), baseline_progressive[5], procs_unskim, ops_1D).Weight(wgt).Tag("ShortName:an_sigkin_baseline_5"); 
  pm.Push<Hist1D>(Axis(100,  100,  200, "llphoton_refit_m",   "m_{ll#gamma} [GeV]",  {180}), baseline_progressive[6], procs_unskim, ops_1D).Weight(wgt).Tag("ShortName:an_sigkin_baseline_6"); 
 

  //Run 2 vs Run 3
  pm.Push<Hist1D>(Axis(40,   80,  100, "ll_m[0]",   "m_{ll} [GeV]",    {}), selection, procs_r2vsr3, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_r2vr3_ll_m");
  pm.Push<Hist1D>(Axis(60,    0,  180, "ll_pt[0]",  "p_{T}(ll) [GeV]", {}), selection, procs_r2vsr3, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_r2vr3_ll_pt");
  pm.Push<Hist1D>(Axis(50, -5.0,  5.0, "ll_eta[0]", "#eta(ll)",        {}), selection, procs_r2vsr3, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_r2vr3_ll_eta");
  pm.Push<Hist1D>(Axis(63,-3.15, 3.15, "ll_phi[0]", "#phi(ll)",        {}), selection, procs_r2vsr3, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_r2vr3_ll_phi");
 
  pm.Push<Hist1D>(Axis(80,    0,   80, "photon_pt[0]", "p_{T}(#gamma) [GeV]", {}), selection, procs_r2vsr3, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_r2vr3_photon_pt");
  pm.Push<Hist1D>(Axis(50, -2.5,  2.5, "photon_eta[0]","#eta(#gamma)",        {}), selection, procs_r2vsr3, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_r2vr3_photon_eta");
  pm.Push<Hist1D>(Axis(63,-3.15, 3.15, "photon_phi[0]","#phi(#gamma)",        {}), selection, procs_r2vsr3, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_r2vr3_photon_phi");
  pm.Push<Hist1D>(Axis(80,    0,    1, "photon_idmva[0]", "#gamma IDMVA",     {}), selection, procs_r2vsr3, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_r2vr3_photon_idmva");
  
  pm.Push<Hist1D>(Axis(80,  110,  140, "llphoton_refit_m",   "m_{ll#gamma} [GeV]",  {}), selection, procs_r2vsr3, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_r2vr3_llphoton_m"); 
  pm.Push<Hist1D>(Axis(60,   0,   180, "llphoton_refit_pt",  "p_{T}(ll#gamma) [GeV]", {}), selection, procs_r2vsr3, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_r2vr3_llphoton_pt");
  pm.Push<Hist1D>(Axis(50, -5.0,  5.0, "llphoton_refit_eta", "#eta(ll#gamma)",        {}), selection, procs_r2vsr3, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_r2vr3_llphoton_eta");
  pm.Push<Hist1D>(Axis(63,-3.15, 3.15, "llphoton_refit_phi", "#phi(ll#gamma)",        {}), selection, procs_r2vsr3, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_r2vr3_llphoton_phi");
 
  pm.Push<Hist1D>(Axis(50,   -1,  1, "llphoton_refit_cosTheta", "cos#Theta",          {}), selection, procs_r2vsr3, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_r2vr3_cosTheta");
  pm.Push<Hist1D>(Axis(50,   -1,  1, "llphoton_refit_costheta", "cos#theta",          {}), selection, procs_r2vsr3, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_r2vr3_costheta");
  pm.Push<Hist1D>(Axis(63,-3.15,3.15,"llphoton_refit_psi",      "#phi",               {}), selection, procs_r2vsr3, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_r2vr3_psi");


  //Signal vs Background
  selection = selection && higgs_window_refit;
  pm.Push<Hist1D>(Axis(40,   80,  100, "ll_m[0]",   "m_{ll} [GeV]",    {}), selection, procs_svb, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_svb_ll_m");
  pm.Push<Hist1D>(Axis(60,    0,  180, "ll_pt[0]",  "p_{T}(ll) [GeV]", {}), selection, procs_svb, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_svb_ll_pt");
  pm.Push<Hist1D>(Axis(50, -5.0,  5.0, "ll_eta[0]", "#eta(ll)",        {}), selection, procs_svb, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_svb_ll_eta");
  pm.Push<Hist1D>(Axis(63,-3.15, 3.15, "ll_phi[0]", "#phi(ll)",        {}), selection, procs_svb, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_svb_ll_phi");
 
  pm.Push<Hist1D>(Axis(80,    0,   80, "photon_pt[0]", "p_{T}(#gamma) [GeV]", {}), selection, procs_svb, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_svb_photon_pt");
  pm.Push<Hist1D>(Axis(50, -2.5,  2.5, "photon_eta[0]","#eta(#gamma)",        {}), selection, procs_svb, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_svb_photon_eta");
  pm.Push<Hist1D>(Axis(63,-3.15, 3.15, "photon_phi[0]","#phi(#gamma)",        {}), selection, procs_svb, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_svb_photon_phi");
  pm.Push<Hist1D>(Axis(80,    0,    1, "photon_idmva[0]", "#gamma IDMVA",     {}), selection, procs_svb, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_svb_photon_idmva");
 
  pm.Push<Hist1D>(Axis(80,  110,  140, "llphoton_refit_m",   "m_{ll#gamma} [GeV]",  {}), selection, procs_svb, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_svb_llphoton_m"); 
  pm.Push<Hist1D>(Axis(60,   0,   180, "llphoton_refit_pt",  "p_{T}(ll#gamma) [GeV]", {}), selection, procs_svb, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_svb_llphoton_pt");
  pm.Push<Hist1D>(Axis(50, -5.0,  5.0, "llphoton_refit_eta", "#eta(ll#gamma)",        {}), selection, procs_svb, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_svb_llphoton_eta");
  pm.Push<Hist1D>(Axis(63,-3.15, 3.15, "llphoton_refit_phi", "#phi(ll#gamma)",        {}), selection, procs_svb, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_svb_llphoton_phi");
 
  pm.Push<Hist1D>(Axis(50,   -1,  1, "llphoton_refit_cosTheta", "cos#Theta",          {}), selection, procs_svb, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_svb_cosTheta");
  pm.Push<Hist1D>(Axis(50,   -1,  1, "llphoton_refit_costheta", "cos#theta",          {}), selection, procs_svb, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_svb_costheta");
  pm.Push<Hist1D>(Axis(63,-3.15,3.15,"llphoton_refit_psi",      "#phi",               {}), selection, procs_svb, ops_shape).Weight(wgt).Tag("ShortName:an_sigkin_svb_psi");


  //This block controls the verbose or non-verbose plotting option and whether or not to use multiple threads
  pm.min_print_ = true;
  pm.multithreaded_ = true;

  //Luminosity option for the plots
  pm.SetLuminosityTag(CRUtilities::R23Lumi).MakePlots(1); 
}


