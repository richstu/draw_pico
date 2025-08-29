//script used to do H->Zgamma trigger studies

#include "core/test.hpp"

#include <algorithm>
#include <bitset>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>
#include <getopt.h>

#include "Math/Vector4D.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TError.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TMath.h"
#include "TPad.h"
#include "TPie.h"
#include "TVector2.h"

#include "core/baby.hpp"
#include "core/cross_sections.hpp"
#include "core/efficiency_plot.hpp"
#include "core/event_scan.hpp"
#include "core/functions.hpp"
#include "core/gamma_params.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/named_func.hpp"
#include "core/named_func_utilities.hpp"
#include "core/palette.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/process.hpp"
#include "core/table.hpp"
#include "core/utilities.hpp"
#include "zgamma/apply_zg_trigeffs.hpp"
#include "zgamma/nano_functions.hpp"
#include "zgamma/zg_functions.hpp"
#include "zgamma/zg_utilities.hpp"

//Update 25.02.2025
//Below this is 6000 lines of several years of trigger studies from the dawn
//of time. Except for the pieces generated from NanoAOD, they are almost
//certainly incompatible with all n-tuples from the last N years. I'm leaving
//them for the record, but be warned very few things are expected to work
//-MO

//This script generates some tables for run3 H->ZGamma trigger studies
//
//Arguments
// --single_thread (-s) to run single thread for debugging
// --unblind (-u) to unblind (not done for AN plots)
// --year (-y) yearname to select data year (2016, 2017, 2018, run2)
// --tag (-t) to add a tag to produced plots
// --string_options (-o) to specify what to plot
//
//possible string options: 
//

using namespace std;
using namespace PlotOptTypes;
using namespace NanoFunctions;
using namespace NamedFuncUtilities;
using namespace ZgUtilities;
using ZgFunctions::w_years;
using ZgFunctions::Ztoee;
using ZgFunctions::ZtoMuMu;

//declare helper functions
//these are defined at the end of the file, after main()
/*! Output efficiencies in C++ format
 
  \param[in] pm pointer to PlotMaker that made efficiency plots
  \param[in] func_name name of function in output file
  \param[in] plot_names vector of names of efficiency plots
  \param[in] out_file_name name of file to append to
 */
void generate_2d_efficiencies(PlotMaker* pm, std::string func_name, 
    std::vector<std::string> plot_names, std::string out_file_name, 
    std::string x_var_name, std::vector<double> x_var_bins,
    std::string y_var_name, std::vector<double> y_var_bins);

//manual bring this in for now
namespace script_utilities { 
  struct ArgStruct {
    bool single_thread;
    std::string year_string;
    std::string string_options;
    std::string tag;
    bool unblind;
  };

  ArgStruct get_options(int argc, char* argv[], std::string default_string_options) {
    //
    //initialize to defaults
    ArgStruct options;
    options.single_thread = false;
    options.year_string = "run2";
    options.string_options = default_string_options;
    options.tag = "";
    options.unblind = false;

    //parse arguments
    while(true){
      static struct option long_options[] = {
        {"single_thread", no_argument, 0, 's'},
        {"year", required_argument, 0, 'y'},
        {"unblind", no_argument, 0, 'u'},
        {"tag", required_argument, 0, 't'},
        {"string_options", required_argument, 0, 'o'},
        {0, 0, 0, 0}
      };

      char opt = -1;
      int option_index;
      opt = getopt_long(argc, argv, "sy:ut:o:", long_options, &option_index);
      if (opt == -1) break;

      std::string optname;
      switch(opt){
      case 's':
        options.single_thread = true;
        break;
      case 'u':
        options.unblind = true;
        break;
      case 'y':
        options.year_string = optarg;
        break;
      case 't':
        options.tag = optarg;
        break;
      case 'o':
        options.string_options = optarg;
        break;
      case 0:
        //handle cases with no short argument form
        std::cout << "Bad option! Found option name " << optname << std::endl;
        break;
      default:
        //std::cout << "Bad option! getopt_long returned character code 0" << opt << std::endl;
        printf("Bad option! getopt_long returned character code 0%o\n", opt);
        break;
      }
    }
    return options;
  }

  std::vector<PlotOpt> plot_lin(bool unblind) {
    PlotOpt lin_norm("txt/plot_styles.txt", "CMSPaper");
    lin_norm.Title(PlotOptTypes::TitleType::info)   
        .Bottom(PlotOptTypes::BottomType::off)
        .YAxis(PlotOptTypes::YAxisType::linear)
        .Stack(PlotOptTypes::StackType::signal_overlay)
        .Overflow(PlotOptTypes::OverflowType::none)
        .LegendColumns(3);
    if (unblind) {
      lin_norm.Stack(PlotOptTypes::StackType::data_norm);
      lin_norm.Bottom(PlotOptTypes::BottomType::ratio);
    }
    return {lin_norm};
  }

  std::vector<PlotOpt> plot_log(bool unblind) {
    PlotOpt log_norm("txt/plot_styles.txt", "CMSPaper");
    log_norm.Title(PlotOptTypes::TitleType::info)   
        .Bottom(PlotOptTypes::BottomType::off)
        .YAxis(PlotOptTypes::YAxisType::log)
        .LogMinimum(.2)
        .Stack(PlotOptTypes::StackType::signal_overlay)
        .Overflow(PlotOptTypes::OverflowType::none)
        .LegendColumns(3);
    if (unblind) {
      log_norm.Stack(PlotOptTypes::StackType::data_norm);
      log_norm.Bottom(PlotOptTypes::BottomType::ratio);
    }
    return {log_norm};
  }

  std::vector<PlotOpt> plot_shapes() {
    PlotOpt shapes_norm("txt/plot_styles.txt", "CMSPaper");
    shapes_norm.Title(PlotOptTypes::TitleType::info)   
        .Bottom(PlotOptTypes::BottomType::off)
        .YAxis(PlotOptTypes::YAxisType::linear)
        .Stack(PlotOptTypes::StackType::shapes)
        .LegendColumns(3);
    return {shapes_norm};
  }

  std::vector<PlotOpt> plot_log_shapes() {
    PlotOpt shapes_norm("txt/plot_styles.txt", "CMSPaper");
    shapes_norm.Title(PlotOptTypes::TitleType::info)   
        .Bottom(PlotOptTypes::BottomType::off)
        .YAxis(PlotOptTypes::YAxisType::log)
        .LogMinimum(.2)
        .Stack(PlotOptTypes::StackType::shapes)
        .LegendColumns(3);
    return {shapes_norm};
  }
}

//returns min(max(0, value),1)
double fix_prob(double value) {
  if (value < 0)
    return 0;
  if (value > 1)
    return 1;
  return value;
}

//count number of times a value appears in a vector
template <class vec_type>
unsigned vector_count(std::vector<vec_type> vec, vec_type value) {
  unsigned counter = 0;
  for (unsigned i = 0; i < vec.size(); i++) {
    if (vec[i]==value)
      counter++;
  }
  return counter;
}

//-----------------------------------------------------------------------------
//                 named funcs to be moved to a common location
//-----------------------------------------------------------------------------

//returns a vector NamedFunc that is the list of indices from list2 matched by deltaR to list1
//can exclude elements of list2 using NamedFunc list2_mask
NamedFunc DeltaRMatchNamedFunc(NamedFunc list1_eta, NamedFunc list1_phi, NamedFunc list2_eta, NamedFunc list2_phi, NamedFunc list2_mask) {
  return NamedFunc("DeltaRMatchNamedFunc("+list1_eta.Name().substr(0,list1_eta.Name().size()-4)+")",[list1_eta, list1_phi, list2_eta, list2_phi, list2_mask](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> list2_indices;
    std::vector<double> list1_eta_ = list1_eta.GetVector(b);
    std::vector<double> list1_phi_ = list1_phi.GetVector(b);
    std::vector<double> list2_eta_ = list2_eta.GetVector(b);
    std::vector<double> list2_phi_ = list2_phi.GetVector(b);
    std::vector<double> list2_mask_ = list2_mask.GetVector(b);
    for (unsigned i1 = 0; i1 < list1_eta_.size(); i1++) {
      double index = -1;
      float min_dr = 999;
      for (unsigned i2 = 0; i2 < list2_eta_.size(); i2++) {
        if (list2_mask_[i2]) {
          float dr = deltaR(list1_eta_[i1],list1_phi_[i1],list2_eta_[i2],list2_phi_[i2]);
          if (dr < min_dr) {
            min_dr = dr;
            index = static_cast<double>(i2);
          }
        }
      }
      list2_indices.push_back(index);
    }
    return list2_indices;
  });
}

////Returns a scalar NamedFunc that is element index of vector_named_func or -999 if no index index
////only considering indices for which mask_named_func is true
//NamedFunc SortNamedFunc(NamedFunc vector_named_func, NamedFunc sort_criteria) {
//  return NamedFunc("SortNamedFunc("+vector_named_func.Name()+","+mask_named_func.Name()+")",[vector_named_func, sort_criteria](const Baby &b) -> NamedFunc::VectorType{
//    std::vector<double> vector_named_func_ = vector_named_func.GetVector(b);
//    std::vector<double> sort_criteria_ = sort_criteria.GetVector(b);
//    std::vector<std::pair<double,double>> sorting_vector;
//    std::vector<double> return_vector;
//    for (unsigned i = 0; i < vector_named_func_.size(); i++) {
//      sorting_vector.push_back(std::pair<double,double>(vector_named_func_, sort_criteria_));
//    }
//    std::sort(sorting_vector.begin(), sorting_vector.end(), [](std::pair<double,double> a, std::pair<double,double> b){a[1]>b[1]});
//    for (unsigned i = 0; i < vector_named_func_.size(); i++) {
//      return_vector.push_back(sorting_vector[i][0]);
//    }
//    return -999;
//  });
//}

int main(int argc, char *argv[]){
  //---------------------------------------------------------------------------
  //                                    constants
  //---------------------------------------------------------------------------
  //can only go as low as Electron/Muon_sig
  //const std::vector<double> el_pt_bins = {7.0, 10.0, 12.0, 15.0, 20.0, 25.0, 27.0, 30.0, 35.0, 40.0, 50.0, 80.0, 120.0, 200.0};
  const std::vector<double> el_pt_bins = {7.0, 12.0, 15.0, 20.0, 25.0, 27.0, 30.0, 35.0, 40.0, 50.0, 80.0, 120.0, 200.0};
  const std::vector<double> el_pt_bins_ptcut = {7.0, 15.0, 20.0, 25.0, 27.0, 30.0, 35.0, 40.0, 50.0, 80.0, 120.0, 200.0};
  const std::vector<double> el_abseta_bins = {0.0, 0.8, 1.4442, 1.566, 2.55};
  const std::vector<double> mu_pt_bins = {5.0, 7.0, 10.0, 15.0, 20.0, 23.0, 25.0, 30.0, 40.0, 50.0, 60.0, 120.0, 300.0};
  const std::vector<double> mu_pt_bins_ptcut = {5.0, 17.0, 20.0, 23.0, 25.0, 30.0, 40.0, 50.0, 60.0, 120.0, 300.0}; //for pt cut efficiencies, can only go down to 10 GeV due to cut on muon TrigObj
  const std::vector<double> mu_abseta_bins = {0.0, 0.9, 1.2, 2.1, 2.4};

  //---------------------------------------------------------------------------
  //                                    initialization
  //---------------------------------------------------------------------------

  //remove 120-130 mass range
  const NamedFunc blind_sr("blind_sr",[](const Baby &b) -> NamedFunc::ScalarType{
    double mllg = HiggsCandidate_mass.GetScalar(b);
    if (mllg>120 && mllg<130) return 0;
    return 1;
  });

  //remove events failing golden json
  GoldenJsonLoader golden_json_loader;
  NamedFunc pass_json = golden_json_loader.pass_json();

  gErrorIgnoreLevel = 6000;
  time_t begtime, endtime;
  time(&begtime);
  Palette colors("txt/colors.txt","default");
  script_utilities::ArgStruct options = script_utilities::get_options(
      argc, argv, "");
  std::vector<PlotOpt> plt_lin = script_utilities::plot_lin(false);
  std::vector<PlotOpt> plt_lin_unblind = script_utilities::plot_lin(true);
  std::vector<PlotOpt> plt_lin_over = script_utilities::plot_lin(false);
  plt_lin_over[0] = plt_lin_over[0].Overflow(OverflowType::both);
  std::vector<PlotOpt> plt_log = script_utilities::plot_log(false);
  std::vector<PlotOpt> plt_shapes = script_utilities::plot_shapes();
  std::vector<PlotOpt> plt_log_shapes = script_utilities::plot_log_shapes();
  PlotOpt style2D("txt/plot_styles.txt", "Scatter");
  std::vector<PlotOpt> twodim_plotopts = {style2D().Title(TitleType::info)
      .YAxis(YAxisType::linear).Overflow(OverflowType::overflow)
      .LogMinimum(0.001)};
  std::vector<PlotOpt> plt_lin_logx = {PlotOpt("txt/plot_styles.txt","CMSPaper")
      .Title(TitleType::info).Bottom(PlotOptTypes::BottomType::off)
      .XAxis(PlotOptTypes::YAxisType::log)
      .Overflow(PlotOptTypes::OverflowType::none).LegendColumns(3)};

  vector<PlotOpt> ops_mc = {PlotOpt("txt/plot_styles.txt","LinLumi")
      .Title(TitleType::private_work)}; 
  vector<PlotOpt> ops_shapes = {PlotOpt("txt/plot_styles.txt","Shapes")
      .Title(TitleType::private_work)}; 

  //set<int> years;
  //HigUtilities::parseYears(options.year_string, years);
  //int lumi_precision = 0;
  //string total_luminosity_string = HigUtilities::getLuminosityString(options.year_string, lumi_precision);
  std::string total_luminosity_string = "41.5";

  std::string base_folder = "/net/cms17/cms17r0/pico/NanoAODv7/nano/2017/signal/";
  std::string base_folder_ul = "/net/cms17/cms17r0/pico/NanoAODv2/nano/2017/";
  std::string base_folder_ulv9 = "/net/cms17/cms17r0/pico/NanoAODv9/nano/2016/";
  std::string base_folder_ulv9_2017 = "/net/cms17/cms17r0/pico/NanoAODv9/nano/2017/";
  std::string base_folder_v7data = "/net/cms17/cms17r0/pico/NanoAODv7/nano/2017/data/";
  std::string base_folder_v7 = "/net/cms17/cms17r0/pico/NanoAODv7/nano/";
  std::string base_folder_v9 = "/net/cms17/cms17r0/pico/NanoAODv9/nano/";

  std::vector<std::shared_ptr<Process>> procs;
  procs.push_back(Process::MakeShared<Baby_nano>("HToZG", 
      Process::Type::signal, kBlack,
      {base_folder+"*GluGluHToZG_ZToLL*"},"1"));

  //2017 triggers
  NamedFunc hlt_el_trigger = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL||HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ";
  NamedFunc hlt_mu_trigger = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8||HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8";
  NamedFunc hlt_single_el_trigger = "HLT_Ele35_WPTight_Gsf||HLT_Ele32_WPTight_Gsf_L1DoubleEG";
  NamedFunc hlt_single_mu_trigger = "HLT_IsoMu27||HLT_IsoMu24";
  //
  NamedFunc hlt_el_trigger_plus = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL||HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ||HLT_Ele35_WPTight_Gsf||HLT_Ele32_WPTight_Gsf_L1DoubleEG";
  NamedFunc hlt_mu_trigger_plus = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8||HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8||HLT_IsoMu27||HLT_IsoMu24";
  NamedFunc hlt_mu_trigger_plusplus = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8||HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8||HLT_IsoMu27||HLT_IsoMu24||HLT_Mu17_Photon30_IsoCaloId";
  NamedFunc l1_el_trigger = "L1_SingleEG24||L1_SingleEG34er2p1||L1_SingleIsoEG24er2p1||L1_SingleIsoEG24||L1_DoubleEG_18_17||L1_DoubleEG_22_10||L1_DoubleEG_LooseIso23_10||L1_DoubleEG_LooseIso24_10";
  NamedFunc l1_mu_trigger = "L1_DoubleMu_12_5||L1_DoubleMu_15_7_SQ||L1_DoubleMu_15_7_SQ_Mass_Min4";
  NamedFunc l1_mu_trigger_singlemu = "L1_SingleMu22||L1_SingleMu25";
  NamedFunc hlt_double_el_noiso_trigger = "HLT_DoubleEle25_CaloIdL_MW";
  NamedFunc hlt_double_mu_noiso_trigger = "HLT_Mu37_TkMu27";
  NamedFunc hlt_single_el_noiso_trigger = "HLT_Ele115_CaloIdVT_GsfTrkIdT";
  NamedFunc hlt_single_mu_noiso_trigger = "HLT_Mu50||HLT_Mu55";
  NamedFunc hlt_single_el_prescale_trigger = "HLT_Ele25_WPTight_Gsf||HLT_Ele27_WPTight_Gsf||HLT_Ele28_WPTight_Gsf||HLT_Ele32_WPTight_Gsf||HLT_Ele20_WPLoose_Gsf||HLT_Ele45_WPLoose_Gsf||"
                                             "HLT_Ele25_eta2p1_WPTight_Gsf||HLT_Ele27_eta2p1_WPTight_Gsf||HLT_Ele20_eta2p1_WPLoose_Gsf||HLT_Ele25_eta2p1_WPLoose_Gsf||HLT_Ele27_eta2p1_WPLoose_Gsf";
  NamedFunc hlt_single_mu_prescale_trigger = "HLT_IsoMu20||HLT_IsoMu22||HLT_IsoMu22_eta2p1||HLT_IsoMu24_eta2p1||HLT_IsoTkMu20||HLT_IsoTkMu22||HLT_IsoTkMu24||HLT_Mu45_eta2p1||HLT_TkMu50";
  NamedFunc hlt_mu_ph_trigger = "HLT_Mu17_Photon30_IsoCaloId";
  NamedFunc hlt_doubleph_trigger = "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55||HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90";
  NamedFunc l1_lep_trigger = l1_el_trigger||l1_mu_trigger;
  NamedFunc hlt_lep_trigger = hlt_el_trigger||hlt_mu_trigger;
  NamedFunc hlt_lep_trigger_plus = hlt_el_trigger_plus||hlt_mu_trigger_plus;
  NamedFunc l1_mu_eg_trigger = "L1_Mu5_EG23||L1_Mu5_LooseIsoEG20";

  //2018 triggers
  NamedFunc hlt_singleel_trigger_2018 = "HLT_Ele32_WPTight_Gsf";
  NamedFunc hlt_singlemu_trigger_2018 = "HLT_IsoMu27";
  NamedFunc hlt_doubleel_trigger_2018 = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL||HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ";
  NamedFunc hlt_doublemu_trigger_2018 = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8";
  NamedFunc hlt_doubleph_trigger_2018 = "HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto||HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90";
  NamedFunc l1_eg_trigger_2018 = "L1_SingleEG26er2p5||L1_SingleIsoEG24er2p1||L1_SingleIsoEG26er2p5||L1_DoubleEG_22_10_er2p5";
  NamedFunc l1_doublemu_trigger_2018 = "L1_DoubleMu_12_5||L1_DoubleMu_15_7";
  NamedFunc l1_mu_eg_trigger_2018 = "L1_Mu5_EG23er2p5||L1_Mu5_LooseIsoEG20er2p5";
  //l1 single mu, hlt mu+ph same

  std::vector<std::shared_ptr<Process>> full_procs;
  //full_procs.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+jets", 
  //    Process::Type::background, TColor::GetColor("#ffb400"),
  //    {base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-"
  //    "pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000*"},"1"));
  //full_procs.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+#gamma", 
  //    Process::Type::background, TColor::GetColor("#16bac5"),
  //    {base_folder_ul+"mc/ZGToLLG_01J_5f_TuneCP5_*"},"1"));
  //full_procs.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+jets", 
  //    Process::Type::background, TColor::GetColor("#ffb400"),
  //    {base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__035CF0D8-DDC5-6D46-B172-3BDC98CC624D*",
  //    base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__223A6F1C-C8C0-9F43-94C3-233AD51ABD32*",
  //    base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__29EF0133-B16A-E747-9E02-89BC6682B669*", //temp?
  //    base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__25118058-0F30-8A41-B82B-C7EAFEE55379*", //temp?
  //    base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__24950148-3722-B84D-9050-C996F5205F0A*"},"1"));
  full_procs.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+jets", 
      Process::Type::background, TColor::GetColor("#ffb400"),
      {base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__0*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__1*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__2*"},"1"));
      //{base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__049988F3-37E5-0346-ABFA-252C0CC3F235.root",
      //base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__04D4C0DC-7E25-784E-B9F9-4FCBD0D02E4A.root",
      //base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__08224897-3E25-064E-A735-5ED118DB12CB.root",
      //base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__0979AFCD-1562-D540-9CF3-513D158A14ED.root",
      //base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__0F6EFF55-2257-4C4E-B26A-DB4FADB06F91.root",
      //base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__135273C9-F6AB-444E-B7C2-8AFE427CB644.root",
      //base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__142F2BF0-9B79-EC42-BDDB-88FD33DDBFB7.root",
      //base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__14C80F56-47FA-CF4A-82AF-ABBA59EA226B.root"},"1"));
  full_procs.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+#gamma", 
      Process::Type::background, TColor::GetColor("#16bac5"),
      {base_folder_ul+"mc/ZGToLLG_01J_5f_TuneCP5*"},"1"));
      //{base_folder_ul+"mc/ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__230000__03E636EE-BAB0-CA42-BFB1-28E772F2D53E.root",
      //base_folder_ul+"mc/ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__230000__1A3F9CD1-6D5C-FA41-968D-F9277B54C454.root",
      //base_folder_ul+"mc/ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__230000__5C9381D1-E437-3C45-8CE4-B1C68F64D70B.root",
      //base_folder_ul+"mc/ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__230000__61113508-9C90-544E-AAD5-2FF5D2C021E7.root",
      //base_folder_ul+"mc/ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__230000__6DAC4396-64CB-5E4C-8674-F89516B35910.root",
      //base_folder_ul+"mc/ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__230000__7A2FE5F0-8F95-784A-8E70-E9487FA641A8.root"},"1"));
  full_procs.push_back(Process::MakeShared<Baby_nano>("gg#rightarrow H#rightarrow ZG (x100)", 
      Process::Type::signal, TColor::GetColor("#ff0000"),
      {base_folder_ul+"signal/GluGluHToZG*"},"1"));

  std::vector<std::shared_ptr<Process>> full_procs_noscale;
  full_procs_noscale.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+jets", 
      Process::Type::background, TColor::GetColor("#ffb400"),
      {base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__0*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__1*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__2*"},"1"));
  full_procs_noscale.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+#gamma", 
      Process::Type::background, TColor::GetColor("#16bac5"),
      {base_folder_ul+"mc/ZGToLLG_01J_5f_TuneCP5*"},"1"));
  full_procs_noscale.push_back(Process::MakeShared<Baby_nano>("gg#rightarrow H#rightarrow ZG", 
      Process::Type::signal, TColor::GetColor("#ff0000"),
      {base_folder_ul+"signal/GluGluHToZG*"},"1"));

  std::vector<std::shared_ptr<Process>> procs_sideband_2018;
  procs_sideband_2018.push_back(Process::MakeShared<Baby_nano>("DYJets", 
      Process::Type::background, TColor::GetColor("#ffb400"),
      {base_folder_v9+"2018/mc/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL18NanoAODv9__106X_upgrade2018_realistic_v16_L1v1-v2__230000__0*",
      base_folder_v9+"2018/mc/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL18NanoAODv9__106X_upgrade2018_realistic_v16_L1v1-v2__230000__1*",
      base_folder_v9+"2018/mc/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL18NanoAODv9__106X_upgrade2018_realistic_v16_L1v1-v2__230000__2*"},"1"));
  procs_sideband_2018.push_back(Process::MakeShared<Baby_nano>("ZGToLLG", 
      Process::Type::background, TColor::GetColor("#16bac5"),
      {base_folder_v9+"2018/mc/ZGToLLG_01J*"},"1"));
  procs_sideband_2018.push_back(Process::MakeShared<Baby_nano>("Data (EGamma)", 
      Process::Type::data, kBlack,
      {base_folder_v9+"2018/data/EGamma*"},"1"));

  std::vector<std::shared_ptr<Process>> procs_data_sideband;
  //for Nanos, have to do overlap removal by hand
  procs_data_sideband.push_back(Process::MakeShared<Baby_nano>("DoubleMuon", 
      Process::Type::data, kBlack,
      {base_folder_v9+"2016/data/DoubleMuon*",
      base_folder_v9+"2016APV/data/DoubleMuon*",
      base_folder_v9+"2017/data/DoubleMuon*",
      base_folder_v9+"2018/data/DoubleMuon*"},pass_json&&blind_sr));
  procs_data_sideband.push_back(Process::MakeShared<Baby_nano>("EGamma", 
      Process::Type::data, kBlack,
      {base_folder_v9+"2018/data/EGamma*"},pass_json&&blind_sr&&!HLT_pass_dimuon));
  procs_data_sideband.push_back(Process::MakeShared<Baby_nano>("DoubleEG", 
      Process::Type::data, kBlack,
      {base_folder_v9+"2016/data/DoubleEG*",
      base_folder_v9+"2016APV/data/DoubleEG*",
      base_folder_v9+"2017/data/DoubleEG*"},pass_json&&blind_sr&&!HLT_pass_dimuon));
  procs_data_sideband.push_back(Process::MakeShared<Baby_nano>("SingleElectron", 
      Process::Type::data, kBlack,
      {base_folder_v9+"2016/data/SingleElectron*",
      base_folder_v9+"2016APV/data/SingleElectron*",
      base_folder_v9+"2017/data/SingleElectron*"},pass_json&&blind_sr&&!HLT_pass_dilepton&&!HLT_pass_diphoton));
  procs_data_sideband.push_back(Process::MakeShared<Baby_nano>("SingleMuon", 
      Process::Type::data, kBlack,
      {base_folder_v9+"2016/data/SingleMuon*",
      base_folder_v9+"2016APV/data/SingleMuon*",
      base_folder_v9+"2017/data/SingleMuon*",
      base_folder_v9+"2018/data/SingleMuon*"},pass_json&&blind_sr&&!HLT_pass_dilepton&&!HLT_pass_diphoton&&!HLT_pass_singleelectron));
  procs_data_sideband.push_back(Process::MakeShared<Baby_nano>("MuonEG", 
      Process::Type::data, kBlack,
      {base_folder_v9+"2016/data/MuonEG*",
      base_folder_v9+"2016APV/data/MuonEG*",
      base_folder_v9+"2017/data/MuonEG*",
      base_folder_v9+"2018/data/MuonEG*"},pass_json&&blind_sr&&!HLT_pass_dilepton&&!HLT_pass_diphoton&&!HLT_pass_singlelepton));

  std::vector<std::shared_ptr<Process>> dy_procs;
  dy_procs.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+jets", 
      Process::Type::background, TColor::GetColor("#ffb400"),
      {base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__049988F3-37E5-0346-ABFA-252C0CC3F235.root",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__04D4C0DC-7E25-784E-B9F9-4FCBD0D02E4A.root",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__08224897-3E25-064E-A735-5ED118DB12CB.root", //temp?
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__0979AFCD-1562-D540-9CF3-513D158A14ED.root"},"1"));
      //{base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__035CF0D8-DDC5-6D46-B172-3BDC98CC624D*",
      //base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__223A6F1C-C8C0-9F43-94C3-233AD51ABD32*",
      //base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__29EF0133-B16A-E747-9E02-89BC6682B669*", //temp?
      //base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__25118058-0F30-8A41-B82B-C7EAFEE55379*", //temp?
      //base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__24950148-3722-B84D-9050-C996F5205F0A*"},"1"));

  std::vector<std::shared_ptr<Process>> dy_procs_2016;
  dy_procs_2016.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+jets", 
      Process::Type::background, TColor::GetColor("#ffb400"),
      {base_folder_ulv9+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL16NanoAODv9__106X_mcRun2_asymptotic_v17-v1__30000__0082C29D-E74C-024A-BE9B-97B29EE7A4A2*",
      base_folder_ulv9+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL16NanoAODv9__106X_mcRun2_asymptotic_v17-v1__30000__0718C107-8960-6B44-B96A-C60D53D52A95*"},"1"));

  std::vector<std::shared_ptr<Process>> signal_trigger_procs;
  signal_trigger_procs.push_back(Process::MakeShared<Baby_nano>("gg#rightarrow H#rightarrow ZG (Dilepton Triggers)", 
      Process::Type::signal, TColor::GetColor("#ff0000"),
      {base_folder_ul+"signal/GluGluHToZG*"},hlt_el_trigger||hlt_mu_trigger));
  signal_trigger_procs.push_back(Process::MakeShared<Baby_nano>("gg#rightarrow H#rightarrow ZG (Only Single Lepton Triggers)", 
      Process::Type::signal, TColor::GetColor("#0000ff"),
      {base_folder_ul+"signal/GluGluHToZG*"},(hlt_single_el_trigger||hlt_single_mu_trigger)&&!(hlt_el_trigger||hlt_mu_trigger)));

  std::vector<std::shared_ptr<Process>> signal_diph_trigger_procs;
  signal_diph_trigger_procs.push_back(Process::MakeShared<Baby_nano>("gg#rightarrow H#rightarrow ZG (Single + Dilepton Triggers)", 
      Process::Type::signal, TColor::GetColor("#ff0000"),
      {base_folder_ul+"signal/GluGluHToZG*"},hlt_single_el_trigger||hlt_single_mu_trigger||hlt_el_trigger||hlt_mu_trigger));
  signal_diph_trigger_procs.push_back(Process::MakeShared<Baby_nano>("gg#rightarrow H#rightarrow ZG (Added by Diphoton Triggers)", 
      Process::Type::signal, TColor::GetColor("#0000ff"),
      {base_folder_ul+"signal/GluGluHToZG*"},(hlt_doubleph_trigger)&&!(hlt_single_el_trigger||hlt_single_mu_trigger||hlt_el_trigger||hlt_mu_trigger)));
  
  std::vector<std::shared_ptr<Process>> signal_procs;
  signal_procs.push_back(Process::MakeShared<Baby_nano>("gg#rightarrow H#rightarrow ZG (x100)", 
      Process::Type::background, TColor::GetColor("#ff0000"),
      {base_folder_ul+"signal/GluGluHToZG*"},"1"));

  std::vector<std::shared_ptr<Process>> signal_procs_noscale;
  signal_procs_noscale.push_back(Process::MakeShared<Baby_nano>("gg#rightarrow H#rightarrow ZG", 
      Process::Type::signal, TColor::GetColor("#ff0000"),
      {base_folder_ul+"signal/GluGluHToZG*"},"1"));

  std::vector<std::shared_ptr<Process>> signal_procs_2018;
  signal_procs_2018.push_back(Process::MakeShared<Baby_nano>("gg#rightarrow H#rightarrow ZG (2018)", 
      Process::Type::signal, TColor::GetColor("#ff0000"),
      {base_folder_v7+"2018/signal/GluGluHToZG*"},"1"));

  std::vector<std::shared_ptr<Process>> signal_procs_noscale_green;
  signal_procs_noscale_green.push_back(Process::MakeShared<Baby_nano>("gg#rightarrow H#rightarrow ZG", 
      Process::Type::signal, TColor::GetColor("#009933"),
      {base_folder_ul+"signal/GluGluHToZG*"},"1"));

  std::vector<std::shared_ptr<Process>> signal_procs_noscale_2d;
  signal_procs_noscale_2d.push_back(Process::MakeShared<Baby_nano>("gg#rightarrow H#rightarrow ZG", 
      Process::Type::background, TColor::GetColor("#ff0000"),
      {base_folder_ul+"signal/GluGluHToZG*"},"1"));

  std::vector<std::shared_ptr<Process>> signal_procs_vbf;
  signal_procs_vbf.push_back(Process::MakeShared<Baby_nano>("VBF H#rightarrow Z#gamma", 
      Process::Type::signal, TColor::GetColor("#ff0000"),
      {base_folder_ul+"signal/VBFHToZG*"},"1"));

  std::vector<std::shared_ptr<Process>> signal_procs_wminush;
  signal_procs_wminush.push_back(Process::MakeShared<Baby_nano>("W^{-}H#rightarrow Z#gamma", 
      Process::Type::signal, TColor::GetColor("#ff0000"),
      {base_folder_ul+"signal/WminusH*"},"1"));

  std::vector<std::shared_ptr<Process>> signal_procs_wplush;
  signal_procs_wplush.push_back(Process::MakeShared<Baby_nano>("W^{+}H#rightarrow Z#gamma", 
      Process::Type::signal, TColor::GetColor("#ff0000"),
      {base_folder_ul+"signal/WplusH*"},"1"));

  std::vector<std::shared_ptr<Process>> signal_procs_zh;
  signal_procs_zh.push_back(Process::MakeShared<Baby_nano>("ZH#rightarrow Z#gamma", 
      Process::Type::signal, TColor::GetColor("#ff0000"),
      {base_folder_ul+"signal/ZH*"},"1"));

  std::vector<std::shared_ptr<Process>> signal_procs_tth;
  signal_procs_tth.push_back(Process::MakeShared<Baby_nano>("ttH#rightarrow Z#gamma", 
      Process::Type::signal, TColor::GetColor("#ff0000"),
      {base_folder_ul+"signal/ttHToZG*"},"1"));

  std::vector<std::shared_ptr<Process>> procs_onelep_data;
  procs_onelep_data.push_back(Process::MakeShared<Baby_nano>("Data", 
      Process::Type::data, kBlack,
      //{base_folder_v7data+"SingleElectron__Run2017B*",
      //base_folder_v7data+"SingleMuon__Run2017B*"},"1"));
      {
      base_folder_ulv9_2017+"/data/SingleElectron__Run2017B*",
      base_folder_ulv9_2017+"/data/SingleMuon__Run2017B*",
      //base_folder_ulv9_2017+"/data/SingleElectron__Run2017B__UL2017_MiniAODv2_NanoAODv9-v1__120000__20867FA7-22D4-CC46-BE86-EE9FA1C846B6.root",
      //base_folder_ulv9_2017+"/data/SingleElectron__Run2017B__UL2017_MiniAODv2_NanoAODv9-v1__120000__20BA2BC5-9F81-1449-98D6-A4DEF081BD4B.root",
      //base_folder_ulv9_2017+"/data/SingleElectron__Run2017B__UL2017_MiniAODv2_NanoAODv9-v1__120000__2C497F12-38D9-CE43-BB34-A85309F9DE9B.root",
      //base_folder_ulv9_2017+"/data/SingleElectron__Run2017B__UL2017_MiniAODv2_NanoAODv9-v1__120000__46E53FF3-D096-C647-83A1-8112BC83D056.root",
      //base_folder_ulv9_2017+"/data/SingleMuon__Run2017B__UL2017_MiniAODv2_NanoAODv9-v1__120000__09FD9FD6-A164-9A45-80BB-F3D1FBF9C462.root",
      //base_folder_ulv9_2017+"/data/SingleMuon__Run2017B__UL2017_MiniAODv2_NanoAODv9-v1__120000__107D61C7-DD24-5848-956B-633AE6499A92.root",
      //base_folder_ulv9_2017+"/data/SingleMuon__Run2017B__UL2017_MiniAODv2_NanoAODv9-v1__120000__12CE0481-77A2-FD48-9FBD-43BE811EAD5D.root",
      //base_folder_ulv9_2017+"/data/SingleMuon__Run2017B__UL2017_MiniAODv2_NanoAODv9-v1__120000__136F7D48-764A-6D49-9AF6-03C79A625F65.root",
      //base_folder_v7data+"SingleElectron__Run2017B__02Apr2020-v1__230000__014FD4FF-A181-2843-B548-BB3198F16E88*",
      //base_folder_v7data+"SingleElectron__Run2017B__02Apr2020-v1__230000__08B5DC81-D780-A54B-80FE-C94EBD267ACA*",
      //base_folder_v7data+"SingleElectron__Run2017B__02Apr2020-v1__230000__0ACA0341-E365-2544-99F0-FBA4C92C0301*",
      //base_folder_v7data+"SingleElectron__Run2017B__02Apr2020-v1__230000__0C083B57-67B2-FF43-AF2E-56C850A094D6*",
      //base_folder_v7data+"SingleMuon__Run2017B__02Apr2020-v1__230000__0C07D411-3FE0-974E-A5A4-25BB763E4038*",
      //base_folder_v7data+"SingleMuon__Run2017B__02Apr2020-v1__230000__0D1DA690-909B-F74F-8DB7-48D6F786F8D8*",
      //base_folder_v7data+"SingleMuon__Run2017B__02Apr2020-v1__230000__12EF443F-9CD8-7842-9563-2BBD627A28FB*",
      //base_folder_v7data+"SingleMuon__Run2017B__02Apr2020-v1__230000__1A14183C-9C11-E947-9D71-6B89994127D9*"
      },"1"));

  std::vector<std::shared_ptr<Process>> procs_run2017d;
  procs_run2017d.push_back(Process::MakeShared<Baby_nano>("Data", 
      Process::Type::data, kBlack,
      {base_folder_v7data+"DoubleEG__Run2017D*",
      base_folder_v7data+"DoubleMuon__Run2017D*",
      base_folder_v7data+"SingleElectron__Run2017D*",
      base_folder_v7data+"SingleMuon__Run2017D*"},"1"));
      
  std::vector<std::shared_ptr<Process>> procs_met_data;
  procs_met_data.push_back(Process::MakeShared<Baby_nano>("Data", 
      Process::Type::data, kBlack,
      {base_folder_v7data+"MET__Run2017B__02Apr2020-v1__230000__00DCCA4E-F5F1-F84D-A6EC-2956ACAB6E02*",
      base_folder_v7data+"MET__Run2017B__02Apr2020-v1__230000__0B4DB0F0-168B-544A-B5F9-A5B3ACFE7F1E*",
      base_folder_v7data+"MET__Run2017B__02Apr2020-v1__230000__0F13F71C-9A13-5641-A5C1-4955D60DE44A*",
      base_folder_v7data+"MET__Run2017B__02Apr2020-v1__230000__1781BB70-1AD7-2F49-8BF6-A6DFE7C58167*",
      base_folder_v7data+"MET__Run2017F__02Apr2020-v1__30000__034FA7B4-A09F-0D47-B752-2897BE4FE43F*",
      base_folder_v7data+"MET__Run2017F__02Apr2020-v1__30000__0D3E9423-8FB8-9548-B63E-B94913999718*",
      base_folder_v7data+"MET__Run2017F__02Apr2020-v1__30000__0EE78B5A-1800-0341-A33E-29EFDCAB3FAA*",
      base_folder_v7data+"MET__Run2017F__02Apr2020-v1__30000__1450C885-647B-074B-BAEB-84763C291B92*"},"1"));

  vector<shared_ptr<Process>> procs_sig_pinnacles = ZgSampleLoader()
      .SetMacro("SKIM_LLG",{"mc/merged_zgmc_llg"})
      .LoadSamples("txt/samples_zgamma.txt","All",{Process::Type::signal});
  procs_sig_pinnacles[0]->cut_ = "1";

  vector<shared_ptr<Process>> procs_sig_pinnacles_unskimmed = ZgSampleLoader()
      .SetMacro("SKIM_LLG",{"mc/unskimmed"})
      .LoadSamples("txt/samples_zgamma.txt","All",{Process::Type::signal});
  procs_sig_pinnacles_unskimmed[0]->cut_ = "1";

  vector<shared_ptr<Process>> procs_sig_pinnacles_compare = ZgSampleLoader()
      .SetMacro("SKIM_LLG",{"mc/merged_zgmc_llg"})
      .LoadSamples("txt/samples_zgamma.txt","All",{Process::Type::signal});
  procs_sig_pinnacles_compare[0]->cut_ = "trig_double_el||trig_double_mu";
  procs_sig_pinnacles_compare[0]->name_ = "H#rightarrowZ#gamma (Dilepton triggers)";
  procs_sig_pinnacles_compare[0]->SetLineColor(TColor::GetColor(63,144,218));
  vector<shared_ptr<Process>> procs_sig_pinnacles_compare2 = ZgSampleLoader()
      .SetMacro("SKIM_LLG",{"mc/merged_zgmc_llg"})
      .LoadSamples("txt/samples_zgamma.txt","All",{Process::Type::signal});
  procs_sig_pinnacles_compare2[0]->cut_ = 
      "trig_double_el||trig_double_mu||trig_single_el||trig_single_mu";
  procs_sig_pinnacles_compare2[0]->name_ = 
      "H#rightarrowZ#gamma (Dilepton or single lepton triggers)";
  //procs_sig_pinnacles_compare2[0]->color_ = TColor::GetColor(189,31,1);
  procs_sig_pinnacles_compare.insert( procs_sig_pinnacles_compare.end(), 
      procs_sig_pinnacles_compare2.begin(),procs_sig_pinnacles_compare2.end());

  vector<shared_ptr<Process>> procs_pinnacles_mc = ZgSampleLoader()
      .SetMacro("SKIM_LLG",{"mc/merged_zgmc_llg"})
      .LoadMCSamples("txt/samples_zgamma.txt","CondensedMC");
  procs_pinnacles_mc.erase(procs_pinnacles_mc.begin()+2);
  for (shared_ptr<Process> &proc : procs_pinnacles_mc) {
    proc->cut_ = "use_event";
    if (proc->type_ == Process::Type::signal)
      proc->name_ = proc->name_+" (x100)";
  }

  //                      named funcs to be moved to a common location
  //---------------------------------------------------------------------------

  const NamedFunc NMuon_trig = ReduceNamedFunc("TrigObj_id==13",reduce_sum);
  const NamedFunc NElectron_trig = ReduceNamedFunc("TrigObj_id==11",
                                                   reduce_sum);
  const NamedFunc NPhoton_trig = ReduceNamedFunc("TrigObj_id==22",reduce_sum);

  const NamedFunc Lead_Electron_abseta("Lead_Electron_abseta",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_sig_[iel]) {
        return abs(b.Electron_eta()->at(iel));
      }
    }
    return -999;
  });

  const NamedFunc Sublead_Electron_abseta("Sublead_Electron_abseta",[](const Baby &b) -> NamedFunc::ScalarType{
    int sig_idx = 0;
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_sig_[iel]) {
        if (sig_idx==1) {
          return abs(b.Electron_eta()->at(iel));
        }
        sig_idx++;
      }
    }
    return -999;
  });

  const NamedFunc Lead_Muon_abseta("Lead_Muon_abseta",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      if (Muon_sig_[imu]) {
        return abs(b.Muon_eta()->at(imu));
      }
    }
    return -999;
  });

  const NamedFunc Sublead_Muon_abseta("Sublead_Muon_abseta",[](const Baby &b) -> NamedFunc::ScalarType{
    int sig_idx = 0;
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      if (Muon_sig_[imu]) {
        if (sig_idx==1)
          return abs(b.Muon_eta()->at(imu));
        sig_idx++;
      }
    }
    return -999;
  });

  //---------------------------------------------------------------------------
  //                                     named funcs
  //---------------------------------------------------------------------------
  
  const NamedFunc weight("weight_nf",[](const Baby &b) -> NamedFunc::ScalarType{
    double h_zg_br = 0.001533;
    double z_ll_br = 0.100974;
    //extra factor for weighting is xs/nevts_effective where effective events are calculated including negative weights
    //0.007519850838359999 GluGluH->ZG->llG xs (pb)
    //967054
    //double xs = 7.7760403e-09; //GluGluHToZG
    double xs = 0.00751985;
    double nents = 400000.0/100.0; //artificial x100
    if (b.FirstFileName().find("VBFHToZG") != std::string::npos) {
      xs = h_zg_br*z_ll_br*3.782;
      nents = 200000;
    }
    if (b.FirstFileName().find("WminusH_HToZG") != std::string::npos) {
      xs = h_zg_br*0.527;
      nents = 299276;
    }
    if (b.FirstFileName().find("WplusH_HToZG") != std::string::npos) {
      xs = h_zg_br*0.831;
      nents = 299978;
    }
    if (b.FirstFileName().find("ZH_HToZG") != std::string::npos) {
      xs = h_zg_br*0.8839;
      nents = 297389;
    }
    if (b.FirstFileName().find("ttHToZG") != std::string::npos) {
      xs = h_zg_br*0.5071;
      nents = 200000;
    }
    if (b.FirstFileName().find("ZGToLLG") != std::string::npos) {
      xs = 117.864;
      //nents = 2.8095231e+09;
      //nents = 1240092;
      //nents = 9315642;
      nents = 29885702;
      if (b.SampleType()==2018) nents = 18750664; //nents(genweight>0)-nents(genweight<0), all samples
    }
    if (b.FirstFileName().find("DYJetsToLL") != std::string::npos) {
      xs = 6077.22;
      //nents = 3.2011619e+11;
      //nents = 3535863;
      //nents = 7680238;
      nents = 19014526;
      if (b.SampleType()==2016) nents = 3865459; 
      if (b.SampleType()==2018) nents = 10052001; //nents(genweight>0)-nents(genweight<0), 230000__0* through 2*
    }
    if (b.FirstFileName().find("EGamma") != std::string::npos) return 1.0; //data unweighted
    float w_lumi = 1.0;
    float w_year = 1.0;
    if (b.Generator_weight()<0) w_lumi = -1.0;
    if (b.SampleType()==2016) w_year = 36.32264; 
    else if (b.SampleType()==2017) w_year = 41.52756;
    else if (b.SampleType()==2018) w_year = 59.67377;
    return w_lumi*w_year*xs*1000.0/nents;
  });

  const NamedFunc weight_noscale("weight_noscale",[](const Baby &b) -> NamedFunc::ScalarType{
    double h_zg_br = 0.001533;
    double z_ll_br = 0.100974;
    //extra factor for weighting is xs/nevts_effective where effective events are calculated including negative weights
    //0.007519850838359999 GluGluH->ZG->llG xs (pb)
    //967054
    //double xs = 7.7760403e-09; //GluGluHToZG
    double xs = 0.00751985;
    double nents = 400000.0;
    if (b.FirstFileName().find("VBFHToZG") != std::string::npos) {
      xs = h_zg_br*z_ll_br*3.782;
      nents = 200000;
    }
    if (b.FirstFileName().find("WminusH_HToZG") != std::string::npos) {
      xs = h_zg_br*0.527;
      nents = 299276;
    }
    if (b.FirstFileName().find("WplusH_HToZG") != std::string::npos) {
      xs = h_zg_br*0.831;
      nents = 299978;
    }
    if (b.FirstFileName().find("ZH_HToZG") != std::string::npos) {
      xs = h_zg_br*0.8839;
      nents = 297389;
    }
    if (b.FirstFileName().find("ttHToZG") != std::string::npos) {
      xs = h_zg_br*0.5071;
      nents = 200000;
    }
    if (b.FirstFileName().find("ZGToLLG") != std::string::npos) {
      xs = 117.864;
      //nents = 2.8095231e+09;
      nents = 1240092;
    }
    if (b.FirstFileName().find("DYJetsToLL") != std::string::npos) {
      xs = 6077.22;
      //nents = 3.2011619e+11;
      //nents = 3535863;
      //nents = 7680238;
      nents = 5056384;
    }
    //float w_lumi = 1.0;
    float w_year = 1.0;
    //if (b.Generator_weight()<0) w_lumi = -1.0;
    if (b.SampleType()==2016) w_year = 36.32264; 
    else if (b.SampleType()==2017) w_year = 41.52756;
    else if (b.SampleType()==2018) w_year = 59.67377;
    return w_year*xs*1000.0/nents;
  });

  const NamedFunc weight_old("weight_old",[](const Baby &b) -> NamedFunc::ScalarType{
    //extra factor for weighting is xs/nevts_effective where effective events are calculated including negative weights
    //0.007519850838359999 GluGluH->ZG->llG xs (pb)
    //967054
    float w_lumi = 1.0;
    float w_year = 1.0;
    if (b.Generator_weight()<0) w_lumi = -1.0;
    if (b.SampleType()==2016) w_year = 36.32264; 
    else if (b.SampleType()==2017) w_year = 41.52756;
    else if (b.SampleType()==2018) w_year = 59.67377;
    float corr_factor = 112.0/104.0; //temp, just to match scale with AN
    return w_lumi*w_year*7.7760403e-09*1000.0*corr_factor;
  });

  const NamedFunc z_decay_pdgid("z_decay_pdgid",[](const Baby &b) -> NamedFunc::ScalarType{
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid != 21 && abs_pdgid != 22 && abs_pdgid != 23) { //skip radiated particles and Z
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==23) {
            int z_mom_idx = b.GenPart_genPartIdxMother()->at(mom_idx);
            while (true) {
              if (z_mom_idx != -1) {
                int z_mom_id = b.GenPart_pdgId()->at(z_mom_idx);
                if (z_mom_id == 23) {
                  z_mom_idx = b.GenPart_genPartIdxMother()->at(z_mom_idx);
                }
                else {
                  if (z_mom_id == 25) { //Z is from H
                    return abs_pdgid;
                  }
                  break;
                }
              }
              else {
                break;
              }
            } //looking for Z parent
          }
        }
      }
    }
    return 0;
  });

  const NamedFunc nh_z_decay_pdgid("nh_z_decay_pdgid",[](const Baby &b) -> NamedFunc::ScalarType{
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid != 21 && abs_pdgid != 22 && abs_pdgid != 23) { //skip radiated particles and Z
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==23) {
            int z_mom_idx = b.GenPart_genPartIdxMother()->at(mom_idx);
            while (true) {
              if (z_mom_idx != -1) {
                int z_mom_id = b.GenPart_pdgId()->at(z_mom_idx);
                if (z_mom_id == 23) {
                  z_mom_idx = b.GenPart_genPartIdxMother()->at(z_mom_idx);
                }
                else {
                  if (z_mom_id != 25) { //not H
                    return abs_pdgid;
                  }
                  break;
                }
              }
              else {
                break;
              }
            } //looking for Z parent
          }
        }
      }
    }
    return 0;
  });

  //number of electrons out of eta acceptance
  const NamedFunc nel_ooeta("nel_ooeta",[](const Baby &b) -> NamedFunc::ScalarType{
    int nel_ooa_ = 0;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==11) {
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==23) {
            if (abs(b.GenPart_eta()->at(imc))>2.5)
              nel_ooa_ += 1;
          }
        }
      }
    }
    return nel_ooa_;
  });

  //true if electrons out of level 1 pt acceptance
  const NamedFunc el_ool1pt("el_ool1pt",[](const Baby &b) -> NamedFunc::ScalarType{
    float el1_pt = 0;
    float el2_pt = 0;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==11) {
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==23) {
            if (b.GenPart_pt()->at(imc) > 24.0) return false;
            if (el1_pt <= 0)
              el1_pt = b.GenPart_pt()->at(imc);
            else
              el2_pt = b.GenPart_pt()->at(imc);
          }
        }
      }
    }
    if ((el1_pt > 18 && el2_pt > 17) || (el1_pt > 17 && el2_pt > 18)) return false;
    if ((el1_pt > 22 && el2_pt > 10) || (el1_pt > 10 && el2_pt > 22)) return false;
    return true;
  });

  //true if muons out of level 1 pt acceptance
  const NamedFunc mu_ool1pt("mu_ool1pt",[](const Baby &b) -> NamedFunc::ScalarType{
    float mu1_pt = 0;
    float mu2_pt = 0;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==13) {
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==23) {
            if (mu1_pt <= 0)
              mu1_pt = b.GenPart_pt()->at(imc);
            else
              mu2_pt = b.GenPart_pt()->at(imc);
          }
        }
      }
    }
    if ((mu1_pt > 12 && mu2_pt > 5) || (mu1_pt > 5 && mu2_pt > 12)) return false;
    return true;
  });

  //true if electrons out of HLT pt acceptance
  const NamedFunc el_oohltpt("el_oohltpt",[](const Baby &b) -> NamedFunc::ScalarType{
    float el1_pt = 0;
    float el2_pt = 0;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==11) {
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==23) {
            if (el1_pt <= 0)
              el1_pt = b.GenPart_pt()->at(imc);
            else
              el2_pt = b.GenPart_pt()->at(imc);
          }
        }
      }
    }
    if ((el1_pt > 23 && el2_pt > 12) || (el1_pt > 12 && el2_pt > 23)) return false;
    return true;
  });

  //true if electrons nearly out of HLT pt acceptance (2 GeV)
  const NamedFunc el_noohltpt("el_noohltpt",[](const Baby &b) -> NamedFunc::ScalarType{
    float el1_pt = 0;
    float el2_pt = 0;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==11) {
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==23) {
            if (el1_pt <= 0)
              el1_pt = b.GenPart_pt()->at(imc);
            else
              el2_pt = b.GenPart_pt()->at(imc);
          }
        }
      }
    }
    if ((el1_pt > 25 && el2_pt > 14) || (el1_pt > 14 && el2_pt > 25)) return false;
    return true;
  });

  const NamedFunc el_hlt_fail_reason("el_hlt_fail_reason",[](const Baby &b) -> NamedFunc::ScalarType{
    int n_pf_el = 0;
    std::vector<float> hlt_el_pt;
    for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
      if (b.TrigObj_id()->at(itrig) == 11) {
        n_pf_el++;
        if ((b.TrigObj_filterBits()->at(itrig) & 0x1)==1) { //CaloIdL_TrackIdL_IsoVL
          hlt_el_pt.push_back(b.TrigObj_pt()->at(itrig));
        }
      }
    }
    if (n_pf_el < 2) return 1; //electron not reconstructed at HLT
    if (hlt_el_pt.size() < 2) return 2; //electron failed CaloIdL_TrackIdL_IsoVL
    std::sort(hlt_el_pt.begin(), hlt_el_pt.end());
    if (hlt_el_pt[hlt_el_pt.size()-1] < 23.0) return 3; //top leg failed pt cut
    if (hlt_el_pt[hlt_el_pt.size()-2] < 12.0) return 4; //bottom leg failed pt cut
    return 0; //unknown
  });

  const NamedFunc nel_id_hlt("nel_id_hlt",[](const Baby &b) -> NamedFunc::ScalarType{
    int nel_hlt = 0;
    for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
      if (b.TrigObj_id()->at(itrig) == 11) {
        if ((b.TrigObj_filterBits()->at(itrig) & 0x1)==1) { //CaloIdL_TrackIdL_IsoVL
          nel_hlt++;
        }
      }
    }
    return nel_hlt;
  });

  const NamedFunc nel_noid_hlt("nel_noid_hlt",[](const Baby &b) -> NamedFunc::ScalarType{
    int nel_hlt = 0;
    for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
      if (b.TrigObj_id()->at(itrig) == 11) {
        nel_hlt++;
      }
    }
    return nel_hlt;
  });

  const NamedFunc hlt_max_el_pt("hlt_max_el_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    float hlt_max_el_pt_ = 0;
    for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
      if (b.TrigObj_id()->at(itrig) == 11) {
        if ((b.TrigObj_filterBits()->at(itrig) & 0x1)==1) { //CaloIdL_TrackIdL_IsoVL
          if (b.TrigObj_pt()->at(itrig) > hlt_max_el_pt_)
            hlt_max_el_pt_ = b.TrigObj_pt()->at(itrig);
        }
      }
    }
    return hlt_max_el_pt_;
  });

  const NamedFunc HLT_Ele27_CaloIdL_TrackIdL_IsoVL_Photon30("HLT_Ele27_CaloIdL_TrackIdL_IsoVL_Photon30",[](const Baby &b) -> NamedFunc::ScalarType{
    bool found_good_el = false;
    bool found_good_ph = false;
    for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
      if (b.TrigObj_id()->at(itrig) == 11) {
        if ((b.TrigObj_filterBits()->at(itrig) & 0x1)==1) { //CaloIdL_TrackIdL_IsoVL
          if (b.TrigObj_pt()->at(itrig) > 27) {
            found_good_el = true;
          }
        }
      }
      if (b.TrigObj_id()->at(itrig) == 22) {
        if (b.TrigObj_pt()->at(itrig) > 30) {
          found_good_ph = true;
        }
      }
    }
    return (found_good_el && found_good_ph);
  });

  const NamedFunc HLT_Ele27_CaloIdL_TrackIdL_IsoVL_Photon50("HLT_Ele27_CaloIdL_TrackIdL_IsoVL_Photon50",[](const Baby &b) -> NamedFunc::ScalarType{
    bool found_good_el = false;
    for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
      if (b.TrigObj_id()->at(itrig) == 11) {
        if ((b.TrigObj_filterBits()->at(itrig) & 0x1)==1) { //CaloIdL_TrackIdL_IsoVL
          if (b.TrigObj_pt()->at(itrig) > 27) {
            found_good_el = true;
          }
        }
      }
    }
    return found_good_el && b.HLT_Photon50();
  });

  const NamedFunc HLT_Ele27_CaloIdL_TrackIdL_IsoVL_Photon33("HLT_Ele27_CaloIdL_TrackIdL_IsoVL_Photon33",[](const Baby &b) -> NamedFunc::ScalarType{
    bool found_good_el = false;
    for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
      if (b.TrigObj_id()->at(itrig) == 11) {
        if ((b.TrigObj_filterBits()->at(itrig) & 0x1)==1) { //CaloIdL_TrackIdL_IsoVL
          if (b.TrigObj_pt()->at(itrig) > 27) {
            found_good_el = true;
          }
        }
      }
    }
    return found_good_el && b.HLT_Photon33();
  });

  const NamedFunc HLT_Ele27_CaloIdL_TrackIdL_IsoVL_Photon25("HLT_Ele27_CaloIdL_TrackIdL_IsoVL_Photon25",[](const Baby &b) -> NamedFunc::ScalarType{
    bool found_good_el = false;
    for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
      if (b.TrigObj_id()->at(itrig) == 11) {
        if ((b.TrigObj_filterBits()->at(itrig) & 0x1)==1) { //CaloIdL_TrackIdL_IsoVL
          if (b.TrigObj_pt()->at(itrig) > 27) {
            found_good_el = true;
          }
        }
      }
    }
    return found_good_el && b.HLT_Photon25();
  });
  
  const NamedFunc HLT_Ele27_CaloIdL_TrackIdL_IsoVL_Photon20_HoverELoose("HLT_Ele27_CaloIdL_TrackIdL_IsoVL_Photon20_HoverELoose",[](const Baby &b) -> NamedFunc::ScalarType{
    bool found_good_el = false;
    for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
      if (b.TrigObj_id()->at(itrig) == 11) {
        if ((b.TrigObj_filterBits()->at(itrig) & 0x1)==1) { //CaloIdL_TrackIdL_IsoVL
          if (b.TrigObj_pt()->at(itrig) > 27) {
            found_good_el = true;
          }
        }
      }
    }
    return found_good_el && b.HLT_Photon20_HoverELoose();
  });

  const NamedFunc HLT_Ele22_CaloIdL_TrackIdL_IsoVL_Photon50("HLT_Ele22_CaloIdL_TrackIdL_IsoVL_Photon50",[](const Baby &b) -> NamedFunc::ScalarType{
    bool found_good_el = false;
    for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
      if (b.TrigObj_id()->at(itrig) == 11) {
        if ((b.TrigObj_filterBits()->at(itrig) & 0x1)==1) { //CaloIdL_TrackIdL_IsoVL
          if (b.TrigObj_pt()->at(itrig) > 22) {
            found_good_el = true;
          }
        }
      }
    }
    return found_good_el && b.HLT_Photon50();
  });

  const NamedFunc HLT_Ele22_CaloIdL_TrackIdL_IsoVL_Photon33("HLT_Ele22_CaloIdL_TrackIdL_IsoVL_Photon33",[](const Baby &b) -> NamedFunc::ScalarType{
    bool found_good_el = false;
    for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
      if (b.TrigObj_id()->at(itrig) == 11) {
        if ((b.TrigObj_filterBits()->at(itrig) & 0x1)==1) { //CaloIdL_TrackIdL_IsoVL
          if (b.TrigObj_pt()->at(itrig) > 22) {
            found_good_el = true;
          }
        }
      }
    }
    return found_good_el && b.HLT_Photon33();
  });

  const NamedFunc HLT_Ele22_CaloIdL_TrackIdL_IsoVL_Photon25("HLT_Ele22_CaloIdL_TrackIdL_IsoVL_Photon25",[](const Baby &b) -> NamedFunc::ScalarType{
    bool found_good_el = false;
    for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
      if (b.TrigObj_id()->at(itrig) == 11) {
        if ((b.TrigObj_filterBits()->at(itrig) & 0x1)==1) { //CaloIdL_TrackIdL_IsoVL
          if (b.TrigObj_pt()->at(itrig) > 22) {
            found_good_el = true;
          }
        }
      }
    }
    return found_good_el && b.HLT_Photon25();
  });

  const NamedFunc HLT_Ele22_CaloIdL_TrackIdL_IsoVL_Photon20_HoverELoose("HLT_Ele22_CaloIdL_TrackIdL_IsoVL_Photon20_HoverELoose",[](const Baby &b) -> NamedFunc::ScalarType{
    bool found_good_el = false;
    for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
      if (b.TrigObj_id()->at(itrig) == 11) {
        if ((b.TrigObj_filterBits()->at(itrig) & 0x1)==1) { //CaloIdL_TrackIdL_IsoVL
          if (b.TrigObj_pt()->at(itrig) > 22) {
            found_good_el = true;
          }
        }
      }
    }
    return found_good_el && b.HLT_Photon20_HoverELoose();
  });

  const NamedFunc HLT_Ele17_Ele8_CaloIdL_TrackIdL_IsoVL_Photon20_HoverELoose("HLT_Ele17_Ele8_CaloIdL_TrackIdL_IsoVL_Photon20_HoverELoose",[](const Baby &b) -> NamedFunc::ScalarType{
    float max_el_pt = 0;
    float subl_el_pt = 0;
    for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
      if (b.TrigObj_id()->at(itrig) == 11) {
        if ((b.TrigObj_filterBits()->at(itrig) & 0x1)==1) { //CaloIdL_TrackIdL_IsoVL
          if (b.TrigObj_pt()->at(itrig) > max_el_pt) {
            subl_el_pt = max_el_pt;
            max_el_pt = b.TrigObj_pt()->at(itrig);
          }
          else if (b.TrigObj_pt()->at(itrig) > subl_el_pt) {
            subl_el_pt = b.TrigObj_pt()->at(itrig);
          }
        }
      }
    }
    return (max_el_pt > 17) && (subl_el_pt > 8) && b.HLT_Photon20_HoverELoose();
  });

  const NamedFunc HLT_Ele17_Ele8_CaloIdL_TrackIdL_IsoVL_Photon25("HLT_Ele17_Ele8_CaloIdL_TrackIdL_IsoVL_Photon25",[](const Baby &b) -> NamedFunc::ScalarType{
    float max_el_pt = 0;
    float subl_el_pt = 0;
    for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
      if (b.TrigObj_id()->at(itrig) == 11) {
        if ((b.TrigObj_filterBits()->at(itrig) & 0x1)==1) { //CaloIdL_TrackIdL_IsoVL
          if (b.TrigObj_pt()->at(itrig) > max_el_pt) {
            subl_el_pt = max_el_pt;
            max_el_pt = b.TrigObj_pt()->at(itrig);
          }
          else if (b.TrigObj_pt()->at(itrig) > subl_el_pt) {
            subl_el_pt = b.TrigObj_pt()->at(itrig);
          }
        }
      }
    }
    return (max_el_pt > 17) && (subl_el_pt > 8) && b.HLT_Photon25();
  });

  NamedFunc hlt_el_trigger_plusplus = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL||HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ||HLT_Ele35_WPTight_Gsf||HLT_Ele32_WPTight_Gsf_L1DoubleEG" || HLT_Ele27_CaloIdL_TrackIdL_IsoVL_Photon33;

  const NamedFunc mu_hlt_fail_reason("mu_hlt_fail_reason",[](const Baby &b) -> NamedFunc::ScalarType{
    int n_pf_mu = 0;
    int n_good_mu = 0;
    std::vector<float> hlt_mu_pt;
    for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
      if (b.TrigObj_id()->at(itrig) == 13) {
        n_pf_mu++;
        if ((b.TrigObj_filterBits()->at(itrig) & 0x1)==1) { //TkIsoVVL
          n_good_mu++;
        }
        hlt_mu_pt.push_back(b.TrigObj_pt()->at(itrig));
      }
    }
    if (n_pf_mu < 2) return 1; //muon not reconstructed at HLT
    //if (hlt_mu_pt.size() < 2) return 2; //muon fails Iso
    std::sort(hlt_mu_pt.begin(), hlt_mu_pt.end());
    if (hlt_mu_pt[hlt_mu_pt.size()-1] < 23.0) return 3; //top leg failed pt cut
    if (hlt_mu_pt[hlt_mu_pt.size()-2] < 12.0) return 4; //bottom leg failed pt cut
    if (n_good_mu < 2) return 2;
    return 0; //unknown
  });

  const NamedFunc nmu_id_hlt("nmu_id_hlt",[](const Baby &b) -> NamedFunc::ScalarType{
    int n_pf_mu = 0;
    for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
      if (b.TrigObj_id()->at(itrig) == 13) {
        if ((b.TrigObj_filterBits()->at(itrig) & 0x1)==1) { //TkIsoVVL
          n_pf_mu++;
        }
      }
    }
    return n_pf_mu;
  });

  const NamedFunc nmu_noid_hlt("nmu_noid_hlt",[](const Baby &b) -> NamedFunc::ScalarType{
    int n_pf_mu = 0;
    for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
      if (b.TrigObj_id()->at(itrig) == 13) {
        n_pf_mu++;
      }
    }
    return n_pf_mu;
  });

  //weights from data that should account for why things fail HLT
  const NamedFunc hlt_fail_weights("hlt_fail_weights",[](const Baby &b) -> NamedFunc::ScalarType{
    float el1_pt = 0;
    float el1_eta = 0;
    float el2_pt = 0;
    float el2_eta = 0;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==11) {
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==23) {
            if (el1_pt <= 0) {
              el1_pt = b.GenPart_pt()->at(imc);
              el1_eta = abs(b.GenPart_eta()->at(imc));
            }
            else {
              el2_pt = b.GenPart_pt()->at(imc);
              el2_eta = abs(b.GenPart_eta()->at(imc));
            }
          }
        }
      }
    }
    float leg1_w = 0.0, leg2_w = 0.0;
    if (el1_pt > 20 && el1_pt < 23 && el1_eta > 0.00 && el1_eta < 0.80) leg1_w = 0.25;
    else if (el1_pt > 20 && el1_pt < 23 && el1_eta > 0.80 && el1_eta < 1.44) leg1_w = 0.18;
    else if (el1_pt > 20 && el1_pt < 23 && el1_eta > 1.44 && el1_eta < 1.56) leg1_w = 0.27;
    else if (el1_pt > 20 && el1_pt < 23 && el1_eta > 1.56 && el1_eta < 2.50) leg1_w = 0.15;
    else if (el1_pt > 23 && el1_pt < 25 && el1_eta > 0.00 && el1_eta < 0.80) leg1_w = 0.80;
    else if (el1_pt > 23 && el1_pt < 25 && el1_eta > 0.80 && el1_eta < 1.44) leg1_w = 0.75;
    else if (el1_pt > 23 && el1_pt < 25 && el1_eta > 1.44 && el1_eta < 1.56) leg1_w = 0.65;
    else if (el1_pt > 23 && el1_pt < 25 && el1_eta > 1.56 && el1_eta < 2.50) leg1_w = 0.68;
    else if (el1_pt > 25 && el1_pt < 30 && !(el1_eta > 1.44 && el1_eta < 1.56)) leg1_w = 0.83;
    else if (el1_pt > 25 && el1_pt < 30 && el1_eta > 1.44 && el1_eta < 1.56) leg1_w = 0.68;
    else if (el1_pt > 30 && el1_pt < 40 && !(el1_eta > 1.44 && el1_eta < 1.56)) leg1_w = 0.90;
    else if (el1_pt > 30 && el1_pt < 40 && el1_eta > 1.44 && el1_eta < 1.56) leg1_w = 0.78;
    else if (el1_pt > 40 && el1_pt < 50 && !(el1_eta > 1.44 && el1_eta < 1.56)) leg1_w = 0.92;
    else if (el1_pt > 40 && el1_pt < 50 && el1_eta > 1.44 && el1_eta < 1.56) leg1_w = 0.83;
    else if (el1_pt > 50 && el1_pt < 100 && !(el1_eta > 1.44 && el1_eta < 1.56)) leg1_w = 0.95;
    else if (el1_pt > 50 && el1_pt < 100 && el1_eta > 1.44 && el1_eta < 1.56) leg1_w = 0.83;
    else if (el1_pt > 100 && !(el1_eta > 1.44 && el1_eta < 1.56)) leg1_w = 0.97;
    else if (el1_pt > 100 && el1_eta > 1.44 && el1_eta < 1.56) leg1_w = 0.91;
    if (el2_pt > 12 && el2_pt < 14 && el2_eta > 0.00 && el2_eta < 0.80) leg2_w = 0.0;
    else if (el2_pt > 12 && el2_pt < 14 && el2_eta > 0.80 && el2_eta < 1.44) leg2_w = 0.0;
    else if (el2_pt > 12 && el2_pt < 14 && el2_eta > 1.44 && el2_eta < 1.56) leg2_w = 0.0;
    else if (el2_pt > 12 && el2_pt < 14 && el2_eta > 1.56 && el2_eta < 2.50) leg2_w = 0.38;
    else if (el2_pt > 14 && el2_pt < 17 && el2_eta > 0.00 && el2_eta < 0.80) leg2_w = 0.78;
    else if (el2_pt > 14 && el2_pt < 17 && el2_eta > 0.80 && el2_eta < 1.44) leg2_w = 0.70;
    else if (el2_pt > 14 && el2_pt < 17 && el2_eta > 1.44 && el2_eta < 1.56) leg2_w = 0.52;
    else if (el2_pt > 14 && el2_pt < 17 && el2_eta > 1.56 && el2_eta < 2.50) leg2_w = 0.57;
    else if (el2_pt > 17 && el2_pt < 20 && el2_eta > 0.00 && el2_eta < 0.80) leg2_w = 0.80;
    else if (el2_pt > 17 && el2_pt < 20 && el2_eta > 0.80 && el2_eta < 1.44) leg2_w = 0.75;
    else if (el2_pt > 17 && el2_pt < 20 && el2_eta > 1.44 && el2_eta < 1.56) leg2_w = 0.54;
    else if (el2_pt > 17 && el2_pt < 20 && el2_eta > 1.56 && el2_eta < 2.50) leg2_w = 0.70;
    else if (el2_pt > 20 && el2_pt < 25 && !(el2_eta > 1.44 && el2_eta < 1.56)) leg2_w = 0.80;
    else if (el2_pt > 20 && el2_pt < 25 && el2_eta > 1.44 && el2_eta < 1.56) leg2_w = 0.65;
    else if (el2_pt > 25 && el2_pt < 27 && !(el2_eta > 1.44 && el2_eta < 1.56)) leg2_w = 0.82;
    else if (el2_pt > 25 && el2_pt < 27 && el2_eta > 1.44 && el2_eta < 1.56) leg2_w = 0.70;
    else if (el2_pt > 27 && el2_pt < 30 && !(el2_eta > 1.44 && el2_eta < 1.56)) leg2_w = 0.85;
    else if (el2_pt > 27 && el2_pt < 30 && el2_eta > 1.44 && el2_eta < 1.56) leg2_w = 0.73;
    else if (el2_pt > 30 && el2_pt < 40 && !(el2_eta > 1.44 && el2_eta < 1.56)) leg2_w = 0.88;
    else if (el2_pt > 30 && el2_pt < 40 && el2_eta > 1.44 && el2_eta < 1.56) leg2_w = 0.76;
    else if (el2_pt > 40 && el2_pt < 60 && !(el2_eta > 1.44 && el2_eta < 1.56)) leg2_w = 0.92;
    else if (el2_pt > 40 && el2_pt < 60 && el2_eta > 1.44 && el2_eta < 1.56) leg2_w = 0.82;
    else if (el2_pt > 60 && !(el2_eta > 1.44 && el2_eta < 1.56)) leg2_w = 0.97;
    else if (el2_pt > 60 && el2_eta > 1.44 && el2_eta < 1.56) leg2_w = 0.91;
    return (1.0-leg1_w*leg2_w);
  });

  //true if muons out of hlt pt acceptance
  const NamedFunc mu_oohltpt("mu_oohltpt",[](const Baby &b) -> NamedFunc::ScalarType{
    float mu1_pt = 0;
    float mu2_pt = 0;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==13) {
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==23) {
            if (mu1_pt <= 0)
              mu1_pt = b.GenPart_pt()->at(imc);
            else
              mu2_pt = b.GenPart_pt()->at(imc);
          }
        }
      }
    }
    if ((mu1_pt > 17 && mu2_pt > 8) || (mu1_pt > 8 && mu2_pt > 17)) return false;
    return true;
  });

  const NamedFunc Electron_hltIndex("Electron_hltIndex",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> el_matchhlt_;
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    double matchhlt = -1.0;
    double dr = 999.0;
    double mindr = 99.0;
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      matchhlt = -1.0;
      mindr = 99.0;
      if (Electron_sig_[iel] > 0.5) {
        //try to match to HLT object
        for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
          if (b.TrigObj_id()->at(itrig) == 11) {
            //if ((b.TrigObj_filterBits()->at(itrig) & 0x1)==1) { //CaloIdL_TrackIdL_IsoVL
            dr = deltaR(b.Electron_eta()->at(iel), b.Electron_phi()->at(iel), b.TrigObj_eta()->at(itrig), b.TrigObj_phi()->at(itrig));
            if (dr < 0.4 && dr < mindr) {
              matchhlt = static_cast<double>(itrig);
              mindr = dr;
            }
            //}
          }
        }
      }
      el_matchhlt_.push_back(matchhlt);
    }
    return el_matchhlt_;
  });

  const NamedFunc Electron_hltId("Electron_hltId",[Electron_hltIndex](const Baby &b) -> NamedFunc::VectorType{
    //-1 - no HLT object, 0 - fails all, 1 - pass CaloIdL_TrackIdL_IsoVL, 2 - pass WPTight
    std::vector<double> el_hltid;
    std::vector<double> hlt_idx = Electron_hltIndex.GetVector(b);
    double id;
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      id = -1;
      if (hlt_idx[iel]>=0) {
        id = 0;
        int itrig = static_cast<int>(hlt_idx[iel]);
        if ((b.TrigObj_filterBits()->at(itrig) & 0x2) == 0x2)
          id = 2;
        else if ((b.TrigObj_filterBits()->at(itrig) & 0x1) == 0x1)
          id = 1;
      }
      el_hltid.push_back(id);
    }
    return el_hltid;
  });

  const NamedFunc Muon_hltIndex("Muon_hltIndex",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> mu_matchhlt_;
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    double matchhlt = -1.0;
    double dr = 999.0;
    double mindr = 99.0;
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      matchhlt = -1.0;
      mindr = 99.0;
      if (Muon_sig_[imu] > 0.5) {
        //try to match to HLT object
        for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
          if (b.TrigObj_id()->at(itrig) == 13) {
            //if ((b.TrigObj_filterBits()->at(itrig) & 0x2)==1) { //Iso, TrkIsoVVL is broken at HLT
            dr = deltaR(b.Muon_eta()->at(imu), b.Muon_phi()->at(imu), b.TrigObj_eta()->at(itrig), b.TrigObj_phi()->at(itrig));
            if (dr < 0.1 && dr < mindr) {
              matchhlt = static_cast<double>(itrig);
              mindr = dr;
            }
            //}
          }
        }
      }
      mu_matchhlt_.push_back(matchhlt);
    }
    return mu_matchhlt_;
  });

  const NamedFunc Muon_hltId("Muon_hltId",[Muon_hltIndex](const Baby &b) -> NamedFunc::VectorType{
    //-1 - no HLT object, 0 - fails all, 1 - pass TrkIsoVVL, 2 - Iso
    std::vector<double> mu_hltid;
    std::vector<double> hlt_idx = Muon_hltIndex.GetVector(b);
    double id;
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      id = -1;
      if (hlt_idx[imu]>=0) {
        id = 0;
        int itrig = static_cast<int>(hlt_idx[imu]);
        if ((b.TrigObj_filterBits()->at(itrig) & 0x2) == 0x2)
          id = 2;
        else if ((b.TrigObj_filterBits()->at(itrig) & 0x1) == 0x1)
          id = 1;
      }
      mu_hltid.push_back(id);
    }
    return mu_hltid;
  });

  const NamedFunc el_gap("el_gap",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_sig_[iel] > 0.5) {
        float etasc = abs(b.Electron_deltaEtaSC()->at(iel) + b.Electron_eta()->at(iel));
        float eta = abs(b.Electron_eta()->at(iel));
        if (eta > 1.44 && eta < 1.56) return true;
        if (etasc > 1.44 && etasc < 1.56) return true;
      }
    }
    return false;
  });

  const NamedFunc Lead_Photon_mvaID("Lead_Photon_mvaID",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> Photon_sig_ = Photon_sig.GetVector(b);
    for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
      if (Photon_sig_[iph]) {
        return b.Photon_mvaID()->at(iph);
      }
    }
    return 0;
  });

  const NamedFunc LeadMuon_absEta("LeadMuon_absEta",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      if (Muon_sig_[imu]) {
        return fabs(b.Muon_eta()->at(imu));
      }
    }
    return 0;
  });

  const NamedFunc TrigObj_filterBit2("TrigObj_filterBit2",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> filterbit2;
    for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
      filterbit2.push_back(static_cast<double>((b.TrigObj_filterBits()->at(itrig) & 0x02) >> 1));
    }
    return filterbit2;
  });

  const NamedFunc Electron_IsoPhotonPt("Electron_IsoPhotonPt",[](const Baby &b) -> NamedFunc::VectorType{
    //-1 - no HLT object, 0 - fails all, 1 - pass CaloIdL_TrackIdL_IsoVL, 2 - pass WPTight
    std::vector<double> el_isophpt;
    std::vector<double> ph_sig_ = Photon_sig.GetVector(b);
    double isophpt = 0;
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      isophpt = 0;
      for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
        if (ph_sig_[iph]>0.5) {
          if (deltaR(b.Electron_eta()->at(iel), b.Electron_phi()->at(iel), b.Photon_eta()->at(iph), b.Photon_phi()->at(iph)) < 0.3) {
            isophpt += b.Photon_pt()->at(iph);
          }
        }
      }
      el_isophpt.push_back(isophpt);
    }
    return el_isophpt;
  });

  const NamedFunc Electron_IsoElectronPt("Electron_IsoElectronPt",[](const Baby &b) -> NamedFunc::VectorType{
    //-1 - no HLT object, 0 - fails all, 1 - pass CaloIdL_TrackIdL_IsoVL, 2 - pass WPTight
    std::vector<double> el_isoelpt;
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    double isoelpt = 0;
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      isoelpt = 0;
      for (unsigned iel2 = 0; iel2 < b.Electron_pt()->size(); iel2++) {
        if (iel == iel2) continue;
        if (Electron_sig_[iel2]>0.5) {
          if (deltaR(b.Electron_eta()->at(iel), b.Electron_phi()->at(iel), b.Electron_eta()->at(iel2), b.Electron_phi()->at(iel2)) < 0.3) {
            isoelpt += b.Electron_pt()->at(iel2);
          }
        }
      }
      el_isoelpt.push_back(isoelpt);
    }
    return el_isoelpt;
  });

  const NamedFunc Electron_IsoJetPt("Electron_IsoJetPt",[](const Baby &b) -> NamedFunc::VectorType{
    //-1 - no HLT object, 0 - fails all, 1 - pass CaloIdL_TrackIdL_IsoVL, 2 - pass WPTight
    std::vector<double> el_isojetpt;
    std::vector<double> jet_sig_ = Jet_sig.GetVector(b);
    double isojetpt = 0;
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      isojetpt = 0;
      for (unsigned ijet = 0; ijet < b.Jet_pt()->size(); ijet++) {
        if (jet_sig_[ijet]>0.5) {
          if (deltaR(b.Electron_eta()->at(iel), b.Electron_phi()->at(iel), b.Jet_eta()->at(ijet), b.Jet_phi()->at(ijet)) < 0.3) {
            isojetpt += b.Jet_pt()->at(ijet);
          }
        }
      }
      el_isojetpt.push_back(isojetpt);
    }
    return el_isojetpt;
  });

  const NamedFunc stitch_dy("stitch_dy",[](const Baby &b) -> NamedFunc::ScalarType{
    //extra factor for weighting is xs/nevts_effective where effective events are calculated including negative weights
    //0.007519850838359999 GluGluH->ZG->llG xs (pb)
    //967054
    std::vector<double> ph_sig_ = Photon_sig.GetVector(b);
    if (b.FirstFileName().find("DYJetsToLL") != std::string::npos) {
      for (int imc = 0; imc < b.nGenPart(); imc++) {
        if (b.GenPart_pdgId()->at(imc)==22 && b.GenPart_status()->at(imc)==1) {
          bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
          if (mc_statusFlags[0]||mc_statusFlags[8]) {
            for (int iph = 0; iph < b.nPhoton(); iph++) {
              if (ph_sig_[iph]) {
                if (deltaR(b.GenPart_eta()->at(imc), b.GenPart_phi()->at(imc), b.Photon_eta()->at(iph), b.Photon_phi()->at(iph)) < 0.1) { 
                  return 0;
                }
              }
            }
          }
        }
      }
    }
    return 1;
  });


  const NamedFunc HiggsCandidate_ptOverMass("HiggsCandidate_ptOverMass",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    std::vector<double> ph_sig_ = Photon_sig.GetVector(b);
    ROOT::Math::PtEtaPhiMVector ll_p;
    double m_ll = -999;
    double pt_llg = 0;
    double m_llg = -999;
    for (unsigned iel = 1; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_sig_[iel]) {
        for (unsigned iel2 = 0; iel2 < iel; iel2++) {
          if (Electron_sig_[iel2]) {
            ROOT::Math::PtEtaPhiMVector p1(b.Electron_pt()->at(iel),b.Electron_eta()->at(iel),b.Electron_phi()->at(iel),0.000511);
            ROOT::Math::PtEtaPhiMVector p2(b.Electron_pt()->at(iel2),b.Electron_eta()->at(iel2),b.Electron_phi()->at(iel2),0.000511);
            double this_mll = (p1+p2).M();
            if (fabs(this_mll-91.2)<fabs(m_ll-91.2)) {
              ll_p = p1+p2;
              m_ll = this_mll;
            }
          }
        }
      }
    }
    for (unsigned imu = 1; imu < b.Muon_pt()->size(); imu++) {
      if (Muon_sig_[imu]) {
        for (unsigned imu2 = 0; imu2 < imu; imu2++) {
          if (Muon_sig_[imu2]) {
            ROOT::Math::PtEtaPhiMVector p1(b.Muon_pt()->at(imu),b.Muon_eta()->at(imu),b.Muon_phi()->at(imu),0.106);
            ROOT::Math::PtEtaPhiMVector p2(b.Muon_pt()->at(imu2),b.Muon_eta()->at(imu2),b.Muon_phi()->at(imu2),0.106);
            double this_mll = (p1+p2).M();
            if (fabs(this_mll-91.2)<fabs(m_ll-91.2)) {
              ll_p = p1+p2; 
              m_ll = this_mll;
            }
          }
        }
      }
    }
    if (m_ll > 0) {
      for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
        if (ph_sig_[iph]) {
          ROOT::Math::PtEtaPhiMVector pph(b.Photon_pt()->at(iph),b.Photon_eta()->at(iph),b.Photon_phi()->at(iph),0.0);
          double this_mass = (pph+ll_p).M();
          if (fabs(this_mass-125.3)<fabs(m_llg-125.3)) {
            pt_llg = (pph+ll_p).Pt();
            m_llg = this_mass;
          }
        }
      }
    }
    return pt_llg/m_llg;
  });

  const NamedFunc HiggsCandidate_mass_cleaned("HiggsCandidate_mass_cleaned",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    std::vector<double> ph_sig_ = Photon_sig.GetVector(b);
    ROOT::Math::PtEtaPhiMVector ll_p;
    double m_ll = -999;
    double m_llg = -999;
    for (unsigned iel = 1; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_sig_[iel]) {
        for (unsigned iel2 = 0; iel2 < iel; iel2++) {
          if (Electron_sig_[iel2]) {
            ROOT::Math::PtEtaPhiMVector p1(b.Electron_pt()->at(iel),b.Electron_eta()->at(iel),b.Electron_phi()->at(iel),0.000511);
            ROOT::Math::PtEtaPhiMVector p2(b.Electron_pt()->at(iel2),b.Electron_eta()->at(iel2),b.Electron_phi()->at(iel2),0.000511);
            double this_mll = (p1+p2).M();
            if (fabs(this_mll-91.2)<fabs(m_ll-91.2)) {
              ll_p = p1+p2;
              m_ll = this_mll;
            }
          }
        }
      }
    }
    for (unsigned imu = 1; imu < b.Muon_pt()->size(); imu++) {
      if (Muon_sig_[imu]) {
        for (unsigned imu2 = 0; imu2 < imu; imu2++) {
          if (Muon_sig_[imu2]) {
            ROOT::Math::PtEtaPhiMVector p1(b.Muon_pt()->at(imu),b.Muon_eta()->at(imu),b.Muon_phi()->at(imu),0.106);
            ROOT::Math::PtEtaPhiMVector p2(b.Muon_pt()->at(imu2),b.Muon_eta()->at(imu2),b.Muon_phi()->at(imu2),0.106);
            double this_mll = (p1+p2).M();
            if (fabs(this_mll-91.2)<fabs(m_ll-91.2)) {
              ll_p = p1+p2; 
              m_ll = this_mll;
            }
          }
        }
      }
    }
    if (m_ll > 0) {
      for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
        if (ph_sig_[iph] && b.Photon_cleanmask()->at(iph)==1) {
          ROOT::Math::PtEtaPhiMVector pph(b.Photon_pt()->at(iph),b.Photon_eta()->at(iph),b.Photon_phi()->at(iph),0.0);
          double this_mass = (pph+ll_p).M();
          if (fabs(this_mass-125.3)<fabs(m_llg-125.3))
            m_llg = this_mass;
        }
      }
    }
    return m_llg;
  });

  const NamedFunc delta_r_cut("delta_r_cut",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    std::vector<double> ph_sig_ = Photon_sig.GetVector(b);
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_sig_[iel]) {
        for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
          if (ph_sig_[iph]) {
            double this_deltar = deltaR(b.Electron_eta()->at(iel), b.Electron_phi()->at(iel), b.Photon_eta()->at(iph), b.Photon_phi()->at(iph));
            if (this_deltar < 0.4)
              return 0;
          }
        }
      }
    }
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      if (Muon_sig_[imu]) {
        for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
          if (ph_sig_[iph]) {
            double this_deltar = deltaR(b.Muon_eta()->at(imu), b.Muon_phi()->at(imu), b.Photon_eta()->at(iph), b.Photon_phi()->at(iph));
            if (this_deltar < 0.4)
              return 0;
          }
        }
      }
    }
    return 1;
  });

  const NamedFunc Min_Lepton_Photon_deltaR("Min_Lepton_Photon_deltaR",[](const Baby &b) -> NamedFunc::ScalarType{
    double min_deltar = 999;
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    std::vector<double> ph_sig_ = Photon_sig.GetVector(b);
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_sig_[iel]) {
        for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
          if (ph_sig_[iph]) {
            double this_deltar = deltaR(b.Electron_eta()->at(iel), b.Electron_phi()->at(iel), b.Photon_eta()->at(iph), b.Photon_phi()->at(iph));
            if (this_deltar < min_deltar)
              min_deltar = this_deltar;
          }
        }
      }
    }
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      if (Muon_sig_[imu]) {
        for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
          if (ph_sig_[iph]) {
            double this_deltar = deltaR(b.Muon_eta()->at(imu), b.Muon_phi()->at(imu), b.Photon_eta()->at(iph), b.Photon_phi()->at(iph));
            if (this_deltar < min_deltar)
              min_deltar = this_deltar;
          }
        }
      }
    }
    return min_deltar;
  });

  const NamedFunc ElectronGenPhotonDeltaR("ElectronGenPhotonDeltaR",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    double min_deltar = 999;
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_sig_[iel]) {
        for (int imc = 0; imc < b.nGenPart(); imc++) {
          if (b.GenPart_pdgId()->at(imc)==22 && b.GenPart_status()->at(imc)==1) {
            bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
            if (mc_statusFlags[0]||mc_statusFlags[8]) {
              if (b.GenPart_pt()->at(imc) > 20) {
                double this_deltar = deltaR(b.Electron_eta()->at(iel), b.Electron_phi()->at(iel), b.GenPart_eta()->at(imc), b.GenPart_phi()->at(imc));
                if (this_deltar < min_deltar)
                  min_deltar = this_deltar;
              }
            }
          }
        }
      }
    }
    return min_deltar;
  });

  //named funcs for efficiency measurements
  const NamedFunc Electron_absEta("Electron_absEta",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> abseta;
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      abseta.push_back(fabs(b.Electron_eta()->at(iel)));
    }
    return abseta;
  });

  const NamedFunc Muon_absEta("Muon_absEta",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> abseta;
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      abseta.push_back(fabs(b.Muon_eta()->at(imu)));
    }
    return abseta;
  });

  const NamedFunc LeadElectron_PassReference("LeadElectron_PassReference",[Electron_hltIndex](const Baby &b) -> NamedFunc::ScalarType{
    double leadel_passref = 0;
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    std::vector<double> hlt_idx = Electron_hltIndex.GetVector(b);
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_sig_[iel] > 0.5) {
        if (hlt_idx[iel]>=0) {
          if (((b.TrigObj_filterBits()->at(hlt_idx[iel]) & 0x2) == 0x2)
              && (b.TrigObj_pt()->at(hlt_idx[iel]) > 35)) {
            leadel_passref = 1;
          }
        }
        break;
      }
    }
    return leadel_passref;
  });

  const NamedFunc SubleadElectron_PassReference("SubleadElectron_PassReference",[Electron_hltIndex](const Baby &b) -> NamedFunc::ScalarType{
    double leadel_passref = 0;
    int Electron_sig_idx = 0;
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    std::vector<double> hlt_idx = Electron_hltIndex.GetVector(b);
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_sig_[iel] > 0.5) {
        if (Electron_sig_idx==1) {
          if (hlt_idx[iel]>=0) {
            if (((b.TrigObj_filterBits()->at(hlt_idx[iel]) & 0x2) == 0x2)
                && (b.TrigObj_pt()->at(hlt_idx[iel]) > 35)) {
              leadel_passref = 1;
            }
          }
          break;
        }
        Electron_sig_idx++;
      }
    }
    return leadel_passref;
  });

  const NamedFunc Electron_OtherPassReference("Electron_OtherPassReference",[LeadElectron_PassReference, SubleadElectron_PassReference](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> otherpassref;
    double ot_passref = 0;
    double leadpassref = LeadElectron_PassReference.GetScalar(b);
    double subleadpassref = SubleadElectron_PassReference.GetScalar(b);
    int Electron_sig_idx = 0;
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      ot_passref = leadpassref;
      if (Electron_sig_[iel] > 0.5) {
        if (Electron_sig_idx==0) {
          ot_passref = subleadpassref;
        }
        Electron_sig_idx++;
      }
      otherpassref.push_back(ot_passref);
    }
    return otherpassref;
  });

  const NamedFunc Electron_OtherPassUpperLeg("Electron_OtherPassUpperLeg",[Electron_hltIndex](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> otherpassupperleg;
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    std::vector<double> hlt_idx = Electron_hltIndex.GetVector(b);
    double this_otherpassupperleg = 0;
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      this_otherpassupperleg = 0;
      if (Electron_sig_[iel] > 0.5) {
        for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
          if (b.TrigObj_id()->at(itrig) == 11) {
            if (itrig != hlt_idx[iel]) {
              if ((((b.TrigObj_filterBits()->at(itrig) & 0x2) == 0x2)
                  || ((b.TrigObj_filterBits()->at(itrig) & 0x1) == 0x1))
                  && (b.TrigObj_pt()->at(itrig) > 23)) {
                this_otherpassupperleg = 1;
              } //TrigObj passes Ele23leg
            } //not matched TrigObj
          } //TrigObj is electron
        } //loop over TrigObj
      } //is signal electron
      otherpassupperleg.push_back(this_otherpassupperleg);
    }
    return otherpassupperleg;
  });

  const NamedFunc Electron_OtherPassLowerLeg("Electron_OtherPassLowerLeg",[Electron_hltIndex](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> otherpassupperleg;
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    std::vector<double> hlt_idx = Electron_hltIndex.GetVector(b);
    double this_otherpassupperleg = 0;
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      this_otherpassupperleg = 0;
      if (Electron_sig_[iel] > 0.5) {
        for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
          if (b.TrigObj_id()->at(itrig) == 11) {
            if (itrig != hlt_idx[iel]) {
              if ((((b.TrigObj_filterBits()->at(itrig) & 0x2) == 0x2)
                  || ((b.TrigObj_filterBits()->at(itrig) & 0x1) == 0x1))
                  && (b.TrigObj_pt()->at(itrig) > 12)) {
                this_otherpassupperleg = 1;
              } //TrigObj passes Ele23leg
            } //not matched TrigObj
          } //TrigObj is electron
        } //loop over TrigObj
      } //is signal electron
      otherpassupperleg.push_back(this_otherpassupperleg);
    }
    return otherpassupperleg;
  });

  const NamedFunc Electron_hltPt("Electron_hltPt",[Electron_hltIndex](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> hlt_pt;
    double this_hlt_pt;
    std::vector<double> hlt_idx = Electron_hltIndex.GetVector(b);
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      this_hlt_pt = -1;
      if (hlt_idx[iel]>=0) {
        this_hlt_pt = b.TrigObj_pt()->at(hlt_idx[iel]);
      }
      hlt_pt.push_back(this_hlt_pt);
    }
    return hlt_pt;
  });

  const NamedFunc LeadMuon_PassReference("LeadMuon_PassReference",[Muon_hltIndex](const Baby &b) -> NamedFunc::ScalarType{
    double leadmu_passref = 0;
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    std::vector<double> hlt_idx = Muon_hltIndex.GetVector(b);
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      if (Muon_sig_[imu] > 0.5) {
        if (hlt_idx[imu]>=0) {
          if (((b.TrigObj_filterBits()->at(hlt_idx[imu]) & 0x2) == 0x2)
              && (b.TrigObj_pt()->at(hlt_idx[imu]) > 27)) {
            leadmu_passref = 1;
          }
        }
        break;
      }
    }
    return leadmu_passref;
  });

  const NamedFunc SubleadMuon_PassReference("SubleadMuon_PassReference",[Muon_hltIndex](const Baby &b) -> NamedFunc::ScalarType{
    double leadmu_passref = 0;
    int Muon_sig_idx = 0;
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    std::vector<double> hlt_idx = Muon_hltIndex.GetVector(b);
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      if (Muon_sig_[imu] > 0.5) {
        if (Muon_sig_idx==1) {
          if (hlt_idx[imu]>=0) {
            if (((b.TrigObj_filterBits()->at(hlt_idx[imu]) & 0x2) == 0x2)
                && (b.TrigObj_pt()->at(hlt_idx[imu]) > 27)) {
              leadmu_passref = 1;
            }
          }
          break;
        }
        Muon_sig_idx++;
      }
    }
    return leadmu_passref;
  });

  const NamedFunc Muon_OtherPassReference("Muon_OtherPassReference",[LeadMuon_PassReference, SubleadMuon_PassReference](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> otherpassref;
    double ot_passref = 0;
    double leadpassref = LeadMuon_PassReference.GetScalar(b);
    double subleadpassref = SubleadMuon_PassReference.GetScalar(b);
    int Muon_sig_idx = 0;
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      ot_passref = leadpassref;
      if (Muon_sig_[imu] > 0.5) {
        if (Muon_sig_idx==0) {
          ot_passref = subleadpassref;
        }
        Muon_sig_idx++;
      }
      otherpassref.push_back(ot_passref);
    }
    return otherpassref;
  });

  const NamedFunc Muon_OtherPassUpperLeg("Muon_OtherPassUpperLeg",[Muon_hltIndex](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> otherpassupperleg;
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    std::vector<double> hlt_idx = Muon_hltIndex.GetVector(b);
    double this_otherpassupperleg = 0;
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      this_otherpassupperleg = 0;
      if (Muon_sig_[imu] > 0.5) {
        for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
          if (b.TrigObj_id()->at(itrig) == 13) {
            if (itrig != hlt_idx[imu]) {
              if ((((b.TrigObj_filterBits()->at(itrig) & 0x2) == 0x2)
                  || ((b.TrigObj_filterBits()->at(itrig) & 0x1) == 0x1))
                  && (b.TrigObj_pt()->at(itrig) > 17)) {
                this_otherpassupperleg = 1;
              } //TrigObj passes Mu17leg
            } //not matched TrigObj
          } //TrigObj is muon
        } //loop over TrigObj
      } //is signal muon
      otherpassupperleg.push_back(this_otherpassupperleg);
    }
    return otherpassupperleg;
  });


  const NamedFunc Muon_OtherPassLowerLeg("Muon_OtherPassLowerLeg",[Muon_hltIndex](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> otherpassupperleg;
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    std::vector<double> hlt_idx = Muon_hltIndex.GetVector(b);
    double this_otherpassupperleg = 0;
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      this_otherpassupperleg = 0;
      if (Muon_sig_[imu] > 0.5) {
        for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
          if (b.TrigObj_id()->at(itrig) == 13) {
            if (itrig != hlt_idx[imu]) {
              if ((((b.TrigObj_filterBits()->at(itrig) & 0x2) == 0x2)
                  || ((b.TrigObj_filterBits()->at(itrig) & 0x1) == 0x1))
                  && (b.TrigObj_pt()->at(itrig) > 8)) {
                this_otherpassupperleg = 1;
              } //TrigObj passes Mu17leg
            } //not matched TrigObj
          } //TrigObj is muon
        } //loop over TrigObj
      } //is signal muon
      otherpassupperleg.push_back(this_otherpassupperleg);
    }
    return otherpassupperleg;
  });

  const NamedFunc Electron_isLeading("Electron_isLeading",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> is_leading;
    double this_is_leading = 0;
    bool first_sig = true;
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      this_is_leading = 0;
      if (Electron_sig_[iel] > 0.5) {
        if (first_sig) {
          first_sig = false;
          this_is_leading = 1;
        }
      }
      is_leading.push_back(this_is_leading);
    }
    return is_leading;
  });

  const NamedFunc Muon_isLeading("Muon_isLeading",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> is_leading;
    double this_is_leading = 0;
    bool first_sig = true;
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      this_is_leading = 0;
      if (Muon_sig_[imu] > 0.5) {
        if (first_sig) {
          first_sig = false;
          this_is_leading = 1;
        }
      }
      is_leading.push_back(this_is_leading);
    }
    return is_leading;
  });

  const NamedFunc Muon_hltPt("Muon_hltPt",[Muon_hltIndex](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> hlt_pt;
    double this_hlt_pt;
    std::vector<double> hlt_idx = Muon_hltIndex.GetVector(b);
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      this_hlt_pt = -1;
      if (hlt_idx[imu]>=0) {
        this_hlt_pt = b.TrigObj_pt()->at(hlt_idx[imu]);
      }
      hlt_pt.push_back(this_hlt_pt);
    }
    return hlt_pt;
  });

  const NamedFunc eff_singlelep("eff_singlelep",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    std::vector<float> lepton_effs;
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_sig_[iel]) {
        lepton_effs.push_back(eff_isoel3235(b.Electron_pt()->at(iel), fabs(b.Electron_eta()->at(iel)))[0]);
      }
    }
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      if (Muon_sig_[imu]) {
        lepton_effs.push_back(eff_isomu2427(b.Muon_pt()->at(imu), fabs(b.Muon_eta()->at(imu)))[0]);
      }
    }
    float eff = 1.0;
    for (float lepton_eff : lepton_effs) {
      eff *= (1.0-lepton_eff);
    }
    eff = 1.0-eff;
    return eff;
  });

  const NamedFunc eff_dilep("eff_dilep",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    float eff_leadel = 0;
    float eff_leadmu = 0;
    float eff_sublel = 0;
    float eff_sublmu = 0;
    float lep_pt = 0;
    float lep_abseta = 0;
    int Electron_sig_idx = 0;
    int Muon_sig_idx = 0;
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_sig_[iel]) {
        lep_pt = b.Electron_pt()->at(iel);
        lep_abseta = fabs(b.Electron_eta()->at(iel));
        if (Electron_sig_idx == 0)
          eff_leadel = eff_dielleg12(lep_pt,lep_abseta)[0]*eff_elhltpt23(lep_pt,lep_abseta)[0];
        else if (Electron_sig_idx == 1)
          eff_sublel = eff_dielleg12(lep_pt,lep_abseta)[0];
        Electron_sig_idx++;
      }
    }
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      if (Muon_sig_[imu]) {
        lep_pt = b.Muon_pt()->at(imu);
        lep_abseta = fabs(b.Muon_eta()->at(imu));
        if (Muon_sig_idx == 0)
          eff_leadmu = eff_dimuleg8(lep_pt,lep_abseta)[0]*eff_muhltpt17(lep_pt,lep_abseta)[0];
        else if (Muon_sig_idx == 1)
          eff_sublmu = eff_dimuleg8(lep_pt,lep_abseta)[0];
        Muon_sig_idx++;
      }
    }
    if (isnan(eff_leadel)) eff_leadel = 0.0;
    if (isnan(eff_sublel)) eff_sublel = 0.0;
    if (isnan(eff_leadmu)) eff_leadmu = 0.0;
    if (isnan(eff_sublmu)) eff_sublmu = 0.0;
    return 1.0-(1.0-eff_leadel*eff_sublel)*(1.0-0.957*eff_leadmu*eff_sublmu);
  });

  const NamedFunc eff_singlelepgivendi("eff_singlelepgivendi",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    float eff_leadel = 0;
    float eff_leadmu = 0;
    float eff_sublel = 0;
    float eff_sublmu = 0;
    float eff_leadel_single = 0;
    float eff_sublel_single = 0;
    float eff_leadmu_single = 0;
    float eff_sublmu_single = 0;
    float lep_pt = 0;
    float lep_abseta = 0;
    int Electron_sig_idx = 0;
    int Muon_sig_idx = 0;
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_sig_[iel]) {
        lep_pt = b.Electron_pt()->at(iel);
        lep_abseta = fabs(b.Electron_eta()->at(iel));
        if (Electron_sig_idx == 0) {
          eff_leadel = eff_dielleg12(lep_pt,lep_abseta)[0]*eff_elhltpt23(lep_pt,lep_abseta)[0];
          eff_leadel_single = eff_isoel3235(lep_pt,lep_abseta)[0];
        }
        else if (Electron_sig_idx == 1) {
          eff_sublel = eff_dielleg12(lep_pt,lep_abseta)[0];
          eff_sublel_single = eff_isoel3235(lep_pt,lep_abseta)[0];
        }
        Electron_sig_idx++;
      }
    }
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      if (Muon_sig_[imu]) {
        lep_pt = b.Muon_pt()->at(imu);
        lep_abseta = fabs(b.Muon_eta()->at(imu));
        if (Muon_sig_idx == 0) {
          eff_leadmu = eff_dimuleg8(lep_pt,lep_abseta)[0]*eff_muhltpt17(lep_pt,lep_abseta)[0];
          eff_leadmu_single = eff_isomu2427(lep_pt,lep_abseta)[0];
        }
        else if (Muon_sig_idx == 1) {
          eff_sublmu = eff_dimuleg8(lep_pt,lep_abseta)[0];
          eff_sublmu_single = eff_isomu2427(lep_pt,lep_abseta)[0];
        }
        Muon_sig_idx++;
      }
    }
    if (isnan(eff_leadel) || eff_leadel<=0) { eff_leadel = 1.0; eff_leadel_single = 0.0; } 
    if (isnan(eff_sublel) || eff_sublel<=0) { eff_sublel = 1.0; eff_sublel_single = 0.0; } 
    if (isnan(eff_leadmu) || eff_leadmu<=0) { eff_leadmu = 1.0; eff_leadmu_single = 0.0; } 
    if (isnan(eff_sublmu) || eff_sublmu<=0) { eff_sublmu = 1.0; eff_sublel_single = 0.0; } 
    return 1.0-(1.0-eff_leadel_single/eff_leadel)*(1.0-eff_sublel_single/eff_sublel)*(1.0-eff_leadmu_single/eff_leadmu)*(1.0-eff_sublmu_single/eff_sublmu);
  });

  //somehow regress Higgs momentum
  //const NamedFunc Sum_SoftActivityJet_pz("Sum_SoftActivityJet_pz",[](const Baby &b) -> NamedFunc::ScalarType{
  //  ROOT::Math::PtEtaPhiMVector total_p;
  //  for (unsigned isj = 0; isj < b.SoftActivityJet_pt->size(); isj++) {
  //    ROOT::Math::PtEtaPhiMVector p;
  //    p.SetCoordinates(b.SoftActivityJet_pt->at(isj),b.SoftActivityJet_eta()->at(isj),b.SoftActivityJet_phi()->at(isj),0);
  //    total_p += p;
  //  }
  //});

  const NamedFunc eff_dilep_hig19014("eff_dilep_hig19014",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    float eff_leadel = 0;
    float eff_leadmu = 0;
    float eff_sublel = 0;
    float eff_sublmu = 0;
    float lep_pt = 0;
    float lep_abseta = 0;
    int Electron_sig_idx = 0;
    int Muon_sig_idx = 0;
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_sig_[iel]) {
        lep_pt = b.Electron_pt()->at(iel);
        lep_abseta = fabs(b.Electron_eta()->at(iel));
        if (Electron_sig_idx == 0)
          eff_leadel = eff_hig19014_ele23(lep_pt,lep_abseta);
        else if (Electron_sig_idx == 1)
          eff_sublel = eff_hig19014_ele12(lep_pt,lep_abseta);
        Electron_sig_idx++;
      }
    }
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      if (Muon_sig_[imu]) {
        lep_pt = b.Muon_pt()->at(imu);
        lep_abseta = fabs(b.Muon_eta()->at(imu));
        if (Muon_sig_idx == 0)
          eff_leadmu = eff_hig19014_mu17(lep_pt, lep_abseta);
        else if (Muon_sig_idx == 1)
          eff_sublmu = eff_hig19014_mu8(lep_pt,lep_abseta);
        Muon_sig_idx++;
      }
    }
    if (isnan(eff_leadel)) eff_leadel = 0.0;
    if (isnan(eff_sublel)) eff_sublel = 0.0;
    if (isnan(eff_leadmu)) eff_leadmu = 0.0;
    if (isnan(eff_sublmu)) eff_sublmu = 0.0;
    return 1.0-(1.0-eff_leadel*eff_sublel)*(1.0-0.957*eff_leadmu*eff_sublmu);
  });

  const NamedFunc eff_singleanddilep("eff_singleanddilep",[eff_singlelep, eff_dilep, eff_singlelepgivendi](const Baby &b) -> NamedFunc::ScalarType{
      return eff_singlelep.GetScalar(b)+eff_dilep.GetScalar(b)-eff_dilep.GetScalar(b)*eff_singlelepgivendi.GetScalar(b);
  });

  const NamedFunc LeptonPhoton_mass("LeptonPhoton_mass",[](const Baby &b) -> NamedFunc::VectorType{
    //get dilepton scale factors
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    std::vector<double> Photon_sig_ = Photon_sig.GetVector(b);
    std::vector<double> lepph_mass;
    ROOT::Math::PtEtaPhiMVector ph, el;
    for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
      if (Photon_sig_[iph]) {
        ph.SetCoordinates(b.Photon_pt()->at(iph),b.Photon_eta()->at(iph),b.Photon_phi()->at(iph),0);
        for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
          if (Electron_sig_[iel]) {
            el.SetCoordinates(b.Electron_pt()->at(iel),b.Electron_eta()->at(iel),b.Electron_phi()->at(iel),0.00511);
            lepph_mass.push_back((el+ph).M());
          }
        }
      }
    }
    return lepph_mass;
  });

  const NamedFunc Max_LeptonPhoton_mass("Max_LeptonPhoton_mass",[LeptonPhoton_mass](const Baby &b) -> NamedFunc::ScalarType{
    float max_lepph_mass = 0;
    //get dilepton scale factors
    std::vector<double> lepph_mass = LeptonPhoton_mass.GetVector(b);
    for (double this_mass : lepph_mass) {
      if (this_mass > max_lepph_mass)
        max_lepph_mass = this_mass;
    }
    return max_lepph_mass;
  });

  const NamedFunc ScaleFactor_Triggers_Old("ScaleFactor_Triggers_Old",[](const Baby &b) -> NamedFunc::ScalarType{
    float total_sf = 1;
    //get dilepton scale factors
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    float probdata_leadel = 0;
    float probdata_leadmu = 0;
    float probdata_sublel = 0;
    float probdata_sublmu = 0;
    float probsimu_leadel = 0;
    float probsimu_leadmu = 0;
    float probsimu_sublel = 0;
    float probsimu_sublmu = 0;
    float lep_pt = 0;
    float lep_abseta = 0;
    int Electron_sig_idx = 0;
    int Muon_sig_idx = 0;
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_sig_[iel]) {
        lep_pt = b.Electron_pt()->at(iel);
        lep_abseta = fabs(b.Electron_eta()->at(iel));
        if (Electron_sig_idx == 0) {
          //current data pt sfs unreliable, assume same between data and mc?
          //small by-hand correction, not sure why
          probsimu_leadel = effsig_dielleg12(lep_pt,lep_abseta)[0]*effsig_elhltpt23(lep_pt,lep_abseta)[0]/0.988;
          if (probsimu_leadel>1.0)
            probsimu_leadel = 1.0;
          probdata_leadel = eff_dielleg12(lep_pt,lep_abseta)[0]*effsig_elhltpt23(lep_pt,lep_abseta)[0];
        }
        else if (Electron_sig_idx == 1) {
          probsimu_sublel = effsig_dielleg12(lep_pt,lep_abseta)[0]/0.988;
          if (probsimu_sublel>1.0)
            probsimu_leadel = 1.0;
          probdata_sublel = eff_dielleg12(lep_pt,lep_abseta)[0];
        }
        Electron_sig_idx++;
      }
    }
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      if (Muon_sig_[imu]) {
        lep_pt = b.Muon_pt()->at(imu);
        lep_abseta = fabs(b.Muon_eta()->at(imu));
        if (Muon_sig_idx == 0) {
          //current data pt sfs unreliable, assume same between data and mc?
          probsimu_leadmu = effsig_dimuleg8(lep_pt,lep_abseta)[0]*effsig_muhltpt17(lep_pt,lep_abseta)[0]/0.998;
          if (probsimu_leadmu>1.0)
            probsimu_leadmu = 1.0;
          probdata_leadmu = eff_dimuleg8(lep_pt,lep_abseta)[0]*effsig_muhltpt17(lep_pt,lep_abseta)[0];
        }
        else if (Muon_sig_idx == 1) {
          probsimu_sublmu = effsig_dimuleg8(lep_pt,lep_abseta)[0]/0.998;
          if (probsimu_sublmu>1.0)
            probsimu_sublmu = 1.0;
          probdata_sublmu = eff_dimuleg8(lep_pt,lep_abseta)[0];
        }
        Muon_sig_idx++;
      }
    }
    probdata_leadel += 0;
    probdata_sublel += 0;
    probdata_leadmu += 0;
    probdata_sublmu += 0;
    if (Electron_sig_idx >= 2) { //if no signal electrons, prob~0 in both data and MC
      //shouldn't need to worry about divide by zeros since same MC was used to generate efficiencies
      if (b.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL()||b.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ()) {
        total_sf *= (probdata_leadel*probdata_sublel)/(probsimu_leadel*probsimu_sublel); 
      }
      else {
        total_sf *= (1.0-(probdata_leadel*probdata_sublel))/(1.0-(probsimu_leadel*probsimu_sublel));
      }
    }
    if (Muon_sig_idx >= 2) { //if no signal electrons, prob~0 in both data and MC
      if (b.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8()
          ||b.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8()) {
        //does simulation not simulate the DZ cut?
        total_sf *= (0.957*probdata_leadmu*probdata_sublmu)/(probsimu_leadmu*probsimu_sublmu);
      }
      else {
        total_sf *= (1.0-(0.957*probdata_leadmu*probdata_sublmu))/(1.0-(probsimu_leadmu*probsimu_sublmu));
      }
    }
    //get single lepton scale factors
    float prob_data = 0;
    float prob_simu = 0;
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_sig_[iel]) {
        lep_pt = b.Electron_pt()->at(iel);
        lep_abseta = fabs(b.Electron_eta()->at(iel));
        prob_simu = 1.0-(1.0-effsig_isoel3235(lep_pt,lep_abseta)[0])*(1.0-prob_simu);
        prob_data = 1.0-(1.0-eff_isoel3235(lep_pt,lep_abseta)[0])*(1.0-prob_data);
      }
    }
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      if (Muon_sig_[imu]) {
        lep_pt = b.Muon_pt()->at(imu);
        lep_abseta = fabs(b.Muon_eta()->at(imu));
        prob_simu = 1.0-(1.0-effsig_isomu2427(lep_pt,lep_abseta)[0])*(1.0-prob_simu);
        prob_data = 1.0-(1.0-eff_isomu2427(lep_pt,lep_abseta)[0])*(1.0-prob_data);
      }
    }
    if (Muon_sig_idx > 0 || Electron_sig_idx > 0) {
      if (b.HLT_Ele35_WPTight_Gsf()||b.HLT_Ele32_WPTight_Gsf_L1DoubleEG()
          ||b.HLT_IsoMu27()||b.HLT_IsoMu24()) {
        if (prob_simu > 0) {
          total_sf *= prob_data/prob_simu;
        }
      }
      else {
        if (prob_simu < 1) {
          total_sf *= (1.0-prob_data)/(1.0-prob_simu);
        }
      }
    }
    return total_sf;
  });

  const NamedFunc ScaleFactor_Triggers("ScaleFactor_Triggers",[](const Baby &b) -> NamedFunc::ScalarType{
    float total_prob_simu = 0;
    float total_prob_data = 0;
    //loop over leptons and construct table of probabilities
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    std::vector<int>    lep_flavor;
    std::vector<double> lep_simu_prob_single;
    std::vector<double> lep_simu_prob_dilep_upper;
    std::vector<double> lep_simu_prob_dilep_lower;
    std::vector<double> lep_simu_prob_fail_all;
    std::vector<double> lep_data_prob_single;
    std::vector<double> lep_data_prob_dilep_upper;
    std::vector<double> lep_data_prob_dilep_lower;
    std::vector<double> lep_data_prob_fail_all;
    float lep_pt = 0;
    float lep_abseta = 0;
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_sig_[iel]) {
        //current data pt sfs unreliable, assume same between data and mc?
        lep_pt = b.Electron_pt()->at(iel);
        lep_abseta = fabs(b.Electron_eta()->at(iel));
        lep_flavor.push_back(11);
        lep_simu_prob_single.push_back(fix_prob(effsig_isoel3235(lep_pt,lep_abseta)[0]));
        lep_simu_prob_dilep_upper.push_back(fix_prob(effsig_dielleg12(lep_pt,lep_abseta)[0]*effsig_elhltpt23(lep_pt,lep_abseta)[0]/0.988-lep_simu_prob_single.back()));
        lep_simu_prob_dilep_lower.push_back(fix_prob(effsig_dielleg12(lep_pt,lep_abseta)[0]*(1.0-effsig_elhltpt23(lep_pt,lep_abseta)[0])/0.988));
        lep_simu_prob_fail_all.push_back(1.0-lep_simu_prob_single.back()-lep_simu_prob_dilep_upper.back()-lep_simu_prob_dilep_lower.back());
        lep_data_prob_single.push_back(fix_prob(eff_isoel3235(lep_pt,lep_abseta)[0]));
        lep_data_prob_dilep_upper.push_back(fix_prob(eff_dielleg12(lep_pt,lep_abseta)[0]*effsig_elhltpt23(lep_pt,lep_abseta)[0]-lep_data_prob_single.back()));
        lep_data_prob_dilep_lower.push_back(fix_prob(eff_dielleg12(lep_pt,lep_abseta)[0]*(1.0-effsig_elhltpt23(lep_pt,lep_abseta)[0])));
        lep_data_prob_fail_all.push_back(1.0-lep_data_prob_single.back()-lep_data_prob_dilep_upper.back()-lep_data_prob_dilep_lower.back());
      }
    }
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      if (Muon_sig_[imu]) {
        //current data pt sfs unreliable, assume same between data and mc?
        lep_pt = b.Muon_pt()->at(imu);
        lep_abseta = fabs(b.Muon_eta()->at(imu));
        lep_flavor.push_back(13);
        lep_simu_prob_single.push_back(fix_prob(effsig_isomu2427(lep_pt,lep_abseta)[0]));
        lep_simu_prob_dilep_upper.push_back(fix_prob(effsig_dimuleg8(lep_pt,lep_abseta)[0]*effsig_muhltpt17(lep_pt,lep_abseta)[0]/0.988-lep_simu_prob_single.back()));
        lep_simu_prob_dilep_lower.push_back(fix_prob(effsig_dimuleg8(lep_pt,lep_abseta)[0]*(1.0-effsig_muhltpt17(lep_pt,lep_abseta)[0])/0.988));
        lep_simu_prob_fail_all.push_back(1.0-lep_simu_prob_single.back()-lep_simu_prob_dilep_upper.back()-lep_simu_prob_dilep_lower.back());
        lep_data_prob_single.push_back(fix_prob(eff_isomu2427(lep_pt,lep_abseta)[0]));
        lep_data_prob_dilep_upper.push_back(fix_prob(eff_dimuleg8(lep_pt,lep_abseta)[0]*effsig_muhltpt17(lep_pt,lep_abseta)[0]-lep_data_prob_single.back()));
        lep_data_prob_dilep_lower.push_back(fix_prob(eff_dimuleg8(lep_pt,lep_abseta)[0]*(1.0-effsig_muhltpt17(lep_pt,lep_abseta)[0])));
        lep_data_prob_fail_all.push_back(1.0-lep_data_prob_single.back()-lep_data_prob_dilep_upper.back()-lep_data_prob_dilep_lower.back());
      }
    }
    //loop over all possible reconstruction scenarios, of ones that trigger, add probability to total
    std::vector<int> combination_idx(lep_flavor.size(), 0); //0 - none, 1 - lower leg of dilepton, 2 - upper leg of dilepton, 3 - single lepton
    bool finished = false;
    while (!finished) {
      float prob_data = 1.0;
      float prob_simu = 1.0;
      bool pass_single = false;
      unsigned n_mu_upper = 0;
      unsigned n_mu_lower = 0;
      unsigned n_el_upper = 0;
      unsigned n_el_lower = 0;
      //loop over leptons
      for (unsigned lepton_idx = 0; lepton_idx < combination_idx.size(); lepton_idx++) {
        //calculate probability and see if we pass triggers for this case
        if (combination_idx[lepton_idx]==3) {
          prob_simu *= lep_simu_prob_single[lepton_idx];
          prob_data *= lep_data_prob_single[lepton_idx];
          pass_single = true;
        }
        else if (combination_idx[lepton_idx]==2) {
          prob_simu *= lep_simu_prob_dilep_upper[lepton_idx];
          prob_data *= lep_data_prob_dilep_upper[lepton_idx];
          if (lep_flavor[lepton_idx]==11) n_el_upper++;
          if (lep_flavor[lepton_idx]==13) n_mu_upper++;
        }
        else if (combination_idx[lepton_idx]==1) {
          prob_simu *= lep_simu_prob_dilep_lower[lepton_idx];
          prob_data *= lep_data_prob_dilep_lower[lepton_idx];
          if (lep_flavor[lepton_idx]==11) n_el_lower++;
          if (lep_flavor[lepton_idx]==13) n_mu_lower++;
        }
        else {
          prob_simu *= lep_simu_prob_fail_all[lepton_idx];
          prob_data *= lep_data_prob_fail_all[lepton_idx];
        }
      }
      //if we pass triggers, add probability to running total
      if (pass_single || (n_mu_upper>=2) || (n_mu_upper==1&&n_mu_lower>=1) || (n_el_upper>=2) || (n_el_upper==1&&n_el_lower>=1)) {
        total_prob_simu += prob_simu;
        total_prob_data += prob_data;
      }
      //next iteration go to next combination
      for (unsigned lepton_idx = 0; lepton_idx < combination_idx.size(); lepton_idx++) {
        if (combination_idx[lepton_idx] != 3) {
          combination_idx[lepton_idx] = combination_idx[lepton_idx]+1;
          break;
        }
        else if (lepton_idx == combination_idx.size()-1) {
          finished = true;
        }
        else {
          combination_idx[lepton_idx] = 0;
        }
      }
    }
    bool pass_triggers = b.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL()
        ||b.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ()||b.HLT_Ele35_WPTight_Gsf()
        ||b.HLT_Ele32_WPTight_Gsf_L1DoubleEG()
        ||b.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8()
        ||b.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8()||b.HLT_IsoMu27()
        ||b.HLT_IsoMu24();
    //return prob(data)/prob(simu) for appropriate case, avoiding divide-by-0
    if (!pass_triggers) {
      if (total_prob_simu>=1.0) {
        return 1.0;
      }
      return (1.0-total_prob_data)/(1.0-total_prob_simu);
    }
    if (total_prob_simu<=0.0)
      return 1.0;
    return total_prob_data/total_prob_simu;
  });

  const NamedFunc trigger_dependent_ptcut("trigger_dependent_ptcut",[hlt_el_trigger,hlt_single_el_trigger,hlt_mu_trigger,hlt_single_mu_trigger,Electron_hltId,Muon_hltId](const Baby &b) -> NamedFunc::ScalarType{
      std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
      std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
      std::vector<double> el_hlt_ = Electron_hltId.GetVector(b);
      std::vector<double> mu_hlt_ = Muon_hltId.GetVector(b);
      if (hlt_el_trigger.GetScalar(b)) {
        if (Lead_SignalElectron_pt.GetScalar(b)>25 && Sublead_SignalElectron_pt.GetScalar(b)>15)
          return 1;
        return 0;
      }
      else if (hlt_mu_trigger.GetScalar(b)) {
        if (Lead_SignalMuon_pt.GetScalar(b)>20 && Sublead_SignalMuon_pt.GetScalar(b)>10)
          return 1;
        return 0;
      }
      else if (hlt_single_el_trigger.GetScalar(b)) {
        //check if lead electron is reco'd and Id'd at HLT
        bool good_lead = false;
        for (unsigned iel = 0; iel < Electron_sig_.size(); iel++) {
          if (Electron_sig_[iel]) {
            if (el_hlt_[iel] >= 2) {
              good_lead = true;
            }
            break;
          }
        }
        if (good_lead) {
          if (Lead_SignalElectron_pt.GetScalar(b)>35)
            return 1;
        }
        else {
          if (Sublead_SignalElectron_pt.GetScalar(b)>35)
            return 1;
        }
        return 0;
      }
      else if (hlt_single_mu_trigger.GetScalar(b)) {
        //check if lead electron is reco'd and Id'd at HLT
        bool good_lead = false;
        for (unsigned imu = 0; imu < Electron_sig_.size(); imu++) {
          if (Muon_sig_[imu]) {
            if (mu_hlt_[imu] >= 2) {
              good_lead = true;
            }
            break;
          }
        }
        if (good_lead) {
          if (Lead_SignalMuon_pt.GetScalar(b)>29)
            return 1;
        }
        else {
          if (Sublead_SignalMuon_pt.GetScalar(b)>29)
            return 1;
        }
        return 0;
      }
      return 0;
  });

  const NamedFunc OneElectronOutOfAcceptance("OneElectronOutOfAcceptance",[](const Baby &b) -> NamedFunc::ScalarType{
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==11) {
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==23) {
            if (abs(b.GenPart_eta()->at(imc))>2.55) {
              return 1;
            }
          }
        }
      }
    }
    return 0;
  });

  const NamedFunc OneMuonOutOfAcceptance("OneMuonOutOfAcceptance",[](const Baby &b) -> NamedFunc::ScalarType{
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==13) {
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==23) {
            if (abs(b.GenPart_eta()->at(imc))>2.45) {
              return 1;
            }
          }
        }
      }
    }
    return 0;
  });

  const NamedFunc GenPart_HZ_Lead_Lepton_pt("GenPart_HZ_Lead_Lepton_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> hz_lepton_pt;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==13||abs_pdgid==11) { //is a lepton
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==23) {//is from Z
            int z_mom_idx = b.GenPart_genPartIdxMother()->at(mom_idx);
            while (true) {
              if (z_mom_idx != -1) {
                int z_mom_id = b.GenPart_pdgId()->at(z_mom_idx);
                if (z_mom_id == 23) {
                  z_mom_idx = b.GenPart_genPartIdxMother()->at(z_mom_idx);
                }
                else {
                  if (z_mom_id == 25) { //Z is from H
                    hz_lepton_pt.push_back(b.GenPart_pt()->at(imc));
                  }
                  break;
                }
              }
              else {
                break;
              }
            } //looking for Z parent
          } //lepton from Z
        } //lepton has parent
      } //particle is a lepton
    }
    if (hz_lepton_pt.size()>0) {
      std::sort(hz_lepton_pt.begin(), hz_lepton_pt.end());
      return hz_lepton_pt[hz_lepton_pt.size()-1];
    }
    return -999;
  });

  const NamedFunc GenPart_HZ_Sublead_Lepton_pt("GenPart_HZ_Sublead_Lepton_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> hz_lepton_pt;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==13||abs_pdgid==11) { //is a lepton
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==23) {//is from Z
            int z_mom_idx = b.GenPart_genPartIdxMother()->at(mom_idx);
            while (true) {
              if (z_mom_idx != -1) {
                int z_mom_id = b.GenPart_pdgId()->at(z_mom_idx);
                if (z_mom_id == 23) {
                  z_mom_idx = b.GenPart_genPartIdxMother()->at(z_mom_idx);
                }
                else {
                  if (z_mom_id == 25) { //Z is from H
                    hz_lepton_pt.push_back(b.GenPart_pt()->at(imc));
                  }
                  break;
                }
              }
              else {
                break;
              }
            } //looking for Z parent
          } //lepton from Z
        } //lepton has parent
      } //particle is a lepton
    }
    if (hz_lepton_pt.size()>1) {
      std::sort(hz_lepton_pt.begin(), hz_lepton_pt.end());
      return hz_lepton_pt[hz_lepton_pt.size()-2];
    }
    return -999;
  });

  const NamedFunc GenPart_HZ_Sublead_Lepton_eta("GenPart_HZ_Sublead_Lepton_eta",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<std::pair<double,double>> hz_lepton_pt;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==13||abs_pdgid==11) { //is a lepton
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==23) {//is from Z
            int z_mom_idx = b.GenPart_genPartIdxMother()->at(mom_idx);
            while (true) {
              if (z_mom_idx != -1) {
                int z_mom_id = b.GenPart_pdgId()->at(z_mom_idx);
                if (z_mom_id == 23) {
                  z_mom_idx = b.GenPart_genPartIdxMother()->at(z_mom_idx);
                }
                else {
                  if (z_mom_id == 25) { //Z is from H
                    hz_lepton_pt.push_back(std::pair<double,double>(b.GenPart_pt()->at(imc),b.GenPart_eta()->at(imc)));
                  }
                  break;
                }
              }
              else {
                break;
              }
            } //looking for Z parent
          } //lepton from Z
        } //lepton has parent
      } //particle is a lepton
    }
    if (hz_lepton_pt.size()>1) {
      std::sort(hz_lepton_pt.begin(), hz_lepton_pt.end(), [](std::pair<double,double> a, std::pair<double,double> c){return a.first<c.first;});
      return hz_lepton_pt[hz_lepton_pt.size()-2].second;
    }
    return -999;
  });

  const NamedFunc GenPart_HZ_Lead_Lepton_eta("GenPart_HZ_Lead_Lepton_eta",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<std::pair<double,double>> hz_lepton_pt;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==13||abs_pdgid==11) { //is a lepton
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==23) {//is from Z
            int z_mom_idx = b.GenPart_genPartIdxMother()->at(mom_idx);
            while (true) {
              if (z_mom_idx != -1) {
                int z_mom_id = b.GenPart_pdgId()->at(z_mom_idx);
                if (z_mom_id == 23) {
                  z_mom_idx = b.GenPart_genPartIdxMother()->at(z_mom_idx);
                }
                else {
                  if (z_mom_id == 25) { //Z is from H
                    hz_lepton_pt.push_back(std::pair<double,double>(b.GenPart_pt()->at(imc),b.GenPart_eta()->at(imc)));
                  }
                  break;
                }
              }
              else {
                break;
              }
            } //looking for Z parent
          } //lepton from Z
        } //lepton has parent
      } //particle is a lepton
    }
    if (hz_lepton_pt.size()>0) {
      std::sort(hz_lepton_pt.begin(), hz_lepton_pt.end(), [](std::pair<double,double> a, std::pair<double,double> c){return a.first<c.first;});
      return hz_lepton_pt[hz_lepton_pt.size()-1].second;
    }
    return -999;
  });

  const NamedFunc GenPart_NZ_Lead_Lepton_pt("GenPart_NZ_Lead_Lepton_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> hz_lepton_pt;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==13||abs_pdgid==11) { //is a lepton
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==23) {//is from Z
            int z_mom_idx = b.GenPart_genPartIdxMother()->at(mom_idx);
            while (true) {
              if (z_mom_idx != -1) {
                int z_mom_id = b.GenPart_pdgId()->at(z_mom_idx);
                if (z_mom_id == 23) {
                  z_mom_idx = b.GenPart_genPartIdxMother()->at(z_mom_idx);
                }
                else {
                  if (z_mom_id != 25) { //not H
                    hz_lepton_pt.push_back(b.GenPart_pt()->at(imc));
                  }
                  break;
                }
              }
              else {
                hz_lepton_pt.push_back(b.GenPart_pt()->at(imc));
                break;
              }
            } //looking for Z parent
          } //lepton from Z
        } //lepton has parent
      } //particle is a lepton
    }
    if (hz_lepton_pt.size()>0) {
      std::sort(hz_lepton_pt.begin(), hz_lepton_pt.end());
      return hz_lepton_pt[hz_lepton_pt.size()-1];
    }
    return -999;
  });

  const NamedFunc GenPart_NZ_Sublead_Lepton_pt("GenPart_NZ_Sublead_Lepton_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> hz_lepton_pt;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==13||abs_pdgid==11) { //is a lepton
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==23) {//is from Z
            int z_mom_idx = b.GenPart_genPartIdxMother()->at(mom_idx);
            while (true) {
              if (z_mom_idx != -1) {
                int z_mom_id = b.GenPart_pdgId()->at(z_mom_idx);
                if (z_mom_id == 23) {
                  z_mom_idx = b.GenPart_genPartIdxMother()->at(z_mom_idx);
                }
                else {
                  if (z_mom_id != 25) { //not H
                    hz_lepton_pt.push_back(b.GenPart_pt()->at(imc));
                  }
                  break;
                }
              }
              else {
                hz_lepton_pt.push_back(b.GenPart_pt()->at(imc));
                break;
              }
            } //looking for Z parent
          } //lepton from Z
        } //lepton has parent
      } //particle is a lepton
    }
    if (hz_lepton_pt.size()>1) {
      std::sort(hz_lepton_pt.begin(), hz_lepton_pt.end());
      return hz_lepton_pt[hz_lepton_pt.size()-2];
    }
    return -999;
  });

  const NamedFunc GenPart_NZ_Sublead_Lepton_eta("GenPart_NZ_Sublead_Lepton_eta",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<std::pair<double,double>> hz_lepton_pt;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==13||abs_pdgid==11) { //is a lepton
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==23) {//is from Z
            int z_mom_idx = b.GenPart_genPartIdxMother()->at(mom_idx);
            while (true) {
              if (z_mom_idx != -1) {
                int z_mom_id = b.GenPart_pdgId()->at(z_mom_idx);
                if (z_mom_id == 23) {
                  z_mom_idx = b.GenPart_genPartIdxMother()->at(z_mom_idx);
                }
                else {
                  if (z_mom_id != 25) { //not H
                    hz_lepton_pt.push_back(std::pair<double,double>(b.GenPart_pt()->at(imc),b.GenPart_eta()->at(imc)));
                  }
                  break;
                }
              }
              else {
                hz_lepton_pt.push_back(std::pair<double,double>(b.GenPart_pt()->at(imc),b.GenPart_eta()->at(imc)));
                break;
              }
            } //looking for Z parent
          } //lepton from Z
        } //lepton has parent
      } //particle is a lepton
    }
    if (hz_lepton_pt.size()>1) {
      std::sort(hz_lepton_pt.begin(), hz_lepton_pt.end(), [](std::pair<double,double> a, std::pair<double,double> c){return a.first<c.first;});
      return hz_lepton_pt[hz_lepton_pt.size()-2].second;
    }
    return -999;
  });

  const NamedFunc GenPart_NZ_Lead_Lepton_eta("GenPart_NZ_Lead_Lepton_eta",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<std::pair<double,double>> hz_lepton_pt;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==13||abs_pdgid==11) { //is a lepton
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==23) {//is from Z
            int z_mom_idx = b.GenPart_genPartIdxMother()->at(mom_idx);
            while (true) {
              if (z_mom_idx != -1) {
                int z_mom_id = b.GenPart_pdgId()->at(z_mom_idx);
                if (z_mom_id == 23) {
                  z_mom_idx = b.GenPart_genPartIdxMother()->at(z_mom_idx);
                }
                else {
                  if (z_mom_id != 25) { //not H
                    hz_lepton_pt.push_back(std::pair<double,double>(b.GenPart_pt()->at(imc),b.GenPart_eta()->at(imc)));
                  }
                  break;
                }
              }
              else {
                hz_lepton_pt.push_back(std::pair<double,double>(b.GenPart_pt()->at(imc),b.GenPart_eta()->at(imc)));
                break;
              }
            } //looking for Z parent
          } //lepton from Z
        } //lepton has parent
      } //particle is a lepton
    }
    if (hz_lepton_pt.size()>0) {
      std::sort(hz_lepton_pt.begin(), hz_lepton_pt.end(), [](std::pair<double,double> a, std::pair<double,double> c){return a.first<c.first;});
      return hz_lepton_pt[hz_lepton_pt.size()-1].second;
    }
    return -999;
  });

  const NamedFunc GenPart_W_Sublead_Electron_eta("GenPart_W_Sublead_Electron_eta",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<std::vector<double>> hz_lepton_pt;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==11 || abs_pdgid==13) { //is a lepton
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==24 || b.GenPart_pdgId()->at(mom_idx)==-24) {//is from W
            hz_lepton_pt.push_back({b.GenPart_pt()->at(imc),b.GenPart_eta()->at(imc),static_cast<double>(abs_pdgid)});
          } //lepton from W
        } //lepton has parent
      } //particle is a lepton
    }
    if (hz_lepton_pt.size()>1) {
      std::sort(hz_lepton_pt.begin(), hz_lepton_pt.end(), [](std::vector<double> a, std::vector<double> c){return a[0]<c[0];});
      if (hz_lepton_pt[hz_lepton_pt.size()-2][2]<=11) {
        return hz_lepton_pt[hz_lepton_pt.size()-2][1];
      }
    }
    return -999;
  });

  const NamedFunc GenPart_W_Lead_Electron_eta("GenPart_W_Lead_Electron_eta",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<std::vector<double>> hz_lepton_pt;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==11 || abs_pdgid==13) { //is a lepton
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==24 || b.GenPart_pdgId()->at(mom_idx)==-24) {//is from W
            hz_lepton_pt.push_back({b.GenPart_pt()->at(imc),b.GenPart_eta()->at(imc),static_cast<double>(abs_pdgid)});
          } //lepton from W
        } //lepton has parent
      } //particle is a lepton
    }
    if (hz_lepton_pt.size()>0) {
      std::sort(hz_lepton_pt.begin(), hz_lepton_pt.end(), [](std::vector<double> a, std::vector<double> c){return a[0]<c[0];});
      if (hz_lepton_pt[hz_lepton_pt.size()-1][2]<=11) {
        return hz_lepton_pt[hz_lepton_pt.size()-1][1];
      }
    }
    return -999;
  });

  const NamedFunc GenPart_W_Sublead_Muon_eta("GenPart_W_Sublead_Muon_eta",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<std::vector<double>> hz_lepton_pt;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==11 || abs_pdgid==13) { //is a lepton
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==24 || b.GenPart_pdgId()->at(mom_idx)==-24) {//is from W
            hz_lepton_pt.push_back({b.GenPart_pt()->at(imc),b.GenPart_eta()->at(imc),static_cast<double>(abs_pdgid)});
          } //lepton from W
        } //lepton has parent
      } //particle is a lepton
    }
    if (hz_lepton_pt.size()>1) {
      std::sort(hz_lepton_pt.begin(), hz_lepton_pt.end(), [](std::vector<double> a, std::vector<double> c){return a[0]<c[0];});
      if (hz_lepton_pt[hz_lepton_pt.size()-2][2]>11) {
        return hz_lepton_pt[hz_lepton_pt.size()-2][1];
      }
    }
    return -999;
  });

  const NamedFunc GenPart_W_Lead_Muon_eta("GenPart_W_Lead_Muon_eta",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<std::vector<double>> hz_lepton_pt;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==11 || abs_pdgid==13) { //is a lepton
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==24 || b.GenPart_pdgId()->at(mom_idx)==-24) {//is from W
            hz_lepton_pt.push_back({b.GenPart_pt()->at(imc),b.GenPart_eta()->at(imc),static_cast<double>(abs_pdgid)});
          } //lepton from W
        } //lepton has parent
      } //particle is a lepton
    }
    if (hz_lepton_pt.size()>0) {
      std::sort(hz_lepton_pt.begin(), hz_lepton_pt.end(), [](std::vector<double> a, std::vector<double> c){return a[0]<c[0];});
      if (hz_lepton_pt[hz_lepton_pt.size()-1][2]>11) {
        return hz_lepton_pt[hz_lepton_pt.size()-1][1];
      }
    }
    return -999;
  });

  const NamedFunc GenPart_W_Sublead_SignalMuon_pt("GenPart_W_Sublead_SignalMuon_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<std::vector<double>> hz_lepton_pt;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==11 || abs_pdgid==13) { //is a lepton
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==24 || b.GenPart_pdgId()->at(mom_idx)==-24) {//is from W
            hz_lepton_pt.push_back({b.GenPart_pt()->at(imc),static_cast<double>(abs_pdgid)});
          } //lepton from W
        } //lepton has parent
      } //particle is a lepton
    }
    if (hz_lepton_pt.size()>1) {
      std::sort(hz_lepton_pt.begin(), hz_lepton_pt.end(), [](std::vector<double> a, std::vector<double> c){return a[0]<c[0];});
      if (hz_lepton_pt[hz_lepton_pt.size()-2][1]>11) {
        return hz_lepton_pt[hz_lepton_pt.size()-2][0];
      }
    }
    return -999;
  });

  const NamedFunc GenPart_W_Lead_SignalMuon_pt("GenPart_W_Lead_SignalMuon_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<std::vector<double>> hz_lepton_pt;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==11 || abs_pdgid==13) { //is a lepton
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==24 || b.GenPart_pdgId()->at(mom_idx)==-24) {//is from W
            hz_lepton_pt.push_back({b.GenPart_pt()->at(imc),static_cast<double>(abs_pdgid)});
          } //lepton from W
        } //lepton has parent
      } //particle is a lepton
    }
    if (hz_lepton_pt.size()>0) {
      std::sort(hz_lepton_pt.begin(), hz_lepton_pt.end(), [](std::vector<double> a, std::vector<double> c){return a[0]<c[0];});
      if (hz_lepton_pt[hz_lepton_pt.size()-1][1]>11) {
        return hz_lepton_pt[hz_lepton_pt.size()-1][0];
      }
    }
    return -999;
  });

  const NamedFunc GenPart_W_Sublead_SignalElectron_pt("GenPart_W_Sublead_SignalElectron_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<std::vector<double>> hz_lepton_pt;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==11 || abs_pdgid==13) { //is a lepton
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==24 || b.GenPart_pdgId()->at(mom_idx)==-24) {//is from W
            hz_lepton_pt.push_back({b.GenPart_pt()->at(imc),static_cast<double>(abs_pdgid)});
          } //lepton from W
        } //lepton has parent
      } //particle is a lepton
    }
    if (hz_lepton_pt.size()>1) {
      std::sort(hz_lepton_pt.begin(), hz_lepton_pt.end(), [](std::vector<double> a, std::vector<double> c){return a[0]<c[0];});
      if (hz_lepton_pt[hz_lepton_pt.size()-2][1]<=11) {
        return hz_lepton_pt[hz_lepton_pt.size()-2][0];
      }
    }
    return -999;
  });

  const NamedFunc GenPart_W_Lead_SignalElectron_pt("GenPart_W_Lead_SignalElectron_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<std::vector<double>> hz_lepton_pt;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==11 || abs_pdgid==13) { //is a lepton
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==24 || b.GenPart_pdgId()->at(mom_idx)==-24) {//is from W
            hz_lepton_pt.push_back({b.GenPart_pt()->at(imc),static_cast<double>(abs_pdgid)});
          } //lepton from W
        } //lepton has parent
      } //particle is a lepton
    }
    if (hz_lepton_pt.size()>0) {
      std::sort(hz_lepton_pt.begin(), hz_lepton_pt.end(), [](std::vector<double> a, std::vector<double> c){return a[0]<c[0];});
      if (hz_lepton_pt[hz_lepton_pt.size()-1][1]<=11) {
        return hz_lepton_pt[hz_lepton_pt.size()-1][0];
      }
    }
    return -999;
  });

  const NamedFunc H_Photon_pt("H_Photon_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==22) { //is a photon
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if (b.GenPart_pdgId()->at(mom_idx)==25) { //parent is H
            return b.GenPart_pt()->at(imc);
          }
        }
      } 
    }
    return -999;
  });

  const NamedFunc H_Photon_eta("H_Photon_eta",[](const Baby &b) -> NamedFunc::ScalarType{
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==22) { //is a photon
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if (b.GenPart_pdgId()->at(mom_idx)==25) { //parent is H
            return b.GenPart_eta()->at(imc);
          }
        }
      } 
    }
    return -999;
  });

  const NamedFunc GenPart_MotherId("GenPart_MotherId",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> ret;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
      if (mom_idx != -1) {
        ret.push_back(static_cast<double>(b.GenPart_pdgId()->at(mom_idx)));
      } 
      else {
        ret.push_back(0);
      }
    }
    return ret;
  });

  //count number of leptonic ws
  const NamedFunc N_Leptonic_W("N_Leptonic_W",[](const Baby &b) -> NamedFunc::ScalarType{
    int n_w_lep = 0;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==13 || abs_pdgid==11 || abs_pdgid==15) {
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if (b.GenPart_pdgId()->at(mom_idx)==24 || b.GenPart_pdgId()->at(mom_idx)==-24) {//is from W
            n_w_lep++;
          } 
        } 
      } 
    }
    return n_w_lep;
  });

  //count number of leptonic zs
  const NamedFunc N_Leptonic_Z("N_Leptonic_Z",[](const Baby &b) -> NamedFunc::ScalarType{
    int n_z_lep = 0;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int pdgid = b.GenPart_pdgId()->at(imc);
      if (pdgid==13 || pdgid==11 || pdgid==15) {
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==23) {//is from Z
            n_z_lep++;
          } 
        } 
      } 
    }
    return n_z_lep;
  });

  //number of leptons in detector acceptance
  const NamedFunc NTruthLeptonAccept("NTruthLeptonAccept",[](const Baby &b) -> NamedFunc::ScalarType{
    double n_truth = 0;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==11||abs_pdgid==13) { //is lepton
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx)>=21 && b.GenPart_pdgId()->at(mom_idx)<=25) || b.GenPart_pdgId()->at(mom_idx)==-24) { //prompt
            float pt_cut = 7;
            float eta_cut = 2.5;
            if (abs_pdgid == 13) {
              pt_cut = 5;
              eta_cut = 2.4;
            }
            if (b.GenPart_pt()->at(imc) > pt_cut && abs(b.GenPart_eta()->at(imc))<eta_cut) {
              n_truth++;
            }
          }
        }
      } 
    }
    return n_truth;
  });

  //number of photons in detector acceptance
  const NamedFunc NTruthPhotonAccept("NTruthPhotonAccept",[](const Baby &b) -> NamedFunc::ScalarType{
    double n_truth = 0;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==22) { //is photon
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if (b.GenPart_pdgId()->at(mom_idx) != 22) { //first copy
            if (b.GenPart_pt()->at(imc) > 15 && abs(b.GenPart_eta()->at(imc))<2.5) {
              n_truth++;
            }
          }
        }
      } 
    }
    return n_truth;
  });

  //hem veto
  const NamedFunc pass_hemveto("pass_hemveto", [](const Baby &b) -> NamedFunc::ScalarType{
      std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
      //only apply for 2018 era C+D and MC
      if (abs(b.SampleType())!=2018) return true; // accept 2016 and 2017 data and mc
      if (b.SampleType()==-2018 && b.run() < 319077) return true; // accept part of 2018 data
      if (b.SampleType()==2018 && (b.event()%1961) >= 1296) return true; //accept part of 2018 mc
      bool pass_hem = true;
      for (unsigned int el_idx = 0; el_idx < b.Electron_pt()->size(); el_idx++) {
        if (b.Electron_miniPFRelIso_chg()->at(el_idx) < 0.1 && -3.0 < b.Electron_eta()->at(el_idx) && b.Electron_eta()->at(el_idx) < -1.4 && -1.57 < b.Electron_phi()->at(el_idx) && b.Electron_phi()->at(el_idx) < -0.87) {
          pass_hem = false;
        }
      }
      for (unsigned int jet_idx = 0; jet_idx < b.Jet_pt()->size(); jet_idx++) {
        if (b.Jet_pt()->at(jet_idx) > 30. && -3.2 < b.Jet_eta()->at(jet_idx) && b.Jet_eta()->at(jet_idx) < -1.2 && -1.77 < b.Jet_phi()->at(jet_idx) && b.Jet_phi()->at(jet_idx) < -0.67) {
          double dphi = fabs(TVector2::Phi_mpi_pi(b.Jet_phi()->at(jet_idx)-b.MET_phi()));
          if (dphi < 0.5) {
            pass_hem = false;
          }
        }
      }
      return pass_hem;
  });

  NamedFunc w_sigx100("w_sigx100",[](const Baby &b) -> NamedFunc::ScalarType{
    if(b.type() >= 200000 && b.type() < 201000)
      return 100.0;
    return 1.0;
  });

  const NamedFunc trig_double_photon = "HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId||HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90";

  const NamedFunc sanitize_weight = NamedFunc("weight/w_lumi<10")
      .Name("sanitize_weight");

  NamedFunc baseline_selection = ((ZCandidate_mass>50)&&((Lead_SignalElectron_pt>25)||(Lead_SignalMuon_pt>20))&&(Lead_SignalPhoton_pt/HiggsCandidate_mass>0.14)&&delta_r_cut&&((HiggsCandidate_mass+ZCandidate_mass)>185));
  NamedFunc baseline_selection_nopt = ((ZCandidate_mass>50)&&(Lead_SignalPhoton_pt/HiggsCandidate_mass>0.14)&&delta_r_cut&&((HiggsCandidate_mass+ZCandidate_mass)>185));
  NamedFunc baseline_selection_nodr = ((ZCandidate_mass>50)&&((Lead_SignalElectron_pt>25)||(Lead_SignalMuon_pt>20))&&(Lead_SignalPhoton_pt/HiggsCandidate_mass>0.14)&&((HiggsCandidate_mass+ZCandidate_mass)>185));
  NamedFunc midleadpt = ((Lead_SignalElectron_pt>25)||(Lead_SignalMuon_pt>20));
  NamedFunc highleadpt = ((Lead_SignalElectron_pt>38)||(Lead_SignalMuon_pt>29));
  NamedFunc highpt = ((Lead_SignalElectron_pt>38&&Sublead_SignalElectron_pt>38)||(Lead_SignalMuon_pt>29&&Sublead_SignalMuon_pt>29));
  NamedFunc dy_selection = ((ZCandidate_mass>81)&&(ZCandidate_mass<101))&&(nSignalElectron>=2||nSignalMuon>=2);

  const NamedFunc zg_baseline_notrig = NamedFunc("(nel>=2||nmu>=2) "
      "&& nphoton>=1 && ((photon_pt[0]/llphoton_m[0])>=15.0/110.0) "
      "&& ll_m[0]>80 && ll_m[0]<100 && ((llphoton_m[0]+ll_m[0])>=185) "
      "&& llphoton_m[0]>100 && llphoton_m[0]<180").Name("baseline_notrig");

  const NamedFunc mc_el_pt = FilterNamedFunc("mc_pt",
      "(mc_id==11||mc_id==-11)&&mc_status==1");
  const NamedFunc mc_mu_pt = FilterNamedFunc("mc_pt",
      "(mc_id==13||mc_id==-13)&&mc_status==1");
  const NamedFunc mc_photon_pt = FilterNamedFunc("mc_pt",
      "(mc_id==22)&&mc_status==1");
  const NamedFunc lead_mc_el_pt = ReduceNamedFunc(mc_el_pt,reduce_max)
                                  .Name("lead_mc_el_pt");
  const NamedFunc sublead_mc_el_pt = ReduceNamedFunc(mc_el_pt,reduce_sublead)
                                     .Name("sublead_mc_el_pt");
  const NamedFunc lead_mc_mu_pt = ReduceNamedFunc(mc_mu_pt,reduce_max)
                                  .Name("lead_mc_mu_pt");
  const NamedFunc sublead_mc_mu_pt = ReduceNamedFunc(mc_mu_pt,reduce_sublead)
                                     .Name("sublead_mc_mu_pt");
  const NamedFunc lead_mc_photon_pt = ReduceNamedFunc(mc_photon_pt,reduce_max)
                                      .Name("lead_mc_photon_pt");

  //---------------------------------------------------------------------------
  //                                     make plots and pie charts
  //---------------------------------------------------------------------------
  
  PlotMaker pm;

  bool plot_fail_dilep_trig = false;
  bool plot_fail_singlelep_trig = false;
  bool plot_more_fail_dilep = false;
  bool plot_cleanmask = false;
  bool plot_single_lep_trig_effects = false;
  bool plot_single_lep_trig_effects_sig = false;
  bool plot_diphoton_trig_effects = false;
  bool make_rereco_plots = false;
  bool make_table_fail_reason = false;
  bool make_compare_preselection_table = false;
  bool plot_trigeffs = false;
  bool plot_trigeffs_signal = false;
  bool make_datamctrigefftable = false;
  bool make_trigcomp_plots = false;
  bool make_mu17ref = false;
  bool make_dy_efftable = false;
  bool plot_jaebak_cutflow = false;
  bool plot_associated_acceptance = false;
  bool plot_ptg = false;
  bool plot_data_sidebands = false;
  bool make_anplots_signal = true;
  bool make_anplots_signal_background = false;

  if (make_anplots_signal) {
    //trigger efficiency plots
    pm.Push<EfficiencyPlot>(
        Axis(30, 0.0, 100.0, "el_pt[0]", "Lead electron p_{T} [GeV]", {}),
        zg_baseline_notrig&&sanitize_weight&&"ll_lepid[0]==11", 
        "trig_double_el",
        procs_sig_pinnacles,true,ops_mc).Weight("weight"*w_years)
        .Tag("zgtrig").YTitle("trigger").YAxisMax(1.4);
    pm.Push<EfficiencyPlot>(
        Axis(30, 0.0, 100.0, "el_pt[0]", "Lead electron p_{T} [GeV]", {}),
        zg_baseline_notrig&&sanitize_weight&&"ll_lepid[0]==11", 
        "trig_single_el||trig_double_el",
        procs_sig_pinnacles,true,ops_mc).Weight("weight"*w_years)
        .Tag("zgtrig").YTitle("trigger").YAxisMax(1.4);
    pm.Push<EfficiencyPlot>(
        Axis(30, 0.0, 100.0, "el_pt[1]", "Sublead electron p_{T} [GeV]", {}),
        zg_baseline_notrig&&sanitize_weight&&"ll_lepid[0]==11", 
        "trig_double_el",
        procs_sig_pinnacles,true,ops_mc).Weight("weight"*w_years)
        .Tag("zgtrig").YTitle("trigger").YAxisMax(1.4);
    pm.Push<EfficiencyPlot>(
        Axis(30, 0.0, 100.0, "el_pt[1]", "Sublead electron p_{T} [GeV]", {}),
        zg_baseline_notrig&&sanitize_weight&&"ll_lepid[0]==11", 
        "trig_single_el||trig_double_el",
        procs_sig_pinnacles,true,ops_mc).Weight("weight"*w_years)
        .Tag("zgtrig").YTitle("trigger").YAxisMax(1.4);
    pm.Push<EfficiencyPlot>(
        Axis(30, 0.0, 100.0, "mu_pt[0]", "Lead muon p_{T} [GeV]", {}),
        zg_baseline_notrig&&sanitize_weight&&"ll_lepid[0]==13", 
        "trig_double_mu",
        procs_sig_pinnacles,true,ops_mc).Weight("weight"*w_years)
        .Tag("zgtrig").YTitle("trigger").YAxisMax(1.4);
    pm.Push<EfficiencyPlot>(
        Axis(30, 0.0, 100.0, "mu_pt[0]", "Lead muon p_{T} [GeV]", {}),
        zg_baseline_notrig&&sanitize_weight&&"ll_lepid[0]==13", 
        "trig_single_mu||trig_double_mu",
        procs_sig_pinnacles,true,ops_mc).Weight("weight"*w_years)
        .Tag("zgtrig").YTitle("trigger").YAxisMax(1.4);
    pm.Push<EfficiencyPlot>(
        Axis(30, 0.0, 100.0, "mu_pt[1]", "Sublead muon p_{T} [GeV]", {}),
        zg_baseline_notrig&&sanitize_weight&&"ll_lepid[0]==13", 
        "trig_double_mu",
        procs_sig_pinnacles,true,ops_mc).Weight("weight"*w_years)
        .Tag("zgtrig").YTitle("trigger").YAxisMax(1.4);
    pm.Push<EfficiencyPlot>(
        Axis(30, 0.0, 100.0, "mu_pt[1]", "Sublead muon p_{T} [GeV]", {}),
        zg_baseline_notrig&&sanitize_weight&&"ll_lepid[0]==13", 
        "trig_single_mu||trig_double_mu",
        procs_sig_pinnacles,true,ops_mc).Weight("weight"*w_years)
        .Tag("zgtrig").YTitle("trigger").YAxisMax(1.4);
    //plots analyzing events added by single lepton triggers for signal
    vector<NamedFunc> channels = {Ztoee, ZtoMuMu};
    vector<NamedFunc> lep_pt = {"el_pt","mu_pt"};
    vector<NamedFunc> lep_eta = {"el_eta","mu_eta"};
    vector<string> lep_name = {"Electron","Muon"};
    vector<NamedFunc> lep_sig = {"el_sig","mu_sig"};
    for (unsigned ich = 0; ich < channels.size(); ich++) {
      pm.Push<Hist1D>(
          Axis(30, 110.0, 140.0, "llphoton_m[0]","m_{ll#gamma} [GeV]"),
          sanitize_weight&&zg_baseline_notrig&&channels[ich],
          procs_sig_pinnacles_compare, ops_shapes).Weight("weight"*w_years)
          .Tag("zgtrig");
      pm.Push<Hist1D>(
          Axis(30, 0.0, 100.0, lep_pt[ich], lep_name[ich]+" p_{T} [GeV]"),
          sanitize_weight&&zg_baseline_notrig&&channels[ich]&&lep_sig[ich],
          procs_sig_pinnacles_compare, ops_shapes).Weight("weight"*w_years)
          .Tag("zgtrig");
      pm.Push<Hist1D>(
          Axis(50, -2.5, 2.5, lep_eta[ich], lep_name[ich]+" #eta [GeV]"),
          sanitize_weight&&zg_baseline_notrig&&channels[ich]&&lep_sig[ich],
          procs_sig_pinnacles_compare, ops_shapes).Weight("weight"*w_years)
          .Tag("zgtrig");
      pm.Push<Hist1D>(
          Axis(30, 0.0, 100.0, "photon_pt[0]", "Photon p_{T} [GeV]"),
          sanitize_weight&&zg_baseline_notrig&&channels[ich],
          procs_sig_pinnacles_compare, ops_shapes).Weight("weight"*w_years)
          .Tag("zgtrig");
      pm.Push<Hist1D>(
          Axis(50, -2.5, 2.5, "photon_eta[0]", "Photon #eta [GeV]"),
          sanitize_weight&&zg_baseline_notrig&&channels[ich],
          procs_sig_pinnacles_compare, ops_shapes).Weight("weight"*w_years)
          .Tag("zgtrig");
    }
    //tables
    pm.Push<Table>("zgtrig_basic", vector<TableRow>{
      TableRow("Baseline (no trigger) events (Z$\\to$ ee)",
          Ztoee&&zg_baseline_notrig,0,0,"weight"*w_years),
      TableRow("Pass dilep trigger",
          Ztoee&&zg_baseline_notrig&&"trig_double_el||trig_double_mu",
          0,0,"weight"*w_years),
      TableRow("Pass single or dilep trigger",
          Ztoee&&zg_baseline_notrig
          &&"trig_double_el||trig_double_mu||trig_single_el||trig_single_mu",
          0,0,"weight"*w_years),
      TableRow("Baseline (no trigger) events (Z$\\to \\mu\\mu$)",
          ZtoMuMu&&zg_baseline_notrig,0,0,"weight"*w_years),
      TableRow("Pass dilep trigger",
          ZtoMuMu&&zg_baseline_notrig&&"trig_double_el||trig_double_mu",
          0,0,"weight"*w_years),
      TableRow("Pass single or dilep trigger",
          ZtoMuMu&&zg_baseline_notrig
          &&"trig_double_el||trig_double_mu||trig_single_el||trig_single_mu",
          0,0,"weight"*w_years),
    },procs_sig_pinnacles,false,true,false,false,false,true).Precision(3);
    pm.Push<Table>("zgtrig_unweighted", vector<TableRow>{
      TableRow("Baseline (no trigger) events (Z$\\to$ ee)",
          Ztoee&&zg_baseline_notrig,0,0,"w_lumi"*w_years),
      TableRow("Pass dilep trigger",
          Ztoee&&zg_baseline_notrig&&"trig_double_el||trig_double_mu",
          0,0,"w_lumi"*w_years),
      TableRow("Pass single lepton trigger",
          Ztoee&&zg_baseline_notrig&&"trig_single_el||trig_single_mu",
          0,0,"w_lumi"*w_years),
      TableRow("Pass single or dilep trigger",
          Ztoee&&zg_baseline_notrig
          &&"trig_double_el||trig_double_mu||trig_single_el||trig_single_mu",
          0,0,"w_lumi"*w_years),
      TableRow("Pass single or dilep trigger or lepton+photon trigger",
          Ztoee&&zg_baseline_notrig
          &&("trig_double_el||trig_double_mu||trig_single_el||trig_single_mu"
          ||trig_double_photon),
          0,0,"w_lumi"*w_years),
      TableRow("Baseline (no trigger) events (Z$\\to \\mu\\mu$)",
          ZtoMuMu&&zg_baseline_notrig,0,0,"w_lumi"*w_years),
      TableRow("Pass dilep trigger",
          ZtoMuMu&&zg_baseline_notrig&&"trig_double_el||trig_double_mu",
          0,0,"w_lumi"*w_years),
      TableRow("Pass single lepton trigger",
          ZtoMuMu&&zg_baseline_notrig&&"trig_single_el||trig_single_mu",
          0,0,"w_lumi"*w_years),
      TableRow("Pass single or dilep trigger",
          ZtoMuMu&&zg_baseline_notrig
          &&("trig_double_el||trig_double_mu||trig_single_el||trig_single_mu"
          ||trig_double_photon),
          0,0,"w_lumi"*w_years),
      TableRow("Pass single or dilep trigger or muon+photon trigger",
          ZtoMuMu&&zg_baseline_notrig
          &&("trig_double_el||trig_double_mu||trig_single_el||trig_single_mu"
          ||NamedFunc("HLT_Mu17_Photon30_IsoCaloId")),
          0,0,"w_lumi"*w_years),
    },procs_sig_pinnacles,false,true,false,false,false,true).Precision(3);

    pm.SetLuminosityTag("138 fb^{-1} (13 TeV) + 62 fb^{-1} (13.6 TeV)");
  }

  if (make_anplots_signal_background) {
    //gen plots
    pm.Push<Hist1D>(
        Axis(30, 0.0, 140.0, lead_mc_el_pt, "Lead Electron p_{T}", {}),
        sanitize_weight&&Ztoee,
        procs_sig_pinnacles_unskimmed, ops_mc).Weight("weight"*w_years)
        .Tag("zgtrig");
    pm.Push<Hist1D>(
        Axis(30, 0.0, 100.0, sublead_mc_el_pt, 
             "Sublead Electron p_{T}", {}),
        sanitize_weight&&Ztoee,
        procs_sig_pinnacles_unskimmed, ops_mc).Weight("weight"*w_years)
        .Tag("zgtrig");
    pm.Push<Hist1D>(
        Axis(30, 0.0, 140.0, lead_mc_mu_pt, "Lead Muon p_{T}", {}),
        sanitize_weight&&ZtoMuMu,
        procs_sig_pinnacles_unskimmed, ops_mc).Weight("weight"*w_years)
        .Tag("zgtrig");
    pm.Push<Hist1D>(
        Axis(30, 0.0, 100.0, sublead_mc_mu_pt, "Sublead Muon p_{T}", {}),
        sanitize_weight&&ZtoMuMu,
        procs_sig_pinnacles_unskimmed, ops_mc).Weight("weight"*w_years)
        .Tag("zgtrig");
    pm.Push<Hist1D>(
        Axis(30, 0.0, 100.0, lead_mc_photon_pt, "Photon p_{T}", {}),
        sanitize_weight,
        procs_sig_pinnacles_unskimmed, ops_mc).Weight("weight"*w_years)
        .Tag("zgtrig");
    //signal/background plots
    pm.Push<Hist1D>(
        Axis(40, 100.0, 180.0, "llphoton_m[0]","m_{ll#gamma} [GeV]"),
        sanitize_weight&&zg_baseline_notrig&&"trig_double_el||trig_double_mu",
        procs_pinnacles_mc, ops_mc)
        .Weight("weight"*w_years*w_sigx100).Tag("zgtrig");
    pm.Push<Hist1D>(
        Axis(40, 100.0, 180.0, "llphoton_m[0]","m_{ll#gamma} [GeV]"),
        sanitize_weight&&zg_baseline_notrig&&
        "trig_double_el||trig_double_mu||trig_single_el||trig_single_mu",
        procs_pinnacles_mc, ops_mc)
        .Weight("weight"*w_years*w_sigx100).Tag("zgtrig");
    pm.SetLuminosityTag("138 fb^{-1} (13 TeV) + 62 fb^{-1} (13.6 TeV)");
  }

  if (plot_associated_acceptance) {
    //ggh
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_HZ_Lead_Lepton_pt, "Lead Electron p_{T}", {}),
        (z_decay_pdgid==11),
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__ggh_zlead_electron_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_HZ_Lead_Lepton_eta, "Lead Electron #eta", {}),
        (z_decay_pdgid==11),
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__ggh_zlead_electron_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_HZ_Sublead_Lepton_pt, "Sublead Electron p_{T}", {}),
        (z_decay_pdgid==11),
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__ggh_zsublead_electron_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_HZ_Sublead_Lepton_eta, "Sublead Electron #eta", {}),
        (z_decay_pdgid==11),
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__ggh_zsublead_electron_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_HZ_Lead_Lepton_pt, "Lead Muon p_{T}", {}),
        (z_decay_pdgid==13),
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__ggh_zlead_muon_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_HZ_Lead_Lepton_eta, "Lead Muon #eta", {}),
        (z_decay_pdgid==13),
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__ggh_zlead_muon_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_HZ_Sublead_Lepton_pt, "Sublead Muon p_{T}", {}),
        (z_decay_pdgid==13),
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__ggh_zsublead_muon_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_HZ_Sublead_Lepton_eta, "Sublead Muon #eta", {}),
        (z_decay_pdgid==13),
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__ggh_zsublead_muon_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, H_Photon_pt, "Photon p_{T}", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11),
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__ggh_photon_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, H_Photon_eta, "Photon #eta", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11),
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__ggh_photon_eta")
        .LuminosityTag(total_luminosity_string);
    //vbfh
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_HZ_Lead_Lepton_pt, "Lead Electron p_{T}", {}),
        (z_decay_pdgid==11),
        signal_procs_vbf, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__vbf_zlead_electron_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_HZ_Lead_Lepton_eta, "Lead Electron #eta", {}),
        (z_decay_pdgid==11),
        signal_procs_vbf, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__vbf_zlead_electron_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_HZ_Sublead_Lepton_pt, "Sublead Electron p_{T}", {}),
        (z_decay_pdgid==11),
        signal_procs_vbf, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__vbf_zsublead_electron_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_HZ_Sublead_Lepton_eta, "Sublead Electron #eta", {}),
        (z_decay_pdgid==11),
        signal_procs_vbf, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__vbf_zsublead_electron_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_HZ_Lead_Lepton_pt, "Lead Muon p_{T}", {}),
        (z_decay_pdgid==13),
        signal_procs_vbf, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__vbf_zlead_muon_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_HZ_Lead_Lepton_eta, "Lead Muon #eta", {}),
        (z_decay_pdgid==13),
        signal_procs_vbf, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__vbf_zlead_muon_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_HZ_Sublead_Lepton_pt, "Sublead Muon p_{T}", {}),
        (z_decay_pdgid==13),
        signal_procs_vbf, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__vbf_zsublead_muon_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_HZ_Sublead_Lepton_eta, "Sublead Muon #eta", {}),
        (z_decay_pdgid==13),
        signal_procs_vbf, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__vbf_zsublead_muon_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, H_Photon_pt, "Photon p_{T}", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11),
        signal_procs_vbf, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__vbf_photon_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, H_Photon_eta, "Photon #eta", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11),
        signal_procs_vbf, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__vbf_photon_eta")
        .LuminosityTag(total_luminosity_string);
    //wplush
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_HZ_Lead_Lepton_pt, "Lead Z Electron p_{T}", {}),
        (z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_wplush, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__wplush_zlead_electron_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_HZ_Lead_Lepton_eta, "Lead Z Electron #eta", {}),
        (z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_wplush, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__wplush_zlead_electron_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_HZ_Sublead_Lepton_pt, "Sublead Z Electron p_{T}", {}),
        (z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_wplush, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__wplush_zsublead_electron_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_HZ_Sublead_Lepton_eta, "Sublead Z Electron #eta", {}),
        (z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_wplush, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__wplush_zsublead_electron_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_HZ_Lead_Lepton_pt, "Lead Z Muon p_{T}", {}),
        (z_decay_pdgid==13)&&(N_Leptonic_W==1),
        signal_procs_wplush, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__wplush_zlead_muon_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_HZ_Lead_Lepton_eta, "Lead Z Muon #eta", {}),
        (z_decay_pdgid==13)&&(N_Leptonic_W==1),
        signal_procs_wplush, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__wplush_zlead_muon_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_HZ_Sublead_Lepton_pt, "Sublead Z Muon p_{T}", {}),
        (z_decay_pdgid==13)&&(N_Leptonic_W==1),
        signal_procs_wplush, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__wplush_zsublead_muon_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_HZ_Sublead_Lepton_eta, "Sublead Z Muon #eta", {}),
        (z_decay_pdgid==13)&&(N_Leptonic_W==1),
        signal_procs_wplush, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__wplush_zsublead_muon_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_W_Lead_SignalMuon_pt, "W Muon p_{T}", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_wplush, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__wplush_w_muon_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_W_Lead_Muon_eta, "W Muon #eta", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_wplush, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__wplush_w_muon_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_W_Lead_SignalElectron_pt, "W Electron p_{T}", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_wplush, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__wplush_w_electron_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_W_Lead_Electron_eta, "W Electron #eta", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_wplush, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__wplush_w_electron_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, H_Photon_pt, "Photon p_{T}", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_wplush, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__wplush_photon_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, H_Photon_eta, "Photon #eta", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_wplush, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__wplush_photon_eta")
        .LuminosityTag(total_luminosity_string);
    //wminush
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_HZ_Lead_Lepton_pt, "Lead Z Electron p_{T}", {}),
        (z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_wminush, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__wminush_zlead_electron_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_HZ_Lead_Lepton_eta, "Lead Z Electron #eta", {}),
        (z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_wminush, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__wminush_zlead_electron_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_HZ_Sublead_Lepton_pt, "Sublead Z Electron p_{T}", {}),
        (z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_wminush, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__wminush_zsublead_electron_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_HZ_Sublead_Lepton_eta, "Sublead Z Electron #eta", {}),
        (z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_wminush, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__wminush_zsublead_electron_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_HZ_Lead_Lepton_pt, "Lead Z Muon p_{T}", {}),
        (z_decay_pdgid==13)&&(N_Leptonic_W==1),
        signal_procs_wminush, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__wminush_zlead_muon_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_HZ_Lead_Lepton_eta, "Lead Z Muon #eta", {}),
        (z_decay_pdgid==13)&&(N_Leptonic_W==1),
        signal_procs_wminush, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__wminush_zlead_muon_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_HZ_Sublead_Lepton_pt, "Sublead Z Muon p_{T}", {}),
        (z_decay_pdgid==13)&&(N_Leptonic_W==1),
        signal_procs_wminush, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__wminush_zsublead_muon_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_HZ_Sublead_Lepton_eta, "Sublead Z Muon #eta", {}),
        (z_decay_pdgid==13)&&(N_Leptonic_W==1),
        signal_procs_wminush, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__wminush_zsublead_muon_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_W_Lead_SignalMuon_pt, "W Muon p_{T}", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_wminush, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__wminush_w_muon_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_W_Lead_Muon_eta, "W Muon #eta", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_wminush, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__wminush_w_muon_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_W_Lead_SignalElectron_pt, "W Electron p_{T}", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_wminush, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__wminush_w_electron_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_W_Lead_Electron_eta, "W Electron #eta", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_wminush, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__wminush_w_electron_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, H_Photon_pt, "Photon p_{T}", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_wminush, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__wminush_photon_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, H_Photon_eta, "Photon #eta", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_wminush, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__wminush_photon_eta")
        .LuminosityTag(total_luminosity_string);
    //zh
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_HZ_Lead_Lepton_pt, "Lead H Z Electron p_{T}", {}),
        (z_decay_pdgid==11)&&(N_Leptonic_Z==2),
        signal_procs_zh, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__zh_hlead_electron_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_HZ_Lead_Lepton_eta, "Lead H Z Electron #eta", {}),
        (z_decay_pdgid==11)&&(N_Leptonic_Z==2),
        signal_procs_zh, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__zh_hlead_electron_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_HZ_Sublead_Lepton_pt, "Sublead H Z Electron p_{T}", {}),
        (z_decay_pdgid==11)&&(N_Leptonic_Z==2),
        signal_procs_zh, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__zh_hsublead_electron_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_HZ_Sublead_Lepton_eta, "Sublead H Z Electron #eta", {}),
        (z_decay_pdgid==11)&&(N_Leptonic_Z==2),
        signal_procs_zh, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__zh_hsublead_electron_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_HZ_Lead_Lepton_pt, "Lead H Z Muon p_{T}", {}),
        (z_decay_pdgid==13)&&(N_Leptonic_Z==2),
        signal_procs_zh, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__zh_hlead_muon_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_HZ_Lead_Lepton_eta, "Lead H Z Muon #eta", {}),
        (z_decay_pdgid==13)&&(N_Leptonic_Z==2),
        signal_procs_zh, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__zh_hlead_muon_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_HZ_Sublead_Lepton_pt, "Sublead H Z Muon p_{T}", {}),
        (z_decay_pdgid==13)&&(N_Leptonic_Z==2),
        signal_procs_zh, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__zh_hsublead_muon_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_HZ_Sublead_Lepton_eta, "Sublead H Z Muon #eta", {}),
        (z_decay_pdgid==13)&&(N_Leptonic_Z==2),
        signal_procs_zh, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__zh_hsublead_muon_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_NZ_Lead_Lepton_pt, "Lead non-H Z Electron p_{T}", {}),
        (nh_z_decay_pdgid==11)&&(N_Leptonic_Z==2),
        signal_procs_zh, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__zh_zlead_electron_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_NZ_Lead_Lepton_eta, "Lead non-H Z Electron #eta", {}),
        (nh_z_decay_pdgid==11)&&(N_Leptonic_Z==2),
        signal_procs_zh, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__zh_zlead_electron_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_NZ_Sublead_Lepton_pt, "Sublead non-H Z Electron p_{T}", {}),
        (nh_z_decay_pdgid==11)&&(N_Leptonic_Z==2),
        signal_procs_zh, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__zh_zsublead_electron_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_NZ_Sublead_Lepton_eta, "Sublead non-H Z Electron #eta", {}),
        (nh_z_decay_pdgid==11)&&(N_Leptonic_Z==2),
        signal_procs_zh, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__zh_zsublead_electron_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_NZ_Lead_Lepton_pt, "Lead non-H Z Muon p_{T}", {}),
        (nh_z_decay_pdgid==13)&&(N_Leptonic_Z==2),
        signal_procs_zh, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__zh_zlead_muon_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_NZ_Lead_Lepton_eta, "Lead non-H Z Muon #eta", {}),
        (nh_z_decay_pdgid==13)&&(N_Leptonic_Z==2),
        signal_procs_zh, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__zh_zlead_muon_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_NZ_Sublead_Lepton_pt, "Sublead non-H Z Muon p_{T}", {}),
        (nh_z_decay_pdgid==13)&&(N_Leptonic_Z==2),
        signal_procs_zh, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__zh_zsublead_muon_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_NZ_Sublead_Lepton_eta, "Sublead non-H Z Muon #eta", {}),
        (nh_z_decay_pdgid==13)&&(N_Leptonic_Z==2),
        signal_procs_zh, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__zh_zsublead_muon_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, H_Photon_pt, "Photon p_{T}", {}),
        (N_Leptonic_Z==2),
        signal_procs_zh, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__zh_photon_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, H_Photon_eta, "Photon #eta", {}),
        (N_Leptonic_Z==2),
        signal_procs_zh, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__zh_photon_eta")
        .LuminosityTag(total_luminosity_string);
    //tth 3l
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_HZ_Lead_Lepton_pt, "Lead Z Electron p_{T}", {}),
        (z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth3l_zlead_electron_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_HZ_Lead_Lepton_eta, "Lead Z Electron #eta", {}),
        (z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth3l_zlead_electron_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_HZ_Sublead_Lepton_pt, "Sublead Z Electron p_{T}", {}),
        (z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth3l_zsublead_electron_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_HZ_Sublead_Lepton_eta, "Sublead Z Electron #eta", {}),
        (z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth3l_zsublead_electron_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_HZ_Lead_Lepton_pt, "Lead Z Muon p_{T}", {}),
        (z_decay_pdgid==13)&&(N_Leptonic_W==1),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth3l_zlead_muon_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_HZ_Lead_Lepton_eta, "Lead Z Muon #eta", {}),
        (z_decay_pdgid==13)&&(N_Leptonic_W==1),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth3l_zlead_muon_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_HZ_Sublead_Lepton_pt, "Sublead Z Muon p_{T}", {}),
        (z_decay_pdgid==13)&&(N_Leptonic_W==1),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth3l_zsublead_muon_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_HZ_Sublead_Lepton_eta, "Sublead Z Muon #eta", {}),
        (z_decay_pdgid==13)&&(N_Leptonic_W==1),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth3l_zsublead_muon_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_W_Lead_SignalMuon_pt, "W Muon p_{T}", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth3l_w_muon_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_W_Lead_Muon_eta, "W Muon #eta", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth3l_w_muon_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_W_Lead_SignalElectron_pt, "W Electron p_{T}", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth3l_w_electron_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_W_Lead_Electron_eta, "W Electron #eta", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth3l_w_electron_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, H_Photon_pt, "Photon p_{T}", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth3l_photon_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, H_Photon_eta, "Photon #eta", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11)&&(N_Leptonic_W==1),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth3l_photon_eta")
        .LuminosityTag(total_luminosity_string);
    //tth 4l
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_HZ_Lead_Lepton_pt, "Lead Z Electron p_{T}", {}),
        (z_decay_pdgid==11)&&(N_Leptonic_W==2),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth4l_zlead_electron_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_HZ_Lead_Lepton_eta, "Lead Z Electron #eta", {}),
        (z_decay_pdgid==11)&&(N_Leptonic_W==2),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth4l_zlead_electron_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_HZ_Sublead_Lepton_pt, "Sublead Z Electron p_{T}", {}),
        (z_decay_pdgid==11)&&(N_Leptonic_W==2),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth4l_zsublead_electron_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_HZ_Sublead_Lepton_eta, "Sublead Z Electron #eta", {}),
        (z_decay_pdgid==11)&&(N_Leptonic_W==2),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth4l_zsublead_electron_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_HZ_Lead_Lepton_pt, "Lead Z Muon p_{T}", {}),
        (z_decay_pdgid==13)&&(N_Leptonic_W==2),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth4l_zlead_muon_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_HZ_Lead_Lepton_eta, "Lead Z Muon #eta", {}),
        (z_decay_pdgid==13)&&(N_Leptonic_W==2),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth4l_zlead_muon_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_HZ_Sublead_Lepton_pt, "Sublead Z Muon p_{T}", {}),
        (z_decay_pdgid==13)&&(N_Leptonic_W==2),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth4l_zsublead_muon_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_HZ_Sublead_Lepton_eta, "Sublead Z Muon #eta", {}),
        (z_decay_pdgid==13)&&(N_Leptonic_W==2),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth4l_zsublead_muon_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_W_Lead_SignalMuon_pt, "W Lead Muon p_{T}", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11)&&(N_Leptonic_W==2),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth4l_w_leadmuon_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_W_Lead_Muon_eta, "W Lead Muon #eta", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11)&&(N_Leptonic_W==2),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth4l_w_leadmuon_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_W_Lead_SignalElectron_pt, "W Lead Electron p_{T}", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11)&&(N_Leptonic_W==2),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth4l_w_leadelectron_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_W_Lead_Electron_eta, "W Lead Electron #eta", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11)&&(N_Leptonic_W==2),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth4l_w_leadelectron_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_W_Sublead_SignalMuon_pt, "W Sublead Muon p_{T}", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11)&&(N_Leptonic_W==2),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth4l_w_subleadmuon_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_W_Sublead_Muon_eta, "W Sublead Muon #eta", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11)&&(N_Leptonic_W==2),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth4l_w_subleadmuon_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, GenPart_W_Sublead_SignalElectron_pt, "W Sublead Electron p_{T}", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11)&&(N_Leptonic_W==2),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth4l_w_subleadelectron_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, GenPart_W_Sublead_Electron_eta, "W Sublead Electron #eta", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11)&&(N_Leptonic_W==2),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth4l_w_subleadelectron_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 120.0, H_Photon_pt, "Photon p_{T}", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11)&&(N_Leptonic_W==2),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth4l_photon_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -4.0, 4.0, H_Photon_eta, "Photon #eta", {}),
        (z_decay_pdgid==13||z_decay_pdgid==11)&&(N_Leptonic_W==2),
        signal_procs_tth, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__tth4l_photon_eta")
        .LuminosityTag(total_luminosity_string);
    //acceptance tables
    pm.Push<Table>("zghlt_ggh_accept", vector<TableRow>{
      TableRow("\\hline $Z\\rightarrow e^{+}e^{-}$ decays", 
          (z_decay_pdgid==11),0,0,weight_noscale),
      TableRow("Electrons in acceptance", 
          (z_decay_pdgid==11)&&(NTruthLeptonAccept>=2),0,0,weight_noscale),
      TableRow("Photon in acceptance", 
          (z_decay_pdgid==11)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1),0,0,weight_noscale),
      TableRow("Electrons reconstructed", 
          (z_decay_pdgid==11)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1)&&(nSignalElectron>=2),0,0,weight_noscale),
      TableRow("Photon reconstructed", 
          (z_decay_pdgid==11)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1)&&(nSignalElectron>=2)&&(nSignalPhoton>=1),0,0,weight_noscale),
      TableRow("L1 Trigger (MC)", 
          (z_decay_pdgid==11)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger,0,0,weight_noscale),
      TableRow("HLT 2l Trigger (MC)", 
          (z_decay_pdgid==11)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&hlt_el_trigger,0,0,weight_noscale),
      TableRow("HLT 2l or 1l Trigger (MC)", 
          (z_decay_pdgid==11)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&hlt_el_trigger_plus,0,0,weight_noscale),
      TableRow("HLT 2l Trigger (Data SFs)", 
          (z_decay_pdgid==11)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&hlt_el_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("HLT 2l or 1l Trigger (Data SFs)", 
          (z_decay_pdgid==11)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&hlt_el_trigger_plus,0,0,weight_noscale*ScaleFactor_Triggers),
      //
      TableRow("\\hline $Z\\rightarrow \\mu^{+}\\mu^{-}$ decays", 
          (z_decay_pdgid==13),0,0,weight_noscale),
      TableRow("Electrons in acceptance", 
          (z_decay_pdgid==13)&&(NTruthLeptonAccept>=2),0,0,weight_noscale),
      TableRow("Photon in acceptance", 
          (z_decay_pdgid==13)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1),0,0,weight_noscale),
      TableRow("Muons reconstructed", 
          (z_decay_pdgid==13)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1)&&(nSignalMuon>=2),0,0,weight_noscale),
      TableRow("Photon reconstructed", 
          (z_decay_pdgid==13)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1)&&(nSignalMuon>=2)&&(nSignalPhoton>=1),0,0,weight_noscale),
      TableRow("L1 Trigger (MC)", 
          (z_decay_pdgid==13)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger,0,0,weight_noscale),
      TableRow("HLT 2l Trigger (MC)", 
          (z_decay_pdgid==13)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&hlt_mu_trigger,0,0,weight_noscale),
      TableRow("HLT 2l or 1l Trigger (MC)", 
          (z_decay_pdgid==13)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&hlt_mu_trigger_plus,0,0,weight_noscale),
      TableRow("HLT 2l Trigger (Data SFs)", 
          (z_decay_pdgid==13)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&hlt_mu_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("HLT 2l or 1l Trigger (Data SFs)", 
          (z_decay_pdgid==13)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&hlt_mu_trigger_plus,0,0,weight_noscale*ScaleFactor_Triggers),
    },signal_procs_noscale,false,true,false,false,false,true).LuminosityTag(total_luminosity_string).Precision(3);
    pm.Push<Table>("zghlt_vbf_accept", vector<TableRow>{
      TableRow("\\hline $Z\\rightarrow e^{+}e^{-}$ decays", 
          (z_decay_pdgid==11),0,0,weight_noscale),
      TableRow("Electrons in acceptance", 
          (z_decay_pdgid==11)&&(NTruthLeptonAccept>=2),0,0,weight_noscale),
      TableRow("Photon in acceptance", 
          (z_decay_pdgid==11)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1),0,0,weight_noscale),
      TableRow("Electrons reconstructed", 
          (z_decay_pdgid==11)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1)&&(nSignalElectron>=2),0,0,weight_noscale),
      TableRow("Photon reconstructed", 
          (z_decay_pdgid==11)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1)&&(nSignalElectron>=2)&&(nSignalPhoton>=1),0,0,weight_noscale),
      TableRow("L1 Trigger (MC)", 
          (z_decay_pdgid==11)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger,0,0,weight_noscale),
      TableRow("HLT 2l Trigger (MC)", 
          (z_decay_pdgid==11)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&hlt_el_trigger,0,0,weight_noscale),
      TableRow("HLT 2l or 1l Trigger (MC)", 
          (z_decay_pdgid==11)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&hlt_el_trigger_plus,0,0,weight_noscale),
      TableRow("HLT 2l Trigger (Data SFs)", 
          (z_decay_pdgid==11)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&hlt_el_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("HLT 2l or 1l Trigger (Data SFs)", 
          (z_decay_pdgid==11)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&hlt_el_trigger_plus,0,0,weight_noscale*ScaleFactor_Triggers),
      //
      TableRow("\\hline $Z\\rightarrow \\mu^{+}\\mu^{-}$ decays", 
          (z_decay_pdgid==13),0,0,weight_noscale),
      TableRow("Electrons in acceptance", 
          (z_decay_pdgid==13)&&(NTruthLeptonAccept>=2),0,0,weight_noscale),
      TableRow("Photon in acceptance", 
          (z_decay_pdgid==13)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1),0,0,weight_noscale),
      TableRow("Muons reconstructed", 
          (z_decay_pdgid==13)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1)&&(nSignalMuon>=2),0,0,weight_noscale),
      TableRow("Photon reconstructed", 
          (z_decay_pdgid==13)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1)&&(nSignalMuon>=2)&&(nSignalPhoton>=1),0,0,weight_noscale),
      TableRow("L1 Trigger (MC)", 
          (z_decay_pdgid==13)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger,0,0,weight_noscale),
      TableRow("HLT 2l Trigger (MC)", 
          (z_decay_pdgid==13)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&hlt_mu_trigger,0,0,weight_noscale),
      TableRow("HLT 2l or 1l Trigger (MC)", 
          (z_decay_pdgid==13)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&hlt_mu_trigger_plus,0,0,weight_noscale),
      TableRow("HLT 2l Trigger (Data SFs)", 
          (z_decay_pdgid==13)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&hlt_mu_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("HLT 2l or 1l Trigger (Data SFs)", 
          (z_decay_pdgid==13)&&(NTruthLeptonAccept>=2)&&(NTruthPhotonAccept>=1)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&hlt_mu_trigger_plus,0,0,weight_noscale*ScaleFactor_Triggers),
    },signal_procs_vbf,false,true,false,false,false,true).LuminosityTag(total_luminosity_string).Precision(3);
    pm.Push<Table>("zghlt_wplush_accept", vector<TableRow>{
      TableRow("\\hline $Z\\rightarrow e^{+}e^{-}$ decays", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1),0,0,weight_noscale),
      TableRow("Electrons in acceptance", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3),0,0,weight_noscale),
      TableRow("Photon in acceptance", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1),0,0,weight_noscale),
      TableRow("Electrons reconstructed", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3),0,0,weight_noscale),
      TableRow("Photon reconstructed", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1),0,0,weight_noscale),
      TableRow("L1 Trigger (MC)", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger,0,0,weight_noscale),
      TableRow("HLT 2l Trigger (MC)", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger,0,0,weight_noscale),
      TableRow("HLT 2l or 1l Trigger (MC)", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger_plus,0,0,weight_noscale),
      TableRow("HLT 2l Trigger (Data SFs)", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("HLT 2l or 1l Trigger (Data SFs)", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger_plus,0,0,weight_noscale*ScaleFactor_Triggers),
      //
      TableRow("\\hline $Z\\rightarrow \\mu^{+}\\mu^{-}$ decays", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1),0,0,weight_noscale),
      TableRow("Muons in acceptance", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3),0,0,weight_noscale),
      TableRow("Photon in acceptance", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1),0,0,weight_noscale),
      TableRow("Muons reconstructed", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3),0,0,weight_noscale),
      TableRow("Photon reconstructed", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1),0,0,weight_noscale),
      TableRow("L1 Trigger (MC)", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger,0,0,weight_noscale),
      TableRow("HLT 2l Trigger (MC)", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger,0,0,weight_noscale),
      TableRow("HLT 2l or 1l Trigger (MC)", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger_plus,0,0,weight_noscale),
      TableRow("HLT 2l Trigger (Data SFs)", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("HLT 2l or 1l Trigger (Data SFs)", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger_plus,0,0,weight_noscale*ScaleFactor_Triggers),
    },signal_procs_wplush,false,true,false,false,false,true).LuminosityTag(total_luminosity_string).Precision(3);
    pm.Push<Table>("zghlt_wminush_accept", vector<TableRow>{
      TableRow("\\hline $Z\\rightarrow e^{+}e^{-}$ decays", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1),0,0,weight_noscale),
      TableRow("Electrons in acceptance", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3),0,0,weight_noscale),
      TableRow("Photon in acceptance", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1),0,0,weight_noscale),
      TableRow("Electrons reconstructed", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3),0,0,weight_noscale),
      TableRow("Photon reconstructed", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1),0,0,weight_noscale),
      TableRow("L1 Trigger (MC)", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger,0,0,weight_noscale),
      TableRow("HLT 2l Trigger (MC)", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger,0,0,weight_noscale),
      TableRow("HLT 2l or 1l Trigger (MC)", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger_plus,0,0,weight_noscale),
      TableRow("HLT 2l Trigger (Data SFs)", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("HLT 2l or 1l Trigger (Data SFs)", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger_plus,0,0,weight_noscale*ScaleFactor_Triggers),
      //
      TableRow("\\hline $Z\\rightarrow \\mu^{+}\\mu^{-}$ decays", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1),0,0,weight_noscale),
      TableRow("Muons in acceptance", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3),0,0,weight_noscale),
      TableRow("Photon in acceptance", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1),0,0,weight_noscale),
      TableRow("Muons reconstructed", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3),0,0,weight_noscale),
      TableRow("Photon reconstructed", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1),0,0,weight_noscale),
      TableRow("L1 Trigger (MC)", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger,0,0,weight_noscale),
      TableRow("HLT 2l Trigger (MC)", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger,0,0,weight_noscale),
      TableRow("HLT 2l or 1l Trigger (MC)", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger_plus,0,0,weight_noscale),
      TableRow("HLT 2l Trigger (Data SFs)", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("HLT 2l or 1l Trigger (Data SFs)", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger_plus,0,0,weight_noscale*ScaleFactor_Triggers),
    },signal_procs_wminush,false,true,false,false,false,true).LuminosityTag(total_luminosity_string).Precision(3);
    pm.Push<Table>("zghlt_zh_accept", vector<TableRow>{
      TableRow("\\hline $Z\\rightarrow e^{+}e^{-}$ decays", 
          (z_decay_pdgid==11)&&(N_Leptonic_Z==2),0,0,weight_noscale),
      TableRow("Electrons in acceptance", 
          (z_decay_pdgid==11)&&(N_Leptonic_Z==2)&&(NTruthLeptonAccept>=4),0,0,weight_noscale),
      TableRow("Photon in acceptance", 
          (z_decay_pdgid==11)&&(N_Leptonic_Z==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1),0,0,weight_noscale),
      TableRow("Electrons reconstructed", 
          (z_decay_pdgid==11)&&(N_Leptonic_Z==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=4),0,0,weight_noscale),
      TableRow("Photon reconstructed", 
          (z_decay_pdgid==11)&&(N_Leptonic_Z==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=4)&&(nSignalPhoton>=1),0,0,weight_noscale),
      TableRow("L1 Trigger (MC)", 
          (z_decay_pdgid==11)&&(N_Leptonic_Z==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=4)&&(nSignalPhoton>=1)&&l1_lep_trigger,0,0,weight_noscale),
      TableRow("HLT 2l Trigger (MC)", 
          (z_decay_pdgid==11)&&(N_Leptonic_Z==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=4)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger,0,0,weight_noscale),
      TableRow("HLT 2l or 1l Trigger (MC)", 
          (z_decay_pdgid==11)&&(N_Leptonic_Z==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=4)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger_plus,0,0,weight_noscale),
      TableRow("HLT 2l Trigger (Data SFs)", 
          (z_decay_pdgid==11)&&(N_Leptonic_Z==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=4)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("HLT 2l or 1l Trigger (Data SFs)", 
          (z_decay_pdgid==11)&&(N_Leptonic_Z==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=4)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger_plus,0,0,weight_noscale*ScaleFactor_Triggers),
      //
      TableRow("\\hline $Z\\rightarrow \\mu^{+}\\mu^{-}$ decays", 
          (z_decay_pdgid==13)&&(N_Leptonic_Z==2),0,0,weight_noscale),
      TableRow("Electrons in acceptance", 
          (z_decay_pdgid==13)&&(N_Leptonic_Z==2)&&(NTruthLeptonAccept>=4),0,0,weight_noscale),
      TableRow("Photon in acceptance", 
          (z_decay_pdgid==13)&&(N_Leptonic_Z==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1),0,0,weight_noscale),
      TableRow("Muons reconstructed", 
          (z_decay_pdgid==13)&&(N_Leptonic_Z==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=4),0,0,weight_noscale),
      TableRow("Photon reconstructed", 
          (z_decay_pdgid==13)&&(N_Leptonic_Z==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=4)&&(nSignalPhoton>=1),0,0,weight_noscale),
      TableRow("L1 Trigger (MC)", 
          (z_decay_pdgid==13)&&(N_Leptonic_Z==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=4)&&(nSignalPhoton>=1)&&l1_lep_trigger,0,0,weight_noscale),
      TableRow("HLT 2l Trigger (MC)", 
          (z_decay_pdgid==13)&&(N_Leptonic_Z==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=4)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger,0,0,weight_noscale),
      TableRow("HLT 2l or 1l Trigger (MC)", 
          (z_decay_pdgid==13)&&(N_Leptonic_Z==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=4)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger_plus,0,0,weight_noscale),
      TableRow("HLT 2l Trigger (Data SFs)", 
          (z_decay_pdgid==13)&&(N_Leptonic_Z==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=4)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("HLT 2l or 1l Trigger (Data SFs)", 
          (z_decay_pdgid==13)&&(N_Leptonic_Z==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=4)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger_plus,0,0,weight_noscale*ScaleFactor_Triggers),
    },signal_procs_zh,false,true,false,false,false,true).LuminosityTag(total_luminosity_string).Precision(3);
    pm.Push<Table>("zghlt_tth3l_accept", vector<TableRow>{
      TableRow("\\hline $Z\\rightarrow e^{+}e^{-}$ decays", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1),0,0,weight_noscale),
      TableRow("Electrons in acceptance", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3),0,0,weight_noscale),
      TableRow("Photon in acceptance", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1),0,0,weight_noscale),
      TableRow("Electrons reconstructed", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3),0,0,weight_noscale),
      TableRow("Photon reconstructed", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1),0,0,weight_noscale),
      TableRow("L1 Trigger (MC)", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger,0,0,weight_noscale),
      TableRow("HLT 2l Trigger (MC)", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger,0,0,weight_noscale),
      TableRow("HLT 2l or 1l Trigger (MC)", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger_plus,0,0,weight_noscale),
      TableRow("HLT 2l Trigger (Data SFs)", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("HLT 2l or 1l Trigger (Data SFs)", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger_plus,0,0,weight_noscale*ScaleFactor_Triggers),
      //
      TableRow("\\hline $Z\\rightarrow \\mu^{+}\\mu^{-}$ decays", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1),0,0,weight_noscale),
      TableRow("Electrons in acceptance", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3),0,0,weight_noscale),
      TableRow("Photon in acceptance", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1),0,0,weight_noscale),
      TableRow("Muons reconstructed", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3),0,0,weight_noscale),
      TableRow("Photon reconstructed", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1),0,0,weight_noscale),
      TableRow("L1 Trigger (MC)", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger,0,0,weight_noscale),
      TableRow("HLT 2l Trigger (MC)", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger,0,0,weight_noscale),
      TableRow("HLT 2l or 1l Trigger (MC)", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger_plus,0,0,weight_noscale),
      TableRow("HLT 2l Trigger (Data SFs)", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("HLT 2l or 1l Trigger (Data SFs)", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==1)&&(NTruthLeptonAccept>=3)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=3)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger_plus,0,0,weight_noscale*ScaleFactor_Triggers),
    },signal_procs_tth,false,true,false,false,false,true).LuminosityTag(total_luminosity_string).Precision(3);
    pm.Push<Table>("zghlt_tth4l_accept", vector<TableRow>{
      TableRow("\\hline $Z\\rightarrow e^{+}e^{-}$ decays", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==2),0,0,weight_noscale),
      TableRow("Electrons in acceptance", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==2)&&(NTruthLeptonAccept>=4),0,0,weight_noscale),
      TableRow("Photon in acceptance", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1),0,0,weight_noscale),
      TableRow("Electrons reconstructed", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=4),0,0,weight_noscale),
      TableRow("Photon reconstructed", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=4)&&(nSignalPhoton>=1),0,0,weight_noscale),
      TableRow("L1 Trigger (MC)", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=4)&&(nSignalPhoton>=1)&&l1_lep_trigger,0,0,weight_noscale),
      TableRow("HLT 2l Trigger (MC)", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=4)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger,0,0,weight_noscale),
      TableRow("HLT 2l or 1l Trigger (MC)", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=4)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger_plus,0,0,weight_noscale),
      TableRow("HLT 2l Trigger (Data SFs)", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=4)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("HLT 2l or 1l Trigger (Data SFs)", 
          (z_decay_pdgid==11)&&(N_Leptonic_W==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=4)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger_plus,0,0,weight_noscale*ScaleFactor_Triggers),
      //
      TableRow("\\hline $Z\\rightarrow \\mu^{+}\\mu^{-}$ decays", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==2),0,0,weight_noscale),
      TableRow("Electrons in acceptance", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==2)&&(NTruthLeptonAccept>=4),0,0,weight_noscale),
      TableRow("Photon in acceptance", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1),0,0,weight_noscale),
      TableRow("Muons reconstructed", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=4),0,0,weight_noscale),
      TableRow("Photon reconstructed", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=4)&&(nSignalPhoton>=1),0,0,weight_noscale),
      TableRow("L1 Trigger (MC)", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=4)&&(nSignalPhoton>=1)&&l1_lep_trigger,0,0,weight_noscale),
      TableRow("HLT 2l Trigger (MC)", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=4)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger,0,0,weight_noscale),
      TableRow("HLT 2l or 1l Trigger (MC)", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=4)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger_plus,0,0,weight_noscale),
      TableRow("HLT 2l Trigger (Data SFs)", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=4)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("HLT 2l or 1l Trigger (Data SFs)", 
          (z_decay_pdgid==13)&&(N_Leptonic_W==2)&&(NTruthLeptonAccept>=4)&&(NTruthPhotonAccept>=1)&&(nSignalLepton>=4)&&(nSignalPhoton>=1)&&l1_lep_trigger&&hlt_lep_trigger_plus,0,0,weight_noscale*ScaleFactor_Triggers),
    },signal_procs_tth,false,true,false,false,false,true).LuminosityTag(total_luminosity_string).Precision(3);

  }

  if (plot_fail_dilep_trig) {
    pm.Push<Hist1D>(Axis(4, -1.5, 2.5, Electron_hltId, "Electron HLT Status", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__elhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(4, -1.5, 2.5, Muon_hltId, "Muon HLT Status", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__muhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(4, -1.5, 2.5, Electron_hltId, "Electron HLT Status", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_elhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(4, -1.5, 2.5, Muon_hltId, "Muon HLT Status", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_muhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 100.0, "Electron_pt", "Electron p_{T} (no HLT match) [GeV]", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId<0.0&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecoel_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -2.5, 2.5, "Electron_eta", "Electron #eta (no HLT match)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId<0.0&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecoel_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -3.16, 3.16, "Electron_phi", "Electron #phi (no HLT match)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId<0.0&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecoel_phi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "Electron_eta", "Electron #eta (no HLT match)", {}),
        Axis(12, -3.16, 3.16, "Electron_phi", "Electron #phi (no HLT match)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId<0.0&&Electron_sig==1.0,
        signal_procs_noscale_2d, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecoel_etaphi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 2.0, "Electron_pfRelIso03_all", "Electron I_{rel,all} (no HLT match)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId<0.0&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecoel_iso")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 2.0, "Electron_pfRelIso03_chg", "Electron I_{rel,chg} (no HLT match)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId<0.0&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecoel_isochg")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 100.0, "Electron_pt", "Electron p_{T} (no HLT ID) [GeV]", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId==0.0&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidel_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -2.5, 2.5, "Electron_eta", "Electron #eta (no HLT ID)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId==0.0&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidel_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -3.16, 3.16, "Electron_phi", "Electron #phi (no HLT ID)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId==0.0&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidel_phi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "Electron_eta", "Electron #eta (no HLT ID)", {}),
        Axis(12, -3.16, 3.16, "Electron_phi", "Electron #phi (no HLT ID)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId==0.0&&Electron_sig==1.0,
        signal_procs_noscale_2d, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidel_etaphi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 2.0, "Electron_pfRelIso03_all", "Electron I_{rel,all} (no HLT ID)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId==0.0&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidel_iso")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 2.0, "Electron_pfRelIso03_chg", "Electron I_{rel,chg} (no HLT ID)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId==0.0&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidel_isochg")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 100.0, "Muon_pt", "Muon p_{T} (no HLT match) [GeV]", {}),
        (l1_mu_trigger||l1_mu_trigger_singlemu)&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId<0.0&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -2.5, 2.5, "Muon_eta", "Muon #eta (no HLT match)", {}),
        (l1_mu_trigger||l1_mu_trigger_singlemu)&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId<0.0&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -3.16, 3.16, "Muon_phi", "Muon #phi (no HLT match)", {}),
        (l1_mu_trigger||l1_mu_trigger_singlemu)&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId<0.0&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_phi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "Muon_eta", "Muon #eta (no HLT match)", {}),
        Axis(12, -3.16, 3.16, "Muon_phi", "Muon #phi (no HLT match)", {}),
        (l1_mu_trigger||l1_mu_trigger_singlemu)&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId<0.0&&Muon_sig==1.0,
        signal_procs_noscale_2d, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_etaphi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 2.0, "Muon_pfRelIso03_all", "Muon I_{rel,all} (no HLT match)", {}),
        (l1_mu_trigger||l1_mu_trigger_singlemu)&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId<0.0&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_iso")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 2.0, "Muon_pfRelIso03_chg", "Muon I_{rel,chg} (no HLT match)", {}),
        (l1_mu_trigger||l1_mu_trigger_singlemu)&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId<0.0&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_isochg")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 100.0, "Muon_pt", "Muon p_{T} (no HLT ID) [GeV]", {}),
        (l1_mu_trigger||l1_mu_trigger_singlemu)&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId==0.0&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidmu_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -2.5, 2.5, "Muon_eta", "Muon #eta (no HLT ID)", {}),
        (l1_mu_trigger||l1_mu_trigger_singlemu)&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId==0.0&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidmu_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -3.16, 3.16, "Muon_phi", "Muon #phi (no HLT ID)", {}),
        (l1_mu_trigger||l1_mu_trigger_singlemu)&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId==0.0&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidmu_phi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "Muon_eta", "Muon #eta (no HLT ID)", {}),
        Axis(12, -3.16, 3.16, "Muon_phi", "Muon #phi (no HLT ID)", {}),
        (l1_mu_trigger||l1_mu_trigger_singlemu)&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId==0.0&&Muon_sig==1.0,
        signal_procs_noscale_2d, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidmu_etaphi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 2.0, "Muon_pfRelIso03_all", "Muon I_{rel,all} (no HLT ID)", {}),
        (l1_mu_trigger||l1_mu_trigger_singlemu)&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId==0.0&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidmu_iso")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 2.0, "Muon_pfRelIso03_chg", "Muon I_{rel,chg} (no HLT ID)", {}),
        (l1_mu_trigger||l1_mu_trigger_singlemu)&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId==0.0&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidmu_isochg")
        .LuminosityTag(total_luminosity_string);
    //muon eta stuff
    pm.Push<Hist1D>(Axis(40, -3.0, 3.0, "Muon_eta", "Muon #eta (no HLT match)", {}),
        (l1_mu_trigger||l1_mu_trigger_singlemu)&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId<0.0&&Muon_sig==1.0&&baseline_selection_nopt,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -3.0, 3.0, "Muon_eta", "Muon #eta (no HLT ID)", {}),
        (l1_mu_trigger||l1_mu_trigger_singlemu)&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId=0.0&&Muon_sig==1.0&&baseline_selection_nopt,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidmu_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -3.0, 3.0, "Muon_eta", "Muon #eta (HLT ID)", {}),
        (l1_mu_trigger||l1_mu_trigger_singlemu)&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId>=1.0&&Muon_sig==1.0&&baseline_selection_nopt,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_idmu_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -3.0, 3.0, "Muon_eta", "Muon #eta (no HLT req, 1#mu and not 2#mu L1 trigger)", {}),
        (!l1_mu_trigger&&l1_mu_trigger_singlemu)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_sig==1.0&&baseline_selection_nopt,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_newl1_mu_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -3.0, 3.0, "Muon_eta", "Muon #eta (no HLT req, 2#mu L1 trigger)", {}),
        (l1_mu_trigger)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_sig==1.0&&baseline_selection_nopt,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_oldl1_mu_eta")
        .LuminosityTag(total_luminosity_string);
    // /muon eta stuff
    // //
  }

  if (plot_more_fail_dilep) {
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "TrigObj_eta", "HLT Muon #eta", {}),
        Axis(12, -3.16, 3.16, "TrigObj_phi", "HLT Muon #phi", {}),
        "TrigObj_id==11" && (nSignalMuon>=2)&&(nSignalPhoton>=1),
        signal_procs_noscale_2d, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_hltmu_etaphi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(12, -3.16, 3.16, "TrigObj_phi", "HLT Muon #phi", {}),
        "TrigObj_id==11" && (nSignalMuon>=2)&&(nSignalPhoton>=1),
        signal_procs_noscale_2d, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_hltmu_phi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "TrigObj_eta", "HLT Muon #eta", {}),
        Axis(12, -3.16, 3.16, "TrigObj_phi", "HLT Muon #phi", {}),
        "TrigObj_id==11",
        signal_procs_noscale_2d, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_hltmu_etaphi_nosel")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "TrigObj_eta", "HLT Muon #eta", {}),
        Axis(12, -3.16, 3.16, "TrigObj_phi", "HLT Muon #phi", {}),
        "TrigObj_id==11" && TrigObj_filterBit2,
        signal_procs_noscale_2d, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_hltmu_etaphi_id")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "Muon_eta", "Offline Muon #eta", {}),
        Axis(12, -3.16, 3.16, "Muon_phi", "Offline Muon #phi", {}),
        (nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_sig,
        signal_procs_noscale_2d, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_offmu_etaphi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "TrigObj_eta", "HLT Muon #eta", {}),
        Axis(12, -3.16, 3.16, "TrigObj_phi", "HLT Muon #phi", {}),
        "TrigObj_id==11",
        dy_procs, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_hltmu_etaphi_nosel_dy")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "Muon_eta", "Offline Muon #eta", {}),
        Axis(12, -3.16, 3.16, "Muon_phi", "Offline Muon #phi", {}),
        (nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_sig,
        dy_procs, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_offmu_etaphi_dy")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "TrigObj_eta", "HLT Muon #eta", {}),
        Axis(12, -3.16, 3.16, "TrigObj_phi", "HLT Muon #phi", {}),
        "TrigObj_id==11",
        dy_procs_2016, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_hltmu_etaphi_nosel_dy16")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "Muon_eta", "Offline Muon #eta", {}),
        Axis(12, -3.16, 3.16, "Muon_phi", "Offline Muon #phi", {}),
        (nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_sig,
        dy_procs_2016, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_offmu_etaphi_dy16")
        .LuminosityTag(total_luminosity_string);
    
    pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "Muon_nStations", "Muon nStations (No HLT Match)", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId<0.0&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_nstations")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(21, -0.5, 20.5, "Muon_nTrackerLayers", "Muon nTrackerLayers (No HLT Match)", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId<0.0&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_ntrackerlayers")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(2, -0.5, 1.5, "Muon_isGlobal", "Muon is global muon (No HLT Match)", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId<0.0&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_isglobal")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "Muon_nStations", "Muon nStations (No HLT Match) (not 0.75 < |#eta| < 1.5)", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId<0.0&&Muon_sig==1.0&&"(Muon_eta>1.5||Muon_eta<-1.5||(Muon_eta>-0.75&&Muon_eta<0.75))",
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_etapeak_nstations")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(21, -0.5, 20.5, "Muon_nTrackerLayers", "Muon nTrackerLayers (No HLT Match) (not 0.75 < |#eta| < 1.5)", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId<0.0&&Muon_sig==1.0&&"(Muon_eta>1.5||Muon_eta<-1.5||(Muon_eta>-0.75&&Muon_eta<0.75))",
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_etapeak_ntrackerlayers")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(2, -0.5, 1.5, "Muon_isGlobal", "Muon is global muon (No HLT Match) (not 0.75 < |#eta| < 1.5)", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId<0.0&&Muon_sig==1.0&&"(Muon_eta>1.5||Muon_eta<-1.5||(Muon_eta>-0.75&&Muon_eta<0.75))",
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_etapeak_isglobal")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "Muon_nStations", "Muon nStations (No HLT Match) (0.75 < |#eta| < 1.5)", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId<0.0&&Muon_sig==1.0&&"((Muon_eta>-1.5&&Muon_eta<-0.75)||(Muon_eta>0.75&&Muon_eta<1.5))",
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_offetapeak_nstations")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(21, -0.5, 20.5, "Muon_nTrackerLayers", "Muon nTrackerLayers (No HLT Match) (0.75 < |#eta| < 1.5)", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId<0.0&&Muon_sig==1.0&&"((Muon_eta>-1.5&&Muon_eta<-0.75)||(Muon_eta>0.75&&Muon_eta<1.5))",
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_offetapeak_ntrackerlayers")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(2, -0.5, 1.5, "Muon_isGlobal", "Muon is global muon (No HLT Match) (0.75 < |#eta| < 1.5)", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId<0.0&&Muon_sig==1.0&&"((Muon_eta>-1.5&&Muon_eta<-0.75)||(Muon_eta>0.75&&Muon_eta<1.5))",
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_offetapeak_isglobal")
        .LuminosityTag(total_luminosity_string);

    pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "Muon_nStations", "Offline Muon nStations", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_offmu_nstations")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(21, -0.5, 20.5, "Muon_nTrackerLayers", "Offline Muon nTrackerLayers", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_offmu_ntrackerlayers")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(2, -0.5, 1.5, "Muon_isGlobal", "Offline Muon is global muon", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_offmu_isglobal")
        .LuminosityTag(total_luminosity_string);

    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "TrigObj_eta", "HLT Electron #eta", {}),
        Axis(12, -3.16, 3.16, "TrigObj_phi", "HLT Electron #phi", {}),
        "TrigObj_id==11",
        signal_procs_noscale_2d, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_hltel_etaphi_nosel")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "TrigObj_eta", "HLT Electron #eta", {}),
        Axis(12, -3.16, 3.16, "TrigObj_phi", "HLT Electron #phi", {}),
        "TrigObj_id==11" && TrigObj_filterBit2,
        signal_procs_noscale_2d, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_hltel_etaphi_id")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "Electron_eta", "Offline Electron #eta", {}),
        Axis(12, -3.16, 3.16, "Electron_phi", "Offline Electron #phi", {}),
        (nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_sig,
        signal_procs_noscale_2d, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_offel_etaphi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "TrigObj_eta", "HLT Electron #eta", {}),
        Axis(12, -3.16, 3.16, "TrigObj_phi", "HLT Electron #phi", {}),
        "TrigObj_id==11",
        dy_procs, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_hltel_etaphi_nosel_dy")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "TrigObj_eta", "HLT Electron #eta", {}),
        Axis(12, -3.16, 3.16, "TrigObj_phi", "HLT Electron #phi", {}),
        "TrigObj_id==11" && TrigObj_filterBit2,
        dy_procs, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_hltel_etaphi_id_dy")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "Electron_eta", "Offline Electron #eta", {}),
        Axis(12, -3.16, 3.16, "Electron_phi", "Offline Electron #phi", {}),
        (nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_sig,
        dy_procs, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_offel_etaphi_dy")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "TrigObj_eta", "HLT Electron #eta", {}),
        Axis(12, -3.16, 3.16, "TrigObj_phi", "HLT Electron #phi", {}),
        "TrigObj_id==11",
        dy_procs_2016, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_hltel_etaphi_nosel_dy16")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "TrigObj_eta", "HLT Electron #eta", {}),
        Axis(12, -3.16, 3.16, "TrigObj_phi", "HLT Electron #phi", {}),
        "TrigObj_id==11" && TrigObj_filterBit2,
        dy_procs_2016, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_hltel_etaphi_id_dy16")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "Electron_eta", "Offline Electron #eta", {}),
        Axis(12, -3.16, 3.16, "Electron_phi", "Offline Electron #phi", {}),
        (nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_sig,
        dy_procs_2016, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_offel_etaphi_dy16")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(8, -0.5, 7.5, "Electron_lostHits", "Electron Lost Hits (No HLT Match)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId<0.0&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecoel_losthits")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, -0.1, 0.1, "Electron_deltaEtaSC", "Electron #Delta#eta_{SC} (No HLT Match)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId<0.0&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecoel_deltaeta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(15, 0.0, 0.03, "Electron_sieie", "Electron #sigma_{i#etai#eta} (No HLT Match)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId<0.0&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecoel_sieie")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 1.0, "Electron_r9", "Electron R_{9} (No HLT Match)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId<0.0&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecoel_r9")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 0.2, "Electron_hoe", "Electron H/E (No HLT Match)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId<0.0&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecoel_hoe")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(8, -0.5, 7.5, "Electron_lostHits", "Electron Lost Hits (No HLT ID)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId==0.0&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidel_losthits")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, -0.1, 0.1, "Electron_deltaEtaSC", "Electron #Delta#eta_{SC} (No HLT ID)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId==0.0&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidel_deltaeta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(15, 0.0, 0.03, "Electron_sieie", "Electron #sigma_{i#etai#eta} (No HLT ID)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId==0.0&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidel_sieie")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 1.0, "Electron_r9", "Electron R_{9} (No HLT ID)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId==0.0&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidel_r9")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 0.2, "Electron_hoe", "Electron H/E (No HLT ID)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId==0.0&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidel_hoe")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(8, -0.5, 7.5, "Electron_lostHits", "Offline Electron Lost Hits", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_offel_losthits")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, -0.1, 0.1, "Electron_deltaEtaSC", "Offline Electron #Delta#eta_{SC}", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_offel_deltaeta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(15, 0.0, 0.03, "Electron_sieie", "Offline Electron #sigma_{i#etai#eta}", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_offel_sieie")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 1.0, "Electron_r9", "Offline Electron R_{9}", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_offel_r9")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 0.2, "Electron_hoe", "Offline Electron H/E", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_offel_hoe")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 50.0, Electron_IsoPhotonPt, "Sum photon p_{T} in electron isolation cone [GeV]", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId==0.0&&Electron_sig==1.0&&"Electron_pfRelIso03_all>0.4",
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidel_isophpt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 50.0, Electron_IsoElectronPt, "Sum electron p_{T} in electron isolation cone [GeV]", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId==0.0&&Electron_sig==1.0&&"Electron_pfRelIso03_all>0.4",
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidel_isoelpt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 50.0, Electron_IsoJetPt, "Sum jet p_{T} in electron isolation cone [GeV]", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId==0.0&&Electron_sig==1.0&&"Electron_pfRelIso03_all>0.4",
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidel_isojetpt")
        .LuminosityTag(total_luminosity_string);
  }

  if (plot_fail_singlelep_trig) {
    pm.Push<Hist1D>(Axis(4, -1.5, 2.5, Electron_hltId, "Electron HLT Status", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_elhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(4, -1.5, 2.5, Muon_hltId, "Muon HLT Status", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_muhlt")
        .LuminosityTag(total_luminosity_string);

    pm.Push<Hist1D>(Axis(20, 0.0, 100.0, "Electron_pt", "Electron p_{T} (no HLT match) [GeV]", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId<0.0&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecoel_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -2.5, 2.5, "Electron_eta", "Electron #eta (no HLT match)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId<0.0&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecoel_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -3.16, 3.16, "Electron_phi", "Electron #phi (no HLT match)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId<0.0&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecoel_phi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "Electron_eta", "Electron #eta (no HLT match)", {}),
        Axis(12, -3.16, 3.16, "Electron_phi", "Electron #phi (no HLT match)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId<0.0&&Electron_sig==1.0,
        signal_procs_noscale_2d, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecoel_etaphi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 2.0, "Electron_pfRelIso03_all", "Electron I_{rel,all} (no HLT match)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId<0.0&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecoel_iso")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(8, -0.5, 7.5, "Electron_lostHits", "Electron Lost Hits (no HLT match)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId<0.0&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecoel_losthits")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, -0.1, 0.1, "Electron_deltaEtaSC", "Electron #Delta#eta_{SC} (no HLT match)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId<0.0&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecoel_deltaeta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(15, 0.0, 0.03, "Electron_sieie", "Electron #sigma_{i#etai#eta} (no HLT match)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId<0.0&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecoel_sieie")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 1.0, "Electron_r9", "Electron R_{9} (no HLT match)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId<0.0&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecoel_r9")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 0.2, "Electron_hoe", "Electron H/E (no HLT match)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_hltId<0.0&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecoel_hoe")
        .LuminosityTag(total_luminosity_string);

    pm.Push<Hist1D>(Axis(20, 0.0, 100.0, "Electron_pt", "Electron p_{T} (no HLT ID) [GeV]", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&(Electron_hltId==0.0||Electron_hltId==1.0)&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidel_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -2.5, 2.5, "Electron_eta", "Electron #eta (no HLT ID)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&(Electron_hltId==0.0||Electron_hltId==1.0)&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidel_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -3.16, 3.16, "Electron_phi", "Electron #phi (no HLT ID)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&(Electron_hltId==0.0||Electron_hltId==1.0)&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidel_phi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "Electron_eta", "Electron #eta (no HLT ID)", {}),
        Axis(12, -3.16, 3.16, "Electron_phi", "Electron #phi (no HLT ID)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&(Electron_hltId==0.0||Electron_hltId==1.0)&&Electron_sig==1.0,
        signal_procs_noscale_2d, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidel_etaphi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 2.0, "Electron_pfRelIso03_all", "Electron I_{rel,all} (no HLT ID)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&(Electron_hltId==0.0||Electron_hltId==1.0)&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidel_iso")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(8, -0.5, 7.5, "Electron_lostHits", "Electron Lost Hits (no HLT ID)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&(Electron_hltId==0.0||Electron_hltId==1.0)&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidel_losthits")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, -0.1, 0.1, "Electron_deltaEtaSC", "Electron #Delta#eta_{SC} (no HLT ID)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&(Electron_hltId==0.0||Electron_hltId==1.0)&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidel_deltaeta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(15, 0.0, 0.03, "Electron_sieie", "Electron #sigma_{i#etai#eta} (no HLT ID)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&(Electron_hltId==0.0||Electron_hltId==1.0)&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidel_sieie")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 1.0, "Electron_r9", "Electron R_{9} (no HLT ID)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&(Electron_hltId==0.0||Electron_hltId==1.0)&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidel_r9")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 0.2, "Electron_hoe", "Electron H/E (no HLT ID)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&(Electron_hltId==0.0||Electron_hltId==1.0)&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidel_hoe")
        .LuminosityTag(total_luminosity_string);

    pm.Push<Hist1D>(Axis(8, -0.5, 7.5, "Electron_lostHits", "Offline Electron Lost Hits", {}),
        (nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_offel_lostHits")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, -0.1, 0.1, "Electron_deltaEtaSC", "Offline Electron #Delta#eta_{SC}", {}),
        (nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_offel_deltaeta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(15, 0.0, 0.03, "Electron_sieie", "Offline Electron #sigma_{i#etai#eta}", {}),
        (nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_offel_sieie")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 1.0, "Electron_r9", "Offline Electron R_{9}", {}),
        (nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_offel_r9")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 0.2, "Electron_hoe", "Offline Electron H/E", {}),
        (nSignalElectron>=2)&&(nSignalPhoton>=1)&&Electron_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_offel_hoe")
        .LuminosityTag(total_luminosity_string);

    pm.Push<Hist1D>(Axis(20, 0.0, 100.0, "Muon_pt", "Muon p_{T} (no HLT match) [GeV]", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId<0.0&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecomu_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -2.5, 2.5, "Muon_eta", "Muon #eta (no HLT match)", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId<0.0&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecomu_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -3.16, 3.16, "Muon_phi", "Muon #phi (no HLT match)", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId<0.0&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecomu_phi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "Muon_eta", "Muon #eta (no HLT match)", {}),
        Axis(12, -3.16, 3.16, "Muon_phi", "Muon #phi (no HLT match)", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId<0.0&&Muon_sig==1.0,
        signal_procs_noscale_2d, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecomu_etaphi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 2.0, "Muon_pfRelIso03_all", "Muon I_{rel,all} (no HLT match)", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId<0.0&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecomu_iso")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "Muon_nStations", "Muon nStations (no HLT match)", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId<0.0&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecomu_nstations")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(21, -0.5, 20.5, "Muon_nTrackerLayers", "Muon nTrackerLayers (no HLT match)", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId<0.0&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecomu_ntrackerlayers")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(2, -0.5, 1.5, "Muon_isGlobal", "Muon is global muon (no HLT match)", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_hltId<0.0&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecomu_isglobal")
        .LuminosityTag(total_luminosity_string);

    pm.Push<Hist1D>(Axis(20, 0.0, 100.0, "Muon_pt", "Muon p_{T} (no HLT ID) [GeV]", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&(Muon_hltId==0.0||Muon_hltId==1.0)&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidmu_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -2.5, 2.5, "Muon_eta", "Muon #eta (no HLT ID)", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&(Muon_hltId==0.0||Muon_hltId==1.0)&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidmu_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -3.16, 3.16, "Muon_phi", "Muon #phi (no HLT ID)", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&(Muon_hltId==0.0||Muon_hltId==1.0)&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidmu_phi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "Muon_eta", "Muon #eta (no HLT ID)", {}),
        Axis(12, -3.16, 3.16, "Muon_phi", "Muon #phi (no HLT ID)", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&(Muon_hltId==0.0||Muon_hltId==1.0)&&Muon_sig==1.0,
        signal_procs_noscale_2d, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidmu_etaphi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 2.0, "Muon_pfRelIso03_all", "Muon I_{rel,all} (no HLT ID)", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&(Muon_hltId==0.0||Muon_hltId==1.0)&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidmu_iso")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "Muon_nStations", "Muon nStations (no HLT ID)", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&(Muon_hltId==0.0||Muon_hltId==1.0)&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidmu_nstations")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(21, -0.5, 20.5, "Muon_nTrackerLayers", "Muon nTrackerLayers (no HLT ID)", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&(Muon_hltId==0.0||Muon_hltId==1.0)&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidmu_ntrackerlayers")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(2, -0.5, 1.5, "Muon_isGlobal", "Muon is global muon (no HLT ID)", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&(Muon_hltId==0.0||Muon_hltId==1.0)&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidmu_isglobal")
        .LuminosityTag(total_luminosity_string);

    pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "Muon_nStations", "Offline Muon nStations", {}),
        (nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_offmu_nstations")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(21, -0.5, 20.5, "Muon_nTrackerLayers", "Offline Muon nTrackerLayers", {}),
        (nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_offmu_ntrackerlayers")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(2, -0.5, 1.5, "Muon_isGlobal", "Offline Muon is global muon", {}),
        (nSignalMuon>=2)&&(nSignalPhoton>=1)&&Muon_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_offmu_isglobal")
        .LuminosityTag(total_luminosity_string);
  }

  if (plot_diphoton_trig_effects) {
    pm.Push<Hist1D>(Axis(40, 0.0, 160.0, LeptonPhoton_mass, "m_{e#gamma} [GeV]", {}),
        (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,
        signal_procs_noscale, plt_lin_over).Weight(weight*ScaleFactor_Triggers)
        .Tag("FixName:zghlt__diph_meg")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, 0.0, 160.0, ZCandidate_mass, "m_{ee} [GeV]", {}),
        (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,
        signal_procs_noscale, plt_lin_over).Weight(weight*ScaleFactor_Triggers)
        .Tag("FixName:zghlt__diph_mee")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, 0.0, 160.0, Max_LeptonPhoton_mass, "Max m_{e#gamma} [GeV]", {}),
        (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,
        signal_procs_noscale, plt_lin_over).Weight(weight*ScaleFactor_Triggers)
        .Tag("FixName:zghlt__diph_maxmeg")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, 0.0, 160.0, Max_LeptonPhoton_mass, "Max m_{e#gamma} [GeV]", {}),
        (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!hlt_el_trigger_plus&&"HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90",
        signal_procs_noscale, plt_lin_over).Weight(weight*ScaleFactor_Triggers)
        .Tag("FixName:zghlt__diph_maxmeg_addedbydiph90")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, 0.0, 160.0, ZCandidate_mass, "m_{ee} [GeV]", {}),
        (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!hlt_el_trigger_plus&&"HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90",
        signal_procs_noscale, plt_lin_over).Weight(weight*ScaleFactor_Triggers)
        .Tag("FixName:zghlt__diph_mee_addedbydiph90")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(15, 110.0, 140.0, HiggsCandidate_mass, "m_{ee#gamma} [GeV]", {}),
        (nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,
        signal_diph_trigger_procs, plt_shapes).Weight(weight*ScaleFactor_Triggers)
        .Tag("FixName:zghlt__mllg_signal_shape_trigs_baseline_el")
        .LuminosityTag(total_luminosity_string);
  }

  if (plot_cleanmask) {
    pm.Push<Hist1D>(Axis(16, 100.0, 180.0, HiggsCandidate_mass, "m_{ee#gamma} [GeV]", {}),
        stitch_dy&&hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nodr,
        full_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__clean_mllg_diel")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 100.0, 180.0, HiggsCandidate_mass, "m_{#mu#mu#gamma} [GeV]", {}),
        stitch_dy&&hlt_mu_trigger&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nodr,
        full_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__clean_mllg_dimu")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 100.0, 180.0, HiggsCandidate_mass_cleaned, "m_{ee#gamma} [GeV]", {}),
        stitch_dy&&hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nodr,
        full_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__clean_mllgclean_diel")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 100.0, 180.0, HiggsCandidate_mass_cleaned, "m_{#mu#mu#gamma} [GeV]", {}),
        stitch_dy&&hlt_mu_trigger&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nodr,
        full_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__clean_mllgclean_dimu")
        .LuminosityTag(total_luminosity_string);
  }

  //NamedFunc baseline_selection = ((ZCandidate_mass>50)&&((Lead_SignalElectron_pt>25)||(Lead_SignalMuon_pt>20))&&(Lead_SignalPhoton_pt/HiggsCandidate_mass>0.14)&&delta_r_cut&&((HiggsCandidate_mass+ZCandidate_mass)>185));
  NamedFunc higgs_window = (HiggsCandidate_mass>120)&&(HiggsCandidate_mass<130);
  if (plot_ptg) {
    pm.Push<Hist1D>(Axis(15, 15.0, 60.0, Lead_SignalPhoton_pt, "p_{T#gamma} [GeV]", {}),
        stitch_dy&&hlt_el_trigger_plus&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection&&higgs_window,
        full_procs, plt_lin).Weight(weight*138.0/41.0)
        .Tag("FixName:zghlt__clean_mllg_diel")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(15, 15.0, 60.0, Lead_SignalPhoton_pt, "p_{T#gamma} [GeV]", {}),
        stitch_dy&&hlt_mu_trigger_plus&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection&&higgs_window,
        full_procs, plt_lin).Weight(weight*138.0/41.0)
        .Tag("FixName:zghlt__clean_mllg_dimu")
        .LuminosityTag(total_luminosity_string);
  }

  if (plot_jaebak_cutflow) {
    pm.Push<Table>("zgtrig_sync_cutflow_"+options.year_string, vector<TableRow>{
      TableRow("Total events", 
          "1",0,0,weight),
      TableRow("Overlap removal (old stitch)", 
          stitch_dy,0,0,weight),
      TableRow("\\hline Dielectron Trigger", 
          stitch_dy&&hlt_el_trigger,0,0,weight),
      TableRow(">=2 Signal Electrons", 
          stitch_dy&&hlt_el_trigger&&(nSignalElectron==2),0,0,weight),
      TableRow(">=1 Signal Photon", 
          stitch_dy&&hlt_el_trigger&&(nSignalElectron==2)&&(nSignalPhoton>=1),0,0,weight),
      TableRow("Z mass > 50 GeV", 
          stitch_dy&&hlt_el_trigger&&(nSignalElectron==2)&&(nSignalPhoton>=1)&&(ZCandidate_mass>50),0,0,weight),
      TableRow("Lead electron pT > 25 GeV", 
          stitch_dy&&hlt_el_trigger&&(nSignalElectron==2)&&(nSignalPhoton>=1)&&(ZCandidate_mass>50)&&(Lead_SignalElectron_pt>25),0,0,weight),
      TableRow("Lead photon pT/Higgs mass > 0.14", 
          stitch_dy&&hlt_el_trigger&&(nSignalElectron==2)&&(nSignalPhoton>=1)&&(ZCandidate_mass>50)&&(Lead_SignalElectron_pt>25)&&(Lead_SignalPhoton_pt/HiggsCandidate_mass>0.14),0,0,weight),
      TableRow("Electron-photon Delta R > 0.4", 
          stitch_dy&&hlt_el_trigger&&(nSignalElectron==2)&&(nSignalPhoton>=1)&&(ZCandidate_mass>50)&&(Lead_SignalElectron_pt>25)&&(Lead_SignalPhoton_pt/HiggsCandidate_mass>0.14)&&delta_r_cut,0,0,weight),
      TableRow("H mass + Z mass > 185 GeV", 
          stitch_dy&&hlt_el_trigger&&(nSignalElectron==2)&&(nSignalPhoton>=1)&&(ZCandidate_mass>50)&&(Lead_SignalElectron_pt>25)&&(Lead_SignalPhoton_pt/HiggsCandidate_mass>0.14)&&delta_r_cut&&((HiggsCandidate_mass+ZCandidate_mass)>185),0,0,weight),
      TableRow("\\hline Dimuon Trigger", 
          stitch_dy&&hlt_mu_trigger,0,0,weight),
      TableRow(">=2 Signal Muons", 
          stitch_dy&&hlt_mu_trigger&&(nSignalMuon==2),0,0,weight),
      TableRow(">=1 Signal Photon", 
          stitch_dy&&hlt_mu_trigger&&(nSignalMuon==2)&&(nSignalPhoton>=1),0,0,weight),
      TableRow("Z mass > 50 GeV", 
          stitch_dy&&hlt_mu_trigger&&(nSignalMuon==2)&&(nSignalPhoton>=1)&&(ZCandidate_mass>50),0,0,weight),
      TableRow("Lead muon pT > 20 GeV", 
          stitch_dy&&hlt_mu_trigger&&(nSignalMuon==2)&&(nSignalPhoton>=1)&&(ZCandidate_mass>50)&&(Lead_SignalMuon_pt>20),0,0,weight),
      TableRow("Lead photon pT/Higgs mass > 0.14", 
          stitch_dy&&hlt_mu_trigger&&(nSignalMuon==2)&&(nSignalPhoton>=1)&&(ZCandidate_mass>50)&&(Lead_SignalMuon_pt>20)&&(Lead_SignalPhoton_pt/HiggsCandidate_mass>0.14),0,0,weight),
      TableRow("Muon-photon Delta R > 0.4", 
          stitch_dy&&hlt_mu_trigger&&(nSignalMuon==2)&&(nSignalPhoton>=1)&&(ZCandidate_mass>50)&&(Lead_SignalMuon_pt>20)&&(Lead_SignalPhoton_pt/HiggsCandidate_mass>0.14)&&delta_r_cut,0,0,weight),
      TableRow("H mass + Z mass > 185 GeV", 
          stitch_dy&&hlt_mu_trigger&&(nSignalMuon==2)&&(nSignalPhoton>=1)&&(ZCandidate_mass>50)&&(Lead_SignalMuon_pt>20)&&(Lead_SignalPhoton_pt/HiggsCandidate_mass>0.14)&&delta_r_cut&&((HiggsCandidate_mass+ZCandidate_mass)>185),0,0,weight),
    },full_procs,false,true,false,false,false,true).LuminosityTag(total_luminosity_string).Precision(2);
  }

  if (plot_single_lep_trig_effects) {
    pm.Push<Hist1D>(Axis(16, 100.0, 180.0, HiggsCandidate_mass, "m_{ee#gamma} [GeV]", {}),
        stitch_dy&&hlt_el_trigger&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection,
        full_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__mllg_diel")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 100.0, 180.0, HiggsCandidate_mass, "m_{#mu#mu#gamma} [GeV]", {}),
        stitch_dy&&hlt_mu_trigger&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection,
        full_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__mllg_dimu")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 100.0, 180.0, HiggsCandidate_mass, "m_{ee#gamma} [GeV]", {}),
        stitch_dy&&(hlt_el_trigger||hlt_single_el_trigger)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection,
        full_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__mllg_eldiel")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 100.0, 180.0, HiggsCandidate_mass, "m_{#mu#mu#gamma} [GeV]", {}),
        stitch_dy&&(hlt_mu_trigger||hlt_single_mu_trigger)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection,
        full_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__mllg_mudimu")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 100.0, 180.0, HiggsCandidate_mass, "m_{ee#gamma} [GeV]", {}),
        stitch_dy&&(!hlt_el_trigger&&hlt_single_el_trigger)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection,
        full_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__mllg_onlyel")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 100.0, 180.0, HiggsCandidate_mass, "m_{#mu#mu#gamma} [GeV]", {}),
        stitch_dy&&(!hlt_mu_trigger&&hlt_single_mu_trigger)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection,
        full_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__mllg_onlymu")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 100.0, 180.0, HiggsCandidate_mass, "m_{ee#gamma} [GeV]", {}),
        stitch_dy&&(hlt_el_trigger&&!hlt_single_el_trigger)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection,
        full_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__mllg_onlydiel")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 100.0, 180.0, HiggsCandidate_mass, "m_{#mu#mu#gamma} [GeV]", {}),
        stitch_dy&&(hlt_mu_trigger&&!hlt_single_mu_trigger)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection,
        full_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__mllg_onlydimu")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 100.0, 180.0, HiggsCandidate_mass, "m_{ee#gamma} [GeV]", {}),
        stitch_dy&&(hlt_el_trigger&&hlt_single_el_trigger)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection,
        full_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__mllg_botheldiel")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 100.0, 180.0, HiggsCandidate_mass, "m_{#mu#mu#gamma} [GeV]", {}),
        stitch_dy&&(hlt_mu_trigger&&hlt_single_mu_trigger)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection,
        full_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__mllg_bothmudimu")
        .LuminosityTag(total_luminosity_string);
  }

  if (plot_single_lep_trig_effects_sig) {
    //single trigger plots
    //pm.Push<Hist1D>(Axis(45, 50.0, 140.0, HiggsCandidate_mass, "m_{ee#gamma} [GeV]", {}),
    //    (nSignalElectron>=2)&&(nSignalPhoton>=1),
    //    signal_trigger_procs, plt_shapes).Weight(weight)
    //    .Tag("FixName:zghlt__mllg_signal_shape_trigs_el")
    //    .LuminosityTag(total_luminosity_string);
    //pm.Push<Hist1D>(Axis(15, 110.0, 140.0, HiggsCandidate_mass, "m_{#mu#mu#gamma} [GeV]", {}),
    //    (nSignalMuon>=2)&&(nSignalPhoton>=1),
    //    signal_trigger_procs, plt_shapes).Weight(weight)
    //    .Tag("FixName:zghlt__mllg_signal_shape_trigs_mu")
    //    .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(15, 110.0, 140.0, HiggsCandidate_mass, "m_{ee#gamma} [GeV]", {}),
        (nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,
        signal_trigger_procs, plt_shapes).Weight(weight*ScaleFactor_Triggers)
        .Tag("FixName:zghlt__mllg_signal_shape_trigs_baseline_el")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(15, 110.0, 140.0, HiggsCandidate_mass, "m_{#mu#mu#gamma} [GeV]", {}),
        (nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,
        signal_trigger_procs, plt_shapes).Weight(weight*ScaleFactor_Triggers)
        .Tag("FixName:zghlt__mllg_signal_shape_trigs_baseline_mu")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -0.5, 1.0, Lead_Photon_mvaID, "Photon ID MVA score", {}),
        (nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,
        signal_trigger_procs, plt_shapes).Weight(weight*ScaleFactor_Triggers)
        .Tag("FixName:zghlt__phmva_signal_shape_trigs_baseline_el")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -0.5, 1.0, Lead_Photon_mvaID, "Photon ID MVA score", {}),
        (nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,
        signal_trigger_procs, plt_shapes).Weight(weight*ScaleFactor_Triggers)
        .Tag("FixName:zghlt__phmva_signal_shape_trigs_baseline_mu")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 4.0, Min_Lepton_Photon_deltaR, "Min Lepton-Photon #DeltaR", {}),
        (nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,
        signal_trigger_procs, plt_shapes).Weight(weight*ScaleFactor_Triggers)
        .Tag("FixName:zghlt__deltar_signal_shape_trigs_baseline_el")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 4.0, Min_Lepton_Photon_deltaR, "Min Lepton-Photon #DeltaR", {}),
        (nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,
        signal_trigger_procs, plt_shapes).Weight(weight*ScaleFactor_Triggers)
        .Tag("FixName:zghlt__deltar_signal_shape_trigs_baseline_mu")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 2.5, HiggsCandidate_ptOverMass, "p_{Tll#gamma}/m_{ll#gamma}", {}),
        (nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,
        signal_trigger_procs, plt_shapes).Weight(weight*ScaleFactor_Triggers)
        .Tag("FixName:zghlt__ptovermllg_signal_shape_trigs_baseline_el")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 2.5, HiggsCandidate_ptOverMass, "p_{Tll#gamma}/m_{ll#gamma}", {}),
        (nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,
        signal_trigger_procs, plt_shapes).Weight(weight*ScaleFactor_Triggers)
        .Tag("FixName:zghlt__ptovermllg_signal_shape_trigs_baseline_mu")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -3.0, 3.0, "Electron_eta", "Electron #eta", {}),
        (nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&Electron_sig,
        signal_trigger_procs, plt_shapes).Weight(weight*ScaleFactor_Triggers)
        .Tag("FixName:zghlt__eleta_signal_shape_trigs_baseline_el")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -3.0, 3.0, "Muon_eta", "Muon #eta", {}),
        (nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&Muon_sig,
        signal_trigger_procs, plt_shapes).Weight(weight*ScaleFactor_Triggers)
        .Tag("FixName:zghlt__mueta_signal_shape_trigs_baseline_mu")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 100.0, "Electron_pt", "Electron p_{T} [GeV]", {}),
        (nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&Electron_sig,
        signal_trigger_procs, plt_shapes).Weight(weight*ScaleFactor_Triggers)
        .Tag("FixName:zghlt__elpt_signal_shape_trigs_baseline_el")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 100.0, "Muon_pt", "Muon p_{T} [GeV]", {}),
        (nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&Muon_sig,
        signal_trigger_procs, plt_shapes).Weight(weight*ScaleFactor_Triggers)
        .Tag("FixName:zghlt__mupt_signal_shape_trigs_baseline_mu")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -3.0, 3.0, "Photon_eta", "Photon #eta", {}),
        (nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&Photon_sig,
        signal_trigger_procs, plt_shapes).Weight(weight*ScaleFactor_Triggers)
        .Tag("FixName:zghlt__pheta_signal_shape_trigs_baseline_el")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -3.0, 3.0, "Photon_eta", "Photon #eta", {}),
        (nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&Photon_sig,
        signal_trigger_procs, plt_shapes).Weight(weight*ScaleFactor_Triggers)
        .Tag("FixName:zghlt__pheta_signal_shape_trigs_baseline_mu")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 100.0, "Photon_pt", "Photon p_{T} [GeV]", {}),
        (nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&Photon_sig,
        signal_trigger_procs, plt_shapes).Weight(weight*ScaleFactor_Triggers)
        .Tag("FixName:zghlt__phpt_signal_shape_trigs_baseline_el")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 100.0, "Photon_pt", "Photon p_{T} [GeV]", {}),
        (nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&Photon_sig,
        signal_trigger_procs, plt_shapes).Weight(weight*ScaleFactor_Triggers)
        .Tag("FixName:zghlt__phpt_signal_shape_trigs_baseline_mu")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(20, 50.0, 150.0, HiggsCandidate_mass, "m_{ee#gamma} [GeV]", {}),
        Axis(20, 0.0, 4.0, ElectronGenPhotonDeltaR, "#Delta R_{e#gamma}^{min}", {}),
        stitch_dy&&(!hlt_el_trigger&&hlt_single_el_trigger)&&(nSignalElectron>=2)&&(nSignalPhoton>=1),
        signal_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__mllg_deltaregamma")
        .LuminosityTag(total_luminosity_string);
  }

  if (make_rereco_plots) {
    pm.Push<Hist1D>(Axis(20, 0.0, 100.0, hlt_max_el_pt, "HLT Electron p_{T}", {}),
        (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&!hlt_el_trigger_plus,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__hltelpt__failhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(27, 0.0, 90.0, "Electron_pt", "Electron p_{T}", {}),
        (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&hlt_el_trigger&&Electron_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__elpt__passhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(27, 0.0, 90.0, "Electron_pt", "Electron p_{T}", {}),
        (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&!hlt_el_trigger&&Electron_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__elpt__failhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, -2.5, 2.5, "Electron_eta", "Electron #eta", {}),
        (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&hlt_el_trigger&&Electron_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__eleta__passhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, -2.5, 2.5, "Electron_eta", "Electron #eta", {}),
        (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&!hlt_el_trigger&&Electron_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__eleta__failhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, -2.5, 2.5, "Electron_eta", "Electron #eta", {}),
        (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==1)&&Electron_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__eleta__failhltreco")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 2.0, "Electron_pfRelIso03_all", "Electron I_{rel}", {}),
        (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&hlt_el_trigger&&Electron_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__elreliso__passhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 2.0, "Electron_pfRelIso03_all", "Electron I_{rel}", {}),
        (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&!hlt_el_trigger&&Electron_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__elreliso__failhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 2.0, "Electron_pfRelIso03_all", "Electron I_{rel}", {}),
        (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==2)&&Electron_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__elreliso__failhltiso")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 2.0, "Electron_pfRelIso03_chg", "Electron I_{rel,chg}", {}),
        (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&hlt_el_trigger&&Electron_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__elrelisochg__passhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 2.0, "Electron_pfRelIso03_chg", "Electron I_{rel,chg}", {}),
        (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&!hlt_el_trigger&&Electron_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__elrelisochg__failhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 2.0, "Muon_pfRelIso03_all", "Muon I_{rel}", {}),
        (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&hlt_mu_trigger&&Muon_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__mureliso__passhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 2.0, "Muon_pfRelIso03_all", "Muon I_{rel}", {}),
        (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&!hlt_mu_trigger&&Muon_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__mureliso__failhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 2.0, "Muon_pfRelIso03_all", "Muon I_{rel}", {}),
        (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==2)&&Muon_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__mureliso__failhltiso")
        .LuminosityTag(total_luminosity_string);


    pm.Push<Hist2D>(Axis(15, -2.5, 2.5, "Electron_eta", "Electron #eta", {}), Axis(10, -3.14, 3.14, "Electron_phi", "Electron #phi", {}),
        (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&hlt_el_trigger&&Electron_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__eletaphi__passhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(15, -2.5, 2.5, "Electron_eta", "Electron #eta", {}), Axis(10, -3.14, 3.14, "Electron_phi", "Electron #phi", {}),
        (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&!hlt_el_trigger&&Electron_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__eletaphi__failhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(15, -2.5, 2.5, "Electron_eta", "Electron #eta", {}), Axis(10, -3.14, 3.14, "Electron_phi", "Electron #phi", {}),
        (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==1)&&Electron_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__eletaphi__failhltreco")
        .LuminosityTag(total_luminosity_string);
  }
  
  if (make_table_fail_reason) {
    pm.Push<Table>("zghlt__2l_fail_reason"+options.year_string, vector<TableRow>{
      TableRow("\\hline \\hline\n $Z\\rightarrow e^{+}e^{-}$ decays", 
          (z_decay_pdgid==11),0,0,weight_noscale),
      TableRow("Offline leptons and photon in acceptance", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1),0,0,weight_noscale),
      TableRow("Baseline Selection (No pT cut)", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,0,0,weight_noscale),
      TableRow("L1 Trigger", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger,0,0,weight_noscale),
      TableRow("Dilepton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&hlt_el_trigger,0,0,weight_noscale),
      TableRow("Dilepton Triggers (Data SFs)", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&hlt_el_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Fail Dilepton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&!hlt_el_trigger,0,0,weight_noscale),
      TableRow("Fail Dilepton Triggers (Data SFs)", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&!hlt_el_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Fail Dilepton Triggers - 1 Trigger Object", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==1),0,0,weight_noscale),
      TableRow("Fail Dilepton Triggers (Data SFs) - 1 Trigger Object", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==1),0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Fail Dilepton Triggers - Failed ID/iso", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==2),0,0,weight_noscale),
      TableRow("Fail Dilepton Triggers (Data SFs) - Failed ID/iso", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==2),0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Fail Dilepton Triggers - Failed top leg pT cut", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==3),0,0,weight_noscale),
      TableRow("Fail Dilepton Triggers (Data SFs) - Failed top leg pT cut", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==3),0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Fail Dilepton Triggers - Failed bottom leg pT cut", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==4),0,0,weight_noscale),
      TableRow("Fail Dilepton Triggers (Data SFs) - Failed bottom leg pT cut", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==4),0,0,weight_noscale*ScaleFactor_Triggers),
      //TableRow("\\hline Fail Single and Dilepton Triggers (Data SFs)", 
      //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&!hlt_el_trigger_plus,0,0,weight_noscale*ScaleFactor_Triggers),
      //TableRow("Fail Single and Dilepton Triggers (Data SFs) 0e0ph", 
      //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&!hlt_el_trigger_plus&&(NElectron_trig==0.0&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale*ScaleFactor_Triggers),
      //TableRow("Fail Single and Dilepton Triggers (Data SFs) 0e1ph", 
      //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&!hlt_el_trigger_plus&&(NElectron_trig==0.0&&"(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale*ScaleFactor_Triggers),
      //TableRow("Fail Single and Dilepton Triggers (Data SFs) 1e0ph", 
      //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&!hlt_el_trigger_plus&&(NElectron_trig==1&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale*ScaleFactor_Triggers),
      //TableRow("Fail Single and Dilepton Triggers (Data SFs) 1e1ph", 
      //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&!hlt_el_trigger_plus&&(NElectron_trig==1&&"(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale*ScaleFactor_Triggers),
      //TableRow("Fail Single and Dilepton Triggers (Data SFs) 2e0ph", 
      //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&!hlt_el_trigger_plus&&(NElectron_trig==2&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale*ScaleFactor_Triggers),
      //TableRow("Fail Single and Dilepton Triggers (Data SFs) 2e1ph", 
      //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&!hlt_el_trigger_plus&&(NElectron_trig==2&&"(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("\\hline Fail Single and Dilepton Triggers (No L1) (Data SFs)", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!hlt_el_trigger_plus,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Fail Single and Dilepton Triggers (Data SFs) 0e0ph", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!hlt_el_trigger_plus&&(NElectron_trig==0.0&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Fail Single and Dilepton Triggers (Data SFs) 0e1ph", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!hlt_el_trigger_plus&&(NElectron_trig==0.0&&"(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Fail Single and Dilepton Triggers (Data SFs) 1e0ph", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!hlt_el_trigger_plus&&(NElectron_trig==1&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Fail Single and Dilepton Triggers (Data SFs) 1e1ph", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!hlt_el_trigger_plus&&(NElectron_trig==1&&"(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Fail Single and Dilepton Triggers (Data SFs) 2e0ph", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!hlt_el_trigger_plus&&(NElectron_trig==2&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Fail Single and Dilepton Triggers (Data SFs) 2e1ph", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!hlt_el_trigger_plus&&(NElectron_trig==2&&"(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("\\hline Fail e, ee, phph Triggers (No L1) (Data SFs)", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!(hlt_el_trigger_plus||hlt_doubleph_trigger),0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Fail e, ee, phph Triggers (Data SFs) 0e0ph", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!(hlt_el_trigger_plus||hlt_doubleph_trigger)&&(NElectron_trig==0.0&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Fail e, ee, phph Triggers (Data SFs) 0e1ph", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!(hlt_el_trigger_plus||hlt_doubleph_trigger)&&(NElectron_trig==0.0&&"(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Fail e, ee, phph Triggers (Data SFs) 1e0ph", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!(hlt_el_trigger_plus||hlt_doubleph_trigger)&&(NElectron_trig==1&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Fail e, ee, phph Triggers (Data SFs) 1e1ph", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!(hlt_el_trigger_plus||hlt_doubleph_trigger)&&(NElectron_trig==1&&"(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Fail e, ee, phph Triggers (Data SFs) 2e0ph", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!(hlt_el_trigger_plus||hlt_doubleph_trigger)&&(NElectron_trig==2&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Fail e, ee, phph Triggers (Data SFs) 2e1ph", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!(hlt_el_trigger_plus||hlt_doubleph_trigger)&&(NElectron_trig==2&&"(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale*ScaleFactor_Triggers),

      TableRow("\\hline \\hline\n $Z\\rightarrow \\mu^{+}\\mu^{-}$ decays", 
          (z_decay_pdgid==13),0,0,weight_noscale),
      TableRow("Offline leptons and photon in acceptance", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1),0,0,weight_noscale),
      TableRow("Baseline Selection (No pT cut)", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,0,0,weight_noscale),
      TableRow("L1 Trigger (only 2l)", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger),0,0,weight_noscale),
      TableRow("Dilepton Triggers", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&hlt_mu_trigger,0,0,weight_noscale),
      TableRow("Dilepton Triggers (Sata SFs)", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&hlt_mu_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Fail Dilepton Triggers", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_mu_trigger&&!hlt_mu_trigger,0,0,weight_noscale),
      TableRow("Fail Dilepton Triggers (Sata SFs)", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_mu_trigger&&!hlt_mu_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Fail Dilepton Triggers - 1 Trigger Object", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==1),0,0,weight_noscale),
      TableRow("Fail Dilepton Triggers (Sata SFs) - 1 Trigger Object", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==1),0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Fail Dilepton Triggers - Failed ID/iso", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==2),0,0,weight_noscale),
      TableRow("Fail Dilepton Triggers (Sata SFs) - Failed ID/iso", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==2),0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Fail Dilepton Triggers - Failed top leg pT cut", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==3),0,0,weight_noscale),
      TableRow("Fail Dilepton Triggers (Sata SFs) - Failed top leg pT cut", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==3),0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Fail Dilepton Triggers - Failed bottom leg pT cut", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==4),0,0,weight_noscale),
      TableRow("Fail Dilepton Triggers (Sata SFs) - Failed bottom leg pT cut", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==4),0,0,weight_noscale*ScaleFactor_Triggers),
      //TableRow("Fail Single and Dimuon Triggers (Data SFs)", 
      //    (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger||l1_mu_trigger_singlemu)&&!hlt_mu_trigger_plus,0,0,weight_noscale),
      //TableRow("Fail Single and Dimuon Triggers (Data SFs) 0mu0ph", 
      //    (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger||l1_mu_trigger_singlemu)&&!hlt_mu_trigger_plus&&(NMuon_trig==0.0&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale),
      //TableRow("Fail Single and Dimuon Triggers (Data SFs) 0mu1ph", 
      //    (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger||l1_mu_trigger_singlemu)&&!hlt_mu_trigger_plus&&(NMuon_trig==0.0&&"(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale),
      //TableRow("Fail Single and Dimuon Triggers (Data SFs) 1mu0ph", 
      //    (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger||l1_mu_trigger_singlemu)&&!hlt_mu_trigger_plus&&(NMuon_trig==1&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale),
      //TableRow("Fail Single and Dimuon Triggers (Data SFs) 1mu1ph", 
      //    (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger||l1_mu_trigger_singlemu)&&!hlt_mu_trigger_plus&&(NMuon_trig==1&&"(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale),
      //TableRow("Fail Single and Dimuon Triggers (Data SFs) 2mu0ph", 
      //    (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger||l1_mu_trigger_singlemu)&&!hlt_mu_trigger_plus&&(NMuon_trig==2&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale),
      //TableRow("Fail Single and Dimuon Triggers (Data SFs) 2mu1ph", 
      //    (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger||l1_mu_trigger_singlemu)&&!hlt_mu_trigger_plus&&(NMuon_trig==2&&"(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale),
      TableRow("Fail Single and Dimuon Triggers (No L1) (Data SFs)", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!hlt_mu_trigger_plus,0,0,weight_noscale),
      TableRow("Fail Single and Dimuon Triggers (Data SFs) 0mu0ph", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!hlt_mu_trigger_plus&&(NMuon_trig==0.0&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale),
      TableRow("Fail Single and Dimuon Triggers (Data SFs) 0mu1ph", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!hlt_mu_trigger_plus&&(NMuon_trig==0.0&&"(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale),
      TableRow("Fail Single and Dimuon Triggers (Data SFs) 1mu0ph", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!hlt_mu_trigger_plus&&(NMuon_trig==1&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale),
      TableRow("Fail Single and Dimuon Triggers (Data SFs) 1mu1ph", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!hlt_mu_trigger_plus&&(NMuon_trig==1&&"(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale),
      TableRow("Fail Single and Dimuon Triggers (Data SFs) 2mu0ph", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!hlt_mu_trigger_plus&&(NMuon_trig==2&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale),
      TableRow("Fail Single and Dimuon Triggers (Data SFs) 2mu1ph", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!hlt_mu_trigger_plus&&(NMuon_trig==2&&"(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale),
      TableRow("Fail mu, mumu, muph Triggers (No L1) (Data SFs)", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!(hlt_mu_trigger_plus||hlt_mu_ph_trigger),0,0,weight_noscale),
      TableRow("Fail mu, mumu, muph Triggers (Data SFs) 0mu0ph", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!(hlt_mu_trigger_plus||hlt_mu_ph_trigger)&&(NMuon_trig==0.0&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale),
      TableRow("Fail mu, mumu, muph Triggers (Data SFs) 0mu1ph", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!(hlt_mu_trigger_plus||hlt_mu_ph_trigger)&&(NMuon_trig==0.0&&"(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale),
      TableRow("Fail mu, mumu, muph Triggers (Data SFs) 1mu0ph", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!(hlt_mu_trigger_plus||hlt_mu_ph_trigger)&&(NMuon_trig==1&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale),
      TableRow("Fail mu, mumu, muph Triggers (Data SFs) 1mu1ph", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!(hlt_mu_trigger_plus||hlt_mu_ph_trigger)&&(NMuon_trig==1&&"(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale),
      TableRow("Fail mu, mumu, muph Triggers (Data SFs) 2mu0ph", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!(hlt_mu_trigger_plus||hlt_mu_ph_trigger)&&(NMuon_trig==2&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale),
      TableRow("Fail mu, mumu, muph Triggers (Data SFs) 2mu1ph", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&!(hlt_mu_trigger_plus||hlt_mu_ph_trigger)&&(NMuon_trig==2&&"(HLT_Photon25||HLT_Photon20_HoverELoose)"),0,0,weight_noscale),
    },signal_procs_noscale,false,true,false,false,false,true).LuminosityTag(total_luminosity_string).Precision(2);

    //pm.Push<Table>("zgtrig_el_"+options.year_string, vector<TableRow>{
    //  //TableRow("L1 Triggers", 
    //  //    (z_decay_pdgid==11)&&l1_el_trigger,0,0,weight),
    //  //TableRow("HLT Triggers", 
    //  //    (z_decay_pdgid==11)&&l1_el_trigger&&hlt_el_trigger,0,0,weight),

    //  TableRow("$Z\\rightarrow e^{+}e^{-}$ decays", 
    //      (z_decay_pdgid==11),0,0,weight),
    //  TableRow("Offline electrons in acceptance", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2),0,0,weight),
    //  TableRow("Offline photon in acceptance", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1),0,0,weight),
    //  TableRow("L1 Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger,0,0,weight),
    //  TableRow("HLT Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&hlt_el_trigger,0,0,weight),
    //  TableRow("HLT Triggers (incl el)", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&hlt_el_trigger_plus,0,0,weight),
    //  TableRow("Events that fail HLT - $0\\gamma 0e$", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&!hlt_el_trigger_plus&&(nel_id_hlt<1)&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)",0,0,weight),
    //  TableRow("Events that fail HLT - $0\\gamma 1e$", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&!hlt_el_trigger_plus&&(nel_id_hlt==1)&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)",0,0,weight),
    //  TableRow("Events that fail HLT - $0\\gamma \\geq 1e$", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&!hlt_el_trigger_plus&&(nel_id_hlt>=2)&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)",0,0,weight),
    //  TableRow("Events that fail HLT - $1\\gamma 0e$", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&!hlt_el_trigger_plus&&(nel_id_hlt<1)&&"(HLT_Photon25||HLT_Photon20_HoverELoose)",0,0,weight),
    //  TableRow("Events that fail HLT - $1\\gamma 1e$", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&!hlt_el_trigger_plus&&(nel_id_hlt==1)&&"(HLT_Photon25||HLT_Photon20_HoverELoose)",0,0,weight),
    //  TableRow("Events that fail HLT - $1\\gamma \\geq 1e$", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&!hlt_el_trigger_plus&&(nel_id_hlt>=2)&&"(HLT_Photon25||HLT_Photon20_HoverELoose)",0,0,weight),
    //  TableRow("HLT Triggers (incl Ele27Ph50)", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&(hlt_el_trigger_plus||HLT_Ele27_CaloIdL_TrackIdL_IsoVL_Photon50),0,0,weight),
    //  TableRow("HLT Triggers (incl Ele27Ph33)", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&(hlt_el_trigger_plus||HLT_Ele27_CaloIdL_TrackIdL_IsoVL_Photon33),0,0,weight),
    //  TableRow("HLT Triggers (incl Ele27Ph25)", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&(hlt_el_trigger_plus||HLT_Ele27_CaloIdL_TrackIdL_IsoVL_Photon25),0,0,weight),
    //  TableRow("HLT Triggers (incl Ele27Ph20)", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&(hlt_el_trigger_plus||HLT_Ele27_CaloIdL_TrackIdL_IsoVL_Photon20_HoverELoose),0,0,weight),
    //  TableRow("HLT Triggers (incl Ele22Ph50)", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&(hlt_el_trigger_plus||HLT_Ele22_CaloIdL_TrackIdL_IsoVL_Photon50),0,0,weight),
    //  TableRow("HLT Triggers (incl Ele22Ph33)", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&(hlt_el_trigger_plus||HLT_Ele22_CaloIdL_TrackIdL_IsoVL_Photon33),0,0,weight),
    //  TableRow("HLT Triggers (incl Ele22Ph25)", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&(hlt_el_trigger_plus||HLT_Ele22_CaloIdL_TrackIdL_IsoVL_Photon25),0,0,weight),
    //  TableRow("HLT Triggers (incl Ele22Ph20)", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&(hlt_el_trigger_plus||HLT_Ele22_CaloIdL_TrackIdL_IsoVL_Photon20_HoverELoose),0,0,weight),
    //  TableRow("HLT Triggers (incl Ele17Ele8Ph25)", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&(hlt_el_trigger_plus||HLT_Ele17_Ele8_CaloIdL_TrackIdL_IsoVL_Photon25),0,0,weight),
    //  TableRow("HLT Triggers (incl Ele17Ele8Ph20)", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&(hlt_el_trigger_plus||HLT_Ele17_Ele8_CaloIdL_TrackIdL_IsoVL_Photon20_HoverELoose),0,0,weight),
    //  //TableRow("HLT Triggers (incl el, el+ph)", 
    //  //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&hlt_el_trigger_plusplus,0,0,weight),
    //  //TableRow("Electrons out of HLT p_{T} acceptance", 
    //  //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&el_oohltpt,0,0,weight),
    //  //TableRow("Electrons nearly out of HLT p_{T} acceptance", 
    //  //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&el_noohltpt,0,0,weight),
    //  //TableRow("Expected events with electrons failing HLT", 
    //  //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&el_oohltpt,0,0,weight*hlt_fail_weights),
    //  //TableRow("HLT fail reason: unknown", 
    //  //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==0.0),0,0,weight),
    //  //TableRow("HLT fail reason: reconstruction", 
    //  //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==1),0,0,weight),
    //  //TableRow("HLT fail reason: reconstruction, e in gap", 
    //  //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==1)&&el_gap,0,0,weight),
    //  //TableRow("HLT fail reason: id/iso", 
    //  //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==2),0,0,weight),
    //  //TableRow("HLT fail reason: top leg", 
    //  //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==3),0,0,weight),
    //  //TableRow("HLT fail reason: bottom leg", 
    //  //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==4),0,0,weight),

    //  TableRow("$Z\\rightarrow \\mu^{+}\\mu^{-}$ decays", 
    //      (z_decay_pdgid==13),0,0,weight),
    //  TableRow("Offline muons in acceptance", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2),0,0,weight),
    //  TableRow("Offline photon in acceptance", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1),0,0,weight),
    //  TableRow("L1 Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger,0,0,weight),
    //  TableRow("HLT Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&hlt_mu_trigger,0,0,weight),
    //  TableRow("HLT Triggers (incl mu)", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&hlt_mu_trigger_plus,0,0,weight),
    //  TableRow("Events that fail HLT - $0\\gamma 0 \\mu$", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&!hlt_mu_trigger_plus&&nmu_id_hlt<1&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)",0,0,weight),
    //  TableRow("Events that fail HLT - $0\\gamma 1 \\mu$", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&!hlt_mu_trigger_plus&&nmu_id_hlt==1&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)",0,0,weight),
    //  TableRow("Events that fail HLT - $0\\gamma 2 \\mu$", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&!hlt_mu_trigger_plus&&nmu_id_hlt>=2&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)",0,0,weight),
    //  TableRow("Events that fail HLT - $1\\gamma 0 \\mu$", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&!hlt_mu_trigger_plus&&nmu_id_hlt<1&&"(HLT_Photon25||HLT_Photon20_HoverELoose)",0,0,weight),
    //  TableRow("Events that fail HLT - $1\\gamma 1 \\mu$", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&!hlt_mu_trigger_plus&&nmu_id_hlt==1&&"(HLT_Photon25||HLT_Photon20_HoverELoose)",0,0,weight),
    //  TableRow("Events that fail HLT - $1\\gamma 2 \\mu$", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&!hlt_mu_trigger_plus&&nmu_id_hlt>=2&&"(HLT_Photon25||HLT_Photon20_HoverELoose)",0,0,weight),
    //  TableRow("HLT Triggers (incl mu, mu+ph)", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&hlt_mu_trigger_plusplus,0,0,weight),
    //  //TableRow("Muons out of HLT acceptance", 
    //  //    (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&hlt_mu_trigger&&mu_oohltpt,0,0,weight),
    //  //TableRow("HLT fail reason: unknown", 
    //  //    (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==0.0),0,0,weight),
    //  //TableRow("HLT fail reason: reconstruction", 
    //  //    (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==1),0,0,weight),
    //  //TableRow("HLT fail reason: iso", 
    //  //    (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==2),0,0,weight),
    //  //TableRow("HLT fail reason: top leg", 
    //  //    (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==3),0,0,weight),
    //  //TableRow("HLT fail reason: bottom leg", 
    //  //    (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==4),0,0,weight),

    //  TableRow("$Z\\rightarrow e^{+}e^{-}$ decays", 
    //      (z_decay_pdgid==11),0,0,weight),
    //  TableRow("L1 Triggers", 
    //      (z_decay_pdgid==11)&&l1_el_trigger,0,0,weight),
    //  TableRow("HLT Triggers", 
    //      (z_decay_pdgid==11)&&l1_el_trigger&&hlt_el_trigger,0,0,weight),
    //  TableRow("HLT fail reason: unknown", 
    //      (z_decay_pdgid==11)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==0.0),0,0,weight),
    //  TableRow("HLT fail reason: reconstruction", 
    //      (z_decay_pdgid==11)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==1),0,0,weight),
    //  TableRow("HLT fail reason: reconstruction, e in gap", 
    //      (z_decay_pdgid==11)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==1)&&el_gap,0,0,weight),
    //  TableRow("HLT fail reason: id/iso", 
    //      (z_decay_pdgid==11)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==2),0,0,weight),
    //  TableRow("HLT fail reason: top leg", 
    //      (z_decay_pdgid==11)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==3),0,0,weight),
    //  TableRow("HLT fail reason: bottom leg", 
    //      (z_decay_pdgid==11)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==4),0,0,weight),

    //  TableRow("$Z\\rightarrow \\mu^{+}\\mu^{-}$ decays", 
    //      (z_decay_pdgid==13),0,0,weight),
    //  TableRow("L1 Triggers", 
    //      (z_decay_pdgid==13)&&l1_mu_trigger,0,0,weight),
    //  TableRow("HLT Triggers", 
    //      (z_decay_pdgid==13)&&l1_mu_trigger&&hlt_mu_trigger,0,0,weight),
    //  TableRow("HLT fail reason: unknwon", 
    //      (z_decay_pdgid==13)&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==0.0),0,0,weight),
    //  TableRow("HLT fail reason: reconstruction", 
    //      (z_decay_pdgid==13)&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==1),0,0,weight),
    //  TableRow("HLT fail reason: iso", 
    //      (z_decay_pdgid==13)&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==2),0,0,weight),
    //  TableRow("HLT fail reason: top leg", 
    //      (z_decay_pdgid==13)&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==3),0,0,weight),
    //  TableRow("HLT fail reason: bottom leg", 
    //      (z_decay_pdgid==13)&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==4),0,0,weight),

    //  //TableRow("$Z\\rightarrow e^{+}e^{-}$ decays, not in l1 acceptance", 
    //  //    (z_decay_pdgid==11)&&((nel_ooeta>=2)||(nel_ool1pt)),0,0,weight),
    //  //TableRow("$Z\\rightarrow e^{+}e^{-}$ decays, not eta acceptance", 
    //  //    (z_decay_pdgid==11)&&((nel_ooeta>=2)),0,0,weight),
    //  //TableRow("$Z\\rightarrow e^{+}e^{-}$ decays, not l1 pt acceptance", 
    //  //    (z_decay_pdgid==11)&&((nel_ool1pt)),0,0,weight),
    //  //TableRow("$Z\\rightarrow e^{+}e^{-}$ decays, offline+L1 but not in hlt pt acceptance", 
    //  //    (z_decay_pdgid==11)&&offline_selection&&l1_el_trigger&&nel_oohltpt,0,0,weight),
    //  //TableRow("HLT Triggers", 
    //  //    (z_decay_pdgid==11)&&l1_el_trigger&&hlt_el_trigger,0,0,weight),
    //},procs,false,true,false,false,false,true).LuminosityTag(total_luminosity_string).Precision(2);
  }
  
  if (make_compare_preselection_table) {
    pm.Push<Table>("ul_trig_eff"+options.year_string, vector<TableRow>{
      TableRow("\\hline\n $Z\\rightarrow e^{+}e^{-}$ decays", 
          (z_decay_pdgid==11),0,0,weight_noscale),
      TableRow("Offline electrons in acceptance", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2),0,0,weight_noscale),
      TableRow("Offline photon in acceptance", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1),0,0,weight_noscale),
      TableRow("L1 Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger,0,0,weight_noscale),
      TableRow("HLT Single or Dilepton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&hlt_el_trigger_plus,0,0,weight_noscale),
      TableRow("HLT Dilepton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&hlt_el_trigger,0,0,weight_noscale),
      TableRow("\\hline\n Baseline Selection", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection,0,0,weight_noscale),
      TableRow("L1 Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection&&l1_el_trigger,0,0,weight_noscale),
      TableRow("HLT Single or Dilepton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection&&l1_el_trigger&&hlt_el_trigger_plus,0,0,weight_noscale),
      TableRow("HLT Dilepton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection&&l1_el_trigger&&hlt_el_trigger,0,0,weight_noscale),
      TableRow("\\hline\n Lepton $p_\\text{T}$ cuts", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&(Lead_SignalElectron_pt>23)&&(Sublead_SignalElectron_pt>12),0,0,weight_noscale),
      TableRow("L1 Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&(Lead_SignalElectron_pt>23)&&(Sublead_SignalElectron_pt>12)&&l1_el_trigger,0,0,weight_noscale),
      TableRow("HLT Single or Dilepton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&(Lead_SignalElectron_pt>23)&&(Sublead_SignalElectron_pt>12)&&l1_el_trigger&&hlt_el_trigger_plus,0,0,weight_noscale),
      TableRow("HLT Dilepton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&(Lead_SignalElectron_pt>23)&&(Sublead_SignalElectron_pt>12)&&l1_el_trigger&&hlt_el_trigger,0,0,weight_noscale),
      TableRow("\\hline\n No Selection", 
          (z_decay_pdgid==11),0,0,weight_noscale),
      TableRow("L1 Triggers", 
          (z_decay_pdgid==11)&&l1_el_trigger,0,0,weight_noscale),
      TableRow("HLT Single or Dilepton Triggers", 
          (z_decay_pdgid==11)&&l1_el_trigger&&hlt_el_trigger_plus,0,0,weight_noscale),
      TableRow("HLT Dilepton Triggers", 
          (z_decay_pdgid==11)&&l1_el_trigger&&hlt_el_trigger,0,0,weight_noscale),

      TableRow("\\hline\n $Z\\rightarrow \\mu^{+}\\mu^{-}$ decays", 
          (z_decay_pdgid==13),0,0,weight_noscale),
      TableRow("Offline muons in acceptance", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2),0,0,weight_noscale),
      TableRow("Offline photon in acceptance", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1),0,0,weight_noscale),
      TableRow("L1 Triggers", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger,0,0,weight_noscale),
      TableRow("HLT Single or Dilepton Triggers", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&hlt_mu_trigger_plus,0,0,weight_noscale),
      TableRow("HLT Dilepton Triggers", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&hlt_mu_trigger,0,0,weight_noscale),
      TableRow("\\hline\n Baseline Selection", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection,0,0,weight_noscale),
      TableRow("L1 Triggers", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection&&l1_mu_trigger,0,0,weight_noscale),
      TableRow("HLT Single or Dilepton Triggers", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection&&l1_mu_trigger&&hlt_mu_trigger_plus,0,0,weight_noscale),
      TableRow("HLT Dilepton Triggers", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection&&l1_mu_trigger&&hlt_mu_trigger,0,0,weight_noscale),
      TableRow("\\hline\n Lepton $p_\\text{T}$ cuts", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&(Lead_SignalMuon_pt>17)&&(Sublead_SignalMuon_pt>10),0,0,weight_noscale),
      TableRow("L1 Triggers", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&(Lead_SignalMuon_pt>17)&&(Sublead_SignalMuon_pt>10)&&l1_mu_trigger,0,0,weight_noscale),
      TableRow("HLT Single or Dilepton Triggers", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&(Lead_SignalMuon_pt>17)&&(Sublead_SignalMuon_pt>10)&&l1_mu_trigger&&hlt_mu_trigger_plus,0,0,weight_noscale),
      TableRow("HLT Dilepton Triggers", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&(Lead_SignalMuon_pt>17)&&(Sublead_SignalMuon_pt>10)&&l1_mu_trigger&&hlt_mu_trigger,0,0,weight_noscale),
      TableRow("\\hline\n No Selection", 
          (z_decay_pdgid==13),0,0,weight_noscale),
      TableRow("L1 Triggers", 
          (z_decay_pdgid==13)&&l1_mu_trigger,0,0,weight_noscale),
      TableRow("HLT Single or Dilepton Triggers", 
          (z_decay_pdgid==13)&&l1_mu_trigger&&hlt_mu_trigger_plus,0,0,weight_noscale),
      TableRow("HLT Dilepton Triggers", 
          (z_decay_pdgid==13)&&l1_mu_trigger&&hlt_mu_trigger,0,0,weight_noscale),
    },signal_procs_noscale,false,true,false,false,false,true).LuminosityTag(total_luminosity_string).Precision(2);
  }
  
  std::vector<std::string> histnames_elhltpt12;
  std::vector<std::string> histnames_elhltpt23;
  std::vector<std::string> histnames_dieleleg12;
  std::vector<std::string> histnames_muhltpt8;
  std::vector<std::string> histnames_muhltpt17;
  std::vector<std::string> histnames_dimuleg8;
  std::vector<std::string> histnames_isoel3235;
  std::vector<std::string> histnames_isomu2427;
  std::vector<std::string> histnames_mu17;
  if (plot_trigeffs) {
    for (unsigned ieta = 0; ieta < el_abseta_bins.size()-1; ieta++) {
      histnames_elhltpt12.push_back("zghlt_dataeff_elhltpt12_eta"+std::to_string(ieta));
      histnames_elhltpt23.push_back("zghlt_dataeff_elhltpt23_eta"+std::to_string(ieta));
      histnames_dieleleg12.push_back("zghlt_dataeff_dieleleg12_eta"+std::to_string(ieta));
      histnames_isoel3235.push_back("zghlt_dataeff_isoel3235_eta"+std::to_string(ieta));
      //pm.Push<EfficiencyPlot>(Axis(el_pt_bins, "Electron_pt", "Offline Electron p_{T} [GeV]", {}),
      //    "HLT_Ele35_WPTight_Gsf"&&(nSignalElectron==2)&&(ZCandidate_mass>81&&ZCandidate_mass<101)
      //    &&Electron_OtherPassReference&&Electron_sig&&Electron_absEta>el_abseta_bins[ieta]
      //    &&Electron_absEta<el_abseta_bins[ieta+1]&&Electron_hltId>0.5,
      //    Electron_hltPt>12,
      //    procs_onelep_data,true,plt_lin_logx).Weight("1").Tag("FixName:zghlt_dataeff_elhltpt12_eta"+std::to_string(ieta)).YTitle("HLT p_{T} > 12 GeV").LuminosityTag(total_luminosity_string);
      pm.Push<EfficiencyPlot>(Axis(el_pt_bins_ptcut, "Electron_pt", "Offline Electron p_{T} [GeV]", {}),
          "HLT_Ele35_WPTight_Gsf"&&(nSignalElectron==2)&&(ZCandidate_mass>81&&ZCandidate_mass<101)
          &&Electron_OtherPassReference&&Electron_sig&&Electron_absEta>el_abseta_bins[ieta]
          &&Electron_absEta<el_abseta_bins[ieta+1]&&Electron_hltId>-0.5,
          Electron_hltPt>23,
          procs_onelep_data,true,plt_lin_logx).Weight("1").Tag("FixName:zghlt_dataeff_elhltpt23_eta"+std::to_string(ieta)).YTitle("HLT p_{T} > 23 GeV").LuminosityTag(total_luminosity_string);
      pm.Push<EfficiencyPlot>(Axis(el_pt_bins, "Electron_pt", "Offline Electron p_{T} [GeV]", {}),
          "HLT_Ele35_WPTight_Gsf"&&(nSignalElectron==2)&&(ZCandidate_mass>81&&ZCandidate_mass<101)
          &&Electron_OtherPassReference&&Electron_OtherPassUpperLeg&&Electron_sig&&!Electron_isLeading
          &&Electron_absEta>el_abseta_bins[ieta]&&Electron_absEta<el_abseta_bins[ieta+1],
          hlt_el_trigger,
          procs_onelep_data,true,plt_lin_logx).Weight("1").Tag("FixName:zghlt_dataeff_dieleleg12_eta"+std::to_string(ieta)).YTitle("Ele12_CaloIdL_TrackIdL_IsoVL").LuminosityTag(total_luminosity_string);
      pm.Push<EfficiencyPlot>(Axis(el_pt_bins, "Electron_pt", "Offline Electron p_{T} [GeV]", {}),
          "HLT_PFMET120_PFMHT120_IDTight&&MET_pt>150"&&(nSignalElectron==1)&&(Electron_sig>0.5)&&"Electron_mvaFall17V2Iso_WP90"&&Electron_absEta>el_abseta_bins[ieta]
          &&Electron_absEta<el_abseta_bins[ieta+1],
          "HLT_Ele32_WPTight_Gsf_L1DoubleEG||HLT_Ele35_WPTight_Gsf",
          procs_met_data,true,plt_lin_logx).Weight("1").Tag("FixName:zghlt_dataeff_isoel3235_eta"+std::to_string(ieta)).YTitle("Single Electron Triggers").LuminosityTag(total_luminosity_string);
    }
    for (unsigned ieta = 0; ieta < mu_abseta_bins.size()-1; ieta++) {
      histnames_muhltpt8.push_back("zghlt_dataeff_muhltpt8_eta"+std::to_string(ieta));
      histnames_muhltpt17.push_back("zghlt_dataeff_muhltpt17_eta"+std::to_string(ieta));
      histnames_dimuleg8.push_back("zghlt_dataeff_dimuleg8_eta"+std::to_string(ieta));
      histnames_isomu2427.push_back("zghlt_dataeff_isomu2427_eta"+std::to_string(ieta));
      histnames_mu17.push_back("zghlt_dataeff_mu17_eta"+std::to_string(ieta));
      pm.Push<EfficiencyPlot>(Axis(mu_pt_bins_ptcut, "Muon_pt", "Offline Muon p_{T} [GeV]", {}),
          "HLT_IsoMu27"&&(nSignalMuon==2)&&(ZCandidate_mass>81&&ZCandidate_mass<101)
          &&Muon_OtherPassReference&&Muon_sig&&Muon_absEta>mu_abseta_bins[ieta]
          &&Muon_absEta<mu_abseta_bins[ieta+1]&&Muon_hltId>-0.5,
          Muon_hltPt>17,
          procs_onelep_data,true,plt_lin_logx).Weight("1").Tag("FixName:zghlt_dataeff_muhltpt17_eta"+std::to_string(ieta)).YTitle("HLT p_{T} > 17 GeV").LuminosityTag(total_luminosity_string);
      pm.Push<EfficiencyPlot>(Axis(mu_pt_bins, "Muon_pt", "Offline Muon p_{T} [GeV]", {}),
          "HLT_IsoMu27"&&(nSignalMuon==2)&&(ZCandidate_mass>81&&ZCandidate_mass<101)
          &&Muon_OtherPassReference&&Muon_OtherPassUpperLeg&&Muon_sig&&!Muon_isLeading
          &&Muon_absEta>mu_abseta_bins[ieta]&&Muon_absEta<mu_abseta_bins[ieta+1],
          hlt_mu_trigger,
          procs_onelep_data,true,plt_lin_logx).Weight("1").Tag("FixName:zghlt_dataeff_dimuleg8_eta"+std::to_string(ieta)).YTitle("Mu8_TrkIsoVVVL").LuminosityTag(total_luminosity_string);
      pm.Push<EfficiencyPlot>(Axis(mu_pt_bins, "Muon_pt", "Offline Muon p_{T} [GeV]", {}),
          "HLT_PFMET120_PFMHT120_IDTight&&MET_pt>150"&&(nSignalMuon==1)&&Muon_sig&&"Muon_mediumId"&&Muon_absEta>mu_abseta_bins[ieta]
          &&Muon_absEta<mu_abseta_bins[ieta+1],
          "HLT_IsoMu27||HLT_IsoMu24",
          procs_met_data,true,plt_lin_logx).Weight("1").Tag("FixName:zghlt_dataeff_isomu2427_eta"+std::to_string(ieta)).YTitle("Single Muon Triggers").LuminosityTag(total_luminosity_string);
      //pm.Push<EfficiencyPlot>(Axis(mu_pt_bins, "Muon_pt", "Offline Muon p_{T} [GeV]", {}),
      //    "HLT_IsoMu27"&&(nSignalMuon==2)&&(ZCandidate_mass>81&&ZCandidate_mass<101)
      //    &&Muon_OtherPassReference&&Muon_sig&&Muon_absEta>mu_abseta_bins[ieta]
      //    &&Muon_absEta<mu_abseta_bins[ieta+1]&&Muon_hltId>0.5,
      //    Muon_hltPt>8,
      //    procs_onelep_data,true,plt_lin_logx).Weight("1").Tag("FixName:zghlt_dataeff_muhltpt8_eta"+std::to_string(ieta)).YTitle("HLT p_{T} > 8 GeV").LuminosityTag(total_luminosity_string);
      //pm.Push<EfficiencyPlot>(Axis(mu_pt_bins, "Muon_pt", "Offline Muon p_{T} [GeV]", {}),
      //    "HLT_PFMET120_PFMHT120_IDTight&&MET_pt>150"&&(nSignalMuon==1)&&Muon_sig&&"Muon_mediumId"&&Muon_absEta>mu_abseta_bins[ieta]
      //    &&Muon_absEta<mu_abseta_bins[ieta+1],
      //    "HLT_Mu17",
      //    procs_met_data,true,plt_lin_logx).Weight("1").Tag("FixName:zghlt_dataeff_mu17_eta"+std::to_string(ieta)).YTitle("Single Muon Triggers").LuminosityTag(total_luminosity_string);
    }
  }

  std::vector<std::string> histnames_sigelhltpt23;
  std::vector<std::string> histnames_sigdieleleg12;
  std::vector<std::string> histnames_sigdieleleg23;
  std::vector<std::string> histnames_sigmuhltpt17;
  std::vector<std::string> histnames_sigdimuleg8;
  std::vector<std::string> histnames_sigdimuleg17;
  std::vector<std::string> histnames_sigisoel3235;
  std::vector<std::string> histnames_sigisomu2427;
  if (plot_trigeffs_signal) {
    for (unsigned ieta = 0; ieta < el_abseta_bins.size()-1; ieta++) {
      histnames_sigelhltpt23.push_back("zghlt_sigeff_elhltpt23_eta"+std::to_string(ieta));
      histnames_sigdieleleg12.push_back("zghlt_sigeff_dieleleg12_eta"+std::to_string(ieta));
      histnames_sigdieleleg23.push_back("zghlt_sigeff_dieleleg23_eta"+std::to_string(ieta));
      histnames_sigisoel3235.push_back("zghlt_sigeff_isoel3235_eta"+std::to_string(ieta));
      pm.Push<EfficiencyPlot>(Axis(el_pt_bins_ptcut, "Electron_pt", "Offline Electron p_{T} [GeV]", {}),
          (nSignalElectron==2)&&Electron_sig&&Electron_absEta>el_abseta_bins[ieta]
          &&Electron_absEta<el_abseta_bins[ieta+1]&&Electron_hltId>-0.5,
          Electron_hltPt>23,
          signal_procs_noscale,true,plt_lin_logx).Weight("1").Tag("FixName:zghlt_sigeff_elhltpt23_eta"+std::to_string(ieta)).YTitle("HLT p_{T} > 23 GeV").LuminosityTag(total_luminosity_string);
      //for statistics at high pT allow leading leptons where pT is not a concern
      pm.Push<EfficiencyPlot>(Axis(el_pt_bins, "Electron_pt", "Offline Electron p_{T} [GeV]", {}),
          (nSignalElectron==2)&&Electron_sig&&delta_r_cut
          &&((Electron_OtherPassUpperLeg&&!Electron_isLeading)||("Electron_pt>30"&&Electron_OtherPassLowerLeg&&Electron_isLeading))
          &&Electron_absEta>el_abseta_bins[ieta]&&Electron_absEta<el_abseta_bins[ieta+1],
          hlt_el_trigger,
          signal_procs_noscale,true,plt_lin_logx).Weight("1").Tag("FixName:zghlt_sigeff_dieleleg12_eta"+std::to_string(ieta)).YTitle("Ele12_CaloIdL_TrackIdL_IsoVL").LuminosityTag(total_luminosity_string);
      pm.Push<EfficiencyPlot>(Axis(el_pt_bins, "Electron_pt", "Offline Electron p_{T} [GeV]", {}),
          (nSignalElectron==2)&&Electron_OtherPassLowerLeg&&Electron_sig&&Electron_isLeading&&delta_r_cut
          &&Electron_absEta>el_abseta_bins[ieta]&&Electron_absEta<el_abseta_bins[ieta+1],
          hlt_el_trigger,
          signal_procs_noscale,true,plt_lin_logx).Weight("1").Tag("FixName:zghlt_sigeff_dieleleg23_eta"+std::to_string(ieta)).YTitle("Ele23_CaloIdL_TrackIdL_IsoVL").LuminosityTag(total_luminosity_string);
      pm.Push<EfficiencyPlot>(Axis(el_pt_bins, "Electron_pt", "Offline Electron p_{T} [GeV]", {}),
          OneElectronOutOfAcceptance&&(nSignalElectron==1)&&Electron_sig&&delta_r_cut
          &&Electron_absEta>el_abseta_bins[ieta]&&Electron_absEta<el_abseta_bins[ieta+1],
          "HLT_Ele32_WPTight_Gsf_L1DoubleEG||HLT_Ele35_WPTight_Gsf",
          signal_procs_noscale,true,plt_lin_logx).Weight("1").Tag("FixName:zghlt_sigeff_isoel3235_eta"+std::to_string(ieta)).YTitle("Single Electron Triggers").LuminosityTag(total_luminosity_string);
    }
    for (unsigned ieta = 0; ieta < mu_abseta_bins.size()-1; ieta++) {
      histnames_sigmuhltpt17.push_back("zghlt_sigeff_muhltpt17_eta"+std::to_string(ieta));
      histnames_sigdimuleg8.push_back("zghlt_sigeff_dimuleg8_eta"+std::to_string(ieta));
      histnames_sigdimuleg17.push_back("zghlt_sigeff_dimuleg17_eta"+std::to_string(ieta));
      histnames_sigisomu2427.push_back("zghlt_sigeff_isomu2427_eta"+std::to_string(ieta));
      pm.Push<EfficiencyPlot>(Axis(mu_pt_bins_ptcut, "Muon_pt", "Offline Muon p_{T} [GeV]", {}),
          (nSignalMuon==2)&&Muon_sig&&delta_r_cut&&Muon_absEta>mu_abseta_bins[ieta]
          &&Muon_absEta<mu_abseta_bins[ieta+1]&&Muon_hltId>-0.5,
          Muon_hltPt>17,
          signal_procs_noscale,true,plt_lin_logx).Weight("1").Tag("FixName:zghlt_sigeff_muhltpt17_eta"+std::to_string(ieta)).YTitle("HLT p_{T} > 17 GeV").LuminosityTag(total_luminosity_string);
      //for statistics at high pT allow leading leptons where pT is not a concern
      pm.Push<EfficiencyPlot>(Axis(mu_pt_bins, "Muon_pt", "Offline Muon p_{T} [GeV]", {}),
          (nSignalMuon==2)&&Muon_sig&&delta_r_cut
          &&((Muon_OtherPassUpperLeg&&!Muon_isLeading)||("Muon_pt>30"&&Muon_OtherPassLowerLeg&&Muon_isLeading))
          &&Muon_absEta>mu_abseta_bins[ieta]&&Muon_absEta<mu_abseta_bins[ieta+1],
          hlt_mu_trigger,
          signal_procs_noscale,true,plt_lin_logx).Weight("1").Tag("FixName:zghlt_sigeff_dimuleg8_eta"+std::to_string(ieta)).YTitle("Mu8_TrkIsoVVVL").LuminosityTag(total_luminosity_string);
      pm.Push<EfficiencyPlot>(Axis(mu_pt_bins, "Muon_pt", "Offline Muon p_{T} [GeV]", {}),
          (nSignalMuon==2)&&Muon_OtherPassLowerLeg&&Muon_sig&&Muon_isLeading&&delta_r_cut
          &&Muon_absEta>mu_abseta_bins[ieta]&&Muon_absEta<mu_abseta_bins[ieta+1],
          hlt_mu_trigger,
          signal_procs_noscale,true,plt_lin_logx).Weight("1").Tag("FixName:zghlt_sigeff_dimuleg17_eta"+std::to_string(ieta)).YTitle("Mu17_TrkIsoVVVL").LuminosityTag(total_luminosity_string);
      pm.Push<EfficiencyPlot>(Axis(mu_pt_bins, "Muon_pt", "Offline Muon p_{T} [GeV]", {}),
          OneMuonOutOfAcceptance&&(nSignalMuon==1)&&Muon_sig&&delta_r_cut&&Muon_absEta>mu_abseta_bins[ieta]
          &&Muon_absEta<mu_abseta_bins[ieta+1],
          "HLT_IsoMu27||HLT_IsoMu24",
          signal_procs_noscale,true,plt_lin_logx).Weight("1").Tag("FixName:zghlt_sigeff_isomu2427_eta"+std::to_string(ieta)).YTitle("Single Muon Triggers").LuminosityTag(total_luminosity_string);
    }
  }

  //bool plot_diph_trigeffs = true;
  //if (plot_diph_trigeffs) {
  //  for (unsigned ieta = 0; ieta < el_abseta_bins.size()-1; ieta++) {
  //    histnames_elhltpt12.push_back("zghlt_dataeff_elhltpt12_eta"+std::to_string(ieta));
  //    histnames_elhltpt23.push_back("zghlt_dataeff_elhltpt23_eta"+std::to_string(ieta));
  //    histnames_dieleleg12.push_back("zghlt_dataeff_dieleleg12_eta"+std::to_string(ieta));
  //    histnames_isoel3235.push_back("zghlt_dataeff_isoel3235_eta"+std::to_string(ieta));
  //    //scEt = pt*(scetoverpt+1)
  //    pm.Push<EfficiencyPlot>(Axis(el_pt_bins_ptcut, Electron_scEt, "Offline Electron SC E_{T} [GeV]", {}),
  //        "HLT_Ele35_WPTight_Gsf"&&(nSignalElectron==2)&&(ZCandidate_mass>81&&ZCandidate_mass<101)
  //        &&Electron_OtherPassReference&&Electron_sig&&Electron_absEta>el_abseta_bins[ieta]
  //        &&Electron_absEta<el_abseta_bins[ieta+1]&&Electron_hltId>-0.5,
  //        hlt_doubleph_trigger,
  //        procs_onelep_data,true,plt_lin_logx).Weight("1").Tag("FixName:zghlt_dataeff_elhltpt23_eta"+std::to_string(ieta)).YTitle("HLT p_{T} > 23 GeV").LuminosityTag(total_luminosity_string);
  //    pm.Push<EfficiencyPlot>(Axis(el_pt_bins, "Electron_pt", "Offline Electron p_{T} [GeV]", {}),
  //        "HLT_Ele35_WPTight_Gsf"&&(nSignalElectron==2)&&(ZCandidate_mass>81&&ZCandidate_mass<101)
  //        &&Electron_OtherPassReference&&Electron_OtherPassUpperLeg&&Electron_sig&&!Electron_isLeading
  //        &&Electron_absEta>el_abseta_bins[ieta]&&Electron_absEta<el_abseta_bins[ieta+1],
  //        hlt_el_trigger,
  //        procs_onelep_data,true,plt_lin_logx).Weight("1").Tag("FixName:zghlt_dataeff_dieleleg12_eta"+std::to_string(ieta)).YTitle("Ele12_CaloIdL_TrackIdL_IsoVL").LuminosityTag(total_luminosity_string);
  //    pm.Push<EfficiencyPlot>(Axis(el_pt_bins, "Electron_pt", "Offline Electron p_{T} [GeV]", {}),
  //        "HLT_PFMET120_PFMHT120_IDTight&&MET_pt>150"&&(nSignalElectron==1)&&(Electron_sig>0.5)&&"Electron_mvaFall17V2Iso_WP90"&&Electron_absEta>el_abseta_bins[ieta]
  //        &&Electron_absEta<el_abseta_bins[ieta+1],
  //        "HLT_Ele32_WPTight_Gsf_L1DoubleEG||HLT_Ele35_WPTight_Gsf",
  //        procs_met_data,true,plt_lin_logx).Weight("1").Tag("FixName:zghlt_dataeff_isoel3235_eta"+std::to_string(ieta)).YTitle("Single Electron Triggers").LuminosityTag(total_luminosity_string);
  //  }
  //}

  //pm.Push<Hist1D>(Axis(30, 0.0, 100.0, "Muon_pt", "Muon p_{T}", {}),
  //    "HLT_IsoMu27",
  //    procs_onelep_data, plt_lin).Weight("1")
  //    .Tag("FixName:zghlt__debugcf0")
  //    .LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(30, 0.0, 100.0, "Muon_pt", "Muon p_{T}", {}),
  //    "HLT_IsoMu27"&&(nSignalMuon==2),
  //    procs_onelep_data, plt_lin).Weight("1")
  //    .Tag("FixName:zghlt__debugcf1")
  //    .LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(30, 0.0, 100.0, "Muon_pt", "Muon p_{T}", {}),
  //    "HLT_IsoMu27"&&(nSignalMuon==2)&&(ZCandidate_mass>81&&ZCandidate_mass<101),
  //    procs_onelep_data, plt_lin).Weight("1")
  //    .Tag("FixName:zghlt__debugcf2")
  //    .LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(30, 0.0, 100.0, "Muon_pt", "Muon p_{T}", {}),
  //    "HLT_IsoMu27"&&(nSignalMuon==2)&&(ZCandidate_mass>81&&ZCandidate_mass<101)&&Muon_OtherPassReference,
  //    procs_onelep_data, plt_lin).Weight("1")
  //    .Tag("FixName:zghlt__debugcf3")
  //    .LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(30, 0.0, 100.0, "Muon_pt", "Muon p_{T}", {}),
  //    "HLT_IsoMu27"&&(nSignalMuon==2)&&(ZCandidate_mass>81&&ZCandidate_mass<101)&&Muon_OtherPassReference&&Muon_hltIndex>=0.0&&Muon_hltPt>8,
  //    procs_onelep_data, plt_lin).Weight("1")
  //    .Tag("FixName:zghlt__debugcf3hipt")
  //    .LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(30, 0.0, 100.0, "Muon_pt", "Muon p_{T}", {}),
  //    "HLT_IsoMu27"&&(nSignalMuon==2)&&(ZCandidate_mass>81&&ZCandidate_mass<101)&&Muon_OtherPassReference&&Muon_hltIndex>=0.0&&Muon_hltPt<8,
  //    procs_onelep_data, plt_lin).Weight("1")
  //    .Tag("FixName:zghlt__debugcf3lopt")
  //    .LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(30, 0.0, 100.0, "Muon_pt", "Muon p_{T}", {}),
  //    "HLT_IsoMu27"&&(nSignalMuon==2)&&(ZCandidate_mass>81&&ZCandidate_mass<101)&&Muon_OtherPassReference&&Muon_sig&&Muon_hltIndex>=0.0&&Muon_hltPt>8,
  //    procs_onelep_data, plt_lin).Weight("1")
  //    .Tag("FixName:zghlt__debugcf4hipt")
  //    .LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(30, 0.0, 100.0, "Muon_pt", "Muon p_{T}", {}),
  //    "HLT_IsoMu27"&&(nSignalMuon==2)&&(ZCandidate_mass>81&&ZCandidate_mass<101)&&Muon_OtherPassReference&&Muon_sig&&Muon_hltIndex>=0.0&&Muon_hltPt<8,
  //    procs_onelep_data, plt_lin).Weight("1")
  //    .Tag("FixName:zghlt__debugcf4lopt")
  //    .LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(30, 0.0, 100.0, "Muon_pt", "Muon p_{T}", {}),
  //    "HLT_IsoMu27"&&(nSignalMuon==2)&&(ZCandidate_mass>81&&ZCandidate_mass<101)&&Muon_OtherPassReference&&Muon_sig&&Muon_OtherPassUpperLeg&&Muon_hltIndex>=0.0&&Muon_hltPt>8,
  //    procs_onelep_data, plt_lin).Weight("1")
  //    .Tag("FixName:zghlt__debugcf5hipt")
  //    .LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(30, 0.0, 100.0, "Muon_pt", "Muon p_{T}", {}),
  //    "HLT_IsoMu27"&&(nSignalMuon==2)&&(ZCandidate_mass>81&&ZCandidate_mass<101)&&Muon_OtherPassReference&&Muon_sig&&Muon_OtherPassUpperLeg&&Muon_hltIndex>=0.0&&Muon_hltPt<8,
  //    procs_onelep_data, plt_lin).Weight("1")
  //    .Tag("FixName:zghlt__debugcf5lopt")
  //    .LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(30, 0.0, 100.0, "Muon_pt", "Muon p_{T}", {}),
  //    "HLT_IsoMu27"&&(nSignalMuon==2)&&(ZCandidate_mass>81&&ZCandidate_mass<101)&&Muon_OtherPassReference&&Muon_sig&&Muon_OtherPassUpperLeg&&Muon_hltIndex>=0.0&&Muon_hltPt>8&&Muon_hltId>0.5,
  //    procs_onelep_data, plt_lin).Weight("1")
  //    .Tag("FixName:zghlt__debugcf6hipt")
  //    .LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(30, 0.0, 100.0, "Muon_pt", "Muon p_{T}", {}),
  //    "HLT_IsoMu27"&&(nSignalMuon==2)&&(ZCandidate_mass>81&&ZCandidate_mass<101)&&Muon_OtherPassReference&&Muon_sig&&Muon_OtherPassUpperLeg&&Muon_hltIndex>=0.0&&Muon_hltPt<8&&Muon_hltId>0.5,
  //    procs_onelep_data, plt_lin).Weight("1")
  //    .Tag("FixName:zghlt__debugcf6lopt")
  //    .LuminosityTag(total_luminosity_string);

  if (make_datamctrigefftable) {
    //pm.Push<Table>("datamc_trig_eff"+options.year_string, vector<TableRow>{
    //  TableRow("\\hline \\hline\n $Z\\rightarrow e^{+}e^{-}$ decays", 
    //      (z_decay_pdgid==11),0,0,weight_noscale),
    //  TableRow("Offline leptons and photon in acceptance", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1),0,0,weight_noscale),
    //  TableRow("MC L1 Trigger", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger,0,0,weight_noscale),
    //  TableRow("\\hline\n MC Single Lepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&hlt_single_el_trigger,0,0,weight_noscale),
    //  TableRow("MC Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&hlt_el_trigger,0,0,weight_noscale),
    //  TableRow("MC Single or Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger&&hlt_el_trigger_plus,0,0,weight_noscale),
    //  TableRow("\\hline\n Data Single Lepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger,0,0,weight_noscale*eff_singlelep),
    //  TableRow("Data Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger,0,0,weight_noscale*eff_dilep),
    //  TableRow("Data Dilepton Triggers (HIG-19-014)", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger,0,0,weight_noscale*eff_dilep_hig19014),
    //  TableRow("Data Single or Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&l1_el_trigger,0,0,weight_noscale*eff_singleanddilep),
    //  TableRow("\\hline \\hline\n Baseline Selection", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection,0,0,weight_noscale),
    //  TableRow("MC L1 Trigger", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection&&l1_el_trigger,0,0,weight_noscale),
    //  TableRow("\\hline\n MC Single Lepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection&&l1_el_trigger&&hlt_single_el_trigger,0,0,weight_noscale),
    //  TableRow("MC Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection&&l1_el_trigger&&hlt_el_trigger,0,0,weight_noscale),
    //  TableRow("MC Single or Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection&&l1_el_trigger&&hlt_el_trigger_plus,0,0,weight_noscale),
    //  TableRow("\\hline\n Data Single Lepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection&&l1_el_trigger,0,0,weight_noscale*eff_singlelep),
    //  TableRow("Data Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection&&l1_el_trigger,0,0,weight_noscale*eff_dilep),
    //  TableRow("Data Dilepton Triggers (HIG-19-014)", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection&&l1_el_trigger,0,0,weight_noscale*eff_dilep_hig19014),
    //  TableRow("Data Single or Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection&&l1_el_trigger,0,0,weight_noscale*eff_singleanddilep),
    //  TableRow("\\hline \\hline\n High Lead p_{T} MC L1 Trigger", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection&&l1_el_trigger&&highleadpt,0,0,weight_noscale),
    //  TableRow("\\hline\n MC Single Lepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection&&highleadpt&&l1_el_trigger&&hlt_single_el_trigger,0,0,weight_noscale),
    //  TableRow("MC Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection&&highleadpt&&l1_el_trigger&&hlt_el_trigger,0,0,weight_noscale),
    //  TableRow("MC Single or Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection&&highleadpt&&l1_el_trigger&&hlt_el_trigger_plus,0,0,weight_noscale),
    //  TableRow("\\hline\n Data Single Lepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection&&highleadpt&&l1_el_trigger,0,0,weight_noscale*eff_singlelep),
    //  TableRow("Data Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection&&highleadpt&&l1_el_trigger,0,0,weight_noscale*eff_dilep),
    //  TableRow("Data Dilepton Triggers (HIG-19-014)", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection&&highleadpt&&l1_el_trigger,0,0,weight_noscale*eff_dilep_hig19014),
    //  TableRow("Data Single or Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection&&highleadpt&&l1_el_trigger,0,0,weight_noscale*eff_singleanddilep),

    //  TableRow("\\hline \\hline\n $Z\\rightarrow \\mu^{+}\\mu^{-}$ decays", 
    //      (z_decay_pdgid==13),0,0,weight_noscale),
    //  TableRow("Offline leptons and photon in acceptance", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1),0,0,weight_noscale),
    //  TableRow("MC L1 Trigger", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger,0,0,weight_noscale),
    //  TableRow("\\hline\n MC Single Lepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&hlt_single_mu_trigger,0,0,weight_noscale),
    //  TableRow("MC Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&hlt_mu_trigger,0,0,weight_noscale),
    //  TableRow("MC Single or Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger&&hlt_mu_trigger_plus,0,0,weight_noscale),
    //  TableRow("\\hline\n Data Single Lepton Triggers ", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger,0,0,weight_noscale*eff_singlelep),
    //  TableRow("Data Dilepton Lepton Triggers ", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger,0,0,weight_noscale*eff_dilep),
    //  TableRow("Data Dilepton Lepton Triggers (HIG-19-014)", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger,0,0,weight_noscale*eff_dilep_hig19014),
    //  TableRow("Data Single or Dilepton Lepton Triggers ", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&l1_mu_trigger,0,0,weight_noscale*eff_singleanddilep),
    //  TableRow("\\hline \\hline\n Baseline Selection", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection,0,0,weight_noscale),
    //  TableRow("MC L1 Trigger", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection&&l1_mu_trigger,0,0,weight_noscale),
    //  TableRow("\\hline\n MC Single Lepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection&&l1_mu_trigger&&hlt_single_mu_trigger,0,0,weight_noscale),
    //  TableRow("MC Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection&&l1_mu_trigger&&hlt_mu_trigger,0,0,weight_noscale),
    //  TableRow("MC Single or Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection&&l1_mu_trigger&&hlt_mu_trigger_plus,0,0,weight_noscale),
    //  TableRow("\\hline\n Data Single Lepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection&&l1_mu_trigger,0,0,weight_noscale*eff_singlelep),
    //  TableRow("Data Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection&&l1_mu_trigger,0,0,weight_noscale*eff_dilep),
    //  TableRow("Data Dilepton Triggers (HIG-19-014)", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection&&l1_mu_trigger,0,0,weight_noscale*eff_dilep_hig19014),
    //  TableRow("Data Single or Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection&&l1_mu_trigger,0,0,weight_noscale*eff_singleanddilep),
    //  TableRow("\\hline \\hline\n High lead p_{T} MC L1 Trigger", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection&&highleadpt&&l1_mu_trigger,0,0,weight_noscale),
    //  TableRow("\\hline\n MC Single Lepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection&&highleadpt&&l1_mu_trigger&&hlt_single_mu_trigger,0,0,weight_noscale),
    //  TableRow("MC Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection&&highleadpt&&l1_mu_trigger&&hlt_mu_trigger,0,0,weight_noscale),
    //  TableRow("MC Single or Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection&&highleadpt&&l1_mu_trigger&&hlt_mu_trigger_plus,0,0,weight_noscale),
    //  TableRow("\\hline\n Data Single Lepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection&&highleadpt&&l1_mu_trigger,0,0,weight_noscale*eff_singlelep),
    //  TableRow("Data Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection&&highleadpt&&l1_mu_trigger,0,0,weight_noscale*eff_dilep),
    //  TableRow("Data Dilepton Triggers (HIG-19-014)", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection&&highleadpt&&l1_mu_trigger,0,0,weight_noscale*eff_dilep_hig19014),
    //  TableRow("Data Single or Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection&&highleadpt&&l1_mu_trigger,0,0,weight_noscale*eff_singleanddilep),
    //},signal_procs_noscale,false,true,false,false,false,true).LuminosityTag(total_luminosity_string).Precision(2);

    //pm.Push<Table>("datamc_trig_eff"+options.year_string, vector<TableRow>{
    //  TableRow("\\hline \\hline\n $Z\\rightarrow e^{+}e^{-}$ decays", 
    //      (z_decay_pdgid==11),0,0,weight_noscale),
    //  TableRow("Offline leptons and photon in acceptance", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1),0,0,weight_noscale),
    //  TableRow("Baseline Selection", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,0,0,weight_noscale),
    //  TableRow("MC L1 Trigger", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger,0,0,weight_noscale),
    //  TableRow("\\hline \n MC Single Lepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&hlt_single_el_trigger,0,0,weight_noscale),
    //  TableRow("MC Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&hlt_el_trigger,0,0,weight_noscale),
    //  TableRow("MC Single or Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&hlt_el_trigger_plus,0,0,weight_noscale),
    //  TableRow("MC Single or Dilepton (incl noiso) Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&(hlt_el_trigger_plus||hlt_double_el_noiso_trigger),0,0,weight_noscale),
    //  TableRow("MC Single (incl noiso) or Dilepton (incl noiso) Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&(hlt_el_trigger_plus||hlt_double_el_noiso_trigger||hlt_single_el_noiso_trigger),0,0,weight_noscale),
    //  TableRow("MC Single (incl noiso, presc.) or Dilepton (incl noiso) Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&(hlt_el_trigger_plus||hlt_double_el_noiso_trigger||hlt_single_el_noiso_trigger||hlt_single_el_prescale_trigger),0,0,weight_noscale),
    //  TableRow("Data Single Lepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger,0,0,weight_noscale*eff_singlelep),
    //  TableRow("Data Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger,0,0,weight_noscale*eff_dilep),
    //  TableRow("Data Dilepton Triggers (HIG-19-014)", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger,0,0,weight_noscale*eff_dilep_hig19014),
    //  TableRow("Data Single or Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger,0,0,weight_noscale*eff_singleanddilep),
    //  TableRow("\\hline \n MC Single Lepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&midleadpt&&l1_el_trigger&&hlt_single_el_trigger,0,0,weight_noscale),
    //  TableRow("MC Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&midleadpt&&l1_el_trigger&&hlt_el_trigger,0,0,weight_noscale),
    //  TableRow("MC Single or Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&midleadpt&&l1_el_trigger&&hlt_el_trigger_plus,0,0,weight_noscale),
    //  TableRow("Data Single Lepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&midleadpt&&l1_el_trigger,0,0,weight_noscale*eff_singlelep),
    //  TableRow("Data Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&midleadpt&&l1_el_trigger,0,0,weight_noscale*eff_dilep),
    //  TableRow("Data Dilepton Triggers (HIG-19-014)", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&midleadpt&&l1_el_trigger,0,0,weight_noscale*eff_dilep_hig19014),
    //  TableRow("Data Single or Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&midleadpt&&l1_el_trigger,0,0,weight_noscale*eff_singleanddilep),
    //  TableRow("\\hline \n MC Single Lepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&highleadpt&&l1_el_trigger&&hlt_single_el_trigger,0,0,weight_noscale),
    //  TableRow("MC Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&highleadpt&&l1_el_trigger&&hlt_el_trigger,0,0,weight_noscale),
    //  TableRow("MC Single or Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&highleadpt&&l1_el_trigger&&hlt_el_trigger_plus,0,0,weight_noscale),
    //  TableRow("Data Single Lepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&highleadpt&&l1_el_trigger,0,0,weight_noscale*eff_singlelep),
    //  TableRow("Data Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&highleadpt&&l1_el_trigger,0,0,weight_noscale*eff_dilep),
    //  TableRow("Data Dilepton Triggers (HIG-19-014)", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&highleadpt&&l1_el_trigger,0,0,weight_noscale*eff_dilep_hig19014),
    //  TableRow("Data Single or Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&highleadpt&&l1_el_trigger,0,0,weight_noscale*eff_singleanddilep),
    //  TableRow("\\hline \n MC Single Lepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&highpt&&l1_el_trigger&&hlt_single_el_trigger,0,0,weight_noscale),
    //  TableRow("MC Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&highpt&&l1_el_trigger&&hlt_el_trigger,0,0,weight_noscale),
    //  TableRow("MC Single or Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&highpt&&l1_el_trigger&&hlt_el_trigger_plus,0,0,weight_noscale),
    //  TableRow("Data Single Lepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&highpt&&l1_el_trigger,0,0,weight_noscale*eff_singlelep),
    //  TableRow("Data Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&highpt&&l1_el_trigger,0,0,weight_noscale*eff_dilep),
    //  TableRow("Data Dilepton Triggers (HIG-19-014)", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&highpt&&l1_el_trigger,0,0,weight_noscale*eff_dilep_hig19014),
    //  TableRow("Data Single or Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&highpt&&l1_el_trigger,0,0,weight_noscale*eff_singleanddilep),
    //  TableRow("\\hline \n MC Single Lepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut&&l1_el_trigger&&hlt_single_el_trigger,0,0,weight_noscale),
    //  TableRow("MC Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut&&l1_el_trigger&&hlt_el_trigger,0,0,weight_noscale),
    //  TableRow("MC Single or Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut&&l1_el_trigger&&hlt_el_trigger_plus,0,0,weight_noscale),
    //  TableRow("Data Single Lepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut&&l1_el_trigger,0,0,weight_noscale*eff_singlelep),
    //  TableRow("Data Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut&&l1_el_trigger,0,0,weight_noscale*eff_dilep),
    //  TableRow("Data Dilepton Triggers (HIG-19-014)", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut&&l1_el_trigger,0,0,weight_noscale*eff_dilep_hig19014),
    //  TableRow("Data Single or Dilepton Triggers", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut&&l1_el_trigger,0,0,weight_noscale*eff_singleanddilep),

    //  TableRow("\\hline \\hline\n $Z\\rightarrow \\mu^{+}\\mu^{-}$ decays", 
    //      (z_decay_pdgid==13),0,0,weight_noscale),
    //  TableRow("Offline leptons and photon in acceptance", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1),0,0,weight_noscale),
    //  TableRow("Baseline Selection", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,0,0,weight_noscale),
    //  TableRow("MC L1 Trigger", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_mu_trigger,0,0,weight_noscale),
    //  TableRow("\\hline\n MC Single Lepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_mu_trigger&&hlt_single_mu_trigger,0,0,weight_noscale),
    //  TableRow("MC Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_mu_trigger&&hlt_mu_trigger,0,0,weight_noscale),
    //  TableRow("MC Single or Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_mu_trigger&&hlt_mu_trigger_plus,0,0,weight_noscale),
    //  TableRow("MC Single or Dilepton or Mu+Ph Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_mu_trigger&&(hlt_mu_trigger_plus||hlt_mu_ph_trigger),0,0,weight_noscale),
    //  TableRow("MC Single or Dilepton (incl noiso) Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_mu_trigger&&(hlt_mu_trigger_plus||hlt_double_mu_noiso_trigger),0,0,weight_noscale),
    //  TableRow("MC Single (incl noiso) or Dilepton (incl noiso) Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_mu_trigger&&(hlt_mu_trigger_plus||hlt_double_mu_noiso_trigger||hlt_single_mu_noiso_trigger),0,0,weight_noscale),
    //  TableRow("MC Single (incl noiso) or Dilepton (incl noiso) or Mu+Ph Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_mu_trigger&&(hlt_mu_trigger_plus||hlt_double_mu_noiso_trigger||hlt_single_mu_noiso_trigger||hlt_mu_ph_trigger),0,0,weight_noscale),
    //  TableRow("MC Single (incl noiso,presc.) or Dilepton (incl noiso) or Mu+Ph Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_mu_trigger&&(hlt_mu_trigger_plus||hlt_double_mu_noiso_trigger||hlt_single_mu_noiso_trigger||hlt_mu_ph_trigger||hlt_single_mu_prescale_trigger),0,0,weight_noscale),
    //  TableRow("Data Single Lepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_mu_trigger,0,0,weight_noscale*eff_singlelep),
    //  TableRow("Data Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_mu_trigger,0,0,weight_noscale*eff_dilep),
    //  TableRow("Data Dilepton Triggers (HIG-19-014)", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_mu_trigger,0,0,weight_noscale*eff_dilep_hig19014),
    //  TableRow("Data Single or Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_mu_trigger,0,0,weight_noscale*eff_singleanddilep),
    //  TableRow("\\hline\n MC Single Lepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&midleadpt&&l1_mu_trigger&&hlt_single_mu_trigger,0,0,weight_noscale),
    //  TableRow("MC Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&midleadpt&&l1_mu_trigger&&hlt_mu_trigger,0,0,weight_noscale),
    //  TableRow("MC Single or Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&midleadpt&&l1_mu_trigger&&hlt_mu_trigger_plus,0,0,weight_noscale),
    //  TableRow("Data Single Lepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&midleadpt&&l1_mu_trigger,0,0,weight_noscale*eff_singlelep),
    //  TableRow("Data Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&midleadpt&&l1_mu_trigger,0,0,weight_noscale*eff_dilep),
    //  TableRow("Data Dilepton Triggers (HIG-19-014)", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&midleadpt&&l1_mu_trigger,0,0,weight_noscale*eff_dilep_hig19014),
    //  TableRow("Data Single or Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&midleadpt&&l1_mu_trigger,0,0,weight_noscale*eff_singleanddilep),
    //  TableRow("\\hline\n MC Single Lepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&highleadpt&&l1_mu_trigger&&hlt_single_mu_trigger,0,0,weight_noscale),
    //  TableRow("MC Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&highleadpt&&l1_mu_trigger&&hlt_mu_trigger,0,0,weight_noscale),
    //  TableRow("MC Single or Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&highleadpt&&l1_mu_trigger&&hlt_mu_trigger_plus,0,0,weight_noscale),
    //  TableRow("Data Single Lepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&highleadpt&&l1_mu_trigger,0,0,weight_noscale*eff_singlelep),
    //  TableRow("Data Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&highleadpt&&l1_mu_trigger,0,0,weight_noscale*eff_dilep),
    //  TableRow("Data Dilepton Triggers (HIG-19-014)", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&highleadpt&&l1_mu_trigger,0,0,weight_noscale*eff_dilep_hig19014),
    //  TableRow("Data Single or Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&highleadpt&&l1_mu_trigger,0,0,weight_noscale*eff_singleanddilep),
    //  TableRow("\\hline\n MC Single Lepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&highpt&&l1_mu_trigger&&hlt_single_mu_trigger,0,0,weight_noscale),
    //  TableRow("MC Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&highpt&&l1_mu_trigger&&hlt_mu_trigger,0,0,weight_noscale),
    //  TableRow("MC Single or Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&highpt&&l1_mu_trigger&&hlt_mu_trigger_plus,0,0,weight_noscale),
    //  TableRow("Data Single Lepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&highpt&&l1_mu_trigger,0,0,weight_noscale*eff_singlelep),
    //  TableRow("Data Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&highpt&&l1_mu_trigger,0,0,weight_noscale*eff_dilep),
    //  TableRow("Data Dilepton Triggers (HIG-19-014)", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&highpt&&l1_mu_trigger,0,0,weight_noscale*eff_dilep_hig19014),
    //  TableRow("Data Single or Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&highpt&&l1_mu_trigger,0,0,weight_noscale*eff_singleanddilep),
    //  TableRow("\\hline\n MC Single Lepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut&&l1_mu_trigger&&hlt_single_mu_trigger,0,0,weight_noscale),
    //  TableRow("MC Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut&&l1_mu_trigger&&hlt_mu_trigger,0,0,weight_noscale),
    //  TableRow("MC Single or Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut&&l1_mu_trigger&&hlt_mu_trigger_plus,0,0,weight_noscale),
    //  TableRow("Data Single Lepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut&&l1_mu_trigger,0,0,weight_noscale*eff_singlelep),
    //  TableRow("Data Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut&&l1_mu_trigger,0,0,weight_noscale*eff_dilep),
    //  TableRow("Data Dilepton Triggers (HIG-19-014)", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut&&l1_mu_trigger,0,0,weight_noscale*eff_dilep_hig19014),
    //  TableRow("Data Single or Dilepton Triggers", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut&&l1_mu_trigger,0,0,weight_noscale*eff_singleanddilep),
    //},signal_procs_noscale,false,true,false,false,false,true).LuminosityTag(total_luminosity_string).Precision(2);
    
    //for (unsigned ieta = 0; ieta < el_abseta_bins.size()-1; ieta++) {
    //  //pm.Push<EfficiencyPlot>(Axis(el_pt_bins, "Electron_pt", "Offline Electron p_{T} [GeV]", {}),
    //  //    "HLT_Ele35_WPTight_Gsf"&&(nSignalElectron==2)&&(ZCandidate_mass>81&&ZCandidate_mass<101)
    //  //    &&Electron_OtherPassReference&&Electron_sig&&Electron_absEta>el_abseta_bins[ieta]
    //  //    &&Electron_absEta<el_abseta_bins[ieta+1]&&Electron_hltId>0.5,
    //  //    Electron_hltPt>12,
    //  //    procs_onelep_data,true,plt_lin_logx).Weight("1").Tag("FixName:zghlt_dataeff_elhltpt12_eta"+std::to_string(ieta)).YTitle("HLT p_{T} > 12 GeV").LuminosityTag(total_luminosity_string);
    //  pm.Push<EfficiencyPlot>(Axis(el_pt_bins_ptcut, "Electron_pt", "Offline Electron p_{T} [GeV]", {}),
    //      "HLT_Ele35_WPTight_Gsf"&&(nSignalElectron==2)&&(ZCandidate_mass>81&&ZCandidate_mass<101)
    //      &&Electron_OtherPassReference&&Electron_sig&&Electron_absEta>el_abseta_bins[ieta]
    //      &&Electron_absEta<el_abseta_bins[ieta+1]&&Electron_hltId>-0.5,
    //      Electron_hltPt>23,
    //      procs_onelep_data,true,plt_lin_logx).Weight("1").Tag("FixName:zghlt_dataeff_elhltpt23_eta"+std::to_string(ieta)).YTitle("HLT p_{T} > 23 GeV").LuminosityTag(total_luminosity_string);
    //  pm.Push<EfficiencyPlot>(Axis(el_pt_bins, "Electron_pt", "Offline Electron p_{T} [GeV]", {}),
    //      "HLT_Ele35_WPTight_Gsf"&&(nSignalElectron==2)&&(ZCandidate_mass>81&&ZCandidate_mass<101)
    //      &&Electron_OtherPassReference&&Electron_OtherPassUpperLeg&&Electron_sig&&!Electron_isLeading
    //      &&Electron_absEta>el_abseta_bins[ieta]&&Electron_absEta<el_abseta_bins[ieta+1],
    //      hlt_el_trigger,
    //      procs_onelep_data,true,plt_lin_logx).Weight("1").Tag("FixName:zghlt_dataeff_dieleleg12_eta"+std::to_string(ieta)).YTitle("Ele12_CaloIdL_TrackIdL_IsoVL").LuminosityTag(total_luminosity_string);
    //  pm.Push<EfficiencyPlot>(Axis(el_pt_bins, "Electron_pt", "Offline Electron p_{T} [GeV]", {}),
    //      "HLT_PFMET120_PFMHT120_IDTight&&MET_pt>150"&&(nSignalElectron==1)&&(Electron_sig>0.5)&&"Electron_mvaFall17V2Iso_WP90"&&Electron_absEta>el_abseta_bins[ieta]
    //      &&Electron_absEta<el_abseta_bins[ieta+1],
    //      "HLT_Ele32_WPTight_Gsf_L1DoubleEG||HLT_Ele35_WPTight_Gsf",
    //      procs_met_data,true,plt_lin_logx).Weight("1").Tag("FixName:zghlt_dataeff_isoel3235_eta"+std::to_string(ieta)).YTitle("Single Electron Triggers").LuminosityTag(total_luminosity_string);
    //}
    //for (unsigned ieta = 0; ieta < mu_abseta_bins.size()-1; ieta++) {
    //}
    //
    //pm.Push<Hist1D>(Axis(25, 0.0, 140.0, Lead_SignalElectron_pt, "Lead e pT", {}),
    //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&hlt_el_trigger,
    //    signal_procs_noscale, plt_lin_over).Weight(weight*ScaleFactor_Triggers)
    //    .Tag("FixName:zghlt__debug_newsf_elead")
    //    .LuminosityTag(total_luminosity_string);
    //pm.Push<Hist1D>(Axis(25, 0.0, 140.0, Lead_SignalElectron_pt, "Lead e pT", {}),
    //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,
    //    signal_procs_noscale, plt_lin_over).Weight(weight*eff_dilep)
    //    .Tag("FixName:zghlt__debug_oldsf_elead")
    //    .LuminosityTag(total_luminosity_string);
    //pm.Push<Hist1D>(Axis(25, 0.0, 140.0, Sublead_SignalElectron_pt, "Sublead e pT", {}),
    //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&hlt_el_trigger,
    //    signal_procs_noscale, plt_lin_over).Weight(weight*ScaleFactor_Triggers)
    //    .Tag("FixName:zghlt__debug_newsf_esubl")
    //    .LuminosityTag(total_luminosity_string);
    //pm.Push<Hist1D>(Axis(25, 0.0, 140.0, Sublead_SignalElectron_pt, "Sublead e pT", {}),
    //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,
    //    signal_procs_noscale, plt_lin_over).Weight(weight*eff_dilep)
    //    .Tag("FixName:zghlt__debug_oldsf_esubl")
    //    .LuminosityTag(total_luminosity_string);
    //pm.Push<Hist1D>(Axis(25, 0.0, 140.0, Lead_SignalMuon_pt, "Lead #mu pT", {}),
    //    (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&hlt_mu_trigger,
    //    signal_procs_noscale, plt_lin_over).Weight(weight*ScaleFactor_Triggers)
    //    .Tag("FixName:zghlt__debug_newsf_mulead")
    //    .LuminosityTag(total_luminosity_string);
    //pm.Push<Hist1D>(Axis(25, 0.0, 140.0, Lead_SignalMuon_pt, "Lead #mu pT", {}),
    //    (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,
    //    signal_procs_noscale, plt_lin_over).Weight(weight*eff_dilep)
    //    .Tag("FixName:zghlt__debug_oldsf_mulead")
    //    .LuminosityTag(total_luminosity_string);
    //pm.Push<Hist1D>(Axis(25, 0.0, 140.0, Sublead_SignalMuon_pt, "Sublead #mu pT", {}),
    //    (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&hlt_mu_trigger,
    //    signal_procs_noscale, plt_lin_over).Weight(weight*ScaleFactor_Triggers)
    //    .Tag("FixName:zghlt__debug_newsf_musubl")
    //    .LuminosityTag(total_luminosity_string);
    //pm.Push<Hist1D>(Axis(25, 0.0, 140.0, Sublead_SignalMuon_pt, "Sublead #mu pT", {}),
    //    (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,
    //    signal_procs_noscale, plt_lin_over).Weight(weight*eff_dilep)
    //    .Tag("FixName:zghlt__debug_oldsf_musubl")
    //    .LuminosityTag(total_luminosity_string);

    //pm.Push<Table>("debug_datamc_trig_eff"+options.year_string, vector<TableRow>{
    //  TableRow("\\hline \\hline\n $Z\\rightarrow e^{+}e^{-}$ decays", 
    //      (z_decay_pdgid==11),0,0,weight_noscale),
    //  TableRow("Offline leptons and photon in acceptance", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1),0,0,weight_noscale),
    //  TableRow("Baseline Selection (No pT cut)", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,0,0,weight_noscale),
    //  TableRow("MC L1 Trigger", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger,0,0,weight_noscale),
    //  TableRow("MC HLT and not L1 Trigger", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&hlt_el_trigger&&!l1_el_trigger,0,0,weight_noscale),
    //  TableRow("Baseline no pt (both barrel)", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(Lead_Electron_abseta<0.8)&&(Sublead_Electron_abseta<0.8),0,0,weight_noscale),
    //  TableRow("Baseline no pt (both endcap)", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(Lead_Electron_abseta>1.566)&&(Lead_Electron_abseta<2.55)&&(Sublead_Electron_abseta>1.566)&&(Sublead_Electron_abseta<2.55),0,0,weight_noscale),
    //  TableRow("Baseline no pt (both barrel) (TSFs)", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&hlt_el_trigger&&(Lead_Electron_abseta<0.8)&&(Sublead_Electron_abseta<0.8),0,0,weight_noscale*ScaleFactor_Triggers),
    //  TableRow("Baseline no pt (both endcap) (TSFs)", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&hlt_el_trigger&&(Lead_Electron_abseta>1.566)&&(Lead_Electron_abseta<2.55)&&(Sublead_Electron_abseta>1.566)&&(Sublead_Electron_abseta<2.55),0,0,weight_noscale*ScaleFactor_Triggers),
    //  TableRow("Baseline no pt (sublmupt > 25)", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(Sublead_SignalElectron_pt>25),0,0,weight_noscale),
    //  TableRow("Baseline no pt (sublmupt > 25) (TSFs)", 
    //      (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&hlt_el_trigger&&(Sublead_SignalElectron_pt>25),0,0,weight_noscale*ScaleFactor_Triggers),


    //  TableRow("\\hline \\hline\n $Z\\rightarrow \\mu^{+}\\mu^{-}$ decays", 
    //      (z_decay_pdgid==13),0,0,weight_noscale),
    //  TableRow("Offline leptons and photon in acceptance", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1),0,0,weight_noscale),
    //  TableRow("Baseline Selection (No pT cut)", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,0,0,weight_noscale),
    //  TableRow("MC L1 Trigger", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_mu_trigger,0,0,weight_noscale),
    //  TableRow("MC HLT and not L1 Trigger", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&hlt_mu_trigger&&!l1_mu_trigger,0,0,weight_noscale),
    //  TableRow("Baseline no pt (both barrel)", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(Lead_Muon_abseta<0.9)&&(Sublead_Muon_abseta<0.9),0,0,weight_noscale),
    //  TableRow("Baseline no pt (both endcap)", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(Lead_Muon_abseta>2.1)&&(Lead_Muon_abseta<2.4)&&(Sublead_Muon_abseta>2.1)&&(Sublead_Muon_abseta<2.4),0,0,weight_noscale),
    //  TableRow("Baseline no pt (both barrel) (TSFs)", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&hlt_mu_trigger&&(Lead_Muon_abseta<0.9)&&(Sublead_Muon_abseta<0.9),0,0,weight_noscale*ScaleFactor_Triggers),
    //  TableRow("Baseline no pt (both endcap) (TSFs)", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&hlt_mu_trigger&&(Lead_Muon_abseta>2.1)&&(Lead_Muon_abseta<2.4)&&(Sublead_Muon_abseta>2.1)&&(Sublead_Muon_abseta<2.4),0,0,weight_noscale*ScaleFactor_Triggers),
    //  TableRow("Baseline no pt (sublmupt > 25)", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(Sublead_SignalMuon_pt>25),0,0,weight_noscale),
    //  TableRow("Baseline no pt (sublmupt > 25) (TSFs)", 
    //      (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&hlt_mu_trigger&&(Sublead_SignalMuon_pt>25),0,0,weight_noscale*ScaleFactor_Triggers),
    //},signal_procs_noscale,false,true,false,false,false,true).LuminosityTag(total_luminosity_string).Precision(2);

    pm.Push<Table>("datamc_trig_eff"+options.year_string, vector<TableRow>{
      TableRow("\\hline \\hline\n $Z\\rightarrow e^{+}e^{-}$ decays", 
          (z_decay_pdgid==11),0,0,weight_noscale),
      TableRow("Offline leptons and photon in acceptance", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1),0,0,weight_noscale),
      TableRow("Baseline Selection (No pT cut)", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,0,0,weight_noscale),
      TableRow("MC L1 Trigger", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger,0,0,weight_noscale),
      TableRow("\\hline MC Single Lepton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&hlt_single_el_trigger,0,0,weight_noscale),
      TableRow("MC Dilepton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&hlt_el_trigger,0,0,weight_noscale),
      TableRow("MC Single or Dilepton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&hlt_el_trigger_plus,0,0,weight_noscale),
      TableRow("\\hline Data Single Lepton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&hlt_single_el_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Data Dilepton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&hlt_el_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Data Dilepton Triggers (Not SF)", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,0,0,weight_noscale*eff_dilep),
      TableRow("Data Dilepton Triggers (HIG-19-014)", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,0,0,weight_noscale*eff_dilep_hig19014),
      TableRow("Data Single or Dilepton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&hlt_el_trigger_plus,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("\\hline Baseline (pT cut)", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection,0,0,weight_noscale),
      TableRow("\\hline MC L1 Trigger (pT cut)", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection&&l1_el_trigger,0,0,weight_noscale),
      TableRow("\\hline MC Single Lepton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection&&hlt_single_el_trigger,0,0,weight_noscale),
      TableRow("MC Dilepton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection&&hlt_el_trigger,0,0,weight_noscale),
      TableRow("MC Single or Dilepton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection&&hlt_el_trigger_plus,0,0,weight_noscale),
      TableRow("\\hline Data L1 Trigger", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("\\hline Data Single Lepton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection&&hlt_single_el_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Data Dilepton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection&&hlt_el_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Data Dilepton Triggers (HIG-19-014)", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection,0,0,weight_noscale*eff_dilep_hig19014),
      TableRow("Data Single or Dilepton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection&&hlt_el_trigger_plus,0,0,weight_noscale*ScaleFactor_Triggers),

      TableRow("\\hline \\hline\n $Z\\rightarrow \\mu^{+}\\mu^{-}$ decays", 
          (z_decay_pdgid==13),0,0,weight_noscale),
      TableRow("Offline leptons and photon in acceptance", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1),0,0,weight_noscale),
      TableRow("Baseline Selection (No pT cut)", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,0,0,weight_noscale),
      TableRow("MC L1 Trigger (only 2l)", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger),0,0,weight_noscale),
      TableRow("MC L1 Trigger", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger||l1_mu_trigger_singlemu),0,0,weight_noscale),
      TableRow("\\hline MC Single Lepton Triggers", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&hlt_single_mu_trigger,0,0,weight_noscale),
      TableRow("MC Dilepton Triggers", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&hlt_mu_trigger,0,0,weight_noscale),
      TableRow("MC Single or Dilepton Triggers", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&hlt_mu_trigger_plus,0,0,weight_noscale),
      TableRow("\\hline Data Single Lepton Triggers", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&hlt_single_mu_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Data Dilepton Triggers", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&hlt_mu_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Data Dilepton Triggers (Not SF)", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,0,0,weight_noscale*eff_dilep),
      TableRow("Data Dilepton Triggers (HIG-19-014)", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,0,0,weight_noscale*eff_dilep_hig19014),
      TableRow("Data Single or Dilepton Triggers", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&hlt_mu_trigger_plus,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("\\hline Baseline (pT cut)", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection,0,0,weight_noscale),
      TableRow("MC L1 Trigger)", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection&&(l1_mu_trigger||l1_mu_trigger_singlemu),0,0,weight_noscale),
      TableRow("\\hline MC Single Lepton Triggers", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection&&hlt_single_mu_trigger,0,0,weight_noscale),
      TableRow("MC Dilepton Triggers", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection&&hlt_mu_trigger,0,0,weight_noscale),
      TableRow("MC Single or Dilepton Triggers", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection&&hlt_mu_trigger_plus,0,0,weight_noscale),
      TableRow("\\hline Data L1 Trigger", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("\\hline Data Single Lepton Triggers", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection&&hlt_single_mu_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Data Dilepton Triggers", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection&&hlt_mu_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Data Dilepton Triggers (HIG-19-014)", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection,0,0,weight_noscale*eff_dilep_hig19014),
      TableRow("Data Single or Dilepton Triggers", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection&&hlt_mu_trigger_plus,0,0,weight_noscale*ScaleFactor_Triggers),
    },signal_procs_noscale,false,true,false,false,false,true).LuminosityTag(total_luminosity_string).Precision(2);

    pm.Push<Table>("datamc_trig_eff_pt_"+options.year_string, vector<TableRow>{
      TableRow("\\hline \\hline\n $Z\\rightarrow e^{+}e^{-}$ decays", 
          (z_decay_pdgid==11),0,0,weight_noscale),
      TableRow("Offline leptons and photon in acceptance", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1),0,0,weight_noscale),
      TableRow("Baseline Selection (No pT cut)", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,0,0,weight_noscale),
      TableRow("MC L1 Trigger", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger,0,0,weight_noscale),
      TableRow("\\hline MC Single Lepton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&hlt_single_el_trigger,0,0,weight_noscale),
      TableRow("MC Dilepton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&hlt_el_trigger,0,0,weight_noscale),
      TableRow("MC Single or Dilepton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&hlt_el_trigger_plus,0,0,weight_noscale),
      TableRow("MC Single or Dilepton Triggers or Diphoton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&(hlt_el_trigger_plus||hlt_doubleph_trigger),0,0,weight_noscale),
      TableRow("\\hline Data Single Lepton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&hlt_single_el_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Data Dilepton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&hlt_el_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Data Dilepton Triggers (Not SF)", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,0,0,weight_noscale*eff_dilep),
      TableRow("Data Dilepton Triggers (HIG-19-014)", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,0,0,weight_noscale*eff_dilep_hig19014),
      TableRow("Data Single or Dilepton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&hlt_el_trigger_plus,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Data Single or Dilepton Triggers or Diphoton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&(hlt_el_trigger_plus||hlt_doubleph_trigger),0,0,weight_noscale*ScaleFactor_Triggers),

      TableRow("MC Single or Dilepton Triggers or 2Ph30PV18PV_pixveto_m55 Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&(hlt_el_trigger_plus||"HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55"),0,0,weight_noscale),
      TableRow("MC Single or Dilepton Triggers or 2Ph30_22_m90 Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&(hlt_el_trigger_plus||"HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90"),0,0,weight_noscale),
      TableRow("Data Single or Dilepton Triggers or 2Ph30PV18PV_pixveto_m55 Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&(hlt_el_trigger_plus||"HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55"),0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Data Single or Dilepton Triggers or 2Ph30_22_m90 Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger&&(hlt_el_trigger_plus||"HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90"),0,0,weight_noscale*ScaleFactor_Triggers),
      //TableRow("MC Single or Dilepton Triggers or 2Ph30_18_nopixveto_m55 Triggers", 
      //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(hlt_el_trigger_plus||"HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55"),0,0,weight_noscale),
      //TableRow("MC Single or Dilepton Triggers or 2Ph30_18_pixveto_m55 Triggers", 
      //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(hlt_el_trigger_plus||"HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55"),0,0,weight_noscale),
      //TableRow("MC Single or Dilepton Triggers or or 2Ph30_22_m90 or 2Ph30_22_m95 Triggers", 
      //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(hlt_el_trigger_plus||"HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90||HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95"),0,0,weight_noscale),
      //TableRow("Data Single or Dilepton Triggers or 2Ph30_18_nopixveto_m55 Triggers", 
      //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(hlt_el_trigger_plus||"HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55"),0,0,weight_noscale*ScaleFactor_Triggers),
      //TableRow("Data Single or Dilepton Triggers or 2Ph30_18_pixveto_m55 Triggers", 
      //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(hlt_el_trigger_plus||"HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55"),0,0,weight_noscale*ScaleFactor_Triggers),
      //TableRow("Data Single or Dilepton Triggers or 2Ph30_22_m95 or 2Ph30_22_m95 Triggers", 
      //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(hlt_el_trigger_plus||"HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90||HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95"),0,0,weight_noscale*ScaleFactor_Triggers),
      //TableRow("\\hline Baseline (pT cut)", 
      //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut,0,0,weight_noscale),
      //TableRow("\\hline MC L1 Trigger (pT cut)", 
      //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut&&l1_el_trigger,0,0,weight_noscale),
      //TableRow("\\hline MC Single Lepton Triggers", 
      //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut&&hlt_single_el_trigger,0,0,weight_noscale),
      //TableRow("MC Dilepton Triggers", 
      //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut&&hlt_el_trigger,0,0,weight_noscale),
      //TableRow("MC Single or Dilepton Triggers", 
      //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut&&hlt_el_trigger_plus,0,0,weight_noscale),
      //TableRow("\\hline Data L1 Trigger", 
      //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut,0,0,weight_noscale*ScaleFactor_Triggers),
      //TableRow("\\hline Data Single Lepton Triggers", 
      //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut&&hlt_single_el_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      //TableRow("Data Dilepton Triggers", 
      //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut&&hlt_el_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      //TableRow("Data Dilepton Triggers (HIG-19-014)", 
      //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut,0,0,weight_noscale*eff_dilep_hig19014),
      //TableRow("Data Single or Dilepton Triggers", 
      //    (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut&&hlt_el_trigger_plus,0,0,weight_noscale*ScaleFactor_Triggers),

      TableRow("\\hline \\hline\n $Z\\rightarrow \\mu^{+}\\mu^{-}$ decays", 
          (z_decay_pdgid==13),0,0,weight_noscale),
      TableRow("Offline leptons and photon in acceptance", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1),0,0,weight_noscale),
      TableRow("Baseline Selection (No pT cut)", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,0,0,weight_noscale),
      TableRow("MC L1 Trigger", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger||l1_mu_trigger_singlemu||l1_mu_eg_trigger),0,0,weight_noscale),
      TableRow("\\hline MC Single Lepton Triggers", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger||l1_mu_trigger_singlemu||l1_mu_eg_trigger)&&hlt_single_mu_trigger,0,0,weight_noscale),
      TableRow("MC Dilepton Triggers", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger||l1_mu_trigger_singlemu||l1_mu_eg_trigger)&&hlt_mu_trigger,0,0,weight_noscale),
      TableRow("MC Single or Dilepton Triggers", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger||l1_mu_trigger_singlemu||l1_mu_eg_trigger)&&hlt_mu_trigger_plus,0,0,weight_noscale),
      TableRow("MC Single or Dilepton Triggers or Muon+Photon Trigger", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger||l1_mu_trigger_singlemu||l1_mu_eg_trigger)&&(hlt_mu_trigger_plus||hlt_mu_ph_trigger),0,0,weight_noscale),
      TableRow("\\hline Data Single Lepton Triggers", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger||l1_mu_trigger_singlemu||l1_mu_eg_trigger)&&hlt_single_mu_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Data Dilepton Triggers", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger||l1_mu_trigger_singlemu||l1_mu_eg_trigger)&&hlt_mu_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Data Dilepton Triggers (Not SF)", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,0,0,weight_noscale*eff_dilep),
      TableRow("Data Dilepton Triggers (HIG-19-014)", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,0,0,weight_noscale*eff_dilep_hig19014),
      TableRow("Data Single or Dilepton Triggers", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger||l1_mu_trigger_singlemu||l1_mu_eg_trigger)&&hlt_mu_trigger_plus,0,0,weight_noscale*ScaleFactor_Triggers),
      TableRow("Data Single or Dilepton Triggers or Muon+Photon Trigger", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger||l1_mu_trigger_singlemu||l1_mu_eg_trigger)&&(hlt_mu_trigger_plus||hlt_mu_ph_trigger),0,0,weight_noscale*ScaleFactor_Triggers),
      //TableRow("\\hline Baseline (pT cut)", 
      //    (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut,0,0,weight_noscale),
      //TableRow("MC L1 Trigger)", 
      //    (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut&&(l1_mu_trigger||l1_mu_trigger_singlemu),0,0,weight_noscale),
      //TableRow("\\hline MC Single Lepton Triggers", 
      //    (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut&&hlt_single_mu_trigger,0,0,weight_noscale),
      //TableRow("MC Dilepton Triggers", 
      //    (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut&&hlt_mu_trigger,0,0,weight_noscale),
      //TableRow("MC Single or Dilepton Triggers", 
      //    (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut&&hlt_mu_trigger_plus,0,0,weight_noscale),
      //TableRow("\\hline Data L1 Trigger", 
      //    (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut,0,0,weight_noscale*ScaleFactor_Triggers),
      //TableRow("\\hline Data Single Lepton Triggers", 
      //    (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut&&hlt_single_mu_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      //TableRow("Data Dilepton Triggers", 
      //    (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut&&hlt_mu_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
      //TableRow("Data Dilepton Triggers (HIG-19-014)", 
      //    (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut,0,0,weight_noscale*eff_dilep_hig19014),
      //TableRow("Data Single or Dilepton Triggers", 
      //    (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&trigger_dependent_ptcut&&hlt_mu_trigger_plus,0,0,weight_noscale*ScaleFactor_Triggers),
    },signal_procs_noscale,false,true,false,false,false,true).LuminosityTag(total_luminosity_string).Precision(2);

    pm.Push<Table>("datamc_trig_eff_2018", vector<TableRow>{
      TableRow("\\hline \\hline\n $Z\\rightarrow e^{+}e^{-}$ decays", 
          (z_decay_pdgid==11),0,0,weight_noscale),
      TableRow("Offline leptons and photon in acceptance", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1),0,0,weight_noscale),
      TableRow("Baseline Selection (No pT cut)", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,0,0,weight_noscale),
      TableRow("MC L1 Trigger", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_eg_trigger_2018,0,0,weight_noscale),
      TableRow("MC Single Lepton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_eg_trigger_2018&&hlt_singleel_trigger_2018,0,0,weight_noscale),
      TableRow("MC Dilepton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_eg_trigger_2018&&hlt_doubleel_trigger_2018,0,0,weight_noscale),
      TableRow("MC Single or Dilepton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_eg_trigger_2018&&(hlt_singleel_trigger_2018||hlt_doubleel_trigger_2018),0,0,weight_noscale),
      TableRow("MC Single or Dilepton Triggers or Diphoton Triggers", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_eg_trigger_2018&&(hlt_singleel_trigger_2018||hlt_doubleel_trigger_2018||hlt_doubleph_trigger_2018),0,0,weight_noscale),
      TableRow("MC Single or Dilepton Triggers or Low Mass Diphoton Trigger", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_eg_trigger_2018&&(hlt_singleel_trigger_2018||hlt_doubleel_trigger_2018||"HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto"),0,0,weight_noscale),
      TableRow("MC Single or Dilepton Triggers or High Mass Diphoton Trigger", 
          (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_eg_trigger_2018&&(hlt_singleel_trigger_2018||hlt_doubleel_trigger_2018||"HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90"),0,0,weight_noscale),

      TableRow("\\hline \\hline\n $Z\\rightarrow \\mu^{+}\\mu^{-}$ decays", 
          (z_decay_pdgid==13),0,0,weight_noscale),
      TableRow("Offline leptons and photon in acceptance", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1),0,0,weight_noscale),
      TableRow("Baseline Selection (No pT cut)", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt,0,0,weight_noscale),
      TableRow("MC L1 Trigger", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger_singlemu||l1_doublemu_trigger_2018||l1_mu_eg_trigger_2018),0,0,weight_noscale),
      TableRow("MC Single Lepton Trigger", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger_singlemu||l1_doublemu_trigger_2018||l1_mu_eg_trigger_2018)&&hlt_singlemu_trigger_2018,0,0,weight_noscale),
      TableRow("MC Double Lepton Trigger", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger_singlemu||l1_doublemu_trigger_2018||l1_mu_eg_trigger_2018)&&hlt_doublemu_trigger_2018,0,0,weight_noscale),
      TableRow("MC Single or Double Lepton Trigger", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger_singlemu||l1_doublemu_trigger_2018||l1_mu_eg_trigger_2018)&&(hlt_singlemu_trigger_2018||hlt_doublemu_trigger_2018),0,0,weight_noscale),
      TableRow("MC Single or Double Lepton Or Muon+Photon Trigger", 
          (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger_singlemu||l1_doublemu_trigger_2018||l1_mu_eg_trigger_2018)&&(hlt_singlemu_trigger_2018||hlt_doublemu_trigger_2018||hlt_mu_ph_trigger),0,0,weight_noscale),
    },signal_procs_2018,false,true,false,false,false,true).LuminosityTag(total_luminosity_string).Precision(2);
  }

  if (make_trigcomp_plots) {
    pm.Push<EfficiencyPlot>(Axis(25, 0.0, 100.0, Lead_SignalElectron_pt, "Lead Electron p_{T} [GeV]", {}),
        (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger,
        hlt_el_trigger,
        signal_procs_noscale,true,plt_lin).Weight(weight_noscale*ScaleFactor_Triggers).Tag("FixName:zghlt_trigcomp_leade_2e").YTitle("Trigger Efficiency").YAxisMax(1.4).LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(25, 0.0, 100.0, Lead_SignalElectron_pt, "Lead Electron p_{T} [GeV]", {}),
        (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger,
        hlt_el_trigger_plus,
        signal_procs_noscale_green,true,plt_lin).Weight(weight_noscale*ScaleFactor_Triggers).Tag("FixName:zghlt_trigcomp_leade_2eor1e").YTitle("Trigger Efficiency").YAxisMax(1.4).LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(25, 0.0, 100.0, Sublead_SignalElectron_pt, "Sublead Electron p_{T} [GeV]", {}),
        (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger,
        hlt_el_trigger,
        signal_procs_noscale,true,plt_lin).Weight(weight_noscale*ScaleFactor_Triggers).Tag("FixName:zghlt_trigcomp_suble_2e").YTitle("Trigger Efficiency").YAxisMax(1.4).LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(25, 0.0, 100.0, Sublead_SignalElectron_pt, "Sublead Electron p_{T} [GeV]", {}),
        (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger,
        hlt_el_trigger_plus,
        signal_procs_noscale_green,true,plt_lin).Weight(weight_noscale*ScaleFactor_Triggers).Tag("FixName:zghlt_trigcomp_suble_2eor1e").YTitle("Trigger Efficiency").YAxisMax(1.4).LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(25, 0.0, 100.0, Lead_SignalMuon_pt, "Lead Muon p_{T} [GeV]", {}),
        (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger||l1_mu_trigger_singlemu),
        hlt_mu_trigger,
        signal_procs_noscale,true,plt_lin).Weight(weight_noscale*ScaleFactor_Triggers).Tag("FixName:zghlt_trigcomp_leadmu_2m").YTitle("Trigger Efficiency").YAxisMax(1.4).LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(25, 0.0, 100.0, Lead_SignalMuon_pt, "Lead Muon p_{T} [GeV]", {}),
        (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger||l1_mu_trigger_singlemu),
        hlt_mu_trigger_plus,
        signal_procs_noscale_green,true,plt_lin).Weight(weight_noscale*ScaleFactor_Triggers).Tag("FixName:zghlt_trigcomp_leadmu_2mor1m").YTitle("Trigger Efficiency").YAxisMax(1.4).LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(25, 0.0, 100.0, Sublead_SignalMuon_pt, "Sublead Muon p_{T} [GeV]", {}),
        (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger||l1_mu_trigger_singlemu),
        hlt_mu_trigger,
        signal_procs_noscale,true,plt_lin).Weight(weight_noscale*ScaleFactor_Triggers).Tag("FixName:zghlt_trigcomp_sublmu_2m").YTitle("Trigger Efficiency").YAxisMax(1.4).LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(25, 0.0, 100.0, Sublead_SignalMuon_pt, "Sublead Muon p_{T} [GeV]", {}),
        (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger||l1_mu_trigger_singlemu),
        hlt_mu_trigger_plus,
        signal_procs_noscale_green,true,plt_lin).Weight(weight_noscale*ScaleFactor_Triggers).Tag("FixName:zghlt_trigcomp_sublmu_2mor1m").YTitle("Trigger Efficiency").YAxisMax(1.4).LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(15, 110.0, 140.0, HiggsCandidate_mass, "m_{ll#gamma} [GeV]", {}),
        (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger,
        hlt_el_trigger,
        signal_procs_noscale,true,plt_lin).Weight(weight_noscale*ScaleFactor_Triggers).Tag("FixName:zghlt_trigcomp_mllg_2e").YTitle("Trigger Efficiency").YAxisMax(1.4).LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(15, 110.0, 140.0, HiggsCandidate_mass, "m_{ll#gamma} [GeV]", {}),
        (z_decay_pdgid==11)&&(nSignalElectron>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&l1_el_trigger,
        hlt_el_trigger_plus,
        signal_procs_noscale_green,true,plt_lin).Weight(weight_noscale*ScaleFactor_Triggers).Tag("FixName:zghlt_trigcomp_mllg_2eor1e").YTitle("Trigger Efficiency").YAxisMax(1.4).LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(15, 110.0, 140.0, HiggsCandidate_mass, "m_{ll#gamma} [GeV]", {}),
        (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger||l1_mu_trigger_singlemu),
        hlt_mu_trigger,
        signal_procs_noscale,true,plt_lin).Weight(weight_noscale*ScaleFactor_Triggers).Tag("FixName:zghlt_trigcomp_mllg_2m").YTitle("Trigger Efficiency").YAxisMax(1.4).LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(15, 110.0, 140.0, HiggsCandidate_mass, "m_{ll#gamma} [GeV]", {}),
        (z_decay_pdgid==13)&&(nSignalMuon>=2)&&(nSignalPhoton>=1)&&baseline_selection_nopt&&(l1_mu_trigger||l1_mu_trigger_singlemu),
        hlt_mu_trigger_plus,
        signal_procs_noscale_green,true,plt_lin).Weight(weight_noscale*ScaleFactor_Triggers).Tag("FixName:zghlt_trigcomp_mllg_2mor1m").YTitle("Trigger Efficiency").YAxisMax(1.4).LuminosityTag(total_luminosity_string);
  }

  if (make_mu17ref) {
    pm.Push<Table>("mu17_releff"+options.year_string, vector<TableRow>{
      TableRow("Data (arbitrary lumi) events with Z and passing Mu17", 
          "HLT_Mu17"&&(nSignalMuon==2)&&(ZCandidate_mass>81&&ZCandidate_mass<101),0,0,"1"),
      TableRow("Pass Single Muon Trigger", 
          "HLT_Mu17"&&(nSignalMuon==2)&&(ZCandidate_mass>81&&ZCandidate_mass<101)&&hlt_single_mu_trigger,0,0,"1"),
      TableRow("Pass Dimuon Trigger", 
          "HLT_Mu17"&&(nSignalMuon==2)&&(ZCandidate_mass>81&&ZCandidate_mass<101)&&hlt_mu_trigger,0,0,"1"),
      TableRow("Pass Single or Dimuon Trigger", 
          "HLT_Mu17"&&(nSignalMuon==2)&&(ZCandidate_mass>81&&ZCandidate_mass<101)&&(hlt_mu_trigger||hlt_single_mu_trigger),0,0,"1"),
    },procs_onelep_data,false,true,false,false,false,true).LuminosityTag(total_luminosity_string).Precision(2);
  }

  if (make_dy_efftable) {
    pm.Push<Table>("dy_mc_sf_table_"+options.year_string, vector<TableRow>{
      TableRow("\\hline \\hline\n $Z\\rightarrow e^{+}e^{-}$ decays", 
          (z_decay_pdgid==11),0,0,weight_noscale),
      TableRow("DY selection", 
          (z_decay_pdgid==11)&&dy_selection,0,0,weight_noscale),
      TableRow("L1 Trigger", 
          (z_decay_pdgid==11)&&dy_selection&&l1_el_trigger,0,0,weight_noscale),
      TableRow("HLT Dielectron Trigger", 
          (z_decay_pdgid==11)&&dy_selection&&l1_el_trigger,0,0,weight_noscale*eff_dilep),
      TableRow("HLT Dielectron Trigger (HIG-19-014)", 
          (z_decay_pdgid==11)&&dy_selection&&l1_el_trigger,0,0,weight_noscale*eff_dilep_hig19014),
      TableRow("HLT 2 or 1 Electron Trigger", 
          (z_decay_pdgid==11)&&dy_selection&&l1_el_trigger,0,0,weight_noscale*eff_singleanddilep),

      TableRow("\\hline \\hline\n $Z\\rightarrow \\mu^{+}\\mu^{-}$ decays", 
          (z_decay_pdgid==13),0,0,weight_noscale),
      TableRow("DY selection", 
          (z_decay_pdgid==13)&&dy_selection,0,0,weight_noscale),
      TableRow("L1 Trigger", 
          (z_decay_pdgid==13)&&dy_selection&&l1_mu_trigger,0,0,weight_noscale),
      TableRow("HLT Dimuon Trigger", 
          (z_decay_pdgid==13)&&dy_selection&&l1_mu_trigger,0,0,weight_noscale*eff_dilep),
      TableRow("HLT Dimuon Trigger (HIG-19-014)", 
          (z_decay_pdgid==13)&&dy_selection&&l1_mu_trigger,0,0,weight_noscale*eff_dilep_hig19014),
      TableRow("HLT 2 or 1 Muon Trigger", 
          (z_decay_pdgid==13)&&dy_selection&&l1_mu_trigger,0,0,weight_noscale*eff_singleanddilep),
    },dy_procs,false,true,false,false,false,true).LuminosityTag(total_luminosity_string).Precision(2);
    //
    //pm.Push<Table>("dy_data_sf_table_"+options.year_string, vector<TableRow>{
    //  TableRow("DY selection", 
    //      (nSignalElectron==2)&&dy_selection,0,0,"1"),
    //  TableRow("HLT 2e Trigger", 
    //      (nSignalElectron==2)&&dy_selection&&hlt_el_trigger,0,0,"1"),
    //  TableRow("HLT 1e Trigger", 
    //      (nSignalElectron==2)&&dy_selection&&hlt_single_el_trigger,0,0,"1"),
    //  TableRow("HLT 2e OR 1e Trigger", 
    //      (nSignalElectron==2)&&dy_selection&&hlt_el_trigger_plus,0,0,"1"),
    //  TableRow("HLT 2e OR 1e OR 1e prescaled Trigger", 
    //      (nSignalElectron==2)&&dy_selection&&(hlt_el_trigger_plus||hlt_single_el_prescale_trigger),0,0,"1"),
    //  TableRow("HLT 2e (incl non-iso) OR 1e Trigger", 
    //      (nSignalElectron==2)&&dy_selection&&(hlt_el_trigger_plus||hlt_double_el_noiso_trigger),0,0,"1"),
    //  TableRow("HLT 2e (incl non-iso) OR 1e (incl non-iso) Trigger", 
    //      (nSignalElectron==2)&&dy_selection&&(hlt_el_trigger_plus||hlt_double_el_noiso_trigger||hlt_single_el_noiso_trigger),0,0,"1"),
    //  TableRow("HLT 2e (incl non-iso) OR 1e (incl non-iso) OR 1e prescaled Trigger", 
    //      (nSignalElectron==2)&&dy_selection&&(hlt_el_trigger_plus||hlt_double_el_noiso_trigger||hlt_single_el_noiso_trigger||hlt_single_el_prescale_trigger),0,0,"1"),

    //  TableRow("DY selection", 
    //      (nSignalMuon==2)&&dy_selection,0,0,"1"),
    //  TableRow("HLT 2mu Trigger", 
    //      (nSignalMuon==2)&&dy_selection&&hlt_mu_trigger,0,0,"1"),
    //  TableRow("HLT 1mu Trigger", 
    //      (nSignalMuon==2)&&dy_selection&&hlt_single_mu_trigger,0,0,"1"),
    //  TableRow("HLT 2mu OR 1mu Trigger", 
    //      (nSignalMuon==2)&&dy_selection&&hlt_mu_trigger_plus,0,0,"1"),
    //  TableRow("HLT 2mu OR 1mu OR mu+ph Trigger", 
    //      (nSignalMuon==2)&&dy_selection&&(hlt_mu_trigger_plus||hlt_mu_ph_trigger),0,0,"1"),
    //  TableRow("HLT 2mu OR 1mu OR 1mu prescaled Trigger", 
    //      (nSignalMuon==2)&&dy_selection&&(hlt_mu_trigger_plus||hlt_single_mu_prescale_trigger),0,0,"1"),
    //  TableRow("HLT 2mu (incl non-iso) OR 1mu Trigger", 
    //      (nSignalMuon==2)&&dy_selection&&(hlt_mu_trigger_plus||hlt_double_mu_noiso_trigger),0,0,"1"),
    //  TableRow("HLT 2mu (incl non-iso) OR 1mu (incl non-iso) Trigger", 
    //      (nSignalMuon==2)&&dy_selection&&(hlt_mu_trigger_plus||hlt_double_mu_noiso_trigger||hlt_single_mu_noiso_trigger),0,0,"1"),
    //  TableRow("HLT 2mu (incl non-iso) OR 1mu (incl non-iso) OR mu+ph Trigger", 
    //      (nSignalMuon==2)&&dy_selection&&(hlt_mu_trigger_plus||hlt_double_mu_noiso_trigger||hlt_single_mu_noiso_trigger||hlt_mu_ph_trigger),0,0,"1"),
    //  TableRow("HLT 2mu (incl non-iso) OR 1mu (incl non-iso) OR mu+ph OR 1mu prescaled Trigger", 
    //      (nSignalMuon==2)&&dy_selection&&(hlt_mu_trigger_plus||hlt_double_mu_noiso_trigger||hlt_single_mu_noiso_trigger||hlt_mu_ph_trigger||hlt_single_mu_prescale_trigger),0,0,"1"),
    //},procs_run2017d,false,true,false,false,false,true).LuminosityTag(total_luminosity_string).Precision(2);
  }

  if (plot_data_sidebands) {
    //pm.Push<Hist1D>(Axis(30, 100.0, 160.0, HiggsCandidate_mass, "m_{ll#gamma} [GeV]", {}),
    //    pass_hemveto&&nSignalElectron>=2&&nSignalPhoton>=1&&baseline_selection_nopt&&Lead_Photon_mvaID>0.5&&Sublead_SignalElectron_pt>15&&hlt_doubleel_trigger_2018&&blind_sr,
    //    procs_sideband_2018, plt_lin_unblind).Weight(weight)
    //    .Tag("FixName:zghlt_diph__mllg_baseline_doubleel")
    //    .LuminosityTag("59.7");
    //pm.Push<Hist1D>(Axis(30, 100.0, 160.0, HiggsCandidate_mass, "m_{ll#gamma} [GeV]", {}),
    //    pass_hemveto&&nSignalElectron>=2&&nSignalPhoton>=1&&baseline_selection_nopt&&Lead_Photon_mvaID>0.5&&Sublead_SignalElectron_pt>15&&(hlt_doubleel_trigger_2018||hlt_singleel_trigger_2018)&&blind_sr,
    //    procs_sideband_2018, plt_lin_unblind).Weight(weight)
    //    .Tag("FixName:zghlt_diph__mllg_baseline_doubleel_or_singleel")
    //    .LuminosityTag("59.7");
    //pm.Push<Hist1D>(Axis(30, 100.0, 160.0, HiggsCandidate_mass, "m_{ll#gamma} [GeV]", {}),
    //    pass_hemveto&&nSignalElectron>=2&&nSignalPhoton>=1&&baseline_selection_nopt&&Lead_Photon_mvaID>0.5&&Sublead_SignalElectron_pt>15&&(hlt_doubleel_trigger_2018||hlt_singleel_trigger_2018||hlt_doubleph_trigger_2018)&&blind_sr,
    //    procs_sideband_2018, plt_lin_unblind).Weight(weight)
    //    .Tag("FixName:zghlt_diph__mllg_baseline_doubleel_or_singleel_ordoubleph")
    //    .LuminosityTag("59.7");
    //pm.Push<Hist1D>(Axis(30, 100.0, 160.0, HiggsCandidate_mass, "m_{ll#gamma} [GeV]", {}),
    //    nSignalMuon>=2&&nSignalPhoton>=1&&baseline_selection_nopt&&Sublead_SignalMuon_pt>10&&hlt_doublemu_trigger_2018&&blind_sr,
    //    procs_sideband_2018, plt_lin).Weight(weight)
    //    .Tag("FixName:zghlt_diph__mllg_baseline_doublemu")
    //    .LuminosityTag(total_luminosity_string);
    //pm.Push<Hist1D>(Axis(30, 100.0, 160.0, HiggsCandidate_mass, "m_{ll#gamma} [GeV]", {}),
    //    nSignalMuon>=2&&nSignalPhoton>=1&&baseline_selection_nopt&&Sublead_SignalMuon_pt>10&&(hlt_doublemu_trigger_2018||hlt_single_mu_trigger)&&blind_sr,
    //    procs_sideband_2018, plt_lin).Weight(weight)
    //    .Tag("FixName:zghlt_diph__mllg_baseline_doublemu_or_singlemu")
    //    .LuminosityTag(total_luminosity_string);
    //pm.Push<Hist1D>(Axis(30, 100.0, 160.0, HiggsCandidate_mass, "m_{ll#gamma} [GeV]", {}),
    //    nSignalMuon>=2&&nSignalPhoton>=1&&baseline_selection_nopt&&Sublead_SignalMuon_pt>10&&(hlt_doublemu_trigger_2018||hlt_single_mu_trigger||hlt_mu_ph_trigger)&&blind_sr,
    //    procs_sideband_2018, plt_lin).Weight(weight)
    //    .Tag("FixName:zghlt_diph__mllg_baseline_doublemu_or_singlemu_ormuph")
    //    .LuminosityTag(total_luminosity_string);

    NamedFunc zg_baseline_noleppt = NamedFunc( nSignalLepton>=2 && nSignalPhoton>=1 &&
        ((Lead_SignalPhoton_pt/HiggsCandidate_mass)>=15.0/110.0) &&
        ((HiggsCandidate_mass+ZCandidate_mass)>185.0) &&
        (Sublead_SignalElectron_pt > 15 || Sublead_SignalMuon_pt > 10) &&
        pass_hemveto).Name("baseline_noleppt");
    NamedFunc near_mass_window = (HiggsCandidate_mass>110&&HiggsCandidate_mass<140);

    pm.Push<Table>("zghlt_diph__eff", vector<TableRow>{
      //TableRow("Electron baseline (no pt) with 2e trigger", 
      //  pass_hemveto&&nSignalElectron>=2&&nSignalPhoton>=1&&baseline_selection_nopt&&Lead_Photon_mvaID>0.5&&Sublead_SignalElectron_pt>15&&hlt_doubleel_trigger_2018&&blind_sr,0,0,weight),
      //TableRow("Electron baseline (no pt) with 2e or 1e trigger", 
      //  pass_hemveto&&nSignalElectron>=2&&nSignalPhoton>=1&&baseline_selection_nopt&&Lead_Photon_mvaID>0.5&&Sublead_SignalElectron_pt>15&&(hlt_doubleel_trigger_2018||hlt_singleel_trigger_2018)&&blind_sr,0,0,weight),
      //TableRow("Electron baseline (no pt) with 2e or 1e or 2ph trigger", 
      //  pass_hemveto&&nSignalElectron>=2&&nSignalPhoton>=1&&baseline_selection_nopt&&Lead_Photon_mvaID>0.5&&Sublead_SignalElectron_pt>15&&(hlt_doubleel_trigger_2018||hlt_singleel_trigger_2018||hlt_doubleph_trigger_2018)&&blind_sr,0,0,weight),
      //TableRow("Muon baseline (no pt) with 2mu trigger", 
      //  nSignalMuon>=2&&nSignalPhoton>=1&&baseline_selection_nopt&&Sublead_SignalMuon_pt>10&&hlt_doublemu_trigger_2018&&blind_sr,0,0,weight),
      //TableRow("Muon baseline (no pt) with 2mu or 1mu trigger", 
      //  nSignalMuon>=2&&nSignalPhoton>=1&&baseline_selection_nopt&&Sublead_SignalMuon_pt>10&&(hlt_doublemu_trigger_2018||hlt_single_mu_trigger)&&blind_sr,0,0,weight),
      //TableRow("Muon baseline (no pt) with 2mu or 1mu or mu+ph trigger", 
      //  nSignalMuon>=2&&nSignalPhoton>=1&&baseline_selection_nopt&&Sublead_SignalMuon_pt>10&&(hlt_doublemu_trigger_2018||hlt_single_mu_trigger||hlt_mu_ph_trigger)&&blind_sr,0,0,weight),
      
      TableRow("Electron baseline (no pt) with 2e trigger (110-140)", 
          zg_baseline_noleppt&&nSignalElectron>=2&&near_mass_window&&HLT_pass_dilepton,0,0,"1"),
      TableRow("Electron baseline (no pt) with 2e or 1e trigger (110-140)", 
          zg_baseline_noleppt&&nSignalElectron>=2&&near_mass_window&&(HLT_pass_dilepton||HLT_pass_singlelepton),0,0,"1"),
      TableRow("Electron baseline (no pt) with 2e or 1e or 2ph trigger (110-140)", 
          zg_baseline_noleppt&&nSignalElectron>=2&&near_mass_window&&(HLT_pass_dilepton||HLT_pass_singlelepton||HLT_pass_diphoton),0,0,"1"),
      TableRow("Muon baseline (no pt) with 2mu trigger (110-140)", 
          zg_baseline_noleppt&&nSignalMuon>=2&&near_mass_window&&HLT_pass_dilepton,0,0,"1"),
      TableRow("Muon baseline (no pt) with 2mu or 1mu trigger (110-140)", 
          zg_baseline_noleppt&&nSignalMuon>=2&&near_mass_window&&(HLT_pass_dilepton||HLT_pass_singlelepton),0,0,"1"),
      TableRow("Muon baseline (no pt) with 2mu or 1mu or mu+ph trigger (110-140)", 
          zg_baseline_noleppt&&nSignalMuon>=2&&near_mass_window&&(HLT_pass_dilepton||HLT_pass_singlelepton||HLT_pass_muonphoton),0,0,"1")
    },procs_data_sideband,false,true,false,false,false,true).LuminosityTag("138").Precision(3);
  }
  
  pm.multithreaded_ = !options.single_thread;
  pm.min_print_ = true; //debugging time baby
  pm.MakePlots(1.);

  //------------------------------------------------------------------------------------
  //                                     post processing
  //------------------------------------------------------------------------------------

  if (plot_trigeffs) {
    std::string outfile_name = "src/zgamma/apply_zg_trigeffs_test.cpp";
    //generate_2d_efficiencies(&pm, "eff_elhltpt12", histnames_elhltpt12, 
    //    outfile_name, "pt", el_pt_bins, 
    //    "abseta", el_abseta_bins);
    generate_2d_efficiencies(&pm, "eff_elhltpt23", histnames_elhltpt23, 
        outfile_name, "pt", el_pt_bins_ptcut, 
        "abseta", el_abseta_bins);
    generate_2d_efficiencies(&pm, "eff_dielleg12", histnames_dieleleg12, 
        outfile_name, "pt", el_pt_bins, 
        "abseta", el_abseta_bins);
    //generate_2d_efficiencies(&pm, "eff_muhltpt8", histnames_muhltpt8, 
    //    outfile_name, "pt", mu_pt_bins, 
    //    "abseta", mu_abseta_bins);
    generate_2d_efficiencies(&pm, "eff_muhltpt17", histnames_muhltpt17, 
        outfile_name, "pt", mu_pt_bins_ptcut, 
        "abseta", mu_abseta_bins);
    generate_2d_efficiencies(&pm, "eff_dimuleg8", histnames_dimuleg8, 
        outfile_name, "pt", mu_pt_bins, 
        "abseta", mu_abseta_bins);
    generate_2d_efficiencies(&pm, "eff_isoel3235", histnames_isoel3235, 
        outfile_name, "pt", el_pt_bins, 
        "abseta", el_abseta_bins);
    generate_2d_efficiencies(&pm, "eff_isomu2427", histnames_isomu2427, 
        outfile_name, "pt", mu_pt_bins, 
        "abseta", mu_abseta_bins);
    //generate_2d_efficiencies(&pm, "eff_mu17", histnames_mu17, 
    //    "src/zgamma/apply_zg_trigeffs.cpp", "pt", mu_pt_bins, 
    //    "abseta", mu_abseta_bins);
  }

  if (plot_trigeffs_signal) {
    std::string outfile_name = "src/zgamma/apply_zg_trigeffs_test.cpp";
    generate_2d_efficiencies(&pm, "effsig_elhltpt23", histnames_sigelhltpt23, 
        outfile_name, "pt", el_pt_bins_ptcut, 
        "abseta", el_abseta_bins);
    generate_2d_efficiencies(&pm, "effsig_dielleg12", histnames_sigdieleleg12, 
        outfile_name, "pt", el_pt_bins, 
        "abseta", el_abseta_bins);
    //generate_2d_efficiencies(&pm, "effsig_dielleg23", histnames_sigdieleleg23, 
    //    outfile_name, "pt", el_pt_bins, 
    //    "abseta", el_abseta_bins);
    generate_2d_efficiencies(&pm, "effsig_muhltpt17", histnames_sigmuhltpt17, 
        outfile_name, "pt", mu_pt_bins_ptcut, 
        "abseta", mu_abseta_bins);
    generate_2d_efficiencies(&pm, "effsig_dimuleg8", histnames_sigdimuleg8, 
        outfile_name, "pt", mu_pt_bins, 
        "abseta", mu_abseta_bins);
    //generate_2d_efficiencies(&pm, "effsig_dimuleg17", histnames_sigdimuleg17, 
    //    outfile_name, "pt", mu_pt_bins, 
    //    "abseta", mu_abseta_bins);
    generate_2d_efficiencies(&pm, "effsig_isoel3235", histnames_sigisoel3235, 
        outfile_name, "pt", el_pt_bins, 
        "abseta", el_abseta_bins);
    generate_2d_efficiencies(&pm, "effsig_isomu2427", histnames_sigisomu2427, 
        outfile_name, "pt", mu_pt_bins, 
        "abseta", mu_abseta_bins);
    //generate_2d_efficiencies(&pm, "eff_mu17", histnames_mu17, 
    //    "src/zgamma/apply_zg_trigeffs.cpp", "pt", mu_pt_bins, 
    //    "abseta", mu_abseta_bins);
  }

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

//helper function definitions
/*! Output efficiencies in C++ format
 
  \param[in] pm pointer to PlotMaker that made efficiency plots
  \param[in] func_name name of function in output file
  \param[in] plot_names vector of names of efficiency plots
  \param[in] out_file_name name of file to append to
 */
void generate_2d_efficiencies(PlotMaker* pm, std::string func_name, 
    std::vector<std::string> plot_names, std::string out_file_name, 
    std::string x_var_name, std::vector<double> x_var_bins,
    std::string y_var_name, std::vector<double> y_var_bins) {
  ofstream out_file;
  out_file.open(out_file_name, std::ios_base::app);
  out_file << "\n";
  out_file << "std::vector<float> " << func_name << "(float " << x_var_name << ", float " << y_var_name << ") {\n";
  out_file << "  float errup=0.0, errdown=0.0, eff=1.0;\n";
  bool first_y = true;
  for (unsigned iy = 0; iy < y_var_bins.size()-1; iy++) {
    EfficiencyPlot * eff_plot = static_cast<EfficiencyPlot*>(pm->GetFigure(plot_names[iy]).get());
    TGraphAsymmErrors* raw_eff_plot = nullptr;
    if (eff_plot->data_ratio_plots_.size()>0) {
      raw_eff_plot = static_cast<TGraphAsymmErrors*>((eff_plot->data_ratio_plots_[0]).get()->Clone());
    }
    else {
      raw_eff_plot = static_cast<TGraphAsymmErrors*>((eff_plot->signal_ratio_plots_[0]).get()->Clone());
    }
    out_file << "  ";
    if (first_y) first_y = false;
    else  out_file << "} else ";
    out_file << "if (" << y_var_name << " > " << y_var_bins[iy] << " && ";
    out_file << y_var_name << " < " << y_var_bins[iy+1] << ") {\n";
    double* eff_x_vals = raw_eff_plot->GetY();
    bool first_x = true;
    for (unsigned ix = 0; ix < x_var_bins.size()-1; ix++) {
      double x_upper = ix==(x_var_bins.size()-2) ? 9999 : x_var_bins[ix+1];
      out_file << "    ";
      if (first_x) first_x = false;
      else out_file << "else ";
      out_file << "if (" << x_var_name << " > " << x_var_bins[ix] << " && ";
      out_file << x_var_name << " < " << x_upper << ")";
      out_file << " {eff = " << eff_x_vals[ix] << "; errup = " << raw_eff_plot->GetErrorYhigh(ix) << "; errdown = " << raw_eff_plot->GetErrorYlow(ix) << ";}\n";
    }
    delete raw_eff_plot;
  }
  out_file << "  }\n";
  out_file << "  return std::vector<float>({eff, errup, errdown});\n";
  out_file << "}\n";
  out_file.close();
}
