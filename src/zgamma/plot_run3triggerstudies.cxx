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

#include "core/test.hpp"

#include <algorithm>
#include <bitset>
#include <cstdlib>
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
#include "core/palette.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/process.hpp"
#include "core/table.hpp"
#include "core/utilities.hpp"
//#include "higgsino/script_utilities.hpp"
//#include "higgsino/hig_functions.hpp"

using namespace std;
using namespace PlotOptTypes;

//manuall bring this in for now
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

int main(int argc, char *argv[]){
  //------------------------------------------------------------------------------------
  //                                    constants
  //------------------------------------------------------------------------------------
  const std::vector<double> el_pt_bins = {5.0, 7.0, 10.0, 15.0, 17.0, 20.0, 25.0, 27.0, 30.0, 40.0, 50.0, 80.0, 120.0, 200.0, 300.0};
  const std::vector<double> el_abseta_bins = {0.0, 0.8, 1.4442, 1.566, 2.5};
  const std::vector<double> mu_pt_bins = {5.0, 8.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 120.0, 200.0};
  const std::vector<double> mu_abseta_bins = {0.0, 0.9, 1.2, 2.1, 2.4};

  //------------------------------------------------------------------------------------
  //                                    initialization
  //------------------------------------------------------------------------------------
  gErrorIgnoreLevel = 6000;
  time_t begtime, endtime;
  time(&begtime);
  Palette colors("txt/colors.txt","default");
  script_utilities::ArgStruct options = script_utilities::get_options(
      argc, argv, "");
  std::vector<PlotOpt> plt_lin = script_utilities::plot_lin(false);
  std::vector<PlotOpt> plt_log = script_utilities::plot_log(false);
  std::vector<PlotOpt> plt_shapes = script_utilities::plot_shapes();
  std::vector<PlotOpt> plt_log_shapes = script_utilities::plot_log_shapes();
  PlotOpt style2D("txt/plot_styles.txt", "Scatter");
  vector<PlotOpt> twodim_plotopts = {style2D().Title(TitleType::info).YAxis(YAxisType::linear).Overflow(OverflowType::overflow).LogMinimum(0.001)};

  //set<int> years;
  //HigUtilities::parseYears(options.year_string, years);
  //int lumi_precision = 0;
  //string total_luminosity_string = HigUtilities::getLuminosityString(options.year_string, lumi_precision);
  std::string total_luminosity_string = "41.5";

  std::string base_folder = "/net/cms17/cms17r0/pico/NanoAODv7/nano/2017/signal/";
  std::string base_folder_ul = "/net/cms17/cms17r0/pico/NanoAODv2/nano/2017/";
  std::string base_folder_ulv9 = "/net/cms17/cms17r0/pico/NanoAODv9/nano/2016/";
  std::string base_folder_v7data = "/net/cms17/cms17r0/pico/NanoAODv7/nano/2017/data/";

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

  std::vector<std::shared_ptr<Process>> full_procs;
  //full_procs.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+jets", 
  //    Process::Type::background, TColor::GetColor("#ffb400"),
  //    {base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-"
  //    "pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000*"},"1"));
  //full_procs.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+#gamma", 
  //    Process::Type::background, TColor::GetColor("#16bac5"),
  //    {base_folder_ul+"mc/ZGToLLG_01J_5f_TuneCP5_*"},"1"));
  full_procs.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+jets", 
      Process::Type::background, TColor::GetColor("#ffb400"),
      {base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__035CF0D8-DDC5-6D46-B172-3BDC98CC624D*",
      base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__223A6F1C-C8C0-9F43-94C3-233AD51ABD32*",
      base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__29EF0133-B16A-E747-9E02-89BC6682B669*", //temp?
      base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__25118058-0F30-8A41-B82B-C7EAFEE55379*", //temp?
      base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__24950148-3722-B84D-9050-C996F5205F0A*"},"1"));
  full_procs.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+#gamma", 
      Process::Type::background, TColor::GetColor("#16bac5"),
      {base_folder_ul+"mc/ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__230000__03E636EE-BAB0-CA42-BFB1-28E772F2D53E*"},"1"));
  full_procs.push_back(Process::MakeShared<Baby_nano>("gg#rightarrow H#rightarrow ZG (x100)", 
      Process::Type::signal, TColor::GetColor("#ff0000"),
      {base_folder_ul+"signal/GluGluHToZG*"},"1"));

  std::vector<std::shared_ptr<Process>> dy_procs;
  dy_procs.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+jets", 
      Process::Type::background, TColor::GetColor("#ffb400"),
      {base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__035CF0D8-DDC5-6D46-B172-3BDC98CC624D*",
      base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__223A6F1C-C8C0-9F43-94C3-233AD51ABD32*",
      base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__29EF0133-B16A-E747-9E02-89BC6682B669*", //temp?
      base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__25118058-0F30-8A41-B82B-C7EAFEE55379*", //temp?
      base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__24950148-3722-B84D-9050-C996F5205F0A*"},"1"));

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
  
  std::vector<std::shared_ptr<Process>> signal_procs;
  signal_procs.push_back(Process::MakeShared<Baby_nano>("gg#rightarrow H#rightarrow ZG (x100)", 
      Process::Type::background, TColor::GetColor("#ff0000"),
      {base_folder_ul+"signal/GluGluHToZG*"},"1"));

  std::vector<std::shared_ptr<Process>> signal_procs_noscale;
  signal_procs_noscale.push_back(Process::MakeShared<Baby_nano>("gg#rightarrow H#rightarrow ZG", 
      Process::Type::signal, TColor::GetColor("#ff0000"),
      {base_folder_ul+"signal/GluGluHToZG*"},"1"));

  std::vector<std::shared_ptr<Process>> signal_procs_noscale_2d;
  signal_procs_noscale_2d.push_back(Process::MakeShared<Baby_nano>("gg#rightarrow H#rightarrow ZG", 
      Process::Type::background, TColor::GetColor("#ff0000"),
      {base_folder_ul+"signal/GluGluHToZG*"},"1"));

  std::vector<std::shared_ptr<Process>> procs_onelep_data;
  procs_onelep_data.push_back(Process::MakeShared<Baby_nano>("Data", 
      Process::Type::data, kBlack,
      {base_folder_v7data+"SingleElectron__Run2017B__02Apr2020-v1__230000__014FD4FF-A181-2843-B548-BB3198F16E88*",
      base_folder_v7data+"SingleElectron__Run2017B__02Apr2020-v1__230000__08B5DC81-D780-A54B-80FE-C94EBD267ACA*",
      base_folder_v7data+"SingleElectron__Run2017B__02Apr2020-v1__230000__0ACA0341-E365-2544-99F0-FBA4C92C0301*",
      base_folder_v7data+"SingleElectron__Run2017B__02Apr2020-v1__230000__0C083B57-67B2-FF43-AF2E-56C850A094D6*",
      base_folder_v7data+"SingleMuon__Run2017B__02Apr2020-v1__230000__0C07D411-3FE0-974E-A5A4-25BB763E4038*",
      base_folder_v7data+"SingleMuon__Run2017B__02Apr2020-v1__230000__0D1DA690-909B-F74F-8DB7-48D6F786F8D8*",
      base_folder_v7data+"SingleMuon__Run2017B__02Apr2020-v1__230000__12EF443F-9CD8-7842-9563-2BBD627A28FB*",
      base_folder_v7data+"SingleMuon__Run2017B__02Apr2020-v1__230000__1A14183C-9C11-E947-9D71-6B89994127D9*"},"1"));

  //------------------------------------------------------------------------------------
  //                                     named funcs
  //------------------------------------------------------------------------------------
  
  const NamedFunc weight("weight_nf",[](const Baby &b) -> NamedFunc::ScalarType{
    //extra factor for weighting is xs/nevts_effective where effective events are calculated including negative weights
    //0.007519850838359999 GluGluH->ZG->llG xs (pb)
    //967054
    //double xs = 7.7760403e-09; //GluGluHToZG
    double xs = 0.00751985;
    double nents = 400000.0/100.0; //artificial x100
    if (b.FirstFileName().find("ZGToLLG") != std::string::npos) {
      xs = 117.864;
      //nents = 2.8095231e+09;
      nents = 1240092;
    }
    if (b.FirstFileName().find("DYJetsToLL") != std::string::npos) {
      xs = 6077.22;
      //nents = 3.2011619e+11;
      //nents = 3535863;
      nents = 7680238;
      if (b.SampleType()==2016) nents = 3865459; 
    }
    //float w_lumi = 1.0;
    float w_year = 1.0;
    //if (b.Generator_weight()<0) w_lumi = -1.0;
    if (b.SampleType()==2016) w_year = 36.32264; 
    else if (b.SampleType()==2017) w_year = 41.52756;
    else if (b.SampleType()==2018) w_year = 59.67377;
    return w_year*xs*1000.0/nents;
  });

  const NamedFunc weight_noscale("weight_noscale",[](const Baby &b) -> NamedFunc::ScalarType{
    //extra factor for weighting is xs/nevts_effective where effective events are calculated including negative weights
    //0.007519850838359999 GluGluH->ZG->llG xs (pb)
    //967054
    //double xs = 7.7760403e-09; //GluGluHToZG
    double xs = 0.00751985;
    double nents = 400000.0;
    if (b.FirstFileName().find("ZGToLLG") != std::string::npos) {
      xs = 117.864;
      //nents = 2.8095231e+09;
      nents = 1240092;
    }
    if (b.FirstFileName().find("DYJetsToLL") != std::string::npos) {
      xs = 6077.22;
      //nents = 3.2011619e+11;
      //nents = 3535863;
      nents = 7680238;
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
      if (abs_pdgid==11 || abs_pdgid == 13 || abs_pdgid == 15) {
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==23) {
            return abs_pdgid;
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
    std::vector<float> hlt_mu_pt;
    for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
      if (b.TrigObj_id()->at(itrig) == 13) {
        n_pf_mu++;
        if ((b.TrigObj_filterBits()->at(itrig) & 0x1)==1) { //TkIsoVVL
          hlt_mu_pt.push_back(b.TrigObj_pt()->at(itrig));
        }
      }
    }
    if (n_pf_mu < 2) return 1; //muon not reconstructed at HLT
    if (hlt_mu_pt.size() < 2) return 2; //muon fails Iso
    std::sort(hlt_mu_pt.begin(), hlt_mu_pt.end());
    if (hlt_mu_pt[hlt_mu_pt.size()-1] < 23.0) return 3; //top leg failed pt cut
    if (hlt_mu_pt[hlt_mu_pt.size()-2] < 12.0) return 4; //bottom leg failed pt cut
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

  const NamedFunc mu_sig("mu_sig",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> mu_sig_;
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      bool is_sig = true;
      if (b.Muon_pt()->at(imu) < 10) is_sig = false; //was 5
      if (abs(b.Muon_eta()->at(imu)) > 2.4) is_sig = false;
      if (abs(b.Muon_dz()->at(imu))>1.0)  is_sig = false;
      if (abs(b.Muon_dxy()->at(imu))>0.5) is_sig = false; 
      if (!((b.Muon_looseId()->at(imu) || (b.Muon_pt()->at(imu) > 200 && b.Muon_highPtId()->at(imu))) && 
          b.Muon_pfRelIso03_all()->at(imu) < 0.35 &&
          b.Muon_sip3d()->at(imu) < 4)) is_sig = false;
      if (is_sig)
        mu_sig_.push_back(1.0);
      else
        mu_sig_.push_back(0.0);
    }
    return mu_sig_;
  });

  const NamedFunc offline_nmu("offline_nmu",[mu_sig](const Baby &b) -> NamedFunc::ScalarType{
    int nmu_sig = 0;
    std::vector<double> mu_sig_ = mu_sig.GetVector(b);
    for (double imu_sig_ : mu_sig_) {
      if (imu_sig_ > 0.5)
        nmu_sig++;
    }
    return nmu_sig;
  });

  const NamedFunc el_sig("el_sig",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> el_sig_;
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      bool is_sig = true;
      if (b.Electron_pt()->at(iel) < 15) is_sig = false; //was 7
      float etasc = b.Electron_deltaEtaSC()->at(iel) + b.Electron_eta()->at(iel);
      if (abs(etasc) > 2.5) is_sig = false;
      if (abs(b.Electron_dz()->at(iel))>1.0) is_sig = false;
      if (abs(b.Electron_dxy()->at(iel))>0.5) is_sig = false; 
      double wp[2][3] = {{-0.145237, -0.0315746, -0.032173},
                         { 0.604775,  0.628743,   0.896462}};
      int ipt(1), ieta(2);
      if(b.Electron_pt()->at(iel)>10) ipt = 0;
      if(fabs(etasc) < 0.8) ieta = 0;
      else if(fabs(etasc) < 1.479) ieta = 1;
      double mva = b.Electron_mvaFall17V2Iso()->at(iel);
      if (mva <= wp[ipt][ieta])
        is_sig = false;
      if (is_sig)
        el_sig_.push_back(1.0);
      else
        el_sig_.push_back(0.0);
    }
    return el_sig_;
  });

  const NamedFunc Electron_hltIndex("Electron_hltIndex",[el_sig](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> el_matchhlt_;
    std::vector<double> el_sig_ = el_sig.GetVector(b);
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      double matchhlt = -1.0;
      if (el_sig_[iel] > 0.5) {
        //try to match to HLT object
        for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
          if (b.TrigObj_id()->at(itrig) == 11) {
            //if ((b.TrigObj_filterBits()->at(itrig) & 0x1)==1) { //CaloIdL_TrackIdL_IsoVL
            if (deltaR(b.Electron_eta()->at(iel), b.Electron_phi()->at(iel), b.TrigObj_eta()->at(itrig), b.TrigObj_phi()->at(itrig)) < 0.3) {
              matchhlt = static_cast<double>(itrig);
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

  const NamedFunc Muon_hltIndex("Muon_hltIndex",[mu_sig](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> mu_matchhlt_;
    std::vector<double> mu_sig_ = mu_sig.GetVector(b);
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      double matchhlt = -1.0;
      if (mu_sig_[imu] > 0.5) {
        //try to match to HLT object
        for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
          if (b.TrigObj_id()->at(itrig) == 13) {
            //if ((b.TrigObj_filterBits()->at(itrig) & 0x2)==1) { //Iso, TrkIsoVVL is broken at HLT
            if (deltaR(b.Muon_eta()->at(imu), b.Muon_phi()->at(imu), b.TrigObj_eta()->at(itrig), b.TrigObj_phi()->at(itrig)) < 0.3) {
              matchhlt = static_cast<double>(itrig);
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

  const NamedFunc el_gap("el_gap",[el_sig](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> el_sig_ = el_sig.GetVector(b);
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (el_sig_[iel] > 0.5) {
        float etasc = abs(b.Electron_deltaEtaSC()->at(iel) + b.Electron_eta()->at(iel));
        float eta = abs(b.Electron_eta()->at(iel));
        if (eta > 1.44 && eta < 1.56) return true;
        if (etasc > 1.44 && etasc < 1.56) return true;
      }
    }
    return false;
  });

  const NamedFunc offline_nel("offline_nel",[el_sig](const Baby &b) -> NamedFunc::ScalarType{
    int nel_sig = 0;
    std::vector<double> el_sig_ = el_sig.GetVector(b);
    for (double iel_sig_ : el_sig_) {
      if (iel_sig_ > 0.5)
        nel_sig++;
    }
    return nel_sig;
  });

  const NamedFunc Lead_Electron_pt("Lead_Electron_pt",[el_sig](const Baby &b) -> NamedFunc::ScalarType{
    double el_pt = 0;
    std::vector<double> el_sig_ = el_sig.GetVector(b);
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (el_sig_[iel]) {
        if (b.Electron_pt()->at(iel) > el_pt)
          el_pt = b.Electron_pt()->at(iel);
      }
    }
    return el_pt;
  });

  const NamedFunc Sublead_Electron_pt("Sublead_Electron_pt",[el_sig](const Baby &b) -> NamedFunc::ScalarType{
    double sublead_el_pt = 0;
    double lead_el_pt = 0;
    std::vector<double> el_sig_ = el_sig.GetVector(b);
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (el_sig_[iel]) {
        if (b.Electron_pt()->at(iel) > lead_el_pt) {
          sublead_el_pt = lead_el_pt;
          lead_el_pt = b.Electron_pt()->at(iel);
        }
        else if (b.Electron_pt()->at(iel) > sublead_el_pt)
          sublead_el_pt = b.Electron_pt()->at(iel);
      }
    }
    return sublead_el_pt;
  });

  const NamedFunc Lead_Muon_pt("Lead_Muon_pt",[mu_sig](const Baby &b) -> NamedFunc::ScalarType{
    double mu_pt = 0;
    std::vector<double> mu_sig_ = mu_sig.GetVector(b);
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      if (mu_sig_[imu]) {
        if (b.Muon_pt()->at(imu) > mu_pt)
          mu_pt = b.Muon_pt()->at(imu);
      }
    }
    return mu_pt;
  });

  const NamedFunc Sublead_Muon_pt("Sublead_Muon_pt",[mu_sig](const Baby &b) -> NamedFunc::ScalarType{
    double lead_mu_pt = 0;
    double sublead_mu_pt = 0;
    std::vector<double> mu_sig_ = mu_sig.GetVector(b);
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      if (mu_sig_[imu]) {
        if (b.Muon_pt()->at(imu) > lead_mu_pt) {
          sublead_mu_pt = lead_mu_pt;
          lead_mu_pt = b.Muon_pt()->at(imu);
        }
        else if (b.Muon_pt()->at(imu) > sublead_mu_pt)
          sublead_mu_pt = b.Muon_pt()->at(imu);
      }
    }
    return sublead_mu_pt;
  });

  const NamedFunc offline_nph("offline_nph",[](const Baby &b) -> NamedFunc::ScalarType{
    int nph_sig = 0;
    for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
      if (b.Photon_pt()->at(iph) < 15) continue;
      if (abs(b.Photon_eta()->at(iph)) > 2.5) continue;
      if (!(b.Photon_isScEtaEB()->at(iph) || b.Photon_isScEtaEE()->at(iph))) continue;
      // if (b.Photon_pfRelIso03_all()->at(iph)==PhotonRelIsoCut) continue; // no isolation cut in 2016...?
      if (b.Photon_mvaID()->at(iph) > 0.2 && 
          b.Photon_electronVeto()->at(iph))
      nph_sig++;
    }
    return nph_sig;
  });

  const NamedFunc TrigObj_filterBit2("TrigObj_filterBit2",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> filterbit2;
    for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
      filterbit2.push_back(static_cast<double>((b.TrigObj_filterBits()->at(itrig) & 0x02) >> 1));
    }
    return filterbit2;
  });

  const NamedFunc Photon_sig("Photon_sig",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> Photon_sig_;
    for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
      double sig = 1;
      if (b.Photon_pt()->at(iph) < 15) sig = 0;
      if (abs(b.Photon_eta()->at(iph)) > 2.5) sig = 0;
      if (!(b.Photon_isScEtaEB()->at(iph) || b.Photon_isScEtaEE()->at(iph))) sig = 0;
      // if (b.Photon_pfRelIso03_all()->at(iph)==PhotonRelIsoCut) continue; // no isolation cut in 2016...?
      if (!(b.Photon_mvaID()->at(iph) > 0.2 && 
          b.Photon_electronVeto()->at(iph))) sig = 0;
      Photon_sig_.push_back(sig);
    }
    return Photon_sig_;
  });

  const NamedFunc Lead_Photon_pt("Lead_Photon_pt",[Photon_sig](const Baby &b) -> NamedFunc::ScalarType{
    double ph_pt = 0;
    std::vector<double> ph_sig_ = Photon_sig.GetVector(b);
    for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
      if (ph_sig_[iph]) {
        if (b.Photon_pt()->at(iph) > ph_pt)
          ph_pt = b.Photon_pt()->at(iph);
      }
    }
    return ph_pt;
  });

  const NamedFunc Jet_sig("Jet_sig",[el_sig, mu_sig, Photon_sig](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> Jet_sig_;
    std::vector<double> el_sig_ = el_sig.GetVector(b);
    std::vector<double> mu_sig_ = mu_sig.GetVector(b);
    std::vector<double> ph_sig_ = Photon_sig.GetVector(b);
    for (unsigned ijet = 0; ijet < b.Jet_pt()->size(); ijet++) {
      double sig = 1;
      if (b.Jet_pt()->at(ijet) < 30) sig = 0;
      if (abs(b.Jet_eta()->at(ijet)) > 2.4) sig = 0;
      for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
        if (el_sig_[iel] > 0.5) {
          if (deltaR(b.Electron_eta()->at(iel), b.Electron_phi()->at(iel), b.Jet_eta()->at(ijet), b.Jet_phi()->at(ijet)) < 0.4)
            sig = 0;
        }
      }
      for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
        if (mu_sig_[imu] > 0.5) {
          if (deltaR(b.Muon_eta()->at(imu), b.Muon_phi()->at(imu), b.Jet_eta()->at(ijet), b.Jet_phi()->at(ijet)) < 0.4)
            sig = 0;
        }
      }
      for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
        if (ph_sig_[iph] > 0.5) {
          if (deltaR(b.Photon_eta()->at(iph), b.Photon_phi()->at(iph), b.Jet_eta()->at(ijet), b.Jet_phi()->at(ijet)) < 0.4)
            sig = 0;
        }
      }
      Jet_sig_.push_back(sig);
    }
    return Jet_sig_;
  });

  const NamedFunc Electron_IsoPhotonPt("Electron_IsoPhotonPt",[Photon_sig](const Baby &b) -> NamedFunc::VectorType{
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

  const NamedFunc Electron_IsoElectronPt("Electron_IsoElectronPt",[el_sig](const Baby &b) -> NamedFunc::VectorType{
    //-1 - no HLT object, 0 - fails all, 1 - pass CaloIdL_TrackIdL_IsoVL, 2 - pass WPTight
    std::vector<double> el_isoelpt;
    std::vector<double> el_sig_ = el_sig.GetVector(b);
    double isoelpt = 0;
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      isoelpt = 0;
      for (unsigned iel2 = 0; iel2 < b.Electron_pt()->size(); iel2++) {
        if (iel == iel2) continue;
        if (el_sig_[iel2]>0.5) {
          if (deltaR(b.Electron_eta()->at(iel), b.Electron_phi()->at(iel), b.Electron_eta()->at(iel2), b.Electron_phi()->at(iel2)) < 0.3) {
            isoelpt += b.Electron_pt()->at(iel2);
          }
        }
      }
      el_isoelpt.push_back(isoelpt);
    }
    return el_isoelpt;
  });

  const NamedFunc Electron_IsoJetPt("Electron_IsoJetPt",[Jet_sig](const Baby &b) -> NamedFunc::VectorType{
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

  const NamedFunc stitch_dy("stitch_dy",[Photon_sig](const Baby &b) -> NamedFunc::ScalarType{
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

  const NamedFunc ZCand_mass("ZCand_mass",[el_sig, mu_sig](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> el_sig_ = el_sig.GetVector(b);
    std::vector<double> mu_sig_ = mu_sig.GetVector(b);
    double m_ll = -999;
    for (unsigned iel = 1; iel < b.Electron_pt()->size(); iel++) {
      if (el_sig_[iel]) {
        for (unsigned iel2 = 0; iel2 < iel; iel2++) {
          if (el_sig_[iel2]) {
            ROOT::Math::PtEtaPhiMVector p1(b.Electron_pt()->at(iel),b.Electron_eta()->at(iel),b.Electron_phi()->at(iel),0.000511);
            ROOT::Math::PtEtaPhiMVector p2(b.Electron_pt()->at(iel2),b.Electron_eta()->at(iel2),b.Electron_phi()->at(iel2),0.000511);
            double this_mll = (p1+p2).M();
            if (fabs(this_mll-91.2)<fabs(m_ll-91.2))
              m_ll = this_mll;
          }
        }
      }
    }
    for (unsigned imu = 1; imu < b.Muon_pt()->size(); imu++) {
      if (mu_sig_[imu]) {
        for (unsigned imu2 = 0; imu2 < imu; imu2++) {
          if (mu_sig_[imu2]) {
            ROOT::Math::PtEtaPhiMVector p1(b.Muon_pt()->at(imu),b.Muon_eta()->at(imu),b.Muon_phi()->at(imu),0.106);
            ROOT::Math::PtEtaPhiMVector p2(b.Muon_pt()->at(imu2),b.Muon_eta()->at(imu2),b.Muon_phi()->at(imu2),0.106);
            double this_mll = (p1+p2).M();
            if (fabs(this_mll-91.2)<fabs(m_ll-91.2))
              m_ll = this_mll;
          }
        }
      }
    }
    return m_ll;
  });

  const NamedFunc HiggsCand_mass("HiggsCand_mass",[el_sig, mu_sig, Photon_sig](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> el_sig_ = el_sig.GetVector(b);
    std::vector<double> mu_sig_ = mu_sig.GetVector(b);
    std::vector<double> ph_sig_ = Photon_sig.GetVector(b);
    ROOT::Math::PtEtaPhiMVector ll_p;
    double m_ll = -999;
    double m_llg = -999;
    for (unsigned iel = 1; iel < b.Electron_pt()->size(); iel++) {
      if (el_sig_[iel]) {
        for (unsigned iel2 = 0; iel2 < iel; iel2++) {
          if (el_sig_[iel2]) {
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
      if (mu_sig_[imu]) {
        for (unsigned imu2 = 0; imu2 < imu; imu2++) {
          if (mu_sig_[imu2]) {
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
          if (fabs(this_mass-125.3)<fabs(m_llg-125.3))
            m_llg = this_mass;
        }
      }
    }
    return m_llg;
  });

  const NamedFunc HiggsCand_mass_cleaned("HiggsCand_mass_cleaned",[el_sig, mu_sig, Photon_sig](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> el_sig_ = el_sig.GetVector(b);
    std::vector<double> mu_sig_ = mu_sig.GetVector(b);
    std::vector<double> ph_sig_ = Photon_sig.GetVector(b);
    ROOT::Math::PtEtaPhiMVector ll_p;
    double m_ll = -999;
    double m_llg = -999;
    for (unsigned iel = 1; iel < b.Electron_pt()->size(); iel++) {
      if (el_sig_[iel]) {
        for (unsigned iel2 = 0; iel2 < iel; iel2++) {
          if (el_sig_[iel2]) {
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
      if (mu_sig_[imu]) {
        for (unsigned imu2 = 0; imu2 < imu; imu2++) {
          if (mu_sig_[imu2]) {
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

  const NamedFunc delta_r_cut("delta_r_cut",[el_sig, mu_sig, Photon_sig](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> el_sig_ = el_sig.GetVector(b);
    std::vector<double> mu_sig_ = mu_sig.GetVector(b);
    std::vector<double> ph_sig_ = Photon_sig.GetVector(b);
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (el_sig_[iel]) {
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
      if (mu_sig_[imu]) {
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

  const NamedFunc ElectronGenPhotonDeltaR("ElectronGenPhotonDeltaR",[el_sig](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> el_sig_ = el_sig.GetVector(b);
    double min_deltar = 999;
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (el_sig_[iel]) {
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

  const NamedFunc LeadElectron_PassReference("LeadElectron_PassReference",[el_sig, Electron_hltIndex](const Baby &b) -> NamedFunc::ScalarType{
    double leadel_passref = 0;
    std::vector<double> el_sig_ = el_sig.GetVector(b);
    std::vector<double> hlt_idx = Electron_hltIndex.GetVector(b);
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (el_sig_[iel] > 0.5) {
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

  const NamedFunc SubleadElectron_PassReference("SubleadElectron_PassReference",[el_sig, Electron_hltIndex](const Baby &b) -> NamedFunc::ScalarType{
    double leadel_passref = 0;
    int el_sig_idx = 0;
    std::vector<double> el_sig_ = el_sig.GetVector(b);
    std::vector<double> hlt_idx = Electron_hltIndex.GetVector(b);
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (el_sig_[iel] > 0.5) {
        if (el_sig_idx==1) {
          if (hlt_idx[iel]>=0) {
            if (((b.TrigObj_filterBits()->at(hlt_idx[iel]) & 0x2) == 0x2)
                && (b.TrigObj_pt()->at(hlt_idx[iel]) > 35)) {
              leadel_passref = 1;
            }
          }
          break;
        }
        el_sig_idx++;
      }
    }
    return leadel_passref;
  });

  const NamedFunc Electron_OtherPassReference("Electron_OtherPassReference",[el_sig, LeadElectron_PassReference, SubleadElectron_PassReference](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> otherpassref;
    double ot_passref = 0;
    double leadpassref = LeadElectron_PassReference.GetScalar(b);
    double subleadpassref = SubleadElectron_PassReference.GetScalar(b);
    int el_sig_idx = 0;
    std::vector<double> el_sig_ = el_sig.GetVector(b);
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      ot_passref = leadpassref;
      if (el_sig_[iel] > 0.5) {
        if (el_sig_idx==0) {
          ot_passref = subleadpassref;
        }
        el_sig_idx++;
      }
      otherpassref.push_back(ot_passref);
    }
    return otherpassref;
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

  const NamedFunc LeadMuon_PassReference("LeadMuon_PassReference",[mu_sig, Muon_hltIndex](const Baby &b) -> NamedFunc::ScalarType{
    double leadmu_passref = 0;
    std::vector<double> mu_sig_ = mu_sig.GetVector(b);
    std::vector<double> hlt_idx = Muon_hltIndex.GetVector(b);
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      if (mu_sig_[imu] > 0.5) {
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

  const NamedFunc SubleadMuon_PassReference("SubleadMuon_PassReference",[mu_sig, Muon_hltIndex](const Baby &b) -> NamedFunc::ScalarType{
    double leadmu_passref = 0;
    int mu_sig_idx = 0;
    std::vector<double> mu_sig_ = mu_sig.GetVector(b);
    std::vector<double> hlt_idx = Muon_hltIndex.GetVector(b);
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      if (mu_sig_[imu] > 0.5) {
        if (mu_sig_idx==1) {
          if (hlt_idx[imu]>=0) {
            if (((b.TrigObj_filterBits()->at(hlt_idx[imu]) & 0x2) == 0x2)
                && (b.TrigObj_pt()->at(hlt_idx[imu]) > 27)) {
              leadmu_passref = 1;
            }
          }
          break;
        }
        mu_sig_idx++;
      }
    }
    return leadmu_passref;
  });

  const NamedFunc Muon_OtherPassReference("Muon_OtherPassReference",[mu_sig, LeadMuon_PassReference, SubleadMuon_PassReference](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> otherpassref;
    double ot_passref = 0;
    double leadpassref = LeadMuon_PassReference.GetScalar(b);
    double subleadpassref = SubleadMuon_PassReference.GetScalar(b);
    int mu_sig_idx = 0;
    std::vector<double> mu_sig_ = mu_sig.GetVector(b);
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      ot_passref = leadpassref;
      if (mu_sig_[imu] > 0.5) {
        if (mu_sig_idx==0) {
          ot_passref = subleadpassref;
        }
        mu_sig_idx++;
      }
      otherpassref.push_back(ot_passref);
    }
    return otherpassref;
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

  NamedFunc baseline_selection = ((ZCand_mass>50)&&((Lead_Electron_pt>25)||(Lead_Muon_pt>20))&&(Lead_Photon_pt/HiggsCand_mass>0.14)&&delta_r_cut&&((HiggsCand_mass+ZCand_mass)>185));
  NamedFunc baseline_selection_nodr = ((ZCand_mass>50)&&((Lead_Electron_pt>25)||(Lead_Muon_pt>20))&&(Lead_Photon_pt/HiggsCand_mass>0.14)&&((HiggsCand_mass+ZCand_mass)>185));

  //------------------------------------------------------------------------------------
  //                                     make plots and pie charts
  //------------------------------------------------------------------------------------
  
  PlotMaker pm;

  bool plot_fail_dilep_trig = false;
  bool plot_fail_singlelep_trig = false;
  bool plot_more_fail_dilep = false;
  bool plot_cleanmask = false;
  bool plot_single_lep_trig_effects = false;
  bool make_rereco_plots = false;
  bool make_rereco_table_fail_reason = false;
  bool make_compare_preselection_table = false;
  bool plot_trigeffs = true;

  if (plot_fail_dilep_trig) {
    pm.Push<Hist1D>(Axis(4, -1.5, 2.5, Electron_hltId, "Electron HLT Status", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&(offline_nel>=2)&&(offline_nph>=1)&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__elhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(4, -1.5, 2.5, Muon_hltId, "Muon HLT Status", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&(offline_nmu>=2)&&(offline_nph>=1)&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__muhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(4, -1.5, 2.5, Electron_hltId, "Electron HLT Status", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_elhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(4, -1.5, 2.5, Muon_hltId, "Muon HLT Status", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(offline_nmu>=2)&&(offline_nph>=1)&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_muhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 100.0, "Electron_pt", "Electron p_{T} (no HLT match) [GeV]", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId<0.0&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecoel_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -2.5, 2.5, "Electron_eta", "Electron #eta (no HLT match)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId<0.0&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecoel_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -3.16, 3.16, "Electron_phi", "Electron #phi (no HLT match)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId<0.0&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecoel_phi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "Electron_eta", "Electron #eta (no HLT match)", {}),
        Axis(12, -3.16, 3.16, "Electron_phi", "Electron #phi (no HLT match)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId<0.0&&el_sig==1.0,
        signal_procs_noscale_2d, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecoel_etaphi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 2.0, "Electron_pfRelIso03_all", "Electron I_{rel,all} (no HLT match)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId<0.0&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecoel_iso")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 2.0, "Electron_pfRelIso03_chg", "Electron I_{rel,chg} (no HLT match)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId<0.0&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecoel_isochg")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 100.0, "Electron_pt", "Electron p_{T} (no HLT ID) [GeV]", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId==0.0&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidel_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -2.5, 2.5, "Electron_eta", "Electron #eta (no HLT ID)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId==0.0&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidel_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -3.16, 3.16, "Electron_phi", "Electron #phi (no HLT ID)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId==0.0&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidel_phi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "Electron_eta", "Electron #eta (no HLT ID)", {}),
        Axis(12, -3.16, 3.16, "Electron_phi", "Electron #phi (no HLT ID)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId==0.0&&el_sig==1.0,
        signal_procs_noscale_2d, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidel_etaphi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 2.0, "Electron_pfRelIso03_all", "Electron I_{rel,all} (no HLT ID)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId==0.0&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidel_iso")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 2.0, "Electron_pfRelIso03_chg", "Electron I_{rel,chg} (no HLT ID)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId==0.0&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidel_isochg")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 100.0, "Muon_pt", "Muon p_{T} (no HLT match) [GeV]", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(offline_nmu>=2)&&(offline_nph>=1)&&Muon_hltId<0.0&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -2.5, 2.5, "Muon_eta", "Muon #eta (no HLT match)", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(offline_nmu>=2)&&(offline_nph>=1)&&Muon_hltId<0.0&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -3.16, 3.16, "Muon_phi", "Muon #phi (no HLT match)", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(offline_nmu>=2)&&(offline_nph>=1)&&Muon_hltId<0.0&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_phi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "Muon_eta", "Muon #eta (no HLT match)", {}),
        Axis(12, -3.16, 3.16, "Muon_phi", "Muon #phi (no HLT match)", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(offline_nmu>=2)&&(offline_nph>=1)&&Muon_hltId<0.0&&mu_sig==1.0,
        signal_procs_noscale_2d, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_etaphi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 2.0, "Muon_pfRelIso03_all", "Muon I_{rel,all} (no HLT match)", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(offline_nmu>=2)&&(offline_nph>=1)&&Muon_hltId<0.0&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_iso")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 2.0, "Muon_pfRelIso03_chg", "Muon I_{rel,chg} (no HLT match)", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(offline_nmu>=2)&&(offline_nph>=1)&&Muon_hltId<0.0&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_isochg")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 100.0, "Muon_pt", "Muon p_{T} (no HLT ID) [GeV]", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(offline_nmu>=2)&&(offline_nph>=1)&&Muon_hltId==0.0&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidmu_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -2.5, 2.5, "Muon_eta", "Muon #eta (no HLT ID)", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(offline_nmu>=2)&&(offline_nph>=1)&&Muon_hltId==0.0&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidmu_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -3.16, 3.16, "Muon_phi", "Muon #phi (no HLT ID)", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(offline_nmu>=2)&&(offline_nph>=1)&&Muon_hltId==0.0&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidmu_phi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "Muon_eta", "Muon #eta (no HLT ID)", {}),
        Axis(12, -3.16, 3.16, "Muon_phi", "Muon #phi (no HLT ID)", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(offline_nmu>=2)&&(offline_nph>=1)&&Muon_hltId==0.0&&mu_sig==1.0,
        signal_procs_noscale_2d, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidmu_etaphi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 2.0, "Muon_pfRelIso03_all", "Muon I_{rel,all} (no HLT ID)", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(offline_nmu>=2)&&(offline_nph>=1)&&Muon_hltId==0.0&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidmu_iso")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 2.0, "Muon_pfRelIso03_chg", "Muon I_{rel,chg} (no HLT ID)", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(offline_nmu>=2)&&(offline_nph>=1)&&Muon_hltId==0.0&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidmu_isochg")
        .LuminosityTag(total_luminosity_string);
  }

  if (plot_more_fail_dilep) {
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "TrigObj_eta", "HLT Muon #eta", {}),
        Axis(12, -3.16, 3.16, "TrigObj_phi", "HLT Muon #phi", {}),
        "TrigObj_id==11" && (offline_nmu>=2)&&(offline_nph>=1),
        signal_procs_noscale_2d, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_hltmu_etaphi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(12, -3.16, 3.16, "TrigObj_phi", "HLT Muon #phi", {}),
        "TrigObj_id==11" && (offline_nmu>=2)&&(offline_nph>=1),
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
        (offline_nmu>=2)&&(offline_nph>=1)&&mu_sig,
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
        (offline_nmu>=2)&&(offline_nph>=1)&&mu_sig,
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
        (offline_nmu>=2)&&(offline_nph>=1)&&mu_sig,
        dy_procs_2016, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_offmu_etaphi_dy16")
        .LuminosityTag(total_luminosity_string);
    
    pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "Muon_nStations", "Muon nStations (No HLT Match)", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(offline_nmu>=2)&&(offline_nph>=1)&&Muon_hltId<0.0&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_nstations")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(21, -0.5, 20.5, "Muon_nTrackerLayers", "Muon nTrackerLayers (No HLT Match)", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(offline_nmu>=2)&&(offline_nph>=1)&&Muon_hltId<0.0&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_ntrackerlayers")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(2, -0.5, 1.5, "Muon_isGlobal", "Muon is global muon (No HLT Match)", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(offline_nmu>=2)&&(offline_nph>=1)&&Muon_hltId<0.0&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_isglobal")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "Muon_nStations", "Muon nStations (No HLT Match) (not 0.75 < |#eta| < 1.5)", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(offline_nmu>=2)&&(offline_nph>=1)&&Muon_hltId<0.0&&mu_sig==1.0&&"(Muon_eta>1.5||Muon_eta<-1.5||(Muon_eta>-0.75&&Muon_eta<0.75))",
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_etapeak_nstations")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(21, -0.5, 20.5, "Muon_nTrackerLayers", "Muon nTrackerLayers (No HLT Match) (not 0.75 < |#eta| < 1.5)", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(offline_nmu>=2)&&(offline_nph>=1)&&Muon_hltId<0.0&&mu_sig==1.0&&"(Muon_eta>1.5||Muon_eta<-1.5||(Muon_eta>-0.75&&Muon_eta<0.75))",
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_etapeak_ntrackerlayers")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(2, -0.5, 1.5, "Muon_isGlobal", "Muon is global muon (No HLT Match) (not 0.75 < |#eta| < 1.5)", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(offline_nmu>=2)&&(offline_nph>=1)&&Muon_hltId<0.0&&mu_sig==1.0&&"(Muon_eta>1.5||Muon_eta<-1.5||(Muon_eta>-0.75&&Muon_eta<0.75))",
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_etapeak_isglobal")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "Muon_nStations", "Muon nStations (No HLT Match) (0.75 < |#eta| < 1.5)", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(offline_nmu>=2)&&(offline_nph>=1)&&Muon_hltId<0.0&&mu_sig==1.0&&"((Muon_eta>-1.5&&Muon_eta<-0.75)||(Muon_eta>0.75&&Muon_eta<1.5))",
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_offetapeak_nstations")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(21, -0.5, 20.5, "Muon_nTrackerLayers", "Muon nTrackerLayers (No HLT Match) (0.75 < |#eta| < 1.5)", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(offline_nmu>=2)&&(offline_nph>=1)&&Muon_hltId<0.0&&mu_sig==1.0&&"((Muon_eta>-1.5&&Muon_eta<-0.75)||(Muon_eta>0.75&&Muon_eta<1.5))",
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_offetapeak_ntrackerlayers")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(2, -0.5, 1.5, "Muon_isGlobal", "Muon is global muon (No HLT Match) (0.75 < |#eta| < 1.5)", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&!hlt_mu_trigger&&(offline_nmu>=2)&&(offline_nph>=1)&&Muon_hltId<0.0&&mu_sig==1.0&&"((Muon_eta>-1.5&&Muon_eta<-0.75)||(Muon_eta>0.75&&Muon_eta<1.5))",
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecomu_offetapeak_isglobal")
        .LuminosityTag(total_luminosity_string);

    pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "Muon_nStations", "Offline Muon nStations", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&(offline_nmu>=2)&&(offline_nph>=1)&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_offmu_nstations")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(21, -0.5, 20.5, "Muon_nTrackerLayers", "Offline Muon nTrackerLayers", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&(offline_nmu>=2)&&(offline_nph>=1)&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_offmu_ntrackerlayers")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(2, -0.5, 1.5, "Muon_isGlobal", "Offline Muon is global muon", {}),
        l1_mu_trigger&&hlt_mu_trigger_plus&&(offline_nmu>=2)&&(offline_nph>=1)&&mu_sig==1.0,
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
        (offline_nel>=2)&&(offline_nph>=1)&&el_sig,
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
        (offline_nel>=2)&&(offline_nph>=1)&&el_sig,
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
        (offline_nel>=2)&&(offline_nph>=1)&&el_sig,
        dy_procs_2016, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_offel_etaphi_dy16")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(8, -0.5, 7.5, "Electron_lostHits", "Electron Lost Hits (No HLT Match)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId<0.0&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecoel_losthits")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, -0.1, 0.1, "Electron_deltaEtaSC", "Electron #Delta#eta_{SC} (No HLT Match)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId<0.0&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecoel_deltaeta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(15, 0.0, 0.03, "Electron_sieie", "Electron #sigma_{i#etai#eta} (No HLT Match)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId<0.0&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecoel_sieie")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 1.0, "Electron_r9", "Electron R_{9} (No HLT Match)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId<0.0&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecoel_r9")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 0.2, "Electron_hoe", "Electron H/E (No HLT Match)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId<0.0&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_norecoel_hoe")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(8, -0.5, 7.5, "Electron_lostHits", "Electron Lost Hits (No HLT ID)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId==0.0&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidel_losthits")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, -0.1, 0.1, "Electron_deltaEtaSC", "Electron #Delta#eta_{SC} (No HLT ID)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId==0.0&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidel_deltaeta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(15, 0.0, 0.03, "Electron_sieie", "Electron #sigma_{i#etai#eta} (No HLT ID)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId==0.0&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidel_sieie")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 1.0, "Electron_r9", "Electron R_{9} (No HLT ID)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId==0.0&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidel_r9")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 0.2, "Electron_hoe", "Electron H/E (No HLT ID)", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId==0.0&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidel_hoe")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(8, -0.5, 7.5, "Electron_lostHits", "Offline Electron Lost Hits", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_offel_losthits")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, -0.1, 0.1, "Electron_deltaEtaSC", "Offline Electron #Delta#eta_{SC}", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_offel_deltaeta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(15, 0.0, 0.03, "Electron_sieie", "Offline Electron #sigma_{i#etai#eta}", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_offel_sieie")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 1.0, "Electron_r9", "Offline Electron R_{9}", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_offel_r9")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 0.2, "Electron_hoe", "Offline Electron H/E", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_offel_hoe")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 50.0, Electron_IsoPhotonPt, "Sum photon p_{T} in electron isolation cone [GeV]", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId==0.0&&el_sig==1.0&&"Electron_pfRelIso03_all>0.4",
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidel_isophpt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 50.0, Electron_IsoElectronPt, "Sum electron p_{T} in electron isolation cone [GeV]", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId==0.0&&el_sig==1.0&&"Electron_pfRelIso03_all>0.4",
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidel_isoelpt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 50.0, Electron_IsoJetPt, "Sum jet p_{T} in electron isolation cone [GeV]", {}),
        l1_el_trigger&&hlt_el_trigger_plus&&!hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId==0.0&&el_sig==1.0&&"Electron_pfRelIso03_all>0.4",
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__dilepfail_noidel_isojetpt")
        .LuminosityTag(total_luminosity_string);
  }

  if (plot_fail_singlelep_trig) {
    pm.Push<Hist1D>(Axis(4, -1.5, 2.5, Electron_hltId, "Electron HLT Status", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(offline_nel>=2)&&(offline_nph>=1)&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_elhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(4, -1.5, 2.5, Muon_hltId, "Muon HLT Status", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(offline_nmu>=2)&&(offline_nph>=1)&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_muhlt")
        .LuminosityTag(total_luminosity_string);

    pm.Push<Hist1D>(Axis(20, 0.0, 100.0, "Electron_pt", "Electron p_{T} (no HLT match) [GeV]", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId<0.0&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecoel_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -2.5, 2.5, "Electron_eta", "Electron #eta (no HLT match)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId<0.0&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecoel_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -3.16, 3.16, "Electron_phi", "Electron #phi (no HLT match)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId<0.0&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecoel_phi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "Electron_eta", "Electron #eta (no HLT match)", {}),
        Axis(12, -3.16, 3.16, "Electron_phi", "Electron #phi (no HLT match)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId<0.0&&el_sig==1.0,
        signal_procs_noscale_2d, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecoel_etaphi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 2.0, "Electron_pfRelIso03_all", "Electron I_{rel,all} (no HLT match)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId<0.0&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecoel_iso")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(8, -0.5, 7.5, "Electron_lostHits", "Electron Lost Hits (no HLT match)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId<0.0&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecoel_losthits")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, -0.1, 0.1, "Electron_deltaEtaSC", "Electron #Delta#eta_{SC} (no HLT match)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId<0.0&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecoel_deltaeta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(15, 0.0, 0.03, "Electron_sieie", "Electron #sigma_{i#etai#eta} (no HLT match)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId<0.0&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecoel_sieie")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 1.0, "Electron_r9", "Electron R_{9} (no HLT match)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId<0.0&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecoel_r9")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 0.2, "Electron_hoe", "Electron H/E (no HLT match)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(offline_nel>=2)&&(offline_nph>=1)&&Electron_hltId<0.0&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecoel_hoe")
        .LuminosityTag(total_luminosity_string);

    pm.Push<Hist1D>(Axis(20, 0.0, 100.0, "Electron_pt", "Electron p_{T} (no HLT ID) [GeV]", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(offline_nel>=2)&&(offline_nph>=1)&&(Electron_hltId==0.0||Electron_hltId==1.0)&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidel_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -2.5, 2.5, "Electron_eta", "Electron #eta (no HLT ID)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(offline_nel>=2)&&(offline_nph>=1)&&(Electron_hltId==0.0||Electron_hltId==1.0)&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidel_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -3.16, 3.16, "Electron_phi", "Electron #phi (no HLT ID)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(offline_nel>=2)&&(offline_nph>=1)&&(Electron_hltId==0.0||Electron_hltId==1.0)&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidel_phi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "Electron_eta", "Electron #eta (no HLT ID)", {}),
        Axis(12, -3.16, 3.16, "Electron_phi", "Electron #phi (no HLT ID)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(offline_nel>=2)&&(offline_nph>=1)&&(Electron_hltId==0.0||Electron_hltId==1.0)&&el_sig==1.0,
        signal_procs_noscale_2d, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidel_etaphi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 2.0, "Electron_pfRelIso03_all", "Electron I_{rel,all} (no HLT ID)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(offline_nel>=2)&&(offline_nph>=1)&&(Electron_hltId==0.0||Electron_hltId==1.0)&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidel_iso")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(8, -0.5, 7.5, "Electron_lostHits", "Electron Lost Hits (no HLT ID)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(offline_nel>=2)&&(offline_nph>=1)&&(Electron_hltId==0.0||Electron_hltId==1.0)&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidel_losthits")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, -0.1, 0.1, "Electron_deltaEtaSC", "Electron #Delta#eta_{SC} (no HLT ID)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(offline_nel>=2)&&(offline_nph>=1)&&(Electron_hltId==0.0||Electron_hltId==1.0)&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidel_deltaeta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(15, 0.0, 0.03, "Electron_sieie", "Electron #sigma_{i#etai#eta} (no HLT ID)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(offline_nel>=2)&&(offline_nph>=1)&&(Electron_hltId==0.0||Electron_hltId==1.0)&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidel_sieie")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 1.0, "Electron_r9", "Electron R_{9} (no HLT ID)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(offline_nel>=2)&&(offline_nph>=1)&&(Electron_hltId==0.0||Electron_hltId==1.0)&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidel_r9")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 0.2, "Electron_hoe", "Electron H/E (no HLT ID)", {}),
        l1_el_trigger&&!hlt_el_trigger_plus&&(offline_nel>=2)&&(offline_nph>=1)&&(Electron_hltId==0.0||Electron_hltId==1.0)&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidel_hoe")
        .LuminosityTag(total_luminosity_string);

    pm.Push<Hist1D>(Axis(8, -0.5, 7.5, "Electron_lostHits", "Offline Electron Lost Hits", {}),
        (offline_nel>=2)&&(offline_nph>=1)&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_offel_lostHits")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, -0.1, 0.1, "Electron_deltaEtaSC", "Offline Electron #Delta#eta_{SC}", {}),
        (offline_nel>=2)&&(offline_nph>=1)&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_offel_deltaeta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(15, 0.0, 0.03, "Electron_sieie", "Offline Electron #sigma_{i#etai#eta}", {}),
        (offline_nel>=2)&&(offline_nph>=1)&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_offel_sieie")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 1.0, "Electron_r9", "Offline Electron R_{9}", {}),
        (offline_nel>=2)&&(offline_nph>=1)&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_offel_r9")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 0.0, 0.2, "Electron_hoe", "Offline Electron H/E", {}),
        (offline_nel>=2)&&(offline_nph>=1)&&el_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_offel_hoe")
        .LuminosityTag(total_luminosity_string);

    pm.Push<Hist1D>(Axis(20, 0.0, 100.0, "Muon_pt", "Muon p_{T} (no HLT match) [GeV]", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(offline_nmu>=2)&&(offline_nph>=1)&&Muon_hltId<0.0&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecomu_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -2.5, 2.5, "Muon_eta", "Muon #eta (no HLT match)", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(offline_nmu>=2)&&(offline_nph>=1)&&Muon_hltId<0.0&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecomu_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -3.16, 3.16, "Muon_phi", "Muon #phi (no HLT match)", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(offline_nmu>=2)&&(offline_nph>=1)&&Muon_hltId<0.0&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecomu_phi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "Muon_eta", "Muon #eta (no HLT match)", {}),
        Axis(12, -3.16, 3.16, "Muon_phi", "Muon #phi (no HLT match)", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(offline_nmu>=2)&&(offline_nph>=1)&&Muon_hltId<0.0&&mu_sig==1.0,
        signal_procs_noscale_2d, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecomu_etaphi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 2.0, "Muon_pfRelIso03_all", "Muon I_{rel,all} (no HLT match)", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(offline_nmu>=2)&&(offline_nph>=1)&&Muon_hltId<0.0&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecomu_iso")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "Muon_nStations", "Muon nStations (no HLT match)", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(offline_nmu>=2)&&(offline_nph>=1)&&Muon_hltId<0.0&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecomu_nstations")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(21, -0.5, 20.5, "Muon_nTrackerLayers", "Muon nTrackerLayers (no HLT match)", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(offline_nmu>=2)&&(offline_nph>=1)&&Muon_hltId<0.0&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecomu_ntrackerlayers")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(2, -0.5, 1.5, "Muon_isGlobal", "Muon is global muon (no HLT match)", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(offline_nmu>=2)&&(offline_nph>=1)&&Muon_hltId<0.0&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_norecomu_isglobal")
        .LuminosityTag(total_luminosity_string);

    pm.Push<Hist1D>(Axis(20, 0.0, 100.0, "Muon_pt", "Muon p_{T} (no HLT ID) [GeV]", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(offline_nmu>=2)&&(offline_nph>=1)&&(Muon_hltId==0.0||Muon_hltId==1.0)&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidmu_pt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -2.5, 2.5, "Muon_eta", "Muon #eta (no HLT ID)", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(offline_nmu>=2)&&(offline_nph>=1)&&(Muon_hltId==0.0||Muon_hltId==1.0)&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidmu_eta")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, -3.16, 3.16, "Muon_phi", "Muon #phi (no HLT ID)", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(offline_nmu>=2)&&(offline_nph>=1)&&(Muon_hltId==0.0||Muon_hltId==1.0)&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidmu_phi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(12, -2.5, 2.5, "Muon_eta", "Muon #eta (no HLT ID)", {}),
        Axis(12, -3.16, 3.16, "Muon_phi", "Muon #phi (no HLT ID)", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(offline_nmu>=2)&&(offline_nph>=1)&&(Muon_hltId==0.0||Muon_hltId==1.0)&&mu_sig==1.0,
        signal_procs_noscale_2d, twodim_plotopts).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidmu_etaphi")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0.0, 2.0, "Muon_pfRelIso03_all", "Muon I_{rel,all} (no HLT ID)", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(offline_nmu>=2)&&(offline_nph>=1)&&(Muon_hltId==0.0||Muon_hltId==1.0)&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidmu_iso")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "Muon_nStations", "Muon nStations (no HLT ID)", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(offline_nmu>=2)&&(offline_nph>=1)&&(Muon_hltId==0.0||Muon_hltId==1.0)&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidmu_nstations")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(21, -0.5, 20.5, "Muon_nTrackerLayers", "Muon nTrackerLayers (no HLT ID)", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(offline_nmu>=2)&&(offline_nph>=1)&&(Muon_hltId==0.0||Muon_hltId==1.0)&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidmu_ntrackerlayers")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(2, -0.5, 1.5, "Muon_isGlobal", "Muon is global muon (no HLT ID)", {}),
        l1_mu_trigger&&!hlt_mu_trigger_plus&&(offline_nmu>=2)&&(offline_nph>=1)&&(Muon_hltId==0.0||Muon_hltId==1.0)&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_noidmu_isglobal")
        .LuminosityTag(total_luminosity_string);

    pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "Muon_nStations", "Offline Muon nStations", {}),
        (offline_nmu>=2)&&(offline_nph>=1)&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_offmu_nstations")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(21, -0.5, 20.5, "Muon_nTrackerLayers", "Offline Muon nTrackerLayers", {}),
        (offline_nmu>=2)&&(offline_nph>=1)&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_offmu_ntrackerlayers")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(2, -0.5, 1.5, "Muon_isGlobal", "Offline Muon is global muon", {}),
        (offline_nmu>=2)&&(offline_nph>=1)&&mu_sig==1.0,
        signal_procs_noscale, plt_lin).Weight(weight_noscale)
        .Tag("FixName:zghlt__singlelepfail_offmu_isglobal")
        .LuminosityTag(total_luminosity_string);
  }

  if (plot_cleanmask) {
    pm.Push<Hist1D>(Axis(16, 100.0, 180.0, HiggsCand_mass, "m_{ee#gamma} [GeV]", {}),
        stitch_dy&&hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&baseline_selection_nodr,
        full_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__clean_mllg_diel")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 100.0, 180.0, HiggsCand_mass, "m_{#mu#mu#gamma} [GeV]", {}),
        stitch_dy&&hlt_mu_trigger&&(offline_nmu>=2)&&(offline_nph>=1)&&baseline_selection_nodr,
        full_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__clean_mllg_dimu")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 100.0, 180.0, HiggsCand_mass_cleaned, "m_{ee#gamma} [GeV]", {}),
        stitch_dy&&hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&baseline_selection_nodr,
        full_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__clean_mllgclean_diel")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 100.0, 180.0, HiggsCand_mass_cleaned, "m_{#mu#mu#gamma} [GeV]", {}),
        stitch_dy&&hlt_mu_trigger&&(offline_nmu>=2)&&(offline_nph>=1)&&baseline_selection_nodr,
        full_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__clean_mllgclean_dimu")
        .LuminosityTag(total_luminosity_string);
  }

  if (plot_single_lep_trig_effects) {
    pm.Push<Hist1D>(Axis(16, 100.0, 180.0, HiggsCand_mass, "m_{ee#gamma} [GeV]", {}),
        stitch_dy&&hlt_el_trigger&&(offline_nel>=2)&&(offline_nph>=1)&&baseline_selection,
        full_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__mllg_diel")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 100.0, 180.0, HiggsCand_mass, "m_{#mu#mu#gamma} [GeV]", {}),
        stitch_dy&&hlt_mu_trigger&&(offline_nmu>=2)&&(offline_nph>=1)&&baseline_selection,
        full_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__mllg_dimu")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 100.0, 180.0, HiggsCand_mass, "m_{ee#gamma} [GeV]", {}),
        stitch_dy&&(hlt_el_trigger||hlt_single_el_trigger)&&(offline_nel>=2)&&(offline_nph>=1)&&baseline_selection,
        full_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__mllg_eldiel")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 100.0, 180.0, HiggsCand_mass, "m_{#mu#mu#gamma} [GeV]", {}),
        stitch_dy&&(hlt_mu_trigger||hlt_single_mu_trigger)&&(offline_nmu>=2)&&(offline_nph>=1)&&baseline_selection,
        full_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__mllg_mudimu")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 100.0, 180.0, HiggsCand_mass, "m_{ee#gamma} [GeV]", {}),
        stitch_dy&&(!hlt_el_trigger&&hlt_single_el_trigger)&&(offline_nel>=2)&&(offline_nph>=1)&&baseline_selection,
        full_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__mllg_onlyel")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 100.0, 180.0, HiggsCand_mass, "m_{#mu#mu#gamma} [GeV]", {}),
        stitch_dy&&(!hlt_mu_trigger&&hlt_single_mu_trigger)&&(offline_nmu>=2)&&(offline_nph>=1)&&baseline_selection,
        full_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__mllg_onlymu")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 100.0, 180.0, HiggsCand_mass, "m_{ee#gamma} [GeV]", {}),
        stitch_dy&&(hlt_el_trigger&&!hlt_single_el_trigger)&&(offline_nel>=2)&&(offline_nph>=1)&&baseline_selection,
        full_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__mllg_onlydiel")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 100.0, 180.0, HiggsCand_mass, "m_{#mu#mu#gamma} [GeV]", {}),
        stitch_dy&&(hlt_mu_trigger&&!hlt_single_mu_trigger)&&(offline_nmu>=2)&&(offline_nph>=1)&&baseline_selection,
        full_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__mllg_onlydimu")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 100.0, 180.0, HiggsCand_mass, "m_{ee#gamma} [GeV]", {}),
        stitch_dy&&(hlt_el_trigger&&hlt_single_el_trigger)&&(offline_nel>=2)&&(offline_nph>=1)&&baseline_selection,
        full_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__mllg_botheldiel")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(16, 100.0, 180.0, HiggsCand_mass, "m_{#mu#mu#gamma} [GeV]", {}),
        stitch_dy&&(hlt_mu_trigger&&hlt_single_mu_trigger)&&(offline_nmu>=2)&&(offline_nph>=1)&&baseline_selection,
        full_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__mllg_bothmudimu")
        .LuminosityTag(total_luminosity_string);
    //single trigger plots
    pm.Push<Hist1D>(Axis(45, 50.0, 140.0, HiggsCand_mass, "m_{ee#gamma} [GeV]", {}),
        (offline_nel>=2)&&(offline_nph>=1),
        signal_trigger_procs, plt_shapes).Weight(weight)
        .Tag("FixName:zghlt__mllg_signal_shape_trigs_el")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(15, 110.0, 140.0, HiggsCand_mass, "m_{#mu#mu#gamma} [GeV]", {}),
        (offline_nmu>=2)&&(offline_nph>=1),
        signal_trigger_procs, plt_shapes).Weight(weight)
        .Tag("FixName:zghlt__mllg_signal_shape_trigs_mu")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(15, 110.0, 140.0, HiggsCand_mass, "m_{ee#gamma} [GeV]", {}),
        (offline_nel>=2)&&(offline_nph>=1)&&baseline_selection,
        signal_trigger_procs, plt_shapes).Weight(weight)
        .Tag("FixName:zghlt__mllg_signal_shape_trigs_baseline_el")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(15, 110.0, 140.0, HiggsCand_mass, "m_{#mu#mu#gamma} [GeV]", {}),
        (offline_nmu>=2)&&(offline_nph>=1)&&baseline_selection,
        signal_trigger_procs, plt_shapes).Weight(weight)
        .Tag("FixName:zghlt__mllg_signal_shape_trigs_baseline_mu")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(20, 50.0, 150.0, HiggsCand_mass, "m_{ee#gamma} [GeV]", {}),
        Axis(20, 0.0, 4.0, ElectronGenPhotonDeltaR, "#Delta R_{e#gamma}^{min}", {}),
        stitch_dy&&(!hlt_el_trigger&&hlt_single_el_trigger)&&(offline_nel>=2)&&(offline_nph>=1),
        signal_procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__mllg_deltaregamma")
        .LuminosityTag(total_luminosity_string);
  }

  if (make_rereco_plots) {
    pm.Push<Hist1D>(Axis(20, 0.0, 100.0, hlt_max_el_pt, "HLT Electron p_{T}", {}),
        (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&!hlt_el_trigger_plus,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__hltelpt__failhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(27, 0.0, 90.0, "Electron_pt", "Electron p_{T}", {}),
        (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&hlt_el_trigger&&el_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__elpt__passhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(27, 0.0, 90.0, "Electron_pt", "Electron p_{T}", {}),
        (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&!hlt_el_trigger&&el_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__elpt__failhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, -2.5, 2.5, "Electron_eta", "Electron #eta", {}),
        (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&hlt_el_trigger&&el_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__eleta__passhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, -2.5, 2.5, "Electron_eta", "Electron #eta", {}),
        (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&!hlt_el_trigger&&el_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__eleta__failhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, -2.5, 2.5, "Electron_eta", "Electron #eta", {}),
        (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==1)&&el_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__eleta__failhltreco")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 2.0, "Electron_pfRelIso03_all", "Electron I_{rel}", {}),
        (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&hlt_el_trigger&&el_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__elreliso__passhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 2.0, "Electron_pfRelIso03_all", "Electron I_{rel}", {}),
        (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&!hlt_el_trigger&&el_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__elreliso__failhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 2.0, "Electron_pfRelIso03_all", "Electron I_{rel}", {}),
        (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==2)&&el_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__elreliso__failhltiso")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 2.0, "Electron_pfRelIso03_chg", "Electron I_{rel,chg}", {}),
        (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&hlt_el_trigger&&el_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__elrelisochg__passhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 2.0, "Electron_pfRelIso03_chg", "Electron I_{rel,chg}", {}),
        (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&!hlt_el_trigger&&el_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__elrelisochg__failhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 2.0, "Muon_pfRelIso03_all", "Muon I_{rel}", {}),
        (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&l1_mu_trigger&&hlt_mu_trigger&&mu_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__mureliso__passhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 2.0, "Muon_pfRelIso03_all", "Muon I_{rel}", {}),
        (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&l1_mu_trigger&&!hlt_mu_trigger&&mu_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__mureliso__failhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 2.0, "Muon_pfRelIso03_all", "Muon I_{rel}", {}),
        (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==2)&&mu_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__mureliso__failhltiso")
        .LuminosityTag(total_luminosity_string);


    pm.Push<Hist2D>(Axis(15, -2.5, 2.5, "Electron_eta", "Electron #eta", {}), Axis(10, -3.14, 3.14, "Electron_phi", "Electron #phi", {}),
        (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&hlt_el_trigger&&el_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__eletaphi__passhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(15, -2.5, 2.5, "Electron_eta", "Electron #eta", {}), Axis(10, -3.14, 3.14, "Electron_phi", "Electron #phi", {}),
        (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&!hlt_el_trigger&&el_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__eletaphi__failhlt")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(15, -2.5, 2.5, "Electron_eta", "Electron #eta", {}), Axis(10, -3.14, 3.14, "Electron_phi", "Electron #phi", {}),
        (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==1)&&el_sig,
        procs, plt_lin).Weight(weight)
        .Tag("FixName:zghlt__eletaphi__failhltreco")
        .LuminosityTag(total_luminosity_string);
  }
  
  if (make_rereco_table_fail_reason) {
    pm.Push<Table>("zgtrig_el_"+options.year_string, vector<TableRow>{
      //TableRow("L1 Triggers", 
      //    (z_decay_pdgid==11)&&l1_el_trigger,0,0,weight),
      //TableRow("HLT Triggers", 
      //    (z_decay_pdgid==11)&&l1_el_trigger&&hlt_el_trigger,0,0,weight),

      TableRow("$Z\\rightarrow e^{+}e^{-}$ decays", 
          (z_decay_pdgid==11),0,0,weight),
      TableRow("Offline electrons in acceptance", 
          (z_decay_pdgid==11)&&(offline_nel>=2),0,0,weight),
      TableRow("Offline photon in acceptance", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1),0,0,weight),
      TableRow("L1 Triggers", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger,0,0,weight),
      TableRow("HLT Triggers", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&hlt_el_trigger,0,0,weight),
      TableRow("HLT Triggers (incl el)", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&hlt_el_trigger_plus,0,0,weight),
      TableRow("Events that fail HLT - $0\\gamma 0e$", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&!hlt_el_trigger_plus&&(nel_id_hlt<1)&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)",0,0,weight),
      TableRow("Events that fail HLT - $0\\gamma 1e$", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&!hlt_el_trigger_plus&&(nel_id_hlt==1)&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)",0,0,weight),
      TableRow("Events that fail HLT - $0\\gamma \\geq 1e$", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&!hlt_el_trigger_plus&&(nel_id_hlt>=2)&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)",0,0,weight),
      TableRow("Events that fail HLT - $1\\gamma 0e$", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&!hlt_el_trigger_plus&&(nel_id_hlt<1)&&"(HLT_Photon25||HLT_Photon20_HoverELoose)",0,0,weight),
      TableRow("Events that fail HLT - $1\\gamma 1e$", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&!hlt_el_trigger_plus&&(nel_id_hlt==1)&&"(HLT_Photon25||HLT_Photon20_HoverELoose)",0,0,weight),
      TableRow("Events that fail HLT - $1\\gamma \\geq 1e$", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&!hlt_el_trigger_plus&&(nel_id_hlt>=2)&&"(HLT_Photon25||HLT_Photon20_HoverELoose)",0,0,weight),
      TableRow("HLT Triggers (incl Ele27Ph50)", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&(hlt_el_trigger_plus||HLT_Ele27_CaloIdL_TrackIdL_IsoVL_Photon50),0,0,weight),
      TableRow("HLT Triggers (incl Ele27Ph33)", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&(hlt_el_trigger_plus||HLT_Ele27_CaloIdL_TrackIdL_IsoVL_Photon33),0,0,weight),
      TableRow("HLT Triggers (incl Ele27Ph25)", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&(hlt_el_trigger_plus||HLT_Ele27_CaloIdL_TrackIdL_IsoVL_Photon25),0,0,weight),
      TableRow("HLT Triggers (incl Ele27Ph20)", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&(hlt_el_trigger_plus||HLT_Ele27_CaloIdL_TrackIdL_IsoVL_Photon20_HoverELoose),0,0,weight),
      TableRow("HLT Triggers (incl Ele22Ph50)", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&(hlt_el_trigger_plus||HLT_Ele22_CaloIdL_TrackIdL_IsoVL_Photon50),0,0,weight),
      TableRow("HLT Triggers (incl Ele22Ph33)", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&(hlt_el_trigger_plus||HLT_Ele22_CaloIdL_TrackIdL_IsoVL_Photon33),0,0,weight),
      TableRow("HLT Triggers (incl Ele22Ph25)", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&(hlt_el_trigger_plus||HLT_Ele22_CaloIdL_TrackIdL_IsoVL_Photon25),0,0,weight),
      TableRow("HLT Triggers (incl Ele22Ph20)", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&(hlt_el_trigger_plus||HLT_Ele22_CaloIdL_TrackIdL_IsoVL_Photon20_HoverELoose),0,0,weight),
      TableRow("HLT Triggers (incl Ele17Ele8Ph25)", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&(hlt_el_trigger_plus||HLT_Ele17_Ele8_CaloIdL_TrackIdL_IsoVL_Photon25),0,0,weight),
      TableRow("HLT Triggers (incl Ele17Ele8Ph20)", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&(hlt_el_trigger_plus||HLT_Ele17_Ele8_CaloIdL_TrackIdL_IsoVL_Photon20_HoverELoose),0,0,weight),
      //TableRow("HLT Triggers (incl el, el+ph)", 
      //    (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&hlt_el_trigger_plusplus,0,0,weight),
      //TableRow("Electrons out of HLT p_{T} acceptance", 
      //    (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&el_oohltpt,0,0,weight),
      //TableRow("Electrons nearly out of HLT p_{T} acceptance", 
      //    (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&el_noohltpt,0,0,weight),
      //TableRow("Expected events with electrons failing HLT", 
      //    (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&el_oohltpt,0,0,weight*hlt_fail_weights),
      //TableRow("HLT fail reason: unknown", 
      //    (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==0.0),0,0,weight),
      //TableRow("HLT fail reason: reconstruction", 
      //    (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==1),0,0,weight),
      //TableRow("HLT fail reason: reconstruction, e in gap", 
      //    (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==1)&&el_gap,0,0,weight),
      //TableRow("HLT fail reason: id/iso", 
      //    (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==2),0,0,weight),
      //TableRow("HLT fail reason: top leg", 
      //    (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==3),0,0,weight),
      //TableRow("HLT fail reason: bottom leg", 
      //    (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==4),0,0,weight),

      TableRow("$Z\\rightarrow \\mu^{+}\\mu^{-}$ decays", 
          (z_decay_pdgid==13),0,0,weight),
      TableRow("Offline muons in acceptance", 
          (z_decay_pdgid==13)&&(offline_nmu>=2),0,0,weight),
      TableRow("Offline photon in acceptance", 
          (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1),0,0,weight),
      TableRow("L1 Triggers", 
          (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&l1_mu_trigger,0,0,weight),
      TableRow("HLT Triggers", 
          (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&l1_mu_trigger&&hlt_mu_trigger,0,0,weight),
      TableRow("HLT Triggers (incl mu)", 
          (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&l1_mu_trigger&&hlt_mu_trigger_plus,0,0,weight),
      TableRow("Events that fail HLT - $0\\gamma 0 \\mu$", 
          (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&l1_mu_trigger&&!hlt_mu_trigger_plus&&nmu_id_hlt<1&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)",0,0,weight),
      TableRow("Events that fail HLT - $0\\gamma 1 \\mu$", 
          (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&l1_mu_trigger&&!hlt_mu_trigger_plus&&nmu_id_hlt==1&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)",0,0,weight),
      TableRow("Events that fail HLT - $0\\gamma 2 \\mu$", 
          (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&l1_mu_trigger&&!hlt_mu_trigger_plus&&nmu_id_hlt>=2&&"!(HLT_Photon25||HLT_Photon20_HoverELoose)",0,0,weight),
      TableRow("Events that fail HLT - $1\\gamma 0 \\mu$", 
          (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&l1_mu_trigger&&!hlt_mu_trigger_plus&&nmu_id_hlt<1&&"(HLT_Photon25||HLT_Photon20_HoverELoose)",0,0,weight),
      TableRow("Events that fail HLT - $1\\gamma 1 \\mu$", 
          (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&l1_mu_trigger&&!hlt_mu_trigger_plus&&nmu_id_hlt==1&&"(HLT_Photon25||HLT_Photon20_HoverELoose)",0,0,weight),
      TableRow("Events that fail HLT - $1\\gamma 2 \\mu$", 
          (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&l1_mu_trigger&&!hlt_mu_trigger_plus&&nmu_id_hlt>=2&&"(HLT_Photon25||HLT_Photon20_HoverELoose)",0,0,weight),
      TableRow("HLT Triggers (incl mu, mu+ph)", 
          (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&l1_mu_trigger&&hlt_mu_trigger_plusplus,0,0,weight),
      //TableRow("Muons out of HLT acceptance", 
      //    (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&l1_mu_trigger&&hlt_mu_trigger&&mu_oohltpt,0,0,weight),
      //TableRow("HLT fail reason: unknown", 
      //    (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==0.0),0,0,weight),
      //TableRow("HLT fail reason: reconstruction", 
      //    (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==1),0,0,weight),
      //TableRow("HLT fail reason: iso", 
      //    (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==2),0,0,weight),
      //TableRow("HLT fail reason: top leg", 
      //    (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==3),0,0,weight),
      //TableRow("HLT fail reason: bottom leg", 
      //    (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==4),0,0,weight),

      TableRow("$Z\\rightarrow e^{+}e^{-}$ decays", 
          (z_decay_pdgid==11),0,0,weight),
      TableRow("L1 Triggers", 
          (z_decay_pdgid==11)&&l1_el_trigger,0,0,weight),
      TableRow("HLT Triggers", 
          (z_decay_pdgid==11)&&l1_el_trigger&&hlt_el_trigger,0,0,weight),
      TableRow("HLT fail reason: unknown", 
          (z_decay_pdgid==11)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==0.0),0,0,weight),
      TableRow("HLT fail reason: reconstruction", 
          (z_decay_pdgid==11)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==1),0,0,weight),
      TableRow("HLT fail reason: reconstruction, e in gap", 
          (z_decay_pdgid==11)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==1)&&el_gap,0,0,weight),
      TableRow("HLT fail reason: id/iso", 
          (z_decay_pdgid==11)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==2),0,0,weight),
      TableRow("HLT fail reason: top leg", 
          (z_decay_pdgid==11)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==3),0,0,weight),
      TableRow("HLT fail reason: bottom leg", 
          (z_decay_pdgid==11)&&l1_el_trigger&&!hlt_el_trigger&&(el_hlt_fail_reason==4),0,0,weight),

      TableRow("$Z\\rightarrow \\mu^{+}\\mu^{-}$ decays", 
          (z_decay_pdgid==13),0,0,weight),
      TableRow("L1 Triggers", 
          (z_decay_pdgid==13)&&l1_mu_trigger,0,0,weight),
      TableRow("HLT Triggers", 
          (z_decay_pdgid==13)&&l1_mu_trigger&&hlt_mu_trigger,0,0,weight),
      TableRow("HLT fail reason: unknwon", 
          (z_decay_pdgid==13)&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==0.0),0,0,weight),
      TableRow("HLT fail reason: reconstruction", 
          (z_decay_pdgid==13)&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==1),0,0,weight),
      TableRow("HLT fail reason: iso", 
          (z_decay_pdgid==13)&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==2),0,0,weight),
      TableRow("HLT fail reason: top leg", 
          (z_decay_pdgid==13)&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==3),0,0,weight),
      TableRow("HLT fail reason: bottom leg", 
          (z_decay_pdgid==13)&&l1_mu_trigger&&!hlt_mu_trigger&&(mu_hlt_fail_reason==4),0,0,weight),

      //TableRow("$Z\\rightarrow e^{+}e^{-}$ decays, not in l1 acceptance", 
      //    (z_decay_pdgid==11)&&((nel_ooeta>=2)||(nel_ool1pt)),0,0,weight),
      //TableRow("$Z\\rightarrow e^{+}e^{-}$ decays, not eta acceptance", 
      //    (z_decay_pdgid==11)&&((nel_ooeta>=2)),0,0,weight),
      //TableRow("$Z\\rightarrow e^{+}e^{-}$ decays, not l1 pt acceptance", 
      //    (z_decay_pdgid==11)&&((nel_ool1pt)),0,0,weight),
      //TableRow("$Z\\rightarrow e^{+}e^{-}$ decays, offline+L1 but not in hlt pt acceptance", 
      //    (z_decay_pdgid==11)&&offline_selection&&l1_el_trigger&&nel_oohltpt,0,0,weight),
      //TableRow("HLT Triggers", 
      //    (z_decay_pdgid==11)&&l1_el_trigger&&hlt_el_trigger,0,0,weight),
    },procs,false,true,false,false,false,true).LuminosityTag(total_luminosity_string).Precision(2);
  }
  
  if (make_compare_preselection_table) {
    pm.Push<Table>("ul_trig_eff"+options.year_string, vector<TableRow>{
      TableRow("\\hline\n $Z\\rightarrow e^{+}e^{-}$ decays", 
          (z_decay_pdgid==11),0,0,weight_noscale),
      TableRow("Offline electrons in acceptance", 
          (z_decay_pdgid==11)&&(offline_nel>=2),0,0,weight_noscale),
      TableRow("Offline photon in acceptance", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1),0,0,weight_noscale),
      TableRow("L1 Triggers", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger,0,0,weight_noscale),
      TableRow("HLT Single or Dilepton Triggers", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&hlt_el_trigger_plus,0,0,weight_noscale),
      TableRow("HLT Dilepton Triggers", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&l1_el_trigger&&hlt_el_trigger,0,0,weight_noscale),
      TableRow("\\hline\n Baseline Selection", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&baseline_selection,0,0,weight_noscale),
      TableRow("L1 Triggers", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&baseline_selection&&l1_el_trigger,0,0,weight_noscale),
      TableRow("HLT Single or Dilepton Triggers", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&baseline_selection&&l1_el_trigger&&hlt_el_trigger_plus,0,0,weight_noscale),
      TableRow("HLT Dilepton Triggers", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&baseline_selection&&l1_el_trigger&&hlt_el_trigger,0,0,weight_noscale),
      TableRow("\\hline\n Lepton $p_\\text{T}$ cuts", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&(Lead_Electron_pt>23)&&(Sublead_Electron_pt>12),0,0,weight_noscale),
      TableRow("L1 Triggers", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&(Lead_Electron_pt>23)&&(Sublead_Electron_pt>12)&&l1_el_trigger,0,0,weight_noscale),
      TableRow("HLT Single or Dilepton Triggers", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&(Lead_Electron_pt>23)&&(Sublead_Electron_pt>12)&&l1_el_trigger&&hlt_el_trigger_plus,0,0,weight_noscale),
      TableRow("HLT Dilepton Triggers", 
          (z_decay_pdgid==11)&&(offline_nel>=2)&&(offline_nph>=1)&&(Lead_Electron_pt>23)&&(Sublead_Electron_pt>12)&&l1_el_trigger&&hlt_el_trigger,0,0,weight_noscale),
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
          (z_decay_pdgid==13)&&(offline_nmu>=2),0,0,weight_noscale),
      TableRow("Offline photon in acceptance", 
          (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1),0,0,weight_noscale),
      TableRow("L1 Triggers", 
          (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&l1_mu_trigger,0,0,weight_noscale),
      TableRow("HLT Single or Dilepton Triggers", 
          (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&l1_mu_trigger&&hlt_mu_trigger_plus,0,0,weight_noscale),
      TableRow("HLT Dilepton Triggers", 
          (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&l1_mu_trigger&&hlt_mu_trigger,0,0,weight_noscale),
      TableRow("\\hline\n Baseline Selection", 
          (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&baseline_selection,0,0,weight_noscale),
      TableRow("L1 Triggers", 
          (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&baseline_selection&&l1_mu_trigger,0,0,weight_noscale),
      TableRow("HLT Single or Dilepton Triggers", 
          (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&baseline_selection&&l1_mu_trigger&&hlt_mu_trigger_plus,0,0,weight_noscale),
      TableRow("HLT Dilepton Triggers", 
          (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&baseline_selection&&l1_mu_trigger&&hlt_mu_trigger,0,0,weight_noscale),
      TableRow("\\hline\n Lepton $p_\\text{T}$ cuts", 
          (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&(Lead_Muon_pt>17)&&(Sublead_Muon_pt>10),0,0,weight_noscale),
      TableRow("L1 Triggers", 
          (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&(Lead_Muon_pt>17)&&(Sublead_Muon_pt>10)&&l1_mu_trigger,0,0,weight_noscale),
      TableRow("HLT Single or Dilepton Triggers", 
          (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&(Lead_Muon_pt>17)&&(Sublead_Muon_pt>10)&&l1_mu_trigger&&hlt_mu_trigger_plus,0,0,weight_noscale),
      TableRow("HLT Dilepton Triggers", 
          (z_decay_pdgid==13)&&(offline_nmu>=2)&&(offline_nph>=1)&&(Lead_Muon_pt>17)&&(Sublead_Muon_pt>10)&&l1_mu_trigger&&hlt_mu_trigger,0,0,weight_noscale),
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
  
  if (plot_trigeffs) {
    for (unsigned ieta = 0; ieta < el_abseta_bins.size()-1; ieta++) {
      pm.Push<EfficiencyPlot>(Axis(el_pt_bins, "Electron_pt", "Offline Electron p_{T} [GeV]", {}),
          "HLT_Ele35_WPTight_Gsf"&&(offline_nel==2)&&(ZCand_mass>81&&ZCand_mass<101)&&Electron_OtherPassReference&&el_sig&&Electron_absEta>el_abseta_bins[ieta]&&Electron_absEta<el_abseta_bins[ieta+1],
          Electron_hltPt>23&&Electron_hltId>0.5,
          procs_onelep_data,true,plt_lin).Weight("1").Tag("FixName:zghlt_dataeff_dieleleg23_eta"+std::to_string(ieta)).YTitle("Ele23_CaloIdL_TrackIdL_IsoVL").LuminosityTag(total_luminosity_string);
      pm.Push<EfficiencyPlot>(Axis(el_pt_bins, "Electron_pt", "Offline Electron p_{T} [GeV]", {}),
          "HLT_Ele35_WPTight_Gsf"&&(offline_nel==2)&&(ZCand_mass>81&&ZCand_mass<101)&&Electron_OtherPassReference&&el_sig&&Electron_absEta>el_abseta_bins[ieta]&&Electron_absEta<el_abseta_bins[ieta+1],
          Electron_hltPt>12&&Electron_hltId>0.5,
          procs_onelep_data,true,plt_lin).Weight("1").Tag("FixName:zghlt_dataeff_dieleleg12_eta"+std::to_string(ieta)).YTitle("Ele12_CaloIdL_TrackIdL_IsoVL").LuminosityTag(total_luminosity_string);
    }
    for (unsigned ieta = 0; ieta < mu_abseta_bins.size()-1; ieta++) {
      pm.Push<EfficiencyPlot>(Axis(mu_pt_bins, "Muon_pt", "Offline Muon p_{T} [GeV]", {}),
          "HLT_IsoMu27"&&(offline_nmu==2)&&(ZCand_mass>81&&ZCand_mass<101)&&Muon_OtherPassReference&&mu_sig&&Muon_absEta>mu_abseta_bins[ieta]&&Muon_absEta<mu_abseta_bins[ieta+1],
          Muon_hltPt>17&&Muon_hltId>0.5,
          procs_onelep_data,true,plt_lin).Weight("1").Tag("FixName:zghlt_dataeff_dimuleg17_eta"+std::to_string(ieta)).YTitle("Mu17_TrkIsoVVVL").LuminosityTag(total_luminosity_string);
      pm.Push<EfficiencyPlot>(Axis(mu_pt_bins, "Muon_pt", "Offline Muon p_{T} [GeV]", {}),
          "HLT_IsoMu27"&&(offline_nmu==2)&&(ZCand_mass>81&&ZCand_mass<101)&&Muon_OtherPassReference&&mu_sig&&Muon_absEta>mu_abseta_bins[ieta]&&Muon_absEta<mu_abseta_bins[ieta+1],
          Muon_hltPt>8&&Muon_hltId>0.5,
          procs_onelep_data,true,plt_lin).Weight("1").Tag("FixName:zghlt_dataeff_dimuleg8_eta"+std::to_string(ieta)).YTitle("Mu8_TrkIsoVVVL").LuminosityTag(total_luminosity_string);
    }
  }
  
  pm.multithreaded_ = !options.single_thread;
  pm.min_print_ = true;
  pm.MakePlots(1.);

  //------------------------------------------------------------------------------------
  //                                     post processing
  //------------------------------------------------------------------------------------

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

