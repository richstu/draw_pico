//collection of things commonly used in plotting scripts
//maybe should be merged with higutilities.cpp

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <getopt.h>
#include <unistd.h>

#include "TColor.h"
#include "TError.h"
#include "TVector2.h"

#include "core/baby.hpp"
#include "core/cross_sections.hpp"
#include "core/event_scan.hpp"
#include "core/functions.hpp"
#include "core/hist1d.hpp"
#include "core/named_func.hpp"
#include "core/palette.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/process.hpp"
#include "core/table.hpp"
#include "core/utilities.hpp"
#include "higgsino/apply_trigeffs2016.hpp"
#include "higgsino/apply_trigeffs2017.hpp"
#include "higgsino/apply_trigeffs2018.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"
#include "higgsino/script_utilities.hpp"

namespace script_utilities { 

  Palette colors() {
    return Palette("txt/colors.txt","default");
  }

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

  std::map<std::string, std::set<std::string>> mctags() {
    std::map<std::string, std::set<std::string>> mctags_; 
    // Set base tags
    mctags_["tt"]     = std::set<std::string>({"*TTJets_SingleLept*","TTJets_DiLept*",
                                    "*_TTZ*.root", "*_TTW*.root",
                                   "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root"});
    mctags_["single_t"] = std::set<std::string>({"*_ST_*.root"});
    mctags_["vjets"]   = std::set<std::string>({"*_ZJet*.root", "*_WJetsToLNu*.root"});
    mctags_["zjets"]   = std::set<std::string>({"*_ZJet*.root", "*DYJetsToLL*.root"});
    mctags_["wjets"]   = std::set<std::string>({"*_WJetsToLNu*.root"});
    mctags_["qcd"]     = std::set<std::string>({"*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                                     "*_QCD_HT700to1000_*", "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                                     "*_QCD_HT2000toInf_*"});
    mctags_["other"]   = std::set<std::string>({"*_WH*.root", "*_ZH_HToBB*.root",
                                       "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});
    mctags_["other_and_single_t"]   = std::set<std::string>({"*_ST_*.root", "*_WH*.root", "*_ZH_HToBB*.root",
                                       "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});
    // Combine all tags
    mctags_["all"] = std::set<std::string>({"*TTJets_SingleLept*",
                                 "*TTJets_DiLept*",
                                 "*_TTZ*.root", "*_TTW*.root",
                                 "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root", "*_ST_*.root",
                                 "*_WJetsToLNu*.root", "*_ZJet*.root",
                                 "*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                                 "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                                 "*_QCD_HT2000toInf_*",
                                 "*_WH*.root", "*_ZH_HToBB*.root",
                                 "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root", "*DYJetsToLL*.root"});
    return mctags_;
  }

  std::vector<std::shared_ptr<Process>> getall_processes(std::set<int> years,
      std::vector<std::string> signal_models, std::string skim_name, bool unblind,
      std::vector<int> signal_colors, bool simplify_mc, bool no_mc) {
    std::unordered_map<std::string, std::string> mc_folder_dict;
    mc_folder_dict["search"] = search_mc_skim_folder;
    mc_folder_dict["met150"] = met150_mc_skim_folder;
    mc_folder_dict["1l2j"] = skim1l2j_mc_skim_folder;
    mc_folder_dict["ttbar"] = ttbar_mc_skim_folder;
    mc_folder_dict["zll"] = zll_mc_skim_folder;
    mc_folder_dict["qcd"] = qcd_mc_skim_folder;
    mc_folder_dict["unskimmed"] = mc_unskimmed_folder;
    std::unordered_map<std::string, std::string> data_folder_dict;
    data_folder_dict["search"] = search_data_skim_folder;
    data_folder_dict["met150"] = met150_data_skim_folder;
    data_folder_dict["1l2j"] = skim1l2j_data_skim_folder;
    data_folder_dict["ttbar"] = ttbar_data_skim_folder;
    data_folder_dict["zll"] = zll_data_skim_folder;
    data_folder_dict["qcd"] = qcd_data_skim_folder;
    data_folder_dict["unskimmed"] = data_unskimmed_folder;
    std::unordered_map<std::string, std::string> signal_folder_dict;
    signal_folder_dict["search"] = search_signal_skim_folder;
    signal_folder_dict["met150"] = met150_signal_skim_folder;
    signal_folder_dict["1l2j"] = skim1l2j_signal_skim_folder;
    signal_folder_dict["ttbar"] = ttbar_signal_skim_folder;
    signal_folder_dict["zll"] = zll_signal_skim_folder;
    signal_folder_dict["qcd"] = qcd_signal_skim_folder;
    signal_folder_dict["unskimmed"] = signal_unskimmed_folder;
    std::unordered_map<std::string, std::string> gluinofast_folder_dict;
    gluinofast_folder_dict["search"] = search_gluinofast_skim_folder;
    gluinofast_folder_dict["met150"] = met150_gluinofast_skim_folder;
    gluinofast_folder_dict["1l2j"] = skim1l2j_gluinofast_skim_folder;
    gluinofast_folder_dict["ttbar"] = ttbar_gluinofast_skim_folder;
    gluinofast_folder_dict["zll"] = zll_gluinofast_skim_folder;
    gluinofast_folder_dict["qcd"] = qcd_gluinofast_skim_folder;
    gluinofast_folder_dict["unskimmed"] = gluinofast_unskimmed_folder;
    std::unordered_map<std::string, std::string> gluinofull_folder_dict;
    gluinofull_folder_dict["search"] = search_gluinofull_skim_folder;
    gluinofull_folder_dict["met150"] = met150_gluinofull_skim_folder;
    gluinofull_folder_dict["1l2j"] = skim1l2j_gluinofull_skim_folder;
    gluinofull_folder_dict["ttbar"] = ttbar_gluinofull_skim_folder;
    gluinofull_folder_dict["zll"] = zll_gluinofull_skim_folder;
    gluinofull_folder_dict["qcd"] = qcd_gluinofull_skim_folder;
    gluinofull_folder_dict["unskimmed"] = gluinofull_unskimmed_folder;
    //std::unordered_map<std::string, NamedFunc> data_trigger_dict;
    //data_trigger_dict["search"] = Higfuncs::met_trigger;
    //data_trigger_dict["met150"] = Higfuncs::met_trigger;
    //data_trigger_dict["ttbar"] = Higfuncs::met_trigger || Higfuncs::el_trigger || Higfuncs::mu_trigger;
    //data_trigger_dict["zll"] = Higfuncs::el_trigger || Higfuncs::mu_trigger;
    //data_trigger_dict["qcd"] = Higfuncs::met_trigger;
    //data_trigger_dict["unskimmed"] = "1";
    //TODO: check for unsupported skimnames

    std::vector<std::shared_ptr<Process>> procs;
    //parse and add signal MC
    if (signal_colors.size() == 0) {
      if (signal_models.size() < 10) 
        signal_colors = {kGreen+1, kRed, kBlue, kYellow, kCyan, kOrange, kViolet, kMagenta, kBlue+4};
      else {
        std::cout << "ERROR: too many signal models, please recode this function" << std::endl;
        return procs;
      }
    }
    if (signal_models.size() > signal_colors.size()) {
      std::cout << "ERROR: insufficient signal colors" << std::endl;
      return procs;
    }
    for (unsigned sig_idx(0); sig_idx < signal_models.size(); sig_idx++) {
      std::string signal_model_string = signal_models[sig_idx];
      std::size_t open_parentheses_location = signal_model_string.find('(',0);
      std::size_t comma_location = signal_model_string.find(',',open_parentheses_location+1);
      std::size_t close_parentheses_location = signal_model_string.find(')',comma_location+1);
      std::string model_name = signal_model_string.substr(0,open_parentheses_location);
      std::string prod_mass = signal_model_string.substr(open_parentheses_location+1,
          comma_location-open_parentheses_location-1);
      std::string lsp_mass = signal_model_string.substr(comma_location+1,
          close_parentheses_location-comma_location-1);
      std::string sig_production_folder;
      std::string sig_skim_folder;
      std::set<std::string> sig_grep_string;
      if (model_name == "TChiHH" || model_name == "TChiHH-G") {
        sig_production_folder = signal_production_folder;
        sig_skim_folder = signal_folder_dict[skim_name];
        sig_grep_string = {"*TChiHH_mChi-"+prod_mass+"_mLSP-"+lsp_mass+"*.root"};
      }
      else if (model_name == "T5HH") {
        sig_production_folder = gluinofull_production_folder;
        sig_skim_folder = gluinofull_folder_dict[skim_name];
        sig_grep_string = {"*mGluino-"+prod_mass+"_*_mLSP-"+lsp_mass+"*.root"};
        //TODO: exclude points outside the 2d mass plane
      }
      else if (model_name == "T5HHfast") {
        sig_production_folder = gluinofast_production_folder;
        sig_skim_folder = gluinofast_folder_dict[skim_name];
        sig_grep_string = {"*mGluino-"+prod_mass+"_*_mLSP-"+lsp_mass+"*.root"};
        //TODO: exclude points outside the 2d mass plane
      }
      else {
        std::cout << "ERROR: unknown signal model name" << std::endl;
        return procs;
      }
      procs.push_back(Process::MakeShared<Baby_pico>(signal_model_string, 
          Process::Type::signal, signal_colors[sig_idx], attach_folder(
          sig_production_folder, years, sig_skim_folder, sig_grep_string), "stitch"));
    }

    //now add background mc
    if (!no_mc) {
      Palette mc_colors("txt/colors.txt","default");
      std::map<std::string, std::set<std::string>> mc_tags = mctags();
      if (!simplify_mc) {
        procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", 
            Process::Type::background, mc_colors("tt_htau"), attach_folder(
            mc_production_folder, years, mc_folder_dict[skim_name], mc_tags["tt"]),
            "stitch&&ntrutauh>0"));
        procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", 
            Process::Type::background, mc_colors("tt_1l"), attach_folder(
            mc_production_folder, years, mc_folder_dict[skim_name], mc_tags["tt"]),
            "stitch&&ntrutauh==0"));
        procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", 
            Process::Type::background, kOrange+1, attach_folder(mc_production_folder, 
            years, mc_folder_dict[skim_name], mc_tags["zjets"]), "stitch"));
        procs.push_back(Process::MakeShared<Baby_pico>("W+jets", 
            Process::Type::background, kGreen+1, attach_folder(mc_production_folder, 
            years, mc_folder_dict[skim_name], mc_tags["wjets"]), "stitch"));
      procs.push_back(Process::MakeShared<Baby_pico>("Single t", 
          Process::Type::background, mc_colors("single_t"), attach_folder(
          mc_production_folder, years, mc_folder_dict[skim_name], mc_tags["single_t"]),
          "stitch"));
      procs.push_back(Process::MakeShared<Baby_pico>("Other", 
          Process::Type::background, kGray+2, attach_folder(mc_production_folder, 
          years, mc_folder_dict[skim_name], mc_tags["other"]), "stitch"));
      }
      else {
        procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", 
            Process::Type::background, mc_colors("tt_htau"), attach_folder(
            mc_production_folder, years, mc_folder_dict[skim_name], mc_tags["tt"]),
            "stitch"));
        procs.push_back(Process::MakeShared<Baby_pico>("V+jets", 
            Process::Type::background, kOrange+1, attach_folder(mc_production_folder, 
            years, mc_folder_dict[skim_name], mc_tags["vjets"]), "stitch"));
      procs.push_back(Process::MakeShared<Baby_pico>("Other", 
          Process::Type::background, kGray+2, attach_folder(mc_production_folder, 
          years, mc_folder_dict[skim_name], mc_tags["other_and_single_t"]), "stitch"));
      }
      procs.push_back(Process::MakeShared<Baby_pico>("QCD", 
          Process::Type::background, mc_colors("other"), attach_folder(
          mc_production_folder, years, mc_folder_dict[skim_name], mc_tags["qcd"]),
          "stitch")); 
    }

    //add data
    if (unblind) {
      NamedFunc triggers = Higfuncs::met_trigger;
      if (skim_name == "ttbar") {
        triggers = Higfuncs::met_trigger || Higfuncs::el_trigger || Higfuncs::mu_trigger;
      }
      else if (skim_name == "zll") {
        triggers = Higfuncs::el_trigger || Higfuncs::mu_trigger;
      }
      procs.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, 
          kBlack, attach_folder(data_production_folder, years, data_folder_dict[skim_name], 
          {"*.root"}),triggers));
    }
    return procs;
  }

  //returns normal weight for everything other than TChiHH(mlsp!=0) which it switches to CN cross
  //section
  const NamedFunc mixed_model_weight("mixed_model_weight",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type() / 1000 != 106) return Higfuncs::final_weight.GetScalar(b); //not TChiHH
    if (b.mlsp() <= 1) return Higfuncs::final_weight.GetScalar(b); //N1N2
    double xsec1d, xsec2d, xsec1d_unc, xsec2d_unc;
    xsec::higgsinoCrossSection(b.mprod(),xsec1d,xsec1d_unc);
    xsec::higgsino2DCrossSection(b.mprod(),xsec2d,xsec2d_unc);
    return Higfuncs::final_weight.GetScalar(b)/xsec1d*xsec2d;
  });

  const NamedFunc hig_nb("hig_nb",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.nbt() < 2) return b.nbt();
    if (b.nbm() < 3) return b.nbm();
    if (b.nbl() == 3) return 3;
    return 4;
  });

  const NamedFunc search_resolved = 
                         "met/mht<2 && met/met_calo<2&&weight<1.5&&"
                         "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&nbt>=2&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40";
  const NamedFunc ttbar_resolved =
                         "met/met_calo<5&&weight<1.5&&"
                         "nlep==1&&mt<=100&&njet>=4&&njet<=5&&nbt>=2&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                         "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))"
                         &&Higfuncs::lead_signal_lepton_pt>30;

  //leppt cut implicit in skim
  const NamedFunc zll_resolved =
                         "met/met_calo<5&&weight<1.5&&"
                         "nlep==2&&njet>=4&&njet<=5&&met<50&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40";
  const NamedFunc qcd_resolved =
                         "met/mht<2 && met/met_calo<2&&"
                         "low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40";

}
