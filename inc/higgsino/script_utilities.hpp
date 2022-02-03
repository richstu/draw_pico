#ifndef H_SCRIPT_UTILITIES
#define H_SCRIPT_UTILITIES

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
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

namespace script_utilities { 

  //struct for holding options passed via command line arguments
  struct ArgStruct {
    bool single_thread;
    std::string year_string;
    std::string string_options;
    std::string tag;
    bool unblind;
  };

  //functions
  
  //returns the palette typically used for the Higgsino analysis
  Palette colors();

  //function to parse command line arguments and return an ArgStruct
  //generic blurb to add to beginning of scripts
  //Arguments
  // --single_thread (-s) to run single thread for debugging
  // --unblind (-u) to unblind
  // --year (-y) yearname to select data year (2016, 2017, 2018, run2)
  // --tag (-t) to add a tag to produced plots
  // --string_options (-o) to specify what to plot
  ArgStruct get_options(int argc, char* argv[], std::string default_string_options="");

  //returns a vector containing a single plotopt corresponding to a data title linear plot
  //(with a ratio plot if unblind==true)
  std::vector<PlotOpt> plot_lin(bool unblind=false);

  //returns a vector containing a single plotopt corresponding to a data title log plot
  //(with a ratio plot if unblind==true)
  std::vector<PlotOpt> plot_log(bool unblind=false);

  //returns a vector containing a single plotopt corresponding to a data title linear shapes plot
  std::vector<PlotOpt> plot_shapes();

  //returns a vector containing a single plotopt corresponding to a data title log shapes plot
  std::vector<PlotOpt> plot_log_shapes();

  //returns a collection of tags, useful for defining processes
  std::map<std::string, std::set<std::string>> mctags();

  //returns a vector of HH+MET processes for use with plotmaker
  std::vector<std::shared_ptr<Process>> getall_processes(std::set<int> years,
      std::vector<std::string> signal_models, std::string skim_name, bool unblind=false,
      std::vector<int> signal_colors={}, bool simplify_mc=false, bool no_mc=false);

  //named func that applies CN weights to higgsino models with mlsp!=0
  extern const NamedFunc mixed_model_weight;

  //nb variable as used in higgsino analysis
  extern const NamedFunc hig_nb;

  //constants
  const std::string mc_production_folder = "/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath_v3/";
  const std::string search_mc_skim_folder = "mc/merged_higmc_higloose/";
  const std::string met150_mc_skim_folder = "mc/skim_met150/";
  const std::string skim1l2j_mc_skim_folder = "mc/skim_1l2j/";
  const std::string ttbar_mc_skim_folder = "mc/merged_higmc_higlep1T/";
  const std::string zll_mc_skim_folder = "mc/merged_higmc_higlep2T/";
  const std::string qcd_mc_skim_folder = "mc/merged_higmc_higqcd/";
  const std::string mc_unskimmed_folder = "mc/unskimmed/";

  const std::string data_production_folder = "/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath_v3/";
  const std::string search_data_skim_folder = "data/merged_higdata_higloose/";
  const std::string met150_data_skim_folder = "data/skim_met150/";
  const std::string skim1l2j_data_skim_folder = "data/skim_1l2j/";
  const std::string ttbar_data_skim_folder = "data/merged_higdata_higlep1T/";
  const std::string zll_data_skim_folder = "data/merged_higdata_higlep2T/";
  const std::string qcd_data_skim_folder = "data/merged_higdata_higqcd/";
  const std::string data_unskimmed_folder = "data/raw_pico/";

  const std::string signal_production_folder = "/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath_v3/";
  const std::string search_signal_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higloose/";
  const std::string met150_signal_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/skim_met150/";
  const std::string skim1l2j_signal_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/skim_1l2j/";
  const std::string ttbar_signal_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higlep1T/";
  const std::string zll_signal_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higlep2T/";
  const std::string qcd_signal_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higqcd/";
  const std::string signal_unskimmed_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/unskimmed/";

  const std::string gluinofast_production_folder = "/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath_v3/";
  const std::string search_gluinofast_skim_folder = "SMS-T5qqqqZH_fastSimJmeCorrection/merged_higmc_higloose/";
  const std::string met150_gluinofast_skim_folder = "SMS-T5qqqqZH_fastSimJmeCorrection/skim_met150/";
  const std::string skim1l2j_gluinofast_skim_folder = "SMS-T5qqqqZH_fastSimJmeCorrection/skim_1l2j/";
  const std::string ttbar_gluinofast_skim_folder = "SMS-T5qqqqZH_fastSimJmeCorrection/merged_higmc_higlep1T/";
  const std::string zll_gluinofast_skim_folder = "SMS-T5qqqqZH_fastSimJmeCorrection/merged_higmc_higlep2T/";
  const std::string qcd_gluinofast_skim_folder = "SMS-T5qqqqZH_fastSimJmeCorrection/merged_higmc_higqcd/";
  const std::string gluinofast_unskimmed_folder = "SMS-T5qqqqZH_fastSimJmeCorrection/unskimmed/";

  const std::string gluinofull_production_folder = "/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath_v3/";
  const std::string search_gluinofull_skim_folder = "SMS-T5qqqqZH_FullSimJmeVariations/merged_higmc_higloose/";
  const std::string met150_gluinofull_skim_folder = "SMS-T5qqqqZH_FullSimJmeVariations/skim_met150/";
  const std::string skim1l2j_gluinofull_skim_folder = "SMS-T5qqqqZH_FullSimJmeVariations/skim_1l2j/";
  const std::string ttbar_gluinofull_skim_folder = "SMS-T5qqqqZH_FullSimJmeVariations/merged_higmc_higlep1T/";
  const std::string zll_gluinofull_skim_folder = "SMS-T5qqqqZH_FullSimJmeVariations/merged_higmc_higlep2T/";
  const std::string qcd_gluinofull_skim_folder = "SMS-T5qqqqZH_FullSimJmeVariations/merged_higmc_higqcd/";
  const std::string gluinofull_unskimmed_folder = "SMS-T5qqqqZH_FullSimJmeVariations/unskimmed/";

  extern const NamedFunc search_resolved;
  extern const NamedFunc ttbar_resolved;
  extern const NamedFunc zll_resolved;
  extern const NamedFunc qcd_resolved;
}

#endif
