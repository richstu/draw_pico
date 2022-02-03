//This script generates the plots seen in the selection section of the AN.
//
//Arguments
// --single_thread (-s) to run single thread for debugging
// --unblind (-u) to unblind (not done for AN plots)
// --year (-y) yearname to select data year (2016, 2017, 2018, run2)
// --tag (-t) to add a tag to produced plots
// --string_options (-o) to specify what to plot
//
//possible string options (default indicates it is an AN plot): 
// plot_variables (default) - generate plots for resolved variables section (7.2.1)
// plot_nminus1 (default) - generates n-1 plots for selection section and appendix
// make_cutflow (default) - generates resolved MC cutflow 
// make_pies (default) - generate pie charts for 3 nb categories x 2 drmax bins x 4 met bins
// plot_nminus1_full - generates more versions of n-1 plots
// plot_am - generates am plots with various nb, met, and drmax selections 

#include "core/test.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>
#include <getopt.h>

#include "TError.h"
#include "TColor.h"
#include "TVector2.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/event_scan.hpp"
#include "core/utilities.hpp"
#include "core/functions.hpp"
#include "higgsino/apply_trigeffs2016.hpp"
#include "higgsino/apply_trigeffs2017.hpp"
#include "higgsino/apply_trigeffs2018.hpp"
#include "core/cross_sections.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"
#include "higgsino/script_utilities.hpp"

using namespace std;

int main(int argc, char *argv[]){
  //------------------------------------------------------------------------------------
  //                                    initialization
  //------------------------------------------------------------------------------------
  gErrorIgnoreLevel = 6000;
  time_t begtime, endtime;
  time(&begtime);
  Palette colors("txt/colors.txt","default");
  script_utilities::ArgStruct options = script_utilities::get_options(
      argc, argv, "plot_variables,plot_nminus1,make_cutflow,make_pies");
  std::vector<PlotOpt> plt_lin = script_utilities::plot_lin(false);
  std::vector<PlotOpt> plt_log = script_utilities::plot_log(false);
  std::vector<PlotOpt> plt_shapes = script_utilities::plot_shapes();
  std::vector<PlotOpt> plt_log_shapes = script_utilities::plot_log_shapes();
  set<int> years;
  HigUtilities::parseYears(options.year_string, years);
  string total_luminosity_string = HigUtilities::getLuminosityString(options.year_string);

  std::vector<std::shared_ptr<Process>> search_signal_procs = 
      script_utilities::getall_processes(
      years, {"TChiHH(175,0)","TChiHH(500,0)","TChiHH(900,0)",
      "TChiHH(350,200)","TChiHH(450,100)"}, "met150", false, 
      {kGreen+1, kRed, kBlue, kYellow, kCyan}, false, true);
  std::vector<std::shared_ptr<Process>> search_procs = 
      script_utilities::getall_processes(
      years, {"TChiHH(175,0)","TChiHH(500,0)","TChiHH(900,0)",
      "TChiHH(350,200)","TChiHH(450,100)"}, "met150");
  std::vector<std::shared_ptr<Process>> cutflow_procs = 
      script_utilities::getall_processes(
      years, {"TChiHH(175,0)","TChiHH(500,0)","TChiHH(900,0)",
      "TChiHH(350,200)","TChiHH(450,100)"}, "met150", false, 
      {kGreen+1, kRed, kBlue, kYellow, kCyan}, true);
  std::vector<std::shared_ptr<Process>> search_procs_higloose = 
      script_utilities::getall_processes(
      years, {"TChiHH(175,0)","TChiHH(500,0)","TChiHH(900,0)",
      "TChiHH(350,200)","TChiHH(450,100)"}, "search");

  //------------------------------------------------------------------------------------
  //                                     named funcs
  //------------------------------------------------------------------------------------
  
  NamedFunc filters = Higfuncs::final_pass_filters;
  NamedFunc base_filters = Higfuncs::final_pass_filters;
  NamedFunc weight = Higfuncs::final_weight;
  NamedFunc sr_baseline = filters&&script_utilities::search_resolved;
  NamedFunc sr_baseline_nob = filters&&
                         "met/mht<2 && met/met_calo<2&&weight<1.5&&"
                         "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40";
  NamedFunc sr_baseline_nonum = filters&&
                         "met/mht<2 && met/met_calo<2&&weight<1.5&&"
                         "njet>=4&&!low_dphi_met&&met>150&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40";
  NamedFunc sr_baseline_noqcd = filters&&
                         "weight<1.5&&ntk==0&&nvlep==0&&met>150&&njet>=4&&njet<=5&&nbt>=2&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40";
  NamedFunc sr_baseline_nohig = filters&&
                         "weight<1.5&&met/mht<2 && met/met_calo<2&&weight<1.5&&"
                         "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&nbt>=2";
  NamedFunc sr_baseline_onlynum = filters&&
                         "weight<1.5&&ntk==0&&nvlep==0&&met>150&&njet>=4&&njet<=5&&nbt>=2";
  NamedFunc mixed_model_weight = script_utilities::mixed_model_weight;
  NamedFunc search_resolved = script_utilities::search_resolved;
  NamedFunc hig_nb = script_utilities::hig_nb;

  //------------------------------------------------------------------------------------
  //                                     make plots
  //------------------------------------------------------------------------------------
  
  PlotMaker pm;

  // 2b3b4b 
  //pm.Push<Table>("FixName:selection__search_pies__2b3b4b_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved, 0, 0, weight)}), search_procs, true, true, true);
  
  if (HigUtilities::is_in_string_options(options.string_options,"plot_variables")) {
    //btag possibilities plots
    pm.Push<Hist1D>(Axis(5, 0.5, 5.5, hig_nb, "b-tag Category (TTML)", {}),
        sr_baseline,
        search_procs, plt_log).Weight(mixed_model_weight)
        .Tag("FixName:kinematicvars__signalbackground_nb_ttml_log_"+options.year_string)
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(5, 0.5, 5.5, "nbm", "b-tag Category (MMMM)", {}),
        sr_baseline_nob && "nbm>=2&&nbm<=4",
        search_procs, plt_log).Weight(mixed_model_weight)
        .Tag("FixName:kinematicvars__signalbackground_nb_mmmm_log_"+options.year_string)
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(5, 0.5, 5.5, "nbl", "b-tag Category (TTLL)", {}),
        sr_baseline && "nbl<=4",
        search_procs, plt_log).Weight(mixed_model_weight)
        .Tag("FixName:kinematicvars__signalbackground_nb_ttll_log_"+options.year_string)
        .LuminosityTag(total_luminosity_string);
    //hig candidate variable signal only plots
    pm.Push<Hist1D>(Axis(25, 0, 250, "hig_cand_am[0]", "#LT m_{bb} #GT", {100,140}),
        Higfuncs::pass_filters && "njet>=4&&nbt>=2&&nbm>=3&&nbl>=4",
        search_signal_procs, plt_shapes).Weight(mixed_model_weight)
        .Tag("FixName:kinematicvars__signal_am_shapes_"+options.year_string)
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0, 150, "hig_cand_dm[0]", "#Delta m_{HH}", {40}),
        Higfuncs::pass_filters && "njet>=4&&nbt>=2&&nbm>=3&&nbl>=4",
        search_signal_procs, plt_shapes).Weight(mixed_model_weight)
        .Tag("FixName:kinematicvars__signal_dm_shapes_"+options.year_string)
        .LuminosityTag(total_luminosity_string);
    //two more drmax plots
    pm.Push<Hist1D>(Axis(40, 0, 4, "hig_cand_drmax[0]", "#Delta R_{max}", {2.2}),
        Higfuncs::pass_filters && "njet>=4&&nbt>=2&&nbm>=3&&nbl>=4",
        search_signal_procs, plt_shapes).Weight(mixed_model_weight)
        .Tag("FixName:kinematicvars__signal_drmax_shapes_"+options.year_string)
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 4, "hig_cand_drmax[0]", "#Delta R_{max}", {2.2}),
        sr_baseline_nohig && "weight<1.5&&njet>=4&&nbt>=2&&nbm>=3&&nbl>=4&&hig_cand_am[0]>=100&&hig_cand_am[0]<=140&&hig_cand_dm[0]<40",
        search_procs, plt_lin).Weight(mixed_model_weight)
        .Tag("FixName:kinematicvars__signalbackground_drmax_lin_"+options.year_string)
        .LuminosityTag(total_luminosity_string);
    //bin variables signal only plots
    pm.Push<Hist1D>(Axis(40, 150, 800., "met", "p_{T}^{miss} [GeV]", {200., 300., 400.}),
        filters&&script_utilities::search_resolved,
        search_signal_procs, plt_shapes).Weight(mixed_model_weight)
        .Tag("FixName:selection__signalcomp_met_shapes_"+options.year_string)
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 2.2, "hig_cand_drmax[0]", "#Delta R_{max}", {1.1}),
        filters&&script_utilities::search_resolved,
        search_signal_procs, plt_shapes).Weight(mixed_model_weight)
        .Tag("FixName:selection__signalcomp_higcanddrmax_shapes_"+options.year_string)
        .LuminosityTag(total_luminosity_string);
    //pm.Push<Hist1D>(Axis(14, 150, 850., "met", "p_{T}^{miss} [GeV]", {200., 300., 400.}),
    //  base_filters&&search_resolved,
    //  search_procs, plt_log).Weight(weight).Tag("FixName:selection__search_met").LuminosityTag(total_luminosity_string);
  }

  if (HigUtilities::is_in_string_options(options.string_options,"plot_nminus1")) {
     //nb: special log plot since no nb=4 req.
     pm.Push<Hist1D>(Axis(5, -0.5, 4.5, hig_nb, "N_{b}", {1.5,2.5,3.5}),
         sr_baseline_nob,
         search_procs, plt_log).Weight(mixed_model_weight)
         .Tag("FixName:selection__nminus1_nb_log_"+options.year_string)
         .LuminosityTag(total_luminosity_string);
     pm.Push<Hist1D>(Axis(40, 150, 800, "met", "MET [GeV]", {200, 300, 400}),
         sr_baseline && "nbm>=3&&nbl>=4",
         search_procs, plt_lin).Weight(mixed_model_weight)
         .Tag("FixName:selection__nminus1_met_nb4_lin_"+options.year_string)
         .LuminosityTag(total_luminosity_string);
     pm.Push<Hist1D>(Axis(8, 3.5, 11.5, "njet", "N_{jets}", {3.5, 5.5}),
         sr_baseline_nonum && "nbt>=2&&nbm>=3&&nbl>=4&&nvlep==0&&ntk==0",
         search_procs, plt_lin).Weight(mixed_model_weight)
         .Tag("FixName:selection__nminus1_njet_nb4_lin_"+options.year_string)
         .LuminosityTag(total_luminosity_string);
     pm.Push<Hist1D>(Axis(4, -0.5, 3.5, "nvlep", "N_{vlep}", {0.5}),
         sr_baseline_nonum && "nbt>=2&&nbm>=3&&nbl>=4&&ntk==0&&njet>=4&&njet<=5",
         search_procs, plt_lin).Weight(mixed_model_weight)
         .Tag("FixName:selection__nminus1_nvlep_nb4_lin_"+options.year_string)
         .LuminosityTag(total_luminosity_string);
     pm.Push<Hist1D>(Axis(4, -0.5, 3.5, "ntk", "N_{isoTk}", {0.5}),
         sr_baseline_nonum && "nbt>=2&&nbm>=3&&nbl>=4&&nvlep==0&&njet>=4&&njet<=5",
         search_procs, plt_lin).Weight(mixed_model_weight)
         .Tag("FixName:selection__nminus1_nisotk_nb4_lin_"+options.year_string)
         .LuminosityTag(total_luminosity_string);
     pm.Push<Hist1D>(Axis(40, 0, 3.2, "jet_met_dphi[0]", "#Delta#Phi_{1}", {0.5}),
         sr_baseline_noqcd && "nbm>=3&&nbl>=4&&jet_met_dphi[1]>0.5 && jet_met_dphi[2]>0.3 &&" 
         "jet_met_dphi[3]>0.3 && (met/met_calo)<2 && (met/mht)<2",
         search_procs, plt_lin).Weight(mixed_model_weight)
         .Tag("FixName:selection__nminus1_dphi1_nb4_lin_"+options.year_string)
         .LuminosityTag(total_luminosity_string);
     pm.Push<Hist1D>(Axis(40, 0, 3.2, "jet_met_dphi[1]", "#Delta#Phi_{2}", {0.5}),
         sr_baseline_noqcd && "nbm>=3&&nbl>=4&&jet_met_dphi[0]>0.5 && jet_met_dphi[2]>0.3 &&" 
         "jet_met_dphi[3]>0.3 && (met/met_calo)<2 && (met/mht)<2",
         search_procs, plt_lin).Weight(mixed_model_weight)
         .Tag("FixName:selection__nminus1_dphi2_nb4_lin_"+options.year_string)
         .LuminosityTag(total_luminosity_string);
     pm.Push<Hist1D>(Axis(40, 0, 3.2, "jet_met_dphi[2]", "#Delta#Phi_{3}", {0.3}),
         sr_baseline_noqcd && "nbm>=3&&nbl>=4&&jet_met_dphi[0]>0.5 && jet_met_dphi[1]>0.5 &&" 
         "jet_met_dphi[3]>0.3 && (met/met_calo)<2 && (met/mht)<2",
         search_procs, plt_lin).Weight(mixed_model_weight)
         .Tag("FixName:selection__nminus1_dphi3_nb4_lin_"+options.year_string)
         .LuminosityTag(total_luminosity_string);
     pm.Push<Hist1D>(Axis(40, 0, 3.2, "jet_met_dphi[3]", "#Delta#Phi_{4}", {0.3}),
         sr_baseline_noqcd && "nbm>=3&&nbl>=4&& jet_met_dphi[0]>0.5 && jet_met_dphi[1]>0.5 &&"
         "jet_met_dphi[2]>0.3 && (met/met_calo)<2 && (met/mht)<2",
         search_procs, plt_lin).Weight(mixed_model_weight)
         .Tag("FixName:selection__nminus1_dphi4_nb4_lin_"+options.year_string)
         .LuminosityTag(total_luminosity_string);
     pm.Push<Hist1D>(Axis(40, 0, 8, "(met/mht)", "p^{miss}_{T}/H^{miss}_{T}", {2}),
         sr_baseline_noqcd && "nbm>=3&&nbl>=4&&!low_dphi_met&&(met/met_calo)<2",
         search_procs, plt_lin).Weight(mixed_model_weight)
         .Tag("FixName:selection__nminus1_metmht_nb4_lin_"+options.year_string)
         .LuminosityTag(total_luminosity_string);
     pm.Push<Hist1D>(Axis(40, 0, 8, "(met/met_calo)", "p^{miss}_{T}/p^{miss}_{Tcalo}", {2}),
         sr_baseline_noqcd && "nbm>=3&&nbl>=4&&!low_dphi_met&&(met/mht)<2",
         search_procs, plt_lin).Weight(mixed_model_weight)
         .Tag("FixName:selection__nminus1_metmetcalo_nb4_lin_"+options.year_string)
         .LuminosityTag(total_luminosity_string);
     pm.Push<Hist1D>(Axis(40, 0, 120, "hig_cand_dm[0]", "#Delta m_{HH} [GeV]", {40}),
         sr_baseline_nohig && "nbm>=3&&nbl>=4&&hig_cand_am[0]<200&&hig_cand_drmax[0]<2.2",
         search_procs, plt_lin).Weight(mixed_model_weight)
         .Tag("FixName:selection__nminus1_higcanddm_nb4_lin_"+options.year_string)
         .LuminosityTag(total_luminosity_string);
     pm.Push<Hist1D>(Axis(40, 0, 4.0, "hig_cand_drmax[0]", "#Delta R_{max}", {1.1,2.2}),
         sr_baseline_nohig && "nbm>=3&&nbl>=4&&hig_cand_am[0]<200&&hig_cand_dm[0]<40",
         search_procs, plt_lin).Weight(mixed_model_weight)
         .Tag("FixName:selection__nminus1_higcanddrmax_nb4_lin_"+options.year_string)
         .LuminosityTag(total_luminosity_string);
     pm.Push<Hist1D>(Axis(40, 0, 200, "hig_cand_am[0]", "#LT m_{bb} #GT [GeV]", {100,140}),
         sr_baseline && "nbm>=3&&nbl>=4",
         search_procs, plt_lin).Weight(mixed_model_weight)
         .Tag("FixName:selection__nminus1_higcandam_nb4_lin_"+options.year_string)
         .LuminosityTag(total_luminosity_string);
  }

  if (HigUtilities::is_in_string_options(options.string_options,"plot_nminus1_full")) {
    //selection section plots
    for (int plot_type_idx = 0; plot_type_idx < 3; plot_type_idx++) {
      vector<PlotOpt> plt_type = plt_log;
      vector<PlotOpt> plt_type_signal = plt_log;
      string plt_type_string = "log";
      if (plot_type_idx == 1) {
        plt_type = plt_shapes;
        plt_type_signal = plt_shapes;
        plt_type_string = "shapes";
      }
      else if (plot_type_idx == 2) {
        plt_type = plt_lin;
        plt_type_signal = plt_lin;
        plt_type_string = "lin";
      }

      if (plot_type_idx == 1) {
        //only make shape plots right now
        //signal only plots
        pm.Push<Hist1D>(Axis(40, 150, 800, "met", "MET [GeV]", {200, 300, 400}),
          sr_baseline,
          search_signal_procs, plt_type_signal).Weight(mixed_model_weight).Tag("FixName:selection__signalcomp_met_"+plt_type_string+"_"+options.year_string).LuminosityTag(total_luminosity_string);
        pm.Push<Hist1D>(Axis(20, 0, 2.2, Higfuncs::jetid_hig_cand_drmax, "#Delta R_{max}", {1.1}),
          sr_baseline_nohig && (Higfuncs::jetid_hig_cand_dm<40) && (Higfuncs::jetid_hig_cand_am<200),
          search_signal_procs, plt_type_signal).Weight(mixed_model_weight).Tag("FixName:selection__signalcomp_higcanddrmax_"+plt_type_string+"_"+options.year_string).LuminosityTag(total_luminosity_string);
      }

      for (int require_nb4 = 0; require_nb4 < 2; require_nb4++) {
        NamedFunc nb_cut = "1";
        string nb_desc = "";
        if (require_nb4 == 0 && plot_type_idx == 2) {
            //don't bother plotting linear plots w/o nb=4; signal too small to be visible
            continue;
        }
        if (require_nb4 == 1) {
            nb_cut = "nbm>=3&&nbt>=4";
            nb_desc = "nb4_";
        }
        
        // n-1 plots
        if (require_nb4==0 && plot_type_idx == 0) {
          //don't plot nb if we are requiring nb=4, right now only make log plot
          pm.Push<Hist1D>(Axis(5, -0.5, 4.5, Higfuncs::jetid_nb, "N_{b}", {1.5}),
            sr_baseline_nob,
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminus1_nb_"+nb_desc+plt_type_string+"_"+options.year_string).LuminosityTag(total_luminosity_string);
        }
        if (plot_type_idx == 2 && require_nb4==1) {
          //right now only make linear nb=4 plots

          pm.Push<Hist1D>(Axis(40, 150, 800, "met", "MET [GeV]", {200, 300, 400}),
            sr_baseline && nb_cut,
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminus1_met_"+nb_desc+plt_type_string+"_"+options.year_string).LuminosityTag(total_luminosity_string);

          pm.Push<Hist1D>(Axis(8, 3.5, 11.5, Higfuncs::jetid_njet, "N_{jets}", {3.5, 5.5}),
            sr_baseline_nonum && nb_cut && "nvlep==0&&ntk==0" && (Higfuncs::jetid_nb>=2),
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminus1_njet_"+nb_desc+plt_type_string+"_"+options.year_string).LuminosityTag(total_luminosity_string);

          pm.Push<Hist1D>(Axis(4, -0.5, 3.5, "nvlep", "N_{vlep}", {0.5}),
            sr_baseline_nonum && nb_cut && "ntk==0" && (Higfuncs::jetid_njet>=4) && (Higfuncs::jetid_njet<=5) && (Higfuncs::jetid_nb>=2),
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminus1_nvlep_"+nb_desc+plt_type_string+"_"+options.year_string).LuminosityTag(total_luminosity_string);

          pm.Push<Hist1D>(Axis(4, -0.5, 3.5, "ntk", "N_{isoTk}", {0.5}),
            sr_baseline_nonum && nb_cut && "nvlep==0" && (Higfuncs::jetid_njet>=4) && (Higfuncs::jetid_njet<=5) && (Higfuncs::jetid_nb>=2),
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminus1_nisotk_"+nb_desc+plt_type_string+"_"+options.year_string).LuminosityTag(total_luminosity_string);

          //pm.Push<Hist1D>(Axis(40, 0, 3.2, "jet_met_dphi[0]", "#Delta#Phi_{1}", {0.5}),
          //  HigUtilities::pass_2016 && nb_cut &&
          //  "nvlep==0 && 4 <= njet && njet <= 5 && nbt>=2 && met>150 && ntk==0 && (met/met_calo)<2 && (met/mht)<2 && hig_cand_dm[0]<=40 && hig_cand_drmax[0]<=2.2 && hig_cand_am[0]<=200",
          //  search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminusphi_dphi1_"+nb_desc+plt_type_string+"_"+options.year_string).LuminosityTag(total_luminosity_string);
          //pm.Push<Hist1D>(Axis(40, 0, 3.2, "jet_met_dphi[1]", "#Delta#Phi_{2}", {0.5}),
          //  HigUtilities::pass_2016 && nb_cut &&
          //  "nvlep==0 && 4 <= njet && njet <= 5 && nbt>=2 && met>150 && ntk==0 && (met/met_calo)<2 && (met/mht)<2 && hig_cand_dm[0]<=40 && hig_cand_drmax[0]<=2.2 && hig_cand_am[0]<=200",
          //  search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminusphi_dphi2_"+nb_desc+plt_type_string+"_"+options.year_string).LuminosityTag(total_luminosity_string);
          //pm.Push<Hist1D>(Axis(40, 0, 3.2, "jet_met_dphi[2]", "#Delta#Phi_{3}", {0.3}),
          //  HigUtilities::pass_2016 && nb_cut &&
          //  "nvlep==0 && 4 <= njet && njet <= 5 && nbt>=2 && met>150 && ntk==0 && (met/met_calo)<2 && (met/mht)<2 && hig_cand_dm[0]<=40 && hig_cand_drmax[0]<=2.2 && hig_cand_am[0]<=200",
          //  search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminusphi_dphi3_"+nb_desc+plt_type_string+"_"+options.year_string).LuminosityTag(total_luminosity_string);
          //pm.Push<Hist1D>(Axis(40, 0, 3.2, "jet_met_dphi[3]", "#Delta#Phi_{4}", {0.3}),
          //  HigUtilities::pass_2016 && nb_cut &&
          //  "nvlep==0 && 4 <= njet && njet <= 5 && nbt>=2 && met>150 && ntk==0 && (met/met_calo)<2 && (met/mht)<2 && hig_cand_dm[0]<=40 && hig_cand_drmax[0]<=2.2 && hig_cand_am[0]<=200",
          //  search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminusphi_dphi4_"+nb_desc+plt_type_string+"_"+options.year_string).LuminosityTag(total_luminosity_string);

          pm.Push<Hist1D>(Axis(40, 0, 3.2, "jet_met_dphi[0]", "#Delta#Phi_{1}", {0.5}),
            sr_baseline_noqcd && nb_cut && "pass_jets && jet_met_dphi[1]>0.5 && jet_met_dphi[2]>0.3 && jet_met_dphi[3]>0.3 && (met/met_calo)<2 && (met/mht)<2",
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminus1_dphi1_"+nb_desc+plt_type_string+"_"+options.year_string).LuminosityTag(total_luminosity_string);

          pm.Push<Hist1D>(Axis(40, 0, 3.2, "jet_met_dphi[1]", "#Delta#Phi_{2}", {0.5}),
            sr_baseline_noqcd && nb_cut && "pass_jets && jet_met_dphi[0]>0.5 && jet_met_dphi[2]>0.3 && jet_met_dphi[3]>0.3 && (met/met_calo)<2 && (met/mht)<2",
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminus1_dphi2_"+nb_desc+plt_type_string+"_"+options.year_string).LuminosityTag(total_luminosity_string);

          pm.Push<Hist1D>(Axis(40, 0, 3.2, "jet_met_dphi[2]", "#Delta#Phi_{3}", {0.3}),
            sr_baseline_noqcd && nb_cut && "pass_jets && jet_met_dphi[0]>0.5 && jet_met_dphi[1]>0.5 && jet_met_dphi[3]>0.3 && (met/met_calo)<2 && (met/mht)<2",
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminus1_dphi3_"+nb_desc+plt_type_string+"_"+options.year_string).LuminosityTag(total_luminosity_string);

          pm.Push<Hist1D>(Axis(40, 0, 3.2, "jet_met_dphi[3]", "#Delta#Phi_{4}", {0.3}),
            sr_baseline_noqcd && nb_cut && "pass_jets && jet_met_dphi[0]>0.5 && jet_met_dphi[1]>0.5 && jet_met_dphi[2]>0.3 && (met/met_calo)<2 && (met/mht)<2",
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminus1_dphi4_"+nb_desc+plt_type_string+"_"+options.year_string).LuminosityTag(total_luminosity_string);

          pm.Push<Hist1D>(Axis(40, 0, 8, "(met/mht)", "p^{miss}_{T}/H^{miss}_{T}", {2}),
            sr_baseline_noqcd && nb_cut && !Higfuncs::jetid_low_dphi_met && "(met/met_calo)<2",
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminus1_metmht_"+nb_desc+plt_type_string+"_"+options.year_string).LuminosityTag(total_luminosity_string);

          pm.Push<Hist1D>(Axis(40, 0, 8, "(met/met_calo)", "p^{miss}_{T}/p^{miss}_{Tcalo}", {2}),
            sr_baseline_noqcd && nb_cut && !Higfuncs::jetid_low_dphi_met && "(met/mht)<2",
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminus1_metmetcalo_"+nb_desc+plt_type_string+"_"+options.year_string).LuminosityTag(total_luminosity_string);

          pm.Push<Hist1D>(Axis(40, 0, 120, Higfuncs::jetid_hig_cand_dm, "#Delta m_{HH} [GeV]", {40}),
            sr_baseline_nohig && nb_cut && (Higfuncs::jetid_hig_cand_am<200) && (Higfuncs::jetid_hig_cand_drmax<2.2),
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminus1_higcanddm_"+nb_desc+plt_type_string+"_"+options.year_string).LuminosityTag(total_luminosity_string);

          pm.Push<Hist1D>(Axis(40, 0, 4.0, Higfuncs::jetid_hig_cand_drmax, "#Delta R_{max}", {2.2}),
            sr_baseline_nohig && nb_cut && (Higfuncs::jetid_hig_cand_dm<40) && (Higfuncs::jetid_hig_cand_am<200),
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminus1_higcanddrmax_"+nb_desc+plt_type_string+"_"+options.year_string).LuminosityTag(total_luminosity_string);

          pm.Push<Hist1D>(Axis(40, 0, 200, Higfuncs::jetid_hig_cand_am, "#LT m_{bb} #GT [GeV]", {100,140}),
            sr_baseline && nb_cut,
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminus1_higcandam_"+nb_desc+plt_type_string+"_"+options.year_string).LuminosityTag(total_luminosity_string);

          //non N-1 variables
          pm.Push<Hist1D>(Axis(40, 0, 1200, "ht", "HT [GeV]", {}),
            sr_baseline && nb_cut,
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__baseline_ht_"+nb_desc+plt_type_string+"_"+options.year_string).LuminosityTag(total_luminosity_string);
        }
      }
    }
  }

  if (HigUtilities::is_in_string_options(options.string_options,"make_cutflow")) {
    pm.Push<Table>("sr_cutflow_"+options.year_string, vector<TableRow>{
      TableRow("Filters, $p_\\text{t}^\\text{miss}>150 \\text{ GeV}$", 
          Higfuncs::pass_filters && "met>150",0,0,mixed_model_weight),
      TableRow("$N_\\text{vl}=N_\\text{tk}=0$", 
          Higfuncs::pass_filters && "met>150 && nvlep==0 && ntk==0",0,0,mixed_model_weight),
      TableRow("$4 \\leq N_\\text{jet} \\leq 5$", 
          Higfuncs::pass_filters && "met>150 && nvlep==0 && ntk==0 && njet >= 4 && njet <= 5",
          0,0,mixed_model_weight),
      TableRow("$N_\\text{b}\\geq 2$", 
          sr_baseline_onlynum,0,0,mixed_model_weight),
      TableRow("Fake MET Cuts",        
          sr_baseline_nohig,0,0,mixed_model_weight),
      TableRow("$\\Delta m<40 \\text{ GeV}$, $\\langle m_\\text{bb} \\rangle<200 \\text{ GeV}$",        
          sr_baseline_nohig && "hig_cand_dm[0]<40&&hig_cand_am[0]<200",0,0,mixed_model_weight),
      TableRow("$\\Delta R_\\text{max}<2.2$",        
          sr_baseline,0,0,mixed_model_weight),
      TableRow("\\hline\n $100 \\leq \\langle m\\rangle < 140 \\text{ GeV}$",        
          sr_baseline && "hig_cand_am[0]>100&&hig_cand_am[0]<140",0,0,mixed_model_weight),
      TableRow("$N_\\text{b}\\geq 3$",        
          sr_baseline && "hig_cand_am[0]>100&&hig_cand_am[0]<140&&nbm>=3",0,0,mixed_model_weight),
      TableRow("$N_\\text{b}=4$",        
          sr_baseline && "hig_cand_am[0]>100&&hig_cand_am[0]<140&&nbm>=3&&nbl>=4",
          0,0,mixed_model_weight),
      TableRow("$p_\\text{T}^\\text{miss}>200 \\text{ GeV}$",        
          sr_baseline && "hig_cand_am[0]>100&&hig_cand_am[0]<140&&nbm>=3&&nbl>=4&&met>200",
          0,0,mixed_model_weight),
      TableRow("$p_\\text{T}^\\text{miss}>300 \\text{ GeV}$",        
          sr_baseline && "hig_cand_am[0]>100&&hig_cand_am[0]<140&&nbm>=3&&nbl>=4&&met>300",
          0,0,mixed_model_weight),
      TableRow("$p_\\text{T}^\\text{miss}>400 \\text{ GeV}$",        
          sr_baseline && "hig_cand_am[0]>100&&hig_cand_am[0]<140&&nbm>=3&&nbl>=4&&met>400",
          0,0,mixed_model_weight),
    //TableRow("$\\Delta R_\\text{max}>1.1$",        
    //  "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>450 && ntk==0 && !low_dphi_met && (met/met_calo)<2 && (met/mht)<2 && hig_cand_dm[0]<=40 && hig_cand_drmax[0]>=1.1 && hig_cand_drmax[0]<=2.2 && hig_cand_am[0]>100 && hig_cand_am[0]<=140 && nbm>=3 && nbl>=4",0,0,mixed_model_weight),
    },cutflow_procs,false,true,false,true,false,true).LuminosityTag(total_luminosity_string);
  }

  //pm.Push<Hist1D>(Axis(10,0,100,"hig_cand_dm[0]", "#Deltam [GeV]", {40.}),
  //  base_filters &&
  //  "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
  //  "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&"
  //  "((nbt>=2&&nbm>=3&&nbl>=4))",
  //  search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_hig_cand_dm").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters &&
  //  "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
  //  "hig_cand_drmax[0]<2.2&&hig_cand_dm[0]<40&&"
  //  "((nbt>=2&&nbm>=3&&nbl>=4))",
  //  search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_hig_cand_am").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(20,0,4,"hig_cand_drmax[0]", "#DeltaR_{max}", {1.1, 2.2}),
  //  base_filters &&
  //  "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
  //  "hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
  //  "((nbt>=2&&nbm>=3&&nbl>=4))",
  //  search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_hig_cand_drmax").LuminosityTag(total_luminosity_string);

  //generate pie-charts

  if (HigUtilities::is_in_string_options(options.string_options,"make_pies")) {
    //// 2b (met: 150, 200, 300, 400) low drmax
    pm.Push<Table>("FixName:selection__search_pies__2b_met150_lowdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt==2&&nbm==2)&&met>150 &&met<=200 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
    pm.Push<Table>("FixName:selection__search_pies__2b_met200_lowdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt==2&&nbm==2)&&met>200 &&met<=250 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
    pm.Push<Table>("FixName:selection__search_pies__2b_met300_lowdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt==2&&nbm==2)&&met>300 &&met<=400 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
    pm.Push<Table>("FixName:selection__search_pies__2b_met400_lowdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt==2&&nbm==2)&&met>400            &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
    // 3b (met: 150, 200, 300, 400) low drmax
    pm.Push<Table>("FixName:selection__search_pies__3b_met150_lowdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>150 &&met<=200 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
    pm.Push<Table>("FixName:selection__search_pies__3b_met200_lowdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>200 &&met<=250 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
    pm.Push<Table>("FixName:selection__search_pies__3b_met300_lowdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>300 &&met<=400 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
    pm.Push<Table>("FixName:selection__search_pies__3b_met400_lowdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>400            &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
    // 4b (met: 150, 200, 300, 400) low drmax
    pm.Push<Table>("FixName:selection__search_pies__4b_met150_lowdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>150 &&met<=200 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
    pm.Push<Table>("FixName:selection__search_pies__4b_met200_lowdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>200 &&met<=250 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
    pm.Push<Table>("FixName:selection__search_pies__4b_met300_lowdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>300 &&met<=400 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
    pm.Push<Table>("FixName:selection__search_pies__4b_met400_lowdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>400            &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
    // 2b (met: 150, 200, 300, 400) high drmax
    pm.Push<Table>("FixName:selection__search_pies__2b_met150_highdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt==2&&nbm==2)&&met>150 &&met<=200 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
    pm.Push<Table>("FixName:selection__search_pies__2b_met200_highdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt==2&&nbm==2)&&met>200 &&met<=250 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
    pm.Push<Table>("FixName:selection__search_pies__2b_met300_highdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt==2&&nbm==2)&&met>300 &&met<=400 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
    pm.Push<Table>("FixName:selection__search_pies__2b_met400_highdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt==2&&nbm==2)&&met>400            &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
    // 3b (met: 150, 200, 300, 400) high drmax
    pm.Push<Table>("FixName:selection__search_pies__3b_met150_highdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>150 &&met<=200 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
    pm.Push<Table>("FixName:selection__search_pies__3b_met200_highdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>200 &&met<=250 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
    pm.Push<Table>("FixName:selection__search_pies__3b_met300_highdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>300 &&met<=400 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
    pm.Push<Table>("FixName:selection__search_pies__3b_met400_highdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>400            &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
    // 4b (met: 150, 200, 300, 400) high drmax
    pm.Push<Table>("FixName:selection__search_pies__4b_met150_highdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>150 &&met<=200 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
    pm.Push<Table>("FixName:selection__search_pies__4b_met200_highdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>200 &&met<=250 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
    pm.Push<Table>("FixName:selection__search_pies__4b_met300_highdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>300 &&met<=400 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
    pm.Push<Table>("FixName:selection__search_pies__4b_met400_highdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>400            &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
    // 2b3b4b (met: 150, 200, 300, 400) low drmax
    pm.Push<Table>("FixName:selection__search_pies__2b3b4b_met150_lowdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"nbt>=2&&met>150 &&met<=200 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
    pm.Push<Table>("FixName:selection__search_pies__2b3b4b_met200_lowdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"nbt>=2&&met>200 &&met<=250 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
    pm.Push<Table>("FixName:selection__search_pies__2b3b4b_met300_lowdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"nbt>=2&&met>300 &&met<=400 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
    pm.Push<Table>("FixName:selection__search_pies__2b3b4b_met400_lowdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"nbt>=2&&met>400            &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
    // 2b3b4b (met: 150, 200, 300, 400) high drmax
    pm.Push<Table>("FixName:selection__search_pies__2b3b4b_met150_highdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"nbt>=2&&met>150 &&met<=200 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
    pm.Push<Table>("FixName:selection__search_pies__2b3b4b_met200_highdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"nbt>=2&&met>200 &&met<=250 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
    pm.Push<Table>("FixName:selection__search_pies__2b3b4b_met300_highdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"nbt>=2&&met>300 &&met<=400 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
    pm.Push<Table>("FixName:selection__search_pies__2b3b4b_met400_highdrmax_"+options.year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"nbt>=2&&met>400            &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
  }

  if (HigUtilities::is_in_string_options(options.string_options,"plot_am")) {
    //SR <m> plots
    // 3b4b [<m>] (met: 150, 200, 300, 400) low drmax
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt>=2&&nbm>=3&&hig_cand_drmax[0]<=1.1 && met>150 && met<=200", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_3b4b_met150_lowdrmax_"+options.year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt>=2&&nbm>=3&&hig_cand_drmax[0]<=1.1 && met>200 && met<=300", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_3b4b_met200_lowdrmax_"+options.year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt>=2&&nbm>=3&&hig_cand_drmax[0]<=1.1 && met>300 && met<=400", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_3b4b_met300_lowdrmax_"+options.year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt>=2&&nbm>=3&&hig_cand_drmax[0]<=1.1 && met>400            ", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_3b4b_met400_lowdrmax_"+options.year_string).LuminosityTag(total_luminosity_string);
    // 3b4b [<m>] (met: 150, 200, 300, 400) high drmax
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt>=2&&nbm>=3&&hig_cand_drmax[0]>1.1 && met>150 && met<=200", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_3b4b_met150_highdrmax_"+options.year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt>=2&&nbm>=3&&hig_cand_drmax[0]>1.1 && met>200 && met<=300", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_3b4b_met200_highdrmax_"+options.year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt>=2&&nbm>=3&&hig_cand_drmax[0]>1.1 && met>300 && met<=400", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_3b4b_met300_highdrmax_"+options.year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt>=2&&nbm>=3&&hig_cand_drmax[0]>1.1 && met>400            ", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_3b4b_met400_highdrmax_"+options.year_string).LuminosityTag(total_luminosity_string);

    // 2b [<m>] (met: 150, 200, 300, 400) low drmax
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt==2&&nbm==2&&hig_cand_drmax[0]<=1.1 && met>150 && met<=200", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_2b_met150_lowdrmax_"+options.year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt==2&&nbm==2&&hig_cand_drmax[0]<=1.1 && met>200 && met<=300", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_2b_met200_lowdrmax_"+options.year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt==2&&nbm==2&&hig_cand_drmax[0]<=1.1 && met>300 && met<=400", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_2b_met300_lowdrmax_"+options.year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt==2&&nbm==2&&hig_cand_drmax[0]<=1.1 && met>400            ", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_2b_met400_lowdrmax_"+options.year_string).LuminosityTag(total_luminosity_string);
    // 2b [<m>] (met: 150, 200, 300, 400) high drmax
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt==2&&nbm==2&&hig_cand_drmax[0]>1.1 && met>150 && met<=200", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_2b_met150_highdrmax_"+options.year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt==2&&nbm==2&&hig_cand_drmax[0]>1.1 && met>200 && met<=300", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_2b_met200_highdrmax_"+options.year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt==2&&nbm==2&&hig_cand_drmax[0]>1.1 && met>300 && met<=400", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_2b_met300_highdrmax_"+options.year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt==2&&nbm==2&&hig_cand_drmax[0]>1.1 && met>400            ", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_2b_met400_highdrmax_"+options.year_string).LuminosityTag(total_luminosity_string);

  }

  pm.multithreaded_ = !options.single_thread;
  pm.min_print_ = true;
  pm.MakePlots(1.);

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

