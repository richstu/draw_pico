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
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"
#include "core/ordered_dict.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace Higfuncs;

namespace{
  bool single_thread = false;
  string year_string = "run2";
  // string_options is split by comma. ex) option1,option2 
  // Use HigUtilities::is_in_string_options(string_options, "option2") to check if in string_options.
  // Options: use_old_trigger,split_ttbar_met,search_ttbar_same_bins,save_entries_weights_to_file,print_entries_weights,do_zbi,do_zbi_signal
  string string_options = "";
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel = 6000;
  time_t begtime, endtime;
  time(&begtime);
  GetOptions(argc, argv);

  Palette colors("txt/colors.txt", "default");

  PlotOpt lin_norm_info("txt/plot_styles.txt", "CMSPaper");
  lin_norm_info.Title(TitleType::info)   
    .Bottom(BottomType::off)
    .YAxis(YAxisType::linear)
    .Stack(StackType::data_norm).LegendColumns(3);
  PlotOpt log_norm_info = lin_norm_info().YAxis(YAxisType::log);
  PlotOpt log_norm = lin_norm_info().YAxis(YAxisType::log).Title(TitleType::info).LogMinimum(.2);
  PlotOpt lin_norm = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info);
  PlotOpt lin_norm_stack = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info).Stack(StackType::lumi_shapes);
  PlotOpt log_norm_stack = log_norm_info().Title(TitleType::info).Stack(StackType::lumi_shapes);
  PlotOpt lin_shapes = lin_norm().Stack(StackType::shapes).Bottom(BottomType::ratio);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info).Bottom(BottomType::off);
  PlotOpt log_norm_data = lin_norm_info().YAxis(YAxisType::log).Title(TitleType::info).LogMinimum(.2).Bottom(BottomType::ratio).PrintVals(true);
  PlotOpt lin_norm_data = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info).Bottom(BottomType::ratio).PrintVals(true);
  PlotOpt lin_lumi = lin_norm_info.Title(TitleType::data).Bottom(BottomType::ratio).YAxis(YAxisType::linear).Stack(StackType::data_norm).RatioMaximum(1.86);

  vector<PlotOpt> plt_norm_info = {lin_norm_info, log_norm_info};
  vector<PlotOpt> plt_lin = {lin_norm};
  vector<PlotOpt> plt_lin_stack = {lin_norm_stack};
  vector<PlotOpt> plt_log_stack = {log_norm_stack};
  vector<PlotOpt> plt_log = {log_norm};
  vector<PlotOpt> plt_shapes = {lin_shapes};
  vector<PlotOpt> plt_shapes_info = {lin_shapes_info};
  vector<PlotOpt> plt_lin_lumi = {lin_lumi};

  set<int> years;
  HigUtilities::parseYears(year_string, years);
  string total_luminosity_string = HigUtilities::getLuminosityString(year_string);

  // Set folders according
  // production, nanoAODFolder, sample_name, year_string, 
  // Set procs
  // Set baseline, filter according: sample_name
  //string production = "higgsino_inyo"; 
  //string nanoAodFolder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r5/pico/NanoAODv7";
  string production = "higgsino_klamath"; 
  string nanoAodFolder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7";

  // Set folders
  string mc_production_folder = nanoAodFolder+"/"+production;
  string data_production_folder = nanoAodFolder+"/"+production;
  string signal_production_folder = nanoAodFolder+"/"+production;
  // preselect:
  //   ((nbt>=2 && njet>=4 && njet<=5)||(Sum$(fjet_pt>300 && fjet_msoftdrop>50)>1))
  //   nvlep==0 && ntk==0 && !low_dphi_met && met>150 && 
  // higloose: 
  //   (nbt>=2 || nbdft>=2 || Sum$(fjet_pt>300 && fjet_msoftdrop>50)>0)&&
  //   met>150 && nvlep==0
  // higlep1T:
  //   (Sum$(fjet_pt>300 && fjet_msoftdrop>50)>1 || ((nbt>=2 || nbdft>=2) && njet>=4 && njet<=5)) &&
  //   nlep==1 && 
  //   (Max$(el_pt*el_sig)>30 || Max$(mu_pt*mu_sig)>30) 
  // higlep2T:
  //   (Sum$(fjet_pt>300 && fjet_msoftdrop>50)>1 || (njet>=4 && njet<=5))
  //   nlep==2 && 
  //   @ll_m.size()>=1 && Sum$(ll_m>80 && ll_m<100)>=1
  //   (Max$(el_pt*el_sig)>30 || Max$(mu_pt*mu_sig)>30) // pass_2l_trig30
  // higqcd with met150:
  //   (Sum$(fjet_pt>300 && fjet_msoftdrop>50)>1 || (njet>=4 && njet<=5))
  //   nvlep==0 && ntk==0 && low_dphi_met &&
  //   met>150  // Since applied to met150 skim
  //string mc_skim_folder = "mc/merged_higmc_preselect/";
  string mc_skim_folder = "mc/unskimmed/";
  string data_skim_folder = "data/merged_higdata_preselect/";
  string signal_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_preselect/";


  // Data cuts
  NamedFunc lepton_triggers = Higfuncs::el_trigger || Higfuncs::mu_trigger;
  NamedFunc met_triggers = Higfuncs::met_trigger;;
  NamedFunc triggers_data = "1";
  //if (sample_name == "zll") triggers_data = lepton_triggers;
  //else if (sample_name == "ttbar") triggers_data = lepton_triggers || met_triggers;
  //else if (sample_name == "qcd") triggers_data = met_triggers;
  //else  triggers_data = met_triggers;
  triggers_data = met_triggers;

  // resolved cuts
  NamedFunc search_resolved_cuts = 
                         "met/mht<2 && met/met_calo<2&&weight<1.5&&"
                         "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
                         "hig_cand_drmax[0]<=2.2&&hig_cand_am[0]<=200&&hig_cand_dm[0]<=40&&"
                         "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))";
  NamedFunc ttbar_resolved_cuts = 
                         "met/met_calo<5&&weight<1.5&&"
                         "nlep==1&&mt<=100&&njet>=4&&njet<=5&&"
                         "hig_cand_drmax[0]<=2.2&&hig_cand_am[0]<=200&&hig_cand_dm[0]<=40&&"
                         "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))"
                         && Higfuncs::lead_signal_lepton_pt>30;
  NamedFunc zll_resolved_cuts =
                         "met/met_calo<5&&weight<1.5&&"
                         "nlep==2&&njet>=4&&njet<=5&&met<=50&&"
                         "hig_cand_drmax[0]<=2.2&&hig_cand_am[0]<=200&&hig_cand_dm[0]<=40&&"
                         "(nbm==0||nbm==1||nbm==2||nbm>=3)";
  NamedFunc qcd_resolved_cuts =
                         "met/mht<2 && met/met_calo<2&&"
                         "low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
                         "hig_cand_drmax[0]<=2.2&&hig_cand_am[0]<=200&&hig_cand_dm[0]<=40&&"
                         "(nbm==0||nbm==1||nbm==2||nbm>=3)";

  // Weights
  NamedFunc weight = Higfuncs::final_weight; //"weight"*eff_higtrig_run2*w_years*Functions::w_pileup;

  // Filters for each sample
  NamedFunc search_filters = Higfuncs::final_pass_filters; //pass_filters&& "met/mht<2 && met/met_calo<2&&weight<1.5"
  NamedFunc ttbar_filters = Higfuncs::final_ttbar_pass_filters; //pass_filters&& "met/met_calo<5&&weight<1.5"
  NamedFunc zll_filters = Higfuncs::final_zll_pass_filters; //pass_filters&& "met/met_calo<5&&weight<1.5"
  NamedFunc qcd_filters = Higfuncs::final_qcd_pass_filters; //pass_filters&& "met/mht<2 && met/met_calo<2"


  // Set procs
  map<string, vector<shared_ptr<Process> > > procsDict;
  procsDict["mc_ttbars"].push_back(Process::MakeShared<Baby_pico>("TTJets SingleLept", Process::Type::background, 1,
                  attach_folder(mc_production_folder, years, mc_skim_folder, 
                  {"*TTJets_*SingleLept*"}
                  ),"stitch"));
  procsDict["mc_ttbars"].push_back(Process::MakeShared<Baby_pico>("TTJets DiLept", Process::Type::background, 2,
                  attach_folder(mc_production_folder, years, mc_skim_folder, 
                  {"*TTJets_*DiLept*"}
                  ),"stitch"));
  procsDict["mc_ttbars"].push_back(Process::MakeShared<Baby_pico>("TTZToLLNuNu", Process::Type::background, 3,
                  attach_folder(mc_production_folder, years, mc_skim_folder, 
                  {"*_TTZToLLNuNu*.root"}
                  ),"stitch"));
  procsDict["mc_ttbars"].push_back(Process::MakeShared<Baby_pico>("TTG", Process::Type::background, 4,
                  attach_folder(mc_production_folder, years, mc_skim_folder, 
                  {"*_TTGJets*.root"}
                  ),"stitch"));
  procsDict["mc_ttbars"].push_back(Process::MakeShared<Baby_pico>("ttHTobb", Process::Type::background, 5,
                  attach_folder(mc_production_folder, years, mc_skim_folder, 
                  {"*ttHTobb*.root"}
                  ),"stitch"));
  procsDict["mc_ttbars"].push_back(Process::MakeShared<Baby_pico>("TTWJetsToLNu", Process::Type::background, 6,
                  attach_folder(mc_production_folder, years, mc_skim_folder, 
                  {"*_TTWJetsToLNu*.root"}
                  ),"stitch"));
  procsDict["mc_ttbars"].push_back(Process::MakeShared<Baby_pico>("TTTT", Process::Type::background, 7,
                  attach_folder(mc_production_folder, years, mc_skim_folder, 
                  {"*_TTTT*.root"}
                  ),"stitch"));
  procsDict["mc_ttbars"].push_back(Process::MakeShared<Baby_pico>("TTZToQQ", Process::Type::background, 8,
                  attach_folder(mc_production_folder, years, mc_skim_folder, 
                  {"*_TTZToQQ*.root"}
                  ),"stitch"));
  procsDict["mc_ttbars"].push_back(Process::MakeShared<Baby_pico>("TTWJetsQQ", Process::Type::background, 9,
                  attach_folder(mc_production_folder, years, mc_skim_folder, 
                  {"*_TTWJetsToQQ*.root"}
                  ),"stitch"));

  procsDict["mc_vjets"].push_back(Process::MakeShared<Baby_pico>("TTJets SingleLept", Process::Type::background, 1,
                  attach_folder(mc_production_folder, years, mc_skim_folder, 
                  {"*TTJets_*SingleLept*"}
                  ),"stitch"));
  procsDict["mc_vjets"].push_back(Process::MakeShared<Baby_pico>("ZJetsToNuNu", Process::Type::background, 2,
                  attach_folder(mc_production_folder, years, mc_skim_folder, 
                  {"*_ZJetsToNuNu_*"}
                  ),"stitch"));
  procsDict["mc_vjets"].push_back(Process::MakeShared<Baby_pico>("WJetsToLNu", Process::Type::background, 3,
                  attach_folder(mc_production_folder, years, mc_skim_folder, 
                  {"*_WJetsToLNu*"}
                  ),"stitch"));
  procsDict["mc_vjets"].push_back(Process::MakeShared<Baby_pico>("DYJetsToLL", Process::Type::background, 4,
                  attach_folder(mc_production_folder, years, mc_skim_folder, 
                  {"*_DYJetsToLL_*"}
                  ),"stitch"));

  procsDict["mc_vh"].push_back(Process::MakeShared<Baby_pico>("TTJets SingleLept", Process::Type::background, 1,
                  attach_folder(mc_production_folder, years, mc_skim_folder, 
                  {"*TTJets_*SingleLept*"}
                  ),"stitch"));
  procsDict["mc_vh"].push_back(Process::MakeShared<Baby_pico>("ZH HToBB ZToNuNu", Process::Type::background, 3,
                  attach_folder(mc_production_folder, years, mc_skim_folder, 
                  {"*_ZH_HToBB*.root"}
                  ),"stitch"));
  procsDict["mc_vh"].push_back(Process::MakeShared<Baby_pico>("WH HToBB WToLNu", Process::Type::background, 2,
                  attach_folder(mc_production_folder, years, mc_skim_folder, 
                  {"*_WH*.root"}
                  ),"stitch"));

  procsDict["mc_vv"].push_back(Process::MakeShared<Baby_pico>("TTJets SingleLept", Process::Type::background, 1,
                  attach_folder(mc_production_folder, years, mc_skim_folder, 
                  {"*TTJets_*SingleLept*"}
                  ),"stitch"));
  procsDict["mc_vv"].push_back(Process::MakeShared<Baby_pico>("WZTo1L1Nu2Q", Process::Type::background, 2,
                  attach_folder(mc_production_folder, years, mc_skim_folder, 
                  {"*_WZTo1L1Nu2Q*.root"}
                  ),"stitch"));
  procsDict["mc_vv"].push_back(Process::MakeShared<Baby_pico>("WWToLNuQQ", Process::Type::background, 3,
                  attach_folder(mc_production_folder, years, mc_skim_folder, 
                  {"*_WWToLNuQQ*.root"}
                  ),"stitch"));
  procsDict["mc_vv"].push_back(Process::MakeShared<Baby_pico>("WWTo2L2Nu", Process::Type::background, 4,
                  attach_folder(mc_production_folder, years, mc_skim_folder, 
                  {"*_WWTo2L2Nu*.root"}
                  ),"stitch"));
  procsDict["mc_vv"].push_back(Process::MakeShared<Baby_pico>("WZTo1L3Nu", Process::Type::background, 5,
                  attach_folder(mc_production_folder, years, mc_skim_folder, 
                  {"*_WZTo1L3Nu*.root"}
                  ),"stitch"));
  procsDict["mc_vv"].push_back(Process::MakeShared<Baby_pico>("WZTo2L2Q", Process::Type::background, 6,
                  attach_folder(mc_production_folder, years, mc_skim_folder, 
                  {"*_WZTo2L2Q*.root"}
                  ),"stitch"));
  procsDict["mc_vv"].push_back(Process::MakeShared<Baby_pico>("WZTo3LNu", Process::Type::background, 7,
                  attach_folder(mc_production_folder, years, mc_skim_folder, 
                  {"*_WZTo3LNu*.root"}
                  ),"stitch"));
  procsDict["mc_vv"].push_back(Process::MakeShared<Baby_pico>("ZZ", Process::Type::background, 8,
                  attach_folder(mc_production_folder, years, mc_skim_folder, 
                  {"*_ZZ_*.root"}
                  ),"stitch"));

  vector<NamedFunc> bin_cuts_search;
  bin_cuts_search.push_back("met>150&&met<=200 && hig_cand_drmax[0]<=1.1");
  bin_cuts_search.push_back("met>200&&met<=300 && hig_cand_drmax[0]<=1.1");
  bin_cuts_search.push_back("met>300&&met<=400 && hig_cand_drmax[0]<=1.1");
  bin_cuts_search.push_back("met>400           && hig_cand_drmax[0]<=1.1");
  bin_cuts_search.push_back("met>150&&met<=200 && hig_cand_drmax[0]>1.1");
  bin_cuts_search.push_back("met>200&&met<=300 && hig_cand_drmax[0]>1.1");
  bin_cuts_search.push_back("met>300&&met<=400 && hig_cand_drmax[0]>1.1");
  bin_cuts_search.push_back("met>400           && hig_cand_drmax[0]>1.1");

  vector<NamedFunc> bin_cuts_ttbar;
  bin_cuts_ttbar.push_back("met>0&&met<=75 && hig_cand_drmax[0]<=1.1");
  bin_cuts_ttbar.push_back("met>75&&met<=150 && hig_cand_drmax[0]<=1.1");
  bin_cuts_ttbar.push_back("met>150&&met<=200 && hig_cand_drmax[0]<=1.1");
  bin_cuts_ttbar.push_back("met>200           && hig_cand_drmax[0]<=1.1");
  bin_cuts_ttbar.push_back("met>0&&met<=75 && hig_cand_drmax[0]>1.1");
  bin_cuts_ttbar.push_back("met>75&&met<=150 && hig_cand_drmax[0]>1.1");
  bin_cuts_ttbar.push_back("met>150&&met<=200 && hig_cand_drmax[0]>1.1");
  bin_cuts_ttbar.push_back("met>200           && hig_cand_drmax[0]>1.1");

  PlotMaker pm;

  // MET plots for different processes
  pm.Push<Hist1D>(Axis(25, 0, 800., "met", "p_{T}^{miss} [GeV]", {150, 200., 300., 400.}),
    search_filters&&search_resolved_cuts,
    procsDict["mc_ttbars"], plt_log_stack).Weight(weight).Tag("FixName:met_phase_space_ttbars_met_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
  pm.Push<Hist1D>(Axis(25, 0, 800., "met", "p_{T}^{miss} [GeV]", {150, 200., 300., 400.}),
    search_filters&&search_resolved_cuts,
    procsDict["mc_vjets"], plt_log_stack).Weight(weight).Tag("FixName:met_phase_space_vjets_met_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
  pm.Push<Hist1D>(Axis(25, 0, 800., "met", "p_{T}^{miss} [GeV]", {150, 200., 300., 400.}),
    search_filters&&search_resolved_cuts,
    procsDict["mc_vh"], plt_log_stack).Weight(weight).Tag("FixName:met_phase_space_vh_met_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
  pm.Push<Hist1D>(Axis(25, 0, 800., "met", "p_{T}^{miss} [GeV]", {150, 200., 300., 400.}),
    search_filters&&search_resolved_cuts,
    procsDict["mc_vv"], plt_log_stack).Weight(weight).Tag("FixName:met_phase_space_vv_met_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);

  // Table for different processes
  vector<TableRow> cutflow;
  std::size_t lines_before = 0;
  std::size_t lines_after = 0;
  bool do_zbi = false;
  cutflow.push_back(TableRow("No cut", search_filters&&"1", lines_before, lines_after, weight));
  cutflow.push_back(TableRow("$p_T^{\\rm{miss}}>$150 GeV", search_filters&&"met>150", lines_before, lines_after, weight));
  cutflow.push_back(TableRow("Baseline", search_filters&&search_resolved_cuts, lines_before, lines_after, weight));
  pm.Push<Table>("FixName:table_met_phase_space_ttbars", cutflow, procsDict["mc_ttbars"], do_zbi).LuminosityTag(total_luminosity_string);
  pm.Push<Table>("FixName:table_met_phase_space_vjets", cutflow, procsDict["mc_vjets"], do_zbi).LuminosityTag(total_luminosity_string);
  pm.Push<Table>("FixName:table_met_phase_space_vh", cutflow, procsDict["mc_vh"], do_zbi).LuminosityTag(total_luminosity_string);
  pm.Push<Table>("FixName:table_met_phase_space_vv", cutflow, procsDict["mc_vv"], do_zbi).LuminosityTag(total_luminosity_string);

  pm.multithreaded_ = !single_thread;
  pm.min_print_ = true;
  pm.MakePlots(1.);

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"single_thread", no_argument, 0, 's'},
      {"year", required_argument, 0, 0},
      {"string_options", required_argument, 0, 'o'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "so:", long_options, &option_index);

    if( opt == -1) break;

    string optname;
    switch(opt){
    case 's':
      single_thread = true;
      break;
    case 'o':
      string_options = optarg;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "year"){
        year_string = optarg;
      }else{
        printf("Bad option! Found option name %s\n", optname.c_str());
      }
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
