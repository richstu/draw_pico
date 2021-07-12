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
#include "higgsino/rdf_hig_functions.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"
#include "higgsino/ordered_dict.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace Higfuncs;

namespace{
  bool single_thread = false;
  string year_string = "run2";
  bool unblind = false;
  bool unblind_signalregion = false;
  // string_options is split by comma. ex) option1,option2 
  // Use HigUtilities::is_in_string_options(string_options, "option2") to check if in string_options.
  string string_options = "plot_in_btags,plot_in_btags_with_met_split,plot_baseline";
  // Other options: plot_in_btags_for_mc,plot_in_btags_split,plot_in_planes,paper_style
  // paper_style - sets plot style and MC categories to match boosted
}

vector<pair<int, int> > parseMassPoints(string const & mass_points_string) {
  vector<string> v_mass_point;
  HigUtilities::stringToVectorString(mass_points_string, v_mass_point, ",");
  vector<pair<int, int> > mass_points;
  for (auto const & mass_point : v_mass_point) {
    vector<string> v_mass;
    HigUtilities::stringToVectorString(mass_point, v_mass, "_");
    int nlsp_mass = stoi(v_mass[0]); int lsp_mass = stoi(v_mass[1]);
    mass_points.push_back({nlsp_mass, lsp_mass});
  }
  return mass_points;
}

void setProcsDict(string const & production, string const & nanoAodFolder, string const & year_string, string const & sample_name, map<string, vector<shared_ptr<Process> > > & procsDict) {
  set<int> years;
  //years_a = {2016};
  HigUtilities::parseYears(year_string, years);

  // Data cuts
  NamedFunc lepton_triggers = Higfuncs::el_trigger || Higfuncs::mu_trigger;
  NamedFunc met_triggers = Higfuncs::met_trigger;;
  NamedFunc triggers_data = "1";
  if (sample_name == "zll") triggers_data = lepton_triggers;
  else if (sample_name == "ttbar") triggers_data = lepton_triggers || met_triggers;
  else if (sample_name == "qcd") triggers_data = met_triggers;
  else  triggers_data = met_triggers;

  // Set folders
  string mc_production_folder = nanoAodFolder+"/"+production;
  string data_production_folder = nanoAodFolder+"/"+production;
  string signal_production_folder = nanoAodFolder+"/"+production;
  string mc_skim_folder, data_skim_folder, signal_skim_folder;
  if (sample_name == "search") {
    //mc_skim_folder = "mc/merged_higmc_higloose/";
    //data_skim_folder = "data/merged_higdata_higloose/";
    //signal_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higloose/";
    mc_skim_folder = "mc/merged_higmc_preselect/";
    data_skim_folder = "data/merged_higdata_preselect/";
    signal_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_preselect/";
  } else if (sample_name == "ttbar") {
    mc_skim_folder = "mc/merged_higmc_higlep1T/";
    data_skim_folder = "data/merged_higdata_higlep1T/";
    signal_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higlep1T/";
  } else if (sample_name == "zll") {
    mc_skim_folder = "mc/merged_higmc_higlep2T/";
    data_skim_folder = "data/merged_higdata_higlep2T/";
    signal_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higlep2T/";
  } else if (sample_name == "qcd") {
    mc_skim_folder = "mc/merged_higmc_higqcd/";
    data_skim_folder = "data/merged_higdata_higqcd/";
    signal_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higqcd/";
  }

  // Set mc file names
  map<string, set<string>> mc_filenames; 
  mc_filenames["tt"]     = set<string>({"*TTJets_*Lept*",
                                  "*_TTZ*.root", "*_TTW*.root",
                                 "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root"});
  mc_filenames["single_t"] = set<string>({"*_ST_*.root"});
  mc_filenames["zjets"]   = set<string>({"*_ZJet*.root", "*DYJetsToLL*.root"});
  mc_filenames["wjets"]   = set<string>({"*_WJetsToLNu*.root"});
  mc_filenames["qcd"]     = set<string>({"*_QCD_HT*"});
  mc_filenames["other"]   = set<string>({"*_WH*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});
  mc_filenames["other_and_single_t"]   = set<string>({"*_WH*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root","*_ST_*.root"});
  mc_filenames["all"] = set<string>({"*TTJets_SingleLept*",
                               "*TTJets_DiLept*",
                               "*_TTZ*.root", "*_TTW*.root",
                               "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root", "*_ST_*.root",
                               "*_WJetsToLNu*.root", "*_ZJet*.root",
                               "*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                               "*_QCD_HT700to1000_*", "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                               "*_QCD_HT2000toInf_*",
                               "*_WH*.root", "*_ZH_HToBB*.root",
                               "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root", "*DYJetsToLL*.root"
  });
  // Set signal file names
  string mass_points_string = "175_0, 500_0, 950_0"; // Original
  //string mass_points_string = "225_0, 700_0, 950_0"; // Below looks good 
  //string mass_points_string = "225_0, 400_0"; // 2016 AN
  //if (unblind && !unblind_signalregion) mass_points_string = "";
  vector<int> sig_colors = {kGreen-6, kRed, kBlue, kOrange}; // Requires mass_points.size() >= sig_colors.size()
  // signal_filenames["TChiHH(nlsp,lsp)"] = {filename}
  torch::OrderedDict<string, set<string> > signal_filenames;
  // mass_points = [[nlsp, lsp]]
  vector<pair<int, int> > mass_points = parseMassPoints(mass_points_string);
  for (auto const & mass_point : mass_points) {
    signal_filenames.insert("TChiHH("+to_string(mass_point.first)+","+to_string(mass_point.second)+")", {"*TChiHH_mChi-"+to_string(mass_point.first)+"_mLSP-"+to_string(mass_point.second)+"_*.root"});
  }

  Palette colors("txt/colors.txt", "default");

  //procsDict["mc_and_sig_and_data"];
  //procsDict["mc_and_sig_and_data"].push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
  //                attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch&&ntrutauh>0"));
  //procsDict["mc_and_sig_and_data"].push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
  //                attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch&&ntrutauh==0"));
  //procsDict["mc_and_sig_and_data"].push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
  //                attach_folder(mc_production_folder, years, mc_skim_folder,mc_filenames["zjets"]),"stitch"));
  //procsDict["mc_and_sig_and_data"].push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
  //                attach_folder(mc_production_folder, years, mc_skim_folder,mc_filenames["wjets"]),"stitch"));
  //procsDict["mc_and_sig_and_data"].push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
  //                attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["single_t"]),"stitch"));
  //procsDict["mc_and_sig_and_data"].push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
  //                attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["qcd"]),"stitch")); 
  //procsDict["mc_and_sig_and_data"].push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
  //                attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["other"]),"stitch"));
  //for (unsigned iMassPoint = 0; iMassPoint < signal_filenames.size(); ++iMassPoint) {
  //  procsDict["mc_and_sig_and_data"].push_back(Process::MakeShared<Baby_pico>(signal_filenames[iMassPoint].key(), Process::Type::signal, 
  //    sig_colors[iMassPoint], attach_folder(signal_production_folder, years, signal_skim_folder, signal_filenames[iMassPoint].value() ), "stitch"));
  //}
  //procsDict["mc_and_sig_and_data"].push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack, 
  //                 attach_folder(data_production_folder, years, data_skim_folder, {"*.root"}), "HLT_PFMET120_PFMHT120_IDTight"));

  //dictionary combines some MC categories for paper
  procsDict["paper_mc_and_sig_and_data"];
  procsDict["paper_mc_and_sig_and_data"].push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background,colors("tt_htau"),
                  attach_folder_years(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch"));
  procsDict["paper_mc_and_sig_and_data"].push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder_years(mc_production_folder, years, mc_skim_folder,mc_filenames["zjets"]),"stitch"));
  procsDict["paper_mc_and_sig_and_data"].push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder_years(mc_production_folder, years, mc_skim_folder,mc_filenames["wjets"]),"stitch"));
  procsDict["paper_mc_and_sig_and_data"].push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder_years(mc_production_folder, years, mc_skim_folder, mc_filenames["qcd"]),"stitch")); 
  procsDict["paper_mc_and_sig_and_data"].push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder_years(mc_production_folder, years, mc_skim_folder, mc_filenames["other_and_single_t"]),"stitch"));
  for (unsigned iMassPoint = 0; iMassPoint < signal_filenames.size(); ++iMassPoint) {
    procsDict["paper_mc_and_sig_and_data"].push_back(Process::MakeShared<Baby_pico>(signal_filenames[iMassPoint].key(), Process::Type::signal, 
      sig_colors[iMassPoint], attach_folder_years(signal_production_folder, years, signal_skim_folder, signal_filenames[iMassPoint].value() ), "stitch"));
  }
  procsDict["paper_mc_and_sig_and_data"].push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack, 
                   attach_folder_years(data_production_folder, years, data_skim_folder, {"*.root"}), "met_trigger"));  //data_trigger

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
    .Stack(StackType::data_norm).LegendColumns(3).ErrorOnZeroData(true).PrintVals(true);
  PlotOpt log_norm_info = lin_norm_info().YAxis(YAxisType::log);
  PlotOpt log_norm = lin_norm_info().YAxis(YAxisType::log).Title(TitleType::info).LogMinimum(.2);
  PlotOpt lin_norm = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info);
  PlotOpt lin_shapes = lin_norm().Stack(StackType::shapes).Bottom(BottomType::ratio);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info).Bottom(BottomType::off);
  PlotOpt log_norm_data = lin_norm_info().YAxis(YAxisType::log).Title(TitleType::info).LogMinimum(.2).Bottom(BottomType::ratio).PrintVals(true);
  PlotOpt lin_norm_data = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info).Bottom(BottomType::ratio).PrintVals(true);
  PlotOpt lin_lumi = lin_norm_info.Title(TitleType::data).Bottom(BottomType::ratio).YAxis(YAxisType::linear).Stack(StackType::data_norm).RatioMaximum(1.86);
  PlotOpt lin_norm_paper = lin_norm().Title(TitleType::data).LegendColumns(2).TitleInFrame(true);
  PlotOpt lin_norm_data_paper = lin_norm_data().Title(TitleType::data).LegendColumns(2).TitleInFrame(true).RatioMaximum(1.9).RatioMinimum(0.1);
  PlotOpt log_norm_paper = log_norm().Title(TitleType::data).LegendColumns(2).TitleInFrame(true);
  PlotOpt log_norm_data_paper = log_norm_data().Title(TitleType::data).LegendColumns(2).TitleInFrame(true).RatioMaximum(1.9).RatioMinimum(0.1);

  vector<PlotOpt> plt_norm_info = {lin_norm_info, log_norm_info};
  vector<PlotOpt> plt_lin = {lin_norm};
  vector<PlotOpt> plt_log = {log_norm};
  vector<PlotOpt> plt_shapes = {lin_shapes};
  vector<PlotOpt> plt_shapes_info = {lin_shapes_info};
  vector<PlotOpt> plt_lin_lumi = {lin_lumi};
  unblind = true;
  if (unblind) plt_lin = {lin_norm_data};
  if (unblind) plt_log = {log_norm_data};
  plt_lin = {lin_norm_paper};
  plt_log = {log_norm_paper};
  if (unblind) {
    plt_lin = {lin_norm_data_paper};
    plt_log = {log_norm_data_paper};
  }

  set<int> years;
  HigUtilities::parseYears(year_string, years);
  int lumi_precision = 1;
  if (HigUtilities::is_in_string_options(string_options, "paper_style"))
    lumi_precision = 0;
  string total_luminosity_string = HigUtilities::getLuminosityString(year_string, lumi_precision);

  // Set folders according
  // production, nanoAODFolder, sample_name, year_string, 
  // Set procs
  // Set baseline, filter according: sample_name
  string production = "higgsino_klamath"; 
  string nanoAodFolder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7";

  // [mc, mcTtbar, mc_and_sig, mc_and_sig_and_data, mc_btag, mc_nisr, mcTtbar_met, mcTtbar_lowmet, mc_and_data, data_3btag]
  map<string, vector<shared_ptr<Process> > > procs_search;
  setProcsDict(production, nanoAodFolder, year_string, "search", procs_search); 

  NamedFunc weight = "final_weight_run2"; 

  // Filters for each sample
  NamedFunc search_filters = "final_pass_filters"; 

  // resolved cuts
  NamedFunc search_resolved_cuts = 
                         "met/mht<2 && met/met_calo<2&&weight<1.5&&"
                         "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
                         "hig_cand_drmax[0]<=2.2&&hig_cand_am[0]<=200&&hig_cand_dm[0]<=40&&"
                         "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))";

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

  string proc_name = "paper_mc_and_sig_and_data";

  NamedFunc extra_cut = "1";
  NamedFunc signalReject = "njet>=4&&njet<=5&&!(((nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))&&(hig_cand_am[0]>100&&hig_cand_am[0]<=140))";

  //pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "njet", "N_{j}", {}),
  //  "met>250"
  //  ,procs_search[proc_name], plt_log).Weight("weight").Tag("FixName:debug__nj").LuminosityTag(total_luminosity_string);
  //// nb
  pm.Push<Hist1D>(Axis(3, 1.5, 4.5, "hig_bcat", "N_{b}", {2.5}),
    search_filters&&search_resolved_cuts
    ,procs_search[proc_name], plt_log).Weight(weight).Tag("FixName:nb__search").LuminosityTag(total_luminosity_string);
  // Plot according to btags
  //vector<pair<string, NamedFunc> > nbCuts = {{"2b3b4b","hig_bcat==2||hig_bcat==3||hig_bcat==4"}, {"3b4b","hig_bcat==3||hig_bcat==4"}, {"4b","hig_bcat==4"}};
  //vector<pair<string, NamedFunc> > nbCuts = {};
  vector<pair<string, NamedFunc> > nbCuts = {{"3b4b","hig_bcat==3||hig_bcat==4"}};
  for (unsigned iCut = 0; iCut < nbCuts.size(); ++iCut) {
    string nbCutName = nbCuts[iCut].first;
    NamedFunc nbCut = nbCuts[iCut].second;
    // am
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am_0", "<m_{bb}> [GeV]", {100, 140}),
      search_filters&&
      "met/mht<2 && met/met_calo<2&&weight<1.5&&"
      "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
      "hig_cand_drmax[0]<=2.2&&hig_cand_dm[0]<=40"&&nbCut
      , procs_search[proc_name], plt_lin).Weight(weight).Tag("FixName:amjj__"+nbCutName+"__search").LuminosityTag(total_luminosity_string);
    // dm
    pm.Push<Hist1D>(Axis(10,0,100,"hig_cand_dm_0", "#Deltam [GeV]", {40.}),
      search_filters&&
      "met/mht<2 && met/met_calo<2&&weight<1.5&&"
      "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
      "hig_cand_drmax[0]<=2.2&&hig_cand_am[0]<=200"&&nbCut
      , procs_search[proc_name], plt_lin).Weight(weight).Tag("FixName:dmjj__"+nbCutName+"__search").LuminosityTag(total_luminosity_string);
    // Delta R max
    pm.Push<Hist1D>(Axis(20,0,4,"hig_cand_drmax_0", "#DeltaR_{max}", {1.1, 2.2}),
      search_filters&&
      "met/mht<2 && met/met_calo<2&&weight<1.5&&"
      "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
      "hig_cand_dm[0]<=40&&hig_cand_am[0]<=200"&&nbCut
      ,procs_search[proc_name], plt_log).Weight(weight).Tag("FixName:drmax__"+nbCutName+"__search_log").LuminosityTag(total_luminosity_string);
    // met
    pm.Push<Hist1D>(Axis(14, 150, 850., "met", "p_{T}^{miss} [GeV]", {200., 300., 400.}),
    search_filters&&search_resolved_cuts&&nbCut
    ,procs_search[proc_name], plt_log).Weight(weight).Tag("FixName:met__"+nbCutName+"__search").LuminosityTag(total_luminosity_string);
  }

  pm.multithreaded_ = !single_thread;
  pm.min_print_ = true;
  pm.DefineRdfColumns(RdfHigfuncs::rdf_hig_functions_slim);
  pm.MakePlotsRdf(1.0);

  //// Get Figures
  //vector<std::unique_ptr<Figure> > const & figures = pm.Figures();
  //for (auto figurePtr = figures.begin(); figurePtr != figures.end(); ++figurePtr) {
  //  Hist1D * figure = dynamic_cast<Hist1D *>(figurePtr->get());
  //  cout<<figure->Name()<<endl;
  //  // Processes
  //}

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
  return 0;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"single_thread", no_argument, 0, 's'},
      {"unblind_signalregion", no_argument, 0, 'a'},
      {"unblind", no_argument, 0, 'u'},
      {"year", required_argument, 0, 'y'},
      {"string_options", required_argument, 0, 'o'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "sauy:o:", long_options, &option_index);

    if( opt == -1) break;

    string optname;
    switch(opt){
    case 's':
      single_thread = true;
      break;
    case 'u':
      unblind = true;
      break;
    case 'a':
      unblind_signalregion = true;
      break;
    case 'y':
      year_string = optarg;
      break;
    case 'o':
      string_options = optarg;
      break;
    case 0:
      optname = long_options[option_index].name;
      if (0) {
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
