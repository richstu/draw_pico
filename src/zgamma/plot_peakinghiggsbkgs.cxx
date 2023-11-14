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
#include "core/utilities.hpp"
#include "zgamma/zg_utilities.hpp"
#include "zgamma/zg_functions.hpp"
using namespace std;
using namespace PlotOptTypes;
using namespace ZgUtilities;
using namespace ZgFunctions;


int main() {

  gErrorIgnoreLevel = 6000;
  Palette col("txt/colors.txt","default");
  Process::Type back =  Process::Type::background;
  //Process::Type sig  =  Process::Type::signal;

  //Folders containing each higgs sample
  string bfolder("/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v0/");
  string spath_2016(bfolder+"/2016/mc/unskimmed/");
  string spath_2016APV(bfolder+"/2016APV/mc/unskimmed/");
  string spath_2017(bfolder+"/2017/mc/unskimmed/");
  string spath_2018(bfolder+"/2018/mc/unskimmed/");

  //Procs for each higgs mode that we want to check (HtoMuMu does not appear for NanoAODv9 it seems)
  auto proc_hmumu   = Process::MakeShared<Baby_pico>("H #rightarrow #mu#mu",       back, TColor::GetColor("#c39de2"), 
                      {spath_2016+"*HToMuMu*.root", spath_2016APV+"*HToMuMu*.root", spath_2017+"*HToMuMu*.root", spath_2018+"*HToMuMu*.root"},   "1");

  auto proc_hgg     = Process::MakeShared<Baby_pico>("H #rightarrow #gamma#gamma", back, TColor::GetColor("#f3ff00"),
                      {spath_2016+"*HToGG*.root", spath_2016APV+"*HToGG*.root", spath_2017+"*HToGG*.root", spath_2018+"*HToGG*.root"},     "1");

  auto proc_htt     = Process::MakeShared<Baby_pico>("H #rightarrow #tau#tau",     back, TColor::GetColor("#ff9632"),
                      {spath_2016+"*HToTauTau*.root", spath_2016APV+"*HToTauTau*.root", spath_2017+"*HToTauTau*.root", spath_2018+"*HToTauTau*.root"}, "1");

  auto proc_hWW     = Process::MakeShared<Baby_pico>("H #rightarrow WW",           back, TColor::GetColor("#1a76a9"), 
                      {spath_2016+"*HToWW*.root", spath_2016APV+"*HToWW*.root", spath_2017+"*HToWW*.root", spath_2018+"*HToWW*.root"},     "1");

  auto proc_hZZ     = Process::MakeShared<Baby_pico>("H #rightarrow ZZ",           back, TColor::GetColor("#08bab7"), 
                      {spath_2016+"*HToZZ*.root", spath_2016APV+"*HToZZ*.root", spath_2017+"*HToZZ*.root", spath_2018+"*HToZZ*.root"},     "1");

  //vector containing the procs
  vector<shared_ptr<Process>> procs = {proc_hmumu,proc_hgg,proc_htt,proc_hWW,proc_hZZ};


  vector<NamedFunc> lep  =  {"1","ll_lepid[0]==11","ll_lepid[0]==13"};
  vector<string>    slep =  {"ll","el","mu"};

  //Baseline selection for analysis
  vector<NamedFunc> cutflow = {"1", pass_trigs_and_pt, "nll>0", "nllphoton > 0","ll_m[0] > 50", "photon_pt[0] > 15", "photon_pt[0]/llphoton_m[0] > 15.0/110.0", "llphoton_m[0] > 100 && llphoton_m[0] < 180", 
                               "ll_m[0] + llphoton_m[0] > 185", "photon_id80[0]" ,"ll_m[0] > 80 && ll_m[0] < 100" };
  vector<string> cutflow_lab = {"All events", "Triggers and $p_{T}(l)$ selections", "$N_{ll} > 0$", "$m_{ll} >$ 50 GeV" ,"$N_{ll\\gamma} > 0$", "$p_{T}(\\gamma) > 15$ GeV", 
                                "$p_{T}(\\gamma)/m_{ll\\gamma} > 15/110$", "100 GeV $< m_{ll\\gamma} < 180$ GeV", "$m_{ll} + m_{ll\\gamma} > 185$ GeV", "\\gamma IDMVA WP80" ,"80 GeV $< m_{ll} < 100$ GeV "};

  for(unsigned int idx_i = 2; idx_i < cutflow.size(); idx_i++){
    cutflow[idx_i] = cutflow[idx_i - 1] && cutflow[idx_i];
  }
  NamedFunc r2_baseline    = r2_baseline;
  NamedFunc tight_baseline = tightened_baseline;


  //Plot options used in this code
  PlotOpt lin_stack("txt/plot_styles.txt","CMSPaper");
  lin_stack.Title(TitleType::info).Stack(StackType::signal_overlay)
           .UseCMYK(false) .FileExtensions({"pdf"});
  vector<PlotOpt> ops = {lin_stack};
  vector<PlotOpt> ops_ls = {lin_stack().Stack(StackType::lumi_shapes)};
  vector<PlotOpt> logops = {lin_stack().YAxis(YAxisType::log)};
 
  //Start making plots
  PlotMaker pm;

  //Loop containing which plots are made
  for(unsigned int idx_i(0); idx_i < lep.size(); idx_i++) { 
    pm.Push<Hist1D>(Axis(80,100,180, "llphoton_m[0]",    "m_{ll#gamma} [GeV]", {122,128}), r2_baseline && lep.at(idx_i), procs, ops).Weight(wgt_1invfb).Tag("ShortName:hpb_");
    pm.Push<Hist1D>(Axis(80,100,180, "llphoton_m[0]",    "m_{ll#gamma} [GeV]", {122,128}), r2_baseline && lep.at(idx_i), procs, logops).Weight(wgt_1invfb).Tag("ShortName:hpb_");
    pm.Push<Hist1D>(Axis(80,100,180, "llphoton_m[0]",    "m_{ll#gamma} [GeV]", {122,128}), r2_baseline && lep.at(idx_i), procs, ops_ls).Weight(wgt_1invfb).Tag("ShortName:hpb_");

    pm.Push<Hist1D>(Axis(80,100,180, "llphoton_m[0]",    "m_{ll#gamma} [GeV]", {122,128}), tight_baseline && lep.at(idx_i), procs, ops).Weight(wgt_1invfb).Tag("ShortName:hpb_tight");
    pm.Push<Hist1D>(Axis(80,100,180, "llphoton_m[0]",    "m_{ll#gamma} [GeV]", {122,128}), tight_baseline && lep.at(idx_i), procs, logops).Weight(wgt_1invfb).Tag("ShortName:hpb_tight");
    pm.Push<Hist1D>(Axis(80,100,180, "llphoton_m[0]",    "m_{ll#gamma} [GeV]", {122,128}), tight_baseline && lep.at(idx_i), procs, ops_ls).Weight(wgt_1invfb).Tag("ShortName:hpb_tight");

    pm.Push<Hist1D>(Axis(80,100,180, "ll_m[0]",          "m_{ll} [GeV]" ,{}),  r2_baseline && lep.at(idx_i), procs, ops).Weight(wgt_1invfb).Tag("ShortName:hpb_");
    pm.Push<Hist1D>(Axis(20,0,1,     "photon_reliso[0]", "Photon reliso" ,{}), r2_baseline && lep.at(idx_i), procs, logops).Weight(wgt_1invfb).Tag("ShortName:hpb_");
    pm.Push<Hist1D>(Axis(20,0,1,     "photon_r9[0]",     "Photon R9" ,{}),     r2_baseline && lep.at(idx_i), procs, logops).Weight(wgt_1invfb).Tag("ShortName:hpb_");
    pm.Push<Hist1D>(Axis(20,0,0.5,   "photon_hoe[0]",    "Photon H/E" ,{}),    r2_baseline && lep.at(idx_i), procs, ops).Weight(wgt_1invfb).Tag("ShortName:hpb_");
    pm.Push<Hist1D>(Axis(16,0.2,1,   "photon_idmva[0]",  "Photon MVA" ,{}),    r2_baseline && lep.at(idx_i), procs, logops).Weight(wgt_1invfb).Tag("ShortName:hpb_");

  }

  NamedFunc hwin = "llphoton_m[0] > 122 && llphoton_m[0] < 128";
  const NamedFunc pass_trigs_and_pt =( (HLT_pass_singlemuon && "nmu > 0 && mu_pt[0] > 30") || (HLT_pass_singleelectron && "nel > 0 && el_pt[0] > 37")) ||
                                ( (HLT_pass_dimuon && "nmu > 1 && mu_pt[0] > 20 && mu_pt[1] > 10") || (HLT_pass_dielectron && "nel > 1 && el_pt[0] > 25 && el_pt[1] > 15"));

  NamedFunc trigs_ee   = (HLT_pass_singleelectron && "nel > 0 && el_pt[0] > 37") || (HLT_pass_dielectron && "nel > 1 && el_pt[0] > 25 && el_pt[1] > 15");
  NamedFunc trigs_mumu = (HLT_pass_singlemuon && "nmu > 0 && mu_pt[0] > 30") || (HLT_pass_dimuon && "nmu > 1 && mu_pt[0] > 20 && mu_pt[1] > 10");

  //Declare rows used for each table
  vector<vector<TableRow>> tableRows = {{},{},{}};
  for(unsigned int idx_fl = 0; idx_fl < lep.size(); idx_fl++){
    for(unsigned int idx_tab = 0; idx_tab < cutflow.size(); idx_tab++){

      //Add row for all events 
      if(idx_tab==0){
        tableRows[idx_fl].push_back(TableRow(cutflow_lab[idx_tab], cutflow[idx_tab], 0,0, wgt_1invfb)); continue;
      }

      //Add row that only includes the triggers
      if(idx_tab==1){
        if(idx_fl==1){
          tableRows[idx_fl].push_back(TableRow(cutflow_lab[idx_tab], trigs_ee, 0,0, wgt_1invfb)); continue;
        } else if(idx_fl==2){
          tableRows[idx_fl].push_back(TableRow(cutflow_lab[idx_tab], trigs_mumu, 0,0, wgt_1invfb)); continue;
        }
      }

      //Add the other rows for each table
      tableRows[idx_fl].push_back(TableRow(cutflow_lab[idx_tab], cutflow[idx_tab] && lep[idx_fl], 0,0, wgt_1invfb));
      if(idx_tab==cutflow.size()-1){
        tableRows[idx_fl].push_back(TableRow("Run 2: 122 GeV $ < m_{ll\\gamma} < 128 GeV", cutflow[idx_tab-2] && lep[idx_fl] && hwin, 0,0, wgt_1invfb));
        tableRows[idx_fl].push_back(TableRow("Tightened Baseline: 122 GeV $ < m_{ll\\gamma} < 128 GeV", cutflow[idx_tab] && lep[idx_fl] && hwin, 0,0, wgt_1invfb));
      }
    }
  }

  //Pushing the tables with split on lepton flavor
  pm.Push<Table>("HiggsPeakingBkg_All_",   tableRows[0], procs).Precision(3);
  pm.Push<Table>("HiggsPeakingBkg_Zee_",   tableRows[1], procs).Precision(3);
  pm.Push<Table>("HiggsPeakingBkg_Zmumu_", tableRows[2], procs).Precision(3);

  //Code that pushes the tables
  pm.min_print_ = true;
  pm.MakePlots(137.61);
}
