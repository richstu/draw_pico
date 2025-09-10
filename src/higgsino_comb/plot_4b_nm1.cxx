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
#include "TSystem.h"
#include "TROOT.h"
#include "TMVA/Reader.h"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/slide_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/utilities.hpp"
#include "core/sample_loader.hpp"
#include "higgsino_comb/vardef.hpp"
#include "higgsino_comb/regions.hpp"
#include "higgsino_comb/weights.hpp"
#include "higgsino_comb/paths.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace vardef;
using namespace regions_4b;
using namespace weights;

// temporary fix for stitch problems (will be removed after new n2p production)
  NamedFunc TTJ_stitch("TTJ_stitch", [](const Baby &b) -> NamedFunc::ScalarType{
    std::string file = *b.FileNames().begin();
    bool isTTJets_incl = false;
    bool isgenmet = false;
    bool use_event = true;
    if (file.find("TTJets_") != file.npos) isTTJets_incl = true;
    if (file.find("genMET") != file.npos) isgenmet = true;
    if (isTTJets_incl && !isgenmet) {
      if (b.met_tru() > 150.0) {
        use_event = false;
      }
     }
    return use_event;
  });

  NamedFunc WJets_stitch("WJets_stitch", [](const Baby &b) -> NamedFunc::ScalarType{
    std::string file = *b.FileNames().begin();
    bool isWJets_incl = false;
    bool use_event = true;
    if (file.find("_WJetsToLNu_TuneCP5") != file.npos) isWJets_incl = true;
    if (isWJets_incl) { // targeting TTJets inclusive samples
      if (b.ht_isr_me() > 70.0) {
        use_event = false;
      }
     }
    return use_event;
  });


int main() {
  gErrorIgnoreLevel = 6000;

  //trigger
  const NamedFunc ptmiss_htmiss_trig = triggers::ptmiss_htmiss_trig;
  const NamedFunc pass_filters = filters::pass_filters;

  //Load sample paths
  auto signal_info = paths_4b::signal_info;
  auto sample_info_met150 = paths_4b::sample_info_met150;

  //Define processes
  // Iterate over the sample info map
  vector<shared_ptr<Process>> procs_vec;

  for (const auto& pair : sample_info_met150) {
    const auto& searchstring = pair.first;
    const auto& ProcInfo = pair.second;
 
   //Pick processes to ignore
    if (ProcInfo.proclegend.find("Data") != ProcInfo.proclegend.npos) continue;

    std::set<string> fullpath;
    for (const auto& proc : searchstring){ 
      string path(ProcInfo.basepath + "*" + proc + "*");
      fullpath.insert(path);
      }
   
    //Other processes
      auto proc = Process::MakeShared<Baby_pico>(ProcInfo.proclegend, ProcInfo.type, TColor::GetColor(ProcInfo.proccolor.c_str()), {fullpath}, ProcInfo.special_proccuts);
      procs_vec.push_back(proc);
  }
 
 for (const auto& pair : signal_info) {
    
    const auto& searchstring = pair.first;
    const auto& ProcInfo = pair.second;
    
    string fullpath(ProcInfo.basepath + "*" + searchstring + "*");    
    //Other processes
      auto proc = Process::MakeShared<Baby_pico>(ProcInfo.proclegend, ProcInfo.type, TColor::GetColor(ProcInfo.proccolor.c_str()), {fullpath}, ProcInfo.special_proccuts);
      procs_vec.push_back(proc);
  }
  
  //Define plot options
  PlotOpt cmspaper("txt/plot_styles.txt","CMSPaper");
  cmspaper.Title(TitleType::info) //simulation, preliminary
          .Overflow(OverflowType::none)
          .FileExtensions({"pdf"})
          .Stack(StackType::signal_overlay); //other options: shapes, data_norm
//  vector<PlotOpt> ops = {cmspaper.YAxis(YAxisType::linear), cmspaper.YAxis(YAxisType::log)};
  vector<PlotOpt> ops = {cmspaper.YAxis(YAxisType::linear)};
  vector<PlotOpt> ops_log = {cmspaper.YAxis(YAxisType::log)};

  //Plots
  PlotMaker pm;
  
  pm.Push<Hist1D>(Axis(45,150,600,"met", "p_{T}^{miss} [GeV]", {}),		   res_baseline && nb4 && ptmiss_htmiss_trig && pass_filters && TTJ_stitch && WJets_stitch,    procs_vec, ops).Weight(w_lumi_hb).Tag("ShortName:res_met_nm1_nb4_ny0");
  pm.Push<Hist1D>(Axis(4,0,4,"nvlep", "n_{veto leptons}", {1}), 		   res_nm_nvl && nb4 && ptmiss_htmiss_trig && pass_filters && TTJ_stitch && WJets_stitch,      procs_vec, ops).Weight(w_lumi_hb).Tag("ShortName:res_nvl_nm1_nb4_ny0");
  pm.Push<Hist1D>(Axis(6,4,10,"njet", "n_{jets}", {4,6}),       		   res_nm_njet && nb4 && ptmiss_htmiss_trig && pass_filters && TTJ_stitch && WJets_stitch,     procs_vec, ops).Weight(w_lumi_hb).Tag("ShortName:res_njet_nm1_nb4_ny0");
  pm.Push<Hist1D>(Axis(5,0,5, num_b, "n_{b}", {2}),		    		   res_nm_nb  && ptmiss_htmiss_trig && pass_filters && TTJ_stitch && WJets_stitch,             procs_vec, ops).Weight(w_lumi_hb).Tag("ShortName:res_nb_nm1_ny0");
  pm.Push<Hist1D>(Axis(5,0,5, "nphoton", "n_{photons}", {1}),			   res_nm_nphoton && nb4 && ptmiss_htmiss_trig && pass_filters && TTJ_stitch && WJets_stitch,  procs_vec, ops).Weight(w_lumi_hb).Tag("ShortName:res_ny_nm1");
  pm.Push<Hist1D>(Axis(25,0,400, "hig_df_cand_am[0]", "<m_{bb}> [GeV]", {200}),    res_nm_higam && nb4 && ptmiss_htmiss_trig && pass_filters && TTJ_stitch && WJets_stitch,    procs_vec, ops).Weight(w_lumi_hb).Tag("ShortName:res_am_nm1_nb4_ny0");
  pm.Push<Hist1D>(Axis(15,0,120, "hig_df_cand_dm[0]", "#Deltam_{HH} [GeV]", {40}), res_nm_higdm && nb4 && ptmiss_htmiss_trig && pass_filters && TTJ_stitch && WJets_stitch,    procs_vec, ops).Weight(w_lumi_hb).Tag("ShortName:res_dm_nm1_nb4_ny0");
  pm.Push<Hist1D>(Axis(20,0,4, "hig_df_cand_drmax[0]", "#DeltaR_{max}", {2.2}),    res_nm_higdrmax && nb4 && ptmiss_htmiss_trig && pass_filters && TTJ_stitch && WJets_stitch, procs_vec, ops).Weight(w_lumi_hb).Tag("ShortName:res_drmax_nm1_nb4_ny0");

  pm.Push<Hist1D>(Axis(45,150,600,"met", "p_{T}^{miss} [GeV]", {}),		   res_baseline && nb4 && ptmiss_htmiss_trig && pass_filters && TTJ_stitch && WJets_stitch,    procs_vec, ops_log).Weight(w_lumi_hb).Tag("ShortName:res_log_met_nm1_nb4_ny0");
  pm.Push<Hist1D>(Axis(4,0,4,"nvlep", "Number of veto leptons", {1}), 		   res_nm_nvl && nb4 && ptmiss_htmiss_trig && pass_filters && TTJ_stitch && WJets_stitch,      procs_vec, ops_log).Weight(w_lumi_hb).Tag("ShortName:res_log_nvl_nm1_nb4_ny0");
  pm.Push<Hist1D>(Axis(6,4,10,"njet", "Number of jets", {4,6}),       		   res_nm_njet && nb4 && ptmiss_htmiss_trig && pass_filters && TTJ_stitch && WJets_stitch,     procs_vec, ops_log).Weight(w_lumi_hb).Tag("ShortName:res_log_njet_nm1_nb4_ny0");
  pm.Push<Hist1D>(Axis(5,0,5, num_b, "n_{b}", {2}),		    		   res_nm_nb  && ptmiss_htmiss_trig && pass_filters && TTJ_stitch && WJets_stitch,             procs_vec, ops_log).Weight(w_lumi_hb).Tag("ShortName:res_log_nb_nm1_ny0");
  pm.Push<Hist1D>(Axis(5,0,5, "nphoton", "n_{photons}", {1}),			   res_nm_nphoton && nb4 && ptmiss_htmiss_trig && pass_filters && TTJ_stitch && WJets_stitch,  procs_vec, ops_log).Weight(w_lumi_hb).Tag("ShortName:res_log_ny_nm1");
  pm.Push<Hist1D>(Axis(25,0,400, "hig_df_cand_am[0]", "<m_{bb}> [GeV]", {200}),    res_nm_higam && nb4 && ptmiss_htmiss_trig && pass_filters && TTJ_stitch && WJets_stitch,    procs_vec, ops_log).Weight(w_lumi_hb).Tag("ShortName:res_log_am_nm1_nb4_ny0");
  pm.Push<Hist1D>(Axis(15,0,120, "hig_df_cand_dm[0]", "#Deltam_{HH} [GeV]", {40}), res_nm_higdm && nb4 && ptmiss_htmiss_trig && pass_filters && TTJ_stitch && WJets_stitch,    procs_vec, ops_log).Weight(w_lumi_hb).Tag("ShortName:res_log_dm_nm1_nb4_ny0");
  pm.Push<Hist1D>(Axis(20,0,4, "hig_df_cand_drmax[0]", "#DeltaR_{max}", {2.2}),    res_nm_higdrmax && nb4 && ptmiss_htmiss_trig && pass_filters && TTJ_stitch && WJets_stitch, procs_vec, ops_log).Weight(w_lumi_hb).Tag("ShortName:res_log_drmax_nm1_nb4_ny0");



//  pm.multithreaded_ = true;
  pm.min_print_ = true;
  pm.MakePlots(1);
  
}






