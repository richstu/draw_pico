  
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
using namespace regions_4b;
using namespace weights;
using namespace vardef;

vector<pair<string, Plotter>> cutflow_plot_params = {
//    {"met",	  { "met", 		    "p_{T}^{miss} [GeV]", {50, 150, 650}, {150}} },
    {"nvl",	  { "nvlep",		    "n_{veto leptons}",   {6, 0, 6},	  {1}} },
    {"ntk", 	  { "ntk", 		    "n_{tracks}",	  {6, 0, 6},	  {1}} },
    {"njet",	  { "njet", 		    "n_{jets}",	  	  {8, 0, 8},	  {4,6}} },
    {"nb", 	  { num_b, 		    "n_{b}",	  	  {6, 0, 6}, 	  {2}} },
//    {"fakemet",   { {dphi1, dphi2, dphi3, dphi4}, {25, 0, 5}},     {0.5, 0.3} },
    {"nphoton",   { "nphoton", 		    "n_{#gamma#}",	  {3, 0, 3},      {1}} },
    {"hig_am",    { "hig_df_cand_am[0]",    "<m_{bb}> [GeV]", 	  {40, 0, 400},   {200}} },
    {"hig_dm",    { "hig_df_cand_dm[0]",    "#Deltam_{HH} [GeV]", {25, 0, 200},   {40}} },
    {"hig_drmax", { "hig_df_cand_drmax[0]", "#DeltaR_{max}", 	  {20, 0, 4},     {2.2}} },
};



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
  //using ProcInfo = bbgg_paths::ProcInfo; 
  auto signal_info = paths_4b::signal_info;
  auto sample_info_met150 = paths_4b::sample_info_met150;

  const NamedFunc basic_cuts = sig_decay_4b && pass_filters && ptmiss_htmiss_trig && TTJ_stitch && WJets_stitch;
 //Group processes
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
  cmspaper.Title(TitleType::info)
          .Overflow(OverflowType::none)
          .FileExtensions({"pdf"})
          .Stack(StackType::signal_overlay); //other options: shapes, data_norm
  vector<PlotOpt> ops = {cmspaper.YAxis(YAxisType::linear), cmspaper.YAxis(YAxisType::log)};

  //Cutflow table using the cuts vector defined in the regions.cpp file
  PlotMaker pm;
  vector<TableRow> tablevec;
  stringstream cuts;
  int index = 0;
  bool at_fakemet = false;
  //table for hh signal region
  tablevec.push_back(TableRow("Filters, p_{T}^{miss} > 150 GeV", basic_cuts, 0, 0, w_lumi_hb_mettrig));
  // Iterate over the sample info map
  for (const auto& pair : cuts_4b_res) {

    const auto& cutlabel = pair.first;
    const auto& cut = pair.second; 
  // make cutflow plots before new cut is added
    for (const auto& pair2 : cutflow_plot_params){
      const auto& label = pair2.first;
      const auto& plotparams = pair2.second;
      if (label == cutlabel){
        const auto var = plotparams.params;
        const auto axes_label = plotparams.axes_name;
        const auto axes = plotparams.axes_limits;
        const auto lines = plotparams.vert_lines;
        if (at_fakemet){
          pm.Push<Hist1D>(Axis(axes[0],axes[1],axes[2],var,axes_label,lines), cuts.str() && basic_cuts && dphi_res, procs_vec, ops).Weight(w_lumi_hb_mettrig).Tag("ShortName:cutflow_res_"+label);
          }
        else{
          pm.Push<Hist1D>(Axis(axes[0],axes[1],axes[2],var,axes_label,lines), cuts.str() && basic_cuts, procs_vec, ops).Weight(w_lumi_hb_mettrig).Tag("ShortName:cutflow_res_"+label);
          }
        }
      }
  // fakemet cut has four plots, make them separately here
    if (cutlabel == "fakemet"){
      pm.Push<Hist1D>(Axis(25,0,5,dphi1,"#Delta#phi_{j_{1},MET}",{0.5}), cuts.str()&&basic_cuts, procs_vec, ops).Weight(w_lumi_hb_mettrig).Tag("ShortName:cutflow_res_dphi1");
      pm.Push<Hist1D>(Axis(25,0,5,dphi2,"#Delta#phi_{j_{2},MET}",{0.5}), cuts.str()&&basic_cuts, procs_vec, ops).Weight(w_lumi_hb_mettrig).Tag("ShortName:cutflow_res_dphi2");
      pm.Push<Hist1D>(Axis(25,0,5,dphi3,"#Delta#phi_{j_{3},MET}",{0.3}), cuts.str()&&basic_cuts, procs_vec, ops).Weight(w_lumi_hb_mettrig).Tag("ShortName:cutflow_res_dphi3");
      pm.Push<Hist1D>(Axis(25,0,5,dphi4,"#Delta#phi_{j_{4},MET}",{0.3}), cuts.str()&&basic_cuts, procs_vec, ops).Weight(w_lumi_hb_mettrig).Tag("ShortName:cutflow_res_dphi4");
      }
    if (index == 0) { cuts << cut; }
    else { cuts << "&&" << cut; }
    index += 1;  
    
    if (cutlabel == "fakemet"){at_fakemet = true;}
    if (at_fakemet){ // beyond fakemet step, add dphi_res cut
        tablevec.push_back(TableRow(cutlabel, cuts.str() && basic_cuts && dphi_res, 0, 0, w_lumi_hb_mettrig));
      }
    else {
        tablevec.push_back(TableRow(cutlabel, cuts.str() && basic_cuts, 0, 0, w_lumi_hb_mettrig));
      }

    }
// add two more rows in table to show effect of higher ptmiss cut and 4b requirement
  tablevec.push_back(TableRow("met > 250 GeV", cuts.str() && basic_cuts && dphi_res && "met > 250", 0, 0, w_lumi_hb_mettrig));
  tablevec.push_back(TableRow("nb = 4", cuts.str() && basic_cuts && dphi_res && "met>250" && nb4, 0, 0, w_lumi_hb_mettrig));

  pm.Push<Hist1D>(Axis(50,150,650,"met","p_{T}^{miss} [GeV]",{250}), cuts.str() && basic_cuts && dphi_res && "met>250", procs_vec, ops).Weight(w_lumi_hb_mettrig).Tag("ShortName:cutflow_res_met250");
  pm.Push<Hist1D>(Axis(4,2,6, num_b, "n_{b}", {4}), cuts.str() && basic_cuts && dphi_res && "met>250" && nb4, procs_vec, ops).Weight(w_lumi_hb_mettrig).Tag("ShortName:cutflow_res_nb4");

  pm.Push<Table>("Cutflow_resolved_4b_sigbkgmc_run2_cutflow_resolved_deepflavb_2dxsec_ptmiss_trig_ttjstitch_dyweight_wjstitch", tablevec, procs_vec, false, true, false, false, false, false); //do_zbi, print_table, print_pie, print_titlepie, do_eff, do_unc
  pm.Push<Table>("Cutflow_eff_resolved_4b_sigbkgmc_run2_cutflow_resolved_deepflavb_2dxsec_ptmiss_trig_ttjstitch_dyweight_wjstitch", tablevec, procs_vec, false, true, false, false, true, false); //do_zbi, print_table, print_pie, print_titlepie, do_eff, do_unc

//  pm.multithreaded_ = true;
  pm.min_print_ = true;
  pm.MakePlots(1);
  
}






