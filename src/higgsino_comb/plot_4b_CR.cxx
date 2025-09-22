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

 NamedFunc met = "met";
 NamedFunc mt = "mt";
 NamedFunc hig_am = "hig_df_cand_am[0]";
 NamedFunc hig_dm = "hig_df_cand_dm[0]";
 NamedFunc hig_drmax = "hig_df_cand_drmax[0]";

 NamedFunc lep1_pt = "lep_pt[0]";
 NamedFunc lep1_eta = "lep_eta[0]";
 NamedFunc lep1_phi = "lep_phi[0]";

 NamedFunc lep_pt1 = "lep_pt[0]";
 NamedFunc lep_pt2 = "lep_pt[1]";
 NamedFunc lep_eta1 = "lep_eta[0]";
 NamedFunc lep_eta2 = "lep_eta[1]";
 NamedFunc lep_phi1 = "lep_phi[0]";
 NamedFunc lep_phi2 = "lep_phi[1]";

 NamedFunc ll_pt = "ll_pt[0]";
 NamedFunc ll_dr = "ll_dr[0]";
 NamedFunc ll_dphi = "ll_dphi[0]";
 NamedFunc ll_deta = "ll_deta[0]";
 NamedFunc ll_m = "ll_m[0]";


int main() {
  gErrorIgnoreLevel = 6000;

  //trigger
  const NamedFunc ptmiss_htmiss_trig = triggers::ptmiss_htmiss_trig;
  const NamedFunc single_ele_trig = triggers::single_ele_trig;
  const NamedFunc single_muon_trig = triggers::single_muon_trig;
  const NamedFunc pass_filters = filters::pass_filters;

  //Load sample paths
  auto sample_info_1lep = paths_4b::sample_info_1lep;
  auto sample_info_2lep = paths_4b::sample_info_2lep;
  auto sample_info_met150 = paths_4b::sample_info_met150;

  //define regions with triggers
  const NamedFunc CR1_el_trig = CR1_el && single_ele_trig;
  const NamedFunc CR1_mu_trig = CR1_mu && single_muon_trig;
  const NamedFunc CR2_el_trig = CR2_el && single_ele_trig;
  const NamedFunc CR2_mu_trig = CR2_mu && single_muon_trig;
  const NamedFunc CR3_met_trig = CR3 && ptmiss_htmiss_trig;

  //Define processes
  // Iterate over the sample info map
  vector<shared_ptr<Process>> procs_vec_met150;
  vector<shared_ptr<Process>> procs_vec_1lep;
  vector<shared_ptr<Process>> procs_vec_2lep;

  for (const auto& pair : sample_info_met150) {
    const auto& searchstring = pair.first;
    const auto& ProcInfo = pair.second;
    std::set<string> fullpath;
    for (const auto& proc : searchstring){ 
      string path(ProcInfo.basepath + "*" + proc + "*");
      fullpath.insert(path);
      }
    //Other processes
      auto proc = Process::MakeShared<Baby_pico>(ProcInfo.proclegend, ProcInfo.type, TColor::GetColor(ProcInfo.proccolor.c_str()), {fullpath}, ProcInfo.special_proccuts && TTJ_stitch && WJets_stitch);
      procs_vec_met150.push_back(proc);
  }

  for (const auto& pair : sample_info_1lep) {
    const auto& searchstring = pair.first;
    const auto& ProcInfo = pair.second;
    std::set<string> fullpath;
    for (const auto& proc : searchstring){ 
      string path(ProcInfo.basepath + "*" + proc + "*");
      fullpath.insert(path);
      }
    //Other processes
      auto proc = Process::MakeShared<Baby_pico>(ProcInfo.proclegend, ProcInfo.type, TColor::GetColor(ProcInfo.proccolor.c_str()), {fullpath}, ProcInfo.special_proccuts && TTJ_stitch && WJets_stitch);
      procs_vec_1lep.push_back(proc);
  }

  for (const auto& pair : sample_info_2lep) {
    const auto& searchstring = pair.first;
    const auto& ProcInfo = pair.second;
    std::set<string> fullpath;
    for (const auto& proc : searchstring){ 
      string path(ProcInfo.basepath + "*" + proc + "*");
      fullpath.insert(path);
      }
    //Other processes
      auto proc = Process::MakeShared<Baby_pico>(ProcInfo.proclegend, ProcInfo.type, TColor::GetColor(ProcInfo.proccolor.c_str()), {fullpath}, ProcInfo.special_proccuts && TTJ_stitch && WJets_stitch);
      procs_vec_2lep.push_back(proc);
  }
 
 
  //Define plot options
  PlotOpt cmspaper("txt/plot_styles.txt","CMSPaper");
  cmspaper.Title(TitleType::info) //simulation, preliminary
          .Overflow(OverflowType::none)
          .FileExtensions({"pdf"})
          .Stack(StackType::data_norm); //other options: shapes, data_norm, signal_overlay
//  vector<PlotOpt> ops = {cmspaper.YAxis(YAxisType::linear), cmspaper.YAxis(YAxisType::log)};
  vector<PlotOpt> ops = {cmspaper.YAxis(YAxisType::linear)};
//  vector<PlotOpt> ops_log = {cmspaper.YAxis(YAxisType::log)};

  //Plots
  PlotMaker pm;
  
  pm.Push<Table>("higgsino4b_CR1_wrun2", vector<TableRow>{
    TableRow("Single lepton CR (electron)", CR1_el && single_ele_trig, 0, 0, w_lumi_gmsb),
    TableRow("Single lepton CR (muon)", CR1_mu && single_muon_trig, 0, 0, w_lumi_gmsb),
  }, procs_vec_1lep, false, true, false, false, false, false);

  pm.Push<Table>("higgsino4b_CR2_wrun2", vector<TableRow>{
    TableRow("Double lepton CR (electron)", CR2_el && single_ele_trig, 0, 0, w_lumi_gmsb),
    TableRow("Double lepton CR (muon)", CR2_mu && single_muon_trig, 0, 0, w_lumi_gmsb),
  }, procs_vec_2lep, false, true, false, false, false, false);

  pm.Push<Table>("higgsino4b_CR3resolved_wrun2", vector<TableRow>{
    TableRow("QCD CR", CR3 && ptmiss_htmiss_trig, 0, 0, w_lumi_gmsb_mettrig),
  }, procs_vec_met150, false, true, false, false, false, false);

//CR1 electron
  pm.Push<Hist1D>(Axis(24,0,240, hig_am, "<m_{bb}> [GeV]", {200}), CR1_el_nm_higam && single_ele_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_el_plot_hig_am_wrun2");
  pm.Push<Hist1D>(Axis(20,0,100, hig_dm, "#Deltam_{HH} [GeV]", {40}), CR1_el_nm_higdm && single_ele_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_el_plot_hig_dm_wrun2");
  pm.Push<Hist1D>(Axis(20,0,4, hig_drmax, "#DeltaR_{max}", {2.2}), CR1_el_nm_higdrmax && single_ele_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_el_plot_hig_drmax_wrun2");

  pm.Push<Hist1D>(Axis(20,0,400, jets_hh_pt, "Jet p_{T} [GeV]", {}), CR1_el_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_el_plot_jets_hh_pt_wrun2");
  pm.Push<Hist1D>(Axis(20,-3.14,3.14, jets_hh_eta, "Jet #eta", {}), CR1_el_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_el_plot_jets_hh_eta_wrun2");
  pm.Push<Hist1D>(Axis(20,-3.14,3.14, jets_hh_phi, "Jet #phi", {}), CR1_el_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_el_plot_jets_hh_phi_wrun2");
  pm.Push<Hist1D>(Axis(30,0,300, met, "p_{T}^{miss} [GeV]", {}), CR1_el_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_el_plot_met_wrun2");
  pm.Push<Hist1D>(Axis(20,0,200, lep1_pt, "Electron p_{T} [GeV]", {}), CR1_el_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_el_plot_el1pt_wrun2");
  pm.Push<Hist1D>(Axis(20,-2.5,2.5, lep1_eta, "Electron #eta", {}), CR1_el_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_el_plot_el1eta_wrun2");
  pm.Push<Hist1D>(Axis(20,-3.14,3.14, lep1_phi, "Electron #phi", {}), CR1_el_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_el_plot_el1phi_wrun2");
  pm.Push<Hist1D>(Axis(20,0,100,mt, "m_{T} [GeV]", {}), CR1_el_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_el_plot_mt_wrun2");

  pm.Push<Hist1D>(Axis(20,0,4, drjl1_min, "min. #DeltaR(jet, lepton)", {}), CR1_el_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_el_plot_mindrjl_wrun2");
  pm.Push<Hist1D>(Axis(25,0,5, dphi1, "#Delta#phi_{j_{1},MET}", {}), CR1_el_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_el_plot_dphi1_wrun2");
  pm.Push<Hist1D>(Axis(25,0,5, dphi2, "#Delta#phi_{j_{2},MET}", {}), CR1_el_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_el_plot_dphi2_wrun2");
  pm.Push<Hist1D>(Axis(25,0,5, dphi3, "#Delta#phi_{j_{3},MET}", {}), CR1_el_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_el_plot_dphi3_wrun2");
  pm.Push<Hist1D>(Axis(25,0,5, dphi4, "#Delta#phi_{j_{4},MET}", {}), CR1_el_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_el_plot_dphi4_wrun2");

// CR1 muon
  pm.Push<Hist1D>(Axis(24,0,240, hig_am, "<m_{bb}> [GeV]", {200}), CR1_mu_nm_higam && single_muon_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_mu_plot_hig_am_wrun2");
  pm.Push<Hist1D>(Axis(20,0,100, hig_dm, "#Deltam_{HH} [GeV]", {40}), CR1_mu_nm_higdm && single_muon_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_mu_plot_hig_dm_wrun2");
  pm.Push<Hist1D>(Axis(20,0,4, hig_drmax, "#DeltaR_{max}", {2.2}), CR1_mu_nm_higdrmax && single_muon_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_mu_plot_hig_drmax_wrun2");

  pm.Push<Hist1D>(Axis(20,0,400, jets_hh_pt, "Jet p_{T} [GeV]", {}), CR1_mu_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_mu_plot_jets_hh_pt_wrun2");
  pm.Push<Hist1D>(Axis(20,-3.14,3.14, jets_hh_eta, "Jet #eta", {}), CR1_mu_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_mu_plot_jets_hh_eta_wrun2");
  pm.Push<Hist1D>(Axis(20,-3.14,3.14, jets_hh_phi, "Jet #phi", {}), CR1_mu_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_mu_plot_jets_hh_phi_wrun2");
  pm.Push<Hist1D>(Axis(30,0,300, met, "p_{T}^{miss} [GeV]", {}), CR1_mu_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_mu_plot_met_wrun2");
  pm.Push<Hist1D>(Axis(20,0,200, lep1_pt, "Muon p_{T} [GeV]", {}), CR1_mu_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_mu_plot_mu1pt_wrun2");
  pm.Push<Hist1D>(Axis(20,-2.5,2.5, lep1_eta, "Muon #eta", {}), CR1_mu_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_mu_plot_mu1eta_wrun2");
  pm.Push<Hist1D>(Axis(20,-3.14,3.14, lep1_phi, "Muon #phi", {}), CR1_mu_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_mu_plot_mu1phi_wrun2");
  pm.Push<Hist1D>(Axis(20,0,100,mt, "m_{T} [GeV]", {}), CR1_mu_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_mu_plot_mt_wrun2");

  pm.Push<Hist1D>(Axis(20,0,4, drjl1_min, "min. #DeltaR(jet, lepton)", {}), CR1_mu_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_mu_plot_mindrjl_wrun2");
  pm.Push<Hist1D>(Axis(25,0,5, dphi1, "#Delta#phi_{j_{1},MET}", {}), CR1_mu_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_mu_plot_dphi1_wrun2");
  pm.Push<Hist1D>(Axis(25,0,5, dphi2, "#Delta#phi_{j_{2},MET}", {}), CR1_mu_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_mu_plot_dphi2_wrun2");
  pm.Push<Hist1D>(Axis(25,0,5, dphi3, "#Delta#phi_{j_{3},MET}", {}), CR1_mu_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_mu_plot_dphi3_wrun2");
  pm.Push<Hist1D>(Axis(25,0,5, dphi4, "#Delta#phi_{j_{4},MET}", {}), CR1_mu_trig, procs_vec_1lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR1_mu_plot_dphi4_wrun2");

// CR2 electron
  pm.Push<Hist1D>(Axis(24,0,240, hig_am, "<m_{bb}> [GeV]", {200}), CR2_el_nm_higam && single_ele_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_el_plot_hig_am_wrun2");
  pm.Push<Hist1D>(Axis(20,0,100, hig_dm, "#Deltam_{HH} [GeV]", {40}), CR2_el_nm_higdm && single_ele_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_el_plot_hig_dm_wrun2");
  pm.Push<Hist1D>(Axis(20,0,4, hig_drmax, "#DeltaR_{max}", {2.2}), CR2_el_nm_higdrmax && single_ele_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_el_plot_hig_drmax_wrun2");
  pm.Push<Hist1D>(Axis(20,0,400, jets_hh_pt, "Jet p_{T} [GeV]", {}), CR2_el_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_el_plot_jets_hh_pt_wrun2");
  pm.Push<Hist1D>(Axis(20,-3.14,3.14, jets_hh_eta, "Jet #eta",{}), CR2_el_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_el_plot_jets_hh_eta_wrun2");
  pm.Push<Hist1D>(Axis(20,-3.14,3.14, jets_hh_phi, "Jet #phi", {}), CR2_el_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_el_plot_jets_hh_phi_wrun2");
  pm.Push<Hist1D>(Axis(10,0,50, met, "p_{T}^{miss} [GeV]", {}), CR2_el_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_el_plot_met_wrun2");
  pm.Push<Hist1D>(Axis(20,0,200, lep_pt1, "Lead electron p_{T} [GeV]", {}), CR2_el_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_el_plot_leadelpt_wrun2");
  pm.Push<Hist1D>(Axis(20,-2.5,2.5, lep_eta1, "Lead electron #eta", {}), CR2_el_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_el_plot_leadeleta_wrun2");
  pm.Push<Hist1D>(Axis(20,-3.14,3.14, lep_phi1, "Lead electron #phi", {}), CR2_el_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_el_plot_leadelphi_wrun2");
  pm.Push<Hist1D>(Axis(20,0,200, lep_pt2, "Sublead electron p_{T} [GeV]", {}), CR2_el_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_el_plot_subelpt_wrun2");
  pm.Push<Hist1D>(Axis(20,-2.5,2.5, lep_eta2, "Sublead electron #eta", {}), CR2_el_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_el_plot_subeleta_wrun2");
  pm.Push<Hist1D>(Axis(20,-3.14,3.14, lep_phi2, "Sublead electron #phi", {}), CR2_el_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_el_plot_subelphi_wrun2");
  pm.Push<Hist1D>(Axis(20,0,400, ll_pt, "p_{T} (ee) [GeV]", {}), CR2_el_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_el_plot_elelpt_run2");
  pm.Push<Hist1D>(Axis(20,0,4, ll_dr, "#DeltaR(ee)", {}), CR2_el_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_el_plot_eleldr_wrun2");
  pm.Push<Hist1D>(Axis(20,0,3.14, ll_dphi, "#Delta#phi(ee)", {}), CR2_el_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_el_plot_eleldphi_wrun2");
  pm.Push<Hist1D>(Axis(20,0,4, ll_deta, "#Delta#eta(ee)", {}), CR2_el_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_el_plot_eleldeta_wrun2");
  pm.Push<Hist1D>(Axis(10,80,100, ll_m, "m(ee) [GeV]", {}), CR2_el_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_el_plot_elelm_wrun2");
  pm.Push<Hist1D>(Axis(20,0,4, drjl1_min, "min. #DeltaR(jet, lead lepton)", {}), CR2_el_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_el_plot_mindrjl1_wrun2");
  pm.Push<Hist1D>(Axis(20,0,4, drjl2_min, "min. #DeltaR(jet, sublead lepton)", {}), CR2_el_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_el_plot_mindrjl2_wrun2");
  pm.Push<Hist1D>(Axis(25,0,5, dphi1, "#Delta#phi_{j_{1},MET}", {}), CR2_el_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_el_plot_dphi1_wrun2");
  pm.Push<Hist1D>(Axis(25,0,5, dphi2, "#Delta#phi_{j_{2},MET}", {}), CR2_el_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_el_plot_dphi2_wrun2");
  pm.Push<Hist1D>(Axis(25,0,5, dphi3, "#Delta#phi_{j_{3},MET}", {}), CR2_el_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_el_plot_dphi3_wrun2");
  pm.Push<Hist1D>(Axis(25,0,5, dphi4, "#Delta#phi_{j_{4},MET}", {}), CR2_el_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_el_plot_dphi4_wrun2");

// CR2 muon
  pm.Push<Hist1D>(Axis(24,0,240, hig_am, "<m_{bb}> [GeV]", {200}), CR2_mu_nm_higam && single_muon_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_mu_plot_hig_am_wrun2");
  pm.Push<Hist1D>(Axis(20,0,100, hig_dm, "#Deltam_{HH} [GeV]", {40}), CR2_mu_nm_higdm && single_muon_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_mu_plot_hig_dm_wrun2");
  pm.Push<Hist1D>(Axis(20,0,4, hig_drmax, "#DeltaR_{max}", {2.2}), CR2_mu_nm_higdrmax && single_muon_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_mu_plot_hig_drmax_wrun2");
  pm.Push<Hist1D>(Axis(20,0,400, jets_hh_pt, "Jet p_{T} [GeV]", {}), CR2_mu_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_mu_plot_jets_hh_pt_wrun2");
  pm.Push<Hist1D>(Axis(20,-3.14,3.14, jets_hh_eta, "Jet #eta", {}), CR2_mu_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_mu_plot_jets_hh_eta_wrun2");
  pm.Push<Hist1D>(Axis(20,-3.14,3.14, jets_hh_phi, "Jet #phi", {}), CR2_mu_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_mu_plot_jets_hh_phi_wrun2");
  pm.Push<Hist1D>(Axis(10,0,50, met, "p_{T}^{miss} [GeV]", {}), CR2_mu_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_mu_plot_met_wrun2");
  pm.Push<Hist1D>(Axis(20,0,200, lep_pt1, "Lead muon p_{T} [GeV]", {}), CR2_mu_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_mu_plot_leadelpt_wrun2");
  pm.Push<Hist1D>(Axis(20,-2.5,2.5, lep_eta1, "Lead muon #eta", {}), CR2_mu_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_mu_plot_leadeleta_wrun2");
  pm.Push<Hist1D>(Axis(20,-3.14,3.14, lep_phi1, "Lead muon #phi", {}), CR2_mu_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_mu_trig_plot_leadelphi_wrun2");
  pm.Push<Hist1D>(Axis(20,0,200, lep_pt2, "Sublead muon p_{T} [GeV]", {}), CR2_mu_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_mu_trig_plot_subelpt_wrun2");
  pm.Push<Hist1D>(Axis(20,-2.5,2.5, lep_eta2, "Sublead muon #eta", {}), CR2_mu_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_mu_plot_subeleta_wrun2");
  pm.Push<Hist1D>(Axis(20,-3.14,3.14, lep_phi2, "Sublead muon #phi", {}), CR2_mu_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_mu_plot_subelphi_wrun2");
  pm.Push<Hist1D>(Axis(20,0,400, ll_pt, "p_{T} (#mu#mu) [GeV]", {}), CR2_mu_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_mu_plot_mumupt_run2");
  pm.Push<Hist1D>(Axis(20,0,4, ll_dr, "#DeltaR(#mu#mu)", {}), CR2_mu_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_mu_plot_mumudr_wrun2");
  pm.Push<Hist1D>(Axis(20,0,3.14, ll_dphi, "#Delta#phi(#mu#mu)", {}), CR2_mu_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_mu_plot_mumudphi_wrun2");
  pm.Push<Hist1D>(Axis(20,0,4, ll_deta, "#Delta#eta(#mu#mu)", {}), CR2_mu_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_mu_plot_mumudeta_wrun2");
  pm.Push<Hist1D>(Axis(10,80,100, ll_m, "m(#mu#mu) [GeV]", {}), CR2_mu_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_mu_plot_mumum_wrun2");
  pm.Push<Hist1D>(Axis(20,0,4, drjl1_min, "min. #DeltaR(jet, lead lepton)", {}), CR2_mu_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_mu_plot_mindrjl1_wrun2");
  pm.Push<Hist1D>(Axis(20,0,4, drjl2_min, "min. #DeltaR(jet, sublead lepton)", {}), CR2_mu_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_mu_plot_mindrjl2_wrun2");
  pm.Push<Hist1D>(Axis(25,0,5, dphi1, "#Delta#phi_{j_{1},MET}", {}), CR2_mu_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_mu_plot_dphi1_wrun2");
  pm.Push<Hist1D>(Axis(25,0,5, dphi2, "#Delta#phi_{j_{2},MET}", {}), CR2_mu_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_mu_plot_dphi2_wrun2");
  pm.Push<Hist1D>(Axis(25,0,5, dphi3, "#Delta#phi_{j_{3},MET}", {}), CR2_mu_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_mu_plot_dphi3_wrun2");
  pm.Push<Hist1D>(Axis(25,0,5, dphi4, "#Delta#phi_{j_{4},MET}", {}), CR2_mu_trig, procs_vec_2lep, ops).Weight(w_lumi_gmsb).Tag("ShortName:CR2_mu_plot_dphi4_wrun2");

// CR3
  pm.Push<Hist1D>(Axis(30,0,600, hig_am, "<m_{bb}> [GeV]", {200}), CR3_nm_higam && ptmiss_htmiss_trig, procs_vec_met150, ops).Weight(w_lumi_gmsb_mettrig).Tag("ShortName:CR3_plot_hig_am_wrun2");
  pm.Push<Hist1D>(Axis(30,0,300, hig_dm, "#Deltam_{HH} [GeV]", {40}), CR3_nm_higdm && ptmiss_htmiss_trig, procs_vec_met150, ops).Weight(w_lumi_gmsb_mettrig).Tag("ShortName:CR3_plot_hig_dm_wrun2");
  pm.Push<Hist1D>(Axis(20,0,4, hig_drmax, "#DeltaR_{max}", {2.2}), CR3_nm_higdrmax && ptmiss_htmiss_trig, procs_vec_met150, ops).Weight(w_lumi_gmsb_mettrig).Tag("ShortName:CR3_plot_hig_drmax_wrun2");
  pm.Push<Hist1D>(Axis(20,0,400, jets_hh_pt, "Jet p_{T} [GeV]", {}), CR3_met_trig, procs_vec_met150, ops).Weight(w_lumi_gmsb_mettrig).Tag("ShortName:CR3_plot_jets_hh_pt_wrun2");
  pm.Push<Hist1D>(Axis(20,-3.14,3.14, jets_hh_eta, "Jet #eta", {}), CR3_met_trig, procs_vec_met150, ops).Weight(w_lumi_gmsb_mettrig).Tag("ShortName:CR3_plot_jets_hh_eta_wrun2");
  pm.Push<Hist1D>(Axis(20,-3.14,3.14, jets_hh_phi, "Jet #phi", {}), CR3_met_trig, procs_vec_met150, ops).Weight(w_lumi_gmsb_mettrig).Tag("ShortName:CR3_plot_jets_hh_phi_wrun2");
  pm.Push<Hist1D>(Axis(30,200,500, met, "p_{T}^{miss} [GeV]", {}), CR3_met_trig, procs_vec_met150, ops).Weight(w_lumi_gmsb_mettrig).Tag("ShortName:CR3_plot_met_wrun2");

  pm.Push<Hist1D>(Axis(25,0,5, dphi1, "#Delta#phi_{j_{1},MET}", {}), CR3_met_trig, procs_vec_met150, ops).Weight(w_lumi_gmsb_mettrig).Tag("ShortName:CR3_plot_dphi1_wrun2");
  pm.Push<Hist1D>(Axis(25,0,5, dphi2, "#Delta#phi_{j_{2},MET}", {}), CR3_met_trig, procs_vec_met150, ops).Weight(w_lumi_gmsb_mettrig).Tag("ShortName:CR3_plot_dphi2_wrun2");
  pm.Push<Hist1D>(Axis(25,0,5, dphi3, "#Delta#phi_{j_{3},MET}", {}), CR3_met_trig, procs_vec_met150, ops).Weight(w_lumi_gmsb_mettrig).Tag("ShortName:CR3_plot_dphi3_wrun2");
  pm.Push<Hist1D>(Axis(25,0,5, dphi4, "#Delta#phi_{j_{4},MET}", {}), CR3_met_trig, procs_vec_met150, ops).Weight(w_lumi_gmsb_mettrig).Tag("ShortName:CR3_plot_dphi4_wrun2");


//  pm.multithreaded_ = true;
  pm.min_print_ = true;
  pm.MakePlots(1);
  
}






