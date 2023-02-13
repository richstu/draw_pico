//This script generates some tables and plots for gg->H->Zg studies

#include "core/test.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>
#include <getopt.h>

#include "TColor.h"
#include "TError.h"

#include "core/baby.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/named_func.hpp"
#include "core/named_func_utilities.hpp"
#include "core/palette.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/process.hpp"
#include "core/table.hpp"
#include "core/utilities.hpp"
#include "zgamma/apply_zg_trigeffs.hpp"
#include "zgamma/zg_functions.hpp"
#include "zgamma/zg_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace ZgFunctions;
using namespace ZgUtilities;

int main() {

  //--------------------------------------------------------------------------
  //                          define NamedFuncs
  //--------------------------------------------------------------------------

  NamedFunc w_sigx20("w_sigx20",[](const Baby &b) -> NamedFunc::ScalarType{
    if(b.type() >= 200000 && b.type() <= 205000)
      return 20.0;
    return 1.0;
  });

  NamedFunc el1_pt = "el_pt[ll_i1[0]]";
  NamedFunc el2_pt = "el_pt[ll_i2[0]]";
  NamedFunc mu1_pt = "mu_pt[ll_i1[0]]";
  NamedFunc mu2_pt = "mu_pt[ll_i2[0]]";
  NamedFunc el1_eta = "el_eta[ll_i1[0]]";
  NamedFunc el2_eta = "el_eta[ll_i2[0]]";
  NamedFunc mu1_eta = "mu_eta[ll_i1[0]]";
  NamedFunc mu2_eta = "mu_eta[ll_i2[0]]";
  NamedFunc photon_pt = "photon_pt[0]";
  NamedFunc llphoton_cosTheta = "llphoton_cosTheta[0]";
  llphoton_cosTheta.Name("llphoton_cosCapTheta0");

  NamedFunc lep1_eta("lep1_eta",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.ll_lepid()->at(0)==11) return b.el_eta()->at(b.ll_i1()->at(0));
    return b.mu_eta()->at(b.ll_i1()->at(0));
  });

  NamedFunc lep2_eta("lep2_eta",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.ll_lepid()->at(0)==11) return b.el_eta()->at(b.ll_i2()->at(0));
    return b.mu_eta()->at(b.ll_i2()->at(0));
  });

  NamedFunc zg_el_cuts = "(ll_lepid[0]==11) && (el_pt[ll_i1[0]]>25) && (el_pt[ll_i2[0]]>15)";
  zg_el_cuts.Name("zg_el_cuts");
  NamedFunc zg_mu_cuts = "(ll_lepid[0]==13) && (mu_pt[ll_i1[0]]>20) && (mu_pt[ll_i2[0]]>10)";
  zg_mu_cuts.Name("zg_mu_cuts");
  NamedFunc higgs_window = "llphoton_m[0]>120&&llphoton_m[0]<130";

  //--------------------------------------------------------------------------
  //                          setup
  //--------------------------------------------------------------------------
  
  gErrorIgnoreLevel = 6000;
  Palette mc_colors("txt/colors_zgamma.txt","default");
  Process::Type back =  Process::Type::background;
  Process::Type sig =  Process::Type::signal;

  bool use_local = true;
  string prod_folder("/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v2/");
  std::set<int> years = {2016,2017,2018};
  std::string lumi_string = "138";
  if (use_local) {
    prod_folder = "/data"+prod_folder;
  }

  bool multiply_sig = true;
  std::string sig_name = "gg#rightarrow H#rightarrow Z#gamma";
  NamedFunc weight = "w_lumi"*w_years;
  if (multiply_sig) {
    sig_name = "gg#rightarrow H#rightarrow Z#gamma (x20)";
    weight = "w_lumi"*w_sigx20*w_years;
  }

  std::cout << "Loading samples:\n";
  auto proc_smzg      = Process::MakeShared<Baby_pico>("Z+#gamma", back, mc_colors("zgtollg"),
                           attach_folder(prod_folder,years,"mc/merged_zgmc_llg",{"*ZGToLLG*"}), 
                           HLT_pass_dilepton);
  auto proc_dy        = Process::MakeShared<Baby_pico>("Z+Fake Photon", back, mc_colors("dyjets"),
                           attach_folder(prod_folder,years,"mc/merged_zgmc_llg",{"*DYJets*"}), 
                           "stitch_dy"&&HLT_pass_dilepton);
  auto proc_tt        = Process::MakeShared<Baby_pico>("tt", back, mc_colors("tt"),
                           attach_folder(prod_folder,years,"mc/merged_zgmc_llg",{"*TTTo2L2Nu*"}), 
                           HLT_pass_dilepton);
  auto proc_tt1l      = Process::MakeShared<Baby_pico>("Fake lepton", back, mc_colors("fakelep"),
                           attach_folder(prod_folder,years,"mc/merged_zgmc_llg",{"*TTG*"}), 
                           "(ntruel+ntrumu+ntrutaul)==1"&&HLT_pass_dilepton);
  //auto proc_hzg       = Process::MakeShared<Baby_pico>(sig_name, sig, kRed,
  //                         attach_folder(prod_folder,years,"mc/merged_zgmc_llg",
  //                         {"*HToZG*M-125*","*HToZG_M125*"}), HLT_pass_dilepton);
  auto proc_hzg       = Process::MakeShared<Baby_pico>(sig_name, sig, kRed,
                           attach_folder(prod_folder,years,"mc/merged_zgmc_llg",
                           {"*GluGluHToZG*M-125*"}), HLT_pass_dilepton);
  auto proc_hzg_bak   = Process::MakeShared<Baby_pico>(sig_name, back, kRed,
                           attach_folder(prod_folder,years,"mc/merged_zgmc_llg",
                           {"*GluGluHToZG*M-125*"}), HLT_pass_dilepton);

  vector<shared_ptr<Process>> procs = {proc_dy, proc_tt, proc_smzg, proc_tt1l, proc_hzg};
  vector<shared_ptr<Process>> procs_bak = {proc_dy, proc_tt, proc_smzg, proc_tt1l};
  vector<shared_ptr<Process>> procs_sig = {proc_hzg_bak};

  PlotOpt lin_lumi("txt/plot_styles.txt","CMSPaper");
  lin_lumi.Title(TitleType::info)
          .Stack(StackType::signal_overlay)
          .Overflow(OverflowType::none)
          .Bottom(BottomType::off)
          .UseCMYK(false)
          .LegendColumns(1)
          .ShowBackgroundError(false)
          .FileExtensions({"pdf"});
  PlotOpt lin_sorb = lin_lumi;
  lin_sorb.Overflow(OverflowType::both).Bottom(BottomType::sorb).FileExtensions({"pdf","root"});
  PlotOpt lin_sorb_upper = lin_lumi;
  lin_sorb_upper.Overflow(OverflowType::both).Bottom(BottomType::sorb_cut_upper).FileExtensions({"pdf","root"});
  PlotOpt lin_shapes = lin_lumi;
  lin_shapes.Stack(StackType::shapes);
  PlotOpt log_lumi2d("txt/plot_styles.txt", "Scatter");
  log_lumi2d.Title(TitleType::info)
                 .YAxis(YAxisType::log)
                 .Overflow(OverflowType::overflow)
                 .LogMinimum(0.001);
  std::vector<PlotOpt> ops = {lin_lumi, lin_shapes};
  std::vector<PlotOpt> ops_sorb = {lin_shapes,lin_sorb};
  std::vector<PlotOpt> ops_sorb_upper = {lin_shapes,lin_sorb_upper};
  std::vector<PlotOpt> ops2d = {log_lumi2d};

  //--------------------------------------------------------------------------
  //                         Make plots & tables
  //--------------------------------------------------------------------------

  PlotMaker pm;

  bool plot_pts = false;
  bool plot_baseline = false;
  bool plot_bdtvars = true;

  if (plot_pts) {
    pm.Push<Hist1D>(Axis(40,0,80, el1_pt, "Lead Electron p_{T} [GeV]", {25}), 
        "ll_lepid[0]==11", procs, ops).Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(40,0,80, el2_pt, "Sublead Electron p_{T} [GeV]", {15}), 
        "ll_lepid[0]==11", procs, ops).Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(40,0,80, mu1_pt, "Lead Muon p_{T} [GeV]", {20}), 
        "ll_lepid[0]==13", procs, ops).Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(40,0,80, mu2_pt, "Sublead Muon p_{T} [GeV]", {10}), 
        "ll_lepid[0]==13", procs, ops).Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
  }

  if (plot_baseline) {
    NamedFunc baseline_nomll = (zg_el_cuts||zg_mu_cuts)&&
        "(photon_pt[0]/llphoton_m[0]>=15.0/110.0)&&photon_drmin[0]>0.4&&photon_idmva[0]>0.5";
    baseline_nomll.Name("baseline_nomll");
    NamedFunc baseline_noptg = (zg_el_cuts||zg_mu_cuts)&&
        "(ll_m[0]>50) && ((llphoton_m[0]+ll_m[0])>=185) && photon_drmin[0]>0.4 && photon_idmva[0]>0.5";
    baseline_noptg.Name("baseline_noptg");
    NamedFunc ptg_over_mllg = "photon_pt[0]/llphoton_m[0]";
    //mll plots
    pm.Push<Hist2D>(Axis(20, 40, 120, "ll_m[0]", "m_{ll} [GeV]", {}),
        Axis(20, 80, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
        baseline_nomll, procs_bak, ops2d).Weight(weight).Tag("zgggh_bak").LuminosityTag(lumi_string);
    pm.Push<Hist2D>(Axis(20, 40, 120, "ll_m[0]", "m_{ll} [GeV]", {}),
        Axis(20, 80, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
        baseline_nomll, procs_sig, ops2d).Weight(weight).Tag("zgggh_sig").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30,100,160, "llphoton_m[0]", "m_{ll#gamma} [GeV]"), 
        baseline_nomll, procs, ops).Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30,100,160, "llphoton_m[0]", "m_{ll#gamma} [GeV]"), 
        zg_baseline&&"photon_idmva[0]>0.5", procs, ops).Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(20,50,130, "ll_m[0]", "m_{ll} [GeV]"), 
        baseline_nomll&&higgs_window, procs, ops).Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    //photon pt plots
    pm.Push<Hist1D>(Axis(26,15,80, "photon_pt[0]", "Photon p_{T} [GeV]"), 
        baseline_noptg&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(26,0.115,0.667, ptg_over_mllg, "Photon p_{T}/m_{ll#gamma}"), 
        baseline_noptg&&higgs_window, procs, ops_sorb)
        .Weight("w_lumi"*w_years).Tag("zgggh").LuminosityTag(lumi_string);
  }

  if (plot_bdtvars) {
    NamedFunc photon_res = "photon_pterr[0]/photon_pt[0]";
    NamedFunc llphoton_pt_mass = "llphoton_pt[0]/llphoton_m[0]";
    NamedFunc ll_pt_mass = "ll_pt[0]/llphoton_m[0]";
    //
    pm.Push<Hist1D>(Axis(30,-1.0,1.0, "photon_idmva[0]", "Photon IDMVA"), 
        zg_baseline&&higgs_window, procs, ops_sorb)
        .Weight("w_lumi"*w_years).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30,0.0,0.5, photon_res, "Photon p_{T} relative uncertainty"), 
        zg_baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    //
    pm.Push<Hist1D>(Axis(30,0.4,3.5, "photon_drmin[0]", "Photon minimum #Delta R(l,#gamma)"), 
        zg_baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30,0.4,5.0, photon_drmax, "Photon maximum #Delta R(l,#gamma)"), 
        zg_baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30,0.0,4.0, llphoton_pt_mass, "p_{Tll#gamma}/m_{ll#gamma}"), 
        zg_baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    //
    pm.Push<Hist1D>(Axis(30,-1.0,1.0, llphoton_cosTheta, "cos #Theta"), 
        zg_baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30,-2.5,2.5, "photon_eta[0]", "Photon #eta"), 
        zg_baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30,-2.5,2.5, lep1_eta, "Lead lepton #eta"), 
        zg_baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30,-2.5,2.5, lep2_eta, "Sublead lepton #eta"), 
        zg_baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    //
    pm.Push<Hist1D>(Axis(30,-1.0,1.0, "llphoton_costheta[0]", "cos #theta"), 
        zg_baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30,-3.14,3.14, "llphoton_psi[0]", "#phi"), 
        zg_baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    //
    pm.Push<Hist1D>(Axis(30,0.0,0.667, ll_pt_mass, "p_{Tll}/m_{ll#gamma}"), 
        zg_baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30,0.0,5.0, "llphoton_dr[0]", "#Delta R(ll,#gamma)"), 
        zg_baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    //
    pm.Push<Hist1D>(Axis(30,0.0,80.0, "met", "p_{T}^{miss} [GeV]"), 
        zg_baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30,0.0,80.0, "met", "p_{T}^{miss} [GeV]"), 
        zg_baseline&&higgs_window, procs, ops)
        .Weight(weight).Tag("zgggh").LuminosityTag(lumi_string);
  }

  pm.min_print_ = true;
  pm.MakePlots(1.0);

  return 0;
}
