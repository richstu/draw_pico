/**
 * Script to write datacard for H->Zgamma analysis
 */

#include <cmath>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "TError.h"
#include "TColor.h"

#include "core/axis.hpp"
#include "core/baby_pico.hpp"
#include "core/fastforest.hpp"
#include "core/hist1d.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/process.hpp"
#include "core/utilities.hpp"
#include "zgamma/eventweighter.hpp"
#include "zgamma/scalesmear.hpp"
#include "zgamma/zg_functions.hpp"
#include "zgamma/zg_syst_functions.hpp"
#include "zgamma/zg_utilities.hpp"

using std::cout;
using std::endl;
using std::make_shared;
using std::set;
using std::shared_ptr;
using std::string;
using std::vector;
using fastforest::FastForest;
using ZgUtilities::ZgSampleLoader;
using ZgUtilities::XGBoostBDTs;
using ZgUtilities::VBFXGBoostBDTs;
using ZgUtilities::xgb_ggf_bdt_250923_offsets;
using ZgUtilities::xgb_vbf_bdt_250923_offsets;
using ZgUtilities::XGBoostBDTScoreCached;
using ZgUtilities::category_ggf4;
using ZgUtilities::category_ggf3;
using ZgUtilities::category_ggf2;
using ZgUtilities::category_ggf1;
using ZgUtilities::VbfBdts;
using ZgUtilities::vbf_bdt_score;
using ZgUtilities::category_vbf4;
using ZgUtilities::category_vbf3;
using ZgUtilities::category_vbf2;
using ZgUtilities::category_vbf1;
using ZgUtilities::rename_signal;
using ZgUtilities::get_w_sigscale;
using namespace ZgFunctions;
//const Process::Type data = Process::Type::data;
const Process::Type signal = Process::Type::signal;
const Process::Type background = Process::Type::background;

int main() {

  //setup
  gErrorIgnoreLevel = 6000;

  //Define processes
  vector<shared_ptr<Process>> procsx300 = ZgSampleLoader()
        //.SetMacro("YEARS",{"2023BPix"}) //debug
        .LoadSamples("txt/samples_zgamma.txt","All");
  vector<shared_ptr<Process>> procsx1000 = ZgSampleLoader()
        //.SetMacro("YEARS",{"2023BPix"}) //debug
        .LoadSamples("txt/samples_zgamma.txt","All");
  vector<shared_ptr<Process>> procs_noscale = ZgSampleLoader()
        //.SetMacro("YEARS",{"2023BPix"}) //debug
        .LoadSamples("txt/samples_zgamma.txt","All");
  vector<shared_ptr<Process>> procs_vbf = ZgSampleLoader() 
        //.SetMacro("YEARS",{"2023BPix"}) //debug
        .LoadSamples("txt/samples_zgamma.txt","AllSplitVBF");
  rename_signal(procsx300, 300);
  rename_signal(procsx1000, 1000);
  rename_signal(procs_vbf, 100);
  vector<PlotOpt> ops = {PlotOpt("txt/plot_styles.txt","LinLumiDataRoot")
                         .RatioMinimum(0.5).RatioMaximum(1.5)
                         .Title(PlotOptTypes::TitleType::supplementary)
                         .Stack(PlotOptTypes::StackType::data_norm)};

  //Define NamedFuncs
  const NamedFunc lead_lepton_pt("lead_lepton_pt",[](const Baby &b) 
      -> NamedFunc::ScalarType{
    if (b.ll_lepid()->at(0)==11) {
      return b.el_pt()->at(0);
    }
    return b.mu_pt()->at(0);
  });

  const NamedFunc sublead_lepton_pt("sublead_lepton_pt",[](const Baby &b) 
      -> NamedFunc::ScalarType{
    if (b.ll_lepid()->at(0)==11) {
      return b.el_pt()->at(1);
    }
    return b.mu_pt()->at(1);
  });

  NamedFunc untagged_category_cached = NamedFunc(untagged_category)
      .EnableCaching(true);
  cs_setter();

  NamedFunc w_sigx1000 = get_w_sigscale(1000);
  w_sigx1000.Name("w_sigx1000");
  NamedFunc w_sigx300 = get_w_sigscale(300);
  w_sigx300.Name("w_sigx300");
  NamedFunc w_sigx100 = get_w_sigscale(100);
  w_sigx100.Name("w_sigx100");

  const vector<FastForest> ggf_bdts = XGBoostBDTs();
  const vector<FastForest> vbf_bdts = VBFXGBoostBDTs();
  const NamedFunc ggf_score_default = XGBoostBDTScoreCached(ggf_bdts, 
      xgb_ggf_bdt_250923_offsets, ggf_bdt_inputs_default, "ggf_default");
  const NamedFunc vbf_score_default = XGBoostBDTScoreCached(vbf_bdts, 
      xgb_vbf_bdt_250923_offsets, vbf_bdt_inputs_default, "vbf_default");

  //weight with some regularization
  const NamedFunc weight_reg(
      "weight_reg",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleTypeString().Contains("-")) 
      return 1.; //data
    //fix for WplusHto2Mu in redwood_v1
    float w = b.weight();
    if (isnan(b.weight()) || isinf(b.weight())) {
      w = b.w_lumi()*b.w_lep()*b.w_bhig_df()*b.w_trig()*b.w_isr()*b.w_pu()
          *b.w_prefire()*b.w_photon()*b.w_phshape()*b.w_fakephoton()
          *b.w_nnlo();
    }
    if (abs(b.SampleType())>2020) {
      if (b.type() == 12500)
        return w*0.000218;
    }
    if (fabs(w/b.w_lumi()) > 2) return b.w_lumi()*2.0;
    return w;
  });

  //Define weight
  NamedFunc weight(w_years*weight_reg*w_el); 

  NamedFunc baseline = NamedFunc("nphoton>=1&&nll>=1"&&sys_trig_pt_default&&"ll_m[0]>80&&ll_m[0]<100&&(photon_pt[0]/llphoton_m[0])>=15.0/110.0&&(llphoton_m[0]+ll_m[0])>=185&&pass").Name("baseline").EnableCaching(true);
  NamedFunc higgs_mass_region = NamedFunc("llphoton_refit_m>120&&llphoton_refit_m<130").Name("higgs_mass_region");
  NamedFunc ggf = NamedFunc("nlep==2"&&sys_met_default<90&&sys_njet_default<2).Name("ggF").EnableCaching(true);
  NamedFunc vbf = NamedFunc("nlep==2"&&sys_nbdfm_default==0.0&&sys_njet_default>=2).Name("VBF").EnableCaching(true);
  NamedFunc ggf4 = NamedFunc(ggf&&category_ggf4(ggf_score_default)).Name("ggf4").EnableCaching(true);
  NamedFunc ggf3 = NamedFunc(ggf&&category_ggf3(ggf_score_default)).Name("ggf3").EnableCaching(true);
  NamedFunc ggf2 = NamedFunc(ggf&&category_ggf2(ggf_score_default)).Name("ggf2").EnableCaching(true);
  NamedFunc ggf1 = NamedFunc(ggf&&category_ggf1(ggf_score_default)).Name("ggf1").EnableCaching(true);
  NamedFunc vbf4 = NamedFunc(vbf&&category_vbf4(vbf_score_default)).Name("vbf4").EnableCaching(true);
  NamedFunc vbf3 = NamedFunc(vbf&&category_vbf3(vbf_score_default)).Name("vbf3").EnableCaching(true);
  NamedFunc vbf2 = NamedFunc(vbf&&category_vbf2(vbf_score_default)).Name("vbf2").EnableCaching(true);
  NamedFunc vbf1 = NamedFunc(vbf&&category_vbf1(vbf_score_default)).Name("vbf1").EnableCaching(true);
  NamedFunc vhmet = NamedFunc("nlep==2"&&sys_met_default>90&&sys_njet_default<2&&"llphoton_pt[0]/llphoton_m[0]>0.4").Name("vhmet").EnableCaching(true);
  NamedFunc vh3l = NamedFunc("nlep>=3"&&sys_met_default>30&&sys_nbdfm_default==0.0&&max_lep_miniso<0.15&&"llphoton_pt[0]/llphoton_m[0]>0.3").Name("vh3l").EnableCaching(true);
  NamedFunc tthhad = NamedFunc("nlep==2"&&sys_nbdfm_default>=1&&sys_njet_default>=5&&"ll_m[0]>85&&ll_m[0]<95").Name("tthhad").EnableCaching(true);
  NamedFunc tthlep = NamedFunc((("nlep==3"&&sys_nbdfm_default>=1&&sys_njet_default>=3)||("nlep>=4"&&sys_nbdfm_default>=1&&sys_njet_default>=1))&&max_lep_miniso<0.1).Name("tthlep").EnableCaching(true);
  vector<NamedFunc> categories = {ggf1, ggf2, ggf3, ggf4, vbf1, vbf2, vbf3, vbf4, vh3l, vhmet, tthlep, tthhad, untagged_category_cached};

  //Make plots
  PlotMaker pm;
  //pm.multithreaded_ = false;
  pm.min_print_ = true;
  pm.max_threads_ = 16;

  vector<Axis> baseline_plots = {
      Axis(40,0.0,120.0, lead_lepton_pt, "Lead lepton p_{T} [GeV]", {}), 
      Axis(40,0.0,100.0, sublead_lepton_pt, "Sublead lepton p_{T} [GeV]", {}), 
      Axis(40,0.0,100.0, "photon_pt[0]", "Lead photon p_{T} [GeV]", {}), 
      Axis(40,0.0,100.0, "llphoton_pt[0]", "p_{T}(ll#gamma) [GeV]", {})};

  vector<Axis> ggf_plots = {
      Axis(40,0.0,1.0, sys_lead_photon_idmva_default, "Photon IDMVA", {}), 
      Axis(35,0.0,3.5, sys_lead_photon_drmin_default, "#Delta R_{min}(#gamma, l)", {}), 
      Axis(40,-1.0,1.0, sys_llphoton_cosTheta_default, "cos(#Theta)", {}), 
      Axis(40,0.0,2.5, sys_llphoton_relpt_default, "p_{T}(ll#gamma)/m_{ll#gamma}", {}), 
      Axis(40,0.0,6.0, sys_photon_jet1_dr_default, "#Delta R(#gamma, j)", {}), 
      Axis(40,0.0,0.15, sys_lead_photon_relpterr_default, "#sigma_{E}(#gamma)/E(#gamma)", {}), 
      Axis(50,0.0,5.0, sys_lead_photon_drmax_default, "#Delta R_{max}(#gamma, l)", {}), 
      Axis(40,0.0,3.1416, sys_photon_mht_dphi_default, "#Delta #phi (#gamma, S_{T}^{miss})", {}), 
      Axis(40,0.0,1.0, sys_llphoton_jet_balance_default, "System balance", {}), 
      Axis(40,-2.5,2.5, sys_lead_photon_eta_default, "Photon #eta", {}), 
      Axis(40,-4.7,4.7, sys_lead_jet_eta_default, "Jet #eta", {}), 
      Axis(40,-2.5,2.5, sys_lead_lepton_eta_default, "Lead lepton #eta", {}), 
      Axis(40,0.0,6.0, sys_photon_zeppenfeld_default, "|#eta(#gamma)-#eta(j)|", {}), 
      Axis(40,30.0,150.0, sys_lead_jet_pt_default, "Jet p_{T} [GeV]", {}), 
      Axis(40,0.0,40.0, sys_lead_jet_m_default, "Jet m [GeV]", {}), 
      Axis(40,0.0,3.1416, sys_llphoton_jet_dphi_default, "#Delta #phi(ll#gamma, j)", {}), 
      Axis(40,-1.0,1.0, sys_llphoton_costheta_default, "cos(#theta)", {}), 
      Axis(40,-2.5,2.5, sys_sublead_lepton_eta_default, "Sublead lepton #eta", {}), 
      Axis(40,-3.1416,3.1416, sys_llphoton_psi_default, "#phi", {}), 
      Axis(2,-0.5,1.5, sys_njet_default, "N_{j}", {})
  };

  vector<Axis> vbf_plots = {
      Axis(20,0.0,2.5, sys_llphoton_relpt_default, "p_{T}(ll#gamma)/m_{ll#gamma}", {}), 
      Axis(25,0.0,5.0, sys_lead_photon_drmax_default, "#Delta R_{max}(#gamma, l)", {}), 
      Axis(20,0.0,4.0, sys_lead_photon_drmin_default, "#Delta R_{min}(#gamma, l)", {}), 
      Axis(20,0.0,3.1416, sys_llphoton_dijet_dphi_default, "#Delta #phi(ll#gamma, jj)", {}), 
      Axis(20,0.0,0.15, sys_lead_photon_relpterr_default, "#sigma_{E}(#gamma)/E(#gamma)", {}), 
      Axis(20,0.0,9.0, sys_dijet_deta_default, "#Delta#eta(j_{1},j_{2})", {}), 
      Axis(20,0.0,1.0, sys_lead_photon_idmva_default, "Photon IDMVA", {}), 
      Axis(20,0.0,6.0, sys_photon_jet1_dr_default, "#Delta R(#gamma, j_{1})", {}), 
      Axis(20,0.0,6.0, sys_photon_jet2_dr_default, "#Delta R(#gamma, j_{2})", {}), 
      Axis(20,0.0,600.0, sys_dijet_m_default, "m_{jj} [GeV]", {}), 
      Axis(20,0.0,3.1416, sys_dijet_dphi_default, "#Delta #phi(j_{1},j_{2})", {}), 
      Axis(20,0.0,1.0, sys_llphoton_jet_balance_default, "System balance", {}), 
      Axis(20,-1.0,1.0, sys_llphoton_costheta_default, "cos(#theta)", {}), 
      Axis(20,30.0,150.0, sys_sublead_jet_pt_default, "Sublead jet p_{T} [GeV]", {}), 
      Axis(20,0.0,6.0, sys_photon_zeppenfeld_default, "Photon Zeppenfeld variable", {}), 
      Axis(20,30.0,200.0, sys_lead_jet_pt_default, "Lead jet p_{T} [GeV]", {}), 
      Axis(20,-1.0,1.0, sys_llphoton_cosTheta_default, "cos(#Theta)", {}), 
      Axis(20,0.0,3.1416, sys_photon_mht_dphi_default, "#Delta #phi (#gamma, S_{T}^{miss})", {}), 
      Axis(20,-2.5,2.5, sys_lead_photon_eta_default, "Photon #eta", {}), 
      Axis(20,-2.5,2.5, sys_lead_lepton_eta_default, "Lead lepton #eta", {}), 
      Axis(20,-2.5,2.5, sys_sublead_lepton_eta_default, "Sublead lepton #eta", {}), 
      Axis(3,1.5,4.5, sys_njet_default, "N_{j}", {}), 
      Axis(20,-3.1416,3.1416, sys_llphoton_psi_default, "#phi", {})
  };

  vector<Axis> mllg_plots = {
      Axis(65,95,160, "llphoton_refit_m", "m_{ll#gamma} [GeV]", {}),
      Axis(65,99,164, "llphoton_refit_m", "m_{ll#gamma} [GeV]", {}),
      Axis(65,105,170, "llphoton_refit_m", "m_{ll#gamma} [GeV]", {}),
      Axis(65,106,171, "llphoton_refit_m", "m_{ll#gamma} [GeV]", {}),
      Axis(65,94,159, "llphoton_refit_m", "m_{ll#gamma} [GeV]", {}),
      Axis(65,97,162, "llphoton_refit_m", "m_{ll#gamma} [GeV]", {}),
      Axis(65,99,164, "llphoton_refit_m", "m_{ll#gamma} [GeV]", {}),
      Axis(65,105,170, "llphoton_refit_m", "m_{ll#gamma} [GeV]", {}),
      Axis(13,100,165, "llphoton_refit_m", "m_{ll#gamma} [GeV]", {}),
      Axis(13,100,165, "llphoton_refit_m", "m_{ll#gamma} [GeV]", {}),
      Axis(13,100,165, "llphoton_refit_m", "m_{ll#gamma} [GeV]", {}),
      Axis(13,100,165, "llphoton_refit_m", "m_{ll#gamma} [GeV]", {}),
      Axis(65,105,170, "llphoton_refit_m", "m_{ll#gamma} [GeV]", {})
  };

  for (Axis& axis : baseline_plots)
    pm.Push<Hist1D>(axis, baseline, procsx1000, ops)
        .Weight(weight*w_sigx1000).Tag("zgsupp");
  for (Axis& axis : ggf_plots)
    pm.Push<Hist1D>(axis, baseline&&higgs_mass_region&&ggf, procsx300, ops)
        .Weight(weight*w_sigx300).Tag("zgsupp");
  for (Axis& axis : vbf_plots)
    pm.Push<Hist1D>(axis, baseline&&higgs_mass_region&&vbf, procs_vbf, ops)
        .Weight(weight*w_sigx100).Tag("zgsupp");
  for (unsigned icat = 0; icat < categories.size(); icat++) {
    pm.Push<Hist1D>(mllg_plots[icat], baseline&&categories[icat],
        procs_noscale, ops).Weight(weight).Tag("zgsupp");
  }

  pm.SetEnergyTag("138 fb^{-1} (13 TeV) + 62 fb^{-1} (13.6 TeV)");
  pm.MakePlots(1.0);

  return 0;
}
