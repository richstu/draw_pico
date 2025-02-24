/**
 * We're trying something new here: I will use this script to do quick studies, 
 * then transfer to dedicated files as needed based on maturity
 */

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <bitset>
#include <unistd.h>
#include <getopt.h>

#include "TError.h"
#include "TChain.h"
#include "TColor.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TMVA/Reader.h"
#include "TMVA/Configurable.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/mva_wrapper.hpp"
#include "core/named_func.hpp"
#include "core/named_func_utilities.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/utilities.hpp"
#include "zgamma/nano_functions.hpp"
#include "zgamma/zg_functions.hpp"
#include "zgamma/zg_utilities.hpp"

using std::cout;
using std::isnan;
using std::make_shared;
using std::set;
using std::shared_ptr;
using std::string;
using std::to_string;
using std::vector;
using NamedFuncUtilities::reduce_max;
using NamedFuncUtilities::reduce_maxfirst;
using NamedFuncUtilities::reduce_sum;
using NamedFuncUtilities::FilterNamedFunc;
using NamedFuncUtilities::MapNamedFunc;
using NamedFuncUtilities::MultiReduceNamedFunc;
using NamedFuncUtilities::ReduceNamedFunc;
using ZgFunctions::max_lep_miniso;
using ZgFunctions::trig_plateau_cuts;
using ZgFunctions::photon_relpterr;
using ZgFunctions::lead_lepton_eta;
using ZgFunctions::sublead_lepton_eta;
using ZgFunctions::w_years;
using ZgUtilities::get_phi;
using ZgUtilities::get_costheta;
using ZgUtilities::get_cosTheta;
using ZgUtilities::ZgSampleLoader;
using ZgUtilities::KinematicBdt;
using ZgUtilities::category_ggh4;
using ZgUtilities::category_ggh3;
using ZgUtilities::category_ggh2;
using ZgUtilities::category_ggh1;
using PlotOptTypes::OverflowType;
//const Process::Type data = Process::Type::data;
const Process::Type signal = Process::Type::signal;
const Process::Type background = Process::Type::background;

const float PI = acos(-1.0);

/**
 * Use linear fit to (x1,y1) and (x2,y2) to esetimate y from x
 */
float linear_interpolation(float x1, float y1, float x2, float y2, float x) {
  float dy = y2-y1;
  float dx = x2-x1;
  return dy/dx*(x-x1)+y1;
}

double double_abs(double x) {
  if (x<0) return -1.0*x;
  return x;
}

int main() {

  //--------------------------------------------------------------------------
  //                                    initialization
  //--------------------------------------------------------------------------

  //setup
  gErrorIgnoreLevel = 6000;
  //string lumi_tag = "138";
  string lumi_tag = "41";

  //vector<shared_ptr<Process>> procs_llskim_reduced = ZgSampleLoader()
  //    //.SetMacro("YEARS",{"2018"})
  //    .LoadSamples("txt/samples_zgamma.txt","llcr");

  string prod_folder("/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/");
  Palette mc_colors("txt/colors_zgamma.txt","default");
  std::set<int> years = {2017};
  vector<shared_ptr<Process>> procs_dy;
  procs_dy.push_back(Process::MakeShared<Baby_pico>("Z+Fake Photon", background, mc_colors("dyjets"),
                     attach_folder(prod_folder,years,"mc/merged_zgmc_ll",
                     {"*DYJetsToLL_M-50_TuneCP5_13TeV-amcatnlo*"}),
                     "use_event&&(trig_single_el||trig_single_mu||trig_double_el||trig_double_mu)"));

  years = {2016,2017,2018};
  vector<shared_ptr<Process>> procs_samples;
  procs_samples.push_back(Process::MakeShared<Baby_pico>("DYJets_amcatnlo", background, mc_colors("dyjets"),
                          attach_folder(prod_folder,years,"mc/merged_zgmc_ll",
                          {"*DYJetsToLL_M-50_TuneCP5_13TeV-amcatnlo*"}),
                          "use_event&&(trig_single_el||trig_single_mu||trig_double_el||trig_double_mu)"));
  procs_samples.push_back(Process::MakeShared<Baby_pico>("DYJets_lo", background, mc_colors("dyjets"),
                          attach_folder(prod_folder,years,"mc/merged_zgmc_ll",
                          {"*DYJetsToLL_M-50_TuneCP5_13TeV-madgraph*"}),
                          "use_event&&(trig_single_el||trig_single_mu||trig_double_el||trig_double_mu)"));
  procs_samples.push_back(Process::MakeShared<Baby_pico>("ZGToLLG_lowMLL", background, mc_colors("zgtollg"),
                          attach_folder(prod_folder,years,"mc/merged_zgmc_llg",
                          {"*ZGToLLG_01J_5f_lowMLL*"}),
                          "use_event&&(trig_single_el||trig_single_mu||trig_double_el||trig_double_mu)"));
  procs_samples.push_back(Process::MakeShared<Baby_pico>("ZGToLLG_nominal", background, mc_colors("zgtollg"),
                          attach_folder(prod_folder,years,"mc/merged_zgmc_llg",
                          {"*ZGToLLG_01J_5f_TuneCP5*"}),
                          "use_event&&(trig_single_el||trig_single_mu||trig_double_el||trig_double_mu)"));

  //std::vector<PlotOpt> ops_data = {PlotOpt("txt/plot_styles.txt","LinLumiDataRoot")}; 
  std::vector<PlotOpt> ops_linlumi = {PlotOpt("txt/plot_styles.txt","LinLumi")
                                      .Overflow(OverflowType::none)};
  std::vector<PlotOpt> ops_shapes = {PlotOpt("txt/plot_styles.txt","LinLumi")
                                     .Overflow(OverflowType::none),
                                     PlotOpt("txt/plot_styles.txt","Shapes")
                                     .Overflow(OverflowType::none)
                                     .FileExtensions({"pdf","root"})}; 

  std::vector<PlotOpt> ops_2d = {PlotOpt("txt/plot_styles.txt","LinLumi2DRoot")};

  //------------------------------------------------------------------------------------
  //                                   NamedFuncs
  //------------------------------------------------------------------------------------

  //weight with some clipping/regularization
  const NamedFunc weight_reg("weight_reg",[](const Baby &b) -> NamedFunc::ScalarType{
    if (fabs(b.weight())>5) return 0.0;
    return b.weight();
  });

  const NamedFunc photon_loosesig = NamedFunc("((photon_isScEtaEB&&photon_idmva>-0.9)||(photon_isScEtaEE&&photon_idmva>-0.9))&&photon_elveto&&(photon_pt>15)&&photon_drmin>0.3").Name("photon_drsig"); 
  const NamedFunc nloosephoton = ReduceNamedFunc(photon_loosesig, reduce_sum).Name("nloosephoton");
  const NamedFunc loosephoton_pt = FilterNamedFunc("photon_pt",photon_loosesig).Name("loosephoton_pt");
  const NamedFunc loosephoton_eta = FilterNamedFunc("photon_eta",photon_loosesig).Name("loosephoton_eta");
  const NamedFunc loosephoton_phi = FilterNamedFunc("photon_phi",photon_loosesig).Name("loosephoton_phi");
  const NamedFunc loosephoton_drmin = FilterNamedFunc("photon_drmin",photon_loosesig).Name("loosephoton_drmin");
  const NamedFunc loosephoton_drmax = FilterNamedFunc("photon_drmax",photon_loosesig).Name("loosephoton_drmax");
  const NamedFunc loosephoton_id80 = FilterNamedFunc("photon_id80",photon_loosesig).Name("loosephoton_id80");
  const NamedFunc loosephoton_idmva = FilterNamedFunc("photon_idmva",photon_loosesig).Name("loosephoton_idmva");
  const NamedFunc loosephoton_isScEtaEB = FilterNamedFunc("photon_isScEtaEB",photon_loosesig).Name("loosephoton_isScEtaEB");
  const NamedFunc lead_loosephoton_pt = ReduceNamedFunc(loosephoton_pt,reduce_max).Name("lead_loosephoton_pt");
  const NamedFunc lead_loosephoton_eta = MultiReduceNamedFunc(
      {loosephoton_pt,loosephoton_eta},reduce_maxfirst).Name("lead_loosephoton_eta");
  const NamedFunc lead_loosephoton_phi = MultiReduceNamedFunc(
      {loosephoton_pt,loosephoton_phi},reduce_maxfirst).Name("lead_loosephoton_phi");
  const NamedFunc lead_loosephoton_drmin = MultiReduceNamedFunc(
      {loosephoton_pt,loosephoton_drmin},reduce_maxfirst).Name("lead_loosephoton_drmin");
  const NamedFunc lead_loosephoton_drmax = MultiReduceNamedFunc(
      {loosephoton_pt,loosephoton_drmax},reduce_maxfirst).Name("lead_loosephoton_drmax");
  const NamedFunc lead_loosephoton_id80 = MultiReduceNamedFunc(
      {loosephoton_pt,loosephoton_id80},reduce_maxfirst).Name("lead_loosephoton_id80");
  const NamedFunc lead_loosephoton_idmva = MultiReduceNamedFunc(
      {loosephoton_pt,loosephoton_idmva},reduce_maxfirst).Name("lead_loosephoton_idmva");
  const NamedFunc lead_loosephoton_isScEtaEB = MultiReduceNamedFunc(
      {loosephoton_pt,loosephoton_isScEtaEB},reduce_maxfirst).Name("lead_loosephoton_isScEtaEB");
  const NamedFunc lead_loosephoton_abseta = MapNamedFunc(lead_loosephoton_eta,double_abs).Name("lead_loosephoton_abseta");

  //relative pt uncertainty of lead photon with ID requirement loosened
  const NamedFunc lead_loosephoton_relpterr("lead_loosephoton_relpterr",[photon_loosesig](const Baby &b) -> NamedFunc::ScalarType{
    vector<double> sig = photon_loosesig.GetVector(b);
    float max_pt = 0;
    unsigned max_idx = 0;
    for (unsigned iph = 0; iph < b.photon_pt()->size(); iph++) {
      if (sig[iph] && b.photon_pt()->at(iph) > max_pt) {
        max_pt = b.photon_pt()->at(iph);
        max_idx = iph;
      }
    }
    return b.photon_energyErr()->at(max_idx)/
        (b.photon_pt()->at(max_idx)*TMath::CosH(b.photon_eta()->at(max_idx)));
  });

  //llphoton mass with photon ID reqiurement loosened
  const NamedFunc llloosephoton_m("llloosephoton_m",[lead_loosephoton_pt, lead_loosephoton_eta, lead_loosephoton_phi](const Baby &b) -> NamedFunc::ScalarType{
    TLorentzVector z, ph;
    z.SetPtEtaPhiM(b.ll_pt()->at(0),b.ll_eta()->at(0),b.ll_phi()->at(0),b.ll_m()->at(0));
    ph.SetPtEtaPhiM(lead_loosephoton_pt.GetScalar(b), lead_loosephoton_eta.GetScalar(b),
                    lead_loosephoton_phi.GetScalar(b), 0.0);
    return (z+ph).M();
  });

  //llphoton pt with photon ID reqiurement loosened
  const NamedFunc llloosephoton_pt("llloosephoton_pt",[lead_loosephoton_pt, lead_loosephoton_eta, lead_loosephoton_phi](const Baby &b) -> NamedFunc::ScalarType{
    TLorentzVector z, ph;
    z.SetPtEtaPhiM(b.ll_pt()->at(0),b.ll_eta()->at(0),b.ll_phi()->at(0),b.ll_m()->at(0));
    ph.SetPtEtaPhiM(lead_loosephoton_pt.GetScalar(b), lead_loosephoton_eta.GetScalar(b),
                    lead_loosephoton_phi.GetScalar(b), 0.0);
    return (z+ph).Pt();
  });

  //llphoton costheta with photon dr reqiurement loosened
  const NamedFunc llloosephoton_costheta("llloosephoton_costheta",[lead_loosephoton_pt, lead_loosephoton_eta, lead_loosephoton_phi](const Baby &b) -> NamedFunc::ScalarType{
    TLorentzVector lm, lp, ph;
    if (b.ll_lepid()->at(0)==11) {
      if (b.el_charge()->at(b.ll_i1()->at(0))==1) {
        lm.SetPtEtaPhiM(b.el_pt()->at(b.ll_i2()->at(0)),b.el_eta()->at(b.ll_i2()->at(0)),
                        b.el_phi()->at(b.ll_i2()->at(0)),0.000511);
        lp.SetPtEtaPhiM(b.el_pt()->at(b.ll_i1()->at(0)),b.el_eta()->at(b.ll_i1()->at(0)),
                        b.el_phi()->at(b.ll_i1()->at(0)),0.000511);
      }
      else {
        lp.SetPtEtaPhiM(b.el_pt()->at(b.ll_i2()->at(0)),b.el_eta()->at(b.ll_i2()->at(0)),
                        b.el_phi()->at(b.ll_i2()->at(0)),0.000511);
        lm.SetPtEtaPhiM(b.el_pt()->at(b.ll_i1()->at(0)),b.el_eta()->at(b.ll_i1()->at(0)),
                        b.el_phi()->at(b.ll_i1()->at(0)),0.000511);
      }
    }
    else {
      if (b.mu_charge()->at(b.ll_i1()->at(0))==1) {
        lm.SetPtEtaPhiM(b.mu_pt()->at(b.ll_i2()->at(0)),b.mu_eta()->at(b.ll_i2()->at(0)),
                        b.mu_phi()->at(b.ll_i2()->at(0)),0.106);
        lp.SetPtEtaPhiM(b.mu_pt()->at(b.ll_i1()->at(0)),b.mu_eta()->at(b.ll_i1()->at(0)),
                        b.mu_phi()->at(b.ll_i1()->at(0)),0.106);
      }
      else {
        lp.SetPtEtaPhiM(b.mu_pt()->at(b.ll_i2()->at(0)),b.mu_eta()->at(b.ll_i2()->at(0)),
                        b.mu_phi()->at(b.ll_i2()->at(0)),0.106);
        lm.SetPtEtaPhiM(b.mu_pt()->at(b.ll_i1()->at(0)),b.mu_eta()->at(b.ll_i1()->at(0)),
                        b.mu_phi()->at(b.ll_i1()->at(0)),0.106);
      }
    }
    ph.SetPtEtaPhiM(lead_loosephoton_pt.GetScalar(b), lead_loosephoton_eta.GetScalar(b),
                    lead_loosephoton_phi.GetScalar(b), 0.0);
    return get_costheta(lm,lp,ph);
  });

  //llphoton cosTheta with photon dr reqiurement loosened
  const NamedFunc llloosephoton_cosTheta("llloosephoton_coscaptheta",[lead_loosephoton_pt, lead_loosephoton_eta, lead_loosephoton_phi](const Baby &b) -> NamedFunc::ScalarType{
    TLorentzVector lm, lp, ph;
    if (b.ll_lepid()->at(0)==11) {
      if (b.el_charge()->at(b.ll_i1()->at(0))==1) {
        lm.SetPtEtaPhiM(b.el_pt()->at(b.ll_i2()->at(0)),b.el_eta()->at(b.ll_i2()->at(0)),
                        b.el_phi()->at(b.ll_i2()->at(0)),0.000511);
        lp.SetPtEtaPhiM(b.el_pt()->at(b.ll_i1()->at(0)),b.el_eta()->at(b.ll_i1()->at(0)),
                        b.el_phi()->at(b.ll_i1()->at(0)),0.000511);
      }
      else {
        lp.SetPtEtaPhiM(b.el_pt()->at(b.ll_i2()->at(0)),b.el_eta()->at(b.ll_i2()->at(0)),
                        b.el_phi()->at(b.ll_i2()->at(0)),0.000511);
        lm.SetPtEtaPhiM(b.el_pt()->at(b.ll_i1()->at(0)),b.el_eta()->at(b.ll_i1()->at(0)),
                        b.el_phi()->at(b.ll_i1()->at(0)),0.000511);
      }
    }
    else {
      if (b.mu_charge()->at(b.ll_i1()->at(0))==1) {
        lm.SetPtEtaPhiM(b.mu_pt()->at(b.ll_i2()->at(0)),b.mu_eta()->at(b.ll_i2()->at(0)),
                        b.mu_phi()->at(b.ll_i2()->at(0)),0.106);
        lp.SetPtEtaPhiM(b.mu_pt()->at(b.ll_i1()->at(0)),b.mu_eta()->at(b.ll_i1()->at(0)),
                        b.mu_phi()->at(b.ll_i1()->at(0)),0.106);
      }
      else {
        lp.SetPtEtaPhiM(b.mu_pt()->at(b.ll_i2()->at(0)),b.mu_eta()->at(b.ll_i2()->at(0)),
                        b.mu_phi()->at(b.ll_i2()->at(0)),0.106);
        lm.SetPtEtaPhiM(b.mu_pt()->at(b.ll_i1()->at(0)),b.mu_eta()->at(b.ll_i1()->at(0)),
                        b.mu_phi()->at(b.ll_i1()->at(0)),0.106);
      }
    }
    ph.SetPtEtaPhiM(lead_loosephoton_pt.GetScalar(b), lead_loosephoton_eta.GetScalar(b),
                    lead_loosephoton_phi.GetScalar(b), 0.0);
    return get_cosTheta(lm,lp,ph);
  });

  //llphoton phi with photon dr reqiurement loosened
  const NamedFunc llloosephoton_psi("llloosephoton_psi",[lead_loosephoton_pt, lead_loosephoton_eta, lead_loosephoton_phi](const Baby &b) -> NamedFunc::ScalarType{
    TLorentzVector lm, lp, ph;
    if (b.ll_lepid()->at(0)==11) {
      if (b.el_charge()->at(b.ll_i1()->at(0))==1) {
        lm.SetPtEtaPhiM(b.el_pt()->at(b.ll_i2()->at(0)),b.el_eta()->at(b.ll_i2()->at(0)),
                        b.el_phi()->at(b.ll_i2()->at(0)),0.000511);
        lp.SetPtEtaPhiM(b.el_pt()->at(b.ll_i1()->at(0)),b.el_eta()->at(b.ll_i1()->at(0)),
                        b.el_phi()->at(b.ll_i1()->at(0)),0.000511);
      }
      else {
        lp.SetPtEtaPhiM(b.el_pt()->at(b.ll_i2()->at(0)),b.el_eta()->at(b.ll_i2()->at(0)),
                        b.el_phi()->at(b.ll_i2()->at(0)),0.000511);
        lm.SetPtEtaPhiM(b.el_pt()->at(b.ll_i1()->at(0)),b.el_eta()->at(b.ll_i1()->at(0)),
                        b.el_phi()->at(b.ll_i1()->at(0)),0.000511);
      }
    }
    else {
      if (b.mu_charge()->at(b.ll_i1()->at(0))==1) {
        lm.SetPtEtaPhiM(b.mu_pt()->at(b.ll_i2()->at(0)),b.mu_eta()->at(b.ll_i2()->at(0)),
                        b.mu_phi()->at(b.ll_i2()->at(0)),0.106);
        lp.SetPtEtaPhiM(b.mu_pt()->at(b.ll_i1()->at(0)),b.mu_eta()->at(b.ll_i1()->at(0)),
                        b.mu_phi()->at(b.ll_i1()->at(0)),0.106);
      }
      else {
        lp.SetPtEtaPhiM(b.mu_pt()->at(b.ll_i2()->at(0)),b.mu_eta()->at(b.ll_i2()->at(0)),
                        b.mu_phi()->at(b.ll_i2()->at(0)),0.106);
        lm.SetPtEtaPhiM(b.mu_pt()->at(b.ll_i1()->at(0)),b.mu_eta()->at(b.ll_i1()->at(0)),
                        b.mu_phi()->at(b.ll_i1()->at(0)),0.106);
      }
    }
    ph.SetPtEtaPhiM(lead_loosephoton_pt.GetScalar(b), lead_loosephoton_eta.GetScalar(b),
                    lead_loosephoton_phi.GetScalar(b), 0.0);
    float phi = get_phi(lm,lp,ph);
    if (phi > PI) phi -= 2.0*PI;
    return phi;
  });

  //TRandom3 rng;
  //random value for IDMVA that is over threshold
  const NamedFunc lead_loosephoton_randomidmva("lead_loosephoton_randomidmva",[lead_loosephoton_eta](const Baby &b) -> NamedFunc::ScalarType{
    if (fabs(lead_loosephoton_eta.GetScalar(b))<1.5)
      return static_cast<float>(b.event()%40)*0.58/40.0+0.42;
    return static_cast<float>(b.event()%40)*0.86/40.0*0.14;
    //return rng.Uniform(0.4,1.0);
    //return rng.Uniform(0.14,1.0);
  });

  //const vector<vector<float>> pt_extrap_factors = {{0.143,0.111,0.089,0.087,0.090},
  //                                                 {0.133,0.117,0.108,0.108,0.094},
  //                                                 {0.158,0.122,0.109,0.118,0.117},
  //                                                 {0.278,0.220,0.241,0.253,0.243}};
  //const vector<float> pt_extrap_bin_centers = {16.25,18.75,22.5,27.5,40.0};
  //const vector<float> eta_extrap_bin_centers = {0.4,1.15,1.75,2.25};
  ////weight for extrapolation
  //const NamedFunc w_extrap_pt("w_extrap_pt",[lead_loosephoton_pt,lead_loosephoton_eta,pt_extrap_factors,pt_extrap_bin_centers,eta_extrap_bin_centers](const Baby &b) -> NamedFunc::ScalarType{
  //  float abseta = fabs(lead_loosephoton_eta.GetScalar(b));
  //  float pt = lead_loosephoton_pt.GetScalar(b);
  //  float w = 0;
  //  unsigned pt_bin_prev = 0;
  //  unsigned pt_bin_next = 0;
  //  unsigned eta_bin_prev = 0;
  //  unsigned eta_bin_next = 0;
  //  //find pt, eta bins
  //  if (pt<16.25)         {pt_bin_prev = 0; pt_bin_next = 0;}
  //  else if (pt<18.75)    {pt_bin_prev = 0; pt_bin_next = 1;}
  //  else if (pt<22.5)     {pt_bin_prev = 1; pt_bin_next = 2;}
  //  else if (pt<27.5)     {pt_bin_prev = 2; pt_bin_next = 3;}
  //  else if (pt<40.0)     {pt_bin_prev = 3; pt_bin_next = 4;}
  //  else                  {pt_bin_prev = 4; pt_bin_next = 4;}
  //  if (abseta<0.4)       {eta_bin_prev = 0; eta_bin_next = 1; }
  //  else if (abseta<1.15) {eta_bin_prev = 0; eta_bin_next = 1; }
  //  else if (abseta<1.75) {eta_bin_prev = 1; eta_bin_next = 2; }
  //  else if (abseta<2.25) {eta_bin_prev = 2; eta_bin_next = 3; } 
  //  else                  {eta_bin_prev = 2; eta_bin_next = 3; }
  //  //get values and interpolate
  //  float w_etalo = 0.0, w_etahi = 0.0;
  //  if (pt_bin_prev != pt_bin_next) {
  //    w_etalo = linear_interpolation(pt_extrap_bin_centers[pt_bin_prev],
  //                                   pt_extrap_factors[eta_bin_prev][pt_bin_prev],
  //                                   pt_extrap_bin_centers[pt_bin_next],
  //                                   pt_extrap_factors[eta_bin_prev][pt_bin_next],
  //                                   pt);
  //    w_etahi = linear_interpolation(pt_extrap_bin_centers[pt_bin_prev],
  //                                   pt_extrap_factors[eta_bin_next][pt_bin_prev],
  //                                   pt_extrap_bin_centers[pt_bin_next],
  //                                   pt_extrap_factors[eta_bin_next][pt_bin_next],
  //                                   pt);
  //  }
  //  else {
  //    w_etalo = pt_extrap_factors[eta_bin_prev][pt_bin_prev];
  //    w_etahi = pt_extrap_factors[eta_bin_next][pt_bin_prev];
  //  }
  //  if (eta_bin_prev != eta_bin_next) {
  //    w = linear_interpolation(eta_extrap_bin_centers[eta_bin_prev],
  //                             w_etalo,
  //                             eta_extrap_bin_centers[eta_bin_next],
  //                             w_etahi,
  //                             abseta);
  //  }
  //  else {
  //    w = w_etalo;
  //  }
  //  return w;
  //});

  //const vector<vector<vector<float>>> extrap_factors = 
  //    { { { 0.141, 0.147, 0.091, 0.099, 0.097 } ,
  //        { 0.173, 0.150, 0.098, 0.083, 0.090 } ,
  //        { 0.145, 0.133, 0.091, 0.094, 0.107 } ,
  //        { 0.301, 0.255, 0.227, 0.274, 0.274 } },
  //      { { 0.138, 0.130, 0.099, 0.092, 0.109 } ,
  //        { 0.123, 0.107, 0.121, 0.098, 0.102 } ,
  //        { 0.232, 0.080, 0.111, 0.112, 0.135 } ,
  //        { 0.236, 0.195, 0.230, 0.293, 0.154 } },
  //      { { 0.148, 0.102, 0.101, 0.076, 0.088 } ,
  //        { 0.154, 0.110, 0.107, 0.121, 0.131 } ,
  //        { 0.132, 0.168, 0.145, 0.122, 0.108 } ,
  //        { 0.298, 0.305, 0.289, 0.335, 0.336 } } ,
  //      { { 0.161, 0.117, 0.078, 0.093, 0.097 } ,
  //        { 0.141, 0.124, 0.121, 0.147, 0.094 } ,
  //        { 0.142, 0.085, 0.132, 0.118, 0.110 } ,
  //        { 0.187, 0.219, 0.172, 0.174, 0.166 } } ,
  //      { { 0.107, 0.100, 0.092, 0.074, 0.074 } ,
  //        { 0.110, 0.078, 0.080, 0.126, 0.087 } ,
  //        { 0.164, 0.112, 0.170, 0.093, 0.087 } ,
  //        { 0.181, 0.258, 0.287, 0.181, 0.551 } } ,
  //      { { 0.120, 0.100, 0.076, 0.062, 0.061 } ,
  //        { 0.124, 0.081, 0.114, 0.098, 0.125 } ,
  //        { 0.146, 0.112, 0.105, 0.157, 0.125 } ,
  //        { 0.368, 0.189, 0.136, 0.230, 0.146 } } };
  //const vector<float> pt_extrap_bin_centers = {16.25,18.75,22.5,27.5,40.0};
  //const vector<float> eta_extrap_bin_centers = {0.4,1.15,1.75,2.25};
  //const vector<float> cosTheta_extrap_bin_centers = {0.1,0.3,0.5,0.7,0.85,0.95};
  ////weight to correct pt, eta, and cosTheta after extrapolation
  //const NamedFunc w_extrap_ptetacosT("w_extrap",
  //                                   [lead_loosephoton_pt,lead_loosephoton_eta,
  //                                    llloosephoton_cosTheta,
  //                                    extrap_factors,
  //                                    pt_extrap_bin_centers,
  //                                    eta_extrap_bin_centers,
  //                                    cosTheta_extrap_bin_centers]
  //                                    (const Baby &b) -> NamedFunc::ScalarType{
  //  float abseta = fabs(lead_loosephoton_eta.GetScalar(b));
  //  float pt = lead_loosephoton_pt.GetScalar(b);
  //  float abs_cosTheta = fabs(llloosephoton_cosTheta.GetScalar(b));
  //  float pt_prv = 0, pt_nxt = 0;
  //  float eta_prv = 0, eta_nxt = 0;
  //  float cT_prv = 0, cT_nxt = 0;
  //  unsigned ipt_prv = 0, ipt_nxt = 0;
  //  unsigned ieta_prv = 0, ieta_nxt = 0;
  //  unsigned icT_prv = 0, icT_nxt = 0;
  //  //find bins
  //  if (pt <= 16.25) {pt_prv = 15.0; pt_nxt = 16.25; ipt_prv = 0; ipt_nxt = 0;}
  //  else if (pt >= 40.0) {pt_prv = 40.0; pt_nxt = 50.0; ipt_prv = 4; ipt_nxt = 4;}
  //  else {
  //    for (unsigned ipt = 1; ipt<pt_extrap_bin_centers.size(); ipt++) {
  //      if (pt < pt_extrap_bin_centers[ipt]) {
  //        pt_prv = pt_extrap_bin_centers[ipt-1];
  //        pt_nxt = pt_extrap_bin_centers[ipt];
  //        ipt_prv = ipt-1;
  //        ipt_nxt = ipt;
  //        break;
  //      }
  //    }
  //  }
  //  if (abseta <= 0.4) {eta_prv = 0.4; eta_nxt = 1.15; ieta_prv = 0; ieta_nxt = 1;}
  //  else if (abseta >= 2.25) {eta_prv = 1.75; eta_nxt = 2.25; ieta_prv = 2; ieta_nxt = 3;}
  //  else {
  //    for (unsigned ieta = 1; ieta<eta_extrap_bin_centers.size(); ieta++) {
  //      if (abseta < eta_extrap_bin_centers[ieta]) {
  //        eta_prv = eta_extrap_bin_centers[ieta-1];
  //        eta_nxt = eta_extrap_bin_centers[ieta];
  //        ieta_prv = ieta-1;
  //        ieta_nxt = ieta;
  //        break;
  //      }
  //    }
  //  }
  //  if (abs_cosTheta <= 0.0) {cT_prv = 0.1; cT_nxt = 0.3; icT_prv = 0; icT_nxt = 1;}
  //  else if (abs_cosTheta >= 0.95) {cT_prv = 0.95; cT_nxt = 0.99; icT_prv = 5; icT_nxt = 5;}
  //  else {
  //    for (unsigned icT = 1; icT<cosTheta_extrap_bin_centers.size(); icT++) {
  //      if (abs_cosTheta < cosTheta_extrap_bin_centers[icT]) {
  //        cT_prv = cosTheta_extrap_bin_centers[icT-1];
  //        cT_nxt = cosTheta_extrap_bin_centers[icT];
  //        icT_prv = icT-1;
  //        icT_nxt = icT;
  //        break;
  //      }
  //    }
  //  }
  //  //do pT interpolation
  //  float w_etalo_cTlo = linear_interpolation(pt_prv, extrap_factors[icT_prv][ieta_prv][ipt_prv],
  //                                            pt_nxt, extrap_factors[icT_prv][ieta_prv][ipt_nxt],
  //                                            pt);
  //  float w_etahi_cTlo = linear_interpolation(pt_prv, extrap_factors[icT_prv][ieta_nxt][ipt_prv],
  //                                            pt_nxt, extrap_factors[icT_prv][ieta_nxt][ipt_nxt],
  //                                            pt);
  //  float w_etalo_cThi = linear_interpolation(pt_prv, extrap_factors[icT_nxt][ieta_prv][ipt_prv],
  //                                            pt_nxt, extrap_factors[icT_nxt][ieta_prv][ipt_nxt],
  //                                            pt);
  //  float w_etahi_cThi = linear_interpolation(pt_prv, extrap_factors[icT_nxt][ieta_nxt][ipt_prv],
  //                                            pt_nxt, extrap_factors[icT_nxt][ieta_nxt][ipt_nxt],
  //                                            pt);
  //  //do eta interpolation
  //  float w_cTlo = linear_interpolation(eta_prv, w_etalo_cTlo,
  //                                      eta_nxt, w_etahi_cTlo,
  //                                      abseta);
  //  float w_cThi = linear_interpolation(eta_prv, w_etalo_cThi,
  //                                      eta_nxt, w_etahi_cThi,
  //                                      abseta);
  //  //do cosTheta interpolation
  //  float w = linear_interpolation(cT_prv, w_cTlo, cT_nxt, w_cThi, abs_cosTheta);
  //  return w;
  //});

  vector<vector<vector<float>>> extrap_factors = { { { 1.218, 0.886, 0.457, 0.345, 0.257 } ,
                                                     { 1.457, 0.895, 0.606, 0.608, 0.400 } ,
                                                     { 0.000, 0.000, 0.000, 0.000, 0.000 } ,
                                                     { 0.000, 0.000, 0.000, 0.000, 0.000 } } ,
                                                   { { 0.210, 0.159, 0.107, 0.104, 0.084 } ,
                                                     { 0.358, 0.256, 0.228, 0.181, 0.151 } ,
                                                     { 1.022, 0.716, 0.672, 0.676, 0.575 } ,
                                                     { 2.421, 1.362, 1.481, 1.070, 0.818 } } ,
                                                   { { 0.071, 0.045, 0.029, 0.012, 0.020 } ,
                                                     { 0.101, 0.089, 0.074, 0.043, 0.029 } ,
                                                     { 0.219, 0.149, 0.112, 0.091, 0.056 } ,
                                                     { 0.351, 0.253, 0.232, 0.162, 0.100 } } ,
                                                   { { 0.014, 0.007, 0.006, 0.005, 0.009 } ,
                                                     { 0.019, 0.016, 0.010, 0.015, 0.009 } ,
                                                     { 0.027, 0.013, 0.011, 0.018, 0.010 } ,
                                                     { 0.034, 0.035, 0.025, 0.019, 0.029 } } ,
                                                   { { 0.002, 0.000, 0.001, 0.003, 0.000 } ,
                                                     { 0.000, 0.000, 0.009, 0.000, 0.000 } ,
                                                     { 0.000, 0.000, 0.000, 0.000, 0.000 } ,
                                                     { 0.008, 0.000, 0.000, 0.000, 0.000 } } };

  const vector<float> pt_extrap_bin_centers = {16.25,18.75,22.5,27.5,40.0};
  const vector<float> eta_extrap_bin_centers = {0.4,1.15,1.75,2.25};
  const vector<float> res_extrap_bin_centers = {0.01,0.03,0.06,0.0115,1.075};
  const NamedFunc w_extrap_ptetares("w_extrap_ptetares",
                                     [lead_loosephoton_pt,lead_loosephoton_eta,
                                      lead_loosephoton_relpterr,
                                      extrap_factors,
                                      pt_extrap_bin_centers,
                                      eta_extrap_bin_centers,
                                      res_extrap_bin_centers]
                                      (const Baby &b) -> NamedFunc::ScalarType{
    float abseta = fabs(lead_loosephoton_eta.GetScalar(b));
    float pt = lead_loosephoton_pt.GetScalar(b);
    float relpterr = lead_loosephoton_relpterr.GetScalar(b);
    float pt_prv = 0, pt_nxt = 0;
    float eta_prv = 0, eta_nxt = 0;
    float res_prv = 0, res_nxt = 0;
    unsigned ipt_prv = 0, ipt_nxt = 0;
    unsigned ieta_prv = 0, ieta_nxt = 0;
    unsigned ires_prv = 0, ires_nxt = 0;
    //find bins
    if (pt <= 16.25) {pt_prv = 15.0; pt_nxt = 16.25; ipt_prv = 0; ipt_nxt = 0;}
    else if (pt >= 40.0) {pt_prv = 40.0; pt_nxt = 50.0; ipt_prv = 4; ipt_nxt = 4;}
    else {
      for (unsigned ipt = 1; ipt<pt_extrap_bin_centers.size(); ipt++) {
        if (pt < pt_extrap_bin_centers[ipt]) {
          pt_prv = pt_extrap_bin_centers[ipt-1];
          pt_nxt = pt_extrap_bin_centers[ipt];
          ipt_prv = ipt-1;
          ipt_nxt = ipt;
          break;
        }
      }
    }
    if (abseta <= 0.4) {eta_prv = 0.4; eta_nxt = 1.15; ieta_prv = 0; ieta_nxt = 1;}
    else if (abseta >= 2.25) {eta_prv = 1.75; eta_nxt = 2.25; ieta_prv = 2; ieta_nxt = 3;}
    else {
      for (unsigned ieta = 1; ieta<eta_extrap_bin_centers.size(); ieta++) {
        if (abseta < eta_extrap_bin_centers[ieta]) {
          eta_prv = eta_extrap_bin_centers[ieta-1];
          eta_nxt = eta_extrap_bin_centers[ieta];
          ieta_prv = ieta-1;
          ieta_nxt = ieta;
          break;
        }
      }
    }
    if (relpterr <= 0.01) {res_prv = 0.01; res_nxt = 0.03; ires_prv = 0; ires_nxt = 1;}
    else if (relpterr >= 1.075) {res_prv = 1.075; res_nxt = 2.0; ires_prv = 5; ires_nxt = 5;}
    else {
      for (unsigned ires = 1; ires<res_extrap_bin_centers.size(); ires++) {
        if (relpterr < res_extrap_bin_centers[ires]) {
          res_prv = res_extrap_bin_centers[ires-1];
          res_nxt = res_extrap_bin_centers[ires];
          ires_prv = ires-1;
          ires_nxt = ires;
          break;
        }
      }
    }
    //do pT interpolation
    float w_etalo_reslo = linear_interpolation(pt_prv, extrap_factors[ires_prv][ieta_prv][ipt_prv],
                                              pt_nxt, extrap_factors[ires_prv][ieta_prv][ipt_nxt],
                                              pt);
    float w_etahi_reslo = linear_interpolation(pt_prv, extrap_factors[ires_prv][ieta_nxt][ipt_prv],
                                              pt_nxt, extrap_factors[ires_prv][ieta_nxt][ipt_nxt],
                                              pt);
    float w_etalo_reshi = linear_interpolation(pt_prv, extrap_factors[ires_nxt][ieta_prv][ipt_prv],
                                              pt_nxt, extrap_factors[ires_nxt][ieta_prv][ipt_nxt],
                                              pt);
    float w_etahi_reshi = linear_interpolation(pt_prv, extrap_factors[ires_nxt][ieta_nxt][ipt_prv],
                                              pt_nxt, extrap_factors[ires_nxt][ieta_nxt][ipt_nxt],
                                              pt);
    //do eta interpolation
    float w_reslo = linear_interpolation(eta_prv, w_etalo_reslo,
                                        eta_nxt, w_etahi_reslo,
                                        abseta);
    float w_reshi = linear_interpolation(eta_prv, w_etalo_reshi,
                                        eta_nxt, w_etahi_reshi,
                                        abseta);
    //do resolution interpolation
    float w = linear_interpolation(res_prv, w_reslo, res_nxt, w_reshi, relpterr);
    return w;
  });

  const NamedFunc january_baseline_looseidmva = NamedFunc("(nel>=2||nmu>=2)" && nloosephoton>=1 && "(ll_m[0]>81 && ll_m[0]<101) && (ll_charge[0]==0)" && ((lead_loosephoton_pt/llloosephoton_m)>15.0/110.0) && ((llloosephoton_m+"ll_m[0]")>185) && trig_plateau_cuts).Name("baseline");

  const NamedFunc january_baseline = NamedFunc("(nel>=2||nmu>=2) && nphoton>=1 && (ll_m[0]>81 && ll_m[0]<101) && (ll_charge[0]==0) && ((photon_pt[0]/llphoton_m[0])>15.0/110.0) && ((llphoton_m[0]+ll_m[0])>185)" && trig_plateau_cuts).Name("baseline");

  const NamedFunc blind_sr_data = NamedFunc("type>=1000||(llphoton_m[0]<120||llphoton_m[0]>130)").Name("blind_sr_data"); //n.b. doesn't work for run 3

  const NamedFunc blind_sr = NamedFunc("(llphoton_m[0]<120||llphoton_m[0]>130)").Name("blind_sr"); //n.b. doesn't work for run 3

  const NamedFunc rel_ptllg = "llphoton_pt[0]/llphoton_m[0]";
  const NamedFunc cos_captheta = NamedFunc("llphoton_cosTheta[0]").Name("llphoton_coscaptheta");

  shared_ptr<MVAWrapper> kinematic_bdt = KinematicBdt();
  NamedFunc kinematic_bdt_score = kinematic_bdt->GetDiscriminant();
  NamedFunc gghall = NamedFunc(kinematic_bdt_score>-0.1);
  NamedFunc ggh4 = category_ggh4(kinematic_bdt);
  NamedFunc ggh3 = category_ggh3(kinematic_bdt);
  NamedFunc ggh2 = category_ggh2(kinematic_bdt);
  NamedFunc ggh1 = category_ggh1(kinematic_bdt);

  std::shared_ptr<MVAWrapper> loose_bdt = std::make_shared<MVAWrapper>("loose_bdt");
  loose_bdt->SetVariable("photon_mva",lead_loosephoton_idmva);
  loose_bdt->SetVariable("min_dR",lead_loosephoton_drmin);
  loose_bdt->SetVariable("max_dR",lead_loosephoton_drmax);
  loose_bdt->SetVariable("pt_mass",llloosephoton_pt/llloosephoton_m);
  loose_bdt->SetVariable("cosTheta",llloosephoton_cosTheta);
  loose_bdt->SetVariable("costheta",llloosephoton_costheta);
  loose_bdt->SetVariable("phi",llloosephoton_psi);
  loose_bdt->SetVariable("photon_res",lead_loosephoton_relpterr);
  loose_bdt->SetVariable("photon_rapidity",lead_loosephoton_eta);
  loose_bdt->SetVariable("l1_rapidity",lead_lepton_eta);
  loose_bdt->SetVariable("l2_rapidity",sublead_lepton_eta);
  loose_bdt->BookMVA("/homes/oshiro/analysis/small_phys_utils/dataset/weights/shuffled_phidcomp_post_phidcomp_post_BDT.weights.xml");
  NamedFunc loose_bdt_score = loose_bdt->GetDiscriminant();
  NamedFunc ggh4_loose = category_ggh4(loose_bdt);
  NamedFunc ggh3_loose = category_ggh3(loose_bdt);
  NamedFunc ggh2_loose = category_ggh2(loose_bdt);
  NamedFunc ggh1_loose = category_ggh1(loose_bdt);

  std::shared_ptr<MVAWrapper> fake_bdt = std::make_shared<MVAWrapper>("fake_bdt");
  fake_bdt->SetVariable("photon_mva",lead_loosephoton_randomidmva);
  fake_bdt->SetVariable("min_dR",lead_loosephoton_drmin);
  fake_bdt->SetVariable("max_dR",lead_loosephoton_drmax);
  fake_bdt->SetVariable("pt_mass",llloosephoton_pt/llloosephoton_m);
  fake_bdt->SetVariable("cosTheta",llloosephoton_cosTheta);
  fake_bdt->SetVariable("costheta",llloosephoton_costheta);
  fake_bdt->SetVariable("phi",llloosephoton_psi);
  fake_bdt->SetVariable("photon_res",lead_loosephoton_relpterr);
  fake_bdt->SetVariable("photon_rapidity",lead_loosephoton_eta);
  fake_bdt->SetVariable("l1_rapidity",lead_lepton_eta);
  fake_bdt->SetVariable("l2_rapidity",sublead_lepton_eta);
  fake_bdt->BookMVA("/homes/oshiro/analysis/small_phys_utils/dataset/weights/shuffled_phidcomp_post_phidcomp_post_BDT.weights.xml");
  NamedFunc fake_bdt_score = fake_bdt->GetDiscriminant();
  NamedFunc gghall_fake = NamedFunc(fake_bdt_score>-0.1);
  NamedFunc ggh4_fake = category_ggh4(fake_bdt);
  NamedFunc ggh3_fake = category_ggh3(fake_bdt);
  NamedFunc ggh2_fake = category_ggh2(fake_bdt);
  NamedFunc ggh1_fake = category_ggh1(fake_bdt);
  
  //------------------------------------------------------------------------------------
  //                                   plots and tables
  //------------------------------------------------------------------------------------

  PlotMaker pm;
  pm.multithreaded_ = false; //using BDT

  //check which MVA variables are correlated with mllg
  //pm.Push<Hist2D>(Axis(30,0.0,2.0, rel_ptllg, "p_{T}^{ll#gamma}/m_{ll#gamma}", {}), 
  //                Axis(30,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
  //                january_baseline, procs_dy, ops_linlumi)
  //                .Weight(weight_reg*w_years).Tag("zgfake");
  //pm.Push<Hist2D>(Axis(30,0.5,1.0, "photon_idmva[0]", "Photon IDMVA", {}), 
  //                Axis(30,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
  //                january_baseline, procs_dy, ops_linlumi)
  //                .Weight(weight_reg*w_years).Tag("zgfake");
  //pm.Push<Hist2D>(Axis(30,0.0,5.0, "photon_drmin[0]", "Photon #Delta R_{min}(l)", {}), 
  //                Axis(30,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
  //                january_baseline, procs_dy, ops_linlumi)
  //                .Weight(weight_reg*w_years).Tag("zgfake");
  //pm.Push<Hist2D>(Axis(30,0.0,6.0, "photon_drmax[0]", "Photon #Delta R_{max}(l)", {}), 
  //                Axis(30,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
  //                january_baseline, procs_dy, ops_linlumi)
  //                .Weight(weight_reg*w_years).Tag("zgfake");
  //pm.Push<Hist2D>(Axis(30,-1.0,1.0, cos_captheta, "cos(#Theta)", {}), 
  //                Axis(30,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
  //                january_baseline, procs_dy, ops_linlumi)
  //                .Weight(weight_reg*w_years).Tag("zgfake");
  //pm.Push<Hist2D>(Axis(30,-1.0,1.0, "llphoton_costheta[0]", "cos(#theta)", {}), 
  //                Axis(30,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
  //                january_baseline, procs_dy, ops_linlumi)
  //                .Weight(weight_reg*w_years).Tag("zgfake");
  //pm.Push<Hist2D>(Axis(30,-3.1416,3.1416, "llphoton_psi[0]", "#psi", {}), 
  //                Axis(30,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
  //                january_baseline, procs_dy, ops_linlumi)
  //                .Weight(weight_reg*w_years).Tag("zgfake");
  //pm.Push<Hist2D>(Axis(30,0.0,0.5, photon_relpterr, "Photon relative resolution", {}), 
  //                Axis(30,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
  //                january_baseline, procs_dy, ops_linlumi)
  //                .Weight(weight_reg*w_years).Tag("zgfake");
  //pm.Push<Hist2D>(Axis(30,-2.5,2.5, "photon_eta[0]", "Photon #eta", {}), 
  //                Axis(30,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
  //                january_baseline, procs_dy, ops_linlumi)
  //                .Weight(weight_reg*w_years).Tag("zgfake");
  //pm.Push<Hist2D>(Axis(30,-2.5,2.5, lead_lepton_eta, "Lead lepton #eta", {}), 
  //                Axis(30,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
  //                january_baseline, procs_dy, ops_linlumi)
  //                .Weight(weight_reg*w_years).Tag("zgfake");
  //pm.Push<Hist2D>(Axis(30,-2.5,2.5, sublead_lepton_eta, "Sublead lepton #eta", {}), 
  //                Axis(30,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
  //                january_baseline, procs_dy, ops_linlumi)
  //                .Weight(weight_reg*w_years).Tag("zgfake");

  ////randomize photon idmva and add pt-eta-res weighting
  //pm.Push<Hist1D>(Axis(60,-1.0,1.0, loosephoton_idmva, "Photon IDMVA", {}), 
  //                january_baseline_looseidmva, procs_dy, ops_linlumi)
  //                .Weight(weight_reg*w_extrap_ptetares*w_years).Tag("zgfake_est");
  //pm.Push<Hist1D>(Axis(80,100.0,180.0, llloosephoton_m, "m_{ll#gamma} [GeV]", {}), 
  //                january_baseline_looseidmva&&lead_loosephoton_id80&&ggh1_loose, procs_dy, ops_linlumi)
  //                .Weight(weight_reg*w_years).Tag("zgfake_est");
  //pm.Push<Hist1D>(Axis(80,100.0,180.0, llloosephoton_m, "m_{ll#gamma} [GeV]", {}), 
  //                january_baseline_looseidmva&&lead_loosephoton_id80&&ggh2_loose, procs_dy, ops_linlumi)
  //                .Weight(weight_reg*w_years).Tag("zgfake_est");
  //pm.Push<Hist1D>(Axis(80,100.0,180.0, llloosephoton_m, "m_{ll#gamma} [GeV]", {}), 
  //                january_baseline_looseidmva&&lead_loosephoton_id80&&ggh3_loose, procs_dy, ops_linlumi)
  //                .Weight(weight_reg*w_years).Tag("zgfake_est");
  //pm.Push<Hist1D>(Axis(80,100.0,180.0, llloosephoton_m, "m_{ll#gamma} [GeV]", {}), 
  //                january_baseline_looseidmva&&lead_loosephoton_id80&&ggh4_loose, procs_dy, ops_linlumi)
  //                .Weight(weight_reg*w_years).Tag("zgfake_est");
  //pm.Push<Hist1D>(Axis(80,100.0,180.0, llloosephoton_m, "m_{ll#gamma} [GeV]", {}), 
  //                january_baseline_looseidmva&&lead_loosephoton_idmva<-0.58&&ggh1_fake, procs_dy, ops_linlumi)
  //                .Weight(weight_reg*w_extrap_ptetares*w_years).Tag("zgfake_est");
  //pm.Push<Hist1D>(Axis(80,100.0,180.0, llloosephoton_m, "m_{ll#gamma} [GeV]", {}), 
  //                january_baseline_looseidmva&&lead_loosephoton_idmva<-0.58&&ggh2_fake, procs_dy, ops_linlumi)
  //                .Weight(weight_reg*w_extrap_ptetares*w_years).Tag("zgfake_est");
  //pm.Push<Hist1D>(Axis(80,100.0,180.0, llloosephoton_m, "m_{ll#gamma} [GeV]", {}), 
  //                january_baseline_looseidmva&&lead_loosephoton_idmva<-0.58&&ggh3_fake, procs_dy, ops_linlumi)
  //                .Weight(weight_reg*w_extrap_ptetares*w_years).Tag("zgfake_est");
  //pm.Push<Hist1D>(Axis(80,100.0,180.0, llloosephoton_m, "m_{ll#gamma} [GeV]", {}), 
  //                january_baseline_looseidmva&&lead_loosephoton_idmva<-0.58&&ggh4_fake, procs_dy, ops_linlumi)
  //                .Weight(weight_reg*w_extrap_ptetares*w_years).Tag("zgfake_est");
  ////check variables to examine method closure
  //std::vector<Axis> plot_vars = {Axis(80,0.0,4.0, lead_loosephoton_drmin, "Photon #Delta R_{min}(l)"),
  //                               Axis(80,0.0,8.0, lead_loosephoton_drmax, "Photon #Delta R_{max}(l)"),
  //                               Axis(80,0.0,2.0, llloosephoton_pt/llloosephoton_m, "p_{T}^{ll#gamma}/m_{ll#gamma}"),
  //                               Axis(80,-1.0,1.0, llloosephoton_cosTheta, "cos(#Theta)"),
  //                               Axis(80,-1.0,1.0, llloosephoton_costheta, "cos(#theta)"),
  //                               Axis(80,-3.1416,3.1416, llloosephoton_psi, "#phi"),
  //                               Axis(80,0.0,0.5, lead_loosephoton_relpterr, "Photon #sigma_{E}/E"),
  //                               Axis(80,15.0,50.0, lead_loosephoton_pt, "Photon p_{T} [GeV]"),
  //                               Axis(80,-2.5,2.5, lead_loosephoton_eta, "Photon #eta"),
  //                               Axis(80,-2.5,2.5, lead_lepton_eta, "Lead lepton #eta"),
  //                               Axis(80,-2.5,2.5, sublead_lepton_eta, "Lead lepton #eta")};
  //for (Axis& var_axis : plot_vars) {
  //  pm.Push<Hist1D>(var_axis,
  //                  january_baseline_looseidmva&&lead_loosephoton_idmva<-0.58, procs_dy, ops_shapes)
  //                  .Weight(weight_reg*w_extrap_ptetares*w_years).Tag("zgfake_est");
  //  pm.Push<Hist1D>(var_axis,
  //                  january_baseline_looseidmva&&lead_loosephoton_id80, procs_dy, ops_shapes)
  //                  .Weight(weight_reg*w_years).Tag("zgfake_est");
  //  pm.Push<Hist1D>(var_axis,
  //                  january_baseline_looseidmva&&lead_loosephoton_idmva<-0.58&&ggh4_fake, procs_dy, ops_shapes)
  //                  .Weight(weight_reg*w_extrap_ptetares*w_years).Tag("zgfake_est");
  //  pm.Push<Hist1D>(var_axis,
  //                  january_baseline_looseidmva&&lead_loosephoton_id80&&ggh4_loose, procs_dy, ops_shapes)
  //                  .Weight(weight_reg*w_years).Tag("zgfake_est");
  //  //pm.Push<Hist1D>(var_axis,
  //  //                january_baseline_looseidmva&&lead_loosephoton_idmva<-0.58&&lead_loosephoton_abseta<1.5, procs_dy, ops_shapes)
  //  //                .Weight(weight_reg*w_extrap_ptetares*w_years).Tag("zgfake_est");
  //  //pm.Push<Hist1D>(var_axis,
  //  //                january_baseline_looseidmva&&lead_loosephoton_id80&&lead_loosephoton_abseta<1.5, procs_dy, ops_shapes)
  //  //                .Weight(weight_reg*w_years).Tag("zgfake_est");
  //  //pm.Push<Hist1D>(var_axis,
  //  //                january_baseline_looseidmva&&lead_loosephoton_idmva<-0.58&&lead_loosephoton_abseta>1.5, procs_dy, ops_shapes)
  //  //                .Weight(weight_reg*w_extrap_ptetares*w_years).Tag("zgfake_est");
  //  //pm.Push<Hist1D>(var_axis,
  //  //                january_baseline_looseidmva&&lead_loosephoton_id80&&lead_loosephoton_abseta>1.5, procs_dy, ops_shapes)
  //  //                .Weight(weight_reg*w_years).Tag("zgfake_est");
  //}
  
  ////check photon IDMVA distributions as a function of pt and eta for reweighting
  //const std::vector<double> pt_bin_boundaries = {15.0,17.5,20.0,25.0,30.0,50.0};
  //const std::vector<double> abseta_bin_boundaries = {0.0,0.8,1.5,2.0,2.5};
  //for (unsigned ipt = 0; ipt < pt_bin_boundaries.size()-1; ipt++) {
  //  for (unsigned ieta = 0; ieta < abseta_bin_boundaries.size()-1; ieta++) {
  //    NamedFunc pt_eta_bin = NamedFunc(lead_loosephoton_pt>pt_bin_boundaries[ipt]
  //                                     &&lead_loosephoton_pt<pt_bin_boundaries[ipt+1]
  //                                     &&lead_loosephoton_abseta>abseta_bin_boundaries[ieta]
  //                                     &&lead_loosephoton_abseta<abseta_bin_boundaries[ieta+1])
  //                                     .Name("ptetabin_"+to_string(ipt)+"_"+to_string(ieta));
  //    pm.Push<Hist1D>(Axis(60,-1.0,1.0, lead_loosephoton_idmva, "Photon IDMVA", {}), 
  //                    january_baseline_looseidmva&&pt_eta_bin, procs_dy, ops_shapes)
  //                    .Weight(weight_reg*w_years).Tag("zgfake_corrfactor");
  //  }
  //}
 
  ////pt-eta-cosTheta weight derivations
  //const std::vector<double> pt_bin_boundaries = {15.0,17.5,20.0,25.0,30.0,50.0};
  //const std::vector<double> abseta_bin_boundaries = {0.0,0.8,1.5,2.0,2.5};
  //const std::vector<double> abscostheta_bin_boundaries = {0.0,0.2,0.4,0.6,0.8,0.9,1.0};
  //for (unsigned ipt = 0; ipt < pt_bin_boundaries.size()-1; ipt++) {
  //  for (unsigned ieta = 0; ieta < abseta_bin_boundaries.size()-1; ieta++) {
  //    for (unsigned ict = 0; ict < abscostheta_bin_boundaries.size()-1; ict++) {
  //      NamedFunc pt_eta_ct_bin = NamedFunc(lead_loosephoton_pt>pt_bin_boundaries[ipt]
  //                                          &&lead_loosephoton_pt<pt_bin_boundaries[ipt+1]
  //                                          &&lead_loosephoton_abseta>abseta_bin_boundaries[ieta]
  //                                          &&lead_loosephoton_abseta<abseta_bin_boundaries[ieta+1]
  //                                          &&llloosephoton_costheta>abscostheta_bin_boundaries[ict]
  //                                          &&llloosephoton_costheta<abscostheta_bin_boundaries[ict+1])
  //                                          .Name("ptetactbin_"+to_string(ipt)+"_"+to_string(ieta)+"_"+to_string(ict));
  //      pm.Push<Hist1D>(Axis(60,-1.0,1.0, lead_loosephoton_idmva, "Photon IDMVA", {}), 
  //                      january_baseline_looseidmva&&pt_eta_ct_bin, procs_dy, ops_shapes)
  //                      .Weight(weight_reg*w_years).Tag("zgfake_corrfactor");
  //    }
  //  }
  //}
  
  ////pt-eta-res weight derivations
  //const std::vector<double> pt_bin_boundaries = {15.0,17.5,20.0,25.0,30.0,50.0};
  //const std::vector<double> abseta_bin_boundaries = {0.0,0.8,1.5,2.0,2.5};
  //const std::vector<double> res_bin_boundaries = {0.00,0.02,0.04,0.08,0.15,2.0};
  //for (unsigned ipt = 0; ipt < pt_bin_boundaries.size()-1; ipt++) {
  //  for (unsigned ieta = 0; ieta < abseta_bin_boundaries.size()-1; ieta++) {
  //    for (unsigned ires = 0; ires < res_bin_boundaries.size()-1; ires++) {
  //      NamedFunc pt_eta_ct_bin = NamedFunc(lead_loosephoton_pt>pt_bin_boundaries[ipt]
  //                                          &&lead_loosephoton_pt<pt_bin_boundaries[ipt+1]
  //                                          &&lead_loosephoton_abseta>abseta_bin_boundaries[ieta]
  //                                          &&lead_loosephoton_abseta<abseta_bin_boundaries[ieta+1]
  //                                          &&lead_loosephoton_relpterr>res_bin_boundaries[ires]
  //                                          &&lead_loosephoton_relpterr<res_bin_boundaries[ires+1])
  //                                          .Name("ptetaresbin_"+to_string(ipt)+"_"+to_string(ieta)+"_"+to_string(ires));
  //      pm.Push<Hist1D>(Axis(60,-1.0,1.0, lead_loosephoton_idmva, "Photon IDMVA", {}), 
  //                      january_baseline_looseidmva&&pt_eta_ct_bin, procs_dy, ops_shapes)
  //                      .Weight(weight_reg*w_years).Tag("zgfake_corrfactor");
  //    }
  //  }
  //}
  
  //check stats
  pm.Push<Table>("sample_stats", vector<TableRow>{
    TableRow("Untagged categories (unweighted)", 
        january_baseline&&gghall,0,0,"1"),
    TableRow("Untagged categories (weighted)", 
        january_baseline&&gghall,0,0,weight_reg*w_years),
    TableRow("Untagged category extrapolation (unweighted)", 
        january_baseline_looseidmva&&lead_loosephoton_idmva<-0.58&&gghall_fake,0,0,"1"),
    TableRow("Untagged category extrapolation (weighted)", 
        january_baseline_looseidmva&&lead_loosephoton_idmva<-0.58&&gghall_fake,0,0,weight_reg*w_years),
  },procs_samples,false,true,false,false,false,true).Precision(3);

  pm.min_print_ = true;
  pm.SetLuminosityTag(lumi_tag).MakePlots(1.0);

  return 0;
}
