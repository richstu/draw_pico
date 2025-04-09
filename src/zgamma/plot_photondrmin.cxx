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

using namespace std;
using namespace PlotOptTypes;
using namespace NamedFuncUtilities;
using namespace ZgFunctions;
using namespace ZgUtilities;

int main() {

  //--------------------------------------------------------------------------
  //                                    initialization
  //--------------------------------------------------------------------------

  //setup
  gErrorIgnoreLevel = 6000;
  string lumi_tag = "138";
  set<string> years = {"2016APV","2016","2017","2018"};
  //string lumi_tag = "60";

  vector<shared_ptr<Process>> procs_llskim_reduced = ZgSampleLoader()
      .SetMacro("YEARS",years)
      .LoadSamples("txt/samples_zgamma.txt","llcr");

  vector<shared_ptr<Process>> procs_llskim_reduced_mc = ZgSampleLoader()
      .SetMacro("YEARS",years)
      .LoadMCSamples("txt/samples_zgamma.txt","llcr");

  vector<shared_ptr<Process>> procs_llskim_inchmumu = ZgSampleLoader()
      .SetMacro("YEARS",years)
      .LoadMCSamples("txt/samples_zgamma.txt","HMuMullskim");
  procs_llskim_inchmumu.insert(procs_llskim_inchmumu.end(), 
                               procs_llskim_reduced_mc.begin(), 
                               procs_llskim_reduced_mc.end());

  std::vector<PlotOpt> ops_data = {PlotOpt("txt/plot_styles.txt","LinLumiDataRoot")}; 
  std::vector<PlotOpt> ops_shapes = {PlotOpt("txt/plot_styles.txt","LinLumi").Overflow(OverflowType::none),
                                     PlotOpt("txt/plot_styles.txt","Shapes").Overflow(OverflowType::none)}; 

  //------------------------------------------------------------------------------------
  //                                   NamedFuncs
  //------------------------------------------------------------------------------------

  //weight with some clipping/regularization
  const NamedFunc weight_reg("weight_reg",[](const Baby &b) -> NamedFunc::ScalarType{
    if (fabs(b.weight())>5) return 0.0;
    return b.weight();
  });

  const NamedFunc photon_drsig = NamedFunc("(photon_isScEtaEB||photon_isScEtaEE)&&photon_id80&&photon_elveto&&(photon_pt>15)").Name("photon_drsig"); 

  const NamedFunc ndrphoton = ReduceNamedFunc(photon_drsig, reduce_sum).Name("ndrphoton");
  const NamedFunc drphoton_pt = FilterNamedFunc("photon_pt",photon_drsig).Name("drphoton_pt");
  const NamedFunc drphoton_eta = FilterNamedFunc("photon_eta",photon_drsig).Name("drphoton_eta");
  const NamedFunc drphoton_phi = FilterNamedFunc("photon_phi",photon_drsig).Name("drphoton_phi");
  const NamedFunc drphoton_drmin = FilterNamedFunc("photon_drmin",photon_drsig).Name("drphoton_drmin");
  const NamedFunc drphoton_elidx = FilterNamedFunc("photon_elidx",photon_drsig).Name("drphoton_elidx");
  const NamedFunc drphoton_drmax = FilterNamedFunc("photon_drmax",photon_drsig).Name("drphoton_drmax");
  const NamedFunc drphoton_idmva = FilterNamedFunc("photon_idmva",photon_drsig).Name("drphoton_idmva");
  const NamedFunc lead_drphoton_pt = ReduceNamedFunc(drphoton_pt,reduce_max).Name("lead_drphoton_pt");
  const NamedFunc lead_drphoton_drmin = MultiReduceNamedFunc(
      {drphoton_pt,drphoton_drmin},reduce_maxfirst).Name("lead_drphoton_drmin");
  const NamedFunc lead_drphoton_eta = MultiReduceNamedFunc(
      {drphoton_pt,drphoton_eta},reduce_maxfirst).Name("lead_drphoton_eta");
  const NamedFunc lead_drphoton_phi = MultiReduceNamedFunc(
      {drphoton_pt,drphoton_phi},reduce_maxfirst).Name("lead_drphoton_phi");
  const NamedFunc lead_drphoton_elidx = MultiReduceNamedFunc(
      {drphoton_pt,drphoton_elidx},reduce_maxfirst).Name("lead_drphoton_elidx");
  const NamedFunc lead_drphoton_drmax = MultiReduceNamedFunc(
      {drphoton_pt,drphoton_drmax},reduce_maxfirst).Name("lead_drphoton_drmax");
  const NamedFunc lead_drphoton_idmva = MultiReduceNamedFunc(
      {drphoton_pt,drphoton_idmva},reduce_maxfirst).Name("lead_drphoton_idmva");

  //relative pt uncertainty of lead photon with dr requirement loosened
  const NamedFunc lead_drphoton_relpterr("lead_drphoton_relpterr",[photon_drsig](const Baby &b) -> NamedFunc::ScalarType{
    vector<double> sig = photon_drsig.GetVector(b);
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

  //phi distance to nearest lepton of lead photon with dr requirement loosened
  const NamedFunc lead_drphoton_dphimin("lead_drphoton_dphimin",[photon_drsig](const Baby &b) -> NamedFunc::ScalarType{
    vector<double> sig = photon_drsig.GetVector(b);
    float max_pt = 0;
    unsigned max_idx = 0;
    for (unsigned iph = 0; iph < b.photon_pt()->size(); iph++) {
      if (sig[iph] && b.photon_pt()->at(iph) > max_pt) {
        max_pt = b.photon_pt()->at(iph);
        max_idx = iph;
      }
    }
    float dphimin = 999.0;
    for (unsigned iel = 0; iel < b.el_pt()->size(); iel++) {
      if (b.el_sig()->at(iel)) {
        float this_dphi = deltaPhi(b.photon_phi()->at(max_idx), b.el_phi()->at(iel));
        if (this_dphi < dphimin)
          dphimin = this_dphi;
      }
    }
    for (unsigned imu = 0; imu < b.mu_pt()->size(); imu++) {
      if (b.mu_sig()->at(imu)) {
        float this_dphi = deltaPhi(b.photon_phi()->at(max_idx), b.mu_phi()->at(imu));
        if (this_dphi < dphimin)
          dphimin = this_dphi;
      }
    }
    return dphimin;
  });

  //eta distance to nearest lepton of lead photon with dr requirement loosened
  const NamedFunc lead_drphoton_detamin("lead_drphoton_detamin",[photon_drsig](const Baby &b) -> NamedFunc::ScalarType{
    vector<double> sig = photon_drsig.GetVector(b);
    float max_pt = 0;
    unsigned max_idx = 0;
    for (unsigned iph = 0; iph < b.photon_pt()->size(); iph++) {
      if (sig[iph] && b.photon_pt()->at(iph) > max_pt) {
        max_pt = b.photon_pt()->at(iph);
        max_idx = iph;
      }
    }
    float detamin = 999.0;
    for (unsigned iel = 0; iel < b.el_pt()->size(); iel++) {
      if (b.el_sig()->at(iel)) {
        float this_deta = b.photon_eta()->at(max_idx)-b.el_eta()->at(iel);
        if (fabs(this_deta) < detamin)
          detamin = this_deta;
      }
    }
    for (unsigned imu = 0; imu < b.mu_pt()->size(); imu++) {
      if (b.mu_sig()->at(imu)) {
        float this_deta = b.photon_eta()->at(max_idx)-b.mu_eta()->at(imu);
        if (fabs(this_deta) < detamin)
          detamin = this_deta;
      }
    }
    return detamin;
  });

  //llphoton mass with photon dr reqiurement loosened
  const NamedFunc lldrphoton_m("lldrphoton_m",[lead_drphoton_pt, lead_drphoton_eta, lead_drphoton_phi](const Baby &b) -> NamedFunc::ScalarType{
    TLorentzVector z, ph;
    z.SetPtEtaPhiM(b.ll_pt()->at(0),b.ll_eta()->at(0),b.ll_phi()->at(0),b.ll_m()->at(0));
    ph.SetPtEtaPhiM(lead_drphoton_pt.GetScalar(b), lead_drphoton_eta.GetScalar(b),
                    lead_drphoton_phi.GetScalar(b), 0.0);
    return (z+ph).M();
  });

  //llphoton pt with photon dr reqiurement loosened
  const NamedFunc lldrphoton_pt("lldrphoton_pt",[lead_drphoton_pt, lead_drphoton_eta, lead_drphoton_phi](const Baby &b) -> NamedFunc::ScalarType{
    TLorentzVector z, ph;
    z.SetPtEtaPhiM(b.ll_pt()->at(0),b.ll_eta()->at(0),b.ll_phi()->at(0),b.ll_m()->at(0));
    ph.SetPtEtaPhiM(lead_drphoton_pt.GetScalar(b), lead_drphoton_eta.GetScalar(b),
                    lead_drphoton_phi.GetScalar(b), 0.0);
    return (z+ph).Pt();
  });

  //llphoton costheta with photon dr reqiurement loosened
  const NamedFunc lldrphoton_costheta("lldrphoton_costheta",[lead_drphoton_pt, lead_drphoton_eta, lead_drphoton_phi](const Baby &b) -> NamedFunc::ScalarType{
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
    ph.SetPtEtaPhiM(lead_drphoton_pt.GetScalar(b), lead_drphoton_eta.GetScalar(b),
                    lead_drphoton_phi.GetScalar(b), 0.0);
    return get_costheta(lm,lp,ph);
  });

  //llphoton cosTheta with photon dr reqiurement loosened
  const NamedFunc lldrphoton_cosTheta("lldrphoton_cosTheta",[lead_drphoton_pt, lead_drphoton_eta, lead_drphoton_phi](const Baby &b) -> NamedFunc::ScalarType{
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
    ph.SetPtEtaPhiM(lead_drphoton_pt.GetScalar(b), lead_drphoton_eta.GetScalar(b),
                    lead_drphoton_phi.GetScalar(b), 0.0);
    return get_cosTheta(lm,lp,ph);
  });

  //llphoton phi with photon dr reqiurement loosened
  const NamedFunc lldrphoton_psi("lldrphoton_psi",[lead_drphoton_pt, lead_drphoton_eta, lead_drphoton_phi](const Baby &b) -> NamedFunc::ScalarType{
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
    ph.SetPtEtaPhiM(lead_drphoton_pt.GetScalar(b), lead_drphoton_eta.GetScalar(b),
                    lead_drphoton_phi.GetScalar(b), 0.0);
    return get_phi(lm,lp,ph);
  });

  const NamedFunc blind_sr_data_nodr = NamedFunc("type>=1000"||(lldrphoton_m<120||lldrphoton_m>130)).Name("blind_sr_data"); //n.b. doesn't work for run 3

  const NamedFunc blind_sr_nodr = NamedFunc((lldrphoton_m<120||lldrphoton_m>130)).Name("blind_sr");

  const NamedFunc january_baseline = NamedFunc("(nel>=2||nmu>=2) && nphoton>=1 && photon_id80[0] && (ll_m[0]>81 && ll_m[0]<101) && (ll_charge[0]==0) && ((photon_pt[0]/llphoton_m[0])>15.0/110.0) && ((llphoton_m[0]+ll_m[0])>185)" && trig_plateau_cuts).Name("baseline");

  const NamedFunc january_baseline_nodr = NamedFunc("(nel>=2||nmu>=2)" && (ndrphoton>=1) && "(ll_m[0]>81 && ll_m[0]<101) && (ll_charge[0]==0) " && ((lead_drphoton_pt/lldrphoton_m)>15.0/110.0) && ((lldrphoton_m+"ll_m[0]")>185) && trig_plateau_cuts ).Name("baseline");

  const NamedFunc lead_drphoton_elsig("lead_drphoton_elsig",[lead_drphoton_elidx](const Baby &b) -> NamedFunc::ScalarType{
    int el_idx = static_cast<int>(lead_drphoton_elidx.GetScalar(b));
    if (el_idx == -1) return false;
    return b.el_sig()->at(el_idx);
  });

  const NamedFunc blind_sr_data = NamedFunc("type>=1000||(llphoton_m[0]<120||llphoton_m[0]>130)").Name("blind_sr_data"); //n.b. doesn't work for run 3

  const NamedFunc blind_sr = NamedFunc("(llphoton_m[0]<120||llphoton_m[0]>130)").Name("blind_sr"); //n.b. doesn't work for run 3

  shared_ptr<MVAWrapper> kinematic_bdt = KinematicBdt();
  NamedFunc ggh4 = category_ggh4_old(kinematic_bdt);
  NamedFunc ggh3 = category_ggh3_old(kinematic_bdt);
  NamedFunc ggh2 = category_ggh2_old(kinematic_bdt);
  NamedFunc ggh1 = category_ggh1_old(kinematic_bdt);

  std::shared_ptr<MVAWrapper> drloose_bdt = std::make_shared<MVAWrapper>("drloose_bdt");
  drloose_bdt->SetVariable("photon_mva",lead_drphoton_idmva);
  drloose_bdt->SetVariable("min_dR",lead_drphoton_drmin);
  drloose_bdt->SetVariable("max_dR",lead_drphoton_drmax);
  drloose_bdt->SetVariable("pt_mass",lldrphoton_pt/lldrphoton_m);
  drloose_bdt->SetVariable("cosTheta",lldrphoton_cosTheta);
  drloose_bdt->SetVariable("costheta",lldrphoton_costheta);
  drloose_bdt->SetVariable("phi",lldrphoton_psi);
  drloose_bdt->SetVariable("photon_res",lead_drphoton_relpterr);
  drloose_bdt->SetVariable("photon_rapidity",lead_drphoton_eta);
  drloose_bdt->SetVariable("l1_rapidity",ZgFunctions::lead_lepton_eta);
  drloose_bdt->SetVariable("l2_rapidity",ZgFunctions::sublead_lepton_eta);
  drloose_bdt->BookMVA("/homes/oshiro/analysis/small_phys_utils/dataset/weights/shuffled_phidcomp_post_phidcomp_post_BDT.weights.xml");
  NamedFunc drloose_bdt_score = drloose_bdt->GetDiscriminant();
  NamedFunc ggh4_drloose = category_ggh4_old(drloose_bdt);
  NamedFunc ggh3_drloose = category_ggh3_old(drloose_bdt);
  NamedFunc ggh2_drloose = category_ggh2_old(drloose_bdt);
  NamedFunc ggh1_drloose = category_ggh1_old(drloose_bdt);
  
  //------------------------------------------------------------------------------------
  //                                   plots and tables
  //------------------------------------------------------------------------------------

  PlotMaker pm;
  pm.multithreaded_ = false; //using BDT

  ////dR plots to motivate loosening cuts/show danger
  pm.Push<Hist1D>(Axis(70,0.0,3.5, lead_drphoton_drmin, "#Delta R_{min}(l,#gamma) [GeV]", {}), 
      january_baseline_nodr&&"ll_lepid[0]==11"&&!lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130, 
      procs_llskim_reduced_mc, ops_shapes).Weight(weight_reg*w_years).Tag("zgdr");
  pm.Push<Hist1D>(Axis(70,0.0,3.5, lead_drphoton_drmin, "#Delta R_{min}(l,#gamma) [GeV]", {}), 
      january_baseline_nodr&&"ll_lepid[0]==13"&&!lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130, 
      procs_llskim_reduced_mc, ops_shapes).Weight(weight_reg*w_years).Tag("zgdr");
  pm.Push<Hist1D>(Axis(70,0.0,3.5, lead_drphoton_drmin, "#Delta R_{min}(l,#gamma) [GeV]", {}), 
      january_baseline_nodr&&"ll_lepid[0]==11"&&lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130, 
      procs_llskim_reduced_mc, ops_shapes).Weight(weight_reg*w_years).Tag("zgdr");
  pm.Push<Hist1D>(Axis(70,0.0,3.5, lead_drphoton_drmin, "#Delta R_{min}(l,#gamma) [GeV]", {}), 
      january_baseline_nodr&&"ll_lepid[0]==13"&&lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130, 
      procs_llskim_reduced_mc, ops_shapes).Weight(weight_reg*w_years).Tag("zgdr");
  pm.Push<Hist1D>(Axis(70,0.0,3.5, lead_drphoton_drmin, "#Delta R_{min}(l,#gamma) [GeV]", {}), 
      january_baseline_nodr&&"ll_lepid[0]==11"&&!lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130, 
      procs_llskim_inchmumu, ops_shapes).Weight(weight_reg*w_years).Tag("zgdr_mm");
  pm.Push<Hist1D>(Axis(70,0.0,3.5, lead_drphoton_drmin, "#Delta R_{min}(l,#gamma) [GeV]", {}), 
      january_baseline_nodr&&"ll_lepid[0]==13"&&!lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130, 
      procs_llskim_inchmumu, ops_shapes).Weight(weight_reg*w_years).Tag("zgdr_mm");
  pm.Push<Hist1D>(Axis(70,0.0,3.5, lead_drphoton_drmin, "#Delta R_{min}(l,#gamma) [GeV]", {}), 
      january_baseline_nodr&&"ll_lepid[0]==11"&&lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130, 
      procs_llskim_inchmumu, ops_shapes).Weight(weight_reg*w_years).Tag("zgdr_mm");
  pm.Push<Hist1D>(Axis(70,0.0,3.5, lead_drphoton_drmin, "#Delta R_{min}(l,#gamma) [GeV]", {}), 
      january_baseline_nodr&&"ll_lepid[0]==13"&&lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130, 
      procs_llskim_inchmumu, ops_shapes).Weight(weight_reg*w_years).Tag("zgdr_mm");
  ////dphi, deta, and lepton flavor
  //pm.Push<Hist1D>(Axis(70,0.0,3.5, lead_drphoton_drmin, "#Delta R_{min}(l,#gamma) [GeV]", {}), 
  //    january_baseline_nodr&&!lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&"nel>=2", procs_llskim_reduced_mc, ops_shapes)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(70,0.0,3.5, lead_drphoton_drmin, "#Delta R_{min}(l,#gamma) [GeV]", {}), 
  //    january_baseline_nodr&&lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&"nel>=2", procs_llskim_reduced_mc, ops_shapes)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(70,0.0,3.5, lead_drphoton_drmin, "#Delta R_{min}(l,#gamma) [GeV]", {}), 
  //    january_baseline_nodr&&!lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&"nmu>=2", procs_llskim_reduced_mc, ops_shapes)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(70,0.0,3.5, lead_drphoton_drmin, "#Delta R_{min}(l,#gamma) [GeV]", {}), 
  //    january_baseline_nodr&&lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&"nmu>=2", procs_llskim_reduced_mc, ops_shapes)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&!lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>-0.2&&lead_drphoton_detamin<-0.15&&"nel>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>-0.2&&lead_drphoton_detamin<-0.15&&"nel>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&!lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>-0.15&&lead_drphoton_detamin<-0.1&&"nel>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>-0.15&&lead_drphoton_detamin<-0.1&&"nel>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&!lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>-0.1&&lead_drphoton_detamin<-0.05&&"nel>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>-0.1&&lead_drphoton_detamin<-0.05&&"nel>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&!lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>-0.05&&lead_drphoton_detamin<0.0&&"nel>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>-0.05&&lead_drphoton_detamin<0.0&&"nel>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&!lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>0.0&&lead_drphoton_detamin<0.05&&"nel>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>0.0&&lead_drphoton_detamin<0.05&&"nel>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&!lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>0.05&&lead_drphoton_detamin<0.1&&"nel>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>0.05&&lead_drphoton_detamin<0.1&&"nel>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&!lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>0.1&&lead_drphoton_detamin<0.15&&"nel>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>0.1&&lead_drphoton_detamin<0.15&&"nel>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&!lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>0.15&&lead_drphoton_detamin<0.2&&"nel>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>0.15&&lead_drphoton_detamin<0.2&&"nel>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&!lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>-0.2&&lead_drphoton_detamin<-0.15&&"nmu>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>-0.2&&lead_drphoton_detamin<-0.15&&"nmu>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&!lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>-0.15&&lead_drphoton_detamin<-0.1&&"nmu>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>-0.15&&lead_drphoton_detamin<-0.1&&"nmu>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&!lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>-0.1&&lead_drphoton_detamin<-0.05&&"nmu>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>-0.1&&lead_drphoton_detamin<-0.05&&"nmu>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&!lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>-0.05&&lead_drphoton_detamin<0.0&&"nmu>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>-0.05&&lead_drphoton_detamin<0.0&&"nmu>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&!lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>0.0&&lead_drphoton_detamin<0.05&&"nmu>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>0.0&&lead_drphoton_detamin<0.05&&"nmu>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&!lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>0.05&&lead_drphoton_detamin<0.1&&"nmu>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>0.05&&lead_drphoton_detamin<0.1&&"nmu>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&!lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>0.1&&lead_drphoton_detamin<0.15&&"nmu>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>0.1&&lead_drphoton_detamin<0.15&&"nmu>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&!lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>0.15&&lead_drphoton_detamin<0.2&&"nmu>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////pm.Push<Hist1D>(Axis(70,0.0,3.0, lead_drphoton_dphimin, "#Delta #phi [GeV]", {}), 
  ////    january_baseline_nodr&&lead_drphoton_elsig&&lldrphoton_m>120&&lldrphoton_m<130&&lead_drphoton_detamin>0.15&&lead_drphoton_detamin<0.2&&"nmu>=2", procs_llskim_reduced_mc, ops_shapes)
  ////    .Weight(weight_reg*w_years).Tag("zgdr");
  ////mllg plots to further quantify danger (cf. Anders slides)
  //pm.Push<Hist1D>(Axis(50,100,150, lldrphoton_m, "m_{ll#gamma} [GeV]", {}), 
  //    january_baseline_nodr&&!lead_drphoton_elsig&&blind_sr_data_nodr, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,100,150, lldrphoton_m, "m_{ll#gamma} [GeV]", {}), 
  //    january_baseline_nodr&&lead_drphoton_elsig&&blind_sr_data_nodr, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,100,150, lldrphoton_m, "m_{ll#gamma} [GeV]", {}), 
  //    january_baseline_nodr&&!lead_drphoton_elsig&&blind_sr_data_nodr&&lead_drphoton_drmin<0.1, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,100,150, lldrphoton_m, "m_{ll#gamma} [GeV]", {}), 
  //    january_baseline_nodr&&lead_drphoton_elsig&&blind_sr_data_nodr&&lead_drphoton_drmin<0.1, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,100,150, lldrphoton_m, "m_{ll#gamma} [GeV]", {}), 
  //    january_baseline_nodr&&!lead_drphoton_elsig&&blind_sr_data_nodr&&lead_drphoton_drmin>0.1&&lead_drphoton_drmin<0.2, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,100,150, lldrphoton_m, "m_{ll#gamma} [GeV]", {}), 
  //    january_baseline_nodr&&lead_drphoton_elsig&&blind_sr_data_nodr&&lead_drphoton_drmin>0.1&&lead_drphoton_drmin<0.2, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,100,150, lldrphoton_m, "m_{ll#gamma} [GeV]", {}), 
  //    january_baseline_nodr&&!lead_drphoton_elsig&&blind_sr_data_nodr&&lead_drphoton_drmin>0.2&&lead_drphoton_drmin<0.3, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,100,150, lldrphoton_m, "m_{ll#gamma} [GeV]", {}), 
  //    january_baseline_nodr&&lead_drphoton_elsig&&blind_sr_data_nodr&&lead_drphoton_drmin>0.2&&lead_drphoton_drmin<0.3, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,100,150, lldrphoton_m, "m_{ll#gamma} [GeV]", {}), 
  //    january_baseline_nodr&&!lead_drphoton_elsig&&blind_sr_data_nodr&&lead_drphoton_drmin>0.3&&lead_drphoton_drmin<0.4, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,100,150, lldrphoton_m, "m_{ll#gamma} [GeV]", {}), 
  //    january_baseline_nodr&&lead_drphoton_elsig&&blind_sr_data_nodr&&lead_drphoton_drmin>0.3&&lead_drphoton_drmin<0.4, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,100,150, lldrphoton_m, "m_{ll#gamma} [GeV]", {}), 
  //    january_baseline_nodr&&!lead_drphoton_elsig&&blind_sr_data_nodr&&lead_drphoton_drmin>0.4, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,100,150, lldrphoton_m, "m_{ll#gamma} [GeV]", {}), 
  //    january_baseline_nodr&&lead_drphoton_elsig&&blind_sr_data_nodr&&lead_drphoton_drmin>0.4, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  ////check effect on analysis variables ver 2
  //NamedFunc drcut = NamedFunc(lead_drphoton_drmin>0.3&&lead_drphoton_drmin<0.4).Name("drcut");
  //pm.Push<Hist1D>(Axis(70,100,170, lldrphoton_m, "m_{ll#gamma} [GeV]", {}), 
  //    january_baseline_nodr&&!lead_drphoton_elsig&&blind_sr_nodr&&drcut, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,-1.0,1.0, lead_drphoton_idmva, "Photon IDMVA", {}), 
  //    january_baseline_nodr&&!lead_drphoton_elsig&&blind_sr_nodr&&drcut, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,0.0,3.5, lead_drphoton_drmin, "Photon #Delta R_{min}(l)", {}), 
  //    january_baseline_nodr&&!lead_drphoton_elsig&&blind_sr_nodr&&drcut, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,0.0,5.0, lead_drphoton_drmax, "Photon #Delta R_{max}(l)", {}), 
  //    january_baseline_nodr&&!lead_drphoton_elsig&&blind_sr_nodr&&drcut, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,0.0,2.5, (lldrphoton_pt/lldrphoton_m), "Higgs candidate p_{T}/m", {}), 
  //    january_baseline_nodr&&!lead_drphoton_elsig&&blind_sr_nodr&&drcut, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,-1.0,1.0, lldrphoton_cosTheta, "cos(#Theta)", {}), 
  //    january_baseline_nodr&&!lead_drphoton_elsig&&blind_sr_nodr&&drcut, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,-1.0,1.0, lldrphoton_costheta, "cos(#theta)", {}), 
  //    january_baseline_nodr&&!lead_drphoton_elsig&&blind_sr_nodr&&drcut, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,0.0,2.0*3.1416, lldrphoton_psi, "#phi", {}), 
  //    january_baseline_nodr&&!lead_drphoton_elsig&&blind_sr_nodr&&drcut, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,0.0,0.2, lead_drphoton_relpterr, "Photon #sigma_{pT}/p_{T}", {}), 
  //    january_baseline_nodr&&!lead_drphoton_elsig&&blind_sr_nodr&&drcut, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,-2.5,2.5, lead_drphoton_eta, "Photon #eta", {}), 
  //    january_baseline_nodr&&!lead_drphoton_elsig&&blind_sr_nodr&&drcut, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,-2.5,2.5, ZgFunctions::lead_lepton_eta, "Lead lepton #eta", {}), 
  //    january_baseline_nodr&&!lead_drphoton_elsig&&blind_sr_nodr&&drcut, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,-2.5,2.5, ZgFunctions::sublead_lepton_eta, "Sublead lepton #eta", {}), 
  //    january_baseline_nodr&&!lead_drphoton_elsig&&blind_sr_nodr&&drcut, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(70,-1.0,1.0, drloose_bdt_score, "Kinematic BDT score", {}), 
  //    january_baseline_nodr&&!lead_drphoton_elsig&&blind_sr_nodr&&drcut, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  ////For toy bias studies
  //pm.Push<Hist1D>(Axis(50,100,150, lldrphoton_m, "m_{ll#gamma} [GeV]", {}), 
  //    january_baseline_nodr&&!lead_drphoton_elsig&&blind_sr_data_nodr&&lead_drphoton_drmin>0.4&&ggh4_drloose, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,100,150, lldrphoton_m, "m_{ll#gamma} [GeV]", {}), 
  //    january_baseline_nodr&&lead_drphoton_elsig&&blind_sr_data_nodr&&lead_drphoton_drmin>0.4&&ggh4_drloose, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,100,150, lldrphoton_m, "m_{ll#gamma} [GeV]", {}), 
  //    january_baseline_nodr&&!lead_drphoton_elsig&&blind_sr_data_nodr&&lead_drphoton_drmin>0.4&&ggh3_drloose, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,100,150, lldrphoton_m, "m_{ll#gamma} [GeV]", {}), 
  //    january_baseline_nodr&&lead_drphoton_elsig&&blind_sr_data_nodr&&lead_drphoton_drmin>0.4&&ggh3_drloose, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,100,150, lldrphoton_m, "m_{ll#gamma} [GeV]", {}), 
  //    january_baseline_nodr&&!lead_drphoton_elsig&&blind_sr_data_nodr&&lead_drphoton_drmin>0.4&&ggh2_drloose, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,100,150, lldrphoton_m, "m_{ll#gamma} [GeV]", {}), 
  //    january_baseline_nodr&&lead_drphoton_elsig&&blind_sr_data_nodr&&lead_drphoton_drmin>0.4&&ggh2_drloose, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,100,150, lldrphoton_m, "m_{ll#gamma} [GeV]", {}), 
  //    january_baseline_nodr&&!lead_drphoton_elsig&&blind_sr_data_nodr&&lead_drphoton_drmin>0.4&&ggh1_drloose, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,100,150, lldrphoton_m, "m_{ll#gamma} [GeV]", {}), 
  //    january_baseline_nodr&&lead_drphoton_elsig&&blind_sr_data_nodr&&lead_drphoton_drmin>0.4&&ggh1_drloose, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,100,150, lldrphoton_m, "m_{ll#gamma} [GeV]", {}), 
  //    january_baseline_nodr&&!lead_drphoton_elsig&&blind_sr_data_nodr&&lead_drphoton_drmin>0.3&&ggh4_drloose, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,100,150, lldrphoton_m, "m_{ll#gamma} [GeV]", {}), 
  //    january_baseline_nodr&&lead_drphoton_elsig&&blind_sr_data_nodr&&lead_drphoton_drmin>0.3&&ggh4_drloose, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,100,150, lldrphoton_m, "m_{ll#gamma} [GeV]", {}), 
  //    january_baseline_nodr&&!lead_drphoton_elsig&&blind_sr_data_nodr&&lead_drphoton_drmin>0.3&&ggh3_drloose, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,100,150, lldrphoton_m, "m_{ll#gamma} [GeV]", {}), 
  //    january_baseline_nodr&&lead_drphoton_elsig&&blind_sr_data_nodr&&lead_drphoton_drmin>0.3&&ggh3_drloose, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,100,150, lldrphoton_m, "m_{ll#gamma} [GeV]", {}), 
  //    january_baseline_nodr&&!lead_drphoton_elsig&&blind_sr_data_nodr&&lead_drphoton_drmin>0.3&&ggh2_drloose, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,100,150, lldrphoton_m, "m_{ll#gamma} [GeV]", {}), 
  //    january_baseline_nodr&&lead_drphoton_elsig&&blind_sr_data_nodr&&lead_drphoton_drmin>0.3&&ggh2_drloose, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,100,150, lldrphoton_m, "m_{ll#gamma} [GeV]", {}), 
  //    january_baseline_nodr&&!lead_drphoton_elsig&&blind_sr_data_nodr&&lead_drphoton_drmin>0.3&&ggh1_drloose, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //pm.Push<Hist1D>(Axis(50,100,150, lldrphoton_m, "m_{ll#gamma} [GeV]", {}), 
  //    january_baseline_nodr&&lead_drphoton_elsig&&blind_sr_data_nodr&&lead_drphoton_drmin>0.3&&ggh1_drloose, procs_llskim_reduced, ops_data)
  //    .Weight(weight_reg*w_years).Tag("zgdr");
  //yields to check sensitivity
  pm.Push<Table>("photondryields", vector<TableRow>{
    TableRow("Yield (dR>0.4 cat1)", 
      january_baseline_nodr&&lldrphoton_m>122&&lldrphoton_m<128&&ggh1_drloose&&lead_drphoton_drmin>0.4,0,0,weight_reg*w_years),
    TableRow("Yield (dR>0.4 cat2)", 
      january_baseline_nodr&&lldrphoton_m>122&&lldrphoton_m<128&&ggh2_drloose&&lead_drphoton_drmin>0.4,0,0,weight_reg*w_years),
    TableRow("Yield (dR>0.4 cat3)", 
      january_baseline_nodr&&lldrphoton_m>122&&lldrphoton_m<128&&ggh3_drloose&&lead_drphoton_drmin>0.4,0,0,weight_reg*w_years),
    TableRow("Yield (dR>0.4 cat4)", 
      january_baseline_nodr&&lldrphoton_m>122&&lldrphoton_m<128&&ggh4_drloose&&lead_drphoton_drmin>0.4,0,0,weight_reg*w_years),
    TableRow("Yield (dR>0.35 cat1)", 
      january_baseline_nodr&&lldrphoton_m>122&&lldrphoton_m<128&&ggh1_drloose&&lead_drphoton_drmin>0.35,0,0,weight_reg*w_years),
    TableRow("Yield (dR>0.35 cat2)", 
      january_baseline_nodr&&lldrphoton_m>122&&lldrphoton_m<128&&ggh2_drloose&&lead_drphoton_drmin>0.35,0,0,weight_reg*w_years),
    TableRow("Yield (dR>0.35 cat3)", 
      january_baseline_nodr&&lldrphoton_m>122&&lldrphoton_m<128&&ggh3_drloose&&lead_drphoton_drmin>0.35,0,0,weight_reg*w_years),
    TableRow("Yield (dR>0.35 cat4)", 
      january_baseline_nodr&&lldrphoton_m>122&&lldrphoton_m<128&&ggh4_drloose&&lead_drphoton_drmin>0.35,0,0,weight_reg*w_years),
    TableRow("Yield (dR>0.3 cat1)", 
      january_baseline_nodr&&lldrphoton_m>122&&lldrphoton_m<128&&ggh1_drloose&&lead_drphoton_drmin>0.3,0,0,weight_reg*w_years),
    TableRow("Yield (dR>0.3 cat2)", 
      january_baseline_nodr&&lldrphoton_m>122&&lldrphoton_m<128&&ggh2_drloose&&lead_drphoton_drmin>0.3,0,0,weight_reg*w_years),
    TableRow("Yield (dR>0.3 cat3)", 
      january_baseline_nodr&&lldrphoton_m>122&&lldrphoton_m<128&&ggh3_drloose&&lead_drphoton_drmin>0.3,0,0,weight_reg*w_years),
    TableRow("Yield (dR>0.3 cat4)", 
      january_baseline_nodr&&lldrphoton_m>122&&lldrphoton_m<128&&ggh4_drloose&&lead_drphoton_drmin>0.3,0,0,weight_reg*w_years),
    TableRow("Yield (dR>0.25 cat1)", 
      january_baseline_nodr&&lldrphoton_m>122&&lldrphoton_m<128&&ggh1_drloose&&lead_drphoton_drmin>0.25,0,0,weight_reg*w_years),
    TableRow("Yield (dR>0.25 cat2)", 
      january_baseline_nodr&&lldrphoton_m>122&&lldrphoton_m<128&&ggh2_drloose&&lead_drphoton_drmin>0.25,0,0,weight_reg*w_years),
    TableRow("Yield (dR>0.25 cat3)", 
      january_baseline_nodr&&lldrphoton_m>122&&lldrphoton_m<128&&ggh3_drloose&&lead_drphoton_drmin>0.25,0,0,weight_reg*w_years),
    TableRow("Yield (dR>0.25 cat4)", 
      january_baseline_nodr&&lldrphoton_m>122&&lldrphoton_m<128&&ggh4_drloose&&lead_drphoton_drmin>0.25,0,0,weight_reg*w_years),
    TableRow("Yield (dR>0.2 cat1)", 
      january_baseline_nodr&&lldrphoton_m>122&&lldrphoton_m<128&&ggh1_drloose&&lead_drphoton_drmin>0.2,0,0,weight_reg*w_years),
    TableRow("Yield (dR>0.2 cat2)", 
      january_baseline_nodr&&lldrphoton_m>122&&lldrphoton_m<128&&ggh2_drloose&&lead_drphoton_drmin>0.2,0,0,weight_reg*w_years),
    TableRow("Yield (dR>0.2 cat3)", 
      january_baseline_nodr&&lldrphoton_m>122&&lldrphoton_m<128&&ggh3_drloose&&lead_drphoton_drmin>0.2,0,0,weight_reg*w_years),
    TableRow("Yield (dR>0.2 cat4)", 
      january_baseline_nodr&&lldrphoton_m>122&&lldrphoton_m<128&&ggh4_drloose&&lead_drphoton_drmin>0.2,0,0,weight_reg*w_years),
    TableRow("Yield (dR>0.15 cat1)", 
      january_baseline_nodr&&lldrphoton_m>122&&lldrphoton_m<128&&ggh1_drloose&&lead_drphoton_drmin>0.15,0,0,weight_reg*w_years),
    TableRow("Yield (dR>0.15 cat2)", 
      january_baseline_nodr&&lldrphoton_m>122&&lldrphoton_m<128&&ggh2_drloose&&lead_drphoton_drmin>0.15,0,0,weight_reg*w_years),
    TableRow("Yield (dR>0.15 cat3)", 
      january_baseline_nodr&&lldrphoton_m>122&&lldrphoton_m<128&&ggh3_drloose&&lead_drphoton_drmin>0.15,0,0,weight_reg*w_years),
    TableRow("Yield (dR>0.15 cat4)", 
      january_baseline_nodr&&lldrphoton_m>122&&lldrphoton_m<128&&ggh4_drloose&&lead_drphoton_drmin>0.15,0,0,weight_reg*w_years),
  },procs_llskim_inchmumu,false,true,false,false,false,true).Precision(3);

  pm.min_print_ = true;
  pm.SetLuminosityTag(lumi_tag).MakePlots(1.0);

  return 0;
}
