#include <iostream>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include "TLorentzVector.h"
#include "TMath.h"
#include "TString.h"

#include "core/baby.hpp"
#include "core/fastforest.hpp"
#include "core/named_func.hpp"
#include "core/named_func_utilities.hpp"
#include "zgamma/KinZfitter.hpp"
#include "zgamma/zg_syst_functions.hpp"
#include "zgamma/zg_utilities.hpp"

using std::lock_guard;
using std::make_shared;
using std::mutex;
using std::shared_ptr;
using std::string;
using std::unordered_map;
using std::vector;
using fastforest::FastForest;
using fastforest::FastForest;
using NamedFuncUtilities::FilterNamedFuncCached;
using NamedFuncUtilities::MultiReduceNamedFuncCached;
using NamedFuncUtilities::MultiMapNamedFuncCached;
using NamedFuncUtilities::ReduceNamedFuncCached;
using NamedFuncUtilities::ReduceNamedFunc;
using NamedFuncUtilities::FilterNamedFunc;
using NamedFuncUtilities::reduce_sum;
using NamedFuncUtilities::reduce_max;
using NamedFuncUtilities::reduce_maxfirst;
using NamedFuncUtilities::reduce_sublead;
using NamedFuncUtilities::reduce_subleadfirst;
using ZgUtilities::CalculateAngles;
using ZgUtilities::get_lep_custom_refit;
using ZgUtilities::get_btag_wp_deepjet;
using ZgUtilities::XGBoostBDTScoreCached;

//Be very careful!! zg_functions MUST be linked before this file so that the
//relevant NamedFuncs are defined in advance. The name of this file is made
//to be alphabetically after zg_functions (and so scons will link it after
//zg_functions), but if the linking order is changed, this file must be linked
//second
//alternatively, all relevant NamedFuncs could be copied over to this file...

namespace ZgFunctions {

  //utility functions
  double map_add(vector<double> input) {
    return input[0]+input[1];
  }

  double map_div(vector<double> input) {
    return input[0]/input[1];
  }

  double map_deltaphi(vector<double> phi_values) {
    return deltaPhi(phi_values[0], phi_values[1]);
  }

  double map_deltar(vector<double> inputs) {
    return deltaR(inputs[0], inputs[1], inputs[2], inputs[3]);
  }

  double reduce_index0(vector<double> inputs) {
    return inputs[0];
  }

  double reduce_index1(vector<double> inputs) {
    return inputs[1];
  }

  double reduce_index2(vector<double> inputs) {
    return inputs[2];
  }

  double reduce_index3(vector<double> inputs) {
    return inputs[3];
  }

  double reduce_index4(vector<double> inputs) {
    return inputs[4];
  }

  double reduce_index5(vector<double> inputs) {
    return inputs[5];
  }

  //13TeV theory uncertainties from the following pages:
  //https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt13TeV
  //https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR
  
  //weight implementing variations in alphaS
  const NamedFunc sys_w_alphas("sys_w_alphas",
      [](const Baby &b) -> NamedFunc::ScalarType{
    float unc = 0.0;
    //in linear approximation N=xs*BR implies dN/N=dxs/xs+dBR/BR
    //TODO add 13.6 TeV uncertainties and year check
    if (b.type() == 200000) //ggF
      unc = 0.026+0.006;
    else if (b.type() == 200100) //VBF
      unc = 0.005+0.006;
    else if (b.type() == 200200 || b.type() == 200300) //WH
      unc = 0.009+0.006;
    else if (b.type() == 200400) //ZH
      unc = 0.009+0.006;
    else if (b.type() == 200500) //ttH
      unc = 0.020+0.006;
    else if (b.type() == 28500) //ggF H->mumu
      unc = 0.026+0.006;
    else if (b.type() == 29500) //VBF H->mumu
      unc = 0.005+0.006;
    else if (b.type() == 12500) //WH H->mumu
      unc = 0.009+0.006;
    else if (b.type() == 13500) //ZH H->mumu
      unc = 0.009+0.006;
    else if (b.type() == 9500) //ttH H->mumu
      unc = 0.020+0.006;
    return 1.0+unc;
  });

  //weight implementing variations in PDFs for ggF
  const NamedFunc sys_w_pdf_ggf("sys_w_pdf_ggf",
      [](const Baby &b) -> NamedFunc::ScalarType{
    float unc = 0.0;
    //in linear approximation N=xs*BR implies dN/N=dxs/xs+dBR/BR
    //TODO add 13.6 TeV uncertainties and year check
    if (b.type() == 200000 || b.type() == 28500) //ggF
      unc = 0.019;
    return 1.0+unc;
  });

  //weight implementing variations in PDFs for VBF/WH/ZH
  const NamedFunc sys_w_pdf_qq("sys_w_pdf_qq",
      [](const Baby &b) -> NamedFunc::ScalarType{
    float unc = 0.0;
    if (b.type() == 200100 || b.type() == 29500) //VBF
      unc = 0.021;
    else if (b.type() == 200200 || b.type() == 200300 || b.type() == 12500)
      //WH
      unc = 0.017;
    else if (b.type() == 200400 || b.type() == 13500) //ZH
      unc = 0.013;
    return 1.0+unc;
  });

  //weight implementing variations in PDFs for ttH
  const NamedFunc sys_w_pdf_tth("sys_w_pdf_tth",
      [](const Baby &b) -> NamedFunc::ScalarType{
    float unc = 0.0;
    if (b.type() == 200500 || b.type() == 9500) //ttH
      unc = 0.030;
    return 1.0+unc;
  });

  //weight implementing variations in ggF cross section
  const NamedFunc sys_w_ggf_xs("sys_w_ggf_xs",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type() == 200000 || b.type() == 28500) //ggF
      return 1.039;
    return 1.0;
  });

  //weight implementing up variation in VBF cross section
  const NamedFunc sys_w_vbf_xs_up("sys_w_vbf_xs_up",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type() == 200100 || b.type() == 29500) //VBF
      return 1.004;
    return 1.0;
  });

  //weight implementing down variation in VBF cross section
  const NamedFunc sys_w_vbf_xs_dn("sys_w_vbf_xs_dn",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type() == 200100 || b.type() == 29500) //VBF
      return 0.997;
    return 1.0;
  });

  //weight implementing up variation in WH cross section
  const NamedFunc sys_w_wh_xs_up("sys_w_wh_xs_up",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type() == 200200 || b.type() == 200300 || b.type() == 12500)
      return 1.005;
    return 1.0;
  });

  //weight implementing down variation in WH cross section
  const NamedFunc sys_w_wh_xs_dn("sys_w_wh_xs_dn",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type() == 200200 || b.type() == 200300 || b.type() == 12500)
      return 0.993;
    return 1.0;
  });

  //weight implementing up variation in ZH cross section
  const NamedFunc sys_w_zh_xs_up("sys_w_zh_xs_up",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type() == 200400 || b.type() == 13500) //ZH
      return 1.038;
    return 1.0;
  });

  //weight implementing down variation in ZH cross section
  const NamedFunc sys_w_zh_xs_dn("sys_w_zh_xs_dn",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type() == 200400 || b.type() == 13500) //ZH
      return 0.969;
    return 1.0;
  });

  //weight implementing up variation in ttH cross section
  const NamedFunc sys_w_tth_xs_up("sys_w_tth_xs_up",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type() == 200500 || b.type() == 9500) //ttH
      return 1.058;
    return 1.0;
  });

  //weight implementing down variation in ttH cross section
  const NamedFunc sys_w_tth_xs_dn("sys_w_tth_xs_dn",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type() == 200500 || b.type() == 9500) //ttH
      return 0.908;
    return 1.0;
  });

  //weight implementing variation in H->Zgamma BR
  const NamedFunc sys_w_htozg_br("sys_w_htozg_br",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type() >= 200000 && b.type() < 201000) //H->Zgamma
      return 1.0571;
    return 1.0;
  });

  //weight implementing variation in H->MuMu BR
  const NamedFunc sys_w_htomumu_br("sys_w_htomumu_br",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type() == 28500 || b.type() == 29500 || b.type() == 12500 
        || b.type() == 13500 || b.type() == 9500) //H->MuMu
      return 1.23;
    return 1.0;
  });

  //weight implementing variation in mq (mostly mtop)
  const NamedFunc sys_w_mq("sys_w_mq",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type() >= 200000 && b.type() < 201000) //H->Zgamma
      return 1.01;
    else if (b.type() == 28500 || b.type() == 29500 || b.type() == 12500 
        || b.type() == 13500 || b.type() == 9500) //H->MuMu
      return 1.01;
    return 1.0;
  });

  //weight implementing variation in run 2 lumi
  const NamedFunc sys_w_lumi_run2("sys_w_lumi_run2",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleTypeString()=="2016APV" || b.SampleTypeString()=="2016" 
        || b.SampleTypeString()=="2017" || b.SampleTypeString()=="2018") {
      return 1.016;
    }
    return 1.0;
  });

  //weight implementing variation in 2022 lumi
  const NamedFunc sys_w_lumi_2022("sys_w_lumi_2022",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleTypeString()=="2022" || b.SampleTypeString()=="2022EE") {
      return 1.013;
    }
    return 1.0;
  });

  //weight implementing variation in 2023 lumi
  const NamedFunc sys_w_lumi_2023("sys_w_lumi_2023",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleTypeString()=="2023" || b.SampleTypeString()=="2023BPix") {
      return 1.014;
    }
    return 1.0;
  });

  //weight implementing variation in PS weights
  const NamedFunc sys_w_ps_2023("sys_w_lumi_2023",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleTypeString()=="2023" || b.SampleTypeString()=="2023BPix") {
      return 1.014;
    }
    return 1.0;
  });

  //weight implementing upward variation in pileup weights 
  const NamedFunc sys_w_pu_up("sys_w_pu_up",
      [](const Baby &b) -> NamedFunc::ScalarType{
    return b.sys_pu()->at(0)/b.w_pu();
  });

  //weight implementing downward variation in pileup weights 
  const NamedFunc sys_w_pu_dn("sys_w_pu_dn",
      [](const Baby &b) -> NamedFunc::ScalarType{
    return b.sys_pu()->at(1)/b.w_pu();
  });

  //weight implementing upward variation in prefiring weights 
  const NamedFunc sys_w_prefire_up("sys_w_prefire_up",
      [](const Baby &b) -> NamedFunc::ScalarType{
    return b.sys_prefire()->at(0)/b.w_prefire();
  });

  //weight implementing downward variation in prefiring weights 
  const NamedFunc sys_w_prefire_dn("sys_w_prefire_dn",
      [](const Baby &b) -> NamedFunc::ScalarType{
    return b.sys_prefire()->at(1)/b.w_prefire();
  });

  //weight implementing upward variation in b/c-tagging weights 
  const NamedFunc sys_w_bctag_up("sys_w_bctag_up",
      [](const Baby &b) -> NamedFunc::ScalarType{
    return b.sys_bchig()->at(0)/b.w_bhig_df();
    //NOTE needs to be updated with next production
  });

  //weight implementing downward variation in b/c-tagging weights 
  const NamedFunc sys_w_bctag_dn("sys_w_bctag_dn",
      [](const Baby &b) -> NamedFunc::ScalarType{
    return b.sys_bchig()->at(1)/b.w_bhig_df();
    //NOTE needs to be updated with next production
  });

  //weight implementing upward variation in udsg-mistagging weights 
  const NamedFunc sys_w_udsgtag_up("sys_w_udsgtag_up",
      [](const Baby &b) -> NamedFunc::ScalarType{
    return b.sys_udsghig()->at(0)/b.w_bhig_df();
    //NOTE needs to be updated with next production
  });

  //weight implementing downward variation in udsg-mistagging weights 
  const NamedFunc sys_w_udsgtag_dn("sys_w_udsgtag_dn",
      [](const Baby &b) -> NamedFunc::ScalarType{
    return b.sys_udsghig()->at(1)/b.w_bhig_df();
    //NOTE needs to be updated with next production
  });

  //weight implementing upward variation in electron weights (efficiency)
  const NamedFunc sys_w_el_up("sys_w_el_up",
      [](const Baby &b) -> NamedFunc::ScalarType{
    return b.sys_el()->at(0)/b.w_el();
  });

  //weight implementing downward variation in electron weights (efficiency)
  const NamedFunc sys_w_el_dn("sys_w_el_dn",
      [](const Baby &b) -> NamedFunc::ScalarType{
    return b.sys_el()->at(1)/b.w_el();
  });

  //weight implementing upward variation in muon weights (efficiency)
  const NamedFunc sys_w_mu_up("sys_w_mu_up",
      [](const Baby &b) -> NamedFunc::ScalarType{
    return b.sys_mu()->at(0)/b.w_mu();
  });

  //weight implementing downward variation in muon weights (efficiency)
  const NamedFunc sys_w_mu_dn("sys_w_mu_dn",
      [](const Baby &b) -> NamedFunc::ScalarType{
    return b.sys_mu()->at(1)/b.w_mu();
  });

  //weight implementing upward variation in photon weights (efficiency)
  const NamedFunc sys_w_photon_up("sys_w_photon_up",
      [](const Baby &b) -> NamedFunc::ScalarType{
    return b.sys_photon()->at(0)/b.w_photon();
  });

  //weight implementing downward variation in photon weights (efficiency)
  const NamedFunc sys_w_photon_dn("sys_w_photon_dn",
      [](const Baby &b) -> NamedFunc::ScalarType{
    return b.sys_photon()->at(1)/b.w_photon();
  });

  //weight implementing upward variation in electron trigger weights
  const NamedFunc sys_w_trig_el_up_pinnacles("sys_w_trig_el_up_pinnacles",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.nel()>0)
      return b.sys_trig()->at(0)/b.w_trig();
    return 1.0;
  });

  //weight implementing downward variation in electron trigger weights
  const NamedFunc sys_w_trig_el_dn_pinnacles("sys_w_trig_el_dn_pinnacles",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.nel()>0)
      return b.sys_trig()->at(1)/b.w_trig();
    return 1.0;
  });

  //weight implementing upward variation in muon trigger weights
  const NamedFunc sys_w_trig_mu_up_pinnacles("sys_w_trig_mu_up_pinnacles",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.nmu()>0)
      return b.sys_trig()->at(0)/b.w_trig();
    return 1.0;
  });

  //weight implementing downward variation in muon trigger weights
  const NamedFunc sys_w_trig_mu_dn_pinnacles("sys_w_trig_mu_dn_pinnacles",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.nmu()>0)
      return b.sys_trig()->at(1)/b.w_trig();
    return 1.0;
  });

  //for reference, electrons and muons failing eta, dxy, or dz cuts are dropped from pico lists
  
  //el_sig and variations
  const NamedFunc sys_el_sig_default = NamedFunc("el_sig");
  const NamedFunc sys_el_sig_scaleup = NamedFunc("(sys_el_pt_scaleup>7)"
      "&&el_idLoose").Name("sys_el_sig_scaleup").EnableCaching(true);
  const NamedFunc sys_el_sig_scaledn = NamedFunc("(sys_el_pt_scaledn>7)"
      "&&el_idLoose").Name("sys_el_sig_scaledn").EnableCaching(true);
  const NamedFunc sys_el_sig_resup = NamedFunc("(sys_el_pt_resup>7)"
      "&&el_idLoose").Name("sys_el_sig_resup").EnableCaching(true);
  const NamedFunc sys_el_sig_resdn = NamedFunc("(sys_el_pt_resdn>7)"
      "&&el_idLoose").Name("sys_el_sig_resdn").EnableCaching(true);

  //el_pt and variations
  const NamedFunc sys_el_pt_default = NamedFunc("el_pt");
  const NamedFunc sys_el_pt_scaleup = NamedFunc("sys_el_pt_scaleup");
  const NamedFunc sys_el_pt_scaledn = NamedFunc("sys_el_pt_scaledn");
  const NamedFunc sys_el_pt_resup = NamedFunc("sys_el_pt_resup");
  const NamedFunc sys_el_pt_resdn = NamedFunc("sys_el_pt_resdn");

  //mu_sig and variations
  //TODO update with next production
  const NamedFunc sys_mu_sig_default = NamedFunc("mu_sig");
  const NamedFunc sys_mu_sig_scaleup = NamedFunc("mu_id&&(mu_reliso<0.35)"
      "&&(mu_sip3d<4)&&((mu_pt)>5)").Name("sys_mu_sig_scaleup")
      .EnableCaching(true);
  const NamedFunc sys_mu_sig_scaledn = NamedFunc("mu_id&&(mu_reliso<0.35)"
      "&&(mu_sip3d<4)&&((mu_pt)>5)").Name("sys_mu_sig_scaledn")
      .EnableCaching(true);
  const NamedFunc sys_mu_sig_resup = NamedFunc("mu_id&&(mu_reliso<0.35)"
      "&&(mu_sip3d<4)&&((mu_pt+mu_ptErr)>5)").Name("sys_mu_sig_resup")
      .EnableCaching(true);
  const NamedFunc sys_mu_sig_resdn = NamedFunc("mu_id&&(mu_reliso<0.35)"
      "&&(mu_sip3d<4)&&((mu_pt-mu_ptErr)>5)").Name("sys_mu_sig_resdn")
      .EnableCaching(true);

  //mu_pt and variations
  //TODO update with next production
  const NamedFunc sys_mu_pt_default = NamedFunc("mu_pt");
  const NamedFunc sys_mu_pt_scaleup = NamedFunc("mu_pt");
  const NamedFunc sys_mu_pt_scaledn = NamedFunc("mu_pt");
  const NamedFunc sys_mu_pt_resup = NamedFunc("mu_pt+mu_ptErr");
  const NamedFunc sys_mu_pt_resdn = NamedFunc("mu_pt-mu_ptErr");

  //nel and variations
  const NamedFunc sys_nel_default("nel");
  const NamedFunc sys_nel_scaleup = ReduceNamedFuncCached(sys_el_sig_scaleup,
      reduce_sum).Name("sys_nel_scaleup").EnableCaching(true);
  const NamedFunc sys_nel_scaledn = ReduceNamedFuncCached(sys_el_sig_scaledn,
      reduce_sum).Name("sys_nel_scaledn").EnableCaching(true);
  const NamedFunc sys_nel_resup = ReduceNamedFuncCached(sys_el_sig_resup,
      reduce_sum).Name("sys_nel_resup").EnableCaching(true);
  const NamedFunc sys_nel_resdn = ReduceNamedFuncCached(sys_el_sig_resdn,
      reduce_sum).Name("sys_nel_resdn").EnableCaching(true);

  //nmu and variations
  const NamedFunc sys_nmu_default("nmu");
  const NamedFunc sys_nmu_scaleup = ReduceNamedFuncCached(sys_mu_sig_scaleup,
      reduce_sum).Name("sys_nmu_scaleup").EnableCaching(true);
  const NamedFunc sys_nmu_scaledn = ReduceNamedFuncCached(sys_mu_sig_scaledn,
      reduce_sum).Name("sys_nmu_scaledn").EnableCaching(true);
  const NamedFunc sys_nmu_resup = ReduceNamedFuncCached(sys_mu_sig_resup,
      reduce_sum).Name("sys_nmu_resup").EnableCaching(true);
  const NamedFunc sys_nmu_resdn = ReduceNamedFuncCached(sys_mu_sig_resdn,
      reduce_sum).Name("sys_nmu_resdn").EnableCaching(true);

  //nlep and variations
  const NamedFunc sys_nlep_default("nlep");
  const NamedFunc sys_nlep_elscaleup = MultiMapNamedFuncCached(
      {&sys_nel_scaleup,&sys_nmu_default},map_add);
  const NamedFunc sys_nlep_elscaledn = MultiMapNamedFuncCached(
      {&sys_nel_scaledn,&sys_nmu_default},map_add);
  const NamedFunc sys_nlep_elresup = MultiMapNamedFuncCached(
      {&sys_nel_resup,&sys_nmu_default},map_add);
  const NamedFunc sys_nlep_elresdn = MultiMapNamedFuncCached(
      {&sys_nel_resdn,&sys_nmu_default},map_add);
  const NamedFunc sys_nlep_muscaleup = MultiMapNamedFuncCached(
      {&sys_nel_default,&sys_nmu_scaleup},map_add);
  const NamedFunc sys_nlep_muscaledn = MultiMapNamedFuncCached(
      {&sys_nel_default,&sys_nmu_scaledn},map_add);
  const NamedFunc sys_nlep_muresup = MultiMapNamedFuncCached(
      {&sys_nel_default,&sys_nmu_resup},map_add);
  const NamedFunc sys_nlep_muresdn = MultiMapNamedFuncCached(
      {&sys_nel_default,&sys_nmu_resdn},map_add);

  //Gets NamedFunc that is number of ll candidates with variation
  NamedFunc assign_variation_nll(const NamedFunc &variation_el_sig, 
                                 const NamedFunc &variation_mu_sig,
                                 const string &name) {
    return NamedFunc(("sys_nll_"+name).c_str(),[&variation_el_sig, 
        &variation_mu_sig] (const Baby &b) -> NamedFunc::ScalarType{
          float nll = 0;
          vector<double> el_sig = variation_el_sig.GetVector(b);
          vector<double> mu_sig = variation_mu_sig.GetVector(b);
          for (unsigned iel = 0; iel < b.el_charge()->size(); iel++) {
            if (el_sig[iel]) {
              for (unsigned iel2 = 0; iel2 < iel; iel2++) {
                if (el_sig[iel2]) {
                  if ((b.el_charge()->at(iel)+b.el_charge()->at(iel2))==0) {
                    nll += 1;
                  }
                }
              }
            }
          }
          for (unsigned imu = 0; imu < b.mu_charge()->size(); imu++) {
            if (mu_sig[imu]) {
              for (unsigned imu2 = 0; imu2 < imu; imu2++) {
                if (mu_sig[imu2]) {
                  if ((b.mu_charge()->at(imu)+b.mu_charge()->at(imu2))==0) {
                    nll += 1;
                  }
                }
              }
            }
          }
          return nll;
        }).EnableCaching(true);
  }

  //nll and variations
  const NamedFunc sys_nll_default("nll");
  const NamedFunc sys_nll_elscaleup = assign_variation_nll(sys_el_sig_scaleup, 
      sys_mu_sig_default, "elscaleup");
  const NamedFunc sys_nll_elscaledn = assign_variation_nll(sys_el_sig_scaledn, 
      sys_mu_sig_default, "elscaledn");
  const NamedFunc sys_nll_elresup = assign_variation_nll(sys_el_sig_resup, 
      sys_mu_sig_default, "elresup");
  const NamedFunc sys_nll_elresdn = assign_variation_nll(sys_el_sig_resdn, 
      sys_mu_sig_default, "elresdn");
  const NamedFunc sys_nll_muscaleup = assign_variation_nll(sys_el_sig_default, 
      sys_mu_sig_scaleup, "muscaleup");
  const NamedFunc sys_nll_muscaledn = assign_variation_nll(sys_el_sig_default, 
      sys_mu_sig_scaledn, "muscaledn");
  const NamedFunc sys_nll_muresup = assign_variation_nll(sys_el_sig_default, 
      sys_mu_sig_resup, "muresup");
  const NamedFunc sys_nll_muresdn = assign_variation_nll(sys_el_sig_default, 
      sys_mu_sig_resdn, "muresdn");

  //leading electron pt with electron scale variation up
  const NamedFunc sys_lead_el_pt_default = ReduceNamedFunc(FilterNamedFunc(
      sys_el_pt_default,sys_el_sig_default),reduce_max)
      .Name("sys_lead_el_pt_default").EnableCaching(true);
  const NamedFunc sys_lead_el_pt_scaleup = ReduceNamedFunc(FilterNamedFunc(
      sys_el_pt_scaleup,sys_el_sig_scaleup),reduce_max)
      .Name("sys_lead_el_pt_scaleup").EnableCaching(true);
  const NamedFunc sys_lead_el_pt_scaledn = ReduceNamedFunc(FilterNamedFunc(
      sys_el_pt_scaledn,sys_el_sig_scaledn),reduce_max)
      .Name("sys_lead_el_pt_scaledn").EnableCaching(true);
  const NamedFunc sys_lead_el_pt_resup = ReduceNamedFunc(FilterNamedFunc(
      sys_el_pt_resup,sys_el_sig_resup),reduce_max)
      .Name("sys_lead_el_pt_resup").EnableCaching(true);
  const NamedFunc sys_lead_el_pt_resdn = ReduceNamedFunc(FilterNamedFunc(
      sys_el_pt_resdn,sys_el_sig_resdn),reduce_max)
      .Name("sys_lead_el_pt_resdn").EnableCaching(true);

  //subleading electron pt and variations
  const NamedFunc sys_sublead_el_pt_default = ReduceNamedFunc(FilterNamedFunc(
      sys_el_pt_default,sys_el_sig_default),reduce_sublead)
      .Name("sys_sublead_el_pt_default").EnableCaching(true);
  const NamedFunc sys_sublead_el_pt_scaleup = ReduceNamedFunc(FilterNamedFunc(
      sys_el_pt_scaleup,sys_el_sig_scaleup),reduce_sublead)
      .Name("sys_sublead_el_pt_scaleup").EnableCaching(true);
  const NamedFunc sys_sublead_el_pt_scaledn = ReduceNamedFunc(FilterNamedFunc(
      sys_el_pt_scaledn,sys_el_sig_scaledn),reduce_sublead)
      .Name("sys_sublead_el_pt_scaledn").EnableCaching(true);
  const NamedFunc sys_sublead_el_pt_resup = ReduceNamedFunc(FilterNamedFunc(
      sys_el_pt_resup,sys_el_sig_resup),reduce_sublead)
      .Name("sys_sublead_el_pt_resup").EnableCaching(true);
  const NamedFunc sys_sublead_el_pt_resdn = ReduceNamedFunc(FilterNamedFunc(
      sys_el_pt_resdn,sys_el_sig_resdn),reduce_sublead)
      .Name("sys_sublead_el_pt_resdn").EnableCaching(true);

  //leading muon pt and variations
  const NamedFunc sys_lead_mu_pt_default = ReduceNamedFunc(FilterNamedFunc(
      sys_mu_pt_default,sys_mu_sig_default),
      reduce_max).Name("sys_lead_mu_pt_default").EnableCaching(true);
  const NamedFunc sys_lead_mu_pt_scaleup = ReduceNamedFunc(FilterNamedFunc(
      sys_mu_pt_scaleup,sys_mu_sig_scaleup),
      reduce_max).Name("sys_lead_mu_pt_scaleup").EnableCaching(true);
  const NamedFunc sys_lead_mu_pt_scaledn = ReduceNamedFunc(FilterNamedFunc(
      sys_mu_pt_scaledn,sys_mu_sig_scaledn),
      reduce_max).Name("sys_lead_mu_pt_scaledn").EnableCaching(true);
  const NamedFunc sys_lead_mu_pt_resup = ReduceNamedFunc(FilterNamedFunc(
      sys_mu_pt_resup,sys_mu_sig_resup),
      reduce_max).Name("sys_lead_mu_pt_resup").EnableCaching(true);
  const NamedFunc sys_lead_mu_pt_resdn = ReduceNamedFunc(FilterNamedFunc(
      sys_mu_pt_resdn,sys_mu_sig_resdn),
      reduce_max).Name("sys_lead_mu_pt_resdn").EnableCaching(true);

  //subleading muon pt with muon scale variation up
  const NamedFunc sys_sublead_mu_pt_default = ReduceNamedFunc(FilterNamedFunc(
      sys_mu_pt_default,sys_mu_sig_default),
      reduce_sublead).Name("sys_sublead_mu_pt_default").EnableCaching(true);
  const NamedFunc sys_sublead_mu_pt_scaleup = ReduceNamedFunc(FilterNamedFunc(
      sys_mu_pt_scaleup,sys_mu_sig_scaleup),
      reduce_sublead).Name("sys_sublead_mu_pt_scaleup").EnableCaching(true);
  const NamedFunc sys_sublead_mu_pt_scaledn = ReduceNamedFunc(FilterNamedFunc(
      sys_mu_pt_scaledn,sys_mu_sig_scaledn),
      reduce_sublead).Name("sys_sublead_mu_pt_scaledn").EnableCaching(true);
  const NamedFunc sys_sublead_mu_pt_resup = ReduceNamedFunc(FilterNamedFunc(
      sys_mu_pt_resup,sys_mu_sig_resup),
      reduce_sublead).Name("sys_sublead_mu_pt_resup").EnableCaching(true);
  const NamedFunc sys_sublead_mu_pt_resdn = ReduceNamedFunc(FilterNamedFunc(
      sys_mu_pt_resdn,sys_mu_sig_resdn),
      reduce_sublead).Name("sys_sublead_mu_pt_resdn").EnableCaching(true);

  //Gets NamedFunc that is trigger+pT cut flag with variation
  NamedFunc assign_variation_trig_pt(const NamedFunc &variation_el_pt, 
                                     const NamedFunc &variation_el_sig,
                                     const NamedFunc &variation_mu_pt,
                                     const NamedFunc &variation_mu_sig,
                                     const string &name) {
    return NamedFunc(("sys_trig_pt_"+name).c_str(),[&variation_el_pt, 
        &variation_el_sig, &variation_mu_pt, &variation_mu_sig] (const Baby &b) 
        -> NamedFunc::ScalarType{

          vector<double> el_sig = variation_el_sig.GetVector(b);
          vector<double> mu_sig = variation_mu_sig.GetVector(b);
          vector<double> el_pt = variation_el_pt.GetVector(b);
          vector<double> mu_pt = variation_mu_pt.GetVector(b);
          float lead_ele_pt = -999.0;
          float subl_ele_pt = -999.0;
          float lead_muo_pt = -999.0;
          float subl_muo_pt = -999.0;
          for (unsigned iel = 0; iel < el_sig.size(); iel++) {
            if (el_sig[iel]) {
              if (el_pt[iel] > lead_ele_pt) {
                subl_ele_pt = lead_ele_pt;
                lead_ele_pt = el_pt[iel];
              }
              else if (el_pt[iel] > subl_ele_pt) {
                subl_ele_pt = el_pt[iel];
              }
            }
          }
          for (unsigned imu = 0; imu < mu_sig.size(); imu++) {
            if (mu_sig[imu]) {
              if (mu_pt[imu] > lead_muo_pt) {
                subl_muo_pt = lead_muo_pt;
                lead_muo_pt = mu_pt[imu];
              }
              else if (mu_pt[imu] > subl_muo_pt) {
                subl_muo_pt = mu_pt[imu];
              }
            }
          }
          double el_cut = 35;
          double mu_cut = 25;
          if (b.SampleTypeString()=="2016APV" 
              || b.SampleTypeString()=="-2016APV" 
              || b.SampleTypeString()=="2016" 
              || b.SampleTypeString()=="-2016")
            el_cut = 30;
          else if (b.SampleTypeString()=="2017" 
                   || b.SampleTypeString()=="-2017")
            mu_cut = 28;
          if (b.nel()>=1)
            if (b.trig_single_el() && lead_ele_pt>el_cut) return true;
          if (b.nel()>=2)
            if (b.trig_double_el() && lead_ele_pt>25 
                && subl_ele_pt>15) return true;
          if (b.nmu()>=1)
            if (b.trig_single_mu() && lead_muo_pt>mu_cut) return true;
          if (b.nmu()>=2)
            if (b.trig_double_mu() && lead_muo_pt>20 
                && subl_muo_pt>10) return true;
          return false;
        }).EnableCaching(true);
  }

  //trigger and pT cuts and variations
  const NamedFunc sys_trig_pt_default = assign_variation_trig_pt(
      sys_el_pt_default, sys_el_sig_default, sys_mu_pt_default, 
      sys_mu_sig_default, "default");
  const NamedFunc sys_trig_pt_elscaleup = assign_variation_trig_pt(
      sys_el_pt_scaleup, sys_el_sig_scaleup, sys_mu_pt_default, 
      sys_mu_sig_default, "elscaleup");
  const NamedFunc sys_trig_pt_elscaledn = assign_variation_trig_pt(
      sys_el_pt_scaledn, sys_el_sig_scaledn, sys_mu_pt_default, 
      sys_mu_sig_default, "elscaledn");
  const NamedFunc sys_trig_pt_elresup = assign_variation_trig_pt(
      sys_el_pt_resup, sys_el_sig_resup, sys_mu_pt_default, 
      sys_mu_sig_default, "elresup");
  const NamedFunc sys_trig_pt_elresdn = assign_variation_trig_pt(
      sys_el_pt_resdn, sys_el_sig_resdn, sys_mu_pt_default, 
      sys_mu_sig_default, "elresdn");
  const NamedFunc sys_trig_pt_muscaleup = assign_variation_trig_pt(
      sys_el_pt_default, sys_el_sig_default, sys_mu_pt_scaleup, 
      sys_mu_sig_scaleup, "muscaleup");
  const NamedFunc sys_trig_pt_muscaledn = assign_variation_trig_pt(
      sys_el_pt_default, sys_el_sig_default, sys_mu_pt_scaledn, 
      sys_mu_sig_scaledn, "muscaledn");
  const NamedFunc sys_trig_pt_muresup = assign_variation_trig_pt(
      sys_el_pt_default, sys_el_sig_default, sys_mu_pt_resup, 
      sys_mu_sig_resup, "muresup");
  const NamedFunc sys_trig_pt_muresdn = assign_variation_trig_pt(
      sys_el_pt_default, sys_el_sig_default, sys_mu_pt_resdn, 
      sys_mu_sig_resdn, "muresdn");

  //Gets NamedFunc that is max lep miniso with variation
  NamedFunc assign_variation_max_lep_miniso(const NamedFunc &variation_el_sig,
                                            const NamedFunc &variation_mu_sig,
                                            const string &name) {
    return NamedFunc(("sys_max_lep_miniso_"+name).c_str(),[&variation_el_sig, 
        &variation_mu_sig] (const Baby &b) -> NamedFunc::ScalarType{

          vector<double> el_sig = variation_el_sig.GetVector(b);
          vector<double> mu_sig = variation_mu_sig.GetVector(b);
          float var_max_lep_miniso = -999.0;
          for (unsigned iel = 0; iel < el_sig.size(); iel++) {
            if (el_sig[iel]) {
              if (b.el_miniso()->at(iel) > var_max_lep_miniso)
                var_max_lep_miniso = b.el_miniso()->at(iel);
            }
          }
          for (unsigned imu = 0; imu < mu_sig.size(); imu++) {
            if (mu_sig[imu]) {
              if (b.mu_miniso()->at(imu) > var_max_lep_miniso)
                var_max_lep_miniso = b.mu_miniso()->at(imu); 
            }
          }
          return var_max_lep_miniso;
        }).EnableCaching(true);
  }

  //max lep miniso and variations
  const NamedFunc sys_max_lep_miniso_default = assign_variation_max_lep_miniso(
      sys_el_sig_default, sys_mu_sig_default, "default");
  const NamedFunc sys_max_lep_miniso_elscaleup = 
      assign_variation_max_lep_miniso(sys_el_sig_scaleup, sys_mu_sig_default, 
      "elscaleup");
  const NamedFunc sys_max_lep_miniso_elscaledn = 
      assign_variation_max_lep_miniso(sys_el_sig_scaledn, sys_mu_sig_default, 
      "elscaledn");
  const NamedFunc sys_max_lep_miniso_elresup = assign_variation_max_lep_miniso(
      sys_el_sig_resup, sys_mu_sig_default, "elresup");
  const NamedFunc sys_max_lep_miniso_elresdn = assign_variation_max_lep_miniso(
      sys_el_sig_resdn, sys_mu_sig_default, "elresdn");
  const NamedFunc sys_max_lep_miniso_muscaleup = 
      assign_variation_max_lep_miniso(sys_el_sig_default, sys_mu_sig_scaleup, 
      "muscaleup");
  const NamedFunc sys_max_lep_miniso_muscaledn = 
      assign_variation_max_lep_miniso(sys_el_sig_default, sys_mu_sig_scaledn, 
      "muscaledn");
  const NamedFunc sys_max_lep_miniso_muresup = assign_variation_max_lep_miniso(
      sys_el_sig_default, sys_mu_sig_resup, "muresup");
  const NamedFunc sys_max_lep_miniso_muresdn = assign_variation_max_lep_miniso(
      sys_el_sig_default, sys_mu_sig_resdn, "muresdn");

  //for reference, photons failing origin eta cuts are dropped from pico list
  
  //if photon is also a signal electron
  const NamedFunc sys_photon_sigel("photon_sigel",[](const Baby &b) 
      -> NamedFunc::VectorType{
    std::vector<double> sigel;
    for (unsigned iph = 0; iph < b.photon_elidx()->size(); iph++) {
      double el_sig = 0;
      if (b.photon_elidx()->at(iph)>=0) {
        if (b.el_sig()->at(b.photon_elidx()->at(iph))) {
          el_sig = 1;
        }
      }
      sigel.push_back(el_sig);
    }
    return sigel;
  });

  //photon_sig and variations
  const NamedFunc sys_photon_sig_default("photon_sig");
  const NamedFunc sys_photon_sig_scaleup = NamedFunc(
      ("(sys_photon_pt_scaleup>15)&&"
      "photon_elveto&&(photon_drmin>0.3)&&photon_id80")&&!sys_photon_sigel)
      .Name("sys_photon_sig_elscaleup").EnableCaching(true);
  const NamedFunc sys_photon_sig_scaledn = NamedFunc(
      ("(sys_photon_pt_scaledn>15)&&"
      "photon_elveto&&(photon_drmin>0.3)&&photon_id80")&&!sys_photon_sigel)
      .Name("sys_photon_sig_elscaledn").EnableCaching(true);
  const NamedFunc sys_photon_sig_resup = NamedFunc(
      ("(sys_photon_pt_resup>15)&&"
      "photon_elveto&&(photon_drmin>0.3)&&photon_id80")&&!sys_photon_sigel)
      .Name("sys_photon_sig_elresup").EnableCaching(true);
  const NamedFunc sys_photon_sig_resdn = NamedFunc(
      ("(sys_photon_pt_resdn>15)&&"
      "photon_elveto&&(photon_drmin>0.3)&&photon_id80")&&!sys_photon_sigel)
      .Name("sys_photon_sig_elresdn").EnableCaching(true);

  //photon_pt and variations
  const NamedFunc sys_photon_pt_default("photon_pt");
  const NamedFunc sys_photon_pt_scaleup("sys_photon_pt_scaleup");
  const NamedFunc sys_photon_pt_scaledn("sys_photon_pt_scaledn");
  const NamedFunc sys_photon_pt_resup("sys_photon_pt_resup");
  const NamedFunc sys_photon_pt_resdn("sys_photon_pt_resdn");

  //some defaults
  const NamedFunc sys_photon_eta_default("photon_eta");
  const NamedFunc sys_photon_phi_default("photon_phi");
  const NamedFunc sys_photon_drmin_default("photon_drmin");
  const NamedFunc sys_photon_drmax_default("photon_drmax");
  const NamedFunc sys_photon_idmva_default("photon_idmva");

  //nphoton and variations
  const NamedFunc sys_nphoton_default("nphoton");
  const NamedFunc sys_nphoton_scaleup = ReduceNamedFuncCached(
      sys_photon_sig_scaleup,reduce_sum).Name("sys_nphoton_scaleup");
  const NamedFunc sys_nphoton_scaledn = ReduceNamedFuncCached(
      sys_photon_sig_scaledn,reduce_sum).Name("sys_nphoton_scaledn");
  const NamedFunc sys_nphoton_resup = ReduceNamedFuncCached(
      sys_photon_sig_resup,reduce_sum).Name("sys_nphoton_resup");
  const NamedFunc sys_nphoton_resdn = ReduceNamedFuncCached(
      sys_photon_sig_resdn,reduce_sum).Name("sys_nphoton_resdn");

  //signal photon properties and variations
  const NamedFunc sys_sig_photon_pt_default = FilterNamedFuncCached(
      sys_photon_pt_default, sys_photon_sig_default).Name(
      "sys_sig_photon_pt_default");
  const NamedFunc sys_sig_photon_pt_scaleup = FilterNamedFuncCached(
      sys_photon_pt_scaleup, sys_photon_sig_scaleup).Name(
      "sys_sig_photon_pt_scaleup");
  const NamedFunc sys_sig_photon_pt_scaledn = FilterNamedFuncCached(
      sys_photon_pt_scaledn, sys_photon_sig_scaledn).Name(
      "sys_sig_photon_pt_scaledn");
  const NamedFunc sys_sig_photon_pt_resup = FilterNamedFuncCached(
      sys_photon_pt_resup, sys_photon_sig_resup).Name(
      "sys_sig_photon_pt_resup");
  const NamedFunc sys_sig_photon_pt_resdn = FilterNamedFuncCached(
      sys_photon_pt_resdn, sys_photon_sig_resdn).Name(
      "sys_sig_photon_pt_resdn");
  const NamedFunc sys_sig_photon_eta_default = FilterNamedFuncCached(
      sys_photon_eta_default, sys_photon_sig_default).Name(
      "sys_sig_photon_eta_default");
  const NamedFunc sys_sig_photon_eta_scaleup = FilterNamedFuncCached(
      sys_photon_eta_default, sys_photon_sig_scaleup).Name(
      "sys_sig_photon_eta_scaleup");
  const NamedFunc sys_sig_photon_eta_scaledn = FilterNamedFuncCached(
      sys_photon_eta_default, sys_photon_sig_scaledn).Name(
      "sys_sig_photon_eta_scaledn");
  const NamedFunc sys_sig_photon_eta_resup = FilterNamedFuncCached(
      sys_photon_eta_default, sys_photon_sig_resup).Name(
      "sys_sig_photon_eta_resup");
  const NamedFunc sys_sig_photon_eta_resdn = FilterNamedFuncCached(
      sys_photon_eta_default, sys_photon_sig_resdn).Name(
      "sys_sig_photon_eta_resdn");
  const NamedFunc sys_sig_photon_phi_default = FilterNamedFuncCached(
      sys_photon_phi_default, sys_photon_sig_default).Name(
      "sys_sig_photon_phi_default");
  const NamedFunc sys_sig_photon_phi_scaleup = FilterNamedFuncCached(
      sys_photon_phi_default, sys_photon_sig_scaleup).Name(
      "sys_sig_photon_phi_scaleup");
  const NamedFunc sys_sig_photon_phi_scaledn = FilterNamedFuncCached(
      sys_photon_phi_default, sys_photon_sig_scaledn).Name(
      "sys_sig_photon_phi_scaledn");
  const NamedFunc sys_sig_photon_phi_resup = FilterNamedFuncCached(
      sys_photon_phi_default, sys_photon_sig_resup).Name(
      "sys_sig_photon_phi_resup");
  const NamedFunc sys_sig_photon_phi_resdn = FilterNamedFuncCached(
      sys_photon_phi_default, sys_photon_sig_resdn).Name(
      "sys_sig_photon_phi_resdn");
  const NamedFunc sys_sig_photon_drmin_default = FilterNamedFuncCached(
      sys_photon_drmin_default, sys_photon_sig_default).Name(
      "sys_sig_photon_drmin_default");
  const NamedFunc sys_sig_photon_drmin_scaleup = FilterNamedFuncCached(
      sys_photon_drmin_default, sys_photon_sig_scaleup).Name(
      "sys_sig_photon_drmin_scaleup");
  const NamedFunc sys_sig_photon_drmin_scaledn = FilterNamedFuncCached(
      sys_photon_drmin_default, sys_photon_sig_scaledn).Name(
      "sys_sig_photon_drmin_scaledn");
  const NamedFunc sys_sig_photon_drmin_resup = FilterNamedFuncCached(
      sys_photon_drmin_default, sys_photon_sig_resup).Name(
      "sys_sig_photon_drmin_resup");
  const NamedFunc sys_sig_photon_drmin_resdn = FilterNamedFuncCached(
      sys_photon_drmin_default, sys_photon_sig_resdn).Name(
      "sys_sig_photon_drmin_resdn");
  const NamedFunc sys_sig_photon_drmax_default = FilterNamedFuncCached(
      sys_photon_drmax_default, sys_photon_sig_default).Name(
      "sys_sig_photon_drmax_default");
  const NamedFunc sys_sig_photon_drmax_scaleup = FilterNamedFuncCached(
      sys_photon_drmax_default, sys_photon_sig_scaleup).Name(
      "sys_sig_photon_drmax_scaleup");
  const NamedFunc sys_sig_photon_drmax_scaledn = FilterNamedFuncCached(
      sys_photon_drmax_default, sys_photon_sig_scaledn).Name(
      "sys_sig_photon_drmax_scaledn");
  const NamedFunc sys_sig_photon_drmax_resup = FilterNamedFuncCached(
      sys_photon_drmax_default, sys_photon_sig_resup).Name(
      "sys_sig_photon_drmax_resup");
  const NamedFunc sys_sig_photon_drmax_resdn = FilterNamedFuncCached(
      sys_photon_drmax_default, sys_photon_sig_resdn).Name(
      "sys_sig_photon_drmax_resdn");
  const NamedFunc sys_sig_photon_idmva_default = FilterNamedFuncCached(
      sys_photon_idmva_default, sys_photon_sig_default).Name(
      "sys_sig_photon_idmva_default");
  const NamedFunc sys_sig_photon_idmva_scaleup = FilterNamedFuncCached(
      sys_photon_idmva_default, sys_photon_sig_scaleup).Name(
      "sys_sig_photon_idmva_scaleup");
  const NamedFunc sys_sig_photon_idmva_scaledn = FilterNamedFuncCached(
      sys_photon_idmva_default, sys_photon_sig_scaledn).Name(
      "sys_sig_photon_idmva_scaledn");
  const NamedFunc sys_sig_photon_idmva_resup = FilterNamedFuncCached(
      sys_photon_idmva_default, sys_photon_sig_resup).Name(
      "sys_sig_photon_idmva_resup");
  const NamedFunc sys_sig_photon_idmva_resdn = FilterNamedFuncCached(
      sys_photon_idmva_default, sys_photon_sig_resdn).Name(
      "sys_sig_photon_idmva_resdn");

  //leading photon properties and variations
  const NamedFunc sys_lead_photon_pt_default("photon_pt[0]");
  const NamedFunc sys_lead_photon_pt_scaleup = ReduceNamedFuncCached(
      sys_sig_photon_pt_scaleup, reduce_max)
      .Name("sys_lead_photon_pt_scaleup");
  const NamedFunc sys_lead_photon_pt_scaledn = ReduceNamedFuncCached(
      sys_sig_photon_pt_scaledn, reduce_max)
      .Name("sys_lead_photon_pt_scaledn");
  const NamedFunc sys_lead_photon_pt_resup = ReduceNamedFuncCached(
      sys_sig_photon_pt_resup, reduce_max).Name("sys_lead_photon_pt_resup");
  const NamedFunc sys_lead_photon_pt_resdn = ReduceNamedFuncCached(
      sys_sig_photon_pt_resdn, reduce_max).Name("sys_lead_photon_pt_resdn");
  const NamedFunc sys_lead_photon_eta_default = NamedFunc("photon_eta[0]");
  const NamedFunc sys_lead_photon_eta_scaleup = MultiReduceNamedFuncCached(
      {&sys_sig_photon_pt_scaleup,&sys_sig_photon_eta_scaleup}, reduce_maxfirst)
      .Name("sys_lead_photon_eta_scaleup");
  const NamedFunc sys_lead_photon_eta_scaledn = MultiReduceNamedFuncCached(
      {&sys_sig_photon_pt_scaledn,&sys_sig_photon_eta_scaledn}, reduce_maxfirst)
      .Name("sys_lead_photon_eta_scaledn");
  const NamedFunc sys_lead_photon_eta_resup = MultiReduceNamedFuncCached(
      {&sys_sig_photon_pt_resup,&sys_sig_photon_eta_resup}, reduce_maxfirst)
      .Name("sys_lead_photon_eta_resup");
  const NamedFunc sys_lead_photon_eta_resdn = MultiReduceNamedFuncCached(
      {&sys_sig_photon_pt_resdn,&sys_sig_photon_eta_resdn}, reduce_maxfirst)
      .Name("sys_lead_photon_eta_resdn");
  const NamedFunc sys_lead_photon_phi_default = NamedFunc("photon_phi[0]");
  const NamedFunc sys_lead_photon_phi_scaleup = MultiReduceNamedFuncCached(
      {&sys_sig_photon_pt_scaleup,&sys_sig_photon_phi_scaleup}, reduce_maxfirst)
      .Name("sys_lead_photon_phi_scaleup");
  const NamedFunc sys_lead_photon_phi_scaledn = MultiReduceNamedFuncCached(
      {&sys_sig_photon_pt_scaledn,&sys_sig_photon_phi_scaledn}, reduce_maxfirst)
      .Name("sys_lead_photon_phi_scaledn");
  const NamedFunc sys_lead_photon_phi_resup = MultiReduceNamedFuncCached(
      {&sys_sig_photon_pt_resup,&sys_sig_photon_phi_resup}, reduce_maxfirst)
      .Name("sys_lead_photon_phi_resup");
  const NamedFunc sys_lead_photon_phi_resdn = MultiReduceNamedFuncCached(
      {&sys_sig_photon_pt_resdn,&sys_sig_photon_phi_resdn}, reduce_maxfirst)
      .Name("sys_lead_photon_phi_resdn");
  const NamedFunc sys_lead_photon_drmin_default = NamedFunc("photon_drmin[0]");
  const NamedFunc sys_lead_photon_drmin_scaleup = MultiReduceNamedFuncCached(
      {&sys_sig_photon_pt_scaleup,&sys_sig_photon_drmin_scaleup}, 
      reduce_maxfirst).Name("sys_lead_photon_drmin_scaleup");
  const NamedFunc sys_lead_photon_drmin_scaledn = MultiReduceNamedFuncCached(
      {&sys_sig_photon_pt_scaledn,&sys_sig_photon_drmin_scaledn}, 
      reduce_maxfirst).Name("sys_lead_photon_drmin_scaledn");
  const NamedFunc sys_lead_photon_drmin_resup = MultiReduceNamedFuncCached(
      {&sys_sig_photon_pt_resup,&sys_sig_photon_drmin_resup}, reduce_maxfirst)
      .Name("sys_lead_photon_drmin_resup");
  const NamedFunc sys_lead_photon_drmin_resdn = MultiReduceNamedFuncCached(
      {&sys_sig_photon_pt_resdn,&sys_sig_photon_drmin_resdn}, reduce_maxfirst)
      .Name("sys_lead_photon_drmin_resdn");
  const NamedFunc sys_lead_photon_drmax_default = NamedFunc("photon_drmax[0]");
  const NamedFunc sys_lead_photon_drmax_scaleup = MultiReduceNamedFuncCached(
      {&sys_sig_photon_pt_scaleup,&sys_sig_photon_drmax_scaleup}, 
      reduce_maxfirst).Name("sys_lead_photon_drmax_scaleup");
  const NamedFunc sys_lead_photon_drmax_scaledn = MultiReduceNamedFuncCached(
      {&sys_sig_photon_pt_scaledn,&sys_sig_photon_drmax_scaledn}, 
      reduce_maxfirst).Name("sys_lead_photon_drmax_scaledn");
  const NamedFunc sys_lead_photon_drmax_resup = MultiReduceNamedFuncCached(
      {&sys_sig_photon_pt_resup,&sys_sig_photon_drmax_resup}, reduce_maxfirst)
      .Name("sys_lead_photon_drmax_resup");
  const NamedFunc sys_lead_photon_drmax_resdn = MultiReduceNamedFuncCached(
      {&sys_sig_photon_pt_resdn,&sys_sig_photon_drmax_resdn}, reduce_maxfirst)
      .Name("sys_lead_photon_drmax_resdn");
  const NamedFunc sys_lead_photon_idmva_default = NamedFunc("photon_idmva[0]");
  const NamedFunc sys_lead_photon_idmva_scaleup = MultiReduceNamedFuncCached(
      {&sys_sig_photon_pt_scaleup,&sys_sig_photon_idmva_scaleup}, 
      reduce_maxfirst).Name("sys_lead_photon_idmva_scaleup");
  const NamedFunc sys_lead_photon_idmva_scaledn = MultiReduceNamedFuncCached(
      {&sys_sig_photon_pt_scaledn,&sys_sig_photon_idmva_scaledn}, 
      reduce_maxfirst).Name("sys_lead_photon_idmva_scaledn");
  const NamedFunc sys_lead_photon_idmva_resup = MultiReduceNamedFuncCached(
      {&sys_sig_photon_pt_resup,&sys_sig_photon_idmva_resup}, reduce_maxfirst)
      .Name("sys_lead_photon_idmva_resup");
  const NamedFunc sys_lead_photon_idmva_resdn = MultiReduceNamedFuncCached(
      {&sys_sig_photon_pt_resdn,&sys_sig_photon_idmva_resdn}, reduce_maxfirst)
      .Name("sys_lead_photon_idmva_resdn");

  //get photon rel pt error with variation
  NamedFunc assign_variation_lead_photon_relpterr(const NamedFunc &ph_sig,
      const NamedFunc &ph_pt, const string &name) {
    return NamedFunc(("sys_lead_photon_relpterr_"+name).c_str(),
        [&ph_sig, &ph_pt](const Baby &b) -> NamedFunc::ScalarType{
      vector<double> photon_pt = ph_pt.GetVector(b);
      vector<double> photon_sig = ph_sig.GetVector(b);
      double lead_photon_eta(0.0), lead_photon_pt(-999.0);
      double lead_photon_energyErr(0.0);
      for (unsigned iph = 0; iph < photon_sig.size(); iph++) {
        if (photon_sig[iph]) {
          if (lead_photon_pt < photon_pt[iph]) {
            lead_photon_pt = photon_pt[iph];
            lead_photon_eta = b.photon_eta()->at(iph);
            lead_photon_energyErr = b.photon_energyErr()->at(iph);
          }
        }
      }
      return lead_photon_energyErr/(lead_photon_pt
                                    *TMath::CosH(lead_photon_eta));
    }).EnableCaching(true);
  }

  //photon relpterr variations
  const NamedFunc sys_lead_photon_relpterr_default = 
      assign_variation_lead_photon_relpterr(sys_photon_sig_default, 
      sys_photon_pt_default, "default");
  const NamedFunc sys_lead_photon_relpterr_scaleup = 
      assign_variation_lead_photon_relpterr(sys_photon_sig_scaleup, 
      sys_photon_pt_scaleup, "scaleup");
  const NamedFunc sys_lead_photon_relpterr_scaledn = 
      assign_variation_lead_photon_relpterr(sys_photon_sig_scaledn, 
      sys_photon_pt_scaledn, "scaledn");
  const NamedFunc sys_lead_photon_relpterr_resup = 
      assign_variation_lead_photon_relpterr(sys_photon_sig_resup, 
      sys_photon_pt_resup, "resup");
  const NamedFunc sys_lead_photon_relpterr_resdn = 
      assign_variation_lead_photon_relpterr(sys_photon_sig_resdn, 
      sys_photon_pt_resdn, "resdn");

  //get Z candidate properties with variation
  //returns (pt, eta, phi, m, lepid, i1, i2, idx)
  NamedFunc assign_variation_ll(const NamedFunc &el_pt, 
                                const NamedFunc &el_sig,
                                const NamedFunc &mu_pt, 
                                const NamedFunc &mu_sig,
                                bool is_el,
                                const string &name) {
    return NamedFunc(("sys_ll_"+name).c_str(),[&el_pt, &el_sig, &mu_pt, 
        &mu_sig, is_el] (const Baby &b) -> NamedFunc::VectorType{
          int reeval_pdgid = 11;
          if (!is_el)
            reeval_pdgid = 13;
          int lepid(0), i1(0), i2(0), idx(-1); 
          //idx points to original index in ll lists, or -1 if reevaluated
          float pt(0.0), eta(0.0), phi(0.0), m(0.0);
          TLorentzVector lep1, lep2, z;
          float min_dm = 999;
          for (unsigned ill = 0; ill<b.ll_pt()->size(); ill++) {
            if (b.ll_lepid()->at(ill) != reeval_pdgid) {
              float dm = fabs(b.ll_m()->at(ill)-91.1876);
              if (dm < min_dm) {
                min_dm = dm;
                pt = b.ll_pt()->at(ill);
                eta = b.ll_eta()->at(ill);
                phi = b.ll_phi()->at(ill);
                m = b.ll_m()->at(ill);
                lepid = b.ll_lepid()->at(ill);
                i1 = b.ll_i1()->at(ill);
                i2 = b.ll_i2()->at(ill);
                idx = static_cast<int>(ill);
              }
            }
          }
          if (is_el) {
            std::vector<double> b_el_sig = el_sig.GetVector(b);
            std::vector<double> b_el_pt = el_pt.GetVector(b);
            for (unsigned iel1 = 0; iel1<b.el_pt()->size(); iel1++) {
              if (!b_el_sig.at(iel1)) continue;
              for (unsigned iel2 = iel1+1; iel2<b.el_pt()->size(); iel2++) {
                if (!b_el_sig.at(iel2)) continue;
                if (b.el_charge()->at(iel1) == b.el_charge()->at(iel2)) 
                  continue;
                lep1.SetPtEtaPhiM(b_el_pt.at(iel1), b.el_eta()->at(iel1), 
                                  b.el_phi()->at(iel1), 0.000511);
                lep2.SetPtEtaPhiM(b_el_pt.at(iel2), b.el_eta()->at(iel2), 
                                  b.el_phi()->at(iel2), 0.000511);
                z = lep1 + lep2;
                float dm = fabs(z.M()-91.1876);
                if (dm < min_dm) {
                  min_dm = dm;
                  pt = z.Pt();
                  eta = z.Eta();
                  phi = z.Phi();
                  m = z.M();
                  lepid = 11;
                  i1 = static_cast<int>(iel1);
                  i2 = static_cast<int>(iel2);
                  idx = -1;
                }
              }
            }
          } //is electron variation
          else {
            std::vector<double> b_mu_sig = mu_sig.GetVector(b);
            std::vector<double> b_mu_pt = mu_pt.GetVector(b);
            for (unsigned imu1 = 0; imu1<b.mu_pt()->size(); imu1++) {
              if (!b_mu_sig.at(imu1)) continue;
              for (unsigned imu2 = imu1+1; imu2<b.mu_pt()->size(); imu2++) {
                if (!b_mu_sig.at(imu2)) continue;
                if (b.mu_charge()->at(imu1) == b.mu_charge()->at(imu2)) 
                  continue;
                lep1.SetPtEtaPhiM(b_mu_pt.at(imu1), b.mu_eta()->at(imu1), 
                                  b.mu_phi()->at(imu1), 0.106);
                lep2.SetPtEtaPhiM(b_mu_pt.at(imu2), b.mu_eta()->at(imu2), 
                                  b.mu_phi()->at(imu2), 0.106);
                z = lep1 + lep2;
                float dm = fabs(z.M()-91.1876);
                if (dm < min_dm) {
                  min_dm = dm;
                  pt = z.Pt();
                  eta = z.Eta();
                  phi = z.Phi();
                  m = z.M();
                  lepid = 13;
                  i1 = static_cast<int>(imu1);
                  i2 = static_cast<int>(imu2);
                  idx = -1;
                }
              }
            }
          } //is muon variation
          return {static_cast<double>(pt), static_cast<double>(eta), 
                  static_cast<double>(phi), static_cast<double>(m), 
                  static_cast<double>(lepid), static_cast<double>(i1),
                  static_cast<double>(i2), static_cast<double>(idx)};
        }).EnableCaching(true);
  }

  //dilepton properties and variations
  const NamedFunc sys_ll_default = NamedFunc("sys_ll_default",[](const Baby &b) 
      -> NamedFunc::VectorType{
    return {static_cast<double>(b.ll_pt()->at(0)), 
            static_cast<double>(b.ll_eta()->at(0)), 
            static_cast<double>(b.ll_phi()->at(0)), 
            static_cast<double>(b.ll_m()->at(0)), 
            static_cast<double>(b.ll_lepid()->at(0)),
            static_cast<double>(b.ll_i1()->at(0)),
            static_cast<double>(b.ll_i2()->at(0)),
            0.0};
  }).EnableCaching(true);
  const NamedFunc sys_ll_elscaleup = assign_variation_ll(sys_el_pt_scaleup, 
      sys_el_sig_scaleup, sys_mu_pt_default, sys_mu_sig_default, true, 
      "elscaleup");
  const NamedFunc sys_ll_elscaledn = assign_variation_ll(sys_el_pt_scaledn, 
      sys_el_sig_scaledn, sys_mu_pt_default, sys_mu_sig_default, true, 
      "elscaledn");
  const NamedFunc sys_ll_elresup = assign_variation_ll(sys_el_pt_resup, 
      sys_el_sig_resup, sys_mu_pt_default, sys_mu_sig_default, true, 
      "elresup");
  const NamedFunc sys_ll_elresdn = assign_variation_ll(sys_el_pt_resdn, 
      sys_el_sig_resdn, sys_mu_pt_default, sys_mu_sig_default, true, 
      "elresdn");
  const NamedFunc sys_ll_muscaleup = assign_variation_ll(sys_el_pt_default, 
      sys_el_sig_default, sys_mu_pt_scaleup, sys_mu_sig_scaleup, false, 
      "muscaleup");
  const NamedFunc sys_ll_muscaledn = assign_variation_ll(sys_el_pt_default, 
      sys_el_sig_default, sys_mu_pt_scaledn, sys_mu_sig_scaledn, false, 
      "muscaledn");
  const NamedFunc sys_ll_muresup = assign_variation_ll(sys_el_pt_default, 
      sys_el_sig_default, sys_mu_pt_resup, sys_mu_sig_resup, false, "muresup");
  const NamedFunc sys_ll_muresdn = assign_variation_ll(sys_el_pt_default, 
      sys_el_sig_default, sys_mu_pt_resdn, sys_mu_sig_resdn, false, "muresdn");
  const NamedFunc sys_ll_m_default = ReduceNamedFuncCached(sys_ll_default, 
      reduce_index3).Name("sys_ll_m_default");
  const NamedFunc sys_ll_m_elscaleup = ReduceNamedFuncCached(sys_ll_elscaleup, 
      reduce_index3).Name("sys_ll_m_elscaleup");
  const NamedFunc sys_ll_m_elscaledn = ReduceNamedFuncCached(sys_ll_elscaledn, 
      reduce_index3).Name("sys_ll_m_elscaledn");
  const NamedFunc sys_ll_m_elresup = ReduceNamedFuncCached(sys_ll_elresup, 
      reduce_index3).Name("sys_ll_m_elresup");
  const NamedFunc sys_ll_m_elresdn = ReduceNamedFuncCached(sys_ll_elresdn, 
      reduce_index3).Name("sys_ll_m_elresdn");
  const NamedFunc sys_ll_m_muscaleup = ReduceNamedFuncCached(sys_ll_muscaleup, 
      reduce_index3).Name("sys_ll_m_muscaleup");
  const NamedFunc sys_ll_m_muscaledn = ReduceNamedFuncCached(sys_ll_muscaledn, 
      reduce_index3).Name("sys_ll_m_muscaledn");
  const NamedFunc sys_ll_m_muresup = ReduceNamedFuncCached(sys_ll_muresup, 
      reduce_index3).Name("sys_ll_m_muresup");
  const NamedFunc sys_ll_m_muresdn = ReduceNamedFuncCached(sys_ll_muresdn, 
      reduce_index3).Name("sys_ll_m_muresdn");

  //get lepton eta with variation
  NamedFunc assign_variation_lep_eta(const NamedFunc &ll, const bool lead,
      const string &name) {
    string place = "lead";
    if (!lead) place = "sublead";
    return NamedFunc(("sys_"+place+"_lep_eta_"+name).c_str(),[&ll, lead] 
        (const Baby &b) -> NamedFunc::ScalarType{
      vector<double> ll_prop = ll.GetVector(b);
      int idx = 5;
      if (!lead) idx = 6;
      if (ll_prop[4]==11) {
        return b.el_eta()->at(static_cast<int>(ll_prop[idx]));
      }
      return b.mu_eta()->at(static_cast<int>(ll_prop[idx]));
      }).EnableCaching(true);
  }

  //lepton eta with variations
  const NamedFunc sys_lead_lepton_eta_default = assign_variation_lep_eta(
      sys_ll_default, true, "default");
  const NamedFunc sys_lead_lepton_eta_elscaleup = assign_variation_lep_eta(
      sys_ll_elscaleup, true, "elscaleup");
  const NamedFunc sys_lead_lepton_eta_elscaledn = assign_variation_lep_eta(
      sys_ll_elscaledn, true, "elscaledn");
  const NamedFunc sys_lead_lepton_eta_elresup = assign_variation_lep_eta(
      sys_ll_elresup, true, "elresup");
  const NamedFunc sys_lead_lepton_eta_elresdn = assign_variation_lep_eta(
      sys_ll_elresdn, true, "elresdn");
  const NamedFunc sys_lead_lepton_eta_muscaleup = assign_variation_lep_eta(
      sys_ll_muscaleup, true, "muscaleup");
  const NamedFunc sys_lead_lepton_eta_muscaledn = assign_variation_lep_eta(
      sys_ll_muscaledn, true, "muscaledn");
  const NamedFunc sys_lead_lepton_eta_muresup = assign_variation_lep_eta(
      sys_ll_muresup, true, "muresup");
  const NamedFunc sys_lead_lepton_eta_muresdn = assign_variation_lep_eta(
      sys_ll_muresdn, true, "muresdn");
  const NamedFunc sys_sublead_lepton_eta_default = assign_variation_lep_eta(
      sys_ll_default, false, "default");
  const NamedFunc sys_sublead_lepton_eta_elscaleup = assign_variation_lep_eta(
      sys_ll_elscaleup, false, "elscaleup");
  const NamedFunc sys_sublead_lepton_eta_elscaledn = assign_variation_lep_eta(
      sys_ll_elscaledn, false, "elscaledn");
  const NamedFunc sys_sublead_lepton_eta_elresup = assign_variation_lep_eta(
      sys_ll_elresup, false, "elresup");
  const NamedFunc sys_sublead_lepton_eta_elresdn = assign_variation_lep_eta(
      sys_ll_elresdn, false, "elresdn");
  const NamedFunc sys_sublead_lepton_eta_muscaleup = assign_variation_lep_eta(
      sys_ll_muscaleup, false, "muscaleup");
  const NamedFunc sys_sublead_lepton_eta_muscaledn = assign_variation_lep_eta(
      sys_ll_muscaledn, false, "muscaledn");
  const NamedFunc sys_sublead_lepton_eta_muresup = assign_variation_lep_eta(
      sys_ll_muresup, false, "muresup");
  const NamedFunc sys_sublead_lepton_eta_muresdn = assign_variation_lep_eta(
      sys_ll_muresdn, false, "muresdn");

  //get H candidate properties with variation
  //returns (pt, eta, phi, m)
  NamedFunc assign_variation_llphoton_p4(const NamedFunc &ll_p4, 
                                         const NamedFunc &lead_photon_pt,
                                         const NamedFunc &lead_photon_eta,
                                         const NamedFunc &lead_photon_phi,
                                         const int var_pdgid,
                                         const string &name) {
    return NamedFunc(("sys_lphoton_p4_"+name).c_str(),[&ll_p4, &lead_photon_pt, 
        &lead_photon_eta, &lead_photon_phi, var_pdgid] (const Baby &b) 
        -> NamedFunc::VectorType{
          vector<double> ll_prop = ll_p4.GetVector(b);
          if (var_pdgid != 22 && ll_prop[7]==0) {
            return {static_cast<double>(b.llphoton_pt()->at(0)),
                    static_cast<double>(b.llphoton_eta()->at(0)),
                    static_cast<double>(b.llphoton_phi()->at(0)),
                    static_cast<double>(b.llphoton_m()->at(0))};
          }
          TLorentzVector photon, ll, llphoton;
          ll.SetPtEtaPhiM(ll_prop[0], ll_prop[1], ll_prop[2], ll_prop[3]);
          photon.SetPtEtaPhiM(lead_photon_pt.GetScalar(b), 
              lead_photon_eta.GetScalar(b), lead_photon_phi.GetScalar(b), 0.0);
          llphoton = ll + photon;
          return {llphoton.Pt(), llphoton.Eta(), llphoton.Phi(), llphoton.M()};
        }).EnableCaching(true);
  }

  //H four momentum and variations
  const NamedFunc sys_llphoton_p4_elscaleup = assign_variation_llphoton_p4(
      sys_ll_elscaleup, sys_lead_photon_pt_default, 
      sys_lead_photon_eta_default, sys_lead_photon_phi_default, 11, 
      "elscaleup");
  const NamedFunc sys_llphoton_p4_elscaledn = assign_variation_llphoton_p4(
      sys_ll_elscaledn, sys_lead_photon_pt_default, 
      sys_lead_photon_eta_default, sys_lead_photon_phi_default, 11, 
      "elscaledn");
  const NamedFunc sys_llphoton_p4_elresup = assign_variation_llphoton_p4(
      sys_ll_elresup, sys_lead_photon_pt_default, sys_lead_photon_eta_default, 
      sys_lead_photon_phi_default, 11, "elresup");
  const NamedFunc sys_llphoton_p4_elresdn = assign_variation_llphoton_p4(
      sys_ll_elresdn, sys_lead_photon_pt_default, sys_lead_photon_eta_default, 
      sys_lead_photon_phi_default, 11, "elresdn");
  const NamedFunc sys_llphoton_p4_muscaleup = assign_variation_llphoton_p4(
      sys_ll_muscaleup, sys_lead_photon_pt_default, 
      sys_lead_photon_eta_default, sys_lead_photon_phi_default, 13, 
      "muscaleup");
  const NamedFunc sys_llphoton_p4_muscaledn = assign_variation_llphoton_p4(
      sys_ll_muscaledn, sys_lead_photon_pt_default, 
      sys_lead_photon_eta_default, sys_lead_photon_phi_default, 13, 
      "muscaledn");
  const NamedFunc sys_llphoton_p4_muresup = assign_variation_llphoton_p4(
      sys_ll_muresup, sys_lead_photon_pt_default, sys_lead_photon_eta_default, 
      sys_lead_photon_phi_default, 13, "muresup");
  const NamedFunc sys_llphoton_p4_muresdn = assign_variation_llphoton_p4(
      sys_ll_muresdn, sys_lead_photon_pt_default, sys_lead_photon_eta_default, 
      sys_lead_photon_phi_default, 13, "muresdn");
  const NamedFunc sys_llphoton_p4_phscaleup = assign_variation_llphoton_p4(
      sys_ll_default, sys_lead_photon_pt_scaleup, sys_lead_photon_eta_scaleup, 
      sys_lead_photon_phi_scaleup, 22, "phscaleup");
  const NamedFunc sys_llphoton_p4_phscaledn = assign_variation_llphoton_p4(
      sys_ll_default, sys_lead_photon_pt_scaledn, sys_lead_photon_eta_scaledn, 
      sys_lead_photon_phi_scaledn, 22, "phscaledn");
  const NamedFunc sys_llphoton_p4_phresup = assign_variation_llphoton_p4(
      sys_ll_default, sys_lead_photon_pt_resup, sys_lead_photon_eta_resup, 
      sys_lead_photon_phi_resup, 22, "phresup");
  const NamedFunc sys_llphoton_p4_phresdn = assign_variation_llphoton_p4(
      sys_ll_default, sys_lead_photon_pt_resdn, sys_lead_photon_eta_resdn, 
      sys_lead_photon_phi_resdn, 22, "phresdn");

  //Higgs candidate mass and variations
  const NamedFunc sys_llphoton_m_default = NamedFunc("llphoton_m[0]");
  const NamedFunc sys_llphoton_m_elscaleup = ReduceNamedFuncCached(
      sys_llphoton_p4_elscaleup, reduce_index3)
      .Name("sys_llphoton_m_elscaleup");
  const NamedFunc sys_llphoton_m_elscaledn = ReduceNamedFuncCached(
      sys_llphoton_p4_elscaledn, reduce_index3)
      .Name("sys_llphoton_m_elscaledn");
  const NamedFunc sys_llphoton_m_elresup = ReduceNamedFuncCached(
      sys_llphoton_p4_elresup, reduce_index3)
      .Name("sys_llphoton_m_elresup");
  const NamedFunc sys_llphoton_m_elresdn = ReduceNamedFuncCached(
      sys_llphoton_p4_elresdn, reduce_index3)
      .Name("sys_llphoton_m_elresdn");
  const NamedFunc sys_llphoton_m_muscaleup = ReduceNamedFuncCached(
      sys_llphoton_p4_muscaleup, reduce_index3)
      .Name("sys_llphoton_m_muscaleup");
  const NamedFunc sys_llphoton_m_muscaledn = ReduceNamedFuncCached(
      sys_llphoton_p4_muscaledn, reduce_index3)
      .Name("sys_llphoton_m_muscaledn");
  const NamedFunc sys_llphoton_m_muresup = ReduceNamedFuncCached(
      sys_llphoton_p4_muresup, reduce_index3)
      .Name("sys_llphoton_m_muresup");
  const NamedFunc sys_llphoton_m_muresdn = ReduceNamedFuncCached(
      sys_llphoton_p4_muresdn, reduce_index3)
      .Name("sys_llphoton_m_muresdn");
  const NamedFunc sys_llphoton_m_phscaleup = ReduceNamedFuncCached(
      sys_llphoton_p4_phscaleup, reduce_index3)
      .Name("sys_llphoton_m_phscaleup");
  const NamedFunc sys_llphoton_m_phscaledn = ReduceNamedFuncCached(
      sys_llphoton_p4_phscaledn, reduce_index3)
      .Name("sys_llphoton_m_phscaledn");
  const NamedFunc sys_llphoton_m_phresup = ReduceNamedFuncCached(
      sys_llphoton_p4_phresup, reduce_index3)
      .Name("sys_llphoton_m_phresup");
  const NamedFunc sys_llphoton_m_phresdn = ReduceNamedFuncCached(
      sys_llphoton_p4_phresdn, reduce_index3)
      .Name("sys_llphoton_m_phresdn");

  //Higgs candidate mass and variations
  const NamedFunc sys_llphoton_pt_default = NamedFunc("llphoton_pt[0]");
  const NamedFunc sys_llphoton_pt_elscaleup = ReduceNamedFuncCached(
      sys_llphoton_p4_elscaleup, reduce_index0)
      .Name("sys_llphoton_pt_elscaleup");
  const NamedFunc sys_llphoton_pt_elscaledn = ReduceNamedFuncCached(
      sys_llphoton_p4_elscaledn, reduce_index0)
      .Name("sys_llphoton_pt_elscaledn");
  const NamedFunc sys_llphoton_pt_elresup = ReduceNamedFuncCached(
      sys_llphoton_p4_elresup, reduce_index0)
      .Name("sys_llphoton_pt_elresup");
  const NamedFunc sys_llphoton_pt_elresdn = ReduceNamedFuncCached(
      sys_llphoton_p4_elresdn, reduce_index0)
      .Name("sys_llphoton_pt_elresdn");
  const NamedFunc sys_llphoton_pt_muscaleup = ReduceNamedFuncCached(
      sys_llphoton_p4_muscaleup, reduce_index0)
      .Name("sys_llphoton_pt_muscaleup");
  const NamedFunc sys_llphoton_pt_muscaledn = ReduceNamedFuncCached(
      sys_llphoton_p4_muscaledn, reduce_index0)
      .Name("sys_llphoton_pt_muscaledn");
  const NamedFunc sys_llphoton_pt_muresup = ReduceNamedFuncCached(
      sys_llphoton_p4_muresup, reduce_index0)
      .Name("sys_llphoton_pt_muresup");
  const NamedFunc sys_llphoton_pt_muresdn = ReduceNamedFuncCached(
      sys_llphoton_p4_muresdn, reduce_index0)
      .Name("sys_llphoton_pt_muresdn");
  const NamedFunc sys_llphoton_pt_phscaleup = ReduceNamedFuncCached(
      sys_llphoton_p4_phscaleup, reduce_index0)
      .Name("sys_llphoton_pt_phscaleup");
  const NamedFunc sys_llphoton_pt_phscaledn = ReduceNamedFuncCached(
      sys_llphoton_p4_phscaledn, reduce_index0)
      .Name("sys_llphoton_pt_phscaledn");
  const NamedFunc sys_llphoton_pt_phresup = ReduceNamedFuncCached(
      sys_llphoton_p4_phresup, reduce_index0)
      .Name("sys_llphoton_pt_phresup");
  const NamedFunc sys_llphoton_pt_phresdn = ReduceNamedFuncCached(
      sys_llphoton_p4_phresdn, reduce_index0)
      .Name("sys_llphoton_pt_phresdn");

  //Gets NamedFunc that is (pt1, eta1, phi1, m1, pt2, eta2, phi2, m2) of lepton
  //refit pT with variation
  NamedFunc assign_variation_lep_refit(
      const NamedFunc &el_pt, 
      const NamedFunc &mu_pt, 
      const NamedFunc &ll, 
      const string &name) {
    //require nll>=1 
    return NamedFunc(("sys_lep_refit_"+name).c_str(),[&el_pt, &mu_pt, &ll] 
        (const Baby &b) -> NamedFunc::VectorType{
      vector<double> ll_ = ll.GetVector(b);
      int lepid = static_cast<int>(ll_[4]);
      int ll_i1 = static_cast<int>(ll_[5]);
      int ll_i2 = static_cast<int>(ll_[6]);
      vector<double> lep_refit = get_lep_custom_refit(b, el_pt, mu_pt, 
      lepid, ll_i1, ll_i2);
      return lep_refit;
    }).EnableCaching(true);
  }

  //lepton 4 momentum and variations
  const NamedFunc sys_lep_refit_default = assign_variation_lep_refit(
      sys_el_pt_default, sys_mu_pt_default, sys_ll_default, "default");
  const NamedFunc sys_lep_refit_elscaleup = assign_variation_lep_refit(
      sys_el_pt_scaleup, sys_mu_pt_default, sys_ll_elscaleup, "elscaleup");
  const NamedFunc sys_lep_refit_elscaledn = assign_variation_lep_refit(
      sys_el_pt_scaledn, sys_mu_pt_default, sys_ll_elscaledn, "elscaledn");
  const NamedFunc sys_lep_refit_elresup = assign_variation_lep_refit(
      sys_el_pt_resup, sys_mu_pt_default, sys_ll_elresup, "elresup");
  const NamedFunc sys_lep_refit_elresdn = assign_variation_lep_refit(
      sys_el_pt_resdn, sys_mu_pt_default, sys_ll_elresdn, "elresdn");
  const NamedFunc sys_lep_refit_muscaleup = assign_variation_lep_refit(
      sys_el_pt_default, sys_mu_pt_scaleup, sys_ll_muscaleup, "muscaleup");
  const NamedFunc sys_lep_refit_muscaledn = assign_variation_lep_refit(
      sys_el_pt_default, sys_mu_pt_scaledn, sys_ll_muscaledn, "muscaledn");
  const NamedFunc sys_lep_refit_muresup = assign_variation_lep_refit(
      sys_el_pt_default, sys_mu_pt_resup, sys_ll_muresup, "muresup");
  const NamedFunc sys_lep_refit_muresdn = assign_variation_lep_refit(
      sys_el_pt_default, sys_mu_pt_resdn, sys_ll_muresdn, "muresdn");

  //Gets NamedFunc that is assigns Z candidate four momentum with variation
  NamedFunc assign_variation_ll_refit_p4(const NamedFunc &ll, 
      const NamedFunc &lep_refit, const string &name) {
    return NamedFunc(("sys_ll_refit_p4_"+name).c_str(),[&ll, &lep_refit]
        (const Baby &b) -> NamedFunc::VectorType{
          int ll_idx = ll.GetVector(b)[7];
          if (ll_idx==0)
            return {b.ll_refit_pt(), b.ll_refit_eta(), 
                    b.ll_refit_phi(), b.ll_refit_m()};
          vector<double> lep_refit0 = lep_refit.GetVector(b);
          TLorentzVector l1, l2, dilep;
          l1.SetPtEtaPhiM(lep_refit0[0], lep_refit0[1], lep_refit0[2], 
                          lep_refit0[3]);
          l2.SetPtEtaPhiM(lep_refit0[4], lep_refit0[5], lep_refit0[6], 
                          lep_refit0[7]);
          dilep = l1+l2;
          return {dilep.Pt(), dilep.Eta(), dilep.Phi(), dilep.M()};
        }).EnableCaching(true);
  }

  //Z candidate refit 4 momentum and variations
  const NamedFunc sys_ll_refit_p4_default = assign_variation_ll_refit_p4(
      sys_ll_default, sys_lep_refit_elscaleup, "default");
  const NamedFunc sys_ll_refit_p4_elscaleup = assign_variation_ll_refit_p4(
      sys_ll_elscaleup, sys_lep_refit_elscaleup, "elscaleup");
  const NamedFunc sys_ll_refit_p4_elscaledn = assign_variation_ll_refit_p4(
      sys_ll_elscaledn, sys_lep_refit_elscaledn, "elscaledn");
  const NamedFunc sys_ll_refit_p4_elresup = assign_variation_ll_refit_p4(
      sys_ll_elresup, sys_lep_refit_elresup, "elresup");
  const NamedFunc sys_ll_refit_p4_elresdn = assign_variation_ll_refit_p4(
      sys_ll_elresdn, sys_lep_refit_elresdn, "elresdn");
  const NamedFunc sys_ll_refit_p4_muscaleup = assign_variation_ll_refit_p4(
      sys_ll_muscaleup, sys_lep_refit_muscaleup, "muscaleup");
  const NamedFunc sys_ll_refit_p4_muscaledn = assign_variation_ll_refit_p4(
      sys_ll_muscaledn, sys_lep_refit_muscaledn, "muscaledn");
  const NamedFunc sys_ll_refit_p4_muresup = assign_variation_ll_refit_p4(
      sys_ll_muresup, sys_lep_refit_muresup, "muresup");
  const NamedFunc sys_ll_refit_p4_muresdn = assign_variation_ll_refit_p4(
      sys_ll_muresdn, sys_lep_refit_muresdn, "muresdn");

  //Gets NamedFunc that assigns Higgs four momentum with variation
  NamedFunc assign_variation_llphoton_refit_p4(
      const NamedFunc &ll, const NamedFunc &ll_refit_p4, 
      const NamedFunc &lead_photon_pt, const NamedFunc &lead_photon_eta, 
      const NamedFunc &lead_photon_phi, bool is_phvar, const string &name) {
    return NamedFunc(("sys_llphoton_refit_p4_"+name).c_str(),[&ll, 
        &ll_refit_p4, &lead_photon_pt, &lead_photon_eta, &lead_photon_phi, 
        is_phvar] (const Baby &b) -> NamedFunc::VectorType{
          if (!is_phvar && ll.GetVector(b)[7]==0) {
            return {b.llphoton_refit_pt(), b.llphoton_refit_eta(), 
                    b.llphoton_refit_phi(), b.llphoton_refit_m()};
          }
          TLorentzVector dilep, photon, llphoton;
          vector<double> ll_refit_p = ll_refit_p4.GetVector(b);
          dilep.SetPtEtaPhiM(ll_refit_p[0], ll_refit_p[1], ll_refit_p[2], 
                          ll_refit_p[3]);
          photon.SetPtEtaPhiM(lead_photon_pt.GetScalar(b), 
                              lead_photon_eta.GetScalar(b), 
                              lead_photon_phi.GetScalar(b), 0.0);
          llphoton = dilep+photon;
          return {llphoton.Pt(), llphoton.Eta(), llphoton.Phi(), 
                  llphoton.M()};
        }).EnableCaching(true);
  }

  //H candidate refit 4 momentum and variations
  const NamedFunc sys_llphoton_refit_p4_elscaleup = 
      assign_variation_llphoton_refit_p4(sys_ll_elscaleup, 
      sys_ll_refit_p4_elscaleup, sys_lead_photon_pt_default, 
      sys_lead_photon_eta_default, sys_lead_photon_phi_default, false, 
      "elscaleup");
  const NamedFunc sys_llphoton_refit_p4_elscaledn = 
      assign_variation_llphoton_refit_p4(sys_ll_elscaledn, 
      sys_ll_refit_p4_elscaledn, sys_lead_photon_pt_default, 
      sys_lead_photon_eta_default, sys_lead_photon_phi_default, false, 
      "elscaledn");
  const NamedFunc sys_llphoton_refit_p4_elresup = 
      assign_variation_llphoton_refit_p4(sys_ll_elresup, 
      sys_ll_refit_p4_elresup, sys_lead_photon_pt_default, 
      sys_lead_photon_eta_default, sys_lead_photon_phi_default, false, 
      "elresup");
  const NamedFunc sys_llphoton_refit_p4_elresdn = 
      assign_variation_llphoton_refit_p4(sys_ll_elresdn, 
      sys_ll_refit_p4_elresdn, sys_lead_photon_pt_default, 
      sys_lead_photon_eta_default, sys_lead_photon_phi_default, false, 
      "elresdn");
  const NamedFunc sys_llphoton_refit_p4_muscaleup = 
      assign_variation_llphoton_refit_p4(sys_ll_muscaleup, 
      sys_ll_refit_p4_muscaleup, sys_lead_photon_pt_default, 
      sys_lead_photon_eta_default, sys_lead_photon_phi_default, false, 
      "muscaleup");
  const NamedFunc sys_llphoton_refit_p4_muscaledn = 
      assign_variation_llphoton_refit_p4(sys_ll_muscaledn, 
      sys_ll_refit_p4_muscaledn, sys_lead_photon_pt_default, 
      sys_lead_photon_eta_default, sys_lead_photon_phi_default, false, 
      "muscaledn");
  const NamedFunc sys_llphoton_refit_p4_muresup = 
      assign_variation_llphoton_refit_p4(sys_ll_muresup, 
      sys_ll_refit_p4_muresup, sys_lead_photon_pt_default, 
      sys_lead_photon_eta_default, sys_lead_photon_phi_default, false, 
      "muresup");
  const NamedFunc sys_llphoton_refit_p4_muresdn = 
      assign_variation_llphoton_refit_p4(sys_ll_muresdn, 
      sys_ll_refit_p4_muresdn, sys_lead_photon_pt_default, 
      sys_lead_photon_eta_default, sys_lead_photon_phi_default, false, 
      "muresdn");
  const NamedFunc sys_llphoton_refit_p4_phscaleup = 
      assign_variation_llphoton_refit_p4(sys_ll_default, 
      sys_ll_refit_p4_default, sys_lead_photon_pt_scaleup, 
      sys_lead_photon_eta_scaleup, sys_lead_photon_phi_scaleup, true, 
      "phscaleup");
  const NamedFunc sys_llphoton_refit_p4_phscaledn = 
      assign_variation_llphoton_refit_p4(sys_ll_default, 
      sys_ll_refit_p4_default, sys_lead_photon_pt_scaledn, 
      sys_lead_photon_eta_scaledn, sys_lead_photon_phi_scaledn, true, 
      "phscaledn");
  const NamedFunc sys_llphoton_refit_p4_phresup = 
      assign_variation_llphoton_refit_p4(sys_ll_default, 
      sys_ll_refit_p4_default, sys_lead_photon_pt_resup, 
      sys_lead_photon_eta_resup, sys_lead_photon_phi_resup, true, "phresup");
  const NamedFunc sys_llphoton_refit_p4_phresdn = 
      assign_variation_llphoton_refit_p4(sys_ll_default, 
      sys_ll_refit_p4_default, sys_lead_photon_pt_resdn, 
      sys_lead_photon_eta_resdn, sys_lead_photon_phi_resdn, true, "phresdn");

  //H candidate refit pt 
  const NamedFunc sys_llphoton_refit_pt_default = 
      NamedFunc("llphoton_refit_pt");
  const NamedFunc sys_llphoton_refit_pt_elscaleup = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_elscaleup, reduce_index0).Name(
      "sys_llphoton_refit_pt_elscaleup");
  const NamedFunc sys_llphoton_refit_pt_elscaledn = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_elscaledn, reduce_index0).Name(
      "sys_llphoton_refit_pt_elscaledn");
  const NamedFunc sys_llphoton_refit_pt_elresup = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_elresup, reduce_index0).Name(
      "sys_llphoton_refit_pt_elresup");
  const NamedFunc sys_llphoton_refit_pt_elresdn = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_elresdn, reduce_index0).Name(
      "sys_llphoton_refit_pt_elresdn");
  const NamedFunc sys_llphoton_refit_pt_muscaleup = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_muscaleup, reduce_index0).Name(
      "sys_llphoton_refit_pt_muscaleup");
  const NamedFunc sys_llphoton_refit_pt_muscaledn = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_muscaledn, reduce_index0).Name(
      "sys_llphoton_refit_pt_muscaledn");
  const NamedFunc sys_llphoton_refit_pt_muresup = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_muresup, reduce_index0).Name(
      "sys_llphoton_refit_pt_muresup");
  const NamedFunc sys_llphoton_refit_pt_muresdn = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_muresdn, reduce_index0).Name(
      "sys_llphoton_refit_pt_muresdn");
  const NamedFunc sys_llphoton_refit_pt_phscaleup = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_phscaleup, reduce_index0).Name(
      "sys_llphoton_refit_pt_phscaleup");
  const NamedFunc sys_llphoton_refit_pt_phscaledn = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_phscaledn, reduce_index0).Name(
      "sys_llphoton_refit_pt_phscaledn");
  const NamedFunc sys_llphoton_refit_pt_phresup = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_phresup, reduce_index0).Name(
      "sys_llphoton_refit_pt_phresup");
  const NamedFunc sys_llphoton_refit_pt_phresdn = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_phresdn, reduce_index0).Name(
      "sys_llphoton_refit_pt_phresdn");

  //H candidate refit m
  const NamedFunc sys_llphoton_refit_m_default = 
      NamedFunc("llphoton_refit_m");
  const NamedFunc sys_llphoton_refit_m_elscaleup = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_elscaleup, reduce_index3).Name(
      "sys_llphoton_refit_m_elscaleup");
  const NamedFunc sys_llphoton_refit_m_elscaledn = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_elscaledn, reduce_index3).Name(
      "sys_llphoton_refit_m_elscaledn");
  const NamedFunc sys_llphoton_refit_m_elresup = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_elresup, reduce_index3).Name(
      "sys_llphoton_refit_m_elresup");
  const NamedFunc sys_llphoton_refit_m_elresdn = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_elresdn, reduce_index3).Name(
      "sys_llphoton_refit_m_elresdn");
  const NamedFunc sys_llphoton_refit_m_muscaleup = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_muscaleup, reduce_index3).Name(
      "sys_llphoton_refit_m_muscaleup");
  const NamedFunc sys_llphoton_refit_m_muscaledn = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_muscaledn, reduce_index3).Name(
      "sys_llphoton_refit_m_muscaledn");
  const NamedFunc sys_llphoton_refit_m_muresup = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_muresup, reduce_index3).Name(
      "sys_llphoton_refit_m_muresup");
  const NamedFunc sys_llphoton_refit_m_muresdn = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_muresdn, reduce_index3).Name(
      "sys_llphoton_refit_m_muresdn");
  const NamedFunc sys_llphoton_refit_m_phscaleup = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_phscaleup, reduce_index3).Name(
      "sys_llphoton_refit_m_phscaleup");
  const NamedFunc sys_llphoton_refit_m_phscaledn = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_phscaledn, reduce_index3).Name(
      "sys_llphoton_refit_m_phscaledn");
  const NamedFunc sys_llphoton_refit_m_phresup = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_phresup, reduce_index3).Name(
      "sys_llphoton_refit_m_phresup");
  const NamedFunc sys_llphoton_refit_m_phresdn = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_phresdn, reduce_index3).Name(
      "sys_llphoton_refit_m_phresdn");

  //H candidate refit pt/m and variations
  const NamedFunc sys_llphoton_refit_relpt_default = MultiMapNamedFuncCached(
      {&sys_llphoton_refit_pt_default,&sys_llphoton_refit_m_default}, map_div);
  const NamedFunc sys_llphoton_refit_relpt_elscaleup = MultiMapNamedFuncCached(
      {&sys_llphoton_refit_pt_elscaleup,&sys_llphoton_refit_m_elscaleup}, 
      map_div);
  const NamedFunc sys_llphoton_refit_relpt_elscaledn = MultiMapNamedFuncCached(
      {&sys_llphoton_refit_pt_elscaledn,&sys_llphoton_refit_m_elscaledn}, 
      map_div);
  const NamedFunc sys_llphoton_refit_relpt_elresup = MultiMapNamedFuncCached(
      {&sys_llphoton_refit_pt_elresup,&sys_llphoton_refit_m_elresup}, 
      map_div);
  const NamedFunc sys_llphoton_refit_relpt_elresdn = MultiMapNamedFuncCached(
      {&sys_llphoton_refit_pt_elresdn,&sys_llphoton_refit_m_elresdn}, 
      map_div);
  const NamedFunc sys_llphoton_refit_relpt_muscaleup = MultiMapNamedFuncCached(
      {&sys_llphoton_refit_pt_muscaleup,&sys_llphoton_refit_m_muscaleup}, 
      map_div);
  const NamedFunc sys_llphoton_refit_relpt_muscaledn = MultiMapNamedFuncCached(
      {&sys_llphoton_refit_pt_muscaledn,&sys_llphoton_refit_m_muscaledn}, 
      map_div);
  const NamedFunc sys_llphoton_refit_relpt_muresup = MultiMapNamedFuncCached(
      {&sys_llphoton_refit_pt_muresup,&sys_llphoton_refit_m_muresup}, 
      map_div);
  const NamedFunc sys_llphoton_refit_relpt_muresdn = MultiMapNamedFuncCached(
      {&sys_llphoton_refit_pt_muresdn,&sys_llphoton_refit_m_muresdn}, 
      map_div);
  const NamedFunc sys_llphoton_refit_relpt_phscaleup = MultiMapNamedFuncCached(
      {&sys_llphoton_refit_pt_phscaleup,&sys_llphoton_refit_m_phscaleup}, 
      map_div);
  const NamedFunc sys_llphoton_refit_relpt_phscaledn = MultiMapNamedFuncCached(
      {&sys_llphoton_refit_pt_phscaledn,&sys_llphoton_refit_m_phscaledn}, 
      map_div);
  const NamedFunc sys_llphoton_refit_relpt_phresup = MultiMapNamedFuncCached(
      {&sys_llphoton_refit_pt_phresup,&sys_llphoton_refit_m_phresup}, 
      map_div);
  const NamedFunc sys_llphoton_refit_relpt_phresdn = MultiMapNamedFuncCached(
      {&sys_llphoton_refit_pt_phresdn,&sys_llphoton_refit_m_phresdn}, 
      map_div);

  //H candidate refit phi and variations
  const NamedFunc sys_llphoton_refit_phi_default = 
      NamedFunc("llphoton_refit_phi");
  const NamedFunc sys_llphoton_refit_phi_elscaleup = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_elscaleup, reduce_index2).Name(
      "sys_llphoton_refit_phi_elscaleup");
  const NamedFunc sys_llphoton_refit_phi_elscaledn = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_elscaledn, reduce_index2).Name(
      "sys_llphoton_refit_phi_elscaledn");
  const NamedFunc sys_llphoton_refit_phi_elresup = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_elresup, reduce_index2).Name(
      "sys_llphoton_refit_phi_elresup");
  const NamedFunc sys_llphoton_refit_phi_elresdn = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_elresdn, reduce_index2).Name(
      "sys_llphoton_refit_phi_elresdn");
  const NamedFunc sys_llphoton_refit_phi_muscaleup = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_muscaleup, reduce_index2).Name(
      "sys_llphoton_refit_phi_muscaleup");
  const NamedFunc sys_llphoton_refit_phi_muscaledn = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_muscaledn, reduce_index2).Name(
      "sys_llphoton_refit_phi_muscaledn");
  const NamedFunc sys_llphoton_refit_phi_muresup = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_muresup, reduce_index2).Name(
      "sys_llphoton_refit_phi_muresup");
  const NamedFunc sys_llphoton_refit_phi_muresdn = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_muresdn, reduce_index2).Name(
      "sys_llphoton_refit_phi_muresdn");
  const NamedFunc sys_llphoton_refit_phi_phscaleup = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_phscaleup, reduce_index2).Name(
      "sys_llphoton_refit_phi_phscaleup");
  const NamedFunc sys_llphoton_refit_phi_phscaledn = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_phscaledn, reduce_index2).Name(
      "sys_llphoton_refit_phi_phscaledn");
  const NamedFunc sys_llphoton_refit_phi_phresup = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_phresup, reduce_index2).Name(
      "sys_llphoton_refit_phi_phresup");
  const NamedFunc sys_llphoton_refit_phi_phresdn = ReduceNamedFuncCached(
      sys_llphoton_refit_p4_phresdn, reduce_index2).Name(
      "sys_llphoton_refit_phi_phresdn");

  //Gets NamedFunc that assigns kinematic angles with variation
  NamedFunc assign_variation_llphoton_refit_angles(const NamedFunc &var_ll,
      const NamedFunc &lep_refit, const NamedFunc &var_photon_pt, 
      const NamedFunc &var_photon_eta, const NamedFunc &var_photon_phi, 
      const string name) {
    return NamedFunc(("sys_llphoton_refit_angles_"+name).c_str(),[&var_ll, 
        &lep_refit, &var_photon_pt, &var_photon_eta, &var_photon_phi] 
        (const Baby &b) -> NamedFunc::VectorType{
          TLorentzVector lplus, lminus, ph;
          vector<double> ll_properties = var_ll.GetVector(b);
          bool default_leptons = false;
          bool has_fsr_photon = false;
          //note: only refit pt saved, need to redo refit if FSR photon
          if (ll_properties[7]==0) {
            default_leptons = true;
            for (unsigned ifsr = 0; ifsr < b.fsrphoton_muonidx()->size(); 
                 ifsr++) {
              if (b.ll_lepid()->at(0)==13) {
                if ((b.fsrphoton_muonidx()->at(ifsr)==b.ll_i1()->at(0)) ||
                    (b.fsrphoton_muonidx()->at(ifsr)==b.ll_i2()->at(0)))
                  has_fsr_photon = true;
              }
            }
          }
          if (default_leptons && !has_fsr_photon) {
            if (b.ll_lepid()->at(0)==11) {
              if (b.el_charge()->at(b.ll_i1()->at(0)) == 1) {
                lplus.SetPtEtaPhiM(b.ll_refit_l1_pt(),
                                   b.el_eta()->at(b.ll_i1()->at(0)),
                                   b.el_phi()->at(b.ll_i1()->at(0)),
                                   0.00511);
                lminus.SetPtEtaPhiM(b.ll_refit_l2_pt(),
                                    b.el_eta()->at(b.ll_i2()->at(0)),
                                    b.el_phi()->at(b.ll_i2()->at(0)),
                                    0.00511);
              }
              else {
                lminus.SetPtEtaPhiM(b.ll_refit_l1_pt(),
                                    b.el_eta()->at(b.ll_i1()->at(0)),
                                    b.el_phi()->at(b.ll_i1()->at(0)),
                                    0.00511);
                lplus.SetPtEtaPhiM(b.ll_refit_l2_pt(),
                                   b.el_eta()->at(b.ll_i2()->at(0)),
                                   b.el_phi()->at(b.ll_i2()->at(0)),
                                   0.00511);
              }
            }
            else {
              if (b.mu_charge()->at(b.ll_i1()->at(0)) == 1) {
                lplus.SetPtEtaPhiM(b.ll_refit_l1_pt(),
                                   b.mu_eta()->at(b.ll_i1()->at(0)),
                                   b.mu_phi()->at(b.ll_i1()->at(0)),
                                   0.106);
                lminus.SetPtEtaPhiM(b.ll_refit_l2_pt(),
                                    b.mu_eta()->at(b.ll_i2()->at(0)),
                                    b.mu_phi()->at(b.ll_i2()->at(0)),
                                    0.106);
              }
              else {
                lminus.SetPtEtaPhiM(b.ll_refit_l1_pt(),
                                    b.mu_eta()->at(b.ll_i1()->at(0)),
                                    b.mu_phi()->at(b.ll_i1()->at(0)),
                                    0.106);
                lplus.SetPtEtaPhiM(b.ll_refit_l2_pt(),
                                   b.mu_eta()->at(b.ll_i2()->at(0)),
                                   b.mu_phi()->at(b.ll_i2()->at(0)),
                                   0.106);
              }
            }
          } //no FSR photon
          else {
            vector<double> lep_refit0 = lep_refit.GetVector(b);
            if ((ll_properties[4]==11 
                && b.el_charge()->at(static_cast<int>(ll_properties[5]))==1)
                ||(ll_properties[4]==13 
                && b.mu_charge()->at(static_cast<int>(ll_properties[5]))==1)) {
              lplus.SetPtEtaPhiM(lep_refit0[0], lep_refit0[1], lep_refit0[2], 
                                 lep_refit0[3]);
              lminus.SetPtEtaPhiM(lep_refit0[4], lep_refit0[5], lep_refit0[6], 
                                  lep_refit0[7]);
            }
            else {
              lminus.SetPtEtaPhiM(lep_refit0[0], lep_refit0[1], lep_refit0[2], 
                                  lep_refit0[3]);
              lplus.SetPtEtaPhiM(lep_refit0[4], lep_refit0[5], lep_refit0[6], 
                                 lep_refit0[7]);
            }
          }
          ph.SetPtEtaPhiM(var_photon_pt.GetScalar(b), 
                          var_photon_eta.GetScalar(b), 
                          var_photon_phi.GetScalar(b), 
                          0.0);
          return CalculateAngles(lplus, lminus, ph);
        }).EnableCaching(true);
  }

  //kinematic angles with variations
  const NamedFunc sys_llphoton_refit_angles_default = 
      assign_variation_llphoton_refit_angles(sys_ll_default, 
      sys_lep_refit_default, sys_lead_photon_pt_default, 
      sys_lead_photon_eta_default, 
      sys_lead_photon_phi_default, "default");
  const NamedFunc sys_llphoton_refit_angles_elscaleup = 
      assign_variation_llphoton_refit_angles(sys_ll_elscaleup, 
      sys_lep_refit_elscaleup, sys_lead_photon_pt_default, 
      sys_lead_photon_eta_default, 
      sys_lead_photon_phi_default, "elscaleup");
  const NamedFunc sys_llphoton_refit_angles_elscaledn = 
      assign_variation_llphoton_refit_angles(sys_ll_elscaledn, 
      sys_lep_refit_elscaledn, sys_lead_photon_pt_default, 
      sys_lead_photon_eta_default, 
      sys_lead_photon_phi_default, "elscaledn");
  const NamedFunc sys_llphoton_refit_angles_elresup = 
      assign_variation_llphoton_refit_angles(sys_ll_elresup, 
      sys_lep_refit_elresup, sys_lead_photon_pt_default, 
      sys_lead_photon_eta_default, 
      sys_lead_photon_phi_default, "elresup");
  const NamedFunc sys_llphoton_refit_angles_elresdn = 
      assign_variation_llphoton_refit_angles(sys_ll_elresdn, 
      sys_lep_refit_elresdn, sys_lead_photon_pt_default, 
      sys_lead_photon_eta_default, 
      sys_lead_photon_phi_default, "elresdn");
  const NamedFunc sys_llphoton_refit_angles_muscaleup = 
      assign_variation_llphoton_refit_angles(sys_ll_muscaleup, 
      sys_lep_refit_muscaleup, sys_lead_photon_pt_default, 
      sys_lead_photon_eta_default, 
      sys_lead_photon_phi_default, "muscaleup");
  const NamedFunc sys_llphoton_refit_angles_muscaledn = 
      assign_variation_llphoton_refit_angles(sys_ll_muscaledn, 
      sys_lep_refit_muscaledn, sys_lead_photon_pt_default, 
      sys_lead_photon_eta_default, 
      sys_lead_photon_phi_default, "muscaledn");
  const NamedFunc sys_llphoton_refit_angles_muresup = 
      assign_variation_llphoton_refit_angles(sys_ll_muresup, 
      sys_lep_refit_muresup, sys_lead_photon_pt_default, 
      sys_lead_photon_eta_default, 
      sys_lead_photon_phi_default, "muresup");
  const NamedFunc sys_llphoton_refit_angles_muresdn = 
      assign_variation_llphoton_refit_angles(sys_ll_muresdn, 
      sys_lep_refit_muresdn, sys_lead_photon_pt_default, 
      sys_lead_photon_eta_default, 
      sys_lead_photon_phi_default, "muresdn");
  const NamedFunc sys_llphoton_refit_angles_phscaleup = 
      assign_variation_llphoton_refit_angles(sys_ll_default, 
      sys_lep_refit_default, sys_lead_photon_pt_scaleup, 
      sys_lead_photon_eta_scaleup, sys_lead_photon_phi_scaleup, "phscaleup");
  const NamedFunc sys_llphoton_refit_angles_phscaledn = 
      assign_variation_llphoton_refit_angles(sys_ll_default, 
      sys_lep_refit_default, sys_lead_photon_pt_scaledn, 
      sys_lead_photon_eta_scaledn, sys_lead_photon_phi_scaledn, "phscaledn");
  const NamedFunc sys_llphoton_refit_angles_phresup = 
      assign_variation_llphoton_refit_angles(sys_ll_default, 
      sys_lep_refit_default, sys_lead_photon_pt_resup, 
      sys_lead_photon_eta_resup, sys_lead_photon_phi_resup, "phresup");
  const NamedFunc sys_llphoton_refit_angles_phresdn = 
      assign_variation_llphoton_refit_angles(sys_ll_default, 
      sys_lep_refit_default, sys_lead_photon_pt_resdn, 
      sys_lead_photon_eta_resdn, sys_lead_photon_phi_resdn, "phresdn");
  const NamedFunc sys_llphoton_refit_cosTheta_default = ReduceNamedFuncCached(
      sys_llphoton_refit_angles_default,reduce_index0).Name(
      "sys_llphoton_refit_coscapTheta_default");
  const NamedFunc sys_llphoton_refit_cosTheta_elscaleup = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_elscaleup,reduce_index0)
      .Name("sys_llphoton_refit_coscapTheta_elscaleup");
  const NamedFunc sys_llphoton_refit_cosTheta_elscaledn = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_elscaledn,reduce_index0)
      .Name("sys_llphoton_refit_coscapTheta_elscaledn");
  const NamedFunc sys_llphoton_refit_cosTheta_elresup = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_elresup,reduce_index0)
      .Name("sys_llphoton_refit_coscapTheta_elresup");
  const NamedFunc sys_llphoton_refit_cosTheta_elresdn = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_elresdn,reduce_index0)
      .Name("sys_llphoton_refit_coscapTheta_elresdn");
  const NamedFunc sys_llphoton_refit_cosTheta_muscaleup = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_muscaleup,reduce_index0)
      .Name("sys_llphoton_refit_coscapTheta_muscaleup");
  const NamedFunc sys_llphoton_refit_cosTheta_muscaledn = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_muscaledn,reduce_index0)
      .Name("sys_llphoton_refit_coscapTheta_muscaledn");
  const NamedFunc sys_llphoton_refit_cosTheta_muresup = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_muresup,reduce_index0)
      .Name("sys_llphoton_refit_coscapTheta_muresup");
  const NamedFunc sys_llphoton_refit_cosTheta_muresdn = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_muresdn,reduce_index0)
      .Name("sys_llphoton_refit_coscapTheta_muresdn");
  const NamedFunc sys_llphoton_refit_cosTheta_phscaleup = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_phscaleup,reduce_index0)
      .Name("sys_llphoton_refit_coscapTheta_phscaleup");
  const NamedFunc sys_llphoton_refit_cosTheta_phscaledn = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_phscaledn,reduce_index0)
      .Name("sys_llphoton_refit_coscapTheta_phscaledn");
  const NamedFunc sys_llphoton_refit_cosTheta_phresup = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_phresup,reduce_index0)
      .Name("sys_llphoton_refit_coscapTheta_phresup");
  const NamedFunc sys_llphoton_refit_cosTheta_phresdn = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_phresdn,reduce_index0)
      .Name("sys_llphoton_refit_coscapTheta_phresdn");
  const NamedFunc sys_llphoton_refit_costheta_default = ReduceNamedFuncCached(
      sys_llphoton_refit_angles_default,reduce_index1).Name(
      "sys_llphoton_refit_costheta_default");
  const NamedFunc sys_llphoton_refit_costheta_elscaleup = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_elscaleup,reduce_index1)
      .Name("sys_llphoton_refit_costheta_elscaleup");
  const NamedFunc sys_llphoton_refit_costheta_elscaledn = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_elscaledn,reduce_index1)
      .Name("sys_llphoton_refit_costheta_elscaledn");
  const NamedFunc sys_llphoton_refit_costheta_elresup = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_elresup,reduce_index1)
      .Name("sys_llphoton_refit_costheta_elresup");
  const NamedFunc sys_llphoton_refit_costheta_elresdn = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_elresdn,reduce_index1)
      .Name("sys_llphoton_refit_costheta_elresdn");
  const NamedFunc sys_llphoton_refit_costheta_muscaleup = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_muscaleup,reduce_index1)
      .Name("sys_llphoton_refit_costheta_muscaleup");
  const NamedFunc sys_llphoton_refit_costheta_muscaledn = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_muscaledn,reduce_index1)
      .Name("sys_llphoton_refit_costheta_muscaledn");
  const NamedFunc sys_llphoton_refit_costheta_muresup = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_muresup,reduce_index1)
      .Name("sys_llphoton_refit_costheta_muresup");
  const NamedFunc sys_llphoton_refit_costheta_muresdn = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_muresdn,reduce_index1)
      .Name("sys_llphoton_refit_costheta_muresdn");
  const NamedFunc sys_llphoton_refit_costheta_phscaleup = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_phscaleup,reduce_index1)
      .Name("sys_llphoton_refit_costheta_phscaleup");
  const NamedFunc sys_llphoton_refit_costheta_phscaledn = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_phscaledn,reduce_index1)
      .Name("sys_llphoton_refit_costheta_phscaledn");
  const NamedFunc sys_llphoton_refit_costheta_phresup = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_phresup,reduce_index1)
      .Name("sys_llphoton_refit_costheta_phresup");
  const NamedFunc sys_llphoton_refit_costheta_phresdn = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_phresdn,reduce_index1)
      .Name("sys_llphoton_refit_costheta_phresdn");
  const NamedFunc sys_llphoton_refit_psi_default = ReduceNamedFuncCached(
      sys_llphoton_refit_angles_default,reduce_index2).Name(
      "sys_llphoton_refit_psi_default");
  const NamedFunc sys_llphoton_refit_psi_elscaleup = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_elscaleup,reduce_index2)
      .Name("sys_llphoton_refit_psi_elscaleup");
  const NamedFunc sys_llphoton_refit_psi_elscaledn = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_elscaledn,reduce_index2)
      .Name("sys_llphoton_refit_psi_elscaledn");
  const NamedFunc sys_llphoton_refit_psi_elresup = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_elresup,reduce_index2)
      .Name("sys_llphoton_refit_psi_elresup");
  const NamedFunc sys_llphoton_refit_psi_elresdn = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_elresdn,reduce_index2)
      .Name("sys_llphoton_refit_psi_elresdn");
  const NamedFunc sys_llphoton_refit_psi_muscaleup = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_muscaleup,reduce_index2)
      .Name("sys_llphoton_refit_psi_muscaleup");
  const NamedFunc sys_llphoton_refit_psi_muscaledn = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_muscaledn,reduce_index2)
      .Name("sys_llphoton_refit_psi_muscaledn");
  const NamedFunc sys_llphoton_refit_psi_muresup = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_muresup,reduce_index2)
      .Name("sys_llphoton_refit_psi_muresup");
  const NamedFunc sys_llphoton_refit_psi_muresdn = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_muresdn,reduce_index2)
      .Name("sys_llphoton_refit_psi_muresdn");
  const NamedFunc sys_llphoton_refit_psi_phscaleup = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_phscaleup,reduce_index2)
      .Name("sys_llphoton_refit_psi_phscaleup");
  const NamedFunc sys_llphoton_refit_psi_phscaledn = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_phscaledn,reduce_index2)
      .Name("sys_llphoton_refit_psi_phscaledn");
  const NamedFunc sys_llphoton_refit_psi_phresup = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_phresup,reduce_index2)
      .Name("sys_llphoton_refit_psi_phresup");
  const NamedFunc sys_llphoton_refit_psi_phresdn = 
      ReduceNamedFuncCached(sys_llphoton_refit_angles_phresdn,reduce_index2)
      .Name("sys_llphoton_refit_psi_phresdn");

  //helper function to get H candidate four vector after variation
  //doesn't consider edge cases ex. if selected Z candidate changes after
  //variation
  NamedFunc SimpleAssignVariationH(const NamedFunc &el_pt, 
                                   const NamedFunc &mu_pt,
                                   const NamedFunc &photon_pt,
                                   const std::string &name) {
    //require nll>=1 and nphoton>=0 ahead of this function to avoid issues
    return NamedFunc(("sys_llphoton_m_"+name).c_str(),[el_pt, mu_pt, photon_pt]
        (const Baby &b) -> NamedFunc::ScalarType{
          std::vector<double> b_photon_pt = photon_pt.GetVector(b);
          std::vector<double> b_el_pt = el_pt.GetVector(b);
          std::vector<double> b_mu_pt = mu_pt.GetVector(b);
          TLorentzVector l1, l2, ph;
          ph.SetPtEtaPhiM(b_photon_pt[0], b.photon_eta()->at(0), 
                          b.photon_phi()->at(0), 0.0);
          int l1_idx = b.ll_i1()->at(0);
          int l2_idx = b.ll_i2()->at(0);
          if (b.ll_lepid()->at(0)==11) {
            l1.SetPtEtaPhiM(b_el_pt[l1_idx], b.el_eta()->at(l1_idx), 
                            b.el_phi()->at(l1_idx), 0.000511);
            l2.SetPtEtaPhiM(b_el_pt[l2_idx], b.el_eta()->at(l2_idx), 
                            b.el_phi()->at(l2_idx), 0.000511);
          }
          else {
            l1.SetPtEtaPhiM(b_mu_pt[l1_idx], b.mu_eta()->at(l1_idx), 
                            b.mu_phi()->at(l1_idx), 0.106);
            l2.SetPtEtaPhiM(b_mu_pt[l2_idx], b.mu_eta()->at(l2_idx), 
                            b.mu_phi()->at(l2_idx), 0.106);
          }
          return (l1+l2+ph).M();
        });
  }

  //Gets NamedFunc that assigns dijet properties (pt eta phi m deta dphi)
  NamedFunc assign_variation_dijet(const NamedFunc &jet_sig, 
      const NamedFunc &jet_pt, const NamedFunc &jet_m, const string &name) {
    return NamedFunc(("sys_dijet_"+name).c_str(), [&jet_sig, &jet_pt, &jet_m] 
        (const Baby &b) -> NamedFunc::VectorType{
      vector<double> evt_jet_sig = jet_sig.GetVector(b);
      vector<double> evt_jet_pt = jet_pt.GetVector(b);
      vector<double> evt_jet_m = jet_m.GetVector(b);
      TLorentzVector j1, j2, dijet;
      float evt_lead_jet_pt(-999.0), evt_sublead_jet_pt(-999.0);
      for (unsigned ijet = 0; ijet < b.jet_pt()->size(); ijet++) {
        if (evt_jet_sig[ijet]) {
          if (evt_jet_pt[ijet] > evt_lead_jet_pt) {
            evt_sublead_jet_pt = evt_lead_jet_pt;
            j2 = j1;
            evt_lead_jet_pt = evt_jet_pt[ijet];
            j1.SetPtEtaPhiM(evt_jet_pt[ijet], b.jet_eta()->at(ijet), 
                            b.jet_phi()->at(ijet), evt_jet_m[ijet]);
          }
          else if (evt_jet_pt[ijet] > evt_sublead_jet_pt) {
            evt_sublead_jet_pt = evt_jet_pt[ijet];
            j2.SetPtEtaPhiM(evt_jet_pt[ijet], b.jet_eta()->at(ijet), 
                            b.jet_phi()->at(ijet), evt_jet_m[ijet]);
          }
        }
      }
      dijet = j1+j2;
      return {dijet.Pt(), dijet.Eta(), dijet.Phi(), dijet.M(), 
              fabs(j1.Eta()-j2.Eta()), deltaPhi(j1.Phi(), j2.Phi())};
    }).EnableCaching(true);
  }

  //Gets NamedFunc that is jet_pt with year_specific_variation
  NamedFunc assign_variation_jet_pt_year(const NamedFunc &var_jet_pt,
      const string year, const string &name) {
    return NamedFunc(("sys_jet_pt_"+name+year).c_str(),[var_jet_pt, year]
        (const Baby &b) -> NamedFunc::VectorType{
          if (b.SampleTypeString()==year
              || b.SampleTypeString()==("-"+year))
            return var_jet_pt.GetVector(b);
          vector<double> def_jet_pt;
          for (unsigned ijet = 0; ijet < b.jet_pt()->size(); ijet++) {
            def_jet_pt.push_back(b.jet_pt()->at(ijet));
          }
          return def_jet_pt;
        }).EnableCaching(true);
  }

  //Gets vector NamedFunc that has variation for one year
  //note this copies the supplied namedfuncs
  NamedFunc assign_vec_variation_year_select(const NamedFunc &var_namedfunc,
      const NamedFunc &def_namedfunc, const string &year, const string &name) {
    return NamedFunc((name+year).c_str(),[var_namedfunc, def_namedfunc, year]
        (const Baby &b) -> NamedFunc::VectorType{
          if (b.SampleTypeString()==year
              || b.SampleTypeString()==("-"+year))
            return var_namedfunc.GetVector(b);
          return def_namedfunc.GetVector(b); 
        }).EnableCaching(true);
  }
  
  //Gets scalar NamedFunc that has variation for one year
  //note this copies the supplied namedfuncs
  NamedFunc assign_sca_variation_year_select(const NamedFunc &var_namedfunc,
      const NamedFunc &def_namedfunc, const string &year, const string &name) {
    return NamedFunc((name+year).c_str(),[var_namedfunc, def_namedfunc, year]
        (const Baby &b) -> NamedFunc::ScalarType{
          if (b.SampleTypeString()==year
              || b.SampleTypeString()==("-"+year))
            return var_namedfunc.GetScalar(b);
          return def_namedfunc.GetScalar(b); 
        }).EnableCaching(true);
  }

  //Gets nbdfm with alternate sig
  NamedFunc assign_variation_nbdfm(const NamedFunc &var_jet_isgood, 
      const string &name) {
    return NamedFunc(("sys_nbdfm_"+name).c_str(),[&var_jet_isgood]
        (const Baby &b) -> NamedFunc::ScalarType{
        float nbdfm = 0;
        float wp = get_btag_wp_deepjet(b.SampleTypeString().Data(), 2);
        vector<double> jet_isgood = var_jet_isgood.GetVector(b);
        for (unsigned ijet = 0; ijet < b.jet_pt()->size(); ijet++) {
          if (jet_isgood[ijet] && b.jet_deepflav()->at(ijet) > wp)
            nbdfm += 1;
        }
        return nbdfm;
        }).EnableCaching(true);
  }

  //jet variations separated by era
  const NamedFunc sys_jet_pt_default("jet_pt");
  vector<NamedFunc> sys_jet_pt_scaleup;
  vector<NamedFunc> sys_jet_pt_scaledn;
  vector<NamedFunc> sys_jet_pt_resup;
  vector<NamedFunc> sys_jet_pt_resdn;
  const NamedFunc sys_jet_m_default("jet_m");
  vector<NamedFunc> sys_jet_m_scaleup;
  vector<NamedFunc> sys_jet_m_scaledn;
  vector<NamedFunc> sys_jet_m_resup;
  vector<NamedFunc> sys_jet_m_resdn;
  const NamedFunc sys_jet_isgood_default("jet_isgood");
  vector<NamedFunc> sys_jet_isgood_scaleup;
  vector<NamedFunc> sys_jet_isgood_scaledn;
  vector<NamedFunc> sys_jet_isgood_resup;
  vector<NamedFunc> sys_jet_isgood_resdn;
  const NamedFunc sys_jet_eta_default("jet_eta");
  const NamedFunc sys_jet_phi_default("jet_phi");
  const NamedFunc sys_jet_deepflav_default("jet_deepflav");
  const NamedFunc sys_sig_jet_pt_default = FilterNamedFuncCached(
      sys_jet_pt_default, sys_jet_isgood_default).Name("sig_jet_pt");
  vector<NamedFunc> sys_sig_jet_pt_scaleup;
  vector<NamedFunc> sys_sig_jet_pt_scaledn;
  vector<NamedFunc> sys_sig_jet_pt_resup;
  vector<NamedFunc> sys_sig_jet_pt_resdn;
  const NamedFunc sys_sig_jet_eta_default = FilterNamedFuncCached(
      sys_jet_eta_default, sys_jet_isgood_default).Name("sig_jet_eta");
  vector<NamedFunc> sys_sig_jet_eta_scaleup;
  vector<NamedFunc> sys_sig_jet_eta_scaledn;
  vector<NamedFunc> sys_sig_jet_eta_resup;
  vector<NamedFunc> sys_sig_jet_eta_resdn;
  const NamedFunc sys_sig_jet_phi_default = FilterNamedFuncCached(
      sys_jet_phi_default, sys_jet_isgood_default).Name("sig_jet_phi");
  vector<NamedFunc> sys_sig_jet_phi_scaleup;
  vector<NamedFunc> sys_sig_jet_phi_scaledn;
  vector<NamedFunc> sys_sig_jet_phi_resup;
  vector<NamedFunc> sys_sig_jet_phi_resdn;
  const NamedFunc sys_sig_jet_m_default = FilterNamedFuncCached(
      sys_jet_m_default, sys_jet_isgood_default).Name("sig_jet_m");
  vector<NamedFunc> sys_sig_jet_m_scaleup;
  vector<NamedFunc> sys_sig_jet_m_scaledn;
  vector<NamedFunc> sys_sig_jet_m_resup;
  vector<NamedFunc> sys_sig_jet_m_resdn;
  const NamedFunc sys_sig_jet_deepflav_default = FilterNamedFuncCached(
      sys_jet_deepflav_default, sys_jet_isgood_default).Name(
      "sig_jet_deepflav");
  vector<NamedFunc> sys_sig_jet_deepflav_scaleup;
  vector<NamedFunc> sys_sig_jet_deepflav_scaledn;
  vector<NamedFunc> sys_sig_jet_deepflav_resup;
  vector<NamedFunc> sys_sig_jet_deepflav_resdn;
  const NamedFunc sys_lead_jet_pt_default = ReduceNamedFuncCached(
      sys_sig_jet_pt_default, reduce_max).Name("lead_jet_pt");
  vector<NamedFunc> sys_lead_jet_pt_scaleup;
  vector<NamedFunc> sys_lead_jet_pt_scaledn;
  vector<NamedFunc> sys_lead_jet_pt_resup;
  vector<NamedFunc> sys_lead_jet_pt_resdn;
  const NamedFunc sys_lead_jet_eta_default = MultiReduceNamedFuncCached(
      {&sys_sig_jet_pt_default, &sys_sig_jet_eta_default}, reduce_maxfirst)
      .Name("lead_jet_eta");
  vector<NamedFunc> sys_lead_jet_eta_scaleup;
  vector<NamedFunc> sys_lead_jet_eta_scaledn;
  vector<NamedFunc> sys_lead_jet_eta_resup;
  vector<NamedFunc> sys_lead_jet_eta_resdn;
  const NamedFunc sys_lead_jet_phi_default = MultiReduceNamedFuncCached(
      {&sys_sig_jet_pt_default, &sys_sig_jet_phi_default}, reduce_maxfirst)
      .Name("lead_jet_phi");
  vector<NamedFunc> sys_lead_jet_phi_scaleup;
  vector<NamedFunc> sys_lead_jet_phi_scaledn;
  vector<NamedFunc> sys_lead_jet_phi_resup;
  vector<NamedFunc> sys_lead_jet_phi_resdn;
  const NamedFunc sys_lead_jet_m_default = MultiReduceNamedFuncCached(
      {&sys_sig_jet_pt_default, &sys_sig_jet_m_default}, reduce_maxfirst)
      .Name("lead_jet_m");
  vector<NamedFunc> sys_lead_jet_m_scaleup;
  vector<NamedFunc> sys_lead_jet_m_scaledn;
  vector<NamedFunc> sys_lead_jet_m_resup;
  vector<NamedFunc> sys_lead_jet_m_resdn;
  const NamedFunc sys_sublead_jet_pt_default = ReduceNamedFuncCached(
      sys_sig_jet_pt_default, reduce_sublead).Name("sublead_jet_pt");
  vector<NamedFunc> sys_sublead_jet_pt_scaleup;
  vector<NamedFunc> sys_sublead_jet_pt_scaledn;
  vector<NamedFunc> sys_sublead_jet_pt_resup;
  vector<NamedFunc> sys_sublead_jet_pt_resdn;
  const NamedFunc sys_sublead_jet_eta_default = MultiReduceNamedFuncCached(
      {&sys_sig_jet_pt_default, &sys_sig_jet_eta_default}, reduce_subleadfirst)
      .Name("sublead_jet_eta");
  vector<NamedFunc> sys_sublead_jet_eta_scaleup;
  vector<NamedFunc> sys_sublead_jet_eta_scaledn;
  vector<NamedFunc> sys_sublead_jet_eta_resup;
  vector<NamedFunc> sys_sublead_jet_eta_resdn;
  const NamedFunc sys_sublead_jet_phi_default = MultiReduceNamedFuncCached(
      {&sys_sig_jet_pt_default, &sys_sig_jet_phi_default}, reduce_subleadfirst)
      .Name("sublead_jet_phi");
  vector<NamedFunc> sys_sublead_jet_phi_scaleup;
  vector<NamedFunc> sys_sublead_jet_phi_scaledn;
  vector<NamedFunc> sys_sublead_jet_phi_resup;
  vector<NamedFunc> sys_sublead_jet_phi_resdn;
  const NamedFunc sys_njet_default("njet");
  vector<NamedFunc> sys_njet_scaleup;
  vector<NamedFunc> sys_njet_scaledn;
  vector<NamedFunc> sys_njet_resup;
  vector<NamedFunc> sys_njet_resdn;
  const NamedFunc sys_nbdfm_default("nbdfm");
  vector<NamedFunc> sys_nbdfm_scaleup;
  vector<NamedFunc> sys_nbdfm_scaledn;
  vector<NamedFunc> sys_nbdfm_resup;
  vector<NamedFunc> sys_nbdfm_resdn;
  const NamedFunc sys_met_default("met");
  vector<NamedFunc> sys_met_scaleup;
  vector<NamedFunc> sys_met_scaledn;
  vector<NamedFunc> sys_met_resup;
  vector<NamedFunc> sys_met_resdn;

  //dijet variations
  const NamedFunc sys_dijet_default = assign_variation_dijet(
      sys_jet_isgood_default,
      sys_jet_pt_default,sys_jet_m_default,"default");
  const NamedFunc sys_dijet_phi_default = ReduceNamedFuncCached(
      sys_dijet_default, reduce_index2).Name("sys_dijet_phi_default");
  const NamedFunc sys_dijet_m_default = ReduceNamedFuncCached(
      sys_dijet_default, reduce_index3).Name("sys_dijet_m_default");
  const NamedFunc sys_dijet_deta_default = ReduceNamedFuncCached(
      sys_dijet_default, reduce_index4).Name("sys_dijet_deta_default");
  const NamedFunc sys_dijet_dphi_default = ReduceNamedFuncCached(
      sys_dijet_default, reduce_index5).Name("sys_dijet_dphi_default");
  vector<NamedFunc> sys_dijet_scaleup;
  vector<NamedFunc> sys_dijet_scaledn;
  vector<NamedFunc> sys_dijet_resup;
  vector<NamedFunc> sys_dijet_resdn;
  vector<NamedFunc> sys_dijet_phi_scaleup;
  vector<NamedFunc> sys_dijet_phi_scaledn;
  vector<NamedFunc> sys_dijet_phi_resup;
  vector<NamedFunc> sys_dijet_phi_resdn;
  vector<NamedFunc> sys_dijet_m_scaleup;
  vector<NamedFunc> sys_dijet_m_scaledn;
  vector<NamedFunc> sys_dijet_m_resup;
  vector<NamedFunc> sys_dijet_m_resdn;
  vector<NamedFunc> sys_dijet_deta_scaleup;
  vector<NamedFunc> sys_dijet_deta_scaledn;
  vector<NamedFunc> sys_dijet_deta_resup;
  vector<NamedFunc> sys_dijet_deta_resdn;
  vector<NamedFunc> sys_dijet_dphi_scaleup;
  vector<NamedFunc> sys_dijet_dphi_scaledn;
  vector<NamedFunc> sys_dijet_dphi_resup;
  vector<NamedFunc> sys_dijet_dphi_resdn;

  //NamedFunc for 2017-2018 (i.e. horn veto, no jet veto) pinnacles production
  const NamedFunc pinnacles_run2_jet_isgood_nopt = NamedFunc("!jet_islep"
      "&&!jet_isphoton&&-4.7<jet_eta&&jet_eta<4.7&&jet_id&&!(((jet_eta>-3.0"
      "&&jet_eta<-2.5)||(jet_eta>2.5&&jet_eta<3.0))&&jet_pt<40)")
      .Name("jet_isgood_nopt").EnableCaching(true);

  //isgood variable for pinnacles
  NamedFunc assign_isgood_pinnacles(NamedFunc& var_pt, const string &name) {
    return NamedFunc(("sys_jet_isgood_"+name).c_str(),[&var_pt]
        (const Baby &b) -> NamedFunc::VectorType{
        vector<double> isgood_nopt = 
            pinnacles_run2_jet_isgood_nopt.GetVector(b);
        vector<double> var_pt_ = var_pt.GetVector(b);
        vector<double> isgood;
        for (unsigned ijet = 0; ijet < isgood_nopt.size(); ijet++) {
           if (isgood_nopt[ijet] && var_pt_[ijet])
             isgood.push_back(1.0);
           else
             isgood.push_back(0.0);
        }
        return isgood;
        }).EnableCaching(true);
  }

  //Gets NamedFunc that is untagged category selections with variation
  NamedFunc assign_variation_untagged_category(
      const NamedFunc &var_nlep, const NamedFunc &var_njet, 
      const NamedFunc &var_nbdfm, const NamedFunc &var_met,
      const NamedFunc &var_llphoton_pt, const NamedFunc &var_llphoton_m,
      const NamedFunc &var_max_lep_miniso, const NamedFunc &var_ll_m,
      const string &name) {

    return NamedFunc(("sys_untagged_category_"+name).c_str(),[&var_nlep, 
        &var_njet, &var_nbdfm, &var_met, &var_llphoton_pt, &var_llphoton_m,
        &var_max_lep_miniso, &var_ll_m]
        (const Baby &b) -> NamedFunc::ScalarType{
        double var_nlep_ = var_nlep.GetScalar(b);
        double var_njet_ = var_njet.GetScalar(b);
        double var_nbdfm_ = var_nbdfm.GetScalar(b);
        double var_met_ = var_met.GetScalar(b);
        double var_llphoton_pt_ = var_llphoton_pt.GetScalar(b);
        double var_llphoton_m_ = var_llphoton_m.GetScalar(b);
        double var_max_lep_miniso_ = var_max_lep_miniso.GetScalar(b);
        double var_ll_m_ = var_ll_m.GetScalar(b);

        return (
            //2l+b but not enough jets for tth
            (var_nlep_==2&&var_nbdfm_>=1&&var_njet_<5)
            //3l0b, but not enough met for vh3l
            ||(var_nlep_>=3&&var_nbdfm_==0.0&&var_met_<=30)
            //3l+b, but not enough jets for tth
            ||(var_nlep_==3&&var_nbdfm_>=1&&var_njet_<3)
            //fail additional llphoton pt/m selections in vhmet
            ||(var_nlep_==2&&var_njet_<=1&&var_met_>90
               &&(var_llphoton_pt_/var_llphoton_m_)<=0.4)
            //fail additional llphoton pt/m selections in vh3l
            ||(var_nlep_==3&&var_nbdfm_==0.0&&var_met_>30
               &&(var_llphoton_pt_/var_llphoton_m_)<=0.3)
            //fail additional miniso selection in vh3l
            ||(var_nlep_>=3&&var_nbdfm_==0.0&&var_met_>30
               &&var_max_lep_miniso_>0.15)
            //fail additional mll selection in tthhad
            ||(var_nlep_==2&&var_nbdfm_>=1&&var_njet_>=5
               &&(var_ll_m_<85||var_ll_m_>95))
            //fail additional miniso seleciton in tthlep
            ||(((var_nlep_==3&&var_nbdfm_>=1&&var_njet_>=3)
                ||(var_nlep_>=4&&var_nbdfm_>=1))
               &&var_max_lep_miniso_>=0.1));
        }).EnableCaching(true);

  }

  //untagged category selection with electron scale variation up
  const NamedFunc untagged_category_elscaleup = 
      assign_variation_untagged_category(
      sys_nlep_elscaleup, sys_njet_default, sys_nbdfm_default, sys_met_default,
      sys_llphoton_pt_elscaleup, 
      sys_llphoton_m_elscaleup, sys_max_lep_miniso_elscaleup, 
      sys_ll_m_elscaleup, "elscaleup");

  //untagged category selection with electron scale variation up
  const NamedFunc untagged_category_elscaledn = 
      assign_variation_untagged_category(
      sys_nlep_elscaledn, sys_njet_default, sys_nbdfm_default, 
      sys_met_default, sys_llphoton_pt_elscaledn, 
      sys_llphoton_m_elscaledn, sys_max_lep_miniso_elscaledn, 
      sys_ll_m_elscaledn, "elscaledn");

  //untagged category selection with electron resolution variation up
  const NamedFunc untagged_category_elresup = 
      assign_variation_untagged_category(
      sys_nlep_elresup, sys_njet_default, sys_nbdfm_default, 
      sys_met_default, sys_llphoton_pt_elresup, 
      sys_llphoton_m_elresup, sys_max_lep_miniso_elresup, 
      sys_ll_m_elresup, "elresup");

  //untagged category selection with electron resolution variation up
  const NamedFunc untagged_category_elresdn = 
      assign_variation_untagged_category(
      sys_nlep_elresdn, sys_njet_default, sys_nbdfm_default, sys_met_default, 
      sys_llphoton_pt_elresdn, 
      sys_llphoton_m_elresdn, sys_max_lep_miniso_elresdn, 
      sys_ll_m_elresdn, "elresdn");

  //untagged category selection with muon scale variation up
  const NamedFunc untagged_category_muscaleup = 
      assign_variation_untagged_category(
      sys_nlep_muscaleup, sys_njet_default, sys_nbdfm_default, sys_met_default,
      sys_llphoton_pt_muscaleup, 
      sys_llphoton_m_muscaleup, sys_max_lep_miniso_muscaleup, 
      sys_ll_m_muscaleup, "muscaleup");

  //untagged category selection with muon scale variation up
  const NamedFunc untagged_category_muscaledn = 
      assign_variation_untagged_category(
      sys_nlep_muscaledn, sys_njet_default, sys_nbdfm_default, sys_met_default,
      sys_llphoton_pt_muscaledn, 
      sys_llphoton_m_muscaledn, sys_max_lep_miniso_muscaledn, 
      sys_ll_m_muscaledn, "muscaledn");

  //untagged category selection with muon resolution variation up
  const NamedFunc untagged_category_muresup = 
      assign_variation_untagged_category(
      sys_nlep_muresup, sys_njet_default, sys_nbdfm_default, sys_met_default,
      sys_llphoton_pt_muresup, 
      sys_llphoton_m_muresup, sys_max_lep_miniso_muresup, 
      sys_ll_m_muresup, "muresup");

  //untagged category selection with muon resolution variation up
  const NamedFunc untagged_category_muresdn = 
      assign_variation_untagged_category(
      sys_nlep_muresdn, sys_njet_default, sys_nbdfm_default, sys_met_default,
      sys_llphoton_pt_muresdn, 
      sys_llphoton_m_muresdn, sys_max_lep_miniso_muresdn, 
      sys_ll_m_muresdn, "muresdn");

  //untagged category selection with photon scale variation up
  const NamedFunc untagged_category_phscaleup = 
      assign_variation_untagged_category(
      sys_nlep_default, sys_njet_default, sys_nbdfm_default, sys_met_default,
      sys_llphoton_pt_phscaleup, 
      sys_llphoton_m_phscaleup, sys_max_lep_miniso_default, sys_ll_m_default, 
      "phscaleup");

  //untagged category selection with photon scale variation down
  const NamedFunc untagged_category_phscaledn = 
      assign_variation_untagged_category(
      sys_nlep_default, sys_njet_default, sys_nbdfm_default, sys_met_default,
      sys_llphoton_pt_phscaledn, 
      sys_llphoton_m_phscaledn, sys_max_lep_miniso_default, sys_ll_m_default, 
      "phscaledn");

  //untagged category selection with photon resolution variation up
  const NamedFunc untagged_category_phresup = 
      assign_variation_untagged_category(
      sys_nlep_default, sys_njet_default, sys_nbdfm_default, sys_met_default,
      sys_llphoton_pt_phresup, 
      sys_llphoton_m_phresup, sys_max_lep_miniso_default, sys_ll_m_default, 
      "phresup");

  //untagged category selection with photon resolution variation down
  const NamedFunc untagged_category_phresdn = 
      assign_variation_untagged_category(
      sys_nlep_default, sys_njet_default, sys_nbdfm_default, sys_met_default,
      sys_llphoton_pt_phresdn, 
      sys_llphoton_m_phresdn, sys_max_lep_miniso_default, sys_ll_m_default, 
      "phresdn");
  
  //untagged category selection with various jet variations
  vector<NamedFunc> untagged_category_jetscaleup;
  vector<NamedFunc> untagged_category_jetscaledn;
  vector<NamedFunc> untagged_category_jetresup;
  vector<NamedFunc> untagged_category_jetresdn;

  //Gets NamedFunc that is llphoton jet dphi with variation
  NamedFunc assign_variation_llphoton_refit_jet_dphi(
      const NamedFunc &llphoton_refit_phi, const NamedFunc &jet_sig, 
      const NamedFunc &jet_pt, const string &name) {
    return NamedFunc(("sys_llphoton_refit_jet_dphi_"+name).c_str(),[
        &llphoton_refit_phi, &jet_sig, &jet_pt](const Baby &b) 
        -> NamedFunc::ScalarType{
      vector<double> evt_jet_pt = jet_pt.GetVector(b);
      vector<double> evt_jet_sig = jet_sig.GetVector(b);
      float evt_lead_jet_pt(-999.0), evt_lead_jet_phi(0.0);
      for (unsigned ijet = 0; ijet < evt_jet_sig.size(); ijet++) {
        if (evt_jet_sig[ijet] && evt_jet_pt[ijet] > evt_lead_jet_pt) {
          evt_lead_jet_pt = evt_jet_pt[ijet];
          evt_lead_jet_phi = b.jet_phi()->at(ijet);
        }
      }
      if (evt_lead_jet_pt < 0.0) return -999.0;
      return deltaPhi(llphoton_refit_phi.GetScalar(b), evt_lead_jet_phi);
    }).EnableCaching(true);
  }

  //Gets NamedFunc that is balance with variation
  NamedFunc assign_variation_llphoton_refit_jet_balance(
      const NamedFunc &ll_refit_p4, const NamedFunc &var_lead_photon_pt, 
      const NamedFunc &var_lead_photon_phi, const NamedFunc &jet_sig, 
      const NamedFunc &jet_pt, const string &name) {
    return NamedFunc(("sys_llphoton_refit_jet_balance_"+name).c_str(),[
        &ll_refit_p4, &var_lead_photon_pt, &var_lead_photon_phi, &jet_sig, 
        &jet_pt] (const Baby &b) -> NamedFunc::ScalarType{
      vector<double> evt_ll_refit_p4 = ll_refit_p4.GetVector(b);
      vector<double> evt_jet_pt = jet_pt.GetVector(b);
      vector<double> evt_jet_sig = jet_sig.GetVector(b);
      double evt_photon_pt = var_lead_photon_pt.GetScalar(b);
      double evt_photon_phi = var_lead_photon_phi.GetScalar(b);
      double mht_x(0.0), mht_y(0.0), evt_ht(0.0);
      double evt_lead_jet_pt(-999.0), evt_lead_jet_phi(0.0);
      double evt_subl_jet_pt(-999.0), evt_subl_jet_phi(0.0);
      for (unsigned ijet = 0; ijet < evt_jet_sig.size(); ijet++) {
        if (evt_jet_sig[ijet]) {
          if (evt_jet_pt[ijet] > evt_lead_jet_pt) {
            evt_subl_jet_pt = evt_lead_jet_pt;
            evt_subl_jet_phi = evt_lead_jet_phi;
            evt_lead_jet_pt = evt_jet_pt[ijet];
            evt_lead_jet_phi = b.jet_phi()->at(ijet);
          }
          else if (evt_jet_pt[ijet] > evt_subl_jet_pt) {
            evt_subl_jet_pt = evt_jet_pt[ijet];
            evt_subl_jet_phi = b.jet_phi()->at(ijet);
          }
        }
      }
      if (evt_lead_jet_pt > 0.0) {
        mht_x -= cos(evt_lead_jet_phi)*evt_lead_jet_pt;
        mht_y -= sin(evt_lead_jet_phi)*evt_lead_jet_pt;
        evt_ht += evt_lead_jet_pt;
      }
      if (evt_subl_jet_pt > 0.0) {
        mht_x -= cos(evt_subl_jet_phi)*evt_subl_jet_pt;
        mht_y -= sin(evt_subl_jet_phi)*evt_subl_jet_pt;
        evt_ht += evt_subl_jet_pt;
      }
      mht_x -= cos(evt_photon_phi)*evt_photon_pt;
      mht_y -= sin(evt_photon_phi)*evt_photon_pt;
      evt_ht += evt_photon_pt;
      mht_x -= cos(evt_ll_refit_p4[2])*evt_ll_refit_p4[0];
      mht_y -= sin(evt_ll_refit_p4[2])*evt_ll_refit_p4[0];
      evt_ht += evt_ll_refit_p4[0];
      return sqrt(mht_x*mht_x+mht_y*mht_y)/evt_ht;
    }).EnableCaching(true);
  }

  //Gets NamedFunc that is mht (mht, phi) with variation
  NamedFunc assign_variation_mht_balance(
      const NamedFunc &var_el_sig, const NamedFunc &var_el_pt,
      const NamedFunc &var_mu_sig, const NamedFunc &var_mu_pt,
      const NamedFunc &var_ph_sig, const NamedFunc &var_ph_pt,
      const NamedFunc &var_jet_sig, const NamedFunc &var_jet_pt, 
      const string &name) {
    return NamedFunc(("sys_mht_"+name).c_str(),[&var_el_sig, &var_el_pt, 
        &var_mu_sig, &var_mu_pt, &var_ph_sig, &var_ph_pt, &var_jet_sig, 
        &var_jet_pt](const Baby &b) -> NamedFunc::VectorType{
      vector<double> evt_el_pt = var_el_pt.GetVector(b);
      vector<double> evt_el_sig = var_el_sig.GetVector(b);
      vector<double> evt_mu_pt = var_mu_pt.GetVector(b);
      vector<double> evt_mu_sig = var_mu_sig.GetVector(b);
      vector<double> evt_ph_pt = var_ph_pt.GetVector(b);
      vector<double> evt_ph_sig = var_ph_sig.GetVector(b);
      vector<double> evt_jet_pt = var_jet_pt.GetVector(b);
      vector<double> evt_jet_sig = var_jet_sig.GetVector(b);
      double mht_x(0.0), mht_y(0.0);
      for (unsigned iel = 0; iel < evt_el_sig.size(); iel++) {
        if (evt_el_sig[iel]) {
          mht_x -= cos(b.el_phi()->at(iel))*evt_el_pt[iel];
          mht_y -= sin(b.el_phi()->at(iel))*evt_el_pt[iel];
        }
      }
      for (unsigned imu = 0; imu < evt_mu_sig.size(); imu++) {
        if (evt_mu_sig[imu]) {
          mht_x -= cos(b.mu_phi()->at(imu))*evt_mu_pt[imu];
          mht_y -= sin(b.mu_phi()->at(imu))*evt_mu_pt[imu];
        }
      }
      for (unsigned iph = 0; iph < evt_ph_sig.size(); iph++) {
        if (evt_ph_sig[iph]) {
          mht_x -= cos(b.photon_phi()->at(iph))*evt_ph_pt[iph];
          mht_y -= sin(b.photon_phi()->at(iph))*evt_ph_pt[iph];
        }
      }
      for (unsigned ijet = 0; ijet < evt_jet_sig.size(); ijet++) {
        if (evt_jet_sig[ijet]) {
          mht_x -= cos(b.jet_phi()->at(ijet))*evt_jet_pt[ijet];
          mht_y -= sin(b.jet_phi()->at(ijet))*evt_jet_pt[ijet];
        }
      }
      return {sqrt(mht_x*mht_x+mht_y*mht_y),atan2(mht_y, mht_x)};
    }).EnableCaching(true);
  }

  //Gets NamedFunc that is photon zeppenfeld with variation
  NamedFunc assign_variation_photon_zeppenfeld(
      const NamedFunc &var_lead_photon_eta, const NamedFunc &var_njet, 
      const NamedFunc &var_lead_jet_eta, 
      const NamedFunc &var_sublead_jet_eta, const string &name) {
    return NamedFunc(("sys_photon_zeppenfeld_"+name).c_str(),[
        &var_lead_photon_eta, &var_njet, &var_lead_jet_eta, 
        &var_sublead_jet_eta] (const Baby &b) -> NamedFunc::ScalarType{
      double evt_njet = var_njet.GetScalar(b);
      if (evt_njet < 1)
        return -999.0;
      //TODO fix 2 in next production
      if (evt_njet < 2)
        return fabs(var_lead_photon_eta.GetScalar(b)
                    -var_lead_jet_eta.GetScalar(b)/2.0);
      return fabs(var_lead_photon_eta.GetScalar(b)
          -(var_lead_jet_eta.GetScalar(b)
            +var_sublead_jet_eta.GetScalar(b))
          /2.0);
    }).EnableCaching(true);
  }


  //jet+llphoton variations
  const NamedFunc sys_llphoton_refit_dijet_dphi_default = 
      MultiMapNamedFuncCached(
      {&sys_llphoton_refit_phi_default, &sys_dijet_phi_default}, 
      map_deltaphi).Name("sys_llphoton_dijet_dphi_default");
  const NamedFunc sys_llphoton_refit_dijet_dphi_elscaleup = 
      MultiMapNamedFuncCached(
      {&sys_llphoton_refit_phi_elscaleup, &sys_dijet_phi_default}, 
      map_deltaphi).Name("sys_llphoton_dijet_dphi_elscaleup");
  const NamedFunc sys_llphoton_refit_dijet_dphi_elscaledn = 
      MultiMapNamedFuncCached(
      {&sys_llphoton_refit_phi_elscaledn, &sys_dijet_phi_default}, 
      map_deltaphi).Name("sys_llphoton_dijet_dphi_elscaledn");
  const NamedFunc sys_llphoton_refit_dijet_dphi_elresup = 
      MultiMapNamedFuncCached(
      {&sys_llphoton_refit_phi_elresup, &sys_dijet_phi_default}, 
      map_deltaphi).Name("sys_llphoton_dijet_dphi_elresup");
  const NamedFunc sys_llphoton_refit_dijet_dphi_elresdn = 
      MultiMapNamedFuncCached(
      {&sys_llphoton_refit_phi_elresdn, &sys_dijet_phi_default}, 
      map_deltaphi).Name("sys_llphoton_dijet_dphi_elresdn");
  const NamedFunc sys_llphoton_refit_dijet_dphi_muscaleup = 
      MultiMapNamedFuncCached(
      {&sys_llphoton_refit_phi_muscaleup, &sys_dijet_phi_default}, 
      map_deltaphi).Name("sys_llphoton_dijet_dphi_muscaleup");
  const NamedFunc sys_llphoton_refit_dijet_dphi_muscaledn = 
      MultiMapNamedFuncCached(
      {&sys_llphoton_refit_phi_muscaledn, &sys_dijet_phi_default}, 
      map_deltaphi).Name("sys_llphoton_dijet_dphi_muscaledn");
  const NamedFunc sys_llphoton_refit_dijet_dphi_muresup = 
      MultiMapNamedFuncCached(
      {&sys_llphoton_refit_phi_muresup, &sys_dijet_phi_default}, 
      map_deltaphi).Name("sys_llphoton_dijet_dphi_muresup");
  const NamedFunc sys_llphoton_refit_dijet_dphi_muresdn = 
      MultiMapNamedFuncCached(
      {&sys_llphoton_refit_phi_muresdn, &sys_dijet_phi_default}, 
      map_deltaphi).Name("sys_llphoton_dijet_dphi_muresdn");
  const NamedFunc sys_llphoton_refit_dijet_dphi_phscaleup = 
      MultiMapNamedFuncCached(
      {&sys_llphoton_refit_phi_phscaleup, &sys_dijet_phi_default}, 
      map_deltaphi).Name("sys_llphoton_dijet_dphi_phscaleup");
  const NamedFunc sys_llphoton_refit_dijet_dphi_phscaledn = 
      MultiMapNamedFuncCached(
      {&sys_llphoton_refit_phi_phscaledn, &sys_dijet_phi_default}, 
      map_deltaphi).Name("sys_llphoton_dijet_dphi_phscaledn");
  const NamedFunc sys_llphoton_refit_dijet_dphi_phresup = 
      MultiMapNamedFuncCached(
      {&sys_llphoton_refit_phi_phresup, &sys_dijet_phi_default}, 
      map_deltaphi).Name("sys_llphoton_dijet_dphi_phresup");
  const NamedFunc sys_llphoton_refit_dijet_dphi_phresdn = 
      MultiMapNamedFuncCached(
      {&sys_llphoton_refit_phi_phresdn, &sys_dijet_phi_default}, 
      map_deltaphi).Name("sys_llphoton_dijet_dphi_phresdn");
  vector<NamedFunc> sys_llphoton_refit_dijet_dphi_jetscaleup;
  vector<NamedFunc> sys_llphoton_refit_dijet_dphi_jetscaledn;
  vector<NamedFunc> sys_llphoton_refit_dijet_dphi_jetresup;
  vector<NamedFunc> sys_llphoton_refit_dijet_dphi_jetresdn;

  const NamedFunc sys_llphoton_refit_jet_dphi_default = 
      assign_variation_llphoton_refit_jet_dphi("llphoton_refit_phi", 
      sys_jet_isgood_default, sys_jet_pt_default, "default");
  const NamedFunc sys_llphoton_refit_jet_dphi_default_elscaleup = 
      assign_variation_llphoton_refit_jet_dphi(
      sys_llphoton_refit_phi_elscaleup, sys_jet_isgood_default, 
      sys_jet_pt_default, "elscaleup");
  const NamedFunc sys_llphoton_refit_jet_dphi_default_elscaledn = 
      assign_variation_llphoton_refit_jet_dphi(
      sys_llphoton_refit_phi_elscaledn, sys_jet_isgood_default, 
      sys_jet_pt_default, "elscaledn");
  const NamedFunc sys_llphoton_refit_jet_dphi_default_elresup = 
      assign_variation_llphoton_refit_jet_dphi(
      sys_llphoton_refit_phi_elresup, sys_jet_isgood_default, 
      sys_jet_pt_default, "elresup");
  const NamedFunc sys_llphoton_refit_jet_dphi_default_elresdn = 
      assign_variation_llphoton_refit_jet_dphi(
      sys_llphoton_refit_phi_elresdn, sys_jet_isgood_default, 
      sys_jet_pt_default, "elresdn");
  const NamedFunc sys_llphoton_refit_jet_dphi_default_muscaleup = 
      assign_variation_llphoton_refit_jet_dphi(
      sys_llphoton_refit_phi_muscaleup, sys_jet_isgood_default, 
      sys_jet_pt_default, "muscaleup");
  const NamedFunc sys_llphoton_refit_jet_dphi_default_muscaledn = 
      assign_variation_llphoton_refit_jet_dphi(
      sys_llphoton_refit_phi_muscaledn, sys_jet_isgood_default, 
      sys_jet_pt_default, "muscaledn");
  const NamedFunc sys_llphoton_refit_jet_dphi_default_muresup = 
      assign_variation_llphoton_refit_jet_dphi(
      sys_llphoton_refit_phi_muresup, sys_jet_isgood_default, 
      sys_jet_pt_default, "muresup");
  const NamedFunc sys_llphoton_refit_jet_dphi_default_muresdn = 
      assign_variation_llphoton_refit_jet_dphi(
      sys_llphoton_refit_phi_muresdn, sys_jet_isgood_default, 
      sys_jet_pt_default, "muresdn");
  const NamedFunc sys_llphoton_refit_jet_dphi_default_phscaleup = 
      assign_variation_llphoton_refit_jet_dphi(
      sys_llphoton_refit_phi_phscaleup, sys_jet_isgood_default, 
      sys_jet_pt_default, "phscaleup");
  const NamedFunc sys_llphoton_refit_jet_dphi_default_phscaledn = 
      assign_variation_llphoton_refit_jet_dphi(
      sys_llphoton_refit_phi_phscaledn, sys_jet_isgood_default, 
      sys_jet_pt_default, "phscaledn");
  const NamedFunc sys_llphoton_refit_jet_dphi_default_phresup = 
      assign_variation_llphoton_refit_jet_dphi(
      sys_llphoton_refit_phi_phresup, sys_jet_isgood_default, 
      sys_jet_pt_default, "phresup");
  const NamedFunc sys_llphoton_refit_jet_dphi_default_phresdn = 
      assign_variation_llphoton_refit_jet_dphi(
      sys_llphoton_refit_phi_phresdn, sys_jet_isgood_default, 
      sys_jet_pt_default, "phresdn");
  vector<NamedFunc> sys_llphoton_refit_jet_dphi_jetscaleup;
  vector<NamedFunc> sys_llphoton_refit_jet_dphi_jetscaledn;
  vector<NamedFunc> sys_llphoton_refit_jet_dphi_jetresup;
  vector<NamedFunc> sys_llphoton_refit_jet_dphi_jetresdn;

  const NamedFunc sys_llphoton_refit_jet_balance_default = 
      assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_default, 
      sys_lead_photon_pt_default, sys_lead_photon_phi_default, 
      sys_jet_isgood_default, sys_jet_pt_default, "default");
  const NamedFunc sys_llphoton_refit_jet_balance_elscaleup = 
      assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_elscaleup, 
      sys_lead_photon_pt_default, sys_lead_photon_phi_default, 
      sys_jet_isgood_default, sys_jet_pt_default, "elscaleup");
  const NamedFunc sys_llphoton_refit_jet_balance_elscaledn = 
      assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_elscaledn, 
      sys_lead_photon_pt_default, sys_lead_photon_phi_default, 
      sys_jet_isgood_default, sys_jet_pt_default, "elscaledn");
  const NamedFunc sys_llphoton_refit_jet_balance_elresup = 
      assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_elresup, 
      sys_lead_photon_pt_default, sys_lead_photon_phi_default, 
      sys_jet_isgood_default, sys_jet_pt_default, "elresup");
  const NamedFunc sys_llphoton_refit_jet_balance_elresdn = 
      assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_elresdn, 
      sys_lead_photon_pt_default, sys_lead_photon_phi_default, 
      sys_jet_isgood_default, sys_jet_pt_default, "elresdn");
  const NamedFunc sys_llphoton_refit_jet_balance_muscaleup = 
      assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_muscaleup, 
      sys_lead_photon_pt_default, sys_lead_photon_phi_default, 
      sys_jet_isgood_default, sys_jet_pt_default, "muscaleup");
  const NamedFunc sys_llphoton_refit_jet_balance_muscaledn = 
      assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_muscaledn, 
      sys_lead_photon_pt_default, sys_lead_photon_phi_default, 
      sys_jet_isgood_default, sys_jet_pt_default, "muscaledn");
  const NamedFunc sys_llphoton_refit_jet_balance_muresup = 
      assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_muresup, 
      sys_lead_photon_pt_default, sys_lead_photon_phi_default, 
      sys_jet_isgood_default, sys_jet_pt_default, "muresup");
  const NamedFunc sys_llphoton_refit_jet_balance_muresdn = 
      assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_muresdn, 
      sys_lead_photon_pt_default, sys_lead_photon_phi_default, 
      sys_jet_isgood_default, sys_jet_pt_default, "muresdn");
  const NamedFunc sys_llphoton_refit_jet_balance_phscaleup = 
      assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_default, 
      sys_lead_photon_pt_scaleup, sys_lead_photon_phi_scaleup, 
      sys_jet_isgood_default, 
      sys_jet_pt_default, "phscaleup");
  const NamedFunc sys_llphoton_refit_jet_balance_phscaledn = 
      assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_default, 
      sys_lead_photon_pt_scaledn, sys_lead_photon_phi_scaledn, 
      sys_jet_isgood_default, 
      sys_jet_pt_default, "phscaledn");
  const NamedFunc sys_llphoton_refit_jet_balance_phresup = 
      assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_default, 
      sys_lead_photon_pt_resup, sys_lead_photon_phi_resup, 
      sys_jet_isgood_default, 
      sys_jet_pt_default, "phresup");
  const NamedFunc sys_llphoton_refit_jet_balance_phresdn = 
      assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_default, 
      sys_lead_photon_pt_resdn, sys_lead_photon_phi_resdn, 
      sys_jet_isgood_default, 
      sys_jet_pt_default, "phresdn");
  vector<NamedFunc> sys_llphoton_refit_jet_balance_jetscaleup;
  vector<NamedFunc> sys_llphoton_refit_jet_balance_jetscaledn;
  vector<NamedFunc> sys_llphoton_refit_jet_balance_jetresup;
  vector<NamedFunc> sys_llphoton_refit_jet_balance_jetresdn;

  const NamedFunc sys_mht_default = assign_variation_mht_balance(
      sys_el_sig_default, sys_el_pt_default, sys_mu_sig_default, 
      sys_mu_pt_default, sys_photon_sig_default, sys_photon_pt_default, 
      sys_jet_isgood_default, sys_jet_pt_default, "default");
  const NamedFunc sys_mht_elscaleup = assign_variation_mht_balance(
      sys_el_sig_scaleup, sys_el_pt_scaleup, sys_mu_sig_default, 
      sys_mu_pt_default, sys_photon_sig_default, sys_photon_pt_default, 
      sys_jet_isgood_default, sys_jet_pt_default, "elscaleup");
  const NamedFunc sys_mht_elscaledn = assign_variation_mht_balance(
      sys_el_sig_scaledn, sys_el_pt_scaledn, sys_mu_sig_default, 
      sys_mu_pt_default, sys_photon_sig_default, sys_photon_pt_default, 
      sys_jet_isgood_default, sys_jet_pt_default, "elscaledn");
  const NamedFunc sys_mht_elresup = assign_variation_mht_balance(
      sys_el_sig_resup, sys_el_pt_resup, sys_mu_sig_default, 
      sys_mu_pt_default, sys_photon_sig_default, sys_photon_pt_default, 
      sys_jet_isgood_default, sys_jet_pt_default, "elresup");
  const NamedFunc sys_mht_elresdn = assign_variation_mht_balance(
      sys_el_sig_resdn, sys_el_pt_resdn, sys_mu_sig_default, 
      sys_mu_pt_default, sys_photon_sig_default, sys_photon_pt_default, 
      sys_jet_isgood_default, sys_jet_pt_default, "elresdn");
  const NamedFunc sys_mht_muscaleup = assign_variation_mht_balance(
      sys_el_sig_default, sys_el_pt_default, sys_mu_sig_scaleup, 
      sys_mu_pt_scaleup, sys_photon_sig_default, sys_photon_pt_default, 
      sys_jet_isgood_default, sys_jet_pt_default, "muscaleup");
  const NamedFunc sys_mht_muscaledn = assign_variation_mht_balance(
      sys_el_sig_default, sys_el_pt_default, sys_mu_sig_scaledn, 
      sys_mu_pt_scaledn, sys_photon_sig_default, sys_photon_pt_default, 
      sys_jet_isgood_default, sys_jet_pt_default, "muscaledn");
  const NamedFunc sys_mht_muresup = assign_variation_mht_balance(
      sys_el_sig_default, sys_el_pt_default, sys_mu_sig_resup, 
      sys_mu_pt_resup, sys_photon_sig_default, sys_photon_pt_default, 
      sys_jet_isgood_default, sys_jet_pt_default, "muresup");
  const NamedFunc sys_mht_muresdn = assign_variation_mht_balance(
      sys_el_sig_default, sys_el_pt_default, sys_mu_sig_resdn, sys_mu_pt_resdn, 
      sys_photon_sig_default, sys_photon_pt_default, sys_jet_isgood_default, 
      sys_jet_pt_default, "muresdn");
  const NamedFunc sys_mht_phscaleup = assign_variation_mht_balance(
      sys_el_sig_default, sys_el_pt_default, sys_mu_sig_default, 
      sys_mu_pt_default, sys_photon_sig_scaleup, sys_photon_pt_scaleup, 
      sys_jet_isgood_default, sys_jet_pt_default, "phscaleup");
  const NamedFunc sys_mht_phscaledn = assign_variation_mht_balance(
      sys_el_sig_default, sys_el_pt_default, sys_mu_sig_default, 
      sys_mu_pt_default, sys_photon_sig_scaledn, sys_photon_pt_scaledn, 
      sys_jet_isgood_default, sys_jet_pt_default, "phscaledn");
  const NamedFunc sys_mht_phresup = assign_variation_mht_balance(
      sys_el_sig_default, sys_el_pt_default, sys_mu_sig_default, 
      sys_mu_pt_default, sys_photon_sig_resup, sys_photon_pt_resup, 
      sys_jet_isgood_default, sys_jet_pt_default, "phresup");
  const NamedFunc sys_mht_phresdn = assign_variation_mht_balance(
      sys_el_sig_default, sys_el_pt_default, sys_mu_sig_default, 
      sys_mu_pt_default, sys_photon_sig_resdn, sys_photon_pt_resdn, 
      sys_jet_isgood_default, sys_jet_pt_default, "phresdn");
  vector<NamedFunc> sys_mht_jetscaleup;
  vector<NamedFunc> sys_mht_jetscaledn;
  vector<NamedFunc> sys_mht_jetresup;
  vector<NamedFunc> sys_mht_jetresdn;

  const NamedFunc sys_mht_phi_default = ReduceNamedFuncCached(sys_mht_default,
      reduce_index1).Name("sys_mht_phi_default");
  const NamedFunc sys_mht_phi_elscaleup = ReduceNamedFuncCached(
      sys_mht_elscaleup, reduce_index1).Name("sys_mht_phi_elscaleup");
  const NamedFunc sys_mht_phi_elscaledn = ReduceNamedFuncCached(
      sys_mht_elscaledn, reduce_index1).Name("sys_mht_phi_elscaledn");
  const NamedFunc sys_mht_phi_elresup = ReduceNamedFuncCached(
      sys_mht_elresup, reduce_index1).Name("sys_mht_phi_elresup");
  const NamedFunc sys_mht_phi_elresdn = ReduceNamedFuncCached(
      sys_mht_elresdn, reduce_index1).Name("sys_mht_phi_elresdn");
  const NamedFunc sys_mht_phi_muscaleup = ReduceNamedFuncCached(
      sys_mht_muscaleup, reduce_index1).Name("sys_mht_phi_muscaleup");
  const NamedFunc sys_mht_phi_muscaledn = ReduceNamedFuncCached(
      sys_mht_muscaledn, reduce_index1).Name("sys_mht_phi_muscaledn");
  const NamedFunc sys_mht_phi_muresup = ReduceNamedFuncCached(
      sys_mht_muresup, reduce_index1).Name("sys_mht_phi_muresup");
  const NamedFunc sys_mht_phi_muresdn = ReduceNamedFuncCached(
      sys_mht_muresdn, reduce_index1).Name("sys_mht_phi_muresdn");
  const NamedFunc sys_mht_phi_phscaleup = ReduceNamedFuncCached(
      sys_mht_phscaleup, reduce_index1).Name("sys_mht_phi_phscaleup");
  const NamedFunc sys_mht_phi_phscaledn = ReduceNamedFuncCached(
      sys_mht_phscaledn, reduce_index1).Name("sys_mht_phi_phscaledn");
  const NamedFunc sys_mht_phi_phresup = ReduceNamedFuncCached(
      sys_mht_phresup, reduce_index1).Name("sys_mht_phi_phresup");
  const NamedFunc sys_mht_phi_phresdn = ReduceNamedFuncCached(
      sys_mht_phresdn, reduce_index1).Name("sys_mht_phi_phresdn");
  vector<NamedFunc> sys_mht_phi_jetscaleup;
  vector<NamedFunc> sys_mht_phi_jetscaledn;
  vector<NamedFunc> sys_mht_phi_jetresup;
  vector<NamedFunc> sys_mht_phi_jetresdn;

  const NamedFunc sys_photon_mht_dphi_default = MultiMapNamedFuncCached(
      {&sys_lead_photon_phi_default, &sys_mht_phi_default}, map_deltaphi)
      .Name("sys_photon_mht_dphi_default");
  const NamedFunc sys_photon_mht_dphi_elscaleup = MultiMapNamedFuncCached(
      {&sys_lead_photon_phi_default, &sys_mht_phi_elscaleup}, map_deltaphi)
      .Name("sys_photon_mht_dphi_elscaleup");
  const NamedFunc sys_photon_mht_dphi_elscaledn = MultiMapNamedFuncCached(
      {&sys_lead_photon_phi_default, &sys_mht_phi_elscaledn}, map_deltaphi)
      .Name("sys_photon_mht_dphi_elscaledn");
  const NamedFunc sys_photon_mht_dphi_elresup = MultiMapNamedFuncCached(
      {&sys_lead_photon_phi_default, &sys_mht_phi_elresup}, map_deltaphi)
      .Name("sys_photon_mht_dphi_elresup");
  const NamedFunc sys_photon_mht_dphi_elresdn = MultiMapNamedFuncCached(
      {&sys_lead_photon_phi_default, &sys_mht_phi_elresdn}, map_deltaphi)
      .Name("sys_photon_mht_dphi_elresdn");
  const NamedFunc sys_photon_mht_dphi_muscaleup = MultiMapNamedFuncCached(
      {&sys_lead_photon_phi_default, &sys_mht_phi_muscaleup}, map_deltaphi)
      .Name("sys_photon_mht_dphi_muscaleup");
  const NamedFunc sys_photon_mht_dphi_muscaledn = MultiMapNamedFuncCached(
      {&sys_lead_photon_phi_default, &sys_mht_phi_muscaledn}, map_deltaphi)
      .Name("sys_photon_mht_dphi_muscaledn");
  const NamedFunc sys_photon_mht_dphi_muresup = MultiMapNamedFuncCached(
      {&sys_lead_photon_phi_default, &sys_mht_phi_muresup}, map_deltaphi)
      .Name("sys_photon_mht_dphi_muresup");
  const NamedFunc sys_photon_mht_dphi_muresdn = MultiMapNamedFuncCached(
      {&sys_lead_photon_phi_default, &sys_mht_phi_muresdn}, map_deltaphi)
      .Name("sys_photon_mht_dphi_muresdn");
  const NamedFunc sys_photon_mht_dphi_phscaleup = MultiMapNamedFuncCached(
      {&sys_lead_photon_phi_scaleup, &sys_mht_phi_phscaleup}, map_deltaphi)
      .Name("sys_photon_mht_dphi_phscaleup");
  const NamedFunc sys_photon_mht_dphi_phscaledn = MultiMapNamedFuncCached(
      {&sys_lead_photon_phi_scaledn, &sys_mht_phi_phscaledn}, map_deltaphi)
      .Name("sys_photon_mht_dphi_phscaledn");
  const NamedFunc sys_photon_mht_dphi_phresup = MultiMapNamedFuncCached(
      {&sys_lead_photon_phi_resup, &sys_mht_phi_phresup}, map_deltaphi)
      .Name("sys_photon_mht_dphi_phresup");
  const NamedFunc sys_photon_mht_dphi_phresdn = MultiMapNamedFuncCached(
      {&sys_lead_photon_phi_resdn, &sys_mht_phi_phresdn}, map_deltaphi)
      .Name("sys_photon_mht_dphi_phresdn");
  vector<NamedFunc> sys_photon_mht_dphi_jetscaleup;
  vector<NamedFunc> sys_photon_mht_dphi_jetscaledn;
  vector<NamedFunc> sys_photon_mht_dphi_jetresup;
  vector<NamedFunc> sys_photon_mht_dphi_jetresdn;

  const NamedFunc sys_photon_jet1_dr_default = MultiMapNamedFuncCached(
      {&sys_lead_photon_eta_default, &sys_lead_photon_phi_default, 
      &sys_lead_jet_eta_default, &sys_lead_jet_phi_default}, map_deltar).Name(
      "sys_photon_jet1_dr_default");
  const NamedFunc sys_photon_jet1_dr_phscaleup = MultiMapNamedFuncCached(
      {&sys_lead_photon_eta_scaleup, &sys_lead_photon_phi_scaleup, 
      &sys_lead_jet_eta_default, &sys_lead_jet_phi_default}, map_deltar).Name(
      "sys_photon_jet1_dr_phscaleup");
  const NamedFunc sys_photon_jet1_dr_phscaledn = MultiMapNamedFuncCached(
      {&sys_lead_photon_eta_scaledn, &sys_lead_photon_phi_scaledn, 
      &sys_lead_jet_eta_default, &sys_lead_jet_phi_default}, map_deltar).Name(
      "sys_photon_jet1_dr_phscaledn");
  const NamedFunc sys_photon_jet1_dr_phresup = MultiMapNamedFuncCached(
      {&sys_lead_photon_eta_resup, &sys_lead_photon_phi_resup, 
      &sys_lead_jet_eta_default, &sys_lead_jet_phi_default}, map_deltar).Name(
      "sys_photon_jet1_dr_phresup");
  const NamedFunc sys_photon_jet1_dr_phresdn = MultiMapNamedFuncCached(
      {&sys_lead_photon_eta_resdn, &sys_lead_photon_phi_resdn, 
      &sys_lead_jet_eta_default, &sys_lead_jet_phi_default}, map_deltar).Name(
      "sys_photon_jet1_dr_phresdn");
  vector<NamedFunc> sys_photon_jet1_dr_jetscaleup;
  vector<NamedFunc> sys_photon_jet1_dr_jetscaledn;
  vector<NamedFunc> sys_photon_jet1_dr_jetresup;
  vector<NamedFunc> sys_photon_jet1_dr_jetresdn;

  const NamedFunc sys_photon_jet2_dr_default = MultiMapNamedFuncCached(
      {&sys_lead_photon_eta_default, &sys_lead_photon_phi_default, 
      &sys_sublead_jet_eta_default, &sys_sublead_jet_phi_default}, map_deltar)
      .Name("sys_photon_jet2_dr_default");
  const NamedFunc sys_photon_jet2_dr_phscaleup = MultiMapNamedFuncCached(
      {&sys_lead_photon_eta_scaleup, &sys_lead_photon_phi_scaleup, 
      &sys_sublead_jet_eta_default, &sys_sublead_jet_phi_default}, map_deltar)
      .Name("sys_photon_jet2_dr_phscaleup");
  const NamedFunc sys_photon_jet2_dr_phscaledn = MultiMapNamedFuncCached(
      {&sys_lead_photon_eta_scaledn, &sys_lead_photon_phi_scaledn, 
      &sys_sublead_jet_eta_default, &sys_sublead_jet_phi_default}, map_deltar)
      .Name("sys_photon_jet2_dr_phscaledn");
  const NamedFunc sys_photon_jet2_dr_phresup = MultiMapNamedFuncCached(
      {&sys_lead_photon_eta_resup, &sys_lead_photon_phi_resup, 
      &sys_sublead_jet_eta_default, &sys_sublead_jet_phi_default}, map_deltar)
      .Name("sys_photon_jet2_dr_phresup");
  const NamedFunc sys_photon_jet2_dr_phresdn = MultiMapNamedFuncCached(
      {&sys_lead_photon_eta_resdn, &sys_lead_photon_phi_resdn, 
      &sys_sublead_jet_eta_default, &sys_sublead_jet_phi_default}, map_deltar)
      .Name("sys_photon_jet2_dr_phresdn");
  vector<NamedFunc> sys_photon_jet2_dr_jetscaleup;
  vector<NamedFunc> sys_photon_jet2_dr_jetscaledn;
  vector<NamedFunc> sys_photon_jet2_dr_jetresup;
  vector<NamedFunc> sys_photon_jet2_dr_jetresdn;

  const NamedFunc sys_photon_zeppenfeld_default = 
      assign_variation_photon_zeppenfeld(sys_lead_photon_eta_default,
      sys_njet_default, sys_lead_jet_eta_default, sys_sublead_jet_eta_default, 
      "default");
  const NamedFunc sys_photon_zeppenfeld_phscaleup = 
      assign_variation_photon_zeppenfeld(sys_lead_photon_eta_scaleup,
      sys_njet_default, sys_lead_jet_eta_default, sys_sublead_jet_eta_default, 
      "phscaleup");
  const NamedFunc sys_photon_zeppenfeld_phscaledn = 
      assign_variation_photon_zeppenfeld(sys_lead_photon_eta_scaledn,
      sys_njet_default, sys_lead_jet_eta_default, sys_sublead_jet_eta_default, 
      "phscaledn");
  const NamedFunc sys_photon_zeppenfeld_phresup = 
      assign_variation_photon_zeppenfeld(sys_lead_photon_eta_resup,
      sys_njet_default, sys_lead_jet_eta_default, sys_sublead_jet_eta_default, 
      "phresup");
  const NamedFunc sys_photon_zeppenfeld_phresdn = 
      assign_variation_photon_zeppenfeld(sys_lead_photon_eta_resdn,
      sys_njet_default, sys_lead_jet_eta_default, sys_sublead_jet_eta_default, 
      "phresdn");
  vector<NamedFunc> sys_photon_zeppenfeld_jetscaleup;
  vector<NamedFunc> sys_photon_zeppenfeld_jetscaledn;
  vector<NamedFunc> sys_photon_zeppenfeld_jetresup;
  vector<NamedFunc> sys_photon_zeppenfeld_jetresdn;

  //btag uncorrelated variations
  vector<NamedFunc> sys_bchig_uncorr_up;
  vector<NamedFunc> sys_bchig_uncorr_dn;
  vector<NamedFunc> sys_udsghig_uncorr_up;
  vector<NamedFunc> sys_udsghig_uncorr_dn;

  //isgood requires !invetomap !inhemveto !inetahornveto !islep !isphoton eta 
  //fixedjetid pt
  //isgood_min includes !islep !isphoton eta fixedjetid !inethornveto
  //but not pt, !invetomap, or !inhemveto

  //horrible hack to get around the fact that I can't put a for loop outside
  //of a function to initialize what should be constant vectors...
  void initialize_jetvariations() {
    const vector<string> years = {"2016APV", "2016", "2017", "2018", 
                                  "2022", "2022EE", "2023", "2023BPix"};
    sys_jet_pt_scaleup.reserve(years.size());
    sys_jet_pt_scaledn.reserve(years.size());
    sys_jet_pt_resup.reserve(years.size());
    sys_jet_pt_resdn.reserve(years.size());
    sys_jet_m_scaleup.reserve(years.size());
    sys_jet_m_scaledn.reserve(years.size());
    sys_jet_m_resup.reserve(years.size());
    sys_jet_m_resdn.reserve(years.size());
    sys_jet_isgood_scaleup.reserve(years.size());
    sys_jet_isgood_scaledn.reserve(years.size());
    sys_jet_isgood_resup.reserve(years.size());
    sys_jet_isgood_resdn.reserve(years.size());
    sys_sig_jet_pt_scaleup.reserve(years.size());
    sys_sig_jet_pt_scaledn.reserve(years.size());
    sys_sig_jet_pt_resup.reserve(years.size());
    sys_sig_jet_pt_resdn.reserve(years.size());
    sys_sig_jet_eta_scaleup.reserve(years.size());
    sys_sig_jet_eta_scaledn.reserve(years.size());
    sys_sig_jet_eta_resup.reserve(years.size());
    sys_sig_jet_eta_resdn.reserve(years.size());
    sys_sig_jet_phi_scaleup.reserve(years.size());
    sys_sig_jet_phi_scaledn.reserve(years.size());
    sys_sig_jet_phi_resup.reserve(years.size());
    sys_sig_jet_phi_resdn.reserve(years.size());
    sys_sig_jet_m_scaleup.reserve(years.size());
    sys_sig_jet_m_scaledn.reserve(years.size());
    sys_sig_jet_m_resup.reserve(years.size());
    sys_sig_jet_m_resdn.reserve(years.size());
    sys_sig_jet_deepflav_scaleup.reserve(years.size());
    sys_sig_jet_deepflav_scaledn.reserve(years.size());
    sys_sig_jet_deepflav_resup.reserve(years.size());
    sys_sig_jet_deepflav_resdn.reserve(years.size());
    sys_lead_jet_pt_scaleup.reserve(years.size());
    sys_lead_jet_pt_scaledn.reserve(years.size());
    sys_lead_jet_pt_resup.reserve(years.size());
    sys_lead_jet_pt_resdn.reserve(years.size());
    sys_lead_jet_eta_scaleup.reserve(years.size());
    sys_lead_jet_eta_scaledn.reserve(years.size());
    sys_lead_jet_eta_resup.reserve(years.size());
    sys_lead_jet_eta_resdn.reserve(years.size());
    sys_lead_jet_phi_scaleup.reserve(years.size());
    sys_lead_jet_phi_scaledn.reserve(years.size());
    sys_lead_jet_phi_resup.reserve(years.size());
    sys_lead_jet_phi_resdn.reserve(years.size());
    sys_lead_jet_m_scaleup.reserve(years.size());
    sys_lead_jet_m_scaledn.reserve(years.size());
    sys_lead_jet_m_resup.reserve(years.size());
    sys_lead_jet_m_resdn.reserve(years.size());
    sys_sublead_jet_pt_scaleup.reserve(years.size());
    sys_sublead_jet_pt_scaledn.reserve(years.size());
    sys_sublead_jet_pt_resup.reserve(years.size());
    sys_sublead_jet_pt_resdn.reserve(years.size());
    sys_sublead_jet_eta_scaleup.reserve(years.size());
    sys_sublead_jet_eta_scaledn.reserve(years.size());
    sys_sublead_jet_eta_resup.reserve(years.size());
    sys_sublead_jet_eta_resdn.reserve(years.size());
    sys_sublead_jet_phi_scaleup.reserve(years.size());
    sys_sublead_jet_phi_scaledn.reserve(years.size());
    sys_sublead_jet_phi_resup.reserve(years.size());
    sys_sublead_jet_phi_resdn.reserve(years.size());
    sys_njet_scaleup.reserve(years.size());
    sys_njet_scaledn.reserve(years.size());
    sys_njet_resup.reserve(years.size());
    sys_njet_resdn.reserve(years.size());
    sys_nbdfm_scaleup.reserve(years.size());
    sys_nbdfm_scaledn.reserve(years.size());
    sys_nbdfm_resup.reserve(years.size());
    sys_nbdfm_resdn.reserve(years.size());
    sys_met_scaleup.reserve(years.size());
    sys_met_scaledn.reserve(years.size());
    sys_met_resup.reserve(years.size());
    sys_met_resdn.reserve(years.size());
    sys_dijet_scaleup.reserve(years.size());
    sys_dijet_scaledn.reserve(years.size());
    sys_dijet_resup.reserve(years.size());
    sys_dijet_resdn.reserve(years.size());
    sys_dijet_phi_scaleup.reserve(years.size());
    sys_dijet_phi_scaledn.reserve(years.size());
    sys_dijet_phi_resup.reserve(years.size());
    sys_dijet_phi_resdn.reserve(years.size());
    sys_dijet_m_scaleup.reserve(years.size());
    sys_dijet_m_scaledn.reserve(years.size());
    sys_dijet_m_resup.reserve(years.size());
    sys_dijet_m_resdn.reserve(years.size());
    sys_dijet_deta_scaleup.reserve(years.size());
    sys_dijet_deta_scaledn.reserve(years.size());
    sys_dijet_deta_resup.reserve(years.size());
    sys_dijet_deta_resdn.reserve(years.size());
    sys_dijet_dphi_scaleup.reserve(years.size());
    sys_dijet_dphi_scaledn.reserve(years.size());
    sys_dijet_dphi_resup.reserve(years.size());
    sys_dijet_dphi_resdn.reserve(years.size());
    untagged_category_jetscaleup.reserve(years.size());
    untagged_category_jetscaledn.reserve(years.size());
    untagged_category_jetresup.reserve(years.size());
    untagged_category_jetresdn.reserve(years.size());
    sys_llphoton_refit_dijet_dphi_jetscaleup.reserve(years.size());
    sys_llphoton_refit_dijet_dphi_jetscaledn.reserve(years.size());
    sys_llphoton_refit_dijet_dphi_jetresup.reserve(years.size());
    sys_llphoton_refit_dijet_dphi_jetresdn.reserve(years.size());
    sys_llphoton_refit_jet_dphi_jetscaleup.reserve(years.size());
    sys_llphoton_refit_jet_dphi_jetscaledn.reserve(years.size());
    sys_llphoton_refit_jet_dphi_jetresup.reserve(years.size());
    sys_llphoton_refit_jet_dphi_jetresdn.reserve(years.size());
    sys_llphoton_refit_jet_balance_jetscaleup.reserve(years.size());
    sys_llphoton_refit_jet_balance_jetscaledn.reserve(years.size());
    sys_llphoton_refit_jet_balance_jetresup.reserve(years.size());
    sys_llphoton_refit_jet_balance_jetresdn.reserve(years.size());
    sys_mht_jetscaleup.reserve(years.size());
    sys_mht_jetscaledn.reserve(years.size());
    sys_mht_jetresup.reserve(years.size());
    sys_mht_jetresdn.reserve(years.size());
    sys_mht_phi_jetscaleup.reserve(years.size());
    sys_mht_phi_jetscaledn.reserve(years.size());
    sys_mht_phi_jetresup.reserve(years.size());
    sys_mht_phi_jetresdn.reserve(years.size());
    sys_photon_mht_dphi_jetscaleup.reserve(years.size());
    sys_photon_mht_dphi_jetscaledn.reserve(years.size());
    sys_photon_mht_dphi_jetresup.reserve(years.size());
    sys_photon_mht_dphi_jetresdn.reserve(years.size());
    sys_photon_jet1_dr_jetscaleup.reserve(years.size());
    sys_photon_jet1_dr_jetscaledn.reserve(years.size());
    sys_photon_jet1_dr_jetresup.reserve(years.size());
    sys_photon_jet1_dr_jetresdn.reserve(years.size());
    sys_photon_jet2_dr_jetscaleup.reserve(years.size());
    sys_photon_jet2_dr_jetscaledn.reserve(years.size());
    sys_photon_jet2_dr_jetresup.reserve(years.size());
    sys_photon_jet2_dr_jetresdn.reserve(years.size());
    sys_photon_zeppenfeld_jetscaleup.reserve(years.size());
    sys_photon_zeppenfeld_jetscaledn.reserve(years.size());
    sys_photon_zeppenfeld_jetresup.reserve(years.size());
    sys_photon_zeppenfeld_jetresdn.reserve(years.size());
    sys_bchig_uncorr_up.reserve(years.size());
    sys_bchig_uncorr_dn.reserve(years.size());
    sys_udsghig_uncorr_up.reserve(years.size());
    sys_udsghig_uncorr_dn.reserve(years.size());
    for (unsigned iyear = 0; iyear < years.size(); iyear++) {
      string year = years[iyear];
      sys_jet_pt_scaleup.push_back(assign_vec_variation_year_select(
          "sys_jet_pt_jesup", "jet_pt", year, "sys_jet_pt_scaleup"));
      sys_jet_pt_scaledn.push_back(assign_vec_variation_year_select(
          "sys_jet_pt_jesdn", "jet_pt", year, "sys_jet_pt_scaledn"));
      sys_jet_pt_resup.push_back(assign_vec_variation_year_select(
          "sys_jet_pt_jerup", "jet_pt", year, "sys_jet_pt_resup"));
      sys_jet_pt_resdn.push_back(assign_vec_variation_year_select(
          "sys_jet_pt_jerdn", "jet_pt", year, "sys_jet_pt_resdn"));
      sys_jet_m_scaleup.push_back(assign_vec_variation_year_select(
          "sys_jet_m_jesup", "jet_m", year, "sys_jet_m_scaleup"));
      sys_jet_m_scaledn.push_back(assign_vec_variation_year_select(
          "sys_jet_m_jesdn", "jet_m", year, "sys_jet_m_scaledn"));
      sys_jet_m_resup.push_back(assign_vec_variation_year_select(
          "sys_jet_m_jerup", "jet_m", year, "sys_jet_m_resup"));
      sys_jet_m_resdn.push_back(assign_vec_variation_year_select(
          "sys_jet_m_jerdn", "jet_m", year, "sys_jet_m_resdn"));
      //TODO uncomment real ones for next production
      sys_jet_isgood_scaleup.push_back(assign_isgood_pinnacles(
          sys_jet_pt_scaleup[iyear], "scaleup"+year));
      sys_jet_isgood_scaledn.push_back(assign_isgood_pinnacles(
          sys_jet_pt_scaledn[iyear], "scaledn"+year));
      sys_jet_isgood_resup.push_back(assign_isgood_pinnacles(
          sys_jet_pt_resup[iyear], "resup"+year));
      sys_jet_isgood_resdn.push_back(assign_isgood_pinnacles(
          sys_jet_pt_resdn[iyear], "resdn"+year));
      //sys_jet_isgood_scaleup.push_back(NamedFunc(
      //    "jet_isgood_min&&!jet_isvetomap&&!jet_is_vetohem"
      //    &&sys_jet_pt_scaleup[iyear]>30.0).Name("sys_jet_isgood_scaleup"+year)
      //    .EnableCaching(true));
      //sys_jet_isgood_scaledn.push_back(NamedFunc(
      //    "jet_isgood_min&&!jet_isvetomap&&!jet_is_vetohem"
      //    &&sys_jet_pt_scaledn[iyear]>30.0).Name("sys_jet_isgood_scaledn"+year)
      //    .EnableCaching(true));
      //sys_jet_isgood_resup.push_back(NamedFunc(
      //    "jet_isgood_min&&!jet_isvetomap&&!jet_is_vetohem"
      //    &&sys_jet_pt_resup[iyear]>30.0).Name("sys_jet_isgood_resup"+year)
      //    .EnableCaching(true));
      //sys_jet_isgood_resdn.push_back(NamedFunc(
      //    "jet_isgood_min&&!jet_isvetomap&&!jet_is_vetohem"
      //    &&sys_jet_pt_resdn[iyear]>30.0).Name("sys_jet_isgood_resdn"+year)
      //    .EnableCaching(true));
      sys_sig_jet_pt_scaleup.push_back(FilterNamedFuncCached(
          sys_jet_pt_scaleup[iyear], sys_jet_isgood_scaleup[iyear]).Name(
          "sys_sig_jet_pt_scaleup"+year));
      sys_sig_jet_pt_scaledn.push_back(FilterNamedFuncCached(
          sys_jet_pt_scaledn[iyear], sys_jet_isgood_scaledn[iyear]).Name(
          "sys_sig_jet_pt_scaledn"+year));
      sys_sig_jet_pt_resup.push_back(FilterNamedFuncCached(
          sys_jet_pt_resup[iyear],
          sys_jet_isgood_resup[iyear]).Name("sys_sig_jet_pt_resup"+year));
      sys_sig_jet_pt_resdn.push_back(FilterNamedFuncCached(
          sys_jet_pt_resdn[iyear],
          sys_jet_isgood_resdn[iyear]).Name("sys_sig_jet_pt_resdn"+year));
      sys_sig_jet_eta_scaleup.push_back(FilterNamedFuncCached(
          sys_jet_eta_default,
          sys_jet_isgood_scaleup[iyear]).Name("sys_sig_jet_eta_scaleup"+year));
      sys_sig_jet_eta_scaledn.push_back(FilterNamedFuncCached(
          sys_jet_eta_default,
          sys_jet_isgood_scaledn[iyear]).Name("sys_sig_jet_eta_scaledn"+year));
      sys_sig_jet_eta_resup.push_back(FilterNamedFuncCached(
          sys_jet_eta_default,
          sys_jet_isgood_resup[iyear]).Name("sys_sig_jet_eta_resup"+year));
      sys_sig_jet_eta_resdn.push_back(FilterNamedFuncCached(
          sys_jet_eta_default,
          sys_jet_isgood_resdn[iyear]).Name("sys_sig_jet_eta_resdn"+year));
      sys_sig_jet_phi_scaleup.push_back(FilterNamedFuncCached(
          sys_jet_phi_default,
          sys_jet_isgood_scaleup[iyear]).Name("sys_sig_jet_phi_scaleup"+year));
      sys_sig_jet_phi_scaledn.push_back(FilterNamedFuncCached(
          sys_jet_phi_default,
          sys_jet_isgood_scaledn[iyear]).Name("sys_sig_jet_phi_scaledn"+year));
      sys_sig_jet_phi_resup.push_back(FilterNamedFuncCached(
          sys_jet_phi_default,
          sys_jet_isgood_resup[iyear]).Name("sys_sig_jet_phi_resup"+year));
      sys_sig_jet_phi_resdn.push_back(FilterNamedFuncCached(
          sys_jet_phi_default,
          sys_jet_isgood_resdn[iyear]).Name("sys_sig_jet_phi_resdn"+year));
      sys_sig_jet_m_scaleup.push_back(FilterNamedFuncCached(sys_jet_m_default,
          sys_jet_isgood_scaleup[iyear]).Name("sys_sig_jet_m_scaleup"+year));
      sys_sig_jet_m_scaledn.push_back(FilterNamedFuncCached(sys_jet_m_default,
          sys_jet_isgood_scaledn[iyear]).Name("sys_sig_jet_m_scaledn"+year));
      sys_sig_jet_m_resup.push_back(FilterNamedFuncCached(sys_jet_m_default,
          sys_jet_isgood_resup[iyear]).Name("sys_sig_jet_m_resup"+year));
      sys_sig_jet_m_resdn.push_back(FilterNamedFuncCached(sys_jet_m_default,
          sys_jet_isgood_resdn[iyear]).Name("sys_sig_jet_m_resdn"+year));
      sys_lead_jet_pt_scaleup.push_back(ReduceNamedFuncCached(
          sys_sig_jet_pt_scaleup[iyear],reduce_max).Name(
          "sys_lead_jet_pt_scaleup"+year));
      sys_lead_jet_pt_scaledn.push_back(ReduceNamedFuncCached(
          sys_sig_jet_pt_scaledn[iyear],reduce_max).Name(
          "sys_lead_jet_pt_scaledn"+year));
      sys_lead_jet_pt_resup.push_back(ReduceNamedFuncCached(
          sys_sig_jet_pt_resup[iyear],reduce_max).Name(
          "sys_lead_jet_pt_resup"+year));
      sys_lead_jet_pt_resdn.push_back(ReduceNamedFuncCached(
          sys_sig_jet_pt_resdn[iyear],reduce_max).Name(
          "sys_lead_jet_pt_resdn"+year));
      sys_sublead_jet_pt_scaleup.push_back(ReduceNamedFuncCached(
          sys_sig_jet_pt_scaleup[iyear],reduce_sublead).Name(
          "sys_sublead_jet_pt_scaleup"+year));
      sys_sublead_jet_pt_scaledn.push_back(ReduceNamedFuncCached(
          sys_sig_jet_pt_scaledn[iyear],reduce_sublead).Name(
          "sys_sublead_jet_pt_scaledn"+year));
      sys_sublead_jet_pt_resup.push_back(ReduceNamedFuncCached(
          sys_sig_jet_pt_resup[iyear],reduce_sublead).Name(
          "sys_sublead_jet_pt_resup"+year));
      sys_sublead_jet_pt_resdn.push_back(ReduceNamedFuncCached(
          sys_sig_jet_pt_resdn[iyear],reduce_sublead).Name(
          "sys_sublead_jet_pt_resdn"+year));
      sys_lead_jet_eta_scaleup.push_back(MultiReduceNamedFuncCached(
          {&sys_sig_jet_pt_scaleup[iyear],&sys_sig_jet_eta_scaleup[iyear]}, 
          reduce_maxfirst).Name("sys_lead_jet_eta_scaleup"+year));
      sys_lead_jet_eta_scaledn.push_back(MultiReduceNamedFuncCached(
          {&sys_sig_jet_pt_scaledn[iyear],&sys_sig_jet_eta_scaledn[iyear]}, 
          reduce_maxfirst).Name("sys_lead_jet_eta_scaledn"+year));
      sys_lead_jet_eta_resup.push_back(MultiReduceNamedFuncCached(
          {&sys_sig_jet_pt_resup[iyear],&sys_sig_jet_eta_resup[iyear]}, 
          reduce_maxfirst).Name("sys_lead_jet_eta_resup"+year));
      sys_lead_jet_eta_resdn.push_back(MultiReduceNamedFuncCached(
          {&sys_sig_jet_pt_resdn[iyear],&sys_sig_jet_eta_resdn[iyear]}, 
          reduce_maxfirst).Name("sys_lead_jet_eta_resdn"+year));
      sys_lead_jet_phi_scaleup.push_back(MultiReduceNamedFuncCached(
          {&sys_sig_jet_pt_scaleup[iyear],&sys_sig_jet_phi_scaleup[iyear]}, 
          reduce_maxfirst).Name("sys_lead_jet_phi_scaleup"+year));
      sys_lead_jet_phi_scaledn.push_back(MultiReduceNamedFuncCached(
          {&sys_sig_jet_pt_scaledn[iyear],&sys_sig_jet_phi_scaledn[iyear]}, 
          reduce_maxfirst).Name("sys_lead_jet_phi_scaledn"+year));
      sys_lead_jet_phi_resup.push_back(MultiReduceNamedFuncCached(
          {&sys_sig_jet_pt_resup[iyear],&sys_sig_jet_phi_resup[iyear]}, 
          reduce_maxfirst).Name("sys_lead_jet_phi_resup"+year));
      sys_lead_jet_phi_resdn.push_back(MultiReduceNamedFuncCached(
          {&sys_sig_jet_pt_resdn[iyear],&sys_sig_jet_phi_resdn[iyear]}, 
          reduce_maxfirst).Name("sys_lead_jet_phi_resdn"+year));
      sys_lead_jet_m_scaleup.push_back(MultiReduceNamedFuncCached(
          {&sys_sig_jet_pt_scaleup[iyear],&sys_sig_jet_m_scaleup[iyear]}, 
          reduce_maxfirst).Name("sys_lead_jet_m_scaleup"+year));
      sys_lead_jet_m_scaledn.push_back(MultiReduceNamedFuncCached(
          {&sys_sig_jet_pt_scaledn[iyear],&sys_sig_jet_m_scaledn[iyear]}, 
          reduce_maxfirst).Name("sys_lead_jet_m_scaledn"+year));
      sys_lead_jet_m_resup.push_back(MultiReduceNamedFuncCached(
          {&sys_sig_jet_pt_resup[iyear],&sys_sig_jet_m_resup[iyear]}, 
          reduce_maxfirst).Name("sys_lead_jet_m_resup"+year));
      sys_lead_jet_m_resdn.push_back(MultiReduceNamedFuncCached(
          {&sys_sig_jet_pt_resdn[iyear],&sys_sig_jet_m_resdn[iyear]}, 
          reduce_maxfirst).Name("sys_lead_jet_m_resdn"+year));
      sys_sublead_jet_eta_scaleup.push_back(MultiReduceNamedFuncCached(
          {&sys_sig_jet_pt_scaleup[iyear],&sys_sig_jet_eta_scaleup[iyear]}, 
          reduce_subleadfirst).Name("sys_sublead_jet_eta_scaleup"+year));
      sys_sublead_jet_eta_scaledn.push_back(MultiReduceNamedFuncCached(
          {&sys_sig_jet_pt_scaledn[iyear],&sys_sig_jet_eta_scaledn[iyear]}, 
          reduce_subleadfirst).Name("sys_sublead_jet_eta_scaledn"+year));
      sys_sublead_jet_eta_resup.push_back(MultiReduceNamedFuncCached(
          {&sys_sig_jet_pt_resup[iyear],&sys_sig_jet_eta_resup[iyear]}, 
          reduce_subleadfirst).Name("sys_sublead_jet_eta_resup"+year));
      sys_sublead_jet_eta_resdn.push_back(MultiReduceNamedFuncCached(
          {&sys_sig_jet_pt_resdn[iyear],&sys_sig_jet_eta_resdn[iyear]}, 
          reduce_subleadfirst).Name("sys_sublead_jet_eta_resdn"+year));
      sys_sublead_jet_phi_scaleup.push_back(MultiReduceNamedFuncCached(
          {&sys_sig_jet_pt_scaleup[iyear],&sys_sig_jet_phi_scaleup[iyear]}, 
          reduce_subleadfirst).Name("sys_sublead_jet_phi_scaleup"+year));
      sys_sublead_jet_phi_scaledn.push_back(MultiReduceNamedFuncCached(
          {&sys_sig_jet_pt_scaledn[iyear],&sys_sig_jet_phi_scaledn[iyear]}, 
          reduce_subleadfirst).Name("sys_sublead_jet_phi_scaledn"+year));
      sys_sublead_jet_phi_resup.push_back(MultiReduceNamedFuncCached(
          {&sys_sig_jet_pt_resup[iyear],&sys_sig_jet_phi_resup[iyear]}, 
          reduce_subleadfirst).Name("sys_sublead_jet_phi_resup"+year));
      sys_sublead_jet_phi_resdn.push_back(MultiReduceNamedFuncCached(
          {&sys_sig_jet_pt_resdn[iyear],&sys_sig_jet_phi_resdn[iyear]}, 
          reduce_subleadfirst).Name("sys_sublead_jet_phi_resdn"+year));
      //sys_njet_scaleup.push_back(assign_sca_variation_year_select("sys_njet[2]",
      //    "njet", year, "sys_njet_scaleup"));
      //sys_njet_scaledn.push_back(assign_sca_variation_year_select("sys_njet[3]",
      //    "njet", year, "sys_njet_scaledn"));
      //sys_njet_resup.push_back(assign_sca_variation_year_select("sys_njet[0]",
      //    "njet", year, "sys_njet_resup"));
      //sys_njet_resdn.push_back(assign_sca_variation_year_select("sys_njet[1]",
      //    "njet", year, "sys_njet_resdn"));
      //sys_nbdfm_scaleup.push_back(assign_sca_variation_year_select("sys_nbm[2]",
      //    "nbdfm", year, "sys_nbdfm_scaleup"));
      //sys_nbdfm_scaledn.push_back(assign_sca_variation_year_select("sys_nbm[3]",
      //    "nbdfm", year, "sys_nbdfm_scaledn"));
      //sys_nbdfm_resup.push_back(assign_sca_variation_year_select("sys_nbm[0]",
      //    "nbdfm", year, "sys_nbdfm_resup"));
      //sys_nbdfm_resdn.push_back(assign_sca_variation_year_select("sys_nbm[1]",
      //    "nbdfm", year, "sys_nbdfm_resdn"));
      //
      sys_njet_scaleup.push_back(ReduceNamedFuncCached(
          sys_jet_isgood_scaleup[iyear],reduce_sum).Name(
          "sys_njet_scaleup"+year));
      sys_njet_scaledn.push_back(ReduceNamedFuncCached(
          sys_jet_isgood_scaledn[iyear],reduce_sum).Name(
          "sys_njet_scaledn"+year));
      sys_njet_resup.push_back(ReduceNamedFuncCached(
          sys_jet_isgood_resup[iyear],reduce_sum).Name(
          "sys_njet_resup"+year));
      sys_njet_resdn.push_back(ReduceNamedFuncCached(
          sys_jet_isgood_resdn[iyear],reduce_sum).Name(
          "sys_njet_resdn"+year));
      sys_nbdfm_scaleup.push_back(assign_variation_nbdfm(
          sys_jet_isgood_scaleup[iyear], "scaleup"+year));
      sys_nbdfm_scaledn.push_back(assign_variation_nbdfm(
          sys_jet_isgood_scaledn[iyear], "scaledn"+year));
      sys_nbdfm_resup.push_back(assign_variation_nbdfm(
          sys_jet_isgood_resup[iyear], "resup"+year));
      sys_nbdfm_resdn.push_back(assign_variation_nbdfm(
          sys_jet_isgood_resdn[iyear], "resdn"+year));
      sys_met_scaleup.push_back(assign_sca_variation_year_select("sys_met[2]",
          "met", year, "sys_met_scaleup"));
      sys_met_scaledn.push_back(assign_sca_variation_year_select("sys_met[3]",
          "met", year, "sys_met_scaledn"));
      sys_met_resup.push_back(assign_sca_variation_year_select("sys_met[0]",
          "met", year, "sys_met_resup"));
      sys_met_resdn.push_back(assign_sca_variation_year_select("sys_met[1]",
          "met", year, "sys_met_resdn"));
      sys_dijet_scaleup.push_back(assign_variation_dijet(
          sys_jet_isgood_scaleup[iyear], sys_jet_pt_scaleup[iyear], 
          sys_jet_m_scaleup[iyear], "scaleup"));
      sys_dijet_scaledn.push_back(assign_variation_dijet(
          sys_jet_isgood_scaledn[iyear], sys_jet_pt_scaledn[iyear], 
          sys_jet_m_scaledn[iyear], "scaledn"));
      sys_dijet_resup.push_back(assign_variation_dijet(
          sys_jet_isgood_resup[iyear], sys_jet_pt_resup[iyear], 
          sys_jet_m_resup[iyear], "resup"));
      sys_dijet_resdn.push_back(assign_variation_dijet(
          sys_jet_isgood_resdn[iyear], sys_jet_pt_resdn[iyear], 
          sys_jet_m_resdn[iyear], "resdn"));
      sys_dijet_phi_scaleup.push_back(ReduceNamedFuncCached(
          sys_dijet_scaleup[iyear], reduce_index2)
          .Name("sys_dijet_phi_scaleup"));
      sys_dijet_phi_scaledn.push_back(ReduceNamedFuncCached(
          sys_dijet_scaledn[iyear], reduce_index2)
          .Name("sys_dijet_phi_scaledn"));
      sys_dijet_phi_resup.push_back(ReduceNamedFuncCached(
          sys_dijet_resup[iyear], reduce_index2).Name("sys_dijet_phi_resup"));
      sys_dijet_phi_resdn.push_back(ReduceNamedFuncCached(
          sys_dijet_resdn[iyear], reduce_index2).Name("sys_dijet_phi_resdn"));
      sys_dijet_m_scaleup.push_back(ReduceNamedFuncCached(
          sys_dijet_scaleup[iyear], reduce_index3)
          .Name("sys_dijet_m_scaleup"));
      sys_dijet_m_scaledn.push_back(ReduceNamedFuncCached(
          sys_dijet_scaledn[iyear], reduce_index3)
          .Name("sys_dijet_m_scaledn"));
      sys_dijet_m_resup.push_back(ReduceNamedFuncCached(
          sys_dijet_resup[iyear], reduce_index3).Name("sys_dijet_m_resup"));
      sys_dijet_m_resdn.push_back(ReduceNamedFuncCached(
          sys_dijet_resdn[iyear], reduce_index3).Name("sys_dijet_m_resdn"));
      sys_dijet_deta_scaleup.push_back(ReduceNamedFuncCached(
          sys_dijet_scaleup[iyear], reduce_index4).Name(
          "sys_dijet_deta_scaleup"));
      sys_dijet_deta_scaledn.push_back(ReduceNamedFuncCached(
          sys_dijet_scaledn[iyear], reduce_index4).Name(
          "sys_dijet_deta_scaledn"));
      sys_dijet_deta_resup.push_back(ReduceNamedFuncCached(
          sys_dijet_resup[iyear], reduce_index4).Name(
          "sys_dijet_deta_resup"));
      sys_dijet_deta_resdn.push_back(ReduceNamedFuncCached(
          sys_dijet_resdn[iyear], reduce_index4).Name(
          "sys_dijet_deta_resdn"));
      sys_dijet_dphi_scaleup.push_back(ReduceNamedFuncCached(
          sys_dijet_scaleup[iyear], reduce_index5).Name(
          "sys_dijet_dphi_scaleup"));
      sys_dijet_dphi_scaledn.push_back(ReduceNamedFuncCached(
          sys_dijet_scaledn[iyear], reduce_index5).Name(
          "sys_dijet_dphi_scaledn"));
      sys_dijet_dphi_resup.push_back(ReduceNamedFuncCached(
          sys_dijet_resup[iyear], reduce_index5).Name(
          "sys_dijet_dphi_resup"));
      sys_dijet_dphi_resdn.push_back(ReduceNamedFuncCached(
          sys_dijet_resdn[iyear], reduce_index5).Name(
          "sys_dijet_dphi_resdn"));
      untagged_category_jetscaleup.push_back(
          assign_variation_untagged_category(
          sys_nlep_default, sys_njet_scaleup[iyear], sys_nbdfm_scaleup[iyear], 
          sys_met_scaleup[iyear], sys_llphoton_pt_default, 
          sys_llphoton_m_default, sys_max_lep_miniso_default, sys_ll_m_default, 
          "jetscaleup"));
      untagged_category_jetscaledn.push_back(
          assign_variation_untagged_category(
          sys_nlep_default, sys_njet_scaledn[iyear], sys_nbdfm_scaledn[iyear], 
          sys_met_scaledn[iyear], sys_llphoton_pt_default, 
          sys_llphoton_m_default, sys_max_lep_miniso_default, sys_ll_m_default, 
          "jetscaledn"));
      untagged_category_jetresup.push_back(assign_variation_untagged_category(
          sys_nlep_default, sys_njet_resup[iyear], sys_nbdfm_resup[iyear], 
          sys_met_resup[iyear], sys_llphoton_pt_default, 
          sys_llphoton_m_default, sys_max_lep_miniso_default, sys_ll_m_default, 
          "jetresup"));
      untagged_category_jetresdn.push_back(assign_variation_untagged_category(
          sys_nlep_default, sys_njet_resdn[iyear], sys_nbdfm_resdn[iyear], 
          sys_met_resdn[iyear], sys_llphoton_pt_default, 
          sys_llphoton_m_default, sys_max_lep_miniso_default, sys_ll_m_default, 
          "jetresdn"));
      sys_llphoton_refit_dijet_dphi_jetscaleup.push_back(
          MultiMapNamedFuncCached({&sys_llphoton_refit_phi_default, 
          &sys_dijet_phi_scaleup[iyear]}, map_deltaphi).Name(
          "sys_llphoton_dijet_dphi_jetscaleup"));
      sys_llphoton_refit_dijet_dphi_jetscaledn.push_back(
          MultiMapNamedFuncCached({&sys_llphoton_refit_phi_default, 
          &sys_dijet_phi_scaledn[iyear]}, map_deltaphi).Name(
          "sys_llphoton_dijet_dphi_jetscaledn"));
      sys_llphoton_refit_dijet_dphi_jetresup.push_back(
          MultiMapNamedFuncCached({&sys_llphoton_refit_phi_default, 
          &sys_dijet_phi_resup[iyear]}, map_deltaphi).Name(
          "sys_llphoton_dijet_dphi_jetresup"));
      sys_llphoton_refit_dijet_dphi_jetresdn.push_back(
          MultiMapNamedFuncCached({&sys_llphoton_refit_phi_default, 
          &sys_dijet_phi_resdn[iyear]}, map_deltaphi).Name(
          "sys_llphoton_dijet_dphi_jetresdn"));
      sys_llphoton_refit_jet_dphi_jetscaleup.push_back(
          assign_variation_llphoton_refit_jet_dphi("llphoton_refit_phi", 
          sys_jet_isgood_scaleup[iyear], sys_jet_pt_scaleup[iyear], 
          "jetscaleup"));
      sys_llphoton_refit_jet_dphi_jetscaledn.push_back(
          assign_variation_llphoton_refit_jet_dphi("llphoton_refit_phi", 
          sys_jet_isgood_scaledn[iyear], sys_jet_pt_scaledn[iyear], 
          "jetscaledn"));
      sys_llphoton_refit_jet_dphi_jetresup.push_back(
          assign_variation_llphoton_refit_jet_dphi("llphoton_refit_phi", 
          sys_jet_isgood_resup[iyear], sys_jet_pt_resup[iyear], "jetresup"));
      sys_llphoton_refit_jet_dphi_jetresdn.push_back(
          assign_variation_llphoton_refit_jet_dphi("llphoton_refit_phi", 
          sys_jet_isgood_resdn[iyear], sys_jet_pt_resdn[iyear], "jetresdn"));
      sys_llphoton_refit_jet_balance_jetscaleup.push_back(
          assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_default, 
          sys_lead_photon_pt_default, sys_lead_photon_phi_default, 
          sys_jet_isgood_scaleup[iyear], 
          sys_jet_pt_scaleup[iyear], "jetscaleup"));
      sys_llphoton_refit_jet_balance_jetscaledn.push_back(
          assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_default, 
          sys_lead_photon_pt_default, sys_lead_photon_phi_default, 
          sys_jet_isgood_scaledn[iyear], 
          sys_jet_pt_scaledn[iyear], "jetscaledn"));
      sys_llphoton_refit_jet_balance_jetresup.push_back(
          assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_default, 
          sys_lead_photon_pt_default, sys_lead_photon_phi_default, 
          sys_jet_isgood_resup[iyear], 
          sys_jet_pt_resup[iyear], "jetresup"));
      sys_llphoton_refit_jet_balance_jetresdn.push_back(
          assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_default, 
          sys_lead_photon_pt_default, sys_lead_photon_phi_default, 
          sys_jet_isgood_resdn[iyear], 
          sys_jet_pt_resdn[iyear], "jetresdn"));
      sys_mht_jetscaleup.push_back(assign_variation_mht_balance(
          sys_el_sig_default, sys_el_pt_default, sys_mu_sig_default, 
          sys_mu_pt_default, sys_photon_sig_default, sys_photon_pt_default, 
          sys_jet_isgood_scaleup[iyear], sys_jet_pt_scaleup[iyear], 
          "jetscaleup"));
      sys_mht_jetscaledn.push_back(assign_variation_mht_balance(
          sys_el_sig_default, sys_el_pt_default, sys_mu_sig_default, 
          sys_mu_pt_default, sys_photon_sig_default, sys_photon_pt_default, 
          sys_jet_isgood_scaledn[iyear], sys_jet_pt_scaledn[iyear], 
          "jetscaledn"));
      sys_mht_jetresup.push_back(assign_variation_mht_balance(
          sys_el_sig_default, sys_el_pt_default, sys_mu_sig_default, 
          sys_mu_pt_default, sys_photon_sig_default, sys_photon_pt_default, 
          sys_jet_isgood_resup[iyear], sys_jet_pt_resup[iyear], "jetresup"));
      sys_mht_jetresdn.push_back(assign_variation_mht_balance(
          sys_el_sig_default, sys_el_pt_default, sys_mu_sig_default, 
          sys_mu_pt_default, sys_photon_sig_default, sys_photon_pt_default, 
          sys_jet_isgood_resdn[iyear], sys_jet_pt_resdn[iyear], "jetresdn"));
      sys_mht_phi_jetscaleup.push_back(ReduceNamedFuncCached(
          sys_mht_jetscaleup[iyear], reduce_index1).Name(
          "sys_mht_phi_jetscaleup"+year));
      sys_mht_phi_jetscaledn.push_back(ReduceNamedFuncCached(
          sys_mht_jetscaledn[iyear], reduce_index1).Name(
          "sys_mht_phi_jetscaledn"+year));
      sys_mht_phi_jetresup.push_back(ReduceNamedFuncCached(
          sys_mht_jetresup[iyear], reduce_index1).Name(
          "sys_mht_phi_jetresup"+year));
      sys_mht_phi_jetresdn.push_back(ReduceNamedFuncCached(
          sys_mht_jetresdn[iyear], reduce_index1).Name(
          "sys_mht_phi_jetresdn"+year));
      sys_photon_mht_dphi_jetscaleup.push_back(MultiMapNamedFuncCached(
          {&sys_lead_photon_phi_default, &sys_mht_phi_jetscaleup[iyear]}, 
           map_deltaphi).Name("sys_photon_mht_dphi_jetscaleup"));
      sys_photon_mht_dphi_jetscaledn.push_back(MultiMapNamedFuncCached(
          {&sys_lead_photon_phi_default, &sys_mht_phi_jetscaledn[iyear]}, 
           map_deltaphi).Name("sys_photon_mht_dphi_jetscaledn"));
      sys_photon_mht_dphi_jetresup.push_back(MultiMapNamedFuncCached(
          {&sys_lead_photon_phi_default, &sys_mht_phi_jetresup[iyear]}, 
           map_deltaphi).Name("sys_photon_mht_dphi_jetresup"));
      sys_photon_mht_dphi_jetresdn.push_back(MultiMapNamedFuncCached(
          {&sys_lead_photon_phi_default, &sys_mht_phi_jetresdn[iyear]}, 
           map_deltaphi).Name("sys_photon_mht_dphi_jetresdn"));
      sys_photon_jet1_dr_jetscaleup.push_back(MultiMapNamedFuncCached(
          {&sys_lead_photon_eta_default, &sys_lead_photon_phi_default, 
           &sys_lead_jet_eta_scaleup[iyear], &sys_lead_jet_phi_scaleup[iyear]}, 
          map_deltar).Name("sys_photon_jet1_dr_jetscaleup"));
      sys_photon_jet1_dr_jetscaledn.push_back(MultiMapNamedFuncCached(
          {&sys_lead_photon_eta_default, &sys_lead_photon_phi_default, 
           &sys_lead_jet_eta_scaledn[iyear], &sys_lead_jet_phi_scaledn[iyear]}, 
          map_deltar).Name("sys_photon_jet1_dr_jetscaledn"));
      sys_photon_jet1_dr_jetresup.push_back(MultiMapNamedFuncCached(
          {&sys_lead_photon_eta_default, &sys_lead_photon_phi_default, 
           &sys_lead_jet_eta_resup[iyear], &sys_lead_jet_phi_resup[iyear]}, 
          map_deltar).Name("sys_photon_jet1_dr_jetresup"));
      sys_photon_jet1_dr_jetresdn.push_back(MultiMapNamedFuncCached(
          {&sys_lead_photon_eta_default, &sys_lead_photon_phi_default, 
           &sys_lead_jet_eta_resdn[iyear], &sys_lead_jet_phi_resdn[iyear]}, 
          map_deltar).Name("sys_photon_jet1_dr_jetresdn"));
      sys_photon_jet2_dr_jetscaleup.push_back(MultiMapNamedFuncCached(
          {&sys_lead_photon_eta_default, &sys_lead_photon_phi_default, 
          &sys_sublead_jet_eta_scaleup[iyear], 
          &sys_sublead_jet_phi_scaleup[iyear]}, map_deltar).Name(
          "sys_photon_jet2_dr_jetscaleup"));
      sys_photon_jet2_dr_jetscaledn.push_back(MultiMapNamedFuncCached(
          {&sys_lead_photon_eta_default, &sys_lead_photon_phi_default, 
          &sys_sublead_jet_eta_scaledn[iyear], 
          &sys_sublead_jet_phi_scaledn[iyear]}, map_deltar).Name(
          "sys_photon_jet2_dr_jetscaledn"));
      sys_photon_jet2_dr_jetresup.push_back(MultiMapNamedFuncCached(
          {&sys_lead_photon_eta_default, &sys_lead_photon_phi_default, 
          &sys_sublead_jet_eta_resup[iyear], 
          &sys_sublead_jet_phi_resup[iyear]}, map_deltar).Name(
          "sys_photon_jet2_dr_jetresup"));
      sys_photon_jet2_dr_jetresdn.push_back(MultiMapNamedFuncCached(
          {&sys_lead_photon_eta_default, &sys_lead_photon_phi_default, 
          &sys_sublead_jet_eta_resdn[iyear], 
          &sys_sublead_jet_phi_resdn[iyear]}, map_deltar).Name(
          "sys_photon_jet2_dr_jetresdn"));
      sys_photon_zeppenfeld_jetscaleup.push_back(
          assign_variation_photon_zeppenfeld(sys_lead_photon_eta_default,
          sys_njet_scaleup[iyear],sys_lead_jet_eta_scaleup[iyear], 
          sys_sublead_jet_eta_scaleup[iyear],"jetscaleup"));
      sys_photon_zeppenfeld_jetscaledn.push_back(
          assign_variation_photon_zeppenfeld(sys_lead_photon_eta_default,
          sys_njet_scaledn[iyear], sys_lead_jet_eta_scaledn[iyear], 
          sys_sublead_jet_eta_scaledn[iyear], "jetscaledn"));
      sys_photon_zeppenfeld_jetresup.push_back(
          assign_variation_photon_zeppenfeld(sys_lead_photon_eta_default,
          sys_njet_resup[iyear], sys_lead_jet_eta_resup[iyear], 
          sys_sublead_jet_eta_resup[iyear], "jetresup"));
      sys_photon_zeppenfeld_jetresdn.push_back(
          assign_variation_photon_zeppenfeld(sys_lead_photon_eta_default,
          sys_njet_resdn[iyear], sys_lead_jet_eta_resdn[iyear], 
          sys_sublead_jet_eta_resdn[iyear], "jetresdn"));
      sys_bchig_uncorr_up.push_back(assign_sca_variation_year_select(
          "sys_bchig_uncorr[0]/w_bhig_df","1",year,"sys_bchig_uncorr_up"));
      sys_bchig_uncorr_dn.push_back(assign_sca_variation_year_select(
          "sys_bchig_uncorr[1]/w_bhig_df","1",year,"sys_bchig_uncorr_dn"));
      sys_udsghig_uncorr_up.push_back(assign_sca_variation_year_select(
          "sys_udsghig_uncorr[0]/w_bhig_df","1",year,"sys_udsghig_uncorr_up"));
      sys_udsghig_uncorr_dn.push_back(assign_sca_variation_year_select(
          "sys_udsghig_uncorr[1]/w_bhig_df","1",year,"sys_udsghig_uncorr_dn"));
    }
  }

  //old (25-03) BDTs
  //returns working version of dijet BDT
  vector<shared_ptr<MVAWrapper>> SystVbfBdts() {
    vector<shared_ptr<MVAWrapper>> vbf_bdt_readers = {
      std::make_shared<MVAWrapper>("vbf_bdt0"),
      std::make_shared<MVAWrapper>("vbf_bdt1"),
      std::make_shared<MVAWrapper>("vbf_bdt2"),
      std::make_shared<MVAWrapper>("vbf_bdt3")};
    for (std::shared_ptr<MVAWrapper> &bdt_reader : vbf_bdt_readers) {
      //default
      bdt_reader->SetVariableRef("dijet_deta",sys_dijet_deta_default);
      bdt_reader->SetVariableRef("dijet_dphi",sys_dijet_dphi_default);
      bdt_reader->SetVariableRef("j1_pt",sys_lead_jet_pt_default);
      bdt_reader->SetVariableRef("j2_pt",sys_sublead_jet_pt_default);
      bdt_reader->SetVariableRef("llphoton_dijet_balance",
                                 sys_llphoton_refit_jet_balance_default);
      bdt_reader->SetVariableRef("llphoton_dijet_dphi",
                                 sys_llphoton_refit_dijet_dphi_default);
      bdt_reader->SetVariableRef("phi",sys_llphoton_refit_psi_default);
      bdt_reader->SetVariableRef("min_dR",sys_lead_photon_drmin_default);
      bdt_reader->SetVariableRef("max_dR",sys_lead_photon_drmax_default);
      bdt_reader->SetVariableRef("costheta",
                                 sys_llphoton_refit_costheta_default);
      bdt_reader->SetVariableRef("cosTheta",
                                 sys_llphoton_refit_cosTheta_default);
      bdt_reader->SetVariableRef("pt_mass",sys_llphoton_refit_relpt_default);
      bdt_reader->SetVariableRef("l1_rapidity",sys_lead_lepton_eta_default);
      bdt_reader->SetVariableRef("l2_rapidity",sys_sublead_lepton_eta_default);
      bdt_reader->SetVariableRef("photon_rapidity",
                                 sys_lead_photon_eta_default);
      bdt_reader->SetVariableRef("photon_mva",sys_lead_photon_idmva_default);
      bdt_reader->SetVariableRef("photon_res",
                                 sys_lead_photon_relpterr_default);
      bdt_reader->SetVariableRef("photon_zeppenfeld",
                                 sys_photon_zeppenfeld_default);
      bdt_reader->SetVariableRef("photon_jet1_dr",sys_photon_jet1_dr_default);
      bdt_reader->SetVariableRef("photon_jet2_dr",sys_photon_jet2_dr_default);
      //elscaleup
      bdt_reader->SetAltVariableRef("elscaleup",sys_dijet_deta_default);
      bdt_reader->SetAltVariableRef("elscaleup",sys_dijet_dphi_default);
      bdt_reader->SetAltVariableRef("elscaleup",sys_lead_jet_pt_default);
      bdt_reader->SetAltVariableRef("elscaleup",sys_sublead_jet_pt_default);
      bdt_reader->SetAltVariableRef("elscaleup",
                                    sys_llphoton_refit_jet_balance_elscaleup);
      bdt_reader->SetAltVariableRef("elscaleup",
                                    sys_llphoton_refit_dijet_dphi_elscaleup);
      bdt_reader->SetAltVariableRef("elscaleup",
                                    sys_llphoton_refit_psi_elscaleup);
      bdt_reader->SetAltVariableRef("elscaleup",sys_lead_photon_drmin_default);
      bdt_reader->SetAltVariableRef("elscaleup",sys_lead_photon_drmax_default);
      bdt_reader->SetAltVariableRef("elscaleup",
                                    sys_llphoton_refit_costheta_elscaleup);
      bdt_reader->SetAltVariableRef("elscaleup",
                                    sys_llphoton_refit_cosTheta_elscaleup);
      bdt_reader->SetAltVariableRef("elscaleup",
                                    sys_llphoton_refit_relpt_elscaleup);
      bdt_reader->SetAltVariableRef("elscaleup",sys_lead_lepton_eta_elscaleup);
      bdt_reader->SetAltVariableRef("elscaleup",
                                    sys_sublead_lepton_eta_elscaleup);
      bdt_reader->SetAltVariableRef("elscaleup",sys_lead_photon_eta_default);
      bdt_reader->SetAltVariableRef("elscaleup",sys_lead_photon_idmva_default);
      bdt_reader->SetAltVariableRef("elscaleup",
                                    sys_lead_photon_relpterr_default);
      bdt_reader->SetAltVariableRef("elscaleup",
                                    sys_photon_zeppenfeld_default);
      bdt_reader->SetAltVariableRef("elscaleup",sys_photon_jet1_dr_default);
      bdt_reader->SetAltVariableRef("elscaleup",sys_photon_jet2_dr_default);
      //elscaledn
      bdt_reader->SetAltVariableRef("elscaledn",sys_dijet_deta_default);
      bdt_reader->SetAltVariableRef("elscaledn",sys_dijet_dphi_default);
      bdt_reader->SetAltVariableRef("elscaledn",sys_lead_jet_pt_default);
      bdt_reader->SetAltVariableRef("elscaledn",sys_sublead_jet_pt_default);
      bdt_reader->SetAltVariableRef("elscaledn",
                                    sys_llphoton_refit_jet_balance_elscaledn);
      bdt_reader->SetAltVariableRef("elscaledn",
                                    sys_llphoton_refit_dijet_dphi_elscaledn);
      bdt_reader->SetAltVariableRef("elscaledn",
                                    sys_llphoton_refit_psi_elscaledn);
      bdt_reader->SetAltVariableRef("elscaledn",sys_lead_photon_drmin_default);
      bdt_reader->SetAltVariableRef("elscaledn",sys_lead_photon_drmax_default);
      bdt_reader->SetAltVariableRef("elscaledn",
                                    sys_llphoton_refit_costheta_elscaledn);
      bdt_reader->SetAltVariableRef("elscaledn",
                                    sys_llphoton_refit_cosTheta_elscaledn);
      bdt_reader->SetAltVariableRef("elscaledn",
                                    sys_llphoton_refit_relpt_elscaledn);
      bdt_reader->SetAltVariableRef("elscaledn",sys_lead_lepton_eta_elscaledn);
      bdt_reader->SetAltVariableRef("elscaledn",
                                    sys_sublead_lepton_eta_elscaledn);
      bdt_reader->SetAltVariableRef("elscaledn",sys_lead_photon_eta_default);
      bdt_reader->SetAltVariableRef("elscaledn",sys_lead_photon_idmva_default);
      bdt_reader->SetAltVariableRef("elscaledn",
                                    sys_lead_photon_relpterr_default);
      bdt_reader->SetAltVariableRef("elscaledn",
                                    sys_photon_zeppenfeld_default);
      bdt_reader->SetAltVariableRef("elscaledn",sys_photon_jet1_dr_default);
      bdt_reader->SetAltVariableRef("elscaledn",sys_photon_jet2_dr_default);
      //elresup
      bdt_reader->SetAltVariableRef("elresup",sys_dijet_deta_default);
      bdt_reader->SetAltVariableRef("elresup",sys_dijet_dphi_default);
      bdt_reader->SetAltVariableRef("elresup",sys_lead_jet_pt_default);
      bdt_reader->SetAltVariableRef("elresup",sys_sublead_jet_pt_default);
      bdt_reader->SetAltVariableRef("elresup",
                                    sys_llphoton_refit_jet_balance_elresup);
      bdt_reader->SetAltVariableRef("elresup",
                                    sys_llphoton_refit_dijet_dphi_elresup);
      bdt_reader->SetAltVariableRef("elresup",sys_llphoton_refit_psi_elresup);
      bdt_reader->SetAltVariableRef("elresup",sys_lead_photon_drmin_default);
      bdt_reader->SetAltVariableRef("elresup",sys_lead_photon_drmax_default);
      bdt_reader->SetAltVariableRef("elresup",
                                    sys_llphoton_refit_costheta_elresup);
      bdt_reader->SetAltVariableRef("elresup",
                                    sys_llphoton_refit_cosTheta_elresup);
      bdt_reader->SetAltVariableRef("elresup",
                                    sys_llphoton_refit_relpt_elresup);
      bdt_reader->SetAltVariableRef("elresup",sys_lead_lepton_eta_elresup);
      bdt_reader->SetAltVariableRef("elresup",sys_sublead_lepton_eta_elresup);
      bdt_reader->SetAltVariableRef("elresup",sys_lead_photon_eta_default);
      bdt_reader->SetAltVariableRef("elresup",sys_lead_photon_idmva_default);
      bdt_reader->SetAltVariableRef("elresup",
                                    sys_lead_photon_relpterr_default);
      bdt_reader->SetAltVariableRef("elresup",
                                    sys_photon_zeppenfeld_default);
      bdt_reader->SetAltVariableRef("elresup",sys_photon_jet1_dr_default);
      bdt_reader->SetAltVariableRef("elresup",sys_photon_jet2_dr_default);
      //elresdn
      bdt_reader->SetAltVariableRef("elresdn",sys_dijet_deta_default);
      bdt_reader->SetAltVariableRef("elresdn",sys_dijet_dphi_default);
      bdt_reader->SetAltVariableRef("elresdn",sys_lead_jet_pt_default);
      bdt_reader->SetAltVariableRef("elresdn",sys_sublead_jet_pt_default);
      bdt_reader->SetAltVariableRef("elresdn",
                                    sys_llphoton_refit_jet_balance_elresdn);
      bdt_reader->SetAltVariableRef("elresdn",
                                    sys_llphoton_refit_dijet_dphi_elresdn);
      bdt_reader->SetAltVariableRef("elresdn",sys_llphoton_refit_psi_elresdn);
      bdt_reader->SetAltVariableRef("elresdn",sys_lead_photon_drmin_default);
      bdt_reader->SetAltVariableRef("elresdn",sys_lead_photon_drmax_default);
      bdt_reader->SetAltVariableRef("elresdn",
                                    sys_llphoton_refit_costheta_elresdn);
      bdt_reader->SetAltVariableRef("elresdn",
                                    sys_llphoton_refit_cosTheta_elresdn);
      bdt_reader->SetAltVariableRef("elresdn",
                                    sys_llphoton_refit_relpt_elresdn);
      bdt_reader->SetAltVariableRef("elresdn",sys_lead_lepton_eta_elresdn);
      bdt_reader->SetAltVariableRef("elresdn",sys_sublead_lepton_eta_elresdn);
      bdt_reader->SetAltVariableRef("elresdn",sys_lead_photon_eta_default);
      bdt_reader->SetAltVariableRef("elresdn",sys_lead_photon_idmva_default);
      bdt_reader->SetAltVariableRef("elresdn",
                                    sys_lead_photon_relpterr_default);
      bdt_reader->SetAltVariableRef("elresdn",
                                    sys_photon_zeppenfeld_default);
      bdt_reader->SetAltVariableRef("elresdn",sys_photon_jet1_dr_default);
      bdt_reader->SetAltVariableRef("elresdn",sys_photon_jet2_dr_default);
      //muscaleup
      bdt_reader->SetAltVariableRef("muscaleup",sys_dijet_deta_default);
      bdt_reader->SetAltVariableRef("muscaleup",sys_dijet_dphi_default);
      bdt_reader->SetAltVariableRef("muscaleup",sys_lead_jet_pt_default);
      bdt_reader->SetAltVariableRef("muscaleup",sys_sublead_jet_pt_default);
      bdt_reader->SetAltVariableRef("muscaleup",
                                    sys_llphoton_refit_jet_balance_muscaleup);
      bdt_reader->SetAltVariableRef("muscaleup",
                                    sys_llphoton_refit_dijet_dphi_muscaleup);
      bdt_reader->SetAltVariableRef("muscaleup",
                                    sys_llphoton_refit_psi_muscaleup);
      bdt_reader->SetAltVariableRef("muscaleup",sys_lead_photon_drmin_default);
      bdt_reader->SetAltVariableRef("muscaleup",sys_lead_photon_drmax_default);
      bdt_reader->SetAltVariableRef("muscaleup",
                                    sys_llphoton_refit_costheta_muscaleup);
      bdt_reader->SetAltVariableRef("muscaleup",
                                    sys_llphoton_refit_cosTheta_muscaleup);
      bdt_reader->SetAltVariableRef("muscaleup",
                                    sys_llphoton_refit_relpt_muscaleup);
      bdt_reader->SetAltVariableRef("muscaleup",sys_lead_lepton_eta_muscaleup);
      bdt_reader->SetAltVariableRef("muscaleup",
                                    sys_sublead_lepton_eta_muscaleup);
      bdt_reader->SetAltVariableRef("muscaleup",sys_lead_photon_eta_default);
      bdt_reader->SetAltVariableRef("muscaleup",sys_lead_photon_idmva_default);
      bdt_reader->SetAltVariableRef("muscaleup",
                                    sys_lead_photon_relpterr_default);
      bdt_reader->SetAltVariableRef("muscaleup",
                                    sys_photon_zeppenfeld_default);
      bdt_reader->SetAltVariableRef("muscaleup",sys_photon_jet1_dr_default);
      bdt_reader->SetAltVariableRef("muscaleup",sys_photon_jet2_dr_default);
      //muscaledn
      bdt_reader->SetAltVariableRef("muscaledn",sys_dijet_deta_default);
      bdt_reader->SetAltVariableRef("muscaledn",sys_dijet_dphi_default);
      bdt_reader->SetAltVariableRef("muscaledn",sys_lead_jet_pt_default);
      bdt_reader->SetAltVariableRef("muscaledn",sys_sublead_jet_pt_default);
      bdt_reader->SetAltVariableRef("muscaledn",
                                    sys_llphoton_refit_jet_balance_muscaledn);
      bdt_reader->SetAltVariableRef("muscaledn",
                                    sys_llphoton_refit_dijet_dphi_muscaledn);
      bdt_reader->SetAltVariableRef("muscaledn",sys_llphoton_refit_psi_muscaledn);
      bdt_reader->SetAltVariableRef("muscaledn",sys_lead_photon_drmin_default);
      bdt_reader->SetAltVariableRef("muscaledn",sys_lead_photon_drmax_default);
      bdt_reader->SetAltVariableRef("muscaledn",
                                    sys_llphoton_refit_costheta_muscaledn);
      bdt_reader->SetAltVariableRef("muscaledn",
                                    sys_llphoton_refit_cosTheta_muscaledn);
      bdt_reader->SetAltVariableRef("muscaledn",
                                    sys_llphoton_refit_relpt_muscaledn);
      bdt_reader->SetAltVariableRef("muscaledn",sys_lead_lepton_eta_muscaledn);
      bdt_reader->SetAltVariableRef("muscaledn",sys_sublead_lepton_eta_muscaledn);
      bdt_reader->SetAltVariableRef("muscaledn",sys_lead_photon_eta_default);
      bdt_reader->SetAltVariableRef("muscaledn",sys_lead_photon_idmva_default);
      bdt_reader->SetAltVariableRef("muscaledn",
                                    sys_lead_photon_relpterr_default);
      bdt_reader->SetAltVariableRef("muscaledn",
                                    sys_photon_zeppenfeld_default);
      bdt_reader->SetAltVariableRef("muscaledn",sys_photon_jet1_dr_default);
      bdt_reader->SetAltVariableRef("muscaledn",sys_photon_jet2_dr_default);
      //muresup
      bdt_reader->SetAltVariableRef("muresup",sys_dijet_deta_default);
      bdt_reader->SetAltVariableRef("muresup",sys_dijet_dphi_default);
      bdt_reader->SetAltVariableRef("muresup",sys_lead_jet_pt_default);
      bdt_reader->SetAltVariableRef("muresup",sys_sublead_jet_pt_default);
      bdt_reader->SetAltVariableRef("muresup",
                                    sys_llphoton_refit_jet_balance_muresup);
      bdt_reader->SetAltVariableRef("muresup",
                                    sys_llphoton_refit_dijet_dphi_muresup);
      bdt_reader->SetAltVariableRef("muresup",sys_llphoton_refit_psi_muresup);
      bdt_reader->SetAltVariableRef("muresup",sys_lead_photon_drmin_default);
      bdt_reader->SetAltVariableRef("muresup",sys_lead_photon_drmax_default);
      bdt_reader->SetAltVariableRef("muresup",
                                    sys_llphoton_refit_costheta_muresup);
      bdt_reader->SetAltVariableRef("muresup",
                                    sys_llphoton_refit_cosTheta_muresup);
      bdt_reader->SetAltVariableRef("muresup",
                                    sys_llphoton_refit_relpt_muresup);
      bdt_reader->SetAltVariableRef("muresup",sys_lead_lepton_eta_muresup);
      bdt_reader->SetAltVariableRef("muresup",sys_sublead_lepton_eta_muresup);
      bdt_reader->SetAltVariableRef("muresup",sys_lead_photon_eta_default);
      bdt_reader->SetAltVariableRef("muresup",sys_lead_photon_idmva_default);
      bdt_reader->SetAltVariableRef("muresup",
                                    sys_lead_photon_relpterr_default);
      bdt_reader->SetAltVariableRef("muresup",
                                    sys_photon_zeppenfeld_default);
      bdt_reader->SetAltVariableRef("muresup",sys_photon_jet1_dr_default);
      bdt_reader->SetAltVariableRef("muresup",sys_photon_jet2_dr_default);
      //muresdn
      bdt_reader->SetAltVariableRef("muresdn",sys_dijet_deta_default);
      bdt_reader->SetAltVariableRef("muresdn",sys_dijet_dphi_default);
      bdt_reader->SetAltVariableRef("muresdn",sys_lead_jet_pt_default);
      bdt_reader->SetAltVariableRef("muresdn",sys_sublead_jet_pt_default);
      bdt_reader->SetAltVariableRef("muresdn",
                                    sys_llphoton_refit_jet_balance_muresdn);
      bdt_reader->SetAltVariableRef("muresdn",
                                    sys_llphoton_refit_dijet_dphi_muresdn);
      bdt_reader->SetAltVariableRef("muresdn",sys_llphoton_refit_psi_muresdn);
      bdt_reader->SetAltVariableRef("muresdn",sys_lead_photon_drmin_default);
      bdt_reader->SetAltVariableRef("muresdn",sys_lead_photon_drmax_default);
      bdt_reader->SetAltVariableRef("muresdn",
                                    sys_llphoton_refit_costheta_muresdn);
      bdt_reader->SetAltVariableRef("muresdn",
                                    sys_llphoton_refit_cosTheta_muresdn);
      bdt_reader->SetAltVariableRef("muresdn",
                                    sys_llphoton_refit_relpt_muresdn);
      bdt_reader->SetAltVariableRef("muresdn",sys_lead_lepton_eta_muresdn);
      bdt_reader->SetAltVariableRef("muresdn",sys_sublead_lepton_eta_muresdn);
      bdt_reader->SetAltVariableRef("muresdn",sys_lead_photon_eta_default);
      bdt_reader->SetAltVariableRef("muresdn",sys_lead_photon_idmva_default);
      bdt_reader->SetAltVariableRef("muresdn",
                                    sys_lead_photon_relpterr_default);
      bdt_reader->SetAltVariableRef("muresdn",
                                    sys_photon_zeppenfeld_default);
      bdt_reader->SetAltVariableRef("muresdn",sys_photon_jet1_dr_default);
      bdt_reader->SetAltVariableRef("muresdn",sys_photon_jet2_dr_default);
      //phscaleup
      bdt_reader->SetAltVariableRef("phscaleup",sys_dijet_deta_default);
      bdt_reader->SetAltVariableRef("phscaleup",sys_dijet_dphi_default);
      bdt_reader->SetAltVariableRef("phscaleup",sys_lead_jet_pt_default);
      bdt_reader->SetAltVariableRef("phscaleup",sys_sublead_jet_pt_default);
      bdt_reader->SetAltVariableRef("phscaleup",
                                    sys_llphoton_refit_jet_balance_phscaleup);
      bdt_reader->SetAltVariableRef("phscaleup",
                                    sys_llphoton_refit_dijet_dphi_phscaleup);
      bdt_reader->SetAltVariableRef("phscaleup",
                                    sys_llphoton_refit_psi_phscaleup);
      bdt_reader->SetAltVariableRef("phscaleup",sys_lead_photon_drmin_scaleup);
      bdt_reader->SetAltVariableRef("phscaleup",sys_lead_photon_drmax_scaleup);
      bdt_reader->SetAltVariableRef("phscaleup",
                                    sys_llphoton_refit_costheta_phscaleup);
      bdt_reader->SetAltVariableRef("phscaleup",
                                    sys_llphoton_refit_cosTheta_phscaleup);
      bdt_reader->SetAltVariableRef("phscaleup",
                                    sys_llphoton_refit_relpt_phscaleup);
      bdt_reader->SetAltVariableRef("phscaleup",sys_lead_lepton_eta_default);
      bdt_reader->SetAltVariableRef("phscaleup",sys_sublead_lepton_eta_default);
      bdt_reader->SetAltVariableRef("phscaleup",sys_lead_photon_eta_scaleup);
      bdt_reader->SetAltVariableRef("phscaleup",sys_lead_photon_idmva_scaleup);
      bdt_reader->SetAltVariableRef("phscaleup",
                                    sys_lead_photon_relpterr_scaleup);
      bdt_reader->SetAltVariableRef("phscaleup",
                                    sys_photon_zeppenfeld_phscaleup);
      bdt_reader->SetAltVariableRef("phscaleup",sys_photon_jet1_dr_phscaleup);
      bdt_reader->SetAltVariableRef("phscaleup",sys_photon_jet2_dr_phscaleup);
      //phscaledn
      bdt_reader->SetAltVariableRef("phscaledn",sys_dijet_deta_default);
      bdt_reader->SetAltVariableRef("phscaledn",sys_dijet_dphi_default);
      bdt_reader->SetAltVariableRef("phscaledn",sys_lead_jet_pt_default);
      bdt_reader->SetAltVariableRef("phscaledn",sys_sublead_jet_pt_default);
      bdt_reader->SetAltVariableRef("phscaledn",
                                    sys_llphoton_refit_jet_balance_phscaledn);
      bdt_reader->SetAltVariableRef("phscaledn",
                                    sys_llphoton_refit_dijet_dphi_phscaledn);
      bdt_reader->SetAltVariableRef("phscaledn",
                                    sys_llphoton_refit_psi_phscaledn);
      bdt_reader->SetAltVariableRef("phscaledn",sys_lead_photon_drmin_scaledn);
      bdt_reader->SetAltVariableRef("phscaledn",sys_lead_photon_drmax_scaledn);
      bdt_reader->SetAltVariableRef("phscaledn",
                                    sys_llphoton_refit_costheta_phscaledn);
      bdt_reader->SetAltVariableRef("phscaledn",
                                    sys_llphoton_refit_cosTheta_phscaledn);
      bdt_reader->SetAltVariableRef("phscaledn",
                                    sys_llphoton_refit_relpt_phscaledn);
      bdt_reader->SetAltVariableRef("phscaledn",sys_lead_lepton_eta_default);
      bdt_reader->SetAltVariableRef("phscaledn",sys_sublead_lepton_eta_default);
      bdt_reader->SetAltVariableRef("phscaledn",sys_lead_photon_eta_scaledn);
      bdt_reader->SetAltVariableRef("phscaledn",sys_lead_photon_idmva_scaledn);
      bdt_reader->SetAltVariableRef("phscaledn",
                                    sys_lead_photon_relpterr_scaledn);
      bdt_reader->SetAltVariableRef("phscaledn",
                                    sys_photon_zeppenfeld_phscaledn);
      bdt_reader->SetAltVariableRef("phscaledn",sys_photon_jet1_dr_phscaledn);
      bdt_reader->SetAltVariableRef("phscaledn",sys_photon_jet2_dr_phscaledn);
      //phresup
      bdt_reader->SetAltVariableRef("phresup",sys_dijet_deta_default);
      bdt_reader->SetAltVariableRef("phresup",sys_dijet_dphi_default);
      bdt_reader->SetAltVariableRef("phresup",sys_lead_jet_pt_default);
      bdt_reader->SetAltVariableRef("phresup",sys_sublead_jet_pt_default);
      bdt_reader->SetAltVariableRef("phresup",
                                    sys_llphoton_refit_jet_balance_phresup);
      bdt_reader->SetAltVariableRef("phresup",
                                    sys_llphoton_refit_dijet_dphi_phresup);
      bdt_reader->SetAltVariableRef("phresup",sys_llphoton_refit_psi_phresup);
      bdt_reader->SetAltVariableRef("phresup",sys_lead_photon_drmin_resup);
      bdt_reader->SetAltVariableRef("phresup",sys_lead_photon_drmax_resup);
      bdt_reader->SetAltVariableRef("phresup",
                                    sys_llphoton_refit_costheta_phresup);
      bdt_reader->SetAltVariableRef("phresup",
                                    sys_llphoton_refit_cosTheta_phresup);
      bdt_reader->SetAltVariableRef("phresup",
                                    sys_llphoton_refit_relpt_phresup);
      bdt_reader->SetAltVariableRef("phresup",sys_lead_lepton_eta_default);
      bdt_reader->SetAltVariableRef("phresup",sys_sublead_lepton_eta_default);
      bdt_reader->SetAltVariableRef("phresup",sys_lead_photon_eta_resup);
      bdt_reader->SetAltVariableRef("phresup",sys_lead_photon_idmva_resup);
      bdt_reader->SetAltVariableRef("phresup",
                                    sys_lead_photon_relpterr_resup);
      bdt_reader->SetAltVariableRef("phresup",
                                    sys_photon_zeppenfeld_phresup);
      bdt_reader->SetAltVariableRef("phresup",sys_photon_jet1_dr_phresup);
      bdt_reader->SetAltVariableRef("phresup",sys_photon_jet2_dr_phresup);
      //phresdn
      bdt_reader->SetAltVariableRef("phresdn",sys_dijet_deta_default);
      bdt_reader->SetAltVariableRef("phresdn",sys_dijet_dphi_default);
      bdt_reader->SetAltVariableRef("phresdn",sys_lead_jet_pt_default);
      bdt_reader->SetAltVariableRef("phresdn",sys_sublead_jet_pt_default);
      bdt_reader->SetAltVariableRef("phresdn",
                                    sys_llphoton_refit_jet_balance_phresdn);
      bdt_reader->SetAltVariableRef("phresdn",
                                    sys_llphoton_refit_dijet_dphi_phresdn);
      bdt_reader->SetAltVariableRef("phresdn",sys_llphoton_refit_psi_phresdn);
      bdt_reader->SetAltVariableRef("phresdn",sys_lead_photon_drmin_resdn);
      bdt_reader->SetAltVariableRef("phresdn",sys_lead_photon_drmax_resdn);
      bdt_reader->SetAltVariableRef("phresdn",
                                    sys_llphoton_refit_costheta_phresdn);
      bdt_reader->SetAltVariableRef("phresdn",
                                    sys_llphoton_refit_cosTheta_phresdn);
      bdt_reader->SetAltVariableRef("phresdn",
                                    sys_llphoton_refit_relpt_phresdn);
      bdt_reader->SetAltVariableRef("phresdn",sys_lead_lepton_eta_default);
      bdt_reader->SetAltVariableRef("phresdn",sys_sublead_lepton_eta_default);
      bdt_reader->SetAltVariableRef("phresdn",sys_lead_photon_eta_resdn);
      bdt_reader->SetAltVariableRef("phresdn",sys_lead_photon_idmva_resdn);
      bdt_reader->SetAltVariableRef("phresdn",
                                    sys_lead_photon_relpterr_resdn);
      bdt_reader->SetAltVariableRef("phresdn",
                                    sys_photon_zeppenfeld_phresdn);
      bdt_reader->SetAltVariableRef("phresdn",sys_photon_jet1_dr_phresdn);
      bdt_reader->SetAltVariableRef("phresdn",sys_photon_jet2_dr_phresdn);
      const vector<string> years = {"2016APV", "2016", "2017", "2018", 
                                    "2022", "2022EE", "2023", "2023BPix"};
      for (unsigned iyear = 0; iyear < years.size(); iyear++) {
        string year = years[iyear];
        //jetscaleup
        bdt_reader->SetAltVariableRef("jetscaleup"+year,
                                      sys_dijet_deta_scaleup[iyear]);
        bdt_reader->SetAltVariableRef("jetscaleup"+year,
                                      sys_dijet_dphi_scaleup[iyear]);
        bdt_reader->SetAltVariableRef("jetscaleup"+year,
                                      sys_lead_jet_pt_scaleup[iyear]);
        bdt_reader->SetAltVariableRef("jetscaleup"+year,
                                      sys_sublead_jet_pt_scaleup[iyear]);
        bdt_reader->SetAltVariableRef("jetscaleup"+year,
            sys_llphoton_refit_jet_balance_jetscaleup[iyear]);
        bdt_reader->SetAltVariableRef("jetscaleup"+year,
            sys_llphoton_refit_dijet_dphi_jetscaleup[iyear]);
        bdt_reader->SetAltVariableRef("jetscaleup"+year,
                                      sys_llphoton_refit_psi_default);
        bdt_reader->SetAltVariableRef("jetscaleup"+year,
                                      sys_lead_photon_drmin_default);
        bdt_reader->SetAltVariableRef("jetscaleup"+year,
                                      sys_lead_photon_drmax_default);
        bdt_reader->SetAltVariableRef("jetscaleup"+year,
                                      sys_llphoton_refit_costheta_default);
        bdt_reader->SetAltVariableRef("jetscaleup"+year,
                                      sys_llphoton_refit_cosTheta_default);
        bdt_reader->SetAltVariableRef("jetscaleup"+year,
                                      sys_llphoton_refit_relpt_default);
        bdt_reader->SetAltVariableRef("jetscaleup"+year,
                                      sys_lead_lepton_eta_default);
        bdt_reader->SetAltVariableRef("jetscaleup"+year,
                                      sys_sublead_lepton_eta_default);
        bdt_reader->SetAltVariableRef("jetscaleup"+year,
            sys_lead_photon_eta_default);
        bdt_reader->SetAltVariableRef("jetscaleup"+year,
            sys_lead_photon_idmva_default);
        bdt_reader->SetAltVariableRef("jetscaleup"+year,
                                      sys_lead_photon_relpterr_default);
        bdt_reader->SetAltVariableRef("jetscaleup"+year,
                                      sys_photon_zeppenfeld_jetscaleup[iyear]);
        bdt_reader->SetAltVariableRef("jetscaleup"+year,
                                      sys_photon_jet1_dr_jetscaleup[iyear]);
        bdt_reader->SetAltVariableRef("jetscaleup"+year,
                                      sys_photon_jet2_dr_jetscaleup[iyear]);
        //jetscaledn
        bdt_reader->SetAltVariableRef("jetscaledn"+year,
                                      sys_dijet_deta_scaledn[iyear]);
        bdt_reader->SetAltVariableRef("jetscaledn"+year,
                                      sys_dijet_dphi_scaledn[iyear]);
        bdt_reader->SetAltVariableRef("jetscaledn"+year,
                                      sys_lead_jet_pt_scaledn[iyear]);
        bdt_reader->SetAltVariableRef("jetscaledn"+year,
                                      sys_sublead_jet_pt_scaledn[iyear]);
        bdt_reader->SetAltVariableRef("jetscaledn"+year,
            sys_llphoton_refit_jet_balance_jetscaledn[iyear]);
        bdt_reader->SetAltVariableRef("jetscaledn"+year,
            sys_llphoton_refit_dijet_dphi_jetscaledn[iyear]);
        bdt_reader->SetAltVariableRef("jetscaledn"+year,
                                      sys_llphoton_refit_psi_default);
        bdt_reader->SetAltVariableRef("jetscaledn"+year,
            sys_lead_photon_drmin_default);
        bdt_reader->SetAltVariableRef("jetscaledn"+year,
            sys_lead_photon_drmax_default);
        bdt_reader->SetAltVariableRef("jetscaledn"+year,
                                      sys_llphoton_refit_costheta_default);
        bdt_reader->SetAltVariableRef("jetscaledn"+year,
                                      sys_llphoton_refit_cosTheta_default);
        bdt_reader->SetAltVariableRef("jetscaledn"+year,
                                      sys_llphoton_refit_relpt_default);
        bdt_reader->SetAltVariableRef("jetscaledn"+year,
                                      sys_lead_lepton_eta_default);
        bdt_reader->SetAltVariableRef("jetscaledn"+year,
                                      sys_sublead_lepton_eta_default);
        bdt_reader->SetAltVariableRef("jetscaledn"+year,
            sys_lead_photon_eta_default);
        bdt_reader->SetAltVariableRef("jetscaledn"+year,
            sys_lead_photon_idmva_default);
        bdt_reader->SetAltVariableRef("jetscaledn"+year,
                                      sys_lead_photon_relpterr_default);
        bdt_reader->SetAltVariableRef("jetscaledn"+year,
                                      sys_photon_zeppenfeld_jetscaledn[iyear]);
        bdt_reader->SetAltVariableRef("jetscaledn"+year,
                                      sys_photon_jet1_dr_jetscaledn[iyear]);
        bdt_reader->SetAltVariableRef("jetscaledn"+year,
                                      sys_photon_jet2_dr_jetscaledn[iyear]);
        //jetresup
        bdt_reader->SetAltVariableRef("jetresup"+year,
                                      sys_dijet_deta_resup[iyear]);
        bdt_reader->SetAltVariableRef("jetresup"+year,
                                      sys_dijet_dphi_resup[iyear]);
        bdt_reader->SetAltVariableRef("jetresup"+year,
                                      sys_lead_jet_pt_resup[iyear]);
        bdt_reader->SetAltVariableRef("jetresup"+year,
                                      sys_sublead_jet_pt_resup[iyear]);
        bdt_reader->SetAltVariableRef("jetresup"+year,
            sys_llphoton_refit_jet_balance_jetresup[iyear]);
        bdt_reader->SetAltVariableRef("jetresup"+year,
            sys_llphoton_refit_dijet_dphi_jetresup[iyear]);
        bdt_reader->SetAltVariableRef("jetresup"+year,
                                      sys_llphoton_refit_psi_default);
        bdt_reader->SetAltVariableRef("jetresup"+year,
                                      sys_lead_photon_drmin_default);
        bdt_reader->SetAltVariableRef("jetresup"+year,
                                      sys_lead_photon_drmax_default);
        bdt_reader->SetAltVariableRef("jetresup"+year,
                                      sys_llphoton_refit_costheta_default);
        bdt_reader->SetAltVariableRef("jetresup"+year,
                                      sys_llphoton_refit_cosTheta_default);
        bdt_reader->SetAltVariableRef("jetresup"+year,
                                      sys_llphoton_refit_relpt_default);
        bdt_reader->SetAltVariableRef("jetresup"+year,
                                      sys_lead_lepton_eta_default);
        bdt_reader->SetAltVariableRef("jetresup"+year,
                                      sys_sublead_lepton_eta_default);
        bdt_reader->SetAltVariableRef("jetresup"+year,
                                      sys_lead_photon_eta_default);
        bdt_reader->SetAltVariableRef("jetresup"+year,
                                      sys_lead_photon_idmva_default);
        bdt_reader->SetAltVariableRef("jetresup"+year,
                                      sys_lead_photon_relpterr_default);
        bdt_reader->SetAltVariableRef("jetresup"+year,
                                      sys_photon_zeppenfeld_jetresup[iyear]);
        bdt_reader->SetAltVariableRef("jetresup"+year,
                                      sys_photon_jet1_dr_jetresup[iyear]);
        bdt_reader->SetAltVariableRef("jetresup"+year,
                                      sys_photon_jet2_dr_jetresup[iyear]);
        //jetresdn
        bdt_reader->SetAltVariableRef("jetresdn"+year,
                                      sys_dijet_deta_resdn[iyear]);
        bdt_reader->SetAltVariableRef("jetresdn"+year,
                                      sys_dijet_dphi_resdn[iyear]);
        bdt_reader->SetAltVariableRef("jetresdn"+year,
                                      sys_lead_jet_pt_resdn[iyear]);
        bdt_reader->SetAltVariableRef("jetresdn"+year,
                                      sys_sublead_jet_pt_resdn[iyear]);
        bdt_reader->SetAltVariableRef("jetresdn"+year,
            sys_llphoton_refit_jet_balance_jetresdn[iyear]);
        bdt_reader->SetAltVariableRef("jetresdn"+year,
            sys_llphoton_refit_dijet_dphi_jetresdn[iyear]);
        bdt_reader->SetAltVariableRef("jetresdn"+year,
                                      sys_llphoton_refit_psi_default);
        bdt_reader->SetAltVariableRef("jetresdn"+year,
                                      sys_lead_photon_drmin_default);
        bdt_reader->SetAltVariableRef("jetresdn"+year,
                                      sys_lead_photon_drmax_default);
        bdt_reader->SetAltVariableRef("jetresdn"+year,
                                      sys_llphoton_refit_costheta_default);
        bdt_reader->SetAltVariableRef("jetresdn"+year,
                                      sys_llphoton_refit_cosTheta_default);
        bdt_reader->SetAltVariableRef("jetresdn"+year,
                                      sys_llphoton_refit_relpt_default);
        bdt_reader->SetAltVariableRef("jetresdn"+year,
                                      sys_lead_lepton_eta_default);
        bdt_reader->SetAltVariableRef("jetresdn"+year,
                                      sys_sublead_lepton_eta_default);
        bdt_reader->SetAltVariableRef("jetresdn"+year,
                                      sys_lead_photon_eta_default);
        bdt_reader->SetAltVariableRef("jetresdn"+year,
                                      sys_lead_photon_idmva_default);
        bdt_reader->SetAltVariableRef("jetresdn"+year,
                                      sys_lead_photon_relpterr_default);
        bdt_reader->SetAltVariableRef("jetresdn"+year,
                                      sys_photon_zeppenfeld_jetresdn[iyear]);
        bdt_reader->SetAltVariableRef("jetresdn"+year,
                                      sys_photon_jet1_dr_jetresdn[iyear]);
        bdt_reader->SetAltVariableRef("jetresdn"+year,
                                      sys_photon_jet2_dr_jetresdn[iyear]);
      } //loop over years
    } //loop over folds
    vbf_bdt_readers[0]->BookMVA("/net/cms37/data1/rui/Training/dataset_2JClassic_pinnacles_run2p3/weights/TMVAClassification_BDT_0.weights.xml");
    vbf_bdt_readers[1]->BookMVA("/net/cms37/data1/rui/Training/dataset_2JClassic_pinnacles_run2p3/weights/TMVAClassification_BDT_1.weights.xml");
    vbf_bdt_readers[2]->BookMVA("/net/cms37/data1/rui/Training/dataset_2JClassic_pinnacles_run2p3/weights/TMVAClassification_BDT_2.weights.xml");
    vbf_bdt_readers[3]->BookMVA("/net/cms37/data1/rui/Training/dataset_2JClassic_pinnacles_run2p3/weights/TMVAClassification_BDT_3.weights.xml");
    return vbf_bdt_readers;
  }

  //returns thread safe NamedFunc that returns VBF score
  NamedFunc syst_vbf_bdt_score(string variation) {
    static mutex mutex_;
    static unordered_map<std::thread::id, vector<shared_ptr<MVAWrapper>>> 
      vbf_bdts;
    return NamedFunc("vbf_bdtscore_"+variation,[variation]
        (const Baby &b) -> NamedFunc::ScalarType{
      if (b.njet()<2) return -1.0;
      int bdt_index = (((b.event()%314159)+1)%4);
      std::thread::id thread_id = std::this_thread::get_id();
      shared_ptr<MVAWrapper> vbf_bdt;
      {
        lock_guard<mutex> lock(mutex_);
        if (vbf_bdts.count(thread_id)==0) {
          vbf_bdts[thread_id] = SystVbfBdts();
        }
        vbf_bdt = vbf_bdts[thread_id][bdt_index];
      }
      return vbf_bdt->GetDiscriminantScore(b, variation);
    }).EnableCaching(true);
  }

  NamedFunc ggfbdt2503_score_default(const vector<FastForest> &xgb_bdts) {
     return XGBoostBDTScoreCached(xgb_bdts,
         {&sys_lead_photon_idmva_default, 
          &sys_lead_photon_relpterr_default, &sys_lead_photon_drmax_default, 
          &sys_lead_photon_drmin_default, 
          &sys_llphoton_refit_relpt_default, 
          &sys_llphoton_refit_cosTheta_default, 
          &sys_llphoton_refit_costheta_default, 
          &sys_llphoton_refit_psi_default,
          &sys_lead_lepton_eta_default, &sys_sublead_lepton_eta_default, 
          &sys_lead_photon_eta_default});
  }

  NamedFunc ggfbdt2503_score_elscaleup(const vector<FastForest> &xgb_bdts) {
     return XGBoostBDTScoreCached(xgb_bdts,
         {&sys_lead_photon_idmva_default, 
          &sys_lead_photon_relpterr_default, &sys_lead_photon_drmax_default, 
          &sys_lead_photon_drmin_default, 
          &sys_llphoton_refit_relpt_elscaleup, 
          &sys_llphoton_refit_cosTheta_elscaleup, 
          &sys_llphoton_refit_costheta_elscaleup, 
          &sys_llphoton_refit_psi_elscaleup,
          &sys_lead_lepton_eta_elscaleup, &sys_sublead_lepton_eta_elscaleup, 
          &sys_lead_photon_eta_default});
  }

  NamedFunc ggfbdt2503_score_elscaledn(const vector<FastForest> &xgb_bdts) {
     return XGBoostBDTScoreCached(xgb_bdts,
         {&sys_lead_photon_idmva_default, 
          &sys_lead_photon_relpterr_default, &sys_lead_photon_drmax_default, 
          &sys_lead_photon_drmin_default, 
          &sys_llphoton_refit_relpt_elscaledn, 
          &sys_llphoton_refit_cosTheta_elscaledn, 
          &sys_llphoton_refit_costheta_elscaledn, 
          &sys_llphoton_refit_psi_elscaledn,
          &sys_lead_lepton_eta_elscaledn, &sys_sublead_lepton_eta_elscaledn, 
          &sys_lead_photon_eta_default});
  }

  NamedFunc ggfbdt2503_score_elresup(const vector<FastForest> &xgb_bdts) {
     return XGBoostBDTScoreCached(xgb_bdts,
         {&sys_lead_photon_idmva_default, 
          &sys_lead_photon_relpterr_default, &sys_lead_photon_drmax_default, 
          &sys_lead_photon_drmin_default, 
          &sys_llphoton_refit_relpt_elresup, 
          &sys_llphoton_refit_cosTheta_elresup, 
          &sys_llphoton_refit_costheta_elresup, 
          &sys_llphoton_refit_psi_elresup,
          &sys_lead_lepton_eta_elresup, &sys_sublead_lepton_eta_elresup, 
          &sys_lead_photon_eta_default});
  }

  NamedFunc ggfbdt2503_score_elresdn(const vector<FastForest> &xgb_bdts) {
     return XGBoostBDTScoreCached(xgb_bdts,
         {&sys_lead_photon_idmva_default, 
          &sys_lead_photon_relpterr_default, &sys_lead_photon_drmax_default, 
          &sys_lead_photon_drmin_default, 
          &sys_llphoton_refit_relpt_elresdn, 
          &sys_llphoton_refit_cosTheta_elresdn, 
          &sys_llphoton_refit_costheta_elresdn, 
          &sys_llphoton_refit_psi_elresdn,
          &sys_lead_lepton_eta_elresdn, &sys_sublead_lepton_eta_elresdn, 
          &sys_lead_photon_eta_default});
  }

  NamedFunc ggfbdt2503_score_muscaleup(const vector<FastForest> &xgb_bdts) {
     return XGBoostBDTScoreCached(xgb_bdts,
         {&sys_lead_photon_idmva_default, 
          &sys_lead_photon_relpterr_default, &sys_lead_photon_drmax_default, 
          &sys_lead_photon_drmin_default, 
          &sys_llphoton_refit_relpt_muscaleup, 
          &sys_llphoton_refit_cosTheta_muscaleup, 
          &sys_llphoton_refit_costheta_muscaleup, 
          &sys_llphoton_refit_psi_muscaleup,
          &sys_lead_lepton_eta_muscaleup, &sys_sublead_lepton_eta_muscaleup, 
          &sys_lead_photon_eta_default});
  }

  NamedFunc ggfbdt2503_score_muscaledn(const vector<FastForest> &xgb_bdts) {
     return XGBoostBDTScoreCached(xgb_bdts,
         {&sys_lead_photon_idmva_default, 
          &sys_lead_photon_relpterr_default, &sys_lead_photon_drmax_default, 
          &sys_lead_photon_drmin_default, 
          &sys_llphoton_refit_relpt_muscaledn, 
          &sys_llphoton_refit_cosTheta_muscaledn, 
          &sys_llphoton_refit_costheta_muscaledn, 
          &sys_llphoton_refit_psi_muscaledn,
          &sys_lead_lepton_eta_muscaledn, &sys_sublead_lepton_eta_muscaledn, 
          &sys_lead_photon_eta_default});
  }

  NamedFunc ggfbdt2503_score_muresup(const vector<FastForest> &xgb_bdts) {
     return XGBoostBDTScoreCached(xgb_bdts,
         {&sys_lead_photon_idmva_default, 
          &sys_lead_photon_relpterr_default, &sys_lead_photon_drmax_default, 
          &sys_lead_photon_drmin_default, 
          &sys_llphoton_refit_relpt_muresup, 
          &sys_llphoton_refit_cosTheta_muresup, 
          &sys_llphoton_refit_costheta_muresup, 
          &sys_llphoton_refit_psi_muresup,
          &sys_lead_lepton_eta_muresup, &sys_sublead_lepton_eta_muresup, 
          &sys_lead_photon_eta_default});
  }

  NamedFunc ggfbdt2503_score_muresdn(const vector<FastForest> &xgb_bdts) {
     return XGBoostBDTScoreCached(xgb_bdts,
         {&sys_lead_photon_idmva_default, 
          &sys_lead_photon_relpterr_default, &sys_lead_photon_drmax_default, 
          &sys_lead_photon_drmin_default, 
          &sys_llphoton_refit_relpt_muresdn, 
          &sys_llphoton_refit_cosTheta_muresdn, 
          &sys_llphoton_refit_costheta_muresdn, 
          &sys_llphoton_refit_psi_muresdn,
          &sys_lead_lepton_eta_muresdn, &sys_sublead_lepton_eta_muresdn, 
          &sys_lead_photon_eta_default});
  }

  NamedFunc ggfbdt2503_score_phscaleup(const vector<FastForest> &xgb_bdts) {
     return XGBoostBDTScoreCached(xgb_bdts,
         {&sys_lead_photon_idmva_scaleup, 
          &sys_lead_photon_relpterr_scaleup, &sys_lead_photon_drmax_scaleup, 
          &sys_lead_photon_drmin_scaleup, 
          &sys_llphoton_refit_relpt_phscaleup, 
          &sys_llphoton_refit_cosTheta_phscaleup, 
          &sys_llphoton_refit_costheta_phscaleup, 
          &sys_llphoton_refit_psi_phscaleup,
          &sys_lead_lepton_eta_default, &sys_sublead_lepton_eta_default, 
          &sys_lead_photon_eta_scaleup});
  }

  NamedFunc ggfbdt2503_score_phscaledn(const vector<FastForest> &xgb_bdts) {
     return XGBoostBDTScoreCached(xgb_bdts,
         {&sys_lead_photon_idmva_scaledn, 
          &sys_lead_photon_relpterr_scaledn, &sys_lead_photon_drmax_scaledn, 
          &sys_lead_photon_drmin_scaledn, 
          &sys_llphoton_refit_relpt_phscaledn, 
          &sys_llphoton_refit_cosTheta_phscaledn, 
          &sys_llphoton_refit_costheta_phscaledn, 
          &sys_llphoton_refit_psi_phscaledn,
          &sys_lead_lepton_eta_default, &sys_sublead_lepton_eta_default, 
          &sys_lead_photon_eta_scaledn});
  }

  NamedFunc ggfbdt2503_score_phresup(const vector<FastForest> &xgb_bdts) {
     return XGBoostBDTScoreCached(xgb_bdts,
         {&sys_lead_photon_idmva_resup, 
          &sys_lead_photon_relpterr_resup, &sys_lead_photon_drmax_resup, 
          &sys_lead_photon_drmin_resup, 
          &sys_llphoton_refit_relpt_phresup, 
          &sys_llphoton_refit_cosTheta_phresup, 
          &sys_llphoton_refit_costheta_phresup, 
          &sys_llphoton_refit_psi_phresup,
          &sys_lead_lepton_eta_default, &sys_sublead_lepton_eta_default, 
          &sys_lead_photon_eta_resup});
  }

  NamedFunc ggfbdt2503_score_phresdn(const vector<FastForest> &xgb_bdts) {
     return XGBoostBDTScoreCached(xgb_bdts,
         {&sys_lead_photon_idmva_resdn, 
          &sys_lead_photon_relpterr_resdn, &sys_lead_photon_drmax_resdn, 
          &sys_lead_photon_drmin_resdn, 
          &sys_llphoton_refit_relpt_phresdn, 
          &sys_llphoton_refit_cosTheta_phresdn, 
          &sys_llphoton_refit_costheta_phresdn, 
          &sys_llphoton_refit_psi_phresdn,
          &sys_lead_lepton_eta_default, &sys_sublead_lepton_eta_default, 
          &sys_lead_photon_eta_resdn});
  }

  //new (25-08) BDTs

}
