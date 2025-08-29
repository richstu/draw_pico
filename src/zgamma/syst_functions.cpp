#include <memory>
#include <string>
#include <vector>

#include "TLorentzVector.h"
#include "TMath.h"

#include "core/baby.hpp"
#include "core/fastforest.hpp"
#include "core/named_func.hpp"
#include "core/named_func_utilities.hpp"
#include "zgamma/KinZfitter.hpp"
#include "zgamma/syst_functions.hpp"
#include "zgamma/zg_functions.hpp"
#include "zgamma/zg_utilities.hpp"

using std::make_shared;
using std::shared_ptr;
using std::string;
using std::vector;
using fastforest::FastForest;
using fastforest::FastForest;
using NamedFuncUtilities::FilterNamedFunc;
using NamedFuncUtilities::MultiReduceNamedFunc;
using NamedFuncUtilities::MultiMapNamedFunc;
using NamedFuncUtilities::ReduceNamedFunc;
using NamedFuncUtilities::reduce_sum;
using NamedFuncUtilities::reduce_max;
using NamedFuncUtilities::reduce_maxfirst;
using NamedFuncUtilities::reduce_sublead;
using NamedFuncUtilities::reduce_subleadfirst;
using ZgUtilities::CalculateAngles;
using ZgUtilities::get_lep_custom_refit;
using ZgUtilities::get_btag_wp_deepjet;
using ZgUtilities::XGBoostBDTScore;

namespace ZgFunctions {

  //One KinZfitter to rule them all
  shared_ptr<KinZfitter> syst_kin_z_fitter = make_shared<KinZfitter>();

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
      unc = 2.000+0.006;
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

  //el_sig with electron scale variation up
  const NamedFunc sys_el_sig_scaleup = NamedFunc("(sys_el_pt_scaleup>7)"
      "&&el_idLoose").Name("sys_el_sig_scaleup").EnableCaching(true);

  //el_sig with electron scale variation down
  const NamedFunc sys_el_sig_scaledn = NamedFunc("(sys_el_pt_scaledn>7)"
      "&&el_idLoose").Name("sys_el_sig_scaledn").EnableCaching(true);

  //el_sig with electron resolution variation up
  const NamedFunc sys_el_sig_resup = NamedFunc("(sys_el_pt_resup>7)"
      "&&el_idLoose").Name("sys_el_sig_resup").EnableCaching(true);

  //el_sig with electron resolution variation down
  const NamedFunc sys_el_sig_resdn = NamedFunc("(sys_el_pt_resdn>7)"
      "&&el_idLoose").Name("sys_el_sig_resdn").EnableCaching(true);

  //mu_sig with muon scale variation up
  //TODO fix in next production
  const NamedFunc sys_mu_sig_scaleup = NamedFunc("mu_id&&(mu_reliso<0.35)"
      "&&(mu_sip3d<4)&&((mu_pt+mu_ptErr)>5)").Name("sys_mu_sig_scaleup")
      .EnableCaching(true);

  //mu_sig with muon scale variation down
  //TODO fix in next production
  const NamedFunc sys_mu_sig_scaledn = NamedFunc("mu_id&&(mu_reliso<0.35)"
      "&&(mu_sip3d<4)&&((mu_pt-mu_ptErr)>5)").Name("sys_mu_sig_scaledn")
      .EnableCaching(true);

  //mu_sig with muon resolution variation up
  //TODO fix in next production
  const NamedFunc sys_mu_sig_resup = NamedFunc("mu_id&&(mu_reliso<0.35)"
      "&&(mu_sip3d<4)&&((mu_pt+mu_ptErr)>5)").Name("sys_mu_sig_resup")
      .EnableCaching(true);

  //mu_sig with muon resolution variation down
  //TODO fix in next production
  const NamedFunc sys_mu_sig_resdn = NamedFunc("mu_id&&(mu_reliso<0.35)"
      "&&(mu_sip3d<4)&&((mu_pt-mu_ptErr)>5)").Name("sys_mu_sig_resdn")
      .EnableCaching(true);

  //nel with electron scale variation up
  const NamedFunc sys_nel_scaleup = ReduceNamedFunc(sys_el_sig_scaleup,
      reduce_sum).Name("sys_nel_scaleup").EnableCaching(true);

  //nel with electron scale variation down
  const NamedFunc sys_nel_scaledn = ReduceNamedFunc(sys_el_sig_scaledn,
      reduce_sum).Name("sys_nel_scaledn").EnableCaching(true);

  //nel with electron resolution variation up
  const NamedFunc sys_nel_resup = ReduceNamedFunc(sys_el_sig_resup,reduce_sum)
      .Name("sys_nel_resup").EnableCaching(true);

  //nel with electron resolution variation down
  const NamedFunc sys_nel_resdn = ReduceNamedFunc(sys_el_sig_resdn,reduce_sum)
      .Name("sys_nel_resdn").EnableCaching(true);

  //nmu with muon scale variation up
  const NamedFunc sys_nmu_scaleup = ReduceNamedFunc(sys_mu_sig_scaleup,
      reduce_sum).Name("sys_nmu_scaleup").EnableCaching(true);

  //nmu with muon scale variation down
  const NamedFunc sys_nmu_scaledn = ReduceNamedFunc(sys_mu_sig_scaledn,
      reduce_sum).Name("sys_nmu_scaledn").EnableCaching(true);

  //nmu with muon resolution variation up
  const NamedFunc sys_nmu_resup = ReduceNamedFunc(sys_mu_sig_resup,reduce_sum)
      .Name("sys_nmu_resup").EnableCaching(true);

  //nmu with muon resolution variation down
  const NamedFunc sys_nmu_resdn = ReduceNamedFunc(sys_mu_sig_resdn,reduce_sum)
      .Name("sys_nmu_resdn").EnableCaching(true);

  //nlep with electron scale variation up
  const NamedFunc sys_nlep_elscaleup = NamedFunc(sys_nel_scaleup+"nmu").Name(
      "sys_nlep_elscaleup").EnableCaching(true);

  //nlep with electron scale variation up
  const NamedFunc sys_nlep_elscaledn = NamedFunc(sys_nel_scaledn+"nmu").Name(
      "sys_nlep_elscaledn").EnableCaching(true);

  //nlep with electron resolution variation up
  const NamedFunc sys_nlep_elresup = NamedFunc(sys_nel_resup+"nmu").Name(
      "sys_nlep_elresup").EnableCaching(true);

  //nlep with electron resolution variation up
  const NamedFunc sys_nlep_elresdn = NamedFunc(sys_nel_resdn+"nmu").Name(
      "sys_nlep_elresdn").EnableCaching(true);

  //nlep with muon scale variation up
  const NamedFunc sys_nlep_muscaleup = NamedFunc("nel"+sys_nmu_scaleup).Name(
      "sys_nlep_muscaleup").EnableCaching(true);

  //nlep with muon scale variation down
  const NamedFunc sys_nlep_muscaledn = NamedFunc("nel"+sys_nmu_scaledn).Name(
      "sys_nlep_muscaledn").EnableCaching(true);

  //nlep with muon resolution variation up
  const NamedFunc sys_nlep_muresup = NamedFunc("nel"+sys_nmu_resup).Name(
      "sys_nlep_muresup").EnableCaching(true);

  //nlep with muon resolution variation down
  const NamedFunc sys_nlep_muresdn = NamedFunc("nel"+sys_nmu_resdn).Name(
      "sys_nlep_muresdn").EnableCaching(true);

  //Gets NamedFunc that is number of ll candidates with variation
  NamedFunc assign_variation_nll(const NamedFunc &variation_el_sig, 
                                 const NamedFunc &variation_mu_sig,
                                 const string &name) {
    return NamedFunc(("sys_nll_"+name).c_str(),[variation_el_sig, 
        variation_mu_sig] (const Baby &b) -> NamedFunc::ScalarType{
          float nll = 0;
          vector<double> el_sig = variation_el_sig.GetVector(b);
          vector<double> mu_sig = variation_mu_sig.GetVector(b);
          for (unsigned iel = 0; iel < b.el_charge()->size(); iel++) {
            if (el_sig[iel]) {
              for (unsigned iel2 = 0; iel2 < iel; iel2++) {
                if (el_sig[iel2]) {
                  if ((b.el_charge()->at(iel)+b.el_charge()->at(iel2))==0) {
                    nll += 1.0;
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

  //nll with electron scale variation up
  NamedFunc sys_nll_elscaleup = assign_variation_nll(sys_el_sig_scaleup, 
      "mu_sig", "elscaleup");

  //nll with electron scale variation down
  NamedFunc sys_nll_elscaledn = assign_variation_nll(sys_el_sig_scaledn, 
      "mu_sig", "elscaledn");

  //nll with electron resolution variation up
  NamedFunc sys_nll_elresup = assign_variation_nll(sys_el_sig_resup, 
      "mu_sig", "elresup");

  //nll with electron resolution variation down
  NamedFunc sys_nll_elresdn = assign_variation_nll(sys_el_sig_resdn, 
      "mu_sig", "elresdn");

  //nll with muon scale variation up
  NamedFunc sys_nll_muscaleup = assign_variation_nll("el_sig", 
      sys_mu_sig_scaleup, "muscaleup");

  //nll with muon scale variation down
  NamedFunc sys_nll_muscaledn = assign_variation_nll("el_sig", 
      sys_mu_sig_scaledn, "muscaledn");

  //nll with muon resolution variation up
  NamedFunc sys_nll_muresup = assign_variation_nll("el_sig", 
      sys_mu_sig_resup, "muresup");

  //nll with muon resolution variation down
  NamedFunc sys_nll_muresdn = assign_variation_nll("el_sig", 
      sys_mu_sig_resdn, "muresdn");

  //leading electron pt with electron scale variation up
  const NamedFunc sys_lead_el_pt_scaleup = ReduceNamedFunc(FilterNamedFunc(
      "sys_el_pt_scaleup",sys_el_sig_scaleup),reduce_max)
      .Name("sys_lead_el_pt_scaleup").EnableCaching(true);

  //subleading electron pt with electron scale variation up
  const NamedFunc sys_sublead_el_pt_scaleup = ReduceNamedFunc(FilterNamedFunc(
      "sys_el_pt_scaleup",sys_el_sig_scaleup),reduce_sublead)
      .Name("sys_sublead_el_pt_scaleup").EnableCaching(true);

  //leading electron pt with electron scale variation downn
  const NamedFunc sys_lead_el_pt_scaledn = ReduceNamedFunc(FilterNamedFunc(
      "sys_el_pt_scaledn",sys_el_sig_scaledn),reduce_max)
      .Name("sys_lead_el_pt_scaledn").EnableCaching(true);

  //subleading electron pt with electron scale variation down
  const NamedFunc sys_sublead_el_pt_scaledn = ReduceNamedFunc(FilterNamedFunc(
      "sys_el_pt_scaledn",sys_el_sig_scaledn),reduce_sublead)
      .Name("sys_sublead_el_pt_scaledn").EnableCaching(true);

  //leading electron pt with electron resolution variation up
  const NamedFunc sys_lead_el_pt_resup = ReduceNamedFunc(FilterNamedFunc(
      "sys_el_pt_resup",sys_el_sig_resup),reduce_max)
      .Name("sys_lead_el_pt_resup").EnableCaching(true);

  //subleading electron pt with electron resolution variation up
  const NamedFunc sys_sublead_el_pt_resup = ReduceNamedFunc(FilterNamedFunc(
      "sys_el_pt_resup",sys_el_sig_resup),reduce_sublead)
      .Name("sys_sublead_el_pt_resup").EnableCaching(true);

  //leading electron pt with electron resolution variation down
  const NamedFunc sys_lead_el_pt_resdn = ReduceNamedFunc(FilterNamedFunc(
      "sys_el_pt_resdn",sys_el_sig_resdn),reduce_max)
      .Name("sys_lead_el_pt_resdn").EnableCaching(true);

  //subleading electron pt with electron resolution variation down
  const NamedFunc sys_sublead_el_pt_resdn = ReduceNamedFunc(FilterNamedFunc(
      "sys_el_pt_resdn",sys_el_sig_resdn),reduce_sublead)
      .Name("sys_sublead_el_pt_resdn").EnableCaching(true);

  //leading muon pt with muon scale variation up
  //TODO fix in next production
  const NamedFunc sys_lead_mu_pt_scaleup = ReduceNamedFunc(FilterNamedFunc(
      "mu_pt+mu_ptErr",sys_mu_sig_scaleup),
      reduce_max).Name("sys_lead_mu_pt_scaleup").EnableCaching(true);

  //subleading muon pt with muon scale variation up
  //TODO fix in next production
  const NamedFunc sys_sublead_mu_pt_scaleup = ReduceNamedFunc(FilterNamedFunc(
      "mu_pt+mu_ptErr",sys_mu_sig_scaleup),
      reduce_sublead).Name("sys_sublead_mu_pt_scaleup").EnableCaching(true);

  //leading muon pt with muon scale variation down
  //TODO fix in next production
  const NamedFunc sys_lead_mu_pt_scaledn = ReduceNamedFunc(FilterNamedFunc(
      "mu_pt-mu_ptErr",sys_mu_sig_scaledn),
      reduce_max).Name("sys_lead_mu_pt_scaledn").EnableCaching(true);

  //subleading muon pt with muon scale variation down
  //TODO fix in next production
  const NamedFunc sys_sublead_mu_pt_scaledn = ReduceNamedFunc(FilterNamedFunc(
      "mu_pt-mu_ptErr",sys_mu_sig_scaledn),
      reduce_sublead).Name("sys_sublead_mu_pt_scaledn").EnableCaching(true);

  //leading muon pt with muon resolution variation up
  //TODO fix in next production
  const NamedFunc sys_lead_mu_pt_resup = ReduceNamedFunc(FilterNamedFunc(
      "mu_pt+mu_ptErr",sys_mu_sig_resup),
      reduce_max).Name("sys_lead_mu_pt_resup").EnableCaching(true);

  //subleading muon pt with muon resolution variation up
  //TODO fix in next production
  const NamedFunc sys_sublead_mu_pt_resup = ReduceNamedFunc(FilterNamedFunc(
      "mu_pt+mu_ptErr",sys_mu_sig_resup),
      reduce_sublead).Name("sys_sublead_mu_pt_resup").EnableCaching(true);

  //leading muon pt with muon resolution variation down
  //TODO fix in next production
  const NamedFunc sys_lead_mu_pt_resdn = ReduceNamedFunc(FilterNamedFunc(
      "mu_pt-mu_ptErr",sys_mu_sig_resdn),
      reduce_max).Name("sys_lead_mu_pt_resdn").EnableCaching(true);

  //subleading muon pt with muon resolution variation down
  //TODO fix in next production
  const NamedFunc sys_sublead_mu_pt_resdn = ReduceNamedFunc(FilterNamedFunc(
      "mu_pt-mu_ptErr",sys_mu_sig_resdn),
      reduce_sublead).Name("sys_sublead_mu_pt_resdn").EnableCaching(true);

  //Gets NamedFunc that is trigger+pT cut flag with variation
  NamedFunc assign_variation_trig_pt(const NamedFunc &variation_el_pt, 
                                     const NamedFunc &variation_el_sig,
                                     const NamedFunc &variation_mu_pt,
                                     const NamedFunc &variation_mu_sig,
                                     const string &name) {
    return NamedFunc(("sys_trig_pt_"+name).c_str(),[variation_el_pt, 
        variation_el_sig, variation_mu_pt, variation_mu_sig] (const Baby &b) 
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

  //trigger and pT cuts with electron scale variation up
  NamedFunc sys_trig_pt_elscaleup = assign_variation_trig_pt(
      "sys_el_pt_scaleup", sys_el_sig_scaleup, "mu_pt", "mu_sig", "elscaleup");

  //trigger and pT cuts with electron scale variation down
  NamedFunc sys_trig_pt_elscaledn = assign_variation_trig_pt(
      "sys_el_pt_scaledn", sys_el_sig_scaledn, "mu_pt", "mu_sig", "elscaledn");

  //trigger and pT cuts with electron resolution variation up
  NamedFunc sys_trig_pt_elresup = assign_variation_trig_pt(
      "sys_el_pt_resup", sys_el_sig_resup, "mu_pt", "mu_sig", "elresup");

  //trigger and pT cuts with electron resolution variation down
  NamedFunc sys_trig_pt_elresdn = assign_variation_trig_pt(
      "sys_el_pt_resdn", sys_el_sig_resdn, "mu_pt", "mu_sig", "elresdn");

  //trigger and pT cuts with muon scale variation up
  //TODO fix with next production
  NamedFunc sys_trig_pt_muscaleup = assign_variation_trig_pt("el_pt",
      "el_sig", "mu_pt+mu_ptErr", sys_mu_sig_scaleup, "muscaleup");

  //trigger and pT cuts with muon scale variation down
  //TODO fix with next production
  NamedFunc sys_trig_pt_muscaledn = assign_variation_trig_pt("el_pt",
      "el_sig", "mu_pt-mu_ptErr", sys_mu_sig_scaledn, "muscaledn");

  //trigger and pT cuts with muon resolution variation up
  //TODO fix with next production
  NamedFunc sys_trig_pt_muresup = assign_variation_trig_pt("el_pt",
      "el_sig", "mu_pt+mu_ptErr", sys_mu_sig_resup, "muresup");

  //trigger and pT cuts with muon resolution variation down
  //TODO fix with next production
  NamedFunc sys_trig_pt_muresdn = assign_variation_trig_pt("el_pt",
      "el_sig", "mu_pt-mu_ptErr", sys_mu_sig_resdn, "muresdn");

  //Gets NamedFunc that is max lep miniso with variation
  NamedFunc assign_variation_max_lep_miniso(const NamedFunc &variation_el_sig,
                                            const NamedFunc &variation_mu_sig,
                                            const string &name) {
    return NamedFunc(("sys_max_lep_miniso_"+name).c_str(),[variation_el_sig, 
        variation_mu_sig] (const Baby &b) -> NamedFunc::ScalarType{

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

  //max lep miniso with electron scale variation up
  NamedFunc sys_max_lep_miniso_elscaleup = assign_variation_max_lep_miniso(
      sys_el_sig_scaleup, "mu_sig", "elscaleup");

  //max lep miniso with electron scale variation down
  NamedFunc sys_max_lep_miniso_elscaledn = assign_variation_max_lep_miniso(
      sys_el_sig_scaledn, "mu_sig", "elscaledn");

  //max lep miniso with electron resolution variation up
  NamedFunc sys_max_lep_miniso_elresup = assign_variation_max_lep_miniso(
      sys_el_sig_resup, "mu_sig", "elresup");

  //max lep miniso with electron resolution variation down
  NamedFunc sys_max_lep_miniso_elresdn = assign_variation_max_lep_miniso(
      sys_el_sig_resdn, "mu_sig", "elresdn");

  //max lep miniso with muon scale variation up
  NamedFunc sys_max_lep_miniso_muscaleup = assign_variation_max_lep_miniso(
      "el_sig", sys_mu_sig_scaleup, "muscaleup");

  //max lep miniso with muon scale variation down
  NamedFunc sys_max_lep_miniso_muscaledn = assign_variation_max_lep_miniso(
      "el_sig", sys_mu_sig_scaledn, "muscaledn");

  //max lep miniso with muon resolution variation up
  NamedFunc sys_max_lep_miniso_muresup = assign_variation_max_lep_miniso(
      "el_sig", sys_mu_sig_resup, "muresup");

  //max lep miniso with muon resolution variation down
  NamedFunc sys_max_lep_miniso_muresdn = assign_variation_max_lep_miniso(
      "el_sig", sys_mu_sig_resdn, "muresdn");

  //for reference, photons failing origin eta cuts are dropped from pico list

  //photon_sig with photon scale variation up
  NamedFunc sys_photon_sig_scaleup = NamedFunc(("(sys_photon_pt_scaleup>15)&&"
      "photon_elveto&&(photon_drmin>0.3)&&photon_id80")&&!photon_sigel)
      .Name("sys_photon_sig_elscaleup").EnableCaching(true);

  //photon_sig with photon scale variation down
  NamedFunc sys_photon_sig_scaledn = NamedFunc(("(sys_photon_pt_scaledn>15)&&"
      "photon_elveto&&(photon_drmin>0.3)&&photon_id80")&&!photon_sigel)
      .Name("sys_photon_sig_elscaledn").EnableCaching(true);

  //photon_sig with photon resolution variation up
  NamedFunc sys_photon_sig_resup = NamedFunc(("(sys_photon_pt_resup>15)&&"
      "photon_elveto&&(photon_drmin>0.3)&&photon_id80")&&!photon_sigel)
      .Name("sys_photon_sig_elresup").EnableCaching(true);

  //photon_sig with photon resolution variation down
  NamedFunc sys_photon_sig_resdn = NamedFunc(("(sys_photon_pt_resdn>15)&&"
      "photon_elveto&&(photon_drmin>0.3)&&photon_id80")&&!photon_sigel)
      .Name("sys_photon_sig_elresdn").EnableCaching(true);

  //nphoton with photon scale variation up
  NamedFunc sys_nphoton_scaleup = ReduceNamedFunc(sys_photon_sig_scaleup,
      reduce_sum).Name("sys_nphoton_scaleup").EnableCaching(true);

  //nphoton with photon scale variation down
  NamedFunc sys_nphoton_scaledn = ReduceNamedFunc(sys_photon_sig_scaledn,
      reduce_sum).Name("sys_nphoton_scaledn").EnableCaching(true);

  //nphoton with photon resolution variation up
  NamedFunc sys_nphoton_resup = ReduceNamedFunc(sys_photon_sig_resup,
      reduce_sum).Name("sys_nphoton_resup").EnableCaching(true);

  //nphoton with photon resolution variation down
  NamedFunc sys_nphoton_resdn = ReduceNamedFunc(sys_photon_sig_resdn,
      reduce_sum).Name("sys_nphoton_resdn").EnableCaching(true);

  //signal photon pt with photon scale variation up
  NamedFunc sys_sig_photon_pt_scaleup = FilterNamedFunc("sys_photon_pt_scaleup",
      sys_photon_sig_scaleup).Name("sys_sig_photon_pt_scaleup")
      .EnableCaching(true);

  //signal photon pt with photon scale variation downn
  NamedFunc sys_sig_photon_pt_scaledn = FilterNamedFunc("sys_photon_pt_scaledn",
      sys_photon_sig_scaledn).Name("sys_sig_photon_pt_scaledn")
      .EnableCaching(true);

  //signal photon pt with photon resolution variation up
  NamedFunc sys_sig_photon_pt_resup = FilterNamedFunc("sys_photon_pt_resup",
      sys_photon_sig_resup).Name("sys_sig_photon_pt_resup")
      .EnableCaching(true);

  //signal photon pt with photon resolution variation down
  NamedFunc sys_sig_photon_pt_resdn = FilterNamedFunc("sys_photon_pt_resdn",
      sys_photon_sig_resdn).Name("sys_sig_photon_pt_resdn")
      .EnableCaching(true);

  //signal photon eta with photon scale variation up
  NamedFunc sys_sig_photon_eta_scaleup = FilterNamedFunc("photon_eta",
      sys_photon_sig_scaleup).Name("sys_sig_photon_eta_scaleup")
      .EnableCaching(true);

  //signal photon eta with photon scale variation downn
  NamedFunc sys_sig_photon_eta_scaledn = FilterNamedFunc("photon_eta",
      sys_photon_sig_scaledn).Name("sys_sig_photon_eta_scaledn")
      .EnableCaching(true);

  //signal photon eta with photon resolution variation up
  NamedFunc sys_sig_photon_eta_resup = FilterNamedFunc("photon_eta",
      sys_photon_sig_resup).Name("sys_sig_photon_eta_resup")
      .EnableCaching(true);

  //signal photon eta with photon resolution variation down
  NamedFunc sys_sig_photon_eta_resdn = FilterNamedFunc("photon_eta",
      sys_photon_sig_resdn).Name("sys_sig_photon_eta_resdn")
      .EnableCaching(true);

  //signal photon phi with photon scale variation up
  NamedFunc sys_sig_photon_phi_scaleup = FilterNamedFunc("photon_phi",
      sys_photon_sig_scaleup).Name("sys_sig_photon_phi_scaleup")
      .EnableCaching(true);

  //signal photon phi with photon scale variation downn
  NamedFunc sys_sig_photon_phi_scaledn = FilterNamedFunc("photon_phi",
      sys_photon_sig_scaledn).Name("sys_sig_photon_phi_scaledn")
      .EnableCaching(true);

  //signal photon phi with photon resolution variation up
  NamedFunc sys_sig_photon_phi_resup = FilterNamedFunc("photon_phi",
      sys_photon_sig_resup).Name("sys_sig_photon_phi_resup")
      .EnableCaching(true);

  //signal photon phi with photon resolution variation down
  NamedFunc sys_sig_photon_phi_resdn = FilterNamedFunc("photon_phi",
      sys_photon_sig_resdn).Name("sys_sig_photon_phi_resdn")
      .EnableCaching(true);

  //signal photon drmin with photon scale variation up
  NamedFunc sys_sig_photon_drmin_scaleup = FilterNamedFunc("photon_drmin",
      sys_photon_sig_scaleup).Name("sys_sig_photon_drmin_scaleup")
      .EnableCaching(true);

  //signal photon drmin with photon scale variation downn
  NamedFunc sys_sig_photon_drmin_scaledn = FilterNamedFunc("photon_drmin",
      sys_photon_sig_scaledn).Name("sys_sig_photon_drmin_scaledn")
      .EnableCaching(true);

  //signal photon drmin with photon resolution variation up
  NamedFunc sys_sig_photon_drmin_resup = FilterNamedFunc("photon_drmin",
      sys_photon_sig_resup).Name("sys_sig_photon_drmin_resup")
      .EnableCaching(true);

  //signal photon drmin with photon resolution variation down
  NamedFunc sys_sig_photon_drmin_resdn = FilterNamedFunc("photon_drmin",
      sys_photon_sig_resdn).Name("sys_sig_photon_drmin_resdn")
      .EnableCaching(true);

  //signal photon drmax with photon scale variation up
  NamedFunc sys_sig_photon_drmax_scaleup = FilterNamedFunc("photon_drmax",
      sys_photon_sig_scaleup).Name("sys_sig_photon_drmax_scaleup")
      .EnableCaching(true);

  //signal photon drmax with photon scale variation downn
  NamedFunc sys_sig_photon_drmax_scaledn = FilterNamedFunc("photon_drmax",
      sys_photon_sig_scaledn).Name("sys_sig_photon_drmax_scaledn")
      .EnableCaching(true);

  //signal photon drmax with photon resolution variation up
  NamedFunc sys_sig_photon_drmax_resup = FilterNamedFunc("photon_drmax",
      sys_photon_sig_resup).Name("sys_sig_photon_drmax_resup")
      .EnableCaching(true);

  //signal photon drmax with photon resolution variation down
  NamedFunc sys_sig_photon_drmax_resdn = FilterNamedFunc("photon_drmax",
      sys_photon_sig_resdn).Name("sys_sig_photon_drmax_resdn")
      .EnableCaching(true);

  //signal photon idmva with photon scale variation up
  NamedFunc sys_sig_photon_idmva_scaleup = FilterNamedFunc("photon_idmva",
      sys_photon_sig_scaleup).Name("sys_sig_photon_idmva_scaleup")
      .EnableCaching(true);

  //signal photon idmva with photon scale variation downn
  NamedFunc sys_sig_photon_idmva_scaledn = FilterNamedFunc("photon_idmva",
      sys_photon_sig_scaledn).Name("sys_sig_photon_idmva_scaledn")
      .EnableCaching(true);

  //signal photon idmva with photon resolution variation up
  NamedFunc sys_sig_photon_idmva_resup = FilterNamedFunc("photon_idmva",
      sys_photon_sig_resup).Name("sys_sig_photon_idmva_resup")
      .EnableCaching(true);

  //signal photon idmva with photon resolution variation down
  NamedFunc sys_sig_photon_idmva_resdn = FilterNamedFunc("photon_idmva",
      sys_photon_sig_resdn).Name("sys_sig_photon_idmva_resdn")
      .EnableCaching(true);

  //leading photon pt with photon scale variation up
  NamedFunc sys_lead_photon_pt_scaleup = ReduceNamedFunc(
      sys_sig_photon_pt_scaleup, reduce_max).Name("sys_lead_photon_pt_scaleup")
      .EnableCaching(true);

  //leading photon pt with photon scale variation downn
  NamedFunc sys_lead_photon_pt_scaledn = ReduceNamedFunc(
      sys_sig_photon_pt_scaledn, reduce_max).Name("sys_lead_photon_pt_scaledn")
      .EnableCaching(true);

  //leading photon pt with photon resolution variation up
  NamedFunc sys_lead_photon_pt_resup = ReduceNamedFunc(sys_sig_photon_pt_resup,
      reduce_max).Name("sys_lead_photon_pt_resup").EnableCaching(true);

  //leading photon pt with photon resolution variation down
  NamedFunc sys_lead_photon_pt_resdn = ReduceNamedFunc(sys_sig_photon_pt_resdn,
      reduce_max).Name("sys_lead_photon_pt_resdn").EnableCaching(true);

  //leading photon eta with photon scale variation up
  NamedFunc sys_lead_photon_eta_scaleup = MultiReduceNamedFunc(
      {sys_sig_photon_pt_scaleup,sys_sig_photon_eta_scaleup}, reduce_maxfirst)
      .Name("sys_lead_photon_eta_scaleup").EnableCaching(true);

  //leading photon eta with photon scale variation down
  NamedFunc sys_lead_photon_eta_scaledn = MultiReduceNamedFunc(
      {sys_sig_photon_pt_scaledn,sys_sig_photon_eta_scaledn}, reduce_maxfirst)
      .Name("sys_lead_photon_eta_scaledn").EnableCaching(true);

  //leading photon eta with photon resolution variation up
  NamedFunc sys_lead_photon_eta_resup = MultiReduceNamedFunc(
      {sys_sig_photon_pt_resup,sys_sig_photon_eta_resup}, reduce_maxfirst)
      .Name("sys_lead_photon_eta_resup").EnableCaching(true);

  //leading photon eta with photon resolution variation down
  NamedFunc sys_lead_photon_eta_resdn = MultiReduceNamedFunc(
      {sys_sig_photon_pt_resdn,sys_sig_photon_eta_resdn}, reduce_maxfirst)
      .Name("sys_lead_photon_eta_resdn").EnableCaching(true);

  //leading photon phi with photon scale variation up
  NamedFunc sys_lead_photon_phi_scaleup = MultiReduceNamedFunc(
      {sys_sig_photon_pt_scaleup,sys_sig_photon_phi_scaleup}, reduce_maxfirst)
      .Name("sys_lead_photon_phi_scaleup").EnableCaching(true);

  //leading photon phi with photon scale variation down
  NamedFunc sys_lead_photon_phi_scaledn = MultiReduceNamedFunc(
      {sys_sig_photon_pt_scaledn,sys_sig_photon_phi_scaledn}, reduce_maxfirst)
      .Name("sys_lead_photon_phi_scaledn").EnableCaching(true);

  //leading photon phi with photon resolution variation up
  NamedFunc sys_lead_photon_phi_resup = MultiReduceNamedFunc(
      {sys_sig_photon_pt_resup,sys_sig_photon_phi_resup}, reduce_maxfirst)
      .Name("sys_lead_photon_phi_resup").EnableCaching(true);

  //leading photon phi with photon resolution variation down
  NamedFunc sys_lead_photon_phi_resdn = MultiReduceNamedFunc(
      {sys_sig_photon_pt_resdn,sys_sig_photon_phi_resdn}, reduce_maxfirst)
      .Name("sys_lead_photon_phi_resdn").EnableCaching(true);

  //leading photon drmin with photon scale variation up
  const NamedFunc sys_lead_photon_drmin_scaleup = MultiReduceNamedFunc(
      {sys_sig_photon_pt_scaleup,sys_sig_photon_drmin_scaleup}, reduce_maxfirst)
      .Name("sys_lead_photon_drmin_scaleup").EnableCaching(true);

  //leading photon drmin with photon scale variation down
  const NamedFunc sys_lead_photon_drmin_scaledn = MultiReduceNamedFunc(
      {sys_sig_photon_pt_scaledn,sys_sig_photon_drmin_scaledn}, reduce_maxfirst)
      .Name("sys_lead_photon_drmin_scaledn").EnableCaching(true);

  //leading photon drmin with photon resolution variation up
  const NamedFunc sys_lead_photon_drmin_resup = MultiReduceNamedFunc(
      {sys_sig_photon_pt_resup,sys_sig_photon_drmin_resup}, reduce_maxfirst)
      .Name("sys_lead_photon_drmin_resup").EnableCaching(true);

  //leading photon drmax with photon scale variation up
  const NamedFunc sys_lead_photon_drmax_scaleup = MultiReduceNamedFunc(
      {sys_sig_photon_pt_scaleup,sys_sig_photon_drmax_scaleup}, reduce_maxfirst)
      .Name("sys_lead_photon_drmax_scaleup").EnableCaching(true);

  //leading photon drmax with photon scale variation down
  const NamedFunc sys_lead_photon_drmax_scaledn = MultiReduceNamedFunc(
      {sys_sig_photon_pt_scaledn,sys_sig_photon_drmax_scaledn}, reduce_maxfirst)
      .Name("sys_lead_photon_drmax_scaledn").EnableCaching(true);

  //leading photon drmax with photon resolution variation up
  const NamedFunc sys_lead_photon_drmax_resup = MultiReduceNamedFunc(
      {sys_sig_photon_pt_resup,sys_sig_photon_drmax_resup}, reduce_maxfirst)
      .Name("sys_lead_photon_drmax_resup").EnableCaching(true);

  //leading photon drmax with photon resolution variation down
  const NamedFunc sys_lead_photon_drmax_resdn = MultiReduceNamedFunc(
      {sys_sig_photon_pt_resdn,sys_sig_photon_drmax_resdn}, reduce_maxfirst)
      .Name("sys_lead_photon_drmax_resdn").EnableCaching(true);

  //leading photon idmva with photon scale variation up
  const NamedFunc sys_lead_photon_idmva_scaleup = MultiReduceNamedFunc(
      {sys_sig_photon_pt_scaleup,sys_sig_photon_idmva_scaleup}, reduce_maxfirst)
      .Name("sys_lead_photon_idmva_scaleup").EnableCaching(true);

  //leading photon idmva with photon scale variation down
  const NamedFunc sys_lead_photon_idmva_scaledn = MultiReduceNamedFunc(
      {sys_sig_photon_pt_scaledn,sys_sig_photon_idmva_scaledn}, reduce_maxfirst)
      .Name("sys_lead_photon_idmva_scaledn").EnableCaching(true);

  //leading photon idmva with photon resolution variation up
  const NamedFunc sys_lead_photon_idmva_resup = MultiReduceNamedFunc(
      {sys_sig_photon_pt_resup,sys_sig_photon_idmva_resup}, reduce_maxfirst)
      .Name("sys_lead_photon_idmva_resup").EnableCaching(true);

  //leading photon idmva with photon resolution variation down
  const NamedFunc sys_lead_photon_idmva_resdn = MultiReduceNamedFunc(
      {sys_sig_photon_pt_resdn,sys_sig_photon_idmva_resdn}, reduce_maxfirst)
      .Name("sys_lead_photon_idmva_resdn").EnableCaching(true);

  //leading photon drmin with photon resolution variation down
  const NamedFunc sys_lead_photon_drmin_resdn = MultiReduceNamedFunc(
      {sys_sig_photon_pt_resdn,sys_sig_photon_drmin_resdn}, reduce_maxfirst)
      .Name("sys_lead_photon_drmin_resdn").EnableCaching(true);

  //signal photon ID80 flag with photon scale variation up
  const NamedFunc sys_sig_photon_id80_scaleup = FilterNamedFunc("photon_id80",
      sys_photon_sig_scaleup).Name("sys_sig_photon_id80_scaleup");

  //signal photon ID80 flag with photon scale variation downn
  const NamedFunc sys_sig_photon_id80_scaledn = FilterNamedFunc("photon_id80",
      sys_photon_sig_scaledn).Name("sys_sig_photon_id80_scaledn");

  //signal photon ID80 flag with photon resolution variation up
  const NamedFunc sys_sig_photon_id80_resup = FilterNamedFunc("photon_id80",
      sys_photon_sig_resup).Name("sys_sig_photon_id80_resup");

  //signal photon ID80 flag with photon resolution variation down
  const NamedFunc sys_sig_photon_id80_resdn = FilterNamedFunc("photon_id80",
      sys_photon_sig_resdn).Name("sys_sig_photon_id80_resdn");

  //leading photon ID80 flag with photon scale variation up
  const NamedFunc sys_lead_photon_id80_scaleup = MultiReduceNamedFunc(
      {sys_sig_photon_pt_scaleup, sys_sig_photon_id80_scaleup}, reduce_maxfirst).Name("sys_lead_photon_id80_scaleup");

  //leading photon ID80 flag with photon scale variation down
  const NamedFunc sys_lead_photon_id80_scaledn = MultiReduceNamedFunc(
      {sys_sig_photon_pt_scaledn, sys_sig_photon_id80_scaledn}, reduce_maxfirst).Name("sys_lead_photon_id80_scaledn");

  //leading photon ID80 flag with photon resolution variation up
  const NamedFunc sys_lead_photon_id80_resup = MultiReduceNamedFunc(
      {sys_sig_photon_pt_resup, sys_sig_photon_id80_resup}, reduce_maxfirst).Name("sys_lead_photon_id80_resup");

  //leading photon ID80 flag with photon resolution variation down
  const NamedFunc sys_lead_photon_id80_resdn = MultiReduceNamedFunc(
      {sys_sig_photon_pt_resdn, sys_sig_photon_id80_resdn}, reduce_maxfirst).Name("sys_lead_photon_id80_resdn");

  //get photon rel pt error with variation
  NamedFunc assign_variation_lead_photon_relpterr(const NamedFunc &ph_sig,
      const NamedFunc &ph_pt, const string &name) {
    return NamedFunc(("sys_lead_photon_relpterr_"+name).c_str(),[ph_sig, ph_pt]
        (const Baby &b) -> NamedFunc::ScalarType{
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
  const NamedFunc sys_lead_photon_relpterr_scaleup = 
      assign_variation_lead_photon_relpterr(sys_photon_sig_scaleup, 
      "sys_photon_pt_scaleup", "scaleup");
  const NamedFunc sys_lead_photon_relpterr_scaledn = 
      assign_variation_lead_photon_relpterr(sys_photon_sig_scaledn, 
      "sys_photon_pt_scaledn", "scaledn");
  const NamedFunc sys_lead_photon_relpterr_resup = 
      assign_variation_lead_photon_relpterr(sys_photon_sig_resup, 
      "sys_photon_pt_resup", "resup");
  const NamedFunc sys_lead_photon_relpterr_resdn = 
      assign_variation_lead_photon_relpterr(sys_photon_sig_resdn, 
      "sys_photon_pt_resdn", "resdn");

  //default dilepton 4 momentum
  const NamedFunc ll_refit_4p("ll_refit_4p",
      [](const Baby &b) -> NamedFunc::VectorType{
    vector<double> momentum = {b.ll_refit_pt(), b.ll_refit_eta(), 
                               b.ll_refit_phi(), b.ll_refit_m()};
    return momentum;
  });

  //get Z candidate properties with variation
  //returns (pt, eta, phi, m, lepid, i1, i2, idx)
  NamedFunc assign_variation_ll(const NamedFunc &el_pt, 
                                const NamedFunc &el_sig,
                                const NamedFunc &mu_pt, 
                                const NamedFunc &mu_sig,
                                bool is_el,
                                const string &name) {
    return NamedFunc(("sys_ll_"+name).c_str(),[el_pt, el_sig, mu_pt, mu_sig, 
        is_el] (const Baby &b) -> NamedFunc::VectorType{
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

  //default dilepton properties
  NamedFunc sys_ll_default = NamedFunc("sys_ll_default",[](const Baby &b) 
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

  //dilepton properties with electron scale variation up
  NamedFunc sys_ll_elscaleup = assign_variation_ll("sys_el_pt_scaleup", 
      sys_el_sig_scaleup, "mu_pt", "mu_sig", true, "elscaleup");

  //dilepton properties with electron scale variation down
  NamedFunc sys_ll_elscaledn = assign_variation_ll("sys_el_pt_scaledn", 
      sys_el_sig_scaledn, "mu_pt", "mu_sig", true, "elscaledn");

  //dilepton properties with electron resolution variation up
  NamedFunc sys_ll_elresup = assign_variation_ll("sys_el_pt_resup", 
      sys_el_sig_resup, "mu_pt", "mu_sig", true, "elresup");

  //dilepton properties with electron resolution variation down
  NamedFunc sys_ll_elresdn = assign_variation_ll("sys_el_pt_resdn", 
      sys_el_sig_resdn, "mu_pt", "mu_sig", true, "elresdn");

  //dilepton properties with muon scale variation up
  //TODO: fix in next production
  NamedFunc sys_ll_muscaleup = assign_variation_ll("el_pt", "el_sig", 
      "mu_pt+mu_ptErr", sys_mu_sig_scaleup, false, "muscaleup");

  //dilepton properties with muon scale variation down
  //TODO: fix in next production
  NamedFunc sys_ll_muscaledn = assign_variation_ll("el_pt", "el_sig", 
      "mu_pt-mu_ptErr", sys_mu_sig_scaledn, false, "muscaledn");

  //dilepton properties with muon resolution variation up
  //TODO: fix in next production
  NamedFunc sys_ll_muresup = assign_variation_ll("el_pt", "el_sig", 
      "mu_pt+mu_ptErr", sys_mu_sig_resup, false, "muresup");

  //dilepton properties with muon resolution variation down
  //TODO: fix in next production
  NamedFunc sys_ll_muresdn = assign_variation_ll("el_pt", "el_sig", 
      "mu_pt-mu_ptErr", sys_mu_sig_resdn, false, "muresdn");

  //dilepton mass with electron scale variation up
  NamedFunc sys_ll_m_elscaleup = NamedFunc(sys_ll_elscaleup[3]).Name(
      "sys_ll_m_elscaleup").EnableCaching(true);

  //dilepton mass with electron scale variation down
  NamedFunc sys_ll_m_elscaledn = NamedFunc(sys_ll_elscaledn[3]).Name(
      "sys_ll_m_elscaledn").EnableCaching(true);

  //dilepton mass with electron resolution variation up
  NamedFunc sys_ll_m_elresup = NamedFunc(sys_ll_elresup[3]).Name(
      "sys_ll_m_elresup").EnableCaching(true);

  //dilepton mass with electron resolution variation down
  NamedFunc sys_ll_m_elresdn = NamedFunc(sys_ll_elresdn[3]).Name(
      "sys_ll_m_elresdn").EnableCaching(true);

  //dilepton mass with muon scale variation up
  NamedFunc sys_ll_m_muscaleup = NamedFunc(sys_ll_muscaleup[3]).Name(
      "sys_ll_m_muscaleup").EnableCaching(true);

  //dilepton mass with muon scale variation down
  NamedFunc sys_ll_m_muscaledn = NamedFunc(sys_ll_muscaledn[3]).Name(
      "sys_ll_m_muscaledn").EnableCaching(true);
  
  //dilepton mass with muon resolution variation up
  NamedFunc sys_ll_m_muresup = NamedFunc(sys_ll_muresup[3]).Name(
      "sys_ll_m_muresup").EnableCaching(true);

  //dilepton mass with muon resolution variation down
  NamedFunc sys_ll_m_muresdn = NamedFunc(sys_ll_muresdn[3]).Name(
      "sys_ll_m_muresdn").EnableCaching(true);

  //get lepton eta with variation
  NamedFunc assign_variation_lep_eta(const NamedFunc &ll, const bool lead,
      const string &name) {
    string place = "lead";
    if (!lead) place = "sublead";
    return NamedFunc(("sys_"+place+"_lep_eta_"+name).c_str(),[ll, lead] 
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
    return NamedFunc(("sys_lphoton_p4_"+name).c_str(),[ll_p4, lead_photon_pt, 
        lead_photon_eta, lead_photon_phi, var_pdgid] (const Baby &b) 
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

  //H four momentum with electron scale variation up 
  NamedFunc sys_llphoton_p4_elscaleup = assign_variation_llphoton_p4(
      sys_ll_elscaleup, "photon_pt[0]", "photon_eta[0]", "photon_phi[0]", 11, 
      "elscaleup");

  //H four momentum with electron scale variation down
  NamedFunc sys_llphoton_p4_elscaledn = assign_variation_llphoton_p4(
      sys_ll_elscaledn, "photon_pt[0]", "photon_eta[0]", "photon_phi[0]", 11, 
      "elscaledn");

  //H four momentum with electron resolution variation up 
  NamedFunc sys_llphoton_p4_elresup = assign_variation_llphoton_p4(
      sys_ll_elresup, "photon_pt[0]", "photon_eta[0]", "photon_phi[0]", 11, 
      "elresup");

  //H four momentum with electron resolution variation down
  NamedFunc sys_llphoton_p4_elresdn = assign_variation_llphoton_p4(
      sys_ll_elresdn, "photon_pt[0]", "photon_eta[0]", "photon_phi[0]", 11, 
      "elresdn");

  //H four momentum with muon scale variation up 
  NamedFunc sys_llphoton_p4_muscaleup = assign_variation_llphoton_p4(
      sys_ll_muscaleup, "photon_pt[0]", "photon_eta[0]", "photon_phi[0]", 13, 
      "muscaleup");

  //H four momentum with muon scale variation down
  NamedFunc sys_llphoton_p4_muscaledn = assign_variation_llphoton_p4(
      sys_ll_muscaledn, "photon_pt[0]", "photon_eta[0]", "photon_phi[0]", 13, 
      "muscaledn");

  //H four momentum with muon resolution variation up 
  NamedFunc sys_llphoton_p4_muresup = assign_variation_llphoton_p4(
      sys_ll_muresup, "photon_pt[0]", "photon_eta[0]", "photon_phi[0]", 13, 
      "muresup");

  //H four momentum with muon resolution variation down
  NamedFunc sys_llphoton_p4_muresdn = assign_variation_llphoton_p4(
      sys_ll_muresdn, "photon_pt[0]", "photon_eta[0]", "photon_phi[0]", 13, 
      "muresdn");

  //H four momentum with photon scale variation up 
  NamedFunc sys_llphoton_p4_phscaleup = assign_variation_llphoton_p4(
      sys_ll_default, sys_lead_photon_pt_scaleup, sys_lead_photon_eta_scaleup, 
      sys_lead_photon_phi_scaleup, 22, "phscaleup");

  //H four momentum with photon scale variation down
  NamedFunc sys_llphoton_p4_phscaledn = assign_variation_llphoton_p4(
      sys_ll_default, sys_lead_photon_pt_scaledn, sys_lead_photon_eta_scaledn, 
      sys_lead_photon_phi_scaledn, 22, "phscaledn");

  //H four momentum with photon resolution variation up 
  NamedFunc sys_llphoton_p4_phresup = assign_variation_llphoton_p4(
      sys_ll_default, sys_lead_photon_pt_resup, sys_lead_photon_eta_resup, 
      sys_lead_photon_phi_resup, 22, "phresup");

  //H four momentum with photon resolution variation down
  NamedFunc sys_llphoton_p4_phresdn = assign_variation_llphoton_p4(
      sys_ll_default, sys_lead_photon_pt_resdn, sys_lead_photon_eta_resdn, 
      sys_lead_photon_phi_resdn, 22, "phresdn");

  //Higgs candidate mass with electron scale variation up
  NamedFunc sys_llphoton_m_elscaleup = NamedFunc(sys_llphoton_p4_elscaleup[3])
      .Name("sys_llphoton_m_elscaleup").EnableCaching(true);

  //Higgs candidate mass with electron scale variation down
  NamedFunc sys_llphoton_m_elscaledn = NamedFunc(sys_llphoton_p4_elscaledn[3])
      .Name("sys_llphoton_m_elscaledn").EnableCaching(true);

  //Higgs candidate mass with electron resolution variation up
  NamedFunc sys_llphoton_m_elresup = NamedFunc(sys_llphoton_p4_elresup[3])
      .Name("sys_llphoton_m_elresup").EnableCaching(true);

  //Higgs candidate mass with electron resolution variation down
  NamedFunc sys_llphoton_m_elresdn = NamedFunc(sys_llphoton_p4_elresdn[3])
      .Name("sys_llphoton_m_elresdn").EnableCaching(true);

  //Higgs candidate mass with muon scale variation up
  NamedFunc sys_llphoton_m_muscaleup = NamedFunc(sys_llphoton_p4_muscaleup[3])
      .Name("sys_llphoton_m_muscaleup").EnableCaching(true);

  //Higgs candidate mass with muon scale variation down
  NamedFunc sys_llphoton_m_muscaledn = NamedFunc(sys_llphoton_p4_muscaledn[3])
      .Name("sys_llphoton_m_muscaledn").EnableCaching(true);

  //Higgs candidate mass with muon resolution variation up
  NamedFunc sys_llphoton_m_muresup = NamedFunc(sys_llphoton_p4_muresup[3])
      .Name("sys_llphoton_m_muresup").EnableCaching(true);

  //Higgs candidate mass with muon resolution variation down
  NamedFunc sys_llphoton_m_muresdn = NamedFunc(sys_llphoton_p4_muresdn[3])
      .Name("sys_llphoton_m_muresdn").EnableCaching(true);

  //Higgs candidate mass with photon scale variation up
  NamedFunc sys_llphoton_m_phscaleup = NamedFunc(sys_llphoton_p4_phscaleup[3])
      .Name("sys_llphoton_m_phscaleup").EnableCaching(true);

  //Higgs candidate mass with photon scale variation down
  NamedFunc sys_llphoton_m_phscaledn = NamedFunc(sys_llphoton_p4_phscaledn[3])
      .Name("sys_llphoton_m_phscaledn").EnableCaching(true);

  //Higgs candidate mass with photon resolution variation up
  NamedFunc sys_llphoton_m_phresup = NamedFunc(sys_llphoton_p4_phresup[3])
      .Name("sys_llphoton_m_phresup").EnableCaching(true);

  //Higgs candidate mass with photon resolution variation down
  NamedFunc sys_llphoton_m_phresdn = NamedFunc(sys_llphoton_p4_phresdn[3])
      .Name("sys_llphoton_m_phresdn").EnableCaching(true);

  //Higgs candidate pt with electron scale variation up
  NamedFunc sys_llphoton_pt_elscaleup = NamedFunc(
      sys_llphoton_p4_elscaleup[0.0]).Name("sys_llphoton_pt_elscaleup")
      .EnableCaching(true);

  //Higgs candidate pt with electron scale variation down
  NamedFunc sys_llphoton_pt_elscaledn = NamedFunc(
      sys_llphoton_p4_elscaledn[0.0]).Name("sys_llphoton_pt_elscaledn")
      .EnableCaching(true);

  //Higgs candidate pt with electron resolution variation up
  NamedFunc sys_llphoton_pt_elresup = NamedFunc(sys_llphoton_p4_elresup[0.0])
      .Name("sys_llphoton_pt_elresup").EnableCaching(true);

  //Higgs candidate pt with electron resolution variation down
  NamedFunc sys_llphoton_pt_elresdn = NamedFunc(sys_llphoton_p4_elresdn[0.0])
      .Name("sys_llphoton_pt_elresdn").EnableCaching(true);

  //Higgs candidate pt with muon scale variation up
  NamedFunc sys_llphoton_pt_muscaleup = NamedFunc(
      sys_llphoton_p4_muscaleup[0.0]).Name("sys_llphoton_pt_muscaleup")
      .EnableCaching(true);

  //Higgs candidate pt with muon scale variation down
  NamedFunc sys_llphoton_pt_muscaledn = NamedFunc(
      sys_llphoton_p4_muscaledn[0.0]).Name("sys_llphoton_pt_muscaledn")
      .EnableCaching(true);

  //Higgs candidate pt with muon resolution variation up
  NamedFunc sys_llphoton_pt_muresup = NamedFunc(sys_llphoton_p4_muresup[0.0])
      .Name("sys_llphoton_pt_muresup").EnableCaching(true);

  //Higgs candidate pt with muon resolution variation down
  NamedFunc sys_llphoton_pt_muresdn = NamedFunc(sys_llphoton_p4_muresdn[0.0])
      .Name("sys_llphoton_pt_muresdn").EnableCaching(true);

  //Higgs candidate pt with photon scale variation up
  NamedFunc sys_llphoton_pt_phscaleup = NamedFunc(
      sys_llphoton_p4_phscaleup[0.0]).Name("sys_llphoton_pt_phscaleup")
      .EnableCaching(true);

  //Higgs candidate pt with photon scale variation down
  NamedFunc sys_llphoton_pt_phscaledn = NamedFunc(
      sys_llphoton_p4_phscaledn[0.0]).Name("sys_llphoton_pt_phscaledn")
      .EnableCaching(true);

  //Higgs candidate pt with photon resolution variation up
  NamedFunc sys_llphoton_pt_phresup = NamedFunc(sys_llphoton_p4_phresup[0.0])
      .Name("sys_llphoton_pt_phresup").EnableCaching(true);

  //Higgs candidate pt with photon resolution variation down
  NamedFunc sys_llphoton_pt_phresdn = NamedFunc(sys_llphoton_p4_phresdn[0.0])
      .Name("sys_llphoton_pt_phresdn").EnableCaching(true);

  //Higgs candidate phi with electron scale variation up
  NamedFunc sys_llphoton_phi_elscaleup = NamedFunc(
      sys_llphoton_p4_elscaleup[2]).Name("sys_llphoton_phi_elscaleup")
      .EnableCaching(true);

  //Higgs candidate phi with electron scale variation down
  NamedFunc sys_llphoton_phi_elscaledn = NamedFunc(
      sys_llphoton_p4_elscaledn[2]).Name("sys_llphoton_phi_elscaledn")
      .EnableCaching(true);

  //Higgs candidate phi with electron resolution variation up
  NamedFunc sys_llphoton_phi_elresup = NamedFunc(sys_llphoton_p4_elresup[2])
      .Name("sys_llphoton_phi_elresup").EnableCaching(true);

  //Higgs candidate phi with electron resolution variation down
  NamedFunc sys_llphoton_phi_elresdn = NamedFunc(sys_llphoton_p4_elresdn[2])
      .Name("sys_llphoton_phi_elresdn").EnableCaching(true);

  //Higgs candidate phi with muon scale variation up
  NamedFunc sys_llphoton_phi_muscaleup = NamedFunc(
      sys_llphoton_p4_muscaleup[2]).Name("sys_llphoton_phi_muscaleup")
      .EnableCaching(true);

  //Higgs candidate phi with muon scale variation down
  NamedFunc sys_llphoton_phi_muscaledn = NamedFunc(
      sys_llphoton_p4_muscaledn[2]).Name("sys_llphoton_phi_muscaledn")
      .EnableCaching(true);

  //Higgs candidate phi with muon resolution variation up
  NamedFunc sys_llphoton_phi_muresup = NamedFunc(sys_llphoton_p4_muresup[2])
      .Name("sys_llphoton_phi_muresup").EnableCaching(true);

  //Higgs candidate phi with muon resolution variation down
  NamedFunc sys_llphoton_phi_muresdn = NamedFunc(sys_llphoton_p4_muresdn[2])
      .Name("sys_llphoton_phi_muresdn").EnableCaching(true);

  //Higgs candidate phi with photon scale variation up
  NamedFunc sys_llphoton_phi_phscaleup = NamedFunc(
      sys_llphoton_p4_phscaleup[2]).Name("sys_llphoton_phi_phscaleup")
      .EnableCaching(true);

  //Higgs candidate phi with photon scale variation down
  NamedFunc sys_llphoton_phi_phscaledn = NamedFunc(
      sys_llphoton_p4_phscaledn[2]).Name("sys_llphoton_phi_phscaledn")
      .EnableCaching(true);

  //Higgs candidate phi with photon resolution variation up
  NamedFunc sys_llphoton_phi_phresup = NamedFunc(sys_llphoton_p4_phresup[2])
      .Name("sys_llphoton_phi_phresup").EnableCaching(true);

  //Higgs candidate phi with photon resolution variation down
  NamedFunc sys_llphoton_phi_phresdn = NamedFunc(sys_llphoton_p4_phresdn[2])
      .Name("sys_llphoton_phi_phresdn").EnableCaching(true);

  //Gets NamedFunc that is (pt1, eta1, phi1, m1, pt2, eta2, phi2, m2) of lepton
  //refit pT with variation
  NamedFunc assign_variation_lep_refit(
      const NamedFunc &el_pt, 
      const NamedFunc &mu_pt, 
      const NamedFunc &ll_lepid, 
      const NamedFunc &ll_i1, 
      const NamedFunc &ll_i2, 
      shared_ptr<KinZfitter> kin_z_fitter, 
      const string &name) {
    //require nll>=1 
    return NamedFunc(("sys_lep_refit_"+name).c_str(),[el_pt, mu_pt, ll_lepid,
        ll_i1, ll_i2, kin_z_fitter] (const Baby &b) -> NamedFunc::VectorType{
          vector<double> lep_refit = get_lep_custom_refit(b, kin_z_fitter, 
              el_pt, mu_pt, ll_lepid, ll_i1, ll_i2);
          return lep_refit;
        }).EnableCaching(true);
  }

  //lepton 4 momentum without variation
  NamedFunc sys_lep_refit_default = assign_variation_lep_refit(
      "el_pt", "mu_pt", sys_ll_default[4], sys_ll_default[5],
      sys_ll_default[6], syst_kin_z_fitter, "default");

  //lepton 4 momentum with electron scale variation up
  NamedFunc sys_lep_refit_elscaleup = assign_variation_lep_refit(
      "sys_el_pt_scaleup", "mu_pt", sys_ll_elscaleup[4], sys_ll_elscaleup[5],
      sys_ll_elscaleup[6], syst_kin_z_fitter, "elscaleup");

  //lepton 4 momentum with electron scale variation down
  NamedFunc sys_lep_refit_elscaledn = assign_variation_lep_refit(
      "sys_el_pt_scaledn", "mu_pt", sys_ll_elscaledn[4], sys_ll_elscaledn[5],
      sys_ll_elscaledn[6], syst_kin_z_fitter, "elscaledn");

  //lepton 4 momentum with electron resolution variation up
  NamedFunc sys_lep_refit_elresup = assign_variation_lep_refit(
      "sys_el_pt_resup", "mu_pt", sys_ll_elresup[4], sys_ll_elresup[5],
      sys_ll_elresup[6], syst_kin_z_fitter, "elresup");

  //lepton 4 momentum with electron resolution variation down
  NamedFunc sys_lep_refit_elresdn = assign_variation_lep_refit(
      "sys_el_pt_resdn", "mu_pt", sys_ll_elresdn[4], sys_ll_elresdn[5],
      sys_ll_elresdn[6], syst_kin_z_fitter, "elresdn");

  //lepton 4 momentum with muon scale variation up
  //TODO fix with next production
  NamedFunc sys_lep_refit_muscaleup = assign_variation_lep_refit(
      "el_pt", "mu_pt+mu_ptErr", sys_ll_muscaleup[4], sys_ll_muscaleup[5],
      sys_ll_muscaleup[6], syst_kin_z_fitter, "muscaleup");

  //lepton 4 momentum with muon scale variation down
  //TODO fix with next production
  NamedFunc sys_lep_refit_muscaledn = assign_variation_lep_refit(
      "el_pt", "mu_pt-mu_ptErr", sys_ll_muscaledn[4], sys_ll_muscaledn[5],
      sys_ll_muscaledn[6], syst_kin_z_fitter, "muscaledn");

  //lepton 4 momentum with muon resolution variation up
  //TODO fix with next production
  NamedFunc sys_lep_refit_muresup = assign_variation_lep_refit(
      "el_pt", "mu_pt+mu_ptErr", sys_ll_muresup[4], sys_ll_muresup[5],
      sys_ll_muresup[6], syst_kin_z_fitter, "muresup");

  //lepton 4 momentum with muon resolution variation down
  //TODO fix with next production
  NamedFunc sys_lep_refit_muresdn = assign_variation_lep_refit(
      "el_pt", "mu_pt-mu_ptErr", sys_ll_muresdn[4], sys_ll_muresdn[5],
      sys_ll_muresdn[6], syst_kin_z_fitter, "muresdn");

  //Gets NamedFunc that is assigns Z candidate four momentum with variation
  NamedFunc assign_variation_ll_refit_p4(const NamedFunc &ll_idx, 
      const NamedFunc &lep_refit, const string &name) {
    return NamedFunc(("sys_ll_refit_p4_"+name).c_str(),[ll_idx, lep_refit]
        (const Baby &b) -> NamedFunc::VectorType{
          if (ll_idx.GetScalar(b)==0)
            return {b.ll_refit_pt(), b.ll_refit_eta(), 
                    b.ll_refit_phi(), b.ll_refit_m()};
          vector<double> lep_refit0 = lep_refit.GetVector(b);
          TLorentzVector l1, l2, ll;
          l1.SetPtEtaPhiM(lep_refit0[0], lep_refit0[1], lep_refit0[2], 
                          lep_refit0[3]);
          l2.SetPtEtaPhiM(lep_refit0[4], lep_refit0[5], lep_refit0[6], 
                          lep_refit0[7]);
          ll = l1+l2;
          return {ll.Pt(), ll.Eta(), ll.Phi(), 
                  ll.M()};
        }).EnableCaching(true);
  }

  //Z candidate refit 4 momentum without any variations
  NamedFunc sys_ll_refit_p4_default = assign_variation_ll_refit_p4("0", 
      sys_lep_refit_elscaleup, "default");

  //Z candidate refit 4 momentum with electron scale variation up
  NamedFunc sys_ll_refit_p4_elscaleup = 
      assign_variation_ll_refit_p4(sys_ll_elscaleup[7], 
      sys_lep_refit_elscaleup, "elscaleup");

  //Z candidate refit 4 momentum with electron scale variation down
  NamedFunc sys_ll_refit_p4_elscaledn = 
      assign_variation_ll_refit_p4(sys_ll_elscaledn[7], 
      sys_lep_refit_elscaledn, "elscaledn");

  //Z candidate refit 4 momentum with electron resolution variation up
  NamedFunc sys_ll_refit_p4_elresup = 
      assign_variation_ll_refit_p4(sys_ll_elresup[7], 
      sys_lep_refit_elresup, "elresup");

  //Z candidate refit 4 momentum with electron resolution variation down
  NamedFunc sys_ll_refit_p4_elresdn = 
      assign_variation_ll_refit_p4(sys_ll_elresdn[7], 
      sys_lep_refit_elresdn, "elresdn");

  //Z candidate refit 4 momentum with muon scale variation up
  NamedFunc sys_ll_refit_p4_muscaleup = 
      assign_variation_ll_refit_p4(sys_ll_muscaleup[7], 
      sys_lep_refit_muscaleup, "muscaleup");

  //Z candidate refit 4 momentum with muon scale variation down
  NamedFunc sys_ll_refit_p4_muscaledn = 
      assign_variation_ll_refit_p4(sys_ll_muscaledn[7], 
      sys_lep_refit_muscaledn, "muscaledn");

  //Z candidate refit 4 momentum with muon resolution variation up
  NamedFunc sys_ll_refit_p4_muresup = 
      assign_variation_ll_refit_p4(sys_ll_muresup[7], 
      sys_lep_refit_muresup, "muresup");

  //Z candidate refit 4 momentum with muon resolution variation down
  NamedFunc sys_ll_refit_p4_muresdn = 
      assign_variation_ll_refit_p4(sys_ll_muresdn[7], 
      sys_lep_refit_muresdn, "muresdn");

  //Gets NamedFunc that assigns Higgs four momentum with variation
  NamedFunc assign_variation_llphoton_refit_p4(
      const NamedFunc &ll_idx, const NamedFunc &ll_refit_p4, 
      const NamedFunc &lead_photon_pt, const NamedFunc &lead_photon_eta, 
      const NamedFunc &lead_photon_phi, bool is_phvar, const string &name) {
    return NamedFunc(("sys_llphoton_refit_p4_"+name).c_str(),[ll_idx, 
        ll_refit_p4, lead_photon_pt, lead_photon_eta, lead_photon_phi, 
        is_phvar] (const Baby &b) -> NamedFunc::VectorType{
          if (!is_phvar && ll_idx.GetScalar(b)==0) {
            return {b.llphoton_refit_pt(), b.llphoton_refit_eta(), 
                    b.llphoton_refit_phi(), b.llphoton_refit_m()};
          }
          TLorentzVector ll, photon, llphoton;
          vector<double> ll_refit_p = ll_refit_p4.GetVector(b);
          ll.SetPtEtaPhiM(ll_refit_p[0], ll_refit_p[1], ll_refit_p[2], 
                          ll_refit_p[3]);
          photon.SetPtEtaPhiM(lead_photon_pt.GetScalar(b), 
                              lead_photon_eta.GetScalar(b), 
                              lead_photon_phi.GetScalar(b), 0.0);
          llphoton = ll+photon;
          return {llphoton.Pt(), llphoton.Eta(), llphoton.Phi(), 
                  llphoton.M()};
        }).EnableCaching(true);
  }

  //H candidate refit 4 momentum with electron scale variation up
  NamedFunc sys_llphoton_refit_p4_elscaleup = 
      assign_variation_llphoton_refit_p4(sys_ll_elscaleup[7], 
      sys_ll_refit_p4_elscaleup, "photon_pt[0]", "photon_eta[0]", 
      "photon_phi[0]", false, "elscaleup");
  
  //H candidate refit 4 momentum with electron scale variation down
  NamedFunc sys_llphoton_refit_p4_elscaledn = 
      assign_variation_llphoton_refit_p4(sys_ll_elscaledn[7], 
      sys_ll_refit_p4_elscaledn, "photon_pt[0]", "photon_eta[0]", 
      "photon_phi[0]", false, "elscaledn");

  //H candidate refit 4 momentum with electron resolution variation up
  NamedFunc sys_llphoton_refit_p4_elresup = 
      assign_variation_llphoton_refit_p4(sys_ll_elresup[7], 
      sys_ll_refit_p4_elresup, "photon_pt[0]", "photon_eta[0]", 
      "photon_phi[0]", false, "elresup");
  
  //H candidate refit 4 momentum with electron resolution variation down
  NamedFunc sys_llphoton_refit_p4_elresdn = 
      assign_variation_llphoton_refit_p4(sys_ll_elresdn[7], 
      sys_ll_refit_p4_elresdn, "photon_pt[0]", "photon_eta[0]", 
      "photon_phi[0]", false, "elresdn");

  //H candidate refit 4 momentum with muon scale variation up
  NamedFunc sys_llphoton_refit_p4_muscaleup = 
      assign_variation_llphoton_refit_p4(sys_ll_muscaleup[7], 
      sys_ll_refit_p4_muscaleup, "photon_pt[0]", "photon_eta[0]", 
      "photon_phi[0]", false, "muscaleup");
  
  //H candidate refit 4 momentum with muon scale variation down
  NamedFunc sys_llphoton_refit_p4_muscaledn = 
      assign_variation_llphoton_refit_p4(sys_ll_muscaledn[7], 
      sys_ll_refit_p4_muscaledn, "photon_pt[0]", "photon_eta[0]", 
      "photon_phi[0]", false, "muscaledn");

  //H candidate refit 4 momentum with muon resolution variation up
  NamedFunc sys_llphoton_refit_p4_muresup = 
      assign_variation_llphoton_refit_p4(sys_ll_muresup[7], 
      sys_ll_refit_p4_muresup, "photon_pt[0]", "photon_eta[0]", 
      "photon_phi[0]", false, "muresup");
  
  //H candidate refit 4 momentum with muon resolution variation down
  NamedFunc sys_llphoton_refit_p4_muresdn = 
      assign_variation_llphoton_refit_p4(sys_ll_muresdn[7], 
      sys_ll_refit_p4_muresdn, "photon_pt[0]", "photon_eta[0]", 
      "photon_phi[0]", false, "muresdn");

  //H candidate refit 4 momentum with photon scale variation up
  NamedFunc sys_llphoton_refit_p4_phscaleup = 
      assign_variation_llphoton_refit_p4("0", sys_ll_refit_p4_default, 
      sys_lead_photon_pt_scaleup, sys_lead_photon_eta_scaleup, 
      sys_lead_photon_phi_scaleup, true, "phscaleup");

  //H candidate refit 4 momentum with photon scale variation dn
  NamedFunc sys_llphoton_refit_p4_phscaledn = 
      assign_variation_llphoton_refit_p4("0", sys_ll_refit_p4_default, 
      sys_lead_photon_pt_scaledn, sys_lead_photon_eta_scaledn, 
      sys_lead_photon_phi_scaledn, true, "phscaledn");

  //H candidate refit 4 momentum with photon resolution variation up
  NamedFunc sys_llphoton_refit_p4_phresup = 
      assign_variation_llphoton_refit_p4("0", sys_ll_refit_p4_default, 
      sys_lead_photon_pt_resup, sys_lead_photon_eta_resup, 
      sys_lead_photon_phi_resup, true, "phresup");

  //H candidate refit 4 momentum with photon resolution variation dn
  NamedFunc sys_llphoton_refit_p4_phresdn = 
      assign_variation_llphoton_refit_p4("0", sys_ll_refit_p4_default, 
      sys_lead_photon_pt_resdn, sys_lead_photon_eta_resdn, 
      sys_lead_photon_phi_resdn, true, "phresdn");

  //H candidate refit pt with electron scale variation up
  NamedFunc sys_llphoton_refit_pt_elscaleup = NamedFunc(
      sys_llphoton_refit_p4_elscaleup[0.0]).Name("sys_llphoton_refit_pt_elscaleup")
      .EnableCaching(true);

  //H candidate refit pt with electron scale variation down
  NamedFunc sys_llphoton_refit_pt_elscaledn = NamedFunc(
      sys_llphoton_refit_p4_elscaledn[0.0]).Name("sys_llphoton_refit_pt_elscaledn")
      .EnableCaching(true);

  //H candidate refit pt with electron resolution variation up
  NamedFunc sys_llphoton_refit_pt_elresup = NamedFunc(
      sys_llphoton_refit_p4_elresup[0.0]).Name("sys_llphoton_refit_pt_elresup")
      .EnableCaching(true);

  //H candidate refit pt with electron resolution variation down
  NamedFunc sys_llphoton_refit_pt_elresdn = NamedFunc(
      sys_llphoton_refit_p4_elresdn[0.0]).Name("sys_llphoton_refit_pt_elresdn")
      .EnableCaching(true);

  //H candidate refit pt with muon scale variation up
  NamedFunc sys_llphoton_refit_pt_muscaleup = NamedFunc(
      sys_llphoton_refit_p4_muscaleup[0.0]).Name("sys_llphoton_refit_pt_muscaleup")
      .EnableCaching(true);

  //H candidate refit pt with muon scale variation down
  NamedFunc sys_llphoton_refit_pt_muscaledn = NamedFunc(
      sys_llphoton_refit_p4_muscaledn[0.0]).Name("sys_llphoton_refit_pt_muscaledn")
      .EnableCaching(true);

  //H candidate refit pt with muon resolution variation up
  NamedFunc sys_llphoton_refit_pt_muresup = NamedFunc(
      sys_llphoton_refit_p4_muresup[0.0]).Name("sys_llphoton_refit_pt_muresup")
      .EnableCaching(true);

  //H candidate refit pt with muon resolution variation down
  NamedFunc sys_llphoton_refit_pt_muresdn = NamedFunc(
      sys_llphoton_refit_p4_muresdn[0.0]).Name("sys_llphoton_refit_pt_muresdn")
      .EnableCaching(true);

  //H candidate refit pt with photon scale variation up
  NamedFunc sys_llphoton_refit_pt_phscaleup = NamedFunc(
      sys_llphoton_refit_p4_phscaleup[0.0]).Name("sys_llphoton_refit_pt_phscaleup")
      .EnableCaching(true);

  //H candidate refit pt with photon scale variation down
  NamedFunc sys_llphoton_refit_pt_phscaledn = NamedFunc(
      sys_llphoton_refit_p4_phscaledn[0.0]).Name("sys_llphoton_refit_pt_phscaledn")
      .EnableCaching(true);

  //H candidate refit pt with photon resolution variation up
  NamedFunc sys_llphoton_refit_pt_phresup = NamedFunc(
      sys_llphoton_refit_p4_phresup[0.0]).Name("sys_llphoton_refit_pt_phresup")
      .EnableCaching(true);

  //H candidate refit pt with photon resolution variation down
  NamedFunc sys_llphoton_refit_pt_phresdn = NamedFunc(
      sys_llphoton_refit_p4_phresdn[0.0]).Name("sys_llphoton_refit_pt_phresdn")
      .EnableCaching(true);
  
  //H candidate refit mass with electron scale variation up
  NamedFunc sys_llphoton_refit_m_elscaleup = NamedFunc(
      sys_llphoton_refit_p4_elscaleup[3]).Name("sys_llphoton_refit_m_elscaleup")
      .EnableCaching(true);

  //H candidate refit mass with electron scale variation down
  NamedFunc sys_llphoton_refit_m_elscaledn = NamedFunc(
      sys_llphoton_refit_p4_elscaledn[3]).Name("sys_llphoton_refit_m_elscaledn")
      .EnableCaching(true);

  //H candidate refit mass with electron resolution variation up
  NamedFunc sys_llphoton_refit_m_elresup = NamedFunc(
      sys_llphoton_refit_p4_elresup[3]).Name("sys_llphoton_refit_m_elresup")
      .EnableCaching(true);

  //H candidate refit mass with electron resolution variation down
  NamedFunc sys_llphoton_refit_m_elresdn = NamedFunc(
      sys_llphoton_refit_p4_elresdn[3]).Name("sys_llphoton_refit_m_elresdn")
      .EnableCaching(true);

  //H candidate refit mass with muon scale variation up
  NamedFunc sys_llphoton_refit_m_muscaleup = NamedFunc(
      sys_llphoton_refit_p4_muscaleup[3]).Name("sys_llphoton_refit_m_muscaleup")
      .EnableCaching(true);

  //H candidate refit mass with muon scale variation down
  NamedFunc sys_llphoton_refit_m_muscaledn = NamedFunc(
      sys_llphoton_refit_p4_muscaledn[3]).Name("sys_llphoton_refit_m_muscaledn")
      .EnableCaching(true);

  //H candidate refit mass with muon resolution variation up
  NamedFunc sys_llphoton_refit_m_muresup = NamedFunc(
      sys_llphoton_refit_p4_muresup[3]).Name("sys_llphoton_refit_m_muresup")
      .EnableCaching(true);

  //H candidate refit mass with muon resolution variation down
  NamedFunc sys_llphoton_refit_m_muresdn = NamedFunc(
      sys_llphoton_refit_p4_muresdn[3]).Name("sys_llphoton_refit_m_muresdn")
      .EnableCaching(true);

  //H candidate refit mass with photon scale variation up
  NamedFunc sys_llphoton_refit_m_phscaleup = NamedFunc(
      sys_llphoton_refit_p4_phscaleup[3]).Name("sys_llphoton_refit_m_phscaleup")
      .EnableCaching(true);

  //H candidate refit mass with photon scale variation down
  NamedFunc sys_llphoton_refit_m_phscaledn = NamedFunc(
      sys_llphoton_refit_p4_phscaledn[3]).Name("sys_llphoton_refit_m_phscaledn")
      .EnableCaching(true);

  //H candidate refit mass with photon resolution variation up
  NamedFunc sys_llphoton_refit_m_phresup = NamedFunc(
      sys_llphoton_refit_p4_phresup[3]).Name("sys_llphoton_refit_m_phresup")
      .EnableCaching(true);

  //H candidate refit mass with photon resolution variation down
  NamedFunc sys_llphoton_refit_m_phresdn = NamedFunc(
      sys_llphoton_refit_p4_phresdn[3]).Name("sys_llphoton_refit_m_phresdn")
      .EnableCaching(true);

  //H candidate refit phi with electron scale variation up
  NamedFunc sys_llphoton_refit_phi_elscaleup = NamedFunc(
      sys_llphoton_refit_p4_elscaleup[2])
      .Name("sys_llphoton_refit_phi_elscaleup").EnableCaching(true);

  //H candidate refit phi with electron scale variation down
  NamedFunc sys_llphoton_refit_phi_elscaledn = NamedFunc(
      sys_llphoton_refit_p4_elscaledn[2])
      .Name("sys_llphoton_refit_phi_elscaledn").EnableCaching(true);

  //H candidate refit phi with electron resolution variation up
  NamedFunc sys_llphoton_refit_phi_elresup = NamedFunc(
      sys_llphoton_refit_p4_elresup[2]).Name("sys_llphoton_refit_phi_elresup")
      .EnableCaching(true);

  //H candidate refit phi with electron resolution variation down
  NamedFunc sys_llphoton_refit_phi_elresdn = NamedFunc(
      sys_llphoton_refit_p4_elresdn[2]).Name("sys_llphoton_refit_phi_elresdn")
      .EnableCaching(true);

  //H candidate refit phi with muon scale variation up
  NamedFunc sys_llphoton_refit_phi_muscaleup = NamedFunc(
      sys_llphoton_refit_p4_muscaleup[2])
      .Name("sys_llphoton_refit_phi_muscaleup").EnableCaching(true);

  //H candidate refit phi with muon scale variation down
  NamedFunc sys_llphoton_refit_phi_muscaledn = NamedFunc(
      sys_llphoton_refit_p4_muscaledn[2])
      .Name("sys_llphoton_refit_phi_muscaledn").EnableCaching(true);

  //H candidate refit phi with muon resolution variation up
  NamedFunc sys_llphoton_refit_phi_muresup = NamedFunc(
      sys_llphoton_refit_p4_muresup[2]).Name("sys_llphoton_refit_phi_muresup")
      .EnableCaching(true);

  //H candidate refit phi with muon resolution variation down
  NamedFunc sys_llphoton_refit_phi_muresdn = NamedFunc(
      sys_llphoton_refit_p4_muresdn[2]).Name("sys_llphoton_refit_phi_muresdn")
      .EnableCaching(true);

  //H candidate refit phi with photon scale variation up
  NamedFunc sys_llphoton_refit_phi_phscaleup = NamedFunc(
      sys_llphoton_refit_p4_phscaleup[2])
      .Name("sys_llphoton_refit_phi_phscaleup").EnableCaching(true);

  //H candidate refit phi with photon scale variation down
  NamedFunc sys_llphoton_refit_phi_phscaledn = NamedFunc(
      sys_llphoton_refit_p4_phscaledn[2])
      .Name("sys_llphoton_refit_phi_phscaledn").EnableCaching(true);

  //H candidate refit phi with photon resolution variation up
  NamedFunc sys_llphoton_refit_phi_phresup = NamedFunc(
      sys_llphoton_refit_p4_phresup[2]).Name("sys_llphoton_refit_phi_phresup")
      .EnableCaching(true);

  //H candidate refit phi with photon resolution variation down
  NamedFunc sys_llphoton_refit_phi_phresdn = NamedFunc(
      sys_llphoton_refit_p4_phresdn[2]).Name("sys_llphoton_refit_phi_phresdn")
      .EnableCaching(true);

  //Gets NamedFunc that assigns kinematic angles with variation
  NamedFunc assign_variation_llphoton_refit_angles(const NamedFunc &var_ll,
      const NamedFunc &lep_refit, const NamedFunc &var_photon_pt, 
      const NamedFunc &var_photon_eta, const NamedFunc &var_photon_phi, 
      const string name) {
    return NamedFunc(("sys_llphoton_refit_angles_"+name).c_str(),[var_ll, 
        lep_refit, var_photon_pt, var_photon_eta, var_photon_phi] 
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
      sys_lep_refit_default, "photon_pt[0]", "photon_eta[0]", 
      "photon_phi[0]", "default");
  const NamedFunc sys_llphoton_refit_angles_elscaleup = 
      assign_variation_llphoton_refit_angles(sys_ll_elscaleup, 
      sys_lep_refit_elscaleup, "photon_pt[0]", "photon_eta[0]", 
      "photon_phi[0]", "elscaleup");
  const NamedFunc sys_llphoton_refit_angles_elscaledn = 
      assign_variation_llphoton_refit_angles(sys_ll_elscaledn, 
      sys_lep_refit_elscaledn, "photon_pt[0]", "photon_eta[0]", 
      "photon_phi[0]", "elscaledn");
  const NamedFunc sys_llphoton_refit_angles_elresup = 
      assign_variation_llphoton_refit_angles(sys_ll_elresup, 
      sys_lep_refit_elresup, "photon_pt[0]", "photon_eta[0]", 
      "photon_phi[0]", "elresup");
  const NamedFunc sys_llphoton_refit_angles_elresdn = 
      assign_variation_llphoton_refit_angles(sys_ll_elresdn, 
      sys_lep_refit_elresdn, "photon_pt[0]", "photon_eta[0]", 
      "photon_phi[0]", "elresdn");
  const NamedFunc sys_llphoton_refit_angles_muscaleup = 
      assign_variation_llphoton_refit_angles(sys_ll_muscaleup, 
      sys_lep_refit_muscaleup, "photon_pt[0]", "photon_eta[0]", 
      "photon_phi[0]", "muscaleup");
  const NamedFunc sys_llphoton_refit_angles_muscaledn = 
      assign_variation_llphoton_refit_angles(sys_ll_muscaledn, 
      sys_lep_refit_muscaledn, "photon_pt[0]", "photon_eta[0]", 
      "photon_phi[0]", "muscaledn");
  const NamedFunc sys_llphoton_refit_angles_muresup = 
      assign_variation_llphoton_refit_angles(sys_ll_muresup, 
      sys_lep_refit_muresup, "photon_pt[0]", "photon_eta[0]", 
      "photon_phi[0]", "muresup");
  const NamedFunc sys_llphoton_refit_angles_muresdn = 
      assign_variation_llphoton_refit_angles(sys_ll_muresdn, 
      sys_lep_refit_muresdn, "photon_pt[0]", "photon_eta[0]", 
      "photon_phi[0]", "muresdn");
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
  const NamedFunc sys_llphoton_refit_cosTheta_default = NamedFunc(
      sys_llphoton_refit_angles_default[0.0]).Name(
      "sys_llphoton_refit_coscapTheta_default").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_cosTheta_elscaleup = NamedFunc(
      sys_llphoton_refit_angles_elscaleup[0.0]).Name(
      "sys_llphoton_refit_coscapTheta_elscaleup").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_cosTheta_elscaledn = NamedFunc(
      sys_llphoton_refit_angles_elscaledn[0.0]).Name(
      "sys_llphoton_refit_coscapTheta_elscaledn").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_cosTheta_elresup = NamedFunc(
      sys_llphoton_refit_angles_elresup[0.0]).Name(
      "sys_llphoton_refit_coscapTheta_elresup").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_cosTheta_elresdn = NamedFunc(
      sys_llphoton_refit_angles_elresdn[0.0]).Name(
      "sys_llphoton_refit_coscapTheta_elresdn").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_cosTheta_muscaleup = NamedFunc(
      sys_llphoton_refit_angles_muscaleup[0.0]).Name(
      "sys_llphoton_refit_coscapTheta_muscaleup").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_cosTheta_muscaledn = NamedFunc(
      sys_llphoton_refit_angles_muscaledn[0.0]).Name(
      "sys_llphoton_refit_coscapTheta_muscaledn").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_cosTheta_muresup = NamedFunc(
      sys_llphoton_refit_angles_muresup[0.0]).Name(
      "sys_llphoton_refit_coscapTheta_muresup").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_cosTheta_muresdn = NamedFunc(
      sys_llphoton_refit_angles_muresdn[0.0]).Name(
      "sys_llphoton_refit_coscapTheta_muresdn").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_cosTheta_phscaleup = NamedFunc(
      sys_llphoton_refit_angles_phscaleup[0.0]).Name(
      "sys_llphoton_refit_coscapTheta_phscaleup").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_cosTheta_phscaledn = NamedFunc(
      sys_llphoton_refit_angles_phscaledn[0.0]).Name(
      "sys_llphoton_refit_coscapTheta_phscaledn").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_cosTheta_phresup = NamedFunc(
      sys_llphoton_refit_angles_phresup[0.0]).Name(
      "sys_llphoton_refit_coscapTheta_phresup").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_cosTheta_phresdn = NamedFunc(
      sys_llphoton_refit_angles_phresdn[0.0]).Name(
      "sys_llphoton_refit_coscapTheta_phresdn").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_costheta_default = NamedFunc(
      sys_llphoton_refit_angles_default[1]).Name(
      "sys_llphoton_refit_costheta_default").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_costheta_elscaleup = NamedFunc(
      sys_llphoton_refit_angles_elscaleup[1]).Name(
      "sys_llphoton_refit_costheta_elscaleup").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_costheta_elscaledn = NamedFunc(
      sys_llphoton_refit_angles_elscaledn[1]).Name(
      "sys_llphoton_refit_costheta_elscaledn").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_costheta_elresup = NamedFunc(
      sys_llphoton_refit_angles_elresup[1]).Name(
      "sys_llphoton_refit_costheta_elresup").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_costheta_elresdn = NamedFunc(
      sys_llphoton_refit_angles_elresdn[1]).Name(
      "sys_llphoton_refit_costheta_elresdn").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_costheta_muscaleup = NamedFunc(
      sys_llphoton_refit_angles_muscaleup[1]).Name(
      "sys_llphoton_refit_costheta_muscaleup").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_costheta_muscaledn = NamedFunc(
      sys_llphoton_refit_angles_muscaledn[1]).Name(
      "sys_llphoton_refit_costheta_muscaledn").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_costheta_muresup = NamedFunc(
      sys_llphoton_refit_angles_muresup[1]).Name(
      "sys_llphoton_refit_costheta_muresup").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_costheta_muresdn = NamedFunc(
      sys_llphoton_refit_angles_muresdn[1]).Name(
      "sys_llphoton_refit_costheta_muresdn").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_costheta_phscaleup = NamedFunc(
      sys_llphoton_refit_angles_phscaleup[1]).Name(
      "sys_llphoton_refit_costheta_phscaleup").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_costheta_phscaledn = NamedFunc(
      sys_llphoton_refit_angles_phscaledn[1]).Name(
      "sys_llphoton_refit_costheta_phscaledn").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_costheta_phresup = NamedFunc(
      sys_llphoton_refit_angles_phresup[1]).Name(
      "sys_llphoton_refit_costheta_phresup").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_costheta_phresdn = NamedFunc(
      sys_llphoton_refit_angles_phresdn[1]).Name(
      "sys_llphoton_refit_costheta_phresdn").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_psi_default = NamedFunc(
      sys_llphoton_refit_angles_default[2]).Name(
      "sys_llphoton_refit_psi_default").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_psi_elscaleup = NamedFunc(
      sys_llphoton_refit_angles_elscaleup[2]).Name(
      "sys_llphoton_refit_psi_elscaleup").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_psi_elscaledn = NamedFunc(
      sys_llphoton_refit_angles_elscaledn[2]).Name(
      "sys_llphoton_refit_psi_elscaledn").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_psi_elresup = NamedFunc(
      sys_llphoton_refit_angles_elresup[2]).Name(
      "sys_llphoton_refit_psi_elresup").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_psi_elresdn = NamedFunc(
      sys_llphoton_refit_angles_elresdn[2]).Name(
      "sys_llphoton_refit_psi_elresdn").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_psi_muscaleup = NamedFunc(
      sys_llphoton_refit_angles_muscaleup[2]).Name(
      "sys_llphoton_refit_psi_muscaleup").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_psi_muscaledn = NamedFunc(
      sys_llphoton_refit_angles_muscaledn[2]).Name(
      "sys_llphoton_refit_psi_muscaledn").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_psi_muresup = NamedFunc(
      sys_llphoton_refit_angles_muresup[2]).Name(
      "sys_llphoton_refit_psi_muresup").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_psi_muresdn = NamedFunc(
      sys_llphoton_refit_angles_muresdn[2]).Name(
      "sys_llphoton_refit_psi_muresdn").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_psi_phscaleup = NamedFunc(
      sys_llphoton_refit_angles_phscaleup[2]).Name(
      "sys_llphoton_refit_psi_phscaleup").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_psi_phscaledn = NamedFunc(
      sys_llphoton_refit_angles_phscaledn[2]).Name(
      "sys_llphoton_refit_psi_phscaledn").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_psi_phresup = NamedFunc(
      sys_llphoton_refit_angles_phresup[2]).Name(
      "sys_llphoton_refit_psi_phresup").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_psi_phresdn = NamedFunc(
      sys_llphoton_refit_angles_phresdn[2]).Name(
      "sys_llphoton_refit_psi_phresdn").EnableCaching(true);

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
    return NamedFunc(("sys_dijet_"+name).c_str(), [jet_sig, jet_pt, jet_m] 
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

  double map_deltaphi(vector<double> phi_values) {
    return deltaPhi(phi_values[0], phi_values[1]);
  }

  double map_deltar(vector<double> inputs) {
    return deltaR(inputs[0], inputs[1], inputs[2], inputs[3]);
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
  NamedFunc assign_vec_variation_year_select(const NamedFunc &var_namedfunc,
      const NamedFunc &def_namedfunc, const string year, const string &name) {
    return NamedFunc((name+year).c_str(),[var_namedfunc, def_namedfunc, year]
        (const Baby &b) -> NamedFunc::VectorType{
          if (b.SampleTypeString()==year
              || b.SampleTypeString()==("-"+year))
            return var_namedfunc.GetVector(b);
          return def_namedfunc.GetVector(b); 
        }).EnableCaching(true);
  }
  
  //Gets scalar NamedFunc that has variation for one year
  NamedFunc assign_sca_variation_year_select(const NamedFunc &var_namedfunc,
      const NamedFunc &def_namedfunc, const string year, const string &name) {
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
    return NamedFunc(("sys_nbdfm_"+name).c_str(),[var_jet_isgood]
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
  vector<NamedFunc> sys_jet_pt_scaleup;
  vector<NamedFunc> sys_jet_pt_scaledn;
  vector<NamedFunc> sys_jet_pt_resup;
  vector<NamedFunc> sys_jet_pt_resdn;
  vector<NamedFunc> sys_jet_isgood_scaleup;
  vector<NamedFunc> sys_jet_isgood_scaledn;
  vector<NamedFunc> sys_jet_isgood_resup;
  vector<NamedFunc> sys_jet_isgood_resdn;
  vector<NamedFunc> sys_sig_jet_pt_scaleup;
  vector<NamedFunc> sys_sig_jet_pt_scaledn;
  vector<NamedFunc> sys_sig_jet_pt_resup;
  vector<NamedFunc> sys_sig_jet_pt_resdn;
  vector<NamedFunc> sys_sig_jet_eta_scaleup;
  vector<NamedFunc> sys_sig_jet_eta_scaledn;
  vector<NamedFunc> sys_sig_jet_eta_resup;
  vector<NamedFunc> sys_sig_jet_eta_resdn;
  vector<NamedFunc> sys_sig_jet_phi_scaleup;
  vector<NamedFunc> sys_sig_jet_phi_scaledn;
  vector<NamedFunc> sys_sig_jet_phi_resup;
  vector<NamedFunc> sys_sig_jet_phi_resdn;
  vector<NamedFunc> sys_sig_jet_m_scaleup;
  vector<NamedFunc> sys_sig_jet_m_scaledn;
  vector<NamedFunc> sys_sig_jet_m_resup;
  vector<NamedFunc> sys_sig_jet_m_resdn;
  vector<NamedFunc> sys_sig_jet_deepflav_scaleup;
  vector<NamedFunc> sys_sig_jet_deepflav_scaledn;
  vector<NamedFunc> sys_sig_jet_deepflav_resup;
  vector<NamedFunc> sys_sig_jet_deepflav_resdn;
  vector<NamedFunc> sys_lead_jet_pt_scaleup;
  vector<NamedFunc> sys_lead_jet_pt_scaledn;
  vector<NamedFunc> sys_lead_jet_pt_resup;
  vector<NamedFunc> sys_lead_jet_pt_resdn;
  vector<NamedFunc> sys_lead_jet_eta_scaleup;
  vector<NamedFunc> sys_lead_jet_eta_scaledn;
  vector<NamedFunc> sys_lead_jet_eta_resup;
  vector<NamedFunc> sys_lead_jet_eta_resdn;
  vector<NamedFunc> sys_lead_jet_phi_scaleup;
  vector<NamedFunc> sys_lead_jet_phi_scaledn;
  vector<NamedFunc> sys_lead_jet_phi_resup;
  vector<NamedFunc> sys_lead_jet_phi_resdn;
  vector<NamedFunc> sys_lead_jet_m_scaleup;
  vector<NamedFunc> sys_lead_jet_m_scaledn;
  vector<NamedFunc> sys_lead_jet_m_resup;
  vector<NamedFunc> sys_lead_jet_m_resdn;
  vector<NamedFunc> sys_sublead_jet_pt_scaleup;
  vector<NamedFunc> sys_sublead_jet_pt_scaledn;
  vector<NamedFunc> sys_sublead_jet_pt_resup;
  vector<NamedFunc> sys_sublead_jet_pt_resdn;
  vector<NamedFunc> sys_sublead_jet_eta_scaleup;
  vector<NamedFunc> sys_sublead_jet_eta_scaledn;
  vector<NamedFunc> sys_sublead_jet_eta_resup;
  vector<NamedFunc> sys_sublead_jet_eta_resdn;
  vector<NamedFunc> sys_sublead_jet_phi_scaleup;
  vector<NamedFunc> sys_sublead_jet_phi_scaledn;
  vector<NamedFunc> sys_sublead_jet_phi_resup;
  vector<NamedFunc> sys_sublead_jet_phi_resdn;
  vector<NamedFunc> sys_njet_scaleup;
  vector<NamedFunc> sys_njet_scaledn;
  vector<NamedFunc> sys_njet_resup;
  vector<NamedFunc> sys_njet_resdn;
  vector<NamedFunc> sys_nbdfm_scaleup;
  vector<NamedFunc> sys_nbdfm_scaledn;
  vector<NamedFunc> sys_nbdfm_resup;
  vector<NamedFunc> sys_nbdfm_resdn;
  vector<NamedFunc> sys_met_scaleup;
  vector<NamedFunc> sys_met_scaledn;
  vector<NamedFunc> sys_met_resup;
  vector<NamedFunc> sys_met_resdn;

  //dijet variations
  const NamedFunc sys_dijet_default = assign_variation_dijet("jet_isgood",
      "jet_pt","jet_m","default");
  const NamedFunc sys_dijet_m_default = (sys_dijet_default[3]).Name(
      "sys_dijet_m_default").EnableCaching(true);
  const NamedFunc sys_dijet_deta_default = (sys_dijet_default[4]).Name(
      "sys_dijet_deta_default").EnableCaching(true);
  const NamedFunc sys_dijet_dphi_default = (sys_dijet_default[5]).Name(
      "sys_dijet_dphi_default").EnableCaching(true);
  vector<NamedFunc> sys_dijet_scaleup;
  vector<NamedFunc> sys_dijet_scaledn;
  vector<NamedFunc> sys_dijet_resup;
  vector<NamedFunc> sys_dijet_resdn;
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
  NamedFunc pinnacles_run2_jet_isgood_nopt = NamedFunc("!jet_islep"
      "&&!jet_isphoton&&-4.7<jet_eta&&jet_eta<4.7&&jet_id&&!(((jet_eta>-3.0"
      "&&jet_eta<-2.5)||(jet_eta>2.5&&jet_eta<3.0))&&jet_pt<40)")
      .Name("jet_isgood_nopt").EnableCaching(true);

  //Gets NamedFunc that is untagged category selections with variation
  NamedFunc assign_variation_untagged_category(
      const NamedFunc &var_nlep, const NamedFunc &var_njet, 
      const NamedFunc &var_nbdfm, const NamedFunc &var_met,
      const NamedFunc &var_llphoton_pt, const NamedFunc &var_llphoton_m,
      const NamedFunc &var_max_lep_miniso, const NamedFunc &var_ll_m,
      const string &name) {

    return NamedFunc(
        //2l+b but not enough jets for tth
        (var_nlep==2&&var_nbdfm>=1&&var_njet<5)
        //3l0b, but not enough met for vh3l
        ||(var_nlep>=3&&var_nbdfm==0.0&&var_met<=30)
        //3l+b, but not enough jets for tth
        ||(var_nlep==3&&var_nbdfm>=1&&var_njet<3)
        //fail additional llphoton pt/m selections in vhmet
        ||(var_nlep==2&&var_njet<=1&&var_met>90
           &&(var_llphoton_pt/var_llphoton_m)<=0.4)
        //fail additional llphoton pt/m selections in vh3l
        ||(var_nlep==3&&var_nbdfm==0.0&&var_met>30
           &&(var_llphoton_pt/var_llphoton_m)<=0.3)
        //fail additional miniso selection in vh3l
        ||(var_nlep>=3&&var_nbdfm==0.0&&var_met>30&&var_max_lep_miniso>0.15)
        //fail additional mll selection in tthhad
        ||(var_nlep==2&&var_nbdfm>=1&&var_njet>=5&&(var_ll_m<85||var_ll_m>95))
        //fail additional miniso seleciton in tthlep
        ||(((var_nlep==3&&var_nbdfm>=1&&var_njet>=3)
            ||(var_nlep>=4&&var_nbdfm>=1))
           &&var_max_lep_miniso>=0.1)).Name("sys_untagged_category_"+name)
        .EnableCaching(true);
  }

  //untagged category selection with electron scale variation up
  const NamedFunc untagged_category_elscaleup = 
      assign_variation_untagged_category(
      sys_nlep_elscaleup, "njet", "nbdfm", "met", sys_llphoton_pt_elscaleup, 
      sys_llphoton_m_elscaleup, sys_max_lep_miniso_elscaleup, 
      sys_ll_m_elscaleup, "elscaleup");

  //untagged category selection with electron scale variation up
  const NamedFunc untagged_category_elscaledn = 
      assign_variation_untagged_category(
      sys_nlep_elscaledn, "njet", "nbdfm", "met", sys_llphoton_pt_elscaledn, 
      sys_llphoton_m_elscaledn, sys_max_lep_miniso_elscaledn, 
      sys_ll_m_elscaledn, "elscaledn");

  //untagged category selection with electron resolution variation up
  const NamedFunc untagged_category_elresup = 
      assign_variation_untagged_category(
      sys_nlep_elresup, "njet", "nbdfm", "met", sys_llphoton_pt_elresup, 
      sys_llphoton_m_elresup, sys_max_lep_miniso_elresup, 
      sys_ll_m_elresup, "elresup");

  //untagged category selection with electron resolution variation up
  const NamedFunc untagged_category_elresdn = 
      assign_variation_untagged_category(
      sys_nlep_elresdn, "njet", "nbdfm", "met", sys_llphoton_pt_elresdn, 
      sys_llphoton_m_elresdn, sys_max_lep_miniso_elresdn, 
      sys_ll_m_elresdn, "elresdn");

  //untagged category selection with muon scale variation up
  const NamedFunc untagged_category_muscaleup = 
      assign_variation_untagged_category(
      sys_nlep_muscaleup, "njet", "nbdfm", "met", sys_llphoton_pt_muscaleup, 
      sys_llphoton_m_muscaleup, sys_max_lep_miniso_muscaleup, 
      sys_ll_m_muscaleup, "muscaleup");

  //untagged category selection with muon scale variation up
  const NamedFunc untagged_category_muscaledn = 
      assign_variation_untagged_category(
      sys_nlep_muscaledn, "njet", "nbdfm", "met", sys_llphoton_pt_muscaledn, 
      sys_llphoton_m_muscaledn, sys_max_lep_miniso_muscaledn, 
      sys_ll_m_muscaledn, "muscaledn");

  //untagged category selection with muon resolution variation up
  const NamedFunc untagged_category_muresup = 
      assign_variation_untagged_category(
      sys_nlep_muresup, "njet", "nbdfm", "met", sys_llphoton_pt_muresup, 
      sys_llphoton_m_muresup, sys_max_lep_miniso_muresup, 
      sys_ll_m_muresup, "muresup");

  //untagged category selection with muon resolution variation up
  const NamedFunc untagged_category_muresdn = 
      assign_variation_untagged_category(
      sys_nlep_muresdn, "njet", "nbdfm", "met", sys_llphoton_pt_muresdn, 
      sys_llphoton_m_muresdn, sys_max_lep_miniso_muresdn, 
      sys_ll_m_muresdn, "muresdn");

  //untagged category selection with photon scale variation up
  const NamedFunc untagged_category_phscaleup = 
      assign_variation_untagged_category(
      "nlep", "njet", "nbdfm", "met", sys_llphoton_pt_phscaleup, 
      sys_llphoton_m_phscaleup, max_lep_miniso, "ll_m[0]", "phscaleup");

  //untagged category selection with photon scale variation down
  const NamedFunc untagged_category_phscaledn = 
      assign_variation_untagged_category(
      "nlep", "njet", "nbdfm", "met", sys_llphoton_pt_phscaledn, 
      sys_llphoton_m_phscaledn, max_lep_miniso, "ll_m[0]", "phscaledn");

  //untagged category selection with photon resolution variation up
  const NamedFunc untagged_category_phresup = 
      assign_variation_untagged_category(
      "nlep", "njet", "nbdfm", "met", sys_llphoton_pt_phresup, 
      sys_llphoton_m_phresup, max_lep_miniso, "ll_m[0]", "phresup");

  //untagged category selection with photon resolution variation down
  const NamedFunc untagged_category_phresdn = 
      assign_variation_untagged_category(
      "nlep", "njet", "nbdfm", "met", sys_llphoton_pt_phresdn, 
      sys_llphoton_m_phresdn, max_lep_miniso, "ll_m[0]", "phresdn");
  
  //untagged category selection with various jet variations
  vector<NamedFunc> untagged_category_jetscaleup;
  vector<NamedFunc> untagged_category_jetscaledn;
  vector<NamedFunc> untagged_category_jetresup;
  vector<NamedFunc> untagged_category_jetresdn;

  //Gets NamedFunc that is untagged category selections with variation
  NamedFunc assign_variation_llphoton_refit_jet_dphi(
      const NamedFunc &llphoton_refit_phi, const NamedFunc &jet_sig, 
      const NamedFunc &jet_pt, const string &name) {
    return NamedFunc(("sys_llphoton_refit_jet_dphi_"+name).c_str(),[
        llphoton_refit_phi,jet_sig, jet_pt](const Baby &b) 
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
        ll_refit_p4, var_lead_photon_pt, var_lead_photon_phi, jet_sig, jet_pt]
        (const Baby &b) -> NamedFunc::ScalarType{
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
    return NamedFunc(("sys_mht_"+name).c_str(),[var_el_sig, var_el_pt, 
        var_mu_sig, var_mu_pt, var_ph_sig, var_ph_pt, var_jet_sig, var_jet_pt]
        (const Baby &b) -> NamedFunc::VectorType{
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
        var_lead_photon_eta, var_njet, var_lead_jet_eta, var_sublead_jet_eta]
        (const Baby &b) -> NamedFunc::ScalarType{
      double evt_njet = var_njet.GetScalar(b);
      if (evt_njet < 1)
        return -999.0;
      //TODO fix 2 in next production
      if (evt_njet < 2)
        return fabs(var_lead_photon_eta.GetScalar(b)
                    -var_lead_jet_eta.GetScalar(b)/2.0);
      return fabs(var_lead_photon_eta.GetScalar(b)
          -(var_lead_jet_eta.GetScalar(b)+var_sublead_jet_eta.GetScalar(b))
          /2.0);
    }).EnableCaching(true);
  }

  //jet+llphoton variations
  const NamedFunc sys_llphoton_refit_dijet_dphi_default = MultiMapNamedFunc(
      {"llphoton_refit_phi", sys_dijet_default[2]}, map_deltaphi)
      .Name("sys_llphoton_dijet_dphi_default").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_dijet_dphi_elscaleup = MultiMapNamedFunc(
      {sys_llphoton_refit_phi_elscaleup, sys_dijet_default[2]}, map_deltaphi)
      .Name("sys_llphoton_dijet_dphi_elscaleup").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_dijet_dphi_elscaledn = MultiMapNamedFunc(
      {sys_llphoton_refit_phi_elscaledn, sys_dijet_default[2]}, map_deltaphi)
      .Name("sys_llphoton_dijet_dphi_elscaledn").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_dijet_dphi_elresup = MultiMapNamedFunc(
      {sys_llphoton_refit_phi_elresup, sys_dijet_default[2]}, map_deltaphi)
      .Name("sys_llphoton_dijet_dphi_elresup").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_dijet_dphi_elresdn = MultiMapNamedFunc(
      {sys_llphoton_refit_phi_elresdn, sys_dijet_default[2]}, map_deltaphi)
      .Name("sys_llphoton_dijet_dphi_elresdn").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_dijet_dphi_muscaleup = MultiMapNamedFunc(
      {sys_llphoton_refit_phi_muscaleup, sys_dijet_default[2]}, map_deltaphi)
      .Name("sys_llphoton_dijet_dphi_muscaleup").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_dijet_dphi_muscaledn = MultiMapNamedFunc(
      {sys_llphoton_refit_phi_muscaledn, sys_dijet_default[2]}, map_deltaphi)
      .Name("sys_llphoton_dijet_dphi_muscaledn").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_dijet_dphi_muresup = MultiMapNamedFunc(
      {sys_llphoton_refit_phi_muresup, sys_dijet_default[2]}, map_deltaphi)
      .Name("sys_llphoton_dijet_dphi_muresup").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_dijet_dphi_muresdn = MultiMapNamedFunc(
      {sys_llphoton_refit_phi_muresdn, sys_dijet_default[2]}, map_deltaphi)
      .Name("sys_llphoton_dijet_dphi_muresdn").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_dijet_dphi_phscaleup = MultiMapNamedFunc(
      {sys_llphoton_refit_phi_phscaleup, sys_dijet_default[2]}, map_deltaphi)
      .Name("sys_llphoton_dijet_dphi_phscaleup").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_dijet_dphi_phscaledn = MultiMapNamedFunc(
      {sys_llphoton_refit_phi_phscaledn, sys_dijet_default[2]}, map_deltaphi)
      .Name("sys_llphoton_dijet_dphi_phscaledn").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_dijet_dphi_phresup = MultiMapNamedFunc(
      {sys_llphoton_refit_phi_phresup, sys_dijet_default[2]}, map_deltaphi)
      .Name("sys_llphoton_dijet_dphi_phresup").EnableCaching(true);
  const NamedFunc sys_llphoton_refit_dijet_dphi_phresdn = MultiMapNamedFunc(
      {sys_llphoton_refit_phi_phresdn, sys_dijet_default[2]}, map_deltaphi)
      .Name("sys_llphoton_dijet_dphi_phresdn").EnableCaching(true);
  vector<NamedFunc> sys_llphoton_refit_dijet_dphi_jetscaleup;
  vector<NamedFunc> sys_llphoton_refit_dijet_dphi_jetscaledn;
  vector<NamedFunc> sys_llphoton_refit_dijet_dphi_jetresup;
  vector<NamedFunc> sys_llphoton_refit_dijet_dphi_jetresdn;

  const NamedFunc sys_llphoton_refit_jet_dphi_default = 
      assign_variation_llphoton_refit_jet_dphi("llphoton_refit_phi", "jet_isgood",
      "jet_pt", "default");
  const NamedFunc sys_llphoton_refit_jet_dphi_default_elscaleup = 
      assign_variation_llphoton_refit_jet_dphi(
      sys_llphoton_refit_phi_elscaleup, "jet_isgood", "jet_pt", "default");
  const NamedFunc sys_llphoton_refit_jet_dphi_default_elscaledn = 
      assign_variation_llphoton_refit_jet_dphi(
      sys_llphoton_refit_phi_elscaledn, "jet_isgood", "jet_pt", "default");
  const NamedFunc sys_llphoton_refit_jet_dphi_default_elresup = 
      assign_variation_llphoton_refit_jet_dphi(
      sys_llphoton_refit_phi_elresup, "jet_isgood", "jet_pt", "default");
  const NamedFunc sys_llphoton_refit_jet_dphi_default_elresdn = 
      assign_variation_llphoton_refit_jet_dphi(
      sys_llphoton_refit_phi_elresdn, "jet_isgood", "jet_pt", "default");
  const NamedFunc sys_llphoton_refit_jet_dphi_default_muscaleup = 
      assign_variation_llphoton_refit_jet_dphi(
      sys_llphoton_refit_phi_muscaleup, "jet_isgood", "jet_pt", "default");
  const NamedFunc sys_llphoton_refit_jet_dphi_default_muscaledn = 
      assign_variation_llphoton_refit_jet_dphi(
      sys_llphoton_refit_phi_muscaledn, "jet_isgood", "jet_pt", "default");
  const NamedFunc sys_llphoton_refit_jet_dphi_default_muresup = 
      assign_variation_llphoton_refit_jet_dphi(
      sys_llphoton_refit_phi_muresup, "jet_isgood", "jet_pt", "default");
  const NamedFunc sys_llphoton_refit_jet_dphi_default_muresdn = 
      assign_variation_llphoton_refit_jet_dphi(
      sys_llphoton_refit_phi_muresdn, "jet_isgood", "jet_pt", "default");
  const NamedFunc sys_llphoton_refit_jet_dphi_default_phscaleup = 
      assign_variation_llphoton_refit_jet_dphi(
      sys_llphoton_refit_phi_phscaleup, "jet_isgood", "jet_pt", "default");
  const NamedFunc sys_llphoton_refit_jet_dphi_default_phscaledn = 
      assign_variation_llphoton_refit_jet_dphi(
      sys_llphoton_refit_phi_phscaledn, "jet_isgood", "jet_pt", "default");
  const NamedFunc sys_llphoton_refit_jet_dphi_default_phresup = 
      assign_variation_llphoton_refit_jet_dphi(
      sys_llphoton_refit_phi_phresup, "jet_isgood", "jet_pt", "default");
  const NamedFunc sys_llphoton_refit_jet_dphi_default_phresdn = 
      assign_variation_llphoton_refit_jet_dphi(
      sys_llphoton_refit_phi_phresdn, "jet_isgood", "jet_pt", "default");
  vector<NamedFunc> sys_llphoton_refit_jet_dphi_jetscaleup;
  vector<NamedFunc> sys_llphoton_refit_jet_dphi_jetscaledn;
  vector<NamedFunc> sys_llphoton_refit_jet_dphi_jetresup;
  vector<NamedFunc> sys_llphoton_refit_jet_dphi_jetresdn;

  const NamedFunc sys_llphoton_refit_jet_balance_default = 
      assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_default, 
      "photon_pt[0]", "photon_phi[0]", "jet_isgood", "jet_pt", "default");
  const NamedFunc sys_llphoton_refit_jet_balance_elscaleup = 
      assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_elscaleup, 
      "photon_pt[0]", "photon_phi[0]", "jet_isgood", "jet_pt", "elscaleup");
  const NamedFunc sys_llphoton_refit_jet_balance_elscaledn = 
      assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_elscaledn, 
      "photon_pt[0]", "photon_phi[0]", "jet_isgood", "jet_pt", "elscaledn");
  const NamedFunc sys_llphoton_refit_jet_balance_elresup = 
      assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_elresup, 
      "photon_pt[0]", "photon_phi[0]", "jet_isgood", "jet_pt", "elresup");
  const NamedFunc sys_llphoton_refit_jet_balance_elresdn = 
      assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_elresdn, 
      "photon_pt[0]", "photon_phi[0]", "jet_isgood", "jet_pt", "elresdn");
  const NamedFunc sys_llphoton_refit_jet_balance_muscaleup = 
      assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_muscaleup, 
      "photon_pt[0]", "photon_phi[0]", "jet_isgood", "jet_pt", "muscaleup");
  const NamedFunc sys_llphoton_refit_jet_balance_muscaledn = 
      assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_muscaledn, 
      "photon_pt[0]", "photon_phi[0]", "jet_isgood", "jet_pt", "muscaledn");
  const NamedFunc sys_llphoton_refit_jet_balance_muresup = 
      assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_muresup, 
      "photon_pt[0]", "photon_phi[0]", "jet_isgood", "jet_pt", "muresup");
  const NamedFunc sys_llphoton_refit_jet_balance_muresdn = 
      assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_muresdn, 
      "photon_pt[0]", "photon_phi[0]", "jet_isgood", "jet_pt", "muresdn");
  const NamedFunc sys_llphoton_refit_jet_balance_phscaleup = 
      assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_default, 
      sys_lead_photon_pt_scaleup, sys_lead_photon_phi_scaleup, "jet_isgood", 
      "jet_pt", "phscaleup");
  const NamedFunc sys_llphoton_refit_jet_balance_phscaledn = 
      assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_default, 
      sys_lead_photon_pt_scaledn, sys_lead_photon_phi_scaledn, "jet_isgood", 
      "jet_pt", "phscaledn");
  const NamedFunc sys_llphoton_refit_jet_balance_phresup = 
      assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_default, 
      sys_lead_photon_pt_resup, sys_lead_photon_phi_resup, "jet_isgood", 
      "jet_pt", "phresup");
  const NamedFunc sys_llphoton_refit_jet_balance_phresdn = 
      assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_default, 
      sys_lead_photon_pt_resdn, sys_lead_photon_phi_resdn, "jet_isgood", 
      "jet_pt", "phresdn");
  vector<NamedFunc> sys_llphoton_refit_jet_balance_jetscaleup;
  vector<NamedFunc> sys_llphoton_refit_jet_balance_jetscaledn;
  vector<NamedFunc> sys_llphoton_refit_jet_balance_jetresup;
  vector<NamedFunc> sys_llphoton_refit_jet_balance_jetresdn;

  const NamedFunc sys_mht_default = assign_variation_mht_balance("el_sig", 
      "el_pt", "mu_sig", "mu_pt", "photon_sig", "photon_pt", "jet_isgood",
      "jet_pt", "default");
  const NamedFunc sys_mht_elscaleup = assign_variation_mht_balance(
      sys_el_sig_scaleup, "sys_el_pt_scaleup", "mu_sig", "mu_pt", 
      "photon_sig", "photon_pt", "jet_isgood", "jet_pt", "elscaleup");
  const NamedFunc sys_mht_elscaledn = assign_variation_mht_balance(
      sys_el_sig_scaledn, "sys_el_pt_scaledn", "mu_sig", "mu_pt", 
      "photon_sig", "photon_pt", "jet_isgood", "jet_pt", "elscaledn");
  const NamedFunc sys_mht_elresup = assign_variation_mht_balance(
      sys_el_sig_resup, "sys_el_pt_resup", "mu_sig", "mu_pt", 
      "photon_sig", "photon_pt", "jet_isgood", "jet_pt", "elresup");
  const NamedFunc sys_mht_elresdn = assign_variation_mht_balance(
      sys_el_sig_resdn, "sys_el_pt_resdn", "mu_sig", "mu_pt", 
      "photon_sig", "photon_pt", "jet_isgood", "jet_pt", "elresdn");
  //TODO fix next 4 NamedFuncs in next production
  const NamedFunc sys_mht_muscaleup = assign_variation_mht_balance(
      "el_pt", "el_sig", sys_mu_sig_scaleup, "mu_pt", 
      "photon_sig", "photon_pt", "jet_isgood", "jet_pt", "muscaleup");
  const NamedFunc sys_mht_muscaledn = assign_variation_mht_balance(
      "el_pt", "el_sig", sys_mu_sig_scaledn, "mu_pt", 
      "photon_sig", "photon_pt", "jet_isgood", "jet_pt", "muscaledn");
  const NamedFunc sys_mht_muresup = assign_variation_mht_balance(
      "el_pt", "el_sig", sys_mu_sig_resup, "mu_pt+mu_ptErr", 
      "photon_sig", "photon_pt", "jet_isgood", "jet_pt", "muresup");
  const NamedFunc sys_mht_muresdn = assign_variation_mht_balance(
      "el_pt", "el_sig", sys_mu_sig_resdn, "mu_pt+mu_ptErr", 
      "photon_sig", "photon_pt", "jet_isgood", "jet_pt", "muresdn");
  const NamedFunc sys_mht_phscaleup = assign_variation_mht_balance(
      "el_pt", "el_sig", "mu_pt", "mu_sig", sys_photon_sig_scaleup, 
      "sys_photon_pt_scaleup", "jet_isgood", "jet_pt", "phscaleup");
  const NamedFunc sys_mht_phscaledn = assign_variation_mht_balance(
      "el_pt", "el_sig", "mu_pt", "mu_sig", sys_photon_sig_scaledn, 
      "sys_photon_pt_scaledn", "jet_isgood", "jet_pt", "phscaledn");
  const NamedFunc sys_mht_phresup = assign_variation_mht_balance(
      "el_pt", "el_sig", "mu_pt", "mu_sig", sys_photon_sig_resup, 
      "sys_photon_pt_resup", "jet_isgood", "jet_pt", "phresup");
  const NamedFunc sys_mht_phresdn = assign_variation_mht_balance(
      "el_pt", "el_sig", "mu_pt", "mu_sig", sys_photon_sig_resdn, 
      "sys_photon_pt_resdn", "jet_isgood", "jet_pt", "phresdn");
  vector<NamedFunc> sys_mht_jetscaleup;
  vector<NamedFunc> sys_mht_jetscaledn;
  vector<NamedFunc> sys_mht_jetresup;
  vector<NamedFunc> sys_mht_jetresdn;

  const NamedFunc sys_photon_mht_dphi_default = MultiMapNamedFunc(
      {"photon_phi[0]", sys_mht_default[1]}, map_deltaphi)
      .Name("sys_photon_mht_dphi_default").EnableCaching(true);
  const NamedFunc sys_photon_mht_dphi_elscaleup = MultiMapNamedFunc(
      {"photon_phi[0]", sys_mht_elscaleup[1]}, map_deltaphi)
      .Name("sys_photon_mht_dphi_elscaleup").EnableCaching(true);
  const NamedFunc sys_photon_mht_dphi_elscaledn = MultiMapNamedFunc(
      {"photon_phi[0]", sys_mht_elscaledn[1]}, map_deltaphi)
      .Name("sys_photon_mht_dphi_elscaledn").EnableCaching(true);
  const NamedFunc sys_photon_mht_dphi_elresup = MultiMapNamedFunc(
      {"photon_phi[0]", sys_mht_elresup[1]}, map_deltaphi)
      .Name("sys_photon_mht_dphi_elresup").EnableCaching(true);
  const NamedFunc sys_photon_mht_dphi_elresdn = MultiMapNamedFunc(
      {"photon_phi[0]", sys_mht_elresdn[1]}, map_deltaphi)
      .Name("sys_photon_mht_dphi_elresdn").EnableCaching(true);
  const NamedFunc sys_photon_mht_dphi_muscaleup = MultiMapNamedFunc(
      {"photon_phi[0]", sys_mht_muscaleup[1]}, map_deltaphi)
      .Name("sys_photon_mht_dphi_muscaleup").EnableCaching(true);
  const NamedFunc sys_photon_mht_dphi_muscaledn = MultiMapNamedFunc(
      {"photon_phi[0]", sys_mht_muscaledn[1]}, map_deltaphi)
      .Name("sys_photon_mht_dphi_muscaledn").EnableCaching(true);
  const NamedFunc sys_photon_mht_dphi_muresup = MultiMapNamedFunc(
      {"photon_phi[0]", sys_mht_muresup[1]}, map_deltaphi)
      .Name("sys_photon_mht_dphi_muresup").EnableCaching(true);
  const NamedFunc sys_photon_mht_dphi_muresdn = MultiMapNamedFunc(
      {"photon_phi[0]", sys_mht_muresdn[1]}, map_deltaphi)
      .Name("sys_photon_mht_dphi_muresdn").EnableCaching(true);
  const NamedFunc sys_photon_mht_dphi_phscaleup = MultiMapNamedFunc(
      {sys_lead_photon_phi_scaleup, sys_mht_phscaleup[1]}, map_deltaphi)
      .Name("sys_photon_mht_dphi_phscaleup").EnableCaching(true);
  const NamedFunc sys_photon_mht_dphi_phscaledn = MultiMapNamedFunc(
      {sys_lead_photon_phi_scaledn, sys_mht_phscaledn[1]}, map_deltaphi)
      .Name("sys_photon_mht_dphi_phscaledn").EnableCaching(true);
  const NamedFunc sys_photon_mht_dphi_phresup = MultiMapNamedFunc(
      {sys_lead_photon_phi_resup, sys_mht_phresup[1]}, map_deltaphi)
      .Name("sys_photon_mht_dphi_phresup").EnableCaching(true);
  const NamedFunc sys_photon_mht_dphi_phresdn = MultiMapNamedFunc(
      {sys_lead_photon_phi_resdn, sys_mht_phresdn[1]}, map_deltaphi)
      .Name("sys_photon_mht_dphi_phresdn").EnableCaching(true);
  vector<NamedFunc> sys_photon_mht_dphi_jetscaleup;
  vector<NamedFunc> sys_photon_mht_dphi_jetscaledn;
  vector<NamedFunc> sys_photon_mht_dphi_jetresup;
  vector<NamedFunc> sys_photon_mht_dphi_jetresdn;

  const NamedFunc sys_photon_jet1_dr_default = MultiMapNamedFunc(
      {"photon_eta[0]", "photon_phi[0]", lead_jet_eta, lead_jet_phi}, 
      map_deltar).Name("sys_photon_jet1_dr_default").EnableCaching(true);
  const NamedFunc sys_photon_jet1_dr_phscaleup = MultiMapNamedFunc(
      {sys_lead_photon_eta_scaleup, sys_lead_photon_phi_scaleup, lead_jet_eta, 
      lead_jet_phi}, map_deltar).Name("sys_photon_jet1_dr_phscaleup")
      .EnableCaching(true);
  const NamedFunc sys_photon_jet1_dr_phscaledn = MultiMapNamedFunc(
      {sys_lead_photon_eta_scaledn, sys_lead_photon_phi_scaledn, lead_jet_eta, 
      lead_jet_phi}, map_deltar).Name("sys_photon_jet1_dr_phscaledn")
      .EnableCaching(true);
  const NamedFunc sys_photon_jet1_dr_phresup = MultiMapNamedFunc(
      {sys_lead_photon_eta_resup, sys_lead_photon_phi_resup, lead_jet_eta, 
      lead_jet_phi}, map_deltar).Name("sys_photon_jet1_dr_phresup")
      .EnableCaching(true);
  const NamedFunc sys_photon_jet1_dr_phresdn = MultiMapNamedFunc(
      {sys_lead_photon_eta_resdn, sys_lead_photon_phi_resdn, lead_jet_eta, 
      lead_jet_phi}, map_deltar).Name("sys_photon_jet1_dr_phresdn")
      .EnableCaching(true);
  vector<NamedFunc> sys_photon_jet1_dr_jetscaleup;
  vector<NamedFunc> sys_photon_jet1_dr_jetscaledn;
  vector<NamedFunc> sys_photon_jet1_dr_jetresup;
  vector<NamedFunc> sys_photon_jet1_dr_jetresdn;

  const NamedFunc sys_photon_jet2_dr_default = MultiMapNamedFunc(
      {"photon_eta[0]", "photon_phi[0]", sublead_jet_eta, sublead_jet_phi}, 
      map_deltar).Name("sys_photon_jet2_dr_default").EnableCaching(true);
  const NamedFunc sys_photon_jet2_dr_phscaleup = MultiMapNamedFunc(
      {sys_lead_photon_eta_scaleup, sys_lead_photon_phi_scaleup, 
      sublead_jet_eta, sublead_jet_phi}, map_deltar).Name(
      "sys_photon_jet2_dr_phscaleup").EnableCaching(true);
  const NamedFunc sys_photon_jet2_dr_phscaledn = MultiMapNamedFunc(
      {sys_lead_photon_eta_scaledn, sys_lead_photon_phi_scaledn, 
      sublead_jet_eta, sublead_jet_phi}, map_deltar).Name(
      "sys_photon_jet2_dr_phscaledn").EnableCaching(true);
  const NamedFunc sys_photon_jet2_dr_phresup = MultiMapNamedFunc(
      {sys_lead_photon_eta_resup, sys_lead_photon_phi_resup, sublead_jet_eta, 
      sublead_jet_phi}, map_deltar).Name("sys_photon_jet2_dr_phresup")
      .EnableCaching(true);
  const NamedFunc sys_photon_jet2_dr_phresdn = MultiMapNamedFunc(
      {sys_lead_photon_eta_resdn, sys_lead_photon_phi_resdn, sublead_jet_eta, 
      sublead_jet_phi}, map_deltar).Name("sys_photon_jet2_dr_phresdn")
      .EnableCaching(true);
  vector<NamedFunc> sys_photon_jet2_dr_jetscaleup;
  vector<NamedFunc> sys_photon_jet2_dr_jetscaledn;
  vector<NamedFunc> sys_photon_jet2_dr_jetresup;
  vector<NamedFunc> sys_photon_jet2_dr_jetresdn;

  const NamedFunc sys_photon_zeppenfeld_default = 
      assign_variation_photon_zeppenfeld("photon_eta[0]","njet",lead_jet_eta,
      sublead_jet_eta, "default");
  const NamedFunc sys_photon_zeppenfeld_phscaleup = 
      assign_variation_photon_zeppenfeld(sys_lead_photon_eta_scaleup,"njet",
      lead_jet_eta, sublead_jet_eta, "phscaleup");
  const NamedFunc sys_photon_zeppenfeld_phscaledn = 
      assign_variation_photon_zeppenfeld(sys_lead_photon_eta_scaledn,"njet",
      lead_jet_eta, sublead_jet_eta, "phscaledn");
  const NamedFunc sys_photon_zeppenfeld_phresup = 
      assign_variation_photon_zeppenfeld(sys_lead_photon_eta_resup,"njet",
      lead_jet_eta, sublead_jet_eta, "phresup");
  const NamedFunc sys_photon_zeppenfeld_phresdn = 
      assign_variation_photon_zeppenfeld(sys_lead_photon_eta_resdn,"njet",
      lead_jet_eta, sublead_jet_eta, "phresdn");
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
    for (unsigned iyear = 0; iyear < years.size(); iyear++) {
      string year = years[iyear];
      sys_jet_pt_scaleup.push_back(assign_variation_jet_pt_year(
          "sys_jet_pt_jesup", year, "scaleup"));
      sys_jet_pt_scaledn.push_back(assign_variation_jet_pt_year(
          "sys_jet_pt_jesdn", year, "scaledn"));
      sys_jet_pt_resup.push_back(assign_variation_jet_pt_year(
          "sys_jet_pt_jerup", year, "resup"));
      sys_jet_pt_resdn.push_back(assign_variation_jet_pt_year(
          "sys_jet_pt_jerdn", year, "resdn"));
      //TODO uncomment real ones for next production
      sys_jet_isgood_scaleup.push_back(NamedFunc(
          pinnacles_run2_jet_isgood_nopt
          &&sys_jet_pt_scaleup[iyear]>30.0).Name("sys_jet_isgood_scaleup"+year)
          .EnableCaching(true));
      sys_jet_isgood_scaledn.push_back(NamedFunc(
          pinnacles_run2_jet_isgood_nopt
          &&sys_jet_pt_scaledn[iyear]>30.0).Name("sys_jet_isgood_scaledn"+year)
          .EnableCaching(true));
      sys_jet_isgood_resup.push_back(NamedFunc(
          pinnacles_run2_jet_isgood_nopt
          &&sys_jet_pt_resup[iyear]>30.0).Name("sys_jet_isgood_resup"+year)
          .EnableCaching(true));
      sys_jet_isgood_resdn.push_back(NamedFunc(
          pinnacles_run2_jet_isgood_nopt
          &&sys_jet_pt_resdn[iyear]>30.0).Name("sys_jet_isgood_resdn"+year)
          .EnableCaching(true));
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
      sys_sig_jet_pt_scaleup.push_back(FilterNamedFunc(sys_jet_pt_scaleup[iyear],
          sys_jet_isgood_scaleup[iyear]).Name("sys_sig_jet_pt_scaleup"+year)
          .EnableCaching(true));
      sys_sig_jet_pt_scaledn.push_back(FilterNamedFunc(sys_jet_pt_scaledn[iyear],
          sys_jet_isgood_scaledn[iyear]).Name("sys_sig_jet_pt_scaledn"+year)
          .EnableCaching(true));
      sys_sig_jet_pt_resup.push_back(FilterNamedFunc(sys_jet_pt_resup[iyear],
          sys_jet_isgood_resup[iyear]).Name("sys_sig_jet_pt_resup"+year)
          .EnableCaching(true));
      sys_sig_jet_pt_resdn.push_back(FilterNamedFunc(sys_jet_pt_resdn[iyear],
          sys_jet_isgood_resdn[iyear]).Name("sys_sig_jet_pt_resdn"+year)
          .EnableCaching(true));
      sys_sig_jet_eta_scaleup.push_back(FilterNamedFunc("jet_eta",
          sys_jet_isgood_scaleup[iyear]).Name("sys_sig_jet_eta_scaleup"+year)
          .EnableCaching(true));
      sys_sig_jet_eta_scaledn.push_back(FilterNamedFunc("jet_eta",
          sys_jet_isgood_scaledn[iyear]).Name("sys_sig_jet_eta_scaledn"+year)
          .EnableCaching(true));
      sys_sig_jet_eta_resup.push_back(FilterNamedFunc("jet_eta",
          sys_jet_isgood_resup[iyear]).Name("sys_sig_jet_eta_resup"+year)
          .EnableCaching(true));
      sys_sig_jet_eta_resdn.push_back(FilterNamedFunc("jet_eta",
          sys_jet_isgood_resdn[iyear]).Name("sys_sig_jet_eta_resdn"+year)
          .EnableCaching(true));
      sys_sig_jet_phi_scaleup.push_back(FilterNamedFunc("jet_phi",
          sys_jet_isgood_scaleup[iyear]).Name("sys_sig_jet_phi_scaleup"+year)
          .EnableCaching(true));
      sys_sig_jet_phi_scaledn.push_back(FilterNamedFunc("jet_phi",
          sys_jet_isgood_scaledn[iyear]).Name("sys_sig_jet_phi_scaledn"+year)
          .EnableCaching(true));
      sys_sig_jet_phi_resup.push_back(FilterNamedFunc("jet_phi",
          sys_jet_isgood_resup[iyear]).Name("sys_sig_jet_phi_resup"+year)
          .EnableCaching(true));
      sys_sig_jet_phi_resdn.push_back(FilterNamedFunc("jet_phi",
          sys_jet_isgood_resdn[iyear]).Name("sys_sig_jet_phi_resdn"+year)
          .EnableCaching(true));
      sys_sig_jet_m_scaleup.push_back(FilterNamedFunc("jet_m",
          sys_jet_isgood_scaleup[iyear]).Name("sys_sig_jet_m_scaleup"+year)
          .EnableCaching(true));
      sys_sig_jet_m_scaledn.push_back(FilterNamedFunc("jet_m",
          sys_jet_isgood_scaledn[iyear]).Name("sys_sig_jet_m_scaledn"+year)
          .EnableCaching(true));
      sys_sig_jet_m_resup.push_back(FilterNamedFunc("jet_m",
          sys_jet_isgood_resup[iyear]).Name("sys_sig_jet_m_resup"+year)
          .EnableCaching(true));
      sys_sig_jet_m_resdn.push_back(FilterNamedFunc("jet_m",
          sys_jet_isgood_resdn[iyear]).Name("sys_sig_jet_m_resdn"+year)
          .EnableCaching(true));
      sys_lead_jet_pt_scaleup.push_back(ReduceNamedFunc(
          sys_sig_jet_pt_scaleup[iyear],reduce_max).Name(
          "sys_lead_jet_pt_scaleup"+year).EnableCaching(true));
      sys_lead_jet_pt_scaledn.push_back(ReduceNamedFunc(
          sys_sig_jet_pt_scaledn[iyear],reduce_max).Name(
          "sys_lead_jet_pt_scaledn"+year).EnableCaching(true));
      sys_lead_jet_pt_resup.push_back(ReduceNamedFunc(
          sys_sig_jet_pt_resup[iyear],reduce_max).Name(
          "sys_lead_jet_pt_resup"+year).EnableCaching(true));
      sys_lead_jet_pt_resdn.push_back(ReduceNamedFunc(
          sys_sig_jet_pt_resdn[iyear],reduce_max).Name(
          "sys_lead_jet_pt_resdn"+year).EnableCaching(true));
      sys_sublead_jet_pt_scaleup.push_back(ReduceNamedFunc(
          sys_sig_jet_pt_scaleup[iyear],reduce_sublead).Name(
          "sys_sublead_jet_pt_scaleup"+year).EnableCaching(true));
      sys_sublead_jet_pt_scaledn.push_back(ReduceNamedFunc(
          sys_sig_jet_pt_scaledn[iyear],reduce_sublead).Name(
          "sys_sublead_jet_pt_scaledn"+year).EnableCaching(true));
      sys_sublead_jet_pt_resup.push_back(ReduceNamedFunc(
          sys_sig_jet_pt_resup[iyear],reduce_sublead).Name(
          "sys_sublead_jet_pt_resup"+year).EnableCaching(true));
      sys_sublead_jet_pt_resdn.push_back(ReduceNamedFunc(
          sys_sig_jet_pt_resdn[iyear],reduce_sublead).Name(
          "sys_sublead_jet_pt_resdn"+year).EnableCaching(true));
      sys_lead_jet_eta_scaleup.push_back(MultiReduceNamedFunc(
          {sys_sig_jet_pt_scaleup[iyear],sys_sig_jet_eta_scaleup[iyear]}, 
          reduce_maxfirst).Name("sys_lead_jet_eta_scaleup"+year)
          .EnableCaching(true));
      sys_lead_jet_eta_scaledn.push_back(MultiReduceNamedFunc(
          {sys_sig_jet_pt_scaledn[iyear],sys_sig_jet_eta_scaledn[iyear]}, 
          reduce_maxfirst).Name("sys_lead_jet_eta_scaledn"+year)
          .EnableCaching(true));
      sys_lead_jet_eta_resup.push_back(MultiReduceNamedFunc(
          {sys_sig_jet_pt_resup[iyear],sys_sig_jet_eta_resup[iyear]}, 
          reduce_maxfirst).Name("sys_lead_jet_eta_resup"+year)
          .EnableCaching(true));
      sys_lead_jet_eta_resdn.push_back(MultiReduceNamedFunc(
          {sys_sig_jet_pt_resdn[iyear],sys_sig_jet_eta_resdn[iyear]}, 
          reduce_maxfirst).Name("sys_lead_jet_eta_resdn"+year)
          .EnableCaching(true));
      sys_lead_jet_phi_scaleup.push_back(MultiReduceNamedFunc(
          {sys_sig_jet_pt_scaleup[iyear],sys_sig_jet_phi_scaleup[iyear]}, 
          reduce_maxfirst).Name("sys_lead_jet_phi_scaleup"+year)
          .EnableCaching(true));
      sys_lead_jet_phi_scaledn.push_back(MultiReduceNamedFunc(
          {sys_sig_jet_pt_scaledn[iyear],sys_sig_jet_phi_scaledn[iyear]}, 
          reduce_maxfirst).Name("sys_lead_jet_phi_scaledn"+year)
          .EnableCaching(true));
      sys_lead_jet_phi_resup.push_back(MultiReduceNamedFunc(
          {sys_sig_jet_pt_resup[iyear],sys_sig_jet_phi_resup[iyear]}, 
          reduce_maxfirst).Name("sys_lead_jet_phi_resup"+year)
          .EnableCaching(true));
      sys_lead_jet_phi_resdn.push_back(MultiReduceNamedFunc(
          {sys_sig_jet_pt_resdn[iyear],sys_sig_jet_phi_resdn[iyear]}, 
          reduce_maxfirst).Name("sys_lead_jet_phi_resdn"+year)
          .EnableCaching(true));
      sys_lead_jet_m_scaleup.push_back(MultiReduceNamedFunc(
          {sys_sig_jet_pt_scaleup[iyear],sys_sig_jet_m_scaleup[iyear]}, 
          reduce_maxfirst).Name("sys_lead_jet_m_scaleup"+year)
          .EnableCaching(true));
      sys_lead_jet_m_scaledn.push_back(MultiReduceNamedFunc(
          {sys_sig_jet_pt_scaledn[iyear],sys_sig_jet_m_scaledn[iyear]}, 
          reduce_maxfirst).Name("sys_lead_jet_m_scaledn"+year)
          .EnableCaching(true));
      sys_lead_jet_m_resup.push_back(MultiReduceNamedFunc(
          {sys_sig_jet_pt_resup[iyear],sys_sig_jet_m_resup[iyear]}, 
          reduce_maxfirst).Name("sys_lead_jet_m_resup"+year)
          .EnableCaching(true));
      sys_lead_jet_m_resdn.push_back(MultiReduceNamedFunc(
          {sys_sig_jet_pt_resdn[iyear],sys_sig_jet_m_resdn[iyear]}, 
          reduce_maxfirst).Name("sys_lead_jet_m_resdn"+year)
          .EnableCaching(true));
      sys_sublead_jet_eta_scaleup.push_back(MultiReduceNamedFunc(
          {sys_sig_jet_pt_scaleup[iyear],sys_sig_jet_eta_scaleup[iyear]}, 
          reduce_subleadfirst).Name("sys_sublead_jet_eta_scaleup"+year)
          .EnableCaching(true));
      sys_sublead_jet_eta_scaledn.push_back(MultiReduceNamedFunc(
          {sys_sig_jet_pt_scaledn[iyear],sys_sig_jet_eta_scaledn[iyear]}, 
          reduce_subleadfirst).Name("sys_sublead_jet_eta_scaledn"+year)
          .EnableCaching(true));
      sys_sublead_jet_eta_resup.push_back(MultiReduceNamedFunc(
          {sys_sig_jet_pt_resup[iyear],sys_sig_jet_eta_resup[iyear]}, 
          reduce_subleadfirst).Name("sys_sublead_jet_eta_resup"+year)
          .EnableCaching(true));
      sys_sublead_jet_eta_resdn.push_back(MultiReduceNamedFunc(
          {sys_sig_jet_pt_resdn[iyear],sys_sig_jet_eta_resdn[iyear]}, 
          reduce_subleadfirst).Name("sys_sublead_jet_eta_resdn"+year)
          .EnableCaching(true));
      sys_sublead_jet_phi_scaleup.push_back(MultiReduceNamedFunc(
          {sys_sig_jet_pt_scaleup[iyear],sys_sig_jet_phi_scaleup[iyear]}, 
          reduce_subleadfirst).Name("sys_sublead_jet_phi_scaleup"+year)
          .EnableCaching(true));
      sys_sublead_jet_phi_scaledn.push_back(MultiReduceNamedFunc(
          {sys_sig_jet_pt_scaledn[iyear],sys_sig_jet_phi_scaledn[iyear]}, 
          reduce_subleadfirst).Name("sys_sublead_jet_phi_scaledn"+year)
          .EnableCaching(true));
      sys_sublead_jet_phi_resup.push_back(MultiReduceNamedFunc(
          {sys_sig_jet_pt_resup[iyear],sys_sig_jet_phi_resup[iyear]}, 
          reduce_subleadfirst).Name("sys_sublead_jet_phi_resup"+year)
          .EnableCaching(true));
      sys_sublead_jet_phi_resdn.push_back(MultiReduceNamedFunc(
          {sys_sig_jet_pt_resdn[iyear],sys_sig_jet_phi_resdn[iyear]}, 
          reduce_subleadfirst).Name("sys_sublead_jet_phi_resdn"+year)
          .EnableCaching(true));
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
      sys_njet_scaleup.push_back(ReduceNamedFunc(sys_jet_isgood_scaleup[iyear],
          reduce_sum).Name("sys_njet_scaleup"+year).EnableCaching(true));
      sys_njet_scaledn.push_back(ReduceNamedFunc(sys_jet_isgood_scaledn[iyear],
          reduce_sum).Name("sys_njet_scaledn"+year).EnableCaching(true));
      sys_njet_resup.push_back(ReduceNamedFunc(sys_jet_isgood_resup[iyear],
          reduce_sum).Name("sys_njet_resup"+year).EnableCaching(true));
      sys_njet_resdn.push_back(ReduceNamedFunc(sys_jet_isgood_resdn[iyear],
          reduce_sum).Name("sys_njet_resdn"+year).EnableCaching(true));
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
      untagged_category_jetscaleup.push_back(assign_variation_untagged_category(
          "nlep", sys_njet_scaleup[iyear], sys_nbdfm_scaleup[iyear], 
          sys_met_scaleup[iyear], "llphoton_pt[0]", "llphoton_m[0]", 
          max_lep_miniso, "ll_m[0]", "jetscaleup"));
      untagged_category_jetscaledn.push_back(assign_variation_untagged_category(
          "nlep", sys_njet_scaledn[iyear], sys_nbdfm_scaledn[iyear], 
          sys_met_scaledn[iyear], "llphoton_pt[0]", "llphoton_m[0]", 
          max_lep_miniso, "ll_m[0]", "jetscaledn"));
      untagged_category_jetresup.push_back(assign_variation_untagged_category(
          "nlep", sys_njet_resup[iyear], sys_nbdfm_resup[iyear], 
          sys_met_resup[iyear], "llphoton_pt[0]", "llphoton_m[0]", 
          max_lep_miniso, "ll_m[0]", "jetresup"));
      untagged_category_jetresdn.push_back(assign_variation_untagged_category(
          "nlep", sys_njet_resdn[iyear], sys_nbdfm_resdn[iyear], 
          sys_met_resdn[iyear], "llphoton_pt[0]", "llphoton_m[0]", 
          max_lep_miniso, "ll_m[0]", "jetresdn"));
      sys_dijet_scaleup.push_back(assign_variation_dijet(
          sys_jet_isgood_scaleup[iyear], "sys_jet_pt_jesup", "sys_jet_m_jesup",
          "scaleup"));
      sys_dijet_scaledn.push_back(assign_variation_dijet(
          sys_jet_isgood_scaledn[iyear], "sys_jet_pt_jesdn", "sys_jet_m_jesdn",
          "scaledn"));
      sys_dijet_resup.push_back(assign_variation_dijet(
          sys_jet_isgood_resup[iyear], "sys_jet_pt_jerup", "sys_jet_m_jerup",
          "resup"));
      sys_dijet_resdn.push_back(assign_variation_dijet(
          sys_jet_isgood_resdn[iyear], "sys_jet_pt_jerdn", "sys_jet_m_jerdn",
          "resdn"));
      sys_dijet_m_scaleup.push_back(NamedFunc(sys_dijet_scaleup[iyear][3])
          .Name("sys_dijet_m_scaleup").EnableCaching(true));
      sys_dijet_m_scaledn.push_back(NamedFunc(sys_dijet_scaledn[iyear][3])
          .Name("sys_dijet_m_scaledn").EnableCaching(true));
      sys_dijet_m_resup.push_back(NamedFunc(sys_dijet_resup[iyear][3])
          .Name("sys_dijet_m_resup").EnableCaching(true));
      sys_dijet_m_resdn.push_back(NamedFunc(sys_dijet_resdn[iyear][3])
          .Name("sys_dijet_m_resdn").EnableCaching(true));
      sys_dijet_deta_scaleup.push_back(NamedFunc(sys_dijet_scaleup[iyear][4])
          .Name("sys_dijet_deta_scaleup").EnableCaching(true));
      sys_dijet_deta_scaledn.push_back(NamedFunc(sys_dijet_scaledn[iyear][4])
          .Name("sys_dijet_deta_scaledn").EnableCaching(true));
      sys_dijet_deta_resup.push_back(NamedFunc(sys_dijet_resup[iyear][4])
          .Name("sys_dijet_deta_resup").EnableCaching(true));
      sys_dijet_deta_resdn.push_back(NamedFunc(sys_dijet_resdn[iyear][4])
          .Name("sys_dijet_deta_resdn").EnableCaching(true));
      sys_dijet_dphi_scaleup.push_back(NamedFunc(sys_dijet_scaleup[iyear][5])
          .Name("sys_dijet_dphi_scaleup").EnableCaching(true));
      sys_dijet_dphi_scaledn.push_back(NamedFunc(sys_dijet_scaledn[iyear][5])
          .Name("sys_dijet_dphi_scaledn").EnableCaching(true));
      sys_dijet_dphi_resup.push_back(NamedFunc(sys_dijet_resup[iyear][5])
          .Name("sys_dijet_dphi_resup").EnableCaching(true));
      sys_dijet_dphi_resdn.push_back(NamedFunc(sys_dijet_resdn[iyear][5])
          .Name("sys_dijet_dphi_resdn").EnableCaching(true));
      sys_llphoton_refit_dijet_dphi_jetscaleup.push_back(MultiMapNamedFunc(
          {"llphoton_refit_phi", sys_dijet_scaleup[iyear][2]}, map_deltaphi)
          .Name("sys_llphoton_dijet_dphi_jetscaleup").EnableCaching(true));
      sys_llphoton_refit_dijet_dphi_jetscaledn.push_back(MultiMapNamedFunc(
          {"llphoton_refit_phi", sys_dijet_scaledn[iyear][2]}, map_deltaphi)
          .Name("sys_llphoton_dijet_dphi_jetscaledn").EnableCaching(true));
      sys_llphoton_refit_dijet_dphi_jetresup.push_back(MultiMapNamedFunc(
          {"llphoton_refit_phi", sys_dijet_resup[iyear][2]}, map_deltaphi)
          .Name("sys_llphoton_dijet_dphi_jetresup").EnableCaching(true));
      sys_llphoton_refit_dijet_dphi_jetresdn.push_back(MultiMapNamedFunc(
          {"llphoton_refit_phi", sys_dijet_resdn[iyear][2]}, map_deltaphi)
          .Name("sys_llphoton_dijet_dphi_jetresdn").EnableCaching(true));
      sys_llphoton_refit_jet_dphi_jetscaleup.push_back(
          assign_variation_llphoton_refit_jet_dphi("llphoton_refit_phi", 
          sys_jet_isgood_scaleup[iyear], "sys_jet_pt_jesup", "jetscaleup"));
      sys_llphoton_refit_jet_dphi_jetscaledn.push_back(
          assign_variation_llphoton_refit_jet_dphi("llphoton_refit_phi", 
          sys_jet_isgood_scaledn[iyear], "sys_jet_pt_jesdn", "jetscaledn"));
      sys_llphoton_refit_jet_dphi_jetresup.push_back(
          assign_variation_llphoton_refit_jet_dphi("llphoton_refit_phi", 
          sys_jet_isgood_resup[iyear], "sys_jet_pt_jerup", "jetresup"));
      sys_llphoton_refit_jet_dphi_jetresdn.push_back(
          assign_variation_llphoton_refit_jet_dphi("llphoton_refit_phi", 
          sys_jet_isgood_resdn[iyear], "sys_jet_pt_jerdn", "jetresdn"));
      sys_llphoton_refit_jet_balance_jetscaleup.push_back(
          assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_default, 
          "photon_pt[0]", "photon_phi[0]", sys_jet_isgood_scaleup[iyear], 
          "sys_jet_pt_jesup", "jetscaleup"));
      sys_llphoton_refit_jet_balance_jetscaledn.push_back(
          assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_default, 
          "photon_pt[0]", "photon_phi[0]", sys_jet_isgood_scaledn[iyear], 
          "sys_jet_pt_jesdn", "jetscaledn"));
      sys_llphoton_refit_jet_balance_jetresup.push_back(
          assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_default, 
          "photon_pt[0]", "photon_phi[0]", sys_jet_isgood_resup[iyear], 
          "sys_jet_pt_jerup", "jetresup"));
      sys_llphoton_refit_jet_balance_jetresdn.push_back(
          assign_variation_llphoton_refit_jet_balance(sys_ll_refit_p4_default, 
          "photon_pt[0]", "photon_phi[0]", sys_jet_isgood_resdn[iyear], 
          "sys_jet_pt_jerdn", "jetresdn"));
      sys_mht_jetscaleup.push_back(assign_variation_mht_balance("el_sig", 
          "el_pt", "mu_sig", "mu_pt", "photon_sig", "photon_pt", 
          sys_jet_isgood_scaleup[iyear], "sys_jet_pt_jesup", "jetscaleup"));
      sys_mht_jetscaledn.push_back(assign_variation_mht_balance("el_sig", 
          "el_pt", "mu_sig", "mu_pt", "photon_sig", "photon_pt", 
          sys_jet_isgood_scaledn[iyear], "sys_jet_pt_jesdn", "jetscaledn"));
      sys_mht_jetresup.push_back(assign_variation_mht_balance("el_sig", 
          "el_pt", "mu_sig", "mu_pt", "photon_sig", "photon_pt", 
          sys_jet_isgood_resup[iyear], "sys_jet_pt_jerup", "jetresup"));
      sys_mht_jetresdn.push_back(assign_variation_mht_balance("el_sig", 
          "el_pt", "mu_sig", "mu_pt", "photon_sig", "photon_pt", 
          sys_jet_isgood_resdn[iyear], "sys_jet_pt_jerdn", "jetresdn"));
      sys_photon_mht_dphi_jetscaleup.push_back(MultiMapNamedFunc(
          {"photon_phi[0]", sys_mht_jetscaleup[iyear][1]}, map_deltaphi)
          .Name("sys_photon_mht_dphi_jetscaleup").EnableCaching(true));
      sys_photon_mht_dphi_jetscaledn.push_back(MultiMapNamedFunc(
          {"photon_phi[0]", sys_mht_jetscaledn[iyear][1]}, map_deltaphi)
          .Name("sys_photon_mht_dphi_jetscaledn").EnableCaching(true));
      sys_photon_mht_dphi_jetresup.push_back(MultiMapNamedFunc(
          {"photon_phi[0]", sys_mht_jetresup[iyear][1]}, map_deltaphi)
          .Name("sys_photon_mht_dphi_jetresup").EnableCaching(true));
      sys_photon_mht_dphi_jetresdn.push_back(MultiMapNamedFunc(
          {"photon_phi[0]", sys_mht_jetresdn[iyear][1]}, map_deltaphi)
          .Name("sys_photon_mht_dphi_jetresdn").EnableCaching(true));
      sys_photon_jet1_dr_jetscaleup.push_back(MultiMapNamedFunc(
          {"photon_eta[0]", "photon_phi[0]", sys_lead_jet_eta_scaleup[iyear], 
          sys_lead_jet_phi_scaleup[iyear]}, map_deltar).Name(
          "sys_photon_jet1_dr_jetscaleup").EnableCaching(true));
      sys_photon_jet1_dr_jetscaledn.push_back(MultiMapNamedFunc(
          {"photon_eta[0]", "photon_phi[0]", sys_lead_jet_eta_scaledn[iyear], 
          sys_lead_jet_phi_scaledn[iyear]}, map_deltar).Name(
          "sys_photon_jet1_dr_jetscaledn").EnableCaching(true));
      sys_photon_jet1_dr_jetresup.push_back(MultiMapNamedFunc(
          {"photon_eta[0]", "photon_phi[0]", sys_lead_jet_eta_resup[iyear], 
          sys_lead_jet_phi_resup[iyear]}, map_deltar).Name(
          "sys_photon_jet1_dr_jetresup").EnableCaching(true));
      sys_photon_jet1_dr_jetresdn.push_back(MultiMapNamedFunc(
          {"photon_eta[0]", "photon_phi[0]", sys_lead_jet_eta_resdn[iyear], 
          sys_lead_jet_phi_resdn[iyear]}, map_deltar).Name(
          "sys_photon_jet1_dr_jetresdn").EnableCaching(true));
      sys_photon_jet2_dr_jetscaleup.push_back(MultiMapNamedFunc(
          {"photon_eta[0]", "photon_phi[0]", 
          sys_sublead_jet_eta_scaleup[iyear], 
          sys_sublead_jet_phi_scaleup[iyear]}, map_deltar).Name(
          "sys_photon_jet2_dr_jetscaleup").EnableCaching(true));
      sys_photon_jet2_dr_jetscaledn.push_back(MultiMapNamedFunc(
          {"photon_eta[0]", "photon_phi[0]", 
          sys_sublead_jet_eta_scaledn[iyear], 
          sys_sublead_jet_phi_scaledn[iyear]}, map_deltar).Name(
          "sys_photon_jet2_dr_jetscaledn").EnableCaching(true));
      sys_photon_jet2_dr_jetresup.push_back(MultiMapNamedFunc(
          {"photon_eta[0]", "photon_phi[0]", 
          sys_sublead_jet_eta_resup[iyear], 
          sys_sublead_jet_phi_resup[iyear]}, map_deltar).Name(
          "sys_photon_jet2_dr_jetresup").EnableCaching(true));
      sys_photon_jet2_dr_jetresdn.push_back(MultiMapNamedFunc(
          {"photon_eta[0]", "photon_phi[0]", 
          sys_sublead_jet_eta_resdn[iyear], 
          sys_sublead_jet_phi_resdn[iyear]}, map_deltar).Name(
          "sys_photon_jet2_dr_jetresdn").EnableCaching(true));
      sys_photon_zeppenfeld_jetscaleup.push_back(
          assign_variation_photon_zeppenfeld("photon_eta[0]",
          sys_njet_scaleup[iyear],sys_lead_jet_eta_scaleup[iyear], 
          sys_sublead_jet_eta_scaleup[iyear],"jetscaleup"));
      sys_photon_zeppenfeld_jetscaledn.push_back(
          assign_variation_photon_zeppenfeld("photon_eta[0]",
          sys_njet_scaledn[iyear], sys_lead_jet_eta_scaledn[iyear], 
          sys_sublead_jet_eta_scaledn[iyear], "jetscaledn"));
      sys_photon_zeppenfeld_jetresup.push_back(
          assign_variation_photon_zeppenfeld("photon_eta[0]",
          sys_njet_resup[iyear], sys_lead_jet_eta_resup[iyear], 
          sys_sublead_jet_eta_resup[iyear], "jetresup"));
      sys_photon_zeppenfeld_jetresdn.push_back(
          assign_variation_photon_zeppenfeld("photon_eta[0]",
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
    for (std::shared_ptr<MVAWrapper> bdt_reader : vbf_bdt_readers) {
      //default
      bdt_reader->SetVariable("dijet_deta",sys_dijet_deta_default);
      bdt_reader->SetVariable("dijet_dphi",sys_dijet_dphi_default);
      bdt_reader->SetVariable("j1_pt",lead_jet_pt);
      bdt_reader->SetVariable("j2_pt",sublead_jet_pt);
      bdt_reader->SetVariable("llphoton_dijet_balance",
                              sys_llphoton_refit_jet_balance_default);
      bdt_reader->SetVariable("llphoton_dijet_dphi",
                              sys_llphoton_refit_dijet_dphi_default);
      bdt_reader->SetVariable("phi",sys_llphoton_refit_psi_default);
      bdt_reader->SetVariable("min_dR","photon_drmin[0]");
      bdt_reader->SetVariable("max_dR","photon_drmax[0]");
      bdt_reader->SetVariable("costheta",sys_llphoton_refit_costheta_default);
      bdt_reader->SetVariable("cosTheta",sys_llphoton_refit_cosTheta_default);
      bdt_reader->SetVariable("pt_mass","llphoton_refit_pt/llphoton_refit_m");
      bdt_reader->SetVariable("l1_rapidity",sys_lead_lepton_eta_default);
      bdt_reader->SetVariable("l2_rapidity",sys_sublead_lepton_eta_default);
      bdt_reader->SetVariable("photon_rapidity","photon_eta[0]");
      bdt_reader->SetVariable("photon_mva","photon_idmva[0]");
      bdt_reader->SetVariable("photon_res",photon_relpterr);
      bdt_reader->SetVariable("photon_zeppenfeld",
                              sys_photon_zeppenfeld_default);
      bdt_reader->SetVariable("photon_jet1_dr",sys_photon_jet1_dr_default);
      bdt_reader->SetVariable("photon_jet2_dr",sys_photon_jet2_dr_default);
      //elscaleup
      bdt_reader->SetAltVariable("elscaleup",sys_dijet_deta_default);
      bdt_reader->SetAltVariable("elscaleup",sys_dijet_dphi_default);
      bdt_reader->SetAltVariable("elscaleup",lead_jet_pt);
      bdt_reader->SetAltVariable("elscaleup",sublead_jet_pt);
      bdt_reader->SetAltVariable("elscaleup",
                                 sys_llphoton_refit_jet_balance_elscaleup);
      bdt_reader->SetAltVariable("elscaleup",
                                 sys_llphoton_refit_dijet_dphi_elscaleup);
      bdt_reader->SetAltVariable("elscaleup",sys_llphoton_refit_psi_elscaleup);
      bdt_reader->SetAltVariable("elscaleup","photon_drmin[0]");
      bdt_reader->SetAltVariable("elscaleup","photon_drmax[0]");
      bdt_reader->SetAltVariable("elscaleup",
                                 sys_llphoton_refit_costheta_elscaleup);
      bdt_reader->SetAltVariable("elscaleup",
                                 sys_llphoton_refit_cosTheta_elscaleup);
      bdt_reader->SetAltVariable("elscaleup",
                                 sys_llphoton_refit_pt_elscaleup
                                 /sys_llphoton_refit_m_elscaleup);
      bdt_reader->SetAltVariable("elscaleup",sys_lead_lepton_eta_elscaleup);
      bdt_reader->SetAltVariable("elscaleup",sys_sublead_lepton_eta_elscaleup);
      bdt_reader->SetAltVariable("elscaleup","photon_eta[0]");
      bdt_reader->SetAltVariable("elscaleup","photon_idmva[0]");
      bdt_reader->SetAltVariable("elscaleup",photon_relpterr);
      bdt_reader->SetAltVariable("elscaleup",
                                 sys_photon_zeppenfeld_default);
      bdt_reader->SetAltVariable("elscaleup",sys_photon_jet1_dr_default);
      bdt_reader->SetAltVariable("elscaleup",sys_photon_jet2_dr_default);
      //elscaledn
      bdt_reader->SetAltVariable("elscaledn",sys_dijet_deta_default);
      bdt_reader->SetAltVariable("elscaledn",sys_dijet_dphi_default);
      bdt_reader->SetAltVariable("elscaledn",lead_jet_pt);
      bdt_reader->SetAltVariable("elscaledn",sublead_jet_pt);
      bdt_reader->SetAltVariable("elscaledn",
                                 sys_llphoton_refit_jet_balance_elscaledn);
      bdt_reader->SetAltVariable("elscaledn",
                                 sys_llphoton_refit_dijet_dphi_elscaledn);
      bdt_reader->SetAltVariable("elscaledn",sys_llphoton_refit_psi_elscaledn);
      bdt_reader->SetAltVariable("elscaledn","photon_drmin[0]");
      bdt_reader->SetAltVariable("elscaledn","photon_drmax[0]");
      bdt_reader->SetAltVariable("elscaledn",
                                 sys_llphoton_refit_costheta_elscaledn);
      bdt_reader->SetAltVariable("elscaledn",
                                 sys_llphoton_refit_cosTheta_elscaledn);
      bdt_reader->SetAltVariable("elscaledn",
                                 sys_llphoton_refit_pt_elscaledn
                                 /sys_llphoton_refit_m_elscaledn);
      bdt_reader->SetAltVariable("elscaledn",sys_lead_lepton_eta_elscaledn);
      bdt_reader->SetAltVariable("elscaledn",sys_sublead_lepton_eta_elscaledn);
      bdt_reader->SetAltVariable("elscaledn","photon_eta[0]");
      bdt_reader->SetAltVariable("elscaledn","photon_idmva[0]");
      bdt_reader->SetAltVariable("elscaledn",photon_relpterr);
      bdt_reader->SetAltVariable("elscaledn",
                                 sys_photon_zeppenfeld_default);
      bdt_reader->SetAltVariable("elscaledn",sys_photon_jet1_dr_default);
      bdt_reader->SetAltVariable("elscaledn",sys_photon_jet2_dr_default);
      //elresup
      bdt_reader->SetAltVariable("elresup",sys_dijet_deta_default);
      bdt_reader->SetAltVariable("elresup",sys_dijet_dphi_default);
      bdt_reader->SetAltVariable("elresup",lead_jet_pt);
      bdt_reader->SetAltVariable("elresup",sublead_jet_pt);
      bdt_reader->SetAltVariable("elresup",
                                 sys_llphoton_refit_jet_balance_elresup);
      bdt_reader->SetAltVariable("elresup",
                                 sys_llphoton_refit_dijet_dphi_elresup);
      bdt_reader->SetAltVariable("elresup",sys_llphoton_refit_psi_elresup);
      bdt_reader->SetAltVariable("elresup","photon_drmin[0]");
      bdt_reader->SetAltVariable("elresup","photon_drmax[0]");
      bdt_reader->SetAltVariable("elresup",
                                 sys_llphoton_refit_costheta_elresup);
      bdt_reader->SetAltVariable("elresup",
                                 sys_llphoton_refit_cosTheta_elresup);
      bdt_reader->SetAltVariable("elresup",
                                 sys_llphoton_refit_pt_elresup
                                 /sys_llphoton_refit_m_elresup);
      bdt_reader->SetAltVariable("elresup",sys_lead_lepton_eta_elresup);
      bdt_reader->SetAltVariable("elresup",sys_sublead_lepton_eta_elresup);
      bdt_reader->SetAltVariable("elresup","photon_eta[0]");
      bdt_reader->SetAltVariable("elresup","photon_idmva[0]");
      bdt_reader->SetAltVariable("elresup",photon_relpterr);
      bdt_reader->SetAltVariable("elresup",
                                 sys_photon_zeppenfeld_default);
      bdt_reader->SetAltVariable("elresup",sys_photon_jet1_dr_default);
      bdt_reader->SetAltVariable("elresup",sys_photon_jet2_dr_default);
      //elresdn
      bdt_reader->SetAltVariable("elresdn",sys_dijet_deta_default);
      bdt_reader->SetAltVariable("elresdn",sys_dijet_dphi_default);
      bdt_reader->SetAltVariable("elresdn",lead_jet_pt);
      bdt_reader->SetAltVariable("elresdn",sublead_jet_pt);
      bdt_reader->SetAltVariable("elresdn",
                                 sys_llphoton_refit_jet_balance_elresdn);
      bdt_reader->SetAltVariable("elresdn",
                                 sys_llphoton_refit_dijet_dphi_elresdn);
      bdt_reader->SetAltVariable("elresdn",sys_llphoton_refit_psi_elresdn);
      bdt_reader->SetAltVariable("elresdn","photon_drmin[0]");
      bdt_reader->SetAltVariable("elresdn","photon_drmax[0]");
      bdt_reader->SetAltVariable("elresdn",
                                 sys_llphoton_refit_costheta_elresdn);
      bdt_reader->SetAltVariable("elresdn",
                                 sys_llphoton_refit_cosTheta_elresdn);
      bdt_reader->SetAltVariable("elresdn",
                                 sys_llphoton_refit_pt_elresdn
                                 /sys_llphoton_refit_m_elresdn);
      bdt_reader->SetAltVariable("elresdn",sys_lead_lepton_eta_elresdn);
      bdt_reader->SetAltVariable("elresdn",sys_sublead_lepton_eta_elresdn);
      bdt_reader->SetAltVariable("elresdn","photon_eta[0]");
      bdt_reader->SetAltVariable("elresdn","photon_idmva[0]");
      bdt_reader->SetAltVariable("elresdn",photon_relpterr);
      bdt_reader->SetAltVariable("elresdn",
                                 sys_photon_zeppenfeld_default);
      bdt_reader->SetAltVariable("elresdn",sys_photon_jet1_dr_default);
      bdt_reader->SetAltVariable("elresdn",sys_photon_jet2_dr_default);
      //muscaleup
      bdt_reader->SetAltVariable("muscaleup",sys_dijet_deta_default);
      bdt_reader->SetAltVariable("muscaleup",sys_dijet_dphi_default);
      bdt_reader->SetAltVariable("muscaleup",lead_jet_pt);
      bdt_reader->SetAltVariable("muscaleup",sublead_jet_pt);
      bdt_reader->SetAltVariable("muscaleup",
                                 sys_llphoton_refit_jet_balance_muscaleup);
      bdt_reader->SetAltVariable("muscaleup",
                                 sys_llphoton_refit_dijet_dphi_muscaleup);
      bdt_reader->SetAltVariable("muscaleup",sys_llphoton_refit_psi_muscaleup);
      bdt_reader->SetAltVariable("muscaleup","photon_drmin[0]");
      bdt_reader->SetAltVariable("muscaleup","photon_drmax[0]");
      bdt_reader->SetAltVariable("muscaleup",
                                 sys_llphoton_refit_costheta_muscaleup);
      bdt_reader->SetAltVariable("muscaleup",
                                 sys_llphoton_refit_cosTheta_muscaleup);
      bdt_reader->SetAltVariable("muscaleup",
                                 sys_llphoton_refit_pt_muscaleup
                                 /sys_llphoton_refit_m_muscaleup);
      bdt_reader->SetAltVariable("muscaleup",sys_lead_lepton_eta_muscaleup);
      bdt_reader->SetAltVariable("muscaleup",sys_sublead_lepton_eta_muscaleup);
      bdt_reader->SetAltVariable("muscaleup","photon_eta[0]");
      bdt_reader->SetAltVariable("muscaleup","photon_idmva[0]");
      bdt_reader->SetAltVariable("muscaleup",photon_relpterr);
      bdt_reader->SetAltVariable("muscaleup",
                                 sys_photon_zeppenfeld_default);
      bdt_reader->SetAltVariable("muscaleup",sys_photon_jet1_dr_default);
      bdt_reader->SetAltVariable("muscaleup",sys_photon_jet2_dr_default);
      //muscaledn
      bdt_reader->SetAltVariable("muscaledn",sys_dijet_deta_default);
      bdt_reader->SetAltVariable("muscaledn",sys_dijet_dphi_default);
      bdt_reader->SetAltVariable("muscaledn",lead_jet_pt);
      bdt_reader->SetAltVariable("muscaledn",sublead_jet_pt);
      bdt_reader->SetAltVariable("muscaledn",
                                 sys_llphoton_refit_jet_balance_muscaledn);
      bdt_reader->SetAltVariable("muscaledn",
                                 sys_llphoton_refit_dijet_dphi_muscaledn);
      bdt_reader->SetAltVariable("muscaledn",sys_llphoton_refit_psi_muscaledn);
      bdt_reader->SetAltVariable("muscaledn","photon_drmin[0]");
      bdt_reader->SetAltVariable("muscaledn","photon_drmax[0]");
      bdt_reader->SetAltVariable("muscaledn",
                                 sys_llphoton_refit_costheta_muscaledn);
      bdt_reader->SetAltVariable("muscaledn",
                                 sys_llphoton_refit_cosTheta_muscaledn);
      bdt_reader->SetAltVariable("muscaledn",
                                 sys_llphoton_refit_pt_muscaledn
                                 /sys_llphoton_refit_m_muscaledn);
      bdt_reader->SetAltVariable("muscaledn",sys_lead_lepton_eta_muscaledn);
      bdt_reader->SetAltVariable("muscaledn",sys_sublead_lepton_eta_muscaledn);
      bdt_reader->SetAltVariable("muscaledn","photon_eta[0]");
      bdt_reader->SetAltVariable("muscaledn","photon_idmva[0]");
      bdt_reader->SetAltVariable("muscaledn",photon_relpterr);
      bdt_reader->SetAltVariable("muscaledn",
                                 sys_photon_zeppenfeld_default);
      bdt_reader->SetAltVariable("muscaledn",sys_photon_jet1_dr_default);
      bdt_reader->SetAltVariable("muscaledn",sys_photon_jet2_dr_default);
      //muresup
      bdt_reader->SetAltVariable("muresup",sys_dijet_deta_default);
      bdt_reader->SetAltVariable("muresup",sys_dijet_dphi_default);
      bdt_reader->SetAltVariable("muresup",lead_jet_pt);
      bdt_reader->SetAltVariable("muresup",sublead_jet_pt);
      bdt_reader->SetAltVariable("muresup",
                                 sys_llphoton_refit_jet_balance_muresup);
      bdt_reader->SetAltVariable("muresup",
                                 sys_llphoton_refit_dijet_dphi_muresup);
      bdt_reader->SetAltVariable("muresup",sys_llphoton_refit_psi_muresup);
      bdt_reader->SetAltVariable("muresup","photon_drmin[0]");
      bdt_reader->SetAltVariable("muresup","photon_drmax[0]");
      bdt_reader->SetAltVariable("muresup",
                                 sys_llphoton_refit_costheta_muresup);
      bdt_reader->SetAltVariable("muresup",
                                 sys_llphoton_refit_cosTheta_muresup);
      bdt_reader->SetAltVariable("muresup",
                                 sys_llphoton_refit_pt_muresup
                                 /sys_llphoton_refit_m_muresup);
      bdt_reader->SetAltVariable("muresup",sys_lead_lepton_eta_muresup);
      bdt_reader->SetAltVariable("muresup",sys_sublead_lepton_eta_muresup);
      bdt_reader->SetAltVariable("muresup","photon_eta[0]");
      bdt_reader->SetAltVariable("muresup","photon_idmva[0]");
      bdt_reader->SetAltVariable("muresup",photon_relpterr);
      bdt_reader->SetAltVariable("muresup",
                                 sys_photon_zeppenfeld_default);
      bdt_reader->SetAltVariable("muresup",sys_photon_jet1_dr_default);
      bdt_reader->SetAltVariable("muresup",sys_photon_jet2_dr_default);
      //muresdn
      bdt_reader->SetAltVariable("muresdn",sys_dijet_deta_default);
      bdt_reader->SetAltVariable("muresdn",sys_dijet_dphi_default);
      bdt_reader->SetAltVariable("muresdn",lead_jet_pt);
      bdt_reader->SetAltVariable("muresdn",sublead_jet_pt);
      bdt_reader->SetAltVariable("muresdn",
                                 sys_llphoton_refit_jet_balance_muresdn);
      bdt_reader->SetAltVariable("muresdn",
                                 sys_llphoton_refit_dijet_dphi_muresdn);
      bdt_reader->SetAltVariable("muresdn",sys_llphoton_refit_psi_muresdn);
      bdt_reader->SetAltVariable("muresdn","photon_drmin[0]");
      bdt_reader->SetAltVariable("muresdn","photon_drmax[0]");
      bdt_reader->SetAltVariable("muresdn",
                                 sys_llphoton_refit_costheta_muresdn);
      bdt_reader->SetAltVariable("muresdn",
                                 sys_llphoton_refit_cosTheta_muresdn);
      bdt_reader->SetAltVariable("muresdn",
                                 sys_llphoton_refit_pt_muresdn
                                 /sys_llphoton_refit_m_muresdn);
      bdt_reader->SetAltVariable("muresdn",sys_lead_lepton_eta_muresdn);
      bdt_reader->SetAltVariable("muresdn",sys_sublead_lepton_eta_muresdn);
      bdt_reader->SetAltVariable("muresdn","photon_eta[0]");
      bdt_reader->SetAltVariable("muresdn","photon_idmva[0]");
      bdt_reader->SetAltVariable("muresdn",photon_relpterr);
      bdt_reader->SetAltVariable("muresdn",
                                 sys_photon_zeppenfeld_default);
      bdt_reader->SetAltVariable("muresdn",sys_photon_jet1_dr_default);
      bdt_reader->SetAltVariable("muresdn",sys_photon_jet2_dr_default);
      //default
      bdt_reader->SetVariable("dijet_deta",sys_dijet_deta_default);
      bdt_reader->SetVariable("dijet_dphi",sys_dijet_dphi_default);
      bdt_reader->SetVariable("j1_pt",lead_jet_pt);
      bdt_reader->SetVariable("j2_pt",sublead_jet_pt);
      bdt_reader->SetVariable("llphoton_dijet_balance",
                              sys_llphoton_refit_jet_balance_default);
      bdt_reader->SetVariable("llphoton_dijet_dphi",
                              sys_llphoton_refit_dijet_dphi_default);
      bdt_reader->SetVariable("phi",sys_llphoton_refit_psi_default);
      bdt_reader->SetVariable("min_dR","photon_drmin[0]");
      bdt_reader->SetVariable("max_dR","photon_drmax[0]");
      bdt_reader->SetVariable("costheta",sys_llphoton_refit_costheta_default);
      bdt_reader->SetVariable("cosTheta",sys_llphoton_refit_cosTheta_default);
      bdt_reader->SetVariable("pt_mass","llphoton_refit_pt/llphoton_refit_m");
      bdt_reader->SetVariable("l1_rapidity",sys_lead_lepton_eta_default);
      bdt_reader->SetVariable("l2_rapidity",sys_sublead_lepton_eta_default);
      bdt_reader->SetVariable("photon_rapidity","photon_eta[0]");
      bdt_reader->SetVariable("photon_mva","photon_idmva[0]");
      bdt_reader->SetVariable("photon_res",photon_relpterr);
      bdt_reader->SetVariable("photon_zeppenfeld",
                              sys_photon_zeppenfeld_default);
      bdt_reader->SetVariable("photon_jet1_dr",sys_photon_jet1_dr_default);
      bdt_reader->SetVariable("photon_jet2_dr",sys_photon_jet2_dr_default);
      //phscaleup
      bdt_reader->SetAltVariable("phscaleup",sys_dijet_deta_default);
      bdt_reader->SetAltVariable("phscaleup",sys_dijet_dphi_default);
      bdt_reader->SetAltVariable("phscaleup",lead_jet_pt);
      bdt_reader->SetAltVariable("phscaleup",sublead_jet_pt);
      bdt_reader->SetAltVariable("phscaleup",
                                 sys_llphoton_refit_jet_balance_phscaleup);
      bdt_reader->SetAltVariable("phscaleup",
                                 sys_llphoton_refit_dijet_dphi_phscaleup);
      bdt_reader->SetAltVariable("phscaleup",sys_llphoton_refit_psi_phscaleup);
      bdt_reader->SetAltVariable("phscaleup",sys_lead_photon_drmin_scaleup);
      bdt_reader->SetAltVariable("phscaleup",sys_lead_photon_drmax_scaleup);
      bdt_reader->SetAltVariable("phscaleup",
                                 sys_llphoton_refit_costheta_phscaleup);
      bdt_reader->SetAltVariable("phscaleup",
                                 sys_llphoton_refit_cosTheta_phscaleup);
      bdt_reader->SetAltVariable("phscaleup",
                                 sys_llphoton_refit_pt_phscaleup
                                 /sys_llphoton_refit_m_phscaleup);
      bdt_reader->SetAltVariable("phscaleup",sys_lead_lepton_eta_default);
      bdt_reader->SetAltVariable("phscaleup",sys_sublead_lepton_eta_default);
      bdt_reader->SetAltVariable("phscaleup",sys_lead_photon_eta_scaleup);
      bdt_reader->SetAltVariable("phscaleup",sys_lead_photon_idmva_scaleup);
      bdt_reader->SetAltVariable("phscaleup",
                                 sys_lead_photon_relpterr_scaleup);
      bdt_reader->SetAltVariable("phscaleup",
                                 sys_photon_zeppenfeld_phscaleup);
      bdt_reader->SetAltVariable("phscaleup",sys_photon_jet1_dr_phscaleup);
      bdt_reader->SetAltVariable("phscaleup",sys_photon_jet2_dr_phscaleup);
      //phscaledn
      bdt_reader->SetAltVariable("phscaledn",sys_dijet_deta_default);
      bdt_reader->SetAltVariable("phscaledn",sys_dijet_dphi_default);
      bdt_reader->SetAltVariable("phscaledn",lead_jet_pt);
      bdt_reader->SetAltVariable("phscaledn",sublead_jet_pt);
      bdt_reader->SetAltVariable("phscaledn",
                                 sys_llphoton_refit_jet_balance_phscaledn);
      bdt_reader->SetAltVariable("phscaledn",
                                 sys_llphoton_refit_dijet_dphi_phscaledn);
      bdt_reader->SetAltVariable("phscaledn",sys_llphoton_refit_psi_phscaledn);
      bdt_reader->SetAltVariable("phscaledn",sys_lead_photon_drmin_scaledn);
      bdt_reader->SetAltVariable("phscaledn",sys_lead_photon_drmax_scaledn);
      bdt_reader->SetAltVariable("phscaledn",
                                 sys_llphoton_refit_costheta_phscaledn);
      bdt_reader->SetAltVariable("phscaledn",
                                 sys_llphoton_refit_cosTheta_phscaledn);
      bdt_reader->SetAltVariable("phscaledn",
                                 sys_llphoton_refit_pt_phscaledn
                                 /sys_llphoton_refit_m_phscaledn);
      bdt_reader->SetAltVariable("phscaledn",sys_lead_lepton_eta_default);
      bdt_reader->SetAltVariable("phscaledn",sys_sublead_lepton_eta_default);
      bdt_reader->SetAltVariable("phscaledn",sys_lead_photon_eta_scaledn);
      bdt_reader->SetAltVariable("phscaledn",sys_lead_photon_idmva_scaledn);
      bdt_reader->SetAltVariable("phscaledn",
                                 sys_lead_photon_relpterr_scaledn);
      bdt_reader->SetAltVariable("phscaledn",
                                 sys_photon_zeppenfeld_phscaledn);
      bdt_reader->SetAltVariable("phscaledn",sys_photon_jet1_dr_phscaledn);
      bdt_reader->SetAltVariable("phscaledn",sys_photon_jet2_dr_phscaledn);
      //phresup
      bdt_reader->SetAltVariable("phresup",sys_dijet_deta_default);
      bdt_reader->SetAltVariable("phresup",sys_dijet_dphi_default);
      bdt_reader->SetAltVariable("phresup",lead_jet_pt);
      bdt_reader->SetAltVariable("phresup",sublead_jet_pt);
      bdt_reader->SetAltVariable("phresup",
                                 sys_llphoton_refit_jet_balance_phresup);
      bdt_reader->SetAltVariable("phresup",
                                 sys_llphoton_refit_dijet_dphi_phresup);
      bdt_reader->SetAltVariable("phresup",sys_llphoton_refit_psi_phresup);
      bdt_reader->SetAltVariable("phresup",sys_lead_photon_drmin_resup);
      bdt_reader->SetAltVariable("phresup",sys_lead_photon_drmax_resup);
      bdt_reader->SetAltVariable("phresup",
                                 sys_llphoton_refit_costheta_phresup);
      bdt_reader->SetAltVariable("phresup",
                                 sys_llphoton_refit_cosTheta_phresup);
      bdt_reader->SetAltVariable("phresup",
                                 sys_llphoton_refit_pt_phresup
                                 /sys_llphoton_refit_m_phresup);
      bdt_reader->SetAltVariable("phresup",sys_lead_lepton_eta_default);
      bdt_reader->SetAltVariable("phresup",sys_sublead_lepton_eta_default);
      bdt_reader->SetAltVariable("phresup",sys_lead_photon_eta_resup);
      bdt_reader->SetAltVariable("phresup",sys_lead_photon_idmva_resup);
      bdt_reader->SetAltVariable("phresup",
                                 sys_lead_photon_relpterr_resup);
      bdt_reader->SetAltVariable("phresup",
                                 sys_photon_zeppenfeld_phresup);
      bdt_reader->SetAltVariable("phresup",sys_photon_jet1_dr_phresup);
      bdt_reader->SetAltVariable("phresup",sys_photon_jet2_dr_phresup);
      //phresdn
      bdt_reader->SetAltVariable("phresdn",sys_dijet_deta_default);
      bdt_reader->SetAltVariable("phresdn",sys_dijet_dphi_default);
      bdt_reader->SetAltVariable("phresdn",lead_jet_pt);
      bdt_reader->SetAltVariable("phresdn",sublead_jet_pt);
      bdt_reader->SetAltVariable("phresdn",
                                 sys_llphoton_refit_jet_balance_phresdn);
      bdt_reader->SetAltVariable("phresdn",
                                 sys_llphoton_refit_dijet_dphi_phresdn);
      bdt_reader->SetAltVariable("phresdn",sys_llphoton_refit_psi_phresdn);
      bdt_reader->SetAltVariable("phresdn",sys_lead_photon_drmin_resdn);
      bdt_reader->SetAltVariable("phresdn",sys_lead_photon_drmax_resdn);
      bdt_reader->SetAltVariable("phresdn",
                                 sys_llphoton_refit_costheta_phresdn);
      bdt_reader->SetAltVariable("phresdn",
                                 sys_llphoton_refit_cosTheta_phresdn);
      bdt_reader->SetAltVariable("phresdn",
                                 sys_llphoton_refit_pt_phresdn
                                 /sys_llphoton_refit_m_phresdn);
      bdt_reader->SetAltVariable("phresdn",sys_lead_lepton_eta_default);
      bdt_reader->SetAltVariable("phresdn",sys_sublead_lepton_eta_default);
      bdt_reader->SetAltVariable("phresdn",sys_lead_photon_eta_resdn);
      bdt_reader->SetAltVariable("phresdn",sys_lead_photon_idmva_resdn);
      bdt_reader->SetAltVariable("phresdn",
                                 sys_lead_photon_relpterr_resdn);
      bdt_reader->SetAltVariable("phresdn",
                                 sys_photon_zeppenfeld_phresdn);
      bdt_reader->SetAltVariable("phresdn",sys_photon_jet1_dr_phresdn);
      bdt_reader->SetAltVariable("phresdn",sys_photon_jet2_dr_phresdn);
      const vector<string> years = {"2016APV", "2016", "2017", "2018", 
                                    "2022", "2022EE", "2023", "2023BPix"};
      for (unsigned iyear = 0; iyear < years.size(); iyear++) {
        string year = years[iyear];
        //jetscaleup
        bdt_reader->SetAltVariable("jetscaleup"+year,
                                   sys_dijet_deta_scaleup[iyear]);
        bdt_reader->SetAltVariable("jetscaleup"+year,
                                   sys_dijet_dphi_scaleup[iyear]);
        bdt_reader->SetAltVariable("jetscaleup"+year,
                                   sys_lead_jet_pt_scaleup[iyear]);
        bdt_reader->SetAltVariable("jetscaleup"+year,
                                   sys_sublead_jet_pt_scaleup[iyear]);
        bdt_reader->SetAltVariable("jetscaleup"+year,
            sys_llphoton_refit_jet_balance_jetscaleup[iyear]);
        bdt_reader->SetAltVariable("jetscaleup"+year,
            sys_llphoton_refit_dijet_dphi_jetscaleup[iyear]);
        bdt_reader->SetAltVariable("jetscaleup"+year,
                                   sys_llphoton_refit_psi_default);
        bdt_reader->SetAltVariable("jetscaleup"+year,"photon_drmin[0]");
        bdt_reader->SetAltVariable("jetscaleup"+year,"photon_drmax[0]");
        bdt_reader->SetAltVariable("jetscaleup"+year,
                                   sys_llphoton_refit_costheta_default);
        bdt_reader->SetAltVariable("jetscaleup"+year,
                                   sys_llphoton_refit_cosTheta_default);
        bdt_reader->SetAltVariable("jetscaleup"+year,
                                   "llphoton_refit_pt/llphoton_refit_m");
        bdt_reader->SetAltVariable("jetscaleup"+year,
                                   sys_lead_lepton_eta_default);
        bdt_reader->SetAltVariable("jetscaleup"+year,
                                   sys_sublead_lepton_eta_default);
        bdt_reader->SetAltVariable("jetscaleup"+year,"photon_eta[0]");
        bdt_reader->SetAltVariable("jetscaleup"+year,"photon_idmva[0]");
        bdt_reader->SetAltVariable("jetscaleup"+year,photon_relpterr);
        bdt_reader->SetAltVariable("jetscaleup"+year,
                                   sys_photon_zeppenfeld_jetscaleup[iyear]);
        bdt_reader->SetAltVariable("jetscaleup"+year,
                                   sys_photon_jet1_dr_jetscaleup[iyear]);
        bdt_reader->SetAltVariable("jetscaleup"+year,
                                   sys_photon_jet2_dr_jetscaleup[iyear]);
        //jetscaledn
        bdt_reader->SetAltVariable("jetscaledn"+year,
                                   sys_dijet_deta_scaledn[iyear]);
        bdt_reader->SetAltVariable("jetscaledn"+year,
                                   sys_dijet_dphi_scaledn[iyear]);
        bdt_reader->SetAltVariable("jetscaledn"+year,
                                   sys_lead_jet_pt_scaledn[iyear]);
        bdt_reader->SetAltVariable("jetscaledn"+year,
                                   sys_sublead_jet_pt_scaledn[iyear]);
        bdt_reader->SetAltVariable("jetscaledn"+year,
            sys_llphoton_refit_jet_balance_jetscaledn[iyear]);
        bdt_reader->SetAltVariable("jetscaledn"+year,
            sys_llphoton_refit_dijet_dphi_jetscaledn[iyear]);
        bdt_reader->SetAltVariable("jetscaledn"+year,
                                   sys_llphoton_refit_psi_default);
        bdt_reader->SetAltVariable("jetscaledn"+year,"photon_drmin[0]");
        bdt_reader->SetAltVariable("jetscaledn"+year,"photon_drmax[0]");
        bdt_reader->SetAltVariable("jetscaledn"+year,
                                   sys_llphoton_refit_costheta_default);
        bdt_reader->SetAltVariable("jetscaledn"+year,
                                   sys_llphoton_refit_cosTheta_default);
        bdt_reader->SetAltVariable("jetscaledn"+year,
                                   "llphoton_refit_pt/llphoton_refit_m");
        bdt_reader->SetAltVariable("jetscaledn"+year,
                                   sys_lead_lepton_eta_default);
        bdt_reader->SetAltVariable("jetscaledn"+year,
                                   sys_sublead_lepton_eta_default);
        bdt_reader->SetAltVariable("jetscaledn"+year,"photon_eta[0]");
        bdt_reader->SetAltVariable("jetscaledn"+year,"photon_idmva[0]");
        bdt_reader->SetAltVariable("jetscaledn"+year,photon_relpterr);
        bdt_reader->SetAltVariable("jetscaledn"+year,
                                   sys_photon_zeppenfeld_jetscaledn[iyear]);
        bdt_reader->SetAltVariable("jetscaledn"+year,
                                   sys_photon_jet1_dr_jetscaledn[iyear]);
        bdt_reader->SetAltVariable("jetscaledn"+year,
                                   sys_photon_jet2_dr_jetscaledn[iyear]);
        //jetresup
        bdt_reader->SetAltVariable("jetresup"+year,
                                   sys_dijet_deta_resup[iyear]);
        bdt_reader->SetAltVariable("jetresup"+year,
                                   sys_dijet_dphi_resup[iyear]);
        bdt_reader->SetAltVariable("jetresup"+year,
                                   sys_lead_jet_pt_resup[iyear]);
        bdt_reader->SetAltVariable("jetresup"+year,
                                   sys_sublead_jet_pt_resup[iyear]);
        bdt_reader->SetAltVariable("jetresup"+year,
            sys_llphoton_refit_jet_balance_jetresup[iyear]);
        bdt_reader->SetAltVariable("jetresup"+year,
            sys_llphoton_refit_dijet_dphi_jetresup[iyear]);
        bdt_reader->SetAltVariable("jetresup"+year,
                                   sys_llphoton_refit_psi_default);
        bdt_reader->SetAltVariable("jetresup"+year,"photon_drmin[0]");
        bdt_reader->SetAltVariable("jetresup"+year,"photon_drmax[0]");
        bdt_reader->SetAltVariable("jetresup"+year,
                                   sys_llphoton_refit_costheta_default);
        bdt_reader->SetAltVariable("jetresup"+year,
                                   sys_llphoton_refit_cosTheta_default);
        bdt_reader->SetAltVariable("jetresup"+year,
                                   "llphoton_refit_pt/llphoton_refit_m");
        bdt_reader->SetAltVariable("jetresup"+year,
                                   sys_lead_lepton_eta_default);
        bdt_reader->SetAltVariable("jetresup"+year,
                                   sys_sublead_lepton_eta_default);
        bdt_reader->SetAltVariable("jetresup"+year,"photon_eta[0]");
        bdt_reader->SetAltVariable("jetresup"+year,"photon_idmva[0]");
        bdt_reader->SetAltVariable("jetresup"+year,photon_relpterr);
        bdt_reader->SetAltVariable("jetresup"+year,
                                   sys_photon_zeppenfeld_jetresup[iyear]);
        bdt_reader->SetAltVariable("jetresup"+year,
                                   sys_photon_jet1_dr_jetresup[iyear]);
        bdt_reader->SetAltVariable("jetresup"+year,
                                   sys_photon_jet2_dr_jetresup[iyear]);
        //jetresdn
        bdt_reader->SetAltVariable("jetresdn"+year,
                                   sys_dijet_deta_resdn[iyear]);
        bdt_reader->SetAltVariable("jetresdn"+year,
                                   sys_dijet_dphi_resdn[iyear]);
        bdt_reader->SetAltVariable("jetresdn"+year,
                                   sys_lead_jet_pt_resdn[iyear]);
        bdt_reader->SetAltVariable("jetresdn"+year,
                                   sys_sublead_jet_pt_resdn[iyear]);
        bdt_reader->SetAltVariable("jetresdn"+year,
            sys_llphoton_refit_jet_balance_jetresdn[iyear]);
        bdt_reader->SetAltVariable("jetresdn"+year,
            sys_llphoton_refit_dijet_dphi_jetresdn[iyear]);
        bdt_reader->SetAltVariable("jetresdn"+year,
                                   sys_llphoton_refit_psi_default);
        bdt_reader->SetAltVariable("jetresdn"+year,"photon_drmin[0]");
        bdt_reader->SetAltVariable("jetresdn"+year,"photon_drmax[0]");
        bdt_reader->SetAltVariable("jetresdn"+year,
                                   sys_llphoton_refit_costheta_default);
        bdt_reader->SetAltVariable("jetresdn"+year,
                                   sys_llphoton_refit_cosTheta_default);
        bdt_reader->SetAltVariable("jetresdn"+year,
                                   "llphoton_refit_pt/llphoton_refit_m");
        bdt_reader->SetAltVariable("jetresdn"+year,
                                   sys_lead_lepton_eta_default);
        bdt_reader->SetAltVariable("jetresdn"+year,
                                   sys_sublead_lepton_eta_default);
        bdt_reader->SetAltVariable("jetresdn"+year,"photon_eta[0]");
        bdt_reader->SetAltVariable("jetresdn"+year,"photon_idmva[0]");
        bdt_reader->SetAltVariable("jetresdn"+year,photon_relpterr);
        bdt_reader->SetAltVariable("jetresdn"+year,
                                   sys_photon_zeppenfeld_jetresdn[iyear]);
        bdt_reader->SetAltVariable("jetresdn"+year,
                                   sys_photon_jet1_dr_jetresdn[iyear]);
        bdt_reader->SetAltVariable("jetresdn"+year,
                                   sys_photon_jet2_dr_jetresdn[iyear]);
      }

    }
    vbf_bdt_readers[0]->BookMVA("/net/cms37/data1/rui/Training/dataset_2JClassic_pinnacles_run2p3/weights/TMVAClassification_BDT_0.weights.xml");
    vbf_bdt_readers[1]->BookMVA("/net/cms37/data1/rui/Training/dataset_2JClassic_pinnacles_run2p3/weights/TMVAClassification_BDT_1.weights.xml");
    vbf_bdt_readers[2]->BookMVA("/net/cms37/data1/rui/Training/dataset_2JClassic_pinnacles_run2p3/weights/TMVAClassification_BDT_2.weights.xml");
    vbf_bdt_readers[3]->BookMVA("/net/cms37/data1/rui/Training/dataset_2JClassic_pinnacles_run2p3/weights/TMVAClassification_BDT_3.weights.xml");
    return vbf_bdt_readers;
  }

  NamedFunc ggfbdt2503_score_default(const vector<FastForest> &xgb_bdts) {
     return XGBoostBDTScore(xgb_bdts, {NamedFunc("photon_idmva[0]"), 
         photon_relpterr, NamedFunc("photon_drmax[0]"), 
         NamedFunc("photon_drmin[0]"), 
         NamedFunc("llphoton_refit_pt/llphoton_refit_m"), 
         sys_llphoton_refit_cosTheta_default, 
         sys_llphoton_refit_costheta_default, sys_llphoton_refit_psi_default,
         sys_lead_lepton_eta_default, sys_sublead_lepton_eta_default, 
         NamedFunc("photon_eta[0]")});
  }

  NamedFunc ggfbdt2503_score_elscaleup(const vector<FastForest> &xgb_bdts) {
     return XGBoostBDTScore(xgb_bdts, {NamedFunc("photon_idmva[0]"), 
         photon_relpterr, NamedFunc("photon_drmax[0]"), 
         NamedFunc("photon_drmin[0]"), 
         sys_llphoton_refit_pt_elscaleup/sys_llphoton_refit_m_elscaleup, 
         sys_llphoton_refit_cosTheta_elscaleup, 
         sys_llphoton_refit_costheta_elscaleup, 
         sys_llphoton_refit_psi_elscaleup,
         sys_lead_lepton_eta_elscaleup, sys_sublead_lepton_eta_elscaleup, 
         NamedFunc("photon_eta[0]")});
  }

  NamedFunc ggfbdt2503_score_elscaledn(const vector<FastForest> &xgb_bdts) {
     return XGBoostBDTScore(xgb_bdts, {NamedFunc("photon_idmva[0]"), 
         photon_relpterr, NamedFunc("photon_drmax[0]"), 
         NamedFunc("photon_drmin[0]"), 
         sys_llphoton_refit_pt_elscaledn/sys_llphoton_refit_m_elscaledn, 
         sys_llphoton_refit_cosTheta_elscaledn, 
         sys_llphoton_refit_costheta_elscaledn, 
         sys_llphoton_refit_psi_elscaledn,
         sys_lead_lepton_eta_elscaledn, sys_sublead_lepton_eta_elscaledn, 
         NamedFunc("photon_eta[0]")});
  }

  NamedFunc ggfbdt2503_score_elresup(const vector<FastForest> &xgb_bdts) {
     return XGBoostBDTScore(xgb_bdts, {NamedFunc("photon_idmva[0]"), 
         photon_relpterr, NamedFunc("photon_drmax[0]"), 
         NamedFunc("photon_drmin[0]"), 
         sys_llphoton_refit_pt_elresup/sys_llphoton_refit_m_elresup, 
         sys_llphoton_refit_cosTheta_elresup, 
         sys_llphoton_refit_costheta_elresup, 
         sys_llphoton_refit_psi_elresup,
         sys_lead_lepton_eta_elresup, sys_sublead_lepton_eta_elresup, 
         NamedFunc("photon_eta[0]")});
  }

  NamedFunc ggfbdt2503_score_elresdn(const vector<FastForest> &xgb_bdts) {
     return XGBoostBDTScore(xgb_bdts, {NamedFunc("photon_idmva[0]"), 
         photon_relpterr, NamedFunc("photon_drmax[0]"), 
         NamedFunc("photon_drmin[0]"), 
         sys_llphoton_refit_pt_elresdn/sys_llphoton_refit_m_elresdn, 
         sys_llphoton_refit_cosTheta_elresdn, 
         sys_llphoton_refit_costheta_elresdn, 
         sys_llphoton_refit_psi_elresdn,
         sys_lead_lepton_eta_elresdn, sys_sublead_lepton_eta_elresdn, 
         NamedFunc("photon_eta[0]")});
  }

  NamedFunc ggfbdt2503_score_muscaleup(const vector<FastForest> &xgb_bdts) {
     return XGBoostBDTScore(xgb_bdts, {NamedFunc("photon_idmva[0]"), 
         photon_relpterr, NamedFunc("photon_drmax[0]"), 
         NamedFunc("photon_drmin[0]"), 
         sys_llphoton_refit_pt_muscaleup/sys_llphoton_refit_m_muscaleup, 
         sys_llphoton_refit_cosTheta_muscaleup, 
         sys_llphoton_refit_costheta_muscaleup, 
         sys_llphoton_refit_psi_muscaleup,
         sys_lead_lepton_eta_muscaleup, sys_sublead_lepton_eta_muscaleup, 
         NamedFunc("photon_eta[0]")});
  }

  NamedFunc ggfbdt2503_score_muscaledn(const vector<FastForest> &xgb_bdts) {
     return XGBoostBDTScore(xgb_bdts, {NamedFunc("photon_idmva[0]"), 
         photon_relpterr, NamedFunc("photon_drmax[0]"), 
         NamedFunc("photon_drmin[0]"), 
         sys_llphoton_refit_pt_muscaledn/sys_llphoton_refit_m_muscaledn, 
         sys_llphoton_refit_cosTheta_muscaledn, 
         sys_llphoton_refit_costheta_muscaledn, 
         sys_llphoton_refit_psi_muscaledn,
         sys_lead_lepton_eta_muscaledn, sys_sublead_lepton_eta_muscaledn, 
         NamedFunc("photon_eta[0]")});
  }

  NamedFunc ggfbdt2503_score_muresup(const vector<FastForest> &xgb_bdts) {
     return XGBoostBDTScore(xgb_bdts, {NamedFunc("photon_idmva[0]"), 
         photon_relpterr, NamedFunc("photon_drmax[0]"), 
         NamedFunc("photon_drmin[0]"), 
         sys_llphoton_refit_pt_muresup/sys_llphoton_refit_m_muresup, 
         sys_llphoton_refit_cosTheta_muresup, 
         sys_llphoton_refit_costheta_muresup, 
         sys_llphoton_refit_psi_muresup,
         sys_lead_lepton_eta_muresup, sys_sublead_lepton_eta_muresup, 
         NamedFunc("photon_eta[0]")});
  }

  NamedFunc ggfbdt2503_score_muresdn(const vector<FastForest> &xgb_bdts) {
     return XGBoostBDTScore(xgb_bdts, {NamedFunc("photon_idmva[0]"), 
         photon_relpterr, NamedFunc("photon_drmax[0]"), 
         NamedFunc("photon_drmin[0]"), 
         sys_llphoton_refit_pt_muresdn/sys_llphoton_refit_m_muresdn, 
         sys_llphoton_refit_cosTheta_muresdn, 
         sys_llphoton_refit_costheta_muresdn, 
         sys_llphoton_refit_psi_muresdn,
         sys_lead_lepton_eta_muresdn, sys_sublead_lepton_eta_muresdn, 
         NamedFunc("photon_eta[0]")});
  }

  NamedFunc ggfbdt2503_score_phscaleup(const vector<FastForest> &xgb_bdts) {
     return XGBoostBDTScore(xgb_bdts, {sys_lead_photon_idmva_scaleup, 
         sys_lead_photon_relpterr_scaleup, sys_lead_photon_drmax_scaleup, 
         sys_lead_photon_drmin_scaleup, 
         sys_llphoton_refit_pt_phscaleup/sys_llphoton_refit_m_phscaleup, 
         sys_llphoton_refit_cosTheta_phscaleup, 
         sys_llphoton_refit_costheta_phscaleup, 
         sys_llphoton_refit_psi_phscaleup,
         sys_lead_lepton_eta_default, sys_sublead_lepton_eta_default, 
         sys_lead_photon_eta_scaleup});
  }

  NamedFunc ggfbdt2503_score_phscaledn(const vector<FastForest> &xgb_bdts) {
     return XGBoostBDTScore(xgb_bdts, {sys_lead_photon_idmva_scaleup, 
         sys_lead_photon_relpterr_scaleup, sys_lead_photon_drmax_scaleup, 
         sys_lead_photon_drmin_scaleup, 
         sys_llphoton_refit_pt_phscaledn/sys_llphoton_refit_m_phscaledn, 
         sys_llphoton_refit_cosTheta_phscaledn, 
         sys_llphoton_refit_costheta_phscaledn, 
         sys_llphoton_refit_psi_phscaledn,
         sys_lead_lepton_eta_default, sys_sublead_lepton_eta_default, 
         sys_lead_photon_eta_scaledn});
  }

  NamedFunc ggfbdt2503_score_phresup(const vector<FastForest> &xgb_bdts) {
     return XGBoostBDTScore(xgb_bdts, {sys_lead_photon_idmva_resup, 
         sys_lead_photon_relpterr_resup, sys_lead_photon_drmax_resup, 
         sys_lead_photon_drmin_resup, 
         sys_llphoton_refit_pt_phresup/sys_llphoton_refit_m_phresup, 
         sys_llphoton_refit_cosTheta_phresup, 
         sys_llphoton_refit_costheta_phresup, 
         sys_llphoton_refit_psi_phresup,
         sys_lead_lepton_eta_default, sys_sublead_lepton_eta_default, 
         sys_lead_photon_eta_resup});
  }

  NamedFunc ggfbdt2503_score_phresdn(const vector<FastForest> &xgb_bdts) {
     return XGBoostBDTScore(xgb_bdts, {sys_lead_photon_idmva_resup, 
         sys_lead_photon_relpterr_resup, sys_lead_photon_drmax_resup, 
         sys_lead_photon_drmin_resup, 
         sys_llphoton_refit_pt_phresdn/sys_llphoton_refit_m_phresdn, 
         sys_llphoton_refit_cosTheta_phresdn, 
         sys_llphoton_refit_costheta_phresdn, 
         sys_llphoton_refit_psi_phresdn,
         sys_lead_lepton_eta_default, sys_sublead_lepton_eta_default, 
         sys_lead_photon_eta_resdn});
  }

  //new (25-08) BDTs

}
