#include <memory>
#include <string>
#include <vector>

#include "TLorentzVector.h"

#include "core/baby.hpp"
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
using NamedFuncUtilities::FilterNamedFunc;
using NamedFuncUtilities::MultiReduceNamedFunc;
using NamedFuncUtilities::ReduceNamedFunc;
using NamedFuncUtilities::reduce_sum;
using NamedFuncUtilities::reduce_max;
using NamedFuncUtilities::reduce_maxfirst;
using NamedFuncUtilities::reduce_sublead;
using NamedFuncUtilities::reduce_subleadfirst;
using ZgUtilities::get_lep_custom_refit;

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

  //for reference, electrons and muons failing eta, dxy, or dz cuts are dropped from pico lists

  //el_sig with electron scale variation up
  const NamedFunc sys_el_sig_scaleup = NamedFunc("(sys_el_pt_scaleup>7)&&el_idLoose").Name("sys_el_sig_scaleup").EnableCaching(true);

  //el_sig with electron scale variation down
  const NamedFunc sys_el_sig_scaledn = NamedFunc("(sys_el_pt_scaledn>7)&&el_idLoose").Name("sys_el_sig_scaledn").EnableCaching(true);

  //el_sig with electron resolution variation up
  const NamedFunc sys_el_sig_resup = NamedFunc("(sys_el_pt_resup>7)&&el_idLoose").Name("sys_el_sig_resup").EnableCaching(true);

  //el_sig with electron resolution variation down
  const NamedFunc sys_el_sig_resdn = NamedFunc("(sys_el_pt_resdn>7)&&el_idLoose").Name("sys_el_sig_resdn").EnableCaching(true);

  //mu_sig with muon systematics up
  const NamedFunc sys_mu_sig_up = NamedFunc("mu_id&&(mu_reliso<0.35)&&(mu_sip3d<4)&&((mu_pt+mu_ptErr)>5)").Name("sys_mu_sig_up");

  //mu_sig with muon systematics down
  const NamedFunc sys_mu_sig_dn = NamedFunc("mu_id&&(mu_reliso<0.35)&&(mu_sip3d<4)&&((mu_pt-mu_ptErr)>5)").Name("sys_mu_sig_dn");

  //nel with electron scale variation up
  const NamedFunc sys_nel_scaleup = ReduceNamedFunc(sys_el_sig_scaleup,reduce_sum).Name("sys_nel_scaleup");

  //nel with electron scale variation down
  const NamedFunc sys_nel_scaledn = ReduceNamedFunc(sys_el_sig_scaledn,reduce_sum).Name("sys_nel_scaledn");

  //nel with electron resolution variation up
  const NamedFunc sys_nel_resup = ReduceNamedFunc(sys_el_sig_resup,reduce_sum).Name("sys_nel_resup");

  //nel with electron resolution variation down
  const NamedFunc sys_nel_resdn = ReduceNamedFunc(sys_el_sig_resdn,reduce_sum).Name("sys_nel_resdn");

  //nmu with muon systematics up
  const NamedFunc sys_nmu_up = ReduceNamedFunc(sys_mu_sig_up,reduce_sum).Name("sys_nmu_up");

  //nmu with muon systematics down
  const NamedFunc sys_nmu_dn = ReduceNamedFunc(sys_mu_sig_dn,reduce_sum).Name("sys_nmu_dn");

  //nlep with electron scale variation up
  const NamedFunc sys_nlep_elscaleup = NamedFunc(sys_nel_scaleup+"nmu").Name("sys_nlep_elscaleup");

  //nlep with electron scale variation up
  const NamedFunc sys_nlep_elscaledn = NamedFunc(sys_nel_scaledn+"nmu").Name("sys_nlep_elscaledn");

  //nlep with electron resolution variation up
  const NamedFunc sys_nlep_elresup = NamedFunc(sys_nel_resup+"nmu").Name("sys_nlep_elresup");

  //nlep with electron resolution variation up
  const NamedFunc sys_nlep_elresdn = NamedFunc(sys_nel_resdn+"nmu").Name("sys_nlep_elresdn");

  //nlep with muon systematics up
  const NamedFunc sys_nlep_muup = NamedFunc("nel"+sys_nmu_up).Name("sys_nlep_muup");

  //nlep with muon systematics down
  const NamedFunc sys_nlep_mudn = NamedFunc("nel"+sys_nmu_dn).Name("sys_nlep_mudn");

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
  //TODO fix with appropriate sys in next production
  NamedFunc sys_nll_muscaleup = assign_variation_nll("el_sig", 
      sys_mu_sig_up, "muscaleup");

  //nll with muon scale variation down
  //TODO fix with appropriate sys in next production
  NamedFunc sys_nll_muscaledn = assign_variation_nll("el_sig", 
      sys_mu_sig_dn, "muscaledn");

  //nll with muon resolution variation up
  //TODO fix with appropriate sys in next production
  NamedFunc sys_nll_muresup = assign_variation_nll("el_sig", 
      sys_mu_sig_up, "muresup");

  //nll with muon resolution variation down
  //TODO fix with appropriate sys in next production
  NamedFunc sys_nll_muresdn = assign_variation_nll("el_sig", 
      sys_mu_sig_dn, "muresdn");

  //leading electron pt with electron scale variation up
  const NamedFunc sys_lead_el_pt_scaleup = ReduceNamedFunc(FilterNamedFunc("sys_el_pt_scaleup",
      sys_el_sig_scaleup),reduce_max).Name("sys_lead_el_pt_scaleup");

  //subleading electron pt with electron scale variation up
  const NamedFunc sys_sublead_el_pt_scaleup = ReduceNamedFunc(FilterNamedFunc("sys_el_pt_scaleup",
      sys_el_sig_scaleup),reduce_sublead).Name("sys_sublead_el_pt_scaleup");

  //leading electron pt with electron scale variation downn
  const NamedFunc sys_lead_el_pt_scaledn = ReduceNamedFunc(FilterNamedFunc("sys_el_pt_scaledn",
      sys_el_sig_scaledn),reduce_max).Name("sys_lead_el_pt_scaledn");

  //subleading electron pt with electron scale variation down
  const NamedFunc sys_sublead_el_pt_scaledn = ReduceNamedFunc(FilterNamedFunc("sys_el_pt_scaledn",
      sys_el_sig_scaledn),reduce_sublead).Name("sys_sublead_el_pt_scaledn");

  //leading electron pt with electron resolution variation up
  const NamedFunc sys_lead_el_pt_resup = ReduceNamedFunc(FilterNamedFunc("sys_el_pt_resup",
      sys_el_sig_resup),reduce_max).Name("sys_lead_el_pt_resup");

  //subleading electron pt with electron resolution variation up
  const NamedFunc sys_sublead_el_pt_resup = ReduceNamedFunc(FilterNamedFunc("sys_el_pt_resup",
      sys_el_sig_resup),reduce_sublead).Name("sys_sublead_el_pt_resup");

  //leading electron pt with electron resolution variation down
  const NamedFunc sys_lead_el_pt_resdn = ReduceNamedFunc(FilterNamedFunc("sys_el_pt_resdn",
      sys_el_sig_resdn),reduce_max).Name("sys_lead_el_pt_resdn");

  //subleading electron pt with electron resolution variation down
  const NamedFunc sys_sublead_el_pt_resdn = ReduceNamedFunc(FilterNamedFunc("sys_el_pt_resdn",
      sys_el_sig_resdn),reduce_sublead).Name("sys_sublead_el_pt_resdn");

  //leading muon pt with muon variation up
  const NamedFunc sys_lead_mu_pt_up = ReduceNamedFunc(FilterNamedFunc(
      "mu_pt+mu_ptErr",sys_mu_sig_up),
      reduce_max).Name("sys_lead_mu_pt_up");

  //subleading muon pt with muon variation up
  const NamedFunc sys_sublead_mu_pt_up = ReduceNamedFunc(FilterNamedFunc(
      "mu_pt+mu_ptErr",sys_mu_sig_up),
      reduce_sublead).Name("sys_sublead_mu_pt_up");

  //leading muon pt with muon variation down
  const NamedFunc sys_lead_mu_pt_dn = ReduceNamedFunc(FilterNamedFunc(
      "mu_pt-mu_ptErr",sys_mu_sig_dn),
      reduce_max).Name("sys_lead_mu_pt_dn");

  //subleading muon pt with muon variation down
  const NamedFunc sys_sublead_mu_pt_dn = ReduceNamedFunc(FilterNamedFunc(
      "mu_pt-mu_ptErr",sys_mu_sig_dn),
      reduce_sublead).Name("sys_sublead_mu_pt_dn");

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
      "el_sig", "mu_pt+mu_ptErr", sys_mu_sig_up, "muscaleup");

  //trigger and pT cuts with muon scale variation down
  //TODO fix with next production
  NamedFunc sys_trig_pt_muscaledn = assign_variation_trig_pt("el_pt",
      "el_sig", "mu_pt-mu_ptErr", sys_mu_sig_dn, "muscaledn");

  //trigger and pT cuts with muon resolution variation up
  //TODO fix with next production
  NamedFunc sys_trig_pt_muresup = assign_variation_trig_pt("el_pt",
      "el_sig", "mu_pt+mu_ptErr", sys_mu_sig_up, "muresup");

  //trigger and pT cuts with muon resolution variation down
  //TODO fix with next production
  NamedFunc sys_trig_pt_muresdn = assign_variation_trig_pt("el_pt",
      "el_sig", "mu_pt-mu_ptErr", sys_mu_sig_dn, "muresdn");

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
      "mu_pt+mu_ptErr", sys_mu_sig_up, false, "muscaleup");

  //dilepton properties with muon scale variation down
  //TODO: fix in next production
  NamedFunc sys_ll_muscaledn = assign_variation_ll("el_pt", "el_sig", 
      "mu_pt-mu_ptErr", sys_mu_sig_dn, false, "muscaledn");

  //dilepton properties with muon resolution variation up
  //TODO: fix in next production
  NamedFunc sys_ll_muresup = assign_variation_ll("el_pt", "el_sig", 
      "mu_pt+mu_ptErr", sys_mu_sig_up, false, "muresup");

  //dilepton properties with muon resolution variation down
  //TODO: fix in next production
  NamedFunc sys_ll_muresdn = assign_variation_ll("el_pt", "el_sig", 
      "mu_pt-mu_ptErr", sys_mu_sig_dn, false, "muresdn");

  //dilepton mass with electron scale variation up
  NamedFunc sys_ll_m_elscaleup = NamedFunc(sys_ll_elscaleup[3]).Name(
      "sys_ll_m_elscaleup");

  //dilepton mass with electron scale variation down
  NamedFunc sys_ll_m_elscaledn = NamedFunc(sys_ll_elscaledn[3]).Name(
      "sys_ll_m_elscaledn");

  //dilepton mass with electron resolution variation up
  NamedFunc sys_ll_m_elresup = NamedFunc(sys_ll_elresup[3]).Name(
      "sys_ll_m_elresup");

  //dilepton mass with electron resolution variation down
  NamedFunc sys_ll_m_elresdn = NamedFunc(sys_ll_elresdn[3]).Name(
      "sys_ll_m_elresdn");

  //dilepton mass with muon scale variation up
  NamedFunc sys_ll_m_muscaleup = NamedFunc(sys_ll_muscaleup[3]).Name(
      "sys_ll_m_muscaleup");

  //dilepton mass with muon scale variation down
  NamedFunc sys_ll_m_muscaledn = NamedFunc(sys_ll_muscaledn[3]).Name(
      "sys_ll_m_muscaledn");
  
  //dilepton mass with muon resolution variation up
  NamedFunc sys_ll_m_muresup = NamedFunc(sys_ll_muresup[3]).Name(
      "sys_ll_m_muresup");

  //dilepton mass with muon resolution variation down
  NamedFunc sys_ll_m_muresdn = NamedFunc(sys_ll_muresdn[3]).Name(
      "sys_ll_m_muresdn");

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

  //Gets NamedFunc that is assigns Higgs four momentum with variation
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

}
