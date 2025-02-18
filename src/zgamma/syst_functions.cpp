#include <vector>

#include "TLorentzVector.h"

#include "core/baby.hpp"
#include "core/named_func.hpp"
#include "core/named_func_utilities.hpp"
#include "zgamma/syst_functions.hpp"
#include "zgamma/zg_functions.hpp"

using NamedFuncUtilities::FilterNamedFunc;
using NamedFuncUtilities::MultiReduceNamedFunc;
using NamedFuncUtilities::ReduceNamedFunc;
using NamedFuncUtilities::reduce_sum;
using NamedFuncUtilities::reduce_max;
using NamedFuncUtilities::reduce_maxfirst;
using NamedFuncUtilities::reduce_sublead;
using NamedFuncUtilities::reduce_subleadfirst;

namespace ZgFunctions {

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

  //for reference, electrons and muons failing eta, dxy, or dz cuts are dropped from pico lists

  //el_sig with electron scale variation up
  const NamedFunc sys_el_sig_scaleup = NamedFunc("(sys_el_pt_scaleup>7)&&el_idLoose").Name("sys_el_sig_scaleup");

  //el_sig with electron scale variation down
  const NamedFunc sys_el_sig_scaledn = NamedFunc("(sys_el_pt_scaledn>7)&&el_idLoose").Name("sys_el_sig_scaledn");

  //el_sig with electron resolution variation up
  const NamedFunc sys_el_sig_resup = NamedFunc("(sys_el_pt_resup>7)&&el_idLoose").Name("sys_el_sig_resup");

  //el_sig with electron resolution variation down
  const NamedFunc sys_el_sig_resdn = NamedFunc("(sys_el_pt_resdn>7)&&el_idLoose").Name("sys_el_sig_resdn");

  //mu_sig with muon systematics up
  const NamedFunc sys_mu_sig_up = NamedFunc("mu_id||(mu_pt>200&&mu_highptid)&&(mu_reliso<0.35)&&(mu_sip3d<4)&&((mu_corrected_pt+mu_corrected_ptErr)>5)").Name("sys_mu_sig_up");

  //mu_sig with muon systematics down
  const NamedFunc sys_mu_sig_dn = NamedFunc("mu_id||(mu_pt>200&&mu_highptid)&&(mu_reliso<0.35)&&(mu_sip3d<4)&&((mu_corrected_pt-mu_corrected_ptErr)>5)").Name("sys_mu_sig_dn");

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
      "mu_corrected_pt+mu_corrected_ptErr",sys_mu_sig_up),
      reduce_max).Name("sys_lead_mu_pt_up");

  //subleading muon pt with muon variation up
  const NamedFunc sys_sublead_mu_pt_up = ReduceNamedFunc(FilterNamedFunc(
      "mu_corrected_pt+mu_corrected_ptErr",sys_mu_sig_up),
      reduce_sublead).Name("sys_sublead_mu_pt_up");

  //leading muon pt with muon variation down
  const NamedFunc sys_lead_mu_pt_dn = ReduceNamedFunc(FilterNamedFunc(
      "mu_corrected_pt-mu_corrected_ptErr",sys_mu_sig_dn),
      reduce_max).Name("sys_lead_mu_pt_dn");

  //subleading muon pt with muon variation down
  const NamedFunc sys_sublead_mu_pt_dn = ReduceNamedFunc(FilterNamedFunc(
      "mu_corrected_pt-mu_corrected_ptErr",sys_mu_sig_dn),
      reduce_sublead).Name("sys_sublead_mu_pt_dn");

  //for reference, photons failing origin eta cuts are dropped from pico list

  //photon_sig with photon scale variation up
  const NamedFunc sys_photon_sig_scaleup = NamedFunc("(sys_photon_pt_scaleup>15)&&photon_elveto&&(photon_drmin>0.3)&&photon_id80"&&!photon_sigel).Name("sys_photon_sig_elscaleup");

  //photon_sig with photon scale variation down
  const NamedFunc sys_photon_sig_scaledn = NamedFunc("(sys_photon_pt_scaledn>15)&&photon_elveto&&(photon_drmin>0.3)&&photon_id80"&&!photon_sigel).Name("sys_photon_sig_elscaledn");

  //photon_sig with photon resolution variation up
  const NamedFunc sys_photon_sig_resup = NamedFunc("(sys_photon_pt_resup>15)&&photon_elveto&&(photon_drmin>0.3)&&photon_id80"&&!photon_sigel).Name("sys_photon_sig_elresup");

  //photon_sig with photon resolution variation down
  const NamedFunc sys_photon_sig_resdn = NamedFunc("(sys_photon_pt_resdn>15)&&photon_elveto&&(photon_drmin>0.3)&&photon_id80"&&!photon_sigel).Name("sys_photon_sig_elresdn");

  //nphoton with photon scale variation up
  const NamedFunc sys_nphoton_scaleup = ReduceNamedFunc(sys_photon_sig_scaleup,reduce_sum).Name("sys_nphoton_scaleup");

  //nphoton with photon scale variation down
  const NamedFunc sys_nphoton_scaledn = ReduceNamedFunc(sys_photon_sig_scaledn,reduce_sum).Name("sys_nphoton_scaledn");

  //nphoton with photon resolution variation up
  const NamedFunc sys_nphoton_resup = ReduceNamedFunc(sys_photon_sig_resup,reduce_sum).Name("sys_nphoton_resup");

  //nphoton with photon resolution variation down
  const NamedFunc sys_nphoton_resdn = ReduceNamedFunc(sys_photon_sig_resdn,reduce_sum).Name("sys_nphoton_resdn");

  //signal photon pt with photon scale variation up
  const NamedFunc sys_sig_photon_pt_scaleup = FilterNamedFunc("sys_photon_pt_scaleup",
      sys_photon_sig_scaleup).Name("sys_sig_photon_pt_scaleup");

  //signal photon pt with photon scale variation downn
  const NamedFunc sys_sig_photon_pt_scaledn = FilterNamedFunc("sys_photon_pt_scaledn",
      sys_photon_sig_scaledn).Name("sys_sig_photon_pt_scaledn");

  //signal photon pt with photon resolution variation up
  const NamedFunc sys_sig_photon_pt_resup = FilterNamedFunc("sys_photon_pt_resup",
      sys_photon_sig_resup).Name("sys_sig_photon_pt_resup");

  //signal photon pt with photon resolution variation down
  const NamedFunc sys_sig_photon_pt_resdn = FilterNamedFunc("sys_photon_pt_resdn",
      sys_photon_sig_resdn).Name("sys_sig_photon_pt_resdn");

  //leading photon pt with photon scale variation up
  const NamedFunc sys_lead_photon_pt_scaleup = ReduceNamedFunc(sys_sig_photon_pt_scaleup,
      reduce_max).Name("sys_lead_photon_pt_scaleup");

  //leading photon pt with photon scale variation downn
  const NamedFunc sys_lead_photon_pt_scaledn = ReduceNamedFunc(sys_sig_photon_pt_scaledn,
      reduce_max).Name("sys_lead_photon_pt_scaledn");

  //leading photon pt with photon resolution variation up
  const NamedFunc sys_lead_photon_pt_resup = ReduceNamedFunc(sys_sig_photon_pt_resup,
      reduce_max).Name("sys_lead_photon_pt_resup");

  //leading photon pt with photon resolution variation down
  const NamedFunc sys_lead_photon_pt_resdn = ReduceNamedFunc(sys_sig_photon_pt_resdn,
      reduce_max).Name("sys_lead_photon_pt_resdn");

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

  //dilepton mass with electron scale variation up
  const NamedFunc sys_ll_m_elscaleup("sys_ll_m_elscaleup",[](const Baby &b) -> NamedFunc::ScalarType{
    return AssignElectronVariationZ(b, "sys_el_pt_scaleup", sys_el_sig_scaleup).M();
  });

  //dilepton mass with electron scale variation down
  const NamedFunc sys_ll_m_elscaledn("sys_ll_m_elscaledn",[](const Baby &b) -> NamedFunc::ScalarType{
    return AssignElectronVariationZ(b, "sys_el_pt_scaledn", sys_el_sig_scaledn).M();
  });

  //dilepton mass with electron resolution variation up
  const NamedFunc sys_ll_m_elresup("sys_ll_m_elresup",[](const Baby &b) -> NamedFunc::ScalarType{
    return AssignElectronVariationZ(b, "sys_el_pt_resup", sys_el_sig_resup).M();
  });

  //dilepton mass with electron resolution variation down
  const NamedFunc sys_ll_m_elresdn("sys_ll_m_elresdn",[](const Baby &b) -> NamedFunc::ScalarType{
    return AssignElectronVariationZ(b, "sys_el_pt_resdn", sys_el_sig_resdn).M();
  });

  //dilepton mass with muon variation up
  const NamedFunc sys_ll_m_muup("sys_ll_m_muup",[](const Baby &b) -> NamedFunc::ScalarType{
    return AssignMuonVariationZ(b, "mu_corrected_pt+mu_corrected_ptErr", sys_mu_sig_up).M();
  });

  //dilepton mass with muon variation dn
  const NamedFunc sys_ll_m_mudn("sys_ll_m_mudn",[](const Baby &b) -> NamedFunc::ScalarType{
    return AssignMuonVariationZ(b, "mu_corrected_pt-mu_corrected_ptErr", sys_mu_sig_dn).M();
  });

  //Higgs candidate mass with electron scale variation up
  const NamedFunc sys_llphoton_m_elscaleup("sys_llphoton_m_elscaleup",
      [](const Baby &b) -> NamedFunc::ScalarType{
    return AssignElectronVariationH(b, "sys_el_pt_scaleup", sys_el_sig_scaleup).M();
  });

  //Higgs candidate mass with electron scale variation down
  const NamedFunc sys_llphoton_m_elscaledn("sys_llphoton_m_elscaledn",
      [](const Baby &b) -> NamedFunc::ScalarType{
    return AssignElectronVariationH(b, "sys_el_pt_scaledn", sys_el_sig_scaledn).M();
  });

  //Higgs candidate mass with electron resolution variation up
  const NamedFunc sys_llphoton_m_elresup("sys_llphoton_m_elresup",
      [](const Baby &b) -> NamedFunc::ScalarType{
    return AssignElectronVariationH(b, "sys_el_pt_resup", sys_el_sig_resup).M();
  });

  //Higgs candidate mass with electron resolution variation down
  const NamedFunc sys_llphoton_m_elresdn("sys_llphoton_m_elresdn",
      [](const Baby &b) -> NamedFunc::ScalarType{
    return AssignElectronVariationH(b, "sys_el_pt_resdn", sys_el_sig_resdn).M();
  });

  //Higgs candidate mass with muon variation up
  const NamedFunc sys_llphoton_m_muup("sys_llphoton_m_muup",
      [](const Baby &b) -> NamedFunc::ScalarType{
    return AssignMuonVariationH(b, "mu_corrected_pt+mu_corrected_ptErr", sys_mu_sig_up).M();
  });

  //Higgs candidate mass with muon variation dn
  const NamedFunc sys_llphoton_m_mudn("sys_llphoton_m_mudn",
      [](const Baby &b) -> NamedFunc::ScalarType{
    return AssignMuonVariationH(b, "mu_corrected_pt-mu_corrected_ptErr", sys_mu_sig_dn).M();
  });

  //Higgs candidate mass with photon scale variation up
  const NamedFunc sys_llphoton_m_phscaleup("sys_llphoton_m_phscaleup",
      [](const Baby &b) -> NamedFunc::ScalarType{
    return AssignPhotonVariationH(b, "sys_photon_pt_scaleup", sys_photon_sig_scaleup).M();
  });

  //Higgs candidate mass with photon scale variation down
  const NamedFunc sys_llphoton_m_phscaledn("sys_llphoton_m_phscaledn",
      [](const Baby &b) -> NamedFunc::ScalarType{
    return AssignPhotonVariationH(b, "sys_photon_pt_scaledn", sys_photon_sig_scaledn).M();
  });

  //Higgs candidate mass with photon resolution variation up
  const NamedFunc sys_llphoton_m_phresup("sys_llphoton_m_phresup",
      [](const Baby &b) -> NamedFunc::ScalarType{
    return AssignPhotonVariationH(b, "sys_photon_pt_resup", sys_photon_sig_resup).M();
  });

  //Higgs candidate mass with photon resolution variation down
  const NamedFunc sys_llphoton_m_phresdn("sys_llphoton_m_phresdn",
      [](const Baby &b) -> NamedFunc::ScalarType{
    return AssignPhotonVariationH(b, "sys_photon_pt_resdn", sys_photon_sig_resdn).M();
  });

  //helper function to get Z candidate four vector after varying electron energy
  TLorentzVector AssignElectronVariationZ(const Baby &b, const NamedFunc &el_pt, 
                                          const NamedFunc &el_sig) {
    TLorentzVector el1, el2, zel, z;
    float min_dm = 999;
    z.SetPtEtaPhiM(0.0,0.0,0.0,-999.0);
    for (unsigned ill = 0; ill<b.ll_pt()->size(); ill++) {
      if (b.ll_lepid()->at(ill) == 13) {
        float dm = fabs(b.ll_m()->at(ill)-91.1876);
        if (dm < min_dm) {
          min_dm = dm;
          z.SetPtEtaPhiM(b.ll_pt()->at(ill), b.ll_eta()->at(ill), b.ll_phi()->at(ill), b.ll_m()->at(ill));
        }
      }
    }
    std::vector<double> b_el_sig = el_sig.GetVector(b);
    for (unsigned iel1 = 0; iel1<b.el_pt()->size(); iel1++) {
      if (!b_el_sig.at(iel1)) continue;
      for (unsigned iel2 = iel1+1; iel2<b.el_pt()->size(); iel2++) {
        if (!b_el_sig.at(iel2)) continue;
        std::vector<double> b_el_pt = el_pt.GetVector(b);
        el1.SetPtEtaPhiM(b_el_pt.at(iel1), b.el_eta()->at(iel1), 
                         b.el_phi()->at(iel1), 0.000511);
        el2.SetPtEtaPhiM(b_el_pt.at(iel2), b.el_eta()->at(iel2), 
                         b.el_phi()->at(iel2), 0.000511);
        zel = el1 + el2;
        float dm = fabs(zel.M()-91.1876);
        if (dm < min_dm) {
          z = zel;
          min_dm = dm;
        }
      }
    }
    return z;
  }

  //helper function to get Z candidate four vector after varying muon energy
  TLorentzVector AssignMuonVariationZ(const Baby &b, const NamedFunc &mu_pt, 
                                      const NamedFunc &mu_sig) {
    TLorentzVector mu1, mu2, zmu, z;
    float min_dm = 999;
    z.SetPtEtaPhiM(0.0,0.0,0.0,-999.0);
    for (unsigned ill = 0; ill<b.ll_pt()->size(); ill++) {
      if (b.ll_lepid()->at(ill) == 11) {
        float dm = fabs(b.ll_m()->at(ill)-91.1876);
        if (dm < min_dm) {
          min_dm = dm;
          z.SetPtEtaPhiM(b.ll_pt()->at(ill), b.ll_eta()->at(ill), b.ll_phi()->at(ill), b.ll_m()->at(ill));
        }
      }
    }
    std::vector<double> b_mu_sig = mu_sig.GetVector(b);
    for (unsigned imu1 = 0; imu1<b.mu_pt()->size(); imu1++) {
      if (!b_mu_sig.at(imu1)) continue;
      for (unsigned imu2 = imu1+1; imu2<b.mu_pt()->size(); imu2++) {
        if (!b_mu_sig.at(imu2)) continue;
        std::vector<double> b_mu_pt = mu_pt.GetVector(b);
        mu1.SetPtEtaPhiM(b_mu_pt.at(imu1), b.mu_eta()->at(imu1), 
                         b.mu_phi()->at(imu1), 0.10566);
        mu2.SetPtEtaPhiM(b_mu_pt.at(imu2), b.mu_eta()->at(imu2), 
                         b.mu_phi()->at(imu2), 0.10566);
        zmu = mu1 + mu2;
        float dm = fabs(zmu.M()-91.1876);
        if (dm < min_dm) {
          z = zmu;
          min_dm = dm;
        }
      }
    }
    return z;
  }

  //helper function to get H candidate four vector after varying electron energy
  TLorentzVector AssignElectronVariationH(const Baby &b, const NamedFunc &el_pt, 
                                          const NamedFunc &el_sig) {
    TLorentzVector ph, z;
    z = AssignElectronVariationZ(b,el_pt,el_sig);
    ph.SetPtEtaPhiM(b.photon_pt()->at(0),b.photon_eta()->at(0),b.photon_phi()->at(0),0);
    return (ph+z);
  }

  //helper function to get H candidate four vector after varying electron energy
  TLorentzVector AssignMuonVariationH(const Baby &b, const NamedFunc &mu_pt, 
                                      const NamedFunc &mu_sig) {
    TLorentzVector ph, z;
    z = AssignMuonVariationZ(b,mu_pt,mu_sig);
    ph.SetPtEtaPhiM(b.photon_pt()->at(0),b.photon_eta()->at(0),b.photon_phi()->at(0),0);
    return (ph+z);
  }

  //helper function to get H candidate four vector after varying photon energy
  TLorentzVector AssignPhotonVariationH(const Baby &b, const NamedFunc &photon_pt, 
                                        const NamedFunc &photon_sig) {
    TLorentzVector ph, z;
    float max_ph_pt = -999;
    std::vector<double> b_photon_pt = photon_pt.GetVector(b);
    std::vector<double> b_photon_sig = photon_sig.GetVector(b);
    z.SetPtEtaPhiM(b.ll_pt()->at(0),b.ll_eta()->at(0),b.ll_phi()->at(0),b.ll_m()->at(0));
    for (unsigned iph = 0; iph < b.photon_pt()->size(); iph++) {
      if (!b_photon_sig.at(iph)) continue;
      if (b_photon_pt.at(iph) > max_ph_pt) {
        ph.SetPtEtaPhiM(b_photon_pt.at(iph),b.photon_eta()->at(iph),b.photon_phi()->at(iph),0);
        max_ph_pt = b_photon_pt.at(iph);
      }
    }
    return (ph+z);
  }

}
