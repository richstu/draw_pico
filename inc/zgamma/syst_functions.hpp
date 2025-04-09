#ifndef H_SYST_FUNCTIONS
#define H_SYST_FUNCTIONS

#include "TLorentzVector.h"

#include "core/baby.hpp"
#include "core/named_func.hpp"

namespace ZgFunctions {

  //weight implementing variations in alphaS
  extern const NamedFunc sys_w_alphas;

  //weight implementing variations in PDFs for ggF
  extern const NamedFunc sys_w_pdf_ggf;

  //weight implementing variations in PDFs for VBF/WH/ZH
  extern const NamedFunc sys_w_pdf_qq;

  //weight implementing variations in PDFs for ttH
  extern const NamedFunc sys_w_pdf_tth;

  //weight implementing variations in ggF cross section
  extern const NamedFunc sys_w_ggf_xs;

  //weight implementing up variation in VBF cross section
  extern const NamedFunc sys_w_vbf_xs_up;

  //weight implementing down variation in VBF cross section
  extern const NamedFunc sys_w_vbf_xs_dn;

  //weight implementing up variation in WH cross section
  extern const NamedFunc sys_w_wh_xs_up;

  //weight implementing down variation in WH cross section
  extern const NamedFunc sys_w_wh_xs_dn;

  //weight implementing up variation in ZH cross section
  extern const NamedFunc sys_w_zh_xs_up;

  //weight implementing down variation in ZH cross section
  extern const NamedFunc sys_w_zh_xs_dn;

  //weight implementing up variation in ttH cross section
  extern const NamedFunc sys_w_tth_xs_up;

  //weight implementing down variation in ttH cross section
  extern const NamedFunc sys_w_tth_xs_dn;

  //weight implementing variation in H->Zgamma BR
  extern const NamedFunc sys_w_htozg_br;

  //weight implementing variation in H->MuMu BR
  extern const NamedFunc sys_w_htomumu_br;

  //weight implementing variation in mq (mostly mtop)
  extern const NamedFunc sys_w_mq;

  //weight implementing variation in run 2 lumi
  extern const NamedFunc sys_w_lumi_run2;

  //weight implementing variation in 2022 lumi
  extern const NamedFunc sys_w_lumi_2022;

  //weight implementing variation in 2023 lumi
  extern const NamedFunc sys_w_lumi_2023;

  //weight implementing upward variation in pileup weights 
  extern const NamedFunc sys_w_pu_up;

  //weight implementing downward variation in pileup weights 
  extern const NamedFunc sys_w_pu_dn;

  //weight implementing upward variation in prefire weights 
  extern const NamedFunc sys_w_prefire_up;

  //weight implementing downward variation in prefire weights 
  extern const NamedFunc sys_w_prefire_dn;

  //weight implementing upward variation in b/c-tagging weights 
  extern const NamedFunc sys_w_bctag_up;

  //weight implementing downward variation in b/c-tagging weights 
  extern const NamedFunc sys_w_bctag_dn;

  //weight implementing upward variation in udsg-mistagging weights 
  extern const NamedFunc sys_w_udsgtag_up;

  //weight implementing downward variation in udsg-mistagging weights 
  extern const NamedFunc sys_w_udsgtag_dn;

  //weight implementing upward variation in electron weights (efficiency)
  extern const NamedFunc sys_w_el_up;

  //weight implementing downward variation in electron weights (efficiency)
  extern const NamedFunc sys_w_el_dn;

  //weight implementing upward variation in muon weights (efficiency)
  extern const NamedFunc sys_w_mu_up;

  //weight implementing downward variation in muon weights (efficiency)
  extern const NamedFunc sys_w_mu_dn;
  
  //weight implementing upward variation in photon weights (efficiency)
  extern const NamedFunc sys_w_photon_up;

  //weight implementing downward variation in photon weights (efficiency)
  extern const NamedFunc sys_w_photon_dn;

  //el_sig with electron scale variation up
  extern const NamedFunc sys_el_sig_scaleup;

  //el_sig with electron scale variation down
  extern const NamedFunc sys_el_sig_scaledn;

  //el_sig with electron resolution variation up
  extern const NamedFunc sys_el_sig_resup;

  //el_sig with electron resolution variation down
  extern const NamedFunc sys_el_sig_resdn;

  //mu_sig with muon systematics up
  extern const NamedFunc sys_mu_sig_up;

  //mu_sig with muon systematics down
  extern const NamedFunc sys_mu_sig_dn;

  //nel with electron scale variation up
  extern const NamedFunc sys_nel_scaleup;

  //nel with electron scale variation down
  extern const NamedFunc sys_nel_scaledn;

  //nel with electron resolution variation up
  extern const NamedFunc sys_nel_resup;

  //nel with electron resolution variation down
  extern const NamedFunc sys_nel_resdn;

  //nmu with muon systematics up
  extern const NamedFunc sys_nmu_up;

  //nmu with muon systematics down
  extern const NamedFunc sys_nmu_dn;

  //nlep with electron scale variation up
  extern const NamedFunc sys_nlep_elscaleup;

  //nlep with electron scale variation up
  extern const NamedFunc sys_nlep_elscaledn;

  //nlep with electron resolution variation up
  extern const NamedFunc sys_nlep_elresup;

  //nlep with electron resolution variation up
  extern const NamedFunc sys_nlep_elresdn;

  //nlep with muon systematics up
  extern const NamedFunc sys_nlep_muup;

  //nlep with muon systematics down
  extern const NamedFunc sys_nlep_mudn;

  //leading electron pt with electron scale variation up
  extern const NamedFunc sys_lead_el_pt_scaleup;

  //subleading electron pt with electron scale variation up
  extern const NamedFunc sys_sublead_el_pt_scaleup;

  //leading electron pt with electron scale variation downn
  extern const NamedFunc sys_lead_el_pt_scaledn;

  //subleading electron pt with electron scale variation down
  extern const NamedFunc sys_sublead_el_pt_scaledn;

  //leading electron pt with electron resolution variation up
  extern const NamedFunc sys_lead_el_pt_resup;

  //subleading electron pt with electron resolution variation up
  extern const NamedFunc sys_sublead_el_pt_resup;

  //leading electron pt with electron resolution variation down
  extern const NamedFunc sys_lead_el_pt_resdn;

  //subleading electron pt with electron resolution variation down
  extern const NamedFunc sys_sublead_el_pt_resdn;

  //leading muon pt with muon variation up
  extern const NamedFunc sys_lead_mu_pt_up;

  //subleading muon pt with muon variation up
  extern const NamedFunc sys_sublead_mu_pt_up;

  //leading muon pt with muon variation down
  extern const NamedFunc sys_lead_mu_pt_dn;

  //subleading muon pt with muon variation down
  extern const NamedFunc sys_sublead_mu_pt_dn;

  //photon_sig with photon scale variation down
  extern const NamedFunc sys_photon_sig_scaleup;

  //photon_sig with photon scale variation down
  extern const NamedFunc sys_photon_sig_scaledn;

  //photon_sig with photon resolution variation down
  extern const NamedFunc sys_photon_sig_resup;

  //photon_sig with photon resolution variation down
  extern const NamedFunc sys_photon_sig_resdn;

  //nphoton with photon scale variation up
  extern const NamedFunc sys_nphoton_scaleup;

  //nphoton with photon scale variation down
  extern const NamedFunc sys_nphoton_scaledn;

  //nphoton with photon resolution variation up
  extern const NamedFunc sys_nphoton_resup;

  //nphoton with photon resolution variation down
  extern const NamedFunc sys_nphoton_resdn;

  //signal photon pt with photon scale variation up
  extern const NamedFunc sys_sig_photon_pt_scaleup;

  //signal photon pt with photon scale variation downn
  extern const NamedFunc sys_sig_photon_pt_scaledn;

  //signal photon pt with photon resolution variation up
  extern const NamedFunc sys_sig_photon_pt_resup;

  //signal photon pt with photon resolution variation down
  extern const NamedFunc sys_sig_photon_pt_resdn;

  //leading photon pt with photon scale variation up
  extern const NamedFunc sys_lead_photon_pt_scaleup;

  //leading photon pt with photon scale variation downn
  extern const NamedFunc sys_lead_photon_pt_scaledn;

  //leading photon pt with photon resolution variation up
  extern const NamedFunc sys_lead_photon_pt_resup;

  //leading photon pt with photon resolution variation down
  extern const NamedFunc sys_lead_photon_pt_resdn;

  //signal photon ID80 flag with photon scale variation up
  extern const NamedFunc sys_sig_photon_id80_scaleup;

  //signal photon ID80 flag with photon scale variation downn
  extern const NamedFunc sys_sig_photon_id80_scaledn;

  //signal photon ID80 flag with photon resolution variation up
  extern const NamedFunc sys_sig_photon_id80_resup;

  //signal photon ID80 flag with photon resolution variation down
  extern const NamedFunc sys_sig_photon_id80_resdn;

  //leading photon ID80 flag with photon scale variation up
  extern const NamedFunc sys_lead_photon_id80_scaleup;

  //leading photon ID80 flag with photon scale variation down
  extern const NamedFunc sys_lead_photon_id80_scaledn;

  //leading photon ID80 flag with photon resolution variation up
  extern const NamedFunc sys_lead_photon_id80_resup;

  //leading photon ID80 flag with photon resolution variation down
  extern const NamedFunc sys_lead_photon_id80_resdn;

  //dilepton mass with electron scale variation up
  extern const NamedFunc sys_ll_m_elscaleup;

  //dilepton mass with electron scale variation down
  extern const NamedFunc sys_ll_m_elscaledn;

  //dilepton mass with electron resolution variation up
  extern const NamedFunc sys_ll_m_elresup;

  //dilepton mass with electron resolution variation down
  extern const NamedFunc sys_ll_m_elresdn;

  //dilepton mass with muon variation up
  extern const NamedFunc sys_ll_m_muup;

  //dilepton mass with muon variation dn
  extern const NamedFunc sys_ll_m_mudn;

  //Higgs candidate mass with electron scale variation up
  extern const NamedFunc sys_llphoton_m_elscaleup;

  //Higgs candidate mass with electron scale variation down
  extern const NamedFunc sys_llphoton_m_elscaledn;

  //Higgs candidate mass with electron resolution variation up
  extern const NamedFunc sys_llphoton_m_elresup;

  //Higgs candidate mass with electron resolution variation down
  extern const NamedFunc sys_llphoton_m_elresdn;

  //Higgs candidate mass with muon variation up
  extern const NamedFunc sys_llphoton_m_muup;

  //Higgs candidate mass with muon variation dn
  extern const NamedFunc sys_llphoton_m_mudn;

  //Higgs candidate mass with photon scale variation up
  extern const NamedFunc sys_llphoton_m_phscaleup;

  //Higgs candidate mass with photon scale variation down
  extern const NamedFunc sys_llphoton_m_phscaledn;

  //Higgs candidate mass with photon resolution variation up
  extern const NamedFunc sys_llphoton_m_phresup;

  //Higgs candidate mass with photon resolution variation down
  extern const NamedFunc sys_llphoton_m_phresdn;

  //helper function to get Z candidate four vector after varying electron energy
  TLorentzVector AssignElectronVariationZ(const Baby &b, const NamedFunc &el_pt, 
                                          const NamedFunc &el_sig);

  //helper function to get Z candidate four vector after varying muon energy
  TLorentzVector AssignMuonVariationZ(const Baby &b, const NamedFunc &mu_pt, 
                                      const NamedFunc &mu_sig);

  //helper function to get H candidate four vector after varying electron energy
  TLorentzVector AssignElectronVariationH(const Baby &b, const NamedFunc &el_pt, 
                                          const NamedFunc &el_sig);

  //helper function to get H candidate four vector after varying electron energy
  TLorentzVector AssignMuonVariationH(const Baby &b, const NamedFunc &mu_pt, 
                                      const NamedFunc &mu_sig);

  //helper function to get H candidate four vector after varying photon energy
  TLorentzVector AssignPhotonVariationH(const Baby &b, const NamedFunc &photon_pt, 
                                        const NamedFunc &photon_sig);

  //helper function to get H candidate four vector after variation
  //doesn't consider edge cases ex. if selected Z candidate changes after
  //variation
  NamedFunc SimpleAssignVariationH(const NamedFunc &el_pt, 
                                   const NamedFunc &mu_pt,
                                   const NamedFunc &photon_pt,
                                   const std::string &name);

}

#endif
