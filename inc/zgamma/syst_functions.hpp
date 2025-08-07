#ifndef H_SYST_FUNCTIONS
#define H_SYST_FUNCTIONS

#include <memory>
#include <vector>

#include "TLorentzVector.h"

#include "core/baby.hpp"
#include "core/named_func.hpp"

#include "zgamma/KinZfitter.hpp"

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

  //mu_sig with muon scale variation up
  extern const NamedFunc sys_mu_sig_scaleup;

  //mu_sig with muon scale variation down
  extern const NamedFunc sys_mu_sig_scaledn;

  //mu_sig with muon resolution variation up
  extern const NamedFunc sys_mu_sig_resup;

  //mu_sig with muon resolution variation down
  extern const NamedFunc sys_mu_sig_resdn;

  //nel with electron scale variation up
  extern const NamedFunc sys_nel_scaleup;

  //nel with electron scale variation down
  extern const NamedFunc sys_nel_scaledn;

  //nel with electron resolution variation up
  extern const NamedFunc sys_nel_resup;

  //nel with electron resolution variation down
  extern const NamedFunc sys_nel_resdn;

  //nmu with muon scale variation up
  extern const NamedFunc sys_nmu_scaleup;

  //nmu with muon scale variation down
  extern const NamedFunc sys_nmu_scaledn;

  //nmu with muon resolution variation up
  extern const NamedFunc sys_nmu_resup;

  //nmu with muon resolution variation down
  extern const NamedFunc sys_nmu_resdn;

  //nlep with electron scale variation up
  extern const NamedFunc sys_nlep_elscaleup;

  //nlep with electron scale variation up
  extern const NamedFunc sys_nlep_elscaledn;

  //nlep with electron resolution variation up
  extern const NamedFunc sys_nlep_elresup;

  //nlep with electron resolution variation up
  extern const NamedFunc sys_nlep_elresdn;

  //nlep with muon scale variation up
  extern const NamedFunc sys_nlep_muscaleup;

  //nlep with muon scale variation down
  extern const NamedFunc sys_nlep_muscaledn;

  //nlep with muon resolution variation up
  extern const NamedFunc sys_nlep_muresup;

  //nlep with muon resolution variation down
  extern const NamedFunc sys_nlep_muresdn;

  //Gets NamedFunc that is number of ll candidates with variation
  NamedFunc assign_variation_nll(const NamedFunc &variation_el_sig, 
                                 const NamedFunc &variation_mu_sig,
                                 const std::string &name);

  //nll with electron scale variation up
  extern NamedFunc sys_nll_elscaleup;

  //nll with electron scale variation down
  extern NamedFunc sys_nll_elscaledn;

  //nll with electron resolution variation up
  extern NamedFunc sys_nll_elresup;

  //nll with electron resolution variation down
  extern NamedFunc sys_nll_elresdn;

  //nll with muon scale variation up
  extern NamedFunc sys_nll_muscaleup;

  //nll with muon scale variation down
  extern NamedFunc sys_nll_muscaledn;

  //nll with muon resolution variation up
  extern NamedFunc sys_nll_muresup;

  //nll with muon resolution variation down
  extern NamedFunc sys_nll_muresdn;

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

  //leading muon pt with muon scale variation up
  extern const NamedFunc sys_lead_mu_pt_scaleup;

  //subleading muon pt with muon scale variation up
  extern const NamedFunc sys_sublead_mu_pt_scaleup;

  //leading muon pt with muon scale variation down
  extern const NamedFunc sys_lead_mu_pt_scaledn;

  //subleading muon pt with muon scale variation down
  extern const NamedFunc sys_sublead_mu_pt_scaledn;

  //leading muon pt with muon resolution variation up
  extern const NamedFunc sys_lead_mu_pt_resup;

  //subleading muon pt with muon resolution variation up
  extern const NamedFunc sys_sublead_mu_pt_resup;

  //leading muon pt with muon resolution variation down
  extern const NamedFunc sys_lead_mu_pt_resdn;

  //subleading muon pt with muon resolution variation down
  extern const NamedFunc sys_sublead_mu_pt_resdn;

  //Gets NamedFunc that is trigger+pT cut flag with variation
  NamedFunc assign_variation_trig_pt(const NamedFunc &variation_el_pt, 
                                     const NamedFunc &variation_el_sig,
                                     const NamedFunc &variation_mu_pt,
                                     const NamedFunc &variation_mu_sig,
                                     const std::string &name);

  //trigger and pT cuts with electron scale variation up
  extern NamedFunc sys_trig_pt_elscaleup;

  //trigger and pT cuts with electron scale variation down
  extern NamedFunc sys_trig_pt_elscaledn;

  //trigger and pT cuts with electron resolution variation up
  extern NamedFunc sys_trig_pt_elresup;

  //trigger and pT cuts with electron resolution variation down
  extern NamedFunc sys_trig_pt_elresdn;

  //trigger and pT cuts with muon scale variation up
  extern NamedFunc sys_trig_pt_muscaleup;

  //trigger and pT cuts with muon scale variation down
  extern NamedFunc sys_trig_pt_muscaledn;

  //trigger and pT cuts with muon resolution variation up
  extern NamedFunc sys_trig_pt_muresup;

  //trigger and pT cuts with muon resolution variation down
  extern NamedFunc sys_trig_pt_muresdn;

  //Gets NamedFunc that is max lep miniso with variation
  NamedFunc assign_variation_max_lep_miniso(const NamedFunc &variation_el_sig,
                                            const NamedFunc &variation_mu_sig,
                                            const std::string &name);

  //max lep miniso with electron scale variation up
  extern NamedFunc sys_max_lep_miniso_elscaleup;

  //max lep miniso with electron scale variation down
  extern NamedFunc sys_max_lep_miniso_elscaledn;

  //max lep miniso with electron resolution variation up
  extern NamedFunc sys_max_lep_miniso_elresup;

  //max lep miniso with electron resolution variation down
  extern NamedFunc sys_max_lep_miniso_elresdn;

  //max lep miniso with muon scale variation up
  extern NamedFunc sys_max_lep_miniso_muscaleup;

  //max lep miniso with muon scale variation down
  extern NamedFunc sys_max_lep_miniso_muscaledn;

  //max lep miniso with muon resolution variation up
  extern NamedFunc sys_max_lep_miniso_muresup;

  //max lep miniso with muon resolution variation down
  extern NamedFunc sys_max_lep_miniso_muresdn;

  //photon_sig with photon scale variation down
  extern NamedFunc sys_photon_sig_scaleup;

  //photon_sig with photon scale variation down
  extern NamedFunc sys_photon_sig_scaledn;

  //photon_sig with photon resolution variation down
  extern NamedFunc sys_photon_sig_resup;

  //photon_sig with photon resolution variation down
  extern NamedFunc sys_photon_sig_resdn;

  //nphoton with photon scale variation up
  extern NamedFunc sys_nphoton_scaleup;

  //nphoton with photon scale variation down
  extern NamedFunc sys_nphoton_scaledn;

  //nphoton with photon resolution variation up
  extern NamedFunc sys_nphoton_resup;

  //nphoton with photon resolution variation down
  extern NamedFunc sys_nphoton_resdn;

  //signal photon pt with photon scale variation up
  extern NamedFunc sys_sig_photon_pt_scaleup;

  //signal photon pt with photon scale variation downn
  extern NamedFunc sys_sig_photon_pt_scaledn;

  //signal photon pt with photon resolution variation up
  extern NamedFunc sys_sig_photon_pt_resup;

  //signal photon pt with photon resolution variation down
  extern NamedFunc sys_sig_photon_pt_resdn;

  //signal photon eta with photon scale variation up
  extern NamedFunc sys_sig_photon_eta_scaleup;

  //signal photon eta with photon scale variation downn
  extern NamedFunc sys_sig_photon_eta_scaledn;

  //signal photon eta with photon resolution variation up
  extern NamedFunc sys_sig_photon_eta_resup;

  //signal photon eta with photon resolution variation down
  extern NamedFunc sys_sig_photon_eta_resdn;

  //signal photon phi with photon scale variation up
  extern NamedFunc sys_sig_photon_phi_scaleup;

  //signal photon phi with photon scale variation downn
  extern NamedFunc sys_sig_photon_phi_scaledn;

  //signal photon phi with photon resolution variation up
  extern NamedFunc sys_sig_photon_phi_resup;

  //signal photon phi with photon resolution variation down
  extern NamedFunc sys_sig_photon_phi_resdn;

  //leading photon pt with photon scale variation up
  extern NamedFunc sys_lead_photon_pt_scaleup;

  //leading photon pt with photon scale variation downn
  extern NamedFunc sys_lead_photon_pt_scaledn;

  //leading photon pt with photon resolution variation up
  extern NamedFunc sys_lead_photon_pt_resup;

  //leading photon pt with photon resolution variation down
  extern NamedFunc sys_lead_photon_pt_resdn;

  //leading photon eta with photon scale variation up
  extern NamedFunc sys_lead_photon_eta_scaleup;

  //leading photon eta with photon scale variation downn
  extern NamedFunc sys_lead_photon_eta_scaledn;

  //leading photon eta with photon resolution variation up
  extern NamedFunc sys_lead_photon_eta_resup;

  //leading photon eta with photon resolution variation down
  extern NamedFunc sys_lead_photon_eta_resdn;

  //leading photon phi with photon scale variation up
  extern NamedFunc sys_lead_photon_phi_scaleup;

  //leading photon phi with photon scale variation downn
  extern NamedFunc sys_lead_photon_phi_scaledn;

  //leading photon phi with photon resolution variation up
  extern NamedFunc sys_lead_photon_phi_resup;

  //leading photon phi with photon resolution variation down
  extern NamedFunc sys_lead_photon_phi_resdn;

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

  //get Z candidate properties with variation
  //returns (pt, eta, phi, m, lepid, i1, i2, idx)
  NamedFunc assign_variation_ll(const NamedFunc &el_pt, 
                                const NamedFunc &el_sig,
                                const NamedFunc &mu_pt, 
                                const NamedFunc &mu_sig,
                                bool is_el,
                                const std::string &name);

  //default dilepton properties
  extern NamedFunc sys_ll_default;

  //dilepton properties with electron scale variation up
  extern NamedFunc sys_ll_elscaleup;

  //dilepton properties with electron scale variation down
  extern NamedFunc sys_ll_elscaledn;

  //dilepton properties with electron resolution variation up
  extern NamedFunc sys_ll_elresup;

  //dilepton properties with electron resolution variation down
  extern NamedFunc sys_ll_elresdn;

  //dilepton properties with muon scale variation up
  extern NamedFunc sys_ll_muscaleup;

  //dilepton properties with muon scale variation dn
  extern NamedFunc sys_ll_muscaledn;

  //dilepton properties with muon resolution variation up
  extern NamedFunc sys_ll_muresup;

  //dilepton properties with muon resolution variation dn
  extern NamedFunc sys_ll_muresdn;

  //dilepton mass with electron scale variation up
  extern NamedFunc sys_ll_m_elscaleup;

  //dilepton mass with electron scale variation down
  extern NamedFunc sys_ll_m_elscaledn;

  //dilepton mass with electron resolution variation up
  extern NamedFunc sys_ll_m_elresup;

  //dilepton mass with electron resolution variation down
  extern NamedFunc sys_ll_m_elresdn;

  //dilepton mass with muon scale variation up
  extern NamedFunc sys_ll_m_muscaleup;

  //dilepton mass with muon scale variation dn
  extern NamedFunc sys_ll_m_muscaledn;

  //dilepton mass with muon resolution variation up
  extern NamedFunc sys_ll_m_muresup;

  //dilepton mass with muon resolution variation dn
  extern NamedFunc sys_ll_m_muresdn;

  //get H candidate properties with variation
  //returns (pt, eta, phi, m)
  NamedFunc assign_variation_llphoton_p4(const NamedFunc &ll_p4, 
                                         const NamedFunc &lead_photon_pt,
                                         const NamedFunc &lead_photon_eta,
                                         const NamedFunc &lead_photon_phi,
                                         const int var_pdgid,
                                         const std::string &name);

  //H four momentum with electron scale variation up 
  extern NamedFunc sys_llphoton_p4_elscaleup;

  //H four momentum with electron scale variation down
  extern NamedFunc sys_llphoton_p4_elscaledn;

  //H four momentum with electron resolution variation up 
  extern NamedFunc sys_llphoton_p4_elresup;

  //H four momentum with electron resolution variation down
  extern NamedFunc sys_llphoton_p4_elresdn;

  //H four momentum with muon scale variation up 
  extern NamedFunc sys_llphoton_p4_muscaleup;

  //H four momentum with muon scale variation down
  extern NamedFunc sys_llphoton_p4_muscaledn;

  //H four momentum with muon resolution variation up 
  extern NamedFunc sys_llphoton_p4_muresup;

  //H four momentum with muon resolution variation down
  extern NamedFunc sys_llphoton_p4_muresdn;

  //H four momentum with photon scale variation up 
  extern NamedFunc sys_llphoton_p4_phscaleup;

  //H four momentum with photon scale variation down
  extern NamedFunc sys_llphoton_p4_phscaledn;

  //H four momentum with photon resolution variation up 
  extern NamedFunc sys_llphoton_p4_phresup;

  //H four momentum with photon resolution variation down
  extern NamedFunc sys_llphoton_p4_phresdn;

  //Higgs candidate mass with electron scale variation up
  extern NamedFunc sys_llphoton_m_elscaleup;

  //Higgs candidate mass with electron scale variation down
  extern NamedFunc sys_llphoton_m_elscaledn;

  //Higgs candidate mass with electron resolution variation up
  extern NamedFunc sys_llphoton_m_elresup;

  //Higgs candidate mass with electron resolution variation down
  extern NamedFunc sys_llphoton_m_elresdn;

  //Higgs candidate mass with muon scale variation up
  extern NamedFunc sys_llphoton_m_muscaleup;

  //Higgs candidate mass with muon scale variation dn
  extern NamedFunc sys_llphoton_m_muscaledn;

  //Higgs candidate mass with muon resolution variation up
  extern NamedFunc sys_llphoton_m_muresup;

  //Higgs candidate mass with muon resolution variation dn
  extern NamedFunc sys_llphoton_m_muresdn;

  //Higgs candidate mass with photon scale variation up
  extern NamedFunc sys_llphoton_m_phscaleup;

  //Higgs candidate mass with photon scale variation down
  extern NamedFunc sys_llphoton_m_phscaledn;

  //Higgs candidate mass with photon resolution variation up
  extern NamedFunc sys_llphoton_m_phresup;

  //Higgs candidate mass with photon resolution variation down
  extern NamedFunc sys_llphoton_m_phresdn;

  //Higgs candidate pt with electron scale variation up
  extern NamedFunc sys_llphoton_pt_elscaleup;

  //Higgs candidate pt with electron scale variation down
  extern NamedFunc sys_llphoton_pt_elscaledn;

  //Higgs candidate pt with electron resolution variation up
  extern NamedFunc sys_llphoton_pt_elresup;

  //Higgs candidate pt with electron resolution variation down
  extern NamedFunc sys_llphoton_pt_elresdn;

  //Higgs candidate pt with muon scale variation up
  extern NamedFunc sys_llphoton_pt_muscaleup;

  //Higgs candidate pt with muon scale variation dn
  extern NamedFunc sys_llphoton_pt_muscaledn;

  //Higgs candidate pt with muon resolution variation up
  extern NamedFunc sys_llphoton_pt_muresup;

  //Higgs candidate pt with muon resolution variation dn
  extern NamedFunc sys_llphoton_pt_muresdn;

  //Higgs candidate pt with photon scale variation up
  extern NamedFunc sys_llphoton_pt_phscaleup;

  //Higgs candidate pt with photon scale variation down
  extern NamedFunc sys_llphoton_pt_phscaledn;

  //Higgs candidate pt with photon resolution variation up
  extern NamedFunc sys_llphoton_pt_phresup;

  //Higgs candidate pt with photon resolution variation down
  extern NamedFunc sys_llphoton_pt_phresdn;

  //Gets NamedFunc that is (pt1, eta1, phi1, m1, pt2, eta2, phi2, m2) of lepton
  //refit pT with variation
  NamedFunc assign_variation_lep_refit(const NamedFunc &el_pt, 
      const NamedFunc &mu_pt, const NamedFunc &ll_lepid, 
      const NamedFunc &ll_i1, const NamedFunc &ll_i2, 
      std::shared_ptr<KinZfitter> kin_z_fitter, const std::string &name);

  //lepton 4 momentum with electron scale variation up
  extern NamedFunc sys_lep_refit_elscaleup;

  //lepton 4 momentum with electron scale variation down
  extern NamedFunc sys_lep_refit_elscaledn;

  //lepton 4 momentum with electron resolution variation up
  extern NamedFunc sys_lep_refit_elresup;

  //lepton 4 momentum with electron resolution variation down
  extern NamedFunc sys_lep_refit_elresdn;

  //lepton 4 momentum with muon scale variation up
  extern NamedFunc sys_lep_refit_muscaleup;

  //lepton 4 momentum with muon scale variation down
  extern NamedFunc sys_lep_refit_muscaledn;

  //lepton 4 momentum with muon resolution variation up
  extern NamedFunc sys_lep_refit_muresup;

  //lepton 4 momentum with muon resolution variation down
  extern NamedFunc sys_lep_refit_muresdn;

  //Gets NamedFunc that is assigns Z candidate four momentum with variation
  NamedFunc assign_variation_ll_refit_p4(const NamedFunc &ll_idx, 
      const NamedFunc &lep_refit, const std::string &name);

  //Z candidate refit 4 momentum without any variations
  extern NamedFunc sys_ll_refit_p4_default;

  //Z candidate refit 4 momentum with electron scale variation up
  extern NamedFunc sys_ll_refit_p4_elscaleup; 

  //Z candidate refit 4 momentum with electron scale variation down
  extern NamedFunc sys_ll_refit_p4_elscaledn; 

  //Z candidate refit 4 momentum with electron resolution variation up
  extern NamedFunc sys_ll_refit_p4_elresup; 

  //Z candidate refit 4 momentum with electron resolution variation down
  extern NamedFunc sys_ll_refit_p4_elresdn; 

  //Z candidate refit 4 momentum with muon scale variation up
  extern NamedFunc sys_ll_refit_p4_muscaleup; 

  //Z candidate refit 4 momentum with muon scale variation down
  extern NamedFunc sys_ll_refit_p4_muscaledn; 

  //Z candidate refit 4 momentum with muon resolution variation up
  extern NamedFunc sys_ll_refit_p4_muresup; 

  //Z candidate refit 4 momentum with muon resolution variation down
  extern NamedFunc sys_ll_refit_p4_muresdn; 

  //Gets NamedFunc that is assigns Higgs four momentum with variation
  NamedFunc assign_variation_llphoton_refit_p4(
      const NamedFunc &ll_idx, const NamedFunc &ll_refit_p4, 
      const NamedFunc &lead_photon_pt, const NamedFunc &lead_photon_eta, 
      const NamedFunc &lead_photon_phi, bool is_phvar, 
      const std::string &name);

  //H candidate refit 4 momentum with electron scale variation up
  extern NamedFunc sys_llphoton_refit_p4_elscaleup;
  
  //H candidate refit 4 momentum with electron scale variation down
  extern NamedFunc sys_llphoton_refit_p4_elscaledn;

  //H candidate refit 4 momentum with electron resolution variation up
  extern NamedFunc sys_llphoton_refit_p4_elresup;
  
  //H candidate refit 4 momentum with electron resolution variation down
  extern NamedFunc sys_llphoton_refit_p4_elresdn;

  //H candidate refit 4 momentum with muon scale variation up
  extern NamedFunc sys_llphoton_refit_p4_muscaleup;
  
  //H candidate refit 4 momentum with muon scale variation down
  extern NamedFunc sys_llphoton_refit_p4_muscaledn;

  //H candidate refit 4 momentum with muon resolution variation up
  extern NamedFunc sys_llphoton_refit_p4_muresup;
  
  //H candidate refit 4 momentum with muon resolution variation down
  extern NamedFunc sys_llphoton_refit_p4_muresdn;

  //H candidate refit 4 momentum with photon scale variation up
  extern NamedFunc sys_llphoton_refit_p4_phscaleup;

  //H candidate refit 4 momentum with photon scale variation dn
  extern NamedFunc sys_llphoton_refit_p4_phscaledn;

  //H candidate refit 4 momentum with photon resolution variation up
  extern NamedFunc sys_llphoton_refit_p4_phresup;

  //H candidate refit 4 momentum with photon resolution variation dn
  extern NamedFunc sys_llphoton_refit_p4_phresdn;

  //H candidate refit mass with electron scale variation up
  extern NamedFunc sys_llphoton_refit_m_elscaleup;

  //H candidate refit mass with electron scale variation down
  extern NamedFunc sys_llphoton_refit_m_elscaledn;

  //H candidate refit mass with electron resolution variation up
  extern NamedFunc sys_llphoton_refit_m_elresup;

  //H candidate refit mass with electron resolution variation down
  extern NamedFunc sys_llphoton_refit_m_elresdn;

  //H candidate refit mass with muon scale variation up
  extern NamedFunc sys_llphoton_refit_m_muscaleup;

  //H candidate refit mass with muon scale variation down
  extern NamedFunc sys_llphoton_refit_m_muscaledn;

  //H candidate refit mass with muon resolution variation up
  extern NamedFunc sys_llphoton_refit_m_muresup;

  //H candidate refit mass with muon resolution variation down
  extern NamedFunc sys_llphoton_refit_m_muresdn;

  //H candidate refit mass with photon scale variation up
  extern NamedFunc sys_llphoton_refit_m_phscaleup;

  //H candidate refit mass with photon scale variation down
  extern NamedFunc sys_llphoton_refit_m_phscaledn;

  //H candidate refit mass with photon resolution variation up
  extern NamedFunc sys_llphoton_refit_m_phresup;

  //H candidate refit mass with photon resolution variation down
  extern NamedFunc sys_llphoton_refit_m_phresdn;

  //helper function to get H candidate four vector after variation
  //doesn't consider edge cases ex. if selected Z candidate changes after
  //variation
  NamedFunc SimpleAssignVariationH(const NamedFunc &el_pt, 
                                   const NamedFunc &mu_pt,
                                   const NamedFunc &photon_pt,
                                   const std::string &name);

  //Gets NamedFunc that is jet_pt with year_specific_variation
  NamedFunc assign_variation_jet_pt_year(const NamedFunc &var_jet_pt,
      const std::string year, const std::string &name);

  //Gets vector NamedFunc that has variation for one year
  NamedFunc assign_vec_variation_year_select(const NamedFunc &var_namedfunc,
      const NamedFunc &def_namedfunc, const std::string year, 
      const std::string &name);
  
  //Gets scalar NamedFunc that has variation for one year
  NamedFunc assign_sca_variation_year_select(const NamedFunc &var_namedfunc,
      const NamedFunc &def_namedfunc, const std::string year, 
      const std::string &name);

  //Gets nbdfm with alternate sig
  NamedFunc assign_variation_nbdfm(const NamedFunc &var_jet_isgood,
      const std::string &name);

  //jet variations separated by era
  extern std::vector<NamedFunc> sys_jet_pt_scaleup;
  extern std::vector<NamedFunc> sys_jet_pt_scaledn;
  extern std::vector<NamedFunc> sys_jet_pt_resup;
  extern std::vector<NamedFunc> sys_jet_pt_resdn;
  extern std::vector<NamedFunc> sys_jet_isgood_scaleup;
  extern std::vector<NamedFunc> sys_jet_isgood_scaledn;
  extern std::vector<NamedFunc> sys_jet_isgood_resup;
  extern std::vector<NamedFunc> sys_jet_isgood_resdn;
  extern std::vector<NamedFunc> sys_lead_jet_pt_scaleup;
  extern std::vector<NamedFunc> sys_lead_jet_pt_scaledn;
  extern std::vector<NamedFunc> sys_lead_jet_pt_resup;
  extern std::vector<NamedFunc> sys_lead_jet_pt_resdn;
  extern std::vector<NamedFunc> sys_sublead_jet_pt_scaleup;
  extern std::vector<NamedFunc> sys_sublead_jet_pt_scaledn;
  extern std::vector<NamedFunc> sys_sublead_jet_pt_resup;
  extern std::vector<NamedFunc> sys_sublead_jet_pt_resdn;
  extern std::vector<NamedFunc> sys_lead_jet_eta_scaleup;
  extern std::vector<NamedFunc> sys_lead_jet_eta_scaledn;
  extern std::vector<NamedFunc> sys_lead_jet_eta_resup;
  extern std::vector<NamedFunc> sys_lead_jet_eta_resdn;
  extern std::vector<NamedFunc> sys_lead_jet_m_scaleup;
  extern std::vector<NamedFunc> sys_lead_jet_m_scaledn;
  extern std::vector<NamedFunc> sys_lead_jet_m_resup;
  extern std::vector<NamedFunc> sys_lead_jet_m_resdn;
  extern std::vector<NamedFunc> sys_njet_scaleup;
  extern std::vector<NamedFunc> sys_njet_scaledn;
  extern std::vector<NamedFunc> sys_njet_resup;
  extern std::vector<NamedFunc> sys_njet_resdn;
  extern std::vector<NamedFunc> sys_nbdfm_scaleup;
  extern std::vector<NamedFunc> sys_nbdfm_scaledn;
  extern std::vector<NamedFunc> sys_nbdfm_resup;
  extern std::vector<NamedFunc> sys_nbdfm_resdn;
  extern std::vector<NamedFunc> sys_met_scaleup;
  extern std::vector<NamedFunc> sys_met_scaledn;
  extern std::vector<NamedFunc> sys_met_resup;
  extern std::vector<NamedFunc> sys_met_resdn;
  
  //Gets NamedFunc that is untagged category selections with variation
  NamedFunc assign_variation_untagged_category(
      const NamedFunc &var_nlep, const NamedFunc &var_njet, 
      const NamedFunc &var_nbdfm, const NamedFunc &var_met,
      const NamedFunc &var_llphoton_pt, const NamedFunc &var_llphoton_m,
      const NamedFunc &var_max_lep_miniso, const NamedFunc &var_ll_m,
      const std::string &name);

  //untagged category selection with variations
  extern const NamedFunc untagged_category_elscaleup;
  extern const NamedFunc untagged_category_elscaledn;
  extern const NamedFunc untagged_category_elresup;
  extern const NamedFunc untagged_category_elresdn;
  extern const NamedFunc untagged_category_muscaleup;
  extern const NamedFunc untagged_category_muscaledn;
  extern const NamedFunc untagged_category_muresup;
  extern const NamedFunc untagged_category_muresdn;
  extern const NamedFunc untagged_category_phscaleup;
  extern const NamedFunc untagged_category_phscaledn;
  extern const NamedFunc untagged_category_phresup;
  extern const NamedFunc untagged_category_phresdn;
  extern std::vector<NamedFunc> untagged_category_jetscaleup;
  extern std::vector<NamedFunc> untagged_category_jetscaledn;
  extern std::vector<NamedFunc> untagged_category_jetresup;
  extern std::vector<NamedFunc> untagged_category_jetresdn;

  //horrible hack to get around the fact that I can't put a for loop outside
  //of a function to initialize what should be constant vectors...
  void initialize_jetvariations();
}

#endif
