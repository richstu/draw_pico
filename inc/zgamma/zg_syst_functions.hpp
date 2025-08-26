#ifndef H_SYST_FUNCTIONS
#define H_SYST_FUNCTIONS

#include <memory>
#include <vector>

#include "TLorentzVector.h"

#include "core/baby.hpp"
#include "core/fastforest.hpp"
#include "core/mva_wrapper.hpp"
#include "core/named_func.hpp"

#include "zgamma/KinZfitter.hpp"

namespace ZgFunctions {

  //utility functions
  double map_add(std::vector<double> input);
  double map_div(std::vector<double> input);
  double map_deltaphi(std::vector<double> phi_values);
  double map_deltar(std::vector<double> inputs);
  double reduce_index0(std::vector<double> inputs);
  double reduce_index1(std::vector<double> inputs);
  double reduce_index2(std::vector<double> inputs);
  double reduce_index3(std::vector<double> inputs);
  double reduce_index4(std::vector<double> inputs);
  double reduce_index5(std::vector<double> inputs);

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

  //weight implementing upward variation in electron trigger weights
  extern const NamedFunc sys_w_trig_el_up_pinnacles;

  //weight implementing downward variation in electron trigger weights
  extern const NamedFunc sys_w_trig_el_dn_pinnacles;

  //weight implementing upward variation in muon trigger weights
  extern const NamedFunc sys_w_trig_mu_up_pinnacles;

  //weight implementing downward variation in muon trigger weights
  extern const NamedFunc sys_w_trig_mu_dn_pinnacles;

  //el_sig and variations
  extern const NamedFunc sys_el_sig_scaleup;
  extern const NamedFunc sys_el_sig_scaleup;
  extern const NamedFunc sys_el_sig_scaledn;
  extern const NamedFunc sys_el_sig_resup;
  extern const NamedFunc sys_el_sig_resdn;

  //el_pt and variations
  extern const NamedFunc sys_el_pt_default;
  extern const NamedFunc sys_el_pt_scaleup;
  extern const NamedFunc sys_el_pt_scaledn;
  extern const NamedFunc sys_el_pt_renup;
  extern const NamedFunc sys_el_pt_rendn;

  //mu_sig and variations
  extern const NamedFunc sys_mu_sig_default;
  extern const NamedFunc sys_mu_sig_scaleup;
  extern const NamedFunc sys_mu_sig_scaledn;
  extern const NamedFunc sys_mu_sig_resup;
  extern const NamedFunc sys_mu_sig_resdn;

  //mu_pt and variations
  extern const NamedFunc sys_mu_pt_default;
  extern const NamedFunc sys_mu_pt_scaleup;
  extern const NamedFunc sys_mu_pt_scaledn;
  extern const NamedFunc sys_mu_pt_renup;
  extern const NamedFunc sys_mu_pt_rendn;

  //nel and variations
  extern const NamedFunc sys_nel_default;
  extern const NamedFunc sys_nel_scaleup;
  extern const NamedFunc sys_nel_scaledn;
  extern const NamedFunc sys_nel_resup;
  extern const NamedFunc sys_nel_resdn;

  //nmu and variations
  extern const NamedFunc sys_nmu_default;
  extern const NamedFunc sys_nmu_scaleup;
  extern const NamedFunc sys_nmu_scaledn;
  extern const NamedFunc sys_nmu_resup;
  extern const NamedFunc sys_nmu_resdn;

  //nlep and variations
  extern const NamedFunc sys_nlep_default;
  extern const NamedFunc sys_nlep_elscaleup;
  extern const NamedFunc sys_nlep_elscaledn;
  extern const NamedFunc sys_nlep_elresup;
  extern const NamedFunc sys_nlep_elresdn;
  extern const NamedFunc sys_nlep_muscaleup;
  extern const NamedFunc sys_nlep_muscaledn;
  extern const NamedFunc sys_nlep_muresup;
  extern const NamedFunc sys_nlep_muresdn;

  //Gets NamedFunc that is number of ll candidates with variation
  NamedFunc assign_variation_nll(const NamedFunc &variation_el_sig, 
                                 const NamedFunc &variation_mu_sig,
                                 const std::string &name);

  //nll and variations
  extern const NamedFunc sys_nll_default;
  extern const NamedFunc sys_nll_elscaleup;
  extern const NamedFunc sys_nll_elscaledn;
  extern const NamedFunc sys_nll_elresup;
  extern const NamedFunc sys_nll_elresdn;
  extern const NamedFunc sys_nll_muscaleup;
  extern const NamedFunc sys_nll_muscaledn;
  extern const NamedFunc sys_nll_muresup;
  extern const NamedFunc sys_nll_muresdn;

  //leading electron pt and variations
  extern const NamedFunc sys_lead_el_pt_default;
  extern const NamedFunc sys_lead_el_pt_scaleup;
  extern const NamedFunc sys_lead_el_pt_scaledn;
  extern const NamedFunc sys_lead_el_pt_resup;
  extern const NamedFunc sys_lead_el_pt_resdn;

  //subleading electron pt and variations
  extern const NamedFunc sys_sublead_el_pt_default;
  extern const NamedFunc sys_sublead_el_pt_scaleup;
  extern const NamedFunc sys_sublead_el_pt_scaledn;
  extern const NamedFunc sys_sublead_el_pt_resup;
  extern const NamedFunc sys_sublead_el_pt_resdn;

  //leading muon pt and variations
  extern const NamedFunc sys_lead_mu_pt_default;
  extern const NamedFunc sys_lead_mu_pt_scaleup;
  extern const NamedFunc sys_lead_mu_pt_scaledn;
  extern const NamedFunc sys_lead_mu_pt_resup;
  extern const NamedFunc sys_lead_mu_pt_resdn;

  //subleading muon and variations
  extern const NamedFunc sys_sublead_mu_pt_default;
  extern const NamedFunc sys_sublead_mu_pt_scaleup;
  extern const NamedFunc sys_sublead_mu_pt_scaledn;
  extern const NamedFunc sys_sublead_mu_pt_resup;
  extern const NamedFunc sys_sublead_mu_pt_resdn;

  //Gets NamedFunc that is trigger+pT cut flag with variation
  NamedFunc assign_variation_trig_pt(const NamedFunc &variation_el_pt, 
                                     const NamedFunc &variation_el_sig,
                                     const NamedFunc &variation_mu_pt,
                                     const NamedFunc &variation_mu_sig,
                                     const std::string &name);

  //trigger and pT cuts and variations
  extern const NamedFunc sys_trig_pt_default;
  extern const NamedFunc sys_trig_pt_elscaleup;
  extern const NamedFunc sys_trig_pt_elscaledn;
  extern const NamedFunc sys_trig_pt_elresup;
  extern const NamedFunc sys_trig_pt_elresdn;
  extern const NamedFunc sys_trig_pt_muscaleup;
  extern const NamedFunc sys_trig_pt_muscaledn;
  extern const NamedFunc sys_trig_pt_muresup;
  extern const NamedFunc sys_trig_pt_muresdn;

  //Gets NamedFunc that is max lep miniso with variation
  NamedFunc assign_variation_max_lep_miniso(const NamedFunc &variation_el_sig,
                                            const NamedFunc &variation_mu_sig,
                                            const std::string &name);

  //max lep miniso and variations
  extern const NamedFunc sys_max_lep_miniso_default;
  extern const NamedFunc sys_max_lep_miniso_elscaleup;
  extern const NamedFunc sys_max_lep_miniso_elscaledn;
  extern const NamedFunc sys_max_lep_miniso_elresup;
  extern const NamedFunc sys_max_lep_miniso_elresdn;
  extern const NamedFunc sys_max_lep_miniso_muscaleup;
  extern const NamedFunc sys_max_lep_miniso_muscaledn;
  extern const NamedFunc sys_max_lep_miniso_muresup;
  extern const NamedFunc sys_max_lep_miniso_muresdn;

  //if photon is also a signal electron
  extern const NamedFunc sys_photon_sigel;

  //photon_sig and variations
  extern const NamedFunc sys_photon_sig_default;
  extern const NamedFunc sys_photon_sig_scaleup;
  extern const NamedFunc sys_photon_sig_scaledn;
  extern const NamedFunc sys_photon_sig_resup;
  extern const NamedFunc sys_photon_sig_resdn;

  //photon_pt and variations
  extern const NamedFunc sys_photon_pt_default;
  extern const NamedFunc sys_photon_pt_scaleup;
  extern const NamedFunc sys_photon_pt_scaledn;
  extern const NamedFunc sys_photon_pt_resup;
  extern const NamedFunc sys_photon_pt_resdn;

  //some defaults
  extern const NamedFunc sys_photon_eta_default;
  extern const NamedFunc sys_photon_phi_default;
  extern const NamedFunc sys_photon_drmin_default;
  extern const NamedFunc sys_photon_drmax_default;
  extern const NamedFunc sys_photon_idmva_default;

  //nphoton and variations
  extern const NamedFunc sys_nphoton_default;
  extern const NamedFunc sys_nphoton_scaleup;
  extern const NamedFunc sys_nphoton_scaledn;
  extern const NamedFunc sys_nphoton_resup;
  extern const NamedFunc sys_nphoton_resdn;

  //signal photon properties and variations
  extern const NamedFunc sys_sig_photon_pt_default;
  extern const NamedFunc sys_sig_photon_pt_scaleup;
  extern const NamedFunc sys_sig_photon_pt_scaledn;
  extern const NamedFunc sys_sig_photon_pt_resup;
  extern const NamedFunc sys_sig_photon_pt_resdn;
  extern const NamedFunc sys_sig_photon_eta_default;
  extern const NamedFunc sys_sig_photon_eta_scaleup;
  extern const NamedFunc sys_sig_photon_eta_scaledn;
  extern const NamedFunc sys_sig_photon_eta_resup;
  extern const NamedFunc sys_sig_photon_eta_resdn;
  extern const NamedFunc sys_sig_photon_phi_default;
  extern const NamedFunc sys_sig_photon_phi_scaleup;
  extern const NamedFunc sys_sig_photon_phi_scaledn;
  extern const NamedFunc sys_sig_photon_phi_resup;
  extern const NamedFunc sys_sig_photon_phi_resdn;
  extern const NamedFunc sys_sig_photon_drmin_default;
  extern const NamedFunc sys_sig_photon_drmin_scaleup;
  extern const NamedFunc sys_sig_photon_drmin_scaledn;
  extern const NamedFunc sys_sig_photon_drmin_resup;
  extern const NamedFunc sys_sig_photon_drmin_resdn;
  extern const NamedFunc sys_sig_photon_drmax_default;
  extern const NamedFunc sys_sig_photon_drmax_scaleup;
  extern const NamedFunc sys_sig_photon_drmax_scaledn;
  extern const NamedFunc sys_sig_photon_drmax_resup;
  extern const NamedFunc sys_sig_photon_drmax_resdn;
  extern const NamedFunc sys_sig_photon_idmva_default;
  extern const NamedFunc sys_sig_photon_idmva_scaleup;
  extern const NamedFunc sys_sig_photon_idmva_scaledn;
  extern const NamedFunc sys_sig_photon_idmva_resup;
  extern const NamedFunc sys_sig_photon_idmva_resdn;

  //leading photon properties and variations
  extern const NamedFunc sys_lead_photon_pt_default;
  extern const NamedFunc sys_lead_photon_pt_scaleup;
  extern const NamedFunc sys_lead_photon_pt_scaledn;
  extern const NamedFunc sys_lead_photon_pt_resup;
  extern const NamedFunc sys_lead_photon_pt_resdn;
  extern const NamedFunc sys_lead_photon_eta_default;
  extern const NamedFunc sys_lead_photon_eta_scaleup;
  extern const NamedFunc sys_lead_photon_eta_scaledn;
  extern const NamedFunc sys_lead_photon_eta_resup;
  extern const NamedFunc sys_lead_photon_eta_resdn;
  extern const NamedFunc sys_lead_photon_phi_default;
  extern const NamedFunc sys_lead_photon_phi_scaleup;
  extern const NamedFunc sys_lead_photon_phi_scaledn;
  extern const NamedFunc sys_lead_photon_phi_resup;
  extern const NamedFunc sys_lead_photon_phi_resdn;
  extern const NamedFunc sys_lead_photon_drmin_default;
  extern const NamedFunc sys_lead_photon_drmin_scaleup;
  extern const NamedFunc sys_lead_photon_drmin_scaledn;
  extern const NamedFunc sys_lead_photon_drmin_resup;
  extern const NamedFunc sys_lead_photon_drmin_resdn;
  extern const NamedFunc sys_lead_photon_drmax_default;
  extern const NamedFunc sys_lead_photon_drmax_scaleup;
  extern const NamedFunc sys_lead_photon_drmax_scaledn;
  extern const NamedFunc sys_lead_photon_drmax_resup;
  extern const NamedFunc sys_lead_photon_drmax_resdn;
  extern const NamedFunc sys_lead_photon_idmva_default;
  extern const NamedFunc sys_lead_photon_idmva_scaleup;
  extern const NamedFunc sys_lead_photon_idmva_scaledn;
  extern const NamedFunc sys_lead_photon_idmva_resup;
  extern const NamedFunc sys_lead_photon_idmva_resdn;

  //get photon rel pt error with variation
  NamedFunc assign_variation_lead_photon_relpterr(const NamedFunc &ph_sig,
      const NamedFunc &ph_pt, const std::string &name);

  //photon relpterr variations
  extern const NamedFunc sys_lead_photon_relpterr_default;
  extern const NamedFunc sys_lead_photon_relpterr_scaleup;
  extern const NamedFunc sys_lead_photon_relpterr_scaledn;
  extern const NamedFunc sys_lead_photon_relpterr_resup;
  extern const NamedFunc sys_lead_photon_relpterr_resdn;

  //get Z candidate properties with variation
  //returns (pt, eta, phi, m, lepid, i1, i2, idx)
  NamedFunc assign_variation_ll(const NamedFunc &el_pt, 
                                const NamedFunc &el_sig,
                                const NamedFunc &mu_pt, 
                                const NamedFunc &mu_sig,
                                bool is_el,
                                const std::string &name);

  //dilepton properties and variations
  extern const NamedFunc sys_ll_default;
  extern const NamedFunc sys_ll_elscaleup;
  extern const NamedFunc sys_ll_elscaledn;
  extern const NamedFunc sys_ll_elresup;
  extern const NamedFunc sys_ll_elresdn;
  extern const NamedFunc sys_ll_muscaleup;
  extern const NamedFunc sys_ll_muscaledn;
  extern const NamedFunc sys_ll_muresup;
  extern const NamedFunc sys_ll_muresdn;
  extern const NamedFunc sys_ll_m_default;
  extern const NamedFunc sys_ll_m_elscaleup;
  extern const NamedFunc sys_ll_m_elscaledn;
  extern const NamedFunc sys_ll_m_elresup;
  extern const NamedFunc sys_ll_m_elresdn;
  extern const NamedFunc sys_ll_m_muscaleup;
  extern const NamedFunc sys_ll_m_muscaledn;
  extern const NamedFunc sys_ll_m_muresup;
  extern const NamedFunc sys_ll_m_muresdn;

  //get lepton eta with variation
  NamedFunc assign_variation_lep_eta(const NamedFunc &ll, const bool lead,
      const std::string &name);

  //lepton eta with variations
  extern const NamedFunc sys_lead_lepton_eta_default;
  extern const NamedFunc sys_lead_lepton_eta_elscaleup;
  extern const NamedFunc sys_lead_lepton_eta_elscaledn;
  extern const NamedFunc sys_lead_lepton_eta_elresup;
  extern const NamedFunc sys_lead_lepton_eta_elresdn;
  extern const NamedFunc sys_lead_lepton_eta_muscaleup;
  extern const NamedFunc sys_lead_lepton_eta_muscaledn;
  extern const NamedFunc sys_lead_lepton_eta_muresup;
  extern const NamedFunc sys_lead_lepton_eta_muresdn;
  extern const NamedFunc sys_sublead_lepton_eta_default;
  extern const NamedFunc sys_sublead_lepton_eta_elscaleup;
  extern const NamedFunc sys_sublead_lepton_eta_elscaledn;
  extern const NamedFunc sys_sublead_lepton_eta_elresup;
  extern const NamedFunc sys_sublead_lepton_eta_elresdn;
  extern const NamedFunc sys_sublead_lepton_eta_muscaleup;
  extern const NamedFunc sys_sublead_lepton_eta_muscaledn;
  extern const NamedFunc sys_sublead_lepton_eta_muresup;
  extern const NamedFunc sys_sublead_lepton_eta_muresdn;

  //get H candidate properties with variation
  //returns (pt, eta, phi, m)
  NamedFunc assign_variation_llphoton_p4(const NamedFunc &ll_p4, 
                                         const NamedFunc &lead_photon_pt,
                                         const NamedFunc &lead_photon_eta,
                                         const NamedFunc &lead_photon_phi,
                                         const int var_pdgid,
                                         const std::string &name);

  //H four momentum and variations
  extern const NamedFunc sys_llphoton_p4_elscaleup;
  extern const NamedFunc sys_llphoton_p4_elscaledn;
  extern const NamedFunc sys_llphoton_p4_elresup;
  extern const NamedFunc sys_llphoton_p4_elresdn;
  extern const NamedFunc sys_llphoton_p4_muscaleup;
  extern const NamedFunc sys_llphoton_p4_muscaledn;
  extern const NamedFunc sys_llphoton_p4_muresup;
  extern const NamedFunc sys_llphoton_p4_muresdn;
  extern const NamedFunc sys_llphoton_p4_phscaleup;
  extern const NamedFunc sys_llphoton_p4_phscaledn;
  extern const NamedFunc sys_llphoton_p4_phresup;
  extern const NamedFunc sys_llphoton_p4_phresdn;

  //Higgs candidate mass and variations
  extern const NamedFunc sys_llphoton_m_default;
  extern const NamedFunc sys_llphoton_m_elscaleup;
  extern const NamedFunc sys_llphoton_m_elscaledn;
  extern const NamedFunc sys_llphoton_m_elresup;
  extern const NamedFunc sys_llphoton_m_elresdn;
  extern const NamedFunc sys_llphoton_m_muscaleup;
  extern const NamedFunc sys_llphoton_m_muscaledn;
  extern const NamedFunc sys_llphoton_m_muresup;
  extern const NamedFunc sys_llphoton_m_muresdn;
  extern const NamedFunc sys_llphoton_m_phscaleup;
  extern const NamedFunc sys_llphoton_m_phscaledn;
  extern const NamedFunc sys_llphoton_m_phresup;
  extern const NamedFunc sys_llphoton_m_phresdn;

  //Higgs candidate pt and variations
  extern const NamedFunc sys_llphoton_pt_default;
  extern const NamedFunc sys_llphoton_pt_elscaleup;
  extern const NamedFunc sys_llphoton_pt_elscaledn;
  extern const NamedFunc sys_llphoton_pt_elresup;
  extern const NamedFunc sys_llphoton_pt_elresdn;
  extern const NamedFunc sys_llphoton_pt_muscaleup;
  extern const NamedFunc sys_llphoton_pt_muscaledn;
  extern const NamedFunc sys_llphoton_pt_muresup;
  extern const NamedFunc sys_llphoton_pt_muresdn;
  extern const NamedFunc sys_llphoton_pt_phscaleup;
  extern const NamedFunc sys_llphoton_pt_phscaledn;
  extern const NamedFunc sys_llphoton_pt_phresup;
  extern const NamedFunc sys_llphoton_pt_phresdn;

  //Gets NamedFunc that is (pt1, eta1, phi1, m1, pt2, eta2, phi2, m2) of lepton
  //refit pT with variation
  NamedFunc assign_variation_lep_refit(const NamedFunc &el_pt, 
      const NamedFunc &mu_pt, const NamedFunc &ll, const std::string &name);

  //lepton 4 momentum and variations
  extern const NamedFunc sys_lep_refit_default;
  extern const NamedFunc sys_lep_refit_elscaleup;
  extern const NamedFunc sys_lep_refit_elscaledn;
  extern const NamedFunc sys_lep_refit_elresup;
  extern const NamedFunc sys_lep_refit_elresdn;
  extern const NamedFunc sys_lep_refit_muscaleup;
  extern const NamedFunc sys_lep_refit_muscaledn;
  extern const NamedFunc sys_lep_refit_muresup;
  extern const NamedFunc sys_lep_refit_muresdn;

  //Gets NamedFunc that is assigns Z candidate four momentum with variation
  NamedFunc assign_variation_ll_refit_p4(const NamedFunc &ll_idx, 
      const NamedFunc &lep_refit, const std::string &name);

  //Z candidate refit 4 momentum and variations
  extern const NamedFunc sys_ll_refit_p4_default;
  extern const NamedFunc sys_ll_refit_p4_elscaleup; 
  extern const NamedFunc sys_ll_refit_p4_elscaledn; 
  extern const NamedFunc sys_ll_refit_p4_elresup; 
  extern const NamedFunc sys_ll_refit_p4_elresdn; 
  extern const NamedFunc sys_ll_refit_p4_muscaleup; 
  extern const NamedFunc sys_ll_refit_p4_muscaledn; 
  extern const NamedFunc sys_ll_refit_p4_muresup; 
  extern const NamedFunc sys_ll_refit_p4_muresdn; 

  //Gets NamedFunc that is assigns Higgs four momentum with variation
  NamedFunc assign_variation_llphoton_refit_p4(
      const NamedFunc &ll_idx, const NamedFunc &ll_refit_p4, 
      const NamedFunc &lead_photon_pt, const NamedFunc &lead_photon_eta, 
      const NamedFunc &lead_photon_phi, bool is_phvar, 
      const std::string &name);

  //H candidate refit 4 momentum and variations
  extern const NamedFunc sys_llphoton_refit_p4_elscaleup;
  extern const NamedFunc sys_llphoton_refit_p4_elscaledn;
  extern const NamedFunc sys_llphoton_refit_p4_elresup;
  extern const NamedFunc sys_llphoton_refit_p4_elresdn;
  extern const NamedFunc sys_llphoton_refit_p4_muscaleup;
  extern const NamedFunc sys_llphoton_refit_p4_muscaledn;
  extern const NamedFunc sys_llphoton_refit_p4_muresup;
  extern const NamedFunc sys_llphoton_refit_p4_muresdn;
  extern const NamedFunc sys_llphoton_refit_p4_phscaleup;
  extern const NamedFunc sys_llphoton_refit_p4_phscaledn;
  extern const NamedFunc sys_llphoton_refit_p4_phresup;
  extern const NamedFunc sys_llphoton_refit_p4_phresdn;

  //H candidate refit pt and variations
  extern const NamedFunc sys_llphoton_refit_pt_default;
  extern const NamedFunc sys_llphoton_refit_pt_elscaleup;
  extern const NamedFunc sys_llphoton_refit_pt_elscaledn;
  extern const NamedFunc sys_llphoton_refit_pt_elresup;
  extern const NamedFunc sys_llphoton_refit_pt_elresdn;
  extern const NamedFunc sys_llphoton_refit_pt_muscaleup;
  extern const NamedFunc sys_llphoton_refit_pt_muscaledn;
  extern const NamedFunc sys_llphoton_refit_pt_muresup;
  extern const NamedFunc sys_llphoton_refit_pt_muresdn;
  extern const NamedFunc sys_llphoton_refit_pt_phscaleup;
  extern const NamedFunc sys_llphoton_refit_pt_phscaledn;
  extern const NamedFunc sys_llphoton_refit_pt_phresup;
  extern const NamedFunc sys_llphoton_refit_pt_phresdn;

  //H candidate refit m and variations
  extern const NamedFunc sys_llphoton_refit_m_default;
  extern const NamedFunc sys_llphoton_refit_m_elscaleup;
  extern const NamedFunc sys_llphoton_refit_m_elscaledn;
  extern const NamedFunc sys_llphoton_refit_m_elresup;
  extern const NamedFunc sys_llphoton_refit_m_elresdn;
  extern const NamedFunc sys_llphoton_refit_m_muscaleup;
  extern const NamedFunc sys_llphoton_refit_m_muscaledn;
  extern const NamedFunc sys_llphoton_refit_m_muresup;
  extern const NamedFunc sys_llphoton_refit_m_muresdn;
  extern const NamedFunc sys_llphoton_refit_m_phscaleup;
  extern const NamedFunc sys_llphoton_refit_m_phscaledn;
  extern const NamedFunc sys_llphoton_refit_m_phresup;
  extern const NamedFunc sys_llphoton_refit_m_phresdn;

  //H candidate refit pt/m and variations
  extern const NamedFunc sys_llphoton_refit_relpt_default;
  extern const NamedFunc sys_llphoton_refit_relpt_elscaleup;
  extern const NamedFunc sys_llphoton_refit_relpt_elscaledn;
  extern const NamedFunc sys_llphoton_refit_relpt_elresup;
  extern const NamedFunc sys_llphoton_refit_relpt_elresdn;
  extern const NamedFunc sys_llphoton_refit_relpt_muscaleup;
  extern const NamedFunc sys_llphoton_refit_relpt_muscaledn;
  extern const NamedFunc sys_llphoton_refit_relpt_muresup;
  extern const NamedFunc sys_llphoton_refit_relpt_muresdn;
  extern const NamedFunc sys_llphoton_refit_relpt_phscaleup;
  extern const NamedFunc sys_llphoton_refit_relpt_phscaledn;
  extern const NamedFunc sys_llphoton_refit_relpt_phresup;
  extern const NamedFunc sys_llphoton_refit_relpt_phresdn;

  //H candidate refit phi and variations
  extern const NamedFunc sys_llphoton_refit_phi_elscaleup;
  extern const NamedFunc sys_llphoton_refit_phi_elscaledn;
  extern const NamedFunc sys_llphoton_refit_phi_elresup;
  extern const NamedFunc sys_llphoton_refit_phi_elresdn;
  extern const NamedFunc sys_llphoton_refit_phi_muscaleup;
  extern const NamedFunc sys_llphoton_refit_phi_muscaledn;
  extern const NamedFunc sys_llphoton_refit_phi_muresup;
  extern const NamedFunc sys_llphoton_refit_phi_muresdn;
  extern const NamedFunc sys_llphoton_refit_phi_phscaleup;
  extern const NamedFunc sys_llphoton_refit_phi_phscaledn;
  extern const NamedFunc sys_llphoton_refit_phi_phresup;
  extern const NamedFunc sys_llphoton_refit_phi_phresdn;

  //Gets NamedFunc that assigns kinematic angles with variation
  NamedFunc assign_variation_llphoton_refit_angles(const NamedFunc &var_ll,
      const NamedFunc &lep_refit, const NamedFunc &var_photon_pt, 
      const NamedFunc &var_photon_eta, const NamedFunc &var_photon_phi, 
      const std::string name);

  //kinematic angles with variations
  extern const NamedFunc sys_llphoton_refit_angles_default;
  extern const NamedFunc sys_llphoton_refit_angles_elscaleup;
  extern const NamedFunc sys_llphoton_refit_angles_elscaledn;
  extern const NamedFunc sys_llphoton_refit_angles_elresup;
  extern const NamedFunc sys_llphoton_refit_angles_elresdn;
  extern const NamedFunc sys_llphoton_refit_angles_muscaleup;
  extern const NamedFunc sys_llphoton_refit_angles_muscaledn;
  extern const NamedFunc sys_llphoton_refit_angles_muresup;
  extern const NamedFunc sys_llphoton_refit_angles_muresdn;
  extern const NamedFunc sys_llphoton_refit_angles_phscaleup;
  extern const NamedFunc sys_llphoton_refit_angles_phscaledn;
  extern const NamedFunc sys_llphoton_refit_angles_phresup;
  extern const NamedFunc sys_llphoton_refit_angles_phresdn;
  extern const NamedFunc sys_llphoton_refit_cosTheta_default;
  extern const NamedFunc sys_llphoton_refit_cosTheta_elscaleup;
  extern const NamedFunc sys_llphoton_refit_cosTheta_elscaledn;
  extern const NamedFunc sys_llphoton_refit_cosTheta_elresup;
  extern const NamedFunc sys_llphoton_refit_cosTheta_elresdn;
  extern const NamedFunc sys_llphoton_refit_cosTheta_muscaleup;
  extern const NamedFunc sys_llphoton_refit_cosTheta_muscaledn;
  extern const NamedFunc sys_llphoton_refit_cosTheta_muresup;
  extern const NamedFunc sys_llphoton_refit_cosTheta_muresdn;
  extern const NamedFunc sys_llphoton_refit_cosTheta_phscaleup;
  extern const NamedFunc sys_llphoton_refit_cosTheta_phscaledn;
  extern const NamedFunc sys_llphoton_refit_cosTheta_phresup;
  extern const NamedFunc sys_llphoton_refit_cosTheta_phresdn;
  extern const NamedFunc sys_llphoton_refit_costheta_default;
  extern const NamedFunc sys_llphoton_refit_costheta_elscaleup;
  extern const NamedFunc sys_llphoton_refit_costheta_elscaledn;
  extern const NamedFunc sys_llphoton_refit_costheta_elresup;
  extern const NamedFunc sys_llphoton_refit_costheta_elresdn;
  extern const NamedFunc sys_llphoton_refit_costheta_muscaleup;
  extern const NamedFunc sys_llphoton_refit_costheta_muscaledn;
  extern const NamedFunc sys_llphoton_refit_costheta_muresup;
  extern const NamedFunc sys_llphoton_refit_costheta_muresdn;
  extern const NamedFunc sys_llphoton_refit_costheta_phscaleup;
  extern const NamedFunc sys_llphoton_refit_costheta_phscaledn;
  extern const NamedFunc sys_llphoton_refit_costheta_phresup;
  extern const NamedFunc sys_llphoton_refit_costheta_phresdn;
  extern const NamedFunc sys_llphoton_refit_psi_default;
  extern const NamedFunc sys_llphoton_refit_psi_elscaleup;
  extern const NamedFunc sys_llphoton_refit_psi_elscaledn;
  extern const NamedFunc sys_llphoton_refit_psi_elresup;
  extern const NamedFunc sys_llphoton_refit_psi_elresdn;
  extern const NamedFunc sys_llphoton_refit_psi_muscaleup;
  extern const NamedFunc sys_llphoton_refit_psi_muscaledn;
  extern const NamedFunc sys_llphoton_refit_psi_muresup;
  extern const NamedFunc sys_llphoton_refit_psi_muresdn;
  extern const NamedFunc sys_llphoton_refit_psi_phscaleup;
  extern const NamedFunc sys_llphoton_refit_psi_phscaledn;
  extern const NamedFunc sys_llphoton_refit_psi_phresup;
  extern const NamedFunc sys_llphoton_refit_psi_phresdn;

  //helper function to get H candidate four vector after variation
  //doesn't consider edge cases ex. if selected Z candidate changes after
  //variation
  NamedFunc SimpleAssignVariationH(const NamedFunc &el_pt, 
                                   const NamedFunc &mu_pt,
                                   const NamedFunc &photon_pt,
                                   const std::string &name);

  //Gets NamedFunc that assigns dijet properties (pt eta phi m deta dphi)
  NamedFunc assign_variation_dijet(const NamedFunc &jet_sig, 
      const NamedFunc &jet_pt, const NamedFunc &jet_m, 
      const std::string &name);

  //Gets vector NamedFunc that has variation for one year
  NamedFunc assign_vec_variation_year_select(const NamedFunc &var_namedfunc,
      const NamedFunc &def_namedfunc, const std::string &year, 
      const std::string &name);
  
  //Gets scalar NamedFunc that has variation for one year
  NamedFunc assign_sca_variation_year_select(const NamedFunc &var_namedfunc,
      const NamedFunc &def_namedfunc, const std::string &year, 
      const std::string &name);

  //Gets nbdfm with alternate sig
  NamedFunc assign_variation_nbdfm(const NamedFunc &var_jet_isgood,
      const std::string &name);

  //jet variations separated by era
  extern const NamedFunc sys_jet_pt_default;
  extern std::vector<NamedFunc> sys_jet_pt_scaleup;
  extern std::vector<NamedFunc> sys_jet_pt_scaledn;
  extern std::vector<NamedFunc> sys_jet_pt_resup;
  extern std::vector<NamedFunc> sys_jet_pt_resdn;
  extern const NamedFunc sys_jet_m_default;
  extern std::vector<NamedFunc> sys_jet_m_scaleup;
  extern std::vector<NamedFunc> sys_jet_m_scaledn;
  extern std::vector<NamedFunc> sys_jet_m_resup;
  extern std::vector<NamedFunc> sys_jet_m_resdn;
  extern const NamedFunc sys_jet_isgood_default;
  extern std::vector<NamedFunc> sys_jet_isgood_scaleup;
  extern std::vector<NamedFunc> sys_jet_isgood_scaledn;
  extern std::vector<NamedFunc> sys_jet_isgood_resup;
  extern std::vector<NamedFunc> sys_jet_isgood_resdn;
  extern const NamedFunc sys_sig_jet_pt_default;
  extern std::vector<NamedFunc> sys_sig_jet_pt_scaleup;
  extern std::vector<NamedFunc> sys_sig_jet_pt_scaledn;
  extern std::vector<NamedFunc> sys_sig_jet_pt_resup;
  extern std::vector<NamedFunc> sys_sig_jet_pt_resdn;
  extern const NamedFunc sys_sig_jet_eta_default;
  extern std::vector<NamedFunc> sys_sig_jet_eta_scaleup;
  extern std::vector<NamedFunc> sys_sig_jet_eta_scaledn;
  extern std::vector<NamedFunc> sys_sig_jet_eta_resup;
  extern std::vector<NamedFunc> sys_sig_jet_eta_resdn;
  extern const NamedFunc sys_sig_jet_phi_default;
  extern std::vector<NamedFunc> sys_sig_jet_phi_scaleup;
  extern std::vector<NamedFunc> sys_sig_jet_phi_scaledn;
  extern std::vector<NamedFunc> sys_sig_jet_phi_resup;
  extern std::vector<NamedFunc> sys_sig_jet_phi_resdn;
  extern const NamedFunc sys_sig_jet_m_default;
  extern std::vector<NamedFunc> sys_sig_jet_m_scaleup;
  extern std::vector<NamedFunc> sys_sig_jet_m_scaledn;
  extern std::vector<NamedFunc> sys_sig_jet_m_resup;
  extern std::vector<NamedFunc> sys_sig_jet_m_resdn;
  extern const NamedFunc sys_lead_jet_pt_default;
  extern std::vector<NamedFunc> sys_lead_jet_pt_scaleup;
  extern std::vector<NamedFunc> sys_lead_jet_pt_scaledn;
  extern std::vector<NamedFunc> sys_lead_jet_pt_resup;
  extern std::vector<NamedFunc> sys_lead_jet_pt_resdn;
  extern const NamedFunc sys_lead_jet_eta_default;
  extern std::vector<NamedFunc> sys_lead_jet_eta_scaleup;
  extern std::vector<NamedFunc> sys_lead_jet_eta_scaledn;
  extern std::vector<NamedFunc> sys_lead_jet_eta_resup;
  extern std::vector<NamedFunc> sys_lead_jet_eta_resdn;
  extern const NamedFunc sys_lead_jet_phi_default;
  extern std::vector<NamedFunc> sys_lead_jet_phi_scaleup;
  extern std::vector<NamedFunc> sys_lead_jet_phi_scaledn;
  extern std::vector<NamedFunc> sys_lead_jet_phi_resup;
  extern std::vector<NamedFunc> sys_lead_jet_phi_resdn;
  extern const NamedFunc sys_lead_jet_m_default;
  extern std::vector<NamedFunc> sys_lead_jet_m_scaleup;
  extern std::vector<NamedFunc> sys_lead_jet_m_scaledn;
  extern std::vector<NamedFunc> sys_lead_jet_m_resup;
  extern std::vector<NamedFunc> sys_lead_jet_m_resdn;
  extern const NamedFunc sys_sublead_jet_pt_default;
  extern std::vector<NamedFunc> sys_sublead_jet_pt_scaleup;
  extern std::vector<NamedFunc> sys_sublead_jet_pt_scaledn;
  extern std::vector<NamedFunc> sys_sublead_jet_pt_resup;
  extern std::vector<NamedFunc> sys_sublead_jet_pt_resdn;
  extern const NamedFunc sys_sublead_jet_eta_default;
  extern std::vector<NamedFunc> sys_sublead_jet_eta_scaleup;
  extern std::vector<NamedFunc> sys_sublead_jet_eta_scaledn;
  extern std::vector<NamedFunc> sys_sublead_jet_eta_resup;
  extern std::vector<NamedFunc> sys_sublead_jet_eta_resdn;
  extern const NamedFunc sys_sublead_jet_phi_default;
  extern std::vector<NamedFunc> sys_sublead_jet_phi_scaleup;
  extern std::vector<NamedFunc> sys_sublead_jet_phi_scaledn;
  extern std::vector<NamedFunc> sys_sublead_jet_phi_resup;
  extern std::vector<NamedFunc> sys_sublead_jet_phi_resdn;
  extern const NamedFunc sys_njet_default;
  extern std::vector<NamedFunc> sys_njet_scaleup;
  extern std::vector<NamedFunc> sys_njet_scaledn;
  extern std::vector<NamedFunc> sys_njet_resup;
  extern std::vector<NamedFunc> sys_njet_resdn;
  extern const NamedFunc sys_nbdfm_default;
  extern std::vector<NamedFunc> sys_nbdfm_scaleup;
  extern std::vector<NamedFunc> sys_nbdfm_scaledn;
  extern std::vector<NamedFunc> sys_nbdfm_resup;
  extern std::vector<NamedFunc> sys_nbdfm_resdn;
  extern const NamedFunc sys_met_default;
  extern std::vector<NamedFunc> sys_met_scaleup;
  extern std::vector<NamedFunc> sys_met_scaledn;
  extern std::vector<NamedFunc> sys_met_resup;
  extern std::vector<NamedFunc> sys_met_resdn;

  //dijet variations
  extern const NamedFunc sys_dijet_default;
  extern const NamedFunc sys_dijet_phi_default;
  extern const NamedFunc sys_dijet_m_default;
  extern const NamedFunc sys_dijet_deta_default;
  extern const NamedFunc sys_dijet_dphi_default;
  extern std::vector<NamedFunc> sys_dijet_scaleup;
  extern std::vector<NamedFunc> sys_dijet_scaledn;
  extern std::vector<NamedFunc> sys_dijet_resup;
  extern std::vector<NamedFunc> sys_dijet_resdn;
  extern std::vector<NamedFunc> sys_dijet_phi_scaleup;
  extern std::vector<NamedFunc> sys_dijet_phi_scaledn;
  extern std::vector<NamedFunc> sys_dijet_phi_resup;
  extern std::vector<NamedFunc> sys_dijet_phi_resdn;
  extern std::vector<NamedFunc> sys_dijet_m_scaleup;
  extern std::vector<NamedFunc> sys_dijet_m_scaledn;
  extern std::vector<NamedFunc> sys_dijet_m_resup;
  extern std::vector<NamedFunc> sys_dijet_m_resdn;
  extern std::vector<NamedFunc> sys_dijet_deta_scaleup;
  extern std::vector<NamedFunc> sys_dijet_deta_scaledn;
  extern std::vector<NamedFunc> sys_dijet_deta_resup;
  extern std::vector<NamedFunc> sys_dijet_deta_resdn;
  extern std::vector<NamedFunc> sys_dijet_dphi_scaleup;
  extern std::vector<NamedFunc> sys_dijet_dphi_scaledn;
  extern std::vector<NamedFunc> sys_dijet_dphi_resup;
  extern std::vector<NamedFunc> sys_dijet_dphi_resdn;

  //jet stuff for pinnacles
  extern const NamedFunc pinnacles_run2_jet_isgood_nopt;
  NamedFunc assign_isgood_pinnacles(NamedFunc& var_pt, 
                                    const std::string &name);
  
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

  //Gets NamedFunc that is llphoton jet dphi with variation
  NamedFunc assign_variation_llphoton_refit_jet_dphi(
      const NamedFunc &llphoton_refit_phi, const NamedFunc &jet_sig, 
      const NamedFunc &jet_pt, const std::string &name);

  //Gets NamedFunc that is balance with variation
  NamedFunc assign_variation_llphoton_refit_jet_balance(
      const NamedFunc &ll_refit_p4, const NamedFunc &var_lead_photon_pt, 
      const NamedFunc &var_lead_photon_phi, const NamedFunc &jet_sig, 
      const NamedFunc &jet_pt, const std::string &name);

  //Gets NamedFunc that is mht (mht, phi) with variation
  NamedFunc assign_variation_mht_balance(
      const NamedFunc &var_el_sig, const NamedFunc &var_el_pt,
      const NamedFunc &var_mu_sig, const NamedFunc &var_mu_pt,
      const NamedFunc &var_ph_sig, const NamedFunc &var_ph_pt,
      const NamedFunc &var_jet_sig, const NamedFunc &var_jet_pt, 
      const std::string &name);

  //Gets NamedFunc that is photon zeppenfeld with variation
  NamedFunc assign_variation_photon_zeppenfeld(
      const NamedFunc &var_lead_photon_eta, const NamedFunc &var_njet, 
      const NamedFunc &var_lead_jet_eta, 
      const NamedFunc &var_sublead_jet_eta, const std::string &name);

  //jet+llphoton variations
  extern const NamedFunc sys_llphoton_refit_dijet_dphi_default;
  extern const NamedFunc sys_llphoton_refit_dijet_dphi_elscaleup;
  extern const NamedFunc sys_llphoton_refit_dijet_dphi_elscaledn;
  extern const NamedFunc sys_llphoton_refit_dijet_dphi_elresup;
  extern const NamedFunc sys_llphoton_refit_dijet_dphi_elresdn;
  extern const NamedFunc sys_llphoton_refit_dijet_dphi_muscaleup;
  extern const NamedFunc sys_llphoton_refit_dijet_dphi_muscaledn;
  extern const NamedFunc sys_llphoton_refit_dijet_dphi_muresup;
  extern const NamedFunc sys_llphoton_refit_dijet_dphi_muresdn;
  extern const NamedFunc sys_llphoton_refit_dijet_dphi_phscaleup;
  extern const NamedFunc sys_llphoton_refit_dijet_dphi_phscaledn;
  extern const NamedFunc sys_llphoton_refit_dijet_dphi_phresup;
  extern const NamedFunc sys_llphoton_refit_dijet_dphi_phresdn;
  extern std::vector<NamedFunc> sys_llphoton_refit_dijet_dphi_jetscaleup;
  extern std::vector<NamedFunc> sys_llphoton_refit_dijet_dphi_jetscaledn;
  extern std::vector<NamedFunc> sys_llphoton_refit_dijet_dphi_jetresup;
  extern std::vector<NamedFunc> sys_llphoton_refit_dijet_dphi_jetresdn;

  extern const NamedFunc sys_llphoton_refit_jet_dphi_default;
  extern const NamedFunc sys_llphoton_refit_jet_dphi_default_elscaleup;
  extern const NamedFunc sys_llphoton_refit_jet_dphi_default_elscaledn;
  extern const NamedFunc sys_llphoton_refit_jet_dphi_default_elresup;
  extern const NamedFunc sys_llphoton_refit_jet_dphi_default_elresdn;
  extern const NamedFunc sys_llphoton_refit_jet_dphi_default_muscaleup;
  extern const NamedFunc sys_llphoton_refit_jet_dphi_default_muscaledn;
  extern const NamedFunc sys_llphoton_refit_jet_dphi_default_muresup;
  extern const NamedFunc sys_llphoton_refit_jet_dphi_default_muresdn;
  extern const NamedFunc sys_llphoton_refit_jet_dphi_default_phscaleup;
  extern const NamedFunc sys_llphoton_refit_jet_dphi_default_phscaledn;
  extern const NamedFunc sys_llphoton_refit_jet_dphi_default_phresup;
  extern const NamedFunc sys_llphoton_refit_jet_dphi_default_phresdn;
  extern std::vector<NamedFunc> sys_llphoton_refit_jet_dphi_jetscaleup;
  extern std::vector<NamedFunc> sys_llphoton_refit_jet_dphi_jetscaledn;
  extern std::vector<NamedFunc> sys_llphoton_refit_jet_dphi_jetresup;
  extern std::vector<NamedFunc> sys_llphoton_refit_jet_dphi_jetresdn;

  extern const NamedFunc sys_llphoton_refit_jet_balance_default;
  extern const NamedFunc sys_llphoton_refit_jet_balance_elscaleup;
  extern const NamedFunc sys_llphoton_refit_jet_balance_elscaledn;
  extern const NamedFunc sys_llphoton_refit_jet_balance_elresup;
  extern const NamedFunc sys_llphoton_refit_jet_balance_elresdn;
  extern const NamedFunc sys_llphoton_refit_jet_balance_muscaleup;
  extern const NamedFunc sys_llphoton_refit_jet_balance_muscaledn;
  extern const NamedFunc sys_llphoton_refit_jet_balance_muresup;
  extern const NamedFunc sys_llphoton_refit_jet_balance_muresdn;
  extern const NamedFunc sys_llphoton_refit_jet_balance_phscaleup;
  extern const NamedFunc sys_llphoton_refit_jet_balance_phscaledn;
  extern const NamedFunc sys_llphoton_refit_jet_balance_phresup;
  extern const NamedFunc sys_llphoton_refit_jet_balance_phresdn;
  extern std::vector<NamedFunc> sys_llphoton_refit_jet_balance_jetscaleup;
  extern std::vector<NamedFunc> sys_llphoton_refit_jet_balance_jetscaledn;
  extern std::vector<NamedFunc> sys_llphoton_refit_jet_balance_jetresup;
  extern std::vector<NamedFunc> sys_llphoton_refit_jet_balance_jetresdn;

  extern const NamedFunc sys_mht_default;
  extern const NamedFunc sys_mht_elscaleup;
  extern const NamedFunc sys_mht_elscaledn;
  extern const NamedFunc sys_mht_elresup;
  extern const NamedFunc sys_mht_elresdn;
  extern const NamedFunc sys_mht_muscaleup;
  extern const NamedFunc sys_mht_muscaledn;
  extern const NamedFunc sys_mht_muresup;
  extern const NamedFunc sys_mht_muresdn;
  extern const NamedFunc sys_mht_phscaleup;
  extern const NamedFunc sys_mht_phscaledn;
  extern const NamedFunc sys_mht_phresup;
  extern const NamedFunc sys_mht_phresdn;
  extern std::vector<NamedFunc> sys_mht_jetscaleup;
  extern std::vector<NamedFunc> sys_mht_jetscaledn;
  extern std::vector<NamedFunc> sys_mht_jetresup;
  extern std::vector<NamedFunc> sys_mht_jetresdn;

  extern const NamedFunc sys_mht_phi_default;
  extern const NamedFunc sys_mht_phi_elscaleup;
  extern const NamedFunc sys_mht_phi_elscaledn;
  extern const NamedFunc sys_mht_phi_elresup;
  extern const NamedFunc sys_mht_phi_elresdn;
  extern const NamedFunc sys_mht_phi_muscaleup;
  extern const NamedFunc sys_mht_phi_muscaledn;
  extern const NamedFunc sys_mht_phi_muresup;
  extern const NamedFunc sys_mht_phi_muresdn;
  extern const NamedFunc sys_mht_phi_phscaleup;
  extern const NamedFunc sys_mht_phi_phscaledn;
  extern const NamedFunc sys_mht_phi_phresup;
  extern const NamedFunc sys_mht_phi_phresdn;
  extern std::vector<NamedFunc> sys_mht_phi_jetscaleup;
  extern std::vector<NamedFunc> sys_mht_phi_jetscaledn;
  extern std::vector<NamedFunc> sys_mht_phi_jetresup;
  extern std::vector<NamedFunc> sys_mht_phi_jetresdn;

  extern const NamedFunc sys_photon_mht_dphi_default;
  extern const NamedFunc sys_photon_mht_dphi_elscaleup;
  extern const NamedFunc sys_photon_mht_dphi_elscaledn;
  extern const NamedFunc sys_photon_mht_dphi_elresup;
  extern const NamedFunc sys_photon_mht_dphi_elresdn;
  extern const NamedFunc sys_photon_mht_dphi_muscaleup;
  extern const NamedFunc sys_photon_mht_dphi_muscaledn;
  extern const NamedFunc sys_photon_mht_dphi_muresup;
  extern const NamedFunc sys_photon_mht_dphi_muresdn;
  extern const NamedFunc sys_photon_mht_dphi_phscaleup;
  extern const NamedFunc sys_photon_mht_dphi_phscaledn;
  extern const NamedFunc sys_photon_mht_dphi_phresup;
  extern const NamedFunc sys_photon_mht_dphi_phresdn;
  extern std::vector<NamedFunc> sys_photon_mht_dphi_jetscaleup;
  extern std::vector<NamedFunc> sys_photon_mht_dphi_jetscaledn;
  extern std::vector<NamedFunc> sys_photon_mht_dphi_jetresup;
  extern std::vector<NamedFunc> sys_photon_mht_dphi_jetresdn;

  extern const NamedFunc sys_photon_jet1_dr_default;
  extern const NamedFunc sys_photon_jet1_dr_phscaleup;
  extern const NamedFunc sys_photon_jet1_dr_phscaledn;
  extern const NamedFunc sys_photon_jet1_dr_phresup;
  extern const NamedFunc sys_photon_jet1_dr_phresdn;
  extern std::vector<NamedFunc> sys_photon_jet1_dr_jetscaleup;
  extern std::vector<NamedFunc> sys_photon_jet1_dr_jetscaledn;
  extern std::vector<NamedFunc> sys_photon_jet1_dr_jetresup;
  extern std::vector<NamedFunc> sys_photon_jet1_dr_jetresdn;

  extern const NamedFunc sys_photon_jet2_dr_default;
  extern const NamedFunc sys_photon_jet2_dr_phscaleup;
  extern const NamedFunc sys_photon_jet2_dr_phscaledn;
  extern const NamedFunc sys_photon_jet2_dr_phresup;
  extern const NamedFunc sys_photon_jet2_dr_phresdn;
  extern std::vector<NamedFunc> sys_photon_jet2_dr_jetscaleup;
  extern std::vector<NamedFunc> sys_photon_jet2_dr_jetscaledn;
  extern std::vector<NamedFunc> sys_photon_jet2_dr_jetresup;
  extern std::vector<NamedFunc> sys_photon_jet2_dr_jetresdn;

  extern const NamedFunc sys_photon_zeppenfeld_default;
  extern const NamedFunc sys_photon_zeppenfeld_phscaleup;
  extern const NamedFunc sys_photon_zeppenfeld_phscaledn;
  extern const NamedFunc sys_photon_zeppenfeld_phresup;
  extern const NamedFunc sys_photon_zeppenfeld_phresdn;
  extern std::vector<NamedFunc> sys_photon_zeppenfeld_jetscaleup;
  extern std::vector<NamedFunc> sys_photon_zeppenfeld_jetscaledn;
  extern std::vector<NamedFunc> sys_photon_zeppenfeld_jetresup;
  extern std::vector<NamedFunc> sys_photon_zeppenfeld_jetresdn;

  //btag uncorrelated variations
  extern std::vector<NamedFunc> sys_bchig_uncorr_up;
  extern std::vector<NamedFunc> sys_bchig_uncorr_dn;
  extern std::vector<NamedFunc> sys_udsghig_uncorr_up;
  extern std::vector<NamedFunc> sys_udsghig_uncorr_dn;

  //horrible hack to get around the fact that I can't put a for loop outside
  //of a function to initialize what should be constant vectors...
  void initialize_jetvariations();

  //old (25-03) BDTs
  std::vector<std::shared_ptr<MVAWrapper>> SystVbfBdts();
  //returns thread safe NamedFunc that returns VBF score
  //NamedFunc syst_vbf_bdt_score(const char* variation="");
  NamedFunc syst_vbf_bdt_score(std::string variation="");

  NamedFunc ggfbdt2503_score_default(
      const std::vector<fastforest::FastForest> &xgb_bdts);
  NamedFunc ggfbdt2503_score_elscaleup(
      const std::vector<fastforest::FastForest> &xgb_bdts);
  NamedFunc ggfbdt2503_score_elscaledn(
      const std::vector<fastforest::FastForest> &xgb_bdts);
  NamedFunc ggfbdt2503_score_elresup(
      const std::vector<fastforest::FastForest> &xgb_bdts);
  NamedFunc ggfbdt2503_score_elresdn(
      const std::vector<fastforest::FastForest> &xgb_bdts);
  NamedFunc ggfbdt2503_score_muscaleup(
      const std::vector<fastforest::FastForest> &xgb_bdts);
  NamedFunc ggfbdt2503_score_muscaledn(
      const std::vector<fastforest::FastForest> &xgb_bdts);
  NamedFunc ggfbdt2503_score_muresup(
      const std::vector<fastforest::FastForest> &xgb_bdts);
  NamedFunc ggfbdt2503_score_muresdn(
      const std::vector<fastforest::FastForest> &xgb_bdts);
  NamedFunc ggfbdt2503_score_phscaleup(
      const std::vector<fastforest::FastForest> &xgb_bdts);
  NamedFunc ggfbdt2503_score_phscaledn(
      const std::vector<fastforest::FastForest> &xgb_bdts);
  NamedFunc ggfbdt2503_score_phresup(
      const std::vector<fastforest::FastForest> &xgb_bdts);
  NamedFunc ggfbdt2503_score_phresdn(
      const std::vector<fastforest::FastForest> &xgb_bdts);
}

#endif
