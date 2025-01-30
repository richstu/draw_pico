#ifndef H_ZG_FUNCTIONS
#define H_ZG_FUNCTIONS

#include "core/named_func.hpp"

namespace ZgFunctions {
  //isolated dielectron triggers for run 2
  extern const NamedFunc HLT_pass_dielectron;

  //isolated dimuon triggers for run 2
  extern const NamedFunc HLT_pass_dimuon;

  //isolated dilepton triggers for run 2
  extern const NamedFunc HLT_pass_dilepton;

  //isolated single electron triggers for run 2
  extern const NamedFunc HLT_pass_singleelectron;

  //isolated single muon triggers for run 2
  extern const NamedFunc HLT_pass_singlemuon;

  //isolated single lepton triggers for run 2
  extern const NamedFunc HLT_pass_singlelepton;

  //year integrated lumi weights
  extern const NamedFunc w_years;

  //Run 3 weight-up
  extern const NamedFunc w_run3;

  //leading muon pt
  extern const NamedFunc lead_mu_pt;

  //subleading muon pt
  extern const NamedFunc sublead_mu_pt;

  //leading electron pt
  extern const NamedFunc lead_el_pt;

  //subleading electron pt
  extern const NamedFunc sublead_el_pt;

  //leading jet pt
  extern const NamedFunc lead_jet_pt;

  //subleading jet pt
  extern const NamedFunc sublead_jet_pt;

  //leading jet eta
  extern const NamedFunc lead_jet_eta;

  //subleading jet eta
  extern const NamedFunc sublead_jet_eta;

  //leading jet phi
  extern const NamedFunc lead_jet_phi;

  //subleading jet phi
  extern const NamedFunc sublead_jet_phi;

  //leading jet m
  extern const NamedFunc lead_jet_m;

  //subleading jet m
  extern const NamedFunc sublead_jet_m;

  //trigger-paired lepton pt plateau cuts
  extern const NamedFunc trig_plateau_cuts;

  //common NamedFuncs for run 2 baseline selection
  extern const NamedFunc zg_baseline_el_run2;
  extern const NamedFunc zg_baseline_mu_run2;
  extern const NamedFunc zg_baseline_run2;

  //new working baseline selection
  extern const NamedFunc zg_baseline;

  //master stitch variable
  extern const NamedFunc stitch;

  //master stitch variable, for deathvalley productions
  extern const NamedFunc stitch_deathvalley;

  //drmax of lead photon
  extern const NamedFunc photon_drmax;

  //relative pt uncertainty of lead photon for kingscanyon_v0 productions and earlier
  extern const NamedFunc photon_relpterr_deathvalley;

  //relative pt uncertainty of lead photon
  extern const NamedFunc photon_relpterr;

  //lead lepton eta (=lep_eta[0], but this isn't saved in slims =( )
  //only works for 2 lepton events
  extern const NamedFunc lead_lepton_eta;

  //sublead lepton eta (=lep_eta[0], but this isn't saved in slims =( )
  //only works for 2 lepton events
  extern const NamedFunc sublead_lepton_eta;

  //pT/m of Higgs candidate
  extern const NamedFunc llphoton_rel_pt;

  //modified ptt variable
  extern const NamedFunc llphoton_pttmod;

  //maximum lepton mini isolation
  extern const NamedFunc max_lep_miniso;

  //vector of whether GenParticles are the first copy
  extern const NamedFunc mc_isFirstCopy;

  //vector mother PDGID of mother of GenParticle
  extern const NamedFunc mc_mommom;
}

#endif
