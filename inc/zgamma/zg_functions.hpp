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

  //common NamedFuncs for run 2 baseline selection
  extern const NamedFunc zg_baseline_el;
  extern const NamedFunc zg_baseline_mu;
  extern const NamedFunc zg_baseline;

  //master stitch variable, currently only has DY/ZG stitch
  extern const NamedFunc stitch;

  //drmax of lead photon
  extern const NamedFunc photon_drmax;

  //lead lepton eta (=lep_eta[0], but this isn't saved in slims =( )
  //only works for 2 lepton events
  extern const NamedFunc lead_lepton_eta;

  //sublead lepton eta (=lep_eta[0], but this isn't saved in slims =( )
  //only works for 2 lepton events
  extern const NamedFunc sublead_lepton_eta;

  //pT/m of Higgs candidate
  extern const NamedFunc llphoton_rel_pt;
}

#endif
