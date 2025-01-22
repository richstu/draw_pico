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

  std::vector<NamedFunc> progressive_cuts(std::vector<NamedFunc> vector_NamedFunc);
  NamedFunc Nreverse1(std::vector<NamedFunc> vector_NamedFunc, unsigned int reverse);
  NamedFunc Nreplace1(std::vector<NamedFunc> vector_NamedFunc, NamedFunc replace, unsigned int skip);
  NamedFunc Nminus1(  std::vector<NamedFunc> vector_NamedFunc, unsigned int skip);
  NamedFunc Nminusk(  std::vector<NamedFunc> vector_NamedFunc, std::vector<unsigned int> skip);

  //This is a vector that gives the individual selections in the baseline selection
  extern const NamedFunc pass_trigs;

  extern const NamedFunc pass_trigs_and_pt;
  extern const std::vector<NamedFunc> vector_run2_baseline;
  extern const std::vector<NamedFunc> vector_tightened_baseline;
  extern const std::vector<NamedFunc> vtb_refit;

  //Below is the baseline used by HIG-19-014
  extern const NamedFunc hig19014_baseline;

  //Below is the tightened baseline with the photon wp80 and the 80 GeV < mll < 100 GeV selection
  extern const NamedFunc tightened_baseline;
  extern const NamedFunc tightened_baseline_refit;
  
  extern const NamedFunc tightened_baseline_pinnacles;

  //Below functions are related to weighting of samples or naming based off of weighting
  extern const NamedFunc wgt;
  extern const NamedFunc wgt_pin_fix;
  
  NamedFunc add_cut(NamedFunc & current_cut, NamedFunc additional_cut);

}

#endif
