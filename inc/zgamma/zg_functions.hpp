#ifndef H_ZG_FUNCTIONS
#define H_ZG_FUNCTIONS

#include "core/named_func.hpp"
#include "core/plot_maker.hpp"

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

  //run 3 weight-up
  extern const NamedFunc w_run3;

  //weights all Run 2 MC to all 1 fb^-1
  extern const NamedFunc wgt_1invfb;

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

  //Turns a vector<NamedFunc> into one usable for a cutflow table
  std::vector<NamedFunc> progressive_cuts(std::vector<NamedFunc> vector_NamedFunc);

  //This function adds all selections and the selection at skip
  NamedFunc Nreverse1(std::vector<NamedFunc> vector_NamedFunc, unsigned int reverse);

  //Returns a NamedFunc with all but one selection (marked by skip)
  NamedFunc Nminus1(std::vector<NamedFunc> vector_NamedFunc, unsigned int skip);

  //Returns a NamedFunc without all selections in the vector skip
  NamedFunc Nminusk(std::vector<NamedFunc> vector_NamedFunc, std::vector<unsigned int> skip);

  extern const NamedFunc mlly;
  extern const NamedFunc pTy_mlly;
  extern const NamedFunc mll_mlly;
  extern const NamedFunc cosTheta;
  extern const NamedFunc costheta;
  extern const NamedFunc Phi;


  //adds ggF control region plots to PlotMaker
  void ggF_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels, int ll_flavor=0);

  //adds VBF control region plots to PlotMaker
  void VBF_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> processes, std::vector<PlotOpt> ops, NamedFunc wgt, std::string labels, int ll_flavor=0);

  //adds lep control region plots to PlotMaker
  void Lep_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> processes, std::vector<PlotOpt> ops, NamedFunc wgt, std::string labels, int ll_flavor=0);

  extern const NamedFunc pass_trigs;
  extern const NamedFunc pass_trigs_and_pt;

  //Different forms of the baseline selection for plots or N-1 plots
  extern const std::vector<NamedFunc> vector_run2_baseline;  
  extern const std::vector<NamedFunc> vector_tightened_baseline;
  extern const NamedFunc run2_baseline;
  extern const NamedFunc tightened_baseline;

}

#endif
