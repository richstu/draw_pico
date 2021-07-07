#ifndef H_RDF_HIG_FUNCTIONS
#define H_RDF_HIG_FUNCTIONS

#include <cstddef>
#include <string>
#include <vector>

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RVec.hxx"
#include "core/rdf_named_func.hpp"
 
namespace RdfHigfuncs{

  template <class C>
  using RVec = ROOT::VecOps::RVec<C>;

  using ProcessedDataFrame = ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void>;
  
  //misc
  extern const std::vector<std::string> lead_signal_lepton_pt_args;
  float lead_signal_lepton_pt (RVec<bool> el_sig, RVec<float> el_pt, 
                               RVec<bool> mu_sig, RVec<float> mu_pt);
  extern const std::vector<std::string> hig_bcat_args;
  int hig_bcat(int nbl, int nbm, int nbt);
  extern const std::string hig_cand_am_0;
  extern const std::string hig_cand_dm_0;
  extern const std::string hig_cand_drmax_0;

  //filters
  extern const std::vector<std::string> pass_ecalnoisejetfilter_args;
  bool pass_ecalnoisejetfilter(RVec<float> jet_pt, RVec<float> jet_eta, 
                         RVec<float> Jet_phi, float met_phi);
  extern const std::vector<std::string> pass_hemveto_args;
  bool pass_hemveto(int type, int run, long event, RVec<float> el_pt, RVec<float> el_miniso,
           RVec<float> el_eta, RVec<float> el_phi, RVec<float> jet_pt, RVec<float> jet_eta,
           RVec<float> jet_phi, float met_phi);
  extern const std::vector<std::string> pass_filters_args;
  bool pass_filters(bool pass_goodv, bool pass_hbhe, bool pass_hbheiso, bool pass_ecaldeadcell, 
      bool pass_badpfmu, bool pass_muon_jet, int type, bool pass_eebadsc, 
      bool pass_cschalo_tight, bool pass_low_neutral_jet, bool pass_htratio_dphi_tight,
      bool pass_jets, bool pass_ecalnoisejetfilter, long event, bool pass_hemveto);
  extern const std::string final_pass_filters;
  extern const std::string final_ttbar_pass_filters;
  extern const std::string final_zll_pass_filters;
  extern const std::string final_qcd_pass_filters;

  //weights
  extern const std::vector<std::string> w_pileup_args;
  float w_pileup(int type, int npu_tru_mean);
  extern const std::vector<std::string> w_run2_args;
  float w_run2(int type);
  extern const std::vector<std::string> eff_higtrig_run2_args;
  float eff_higtrig_run2(int type, int nvlep, int nel, int nmu, float met, 
                         float ht, float lead_signal_lepton_pt);
  extern const std::string final_weight_run2;
  extern const std::string final_weight_run2_data;

  //triggers
  extern const std::string jet_trigger;
  extern const std::string met_trigger;
  extern const std::string el_trigger;
  extern const std::string mu_trigger;
  extern const std::string ht_trigger;
  extern const std::string jetht_trigger;

  //functions to add columns to RDataFrame
  ProcessedDataFrame rdf_hig_functions(ProcessedDataFrame &data_frame);
  ProcessedDataFrame rdf_hig_functions_slim_mc(ProcessedDataFrame &data_frame);
  ProcessedDataFrame rdf_hig_functions_slim_data(ProcessedDataFrame &data_frame);

}

#endif //H_RDF_HIG_FUNCTIONS
