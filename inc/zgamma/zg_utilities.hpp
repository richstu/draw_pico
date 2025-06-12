#ifndef H_ZG_UTILITIES
#define H_ZG_UTILITIES

#include <iostream>
#include <string>
#include <utility>
#include <set>
#include <vector>
#include <map>
#include <memory>

#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "TString.h"
#include "TLorentzVector.h"

#include "core/fastforest.hpp"
#include "core/gamma_params.hpp"
#include "core/mva_wrapper.hpp"
#include "core/named_func.hpp"
#include "core/named_func.hpp"
#include "core/process.hpp"
#include "core/plot_maker.hpp"
#include "core/sample_loader.hpp"
#include "core/table.hpp"
#include "core/table_row.hpp"
#include "zgamma/KinZfitter.hpp"

namespace ZgUtilities {
  TLorentzVector AssignL1(const Baby &b, bool gen = false);
  TLorentzVector AssignL2(const Baby &b, bool gen = false);
  TLorentzVector AssignZ (const Baby &b, bool gen = false);
  TLorentzVector AssignH (const Baby &b, bool gen = false);
  TLorentzVector AssignQ1(const Baby &b, bool gen = false);
  TLorentzVector AssignQ2(const Baby &b, bool gen = false);
  TLorentzVector AssignJJ(const Baby &b, bool gen = false);
  TLorentzVector AssignGamma (const Baby &b, bool gen = false);
  TLorentzVector FindGamma (const Baby &b, bool gen = false);
  double Findlly(const Baby &b);
  double lambdaZ(const Baby &b, bool gen = false);
  double cos_theta(const Baby &b, bool gen = false);
  double cos_Theta(const Baby &b, bool gen = false);
  double Getphi(const Baby &b, bool gen = false);
  double pdrmax(const Baby &b);

  //correct angle functions from Jaebak
  TLorentzVector get_q1(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma);
  TLorentzVector get_q2(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma);
  double get_lambdaZ(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma);
  double get_cosTheta(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma);
  double get_costheta(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma);
  double get_phi(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma);

  //returns working version of kinematic BDT
  std::shared_ptr<MVAWrapper> KinematicBdt();
  //returns NamedFunc that selects low BDT score category "ggF/untagged 4"
  NamedFunc category_ggh4_old(std::shared_ptr<MVAWrapper> kinematic_bdt);
  //returns NamedFunc that selects medium BDT score category "ggF/untagged 3"
  NamedFunc category_ggh3_old(std::shared_ptr<MVAWrapper> kinematic_bdt);
  //returns NamedFunc that selects high BDT score category "ggF/untagged 2"
  NamedFunc category_ggh2_old(std::shared_ptr<MVAWrapper> kinematic_bdt);
  //returns NamedFunc that selects very high BDT score category "ggF/untagged 1"
  NamedFunc category_ggh1_old(std::shared_ptr<MVAWrapper> kinematic_bdt);

  //returns working version of dijet BDT
  std::vector<std::shared_ptr<MVAWrapper>> VbfBdts();
  //returns NamedFunc that returns VBF score
  NamedFunc vbf_bdt_score(std::vector<std::shared_ptr<MVAWrapper>> vbf_bdts);
  //returns NamedFunc that selects very high BDT score VBF category 
  NamedFunc category_vbf1(std::vector<std::shared_ptr<MVAWrapper>> vbf_bdts);
  //returns NamedFunc that selects high BDT score VBF category 
  NamedFunc category_vbf2(std::vector<std::shared_ptr<MVAWrapper>> vbf_bdts);
  //returns NamedFunc that selects medium BDT score VBF category
  NamedFunc category_vbf3(std::vector<std::shared_ptr<MVAWrapper>> vbf_bdts);
  //returns NamedFunc that selects low BDT score VBF category 
  NamedFunc category_vbf4(std::vector<std::shared_ptr<MVAWrapper>> vbf_bdts);

  //returns XGBoost BDTs
  const std::vector<fastforest::FastForest> XGBoostBDTs();
  //Returns NamedFunc that returns XGBoost score
  NamedFunc XGBoostBDTScore(
      const std::vector<fastforest::FastForest> &xgb_bdts);
  //returns NamedFunc that selects low BDT score category "ggF 4"
  NamedFunc category_ggf4(const std::vector<fastforest::FastForest> &xgb_bdts);
  //returns NamedFunc that selects medium BDT score category "ggF 3"
  NamedFunc category_ggf3(const std::vector<fastforest::FastForest> &xgb_bdts);
  //returns NamedFunc that selects high BDT score category "ggF 2"
  NamedFunc category_ggf2(const std::vector<fastforest::FastForest> &xgb_bdts);
  //returns NamedFunc that selects very high BDT score category "ggF 1"
  NamedFunc category_ggf1(const std::vector<fastforest::FastForest> &xgb_bdts);

  //returns a sample loader that has the H->Zy colors pre-sets and NamedFuncs loaded
  SampleLoader ZgSampleLoader();

  void rename_signal(std::vector<std::shared_ptr<Process>> &procs, int factor);
  std::vector<std::shared_ptr<Process>> procs_with_sig_scale(const std::string file_name, const std::string loader_name, int sig_factor);

  int get_btag_wp_deepjet(int year, float discriminator_value);

  //Kinematic Refit Functions
  bool isFSRphoton(const Baby &b);
  std::map<unsigned int, TLorentzVector> fsrphoton_ret(const Baby &b);
  double KinRefit(const Baby &b);
  double KinRefit(const Baby &b,TString txtFile);
  std::vector<TLorentzVector> RefitP4(const Baby &b);
  std::vector<TLorentzVector> RefitP4(const Baby &b,TString txtFile);
  double AssignL1Error(const Baby &b);
  double AssignL2Error(const Baby &b);
  double difference_check(const Baby &b,TString txtFile);
  double difference_check_lly(const Baby &b,TString txtFile);

  //Functions to test the mu_correctedPt
  TLorentzVector AssignCorrL1(const Baby &b);
  TLorentzVector AssignCorrL2(const Baby &b);
  double KinRefitCorrected(const Baby &b, TString txtFile);
  std::vector<TLorentzVector> RefitP4Corr(const Baby &b, TString txtFile); 

  //returns WP. 1 is loose, 2 is medium, 3 is tight
  float get_btag_wp_deepjet(const std::string& year, int wp);

  //sets all processes to background, useful for making colz 2D plots of data
  void SetProcessesBackground(
      std::vector<std::shared_ptr<Process>> &processes);

  //returns three-body invariant mass with custom refit
  //TODO update this
  std::vector<double> get_lep_pt_custom_refit(const Baby &b, 
      std::shared_ptr<KinZfitter> kinZfitter, NamedFunc el_pt_var, 
      NamedFunc mu_pt_var);
}
#endif
