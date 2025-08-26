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
  std::vector<double> CalculateAngles(TLorentzVector lplus, 
      TLorentzVector lminus, TLorentzVector ph);

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
  NamedFunc vbf_bdt_score(std::string variation="");
  //returns NamedFunc that selects very high BDT score VBF category 
  NamedFunc category_vbf1(const NamedFunc &bdt_scores);
  //returns NamedFunc that selects high BDT score VBF category 
  NamedFunc category_vbf2(const NamedFunc &bdt_scores);
  //returns NamedFunc that selects medium BDT score VBF category
  NamedFunc category_vbf3(const NamedFunc &bdt_scores);
  //returns NamedFunc that selects low BDT score VBF category 
  NamedFunc category_vbf4(const NamedFunc &bdt_scores);

  //returns XGBoost BDTs
  const std::vector<fastforest::FastForest> XGBoostBDTs();
  //Returns NamedFunc that returns XGBoost score
  NamedFunc OldXGBoostBDTScore(
      const std::vector<fastforest::FastForest> &xgb_bdts);
  //Returns NamedFunc that returns XGBoost score
  NamedFunc XGBoostBDTScore(
      const std::vector<fastforest::FastForest> &xgb_bdts, 
      const std::vector<NamedFunc> &inputs);
  //Returns NamedFunc that returns XGBoost score
  //Be careful: neither NamedFuncs nor xgb_bdt may be allowed to go out of 
  //scope
  NamedFunc XGBoostBDTScoreCached(
      const std::vector<fastforest::FastForest> &xgb_bdts, 
      const std::vector<const NamedFunc*> &inputs);
  //returns NamedFunc that selects low BDT score category "ggF 4"
  NamedFunc category_ggf4(const NamedFunc &bdtscore);
  //returns NamedFunc that selects medium BDT score category "ggF 3"
  NamedFunc category_ggf3(const NamedFunc &bdtscore);
  //returns NamedFunc that selects high BDT score category "ggF 2"
  NamedFunc category_ggf2(const NamedFunc &bdtscore);
  //returns NamedFunc that selects very high BDT score category "ggF 1"
  NamedFunc category_ggf1(const NamedFunc &bdtscore);

  //returns a sample loader that has the H->Zy colors pre-sets and NamedFuncs loaded
  SampleLoader ZgSampleLoader();

  void rename_signal(std::vector<std::shared_ptr<Process>> &procs, int factor);
  std::vector<std::shared_ptr<Process>> procs_with_sig_scale(const std::string file_name, const std::string loader_name, int sig_factor);

  int get_btag_wp_deepjet(int year, float discriminator_value);

  //Kinematic Refit Functions
  bool isFSRphoton(const Baby &b);
  std::map<unsigned int, TLorentzVector> fsrphoton_ret(const Baby &b);
  std::map<unsigned int, TLorentzVector> fsrphoton_ret_customll(const Baby &b,
      int lepid, int i1, int i2);
  double KinRefit(const Baby &b);
  std::vector<TLorentzVector> RefitP4(const Baby &b);
  double AssignL1Error(const Baby &b);
  double AssignL2Error(const Baby &b);
  double difference_check(const Baby &b);
  double difference_check_lly(const Baby &b);

  //Functions to test the mu_correctedPt
  TLorentzVector AssignCorrL1(const Baby &b);
  TLorentzVector AssignCorrL2(const Baby &b);
  double KinRefitCorrected(const Baby &b);
  std::vector<TLorentzVector> RefitP4Corr(const Baby &b); 

  //returns WP. 1 is loose, 2 is medium, 3 is tight
  float get_btag_wp_deepjet(const std::string& year, int wp);

  //sets all processes to background, useful for making colz 2D plots of data
  void SetProcessesBackground(
      std::vector<std::shared_ptr<Process>> &processes);

  //returns lepton (pt1, eta1, phi1, m1, pt2, eta2, phi2, m2) with custom refit
  std::vector<double> get_lep_custom_refit(const Baby &b, 
      NamedFunc el_pt, NamedFunc mu_pt, int ll_lepid, int ll_i1, int ll_i2);
}
#endif
