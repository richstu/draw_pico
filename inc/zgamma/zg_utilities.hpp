#ifndef H_ZG_UTILITIES
#define H_ZG_UTILITIES

#include <iostream>
#include <string>
#include <utility>
#include <set>
#include <vector>
#include <map>
#include <memory>

#include "TString.h"
#include "TLorentzVector.h"

#include "core/named_func.hpp"
#include "core/named_func.hpp"
#include "core/process.hpp"
#include "core/plot_maker.hpp"
#include "core/sample_loader.hpp"
#include "core/table.hpp"
#include "core/table_row.hpp"
#include "core/gamma_params.hpp"

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


}
#endif
