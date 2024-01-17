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

#include "core/gamma_params.hpp"
#include "core/mva_wrapper.hpp"
#include "core/named_func.hpp"
#include "core/named_func.hpp"
#include "core/process.hpp"
#include "core/plot_maker.hpp"
#include "core/sample_loader.hpp"
#include "core/table.hpp"
#include "core/table_row.hpp"

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

  //returns working version of kinematic BDT
  std::shared_ptr<MVAWrapper> KinematicBdt();
  //returns NamedFunc that selects low BDT score category "ggF/untagged 4"
  NamedFunc category_ggh4(std::shared_ptr<MVAWrapper> kinematic_bdt);
  //returns NamedFunc that selects medium BDT score category "ggF/untagged 3"
  NamedFunc category_ggh3(std::shared_ptr<MVAWrapper> kinematic_bdt);
  //returns NamedFunc that selects high BDT score category "ggF/untagged 2"
  NamedFunc category_ggh2(std::shared_ptr<MVAWrapper> kinematic_bdt);
  //returns NamedFunc that selects very high BDT score category "ggF/untagged 1"
  NamedFunc category_ggh1(std::shared_ptr<MVAWrapper> kinematic_bdt);

  //returns a sample loader that has the H->Zy colors pre-sets and NamedFuncs loaded
  SampleLoader ZgSampleLoader();
  int get_btag_wp_deepjet(int year, float discriminator_value);

  //returns a list of just background processes from a list of processes
  std::vector<std::shared_ptr<Process>> GetBackgroundProcesses(
      std::vector<std::shared_ptr<Process>> processes);

  //returns a list of just data processes from a list of processes
  std::vector<std::shared_ptr<Process>> GetDataProcesses(
      std::vector<std::shared_ptr<Process>> processes);

  //sets all processes to background, useful for making colz 2D plots of data
  void SetProcessesBackground(
      std::vector<std::shared_ptr<Process>> &processes);
}
#endif
