#ifndef CAT_UTILITIES
#define CAT_UTILITIES

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
#include "core/process.hpp"
#include "core/gamma_params.hpp"

namespace CatUtilities {

  //Test Zone

  int nlep_dxyz(const Baby &b);
  bool btag_DF_pass(TString year,float deep_flav,int lmt);
  TLorentzVector bjet_cand(const Baby &b,int lmt);
  TLorentzVector AssignJet(const Baby &b,int jet_num);
  TLorentzVector AssignLep(const Baby &b, int flav, int idx);
  float get_dphi(float phi1, float phi2);
  float dphi_Hmet_d(const Baby &b);
  double max_deepflav_d(const Baby &b);
  int el_firstIdx(const Baby &b, std::vector<int> idces);
  int mu_firstIdx(const Baby &b, std::vector<int> idces);
  int el3_idx(const Baby &b);
  int mu3_idx(const Baby &b);
  int Z2_idx(const Baby &b);
  std::vector<int> L3_Idx(const Baby &b);
  TLorentzVector AssignL3(const Baby &b);
  TLorentzVector AssignL4(const Baby &b);
  int L3_flav_i(const Baby &b);
  int L4_flav_i(const Baby &b);
  double l3_pt_d(const Baby &b);
  double l4_pt_d(const Baby &b);
  double l3_ZCheck_d(const Baby &b);
  double l3_Z_DR_d(const Baby &b);
  double l3_idmva_d(const Baby &b);
  double l4_idmva_d(const Baby &b);
  double l3_Imini_d(const Baby &b);
  double l4_Imini_d(const Baby &b);
  double l3_mT(const Baby &b);
  double Z2_pt_d(const Baby &b);
  double Z2_m_d(const Baby &b);
  bool   Z2_flav_b(const Baby &b);
  double Z2y_dR_d(const Baby &b);
  double ZZ_dR_d(const Baby &b);
  double bjetl_drmin3_d(const Baby &b);
  double m_top(const Baby &b);

  extern const NamedFunc mbb;
  extern const NamedFunc mtop;
  extern const NamedFunc lly_pT;
  extern const NamedFunc lly_dR;
  extern const NamedFunc pTlly_mlly;
  extern const NamedFunc l3_pt;
  extern const NamedFunc l3_idmva;
  extern const NamedFunc l3_mini;
  extern const NamedFunc l3_flav;
  extern const NamedFunc l4_pt;
  extern const NamedFunc l4_idmva;
  extern const NamedFunc l4_mini;
  extern const NamedFunc l4_flav;

  extern const NamedFunc l3_ZCheck;
  extern const NamedFunc l3_Z_DR;

  extern const NamedFunc ptj_m4j;
  extern const NamedFunc max_deepflav;
  extern const NamedFunc Z2_pt;
  extern const NamedFunc Z2_m;
  extern const NamedFunc Z2_flav;
  extern const NamedFunc Z2y_dR;
  extern const NamedFunc ZZ_dR;
  extern const NamedFunc dphi_Hmet; 

  //functions of categories' object selections
  int R2_Cats(const Baby &b);
  int cat_obj_sels(const Baby &b);

  //Run 2 categories
  extern const NamedFunc run2_lep;
  extern const NamedFunc run2_VBF;
  extern const NamedFunc run2_ggF;

  //Vector containing the categories used for Run 2 
  extern const std::vector<NamedFunc> run2_categories;

  //working approach to categories with new VH/ttH categories
  extern const NamedFunc cat_ttH_lep_4l;
  extern const NamedFunc cat_ttH_lep_3l2b;
  extern const NamedFunc cat_ttH_lep_3l1b;
  extern const NamedFunc cat_VH_4l;
  extern const NamedFunc cat_VH_3l;
  extern const NamedFunc cat_ttH_had_2b;
  extern const NamedFunc cat_ttH_had_1b;
  extern const NamedFunc cat_ZH_had;
  extern const NamedFunc cat_ZH_met;
  extern const NamedFunc cat_VBF;
  extern const NamedFunc cat_ggF;

  //Vector returning all the categories used for Run 3
  extern const std::vector<NamedFunc>   run3_category_vector;
  extern const std::vector<std::string> run3_category_labels;

  //Each category's sample plots
  void ttH_lep_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels);
  void VH_4l_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels);
  void VH_3l_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels);
  void ttH_had_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels);
  void VH_2bl_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels);
  void VH_MET_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels);

}
#endif
