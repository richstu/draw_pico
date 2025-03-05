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

  TLorentzVector AssignJet(const Baby &b,int jet_num);
  TLorentzVector AssignLep(const Baby &b, int flav, int idx);
  float get_dphi(float phi1, float phi2);
  float dphi_Hmet_d(const Baby &b);
  double max_deepflav_d(const Baby &b);
  int el_firstIdx(const Baby &b, std::vector<int> idces);
  int mu_firstIdx(const Baby &b, std::vector<int> idces);
  int el3_idx(const Baby &b);
  int mu3_idx(const Baby &b);
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
  bool checkBit(int i, int n);

  size_t idx_j1(const Baby &b);
  size_t idx_j2(const Baby &b);
  size_t idx_jN(const Baby &b, int num_jet);

  extern const NamedFunc mbb;
  extern const NamedFunc mtop;
  extern const NamedFunc lly_pT;
  extern const NamedFunc lly_dR;
  extern const NamedFunc pTlly_mlly;
  extern const NamedFunc mlly;
  extern const NamedFunc l3_pt;
  extern const NamedFunc l3_idmva;
  extern const NamedFunc l3_mini;
  extern const NamedFunc l3_flav;
  extern const NamedFunc l4_pt;
  extern const NamedFunc l4_idmva;
  extern const NamedFunc l4_mini;
  extern const NamedFunc l4_flav;
  extern const NamedFunc max_Imini;
  extern const NamedFunc l1_eta;
  extern const NamedFunc l2_eta;

  extern const NamedFunc ptj_m4j;
  extern const NamedFunc max_deepflav;
  extern const NamedFunc ptj_m5j;

  extern const NamedFunc pTy_mlly;
  extern const NamedFunc mlly;
  extern const NamedFunc pT_lly;
  extern const NamedFunc eta_lly;
  extern const NamedFunc phi_lly;
  extern const NamedFunc sigy_pTy;
  extern const NamedFunc diphi_dijet;

  extern const NamedFunc j1_pt;
  extern const NamedFunc j2_pt;
  extern const NamedFunc j3_pt;
  extern const NamedFunc j4_pt;

  extern const NamedFunc j1_eta;
  extern const NamedFunc j2_eta;
  extern const NamedFunc j3_eta;
  extern const NamedFunc j4_eta;

  extern const NamedFunc j1_phi;
  extern const NamedFunc j2_phi;
  extern const NamedFunc j3_phi;
  extern const NamedFunc j4_phi;

  extern const NamedFunc j1_m;
  extern const NamedFunc j2_m;
  extern const NamedFunc j3_m;
  extern const NamedFunc j4_m;

  //functions of categories' object selections
  int R2_Cats(const Baby &b);
  int cat_obj_sels(const Baby &b);

  //Run 2 categories
  extern const NamedFunc run2_lep;
  extern const NamedFunc run2_VBF;
  extern const NamedFunc run2_ggF;

  //Vector containing the categories used for Run 2 
  extern const std::vector<NamedFunc> run2_categories;
  extern const std::vector<std::string> run2_cat_labs;

  extern const NamedFunc cat_ttH_lep;
  extern const NamedFunc cat_WH_3l;
  extern const NamedFunc cat_ttH_had;
  extern const NamedFunc cat_ZH_met;
  extern const NamedFunc cat_VBF;
  extern const NamedFunc cat_ggF;


  //Vector returning all the categories used for Run 3
  extern const std::vector<NamedFunc>   run3_category_vector;
  extern const std::vector<NamedFunc>   run3_catwsel_vector;
  extern const std::vector<NamedFunc>   run3_refit_catwsel_vector;
  extern const std::vector<NamedFunc>   run3_catbs_vector;
  extern const std::vector<NamedFunc>   run3_catbs_nomll_vector; 
  extern const std::vector<NamedFunc>   run3_catwsel_nomll_vector;
  extern const std::vector<std::string> run3_category_labels;
  extern const std::vector<std::vector<NamedFunc>> categories_selections_vector; 

  //Each category's sample plots
  void ttH_lep_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels, bool run3=false);
  void   ZH_4l_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels);
  void   WH_3l_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels);
  void ttH_had_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels, bool run3=false);
  void  ZH_2bl_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels);
  void  ZH_MET_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels);
  void     ggF_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels);
  void             ggF_input_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels);
  void             VBF_input_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels);
  void     VBF_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels);
  void     Lep_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels);
  void       sample_kinrefit_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels);

  void ttH_lep_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels, int nbins, bool run3=false);
  void   WH_3l_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels, int nbins);
  void ttH_had_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels, int nbins, bool run3=false);
  void  ZH_MET_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels, int nbins);
  void     ggF_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels, int nbins);
  void     VBF_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels, int nbins);
  void       sample_kinrefit_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, int nbins, std::string labels);


  //MVP plots
  void                 mvp_objects(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels);
  void   mvp_procs_and_corrections(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels);

  std::vector<std::vector<float>> plotting_bins(bool plotting_run3);
}
#endif
