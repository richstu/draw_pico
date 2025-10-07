/**
 * Script to write datacard for H->Zgamma analysis
 */

#include <cmath>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "TError.h"
#include "TColor.h"

#include "core/axis.hpp"
#include "core/baby_pico.hpp"
#include "core/datacard.hpp"
#include "core/fastforest.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/process.hpp"
#include "core/utilities.hpp"
#include "zgamma/zg_functions.hpp"
#include "zgamma/zg_syst_functions.hpp"
#include "zgamma/zg_utilities.hpp"

using std::make_shared;
using std::set;
using std::shared_ptr;
using std::string;
using std::vector;
using fastforest::FastForest;
using ZgUtilities::ZgSampleLoader;
using ZgUtilities::XGBoostBDTs;
using ZgUtilities::VBFXGBoostBDTs;
using ZgUtilities::xgb_ggf_bdt_250923_offsets;
using ZgUtilities::xgb_vbf_bdt_250923_offsets;
using ZgUtilities::XGBoostBDTScoreCached;
using ZgUtilities::category_ggf4;
using ZgUtilities::category_ggf3;
using ZgUtilities::category_ggf2;
using ZgUtilities::category_ggf1;
using ZgUtilities::VbfBdts;
using ZgUtilities::vbf_bdt_score;
using ZgUtilities::category_vbf4;
using ZgUtilities::category_vbf3;
using ZgUtilities::category_vbf2;
using ZgUtilities::category_vbf1;
using SelectionList = Datacard::SelectionList;
using Systematic = Datacard::Systematic;
using namespace ZgFunctions;
//const Process::Type data = Process::Type::data;
const Process::Type signal = Process::Type::signal;
const Process::Type background = Process::Type::background;

int main() {

  //setup
  gErrorIgnoreLevel = 6000;

  //Define processes
  vector<shared_ptr<Process>> processes = ZgSampleLoader() 
        .LoadSamples("txt/samples_zgamma.txt","Datacard");
  vector<shared_ptr<Process>> processes_tuneup = ZgSampleLoader() 
        .LoadSamples("txt/samples_zgamma.txt","DatacardTuneUp");
  vector<shared_ptr<Process>> processes_tunedn = ZgSampleLoader() 
        .LoadSamples("txt/samples_zgamma.txt","DatacardTuneDown");
  vector<shared_ptr<Process>> processes_aux = ZgSampleLoader() 
        .LoadSamples("txt/samples_zgamma.txt","DatacardAux"); //m120, m130

  //Define NamedFuncs
  NamedFunc mllg = NamedFunc("llphoton_refit_m").Name("mllg");
  NamedFunc untagged_category_cached = NamedFunc(untagged_category)
      .EnableCaching(true);
  initialize_jetvariations();
  const vector<string> years = {"2016APV", "2016", "2017", "2018", 
                                "2022", "2022EE", "2023", "2023BPix"};

  const vector<FastForest> ggf_bdts = XGBoostBDTs();
  const vector<FastForest> vbf_bdts = VBFXGBoostBDTs();
  const NamedFunc ggf_score_default = XGBoostBDTScoreCached(ggf_bdts, 
      xgb_ggf_bdt_250923_offsets, ggf_bdt_inputs_default, "ggf_default");
  const NamedFunc ggf_score_elscaleup = XGBoostBDTScoreCached(ggf_bdts, 
      xgb_ggf_bdt_250923_offsets, ggf_bdt_inputs_elscaleup, "ggf_elscaleup");
  const NamedFunc ggf_score_elscaledn = XGBoostBDTScoreCached(ggf_bdts, 
      xgb_ggf_bdt_250923_offsets, ggf_bdt_inputs_elscaledn, "ggf_elscaledn");
  const NamedFunc ggf_score_elresup = XGBoostBDTScoreCached(ggf_bdts, 
      xgb_ggf_bdt_250923_offsets, ggf_bdt_inputs_elresup, "ggf_elresup");
  const NamedFunc ggf_score_elresdn = XGBoostBDTScoreCached(ggf_bdts, 
      xgb_ggf_bdt_250923_offsets, ggf_bdt_inputs_elresdn, "ggf_elresdn");
  const NamedFunc ggf_score_muscaleup = XGBoostBDTScoreCached(ggf_bdts, 
      xgb_ggf_bdt_250923_offsets, ggf_bdt_inputs_muscaleup, "ggf_muscaleup");
  const NamedFunc ggf_score_muscaledn = XGBoostBDTScoreCached(ggf_bdts, 
      xgb_ggf_bdt_250923_offsets, ggf_bdt_inputs_muscaledn, "ggf_muscaledn");
  const NamedFunc ggf_score_muresup = XGBoostBDTScoreCached(ggf_bdts, 
      xgb_ggf_bdt_250923_offsets, ggf_bdt_inputs_muresup, "ggf_muresup");
  const NamedFunc ggf_score_muresdn = XGBoostBDTScoreCached(ggf_bdts, 
      xgb_ggf_bdt_250923_offsets, ggf_bdt_inputs_muresdn, "ggf_muresdn");
  const NamedFunc ggf_score_phscaleup = XGBoostBDTScoreCached(ggf_bdts, 
      xgb_ggf_bdt_250923_offsets, ggf_bdt_inputs_phscaleup, "ggf_phscaleup");
  const NamedFunc ggf_score_phscaledn = XGBoostBDTScoreCached(ggf_bdts, 
      xgb_ggf_bdt_250923_offsets, ggf_bdt_inputs_phscaledn, "ggf_phscaledn");
  const NamedFunc ggf_score_phresup = XGBoostBDTScoreCached(ggf_bdts, 
      xgb_ggf_bdt_250923_offsets, ggf_bdt_inputs_phresup, "ggf_phresup");
  const NamedFunc ggf_score_phresdn = XGBoostBDTScoreCached(ggf_bdts, 
      xgb_ggf_bdt_250923_offsets, ggf_bdt_inputs_phresdn, "ggf_phresdn");
  const NamedFunc vbf_score_default = XGBoostBDTScoreCached(vbf_bdts, 
      xgb_vbf_bdt_250923_offsets, vbf_bdt_inputs_default, "vbf_default");
  const NamedFunc vbf_score_elscaleup = XGBoostBDTScoreCached(vbf_bdts, 
      xgb_vbf_bdt_250923_offsets, vbf_bdt_inputs_elscaleup, "vbf_elscaleup");
  const NamedFunc vbf_score_elscaledn = XGBoostBDTScoreCached(vbf_bdts, 
      xgb_vbf_bdt_250923_offsets, vbf_bdt_inputs_elscaledn, "vbf_elscaledn");
  const NamedFunc vbf_score_elresup = XGBoostBDTScoreCached(vbf_bdts, 
      xgb_vbf_bdt_250923_offsets, vbf_bdt_inputs_elresup, "vbf_elresup");
  const NamedFunc vbf_score_elresdn = XGBoostBDTScoreCached(vbf_bdts, 
      xgb_vbf_bdt_250923_offsets, vbf_bdt_inputs_elresdn, "vbf_elresdn");
  const NamedFunc vbf_score_muscaleup = XGBoostBDTScoreCached(vbf_bdts, 
      xgb_vbf_bdt_250923_offsets, vbf_bdt_inputs_muscaleup, "vbf_muscaleup");
  const NamedFunc vbf_score_muscaledn = XGBoostBDTScoreCached(vbf_bdts, 
      xgb_vbf_bdt_250923_offsets, vbf_bdt_inputs_muscaledn, "vbf_muscaledn");
  const NamedFunc vbf_score_muresup = XGBoostBDTScoreCached(vbf_bdts, 
      xgb_vbf_bdt_250923_offsets, vbf_bdt_inputs_muresup, "vbf_muresup");
  const NamedFunc vbf_score_muresdn = XGBoostBDTScoreCached(vbf_bdts, 
      xgb_vbf_bdt_250923_offsets, vbf_bdt_inputs_muresdn, "vbf_muresdn");
  const NamedFunc vbf_score_phscaleup = XGBoostBDTScoreCached(vbf_bdts, 
      xgb_vbf_bdt_250923_offsets, vbf_bdt_inputs_phscaleup, "vbf_phscaleup");
  const NamedFunc vbf_score_phscaledn = XGBoostBDTScoreCached(vbf_bdts, 
      xgb_vbf_bdt_250923_offsets, vbf_bdt_inputs_phscaledn, "vbf_phscaledn");
  const NamedFunc vbf_score_phresup = XGBoostBDTScoreCached(vbf_bdts, 
      xgb_vbf_bdt_250923_offsets, vbf_bdt_inputs_phresup, "vbf_phresup");
  const NamedFunc vbf_score_phresdn = XGBoostBDTScoreCached(vbf_bdts, 
      xgb_vbf_bdt_250923_offsets, vbf_bdt_inputs_phresdn, "vbf_phresdn");

  vector<NamedFunc> ggf_score_jetscaleup;
  vector<NamedFunc> ggf_score_jetscaledn;
  vector<NamedFunc> ggf_score_jetresup;
  vector<NamedFunc> ggf_score_jetresdn;
  vector<NamedFunc> vbf_score_jetscaleup;
  vector<NamedFunc> vbf_score_jetscaledn;
  vector<NamedFunc> vbf_score_jetresup;
  vector<NamedFunc> vbf_score_jetresdn;
  for (unsigned iyear = 0; iyear < years.size(); iyear++) {
    string year = years[iyear];
    ggf_score_jetscaleup.push_back(XGBoostBDTScoreCached(ggf_bdts, 
      xgb_ggf_bdt_250923_offsets, ggf_bdt_inputs_jetscaleup[iyear], 
      "ggf_jetscaleup"+year));
    ggf_score_jetscaledn.push_back(XGBoostBDTScoreCached(ggf_bdts, 
      xgb_ggf_bdt_250923_offsets, ggf_bdt_inputs_jetscaledn[iyear], 
      "ggf_jetscaledn"+year));
    ggf_score_jetresup.push_back(XGBoostBDTScoreCached(ggf_bdts, 
      xgb_ggf_bdt_250923_offsets, ggf_bdt_inputs_jetresup[iyear], 
      "ggf_jetresup"+year));
    ggf_score_jetresdn.push_back(XGBoostBDTScoreCached(ggf_bdts, 
      xgb_ggf_bdt_250923_offsets, ggf_bdt_inputs_jetresdn[iyear], 
      "ggf_jetresdn"+year));
    vbf_score_jetscaleup.push_back(XGBoostBDTScoreCached(vbf_bdts, 
      xgb_vbf_bdt_250923_offsets, vbf_bdt_inputs_jetscaleup[iyear], 
      "vbf_jetscaleup"+year));
    vbf_score_jetscaledn.push_back(XGBoostBDTScoreCached(vbf_bdts, 
      xgb_vbf_bdt_250923_offsets, vbf_bdt_inputs_jetscaledn[iyear], 
      "vbf_jetscaledn"+year));
    vbf_score_jetresup.push_back(XGBoostBDTScoreCached(vbf_bdts, 
      xgb_vbf_bdt_250923_offsets, vbf_bdt_inputs_jetresup[iyear], 
      "vbf_jetresup"+year));
    vbf_score_jetresdn.push_back(XGBoostBDTScoreCached(vbf_bdts, 
      xgb_vbf_bdt_250923_offsets, vbf_bdt_inputs_jetresdn[iyear], 
      "vbf_jetresdn"+year));
  }

  //weight with some regularization
  const NamedFunc weight_reg(
      "weight_reg",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleTypeString().Contains("-")) 
      return 1.; //data
    if (fabs(b.weight()/b.w_lumi()) > 10) return b.w_lumi()*10.0;
    return b.weight();
  });

  //Define weight
  NamedFunc weight(w_years*weight_reg); 

  //Define channels
  SelectionList baseline("baseline");
  baseline.AddSelection("objectreq","nphoton>=1&&nll>=1");
  baseline.AddSelection("lepptcuts",sys_trig_pt_default);
  baseline.AddSelection("zmassreq","ll_m[0]>80&&ll_m[0]<100");
  baseline.AddSelection("photonptreq",
                        "(photon_pt[0]/llphoton_m[0])>=15.0/110.0");
  baseline.AddSelection("mllmllgreq","(llphoton_m[0]+ll_m[0])>=185");
  baseline.AddSelection("metfilters","pass");

  SelectionList cat_ggf4("cat_ggf4",baseline);
  SelectionList cat_ggf3("cat_ggf3",baseline);
  SelectionList cat_ggf2("cat_ggf2",baseline);
  SelectionList cat_ggf1("cat_ggf1",baseline);
  SelectionList cat_vbf4("cat_vbf4",baseline);
  SelectionList cat_vbf3("cat_vbf3",baseline);
  SelectionList cat_vbf2("cat_vbf2",baseline);
  SelectionList cat_vbf1("cat_vbf1",baseline);
  SelectionList cat_vhmet("cat_vhmet",baseline);
  SelectionList cat_vh3l("cat_vh3l",baseline);
  SelectionList cat_tthhad("cat_tthhad",baseline);
  SelectionList cat_tthlep("cat_tthlep",baseline);
  SelectionList cat_untagged("cat_untagged",baseline);

  cat_ggf4.AddSelection("ggfobjectreq","nlep==2&&njet<2&&met<90");
  cat_ggf4.AddSelection("ggf4bdtcuts",category_ggf4(ggf_score_default));
  cat_ggf3.AddSelection("ggfobjectreq","nlep==2&&njet<2&&met<90");
  cat_ggf3.AddSelection("ggf3bdtcuts",category_ggf3(ggf_score_default));
  cat_ggf2.AddSelection("ggfobjectreq","nlep==2&&njet<2&&met<90");
  cat_ggf2.AddSelection("ggf2bdtcuts",category_ggf2(ggf_score_default));
  cat_ggf1.AddSelection("ggfobjectreq","nlep==2&&njet<2&&met<90");
  cat_ggf1.AddSelection("ggf1bdtcuts",category_ggf1(ggf_score_default));

  cat_vbf4.AddSelection("vbfobjectreq","nlep==2&&njet>=2&&nbdfm==0");
  cat_vbf4.AddSelection("vbf4bdtcuts",category_vbf4(vbf_score_default));
  cat_vbf3.AddSelection("vbfobjectreq","nlep==2&&njet>=2&&nbdfm==0");
  cat_vbf3.AddSelection("vbf3bdtcuts",category_vbf3(vbf_score_default));
  cat_vbf2.AddSelection("vbfobjectreq","nlep==2&&njet>=2&&nbdfm==0");
  cat_vbf2.AddSelection("vbf2bdtcuts",category_vbf2(vbf_score_default));
  cat_vbf1.AddSelection("vbfobjectreq","nlep==2&&njet>=2&&nbdfm==0");
  cat_vbf1.AddSelection("vbf1bdtcuts",category_vbf1(vbf_score_default));

  cat_vhmet.AddSelection("vhmetobjectreq","nlep==2&&njet<2&&met>90");
  cat_vhmet.AddSelection("vhmetptllgreq","llphoton_pt[0]/llphoton_m[0]>0.4");
  cat_vh3l.AddSelection("vh3lobjectreq","nlep>=3&&nbdfm==0&&met>30");
  cat_vh3l.AddSelection("vh3lminisoreq",max_lep_miniso<0.15);
  cat_vh3l.AddSelection("vh3lptllgreq","llphoton_pt[0]/llphoton_m[0]>0.3");

  cat_tthhad.AddSelection("tthhadobjectreq","nlep==2&&njet>=5&&nbdfm>=1");
  cat_tthhad.AddSelection("tthhadzmassreq","ll_m[0]>85&&ll_m[0]<95");
  cat_tthlep.AddSelection("tthlepobjectreq",
      "(nlep==3&&njet>=3&&nbdfm>=1)||(nlep>=4&&njet>=1&&nbdfm>=1)");
  cat_tthlep.AddSelection("tthlepminisoreq",max_lep_miniso<0.1);

  cat_untagged.AddSelection("untagged",untagged_category_cached);

  vector<SelectionList> channels = {cat_ggf4,cat_ggf3,cat_ggf2,cat_ggf1,
                                    cat_vbf4,cat_vbf3,cat_vbf2,cat_vbf1,
                                    cat_vh3l,cat_vhmet,cat_tthhad,cat_tthlep,
                                    cat_untagged};

  //Define systematics
  //follow CAT naming when possible
  vector<Systematic> systematics;
  systematics.push_back(Systematic("pdf_Higgs_gg",{"weight"},
                                   {weight*sys_w_pdf_ggf}));
  systematics.push_back(Systematic("pdf_Higgs_qqbar",{"weight"},
                                   {weight*sys_w_pdf_qq}));
  systematics.push_back(Systematic("pdf_Higgs_ttH",{"weight"},
                                   {weight*sys_w_pdf_tth}));
  systematics.push_back(Systematic("QCD_scale_ggH",{"weight"},
                                   {weight*sys_w_ggf_xs}));
  systematics.push_back(Systematic("QCD_scale_qqH",{"weight"},
                                   {weight*sys_w_vbf_xs_up},
                                   {weight*sys_w_vbf_xs_dn}));
  systematics.push_back(Systematic("QCD_scale_VH",{"weight"},
                                   {weight*sys_w_wh_xs_up*sys_w_zh_xs_up},
                                   {weight*sys_w_wh_xs_dn*sys_w_zh_xs_dn}));
  systematics.push_back(Systematic("QCD_scale_ttH",{"weight"},
                                   {weight*sys_w_tth_xs_up},
                                   {weight*sys_w_tth_xs_dn}));
  systematics.push_back(Systematic("QCD_scale_diff_ggH",{"weight"},
                                   {weight*sys_diffscale_ggf_up},
                                   {weight*sys_diffscale_ggf_dn}));
  systematics.push_back(Systematic("QCD_scale_diff_qqH",{"weight"},
                                   {weight*sys_diffscale_vbf_up},
                                   {weight*sys_diffscale_vbf_dn}));
  systematics.push_back(Systematic("QCD_scale_diff_VH",{"weight"},
                                   {weight*sys_diffscale_vh_up},
                                   {weight*sys_diffscale_vh_dn}));
  systematics.push_back(Systematic("QCD_scale_diff_ttH",{"weight"},
                                   {weight*sys_diffscale_tth_up},
                                   {weight*sys_diffscale_tth_dn}));
  systematics.push_back(Systematic("BR_hzg",{"weight"},
                                   {weight*sys_w_htozg_br}));
  systematics.push_back(Systematic("BR_hmm",{"weight"},
                                   {weight*sys_w_htomumu_br}));
  //no CAT guidance on parameter naming
  systematics.push_back(Systematic("param_alphas",{"weight"},
                                   {weight*sys_w_alphas}));
  systematics.push_back(Systematic("param_mq",{"weight"},
                                   {weight*sys_w_mq}));
  systematics.push_back(Systematic("ps_isr",{"weight"},
                                   {weight*"sys_ps[0]"},
                                   {weight*"sys_ps[2]"}));
  systematics.push_back(Systematic("ps_fsr",{"weight"},
                                   {weight*"sys_ps[1]"},
                                   {weight*"sys_ps[3]"}));
  systematics.push_back(Systematic("underlying_event", {}, {}, {}, false, 
                                   processes_tuneup, processes_tunedn));
  //no lumi correlation yet, waiting on official recos from LUM
  systematics.push_back(Systematic("lumi_13TeV",{"weight"},
                                   {weight*sys_w_lumi_run2}));
  systematics.push_back(Systematic("lumi_2022",{"weight"},
                                   {weight*sys_w_lumi_2022}));
  systematics.push_back(Systematic("lumi_2023",{"weight"},
                                   {weight*sys_w_lumi_2023}));
  systematics.push_back(Systematic("CMS_pileup",{"weight"},
                                   {weight*"sys_pu[0]/w_pu"},
                                   {weight*"sys_pu[1]/w_pu"}));
  systematics.push_back(Systematic("CMS_l1_ecal_prefiring",{"weight"},
                                   {weight*"sys_prefire[0]/w_prefire"},
                                   {weight*"sys_prefire[1]/w_prefire"}));
  systematics.push_back(Systematic("CMS_eff_e",{"weight"},
                                   {weight*sys_w_el_up},
                                   {weight*sys_w_el_dn}));
  systematics.push_back(Systematic("CMS_eff_m",{"weight"},
                                   {weight*sys_w_mu_up},
                                   {weight*sys_w_mu_dn}));
  systematics.push_back(Systematic("CMS_eff_g",{"weight"},
                                   {weight*"sys_photon[0]/w_photon"},
                                   {weight*"sys_photon[1]/w_photon"}));
  systematics.push_back(Systematic("CMS_trigger_e",{"weight"},
                                   {weight*sys_trig_el_up},
                                   {weight*sys_trig_el_dn}));
  systematics.push_back(Systematic("CMS_trigger_m",{"weight"},
                                   {weight*sys_trig_mu_up},
                                   {weight*sys_trig_mu_dn}));
  systematics.push_back(Systematic("CMS_quality_g",{"weight"},
                                   {weight*"(1.5*w_phshape-0.5)/w_phshape"},
                                   {weight*"(0.5*w_phshape+0.5)/w_phshape"}));
  systematics.push_back(Systematic("CMS_btag_heavy",{"weight"},
                                   {weight*"sys_bchig[0]/w_bhig_df"},
                                   {weight*"sys_bchig[1]/w_bhig_df"}));
  systematics.push_back(Systematic("CMS_btag_light",{"weight"},
                                   {weight*"sys_udsghig[0]/w_bhig_df"},
                                   {weight*"sys_udsghig[1]/w_bhig_df"}));
  for (unsigned iyear = 0; iyear < years.size(); iyear++) {
    string year = years[iyear];
    systematics.push_back(Systematic("CMS_btag_heavy_"+year,{"weight"},
                                     {weight*sys_bchig_uncorr_up[iyear]},
                                     {weight*sys_bchig_uncorr_dn[iyear]}));
    systematics.push_back(Systematic("CMS_btag_light_"+year,{"weight"},
                                     {weight*sys_udsghig_uncorr_up[iyear]},
                                     {weight*sys_udsghig_uncorr_up[iyear]}));
  }
  systematics.push_back(Systematic("CMS_scale_e",
      {"objectreq","lepptcuts","zmassreq","photonptreq","mllmllgreq",
       "ggfobjectreq","vbfobjectreq","vhmetobjectreq","vh3lobjectreq",
       "tthhadobjectreq","tthlepobjectreq","untagged","vhmetptllgreq",
       "vh3lminisoreq","vh3lptllgreq","tthhadzmassreq","tthlepminisoreq",
       "ggf4bdtcuts","ggf3bdtcuts","ggf2bdtcuts","ggf1bdtcuts",
       "vbf4bdtcuts","vbf3bdtcuts","vbf2bdtcuts","vbf1bdtcuts", "fitvar"},
      {"nphoton>=1"&&sys_nll_elscaleup>=1, 
       sys_trig_pt_elscaleup, 
       sys_ll_m_elscaleup>80&&sys_ll_m_elscaleup<100, 
       "photon_pt[0]"/sys_llphoton_m_elscaleup>(15.0/110.0),
       (sys_ll_m_elscaleup+sys_llphoton_m_elscaleup)>185.0,
       sys_nlep_elscaleup==2&&"njet<2&&met<90",
       sys_nlep_elscaleup==2&&"njet>=2&&nbdfm==0",
       sys_nlep_elscaleup==2&&"njet<2&&met>90",
       sys_nlep_elscaleup>=3&&"nbdfm==0&&met>30",
       sys_nlep_elscaleup==2&&"nbdfm>=1&&njet>=5",
       ((sys_nlep_elscaleup==3&&"njet>=3")||(sys_nlep_elscaleup>=4))
         &&"nbdfm>=1",
       untagged_category_elscaleup,
       sys_llphoton_pt_elscaleup/sys_llphoton_m_elscaleup>0.4,
       sys_max_lep_miniso_elscaleup<0.15,
       sys_llphoton_pt_elscaleup/sys_llphoton_m_elscaleup>0.3,
       sys_ll_m_elscaleup>85&&sys_ll_m_elscaleup<95, 
       sys_max_lep_miniso_elscaleup<0.1,
       category_ggf4(ggf_score_elscaleup),
       category_ggf3(ggf_score_elscaleup),
       category_ggf2(ggf_score_elscaleup),
       category_ggf1(ggf_score_elscaleup),
       category_vbf4(vbf_score_elscaleup),
       category_vbf3(vbf_score_elscaleup),
       category_vbf2(vbf_score_elscaleup),
       category_vbf1(vbf_score_elscaleup),
       sys_llphoton_refit_m_elscaleup},
      {"nphoton>=1"&&sys_nll_elscaledn>=1, 
       sys_trig_pt_elscaledn, 
       sys_ll_m_elscaledn>80&&sys_ll_m_elscaledn<100, 
       "photon_pt[0]"/sys_llphoton_m_elscaledn>(15.0/110.0),
       (sys_ll_m_elscaledn+sys_llphoton_m_elscaledn)>185.0,
       sys_nlep_elscaledn==2&&"njet<2&&met<90",
       sys_nlep_elscaledn==2&&"njet>=2&&nbdfm==0",
       sys_nlep_elscaledn==2&&"njet<2&&met>90",
       sys_nlep_elscaledn>=3&&"nbdfm==0&&met>30",
       sys_nlep_elscaledn==2&&"nbdfm>=1&&njet>=5",
       ((sys_nlep_elscaledn==3&&"njet>=3")||(sys_nlep_elscaledn>=4))
         &&"nbdfm>=1",
       untagged_category_elscaledn,
       sys_llphoton_pt_elscaledn/sys_llphoton_m_elscaledn>0.4,
       sys_max_lep_miniso_elscaledn<0.15,
       sys_llphoton_pt_elscaledn/sys_llphoton_m_elscaledn>0.3,
       sys_ll_m_elscaledn>85&&sys_ll_m_elscaledn<95, 
       sys_max_lep_miniso_elscaledn<0.1,
       category_ggf4(ggf_score_elscaledn),
       category_ggf3(ggf_score_elscaledn),
       category_ggf2(ggf_score_elscaledn),
       category_ggf1(ggf_score_elscaledn),
       category_vbf4(vbf_score_elscaledn),
       category_vbf3(vbf_score_elscaledn),
       category_vbf2(vbf_score_elscaledn),
       category_vbf1(vbf_score_elscaledn),
       sys_llphoton_refit_m_elscaledn},true));
  systematics.push_back(Systematic("CMS_res_e",
      {"objectreq","lepptcuts","zmassreq","photonptreq","mllmllgreq",
       "ggfobjectreq","vbfobjectreq","vhmetobjectreq","vh3lobjectreq",
       "tthhadobjectreq","tthlepobjectreq","untagged","vhmetptllgreq",
       "vh3lminisoreq","vh3lptllgreq","tthhadzmassreq","tthlepminisoreq",
       "ggf4bdtcuts","ggf3bdtcuts","ggf2bdtcuts","ggf1bdtcuts",
       "vbf4bdtcuts","vbf3bdtcuts","vbf2bdtcuts","vbf1bdtcuts", "fitvar"},
      {"nphoton>=1"&&sys_nll_elresup>=1, 
       sys_trig_pt_elresup, 
       sys_ll_m_elresup>80&&sys_ll_m_elresup<100, 
       "photon_pt[0]"/sys_llphoton_m_elresup>(15.0/110.0),
       (sys_ll_m_elresup+sys_llphoton_m_elresup)>185.0,
       sys_nlep_elresup==2&&"njet<2&&met<90",
       sys_nlep_elresup==2&&"njet>=2&&nbdfm==0",
       sys_nlep_elresup==2&&"njet<2&&met>90",
       sys_nlep_elresup>=3&&"nbdfm==0&&met>30",
       sys_nlep_elresup==2&&"nbdfm>=1&&njet>=5",
       ((sys_nlep_elresup==3&&"njet>=3")||(sys_nlep_elresup>=4))&&"nbdfm>=1",
       untagged_category_elresup,
       sys_llphoton_pt_elresup/sys_llphoton_m_elresup>0.4,
       sys_max_lep_miniso_elresup<0.15,
       sys_llphoton_pt_elresup/sys_llphoton_m_elresup>0.3,
       sys_ll_m_elresup>85&&sys_ll_m_elresup<95, 
       sys_max_lep_miniso_elresup<0.1,
       category_ggf4(ggf_score_elresup),
       category_ggf3(ggf_score_elresup),
       category_ggf2(ggf_score_elresup),
       category_ggf1(ggf_score_elresup),
       category_vbf4(vbf_score_elresup),
       category_vbf3(vbf_score_elresup),
       category_vbf2(vbf_score_elresup),
       category_vbf1(vbf_score_elresup),
       sys_llphoton_refit_m_elresup},
      {"nphoton>=1"&&sys_nll_elresdn>=1, 
       sys_trig_pt_elresdn, 
       sys_ll_m_elresdn>80&&sys_ll_m_elresdn<100, 
       "photon_pt[0]"/sys_llphoton_m_elresdn>(15.0/110.0),
       (sys_ll_m_elresdn+sys_llphoton_m_elresdn)>185.0,
       sys_nlep_elresdn==2&&"njet<2&&met<90",
       sys_nlep_elresdn==2&&"njet>=2&&nbdfm==0",
       sys_nlep_elresdn==2&&"njet<2&&met>90",
       sys_nlep_elresdn>=3&&"nbdfm==0&&met>30",
       sys_nlep_elresdn==2&&"nbdfm>=1&&njet>=5",
       ((sys_nlep_elresdn==3&&"njet>=3")||(sys_nlep_elresdn>=4))&&"nbdfm>=1",
       untagged_category_elresdn,
       sys_llphoton_pt_elresdn/sys_llphoton_m_elresdn>0.4,
       sys_max_lep_miniso_elresdn<0.15,
       sys_llphoton_pt_elresdn/sys_llphoton_m_elresdn>0.3,
       sys_ll_m_elresdn>85&&sys_ll_m_elresdn<95, 
       sys_max_lep_miniso_elresdn<0.1,
       category_ggf4(ggf_score_elresdn),
       category_ggf3(ggf_score_elresdn),
       category_ggf2(ggf_score_elresdn),
       category_ggf1(ggf_score_elresdn),
       category_vbf4(vbf_score_elresdn),
       category_vbf3(vbf_score_elresdn),
       category_vbf2(vbf_score_elresdn),
       category_vbf1(vbf_score_elresdn),
       sys_llphoton_refit_m_elresdn},true));
  systematics.push_back(Systematic("CMS_scale_m",
      {"objectreq","lepptcuts","zmassreq","photonptreq","mllmllgreq",
       "ggfobjectreq","vbfobjectreq","vhmetobjectreq","vh3lobjectreq",
       "tthhadobjectreq","tthlepobjectreq","untagged","vhmetptllgreq",
       "vh3lminisoreq","vh3lptllgreq","tthhadzmassreq","tthlepminisoreq",
       "ggf4bdtcuts","ggf3bdtcuts","ggf2bdtcuts","ggf1bdtcuts",
       "vbf4bdtcuts","vbf3bdtcuts","vbf2bdtcuts","vbf1bdtcuts", "fitvar"},
      {"nphoton>=1"&&sys_nll_muscaleup>=1, 
       sys_trig_pt_muscaleup, 
       sys_ll_m_muscaleup>80&&sys_ll_m_muscaleup<100, 
       "photon_pt[0]"/sys_llphoton_m_muscaleup>(15.0/110.0),
       (sys_ll_m_muscaleup+sys_llphoton_m_muscaleup)>185.0,
       sys_nlep_muscaleup==2&&"njet<2&&met<90",
       sys_nlep_muscaleup==2&&"njet>=2&&nbdfm==0",
       sys_nlep_muscaleup==2&&"njet<2&&met>90",
       sys_nlep_muscaleup>=3&&"nbdfm==0&&met>30",
       sys_nlep_muscaleup==2&&"nbdfm>=1&&njet>=5",
       ((sys_nlep_muscaleup==3&&"njet>=3")||(sys_nlep_muscaleup>=4))
           &&"nbdfm>=1",
       untagged_category_muscaleup,
       sys_llphoton_pt_muscaleup/sys_llphoton_m_muscaleup>0.4,
       sys_max_lep_miniso_muscaleup<0.15,
       sys_llphoton_pt_muscaleup/sys_llphoton_m_muscaleup>0.3,
       sys_ll_m_muscaleup>85&&sys_ll_m_muscaleup<95, 
       sys_max_lep_miniso_muscaleup<0.1,
       category_ggf4(ggf_score_muscaleup),
       category_ggf3(ggf_score_muscaleup),
       category_ggf2(ggf_score_muscaleup),
       category_ggf1(ggf_score_muscaleup),
       category_vbf4(vbf_score_muscaleup),
       category_vbf3(vbf_score_muscaleup),
       category_vbf2(vbf_score_muscaleup),
       category_vbf1(vbf_score_muscaleup),
       sys_llphoton_refit_m_muscaleup},
      {"nphoton>=1"&&sys_nll_muscaledn>=1, 
       sys_trig_pt_muscaledn, 
       sys_ll_m_muscaledn>80&&sys_ll_m_muscaledn<100, 
       "photon_pt[0]"/sys_llphoton_m_muscaledn>(15.0/110.0),
       (sys_ll_m_muscaledn+sys_llphoton_m_muscaledn)>185.0,
       sys_nlep_muscaledn==2&&"njet<2&&met<90",
       sys_nlep_muscaledn==2&&"njet>=2&&nbdfm==0",
       sys_nlep_muscaledn==2&&"njet<2&&met>90",
       sys_nlep_muscaledn>=3&&"nbdfm==0&&met>30",
       sys_nlep_muscaledn==2&&"nbdfm>=1&&njet>=5",
       ((sys_nlep_muscaledn==3&&"njet>=3")||(sys_nlep_muscaledn>=4))
           &&"nbdfm>=1",
       untagged_category_muscaledn,
       sys_llphoton_pt_muscaledn/sys_llphoton_m_muscaledn>0.4,
       sys_max_lep_miniso_muscaledn<0.15,
       sys_llphoton_pt_muscaledn/sys_llphoton_m_muscaledn>0.3,
       sys_ll_m_muscaledn>85&&sys_ll_m_muscaledn<95, 
       sys_max_lep_miniso_muscaledn<0.1,
       category_ggf4(ggf_score_muscaledn),
       category_ggf3(ggf_score_muscaledn),
       category_ggf2(ggf_score_muscaledn),
       category_ggf1(ggf_score_muscaledn),
       category_vbf4(vbf_score_muscaledn),
       category_vbf3(vbf_score_muscaledn),
       category_vbf2(vbf_score_muscaledn),
       category_vbf1(vbf_score_muscaledn),
       sys_llphoton_refit_m_muscaledn},true));
  systematics.push_back(Systematic("CMS_res_m",
      {"objectreq","lepptcuts","zmassreq","photonptreq","mllmllgreq",
       "ggfobjectreq","vbfobjectreq","vhmetobjectreq","vh3lobjectreq",
       "tthhadobjectreq","tthlepobjectreq","untagged","vhmetptllgreq",
       "vh3lminisoreq","vh3lptllgreq","tthhadzmassreq","tthlepminisoreq",
       "ggf4bdtcuts","ggf3bdtcuts","ggf2bdtcuts","ggf1bdtcuts",
       "vbf4bdtcuts","vbf3bdtcuts","vbf2bdtcuts","vbf1bdtcuts", "fitvar"},
      {"nphoton>=1"&&sys_nll_muresup>=1, 
       sys_trig_pt_muresup, 
       sys_ll_m_muresup>80&&sys_ll_m_muresup<100, 
       "photon_pt[0]"/sys_llphoton_m_muresup>(15.0/110.0),
       (sys_ll_m_muresup+sys_llphoton_m_muresup)>185.0,
       sys_nlep_muresup==2&&"njet<2&&met<90",
       sys_nlep_muresup==2&&"njet>=2&&nbdfm==0",
       sys_nlep_muresup==2&&"njet<2&&met>90",
       sys_nlep_muresup>=3&&"nbdfm==0&&met>30",
       sys_nlep_muresup==2&&"nbdfm>=1&&njet>=5",
       ((sys_nlep_muresup==3&&"njet>=3")||(sys_nlep_muresup>=4))&&"nbdfm>=1",
       untagged_category_muresup,
       sys_llphoton_pt_muresup/sys_llphoton_m_muresup>0.4,
       sys_max_lep_miniso_muresup<0.15,
       sys_llphoton_pt_muresup/sys_llphoton_m_muresup>0.3,
       sys_ll_m_muresup>85&&sys_ll_m_muresup<95, 
       sys_max_lep_miniso_muresup<0.1,
       category_ggf4(ggf_score_muresup),
       category_ggf3(ggf_score_muresup),
       category_ggf2(ggf_score_muresup),
       category_ggf1(ggf_score_muresup),
       category_vbf4(vbf_score_muresup),
       category_vbf3(vbf_score_muresup),
       category_vbf2(vbf_score_muresup),
       category_vbf1(vbf_score_muresup),
       sys_llphoton_refit_m_muresup},
      {"nphoton>=1"&&sys_nll_muresdn>=1, 
       sys_trig_pt_muresdn, 
       sys_ll_m_muresdn>80&&sys_ll_m_muresdn<100, 
       "photon_pt[0]"/sys_llphoton_m_muresdn>(15.0/110.0),
       (sys_ll_m_muresdn+sys_llphoton_m_muresdn)>185.0,
       sys_nlep_muresdn==2&&"njet<2&&met<90",
       sys_nlep_muresdn==2&&"njet>=2&&nbdfm==0",
       sys_nlep_muresdn==2&&"njet<2&&met>90",
       sys_nlep_muresdn>=3&&"nbdfm==0&&met>30",
       sys_nlep_muresdn==2&&"nbdfm>=1&&njet>=5",
       ((sys_nlep_muresdn==3&&"njet>=3")||(sys_nlep_muresdn>=4))&&"nbdfm>=1",
       untagged_category_muresdn,
       sys_llphoton_pt_muresdn/sys_llphoton_m_muresdn>0.4,
       sys_max_lep_miniso_muresdn<0.15,
       sys_llphoton_pt_muresdn/sys_llphoton_m_muresdn>0.3,
       sys_ll_m_muresdn>85&&sys_ll_m_muresdn<95, 
       sys_max_lep_miniso_muresdn<0.1,
       category_ggf4(ggf_score_muresdn),
       category_ggf3(ggf_score_muresdn),
       category_ggf2(ggf_score_muresdn),
       category_ggf1(ggf_score_muresdn),
       category_vbf4(vbf_score_muresdn),
       category_vbf3(vbf_score_muresdn),
       category_vbf2(vbf_score_muresdn),
       category_vbf1(vbf_score_muresdn),
       sys_llphoton_refit_m_muresdn},true));
  systematics.push_back(Systematic("CMS_scale_g",
      {"objectreq","photonptreq","mllmllgreq","untagged","vhmetptllgreq",
       "vh3lptllgreq","ggf4bdtcuts","ggf3bdtcuts","ggf2bdtcuts","ggf1bdtcuts",
       "vbf4bdtcuts","vbf3bdtcuts","vbf2bdtcuts","vbf1bdtcuts", "fitvar"},
      {sys_nphoton_scaleup>=1&&"nll>=1",
       sys_lead_photon_pt_scaleup/sys_llphoton_m_phscaleup>(15.0/110.0),
       ("ll_m[0]"+sys_llphoton_m_phscaleup)>185.0,
       untagged_category_phscaleup,
       sys_llphoton_pt_phscaleup/sys_llphoton_m_phscaleup>0.4,
       sys_llphoton_pt_phscaleup/sys_llphoton_m_phscaleup>0.3,
       category_ggf4(ggf_score_phscaleup),
       category_ggf3(ggf_score_phscaleup),
       category_ggf2(ggf_score_phscaleup),
       category_ggf1(ggf_score_phscaleup),
       category_vbf4(vbf_score_phscaleup),
       category_vbf3(vbf_score_phscaleup),
       category_vbf2(vbf_score_phscaleup),
       category_vbf1(vbf_score_phscaleup),
       sys_llphoton_refit_m_phscaleup},
      {sys_nphoton_scaledn>=1&&"nll>=1",
       sys_lead_photon_pt_scaledn/sys_llphoton_m_phscaledn>(15.0/110.0),
       ("ll_m[0]"+sys_llphoton_m_phscaledn)>185.0,
       untagged_category_phscaledn,
       sys_llphoton_pt_phscaledn/sys_llphoton_m_phscaledn>0.4,
       sys_llphoton_pt_phscaledn/sys_llphoton_m_phscaledn>0.3,
       category_ggf4(ggf_score_phscaledn),
       category_ggf3(ggf_score_phscaledn),
       category_ggf2(ggf_score_phscaledn),
       category_ggf1(ggf_score_phscaledn),
       category_vbf4(vbf_score_phscaledn),
       category_vbf3(vbf_score_phscaledn),
       category_vbf2(vbf_score_phscaledn),
       category_vbf1(vbf_score_phscaledn),
       sys_llphoton_refit_m_phscaledn},true));
  systematics.push_back(Systematic("CMS_res_g",
      {"objectreq","photonptreq","mllmllgreq","untagged","vhmetptllgreq",
       "vh3lptllgreq","ggf4bdtcuts","ggf3bdtcuts","ggf2bdtcuts","ggf1bdtcuts",
       "vbf4bdtcuts","vbf3bdtcuts","vbf2bdtcuts","vbf1bdtcuts", "fitvar"},
      {sys_nphoton_resup>=1&&"nll>=1",
       sys_lead_photon_pt_resup/sys_llphoton_m_phresup>(15.0/110.0),
       ("ll_m[0]"+sys_llphoton_m_phresup)>185.0,
       untagged_category_phresup,
       sys_llphoton_pt_phresup/sys_llphoton_m_phresup>0.4,
       sys_llphoton_pt_phresup/sys_llphoton_m_phresup>0.3,
       category_ggf4(ggf_score_phresup),
       category_ggf3(ggf_score_phresup),
       category_ggf2(ggf_score_phresup),
       category_ggf1(ggf_score_phresup),
       category_vbf4(vbf_score_phresup),
       category_vbf3(vbf_score_phresup),
       category_vbf2(vbf_score_phresup),
       category_vbf1(vbf_score_phresup),
       sys_llphoton_refit_m_phresup},
      {sys_nphoton_resdn>=1&&"nll>=1",
       sys_lead_photon_pt_resdn/sys_llphoton_m_phresdn>(15.0/110.0),
       ("ll_m[0]"+sys_llphoton_m_phresdn)>185.0,
       untagged_category_phresdn,
       sys_llphoton_pt_phresdn/sys_llphoton_m_phresdn>0.4,
       sys_llphoton_pt_phresdn/sys_llphoton_m_phresdn>0.3,
       category_ggf4(ggf_score_phresdn),
       category_ggf3(ggf_score_phresdn),
       category_ggf2(ggf_score_phresdn),
       category_ggf1(ggf_score_phresdn),
       category_vbf4(vbf_score_phresdn),
       category_vbf3(vbf_score_phresdn),
       category_vbf2(vbf_score_phresdn),
       category_vbf1(vbf_score_phresdn),
       sys_llphoton_refit_m_phresdn},true));
  for (unsigned iyear = 0; iyear < years.size(); iyear++) {
    string year = years[iyear];
    systematics.push_back(Systematic("CMS_scale_j_"+year,
        {"ggfobjectreq","vbfobjectreq","vhmetobjectreq","vh3lobjectreq",
         "tthhadobjectreq","tthlepobjectreq","untagged","ggf4bdtcuts",
         "ggf3bdtcuts","ggf2bdtcuts","ggf1bdtcuts","vbf4bdtcuts","vbf3bdtcuts",
         "vbf2bdtcuts","vbf1bdtcuts"},
        {"nlep==2"&&sys_njet_scaleup[iyear]<2&&sys_met_scaleup[iyear]<90,
         "nlep==2"&&sys_njet_scaleup[iyear]>=2&&sys_nbdfm_scaleup[iyear]==0.0,
         "nlep==2"&&sys_njet_scaleup[iyear]<2&&sys_met_scaleup[iyear]>90,
         "nlep>=3"&&sys_nbdfm_scaleup[iyear]==0.0&&sys_met_scaleup[iyear]>30,
         "nlep==2"&&sys_nbdfm_scaleup[iyear]>=1&&sys_njet_scaleup[iyear]>=5,
         (("nlep==3"&&sys_njet_scaleup[iyear]>=3)||"nlep>=4")
           &&sys_nbdfm_scaleup[iyear]>=1,
         untagged_category_jetscaleup[iyear],
         category_ggf4(ggf_score_jetscaleup[iyear]),
         category_ggf3(ggf_score_jetscaleup[iyear]),
         category_ggf2(ggf_score_jetscaleup[iyear]),
         category_ggf1(ggf_score_jetscaleup[iyear]),
         category_vbf4(vbf_score_jetscaleup[iyear]),
         category_vbf3(vbf_score_jetscaleup[iyear]),
         category_vbf2(vbf_score_jetscaleup[iyear]),
         category_vbf1(vbf_score_jetscaleup[iyear])},
        {"nlep==2"&&sys_njet_scaledn[iyear]<2&&sys_met_scaledn[iyear]<90,
         "nlep==2"&&sys_njet_scaledn[iyear]>=2&&sys_nbdfm_scaledn[iyear]==0.0,
         "nlep==2"&&sys_njet_scaledn[iyear]<2&&sys_met_scaledn[iyear]>90,
         "nlep>=3"&&sys_nbdfm_scaledn[iyear]==0.0&&sys_met_scaledn[iyear]>30,
         "nlep==2"&&sys_nbdfm_scaledn[iyear]>=1&&sys_njet_scaledn[iyear]>=5,
         (("nlep==3"&&sys_njet_scaledn[iyear]>=3)||"nlep>=4")
           &&sys_nbdfm_scaledn[iyear]>=1,
         untagged_category_jetscaledn[iyear],
         category_ggf4(ggf_score_jetscaledn[iyear]),
         category_ggf3(ggf_score_jetscaledn[iyear]),
         category_ggf2(ggf_score_jetscaledn[iyear]),
         category_ggf1(ggf_score_jetscaledn[iyear]),
         category_vbf4(vbf_score_jetscaledn[iyear]),
         category_vbf3(vbf_score_jetscaledn[iyear]),
         category_vbf2(vbf_score_jetscaledn[iyear]),
         category_vbf1(vbf_score_jetscaledn[iyear])},false));
    systematics.push_back(Systematic("CMS_res_j_"+year,
        {"ggfobjectreq","vbfobjectreq","vhmetobjectreq","vh3lobjectreq",
         "tthhadobjectreq","tthlepobjectreq","untagged","ggf4bdtcuts",
         "ggf3bdtcuts","ggf2bdtcuts","ggf1bdtcuts","vbf4bdtcuts","vbf3bdtcuts",
         "vbf2bdtcuts","vbf1bdtcuts"},
        {"nlep==2"&&sys_njet_resup[iyear]<2&&sys_met_resup[iyear]<90,
         "nlep==2"&&sys_njet_resup[iyear]>=2&&sys_nbdfm_resup[iyear]==0.0,
         "nlep==2"&&sys_njet_resup[iyear]<2&&sys_met_resup[iyear]>90,
         "nlep>=3"&&sys_nbdfm_resup[iyear]==0.0&&sys_met_resup[iyear]>30,
         "nlep==2"&&sys_nbdfm_resup[iyear]>=1&&sys_njet_resup[iyear]>=5,
         (("nlep==3"&&sys_njet_resup[iyear]>=3)||"nlep>=4")
           &&sys_nbdfm_resup[iyear]>=1,
         untagged_category_jetresup[iyear],
         category_ggf4(ggf_score_jetresup[iyear]),
         category_ggf3(ggf_score_jetresup[iyear]),
         category_ggf2(ggf_score_jetresup[iyear]),
         category_ggf1(ggf_score_jetresup[iyear]),
         category_vbf4(vbf_score_jetresup[iyear]),
         category_vbf3(vbf_score_jetresup[iyear]),
         category_vbf2(vbf_score_jetresup[iyear]),
         category_vbf1(vbf_score_jetresup[iyear])},
        {"nlep==2"&&sys_njet_resdn[iyear]<2&&sys_met_resdn[iyear]<90,
         "nlep==2"&&sys_njet_resdn[iyear]>=2&&sys_nbdfm_resdn[iyear]==0.0,
         "nlep==2"&&sys_njet_resdn[iyear]<2&&sys_met_resdn[iyear]>90,
         "nlep>=3"&&sys_nbdfm_resdn[iyear]==0.0&&sys_met_resdn[iyear]>30,
         "nlep==2"&&sys_nbdfm_resdn[iyear]>=1&&sys_njet_resdn[iyear]>=5,
         (("nlep==3"&&sys_njet_resdn[iyear]>=3)||"nlep>=4")
           &&sys_nbdfm_resdn[iyear]>=1,
         untagged_category_jetresdn[iyear],
         category_ggf4(ggf_score_jetresdn[iyear]),
         category_ggf3(ggf_score_jetresdn[iyear]),
         category_ggf2(ggf_score_jetresdn[iyear]),
         category_ggf1(ggf_score_jetresdn[iyear]),
         category_vbf4(vbf_score_jetresdn[iyear]),
         category_vbf3(vbf_score_jetresdn[iyear]),
         category_vbf2(vbf_score_jetresdn[iyear]),
         category_vbf1(vbf_score_jetresdn[iyear])},false));
  }


  //Make datacard
  PlotMaker pm;
  pm.multithreaded_ = true;
  pm.min_print_ = true;
  pm.max_threads_ = 16;

  //set axis range to be larger than range in any individual category
  pm.Push<Datacard>("hzg_datacard_v1p0", channels, systematics, 
      processes, weight,
      Axis(340, 95.0, 180.0, mllg, "m_{ll#gamma} [GeV]", {}))
      .AddHistOnlyProcesses(processes_aux)
      .AddParametricProcess("background")
      .SaveDataAsHist()
      .IncludeStatUncertainties();

  pm.max_entries_ = 500;
  pm.MakePlots(1.0);

  return 0;
}
