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
#include "zgamma/syst_functions.hpp"
#include "zgamma/zg_functions.hpp"
#include "zgamma/zg_utilities.hpp"

using std::make_shared;
using std::set;
using std::shared_ptr;
using std::string;
using std::vector;
using fastforest::FastForest;
using ZgUtilities::ZgSampleLoader;
using ZgUtilities::XGBoostBDTs;
using ZgUtilities::category_ggf4;
using ZgUtilities::category_ggf3;
using ZgUtilities::category_ggf2;
using ZgUtilities::category_ggf1;
using ZgUtilities::VbfBdts;
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
        .SetMacro("YEARS",{"2018"})
        .LoadSamples("txt/samples_zgamma.txt","Datacard");

  //TODO implement nullptrs in sample loader
  //TODO implement processes in samples.txt
  vector<shared_ptr<Process>> processes_tuneup = ZgSampleLoader() 
        .SetMacro("YEARS",{"2018"})
        .LoadSamples("txt/samples_zgamma.txt","DatacardTuneUp");

  vector<shared_ptr<Process>> processes_tunedn = ZgSampleLoader() 
        .SetMacro("YEARS",{"2018"})
        .LoadSamples("txt/samples_zgamma.txt","DatacardTuneDown");

  vector<shared_ptr<Process>> processes_aux = ZgSampleLoader() 
        .SetMacro("YEARS",{"2018"})
        .LoadSamples("txt/samples_zgamma.txt","DatacardAux");

  //Define NamedFuncs
  NamedFunc mllg = NamedFunc("llphoton_refit_m").Name("mllg");
  NamedFunc untagged_category_cached = NamedFunc(untagged_category)
      .EnableCaching(true);
  initialize_jetvariations();
  //NamedFunc mllg_range_cut = NamedFunc("llphoton_m[0]>100&&llphoton_m[0]<165");

  const vector<FastForest> kinematic_bdt = XGBoostBDTs();
  NamedFunc ggf4 = category_ggf4(kinematic_bdt);
  NamedFunc ggf3 = category_ggf3(kinematic_bdt);
  NamedFunc ggf2 = category_ggf2(kinematic_bdt);
  NamedFunc ggf1 = category_ggf1(kinematic_bdt);

  vector<shared_ptr<MVAWrapper>> vbf_bdt = VbfBdts();
  NamedFunc vbf4 = category_vbf4(vbf_bdt);
  NamedFunc vbf3 = category_vbf3(vbf_bdt);
  NamedFunc vbf2 = category_vbf2(vbf_bdt);
  NamedFunc vbf1 = category_vbf1(vbf_bdt);

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
  baseline.AddSelection("lepptcuts",trig_plateau_cuts);
  baseline.AddSelection("zmassreq","ll_m[0]>80&&ll_m[0]<100");
  baseline.AddSelection("photonptreq",
                        "(photon_pt[0]/llphoton_m[0])>=15.0/110.0");
  baseline.AddSelection("mllmllgreq","(llphoton_m[0]+ll_m[0])>=185");
  //baseline.AddSelection("fitrange",mllg_range_cut);
  baseline.AddSelection("metfilters","pass");

  SelectionList category_ggf4("cat_ggf4",baseline);
  SelectionList category_ggf3("cat_ggf3",baseline);
  SelectionList category_ggf2("cat_ggf2",baseline);
  SelectionList category_ggf1("cat_ggf1",baseline);
  SelectionList category_vbf4("cat_vbf4",baseline);
  SelectionList category_vbf3("cat_vbf3",baseline);
  SelectionList category_vbf2("cat_vbf2",baseline);
  SelectionList category_vbf1("cat_vbf1",baseline);
  SelectionList category_vhmet("cat_vhmet",baseline);
  SelectionList category_vh3l("cat_vh3l",baseline);
  SelectionList category_tthhad("cat_tthhad",baseline);
  SelectionList category_tthlep("cat_tthlep",baseline);
  SelectionList category_untagged("cat_untagged",baseline);

  category_ggf4.AddSelection("ggfobjectreq","nlep==2&&njet<2&&met<90");
  //category_ggf4.AddSelection("fitrange",
  //                           "llphoton_m[0]>107&&llphoton_m[0]<172");
  category_ggf4.AddSelection("ggf4bdtcuts",ggf4);
  category_ggf3.AddSelection("ggfobjectreq","nlep==2&&njet<2&&met<90");
  //category_ggf3.AddSelection("fitrange",
  //                           "llphoton_m[0]>105&&llphoton_m[0]<170");
  category_ggf3.AddSelection("ggf3bdtcuts",ggf3);
  category_ggf2.AddSelection("ggfobjectreq","nlep==2&&njet<2&&met<90");
  //category_ggf2.AddSelection("fitrange",
  //                           "llphoton_m[0]>103&&llphoton_m[0]<168");
  category_ggf2.AddSelection("ggf2bdtcuts",ggf2);
  category_ggf1.AddSelection("ggfobjectreq","nlep==2&&njet<2&&met<90");
  //category_ggf1.AddSelection("fitrange",
  //                           "llphoton_m[0]>97&&llphoton_m[0]<162");
  category_ggf1.AddSelection("ggf1bdtcuts",ggf1);

  category_vbf4.AddSelection("vbfobjectreq","nlep==2&&njet>=2&&nbdfm==0");
  //category_vbf4.AddSelection("fitrange",
  //                           "llphoton_m[0]>105&&llphoton_m[0]<170");
  category_vbf4.AddSelection("vbf4bdtcuts",vbf4);
  category_vbf3.AddSelection("vbfobjectreq","nlep==2&&njet>=2&&nbdfm==0");
  //category_vbf3.AddSelection("fitrange",
  //                           "llphoton_m[0]>100&&llphoton_m[0]<165");
  category_vbf3.AddSelection("vbf3bdtcuts",vbf3);
  category_vbf2.AddSelection("vbfobjectreq","nlep==2&&njet>=2&&nbdfm==0");
  //category_vbf2.AddSelection("fitrange",
  //                           "llphoton_m[0]>96&&llphoton_m[0]<161");
  category_vbf2.AddSelection("vbf2bdtcuts",vbf2);
  category_vbf1.AddSelection("vbfobjectreq","nlep==2&&njet>=2&&nbdfm==0");
  //category_vbf1.AddSelection("fitrange",
  //                           "llphoton_m[0]>95&&llphoton_m[0]<160");
  category_vbf1.AddSelection("vbfbdtcuts",vbf1);

  category_vhmet.AddSelection("vhmetobjectreq","nlep==2&&njet<2&&met>90");
  category_vhmet.AddSelection("vhmetptllgreq",
                              "llphoton_pt[0]/llphoton_m[0]>0.4");
  //category_vhmet.AddSelection("fitrange",
  //                            "llphoton_m[0]>100&&llphoton_m[0]<165");

  category_vh3l.AddSelection("vh3lobjectreq","nlep>=3&&nbdfm==0&&met>30");
  category_vh3l.AddSelection("vh3lminisoreq",max_lep_miniso<0.15);
  category_vh3l.AddSelection("vh3lptllgreq","llphoton_pt[0]/llphoton_m[0]>0.3");
  //category_vh3l.AddSelection("fitrange",
  //                           "llphoton_m[0]>100&&llphoton_m[0]<165");

  category_tthhad.AddSelection("tthhadobjectreq","nlep==2&&njet>=5&&nbdfm>=1");
  category_tthhad.AddSelection("tthhadzmassreq","ll_m[0]>85&&ll_m[0]<95");
  //category_tthhad.AddSelection("fitrange",
  //                             "llphoton_m[0]>100&&llphoton_m[0]<165");

  category_tthlep.AddSelection("tthlepobjectreq",
      "(nlep==3&&njet>=3&&nbdfm>=1)||(nlep>=4&&njet>=1&&nbdfm>=1)");
  category_tthlep.AddSelection("tthlepminisoreq",max_lep_miniso<0.1);
  //category_tthlep.AddSelection("fitrange",
  //                             "llphoton_m[0]>100&&llphoton_m[0]<165");

  category_untagged.AddSelection("untagged",untagged_category_cached);

  vector<SelectionList> channels = {category_ggf4,category_ggf3,category_ggf2,
                                    category_ggf1,category_vbf4,category_vbf3,
                                    category_vbf2,category_vbf1,category_vh3l,
                                    category_vhmet,category_tthhad,
                                    category_tthlep,category_untagged};

  //Define systematics

  vector<Systematic> systematics;
  //no CAT guidance on alphaS naming
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
                                   {weight*"sys_el[0]/w_el"},
                                   {weight*"sys_el[1]/w_el"}));
  systematics.push_back(Systematic("CMS_eff_m",{"weight"},
                                   {weight*"sys_mu[0]/w_mu"},
                                   {weight*"sys_mu[1]/w_mu"}));
  systematics.push_back(Systematic("CMS_eff_g",{"weight"},
                                   {weight*"sys_photon[0]/w_photon"},
                                   {weight*"sys_photon[1]/w_photon"}));
  //currently includes both e and mu
  systematics.push_back(Systematic("CMS_trigger",{"weight"},
                                   {weight*"sys_trig[0]/w_trig"},
                                   {weight*"sys_trig[1]/w_trig"}));
  //needs to be updated with new production
  systematics.push_back(Systematic("CMS_btag_heavy",{"weight"},
                                   {weight*"sys_bchig[0]/w_bhig_df"},
                                   {weight*"sys_bchig[1]/w_bhig_df"}));
  systematics.push_back(Systematic("CMS_btag_light",{"weight"},
                                   {weight*"sys_udsghig[0]/w_bhig_df"},
                                   {weight*"sys_udsghig[1]/w_bhig_df"}));
  systematics.push_back(Systematic("CMS_scale_e",
      {"objectreq","lepptcuts","zmassreq","photonptreq","mllmllgreq",
       "ggfobjectreq","vbfobjectreq","vhmetobjectreq","vh3lobjectreq",
       "tthhadobjectreq","tthlepobjectreq","untagged","vhmetptllgreq",
       "vh3lminisoreq","vh3lptllgreq","tthhadzmassreq","tthlepminisoreq",
       "fitvar"},
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
       sys_llphoton_refit_m_elscaledn},true));
  systematics.push_back(Systematic("CMS_res_e",
      {"objectreq","lepptcuts","zmassreq","photonptreq","mllmllgreq",
       "ggfobjectreq","vbfobjectreq","vhmetobjectreq","vh3lobjectreq",
       "tthhadobjectreq","tthlepobjectreq","untagged","vhmetptllgreq",
       "vh3lminisoreq","vh3lptllgreq","tthhadzmassreq","tthlepminisoreq",
       "fitvar"},
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
       sys_llphoton_refit_m_elresdn},true));
  systematics.push_back(Systematic("CMS_res_m",
      {"objectreq","lepptcuts","zmassreq","photonptreq","mllmllgreq",
       "ggfobjectreq","vbfobjectreq","vhmetobjectreq","vh3lobjectreq",
       "tthhadobjectreq","tthlepobjectreq","untagged","vhmetptllgreq",
       "vh3lminisoreq","vh3lptllgreq","tthhadzmassreq","tthlepminisoreq",
       "fitvar"},
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
       sys_llphoton_refit_m_muresdn},true));
  systematics.push_back(Systematic("CMS_scale_g",
      {"objectreq","photonptreq","mllmllgreq","untagged","vhmetptllgreq",
       "vh3lptllgreq","fitvar"},
      {sys_nphoton_scaleup>=1&&"nll>=1",
       sys_lead_photon_pt_scaleup/sys_llphoton_m_phscaleup>(15.0/110.0),
       ("ll_m[0]"+sys_llphoton_m_phscaleup)>185.0,
       untagged_category_phscaleup,
       sys_llphoton_pt_phscaleup/sys_llphoton_m_phscaleup>0.4,
       sys_llphoton_pt_phscaleup/sys_llphoton_m_phscaleup>0.3,
       sys_llphoton_refit_m_phscaleup},
      {sys_nphoton_scaledn>=1&&"nll>=1",
       sys_lead_photon_pt_scaledn/sys_llphoton_m_phscaledn>(15.0/110.0),
       ("ll_m[0]"+sys_llphoton_m_phscaledn)>185.0,
       untagged_category_phscaledn,
       sys_llphoton_pt_phscaledn/sys_llphoton_m_phscaledn>0.4,
       sys_llphoton_pt_phscaledn/sys_llphoton_m_phscaledn>0.3,
       sys_llphoton_refit_m_phscaledn},true));
  systematics.push_back(Systematic("CMS_res_g",
      {"objectreq","photonptreq","mllmllgreq","untagged","vhmetptllgreq",
       "vh3lptllgreq","fitvar"},
      {sys_nphoton_resup>=1&&"nll>=1",
       sys_lead_photon_pt_resup/sys_llphoton_m_phresup>(15.0/110.0),
       ("ll_m[0]"+sys_llphoton_m_phresup)>185.0,
       untagged_category_phresup,
       sys_llphoton_pt_phresup/sys_llphoton_m_phresup>0.4,
       sys_llphoton_pt_phresup/sys_llphoton_m_phresup>0.3,
       sys_llphoton_refit_m_phresup},
      {sys_nphoton_resdn>=1&&"nll>=1",
       sys_lead_photon_pt_resdn/sys_llphoton_m_phresdn>(15.0/110.0),
       ("ll_m[0]"+sys_llphoton_m_phresdn)>185.0,
       untagged_category_phresdn,
       sys_llphoton_pt_phresdn/sys_llphoton_m_phresdn>0.4,
       sys_llphoton_pt_phresdn/sys_llphoton_m_phresdn>0.3,
       sys_llphoton_refit_m_phresdn},true));
  const vector<string> years = {"2016APV", "2016", "2017", "2018", 
                                "2022", "2022EE", "2023", "2023BPix"};
  for (unsigned iyear = 0; iyear < years.size(); iyear++) {
    string year = years[iyear];
    systematics.push_back(Systematic("CMS_scale_j_"+year,
        {"ggfobjectreq","vbfobjectreq","vhmetobjectreq","vh3lobjectreq",
         "tthhadobjectreq","tthlepobjectreq","untagged"},
        {"nlep==2"&&sys_njet_scaleup[iyear]<2&&sys_met_scaleup[iyear]<90,
         "nlep==2"&&sys_njet_scaleup[iyear]>=2&&sys_nbdfm_scaleup[iyear]==0.0
          &&"njet>=2", //TODO remove once BDT vars propagated
         "nlep==2"&&sys_njet_scaleup[iyear]<2&&sys_met_scaleup[iyear]>90,
         "nlep>=3"&&sys_nbdfm_scaleup[iyear]==0.0&&sys_met_scaleup[iyear]>30,
         "nlep==2"&&sys_nbdfm_scaleup[iyear]>=1&&sys_njet_scaleup[iyear]>=5,
         (("nlep==3"&&sys_njet_scaleup[iyear]>=3)||"nlep>=4")
           &&sys_nbdfm_scaleup[iyear]>=1,
         untagged_category_jetscaleup[iyear]},
        {"nlep==2"&&sys_njet_scaledn[iyear]<2&&sys_met_scaledn[iyear]<90,
         "nlep==2"&&sys_njet_scaledn[iyear]>=2&&sys_nbdfm_scaledn[iyear]==0.0
          &&"njet>=2", //TODO remove once BDT vars propagated
         "nlep==2"&&sys_njet_scaledn[iyear]<2&&sys_met_scaledn[iyear]>90,
         "nlep>=3"&&sys_nbdfm_scaledn[iyear]==0.0&&sys_met_scaledn[iyear]>30,
         "nlep==2"&&sys_nbdfm_scaledn[iyear]>=1&&sys_njet_scaledn[iyear]>=5,
         (("nlep==3"&&sys_njet_scaledn[iyear]>=3)||"nlep>=4")
           &&sys_nbdfm_scaledn[iyear]>=1,
         untagged_category_jetscaledn[iyear]},false));
    systematics.push_back(Systematic("CMS_res_j_"+year,
        {"ggfobjectreq","vbfobjectreq","vhmetobjectreq","vh3lobjectreq",
         "tthhadobjectreq","tthlepobjectreq","untagged"},
        {"nlep==2"&&sys_njet_resup[iyear]<2&&sys_met_resup[iyear]<90,
         "nlep==2"&&sys_njet_resup[iyear]>=2&&sys_nbdfm_resup[iyear]==0.0
          &&"njet>=2", //TODO remove once BDT vars propagated
         "nlep==2"&&sys_njet_resup[iyear]<2&&sys_met_resup[iyear]>90,
         "nlep>=3"&&sys_nbdfm_resup[iyear]==0.0&&sys_met_resup[iyear]>30,
         "nlep==2"&&sys_nbdfm_resup[iyear]>=1&&sys_njet_resup[iyear]>=5,
         (("nlep==3"&&sys_njet_resup[iyear]>=3)||"nlep>=4")
           &&sys_nbdfm_resup[iyear]>=1,
         untagged_category_jetresup[iyear]},
        {"nlep==2"&&sys_njet_resdn[iyear]<2&&sys_met_resdn[iyear]<90,
         "nlep==2"&&sys_njet_resdn[iyear]>=2&&sys_nbdfm_resdn[iyear]==0.0
          &&"njet>=2", //TODO remove once BDT vars propagated
         "nlep==2"&&sys_njet_resdn[iyear]<2&&sys_met_resdn[iyear]>90,
         "nlep>=3"&&sys_nbdfm_resdn[iyear]==0.0&&sys_met_resdn[iyear]>30,
         "nlep==2"&&sys_nbdfm_resdn[iyear]>=1&&sys_njet_resdn[iyear]>=5,
         (("nlep==3"&&sys_njet_resdn[iyear]>=3)||"nlep>=4")
           &&sys_nbdfm_resdn[iyear]>=1,
         untagged_category_jetresdn[iyear]},false));
  }


  //Make datacard
  PlotMaker pm;
  pm.multithreaded_ = false;
  pm.min_print_ = true;

  //set axis range to be larger than range in any individual category
  pm.Push<Datacard>("test_datacard9", channels, systematics, 
      processes, weight,
      Axis(80, 95.0, 180.0, mllg, "m_{ll#gamma} [GeV]", {}))
      .AddHistOnlyProcesses(processes_aux)
      .AddParametricProcess("background")
      .SaveDataAsHist()
      .IncludeStatUncertainties();

  //pm.max_entries_ = 1000;
  pm.MakePlots(1.0);

  return 0;
}
