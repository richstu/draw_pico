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
using ZgFunctions::max_lep_miniso;
using ZgFunctions::trig_plateau_cuts;
using ZgFunctions::w_years;
using ZgFunctions::sys_w_alphas;
using ZgFunctions::sys_w_pdf_ggf;
using ZgFunctions::sys_w_pdf_qq;
using ZgFunctions::sys_w_pdf_tth;
using ZgFunctions::sys_w_ggf_xs;
using ZgFunctions::sys_w_vbf_xs_up;
using ZgFunctions::sys_w_vbf_xs_dn;
using ZgFunctions::sys_w_wh_xs_up;
using ZgFunctions::sys_w_wh_xs_dn;
using ZgFunctions::sys_w_zh_xs_up;
using ZgFunctions::sys_w_zh_xs_dn;
using ZgFunctions::sys_w_tth_xs_up;
using ZgFunctions::sys_w_tth_xs_dn;
using ZgFunctions::sys_w_htozg_br;
using ZgFunctions::sys_w_htomumu_br;
using ZgFunctions::sys_w_lumi_run2;
using ZgFunctions::sys_w_lumi_2022;
using ZgFunctions::sys_w_lumi_2023;
using ZgFunctions::sys_w_mq;
using ZgFunctions::sys_llphoton_m_elscaleup;
using ZgFunctions::sys_llphoton_m_elscaledn;
using ZgFunctions::sys_llphoton_m_elresup;
using ZgFunctions::sys_llphoton_m_elresdn;
using ZgFunctions::sys_llphoton_m_phscaleup;
using ZgFunctions::sys_llphoton_m_phscaledn;
using ZgFunctions::sys_llphoton_m_phresup;
using ZgFunctions::sys_llphoton_m_phresdn;
using ZgFunctions::sys_llphoton_m_muup;
using ZgFunctions::sys_llphoton_m_mudn;
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
using ZgUtilities::AddGaussStepExponential;
using ZgUtilities::AddGaussStepBernstein;
using SelectionList = Datacard::SelectionList;
using Systematic = Datacard::Systematic;
//const Process::Type data = Process::Type::data;
const Process::Type signal = Process::Type::signal;
const Process::Type background = Process::Type::background;

int main() {

  //setup
  gErrorIgnoreLevel = 6000;

  //Define processes
  vector<shared_ptr<Process>> processes = ZgSampleLoader() 
        //.SetMacro("YEARS",{"2018"})
        .LoadSamples("txt/samples_zgamma.txt","Datacard");

  vector<shared_ptr<Process>> processes_aux = ZgSampleLoader() 
        //.SetMacro("YEARS",{"2018"})
        .LoadSamples("txt/samples_zgamma.txt","DatacardAux");

  //Define NamedFuncs
  NamedFunc mllg = NamedFunc("llphoton_m[0]").Name("mllg");
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

  category_ggf4.AddSelection("catobjectreq","nlep==2&&njet<2&&met<90");
  category_ggf4.AddSelection("fitrange",
                             "llphoton_m[0]>107&&llphoton_m[0]<172");
  category_ggf4.AddSelection("bdtcuts",ggf4);
  category_ggf3.AddSelection("catobjectreq","nlep==2&&njet<2&&met<90");
  category_ggf3.AddSelection("fitrange",
                             "llphoton_m[0]>105&&llphoton_m[0]<170");
  category_ggf3.AddSelection("bdtcuts",ggf3);
  category_ggf2.AddSelection("catobjectreq","nlep==2&&njet<2&&met<90");
  category_ggf2.AddSelection("fitrange",
                             "llphoton_m[0]>103&&llphoton_m[0]<168");
  category_ggf2.AddSelection("bdtcuts",ggf2);
  category_ggf1.AddSelection("catobjectreq","nlep==2&&njet<2&&met<90");
  category_ggf1.AddSelection("fitrange",
                             "llphoton_m[0]>97&&llphoton_m[0]<162");
  category_ggf1.AddSelection("bdtcuts",ggf1);

  category_vbf4.AddSelection("catobjectreq","nlep==2&&njet>=2&&nbdfm==0");
  category_vbf4.AddSelection("fitrange",
                             "llphoton_m[0]>105&&llphoton_m[0]<170");
  category_vbf4.AddSelection("bdtcuts",vbf4);
  category_vbf3.AddSelection("catobjectreq","nlep==2&&njet>=2&&nbdfm==0");
  category_vbf3.AddSelection("fitrange",
                             "llphoton_m[0]>100&&llphoton_m[0]<165");
  category_vbf3.AddSelection("bdtcuts",vbf3);
  category_vbf2.AddSelection("catobjectreq","nlep==2&&njet>=2&&nbdfm==0");
  category_vbf2.AddSelection("fitrange","llphoton_m[0]>96&&llphoton_m[0]<161");
  category_vbf2.AddSelection("bdtcuts",vbf2);
  category_vbf1.AddSelection("catobjectreq","nlep==2&&njet>=2&&nbdfm==0");
  category_vbf1.AddSelection("fitrange","llphoton_m[0]>95&&llphoton_m[0]<160");
  category_vbf1.AddSelection("bdtcuts",vbf1);

  category_vhmet.AddSelection("catobjectreq","nlep==2&&njet<2&&met>90");
  category_vhmet.AddSelection("zmassreq","ll_m[0]>85&&ll_m[0]<95");
  category_vhmet.AddSelection("ptllgreq","llphoton_pt[0]/llphoton_m[0]>0.4");
  category_vhmet.AddSelection("fitrange",
                              "llphoton_m[0]>100&&llphoton_m[0]<165");

  category_vh3l.AddSelection("catobjectreq",
                             "nphoton>=1&&nlep>=3&&nbdfm==0&&met>30");
  category_vh3l.AddSelection("zmassreq","ll_m[0]>85&&ll_m[0]<95");
  category_vh3l.AddSelection("minisoreq",max_lep_miniso<0.15);
  category_vh3l.AddSelection("ptllgreq","llphoton_pt[0]/llphoton_m[0]>0.3");
  category_vh3l.AddSelection("fitrange",
                             "llphoton_m[0]>100&&llphoton_m[0]<165");

  category_tthhad.AddSelection("catobjectreq","nlep==2&&njet>=5&&nbdfm>=1");
  category_tthhad.AddSelection("zmassreq","ll_m[0]>85&&ll_m[0]<95");
  category_tthhad.AddSelection("fitrange",
                               "llphoton_m[0]>100&&llphoton_m[0]<165");

  category_tthlep.AddSelection("catobjectreq",
      "(nlep==3&&njet>=3&&nbdfm>=1)||(nlep>=4&&njet>=1&&nbdfm>=1)");
  category_tthlep.AddSelection("zmassreq","ll_m[0]>85&&ll_m[0]<95");
  category_tthlep.AddSelection("minisoreq",max_lep_miniso<0.1);
  category_tthlep.AddSelection("fitrange",
                               "llphoton_m[0]>100&&llphoton_m[0]<165");

  vector<SelectionList> channels = {category_ggf4,category_ggf3,category_ggf2,
                                    category_ggf1,category_vbf4,category_vbf3,
                                    category_vbf2,category_vbf1,category_vh3l,
                                    category_vhmet,category_tthhad,
                                    category_tthlep};

  //Define systematics

  vector<Systematic> systematics;
  //no CAT guidance on alphaS naming
  systematics.push_back(Systematic("alphas",{"weight"},{weight*sys_w_alphas}));
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
  //no CAT guidance on mq naming
  systematics.push_back(Systematic("mq",{"weight"},
                                   {weight*sys_w_mq}));
  systematics.push_back(Systematic("CMS_scale_e",{"fitvar"},
                                   {sys_llphoton_m_elscaleup},
                                   {sys_llphoton_m_elscaledn},true));
  systematics.push_back(Systematic("CMS_res_e",{"fitvar"},
                                   {sys_llphoton_m_elresup},
                                   {sys_llphoton_m_elresdn},true));
  systematics.push_back(Systematic("CMS_scale_g",{"fitvar"},
                                   {sys_llphoton_m_phscaleup},
                                   {sys_llphoton_m_phscaledn},true));
  systematics.push_back(Systematic("CMS_res_g",{"fitvar"},
                                   {sys_llphoton_m_phresup},
                                   {sys_llphoton_m_phresdn},true));
  systematics.push_back(Systematic("CMS_scale_m",{"fitvar"},
                                   {sys_llphoton_m_muup},
                                   {sys_llphoton_m_mudn},true));

  //Make datacard
  PlotMaker pm;
  pm.multithreaded_ = false;
  pm.min_print_ = true;

  //set axis range to be larger than range in any individual category
  pm.Push<Datacard>("test_datacard7", channels, systematics, 
      processes, weight,
      Axis(80, 95.0, 180.0, mllg, "m_{ll#gamma} [GeV]", {}))
      .AddHistOnlyProcesses(processes_aux)
      .AddParametricProcess("background");

  pm.MakePlots(1.0);

  return 0;
}
