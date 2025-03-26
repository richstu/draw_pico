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
using ZgFunctions::sys_w_mq;
using ZgUtilities::ZgSampleLoader;
using ZgUtilities::KinematicBdt;
using ZgUtilities::category_ggh4;
using ZgUtilities::category_ggh3;
using ZgUtilities::category_ggh2;
using ZgUtilities::category_ggh1;
using ZgUtilities::VbfBdt;
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
        .SetMacro("YEARS",{"2018"})
        .LoadSamples("txt/samples_zgamma.txt","Datacard");

  vector<shared_ptr<Process>> processes_aux = ZgSampleLoader() 
        .SetMacro("YEARS",{"2018"})
        .LoadSamples("txt/samples_zgamma.txt","DatacardAux");

  //Define NamedFuncs
  NamedFunc mllg = NamedFunc("llphoton_m[0]").Name("mllg");
  NamedFunc mllg_range_cut = NamedFunc("llphoton_m[0]>100&&llphoton_m[0]<180");

  shared_ptr<MVAWrapper> kinematic_bdt = KinematicBdt();
  NamedFunc ggh4 = category_ggh4(kinematic_bdt);
  NamedFunc ggh3 = category_ggh3(kinematic_bdt);
  NamedFunc ggh2 = category_ggh2(kinematic_bdt);
  NamedFunc ggh1 = category_ggh1(kinematic_bdt);

  shared_ptr<MVAWrapper> vbf_bdt = VbfBdt();
  NamedFunc vbf3 = category_vbf3(vbf_bdt);
  NamedFunc vbf2 = category_vbf2(vbf_bdt);
  NamedFunc vbf1 = category_vbf1(vbf_bdt);

  //weight with some regularization
  const NamedFunc weight_reg(
      "weight_reg",[](const Baby &b) -> NamedFunc::ScalarType{
    if (fabs(b.weight()/b.w_lumi()) > 10) return b.w_lumi()*10.0;
    return b.weight();
  });

  //Define weight
  NamedFunc weight(w_years*weight_reg); 

  //Define channels
  SelectionList baseline("baseline");
  baseline.AddSelection("objectreq","nphoton>=1&&nlep>=2");
  baseline.AddSelection("lepptcuts",trig_plateau_cuts);
  baseline.AddSelection("zmassreq","ll_m[0]>80&&ll_m[0]<100");
  baseline.AddSelection("photonptreq","(photon_pt[0]/llphoton_m[0])>=15.0/110.0");
  baseline.AddSelection("mllmllgreq","(llphoton_m[0]+ll_m[0])>=185");
  baseline.AddSelection("photonidreq","photon_id80[0]");
  baseline.AddSelection("mllgrange",mllg_range_cut);

  SelectionList category_ggh4("cat_ggh4",baseline);
  SelectionList category_ggh3("cat_ggh3",baseline);
  SelectionList category_ggh2("cat_ggh2",baseline);
  SelectionList category_ggh1("cat_ggh1",baseline);
  SelectionList category_vbf3("cat_vbf3",baseline);
  SelectionList category_vbf2("cat_vbf2",baseline);
  SelectionList category_vbf1("cat_vbf1",baseline);
  SelectionList category_vhmet("cat_vhmet",baseline);
  SelectionList category_vh3l("cat_vh3l",baseline);
  SelectionList category_tthhad("cat_tthhad",baseline);
  SelectionList category_tthlep("cat_tthlep",baseline);

  category_ggh4.AddSelection("catobjectreq","nlep==2&&njet<2&&met<90");
  category_ggh4.AddSelection("bdtcuts",ggh4);
  category_ggh3.AddSelection("catobjectreq","nlep==2&&njet<2&&met<90");
  category_ggh3.AddSelection("bdtcuts",ggh3);
  category_ggh2.AddSelection("catobjectreq","nlep==2&&njet<2&&met<90");
  category_ggh2.AddSelection("bdtcuts",ggh2);
  category_ggh1.AddSelection("catobjectreq","nlep==2&&njet<2&&met<90");
  category_ggh1.AddSelection("bdtcuts",ggh1);

  category_vbf3.AddSelection("catobjectreq","nlep==2&&njet>=2&&nbdfm==0");
  category_vbf3.AddSelection("bdtcuts",vbf3);
  category_vbf2.AddSelection("catobjectreq","nlep==2&&njet>=2&&nbdfm==0");
  category_vbf2.AddSelection("bdtcuts",vbf2);
  category_vbf1.AddSelection("catobjectreq","nlep==2&&njet>=2&&nbdfm==0");
  category_vbf1.AddSelection("bdtcuts",vbf1);

  category_vhmet.AddSelection("catobjectreq","nlep==2&&njet<2&&met>90");
  category_vhmet.AddSelection("zmassreq","ll_m[0]>85&&ll_m[0]<95");
  category_vhmet.AddSelection("ptllgreq","llphoton_pt[0]/llphoton_m[0]>0.4");

  category_vh3l.AddSelection("catobjectreq",
                             "nphoton>=1&&nlep>=3&&nbdfm==0&&met>30");
  category_vh3l.AddSelection("zmassreq","ll_m[0]>85&&ll_m[0]<95");
  category_vh3l.AddSelection("minisoreq",max_lep_miniso<0.15);
  category_vh3l.AddSelection("ptllgreq","llphoton_pt[0]/llphoton_m[0]>0.3");

  category_tthhad.AddSelection("catobjectreq","nlep==2&&njet>=5&&nbdfm>=1");
  category_tthhad.AddSelection("zmassreq","ll_m[0]>85&&ll_m[0]<95");

  category_tthlep.AddSelection("catobjectreq",
      "(nlep==3&&njet>=3&&nbdfm>=1)||(nlep>=4&&njet>=1&&nbdfm>=1)");
  category_tthlep.AddSelection("zmassreq","ll_m[0]>85&&ll_m[0]<95");
  category_tthlep.AddSelection("minisoreq",max_lep_miniso<0.1);

  vector<SelectionList> channels = {category_ggh4,category_ggh3,category_ggh2,
                                    category_ggh1,category_vbf3,category_vbf2,
                                    category_vbf1,category_vh3l,category_vhmet,
                                    category_tthhad,category_tthlep};

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
  //no CAT guidance on mq naming
  systematics.push_back(Systematic("mq",{"weight"},
                                   {weight*sys_w_mq}));

  //Make datacard
  PlotMaker pm;
  pm.multithreaded_ = false;
  pm.min_print_ = true;

  pm.Push<Datacard>("test_datacard2", channels, systematics, processes, weight,
      Axis(80, 100.0, 180.0, mllg, "m_{ll#gamma} [GeV]", {}))
      .AddHistOnlyProcesses(processes_aux)
      .AddParametricProcess("background");

  pm.MakePlots(1.0);

  return 0;
}
