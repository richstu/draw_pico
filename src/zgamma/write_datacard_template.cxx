/**
 * Script to write datacard for H->Zgamma analysis
 */

#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "RooArgList.h"
#include "RooGenericPdf.h"
#include "RooRealVar.h"
#include "TError.h"
#include "TColor.h"

#include "core/axis.hpp"
#include "core/baby_pico.hpp"
#include "core/datacard.hpp"
#include "core/mva_wrapper.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/process.hpp"
#include "core/utilities.hpp"
#include "zgamma/zg_functions.hpp"

using std::set;
using std::shared_ptr;
using std::string;
using std::vector;
using ZgFunctions::HLT_pass_dilepton;
using ZgFunctions::HLT_pass_singlelepton;
using ZgFunctions::lead_lepton_eta;
using ZgFunctions::sublead_lepton_eta;
using ZgFunctions::photon_drmax;
using ZgFunctions::stitch;
using ZgFunctions::w_years;
using SelectionList = Datacard::SelectionList;
using Systematic = Datacard::Systematic;
//const Process::Type data = Process::Type::data;
const Process::Type signal = Process::Type::signal;
const Process::Type background = Process::Type::background;

int main() {

  //setup
  gErrorIgnoreLevel = 6000;

  //Define processes
  string prod_folder("/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/");
  set<string> years = {"2016APV","2016","2017","2018"};
  NamedFunc trig_and_stitch = (HLT_pass_dilepton||HLT_pass_singlelepton)&&stitch;
  shared_ptr<Process> proc_pseudodata = Process::MakeShared<Baby_pico>(
      "data_obs", Process::Type::data, kBlack, attach_folder(prod_folder,years,
      "mc/merged_zgmc_llg",{"*DYJets*","*ZGToLLG*Tune*","*HToZG*M-125*"}),
      trig_and_stitch);
  shared_ptr<Process> proc_background = Process::MakeShared<Baby_pico>(
      "background", background, kBlack, attach_folder(prod_folder,years,
      "mc/merged_zgmc_llg",{"*DYJets*","*ZGToLLG*Tune*"}),
      trig_and_stitch);
  shared_ptr<Process> proc_signal = Process::MakeShared<Baby_pico>(
      "htozg", signal, kRed, attach_folder(prod_folder,years,
      "mc/merged_zgmc_llg",{"*HToZG*M-125*"}), trig_and_stitch);
  vector<shared_ptr<Process>> processes = {proc_pseudodata, proc_signal, proc_background};

  //Define NamedFuncs
  NamedFunc mllg = NamedFunc("llphoton_m[0]").Name("mllg");
  NamedFunc zg_el_cuts("(ll_lepid[0]==11) && (el_pt[ll_i1[0]]>25) && (el_pt[ll_i2[0]]>15)");
  NamedFunc zg_mu_cuts("(ll_lepid[0]==13) && (mu_pt[ll_i1[0]]>20) && (mu_pt[ll_i2[0]]>10)");
  MVAWrapper kin_bdt_reader("kinematic_bdt");
  kin_bdt_reader.SetVariable("photon_mva","photon_idmva[0]");
  kin_bdt_reader.SetVariable("min_dR","photon_drmin[0]");
  kin_bdt_reader.SetVariable("max_dR",photon_drmax);
  kin_bdt_reader.SetVariable("pt_mass","llphoton_pt[0]/llphoton_m[0]");
  kin_bdt_reader.SetVariable("cosTheta","llphoton_cosTheta[0]");
  kin_bdt_reader.SetVariable("costheta","llphoton_costheta[0]");
  kin_bdt_reader.SetVariable("phi","llphoton_psi[0]");
  kin_bdt_reader.SetVariable("photon_res","photon_pterr[0]/photon_pt[0]");
  kin_bdt_reader.SetVariable("photon_rapidity","photon_eta[0]");
  kin_bdt_reader.SetVariable("l1_rapidity",lead_lepton_eta);
  kin_bdt_reader.SetVariable("l2_rapidity",sublead_lepton_eta);
  kin_bdt_reader.BookMVA("/homes/oshiro/analysis/small_phys_utils/dataset/weights/shuffled_kinbdt_masscut_run2_BDT.weights.xml");
  NamedFunc bdt_score = kin_bdt_reader.GetDiscriminant();

  //Define weight
  NamedFunc weight(w_years*"w_lumi"); //no ewkzgamma sample for now

  //Define channels
  SelectionList category_ggh1("cat_ggh1");
  category_ggh1.AddSelection("objectreq","nphoton>=1&&nlep>=2");
  category_ggh1.AddSelection("zmassreq","ll_m[0]>81&&ll_m[0]<101");
  category_ggh1.AddSelection("photonptreq","(photon_pt[0]/llphoton_m[0])>=15.0/110.0");
  category_ggh1.AddSelection("mllmllgreq","(llphoton_m[0]+ll_m[0])>=185");
  category_ggh1.AddSelection("photonidcut","photon_id80[0]");
  category_ggh1.AddSelection("lepptcuts",zg_el_cuts||zg_mu_cuts);
  category_ggh1.AddSelection("bdt_selection",bdt_score>0.52);
  SelectionList category_ggh2("cat_ggh2");
  category_ggh2.AddSelection("objectreq","nphoton>=1&&nlep>=2");
  category_ggh2.AddSelection("zmassreq","ll_m[0]>81&&ll_m[0]<101");
  category_ggh2.AddSelection("photonptreq","(photon_pt[0]/llphoton_m[0])>=15.0/110.0");
  category_ggh2.AddSelection("mllmllgreq","(llphoton_m[0]+ll_m[0])>=185");
  category_ggh2.AddSelection("photonidcut","photon_id80[0]");
  category_ggh2.AddSelection("lepptcuts",zg_el_cuts||zg_mu_cuts);
  category_ggh2.AddSelection("bdt_selection",bdt_score>0.28&&bdt_score<0.52);
  SelectionList category_ggh3("cat_ggh3");
  category_ggh3.AddSelection("objectreq","nphoton>=1&&nlep>=2");
  category_ggh3.AddSelection("zmassreq","ll_m[0]>81&&ll_m[0]<101");
  category_ggh3.AddSelection("photonptreq","(photon_pt[0]/llphoton_m[0])>=15.0/110.0");
  category_ggh3.AddSelection("mllmllgreq","(llphoton_m[0]+ll_m[0])>=185");
  category_ggh3.AddSelection("photonidcut","photon_id80[0]");
  category_ggh3.AddSelection("lepptcuts",zg_el_cuts||zg_mu_cuts);
  category_ggh3.AddSelection("bdt_selection",bdt_score>0.02&&bdt_score<0.28);
  SelectionList category_ggh4("cat_ggh4");
  category_ggh4.AddSelection("objectreq","nphoton>=1&&nlep>=2");
  category_ggh4.AddSelection("zmassreq","ll_m[0]>81&&ll_m[0]<101");
  category_ggh4.AddSelection("photonptreq","(photon_pt[0]/llphoton_m[0])>=15.0/110.0");
  category_ggh4.AddSelection("mllmllgreq","(llphoton_m[0]+ll_m[0])>=185");
  category_ggh4.AddSelection("photonidcut","photon_id80[0]");
  category_ggh4.AddSelection("lepptcuts",zg_el_cuts||zg_mu_cuts);
  category_ggh4.AddSelection("bdt_selection",bdt_score<0.02);
  vector<SelectionList> channels = {category_ggh1, category_ggh2, category_ggh3, category_ggh4};

  //Define systematics
  //Systematic syst_altweight("altw",weight*"w_photon");
  //Systematic syst_altselection("altsel","zmassreq","ll_m[0]>80&&ll_m[0]<101"); //can make namedfunc that selects signal only using type
  //vector<Systematic> systematics = {syst_altweight, syst_altselection};
  vector<Systematic> systematics;

  //Define parametric PDFs
  //RooRealVar rrv_mllg("mllg","llgamma invariant mass",114.0,150.0);
  //RooRealVar c0("c0","exponential coefficient",0.0,10.0);
  //c0.setVal(1.0);
  //RooGenericPdf exp_pdf_el("pdf_background_cat_el","exp_pdf","exp(-1.0*c0*(mllg-114.0)/36.0)",RooArgSet(rrv_mllg,c0));
  //RooGenericPdf exp_pdf_mu("pdf_background_cat_mu","exp_pdf","exp(-1.0*c0*(mllg-114.0)/36.0)",RooArgSet(rrv_mllg,c0));
  //vector<RooAbsPdf*> background_pdfs;
  //background_pdfs.push_back(&exp_pdf_el);
  //background_pdfs.push_back(&exp_pdf_mu);

  //Make datacard
  PlotMaker pm;
  pm.multithreaded_ = false; //use MVA
  pm.min_print_ = true;

  pm.Push<Datacard>("signif_datacard", channels, systematics, processes, weight,
      Axis(60, 100.0, 160.0, mllg, "m_{ll#gamma} [GeV]", {}));

  pm.MakePlots(1.0);

  return 0;
}
