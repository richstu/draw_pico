/**
 * Script to write datacard for H->Zgamma analysis
 */

#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "RooArgList.h"
#include "RooCategory.h"
#include "RooGenericPdf.h"
#include "core/RooMultiPdf.hpp"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "TError.h"
#include "TColor.h"

#include "core/axis.hpp"
#include "core/baby_pico.hpp"
#include "core/datacard.hpp"
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
  set<string> years = {"2018"};
  NamedFunc trig_and_stitch = (HLT_pass_dilepton||HLT_pass_singlelepton)&&stitch;
  shared_ptr<Process> proc_pseudodata = Process::MakeShared<Baby_pico>(
      "data_obs", Process::Type::data, kBlack, attach_folder(prod_folder,years,
      "mc/merged_zgmc_llg",{"*DYJets*","*ZGToLLG*Tune*","*HToZG*M-125*"}),
      trig_and_stitch);
  shared_ptr<Process> proc_signal = Process::MakeShared<Baby_pico>(
      "htozg", signal, kRed, attach_folder(prod_folder,years,
      "mc/merged_zgmc_llg",{"*HToZG*M-125*"}), trig_and_stitch);
  vector<shared_ptr<Process>> processes = {proc_pseudodata, proc_signal};

  //Define NamedFuncs
  NamedFunc mllg = NamedFunc("llphoton_m[0]").Name("mllg");
  NamedFunc zg_el_cuts("(ll_lepid[0]==11) && (el_pt[ll_i1[0]]>25) && (el_pt[ll_i2[0]]>15)");
  NamedFunc zg_mu_cuts("(ll_lepid[0]==13) && (mu_pt[ll_i1[0]]>20) && (mu_pt[ll_i2[0]]>10)");

  //Define weight
  NamedFunc weight(w_years*"w_lumi"); //no ewkzgamma sample for now

  //Define channels
  SelectionList category_electron("cat_el");
  category_electron.AddSelection("objectreq","nphoton>=1&&nlep>=2");
  category_electron.AddSelection("zmassreq","ll_m[0]>81&&ll_m[0]<101");
  category_electron.AddSelection("photonptreq","(photon_pt[0]/llphoton_m[0])>=15.0/110.0");
  category_electron.AddSelection("mllmllgreq","(llphoton_m[0]+ll_m[0])>=185");
  category_electron.AddSelection("photonidcut","photon_id80[0]");
  category_electron.AddSelection("lepptcuts",zg_el_cuts);
  category_electron.AddSelection("mllgcuts","llphoton_m[0]>100&&llphoton_m[0]<160");
  SelectionList category_muon("cat_mu");
  category_muon.AddSelection("objectreq","nphoton>=1&&nlep>=2");
  category_muon.AddSelection("zmassreq","ll_m[0]>81&&ll_m[0]<101");
  category_muon.AddSelection("photonptreq","(photon_pt[0]/llphoton_m[0])>=15.0/110.0");
  category_muon.AddSelection("mllmllgreq","(llphoton_m[0]+ll_m[0])>=185");
  category_muon.AddSelection("photonidcut","photon_id80[0]");
  category_muon.AddSelection("lepptcuts",zg_mu_cuts);
  category_muon.AddSelection("mllgcuts","llphoton_m[0]>100&&llphoton_m[0]<160");
  vector<SelectionList> channels = {category_electron, category_muon};

  //Define systematics
  Systematic syst_altweight("altw",weight*"w_photon");
  Systematic syst_altselection("altsel","zmassreq","ll_m[0]>80&&ll_m[0]<101"); //can make namedfunc that selects signal only using type
  vector<Systematic> systematics = {syst_altweight, syst_altselection};

  //Define parametric PDFs
  RooRealVar rrv_mllg("mllg","llgamma invariant mass",114.0,150.0);
  RooRealVar c0("c0","exponential coefficient",0.0,10.0);
  c0.setVal(1.0);
  RooGenericPdf exp_pdf_el("pdf_background_cat_el_EXP1","exp1_pdf","exp(-1.0*c0*(mllg-114.0)/36.0)",RooArgSet(rrv_mllg,c0));
  RooGenericPdf exp_pdf_mu("pdf_background_cat_mu_EXP1","exp1_pdf","exp(-1.0*c0*(mllg-114.0)/36.0)",RooArgSet(rrv_mllg,c0));

  RooRealVar sigma1("sigma1","sigma1",3., 0.0,10.0);
  RooRealVar mean1("mean1","mean1",130., 120.0,160.0);

  RooRealVar sigma2("sigma2","sigma2",3., 0.0,10.0);
  RooRealVar mean2("mean2","mean2",130., 120.0,160.0);

  RooGaussian gau_pdf_el("pdf_background_cat_el_GAU1","gau1_el_pdf", rrv_mllg, mean1, sigma1);
  RooGaussian gau_pdf_mu("pdf_background_cat_mu_GAU1","gau1_mu_pdf", rrv_mllg, mean2, sigma2);

  vector<RooAbsPdf*> background_pdfs;
  background_pdfs.push_back(&exp_pdf_el);
  background_pdfs.push_back(&exp_pdf_mu);
  background_pdfs.push_back(&gau_pdf_el);
  background_pdfs.push_back(&gau_pdf_mu);

  //Make datacard
  PlotMaker pm;
  pm.multithreaded_ = true;
  pm.min_print_ = true;

  pm.Push<Datacard>("test_datacard", channels, systematics, processes, weight,
      Axis(18, 100.0, 160.0, mllg, "m_{ll#gamma} [GeV]", {}))
      .AddParametricProcess("background",background_pdfs)
      .MakeProcessParametric("htozg",signal_pdfs);
      //.SaveDataAsHist();

  pm.MakePlots(1.0);

  return 0;
}
