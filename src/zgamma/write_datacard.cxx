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
#include "RooMultiPdf.h"
#include "RooRealVar.h"
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
using ZgFunctions::stitch_deathvalley;
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
  NamedFunc trig_and_stitch = (HLT_pass_dilepton||HLT_pass_singlelepton)&&stitch_deathvalley;
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
  vector<RooAbsPdf*> background_pdfs;
  vector<RooAbsPdf*> signal_pdfs;
  vector<RooRealVar*> vars;
  
  for (string category : {"cat_el", "cat_mu"}) {
    //RooRealVar rrv_mllg("mllg_cat_el","mllg cat",100.0,160.0);
    //RooRealVar c0("c0_cat_el","exponential coefficient",0.0,10.0);
    //RooRealVar c1("c1_cat_el","sigmoid offset",100.0,115.0);
    //RooRealVar c2("c2_cat_el","sigmoid width",0.05,20.0);
    //vars.push_back(new RooRealVar(("MH").c_str(),"Higgs mass",120.0,130.0));
    vars.push_back(new RooRealVar(("mllg_"+category).c_str(),"mllg cat",100.0,160.0));
    vars.push_back(new RooRealVar(("c0_"+category).c_str(),"exponential coefficient",0.0,10.0));
    vars.push_back(new RooRealVar(("c1_"+category).c_str(),"sigmoid offset",100.0,115.0));
    vars.push_back(new RooRealVar(("c2_"+category).c_str(),"sigmoid width",0.05,20.0));
    unsigned len_vars = vars.size();
    //RooGenericPdf pdf_bkg_cat(("pdf_background_"+category).c_str(),"bkg_pdf",
    //    "exp(-1.0*@1*(@0-114.0)/36.0)/(1.0+exp(-1.0*@2*(@0-@3)))",
    //    RooArgSet(vars[len_vars-4],vars[len_vars-3],vars[len_vars-2],vars[len_vars-1]));
    //RooGenericPdf pdf_sig_cat(("pdf_htozg_"+category).c_str(),"sig_pdf",
    //    "50.0*TMath::Gaus(@0,125.0,1.5)",
    //    RooArgSet(vars[len_vars-4]));
    RooGenericPdf* pdf_bkg_cat = new RooGenericPdf(("pdf_background_"+category).c_str(),"bkg_pdf",
        "exp(-1.0*@1*(@0-100.0)/60.0)/(1.0+exp(-1.0*@3*(@0-@2)))",
        RooArgSet(*vars[len_vars-4],*vars[len_vars-3],*vars[len_vars-2],*vars[len_vars-1]));
    RooGenericPdf* pdf_sig_cat = new RooGenericPdf(("pdf_htozg_"+category).c_str(),"sig_pdf",
        "TMath::Gaus(@0,125.0,1.5)",
        RooArgSet(*vars[len_vars-4]));
    background_pdfs.push_back(pdf_bkg_cat);
    signal_pdfs.push_back(pdf_sig_cat);
    //background_pdfs.push_back(new RooGenericPdf(("pdf_background_"+category).c_str(),"exp_pdf",
    //    "exp(-1.0*@1*(@0-114.0)/36.0)/(1.0+exp(-1.0*@2*(@0-@3)))",
    //    RooArgSet(vars[len_vars-4],vars[len_vars-3],vars[len_vars-2],vars[len_vars-1])));
    //signal_pdfs.push_back(new RooGenericPdf(("pdf_htozg_"+category).c_str(),"sig_pdf",
    //    "50.0*TMath::Gaus(@0,125.0,1.5)",
    //    RooArgSet(vars[len_vars-4])));
  }
  //RooGenericPdf pdf_bac_el("pdf_background_cat_el","exp_pdf",
  //    "exp(-1.0*@1*(@0-114.0)/36.0)/(1.0+exp(-1.0*@2*(@0-@3)))",
  //    RooArgSet(vars[0],vars[1],vars[2],vars[3]));
  //RooGenericPdf pdf_sig_el("pdf_htozg_cat_el","",
  //    "50.0*TMath::Gaus(@0,125.0,1.5)",
  //    RooArgSet(vars[0]));
  //background_pdfs.push_back(&pdf_bac_el);
  //signal_pdfs.push_back(&pdf_sig_el);
  std::cout << background_pdfs[0]->GetName() << std::endl;
  std::cout << signal_pdfs[0]->GetName() << std::endl;
  std::cout << background_pdfs[1]->GetName() << std::endl;
  std::cout << signal_pdfs[1]->GetName() << std::endl;

  //RooCategory cat("pdf_index","Index of Pdf which is active");
  //RooArgList mypdfs;
  //mypdfs.add(exp_pdf_el);
  //mypdfs.add(exp_pdf_mu);
  //RooMultiPdf multipdf("roomultipdf","All Pdfs",cat,mypdfs);
  //multipdf.GetName();

  //Make datacard
  PlotMaker pm;
  pm.multithreaded_ = true;
  pm.min_print_ = true;

  pm.Push<Datacard>("test_datacard", channels, systematics, processes, weight,
      Axis(18, 100.0, 160.0, mllg, "m_{ll#gamma} [GeV]", {}))
      .AddParametricProcess("background",background_pdfs)
      .MakeProcessParametric("htozg",signal_pdfs);

  pm.MakePlots(1.0);

  return 0;
}
