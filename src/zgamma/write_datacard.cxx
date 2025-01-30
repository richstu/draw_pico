/**
 * Script to write datacard for H->Zgamma analysis
 */

#include <cmath>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "RooAbsPdf.h"
#include "RooArgList.h"
#include "RooBernstein.h"
#include "RooCategory.h"
#include "RooGenericPdf.h"
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
#include "core/RooGaussStepBernstein.hpp"
#include "core/RooDoubleCBFast.hpp"
#include "core/RooMultiPdf.hpp"
#include "zgamma/zg_functions.hpp"
#include "zgamma/zg_utilities.hpp"

using std::make_shared;
using std::set;
using std::shared_ptr;
using std::string;
using std::vector;
using ZgFunctions::llphoton_pttmod;
using ZgFunctions::max_lep_miniso;
using ZgFunctions::trig_plateau_cuts;
using ZgFunctions::w_years;
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
  vector<shared_ptr<Process>> processes = ZgSampleLoader() //.SetMacro("YEARS",{"2018"}).Verbose(true)
        .LoadSamples("txt/samples_zgamma.txt","DatacardMC");

  //Define NamedFuncs
  NamedFunc mllg = NamedFunc("llphoton_m[0]").Name("mllg");
  NamedFunc mllg_range_cut = NamedFunc("llphoton_m[0]>105&&llphoton_m[0]<160");

  shared_ptr<MVAWrapper> kinematic_bdt = KinematicBdt();
  NamedFunc ggh4 = category_ggh4(kinematic_bdt);
  NamedFunc ggh3 = category_ggh3(kinematic_bdt);
  NamedFunc ggh2 = category_ggh2(kinematic_bdt);
  NamedFunc ggh1 = category_ggh1(kinematic_bdt);

  shared_ptr<MVAWrapper> vbf_bdt = VbfBdt();
  NamedFunc vbf3 = category_vbf3(vbf_bdt);
  NamedFunc vbf2 = category_vbf2(vbf_bdt);
  NamedFunc vbf1 = category_vbf1(vbf_bdt);

  //weight with some clipping/regularization
  const NamedFunc weight_reg("weight_reg",[](const Baby &b) -> NamedFunc::ScalarType{
    if (fabs(b.weight())>5) return 0.0;
    return b.weight();
  });

  //Define weight
  NamedFunc weight(w_years*weight_reg); 

  //Define channels
  SelectionList category_ggh4("cat_ggh4");
  category_ggh4.AddSelection("objectreq","nphoton>=1&&nlep==2&&njet<2&&met<90");
  category_ggh4.AddSelection("lepptcuts",trig_plateau_cuts);
  category_ggh4.AddSelection("zmassreq","ll_m[0]>81&&ll_m[0]<101");
  category_ggh4.AddSelection("zchargereq","ll_charge[0]==0");
  category_ggh4.AddSelection("photonptreq","(photon_pt[0]/llphoton_m[0])>=15.0/110.0");
  category_ggh4.AddSelection("mllmllgreq","(llphoton_m[0]+ll_m[0])>=185");
  category_ggh4.AddSelection("photonidreq","photon_id80[0]");
  category_ggh4.AddSelection("mllgrange",mllg_range_cut);
  category_ggh4.AddSelection("bdtcuts",ggh4);
  SelectionList category_ggh3("cat_ggh3");
  category_ggh3.AddSelection("objectreq","nphoton>=1&&nlep==2&&njet<2&&met<90");
  category_ggh3.AddSelection("lepptcuts",trig_plateau_cuts);
  category_ggh3.AddSelection("zmassreq","ll_m[0]>81&&ll_m[0]<101");
  category_ggh3.AddSelection("zchargereq","ll_charge[0]==0");
  category_ggh3.AddSelection("photonptreq","(photon_pt[0]/llphoton_m[0])>=15.0/110.0");
  category_ggh3.AddSelection("mllmllgreq","(llphoton_m[0]+ll_m[0])>=185");
  category_ggh3.AddSelection("photonidreq","photon_id80[0]");
  category_ggh3.AddSelection("mllgrange",mllg_range_cut);
  category_ggh3.AddSelection("bdtcuts",ggh3);
  SelectionList category_ggh2("cat_ggh2");
  category_ggh2.AddSelection("objectreq","nphoton>=1&&nlep==2&&njet<2&&met<90");
  category_ggh2.AddSelection("lepptcuts",trig_plateau_cuts);
  category_ggh2.AddSelection("zmassreq","ll_m[0]>81&&ll_m[0]<101");
  category_ggh2.AddSelection("zchargereq","ll_charge[0]==0");
  category_ggh2.AddSelection("photonptreq","(photon_pt[0]/llphoton_m[0])>=15.0/110.0");
  category_ggh2.AddSelection("mllmllgreq","(llphoton_m[0]+ll_m[0])>=185");
  category_ggh2.AddSelection("photonidreq","photon_id80[0]");
  category_ggh2.AddSelection("mllgrange",mllg_range_cut);
  category_ggh2.AddSelection("bdtcuts",ggh2);
  SelectionList category_ggh1("cat_ggh1");
  category_ggh1.AddSelection("objectreq","nphoton>=1&&nlep==2&&njet<2&&met<90");
  category_ggh1.AddSelection("lepptcuts",trig_plateau_cuts);
  category_ggh1.AddSelection("zmassreq","ll_m[0]>81&&ll_m[0]<101");
  category_ggh1.AddSelection("zchargereq","ll_charge[0]==0");
  category_ggh1.AddSelection("photonptreq","(photon_pt[0]/llphoton_m[0])>=15.0/110.0");
  category_ggh1.AddSelection("mllmllgreq","(llphoton_m[0]+ll_m[0])>=185");
  category_ggh1.AddSelection("photonidreq","photon_id80[0]");
  category_ggh1.AddSelection("mllgrange",mllg_range_cut);
  category_ggh1.AddSelection("bdtcuts",ggh1);
  SelectionList category_vbf3("cat_vbf3");
  category_vbf3.AddSelection("objectreq","nphoton>=1&&nlep==2&&njet>=2&&nbm==0");
  category_vbf3.AddSelection("lepptcuts",trig_plateau_cuts);
  category_vbf3.AddSelection("zmassreq","ll_m[0]>81&&ll_m[0]<101");
  category_vbf3.AddSelection("zchargereq","ll_charge[0]==0");
  category_vbf3.AddSelection("photonptreq","(photon_pt[0]/llphoton_m[0])>=15.0/110.0");
  category_vbf3.AddSelection("mllmllgreq","(llphoton_m[0]+ll_m[0])>=185");
  category_vbf3.AddSelection("photonidreq","photon_id80[0]");
  category_vbf3.AddSelection("mllgrange",mllg_range_cut);
  category_vbf3.AddSelection("bdtcuts",vbf3);
  SelectionList category_vbf2("cat_vbf2");
  category_vbf2.AddSelection("objectreq","nphoton>=1&&nlep==2&&njet>=2&&nbm==0");
  category_vbf2.AddSelection("lepptcuts",trig_plateau_cuts);
  category_vbf2.AddSelection("zmassreq","ll_m[0]>81&&ll_m[0]<101");
  category_vbf2.AddSelection("zchargereq","ll_charge[0]==0");
  category_vbf2.AddSelection("photonptreq","(photon_pt[0]/llphoton_m[0])>=15.0/110.0");
  category_vbf2.AddSelection("mllmllgreq","(llphoton_m[0]+ll_m[0])>=185");
  category_vbf2.AddSelection("photonidreq","photon_id80[0]");
  category_vbf2.AddSelection("mllgrange",mllg_range_cut);
  category_vbf2.AddSelection("bdtcuts",vbf2);
  SelectionList category_vbf1("cat_vbf1");
  category_vbf1.AddSelection("objectreq","nphoton>=1&&nlep==2&&njet>=2&&nbm==0");
  category_vbf1.AddSelection("lepptcuts",trig_plateau_cuts);
  category_vbf1.AddSelection("zmassreq","ll_m[0]>81&&ll_m[0]<101");
  category_vbf1.AddSelection("zchargereq","ll_charge[0]==0");
  category_vbf1.AddSelection("photonptreq","(photon_pt[0]/llphoton_m[0])>=15.0/110.0");
  category_vbf1.AddSelection("mllmllgreq","(llphoton_m[0]+ll_m[0])>=185");
  category_vbf1.AddSelection("photonidreq","photon_id80[0]");
  category_vbf1.AddSelection("mllgrange",mllg_range_cut);
  category_vbf1.AddSelection("bdtcuts",vbf1);
  SelectionList category_vhmet("cat_vhmet");
  category_vhmet.AddSelection("objectreq","nphoton>=1&&nlep==2&&njet<2&&met>90");
  category_vhmet.AddSelection("lepptcuts",trig_plateau_cuts);
  category_vhmet.AddSelection("zmassreq","ll_m[0]>81&&ll_m[0]<101");
  category_vhmet.AddSelection("zchargereq","ll_charge[0]==0");
  category_vhmet.AddSelection("photonptreq","(photon_pt[0]/llphoton_m[0])>=15.0/110.0");
  category_vhmet.AddSelection("mllmllgreq","(llphoton_m[0]+ll_m[0])>=185");
  category_vhmet.AddSelection("photonidreq","photon_id80[0]");
  category_vhmet.AddSelection("mllgrange",mllg_range_cut);
  category_vhmet.AddSelection("pttmodreq",llphoton_pttmod<60.0);
  category_vhmet.AddSelection("ptllgreq","llphoton_pt[0]/llphoton_m[0]>0.4");
  SelectionList category_vh3l("cat_vh3l");
  category_vh3l.AddSelection("objectreq","nphoton>=1&&nlep>=3&&nbm==0");
  category_vh3l.AddSelection("lepptcuts",trig_plateau_cuts);
  category_vh3l.AddSelection("zmassreq","ll_m[0]>81&&ll_m[0]<101");
  category_vh3l.AddSelection("zchargereq","ll_charge[0]==0");
  category_vh3l.AddSelection("photonptreq","(photon_pt[0]/llphoton_m[0])>=15.0/110.0");
  category_vh3l.AddSelection("mllmllgreq","(llphoton_m[0]+ll_m[0])>=185");
  category_vh3l.AddSelection("photonidreq","photon_id80[0]");
  category_vh3l.AddSelection("mllgrange",mllg_range_cut);
  category_vh3l.AddSelection("minisoreq",max_lep_miniso<0.15);
  category_vh3l.AddSelection("ptllgreq","llphoton_pt[0]/llphoton_m[0]>0.3");
  SelectionList category_tthhad("cat_tthhad");
  category_tthhad.AddSelection("objectreq","nphoton>=1&&nlep==2&&njet>=5&&nbm>=1");
  category_tthhad.AddSelection("lepptcuts",trig_plateau_cuts);
  category_tthhad.AddSelection("zmassreq","ll_m[0]>81&&ll_m[0]<101");
  category_tthhad.AddSelection("zchargereq","ll_charge[0]==0");
  category_tthhad.AddSelection("photonptreq","(photon_pt[0]/llphoton_m[0])>=15.0/110.0");
  category_tthhad.AddSelection("mllmllgreq","(llphoton_m[0]+ll_m[0])>=185");
  category_tthhad.AddSelection("photonidreq","photon_id80[0]");
  category_tthhad.AddSelection("mllgrange",mllg_range_cut);
  SelectionList category_tthlep("cat_tthlep");
  category_tthlep.AddSelection("objectreq","nphoton>=1&&nlep>=3&&njet>=3&&nbm>=1");
  category_tthlep.AddSelection("lepptcuts",trig_plateau_cuts);
  category_tthlep.AddSelection("zmassreq","ll_m[0]>81&&ll_m[0]<101");
  category_tthlep.AddSelection("zchargereq","ll_charge[0]==0");
  category_tthlep.AddSelection("photonptreq","(photon_pt[0]/llphoton_m[0])>=15.0/110.0");
  category_tthlep.AddSelection("mllmllgreq","(llphoton_m[0]+ll_m[0])>=185");
  category_tthlep.AddSelection("photonidreq","photon_id80[0]");
  category_tthlep.AddSelection("mllgrange",mllg_range_cut);
  category_tthlep.AddSelection("minisoreq",max_lep_miniso<0.15);
  vector<SelectionList> channels = {category_ggh4,category_ggh3,category_ggh2,
                                    category_ggh1,category_vbf3,category_vbf2,
                                    category_vbf1,category_vh3l,category_vhmet,
                                    category_tthhad,category_tthlep};


  //Define systematics
  //Systematic syst_altweight("altw",weight*"w_photon");
  //Systematic syst_altselection("altsel","zmassreq","ll_m[0]>80&&ll_m[0]<101"); 
  //can make namedfunc that selects signal only using type
  //vector<Systematic> systematics = {syst_altweight, syst_altselection};
  vector<Systematic> systematics = {};

  //Define parametric PDFs
  vector<shared_ptr<RooAbsPdf>> background_pdfs;
  vector<shared_ptr<RooAbsPdf>> aux_pdfs;
  //vector<vector<shared_ptr<RooAbsPdf>>> background_pdfs;
  vector<shared_ptr<RooAbsPdf>> signal_pdfs;
  vector<shared_ptr<RooRealVar>> vars;
  
  for (string category : {"cat_ggh4", "cat_ggh3", "cat_ggh2", "cat_ggh1", 
                          "cat_vbf3", "cat_vbf2", "cat_vbf1", "cat_vh3l", 
                          "cat_vhmet", "cat_tthhad", "cat_tthlep"}) {
  //for (string category : {"cat_ggh4"}) {
    shared_ptr<RooRealVar> pdf_mllg = make_shared<RooRealVar>(
        ("mllg_"+category).c_str(),"mllg cat",105.0,160.0); 
    vars.push_back(pdf_mllg);

    AddGaussStepExponential(background_pdfs, aux_pdfs, vars, pdf_mllg, category, 2);

    shared_ptr<RooRealVar> dscb_mu = make_shared<RooRealVar>(
        ("dscb_mu_"+category).c_str(),"Signal DSCB gaussian mean",120.0,130.0); 
    vars.push_back(dscb_mu);
    shared_ptr<RooRealVar> dscb_sigma = make_shared<RooRealVar>(
        ("dscb_sigma_"+category).c_str(),"Signal DSCB gaussian sigma",0.2,4.0); 
    vars.push_back(dscb_sigma);
    shared_ptr<RooRealVar> dscb_alphal = make_shared<RooRealVar>(
        ("dscb_alphal_"+category).c_str(),"Signal DSCB alphaL",0.5,4.0); 
    vars.push_back(dscb_alphal);
    shared_ptr<RooRealVar> dscb_alphar = make_shared<RooRealVar>(
        ("dscb_alphar_"+category).c_str(),"Signal DSCB alphaR",0.5,4.0); 
    vars.push_back(dscb_alphar);
    shared_ptr<RooRealVar> dscb_nl = make_shared<RooRealVar>(
        ("dscb_nl_"+category).c_str(),"Signal DSCB nL",0.2,10.0); 
    vars.push_back(dscb_nl);
    shared_ptr<RooRealVar> dscb_nr = make_shared<RooRealVar>(
        ("dscb_nr_"+category).c_str(),"Signal DSCB nR",0.2,10.0); 
    vars.push_back(dscb_nr);
    shared_ptr<RooDoubleCBFast> pdf_sig_cat = make_shared<RooDoubleCBFast>(
        ("pdf_htozg_"+category).c_str(),"sig_pdf",
        *pdf_mllg, *dscb_mu, *dscb_sigma, *dscb_alphal, *dscb_nl, 
        *dscb_alphar, *dscb_nr);

    //background_pdfs.push_back(vector<RooAbsPdf*>{pdf_bkg_cat_bern3,pdf_bkg_cat_exp1});
    //background_pdfs.push_back(pdf_bkg_cat_bern4);
    signal_pdfs.push_back(pdf_sig_cat);
  }

  //RooCategory cat("pdf_index","Index of Pdf which is active");
  //RooArgList mypdfs;
  //mypdfs.add(exp_pdf_el);
  //mypdfs.add(exp_pdf_mu);
  //RooMultiPdf multipdf("roomultipdf","All Pdfs",cat,mypdfs);
  //multipdf.GetName();

  //Make datacard
  PlotMaker pm;
  pm.multithreaded_ = false;
  pm.min_print_ = true;

  pm.Push<Datacard>("test_datacard2", channels, systematics, processes, weight,
      Axis(55, 105.0, 160.0, mllg, "m_{ll#gamma} [GeV]", {}))
      .SaveDataAsHist()
      .AddParametricProcess("background",background_pdfs)
      .MakeProcessParametric("htozgamma",signal_pdfs);

  pm.MakePlots(1.0);

  return 0;
}
