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
#include "zgamma/zg_utilities.hpp"

using std::make_shared;
using std::set;
using std::shared_ptr;
using std::string;
using std::vector;
using ZgFunctions::HLT_pass_dilepton;
using ZgFunctions::HLT_pass_singlelepton;
using ZgFunctions::stitch_deathvalley;
using ZgFunctions::w_years;
using ZgUtilities::ZgSampleLoader;
using ZgUtilities::KinematicBdt;
using ZgUtilities::category_ggh4;
using ZgUtilities::category_ggh3;
using ZgUtilities::category_ggh2;
using ZgUtilities::category_ggh1;
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
  NamedFunc lep_pt_cuts = NamedFunc("(nel>=2&&el_pt[0]>25&&el_pt[1]>15)||(nmu>=2&&mu_pt[0]>20&&mu_pt[1]>10)").Name("lep_pt_cuts");
  shared_ptr<MVAWrapper> kinematic_bdt = KinematicBdt();
  NamedFunc ggh4 = category_ggh4(kinematic_bdt);
  NamedFunc ggh3 = category_ggh3(kinematic_bdt);
  NamedFunc ggh2 = category_ggh2(kinematic_bdt);
  NamedFunc ggh1 = category_ggh1(kinematic_bdt);

  //Define weight
  NamedFunc weight(w_years*"weight"); 

  //Define channels
  SelectionList category_ggh4("cat_ggh4");
  category_ggh4.AddSelection("objectreq","nphoton>=1&&nlep>=2");
  category_ggh4.AddSelection("zmassreq","ll_m[0]>81&&ll_m[0]<101");
  category_ggh4.AddSelection("photonptreq","(photon_pt[0]/llphoton_m[0])>=15.0/110.0");
  category_ggh4.AddSelection("mllmllgreq","(llphoton_m[0]+ll_m[0])>=185");
  category_ggh4.AddSelection("photonidcut","photon_id80[0]");
  category_ggh4.AddSelection("lepptcuts",lep_pt_cuts);
  category_ggh4.AddSelection("mllgcuts","llphoton_m[0]>110&&llphoton_m[0]<160");
  category_ggh4.AddSelection("bdtcuts",ggh4);
  SelectionList category_ggh3("cat_ggh3");
  category_ggh3.AddSelection("objectreq","nphoton>=1&&nlep>=2");
  category_ggh3.AddSelection("zmassreq","ll_m[0]>81&&ll_m[0]<101");
  category_ggh3.AddSelection("photonptreq","(photon_pt[0]/llphoton_m[0])>=15.0/110.0");
  category_ggh3.AddSelection("mllmllgreq","(llphoton_m[0]+ll_m[0])>=185");
  category_ggh3.AddSelection("photonidcut","photon_id80[0]");
  category_ggh3.AddSelection("lepptcuts",lep_pt_cuts);
  category_ggh3.AddSelection("mllgcuts","llphoton_m[0]>110&&llphoton_m[0]<160");
  category_ggh3.AddSelection("bdtcuts",ggh3);
  SelectionList category_ggh2("cat_ggh2");
  category_ggh2.AddSelection("objectreq","nphoton>=1&&nlep>=2");
  category_ggh2.AddSelection("zmassreq","ll_m[0]>81&&ll_m[0]<101");
  category_ggh2.AddSelection("photonptreq","(photon_pt[0]/llphoton_m[0])>=15.0/110.0");
  category_ggh2.AddSelection("mllmllgreq","(llphoton_m[0]+ll_m[0])>=185");
  category_ggh2.AddSelection("photonidcut","photon_id80[0]");
  category_ggh2.AddSelection("lepptcuts",lep_pt_cuts);
  category_ggh2.AddSelection("mllgcuts","llphoton_m[0]>110&&llphoton_m[0]<160");
  category_ggh2.AddSelection("bdtcuts",ggh2);
  SelectionList category_ggh1("cat_ggh1");
  category_ggh1.AddSelection("objectreq","nphoton>=1&&nlep>=2");
  category_ggh1.AddSelection("zmassreq","ll_m[0]>81&&ll_m[0]<101");
  category_ggh1.AddSelection("photonptreq","(photon_pt[0]/llphoton_m[0])>=15.0/110.0");
  category_ggh1.AddSelection("mllmllgreq","(llphoton_m[0]+ll_m[0])>=185");
  category_ggh1.AddSelection("photonidcut","photon_id80[0]");
  category_ggh1.AddSelection("lepptcuts",lep_pt_cuts);
  category_ggh1.AddSelection("mllgcuts","llphoton_m[0]>110&&llphoton_m[0]<160");
  category_ggh1.AddSelection("bdtcuts",ggh1);
  vector<SelectionList> channels = {category_ggh4,category_ggh3,category_ggh2,category_ggh1};

  //Define systematics
  //Systematic syst_altweight("altw",weight*"w_photon");
  //Systematic syst_altselection("altsel","zmassreq","ll_m[0]>80&&ll_m[0]<101"); 
  //can make namedfunc that selects signal only using type
  //vector<Systematic> systematics = {syst_altweight, syst_altselection};
  vector<Systematic> systematics = {};

  //Define parametric PDFs
  //TODO switch everything to shared_ptr
  vector<RooAbsPdf*> background_pdfs;
  //vector<vector<RooAbsPdf*>> background_pdfs;
  vector<RooAbsPdf*> signal_pdfs;
  vector<shared_ptr<RooRealVar>> vars;
  
  for (string category : {"cat_ggh4", "cat_ggh3", "cat_ggh2", "cat_ggh1"}) {
    shared_ptr<RooRealVar> pdf_mllg = make_shared<RooRealVar>(("mllg_"+category).c_str(),"mllg cat",110.0,160.0); vars.push_back(pdf_mllg);

    //shared_ptr<RooRealVar> b5_c0 = make_shared<RooRealVar>(("b5_c0_"+category).c_str(),"bern 5 coefficient 0",-1.0,1.0); vars.push_back(b5_c0);
    //shared_ptr<RooRealVar> b5_c1 = make_shared<RooRealVar>(("b5_c1_"+category).c_str(),"bern 5 coefficient 1",-1.0,1.0); vars.push_back(b5_c1);
    //shared_ptr<RooRealVar> b5_c2 = make_shared<RooRealVar>(("b5_c2_"+category).c_str(),"bern 5 coefficient 2",-1.0,1.0); vars.push_back(b5_c2);
    //shared_ptr<RooRealVar> b5_c3 = make_shared<RooRealVar>(("b5_c3_"+category).c_str(),"bern 5 coefficient 3",-1.0,1.0); vars.push_back(b5_c3);
    //shared_ptr<RooRealVar> b5_c4 = make_shared<RooRealVar>(("b5_c4_"+category).c_str(),"bern 5 coefficient 4",-1.0,1.0); vars.push_back(b5_c4);
    //shared_ptr<RooRealVar> b5_c5 = make_shared<RooRealVar>(("b5_c5_"+category).c_str(),"bern 5 coefficient 5",-1.0,1.0); vars.push_back(b5_c5);
    //shared_ptr<RooRealVar> b5_so = make_shared<RooRealVar>(("b5_so_"+category).c_str(),"bern 5 sigmoid offset",90.0,125.0); vars.push_back(b5_so);
    //shared_ptr<RooRealVar> b5_sw = make_shared<RooRealVar>(("b5_sw_"+category).c_str(),"bern 5 sigmoid width",0.005,30.0); vars.push_back(b5_sw);
    //b5_c0->setVal(1.0);
    //b5_c1->setVal(0.1);
    //b5_c2->setVal(0.1);
    //b5_c3->setVal(0.1);
    //b5_c4->setVal(0.1);
    //b5_c5->setVal(0.1);
    //b5_so->setVal(108);
    //b5_sw->setVal(0.2);
    //RooGenericPdf* pdf_bkg_cat_bern5 = new RooGenericPdf(("pdf_background_"+category+"_bern5").c_str(),"bkg_pdf",
    //    "(@1*pow((165-@0)/60,5)+@2*(5*(@0-105)/60*pow((165-@0)/60,4))+@3*(10*pow((@0-105)/60,2)*pow((165-@0)/60,3))+@4*(10*pow((@0-105)/60,3)*pow((165-@0)/60,2))+@5*(5*pow((@0-105)/60,4)*(165-@0)/60)+@6*pow((@0-105)/60,5))/(1.0+exp(-1.0*@8*(@0-@7)))",
    //    RooArgSet(*pdf_mllg,*b5_c0,*b5_c1,*b5_c2,*b5_c3,*b5_c4,*b5_c5,*b5_so,*b5_sw));

    RooArgSet bern6_argset;
    bern6_argset.add(*pdf_mllg);
    shared_ptr<RooRealVar> b6_c0 = make_shared<RooRealVar>(("b6_c0_"+category).c_str(),"bern 6 coefficient 0",-1.0,1.0); vars.push_back(b6_c0); bern6_argset.add(*b6_c0);
    shared_ptr<RooRealVar> b6_c1 = make_shared<RooRealVar>(("b6_c1_"+category).c_str(),"bern 6 coefficient 1",-1.0,1.0); vars.push_back(b6_c1); bern6_argset.add(*b6_c1);
    shared_ptr<RooRealVar> b6_c2 = make_shared<RooRealVar>(("b6_c2_"+category).c_str(),"bern 6 coefficient 2",-1.0,1.0); vars.push_back(b6_c2); bern6_argset.add(*b6_c2);
    shared_ptr<RooRealVar> b6_c3 = make_shared<RooRealVar>(("b6_c3_"+category).c_str(),"bern 6 coefficient 3",-1.0,1.0); vars.push_back(b6_c3); bern6_argset.add(*b6_c3);
    shared_ptr<RooRealVar> b6_c4 = make_shared<RooRealVar>(("b6_c4_"+category).c_str(),"bern 6 coefficient 4",-1.0,1.0); vars.push_back(b6_c4); bern6_argset.add(*b6_c4);
    shared_ptr<RooRealVar> b6_c5 = make_shared<RooRealVar>(("b6_c5_"+category).c_str(),"bern 6 coefficient 5",-1.0,1.0); vars.push_back(b6_c5); bern6_argset.add(*b6_c5);
    shared_ptr<RooRealVar> b6_c6 = make_shared<RooRealVar>(("b6_c6_"+category).c_str(),"bern 6 coefficient 6",-1.0,1.0); vars.push_back(b6_c6); bern6_argset.add(*b6_c6);
    shared_ptr<RooRealVar> b6_so = make_shared<RooRealVar>(("b6_so_"+category).c_str(),"bern 6 sigmoid offset",90.0,125.0); vars.push_back(b6_so); bern6_argset.add(*b6_so);
    shared_ptr<RooRealVar> b6_sw = make_shared<RooRealVar>(("b6_sw_"+category).c_str(),"bern 6 sigmoid width",0.005,30.0); vars.push_back(b6_sw); bern6_argset.add(*b6_sw);
    b6_c0->setVal(1.0);
    b6_c1->setVal(0.1);
    b6_c2->setVal(0.1);
    b6_c3->setVal(0.1);
    b6_c4->setVal(0.1);
    b6_c5->setVal(0.0);
    b6_c6->setVal(0.0);
    b6_so->setVal(110);
    b6_sw->setVal(0.25);
    RooGenericPdf* pdf_bkg_cat_bern6 = new RooGenericPdf(("pdf_background_"+category+"_bern6").c_str(),"bkg_pdf",
        "(@1*pow((160-@0)/50,6)+@2*(6*(@0-110)/50*pow((160-@0)/50,5))+@3*(15*pow((@0-110)/50,2)*pow((160-@0)/50,4))+@4*(20*pow((@0-110)/50,3)*pow((160-@0)/50,3))+@5*(15*pow((@0-110)/50,4)*pow((160-@0)/50,2))+@6*(6*pow((@0-110)/50,5)*(160-@0)/50)+@7*pow((@0-110)/50,6))/(1.0+exp(-1.0*@9*(@0-@8)))",
        bern6_argset);

    //RooArgSet laurent2_argset;
    //laurent2_argset.add(*pdf_mllg);
    //shared_ptr<RooRealVar> l2_c0 = make_shared<RooRealVar>(("l2_c0_"+category).c_str(),"Laurent 5 coefficient 0",0.0,1.0); vars.push_back(l2_c0); laurent2_argset.add(*l2_c0);
    //shared_ptr<RooRealVar> l2_c1 = make_shared<RooRealVar>(("l2_c1_"+category).c_str(),"Laurent 5 coefficient 1",0.0,1.0); vars.push_back(l2_c1); laurent2_argset.add(*l2_c1);
    //shared_ptr<RooRealVar> l2_c2 = make_shared<RooRealVar>(("l2_c2_"+category).c_str(),"Laurent 5 coefficient 2",0.0,1.0); vars.push_back(l2_c2); laurent2_argset.add(*l2_c2);
    //shared_ptr<RooRealVar> l2_o1 = make_shared<RooRealVar>(("l2_o1_"+category).c_str(),"Laurent 5 sigmoid 1 offset",90.0,115.0); vars.push_back(l2_o1); laurent2_argset.add(*l2_o1);
    //shared_ptr<RooRealVar> l2_w1 = make_shared<RooRealVar>(("l2_w1_"+category).c_str(),"Laurent 5 sigmoid 1 width",0.005,30.0); vars.push_back(l2_w1); laurent2_argset.add(*l2_w1);
    //shared_ptr<RooRealVar> l2_i1 = make_shared<RooRealVar>(("l2_i1_"+category).c_str(),"Laurent 5 sigmoid 1 pedestal",0.0,1.0); vars.push_back(l2_i1); laurent2_argset.add(*l2_i1);
    //shared_ptr<RooRealVar> l2_o2 = make_shared<RooRealVar>(("l2_o2_"+category).c_str(),"Laurent 5 sigmoid 2 offset",90.0,115.0); vars.push_back(l2_o2); laurent2_argset.add(*l2_o2);
    //shared_ptr<RooRealVar> l2_w2 = make_shared<RooRealVar>(("l2_w2_"+category).c_str(),"Laurent 5 sigmoid 2 width",0.005,30.0); vars.push_back(l2_w2); laurent2_argset.add(*l2_w2);
    //shared_ptr<RooRealVar> l2_i2 = make_shared<RooRealVar>(("l2_i2_"+category).c_str(),"Laurent 5 sigmoid 2 pedestal",0.0,1.0); vars.push_back(l2_i2); laurent2_argset.add(*l2_i2);
    //l2_c0->setVal(1.0);
    //l2_c1->setVal(0.1);
    //l2_c2->setVal(0.1);
    //l2_o1->setVal(105);
    //l2_w1->setVal(0.03);
    //l2_i1->setVal(0.0);
    //l2_o1->setVal(110);
    //l2_w1->setVal(0.03);
    //l2_i1->setVal(0.1);
    //RooGenericPdf* pdf_bkg_cat_laurent2 = new RooGenericPdf(("pdf_background_"+category+"_laurent2").c_str(),"bkg_pdf",
    //    "(@1/pow((@0-45)/60,3)+@2/pow((@0-45)/60,4)+@3/pow((@0-45)/60,5))*(@6+(1.0-@6)/(1.0+exp(-1.0*@5*(@0-@4))))*(@9+(1.0-@9)/(1.0+exp(-1.0*@8*(@0-@7))))",
    //    laurent2_argset);

    //RooArgSet doublelaurent_argset;
    //doublelaurent_argset.add(*pdf_mllg);
    //shared_ptr<RooRealVar> dl_c0 = make_shared<RooRealVar>(("dl_c0_"+category).c_str(),"Double Laurent 1 coefficient 0,0",-1.0,1.0); vars.push_back(dl_c0); doublelaurent_argset.add(*dl_c0);
    //shared_ptr<RooRealVar> dl_c1 = make_shared<RooRealVar>(("dl_c1_"+category).c_str(),"Double Laurent 1 coefficient 0,1",-1.0,1.0); vars.push_back(dl_c1); doublelaurent_argset.add(*dl_c1);
    //shared_ptr<RooRealVar> dl_c2 = make_shared<RooRealVar>(("dl_c2_"+category).c_str(),"Double Laurent 1 coefficient 1,0",-1.0,1.0); vars.push_back(dl_c2); doublelaurent_argset.add(*dl_c2);
    //shared_ptr<RooRealVar> dl_c3 = make_shared<RooRealVar>(("dl_c3_"+category).c_str(),"Double Laurent 1 coefficient 1,1",-1.0,1.0); vars.push_back(dl_c3); doublelaurent_argset.add(*dl_c3);
    //shared_ptr<RooRealVar> dl_o0 = make_shared<RooRealVar>(("dl_o0_"+category).c_str(),"Double Laurent 1 sigmoid offset 0",90.0,115.0); vars.push_back(dl_o0); doublelaurent_argset.add(*dl_o0);
    //shared_ptr<RooRealVar> dl_o1 = make_shared<RooRealVar>(("dl_o1_"+category).c_str(),"Double Laurent 1 sigmoid offset 1",90.0,115.0); vars.push_back(dl_o1); doublelaurent_argset.add(*dl_o1);
    //shared_ptr<RooRealVar> dl_w0 = make_shared<RooRealVar>(("dl_w0_"+category).c_str(),"Double Laurent 1 sigmoid width 0",0.005,20.0); vars.push_back(dl_w0); doublelaurent_argset.add(*dl_w0);
    //shared_ptr<RooRealVar> dl_w1 = make_shared<RooRealVar>(("dl_w1_"+category).c_str(),"Double Laurent 1 sigmoid width 1",0.005,20.0); vars.push_back(dl_w1); doublelaurent_argset.add(*dl_w1);
    //dl_c0->setVal(0.5);
    //dl_c1->setVal(0.5);
    //dl_c2->setVal(0.5);
    //dl_c3->setVal(0.5);
    //dl_o0->setVal(105);
    //dl_o1->setVal(105);
    //dl_w0->setVal(0.03);
    //dl_w1->setVal(0.03);
    //RooGenericPdf* pdf_bkg_cat_doublelaurent1 = new RooGenericPdf(("pdf_background_"+category+"_doublelaurent1").c_str(),"bkg_pdf",
    //    "(@1/pow((@0-45)/60,4)+@2/pow((@0-45)/60,5))/(1.0+exp(-1.0*@7*(@0-@5)))+(@3/pow((@0-45)/60,4)+@4/pow((@0-45)/60,5))/(1.0+exp(-1.0*@8*(@0-@6)))",
    //    doublelaurent_argset);
    
    //RooArgSet doubleexp_argset;
    //doubleexp_argset.add(*pdf_mllg);
    //shared_ptr<RooRealVar> de_c0 = make_shared<RooRealVar>(("de_c0_"+category).c_str(),"Double Exponential 1 coefficient 0,0",-1.0,1.0); vars.push_back(de_c0); doubleexp_argset.add(*de_c0);
    //shared_ptr<RooRealVar> de_c1 = make_shared<RooRealVar>(("de_c1_"+category).c_str(),"Double Exponential 1 coefficient 0,1",-1.0,1.0); vars.push_back(de_c1); doubleexp_argset.add(*de_c1);
    //shared_ptr<RooRealVar> de_p0 = make_shared<RooRealVar>(("de_p0_"+category).c_str(),"Double Exponential 1 width 1,0",0.0,20.0); vars.push_back(de_p0); doubleexp_argset.add(*de_p0);
    //shared_ptr<RooRealVar> de_p1 = make_shared<RooRealVar>(("de_p1_"+category).c_str(),"Double Exponential 1 width 1,1",0.0,20.0); vars.push_back(de_p1); doubleexp_argset.add(*de_p1);
    //shared_ptr<RooRealVar> de_o0 = make_shared<RooRealVar>(("de_o0_"+category).c_str(),"Double Exponential 1 sigmoid offset 0",90.0,120.0); vars.push_back(de_o0); doubleexp_argset.add(*de_o0);
    //shared_ptr<RooRealVar> de_o1 = make_shared<RooRealVar>(("de_o1_"+category).c_str(),"Double Exponential 1 sigmoid offset 1",90.0,120.0); vars.push_back(de_o1); doubleexp_argset.add(*de_o1);
    //shared_ptr<RooRealVar> de_w0 = make_shared<RooRealVar>(("de_w0_"+category).c_str(),"Double Exponential 1 sigmoid width 0",0.001,20.0); vars.push_back(de_w0); doubleexp_argset.add(*de_w0);
    //shared_ptr<RooRealVar> de_w1 = make_shared<RooRealVar>(("de_w1_"+category).c_str(),"Double Exponential 1 sigmoid width 1",0.001,20.0); vars.push_back(de_w1); doubleexp_argset.add(*de_w1);
    //shared_ptr<RooRealVar> de_pd = make_shared<RooRealVar>(("de_pd_"+category).c_str(),"Double Exponential 1 pedestal",-1.0,1.0); vars.push_back(de_pd); doubleexp_argset.add(*de_pd);
    //de_c0->setVal(0.5);
    //de_c1->setVal(0.1);
    //de_p0->setVal(4.0);
    //de_p1->setVal(9.0);
    //de_o0->setVal(105);
    //de_o1->setVal(110);
    //de_w0->setVal(0.8);
    //de_w1->setVal(0.5);
    //de_pd->setVal(0.0);
    //RooGenericPdf* pdf_bkg_cat_doubleexp1 = new RooGenericPdf(("pdf_background_"+category+"_doubleexp1").c_str(),"bkg_pdf",
    //    "(@1*exp(-1.0*@3*(@0-105)/60))/(1.0+exp(-1.0*@7*(@0-@5)))+(@2*exp(-1.0*@4*(@0-105)/60))/(1.0+exp(-1.0*@8*(@0-@6)))+@9",
    //    doubleexp_argset);

    //RooArgSet doubleexp_argset;
    //doubleexp_argset.add(*pdf_mllg);
    //shared_ptr<RooRealVar> de_c0 = make_shared<RooRealVar>(("de_c0_"+category).c_str(),"Double Exponential 1 coefficient 0,0",-1.0,1.0); vars.push_back(de_c0); doubleexp_argset.add(*de_c0);
    //shared_ptr<RooRealVar> de_c1 = make_shared<RooRealVar>(("de_c1_"+category).c_str(),"Double Exponential 1 coefficient 0,1",-1.0,1.0); vars.push_back(de_c1); doubleexp_argset.add(*de_c1);
    //shared_ptr<RooRealVar> de_p0 = make_shared<RooRealVar>(("de_p0_"+category).c_str(),"Double Exponential 1 width 1,0",0.0,20.0); vars.push_back(de_p0); doubleexp_argset.add(*de_p0);
    //shared_ptr<RooRealVar> de_p1 = make_shared<RooRealVar>(("de_p1_"+category).c_str(),"Double Exponential 1 width 1,1",0.0,20.0); vars.push_back(de_p1); doubleexp_argset.add(*de_p1);
    //shared_ptr<RooRealVar> de_o0 = make_shared<RooRealVar>(("de_o0_"+category).c_str(),"Double Exponential 1 sigmoid offset 0",90.0,120.0); vars.push_back(de_o0); doubleexp_argset.add(*de_o0);
    //shared_ptr<RooRealVar> de_o1 = make_shared<RooRealVar>(("de_o1_"+category).c_str(),"Double Exponential 1 sigmoid offset 1",90.0,120.0); vars.push_back(de_o1); doubleexp_argset.add(*de_o1);
    //shared_ptr<RooRealVar> de_w0 = make_shared<RooRealVar>(("de_w0_"+category).c_str(),"Double Exponential 1 sigmoid width 0",0.001,20.0); vars.push_back(de_w0); doubleexp_argset.add(*de_w0);
    //shared_ptr<RooRealVar> de_w1 = make_shared<RooRealVar>(("de_w1_"+category).c_str(),"Double Exponential 1 sigmoid width 1",0.001,20.0); vars.push_back(de_w1); doubleexp_argset.add(*de_w1);
    //shared_ptr<RooRealVar> de_pd = make_shared<RooRealVar>(("de_pd_"+category).c_str(),"Double Exponential 1 pedestal",-1.0,1.0); vars.push_back(de_pd); doubleexp_argset.add(*de_pd);
    //de_c0->setVal(0.1);
    //de_c1->setVal(0.8);
    //de_p0->setVal(4.0);
    //de_p1->setVal(5.0);
    //de_o0->setVal(105);
    //de_o1->setVal(112);
    //de_w0->setVal(0.8);
    //de_w1->setVal(0.8);
    //de_pd->setVal(0.0);
    //RooGenericPdf* pdf_bkg_cat_doubleexp1 = new RooGenericPdf(("pdf_background_"+category+"_doubleexp1").c_str(),"bkg_pdf",
    //    "(@1*exp(-1.0*@3*(@0-105)/60))/(1.0+exp(-1.0*@7*(@0-@5)))+(@2*exp(-1.0*@4*(165-@0)/60))/(1.0+exp(1.0*@8*(@0-@6)))+@9",
    //    doubleexp_argset);
    
    //RooArgSet laurent3_argset;
    //laurent3_argset.add(*pdf_mllg);
    //shared_ptr<RooRealVar> l2_c0 = make_shared<RooRealVar>(("l2_c0_"+category).c_str(),"Laurent 5 coefficient 0",-1.0,1.0); vars.push_back(l2_c0); laurent3_argset.add(*l2_c0);
    //shared_ptr<RooRealVar> l2_c1 = make_shared<RooRealVar>(("l2_c1_"+category).c_str(),"Laurent 5 coefficient 1",-1.0,1.0); vars.push_back(l2_c1); laurent3_argset.add(*l2_c1);
    //shared_ptr<RooRealVar> l2_c2 = make_shared<RooRealVar>(("l2_c2_"+category).c_str(),"Laurent 5 coefficient 2",-1.0,1.0); vars.push_back(l2_c2); laurent3_argset.add(*l2_c2);
    //shared_ptr<RooRealVar> l2_c3 = make_shared<RooRealVar>(("l2_c3_"+category).c_str(),"Laurent 5 coefficient 3",-1.0,1.0); vars.push_back(l2_c3); laurent3_argset.add(*l2_c3);
    //shared_ptr<RooRealVar> l2_o1 = make_shared<RooRealVar>(("l2_o1_"+category).c_str(),"Laurent 5 sigmoid offset",90.0,120.0); vars.push_back(l2_o1); laurent3_argset.add(*l2_o1);
    //shared_ptr<RooRealVar> l2_w1 = make_shared<RooRealVar>(("l2_w1_"+category).c_str(),"Laurent 5 sigmoid width",0.005,30.0); vars.push_back(l2_w1); laurent3_argset.add(*l2_w1);
    //l2_c0->setVal(0.2);
    //l2_c1->setVal(0.2);
    //l2_c2->setVal(0.2);
    //l2_c3->setVal(0.8);
    //l2_o1->setVal(105);
    //l2_w1->setVal(0.1);
    //RooGenericPdf* pdf_bkg_cat_laurent3 = new RooGenericPdf(("pdf_background_"+category+"_laurent3").c_str(),"bkg_pdf",
    //    "(@1/pow((@0-75)/40,3)+@2/pow((@0-75)/40,4)+@3/pow((@0-75)/40,5)+@4/pow((@0-75)/40,6))/(1.0+exp(-1.0*@6*(@0-@5)))",
    //    laurent3_argset);
    
    //RooArgSet exp5_argset;
    //exp5_argset.add(*pdf_mllg);
    //shared_ptr<RooRealVar> e5_c0 = make_shared<RooRealVar>(("e5_c0_"+category).c_str(),"Exponential 5 coefficient 0",-1.0,1.0); vars.push_back(e5_c0); exp5_argset.add(*e5_c0);
    //shared_ptr<RooRealVar> e5_c1 = make_shared<RooRealVar>(("e5_c1_"+category).c_str(),"Exponential 5 coefficient 1",-1.0,1.0); vars.push_back(e5_c1); exp5_argset.add(*e5_c1);
    //shared_ptr<RooRealVar> e5_c2 = make_shared<RooRealVar>(("e5_c2_"+category).c_str(),"Exponential 5 coefficient 2",-1.0,1.0); vars.push_back(e5_c2); exp5_argset.add(*e5_c2);
    //shared_ptr<RooRealVar> e5_c3 = make_shared<RooRealVar>(("e5_c3_"+category).c_str(),"Exponential 5 coefficient 3",-1.0,1.0); vars.push_back(e5_c3); exp5_argset.add(*e5_c3);
    //shared_ptr<RooRealVar> e5_c4 = make_shared<RooRealVar>(("e5_c4_"+category).c_str(),"Exponential 5 coefficient 4",-1.0,1.0); vars.push_back(e5_c4); exp5_argset.add(*e5_c4);
    //shared_ptr<RooRealVar> e5_p0 = make_shared<RooRealVar>(("e5_p0_"+category).c_str(),"Exponential 5 width 0",-20.0,0.0); vars.push_back(e5_p0); exp5_argset.add(*e5_p0);
    //shared_ptr<RooRealVar> e5_p1 = make_shared<RooRealVar>(("e5_p1_"+category).c_str(),"Exponential 5 width 1",-20.0,0.0); vars.push_back(e5_p1); exp5_argset.add(*e5_p1);
    //shared_ptr<RooRealVar> e5_p2 = make_shared<RooRealVar>(("e5_p2_"+category).c_str(),"Exponential 5 width 2",-20.0,0.0); vars.push_back(e5_p2); exp5_argset.add(*e5_p2);
    //shared_ptr<RooRealVar> e5_p3 = make_shared<RooRealVar>(("e5_p3_"+category).c_str(),"Exponential 5 width 3",-20.0,0.0); vars.push_back(e5_p3); exp5_argset.add(*e5_p3);
    //shared_ptr<RooRealVar> e5_p4 = make_shared<RooRealVar>(("e5_p4_"+category).c_str(),"Exponential 5 width 4",-20.0,0.0); vars.push_back(e5_p4); exp5_argset.add(*e5_p4);
    //e5_c0->setVal(0.5);
    //e5_c1->setVal(0.5);
    //e5_c2->setVal(0.5);
    //e5_c3->setVal(0.5);
    //e5_c4->setVal(0.5);
    //e5_p0->setVal(0.5);
    //e5_p1->setVal(1.0);
    //e5_p2->setVal(1.5);
    //e5_p3->setVal(2.0);
    //e5_p4->setVal(2.5);
    //RooGenericPdf* pdf_bkg_cat_exp5 = new RooGenericPdf(("pdf_background_"+category+"_exp5").c_str(),"bkg_pdf",
    //    "@1*exp(@6*(@0-120)/35)+@2*exp(@7*(@0-120)/35)+@3*exp(@8*(@0-120)/35)+@4*exp(@9*(@0-120)/35)+@5*exp(@10*(@0-120)/35)",
    //    exp5_argset);

    //RooArgSet exp7_argset;
    //exp7_argset.add(*pdf_mllg);
    //shared_ptr<RooRealVar> e7_c0 = make_shared<RooRealVar>(("e7_c0_"+category).c_str(),"Exponential 7 coefficient 0",-1.0,1.0); vars.push_back(e7_c0); exp7_argset.add(*e7_c0);
    //shared_ptr<RooRealVar> e7_c1 = make_shared<RooRealVar>(("e7_c1_"+category).c_str(),"Exponential 7 coefficient 1",-1.0,1.0); vars.push_back(e7_c1); exp7_argset.add(*e7_c1);
    //shared_ptr<RooRealVar> e7_c2 = make_shared<RooRealVar>(("e7_c2_"+category).c_str(),"Exponential 7 coefficient 2",-1.0,1.0); vars.push_back(e7_c2); exp7_argset.add(*e7_c2);
    //shared_ptr<RooRealVar> e7_c3 = make_shared<RooRealVar>(("e7_c3_"+category).c_str(),"Exponential 7 coefficient 3",-1.0,1.0); vars.push_back(e7_c3); exp7_argset.add(*e7_c3);
    //shared_ptr<RooRealVar> e7_c4 = make_shared<RooRealVar>(("e7_c4_"+category).c_str(),"Exponential 7 coefficient 4",-1.0,1.0); vars.push_back(e7_c4); exp7_argset.add(*e7_c4);
    //shared_ptr<RooRealVar> e7_c5 = make_shared<RooRealVar>(("e7_c5_"+category).c_str(),"Exponential 7 coefficient 5",-1.0,1.0); vars.push_back(e7_c5); exp7_argset.add(*e7_c5);
    //shared_ptr<RooRealVar> e7_c6 = make_shared<RooRealVar>(("e7_c6_"+category).c_str(),"Exponential 7 coefficient 6",-1.0,1.0); vars.push_back(e7_c6); exp7_argset.add(*e7_c6);
    //shared_ptr<RooRealVar> e7_p0 = make_shared<RooRealVar>(("e7_p0_"+category).c_str(),"Exponential 7 width 0",-20.0,0.0); vars.push_back(e7_p0); exp7_argset.add(*e7_p0);
    //shared_ptr<RooRealVar> e7_p1 = make_shared<RooRealVar>(("e7_p1_"+category).c_str(),"Exponential 7 width 1",-20.0,0.0); vars.push_back(e7_p1); exp7_argset.add(*e7_p1);
    //shared_ptr<RooRealVar> e7_p2 = make_shared<RooRealVar>(("e7_p2_"+category).c_str(),"Exponential 7 width 2",-20.0,0.0); vars.push_back(e7_p2); exp7_argset.add(*e7_p2);
    //shared_ptr<RooRealVar> e7_p3 = make_shared<RooRealVar>(("e7_p3_"+category).c_str(),"Exponential 7 width 3",-20.0,0.0); vars.push_back(e7_p3); exp7_argset.add(*e7_p3);
    //shared_ptr<RooRealVar> e7_p4 = make_shared<RooRealVar>(("e7_p4_"+category).c_str(),"Exponential 7 width 4",-20.0,0.0); vars.push_back(e7_p4); exp7_argset.add(*e7_p4);
    //shared_ptr<RooRealVar> e7_p5 = make_shared<RooRealVar>(("e7_p5_"+category).c_str(),"Exponential 7 width 5",-20.0,0.0); vars.push_back(e7_p5); exp7_argset.add(*e7_p5);
    //shared_ptr<RooRealVar> e7_p6 = make_shared<RooRealVar>(("e7_p6_"+category).c_str(),"Exponential 7 width 6",-20.0,0.0); vars.push_back(e7_p6); exp7_argset.add(*e7_p6);
    //shared_ptr<RooRealVar> e7_so = make_shared<RooRealVar>(("e7_so_"+category).c_str(),"Exponential 7 sigmoid offset",90.0,120.0); vars.push_back(e7_so); exp7_argset.add(*e7_so);
    //shared_ptr<RooRealVar> e7_sw = make_shared<RooRealVar>(("e7_sw_"+category).c_str(),"Exponential 7 sigmoid width",0.005,20.0); vars.push_back(e7_sw); exp7_argset.add(*e7_sw);
    //e7_c0->setVal(0.5);
    //e7_c1->setVal(0.5);
    //e7_c2->setVal(0.5);
    //e7_c3->setVal(0.5);
    //e7_c4->setVal(0.5);
    //e7_c5->setVal(0.5);
    //e7_c6->setVal(0.5);
    //e7_p0->setVal(0.5);
    //e7_p1->setVal(1.0);
    //e7_p2->setVal(1.5);
    //e7_p3->setVal(2.0);
    //e7_p4->setVal(2.5);
    //e7_p5->setVal(3.0);
    //e7_p6->setVal(3.5);
    //e7_so->setVal(105);
    //e7_sw->setVal(0.1);
    //RooGenericPdf* pdf_bkg_cat_exp7 = new RooGenericPdf(("pdf_background_"+category+"_exp7").c_str(),"bkg_pdf",
    //    "(@1*exp(@8*(@0-105)/55)+@2*exp(@9*(@0-105)/55)+@3*exp(@10*(@0-105)/55)+@4*exp(@11*(@0-105)/55)+@5*exp(@12*(@0-105)/55)+@6*exp(@13*(@0-105)/55)+@7*exp(@14*(@0-105)/55))/(1.0+exp(-1.0*@16*(@0-@15)))",
    //    exp7_argset);

    //RooArgSet exp5_argset;
    //exp5_argset.add(*pdf_mllg);
    //shared_ptr<RooRealVar> e5_c0 = make_shared<RooRealVar>(("e5_c0_"+category).c_str(),"Exponential 5 coefficient 0",-1.0,1.0); vars.push_back(e5_c0); exp5_argset.add(*e5_c0);
    //shared_ptr<RooRealVar> e5_c1 = make_shared<RooRealVar>(("e5_c1_"+category).c_str(),"Exponential 5 coefficient 1",-1.0,1.0); vars.push_back(e5_c1); exp5_argset.add(*e5_c1);
    //shared_ptr<RooRealVar> e5_c2 = make_shared<RooRealVar>(("e5_c2_"+category).c_str(),"Exponential 5 coefficient 2",-1.0,1.0); vars.push_back(e5_c2); exp5_argset.add(*e5_c2);
    //shared_ptr<RooRealVar> e5_c3 = make_shared<RooRealVar>(("e5_c3_"+category).c_str(),"Exponential 5 coefficient 3",-1.0,1.0); vars.push_back(e5_c3); exp5_argset.add(*e5_c3);
    //shared_ptr<RooRealVar> e5_c4 = make_shared<RooRealVar>(("e5_c4_"+category).c_str(),"Exponential 5 coefficient 4",-1.0,1.0); vars.push_back(e5_c4); exp5_argset.add(*e5_c4);
    //shared_ptr<RooRealVar> e5_p0 = make_shared<RooRealVar>(("e5_p0_"+category).c_str(),"Exponential 5 width 0",-20.0,0.0); vars.push_back(e5_p0); exp5_argset.add(*e5_p0);
    //shared_ptr<RooRealVar> e5_p1 = make_shared<RooRealVar>(("e5_p1_"+category).c_str(),"Exponential 5 width 1",-20.0,0.0); vars.push_back(e5_p1); exp5_argset.add(*e5_p1);
    //shared_ptr<RooRealVar> e5_p2 = make_shared<RooRealVar>(("e5_p2_"+category).c_str(),"Exponential 5 width 2",-20.0,0.0); vars.push_back(e5_p2); exp5_argset.add(*e5_p2);
    //shared_ptr<RooRealVar> e5_p3 = make_shared<RooRealVar>(("e5_p3_"+category).c_str(),"Exponential 5 width 3",-20.0,0.0); vars.push_back(e5_p3); exp5_argset.add(*e5_p3);
    //shared_ptr<RooRealVar> e5_p4 = make_shared<RooRealVar>(("e5_p4_"+category).c_str(),"Exponential 5 width 4",-20.0,0.0); vars.push_back(e5_p4); exp5_argset.add(*e5_p4);
    //shared_ptr<RooRealVar> e5_so = make_shared<RooRealVar>(("e5_so_"+category).c_str(),"Exponential 5 sigmoid offset",90.0,120.0); vars.push_back(e5_so); exp5_argset.add(*e5_so);
    //shared_ptr<RooRealVar> e5_sw = make_shared<RooRealVar>(("e5_sw_"+category).c_str(),"Exponential 5 sigmoid width",0.005,20.0); vars.push_back(e5_sw); exp5_argset.add(*e5_sw);
    //e5_c0->setVal(0.5);
    //e5_c1->setVal(0.5);
    //e5_c2->setVal(0.5);
    //e5_c3->setVal(0.5);
    //e5_c4->setVal(0.5);
    //e5_p0->setVal(0.5);
    //e5_p1->setVal(1.0);
    //e5_p2->setVal(1.5);
    //e5_p3->setVal(2.0);
    //e5_p4->setVal(2.5);
    //e5_so->setVal(104);
    //e5_sw->setVal(0.25);
    //RooGenericPdf* pdf_bkg_cat_exp5 = new RooGenericPdf(("pdf_background_"+category+"_exp5").c_str(),"bkg_pdf",
    //    "(@1*exp(@6*(@0-105)/55)+@2*exp(@7*(@0-105)/55)+@3*exp(@8*(@0-105)/55)+@4*exp(@9*(@0-105)/55)+@5*exp(@10*(@0-105)/55))/(1.0+exp(-1*@12*(@0-@11)))",
    //    exp5_argset);

    //RooArgSet bern6_argset;
    //bern6_argset.add(*pdf_mllg);
    //shared_ptr<RooRealVar> b6_c0 = make_shared<RooRealVar>(("b6_c0_"+category).c_str(),"bern 6 coefficient 0",0.0,1.0); vars.push_back(b6_c0); bern6_argset.add(*b6_c0);
    //shared_ptr<RooRealVar> b6_c1 = make_shared<RooRealVar>(("b6_c1_"+category).c_str(),"bern 6 coefficient 1",0.0,1.0); vars.push_back(b6_c1); bern6_argset.add(*b6_c1);
    //shared_ptr<RooRealVar> b6_c2 = make_shared<RooRealVar>(("b6_c2_"+category).c_str(),"bern 6 coefficient 2",0.0,1.0); vars.push_back(b6_c2); bern6_argset.add(*b6_c2);
    //shared_ptr<RooRealVar> b6_c3 = make_shared<RooRealVar>(("b6_c3_"+category).c_str(),"bern 6 coefficient 3",0.0,1.0); vars.push_back(b6_c3); bern6_argset.add(*b6_c3);
    //shared_ptr<RooRealVar> b6_c4 = make_shared<RooRealVar>(("b6_c4_"+category).c_str(),"bern 6 coefficient 4",0.0,1.0); vars.push_back(b6_c4); bern6_argset.add(*b6_c4);
    //shared_ptr<RooRealVar> b6_c5 = make_shared<RooRealVar>(("b6_c5_"+category).c_str(),"bern 6 coefficient 5",0.0,1.0); vars.push_back(b6_c5); bern6_argset.add(*b6_c5);
    //shared_ptr<RooRealVar> b6_c6 = make_shared<RooRealVar>(("b6_c6_"+category).c_str(),"bern 6 coefficient 6",0.0,1.0); vars.push_back(b6_c6); bern6_argset.add(*b6_c6);
    //shared_ptr<RooRealVar> b6_so = make_shared<RooRealVar>(("b6_so_"+category).c_str(),"bern 6 sigmoid offset",90.0,115.0); vars.push_back(b6_so); bern6_argset.add(*b6_so);
    //shared_ptr<RooRealVar> b6_sw = make_shared<RooRealVar>(("b6_sw_"+category).c_str(),"bern 6 sigmoid width",0.005,30.0); vars.push_back(b6_sw); bern6_argset.add(*b6_sw);
    //b6_c0->setVal(1.0);
    //b6_c1->setVal(0.1);
    //b6_c2->setVal(0.1);
    //b6_c3->setVal(0.1);
    //b6_c4->setVal(0.1);
    //b6_c5->setVal(0.1);
    //b6_c6->setVal(0.1);
    //b6_so->setVal(105);
    //b6_sw->setVal(0.03);
    //RooGenericPdf* pdf_bkg_cat_bern6 = new RooGenericPdf(("pdf_background_"+category+"_bern6").c_str(),"bkg_pdf",
    //    "(@1*pow((160-@0)/60,6)+@2*(6*(@0-100)/60*pow((160-@0)/60,5))+@3*(15*pow((@0-100)/60,2)*pow((160-@0)/60,4))+@4*(20*pow((@0-100)/60,3)*pow((160-@0)/60,3))+@5*(15*pow((@0-100)/60,4)*pow((160-@0)/60,2))+@6*(5*pow((@0-100)/60,5)*(160-@0)/60)+@7*pow((@0-100)/60,6))/(1.0+exp(-1.0*@9*(@0-@8)))",
    //    bern6_argset);
    
    //RooArgSet bern6_argset;
    //bern6_argset.add(*pdf_mllg);
    //shared_ptr<RooRealVar> b6_pd = make_shared<RooRealVar>(("b6_pd_"+category).c_str(),"bern 6 pedestal",-1.0,1.0); vars.push_back(b6_pd); bern6_argset.add(*b6_pd);
    //shared_ptr<RooRealVar> b6_c0 = make_shared<RooRealVar>(("b6_c0_"+category).c_str(),"bern 6 coefficient 0",-1.0,1.0); vars.push_back(b6_c0); bern6_argset.add(*b6_c0);
    //shared_ptr<RooRealVar> b6_c1 = make_shared<RooRealVar>(("b6_c1_"+category).c_str(),"bern 6 coefficient 1",-1.0,1.0); vars.push_back(b6_c1); bern6_argset.add(*b6_c1);
    //shared_ptr<RooRealVar> b6_c2 = make_shared<RooRealVar>(("b6_c2_"+category).c_str(),"bern 6 coefficient 2",-1.0,1.0); vars.push_back(b6_c2); bern6_argset.add(*b6_c2);
    //shared_ptr<RooRealVar> b6_c3 = make_shared<RooRealVar>(("b6_c3_"+category).c_str(),"bern 6 coefficient 3",-1.0,1.0); vars.push_back(b6_c3); bern6_argset.add(*b6_c3);
    //shared_ptr<RooRealVar> b6_c4 = make_shared<RooRealVar>(("b6_c4_"+category).c_str(),"bern 6 coefficient 4",-1.0,1.0); vars.push_back(b6_c4); bern6_argset.add(*b6_c4);
    //shared_ptr<RooRealVar> b6_c5 = make_shared<RooRealVar>(("b6_c5_"+category).c_str(),"bern 6 coefficient 5",-1.0,1.0); vars.push_back(b6_c5); bern6_argset.add(*b6_c5);
    //shared_ptr<RooRealVar> b6_c6 = make_shared<RooRealVar>(("b6_c6_"+category).c_str(),"bern 6 coefficient 6",-1.0,1.0); vars.push_back(b6_c6); bern6_argset.add(*b6_c6);
    //shared_ptr<RooRealVar> b6_12_c1 = make_shared<RooRealVar>(("b6_12_c1_"+category).c_str(),"bern 12 coefficient 1",-1.0,1.0); vars.push_back(b6_12_c1); bern6_argset.add(*b6_12_c1);
    //shared_ptr<RooRealVar> b6_30_c2 = make_shared<RooRealVar>(("b6_30_c2_"+category).c_str(),"bern 30 coefficient 2",-1.0,1.0); vars.push_back(b6_30_c2); bern6_argset.add(*b6_30_c2);
    //b6_pd->setVal(0.0);
    //b6_c0->setVal(1.0);
    //b6_c1->setVal(0.1);
    //b6_c2->setVal(0.1);
    //b6_c3->setVal(0.1);
    //b6_c4->setVal(0.1);
    //b6_c5->setVal(0.1);
    //b6_c6->setVal(0.1);
    //b6_12_c1->setVal(0.1);
    //b6_30_c2->setVal(0.1);
    //RooGenericPdf* pdf_bkg_cat_bern6 = new RooGenericPdf(("pdf_background_"+category+"_bern6").c_str(),"bkg_pdf",
    //    "@1+@2*pow((165-@0)/60,6)+@3*(6*(@0-105)/60*pow((165-@0)/60,5))+@4*(15*pow((@0-105)/60,2)*pow((165-@0)/60,4))+@5*(20*pow((@0-105)/60,3)*pow((165-@0)/60,3))+@6*(15*pow((@0-105)/60,4)*pow((165-@0)/60,2))+@7*(5*pow((@0-105)/60,5)*(165-@0)/60)+@8*pow((@0-105)/60,6)+@9*12*((@0-105)/60)*pow((165-@0)/60,12)+@10*435*pow((@0-105)/60,2)*pow((165-@0)/60,28)",
    //    bern6_argset);

    //RooArgSet bern8_argset;
    //shared_ptr<RooRealVar> b8_c0 = make_shared<RooRealVar>(("b8_c0_"+category).c_str(),"bern 8 coefficient 0",-1.0,1.0); vars.push_back(b8_c0); bern8_argset.add(*b8_c0);
    //shared_ptr<RooRealVar> b8_c1 = make_shared<RooRealVar>(("b8_c1_"+category).c_str(),"bern 8 coefficient 1",-1.0,1.0); vars.push_back(b8_c1); bern8_argset.add(*b8_c1);
    //shared_ptr<RooRealVar> b8_c2 = make_shared<RooRealVar>(("b8_c2_"+category).c_str(),"bern 8 coefficient 2",-1.0,1.0); vars.push_back(b8_c2); bern8_argset.add(*b8_c2);
    //shared_ptr<RooRealVar> b8_c3 = make_shared<RooRealVar>(("b8_c3_"+category).c_str(),"bern 8 coefficient 3",-1.0,1.0); vars.push_back(b8_c3); bern8_argset.add(*b8_c3);
    //shared_ptr<RooRealVar> b8_c4 = make_shared<RooRealVar>(("b8_c4_"+category).c_str(),"bern 8 coefficient 4",-1.0,1.0); vars.push_back(b8_c4); bern8_argset.add(*b8_c4);
    //shared_ptr<RooRealVar> b8_c5 = make_shared<RooRealVar>(("b8_c5_"+category).c_str(),"bern 8 coefficient 5",-1.0,1.0); vars.push_back(b8_c5); bern8_argset.add(*b8_c5);
    //shared_ptr<RooRealVar> b8_c6 = make_shared<RooRealVar>(("b8_c6_"+category).c_str(),"bern 8 coefficient 6",-1.0,1.0); vars.push_back(b8_c6); bern8_argset.add(*b8_c6);
    //shared_ptr<RooRealVar> b8_c7 = make_shared<RooRealVar>(("b8_c7_"+category).c_str(),"bern 8 coefficient 7",-1.0,1.0); vars.push_back(b8_c7); bern8_argset.add(*b8_c7);
    //shared_ptr<RooRealVar> b8_c8 = make_shared<RooRealVar>(("b8_c8_"+category).c_str(),"bern 8 coefficient 8",-1.0,1.0); vars.push_back(b8_c8); bern8_argset.add(*b8_c8);
    //b8_c0->setVal(0.8);
    //b8_c1->setVal(0.2);
    //b8_c2->setVal(0.1);
    //b8_c3->setVal(0.1);
    //b8_c4->setVal(0.1);
    //b8_c5->setVal(0.1);
    //b8_c6->setVal(0.1);
    //b8_c7->setVal(0.1);
    //b8_c8->setVal(0.1);
    //RooBernstein* pdf_bkg_cat_bern8 = new RooBernstein(("pdf_background_"+category+"_bern8").c_str(),"bkg_pdf",*pdf_mllg,bern8_argset);

    //RooArgSet bern3_argset;
    //bern3_argset.add(*pdf_mllg);
    //shared_ptr<RooRealVar> b3_c0 = make_shared<RooRealVar>(("b3_c0_"+category).c_str(),"bern 3 coefficient 0",0.0,1.0); vars.push_back(b3_c0); bern3_argset.add(*b3_c0);
    //shared_ptr<RooRealVar> b3_c1 = make_shared<RooRealVar>(("b3_c1_"+category).c_str(),"bern 3 coefficient 1",0.0,1.0); vars.push_back(b3_c1); bern3_argset.add(*b3_c1);
    //shared_ptr<RooRealVar> b3_c2 = make_shared<RooRealVar>(("b3_c2_"+category).c_str(),"bern 3 coefficient 2",0.0,1.0); vars.push_back(b3_c2); bern3_argset.add(*b3_c2);
    //shared_ptr<RooRealVar> b3_c3 = make_shared<RooRealVar>(("b3_c3_"+category).c_str(),"bern 3 coefficient 3",0.0,1.0); vars.push_back(b3_c3); bern3_argset.add(*b3_c3);
    //shared_ptr<RooRealVar> b3_o1 = make_shared<RooRealVar>(("b3_o1_"+category).c_str(),"bern 3 sigmoid 1 offset",100.0,115.0); vars.push_back(b3_o1); bern3_argset.add(*b3_o1);
    //shared_ptr<RooRealVar> b3_w1 = make_shared<RooRealVar>(("b3_w1_"+category).c_str(),"bern 3 sigmoid 1 width",0.05,20.0); vars.push_back(b3_w1); bern3_argset.add(*b3_w1);
    //shared_ptr<RooRealVar> b3_i1 = make_shared<RooRealVar>(("b3_i1_"+category).c_str(),"bern 3 sigmoid 1 pedestal",0.0,1.0); vars.push_back(b3_i1); bern3_argset.add(*b3_i1);
    //shared_ptr<RooRealVar> b3_o2 = make_shared<RooRealVar>(("b3_o2_"+category).c_str(),"bern 3 sigmoid 2 offset",100.0,115.0); vars.push_back(b3_o2); bern3_argset.add(*b3_o2);
    //shared_ptr<RooRealVar> b3_w2 = make_shared<RooRealVar>(("b3_w2_"+category).c_str(),"bern 3 sigmoid 2 width",0.05,20.0); vars.push_back(b3_w2); bern3_argset.add(*b3_w2);
    //shared_ptr<RooRealVar> b3_i2 = make_shared<RooRealVar>(("b3_i2_"+category).c_str(),"bern 3 sigmoid 2 pedestal",0.0,1.0); vars.push_back(b3_i2); bern3_argset.add(*b3_i2);
    //b3_c0->setVal(1.0);
    //b3_c1->setVal(0.1);
    //b3_c2->setVal(0.1);
    //b3_c3->setVal(0.1);
    //b3_o1->setVal(105);
    //b3_w1->setVal(0.03);
    //b3_i1->setVal(0.0);
    //b3_o2->setVal(110);
    //b3_w2->setVal(0.03);
    //b3_i2->setVal(0.1);
    //RooGenericPdf* pdf_bkg_cat_bern3 = new RooGenericPdf(("pdf_background_"+category+"_bern3").c_str(),"bkg_pdf",
    //    "(@1*pow(1-(@0-100)/60,3)+@2*(3*(@0-100)/60*pow(1-(@0-100)/60,2))+@3*(3*pow((@0-100)/60,2)*(1-(@0-100)/60))+@4*pow((@0-100)/60,3))*(@7+(1.0-@7)/(1.0+exp(-1.0*@6*(@0-@5))))*(@10+(1.0-@10)/(1.0+exp(-1.0*@9*(@0-@8))))",
    //    bern3_argset);

    //RooArgSet bern3_argset;
    //bern3_argset.add(*pdf_mllg);
    //shared_ptr<RooRealVar> b3_c0 = make_shared<RooRealVar>(("b3_c0_"+category).c_str(),"bern 3 coefficient 0",-1.0,1.0); vars.push_back(b3_c0); bern3_argset.add(*b3_c0);
    //shared_ptr<RooRealVar> b3_c1 = make_shared<RooRealVar>(("b3_c1_"+category).c_str(),"bern 3 coefficient 1",-1.0,1.0); vars.push_back(b3_c1); bern3_argset.add(*b3_c1);
    //shared_ptr<RooRealVar> b3_c2 = make_shared<RooRealVar>(("b3_c2_"+category).c_str(),"bern 3 coefficient 2",-1.0,1.0); vars.push_back(b3_c2); bern3_argset.add(*b3_c2);
    //shared_ptr<RooRealVar> b3_c3 = make_shared<RooRealVar>(("b3_c3_"+category).c_str(),"bern 3 coefficient 3",-1.0,1.0); vars.push_back(b3_c3); bern3_argset.add(*b3_c3);
    //shared_ptr<RooRealVar> b3_o1 = make_shared<RooRealVar>(("b3_o1_"+category).c_str(),"bern 3 sigmoid offset",100.0,115.0); vars.push_back(b3_o1); bern3_argset.add(*b3_o1);
    //shared_ptr<RooRealVar> b3_w1 = make_shared<RooRealVar>(("b3_w1_"+category).c_str(),"bern 3 sigmoid linear width",0.05,20.0); vars.push_back(b3_w1); bern3_argset.add(*b3_w1);
    //shared_ptr<RooRealVar> b3_w2 = make_shared<RooRealVar>(("b3_w2_"+category).c_str(),"bern 3 sigmoid sqrt width",0.0,20.0); vars.push_back(b3_w2); bern3_argset.add(*b3_w2);
    //b3_c0->setVal(1.0);
    //b3_c1->setVal(0.1);
    //b3_c2->setVal(0.1);
    //b3_c3->setVal(0.1);
    //b3_o1->setVal(105);
    //b3_w1->setVal(0.01);
    //b3_w2->setVal(0.01);
    //RooGenericPdf* pdf_bkg_cat_bern3 = new RooGenericPdf(("pdf_background_"+category+"_bern3").c_str(),"bkg_pdf",
    //    "(@1*pow(1-(@0-100)/60,3)+@2*(3*(@0-100)/60*pow(1-(@0-100)/60,2))+@3*(3*pow((@0-100)/60,2)*(1-(@0-100)/60))+@4*pow((@0-100)/60,3))/(1.0+exp(-1.0*@6*(@0-@5)-1.0*@7*TMath::Sign(sqrt(abs(@0-@5)),@0-@5)))",
    //    bern3_argset);

    //shared_ptr<RooRealVar> e1_c0 = make_shared<RooRealVar>(("c0_"+category).c_str(),"exp 1 exponential coefficient",0.0,10.0); vars.push_back(e1_c0);
    //shared_ptr<RooRealVar> e1_so = make_shared<RooRealVar>(("c1_"+category).c_str(),"exp 1 sigmoid offset",100.0,115.0); vars.push_back(e1_so);
    //shared_ptr<RooRealVar> e1_sw = make_shared<RooRealVar>(("c2_"+category).c_str(),"exp 1 sigmoid width",0.05,20.0); vars.push_back(e1_sw);
    //RooGenericPdf* pdf_bkg_cat_exp1 = new RooGenericPdf(("pdf_background_"+category+"_exp1").c_str(),"bkg_pdf",
    //    "exp(-1.0*@1*(@0-100.0)/60.0)/(1.0+exp(-1.0*@3*(@0-@2)))",
    //    RooArgSet(*pdf_mllg,*e1_c0,*e1_so,*e1_sw));

    RooGenericPdf* pdf_sig_cat = new RooGenericPdf(("pdf_htozg_"+category).c_str(),"sig_pdf",
        "TMath::Gaus(@0,125.0,1.5)",
        RooArgSet(*pdf_mllg));

    //background_pdfs.push_back(vector<RooAbsPdf*>{pdf_bkg_cat_bern3,pdf_bkg_cat_exp1});
    background_pdfs.push_back(pdf_bkg_cat_bern6);
    signal_pdfs.push_back(pdf_sig_cat);
  }
  ////std::cout << background_pdfs[0][0]->GetName() << std::endl;
  ////std::cout << background_pdfs[0][1]->GetName() << std::endl;
  //std::cout << background_pdfs[0]->GetName() << std::endl;
  //std::cout << signal_pdfs[0]->GetName() << std::endl;
  ////std::cout << background_pdfs[1][0]->GetName() << std::endl;
  ////std::cout << background_pdfs[1][1]->GetName() << std::endl;
  //std::cout << background_pdfs[1]->GetName() << std::endl;
  //std::cout << signal_pdfs[1]->GetName() << std::endl;

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

  pm.Push<Datacard>("test_datacard", channels, systematics, processes, weight,
      Axis(50, 110.0, 160.0, mllg, "m_{ll#gamma} [GeV]", {}))
      .AddParametricProcess("background",background_pdfs)
      .MakeProcessParametric("htozgamma",signal_pdfs);

  pm.MakePlots(1.0);

  return 0;
}
