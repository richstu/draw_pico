#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <bitset>
#include <unistd.h>
#include <getopt.h>
#include "TError.h"
#include "TColor.h"
#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"
#include "TMVA/Reader.h"
#include "TMVA/Configurable.h"
#include "TLorentzVector.h"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/utilities.hpp"
#include "zgamma/zg_utilities.hpp"
using namespace std;
using namespace PlotOptTypes;
using namespace ZgUtilities;

namespace {

  float get_l1_rapidity(const Baby &b) {
    if (b.ll_lepid()->at(0)==11) {
      return (b.el_pt()->at(b.ll_i1()->at(0)) > b.el_pt()->at(b.ll_i2()->at(0))) ? b.el_eta()->at(b.ll_i1()->at(0)) : b.el_eta()->at(b.ll_i2()->at(0));
    }
    return (b.mu_pt()->at(b.ll_i1()->at(0)) > b.mu_pt()->at(b.ll_i2()->at(0))) ? b.mu_eta()->at(b.ll_i1()->at(0)) : b.mu_eta()->at(b.ll_i2()->at(0));
  }

  float get_l2_rapidity(const Baby &b) {
    if (b.ll_lepid()->at(0)==11) {
      return (b.el_pt()->at(b.ll_i1()->at(0)) < b.el_pt()->at(b.ll_i2()->at(0))) ? b.el_eta()->at(b.ll_i1()->at(0)) : b.el_eta()->at(b.ll_i2()->at(0));
    }
    return (b.mu_pt()->at(b.ll_i1()->at(0)) < b.mu_pt()->at(b.ll_i2()->at(0))) ? b.mu_eta()->at(b.ll_i1()->at(0)) : b.mu_eta()->at(b.ll_i2()->at(0));
  }

  //All global variables that help, Find best place to delete reader when done
  Float_t ct = 0;
  Float_t cT = 0;
  Float_t pdmax = 0;
  Float_t gphi = 0;
  Float_t yrV = 0;
  Float_t yres = 0;
  Float_t yrap = 0;
  Float_t yidmva = 0;
  Float_t pdmin = 0;
  Float_t l1eta = 0;
  Float_t l2eta = 0;
  TMVA::Reader *reader =  new TMVA::Reader();
  
  //Long64_t  evNumber = -1;
  double currentBDTScore;
  
  
  //Function that returns BDT score
  double mva_score(const Baby &b){//, TMVA::Reader *reader){ 
  
    //cout << "top of mva_score" << endl;
  
    //Checks if the value has been evaluated during that event before <--- New as of July 15th, 2022
    //if(evNumber == b.event() ){
    //  return currentBDTScore;
    //} else {
    //  evNumber = b.event();
    //}
  
    //Redundant check to ensure no errors with checking out of the bounds of an array
    if( ( (b.el_pt() -> size() < 2) && (b.mu_pt() -> size() < 2) ) || (b.photon_pt() -> size() < 1) ){ return -2; }
  
  
    //This sets lepton eta (could use lep_eta branch?)
    if( b.ll_lepid() -> at(0) == 11){
      l1eta = static_cast<Float_t>( b.el_eta() -> at( b.ll_i1() -> at(0) ) );
      l2eta = static_cast<Float_t>( b.el_eta() -> at( b.ll_i2() -> at(0) ) );
    } else if( b.ll_lepid() -> at(0) == 13 ){
      l1eta = static_cast<Float_t>( b.mu_eta() -> at( b.ll_i1() -> at(0) ) );
      l2eta = static_cast<Float_t>( b.mu_eta() -> at( b.ll_i2() -> at(0) ) );
    } else {
      return -2; //Another redundancy
    }
   
    //cout << "Past conditional statements" << endl;
  
  
    //Set all the values of the constants to change values of pointers
    //MO - switching definitions to match training
    ct = static_cast<Float_t>(b.llphoton_costheta()->at(0));
    cT = static_cast<Float_t>(b.llphoton_cosTheta()->at(0));
    pdmax = static_cast<Float_t>(pdrmax(b));
    gphi = static_cast<Float_t>(b.llphoton_psi()->at(0)); 
    yrV = static_cast<Float_t>(b.llphoton_pt()->at(0)/b.llphoton_m()->at(0) ); 
    yres = static_cast<Float_t>((b.photon_pterr()->at(0))/(b.photon_pt()->at(0)));
    yrap = static_cast<Float_t>(b.photon_eta()->at(0));
    yidmva = static_cast<Float_t>(b.photon_idmva() -> at(0));
    pdmin = static_cast<Float_t>(b.photon_drmin() -> at(0));
  
  
    //cout << "before evaluate" << endl;
  
    //EvaluateMVA to get BDT score
    currentBDTScore = reader -> TMVA::Reader::EvaluateMVA("BDT");
  
    //cout << "got value " << endl;
  
    //Return new BDT score 
    return currentBDTScore; 
  }

  void set_vars(const Baby& b) {
    ct = static_cast<Float_t>(cos_theta(b));
    cT = static_cast<Float_t>(cos_Theta(b));
    pdmax = static_cast<Float_t>(pdrmax(b));
    gphi = static_cast<Float_t>(Getphi(b)); 
    yrV = static_cast<Float_t>(b.llphoton_pt()->at(0)/b.llphoton_m()->at(0) ); 
    yidmva = static_cast<Float_t>(b.photon_idmva() -> at(0));
    pdmin = static_cast<Float_t>(b.photon_drmin() -> at(0));
    yres = static_cast<Float_t>((b.photon_pterr()->at(0))/(b.photon_pt()->at(0)));
    yrap = static_cast<Float_t>(b.photon_eta()->at(0));
    l1eta = get_l1_rapidity(b);
    l2eta = get_l2_rapidity(b);
  }
  
  //Wrapper around mva_score function so it can be returned easily in draw_pico
  NamedFunc MVAScore("MVAScore",[](const Baby &b) -> NamedFunc::ScalarType{ 
      return mva_score(b);
  });

  TMVA::Reader boosted_tagger_reader;
  TMVA::Reader decorr_reader;
  TMVA::Reader decorrratio_reader;
  TMVA::Reader decorrratio_nomasscut_reader;
  TMVA::Reader phptl35_reader;
  TMVA::Reader phptg35_reader;
  TMVA::Reader topfour_reader;
  TMVA::Reader topfive_reader;
  TMVA::Reader topseven_reader;
  TMVA::Reader topnine_reader;
  TMVA::Reader topeleven_reader;
  TMVA::Reader ptbias_reader;
}

int main() {

  //setup
  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");
  Process::Type back =  Process::Type::background;
  //Process::Type data =  Process::Type::data;
  Process::Type sig =  Process::Type::signal;

  Float_t *ct_ptr = &ct;//,cT,pdmax,gphi,yrV,yres,yrap,yidmva,pdmin,l1eta,l2eta;
  Float_t *cT_ptr = &cT;
  Float_t *pdmax_ptr = &pdmax;
  Float_t *gphi_ptr = &gphi;
  Float_t *yrV_ptr = &yrV;
  Float_t *yres_ptr = &yres;
  Float_t *yrap_ptr = &yrap;
  Float_t *yidmva_ptr = &yidmva;
  Float_t *pdmin_ptr = &pdmin;
  Float_t *l1eta_ptr = &l1eta;
  Float_t *l2eta_ptr = &l2eta;

  reader -> TMVA::Reader::AddVariable("photon_mva",yidmva_ptr);//static_cast<Float_t*>(&(b.photon_idmva() -> at(0))) );
  reader -> TMVA::Reader::AddVariable("min_dR", pdmin_ptr); //static_cast<Float_t*>(&(b.photon_drmin() -> at(0))) );
  reader -> TMVA::Reader::AddVariable("max_dR", pdmax_ptr); //&(static_cast<Float_t>(pdrmax(b))) );
  reader -> TMVA::Reader::AddVariable("pt_mass", yrV_ptr); //(static_cast<Float_t>( y_RatVar(b) )) );
  reader -> TMVA::Reader::AddVariable("cosTheta", cT_ptr);//(static_cast<Float_t>(cos_Theta(b))) );
  reader -> TMVA::Reader::AddVariable("costheta", ct_ptr);//(static_cast<Float_t>( cos_theta(b))) );
  reader -> TMVA::Reader::AddVariable("phi",gphi_ptr);//(static_cast<Float_t>( Getphi(b)))) );
  reader -> TMVA::Reader::AddVariable("photon_res",yres_ptr); //static_cast<Float_t*>( &( (b.photon_pterr()->at(0))/(b.photon_pt()->at(0)) ) ) );
  reader -> TMVA::Reader::AddVariable("photon_rapidity",yrap_ptr); 
  reader -> TMVA::Reader::AddVariable("l1_rapidity", l1eta_ptr); //static_cast<Float_t*>(&(b.lep_eta() -> at(0))) );
  reader -> TMVA::Reader::AddVariable("l2_rapidity", l2eta_ptr); // static_cast<Float_t*>(&(b.lep_eta() -> at(1))) );
  //reader -> TMVA::Reader::BookMVA("BDT","/homes/abarzdukas/ZGamma/kinematicMVAs/TMVAClassification_BDTG_DM.weights.xml");
  //reader -> TMVA::Reader::BookMVA("BDT","/homes/oshiro/analysis/mltools/dataset/weights/zgnodecorr_bdt_BDT.weights.xml");
  //reader -> TMVA::Reader::BookMVA("BDT","/homes/oshiro/analysis/mltools/dataset/weights/zg_idmva_nodecorr_bdt_BDT.weights.xml");
  reader -> TMVA::Reader::BookMVA("BDT","/homes/oshiro/analysis/small_phys_utils/dataset/weights/shuffled_kinbdt_masscut_phpt0_BDT.weights.xml");

  float* boost_tagger_hgpt = new float();
  float* boost_tagger_hgdr = new float();
  float* boost_tagger_zpt = new float();
  float* boost_tagger_phpt = new float();
  boosted_tagger_reader.AddVariable("higgspt",boost_tagger_hgpt);
  boosted_tagger_reader.AddVariable("higgsdr",boost_tagger_hgdr);
  boosted_tagger_reader.AddVariable("zpt", boost_tagger_zpt);
  boosted_tagger_reader.AddVariable("phpt",boost_tagger_phpt);
  boosted_tagger_reader.BookMVA("BDT","/homes/oshiro/analysis/mltools/dataset/weights/zgboost_bdt_BDT.weights.xml");

  NamedFunc BoostScore("BoostScore",[boost_tagger_hgpt,boost_tagger_hgdr,boost_tagger_zpt,boost_tagger_phpt](const Baby &b) -> NamedFunc::ScalarType{ 
    *boost_tagger_hgpt = b.llphoton_pt()->at(0);
    *boost_tagger_hgdr = b.llphoton_dr()->at(0);
    *boost_tagger_zpt = b.ll_pt()->at(0);
    *boost_tagger_phpt = b.photon_pt()->at(0);
    return boosted_tagger_reader.EvaluateMVA("BDT");
  });

  float* decorr_photon_pt_ptr = new float();
  decorr_reader.AddVariable("photon_mva",yidmva_ptr);//static_cast<Float_t*>(&(b.photon_idmva() -> at(0))) );
  decorr_reader.AddVariable("min_dR", pdmin_ptr); //static_cast<Float_t*>(&(b.photon_drmin() -> at(0))) );
  decorr_reader.AddVariable("max_dR", pdmax_ptr); //&(static_cast<Float_t>(pdrmax(b))) );
  decorr_reader.AddVariable("pt_mass", yrV_ptr); //(static_cast<Float_t>( y_RatVar(b) )) );
  decorr_reader.AddVariable("cosTheta", cT_ptr);//(static_cast<Float_t>(cos_Theta(b))) );
  decorr_reader.AddVariable("costheta", ct_ptr);//(static_cast<Float_t>( cos_theta(b))) );
  decorr_reader.AddVariable("phi",gphi_ptr);//(static_cast<Float_t>( Getphi(b)))) );
  decorr_reader.AddVariable("photon_res",yres_ptr); //static_cast<Float_t*>( &( (b.photon_pterr()->at(0))/(b.photon_pt()->at(0)) ) ) );
  decorr_reader.AddVariable("photon_rapidity",yrap_ptr); 
  decorr_reader.AddVariable("l1_rapidity", l1eta_ptr); //static_cast<Float_t*>(&(b.lep_eta() -> at(0))) );
  decorr_reader.AddVariable("l2_rapidity", l2eta_ptr); // static_cast<Float_t*>(&(b.lep_eta() -> at(1))) );
  decorr_reader.AddVariable("decorr_photon_pt", decorr_photon_pt_ptr); // static_cast<Float_t*>(&(b.lep_eta() -> at(1))) );
  decorr_reader.BookMVA("BDT","/homes/oshiro/analysis/mltools/dataset/weights/zg_idmva_lindecorr_bdt_BDT.weights.xml");

  float* photon_pt_mass_ptr = new float();
  decorrratio_reader.AddVariable("photon_mva",yidmva_ptr);//static_cast<Float_t*>(&(b.photon_idmva() -> at(0))) );
  decorrratio_reader.AddVariable("min_dR", pdmin_ptr); //static_cast<Float_t*>(&(b.photon_drmin() -> at(0))) );
  decorrratio_reader.AddVariable("max_dR", pdmax_ptr); //&(static_cast<Float_t>(pdrmax(b))) );
  decorrratio_reader.AddVariable("pt_mass", yrV_ptr); //(static_cast<Float_t>( y_RatVar(b) )) );
  decorrratio_reader.AddVariable("cosTheta", cT_ptr);//(static_cast<Float_t>(cos_Theta(b))) );
  decorrratio_reader.AddVariable("costheta", ct_ptr);//(static_cast<Float_t>( cos_theta(b))) );
  decorrratio_reader.AddVariable("phi",gphi_ptr);//(static_cast<Float_t>( Getphi(b)))) );
  decorrratio_reader.AddVariable("photon_res",yres_ptr); //static_cast<Float_t*>( &( (b.photon_pterr()->at(0))/(b.photon_pt()->at(0)) ) ) );
  decorrratio_reader.AddVariable("photon_rapidity",yrap_ptr); 
  decorrratio_reader.AddVariable("l1_rapidity", l1eta_ptr); //static_cast<Float_t*>(&(b.lep_eta() -> at(0))) );
  decorrratio_reader.AddVariable("l2_rapidity", l2eta_ptr); // static_cast<Float_t*>(&(b.lep_eta() -> at(1))) );
  decorrratio_reader.AddVariable("photon_pt_mass", photon_pt_mass_ptr); // static_cast<Float_t*>(&(b.lep_eta() -> at(1))) );
  //decorrratio_reader.BookMVA("BDT","/homes/oshiro/analysis/mltools/dataset/weights/zg_idmva_ratdecorr_bdt_BDT.weights.xml");
  decorrratio_reader.BookMVA("BDT","/homes/oshiro/analysis/small_phys_utils/dataset/weights/shuffled_kinbdt_masscut_phpt1_BDT.weights.xml");
  
  decorrratio_nomasscut_reader.AddVariable("photon_mva",yidmva_ptr);//static_cast<Float_t*>(&(b.photon_idmva() -> at(0))) );
  decorrratio_nomasscut_reader.AddVariable("min_dR", pdmin_ptr); //static_cast<Float_t*>(&(b.photon_drmin() -> at(0))) );
  decorrratio_nomasscut_reader.AddVariable("max_dR", pdmax_ptr); //&(static_cast<Float_t>(pdrmax(b))) );
  decorrratio_nomasscut_reader.AddVariable("pt_mass", yrV_ptr); //(static_cast<Float_t>( y_RatVar(b) )) );
  decorrratio_nomasscut_reader.AddVariable("cosTheta", cT_ptr);//(static_cast<Float_t>(cos_Theta(b))) );
  decorrratio_nomasscut_reader.AddVariable("costheta", ct_ptr);//(static_cast<Float_t>( cos_theta(b))) );
  decorrratio_nomasscut_reader.AddVariable("phi",gphi_ptr);//(static_cast<Float_t>( Getphi(b)))) );
  decorrratio_nomasscut_reader.AddVariable("photon_res",yres_ptr); //static_cast<Float_t*>( &( (b.photon_pterr()->at(0))/(b.photon_pt()->at(0)) ) ) );
  decorrratio_nomasscut_reader.AddVariable("photon_rapidity",yrap_ptr); 
  decorrratio_nomasscut_reader.AddVariable("l1_rapidity", l1eta_ptr); //static_cast<Float_t*>(&(b.lep_eta() -> at(0))) );
  decorrratio_nomasscut_reader.AddVariable("l2_rapidity", l2eta_ptr); // static_cast<Float_t*>(&(b.lep_eta() -> at(1))) );
  decorrratio_nomasscut_reader.AddVariable("photon_pt_mass", photon_pt_mass_ptr); // static_cast<Float_t*>(&(b.lep_eta() -> at(1))) );
  decorrratio_nomasscut_reader.BookMVA("BDT","/homes/oshiro/analysis/mltools/dataset/weights/zg_idmva_ratdecorr_nomasscut_bdt_BDT.weights.xml");

  phptl35_reader.AddVariable("photon_mva",yidmva_ptr);//static_cast<Float_t*>(&(b.photon_idmva() -> at(0))) );
  phptl35_reader.AddVariable("min_dR", pdmin_ptr); //static_cast<Float_t*>(&(b.photon_drmin() -> at(0))) );
  phptl35_reader.AddVariable("max_dR", pdmax_ptr); //&(static_cast<Float_t>(pdrmax(b))) );
  phptl35_reader.AddVariable("pt_mass", yrV_ptr); //(static_cast<Float_t>( y_RatVar(b) )) );
  phptl35_reader.AddVariable("cosTheta", cT_ptr);//(static_cast<Float_t>(cos_Theta(b))) );
  phptl35_reader.AddVariable("costheta", ct_ptr);//(static_cast<Float_t>( cos_theta(b))) );
  phptl35_reader.AddVariable("phi",gphi_ptr);//(static_cast<Float_t>( Getphi(b)))) );
  phptl35_reader.AddVariable("photon_res",yres_ptr); //static_cast<Float_t*>( &( (b.photon_pterr()->at(0))/(b.photon_pt()->at(0)) ) ) );
  phptl35_reader.AddVariable("photon_rapidity",yrap_ptr); 
  phptl35_reader.AddVariable("l1_rapidity", l1eta_ptr); //static_cast<Float_t*>(&(b.lep_eta() -> at(0))) );
  phptl35_reader.AddVariable("l2_rapidity", l2eta_ptr); // static_cast<Float_t*>(&(b.lep_eta() -> at(1))) );
  //phptl35_reader.BookMVA("BDT","/homes/oshiro/analysis/mltools/dataset/weights/zg_idmva_ratdecorr_phptl35_bdt_BDT.weights.xml");
  phptl35_reader.BookMVA("BDT","/homes/oshiro/analysis/small_phys_utils/dataset/weights/shuffled_kinbdt_masscut_phptl35_phpt2_BDT.weights.xml");

  phptg35_reader.AddVariable("photon_mva",yidmva_ptr);//static_cast<Float_t*>(&(b.photon_idmva() -> at(0))) );
  phptg35_reader.AddVariable("min_dR", pdmin_ptr); //static_cast<Float_t*>(&(b.photon_drmin() -> at(0))) );
  phptg35_reader.AddVariable("max_dR", pdmax_ptr); //&(static_cast<Float_t>(pdrmax(b))) );
  phptg35_reader.AddVariable("pt_mass", yrV_ptr); //(static_cast<Float_t>( y_RatVar(b) )) );
  phptg35_reader.AddVariable("cosTheta", cT_ptr);//(static_cast<Float_t>(cos_Theta(b))) );
  phptg35_reader.AddVariable("costheta", ct_ptr);//(static_cast<Float_t>( cos_theta(b))) );
  phptg35_reader.AddVariable("phi",gphi_ptr);//(static_cast<Float_t>( Getphi(b)))) );
  phptg35_reader.AddVariable("photon_res",yres_ptr); //static_cast<Float_t*>( &( (b.photon_pterr()->at(0))/(b.photon_pt()->at(0)) ) ) );
  phptg35_reader.AddVariable("photon_rapidity",yrap_ptr); 
  phptg35_reader.AddVariable("l1_rapidity", l1eta_ptr); //static_cast<Float_t*>(&(b.lep_eta() -> at(0))) );
  phptg35_reader.AddVariable("l2_rapidity", l2eta_ptr); // static_cast<Float_t*>(&(b.lep_eta() -> at(1))) );
  //phptg35_reader.BookMVA("BDT","/homes/oshiro/analysis/mltools/dataset/weights/zg_idmva_ratdecorr_phptg35_bdt_BDT.weights.xml");
  phptg35_reader.BookMVA("BDT","/homes/oshiro/analysis/small_phys_utils/dataset/weights/shuffled_kinbdt_masscut_phptg35_phpt3_BDT.weights.xml");

  topfour_reader.AddVariable("photon_mva",yidmva_ptr);//static_cast<Float_t*>(&(b.photon_idmva() -> at(0))) );
  topfour_reader.AddVariable("min_dR", pdmin_ptr); //static_cast<Float_t*>(&(b.photon_drmin() -> at(0))) );
  topfour_reader.AddVariable("pt_mass", yrV_ptr); //(static_cast<Float_t>( y_RatVar(b) )) );
  topfour_reader.AddVariable("cosTheta", cT_ptr);//(static_cast<Float_t>(cos_Theta(b))) );
  topfour_reader.BookMVA("BDT","/homes/oshiro/analysis/small_phys_utils/dataset/weights/kinbdt_masscut_idmvacut_topfour_200trees_BDT.weights.xml");

  topfive_reader.AddVariable("photon_mva",yidmva_ptr);//static_cast<Float_t*>(&(b.photon_idmva() -> at(0))) );
  topfive_reader.AddVariable("min_dR", pdmin_ptr); //static_cast<Float_t*>(&(b.photon_drmin() -> at(0))) );
  topfive_reader.AddVariable("pt_mass", yrV_ptr); //(static_cast<Float_t>( y_RatVar(b) )) );
  topfive_reader.AddVariable("cosTheta", cT_ptr);//(static_cast<Float_t>(cos_Theta(b))) );
  topfive_reader.AddVariable("photon_rapidity",yrap_ptr); 
  //topfive_reader.BookMVA("BDT","/homes/oshiro/analysis/mltools/dataset/weights/zg_idmva_topfour_rapid_bdt_BDT.weights.xml");
  topfive_reader.BookMVA("BDT","/homes/oshiro/analysis/small_phys_utils/dataset/weights/kinbdt_masscut_idmvacut_topfive_100trees_BDT.weights.xml");

  topseven_reader.AddVariable("photon_mva",yidmva_ptr);//static_cast<Float_t*>(&(b.photon_idmva() -> at(0))) );
  topseven_reader.AddVariable("min_dR", pdmin_ptr); //static_cast<Float_t*>(&(b.photon_drmin() -> at(0))) );
  topseven_reader.AddVariable("pt_mass", yrV_ptr); //(static_cast<Float_t>( y_RatVar(b) )) );
  topseven_reader.AddVariable("cosTheta", cT_ptr);//(static_cast<Float_t>(cos_Theta(b))) );
  topseven_reader.AddVariable("photon_rapidity",yrap_ptr); 
  topseven_reader.AddVariable("l1_rapidity", l1eta_ptr); //static_cast<Float_t*>(&(b.lep_eta() -> at(0))) );
  topseven_reader.AddVariable("l2_rapidity", l2eta_ptr); // static_cast<Float_t*>(&(b.lep_eta() -> at(1))) );
  topseven_reader.BookMVA("BDT","/homes/oshiro/analysis/small_phys_utils/dataset/weights/kinbdt_masscut_idmvacut_topseven_100trees_BDT.weights.xml");

  topnine_reader.AddVariable("photon_mva",yidmva_ptr);//static_cast<Float_t*>(&(b.photon_idmva() -> at(0))) );
  topnine_reader.AddVariable("min_dR", pdmin_ptr); //static_cast<Float_t*>(&(b.photon_drmin() -> at(0))) );
  topnine_reader.AddVariable("pt_mass", yrV_ptr); //(static_cast<Float_t>( y_RatVar(b) )) );
  topnine_reader.AddVariable("cosTheta", cT_ptr);//(static_cast<Float_t>(cos_Theta(b))) );
  topnine_reader.AddVariable("costheta", ct_ptr);//(static_cast<Float_t>(cos_Theta(b))) );
  topnine_reader.AddVariable("photon_res", yres_ptr);//(static_cast<Float_t>(cos_Theta(b))) );
  topnine_reader.AddVariable("photon_rapidity",yrap_ptr); 
  topnine_reader.AddVariable("l1_rapidity", l1eta_ptr); //static_cast<Float_t*>(&(b.lep_eta() -> at(0))) );
  topnine_reader.AddVariable("l2_rapidity", l2eta_ptr); // static_cast<Float_t*>(&(b.lep_eta() -> at(1))) );
  topnine_reader.BookMVA("BDT","/homes/oshiro/analysis/small_phys_utils/dataset/weights/kinbdt_masscut_idmvacut_topnine_100trees_BDT.weights.xml");

  topeleven_reader.AddVariable("photon_mva",yidmva_ptr);//static_cast<Float_t*>(&(b.photon_idmva() -> at(0))) );
  topeleven_reader.AddVariable("min_dR", pdmin_ptr); //static_cast<Float_t*>(&(b.photon_drmin() -> at(0))) );
  topeleven_reader.AddVariable("max_dR", pdmax_ptr); //static_cast<Float_t*>(&(b.photon_drmin() -> at(0))) );
  topeleven_reader.AddVariable("pt_mass", yrV_ptr); //(static_cast<Float_t>( y_RatVar(b) )) );
  topeleven_reader.AddVariable("cosTheta", cT_ptr);//(static_cast<Float_t>(cos_Theta(b))) );
  topeleven_reader.AddVariable("costheta", ct_ptr);//(static_cast<Float_t>(cos_Theta(b))) );
  topeleven_reader.AddVariable("phi",gphi_ptr);//(static_cast<Float_t>( Getphi(b)))) );
  topeleven_reader.AddVariable("photon_res", yres_ptr);//(static_cast<Float_t>(cos_Theta(b))) );
  topeleven_reader.AddVariable("photon_rapidity",yrap_ptr); 
  topeleven_reader.AddVariable("l1_rapidity", l1eta_ptr); //static_cast<Float_t*>(&(b.lep_eta() -> at(0))) );
  topeleven_reader.AddVariable("l2_rapidity", l2eta_ptr); // static_cast<Float_t*>(&(b.lep_eta() -> at(1))) );
  topeleven_reader.BookMVA("BDT","/homes/oshiro/analysis/small_phys_utils/dataset/weights/kinbdt_masscut_idmvacut_topeleven_200trees_BDT.weights.xml");

  float* photon_pt_ptr = new float();
  ptbias_reader.AddVariable("photon_mva",yidmva_ptr);//static_cast<Float_t*>(&(b.photon_idmva() -> at(0))) );
  ptbias_reader.AddVariable("min_dR", pdmin_ptr); //static_cast<Float_t*>(&(b.photon_drmin() -> at(0))) );
  ptbias_reader.AddVariable("max_dR", pdmax_ptr); //&(static_cast<Float_t>(pdrmax(b))) );
  ptbias_reader.AddVariable("pt_mass", yrV_ptr); //(static_cast<Float_t>( y_RatVar(b) )) );
  ptbias_reader.AddVariable("cosTheta", cT_ptr);//(static_cast<Float_t>(cos_Theta(b))) );
  ptbias_reader.AddVariable("costheta", ct_ptr);//(static_cast<Float_t>( cos_theta(b))) );
  ptbias_reader.AddVariable("phi",gphi_ptr);//(static_cast<Float_t>( Getphi(b)))) );
  ptbias_reader.AddVariable("photon_res",yres_ptr); //static_cast<Float_t*>( &( (b.photon_pterr()->at(0))/(b.photon_pt()->at(0)) ) ) );
  ptbias_reader.AddVariable("photon_rapidity",yrap_ptr); 
  ptbias_reader.AddVariable("l1_rapidity", l1eta_ptr); //static_cast<Float_t*>(&(b.lep_eta() -> at(0))) );
  ptbias_reader.AddVariable("l2_rapidity", l2eta_ptr); // static_cast<Float_t*>(&(b.lep_eta() -> at(1))) );
  ptbias_reader.AddVariable("photon_ptransverse", photon_pt_ptr); // static_cast<Float_t*>(&(b.lep_eta() -> at(1))) );
  ptbias_reader.BookMVA("BDT","/homes/oshiro/analysis/mltools/dataset/weights/zg_idmva_ptbias_bdt_BDT.weights.xml");

  NamedFunc DecorrBdtScore("DecorrBdtScore",[decorr_photon_pt_ptr](const Baby &b) -> NamedFunc::ScalarType{ 
    if( b.ll_lepid() -> at(0) == 11){
      l1eta = static_cast<Float_t>( b.el_eta() -> at( b.ll_i1() -> at(0) ) );
      l2eta = static_cast<Float_t>( b.el_eta() -> at( b.ll_i2() -> at(0) ) );
    } else if( b.ll_lepid() -> at(0) == 13 ){
      l1eta = static_cast<Float_t>( b.mu_eta() -> at( b.ll_i1() -> at(0) ) );
      l2eta = static_cast<Float_t>( b.mu_eta() -> at( b.ll_i2() -> at(0) ) );
    }
    ct = static_cast<Float_t>(cos_theta(b));
    cT = static_cast<Float_t>(cos_Theta(b));
    pdmax = static_cast<Float_t>(pdrmax(b));
    gphi = static_cast<Float_t>(Getphi(b)); 
    yrV = static_cast<Float_t>(b.llphoton_pt()->at(0)/b.llphoton_m()->at(0) ); 
    yres = static_cast<Float_t>((b.photon_pterr()->at(0))/(b.photon_pt()->at(0)));
    yrap = static_cast<Float_t>(b.photon_eta()->at(0));
    yidmva = static_cast<Float_t>(b.photon_idmva() -> at(0));
    pdmin = static_cast<Float_t>(b.photon_drmin() -> at(0));
    *decorr_photon_pt_ptr = b.photon_pt()->at(0)-0.207*b.llphoton_m()->at(0);
    return decorr_reader.EvaluateMVA("BDT");
  });

  NamedFunc DecorrRatioBdtScore("DecorrRatioBdtScore",[photon_pt_mass_ptr](const Baby &b) -> NamedFunc::ScalarType{ 
    if( b.ll_lepid() -> at(0) == 11){
      l1eta = static_cast<Float_t>( b.el_eta() -> at( b.ll_i1() -> at(0) ) );
      l2eta = static_cast<Float_t>( b.el_eta() -> at( b.ll_i2() -> at(0) ) );
    } else if( b.ll_lepid() -> at(0) == 13 ){
      l1eta = static_cast<Float_t>( b.mu_eta() -> at( b.ll_i1() -> at(0) ) );
      l2eta = static_cast<Float_t>( b.mu_eta() -> at( b.ll_i2() -> at(0) ) );
    }
    ct = static_cast<Float_t>(cos_theta(b));
    cT = static_cast<Float_t>(cos_Theta(b));
    pdmax = static_cast<Float_t>(pdrmax(b));
    gphi = static_cast<Float_t>(Getphi(b)); 
    yrV = static_cast<Float_t>(b.llphoton_pt()->at(0)/b.llphoton_m()->at(0) ); 
    yres = static_cast<Float_t>((b.photon_pterr()->at(0))/(b.photon_pt()->at(0)));
    yrap = static_cast<Float_t>(b.photon_eta()->at(0));
    yidmva = static_cast<Float_t>(b.photon_idmva() -> at(0));
    pdmin = static_cast<Float_t>(b.photon_drmin() -> at(0));
    *photon_pt_mass_ptr = b.photon_pt()->at(0)/b.llphoton_m()->at(0);
    return decorrratio_reader.EvaluateMVA("BDT");
  });

  NamedFunc NodecorrBdtScore("NodecorrBdtScore",[photon_pt_ptr](const Baby &b) -> NamedFunc::ScalarType{ 
    if( b.ll_lepid() -> at(0) == 11){
      l1eta = static_cast<Float_t>( b.el_eta() -> at( b.ll_i1() -> at(0) ) );
      l2eta = static_cast<Float_t>( b.el_eta() -> at( b.ll_i2() -> at(0) ) );
    } else if( b.ll_lepid() -> at(0) == 13 ){
      l1eta = static_cast<Float_t>( b.mu_eta() -> at( b.ll_i1() -> at(0) ) );
      l2eta = static_cast<Float_t>( b.mu_eta() -> at( b.ll_i2() -> at(0) ) );
    }
    ct = static_cast<Float_t>(cos_theta(b));
    cT = static_cast<Float_t>(cos_Theta(b));
    pdmax = static_cast<Float_t>(pdrmax(b));
    gphi = static_cast<Float_t>(Getphi(b)); 
    yrV = static_cast<Float_t>(b.llphoton_pt()->at(0)/b.llphoton_m()->at(0) ); 
    yres = static_cast<Float_t>((b.photon_pterr()->at(0))/(b.photon_pt()->at(0)));
    yrap = static_cast<Float_t>(b.photon_eta()->at(0));
    yidmva = static_cast<Float_t>(b.photon_idmva() -> at(0));
    pdmin = static_cast<Float_t>(b.photon_drmin() -> at(0));
    *photon_pt_ptr = b.photon_pt()->at(0);
    return ptbias_reader.EvaluateMVA("BDT");
  });

  NamedFunc DecorrRatioNoMassCutBdtScore("DecorrRatioNoMassCutBdtScore",[photon_pt_mass_ptr](const Baby &b) -> NamedFunc::ScalarType{ 
    if( b.ll_lepid() -> at(0) == 11){
      l1eta = static_cast<Float_t>( b.el_eta() -> at( b.ll_i1() -> at(0) ) );
      l2eta = static_cast<Float_t>( b.el_eta() -> at( b.ll_i2() -> at(0) ) );
    } else if( b.ll_lepid() -> at(0) == 13 ){
      l1eta = static_cast<Float_t>( b.mu_eta() -> at( b.ll_i1() -> at(0) ) );
      l2eta = static_cast<Float_t>( b.mu_eta() -> at( b.ll_i2() -> at(0) ) );
    }
    ct = static_cast<Float_t>(cos_theta(b));
    cT = static_cast<Float_t>(cos_Theta(b));
    pdmax = static_cast<Float_t>(pdrmax(b));
    gphi = static_cast<Float_t>(Getphi(b)); 
    yrV = static_cast<Float_t>(b.llphoton_pt()->at(0)/b.llphoton_m()->at(0) ); 
    yres = static_cast<Float_t>((b.photon_pterr()->at(0))/(b.photon_pt()->at(0)));
    yrap = static_cast<Float_t>(b.photon_eta()->at(0));
    yidmva = static_cast<Float_t>(b.photon_idmva() -> at(0));
    pdmin = static_cast<Float_t>(b.photon_drmin() -> at(0));
    *photon_pt_mass_ptr = b.photon_pt()->at(0)/b.llphoton_m()->at(0);
    return decorrratio_nomasscut_reader.EvaluateMVA("BDT");
  });

  NamedFunc LowPhPtBdtScore("LowPhPtBdtScore",[photon_pt_mass_ptr](const Baby &b) -> NamedFunc::ScalarType{ 
    if( b.ll_lepid() -> at(0) == 11){
      l1eta = static_cast<Float_t>( b.el_eta() -> at( b.ll_i1() -> at(0) ) );
      l2eta = static_cast<Float_t>( b.el_eta() -> at( b.ll_i2() -> at(0) ) );
    } else if( b.ll_lepid() -> at(0) == 13 ){
      l1eta = static_cast<Float_t>( b.mu_eta() -> at( b.ll_i1() -> at(0) ) );
      l2eta = static_cast<Float_t>( b.mu_eta() -> at( b.ll_i2() -> at(0) ) );
    }
    ct = static_cast<Float_t>(cos_theta(b));
    cT = static_cast<Float_t>(cos_Theta(b));
    pdmax = static_cast<Float_t>(pdrmax(b));
    gphi = static_cast<Float_t>(Getphi(b)); 
    yrV = static_cast<Float_t>(b.llphoton_pt()->at(0)/b.llphoton_m()->at(0) ); 
    yres = static_cast<Float_t>((b.photon_pterr()->at(0))/(b.photon_pt()->at(0)));
    yrap = static_cast<Float_t>(b.photon_eta()->at(0));
    yidmva = static_cast<Float_t>(b.photon_idmva() -> at(0));
    pdmin = static_cast<Float_t>(b.photon_drmin() -> at(0));
    return phptl35_reader.EvaluateMVA("BDT");
  });

  NamedFunc HighPhPtBdtScore("HighPhPtBdtScore",[photon_pt_mass_ptr](const Baby &b) -> NamedFunc::ScalarType{ 
    if( b.ll_lepid() -> at(0) == 11){
      l1eta = static_cast<Float_t>( b.el_eta() -> at( b.ll_i1() -> at(0) ) );
      l2eta = static_cast<Float_t>( b.el_eta() -> at( b.ll_i2() -> at(0) ) );
    } else if( b.ll_lepid() -> at(0) == 13 ){
      l1eta = static_cast<Float_t>( b.mu_eta() -> at( b.ll_i1() -> at(0) ) );
      l2eta = static_cast<Float_t>( b.mu_eta() -> at( b.ll_i2() -> at(0) ) );
    }
    ct = static_cast<Float_t>(cos_theta(b));
    cT = static_cast<Float_t>(cos_Theta(b));
    pdmax = static_cast<Float_t>(pdrmax(b));
    gphi = static_cast<Float_t>(Getphi(b)); 
    yrV = static_cast<Float_t>(b.llphoton_pt()->at(0)/b.llphoton_m()->at(0) ); 
    yres = static_cast<Float_t>((b.photon_pterr()->at(0))/(b.photon_pt()->at(0)));
    yrap = static_cast<Float_t>(b.photon_eta()->at(0));
    yidmva = static_cast<Float_t>(b.photon_idmva() -> at(0));
    pdmin = static_cast<Float_t>(b.photon_drmin() -> at(0));
    return phptg35_reader.EvaluateMVA("BDT");
  });

  NamedFunc TopFourBdtScore("TopFourBdtScore",[](const Baby &b) -> NamedFunc::ScalarType{ 
    cT = static_cast<Float_t>(cos_Theta(b));
    yrV = static_cast<Float_t>(b.llphoton_pt()->at(0)/b.llphoton_m()->at(0) ); 
    yidmva = static_cast<Float_t>(b.photon_idmva() -> at(0));
    pdmin = static_cast<Float_t>(b.photon_drmin() -> at(0));
    return topfour_reader.EvaluateMVA("BDT");
  });

  NamedFunc TopFiveBdtScore("TopFiveBdtScore",[](const Baby &b) -> NamedFunc::ScalarType{ 
    set_vars(b);
    return topfive_reader.EvaluateMVA("BDT");
  });

  NamedFunc TopSevenBdtScore("TopSevenBdtScore",[](const Baby &b) -> NamedFunc::ScalarType{ 
    set_vars(b);
    return topseven_reader.EvaluateMVA("BDT");
  });

  NamedFunc TopNineBdtScore("TopNineBdtScore",[](const Baby &b) -> NamedFunc::ScalarType{ 
    set_vars(b);
    return topnine_reader.EvaluateMVA("BDT");
  });

  NamedFunc TopElevenBdtScore("TopElevenBdtScore",[](const Baby &b) -> NamedFunc::ScalarType{ 
    set_vars(b);
    return topeleven_reader.EvaluateMVA("BDT");
  });

  string bfolder("/net/cms17/cms17r0/pico/");
  string mc_path( bfolder+"NanoAODv9/htozgamma_deathvalley_v1/2017/mc/skim_llg/");
  string sig_path(bfolder+"NanoAODv2/zgamma_signal/2017/signal/skim_llg/");
  string mc_slim_path( bfolder+"NanoAODv9/htozgamma_deathvalley_v1/2017/mc/merged_zgmc_llg/");
  string sig_slim_path(bfolder+"NanoAODv2/zgamma_signal/2017/signal/merged_zgmc_llg/");

  //NamedFunc el_trigs = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL||HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ||HLT_Ele35_WPTight_Gsf||HLT_Ele32_WPTight_Gsf_L1DoubleEG";
  //NamedFunc mu_trigs = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8||HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8||HLT_IsoMu27||HLT_IsoMu24";
  NamedFunc el_trigs = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL||HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ";
  NamedFunc mu_trigs = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ";
  //NamedFunc el_trigs("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
  //NamedFunc mu_trigs("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8");
  NamedFunc trigs(el_trigs || mu_trigs);

  //NamedFunc to remove events with no ll or llphoton objects, which can happen in prodcutions from prasanna's earlier branch with ll sign req.
  NamedFunc sanitize("sanitize",[](const Baby &b) -> NamedFunc::ScalarType{ 
    unsigned ll_size = b.ll_m()->size();
    unsigned llg_size = b.llphoton_m()->size();
    if (ll_size == 0 || llg_size == 0) return 0;
    return 1;
  });

  NamedFunc z_to_ll("z_to_ll",[](const Baby &b) -> NamedFunc::ScalarType{ 
    for (unsigned imc = 0; imc < b.mc_pt()->size(); imc++) {
      if ((abs(b.mc_id()->at(imc))==11 || abs(b.mc_id()->at(imc))==13)&&b.mc_mom()->at(imc)==23) {
        int mom_idx = b.mc_momidx()->at(imc);
        if (mom_idx != -1) {
          if (b.mc_mom()->at(mom_idx)==25) {
            return 1;
          }
        }
      }
    }
    return 0;
  });

  auto proc_smzg        = Process::MakeShared<Baby_pico>("Z+#gamma",       back, 
                             TColor::GetColor("#16bac5"),{mc_path+"*ZGToLLG*"}, trigs);
  auto proc_dy          = Process::MakeShared<Baby_pico>("Z+Fake Photon",               back, 
                             TColor::GetColor("#ffb400"),{mc_path+"*DYJets*"},  trigs && "stitch_dy");
  auto proc_tt          = Process::MakeShared<Baby_pico>("tt",               back, 
                             TColor::GetColor("#ff6600"),{mc_path+"*TTTo2L2Nu*"},  trigs);
  //auto proc_mc_back     = Process::MakeShared<Baby_pico>("Pseudodata", data, 
  //                         kBlack, {mc_path+"*ZGToLLG*", mc_path+"*LLAJJ*",mc_path+"*DYJets*"},
  //                         trigs && "(type != 6200 || stitch_dy)");
  auto proc_hzg_gg      = Process::MakeShared<Baby_pico>("gg#rightarrow H#rightarrow Z#gamma", sig, 
                           kRed     ,{sig_path+"*GluGluHToZG*"},   trigs&&sanitize);
  auto proc_hzg         = Process::MakeShared<Baby_pico>("H#rightarrow Z#gamma", sig, 
                           kRed ,{sig_path+"*.root"},   trigs&&sanitize);
  auto proc_hzg_bak     = Process::MakeShared<Baby_pico>("H#rightarrow Z#gamma", back, 
                           kRed ,{sig_path+"*.root"},   trigs&&sanitize);
  auto proc_hzg_gg_x20  = Process::MakeShared<Baby_pico>("gg#rightarrow H#rightarrow Z#gamma (x20)", sig, 
                          kRed     ,{sig_path+"*GluGluHToZG*"},   trigs&&sanitize);
  auto proc_hzg_x20     = Process::MakeShared<Baby_pico>("H#rightarrow Z#gamma (x20)", sig, 
                          kRed ,{sig_path+"*.root"},   trigs&&sanitize);
  auto proc_hzg_unskimmed = Process::MakeShared<Baby_pico>("H#rightarrow Z#gamma", back, 
                           kRed ,{sig_path+"/../unskimmed/*.root"}, z_to_ll);
  auto proc_smzg_run2   = Process::MakeShared<Baby_pico>("Z+#gamma",       back, 
                             TColor::GetColor("#16bac5"),{mc_path+"*ZGToLLG*",
                             mc_path+"../../../2016/mc/skim_llg/*ZGToLLG*",mc_path+"../../../2018/mc/skim_llg/*ZGToLLG*"}, trigs);
  auto proc_dy_run2     = Process::MakeShared<Baby_pico>("Z+Fake Photon",               back, 
                             TColor::GetColor("#ffb400"),{mc_path+"*DYJets*",
                             mc_path+"../../../2016/mc/skim_llg/*DYJets*",mc_path+"../../../2018/mc/skim_llg/*DYJets*"}, trigs&&"stitch_dy");
  auto proc_dy_noskim   = Process::MakeShared<Baby_pico>("Z+#gamma",               back, 
                             TColor::GetColor("#16bac5"),{mc_path+"../unskimmed/*DYJets*",mc_path+"../unskimmed/*ZGToLLG*"},  "1");

  auto proc_hzg_valid   = Process::MakeShared<Baby_pico>("H#rightarrow Z#gamma", sig, 
                           kRed ,{sig_slim_path+"*.root"},   trigs&&sanitize);
  auto proc_dyj_valid   = Process::MakeShared<Baby_pico>("Z+Fake Photon", back, 
                           kRed ,{mc_slim_path+"*DYJets*"},   trigs&&sanitize);
  auto proc_zgm_valid   = Process::MakeShared<Baby_pico>("Z+#gamma", back, 
                           kRed ,{mc_slim_path+"*ZGToLLG*"},   trigs&&sanitize);

  auto proc_hzg_debug   = Process::MakeShared<Baby_pico>("H#rightarrow Z#gamma", sig, 
                           kRed ,{sig_path+"*.root"},   "1");
  auto proc_dyj_debug   = Process::MakeShared<Baby_pico>("Z+Fake Photon", back, 
                           TColor::GetColor("#ffb400") ,{mc_path+"*DYJets*"},   "1");
  auto proc_zgm_debug   = Process::MakeShared<Baby_pico>("Z+#gamma", back, 
                           TColor::GetColor("#16bac5") ,{mc_path+"*ZGToLLG*"},   "1");

  proc_smzg->SetLineWidth(1); proc_dy->SetLineWidth(1); 
  proc_hzg->SetLineWidth(3); proc_hzg_gg->SetLineWidth(3); 
  //proc_mc_back->SetMarkerSize(1); 
  vector<shared_ptr<Process>> procs = {proc_dy, proc_smzg, proc_hzg};
  vector<shared_ptr<Process>> procs_sigx20 = {proc_dy, proc_smzg, proc_hzg_x20};
  vector<shared_ptr<Process>> procs_sigonly = {proc_hzg_bak};
  vector<shared_ptr<Process>> procs_bakonly = {proc_smzg};
  vector<shared_ptr<Process>> procs_allbaks = {proc_tt, proc_smzg, proc_dy};
  vector<shared_ptr<Process>> procs_allprocs = {proc_tt, proc_smzg, proc_dy, proc_hzg};
  vector<shared_ptr<Process>> procs_sig_unskimmed = {proc_hzg_unskimmed};
  vector<shared_ptr<Process>> procs_extrastats = {proc_dy_run2, proc_smzg_run2, proc_hzg};
  vector<shared_ptr<Process>> procs_dy_noskim = {proc_dy_noskim};
  vector<shared_ptr<Process>> procs_debug = {proc_dyj_debug,proc_zgm_debug,proc_hzg_debug};
  PlotOpt lin_lumi("txt/plot_styles.txt","CMSPaper");
  lin_lumi.Title(TitleType::info)
          .Stack(StackType::signal_overlay)
          .Overflow(OverflowType::none)
          .Bottom(BottomType::off)
          .UseCMYK(false)
          .LegendColumns(1)
          .ShowBackgroundError(false)
          .FileExtensions({"pdf"});
  PlotOpt lin_sorb = lin_lumi;
  lin_sorb.Overflow(OverflowType::both).Bottom(BottomType::sorb).FileExtensions({"pdf","root"});
  PlotOpt lin_sorb_upper = lin_lumi;
  lin_sorb_upper.Overflow(OverflowType::both).Bottom(BottomType::sorb_cut_upper).FileExtensions({"pdf","root"});
  PlotOpt lin_shapes = lin_lumi;
  lin_shapes.Stack(StackType::shapes);
  PlotOpt log_lumi = lin_lumi;
  log_lumi.YAxis(YAxisType::log);
  vector<PlotOpt> ops_shapes = {lin_shapes};
  vector<PlotOpt> ops = {lin_lumi};
  vector<PlotOpt> ops_log = {log_lumi};
  vector<PlotOpt> ops_sorb = {lin_sorb};
  vector<PlotOpt> ops_sorb_upper = {lin_sorb_upper};
  PlotOpt twodim_lin_lumi("txt/plot_styles.txt", "Scatter");
  twodim_lin_lumi.Title(TitleType::info)
                 .YAxis(YAxisType::log)
                 .Overflow(OverflowType::overflow)
                 .LogMinimum(0.001);
  vector<PlotOpt> twodim_ops = {twodim_lin_lumi};

  //NamedFuncs
  //NamedFunc pTt("pTt",[](const Baby &b) -> NamedFunc::ScalarType{
  //  TVector3 g = AssignGamma(b).Vect();
  //  TVector3 h = AssignH(b).Vect();
  //  TVector3 z = AssignZ(b).Vect();
  //  g.SetZ(0); h.SetZ(0); z.SetZ(0);
  //  return h.Cross((z-g).Unit()).Mag();
  //});
  //NamedFunc pTt2("pTt2",[](const Baby &b) -> NamedFunc::ScalarType{
  //  TVector3 g = AssignGamma(b).Vect();
  //  TVector3 h = AssignH(b).Vect();
  //  TVector3 z = AssignZ(b).Vect();
  //  g.SetZ(0); h.SetZ(0); z.SetZ(0);
  //  return h.Cross(z-g).Mag()/h.Mag();
  //});
  //NamedFunc el_check("el_check",[](const Baby &b) -> NamedFunc::ScalarType{
  //  if(abs(b.el_eta()->at(b.ll_i1()->at(0)) - b.el_eta()->at(b.ll_i2()->at(0))) < 0.075 &&
  //     abs(b.el_phi()->at(b.ll_i1()->at(0)) - b.el_phi()->at(b.ll_i2()->at(0))) < 0.125)
  //     return false;
  //  return true;
  //});
  //NamedFunc relpT("p_{T}^{#gamma}/m_{Z#gamma}",[](const Baby &b) -> NamedFunc::ScalarType{
  //  if(b.nphoton() == 1) return b.photon_pt()->at(0)/b.llphoton_m()->at(0);
  //  double relp(-1);
  //  if(b.photon_pt()->at(1) > b.photon_pt()->at(0)) 
  //    if(b.photon_drmin()->at(1) > 0.3) 
  //      for(int i = 0; i < b.nllphoton(); i++)
  //        if(b.llphoton_iph()->at(i) == 1) {
  //          relp = b.photon_pt()->at(1)/b.llphoton_m()->at(i);
  //          break;
  //        }
  //  if(relp < 0) relp = b.photon_pt()->at(0)/b.llphoton_m()->at(0);
  //  return relp;
  //});
  //NamedFunc llp_m("llp_m",[](const Baby &b) -> NamedFunc::ScalarType{
  //  if(b.nphoton() == 1) return b.llphoton_m()->at(0);
  //  double m(-1);
  //  if(b.photon_pt()->at(1) > b.photon_pt()->at(0)) 
  //    if(b.photon_drmin()->at(1) > 0.3) 
  //      for(size_t i = 0; i < b.llphoton_pt()->size(); i++)
  //        if(b.llphoton_iph()->at(i) == 1) {
  //          m = b.llphoton_m()->at(i);
  //          break;
  //        }
  //  if(m < 0) m = b.llphoton_m()->at(0);
  //  return m;
  //});
  //NamedFunc els("electron channel",[](const Baby &b) -> NamedFunc::ScalarType{
  //  bool e =  b.ll_lepid()->at(0) == 11 && b.el_pt()->at(b.ll_i1()->at(0)) > 26 && b.el_pt()->at(b.ll_i2()->at(0)) > 17;
  //  return e;
  //});
  //NamedFunc mus("muon channel",[](const Baby &b) -> NamedFunc::ScalarType{
  //  bool mu = b.ll_lepid()->at(0) == 13 && b.mu_pt()->at(b.ll_i1()->at(0)) > 26 && b.mu_pt()->at(b.ll_i2()->at(0)) > 8;
  //  return mu;
  //});
  //NamedFunc vbf("VBF cut",[](const Baby &b) -> NamedFunc::ScalarType{
  //  return b.vbf_mva() > -0.99;
  //});
  //NamedFunc leps("e-or-#mu",[](const Baby &b) -> NamedFunc::ScalarType{
  //  bool e =  b.ll_lepid()->at(0) == 11 && b.el_pt()->at(b.ll_i1()->at(0)) > 26 && b.el_pt()->at(b.ll_i2()->at(0)) > 17;
  //  bool mu = b.ll_lepid()->at(0) == 13 && b.mu_pt()->at(b.ll_i1()->at(0)) > 26 && b.mu_pt()->at(b.ll_i2()->at(0)) > 8;
  //  return e || mu;
  //});
  //NamedFunc wgt("weight",[](const Baby &b) -> NamedFunc::ScalarType{ 
  //  double weight = b.w_lumi();
  //  if(b.type() >= 200000 && b.type() <= 205000)
  //    return weight;
  //  else if(b.type() == 6200)
  //    return 0.823*weight/1.664;
  //  return   0.823*weight;
  //});
  
  NamedFunc mc_higgs_pt("mc_higgs_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
      if (b.mc_id()->at(imc)==25) {
        return b.mc_pt()->at(imc);
      }
    }
    return -999;
  });

  NamedFunc mc_photon_pt("mc_photon_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
      if (b.mc_id()->at(imc)==22 && b.mc_mom()->at(imc)==25) {
        return b.mc_pt()->at(imc);
      }
    }
    return -999;
  });

  NamedFunc mc_deltar_zg("mc_deltar_zg",[](const Baby &b) -> NamedFunc::ScalarType{
    TVector3 ph;
    TVector3 Z;
    for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
      if (b.mc_id()->at(imc)==22 && b.mc_mom()->at(imc)==25) {
        ph.SetPtEtaPhi(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc));
      }
      if (b.mc_id()->at(imc)==23 && b.mc_mom()->at(imc)==25) {
        Z.SetPtEtaPhi(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc));
      }
    }
    return ph.DeltaR(Z);
  });

  NamedFunc photon_drmax("photon_drmax",[](const Baby &b) -> NamedFunc::ScalarType{
      return pdrmax(b);
  });

  NamedFunc mc_llg_in_acceptance("mc_llg_in_acceptance",[](const Baby &b) -> NamedFunc::ScalarType{
    for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
      if (b.mc_id()->at(imc)==22 && b.mc_mom()->at(imc)==25) {
        if (!(b.mc_pt()->at(imc)>15 && abs(b.mc_eta()->at(imc))<2.5))
          return 0;
      }
      if ((abs(b.mc_id()->at(imc))==11 || abs(b.mc_id()->at(imc))==13)&&b.mc_mom()->at(imc)==23) {
        int mom_idx = b.mc_momidx()->at(imc);
        if (mom_idx != -1) {
          if (b.mc_mom()->at(mom_idx)==25) {
            float eta_cut = 2.5;
            float pt_cut = 15;
            if (abs(b.mc_id()->at(imc))==13) {
              eta_cut = 2.4;
              pt_cut = 10;
            }
            if (!(b.mc_pt()->at(imc)>pt_cut && abs(b.mc_eta()->at(imc))<eta_cut))
              return 0;
          }
        }
      }
    }
    return 1;
  });

  NamedFunc mc_llg_reco("mc_llg_reco",[](const Baby &b) -> NamedFunc::ScalarType{
    TVector3 l1, l2, ph;
    TVector3 reco_obj;
    int z_decay_pdgid = 11;
    int lep_idx = 0;
    bool found_ph(false), found_l1(false), found_l2(false);
    for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
      if (b.mc_id()->at(imc)==22 && b.mc_mom()->at(imc)==25) {
        ph.SetPtEtaPhi(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc));
      }
      if ((abs(b.mc_id()->at(imc))==11 || abs(b.mc_id()->at(imc))==13)&&b.mc_mom()->at(imc)==23) {
        int mom_idx = b.mc_momidx()->at(imc);
        if (mom_idx != -1) {
          if (b.mc_mom()->at(mom_idx)==25) {
            z_decay_pdgid = abs(b.mc_id()->at(imc));
            if (lep_idx==0) {
              l1.SetPtEtaPhi(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc));
              lep_idx = 1;
            }
            else {
              l2.SetPtEtaPhi(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc));
            }
          }
        }
      }
    }
    for (unsigned iph = 0; iph < b.photon_pt()->size(); iph++) {
      reco_obj.SetPtEtaPhi(b.photon_pt()->at(iph),b.photon_eta()->at(iph),b.photon_phi()->at(iph));
      if (reco_obj.DeltaR(ph)<0.2) {
        found_ph = true;
      }
    }
    if (z_decay_pdgid==13) {
      for (unsigned imu = 0; imu < b.mu_pt()->size(); imu++) {
        reco_obj.SetPtEtaPhi(b.mu_pt()->at(imu),b.mu_eta()->at(imu),b.mu_phi()->at(imu));
        if (reco_obj.DeltaR(l1)<0.2) {
          found_l1 = true;
        }
        else if (reco_obj.DeltaR(l2)<0.2) {
          found_l2 = true;
        }
      }
    }
    else {
      for (unsigned iel = 0; iel < b.el_pt()->size(); iel++) {
        reco_obj.SetPtEtaPhi(b.el_pt()->at(iel),b.el_eta()->at(iel),b.el_phi()->at(iel));
        if (reco_obj.DeltaR(l1)<0.2) {
          found_l1 = true;
        }
        else if (reco_obj.DeltaR(l2)<0.2) {
          found_l2 = true;
        }
      }
    }
    if (found_ph&&found_l1&&found_l2) return 1;
    return 0;
  });

  NamedFunc nphoton_nodr("nphoton_nodr",[](const Baby &b) -> NamedFunc::ScalarType{
    double nphoton = 0;
    for (unsigned iph = 0; iph < b.photon_pt()->size(); iph++) {
      if (b.photon_pt()->at(iph)<15) continue;
      if (abs(b.photon_eta()->at(iph))>2.5) continue;
      if (abs(b.photon_eta()->at(iph))<1.4442 && b.photon_idmva()->at(iph) < -0.4) continue;
      if (abs(b.photon_eta()->at(iph))>1.4442 && abs(b.photon_eta()->at(iph))<1.566) continue;
      if (abs(b.photon_eta()->at(iph))>1.566 && b.photon_idmva()->at(iph) < -0.58) continue;
      if (!b.photon_elveto()->at(iph)) continue;
      nphoton += 1;
    }
    return nphoton;
  });

  NamedFunc sigphoton_mindr("sigphoton_mindr",[](const Baby &b) -> NamedFunc::ScalarType{
    for (unsigned iph = 0; iph < b.photon_pt()->size(); iph++) {
      if (b.photon_pt()->at(iph)<15) continue;
      if (abs(b.photon_eta()->at(iph))>2.5) continue;
      if (abs(b.photon_eta()->at(iph))<1.4442 && b.photon_idmva()->at(iph) < -0.4) continue;
      if (abs(b.photon_eta()->at(iph))>1.4442 && abs(b.photon_eta()->at(iph))<1.566) continue;
      if (abs(b.photon_eta()->at(iph))>1.566 && b.photon_idmva()->at(iph) < -0.58) continue;
      if (!b.photon_elveto()->at(iph)) continue;
      //use first signal photon
      TVector3 ph, lep;
      ph.SetPtEtaPhi(b.photon_pt()->at(iph),b.photon_eta()->at(iph),b.photon_phi()->at(iph));
      float mindr = 999;
      for (unsigned iel = 0; iel < b.el_pt()->size(); iel++) {
        if (b.el_sig()->at(iel)) {
          lep.SetPtEtaPhi(b.el_pt()->at(iel),b.el_eta()->at(iel),b.el_phi()->at(iel));
          if (ph.DeltaR(lep)<mindr) mindr = ph.DeltaR(lep);
        }
      }
      for (unsigned imu = 0; imu < b.mu_pt()->size(); imu++) {
        if (b.mu_sig()->at(imu)) {
          lep.SetPtEtaPhi(b.mu_pt()->at(imu),b.mu_eta()->at(imu),b.mu_phi()->at(imu));
          if (ph.DeltaR(lep)<mindr) mindr = ph.DeltaR(lep);
        }
      }
      return mindr;
    }
    return 999;
  });

  NamedFunc coslowertheta("costheta",[](const Baby &b) -> NamedFunc::ScalarType{
    return cos_theta(b);
  });

  NamedFunc coscaptheta("coscaptheta",[](const Baby &b) -> NamedFunc::ScalarType{
    return cos_Theta(b);
  });

  NamedFunc llgphi("phi",[](const Baby &b) -> NamedFunc::ScalarType{
    return Getphi(b);
  });

  NamedFunc wgt("weight",[](const Baby &b) -> NamedFunc::ScalarType{
    double w_lumi = b.w_lumi();
    //fix mis-weighted assoc. prod
    if(b.type() >= 202000 && b.type() <= 205000)
      w_lumi /= 0.100974;
    return w_lumi;
  });

  NamedFunc wgt_extrastats("wgt_extrastats",[](const Baby &b) -> NamedFunc::ScalarType{
    double w_lumi = b.w_lumi();
    //fix mis-weighted assoc. prod
    if(b.type() >= 202000 && b.type() <= 205000)
      w_lumi /= 0.100974;
    else if (b.type() < 200000 || b.type() > 205000)
      w_lumi /= 3.0;
    return w_lumi;
  });

  NamedFunc wgt_sigx20("weight",[](const Baby &b) -> NamedFunc::ScalarType{
    double w_lumi = b.w_lumi();
    //fix mis-weighted assoc. prod
    if(b.type() >= 202000 && b.type() <= 205000)
      w_lumi /= 0.100974;
    if(b.type() >= 200000 && b.type() <= 205000)
      w_lumi *= 20.0;
    return w_lumi;
  });

  NamedFunc higgs_pt("higgs_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    if(b.type() >= 200000 && b.type() <= 205000)
      return b.w_lumi();
    return b.w_lumi();
  });

  NamedFunc decorr_photon_pt("decorr_photon_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.photon_pt()->size()==0) return 0;
    if (b.llphoton_m()->size()==0) return 0;
    //return (b.photon_pt()->at(0))-0.202*(b.llphoton_m()->at(0));
    return (b.photon_pt()->at(0))-0.207*(b.llphoton_m()->at(0));
  });

  NamedFunc photon_pt_mass("photon_pt_mass",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.photon_pt()->size()==0) return 0;
    if (b.llphoton_m()->size()==0) return 0;
    //return (b.photon_pt()->at(0))-0.202*(b.llphoton_m()->at(0));
    return (b.photon_pt()->at(0))/(b.llphoton_m()->at(0));
  });

  NamedFunc truth_mllg("truth_mllg",[](const Baby &b) -> NamedFunc::ScalarType{
    //calculate truth mllg for dy (pre-stitch) events with ISR
    int lep_num = 0;
    float max_ph_pt = 0;
    TLorentzVector lep1, lep2, ph;
    for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
      if (abs(b.mc_id()->at(imc))==11) {
        if (lep_num==0)
          lep1.SetPtEtaPhiM(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc),0.000511);
        else
          lep2.SetPtEtaPhiM(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc),0.000511);
        lep_num++;
      }
      if (abs(b.mc_id()->at(imc))==13) {
        if (lep_num==0)
          lep1.SetPtEtaPhiM(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc),0.106);
        else
          lep2.SetPtEtaPhiM(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc),0.106);
        lep_num++;
      }
      if (b.mc_id()->at(imc)==22) {
        if (b.mc_pt()->at(imc) > max_ph_pt) {
          max_ph_pt = b.mc_pt()->at(imc);
          ph.SetPtEtaPhiM(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc),0.0);
        }
      }
    }
    if (lep_num<2 || max_ph_pt <= 0) return -999;
    return (lep1+lep2+ph).M();
  });

  NamedFunc truth_photon_pt("truth_photon_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    float max_ph_pt = 0;
    for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
      if (b.mc_id()->at(imc)==22) {
        if (b.mc_pt()->at(imc) > max_ph_pt) {
          max_ph_pt = b.mc_pt()->at(imc);
        }
      }
    }
    return max_ph_pt;
  });

  NamedFunc truth_photon_abseta("truth_photon_abseta",[](const Baby &b) -> NamedFunc::ScalarType{
    float max_ph_pt = 0;
    float ph_abseta = 0;
    for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
      if (b.mc_id()->at(imc)==22) {
        if (b.mc_pt()->at(imc) > max_ph_pt) {
          max_ph_pt = b.mc_pt()->at(imc);
          ph_abseta = abs(b.mc_eta()->at(imc)); 
        }
      }
    }
    return ph_abseta;
  });

  NamedFunc stitch_zg("stitch_zg",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type()==6200) return b.stitch_dy();
    return 1;
  });

  NamedFunc truth_lepton_minpt("truth_lepton_minpt",[](const Baby &b) -> NamedFunc::ScalarType{
    float min_lep_pt = 999;
    for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
      if (abs(b.mc_id()->at(imc))==11||abs(b.mc_id()->at(imc))==13) {
        if (b.mc_pt()->at(imc) < min_lep_pt) {
          min_lep_pt = b.mc_pt()->at(imc);
        }
      }
    }
    return min_lep_pt;
  });

  NamedFunc truth_lepton_eta_acc("truth_lepton_eta_acc",[](const Baby &b) -> NamedFunc::ScalarType{
    for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
      if (abs(b.mc_id()->at(imc))==11||abs(b.mc_id()->at(imc))==13) {
        if (abs(b.mc_eta()->at(imc))>2.4)
          return 0.0;
      }
    }
    return 1.0;
  });

  NamedFunc truth_min_dr("truth_min_dr",[](const Baby &b) -> NamedFunc::ScalarType{
    float ph_eta = -999;
    float ph_phi = 0;
    float max_ph_pt = 0;
    float min_dr = 999;
    float dr = 999;
    for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
      if (b.mc_id()->at(imc)==22) {
        if (b.mc_pt()->at(imc) > max_ph_pt) {
          max_ph_pt = b.mc_pt()->at(imc);
          ph_eta = b.mc_eta()->at(imc);
          ph_phi = b.mc_phi()->at(imc);
        }
      }
    }
    for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
      if (abs(b.mc_id()->at(imc))==11||abs(b.mc_id()->at(imc))==13) {
        dr = deltaR(ph_eta, ph_phi, b.mc_eta()->at(imc), b.mc_phi()->at(imc));
        if (dr < min_dr)
          min_dr = dr;
      }
    }
    return min_dr;
  });

  NamedFunc abs_cos_theta("abs_cos_theta",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.llphoton_cosTheta()->size()==0) return -999;
    return abs(b.llphoton_cosTheta()->at(0));
  });

  NamedFunc photon_res = "photon_pterr[0]/photon_pt[0]";
  NamedFunc pt_mass = "llphoton_pt[0]/llphoton_m[0]";

  NamedFunc baseline_nolep = "(ll_m[0]>50) && (photon_pt[0]/llphoton_m[0]>=15.0/110.0) && ((llphoton_m[0]+ll_m[0])>=185) && (photon_drmin[0]>0.4)";
  NamedFunc baseline_el = "(ll_lepid[0]==11) && (el_pt[ll_i1[0]]>25) && (el_pt[ll_i2[0]]>15)" && baseline_nolep;
  NamedFunc baseline_mu = "(ll_lepid[0]==13) && (mu_pt[ll_i1[0]]>20) && (mu_pt[ll_i2[0]]>10)" && baseline_nolep;
  baseline_el.Name("electron_baseline");
  baseline_mu.Name("muon_baseline");
  NamedFunc baseline = (baseline_el||baseline_mu);
  baseline.Name("baseline");
  std::vector<NamedFunc> baseline_lep = {baseline_el, baseline_mu, baseline_el||baseline_mu};

  NamedFunc baseline_nolep_nomll = "(photon_pt[0]/llphoton_m[0]>=15.0/110.0) && (photon_drmin[0]>0.4)";
  NamedFunc baseline_el_nomll = "(ll_lepid[0]==11) && (el_pt[ll_i1[0]]>25) && (el_pt[ll_i2[0]]>15)" && baseline_nolep_nomll;
  NamedFunc baseline_mu_nomll = "(ll_lepid[0]==13) && (mu_pt[ll_i1[0]]>20) && (mu_pt[ll_i2[0]]>10)" && baseline_nolep_nomll;
  baseline_el_nomll.Name("electron_baseline");
  baseline_mu_nomll.Name("muon_baseline");
  NamedFunc baseline_nomll = (baseline_el_nomll||baseline_mu_nomll);
  baseline_nomll.Name("baseline (no m_{ll})");

  NamedFunc baseline_nolep_noph = "(ll_m[0]>50) && ((llphoton_m[0]+ll_m[0])>=185) && (photon_drmin[0]>0.4)";
  NamedFunc baseline_el_noph = "(ll_lepid[0]==11) && (el_pt[ll_i1[0]]>25) && (el_pt[ll_i2[0]]>15)" && baseline_nolep_noph;
  NamedFunc baseline_mu_noph = "(ll_lepid[0]==13) && (mu_pt[ll_i1[0]]>20) && (mu_pt[ll_i2[0]]>10)" && baseline_nolep_noph;
  baseline_el_noph.Name("electron_baseline");
  baseline_mu_noph.Name("muon_baseline");
  NamedFunc baseline_noph = (baseline_el_noph||baseline_mu_noph);
  baseline_noph.Name("baseline (no photon p_{T})");

  //sanitize separate
  NamedFunc baseline_debug_training = 
  sanitize&&"(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL||HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ||HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ)&&(stitch_dy||(type/1000>=6&&type/1000<7))&&((ll_m[0]>50)&&(photon_pt[0]/llphoton_m[0]>=15.0/110.0)&&((llphoton_m[0]+ll_m[0])>185)&&(photon_drmin[0]>0.4))&&(ll_lepid[0]==11&&el_pt[ll_i1[0]]>25&&el_pt[ll_i2[0]]>15)||(ll_lepid[0]==13&&mu_pt[ll_i1[0]]>20&&mu_pt[ll_i2[0]]>10)&&(photon_idmva[0]>0.5)&&(llphoton_m[0]>120&&llphoton_m[0]<130)";

  NamedFunc higgs_window = "llphoton_m[0]>122&&llphoton_m[0]<128";
  NamedFunc photon_pt = "photon_pt[0]";
  std::vector<float> photon_pt_bins = {15,20,25,30,40,50,200};

  //note skim_llg contains nlep>1 && nphoton>0 (delta r in nphoton?)
  //baseline include mZ>50 leadleppt>25(e)20(mu) ptgam/mllg>0.14 mll+mllg>185
  PlotMaker pm;

  bool plot_truth = false;
  bool plot_boostbdt = false;
  bool plot_kinematics = false;
  bool plot_basic = false;
  bool plot_corr = false;
  bool plot_yieldtable = false;
  bool plot_decorrelate = false;
  bool plot_photonpt_story = false;
  bool plot_ggf_general = false;
  bool plot_photonpt_compare = false;
  bool plot_mva_debug = false;
  bool plot_final_bdt_mllg = true;
  //bool plot_fitshapes = false;

  if (plot_truth) {
    pm.Push<Hist1D>(Axis(12,40,400, "mc_pt", "Higgs p_{T} [GeV]", {}), 
        "mc_id==25", procs_sig_unskimmed, ops).Weight(wgt).Tag("zgboost_acc");
    pm.Push<Hist2D>(Axis(20, 0, 400, mc_higgs_pt, "Higgs p_{T} [GeV]", {}),
        Axis(24, 0, 120.0, mc_photon_pt, "Photon p_{T} [GeV]", {}),
        "1", procs_sig_unskimmed, twodim_ops).Weight(wgt).Tag("zgboost_acc");
    pm.Push<Hist2D>(Axis(20, 0, 400, mc_higgs_pt, "Higgs p_{T} [GeV]", {}),
        Axis(20, 0, 2.0, mc_deltar_zg, "#Delta R(Z #gamma)", {}),
        "1", procs_sig_unskimmed, twodim_ops).Weight(wgt).Tag("zgboost_acc");
    pm.Push<Hist1D>(Axis(20,0,2.0, mc_deltar_zg, "#Delta R(Z,#gamma)", {}), 
        mc_higgs_pt>160, procs_sig_unskimmed, ops).Weight(wgt).Tag("zgboost_acc");
    pm.Push<Hist1D>(Axis(20,0,200, mc_photon_pt, "Photon p_{T} [GeV]", {}), 
        mc_higgs_pt>160, procs_sig_unskimmed, ops).Weight(wgt).Tag("zgboost_acc");
    pm.Push<Hist1D>(Axis(40,0,400, mc_higgs_pt, "Higgs p_{T} [GeV]", {}), 
        mc_photon_pt>50, procs_sig_unskimmed, ops).Weight(wgt).Tag("zgboost_acc");
    pm.Push<Table>("truth_cutflow", vector<TableRow>{
      TableRow("Inclusive", 
          "1",0,0,wgt),
      TableRow("High p_{T} Higgs", 
          mc_higgs_pt>160,0,0,wgt),
      TableRow("Objects in acceptance", 
          mc_higgs_pt>160&&mc_llg_in_acceptance,0,0,wgt),
      TableRow("Objects reconstructed", 
          mc_higgs_pt>160&&mc_llg_in_acceptance&&mc_llg_reco,0,0,wgt),
      TableRow("Objects pass ID", 
          mc_higgs_pt>160&&mc_llg_in_acceptance&&mc_llg_reco&&nphoton_nodr>=1&&"nlep>=2",0,0,wgt),
      TableRow("Min $\\Delta R(l,\\gamma)>0.4$", 
          mc_higgs_pt>160&&mc_llg_in_acceptance&&mc_llg_reco&&nphoton_nodr>=1&&"nlep>=2"&&sigphoton_mindr>0.4,0,0,wgt),
    },procs_sig_unskimmed,false,true,false,false,false,true).Precision(3);

  }

  std::vector<float> boost_bins = {-1.0,-0.10,0.05,1.0};
  if (plot_boostbdt) {
    for (unsigned iscale = 0; iscale < 2; iscale++) {
      vector<shared_ptr<Process>> plot_procs = procs;
      NamedFunc plot_weight = wgt;
      std::string base_tag_string = "zgboost_bdt_noscale";
      if (iscale==1) {
        plot_procs = procs_sigx20;
        plot_weight = wgt_sigx20;
        base_tag_string = "zgboost_bdt_signalx20";
      }
      for (unsigned ilep = 0; ilep < 3; ilep++) {
        pm.Push<Hist1D>(Axis(40,-1.0,1.0, BoostScore, "Boosted ll#gamma BDT Score", {}), 
            baseline_lep[ilep]&&"photon_idmva[0]>0.5", plot_procs, ops_sorb).Weight(plot_weight).Tag(base_tag_string);
        pm.Push<Hist1D>(Axis(40,-1.0,1.0, BoostScore, "Boosted ll#gamma BDT Score", {}), 
            baseline_lep[ilep]&&"photon_idmva[0]>0.5"&&higgs_window, plot_procs, ops_sorb).Weight(plot_weight).Tag(base_tag_string);
      }
      //for (unsigned iboost = 0; iboost < 3; iboost++) {
      //  pm.Push<Hist1D>(Axis(30,100,160, "llphoton_m[0]", "Higgs Candidate m [GeV]", {}), 
      //      baseline_lep[2]&&"photon_idmva[0]>0.5"&&BoostScore>boost_bins[iboost]&&BoostScore<boost_bins[iboost+1], plot_procs, ops).Weight(plot_weight).Tag(base_tag_string);
      //  pm.Push<Hist1D>(Axis(30,0,240, "llphoton_pt[0]", "Higgs Candidate p_{T} [GeV]", {}), 
      //      baseline_lep[2]&&"photon_idmva[0]>0.5"&&BoostScore>boost_bins[iboost]&&BoostScore<boost_bins[iboost+1], plot_procs, ops).Weight(plot_weight).Tag(base_tag_string);
      //  pm.Push<Hist1D>(Axis(40,0,4.0, "llphoton_dr[0]", "#Delta R(Z,#gamma)", {}), 
      //      baseline_lep[2]&&"photon_idmva[0]>0.5"&&BoostScore>boost_bins[iboost]&&BoostScore<boost_bins[iboost+1], plot_procs, ops).Weight(plot_weight).Tag(base_tag_string);
      //  pm.Push<Hist1D>(Axis(30,15,75, "photon_pt[0]", "Photon p_{T} [GeV]", {}), 
      //      baseline_lep[2]&&"photon_idmva[0]>0.5"&&BoostScore>boost_bins[iboost]&&BoostScore<boost_bins[iboost+1], plot_procs, ops).Weight(plot_weight).Tag(base_tag_string);
      //  pm.Push<Hist1D>(Axis(30,0,210, "ll_pt[0]", "Z Candidate p_{T} [GeV]", {}), 
      //      baseline_lep[2]&&"photon_idmva[0]>0.5"&&BoostScore>boost_bins[iboost]&&BoostScore<boost_bins[iboost+1], plot_procs, ops).Weight(plot_weight).Tag(base_tag_string);
      //}
      //pm.Push<Hist1D>(Axis(12,100,160, "llphoton_m[0]", "Higgs Candidate m [GeV]", {}), 
      //    baseline_lep[2]&&"photon_idmva[0]>0.5"&&BoostScore>0.35, plot_procs, ops).Weight(plot_weight).Tag(base_tag_string);
      //pm.Push<Hist1D>(Axis(30,0,240, "llphoton_pt[0]", "Higgs Candidate p_{T} [GeV]", {}), 
      //    baseline_lep[2]&&"photon_idmva[0]>0.5"&&BoostScore>0.35, plot_procs, ops).Weight(plot_weight).Tag(base_tag_string);
      //pm.Push<Hist1D>(Axis(40,0,4.0, "llphoton_dr[0]", "#Delta R(Z,#gamma)", {}), 
      //    baseline_lep[2]&&"photon_idmva[0]>0.5"&&BoostScore>0.35, plot_procs, ops).Weight(plot_weight).Tag(base_tag_string);
      //pm.Push<Hist1D>(Axis(30,15,75, "photon_pt[0]", "Photon p_{T} [GeV]", {}), 
      //    baseline_lep[2]&&"photon_idmva[0]>0.5"&&BoostScore>0.35, plot_procs, ops).Weight(plot_weight).Tag(base_tag_string);
      //pm.Push<Hist1D>(Axis(30,0,210, "ll_pt[0]", "Z Candidate p_{T} [GeV]", {}), 
      //    baseline_lep[2]&&"photon_idmva[0]>0.5"&&BoostScore>0.35, plot_procs, ops).Weight(plot_weight).Tag(base_tag_string);
    }
  }

  if (plot_kinematics) {
    for (unsigned iscale = 0; iscale < 2; iscale++) {
      vector<shared_ptr<Process>> plot_procs = procs;
      NamedFunc plot_weight = wgt;
      std::string base_tag_string = "zgboost_kin_noscale";
      if (iscale==1) {
        plot_procs = procs_sigx20;
        plot_weight = wgt_sigx20;
        base_tag_string = "zgboost_kin_signalx20";
      }

      pm.Push<Hist1D>(Axis(40,-1.0,1.0, "photon_idmva[0]", "Photon IDMVA", {}), 
          baseline_lep[2]&&higgs_window, plot_procs, ops_sorb).Weight(plot_weight).Tag(base_tag_string);
      pm.Push<Hist1D>(Axis(30,0,240, "llphoton_pt[0]", "Higgs Candidate p_{T} [GeV]", {}), 
          baseline_lep[2]&&higgs_window&&"photon_idmva[0]>0.5", plot_procs, ops_sorb).Weight(plot_weight).Tag(base_tag_string);
      pm.Push<Hist1D>(Axis(40,0,4.0, "llphoton_dr[0]", "#Delta R(Z,#gamma)", {}), 
          baseline_lep[2]&&higgs_window&&"photon_idmva[0]>0.5", procs_extrastats, ops_sorb_upper).Weight(wgt_extrastats).Tag(base_tag_string);
      pm.Push<Hist1D>(Axis(30,15,75, "photon_pt[0]", "Photon p_{T} [GeV]", {}), 
          baseline_lep[2]&&higgs_window&&"photon_idmva[0]>0.5", plot_procs, ops_sorb).Weight(plot_weight).Tag(base_tag_string);
      pm.Push<Hist1D>(Axis(30,0,210, "ll_pt[0]", "Z Candidate p_{T} [GeV]", {}), 
          baseline_lep[2]&&higgs_window&&"photon_idmva[0]>0.5", plot_procs, ops_sorb).Weight(plot_weight).Tag(base_tag_string);

      pm.Push<Hist1D>(Axis(20,100,180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
          baseline_lep[2]&&"llphoton_pt[0]>150", plot_procs, ops).Weight(plot_weight).Tag(base_tag_string);
      pm.Push<Hist1D>(Axis(20,100,180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
          baseline_lep[2]&&"photon_pt[0]>45", plot_procs, ops).Weight(plot_weight).Tag(base_tag_string);
      pm.Push<Hist1D>(Axis(20,100,180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
          baseline_lep[2]&&"llphoton_dr[0]<0.8", plot_procs, ops).Weight(plot_weight).Tag(base_tag_string);
      pm.Push<Hist1D>(Axis(20,100,180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
          baseline_lep[2]&&"ll_pt[0]>125", plot_procs, ops).Weight(plot_weight).Tag(base_tag_string);
      pm.Push<Hist1D>(Axis(20,100,180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
          baseline_lep[2]&&"llphoton_pt[0]/llphoton_m[0]>1.2", plot_procs, ops).Weight(plot_weight).Tag(base_tag_string);
      pm.Push<Hist1D>(Axis(20,100,180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
          baseline_lep[2]&&"photon_pt[0]/llphoton_m[0]>0.36", plot_procs, ops).Weight(plot_weight).Tag(base_tag_string);
      pm.Push<Hist1D>(Axis(20,100,180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
          baseline_lep[2]&&"ll_pt[0]/llphoton_m[0]>1.0", plot_procs, ops).Weight(plot_weight).Tag(base_tag_string);

      for (unsigned in_signal_window = 0; in_signal_window < 4; in_signal_window++) {
        NamedFunc additional_cuts = "1";
        std::string tag_string = base_tag_string;
        if (in_signal_window==1) {
          additional_cuts = higgs_window;
          tag_string += "_higgswindow";
        }
        else if (in_signal_window==2) {
          additional_cuts = "photon_idmva[0]>0.6";
          tag_string += "_highphmva";
        }
        else if (in_signal_window==3) {
          additional_cuts = higgs_window&&"photon_idmva[0]>0.6";
          tag_string += "_highphmva_higgswindow";
        }

        for (unsigned ilep = 0; ilep < 3; ilep++) {
          //pt variables
          pm.Push<Hist1D>(Axis(16,0,160, "photon_pt[0]", "Photon p_{T} [GeV]", {}), 
              additional_cuts&&baseline_lep[ilep]&&"llphoton_pt[0]>160", plot_procs, ops).Weight(plot_weight).Tag(tag_string);
          pm.Push<Hist1D>(Axis(16,0,240, "lep_pt[0]", "Lead Lepton p_{T} [GeV]", {}), 
              additional_cuts&&baseline_lep[ilep]&&"llphoton_pt[0]>160", plot_procs, ops).Weight(plot_weight).Tag(tag_string);
          pm.Push<Hist1D>(Axis(16,0,160, "lep_pt[1]", "Sublead Lepton p_{T} [GeV]", {}), 
              additional_cuts&&baseline_lep[ilep]&&"llphoton_pt[0]>160", plot_procs, ops).Weight(plot_weight).Tag(tag_string);
          pm.Push<Hist1D>(Axis(12,160,400, "llphoton_pt[0]", "Higgs Candidate p_{T} [GeV]", {}), 
              additional_cuts&&baseline_lep[ilep]&&"llphoton_pt[0]>160", plot_procs, ops).Weight(plot_weight).Tag(tag_string);
          pm.Push<Hist1D>(Axis(12,0,240, "ll_pt[0]", "Z Candidate p_{T} [GeV]", {}), 
              additional_cuts&&baseline_lep[ilep]&&"llphoton_pt[0]>160", plot_procs, ops).Weight(plot_weight).Tag(tag_string);
          //eta variables
          pm.Push<Hist1D>(Axis(13,-2.6,2.6, "photon_eta[0]", "Photon #eta", {}), 
              additional_cuts&&baseline_lep[ilep]&&"llphoton_pt[0]>160", plot_procs, ops).Weight(plot_weight).Tag(tag_string);
          pm.Push<Hist1D>(Axis(13,-2.6,2.6, "lep_eta[0]", "Lead Lepton #eta", {}), 
              additional_cuts&&baseline_lep[ilep]&&"llphoton_pt[0]>160", plot_procs, ops).Weight(plot_weight).Tag(tag_string);
          pm.Push<Hist1D>(Axis(13,-2.6,2.6, "lep_eta[1]", "Sublead Lepton #eta", {}), 
              additional_cuts&&baseline_lep[ilep]&&"llphoton_pt[0]>160", plot_procs, ops).Weight(plot_weight).Tag(tag_string);
          //angular variables
          pm.Push<Hist1D>(Axis(10,0,4.0, "photon_drmin[0]", "Min #Delta R(l,#gamma)", {}), 
              additional_cuts&&baseline_lep[ilep]&&"llphoton_pt[0]>160", plot_procs, ops).Weight(plot_weight).Tag(tag_string);
          pm.Push<Hist1D>(Axis(10,0,4.0, photon_drmax, "Max #Delta R(l,#gamma)", {}), 
              additional_cuts&&baseline_lep[ilep]&&"llphoton_pt[0]>160", plot_procs, ops).Weight(plot_weight).Tag(tag_string);
          pm.Push<Hist1D>(Axis(10,0,4.0, "llphoton_dr[0]", "#Delta R(Z,#gamma)", {}), 
              additional_cuts&&baseline_lep[ilep]&&"llphoton_pt[0]>160", plot_procs, ops).Weight(plot_weight).Tag(tag_string);
          pm.Push<Hist1D>(Axis(10,-1.0,1.0, coslowertheta, "cos#theta", {}), 
              additional_cuts&&baseline_lep[ilep]&&"llphoton_pt[0]>160", plot_procs, ops).Weight(plot_weight).Tag(tag_string);
          pm.Push<Hist1D>(Axis(10,-1.0,1.0, coscaptheta, "cos#Theta", {}), 
              additional_cuts&&baseline_lep[ilep]&&"llphoton_pt[0]>160", plot_procs, ops).Weight(plot_weight).Tag(tag_string);
          pm.Push<Hist1D>(Axis(10,-3.15,3.15, llgphi, "#phi", {}), 
              additional_cuts&&baseline_lep[ilep]&&"llphoton_pt[0]>160", plot_procs, ops).Weight(plot_weight).Tag(tag_string);
          //photon quality variables
          pm.Push<Hist1D>(Axis(10,-1.0,1.0, "photon_idmva[0]", "Photon ID MVA", {}), 
              additional_cuts&&baseline_lep[ilep]&&"llphoton_pt[0]>160", plot_procs, ops).Weight(plot_weight).Tag(tag_string);
          pm.Push<Hist1D>(Axis(10,0.0,0.2, photon_res, "Photon Resolution", {}), 
              additional_cuts&&baseline_lep[ilep]&&"llphoton_pt[0]>160", plot_procs, ops).Weight(plot_weight).Tag(tag_string);
          //mass variables
          pm.Push<Hist1D>(Axis(20,80,160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
              additional_cuts&&baseline_lep[ilep]&&"llphoton_pt[0]>160", plot_procs, ops).Weight(plot_weight).Tag(tag_string);
          pm.Push<Hist1D>(Axis(15,50,125, "ll_m[0]", "m_{ll} [GeV]", {}), 
              additional_cuts&&baseline_lep[ilep]&&"llphoton_pt[0]>160", plot_procs, ops).Weight(plot_weight).Tag(tag_string);
        }
      }
    }

    pm.Push<Hist1D>(Axis(1,122,128, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        higgs_window&&"photon_idmva[0]>0.5"&&baseline_lep[2]&&"llphoton_pt[0]>160", procs, ops).Weight(wgt).Tag("zgboost_kin_yield");
    pm.Push<Hist1D>(Axis(10,0,3.0, "photon_drmin[0]", "Min #Delta R(l,#gamma)", {}), 
        higgs_window&&"photon_idmva[0]>0.5"&&baseline_lep[2]&&"llphoton_pt[0]>160", procs, ops).Weight(wgt).Tag("zgboost_kin_yield");
    pm.Push<Hist1D>(Axis(1,122,128, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        "photon_drmin[0]<1.0"&&higgs_window&&"photon_idmva[0]>0.5"&&baseline_lep[2]&&"llphoton_pt[0]>160", procs, ops).Weight(wgt).Tag("zgboost_kin_yield");
    pm.Push<Hist1D>(Axis(1,122,128, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        "photon_pt[0]>40&&photon_drmin[0]<1.0"&&higgs_window&&"photon_idmva[0]>0.5"&&baseline_lep[2]&&"llphoton_pt[0]>160", procs, ops).Weight(wgt).Tag("zgboost_kin_yield");
    pm.Push<Hist1D>(Axis(10,110,170, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        "photon_drmin[0]<1.0&&photon_idmva[0]>0.5"&&baseline_lep[2]&&"llphoton_pt[0]>160", procs, ops).Weight(wgt).Tag("zgboost_kin_yield");
    pm.Push<Hist1D>(Axis(10,110,170, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        "photon_pt[0]>40&&photon_drmin[0]<1.0&&photon_idmva[0]>0.5"&&baseline_lep[2]&&"llphoton_pt[0]>160", procs, ops).Weight(wgt).Tag("zgboost_kin_yield");
  }

  if (plot_basic) {
    for (unsigned iscale = 0; iscale < 2; iscale++) {
      vector<shared_ptr<Process>> plot_procs = procs;
      NamedFunc plot_weight = wgt;
      std::string tag_string = "zgboost_noscale";
      if (iscale==1) {
        plot_procs = procs_sigx20;
        plot_weight = wgt_sigx20;
        tag_string = "zgboost_signalx20";
      }

      for (unsigned ilep = 0; ilep < 2; ilep++) {
        pm.Push<Hist1D>(Axis(12,40,400, "llphoton_pt[0]", "Higgs candidate p_{T} [GeV]", {}), 
            baseline_lep[ilep], plot_procs, ops).Weight(plot_weight).Tag(tag_string);
        pm.Push<Hist1D>(Axis(12,40,400, "llphoton_pt[0]", "Higgs candidate p_{T} [GeV]", {}), 
            baseline_lep[ilep]&&higgs_window, plot_procs, ops).Weight(plot_weight).Tag(tag_string);
        pm.Push<Hist1D>(Axis(20,100,180, "llphoton_m[0]", "m_{Z#gamma} [GeV]", {}), 
            baseline_lep[ilep], plot_procs, ops).Weight(plot_weight).Tag(tag_string);
        pm.Push<Hist1D>(Axis(45,15,60, "photon_pt[0]", "Photon p_{T} [GeV]", {}), 
            baseline_lep[ilep]&&higgs_window, plot_procs, ops).Weight(plot_weight).Tag(tag_string+"_ext");
        pm.Push<Hist1D>(Axis(20,40,60, "photon_pt[0]", "Photon p_{T} [GeV]", {}), 
            baseline_lep[ilep]&&higgs_window, plot_procs, ops).Weight(plot_weight).Tag(tag_string);
        for (unsigned iphpt = 0; iphpt < (photon_pt_bins.size()-1); iphpt++) {
          pm.Push<Hist1D>(Axis(20,100,180, "llphoton_m[0]", "m_{Z#gamma} [GeV]", {}), 
              baseline_lep[ilep]&&(photon_pt>=photon_pt_bins[iphpt]&&photon_pt<photon_pt_bins[iphpt+1]), 
              plot_procs, ops).Weight(plot_weight).Tag(tag_string);
        }
      }

      //now with high mva score
      for (unsigned ilep = 0; ilep < 2; ilep++) {
        pm.Push<Hist1D>(Axis(12,40,400, "llphoton_pt[0]", "Higgs candidate p_{T} [GeV]", {}), 
            baseline_lep[ilep]&&MVAScore>0.683, plot_procs, ops).Weight(plot_weight).Tag(tag_string);
        pm.Push<Hist1D>(Axis(12,40,400, "llphoton_pt[0]", "Higgs candidate p_{T} [GeV]", {}), 
            baseline_lep[ilep]&&higgs_window&&MVAScore>0.683, plot_procs, ops).Weight(plot_weight).Tag(tag_string);
        pm.Push<Hist1D>(Axis(20,100,180, "llphoton_m[0]", "m_{Z#gamma} [GeV]", {}), 
            baseline_lep[ilep]&&MVAScore>0.683, plot_procs, ops).Weight(plot_weight).Tag(tag_string);
        pm.Push<Hist1D>(Axis(45,15,60, "photon_pt[0]", "Photon p_{T} [GeV]", {}), 
            baseline_lep[ilep]&&higgs_window&&MVAScore>0.683, plot_procs, ops).Weight(plot_weight).Tag(tag_string+"_ext");
        pm.Push<Hist1D>(Axis(20,40,60, "photon_pt[0]", "Photon p_{T} [GeV]", {}), 
            baseline_lep[ilep]&&higgs_window&&MVAScore>0.683, plot_procs, ops).Weight(plot_weight).Tag(tag_string);
        for (unsigned iphpt = 0; iphpt < (photon_pt_bins.size()-1); iphpt++) {
          pm.Push<Hist1D>(Axis(20,100,180, "llphoton_m[0]", "m_{Z#gamma} [GeV]", {}), 
              baseline_lep[ilep]&&MVAScore>0.683&&(photon_pt>=photon_pt_bins[iphpt]&&photon_pt<photon_pt_bins[iphpt+1]), 
              plot_procs, ops).Weight(plot_weight).Tag(tag_string);
        }
      }
    }
  }
  if (plot_corr) {
    //2D plots of boost proxy variables
    pm.Push<Hist2D>(Axis(30, 15, 75, "photon_pt[0]", "Photon p_{T} [GeV]", {}),
        Axis(40, 0, 4.0, "llphoton_dr[0]", "#Delta R(Z,#gamma)", {}),
        baseline, procs_sigonly, twodim_ops).Weight(wgt).Tag("zgboost_sig");
    pm.Push<Hist2D>(Axis(30, 15, 75, "photon_pt[0]", "Photon p_{T} [GeV]", {}),
        Axis(30, 0, 240, "llphoton_pt[0]", "Higgs Candidate p_{T} [GeV]", {}),
        baseline, procs_sigonly, twodim_ops).Weight(wgt).Tag("zgboost_sig");
    pm.Push<Hist2D>(Axis(30, 15, 75, "photon_pt[0]", "Photon p_{T} [GeV]", {}),
        Axis(30, 0, 210, "ll_pt[0]", "Z Candidate p_{T} [GeV]", {}),
        baseline, procs_sigonly, twodim_ops).Weight(wgt).Tag("zgboost_sig");
    pm.Push<Hist2D>(Axis(30, 0, 240, "llphoton_pt[0]", "Higgs Candidate p_{T} [GeV]", {}),
        Axis(40, 0, 4.0, "llphoton_dr[0]", "#Delta R(Z,#gamma)", {}),
        baseline, procs_sigonly, twodim_ops).Weight(wgt).Tag("zgboost_sig");
    pm.Push<Hist2D>(Axis(30, 0, 240, "llphoton_pt[0]", "Higgs Candidate p_{T} [GeV]", {}),
        Axis(30, 0, 210, "ll_pt[0]", "Z Candidate p_{T} [GeV]", {}),
        baseline, procs_sigonly, twodim_ops).Weight(wgt).Tag("zgboost_sig");
    pm.Push<Hist2D>(Axis(30, 0, 210, "ll_pt[0]", "Z Candidate p_{T} [GeV]", {}),
        Axis(40, 0, 4.0, "llphoton_dr[0]", "#Delta R(Z,#gamma)", {}),
        baseline, procs_sigonly, twodim_ops).Weight(wgt).Tag("zgboost_sig");

    pm.Push<Hist2D>(Axis(30, 15, 75, "photon_pt[0]", "Photon p_{T} [GeV]", {}),
        Axis(40, 0, 4.0, "llphoton_dr[0]", "#Delta R(Z,#gamma)", {}),
        baseline&&higgs_window, procs_bakonly, twodim_ops).Weight(wgt).Tag("zgboost_bak");
    pm.Push<Hist2D>(Axis(30, 15, 75, "photon_pt[0]", "Photon p_{T} [GeV]", {}),
        Axis(30, 0, 240, "llphoton_pt[0]", "Higgs Candidate p_{T} [GeV]", {}),
        baseline&&higgs_window, procs_bakonly, twodim_ops).Weight(wgt).Tag("zgboost_bak");
    pm.Push<Hist2D>(Axis(30, 15, 75, "photon_pt[0]", "Photon p_{T} [GeV]", {}),
        Axis(30, 0, 210, "ll_pt[0]", "Z Candidate p_{T} [GeV]", {}),
        baseline&&higgs_window, procs_bakonly, twodim_ops).Weight(wgt).Tag("zgboost_bak");
    pm.Push<Hist2D>(Axis(30, 0, 240, "llphoton_pt[0]", "Higgs Candidate p_{T} [GeV]", {}),
        Axis(40, 0, 4.0, "llphoton_dr[0]", "#Delta R(Z,#gamma)", {}),
        baseline&&higgs_window, procs_bakonly, twodim_ops).Weight(wgt).Tag("zgboost_bak");
    pm.Push<Hist2D>(Axis(30, 0, 240, "llphoton_pt[0]", "Higgs Candidate p_{T} [GeV]", {}),
        Axis(30, 0, 210, "ll_pt[0]", "Z Candidate p_{T} [GeV]", {}),
        baseline&&higgs_window, procs_bakonly, twodim_ops).Weight(wgt).Tag("zgboost_bak");
    pm.Push<Hist2D>(Axis(30, 0, 210, "ll_pt[0]", "Z Candidate p_{T} [GeV]", {}),
        Axis(40, 0, 4.0, "llphoton_dr[0]", "#Delta R(Z,#gamma)", {}),
        baseline&&higgs_window, procs_bakonly, twodim_ops).Weight(wgt).Tag("zgboost_bak");

    //old corr plots
    pm.Push<Hist2D>(Axis(45, 15, 60, "photon_pt[0]", "Photon p_{T} [GeV]", {}),
        Axis(36, 0.4, 4.0, "photon_drmin[0]", "Minimum #Delta R_{l #gamma}", {}),
        baseline, procs_sigonly, twodim_ops).Weight(wgt).Tag("zgboost_sig");
    pm.Push<Hist2D>(Axis(45, 15, 60, "photon_pt[0]", "Photon p_{T} [GeV]", {}),
        Axis(20, -1.0, 1.0, coslowertheta, "ll#gamma cos#theta", {}),
        baseline, procs_sigonly, twodim_ops).Weight(wgt).Tag("zgboost_sig");
    pm.Push<Hist2D>(Axis(45, 15, 60, "photon_pt[0]", "Photon p_{T} [GeV]", {}),
        Axis(20, -1.0, 1.0, coscaptheta, "cos#Theta", {}),
        baseline, procs_sigonly, twodim_ops).Weight(wgt).Tag("zgboost_sig");
    pm.Push<Hist2D>(Axis(45, 15, 60, "photon_pt[0]", "Photon p_{T} [GeV]", {}),
        Axis(20, -3.15, 3.15, llgphi, "#phi", {}),
        baseline, procs_sigonly, twodim_ops).Weight(wgt).Tag("zgboost_sig");
    pm.Push<Hist2D>(Axis(26, 40, 300, "llphoton_pt[0]", "Higgs Candidate p_{T} [GeV]", {}),
        Axis(20, 100, 160, "llphoton_m[0]", "Higgs Candidate m [GeV]", {}),
        baseline, procs_sigonly, twodim_ops).Weight(wgt).Tag("zgboost_sig");

    pm.Push<Hist2D>(Axis(45, 15, 60, "photon_pt[0]", "Photon p_{T} [GeV]", {}),
        Axis(36, 0.4, 4.0, "photon_drmin[0]", "Minimum #Delta R_{l #gamma}", {}),
        baseline, procs_bakonly, twodim_ops).Weight(wgt).Tag("zgboost_bak");
    pm.Push<Hist2D>(Axis(45, 15, 60, "photon_pt[0]", "Photon p_{T} [GeV]", {}),
        Axis(20, -1.0, 1.0, coslowertheta, "ll#gamma cos#theta", {}),
        baseline, procs_bakonly, twodim_ops).Weight(wgt).Tag("zgboost_bak");
    pm.Push<Hist2D>(Axis(45, 15, 60, "photon_pt[0]", "Photon p_{T} [GeV]", {}),
        Axis(20, -1.0, 1.0, coscaptheta, "cos#Theta", {}),
        baseline, procs_bakonly, twodim_ops).Weight(wgt).Tag("zgboost_bak");
    pm.Push<Hist2D>(Axis(45, 15, 60, "photon_pt[0]", "Photon p_{T} [GeV]", {}),
        Axis(20, -3.15, 3.15, llgphi, "#phi", {}),
        baseline, procs_bakonly, twodim_ops).Weight(wgt).Tag("zgboost_bak");
    pm.Push<Hist2D>(Axis(26, 40, 300, "llphoton_pt[0]", "Higgs Candidate p_{T} [GeV]", {}),
        Axis(20, 100, 160, "llphoton_m[0]", "Higgs Candidate m [GeV]", {}),
        baseline, procs_bakonly, twodim_ops).Weight(wgt).Tag("zgboost_bak");
  }

  if (plot_decorrelate) {
    pm.Push<Hist1D>(Axis(20, 0, 60, "photon_pt[0]", "Photon p_{T} [GeV]", {}), 
        baseline_lep[2]&&higgs_window, procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");
    pm.Push<Hist1D>(Axis(20, -20, 60, decorr_photon_pt, "Decorrelated Photon p_{T} [GeV]", {}), 
        baseline_lep[2]&&higgs_window, procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");
    pm.Push<Hist1D>(Axis(20, 0, 1.0, photon_pt_mass, "Photon p_{T}/Higgs candidate m", {}), 
        baseline_lep[2]&&higgs_window, procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");
    pm.Push<Hist2D>(Axis(20, 0, 60, "photon_pt[0]", "Photon p_{T} [GeV]", {}),
        Axis(20, 100, 160, "llphoton_m[0]", "Higgs Candidate m [GeV]", {}),
        baseline, procs_sigonly, twodim_ops).Weight(wgt).Tag("FixName:zgboost_decorr__sig_ptmcorr");
    pm.Push<Hist2D>(Axis(20, 0, 60, "photon_pt[0]", "Photon p_{T} [GeV]", {}),
        Axis(20, 100, 160, "llphoton_m[0]", "Higgs Candidate m [GeV]", {}),
        baseline, procs_bakonly, twodim_ops).Weight(wgt).Tag("FixName:zgboost_decorr__bak_ptmcorr");
    pm.Push<Hist2D>(Axis(20, -20, 60, decorr_photon_pt, "Decorrelated Photon p_{T} [GeV]", {}),
        Axis(20, 100, 160, "llphoton_m[0]", "Higgs Candidate m [GeV]", {}),
        baseline, procs_sigonly, twodim_ops).Weight(wgt).Tag("FixName:zgboost_decorr__sig_ptmcorr_decorr");
    pm.Push<Hist2D>(Axis(20, -20, 60, decorr_photon_pt, "Decorrelated Photon p_{T} [GeV]", {}),
        Axis(20, 100, 160, "llphoton_m[0]", "Higgs Candidate m [GeV]", {}),
        baseline, procs_bakonly, twodim_ops).Weight(wgt).Tag("FixName:zgboost_decorr__bak_ptmcorr_decorr");
    pm.Push<Hist2D>(Axis(20, 0, 1.0, photon_pt_mass, "Photon p_{T}/Higgs candidate m", {}),
        Axis(20, 100, 160, "llphoton_m[0]", "Higgs Candidate m [GeV]", {}),
        baseline, procs_sigonly, twodim_ops).Weight(wgt).Tag("FixName:zgboost_decorr__sig_ptmcorr_decorrratio");
    pm.Push<Hist2D>(Axis(20, 0, 1.0, photon_pt_mass, "Photon p_{T}/Higgs candidate m", {}),
        Axis(20, 100, 160, "llphoton_m[0]", "Higgs Candidate m [GeV]", {}),
        baseline, procs_bakonly, twodim_ops).Weight(wgt).Tag("FixName:zgboost_decorr__bak_ptmcorr_decorrratio");
    //pm.Push<Hist1D>(Axis(60,-0.5,0.5, MVAScore, "ll#gamma BDT Score", {}), 
    pm.Push<Hist1D>(Axis(60,-0.3,0.3, MVAScore, "ll#gamma BDT Score", {}), 
        baseline_lep[2]&&higgs_window&&"photon_idmva[0]>0.5", procs, ops_sorb).Weight(wgt).Tag("FixName:zgboost_decorr__nodecorrbdt");
    pm.Push<Hist1D>(Axis(60,-0.3,0.3, TopFourBdtScore, "ll#gamma BDT Score (top four only)", {}), 
        baseline_lep[2]&&higgs_window&&"photon_idmva[0]>0.5", procs, ops_sorb).Weight(wgt).Tag("FixName:zgboost_decorr__topfourbdt");
    pm.Push<Hist1D>(Axis(60,-0.3,0.3, DecorrBdtScore, "ll#gamma BDT Score with Decorrelated Photon pt", {}), 
        baseline_lep[2]&&higgs_window&&"photon_idmva[0]>0.5", procs, ops_sorb).Weight(wgt).Tag("FixName:zgboost_decorr__lindecorrbdt");
    pm.Push<Hist1D>(Axis(60,-0.3,0.3, DecorrRatioBdtScore, "ll#gamma BDT Score with Photon pt Ratio", {}), 
        baseline_lep[2]&&higgs_window&&"photon_idmva[0]>0.5", procs, ops_sorb).Weight(wgt).Tag("FixName:zgboost_decorr__ratdecorrbdt");
    //
    //
    pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        baseline_lep[2]&&MVAScore<0.0&&"photon_idmva[0]>0.5", procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");
    pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        baseline_lep[2]&&MVAScore>0.0&&MVAScore<0.05&&"photon_idmva[0]>0.5", procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");
    pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        baseline_lep[2]&&MVAScore>0.05&&MVAScore<0.1&&"photon_idmva[0]>0.5", procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");
    pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        baseline_lep[2]&&MVAScore>0.1&&"photon_idmva[0]>0.5", procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");

    pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        baseline_lep[2]&&DecorrBdtScore<0.033&&"photon_idmva[0]>0.5", procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");
    pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        baseline_lep[2]&&DecorrBdtScore>0.033&&DecorrBdtScore<0.117&&"photon_idmva[0]>0.5", procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");
    pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        baseline_lep[2]&&DecorrBdtScore>0.117&&DecorrBdtScore<0.2&&"photon_idmva[0]>0.5", procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");
    pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        baseline_lep[2]&&DecorrBdtScore>0.2&&"photon_idmva[0]>0.5", procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");

    pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        baseline_lep[2]&&DecorrRatioBdtScore<0.0&&"photon_idmva[0]>0.5", procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");
    pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        baseline_lep[2]&&DecorrRatioBdtScore>0.0&&DecorrRatioBdtScore<0.066&&"photon_idmva[0]>0.5", procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");
    pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        baseline_lep[2]&&DecorrRatioBdtScore>0.066&&DecorrRatioBdtScore<0.133&&"photon_idmva[0]>0.5", procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");
    pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        baseline_lep[2]&&DecorrRatioBdtScore>0.133&&"photon_idmva[0]>0.5", procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");

    pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        baseline_lep[2]&&DecorrRatioBdtScore<0.0&&"photon_idmva[0]>0.5", procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");
    pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        baseline_lep[2]&&DecorrRatioBdtScore>0.0&&DecorrRatioBdtScore<0.066&&"photon_idmva[0]>0.5", procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");
    pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        baseline_lep[2]&&DecorrRatioBdtScore>0.066&&DecorrRatioBdtScore<0.133&&"photon_idmva[0]>0.5", procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");
    pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        baseline_lep[2]&&DecorrRatioBdtScore>0.133&&"photon_idmva[0]>0.5", procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");

    pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        baseline_lep[2]&&DecorrRatioNoMassCutBdtScore<0.0&&"photon_idmva[0]>0.5", procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");
    pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        baseline_lep[2]&&DecorrRatioNoMassCutBdtScore>0.0&&DecorrRatioBdtScore<0.066&&"photon_idmva[0]>0.5", procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");
    pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        baseline_lep[2]&&DecorrRatioNoMassCutBdtScore>0.066&&DecorrRatioBdtScore<0.133&&"photon_idmva[0]>0.5", procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");
    pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        baseline_lep[2]&&DecorrRatioNoMassCutBdtScore>0.133&&"photon_idmva[0]>0.5", procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");
  }

  if (plot_final_bdt_mllg) {
    pm.multithreaded_ = false;
    std::vector<float> run2_bdt_bins = {-1.1,0.0,0.1,0.2,1.1};
    std::vector<float> phpt_bdt_bins = {-1.1,0.04,0.14,0.24,1.1};
    std::vector<float> splt_bdt_bins = {-1.1,0.0,0.08,1.1};
    for (unsigned imva = 0; imva < 4; imva++) {
      //pm.Push<Hist1D>(Axis(30,100,160, "llphoton_m[0]", "m_{ll#gamma}", {}), 
      //    baseline&&MVAScore>run2_bdt_bins[imva]&&MVAScore<=run2_bdt_bins[imva+1], procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_phpt");
      //pm.Push<Hist1D>(Axis(30,100,160, "llphoton_m[0]", "m_{ll#gamma}", {}), 
      //    baseline&&DecorrRatioBdtScore>phpt_bdt_bins[imva]&&DecorrRatioBdtScore<=phpt_bdt_bins[imva+1], procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_phpt");
    }
    for (unsigned imva = 0; imva < 3; imva++) {
      //pm.Push<Hist1D>(Axis(30,100,160, "llphoton_m[0]", "m_{ll#gamma}", {}), 
      //    baseline&&"photon_pt[0]<35"&&LowPhPtBdtScore>splt_bdt_bins[imva]&&LowPhPtBdtScore<=splt_bdt_bins[imva+1], procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_phpt");
      pm.Push<Hist1D>(Axis(30,100,160, "llphoton_m[0]", "m_{ll#gamma}", {}), 
          baseline&&"photon_pt[0]>35"&&HighPhPtBdtScore>splt_bdt_bins[imva]&&HighPhPtBdtScore<=splt_bdt_bins[imva+1], procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_phpt");
    }
  }

  if (plot_photonpt_story) {
    pm.multithreaded_ = false;
    pm.Push<Hist1D>(Axis(20,15,60, "photon_pt[0]", "Photon p_{T} [GeV]", {}), 
        baseline_noph&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops_shapes).Weight(wgt).Tag("zgboost_phpt");
    pm.Push<Hist1D>(Axis(26,0.1,0.5, photon_pt_mass, "Photon p_{T}/m_{ll#gamma}", {}), 
        baseline_noph&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops_shapes).Weight(wgt).Tag("zgboost_phpt");
    pm.Push<Hist1D>(Axis(20,15,60, "photon_pt[0]", "Photon p_{T} [GeV]", {}), 
        baseline_noph&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops_sorb).Weight(wgt).Tag("zgboost_phpt_sorb");
    pm.Push<Hist1D>(Axis(20,0.1,0.5, photon_pt_mass, "Photon p_{T}/m_{ll#gamma}", {}), 
        baseline_noph&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops_sorb).Weight(wgt).Tag("zgboost_phpt_sorb");

    pm.Push<Hist1D>(Axis(60,-0.3,0.3, MVAScore, "ll#gamma BDT Score", {}), 
        baseline_lep[2]&&higgs_window&&"photon_idmva[0]>0.5", procs, ops_sorb).Weight(wgt).Tag("FixName:zgboost_decorr__nodecorrbdt");
    pm.Push<Hist1D>(Axis(60,-0.3,0.3, NodecorrBdtScore, "ll#gamma BDT Score with Photon pt", {}), 
        baseline_lep[2]&&higgs_window&&"photon_idmva[0]>0.5", procs, ops_sorb).Weight(wgt).Tag("FixName:zgboost_decorr__decorrlessbdt");
    pm.Push<Hist1D>(Axis(60,-0.3,0.3, DecorrRatioNoMassCutBdtScore, "ll#gamma BDT Score with Photon pt Ratio (wide window)", {}), 
        baseline_lep[2]&&higgs_window&&"photon_idmva[0]>0.5", procs, ops_sorb).Weight(wgt).Tag("FixName:zgboost_decorr__widedecorrbdt");
    pm.Push<Hist1D>(Axis(60,-0.3,0.3, DecorrRatioBdtScore, "ll#gamma BDT Score with Photon pt Ratio", {}), 
        baseline_lep[2]&&higgs_window&&"photon_idmva[0]>0.5", procs, ops_sorb).Weight(wgt).Tag("FixName:zgboost_decorr__ratdecorrbdt");
    pm.Push<Hist1D>(Axis(60,-0.3,0.3, LowPhPtBdtScore, "ll#gamma BDT Score low photon pT", {}), 
        baseline_lep[2]&&higgs_window&&"photon_idmva[0]>0.5&&photon_pt[0]<35", procs, ops_sorb).Weight(wgt).Tag("FixName:zgboost_decorr__lowphptbdt");
    pm.Push<Hist1D>(Axis(60,-0.3,0.3, HighPhPtBdtScore, "ll#gamma BDT Score low photon pT", {}), 
        baseline_lep[2]&&higgs_window&&"photon_idmva[0]>0.5&&photon_pt[0]>35", procs, ops_sorb).Weight(wgt).Tag("FixName:zgboost_decorr__highphptbdt");

    pm.Push<Hist1D>(Axis(60,-0.3,0.3, MVAScore, "ll#gamma BDT Score", {}), 
        baseline_lep[2]&&higgs_window&&"photon_idmva[0]>0.5", procs, ops_shapes).Weight(wgt).Tag("zgboost_phpt");
    pm.Push<Hist1D>(Axis(60,-0.3,0.3, NodecorrBdtScore, "ll#gamma BDT Score with Photon pt", {}), 
        baseline_lep[2]&&higgs_window&&"photon_idmva[0]>0.5", procs, ops_shapes).Weight(wgt).Tag("zgboost_phpt");
    pm.Push<Hist1D>(Axis(60,-0.3,0.3, DecorrRatioNoMassCutBdtScore, "ll#gamma BDT Score with Photon pt Ratio (wide window)", {}), 
        baseline_lep[2]&&higgs_window&&"photon_idmva[0]>0.5", procs, ops_shapes).Weight(wgt).Tag("zgboost_phpt");
    pm.Push<Hist1D>(Axis(60,-0.3,0.3, DecorrRatioBdtScore, "ll#gamma BDT Score with Photon pt Ratio", {}), 
        baseline_lep[2]&&higgs_window&&"photon_idmva[0]>0.5", procs, ops_shapes).Weight(wgt).Tag("zgboost_phpt");
    pm.Push<Hist1D>(Axis(60,-0.3,0.3, LowPhPtBdtScore, "ll#gamma BDT Score low photon pT", {}), 
        baseline_lep[2]&&higgs_window&&"photon_idmva[0]>0.5&&photon_pt[0]<35", procs, ops_shapes).Weight(wgt).Tag("zgboost_phpt");
    pm.Push<Hist1D>(Axis(60,-0.3,0.3, HighPhPtBdtScore, "ll#gamma BDT Score low photon pT", {}), 
        baseline_lep[2]&&higgs_window&&"photon_idmva[0]>0.5&&photon_pt[0]>35", procs, ops_shapes).Weight(wgt).Tag("zgboost_phpt");

    std::vector<float> mva_bins = {-1.0,-0.01,0.04,0.10,1.0};
    std::vector<float> nondecorr_mva_bins = {-1.0,0.01,0.07,0.14,1.0};
    std::vector<float> large_ratio_mva_bins = {-1.0,0.02,0.09,0.14,1.0};
    std::vector<float> ratio_mva_bins = {-1.0,0.03,0.11,0.16,1.0};

    std::vector<TableRow> hzgyields_rows;
    hzgyields_rows.push_back(TableRow("Inclusive baseline+Higgs window", 
        (baseline)&&"photon_idmva[0]>0.5"&&higgs_window,0,0,wgt));
    //default MVA bins
    for (unsigned imva = 0; imva < mva_bins.size()-1; imva++) {
      hzgyields_rows.push_back(TableRow("MVA bin "+std::to_string(imva),
          (baseline)&&"photon_idmva[0]>0.5"&&higgs_window&&MVAScore>mva_bins[imva]&&MVAScore<=mva_bins[imva+1],0,0,wgt));
      pm.Push<Hist1D>(Axis(30,100,160, "llphoton_m[0]", "m_{ll#gamma}", {}), 
          baseline&&"photon_idmva[0]>0.5"&&MVAScore>mva_bins[imva]&&MVAScore<=mva_bins[imva+1], procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_phpt");
    }
    //MVA bins with photon pt cut
    for (unsigned imva = 0; imva < mva_bins.size()-1; imva++) {
      hzgyields_rows.push_back(TableRow("MVA bin "+std::to_string(imva),
          (baseline)&&"photon_idmva[0]>0.5"&&higgs_window&&MVAScore>mva_bins[imva]&&MVAScore<=mva_bins[imva+1]&&photon_pt_mass>0.18,0,0,wgt));
      pm.Push<Hist1D>(Axis(30,100,160, "llphoton_m[0]", "m_{ll#gamma}", {}), 
          baseline&&"photon_idmva[0]>0.5"&&MVAScore>mva_bins[imva]&&MVAScore<=mva_bins[imva+1]&&photon_pt_mass>0.18, procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_phpt");
    }
    //MVA bins with photon pt
    for (unsigned imva = 0; imva < nondecorr_mva_bins.size()-1; imva++) {
      hzgyields_rows.push_back(TableRow("MVA (Photon pt) bin "+std::to_string(imva),
          (baseline)&&"photon_idmva[0]>0.5"&&higgs_window&&NodecorrBdtScore>nondecorr_mva_bins[imva]&&NodecorrBdtScore<=nondecorr_mva_bins[imva+1],0,0,wgt));
      pm.Push<Hist1D>(Axis(30,100,160, "llphoton_m[0]", "m_{ll#gamma}", {}), 
          baseline&&"photon_idmva[0]>0.5"&&NodecorrBdtScore>nondecorr_mva_bins[imva]&&NodecorrBdtScore<=nondecorr_mva_bins[imva+1], procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_phpt");
    }
    //MVA bins with photon pt ratio (wide)
    for (unsigned imva = 0; imva < mva_bins.size()-1; imva++) {
      hzgyields_rows.push_back(TableRow("MVA (Photon pt ratio large window) bin "+std::to_string(imva),
          (baseline)&&"photon_idmva[0]>0.5"&&higgs_window&&DecorrRatioNoMassCutBdtScore>large_ratio_mva_bins[imva]&&DecorrRatioNoMassCutBdtScore<=large_ratio_mva_bins[imva+1],0,0,wgt));
      pm.Push<Hist1D>(Axis(30,100,160, "llphoton_m[0]", "m_{ll#gamma}", {}), 
          baseline&&"photon_idmva[0]>0.5"&&DecorrRatioNoMassCutBdtScore>large_ratio_mva_bins[imva]&&DecorrRatioNoMassCutBdtScore<=large_ratio_mva_bins[imva+1], procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_phpt");
    }
    //MVA bins with photon pt ratio
    for (unsigned imva = 0; imva < mva_bins.size()-1; imva++) {
      hzgyields_rows.push_back(TableRow("MVA (Photon pt ratio) bin "+std::to_string(imva),
          (baseline)&&"photon_idmva[0]>0.5"&&higgs_window&&DecorrRatioBdtScore>ratio_mva_bins[imva]&&DecorrRatioBdtScore<=ratio_mva_bins[imva+1],0,0,wgt));
      pm.Push<Hist1D>(Axis(30,100,160, "llphoton_m[0]", "m_{ll#gamma}", {}), 
          baseline&&"photon_idmva[0]>0.5"&&DecorrRatioBdtScore>ratio_mva_bins[imva]&&DecorrRatioBdtScore<=ratio_mva_bins[imva+1], procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_phpt");
    }
    //pt-binned MVA bins
    hzgyields_rows.push_back(TableRow("Low pt low MVA bin ",
        baseline&&"photon_idmva[0]>0.5"&&higgs_window&&"photon_pt[0]<35"&&LowPhPtBdtScore<0.02,0,0,wgt));
    hzgyields_rows.push_back(TableRow("Low pt high MVA bin ",
        baseline&&"photon_idmva[0]>0.5"&&higgs_window&&"photon_pt[0]<35"&&LowPhPtBdtScore>0.02,0,0,wgt));
    hzgyields_rows.push_back(TableRow("High pt low MVA bin ",
        baseline&&"photon_idmva[0]>0.5"&&higgs_window&&"photon_pt[0]>35"&&HighPhPtBdtScore<0.04,0,0,wgt));
    hzgyields_rows.push_back(TableRow("High pt high MVA bin ",
        baseline&&"photon_idmva[0]>0.5"&&higgs_window&&"photon_pt[0]>35"&&HighPhPtBdtScore>0.04,0,0,wgt));
    pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        baseline&&"photon_pt[0]<35"&&LowPhPtBdtScore<0.02&&"photon_idmva[0]>0.5", procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");
    pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        baseline&&"photon_pt[0]<35"&&LowPhPtBdtScore>0.02&&"photon_idmva[0]>0.5", procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");
    pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        baseline&&"photon_pt[0]>35"&&HighPhPtBdtScore<0.04&&"photon_idmva[0]>0.5", procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");
    pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        baseline&&"photon_pt[0]>35"&&HighPhPtBdtScore>0.04&&"photon_idmva[0]>0.5", procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");

    pm.Push<Table>("zgboost_phpt_yields",hzgyields_rows,procs,false,true,false,false,false,true).Precision(3);
  }

  if (plot_mva_debug) {
    pm.multithreaded_ = false;
    pm.Push<Hist1D>(Axis(100,-1.0,1.0, TopFourBdtScore, "ll#gamma BDT Score top 4", {}), 
        baseline_debug_training, procs_debug, ops_sorb).Weight("w_lumi").Tag("FixName:zgboost_debug__topfourbdt");
    pm.Push<Hist1D>(Axis(100,-1.0,1.0, TopFiveBdtScore, "ll#gamma BDT Score top 5", {}), 
        baseline_debug_training, procs_debug, ops_sorb).Weight("w_lumi").Tag("FixName:zgboost_debug__topfivebdt");
    pm.Push<Hist1D>(Axis(100,-1.0,1.0, TopSevenBdtScore, "ll#gamma BDT Score top 7", {}), 
        baseline_debug_training, procs_debug, ops_sorb).Weight("w_lumi").Tag("FixName:zgboost_debug__topsevenbdt");
    pm.Push<Hist1D>(Axis(100,-1.0,1.0, TopNineBdtScore, "ll#gamma BDT Score top 9", {}), 
        baseline_debug_training, procs_debug, ops_sorb).Weight("w_lumi").Tag("FixName:zgboost_debug__topninebdt");
    pm.Push<Hist1D>(Axis(100,-1.0,1.0, TopElevenBdtScore, "ll#gamma BDT Score top 11", {}), 
        baseline_debug_training, procs_debug, ops_sorb).Weight("w_lumi").Tag("FixName:zgboost_debug__topelevenbdt");
    
    //pm.Push<Hist1D>(Axis(100,-1.0,1.0, TopFourBdtScore, "ll#gamma BDT Score top 4", {}), 
    //    baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops_sorb).Weight(wgt).Tag("FixName:zgboost_debug__topfourbdt");
    //pm.Push<Hist1D>(Axis(100,-1.0,1.0, TopFiveBdtScore, "ll#gamma BDT Score top 5", {}), 
    //    baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops_sorb).Weight(wgt).Tag("FixName:zgboost_debug__topfivebdt");
    //pm.Push<Hist1D>(Axis(100,-1.0,1.0, TopSevenBdtScore, "ll#gamma BDT Score top 7", {}), 
    //    baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops_sorb).Weight(wgt).Tag("FixName:zgboost_debug__topsevenbdt");
    //pm.Push<Hist1D>(Axis(100,-1.0,1.0, TopNineBdtScore, "ll#gamma BDT Score top 9", {}), 
    //    baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops_sorb).Weight(wgt).Tag("FixName:zgboost_debug__topninebdt");
    //pm.Push<Hist1D>(Axis(100,-1.0,1.0, TopElevenBdtScore, "ll#gamma BDT Score top 11", {}), 
    //    baseline&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops_sorb).Weight(wgt).Tag("FixName:zgboost_debug__topelevenbdt");
  }

  if (plot_ggf_general) {
    pm.multithreaded_ = false;
    //baseline vars
    pm.Push<Hist1D>(Axis(26,15,80, "photon_pt[0]", "Photon p_{T} [GeV]", {}), 
        baseline_noph&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops_shapes).Weight(wgt).Tag("zgboost_phpt");
    pm.Push<Hist1D>(Axis(26,0.0,1.0, photon_pt_mass, "Photon p_{T}/m_{ll#gamma} [GeV]", {}), 
        baseline_noph&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops_shapes).Weight(wgt).Tag("zgboost_phpt");
    pm.Push<Hist1D>(Axis(26,15,80, "photon_pt[0]", "Photon p_{T} [GeV]", {}), 
        baseline_noph&&"photon_idmva[0]>0.5"&&higgs_window, procs_allprocs, ops_sorb).Weight(wgt).Tag("zgboost_ggf_sorb");
    pm.Push<Hist1D>(Axis(26,0.0,1.0, photon_pt_mass, "Photon p_{T}/m_{ll#gamma} [GeV]", {}), 
        baseline_noph&&"photon_idmva[0]>0.5"&&higgs_window, procs_allprocs, ops_sorb).Weight(wgt).Tag("zgboost_ggf_sorb");
    pm.Push<Hist1D>(Axis(20,50,130, "ll_m[0]", "m_{ll} [GeV]", {}), 
        baseline_nomll&&"photon_idmva[0]>0.5"&&higgs_window, procs, ops_shapes).Weight(wgt).Tag("zgboost_phpt");
    pm.Push<Hist1D>(Axis(20,50,130, "ll_m[0]", "m_{ll} [GeV]", {}), 
        baseline_nomll&&"photon_idmva[0]>0.5"&&higgs_window, procs_allprocs, ops).Weight(wgt).Tag("zgboost_phpt");
    pm.Push<Hist2D>(Axis(20, 60, 120, "ll_m[0]", "m_{ll} [GeV]", {}),
        Axis(20, 80, 160, "llphoton_m[0]", "m_{ll#gamma} Candidate m [GeV]", {}),
        baseline_nomll&&"photon_idmva[0]>0.5", procs_sigonly, twodim_ops).Weight(wgt).Tag("zgboost_phpt_sig");
    pm.Push<Hist2D>(Axis(20, 60, 120, "ll_m[0]", "m_{ll} [GeV]", {}),
        Axis(20, 80, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
        baseline_nomll&&"photon_idmva[0]>0.5", procs_allbaks, twodim_ops).Weight(wgt).Tag("zgboost_phpt_bak");
    //BDT vars
    pm.Push<Hist1D>(Axis(20,-1.0,1.0, "photon_idmva[0]", "Photon IDMVA", {}), 
        baseline&&higgs_window, procs_allprocs, ops_sorb).Weight(wgt).Tag("zgboost_phpt_sorb");
    pm.Push<Hist1D>(Axis(20,-1.0,1.0, "photon_idmva[0]", "Photon IDMVA", {}), 
        baseline&&higgs_window, procs, ops_shapes).Weight(wgt).Tag("zgboost_phpt");
    pm.Push<Hist1D>(Axis(20,0.0,1.0, abs_cos_theta, "|cos(#Theta)|", {}), 
        baseline&&higgs_window&&"photon_idmva[0]>0.5", procs_allprocs, ops_sorb_upper).Weight(wgt).Tag("zgboost_phpt_sorb");
    pm.Push<Hist1D>(Axis(20,0.0,1.0, abs_cos_theta, "|cos(#Theta)|", {}), 
        baseline&&higgs_window&&"photon_idmva[0]>0.5", procs, ops_shapes).Weight(wgt).Tag("zgboost_phpt");
    pm.Push<Hist1D>(Axis(20,0.4,3.0, "photon_drmin[0]", "Min #Delta R(l,#gamma)", {}), 
        baseline&&higgs_window&&"photon_idmva[0]>0.5", procs_allprocs, ops_sorb_upper).Weight(wgt).Tag("zgboost_phpt_sorb");
    pm.Push<Hist1D>(Axis(20,0.4,3.0, "photon_drmin[0]", "Min #Delta R(l,#gamma)", {}), 
        baseline&&higgs_window&&"photon_idmva[0]>0.5", procs, ops_shapes).Weight(wgt).Tag("zgboost_phpt");
    pm.Push<Hist1D>(Axis(20,0.0,2.0, pt_mass, "p_{Tll#gamma}/m_{ll#gamma}", {}), 
        baseline&&higgs_window, procs_allprocs, ops_sorb).Weight(wgt).Tag("zgboost_phpt_sorb");
    pm.Push<Hist1D>(Axis(20,0.0,2.0, pt_mass, "p_{Tll#gamma}/m_{ll#gamma}", {}), 
        baseline&&higgs_window, procs, ops_shapes).Weight(wgt).Tag("zgboost_phpt");
  }

  if (plot_photonpt_compare) {

    pm.Push<Hist1D>(Axis(60,-0.5,0.5, DecorrRatioBdtScore, "ll#gamma BDT Score with Photon pt Ratio", {}), 
        baseline_lep[2]&&higgs_window&&"photon_idmva[0]>0.5", procs, ops_sorb).Weight(wgt).Tag("FixName:zgboost_decorr__ratdecorrbdt");
    pm.Push<Hist1D>(Axis(40,0.0,100.0, "photon_pt[0]", "Photon p_{T} [GeV]", {}), 
        baseline_lep[2]&&higgs_window&&"photon_idmva[0]>0.5", procs, ops_sorb).Weight(wgt).Tag("FixName:zgboost_decorr__photonpt");
    pm.Push<Hist1D>(Axis(60,-0.5,0.5, LowPhPtBdtScore, "BDT Score in low photon pt region", {}), 
        baseline_lep[2]&&higgs_window&&"photon_idmva[0]>0.5&&photon_pt[0]<35", procs, ops_sorb).Weight(wgt).Tag("FixName:zgboost_decorr__lowphptbdt");
    pm.Push<Hist1D>(Axis(60,-0.5,0.5, HighPhPtBdtScore, "BDT Score in high photon pt region", {}), 
        baseline_lep[2]&&higgs_window&&"photon_idmva[0]>0.5&&photon_pt[0]>35", procs, ops_sorb).Weight(wgt).Tag("FixName:zgboost_decorr__highphptbdt");

    pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        baseline_lep[2]&&"photon_pt[0]<35"&&LowPhPtBdtScore<0.0&&"photon_idmva[0]>0.5", procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");
    pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        baseline_lep[2]&&"photon_pt[0]<35"&&LowPhPtBdtScore>0.0&&"photon_idmva[0]>0.5", procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");
    pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        baseline_lep[2]&&"photon_pt[0]>35"&&HighPhPtBdtScore<0.05&&"photon_idmva[0]>0.5", procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");
    pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        baseline_lep[2]&&"photon_pt[0]>35"&&HighPhPtBdtScore>0.05&&"photon_idmva[0]>0.5", procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_decorr");

    std::vector<float> mva_bins = {-1.0,0.0,0.066,0.133,1.0};
    std::vector<TableRow> hzgyields_rows;
    hzgyields_rows.push_back(TableRow("Inclusive baseline+Higgs window", 
        baseline_lep[2]&&"photon_idmva[0]>0.5"&&higgs_window,0,0,wgt));
    for (unsigned imva = 0; imva < mva_bins.size()-1; imva++) {
      hzgyields_rows.push_back(TableRow("MVA bin "+std::to_string(imva),
          baseline_lep[2]&&"photon_idmva[0]>0.5"&&higgs_window&&DecorrRatioBdtScore>mva_bins[imva]&&DecorrRatioBdtScore<=mva_bins[imva+1],0,0,wgt));
    }
    hzgyields_rows.push_back(TableRow("Low pt low MVA bin ",
        baseline_lep[2]&&"photon_idmva[0]>0.5"&&higgs_window&&"photon_pt[0]<35"&&LowPhPtBdtScore<0.0,0,0,wgt));
    hzgyields_rows.push_back(TableRow("Low pt high MVA bin ",
        baseline_lep[2]&&"photon_idmva[0]>0.5"&&higgs_window&&"photon_pt[0]<35"&&LowPhPtBdtScore>0.0,0,0,wgt));
    hzgyields_rows.push_back(TableRow("High pt low MVA bin ",
        baseline_lep[2]&&"photon_idmva[0]>0.5"&&higgs_window&&"photon_pt[0]>35"&&HighPhPtBdtScore<0.05,0,0,wgt));
    hzgyields_rows.push_back(TableRow("High pt high MVA bin ",
        baseline_lep[2]&&"photon_idmva[0]>0.5"&&higgs_window&&"photon_pt[0]>35"&&HighPhPtBdtScore>0.05,0,0,wgt));
    pm.Push<Table>("hzg_phpt_yields",hzgyields_rows,procs,false,true,false,false,false,true).Precision(3);
  }

  //if (plot_fitshapes) {
  //  pm.Push<Hist1D>(Axis(40, 80, 160, truth_mllg, "Truth m_{ll#gamma} [GeV]", {}), 
  //      stitch_zg, procs_dy_noskim, ops).Weight(wgt_sigx20).Tag("zgboost_fitshapes");
  //  pm.Push<Hist1D>(Axis(40, 0, 100, truth_photon_pt, "Truth photon p_{T} [GeV]", {}), 
  //      stitch_zg, procs_dy_noskim, ops).Weight(wgt_sigx20).Tag("zgboost_fitshapes");
  //  pm.Push<Hist1D>(Axis(40, 0, 100, truth_photon_pt, "Truth photon p_{T} [GeV]", {}), 
  //      stitch_zg, procs_dy_noskim, ops_log).Weight(wgt_sigx20).Tag("zgboost_fitshapes_log");
  //  //std::vector<float> phpt_bins = {0,5,10,15,20,25,35,999};
  //  //for (unsigned iphpt = 0; iphpt < phpt_bins.size()-1; iphpt++) {
  //  //  pm.Push<Hist1D>(Axis(40, 80, 160, truth_mllg, "Truth m_{ll#gamma} [GeV]", {}), 
  //  //      stitch_zg&&truth_photon_pt>phpt_bins[iphpt]&&truth_photon_pt<phpt_bins[iphpt+1], procs_dy_noskim, ops).Weight(wgt_sigx20).Tag("zgboost_fitshapes");
  //  //  pm.Push<Hist1D>(Axis(40, 80, 160, truth_mllg, "Truth m_{ll#gamma} [GeV]", {}), 
  //  //      stitch_zg&&truth_photon_pt>phpt_bins[iphpt], procs_dy_noskim, ops).Weight(wgt_sigx20).Tag("zgboost_fitshapes");
  //  //}
  //  pm.Push<Hist1D>(Axis(40, 80, 160, truth_mllg, "Truth m_{ll#gamma} [GeV]", {}), 
  //      stitch_zg&&truth_photon_pt>15, procs_dy_noskim, ops).Weight(wgt_sigx20).Tag("zgboost_fitshapes");
  //  pm.Push<Hist1D>(Axis(40, 80, 160, truth_mllg, "Truth m_{ll#gamma} [GeV]", {}), 
  //      stitch_zg&&truth_photon_pt>15&&truth_photon_abseta<2.5, procs_dy_noskim, ops).Weight(wgt_sigx20).Tag("zgboost_fitshapes");
  //  pm.Push<Hist1D>(Axis(40, 80, 160, truth_mllg, "Truth m_{ll#gamma} [GeV]", {}), 
  //      stitch_zg&&truth_photon_pt>15&&truth_photon_abseta<2.5&&truth_lepton_minpt>15, procs_dy_noskim, ops).Weight(wgt_sigx20).Tag("zgboost_fitshapes");
  //  pm.Push<Hist1D>(Axis(40, 80, 160, truth_mllg, "Truth m_{ll#gamma} [GeV]", {}), 
  //      stitch_zg&&truth_photon_pt>15&&truth_photon_abseta<2.5&&truth_lepton_minpt>15&&truth_lepton_eta_acc, procs_dy_noskim, ops).Weight(wgt_sigx20).Tag("zgboost_fitshapes");
  //  pm.Push<Hist1D>(Axis(40, 80, 160, truth_mllg, "Truth m_{ll#gamma} [GeV]", {}), 
  //      stitch_zg&&truth_photon_pt>15&&truth_photon_abseta<2.5&&truth_lepton_minpt>15&&truth_lepton_eta_acc&&truth_min_dr>0.4, procs_dy_noskim, ops).Weight(wgt_sigx20).Tag("zgboost_fitshapes");

  //  //pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
  //  //    baseline_lep[2], procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_fitshapes");
  //  //pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
  //  //    baseline_lep[2]&&"photon_idmva[0]>0.5", procs_sigx20, ops).Weight(wgt_sigx20).Tag("zgboost_fitshapes");
  //  //pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
  //  //    baseline_lep[2], procs_sigx20, ops_shapes).Weight(wgt_sigx20).Tag("zgboost_fitshapes_shapes");
  //  //pm.Push<Hist1D>(Axis(30, 100, 160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
  //  //    baseline_lep[2]&&"photon_idmva[0]>0.5", procs_sigx20, ops_shapes).Weight(wgt_sigx20).Tag("zgboost_fitshapes_shapes");
  //}

  if (plot_yieldtable) {
    NamedFunc photon_pt0 = "photon_pt[0]";
    NamedFunc llphoton_dr0 = "llphoton_dr[0]";

    std::vector<float> mva_bins = {-1.0,0.0,0.55,0.85,1.0};
    std::vector<float> phpt_bins = {15,25,35,61,999};
    std::vector<float> hgdr_bins = {0,1.3,2.3,3.0,6.0};
    std::vector<float> boostbdt_bins = {-1.0,0.0,0.1,0.3,1.0};

    std::vector<TableRow> hzgyields_rows;
    hzgyields_rows.push_back(TableRow("Inclusive baseline+Higgs window", 
        (baseline)&&higgs_window,0,0,wgt));
    hzgyields_rows.push_back(TableRow("Photon pt>49 GeV", 
        (baseline)&&higgs_window&&photon_pt0>49,0,0,wgt));
    hzgyields_rows.push_back(TableRow("Boost BDT score>0.3", 
        (baseline)&&higgs_window&&BoostScore>0.3,0,0,wgt));
    for (unsigned imva = 0; imva < mva_bins.size()-1; imva++) {
      hzgyields_rows.push_back(TableRow("MVA bin "+std::to_string(imva),
          (baseline)&&higgs_window&&MVAScore>mva_bins[imva]&&MVAScore<=mva_bins[imva+1],0,0,wgt));
    }
    for (unsigned imva = 0; imva < mva_bins.size()-1; imva++) {
      hzgyields_rows.push_back(TableRow("MVA bin "+std::to_string(imva)+", Veto photon pt bin 49",
          (baseline)&&higgs_window&&MVAScore>mva_bins[imva]&&MVAScore<=mva_bins[imva+1]&&photon_pt0<49,0,0,wgt));
    }
    for (unsigned imva = 0; imva < mva_bins.size()-1; imva++) {
      hzgyields_rows.push_back(TableRow("MVA bin "+std::to_string(imva)+", Veto Boost BDT bin 0.3",
          (baseline)&&higgs_window&&MVAScore>mva_bins[imva]&&MVAScore<=mva_bins[imva+1]&&BoostScore<0.3,0,0,wgt));
    }
    for (unsigned imva = 0; imva < mva_bins.size()-1; imva++) {
      for (unsigned iphpt = 0; iphpt < phpt_bins.size()-1; iphpt++) {
        hzgyields_rows.push_back(TableRow("MVA bin "+std::to_string(imva)+", Photon pt bin "+std::to_string(iphpt),
            (baseline)&&higgs_window&&MVAScore>mva_bins[imva]&&MVAScore<=mva_bins[imva+1]&&photon_pt0>phpt_bins[iphpt]&&photon_pt0<=phpt_bins[iphpt+1],0,0,wgt));
      }
    }
    for (unsigned imva = 0; imva < mva_bins.size()-1; imva++) {
      for (unsigned ihgdr = 0; ihgdr < hgdr_bins.size()-1; ihgdr++) {
        hzgyields_rows.push_back(TableRow("MVA bin "+std::to_string(imva)+", Higgs #Delta R bin "+std::to_string(ihgdr),
            (baseline)&&higgs_window&&MVAScore>mva_bins[imva]&&MVAScore<=mva_bins[imva+1]&&llphoton_dr0>hgdr_bins[ihgdr]&&llphoton_dr0<=hgdr_bins[ihgdr+1],0,0,wgt));
      }
    }
    for (unsigned imva = 0; imva < mva_bins.size()-1; imva++) {
      for (unsigned iboostbdt = 0; iboostbdt < boostbdt_bins.size()-1; iboostbdt++) {
        hzgyields_rows.push_back(TableRow("MVA bin "+std::to_string(imva)+", Boost BDT bin "+std::to_string(iboostbdt),
            (baseline)&&higgs_window&&MVAScore>mva_bins[imva]&&MVAScore<=mva_bins[imva+1]&&BoostScore>boostbdt_bins[iboostbdt]&&BoostScore<=boostbdt_bins[iboostbdt+1],0,0,wgt));
      }
    }
    for (unsigned iboostbdt = 0; iboostbdt < boostbdt_bins.size()-1; iboostbdt++) {
      hzgyields_rows.push_back(TableRow("Boost BDT bin "+std::to_string(iboostbdt),
          (baseline)&&higgs_window&&BoostScore>boostbdt_bins[iboostbdt]&&BoostScore<=boostbdt_bins[iboostbdt+1],0,0,wgt));
    }

    pm.Push<Table>("hzgyields",hzgyields_rows,procs,false,true,false,false,false,true).Precision(3);
      //std::vector<TableRow>{

      ////MVA binning
      //TableRow("Baseline+Higgs window, MVA bin 1", 
      //    (baseline)&&higgs_window&&MVAScore<-0.159,0,0,wgt),
      //TableRow("Baseline+Higgs window, MVA bin 2", 
      //    (baseline)&&higgs_window&&MVAScore>-0.159&&MVAScore<0.285,0,0,wgt),
      //TableRow("Baseline+Higgs window, MVA bin 3", 
      //    (baseline)&&higgs_window&&MVAScore>0.285&&MVAScore<0.683,0,0,wgt),
      //TableRow("Baseline+Higgs window, MVA bin 4", 
      //    (baseline)&&higgs_window&&MVAScore>0.683,0,0,wgt),

      ////MVA binning+Photon p_{T} binning
      //TableRow("Baseline+Higgs window, MVA bin 1", 
      //    (baseline)&&higgs_window&&MVAScore<-0.159,0,0,wgt),
      //TableRow("Baseline+Higgs window, MVA bin 2", 
      //    (baseline)&&higgs_window&&MVAScore>-0.159&&MVAScore<0.285,0,0,wgt),

      //TableRow("Baseline+Higgs window, MVA bin 3, $p_{T\\gamma}>61 GeV$", 
      //    (baseline)&&higgs_window&&MVAScore>0.285&&MVAScore<0.683&&"photon_pt[0]>61",0,0,wgt),

      //TableRow("Baseline+Higgs window, MVA bin 4, $p_{T\\gamma}<25$ GeV", 
      //    (baseline)&&higgs_window&&MVAScore>0.683&&"photon_pt[0]<25",0,0,wgt),
      //TableRow("Baseline+Higgs window, MVA bin 4, $25<p_{T\\gamma}<35$ GeV", 
      //    (baseline)&&higgs_window&&MVAScore>0.683&&"photon_pt[0]>25&&photon_pt[0]<35",0,0,wgt),
      //TableRow("Baseline+Higgs window, MVA bin 4, $35<p_{T\\gamma}<61$ GeV", 
      //    (baseline)&&higgs_window&&MVAScore>0.683&&"photon_pt[0]>35&&photon_pt[0]<61",0,0,wgt),
      //TableRow("Baseline+Higgs window, MVA bin 4, $p_{T\\gamma}>61$ GeV", 
      //    (baseline)&&higgs_window&&MVAScore>0.683&&"photon_pt[0]>61",0,0,wgt),

      ////MVA binning+Delta R binning
      //TableRow("Baseline+Higgs window, MVA bin 1", 
      //    (baseline)&&higgs_window&&MVAScore<-0.159,0,0,wgt),
      //TableRow("Baseline+Higgs window, MVA bin 2", 
      //    (baseline)&&higgs_window&&MVAScore>-0.159&&MVAScore<0.285,0,0,wgt),
      //TableRow("Baseline+Higgs window, MVA bin 3", 
      //    (baseline)&&higgs_window&&MVAScore>0.285&&MVAScore<0.683,0,0,wgt),
      //TableRow("Baseline+Higgs window, MVA bin 4", 
      //    (baseline)&&higgs_window&&MVAScore>0.683,0,0,wgt),

      ////TableRow("Baseline+Higgs window, $p_{T\\gamma}<23$ GeV", 
      ////    (baseline)&&higgs_window&&"photon_pt[0]<23",0,0,wgt),
      ////TableRow("Baseline+Higgs window, $23<p_{T\\gamma}<28$ GeV", 
      ////    (baseline)&&higgs_window&&"photon_pt[0]>23&&photon_pt[0]<28",0,0,wgt),
      ////TableRow("Baseline+Higgs window, $28<p_{T\\gamma}<35$ GeV", 
      ////    (baseline)&&higgs_window&&"photon_pt[0]>28&&photon_pt[0]<35",0,0,wgt),
      ////TableRow("Baseline+Higgs window, $35<p_{T\\gamma}<45$ GeV", 
      ////    (baseline)&&higgs_window&&"photon_pt[0]>35&&photon_pt[0]<45",0,0,wgt),
      ////TableRow("Baseline+Higgs window, $45<p_{T\\gamma}<60$ GeV", 
      ////    (baseline)&&higgs_window&&"photon_pt[0]>45&&photon_pt[0]<60",0,0,wgt),
      ////TableRow("Baseline+Higgs window, $p_{T\\gamma}>60$ GeV", 
      ////    (baseline)&&higgs_window&&"photon_pt[0]>60",0,0,wgt),

      ////TableRow("Baseline+Higgs window, $p_{Tll}<30$ GeV", 
      ////    (baseline)&&higgs_window&&"ll_pt[0]<30",0,0,wgt),
      ////TableRow("Baseline+Higgs window, $30<p_{Tll}<50$ GeV", 
      ////    (baseline)&&higgs_window&&"ll_pt[0]>30&&ll_pt[0]<50",0,0,wgt),
      ////TableRow("Baseline+Higgs window, $50<p_{Tll}<75$ GeV", 
      ////    (baseline)&&higgs_window&&"ll_pt[0]>50&&ll_pt[0]<75",0,0,wgt),
      ////TableRow("Baseline+Higgs window, $75<p_{Tll}<115$ GeV", 
      ////    (baseline)&&higgs_window&&"ll_pt[0]>75&&ll_pt[0]<115",0,0,wgt),
      ////TableRow("Baseline+Higgs window, $115<p_{Tll}<130$ GeV", 
      ////    (baseline)&&higgs_window&&"ll_pt[0]>115&&ll_pt[0]<130",0,0,wgt),
      ////TableRow("Baseline+Higgs window, $p_{Tll}>130$ GeV", 
      ////    (baseline)&&higgs_window&&"ll_pt[0]>130",0,0,wgt),

      ////TableRow("Baseline+Higgs window, $p_{Tll\\gamma}<25$ GeV", 
      ////    (baseline)&&higgs_window&&"llphoton_pt[0]<25",0,0,wgt),
      ////TableRow("Baseline+Higgs window, $25<p_{Tll\\gamma}<45$ GeV", 
      ////    (baseline)&&higgs_window&&"llphoton_pt[0]>25&&llphoton_pt[0]<45",0,0,wgt),
      ////TableRow("Baseline+Higgs window, $45<p_{Tll\\gamma}<70$ GeV", 
      ////    (baseline)&&higgs_window&&"llphoton_pt[0]>45&&llphoton_pt[0]<70",0,0,wgt),
      ////TableRow("Baseline+Higgs window, $70<p_{Tll\\gamma}<130$ GeV", 
      ////    (baseline)&&higgs_window&&"llphoton_pt[0]>70&&llphoton_pt[0]<130",0,0,wgt),
      ////TableRow("Baseline+Higgs window, $130<p_{Tll\\gamma}<160$ GeV", 
      ////    (baseline)&&higgs_window&&"llphoton_pt[0]>130&&llphoton_pt[0]<160",0,0,wgt),
      ////TableRow("Baseline+Higgs window, $p_{Tll\\gamma}>160$ GeV", 
      ////    (baseline)&&higgs_window&&"llphoton_pt[0]>160",0,0,wgt),

      ////TableRow("Baseline+Higgs window, $\\Delta R<3.0$ GeV", 
      ////    (baseline)&&higgs_window&&"llphoton_dr[0]<3.0",0,0,wgt),
      ////TableRow("Baseline+Higgs window, $2.5<\\Delta R<3.0$ GeV", 
      ////    (baseline)&&higgs_window&&"llphoton_dr[0]>2.5&&llphoton_dr[0]<3.0",0,0,wgt),
      ////TableRow("Baseline+Higgs window, $1.7<\\Delta R<2.5$ GeV", 
      ////    (baseline)&&higgs_window&&"llphoton_dr[0]>1.7&&llphoton_dr[0]<2.5",0,0,wgt),
      ////TableRow("Baseline+Higgs window, $0.8<\\Delta R<1.7$ GeV", 
      ////    (baseline)&&higgs_window&&"llphoton_dr[0]>0.8&&llphoton_dr[0]<1.7",0,0,wgt),
      ////TableRow("Baseline+Higgs window, $0.4<\\Delta R<0.8$ GeV", 
      ////    (baseline)&&higgs_window&&"llphoton_dr[0]>0.4&&llphoton_dr[0]<0.8",0,0,wgt),
      ////TableRow("Baseline+Higgs window, $\\Delta R<0.4$ GeV", 
      ////    (baseline)&&higgs_window&&"llphoton_dr[0]<0.4",0,0,wgt),

      ////TableRow("Baseline+Higgs window, MVA bin 4, $p_{T\\gamma}<20$ GeV", 
      ////    (baseline)&&higgs_window&&MVAScore>0.683&&"photon_pt[0]<20",0,0,wgt),
      ////TableRow("Baseline+Higgs window, MVA bin 4, $20<p_{T\\gamma}<25$ GeV", 
      ////    (baseline)&&higgs_window&&MVAScore>0.683&&"photon_pt[0]>20&&photon_pt[0]<25",0,0,wgt),
      ////TableRow("Baseline+Higgs window, MVA bin 4, $25<p_{T\\gamma}<30$ GeV", 
      ////    (baseline)&&higgs_window&&MVAScore>0.683&&"photon_pt[0]>25&&photon_pt[0]<30",0,0,wgt),
      ////TableRow("Baseline+Higgs window, MVA bin 4, $30<p_{T\\gamma}<40$ GeV", 
      ////    (baseline)&&higgs_window&&MVAScore>0.683&&"photon_pt[0]>30&&photon_pt[0]<40",0,0,wgt),
      ////TableRow("Baseline+Higgs window, MVA bin 4, $40<p_{T\\gamma}<50$ GeV", 
      ////    (baseline)&&higgs_window&&MVAScore>0.683&&"photon_pt[0]>40&&photon_pt[0]<50",0,0,wgt),
      ////TableRow("Baseline+Higgs window, MVA bin 4, $p_{T\\gamma}>50$ GeV", 
      ////    (baseline)&&higgs_window&&MVAScore>0.683&&"photon_pt[0]>50",0,0,wgt),

      ////TableRow("Baseline+Higgs window, MVA bin 3, $p_{T\\gamma}<20$ GeV", 
      ////    (baseline)&&higgs_window&&MVAScore>0.285&&MVAScore<0.683&&"photon_pt[0]<20",0,0,wgt),
      ////TableRow("Baseline+Higgs window, MVA bin 3, $20<p_{T\\gamma}<25$ GeV", 
      ////    (baseline)&&higgs_window&&MVAScore>0.285&&MVAScore<0.683&&"photon_pt[0]>20&&photon_pt[0]<25",0,0,wgt),
      ////TableRow("Baseline+Higgs window, MVA bin 3, $25<p_{T\\gamma}<30$ GeV", 
      ////    (baseline)&&higgs_window&&MVAScore>0.285&&MVAScore<0.683&&"photon_pt[0]>25&&photon_pt[0]<30",0,0,wgt),
      ////TableRow("Baseline+Higgs window, MVA bin 3, $30<p_{T\\gamma}<40$ GeV", 
      ////    (baseline)&&higgs_window&&MVAScore>0.285&&MVAScore<0.683&&"photon_pt[0]>30&&photon_pt[0]<40",0,0,wgt),
      ////TableRow("Baseline+Higgs window, MVA bin 3, $40<p_{T\\gamma}<50$ GeV", 
      ////    (baseline)&&higgs_window&&MVAScore>0.285&&MVAScore<0.683&&"photon_pt[0]>40&&photon_pt[0]<50",0,0,wgt),
      ////TableRow("Baseline+Higgs window, MVA bin 3, $p_{T\\gamma}>50$ GeV", 
      ////    (baseline)&&higgs_window&&MVAScore>0.285&&MVAScore<0.683&&"photon_pt[0]>50",0,0,wgt),

      ////TableRow("Baseline+Higgs window, MVA bin 2, $p_{T\\gamma}<20$ GeV", 
      ////    (baseline)&&higgs_window&&MVAScore>-0.159&&MVAScore<0.285&&"photon_pt[0]<20",0,0,wgt),
      ////TableRow("Baseline+Higgs window, MVA bin 2, $20<p_{T\\gamma}<25$ GeV", 
      ////    (baseline)&&higgs_window&&MVAScore>-0.159&&MVAScore<0.285&&"photon_pt[0]>20&&photon_pt[0]<25",0,0,wgt),
      ////TableRow("Baseline+Higgs window, MVA bin 2, $25<p_{T\\gamma}<30$ GeV", 
      ////    (baseline)&&higgs_window&&MVAScore>-0.159&&MVAScore<0.285&&"photon_pt[0]>25&&photon_pt[0]<30",0,0,wgt),
      ////TableRow("Baseline+Higgs window, MVA bin 2, $30<p_{T\\gamma}<40$ GeV", 
      ////    (baseline)&&higgs_window&&MVAScore>-0.159&&MVAScore<0.285&&"photon_pt[0]>30&&photon_pt[0]<40",0,0,wgt),
      ////TableRow("Baseline+Higgs window, MVA bin 2, $40<p_{T\\gamma}<50$ GeV", 
      ////    (baseline)&&higgs_window&&MVAScore>-0.159&&MVAScore<0.285&&"photon_pt[0]>40&&photon_pt[0]<50",0,0,wgt),
      ////TableRow("Baseline+Higgs window, MVA bin 2, $p_{T\\gamma}>50$ GeV", 
      ////    (baseline)&&higgs_window&&MVAScore>-0.159&&MVAScore<0.285&&"photon_pt[0]>50",0,0,wgt),
    //},procs,false,true,false,false,false,true).Precision(3);
  }

  //pm.multithreaded_ = false;
  //if (plot_fitshapes) {
  //  //enable for plotting that doesn't use BDTs
  //  pm.multithreaded_ = true;
  //}
  pm.min_print_ = true;
  pm.MakePlots(340.0);

  if (plot_decorrelate) {
    Hist2D * ptmcorr = static_cast<Hist2D*>(pm.GetFigure("zgboost_decorr__bak_ptmcorr").get());
    TH2D plot_ptmcorr = ptmcorr->GetBkgHist(true);
    std::cout << "Mean of var 1 is " << plot_ptmcorr.GetMean(1) << std::endl;
    std::cout << "Mean of var 2 is " << plot_ptmcorr.GetMean(2) << std::endl;
    std::cout << "Variance of var 1 is " << plot_ptmcorr.GetStdDev(1) << std::endl;
    std::cout << "Variance of var 2 is " << plot_ptmcorr.GetStdDev(2) << std::endl;
    std::cout << "Covariance is " << plot_ptmcorr.GetCovariance() << std::endl;
    Hist2D * ptmcorr2 = static_cast<Hist2D*>(pm.GetFigure("zgboost_decorr__bak_ptmcorr_decorr").get());
    TH2D plot_ptmcorr2 = ptmcorr2->GetBkgHist(true);
    std::cout << "Mean of var 1 is " << plot_ptmcorr2.GetMean(1) << std::endl;
    std::cout << "Mean of var 2 is " << plot_ptmcorr2.GetMean(2) << std::endl;
    std::cout << "Variance of var 1 is " << plot_ptmcorr2.GetStdDev(1) << std::endl;
    std::cout << "Variance of var 2 is " << plot_ptmcorr2.GetStdDev(2) << std::endl;
    std::cout << "Covariance is " << plot_ptmcorr2.GetCovariance() << std::endl;
  }

  delete reader;
  delete boost_tagger_hgpt;
  delete boost_tagger_hgdr;
  delete boost_tagger_zpt;
  delete boost_tagger_phpt;
  delete decorr_photon_pt_ptr;
  delete photon_pt_mass_ptr;

  return 0;
}
