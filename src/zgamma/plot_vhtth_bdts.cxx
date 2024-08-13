#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <bitset>
#include <unistd.h>
#include <getopt.h>
#include "TError.h"
#include "TColor.h"
#include "TVector3.h"
#include "TMVA/Reader.h"
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
#include "zgamma/zg_functions.hpp"
#include "zgamma/categorization_utilities.hpp"
using namespace std;
using namespace PlotOptTypes;
using namespace ZgUtilities;
using namespace ZgFunctions;
using namespace CatUtilities;

//This file is used to plot the different BDTs created for the ttH and VH categories

double wgt_tree_select(const Baby &b, int tree_select){ 
  if(b.SampleTypeString().Contains("-")){return 1;}
  return b.event()%314159%3==tree_select ? wgt.GetScalar(b)*3 : 0; 
}


NamedFunc wgt_train("wgt_train",[](const Baby &b) -> NamedFunc::ScalarType{ return wgt_tree_select(b,0); });
NamedFunc wgt_valid("wgt_valid",[](const Baby &b) -> NamedFunc::ScalarType{ return wgt_tree_select(b,1); });
NamedFunc wgt_test( "wgt_test", [](const Baby &b) -> NamedFunc::ScalarType{ return wgt_tree_select(b,2); });



//All global variables that help, Find best place to delete reader when done
Float_t pdmax_ft  = 0;
Float_t yrV_ft    = 0;
Float_t yres_ft   = 0;
Float_t yrap_ft   = 0;
Float_t yidmva_ft = 0;
Float_t pdmin_ft  = 0;
Float_t l1eta_ft  = 0;
Float_t l2eta_ft  = 0;
Float_t ll_pt_ft  = 0;
Float_t ll_m_ft   = 0;

//WH 3l and ZH MET bdt variables
Float_t l3_pt_ft   = 0;
Float_t l3_mini_ft = 0;
Float_t met_ft     = 0;
Float_t dphi_met_lly_ft = 0;

//ttH bdt variables
Float_t njet_ft = 0;
Float_t mtop_ft = 0;
Float_t ht_ft   = 0;
Float_t max_deepflavor_ft = 0;
Float_t max_imini_ft = 0;
Float_t sumpt4jom_ft = 0;

TMVA::Reader *bdt_reader = nullptr;
int bdt_select = 0;
int event_cache = -1;
double BDT_score_cache = -1;
string bdt_name = "";

//Function that returns BDT score
Float_t mva_score(const Baby &b){
  //Attempt at a cache
  if(event_cache == b.event()){return BDT_score_cache;}
  
  pdmax_ft  = static_cast<Float_t>(b.photon_drmax() -> at(0));
  pdmin_ft  = static_cast<Float_t>(b.photon_drmin() -> at(0));
  yidmva_ft = static_cast<Float_t>(b.photon_idmva() -> at(0));
  yrap_ft   = static_cast<Float_t>(b.photon_eta() -> at(0));
  yres_ft   = static_cast<Float_t>((b.photon_energyErr()->at(0))/(b.photon_pt() ->at(0)));
  ll_pt_ft  = static_cast<Float_t>(b.ll_pt() -> at(0));
  ll_m_ft   = static_cast<Float_t>(b.ll_m() -> at(0));

  if(bdt_select==3){
    // 'l3_pt', 'l3_imini', 'met', 'lly_met_dphi'
    l3_pt_ft   = static_cast<Float_t>( l3_pt_d(b) );
    l3_mini_ft = static_cast<Float_t>( l3_mini.GetScalar(b) );
    met_ft = static_cast<Float_t>(b.met());
    dphi_met_lly_ft = static_cast<Float_t>(dphi_Hmet_d(b));

  } else if(bdt_select==2){
    //'met', 'lly_met_dphi'
    dphi_met_lly_ft = static_cast<Float_t>(dphi_Hmet_d(b));
    met_ft = static_cast<Float_t>(b.met());

  } else if(bdt_select==1){
    //'mtop', 'njet', 'ht', 'sumjptom_4j', 'max_deepflavor'
    njet_ft = static_cast<Float_t>(b.njet());
    mtop_ft = static_cast<Float_t>(mtop.GetScalar(b));
    ht_ft = static_cast<Float_t>(b.ht());
    max_deepflavor_ft = static_cast<Float_t>(max_deepflav_d(b));
    sumpt4jom_ft = static_cast<Float_t>(ptj_m4j.GetScalar(b));

  } else if(bdt_select==0){
    //'max_imini', 'max_deepflavor', 'met', 'ht'
    max_imini_ft = static_cast<Float_t>(max_Imini.GetScalar(b));
    max_deepflavor_ft = static_cast<Float_t>(max_deepflav_d(b));
    met_ft = static_cast<Float_t>(b.met());
    ht_ft = static_cast<Float_t>(b.ht());
  }
 
  //EvaluateMVA to get BDT score <-- Without Caching this is not needed. Can directly return value of BDT
  event_cache = b.event(); 
  BDT_score_cache = (bdt_reader -> TMVA::Reader::EvaluateMVA(bdt_name));
  return BDT_score_cache;
}


//Wrapper around mva_score function so it can be returned easily in draw_pico
NamedFunc SignalBDT("SignalBDT",[](const Baby &b) -> NamedFunc::ScalarType{ 
    return mva_score(b);
});


//map<string,vector<TString>> variables_bdt= {}
/*
ifcategory=='WH_3l':
features=['y_mva','yl_drmin','yl_drmax','lly_ptom','y_res','y_eta','l3_pt','l3_imini','leplep_pt','leplep_m','met','lly_met_dphi']
elifcategory=='ZH_MET':
features=['y_eta',]
elifcategory=='ttH_had':
features=['y_mva','yl_drmin','yl_drmax','lly_ptom',]
elifcategory=='ttH_lep':
features=['y_mva','yl_drmin','yl_drmax','lly_ptom','y_res','y_eta',]
else:


*/

//Wrapper around mva_score function so it can be returned easily in draw_pico
NamedFunc mlly_datacut("mlly_datacut",[](const Baby &b) -> NamedFunc::ScalarType{ 
  if( !(b.SampleTypeString().Contains("-")) ){ return true;}
  if( b.llphoton_m() -> at(0) < 122 || b.llphoton_m() -> at(0) > 128){ return true; }
  return false;
});


int main(int argc, char *argv[]) {
  string bdt_folder = "/net/cms27/cms27r0/abarzdukas/bdttools/htozgamma_mva/vhtthbdts/";

  //int category = -1;
  string bdt = "";
  if(argc<2){
    cout << "Please give the option of a category." << endl;
    return -1;
  } else {
    //category = std::atoi(argv[1]);
    bdt_select = std::atoi(argv[1]);

  }

/*
  if(argv<3){
    cout << "Please give the bdt selection option" << endl;
    return -1;
  } else {
    bdt_select = argv[2];
  }
*/
  cout << "debug 0" << endl;

  Float_t *pdmax_ptr = &pdmax_ft;
  Float_t *yrV_ptr = &yrV_ft;
  Float_t *yres_ptr = &yres_ft;
  Float_t *yrap_ptr = &yrap_ft;
  Float_t *yidmva_ptr = &yidmva_ft;
  Float_t *pdmin_ptr = &pdmin_ft;
  //Float_t *l1eta_ptr = &l1eta_ft;
  //Float_t *l2eta_ptr = &l2eta_ft;
  Float_t *ll_pt_ptr = &ll_pt_ft;
  Float_t *ll_m_ptr = &ll_m_ft;

  //WH 3l and ZH MET bdt variables
  Float_t *l3_pt_ptr = &l3_pt_ft;
  Float_t *l3_mini_ptr = &l3_mini_ft;
  Float_t *met_ptr = &met_ft;
  Float_t *dphi_met_lly_ptr = &dphi_met_lly_ft;

  //ttH bdt variables
  Float_t *njet_ptr = &njet_ft;
  Float_t *mtop_ptr = &mtop_ft;
  Float_t *ht_ptr = &ht_ft;
  Float_t *max_deepflavor_ptr = &max_deepflavor_ft;
  Float_t *max_imini_ptr = &max_imini_ft;
  Float_t *sumjpt4jom_ptr = &sumpt4jom_ft;
 
  cout << "debug 1" << endl;
  bdt_reader = new TMVA::Reader;

  float* lly_m_sp = 0; float* w1 = 0; float* w2 = 0; float* w3 = 0;
  bdt_reader -> TMVA::Reader::AddSpectator("lly_m", lly_m_sp );
  bdt_reader -> TMVA::Reader::AddSpectator("weightXyear", w1 );
  bdt_reader -> TMVA::Reader::AddSpectator("w_lumiXyear", w2 );
  bdt_reader -> TMVA::Reader::AddSpectator("sampleID", w3 );

  bdt_reader -> TMVA::Reader::AddVariable("y_mva", yidmva_ptr);
  bdt_reader -> TMVA::Reader::AddVariable("yl_drmin", pdmin_ptr);
  bdt_reader -> TMVA::Reader::AddVariable("yl_drmax", pdmax_ptr);
  bdt_reader -> TMVA::Reader::AddVariable("lly_ptom", yrV_ptr);
  bdt_reader -> TMVA::Reader::AddVariable("y_res", yres_ptr);

  cout << "debug 1" << endl;

  double bdt_selections = -1;
  if(bdt_select==3){
    bdt_name = "tmva_bdta_vhtth_WH_3l_mt_1_win_1_model";
    bdt_selections = 0.0;

    bdt_reader -> TMVA::Reader::AddVariable("y_eta", yrap_ptr);
    bdt_reader -> TMVA::Reader::AddVariable("l3_pt", l3_pt_ptr);
    bdt_reader -> TMVA::Reader::AddVariable("l3_imini", l3_mini_ptr);
    bdt_reader -> TMVA::Reader::AddVariable("leplep_pt", ll_pt_ptr);
    bdt_reader -> TMVA::Reader::AddVariable("leplep_m", ll_m_ptr);
    bdt_reader -> TMVA::Reader::AddVariable("met", met_ptr);
    bdt_reader -> TMVA::Reader::AddVariable("lly_met_dphi", dphi_met_lly_ptr);
    bdt_reader -> TMVA::Reader::BookMVA(bdt_name, bdt_folder + bdt_name + "/weights/TMVAClassification_BDT.weights.xml"); 

  } else if(bdt_select==2){
    bdt_name = "tmva_bdta_vhtth_ZH_MET_mt_4_win_0_model";
    bdt_selections = 0.03;

    bdt_reader -> TMVA::Reader::AddVariable("y_eta", yrap_ptr);
    bdt_reader -> TMVA::Reader::AddVariable("met", met_ptr);
    bdt_reader -> TMVA::Reader::AddVariable("lly_met_dphi", dphi_met_lly_ptr);
    bdt_reader -> TMVA::Reader::AddVariable("leplep_pt", ll_pt_ptr);
    bdt_reader -> TMVA::Reader::AddVariable("leplep_m", ll_m_ptr);
    bdt_reader -> TMVA::Reader::BookMVA(bdt_name, bdt_folder + bdt_name + "/weights/TMVAClassification_BDT.weights.xml"); 

  } else if(bdt_select==1){
    bdt_name = "tmva_bdta_vhtth_ttH_had_mt_2_win_1_model";
    bdt_selections = 0.0;

    bdt_reader -> TMVA::Reader::AddVariable("mtop", mtop_ptr);
    bdt_reader -> TMVA::Reader::AddVariable("njet", njet_ptr);
    bdt_reader -> TMVA::Reader::AddVariable("ht", ht_ptr);
    bdt_reader -> TMVA::Reader::AddVariable("sumjptom_4j", sumjpt4jom_ptr);
    bdt_reader -> TMVA::Reader::AddVariable("max_deepflavor", max_deepflavor_ptr);
    bdt_reader -> TMVA::Reader::AddVariable("leplep_pt",  ll_pt_ptr);
    bdt_reader -> TMVA::Reader::AddVariable("leplep_m", ll_m_ptr);
    bdt_reader -> TMVA::Reader::BookMVA(bdt_name, bdt_folder + bdt_name + "/weights/TMVAClassification_BDT.weights.xml"); 

  } else if(bdt_select==0){
    bdt_name = "tmva_bdta_vhtth_ttH_lep_mt_1_win_1_model";
    bdt_selections = 0.2;
 
    bdt_reader -> TMVA::Reader::AddVariable("y_eta", yrap_ptr);
    bdt_reader -> TMVA::Reader::AddVariable("max_imini", max_imini_ptr);
    bdt_reader -> TMVA::Reader::AddVariable("max_deepflavor", max_deepflavor_ptr);
    bdt_reader -> TMVA::Reader::AddVariable("met", met_ptr);
    bdt_reader -> TMVA::Reader::AddVariable("leplep_pt",  ll_pt_ptr);
    bdt_reader -> TMVA::Reader::AddVariable("leplep_m", ll_m_ptr);
    bdt_reader -> TMVA::Reader::AddVariable("ht", ht_ptr);
    bdt_reader -> TMVA::Reader::BookMVA(bdt_name, bdt_folder + bdt_name + "/weights/TMVAClassification_BDT.weights.xml"); 
  }



  cout << "debug 2: Book" << endl;

  


  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");

  //Defines the plot options for 1D and 2D plots
  PlotOpt style("txt/plot_styles.txt","Eff2D");
  vector<PlotOpt> bkg_hist = {style().Stack(StackType::data_norm)
                                     .YTitleOffset(1.)
                                     .Overflow(OverflowType::both)
                                     .LabelSize(0.04)
                                     .UseCMYK(true)  
                                     .YAxis(YAxisType::log)
                                     .Title(TitleType::info)};
  vector<PlotOpt> ops_2D = {bkg_hist};

  PlotOpt lin_lumi("txt/plot_styles.txt","CMSPaper");
  lin_lumi.YAxis(YAxisType::linear)
          .Stack(StackType::signal_overlay)
          .YTitleOffset(1.75)
          .AutoYAxis(false)
          .UseCMYK(false)
          .LeftMargin(0.17)
          .CanvasWidth(800)
          .CanvasHeight(800)
          .Title(TitleType::info)
          .FileExtensions({"pdf"});
  vector<PlotOpt> ops = {lin_lumi()};


  cout << "debug 2.5" << endl;

  //Uses txt/zg_samples.txt to define needed processes
  //vector<shared_ptr<Process>> procs =  ZgUtilities::procs_with_sig_scale("txt/samples_zgamma.txt","AllMoreTrigs",10);
  //vector<shared_ptr<Process>> procs = ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_zgamma.txt","AllMoreTrigsData");
  vector<shared_ptr<Process>> procs = ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_zgamma.txt","ttXonly");

  //This defines the splits based on lepton flavor that will be plotted
  vector<NamedFunc> lep = { "ll_lepid[0]==11","ll_lepid[0]==13","1"};
  vector<string> lep_lab= {"_ee", "_mumu", "_ll"};

  //This defines the categorization that will be used for plotting
  vector<NamedFunc> cat_vec = CatUtilities::run3_category_vector;
  vector<string> cat_vec_str= CatUtilities::run3_category_labels;                          

  //These are all the selections used as a part of the baseline selection from the zg_functions.cpp file
  NamedFunc r2_baseline    = ZgFunctions::run2_baseline;
  NamedFunc tight_baseline = ZgFunctions::tightened_baseline;

  //This block of code loops through to create the NamedFuncs for each category
  //The purpose is to only require one loop while making plots
  vector<NamedFunc> NamedFunc_loop_tight  = {};
  vector<string>    string_loop_label_tight = {};
  for(unsigned int idx_i = 0; idx_i < lep.size(); idx_i++){
    //Vectors used for plots with the tightened baseline selection
    NamedFunc_loop_tight.push_back( tight_baseline && cat_vec[bdt_select] && lep[idx_i]);
    string_loop_label_tight.push_back(cat_vec_str[bdt_select] + lep_lab[idx_i] + "_TIGHT_BASELINE");
  }

  //This block of code handles all of the plot making
  PlotMaker pm;
  string labels = "";
  NamedFunc selection = "1";
  vector<TableRow> higgs_window_yields = {};

  cout << "debug 3: Pls wurk" << endl;

  //Makes plots with the tightened baseline selection
  //These are not necessary but sometimes make it a bit more concise when making plots
  for(unsigned int idx_lep = 0; idx_lep < lep.size(); idx_lep++){

    labels = string_loop_label_tight[idx_lep];
    selection  = NamedFunc_loop_tight[idx_lep] && mlly_datacut;
    pm.Push<Hist1D>(Axis(50, -1, 1, SignalBDT, "BDT score", {}), selection, procs, ops).Weight(wgt_test).Tag("ShortName:cat_vhtth_bdts_" + labels + "_bdt_score");

    higgs_window_yields.push_back(TableRow(labels, selection && "llphoton_m[0] > 122 && llphoton_m[0] < 128", 0,0, wgt_test));

    labels = string_loop_label_tight[idx_lep] + "_lbdtselection";
    selection  = NamedFunc_loop_tight[idx_lep] && mlly_datacut && SignalBDT < bdt_selections;

    pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_refit_m", "m_{ll#gamma}", {}), selection, procs, ops).Weight(wgt_test).Tag("ShortName:cat_vhtth_bdts_" + labels + "_refit_mlly");


    if(bdt_select==0){ CatUtilities::ttH_lep_controlregion_plots(pm,selection,procs,ops,wgt_test,labels); }
    if(bdt_select==1){ CatUtilities::ttH_had_controlregion_plots(pm,selection,procs,ops,wgt_test,labels); }
    if(bdt_select==2){  CatUtilities::ZH_MET_controlregion_plots(pm,selection,procs,ops,wgt_test,labels); }
    if(bdt_select==3){   CatUtilities::WH_3l_controlregion_plots(pm,selection,procs,ops,wgt_test,labels); }
    if(bdt_select==4){     CatUtilities::VBF_controlregion_plots(pm,selection,procs,ops,wgt_test,labels); }
    if(bdt_select==5){     CatUtilities::ggF_controlregion_plots(pm,selection,procs,ops,wgt_test,labels); }

    higgs_window_yields.push_back(TableRow(labels, selection && "llphoton_m[0] > 122 && llphoton_m[0] < 128", 0,0, wgt_test));

    labels = string_loop_label_tight[idx_lep] + "_gbdtselection";
    selection  = NamedFunc_loop_tight[idx_lep] && mlly_datacut && SignalBDT > bdt_selections;

    pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_refit_m", "m_{ll#gamma}", {}), selection, procs, ops).Weight(wgt_test).Tag("ShortName:cat_vhtth_bdts_" + labels + "_refit_mlly");

    if(bdt_select==0){ CatUtilities::ttH_lep_controlregion_plots(pm,selection,procs,ops,wgt_test,labels); }
    if(bdt_select==1){ CatUtilities::ttH_had_controlregion_plots(pm,selection,procs,ops,wgt_test,labels); }
    if(bdt_select==2){  CatUtilities::ZH_MET_controlregion_plots(pm,selection,procs,ops,wgt_test,labels); }
    if(bdt_select==3){   CatUtilities::WH_3l_controlregion_plots(pm,selection,procs,ops,wgt_test,labels); }
    if(bdt_select==4){     CatUtilities::VBF_controlregion_plots(pm,selection,procs,ops,wgt_test,labels); }
    if(bdt_select==5){     CatUtilities::ggF_controlregion_plots(pm,selection,procs,ops,wgt_test,labels); }

    higgs_window_yields.push_back(TableRow(labels, selection && "llphoton_m[0] > 122 && llphoton_m[0] < 128", 0,0, wgt_test));
  }

  //Tables of yields for each final state signal MC
  pm.Push<Table>("vhtth_bdts_" + to_string(bdt_select) , higgs_window_yields, procs, 0, 1, 0, 1, 0, 1).Precision(3);


  //This block controls the verbose or non-verbose plotting option and whether or not to use multiple threads
  pm.min_print_ = true;
  pm.multithreaded_ = true;

  //Luminosity option for the plots
  pm.SetLuminosityTag("137.61").MakePlots(1);

}
