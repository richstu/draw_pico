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
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMVA/Reader.h"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/slide_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/utilities.hpp"
#include "core/mva_wrapper.hpp"
#include "zgamma/zg_utilities.hpp"
#include "zgamma/zg_functions.hpp"
#include "zgamma/categorization_utilities.hpp"
using namespace std;
using namespace PlotOptTypes;
using namespace ZgUtilities;
using namespace ZgFunctions;
using namespace CatUtilities;

//Defining NamedFuncs used for plotting
int l1_eta_d(const Baby &b){
  int l1 = b.ll_i1() -> at(0); int l2 = b.ll_i2() -> at(0);

  if(b.ll_lepid()->at(0)==11){
    return b.el_pt() -> at(l1) > b.el_pt() -> at(l2) ? b.el_eta() -> at(l1) : b.el_eta() -> at(l2);
  }
  return b.mu_pt() -> at(l1) > b.mu_pt() -> at(l2) ? b.mu_eta() -> at(l1) : b.mu_eta() -> at(l2);
}


int l2_eta_d(const Baby &b){
  int l1 = b.ll_i1() -> at(0); int l2 = b.ll_i2() -> at(0);

  if(b.ll_lepid()->at(0)==11){
    return b.el_pt() -> at(l1) < b.el_pt() -> at(l2) ? b.el_eta() -> at(l1) : b.el_eta() -> at(l2);
  }
  return b.mu_pt() -> at(l1) < b.mu_pt() -> at(l2) ? b.mu_eta() -> at(l1) : b.mu_eta() -> at(l2);
}


NamedFunc l1_eta("l1_eta", [](const Baby &b) -> NamedFunc::ScalarType{ return l1_eta_d(b); });

NamedFunc l2_eta("l2_eta", [](const Baby &b) -> NamedFunc::ScalarType{ return l2_eta_d(b); });

NamedFunc yj1_dr("yj1_dr", [](const Baby &b) -> NamedFunc::ScalarType{
  int ij1 = idx_j1(b);
  if(ij1<0){return -999;}
  TLorentzVector ph = AssignGamma(b);
  TLorentzVector j1; j1.SetPtEtaPhiM(b.jet_pt()->at(ij1), b.jet_eta()->at(ij1), b.jet_phi()->at(ij1), b.jet_m()->at(ij1));
  return ph.DeltaR(j1); 
});

NamedFunc yj2_dr("yj2_dr", [](const Baby &b) -> NamedFunc::ScalarType{
  int ij2 = idx_j2(b);
  if(ij2<0){return -999;}
  TLorentzVector ph = AssignGamma(b);
  TLorentzVector j2; j2.SetPtEtaPhiM(b.jet_pt()->at(ij2), b.jet_eta()->at(ij2), b.jet_phi()->at(ij2), b.jet_m()->at(ij2));
  return ph.DeltaR(j2); 
});

NamedFunc photon_energyerr("photon_energyerr", [](const Baby &b) -> NamedFunc::ScalarType{
  return (b.photon_energyErr() ->at(0))/AssignGamma(b).E(); 
});

NamedFunc blind_data("blind_data", [](const Baby &b) -> NamedFunc::ScalarType{
  if(b.SampleTypeString().Contains("-") && b.nllphoton() > 0 && b.llphoton_m()->at(0) > 120 && b.llphoton_m()->at(0) < 130){return false;}
  return true;
});

int main(int argc, char *argv[]) {
  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");
  NamedFunc type("rtype",  [](const Baby &b) -> NamedFunc::ScalarType{ return b.type(); });

  //Uses first input to determine whether or not to use run2 or run3 samples.
  //The second input (if there) will instead plot VBF plots
  string plot_option = "";
  //bool vbf_cat = true;
  if(argc>1){ plot_option=argv[1]; }
  //if(argc>2){ vbf_cat = false;}

  //Declares the samples
  cout << "Rui" << endl;
  //  vector<shared_ptr<Process>> procs = ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_zgamma_pinnacles.txt","BDTMCData");//"MinimalMoreTrigsData");//"AllMoreTrigsData");
  string plot_names = "_run2_1";
  string lumi_label = "137.61";
  vector<shared_ptr<Process>> procs_mc = ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_zgamma_pinnacles.txt","BDTMC");//"AllSplit");
  cout << "Rui" << endl;
  if(plot_option=="run3"){
    //    procs_mc      = ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_zgamma_pinnacles.txt","BDTMCDataRun3");
    procs_mc      = ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_zgamma_pinnacles.txt","BDTMCRun3");
    plot_names = "_run3";
    lumi_label = "62.32";
    cout << "Plotting Run 3" << endl;
  }else{
    cout << "Plotting Run 2" << endl;
  }
 
  //Define and name the BDT that will be used. Change the path of BookMVA to your own BDT and change the vbf_bin_boundaries to the optimized values.
  //One step I will likely take is find a way to condense this down.
  // MVAWrapper vbf_bdt_reader("vbf_01j_bdt");
  // vbf_bdt_reader.SetVariable("y_mva","photon_idmva[0]"); //1
  // vbf_bdt_reader.SetVariable("yl_drmin","photon_drmin[0]"); //2
  // vbf_bdt_reader.SetVariable("yl_drmax","photon_drmax[0]"); //3
  // vbf_bdt_reader.SetVariable("cosTheta","llphoton_cosTheta[0]"); //4
  // vbf_bdt_reader.SetVariable("costheta","llphoton_costheta[0]"); //5
  // vbf_bdt_reader.SetVariable("phi","llphoton_psi[0]"); //6 
  // vbf_bdt_reader.SetVariable("lly_ptmass","llphoton_pt[0]/llphoton_m[0]"); //7
  // vbf_bdt_reader.SetVariable("y_res",photon_energyerr); //8
  // vbf_bdt_reader.SetVariable("y_eta","photon_eta[0]"); //9 
  // vbf_bdt_reader.SetVariable("l1_eta",l1_eta); //10
  // vbf_bdt_reader.SetVariable("l2_eta",l2_eta); //11
  // vbf_bdt_reader.SetSpectator("sampleID","type");
  // vbf_bdt_reader.BookMVA("/net/cms27/cms27r0/abarzdukas/bdttools/htozgamma_mva/vbfvbf/tmva_bdta_vbfvbf_tight_01j_higwindow_var11_model/weights/TMVAClassification_BDT.weights.xml");
  //  NamedFunc vbf_bdt_score = vbf_bdt_reader.GetDiscriminant();
  //  vector<double> vbf_bin_boundaries={0.042, 0.094, 0.29};
  //  vector<NamedFunc> vbf_categories = {vbf_bdt_score < vbf_bin_boundaries[0], vbf_bdt_score > vbf_bin_boundaries[0] && vbf_bdt_score < vbf_bin_boundaries[1], 
  //                                      vbf_bdt_score > vbf_bin_boundaries[1] && vbf_bdt_score < vbf_bin_boundaries[2], vbf_bdt_score > vbf_bin_boundaries[2]}; 


  //VBF BDT declaration. Commented out but feel free to use.
  /*
  MVAWrapper vbf_bdt_reader("vbf_2j_bdt");
  vbf_bdt_reader.SetVariable("y_mva","photon_idmva[0]"); //1
  vbf_bdt_reader.SetVariable("yl_drmin","photon_drmin[0]"); //2
  vbf_bdt_reader.SetVariable("yl_drmax","photon_drmax[0]"); //3
  vbf_bdt_reader.SetVariable("cosTheta","llphoton_cosTheta[0]"); //4
  vbf_bdt_reader.SetVariable("costheta","llphoton_costheta[0]"); //5
  vbf_bdt_reader.SetVariable("phi","llphoton_psi[0]"); //6
  vbf_bdt_reader.SetVariable("y_res",photon_energyerr); //7
  vbf_bdt_reader.SetVariable("y_eta","photon_eta[0]"); //8
  vbf_bdt_reader.SetVariable("l1_eta",l1_eta); //9
  vbf_bdt_reader.SetVariable("l2_eta",l2_eta); //10
  vbf_bdt_reader.SetVariable("lly_ptmass","llphoton_pt[0]/llphoton_m[0]"); //11

  vbf_bdt_reader.SetVariable("lly_ptt","llphoton_pTt[0]"); //12
  vbf_bdt_reader.SetVariable("jj_deta","dijet_deta"); //13
  vbf_bdt_reader.SetVariable("jj_dphi","dijet_dphi"); //14
  vbf_bdt_reader.SetVariable("yj1_dr",yj1_dr); //15
  vbf_bdt_reader.SetVariable("yj2_dr",yj2_dr); //16
  vbf_bdt_reader.SetVariable("llyjj_dphi","llphoton_dijet_dphi[0]"); //17
  vbf_bdt_reader.SetVariable("j1_pt",j1_pt); //18
  vbf_bdt_reader.SetVariable("j2_pt",j2_pt); //19
  vbf_bdt_reader.SetVariable("llyjj_ptbal","llphoton_dijet_balance[0]"); //20
  vbf_bdt_reader.SetVariable("yjj_zep","photon_zeppenfeld[0]"); //21
  vbf_bdt_reader.SetSpectator("sampleID","type");
  vbf_bdt_reader.BookMVA("/net/cms27/cms27r0/abarzdukas/bdttools/htozgamma_mva/vbfvbf/tmva_bdta_vbfvbf_tight_2j_higwindow_var21_model/weights/TMVAClassification_BDT.weights.xml");
  NamedFunc vbf_bdt_score = vbf_bdt_reader.GetDiscriminant();
  vector<double> vbf_bin_boundaries={0.057,0.165,0.317};
  vector<NamedFunc> vbf_categories = {vbf_bdt_score < vbf_bin_boundaries[0], vbf_bdt_score > vbf_bin_boundaries[0] && vbf_bdt_score < vbf_bin_boundaries[1], 
                                      vbf_bdt_score > vbf_bin_boundaries[1] && vbf_bdt_score < vbf_bin_boundaries[2], vbf_bdt_score > vbf_bin_boundaries[2]}; 
  */

  //PlotOpt for control regions plots. This should (and soon will) be defined in controlregion_utilities.cpp when I push that.
  //  PlotOpt lin_lumi("txt/plot_styles.txt","Std1D");
  PlotOpt lin_lumi("txt/plot_styles.txt","CMSPaper");//2022");
  lin_lumi.Title(TitleType::preliminary)
          .YAxis(YAxisType::linear)
          .Overflow(OverflowType::both)
          .YTitleOffset(1.)
          .LogMinimum(0.001)
          .AutoYAxis(true)
          .UseCMYK(false)
          .LeftMargin(0.17)
          .LegendColumns(1)
          .CanvasWidth(1077)
          .CanvasWidth(900)
          .Stack(StackType::data_norm) //<---STACK TYPE TO NORMALIZE MC TO DATA. MAKE SURE TO ALSO not_hwin and NOT blind_data. comment out if you want to plot out of box
          .Bottom(PlotOptTypes::BottomType::ratio) 
          .RatioMinimum(0.5)
          .RatioMaximum(1.5)
          .FileExtensions({"pdf"});


  PlotOpt lin_lumi_shape("txt/plot_styles.txt","CMSPaper2022");
  lin_lumi_shape.Title(TitleType::preliminary)
    .YAxis(YAxisType::linear)
    .Overflow(OverflowType::both)
    .YTitleOffset(1.)
    .LogMinimum(0.001)
    .AutoYAxis(true)
    .UseCMYK(false)
    .LeftMargin(0.1)
    .LegendColumns(1)
    .CanvasWidth(1077)
    .CanvasWidth(900)
    .Stack(StackType::prop_shape_stack)
    .FileExtensions({"pdf"});


  //Ops vector used for plots
  vector<PlotOpt> ops = {lin_lumi};
  vector<PlotOpt> ops_shape = {lin_lumi_shape};

  //Baseline selection used for the plots
  //  NamedFunc baseline = ZgFunctions::tightened_baseline_pinnacles;
  //NamedFunc baseline = ZgFunctions::tightened_baseline; //<--- IF USING LASSEN_V0 OR EARLIER USE THIS


  //  vector<string>    run3_cat_labels = CatUtilities::run3_category_labels;                          
  //  vector<NamedFunc> run3_catwsel    = CatUtilities::run3_catwsel_vector;

  //If you want to plot electrons or muons this vector can be used in the loops below.
  vector<NamedFunc> leps = {"ll_lepid[0] == 11","ll_lepid[0] == 13","1"};
  vector<string> leps_str= {"_ee","_mumu","_ll"};

  //Declaring & filling vectors used to loop over during plotting
  std::vector<NamedFunc> vbf_categories_vector    = {};
  std::vector<string> vbf_categories_name_vector  = {};
  std::vector<string> string_vector               = {"_bin1","_bin2","_bin3","_bin4"};
  //  for(unsigned int idx_i = 0; idx_i < vbf_categories.size(); idx_i++){
  //      vbf_categories_vector.push_back(baseline && cat_vbf && vbf_categories[idx_i] );
  //      vbf_categories_name_vector.push_back("vbf" + string_vector[idx_i] + plot_names );

      //Alternative looping over lepton flavor
      //for(unsigned int idx_lep = 0; idx_lep < leps.size(); idx_lep++){
      //  vbf_categories_vector.push_back(baseline && cat_vbf && leps[idx_lep] && vbf_categories[idx_i] );
      //  vbf_categories_name_vector.push_back("vbf" + string_vector[idx_i] + leps_str[idx_lep] + plot_names );
      //}
  //  }

  //VBF code is also commented out here
  /*
  std::vector<NamedFunc> vbf_categories_vector = {};
  std::vector<string> vbf_categories_name_vector  = {};

  for(unsigned int idx_i = 0; idx_i < vbf_categories.size(); idx_i++){
    for(unsigned int idx_lep = 0; idx_lep < leps.size(); idx_lep++){
      vbf_categories_vector.push_back(tightbase && cat_VBF && leps[idx_lep] && vbf_categories[idx_i] );
      vbf_categories_name_vector.push_back("VBF" + string_vector[idx_i] + leps_str[idx_lep] + plot_names );
    }
  }
  */

  //Declaring PlotMaker
  PlotMaker pm;

  //Initial definitions of variables used in the loop
  NamedFunc selection_blind = "1"; string labels = "VBF"+plot_names;
  NamedFunc selection_hmass = "1";
  NamedFunc not_hwin = mlly < 120 || mlly > 130; //<------ ALTERNATIVE TO BLIND_DATA. NEEDED IF USING StackType("data_norm")
  NamedFunc hmasswindow = mlly > 120 && mlly < 130; 
  
  //  selection = baseline && cat_VBF && hmasswindow; //&& not_hwin;// blind_data;
  selection_blind = cat_VBF && not_hwin;//phmasswindow;
  selection_hmass = cat_VBF && hmasswindow;
  cout << "Rui" << endl;
  //  CatUtilities::VBF_input_plots(pm,selection_blind,procs,ops,wgt,labels);
  //  for(unsigned int idx_sel = 0; idx_sel < vbf_categories_vector.size(); idx_sel++){
    //    selection = vbf_categories_vector[idx_sel] && blind_data;
    //    labels = vbf_categories_name_vector[idx_sel];  
  //CatUtilities::vbf_controlregion_plots(pm,selection,procs,ops,wgt_pin_fix,labels);
  //  CatUtilities::vbf_input_plots(pm,selection_blind,procs,ops,wgt,labels);
  cout << "Rui" << endl;
  CatUtilities::VBF_input_plots(pm,selection_hmass,procs_mc,ops_shape,wgt,labels);  
  cout << "Rui" << endl;
    //CatUtilities::vbf_controlregion_plots(pm,selection,procs,ops,wgt_pin_fix,labels,bins); //<--- Option that allows you to set the bins for ALL plots
    //  }
  //  CatUtilities::vbf_input_plots(pm,selection,procs,ops,wgt,labels);  
  //Plot options and MakePlots command
  //pm.multithreaded_ = false;
  pm.min_print_ = true;
  pm.SetLuminosityTag(lumi_label).MakePlots(1);

}

  //Implementation using the vbf_cat variable set by the input arguments to this file
  /*
  if(vbf_cat){
    for(unsigned int idx_sel = 0; idx_sel < vbf_categories_vector.size(); idx_sel++){
      selection = vbf_categories_vector[idx_sel];
      labels = vbf_categories_name_vector[idx_sel];  
      CatUtilities::vbf_controlregion_plots(pm,selection,procs,ops,wgt,labels);
      //CatUtilities::vbf_controlregion_plots(pm,selection,procs,ops,wgt,labels,bins); //<--- Option that allows you to set the bins for ALL plots
    }
  } else {
    for(unsigned int idx_sel = 0; idx_sel < vbf_categories.size(); idx_sel++){
      selection = vbf_categories_vector[idx_sel];
      labels = vbf_categories_name_vector[idx_sel]; 
      CatUtilities::VBF_controlregion_plots(pm,selection,procs,ops,wgt,labels);
      //CatUtilities::VBF_controlregion_plots(pm,selection,procs,ops,wgt,labels,bins); //<--- Option that allows you to set the bins for ALL plots

    }
  }
  */


