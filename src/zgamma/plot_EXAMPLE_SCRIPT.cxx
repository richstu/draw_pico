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

//temp
#include <Math/Rotation3D.h>
#include <Math/RotationX.h>
#include <Math/RotationY.h>
#include <Math/RotationZ.h>
#include <Math/Vector4D.h>
#include <Math/Boost.h>
#include <Math/BoostX.h>
#include <Math/BoostY.h>
#include <Math/BoostZ.h>

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
#include "zgamma/zg_utilities.hpp"
#include "zgamma/zg_functions.hpp"
using namespace std;
using namespace PlotOptTypes;
using namespace ZgUtilities;
using namespace ZgFunctions;

//------------------------------------GIST--------------------------------------//
//This plot file makes a number of plots using a variety of draw_pico methods
//------------------------------------------------------------------------------//

//NamedFuncs. Use NamedFuncs to perform operations beyond simple cuts on pico branches.
                     
NamedFunc cutbitcheck(int bit){//NamedFunc to check cutflow bitmap
  return NamedFunc("comp",[bit](const Baby &b)->NamedFunc::ScalarType{
    return (bit & (b.zg_cutBitMap()))>0;
  });
}

NamedFunc catbitcheck(int bit){//NamedFunc to check categorization bitmap
  return NamedFunc("comp",[bit](const Baby &b)->NamedFunc::ScalarType{
    return (bit == (b.zg_categorizationBitMap()));
  });
}

NamedFunc hastz("hastz",[](const Baby &b) -> NamedFunc::ScalarType{//Example NamedFunc that identifies generator Zs in picos
  if(b.mc_pt()->size()==0) return false;
  bool foundz = false;
  for(unsigned int imc(0); imc<b.mc_pt()->size(); imc++){
    if(b.mc_id()->at(imc)==23) foundz = true;
  }
  return foundz;
});

NamedFunc wgt_std("wgt_std",[](const Baby &b) -> NamedFunc::ScalarType{//NamedFunc which applies weights to events
  if(b.SampleTypeString().Contains("-")) {
    return 1;
  }
  //if( b.type() >= 200000 && b.type() <= 200500 ){return b.weight()*w_years.GetScalar(b)*100;}
  return b.weight()*w_years.GetScalar(b);
});

//Main plotting function
int main(){
 
  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");

//===========================Define Processes=============================//

  //Manual sample loading:
  //vector<shared_ptr<Process>> internal_name       = Process::MakeShared<Baby_pico>("Plot legend name", Process::Type::[signal,background,data], [TColor of choice], {"pico_path", "pico_path_2",..., "pico_path_n"}, {"simple_cut_1",...,"simple_cut_n"});
  auto ggF_2022     = Process::MakeShared<Baby_pico>("ggF 2022", Process::Type::background, kBlack, {"/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022/mc/skim_ll/*GluGluHtoZG_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8*.root"},"1");
  auto vbf_2022     = Process::MakeShared<Baby_pico>("vbf 2022", Process::Type::background, kBlue, {"/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022/mc/skim_ll/*VBFHtoZG_Zto2L_M-125_TuneCP5_withDipoleRecoil_13p6TeV_powheg-pythia8*.root"},"1");
 
  //Automated sample loading: One can define much more complicated sample loading patterns in the txt files.
  //vector<shared_ptr<Process>> internal_name       = ZgUtilities::ZgSampleLoader().LoadSamples("sample_loading_txt_file","key");
  vector<shared_ptr<Process>> MC_test = ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_zgamma.txt","MCtest");
  //Choose samples to plot
  vector<string> procs_name = {"_ggF_2022_","_vbf_2022","_mc_test_","both_ggf_vbf"};//Names to append to plot files
  vector<vector<shared_ptr<Process>>> procs_test = {{ggF_2022},{vbf_2022},{MC_test},{ggF_2022, vbf_2022}};//combinations of samples to plot
  vector<string> e_label = {"137.61 fb^{-1} (13 TeV) + 62.32 fb^{-1} (13.6 TeV)"};//A custom energy and luminosity string can be used as desired
  vector<string> lumi_vec_str = {"7.98","7.98","7.98","7.98"};//luminosities for display


  //This code creates the plotting styles for the 1D and 2D histograms
  //One can either take plotting styles direct from input txt 
  //or take a style and modify it.
  vector<PlotOpt> style_1D = {PlotOpt("txt/plot_styles.txt","LinLumi").Title(TitleType::preliminary).TitleSize(0.032).DisplayMCNorm(true)};
  vector<PlotOpt> style_2D = {PlotOpt("txt/plot_styles.txt","LogLumi2D").Title(TitleType::preliminary).FileExtensions({"pdf", "root"})};
  //This defines the plotting object used to create plots
  PlotMaker pm;

  //Define additional namedfunctions if necessary
  NamedFunc init = "pass==1 && use_event";
  NamedFunc base_half = cutbitcheck(0b100000000000) && cutbitcheck(0b000100000000) && cutbitcheck(0b000010000000) && cutbitcheck(0b000001000000);//half baseline selection


//===========================Generate Plots=============================//


for(unsigned int i(0); i < procs_test.size();  i++){//loop over combinations of samples to plot
  //2D histogram example
  //pm.Push<Hist2D>(Axis(bins, low edge, high edge, "pico_branch" or NamedFunc, "Axis label", {Paired values to create simple cut window displayed on plot},
  //                Axis(bins, low edge, high edge, "pico_branch" or NamedFunc, "Axis label", {Paired values to create simple cut window displayed on plot},
  //                "pico_branch" or NamedFunc cuts, sample proc to use, plotting style).Weight(weight NamedFunc).LuminosityTag(string with luminosity value).Tag("ShortName:description to include in plot name");
  pm.Push<Hist2D>(Axis(20,-4.7,4.7, "ll_eta[0]", "#eta_{ll}", {}),
                  Axis(20,-3.14,3.14, "ll_phi[0]", "#phi_{ll}", {}),
                  init && base_half, procs_test[i], style_2D).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot_example2D_ll_etaphi"+ procs_name[i]);

  //1D histogram example
  //pm.Push<Hist1D>(Axis(bins,low edge,high edge , "pico_branch" or NamedFunc, "Axis label" ,{Paired values to create simple cut window displayed on plot}), "pico_branch" or NamedFunc cuts, sample proc to use, plotting style).Weight(weight NamedFunc).LuminosityTag(string with luminosity value).Tag("ShortName:description to include in plot name");
  pm.Push<Hist1D>(Axis(40,70,110 , "ll_m[0]"   ,"m_{ll} [GeV]"   ,{80,100}), init && base_half, procs_test[i], style_1D).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot_example1D_mll" + procs_name[i]);
  pm.Push<Hist1D>(Axis(40,70,110 , "ll_m[0]"   ,"m_{ll} [GeV]"   ,{80,100}), init && base_half, procs_test[i], style_1D).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).EnergyLabel(e_label[0]).Tag("ShortName:plot_example1D_mll_label" + procs_name[i]); //Same plots but with custom energy and lumi label
}

  //This chunk of code is needed to actually create the plots. The 1.0 can be changed to a luminosity value. 
  pm.min_print_ = true;
  pm.MakePlots(1.01);
}


