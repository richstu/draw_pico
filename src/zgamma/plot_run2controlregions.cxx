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
#include "zgamma/zg_utilities.hpp"
#include "zgamma/zg_functions.hpp"
#include "zgamma/categorization_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace ZgUtilities;
using namespace ZgFunctions;

NamedFunc hasNanoPhoton("hasNanoPhoton",[](const Baby &b) -> NamedFunc::ScalarType{ return b.photon_pt() -> size() > 0; });

int main() {
  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");

  //Creating the vector of samples to process1
  vector<shared_ptr<Process>> procs = ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_zgamma.txt","AllUnskimmedTrigsData");

  //Declaring the plot option we want to use
  PlotOpt lin_dnorm("txt/plot_styles.txt","CMSPaper");
  lin_dnorm.Title(TitleType::info)
           .YAxis(YAxisType::linear)
           .Stack(StackType::data_norm)
           .Overflow(OverflowType::both)
           .YTitleOffset(1.)
           .AutoYAxis(false)
           .UseCMYK(false)
           .LeftMargin(0.17)
           .LegendColumns(1)
           .CanvasWidth(1077)
           .CanvasWidth(900)
           .FileExtensions({"pdf"})
           .Bottom(BottomType::ratio);

  vector<PlotOpt> ops_dwrat = { lin_dnorm.Bottom(BottomType::ratio) };


  cout << "db1" << endl;

  //Required selections for Control regions
  NamedFunc r2base = ZgFunctions::run2_baseline && "llphoton_m[0] < 122 || llphoton_m[0] > 128";
  NamedFunc tightbase = ZgFunctions::tightened_baseline && "llphoton_m[0] < 122 || llphoton_m[0] > 128";
  NamedFunc nqualityCut = "(photon_idmva[0] < -0.4 && (photon_eta[0] < 1.4442 && photon_eta[0] >-1.4442) ) || (photon_idmva[0] < -0.58 && (photon_eta[0] < -1.566 || photon_eta[0] > 1.566) && (photon_eta[0] < 2.5&& photon_eta[0] > -2.5 ) )";


  //Various control regions
  //const std::vector<NamedFunc> vector_run2_baseline      = {"nll>0", pass_trigs_and_pt, "ll_m[0] > 50", "nphoton > 0",  pTy_mlly > 15.0/110, mll_mlly > 185, mlly > 100 && mlly < 180};
  vector<NamedFunc> ControlRegions = { r2base, Nreverse1(ZgFunctions::vector_run2_baseline, 4),
                                       Nminusk(ZgFunctions::vector_run2_baseline,{5,6}) && mll_mlly < 185 && mlly < 100,
                                       r2base,//hasNanoPhoton && Nminus1(ZgFunctions::vector_run2_baseline,3) && "photon_pt[0] > 15" && nqualityCut,
                                       r2base && CatUtilities::run2_ggF, r2base && CatUtilities::run2_VBF, r2base && CatUtilities::run2_lep, 
                                       tightbase && CatUtilities::run2_ggF, tightbase && CatUtilities::run2_VBF, tightbase && CatUtilities::run2_lep
                                     };
  vector<string> CRString  = {"_sideband","_inv_pTgamma_mlly","_fsr_Zpeak","_low_idmva","_sideband_ggF","_sideband_VBF" , "_sideband_lep",
                             "_tight_sb_ggF","_tight_sb_VBF","_tight_sb_lep"};

  cout << "db2" << endl;

  //Adds flavor splits to the control regions
  vector<NamedFunc> FlavorCut = {"ll_lepid[0]==11" , "ll_lepid[0]==13", "1"};
  vector<int>       FlavorInt = {11,13,0};
  vector<string>    FlavString= {"_ee","_mumu","_ll"};

  //Loops through flavor and control regions and makes respective plots
  //List of CR plots can be found in zg_functions.cpp
  PlotMaker pm;
  NamedFunc curr_sel = "1"; string curr_string = ""; int ll_flavor = 0;
  for(unsigned int idx_CR(0); idx_CR < ControlRegions.size(); idx_CR++){
    for(unsigned int idx_fl(0); idx_fl < FlavorCut.size(); idx_fl++){
      //Sets strings and parameters to make function calls more clear
      curr_sel = ControlRegions[idx_CR] && FlavorCut[idx_fl];
      curr_string = CRString[idx_CR] + FlavString[idx_fl];
      ll_flavor = FlavorInt[idx_fl];

      //Adds plots to plotmaker via functions in zg_functions
      if(idx_CR<5 || idx_CR==7){  ggF_controlregion_plots(pm, curr_sel, procs, ops_dwrat, wgt_1invfb, curr_string, ll_flavor);}
      if(idx_CR==5 || idx_CR==8){ VBF_controlregion_plots(pm, curr_sel, procs, ops_dwrat, wgt_1invfb, curr_string, ll_flavor);}
      if(idx_CR==6 || idx_CR==9){ Lep_controlregion_plots(pm, curr_sel, procs, ops_dwrat, wgt_1invfb, curr_string, ll_flavor);}
    }
  }

  cout << "db3" << endl;


  pm.min_print_ = true;
  pm.MakePlots(137.61);
}


