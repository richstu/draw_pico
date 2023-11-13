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
//#include "TMVA/Reader.h"
#include "RooGenericPdf.h"
#include "RooRealVar.h"
#include "RooArgList.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/datacard.hpp"
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
using namespace CatUtilities;
using SelectionList = Datacard::SelectionList;
using Systematic = Datacard::Systematic;
  
int main() {
  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");
  NamedFunc type("rtype",  [](const Baby &b) -> NamedFunc::ScalarType{ return b.type(); });

  //Processes with dilepton triggers
  vector<shared_ptr<Process>> procs     = ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_zgamma.txt","All");

  //Plot options 
  PlotOpt lin_lumi("txt/plot_styles.txt","Std1D");
  lin_lumi.Title(TitleType::info)
          .YAxis(YAxisType::linear)
          .Stack(StackType::signal_overlay)
          .Overflow(OverflowType::both)
          .YTitleOffset(1.)
          .LogMinimum(0.001)
          .AutoYAxis(true)
          .UseCMYK(false)
          .LeftMargin(0.17)
          .LegendColumns(1)
          .CanvasWidth(1077)
          .CanvasWidth(900)
          .FileExtensions({"pdf"});
  vector<PlotOpt> ops = {lin_lumi, lin_lumi().YAxis(YAxisType::log)};

  PlotOpt style("txt/plot_styles.txt","Eff2D");
  vector<PlotOpt> bkg_hist = {style().Stack(StackType::data_norm)
                                     .YTitleOffset(1.)
                                     .LogMinimum(0.1)
                                     //.LogMaximum(2000)
                                     .Overflow(OverflowType::both)
                                     .CanvasWidth(600)
                                     .LabelSize(0.04)
                                     .UseCMYK(true)  
                                     .YAxis(YAxisType::log)
                                     .Title(TitleType::info)};
                           
  //Defining the categories that will be plotted 
  vector<string>    run3_cat_labels = CatUtilities::run3_category_labels;                          
  vector<NamedFunc> run3_cats       = CatUtilities::run3_category_vector;
     
  //Defining if we are plotting based on lepton splits 
  vector<NamedFunc> flavor_NamedFunc_vector = {"1"};
  vector<string> flavor_string_vector = {"_ll"};

  //Loops through categories and flavors to create a 1D vector
  vector<NamedFunc> NamedFunc_loop = {};
  vector<string> labels_loop = {};
  for(unsigned int idx_fl=0; idx_fl < flavor_NamedFunc_vector.size(); idx_fl++){
    for(unsigned int idx_cat = 0; idx_cat < run3_cats.size(); idx_cat++){
      NamedFunc_loop.push_back(ZgFunctions::tightened_baseline && run3_cats[idx_cat] && flavor_NamedFunc_vector[idx_fl]);        
      labels_loop.push_back(run3_cat_labels[idx_cat] + flavor_string_vector[idx_fl]);
    }
  }


PlotMaker pm;

//Some variables useful for plotting
int Ncategories = run3_cats.size();
NamedFunc selection = "1";
string labels = "";

//This loops through the categories and adds the respective categories' plots
for(unsigned int idx_sel = 0; idx_sel < NamedFunc_loop.size(); idx_sel++){
  //Selections and labels for each category
  selection = NamedFunc_loop[idx_sel];
  labels =    labels_loop[idx_sel] + "_MCONLY";

  //Each category with the respective plots
  if(idx_sel%Ncategories==0){ CatUtilities::ttH_lep_controlregion_plots(pm,selection,procs,ops,wgt_1invfb,labels);}
  if(idx_sel%Ncategories==1){   CatUtilities::VH_4l_controlregion_plots(pm,selection,procs,ops,wgt_1invfb,labels);}
  if(idx_sel%Ncategories==2){   CatUtilities::VH_3l_controlregion_plots(pm,selection,procs,ops,wgt_1invfb,labels);}
  if(idx_sel%Ncategories==3){ CatUtilities::ttH_had_controlregion_plots(pm,selection,procs,ops,wgt_1invfb,labels);}
  if(idx_sel%Ncategories==4){  CatUtilities::VH_2bl_controlregion_plots(pm,selection,procs,ops,wgt_1invfb,labels);}
  if(idx_sel%Ncategories==5){  CatUtilities::VH_MET_controlregion_plots(pm,selection,procs,ops,wgt_1invfb,labels);}
  if(idx_sel%Ncategories==6){  ZgFunctions::VBF_controlregion_plots(pm,selection,procs,ops,wgt_1invfb,labels);}
  if(idx_sel%Ncategories==7){  ZgFunctions::VBF_controlregion_plots(pm,selection,procs,ops,wgt_1invfb,labels);}
  if(idx_sel%Ncategories==8){  ZgFunctions::ggF_controlregion_plots(pm,selection,procs,ops,wgt_1invfb,labels);}

}

  pm.min_print_ = true;
  pm.MakePlots(200.0);
}

