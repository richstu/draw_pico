#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <bitset>
#include <unordered_map>
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
using namespace CatUtilities;

 NamedFunc hasNanoPhoton("hasNanoPhoton",[](const Baby &b) -> NamedFunc::ScalarType{ return b.photon_pt() -> size() > 0; });
 NamedFunc hasNll("hasNll",[](const Baby &b) -> NamedFunc::ScalarType{ return b.nll() > 0; });

  bool contains(string string1,string string2){
    return string1.find(string2,0)!=std::string::npos ? true : false;
  }

int main(int argc, char *argv[]) {

  //Determines whether or not to plot all control regions. 
  //To plot all give the input of "all" as the first input argument
  bool   plot_all    = false;
  string plot_option = "";

  if(argc > 1){
    plot_option = argv[1];
  }
  if(plot_option == "all"){
    plot_all = true;
    std::cout << "Selected option to plot all control regions. Note this will take longer than just plotting the sideband control region." << std::endl; 
  } else {
    std::cout << "Selected default to print just sideband control region. Using skimmed picos." << std::endl;
  }

  //Select a single category to run the plots over.
  int category = -1;
  if(argc > 2){
    category = atoi(argv[2]);
  }

  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");
  NamedFunc type("rtype",  [](const Baby &b) -> NamedFunc::ScalarType{ return b.type(); });

  //Declares the samples  
  vector<shared_ptr<Process>> procs = ZgUtilities::procs_with_sig_scale("txt/samples_zgamma.txt","AllMoreTrigsData",10);
  //vector<shared_ptr<Process>> procs = ZgUtilities::procs_with_sig_scale("txt/samples_zgamma.txt","MinimalMoreTrigsData",10);

  if(plot_all){
    procs = ZgUtilities::procs_with_sig_scale("txt/samples_zgamma.txt","AllSkimLLMoreTrigsData",10);
  }

  PlotOpt lin_lumi("txt/plot_styles.txt","Std1D");
  lin_lumi.Title(TitleType::info)
          .YAxis(YAxisType::linear)
          .Stack(StackType::data_norm)
          .Overflow(OverflowType::both)
          .YTitleOffset(1.)
          .LogMinimum(0.001)
          .AutoYAxis(true)
          .UseCMYK(false)
          .LeftMargin(0.17)
          .LegendColumns(1)
          .CanvasWidth(1077)
          .CanvasWidth(900)
          .FileExtensions({"pdf"})
          .Bottom(BottomType::ratio);

  PlotOpt style("txt/plot_styles.txt","Eff2D");
  vector<PlotOpt> bkg_hist = {style().Stack(StackType::data_norm)
                                     .YTitleOffset(1.)
                                     .LogMinimum(0.1)
                                     .Overflow(OverflowType::both)
                                     .CanvasWidth(600)
                                     .LabelSize(0.04)
                                     .UseCMYK(true)  
                                     .YAxis(YAxisType::log)
                                     .Title(TitleType::info)};
  vector<PlotOpt> ops = {lin_lumi};



  NamedFunc tightbase  = ZgFunctions::tightened_baseline && "llphoton_m[0] < 122 || llphoton_m[0] > 128";
  NamedFunc Nllgamma_g0_minsel = hasNll && hasNanoPhoton;
  //Various control regions. Default option only plots the sideband region
  vector<NamedFunc> controlregions_vec = { tightbase };
  vector<string> cr_str_vec        = {"_sideband"};
  if(plot_all==true){
    controlregions_vec = {tightbase, Nllgamma_g0_minsel && Nreverse1(ZgFunctions::vector_tightened_baseline, 4), 
                      Nllgamma_g0_minsel && Nminusk(ZgFunctions::vector_tightened_baseline,{5,6,7}) && "ll_m[0] < 80 && llphoton_m[0] < 100",
                      Nllgamma_g0_minsel && Nreplace1(ZgFunctions::vector_tightened_baseline, hasNanoPhoton && "photon_pt[0] > 15 && photon_drmin[0] > 0.3 && !photon_id80[0]", 3)};  
    cr_str_vec     = {"_sideband","_ptyrat_l15o110","_Zfsrpeak","_lowIDMVA"};

  }


  //Map is used to see if a category has a specific control region binning. The keys here are 10*category + idx_of_control_region
  //The map itself does not contain the selections because a map<int,NamedFunc> was complaining about trying to access a private contructor

  std::vector<NamedFunc> map_pair_selections = {"met > 50", (pTlly_mlly < 0.4), (l3_mini > 0.15), "dijet_m > 800"};

  //move this to tuple
  std::unordered_map<int,int>    category_specific_control_regions = {{2, 0}, {0, 0}, {12, 0}, {22, 1}, {32, 2}, {42, 3}, {40, 3}};
  std::unordered_map<int,int>    category_specific_control_regions_whbs = {{2, -1}, {0, -1}, {12, -1}, {22, 0}, {32, 0}, {42, -1}, {40, -1}};
  std::unordered_map<int,string> category_specific_control_regions_names = {{2, "metg50"},
                                                             {0,  "metg50"},
                                                             {12, "metg50"},
                                                             {22, "pTlly_mlly_l0p4"},
                                                             {32, "l3_mini_g0p15"},
                                                             {42, "dijet_m_g800"},
                                                             {40, "dijet_m_g800"}};
 

  //Standard lepton flavor splits and string for each split
  vector<NamedFunc> lep_flavor = {"ll_lepid[0] == 11","ll_lepid[0] == 13","1"};
  vector<string>    lep_str    = {"_ee","_mumu","_ll"};

  cout << "Making vector with each control region" << endl;

  PlotMaker pm;

  //Makes 1D vector for plotting. As long as inner most loop is the categories can take %Ncategories to get correct plots for each category
  std::vector<NamedFunc> NamedFunc_vector = {};
  std::vector<string>    string_vector    = {};
  NamedFunc control_region_selection = "1";
  NamedFunc category_specific_control_region_selection = "1";

  string control_region_name = "";
  string category_specific_control_region_name = "";

  //Loop that handles the addition of all the different control regions to plot
  for(unsigned int idx_cr = 0; idx_cr < controlregions_vec.size(); idx_cr++){
    for(unsigned int idx_lep = 0; idx_lep < lep_flavor.size(); idx_lep++){

      //Some initial variables to set. Basically setting general control region selections  
      //Note: The control_region_name portion here is actually the ending of the string that is used. 
      //This is because the category should be before the control region.
      control_region_selection = controlregions_vec[idx_cr] && lep_flavor[idx_lep];
      control_region_name =  lep_str[idx_lep] + cr_str_vec[idx_cr];

      //---------------------------------------------------------------------------------------------------------//
      //------------------- Portion for plotting all categories' general control regions ------------------------//
      //---------------------------------------------------------------------------------------------------------//

      if(category == -1){
        for(unsigned int idx_cat = 0; idx_cat < run3_category_vector.size(); idx_cat++){
          control_region_name = "an_controlregions_" + run3_category_labels[idx_cat] + control_region_name; 

          //FSR photon control region requires the categories to have the mll selections not applied.
          if(idx_cr==2){
            control_region_selection = control_region_selection && CatUtilities::run3_catwsel_nomll_vector[idx_cat];
          //Otherwise set the selection to the full category specific selection
          } else {
            control_region_selection = control_region_selection && CatUtilities::run3_catwsel_vector[idx_cat]; 
          }

          //Adding this control region to the ones that we want to plot
          NamedFunc_vector.push_back(control_region_selection);
          string_vector.push_back(control_region_name);

          //Some messages for debugging
          //cout << "NAME: " << control_region_name << endl;
          //cout << "SELECTION: " << control_region_selection << endl;

          //Have to reset the name in order to addition of previous category names to current category
          control_region_name = lep_str[idx_lep] + cr_str_vec[idx_cr];
          control_region_selection = controlregions_vec[idx_cr] && lep_flavor[idx_lep];
        }
        
        //After adding each control region's plots continuing to not add any other control regions
        continue;
      }

      //---------------------------------------------------------------------------------------------------------//
      //-------------------- Portion for plotting a specific category's control regions -------------------------//
      //---------------------------------------------------------------------------------------------------------//

      //FSR photon control region requires the categories to have the mll selections not applied.
      if(idx_cr==2){
        control_region_selection = control_region_selection && CatUtilities::run3_catwsel_nomll_vector[category]; 
      //Otherwise set the selection to the full category specific selection
      } else {
        control_region_selection = control_region_selection && CatUtilities::run3_catwsel_vector[category]; 
      }
      control_region_name = "an_controlregions_" + run3_category_labels[category] + control_region_name; 

      //Add the general control region
      NamedFunc_vector.push_back(control_region_selection);
      string_vector.push_back(control_region_name);

      //Some messages for debugging
      //cout << "NAME: " << control_region_name << endl;
      //cout << "SELECTION: " << control_region_selection << endl;

      //Checking if there is a category specific control region split
      //NOTE: This will be run ONLY IF a specific category has been specified AND the plot_all option has been selected
      //Running in this way would look like /run/zgamma/plot_an_controlregions.exe all 3
      int control_region_key = 10*category + static_cast<int>(idx_cr);
      if(plot_all && category_specific_control_regions.count(control_region_key)){

        //Check if control_region_key matches that one selection was reversed
        if(category_specific_control_regions_whbs[control_region_key] > -1){
          control_region_selection = controlregions_vec[idx_cr] && CatUtilities::run3_category_vector[category];
          control_region_selection = control_region_selection && ZgFunctions::Nminus1(categories_selections_vector[category],category_specific_control_regions_whbs[control_region_key]);
        }

        //Properly sets the additions to the control region for each selection
        category_specific_control_region_selection = map_pair_selections[category_specific_control_regions[control_region_key]];
        category_specific_control_region_name      = category_specific_control_regions_names[control_region_key];

        //Adding the control regions based on the cat. specific selection
        NamedFunc_vector.push_back(control_region_selection && category_specific_control_region_selection);
        string_vector.push_back(control_region_name + category_specific_control_region_name);

        NamedFunc_vector.push_back(control_region_selection && !category_specific_control_region_selection);
        string_vector.push_back(control_region_name + "_NOT_" + category_specific_control_region_name);

        //Some messages for debugging
        cout << "NAME: " << control_region_name + category_specific_control_region_name << endl;
        cout << "SELECTION: " << (control_region_selection && category_specific_control_region_selection) << endl;
        cout << "NAME: " << control_region_name + "_NOT_" + category_specific_control_region_name << endl;
        cout << "SELECTION: " << (control_region_selection && !category_specific_control_region_selection) << endl;

        control_region_selection = controlregions_vec[idx_cr] && lep_flavor[idx_lep];
        control_region_name =  lep_str[idx_lep] + cr_str_vec[idx_cr];

      }

    }
  }

  cout << "Adding plots to PlotMaker" << endl;

  //Main loop that handles all of the plot making
  NamedFunc selection   = "1"; 
  string    labels      = "";
  int       Ncategories = 6;
  if(category==-1){
    for(unsigned int idx_sel = 0; idx_sel < NamedFunc_vector.size(); idx_sel++){
      selection = NamedFunc_vector[idx_sel];
      labels = string_vector[idx_sel];
 
      if(idx_sel%Ncategories==0){ CatUtilities::ttH_lep_controlregion_plots(pm,selection,procs,ops,wgt,labels); continue;}
      if(idx_sel%Ncategories==1){ CatUtilities::ttH_had_controlregion_plots(pm,selection,procs,ops,wgt,labels); continue;}
      if(idx_sel%Ncategories==2){  CatUtilities::ZH_MET_controlregion_plots(pm,selection,procs,ops,wgt,labels); continue;}
      if(idx_sel%Ncategories==3){   CatUtilities::WH_3l_controlregion_plots(pm,selection,procs,ops,wgt,labels); continue;}
      if(idx_sel%Ncategories==4){     CatUtilities::VBF_controlregion_plots(pm,selection,procs,ops,wgt,labels); continue;}
      if(idx_sel%Ncategories==5){     CatUtilities::ggF_controlregion_plots(pm,selection,procs,ops,wgt,labels); continue;}
    }
  } else {
     for(unsigned int idx_sel = 0; idx_sel < NamedFunc_vector.size(); idx_sel++){
      selection = NamedFunc_vector[idx_sel];
      labels = string_vector[idx_sel];
 
      if(category==0){ CatUtilities::ttH_lep_controlregion_plots(pm,selection,procs,ops,wgt,labels); continue;}
      if(category==1){ CatUtilities::ttH_had_controlregion_plots(pm,selection,procs,ops,wgt,labels); continue;}
      if(category==2){  CatUtilities::ZH_MET_controlregion_plots(pm,selection,procs,ops,wgt,labels); continue;}
      if(category==3){   CatUtilities::WH_3l_controlregion_plots(pm,selection,procs,ops,wgt,labels); continue;}
      if(category==4){     CatUtilities::VBF_controlregion_plots(pm,selection,procs,ops,wgt,labels); continue;}
      if(category==5){     CatUtilities::ggF_controlregion_plots(pm,selection,procs,ops,wgt,labels); continue;}
    }
  }
  cout << "Running over event loop" << endl;

  //Setting luminosity equal to Run 2. Right now this makes control regions only for inclusive Run 2 data
  pm.min_print_ = true;
  pm.SetLuminosityTag("137.61").MakePlots(1);
}


