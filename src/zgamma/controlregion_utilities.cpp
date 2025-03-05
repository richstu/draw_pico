#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <bitset>
#include <unistd.h>
#include <getopt.h>

#include "TError.h"
#include "TColor.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "core/baby.hpp"
#include "core/named_func.hpp"
#include "core/table_row.hpp"
#include "core/hist1d.hpp"
#include "core/utilities.hpp"
#include "core/named_func_utilities.hpp"
#include "core/plot_opt.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "zgamma/zg_utilities.hpp"
#include "zgamma/zg_functions.hpp"
#include "zgamma/categorization_utilities.hpp"
#include "zgamma/controlregion_utilities.hpp"

using namespace NamedFuncUtilities;
using namespace ZgUtilities;
using namespace CatUtilities;
using namespace ZgFunctions;

namespace CRUtilities{

  //Boolean used to select which era to use for plotting. False == Run2, True == Run3
  bool select_era(int argc, char *argv[]){
    bool era = false;
    if(argc<2){return era;}
    
    std::string era_sel = argv[1];
    if(era_sel == "run3" || era_sel =="r3"){era=true;}

    return era; 
  }

  const NamedFunc hasNanoPhoton("hasNanoPhoton",[](const Baby &b) -> NamedFunc::ScalarType{ return b.photon_pt() -> size() > 0; });
  const NamedFunc hasNll("hasNll",[](const Baby &b) -> NamedFunc::ScalarType{ return b.nll() > 0; });
  const NamedFunc Nllgamma_g0_minsel = hasNll && hasNanoPhoton;


  const NamedFunc higgs_window = "llphoton_m[0] < 120 || llphoton_m[0] > 130";
  const NamedFunc cr_sideband = ZgFunctions::tightened_baseline && !higgs_window;
  const NamedFunc cr_fake_photon = Nllgamma_g0_minsel && Nreplace1(ZgFunctions::vector_tightened_baseline, hasNanoPhoton && "photon_pt[0] > 15 && photon_drmin[0] > 0.3 && !photon_id80[0]",3);
  const NamedFunc cr_Zfsr_photon = Nllgamma_g0_minsel && Nminusk(ZgFunctions::vector_tightened_baseline,{5,6,7}) && "ll_m[0] < 80 && llphoton_m[0] < 100";
  const NamedFunc cr_soft_photon = Nllgamma_g0_minsel && Nreverse1(ZgFunctions::vector_tightened_baseline, 4);

  const std::vector<NamedFunc>   controlregions_vec     = {cr_sideband,   cr_soft_photon, cr_Zfsr_photon, cr_fake_photon};  
  const std::vector<std::string> controlregions_str_vec = {"_sideband","_ptyrat_l15o110",    "_Zfsrpeak", "_lowIDMVA"};

  const NamedFunc cr_sideband_refit = ZgFunctions::tightened_baseline_refit && !higgs_window;
  const NamedFunc cr_fake_photon_refit = Nllgamma_g0_minsel && Nreplace1(vtb_refit, hasNanoPhoton && "photon_pt[0] > 15 && photon_drmin[0] > 0.3 && !photon_id80[0]",3);
  const NamedFunc cr_Zfsr_photon_refit = Nllgamma_g0_minsel && Nminusk(vtb_refit,{5,6,7}) && "ll_m[0] < 80 && llphoton_refit_m < 100";
  const NamedFunc cr_soft_photon_refit = Nllgamma_g0_minsel && Nreverse1(vtb_refit, 4);

  const std::vector<NamedFunc> controlregions_refit_vec = {cr_sideband_refit, cr_soft_photon_refit, cr_Zfsr_photon_refit, cr_fake_photon_refit};  




  //This code defines additional control regions split by the selections listed below
  const std::vector<NamedFunc> map_pair_selections = {"met > 50", (pTlly_mlly < 0.4), (l3_mini > 0.15), "dijet_m > 800"};
  const std::unordered_map<int,int>    category_specific_control_regions = {{2, 0}, {0, 0}, {12, 0}, {22, 1}, {32, 2}, {42, 3}, {40, 3}};
  const std::unordered_map<int,int>    category_specific_control_regions_whbs = {{2, -1}, {0, -1}, {12, -1}, {22, 0}, {32, 0}, {42, -1}, {40, -1}};
  const std::unordered_map<int,std::string> category_specific_control_regions_names = {{2, "metg50"},         {0,  "metg50"},       {12, "metg50"},       {22, "pTlly_mlly_l0p4"}, 
                                                                            {32, "l3_mini_g0p15"}, {42, "dijet_m_g800"}, {40, "dijet_m_g800"}, {10, "metg50"}}; 

  std::vector<std::shared_ptr<Process>> control_region_procs(std::string year_select, int category){
    std::vector<std::shared_ptr<Process>> procs = ZgUtilities::procs_with_sig_scale("txt/samples_zgamma.txt","AllMoreTrigsData",10); 
    std::unordered_map<int,std::string> run2_skims = {{0, "Run2SkimttHlep"}, {1, "Run2SkimttHhad"}, {2, "Run2SkimZHMET"}, {3, "Run2SkimWH3l"}, {4, "Run2SkimVBF"},{5, "Run2SkimggF"}};
    std::unordered_map<int,std::string> run3_skims = {{0, "Run3SkimttHlep"}, {1, "Run3SkimttHhad"}, {2, "Run3SkimZHMET"}, {3, "Run3SkimWH3l"}, {4, "Run3SkimVBF"},{5, "Run3SkimggF"}};
    if(year_select=="run3"){
      if(category>-1 && category<6){ 
       procs = ZgUtilities::procs_with_sig_scale("txt/samples_zgamma_skims.txt", run3_skims[category],10);
       std::cout << "Using Run 3 category skims" << std::endl; 
      } else {
       procs = ZgUtilities::procs_with_sig_scale("txt/samples_zgamma.txt","AllMoreTrigsDataRun3",10);
       std::cout << "Using Run 3 LL skims" << std::endl; 
      }
    } else if(Contains(year_select,"20")){
      procs = ZgUtilities::procs_with_sig_scale("txt/samples_zgamma.txt","AllMoreTrigsData" + year_select,10); 
    } else if(category > -1){
      procs = ZgUtilities::procs_with_sig_scale("txt/samples_zgamma_skims.txt", run2_skims[category],10);
      std::cout << "Using Run 2 category skims" << std::endl; 
    } else{
       std::cout << "Using Run 2 LL skims" << std::endl; 
    }

    return procs;
  }

  std::string lumi_label(std::string year_select){
    std::unordered_map<std::string,std::string> lumi_select = {{"2016APV","16.80"},{"2016","19.51"},{"2017","41.48"},{"2018","8.17"},{"2022","8.17"},{"2022EE","27.01"},{"2023","17.61"},
                                                             {"2023BPix","9.53"},{"run2","137.61"},{"run3","62.32"}};
    return lumi_select[year_select];
  }

  //Function that can be called to return the plotting options you want from preset list
  std::vector<PlotOpt> an_plotting_options(std::string select_options){
    std::vector<PlotOpt> ops = {};

    //WAIT WHY AM I DOING THIS JUST PUT IT IN PLOT_STYLES
    if(select_options == "control_region"){
        PlotOpt an_cr_plot("txt/plot_styles.txt","CMSPaper");
        an_cr_plot.YAxis(PlotOptTypes::YAxisType::linear)
                .Stack(PlotOptTypes::StackType::data_norm)
                .YTitleOffset(1.75)
                .AutoYAxis(false)
                .UseCMYK(false)
                .LeftMargin(0.17)
                .CanvasWidth(800)
                .CanvasHeight(800)
                .FileExtensions({"pdf"});
        ops = {an_cr_plot};
    } else if(select_options=="2D"){
        PlotOpt an_2D_hist("txt/plot_styles.txt","Eff2D");
        an_2D_hist().Stack(PlotOptTypes::StackType::data_norm)
             .YTitleOffset(1.)
             .LabelSize(0.04)
             .UseCMYK(true); 
        ops = {an_2D_hist};
    } else if(select_options == "1D root"){
        PlotOpt an_cr_plot("txt/plot_styles.txt","CMSPaper");
        an_cr_plot.YAxis(PlotOptTypes::YAxisType::linear)
                .YTitleOffset(1.75)
                .AutoYAxis(false)
                .UseCMYK(false)
                .LeftMargin(0.17)
                .CanvasWidth(800)
                .CanvasHeight(800)
                .FileExtensions({"pdf","root"});
        ops = {an_cr_plot};
    } else if(select_options == "soverb"){
      PlotOpt soverb("txt/plot_styles.txt","CMSPaper");
      soverb.YAxis(PlotOptTypes::YAxisType::log)
          .Stack(PlotOptTypes::StackType::signal_overlay)
          .Overflow(PlotOptTypes::OverflowType::both)
          .YTitleOffset(1.)
          .LogMinimum(0.001)
          .AutoYAxis(true)
          .UseCMYK(false)
          .LeftMargin(0.17)
          .LegendColumns(1)
          .CanvasWidth(900)
          .Bottom(PlotOptTypes::BottomType::sorb)
          .FileExtensions({"pdf"});
      ops = {soverb};
    } else if(select_options == "soverb upper"){
      PlotOpt soverb_upper("txt/plot_styles.txt","CMSPaper");
      soverb_upper.Title(PlotOptTypes::TitleType::info)
          .YAxis(PlotOptTypes::YAxisType::log)
          .Stack(PlotOptTypes::StackType::signal_overlay)
          .Overflow(PlotOptTypes::OverflowType::both)
          .YTitleOffset(1.)
          .LogMinimum(0.001)
          .AutoYAxis(true)
          .UseCMYK(false)
          .LeftMargin(0.17)
          .LegendColumns(1)
          .CanvasWidth(900)
          .Bottom(PlotOptTypes::BottomType::sorb_cut_upper)
          .FileExtensions({"pdf"});
      ops = {soverb_upper};
    } else {
        PlotOpt an_1D_plot("txt/plot_styles.txt","CMSPaper");
        an_1D_plot.YAxis(PlotOptTypes::YAxisType::linear)
                .YTitleOffset(1.75)
                .AutoYAxis(false)
                .UseCMYK(false)
                .LeftMargin(0.17)
                .CanvasWidth(800)
                .CanvasHeight(800)
                .FileExtensions({"pdf"});
         ops = {an_1D_plot};
    }
 
    return ops;
 }

/*
std::tuple<std::vector<NamedFunc>,std::vector<std::string>> create_control_region_selections_and_labels(std::string label, bool plot_all, int category){
  //Standard lepton flavor splits and string for each split
  std::vector<NamedFunc> lep_flavor = {"ll_lepid[0] == 11","ll_lepid[0] == 13","1"};
  std::vector<std::string>    lep_str    = {"_ee","_mumu","_ll"};

  std::cout << "Making vector with each control region" << std::endl;

  std::vector<std::string> temp_controlregions_str_vec = controlregions_str_vec;
  for(unsigned int idx = 0; idx < temp_controlregions_str_vec.size(); idx++){
    temp_controlregions_str_vec[idx] = temp_controlregions_str_vec[idx] + label;
  }

  //Makes 1D vector for plotting. As long as inner most loop is the categories can take %Ncategories to get correct plots for each category
  std::vector<NamedFunc> NamedFunc_vector = {};
  std::vector<std::string>    string_vector    = {};
  NamedFunc control_region_selection = "1";
  NamedFunc category_specific_control_region_selection = "1";

  std::string control_region_name = "";
  std::string category_specific_control_region_name = "";

  //Loop that handles the addition of all the different control regions to plot
  for(unsigned int idx_cr = 0; idx_cr < controlregions_vec.size(); idx_cr++){
    for(unsigned int idx_lep = 0; idx_lep < lep_flavor.size(); idx_lep++){

      //Some initial variables to set. Basically setting general control region selections  
      //Note: The control_region_name portion here is actually the ending of the string that is used. 
      //This is because the category should be before the control region.
      control_region_selection = controlregions_vec[idx_cr] && lep_flavor[idx_lep];
      control_region_name =  lep_str[idx_lep] + temp_controlregions_str_vec[idx_cr];

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
            //if(idx_cr>0){
            //  control_region_selection = control_region_selection && CatUtilities::run3_category_vector[idx_cat]; 
            //} else {
            //  control_region_selection = control_region_selection && CatUtilities::run3_catwsel_vector[idx_cat];
            //}
            //control_region_selection = control_region_selection && CatUtilities::run3_category_vector[idx_cat]; 
          }

          //Adding this control region to the ones that we want to plot
          NamedFunc_vector.push_back(control_region_selection);
          string_vector.push_back(control_region_name);

          //Some messages for debugging
          //std::cout << "NAME: " << control_region_name << std::endl;
          //std::cout << "SELECTION: " << control_region_selection << std::endl;

          //Have to reset the name in order to addition of previous category names to current category
          control_region_name = lep_str[idx_lep] + temp_controlregions_str_vec[idx_cr];
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
        if(category>-1){
          control_region_selection = control_region_selection && CatUtilities::run3_catbs_nomll_vector[category]; 
        } else {
          control_region_selection = control_region_selection && CatUtilities::run3_catwsel_nomll_vector[category]; 
        }
      //Otherwise set the selection to the full category specific selection
      } else {
        if(category>-1){
          control_region_selection = control_region_selection && CatUtilities::run3_catbs_vector[category]; 
        } else {
          control_region_selection = control_region_selection && CatUtilities::run3_catwsel_vector[category]; 
        }
      }
      control_region_name = "an_controlregions_" + run3_category_labels[category] + control_region_name; 

      //Add the general control region
      NamedFunc_vector.push_back(control_region_selection);
      string_vector.push_back(control_region_name);

      //Some messages for debugging
      //std::cout << "NAME: " << control_region_name << std::endl;
      //std::cout << "SELECTION: " << control_region_selection << std::endl;

      //Checking if there is a category specific control region split
      //NOTE: This will be run ONLY IF a specific category has been specified AND the plot_all option has been selected
      //Running in this way would look like /run/zgamma/plot_an_controlregions.exe all 3
      int control_region_key = 10*category + static_cast<int>(idx_cr);
      if(plot_all && category_specific_control_regions.count(control_region_key)){

        //Check if control_region_key matches that one selection was reversed
        if(category_specific_control_regions_whbs[control_region_key] > -1){
          control_region_selection = controlregions_vec[idx_cr] && CatUtilities::run3_category_vector[category];
          control_region_selection = control_region_selection && NamedFuncUtilities::Nminus1(categories_selections_vector[category],category_specific_control_regions_whbs[control_region_key]);
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
        //std::cout << "NAME: " << control_region_name + category_specific_control_region_name << std::endl;
        //std::cout << "SELECTION: " << (control_region_selection && category_specific_control_region_selection) << std::endl;
        //std::cout << "NAME: " << control_region_name + "_NOT_" + category_specific_control_region_name << std::endl;
        //std::cout << "SELECTION: " << (control_region_selection && !category_specific_control_region_selection) << std::endl;

        control_region_selection = controlregions_vec[idx_cr] && lep_flavor[idx_lep];
        control_region_name =  lep_str[idx_lep] + temp_controlregions_str_vec[idx_cr];

      }

    }
  }
  return make_tuple(NamedFunc_vector,string_vector);
}
*/


}
