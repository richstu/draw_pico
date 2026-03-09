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
#include "zgamma/controlregion_utilities.hpp"
#include "core/named_func_utilities.hpp"
using namespace std;
using namespace PlotOptTypes;
using namespace ZgUtilities;
using namespace ZgFunctions;
using namespace CatUtilities;
using namespace CRUtilities;
using namespace NamedFuncUtilities;

// NamedFunc hasNanoPhoton("hasNanoPhoton",[](const Baby &b) -> NamedFunc::ScalarType{ return b.photon_pt() -> size() > 0; });
// NamedFunc hasNll("hasNll",[](const Baby &b) -> NamedFunc::ScalarType{ return b.nll() > 0; });

  bool contains(string string1,string string2){
    return string1.find(string2,0)!=std::string::npos ? true : false;
  }

  NamedFunc jet_veto_event("jet_veto_event",[](const Baby &b) -> NamedFunc::ScalarType{
    for(unsigned int idx = 0; idx < b.jet_pt() -> size(); idx++){
      if(b.jet_isgood() -> at(idx) && fabs(b.jet_eta() -> at(idx)) > 2.5 && fabs(b.jet_eta() -> at(idx)) < 3.139 && b.jet_pt() ->at(idx) < 50){return false; }
    }
    return true;
  });

  NamedFunc njetcorr("njetcorr",[](const Baby &b) -> NamedFunc::ScalarType{
    int njet = 0;
    for(unsigned int idx = 0; idx < b.jet_pt() -> size(); idx++){
      if(b.jet_isgood() -> at(idx) && b.jet_pt() -> at(idx) < 50){njet++; }
    }
    return njet;
  });


int main(int argc, char *argv[]) {

  //Determines whether or not to plot all control regions. 
  //To plot all give the input of "all" as the first input argument
  bool   plot_all    = false;
  bool   run3_plots  = false;
  string plot_option = "";
  string year_select = "";

  if(argc > 1){
    plot_option = argv[1];
  }

  //Parse input options
  if(plot_option == "all"){
    plot_all = true;
    std::cout << "Selected option to plot all control regions. Note this will take longer than just plotting the sideband control region." << std::endl;
  } else if(plot_option == "all_run3"){
    plot_all = true;
    run3_plots = true;
    year_select = "run3";
    std::cout << "Selected plot option for all control regions with Run 3 Data/MC." << endl;
  } else if(plot_option == "year") {
    if(argc < 2){std::cout <<  "Please provide a year as the second option." << std::endl;return -1; }
    year_select = argv[2];
    if(contains(year_select,"2022") || contains(year_select,"2023")){run3_plots = true;}
    std::cout << "Selected the year option: " << year_select << endl;
  } else if(plot_option == "sideband_r2"){
    std::cout << "Selected default to print just sideband control region. Using skimmed picos." << std::endl; 
  } else if(plot_option == "sideband_r3"){
    run3_plots = true;
    year_select = "run3";
    std::cout << "Selected default to print just sideband control region for Run 3. Using skimmed picos." << std::endl;
  } else {
    std::cout << "Selected default to print just sideband control region. Using skimmed picos." << std::endl;
    //return -1;
  }

  //Select a single category to run the plots over.
  int category = -1;
  if(argc > 2 && plot_option!="year"){
    category = atoi(argv[2]);
  }

  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");
  NamedFunc type("rtype",  [](const Baby &b) -> NamedFunc::ScalarType{ return b.type(); });

  string lumi_label = CRUtilities::R2Lumi;
  //Declares the samples
  string samples_text_file = CRUtilities::SamplesFile;
  string samples_selector  = "AllMoreTrigsData";
  int signal_multiplier    = 10;
  if(year_select=="2016"){           samples_selector = "AllMoreTrigsData2016";    lumi_label = "19.51";
  } else if(year_select=="2016APV"){ samples_selector = "AllMoreTrigsData2016APV"; lumi_label = "16.80";
  } else if(year_select=="2017"){    samples_selector = "AllMoreTrigsData2017";    lumi_label = "41.48";
  } else if(year_select=="2018"){    samples_selector = "AllMoreTrigsData2018";    lumi_label = "59.83";
  } else if(year_select=="2022"){    samples_selector = "AllMoreTrigsData2022";    lumi_label = "8.17";
  } else if(year_select=="2022EE"){  samples_selector = "AllMoreTrigsData2022EE";  lumi_label = "27.01";
  } else if(year_select=="2023"){    samples_selector = "AllMoreTrigsData2023";    lumi_label = "17.61";
  } else if(year_select=="2023BPix"){samples_selector = "AllMoreTrigsData2023BPix";lumi_label = "9.53";
  } else if(year_select=="run3"){    samples_selector = "AllMoreTrigsDataRun3";    lumi_label = CRUtilities::R3Lumi;
  } else if(plot_option=="year"){
     std::cout << "Please enter a valid year (2016,2016APV,2017,2018)" << std::endl;
     return -1;
  } else {
      cout << "Using Run 2 samples" << endl; 
  }

  vector<shared_ptr<Process>> procs = ZgUtilities::procs_with_sig_scale(samples_text_file,samples_selector,signal_multiplier); 
 
 std::vector<PlotOpt> ops = an_plotting_options("control_region");
 

  NamedFunc tightbase  = ZgFunctions::tightened_baseline_pinnacles && "llphoton_m[0] < 120 || llphoton_m[0] > 130";
  //NamedFunc Nllgamma_g0_minsel = hasNll && hasNanoPhoton;
  //&& jet_veto_event
  // && "(nlep==2 && njet>=2 && nbdfm==0)" && jet_veto_event
  //Various control regions. Default option only plots the sideband region
  /*vector<NamedFunc> controlregions_vec = { tightbase};
  vector<string> cr_str_vec        = {"_sideband"};
 
  if(plot_all==true){ 
    controlregions_vec = CRUtilities::controlregions_vec; 
    cr_str_vec = CRUtilities::controlregions_str_vec;
  }
  
  if(plot_option=="year"){
    for(unsigned int idx_str = 0; idx_str < cr_str_vec.size(); idx_str++){
      cr_str_vec[idx_str] = cr_str_vec[idx_str] + "_" + year_select;  
    }
  } else if(run3_plots){
    for(unsigned int idx_str = 0; idx_str < cr_str_vec.size(); idx_str++){
      cr_str_vec[idx_str] = cr_str_vec[idx_str] + "_" + year_select;  
    }
  }
  */
  //
   if(plot_all==true){ cout << "Wow a useless boolean needs to be used to not throw an error. I wish I wasn't so lazy and would take the 30 seconds to remove all references of this bool instead!" << endl; } 

  cout << "Making vector with each control region" << endl;
  std::tuple<std::vector<NamedFunc>,std::vector<std::string>> tuple_control_regions = create_control_region_selections_and_labels(category);
  vector<NamedFunc> NamedFunc_vector = {};//get<0>(tuple_control_regions);
  vector<string> string_vector = {};//get<1>(tuple_control_regions); 
  
  std::vector<NamedFunc>   cr_vec     = {ZgFunctions::tightened_baseline && "llphoton_refit_m > 100 && llphoton_refit_m < 180" && !higgs_window_refit,
  //Nllgamma_g0_minsel && Nreplace1(ZgFunctions::vector_tightened_baseline, hasNanoPhoton && "photon_pt[0] > 15 && photon_drmin[0] > 0.3 && !photon_id80[0] && photon_id[0]",3),
  Nllgamma_g0_minsel && Nminusk(ZgFunctions::vector_tightened_baseline,{5,6,7}) && "ll_m[0] < 80 && llphoton_m[0] < 100",
  Nllgamma_g0_minsel && Nreverse1(ZgFunctions::vector_tightened_baseline, 4)};
  
  std::vector<std::string> controlregions_str_vec = {"_sideband", //"_lowIDMVA", 
                                                     "_Zfsrpeak", "_ptyrat_l15o110"};
  std::vector<NamedFunc>   lep_flavor = {"ll_lepid[0] == 11", "ll_lepid[0] == 13",   "1"};
  std::vector<std::string> lep_str    = {              "_ee",             "_mumu", "_ll"};

  string r3_str = run3_plots ? "_run3" : "";

  for(unsigned int idx_cr = 0; idx_cr < cr_vec.size(); idx_cr++){
    for(unsigned int idx_lep = 0; idx_lep < lep_flavor.size(); idx_lep++){
      if(idx_cr!=1){
       NamedFunc_vector.push_back(cr_vec[idx_cr] && lep_flavor[idx_lep] && run3_catwsel_vector[category] ); 
      } else{
        NamedFunc_vector.push_back(cr_vec[idx_cr] && lep_flavor[idx_lep] && run3_catwsel_nomll_vector[category] );
 
      }
      string_vector.push_back("an_controlregions_"  + run3_category_labels[category] + lep_str[idx_lep] + controlregions_str_vec[idx_cr] + r3_str ); 
    }
  }
  PlotMaker pm;

  cout << "Adding plots to PlotMaker" << endl;
  std::vector<std::vector<float>> vec_nbins = plotting_bins(run3_plots);

  //Main loop that handles all of the plot making
  NamedFunc selection   = "1"; 
  string    labels      = "";
  int       Ncategories = 7;
  vector<TableRow> event_yields= {};

  //pm.Push<Hist1D>(Axis(100, 100, 180, "llphoton_refit_m", "m_{ll#gamma} [GeV]", {}), Nllgamma_g0_minsel && Nreplace1(ZgFunctions::vector_tightened_baseline, hasNanoPhoton && "photon_pt[0] > 15 && photon_drmin[0] > 0.3 && !photon_id80[0]",3), procs, ops).Weight(wgt).Tag("ShortName:test_" + labels);
  //pm.Push<Hist1D>(Axis(100, 100, 180, "llphoton_refit_m", "m_{ll#gamma} [GeV]", {}), controlregions_vec[0], procs, ops).Weight(wgt).Tag("ShortName:test_" + labels);



   if(category==-1){
    for(unsigned int idx_sel = 0; idx_sel < NamedFunc_vector.size(); idx_sel++){
      selection = NamedFunc_vector[idx_sel];
      labels = string_vector[idx_sel];
      int nbins = vec_nbins[idx_sel%7][idx_sel/7];

      //Debugging message checking binning 
      //cout<< "for control region: " << labels << " using nbins= " << to_string(nbins) << endl;

      //CatUtilities::sample_kinrefit_plots(pm, selection, procs, ops, wgt, nbins, labels);
      CatUtilities::mvp_objects(          pm, selection, procs, ops, wgt, labels + run3_category_labels[idx_sel/3]);
      if(idx_sel%Ncategories==0){ CatUtilities::ttH_lep_controlregion_plots(pm,selection,procs,ops,wgt,labels, nbins, run3_plots); continue;}
      if(idx_sel%Ncategories==1){ CatUtilities::ttH_had_controlregion_plots(pm,selection,procs,ops,wgt,labels, nbins, run3_plots); continue;}
      if(idx_sel%Ncategories==2){  CatUtilities::ZH_MET_controlregion_plots(pm,selection,procs,ops,wgt,labels, nbins); continue;}
      if(idx_sel%Ncategories==3){   CatUtilities::WH_3l_controlregion_plots(pm,selection,procs,ops,wgt,labels, nbins); continue;}
      if(idx_sel%Ncategories==4){     CatUtilities::VBF_controlregion_plots(pm,selection,procs,ops,wgt,labels, nbins); continue;}
      if(idx_sel%Ncategories==5){     CatUtilities::ggF_controlregion_plots(pm,selection,procs,ops,wgt,labels, nbins); continue;}
      
      //Untagged CR
      if(idx_sel%Ncategories==6){     CatUtilities::ggF_controlregion_plots(pm,selection,procs,ops,wgt,labels, nbins); continue;}
 
  }

    event_yields.push_back(TableRow(labels, selection, 0, 0, wgt));    

  } else {
     for(unsigned int idx_sel = 0; idx_sel < NamedFunc_vector.size(); idx_sel++){
      selection = NamedFunc_vector[idx_sel];
      labels = string_vector[idx_sel];
      int nbins = vec_nbins[category][idx_sel%3];
      //cout<< "for control region: " << labels << " using nbins= " << to_string(nbins) << endl;

      CatUtilities::sample_kinrefit_plots(pm, selection, procs, ops, wgt, nbins, labels);
      CatUtilities::mvp_objects(          pm, selection, procs, ops, wgt, labels);
      if(category==0){ CatUtilities::ttH_lep_controlregion_plots(pm,selection,procs,ops,wgt,labels,nbins,run3_plots); continue;}
      if(category==1){ CatUtilities::ttH_had_controlregion_plots(pm,selection,procs,ops,wgt,labels,nbins,run3_plots); continue;}
      if(category==2){  CatUtilities::ZH_MET_controlregion_plots(pm,selection,procs,ops,wgt,labels,nbins); continue;}
      if(category==3){   CatUtilities::WH_3l_controlregion_plots(pm,selection,procs,ops,wgt,labels,nbins); continue;}
      if(category==4){     CatUtilities::VBF_controlregion_plots(pm,selection,procs,ops,wgt,labels,nbins); continue;}
      if(category==5){     CatUtilities::ggF_controlregion_plots(pm,selection,procs,ops,wgt,labels,nbins); continue;}

      //Untagged CR
      if(category==6){     CatUtilities::ggF_controlregion_plots(pm,selection,procs,ops,wgt,labels,nbins); continue;}
     }
  }
  cout << "Running over event loop" << endl;

  pm.Push<Table>("event_yields_" + year_select, event_yields, procs, false,true,false,false,false,true);

  //Setting luminosity based on lumi_labels.
  pm.min_print_ = true;
  pm.SetLuminosityTag(lumi_label).MakePlots(1);
}


  /*vec_nbins_cr_all= {
  {{100, 100, 100}, { , , }, { , , }, { , , }},
  {{100, 100, 100}, { , , }, { , , }, { , , }},
  {{12,   12,  20}, { , , }, { , , }, { , , }},
  {{12,   12,  20}, { , , }, { , , }, { , , }},
  {{12,   12,  20}, { , , }, { , , }, { , , }},
  {{5,     5,  10}, { , , }, { , , }, { , , }}};
  */


/*
 *
 *
 *
 *
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
        //cout << "NAME: " << control_region_name + category_specific_control_region_name << endl;
        //cout << "SELECTION: " << (control_region_selection && category_specific_control_region_selection) << endl;
        //cout << "NAME: " << control_region_name + "_NOT_" + category_specific_control_region_name << endl;
        //cout << "SELECTION: " << (control_region_selection && !category_specific_control_region_selection) << endl;

        control_region_selection = controlregions_vec[idx_cr] && lep_flavor[idx_lep];
        control_region_name =  lep_str[idx_lep] + cr_str_vec[idx_cr];

      }

    }
  }



 */
