//-----------------------------README-------------------------------//
//This script performs the backend work of producing tables. The
//script takes as input a json file (use relative path from current
//directory being used), containing information on the samples used, 
//years, and various other options.
//
//The layout of this file is:
//--Various functions not available in ZgUtilities and ZgFunctions
//--Cutflow constructor, creating the table rows in order
//--Table Generator, creates the table with various options
//--Main interprets the json and calls the other functions. Saves
//the tables to subfolder within the tables/ directory in draw_pico

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <getopt.h>
#include <regex>
#include <experimental/filesystem>

#include "TError.h" // Controls error level reporting
#include "TColor.h" // Controls error level reporting
#include "TLegend.h" // Controls error level reporting
#include "TVector2.h"

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/event_scan.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/styles.hpp"
#include "core/hist1d.hpp"
#include "core/functions.hpp"
#include "zgamma/zg_utilities.hpp"
#include "zgamma/zg_functions.hpp"
#include "json.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace ZgFunctions;
using namespace ZgUtilities;

using json = nlohmann::json;
namespace fs = std::experimental::filesystem;
using namespace fs;

set<string> vtos(vector<string> v){
  set<string> out(v.begin(), v.end());
  return out;
}

bool Contains(const string& text, const string& pattern){
  return text.find(pattern) != string::npos;
}

string getLuminosityString(vector<string> years) {
  float total_luminosity = 0;
  for (auto const & year : years) {
    // https://twiki.cern.ch/twiki/bin/view/CMS/PdmVDatasetsUL2016
    if (year == "2016APV") total_luminosity += 19.51;
    if (year == "2016") total_luminosity += 16.80;
    if (year == "2017") total_luminosity += 41.48;
    if (year == "2018") total_luminosity += 59.83;
    if (year == "2022") total_luminosity +=7.9804;
    if (year == "2022EE") total_luminosity +=26.6717;
    if (year == "2023") total_luminosity +=17.794;
    if (year == "2023BPIX") total_luminosity +=9.451;
 }
  string total_luminosity_string = RoundNumber(total_luminosity, 1, 1).Data();
  return total_luminosity_string;
}

const NamedFunc w_lumi("w_lumi", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleTypeString().Contains("-")) return 1.;
  else return b.w_lumi();
});

const NamedFunc fullwgt("fullwgt",[](const Baby &b) -> NamedFunc::ScalarType{
  if(b.SampleTypeString().Contains("-")) return 1.;
  return b.weight();
 });

const NamedFunc combwgt("combwgt", [](const Baby &b) -> NamedFunc::ScalarType{ //customizable weight application
  if (b.SampleTypeString().Contains("-")) return 1.;
  return b.w_lumi()*b.w_prefire()*b.w_pu()*b.w_bhig_df()*b.w_el()*b.w_mu();
});

const NamedFunc sideband("sideband",[](const Baby &b)->NamedFunc::ScalarType{
  bool out = true;
  if(b.llphoton_m()->at(0)>120.f && b.llphoton_m()->at(0)<130.f) out=false;
  return out;
});

const NamedFunc higgs_window("higgs_window",[](const Baby &b)->NamedFunc::ScalarType{
  bool out = false;
  if(b.llphoton_m()->at(0)>120.f && b.llphoton_m()->at(0)<130.f) out=true;
  return out;
});

NamedFunc cutbitcheck(int bit){
  return NamedFunc("comp",[bit](const Baby &b)->NamedFunc::ScalarType{
    return (bit & (b.zg_cutBitMap()))>0;
  });
}

NamedFunc catbitcheck(int bit){
  return NamedFunc("comp",[bit](const Baby &b)->NamedFunc::ScalarType{
    return (bit == (b.zg_categorizationBitMap()));
  });
}



void constructCutflowTable(vector<TableRow> & tablerows, NamedFunc weight, const string window, int lepton, bool isMC) {
  NamedFunc current_cut = "1";
  if(isMC == true) current_cut = "use_event";

  tablerows = vector<TableRow>{
    TableRow("Initial events", add_cut(current_cut, "1"),0,0, weight),
  };

  if(lepton == 0) { // el specific cuts
    tablerows.push_back(TableRow("$N_e \\geq 2$", add_cut(current_cut, cutbitcheck(0b010000000000)),0,0, weight));
    tablerows.push_back(TableRow("e, ee trigger", add_cut(current_cut, cutbitcheck(0b000100000000)),0,0,weight));
    tablerows.push_back(TableRow("el trigger $p_{T}$ cuts", add_cut(current_cut, cutbitcheck(0b000010000000)),0,0, weight));
  }else if (lepton == 1) { // mu specific cuts
    tablerows.push_back(TableRow("$N_{\\mu} \\geq 2 $", add_cut(current_cut, cutbitcheck(0b001000000000)),0,0, weight));
    tablerows.push_back(TableRow("$\\mu, \\mu\\mu$ trigger", add_cut(current_cut, cutbitcheck(0b000100000000)),0,0, weight));
    tablerows.push_back(TableRow("mu trigger $p_{T}$ cuts", add_cut(current_cut, cutbitcheck(0b000010000000)),0,0, weight));
  }else{ // mu + el cuts
    tablerows.push_back(TableRow("$N_e \\geq 2 || N_{\\mu} \\geq 2 $", add_cut(current_cut, cutbitcheck(0b100000000000)),0,0, weight));
    tablerows.push_back(TableRow("e, ee trigger $||$ $\\mu, \\mu\\mu$ trigger", add_cut(current_cut, cutbitcheck(0b000100000000)),0,0, weight));
    tablerows.push_back(TableRow("mu and el trigger $p_{T}$ cuts", add_cut(current_cut, cutbitcheck(0b000010000000)),0,0, weight));
  }

  tablerows.push_back(TableRow("$N_{\\gamma} \\geq 1$", add_cut(current_cut, cutbitcheck(0b000001000000)),0,0, weight));
  tablerows.push_back(TableRow("$80 \\text{ GeV} < m_{ll} < 100 \\text{ GeV}$", add_cut(current_cut, cutbitcheck(0b000000100000)),0,0, weight));
  tablerows.push_back(TableRow("$p_{T}^{\\gamma}/m_{ll\\gamma} > 15./110$", add_cut(current_cut, cutbitcheck(0b000000010000)),0,0, weight));
  tablerows.push_back(TableRow("$m_{ll\\gamma}+m_{ll} > 185 \\text{ GeV}$", add_cut(current_cut, cutbitcheck(0b000000001000)),0,0, weight)); 
  tablerows.push_back(TableRow("Event filters", add_cut(current_cut, cutbitcheck(0b000000000010)),0,0, weight));
  if(window=="higgs") tablerows.push_back(TableRow("$120 < m_{ll\\gamma} < 130 \\text{ GeV}$", add_cut(current_cut, higgs_window),0,0, weight));
  if(window=="sideband") tablerows.push_back(TableRow("$m_{ll\\gamma}<=120 \\text{ GeV} || 130\\text{ GeV} <=m_{ll\\gamma}$", add_cut(current_cut, sideband),0,0, weight));

  tablerows.push_back(TableRow("ttH lep.", current_cut && catbitcheck(0b00100001),0,0, weight));
  tablerows.push_back(TableRow("ttH had.", current_cut && catbitcheck(0b00001001),0,0, weight));
  tablerows.push_back(TableRow("ZH $p_{T}^{miss}$", current_cut && catbitcheck(0b00000101),0,0, weight));
  tablerows.push_back(TableRow("WH/ZH $N_{l}\\geq 3$", current_cut && catbitcheck(0b00010001),0,0, weight));
  tablerows.push_back(TableRow("VBF", current_cut && catbitcheck(0b01000000),0,0, weight));
  tablerows.push_back(TableRow("ggF", current_cut && catbitcheck(0b10000000),0,0, weight));
}



void generateTable(PlotMaker& pm, Palette colors, vector<string> in_years, vector<string> in_files, const string now_str, 
  const string prod_name, const string path_to_production, const string label, const string table_luminosity, const string skim_t, 
  const string window, const string sample_t, NamedFunc weight, TString tablename,
  const float lumi, unsigned int lepton, int precision_level, bool print_unc = false, bool do_eff = false, bool rel_prev_cut = false){//Performs the details of cutflow generation

  time_t begtime, endtime;
  time(&begtime);
  vector<shared_ptr<Process> > procs_dat;
  set<string> pathNames;
  string title, color, loc;
  bool isMC = false;

  if(sample_t == "bkg"){
    color = "dy"; title = "Background"; isMC = true; loc = "mc";
    pathNames = attach_folder(path_to_production, vtos(in_years), loc + "/" + skim_t, vtos(in_files));
    procs_dat.push_back(Process::MakeShared<Baby_pico>(title, Process::Type::background, colors(color), pathNames, "1"));
  } else if(sample_t == "sig"){
    color = "other"; title = "Signal"; isMC = true; loc = "mc";
    pathNames = attach_folder(path_to_production, vtos(in_years), loc + "/"  + skim_t, vtos(in_files));
    procs_dat.push_back(Process::MakeShared<Baby_pico>(title, Process::Type::signal, colors(color), pathNames, "1"));
  } else if(sample_t == "data"){
    color = "data"; title = "Data"; loc = "data";
    pathNames = attach_folder(path_to_production, vtos(in_years), loc + "/" + skim_t, vtos(in_files));
    procs_dat.push_back(Process::MakeShared<Baby_pico>(title, Process::Type::data, colors(color), pathNames, "1"));
  }

  vector<TableRow> tablerows;
  constructCutflowTable(tablerows, weight, window, lepton, isMC);
  pm.Push<Table>(tablename.Data(), tablerows, procs_dat, 0,1,0,1,do_eff,print_unc,rel_prev_cut).LuminosityTag(table_luminosity).Precision(precision_level);
  pm.min_print_ = true;
  pm.MakePlots(lumi, (prod_name+"_"+now_str+"/"+label));
  pm.Clear();

  time(&endtime);
  cout<<endl<<"Making cutflow took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}




int main(int argc, char* argv[]){//This script creates cutflows of the Run 2 and Run 3 data for validation purposes using json inputs
  if(argc != 2){
    cout<<"Input should be the path to one json file containing table options"<<endl;
    return 1;
  }
  string infile = argv[1];
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  Palette colors("txt/colors.txt","default");
  PlotMaker pm;
  TString tablename;
  
  //Load instructions from json
  ifstream in_json_data(infile);
  json inputs = json::parse(in_json_data);
  cout<<"Creating cutflows using "<<infile<<" instructions :)"<<endl;

  // Initializing paths and selections
  float lumi = 1;
  string prod_name, label_s, output_folder, skim_t, sample_t, weight, path_base, tag, nanov, path_to_production;
  string window = "full";
  vector<string> years, filenames, label, lepton_t, in_years, in_files;
  int n_year_tables, n_sample_tables;
  int precision_level = 0;
  int lep = 2;

  bool combine_years = false;
  bool combine_samples = false;
  bool print_unc = false;
  bool do_eff = false;
  bool rel_prev_cut = false;

  NamedFunc weight_i = "1";

  time_t now;
  time(&now);
  stringstream now_s;
  now_s << now;
  string now_str = now_s.str();

  auto jlist = inputs["set_of_cutflows"];
  prod_name = inputs["production_name"];
  path_base = inputs["path_to_production_base"];
  tag = inputs["production_tag"];
  cout<<"Producing "<<inputs["set_of_cutflows"].size()<<" sets of cutflows:"<<endl;

  for(unsigned int i = 0;i<inputs["set_of_cutflows"].size();i++){//Loop over cutflow sets
    auto set_i = jlist[i];
    label_s =                       set_i["short_name"];
    nanov =                         set_i["NanoAOD_tag"];
    output_folder =                 "tables/"+prod_name+"_"+now_s.str()+"/"+label_s;
    fs::create_directories(output_folder);
    years =                         set_i["years"].get<vector<string>>();
    combine_years =                 set_i["bool_combine_years"];
    combine_samples =               set_i["bool_combine_samples"];
    skim_t =                        set_i["skim_type"];
    sample_t =                      set_i["sample_type"];
    filenames =                     set_i["samples"].get<vector<string>>();
    label =                         set_i["sample_names"].get<vector<string>>();
    weight =                        set_i["weight"];
    window =                        set_i["mlly_window"];
    lepton_t =                      set_i["lepton_type"].get<vector<string>>();
    print_unc =                     set_i["bool_print_unc"];
    do_eff =                        set_i["bool_do_eff"];
    rel_prev_cut =                  set_i["bool_rel_prev"];
    precision_level =               set_i["decimal_precision"];

    if(combine_years == true) n_year_tables = 1;
    else n_year_tables = years.size();
    if(combine_samples == true) n_sample_tables = 1;
    else n_sample_tables = filenames.size();

    if(weight == "1") weight_i = "1";//defined weighting options
    else if(weight == "lumi") weight_i = w_years*w_lumi;
    else if(weight == "custom") weight_i = combwgt*w_years;
    else if(weight == "weight") weight_i = fullwgt*w_years;

    for(int yr = 0;yr<n_year_tables;yr++){//loop over years (if not combining)
      for(int samp = 0;samp<n_sample_tables;samp++){//loop over samples (if not combining)
        string table_luminosity;
        if(combine_years==true){
          path_to_production = path_base + "/" + nanov + "/" + tag;
          table_luminosity = getLuminosityString(years);
          in_years = years;
        } else{
          path_to_production = path_base + "/" + nanov + "/" + tag;
          table_luminosity = getLuminosityString({years[yr]});
          in_years = {years[yr]};
        }
        if(combine_samples==true) in_files = filenames;
        else in_files = {filenames[samp]};

        for(unsigned int l = 0;l<lepton_t.size();l++){//loop over leptons
          if (lepton_t[l]=="electron") lep = 0;
          else if(lepton_t[l]=="muon") lep = 1;
          else if(lepton_t[l]=="lepton") lep = 2;
          tablename = "FixName:cutflow_table_";
          if(combine_samples==false) tablename += label[samp];
          else tablename += label_s;
          tablename += lepton_t[l];
          if(combine_years==false) tablename += in_years[0];//If doing combined R2 or R3 make sure to mention in the short_name
          generateTable(pm, colors, in_years, in_files, now_str, prod_name, path_to_production, 
            label_s, table_luminosity, skim_t, window, sample_t, weight_i, tablename, lumi, 
            lep, precision_level, print_unc, do_eff, rel_prev_cut);
        }
      }//sample loop
    }//year loop
  }//cutflow production loop

  fs::copy(__FILE__, "tables/"+prod_name+"_"+now_str+"/cutflow_BE_strict_"+now_str +".cxx", fs::copy_options::update_existing);
  fs::copy(infile, "tables/"+prod_name+"_"+now_str+"/cutflow_inputs_strict_"+now_str+".json", fs::copy_options::update_existing);
  //Saves a copy of FE and BE to the output directory. Append date and time to the name of directory


}
