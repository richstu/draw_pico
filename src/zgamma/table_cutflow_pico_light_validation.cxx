#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <getopt.h>
#include <regex>

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

using namespace std;
using namespace PlotOptTypes;


namespace{
  const string path_to_production_base = "/net/cms11/cms11r0/pico/"; //Some things may need to be changed in main() if path structure changes significantly
  const string tag_a= "htozgamma_kingscanyon_v1";
  const string skim_a = "merged_*_llg";
  float lumi = 1;
  string tag = "";
}

bool Contains(const string& text, const string& pattern){
  return text.find(pattern) != string::npos;
}

NamedFunc add_cut(NamedFunc & current_cut, NamedFunc additional_cut) {
  current_cut = current_cut && additional_cut;
  return current_cut;
}
/*
string getLuminosityString(string const & years_string_a) {
  float total_luminosity = 0;
  for (auto const & year : years_string_a) {
    // https://twiki.cern.ch/twiki/bin/view/CMS/PdmVDatasetsUL2016
    if (strcmp(&year, "2016APV") == 0) total_luminosity += 19.51;
    if (strcmp(&year, "2016") == 0) total_luminosity += 16.80;
    if (strcmp(&year, "2017") == 0) total_luminosity += 41.48;
    if (strcmp(&year, "2018") == 0) total_luminosity += 59.83;
    if (strcmp(&year, "2022") == 0) total_luminosity +=8.17;
    if (strcmp(&year, "2022EE") == 0) total_luminosity +=27.01;
    if (strcmp(&year, "2023") == 0) total_luminosity +=17.61;
  }
  
  string total_luminosity_string = RoundNumber(total_luminosity, 1, 1).Data();
  std::cout<<total_luminosity_string<<endl;
  return total_luminosity_string;
}
*/

string getLuminosityString(set<string> years) {
  float total_luminosity = 0;
  for (auto const & year : years) {
    // https://twiki.cern.ch/twiki/bin/view/CMS/PdmVDatasetsUL2016
    if (year == "2016APV") total_luminosity += 19.51;
    if (year == "2016") total_luminosity += 16.80;
    if (year == "2017") total_luminosity += 41.48;
    if (year == "2018") total_luminosity += 59.83;
    if (year == "2022") total_luminosity += 8.17;
    if (year == "2022EE") total_luminosity += 27.01;
    if (year == "2023") total_luminosity +=17.61;
    if (year == "2023BPIX") total_luminosity += 9.53;
  }
  string total_luminosity_string = RoundNumber(total_luminosity, 1, 1).Data();
  return total_luminosity_string;
}


// Luminosities
const NamedFunc w_years("w_years", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleTypeString().Contains("-")) return 1.;
  if (b.SampleTypeString()=="2016APV") return 19.51;
  else if (b.SampleTypeString()=="2016") return 16.80;
  else if (b.SampleTypeString()=="2017") return 41.48;
  else if (b.SampleTypeString()=="2018") return 59.83;
  else if (b.SampleTypeString()=="2022") return 8.17;
  else if (b.SampleTypeString()=="2022EE") return 27.01;
  else if (b.SampleTypeString()=="2023") return 17.61; 
  else if (b.SampleTypeString()=="2023BPIX") return 9.53;
  else return 0.;
});

const NamedFunc w_lumi("w_lumi", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleTypeString().Contains("-")) return 1.;
  else return b.w_lumi();
});

const NamedFunc wgt("wgt",[](const Baby &b) -> NamedFunc::ScalarType{
  if(b.SampleTypeString().Contains("-")) {
    return 1;
  }
  return b.weight();
 });

const NamedFunc signal_lead_muon_pt("signal_lead_muon_pt",[](const Baby &b) -> NamedFunc::ScalarType{
  for (unsigned iPart = 0; iPart<b.mu_pt()->size(); iPart++) {
    if (b.mu_sig()->at(iPart)) return b.mu_pt()->at(iPart);
  }
  return -999;
});
const NamedFunc signal_lead_electron_pt("signal_lead_electron_pt",[](const Baby &b) -> NamedFunc::ScalarType{
  for (unsigned iPart = 0; iPart<b.el_pt()->size(); iPart++) {
    if (b.el_sig()->at(iPart)) return b.el_pt()->at(iPart);
  }
  return -999;
});

const NamedFunc signal_sublead_muon_pt("signal_sublead_muon_pt",[](const Baby &b) -> NamedFunc::ScalarType{
  bool sublead = false;
  for (unsigned iPart = 0; iPart<b.mu_pt()->size(); iPart++) {
    if (b.mu_sig()->at(iPart)) {
      if (sublead == false) sublead = true;
      else return b.mu_pt()->at(iPart);
    }
  }
  return -999;
});
const NamedFunc signal_sublead_electron_pt("signal_sublead_electron_pt",[](const Baby &b) -> NamedFunc::ScalarType{
  bool sublead = false;
  for (unsigned iPart = 0; iPart<b.el_pt()->size(); iPart++) {
    if (b.el_sig()->at(iPart)) {
      if (sublead == false) sublead = true;
      else return b.el_pt()->at(iPart);
    }
  }
  return -999;
});


void constructCutflowTable(vector<TableRow> & tablerows, NamedFunc weight, int electron_or_muon, bool isMC = false) {
    NamedFunc current_cut = "1";
    if (isMC == true){
      current_cut = "use_event";
    }
    NamedFunc el_trigger = "ll_lepid[0]==11 && (HLT_Ele27_WPTight_Gsf || HLT_Ele30_WPTight_Gsf || HLT_Ele35_WPTight_Gsf || HLT_Ele32_WPTight_Gsf || HLT_Ele32_WPTight_Gsf_L1DoubleEG || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL)"; //run 2 and run 3 e, ee triggers. Add in Run 3
    NamedFunc mu_trigger = "ll_lepid[0]==13 && (HLT_IsoMu24 || HLT_IsoTkMu24 || HLT_IsoMu27 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8)"; //run 2 and run 3 mu, mumu triggers. Add in Run 3
    if (electron_or_muon == 0) { // el
      tablerows = vector<TableRow>{
        TableRow("No selection", add_cut(current_cut,"1"),0,0, weight),
        TableRow("$N_e \\geq 2$", add_cut(current_cut, "nel>=2 && ll_lepid[0]==11"),0,0, weight),
        TableRow("e, ee trigger", add_cut(current_cut, el_trigger),0,0, weight),
        TableRow("$p_{T}^{\\text{lead }e} \\geq 25\\text{ GeV}$", add_cut(current_cut, signal_lead_electron_pt > 25),0,0, weight),
        TableRow("$p_{T}^{\\text{sublead }e} \\geq 15\\text{ GeV}$", add_cut(current_cut, signal_sublead_electron_pt > 15),0,1, weight),
        TableRow("$N_{\\gamma} \\geq 1$", add_cut(current_cut, "nphoton>=1"),0,0, weight),
        TableRow("$Charge_{ll} = 0$", add_cut(current_cut,"ll_charge[0] == 0"),0,0, weight),
        TableRow("$80 \\text{ GeV} \\leq m_{ee} \\leq 100 \\text{ GeV}$", add_cut(current_cut, "ll_m[0] > 80 && ll_m[0]<100"),0,0, weight),
        TableRow("$p_{T}^{\\gamma}/m_{ee\\gamma} > 15./110$", add_cut(current_cut, "(photon_pt[0]/llphoton_m[0])>(15./110)"),0,0, weight),
        TableRow("$m_{ee} + m_{ee\\gamma} > 185 \\text{ GeV}$", add_cut(current_cut, "(llphoton_m[0] + ll_m[0]) > 185"),0,0, weight),
        TableRow("$100 < m_{ee\\gamma} < 180 \\text{ GeV}$", add_cut(current_cut, "llphoton_m[0] > 100 && llphoton_m[0] < 180"),0,0, weight),
      };
    } else if (electron_or_muon == 1) { // mu
      tablerows = vector<TableRow>{
        TableRow("No selection", add_cut(current_cut,"1"),0,0, weight),
        TableRow("$N_{\\mu} \\geq 2 $", add_cut(current_cut, "nmu>=2 && ll_lepid[0]==13"),0,0, weight),
        TableRow("$\\mu, \\mu\\mu$ trigger", add_cut(current_cut, mu_trigger),0,0, weight),
        TableRow("$p_{T}^{\\text{lead }\\mu} \\geq 20\\text{ GeV}$", add_cut(current_cut, signal_lead_muon_pt > 20),0,0, weight),
        TableRow("$p_{T}^{\\text{sublead }\\mu} \\geq 10\\text{ GeV}$", add_cut(current_cut, signal_sublead_muon_pt > 10),0,1, weight),
        TableRow("$N_{\\gamma} \\geq 1$", add_cut(current_cut, "nphoton>=1"),0,0, weight),
        TableRow("$Charge_{ll} = 0$", add_cut(current_cut,"ll_charge[0] == 0"),0,0, weight),
        TableRow("$80 \\text{ GeV} \\leq m_{\\mu\\mu} \\leq 100 \\text{ GeV}$", add_cut(current_cut, "ll_m[0] > 80 && ll_m[0]<100"),0,0, weight),
        TableRow("$p_{T}^{\\gamma}/m_{\\mu\\mu\\gamma} > 15./110$", add_cut(current_cut, "(photon_pt[0]/llphoton_m[0])>(15./110)"),0,0, weight),
        TableRow("$m_{\\mu\\mu\\gamma}+m_{\\mu\\mu} > 185 \\text{ GeV}$", add_cut(current_cut, "(llphoton_m[0]+ll_m[0]) > 185"),0,0, weight),
        TableRow("$100 < m_{\\mu\\mu\\gamma} < 180 \\text{ GeV}$", add_cut(current_cut, "llphoton_m[0] > 100 && llphoton_m[0] < 180"),0,0, weight),
      };
    } else { // mu + el
      tablerows = vector<TableRow>{
        TableRow("No selection", add_cut(current_cut,"1"),0,0, weight),
        TableRow("$N_e \\geq 2 || N_{\\mu} \\geq 2 $", add_cut(current_cut, "(nel>=2 && ll_lepid[0]==11) || (nmu>=2 && ll_lepid[0]==13)"),0,0, weight),
        TableRow("e, ee trigger $||$ $\\mu, \\mu\\mu$ trigger", add_cut(current_cut, el_trigger || mu_trigger),0,0, weight),
        TableRow("$p_{T}^{\\text{lead }l} \\geq 25 (20)\\text{ GeV}$", add_cut(current_cut, ("ll_lepid[0]==13"&&signal_lead_muon_pt > 20) || ("ll_lepid[0]==11"&&signal_lead_electron_pt > 25)),0,0, weight),
        TableRow("$p_{T}^{\\text{sublead }l} \\geq 15 (10)\\text{ GeV} $", add_cut(current_cut, ("ll_lepid[0]==13"&&signal_sublead_muon_pt > 10) || ("ll_lepid[0]==11"&&signal_sublead_electron_pt > 15)),0,1, weight),
        TableRow("$N_{\\gamma} \\geq 1$", add_cut(current_cut, "nphoton>=1"),0,0, weight),
        TableRow("$Charge_{ll} = 0$", add_cut(current_cut,"ll_charge[0] == 0"),0,0, weight),
        TableRow("$80 \\text{ GeV} \\leq m_{ll} \\leq 100 \\text{ GeV}$", add_cut(current_cut, "ll_m[0] > 80 && ll_m[0]<100"),0,0, weight),
        TableRow("$p_{T}^{\\gamma}/m_{ll\\gamma} > 15./110$", add_cut(current_cut, "(photon_pt[0]/llphoton_m[0])>(15./110)"),0,0, weight),
        TableRow("$m_{ll\\gamma}+m_{ll} > 185 \\text{ GeV}$", add_cut(current_cut, "(llphoton_m[0]+ll_m[0]) > 185"),0,0, weight),
        TableRow("$100 < m_{ll\\gamma} < 180 \\text{ GeV}$", add_cut(current_cut, "llphoton_m[0] > 100 && llphoton_m[0] < 180"),0,0, weight),
      };
    }
}

void generateTable(int run_no, PlotMaker& pm, Palette colors, set<string> year_t, const string skim_a, NamedFunc weight_a, TString tablename_a, 
      const string total_luminosity_string_a, const float lumi, unsigned int lepton, bool isMC = false, bool isSig = false){
  //Performs the details of cutflow generation
  bool print_uncertainty = false;
  time_t begtime, endtime;
  time(&begtime);
  vector<shared_ptr<Process> > procs_dat_a;
  set<string> filenames;
  set<string> pathNames;
  string type;
  string label;
  string color;


  if(isMC == true && isSig == false){
    filenames = {"*DYJetsToLL*amcatnloFXFX*.root", "*ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX*.root"}; //Check on ZGToLLG
    if(run_no == 3){ //Change how this check is done
      filenames = {"*DYto2L*amcatnloFXFX*.root", "*DYGto2LG-1Jets_MLL-50*.root"}; //Double check on these
    }
    type = "mc/"; label = "DYJets and SMZy"; color = "t1tttt"; print_uncertainty = true;
    pathNames = attach_folder(path_to_production_base, year_t, type + skim_a, filenames); 
    procs_dat_a.push_back(Process::MakeShared<Baby_pico>(label, Process::Type::signal, colors(color), pathNames, "1")); //In if statement because I could not find way to dynamically choose the process type.
  } else if(isMC == true && isSig == true){
    //filenames = {"*GluGlu*HToZG*ZToLL_M-125_*.root","VBF*HToZG*ZToLL_M-125_*.root","*WplusH*HToZG*ZToLL_M-125_*.root", "*WminusH*HToZG*ZToLL_M-125_*.root", "*tt*HToZG*ZToLL_M-125_*.root", "*ZH*HToZG_ZToAll_M-125_*.root"};
    filenames = {"*GluGlu*HToZG*ZToLL_M-125_*.root","VBF*HToZG*ZToLL_M-125_*.root", "*tt*HToZG*ZToLL_M-125_*.root", "*ZH*HToZG_ZToAll_M-125_*.root"};//Add back in WplusH and WminusH later
    type = "mc/"; label = "Signal MC"; color = "dy"; print_uncertainty = true;
    pathNames = attach_folder(path_to_production_base, year_t, type + skim_a, filenames); 
    procs_dat_a.push_back(Process::MakeShared<Baby_pico>(label, Process::Type::background, colors(color), pathNames, "1"));
  } else if(isMC == false){
    filenames = {"*.root"};
    type = "data/"; label = "Data"; color = "data"; print_uncertainty = false;
    pathNames = attach_folder(path_to_production_base, year_t, type + skim_a, filenames); 
    procs_dat_a.push_back(Process::MakeShared<Baby_pico>(label, Process::Type::data, colors(color), pathNames, "1"));
  }

  vector<TableRow> tablerows_a;
  constructCutflowTable(tablerows_a, weight_a, lepton, isMC);
  pm.Push<Table>(tablename_a.Data(), tablerows_a, procs_dat_a, 0, 1, 0, 1, 0, print_uncertainty).LuminosityTag(total_luminosity_string_a).Precision(3);
  pm.min_print_ = true;
  pm.MakePlots(lumi);
  pm.Clear();
  time(&endtime);
  cout<<endl<<"Making cutflow took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}


int main(){
  //Minimal validation script: Does Z->ll only for run 2, run 3, data, background, and signal: 6 cutflow tables
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches


  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////// Initializing paths and selections ////////////////////////////


  string lep_tag;

  NamedFunc weight_a = "1";
  Palette colors("txt/colors.txt", "default");
  PlotMaker pm; //Define just once

  set<string> run2_year_paths = {"NanoAODv9/" + tag_a + "/2016APV", "NanoAODv9/" + tag_a + "/2016", "NanoAODv9/" + tag_a + "/2017", "NanoAODv9/" + tag_a + "/2018"};
  set<string> run3_year_paths = {"NanoAODv11/" + tag_a + "/2022", "NanoAODv11/" + tag_a + "/2022EE", "NanoAODv11p9/" + tag_a + "/2023"};
  set<string> run2_years = {"2016APV","2016","2017","2018"};
  //set<string> run3_years = {"2022", "2022EE", "2023"};
  set<string> run3_years = {"2022", "2022EE"}; //2023 removed temporarily
  /////////////////////////////////////// Data and MC Validation ////////////////////////////////////
  //////////////////////////////////// Loop over years and leptons //////////////////////////////////

  lep_tag = "_lepton";
  TString tablename_a = "FixName:validation_table_"+tag_a+lep_tag+"_"+"run2";
  string total_luminosity_string_a = getLuminosityString(run2_years);

  //    Generate run2 data cutflow table
  //--------------------------------------------------
  tablename_a = "FixName:validation_table_"+tag_a+lep_tag+"_"+"run2"+"_data";
  weight_a = "1";
  generateTable(2, pm, colors, run2_year_paths, skim_a, weight_a, tablename_a, total_luminosity_string_a, lumi, 2, false, false);

  //    Generate run2 background MC (w/ SF) cutflow table
  //--------------------------------------------------
  tablename_a = "FixName:validation_table_"+tag_a+lep_tag+"_"+"run2"+"_backmc_sf";
  weight_a = wgt*w_years;
  generateTable(2, pm, colors, run2_year_paths, skim_a, weight_a, tablename_a, total_luminosity_string_a, lumi, 2, true, false);

  //    Generate run2 signal MC (w/ SF) cutflow table
  //--------------------------------------------------
  tablename_a = "FixName:validation_table_"+tag_a+lep_tag+"_"+"run2"+"_sigmc";
  weight_a = wgt*w_years;
  generateTable(2, pm, colors, run2_year_paths, skim_a, weight_a, tablename_a, total_luminosity_string_a, lumi, 2, true, true);



  total_luminosity_string_a = getLuminosityString(run3_years);
  //    Generate run3 data cutflow table
  //--------------------------------------------------
  tablename_a = "FixName:validation_table_"+tag_a+lep_tag+"_"+"run3"+"_data";
  weight_a = "1";
  generateTable(3, pm, colors, run3_year_paths, skim_a, weight_a, tablename_a, total_luminosity_string_a, lumi, 2, false, false);
  
  //    Generate run3 background MC (w/ SF) cutflow table
  //--------------------------------------------------
  tablename_a = "FixName:validation_table_"+tag_a+lep_tag+"_"+"run3"+"_backmc_sf";
  weight_a = wgt*w_years;
  generateTable(3, pm, colors, run3_year_paths, skim_a, weight_a, tablename_a, total_luminosity_string_a, lumi, 2, true, false);

  //    Generate run3 signal MC (w/ SF) cutflow table
  //--------------------------------------------------
  tablename_a = "FixName:validation_table_"+tag_a+lep_tag+"_"+"run3"+"_sigmc_sf";
  weight_a = wgt*w_years;
  generateTable(3, pm, colors, run3_year_paths, skim_a, weight_a, tablename_a, total_luminosity_string_a, lumi, 2, true, true);
  

}

