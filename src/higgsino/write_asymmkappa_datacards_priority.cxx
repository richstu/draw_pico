#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <sstream>
#include <ctime>
#include <algorithm>
#include <unistd.h> // getopt in Macs
#include <getopt.h>
#include <dirent.h>
#include <regex>

#include "TSystem.h"
#include "TString.h"
#include "TMath.h"
#include "TError.h" // Controls error level reporting
#include "TVector2.h"

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/table.hpp"
#include "core/palette.hpp"
#include "core/config_parser.hpp"
#include "core/functions.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/write_kappa_datacards.hpp"
#include "higgsino/hig_utilities.hpp"

using namespace std;
namespace 
{
  string outFolder = getenv("PWD");
  string mass_points_string = ""; // run over all the points, otherwise specify, e.g. "175_1,200_1"
  string mass_point_glob = "";
  string years_string = "2016,2017,2018";
  float luminosity = 1.;
  //float luminosity = 35.9;
  string dimensionFilePath = "";
  bool unblind = false;
  string tag = "resolved";
  bool do_met_average = true;
  string higgsino_model = "N1N2";
  // 1. RSR, RCR, BSR, BCR
  // 2. BSR, BCR, RSR, RCR
  // 3. RSR, BSR, RCR, BCR
  // 4. BSR, RSR, BCR, RCR
  int priority = 1;
  bool is1D = true;
}

const NamedFunc min_jet_dphi("min_jet_dphi", [](const Baby &b) -> NamedFunc::ScalarType{
  float min_dphi = 4;
  for (size_t ijet = 0; ijet < (*b.jet_phi()).size(); ++ijet) {
    float dphi = fabs(TVector2::Phi_mpi_pi((*b.jet_phi())[ijet]-b.met_phi()));
    if (dphi < min_dphi) min_dphi = dphi;
  }
  return min_dphi;
});

string regSearch(string inString, string const & regExPattern) {
    std::smatch match;
    std::regex pattern(regExPattern);
    std::regex_search (inString, match, pattern);
    return match.str();
}

string getPartialDatasetName(string filename) {
  string baseFilename = filename.substr(filename.find_last_of("/\\")+1);
  string fullDatasetName = regSearch(baseFilename, "[A-Z].*|ttHTobb.*");
  //cout<<datasetName.substr(0,datasetName.find_last_of("Tune"))<<endl;;
  string datasetName = fullDatasetName;
  if (datasetName.find("Run") != string::npos) datasetName = "DATA";
  if (datasetName.find("SMS-TChiHH") != string::npos) datasetName = "TChiHH";
  if (datasetName.find("Tune") != string::npos) datasetName = datasetName.substr(0, datasetName.find("Tune")+4);
  if (datasetName.find("_13TeV") != string::npos) datasetName = datasetName.substr(0, datasetName.find("_13TeV"));
  return datasetName;
}

bool findStringIC(const std::string & strHaystack, const std::string & strNeedle)
{
  auto it = std::search(
    strHaystack.begin(), strHaystack.end(),
    strNeedle.begin(),   strNeedle.end(),
    [](char ch1, char ch2) { return std::toupper(ch1) == std::toupper(ch2); }
  );
  return (it != strHaystack.end() );
}

pair<int, int> getSignalMassValues(string filename) {
  //pair<int, int> nlsp_lsp_mass;
  string baseFilename = filename.substr(filename.find_last_of("/\\")+1);
  string datasetName = regSearch(baseFilename, "[A-Z].*|ttHTobb.*");
  if (datasetName.find("SMS-TChiHH") == string::npos) return {-1,-1};
  string NLSPMassString = regSearch(datasetName, "mChi-.*?_").substr(5);
  string LSPMassString = regSearch(datasetName, "mLSP-.*?_").substr(5);
  //cout<<NLSPMassString<<" "<<LSPMassString<<endl;
  return {stoi(NLSPMassString),stoi(LSPMassString)};
}

bool regionCut(const Baby & b, int regionIndex) {
  bool inRegion = false;
  vector<map<Long64_t, set<tuple<string, int, int, int, int, int, int> > > > * eventNumberData = static_cast<vector<map<Long64_t, set< tuple<string, int, int, int, int, int, int> > > > *> (b.EventVetoData());
  // 0: sampleType, 1: year, 2: run, 3: lumiblock, 4: met, 5: nlsp_mass, 6: lsp_mass
  // Check event nuber
  if ((*eventNumberData)[regionIndex].count(b.event())==1) {
    for (auto const & it : (*eventNumberData)[regionIndex][b.event()]) {
      // Check year
      if (get<1>(it) != abs(b.SampleType())) continue;
      // Check run
      if (get<2>(it) != b.run()) continue;
      bool isSignal = ((*b.FileNames().begin()).find("TChiHH") != string::npos) ? true : false;
      // Check met out of 20%
      if (isSignal) {if (get<4>(it) < b.met()*0.8 || get<4>(it) > b.met()*1.2) continue;}
      else if (get<4>(it) != round(b.met())) continue;
      pair<int, int> signal_mass = getSignalMassValues(*b.FileNames().begin());
      // Check nlsp_mass
      if (signal_mass.first != get<5>(it)) continue;
      // Check lsp_mass
      int regionLsp = get<6>(it)==1? 0:get<6>(it);
      if (signal_mass.second != regionLsp) continue;
      // Check sample name
      //if (get<5>(it) != -1) cout<<*b.FileNames().begin()<<" "<<get<0>(it)<<endl;
      if ((*b.FileNames().begin()).find("TChiHH") != string::npos && get<0>(it).find("TChiHH") != string::npos)  {
        if ((*b.FileNames().begin()).find("HToBB_2D") != get<0>(it).find("2D")) continue;
      } else if (!findStringIC(*b.FileNames().begin(), get<0>(it))) continue;

      //if (b.met()<300) cout<<get<0>(it)<<" "<<get<1>(it)<<" "<<get<2>(it)<<" "<<get<3>(it)<<" "<<b.event()<<endl;
      //cout<<"sampleType: "<<get<0>(it)<<" year: "<<get<1>(it)<<" run: "<<get<2>(it)<<" lumiblock: "<<get<3>(it)<<" met: "<<b.event()<<get<4>(it)<<" nlsp_mass: "<<b.event()<<get<5>(it)<<" lsp_mass: "<<b.event()<<get<6>(it)<<endl;

      inRegion = true;
      break;
    }
  }
  return inRegion;
}

const NamedFunc boostSignalRegion("boostSignalRegion",[](const Baby &b) -> NamedFunc::ScalarType{
  return regionCut(b, 0);
});
const NamedFunc boostControlRegion("boostControlRegion",[](const Baby &b) -> NamedFunc::ScalarType{
  return regionCut(b, 1);
});

void addEventNumberData(set<string> eventNumberFilenames, vector<map<Long64_t, set<tuple<string, int, int, int, int, int, int> > > > * eventNumberData) {
  eventNumberData->push_back(map<Long64_t, set<tuple<string, int, int, int, int, int, int> > > ());
  for (string const & eventNumberFilename : eventNumberFilenames) {
    ifstream eventNumberFile(eventNumberFilename);
    string line;
    while (getline(eventNumberFile, line)) {
      if(line[0]=='#') continue;
      if(line.find("tree") != string::npos) continue;
      vector<string> lineSplit;
      HigUtilities::stringToVectorString(line, lineSplit, ",");
      string const & sampleType = lineSplit[0];
      int year = stoi(lineSplit[1]);
      int run = stoi(lineSplit[2]);
      int lumiblock = stoi(lineSplit[3]);
      Long64_t event = stoll(lineSplit[4]);
      int met = stoi(lineSplit[5]);
      int nlsp_mass = -1; int lsp_mass = -1;
      bool isSignal = false; if (sampleType.find("TChiHH") != string::npos) isSignal = true;
      if (isSignal) {nlsp_mass = stoi(lineSplit[6]); lsp_mass = stoi(lineSplit[7]);}
      // Insert data
      if ((*eventNumberData).back().count(event) == 0) {
        //cout<<"sampleType: "<<sampleType<<" year: "<<year<<" run: "<<run<<" lumiblock: "<<lumiblock<<" met: "<<met<<" nlsp_mass: "<<nlsp_mass<<" lsp_mass: "<<lsp_mass<<endl;
        (*eventNumberData).back()[event] = {{sampleType, year, run, lumiblock, met, nlsp_mass, lsp_mass}};
      } else {
        if ((*eventNumberData).back()[event].count({sampleType, year, run, lumiblock, met, nlsp_mass, lsp_mass}) != 0) cout<<"Duplicate: "<<sampleType<<" "<<year<<" "<<run<<" "<<lumiblock<<" "<<event<<" "<<met<<" "<<nlsp_mass<<" "<<lsp_mass<<endl;
        (*eventNumberData).back()[event].insert({sampleType, year, run, lumiblock, met, nlsp_mass, lsp_mass});
      }
    }
    eventNumberFile.close();
  }

  //// Read out map
  //int regionIndex = 1;
  //for (auto const & regionIt : (*eventNumberData)[regionIndex]) {
  //  Long64_t const & event = regionIt.first;
  //  for (auto const & it : (*eventNumberData)[regionIndex][event]) {
  //      cout<<"sampleType: "<<get<0>(it)<<" year: "<<get<1>(it)<<" run: "<<get<2>(it)<<" lumiblock: "<<get<3>(it)<<" met: "<<get<4>(it)<<" nlsp_mass: "<<get<5>(it)<<" lsp_mass: "<<get<6>(it)<<endl;
  //  }
  //}

}

int main(int argc, char *argv[])
{
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  time_t begtime, endtime;
  time(&begtime);
  HigWriteDataCards::GetOptions(argc, argv);

  string baseFolder = ""; string hostName = execute("echo $HOSTNAME");
  if(Contains(hostName, "cms") || Contains(hostName, "compute-")) baseFolder = "/net/cms29";
  //baseFolder = "";

  cout<<"INFO:: Systematics are ON. Make sure to run on appropriate babies, i.e. unskimmed or skim_sys_abcd!!"<<endl;
  gSystem->mkdir(outFolder.c_str(), kTRUE);

  set<int> years;
  HigUtilities::parseYears(years_string, years);
  float total_luminosity = 0;
  for (auto const & year : years) {
    if (year == 2016) total_luminosity += 35.9;
    if (year == 2017) total_luminosity += 41.5;
    if (year == 2018) total_luminosity += 60;
  }
  string total_luminosity_string = RoundNumber(total_luminosity, 1, 1).Data();

  map<string, string> samplePaths;
  ////samplePaths["mc_2016"] = baseFolder + "/cms29r0/pico/NanoAODv5/higgsino_eldorado/2016/mc/merged_higmc_preselect/";
  ////samplePaths["signal_2016"] = baseFolder + "/cms29r0/pico/NanoAODv5/higgsino_eldorado/2016/SMS-TChiHH_2D/merged_higmc_preselect/";
  ////samplePaths["data_2016"] = baseFolder + "/cms2r0/babymaker/babies/2017_02_14/data/merged_higdata_higloose/";
  //samplePaths["mc_2016"] = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/2016/mc/merged_higmc_preselect/";
  //samplePaths["signal_2016"] = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/2016/SMS-TChiHH_2D/merged_higmc_preselect/";
  //samplePaths["mc_2017"] = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldtv3/2017/mc/merged_higmc_preselect/";
  //samplePaths["signal_2017"] = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldtv3/2017/SMS-TChiHH_2D/merged_higmc_preselect/";
  //samplePaths["mc_2018"] = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/2018/mc/merged_higmc_preselect/";
  //samplePaths["signal_2018"] = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/2018/SMS-TChiHH_2D/merged_higmc_preselect/";

  //samplePaths["mc_2016"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/2016/mc/merged_higmc_preselect/";
  //samplePaths["signal_2016"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/2016/SMS-TChiHH_2D/merged_higmc_preselect/";
  //samplePaths["mc_2017"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/2017/mc/merged_higmc_preselect/";
  //samplePaths["signal_2017"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/2017/SMS-TChiHH_2D/merged_higmc_preselect/";
  //samplePaths["mc_2018"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/2018/mc/merged_higmc_preselect/";
  //samplePaths["signal_2018"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/2018/SMS-TChiHH_2D/merged_higmc_preselect/";

  samplePaths["mc_2016"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/2016/mc/merged_higmc_preselect/";
  samplePaths["signal_2016"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/2016/SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_preselect/";
  samplePaths["data_2016"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/2016/data/merged_higdata_preselect/";
  samplePaths["mc_2017"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/2017/mc/merged_higmc_preselect/";
  samplePaths["signal_2017"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/2017/SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_preselect/";
  samplePaths["data_2017"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/2017/data/merged_higdata_preselect/";
  samplePaths["mc_2018"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/2018/mc/merged_higmc_preselect/";
  samplePaths["signal_2018"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/2018/SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_preselect/";
  samplePaths["data_2018"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/2018/data/merged_higdata_preselect/";

  //// massPoints = { {"1000","1"} }
  vector<pair<string, string> > massPoints;
  //if (higgsino_model == "N1N2") {
  //  string signal_folder = samplePaths["signal_2016"];
  //  set<string> signal_files = Glob(signal_folder+"/*.root");
  //  string mLSP, mChi;
  //  for (string signal_file : signal_files) {
  //    HigUtilities::filenameToMassPoint(signal_file, mChi, mLSP);
  //    if (stoi(mChi)>800) continue; // Ingore points above 800 GeV
  //    massPoints.push_back({mChi, mLSP});
  //  }
  //}
  //else if (mass_points_string == "") HigUtilities::findMassPoints(samplePaths["signal_2016"], massPoints);
  //else HigUtilities::parseMassPoints(mass_points_string, massPoints);

  if (mass_points_string != "") {
    HigUtilities::parseMassPoints(mass_points_string, massPoints);
  } else if (mass_point_glob != "") {
    string signal_folder = samplePaths["signal_2016"];
    set<string> signal_files;
    signal_files = Glob(signal_folder+"/"+mass_point_glob);
    string mLSP, mChi;
    for (string signal_file : signal_files) {
      HigUtilities::filenameToMassPoint(signal_file, mChi, mLSP);
      if (!is1D && stoi(mChi)>800) continue; // Ingore points above 800 GeV for 2D
      massPoints.push_back({mChi, mLSP});
    }
  } else {
    string signal_folder = samplePaths["signal_2016"];
    set<string> signal_files;
    if (is1D) {
      signal_files = Glob(signal_folder+"/*mLSP-0_*.root");
    } else {
      signal_files = Glob(signal_folder+"/*.root");
    }
    string mLSP, mChi;
    for (string signal_file : signal_files) {
      HigUtilities::filenameToMassPoint(signal_file, mChi, mLSP);
      if (!is1D && stoi(mChi)>800) continue; // Ingore points above 800 GeV for 2D
      massPoints.push_back({mChi, mLSP});
    }
  }

  //NamedFunc filters = HigUtilities::pass_2016;
  //NamedFunc filters = Functions::hem_veto && "pass && met/mht<2 && met/met_calo<2";
  NamedFunc filters = Higfuncs::final_pass_filters;

  // sampleProcesses[mc, data, signal]
  map<string, vector<shared_ptr<Process> > > sampleProcesses;
  HigUtilities::setMcProcesses(years, samplePaths, filters && "stitch", sampleProcesses);
  // Higfuncs::trig_hig will only work for 2016
  //HigUtilities::setDataProcesses(years, samplePaths, filters&&Higfuncs::trig_hig>0., sampleProcesses);
  NamedFunc met_triggers = Higfuncs::met_trigger;
  if(unblind) HigUtilities::setDataProcesses(years, samplePaths, filters&&met_triggers, sampleProcesses);
  HigUtilities::setSignalProcesses(massPoints, years, samplePaths, filters, sampleProcesses);
  
  //NamedFunc weight = "weight"*Higfuncs::eff_higtrig_run2*Higfuncs::w_years;
  //NamedFunc weight = "weight"*Higfuncs::eff_higtrig_run2*Higfuncs::w_years*Functions::w_pileup;
  NamedFunc weight = Higfuncs::final_weight;
  //NamedFunc weight = "weight"*Higfuncs::eff_higtrig_run2*Higfuncs::w_years;
  //NamedFunc weight = "w_lumi*w_isr"*Higfuncs::eff_higtrig*Higfuncs::w_years;

  if (higgsino_model=="N1N2") weight *= HigUtilities::w_CNToN1N2;
  string baseline = "!low_dphi_met && nvlep==0 && ntk==0";
  string higtrim = "hig_cand_drmax[0]<=2.2 && hig_cand_dm[0] <= 40 && hig_cand_am[0]<=200 && met/mht<2 && met/met_calo<2&& weight<1.5";
  if (tag=="resolved") baseline += "&& njet>=4 && njet<=5 && nbt>=2 && "+higtrim;
  else if (tag=="boosted") {
    baseline += " && ht>600 && nfjet>1 && fjet_pt[0]>300 && fjet_pt[1]>300";
    baseline += " && fjet_msoftdrop[0]>50 && fjet_msoftdrop[0]<=250 && fjet_msoftdrop[1]>50 && fjet_msoftdrop[1]<=250";
  }

  // Set ABCDbins
  // "bkg" and "sig" should be kept for setDataCardBackground()
  map<string, string> xBins;
  if (tag=="resolved") {
    xBins["bkg"] = "nbm==2";
    xBins["sig0"] = "nbm==3&&nbl==3";
    xBins["sig1"] = "nbm>=3&&nbl>=4";
  } else if (tag=="boosted") {
    string j1bb = "((met<=300 && fjet_mva_hbb_btv[0]>0.6)||(met>300 && fjet_mva_hbb_btv[0]>0.3))";
    string j2bb = "((met<=300 && fjet_mva_hbb_btv[1]>0.6)||(met>300 && fjet_mva_hbb_btv[1]>0.3))";
    xBins["bkg"] = "(!"+j1bb+"&& !"+j2bb+")";
    xBins["sig0"] = "((!"+j1bb+"&&"+j2bb+")||("+j1bb+"&& !"+j2bb+"))";
    xBins["sig1"] = j1bb+"&&"+j2bb;
  }
  map<string, string> yBins; // Shares sideband
  if (tag=="resolved") {
    yBins["sig"] = "hig_cand_am[0]>100 && hig_cand_am[0]<=140";
    yBins["bkg"] = "!("+yBins["sig"]+")";
  } else if (tag=="boosted") {
    yBins["sig"] = "fjet_msoftdrop[0]>95 && fjet_msoftdrop[0]<=145 && fjet_msoftdrop[1]>95 && fjet_msoftdrop[1]<=145";
    yBins["bkg"] = "!("+yBins["sig"]+")";
  }
  
  map<string, vector<pair<string, string> > > dimensionBins;
  if (dimensionFilePath==""){
    if (tag=="resolved") {
      //// Old
      //dimensionBins["met"].push_back({"met0", "met>150 && met<=200"});
      //dimensionBins["met"].push_back({"met1", "met>200 && met<=300"});
      //dimensionBins["met"].push_back({"met2", "met>300 && met<=450"});
      //dimensionBins["met"].push_back({"met3", "met>450"});

      // Original
      dimensionBins["met"].push_back({"met0", "met>150 && met<=200"});
      dimensionBins["met"].push_back({"met1", "met>200 && met<=300"});
      dimensionBins["met"].push_back({"met2", "met>300 && met<=400"});
      dimensionBins["met"].push_back({"met3", "met>400"});
      dimensionBins["drmax"].push_back({"drmax0", "hig_cand_drmax[0]<=1.1"});
      dimensionBins["drmax"].push_back({"drmax1", "hig_cand_drmax[0]>1.1"});

      //// Slightly lower met binning
      //dimensionBins["met"].push_back({"met0", "met>150 && met<=200"});
      //dimensionBins["met"].push_back({"met1", "met>200 && met<=250"});
      //dimensionBins["met"].push_back({"met2", "met>250 && met<=325"});
      //dimensionBins["met"].push_back({"met3", "met>325"});
      //dimensionBins["drmax"].push_back({"drmax0", "hig_cand_drmax[0]<=1.1"});
      //dimensionBins["drmax"].push_back({"drmax1", "hig_cand_drmax[0]>1.1"});

      //// 3 met bin
      //dimensionBins["met"].push_back({"met0", "met>150 && met<=200"});
      //dimensionBins["met"].push_back({"met1", "met>200 && met<=300"});
      //dimensionBins["met"].push_back({"met2", "met>300"});
      //dimensionBins["drmax"].push_back({"drmax0", "hig_cand_drmax[0]<=1.1"});
      //dimensionBins["drmax"].push_back({"drmax1", "hig_cand_drmax[0]>1.1"});

      //// mix binning
      //dimensionBins["mix"].push_back({"met0_drmax1", "met>150 && met<=200 &&hig_cand_drmax[0]>1.1"});
      //dimensionBins["mix"].push_back({"met1_drmax1", "met>200 && met<=300 &&hig_cand_drmax[0]>1.1"});
      //dimensionBins["mix"].push_back({"met2_drmax1", "met>300 && met<=400 &&hig_cand_drmax[0]>1.1"});
      //dimensionBins["mix"].push_back({"met3_drmax1", "met>400 && hig_cand_drmax[0]>1.1"});
      //dimensionBins["mix"].push_back({"met0_drmax0", "met>150 && met<=200 &&hig_cand_drmax[0]<=1.1"});
      //dimensionBins["mix"].push_back({"met1_drmax0", "met>200 && met<=300 &&hig_cand_drmax[0]<=1.1"});
      //dimensionBins["mix"].push_back({"met2_drmax0", "met>300 &&hig_cand_drmax[0]<=1.1"});

    } else {
      dimensionBins["met"].push_back({"met0", "met>150 && met<=200"});
      dimensionBins["met"].push_back({"met1", "met>200 && met<=300"});
      dimensionBins["met"].push_back({"met2", "met>300 && met<=500"});
      dimensionBins["met"].push_back({"met3", "met>500 && met<=700"});
      dimensionBins["met"].push_back({"met4", "met>700"});
    }
  } else {
    HigWriteDataCards::readDimensionFile(dimensionFilePath, dimensionBins);
    for (auto & dimension : dimensionBins){
      cout<<"dimension bins"<<endl;
      cout<<dimension.first<<endl;
      for (auto & entry : dimension.second){
        cout<<" "<<entry.first<<" "<<entry.second<<endl;
      }
    }
  }

  // sampleBins = { {label, cut} }
  vector<pair<string, string> > sampleBins;
  // 1. RSR, RCR, BSR, BCR
  // 2. BSR, BCR, RSR, RCR
  // 3. RSR, BSR, RCR, BCR
  // 4. BSR, RSR, BCR, RCR
  HigUtilities::setABCDBinsPriority(xBins, yBins, dimensionBins, sampleBins, priority);
  //sampleBins = {{"test","1"}};

  // Set systematics somehow..
  // controlSystematics[xsig0_ybkg_met0_drmax0]["ttbar"] = systematic value
  map<string, map<string, float> > controlSystematics;
  HigWriteDataCards::setControlSystematics(controlSystematics);
  // Consider weight for nonHH somehow..

  // cuts[mc,data,signal] = RowInformation(labels, tableRows, yields)
  map<string, HigUtilities::RowInformation > cutTable;
  if(unblind) HigUtilities::addBinCuts(sampleBins, baseline, weight, "data", cutTable["data"]);
  HigUtilities::addBinCuts(sampleBins, baseline, weight, "signal", cutTable["signal"]);
  HigUtilities::addBinCuts(sampleBins, baseline, weight, "signalGenMet", HigUtilities::nom2genmet, cutTable["signal"]);
  HigUtilities::addBinCuts(sampleBins, baseline, weight, "mc", cutTable["mc"]);

  PlotMaker pm;
  set<string> boostedSignalRegion = {
                                     "BoostedEvents/boostedEvts_MC2016bkg_SR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2017bkg_SR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2018bkg_SR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2016sig1D_SR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2017sig1D_SR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2018sig1D_SR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2016sig2D_SR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2017sig2D_SR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2018sig2D_SR_noVeto.txt",
                                     };
  set<string> boostedControlRegion = {
                                     "BoostedEvents/boostedEvts_MC2016bkg_CR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2017bkg_CR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2018bkg_CR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2016sig1D_CR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2017sig1D_CR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2018sig1D_CR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2016sig2D_CR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2017sig2D_CR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2018sig2D_CR_noVeto.txt",
                                     };
  //set<string> resolvedSignalRegion = {
  //                                    "eventlist/processed_resolved_list_SR_SCAN_MC.txt",
  //                                    "eventlist/processed_resolved_list_SR_SCAN_TChiHH1D.txt",
  //                                    //"eventlist/processed_resolved_list_SR_SCAN_TChiHH2D.txt",
  //                                   };
  //set<string> resolvedControlRegion = {
  //                                    "eventlist/processed_resolved_list_CR_SCAN_MC.txt",
  //                                    "eventlist/processed_resolved_list_CR_SCAN_TChiHH1D.txt",
  //                                    //"eventlist/processed_resolved_list_CR_SCAN_TChiHH2D.txt",
  //                                    };
  //eventNumberData[0][event] = {[sampleType, year, run, lumiblock, met, NLSP, LSP]} // Boosted SR
  //eventNumberData[1][event] = {[sampleType, year, run, lumiblock, met, NSLP, LSP]} // Boosted CR
  //eventNumberData[2][event] = {[sampleType, year, run, lumiblock, met, NLSP, LSP]} // Resolved SR
  //eventNumberData[3][event] = {[sampleType, year, run, lumiblock, met, NLSP, LSP]} // Resolved CR
  vector<map<Long64_t, set<tuple<string, int, int, int, int, int, int> > > > * eventNumberData = new vector<map<Long64_t, set<tuple<string, int, int, int, int, int, int> > > >();
  addEventNumberData(boostedSignalRegion, eventNumberData);
  addEventNumberData(boostedControlRegion, eventNumberData);
  //addEventNumberData(resolvedSignalRegion, eventNumberData);
  //addEventNumberData(resolvedControlRegion, eventNumberData);
  pm.SetEventVetoData(static_cast<void * >(eventNumberData));
  // Luminosity used for labeling for table
  // Luminosity used for scaling for hist1d
  bool verbose = false;
  HigUtilities::makePlots(cutTable, sampleProcesses, luminosity, pm, verbose);

  // mYields[process_tag_sampleBinLabel] = GammaParams, TableRow
  map<string, pair<GammaParams, TableRow> > mYields;
  if(unblind) HigUtilities::fillDataYields(pm, cutTable["data"], mYields, true);
  // Luminosity used for scaling
  HigUtilities::fillMcYields(pm, luminosity, cutTable["mc"], mYields, true);
  // Luminosity used for scaling
  HigUtilities::fillSignalYieldsProcesses(pm, luminosity, sampleProcesses["signal"], cutTable["signal"], mYields, false);
  HigUtilities::fillAverageGenMetYields(sampleProcesses["signal"], sampleBins, "signal", "signalGenMet", "signalAverageGenMet", mYields);

  // Calculate kappas from mc
  // kappas[xsigX_ysigY_metA_drmaxB]
  map<string, double> kappas;
  // kappa_uncertainties[metX_drmaxY]
  map<string, pair<double, double> > kappa_uncertainties;
  HigWriteDataCards::calculateKappas(mYields, dimensionBins, kappas, kappa_uncertainties);
  //// For testing
  //for (auto & item : kappas) item.second = 2;
  //for (auto & item : kappa_uncertainties) item.second = 0.001;

  //cout<<"cut name: "<<cutTable["signal"].tableRows[5].cut_.Name()<<endl;

  for (auto & process : sampleProcesses["signal"]){
    cout<<"Process name: "<<process->name_<<endl;
    string model;
    int mGluino=0, mLSP=0;
    HigUtilities::getInfoFromProcessName(process->name_, model, mGluino, mLSP);
    TString outPath = outFolder+"/datacard-"+HigUtilities::setProcessNameLong(model, mGluino, mLSP)+"_"+years_string+"_priority"+to_string(priority)+"_"+tag+".txt";
    cout<<"open "<<outPath<<endl;
    ofstream cardFile(outPath);
    HigWriteDataCards::writeDataCardHeader(sampleBins,cardFile);

    vector<vector<string> > tableValues;
    if (unblind) HigWriteDataCards::setDataCardObserved(mYields, sampleBins, "data", tableValues);
    else HigWriteDataCards::setDataCardObserved(mYields, sampleBins, "mc", tableValues);
    if (do_met_average) {
      HigWriteDataCards::setDataCardSignalBackground(process->name_, "signalAverageGenMet", mYields, sampleBins, tableValues);
      HigWriteDataCards::setDataCardSignalStatistics(process->name_, "signalAverageGenMet", mYields, sampleBins, tableValues);
    } else {
      HigWriteDataCards::setDataCardSignalBackground(process->name_, "signal", mYields, sampleBins, tableValues);
      HigWriteDataCards::setDataCardSignalStatistics(process->name_, "signal", mYields, sampleBins, tableValues);      
    }
    HigWriteDataCards::setDataCardControlSystematics(controlSystematics, sampleBins, tableValues);
    HigWriteDataCards::writeTableValues(tableValues,cardFile);
    tableValues.clear();

    if (unblind) HigWriteDataCards::setDataCardBackground(mYields, sampleBins, "data", tableValues);
    else HigWriteDataCards::setDataCardBackground(mYields, sampleBins, "mc", tableValues);
    HigWriteDataCards::writeTableValues(tableValues,cardFile, true);
    tableValues.clear();
    HigWriteDataCards::setDataCardKappa(kappas, kappa_uncertainties, dimensionBins, tableValues);
    HigWriteDataCards::writeTableValues(tableValues, cardFile);
    tableValues.clear();
    //writeDataCardClosure()
    //writeDataCardUncertainties()
    //writeDataCardParam()
  }

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

namespace HigWriteDataCards{
  void calculateKappas(std::map<std::string, std::pair<GammaParams, TableRow> > & mYields, std::map<std::string, std::vector<std::pair<std::string, std::string> > > & dimensionBins, std::map<std::string, double> & kappas, std::map<std::string, pair<double, double> > & kappa_uncertainties) {
    map<string, vector<pair<string, string> > > t_dimensionBins = dimensionBins;
    // Combine dimensions to one list of dimensions.
    HigUtilities::combineDimensionBins(t_dimensionBins);

    auto dimension = t_dimensionBins.begin();
    for (unsigned iBin =0; iBin < dimension->second.size(); ++iBin) {
      for (unsigned iXbin = 0; iXbin < 2; iXbin++) {
        GammaParams xsig_ysig = mYields.at("mc_xsig"+to_string(iXbin)+"_ysig_"+(dimension->second)[iBin].first).first;
        GammaParams xsig_ybkg = mYields.at("mc_xsig"+to_string(iXbin)+"_ybkg_"+(dimension->second)[iBin].first).first;
        GammaParams xbkg_ysig = mYields.at("mc_xbkg_ysig_"+(dimension->second)[iBin].first).first;
        GammaParams xbkg_ybkg = mYields.at("mc_xbkg_ybkg_"+(dimension->second)[iBin].first).first;
        vector<vector<float> > entries = {{static_cast<float>(xbkg_ybkg.NEffective())}, {static_cast<float>(xsig_ybkg.NEffective())}, {static_cast<float>(xbkg_ysig.NEffective())}, {static_cast<float>(xsig_ysig.NEffective())}};
        vector<vector<float> > entry_weights = {{static_cast<float>(xbkg_ybkg.Weight())}, {static_cast<float>(xsig_ybkg.Weight())}, {static_cast<float>(xbkg_ysig.Weight())}, {static_cast<float>(xsig_ysig.Weight())}};
        vector<float> powers = {1, -1, -1, 1};
        float kappa_up = 0, kappa_down = 0;
        double kappa = calcKappa(entries, entry_weights, powers, kappa_down, kappa_up);
        //double kappa_uncertainty = kappa_up>kappa_down ? kappa_up:kappa_down;
        //cout<<xsig_ysig.Yield()<<" "<<xsig_ybkg.Yield()<<" "<<xbkg_ysig.Yield()<<" "<<xbkg_ybkg.Yield()<<endl;
        //cout<<"  kappa: "<<kappa<<" stat. uncertainty: "<<kappa_uncertainty<<endl;
        string kappaName = "xsig"+to_string(iXbin)+"_ysig_"+(dimension->second)[iBin].first;
        kappas[kappaName] = kappa;
        kappa_uncertainties[kappaName] = {kappa_down, kappa_up};
      }
    }
  }

  void setControlSystematics(map<string, map<string, float> > & controlSystematics) {
    controlSystematics["xsig0_ybkg_met0_drmax0"]["ttbar"] = 1.12; 
    controlSystematics["xsig1_ybkg_met0_drmax0"]["ttbar"] = 1.15;
    controlSystematics["xsig0_ybkg_met0_drmax1"]["ttbar"] = 1.02;
    controlSystematics["xsig1_ybkg_met0_drmax1"]["ttbar"] = 1.06;
    controlSystematics["xsig0_ybkg_met1_drmax0"]["ttbar"] = 1.10;
    controlSystematics["xsig1_ybkg_met1_drmax0"]["ttbar"] = 1.14;
    controlSystematics["xsig0_ybkg_met1_drmax1"]["ttbar"] = 1.02;
    controlSystematics["xsig1_ybkg_met1_drmax1"]["ttbar"] = 1.06;
    controlSystematics["xsig0_ybkg_met2_drmax0"]["ttbar"] = 1.06;
    controlSystematics["xsig1_ybkg_met2_drmax0"]["ttbar"] = 1.11;
    controlSystematics["xsig0_ybkg_met2_drmax1"]["ttbar"] = 1.01;
    controlSystematics["xsig1_ybkg_met2_drmax1"]["ttbar"] = 1.05;
    controlSystematics["xsig0_ybkg_met3_drmax0"]["ttbar"] = 1.06;
    controlSystematics["xsig1_ybkg_met3_drmax0"]["ttbar"] = 1.09;
    controlSystematics["xsig0_ybkg_met3_drmax1"]["ttbar"] = 1.01;
    controlSystematics["xsig1_ybkg_met3_drmax1"]["ttbar"] = 1.05;

    controlSystematics["xsig0_ybkg_met0_drmax0"]["vjets"] = 1.02; 
    controlSystematics["xsig1_ybkg_met0_drmax0"]["vjets"] = 1.01;
    controlSystematics["xsig0_ybkg_met0_drmax1"]["vjets"] = 1.00;
    controlSystematics["xsig1_ybkg_met0_drmax1"]["vjets"] = 1.00;
    controlSystematics["xsig0_ybkg_met1_drmax0"]["vjets"] = 1.04;
    controlSystematics["xsig1_ybkg_met1_drmax0"]["vjets"] = 1.05;
    controlSystematics["xsig0_ybkg_met1_drmax1"]["vjets"] = 1.01;
    controlSystematics["xsig1_ybkg_met1_drmax1"]["vjets"] = 1.00;
    controlSystematics["xsig0_ybkg_met2_drmax0"]["vjets"] = 1.09;
    controlSystematics["xsig1_ybkg_met2_drmax0"]["vjets"] = 1.06;
    controlSystematics["xsig0_ybkg_met2_drmax1"]["vjets"] = 1.01;
    controlSystematics["xsig1_ybkg_met2_drmax1"]["vjets"] = 1.01;
    controlSystematics["xsig0_ybkg_met3_drmax0"]["vjets"] = 1.07;
    controlSystematics["xsig1_ybkg_met3_drmax0"]["vjets"] = 1.08;
    controlSystematics["xsig0_ybkg_met3_drmax1"]["vjets"] = 1.02;
    controlSystematics["xsig1_ybkg_met3_drmax1"]["vjets"] = 1.00;

    controlSystematics["xsig0_ybkg_met0_drmax0"]["qcd"] = 1.00; 
    controlSystematics["xsig1_ybkg_met0_drmax0"]["qcd"] = 1.00;
    controlSystematics["xsig0_ybkg_met0_drmax1"]["qcd"] = 1.00;
    controlSystematics["xsig1_ybkg_met0_drmax1"]["qcd"] = 1.00;
    controlSystematics["xsig0_ybkg_met1_drmax0"]["qcd"] = 1.00;
    controlSystematics["xsig1_ybkg_met1_drmax0"]["qcd"] = 1.02;
    controlSystematics["xsig0_ybkg_met1_drmax1"]["qcd"] = 1.00;
    controlSystematics["xsig1_ybkg_met1_drmax1"]["qcd"] = 1.00;
    controlSystematics["xsig0_ybkg_met2_drmax0"]["qcd"] = 1.00;
    controlSystematics["xsig1_ybkg_met2_drmax0"]["qcd"] = 1.00;
    controlSystematics["xsig0_ybkg_met2_drmax1"]["qcd"] = 1.00;
    controlSystematics["xsig1_ybkg_met2_drmax1"]["qcd"] = 1.00;
    controlSystematics["xsig0_ybkg_met3_drmax0"]["qcd"] = 1.00;
    controlSystematics["xsig1_ybkg_met3_drmax0"]["qcd"] = 1.00;
    controlSystematics["xsig0_ybkg_met3_drmax1"]["qcd"] = 1.00;
    controlSystematics["xsig1_ybkg_met3_drmax1"]["qcd"] = 1.01;
  }

  void readDimensionFile(std::string const & dimensionFilePath, std::map<std::string, std::vector<std::pair<std::string, std::string> > > & dimensionBins){
    map<string, int> dimensionIndex;
    ifstream dimensionFile(dimensionFilePath.c_str());
    if (dimensionFile.is_open()){
      string dimension, dimensionCut;
      while (dimensionFile >> dimension >> dimensionCut){
        dimensionBins[dimension].push_back(std::make_pair(dimension+to_string(dimensionIndex[dimension]), dimensionCut));
        dimensionIndex[dimension]++;
      }
      dimensionFile.close();
    }
  }

  void writeDataCardHeader(vector<pair<string, string> > sampleBins, ofstream & cardFile)
  {
    cardFile<<"imax "<<sampleBins.size()<<"  number of channels\n";
    cardFile<<"jmax 1  number of backgrounds\n";
    cardFile<<"kmax *  number of nuisance parameters\n";
    cardFile<<"shapes * * FAKE\n\n";
  }

  void setDataCardObserved(map<string, pair<GammaParams, TableRow> > & mYields, vector<pair<string, string> > sampleBins, string const & dataTag, vector<vector<string> > & tableValues){
    // title + type + nBins * 2
    vector<string> row(2+2*sampleBins.size());
    row[0] = "bin";
    for (unsigned iBin=0;iBin<sampleBins.size();++iBin){
      row[iBin*2+3] = sampleBins[iBin].first;
    }
    tableValues.push_back(row);
    setRow(row,"");
  
    row[0] = "Observation";
    for (unsigned iBin=0;iBin<sampleBins.size();++iBin){
      //row[iBin*2+3] = to_string(int(mYields.at(dataTag+"_"+sampleBins[iBin].first).first.Yield())); // Integerize. Round down
      //row[iBin*2+3] = to_string(int(mYields.at(dataTag+"_"+sampleBins[iBin].first).first.Yield())+1); // Integerize. Round up
      //row[iBin*2+3] = RoundNumber(mYields.at(dataTag+"_"+sampleBins[iBin].first).first.Yield(), 0, 1).Data(); // Integerize. Round
      row[iBin*2+3] = to_string(mYields.at(dataTag+"_"+sampleBins[iBin].first).first.Yield());
    }
    tableValues.push_back(row);
    setRow(row,"");
  
    tableValues.push_back(row);
  }

  void setDataCardSignalBackground(string const & processName, string const & signalAverageGenMetTag, map<string, pair<GammaParams, TableRow> > & mYields, vector<pair<string, string> > sampleBins, vector<vector<string> > & tableValues){
    // title + type + nBins * 2
    vector<string> row(2+2*sampleBins.size());
    row[0] = "bin";
    for (unsigned iBin=0;iBin<sampleBins.size();++iBin){
      row[iBin*2+2] = sampleBins[iBin].first;
      row[iBin*2+3] = sampleBins[iBin].first;
    }
    tableValues.push_back(row);
    setRow(row,"");
  
    row[0] = "process";
    for (unsigned iBin=0;iBin<sampleBins.size();++iBin){
      row[iBin*2+2] = "sig";
      row[iBin*2+3] = "bkg";
    }
    tableValues.push_back(row);
    setRow(row,"");
  
    row[0] = "process";
    for (unsigned iBin=0;iBin<sampleBins.size();++iBin){
      row[iBin*2+2] = "0";
      row[iBin*2+3] = "1";
    }
    tableValues.push_back(row);
    setRow(row,"");
  
    row[0] = "rate";
    for (unsigned iBin=0;iBin<sampleBins.size();++iBin){
      string label = processName + "_"+signalAverageGenMetTag+"_" +sampleBins[iBin].first;
      row[iBin*2+2] = to_string(mYields.at(label).first.Yield());
      row[iBin*2+3] = "1";
    }
    tableValues.push_back(row);
    setRow(row,"");
  
    tableValues.push_back(row);
  }
  
  void setDataCardSignalStatistics(string const & processName, string const & signalAverageGenMetTag, map<string, pair<GammaParams, TableRow> > & mYields, vector<pair<string, string> > sampleBins, vector<vector<string> > & tableValues)
  {
    // title + type + nBins * 2
    vector<string> row(2+2*sampleBins.size());
    string label;
  
    for (unsigned rBin=0;rBin<sampleBins.size();++rBin)
    {
      setRow(row,"-");
      row[0] = "stat_"+sampleBins[rBin].first;
      row[1] = "lnN";
      label = processName + "_"+signalAverageGenMetTag+"_" +sampleBins[rBin].first;
      string uncertainty;
      if (mYields.at(label).first.Yield() == 0) uncertainty = "2.00";
      else uncertainty = to_string(1+mYields.at(label).first.Uncertainty()/mYields.at(label).first.Yield());
      row[2+2*rBin] = uncertainty;
      tableValues.push_back(row);
    }
  
    setRow(row,"");
    tableValues.push_back(row);
  }

  void setDataCardControlSystematics(map<string, map<string, float> > controlSystematics, vector<pair<string, string> > sampleBins, vector<vector<string> > & tableValues) 
  {
    // controlSystematics[xsig0_ybkg_met0_drmax0]["ttbar"] = 1.03

    // title + type + nBins * 2
    vector<string> row(2+2*sampleBins.size());
    
    // controlSystematicNames[0] = ttbar
    vector<string> controlSystematicNames;
    for (auto it = controlSystematics.begin()->second.begin(); it != controlSystematics.begin()->second.end(); ++it) controlSystematicNames.push_back(it->first);

    // Set control sample systematics
    for (unsigned iControlSystematic = 0; iControlSystematic < controlSystematicNames.size(); ++iControlSystematic) {
      string systematicName = controlSystematicNames[iControlSystematic]; // ex: ttbar
      setRow(row,"-");
      row[0] = "cr_"+systematicName;
      row[1] = "lnN";
      for (unsigned rBin=0;rBin<sampleBins.size();++rBin) {
        string regionName = sampleBins[rBin].first;
        if (controlSystematics.find(regionName) != controlSystematics.end()) {
          row[2+2*rBin+1] = to_string(controlSystematics[regionName][systematicName]);
        }
      }
      tableValues.push_back(row);
    }

    setRow(row, "");
    tableValues.push_back(row);
  }

  
  // returns count of non-overlapping occurrences of 'sub' in 'str'
  int countSubstring(const std::string& str, const std::string& sub)
  {
    if (sub.length() == 0) return 0;
    int count = 0;
    for (size_t offset = str.find(sub); offset != std::string::npos; offset = str.find(sub, offset + sub.length()))
    {
        ++count;
    }
    return count;
  }

  void setDataCardBackground(map<string, pair<GammaParams, TableRow> > & mYields, vector<pair<string, string> > sampleBins, string const & mcTag, vector<vector<string> > & tableValues)
  {
    // title + type + tag + tagType + (yield,equation), equation arguments
    vector<string> row(6);
    for (unsigned rBin=0;rBin<sampleBins.size();++rBin)
    {
      row[0] = "rp_"+sampleBins[rBin].first;
      row[1] = "rateParam";
      row[2] = sampleBins[rBin].first;
      row[3] = "bkg";
      if (countSubstring(sampleBins[rBin].first,"sig") != 2) 
      {
        //row[4] = to_string(int(mYields.at(mcTag+"_"+sampleBins[rBin].first).first.Yield())); //Integerize. Round down
        //row[4] = to_string(int(mYields.at(mcTag+"_"+sampleBins[rBin].first).first.Yield())+1); //Integerize. Round up
        //row[4] = RoundNumber(mYields.at(mcTag+"_"+sampleBins[rBin].first).first.Yield(), 0, 1).Data(); //Integerize. Round
        row[4] = to_string(mYields.at(mcTag+"_"+sampleBins[rBin].first).first.Yield());
        // Infinity defined in https://root.cern.ch/doc/master/RooNumber_8cxx_source.html
        row[5] = "[0,1.0e30]";
      }
      else 
      {
        //row[4] = "(@0*@1/@2)";
        //vector<string> xydimension;
        //HigUtilities::stringToVectorString(sampleBins[rBin].first, xydimension, "_");
        //vector<string> aNameVector = xydimension;
        //aNameVector[0] = "xbkg";
        //aNameVector[1] = "ybkg";
        //string aName;
        //HigUtilities::vectorStringToString(aNameVector, aName, "_");
        //vector<string> bNameVector = xydimension;
        //bNameVector[1] = "ybkg";
        //string bName;
        //HigUtilities::vectorStringToString(bNameVector, bName, "_");
        //vector<string> cNameVector = xydimension;
        //cNameVector[0] = "xbkg";
        //string cName;
        //HigUtilities::vectorStringToString(cNameVector, cName, "_");
        ////cout<<aName<<" "<<bName<<" "<<cName<<endl;
        //row[5] = "rp_"+cName+",rp_"+bName+",rp_"+aName;

        row[4] = "(@0*@1/@2)*@3";
        vector<string> xydimension;
        HigUtilities::stringToVectorString(sampleBins[rBin].first, xydimension, "_");
        vector<string> aNameVector = xydimension;
        aNameVector[0] = "xbkg";
        aNameVector[1] = "ybkg";
        string aName;
        HigUtilities::vectorStringToString(aNameVector, aName, "_");
        vector<string> bNameVector = xydimension;
        bNameVector[1] = "ybkg";
        string bName;
        HigUtilities::vectorStringToString(bNameVector, bName, "_");
        vector<string> cNameVector = xydimension;
        cNameVector[0] = "xbkg";
        string cName;
        HigUtilities::vectorStringToString(cNameVector, cName, "_");
        string kappaName;
        vector<string> kappaNameVector = xydimension;
        HigUtilities::vectorStringToString(kappaNameVector, kappaName, "_");
        //cout<<aName<<" "<<bName<<" "<<cName<<endl;
        row[5] = "rp_"+cName+",rp_"+bName+",rp_"+aName+",kappa_"+kappaName;
      }
      tableValues.push_back(row);
      setRow(row,"");
    }
  
  }

  void setDataCardKappa(map<string, double > & kappas, map<string, pair<double,double> > & kappa_uncertainties, map<string, vector<pair<string, string> > > & dimensionBins, vector<vector<string> > & tableValues)
  {
    map<string, vector<pair<string, string> > > t_dimensionBins = dimensionBins;
    // Combine dimensions to one list of dimensions.
    HigUtilities::combineDimensionBins(t_dimensionBins);

    auto dimension = t_dimensionBins.begin();
    vector<string> row(4);
    for (unsigned iBin =0; iBin < dimension->second.size(); ++iBin) {
      for (unsigned iXbin = 0; iXbin < 2; iXbin++) {
        //string kappaName = "xsig"+to_string(iXbin)+"_ysig_met"+to_string(iMet)+"_drmax"+to_string(iDrmax);
        string kappaName = "xsig"+to_string(iXbin)+"_ysig_"+(dimension->second)[iBin].first;
        row[0] = "kappa_"+kappaName;
        row[1] = "param";
        row[2] = to_string(kappas[kappaName]);
        row[3] = "-"+to_string(kappa_uncertainties[kappaName].first)+"/+"+to_string(kappa_uncertainties[kappaName].second);
        tableValues.push_back(row);
        setRow(row,"");
      }
    }
    //cout<<"End"<<endl;

    //for (unsigned iMet = 0; iMet < dimensionBins["met"].size(); iMet++) {
    //  for (unsigned iDrmax = 0 ; iDrmax < dimensionBins["drmax"].size(); iDrmax++) {
    //    for (unsigned iXbin = 0; iXbin < 2; iXbin++) {
    //      string kappaName = "xsig"+to_string(iXbin)+"_ysig_met"+to_string(iMet)+"_drmax"+to_string(iDrmax);
    //      row[0] = "kappa_"+kappaName;
    //      row[1] = "param";
    //      row[2] = to_string(kappas[kappaName]);
    //      row[3] = to_string(kappa_uncertainties[kappaName]);
    //      tableValues.push_back(row);
    //      setRow(row,"");
    //    }
    //  }
    //}
  }
  
  void writeTableValues(vector<vector<string> > & tableValues, ofstream & cardFile, bool alignLeft)
  {
    // Set space values
    vector<unsigned> xSpace(tableValues.back().size());
    for (unsigned yIndex = 0; yIndex < tableValues.size(); ++yIndex)
    {
      for (unsigned xIndex = 0; xIndex < tableValues[yIndex].size(); ++xIndex)
      {
        if (xSpace[xIndex] < tableValues[yIndex][xIndex].size()) xSpace[xIndex] = tableValues[yIndex][xIndex].size();
      }
    }
  
    // Write values
    for (unsigned yIndex = 0; yIndex < tableValues.size(); ++yIndex)
    {
      for (unsigned xIndex = 0; xIndex < tableValues[yIndex].size(); ++xIndex)
      {
        if (alignLeft) cardFile << setw(xSpace[xIndex]) << std::left << tableValues[yIndex][xIndex] << " ";
        else cardFile << setw(xSpace[xIndex]) << tableValues[yIndex][xIndex] << " ";
      }
      cardFile<<endl;
    }
  }
  
  void setRow(vector<string> & row, string const & value)
  {
    for(auto & item : row) item = value;
  }

  void GetOptions(int argc, char *argv[])
  {
    while(true){
      static struct option long_options[] = {
        {"output_folder", required_argument, 0, 'o'},
        {"mass_points", required_argument, 0, 'p'},
        {"years", required_argument, 0, 'y'},
        {"luminosity", required_argument, 0, 'l'},
        {"dimension", required_argument, 0, 'd'},
        {"tag", required_argument, 0, 't'},
        {"unblind", no_argument, 0, 'u'},
        {"recomet", no_argument, 0, 0},
        {"higgsino_model", required_argument, 0, 'm'},
        {"priority", required_argument, 0, 'r'},
        {0, 0, 0, 0}
      };
  
      char opt = -1;
      int option_index;
      opt = getopt_long(argc, argv, "o:p:g:y:l:d:t:m:r:nu12", long_options, &option_index);
      if( opt == -1) break;

      string optname;
      switch(opt){
        case 'o': outFolder = optarg; break;
        case 'p': mass_points_string = optarg; break;
        case 'g': mass_point_glob = optarg; break;
        case 'y': years_string = optarg; break;
        case 'l': luminosity = atof(optarg); break;
        case 't': tag = optarg; break;
        case 'm': higgsino_model = optarg; break;
        case 'r': priority = atoi(optarg); break;
        case 'd': 
          dimensionFilePath = optarg; 
          if (!FileExists(dimensionFilePath)) 
          {
            cout<<"[Error] No file caled "<<dimensionFilePath<<endl;
            exit(EXIT_FAILURE);
          }
          break;
        case 'u':
          unblind = true;
          break;
        case '1':
          is1D = true;
          break;
        case '2':
          is1D = false;
          break;
        case 0:
          optname = long_options[option_index].name;
          if(optname == "recomet"){
            do_met_average = false;
          }else{
            printf("Bad option! Found option name %s\n", optname.c_str());
          }
          break;
        default: 
          printf("Bad option! getopt_long returned character code 0%o\n", opt); 
          break;
      }
    }
  }

}
