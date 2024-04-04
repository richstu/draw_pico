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
  const string tag_a = "htozgamma_kingscanyon_v1";
  const string path_to_production_base = "/net/cms11/cms11r0/pico/";
  const string years_string_a[7] = {"2016APV", "2016", "2017", "2018", "2022", "2022EE", "2023"};//add in 2023 BPIX eventually
  //const string years_string_a[4]={"2016APV","2016","2017","2018"};
  //const string years_string_a[1]={"2017"};
  const string skim_a = "merged_*_llg";
  //const string skim_a = "unskimmed";
  float lumi = 1;
}

bool Contains(const string& text, const string& pattern){
  return text.find(pattern) != string::npos;
}

NamedFunc add_cut(NamedFunc & current_cut, NamedFunc additional_cut) {
  current_cut = current_cut && additional_cut;
  return current_cut;
}

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
    if (strcmp(&year, "2023BPIX") == 0) total_luminosity +=9.53;
 }
  
  string total_luminosity_string = RoundNumber(total_luminosity, 1, 1).Data();
  std::cout<<total_luminosity_string<<endl;
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


//This code is for chasing down indices. For nano2pico productions PAST kingscanyon_v1 this should no longer be necessary. Remove if desired.
/*
const NamedFunc get_Zid("get_Zid",[](const Baby &b)->NamedFunc::ScalarType{
  int Zid = -999;
  double dm = 99999.;
  double zmass = 91.1876;
  for (int i = 0; i < b.nll(); i++){
    if((b.ll_charge()->at(i) == 0) && fabs(b.ll_m()->at(i) - zmass) < dm){
      Zid = i;
      dm = fabs(b.ll_m()->at(i) - zmass);
    }
  }
  return Zid;
});

const NamedFunc get_flav("get_flav",[](const Baby &b)->NamedFunc::ScalarType{
  int out = 0;
  if (get_Zid.GetScalar(b) < 0) return out;
  else{
    return b.ll_lepid()->at(get_Zid.GetScalar(b));
  }
});

const NamedFunc get_mll("get_mll",[] (const Baby &b)->NamedFunc::ScalarType{
  double out = -1;
  if (get_Zid.GetScalar(b) < 0) return out;
  else{
    return b.ll_m()->at(get_Zid.GetScalar(b));
  }
});

const NamedFunc get_yid("get_yid",[](const Baby &b)->NamedFunc::ScalarType{
  int yid = -999;
  double ypt = -1.;
  for (int i = 0; i < b.nphoton(); i++){
    if(b.photon_id80()->at(i) && b.photon_pt()->at(i)>ypt){
      yid = i;
      ypt = b.photon_pt()->at(i);
      //break;
    }
  }
  return yid;
});

const NamedFunc wp80fix("wp80fix",[](const Baby &b)->NamedFunc::ScalarType{
  bool WP = false;
  if (get_yid.GetScalar(b) < 0) return WP;
  else if (b.photon_id80()->at(get_yid.GetScalar(b))){
    WP = true;
  }
  return WP;
});

const NamedFunc get_Hid("get_Hid",[](const Baby &b)->NamedFunc::ScalarType{
  int Hid = -999;
  if (get_Zid.GetScalar(b) < 0 || get_yid.GetScalar(b) < 0) return Hid;
  else{
    for (int i = 0; i < b.nllphoton(); i++){
      if(b.llphoton_ill()->at(i) == get_Zid.GetScalar(b) && b.llphoton_iph()->at(i) == get_yid.GetScalar(b)){
        Hid = i;
        break;
      }
    }
    return Hid;
  }
});

const NamedFunc get_mlly("get_mlly",[](const Baby &b)->NamedFunc::ScalarType{
  double out = -1;
  if(get_Hid.GetScalar(b) < 0) return out;
  else{
    return b.llphoton_m()->at(get_Hid.GetScalar(b));
  }
});

const NamedFunc get_ratio("get_ratio",[](const Baby &b)->NamedFunc::ScalarType{
  double out = -1;
  if(get_Hid.GetScalar(b) < 0) return out;
  else{
    return b.photon_pt()->at(get_yid.GetScalar(b))/b.llphoton_m()->at(get_Hid.GetScalar(b));
  }
});

const NamedFunc get_sum("get_sum",[](const Baby &b)->NamedFunc::ScalarType{
  double out = -1;
  if(get_Hid.GetScalar(b) < 0) return out;
  else{
    return b.ll_m()->at(get_Zid.GetScalar(b)) + b.llphoton_m()->at(get_Hid.GetScalar(b));
  }
});
*/

  //This code defines the triggers and trigger pt selections. The muon triggers may be incorrect, but I tried to add all the triggers it could be
  const NamedFunc trigs_2016_e = "HLT_Ele27_WPTight_Gsf";                  const NamedFunc trigs_2016_ee = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ";
  //const NamedFunc trigs_2017_e = "HLT_Ele32_WPTight_Gsf_L1DoubleEG";       const NamedFunc trigs_2017_ee = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"; //Use for kingscanyon_v1
  const NamedFunc trigs_2017_e = "HLT_Ele32_WPTight_Gsf_Emu";              const NamedFunc trigs_2017_ee = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"; //use for all versions past kingscanyon_v1
  const NamedFunc trigs_2018_e = "HLT_Ele32_WPTight_Gsf";                  const NamedFunc trigs_2018_ee = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL";
  const NamedFunc trigs_2022_e = "HLT_Ele30_WPTight_Gsf";                  const NamedFunc trigs_2022_ee = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL";
  const NamedFunc trigs_2023_e = "HLT_Ele30_WPTight_Gsf";                  const NamedFunc trigs_2023_ee = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL";

  const NamedFunc trigs_2016_mu = "HLT_IsoMu24 || HLT_IsoTkMu24";          const NamedFunc trigs_2016_mumu = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ";
  const NamedFunc trigs_2017_mu = "HLT_IsoMu27";                           const NamedFunc trigs_2017_mumu = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8";
  const NamedFunc trigs_2018_mu = "HLT_IsoMu24";                           const NamedFunc trigs_2018_mumu = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8";
  const NamedFunc trigs_2022_mu = "HLT_IsoMu24";                           const NamedFunc trigs_2022_mumu = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8";
  const NamedFunc trigs_2023_mu = "HLT_IsoMu24";                           const NamedFunc trigs_2023_mumu = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8";

  const NamedFunc trigs_pT_2016_el = (trigs_2016_e && signal_lead_electron_pt > 30) || (trigs_2016_ee && signal_lead_electron_pt > 25 && signal_sublead_electron_pt > 15);
  const NamedFunc trigs_pT_2016_mu = (trigs_2016_mu && signal_lead_muon_pt > 25) || (trigs_2016_mumu && signal_lead_muon_pt > 20 && signal_sublead_muon_pt > 10);
  const NamedFunc trigs_pT_2016 = trigs_pT_2016_el || trigs_pT_2016_mu;

  const NamedFunc trigs_pT_2017_el = (trigs_2017_e && signal_lead_electron_pt > 35) || (trigs_2017_ee && signal_lead_electron_pt > 25 && signal_sublead_electron_pt > 15); 
  const NamedFunc trigs_pT_2017_mu = (trigs_2017_mu && signal_lead_muon_pt > 28) || (trigs_2017_mumu && signal_lead_muon_pt > 20 && signal_sublead_muon_pt > 10);
  const NamedFunc trigs_pT_2017 = trigs_pT_2017_el || trigs_pT_2017_mu;

  const NamedFunc trigs_pT_2018_el = (trigs_2018_e && signal_lead_electron_pt > 35) || (trigs_2018_ee && signal_lead_electron_pt > 25 && signal_sublead_electron_pt > 15);
  const NamedFunc trigs_pT_2018_mu = (trigs_2018_mu && signal_lead_muon_pt > 25) || (trigs_2018_mumu && signal_lead_muon_pt > 20 && signal_sublead_muon_pt > 10);
  const NamedFunc trigs_pT_2018 = trigs_pT_2018_el || trigs_pT_2018_mu;

  const NamedFunc trigs_pT_2022_el = (trigs_2022_e && signal_lead_electron_pt > 35) || (trigs_2022_ee && signal_lead_electron_pt > 25 && signal_sublead_electron_pt > 15);
  const NamedFunc trigs_pT_2022_mu = (trigs_2022_mu && signal_lead_muon_pt > 25) || (trigs_2022_mumu && signal_lead_muon_pt > 20 && signal_sublead_muon_pt > 10);
  const NamedFunc trigs_pT_2022 = trigs_pT_2022_el || trigs_pT_2022_mu;

  const NamedFunc trigs_pT_2023_el = (trigs_2023_e && signal_lead_electron_pt > 35) || (trigs_2023_ee && signal_lead_electron_pt > 25 && signal_sublead_electron_pt > 15);
  const NamedFunc trigs_pT_2023_mu = (trigs_2023_mu && signal_lead_muon_pt > 25) || (trigs_2023_mumu && signal_lead_muon_pt > 20 && signal_sublead_muon_pt > 10);
  const NamedFunc trigs_pT_2023 = trigs_pT_2023_el || trigs_pT_2023_mu;

  //Electron trigs for run 2 and run 3, no pT cut
  const NamedFunc trigs_el("trigs_el", [](const Baby &b)->NamedFunc::ScalarType{
      bool res = false;
      if (abs(b.SampleType()) == 2016){
          res = trigs_2016_e.GetScalar(b) || trigs_2016_ee.GetScalar(b);
      } else if (abs(b.SampleType()) == 2017){
          res = trigs_2017_e.GetScalar(b) || trigs_2017_ee.GetScalar(b);
      } else if (abs(b.SampleType()) == 2018){
          res = trigs_2018_e.GetScalar(b) || trigs_2018_ee.GetScalar(b);
      } else if (abs(b.SampleType()) == 2022){
          res = trigs_2022_e.GetScalar(b) || trigs_2022_ee.GetScalar(b);
      } else if (abs(b.SampleType()) == 2023){
          res = trigs_2023_e.GetScalar(b) || trigs_2023_ee.GetScalar(b);
      };
      return res;
  });

  //Electron trigs for run 2 and 3, pT cut included
  const NamedFunc trigs_el_pT("trigs_el_pT", [](const Baby &b)->NamedFunc::ScalarType{
      bool res = false;
      if (abs(b.SampleType()) == 2016){
          res = trigs_pT_2016_el.GetScalar(b);
      } else if (abs(b.SampleType()) == 2017){
          res = trigs_pT_2017_el.GetScalar(b);
      } else if (abs(b.SampleType()) == 2018){
          res = trigs_pT_2018_el.GetScalar(b);
      } else if (abs(b.SampleType()) == 2022){
          res = trigs_pT_2022_el.GetScalar(b);
      } else if (abs(b.SampleType()) == 2023){
          res = trigs_pT_2023_el.GetScalar(b);
      };
      return res;
  });

  //Muon trigs for run 2 and run 3, no pT cut
  const NamedFunc trigs_mu("trigs_mu", [](const Baby &b)->NamedFunc::ScalarType{
      bool res = false;
      if (abs(b.SampleType()) == 2016){
          res = trigs_2016_mu.GetScalar(b) || trigs_2016_mumu.GetScalar(b);
      } else if (abs(b.SampleType()) == 2017){
          res = trigs_2017_mu.GetScalar(b) || trigs_2017_mumu.GetScalar(b);
      } else if (abs(b.SampleType()) == 2018){
          res = trigs_2018_mu.GetScalar(b) || trigs_2018_mumu.GetScalar(b);
      } else if (abs(b.SampleType()) == 2022){
          res = trigs_2022_mu.GetScalar(b) || trigs_2022_mumu.GetScalar(b);
      } else if (abs(b.SampleType()) == 2023){
          res = trigs_2023_mu.GetScalar(b) || trigs_2023_mumu.GetScalar(b);
      };
      return res;
  });

  //Muon trigs for run 2 and run 3, pT cut included
  const NamedFunc trigs_mu_pT("trigs_mu_pT", [](const Baby &b)->NamedFunc::ScalarType{
      bool res = false;
      if (abs(b.SampleType()) == 2016){
          res = trigs_pT_2016_mu.GetScalar(b);
      } else if (abs(b.SampleType()) == 2017){
          res = trigs_pT_2017_mu.GetScalar(b);
      } else if (abs(b.SampleType()) == 2018){
          res = trigs_pT_2018_mu.GetScalar(b);
      } else if (abs(b.SampleType()) == 2022){
          res = trigs_pT_2022_mu.GetScalar(b);
      } else if (abs(b.SampleType()) == 2023){
          res = trigs_pT_2023_mu.GetScalar(b);
      };
      return res;
  });

  const NamedFunc trigs_Run2_pT("trigs_Run2_pT", [](const Baby &b)->NamedFunc::ScalarType{
      bool res = false;
      if (abs(b.SampleType()) == 2016){
          res = trigs_pT_2016.GetScalar(b);
      } else if (abs(b.SampleType()) == 2017){
          res = trigs_pT_2017.GetScalar(b);
      } else if (abs(b.SampleType()) == 2018){
          res = trigs_pT_2018.GetScalar(b);
      };
      return res;
  });

  const NamedFunc trigs_Run3_pT("trigs_Run3_pT", [](const Baby &b)->NamedFunc::ScalarType{
      bool res = false;
      if (abs(b.SampleType()) == 2022){
          res = trigs_pT_2022.GetScalar(b);
      } else if (abs(b.SampleType()) == 2023){
          res = trigs_pT_2023.GetScalar(b);
      };
      return res;
  });


void constructCutflowTable(vector<TableRow> & tablerows, NamedFunc weight, int electron_or_muon, bool isMC = false) {
    NamedFunc current_cut = "1";
    if (isMC == true){
      current_cut = "use_event";
    }
    //NamedFunc el_trigger = "ll_lepid[0]==11 && (HLT_Ele27_WPTight_Gsf || HLT_Ele30_WPTight_Gsf || HLT_Ele35_WPTight_Gsf || HLT_Ele32_WPTight_Gsf || HLT_Ele32_WPTight_Gsf_L1DoubleEG || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL)"; //run 2 and run 3 e, ee triggers. Add in Run 3
    //NamedFunc mu_trigger = "ll_lepid[0]==13 && (HLT_IsoMu24 || HLT_IsoTkMu24 || HLT_IsoMu27 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8)"; //run 2 and run 3 mu, mumu triggers. Add in Run 3
    if (electron_or_muon == 0) { // el
      tablerows = vector<TableRow>{
        TableRow("Initial events", add_cut(current_cut,"1"),0,0, weight),
        TableRow("$N_e \\geq 2$", add_cut(current_cut, "nel>=2 && ll_lepid[0]==11"),0,0, weight),
        TableRow("e, ee trigger", add_cut(current_cut, trigs_el),0,0, weight),
        TableRow("el trigger $p_{T}$ cuts", add_cut(current_cut, trigs_el_pT),0,0, weight),
        TableRow("$N_{\\gamma} \\geq 1$", add_cut(current_cut, "nphoton>=1"),0,0, weight),
        TableRow("$80 \\text{ GeV} < m_{ee} < 100 \\text{ GeV}$", add_cut(current_cut, "ll_m[0] > 80 && ll_m[0]<100"),0,0, weight),
        TableRow("$p_{T}^{\\gamma}/m_{ee\\gamma} > 15./110$", add_cut(current_cut, "(photon_pt[0]/llphoton_m[0])>(15./110)"),0,0, weight),
        TableRow("$m_{ee} + m_{ee\\gamma} > 185 \\text{ GeV}$", add_cut(current_cut, "(llphoton_m[0] + ll_m[0]) > 185"),0,0, weight),
        TableRow("$100 < m_{ee\\gamma} < 180 \\text{ GeV}$", add_cut(current_cut, "llphoton_m[0] > 100 && llphoton_m[0] < 180"),0,0, weight),
        TableRow("ggF cat",current_cut && "nlep == 2 && njet <= 1 && met < 90",0,0,weight),
        TableRow("VBF cat",current_cut && "nlep == 2 && njet >= 2 && nbdfm == 0",0,0,weight),
        TableRow("ZH/WH cat",current_cut && "nlep >= 3 && nbdfm == 0",0,0,weight),
        TableRow("ZH pt miss",current_cut && "nlep == 2 && njet <= 1 && met > 90",0,0,weight),
        TableRow("ttH hadronic",current_cut && "nlep == 2 && njet >= 5 && nbdfm >=1",0,0,weight),
        TableRow("ttH leptonic",current_cut && "(nlep == 3 && njet >=3 && nbdfm >=1) || (nlep>=4 && njet >=1 && nbdfm>=1)",0,0,weight),

      };
    } else if (electron_or_muon == 1) { // mu
      tablerows = vector<TableRow>{
        TableRow("Initial events", add_cut(current_cut,"1"),0,0, weight),
        TableRow("$N_{\\mu} \\geq 2 $", add_cut(current_cut, "nmu>=2 && ll_lepid[0]==13"),0,0, weight),
        TableRow("$\\mu, \\mu\\mu$ trigger", add_cut(current_cut, trigs_mu),0,0, weight),
        TableRow("mu trigger $p_{T}$ cuts", add_cut(current_cut, trigs_mu_pT),0,0, weight),
        TableRow("$N_{\\gamma} \\geq 1$", add_cut(current_cut, "nphoton>=1"),0,0, weight),
        TableRow("$80 \\text{ GeV} < m_{ee} < 100 \\text{ GeV}$", add_cut(current_cut, "ll_m[0] > 80 && ll_m[0]<100"),0,0, weight),
        TableRow("$p_{T}^{\\gamma}/m_{ee\\gamma} > 15./110$", add_cut(current_cut, "(photon_pt[0]/llphoton_m[0])>(15./110)"),0,0, weight),
        TableRow("$m_{ee} + m_{ee\\gamma} > 185 \\text{ GeV}$", add_cut(current_cut, "(llphoton_m[0] + ll_m[0]) > 185"),0,0, weight),
        TableRow("$100 < m_{ee\\gamma} < 180 \\text{ GeV}$", add_cut(current_cut, "llphoton_m[0] > 100 && llphoton_m[0] < 180"),0,0, weight),
        TableRow("ggF cat",current_cut && "nlep == 2 && njet <= 1 && met < 90",0,0,weight),
        TableRow("VBF cat",current_cut && "nlep == 2 && njet >= 2 && nbdfm == 0",0,0,weight),
        TableRow("ZH/WH cat",current_cut && "nlep >= 3 && nbdfm == 0",0,0,weight),
        TableRow("ZH pt miss",current_cut && "nlep == 2 && njet <= 1 && met > 90",0,0,weight),
        TableRow("ttH hadronic",current_cut && "nlep == 2 && njet >= 5 && nbdfm >=1",0,0,weight),
        TableRow("ttH leptonic",current_cut && "(nlep == 3 && njet >=3 && nbdfm >=1) || (nlep>=4 && njet >=1 && nbdfm>=1)",0,0,weight),

      };
    } else { // mu + el
      tablerows = vector<TableRow>{
        TableRow("Initial events", add_cut(current_cut,"1"),0,0, weight),
        TableRow("$N_e \\geq 2 || N_{\\mu} \\geq 2 $", add_cut(current_cut, "(nel>=2 && ll_lepid[0]==11) || (nmu>=2 && ll_lepid[0]==13)"),0,0, weight),
        TableRow("e, ee trigger $||$ $\\mu, \\mu\\mu$ trigger", add_cut(current_cut, trigs_mu || trigs_el),0,0, weight),
        TableRow("mu and el trigger $p_{T}$ cuts", add_cut(current_cut, (trigs_mu_pT || trigs_el_pT)),0,0, weight),
        TableRow("$N_{\\gamma} \\geq 1$", add_cut(current_cut, "nphoton>=1"),0,0, weight),
        TableRow("$80 \\text{ GeV} < m_{ee} < 100 \\text{ GeV}$", add_cut(current_cut, "ll_m[0] > 80 && ll_m[0]<100"),0,0, weight),
        TableRow("$p_{T}^{\\gamma}/m_{ee\\gamma} > 15./110$", add_cut(current_cut, "(photon_pt[0]/llphoton_m[0])>(15./110)"),0,0, weight),
        TableRow("$m_{ee} + m_{ee\\gamma} > 185 \\text{ GeV}$", add_cut(current_cut, "(llphoton_m[0] + ll_m[0]) > 185"),0,0, weight),
        TableRow("$100 < m_{ee\\gamma} < 180 \\text{ GeV}$", add_cut(current_cut, "llphoton_m[0] > 100 && llphoton_m[0] < 180"),0,0, weight),
        TableRow("ggF cat",current_cut && "nlep == 2 && njet <= 1 && met < 90",0,0,weight),
        TableRow("VBF cat",current_cut && "nlep == 2 && njet >= 2 && nbdfm == 0",0,0,weight),
        TableRow("ZH/WH cat",current_cut && "nlep >= 3 && nbdfm == 0",0,0,weight),
        TableRow("ZH pt miss",current_cut && "nlep == 2 && njet <= 1 && met > 90",0,0,weight),
        TableRow("ttH hadronic",current_cut && "nlep == 2 && njet >= 5 && nbdfm >=1",0,0,weight),
        TableRow("ttH leptonic",current_cut && "(nlep == 3 && njet >=3 && nbdfm >=1) || (nlep>=4 && njet >=1 && nbdfm>=1)",0,0,weight),
      };
    }
}

void generateTable(const string path_to_production_a, PlotMaker& pm, Palette colors, set<string> year_t, const string skim_a, NamedFunc weight_a, TString tablename_a, 
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
    //filenames = {"*DYJetsToLL*amcatnloFXFX*.root"};
    if(path_to_production_a == "/net/cms11/cms11r0/pico/NanoAODv11/htozgamma_kingscanyon_v1/"|| path_to_production_a == "/net/cms11/cms11r0/pico/NanoAODv11p9/htozgamma_kingscanyon_v1/"){
      filenames = {"*DYto2L*amcatnloFXFX*.root", "*DYGto2LG-1Jets_MLL-50*.root"}; //
    }
    type = "mc/"; label = "DYJets and SMZy"; color = "t1tttt"; print_uncertainty = true;
    pathNames = attach_folder(path_to_production_a, year_t, type + skim_a, filenames); 
    procs_dat_a.push_back(Process::MakeShared<Baby_pico>(label, Process::Type::signal, colors(color), pathNames, "1")); 
  } else if(isMC == true && isSig == true){
    //filenames = {"*GluGlu*HToZG*ZToLL_M-125_*.root","VBF*HToZG*ZToLL_M-125_*.root","*WplusH*HToZG*ZToLL_M-125_*.root", "*WminusH*HToZG*ZToLL_M-125_*.root", "*tt*HToZG*ZToLL_M-125_*.root", "*ZH*HToZG_ZToAll_M-125_*.root"};
    filenames = {"*GluGlu*HToZG*ZToLL_M-125_*.root","VBF*HToZG*ZToLL_M-125_*.root", "*tt*HToZG*ZToLL_M-125_*.root", "*ZH*HToZG_ZToAll_M-125_*.root"};//add back in Wplus and Wminus later
    //filenames = {"*GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8__RunIISummer20UL17NanoAODv9__106X_mc2017_realistic_v9-v1*.root"};
    //filenames = {"*GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV*.root"};
    type = "mc/"; label = "Signal MC"; color = "dy"; print_uncertainty = true;
    pathNames = attach_folder(path_to_production_a, year_t, type + skim_a, filenames); 
    procs_dat_a.push_back(Process::MakeShared<Baby_pico>(label, Process::Type::background, colors(color), pathNames, "1"));
  } else if(isMC == false){
    filenames = {"*DoubleEG*.root","*DoubleMuon*.root","*SingleElectron*.root","*SingleMuon*.root","*EGamma*.root","*_Muon*.root"};
    type = "data/"; label = "Data"; color = "data"; print_uncertainty = false;
    pathNames = attach_folder(path_to_production_a, year_t, type + skim_a, filenames); 
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
  //This script creates cutflows of the Run 2 and Run 3 data for validation purposes using the full standard baseline selection. Add in MC later. 
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches


  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////// Initializing paths and selections ////////////////////////////

  string path_to_production_3p = path_to_production_base +"NanoAODv11p9/" + tag_a + "/";
  string path_to_production_3 = path_to_production_base +"NanoAODv11/" + tag_a + "/";
  string path_to_production_2 = path_to_production_base +"NanoAODv9/" + tag_a +"/";
  string path_to_production_a;
  string lep_tag;
  set<string> year_t; 

  NamedFunc weight_a = "1";
  Palette colors("txt/colors.txt", "default");
  PlotMaker pm; //Define just once


  /////////////////////////////////////// Data and MC Validation ////////////////////////////////////
  //////////////////////////////////// Loop over years and leptons //////////////////////////////////
  for (unsigned int k = 0; k < size(years_string_a); k++){
    //Loop over all datasets, and make separate cutflow table for each.
    if (Contains(years_string_a[k], "2023")){
      path_to_production_a = path_to_production_3p;
    } else if (Contains(years_string_a[k], "2022")){
      path_to_production_a = path_to_production_3;
    } else {
      path_to_production_a = path_to_production_2;
    }
    for (unsigned int l = 0; l<=2; l++){
      //Loop over lepton types, 0 = electron, 1 = muon, 2 = electron and muon

        if (l==0) lep_tag = "_electron";
        else if (l==1) lep_tag = "_muon";
        else  lep_tag = "_lepton";
        TString tablename_a = "FixName:validation_table_"+tag_a+lep_tag+"_"+years_string_a[k];
        string total_luminosity_string_a = getLuminosityString(years_string_a[k]);
        year_t = {years_string_a[k]};


        //    Generate data cutflow table
        //--------------------------------------------------
        tablename_a = "FixName:validation_table_"+tag_a+lep_tag+"_"+years_string_a[k]+"_data";
        weight_a = "1";
        generateTable(path_to_production_a, pm, colors, year_t, skim_a, weight_a, tablename_a, total_luminosity_string_a, lumi, l, false, false);

        //    Generate background MC (no SF) cutflow table
        //--------------------------------------------------
        tablename_a = "FixName:validation_table_"+tag_a+lep_tag+"_"+years_string_a[k]+"_backmc";
        weight_a = "1";
        generateTable(path_to_production_a, pm, colors, year_t, skim_a, weight_a, tablename_a, total_luminosity_string_a, lumi, l, true, false);

        //    Generate background MC (w/ SF) cutflow table
        //--------------------------------------------------
        tablename_a = "FixName:validation_table_"+tag_a+lep_tag+"_"+years_string_a[k]+"_backmc_sf";
        weight_a = wgt*w_years;//can change to w_lumi for testing
        generateTable(path_to_production_a, pm, colors, year_t, skim_a, weight_a, tablename_a, total_luminosity_string_a, lumi, l, true, false);

        //    Generate signal MC (no SF) cutflow table
        //--------------------------------------------------
        tablename_a = "FixName:validation_table_"+tag_a+lep_tag+"_"+years_string_a[k]+"_sigmc";
        weight_a = "1";
        generateTable(path_to_production_a, pm, colors, year_t, skim_a, weight_a, tablename_a, total_luminosity_string_a, lumi, l, true, true);

        //    Generate signal MC (no SF) cutflow table
        //--------------------------------------------------
        tablename_a = "FixName:validation_table_"+tag_a+lep_tag+"_"+years_string_a[k]+"_sigmc_sf";
        weight_a = wgt*w_years;//can change to w_lumi for testing
        generateTable(path_to_production_a, pm, colors, year_t, skim_a, weight_a, tablename_a, total_luminosity_string_a, lumi, l, true, true);

    }
  }



}

