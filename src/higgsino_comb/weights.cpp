#include <algorithm>
#include <stdlib.h>
#include <regex>
#include <string>
#include "TVector2.h"
#include "core/hist1d.hpp"
#include "core/utilities.hpp"
#include "core/palette.hpp"
#include "core/plot_opt.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/utilities.hpp"
#include "core/functions.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"
#include "higgsino_comb/vardef.hpp"
#include "higgsino_comb/regions.hpp"
#include "higgsino_comb/weights.hpp"

using namespace std;

namespace weights{

  const NamedFunc w_run2("w_run2", [](const Baby &b) -> NamedFunc::ScalarType{
    double w_year = 1;
    if(b.SampleTypeString().Contains("-")) w_year = 1;
    else{
      bool isdiph = false;
      bool isQCD2 = false;
      bool isttWToLNu = false;
      bool isSMS = false;
      std::string file = *b.FileNames().begin();
      if (file.find("DiPhotonJetsBox") != file.npos) isdiph = true;
      if (file.find("QCD_Pt-40") != file.npos) isQCD2 = true;
      if (file.find("TTWJetsToLNu") != file.npos) isttWToLNu = true;
      if (file.find("TChiHH") != file.npos) isSMS = true;
      if (b.SampleTypeString()=="2016") {
          w_year = 19.5;
      } else if (b.SampleTypeString()=="2016APV"){
          w_year = 16.8;
      } else if (b.SampleTypeString()=="2017"){
          w_year = 41.48;
      } else if (b.SampleTypeString()=="2018"){
          w_year = 59.83;
      }
      if (isSMS) w_year = 137.61*vardef::xsec_ratio.GetScalar(b); //since we have only 2016 v7 MC
      w_year = w_year*b.w_lumi();//*b.w_trig();
      if(isdiph) w_year *= b.mc_weight()>0 ? b.mc_weight():-b.mc_weight(); /*w_year = w_year*b.mc_weight();*/ //w_lumi already knows the sign
      if(isQCD2) w_year *= -0.1131001131; //due to it being unable to find cross section in nano2pico
      if (isttWToLNu) w_year *= 0.2163/66680; // remove after 4b background has been processed again
    }
    return w_year;
  });

  const NamedFunc w_lumi("w_lumi", [](const Baby &b) -> NamedFunc::ScalarType{
    double w_year = 1;
    if(b.SampleTypeString().Contains("-")) w_year = 1;
    else{
      bool isdiph = false;
      bool isQCD2 = false;
      bool isttWToLNu = false;
      bool isSMS = false;
      std::string file = *b.FileNames().begin();
      if (file.find("DiPhotonJetsBox") != file.npos) isdiph = true;
      if (file.find("QCD_Pt-40") != file.npos) isQCD2 = true;
      if (file.find("TTWJetsToLNu") != file.npos) isttWToLNu = true;
      if (file.find("TChiHH") != file.npos) isSMS = true;
      if (b.SampleTypeString()=="2016") {
          w_year = 19.5;
      } else if (b.SampleTypeString()=="2016APV"){
          w_year = 16.8;
      } else if (b.SampleTypeString()=="2017"){
          w_year = 41.48;
      } else if (b.SampleTypeString()=="2018"){
          w_year = 59.83;
      }
      if (isSMS) w_year = 137.61*vardef::xsec_ratio.GetScalar(b); //since we have only 2016 v7 MC
      w_year = w_year*b.w_lumi();//*b.w_trig();
      if(isdiph) w_year *= b.mc_weight()>0 ? b.mc_weight():-b.mc_weight(); /*w_year = w_year*b.mc_weight();*/ //w_lumi already knows the sign
      if(isQCD2) w_year *= -0.1131001131; //due to it being unable to find cross section in nano2pico
      if (isttWToLNu) w_year *= 0.2163/66680; // remove after 4b background has been processed again
      w_year *= 200/137.61; // scale for 200fb-1
    }
    return w_year;
  });

}

namespace triggers{
  const NamedFunc single_ele_trig("single_ele_trig", [](const Baby &b) -> NamedFunc::ScalarType{
    bool trig = b.HLT_Ele25_WPTight_Gsf() || b.HLT_Ele27_WPTight_Gsf() || b.HLT_Ele32_WPTight_Gsf() || b.HLT_Ele35_WPTight_Gsf() || b.HLT_Ele32_WPTight_Gsf_L1DoubleEG() || b.HLT_Ele20_WPLoose_Gsf() || b.HLT_Ele45_WPLoose_Gsf() || b.HLT_Ele105_CaloIdVT_GsfTrkIdT() || b.HLT_Ele115_CaloIdVT_GsfTrkIdT() || b.HLT_Ele135_CaloIdVT_GsfTrkIdT() || b.HLT_Ele145_CaloIdVT_GsfTrkIdT() || b.HLT_Ele25_eta2p1_WPTight_Gsf() || b.HLT_Ele27_eta2p1_WPTight_Gsf() || b.HLT_Ele20_eta2p1_WPLoose_Gsf() || b.HLT_Ele25_eta2p1_WPLoose_Gsf() || b.HLT_Ele27_eta2p1_WPLoose_Gsf() || b.HLT_Ele15_IsoVVVL_PFHT350() || b.HLT_Ele15_IsoVVVL_PFHT400() || b.HLT_Ele15_IsoVVVL_PFHT450() || b.HLT_Ele15_IsoVVVL_PFHT600();
    return trig;
  });

  const NamedFunc single_muon_trig("single_muon_trig", [](const Baby &b) -> NamedFunc::ScalarType{
    bool trig = b.HLT_IsoMu20() || b.HLT_IsoMu22() || b.HLT_IsoMu24() || b.HLT_IsoMu27() || b.HLT_IsoTkMu20() || b.HLT_IsoTkMu22() || b.HLT_IsoTkMu24() || b.HLT_Mu50() || b.HLT_Mu55() || b.HLT_TkMu50() || b.HLT_IsoMu22_eta2p1() || b.HLT_IsoMu24_eta2p1() ||b.HLT_Mu45_eta2p1() || b.HLT_Mu15_IsoVVVL_PFHT350() || b.HLT_Mu15_IsoVVVL_PFHT400() || b.HLT_Mu15_IsoVVVL_PFHT450()|| b.HLT_Mu15_IsoVVVL_PFHT600() || b.HLT_Mu50_IsoVVVL_PFHT400();
    return trig;
  });

  const NamedFunc ptmiss_htmiss_trig("ptmiss_htmiss_trig", [](const Baby &b) -> NamedFunc::ScalarType{
    bool trig = b.HLT_PFMET90_PFMHT90_IDTight() || b.HLT_PFMET100_PFMHT100_IDTight() || b.HLT_PFMET110_PFMHT110_IDTight() || b.HLT_PFMET120_PFMHT120_IDTight() || b.HLT_PFMET130_PFMHT130_IDTight() || b.HLT_PFMET140_PFMHT140_IDTight() || b.HLT_PFMETNoMu90_PFMHTNoMu90_IDTight() || b.HLT_PFMETNoMu100_PFMHTNoMu100_IDTight() || b.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight() || b.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight() || b.HLT_PFMETNoMu130_PFMHTNoMu130_IDTight() || b.HLT_PFMETNoMu140_PFMHTNoMu140_IDTight() || b.HLT_PFMET100_PFMHT100_IDTight_PFHT60() || b.HLT_PFMET110_PFMHT110_IDTight_PFHT60() || b.HLT_PFMET120_PFMHT120_IDTight_PFHT60() || b.HLT_PFMET130_PFMHT130_IDTight_PFHT60() || b.HLT_PFMET140_PFMHT140_IDTight_PFHT60() || b.HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60() || b.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_PFHT60() || b.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60() || b.HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_PFHT60() || b.HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_PFHT60() || b.HLT_PFMET120_PFMHT120_IDTight_HFCleaned() || b.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned() || b.HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned();
    return trig;
  });

  const NamedFunc jet_ht_trig("jet_ht_trig", [](const Baby &b) -> NamedFunc::ScalarType{
    bool trig = b.HLT_PFJet500() || b.HLT_PFHT125() || b.HLT_PFHT200() || b.HLT_PFHT300() || b.HLT_PFHT400() || b.HLT_PFHT475() || b.HLT_PFHT600() || b.HLT_PFHT650() || b.HLT_PFHT800() || b.HLT_PFHT900() || b.HLT_PFHT180() || b.HLT_PFHT370() || b.HLT_PFHT430() || b.HLT_PFHT510() || b.HLT_PFHT590() || b.HLT_PFHT680() || b.HLT_PFHT780() || b.HLT_PFHT890() || b.HLT_PFHT1050() || b.HLT_PFHT250() || b.HLT_PFHT350();
    return trig;
  });

}

namespace filters{
  const NamedFunc pass_hemveto("pass_hemveto", [](const Baby &b) -> NamedFunc::ScalarType{
    //only apply for 2018 era C+D and MC
    if (abs(b.SampleType())!=2018) return true; // accept 2016 and 2017 data and mc
    if (b.SampleType()==-2018 && b.run() < 319077) return true; // accept part of 2018 data
    if (b.SampleType()==2018 && (b.event()%1961) >= 1296) return true; //accept part of 2018 mc
    bool pass_hem = true;
    for (unsigned int el_idx = 0; el_idx < b.el_pt()->size(); el_idx++) {
      if (b.el_miniso()->at(el_idx) < 0.1 && -3.0 < b.el_eta()->at(el_idx) && b.el_eta()->at(el_idx) < -1.4 && -1.57 < b.el_phi()->at(el_idx) && b.el_phi()->at(el_idx) < -0.87) {
        pass_hem = false;
      }
    }
    for (unsigned int jet_idx = 0; jet_idx < b.jet_pt()->size(); jet_idx++) {
      if (b.jet_pt()->at(jet_idx) > 30. && -3.2 < b.jet_eta()->at(jet_idx) && b.jet_eta()->at(jet_idx) < -1.2 && -1.77 < b.jet_phi()->at(jet_idx) && b.jet_phi()->at(jet_idx) < -0.67) {
        double dphi = fabs(TVector2::Phi_mpi_pi(b.jet_phi()->at(jet_idx)-b.met_phi()));
        if (dphi < 0.5) {
          pass_hem = false;
        }
      }
    }
    return pass_hem;
  });

  const NamedFunc pass_filters("pass_filters", [](const Baby &b) -> NamedFunc::ScalarType{
    if (!b.pass_goodv() || !b.pass_hbhe() || !b.pass_hbheiso() || !b.pass_ecaldeadcell() || !b.pass_badpfmu() || !b.pass_muon_jet()) return false;
    if (b.type()/1000 == 0 && !b.pass_eebadsc()) return false; //only apply eebadsc fiter for data
    if (!b.FastSim() && !b.pass_cschalo_tight()) return false; //not for fastsim
    if (!b.pass_low_neutral_jet()) return false;
    if (!b.pass_htratio_dphi_tight()) return false;
    if (b.FastSim() && !b.pass_jets()) return false; //back to only for fastsim
    //if (!b.pass_jets()) return false; //was modified
    if (!b.pass_ecalnoisejet()) return false;
    if (!pass_hemveto.GetScalar(b)) return false;
    return true;
  });
 

}
