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

  const NamedFunc w_run2_datavmc("w_run2_datavmc", [](const Baby &b) -> NamedFunc::ScalarType{
    double w_year = 1;
    if(b.SampleTypeString().Contains("-")) w_year = 1;
    else{
      bool isdiph = false;
      bool isQCD2 = false;
      bool isttWToLNu = false;
      bool isDYJetsHT = false;
      bool isSMS = false;
      std::string file = *b.FileNames().begin();
      if (file.find("DiPhotonJetsBox") != file.npos) isdiph = true;
      if (file.find("QCD_Pt-40") != file.npos) isQCD2 = true;
      if (file.find("TTWJetsToLNu") != file.npos) isttWToLNu = true;
      if (file.find("DYJetsToLL_M-50_HT") != file.npos) isDYJetsHT = true;
      if (file.find("TChiHH") != file.npos || file.find("TChiHZ") != file.npos) isSMS = true;
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
      if (isDYJetsHT) w_year *= 6077.22/5398.0; // NNLO/LO ratio for DYJets HT-binned samples
    }
    return w_year;
  });

  const NamedFunc w_lumi_hb("w_lumi_hb", [](const Baby &b) -> NamedFunc::ScalarType{
    double w_year = 1;
    if(b.SampleTypeString().Contains("-")) w_year = 1;
    else{
      bool isdiph = false;
      bool isQCD2 = false;
      bool isttWToLNu = false;
      bool isDYJetsHT = false;
      bool isSMS = false;
      std::string file = *b.FileNames().begin();
      if (file.find("DiPhotonJetsBox") != file.npos) isdiph = true;
      if (file.find("QCD_Pt-40") != file.npos) isQCD2 = true;
      if (file.find("TTWJetsToLNu") != file.npos) isttWToLNu = true;
      if (file.find("DYJetsToLL_M-50_HT") != file.npos) isDYJetsHT = true;
      if (file.find("TChiHH") != file.npos || file.find("TChiHZ") != file.npos) isSMS = true;
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
      if (isDYJetsHT) w_year *= 6077.22/5398.0; // NNLO/LO ratio for DYJets HT-binned samples
      //w_year *= 200/137.61; // scale for 200fb-1
    }
    return w_year;
  });

  const NamedFunc w_lumi_hb_mettrig("w_lumi_hb_mettrig", [](const Baby &b) -> NamedFunc::ScalarType{
    double w_year = 1;
    if(b.SampleTypeString().Contains("-")) w_year = 1;
    else{
      bool isdiph = false;
      bool isQCD2 = false;
      bool isttWToLNu = false;
      bool isDYJetsHT = false;
      bool isSMS = false;
      std::string file = *b.FileNames().begin();
      if (file.find("DiPhotonJetsBox") != file.npos) isdiph = true;
      if (file.find("QCD_Pt-40") != file.npos) isQCD2 = true;
      if (file.find("TTWJetsToLNu") != file.npos) isttWToLNu = true;
      if (file.find("DYJetsToLL_M-50_HT") != file.npos) isDYJetsHT = true;
      if (file.find("TChiHH") != file.npos || file.find("TChiHZ") != file.npos) isSMS = true;
      if (b.SampleTypeString()=="2016") {
          w_year = 19.5;
      } else if (b.SampleTypeString()=="2016APV"){
          w_year = 16.8;
      } else if (b.SampleTypeString()=="2017"){
          w_year = 40.6;
      } else if (b.SampleTypeString()=="2018"){
          w_year = 59.83;
      }
      if (isSMS) w_year = 136.71*vardef::xsec_ratio.GetScalar(b); //since we have only 2016 v7 MC
      w_year = w_year*b.w_lumi();//*b.w_trig();
      if(isdiph) w_year *= b.mc_weight()>0 ? b.mc_weight():-b.mc_weight(); /*w_year = w_year*b.mc_weight();*/ //w_lumi already knows the sign
      if(isQCD2) w_year *= -0.1131001131; //due to it being unable to find cross section in nano2pico
      if (isttWToLNu) w_year *= 0.2163/66680; // remove after 4b background has been processed again
      if (isDYJetsHT) w_year *= 6077.22/5398.0; // NNLO/LO ratio for DYJets HT-binned samples
      //w_year *= 200/137.61; // scale for 200fb-1
    }
    return w_year;
  });


  const NamedFunc w_lumi_gmsb("w_lumi_gmsb", [](const Baby &b) -> NamedFunc::ScalarType{
    double w_year = 1;
    if(b.SampleTypeString().Contains("-")) w_year = 1;
    else{
      bool isdiph = false;
      bool isQCD2 = false;
      bool isttWToLNu = false;
      bool isDYJetsHT = false;
      bool isSMS = false;
      std::string file = *b.FileNames().begin();
      if (file.find("DiPhotonJetsBox") != file.npos) isdiph = true;
      if (file.find("QCD_Pt-40") != file.npos) isQCD2 = true;
      if (file.find("TTWJetsToLNu") != file.npos) isttWToLNu = true;
      if (file.find("DYJetsToLL_M-50_HT") != file.npos) isDYJetsHT = true;
      if (file.find("TChiHH") != file.npos || file.find("TChiHZ") != file.npos) isSMS = true;
      if (b.SampleTypeString()=="2016") {
          w_year = 19.5;
      } else if (b.SampleTypeString()=="2016APV"){
          w_year = 16.8;
      } else if (b.SampleTypeString()=="2017"){
          w_year = 41.48;
      } else if (b.SampleTypeString()=="2018"){
          w_year = 59.83;
      }
      if (isSMS) w_year = 137.61; //since we have only 2016 v7 MC
      w_year = w_year*b.w_lumi();//*b.w_trig();
      if(isdiph) w_year *= b.mc_weight()>0 ? b.mc_weight():-b.mc_weight(); /*w_year = w_year*b.mc_weight();*/ //w_lumi already knows the sign
      if(isQCD2) w_year *= -0.1131001131; //due to it being unable to find cross section in nano2pico
      if (isttWToLNu) w_year *= 0.2163/66680; // remove after 4b background has been processed again
      if (isDYJetsHT) w_year *= 6077.22/5398.0; // NNLO/LO ratio for DYJets HT-binned samples
      //w_year *= 200/137.61; // scale for 200fb-1
    }
    return w_year;
  });

  const NamedFunc w_lumi_gmsb_mettrig("w_lumi_gmsb_mettrig", [](const Baby &b) -> NamedFunc::ScalarType{
    double w_year = 1;
    if(b.SampleTypeString().Contains("-")) w_year = 1;
    else{
      bool isdiph = false;
      bool isQCD2 = false;
      bool isttWToLNu = false;
      bool isDYJetsHT = false;
      bool isSMS = false;
      std::string file = *b.FileNames().begin();
      if (file.find("DiPhotonJetsBox") != file.npos) isdiph = true;
      if (file.find("QCD_Pt-40") != file.npos) isQCD2 = true;
      if (file.find("TTWJetsToLNu") != file.npos) isttWToLNu = true;
      if (file.find("DYJetsToLL_M-50_HT") != file.npos) isDYJetsHT = true;
      if (file.find("TChiHH") != file.npos || file.find("TChiHZ") != file.npos) isSMS = true;
      if (b.SampleTypeString()=="2016") {
          w_year = 19.5;
      } else if (b.SampleTypeString()=="2016APV"){
          w_year = 16.8;
      } else if (b.SampleTypeString()=="2017"){
          w_year = 40.6; // scaling down 2017 luminosity to account for trigger
      } else if (b.SampleTypeString()=="2018"){
          w_year = 59.83;
      }
      if (isSMS) w_year = 136.71; // scaling down total luminosity to account for trigger
      w_year = w_year*b.w_lumi();//*b.w_trig();
      if(isdiph) w_year *= b.mc_weight()>0 ? b.mc_weight():-b.mc_weight(); /*w_year = w_year*b.mc_weight();*/ //w_lumi already knows the sign
      if(isQCD2) w_year *= -0.1131001131; //due to it being unable to find cross section in nano2pico
      if (isttWToLNu) w_year *= 0.2163/66680; // remove after 4b background has been processed again
      if (isDYJetsHT) w_year *= 6077.22/5398.0; // NNLO/LO ratio for DYJets HT-binned samples
      //w_year *= 200/137.61; // scale for 200fb-1
    }
    return w_year;
  });


}

namespace triggers{
  const NamedFunc single_ele_trig("single_ele_trig", [](const Baby &b) -> NamedFunc::ScalarType{
    bool trig = false;
    if (b.SampleTypeString().Contains("2016") || b.SampleTypeString().Contains("2016APV")){
        trig = b.HLT_Ele27_WPTight_Gsf();
    } else if (b.SampleTypeString().Contains("2017")){
        trig = b.HLT_Ele32_WPTight_Gsf_L1DoubleEG();
    } else if (b.SampleTypeString().Contains("2018")){
        trig = b.HLT_Ele32_WPTight_Gsf();
    }
    return trig;
  });

  const NamedFunc single_muon_trig("single_muon_trig", [](const Baby &b) -> NamedFunc::ScalarType{
    bool trig = false;
    if (b.SampleTypeString().Contains("2016") || b.SampleTypeString().Contains("2016APV") || b.SampleTypeString().Contains("2018")){
        trig = b.HLT_IsoMu24();
    } else if (b.SampleTypeString().Contains("2017")){
        trig = b.HLT_IsoMu27();
    }   
    return trig;
  });

  const NamedFunc ptmiss_htmiss_trig("ptmiss_htmiss_trig", [](const Baby &b) -> NamedFunc::ScalarType{
    bool trig = false;
    trig = b.HLT_PFMET120_PFMHT120_IDTight(); 
    return trig;
  });


  const NamedFunc diphoton_trig("diphoton_trig",[](const Baby &b) -> NamedFunc::ScalarType{

      bool trig_year = false;
      if (b.SampleTypeString().Contains("2016")|| b.SampleTypeString().Contains("2016APV")) {
          trig_year = b.HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90();
      } else if (b.SampleTypeString().Contains("2017") || b.SampleTypeString().Contains("2018")){
          trig_year = b.HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90();
      }
      return trig_year;
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
