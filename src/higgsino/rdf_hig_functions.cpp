#include "higgsino/rdf_hig_functions.hpp"
#include "higgsino/apply_trigeffs_v0.hpp"
#include "higgsino/apply_trigeffs2016.hpp"
#include "higgsino/apply_trigeffs2017.hpp"
#include "higgsino/apply_trigeffs2018.hpp"

#include "ROOT/RVec.hxx"
#include "TVector2.h"
#include "TMath.h"
#include "Math/Vector4D.h"

#include "core/utilities.hpp"
#include "core/config_parser.hpp"

namespace RdfHigfuncs{

//misc

const std::vector<std::string> lead_signal_lepton_pt_args =  {"el_sig", "el_pt", "mu_sig", "mu_pt"};
float lead_signal_lepton_pt (RVec<bool> el_sig, RVec<float> el_pt, 
                             RVec<bool> mu_sig, RVec<float> mu_pt) {
  // Search for signal electrons
  float lead_electron_pt = -1;
  for (unsigned iEl = 0; iEl < el_sig.size(); ++iEl) {
    if (el_sig.at(iEl)) {
      lead_electron_pt = el_pt.at(iEl); 
      break;
    }
  }
  // Search for signal muons
  float lead_muon_pt = -1;
  for (unsigned iMu = 0; iMu < mu_sig.size(); ++iMu) {
    if (mu_sig.at(iMu)) {
      lead_muon_pt = mu_pt.at(iMu); 
      break;
    }
  }
  // Warnings
  if (lead_electron_pt==-1 && lead_muon_pt==-1) {
    return -1;
  } else if (lead_electron_pt != -1 && lead_muon_pt != -1) {
    // Return max pt
    return std::max(lead_electron_pt, lead_muon_pt);
  } else if (lead_electron_pt != -1 && lead_muon_pt == -1) { // Electron case
    return lead_electron_pt;
  } else { // Muon case
    return lead_muon_pt;
  }
}

const std::vector<std::string> hig_bcat_args = {"nbl", "nbm", "nbt"};
int hig_bcat(int nbl, int nbm, int nbt) {
  if (nbt < 2) return nbt;
  if (nbt == 2 && nbm == 2) return 2;
  if (nbt >= 2 && nbm == 3 && nbl == 3) return 3;
  if (nbt >= 2 && nbm >= 3 && nbl >= 4) return 4;
  return 0;
}

const std::string hig_cand_am_0 = "hig_cand_am[0]";
const std::string hig_cand_dm_0 = "hig_cand_dm[0]";
const std::string hig_cand_drmax_0 = "hig_cand_drmax[0]";

//filters
 
const std::vector<std::string> pass_ecalnoisejetfilter_args = {"jet_pt","jet_eta","jet_phi","met_phi"};
bool pass_ecalnoisejetfilter(RVec<float> jet_pt, RVec<float> jet_eta, 
                       RVec<float> jet_phi, float met_phi) {
  //check top two highest pt jets in eta 2.4 to 5.0 region, if either has pt>250 and is closely aligned or anti-aligned with MET, then the event is vetoed
  int counter = 0;
  bool r_pass_ecalnoisejet;
  bool goodjet[2] = {true, true};
  double dphi = 0.;
  for (unsigned int jet_idx = 0; jet_idx < jet_pt.size(); jet_idx++) {
    if (counter >= 2) break;
    if (jet_pt.at(jet_idx)>30 && fabs(jet_eta.at(jet_idx))>2.4 && fabs(jet_eta.at(jet_idx))<5.0) {
      dphi = fabs(TVector2::Phi_mpi_pi(jet_phi.at(jet_idx)-met_phi));
      if (jet_pt.at(jet_idx)>250 && (dphi > 2.6 || dphi < 0.1)) goodjet[counter] = false;
      ++counter;
    }
  }
  r_pass_ecalnoisejet = goodjet[0] && goodjet[1];
  return r_pass_ecalnoisejet;
}

const std::vector<std::string> pass_hemveto_args = {"type","run","event",
    "el_pt","el_miniso","el_eta","el_phi","jet_pt","jet_eta","jet_phi","met_phi"};
bool pass_hemveto(int type, int run, long event, RVec<float> el_pt, RVec<float> el_miniso, 
         RVec<float> el_eta, RVec<float> el_phi, RVec<float> jet_pt, RVec<float> jet_eta,
         RVec<float> jet_phi, float met_phi) {
    //only apply for 2018 era C+D and some 2018 MC
    if ((type/1000 == 0) && run < 319077) return true; //data before 2018C
    //for now veto fraction across all years in MC since no way to distinguish
    if ((type/1000 != 0) && (event%4501) >= 1296) return true;
    bool pass_hem = true;
    for (unsigned int el_idx = 0; el_idx < el_pt.size(); el_idx++) {
      if (el_miniso.at(el_idx) < 0.1 && -3.0 < el_eta.at(el_idx) && el_eta.at(el_idx) < -1.4 && -1.57 < el_phi.at(el_idx) && el_phi.at(el_idx) < -0.87) {
        pass_hem = false;
      }
    }
    for (unsigned int jet_idx = 0; jet_idx < jet_pt.size(); jet_idx++) {
      if (jet_pt.at(jet_idx) > 30. && -3.2 < jet_eta.at(jet_idx) && jet_eta.at(jet_idx) < -1.2 && -1.77 < jet_phi.at(jet_idx) && jet_phi.at(jet_idx) < -0.67) {
        double dphi = fabs(TVector2::Phi_mpi_pi(jet_phi.at(jet_idx)-met_phi));
        if (dphi < 0.5) {
          pass_hem = false;
        }
      }
    }
    return pass_hem;
}

const std::vector<std::string> pass_filters_args = {"pass_goodv","pass_hbhe",
    "pass_hbheiso","pass_ecaldeadcell","pass_badpfmu",
    "pass_muon_jet","type","pass_eebadsc","pass_cschalo_tight","pass_low_neutral_jet",
    "pass_htratio_dphi_tight","pass_jets","pass_ecalnoisejetfilter","event",
    "pass_hemveto"};
bool pass_filters(bool pass_goodv, bool pass_hbhe, bool pass_hbheiso, bool pass_ecaldeadcell, 
    bool pass_badpfmu, bool pass_muon_jet, int type, bool pass_eebadsc, 
    bool pass_cschalo_tight, bool pass_low_neutral_jet, bool pass_htratio_dphi_tight,
    bool pass_jets, bool pass_ecalnoisejetfilter, long event, bool pass_hemveto){
  //note this treats all TChiHH as fastsim and all T5HH as fullsim
  if (!pass_goodv || !pass_hbhe || !pass_hbheiso || !pass_ecaldeadcell) return false; 
  if (!pass_badpfmu || !pass_muon_jet) return false;
  if ((type/1000 == 0) && !pass_eebadsc) return false;
  if (!(type/1000==106) && !pass_cschalo_tight) return false;
  if (!pass_low_neutral_jet) return false;
  if (!pass_htratio_dphi_tight) return false;
  if ((type/1000==106) && !pass_jets) return false;
  if (type/1000 == 0) {
    if (((type%100) > 26*2) && !pass_ecalnoisejetfilter) return false; 
  }
  else {
    if ((event%1375 > 363) && !pass_ecalnoisejetfilter) return false; 
  }
  if (!pass_hemveto) return false;
  return true;
}

const std::string final_pass_filters = "pass_filters&&met/mht<2&&met/met_calo<2&&weight<1.5";
const std::string final_ttbar_pass_filters = "pass_filters&&met/met_calo<5&&weight<1.5";
const std::string final_zll_pass_filters = "pass_filters&&met/met_calo<5&&weight<1.5";
const std::string final_qcd_pass_filters = "pass_filters&&met/mht<2&&met/met_calo<2";

//weights except trigger weights

const std::vector<std::string> w_pileup_args = {"type", "npu_tru_mean"};
float w_pileup(int type, int npu_tru_mean){
  if (type/1000 == 0) return 1.0;

  // weights are from get_puweights.py
  static std::vector<double> weights2016 = std::vector<double>({3.661e-01, 8.939e-01, 1.198e+00, 9.627e-01, 1.121e+00, 1.165e+00, 7.956e-01, 4.958e-01, 7.422e-01, 8.789e-01, 9.642e-01, 1.072e+00, 1.125e+00, 1.176e+00, 1.202e+00, 1.208e+00, 1.200e+00, 1.183e+00, 1.144e+00, 1.097e+00, 1.066e+00, 1.051e+00, 1.052e+00, 1.051e+00, 1.050e+00, 1.058e+00, 1.072e+00, 1.083e+00, 1.096e+00, 1.108e+00, 1.095e+00, 1.083e+00, 1.041e+00, 9.858e-01, 9.108e-01, 8.209e-01, 7.168e-01, 6.100e-01, 5.031e-01, 4.048e-01, 3.092e-01, 2.279e-01, 1.637e-01, 1.132e-01, 7.730e-02, 5.092e-02, 3.189e-02, 2.009e-02, 1.226e-02, 7.426e-03, 4.380e-03, 2.608e-03, 1.566e-03, 9.714e-04, 7.292e-04, 6.727e-04, 7.305e-04, 9.488e-04, 1.355e-03, 1.894e-03, 3.082e-03, 4.097e-03, 4.874e-03, 5.256e-03, 5.785e-03, 5.515e-03, 5.000e-03, 4.410e-03, 4.012e-03, 3.548e-03, 3.108e-03, 2.702e-03, 2.337e-03, 2.025e-03, 1.723e-03});
  static std::vector<double> weights2017 = std::vector<double>({1.835e-01, 3.933e+00, 3.471e+00, 2.492e+00, 1.625e+00, 1.515e+00, 1.288e+00, 1.278e+00, 6.158e-01, 1.452e+00, 1.498e+00, 1.487e+00, 1.331e+00, 1.164e+00, 1.079e+00, 1.054e+00, 1.081e+00, 1.129e+00, 1.166e+00, 1.189e+00, 1.213e+00, 1.238e+00, 1.260e+00, 1.271e+00, 1.273e+00, 1.271e+00, 1.271e+00, 1.267e+00, 1.274e+00, 1.252e+00, 1.221e+00, 1.170e+00, 1.109e+00, 1.038e+00, 9.694e-01, 9.119e-01, 8.668e-01, 8.345e-01, 7.884e-01, 7.508e-01, 7.592e-01, 7.932e-01, 8.583e-01, 9.588e-01, 1.094e+00, 1.256e+00, 1.420e+00, 1.496e+00, 1.531e+00, 1.462e+00, 1.337e+00, 1.154e+00, 9.506e-01, 7.497e-01, 5.698e-01, 4.107e-01, 2.900e-01, 1.987e-01, 1.376e-01, 9.664e-02, 6.923e-02, 5.086e-02, 3.844e-02, 2.996e-02, 2.409e-02, 1.712e-02, 1.248e-02, 1.077e-02, 9.597e-03, 8.812e-03, 8.310e-03, 8.018e-03, 7.885e-03, 7.875e-03, 6.337e-03, 5.334e-03, 5.444e-03, 5.584e-03, 5.744e-03, 5.916e-03, 6.090e-03, 6.257e-03, 6.408e-03, 5.095e-03, 4.228e-03, 4.259e-03, 4.266e-03, 4.247e-03, 4.203e-03, 4.137e-03, 4.052e-03, 3.950e-03, 2.943e-03, 2.296e-03, 2.200e-03, 2.101e-03, 2.002e-03, 1.902e-03, 1.804e-03});
  static std::vector<double> weights2018 = std::vector<double>({9.687e+03, 1.327e+01, 4.377e+01, 1.869e+01, 1.251e+01, 9.042e+00, 6.579e+00, 4.873e+00, 3.629e+00, 2.763e+00, 2.224e+00, 1.900e+00, 1.702e+00, 1.579e+00, 1.504e+00, 1.465e+00, 1.450e+00, 1.451e+00, 1.459e+00, 1.464e+00, 1.458e+00, 1.437e+00, 1.401e+00, 1.353e+00, 1.301e+00, 1.250e+00, 1.205e+00, 1.167e+00, 1.137e+00, 1.115e+00, 1.100e+00, 1.089e+00, 1.082e+00, 1.077e+00, 1.073e+00, 1.069e+00, 1.063e+00, 1.052e+00, 1.036e+00, 1.012e+00, 9.794e-01, 9.373e-01, 8.862e-01, 8.272e-01, 7.618e-01, 6.923e-01, 6.211e-01, 5.505e-01, 4.827e-01, 4.194e-01, 3.618e-01, 3.105e-01, 2.659e-01, 2.276e-01, 1.952e-01, 1.679e-01, 1.450e-01, 1.259e-01, 1.097e-01, 9.590e-02, 8.401e-02, 7.361e-02, 6.442e-02, 5.623e-02, 4.889e-02, 4.231e-02, 3.643e-02, 3.122e-02, 2.662e-02, 2.262e-02, 1.915e-02, 1.618e-02, 1.363e-02, 1.147e-02, 9.633e-03, 8.075e-03, 6.750e-03, 5.623e-03, 4.661e-03, 3.838e-03, 3.134e-03, 2.531e-03, 2.016e-03, 1.578e-03, 1.208e-03, 8.991e-04, 6.465e-04, 4.458e-04, 2.929e-04, 1.824e-04, 1.094e-04, 6.258e-05, 3.428e-05, 1.807e-05, 9.221e-06, 4.578e-06, 2.221e-06, 1.058e-06, 4.964e-07, 2.308e-07});
  unsigned npu_tru_mean_uns = static_cast<unsigned>(npu_tru_mean);
  unsigned wgt2016_idx = npu_tru_mean_uns;
  unsigned wgt2017_idx = npu_tru_mean_uns;
  unsigned wgt2018_idx = npu_tru_mean_uns;
  if (npu_tru_mean_uns>=weights2016.size()) wgt2016_idx = weights2016.size()-1;
  if (npu_tru_mean_uns>=weights2017.size()) wgt2017_idx = weights2016.size()-1;
  if (npu_tru_mean_uns>=weights2018.size()) wgt2018_idx = weights2016.size()-1;
  
  return (36.3*weights2016.at(wgt2016_idx)+
      41.5*weights2017.at(wgt2017_idx)+59.7*weights2018.at(wgt2018_idx))/137.5;
}

const std::vector<std::string> w_run2_args = {"type"};
float w_run2(int type) {
  if (type/1000==0) return 1.0;
  return (36.32264+41.52756+59.67377);
}

const std::string final_weight_run2 = "weight*eff_higtrig_run2*w_run2*w_pileup";
const std::string final_weight_run2_data = "weight*eff_higtrig_run2*w_run2";

//triggers

const std::string jet_trigger = "HLT_PFJet500";
const std::string met_trigger = "HLT_PFMET90_PFMHT90_IDTight||HLT_PFMETNoMu90_PFMHTNoMu90_IDTight||HLT_PFMET100_PFMHT100_IDTight||HLT_PFMETNoMu100_PFMHTNoMu100_IDTight||HLT_PFMET110_PFMHT110_IDTight||HLT_PFMETNoMu110_PFMHTNoMu110_IDTight||HLT_PFMET120_PFMHT120_IDTight||HLT_PFMETNoMu120_PFMHTNoMu120_IDTight||HLT_PFMET130_PFMHT130_IDTight||HLT_PFMETNoMu130_PFMHTNoMu130_IDTight||HLT_PFMET140_PFMHT140_IDTight||HLT_PFMETNoMu140_PFMHTNoMu140_IDTight||HLT_PFMET100_PFMHT100_IDTight_PFHT60||HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60||HLT_PFMET110_PFMHT110_IDTight_PFHT60||HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_PFHT60||HLT_PFMET120_PFMHT120_IDTight_PFHT60||HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60||HLT_PFMET130_PFMHT130_IDTight_PFHT60||HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_PFHT60||HLT_PFMET140_PFMHT140_IDTight_PFHT60||HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_PFHT60||HLT_PFMET120_PFMHT120_IDTight_HFCleaned||HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned||HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned";
const std::string el_trigger =  "HLT_Ele25_WPTight_Gsf||HLT_Ele27_WPTight_Gsf||HLT_Ele28_WPTight_Gsf||HLT_Ele32_WPTight_Gsf_L1DoubleEG||HLT_Ele32_WPTight_Gsf||HLT_Ele35_WPTight_Gsf||HLT_Ele45_WPLoose_Gsf||HLT_Ele105_CaloIdVT_GsfTrkIdT||HLT_Ele115_CaloIdVT_GsfTrkIdT||HLT_Ele135_CaloIdVT_GsfTrkIdT||HLT_Ele145_CaloIdVT_GsfTrkIdT||HLT_Ele25_eta2p1_WPTight_Gsf||HLT_Ele27_eta2p1_WPTight_Gsf||HLT_Ele27_eta2p1_WPLoose_Gsf||HLT_Ele20_WPLoose_Gsf||HLT_Ele20_eta2p1_WPLoose_Gsf||HLT_Ele25_eta2p1_WPLoose_Gsf||HLT_Ele15_IsoVVVL_PFHT350||HLT_Ele15_IsoVVVL_PFHT400||HLT_Ele15_IsoVVVL_PFHT450||HLT_Ele15_IsoVVVL_PFHT600||HLT_Ele50_IsoVVVL_PFHT450";
const std::string mu_trigger = "HLT_IsoMu20||HLT_IsoMu22||HLT_IsoMu24||HLT_IsoMu27||HLT_IsoTkMu20||HLT_IsoTkMu22||HLT_IsoTkMu24||HLT_Mu50||HLT_Mu55||HLT_TkMu50||HLT_IsoMu22_eta2p1||HLT_IsoMu24_eta2p1||HLT_Mu45_eta2p1||HLT_Mu15_IsoVVVL_PFHT350||HLT_Mu15_IsoVVVL_PFHT400||HLT_Mu15_IsoVVVL_PFHT450||HLT_Mu15_IsoVVVL_PFHT600||HLT_Mu50_IsoVVVL_PFHT400||HLT_Mu50_IsoVVVL_PFHT450";
const std::string ht_trigger = "HLT_PFHT125||HLT_PFHT200||HLT_PFHT300||HLT_PFHT400||HLT_PFHT475||HLT_PFHT600||HLT_PFHT650||HLT_PFHT800||HLT_PFHT900||HLT_PFHT180||HLT_PFHT370||HLT_PFHT430||HLT_PFHT510||HLT_PFHT590||HLT_PFHT680||HLT_PFHT780||HLT_PFHT890||HLT_PFHT1050||HLT_PFHT250||HLT_PFHT350";
const std::string jetht_trigger = "jet_trigger||ht_trigger";

//function to add columns to RDataFrame
ProcessedDataFrame rdf_hig_functions(ProcessedDataFrame &data_frame) {
  data_frame = data_frame.Define("lead_signal_lepton_pt",lead_signal_lepton_pt,
      lead_signal_lepton_pt_args);
  data_frame = data_frame.Define("hig_bcat",hig_bcat,hig_bcat_args);
  data_frame = data_frame.Define("hig_cand_am_0",hig_cand_am_0);
  data_frame = data_frame.Define("hig_cand_dm_0",hig_cand_dm_0);
  data_frame = data_frame.Define("hig_cand_drmax_0",hig_cand_drmax_0);
  data_frame = data_frame.Define("pass_hemveto",pass_hemveto,pass_hemveto_args);
  data_frame = data_frame.Define("pass_ecalnoisejetfilter",pass_ecalnoisejetfilter,pass_ecalnoisejetfilter_args);
  data_frame = data_frame.Define("pass_filters",pass_filters,pass_filters_args);
  data_frame = data_frame.Define("final_pass_filters",final_pass_filters);
  data_frame = data_frame.Define("final_ttbar_pass_filters",final_ttbar_pass_filters);
  data_frame = data_frame.Define("final_zll_pass_filters",final_zll_pass_filters);
  data_frame = data_frame.Define("final_qcd_pass_filters",final_qcd_pass_filters);
  data_frame = data_frame.Define("w_pileup",w_pileup,w_pileup_args);
  data_frame = data_frame.Define("w_run2",w_run2,w_run2_args);
  data_frame = data_frame.Define("eff_higtrig_run2",eff_higtrig_run2,eff_higtrig_run2_args);
  data_frame = data_frame.Define("final_weight_run2",final_weight_run2);
  data_frame = data_frame.Define("jet_trigger",jet_trigger);
  data_frame = data_frame.Define("met_trigger",met_trigger);
  data_frame = data_frame.Define("el_trigger",el_trigger);
  data_frame = data_frame.Define("mu_trigger",mu_trigger);
  data_frame = data_frame.Define("ht_trigger",ht_trigger);
  data_frame = data_frame.Define("jetht_trigger",jetht_trigger);
  return data_frame;
}

//function to add columns to RDataFrame
//note that slim doesn't have the other trigger branches
ProcessedDataFrame rdf_hig_functions_slim_mc(ProcessedDataFrame &data_frame) {
  data_frame = data_frame.Define("lead_signal_lepton_pt",lead_signal_lepton_pt,
      lead_signal_lepton_pt_args);
  data_frame = data_frame.Define("hig_bcat",hig_bcat,hig_bcat_args);
  data_frame = data_frame.Define("hig_cand_am_0",hig_cand_am_0);
  data_frame = data_frame.Define("hig_cand_dm_0",hig_cand_dm_0);
  data_frame = data_frame.Define("hig_cand_drmax_0",hig_cand_drmax_0);
  data_frame = data_frame.Define("pass_hemveto",pass_hemveto,pass_hemveto_args);
  data_frame = data_frame.Define("pass_ecalnoisejetfilter",pass_ecalnoisejetfilter,pass_ecalnoisejetfilter_args);
  data_frame = data_frame.Define("pass_filters",pass_filters,pass_filters_args);
  data_frame = data_frame.Define("final_pass_filters",final_pass_filters);
  data_frame = data_frame.Define("final_ttbar_pass_filters",final_ttbar_pass_filters);
  data_frame = data_frame.Define("final_zll_pass_filters",final_zll_pass_filters);
  data_frame = data_frame.Define("final_qcd_pass_filters",final_qcd_pass_filters);
  data_frame = data_frame.Define("w_pileup",w_pileup,w_pileup_args);
  data_frame = data_frame.Define("w_run2",w_run2,w_run2_args);
  data_frame = data_frame.Define("eff_higtrig_run2",eff_higtrig_run2,eff_higtrig_run2_args);
  data_frame = data_frame.Define("final_weight_run2",final_weight_run2);
  data_frame = data_frame.Define("met_trigger",met_trigger);
  return data_frame;
}

ProcessedDataFrame rdf_hig_functions_slim_data(ProcessedDataFrame &data_frame) {
  data_frame = data_frame.Define("lead_signal_lepton_pt",lead_signal_lepton_pt,
      lead_signal_lepton_pt_args);
  data_frame = data_frame.Define("hig_bcat",hig_bcat,hig_bcat_args);
  data_frame = data_frame.Define("hig_cand_am_0",hig_cand_am_0);
  data_frame = data_frame.Define("hig_cand_dm_0",hig_cand_dm_0);
  data_frame = data_frame.Define("hig_cand_drmax_0",hig_cand_drmax_0);
  data_frame = data_frame.Define("pass_hemveto",pass_hemveto,pass_hemveto_args);
  data_frame = data_frame.Define("pass_ecalnoisejetfilter",pass_ecalnoisejetfilter,pass_ecalnoisejetfilter_args);
  data_frame = data_frame.Define("pass_filters",pass_filters,pass_filters_args);
  data_frame = data_frame.Define("final_pass_filters",final_pass_filters);
  data_frame = data_frame.Define("final_ttbar_pass_filters",final_ttbar_pass_filters);
  data_frame = data_frame.Define("final_zll_pass_filters",final_zll_pass_filters);
  data_frame = data_frame.Define("final_qcd_pass_filters",final_qcd_pass_filters);
  data_frame = data_frame.Define("w_run2",w_run2,w_run2_args);
  data_frame = data_frame.Define("eff_higtrig_run2",eff_higtrig_run2,eff_higtrig_run2_args);
  data_frame = data_frame.Define("final_weight_run2",final_weight_run2_data);
  data_frame = data_frame.Define("met_trigger",met_trigger);
  return data_frame;
}

//trigger weights

const std::vector<std::string> eff_higtrig_run2_args = {"type","nvlep","nel",
    "nmu","met","ht","lead_signal_lepton_pt"};
float eff_higtrig_run2(int type, int nvlep, int nel, int nmu, float met, 
                       float ht, float lead_signal_lepton_pt){
  if (type/1000 == 0) return 1.0; //data
  float errup=0., errdown=0.; // Not used, but for reference
  errup += errdown; //suppresses unused warning
  float eff2016=1.0,eff2017=1.0,eff2018=1.0;
  float lep_pt = lead_signal_lepton_pt;
  if (nvlep==0) { //search/QCD CR
    if(type/1000 != 7) { // REAL MET MC
      if (ht> 0 && ht<= 350 && met> 150 && met<= 160) {eff2016 = 0.404605; errup = 0.0281502; errdown = 0.0281502;}
      else if (ht> 0 && ht<= 350 && met> 160 && met<= 170) {eff2016 = 0.571429; errup = 0.0374088; errdown = 0.0374088;}
      else if (ht> 0 && ht<= 350 && met> 170 && met<= 180) {eff2016 = 0.613636; errup = 0.0423806; errdown = 0.0423806;}
      else if (ht> 0 && ht<= 350 && met> 180 && met<= 190) {eff2016 = 0.716418; errup = 0.0550662; errdown = 0.0550662;}
      else if (ht> 0 && ht<= 350 && met> 190 && met<= 200) {eff2016 = 0.666667; errup = 0.0727393; errdown = 0.0727393;}
      else if (ht> 0 && ht<= 350 && met> 200 && met<= 225) {eff2016 = 0.656716; errup = 0.0580067; errdown = 0.0580067;}
      else if (ht> 0 && ht<= 350 && met> 225 && met<= 250) {eff2016 = 0.846154; errup = 0.0707589; errdown = 0.0707589;}
      else if (ht> 0 && ht<= 350 && met> 250 && met<= 9999) {eff2016 = 0.73913; errup = 0.0915605; errdown = 0.0915605;}
      else if (ht> 350 && ht<= 450 && met> 150 && met<= 155) {eff2016 = 0.457086; errup = 0.0222559; errdown = 0.0222559;}
      else if (ht> 350 && ht<= 450 && met> 155 && met<= 160) {eff2016 = 0.509589; errup = 0.0261664; errdown = 0.0261664;}
      else if (ht> 350 && ht<= 450 && met> 160 && met<= 165) {eff2016 = 0.52795; errup = 0.0278203; errdown = 0.0278203;}
      else if (ht> 350 && ht<= 450 && met> 165 && met<= 170) {eff2016 = 0.60084; errup = 0.0317442; errdown = 0.0317442;}
      else if (ht> 350 && ht<= 450 && met> 170 && met<= 175) {eff2016 = 0.636364; errup = 0.0332746; errdown = 0.0332746;}
      else if (ht> 350 && ht<= 450 && met> 175 && met<= 180) {eff2016 = 0.517647; errup = 0.0383244; errdown = 0.0383244;}
      else if (ht> 350 && ht<= 450 && met> 180 && met<= 185) {eff2016 = 0.661417; errup = 0.0419922; errdown = 0.0419922;}
      else if (ht> 350 && ht<= 450 && met> 185 && met<= 190) {eff2016 = 0.694444; errup = 0.0443253; errdown = 0.0443253;}
      else if (ht> 350 && ht<= 450 && met> 190 && met<= 195) {eff2016 = 0.621053; errup = 0.0497728; errdown = 0.0497728;}
      else if (ht> 350 && ht<= 450 && met> 195 && met<= 200) {eff2016 = 0.659091; errup = 0.0505302; errdown = 0.0505302;}
      else if (ht> 350 && ht<= 450 && met> 200 && met<= 210) {eff2016 = 0.778626; errup = 0.0362737; errdown = 0.0362737;}
      else if (ht> 350 && ht<= 450 && met> 210 && met<= 220) {eff2016 = 0.738636; errup = 0.0468378; errdown = 0.0468378;}
      else if (ht> 350 && ht<= 450 && met> 220 && met<= 230) {eff2016 = 0.746269; errup = 0.0531615; errdown = 0.0531615;}
      else if (ht> 350 && ht<= 450 && met> 230 && met<= 240) {eff2016 = 0.758621; errup = 0.0561886; errdown = 0.0561886;}
      else if (ht> 350 && ht<= 450 && met> 240 && met<= 250) {eff2016 = 0.744186; errup = 0.0665378; errdown = 0.0665378;}
      else if (ht> 350 && ht<= 450 && met> 250 && met<= 300) {eff2016 = 0.754098; errup = 0.0389865; errdown = 0.0389865;}
      else if (ht> 350 && ht<= 450 && met> 300 && met<= 9999) {eff2016 = 0.755556; errup = 0.0640644; errdown = 0.0640644;}
      else if (ht> 450 && ht<= 550 && met> 150 && met<= 155) {eff2016 = 0.432544; errup = 0.0158991; errdown = 0.0158991;}
      else if (ht> 450 && ht<= 550 && met> 155 && met<= 160) {eff2016 = 0.467866; errup = 0.0178888; errdown = 0.0178888;}
      else if (ht> 450 && ht<= 550 && met> 160 && met<= 165) {eff2016 = 0.510542; errup = 0.0193994; errdown = 0.0193994;}
      else if (ht> 450 && ht<= 550 && met> 165 && met<= 170) {eff2016 = 0.525773; errup = 0.0206981; errdown = 0.0206981;}
      else if (ht> 450 && ht<= 550 && met> 170 && met<= 175) {eff2016 = 0.569476; errup = 0.0236322; errdown = 0.0236322;}
      else if (ht> 450 && ht<= 550 && met> 175 && met<= 180) {eff2016 = 0.604369; errup = 0.0240906; errdown = 0.0240906;}
      else if (ht> 450 && ht<= 550 && met> 180 && met<= 185) {eff2016 = 0.629747; errup = 0.0271637; errdown = 0.0271637;}
      else if (ht> 450 && ht<= 550 && met> 185 && met<= 190) {eff2016 = 0.620567; errup = 0.028896; errdown = 0.028896;}
      else if (ht> 450 && ht<= 550 && met> 190 && met<= 195) {eff2016 = 0.626609; errup = 0.0316885; errdown = 0.0316885;}
      else if (ht> 450 && ht<= 550 && met> 195 && met<= 200) {eff2016 = 0.633663; errup = 0.0338995; errdown = 0.0338995;}
      else if (ht> 450 && ht<= 550 && met> 200 && met<= 210) {eff2016 = 0.680645; errup = 0.0264799; errdown = 0.0264799;}
      else if (ht> 450 && ht<= 550 && met> 210 && met<= 220) {eff2016 = 0.744589; errup = 0.0286928; errdown = 0.0286928;}
      else if (ht> 450 && ht<= 550 && met> 220 && met<= 230) {eff2016 = 0.658163; errup = 0.0338804; errdown = 0.0338804;}
      else if (ht> 450 && ht<= 550 && met> 230 && met<= 240) {eff2016 = 0.706667; errup = 0.0371743; errdown = 0.0371743;}
      else if (ht> 450 && ht<= 550 && met> 240 && met<= 250) {eff2016 = 0.700855; errup = 0.0423314; errdown = 0.0423314;}
      else if (ht> 450 && ht<= 550 && met> 250 && met<= 300) {eff2016 = 0.642173; errup = 0.0270951; errdown = 0.0270951;}
      else if (ht> 450 && ht<= 550 && met> 300 && met<= 400) {eff2016 = 0.609626; errup = 0.035674; errdown = 0.035674;}
      else if (ht> 450 && ht<= 550 && met> 400 && met<= 9999) {eff2016 = 0.926471; errup = 0.0316513; errdown = 0.0316513;}
      else if (ht> 550 && ht<= 650 && met> 150 && met<= 155) {eff2016 = 0.404868; errup = 0.0109406; errdown = 0.0109406;}
      else if (ht> 550 && ht<= 650 && met> 155 && met<= 160) {eff2016 = 0.437424; errup = 0.0122272; errdown = 0.0122272;}
      else if (ht> 550 && ht<= 650 && met> 160 && met<= 165) {eff2016 = 0.488841; errup = 0.0134125; errdown = 0.0134125;}
      else if (ht> 550 && ht<= 650 && met> 165 && met<= 170) {eff2016 = 0.510135; errup = 0.014528; errdown = 0.014528;}
      else if (ht> 550 && ht<= 650 && met> 170 && met<= 175) {eff2016 = 0.534884; errup = 0.0158603; errdown = 0.0158603;}
      else if (ht> 550 && ht<= 650 && met> 175 && met<= 180) {eff2016 = 0.554398; errup = 0.0169094; errdown = 0.0169094;}
      else if (ht> 550 && ht<= 650 && met> 180 && met<= 185) {eff2016 = 0.565807; errup = 0.0182575; errdown = 0.0182575;}
      else if (ht> 550 && ht<= 650 && met> 185 && met<= 190) {eff2016 = 0.625; errup = 0.0189018; errdown = 0.0189018;}
      else if (ht> 550 && ht<= 650 && met> 190 && met<= 195) {eff2016 = 0.590909; errup = 0.0209647; errdown = 0.0209647;}
      else if (ht> 550 && ht<= 650 && met> 195 && met<= 200) {eff2016 = 0.655031; errup = 0.0215405; errdown = 0.0215405;}
      else if (ht> 550 && ht<= 650 && met> 200 && met<= 210) {eff2016 = 0.616162; errup = 0.0172806; errdown = 0.0172806;}
      else if (ht> 550 && ht<= 650 && met> 210 && met<= 220) {eff2016 = 0.635556; errup = 0.0185242; errdown = 0.0185242;}
      else if (ht> 550 && ht<= 650 && met> 220 && met<= 230) {eff2016 = 0.642147; errup = 0.021374; errdown = 0.021374;}
      else if (ht> 550 && ht<= 650 && met> 230 && met<= 240) {eff2016 = 0.663529; errup = 0.0229197; errdown = 0.0229197;}
      else if (ht> 550 && ht<= 650 && met> 240 && met<= 250) {eff2016 = 0.617391; errup = 0.0261666; errdown = 0.0261666;}
      else if (ht> 550 && ht<= 650 && met> 250 && met<= 275) {eff2016 = 0.606112; errup = 0.0201328; errdown = 0.0201328;}
      else if (ht> 550 && ht<= 650 && met> 275 && met<= 300) {eff2016 = 0.602532; errup = 0.0246231; errdown = 0.0246231;}
      else if (ht> 550 && ht<= 650 && met> 300 && met<= 400) {eff2016 = 0.704819; errup = 0.017701; errdown = 0.017701;}
      else if (ht> 550 && ht<= 650 && met> 400 && met<= 9999) {eff2016 = 0.955; errup = 0.00732931; errdown = 0.00732931;}
      else if (ht> 650 && ht<= 800 && met> 150 && met<= 155) {eff2016 = 0.403256; errup = 0.00518037; errdown = 0.00518037;}
      else if (ht> 650 && ht<= 800 && met> 155 && met<= 160) {eff2016 = 0.433503; errup = 0.00573216; errdown = 0.00573216;}
      else if (ht> 650 && ht<= 800 && met> 160 && met<= 165) {eff2016 = 0.476213; errup = 0.00629978; errdown = 0.00629978;}
      else if (ht> 650 && ht<= 800 && met> 165 && met<= 170) {eff2016 = 0.504982; errup = 0.00679124; errdown = 0.00679124;}
      else if (ht> 650 && ht<= 800 && met> 170 && met<= 175) {eff2016 = 0.529961; errup = 0.00735406; errdown = 0.00735406;}
      else if (ht> 650 && ht<= 800 && met> 175 && met<= 180) {eff2016 = 0.56639; errup = 0.00786624; errdown = 0.00786624;}
      else if (ht> 650 && ht<= 800 && met> 180 && met<= 185) {eff2016 = 0.577905; errup = 0.00836627; errdown = 0.00836627;}
      else if (ht> 650 && ht<= 800 && met> 185 && met<= 190) {eff2016 = 0.600321; errup = 0.00877644; errdown = 0.00877644;}
      else if (ht> 650 && ht<= 800 && met> 190 && met<= 195) {eff2016 = 0.619677; errup = 0.00930155; errdown = 0.00930155;}
      else if (ht> 650 && ht<= 800 && met> 195 && met<= 200) {eff2016 = 0.623563; errup = 0.0098163; errdown = 0.0098163;}
      else if (ht> 650 && ht<= 800 && met> 200 && met<= 210) {eff2016 = 0.66071; errup = 0.0074079; errdown = 0.0074079;}
      else if (ht> 650 && ht<= 800 && met> 210 && met<= 220) {eff2016 = 0.683814; errup = 0.00823533; errdown = 0.00823533;}
      else if (ht> 650 && ht<= 800 && met> 220 && met<= 230) {eff2016 = 0.702176; errup = 0.00893583; errdown = 0.00893583;}
      else if (ht> 650 && ht<= 800 && met> 230 && met<= 240) {eff2016 = 0.724342; errup = 0.00952245; errdown = 0.00952245;}
      else if (ht> 650 && ht<= 800 && met> 240 && met<= 250) {eff2016 = 0.745802; errup = 0.0104774; errdown = 0.0104774;}
      else if (ht> 650 && ht<= 800 && met> 250 && met<= 275) {eff2016 = 0.775077; errup = 0.00732399; errdown = 0.00732399;}
      else if (ht> 650 && ht<= 800 && met> 275 && met<= 300) {eff2016 = 0.808443; errup = 0.00838429; errdown = 0.00838429;}
      else if (ht> 650 && ht<= 800 && met> 300 && met<= 350) {eff2016 = 0.844642; errup = 0.00697529; errdown = 0.00697529;}
      else if (ht> 650 && ht<= 800 && met> 350 && met<= 400) {eff2016 = 0.904545; errup = 0.00748778; errdown = 0.00748778;}
      else if (ht> 650 && ht<= 800 && met> 400 && met<= 450) {eff2016 = 0.953261; errup = 0.00695909; errdown = 0.00695909;}
      else if (ht> 650 && ht<= 800 && met> 450 && met<= 500) {eff2016 = 0.969849; errup = 0.00699865; errdown = 0.00699865;}
      else if (ht> 650 && ht<= 800 && met> 500 && met<= 9999) {eff2016 = 0.990217; errup = 0.00324488; errdown = 0.00324488;}
      else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff2016 = 0.393245; errup = 0.00244599; errdown = 0.00244599;}
      else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff2016 = 0.43438; errup = 0.00271824; errdown = 0.00271824;}
      else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff2016 = 0.480271; errup = 0.00299093; errdown = 0.00299093;}
      else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff2016 = 0.51571; errup = 0.00327196; errdown = 0.00327196;}
      else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff2016 = 0.553577; errup = 0.00351606; errdown = 0.00351606;}
      else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff2016 = 0.592697; errup = 0.00377356; errdown = 0.00377356;}
      else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff2016 = 0.619873; errup = 0.00400325; errdown = 0.00400325;}
      else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff2016 = 0.651545; errup = 0.00424702; errdown = 0.00424702;}
      else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff2016 = 0.682976; errup = 0.00444573; errdown = 0.00444573;}
      else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff2016 = 0.704728; errup = 0.00462451; errdown = 0.00462451;}
      else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff2016 = 0.734168; errup = 0.00347758; errdown = 0.00347758;}
      else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff2016 = 0.777321; errup = 0.00366278; errdown = 0.00366278;}
      else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff2016 = 0.803227; errup = 0.00394339; errdown = 0.00394339;}
      else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff2016 = 0.827035; errup = 0.00416252; errdown = 0.00416252;}
      else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff2016 = 0.84477; errup = 0.00445204; errdown = 0.00445204;}
      else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff2016 = 0.859084; errup = 0.00316684; errdown = 0.00316684;}
      else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff2016 = 0.890933; errup = 0.00359947; errdown = 0.00359947;}
      else if (ht> 800 && ht<= 1000 && met> 300 && met<= 350) {eff2016 = 0.905381; errup = 0.00318308; errdown = 0.00318308;}
      else if (ht> 800 && ht<= 1000 && met> 350 && met<= 400) {eff2016 = 0.927529; errup = 0.00414893; errdown = 0.00414893;}
      else if (ht> 800 && ht<= 1000 && met> 400 && met<= 450) {eff2016 = 0.947605; errup = 0.00497749; errdown = 0.00497749;}
      else if (ht> 800 && ht<= 1000 && met> 450 && met<= 500) {eff2016 = 0.967402; errup = 0.00549868; errdown = 0.00549868;}
      else if (ht> 800 && ht<= 1000 && met> 500 && met<= 9999) {eff2016 = 0.985313; errup = 0.0030399; errdown = 0.0030399;}
      else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff2016 = 0.294527; errup = 0.00165956; errdown = 0.00165956;}
      else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff2016 = 0.323816; errup = 0.00185338; errdown = 0.00185338;}
      else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff2016 = 0.352023; errup = 0.00204674; errdown = 0.00204674;}
      else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff2016 = 0.378594; errup = 0.00225679; errdown = 0.00225679;}
      else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff2016 = 0.402409; errup = 0.00245637; errdown = 0.00245637;}
      else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff2016 = 0.433744; errup = 0.00267275; errdown = 0.00267275;}
      else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff2016 = 0.463049; errup = 0.00289395; errdown = 0.00289395;}
      else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff2016 = 0.48537; errup = 0.00310318; errdown = 0.00310318;}
      else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff2016 = 0.513993; errup = 0.00332597; errdown = 0.00332597;}
      else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff2016 = 0.533194; errup = 0.00356519; errdown = 0.00356519;}
      else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff2016 = 0.560839; errup = 0.00275238; errdown = 0.00275238;}
      else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff2016 = 0.599586; errup = 0.00306054; errdown = 0.00306054;}
      else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff2016 = 0.629718; errup = 0.00337431; errdown = 0.00337431;}
      else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff2016 = 0.669613; errup = 0.00367228; errdown = 0.00367228;}
      else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff2016 = 0.69982; errup = 0.00397398; errdown = 0.00397398;}
      else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff2016 = 0.733677; errup = 0.00285148; errdown = 0.00285148;}
      else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff2016 = 0.775371; errup = 0.00333625; errdown = 0.00333625;}
      else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 350) {eff2016 = 0.805159; errup = 0.00296283; errdown = 0.00296283;}
      else if (ht> 1000 && ht<= 9999 && met> 350 && met<= 400) {eff2016 = 0.847872; errup = 0.00376673; errdown = 0.00376673;}
      else if (ht> 1000 && ht<= 9999 && met> 400 && met<= 450) {eff2016 = 0.875209; errup = 0.00478206; errdown = 0.00478206;}
      else if (ht> 1000 && ht<= 9999 && met> 450 && met<= 500) {eff2016 = 0.888848; errup = 0.00603014; errdown = 0.00603014;}
      else if (ht> 1000 && ht<= 9999 && met> 500 && met<= 9999) {eff2016 = 0.92654; errup = 0.00382917; errdown = 0.00382917;}
      //2017
      if (ht> 0 && ht<= 200 && met> 150 && met<= 155) {eff2017 = 0.123457; errup = 0.0511779; errdown = 0.0511779;}
      else if (ht> 0 && ht<= 200 && met> 155 && met<= 160) {eff2017 = 0.151724; errup = 0.0624119; errdown = 0.0624119;}
      else if (ht> 0 && ht<= 200 && met> 160 && met<= 165) {eff2017 = 0.228311; errup = 0.0857248; errdown = 0.0857248;}
      else if (ht> 0 && ht<= 200 && met> 165 && met<= 170) {eff2017 = 0.228395; errup = 0.0873896; errdown = 0.0873896;}
      else if (ht> 0 && ht<= 200 && met> 170 && met<= 180) {eff2017 = 0.251366; errup = 0.0946626; errdown = 0.0946626;}
      else if (ht> 0 && ht<= 200 && met> 180 && met<= 190) {eff2017 = 0.285714; errup = 0.103601; errdown = 0.103601;}
      else if (ht> 0 && ht<= 200 && met> 190 && met<= 200) {eff2017 = 0.333333; errup = 0.127092; errdown = 0.127092;}
      else if (ht> 0 && ht<= 200 && met> 200 && met<= 9999) {eff2017 = 0.574468; errup = 0.130371; errdown = 0.130371;}
      else if (ht> 200 && ht<= 300 && met> 150 && met<= 155) {eff2017 = 0.297837; errup = 0.116824; errdown = 0.116824;}
      else if (ht> 200 && ht<= 300 && met> 155 && met<= 160) {eff2017 = 0.347426; errup = 0.136067; errdown = 0.136067;}
      else if (ht> 200 && ht<= 300 && met> 160 && met<= 165) {eff2017 = 0.385093; errup = 0.138234; errdown = 0.138234;}
      else if (ht> 200 && ht<= 300 && met> 165 && met<= 170) {eff2017 = 0.454955; errup = 0.162926; errdown = 0.162926;}
      else if (ht> 200 && ht<= 300 && met> 170 && met<= 175) {eff2017 = 0.497354; errup = 0.178092; errdown = 0.178092;}
      else if (ht> 200 && ht<= 300 && met> 175 && met<= 180) {eff2017 = 0.582192; errup = 0.208295; errdown = 0.208295;}
      else if (ht> 200 && ht<= 300 && met> 180 && met<= 185) {eff2017 = 0.656934; errup = 0.223291; errdown = 0.223291;}
      else if (ht> 200 && ht<= 300 && met> 185 && met<= 190) {eff2017 = 0.663551; errup = 0.225992; errdown = 0.225992;}
      else if (ht> 200 && ht<= 300 && met> 190 && met<= 195) {eff2017 = 0.716495; errup = 0.243676; errdown = 0.243676;}
      else if (ht> 200 && ht<= 300 && met> 195 && met<= 200) {eff2017 = 0.789189; errup = 0.210811; errdown = 0.267707;}
      else if (ht> 200 && ht<= 300 && met> 200 && met<= 210) {eff2017 = 0.857143; errup = 0.142857; errdown = 0.180304;}
      else if (ht> 200 && ht<= 300 && met> 210 && met<= 220) {eff2017 = 0.895028; errup = 0.104972; errdown = 0.188319;}
      else if (ht> 200 && ht<= 300 && met> 220 && met<= 230) {eff2017 = 0.89726; errup = 0.10274; errdown = 0.189079;}
      else if (ht> 200 && ht<= 300 && met> 230 && met<= 240) {eff2017 = 0.932039; errup = 0.0679612; errdown = 0.0964633;}
      else if (ht> 200 && ht<= 300 && met> 240 && met<= 250) {eff2017 = 0.914286; errup = 0.0857143; errdown = 0.0973746;}
      else if (ht> 200 && ht<= 300 && met> 250 && met<= 275) {eff2017 = 0.9375; errup = 0.0331599; errdown = 0.0331599;}
      else if (ht> 200 && ht<= 300 && met> 275 && met<= 9999) {eff2017 = 0.977778; errup = 0.0222222; errdown = 0.0297024;}
      else if (ht> 300 && ht<= 400 && met> 150 && met<= 155) {eff2017 = 0.419244; errup = 0.197077; errdown = 0.197077;}
      else if (ht> 300 && ht<= 400 && met> 155 && met<= 160) {eff2017 = 0.440789; errup = 0.206929; errdown = 0.206929;}
      else if (ht> 300 && ht<= 400 && met> 160 && met<= 165) {eff2017 = 0.565657; errup = 0.0858744; errdown = 0.0858744;}
      else if (ht> 300 && ht<= 400 && met> 165 && met<= 170) {eff2017 = 0.5625; errup = 0.0866008; errdown = 0.0866008;}
      else if (ht> 300 && ht<= 400 && met> 170 && met<= 175) {eff2017 = 0.609649; errup = 0.0929996; errdown = 0.0929996;}
      else if (ht> 300 && ht<= 400 && met> 175 && met<= 180) {eff2017 = 0.696335; errup = 0.105018; errdown = 0.105018;}
      else if (ht> 300 && ht<= 400 && met> 180 && met<= 185) {eff2017 = 0.73913; errup = 0.178528; errdown = 0.178528;}
      else if (ht> 300 && ht<= 400 && met> 185 && met<= 190) {eff2017 = 0.783626; errup = 0.188336; errdown = 0.188336;}
      else if (ht> 300 && ht<= 400 && met> 190 && met<= 195) {eff2017 = 0.804348; errup = 0.193564; errdown = 0.193564;}
      else if (ht> 300 && ht<= 400 && met> 195 && met<= 200) {eff2017 = 0.845161; errup = 0.154839; errdown = 0.202364;}
      else if (ht> 300 && ht<= 400 && met> 200 && met<= 210) {eff2017 = 0.88664; errup = 0.0943624; errdown = 0.0943624;}
      else if (ht> 300 && ht<= 400 && met> 210 && met<= 220) {eff2017 = 0.925532; errup = 0.0744681; errdown = 0.098111;}
      else if (ht> 300 && ht<= 400 && met> 220 && met<= 230) {eff2017 = 0.956522; errup = 0.0434783; errdown = 0.100737;}
      else if (ht> 300 && ht<= 400 && met> 230 && met<= 240) {eff2017 = 0.978261; errup = 0.0217391; errdown = 0.0338717;}
      else if (ht> 300 && ht<= 400 && met> 240 && met<= 250) {eff2017 = 0.964912; errup = 0.0350877; errdown = 0.0355423;}
      else if (ht> 300 && ht<= 400 && met> 250 && met<= 275) {eff2017 = 0.970588; errup = 0.0227573; errdown = 0.0227573;}
      else if (ht> 300 && ht<= 400 && met> 275 && met<= 300) {eff2017 = 0.988235; errup = 0.0117647; errdown = 0.0229914;}
      else if (ht> 300 && ht<= 400 && met> 300 && met<= 9999) {eff2017 = 0.984848; errup = 0.0151515; errdown = 0.0912803;}
      else if (ht> 400 && ht<= 600 && met> 150 && met<= 155) {eff2017 = 0.407002; errup = 0.0642462; errdown = 0.0642462;}
      else if (ht> 400 && ht<= 600 && met> 155 && met<= 160) {eff2017 = 0.495192; errup = 0.0770015; errdown = 0.0770015;}
      else if (ht> 400 && ht<= 600 && met> 160 && met<= 165) {eff2017 = 0.532374; errup = 0.0815655; errdown = 0.0815655;}
      else if (ht> 400 && ht<= 600 && met> 165 && met<= 170) {eff2017 = 0.608808; errup = 0.0923941; errdown = 0.0923941;}
      else if (ht> 400 && ht<= 600 && met> 170 && met<= 175) {eff2017 = 0.593548; errup = 0.0911365; errdown = 0.0911365;}
      else if (ht> 400 && ht<= 600 && met> 175 && met<= 180) {eff2017 = 0.716981; errup = 0.107806; errdown = 0.107806;}
      else if (ht> 400 && ht<= 600 && met> 180 && met<= 185) {eff2017 = 0.716475; errup = 0.0803522; errdown = 0.0803522;}
      else if (ht> 400 && ht<= 600 && met> 185 && met<= 190) {eff2017 = 0.749049; errup = 0.0831922; errdown = 0.0831922;}
      else if (ht> 400 && ht<= 600 && met> 190 && met<= 195) {eff2017 = 0.833333; errup = 0.0908854; errdown = 0.0908854;}
      else if (ht> 400 && ht<= 600 && met> 195 && met<= 200) {eff2017 = 0.851282; errup = 0.0930868; errdown = 0.0930868;}
      else if (ht> 400 && ht<= 600 && met> 200 && met<= 210) {eff2017 = 0.879795; errup = 0.089772; errdown = 0.089772;}
      else if (ht> 400 && ht<= 600 && met> 210 && met<= 220) {eff2017 = 0.959375; errup = 0.040625; errdown = 0.0968661;}
      else if (ht> 400 && ht<= 600 && met> 220 && met<= 230) {eff2017 = 0.954064; errup = 0.0459364; errdown = 0.0965082;}
      else if (ht> 400 && ht<= 600 && met> 230 && met<= 240) {eff2017 = 0.953975; errup = 0.0213923; errdown = 0.0213923;}
      else if (ht> 400 && ht<= 600 && met> 240 && met<= 250) {eff2017 = 0.969543; errup = 0.0208045; errdown = 0.0208045;}
      else if (ht> 400 && ht<= 600 && met> 250 && met<= 275) {eff2017 = 0.978261; errup = 0.0217391; errdown = 0.0281273;}
      else if (ht> 400 && ht<= 600 && met> 275 && met<= 300) {eff2017 = 0.98773; errup = 0.0122699; errdown = 0.0281309;}
      else if (ht> 400 && ht<= 600 && met> 300 && met<= 9999) {eff2017 = 0.998663; errup = 0.0013369; errdown = 0.0840732;}
      else if (ht> 600 && ht<= 950 && met> 150 && met<= 155) {eff2017 = 0.415385; errup = 0.0737943; errdown = 0.0737943;}
      else if (ht> 600 && ht<= 950 && met> 155 && met<= 160) {eff2017 = 0.518828; errup = 0.0899055; errdown = 0.0899055;}
      else if (ht> 600 && ht<= 950 && met> 160 && met<= 165) {eff2017 = 0.555556; errup = 0.0702944; errdown = 0.0702944;}
      else if (ht> 600 && ht<= 950 && met> 165 && met<= 170) {eff2017 = 0.580488; errup = 0.0726663; errdown = 0.0726663;}
      else if (ht> 600 && ht<= 950 && met> 170 && met<= 175) {eff2017 = 0.615385; errup = 0.0768102; errdown = 0.0768102;}
      else if (ht> 600 && ht<= 950 && met> 175 && met<= 180) {eff2017 = 0.716667; errup = 0.0858251; errdown = 0.0858251;}
      else if (ht> 600 && ht<= 950 && met> 180 && met<= 185) {eff2017 = 0.675497; errup = 0.0734554; errdown = 0.0734554;}
      else if (ht> 600 && ht<= 950 && met> 185 && met<= 190) {eff2017 = 0.795775; errup = 0.0813517; errdown = 0.0813517;}
      else if (ht> 600 && ht<= 950 && met> 190 && met<= 195) {eff2017 = 0.742857; errup = 0.0783216; errdown = 0.0783216;}
      else if (ht> 600 && ht<= 950 && met> 195 && met<= 200) {eff2017 = 0.843284; errup = 0.0844566; errdown = 0.0844566;}
      else if (ht> 600 && ht<= 950 && met> 200 && met<= 210) {eff2017 = 0.885593; errup = 0.0605422; errdown = 0.0605422;}
      else if (ht> 600 && ht<= 950 && met> 210 && met<= 220) {eff2017 = 0.921053; errup = 0.0623144; errdown = 0.0623144;}
      else if (ht> 600 && ht<= 950 && met> 220 && met<= 230) {eff2017 = 0.940541; errup = 0.0594595; errdown = 0.0628678;}
      else if (ht> 600 && ht<= 950 && met> 230 && met<= 240) {eff2017 = 0.947368; errup = 0.0261934; errdown = 0.0261934;}
      else if (ht> 600 && ht<= 950 && met> 240 && met<= 250) {eff2017 = 0.983871; errup = 0.016129; errdown = 0.0215317;}
      else if (ht> 600 && ht<= 950 && met> 250 && met<= 275) {eff2017 = 0.979167; errup = 0.0208333; errdown = 0.0244204;}
      else if (ht> 600 && ht<= 950 && met> 275 && met<= 300) {eff2017 = 0.990148; errup = 0.00985222; errdown = 0.0241957;}
      else if (ht> 600 && ht<= 950 && met> 300 && met<= 9999) {eff2017 = 0.997106; errup = 0.00289436; errdown = 0.0843903;}
      else if (ht> 950 && ht<= 9999 && met> 150 && met<= 160) {eff2017 = 0.539216; errup = 0.0849513; errdown = 0.0849513;}
      else if (ht> 950 && ht<= 9999 && met> 160 && met<= 170) {eff2017 = 0.567901; errup = 0.0880225; errdown = 0.0880225;}
      else if (ht> 950 && ht<= 9999 && met> 170 && met<= 180) {eff2017 = 0.714286; errup = 0.101882; errdown = 0.101882;}
      else if (ht> 950 && ht<= 9999 && met> 180 && met<= 190) {eff2017 = 0.825397; errup = 0.16131; errdown = 0.16131;}
      else if (ht> 950 && ht<= 9999 && met> 190 && met<= 200) {eff2017 = 0.724138; errup = 0.147348; errdown = 0.147348;}
      else if (ht> 950 && ht<= 9999 && met> 200 && met<= 210) {eff2017 = 0.862745; errup = 0.119452; errdown = 0.119452;}
      else if (ht> 950 && ht<= 9999 && met> 210 && met<= 220) {eff2017 = 0.948718; errup = 0.0512821; errdown = 0.125276;}
      else if (ht> 950 && ht<= 9999 && met> 220 && met<= 230) {eff2017 = 0.875; errup = 0.122569; errdown = 0.122569;}
      else if (ht> 950 && ht<= 9999 && met> 230 && met<= 240) {eff2017 = 0.975; errup = 0.025; errdown = 0.032736;}
      else if (ht> 950 && ht<= 9999 && met> 240 && met<= 250) {eff2017 = 1; errup = 0; errdown = 0.0220517;}
      else if (ht> 950 && ht<= 9999 && met> 250 && met<= 275) {eff2017 = 1; errup = 0; errdown = 0.0198667;}
      else if (ht> 950 && ht<= 9999 && met> 275 && met<= 300) {eff2017 = 0.972222; errup = 0.0277778; errdown = 0.0335147;}
      else if (ht> 950 && ht<= 9999 && met> 300 && met<= 9999) {eff2017 = 1; errup = 0; errdown = 0.089518;}
      //2018
      if (ht> 0 && ht<= 200 && met> 150 && met<= 155) {eff2018 = 0.0519126; errup = 0.0266522; errdown = 0.0266522;}
      else if (ht> 0 && ht<= 200 && met> 155 && met<= 160) {eff2018 = 0.102473; errup = 0.0506841; errdown = 0.0506841;}
      else if (ht> 0 && ht<= 200 && met> 160 && met<= 165) {eff2018 = 0.111702; errup = 0.0336122; errdown = 0.0336122;}
      else if (ht> 0 && ht<= 200 && met> 165 && met<= 170) {eff2018 = 0.111702; errup = 0.0336122; errdown = 0.0336122;}
      else if (ht> 0 && ht<= 200 && met> 170 && met<= 180) {eff2018 = 0.177866; errup = 0.0458729; errdown = 0.0458729;}
      else if (ht> 0 && ht<= 200 && met> 180 && met<= 190) {eff2018 = 0.221311; errup = 0.0614236; errdown = 0.0614236;}
      else if (ht> 0 && ht<= 200 && met> 190 && met<= 200) {eff2018 = 0.272727; errup = 0.0747522; errdown = 0.0747522;}
      else if (ht> 0 && ht<= 200 && met> 200 && met<= 9999) {eff2018 = 0.336957; errup = 0.0757904; errdown = 0.0757904;}
      else if (ht> 200 && ht<= 300 && met> 150 && met<= 155) {eff2018 = 0.200972; errup = 0.0942924; errdown = 0.0942924;}
      else if (ht> 200 && ht<= 300 && met> 155 && met<= 160) {eff2018 = 0.259058; errup = 0.121196; errdown = 0.121196;}
      else if (ht> 200 && ht<= 300 && met> 160 && met<= 165) {eff2018 = 0.283133; errup = 0.0653851; errdown = 0.0653851;}
      else if (ht> 200 && ht<= 300 && met> 165 && met<= 170) {eff2018 = 0.331126; errup = 0.076019; errdown = 0.076019;}
      else if (ht> 200 && ht<= 300 && met> 170 && met<= 175) {eff2018 = 0.451128; errup = 0.102174; errdown = 0.102174;}
      else if (ht> 200 && ht<= 300 && met> 175 && met<= 180) {eff2018 = 0.41875; errup = 0.0960246; errdown = 0.0960246;}
      else if (ht> 200 && ht<= 300 && met> 180 && met<= 185) {eff2018 = 0.515504; errup = 0.117364; errdown = 0.117364;}
      else if (ht> 200 && ht<= 300 && met> 185 && met<= 190) {eff2018 = 0.579439; errup = 0.1316; errdown = 0.1316;}
      else if (ht> 200 && ht<= 300 && met> 190 && met<= 195) {eff2018 = 0.568807; errup = 0.129293; errdown = 0.129293;}
      else if (ht> 200 && ht<= 300 && met> 195 && met<= 200) {eff2018 = 0.670968; errup = 0.152051; errdown = 0.152051;}
      else if (ht> 200 && ht<= 300 && met> 200 && met<= 210) {eff2018 = 0.73516; errup = 0.129121; errdown = 0.129121;}
      else if (ht> 200 && ht<= 300 && met> 210 && met<= 220) {eff2018 = 0.774038; errup = 0.135416; errdown = 0.135416;}
      else if (ht> 200 && ht<= 300 && met> 220 && met<= 230) {eff2018 = 0.867089; errup = 0.132911; errdown = 0.150618;}
      else if (ht> 200 && ht<= 300 && met> 230 && met<= 240) {eff2018 = 0.87156; errup = 0.0909439; errdown = 0.0909439;}
      else if (ht> 200 && ht<= 300 && met> 240 && met<= 250) {eff2018 = 0.886076; errup = 0.093621; errdown = 0.093621;}
      else if (ht> 200 && ht<= 300 && met> 250 && met<= 275) {eff2018 = 0.935484; errup = 0.0645161; errdown = 0.0808035;}
      else if (ht> 200 && ht<= 300 && met> 275 && met<= 9999) {eff2018 = 0.972222; errup = 0.0277778; errdown = 0.0842695;}
      else if (ht> 300 && ht<= 400 && met> 150 && met<= 155) {eff2018 = 0.337349; errup = 0.163848; errdown = 0.163848;}
      else if (ht> 300 && ht<= 400 && met> 155 && met<= 160) {eff2018 = 0.360269; errup = 0.175003; errdown = 0.175003;}
      else if (ht> 300 && ht<= 400 && met> 160 && met<= 165) {eff2018 = 0.422925; errup = 0.0939919; errdown = 0.0939919;}
      else if (ht> 300 && ht<= 400 && met> 165 && met<= 170) {eff2018 = 0.466926; errup = 0.102767; errdown = 0.102767;}
      else if (ht> 300 && ht<= 400 && met> 170 && met<= 175) {eff2018 = 0.553398; errup = 0.121137; errdown = 0.121137;}
      else if (ht> 300 && ht<= 400 && met> 175 && met<= 180) {eff2018 = 0.569832; errup = 0.125124; errdown = 0.125124;}
      else if (ht> 300 && ht<= 400 && met> 180 && met<= 185) {eff2018 = 0.701493; errup = 0.134049; errdown = 0.134049;}
      else if (ht> 300 && ht<= 400 && met> 185 && met<= 190) {eff2018 = 0.700535; errup = 0.134175; errdown = 0.134175;}
      else if (ht> 300 && ht<= 400 && met> 190 && met<= 195) {eff2018 = 0.724138; errup = 0.13934; errdown = 0.13934;}
      else if (ht> 300 && ht<= 400 && met> 195 && met<= 200) {eff2018 = 0.826667; errup = 0.156406; errdown = 0.156406;}
      else if (ht> 300 && ht<= 400 && met> 200 && met<= 210) {eff2018 = 0.862069; errup = 0.0880559; errdown = 0.0880559;}
      else if (ht> 300 && ht<= 400 && met> 210 && met<= 220) {eff2018 = 0.867816; errup = 0.0890126; errdown = 0.0890126;}
      else if (ht> 300 && ht<= 400 && met> 220 && met<= 230) {eff2018 = 0.940541; errup = 0.0594595; errdown = 0.0939934;}
      else if (ht> 300 && ht<= 400 && met> 230 && met<= 240) {eff2018 = 0.950413; errup = 0.0401099; errdown = 0.0401099;}
      else if (ht> 300 && ht<= 400 && met> 240 && met<= 250) {eff2018 = 0.982609; errup = 0.0173913; errdown = 0.0381041;}
      else if (ht> 300 && ht<= 400 && met> 250 && met<= 275) {eff2018 = 0.985437; errup = 0.0145631; errdown = 0.0236853;}
      else if (ht> 300 && ht<= 400 && met> 275 && met<= 300) {eff2018 = 0.985507; errup = 0.0144928; errdown = 0.0243905;}
      else if (ht> 300 && ht<= 400 && met> 300 && met<= 9999) {eff2018 = 0.99; errup = 0.01; errdown = 0.0169057;}
      else if (ht> 400 && ht<= 600 && met> 150 && met<= 155) {eff2018 = 0.381271; errup = 0.107837; errdown = 0.107837;}
      else if (ht> 400 && ht<= 600 && met> 155 && met<= 160) {eff2018 = 0.43619; errup = 0.123176; errdown = 0.123176;}
      else if (ht> 400 && ht<= 600 && met> 160 && met<= 165) {eff2018 = 0.446581; errup = 0.0874369; errdown = 0.0874369;}
      else if (ht> 400 && ht<= 600 && met> 165 && met<= 170) {eff2018 = 0.558411; errup = 0.108185; errdown = 0.108185;}
      else if (ht> 400 && ht<= 600 && met> 170 && met<= 175) {eff2018 = 0.581047; errup = 0.112496; errdown = 0.112496;}
      else if (ht> 400 && ht<= 600 && met> 175 && met<= 180) {eff2018 = 0.677778; errup = 0.130386; errdown = 0.130386;}
      else if (ht> 400 && ht<= 600 && met> 180 && met<= 185) {eff2018 = 0.722071; errup = 0.0808773; errdown = 0.0808773;}
      else if (ht> 400 && ht<= 600 && met> 185 && met<= 190) {eff2018 = 0.780899; errup = 0.0865531; errdown = 0.0865531;}
      else if (ht> 400 && ht<= 600 && met> 190 && met<= 195) {eff2018 = 0.780488; errup = 0.0871803; errdown = 0.0871803;}
      else if (ht> 400 && ht<= 600 && met> 195 && met<= 200) {eff2018 = 0.857143; errup = 0.0943144; errdown = 0.0943144;}
      else if (ht> 400 && ht<= 600 && met> 200 && met<= 210) {eff2018 = 0.858427; errup = 0.0715949; errdown = 0.0715949;}
      else if (ht> 400 && ht<= 600 && met> 210 && met<= 220) {eff2018 = 0.925926; errup = 0.0740741; errdown = 0.0762577;}
      else if (ht> 400 && ht<= 600 && met> 220 && met<= 230) {eff2018 = 0.955844; errup = 0.0441558; errdown = 0.0782704;}
      else if (ht> 400 && ht<= 600 && met> 230 && met<= 240) {eff2018 = 0.978723; errup = 0.0212766; errdown = 0.0309743;}
      else if (ht> 400 && ht<= 600 && met> 240 && met<= 250) {eff2018 = 0.971631; errup = 0.0283688; errdown = 0.0313196;}
      else if (ht> 400 && ht<= 600 && met> 250 && met<= 275) {eff2018 = 0.994444; errup = 0.00555556; errdown = 0.022177;}
      else if (ht> 400 && ht<= 600 && met> 275 && met<= 300) {eff2018 = 0.987981; errup = 0.0120192; errdown = 0.0224476;}
      else if (ht> 400 && ht<= 600 && met> 300 && met<= 9999) {eff2018 = 1; errup = 0; errdown = 0.00358039;}
      else if (ht> 600 && ht<= 950 && met> 150 && met<= 155) {eff2018 = 0.400491; errup = 0.0744788; errdown = 0.0744788;}
      else if (ht> 600 && ht<= 950 && met> 155 && met<= 160) {eff2018 = 0.504021; errup = 0.0923124; errdown = 0.0923124;}
      else if (ht> 600 && ht<= 950 && met> 160 && met<= 165) {eff2018 = 0.507003; errup = 0.103808; errdown = 0.103808;}
      else if (ht> 600 && ht<= 950 && met> 165 && met<= 170) {eff2018 = 0.551155; errup = 0.1128; errdown = 0.1128;}
      else if (ht> 600 && ht<= 950 && met> 170 && met<= 175) {eff2018 = 0.640927; errup = 0.130348; errdown = 0.130348;}
      else if (ht> 600 && ht<= 950 && met> 175 && met<= 180) {eff2018 = 0.743191; errup = 0.149643; errdown = 0.149643;}
      else if (ht> 600 && ht<= 950 && met> 180 && met<= 185) {eff2018 = 0.725806; errup = 0.0445408; errdown = 0.0445408;}
      else if (ht> 600 && ht<= 950 && met> 185 && met<= 190) {eff2018 = 0.835498; errup = 0.046481; errdown = 0.046481;}
      else if (ht> 600 && ht<= 950 && met> 190 && met<= 195) {eff2018 = 0.797297; errup = 0.046407; errdown = 0.046407;}
      else if (ht> 600 && ht<= 950 && met> 195 && met<= 200) {eff2018 = 0.850299; errup = 0.0488229; errdown = 0.0488229;}
      else if (ht> 600 && ht<= 950 && met> 200 && met<= 210) {eff2018 = 0.886628; errup = 0.0654922; errdown = 0.0654922;}
      else if (ht> 600 && ht<= 950 && met> 210 && met<= 220) {eff2018 = 0.926357; errup = 0.068027; errdown = 0.068027;}
      else if (ht> 600 && ht<= 950 && met> 220 && met<= 230) {eff2018 = 0.933333; errup = 0.0666667; errdown = 0.0684725;}
      else if (ht> 600 && ht<= 950 && met> 230 && met<= 240) {eff2018 = 0.982759; errup = 0.0172414; errdown = 0.0310059;}
      else if (ht> 600 && ht<= 950 && met> 240 && met<= 250) {eff2018 = 0.980583; errup = 0.0194175; errdown = 0.0312543;}
      else if (ht> 600 && ht<= 950 && met> 250 && met<= 275) {eff2018 = 0.989333; errup = 0.0106667; errdown = 0.0226867;}
      else if (ht> 600 && ht<= 950 && met> 275 && met<= 300) {eff2018 = 0.996753; errup = 0.00324675; errdown = 0.0224584;}
      else if (ht> 600 && ht<= 950 && met> 300 && met<= 9999) {eff2018 = 0.996324; errup = 0.00367647; errdown = 0.00711638;}
      else if (ht> 950 && ht<= 9999 && met> 150 && met<= 160) {eff2018 = 0.493151; errup = 0.110839; errdown = 0.110839;}
      else if (ht> 950 && ht<= 9999 && met> 160 && met<= 170) {eff2018 = 0.644444; errup = 0.103253; errdown = 0.103253;}
      else if (ht> 950 && ht<= 9999 && met> 170 && met<= 180) {eff2018 = 0.681416; errup = 0.109285; errdown = 0.109285;}
      else if (ht> 950 && ht<= 9999 && met> 180 && met<= 190) {eff2018 = 0.863158; errup = 0.100786; errdown = 0.100786;}
      else if (ht> 950 && ht<= 9999 && met> 190 && met<= 200) {eff2018 = 0.847826; errup = 0.100015; errdown = 0.100015;}
      else if (ht> 950 && ht<= 9999 && met> 200 && met<= 210) {eff2018 = 0.903614; errup = 0.0574661; errdown = 0.0574661;}
      else if (ht> 950 && ht<= 9999 && met> 210 && met<= 220) {eff2018 = 0.904762; errup = 0.0602203; errdown = 0.0602203;}
      else if (ht> 950 && ht<= 9999 && met> 220 && met<= 230) {eff2018 = 0.953125; errup = 0.046875; errdown = 0.0566106;}
      else if (ht> 950 && ht<= 9999 && met> 230 && met<= 240) {eff2018 = 0.979167; errup = 0.0208333; errdown = 0.0411332;}
      else if (ht> 950 && ht<= 9999 && met> 240 && met<= 250) {eff2018 = 1; errup = 0; errdown = 0.0363517;}
      else if (ht> 950 && ht<= 9999 && met> 250 && met<= 275) {eff2018 = 0.965909; errup = 0.0291201; errdown = 0.0291201;}
      else if (ht> 950 && ht<= 9999 && met> 275 && met<= 300) {eff2018 = 1; errup = 0; errdown = 0.0225349;}
      else if (ht> 950 && ht<= 9999 && met> 300 && met<= 9999) {eff2018 = 1; errup = 0; errdown = 0.021753;}
    }
    else { //Fake MET MC
      if (ht> 0 && ht<= 200 && met> 150 && met<= 155) {eff2016 = 0.546973; errup = 0.114682; errdown = 0.114682;}
      else if (ht> 0 && ht<= 200 && met> 155 && met<= 160) {eff2016 = 0.57971; errup = 0.121576; errdown = 0.121576;}
      else if (ht> 0 && ht<= 200 && met> 160 && met<= 165) {eff2016 = 0.656766; errup = 0.134811; errdown = 0.134811;}
      else if (ht> 0 && ht<= 200 && met> 165 && met<= 170) {eff2016 = 0.719149; errup = 0.147506; errdown = 0.147506;}
      else if (ht> 0 && ht<= 200 && met> 170 && met<= 180) {eff2016 = 0.758621; errup = 0.154213; errdown = 0.154213;}
      else if (ht> 0 && ht<= 200 && met> 180 && met<= 190) {eff2016 = 0.817734; errup = 0.0773914; errdown = 0.0773914;}
      else if (ht> 0 && ht<= 200 && met> 190 && met<= 200) {eff2016 = 0.877193; errup = 0.0836193; errdown = 0.0836193;}
      else if (ht> 0 && ht<= 200 && met> 200 && met<= 9999) {eff2016 = 0.927419; errup = 0.0725806; errdown = 0.0731884;}
      else if (ht> 200 && ht<= 300 && met> 150 && met<= 155) {eff2016 = 0.731915; errup = 0.151102; errdown = 0.151102;}
      else if (ht> 200 && ht<= 300 && met> 155 && met<= 160) {eff2016 = 0.777531; errup = 0.160385; errdown = 0.160385;}
      else if (ht> 200 && ht<= 300 && met> 160 && met<= 165) {eff2016 = 0.797994; errup = 0.161131; errdown = 0.161131;}
      else if (ht> 200 && ht<= 300 && met> 165 && met<= 170) {eff2016 = 0.857357; errup = 0.142643; errdown = 0.172877;}
      else if (ht> 200 && ht<= 300 && met> 170 && met<= 175) {eff2016 = 0.902156; errup = 0.0978441; errdown = 0.181754;}
      else if (ht> 200 && ht<= 300 && met> 175 && met<= 180) {eff2016 = 0.915905; errup = 0.0840951; errdown = 0.184497;}
      else if (ht> 200 && ht<= 300 && met> 180 && met<= 185) {eff2016 = 0.926174; errup = 0.0738255; errdown = 0.0830325;}
      else if (ht> 200 && ht<= 300 && met> 185 && met<= 190) {eff2016 = 0.930233; errup = 0.0697674; errdown = 0.083602;}
      else if (ht> 200 && ht<= 300 && met> 190 && met<= 195) {eff2016 = 0.942779; errup = 0.0572207; errdown = 0.0844531;}
      else if (ht> 200 && ht<= 300 && met> 195 && met<= 200) {eff2016 = 0.969595; errup = 0.0304054; errdown = 0.0865329;}
      else if (ht> 200 && ht<= 300 && met> 200 && met<= 210) {eff2016 = 0.963039; errup = 0.036961; errdown = 0.072551;}
      else if (ht> 200 && ht<= 300 && met> 210 && met<= 220) {eff2016 = 0.976331; errup = 0.0236686; errdown = 0.0735065;}
      else if (ht> 200 && ht<= 300 && met> 220 && met<= 230) {eff2016 = 0.991803; errup = 0.00819672; errdown = 0.0744216;}
      else if (ht> 200 && ht<= 300 && met> 230 && met<= 240) {eff2016 = 0.98773; errup = 0.0122699; errdown = 0.0184646;}
      else if (ht> 200 && ht<= 300 && met> 240 && met<= 250) {eff2016 = 1; errup = 0; errdown = 0.0165304;}
      else if (ht> 200 && ht<= 300 && met> 250 && met<= 275) {eff2016 = 1; errup = 0; errdown = 0.00305686;}
      else if (ht> 200 && ht<= 300 && met> 275 && met<= 9999) {eff2016 = 1; errup = 0; errdown = 0.00305686;}
      else if (ht> 300 && ht<= 400 && met> 150 && met<= 155) {eff2016 = 0.747956; errup = 0.126351; errdown = 0.126351;}
      else if (ht> 300 && ht<= 400 && met> 155 && met<= 160) {eff2016 = 0.813636; errup = 0.137177; errdown = 0.137177;}
      else if (ht> 300 && ht<= 400 && met> 160 && met<= 165) {eff2016 = 0.847368; errup = 0.123161; errdown = 0.123161;}
      else if (ht> 300 && ht<= 400 && met> 165 && met<= 170) {eff2016 = 0.873264; errup = 0.126732; errdown = 0.126732;}
      else if (ht> 300 && ht<= 400 && met> 170 && met<= 175) {eff2016 = 0.918142; errup = 0.0818584; errdown = 0.133072;}
      else if (ht> 300 && ht<= 400 && met> 175 && met<= 180) {eff2016 = 0.928571; errup = 0.0714286; errdown = 0.13447;}
      else if (ht> 300 && ht<= 400 && met> 180 && met<= 185) {eff2016 = 0.93932; errup = 0.0450763; errdown = 0.0450763;}
      else if (ht> 300 && ht<= 400 && met> 185 && met<= 190) {eff2016 = 0.95599; errup = 0.0440098; errdown = 0.0454335;}
      else if (ht> 300 && ht<= 400 && met> 190 && met<= 195) {eff2016 = 0.96134; errup = 0.0386598; errdown = 0.0455975;}
      else if (ht> 300 && ht<= 400 && met> 195 && met<= 200) {eff2016 = 0.966777; errup = 0.0332226; errdown = 0.0459625;}
      else if (ht> 300 && ht<= 400 && met> 200 && met<= 210) {eff2016 = 0.970018; errup = 0.0299824; errdown = 0.0657562;}
      else if (ht> 300 && ht<= 400 && met> 210 && met<= 220) {eff2016 = 0.984716; errup = 0.0152838; errdown = 0.0666026;}
      else if (ht> 300 && ht<= 400 && met> 220 && met<= 230) {eff2016 = 0.994819; errup = 0.00518135; errdown = 0.0671357;}
      else if (ht> 300 && ht<= 400 && met> 230 && met<= 240) {eff2016 = 0.990881; errup = 0.00911854; errdown = 0.0193026;}
      else if (ht> 300 && ht<= 400 && met> 240 && met<= 250) {eff2016 = 1; errup = 0; errdown = 0.0187486;}
      else if (ht> 300 && ht<= 400 && met> 250 && met<= 275) {eff2016 = 0.994012; errup = 0.00598802; errdown = 0.00626348;}
      else if (ht> 300 && ht<= 400 && met> 275 && met<= 300) {eff2016 = 1; errup = 0; errdown = 0.00526128;}
      else if (ht> 300 && ht<= 400 && met> 300 && met<= 9999) {eff2016 = 1; errup = 0; errdown = 0;}
      else if (ht> 400 && ht<= 600 && met> 150 && met<= 155) {eff2016 = 0.731836; errup = 0.0641633; errdown = 0.0641633;}
      else if (ht> 400 && ht<= 600 && met> 155 && met<= 160) {eff2016 = 0.776899; errup = 0.0679826; errdown = 0.0679826;}
      else if (ht> 400 && ht<= 600 && met> 160 && met<= 165) {eff2016 = 0.798233; errup = 0.0841752; errdown = 0.0841752;}
      else if (ht> 400 && ht<= 600 && met> 165 && met<= 170) {eff2016 = 0.854839; errup = 0.0898686; errdown = 0.0898686;}
      else if (ht> 400 && ht<= 600 && met> 170 && met<= 175) {eff2016 = 0.882243; errup = 0.0925193; errdown = 0.0925193;}
      else if (ht> 400 && ht<= 600 && met> 175 && met<= 180) {eff2016 = 0.921971; errup = 0.0780287; errdown = 0.0963521;}
      else if (ht> 400 && ht<= 600 && met> 180 && met<= 185) {eff2016 = 0.933638; errup = 0.0605031; errdown = 0.0605031;}
      else if (ht> 400 && ht<= 600 && met> 185 && met<= 190) {eff2016 = 0.938875; errup = 0.0608173; errdown = 0.0608173;}
      else if (ht> 400 && ht<= 600 && met> 190 && met<= 195) {eff2016 = 0.944444; errup = 0.0555556; errdown = 0.0611521;}
      else if (ht> 400 && ht<= 600 && met> 195 && met<= 200) {eff2016 = 0.967655; errup = 0.032345; errdown = 0.0621634;}
      else if (ht> 400 && ht<= 600 && met> 200 && met<= 210) {eff2016 = 0.977124; errup = 0.0228758; errdown = 0.067068;}
      else if (ht> 400 && ht<= 600 && met> 210 && met<= 220) {eff2016 = 0.975701; errup = 0.0242991; errdown = 0.0670292;}
      else if (ht> 400 && ht<= 600 && met> 220 && met<= 230) {eff2016 = 0.995726; errup = 0.0042735; errdown = 0.0681335;}
      else if (ht> 400 && ht<= 600 && met> 230 && met<= 240) {eff2016 = 0.995316; errup = 0.00468384; errdown = 0.0153583;}
      else if (ht> 400 && ht<= 600 && met> 240 && met<= 250) {eff2016 = 0.994505; errup = 0.00549451; errdown = 0.0154792;}
      else if (ht> 400 && ht<= 600 && met> 250 && met<= 275) {eff2016 = 0.994413; errup = 0.00459548; errdown = 0.00459548;}
      else if (ht> 400 && ht<= 600 && met> 275 && met<= 300) {eff2016 = 1; errup = 0; errdown = 0.0036756;}
      else if (ht> 400 && ht<= 600 && met> 300 && met<= 9999) {eff2016 = 1; errup = 0; errdown = 0;}
      else if (ht> 600 && ht<= 950 && met> 150 && met<= 155) {eff2016 = 0.762332; errup = 0.0497612; errdown = 0.0497612;}
      else if (ht> 600 && ht<= 950 && met> 155 && met<= 160) {eff2016 = 0.72; errup = 0.0487859; errdown = 0.0487859;}
      else if (ht> 600 && ht<= 950 && met> 160 && met<= 165) {eff2016 = 0.813725; errup = 0.0946634; errdown = 0.0946634;}
      else if (ht> 600 && ht<= 950 && met> 165 && met<= 170) {eff2016 = 0.832402; errup = 0.0968458; errdown = 0.0968458;}
      else if (ht> 600 && ht<= 950 && met> 170 && met<= 175) {eff2016 = 0.870748; errup = 0.100876; errdown = 0.100876;}
      else if (ht> 600 && ht<= 950 && met> 175 && met<= 180) {eff2016 = 0.927152; errup = 0.0728477; errdown = 0.105433;}
      else if (ht> 600 && ht<= 950 && met> 180 && met<= 185) {eff2016 = 0.916667; errup = 0.0336675; errdown = 0.0336675;}
      else if (ht> 600 && ht<= 950 && met> 185 && met<= 190) {eff2016 = 0.929078; errup = 0.0329663; errdown = 0.0329663;}
      else if (ht> 600 && ht<= 950 && met> 190 && met<= 195) {eff2016 = 0.883929; errup = 0.0384289; errdown = 0.0384289;}
      else if (ht> 600 && ht<= 950 && met> 195 && met<= 200) {eff2016 = 0.960784; errup = 0.0321224; errdown = 0.0321224;}
      else if (ht> 600 && ht<= 950 && met> 200 && met<= 210) {eff2016 = 0.962791; errup = 0.0372093; errdown = 0.0750713;}
      else if (ht> 600 && ht<= 950 && met> 210 && met<= 220) {eff2016 = 0.972222; errup = 0.0277778; errdown = 0.0756755;}
      else if (ht> 600 && ht<= 950 && met> 220 && met<= 230) {eff2016 = 0.964286; errup = 0.0357143; errdown = 0.0754391;}
      else if (ht> 600 && ht<= 950 && met> 230 && met<= 240) {eff2016 = 1; errup = 0; errdown = 0.0148911;}
      else if (ht> 600 && ht<= 950 && met> 240 && met<= 250) {eff2016 = 1; errup = 0; errdown = 0.0148911;}
      else if (ht> 600 && ht<= 950 && met> 250 && met<= 275) {eff2016 = 1; errup = 0; errdown = 0.00305686;}
      else if (ht> 600 && ht<= 950 && met> 275 && met<= 300) {eff2016 = 1; errup = 0; errdown = 0.00305686;}
      else if (ht> 600 && ht<= 950 && met> 300 && met<= 9999) {eff2016 = 1; errup = 0; errdown = 9.2029e-05;}
      else if (ht> 950 && ht<= 9999 && met> 150 && met<= 160) {eff2016 = 0.686275; errup = 0.154105; errdown = 0.154105;}
      else if (ht> 950 && ht<= 9999 && met> 160 && met<= 170) {eff2016 = 0.698413; errup = 0.119721; errdown = 0.119721;}
      else if (ht> 950 && ht<= 9999 && met> 170 && met<= 180) {eff2016 = 0.796875; errup = 0.129753; errdown = 0.129753;}
      else if (ht> 950 && ht<= 9999 && met> 180 && met<= 190) {eff2016 = 0.836735; errup = 0.0579891; errdown = 0.0579891;}
      else if (ht> 950 && ht<= 9999 && met> 190 && met<= 200) {eff2016 = 0.90625; errup = 0.0447411; errdown = 0.0447411;}
      else if (ht> 950 && ht<= 9999 && met> 200 && met<= 210) {eff2016 = 0.9375; errup = 0.0625; errdown = 0.0722521;}
      else if (ht> 950 && ht<= 9999 && met> 210 && met<= 220) {eff2016 = 0.925; errup = 0.075; errdown = 0.0750206;}
      else if (ht> 950 && ht<= 9999 && met> 220 && met<= 230) {eff2016 = 0.882353; errup = 0.0812162; errdown = 0.0812162;}
      else if (ht> 950 && ht<= 9999 && met> 230 && met<= 240) {eff2016 = 1; errup = 0; errdown = 0.0340618;}
      else if (ht> 950 && ht<= 9999 && met> 240 && met<= 250) {eff2016 = 1; errup = 0; errdown = 0.0340618;}
      else if (ht> 950 && ht<= 9999 && met> 250 && met<= 275) {eff2016 = 1; errup = 0; errdown = 0.00305732;}
      else if (ht> 950 && ht<= 9999 && met> 275 && met<= 300) {eff2016 = 1; errup = 0; errdown = 0.00305732;}
      else if (ht> 950 && ht<= 9999 && met> 300 && met<= 9999) {eff2016 = 1; errup = 0; errdown = 2.32418e-05;}
      //2017
      if (ht> 0 && ht<= 350 && met> 150 && met<= 160) {eff2017 = 0.0974304; errup = 0.00970318; errdown = 0.00970318;}
      else if (ht> 0 && ht<= 350 && met> 160 && met<= 170) {eff2017 = 0.14405; errup = 0.016044; errdown = 0.016044;}
      else if (ht> 0 && ht<= 350 && met> 170 && met<= 180) {eff2017 = 0.184874; errup = 0.025163; errdown = 0.025163;}
      else if (ht> 0 && ht<= 350 && met> 180 && met<= 190) {eff2017 = 0.321101; errup = 0.0447209; errdown = 0.0447209;}
      else if (ht> 0 && ht<= 350 && met> 190 && met<= 200) {eff2017 = 0.348837; errup = 0.0513934; errdown = 0.0513934;}
      else if (ht> 0 && ht<= 350 && met> 200 && met<= 225) {eff2017 = 0.545455; errup = 0.0567443; errdown = 0.0567443;}
      else if (ht> 0 && ht<= 350 && met> 225 && met<= 250) {eff2017 = 0.62069; errup = 0.0901022; errdown = 0.0901022;}
      else if (ht> 0 && ht<= 350 && met> 250 && met<= 9999) {eff2017 = 0.818182; errup = 0.116291; errdown = 0.116291;}
      else if (ht> 350 && ht<= 450 && met> 150 && met<= 155) {eff2017 = 0.148239; errup = 0.0101691; errdown = 0.0101691;}
      else if (ht> 350 && ht<= 450 && met> 155 && met<= 160) {eff2017 = 0.170517; errup = 0.012474; errdown = 0.012474;}
      else if (ht> 350 && ht<= 450 && met> 160 && met<= 165) {eff2017 = 0.214286; errup = 0.0150635; errdown = 0.0150635;}
      else if (ht> 350 && ht<= 450 && met> 165 && met<= 170) {eff2017 = 0.241993; errup = 0.0180663; errdown = 0.0180663;}
      else if (ht> 350 && ht<= 450 && met> 170 && met<= 175) {eff2017 = 0.277533; errup = 0.0210154; errdown = 0.0210154;}
      else if (ht> 350 && ht<= 450 && met> 175 && met<= 180) {eff2017 = 0.369444; errup = 0.0254381; errdown = 0.0254381;}
      else if (ht> 350 && ht<= 450 && met> 180 && met<= 185) {eff2017 = 0.395833; errup = 0.0315667; errdown = 0.0315667;}
      else if (ht> 350 && ht<= 450 && met> 185 && met<= 190) {eff2017 = 0.424731; errup = 0.036244; errdown = 0.036244;}
      else if (ht> 350 && ht<= 450 && met> 190 && met<= 195) {eff2017 = 0.486188; errup = 0.0371505; errdown = 0.0371505;}
      else if (ht> 350 && ht<= 450 && met> 195 && met<= 200) {eff2017 = 0.51049; errup = 0.0418029; errdown = 0.0418029;}
      else if (ht> 350 && ht<= 450 && met> 200 && met<= 210) {eff2017 = 0.622642; errup = 0.0332911; errdown = 0.0332911;}
      else if (ht> 350 && ht<= 450 && met> 210 && met<= 220) {eff2017 = 0.669173; errup = 0.0407985; errdown = 0.0407985;}
      else if (ht> 350 && ht<= 450 && met> 220 && met<= 230) {eff2017 = 0.728261; errup = 0.0463795; errdown = 0.0463795;}
      else if (ht> 350 && ht<= 450 && met> 230 && met<= 240) {eff2017 = 0.709677; errup = 0.0576468; errdown = 0.0576468;}
      else if (ht> 350 && ht<= 450 && met> 240 && met<= 250) {eff2017 = 0.757576; errup = 0.0746009; errdown = 0.0746009;}
      else if (ht> 350 && ht<= 450 && met> 250 && met<= 300) {eff2017 = 0.916031; errup = 0.0242315; errdown = 0.0242315;}
      else if (ht> 350 && ht<= 450 && met> 300 && met<= 9999) {eff2017 = 0.954545; errup = 0.0314022; errdown = 0.0314022;}
      else if (ht> 450 && ht<= 550 && met> 150 && met<= 155) {eff2017 = 0.189972; errup = 0.00841328; errdown = 0.00841328;}
      else if (ht> 450 && ht<= 550 && met> 155 && met<= 160) {eff2017 = 0.228241; errup = 0.0103543; errdown = 0.0103543;}
      else if (ht> 450 && ht<= 550 && met> 160 && met<= 165) {eff2017 = 0.243159; errup = 0.0119953; errdown = 0.0119953;}
      else if (ht> 450 && ht<= 550 && met> 165 && met<= 170) {eff2017 = 0.323308; errup = 0.0153295; errdown = 0.0153295;}
      else if (ht> 450 && ht<= 550 && met> 170 && met<= 175) {eff2017 = 0.377577; errup = 0.0174026; errdown = 0.0174026;}
      else if (ht> 450 && ht<= 550 && met> 175 && met<= 180) {eff2017 = 0.407874; errup = 0.0195022; errdown = 0.0195022;}
      else if (ht> 450 && ht<= 550 && met> 180 && met<= 185) {eff2017 = 0.418468; errup = 0.0218655; errdown = 0.0218655;}
      else if (ht> 450 && ht<= 550 && met> 185 && met<= 190) {eff2017 = 0.486339; errup = 0.0261257; errdown = 0.0261257;}
      else if (ht> 450 && ht<= 550 && met> 190 && met<= 195) {eff2017 = 0.454286; errup = 0.0266142; errdown = 0.0266142;}
      else if (ht> 450 && ht<= 550 && met> 195 && met<= 200) {eff2017 = 0.582143; errup = 0.0294747; errdown = 0.0294747;}
      else if (ht> 450 && ht<= 550 && met> 200 && met<= 210) {eff2017 = 0.590571; errup = 0.0244947; errdown = 0.0244947;}
      else if (ht> 450 && ht<= 550 && met> 210 && met<= 220) {eff2017 = 0.671141; errup = 0.0272147; errdown = 0.0272147;}
      else if (ht> 450 && ht<= 550 && met> 220 && met<= 230) {eff2017 = 0.795349; errup = 0.0275148; errdown = 0.0275148;}
      else if (ht> 450 && ht<= 550 && met> 230 && met<= 240) {eff2017 = 0.823944; errup = 0.0319617; errdown = 0.0319617;}
      else if (ht> 450 && ht<= 550 && met> 240 && met<= 250) {eff2017 = 0.796875; errup = 0.0355608; errdown = 0.0355608;}
      else if (ht> 450 && ht<= 550 && met> 250 && met<= 300) {eff2017 = 0.911972; errup = 0.0168129; errdown = 0.0168129;}
      else if (ht> 450 && ht<= 550 && met> 300 && met<= 400) {eff2017 = 0.975309; errup = 0.0121923; errdown = 0.0121923;}
      else if (ht> 450 && ht<= 550 && met> 400 && met<= 9999) {eff2017 = 0.966667; errup = 0.0231741; errdown = 0.0231741;}
      else if (ht> 550 && ht<= 650 && met> 150 && met<= 155) {eff2017 = 0.245086; errup = 0.00879666; errdown = 0.00879666;}
      else if (ht> 550 && ht<= 650 && met> 155 && met<= 160) {eff2017 = 0.279679; errup = 0.0103794; errdown = 0.0103794;}
      else if (ht> 550 && ht<= 650 && met> 160 && met<= 165) {eff2017 = 0.319121; errup = 0.0118475; errdown = 0.0118475;}
      else if (ht> 550 && ht<= 650 && met> 165 && met<= 170) {eff2017 = 0.363946; errup = 0.0140301; errdown = 0.0140301;}
      else if (ht> 550 && ht<= 650 && met> 170 && met<= 175) {eff2017 = 0.376166; errup = 0.0155941; errdown = 0.0155941;}
      else if (ht> 550 && ht<= 650 && met> 175 && met<= 180) {eff2017 = 0.445813; errup = 0.0174432; errdown = 0.0174432;}
      else if (ht> 550 && ht<= 650 && met> 180 && met<= 185) {eff2017 = 0.507184; errup = 0.0189505; errdown = 0.0189505;}
      else if (ht> 550 && ht<= 650 && met> 185 && met<= 190) {eff2017 = 0.532957; errup = 0.021651; errdown = 0.021651;}
      else if (ht> 550 && ht<= 650 && met> 190 && met<= 195) {eff2017 = 0.599109; errup = 0.0231283; errdown = 0.0231283;}
      else if (ht> 550 && ht<= 650 && met> 195 && met<= 200) {eff2017 = 0.631043; errup = 0.02434; errdown = 0.02434;}
      else if (ht> 550 && ht<= 650 && met> 200 && met<= 210) {eff2017 = 0.685155; errup = 0.0187591; errdown = 0.0187591;}
      else if (ht> 550 && ht<= 650 && met> 210 && met<= 220) {eff2017 = 0.694382; errup = 0.0218378; errdown = 0.0218378;}
      else if (ht> 550 && ht<= 650 && met> 220 && met<= 230) {eff2017 = 0.740413; errup = 0.0238111; errdown = 0.0238111;}
      else if (ht> 550 && ht<= 650 && met> 230 && met<= 240) {eff2017 = 0.8; errup = 0.0245718; errdown = 0.0245718;}
      else if (ht> 550 && ht<= 650 && met> 240 && met<= 250) {eff2017 = 0.856436; errup = 0.0246715; errdown = 0.0246715;}
      else if (ht> 550 && ht<= 650 && met> 250 && met<= 275) {eff2017 = 0.880342; errup = 0.0173238; errdown = 0.0173238;}
      else if (ht> 550 && ht<= 650 && met> 275 && met<= 300) {eff2017 = 0.903382; errup = 0.0205343; errdown = 0.0205343;}
      else if (ht> 550 && ht<= 650 && met> 300 && met<= 400) {eff2017 = 0.933333; errup = 0.0117589; errdown = 0.0117589;}
      else if (ht> 550 && ht<= 650 && met> 400 && met<= 9999) {eff2017 = 0.978102; errup = 0.00510454; errdown = 0.00510454;}
      else if (ht> 650 && ht<= 800 && met> 150 && met<= 155) {eff2017 = 0.322707; errup = 0.00710557; errdown = 0.00710557;}
      else if (ht> 650 && ht<= 800 && met> 155 && met<= 160) {eff2017 = 0.366185; errup = 0.00810742; errdown = 0.00810742;}
      else if (ht> 650 && ht<= 800 && met> 160 && met<= 165) {eff2017 = 0.406324; errup = 0.00905657; errdown = 0.00905657;}
      else if (ht> 650 && ht<= 800 && met> 165 && met<= 170) {eff2017 = 0.439822; errup = 0.00997531; errdown = 0.00997531;}
      else if (ht> 650 && ht<= 800 && met> 170 && met<= 175) {eff2017 = 0.463355; errup = 0.0110594; errdown = 0.0110594;}
      else if (ht> 650 && ht<= 800 && met> 175 && met<= 180) {eff2017 = 0.515652; errup = 0.0119226; errdown = 0.0119226;}
      else if (ht> 650 && ht<= 800 && met> 180 && met<= 185) {eff2017 = 0.554072; errup = 0.0128428; errdown = 0.0128428;}
      else if (ht> 650 && ht<= 800 && met> 185 && met<= 190) {eff2017 = 0.598348; errup = 0.0134323; errdown = 0.0134323;}
      else if (ht> 650 && ht<= 800 && met> 190 && met<= 195) {eff2017 = 0.642009; errup = 0.0144877; errdown = 0.0144877;}
      else if (ht> 650 && ht<= 800 && met> 195 && met<= 200) {eff2017 = 0.696538; errup = 0.0146713; errdown = 0.0146713;}
      else if (ht> 650 && ht<= 800 && met> 200 && met<= 210) {eff2017 = 0.715486; errup = 0.0110539; errdown = 0.0110539;}
      else if (ht> 650 && ht<= 800 && met> 210 && met<= 220) {eff2017 = 0.781983; errup = 0.0113604; errdown = 0.0113604;}
      else if (ht> 650 && ht<= 800 && met> 220 && met<= 230) {eff2017 = 0.827618; errup = 0.0114987; errdown = 0.0114987;}
      else if (ht> 650 && ht<= 800 && met> 230 && met<= 240) {eff2017 = 0.839763; errup = 0.0115368; errdown = 0.0115368;}
      else if (ht> 650 && ht<= 800 && met> 240 && met<= 250) {eff2017 = 0.884467; errup = 0.0114531; errdown = 0.0114531;}
      else if (ht> 650 && ht<= 800 && met> 250 && met<= 275) {eff2017 = 0.898712; errup = 0.00730038; errdown = 0.00730038;}
      else if (ht> 650 && ht<= 800 && met> 275 && met<= 300) {eff2017 = 0.930855; errup = 0.00691768; errdown = 0.00691768;}
      else if (ht> 650 && ht<= 800 && met> 300 && met<= 350) {eff2017 = 0.941468; errup = 0.00522821; errdown = 0.00522821;}
      else if (ht> 650 && ht<= 800 && met> 350 && met<= 400) {eff2017 = 0.959596; errup = 0.00528902; errdown = 0.00528902;}
      else if (ht> 650 && ht<= 800 && met> 400 && met<= 450) {eff2017 = 0.976963; errup = 0.00485453; errdown = 0.00485453;}
      else if (ht> 650 && ht<= 800 && met> 450 && met<= 500) {eff2017 = 0.991124; errup = 0.00360739; errdown = 0.00360739;}
      else if (ht> 650 && ht<= 800 && met> 500 && met<= 9999) {eff2017 = 0.973333; errup = 0.00515956; errdown = 0.00515956;}
      else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff2017 = 0.489496; errup = 0.00369257; errdown = 0.00369257;}
      else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff2017 = 0.522675; errup = 0.00395583; errdown = 0.00395583;}
      else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff2017 = 0.568861; errup = 0.00421634; errdown = 0.00421634;}
      else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff2017 = 0.592914; errup = 0.00450707; errdown = 0.00450707;}
      else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff2017 = 0.635472; errup = 0.00469609; errdown = 0.00469609;}
      else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff2017 = 0.674262; errup = 0.00491711; errdown = 0.00491711;}
      else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff2017 = 0.686402; errup = 0.00517974; errdown = 0.00517974;}
      else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff2017 = 0.726592; errup = 0.00534945; errdown = 0.00534945;}
      else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff2017 = 0.737255; errup = 0.00551234; errdown = 0.00551234;}
      else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff2017 = 0.775837; errup = 0.00559582; errdown = 0.00559582;}
      else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff2017 = 0.802235; errup = 0.00410939; errdown = 0.00410939;}
      else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff2017 = 0.83659; errup = 0.00423424; errdown = 0.00423424;}
      else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff2017 = 0.849854; errup = 0.00442152; errdown = 0.00442152;}
      else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff2017 = 0.880527; errup = 0.00454664; errdown = 0.00454664;}
      else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff2017 = 0.893832; errup = 0.00466479; errdown = 0.00466479;}
      else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff2017 = 0.916677; errup = 0.00306379; errdown = 0.00306379;}
      else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff2017 = 0.936391; errup = 0.00323515; errdown = 0.00323515;}
      else if (ht> 800 && ht<= 1000 && met> 300 && met<= 350) {eff2017 = 0.941998; errup = 0.0028653; errdown = 0.0028653;}
      else if (ht> 800 && ht<= 1000 && met> 350 && met<= 400) {eff2017 = 0.952878; errup = 0.00376834; errdown = 0.00376834;}
      else if (ht> 800 && ht<= 1000 && met> 400 && met<= 450) {eff2017 = 0.963603; errup = 0.00465149; errdown = 0.00465149;}
      else if (ht> 800 && ht<= 1000 && met> 450 && met<= 500) {eff2017 = 0.954495; errup = 0.00694311; errdown = 0.00694311;}
      else if (ht> 800 && ht<= 1000 && met> 500 && met<= 9999) {eff2017 = 0.983751; errup = 0.00328979; errdown = 0.00328979;}
      else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff2017 = 0.387888; errup = 0.001472; errdown = 0.001472;}
      else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff2017 = 0.417927; errup = 0.00160609; errdown = 0.00160609;}
      else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff2017 = 0.451373; errup = 0.0017459; errdown = 0.0017459;}
      else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff2017 = 0.475537; errup = 0.00189276; errdown = 0.00189276;}
      else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff2017 = 0.512738; errup = 0.00204099; errdown = 0.00204099;}
      else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff2017 = 0.53824; errup = 0.00218953; errdown = 0.00218953;}
      else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff2017 = 0.567544; errup = 0.0023476; errdown = 0.0023476;}
      else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff2017 = 0.598143; errup = 0.00249669; errdown = 0.00249669;}
      else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff2017 = 0.624385; errup = 0.00264351; errdown = 0.00264351;}
      else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff2017 = 0.649303; errup = 0.00279531; errdown = 0.00279531;}
      else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff2017 = 0.681851; errup = 0.002131; errdown = 0.002131;}
      else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff2017 = 0.729521; errup = 0.00231462; errdown = 0.00231462;}
      else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff2017 = 0.762535; errup = 0.00251095; errdown = 0.00251095;}
      else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff2017 = 0.798438; errup = 0.00267273; errdown = 0.00267273;}
      else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff2017 = 0.821241; errup = 0.00285038; errdown = 0.00285038;}
      else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff2017 = 0.856744; errup = 0.00196232; errdown = 0.00196232;}
      else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff2017 = 0.886954; errup = 0.00223653; errdown = 0.00223653;}
      else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 350) {eff2017 = 0.916497; errup = 0.00185908; errdown = 0.00185908;}
      else if (ht> 1000 && ht<= 9999 && met> 350 && met<= 400) {eff2017 = 0.936176; errup = 0.00234582; errdown = 0.00234582;}
      else if (ht> 1000 && ht<= 9999 && met> 400 && met<= 450) {eff2017 = 0.947549; errup = 0.00292349; errdown = 0.00292349;}
      else if (ht> 1000 && ht<= 9999 && met> 450 && met<= 500) {eff2017 = 0.954474; errup = 0.0037065; errdown = 0.0037065;}
      else if (ht> 1000 && ht<= 9999 && met> 500 && met<= 9999) {eff2017 = 0.967807; errup = 0.00238727; errdown = 0.00238727;}
      //2018
      if (ht> 0 && ht<= 350 && met> 150 && met<= 160) {eff2018 = 0.078; errup = 0.011993; errdown = 0.011993;}
      else if (ht> 0 && ht<= 350 && met> 160 && met<= 170) {eff2018 = 0.102837; errup = 0.0180878; errdown = 0.0180878;}
      else if (ht> 0 && ht<= 350 && met> 170 && met<= 180) {eff2018 = 0.186207; errup = 0.0323274; errdown = 0.0323274;}
      else if (ht> 0 && ht<= 350 && met> 180 && met<= 190) {eff2018 = 0.253731; errup = 0.0531615; errdown = 0.0531615;}
      else if (ht> 0 && ht<= 350 && met> 190 && met<= 200) {eff2018 = 0.363636; errup = 0.0725204; errdown = 0.0725204;}
      else if (ht> 0 && ht<= 350 && met> 200 && met<= 225) {eff2018 = 0.428571; errup = 0.0763604; errdown = 0.0763604;}
      else if (ht> 0 && ht<= 350 && met> 225 && met<= 250) {eff2018 = 0.529412; errup = 0.121058; errdown = 0.121058;}
      else if (ht> 0 && ht<= 350 && met> 250 && met<= 9999) {eff2018 = 0.823529; errup = 0.0924594; errdown = 0.0924594;}
      else if (ht> 350 && ht<= 450 && met> 150 && met<= 155) {eff2018 = 0.153846; errup = 0.0114786; errdown = 0.0114786;}
      else if (ht> 350 && ht<= 450 && met> 155 && met<= 160) {eff2018 = 0.164521; errup = 0.014023; errdown = 0.014023;}
      else if (ht> 350 && ht<= 450 && met> 160 && met<= 165) {eff2018 = 0.227273; errup = 0.0178692; errdown = 0.0178692;}
      else if (ht> 350 && ht<= 450 && met> 165 && met<= 170) {eff2018 = 0.239374; errup = 0.0201823; errdown = 0.0201823;}
      else if (ht> 350 && ht<= 450 && met> 170 && met<= 175) {eff2018 = 0.285294; errup = 0.024489; errdown = 0.024489;}
      else if (ht> 350 && ht<= 450 && met> 175 && met<= 180) {eff2018 = 0.309353; errup = 0.0277225; errdown = 0.0277225;}
      else if (ht> 350 && ht<= 450 && met> 180 && met<= 185) {eff2018 = 0.402655; errup = 0.0326231; errdown = 0.0326231;}
      else if (ht> 350 && ht<= 450 && met> 185 && met<= 190) {eff2018 = 0.411111; errup = 0.0366741; errdown = 0.0366741;}
      else if (ht> 350 && ht<= 450 && met> 190 && met<= 195) {eff2018 = 0.503876; errup = 0.0440212; errdown = 0.0440212;}
      else if (ht> 350 && ht<= 450 && met> 195 && met<= 200) {eff2018 = 0.473684; errup = 0.0512278; errdown = 0.0512278;}
      else if (ht> 350 && ht<= 450 && met> 200 && met<= 210) {eff2018 = 0.649351; errup = 0.0384517; errdown = 0.0384517;}
      else if (ht> 350 && ht<= 450 && met> 210 && met<= 220) {eff2018 = 0.666667; errup = 0.0488824; errdown = 0.0488824;}
      else if (ht> 350 && ht<= 450 && met> 220 && met<= 230) {eff2018 = 0.860759; errup = 0.0389502; errdown = 0.0389502;}
      else if (ht> 350 && ht<= 450 && met> 230 && met<= 240) {eff2018 = 0.722222; errup = 0.0609519; errdown = 0.0609519;}
      else if (ht> 350 && ht<= 450 && met> 240 && met<= 250) {eff2018 = 0.837209; errup = 0.0562986; errdown = 0.0562986;}
      else if (ht> 350 && ht<= 450 && met> 250 && met<= 300) {eff2018 = 0.93617; errup = 0.025213; errdown = 0.025213;}
      else if (ht> 350 && ht<= 450 && met> 300 && met<= 9999) {eff2018 = 0.973684; errup = 0.0259672; errdown = 0.0259672;}
      else if (ht> 450 && ht<= 550 && met> 150 && met<= 155) {eff2018 = 0.185145; errup = 0.0090772; errdown = 0.0090772;}
      else if (ht> 450 && ht<= 550 && met> 155 && met<= 160) {eff2018 = 0.224373; errup = 0.0111693; errdown = 0.0111693;}
      else if (ht> 450 && ht<= 550 && met> 160 && met<= 165) {eff2018 = 0.253235; errup = 0.0132203; errdown = 0.0132203;}
      else if (ht> 450 && ht<= 550 && met> 165 && met<= 170) {eff2018 = 0.296253; errup = 0.0156247; errdown = 0.0156247;}
      else if (ht> 450 && ht<= 550 && met> 170 && met<= 175) {eff2018 = 0.351145; errup = 0.0186508; errdown = 0.0186508;}
      else if (ht> 450 && ht<= 550 && met> 175 && met<= 180) {eff2018 = 0.397363; errup = 0.0212361; errdown = 0.0212361;}
      else if (ht> 450 && ht<= 550 && met> 180 && met<= 185) {eff2018 = 0.430233; errup = 0.0238763; errdown = 0.0238763;}
      else if (ht> 450 && ht<= 550 && met> 185 && met<= 190) {eff2018 = 0.512894; errup = 0.0267555; errdown = 0.0267555;}
      else if (ht> 450 && ht<= 550 && met> 190 && met<= 195) {eff2018 = 0.553191; errup = 0.0296056; errdown = 0.0296056;}
      else if (ht> 450 && ht<= 550 && met> 195 && met<= 200) {eff2018 = 0.571429; errup = 0.033065; errdown = 0.033065;}
      else if (ht> 450 && ht<= 550 && met> 200 && met<= 210) {eff2018 = 0.626263; errup = 0.0243116; errdown = 0.0243116;}
      else if (ht> 450 && ht<= 550 && met> 210 && met<= 220) {eff2018 = 0.731707; errup = 0.0282492; errdown = 0.0282492;}
      else if (ht> 450 && ht<= 550 && met> 220 && met<= 230) {eff2018 = 0.797753; errup = 0.0301069; errdown = 0.0301069;}
      else if (ht> 450 && ht<= 550 && met> 230 && met<= 240) {eff2018 = 0.758865; errup = 0.0360249; errdown = 0.0360249;}
      else if (ht> 450 && ht<= 550 && met> 240 && met<= 250) {eff2018 = 0.824074; errup = 0.0366384; errdown = 0.0366384;}
      else if (ht> 450 && ht<= 550 && met> 250 && met<= 300) {eff2018 = 0.90458; errup = 0.0181507; errdown = 0.0181507;}
      else if (ht> 450 && ht<= 550 && met> 300 && met<= 400) {eff2018 = 0.948387; errup = 0.0177708; errdown = 0.0177708;}
      else if (ht> 450 && ht<= 550 && met> 400 && met<= 9999) {eff2018 = 0.959459; errup = 0.0229267; errdown = 0.0229267;}
      else if (ht> 550 && ht<= 650 && met> 150 && met<= 155) {eff2018 = 0.2286; errup = 0.00942298; errdown = 0.00942298;}
      else if (ht> 550 && ht<= 650 && met> 155 && met<= 160) {eff2018 = 0.262344; errup = 0.0107295; errdown = 0.0107295;}
      else if (ht> 550 && ht<= 650 && met> 160 && met<= 165) {eff2018 = 0.332286; errup = 0.0132019; errdown = 0.0132019;}
      else if (ht> 550 && ht<= 650 && met> 165 && met<= 170) {eff2018 = 0.356405; errup = 0.0153936; errdown = 0.0153936;}
      else if (ht> 550 && ht<= 650 && met> 170 && met<= 175) {eff2018 = 0.418685; errup = 0.0167548; errdown = 0.0167548;}
      else if (ht> 550 && ht<= 650 && met> 175 && met<= 180) {eff2018 = 0.471642; errup = 0.0192856; errdown = 0.0192856;}
      else if (ht> 550 && ht<= 650 && met> 180 && met<= 185) {eff2018 = 0.472222; errup = 0.0222374; errdown = 0.0222374;}
      else if (ht> 550 && ht<= 650 && met> 185 && met<= 190) {eff2018 = 0.543033; errup = 0.02255; errdown = 0.02255;}
      else if (ht> 550 && ht<= 650 && met> 190 && met<= 195) {eff2018 = 0.584699; errup = 0.0257577; errdown = 0.0257577;}
      else if (ht> 550 && ht<= 650 && met> 195 && met<= 200) {eff2018 = 0.635328; errup = 0.0256919; errdown = 0.0256919;}
      else if (ht> 550 && ht<= 650 && met> 200 && met<= 210) {eff2018 = 0.669922; errup = 0.0207819; errdown = 0.0207819;}
      else if (ht> 550 && ht<= 650 && met> 210 && met<= 220) {eff2018 = 0.737113; errup = 0.0223478; errdown = 0.0223478;}
      else if (ht> 550 && ht<= 650 && met> 220 && met<= 230) {eff2018 = 0.744409; errup = 0.0246551; errdown = 0.0246551;}
      else if (ht> 550 && ht<= 650 && met> 230 && met<= 240) {eff2018 = 0.791304; errup = 0.0267957; errdown = 0.0267957;}
      else if (ht> 550 && ht<= 650 && met> 240 && met<= 250) {eff2018 = 0.828125; errup = 0.0272272; errdown = 0.0272272;}
      else if (ht> 550 && ht<= 650 && met> 250 && met<= 275) {eff2018 = 0.881459; errup = 0.0178212; errdown = 0.0178212;}
      else if (ht> 550 && ht<= 650 && met> 275 && met<= 300) {eff2018 = 0.88587; errup = 0.023441; errdown = 0.023441;}
      else if (ht> 550 && ht<= 650 && met> 300 && met<= 400) {eff2018 = 0.887029; errup = 0.014479; errdown = 0.014479;}
      else if (ht> 550 && ht<= 650 && met> 400 && met<= 9999) {eff2018 = 0.977528; errup = 0.00473689; errdown = 0.00473689;}
      else if (ht> 650 && ht<= 800 && met> 150 && met<= 155) {eff2018 = 0.312199; errup = 0.00737961; errdown = 0.00737961;}
      else if (ht> 650 && ht<= 800 && met> 155 && met<= 160) {eff2018 = 0.355075; errup = 0.00835433; errdown = 0.00835433;}
      else if (ht> 650 && ht<= 800 && met> 160 && met<= 165) {eff2018 = 0.397912; errup = 0.00962519; errdown = 0.00962519;}
      else if (ht> 650 && ht<= 800 && met> 165 && met<= 170) {eff2018 = 0.446444; errup = 0.0104157; errdown = 0.0104157;}
      else if (ht> 650 && ht<= 800 && met> 170 && met<= 175) {eff2018 = 0.505056; errup = 0.0115341; errdown = 0.0115341;}
      else if (ht> 650 && ht<= 800 && met> 175 && met<= 180) {eff2018 = 0.540892; errup = 0.012404; errdown = 0.012404;}
      else if (ht> 650 && ht<= 800 && met> 180 && met<= 185) {eff2018 = 0.593127; errup = 0.0128787; errdown = 0.0128787;}
      else if (ht> 650 && ht<= 800 && met> 185 && met<= 190) {eff2018 = 0.609113; errup = 0.0137958; errdown = 0.0137958;}
      else if (ht> 650 && ht<= 800 && met> 190 && met<= 195) {eff2018 = 0.628141; errup = 0.0153217; errdown = 0.0153217;}
      else if (ht> 650 && ht<= 800 && met> 195 && met<= 200) {eff2018 = 0.692389; errup = 0.0150048; errdown = 0.0150048;}
      else if (ht> 650 && ht<= 800 && met> 200 && met<= 210) {eff2018 = 0.703607; errup = 0.011694; errdown = 0.011694;}
      else if (ht> 650 && ht<= 800 && met> 210 && met<= 220) {eff2018 = 0.790055; errup = 0.0114418; errdown = 0.0114418;}
      else if (ht> 650 && ht<= 800 && met> 220 && met<= 230) {eff2018 = 0.805634; errup = 0.0121256; errdown = 0.0121256;}
      else if (ht> 650 && ht<= 800 && met> 230 && met<= 240) {eff2018 = 0.862187; errup = 0.0116332; errdown = 0.0116332;}
      else if (ht> 650 && ht<= 800 && met> 240 && met<= 250) {eff2018 = 0.868639; errup = 0.0116205; errdown = 0.0116205;}
      else if (ht> 650 && ht<= 800 && met> 250 && met<= 275) {eff2018 = 0.878536; errup = 0.00769318; errdown = 0.00769318;}
      else if (ht> 650 && ht<= 800 && met> 275 && met<= 300) {eff2018 = 0.922261; errup = 0.00711815; errdown = 0.00711815;}
      else if (ht> 650 && ht<= 800 && met> 300 && met<= 350) {eff2018 = 0.934977; errup = 0.00508305; errdown = 0.00508305;}
      else if (ht> 650 && ht<= 800 && met> 350 && met<= 400) {eff2018 = 0.95888; errup = 0.00484746; errdown = 0.00484746;}
      else if (ht> 650 && ht<= 800 && met> 400 && met<= 450) {eff2018 = 0.969673; errup = 0.00484455; errdown = 0.00484455;}
      else if (ht> 650 && ht<= 800 && met> 450 && met<= 500) {eff2018 = 0.979543; errup = 0.00491061; errdown = 0.00491061;}
      else if (ht> 650 && ht<= 800 && met> 500 && met<= 9999) {eff2018 = 0.986716; errup = 0.00311023; errdown = 0.00311023;}
      else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff2018 = 0.539823; errup = 0.0034145; errdown = 0.0034145;}
      else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff2018 = 0.571429; errup = 0.00358377; errdown = 0.00358377;}
      else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff2018 = 0.61682; errup = 0.00383997; errdown = 0.00383997;}
      else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff2018 = 0.65604; errup = 0.00399833; errdown = 0.00399833;}
      else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff2018 = 0.675896; errup = 0.00424876; errdown = 0.00424876;}
      else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff2018 = 0.719656; errup = 0.00434348; errdown = 0.00434348;}
      else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff2018 = 0.735982; errup = 0.00455557; errdown = 0.00455557;}
      else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff2018 = 0.765251; errup = 0.00464919; errdown = 0.00464919;}
      else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff2018 = 0.773482; errup = 0.00482465; errdown = 0.00482465;}
      else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff2018 = 0.806076; errup = 0.00486078; errdown = 0.00486078;}
      else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff2018 = 0.833808; errup = 0.00351213; errdown = 0.00351213;}
      else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff2018 = 0.858087; errup = 0.00366153; errdown = 0.00366153;}
      else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff2018 = 0.874377; errup = 0.00379521; errdown = 0.00379521;}
      else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff2018 = 0.891151; errup = 0.00394335; errdown = 0.00394335;}
      else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff2018 = 0.909908; errup = 0.0039598; errdown = 0.0039598;}
      else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff2018 = 0.920343; errup = 0.00270517; errdown = 0.00270517;}
      else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff2018 = 0.940531; errup = 0.00282431; errdown = 0.00282431;}
      else if (ht> 800 && ht<= 1000 && met> 300 && met<= 350) {eff2018 = 0.943467; errup = 0.0025356; errdown = 0.0025356;}
      else if (ht> 800 && ht<= 1000 && met> 350 && met<= 400) {eff2018 = 0.952579; errup = 0.00332293; errdown = 0.00332293;}
      else if (ht> 800 && ht<= 1000 && met> 400 && met<= 450) {eff2018 = 0.968186; errup = 0.00382436; errdown = 0.00382436;}
      else if (ht> 800 && ht<= 1000 && met> 450 && met<= 500) {eff2018 = 0.969146; errup = 0.00486385; errdown = 0.00486385;}
      else if (ht> 800 && ht<= 1000 && met> 500 && met<= 9999) {eff2018 = 0.987481; errup = 0.00248803; errdown = 0.00248803;}
      else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff2018 = 0.432701; errup = 0.00130445; errdown = 0.00130445;}
      else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff2018 = 0.468416; errup = 0.00141951; errdown = 0.00141951;}
      else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff2018 = 0.500201; errup = 0.00153028; errdown = 0.00153028;}
      else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff2018 = 0.534339; errup = 0.00164201; errdown = 0.00164201;}
      else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff2018 = 0.565809; errup = 0.00176191; errdown = 0.00176191;}
      else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff2018 = 0.597802; errup = 0.00187686; errdown = 0.00187686;}
      else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff2018 = 0.624353; errup = 0.00198579; errdown = 0.00198579;}
      else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff2018 = 0.657567; errup = 0.00209913; errdown = 0.00209913;}
      else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff2018 = 0.683468; errup = 0.00220487; errdown = 0.00220487;}
      else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff2018 = 0.707921; errup = 0.00231018; errdown = 0.00231018;}
      else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff2018 = 0.736016; errup = 0.00174943; errdown = 0.00174943;}
      else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff2018 = 0.77612; errup = 0.0018808; errdown = 0.0018808;}
      else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff2018 = 0.81054; errup = 0.00199175; errdown = 0.00199175;}
      else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff2018 = 0.831745; errup = 0.00215975; errdown = 0.00215975;}
      else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff2018 = 0.855244; errup = 0.00226767; errdown = 0.00226767;}
      else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff2018 = 0.882622; errup = 0.00156784; errdown = 0.00156784;}
      else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff2018 = 0.90554; errup = 0.00178883; errdown = 0.00178883;}
      else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 350) {eff2018 = 0.927802; errup = 0.00151072; errdown = 0.00151072;}
      else if (ht> 1000 && ht<= 9999 && met> 350 && met<= 400) {eff2018 = 0.941449; errup = 0.00195205; errdown = 0.00195205;}
      else if (ht> 1000 && ht<= 9999 && met> 400 && met<= 450) {eff2018 = 0.944989; errup = 0.00260629; errdown = 0.00260629;}
      else if (ht> 1000 && ht<= 9999 && met> 450 && met<= 500) {eff2018 = 0.954409; errup = 0.00317334; errdown = 0.00317334;}
      else if (ht> 1000 && ht<= 9999 && met> 500 && met<= 9999) {eff2018 = 0.960595; errup = 0.00225241; errdown = 0.00225241;}
    }
  }
  else if (nel==1 && nmu==0) { // 1 electron CR
    if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 0 && met<= 50) {eff2016 = 0.10209; errup = 0.00591574; errdown = 0.00563482;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 0 && met<= 50) {eff2016 = 0.508521; errup = 0.0133826; errdown = 0.0133946;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 0 && met<= 50) {eff2016 = 0.540469; errup = 0.00506853; errdown = 0.00507683;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 0 && met<= 50) {eff2016 = 0.197796; errup = 0.00954788; errdown = 0.00922175;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 0 && met<= 50) {eff2016 = 0.522956; errup = 0.0172376; errdown = 0.0172907;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 0 && met<= 50) {eff2016 = 0.592268; errup = 0.00609976; errdown = 0.00612797;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 40 && met> 0 && met<= 50) {eff2016 = 0.343518; errup = 0.00991985; errdown = 0.0097879;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 40 && met> 0 && met<= 50) {eff2016 = 0.632466; errup = 0.0136146; errdown = 0.0138232;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 40 && met> 0 && met<= 50) {eff2016 = 0.676785; errup = 0.00477725; errdown = 0.00481403;}
    else if (ht> 0 && ht<= 400 && lep_pt> 40 && lep_pt<= 110 && met> 0 && met<= 50) {eff2016 = 0.51545; errup = 0.00513898; errdown = 0.00514221;}
    else if (ht> 400 && ht<= 600 && lep_pt> 40 && lep_pt<= 110 && met> 0 && met<= 50) {eff2016 = 0.734544; errup = 0.0056898; errdown = 0.00576763;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 40 && lep_pt<= 110 && met> 0 && met<= 50) {eff2016 = 0.756135; errup = 0.00212987; errdown = 0.00214247;}
    else if (ht> 0 && ht<= 400 && lep_pt> 110 && lep_pt<= 120 && met> 0 && met<= 50) {eff2016 = 0.76007; errup = 0.0128345; errdown = 0.0133055;}
    else if (ht> 400 && ht<= 600 && lep_pt> 110 && lep_pt<= 120 && met> 0 && met<= 50) {eff2016 = 0.85582; errup = 0.0129388; errdown = 0.0139281;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 110 && lep_pt<= 120 && met> 0 && met<= 50) {eff2016 = 0.860034; errup = 0.00541667; errdown = 0.00559399;}
    else if (ht> 0 && ht<= 400 && lep_pt> 120 && lep_pt<= 9999 && met> 0 && met<= 50) {eff2016 = 0.851656; errup = 0.00672062; errdown = 0.00697494;}
    else if (ht> 400 && ht<= 600 && lep_pt> 120 && lep_pt<= 9999 && met> 0 && met<= 50) {eff2016 = 0.884747; errup = 0.00665368; errdown = 0.00699398;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 120 && lep_pt<= 9999 && met> 0 && met<= 50) {eff2016 = 0.90492; errup = 0.00282916; errdown = 0.00290529;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 50 && met<= 75) {eff2016 = 0.207113; errup = 0.0202127; errdown = 0.0189198;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 50 && met<= 75) {eff2016 = 0.564516; errup = 0.0268312; errdown = 0.0271955;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 50 && met<= 75) {eff2016 = 0.564735; errup = 0.00789192; errdown = 0.00792435;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 50 && met<= 75) {eff2016 = 0.371345; errup = 0.0279443; errdown = 0.0271513;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 50 && met<= 75) {eff2016 = 0.664093; errup = 0.0305449; errdown = 0.0318921;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 50 && met<= 75) {eff2016 = 0.651934; errup = 0.00869838; errdown = 0.00879899;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 40 && met> 50 && met<= 75) {eff2016 = 0.523629; errup = 0.0225912; errdown = 0.0226843;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 40 && met> 50 && met<= 75) {eff2016 = 0.722334; errup = 0.0206066; errdown = 0.021545;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 40 && met> 50 && met<= 75) {eff2016 = 0.728115; errup = 0.00644023; errdown = 0.0065357;}
    else if (ht> 0 && ht<= 400 && lep_pt> 40 && lep_pt<= 110 && met> 50 && met<= 75) {eff2016 = 0.656417; errup = 0.0101746; errdown = 0.0103171;}
    else if (ht> 400 && ht<= 600 && lep_pt> 40 && lep_pt<= 110 && met> 50 && met<= 75) {eff2016 = 0.792651; errup = 0.00856531; errdown = 0.00882788;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 40 && lep_pt<= 110 && met> 50 && met<= 75) {eff2016 = 0.802124; errup = 0.00282306; errdown = 0.00285346;}
    else if (ht> 0 && ht<= 400 && lep_pt> 110 && lep_pt<= 120 && met> 50 && met<= 75) {eff2016 = 0.829493; errup = 0.0261432; errdown = 0.0294446;}
    else if (ht> 400 && ht<= 600 && lep_pt> 110 && lep_pt<= 120 && met> 50 && met<= 75) {eff2016 = 0.852459; errup = 0.0231571; errdown = 0.0262989;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 110 && lep_pt<= 120 && met> 50 && met<= 75) {eff2016 = 0.879075; errup = 0.00753224; errdown = 0.00794474;}
    else if (ht> 0 && ht<= 400 && lep_pt> 120 && lep_pt<= 9999 && met> 50 && met<= 75) {eff2016 = 0.866873; errup = 0.013536; errdown = 0.0147365;}
    else if (ht> 400 && ht<= 600 && lep_pt> 120 && lep_pt<= 9999 && met> 50 && met<= 75) {eff2016 = 0.90665; errup = 0.0104923; errdown = 0.0115955;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 120 && lep_pt<= 9999 && met> 50 && met<= 75) {eff2016 = 0.91716; errup = 0.00409646; errdown = 0.0042842;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 75 && met<= 100) {eff2016 = 0.283465; errup = 0.0456264; errdown = 0.0418996;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 75 && met<= 100) {eff2016 = 0.75; errup = 0.0330874; errdown = 0.036026;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 75 && met<= 100) {eff2016 = 0.664792; errup = 0.00995901; errdown = 0.0101047;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 75 && met<= 100) {eff2016 = 0.578512; errup = 0.0481009; errdown = 0.0495105;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 75 && met<= 100) {eff2016 = 0.741722; errup = 0.0370708; errdown = 0.0405552;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 75 && met<= 100) {eff2016 = 0.735936; errup = 0.0107074; errdown = 0.0109856;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 40 && met> 75 && met<= 100) {eff2016 = 0.703883; errup = 0.033095; errdown = 0.035219;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 40 && met> 75 && met<= 100) {eff2016 = 0.814545; errup = 0.0239949; errdown = 0.0264544;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 40 && met> 75 && met<= 100) {eff2016 = 0.78274; errup = 0.00762671; errdown = 0.00782077;}
    else if (ht> 0 && ht<= 400 && lep_pt> 40 && lep_pt<= 110 && met> 75 && met<= 100) {eff2016 = 0.784091; errup = 0.0140963; errdown = 0.014768;}
    else if (ht> 400 && ht<= 600 && lep_pt> 40 && lep_pt<= 110 && met> 75 && met<= 100) {eff2016 = 0.871389; errup = 0.0103228; errdown = 0.0110463;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 40 && lep_pt<= 110 && met> 75 && met<= 100) {eff2016 = 0.849313; errup = 0.00324759; errdown = 0.00330547;}
    else if (ht> 0 && ht<= 400 && lep_pt> 110 && lep_pt<= 120 && met> 75 && met<= 100) {eff2016 = 0.870968; errup = 0.0354651; errdown = 0.0446107;}
    else if (ht> 400 && ht<= 600 && lep_pt> 110 && lep_pt<= 120 && met> 75 && met<= 100) {eff2016 = 0.966942; errup = 0.015755; errdown = 0.0253664;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 110 && lep_pt<= 120 && met> 75 && met<= 100) {eff2016 = 0.91099; errup = 0.00864424; errdown = 0.00942993;}
    else if (ht> 0 && ht<= 400 && lep_pt> 120 && lep_pt<= 9999 && met> 75 && met<= 100) {eff2016 = 0.893333; errup = 0.0208625; errdown = 0.0247257;}
    else if (ht> 400 && ht<= 600 && lep_pt> 120 && lep_pt<= 9999 && met> 75 && met<= 100) {eff2016 = 0.926509; errup = 0.013463; errdown = 0.0159179;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 120 && lep_pt<= 9999 && met> 75 && met<= 100) {eff2016 = 0.944102; errup = 0.00452597; errdown = 0.00488294;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 100 && met<= 125) {eff2016 = 0.5; errup = 0.0780973; errdown = 0.0780973;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 100 && met<= 125) {eff2016 = 0.739726; errup = 0.0541361; errdown = 0.0615126;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 100 && met<= 125) {eff2016 = 0.746437; errup = 0.0120965; errdown = 0.0124777;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 100 && met<= 125) {eff2016 = 0.54; errup = 0.0786763; errdown = 0.080479;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 100 && met<= 125) {eff2016 = 0.794872; errup = 0.0475795; errdown = 0.0561273;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 100 && met<= 125) {eff2016 = 0.819343; errup = 0.0117716; errdown = 0.0123769;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 40 && met> 100 && met<= 125) {eff2016 = 0.843478; errup = 0.0347781; errdown = 0.0414801;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 40 && met> 100 && met<= 125) {eff2016 = 0.908397; errup = 0.0254458; errdown = 0.0325889;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 40 && met> 100 && met<= 125) {eff2016 = 0.842222; errup = 0.00867254; errdown = 0.00906496;}
    else if (ht> 0 && ht<= 400 && lep_pt> 40 && lep_pt<= 110 && met> 100 && met<= 125) {eff2016 = 0.864662; errup = 0.0173909; errdown = 0.0193496;}
    else if (ht> 400 && ht<= 600 && lep_pt> 40 && lep_pt<= 110 && met> 100 && met<= 125) {eff2016 = 0.926448; errup = 0.0103967; errdown = 0.0118313;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 40 && lep_pt<= 110 && met> 100 && met<= 125) {eff2016 = 0.89891; errup = 0.00351003; errdown = 0.00361948;}
    else if (ht> 0 && ht<= 400 && lep_pt> 110 && lep_pt<= 120 && met> 100 && met<= 125) {eff2016 = 0.966667; errup = 0.0275914; errdown = 0.072517;}
    else if (ht> 400 && ht<= 600 && lep_pt> 110 && lep_pt<= 120 && met> 100 && met<= 125) {eff2016 = 0.942529; errup = 0.0246014; errdown = 0.0370083;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 110 && lep_pt<= 120 && met> 100 && met<= 125) {eff2016 = 0.9377; errup = 0.00971033; errdown = 0.011223;}
    else if (ht> 0 && ht<= 400 && lep_pt> 120 && lep_pt<= 9999 && met> 100 && met<= 125) {eff2016 = 0.966292; errup = 0.0182896; errdown = 0.0317032;}
    else if (ht> 400 && ht<= 600 && lep_pt> 120 && lep_pt<= 9999 && met> 100 && met<= 125) {eff2016 = 0.950495; errup = 0.0152082; errdown = 0.0203749;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 120 && lep_pt<= 9999 && met> 100 && met<= 125) {eff2016 = 0.965256; errup = 0.00460922; errdown = 0.0052372;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 125 && met<= 150) {eff2016 = 0.851852; errup = 0.0695101; errdown = 0.101731;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 125 && met<= 150) {eff2016 = 0.936508; errup = 0.0301398; errdown = 0.0473598;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 125 && met<= 150) {eff2016 = 0.853598; errup = 0.0126086; errdown = 0.013529;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 125 && met<= 150) {eff2016 = 0.909091; errup = 0.0490632; errdown = 0.0805904;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 125 && met<= 150) {eff2016 = 1; errup = 0; errdown = 0.0472931;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 125 && met<= 150) {eff2016 = 0.871069; errup = 0.0134541; errdown = 0.0146888;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 40 && met> 125 && met<= 150) {eff2016 = 0.883333; errup = 0.0420083; errdown = 0.0571758;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 40 && met> 125 && met<= 150) {eff2016 = 0.929825; errup = 0.0332829; errdown = 0.0520157;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 40 && met> 125 && met<= 150) {eff2016 = 0.876818; errup = 0.00970247; errdown = 0.0103755;}
    else if (ht> 0 && ht<= 400 && lep_pt> 40 && lep_pt<= 110 && met> 125 && met<= 150) {eff2016 = 0.920354; errup = 0.0181519; errdown = 0.0223173;}
    else if (ht> 400 && ht<= 600 && lep_pt> 40 && lep_pt<= 110 && met> 125 && met<= 150) {eff2016 = 0.948276; errup = 0.0103093; errdown = 0.0124414;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 40 && lep_pt<= 110 && met> 125 && met<= 150) {eff2016 = 0.935926; errup = 0.00348682; errdown = 0.00366757;}
    else if (ht> 0 && ht<= 400 && lep_pt> 110 && lep_pt<= 120 && met> 125 && met<= 150) {eff2016 = 0.952381; errup = 0.0394264; errdown = 0.101134;}
    else if (ht> 400 && ht<= 600 && lep_pt> 110 && lep_pt<= 120 && met> 125 && met<= 150) {eff2016 = 0.981481; errup = 0.0153245; errdown = 0.0412969;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 110 && lep_pt<= 120 && met> 125 && met<= 150) {eff2016 = 0.971074; errup = 0.00757748; errdown = 0.00978541;}
    else if (ht> 0 && ht<= 400 && lep_pt> 120 && lep_pt<= 9999 && met> 125 && met<= 150) {eff2016 = 0.981481; errup = 0.0153245; errdown = 0.0412969;}
    else if (ht> 400 && ht<= 600 && lep_pt> 120 && lep_pt<= 9999 && met> 125 && met<= 150) {eff2016 = 0.951613; errup = 0.0190274; errdown = 0.0277831;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 120 && lep_pt<= 9999 && met> 125 && met<= 150) {eff2016 = 0.97007; errup = 0.00505453; errdown = 0.00595382;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 150 && met<= 175) {eff2016 = 0.857143; errup = 0.0917089; errdown = 0.158271;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 150 && met<= 175) {eff2016 = 1; errup = 0; errdown = 0.0595223;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 150 && met<= 175) {eff2016 = 0.925703; errup = 0.0118357; errdown = 0.0136893;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 150 && met<= 175) {eff2016 = 1; errup = 0; errdown = 0.205568;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 150 && met<= 175) {eff2016 = 1; errup = 0; errdown = 0.0709947;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 150 && met<= 175) {eff2016 = 0.932367; errup = 0.0124155; errdown = 0.0147054;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 40 && met> 150 && met<= 175) {eff2016 = 0.810811; errup = 0.0670146; errdown = 0.0869581;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 40 && met> 150 && met<= 175) {eff2016 = 0.980769; errup = 0.0159141; errdown = 0.0428344;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 40 && met> 150 && met<= 175) {eff2016 = 0.950556; errup = 0.00764631; errdown = 0.00884972;}
    else if (ht> 0 && ht<= 400 && lep_pt> 40 && lep_pt<= 110 && met> 150 && met<= 175) {eff2016 = 0.972603; errup = 0.0130667; errdown = 0.0211322;}
    else if (ht> 400 && ht<= 600 && lep_pt> 40 && lep_pt<= 110 && met> 150 && met<= 175) {eff2016 = 0.977564; errup = 0.00823774; errdown = 0.0118755;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 40 && lep_pt<= 110 && met> 150 && met<= 175) {eff2016 = 0.960628; errup = 0.00329173; errdown = 0.00356596;}
    else if (ht> 0 && ht<= 400 && lep_pt> 110 && lep_pt<= 120 && met> 150 && met<= 175) {eff2016 = 1; errup = 0; errdown = 0.115502;}
    else if (ht> 400 && ht<= 600 && lep_pt> 110 && lep_pt<= 120 && met> 150 && met<= 175) {eff2016 = 1; errup = 0; errdown = 0.0409781;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 110 && lep_pt<= 120 && met> 150 && met<= 175) {eff2016 = 0.975385; errup = 0.00847391; errdown = 0.0119156;}
    else if (ht> 0 && ht<= 400 && lep_pt> 120 && lep_pt<= 9999 && met> 150 && met<= 175) {eff2016 = 1; errup = 0; errdown = 0.0472931;}
    else if (ht> 400 && ht<= 600 && lep_pt> 120 && lep_pt<= 9999 && met> 150 && met<= 175) {eff2016 = 0.990654; errup = 0.00773258; errdown = 0.0211614;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 120 && lep_pt<= 9999 && met> 150 && met<= 175) {eff2016 = 0.988032; errup = 0.00390372; errdown = 0.00541802;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 175 && met<= 215) {eff2016 = 1; errup = 0; errdown = 0.123222;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 175 && met<= 215) {eff2016 = 1; errup = 0; errdown = 0.0709947;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 175 && met<= 215) {eff2016 = 0.95572; errup = 0.00884752; errdown = 0.0107028;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 175 && met<= 215) {eff2016 = 1; errup = 0; errdown = 0.115502;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 175 && met<= 215) {eff2016 = 0.96875; errup = 0.025866; errdown = 0.0682225;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 175 && met<= 215) {eff2016 = 0.974895; errup = 0.00709912; errdown = 0.00937376;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 40 && met> 175 && met<= 215) {eff2016 = 1; errup = 0; errdown = 0.102638;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 40 && met> 175 && met<= 215) {eff2016 = 1; errup = 0; errdown = 0.0376284;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 40 && met> 175 && met<= 215) {eff2016 = 0.982256; errup = 0.00466259; errdown = 0.00604897;}
    else if (ht> 0 && ht<= 400 && lep_pt> 40 && lep_pt<= 110 && met> 175 && met<= 215) {eff2016 = 0.993333; errup = 0.00551564; errdown = 0.0151622;}
    else if (ht> 400 && ht<= 600 && lep_pt> 40 && lep_pt<= 110 && met> 175 && met<= 215) {eff2016 = 0.990291; errup = 0.00527929; errdown = 0.00935364;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 40 && lep_pt<= 110 && met> 175 && met<= 215) {eff2016 = 0.981423; errup = 0.00219908; errdown = 0.00247001;}
    else if (ht> 0 && ht<= 400 && lep_pt> 110 && lep_pt<= 120 && met> 175 && met<= 215) {eff2016 = 1; errup = 0; errdown = 0.205568;}
    else if (ht> 400 && ht<= 600 && lep_pt> 110 && lep_pt<= 120 && met> 175 && met<= 215) {eff2016 = 1; errup = 0; errdown = 0.0419109;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 110 && lep_pt<= 120 && met> 175 && met<= 215) {eff2016 = 0.994505; errup = 0.00354816; errdown = 0.00720076;}
    else if (ht> 0 && ht<= 400 && lep_pt> 120 && lep_pt<= 9999 && met> 175 && met<= 215) {eff2016 = 1; errup = 0; errdown = 0.0879414;}
    else if (ht> 400 && ht<= 600 && lep_pt> 120 && lep_pt<= 9999 && met> 175 && met<= 215) {eff2016 = 1; errup = 0; errdown = 0.0180628;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 120 && lep_pt<= 9999 && met> 175 && met<= 215) {eff2016 = 0.992408; errup = 0.00279574; errdown = 0.00406535;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 215 && met<= 9999) {eff2016 = 1; errup = 0; errdown = 0.168149;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 215 && met<= 9999) {eff2016 = 1; errup = 0; errdown = 0.108691;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 215 && met<= 9999) {eff2016 = 0.991718; errup = 0.00395933; errdown = 0.00649962;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 215 && met<= 9999) {eff2016 = 1; errup = 0; errdown = 0.264229;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 215 && met<= 9999) {eff2016 = 1; errup = 0; errdown = 0.0839348;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 215 && met<= 9999) {eff2016 = 0.986143; errup = 0.00548298; errdown = 0.00818471;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 40 && met> 215 && met<= 9999) {eff2016 = 1; errup = 0; errdown = 0.115502;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 40 && met> 215 && met<= 9999) {eff2016 = 1; errup = 0; errdown = 0.0559083;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 40 && met> 215 && met<= 9999) {eff2016 = 0.99409; errup = 0.00255052; errdown = 0.00397839;}
    else if (ht> 0 && ht<= 400 && lep_pt> 40 && lep_pt<= 110 && met> 215 && met<= 9999) {eff2016 = 1; errup = 0; errdown = 0.0173807;}
    else if (ht> 400 && ht<= 600 && lep_pt> 40 && lep_pt<= 110 && met> 215 && met<= 9999) {eff2016 = 0.996441; errup = 0.00294413; errdown = 0.00813542;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 40 && lep_pt<= 110 && met> 215 && met<= 9999) {eff2016 = 0.994857; errup = 0.0011116; errdown = 0.00138053;}
    else if (ht> 0 && ht<= 400 && lep_pt> 110 && lep_pt<= 120 && met> 215 && met<= 9999) {eff2016 = 1; errup = 0; errdown = 0.458642;}
    else if (ht> 400 && ht<= 600 && lep_pt> 110 && lep_pt<= 120 && met> 215 && met<= 9999) {eff2016 = 1; errup = 0; errdown = 0.0498539;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 110 && lep_pt<= 120 && met> 215 && met<= 9999) {eff2016 = 0.995455; errup = 0.00293541; errdown = 0.00596358;}
    else if (ht> 0 && ht<= 400 && lep_pt> 120 && lep_pt<= 9999 && met> 215 && met<= 9999) {eff2016 = 1; errup = 0; errdown = 0.142229;}
    else if (ht> 400 && ht<= 600 && lep_pt> 120 && lep_pt<= 9999 && met> 215 && met<= 9999) {eff2016 = 1; errup = 0; errdown = 0.0149771;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 120 && lep_pt<= 9999 && met> 215 && met<= 9999) {eff2016 = 0.998031; errup = 0.00127137; errdown = 0.0025904;}
    //2017
    if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 0 && met<= 50) {eff2017 = 0.00606061; errup = 0.00297522; errdown = 0.0020946;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 0 && met<= 50) {eff2017 = 0.252727; errup = 0.0137849; errdown = 0.0133199;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 0 && met<= 50) {eff2017 = 0.47488; errup = 0.0101789; errdown = 0.0101585;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 0 && met<= 50) {eff2017 = 0.0580311; errup = 0.00853854; errdown = 0.00756;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 0 && met<= 50) {eff2017 = 0.273902; errup = 0.0169699; errdown = 0.0163626;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 0 && met<= 50) {eff2017 = 0.512061; errup = 0.0119679; errdown = 0.0119814;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 40 && met> 0 && met<= 50) {eff2017 = 0.239788; errup = 0.0123203; errdown = 0.0119141;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 40 && met> 0 && met<= 50) {eff2017 = 0.432836; errup = 0.0161873; errdown = 0.0160495;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 40 && met> 0 && met<= 50) {eff2017 = 0.610241; errup = 0.00967965; errdown = 0.00976557;}
    else if (ht> 0 && ht<= 400 && lep_pt> 40 && lep_pt<= 110 && met> 0 && met<= 50) {eff2017 = 0.368855; errup = 0.0057901; errdown = 0.00575302;}
    else if (ht> 400 && ht<= 600 && lep_pt> 40 && lep_pt<= 110 && met> 0 && met<= 50) {eff2017 = 0.58404; errup = 0.00788592; errdown = 0.0079285;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 40 && lep_pt<= 110 && met> 0 && met<= 50) {eff2017 = 0.718103; errup = 0.00425673; errdown = 0.00429573;}
    else if (ht> 0 && ht<= 400 && lep_pt> 110 && lep_pt<= 120 && met> 0 && met<= 50) {eff2017 = 0.358079; errup = 0.0165383; errdown = 0.0162179;}
    else if (ht> 400 && ht<= 600 && lep_pt> 110 && lep_pt<= 120 && met> 0 && met<= 50) {eff2017 = 0.574176; errup = 0.0270377; errdown = 0.0274661;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 110 && lep_pt<= 120 && met> 0 && met<= 50) {eff2017 = 0.721596; errup = 0.01343; errdown = 0.0138268;}
    else if (ht> 0 && ht<= 400 && lep_pt> 120 && lep_pt<= 9999 && met> 0 && met<= 50) {eff2017 = 0.727555; errup = 0.00887959; errdown = 0.00906034;}
    else if (ht> 400 && ht<= 600 && lep_pt> 120 && lep_pt<= 9999 && met> 0 && met<= 50) {eff2017 = 0.776404; errup = 0.014199; errdown = 0.0148447;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 120 && lep_pt<= 9999 && met> 0 && met<= 50) {eff2017 = 0.830705; errup = 0.00707675; errdown = 0.00731432;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 50 && met<= 75) {eff2017 = 0.0128205; errup = 0.00857983; errdown = 0.00552679;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 50 && met<= 75) {eff2017 = 0.251185; errup = 0.0228961; errdown = 0.0216526;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 50 && met<= 75) {eff2017 = 0.515088; errup = 0.0135744; errdown = 0.0135961;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 50 && met<= 75) {eff2017 = 0.07; errup = 0.0180052; errdown = 0.0148206;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 50 && met<= 75) {eff2017 = 0.359756; errup = 0.0284325; errdown = 0.0275299;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 50 && met<= 75) {eff2017 = 0.58209; errup = 0.0154441; errdown = 0.0156019;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 40 && met> 50 && met<= 75) {eff2017 = 0.313953; errup = 0.0239721; errdown = 0.0230631;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 40 && met> 50 && met<= 75) {eff2017 = 0.564612; errup = 0.0229432; errdown = 0.0232113;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 40 && met> 50 && met<= 75) {eff2017 = 0.668443; errup = 0.0116455; errdown = 0.0118501;}
    else if (ht> 0 && ht<= 400 && lep_pt> 40 && lep_pt<= 110 && met> 50 && met<= 75) {eff2017 = 0.478528; errup = 0.0103138; errdown = 0.0102959;}
    else if (ht> 400 && ht<= 600 && lep_pt> 40 && lep_pt<= 110 && met> 50 && met<= 75) {eff2017 = 0.690422; errup = 0.0112097; errdown = 0.0114324;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 40 && lep_pt<= 110 && met> 50 && met<= 75) {eff2017 = 0.76283; errup = 0.00509721; errdown = 0.00517276;}
    else if (ht> 0 && ht<= 400 && lep_pt> 110 && lep_pt<= 120 && met> 50 && met<= 75) {eff2017 = 0.404682; errup = 0.0303426; errdown = 0.0296689;}
    else if (ht> 400 && ht<= 600 && lep_pt> 110 && lep_pt<= 120 && met> 50 && met<= 75) {eff2017 = 0.722973; errup = 0.0383999; errdown = 0.0416756;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 110 && lep_pt<= 120 && met> 50 && met<= 75) {eff2017 = 0.766862; errup = 0.0165056; errdown = 0.0173229;}
    else if (ht> 0 && ht<= 400 && lep_pt> 120 && lep_pt<= 9999 && met> 50 && met<= 75) {eff2017 = 0.73375; errup = 0.0159396; errdown = 0.0165468;}
    else if (ht> 400 && ht<= 600 && lep_pt> 120 && lep_pt<= 9999 && met> 50 && met<= 75) {eff2017 = 0.814136; errup = 0.0203184; errdown = 0.0220684;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 120 && lep_pt<= 9999 && met> 50 && met<= 75) {eff2017 = 0.843819; errup = 0.00860726; errdown = 0.0089989;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 75 && met<= 100) {eff2017 = 0.0526316; errup = 0.0271984; errdown = 0.019208;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 75 && met<= 100) {eff2017 = 0.295082; errup = 0.0375578; errdown = 0.0351452;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 75 && met<= 100) {eff2017 = 0.577465; errup = 0.0167048; errdown = 0.0168781;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 75 && met<= 100) {eff2017 = 0.196429; errup = 0.0448753; errdown = 0.0388449;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 75 && met<= 100) {eff2017 = 0.467105; errup = 0.0438353; errdown = 0.0433692;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 75 && met<= 100) {eff2017 = 0.620134; errup = 0.0182735; errdown = 0.018608;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 40 && met> 75 && met<= 100) {eff2017 = 0.495098; errup = 0.037384; errdown = 0.0373327;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 40 && met> 75 && met<= 100) {eff2017 = 0.59375; errup = 0.03452; errdown = 0.0354114;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 40 && met> 75 && met<= 100) {eff2017 = 0.738619; errup = 0.0125052; errdown = 0.0128916;}
    else if (ht> 0 && ht<= 400 && lep_pt> 40 && lep_pt<= 110 && met> 75 && met<= 100) {eff2017 = 0.584366; errup = 0.0156025; errdown = 0.0157683;}
    else if (ht> 400 && ht<= 600 && lep_pt> 40 && lep_pt<= 110 && met> 75 && met<= 100) {eff2017 = 0.737811; errup = 0.0147447; errdown = 0.0152789;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 40 && lep_pt<= 110 && met> 75 && met<= 100) {eff2017 = 0.796101; errup = 0.00566351; errdown = 0.00578094;}
    else if (ht> 0 && ht<= 400 && lep_pt> 110 && lep_pt<= 120 && met> 75 && met<= 100) {eff2017 = 0.607843; errup = 0.0518085; errdown = 0.0541227;}
    else if (ht> 400 && ht<= 600 && lep_pt> 110 && lep_pt<= 120 && met> 75 && met<= 100) {eff2017 = 0.77551; errup = 0.0439117; errdown = 0.050173;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 110 && lep_pt<= 120 && met> 75 && met<= 100) {eff2017 = 0.801205; errup = 0.0182366; errdown = 0.0195126;}
    else if (ht> 0 && ht<= 400 && lep_pt> 120 && lep_pt<= 9999 && met> 75 && met<= 100) {eff2017 = 0.767516; errup = 0.0244905; errdown = 0.0263047;}
    else if (ht> 400 && ht<= 600 && lep_pt> 120 && lep_pt<= 9999 && met> 75 && met<= 100) {eff2017 = 0.791878; errup = 0.0297937; errdown = 0.0330066;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 120 && lep_pt<= 9999 && met> 75 && met<= 100) {eff2017 = 0.863636; errup = 0.00987432; errdown = 0.0104889;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 100 && met<= 125) {eff2017 = 0.0338983; errup = 0.0429671; errdown = 0.021865;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 100 && met<= 125) {eff2017 = 0.357895; errup = 0.0557905; errdown = 0.0525014;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 100 && met<= 125) {eff2017 = 0.63606; errup = 0.0202397; errdown = 0.0207128;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 100 && met<= 125) {eff2017 = 0.230769; errup = 0.0733859; errdown = 0.0614831;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 100 && met<= 125) {eff2017 = 0.532468; errup = 0.0624335; errdown = 0.0633651;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 100 && met<= 125) {eff2017 = 0.701613; errup = 0.0211067; errdown = 0.0219585;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 40 && met> 100 && met<= 125) {eff2017 = 0.645161; errup = 0.0455298; errdown = 0.0480775;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 40 && met> 100 && met<= 125) {eff2017 = 0.700787; errup = 0.0426733; errdown = 0.0461241;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 40 && met> 100 && met<= 125) {eff2017 = 0.807339; errup = 0.0135565; errdown = 0.0142912;}
    else if (ht> 0 && ht<= 400 && lep_pt> 40 && lep_pt<= 110 && met> 100 && met<= 125) {eff2017 = 0.716075; errup = 0.0211522; errdown = 0.0220989;}
    else if (ht> 400 && ht<= 600 && lep_pt> 40 && lep_pt<= 110 && met> 100 && met<= 125) {eff2017 = 0.797227; errup = 0.0170527; errdown = 0.0181352;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 40 && lep_pt<= 110 && met> 100 && met<= 125) {eff2017 = 0.84023; errup = 0.0061067; errdown = 0.00629736;}
    else if (ht> 0 && ht<= 400 && lep_pt> 110 && lep_pt<= 120 && met> 100 && met<= 125) {eff2017 = 0.785714; errup = 0.0664487; errdown = 0.0823485;}
    else if (ht> 400 && ht<= 600 && lep_pt> 110 && lep_pt<= 120 && met> 100 && met<= 125) {eff2017 = 0.891892; errup = 0.0510106; errdown = 0.077259;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 110 && lep_pt<= 120 && met> 100 && met<= 125) {eff2017 = 0.827922; errup = 0.0219648; errdown = 0.0242503;}
    else if (ht> 0 && ht<= 400 && lep_pt> 120 && lep_pt<= 9999 && met> 100 && met<= 125) {eff2017 = 0.841667; errup = 0.0342095; errdown = 0.0405813;}
    else if (ht> 400 && ht<= 600 && lep_pt> 120 && lep_pt<= 9999 && met> 100 && met<= 125) {eff2017 = 0.893443; errup = 0.0283253; errdown = 0.0356815;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 120 && lep_pt<= 9999 && met> 100 && met<= 125) {eff2017 = 0.905502; errup = 0.0102026; errdown = 0.0112294;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 125 && met<= 150) {eff2017 = 0.2; errup = 0.130119; errdown = 0.0931225;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 125 && met<= 150) {eff2017 = 0.604167; errup = 0.0776099; errdown = 0.0825214;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 125 && met<= 150) {eff2017 = 0.801508; errup = 0.0204268; errdown = 0.0220347;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 125 && met<= 150) {eff2017 = 0.222222; errup = 0.109068; errdown = 0.0843892;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 125 && met<= 150) {eff2017 = 0.575; errup = 0.0872414; errdown = 0.0915181;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 125 && met<= 150) {eff2017 = 0.788406; errup = 0.0225238; errdown = 0.0243027;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 40 && met> 125 && met<= 150) {eff2017 = 0.745098; errup = 0.0646476; errdown = 0.0756564;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 40 && met> 125 && met<= 150) {eff2017 = 0.813333; errup = 0.0466366; errdown = 0.0561452;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 40 && met> 125 && met<= 150) {eff2017 = 0.85; errup = 0.0150446; errdown = 0.0163204;}
    else if (ht> 0 && ht<= 400 && lep_pt> 40 && lep_pt<= 110 && met> 125 && met<= 150) {eff2017 = 0.815686; errup = 0.0248726; errdown = 0.027542;}
    else if (ht> 400 && ht<= 600 && lep_pt> 40 && lep_pt<= 110 && met> 125 && met<= 150) {eff2017 = 0.876254; errup = 0.0193393; errdown = 0.0220715;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 40 && lep_pt<= 110 && met> 125 && met<= 150) {eff2017 = 0.895428; errup = 0.00633502; errdown = 0.00668102;}
    else if (ht> 0 && ht<= 400 && lep_pt> 110 && lep_pt<= 120 && met> 125 && met<= 150) {eff2017 = 0.666667; errup = 0.106107; errdown = 0.122517;}
    else if (ht> 400 && ht<= 600 && lep_pt> 110 && lep_pt<= 120 && met> 125 && met<= 150) {eff2017 = 0.892857; errup = 0.0577336; errdown = 0.0933526;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 110 && lep_pt<= 120 && met> 125 && met<= 150) {eff2017 = 0.885965; errup = 0.0213678; errdown = 0.0250959;}
    else if (ht> 0 && ht<= 400 && lep_pt> 120 && lep_pt<= 9999 && met> 125 && met<= 150) {eff2017 = 0.982456; errup = 0.0145177; errdown = 0.0391869;}
    else if (ht> 400 && ht<= 600 && lep_pt> 120 && lep_pt<= 9999 && met> 125 && met<= 150) {eff2017 = 0.970588; errup = 0.0189746; errdown = 0.0374789;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 120 && lep_pt<= 9999 && met> 125 && met<= 150) {eff2017 = 0.93956; errup = 0.0102474; errdown = 0.0120003;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 150 && met<= 175) {eff2017 = 0.461538; errup = 0.172004; errdown = 0.164847;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 150 && met<= 175) {eff2017 = 0.612903; errup = 0.0974957; errdown = 0.105934;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 150 && met<= 175) {eff2017 = 0.865546; errup = 0.0225156; errdown = 0.0258705;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 150 && met<= 175) {eff2017 = 0.454545; errup = 0.189662; errdown = 0.179582;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 150 && met<= 175) {eff2017 = 0.666667; errup = 0.106107; errdown = 0.122517;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 150 && met<= 175) {eff2017 = 0.883178; errup = 0.0223037; errdown = 0.0262543;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 40 && met> 150 && met<= 175) {eff2017 = 0.75; errup = 0.0943497; errdown = 0.119341;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 40 && met> 150 && met<= 175) {eff2017 = 0.906977; errup = 0.0439838; errdown = 0.0674608;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 40 && met> 150 && met<= 175) {eff2017 = 0.914286; errup = 0.0132412; errdown = 0.0152122;}
    else if (ht> 0 && ht<= 400 && lep_pt> 40 && lep_pt<= 110 && met> 150 && met<= 175) {eff2017 = 0.873418; errup = 0.026946; errdown = 0.0322169;}
    else if (ht> 400 && ht<= 600 && lep_pt> 40 && lep_pt<= 110 && met> 150 && met<= 175) {eff2017 = 0.929648; errup = 0.0182155; errdown = 0.0231163;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 40 && lep_pt<= 110 && met> 150 && met<= 175) {eff2017 = 0.929748; errup = 0.00601652; errdown = 0.00650973;}
    else if (ht> 0 && ht<= 400 && lep_pt> 110 && lep_pt<= 120 && met> 150 && met<= 175) {eff2017 = 0.809524; errup = 0.0888157; errdown = 0.125184;}
    else if (ht> 400 && ht<= 600 && lep_pt> 110 && lep_pt<= 120 && met> 150 && met<= 175) {eff2017 = 1; errup = 0; errdown = 0.0709947;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 110 && lep_pt<= 120 && met> 150 && met<= 175) {eff2017 = 0.907801; errup = 0.0246102; errdown = 0.0312058;}
    else if (ht> 0 && ht<= 400 && lep_pt> 120 && lep_pt<= 9999 && met> 150 && met<= 175) {eff2017 = 0.942857; errup = 0.0368224; errdown = 0.0704444;}
    else if (ht> 400 && ht<= 600 && lep_pt> 120 && lep_pt<= 9999 && met> 150 && met<= 175) {eff2017 = 0.975; errup = 0.0206905; errdown = 0.0551516;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 120 && lep_pt<= 9999 && met> 150 && met<= 175) {eff2017 = 0.951351; errup = 0.011195; errdown = 0.0139251;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 175 && met<= 215) {eff2017 = 0.8; errup = 0.106751; errdown = 0.157061;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 175 && met<= 215) {eff2017 = 0.916667; errup = 0.0536391; errdown = 0.0995072;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 175 && met<= 215) {eff2017 = 0.955645; errup = 0.0130229; errdown = 0.017253;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 175 && met<= 215) {eff2017 = 0.769231; errup = 0.122762; errdown = 0.174724;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 175 && met<= 215) {eff2017 = 0.869565; errup = 0.0701228; errdown = 0.110814;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 175 && met<= 215) {eff2017 = 0.920833; errup = 0.0175659; errdown = 0.0214819;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 40 && met> 175 && met<= 215) {eff2017 = 0.944444; errup = 0.046004; errdown = 0.116415;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 40 && met> 175 && met<= 215) {eff2017 = 0.97619; errup = 0.0197048; errdown = 0.0526298;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 40 && met> 175 && met<= 215) {eff2017 = 0.961039; errup = 0.00899222; errdown = 0.0112252;}
    else if (ht> 0 && ht<= 400 && lep_pt> 40 && lep_pt<= 110 && met> 175 && met<= 215) {eff2017 = 0.950413; errup = 0.0194949; errdown = 0.0284436;}
    else if (ht> 400 && ht<= 600 && lep_pt> 40 && lep_pt<= 110 && met> 175 && met<= 215) {eff2017 = 0.978022; errup = 0.0104893; errdown = 0.0170363;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 40 && lep_pt<= 110 && met> 175 && met<= 215) {eff2017 = 0.969503; errup = 0.00381731; errdown = 0.00430889;}
    else if (ht> 0 && ht<= 400 && lep_pt> 110 && lep_pt<= 120 && met> 175 && met<= 215) {eff2017 = 0.916667; errup = 0.0690403; errdown = 0.16652;}
    else if (ht> 400 && ht<= 600 && lep_pt> 110 && lep_pt<= 120 && met> 175 && met<= 215) {eff2017 = 1; errup = 0; errdown = 0.0802771;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 110 && lep_pt<= 120 && met> 175 && met<= 215) {eff2017 = 0.964706; errup = 0.013912; errdown = 0.0204859;}
    else if (ht> 0 && ht<= 400 && lep_pt> 120 && lep_pt<= 9999 && met> 175 && met<= 215) {eff2017 = 0.966667; errup = 0.0275914; errdown = 0.072517;}
    else if (ht> 400 && ht<= 600 && lep_pt> 120 && lep_pt<= 9999 && met> 175 && met<= 215) {eff2017 = 1; errup = 0; errdown = 0.0384134;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 120 && lep_pt<= 9999 && met> 175 && met<= 215) {eff2017 = 0.980815; errup = 0.00661202; errdown = 0.00932515;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 215 && met<= 9999) {eff2017 = 1; errup = 0; errdown = 0.205568;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 215 && met<= 9999) {eff2017 = 0.913043; errup = 0.0559625; errdown = 0.103371;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 215 && met<= 9999) {eff2017 = 0.981651; errup = 0.00725455; errdown = 0.0107985;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 215 && met<= 9999) {eff2017 = 1; errup = 0; errdown = 0.205568;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 215 && met<= 9999) {eff2017 = 1; errup = 0; errdown = 0.0769247;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 215 && met<= 9999) {eff2017 = 0.984674; errup = 0.00732055; errdown = 0.0119517;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 40 && met> 215 && met<= 9999) {eff2017 = 1; errup = 0; errdown = 0.0659133;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 40 && met> 215 && met<= 9999) {eff2017 = 1; errup = 0; errdown = 0.042887;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 40 && met> 215 && met<= 9999) {eff2017 = 0.997938; errup = 0.00170573; errdown = 0.00472518;}
    else if (ht> 0 && ht<= 400 && lep_pt> 40 && lep_pt<= 110 && met> 215 && met<= 9999) {eff2017 = 0.990099; errup = 0.00819202; errdown = 0.0223979;}
    else if (ht> 400 && ht<= 600 && lep_pt> 40 && lep_pt<= 110 && met> 215 && met<= 9999) {eff2017 = 1; errup = 0; errdown = 0.00964279;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 40 && lep_pt<= 110 && met> 215 && met<= 9999) {eff2017 = 0.992439; errup = 0.00171528; errdown = 0.00215231;}
    else if (ht> 0 && ht<= 400 && lep_pt> 110 && lep_pt<= 120 && met> 215 && met<= 9999) {eff2017 = 1; errup = 0; errdown = 0.205568;}
    else if (ht> 400 && ht<= 600 && lep_pt> 110 && lep_pt<= 120 && met> 215 && met<= 9999) {eff2017 = 1; errup = 0; errdown = 0.115502;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 110 && lep_pt<= 120 && met> 215 && met<= 9999) {eff2017 = 1; errup = 0; errdown = 0.00804214;}
    else if (ht> 0 && ht<= 400 && lep_pt> 120 && lep_pt<= 9999 && met> 215 && met<= 9999) {eff2017 = 1; errup = 0; errdown = 0.108691;}
    else if (ht> 400 && ht<= 600 && lep_pt> 120 && lep_pt<= 9999 && met> 215 && met<= 9999) {eff2017 = 1; errup = 0; errdown = 0.0283562;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 120 && lep_pt<= 9999 && met> 215 && met<= 9999) {eff2017 = 0.998397; errup = 0.00132575; errdown = 0.0036754;}
    //2018
    if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 0 && met<= 50) {eff2018 = 0.0097629; errup = 0.00521891; errdown = 0.00359356;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 0 && met<= 50) {eff2018 = 0.238559; errup = 0.0140446; errdown = 0.0135171;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 0 && met<= 50) {eff2018 = 0.523018; errup = 0.00945585; errdown = 0.00947208;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 0 && met<= 50) {eff2018 = 0.017301; errup = 0.00729105; errdown = 0.00535648;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 0 && met<= 50) {eff2018 = 0.297573; errup = 0.0172358; errdown = 0.0166989;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 0 && met<= 50) {eff2018 = 0.583857; errup = 0.0110342; errdown = 0.011117;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 40 && met> 0 && met<= 50) {eff2018 = 0.260109; errup = 0.0153141; errdown = 0.0147704;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 40 && met> 0 && met<= 50) {eff2018 = 0.444776; errup = 0.0162215; errdown = 0.0161081;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 40 && met> 0 && met<= 50) {eff2018 = 0.650049; errup = 0.00869262; errdown = 0.00879159;}
    else if (ht> 0 && ht<= 400 && lep_pt> 40 && lep_pt<= 110 && met> 0 && met<= 50) {eff2018 = 0.387642; errup = 0.0066101; errdown = 0.00656955;}
    else if (ht> 400 && ht<= 600 && lep_pt> 40 && lep_pt<= 110 && met> 0 && met<= 50) {eff2018 = 0.596085; errup = 0.00792443; errdown = 0.00797407;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 40 && lep_pt<= 110 && met> 0 && met<= 50) {eff2018 = 0.751802; errup = 0.00372287; errdown = 0.00376028;}
    else if (ht> 0 && ht<= 400 && lep_pt> 110 && lep_pt<= 120 && met> 0 && met<= 50) {eff2018 = 0.386105; errup = 0.0171245; errdown = 0.0168562;}
    else if (ht> 400 && ht<= 600 && lep_pt> 110 && lep_pt<= 120 && met> 0 && met<= 50) {eff2018 = 0.643836; errup = 0.0259838; errdown = 0.0268136;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 110 && lep_pt<= 120 && met> 0 && met<= 50) {eff2018 = 0.754565; errup = 0.0118009; errdown = 0.0121843;}
    else if (ht> 0 && ht<= 400 && lep_pt> 120 && lep_pt<= 9999 && met> 0 && met<= 50) {eff2018 = 0.761468; errup = 0.00860319; errdown = 0.00881662;}
    else if (ht> 400 && ht<= 600 && lep_pt> 120 && lep_pt<= 9999 && met> 0 && met<= 50) {eff2018 = 0.800239; errup = 0.0140459; errdown = 0.0147947;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 120 && lep_pt<= 9999 && met> 0 && met<= 50) {eff2018 = 0.87316; errup = 0.00563279; errdown = 0.0058494;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 50 && met<= 75) {eff2018 = 0.00446429; errup = 0.0101903; errdown = 0.00369336;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 50 && met<= 75) {eff2018 = 0.267442; errup = 0.0230545; errdown = 0.0219156;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 50 && met<= 75) {eff2018 = 0.551826; errup = 0.0123281; errdown = 0.0123906;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 50 && met<= 75) {eff2018 = 0.0408163; errup = 0.0195138; errdown = 0.0140021;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 50 && met<= 75) {eff2018 = 0.307937; errup = 0.0282018; errdown = 0.0269105;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 50 && met<= 75) {eff2018 = 0.622134; errup = 0.0139227; errdown = 0.0141213;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 40 && met> 50 && met<= 75) {eff2018 = 0.33518; errup = 0.0266756; errdown = 0.0257133;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 40 && met> 50 && met<= 75) {eff2018 = 0.524887; errup = 0.024796; errdown = 0.0249138;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 40 && met> 50 && met<= 75) {eff2018 = 0.689604; errup = 0.010443; errdown = 0.0106352;}
    else if (ht> 0 && ht<= 400 && lep_pt> 40 && lep_pt<= 110 && met> 50 && met<= 75) {eff2018 = 0.490962; errup = 0.0112954; errdown = 0.0112864;}
    else if (ht> 400 && ht<= 600 && lep_pt> 40 && lep_pt<= 110 && met> 50 && met<= 75) {eff2018 = 0.663012; errup = 0.0115434; errdown = 0.0117363;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 40 && lep_pt<= 110 && met> 50 && met<= 75) {eff2018 = 0.785698; errup = 0.0044274; errdown = 0.00449406;}
    else if (ht> 0 && ht<= 400 && lep_pt> 110 && lep_pt<= 120 && met> 50 && met<= 75) {eff2018 = 0.470199; errup = 0.0304274; errdown = 0.0302191;}
    else if (ht> 400 && ht<= 600 && lep_pt> 110 && lep_pt<= 120 && met> 50 && met<= 75) {eff2018 = 0.636364; errup = 0.0394636; errdown = 0.0412433;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 110 && lep_pt<= 120 && met> 50 && met<= 75) {eff2018 = 0.793506; errup = 0.0148321; errdown = 0.0156277;}
    else if (ht> 0 && ht<= 400 && lep_pt> 120 && lep_pt<= 9999 && met> 50 && met<= 75) {eff2018 = 0.790724; errup = 0.0138995; errdown = 0.014584;}
    else if (ht> 400 && ht<= 600 && lep_pt> 120 && lep_pt<= 9999 && met> 50 && met<= 75) {eff2018 = 0.872; errup = 0.0175109; errdown = 0.0196448;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 120 && lep_pt<= 9999 && met> 50 && met<= 75) {eff2018 = 0.897516; errup = 0.00642761; errdown = 0.00679249;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 75 && met<= 100) {eff2018 = 0.0361446; errup = 0.0339112; errdown = 0.0196074;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 75 && met<= 100) {eff2018 = 0.278607; errup = 0.0352149; errdown = 0.0328462;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 75 && met<= 100) {eff2018 = 0.65252; errup = 0.0144558; errdown = 0.0147338;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 75 && met<= 100) {eff2018 = 0.045977; errup = 0.0348613; errdown = 0.0218756;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 75 && met<= 100) {eff2018 = 0.409639; errup = 0.0416348; errdown = 0.0404647;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 75 && met<= 100) {eff2018 = 0.704675; errup = 0.0157244; errdown = 0.0162081;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 40 && met> 75 && met<= 100) {eff2018 = 0.453988; errup = 0.0422296; errdown = 0.041623;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 40 && met> 75 && met<= 100) {eff2018 = 0.61597; errup = 0.0313617; errdown = 0.0322973;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 40 && met> 75 && met<= 100) {eff2018 = 0.763399; errup = 0.0110119; errdown = 0.0113664;}
    else if (ht> 0 && ht<= 400 && lep_pt> 40 && lep_pt<= 110 && met> 75 && met<= 100) {eff2018 = 0.626178; errup = 0.0160347; errdown = 0.0163076;}
    else if (ht> 400 && ht<= 600 && lep_pt> 40 && lep_pt<= 110 && met> 75 && met<= 100) {eff2018 = 0.76572; errup = 0.0137106; errdown = 0.0142695;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 40 && lep_pt<= 110 && met> 75 && met<= 100) {eff2018 = 0.825713; errup = 0.00481489; errdown = 0.00492038;}
    else if (ht> 0 && ht<= 400 && lep_pt> 110 && lep_pt<= 120 && met> 75 && met<= 100) {eff2018 = 0.610687; errup = 0.0453152; errdown = 0.0471475;}
    else if (ht> 400 && ht<= 600 && lep_pt> 110 && lep_pt<= 120 && met> 75 && met<= 100) {eff2018 = 0.775; errup = 0.0487807; errdown = 0.0565118;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 110 && lep_pt<= 120 && met> 75 && met<= 100) {eff2018 = 0.859083; errup = 0.0145352; errdown = 0.0158254;}
    else if (ht> 0 && ht<= 400 && lep_pt> 120 && lep_pt<= 9999 && met> 75 && met<= 100) {eff2018 = 0.850498; errup = 0.0209401; errdown = 0.0234521;}
    else if (ht> 400 && ht<= 600 && lep_pt> 120 && lep_pt<= 9999 && met> 75 && met<= 100) {eff2018 = 0.859459; errup = 0.0260778; errdown = 0.0303603;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 120 && lep_pt<= 9999 && met> 75 && met<= 100) {eff2018 = 0.914036; errup = 0.00731012; errdown = 0.00789175;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 100 && met<= 125) {eff2018 = 0.111111; errup = 0.0791717; errdown = 0.0524058;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 100 && met<= 125) {eff2018 = 0.3125; errup = 0.0543931; errdown = 0.0500862;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 100 && met<= 125) {eff2018 = 0.710778; errup = 0.0171191; errdown = 0.0177172;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 100 && met<= 125) {eff2018 = 0; errup = 0.0683597; errdown = 0;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 100 && met<= 125) {eff2018 = 0.363636; errup = 0.0583279; errdown = 0.0549108;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 100 && met<= 125) {eff2018 = 0.755848; errup = 0.0167559; errdown = 0.0175366;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 40 && met> 100 && met<= 125) {eff2018 = 0.588235; errup = 0.0577422; errdown = 0.0600294;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 40 && met> 100 && met<= 125) {eff2018 = 0.702479; errup = 0.0436783; errdown = 0.0473382;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 40 && met> 100 && met<= 125) {eff2018 = 0.807843; errup = 0.012509; errdown = 0.0131362;}
    else if (ht> 0 && ht<= 400 && lep_pt> 40 && lep_pt<= 110 && met> 100 && met<= 125) {eff2018 = 0.705495; errup = 0.0219711; errdown = 0.0229195;}
    else if (ht> 400 && ht<= 600 && lep_pt> 40 && lep_pt<= 110 && met> 100 && met<= 125) {eff2018 = 0.808732; errup = 0.018281; errdown = 0.0196374;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 40 && lep_pt<= 110 && met> 100 && met<= 125) {eff2018 = 0.868353; errup = 0.0052035; errdown = 0.00537996;}
    else if (ht> 0 && ht<= 400 && lep_pt> 110 && lep_pt<= 120 && met> 100 && met<= 125) {eff2018 = 0.865385; errup = 0.0482807; errdown = 0.0649611;}
    else if (ht> 400 && ht<= 600 && lep_pt> 110 && lep_pt<= 120 && met> 100 && met<= 125) {eff2018 = 0.942308; errup = 0.0312343; errdown = 0.0529489;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 110 && lep_pt<= 120 && met> 100 && met<= 125) {eff2018 = 0.876147; errup = 0.0159935; errdown = 0.0178418;}
    else if (ht> 0 && ht<= 400 && lep_pt> 120 && lep_pt<= 9999 && met> 100 && met<= 125) {eff2018 = 0.898438; errup = 0.0270367; errdown = 0.0341365;}
    else if (ht> 400 && ht<= 600 && lep_pt> 120 && lep_pt<= 9999 && met> 100 && met<= 125) {eff2018 = 0.9; errup = 0.0338535; errdown = 0.0456142;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 120 && lep_pt<= 9999 && met> 100 && met<= 125) {eff2018 = 0.933071; errup = 0.00788272; errdown = 0.00878826;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 125 && met<= 150) {eff2018 = 0.1875; errup = 0.149399; errdown = 0.100212;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 125 && met<= 150) {eff2018 = 0.583333; errup = 0.0918854; errdown = 0.0971962;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 125 && met<= 150) {eff2018 = 0.799163; errup = 0.0186952; errdown = 0.0200167;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 125 && met<= 150) {eff2018 = 0.238095; errup = 0.128988; errdown = 0.0987131;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 125 && met<= 150) {eff2018 = 0.673077; errup = 0.0701338; errdown = 0.0776765;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 125 && met<= 150) {eff2018 = 0.842784; errup = 0.0188114; errdown = 0.0206985;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 40 && met> 125 && met<= 150) {eff2018 = 0.522727; errup = 0.0849595; errdown = 0.0861304;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 40 && met> 125 && met<= 150) {eff2018 = 0.786885; errup = 0.0548156; errdown = 0.0655726;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 40 && met> 125 && met<= 150) {eff2018 = 0.887324; errup = 0.0119924; errdown = 0.0131472;}
    else if (ht> 0 && ht<= 400 && lep_pt> 40 && lep_pt<= 110 && met> 125 && met<= 150) {eff2018 = 0.865546; errup = 0.0225156; errdown = 0.0258705;}
    else if (ht> 400 && ht<= 600 && lep_pt> 40 && lep_pt<= 110 && met> 125 && met<= 150) {eff2018 = 0.885621; errup = 0.0184512; errdown = 0.0211925;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 40 && lep_pt<= 110 && met> 125 && met<= 150) {eff2018 = 0.908734; errup = 0.00523456; errdown = 0.00551032;}
    else if (ht> 0 && ht<= 400 && lep_pt> 110 && lep_pt<= 120 && met> 125 && met<= 150) {eff2018 = 0.857143; errup = 0.0670785; errdown = 0.0986253;}
    else if (ht> 400 && ht<= 600 && lep_pt> 110 && lep_pt<= 120 && met> 125 && met<= 150) {eff2018 = 0.9; errup = 0.0643201; errdown = 0.116971;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 110 && lep_pt<= 120 && met> 125 && met<= 150) {eff2018 = 0.910448; errup = 0.0176206; errdown = 0.0210035;}
    else if (ht> 0 && ht<= 400 && lep_pt> 120 && lep_pt<= 9999 && met> 125 && met<= 150) {eff2018 = 0.915254; errup = 0.0361154; errdown = 0.0532642;}
    else if (ht> 400 && ht<= 600 && lep_pt> 120 && lep_pt<= 9999 && met> 125 && met<= 150) {eff2018 = 0.940299; errup = 0.0283546; errdown = 0.0446914;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 120 && lep_pt<= 9999 && met> 125 && met<= 150) {eff2018 = 0.953344; errup = 0.00833639; errdown = 0.00987704;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 150 && met<= 175) {eff2018 = 0.5; errup = 0.195182; errdown = 0.195182;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 150 && met<= 175) {eff2018 = 0.666667; errup = 0.106107; errdown = 0.122517;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 150 && met<= 175) {eff2018 = 0.871875; errup = 0.0189792; errdown = 0.021493;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 150 && met<= 175) {eff2018 = 0.666667; errup = 0.17521; errdown = 0.221361;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 150 && met<= 175) {eff2018 = 0.772727; errup = 0.0944237; errdown = 0.12452;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 150 && met<= 175) {eff2018 = 0.9375; errup = 0.014321; errdown = 0.0177218;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 40 && met> 150 && met<= 175) {eff2018 = 0.764706; errup = 0.108959; errdown = 0.147312;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 40 && met> 150 && met<= 175) {eff2018 = 0.913043; errup = 0.0411491; errdown = 0.0634308;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 40 && met> 150 && met<= 175) {eff2018 = 0.924335; errup = 0.0120474; errdown = 0.013929;}
    else if (ht> 0 && ht<= 400 && lep_pt> 40 && lep_pt<= 110 && met> 150 && met<= 175) {eff2018 = 0.894309; errup = 0.028102; errdown = 0.0354144;}
    else if (ht> 400 && ht<= 600 && lep_pt> 40 && lep_pt<= 110 && met> 150 && met<= 175) {eff2018 = 0.94; errup = 0.0168114; errdown = 0.0218501;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 40 && lep_pt<= 110 && met> 150 && met<= 175) {eff2018 = 0.952295; errup = 0.00453493; errdown = 0.00496221;}
    else if (ht> 0 && ht<= 400 && lep_pt> 110 && lep_pt<= 120 && met> 150 && met<= 175) {eff2018 = 1; errup = 0; errdown = 0.142229;}
    else if (ht> 400 && ht<= 600 && lep_pt> 110 && lep_pt<= 120 && met> 150 && met<= 175) {eff2018 = 1; errup = 0; errdown = 0.15411;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 110 && lep_pt<= 120 && met> 150 && met<= 175) {eff2018 = 0.969543; errup = 0.0120157; errdown = 0.0177485;}
    else if (ht> 0 && ht<= 400 && lep_pt> 120 && lep_pt<= 9999 && met> 150 && met<= 175) {eff2018 = 0.969697; errup = 0.0250817; errdown = 0.0662602;}
    else if (ht> 400 && ht<= 600 && lep_pt> 120 && lep_pt<= 9999 && met> 150 && met<= 175) {eff2018 = 0.972973; errup = 0.0223689; errdown = 0.0594217;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 120 && lep_pt<= 9999 && met> 150 && met<= 175) {eff2018 = 0.970297; errup = 0.00752031; errdown = 0.00962462;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 175 && met<= 215) {eff2018 = 0.666667; errup = 0.277375; errdown = 0.414535;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 175 && met<= 215) {eff2018 = 0.8; errup = 0.0835235; errdown = 0.112668;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 175 && met<= 215) {eff2018 = 0.952381; errup = 0.0116187; errdown = 0.0146509;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 175 && met<= 215) {eff2018 = 0.777778; errup = 0.142118; errdown = 0.221429;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 175 && met<= 215) {eff2018 = 0.8; errup = 0.0835235; errdown = 0.112668;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 175 && met<= 215) {eff2018 = 0.983444; errup = 0.00713312; errdown = 0.0110449;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 40 && met> 175 && met<= 215) {eff2018 = 0.9; errup = 0.0643201; errdown = 0.116971;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 40 && met> 175 && met<= 215) {eff2018 = 0.918919; errup = 0.0438002; errdown = 0.0726265;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 40 && met> 175 && met<= 215) {eff2018 = 0.965649; errup = 0.00793934; errdown = 0.00992767;}
    else if (ht> 0 && ht<= 400 && lep_pt> 40 && lep_pt<= 110 && met> 175 && met<= 215) {eff2018 = 0.957983; errup = 0.0180304; errdown = 0.027424;}
    else if (ht> 400 && ht<= 600 && lep_pt> 40 && lep_pt<= 110 && met> 175 && met<= 215) {eff2018 = 0.960199; errup = 0.0136568; errdown = 0.0190433;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 40 && lep_pt<= 110 && met> 175 && met<= 215) {eff2018 = 0.984898; errup = 0.00245839; errdown = 0.00288698;}
    else if (ht> 0 && ht<= 400 && lep_pt> 110 && lep_pt<= 120 && met> 175 && met<= 215) {eff2018 = 1; errup = 0; errdown = 0.168149;}
    else if (ht> 400 && ht<= 600 && lep_pt> 110 && lep_pt<= 120 && met> 175 && met<= 215) {eff2018 = 0.95; errup = 0.0413995; errdown = 0.105764;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 110 && lep_pt<= 120 && met> 175 && met<= 215) {eff2018 = 0.975845; errup = 0.0103945; errdown = 0.0160098;}
    else if (ht> 0 && ht<= 400 && lep_pt> 120 && lep_pt<= 9999 && met> 175 && met<= 215) {eff2018 = 1; errup = 0; errdown = 0.0576587;}
    else if (ht> 400 && ht<= 600 && lep_pt> 120 && lep_pt<= 9999 && met> 175 && met<= 215) {eff2018 = 1; errup = 0; errdown = 0.0472931;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 120 && lep_pt<= 9999 && met> 175 && met<= 215) {eff2018 = 0.994718; errup = 0.00287314; errdown = 0.0051109;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 215 && met<= 9999) {eff2018 = 1; errup = 0; errdown = 0.15411;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 215 && met<= 9999) {eff2018 = 1; errup = 0; errdown = 0.0879414;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 215 && met<= 9999) {eff2018 = 0.991758; errup = 0.00448218; errdown = 0.00795189;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 215 && met<= 9999) {eff2018 = 1; errup = 0; errdown = 0.264229;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 215 && met<= 9999) {eff2018 = 0.947368; errup = 0.0435805; errdown = 0.110836;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 215 && met<= 9999) {eff2018 = 0.988201; errup = 0.00563868; errdown = 0.00923117;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 40 && met> 215 && met<= 9999) {eff2018 = 1; errup = 0; errdown = 0.184992;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 40 && met> 215 && met<= 9999) {eff2018 = 0.979592; errup = 0.0168888; errdown = 0.0453679;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 40 && met> 215 && met<= 9999) {eff2018 = 0.995208; errup = 0.00260705; errdown = 0.00463962;}
    else if (ht> 0 && ht<= 400 && lep_pt> 40 && lep_pt<= 110 && met> 215 && met<= 9999) {eff2018 = 1; errup = 0; errdown = 0.0200277;}
    else if (ht> 400 && ht<= 600 && lep_pt> 40 && lep_pt<= 110 && met> 215 && met<= 9999) {eff2018 = 1; errup = 0; errdown = 0.0104058;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 40 && lep_pt<= 110 && met> 215 && met<= 9999) {eff2018 = 0.996068; errup = 0.00111791; errdown = 0.00148999;}
    else if (ht> 0 && ht<= 400 && lep_pt> 110 && lep_pt<= 120 && met> 215 && met<= 9999) {eff2018 = 1; errup = 0; errdown = 0.308024;}
    else if (ht> 400 && ht<= 600 && lep_pt> 110 && lep_pt<= 120 && met> 215 && met<= 9999) {eff2018 = 1; errup = 0; errdown = 0.0802771;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 110 && lep_pt<= 120 && met> 215 && met<= 9999) {eff2018 = 0.993266; errup = 0.00434837; errdown = 0.00881245;}
    else if (ht> 0 && ht<= 400 && lep_pt> 120 && lep_pt<= 9999 && met> 215 && met<= 9999) {eff2018 = 1; errup = 0; errdown = 0.102638;}
    else if (ht> 400 && ht<= 600 && lep_pt> 120 && lep_pt<= 9999 && met> 215 && met<= 9999) {eff2018 = 1; errup = 0; errdown = 0.0271039;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 120 && lep_pt<= 9999 && met> 215 && met<= 9999) {eff2018 = 0.99867; errup = 0.00110009; errdown = 0.00305118;}
  }
  else if (nmu==1 && nel==0) { // 1 muon CR
    if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 0 && met<= 50) {eff2016 = 0.309278; errup = 0.012177; errdown = 0.0119246;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 0 && met<= 50) {eff2016 = 0.896958; errup = 0.00960512; errdown = 0.0104239;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 0 && met<= 50) {eff2016 = 0.902498; errup = 0.00336926; errdown = 0.00347436;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 0 && met<= 50) {eff2016 = 0.433962; errup = 0.0171325; errdown = 0.016981;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 0 && met<= 50) {eff2016 = 0.867868; errup = 0.0132858; errdown = 0.0144527;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 0 && met<= 50) {eff2016 = 0.908151; errup = 0.00408863; errdown = 0.00425486;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 50 && met> 0 && met<= 50) {eff2016 = 0.631429; errup = 0.0117375; errdown = 0.0118913;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 50 && met> 0 && met<= 50) {eff2016 = 0.929134; errup = 0.00660721; errdown = 0.00719829;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 50 && met> 0 && met<= 50) {eff2016 = 0.935745; errup = 0.0022462; errdown = 0.00232037;}
    else if (ht> 0 && ht<= 400 && lep_pt> 50 && lep_pt<= 9999 && met> 0 && met<= 50) {eff2016 = 0.938377; errup = 0.00596572; errdown = 0.00652796;}
    else if (ht> 400 && ht<= 600 && lep_pt> 50 && lep_pt<= 9999 && met> 0 && met<= 50) {eff2016 = 0.962165; errup = 0.00418374; errdown = 0.00465172;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 50 && lep_pt<= 9999 && met> 0 && met<= 50) {eff2016 = 0.961258; errup = 0.00142423; errdown = 0.00147532;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 50 && met<= 75) {eff2016 = 0.379603; errup = 0.0275657; errdown = 0.0268474;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 50 && met<= 75) {eff2016 = 0.91954; errup = 0.0147031; errdown = 0.0173473;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 50 && met<= 75) {eff2016 = 0.898394; errup = 0.00492584; errdown = 0.00514117;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 50 && met<= 75) {eff2016 = 0.595238; errup = 0.0324355; errdown = 0.0332379;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 50 && met<= 75) {eff2016 = 0.91791; errup = 0.0169135; errdown = 0.0203713;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 50 && met<= 75) {eff2016 = 0.928719; errup = 0.00479442; errdown = 0.00510027;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 50 && met> 50 && met<= 75) {eff2016 = 0.734667; errup = 0.0164532; errdown = 0.0171043;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 50 && met> 50 && met<= 75) {eff2016 = 0.927536; errup = 0.00906659; errdown = 0.0101678;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 50 && met> 50 && met<= 75) {eff2016 = 0.94032; errup = 0.0026961; errdown = 0.00281232;}
    else if (ht> 0 && ht<= 400 && lep_pt> 50 && lep_pt<= 9999 && met> 50 && met<= 75) {eff2016 = 0.948035; errup = 0.00793019; errdown = 0.00915617;}
    else if (ht> 400 && ht<= 600 && lep_pt> 50 && lep_pt<= 9999 && met> 50 && met<= 75) {eff2016 = 0.965458; errup = 0.00537137; errdown = 0.00623961;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 50 && lep_pt<= 9999 && met> 50 && met<= 75) {eff2016 = 0.963843; errup = 0.00173444; errdown = 0.00181629;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 75 && met<= 100) {eff2016 = 0.553333; errup = 0.0433853; errdown = 0.0441515;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 75 && met<= 100) {eff2016 = 0.930233; errup = 0.0194967; errdown = 0.0252262;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 75 && met<= 100) {eff2016 = 0.915444; errup = 0.00580904; errdown = 0.00618079;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 75 && met<= 100) {eff2016 = 0.675214; errup = 0.0457257; errdown = 0.0489983;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 75 && met<= 100) {eff2016 = 0.951515; errup = 0.016605; errdown = 0.0230425;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 75 && met<= 100) {eff2016 = 0.94313; errup = 0.00521512; errdown = 0.00568264;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 50 && met> 75 && met<= 100) {eff2016 = 0.848397; errup = 0.019717; errdown = 0.0218975;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 50 && met> 75 && met<= 100) {eff2016 = 0.973843; errup = 0.00711092; errdown = 0.00928315;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 50 && met> 75 && met<= 100) {eff2016 = 0.950236; errup = 0.00288025; errdown = 0.00304225;}
    else if (ht> 0 && ht<= 400 && lep_pt> 50 && lep_pt<= 9999 && met> 75 && met<= 100) {eff2016 = 0.946067; errup = 0.0107415; errdown = 0.0129539;}
    else if (ht> 400 && ht<= 600 && lep_pt> 50 && lep_pt<= 9999 && met> 75 && met<= 100) {eff2016 = 0.980488; errup = 0.00479992; errdown = 0.006119;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 50 && lep_pt<= 9999 && met> 75 && met<= 100) {eff2016 = 0.967188; errup = 0.00193449; errdown = 0.00204772;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 100 && met<= 125) {eff2016 = 0.671233; errup = 0.0587767; errdown = 0.0640049;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 100 && met<= 125) {eff2016 = 0.914634; errup = 0.0309441; errdown = 0.0429513;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 100 && met<= 125) {eff2016 = 0.933731; errup = 0.00643465; errdown = 0.00703905;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 100 && met<= 125) {eff2016 = 0.775862; errup = 0.0574412; errdown = 0.068324;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 100 && met<= 125) {eff2016 = 0.948454; errup = 0.0220862; errdown = 0.0333658;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 100 && met<= 125) {eff2016 = 0.958838; errup = 0.00565631; errdown = 0.00644963;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 50 && met> 100 && met<= 125) {eff2016 = 0.918478; errup = 0.0203311; errdown = 0.0254746;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 50 && met> 100 && met<= 125) {eff2016 = 0.947368; errup = 0.0132382; errdown = 0.0167863;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 50 && met> 100 && met<= 125) {eff2016 = 0.966378; errup = 0.00292606; errdown = 0.00318206;}
    else if (ht> 0 && ht<= 400 && lep_pt> 50 && lep_pt<= 9999 && met> 100 && met<= 125) {eff2016 = 0.983264; errup = 0.00799298; errdown = 0.0130352;}
    else if (ht> 400 && ht<= 600 && lep_pt> 50 && lep_pt<= 9999 && met> 100 && met<= 125) {eff2016 = 0.985102; errup = 0.00513915; errdown = 0.00726487;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 50 && lep_pt<= 9999 && met> 100 && met<= 125) {eff2016 = 0.97674; errup = 0.00198722; errdown = 0.00215996;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 125 && met<= 150) {eff2016 = 0.828571; errup = 0.0658118; errdown = 0.0883336;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 125 && met<= 150) {eff2016 = 0.955556; errup = 0.0286549; errdown = 0.0556317;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 125 && met<= 150) {eff2016 = 0.970874; errup = 0.00551699; errdown = 0.00663155;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 125 && met<= 150) {eff2016 = 0.973684; errup = 0.02178; errdown = 0.0579268;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 125 && met<= 150) {eff2016 = 0.964286; errup = 0.0230346; errdown = 0.0451714;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 125 && met<= 150) {eff2016 = 0.976055; errup = 0.00514512; errdown = 0.00635002;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 50 && met> 125 && met<= 150) {eff2016 = 0.94382; errup = 0.0240536; errdown = 0.0362177;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 50 && met> 125 && met<= 150) {eff2016 = 0.987952; errup = 0.00777825; errdown = 0.0156694;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 50 && met> 125 && met<= 150) {eff2016 = 0.983293; errup = 0.00258303; errdown = 0.0030079;}
    else if (ht> 0 && ht<= 400 && lep_pt> 50 && lep_pt<= 9999 && met> 125 && met<= 150) {eff2016 = 0.987261; errup = 0.0082239; errdown = 0.0165543;}
    else if (ht> 400 && ht<= 600 && lep_pt> 50 && lep_pt<= 9999 && met> 125 && met<= 150) {eff2016 = 0.988604; errup = 0.00544618; errdown = 0.00891882;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 50 && lep_pt<= 9999 && met> 125 && met<= 150) {eff2016 = 0.982685; errup = 0.00206541; errdown = 0.00232236;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 150 && met<= 175) {eff2016 = 1; errup = 0; errdown = 0.0923495;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 150 && met<= 175) {eff2016 = 1; errup = 0; errdown = 0.0683597;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 150 && met<= 175) {eff2016 = 0.987118; errup = 0.0044459; errdown = 0.00629175;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 150 && met<= 175) {eff2016 = 1; errup = 0; errdown = 0.142229;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 150 && met<= 175) {eff2016 = 1; errup = 0; errdown = 0.0659133;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 150 && met<= 175) {eff2016 = 0.978873; errup = 0.00598041; errdown = 0.00791061;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 50 && met> 150 && met<= 175) {eff2016 = 0.952381; errup = 0.0258048; errdown = 0.044158;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 50 && met> 150 && met<= 175) {eff2016 = 1; errup = 0; errdown = 0.0149771;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 50 && met> 150 && met<= 175) {eff2016 = 0.987744; errup = 0.00258316; errdown = 0.00318481;}
    else if (ht> 0 && ht<= 400 && lep_pt> 50 && lep_pt<= 9999 && met> 150 && met<= 175) {eff2016 = 0.988764; errup = 0.00929678; errdown = 0.0253616;}
    else if (ht> 400 && ht<= 600 && lep_pt> 50 && lep_pt<= 9999 && met> 150 && met<= 175) {eff2016 = 0.986425; errup = 0.0073789; errdown = 0.013028;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 50 && lep_pt<= 9999 && met> 150 && met<= 175) {eff2016 = 0.994496; errup = 0.00135936; errdown = 0.00174222;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 175 && met<= 215) {eff2016 = 1; errup = 0; errdown = 0.123222;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 175 && met<= 215) {eff2016 = 0.973684; errup = 0.02178; errdown = 0.0579268;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 175 && met<= 215) {eff2016 = 0.982704; errup = 0.00511276; errdown = 0.00685994;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 175 && met<= 215) {eff2016 = 1; errup = 0; errdown = 0.0802771;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 175 && met<= 215) {eff2016 = 1; errup = 0; errdown = 0.0512411;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 175 && met<= 215) {eff2016 = 0.998138; errup = 0.00154055; errdown = 0.00426903;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 50 && met> 175 && met<= 215) {eff2016 = 1; errup = 0; errdown = 0.0361508;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 50 && met> 175 && met<= 215) {eff2016 = 0.991379; errup = 0.00713254; errdown = 0.019543;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 50 && met> 175 && met<= 215) {eff2016 = 0.997409; errup = 0.00111857; errdown = 0.00174877;}
    else if (ht> 0 && ht<= 400 && lep_pt> 50 && lep_pt<= 9999 && met> 175 && met<= 215) {eff2016 = 0.985294; errup = 0.0121686; errdown = 0.0330032;}
    else if (ht> 400 && ht<= 600 && lep_pt> 50 && lep_pt<= 9999 && met> 175 && met<= 215) {eff2016 = 1; errup = 0; errdown = 0.00811302;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 50 && lep_pt<= 9999 && met> 175 && met<= 215) {eff2016 = 0.996846; errup = 0.000979628; errdown = 0.00134259;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 215 && met<= 9999) {eff2016 = 1; errup = 0; errdown = 0.205568;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 215 && met<= 9999) {eff2016 = 1; errup = 0; errdown = 0.0738409;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 215 && met<= 9999) {eff2016 = 0.996764; errup = 0.00209004; errdown = 0.00425238;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 215 && met<= 9999) {eff2016 = 1; errup = 0; errdown = 0.23126;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 215 && met<= 9999) {eff2016 = 1; errup = 0; errdown = 0.0559083;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 215 && met<= 9999) {eff2016 = 0.994709; errup = 0.0028782; errdown = 0.00511987;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 50 && met> 215 && met<= 9999) {eff2016 = 1; errup = 0; errdown = 0.0498539;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 50 && met> 215 && met<= 9999) {eff2016 = 1; errup = 0; errdown = 0.0172182;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 50 && met> 215 && met<= 9999) {eff2016 = 0.997996; errup = 0.000958848; errdown = 0.0015817;}
    else if (ht> 0 && ht<= 400 && lep_pt> 50 && lep_pt<= 9999 && met> 215 && met<= 9999) {eff2016 = 1; errup = 0; errdown = 0.0335184;}
    else if (ht> 400 && ht<= 600 && lep_pt> 50 && lep_pt<= 9999 && met> 215 && met<= 9999) {eff2016 = 1; errup = 0; errdown = 0.00783675;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 50 && lep_pt<= 9999 && met> 215 && met<= 9999) {eff2016 = 0.999413; errup = 0.000378822; errdown = 0.00077304;}
    //2017
    if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 0 && met<= 50) {eff2017 = 0.0708592; errup = 0.00848378; errdown = 0.00767985;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 0 && met<= 50) {eff2017 = 0.386643; errup = 0.0149712; errdown = 0.0147661;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 0 && met<= 50) {eff2017 = 0.79144; errup = 0.00809149; errdown = 0.00832374;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 0 && met<= 50) {eff2017 = 0.271298; errup = 0.0170541; errdown = 0.0164307;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 0 && met<= 50) {eff2017 = 0.550075; errup = 0.0198882; errdown = 0.0200436;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 0 && met<= 50) {eff2017 = 0.843672; errup = 0.0091337; errdown = 0.00957456;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 50 && met> 0 && met<= 50) {eff2017 = 0.431312; errup = 0.012669; errdown = 0.0125819;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 50 && met> 0 && met<= 50) {eff2017 = 0.704082; errup = 0.0120989; errdown = 0.0123844;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 50 && met> 0 && met<= 50) {eff2017 = 0.886586; errup = 0.00491944; errdown = 0.00510812;}
    else if (ht> 0 && ht<= 400 && lep_pt> 50 && lep_pt<= 9999 && met> 0 && met<= 50) {eff2017 = 0.887264; errup = 0.00786673; errdown = 0.00835697;}
    else if (ht> 400 && ht<= 600 && lep_pt> 50 && lep_pt<= 9999 && met> 0 && met<= 50) {eff2017 = 0.92204; errup = 0.00652681; errdown = 0.00704377;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 50 && lep_pt<= 9999 && met> 0 && met<= 50) {eff2017 = 0.936299; errup = 0.00313807; errdown = 0.00328508;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 50 && met<= 75) {eff2017 = 0.0614251; errup = 0.0143316; errdown = 0.0119574;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 50 && met<= 75) {eff2017 = 0.414583; errup = 0.0236888; errdown = 0.0233169;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 50 && met<= 75) {eff2017 = 0.814067; errup = 0.00972747; errdown = 0.0101237;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 50 && met<= 75) {eff2017 = 0.325843; errup = 0.0311932; errdown = 0.0298068;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 50 && met<= 75) {eff2017 = 0.614198; errup = 0.0281699; errdown = 0.0289136;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 50 && met<= 75) {eff2017 = 0.857381; errup = 0.0102039; errdown = 0.0108243;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 50 && met> 50 && met<= 75) {eff2017 = 0.543956; errup = 0.0190702; errdown = 0.0191953;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 50 && met> 50 && met<= 75) {eff2017 = 0.780458; errup = 0.0146207; errdown = 0.0153252;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 50 && met> 50 && met<= 75) {eff2017 = 0.893736; errup = 0.0056122; errdown = 0.00587802;}
    else if (ht> 0 && ht<= 400 && lep_pt> 50 && lep_pt<= 9999 && met> 50 && met<= 75) {eff2017 = 0.921429; errup = 0.0093482; errdown = 0.0104152;}
    else if (ht> 400 && ht<= 600 && lep_pt> 50 && lep_pt<= 9999 && met> 50 && met<= 75) {eff2017 = 0.952077; errup = 0.00699194; errdown = 0.00802783;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 50 && lep_pt<= 9999 && met> 50 && met<= 75) {eff2017 = 0.940981; errup = 0.0035208; errdown = 0.00372265;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 75 && met<= 100) {eff2017 = 0.15; errup = 0.0337209; errdown = 0.0288809;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 75 && met<= 100) {eff2017 = 0.45098; errup = 0.0332565; errdown = 0.0328488;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 75 && met<= 100) {eff2017 = 0.842715; errup = 0.0105919; errdown = 0.0111815;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 75 && met<= 100) {eff2017 = 0.520325; errup = 0.0487155; errdown = 0.0490739;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 75 && met<= 100) {eff2017 = 0.647668; errup = 0.0360589; errdown = 0.037699;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 75 && met<= 100) {eff2017 = 0.861354; errup = 0.0115486; errdown = 0.0123748;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 50 && met> 75 && met<= 100) {eff2017 = 0.666667; errup = 0.025072; errdown = 0.0260004;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 50 && met> 75 && met<= 100) {eff2017 = 0.806922; errup = 0.0171572; errdown = 0.0183345;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 50 && met> 75 && met<= 100) {eff2017 = 0.903462; errup = 0.00599439; errdown = 0.00633423;}
    else if (ht> 0 && ht<= 400 && lep_pt> 50 && lep_pt<= 9999 && met> 75 && met<= 100) {eff2017 = 0.941414; errup = 0.0106018; errdown = 0.0125537;}
    else if (ht> 400 && ht<= 600 && lep_pt> 50 && lep_pt<= 9999 && met> 75 && met<= 100) {eff2017 = 0.954397; errup = 0.00843586; errdown = 0.0100577;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 50 && lep_pt<= 9999 && met> 75 && met<= 100) {eff2017 = 0.95625; errup = 0.00337964; errdown = 0.00363757;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 100 && met<= 125) {eff2017 = 0.14; errup = 0.0672437; errdown = 0.0501519;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 100 && met<= 125) {eff2017 = 0.556452; errup = 0.0479412; errdown = 0.0489288;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 100 && met<= 125) {eff2017 = 0.868726; errup = 0.0122566; errdown = 0.0132553;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 100 && met<= 125) {eff2017 = 0.619048; errup = 0.0664109; errdown = 0.0706357;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 100 && met<= 125) {eff2017 = 0.816514; errup = 0.0382778; errdown = 0.044765;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 100 && met<= 125) {eff2017 = 0.919672; errup = 0.0110898; errdown = 0.0125669;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 50 && met> 100 && met<= 125) {eff2017 = 0.73262; errup = 0.0336118; errdown = 0.0362961;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 50 && met> 100 && met<= 125) {eff2017 = 0.855596; errup = 0.0215151; errdown = 0.0242958;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 50 && met> 100 && met<= 125) {eff2017 = 0.931906; errup = 0.00593101; errdown = 0.0064273;}
    else if (ht> 0 && ht<= 400 && lep_pt> 50 && lep_pt<= 9999 && met> 100 && met<= 125) {eff2017 = 0.956364; errup = 0.0122799; errdown = 0.01608;}
    else if (ht> 400 && ht<= 600 && lep_pt> 50 && lep_pt<= 9999 && met> 100 && met<= 125) {eff2017 = 0.972569; errup = 0.00808855; errdown = 0.0108016;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 50 && lep_pt<= 9999 && met> 100 && met<= 125) {eff2017 = 0.968242; errup = 0.00341417; errdown = 0.00378781;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 125 && met<= 150) {eff2017 = 0.423077; errup = 0.116981; errdown = 0.110076;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 125 && met<= 150) {eff2017 = 0.666667; errup = 0.0654592; errdown = 0.0717057;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 125 && met<= 150) {eff2017 = 0.912046; errup = 0.0124958; errdown = 0.0141913;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 125 && met<= 150) {eff2017 = 0.724138; errup = 0.089299; errdown = 0.107524;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 125 && met<= 150) {eff2017 = 0.850746; errup = 0.044684; errdown = 0.0568246;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 125 && met<= 150) {eff2017 = 0.943262; errup = 0.0112896; errdown = 0.0136025;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 50 && met> 125 && met<= 150) {eff2017 = 0.94186; errup = 0.0248848; errdown = 0.0374166;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 50 && met> 125 && met<= 150) {eff2017 = 0.928571; errup = 0.0191788; errdown = 0.024547;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 50 && met> 125 && met<= 150) {eff2017 = 0.956884; errup = 0.00560047; errdown = 0.0063382;}
    else if (ht> 0 && ht<= 400 && lep_pt> 50 && lep_pt<= 9999 && met> 125 && met<= 150) {eff2017 = 0.980392; errup = 0.0106526; errdown = 0.0187052;}
    else if (ht> 400 && ht<= 600 && lep_pt> 50 && lep_pt<= 9999 && met> 125 && met<= 150) {eff2017 = 0.973913; errup = 0.0102998; errdown = 0.0152564;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 50 && lep_pt<= 9999 && met> 125 && met<= 150) {eff2017 = 0.981845; errup = 0.00295218; errdown = 0.0034642;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 150 && met<= 175) {eff2017 = 0.533333; errup = 0.152935; errdown = 0.15827;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 150 && met<= 175) {eff2017 = 0.882353; errup = 0.0554381; errdown = 0.0832913;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 150 && met<= 175) {eff2017 = 0.929412; errup = 0.01398; errdown = 0.0167685;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 150 && met<= 175) {eff2017 = 0.769231; errup = 0.122762; errdown = 0.174724;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 150 && met<= 175) {eff2017 = 0.904762; errup = 0.0450174; errdown = 0.0689193;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 150 && met<= 175) {eff2017 = 0.981308; errup = 0.0073897; errdown = 0.0109973;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 50 && met> 150 && met<= 175) {eff2017 = 1; errup = 0; errdown = 0.0283562;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 50 && met> 150 && met<= 175) {eff2017 = 0.94898; errup = 0.0218627; errdown = 0.0330405;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 50 && met> 150 && met<= 175) {eff2017 = 0.973348; errup = 0.00524904; errdown = 0.00636007;}
    else if (ht> 0 && ht<= 400 && lep_pt> 50 && lep_pt<= 9999 && met> 150 && met<= 175) {eff2017 = 0.988636; errup = 0.00940245; errdown = 0.0256444;}
    else if (ht> 400 && ht<= 600 && lep_pt> 50 && lep_pt<= 9999 && met> 150 && met<= 175) {eff2017 = 0.993711; errup = 0.0052034; errdown = 0.0143129;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 50 && lep_pt<= 9999 && met> 150 && met<= 175) {eff2017 = 0.985469; errup = 0.00306041; errdown = 0.00377044;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 175 && met<= 215) {eff2017 = 1; errup = 0; errdown = 0.142229;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 175 && met<= 215) {eff2017 = 0.888889; errup = 0.0598486; errdown = 0.0963981;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 175 && met<= 215) {eff2017 = 0.98374; errup = 0.00643122; errdown = 0.00958564;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 175 && met<= 215) {eff2017 = 0.947368; errup = 0.0435805; errdown = 0.110836;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 175 && met<= 215) {eff2017 = 1; errup = 0; errdown = 0.0738409;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 175 && met<= 215) {eff2017 = 0.988473; errup = 0.00550887; errdown = 0.00902056;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 50 && met> 175 && met<= 215) {eff2017 = 1; errup = 0; errdown = 0.0392319;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 50 && met> 175 && met<= 215) {eff2017 = 1; errup = 0; errdown = 0.0182418;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 50 && met> 175 && met<= 215) {eff2017 = 0.990955; errup = 0.00295229; errdown = 0.00410363;}
    else if (ht> 0 && ht<= 400 && lep_pt> 50 && lep_pt<= 9999 && met> 175 && met<= 215) {eff2017 = 1; errup = 0; errdown = 0.0249041;}
    else if (ht> 400 && ht<= 600 && lep_pt> 50 && lep_pt<= 9999 && met> 175 && met<= 215) {eff2017 = 0.993711; errup = 0.0052034; errdown = 0.0143129;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 50 && lep_pt<= 9999 && met> 175 && met<= 215) {eff2017 = 0.995131; errup = 0.00168324; errdown = 0.00239245;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 215 && met<= 9999) {eff2017 = 1; errup = 0; errdown = 0.205568;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 215 && met<= 9999) {eff2017 = 1; errup = 0; errdown = 0.0709947;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 215 && met<= 9999) {eff2017 = 0.992788; errup = 0.00392226; errdown = 0.00696502;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 215 && met<= 9999) {eff2017 = 1; errup = 0; errdown = 0.15411;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 215 && met<= 9999) {eff2017 = 1; errup = 0; errdown = 0.0972223;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 215 && met<= 9999) {eff2017 = 0.99505; errup = 0.00319693; errdown = 0.00649192;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 50 && met> 215 && met<= 9999) {eff2017 = 1; errup = 0; errdown = 0.0542609;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 50 && met> 215 && met<= 9999) {eff2017 = 1; errup = 0; errdown = 0.0224723;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 50 && met> 215 && met<= 9999) {eff2017 = 0.994958; errup = 0.00199818; errdown = 0.00299934;}
    else if (ht> 0 && ht<= 400 && lep_pt> 50 && lep_pt<= 9999 && met> 215 && met<= 9999) {eff2017 = 1; errup = 0; errdown = 0.0461088;}
    else if (ht> 400 && ht<= 600 && lep_pt> 50 && lep_pt<= 9999 && met> 215 && met<= 9999) {eff2017 = 1; errup = 0; errdown = 0.012972;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 50 && lep_pt<= 9999 && met> 215 && met<= 9999) {eff2017 = 0.997174; errup = 0.00112047; errdown = 0.0016842;}
    //2018
    if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 0 && met<= 50) {eff2018 = 0.0509554; errup = 0.0103616; errdown = 0.00880306;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 0 && met<= 50) {eff2018 = 0.406528; errup = 0.0160299; errdown = 0.0158392;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 0 && met<= 50) {eff2018 = 0.920266; errup = 0.00495963; errdown = 0.00524822;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 0 && met<= 50) {eff2018 = 0.296203; errup = 0.024754; errdown = 0.0236669;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 0 && met<= 50) {eff2018 = 0.57037; errup = 0.019671; errdown = 0.0198874;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 0 && met<= 50) {eff2018 = 0.918153; errup = 0.00631369; errdown = 0.00677068;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 50 && met> 0 && met<= 50) {eff2018 = 0.45394; errup = 0.0171845; errdown = 0.017079;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 50 && met> 0 && met<= 50) {eff2018 = 0.734914; errup = 0.0120119; errdown = 0.0123595;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 50 && met> 0 && met<= 50) {eff2018 = 0.94219; errup = 0.00339353; errdown = 0.0035852;}
    else if (ht> 0 && ht<= 400 && lep_pt> 50 && lep_pt<= 9999 && met> 0 && met<= 50) {eff2018 = 0.904807; errup = 0.00908061; errdown = 0.00988317;}
    else if (ht> 400 && ht<= 600 && lep_pt> 50 && lep_pt<= 9999 && met> 0 && met<= 50) {eff2018 = 0.939845; errup = 0.00607334; errdown = 0.00667259;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 50 && lep_pt<= 9999 && met> 0 && met<= 50) {eff2018 = 0.953995; errup = 0.00248988; errdown = 0.00262123;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 50 && met<= 75) {eff2018 = 0.0708333; errup = 0.0206706; errdown = 0.0166547;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 50 && met<= 75) {eff2018 = 0.462222; errup = 0.0246733; errdown = 0.0244977;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 50 && met<= 75) {eff2018 = 0.905621; errup = 0.00668024; errdown = 0.00711475;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 50 && met<= 75) {eff2018 = 0.358382; errup = 0.0401117; errdown = 0.0383512;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 50 && met<= 75) {eff2018 = 0.66548; errup = 0.0292488; errdown = 0.0304982;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 50 && met<= 75) {eff2018 = 0.924117; errup = 0.00730077; errdown = 0.00797078;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 50 && met> 50 && met<= 75) {eff2018 = 0.5629; errup = 0.0238039; errdown = 0.0240843;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 50 && met> 50 && met<= 75) {eff2018 = 0.770013; errup = 0.015271; errdown = 0.0159857;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 50 && met> 50 && met<= 75) {eff2018 = 0.950936; errup = 0.0036475; errdown = 0.00391321;}
    else if (ht> 0 && ht<= 400 && lep_pt> 50 && lep_pt<= 9999 && met> 50 && met<= 75) {eff2018 = 0.923333; errup = 0.010941; errdown = 0.0124589;}
    else if (ht> 400 && ht<= 600 && lep_pt> 50 && lep_pt<= 9999 && met> 50 && met<= 75) {eff2018 = 0.946875; errup = 0.00726638; errdown = 0.00826383;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 50 && lep_pt<= 9999 && met> 50 && met<= 75) {eff2018 = 0.95772; errup = 0.00278904; errdown = 0.00297018;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 75 && met<= 100) {eff2018 = 0.105263; errup = 0.0415156; errdown = 0.0318852;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 75 && met<= 100) {eff2018 = 0.44898; errup = 0.0339661; errdown = 0.033524;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 75 && met<= 100) {eff2018 = 0.924188; errup = 0.0071534; errdown = 0.0077968;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 75 && met<= 100) {eff2018 = 0.410256; errup = 0.0629437; errdown = 0.0603993;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 75 && met<= 100) {eff2018 = 0.71978; errup = 0.0346421; errdown = 0.0372475;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 75 && met<= 100) {eff2018 = 0.93621; errup = 0.00752273; errdown = 0.00839199;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 50 && met> 75 && met<= 100) {eff2018 = 0.681159; errup = 0.0291135; errdown = 0.0305079;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 50 && met> 75 && met<= 100) {eff2018 = 0.833676; errup = 0.0171643; errdown = 0.0186162;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 50 && met> 75 && met<= 100) {eff2018 = 0.953016; errup = 0.00387315; errdown = 0.00418804;}
    else if (ht> 0 && ht<= 400 && lep_pt> 50 && lep_pt<= 9999 && met> 75 && met<= 100) {eff2018 = 0.918539; errup = 0.0146216; errdown = 0.0171965;}
    else if (ht> 400 && ht<= 600 && lep_pt> 50 && lep_pt<= 9999 && met> 75 && met<= 100) {eff2018 = 0.956954; errup = 0.00826908; errdown = 0.00993237;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 50 && lep_pt<= 9999 && met> 75 && met<= 100) {eff2018 = 0.965588; errup = 0.00275628; errdown = 0.00297726;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 100 && met<= 125) {eff2018 = 0.111111; errup = 0.0791717; errdown = 0.0524058;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 100 && met<= 125) {eff2018 = 0.519685; errup = 0.0478966; errdown = 0.0482324;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 100 && met<= 125) {eff2018 = 0.938907; errup = 0.00787908; errdown = 0.00888368;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 100 && met<= 125) {eff2018 = 0.68; errup = 0.102272; errdown = 0.119276;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 100 && met<= 125) {eff2018 = 0.788462; errup = 0.0415895; errdown = 0.0477639;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 100 && met<= 125) {eff2018 = 0.947566; errup = 0.00790478; errdown = 0.00911009;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 50 && met> 100 && met<= 125) {eff2018 = 0.703125; errup = 0.042389; errdown = 0.0458523;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 50 && met> 100 && met<= 125) {eff2018 = 0.882759; errup = 0.0191698; errdown = 0.0220444;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 50 && met> 100 && met<= 125) {eff2018 = 0.965564; errup = 0.00386237; errdown = 0.00430287;}
    else if (ht> 0 && ht<= 400 && lep_pt> 50 && lep_pt<= 9999 && met> 100 && met<= 125) {eff2018 = 0.976303; errup = 0.0101982; errdown = 0.0157124;}
    else if (ht> 400 && ht<= 600 && lep_pt> 50 && lep_pt<= 9999 && met> 100 && met<= 125) {eff2018 = 0.971354; errup = 0.00844407; errdown = 0.01127;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 50 && lep_pt<= 9999 && met> 100 && met<= 125) {eff2018 = 0.974014; errup = 0.00275157; errdown = 0.00305009;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 125 && met<= 150) {eff2018 = 0.1; errup = 0.116971; errdown = 0.0643201;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 125 && met<= 150) {eff2018 = 0.627119; errup = 0.0683561; errdown = 0.0731904;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 125 && met<= 150) {eff2018 = 0.973875; errup = 0.00605301; errdown = 0.00759168;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 125 && met<= 150) {eff2018 = 0.75; errup = 0.115499; errdown = 0.153966;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 125 && met<= 150) {eff2018 = 0.869565; errup = 0.0504938; errdown = 0.0697763;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 125 && met<= 150) {eff2018 = 0.97757; errup = 0.00634716; errdown = 0.00839086;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 50 && met> 125 && met<= 150) {eff2018 = 0.923077; errup = 0.032824; errdown = 0.0486878;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 50 && met> 125 && met<= 150) {eff2018 = 0.923858; errup = 0.01902; errdown = 0.0238853;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 50 && met> 125 && met<= 150) {eff2018 = 0.980064; errup = 0.00353709; errdown = 0.0042104;}
    else if (ht> 0 && ht<= 400 && lep_pt> 50 && lep_pt<= 9999 && met> 125 && met<= 150) {eff2018 = 0.963964; errup = 0.0171678; errdown = 0.0275761;}
    else if (ht> 400 && ht<= 600 && lep_pt> 50 && lep_pt<= 9999 && met> 125 && met<= 150) {eff2018 = 0.995305; errup = 0.00388411; errdown = 0.0107125;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 50 && lep_pt<= 9999 && met> 125 && met<= 150) {eff2018 = 0.984236; errup = 0.00249933; errdown = 0.00292235;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 150 && met<= 175) {eff2018 = 0.875; errup = 0.103637; errdown = 0.23225;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 150 && met<= 175) {eff2018 = 0.723404; errup = 0.0696186; errdown = 0.0805167;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 150 && met<= 175) {eff2018 = 0.987923; errup = 0.00520699; errdown = 0.00808753;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 150 && met<= 175) {eff2018 = 0.941176; errup = 0.048713; errdown = 0.12258;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 150 && met<= 175) {eff2018 = 0.962963; errup = 0.0306592; errdown = 0.080075;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 150 && met<= 175) {eff2018 = 0.980978; errup = 0.00698888; errdown = 0.0100953;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 50 && met> 150 && met<= 175) {eff2018 = 0.971429; errup = 0.0236478; errdown = 0.0626552;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 50 && met> 150 && met<= 175) {eff2018 = 0.957447; errup = 0.0202555; errdown = 0.0323679;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 50 && met> 150 && met<= 175) {eff2018 = 0.981584; errup = 0.00406048; errdown = 0.00504708;}
    else if (ht> 0 && ht<= 400 && lep_pt> 50 && lep_pt<= 9999 && met> 150 && met<= 175) {eff2018 = 0.967742; errup = 0.0208084; errdown = 0.0409676;}
    else if (ht> 400 && ht<= 600 && lep_pt> 50 && lep_pt<= 9999 && met> 150 && met<= 175) {eff2018 = 0.986667; errup = 0.00860748; errdown = 0.0173148;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 50 && lep_pt<= 9999 && met> 150 && met<= 175) {eff2018 = 0.992337; errup = 0.00201909; errdown = 0.00263025;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 175 && met<= 215) {eff2018 = 0.9; errup = 0.082873; errdown = 0.194135;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 175 && met<= 215) {eff2018 = 0.9; errup = 0.0539222; errdown = 0.0877974;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 175 && met<= 215) {eff2018 = 0.989035; errup = 0.00472821; errdown = 0.00734954;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 175 && met<= 215) {eff2018 = 0.857143; errup = 0.11848; errdown = 0.257124;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 175 && met<= 215) {eff2018 = 1; errup = 0; errdown = 0.0659133;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 175 && met<= 215) {eff2018 = 0.994819; errup = 0.00334597; errdown = 0.00679283;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 50 && met> 175 && met<= 215) {eff2018 = 0.971429; errup = 0.0236478; errdown = 0.0626552;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 50 && met> 175 && met<= 215) {eff2018 = 0.98913; errup = 0.00899357; errdown = 0.0245495;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 50 && met> 175 && met<= 215) {eff2018 = 0.989756; errup = 0.00279678; errdown = 0.00367595;}
    else if (ht> 0 && ht<= 400 && lep_pt> 50 && lep_pt<= 9999 && met> 175 && met<= 215) {eff2018 = 0.967213; errup = 0.0211491; errdown = 0.0416131;}
    else if (ht> 400 && ht<= 600 && lep_pt> 50 && lep_pt<= 9999 && met> 175 && met<= 215) {eff2018 = 1; errup = 0; errdown = 0.0133482;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 50 && lep_pt<= 9999 && met> 175 && met<= 215) {eff2018 = 0.99708; errup = 0.00115753; errdown = 0.0017398;}
    else if (ht> 0 && ht<= 400 && lep_pt> 20 && lep_pt<= 25 && met> 215 && met<= 9999) {eff2018 = 1; errup = 0; errdown = 0.601684;}
    else if (ht> 400 && ht<= 600 && lep_pt> 20 && lep_pt<= 25 && met> 215 && met<= 9999) {eff2018 = 1; errup = 0; errdown = 0.0769247;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 20 && lep_pt<= 25 && met> 215 && met<= 9999) {eff2018 = 1; errup = 0; errdown = 0.00350724;}
    else if (ht> 0 && ht<= 400 && lep_pt> 25 && lep_pt<= 30 && met> 215 && met<= 9999) {eff2018 = 1; errup = 0; errdown = 0.601684;}
    else if (ht> 400 && ht<= 600 && lep_pt> 25 && lep_pt<= 30 && met> 215 && met<= 9999) {eff2018 = 1; errup = 0; errdown = 0.0802771;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 25 && lep_pt<= 30 && met> 215 && met<= 9999) {eff2018 = 0.995204; errup = 0.00309728; errdown = 0.00629067;}
    else if (ht> 0 && ht<= 400 && lep_pt> 30 && lep_pt<= 50 && met> 215 && met<= 9999) {eff2018 = 1; errup = 0; errdown = 0.0636358;}
    else if (ht> 400 && ht<= 600 && lep_pt> 30 && lep_pt<= 50 && met> 215 && met<= 9999) {eff2018 = 1; errup = 0; errdown = 0.016449;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 30 && lep_pt<= 50 && met> 215 && met<= 9999) {eff2018 = 0.998025; errup = 0.00107467; errdown = 0.00191738;}
    else if (ht> 0 && ht<= 400 && lep_pt> 50 && lep_pt<= 9999 && met> 215 && met<= 9999) {eff2018 = 1; errup = 0; errdown = 0.0439098;}
    else if (ht> 400 && ht<= 600 && lep_pt> 50 && lep_pt<= 9999 && met> 215 && met<= 9999) {eff2018 = 1; errup = 0; errdown = 0.01428;}
    else if (ht> 600 && ht<= 9999 && lep_pt> 50 && lep_pt<= 9999 && met> 215 && met<= 9999) {eff2018 = 0.998493; errup = 0.000720896; errdown = 0.00118964;}
  }
  else if (nel==2 && nmu==0) { // 2 electron CR
    if (lep_pt> 40 && lep_pt<= 45) {eff2016 = 0.938144; errup = 0.0141214; errdown = 0.0141214;}
    else if (lep_pt> 45 && lep_pt<= 50) {eff2016 = 0.949468; errup = 0.0112961; errdown = 0.0112961;}
    else if (lep_pt> 50 && lep_pt<= 55) {eff2016 = 0.966443; errup = 0.00851778; errdown = 0.00851778;}
    else if (lep_pt> 55 && lep_pt<= 60) {eff2016 = 0.962617; errup = 0.0082014; errdown = 0.0082014;}
    else if (lep_pt> 60 && lep_pt<= 65) {eff2016 = 0.972632; errup = 0.00748604; errdown = 0.00748604;}
    else if (lep_pt> 65 && lep_pt<= 70) {eff2016 = 0.973626; errup = 0.00751234; errdown = 0.00751234;}
    else if (lep_pt> 70 && lep_pt<= 75) {eff2016 = 0.960688; errup = 0.00963289; errdown = 0.00963289;}
    else if (lep_pt> 75 && lep_pt<= 80) {eff2016 = 0.983213; errup = 0.00629125; errdown = 0.00629125;}
    else if (lep_pt> 80 && lep_pt<= 85) {eff2016 = 0.987531; errup = 0.00554136; errdown = 0.00554136;}
    else if (lep_pt> 85 && lep_pt<= 90) {eff2016 = 0.973988; errup = 0.00855701; errdown = 0.00855701;}
    else if (lep_pt> 90 && lep_pt<= 95) {eff2016 = 0.985876; errup = 0.00627181; errdown = 0.00627181;}
    else if (lep_pt> 95 && lep_pt<= 100) {eff2016 = 0.974522; errup = 0.00889224; errdown = 0.00889224;}
    else if (lep_pt> 100 && lep_pt<= 105) {eff2016 = 0.977199; errup = 0.00851926; errdown = 0.00851926;}
    else if (lep_pt> 105 && lep_pt<= 110) {eff2016 = 0.985294; errup = 0.00729868; errdown = 0.00729868;}
    else if (lep_pt> 110 && lep_pt<= 9999) {eff2016 = 0.995914; errup = 0.000770677; errdown = 0.000770677;}
    //2017
    if (lep_pt> 40 && lep_pt<= 45) {eff2017 = 0.897959; errup = 0.0249664; errdown = 0.0249664;}
    else if (lep_pt> 45 && lep_pt<= 50) {eff2017 = 0.949309; errup = 0.0148916; errdown = 0.0148916;}
    else if (lep_pt> 50 && lep_pt<= 55) {eff2017 = 0.939502; errup = 0.0142222; errdown = 0.0142222;}
    else if (lep_pt> 55 && lep_pt<= 60) {eff2017 = 0.930233; errup = 0.0158603; errdown = 0.0158603;}
    else if (lep_pt> 60 && lep_pt<= 65) {eff2017 = 0.919014; errup = 0.0161885; errdown = 0.0161885;}
    else if (lep_pt> 65 && lep_pt<= 70) {eff2017 = 0.931818; errup = 0.0155131; errdown = 0.0155131;}
    else if (lep_pt> 70 && lep_pt<= 75) {eff2017 = 0.94382; errup = 0.0140922; errdown = 0.0140922;}
    else if (lep_pt> 75 && lep_pt<= 80) {eff2017 = 0.973913; errup = 0.0105101; errdown = 0.0105101;}
    else if (lep_pt> 80 && lep_pt<= 85) {eff2017 = 0.941704; errup = 0.01569; errdown = 0.01569;}
    else if (lep_pt> 85 && lep_pt<= 90) {eff2017 = 0.960591; errup = 0.0136558; errdown = 0.0136558;}
    else if (lep_pt> 90 && lep_pt<= 95) {eff2017 = 0.919075; errup = 0.0207345; errdown = 0.0207345;}
    else if (lep_pt> 95 && lep_pt<= 100) {eff2017 = 0.941176; errup = 0.0180462; errdown = 0.0180462;}
    else if (lep_pt> 100 && lep_pt<= 105) {eff2017 = 0.972222; errup = 0.0122488; errdown = 0.0122488;}
    else if (lep_pt> 105 && lep_pt<= 110) {eff2017 = 0.934211; errup = 0.0201085; errdown = 0.0201085;}
    else if (lep_pt> 110 && lep_pt<= 9999) {eff2017 = 0.983325; errup = 0.00197631; errdown = 0.00197631;}
    //2018
    if (lep_pt> 40 && lep_pt<= 45) {eff2018 = 0.897959; errup = 0.0249664; errdown = 0.0249664;}
    else if (lep_pt> 45 && lep_pt<= 50) {eff2018 = 0.950413; errup = 0.0139551; errdown = 0.0139551;}
    else if (lep_pt> 50 && lep_pt<= 55) {eff2018 = 0.968085; errup = 0.0104672; errdown = 0.0104672;}
    else if (lep_pt> 55 && lep_pt<= 60) {eff2018 = 0.95203; errup = 0.0129816; errdown = 0.0129816;}
    else if (lep_pt> 60 && lep_pt<= 65) {eff2018 = 0.92926; errup = 0.0145385; errdown = 0.0145385;}
    else if (lep_pt> 65 && lep_pt<= 70) {eff2018 = 0.953846; errup = 0.0130124; errdown = 0.0130124;}
    else if (lep_pt> 70 && lep_pt<= 75) {eff2018 = 0.962963; errup = 0.0114932; errdown = 0.0114932;}
    else if (lep_pt> 75 && lep_pt<= 80) {eff2018 = 0.960937; errup = 0.012109; errdown = 0.012109;}
    else if (lep_pt> 80 && lep_pt<= 85) {eff2018 = 0.952381; errup = 0.0140117; errdown = 0.0140117;}
    else if (lep_pt> 85 && lep_pt<= 90) {eff2018 = 0.961353; errup = 0.0133973; errdown = 0.0133973;}
    else if (lep_pt> 90 && lep_pt<= 95) {eff2018 = 0.975248; errup = 0.0109318; errdown = 0.0109318;}
    else if (lep_pt> 95 && lep_pt<= 100) {eff2018 = 0.978022; errup = 0.0108676; errdown = 0.0108676;}
    else if (lep_pt> 100 && lep_pt<= 105) {eff2018 = 0.979899; errup = 0.00994873; errdown = 0.00994873;}
    else if (lep_pt> 105 && lep_pt<= 110) {eff2018 = 0.987578; errup = 0.00872921; errdown = 0.00872921;}
    else if (lep_pt> 110 && lep_pt<= 9999) {eff2018 = 0.996761; errup = 0.000784379; errdown = 0.000784379;}
  }
  else if (nmu==2 && nel==0) { // 2 muon CR
    if (lep_pt> 40 && lep_pt<= 45) {eff2016 = 0.978495; errup = 0.0075211; errdown = 0.0075211;}
    else if (lep_pt> 45 && lep_pt<= 50) {eff2016 = 0.993151; errup = 0.0030526; errdown = 0.0030526;}
    else if (lep_pt> 50 && lep_pt<= 9999) {eff2016 = 0.99209; errup = 0.000276502; errdown = 0.000276502;}
    //2017
    if (lep_pt> 40 && lep_pt<= 45) {eff2017 = 0.948413; errup = 0.0139338; errdown = 0.0139338;}
    else if (lep_pt> 45 && lep_pt<= 50) {eff2017 = 0.963134; errup = 0.0090451; errdown = 0.0090451;}
    else if (lep_pt> 50 && lep_pt<= 9999) {eff2017 = 0.983369; errup = 0.000433117; errdown = 0.000433117;}
    //2018
    if (lep_pt> 40 && lep_pt<= 45) {eff2018 = 0.984496; errup = 0.00769161; errdown = 0.00769161;}
    else if (lep_pt> 45 && lep_pt<= 50) {eff2018 = 0.960106; errup = 0.0100929; errdown = 0.0100929;}
    else if (lep_pt> 50 && lep_pt<= 9999) {eff2018 = 0.988916; errup = 0.000322656; errdown = 0.000322656;}
  }
  return (36.3*eff2016+41.5*eff2017+59.7*eff2018)/137.5;
}


}
