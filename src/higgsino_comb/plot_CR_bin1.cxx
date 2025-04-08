#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <bitset>
#include <unistd.h>
#include <regex>
#include <getopt.h>
#include "TError.h"
#include "TColor.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/utilities.hpp"
#include "bbgg/lester_mt2_bisect.h"
#include <iomanip>
#include <sstream>
using namespace std;
using namespace PlotOptTypes;
using namespace ROOT::Math;

long double DeltaPhi(long double phi1, long double phi2){
  long double dphi = fmod(fabs(phi2-phi1), 2.L*3.14159);
  return dphi>3.14159 ? 2.L*3.14159-dphi : dphi;
}

long double AddInQuadrature(long double x, long double y){
  if(fabs(y)>fabs(x)){
    const long double temp = y;
    y=x;
    x=temp;
  }
  if(x==0.) return y;
  const long double rat=y/x;
  return fabs(x)*sqrt(1.0L+rat*rat);
}

float dR(float eta1, float eta2, float phi1, float phi2) {
  return AddInQuadrature(eta1-eta2, DeltaPhi(phi1,phi2));
}

double roundDouble( double x, int n )
{
   stringstream ss;
   ss << scientific << setprecision( n - 1 ) << x; 
   return stod( ss.str() );
}


std::vector<std::string> BkgSampleNames = {"DiPhotonJetsBox_", 
                                          "DiPhotonJetsBox1BJet", 
                                          "DiPhotonJetsBox2BJets",
                                          "GJet_Pt-20to40",
                                          "GJet_Pt-40toInf",
                                          "TTGG_0Jets",
                                          "TTGJets",
                                          "TTTo2L2Nu",
                                          "QCD_Pt-30to40",
                                          "QCD_Pt-40"};

std::vector<std::string> SigSampleNames = {"SMS-TChiHH_mChi-150_mLSP-1_HToGG", 
                                          "SMS-TChiHH_mChi-300_mLSP-1_HToGG", 
                                          "SMS-TChiHH_mChi-500_mLSP-1_HToGG"};

//Sept production:
//string bfolderbkg("/net/cms11/cms11r0/pico/NanoAODv7/bbgg_higgsino_modifiedprod_27sep23_2/2016/mc/unskimmed/");
//string bfolderbkg("/net/cms11/cms11r0/pico/NanoAODv7/bbgg_higgsino_signal_prod_10/2016/mc/unskimmed/");
string bfoldersig("/net/cms17/cms17r0/pico/NanoAODv7/bbgg_higgsino_signal_prod_11/2016/mc/unskimmed/");
//string bfolderdata("/net/cms17/cms17r0/pico/NanoAODv7/bbgg_higgsino_signal_prod_10/2016/data/unskimmed/");

string debugsigfile("/net/cms17/cms17r0/pico/NanoAODv7/bbgg_higgsino_signal_prod_11/2016/TChiHH/unskimmed/pico_SMS-TChiHH_mChi-150_mLSP-1_HToGG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISummer16NanoAODv7__PUSummer16v3Fast_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1.root");

string diphlessmcpath("/net/cms37/data1/mhussain/HH-MET/DataSample/Background/picos/");
//Oct production: bb candidate medium cut placed, changed 2D cross sec for the signal, using photon_pt thresh in draw_pico 
//string bfolderbkg("/net/cms11/cms11r0/pico/NanoAODv7/bbgg_higgsino_bkg_prod_13Oct23/2016/mc/unskimmed/");
//string bfolderbkg("/net/cms17/cms17r0/pico/NanoAODv7/bbgg_higgsino_bkg_prod_13Oct23/2016/mc/unskimmed/");
//string bfoldersig("/net/cms17/cms17r0/pico/NanoAODv7/bbgg_higgsino_bkg_prod_13Oct23/2016/TChiHH/unskimmed/");

string Debug_file_diphoton("/net/cms11/cms11r0/pico/NanoAODv9/bbgg_higgsino_signal_prod_11/2016/mc/unskimmed/pico_DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa__RunIISummer20UL16NanoAODv9__106X_mcRun2_asymptotic_v17-v2__270000__E5C3E208-DDB5-8E4A-B3A7-61F981802CF4.root");
//string Debug_file_diphoton("/net/cms37/data1/mhussain/HH-MET/nano2pico/output/unit_test6/mc/unskimmed/pico_DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa__RunIISummer20UL16NanoAODv9__106X_mcRun2_asymptotic_v17-v2__270000__E5C3E208-DDB5-8E4A-B3A7-61F981802CF4.root");
//March 2024 Run 2 UL production
string bfolderbkg("/net/cms11/cms11r0/pico/NanoAODv9/bbgg_higgsino_signal_prod_11/20*/mc/unskimmed/");
string bfolderdata_16("/net/cms11/cms11r0/pico/NanoAODv9/bbgg_higgsino_signal_prod_11/2016*/data/unskimmed/");
string bfolderdata_17("/net/cms11/cms11r0/pico/NanoAODv9/bbgg_higgsino_signal_prod_11/2017/data/unskimmed/");
string bfolderdata_18("/net/cms11/cms11r0/pico/NanoAODv9/bbgg_higgsino_signal_prod_11/2018/data/unskimmed/");
//string bfolderdata("/net/cms11/cms11r0/pico/NanoAODv9/bbgg_higgsino_signal_prod_11/20*/data/unskimmed/");
string bfolderbkg_gjet2("/net/cms11/cms11r0/pico/NanoAODv2/bbgg_higgsino_signal_prod_11/201*/mc/unskimmed/");

string mchi_150(bfoldersig + "*" + SigSampleNames[0] + "*");
string mchi_500(bfoldersig + "*" + SigSampleNames[1] + "*");
string mchi_1000(bfoldersig + "*" + SigSampleNames[2] + "*");

//Alternative Diphoton (madgraph sample)
string DiPhoton_mad("/net/cms17/cms17r0/pico/NanoAODv9/bbgg_higgsino_signal_prod_10/2016/mc/unskimmed/*");
string DiPhoton(bfolderbkg + "*" + BkgSampleNames[0] + "*");
string DiPhoton1B(bfolderbkg + "*" + BkgSampleNames[1] + "*");
string DiPhoton2B(bfolderbkg + "*" + BkgSampleNames[2] + "*");
string GJet1(bfolderbkg + "*" + BkgSampleNames[3] + "*");
string GJet2(bfolderbkg_gjet2 + "*" + BkgSampleNames[4] + "*");
string TTGG(bfolderbkg + "*" + BkgSampleNames[5] + "*");
string TTGJets(bfolderbkg + "*" + BkgSampleNames[6] + "*");
string TTTo2L2Nu(bfolderbkg + "*" + BkgSampleNames[7] + "*");
string QCD_Pt1(bfolderbkg + "*" + BkgSampleNames[8] + "*");
string QCD_Pt2(bfolderbkg + "*" + BkgSampleNames[9] + "*");

string data_16(bfolderdata_16 + "*DoubleEG*");
string data_17(bfolderdata_17 + "*DoubleEG*");
string data_18(bfolderdata_18 + "*EGamma*");

//float Luminosity = 36.3;
float Luminosity = 1.0; //If running the full UL Run 2


//////////////////// Overlap removal code

const NamedFunc num_gengamma_basic("num_gengamma_basic",[](const Baby &b) -> NamedFunc::ScalarType{ 
  
  int number_of_promptgamma = 0;

  for (size_t mc_idx = 0; mc_idx<b.mc_pt()->size(); mc_idx++) {

    bitset<15> mc_statusFlags(b.mc_statusflag()->at(mc_idx));
    
    if(b.mc_id()->at(mc_idx) == 22 && mc_statusFlags[0] && b.mc_pt()->at(mc_idx) >= 10){
      
      if (b.mc_mom()->at(mc_idx) == 5 || b.mc_mom()->at(mc_idx) == -5 || b.mc_mom()->at(mc_idx) == 4 || b.mc_mom()->at(mc_idx) == -4 || b.mc_mom()->at(mc_idx) == 3 || b.mc_mom()->at(mc_idx) == -3 || b.mc_mom()->at(mc_idx) == 2 || b.mc_mom()->at(mc_idx) == -2 || b.mc_mom()->at(mc_idx) == 1 || b.mc_mom()->at(mc_idx) == -1 || b.mc_mom()->at(mc_idx) == 21){
        number_of_promptgamma += 1;
      }

    }

  }

  return number_of_promptgamma;

});


const NamedFunc num_gen_b("num_gen_b",[](const Baby &b) -> NamedFunc::ScalarType{
  
  int number_of_b = 0;
  
  for (size_t mc_idx = 0; mc_idx<b.mc_pt()->size(); mc_idx++) {

    if(b.mc_id()->at(mc_idx) == 5 || b.mc_id()->at(mc_idx) == -5){
      
      bitset<15> mc_statusFlags(b.mc_statusflag()->at(mc_idx));

      if(mc_statusFlags[7] && b.mc_pt()->at(mc_idx) >= 10){  // isHardProcess and pt >= 10

        TVector3 compPart,genb;
        genb.SetPtEtaPhi(b.mc_pt()->at(mc_idx), 
                            b.mc_eta()->at(mc_idx), 
                            b.mc_phi()->at(mc_idx));      
          
        bool found_other_particles = false;
        
        for (size_t mc_idx_2 = 0; mc_idx_2 < b.mc_pt()->size(); mc_idx_2++) {     
          
          double isocone_bbgg_b = 0.1;
          bitset<15> mc_statusFlags2(b.mc_statusflag()->at(mc_idx_2));
          compPart.SetPtEtaPhi(b.mc_pt()->at(mc_idx_2), 
                               b.mc_eta()->at(mc_idx_2), 
                               b.mc_phi()->at(mc_idx_2));

          if (b.mc_id()->at(mc_idx_2) == 22) isocone_bbgg_b = 0.29; // 0.29 set instead of 0.3 because a lot of b were not meeting the criteria just because a photon was 0.29.. delta r away.

          if ((compPart.Pt() >= 10.0) && (mc_idx != mc_idx_2) && (mc_statusFlags2[7] && genb.DeltaR(compPart) < isocone_bbgg_b)) {
            found_other_particles = true;
            break;
          }     
        }

        if (!found_other_particles) number_of_b += 1;
 
      } 
    }
  }

  return number_of_b;

});


const NamedFunc weight_numgenpart("weight_numgenpart",[](const Baby &b) -> NamedFunc::ScalarType{
  
  std::smatch matches_bkg;
  int number_of_jets = 0;
  int number_of_photons= 0;
  double weight = 1;
  bool isdiph = false;

  std::regex_search(*b.FileNames().begin(),matches_bkg,std::regex("unskimmed/pico_(.*)JetsBox"));
  
  if (matches_bkg[1].str() == "DiPhoton"){
    isdiph = true;
  }

  if (isdiph){

    for (size_t mc_idx = 0; mc_idx<b.mc_pt()->size(); mc_idx++) {
      bitset<15> mc_statusFlags(b.mc_statusflag()->at(mc_idx));
      if((mc_statusFlags[7]) && b.mc_pt()->at(mc_idx) > 10){  // Which are isHardProcess
        if (b.mc_id()->at(mc_idx) == 5 || b.mc_id()->at(mc_idx) == -5 || b.mc_id()->at(mc_idx) == 4 || b.mc_id()->at(mc_idx) == -4 || b.mc_id()->at(mc_idx) == 3 || b.mc_id()->at(mc_idx) == -3 || b.mc_id()->at(mc_idx) == 2 || b.mc_id()->at(mc_idx) == -2 || b.mc_id()->at(mc_idx) == 1 || b.mc_id()->at(mc_idx) == -1 || b.mc_id()->at(mc_idx) == 21 )  number_of_jets += 1;
        else if (b.mc_id()->at(mc_idx) == 22) number_of_photons += 1;
      } 
    }
    
    if (number_of_photons == 2 && number_of_jets == 0) {
      weight = 1;  
    }
    else if (number_of_photons == 2 && number_of_jets == 1){ 
      weight = 0.5; 
    }
    else if (number_of_photons == 2 && number_of_jets == 2){ 
      weight = 0.2; 
    }
    else if (number_of_photons == 2 && number_of_jets == 3){
      weight = 0.1; 
    }
    else {
      weight = 1; 
    }

  }

  return weight;
});



map<int, vector<float>> btag_df_wpts{
  {2016, vector<float>({0.0614, 0.3093, 0.7221})},
  {2017, vector<float>({0.0521, 0.3033, 0.7489})},
  {2018, vector<float>({0.0494, 0.2770, 0.7264})},
  {2022, vector<float>({0.0494, 0.2770, 0.7264})},
  {2023, vector<float>({0.0494, 0.2770, 0.7264})}
};

const NamedFunc DeepFlavBool_med("DeepFlavBool_med",[](const Baby &b) -> NamedFunc::ScalarType{
  
  std::string year_string = "";

  if (b.SampleTypeString().Contains("2016") || b.SampleTypeString().Contains("2016APV")) year_string = "2016";
  if (b.SampleTypeString().Contains("2017")) year_string = "2017";
  if (b.SampleTypeString().Contains("2018")) year_string = "2018";

  double deepflav_med = btag_df_wpts[stoi(year_string)][1]; 

  bool DeepFlav = false;
  if (b.bb_ileadscorejet()->size() > 0 && b.bb_isubscorejet()->size() > 0) 
    DeepFlav = b.jet_deepflav()->at(b.bb_ileadscorejet()->at(0)) >= deepflav_med && b.jet_deepflav()->at(b.bb_isubscorejet()->at(0)) >= deepflav_med;
  return DeepFlav;
});

const NamedFunc DeepFlavBool_loose("DeepFlavBool_loose",[](const Baby &b) -> NamedFunc::ScalarType{
  std::string year_string = "";

  if (b.SampleTypeString().Contains("2016") || b.SampleTypeString().Contains("2016APV")) year_string = "2016";
  if (b.SampleTypeString().Contains("2017")) year_string = "2017";
  if (b.SampleTypeString().Contains("2018")) year_string = "2018";
  double deepflav_loose = btag_df_wpts[stoi(year_string)][0];
  bool DeepFlav = false;
  if (b.bb_ileadscorejet()->size() > 0 && b.bb_isubscorejet()->size() > 0) 
    DeepFlav = b.jet_deepflav()->at(b.bb_ileadscorejet()->at(0)) >= deepflav_loose && b.jet_deepflav()->at(b.bb_isubscorejet()->at(0)) >= deepflav_loose;
  return DeepFlav;
});

const NamedFunc DeepFlavBool_tight("DeepFlavBool_tight",[](const Baby &b) -> NamedFunc::ScalarType{
  std::string year_string = "";

  if (b.SampleTypeString().Contains("2016") || b.SampleTypeString().Contains("2016APV")) year_string = "2016";
  if (b.SampleTypeString().Contains("2017")) year_string = "2017";
  if (b.SampleTypeString().Contains("2018")) year_string = "2018";
  double deepflav_tight = btag_df_wpts[stoi(year_string)][2];
  bool DeepFlav = false;
  if (b.bb_ileadscorejet()->size() > 0 && b.bb_isubscorejet()->size() > 0) 
    DeepFlav = b.jet_deepflav()->at(b.bb_ileadscorejet()->at(0)) >= deepflav_tight && b.jet_deepflav()->at(b.bb_isubscorejet()->at(0)) >= deepflav_tight;
  return DeepFlav;
});

const NamedFunc DeepFlavBool_tightmed("DeepFlavBool_tightmed",[](const Baby &b) -> NamedFunc::ScalarType{
  std::string year_string = "";

  if (b.SampleTypeString().Contains("2016") || b.SampleTypeString().Contains("2016APV")) year_string = "2016";
  if (b.SampleTypeString().Contains("2017")) year_string = "2017";
  if (b.SampleTypeString().Contains("2018")) year_string = "2018";
  double deepflav_tight = btag_df_wpts[stoi(year_string)][2];
  double deepflav_med = btag_df_wpts[stoi(year_string)][1]; 
  bool DeepFlav = false;
  if (b.bb_ileadscorejet()->size() > 0 && b.bb_isubscorejet()->size() > 0) 
    DeepFlav = (b.jet_deepflav()->at(b.bb_ileadscorejet()->at(0)) >= deepflav_tight && b.jet_deepflav()->at(b.bb_isubscorejet()->at(0)) >= deepflav_med);
  return DeepFlav;
});

const NamedFunc DeepFlavBool_tightloose("DeepFlavBool_tightloose",[](const Baby &b) -> NamedFunc::ScalarType{
  std::string year_string = "";

  if (b.SampleTypeString().Contains("2016") || b.SampleTypeString().Contains("2016APV")) year_string = "2016";
  if (b.SampleTypeString().Contains("2017")) year_string = "2017";
  if (b.SampleTypeString().Contains("2018")) year_string = "2018";
  double deepflav_tight = btag_df_wpts[stoi(year_string)][2];
  double deepflav_med = btag_df_wpts[stoi(year_string)][1]; 
  bool DeepFlav = false;
  if (b.bb_ileadscorejet()->size() > 0 && b.bb_isubscorejet()->size() > 0) 
    DeepFlav = (b.jet_deepflav()->at(b.bb_ileadscorejet()->at(0)) >= deepflav_tight && b.jet_deepflav()->at(b.bb_isubscorejet()->at(0)) < deepflav_med);
  return DeepFlav;
});

const NamedFunc DeepFlavBool_leadnm1("DeepFlavBool_leadnm1",[](const Baby &b) -> NamedFunc::ScalarType{
  std::string year_string = "";

  if (b.SampleTypeString().Contains("2016") || b.SampleTypeString().Contains("2016APV")) year_string = "2016";
  if (b.SampleTypeString().Contains("2017")) year_string = "2017";
  if (b.SampleTypeString().Contains("2018")) year_string = "2018";
  double deepflav_med = btag_df_wpts[stoi(year_string)][1]; 
  bool DeepFlav = false;
  if (b.bb_ileadscorejet()->size() > 0 && b.bb_isubscorejet()->size() > 0) 
    DeepFlav = (b.jet_deepflav()->at(b.bb_isubscorejet()->at(0)) >= deepflav_med);
  return DeepFlav;
});

const NamedFunc DeepFlavBool_subleadnm1("DeepFlavBool_subleadnm1",[](const Baby &b) -> NamedFunc::ScalarType{
  std::string year_string = "";

  if (b.SampleTypeString().Contains("2016") || b.SampleTypeString().Contains("2016APV")) year_string = "2016";
  if (b.SampleTypeString().Contains("2017")) year_string = "2017";
  if (b.SampleTypeString().Contains("2018")) year_string = "2018";
  double deepflav_med = btag_df_wpts[stoi(year_string)][1]; 
  bool DeepFlav = false;
  if (b.bb_ileadscorejet()->size() > 0 && b.bb_isubscorejet()->size() > 0) 
    DeepFlav = b.jet_deepflav()->at(b.bb_ileadscorejet()->at(0)) >= deepflav_med;
  return DeepFlav;
});

/* ///Old deepflav functions for only 2016 mc data comp
const NamedFunc DeepFlavBool_med("DeepFlavBool_med",[](const Baby &b) -> NamedFunc::ScalarType{
  bool DeepFlav = false;
  if (b.bb_ileadscorejet()->size() > 0 && b.bb_isubscorejet()->size() > 0) 
    DeepFlav = b.jet_deepflav()->at(b.bb_ileadscorejet()->at(0)) >= 0.3093 && b.jet_deepflav()->at(b.bb_isubscorejet()->at(0)) >= 0.3093;
  return DeepFlav;
});

const NamedFunc DeepFlavBool_loose("DeepFlavBool_loose",[](const Baby &b) -> NamedFunc::ScalarType{
  bool DeepFlav = false;
  if (b.bb_ileadscorejet()->size() > 0 && b.bb_isubscorejet()->size() > 0) 
    DeepFlav = b.jet_deepflav()->at(b.bb_ileadscorejet()->at(0)) >= 0.0614 && b.jet_deepflav()->at(b.bb_isubscorejet()->at(0)) >= 0.0614;
  return DeepFlav;
});

const NamedFunc DeepFlavBool_tight("DeepFlavBool_tight",[](const Baby &b) -> NamedFunc::ScalarType{
  bool DeepFlav = false;
  if (b.bb_ileadscorejet()->size() > 0 && b.bb_isubscorejet()->size() > 0) 
    DeepFlav = b.jet_deepflav()->at(b.bb_ileadscorejet()->at(0)) >= 0.7221 && b.jet_deepflav()->at(b.bb_isubscorejet()->at(0)) >= 0.7221;
  return DeepFlav;
});

const NamedFunc DeepFlavBool_tightmed("DeepFlavBool_tightmed",[](const Baby &b) -> NamedFunc::ScalarType{
  bool DeepFlav = false;
  if (b.bb_ileadscorejet()->size() > 0 && b.bb_isubscorejet()->size() > 0) 
    DeepFlav = (b.jet_deepflav()->at(b.bb_ileadscorejet()->at(0)) >= 0.7221 && b.jet_deepflav()->at(b.bb_isubscorejet()->at(0)) >= 0.3093) || (b.jet_deepflav()->at(b.bb_ileadscorejet()->at(0)) >= 0.3093 && b.jet_deepflav()->at(b.bb_isubscorejet()->at(0)) >= 0.7221);
  return DeepFlav;
});

*/

const NamedFunc DeepFlavScore_lead("DeepFlavScore_lead",[](const Baby &b) -> NamedFunc::ScalarType{
  return b.jet_deepflav()->at(b.bb_ileadscorejet()->at(0));
});

const NamedFunc DeepFlavScore_sublead("DeepFlavScore_sublead",[](const Baby &b) -> NamedFunc::ScalarType{
  return b.jet_deepflav()->at(b.bb_isubscorejet()->at(0));
});

const NamedFunc jet_pt_lead("jet_pt_lead",[](const Baby &b) -> NamedFunc::ScalarType{
  return b.jet_pt()->at(b.bb_ileadscorejet()->at(0));
});

const NamedFunc jet_pt_sublead("jet_pt_sublead",[](const Baby &b) -> NamedFunc::ScalarType{
  return b.jet_pt()->at(b.bb_isubscorejet()->at(0));
});

const NamedFunc jet_eta_lead("jet_eta_lead",[](const Baby &b) -> NamedFunc::ScalarType{
  return b.jet_eta()->at(b.bb_ileadscorejet()->at(0));
});

const NamedFunc jet_eta_sublead("jet_eta_sublead",[](const Baby &b) -> NamedFunc::ScalarType{
  return b.jet_eta()->at(b.bb_isubscorejet()->at(0));
});

const NamedFunc jet_phi_lead("jet_phi_lead",[](const Baby &b) -> NamedFunc::ScalarType{
  return b.jet_phi()->at(b.bb_ileadscorejet()->at(0));
});

const NamedFunc jet_phi_sublead("jet_phi_sublead",[](const Baby &b) -> NamedFunc::ScalarType{
  return b.jet_phi()->at(b.bb_isubscorejet()->at(0));
});

const NamedFunc DeepFlavScore_sum("DeepFlavScore_sum",[](const Baby &b) -> NamedFunc::ScalarType{
  return b.jet_deepflav()->at(b.bb_ileadscorejet()->at(0)) + b.jet_deepflav()->at(b.bb_isubscorejet()->at(0));
});

const NamedFunc Photon_pt_lead("Photon_pt_lead",[](const Baby &b) -> NamedFunc::ScalarType{
  return b.photon_pt()->at(0);
  //return b.photon_pt()->at(0)/b.photonphoton_m()->at(0);
});

const NamedFunc Photon_pt_sublead("Photon_pt_sublead",[](const Baby &b) -> NamedFunc::ScalarType{
  return b.photon_pt()->at(1);
  //return b.photon_pt()->at(1)/b.photonphoton_m()->at(0);
});

const NamedFunc Photon_eta_lead("Photon_eta_lead",[](const Baby &b) -> NamedFunc::ScalarType{
  return b.photon_eta()->at(0);
  
});

const NamedFunc Photon_eta_sublead("Photon_eta_sublead",[](const Baby &b) -> NamedFunc::ScalarType{
  return b.photon_eta()->at(1);
  
});

const NamedFunc Photon_phi_lead("Photon_phi_lead",[](const Baby &b) -> NamedFunc::ScalarType{
  return b.photon_phi()->at(0);
  
});

const NamedFunc Photon_phi_sublead("Photon_phi_sublead",[](const Baby &b) -> NamedFunc::ScalarType{
  return b.photon_phi()->at(1);
  
});

const NamedFunc Photon_id_lead("Photon_id_lead",[](const Baby &b) -> NamedFunc::ScalarType{
  return b.photon_idmva()->at(0);
  //return b.photon_pt()->at(0)/b.photonphoton_m()->at(0);
});

const NamedFunc Photon_id_sublead("Photon_id_sublead",[](const Baby &b) -> NamedFunc::ScalarType{
  return b.photon_idmva()->at(1);
  //return b.photon_pt()->at(0)/b.photonphoton_m()->at(0);
});


const NamedFunc Photon_p_lead("Photon_p_lead",[](const Baby &b) -> NamedFunc::ScalarType{
  return b.photon_pt()->at(0) * sinh(b.photon_eta()->at(0));
  //return b.photon_pt()->at(0)/b.photonphoton_m()->at(0);
});

const NamedFunc Photon_p_sublead("Photon_p_sublead",[](const Baby &b) -> NamedFunc::ScalarType{
  return b.photon_pt()->at(1) * sinh(b.photon_eta()->at(1));
  //return b.photon_pt()->at(0)/b.photonphoton_m()->at(0);
});



const NamedFunc Photon_pt_thresh("Photon_pt_thresh",[](const Baby &b) -> NamedFunc::ScalarType{
  bool thresh = false;
  if (b.nphoton() >= 2) thresh = (b.photon_pt()->at(0) > 35 && b.photon_pt()->at(1) > 25 /*22*/ );
  return thresh;
});

const NamedFunc Photon_id80("Photon_id80",[](const Baby &b) -> NamedFunc::ScalarType{
  bool thresh = false;
  if (b.nphoton() >= 2) thresh = (b.photon_id80()->at(0) && b.photon_id80()->at(1));
  return thresh;
});

const NamedFunc Photon_id80_lead("Photon_id80_lead",[](const Baby &b) -> NamedFunc::ScalarType{
  bool thresh = false;
  if (b.nphoton() >= 2) thresh = (b.photon_id80()->at(0) && !b.photon_id80()->at(1));
  return thresh;
});

const NamedFunc Photon_idlessthan80("Photon_idlessthan80",[](const Baby &b) -> NamedFunc::ScalarType{
  bool thresh = false;
  if (b.nphoton() >= 2) thresh = (!b.photon_id80()->at(0) && !b.photon_id80()->at(1));
  return thresh;
});

const NamedFunc bjet_pt_thresh("bjet_pt_thresh",[](const Baby &b) -> NamedFunc::ScalarType{
  bool thresh = false;
  if (b.nbb() == 1) thresh = (b.jet_pt()->at(b.bb_ileadscorejet()->at(0)) > 40 && b.jet_pt()->at(b.bb_isubscorejet()->at(0)) > 40);
  return thresh;
});

const NamedFunc nearestjet_met_dphi("nearestjet_met_dphi",[](const Baby &b) -> NamedFunc::ScalarType{
  double dphi = 5;
  for (size_t iIdx=0; iIdx < b.jet_pt()->size(); ++iIdx) {     
    if (b.jet_isgood()->at(iIdx) && b.jet_deepflav()->at(iIdx) > 0.3093) {
      if (b.jet_met_dphi()->at(iIdx) < dphi) dphi = b.jet_met_dphi()->at(iIdx);        
    } 
  }
  return dphi;
});

const NamedFunc nearestjet_met_pt("nearestjet_met_pt",[](const Baby &b) -> NamedFunc::ScalarType{
  double dphi = 5;
  double pt = 0;
  for (size_t iIdx=0; iIdx < b.jet_pt()->size(); ++iIdx) {     
    if (b.jet_isgood()->at(iIdx) && b.jet_deepflav()->at(iIdx) > 0.3093) {
      if (b.jet_met_dphi()->at(iIdx) < dphi) {
       dphi = b.jet_met_dphi()->at(iIdx);
       pt = b.jet_pt()->at(iIdx);
      }        
    } 
  }
  return pt;
});
  

const NamedFunc jet_pom_lead("jet_pom_lead",[](const Baby &b) -> NamedFunc::ScalarType{
  return b.jet_pt()->at(b.bb_ileadscorejet()->at(0))/b.bb_m()->at(0);
});

const NamedFunc jet_pom_sublead("jet_pom_sublead",[](const Baby &b) -> NamedFunc::ScalarType{
  return b.jet_pt()->at(b.bb_isubscorejet()->at(0))/b.bb_m()->at(0);
});

const NamedFunc photon_pom_lead("photon_pom_lead",[](const Baby &b) -> NamedFunc::ScalarType{
  return b.photon_pt()->at(0)/b.photonphoton_m()->at(0);
});

const NamedFunc photon_pom_sublead("photon_pom_sublead",[](const Baby &b) -> NamedFunc::ScalarType{
  return b.photon_pt()->at(1)/b.photonphoton_m()->at(0);
});

const NamedFunc dau_jet_met_dphi_sublead("dau_jet_met_dphi_sublead",[](const Baby &b) -> NamedFunc::ScalarType{
  return DeltaPhi(b.jet_phi()->at(b.bb_isubscorejet()->at(0)), b.met_phi());
});
const NamedFunc dau_jet_met_dphi_lead("dau_jet_met_dphi_lead",[](const Baby &b) -> NamedFunc::ScalarType{
  return DeltaPhi(b.jet_phi()->at(b.bb_ileadscorejet()->at(0)), b.met_phi());
});

const NamedFunc DeltaPhi_METak4blead("DeltaPhi_METak4b",[](const Baby &b) -> NamedFunc::ScalarType{
  return DeltaPhi(b.jet_phi()->at(b.bb_ileadscorejet()->at(0)), b.met_phi());
});

const NamedFunc DeltaPhi_METak4bsublead("DeltaPhi_METak4b",[](const Baby &b) -> NamedFunc::ScalarType{
  return DeltaPhi(b.jet_phi()->at(b.bb_isubscorejet()->at(0)), b.met_phi());
});

const NamedFunc DeltaPhi_HH("DeltaPhi_HH",[](const Baby &b) -> NamedFunc::ScalarType{
  return DeltaPhi(b.photonphoton_phi()->at(0), b.bb_phi()->at(0));
});

const NamedFunc DeltaR_HH("DeltaR_HH",[](const Baby &b) -> NamedFunc::ScalarType{
  return dR(b.photonphoton_eta()->at(0), b.bb_eta()->at(0), b.photonphoton_phi()->at(0), b.bb_phi()->at(0));
});


const NamedFunc mt2("mt2",[](const Baby &b) -> NamedFunc::ScalarType{
  double MT2 =  asymm_mt2_lester_bisect::get_mT2(
             b.photonphoton_m()->at(0), b.photonphoton_pt()->at(0)*cos(b.photonphoton_phi()->at(0)), b.photonphoton_pt()->at(0)*sin(b.photonphoton_phi()->at(0)),
             b.bb_m()->at(0), b.bb_pt()->at(0)*cos(b.bb_phi()->at(0)), b.bb_pt()->at(0)*sin(b.bb_phi()->at(0)),
             b.met()*cos(b.met_phi()), b.met()*sin(b.met_phi()),
             1, 1,
             0);
  return MT2;
});

const NamedFunc higgs_gg_lab_boost("higgs_gg_lab_boost",[](const Baby &b) -> NamedFunc::ScalarType{
  double gg_pt = b.photonphoton_pt()->at(0);
  double gg_pz = b.photonphoton_pt()->at(0)*tan(b.photonphoton_eta()->at(0));
  double gg_m = b.photonphoton_m()->at(0);
  double gg_boost = sqrt(gg_pt*gg_pt + gg_pz*gg_pz)/sqrt(gg_pt*gg_pt + gg_pz*gg_pz + gg_m*gg_m);
  return gg_boost;
});

const NamedFunc higgs_bb_lab_boost("higgs_bb_lab_boost",[](const Baby &b) -> NamedFunc::ScalarType{
  double bb_pt = b.bb_pt()->at(0);
  double bb_pz = b.bb_pt()->at(0)*tan(b.bb_eta()->at(0));
  double bb_m = b.bb_m()->at(0);
  double bb_boost = sqrt(bb_pt*bb_pt + bb_pz*bb_pz)/sqrt(bb_pt*bb_pt + bb_pz*bb_pz + bb_m*bb_m);
  return bb_boost;
});


const NamedFunc gg_pz("gg_pz",[](const Baby &b) -> NamedFunc::ScalarType{
  double h_pz = b.photonphoton_pt()->at(0)*tan(b.photonphoton_eta()->at(0));
  return h_pz;
});

const NamedFunc bb_pz("bb_pz",[](const Baby &b) -> NamedFunc::ScalarType{
  double h_pz = b.bb_pt()->at(0)*tan(b.bb_eta()->at(0));
  return h_pz;
});

const NamedFunc gg_p("gg_p",[](const Baby &b) -> NamedFunc::ScalarType{
  double h_p = b.photonphoton_pt()->at(0)/cos(b.photonphoton_eta()->at(0));
  return h_p;
});

const NamedFunc bb_p("bb_p",[](const Baby &b) -> NamedFunc::ScalarType{
  double h_p = b.bb_pt()->at(0)/cos(b.bb_eta()->at(0));
  return h_p;
});

const NamedFunc bb_pmet("bb_pmet",[](const Baby &b) -> NamedFunc::ScalarType{
  double tot_x = b.bb_pt()->at(0)*cos(b.bb_phi()->at(0)) + b.met()*cos(b.met_phi());
  double tot_y = b.bb_pt()->at(0)*sin(b.bb_phi()->at(0)) + b.met()*sin(b.met_phi());
  return sqrt(tot_x*tot_x + tot_y*tot_y);
});



const NamedFunc cos_theta_CS("cos_theta_CS",[](const Baby &b) -> NamedFunc::ScalarType{

  // daughter b jet 1
  double b1_pt = b.photon_pt()->at(0);
  double b1_pz = b1_pt*sinh(b.photon_eta()->at(0));
  double b1_m = 0;
  double b1_p = sqrt(b1_pt*b1_pt + b1_pz * b1_pz);
  double b1_E = sqrt(b1_p*b1_p + b1_m * b1_m);
 
  double b2_pt = b.photon_pt()->at(1);
  double b2_pz = b2_pt*sinh(b.photon_eta()->at(1));
  double b2_m = 0;
  double b2_p = sqrt(b2_pt*b2_pt + b2_pz * b2_pz);
  double b2_E = sqrt(b2_p*b2_p + b2_m * b2_m);

  // higgs 
  double higgs_pt = b.photonphoton_pt()->at(0);
  double higgs_m = b.photonphoton_m()->at(0);


  double costhetCS = 2*(b1_pz*b2_E - b2_pz*b1_E)/(higgs_m*sqrt(higgs_m*higgs_m + higgs_pt*higgs_pt));
  return costhetCS;

});




const NamedFunc cos_theta_Jeff_lead("cos_theta_Jeff",[](const Baby &b) -> NamedFunc::ScalarType{

  // daughter b jet 1
  double b1_pt = b.photon_pt()->at(0);
  double b1_px = b1_pt*cos(b.photon_phi()->at(0));
  double b1_py = b1_pt*sin(b.photon_phi()->at(0));
  double b1_pz = b1_pt*sinh(b.photon_eta()->at(0));
  double b1_m = 0;
  double b1_p = sqrt(b1_pt*b1_pt + b1_pz * b1_pz);
  double b1_E = sqrt(b1_p*b1_p + b1_m * b1_m);
 
  // higgs 
  double higgs_pt = b.photonphoton_pt()->at(0);
  double higgs_phi = b.photonphoton_phi()->at(0);  
  double higgs_eta = b.photonphoton_eta()->at(0); 
  double higgs_px = higgs_pt*cos(higgs_phi);
  double higgs_py = higgs_pt*sin(higgs_phi);
  double higgs_pz = higgs_pt*sinh(higgs_eta);
  double higgs_m = b.photonphoton_m()->at(0);
  double higgs_p = sqrt(higgs_pt*higgs_pt + higgs_pz * higgs_pz);
  double higgs_E = sqrt(higgs_p*higgs_p + higgs_m * higgs_m);

  // calculating boost from lab to rest frame
  double b_labtorest = higgs_p/higgs_E;
  double gamma = 1/sqrt(1 - b_labtorest * b_labtorest);

  // Define this boost direction to be z axis, new components for b jets along this direction
  double polarangle = acos((b1_px*higgs_px + b1_py*higgs_py + b1_pz*higgs_pz)/(b1_p*higgs_p));
  double b1_pt_ah = b1_p * sin(polarangle);
  double b1_pz_ah = b1_p * cos(polarangle);  

  // boosting daughter b jets to higgs rest frame along this new z direction
  double b1_pt_rest = b1_pt_ah;
  double b1_pz_rest = gamma * b1_pz_ah - b_labtorest * gamma * b1_E;

  double costhet = b1_pz_rest/sqrt(b1_pz_rest*b1_pz_rest + b1_pt_rest*b1_pt_rest);
  //(void)costhet;
  return costhet;
});

const NamedFunc cos_theta_Jeff_sublead("cos_theta_Jeff",[](const Baby &b) -> NamedFunc::ScalarType{

  // daughter b jet 1
  double b1_pt = b.photon_pt()->at(1);
  double b1_px = b1_pt*cos(b.photon_phi()->at(1));
  double b1_py = b1_pt*sin(b.photon_phi()->at(1));
  double b1_pz = b1_pt*sinh(b.photon_eta()->at(1));
  double b1_m = 0;
  double b1_p = sqrt(b1_pt*b1_pt + b1_pz * b1_pz);
  double b1_E = sqrt(b1_p*b1_p + b1_m * b1_m);
 
  // higgs 
  double higgs_pt = b.photonphoton_pt()->at(0);
  double higgs_phi = b.photonphoton_phi()->at(0);  
  double higgs_eta = b.photonphoton_eta()->at(0); 
  double higgs_px = higgs_pt*cos(higgs_phi);
  double higgs_py = higgs_pt*sin(higgs_phi);
  double higgs_pz = higgs_pt*sinh(higgs_eta);
  double higgs_m = b.photonphoton_m()->at(0);
  double higgs_p = sqrt(higgs_pt*higgs_pt + higgs_pz * higgs_pz);
  double higgs_E = sqrt(higgs_p*higgs_p + higgs_m * higgs_m);

  // calculating boost from lab to rest frame
  double b_labtorest = higgs_p/higgs_E;
  double gamma = 1/sqrt(1 - b_labtorest * b_labtorest);

  // Define this boost direction to be z axis, new components for b jets along this direction
  double polarangle = acos((b1_px*higgs_px + b1_py*higgs_py + b1_pz*higgs_pz)/(b1_p*higgs_p));
  double b1_pt_ah = b1_p * sin(polarangle);
  double b1_pz_ah = b1_p * cos(polarangle);  

  // boosting daughter b jets to higgs rest frame along this new z direction
  double b1_pt_rest = b1_pt_ah;
  double b1_pz_rest = gamma * b1_pz_ah - b_labtorest * gamma * b1_E;

  double costhet = b1_pz_rest/sqrt(b1_pz_rest*b1_pz_rest + b1_pt_rest*b1_pt_rest);
  //(void)costhet;
  return costhet;
});


const NamedFunc cos_theta_Jeff("cos_theta_Jeff",[](const Baby &b) -> NamedFunc::VectorType{

  double b1_pt = b.photon_pt()->at(0);
  double b1_px = b1_pt*cos(b.photon_phi()->at(0));
  double b1_py = b1_pt*sin(b.photon_phi()->at(0));
  double b1_pz = b1_pt*sinh(b.photon_eta()->at(0));
  double b1_m = 0;
  double b1_p = sqrt(b1_pt*b1_pt + b1_pz * b1_pz);
  double b1_E = sqrt(b1_p*b1_p + b1_m * b1_m);
 
  // daughter b jet 2
  double b2_pt = b.photon_pt()->at(1);
  double b2_px = b2_pt*cos(b.photon_phi()->at(1));
  double b2_py = b2_pt*sin(b.photon_phi()->at(1));
  double b2_pz = b2_pt*sinh(b.photon_eta()->at(1));
  double b2_m = 0;
  double b2_p = sqrt(b2_pt*b2_pt + b2_pz * b2_pz);
  double b2_E = sqrt(b2_p*b2_p + b2_m * b2_m);

  // higgs 
  double higgs_pt = b.photonphoton_pt()->at(0);
  double higgs_phi = b.photonphoton_phi()->at(0);  
  double higgs_eta = b.photonphoton_eta()->at(0); 
  double higgs_px = higgs_pt*cos(higgs_phi);
  double higgs_py = higgs_pt*sin(higgs_phi);
  double higgs_pz = higgs_pt*sinh(higgs_eta);
  double higgs_m = b.photonphoton_m()->at(0);
  double higgs_p = sqrt(higgs_pt*higgs_pt + higgs_pz * higgs_pz);
  double higgs_E = sqrt(higgs_p*higgs_p + higgs_m * higgs_m);

  // calculating boost from lab to rest frame
  double b_labtorest = higgs_p/higgs_E;
  double gamma = 1/sqrt(1 - b_labtorest * b_labtorest);



  // Define this boost direction to be z axis, new components for b jets along this direction
  double polarangle1 = acos((b1_px*higgs_px + b1_py*higgs_py + b1_pz*higgs_pz)/(b1_p*higgs_p));
  double b1_pt_ah = b1_p * sin(polarangle1);
  double b1_pz_ah = b1_p * cos(polarangle1);  

  // boosting daughter b jets to higgs rest frame along this new z direction
  double b1_pt_rest = b1_pt_ah;
  double b1_pz_rest = gamma * b1_pz_ah - b_labtorest * gamma * b1_E;

  double costhet_lead = b1_pz_rest/sqrt(b1_pz_rest*b1_pz_rest + b1_pt_rest*b1_pt_rest);


  // Define this boost direction to be z axis, new components for b jets along this direction
  double polarangle2 = acos((b2_px*higgs_px + b2_py*higgs_py + b2_pz*higgs_pz)/(b2_p*higgs_p));
  double b2_pt_ah = b2_p * sin(polarangle2);
  double b2_pz_ah = b2_p * cos(polarangle2);  

  // boosting daughter b jets to higgs rest frame along this new z direction
  double b2_pt_rest = b2_pt_ah;
  double b2_pz_rest = gamma * b2_pz_ah - b_labtorest * gamma * b2_E;

  double costhet_sublead = b2_pz_rest/sqrt(b2_pz_rest*b2_pz_rest + b2_pt_rest*b2_pt_rest);



  std::vector<double> costhet {costhet_lead,costhet_sublead};
  //(void)costhet;
  return costhet;
});


const NamedFunc cos_theta_Jeff_truth("cos_theta_Jeff_truth",[](const Baby &b) -> NamedFunc::ScalarType{

  double costhet = -2;

  for (size_t iIdx=0; iIdx < b.mc_id()->size(); ++iIdx) {
    
    if ((b.mc_id()->at(iIdx) == 5 && b.mc_momidx()->at(iIdx) != -1 && b.mc_id()->at(b.mc_momidx()->at(iIdx)) == 25)) {
    
      double b1_pt = b.mc_pt()->at(iIdx);
      double b1_px = b1_pt*cos(b.mc_phi()->at(iIdx));
      double b1_py = b1_pt*sin(b.mc_phi()->at(iIdx));
      double b1_pz = b1_pt*sinh(b.mc_eta()->at(iIdx));
      double b1_m = b.mc_mass()->at(iIdx);
      double b1_p = sqrt(b1_pt*b1_pt + b1_pz * b1_pz);
      double b1_E = sqrt(b1_p*b1_p + b1_m * b1_m);  

      double higgs_pt = b.mc_pt()->at(b.mc_momidx()->at(iIdx));
      double higgs_px = higgs_pt*cos(b.mc_phi()->at(b.mc_momidx()->at(iIdx)));
      double higgs_py = higgs_pt*sin(b.mc_phi()->at(b.mc_momidx()->at(iIdx)));
      double higgs_pz = higgs_pt*sinh(b.mc_eta()->at(b.mc_momidx()->at(iIdx)));
      double higgs_m = b.mc_mass()->at(b.mc_momidx()->at(iIdx));
      double higgs_p = sqrt(higgs_pt*higgs_pt + higgs_pz * higgs_pz);
      double higgs_E = sqrt(higgs_p*higgs_p + higgs_m * higgs_m);

      // calculating boost from lab to rest frame
      double b_labtorest = higgs_p/higgs_E;
      double gamma = 1/sqrt(1 - b_labtorest * b_labtorest);

      // Define this boost direction to be z axis, new components for b jets along this direction
      double polarangle = acos((b1_px*higgs_px + b1_py*higgs_py + b1_pz*higgs_pz)/(b1_p*higgs_p));
      double b1_pt_ah = b1_p * sin(polarangle);
      double b1_pz_ah = b1_p * cos(polarangle);  

      // boosting daughter b jets to higgs rest frame along this new z direction
      double b1_pt_rest = b1_pt_ah;
      double b1_pz_rest = gamma * b1_pz_ah - b_labtorest * gamma * b1_E;

      costhet = b1_pz_rest/sqrt(b1_pz_rest*b1_pz_rest + b1_pt_rest*b1_pt_rest);
    
    }
  
  }

  return costhet;
});

const NamedFunc genb_counter("genb_counter",[](const Baby &b) -> NamedFunc::ScalarType{
  int n = 0;
  if(b.no_genb()) n =0;
  else if (b.one_genb()) n=1;
  else if (b.two_genb()) n=2;
  else n=3;
  return n;
});

const NamedFunc process_bin("process_bin",[](const Baby &b) -> NamedFunc::ScalarType{
  int n = 0;
  double weight = roundDouble(b.mc_weight(), 2);
  if(weight == 1) n =2;
  else if (weight == 0.5) n=3;
  else if (weight == 0.2) n=4;
  else if (weight == 0.1) n=5;

  return n;
});

const NamedFunc geng_counter("geng_counter",[](const Baby &b) -> NamedFunc::ScalarType{
  int n = 0;
  if(b.use_event()) n =1;
  else n=2;

  return n;
});


const NamedFunc overlapDiph_cutflow("overlap_cutflow",[](const Baby &b) -> NamedFunc::ScalarType{
  
  std::string file = *b.FileNames().begin();
  if (file.find("DiPhotonJetsBox_") != file.npos) return b.no_genb();
  else if (file.find("DiPhotonJetsBox1BJet") != file.npos) return b.one_genb();
  else if (file.find("DiPhotonJetsBox2BJets") != file.npos) return b.two_genb();
  else return true;

});

const NamedFunc overlapGjet_cutflow("overlap_cutflow",[](const Baby &b) -> NamedFunc::ScalarType{
  
  std::string file = *b.FileNames().begin();
  if (file.find("GJet_Pt") != file.npos) return b.use_event();
  else return true;

});



int main(){

  gErrorIgnoreLevel = 6000;

  /*
  const NamedFunc w_lumi("w_lumi", [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleTypeString().Contains("-")) {
      return 1.;
    }
    else {
      return b.w_lumi();
    } 
  });
  */
  const NamedFunc wgt("wgt", [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleTypeString().Contains("-")) {
      return 1.;
    }
    else {
      return b.weight();
    } 
  });
  /*
  NamedFunc trig("trig",[](const Baby &b) -> NamedFunc::ScalarType{

      bool trig_year = false;

      if (b.SampleTypeString()=="2016" || b.SampleTypeString()=="2016APV") {
          trig_year = b.HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90();
      } else if (b.SampleTypeString()=="2017" || b.SampleTypeString()=="2018"){
          trig_year = b.HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90();
      }


      return trig_year;
  });
  */

  NamedFunc trig("trig",[](const Baby &b) -> NamedFunc::ScalarType{

      bool trig_year = false;
      if (b.SampleTypeString().Contains("2016")|| b.SampleTypeString().Contains("2016APV")) {
          trig_year = b.HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90();
      } else if (b.SampleTypeString().Contains("2017") || b.SampleTypeString().Contains("2018")){
          trig_year = b.HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90();
      }

      return trig_year;
  });

  NamedFunc w_lumi("w_lumi",[](const Baby &b) -> NamedFunc::ScalarType{
    double w_year = 1;

    if(b.SampleTypeString().Contains("-")) w_year = 1;
    else{

      bool isdiph = false;
      bool isQCD2 = false;
      bool isSMS = false;

      //Regex search takes too long, switching to find()
      //std::smatch matches_bkg;
      //cout<<*b.FileNames().begin()<<"  "<<typeid(*b.FileNames().begin()).name()<<endl;
      //std::regex_search(*b.FileNames().begin(),matches_bkg,std::regex("unskimmed/pico_(.*)JetsBox"));
      //if (matches_bkg[1].str() == "DiPhoton") isdiph = true;
      
      std::string file = *b.FileNames().begin();

      if (file.find("DiPhotonJetsBox") != file.npos) isdiph = true;
      if (file.find("QCD_Pt-40") != file.npos) isQCD2 = true;
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

      if (isSMS) w_year = 137.61; //since we have only 2016 v7 MC

      w_year = w_year*b.w_lumi();
      
      if(isdiph) w_year *= b.mc_weight()>0 ? b.mc_weight():-b.mc_weight(); /*w_year = w_year*b.mc_weight();*/ //w_lumi already knows the sign
      if(isQCD2) w_year *= -0.1131001131; //due to it being unable to find cross section in nano2pico
    }

    return w_year;
  
  });
  
  //w_lumi function without undoing enhance factors
  NamedFunc w_lumi_ef("w_lumi_ef",[](const Baby &b) -> NamedFunc::ScalarType{

    double w_year = 1;
    
    if(b.SampleTypeString().Contains("-")) w_year = 1;
    else{

      bool isdiph = false;
      bool isQCD2 = false;
      bool isSMS = false;
      double neff = 0;
      double nent = 0;

      std::string file = *b.FileNames().begin();
      if (file.find("QCD_Pt-40") != file.npos) isQCD2 = true;
      if (file.find("TChiHH") != file.npos) isSMS = true;

      if (file.find("DiPhotonJetsBox_") != file.npos) {
        isdiph = true;
        if (b.SampleTypeString()=="2016") {
           neff = 6167042.9;
           nent = 14263911;
        } else if (b.SampleTypeString()=="2016APV"){
           neff = 6198395.7;
           nent = 14333896;
        } else if (b.SampleTypeString()=="2017"){
           neff = 13897118;
           nent = 29676800;
        } else if (b.SampleTypeString()=="2018"){
           neff = 14028260;
           nent = 29948772;
        }   
      }

      if (file.find("DiPhotonJetsBox1BJet") != file.npos) {
        isdiph = true;
        if (b.SampleTypeString()=="2016") {
           neff = 171530.9;
           nent = 953898;
        } else if (b.SampleTypeString()=="2016APV"){
           neff = 179747.12;
           nent = 999924;
        } else if (b.SampleTypeString()=="2017"){
           neff = 386021.5;
           nent = 1999832;
        } else if (b.SampleTypeString()=="2018"){
           neff = 386147.4;
           nent = 1999794;
        }   
      }

      if (file.find("DiPhotonJetsBox2BJets") != file.npos) {
        isdiph = true;
        if (b.SampleTypeString()=="2016") {
           neff = 168499.12;
           nent = 999897;
        } else if (b.SampleTypeString()=="2016APV"){
           neff = 165321.67;
           nent = 981922;
        } else if (b.SampleTypeString()=="2017"){
           neff = 340044.31;
           nent = 1999826;
        } else if (b.SampleTypeString()=="2018"){
           neff = 340043.5;
           nent = 1999823;
        }   
      }

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

      w_year = w_year*b.w_lumi();
      
      if(isdiph) w_year *= neff/nent; 
      if(isQCD2) w_year *= -0.1131001131; //due to it being unable to find cross section in nano2pico

    }

    return w_year;
  
  });


  //NamedFunc trigs_2016_gg   = "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90";
  // Setting up Plotmaker

  Palette colors("txt/colors.txt","default");
  Process::Type sig  = Process::Type::signal;
  Process::Type bkg  = Process::Type::background;

  PlotOpt log_lumi("txt/plot_styles.txt","CMSPaper");
  log_lumi.Title(TitleType::info)
          .YAxis(YAxisType::log)
          .Stack(StackType::signal_overlay)
          .Overflow(OverflowType::overflow)
          .YTitleOffset(1.75)
          .AutoYAxis(true)
          .UseCMYK(false)
          .LeftMargin(0.17)
          .LegendColumns(1)
          .CanvasWidth(800)
          .CanvasHeight(800)
          .FileExtensions({"pdf"});

  PlotOpt lin_lumi("txt/plot_styles.txt","CMSPaper");
  lin_lumi.Title(TitleType::info)
          .YAxis(YAxisType::linear)
          .Overflow(OverflowType::none)
          .YTitleOffset(1.75)
          .AutoYAxis(true)
          .UseCMYK(false)
          .LeftMargin(0.17)
          .LegendColumns(1)
          .CanvasWidth(800)
          .CanvasHeight(800)
          .FileExtensions({"pdf"});
  
  PlotOpt style("txt/plot_styles.txt","Eff2D");
  vector<PlotOpt> bkg_hist = {style().Stack(StackType::data_norm)
                                     .YTitleOffset(1.)
                                     .LogMinimum(0.0001)
                                     //.LogMaximum(2000)
                                     .Overflow(OverflowType::both)
                                     .CanvasWidth(600)
                                     .LabelSize(0.04)
                                     .UseCMYK(false)
                                     //.UseROOTColors(true)
                                     //.YAxis(YAxisType::log)
                                     .Title(TitleType::info)};


  PlotOpt style2D("txt/plot_styles.txt","Scatter");
  vector<PlotOpt> bkg_hist_scat = {style2D().Stack(StackType::lumi_shapes)
                                     .Title(TitleType::info)};



  PlotOpt log_stack = log_lumi().Stack(StackType::signal_overlay);
  PlotOpt lin_stack = lin_lumi().Stack(StackType::signal_overlay);
  PlotOpt lin_dnorm = lin_lumi().Stack(StackType::data_norm)
                                .Bottom(BottomType::ratio);
  PlotOpt log_dnorm = log_lumi().Stack(StackType::data_norm)
                                .Bottom(BottomType::ratio);
  
  PlotOpt lin_sigtop = lin_lumi().Stack(StackType::signal_overlay);
  vector<PlotOpt> ops_log = {log_stack};
  vector<PlotOpt> ops_lin = {lin_stack};
  vector<PlotOpt> ops_linlog = {lin_stack, log_stack};
  vector<PlotOpt> ops_lin_dnorm = {lin_dnorm};
  vector<PlotOpt> ops_log_dnorm = {log_dnorm};
  vector<PlotOpt> ops_lin_nodnorm = {lin_sigtop};
  vector<PlotOpt> ops_2D = {bkg_hist};
  vector<PlotOpt> ops_2D_scat = {bkg_hist_scat};
  

  // Sept production baseline
  NamedFunc bin1 = "met >=20 && met <=70";
  NamedFunc bin2 = "met <=130 && met >70";
  NamedFunc bin3 = "met >130";
  NamedFunc pmet = "met >= 0";
  NamedFunc bin1OR2 = bin1 || bin2;
  NamedFunc bin2OR3 = bin2 || bin3;
  NamedFunc novetolepton = "nvmu == 0 && nvel == 0";
  NamedFunc njetcut = "njet < 4";
  NamedFunc mbbcut = "bb_m < 140 && bb_m > 90";
  NamedFunc flavmedcut = DeepFlavBool_med;
  NamedFunc flavloosecut = DeepFlavBool_loose;
  NamedFunc flavtightcut = DeepFlavBool_tight;
  NamedFunc mggcut = "photonphoton_m >= 100 && photonphoton_m <= 200";
  NamedFunc ggdrcut = "photonphoton_dr < 2.5";
  NamedFunc ggptcut = "photonphoton_pt > 75"; 
  NamedFunc ggdrcut_bin2 = "photonphoton_dr < 2.2";
  NamedFunc ggptcut_bin2 = "photonphoton_pt > 110";
  NamedFunc ggdrcut_bin3 = "photonphoton_dr < 2";
  NamedFunc ggptcut_bin3 = "photonphoton_pt > 130";
  NamedFunc costhetacut = cos_theta_Jeff_lead < 0.75;
  NamedFunc mbbwindow = "bb_m <= 140 && bb_m >= 90";
  NamedFunc mbbloosewindow = "bb_m <= 140 && bb_m >= 70";
  NamedFunc mggwindow = "photonphoton_m >= 115 && photonphoton_m <= 135";
  NamedFunc mggwindow_broad = "photonphoton_m >= 110 && photonphoton_m <= 140";
  NamedFunc mggloosencut = "photonphoton_m >= 80 && photonphoton_m <= 200";
  NamedFunc mggloosencut_2 = "photonphoton_m >= 110 && photonphoton_m <= 200";
  NamedFunc mbbloosencut = "50 < bb_m < 150";
  NamedFunc bin1sr = bin1 && novetolepton && Photon_pt_thresh && flavmedcut && njetcut && mbbcut && mggcut && ggdrcut && ggptcut;
  NamedFunc Sept_baseline = bin1 && novetolepton && njetcut && mbbcut && flavmedcut && mggcut && ggdrcut && ggptcut && costhetacut;
  NamedFunc newnjetcut = "njet == 2";
  NamedFunc njetcut_jetflav = "njet >= 2";
  NamedFunc newreverseggdrcut = "photonphoton_dr > 3";
  NamedFunc nphotoncut = "nphoton == 2";

  NamedFunc no_genb_dp = num_gen_b < 1;
  NamedFunc one_genb_dp = num_gen_b >= 1 && num_gen_b < 2;
  NamedFunc two_genb_dp = num_gen_b >= 2 && num_gen_b < 3;
  NamedFunc isonegamma = num_gengamma_basic >= 1 && num_gengamma_basic < 2;
  NamedFunc no_genb = "no_genb";
  NamedFunc one_genb = "one_genb";
  NamedFunc two_genb = "two_genb";
  NamedFunc use_event = "use_event";


  auto proc_mchi_150 = Process::MakeShared<Baby_pico>(SigSampleNames[0], sig, TColor::GetColor("#ff0000"), {mchi_150}, "1");
  auto proc_mchi_500 = Process::MakeShared<Baby_pico>(SigSampleNames[1], sig, TColor::GetColor("#00ff00"), {mchi_500}, "1");
  auto proc_mchi_1000 = Process::MakeShared<Baby_pico>(SigSampleNames[2], sig, TColor::GetColor("#0000ff"), {mchi_1000}, "1");

  auto proc_mchi_150_bkg= Process::MakeShared<Baby_pico>(SigSampleNames[0], bkg, TColor::GetColor("#ff0000"), {mchi_150}, "1");
  auto proc_mchi_500_bkg= Process::MakeShared<Baby_pico>(SigSampleNames[1], bkg, TColor::GetColor("#00ff00"), {mchi_500}, "1");
  auto proc_mchi_1000_bkg = Process::MakeShared<Baby_pico>(SigSampleNames[2], bkg, TColor::GetColor("#0000ff"), {mchi_1000}, "1");

  auto proc_sig150_debug = Process::MakeShared<Baby_pico>("debug file", sig, TColor::GetColor("#006600"), {debugsigfile}, "1");
  auto proc_DiPhoton_mad = Process::MakeShared<Baby_pico>("Diphoton madgraph", bkg, TColor::GetColor("#008080"), {DiPhoton_mad}, "1");

  auto proc_DiPhoton = Process::MakeShared<Baby_pico>("DiPhotonJetsBox", bkg, TColor::GetColor("#006600"), {DiPhoton}, trig && no_genb);
  auto proc_DiPhoton1B = Process::MakeShared<Baby_pico>(BkgSampleNames[1], bkg, TColor::GetColor("#ff6600"), {DiPhoton1B}, trig && one_genb);
  auto proc_DiPhoton2B = Process::MakeShared<Baby_pico>(BkgSampleNames[2], bkg, TColor::GetColor("#ffcc00"), {DiPhoton2B}, trig && two_genb);
  auto proc_GJet1 = Process::MakeShared<Baby_pico>(BkgSampleNames[3], bkg, TColor::GetColor("#99dcff"), {GJet1}, trig && use_event);
  auto proc_GJet2 = Process::MakeShared<Baby_pico>(BkgSampleNames[4], bkg, TColor::GetColor("#9999ff"), {GJet2}, trig && use_event);
  auto proc_TTGG = Process::MakeShared<Baby_pico>(BkgSampleNames[5], bkg, TColor::GetColor("#401641"), {TTGG}, trig);
  auto proc_TTGJets = Process::MakeShared<Baby_pico>(BkgSampleNames[6], bkg, TColor::GetColor("#09ba01"), {TTGJets}, trig);
  auto proc_TTTo2L2Nu = Process::MakeShared<Baby_pico>(BkgSampleNames[7], bkg, TColor::GetColor("#deba87"), {TTTo2L2Nu}, trig);
  auto proc_QCDPt1 = Process::MakeShared<Baby_pico>(BkgSampleNames[8], bkg, TColor::GetColor("#fe6200"), {QCD_Pt1}, trig);
  auto proc_QCDPt2 = Process::MakeShared<Baby_pico>("QCD_Pt-40toInf", bkg, TColor::GetColor("#0194da"), {QCD_Pt2}, trig);

  auto proc_DiPhoton_notrigs = Process::MakeShared<Baby_pico>("DiPhotonJetsBox", bkg, TColor::GetColor("#006600"), {DiPhoton}, "1");
  auto proc_DiPhoton1B_notrigs = Process::MakeShared<Baby_pico>(BkgSampleNames[1], bkg, TColor::GetColor("#ff6600"), {DiPhoton1B}, "1");
  auto proc_DiPhoton2B_notrigs = Process::MakeShared<Baby_pico>(BkgSampleNames[2], bkg, TColor::GetColor("#ffcc00"), {DiPhoton2B}, "1");
  auto proc_GJet1_notrigs = Process::MakeShared<Baby_pico>(BkgSampleNames[3], bkg, TColor::GetColor("#99dcff"), {GJet1}, "1");
  auto proc_GJet2_notrigs = Process::MakeShared<Baby_pico>(BkgSampleNames[4], bkg, TColor::GetColor("#9999ff"), {GJet2}, "1");
  auto proc_TTGG_notrigs = Process::MakeShared<Baby_pico>(BkgSampleNames[5], bkg, TColor::GetColor("#401641"), {TTGG}, "1");
  auto proc_TTGJets_notrigs = Process::MakeShared<Baby_pico>(BkgSampleNames[6], bkg, TColor::GetColor("#09ba01"), {TTGJets}, "1");
  auto proc_TTTo2L2Nu_notrigs = Process::MakeShared<Baby_pico>(BkgSampleNames[7], bkg, TColor::GetColor("#deba87"), {TTTo2L2Nu}, "1");
  auto proc_QCDPt1_notrigs = Process::MakeShared<Baby_pico>(BkgSampleNames[8], bkg, TColor::GetColor("#fe6200"), {QCD_Pt1}, "1");
  auto proc_QCDPt2_notrigs = Process::MakeShared<Baby_pico>("QCD_Pt-40toInf", bkg, TColor::GetColor("#0194da"), {QCD_Pt2}, "1");

  auto proc_DiPhoton_debug = Process::MakeShared<Baby_pico>("debug file", sig, TColor::GetColor("#006600"), {Debug_file_diphoton}, trig);
  auto proc_DiPhoton_sig = Process::MakeShared<Baby_pico>(BkgSampleNames[0], sig, TColor::GetColor("#006600"), {DiPhoton}, trig);
  auto proc_DiPhoton1B_sig = Process::MakeShared<Baby_pico>(BkgSampleNames[1], sig, TColor::GetColor("#ff6600"), {DiPhoton1B}, trig);
  auto proc_DiPhoton2B_sig = Process::MakeShared<Baby_pico>(BkgSampleNames[2], sig, TColor::GetColor("#ffcc00"), {DiPhoton2B}, trig);
  auto proc_GJet1_sig = Process::MakeShared<Baby_pico>(BkgSampleNames[3], sig, TColor::GetColor("#99dcff"), {GJet1}, trig);
  auto proc_GJet2_sig = Process::MakeShared<Baby_pico>(BkgSampleNames[4], sig, TColor::GetColor("#9999ff"), {GJet2}, trig);


  auto proc_data_blind = Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack, {data_16, data_17,data_18}, trig && !mggwindow);
  auto proc_data = Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack, {data_16, data_17,data_18}, trig);
  //auto proc_data = Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack, {data_all}, trig);

  auto proc_DiPhotonAnd1bOverlap = Process::MakeShared<Baby_pico>(BkgSampleNames[0], sig, TColor::GetColor("#006600"), {DiPhoton}, trig);
  auto proc_DiPhotonAnd2bOverlap = Process::MakeShared<Baby_pico>(BkgSampleNames[0], sig, TColor::GetColor("#006600"), {DiPhoton}, trig);
  //vector<shared_ptr<Process>> procs_DiPhoton1bcomp = {proc_DiPhoton1B, proc_DiPhotonAnd1bOverlap};
  //vector<shared_ptr<Process>> procs_DiPhoton2bcomp = {proc_DiPhoton2B, proc_DiPhotonAnd2bOverlap};

  auto proc_DiPhoton1bAnd0bOverlap = Process::MakeShared<Baby_pico>(BkgSampleNames[0], sig, TColor::GetColor("#ff6600"), {DiPhoton1B}, trig);
  auto proc_DiPhoton1bAnd2bOverlap = Process::MakeShared<Baby_pico>(BkgSampleNames[0], sig, TColor::GetColor("#ff6600"), {DiPhoton1B}, trig);
  vector<shared_ptr<Process>> procs_DiPhoton1bcomp = {proc_DiPhoton, proc_DiPhoton1bAnd0bOverlap};
  vector<shared_ptr<Process>> procs_DiPhoton2bcomp = {proc_DiPhoton2B, proc_DiPhoton1bAnd2bOverlap};

  proc_mchi_150->SetLineWidth(3);
  proc_mchi_500->SetLineWidth(3);
  proc_mchi_1000->SetLineWidth(3);

  vector<shared_ptr<Process>> procs_sig_debug = {proc_sig150_debug};
  vector<shared_ptr<Process>> procs_debug_diphoton = {proc_DiPhoton_debug};
  vector<shared_ptr<Process>> procs_sig = {proc_mchi_150, proc_mchi_500, proc_mchi_1000};
  vector<shared_ptr<Process>> procs_bkg = {proc_TTTo2L2Nu, proc_QCDPt2, proc_QCDPt1,  proc_DiPhoton, proc_GJet2, proc_DiPhoton1B, proc_DiPhoton2B, proc_TTGJets, proc_GJet1, proc_TTGG};
  vector<shared_ptr<Process>> procs = {proc_mchi_150, proc_mchi_500, proc_mchi_1000, proc_TTTo2L2Nu, proc_QCDPt2, proc_QCDPt1, proc_DiPhoton, proc_GJet2, proc_DiPhoton1B, proc_DiPhoton2B, proc_TTGJets, proc_GJet1, proc_TTGG};

  vector<shared_ptr<Process>> procs_bkg_notrigs = {proc_TTTo2L2Nu_notrigs, proc_QCDPt2_notrigs, proc_QCDPt1_notrigs,  proc_DiPhoton_notrigs, proc_GJet2_notrigs, proc_DiPhoton1B_notrigs, proc_DiPhoton2B_notrigs, proc_TTGJets_notrigs, proc_GJet1_notrigs, proc_TTGG_notrigs};

  vector<shared_ptr<Process>> procs_sig150_bkg = {proc_mchi_150_bkg};
  vector<shared_ptr<Process>> procs_sig500_bkg = {proc_mchi_500_bkg};
  vector<shared_ptr<Process>> procs_sig1000_bkg = {proc_mchi_1000_bkg};
  vector<shared_ptr<Process>> procs_sig_bkg = {proc_mchi_150_bkg, proc_mchi_500_bkg, proc_mchi_1000_bkg};
  vector<shared_ptr<Process>> procs__CR = {proc_mchi_150_bkg, proc_TTTo2L2Nu, proc_QCDPt2, proc_QCDPt1, proc_DiPhoton, proc_GJet2, proc_DiPhoton1B, proc_DiPhoton2B, proc_TTGJets, proc_GJet1, proc_TTGG};
  vector<shared_ptr<Process>> procs__CR_noqcd = {proc_mchi_150_bkg, proc_TTTo2L2Nu, proc_DiPhoton, proc_GJet2, proc_DiPhoton1B, proc_DiPhoton2B, proc_TTGJets, proc_GJet1, proc_TTGG};
  
  vector<shared_ptr<Process>> procs_SOB_bin1 = {proc_mchi_150_bkg, proc_TTTo2L2Nu, /*proc_QCDPt2, proc_QCDPt1,*/ proc_DiPhoton, proc_GJet2, proc_DiPhoton1B, proc_DiPhoton2B, proc_TTGJets, proc_GJet1, proc_TTGG};
  vector<shared_ptr<Process>> procs_SOB_bin1_sherpa = {proc_mchi_150_bkg, proc_TTTo2L2Nu, proc_QCDPt2, proc_QCDPt1, proc_DiPhoton, proc_GJet2, proc_TTGJets, proc_GJet1, proc_TTGG};
  vector<shared_ptr<Process>> procs_SOB_bin1_mad = {proc_mchi_150_bkg, proc_TTTo2L2Nu, proc_QCDPt2, proc_QCDPt1, proc_DiPhoton_mad, proc_GJet2, proc_TTGJets, proc_GJet1, proc_TTGG};


  vector<shared_ptr<Process>> procs_SOB_bin2 = {proc_mchi_150_bkg, proc_TTTo2L2Nu, proc_QCDPt2, proc_QCDPt1, proc_DiPhoton, proc_GJet2, proc_DiPhoton1B, proc_DiPhoton2B, proc_TTGJets, proc_GJet1, proc_TTGG};
  vector<shared_ptr<Process>> procs_SOB_bin3 = {proc_mchi_500_bkg, proc_TTTo2L2Nu, proc_QCDPt2, proc_QCDPt1, proc_DiPhoton, proc_GJet2, proc_DiPhoton1B, proc_DiPhoton2B, proc_TTGJets, proc_GJet1, proc_TTGG};

  vector<shared_ptr<Process>> procs_DiPhotonsamples_comp = {proc_DiPhoton,proc_DiPhoton1B,proc_DiPhoton2B};
  vector<shared_ptr<Process>> procs_DiPhotonsamples_comp_sig = {proc_DiPhoton_sig,proc_DiPhoton1B_sig,proc_DiPhoton2B_sig};
  vector<shared_ptr<Process>> procs_DiPhoton_sig = {proc_DiPhoton_sig};
  vector<shared_ptr<Process>> procs_DiPhoton1B_sig = {proc_DiPhoton_sig,proc_DiPhoton1B_sig, proc_DiPhoton2B_sig};
  vector<shared_ptr<Process>> procs_DiPhoton2B_sig = {proc_DiPhoton2B_sig};
  vector<shared_ptr<Process>> procs_GJet_sig = {proc_GJet1_sig, proc_GJet2_sig};

  vector<shared_ptr<Process>> procs_DiPhoton_compare = {proc_DiPhoton, proc_DiPhoton_mad};

  vector<shared_ptr<Process>> procs_data = {proc_data};


  vector<shared_ptr<Process>> procs_TTTo2L2N = {proc_TTTo2L2Nu};
  vector<shared_ptr<Process>> procs_DiPhoton = {proc_DiPhoton};
  vector<shared_ptr<Process>> procs_GJet = {proc_GJet1,proc_GJet2};
  vector<shared_ptr<Process>> procs_DiPhoton1B = {proc_DiPhoton1B};
  vector<shared_ptr<Process>> procs_DiPhoton2B = {proc_DiPhoton2B};
  vector<shared_ptr<Process>> procs_TTGJets = {proc_TTGJets};
  vector<shared_ptr<Process>> procs_TTGG = {proc_TTGG};
  vector<shared_ptr<Process>> procs_QCD = {proc_QCDPt2, proc_QCDPt1};

  vector<shared_ptr<Process>> procs_CR_dataMC = {proc_data, proc_TTTo2L2Nu, proc_DiPhoton, proc_GJet2, proc_DiPhoton1B, proc_DiPhoton2B, proc_TTGJets, proc_GJet1, proc_TTGG};
  vector<shared_ptr<Process>> procs_CR_dataMC_withqcd_datablind = {proc_data_blind, proc_TTTo2L2Nu, proc_QCDPt2, proc_QCDPt1, proc_DiPhoton, proc_GJet2, proc_DiPhoton1B, proc_DiPhoton2B, proc_TTGJets, proc_GJet1, proc_TTGG};
  vector<shared_ptr<Process>> procs_CR_dataMC_withqcd = {proc_data, proc_TTTo2L2Nu, proc_QCDPt2, proc_QCDPt1, proc_DiPhoton, proc_GJet2, proc_DiPhoton1B, proc_DiPhoton2B, proc_TTGJets, proc_GJet1, proc_TTGG};

  vector<shared_ptr<Process>> procs_CR_dataMC_withqcd_withsig_datablind = {proc_mchi_150, proc_mchi_500, proc_mchi_1000, proc_data_blind, proc_TTTo2L2Nu, proc_QCDPt2, proc_QCDPt1, proc_DiPhoton, proc_GJet2, proc_DiPhoton1B, proc_DiPhoton2B, proc_TTGJets, proc_GJet1, proc_TTGG};
  vector<shared_ptr<Process>> procs_CR_dataMC_withqcd_withsig = {proc_mchi_150, proc_mchi_500, proc_mchi_1000, proc_data, proc_TTTo2L2Nu, proc_QCDPt2, proc_QCDPt1, proc_DiPhoton, proc_GJet2, proc_DiPhoton1B, proc_DiPhoton2B, proc_TTGJets, proc_GJet1, proc_TTGG};
 
  //vector<shared_ptr<Process>> procs_CR_dataMC_withqcd_withsig_bin2 = {proc_mchi_500, proc_data, proc_TTTo2L2Nu, proc_QCDPt2, proc_QCDPt1, proc_DiPhoton, proc_GJet2, proc_DiPhoton1B, proc_DiPhoton2B, proc_TTGJets, proc_GJet1, proc_TTGG};
  //vector<shared_ptr<Process>> procs_CR_dataMC_withqcd_withsig_bin3 = {proc_mchi_1000, proc_data, proc_TTTo2L2Nu, proc_QCDPt2, proc_QCDPt1, proc_DiPhoton, proc_GJet2, proc_DiPhoton1B, proc_DiPhoton2B, proc_TTGJets, proc_GJet1, proc_TTGG};


  vector<shared_ptr<Process>> procs_CR_dataMC_withqcd_nodata = {proc_TTTo2L2Nu, proc_QCDPt2, proc_QCDPt1, proc_DiPhoton, proc_GJet2, proc_DiPhoton1B, proc_DiPhoton2B, proc_TTGJets, proc_GJet1, proc_TTGG};

  vector<shared_ptr<Process>> procs_CR_dataMC_2 = {proc_data, proc_TTTo2L2Nu, /*proc_QCDPt2, proc_QCDPt1,*/ proc_DiPhoton, proc_GJet2, /*proc_DiPhoton1B, proc_DiPhoton2B,*/ proc_TTGJets, proc_GJet1, proc_TTGG};
  vector<shared_ptr<Process>> procs_CR_dataMC_compareDiph_mad = {proc_data, proc_TTTo2L2Nu, proc_QCDPt2, proc_QCDPt1, proc_DiPhoton_mad, proc_GJet2, proc_TTGJets, proc_GJet1, proc_TTGG};
  vector<shared_ptr<Process>> procs_CR_dataMC_compareDiph_sherpa = {proc_data, proc_TTTo2L2Nu, proc_QCDPt2, proc_QCDPt1, proc_DiPhoton, proc_GJet2, proc_TTGJets, proc_GJet1, proc_TTGG};
  vector<shared_ptr<Process>> procs_CR_dataMC_compareDiph_sherpa_With1b2b = {proc_data, proc_TTTo2L2Nu, proc_QCDPt2, proc_QCDPt1, proc_DiPhoton, proc_GJet2, proc_DiPhoton1B, proc_DiPhoton2B, proc_TTGJets, proc_GJet1, proc_TTGG};

  vector<shared_ptr<Process>> procs_CR_dataMC_compareDiph_noqcd_mad = {proc_data, proc_TTTo2L2Nu, /*proc_QCDPt2, proc_QCDPt1,*/ proc_DiPhoton_mad, proc_GJet2, proc_TTGJets, proc_GJet1, proc_TTGG};
  vector<shared_ptr<Process>> procs_CR_dataMC_compareDiph_noqcd_sherpa = {proc_data, proc_TTTo2L2Nu, /*proc_QCDPt2, proc_QCDPt1,*/ proc_DiPhoton, proc_GJet2, proc_TTGJets, proc_GJet1, proc_TTGG};



  PlotMaker pm;


  ///////////////////////////////
  ////// DiPhoton CR 1 (CR 1)//////////
  ///////////////////////////////

  //gg distributions
  pm.Push<Hist1D>(Axis(24, 80, 200, "photonphoton_m", "M_{\\gamma\\gamma} [GeV]", {115,135}), bin1 && novetolepton && Photon_pt_thresh && !njetcut && !mbbloosewindow && flavloosecut && !mggwindow && ggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:mgg_DiPhotonCR1");
  pm.Push<Hist1D>(Axis(30, 0, 6, "photonphoton_dr", "\\DeltaR_{\\gamma\\gamma}", {2.5}), bin1 && novetolepton && Photon_pt_thresh && !njetcut && !mbbloosewindow && flavloosecut && mggloosencut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:ggdr_DiPhotonCR1");
  pm.Push<Hist1D>(Axis(30, 0, 300, "photonphoton_pt", "p_{T}^{\\gamma\\gamma} [GeV]", {75}), bin1 && novetolepton && Photon_pt_thresh && !njetcut && !mbbloosewindow && flavloosecut && mggloosencut && ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:ggpt_DiPhotonCR1");

  //bb distributions 
  pm.Push<Hist1D>(Axis(40, 0, 300, "bb_m", "M_{jj} [GeV]", {70,140}), bin1 && novetolepton && Photon_pt_thresh && !njetcut && !mbbloosewindow && flavloosecut && mggloosencut && ggptcut && ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:mbb_DiPhotonCR1");
  pm.Push<Hist1D>(Axis(30, 0, 6, "bb_dr", "\\DeltaR_{jj}", {}), bin1 && novetolepton && Photon_pt_thresh && !njetcut && !mbbloosewindow && flavloosecut && mggloosencut && ggptcut && ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:bbdr_DiPhotonCR1");
  pm.Push<Hist1D>(Axis(30, 0, 300, "bb_pt", "p_{T}^{jj} [GeV]", {}), bin1 && novetolepton && Photon_pt_thresh && !njetcut && !mbbloosewindow && flavloosecut && mggloosencut && ggptcut && ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:bbpt_DiPhotonCR1");
  
  //Lead photon distributions
  pm.Push<Hist1D>(Axis(30, 0, 300, Photon_pt_lead, "p_{T}^{lead \\gamma} [GeV]", {35}), bin1 && novetolepton && Photon_pt_thresh && !njetcut && !mbbloosewindow && flavloosecut && mggloosencut && ggptcut && ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadphotonpt_DiPhotonCR1");
  pm.Push<Hist1D>(Axis(20, -1, 1, Photon_id_lead, "lead photon id", {}), bin1 && novetolepton && Photon_pt_thresh && !njetcut && !mbbloosewindow && flavloosecut && mggloosencut && ggptcut && ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadphotonid_DiPhotonCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, Photon_eta_lead, "\\eta^{lead \\gamma}", {}), bin1 && novetolepton && Photon_pt_thresh && !njetcut && !mbbloosewindow && flavloosecut && mggloosencut && ggptcut && ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadphotoneta_DiPhotonCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, Photon_phi_lead, "\\phi^{lead \\gamma}", {}), bin1 && novetolepton && Photon_pt_thresh && !njetcut && !mbbloosewindow && flavloosecut && mggloosencut && ggptcut && ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadphotonphi_DiPhotonCR1");
  
  //Sublead photon distributions
  pm.Push<Hist1D>(Axis(30, 0, 300, Photon_pt_sublead, "p_{T}^{sublead \\gamma} [GeV]", {25}), bin1 && novetolepton && Photon_pt_thresh && !njetcut && !mbbloosewindow && flavloosecut && mggloosencut && ggptcut && ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadphotonpt_DiPhotonCR1");
  pm.Push<Hist1D>(Axis(20, -1, 1, Photon_id_sublead, "sublead photon id", {}), bin1 && novetolepton && Photon_pt_thresh && !njetcut && !mbbloosewindow && flavloosecut && mggloosencut && ggptcut && ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadphotonid_DiPhotonCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, Photon_eta_sublead, "\\eta^{sublead \\gamma}", {}), bin1 && novetolepton && Photon_pt_thresh && !njetcut && !mbbloosewindow && flavloosecut && mggloosencut && ggptcut && ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadphotoneta_DiPhotonCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, Photon_phi_sublead, "\\phi^{sublead \\gamma}", {}), bin1 && novetolepton && Photon_pt_thresh && !njetcut && !mbbloosewindow && flavloosecut && mggloosencut && ggptcut && ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadphotonphi_DiPhotonCR1");
  
  //Lead jet distributions
  pm.Push<Hist1D>(Axis(30, 0, 300, jet_pt_lead, "p_{T}^{lead jet} [GeV]", {}), bin1 && novetolepton && Photon_pt_thresh && !njetcut && !mbbloosewindow && flavloosecut && mggloosencut && ggptcut && ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadjetpt_DiPhotonCR1");
  pm.Push<Hist1D>(Axis(20, 0, 1, DeepFlavScore_lead, "Lead jet DeepFlav", {}), bin1 && novetolepton && Photon_pt_thresh && !njetcut && !mbbloosewindow && mggloosencut && ggptcut && ggdrcut, procs_CR_dataMC_withqcd, ops_log_dnorm).Weight(w_lumi).Tag("ShortName:leadjetid_DiPhotonCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, jet_eta_lead, "\\eta^{lead jet}", {}), bin1 && novetolepton && Photon_pt_thresh && !njetcut && !mbbloosewindow && flavloosecut && mggloosencut && ggptcut && ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadjeteta_DiPhotonCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, jet_phi_lead, "\\phi^{lead jet}", {}), bin1 && novetolepton && Photon_pt_thresh && !njetcut && !mbbloosewindow && flavloosecut && mggloosencut && ggptcut && ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadjetphi_DiPhotonCR1");
  
  //Sublead jet distributions
  pm.Push<Hist1D>(Axis(30, 0, 300, jet_pt_sublead, "p_{T}^{sublead jet} [GeV]", {}), bin1 && novetolepton && Photon_pt_thresh && !njetcut && !mbbloosewindow && flavloosecut && mggloosencut && ggptcut && ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadjetpt_DiPhotonCR1");
  pm.Push<Hist1D>(Axis(20, 0, 1, DeepFlavScore_sublead, "Sublead jet DeepFlav", {}), bin1 && novetolepton && Photon_pt_thresh && !njetcut && !mbbloosewindow && mggloosencut && ggptcut && ggdrcut, procs_CR_dataMC_withqcd, ops_log_dnorm).Weight(w_lumi).Tag("ShortName:subleadjetid_DiPhotonCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, jet_eta_sublead, "\\eta^{sublead jet}", {}), bin1 && novetolepton && Photon_pt_thresh && !njetcut && !mbbloosewindow && flavloosecut && mggloosencut && ggptcut && ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadjeteta_DiPhotonCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, jet_phi_sublead, "\\phi^{sublead jet}", {}), bin1 && novetolepton && Photon_pt_thresh && !njetcut && !mbbloosewindow && flavloosecut && mggloosencut && ggptcut && ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadjetphi_DiPhotonCR1");

  //other distributions
  pm.Push<Hist1D>(Axis(30, 0, 200, "met", "p_{T}^{miss} [GeV]", {20,70}), novetolepton && Photon_pt_thresh && !njetcut && !mbbloosewindow && flavloosecut && mggloosencut && ggptcut && ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:met_DiPhotonCR1");
  pm.Push<Hist1D>(Axis(8, 0, 8, "njet", "Number of jets", {4}), bin1 && novetolepton && Photon_pt_thresh && !mbbloosewindow && flavloosecut && mggloosencut && ggptcut && ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:njet_DiPhotonCR1");
  

  ///////////////////////////////////
  ////// DiPhoton1b2b CR 1 (CR2) //////////
  ///////////////////////////////////

  //gg distributions
  pm.Push<Hist1D>(Axis(24, 80, 200, "photonphoton_m", "M_{\\gamma\\gamma} [GeV]", {115,135}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && !mggwindow && flavtightcut && !ggdrcut && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:mgg_DiPhoton1b2bCR1");
  pm.Push<Hist1D>(Axis(30, 0, 6, "photonphoton_dr", "\\DeltaR_{\\gamma\\gamma}", {2.5}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:ggdr_DiPhoton1b2bCR1");
  pm.Push<Hist1D>(Axis(30, 0, 300, "photonphoton_pt", "p_{T}^{\\gamma\\gamma} [GeV]", {75}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavloosecut && mggloosencut && !ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:ggpt_DiPhoton1b2bCR1");

  //bb distributions 
  pm.Push<Hist1D>(Axis(40, 0, 300, "bb_m", "M_{jj} [GeV]", {70,140}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && !ggdrcut && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:mbb_DiPhoton1b2bCR1");
  pm.Push<Hist1D>(Axis(30, 0, 6, "bb_dr", "\\DeltaR_{jj}", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && !ggdrcut && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:bbdr_DiPhoton1b2bCR1");
  pm.Push<Hist1D>(Axis(30, 0, 300, "bb_pt", "p_{T}^{jj} [GeV]", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && !ggdrcut && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:bbpt_DiPhoton1b2bCR1");
  
  //Lead photon distributions
  pm.Push<Hist1D>(Axis(30, 0, 300, Photon_pt_lead, "p_{T}^{lead \\gamma} [GeV]", {35}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && !ggdrcut && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadphotonpt_DiPhoton1b2bCR1");
  pm.Push<Hist1D>(Axis(20, -1, 1, Photon_id_lead, "lead photon id", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && !ggdrcut && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadphotonid_DiPhoton1b2bCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, Photon_eta_lead, "\\eta^{lead \\gamma}", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && !ggdrcut && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadphotoneta_DiPhoton1b2bCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, Photon_phi_lead, "\\phi^{lead \\gamma}", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && !ggdrcut && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadphotonphi_DiPhoton1b2bCR1");
  
  //Sublead photon distributions
  pm.Push<Hist1D>(Axis(30, 0, 300, Photon_pt_sublead, "p_{T}^{sublead \\gamma} [GeV]", {25}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && !ggdrcut && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadphotonpt_DiPhoton1b2bCR1");
  pm.Push<Hist1D>(Axis(20, -1, 1, Photon_id_sublead, "sublead photon id", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && !ggdrcut && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadphotonid_DiPhoton1b2bCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, Photon_eta_sublead, "\\eta^{sublead \\gamma}", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && !ggdrcut && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadphotoneta_DiPhoton1b2bCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, Photon_phi_sublead, "\\phi^{sublead \\gamma}", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && !ggdrcut && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadphotonphi_DiPhoton1b2bCR1");

  //Lead jet distributions
  pm.Push<Hist1D>(Axis(30, 0, 300, jet_pt_lead, "p_{T}^{lead jet} [GeV]", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && !ggdrcut && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadjetpt_DiPhoton1b2bCR1");
  pm.Push<Hist1D>(Axis(20, 0, 1, DeepFlavScore_lead, "Lead jet DeepFlav", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && !ggdrcut && !ggptcut, procs_CR_dataMC_withqcd, ops_log_dnorm).Weight(w_lumi).Tag("ShortName:leadjetid_DiPhoton1b2bCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, jet_eta_lead, "\\eta^{lead jet}", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && !ggdrcut && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadjeteta_DiPhoton1b2bCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, jet_phi_lead, "\\phi^{lead jet}", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && !ggdrcut && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadjetphi_DiPhoton1b2bCR1");
  
  //Sublead jet distributions
  pm.Push<Hist1D>(Axis(30, 0, 300, jet_pt_sublead, "p_{T}^{sublead jet} [GeV]", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && !ggdrcut && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadjetpt_DiPhoton1b2bCR1");
  pm.Push<Hist1D>(Axis(20, 0, 1, DeepFlavScore_sublead, "Sublead jet DeepFlav", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && !ggdrcut && !ggptcut, procs_CR_dataMC_withqcd, ops_log_dnorm).Weight(w_lumi).Tag("ShortName:subleadjetid_DiPhoton1b2bCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, jet_eta_sublead, "\\eta^{sublead jet}", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && !ggdrcut && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadjeteta_DiPhoton1b2bCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, jet_phi_sublead, "\\phi^{sublead jet}", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && !ggdrcut && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadjetphi_DiPhoton1b2bCR1");

  //other distributions
  pm.Push<Hist1D>(Axis(30, 0, 200, "met", "p_{T}^{miss} [GeV]", {20,70}), novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && !ggdrcut && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:met_DiPhoton1b2bCR1");
  pm.Push<Hist1D>(Axis(8, 0, 8, "njet", "Number of jets", {4}), bin1 && novetolepton && Photon_pt_thresh && !mbbloosewindow && mggloosencut && flavtightcut && !ggdrcut && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:njet_DiPhoton1b2bCR1");


  ///////////////////////////////////
  ////// DiPhoton1b2b CR 2 (CR 5) //////////
  ///////////////////////////////////

  //gg distributions
  pm.Push<Hist1D>(Axis(24, 80, 200, "photonphoton_m", "M_{\\gamma\\gamma} [GeV]", {115,135}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && !mggwindow && flavtightcut && Photon_id80 && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:mgg_DiPhoton1b2bCR2");
  pm.Push<Hist1D>(Axis(30, 0, 6, "photonphoton_dr", "\\DeltaR_{\\gamma\\gamma}", {2.5}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && Photon_id80 && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:ggdr_DiPhoton1b2bCR2");
  pm.Push<Hist1D>(Axis(30, 0, 300, "photonphoton_pt", "p_{T}^{\\gamma\\gamma} [GeV]", {75}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavloosecut && mggloosencut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:ggpt_DiPhoton1b2bCR2");

  //bb distributions 
  pm.Push<Hist1D>(Axis(40, 0, 300, "bb_m", "M_{jj} [GeV]", {70,140}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && Photon_id80 && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:mbb_DiPhoton1b2bCR2");
  pm.Push<Hist1D>(Axis(30, 0, 6, "bb_dr", "\\DeltaR_{jj}", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && Photon_id80 && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:bbdr_DiPhoton1b2bCR2");
  pm.Push<Hist1D>(Axis(30, 0, 300, "bb_pt", "p_{T}^{jj} [GeV]", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && Photon_id80 && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:bbpt_DiPhoton1b2bCR2");
  
  //Lead photon distributions
  pm.Push<Hist1D>(Axis(30, 0, 300, Photon_pt_lead, "p_{T}^{lead \\gamma} [GeV]", {35}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && Photon_id80 && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadphotonpt_DiPhoton1b2bCR2");
  pm.Push<Hist1D>(Axis(20, -1, 1, Photon_id_lead, "lead photon id", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && Photon_id80 && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadphotonid_DiPhoton1b2bCR2");
  pm.Push<Hist1D>(Axis(30, -4, 4, Photon_eta_lead, "\\eta^{lead \\gamma}", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && Photon_id80 && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadphotoneta_DiPhoton1b2bCR2");
  pm.Push<Hist1D>(Axis(30, -4, 4, Photon_phi_lead, "\\phi^{lead \\gamma}", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && Photon_id80 && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadphotonphi_DiPhoton1b2bCR2");
  
  //Sublead photon distributions
  pm.Push<Hist1D>(Axis(30, 0, 300, Photon_pt_sublead, "p_{T}^{sublead \\gamma} [GeV]", {25}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && Photon_id80 && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadphotonpt_DiPhoton1b2bCR2");
  pm.Push<Hist1D>(Axis(20, -1, 1, Photon_id_sublead, "sublead photon id", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && Photon_id80 && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadphotonid_DiPhoton1b2bCR2");
  pm.Push<Hist1D>(Axis(30, -4, 4, Photon_eta_sublead, "\\eta^{sublead \\gamma}", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && Photon_id80 && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadphotoneta_DiPhoton1b2bCR2");
  pm.Push<Hist1D>(Axis(30, -4, 4, Photon_phi_sublead, "\\phi^{sublead \\gamma}", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && Photon_id80 && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadphotonphi_DiPhoton1b2bCR2");

  //Lead jet distributions
  pm.Push<Hist1D>(Axis(30, 0, 300, jet_pt_lead, "p_{T}^{lead jet} [GeV]", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && Photon_id80 && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadjetpt_DiPhoton1b2bCR2");
  pm.Push<Hist1D>(Axis(20, 0, 1, DeepFlavScore_lead, "Lead jet DeepFlav", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && Photon_id80 && !ggptcut, procs_CR_dataMC_withqcd, ops_log_dnorm).Weight(w_lumi).Tag("ShortName:leadjetid_DiPhoton1b2bCR2");
  pm.Push<Hist1D>(Axis(30, -4, 4, jet_eta_lead, "\\eta^{lead jet}", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && Photon_id80 && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadjeteta_DiPhoton1b2bCR2");
  pm.Push<Hist1D>(Axis(30, -4, 4, jet_phi_lead, "\\phi^{lead jet}", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && Photon_id80 && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadjetphi_DiPhoton1b2bCR2");
  
  //Sublead jet distributions
  pm.Push<Hist1D>(Axis(30, 0, 300, jet_pt_sublead, "p_{T}^{sublead jet} [GeV]", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && Photon_id80 && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadjetpt_DiPhoton1b2bCR2");
  pm.Push<Hist1D>(Axis(20, 0, 1, DeepFlavScore_sublead, "Sublead jet DeepFlav", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && Photon_id80 && !ggptcut, procs_CR_dataMC_withqcd, ops_log_dnorm).Weight(w_lumi).Tag("ShortName:subleadjetid_DiPhoton1b2bCR2");
  pm.Push<Hist1D>(Axis(30, -4, 4, jet_eta_sublead, "\\eta^{sublead jet}", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && Photon_id80 && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadjeteta_DiPhoton1b2bCR2");
  pm.Push<Hist1D>(Axis(30, -4, 4, jet_phi_sublead, "\\phi^{sublead jet}", {}), bin1 && novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && Photon_id80 && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadjetphi_DiPhoton1b2bCR2");

  //other distributions
  pm.Push<Hist1D>(Axis(30, 0, 200, "met", "p_{T}^{miss} [GeV]", {20,70}), novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && mggloosencut && flavtightcut && Photon_id80 && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:met_DiPhoton1b2bCR2");
  pm.Push<Hist1D>(Axis(8, 0, 8, "njet", "Number of jets", {4}), bin1 && novetolepton && Photon_pt_thresh && !mbbloosewindow && mggloosencut && flavtightcut && Photon_id80 && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:njet_DiPhoton1b2bCR2");
  

  ///////////////////////////////////
  ////// GJet CR 1 (CR 3)//////////////////
  ///////////////////////////////////

  //gg distributions
  pm.Push<Hist1D>(Axis(24, 80, 200, "photonphoton_m", "M_{\\gamma\\gamma} [GeV]", {115,135}), bin1 && novetolepton && Photon_pt_thresh && newnjetcut && mbbloosencut && flavloosecut && !mggwindow && newreverseggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:mgg_GJetCR1");
  pm.Push<Hist1D>(Axis(30, 0, 6, "photonphoton_dr", "\\DeltaR_{\\gamma\\gamma}", {3}), bin1 && novetolepton && Photon_pt_thresh && newnjetcut && mbbloosencut && flavloosecut && mggloosencut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:ggdr_GJetCR1");
  pm.Push<Hist1D>(Axis(30, 0, 300, "photonphoton_pt", "p_{T}^{\\gamma\\gamma} [GeV]", {75}), bin1 && novetolepton && Photon_pt_thresh && newnjetcut && mbbloosencut && flavloosecut && mggloosencut && newreverseggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:ggpt_GJetCR1");

  //bb distributions 
  pm.Push<Hist1D>(Axis(40, 0, 300, "bb_m", "M_{jj} [GeV]", {70,140}), bin1 && novetolepton && Photon_pt_thresh && newnjetcut && !mbbloosewindow && flavloosecut && mggloosencut && newreverseggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:mbb_GJetCR1");
  pm.Push<Hist1D>(Axis(30, 0, 6, "bb_dr", "\\DeltaR_{jj}", {}), bin1 && novetolepton && Photon_pt_thresh && newnjetcut && mbbloosencut && flavloosecut && mggloosencut && newreverseggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:bbdr_GJetCR1");
  pm.Push<Hist1D>(Axis(30, 0, 300, "bb_pt", "p_{T}^{jj} [GeV]", {}), bin1 && novetolepton && Photon_pt_thresh && newnjetcut && mbbloosencut && flavloosecut && mggloosencut && newreverseggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:bbpt_GJetCR1");
  
  //Lead photon distributions
  pm.Push<Hist1D>(Axis(30, 0, 300, Photon_pt_lead, "p_{T}^{lead \\gamma} [GeV]", {35}), bin1 && novetolepton && Photon_pt_thresh && newnjetcut && mbbloosencut && flavloosecut && mggloosencut && newreverseggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadphotonpt_GJetCR1");
  pm.Push<Hist1D>(Axis(20, -1, 1, Photon_id_lead, "lead photon id", {}), bin1 && novetolepton && Photon_pt_thresh && newnjetcut && mbbloosencut && flavloosecut && mggloosencut && newreverseggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadphotonid_GJetCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, Photon_eta_lead, "\\eta^{lead \\gamma}", {}), bin1 && novetolepton && Photon_pt_thresh && newnjetcut && mbbloosencut && flavloosecut && mggloosencut && newreverseggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadphotoneta_GJetCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, Photon_phi_lead, "\\phi^{lead \\gamma}", {}), bin1 && novetolepton && Photon_pt_thresh && newnjetcut && mbbloosencut && flavloosecut && mggloosencut && newreverseggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadphotonphi_GJetCR1");
  
  //Sublead photon distributions
  pm.Push<Hist1D>(Axis(30, 0, 300, Photon_pt_sublead, "p_{T}^{sublead \\gamma} [GeV]", {25}), bin1 && novetolepton && Photon_pt_thresh && newnjetcut && mbbloosencut && flavloosecut && mggloosencut && newreverseggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadphotonpt_GJetCR1");
  pm.Push<Hist1D>(Axis(20, -1, 1, Photon_id_sublead, "sublead photon id", {}), bin1 && novetolepton && Photon_pt_thresh && newnjetcut && mbbloosencut && flavloosecut && mggloosencut && newreverseggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadphotonid_GJetCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, Photon_eta_sublead, "\\eta^{sublead \\gamma}", {}), bin1 && novetolepton && Photon_pt_thresh && newnjetcut && mbbloosencut && flavloosecut && mggloosencut && newreverseggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadphotoneta_GJetCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, Photon_phi_sublead, "\\phi^{sublead \\gamma}", {}), bin1 && novetolepton && Photon_pt_thresh && newnjetcut && mbbloosencut && flavloosecut && mggloosencut && newreverseggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadphotonphi_GJetCR1");

  //Lead jet distributions
  pm.Push<Hist1D>(Axis(30, 0, 300, jet_pt_lead, "p_{T}^{lead jet} [GeV]", {}), bin1 && novetolepton && Photon_pt_thresh && newnjetcut && mbbloosencut && flavloosecut && mggloosencut && newreverseggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadjetpt_GJetCR1");
  pm.Push<Hist1D>(Axis(20, 0, 1, DeepFlavScore_lead, "Lead jet DeepFlav", {}), bin1 && novetolepton && Photon_pt_thresh && newnjetcut && mbbloosencut && mggloosencut && newreverseggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_log_dnorm).Weight(w_lumi).Tag("ShortName:leadjetid_GJetCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, jet_eta_lead, "\\eta^{lead jet}", {}), bin1 && novetolepton && Photon_pt_thresh && newnjetcut && mbbloosencut && flavloosecut && mggloosencut && newreverseggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadjeteta_GJetCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, jet_phi_lead, "\\phi^{lead jet}", {}), bin1 && novetolepton && Photon_pt_thresh && newnjetcut && mbbloosencut && flavloosecut && mggloosencut && newreverseggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadjetphi_GJetCR1");
  
  //Sublead jet distributions
  pm.Push<Hist1D>(Axis(30, 0, 300, jet_pt_sublead, "p_{T}^{sublead jet} [GeV]", {}), bin1 && novetolepton && Photon_pt_thresh && newnjetcut && mbbloosencut && flavloosecut && mggloosencut && newreverseggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadjetpt_GJetCR1");
  pm.Push<Hist1D>(Axis(20, 0, 1, DeepFlavScore_sublead, "Sublead jet DeepFlav", {}), bin1 && novetolepton && Photon_pt_thresh && newnjetcut && mbbloosencut && mggloosencut && newreverseggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_log_dnorm).Weight(w_lumi).Tag("ShortName:subleadjetid_GJetCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, jet_eta_sublead, "\\eta^{sublead jet}", {}), bin1 && novetolepton && Photon_pt_thresh && newnjetcut && mbbloosencut && flavloosecut && mggloosencut && newreverseggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadjeteta_GJetCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, jet_phi_sublead, "\\phi^{sublead jet}", {}), bin1 && novetolepton && Photon_pt_thresh && newnjetcut && mbbloosencut && flavloosecut && mggloosencut && newreverseggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadjetphi_GJetCR1");

  //other distributions
  pm.Push<Hist1D>(Axis(30, 0, 200, "met", "p_{T}^{miss} [GeV]", {20,70}), novetolepton && Photon_pt_thresh && newnjetcut && mbbloosencut && flavloosecut && mggloosencut && newreverseggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:met_GJetCR1");
  pm.Push<Hist1D>(Axis(8, 0, 8, "njet", "Number of jets", {3}), bin1 && novetolepton && Photon_pt_thresh && mbbloosencut && flavloosecut && mggloosencut && newreverseggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:njet_GJetCR1");
  

  ///////////////////////////////////
  ////// GJet CR 2 (CR 6)//////////////////
  ///////////////////////////////////

  //gg distributions
  pm.Push<Hist1D>(Axis(24, 80, 200, "photonphoton_m", "M_{\\gamma\\gamma} [GeV]", {115,135}), bin1 && novetolepton && Photon_pt_thresh && newnjetcut && Photon_id80_lead && DeepFlavBool_tightloose && !mggwindow && !ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:mgg_GJetCR2");
  pm.Push<Hist1D>(Axis(30, 0, 300, "photonphoton_pt", "p_{T}^{\\gamma\\gamma} [GeV]", {75}), bin1 && novetolepton && Photon_pt_thresh && newnjetcut && Photon_id80_lead && DeepFlavBool_tightloose && mggloosencut && !ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:ggpt_binsws_GJetCR2");
  
  //bb distributions 
  pm.Push<Hist1D>(Axis(40, 0, 300, "bb_m", "M_{jj} [GeV]", {70,140}), bin1 && novetolepton && Photon_pt_thresh && newnjetcut && Photon_id80_lead && DeepFlavBool_tightloose && mggloosencut && !mbbloosewindow && !ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:mbb_GJetCR2");
  pm.Push<Hist1D>(Axis(30, 0, 6, "bb_dr", "\\DeltaR_{jj}", {}), bin1 && novetolepton && Photon_pt_thresh && newnjetcut && Photon_id80_lead && DeepFlavBool_tightloose && mggloosencut && !ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:bbdr_GJetCR2");
  pm.Push<Hist1D>(Axis(30, 0, 300, "bb_pt", "p_{T}^{jj} [GeV]", {}), bin1 && novetolepton && Photon_pt_thresh && newnjetcut && Photon_id80_lead && DeepFlavBool_tightloose && mggloosencut && !ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:bbpt_GJetCR2");
  
  //other distributions
  pm.Push<Hist1D>(Axis(30, 0, 200, "met", "p_{T}^{miss} [GeV]", {20,70}), novetolepton && Photon_pt_thresh && newnjetcut && Photon_id80_lead && DeepFlavBool_tightloose && mggloosencut && !ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:met_GJetCR2");
  pm.Push<Hist1D>(Axis(8, 0, 8, "njet", "Number of jets", {3}), bin1 && novetolepton && Photon_pt_thresh && Photon_id80_lead && DeepFlavBool_tightloose && mggloosencut && !ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:njet_GJetCR2");


  ///////////////////////////////////
  ////// QCD CR 1 (CR 7)///////////////////
  ///////////////////////////////////

  //gg distributions
  pm.Push<Hist1D>(Axis(24, 80, 200, "photonphoton_m", "M_{\\gamma\\gamma} [GeV]", {115,135}), bin1 && novetolepton && Photon_pt_thresh && !mggwindow && !ggdrcut && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:mgg_QCDCR1");
  pm.Push<Hist1D>(Axis(30, 0, 6, "photonphoton_dr", "\\DeltaR_{\\gamma\\gamma}", {2.5}), bin1 && novetolepton && Photon_pt_thresh && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:ggdr_QCDCR1");
  pm.Push<Hist1D>(Axis(30, 0, 300, "photonphoton_pt", "p_{T}^{\\gamma\\gamma} [GeV]", {75}), bin1 && novetolepton && Photon_pt_thresh && !ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:ggpt_binsws_QCDCR1");
  
  //bb distributions 
  pm.Push<Hist1D>(Axis(40, 0, 300, "bb_m", "M_{jj} [GeV]", {70,140}), bin1 && novetolepton && Photon_pt_thresh && !mbbloosewindow && !ggdrcut && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:mbb_QCDCR1");
  pm.Push<Hist1D>(Axis(30, 0, 6, "bb_dr", "\\DeltaR_{jj}", {}), bin1 && novetolepton && Photon_pt_thresh && !ggdrcut && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:bbdr_QCDCR1");
  pm.Push<Hist1D>(Axis(30, 0, 300, "bb_pt", "p_{T}^{jj} [GeV]", {}), bin1 && novetolepton && Photon_pt_thresh && !ggdrcut && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:bbpt_QCDCR1");
  
  //other distributions
  pm.Push<Hist1D>(Axis(30, 0, 200, "met", "p_{T}^{miss} [GeV]", {20,70}), novetolepton && Photon_pt_thresh && !ggdrcut && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:met_QCDCR1");
  pm.Push<Hist1D>(Axis(8, 0, 8, "njet", "Number of jets", {4}), bin1 && novetolepton && Photon_pt_thresh && !ggdrcut && !ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:njet_QCDCR1"); 


  ///////////////////////////////////
  ////// TTT CR 1 (CR 4) ///////////////////
  ///////////////////////////////////

  //gg distributions
  pm.Push<Hist1D>(Axis(24, 80, 200, "photonphoton_m", "M_{\\gamma\\gamma} [GeV]", {115,135}), bin1 && !novetolepton && Photon_pt_thresh && njetcut && mbbloosencut && flavmedcut && !mggwindow && !ggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:mgg_TTTCR1");
  pm.Push<Hist1D>(Axis(30, 0, 6, "photonphoton_dr", "\\DeltaR_{\\gamma\\gamma}", {2.5}), bin1 && !novetolepton && Photon_pt_thresh && njetcut && mbbloosencut && flavmedcut && mggloosencut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:ggdr_TTTCR1");
  pm.Push<Hist1D>(Axis(30, 0, 300, "photonphoton_pt", "p_{T}^{\\gamma\\gamma} [GeV]", {75}), bin1 && !novetolepton && Photon_pt_thresh && njetcut && mbbloosencut && flavmedcut && mggloosencut && !ggdrcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:ggpt_TTTCR1");

  //bb distributions 
  pm.Push<Hist1D>(Axis(40, 0, 300, "bb_m", "M_{jj} [GeV]", {70,140}), bin1 && !novetolepton && Photon_pt_thresh && njetcut && !mbbloosewindow && flavmedcut && mggloosencut && !ggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:mbb_TTTCR1");
  pm.Push<Hist1D>(Axis(30, 0, 6, "bb_dr", "\\DeltaR_{jj}", {}), bin1 && !novetolepton && Photon_pt_thresh && njetcut && mbbloosencut && flavmedcut && mggloosencut && !ggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:bbdr_TTTCR1");
  pm.Push<Hist1D>(Axis(30, 0, 300, "bb_pt", "p_{T}^{jj} [GeV]", {}), bin1 && !novetolepton && Photon_pt_thresh && njetcut && mbbloosencut && flavmedcut && mggloosencut && !ggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:bbpt_TTTCR1");
  
  //Lead photon distributions
  pm.Push<Hist1D>(Axis(30, 0, 300, Photon_pt_lead, "p_{T}^{lead \\gamma} [GeV]", {35}), bin1 && !novetolepton && Photon_pt_thresh && njetcut && mbbloosencut && flavmedcut && mggloosencut && !ggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadphotonpt_TTTCR1");
  pm.Push<Hist1D>(Axis(20, -1, 1, Photon_id_lead, "lead photon id", {}), bin1 && !novetolepton && Photon_pt_thresh && njetcut && mbbloosencut && flavmedcut && mggloosencut && !ggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadphotonid_TTTCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, Photon_eta_lead, "\\eta^{lead \\gamma}", {}), bin1 && !novetolepton && Photon_pt_thresh && njetcut && mbbloosencut && flavmedcut && mggloosencut && !ggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadphotoneta_TTTCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, Photon_phi_lead, "\\phi^{lead \\gamma}", {}), bin1 && !novetolepton && Photon_pt_thresh && njetcut && mbbloosencut && flavmedcut && mggloosencut && !ggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadphotonphi_TTTCR1");
  
  //Sublead photon distributions
  pm.Push<Hist1D>(Axis(30, 0, 300, Photon_pt_sublead, "p_{T}^{sublead \\gamma} [GeV]", {25}), bin1 && !novetolepton && Photon_pt_thresh && njetcut && mbbloosencut && flavmedcut && mggloosencut && !ggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadphotonpt_TTTCR1");
  pm.Push<Hist1D>(Axis(20, -1, 1, Photon_id_sublead, "sublead photon id", {}), bin1 && !novetolepton && Photon_pt_thresh && njetcut && mbbloosencut && flavmedcut && mggloosencut && !ggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadphotonid_TTTCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, Photon_eta_sublead, "\\eta^{sublead \\gamma}", {}), bin1 && !novetolepton && Photon_pt_thresh && njetcut && mbbloosencut && flavmedcut && mggloosencut && !ggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadphotoneta_TTTCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, Photon_phi_sublead, "\\phi^{sublead \\gamma}", {}), bin1 && !novetolepton && Photon_pt_thresh && njetcut && mbbloosencut && flavmedcut && mggloosencut && !ggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadphotonphi_TTTCR1");

  //Lead jet distributions
  pm.Push<Hist1D>(Axis(30, 0, 300, jet_pt_lead, "p_{T}^{lead jet} [GeV]", {}), bin1 && !novetolepton && Photon_pt_thresh && njetcut && mbbloosencut && flavmedcut && mggloosencut && !ggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadjetpt_TTTCR1");
  pm.Push<Hist1D>(Axis(20, 0, 1, DeepFlavScore_lead, "Lead jet DeepFlav", {}), bin1 && !novetolepton && Photon_pt_thresh && njetcut && mbbloosencut && mggloosencut && !ggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_log_dnorm).Weight(w_lumi).Tag("ShortName:leadjetid_TTTCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, jet_eta_lead, "\\eta^{lead jet}", {}), bin1 && !novetolepton && Photon_pt_thresh && njetcut && mbbloosencut && flavmedcut && mggloosencut && !ggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadjeteta_TTTCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, jet_phi_lead, "\\phi^{lead jet}", {}), bin1 && !novetolepton && Photon_pt_thresh && njetcut && mbbloosencut && flavmedcut && mggloosencut && !ggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:leadjetphi_TTTCR1");
  
  //Sublead jet distributions
  pm.Push<Hist1D>(Axis(30, 0, 300, jet_pt_sublead, "p_{T}^{sublead jet} [GeV]", {}), bin1 && !novetolepton && Photon_pt_thresh && njetcut && mbbloosencut && flavmedcut && mggloosencut && !ggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadjetpt_TTTCR1");
  pm.Push<Hist1D>(Axis(20, 0, 1, DeepFlavScore_sublead, "Sublead jet DeepFlav", {}), bin1 && !novetolepton && Photon_pt_thresh && njetcut && mbbloosencut && mggloosencut && !ggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_log_dnorm).Weight(w_lumi).Tag("ShortName:subleadjetid_TTTCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, jet_eta_sublead, "\\eta^{sublead jet}", {}), bin1 && !novetolepton && Photon_pt_thresh && njetcut && mbbloosencut && flavmedcut && mggloosencut && !ggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadjeteta_TTTCR1");
  pm.Push<Hist1D>(Axis(30, -4, 4, jet_phi_sublead, "\\phi^{sublead jet}", {}), bin1 && !novetolepton && Photon_pt_thresh && njetcut && mbbloosencut && flavmedcut && mggloosencut && !ggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:subleadjetphi_TTTCR1");

  //other distributions
  pm.Push<Hist1D>(Axis(30, 0, 200, "met", "p_{T}^{miss} [GeV]", {20,70}), !novetolepton && Photon_pt_thresh && njetcut && mbbloosencut && flavmedcut && mggloosencut && !ggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:met_TTTCR1");
  pm.Push<Hist1D>(Axis(8, 0, 8, "njet", "Number of jets", {4}), bin1 && !novetolepton && Photon_pt_thresh && mbbloosencut && flavmedcut && mggloosencut && !ggdrcut && ggptcut, procs_CR_dataMC_withqcd, ops_lin_dnorm).Weight(w_lumi).Tag("ShortName:njet_TTTCR1");


  pm.min_print_ = true;
  pm.MakePlots(Luminosity);
}

