//-----------------------------README-------------------------------//
//This script creates one off cutflow tables. Can also be used for
//printing event details. Extremely inefficient for making more than
//one cutflow at a time. But very flexible. All inputs defined in 
//this file, not from command line.
//
//This file contains a fair number of obsolete functions at the start
//I don't particularly care to clean them right now.

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <getopt.h>
#include <regex>

#include "TError.h" // Controls error level reporting
#include "TColor.h"
#include "TLegend.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TMath.h"

//temp
#include <Math/Rotation3D.h>
#include <Math/RotationX.h>
#include <Math/RotationY.h>
#include <Math/RotationZ.h>
#include <Math/Vector4D.h>
#include <Math/Boost.h>
#include <Math/BoostX.h>
#include <Math/BoostY.h>
#include <Math/BoostZ.h>

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
#include "zgamma/zg_functions.hpp"
#include "zgamma/zg_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace ZgFunctions;
using namespace ZgUtilities;


namespace{
  const string tag_a = "htozgamma_pinnacles_v0";
  const string path_to_production_base = "/net/cms11/cms11r0/pico/";
  const string years_string_a[1]={"2017"};
  const string skim_a = "unskimmed";
  float lumi = 1;
}

bool Contains(const string& text, const string& pattern){
  return text.find(pattern) != string::npos;
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

const NamedFunc w_lumi("w_lumi", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleTypeString().Contains("-")) return 1.;
  else return b.w_lumi();
});

const NamedFunc fullwgt("fullwgt",[](const Baby &b) -> NamedFunc::ScalarType{
  if(b.SampleTypeString().Contains("-")) return 1;
  return b.weight();
 });

const NamedFunc combwgt("combwgt", [](const Baby &b)-> NamedFunc::ScalarType{//customizable weight application
  return b.w_prefire()*b.w_pu()*b.w_bhig_df();
});

const NamedFunc get_Zid("get_Zid",[](const Baby &b)->NamedFunc::ScalarType{//Utilities to select Z, y from scratch
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

const NamedFunc get_yid("get_yid",[](const Baby &b)->NamedFunc::ScalarType{//using this to get highest pT
  int yid = -999;
  double ypt = -10.;
  for (int i = 0; i < b.nphoton(); i++){
    if(b.photon_pt()->at(i)>ypt){//id80 req now in n2p
      yid = i;
      ypt = b.photon_pt()->at(i);
      //break;
    }
  }
  return yid;
});

const NamedFunc sig_g_fix("sig_g_fix",[](const Baby &b)->NamedFunc::ScalarType{
  bool pass = false;
  if (get_yid.GetScalar(b) < 0) return pass;
  else {
    pass = true;
  }
  return pass;
});

const NamedFunc get_Hid("get_Hid",[](const Baby &b)->NamedFunc::ScalarType{
  int Hid = -999;
  if (b.nllphoton() < 1 || get_yid.GetScalar(b) < 0) return Hid;
  else{
    for (int i = 0; i < b.nllphoton(); i++){
      if(b.llphoton_ill()->at(i) == 0  && b.llphoton_iph()->at(i) == get_yid.GetScalar(b)){
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
    return b.ll_m()->at(0) + b.llphoton_m()->at(get_Hid.GetScalar(b));
  }
});

const NamedFunc get_mllyratio("get_mllyratio",[](const Baby &b)->NamedFunc::ScalarType{
  double out = -1;
  if(get_Hid.GetScalar(b) < 0) return out;
  else{
    return b.llphoton_pt()->at(get_Hid.GetScalar(b))/b.llphoton_m()->at(get_Hid.GetScalar(b));
  }
});

const NamedFunc max_mini_iso("max_mini_iso",[](const Baby &b)->NamedFunc::ScalarType{
  double out = -1;
  for (unsigned int i = 0; i<b.el_pt()->size(); i++){
    if(b.el_sig()->at(i)==1 && b.el_miniso()->at(i) > out){
      out=b.el_miniso()->at(i);
    }
  }
  for (unsigned int j = 0; j<b.mu_pt()->size(); j++){
    if(b.mu_sig()->at(j)==1 && b.mu_miniso()->at(j) > out){
      out=b.mu_miniso()->at(j);
    }
  }
  return out;
});

const NamedFunc sideband("sideband",[](const Baby &b)->NamedFunc::ScalarType{//fix to use proper mlly
  bool out = true;
  if(get_mlly.GetScalar(b)>120.f && get_mlly.GetScalar(b)<130.f) return true;//changed for now, false for blinding
  else{
    return out;
  }
});

const NamedFunc sideband2("sideband2",[](const Baby &b)->NamedFunc::ScalarType{//Fix to use proper mlly
  bool out = true;
  //if(b.llphoton_m()->at(0)>122.f && b.llphoton_m()->at(0)<128.f) out=false;
  if(b.llphoton_m()->at(0)>122.f && b.llphoton_m()->at(0)<128.f) out=true;
  return out;
});

const NamedFunc getJets("getJets",[](const Baby &b)->NamedFunc::ScalarType{//gets number of all jets, not just signal
  int out = 0;
  for(unsigned int i = 0;i < b.jet_pt()->size();i++){
    if (b.jet_pt()->at(i)>30) out+=1;
  }
  return out;
});

const NamedFunc getjID("getjID",[](const Baby &b)->NamedFunc::ScalarType{//gets number of all jets, not just signal
  bool out = false;
  for(unsigned int i = 0;i < b.jet_pt()->size();i++){
    if (b.jet_pt()->at(i)>30 && b.jet_id()->at(i)==true) out = true;
  }
  return out;
});

const long double PI = acos(-1.L);

const NamedFunc printer("printer",[](const Baby &b)->NamedFunc::ScalarType{
  ofstream myfile;
  myfile.open ("event_print.txt", std::ios_base::app);
  myfile<<b.run()<<" ";
  myfile<<b.lumiblock()<<" ";
  myfile<<b.event()<<endl;
  myfile.close();
  return true;
});

const NamedFunc printer2("printer2",[](const Baby &b)->NamedFunc::ScalarType{
  cout<<b.event()<<endl;
  for (unsigned int i(0); i<b.mu_pt()->size();i++){
    cout<<"Muon "<<i<<endl;
    cout<<" muon_pt: "<<b.mu_pt()->at(i)<<endl;
    cout<<" muon_eta: "<<b.mu_eta()->at(i)<<endl;
    cout<<" muon_phi: "<<b.mu_phi()->at(i)<<endl;
  }
  for(unsigned int i(0);i<b.ll_m()->size();i++) cout<<"mll: "<<b.ll_m()->at(i)<<endl;
  cout<<"mllg: "<<b.llphoton_m()->at(0)<<endl;
  cout<<"ratio: "<<get_ratio.GetScalar(b)<<endl;
  return true;
});

const NamedFunc cosTheta_daje("cosTheta_daje",[](const Baby &b)->NamedFunc::ScalarType{

  ROOT::Math::PtEtaPhiMVector llg_(b.llphoton_pt()->at(get_Hid.GetScalar(b)), b.llphoton_eta()->at(get_Hid.GetScalar(b)), b.llphoton_phi()->at(get_Hid.GetScalar(b)), b.llphoton_m()->at(get_Hid.GetScalar(b)));//lab frame
  ROOT::Math::PtEtaPhiMVector ll_(b.ll_pt()->at(0), b.ll_eta()->at(0), b.ll_phi()->at(0), b.ll_m()->at(0));//lab frame

  ROOT::Math::PxPyPzEVector llg, ll;
  llg = llg_;
  ll = ll_;
  float llg_theta = llg.Theta();
  float llg_phi = llg.Phi();

  ROOT::Math::RotationZ r_z_llg(-llg_phi);//rotation around z axis
  ROOT::Math::RotationY r_y_llg(-llg_theta);//rotation around y axis
  ROOT::Math::PxPyPzEVector llg_derotated = r_y_llg(r_z_llg(llg));//H in rotated frame with Higgs along Z

  ROOT::Math::PxPyPzEVector ll_derotated = r_y_llg(r_z_llg(ll));//Z in rotated frame with Higgs along Z

  ROOT::Math::BoostZ b_llg(-llg_derotated.P()/llg_derotated.E());//boost factor for going to higgs rest frame
  ROOT::Math::PxPyPzEVector ll_in_llg_rest = b_llg(ll_derotated);

  float theta_c = ll_in_llg_rest.Theta();
  float cos_Theta = cos(theta_c);
  return cos_Theta;
});

float get_dr(float eta1, float phi1, float eta2, float phi2) {
  //const double PI = 3.1415;
  double dphi = fmod(fabs(phi2-phi1), 2.*PI);
  dphi = dphi>PI ? 2.*PI-dphi : dphi;
  double deta = fabs(eta1-eta2);
  return sqrt(deta*deta+dphi*dphi);
} 

float get_max_dr(vector<float> photon_eta, vector<float> photon_phi, 
    vector<float> el_eta, vector<float> el_phi, vector<float> mu_eta,
    vector<float> mu_phi, vector<int> ll_lepid, vector<int> ll_i1,
    vector<int> ll_i2) {
  float dr1, dr2;
  if (ll_lepid[0]==11) {
    dr1 = get_dr(photon_eta[0],photon_phi[0],el_eta[ll_i1[0]],el_phi[ll_i1[0]]);
    dr2 = get_dr(photon_eta[0],photon_phi[0],el_eta[ll_i2[0]],el_phi[ll_i2[0]]);
    return dr1 > dr2 ? dr1 : dr2;
  }
  dr1 = get_dr(photon_eta[0],photon_phi[0],mu_eta[ll_i1[0]],mu_phi[ll_i1[0]]);
  dr2 = get_dr(photon_eta[0],photon_phi[0],mu_eta[ll_i2[0]],mu_phi[ll_i2[0]]);
  return dr1 > dr2 ? dr1 : dr2;
}

float get_l1_eta(vector<float> el_pt, vector<float> el_eta, 
    vector<float> mu_pt, vector<float> mu_eta, vector<int> ll_lepid, 
    vector<int> ll_i1, vector<int> ll_i2) {
  if (ll_lepid[0]==11) {
    return (el_pt[ll_i1[0]] > el_pt[ll_i2[0]]) ? el_eta[ll_i1[0]] : el_eta[ll_i2[0]];
  }
  return (mu_pt[ll_i1[0]] > mu_pt[ll_i2[0]]) ? mu_eta[ll_i1[0]] : mu_eta[ll_i2[0]];
}

float get_l2_eta(vector<float> el_pt, vector<float> el_eta, 
    vector<float> mu_pt, vector<float> mu_eta, vector<int> ll_lepid, 
    vector<int> ll_i1, vector<int> ll_i2) {
  if (ll_lepid[0]==11) {
    return (el_pt[ll_i1[0]] > el_pt[ll_i2[0]]) ? el_eta[ll_i2[0]] : el_eta[ll_i1[0]];
  }
  return (mu_pt[ll_i1[0]] > mu_pt[ll_i2[0]]) ? mu_eta[ll_i2[0]] : mu_eta[ll_i1[0]];
}

float get_l1_pt(vector<float> el_pt,
    vector<float> mu_pt, vector<int> ll_lepid, 
    vector<int> ll_i1, vector<int> ll_i2) {
  if (ll_lepid[0]==11) {
    return (el_pt[ll_i1[0]] > el_pt[ll_i2[0]]) ? el_pt[ll_i1[0]] : el_pt[ll_i2[0]];
  }
  return (mu_pt[ll_i1[0]] > mu_pt[ll_i2[0]]) ? mu_pt[ll_i1[0]] : mu_pt[ll_i2[0]];
}

float get_l2_pt(vector<float> el_pt,
    vector<float> mu_pt, vector<int> ll_lepid, 
    vector<int> ll_i1, vector<int> ll_i2) {
  if (ll_lepid[0]==11) {
    return (el_pt[ll_i1[0]] > el_pt[ll_i2[0]]) ? el_pt[ll_i2[0]] : el_pt[ll_i1[0]];
  }
  return (mu_pt[ll_i1[0]] > mu_pt[ll_i2[0]]) ? mu_pt[ll_i2[0]] : mu_pt[ll_i1[0]];
}

float get_l1_phi(vector<float> el_pt, vector<float> el_phi, 
    vector<float> mu_pt, vector<float> mu_phi, vector<int> ll_lepid, 
    vector<int> ll_i1, vector<int> ll_i2) {
  if (ll_lepid[0]==11) {
    return (el_pt[ll_i1[0]] > el_pt[ll_i2[0]]) ? el_phi[ll_i1[0]] : el_phi[ll_i2[0]];
  }
  return (mu_pt[ll_i1[0]] > mu_pt[ll_i2[0]]) ? mu_phi[ll_i1[0]] : mu_phi[ll_i2[0]];
}

float get_l2_phi(vector<float> el_pt, vector<float> el_phi, 
    vector<float> mu_pt, vector<float> mu_phi, vector<int> ll_lepid, 
    vector<int> ll_i1, vector<int> ll_i2) {
  if (ll_lepid[0]==11) {
    return (el_pt[ll_i1[0]] > el_pt[ll_i2[0]]) ? el_phi[ll_i2[0]] : el_phi[ll_i1[0]];
  }
  return (mu_pt[ll_i1[0]] > mu_pt[ll_i2[0]]) ? mu_phi[ll_i2[0]] : mu_phi[ll_i1[0]];
}

float H_t(vector<float> photon_pt,vector<float> photon_eta,vector<float> photon_phi,
          vector<float> el_pt,vector<float> el_eta,vector<float> el_phi,
          vector<float> mu_pt,vector<float> mu_eta,vector<float> mu_phi,
          vector<float> jet_pt,vector<float> jet_eta,vector<float> jet_phi,vector<float> jet_m){

    TLorentzVector tot,temp;
    tot.SetPtEtaPhiM(photon_pt[0],photon_eta[0],photon_phi[0],0.0);

    for(unsigned int idx = 0; idx < jet_pt.size(); idx++){
      temp.SetPtEtaPhiM(jet_pt[idx],jet_eta[idx],jet_phi[idx],jet_m[idx]);
      tot = tot + temp;
    }

    for(unsigned int idx = 0; idx < el_pt.size(); idx++){
      temp.SetPtEtaPhiM(el_pt[idx],el_eta[idx],el_phi[idx],0.000511);
      tot = tot + temp;
    }

    for(unsigned int idx = 0; idx < mu_pt.size(); idx++){
      temp.SetPtEtaPhiM(mu_pt[idx],mu_eta[idx],mu_phi[idx],0.1057);
      tot = tot + temp;
    }

    return tot.Et();

}

float get_st(vector<float> photon_pt,vector<bool> photon_sig,
          vector<float> el_pt,vector<float> el_sig,
          vector<float> mu_pt,vector<float> mu_sig,
          vector<float> jet_pt,vector<bool> jet_isgood){

    float result = 0;

    for(unsigned int idx=0; idx < photon_pt.size(); ++idx) {
      if (!photon_sig[idx]) continue;
      result += photon_pt[idx];
    }

    for(unsigned int idx = 0; idx < jet_pt.size(); idx++){
      if (!jet_isgood[idx]) continue;
      result += jet_pt[idx];
    }

    for(unsigned int idx = 0; idx < el_pt.size(); idx++){
      if (!el_sig[idx]) continue;
      result += el_pt[idx];
    }

    for(unsigned int idx = 0; idx < mu_pt.size(); idx++){
      if (!mu_sig[idx]) continue;
      result += mu_pt[idx];
    }

    return result;
}
/*
float get_weight(float w_lumi ,float w_year, vector<float> llphoton_l1_masserr,vector<float> llphoton_l2_masserr,vector<float> llphoton_ph_masserr, bool isNotSig){
    float dm = 1.0;
    if(isNotSig){dm=1;} else {
      float dml1,dml2,dmph;
      dml1 = llphoton_l1_masserr[0];
      dml2 = llphoton_l2_masserr[0];
      dmph = llphoton_ph_masserr[0];
      dm = sqrt(dml1 * dml1 + dml2 * dml2 + dmph * dmph);
    }
}
*/

float get_j1_pt(int njet, vector<int> jet_isgood, vector<float> jet_pt){
  if (njet<1) return -999;
  for (unsigned iPart = 0; iPart<jet_pt.size(); iPart++) {
      if (jet_isgood.at(iPart)) return jet_pt.at(iPart);
  }
  return -999;
}
float get_j1_eta(int njet, vector<int> jet_isgood, vector<float> jet_eta){
  if (njet<1) return -999;
  for (unsigned iPart = 0; iPart<jet_eta.size(); iPart++) {
      if (jet_isgood.at(iPart)) return jet_eta.at(iPart);
  }
  return -999;
}
float get_j1_phi(int njet, vector<int> jet_isgood, vector<float> jet_phi){
  if (njet<1) return -999;
  for (unsigned iPart = 0; iPart<jet_phi.size(); iPart++) {
      if (jet_isgood.at(iPart)) return jet_phi.at(iPart);
  }
  return -999;
}
float get_llyj_dphi(float llphoton_phi, int njet, vector<int> jet_isgood, vector<float> jet_phi){
  if (njet<1) return -999;
  float signal_jet_phi = -999;
  for (unsigned iPart = 0; iPart<jet_phi.size(); iPart++) {
      if (jet_isgood.at(iPart)) {
        signal_jet_phi = jet_phi.at(iPart);
        break;
      }
  }
  if (signal_jet_phi<-998) return signal_jet_phi;
  float dphi_lead = TVector2::Phi_mpi_pi(llphoton_phi - signal_jet_phi);
  return dphi_lead;
}

float get_yj1_deta(float photon_eta, int njet, vector<int> jet_isgood, vector<float> jet_eta){
  if (njet<1) return -999;
  float signal_jet_eta = -999;
  for (unsigned iPart = 0; iPart<jet_eta.size(); iPart++) {
      if (jet_isgood.at(iPart)) {
        signal_jet_eta = jet_eta.at(iPart);
        break;
      }
  }
  if (signal_jet_eta<-998) return signal_jet_eta;
  float deta_lead = photon_eta - signal_jet_eta;
  return deta_lead;
}
float get_yj1_dr(float photon_eta, float photon_phi, int njet, vector<int> jet_isgood, vector<float> jet_eta, vector<float> jet_phi){
  if (njet<1) return -999;
  float signal_jet_eta = -999;
  float signal_jet_phi = -999;
  for (unsigned iPart = 0; iPart<jet_eta.size(); iPart++) {
      if (jet_isgood.at(iPart)) {
        signal_jet_eta = jet_eta.at(iPart);
        signal_jet_phi = jet_phi.at(iPart);
        break;
      }
  }
  if (signal_jet_eta<-998) return signal_jet_eta;
  return get_dr(photon_eta, photon_phi, signal_jet_eta, signal_jet_phi);
}
float get_llyj1_ptbal(float ll_pt, float ll_eta, float ll_phi, float photon_pt, float photon_eta, float photon_phi, int njet, vector<int> jet_isgood, vector<float> jet_pt, vector<float> jet_eta, vector<float> jet_phi){
  if (njet<1) return -999;
  float signal_jet_eta = -999;
  float signal_jet_phi = -999;
  float signal_jet_pt = -999;
  for (unsigned iPart = 0; iPart<jet_eta.size(); iPart++) {
      if (jet_isgood.at(iPart)) {
        signal_jet_pt = jet_pt.at(iPart);
        signal_jet_eta = jet_eta.at(iPart);
        signal_jet_phi = jet_phi.at(iPart);
        break;
      }
  }
  if (signal_jet_eta<-998) return signal_jet_eta;
  TVector3 zboson; zboson.SetPtEtaPhi(ll_pt, ll_eta, ll_phi);
  TVector3 gamma; gamma.SetPtEtaPhi(photon_pt, photon_eta, photon_phi);
  TVector3 jet; jet.SetPtEtaPhi(signal_jet_pt, signal_jet_eta, signal_jet_phi);
  return (zboson+gamma+jet).Pt()/(zboson.Pt()+gamma.Pt()+jet.Pt());
}
float get_j2_pt(int njet, vector<int> jet_isgood, vector<float> jet_pt){
  if (njet<2) return -999;
  bool sublead = false;
  for (unsigned iPart = 0; iPart<jet_pt.size(); iPart++) {
    if (jet_isgood.at(iPart) == 0) continue;
    if (sublead == true) return jet_pt.at(iPart);
    sublead = true;
  }
  return -999;
}

float get_llg_ptt(vector<float> photon_pt, vector<float> photon_eta, vector<float> photon_phi, 
                  vector<float> llphoton_pt, vector<float> llphoton_eta, vector<float> llphoton_phi,
                  vector<float> ll_pt, vector<float> ll_eta, vector<float> ll_phi) {
  TVector3 gamma; gamma.SetPtEtaPhi(photon_pt[0], photon_eta[0], photon_phi[0]);
  TVector3 higgs; higgs.SetPtEtaPhi(llphoton_pt[0], llphoton_eta[0], llphoton_phi[0]);
  TVector3 zboson; zboson.SetPtEtaPhi(ll_pt[0], ll_eta[0], ll_phi[0]);
  gamma.SetZ(0); higgs.SetZ(0); zboson.SetZ(0);
  return higgs.Cross((zboson-gamma).Unit()).Mag();
}

float get_tru_leplep_m(vector<float> mc_id, vector<float> mc_status, vector<float> mc_pt, vector<float> mc_eta, vector<float> mc_phi, vector<float> mc_mass) {
  float leplep_m = -999;
  TLorentzVector lep_plus;
  TLorentzVector lep_minus;
  bool set_lep_plus = false;
  bool set_lep_minus = false;
  for (unsigned imc = 0; imc < mc_id.size(); ++imc) {
    if (mc_id[imc]==23) {
      leplep_m = mc_mass[imc];
      break;
    }
    if (mc_status[imc]==23)  {
      if (mc_id[imc]==11 || mc_id[imc]==13) {
        lep_plus.SetPtEtaPhiM(mc_pt[imc],mc_eta[imc],mc_phi[imc],mc_mass[imc]);
        set_lep_plus = true;
      } else if (mc_id[imc]==-11 || mc_id[imc]==-13) {
        lep_minus.SetPtEtaPhiM(mc_pt[imc],mc_eta[imc],mc_phi[imc],mc_mass[imc]);
        set_lep_minus = true;
      }
    }
    if (set_lep_plus && set_lep_minus) break;
  }
  if (leplep_m>-998) return leplep_m;
  if (set_lep_minus && set_lep_plus) return (lep_plus+lep_minus).M();
  return -999;
}






  //This code defines the triggers and trigger pt selections.
  const NamedFunc trigs_2016_e = "HLT_Ele27_WPTight_Gsf";                  const NamedFunc trigs_2016_ee = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ";
  const NamedFunc trigs_2017_e = "HLT_Ele32_WPTight_Gsf_Emu";              const NamedFunc trigs_2017_ee = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"; //Use for kingscanyon_v1
  const NamedFunc trigs_2018_e = "HLT_Ele32_WPTight_Gsf";                  const NamedFunc trigs_2018_ee = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL";
  const NamedFunc trigs_2022_e = "HLT_Ele30_WPTight_Gsf";                  const NamedFunc trigs_2022_ee = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL";
  const NamedFunc trigs_2023_e = "HLT_Ele30_WPTight_Gsf";                  const NamedFunc trigs_2023_ee = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL";

  const NamedFunc trigs_2016_mu = "HLT_IsoMu24 || HLT_IsoTkMu24";          const NamedFunc trigs_2016_mumu = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ";
  const NamedFunc trigs_2017_mu = "HLT_IsoMu27";                           const NamedFunc trigs_2017_mumu = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8";
  const NamedFunc trigs_2018_mu = "HLT_IsoMu24";                           const NamedFunc trigs_2018_mumu = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8";
  const NamedFunc trigs_2022_mu = "HLT_IsoMu24";                           const NamedFunc trigs_2022_mumu = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8";
  const NamedFunc trigs_2023_mu = "HLT_IsoMu24";                           const NamedFunc trigs_2023_mumu = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8";

  const NamedFunc trigs_pT_2016_el = (trigs_2016_e && lead_el_pt > 30) || (trigs_2016_ee && lead_el_pt > 25 && sublead_el_pt > 15);
  const NamedFunc trigs_pT_2016_mu = (trigs_2016_mu && lead_mu_pt > 25) || (trigs_2016_mumu && lead_mu_pt > 20 && sublead_mu_pt > 10);
  const NamedFunc trigs_pT_2016 = trigs_pT_2016_el || trigs_pT_2016_mu;

  const NamedFunc trigs_pT_2017_el = (trigs_2017_e && lead_el_pt > 35) || (trigs_2017_ee && lead_el_pt > 25 && sublead_el_pt > 15); 
  const NamedFunc trigs_pT_2017_mu = (trigs_2017_mu && lead_mu_pt > 28) || (trigs_2017_mumu && lead_mu_pt > 20 && sublead_mu_pt > 10);
  const NamedFunc trigs_pT_2017 = trigs_pT_2017_el || trigs_pT_2017_mu;

  const NamedFunc trigs_pT_2018_el = (trigs_2018_e && lead_el_pt > 35) || (trigs_2018_ee && lead_el_pt > 25 && sublead_el_pt > 15);
  const NamedFunc trigs_pT_2018_mu = (trigs_2018_mu && lead_mu_pt > 25) || (trigs_2018_mumu && lead_mu_pt > 20 && sublead_mu_pt > 10);
  const NamedFunc trigs_pT_2018 = trigs_pT_2018_el || trigs_pT_2018_mu;

  const NamedFunc trigs_pT_2022_el = (trigs_2022_e && lead_el_pt > 35) || (trigs_2022_ee && lead_el_pt > 25 && sublead_el_pt > 15);
  const NamedFunc trigs_pT_2022_mu = (trigs_2022_mu && lead_mu_pt > 25) || (trigs_2022_mumu && lead_mu_pt > 20 && sublead_mu_pt > 10);
  const NamedFunc trigs_pT_2022 = trigs_pT_2022_el || trigs_pT_2022_mu;

  const NamedFunc trigs_pT_2023_el = (trigs_2023_e && lead_el_pt > 35) || (trigs_2023_ee && lead_el_pt > 25 && sublead_el_pt > 15);
  const NamedFunc trigs_pT_2023_mu = (trigs_2023_mu && lead_mu_pt > 25) || (trigs_2023_mumu && lead_mu_pt > 20 && sublead_mu_pt > 10);
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
    if (electron_or_muon == 0) { // el
      tablerows = vector<TableRow>{
        
        TableRow("Initial events", add_cut(current_cut, "1"),0,0, weight),
        //TableRow("printer", current_cut && printer2,0,0,weight),
        TableRow("$N_e \\geq 2$", add_cut(current_cut, "nll>0 && nel>=2 && ll_lepid[0]==11"),0,0, weight),
        TableRow("e, ee trigger", add_cut(current_cut, "trig_double_el==1 || trig_single_el==1"),0,0, weight),
        TableRow("el trigger $p_{T}$ cuts", add_cut(current_cut, "trig_el_pt==1"),0,0, weight),
        TableRow("$N_{\\gamma} \\geq 1$", add_cut(current_cut, "nphoton>=1"),0,0, weight),
/*
        TableRow("cosTheta>0.5f", current_cut && cosTheta_n2p, 0,0, weight),
        TableRow("costheta>0.5f", current_cut && costheta_n2p, 0,0, weight),
        TableRow("phi (psi)>0.5f", current_cut && phipsi, 0,0, weight),
        TableRow("Photon_mva>0.5f", current_cut && "photon_idmva>0.5f", 0,0, weight),
        TableRow("min_dR>0.3"
*/
        TableRow("$80 \\text{ GeV} < m_{ee} < 100 \\text{ GeV}$", add_cut(current_cut, "ll_m[0]>80 && ll_m[0]<100"),0,0, weight),
        TableRow("$p_{T}^{\\gamma}/m_{ee\\gamma} > 15./110$", add_cut(current_cut, "photon_pt[llphoton_iph[0]]/llphoton_m[0]>(15./110)"),0,0, weight),
        TableRow("$m_{ee} + m_{ee\\gamma} > 185 \\text{ GeV}$", add_cut(current_cut, "ll_m[0]+llphoton_m[0] > 185"),0,0, weight),
        TableRow("$100 < m_{ee\\gamma} < 180 \\text{ GeV}$", add_cut(current_cut, "llphoton_m[0] > 100 && llphoton_m[0] < 180"),0,0, weight),
        TableRow("Event filters", add_cut(current_cut, "pass==1"),0,0, weight),
        //Category selections
        TableRow("ttH lep.", current_cut && max_mini_iso<0.1 && "njet>=3" && "ll_m[0]>85 && ll_m[0]<95" && "(nlep==3 && nbdfm>=1) || (nlep>=4 && njet>=1 && nbdfm>=1)",0,0, weight),
        TableRow("ttH had.", current_cut && "ll_m[0]>85 && ll_m[0]<95" && "njet>=5" &&"nlep==2 && nbdfm>=1",0,0, weight),
        TableRow("ZH $p_{T}^{miss}$", current_cut && "llphoton_pt[0]/llphoton_m[0]>0.4" && "ll_m[0]>85 && ll_m[0]<95" && "njet<=1" && "nlep==2 && met>90",0,0, weight),
        TableRow("WH/ZH $N_{l}\\geq 3$", current_cut && max_mini_iso<0.15 && "llphoton_pt[0]/llphoton_m[0]>0.3" && "ll_m[0]>85 && ll_m[0]<95" && "nlep>=3 && nbdfm==0 && met>30",0,0, weight),
        TableRow("VBF", add_cut(current_cut, ("njet>=2" && "nlep==2 && nbdfm==0")),0,0, weight),
        TableRow("ggF", current_cut && "njet<=1" && "nlep==2 && met<90",0,0, weight),

      };
    } else if (electron_or_muon == 1) { // mu
      tablerows = vector<TableRow>{
        TableRow("Initial events", add_cut(current_cut, "1"),0,0, weight),
        TableRow("$N_{\\mu} \\geq 2 $", add_cut(current_cut, "nll>0 && nmu>=2 && ll_lepid[0]==13"),0,0, weight),
        TableRow("$\\mu, \\mu\\mu$ trigger", add_cut(current_cut, "trig_double_mu==1 || trig_single_mu==1"),0,0, weight),
        TableRow("mu trigger $p_{T}$ cuts", add_cut(current_cut, "trig_mu_pt==1"),0,0, weight),
        TableRow("$N_{\\gamma} \\geq 1$", add_cut(current_cut, "nphoton>=1"),0,0, weight),
        TableRow("$80 \\text{ GeV} < m_{\\mu\\mu} < 100 \\text{ GeV}$", add_cut(current_cut, "ll_m[0]>80 && ll_m[0]<100"),0,0, weight),
        TableRow("$p_{T}^{\\gamma}/m_{\\mu\\mu\\gamma} > 15./110$", add_cut(current_cut, "photon_pt[llphoton_iph[0]]/llphoton_m[0]>(15./110)"),0,0, weight),
        TableRow("$m_{\\mu\\mu\\gamma}+m_{\\mu\\mu} > 185 \\text{ GeV}$", add_cut(current_cut, "ll_m[0]+llphoton_m[0] > 185"),0,0, weight),
        TableRow("$100 < m_{\\mu\\mu\\gamma} < 180 \\text{ GeV}$", add_cut(current_cut, "llphoton_m[0] > 100 && llphoton_m[0] < 180"),0,0, weight),
        TableRow("Event filters", add_cut(current_cut, "pass==1"),0,0, weight),
        //Category selections
        TableRow("ttH lep.", current_cut && max_mini_iso<0.1 && "ll_m[0]>85 && ll_m[0]<95" && "(nlep==3 && njet>=3 && nbdfm>=1) || (nlep>=4 && njet>=1 && nbdfm>=1)",0,0, weight),
        TableRow("ttH had.", current_cut && "ll_m[0]>85 && ll_m[0]<95" && "nlep==2 && njet>=5 && nbdfm>=1",0,0, weight),
        TableRow("ZH $p_{T}^{miss}$", current_cut && "llphoton_pt[0]/llphoton_m[0]>0.4" && "ll_m[0]>85 && ll_m[0]<95" && "nlep==2 && njet<=1 && met>90",0,0, weight),
        TableRow("WH/ZH $N_{l}\\geq 3$", current_cut && max_mini_iso<0.15 && "llphoton_pt[0]/llphoton_m[0]>0.3" && "ll_m[0]>85 && ll_m[0]<95" && "nlep>=3 && nbdfm==0 && met>30",0,0, weight),
        TableRow("VBF", current_cut && "nlep==2 && njet>=2 && nbdfm==0",0,0, weight),
        TableRow("ggF", current_cut && "nlep==2 && njet<=1 && met<90",0,0, weight),

      };
    } else { // mu + el
      tablerows = vector<TableRow>{
        TableRow("Initial events", add_cut(current_cut,"1"),0,0, weight),
        TableRow("$N_e \\geq 2 || N_{\\mu} \\geq 2 $", add_cut(current_cut, "nll>0 && ((nel>=2 && ll_lepid[0]==11) || (nmu>=2 && ll_lepid[0]==13))"),0,0, weight),
        TableRow("e, ee trigger $||$ $\\mu, \\mu\\mu$ trigger", add_cut(current_cut, "trig_double_el || trig_single_el || trig_double_mu==1 || trig_single_mu==1"),0,0, weight),
        TableRow("mu and el trigger $p_{T}$ cuts", add_cut(current_cut, "trig_el_pt==1 || trig_mu_pt==1"),0,0, weight),
        TableRow("$N_{\\gamma} \\geq 1$", add_cut(current_cut, "nphoton>=1"),0,0, weight),
        TableRow("$80 \\text{ GeV} < m_{ll} < 100 \\text{ GeV}$", add_cut(current_cut, "ll_m[0]>80 && ll_m[0]<100"),0,0, weight),
        TableRow("$p_{T}^{\\gamma}/m_{ll\\gamma} > 15./110$", add_cut(current_cut, "photon_pt[0]/llphoton_m[0]>(15./110)"),0,0, weight),
        TableRow("$m_{ll\\gamma}+m_{ll} > 185 \\text{ GeV}$", add_cut(current_cut, "ll_m[0]+llphoton_m[0]>185"),0,0, weight),
        TableRow("$100 < m_{ll\\gamma} < 180 \\text{ GeV}$", add_cut(current_cut, "llphoton_m[0]>100 && llphoton_m[0]<180"),0,0, weight),
        TableRow("Event filters", add_cut(current_cut, "pass==1"),0,0, weight),
        //Category selections
        TableRow("ttH lep.", current_cut && max_mini_iso<0.1 && "ll_m[0]>85 && ll_m[0]<95" && "(nlep==3 && njet>=3 && nbdfm>=1) || (nlep>=4 && njet>=1 && nbdfm>=1)",0,0, weight),
        TableRow("ttH had.", current_cut && "ll_m[0]>85 && ll_m[0]<95" && "nlep==2 && njet>=5 && nbdfm>=1",0,0, weight),
        TableRow("ZH $p_{T}^{miss}$", current_cut && "llphoton_pt[0]/llphoton_m[0]>0.4" && "ll_m[0]>85 && ll_m[0]<95" && "nlep==2 && njet<=1 && met>90",0,0, weight),
        TableRow("WH/ZH $N_{l}\\geq 3$", current_cut && max_mini_iso<0.15 && "llphoton_pt[0]/llphoton_m[0]>0.3" && "ll_m[0]>85 && ll_m[0]<95" && "nlep>=3 && nbdfm==0 && met>30",0,0, weight),
        TableRow("VBF", current_cut && "nlep==2 && njet>=2 && nbdfm==0",0,0, weight),
        TableRow("ggF", current_cut && "nlep==2 && njet<=1 && met<90",0,0, weight),

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
    filenames = {"*ZGToLLG_01J_5f_lowMLL_lowGPt_TuneCP5_13TeV*0F18644F-6B57-E543-B9F3-862E4E747E90.root"};
    type = "mc/"; label = "DYJets and SMZy"; color = "t1tttt"; print_uncertainty = true;
    pathNames = attach_folder(path_to_production_a, year_t, type + skim_a, filenames); 
    procs_dat_a.push_back(Process::MakeShared<Baby_pico>(label, Process::Type::signal, colors(color), pathNames, "1")); 
  } else if(isMC == true && isSig == true){
    filenames = {"*GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV*.root"};
    type = "mc/"; label = "Signal MC"; color = "dy"; print_uncertainty = true;
    pathNames = attach_folder(path_to_production_a, year_t, type + skim_a, filenames); 
    procs_dat_a.push_back(Process::MakeShared<Baby_pico>(label, Process::Type::background, colors(color), pathNames, "1"));
  } else if(isMC == false){
    filenames = {"*SingleElectron*.root"};
    type = "data/"; label = "Data"; color = "data"; print_uncertainty = false;
    pathNames = attach_folder(path_to_production_a, year_t, type + skim_a, filenames); 
    procs_dat_a.push_back(Process::MakeShared<Baby_pico>(label, Process::Type::data, colors(color), pathNames, "1"));
  }
  vector<TableRow> tablerows_a;
  constructCutflowTable(tablerows_a, weight_a, lepton, isMC);
  pm.Push<Table>(tablename_a.Data(), tablerows_a, procs_dat_a, 0, 1, 0, 1, 0, print_uncertainty, 0).LuminosityTag(total_luminosity_string_a).Precision(3);
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
  string path_to_production_3 = path_to_production_base +"NanoAODv12/" + tag_a + "/";
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
    if (Contains(years_string_a[k], "2023") || Contains(years_string_a[k], "2022")){
      path_to_production_a = path_to_production_3;
    } else {
      path_to_production_a = path_to_production_2;
    }
    for (unsigned int l = 0; l<=2; l++){
      //Loop over lepton types, 0 = electron, 1 = muon, 2 = electron and muon
        if (l==0){
          lep_tag = "_electron";
          continue;
        }
        else if (l==1){
          lep_tag = "_muon";
        }
        else{
          lep_tag = "_lepton";
          continue;
        }

        TString tablename_a = "FixName:validation_table_"+tag_a+lep_tag+"_"+years_string_a[k];
        string total_luminosity_string_a = getLuminosityString(years_string_a[k]);
        year_t = {years_string_a[k]};

/*
        //    Generate data cutflow table
        //--------------------------------------------------
        tablename_a = "FixName:validation_table_"+tag_a+lep_tag+"_"+years_string_a[k]+"_data";
        weight_a = "1";
        generateTable(path_to_production_a, pm, colors, year_t, skim_a, weight_a, tablename_a, total_luminosity_string_a, lumi, l, false, false);


        //    Generate background MC (no SF) cutflow table
        //--------------------------------------------------
        tablename_a = "FixName:validation_table_"+tag_a+lep_tag+"_"+years_string_a[k]+"_backmc";
        //weight_a = w_lumi*w_years;
        weight_a = "1";
        generateTable(path_to_production_a, pm, colors, year_t, skim_a, weight_a, tablename_a, total_luminosity_string_a, lumi, l, true, false);


        //    Generate background MC (w/ SF) cutflow table
        //--------------------------------------------------
        tablename_a = "FixName:validation_table_"+tag_a+lep_tag+"_"+years_string_a[k]+"_backmc_sf";
        weight_a = wgt*w_years;
        generateTable(path_to_production_a, pm, colors, year_t, skim_a, weight_a, tablename_a, total_luminosity_string_a, lumi, l, true, false);
*/
        //    Generate signal MC (no SF) cutflow table
        //--------------------------------------------------
        tablename_a = "FixName:validation_table_"+tag_a+lep_tag+"_"+years_string_a[k]+"_sigmc";
        weight_a = "1";
        generateTable(path_to_production_a, pm, colors, year_t, skim_a, weight_a, tablename_a, total_luminosity_string_a, lumi, l, true, true);
/*
        //    Generate signal MC (no SF) cutflow table
        //--------------------------------------------------
        tablename_a = "FixName:validation_table_"+tag_a+lep_tag+"_"+years_string_a[k]+"_sigmc_sf";
        weight_a = w_lumi*w_years;
        generateTable(path_to_production_a, pm, colors, year_t, skim_a, weight_a, tablename_a, total_luminosity_string_a, lumi, l, true, true);
*/
    }
  }



}

