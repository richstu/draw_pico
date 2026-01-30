//This script runs over Nanos to generate some plots to investigate electron-photon shower shape differences
//

#include "core/test.hpp"

#include <algorithm>
#include <bitset>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <unistd.h>
#include <getopt.h>

#include "TColor.h"
#include "TError.h"
#include "TLorentzVector.h"

#include "core/baby.hpp"
#include "core/efficiency_plot.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/named_func.hpp"
#include "core/named_func_utilities.hpp"
#include "core/palette.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/process.hpp"
#include "core/table.hpp"
#include "core/utilities.hpp"

#include "zgamma/nano_functions.hpp"
#include "zgamma/zg_functions.hpp"
#include "zgamma/zg_utilities.hpp"

using std::vector;
using std::string;
using std::shared_ptr;
using std::to_string;
using NamedFuncUtilities::reduce_max;
using NamedFuncUtilities::reduce_maxfirst;
using NamedFuncUtilities::reduce_sum;
using NamedFuncUtilities::FilterNamedFunc;
using NamedFuncUtilities::MapNamedFunc;
using NamedFuncUtilities::MultiReduceNamedFunc;
using NamedFuncUtilities::ReduceNamedFunc;
using NanoFunctions::Electron_sig;
using NanoFunctions::nSignalElectron;
using NanoFunctions::Muon_sig;
using NanoFunctions::Lead_SignalMuon_pt;
using NanoFunctions::Lead_SignalMuon_eta;
using NanoFunctions::Lead_SignalMuon_phi;
using NanoFunctions::Sublead_SignalMuon_pt;
using NanoFunctions::Sublead_SignalMuon_eta;
using NanoFunctions::Sublead_SignalMuon_phi;
using NanoFunctions::nSignalMuon;
using NanoFunctions::Photon_sig;
using NanoFunctions::Photon_abseta;
using NanoFunctions::SignalPhoton_pt;
using NanoFunctions::Lead_SignalPhoton_pt;
using NanoFunctions::Lead_SignalPhoton_abseta;
using NanoFunctions::Lead_SignalPhoton_eta;
using NanoFunctions::Lead_SignalPhoton_phi;
using NanoFunctions::Lead_SignalPhoton_mvaID;
using NanoFunctions::nSignalPhoton;
using NanoFunctions::HiggsCandidate_mass;
const Process::Type background = Process::Type::background;

int main(){
  //------------------------------------------------------------------------------------
  //                              constants and namedfuncs
  //------------------------------------------------------------------------------------
  
  //Signal flag for photons indicating HIG-19-014 criteria except electron veto
  //and deltaR from negative electrons
  const NamedFunc Photon_looseSig("Photon_looseSig",[](const Baby &b) -> NamedFunc::VectorType{
    vector<double> ph_sig;
    vector<double> Electron_signal = Electron_sig.GetVector(b);
    for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
      if (b.Photon_pt()->at(iph) > 15 && 
          ((b.Photon_isScEtaEB()->at(iph) && b.Photon_mvaID()->at(iph) > -0.4) ||
           (b.Photon_isScEtaEE()->at(iph) && b.Photon_mvaID()->at(iph) > -0.58))) {
        bool nearby_el = false;
        for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
          if (Electron_signal[iel] != 0 && b.Electron_charge()->at(iel)==1) {
            if (deltaR(b.Electron_eta()->at(iel), b.Electron_phi()->at(iel),
                       b.Photon_eta()->at(iph), b.Photon_phi()->at(iph))<0.4) {
              nearby_el = true;
            }
          }
        } //loop over electrons
        if (!nearby_el) {
          ph_sig.push_back(1);
        }
        else {
          ph_sig.push_back(0);
        }
      }
      else {
        ph_sig.push_back(0);
      }
    } //loop over photons
    return ph_sig;
  });

  //number of photons passing HIG-19-014 criteria except electron veto and deltaR 
  //from negative electron
  const NamedFunc nLoosePhoton = ReduceNamedFunc(Photon_looseSig, reduce_sum).Name("nLoosePhoton");
  
  //NamedFunc that is true when there is a signal positron passing trigger and pT cuts as well as
  //a photon passing HIG-19-014 criteria except the electron veto and deltaR from negative electrons
  const NamedFunc Pass_ee_tnp("Pass_ee_tnp",[nLoosePhoton](const Baby &b) -> NamedFunc::ScalarType{

    if (!b.HLT_Ele32_WPTight_Gsf()) return false;

    //check for signal positron passing trigger and pT cuts
    vector<double> Electron_signal = Electron_sig.GetVector(b);
    bool found_good_el = false;
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_signal[iel] != 0 && b.Electron_charge()->at(iel)==1 &&
          b.Electron_pt()->at(iel) > 35) {
        for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
          if (b.TrigObj_id()->at(itrig) == 11 && 
              (b.TrigObj_filterBits()->at(itrig) & 0x2) != 0) {
            if (deltaR(b.Electron_eta()->at(iel), b.Electron_phi()->at(iel),
                       b.TrigObj_eta()->at(itrig), b.TrigObj_phi()->at(itrig))<0.2) {

              found_good_el = true;
            }
          }
        } //loop over TrigObjs
      }
    } //loop over electrons

    if (!found_good_el) return false;
    return (nLoosePhoton.GetScalar(b) > 0);
  });

  //loose photon pt
  const NamedFunc LoosePhoton_pt = FilterNamedFunc("Photon_pt",Photon_looseSig)
                                   .Name("LoosePhoton_pt");

  //lead loose photon pt
  const NamedFunc Lead_LoosePhoton_pt = ReduceNamedFunc(LoosePhoton_pt, reduce_max).Name("Lead_LoosePhoton_pt");

  //loose photon abs eta
  const NamedFunc LoosePhoton_abseta = FilterNamedFunc(Photon_abseta,Photon_looseSig)
                                   .Name("LoosePhoton_abseta");

  //Lead loose photon abseta
  const NamedFunc Lead_LoosePhoton_abseta = MultiReduceNamedFunc(
      {LoosePhoton_pt,LoosePhoton_abseta}, reduce_maxfirst).Name("Lead_LoosePhoton_abseta");

  //loose photon mvaID
  const NamedFunc LoosePhoton_mvaID = FilterNamedFunc("Photon_mvaID",Photon_looseSig)
                                      .Name("LoosePhoton_mvaID");

  //Lead loose photon IDMVA
  const NamedFunc Lead_LoosePhoton_mvaID = MultiReduceNamedFunc(
      {LoosePhoton_pt,LoosePhoton_mvaID}, reduce_maxfirst).Name("Lead_LoosePhoton_mvaID");

  //loose photon r9
  const NamedFunc LoosePhoton_r9 = FilterNamedFunc("Photon_r9",Photon_looseSig)
                                      .Name("LoosePhoton_r9");

  //Lead loose photon r9
  const NamedFunc Lead_LoosePhoton_r9 = MultiReduceNamedFunc(
      {LoosePhoton_pt,LoosePhoton_r9}, reduce_maxfirst).Name("Lead_LoosePhoton_r9");

  //loose photon sieie
  const NamedFunc LoosePhoton_sieie = FilterNamedFunc("Photon_sieie",Photon_looseSig)
                                      .Name("LoosePhoton_sieie");

  //Lead loose photon sieie
  const NamedFunc Lead_LoosePhoton_sieie = MultiReduceNamedFunc(
      {LoosePhoton_pt,LoosePhoton_sieie}, reduce_maxfirst).Name("Lead_LoosePhoton_sieie");

  //variation on H->yy preselection
  const NamedFunc Photon_hggPresel = NamedFunc("Photon_hoe<0.08&&Photon_sieie<0.015&&Photon_r9>0.5&&Photon_pfRelIso03_all<0.035").Name("Photon_hggPresel");

  //loose photon hggPresel
  const NamedFunc LoosePhoton_hggPresel = FilterNamedFunc(Photon_hggPresel,Photon_sig).Name("LoosePhoton_hggPresel");

  //lead loose photon hggPresel
  const NamedFunc Lead_LoosePhoton_hggPresel = MultiReduceNamedFunc(
      {LoosePhoton_pt,LoosePhoton_hggPresel}, reduce_maxfirst).Name("Lead_LoosePhoton_hggPresel");

  //Electron+"Photon" (i.e. Z candidate) invariant mass
  const NamedFunc ElectronPhoton_mass("ElectronPhoton_mass",[Photon_looseSig](const Baby &b) -> NamedFunc::ScalarType{
  
    vector<double> Electron_signal = Electron_sig.GetVector(b);
    vector<double> Photon_signal = Photon_looseSig.GetVector(b);
    float lead_el_pt = 0;
    float lead_el_eta = 0;
    float lead_el_phi = 0;
    float lead_ph_pt = 0;
    float lead_ph_eta = 0;
    float lead_ph_phi = 0;
    TLorentzVector el,ph;
  
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_signal[iel] != 0 && b.Electron_charge()->at(iel)==1) {
        if (b.Electron_pt()->at(iel) > lead_el_pt) {
          lead_el_pt = b.Electron_pt()->at(iel);
          lead_el_eta = b.Electron_eta()->at(iel);
          lead_el_phi = b.Electron_phi()->at(iel);
        }
      }
    } //loop over electrons

    for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
      if (Photon_signal[iph] != 0) {
        if (b.Photon_pt()->at(iph) > lead_ph_pt) {
          lead_ph_pt = b.Photon_pt()->at(iph);
          lead_ph_eta = b.Photon_eta()->at(iph);
          lead_ph_phi = b.Photon_phi()->at(iph);
        }
      }
    } //loop over photons

    el.SetPtEtaPhiM(lead_el_pt,lead_el_eta,lead_el_phi,0.000511);
    ph.SetPtEtaPhiM(lead_ph_pt,lead_ph_eta,lead_ph_phi,0.0);
    return (el+ph).M();
  });

  //signal photon hggPresel
  const NamedFunc SignalPhoton_hggPresel = FilterNamedFunc(Photon_hggPresel,Photon_sig).Name("SignalPhoton_hggPresel");

  //lead signal photon hggPresel
  const NamedFunc Lead_SignalPhoton_hggPresel = MultiReduceNamedFunc(
      {SignalPhoton_pt,SignalPhoton_hggPresel}, reduce_maxfirst).Name("Lead_SignalPhoton_hggPresel");

  //signal photon r9
  const NamedFunc SignalPhoton_r9 = FilterNamedFunc("Photon_r9",Photon_sig).Name("SignalPhoton_r9");

  //lead signal photon r9
  const NamedFunc Lead_SignalPhoton_r9 = MultiReduceNamedFunc(
      {SignalPhoton_pt,SignalPhoton_r9}, reduce_maxfirst).Name("Lead_SignalPhoton_r9");

  //signal photon sieie
  const NamedFunc SignalPhoton_sieie = FilterNamedFunc("Photon_sieie",Photon_sig).Name("SignalPhoton_sieie");

  //lead signal photon sieie
  const NamedFunc Lead_SignalPhoton_sieie = MultiReduceNamedFunc(
      {SignalPhoton_pt,SignalPhoton_sieie}, reduce_maxfirst).Name("Lead_SignalPhoton_sieie");

  //True when there are two signal muons, double muon trigger is fired, and a
  //signal photon
  const NamedFunc Pass_mumuph_tnp = NamedFunc(nSignalMuon==2 && nSignalPhoton==1 &&
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8").Name("Pass_mumuph_tnp");

  const NamedFunc w_lumi_dy = NamedFunc("Generator_weight/2.735e10*6077.22*1000.0")
                              .Name("w_lumi_dy");

  const NamedFunc w_lumi_zg = NamedFunc("Generator_weight/5.803e8*172.4*1000.0")
                              .Name("w_lumi_zg");

  const NamedFunc weight_dy = NamedFunc(w_lumi_dy*60.0).Name("weight_dy");
  const NamedFunc weight_zg = NamedFunc(w_lumi_zg*60.0).Name("weight_zg");

  //------------------------------------------------------------------------------------
  //                                 initialization
  //------------------------------------------------------------------------------------
  
  gErrorIgnoreLevel = 6000;

  string lumi_string = "60";
  string base_folder = "/net/cms11/cms11r0/pico/NanoAODv9/nano/2018/mc/";

  vector<shared_ptr<Process>> procs_ee;
  procs_ee.push_back(Process::MakeShared<Baby_nano>("Electron (Z#rightarrow ee)",
      Process::Type::background, TColor::GetColor("#5790fc"),
      {base_folder+"DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX*.root"}, Pass_ee_tnp));
      //{base_folder+"DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL18NanoAODv9__106X_upgrade2018_realistic_v16_L1v1-v2__230000__A9969A0B-4994-7943-AEB7-301174DE97AF.root"}, Pass_ee_tnp));

  vector<shared_ptr<Process>> procs_mumuph;
  procs_mumuph.push_back(Process::MakeShared<Baby_nano>("Photon (Z#rightarrow #mu#mu#gamma)",
      Process::Type::background, TColor::GetColor("#f89c20"),
      {base_folder+"ZGToLLG_01J_5f_lowMLL*.root"}, Pass_mumuph_tnp));
      //{base_folder+"ZGToLLG_01J_5f_lowMLL_lowGPt_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL18NanoAODv9__106X_upgrade2018_realistic_v16_L1v1-v1__2520000__231FBC7E-D6F7-EB4E-AA25-EF1C2871DA29.root"}, Pass_mumuph_tnp));

  std::vector<PlotOpt> ops_shapes = {PlotOpt("txt/plot_styles.txt","Shapes")};

  //------------------------------------------------------------------------------------
  //                                     make plots and pie charts
  //------------------------------------------------------------------------------------
  
  PlotMaker pm;

  //Electrons
  pm.Push<Hist1D>(Axis(40, -0.58, 1.0, Lead_LoosePhoton_mvaID, "Photon IDMVA", {}), 
      ElectronPhoton_mass>80&&ElectronPhoton_mass<100&&Lead_LoosePhoton_hggPresel, procs_ee, ops_shapes).Weight(weight_dy).Tag("zgshower");
  pm.Push<Hist1D>(Axis(40, -0.58, 1.0, Lead_LoosePhoton_mvaID, "Photon IDMVA", {}), 
      ElectronPhoton_mass>80&&ElectronPhoton_mass<100&&Lead_LoosePhoton_pt>27.5&&Lead_LoosePhoton_pt<30&&Lead_LoosePhoton_abseta<0.8&&Lead_LoosePhoton_hggPresel, procs_ee, ops_shapes).Weight(weight_dy).Tag("zgshower");
  pm.Push<Hist1D>(Axis(40, 0.5, 1.0, Lead_LoosePhoton_r9, "Photon R_{9}", {}), 
      ElectronPhoton_mass>80&&ElectronPhoton_mass<100&&Lead_LoosePhoton_hggPresel, procs_ee, ops_shapes).Weight(weight_dy).Tag("zgshower");
  pm.Push<Hist1D>(Axis(40, 0.0, 0.02, Lead_LoosePhoton_sieie, "Photon #sigma_{i#eta i#eta}", {}), 
      ElectronPhoton_mass>80&&ElectronPhoton_mass<100&&Lead_LoosePhoton_hggPresel, procs_ee, ops_shapes).Weight(weight_dy).Tag("zgshower");
  pm.Push<Hist1D>(Axis(40, 70.0, 110.0, ElectronPhoton_mass, "m_{e#gamma} [GeV]", {}), 
      Lead_LoosePhoton_hggPresel, procs_ee, ops_shapes).Weight(weight_dy).Tag("zgshower");
  pm.Push<Hist1D>(Axis(40, 70.0, 110.0, ElectronPhoton_mass, "m_{e#gamma} [GeV]", {}), 
      Lead_LoosePhoton_pt<20&&Lead_LoosePhoton_hggPresel, procs_ee, ops_shapes).Weight(weight_dy).Tag("zgshower");

  //Photons
  pm.Push<Hist1D>(Axis(40, -0.58, 1.0, Lead_SignalPhoton_mvaID, "Photon IDMVA", {}), 
      HiggsCandidate_mass>80&&HiggsCandidate_mass<100&&Lead_SignalPhoton_hggPresel, procs_mumuph, ops_shapes).Weight(weight_zg).Tag("zgshower");
  pm.Push<Hist1D>(Axis(40, 0.5, 1.0, Lead_SignalPhoton_r9, "Photon R_{9}", {}), 
      HiggsCandidate_mass>80&&HiggsCandidate_mass<100&&Lead_SignalPhoton_hggPresel, procs_mumuph, ops_shapes).Weight(weight_zg).Tag("zgshower");
  pm.Push<Hist1D>(Axis(40, 0.0, 0.02, Lead_SignalPhoton_sieie, "Photon #sigma_{i#eta i#eta}", {}), 
      HiggsCandidate_mass>80&&HiggsCandidate_mass<100&&Lead_SignalPhoton_hggPresel, procs_mumuph, ops_shapes).Weight(weight_zg).Tag("zgshower");
  pm.Push<Hist1D>(Axis(40, 70.0, 110.0, HiggsCandidate_mass, "m_{#mu#mu#gamma} [GeV]", {}), 
      Lead_SignalPhoton_hggPresel, procs_mumuph, ops_shapes).Weight(weight_zg).Tag("zgshower");
  pm.Push<Hist1D>(Axis(40, 70.0, 110.0, HiggsCandidate_mass, "m_{#mu#mu#gamma} [GeV]", {}), 
      Lead_SignalPhoton_pt<20&&Lead_SignalPhoton_hggPresel, procs_mumuph, ops_shapes).Weight(weight_zg).Tag("zgshower");

  vector<double> pt_bins = {15.0,17.5,20.0,22.5,25.0,27.5,30.0,35.0,40.0,45.0,50.0,60.0,70.0};
  vector<double> abseta_bins = {0,0.4,0.8,1.5,2.0,2.5};
  for (unsigned ipt = 0; ipt < pt_bins.size()-1; ipt++) {
    for (unsigned ieta = 0; ieta < abseta_bins.size()-1; ieta++) {
      NamedFunc pteta_bin_ee = NamedFunc(Lead_LoosePhoton_pt>pt_bins[ipt]&&Lead_LoosePhoton_pt<pt_bins[ipt+1]&&
                                         Lead_LoosePhoton_abseta>abseta_bins[ieta]&&Lead_LoosePhoton_abseta<abseta_bins[ieta+1])
                                        .Name("bin"+to_string(ipt)+"_"+to_string(ieta));
      NamedFunc pteta_bin_mumuph = NamedFunc(Lead_SignalPhoton_pt>pt_bins[ipt]&&Lead_SignalPhoton_pt<pt_bins[ipt+1]&&
                                             Lead_SignalPhoton_abseta>abseta_bins[ieta]&&Lead_SignalPhoton_abseta<abseta_bins[ieta+1])
                                            .Name("bin"+to_string(ipt)+"_"+to_string(ieta));

      pm.Push<Hist1D>(Axis(40, -0.58, 1.0, Lead_LoosePhoton_mvaID, "Photon IDMVA", {}), 
          ElectronPhoton_mass>80&&ElectronPhoton_mass<100&&Lead_LoosePhoton_hggPresel&&pteta_bin_ee, procs_ee, ops_shapes).Weight(weight_dy).Tag("zgshower");
      pm.Push<Hist1D>(Axis(40, -0.58, 1.0, Lead_SignalPhoton_mvaID, "Photon IDMVA", {}), 
          HiggsCandidate_mass>80&&HiggsCandidate_mass<100&&Lead_SignalPhoton_hggPresel&&pteta_bin_mumuph, procs_mumuph, ops_shapes).Weight(weight_zg).Tag("zgshower");
      pm.Push<Hist1D>(Axis(40, 0.5, 1.0, Lead_LoosePhoton_r9, "Photon R_{9}", {}), 
          ElectronPhoton_mass>80&&ElectronPhoton_mass<100&&Lead_LoosePhoton_hggPresel&&pteta_bin_ee, procs_ee, ops_shapes).Weight(weight_dy).Tag("zgshower");
      pm.Push<Hist1D>(Axis(40, 0.5, 1.0, Lead_SignalPhoton_r9, "Photon R_{9}", {}), 
          HiggsCandidate_mass>80&&HiggsCandidate_mass<100&&Lead_SignalPhoton_hggPresel&&pteta_bin_mumuph, procs_mumuph, ops_shapes).Weight(weight_zg).Tag("zgshower");
      pm.Push<Hist1D>(Axis(40, 0.0, 0.02, Lead_LoosePhoton_sieie, "Photon #sigma_{i#eta i#eta}", {}), 
          ElectronPhoton_mass>80&&ElectronPhoton_mass<100&&Lead_LoosePhoton_hggPresel&&pteta_bin_ee, procs_ee, ops_shapes).Weight(weight_dy).Tag("zgshower");
      pm.Push<Hist1D>(Axis(40, 0.0, 0.02, Lead_SignalPhoton_sieie, "Photon #sigma_{i#eta i#eta}", {}), 
          HiggsCandidate_mass>80&&HiggsCandidate_mass<100&&Lead_SignalPhoton_hggPresel&&pteta_bin_mumuph, procs_mumuph, ops_shapes).Weight(weight_zg).Tag("zgshower");
    }
  }

  pm.multithreaded_ = true; //no MVA
  pm.min_print_ = true;
  pm.SetLuminosityTag(lumi_string).MakePlots(1.0);

}
