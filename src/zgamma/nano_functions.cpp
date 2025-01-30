#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "Math/Vector4D.h"

#include "core/baby.hpp"
#include "core/named_func.hpp"
#include "core/named_func_utilities.hpp"
#include "core/utilities.hpp"
#include "zgamma/nano_functions.hpp"
#include "zgamma/zg_utilities.hpp"

using namespace NamedFuncUtilities;

//have to have a wrapper because fabs overloading confuses MapNamedFunc
double abs_fn(double x){
  return fabs(x);
}

//Namedfuncs that replicate standard nano2pico behavior
namespace NanoFunctions {

  //H->Zg standard signal electron criteria
  const NamedFunc Electron_sig("Electron_sig",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> Electron_sig_;
    double wp[2][3] = {{-0.145237, -0.0315746, -0.032173},
                       { 0.604775,  0.628743,   0.896462}};
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      bool is_sig = true;
      float etasc = b.Electron_deltaEtaSC()->at(iel) + b.Electron_eta()->at(iel);
      int ipt(1), ieta(2);
      double mva = b.Electron_mvaFall17V2Iso()->at(iel);
      if (b.Electron_pt()->at(iel)>10) ipt = 0;
      if (fabs(etasc) < 0.8) ieta = 0;
      else if (fabs(etasc) < 1.479) ieta = 1;
      //check signal criteria
      if (b.Electron_pt()->at(iel) < 7) is_sig = false; //was 15
      else if (abs(etasc) > 2.5) is_sig = false;
      else if (abs(b.Electron_dz()->at(iel))>1.0) is_sig = false;
      else if (abs(b.Electron_dxy()->at(iel))>0.5) is_sig = false; 
      else if (mva <= wp[ipt][ieta]) is_sig = false;
      if (is_sig)
        Electron_sig_.push_back(1.0);
      else
        Electron_sig_.push_back(0.0);
    }
    return Electron_sig_;
  });

  //number of signal electrons
  const NamedFunc nSignalElectron = ReduceNamedFunc(Electron_sig, reduce_sum).Name("nSignalElectron");

  //signal electron pt
  const NamedFunc SignalElectron_pt = FilterNamedFunc("Electron_pt", Electron_sig).Name("SignalElectron_pt");

  //lead signal electron pt
  const NamedFunc Lead_SignalElectron_pt = ReduceNamedFunc(SignalElectron_pt, reduce_max).Name("Lead_SignalElectron_pt");

  //sublead signal electron pt
  const NamedFunc Sublead_SignalElectron_pt = ReduceNamedFunc(SignalElectron_pt, reduce_sublead).Name("Sublead_SignalElectron_pt");

  //H->Zg standard signal muon criteria
  const NamedFunc Muon_sig("Muon_sig",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> Muon_sig_;
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      bool is_sig = true;
      if (b.Muon_pt()->at(imu) < 5) is_sig = false; //was 10
      else if (abs(b.Muon_eta()->at(imu)) > 2.4) is_sig = false;
      else if (abs(b.Muon_dz()->at(imu))>1.0)  is_sig = false;
      else if (abs(b.Muon_dxy()->at(imu))>0.5) is_sig = false; 
      else if (!((b.Muon_looseId()->at(imu) || 
               (b.Muon_pt()->at(imu) > 200 && b.Muon_highPtId()->at(imu))) && 
               b.Muon_pfRelIso03_all()->at(imu) < 0.35 &&
               b.Muon_sip3d()->at(imu) < 4)) is_sig = false;
      if (is_sig)
        Muon_sig_.push_back(1.0);
      else
        Muon_sig_.push_back(0.0);
    }
    return Muon_sig_;
  });

  //number of signal muons
  const NamedFunc nSignalMuon = ReduceNamedFunc(Muon_sig, reduce_sum).Name("nSignalMuon");

  //signal muon pt
  const NamedFunc SignalMuon_pt = FilterNamedFunc("Muon_pt", Muon_sig).Name("SignalMuon_pt");

  //signal muon eta
  const NamedFunc SignalMuon_eta = FilterNamedFunc("Muon_eta", Muon_sig).Name("SignalMuon_eta");

  //signal muon phi
  const NamedFunc SignalMuon_phi = FilterNamedFunc("Muon_phi", Muon_sig).Name("SignalMuon_phi");

  //lead signal muon pt
  const NamedFunc Lead_SignalMuon_pt = ReduceNamedFunc(SignalMuon_pt, reduce_max).Name("Lead_SignalMuon_pt");

  //lead signal muon eta
  const NamedFunc Lead_SignalMuon_eta = MultiReduceNamedFunc(
      {SignalMuon_pt,SignalMuon_eta}, reduce_maxfirst).Name("Lead_SignalMuon_eta");

  //lead signal muon phi
  const NamedFunc Lead_SignalMuon_phi = MultiReduceNamedFunc(
      {SignalMuon_pt,SignalMuon_phi}, reduce_maxfirst).Name("Lead_SignalMuon_phi");

  //sublead signal muon pt
  const NamedFunc Sublead_SignalMuon_pt = ReduceNamedFunc(SignalMuon_pt, reduce_sublead).Name("Sublead_SignalMuon_pt");

  //sublead signal muon eta
  const NamedFunc Sublead_SignalMuon_eta = MultiReduceNamedFunc(
      {SignalMuon_pt,SignalMuon_eta}, reduce_subleadfirst).Name("Sublead_SignalMuon_eta");

  //sublead signal muon phi
  const NamedFunc Sublead_SignalMuon_phi = MultiReduceNamedFunc(
      {SignalMuon_pt,SignalMuon_phi}, reduce_subleadfirst).Name("Sublead_SignalMuon_phi");

  //Minimum delta r between photon and a signal lepton
  const NamedFunc Photon_drmin("Photon_drmin",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    std::vector<double> Photon_drmin_;
    for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
      float drmin = 999;
      for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
        if (!Electron_sig_[iel]) continue;
        float this_dr = deltaR(b.Photon_eta()->at(iph),b.Photon_phi()->at(iph),
                               b.Electron_eta()->at(iel),b.Electron_phi()->at(iel));
        if (this_dr < drmin) drmin = this_dr;
      }
      for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
        if (!Muon_sig_[imu]) continue;
        float this_dr = deltaR(b.Photon_eta()->at(iph),b.Photon_phi()->at(iph),
                               b.Muon_eta()->at(imu),b.Muon_phi()->at(imu));
        if (this_dr < drmin) drmin = this_dr;
      }
      Photon_drmin_.push_back(drmin);
    }
    return Photon_drmin_;
  });

  //H->Zg standard signal photon criteria
  const NamedFunc Photon_sig("Photon_sig",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> Photon_drmin_ = Photon_drmin.GetVector(b);
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    std::vector<double> Photon_sig_;
    for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
      double sig = 1;
      if (b.Photon_pt()->at(iph) < 15) sig = 0;
      else if (!(b.Photon_isScEtaEB()->at(iph) || b.Photon_isScEtaEE()->at(iph))) sig = 0;
      else if (!(b.Photon_mvaID_WP80()->at(iph))) sig = 0;
      //else if (abs(b.Photon_eta()->at(iph))<1.4442 && b.Photon_mvaID()->at(iph)<-0.4) sig = 0;
      //else if (abs(b.Photon_eta()->at(iph))>1.566 && b.Photon_mvaID()->at(iph)<-0.58) sig = 0;
      else if (!b.Photon_electronVeto()->at(iph)) sig = 0;
      else if (Photon_drmin_[iph]<0.3) sig = 0;
      else if (b.Photon_electronIdx()->at(iph) != -1) {
        if (Electron_sig_[b.Photon_electronIdx()->at(iph)]) sig = 0;
      }
      Photon_sig_.push_back(sig);
    }
    return Photon_sig_;
  });

  //signal photon pt
  const NamedFunc SignalPhoton_pt = FilterNamedFunc("Photon_pt",Photon_sig).Name("SignalPhoton_pt");

  //signal photon eta
  const NamedFunc SignalPhoton_eta = FilterNamedFunc("Photon_eta",Photon_sig).Name("SignalPhoton_eta");

  //photon abs eta
  const NamedFunc Photon_abseta = MapNamedFunc("Photon_eta",abs_fn).Name("Photon_abseta");

  //signal photon abs eta
  const NamedFunc SignalPhoton_abseta = FilterNamedFunc(Photon_abseta,Photon_sig).Name("SignalPhoton_abseta");

  //signal photon phi
  const NamedFunc SignalPhoton_phi = FilterNamedFunc("Photon_phi",Photon_sig).Name("SignalPhoton_phi");

  //signal photon mvaID
  const NamedFunc SignalPhoton_mvaID = FilterNamedFunc("Photon_mvaID",Photon_sig).Name("SignalPhoton_mvaID");

  //signal photon parton flavor (1= photon, 11=electron, 0=other)
  const NamedFunc SignalPhoton_genPartFlav = FilterNamedFunc("Photon_genPartFlav",Photon_sig).Name("SignalPhoton_genPartFlav");

  //number of signal photons
  const NamedFunc nSignalPhoton = ReduceNamedFunc(Photon_sig, reduce_sum).Name("nSignalPhoton");

  //lead signal photon pt
  const NamedFunc Lead_SignalPhoton_pt = ReduceNamedFunc(SignalPhoton_pt, reduce_max).Name("Lead_SignalPhoton_pt");

  //lead signal photon eta
  const NamedFunc Lead_SignalPhoton_eta = MultiReduceNamedFunc(
      {SignalPhoton_pt,SignalPhoton_eta}, reduce_maxfirst).Name("Lead_SignalPhoton_eta");

  //lead signal photon abs eta
  const NamedFunc Lead_SignalPhoton_abseta = MultiReduceNamedFunc(
      {SignalPhoton_pt,SignalPhoton_abseta}, reduce_maxfirst).Name("Lead_SignalPhoton_abseta");

  //lead signal photon phi
  const NamedFunc Lead_SignalPhoton_phi = MultiReduceNamedFunc(
      {SignalPhoton_pt,SignalPhoton_phi}, reduce_maxfirst).Name("Lead_SignalPhoton_phi");

  //lead signal photon idmva
  const NamedFunc Lead_SignalPhoton_mvaID = MultiReduceNamedFunc(
      {SignalPhoton_pt,SignalPhoton_mvaID}, reduce_maxfirst).Name("Lead_SignalPhoton_mvaID");

  //number of signal leptons
  const NamedFunc nSignalLepton = (nSignalElectron+nSignalMuon).Name("nSignalLepton");

  //signal jet criteria
  const NamedFunc Jet_sig("Jet_sig",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> Jet_sig_;
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    std::vector<double> ph_sig_ = Photon_sig.GetVector(b);
    for (unsigned ijet = 0; ijet < b.Jet_pt()->size(); ijet++) {
      double sig = 1;
      if (b.Jet_pt()->at(ijet) < 30) sig = 0;
      if (abs(b.Jet_eta()->at(ijet)) > 4.7) sig = 0;
      if (b.Jet_jetId()->at(ijet) < 1) sig = 0;
      for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
        if (Electron_sig_[iel] > 0.5) {
          if (deltaR(b.Electron_eta()->at(iel), b.Electron_phi()->at(iel), b.Jet_eta()->at(ijet), b.Jet_phi()->at(ijet)) < 0.4)
            sig = 0;
        }
      }
      for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
        if (Muon_sig_[imu] > 0.5) {
          if (deltaR(b.Muon_eta()->at(imu), b.Muon_phi()->at(imu), b.Jet_eta()->at(ijet), b.Jet_phi()->at(ijet)) < 0.4)
            sig = 0;
        }
      }
      for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
        if (ph_sig_[iph] > 0.5) {
          if (deltaR(b.Photon_eta()->at(iph), b.Photon_phi()->at(iph), b.Jet_eta()->at(ijet), b.Jet_phi()->at(ijet)) < 0.4)
            sig = 0;
        }
      }
      Jet_sig_.push_back(sig);
    }
    return Jet_sig_;
  });

  //number of signal jets
  const NamedFunc nSignalJet = ReduceNamedFunc(Jet_sig, reduce_sum).Name("nSignalJet");

  //signal jet eta
  const NamedFunc SignalJet_eta = FilterNamedFunc("Jet_eta", Jet_sig).Name("SignalJet_eta");

  //number of deep jet/flavor medium-tagged jets
  const NamedFunc nJet_bdfm("nJet_bdfm",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> Jet_sig_ = Jet_sig.GetVector(b);
    float nbdfm = 0;
    float medium_wp = 0.3093; //2016
    if (b.SampleType()==2017) medium_wp = 0.3033;
    else if (b.SampleType()==2018) medium_wp = 0.2770;
    for (unsigned ijet = 0; ijet < b.Jet_pt()->size(); ijet++) {
      if (!Jet_sig_[ijet]) continue;
      if (b.Jet_btagDeepFlavB()->at(ijet) > medium_wp) nbdfm += 1;
    }
    return nbdfm;
  });

  //Z candidate mass
  const NamedFunc ZCandidate_mass("ZCandidate_mass",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    double m_ll = -999;
    for (unsigned iel = 1; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_sig_[iel]) {
        for (unsigned iel2 = 0; iel2 < iel; iel2++) {
          if (Electron_sig_[iel2]) {
            ROOT::Math::PtEtaPhiMVector p1(b.Electron_pt()->at(iel),b.Electron_eta()->at(iel),b.Electron_phi()->at(iel),0.000511);
            ROOT::Math::PtEtaPhiMVector p2(b.Electron_pt()->at(iel2),b.Electron_eta()->at(iel2),b.Electron_phi()->at(iel2),0.000511);
            double this_mll = (p1+p2).M();
            if (fabs(this_mll-91.2)<fabs(m_ll-91.2))
              m_ll = this_mll;
          }
        }
      }
    }
    for (unsigned imu = 1; imu < b.Muon_pt()->size(); imu++) {
      if (Muon_sig_[imu]) {
        for (unsigned imu2 = 0; imu2 < imu; imu2++) {
          if (Muon_sig_[imu2]) {
            ROOT::Math::PtEtaPhiMVector p1(b.Muon_pt()->at(imu),b.Muon_eta()->at(imu),b.Muon_phi()->at(imu),0.106);
            ROOT::Math::PtEtaPhiMVector p2(b.Muon_pt()->at(imu2),b.Muon_eta()->at(imu2),b.Muon_phi()->at(imu2),0.106);
            double this_mll = (p1+p2).M();
            if (fabs(this_mll-91.2)<fabs(m_ll-91.2))
              m_ll = this_mll;
          }
        }
      }
    }
    return m_ll;
  });

  //Higgs candidate mass
  const NamedFunc HiggsCandidate_mass("HiggsCandidate_mass",[](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    std::vector<double> ph_sig_ = Photon_sig.GetVector(b);
    ROOT::Math::PtEtaPhiMVector ll_p;
    double m_ll = -999;
    double m_llg = -999;
    for (unsigned iel = 1; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_sig_[iel]) {
        for (unsigned iel2 = 0; iel2 < iel; iel2++) {
          if (Electron_sig_[iel2]) {
            ROOT::Math::PtEtaPhiMVector p1(b.Electron_pt()->at(iel),b.Electron_eta()->at(iel),b.Electron_phi()->at(iel),0.000511);
            ROOT::Math::PtEtaPhiMVector p2(b.Electron_pt()->at(iel2),b.Electron_eta()->at(iel2),b.Electron_phi()->at(iel2),0.000511);
            double this_mll = (p1+p2).M();
            if (fabs(this_mll-91.2)<fabs(m_ll-91.2)) {
              ll_p = p1+p2;
              m_ll = this_mll;
            }
          }
        }
      }
    }
    for (unsigned imu = 1; imu < b.Muon_pt()->size(); imu++) {
      if (Muon_sig_[imu]) {
        for (unsigned imu2 = 0; imu2 < imu; imu2++) {
          if (Muon_sig_[imu2]) {
            ROOT::Math::PtEtaPhiMVector p1(b.Muon_pt()->at(imu),b.Muon_eta()->at(imu),b.Muon_phi()->at(imu),0.106);
            ROOT::Math::PtEtaPhiMVector p2(b.Muon_pt()->at(imu2),b.Muon_eta()->at(imu2),b.Muon_phi()->at(imu2),0.106);
            double this_mll = (p1+p2).M();
            if (fabs(this_mll-91.2)<fabs(m_ll-91.2)) {
              ll_p = p1+p2; 
              m_ll = this_mll;
            }
          }
        }
      }
    }
    if (m_ll > 0) {
      for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
        if (ph_sig_[iph]) {
          ROOT::Math::PtEtaPhiMVector pph(b.Photon_pt()->at(iph),b.Photon_eta()->at(iph),b.Photon_phi()->at(iph),0.0);
          double this_mass = (pph+ll_p).M();
          if (fabs(this_mass-125.3)<fabs(m_llg-125.3))
            m_llg = this_mass;
        }
      }
    }
    return m_llg;
  });

  //isolated dielectron triggers for run 2
  const NamedFunc HLT_pass_dielectron("dielectron triggers",[](const Baby &b) -> NamedFunc::ScalarType{
    return b.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL()||b.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ();
  });

  //isolated dimuon triggers for run 2
  const NamedFunc HLT_pass_dimuon("dimuon triggers",[](const Baby &b) -> NamedFunc::ScalarType{
    if (abs(b.SampleType())==2016)
      return b.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ()||b.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ();
    return b.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8()||b.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8();
  });

  //isolated dilepton triggers for run 2
  const NamedFunc HLT_pass_dilepton = HLT_pass_dielectron||HLT_pass_dimuon;

  //isolated single electron triggers for run 2
  const NamedFunc HLT_pass_singleelectron("single electron triggers",[](const Baby &b) -> NamedFunc::ScalarType{
    if (abs(b.SampleType())==2016)
      return b.HLT_Ele27_WPTight_Gsf();
    if (abs(b.SampleType())==2017)
      return b.HLT_Ele35_WPTight_Gsf()||b.HLT_Ele32_WPTight_Gsf_L1DoubleEG();
    return b.HLT_Ele32_WPTight_Gsf();
  });

  //isolated single muon triggers for run 2
  const NamedFunc HLT_pass_singlemuon("single muon triggers",[](const Baby &b) -> NamedFunc::ScalarType{
    if (abs(b.SampleType())==2016)
      return b.HLT_IsoMu24();
    if (abs(b.SampleType())==2017)
      return b.HLT_IsoMu27()||b.HLT_IsoMu24();
    return b.HLT_IsoMu24();
  });

  //isolated single lepton triggers for run 2
  const NamedFunc HLT_pass_singlelepton = HLT_pass_singleelectron||HLT_pass_singlemuon;

  //isolated diphoton triggers for run 2
  const NamedFunc HLT_pass_diphoton("diphoton triggers",[](const Baby &b) -> NamedFunc::ScalarType{
    //Note: new EXO trigger for 2022??
    //same high-mass trigger for all years
    if (b.HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90())
      return 1.0;
    if (abs(b.SampleType())==2016) 
      return b.HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55() ||
             b.HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55();
    if (abs(b.SampleType())==2017)
      return b.HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55();
    return b.HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto();
  });

  //isolated muon+photon triggers for run 2
  const NamedFunc HLT_pass_muonphoton = NamedFunc("HLT_Mu17_Photon30_IsoCaloId")
      .Name("muon+photon triggers");

  //standard triggers and pT cuts
  const NamedFunc Selection_HLT_pt("Selection_HLT_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    double el_cut = 35;
    double mu_cut = 25;
    if (b.SampleTypeString()=="2016APV" || b.SampleTypeString()=="-2016APV" ||
        b.SampleTypeString()=="2016" || b.SampleTypeString()=="-2016")
      el_cut = 30;
    else if (b.SampleTypeString()=="2017" || b.SampleTypeString()=="-2017")
      mu_cut = 28;
    if (nSignalElectron.GetScalar(b)>0) {
      if (HLT_pass_singleelectron.GetScalar(b) && (Lead_SignalElectron_pt.GetScalar(b) > el_cut))
        return true;
    }
    if (nSignalMuon.GetScalar(b)>0) {
      if (HLT_pass_singlemuon.GetScalar(b) && (Lead_SignalMuon_pt.GetScalar(b) > mu_cut))
        return true;
    }
    if (nSignalElectron.GetScalar(b)>1) {
      if (HLT_pass_dielectron.GetScalar(b) && (Lead_SignalElectron_pt.GetScalar(b) > 25) 
          && Sublead_SignalElectron_pt.GetScalar(b)>15)
        return true;
    }
    if (nSignalMuon.GetScalar(b)>1) {
      if (HLT_pass_dimuon.GetScalar(b) && (Lead_SignalMuon_pt.GetScalar(b) > 20) 
          && Sublead_SignalMuon_pt.GetScalar(b)>10)
        return true;
    }
    return false;
  });

  //baseline for H->Zgamma analysis
  const NamedFunc zg_baseline = NamedFunc( 
      (nSignalElectron>=2 || nSignalMuon>=2) && nSignalPhoton>=1 &&
      Selection_HLT_pt &&
      ((Lead_SignalPhoton_pt/HiggsCandidate_mass)>=15.0/110.0) &&
      (ZCandidate_mass>80.0 && ZCandidate_mass<100.0) &&
      (HiggsCandidate_mass>100.0 && HiggsCandidate_mass<180.0) &&
      ((HiggsCandidate_mass+ZCandidate_mass)>185.0)).Name("baseline");

  //poor man's stitch variable
  const NamedFunc stitch("stitch",[](const Baby &b) -> NamedFunc::ScalarType{
    //skip data
    if (b.SampleTypeString().Contains("-")) return 1;
    //for MC, just use existence of true photon
    std::vector<double> photon_pflavor = SignalPhoton_genPartFlav.GetVector(b);
    bool is_real_photon = false;
    if (photon_pflavor.size()>0)
      if (photon_pflavor[0] == 1)
        is_real_photon = true;
    if ((b.FirstFileName().find("DYJets") != std::string::npos) && is_real_photon) return 0;
    if ((b.FirstFileName().find("TTJets") != std::string::npos) && is_real_photon) return 0;
    if ((b.FirstFileName().find("TTTo2L") != std::string::npos) && is_real_photon) return 0;
    //can't use VVG samples because ZZ only has leptonic decays
    //if ((b.FirstFileName().find("WW_Tune") == 0) && is_real_photon) return 0;
    //if ((b.FirstFileName().find("WZ_Tune") == 0) && is_real_photon) return 0;
    //if ((b.FirstFileName().find("ZZ_Tune") == 0) && is_real_photon) return 0;
    return 1;
  });

  //lumi weight for Nanos
  const NamedFunc lumiWeight("lumiWeight",[](const Baby &b) -> NamedFunc::ScalarType{
    //adding in samples as needed
    if (b.SampleTypeString()=="-2018") return 1.0;
    if (b.FirstFileName().find("ZGToLLG_01J_5f_lowMLL") != std::string::npos) 
      if (b.SampleTypeString()=="2018") return 0.0032085*b.genWeight()/fabs(b.genWeight());
    return b.genWeight();
  });

  //golden json loader
  GoldenJsonLoader::GoldenJsonLoader() {
    std::vector<std::string> golden_filenames = {
      "txt/golden/golden_Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt",
      "txt/golden/golden_Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt",
      "txt/golden/golden_Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"};
    for (std::string golden_filename : golden_filenames) {
      std::ifstream orgJSON;
      orgJSON.open(golden_filename.c_str());
      std::vector<int> VRunLumi;
      if(orgJSON.is_open()){
        char inChar;
        int inInt;
        std::string str;
        while(!orgJSON.eof()){
          char next = orgJSON.peek();
          if( next == '1' || next == '2' || next == '3' ||
              next == '4' || next == '5' || next == '6' ||
              next == '7' || next == '8' || next == '9' || 
              next == '0'){     
            orgJSON >> inInt;
            VRunLumi.push_back(inInt);        
          }
          else if(next == ' '){
            getline(orgJSON,str,' ');
          }
          else{
            orgJSON>>inChar;
          }
        }
      }//check if the file opened.
      else{
        std::cout<<"Invalid JSON File:"<<golden_filename<<"!\n";
      }
      orgJSON.close();
      if(VRunLumi.size() == 0){
        std::cout<<"No Lumiblock found in JSON file\n";
      }
      for(unsigned int i = 0; i+2 < VRunLumi.size();){
        if(VRunLumi[i] > 130000){
          std::vector<int> RunLumi;
          RunLumi.push_back(VRunLumi[i]);
          while(VRunLumi[i+1] < 130000 && i+1 < VRunLumi.size()){
            RunLumi.push_back(VRunLumi[i+1]);
            ++i;
          }
          VVRunLumi.push_back(RunLumi);
          ++i;
        }
      }
    }
  }

  NamedFunc GoldenJsonLoader::pass_json() {
    return NamedFunc("pass_json",[&](const Baby &b) -> NamedFunc::ScalarType{
      bool answer = false;
      if(b.run() < 120000){
        answer = true;
      }
      else{
        for(unsigned int i = 0; i < VVRunLumi.size();++i){
          if(b.run() == VVRunLumi[i][0]){
            for(unsigned int j = 1; j+1 < VVRunLumi[i].size();j=j+2){
              if(b.luminosityBlock() >= VVRunLumi[i][j] && b.luminosityBlock() <= VVRunLumi[i][j+1]){
                answer = true;
              }
            }
          }
        }
      }
      return answer;
    });
  }

}







