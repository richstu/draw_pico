#ifndef H_NANO_FUNCTIONS
#define H_NANO_FUNCTIONS

#include <vector>

#include "Math/Vector4D.h"

#include "core/baby.hpp"
#include "core/named_func.hpp"

//Namedfuncs that replicate standard nano2pico behavior
namespace NanoFunctions {

  //H->Zg standard signal electron criteria
  extern const NamedFunc Electron_sig;

  //number of signal electrons
  extern const NamedFunc nSignalElectron;

  //signal electron pt
  extern const NamedFunc SignalElectron_pt;

  //lead signal electron pt
  extern const NamedFunc Lead_SignalElectron_pt;

  //sublead signal electron pt
  extern const NamedFunc Sublead_SignalElectron_pt;

  //H->Zg standard signal muon criteria
  extern const NamedFunc Muon_sig;

  //number of signal muons
  extern const NamedFunc nSignalMuon;

  //signal muon pt
  extern const NamedFunc SignalMuon_pt;

  //signal muon eta
  extern const NamedFunc SignalMuon_eta;

  //signal muon phi
  extern const NamedFunc SignalMuon_phi;

  //lead signal muon pt
  extern const NamedFunc Lead_SignalMuon_pt;

  //lead signal muon eta
  extern const NamedFunc Lead_SignalMuon_eta;

  //lead signal muon phi
  extern const NamedFunc Lead_SignalMuon_phi;

  //sublead signal muon pt
  extern const NamedFunc Sublead_SignalMuon_pt;

  //sublead signal muon eta
  extern const NamedFunc Sublead_SignalMuon_eta;

  //sublead signal muon phi
  extern const NamedFunc Sublead_SignalMuon_phi;

  //Minimum delta r between photon and a signal lepton
  extern const NamedFunc Photon_drmin;

  //H->Zg standard signal photon criteria
  extern const NamedFunc Photon_sig;

  //photon abs eta
  extern const NamedFunc Photon_abseta;

  //signal photon pt
  extern const NamedFunc SignalPhoton_pt;

  //signal photon eta
  extern const NamedFunc SignalPhoton_eta;

  //signal photon abs eta
  extern const NamedFunc SignalPhoton_abseta;

  //signal photon phi
  extern const NamedFunc SignalPhoton_phi;

  //number of signal photons
  extern const NamedFunc nSignalPhoton;

  //lead signal photon pt
  extern const NamedFunc Lead_SignalPhoton_pt;

  //lead signal photon eta
  extern const NamedFunc Lead_SignalPhoton_eta;

  //lead signal photon abs eta
  extern const NamedFunc Lead_SignalPhoton_abseta;

  //lead signal photon phi
  extern const NamedFunc Lead_SignalPhoton_phi;

  //lead signal photon idmva
  extern const NamedFunc Lead_SignalPhoton_mvaID;

  //number of signal leptons
  extern const NamedFunc nSignalLepton;

  //signal jet criteria
  extern const NamedFunc Jet_sig;

  //number of signal jets
  extern const NamedFunc nSignalJet;

  //signal jet eta
  extern const NamedFunc SignalJet_eta;

  //signal jet phi
  extern const NamedFunc SignalJet_phi;

  //number of deep jet/flavor medium-tagged jets
  extern const NamedFunc nJet_bdfm;

  //Z candidate mass
  extern const NamedFunc ZCandidate_mass;

  //Higgs candidate mass
  extern const NamedFunc HiggsCandidate_mass;

  //isolated dielectron triggers for run 2
  extern const NamedFunc HLT_pass_dielectron;

  //isolated dimuon triggers for run 2
  extern const NamedFunc HLT_pass_dimuon;

  //isolated dilepton triggers for run 2
  extern const NamedFunc HLT_pass_dilepton;

  //isolated single electron triggers for run 2
  extern const NamedFunc HLT_pass_singleelectron;

  //isolated single muon triggers for run 2
  extern const NamedFunc HLT_pass_singlemuon;

  //isolated single lepton triggers for run 2
  extern const NamedFunc HLT_pass_singlelepton;

  //isolated diphoton triggers for run 2
  extern const NamedFunc HLT_pass_diphoton;

  //isolated muon+photon triggers for run 2
  extern const NamedFunc HLT_pass_muonphoton;

  //standard triggers and pT cuts
  extern const NamedFunc Selection_HLT_pt;

  //baseline for H->Zgamma analysis
  extern const NamedFunc zg_baseline;

  //poor man's stitch variable
  extern const NamedFunc stitch;

  //lumi weight for Nanos
  extern const NamedFunc lumiWeight;

  //golden json loader
  class GoldenJsonLoader {
  public:
    GoldenJsonLoader();
    GoldenJsonLoader(GoldenJsonLoader &&) = default;
    GoldenJsonLoader& operator=(GoldenJsonLoader &&) = default;
    ~GoldenJsonLoader() = default;
    NamedFunc pass_json();
  private:
    std::vector<std::vector<int>> VVRunLumi;
  };

  float ConvertHZZMVA(float mva_mini);

  bool HzzId_WP2022(float pt, float etasc, float hzzmvaid);

  ROOT::Math::PtEtaPhiMVector GetZCandidateP4(const Baby &b);

}

#endif
