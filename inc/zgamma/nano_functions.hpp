#ifndef H_NANO_FUNCTIONS
#define H_NANO_FUNCTIONS

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

  //lead signal muon pt
  extern const NamedFunc Lead_SignalMuon_pt;

  //sublead signal muon pt
  extern const NamedFunc Sublead_SignalMuon_pt;

  //Minimum delta r between photon and a signal lepton
  extern const NamedFunc Photon_drmin;

  //H->Zg standard signal photon criteria
  extern const NamedFunc Photon_sig;

  //signal photon pt
  extern const NamedFunc SignalPhoton_pt;

  //number of signal photons
  extern const NamedFunc nSignalPhoton;

  //lead signal photon pt
  extern const NamedFunc Lead_SignalPhoton_pt;

  //lead signal photon idmva
  extern const NamedFunc Lead_SignalPhoton_mvaID;

  //number of signal leptons
  extern const NamedFunc nSignalLepton;

  //signal jet criteria
  extern const NamedFunc Jet_sig;

  //number of signal jets
  extern const NamedFunc nSignalJet;

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

  //baseline for H->Zgamma analysis
  extern const NamedFunc zg_baseline;

  //poor man's stitch variable
  extern const NamedFunc stitch;

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

}

#endif
