/**
 * Script to validate physics objects between pico productions
 */

#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <memory>

#include "TError.h"
#include "TLorentzVector.h"
#include "TMath.h"

#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/plot_maker.hpp"
#include "core/process.hpp"
#include "zgamma/zg_functions.hpp"
#include "zgamma/zg_utilities.hpp"
#include "zgamma/photon_weighter.hpp"

#include "core/baby.hpp"
#include "core/mva_wrapper.hpp"
#include "core/named_func.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/utilities.hpp"

using PlotOptTypes::OverflowType;
using std::invalid_argument;
using std::shared_ptr;
using std::string;
using std::set;
using std::vector;
using ZgFunctions::lead_jet_pt;
using ZgFunctions::lead_jet_eta;
using ZgFunctions::lead_jet_phi;
using ZgFunctions::lead_jet_m;
using ZgFunctions::sublead_jet_pt;
using ZgFunctions::sublead_jet_eta;
using ZgFunctions::sublead_jet_phi;
using ZgFunctions::sublead_jet_m;
using ZgFunctions::w_years;
using ZgFunctions::photon_relpterr;
using ZgFunctions::photon_isjet;
using ZgFunctions::llphoton_rel_pt;
using ZgFunctions::zg_baseline;
using ZgUtilities::SetProcessesBackground;
using ZgUtilities::ZgSampleLoader;

//alias to avoid naming confusion
const NamedFunc llphoton_coscaptheta0 = NamedFunc("llphoton_cosTheta[0]")
    .Name("llphoton_coscaptheta0");

//makes plots for electrons: pT eta phi (6 total)
void AddElectronPlots(PlotMaker& pm, const NamedFunc& selection, 
                      const vector<shared_ptr<Process>>& procs, 
                      const vector<PlotOpt>& ops,
                      const NamedFunc& weight, const string tag) {
  pm.Push<Hist1D>(
      Axis(25,0.0,100.0, "el_pt[0]", "Lead electron p_{T} [GeV]", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,0.0,100.0, "el_pt[1]", "Sublead electron p_{T} [GeV]", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,-2.5,2.5, "el_eta[0]", "Lead electron #eta", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,-2.5,2.5, "el_eta[1]", "Sublead electron #eta", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,-3.1416,3.1416, "el_phi[0]", "Lead electron #phi", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,-3.1416,3.1416, "el_phi[1]", "Sublead electron #phi", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  //pm.Push<Hist1D>(
  //    Axis(25,0.0,0.5, "el_reliso[0]", "Lead electron I_{rel}", {}), 
  //    selection, procs, ops).Weight(weight).Tag(tag);
  //pm.Push<Hist1D>(
  //    Axis(25,0.0,0.5, "el_reliso[1]", "Sublead electron I_{rel}", {}), 
  //    selection, procs, ops).Weight(weight).Tag(tag);
  //pm.Push<Hist1D>(
  //    Axis(25,-1.0,1.0, "el_dz[0]", "Lead electron d_{z} [cm]", {}), 
  //    selection, procs, ops_log).Weight(weight).Tag(tag);
  //pm.Push<Hist1D>(
  //    Axis(25,-1.0,1.0, "el_dz[1]", "Sublead electron d_{z} [cm]", {}), 
  //    selection, procs, ops_log).Weight(weight).Tag(tag);
  //pm.Push<Hist1D>(
  //    Axis(25,-0.5,0.5, "el_dxy[0]", "Lead electron d_{xy} [cm]", {}), 
  //    selection, procs, ops_log).Weight(weight).Tag(tag);
  //pm.Push<Hist1D>(
  //    Axis(25,-0.5,0.5, "el_dxy[1]", "Sublead electron d_{xy} [cm]", {}), 
  //    selection, procs, ops_log).Weight(weight).Tag(tag);
  //pm.Push<Hist2D>(
  //    Axis(25,-2.5,2.5, "el_eta[0]", "Lead electron #eta", {}), 
  //    Axis(25,-3.1416,3.1416, "el_phi[0]", "Lead electron #phi", {}), 
  //    selection, procs_mc, ops2d).Weight(weight)
  //    .Tag(tag+"_mc");
  //pm.Push<Hist2D>(
  //    Axis(25,-2.5,2.5, "el_eta[1]", "Sublead electron #eta", {}), 
  //    Axis(25,-3.1416,3.1416, "el_phi[1]", "Sublead electron #phi", {}), 
  //    selection, procs_mc, ops2d).Weight(weight)
  //    .Tag(tag+"_mc");
  //pm.Push<Hist2D>(
  //    Axis(25,-2.5,2.5, "el_eta[0]", "Lead electron #eta", {}), 
  //    Axis(25,-3.1416,3.1416, "el_phi[0]", "Lead electron #phi", {}), 
  //    selection, procs_data, ops2d).Weight(weight)
  //    .Tag(tag+"_data");
  //pm.Push<Hist2D>(
  //    Axis(25,-2.5,2.5, "el_eta[1]", "Sublead electron #eta", {}), 
  //    Axis(25,-3.1416,3.1416, "el_phi[1]", "Sublead electron #phi", {}), 
  //    selection, procs_data, ops2d).Weight(weight)
  //    .Tag(tag+"_data");
}

//makes plots for muons: pT eta phi (6 total)
void AddMuonPlots(PlotMaker& pm, const NamedFunc& selection, 
                  const vector<shared_ptr<Process>>& procs, 
                  const vector<PlotOpt>& ops,
                  const NamedFunc& weight, const string tag) {
  pm.Push<Hist1D>(
      Axis(25,0.0,100.0, "mu_pt[0]", "Lead muon p_{T} [GeV]", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,0.0,100.0, "mu_pt[1]", "Sublead muon p_{T} [GeV]", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,-2.5,2.5, "mu_eta[0]", "Lead muon #eta", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,-2.5,2.5, "mu_eta[1]", "Sublead muon #eta", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,-3.1416,3.1416, "mu_phi[0]", "Lead muon #phi", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,-3.1416,3.1416, "mu_phi[1]", "Sublead muon #phi", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  //pm.Push<Hist1D>(
  //    Axis(25,0.0,0.5, "mu_reliso[0]", "Lead muon I_{rel}", {}), 
  //    selection, procs, ops).Weight(weight).Tag(tag);
  //pm.Push<Hist1D>(
  //    Axis(25,0.0,0.5, "mu_reliso[1]", "Sublead muon I_{rel}", {}), 
  //    selection, procs, ops).Weight(weight).Tag(tag);
  //pm.Push<Hist1D>(
  //    Axis(25,-1.0,1.0, "mu_dz[0]", "Lead muon d_{z} [cm]", {}), 
  //    selection, procs, ops_log).Weight(weight).Tag(tag);
  //pm.Push<Hist1D>(
  //    Axis(25,-1.0,1.0, "mu_dz[1]", "Sublead muon d_{z} [cm]", {}), 
  //    selection, procs, ops_log).Weight(weight).Tag(tag);
  //pm.Push<Hist1D>(
  //    Axis(25,-0.5,0.5, "mu_dxy[0]", "Lead muon d_{xy} [cm]", {}), 
  //    selection, procs, ops_log).Weight(weight).Tag(tag);
  //pm.Push<Hist1D>(
  //    Axis(25,-0.5,0.5, "mu_dxy[1]", "Sublead muon d_{xy} [cm]", {}), 
  //    selection, procs, ops_log).Weight(weight).Tag(tag);
  //pm.Push<Hist2D>(
  //    Axis(25,-2.5,2.5, "mu_eta[0]", "Lead muon #eta", {}), 
  //    Axis(25,-3.1416,3.1416, "mu_phi[0]", "Lead muon #phi", {}), 
  //    selection, procs_mc, ops2d).Weight(weight)
  //    .Tag(tag+"_mc");
  //pm.Push<Hist2D>(
  //   Axis(25,-2.5,2.5, "mu_eta[1]", "Sublead muon #eta", {}), 
  //    Axis(25,-3.1416,3.1416, "mu_phi[1]", "Sublead muon #phi", {}), 
  //    selection, procs_mc, ops2d).Weight(weight)
  //    .Tag(tag+"_mc");
  //pm.Push<Hist2D>(
  //    Axis(25,-2.5,2.5, "mu_eta[0]", "Lead muon #eta", {}), 
  //    Axis(25,-3.1416,3.1416, "mu_phi[0]", "Lead muon #phi", {}), 
  //    selection, procs_data, ops2d).Weight(weight)
  //    .Tag(tag+"_data");
  //pm.Push<Hist2D>(
  //    Axis(25,-2.5,2.5, "mu_eta[1]", "Sublead muon #eta", {}), 
  //    Axis(25,-3.1416,3.1416, "mu_phi[1]", "Sublead muon #phi", {}), 
  //    selection, procs_data, ops2d).Weight(weight)
  //    .Tag(tag+"_data");

  //deprecated in newer productions
  //pm.Push<Hist1D>(
  //    Axis(25,0.0,100.0, "mu_corrected_pt[0]", 
  //         "Lead muon corrected p_{T} [GeV]", {}), 
  //    selection, procs, ops).Weight(weight).Tag(tag);
  //pm.Push<Hist1D>(
  //    Axis(25,0.0,100.0, "mu_corrected_pt[1]", 
  //         "Sublead muon corrected p_{T} [GeV]", {}), 
  //    selection, procs, ops).Weight(weight).Tag(tag);
}

//plots pT eta phi m of dilepton system (4 plots)
void AddDilepPlots(PlotMaker& pm, const NamedFunc& selection, 
                   const vector<shared_ptr<Process>>& procs, 
                   const vector<PlotOpt>& ops,
                   const NamedFunc& weight, const string tag) {
  pm.Push<Hist1D>(
      Axis(25,0.0,100.0, "ll_pt[0]", "p_{T}^{ll} [GeV]", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,-6.0,6.0, "ll_eta[0]", "#eta_{ll}", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,-3.1416,3.1416, "ll_phi[0]", "#phi_{ll}", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,80,100, "ll_m[0]", "m_{ll} [GeV]", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
}

//plots pT eta phi m of llphoton system (4 plots)
void AddDilepphotonPlots(PlotMaker& pm, const NamedFunc& selection, 
                         const vector<shared_ptr<Process>>& procs, 
                         const vector<PlotOpt>& ops,
                         const NamedFunc& weight, const string tag,
                         const bool is_sideband) {
  pm.Push<Hist1D>(
      Axis(25,0.0,100.0, "llphoton_pt[0]", "p_{T}^{ll#gamma} [GeV]", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,-6.0,6.0, "llphoton_eta[0]", "#eta_{ll#gamma}", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,-3.1416,3.1416, "llphoton_phi[0]", "#phi_{ll#gamma}", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  if (!is_sideband) {
    pm.Push<Hist1D>(
        Axis(25,80,100, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        selection, procs, ops).Weight(weight).Tag(tag);
  }
  else {
    pm.Push<Hist1D>(
        Axis(25,100,180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        selection, procs, ops).Weight(weight).Tag(tag);
  }
}

//Adds pT eta phi plots for photons (3 plots)
void AddPhotonPlots(PlotMaker& pm, const NamedFunc& selection, 
                    const vector<shared_ptr<Process>>& procs, 
                    const vector<PlotOpt>& ops,
                    const NamedFunc& weight, const string tag) {
  pm.Push<Hist1D>(
      Axis(25,0.0,80.0, "photon_pt[0]", "Photon p_{T} [GeV]", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,-2.5,2.5, "photon_eta[0]", "Photon #eta", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,-3.1416,3.1416, "photon_phi[0]", "Photon #phi", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  //pm.Push<Hist2D>(
  //    Axis(25,-2.5,2.5, "photon_eta[0]", "Photon #eta", {}), 
  //    Axis(25,-3.1416,3.1416, "photon_phi[0]", "Photon #phi", {}), 
  //    selection, procs_mc, ops2d).Weight(weight)
  //    .Tag(tag+"_mc");
  //pm.Push<Hist2D>(
  //    Axis(25,-2.5,2.5, "photon_eta[0]", "Photon #eta", {}), 
  //    Axis(25,-3.1416,3.1416, "photon_phi[0]", "Photon #phi", {}), 
  //    selection, procs_data, ops2d).Weight(weight)
  //    .Tag(tag+"_data");
}

//Adds physics object multiplicity plots (4 total): njet, nb, nlep, met
void AddMultiplicityPlots(PlotMaker& pm, const NamedFunc& selection, 
                          const vector<shared_ptr<Process>>& procs, 
                          const vector<PlotOpt>& ops,
                          const NamedFunc& weight, const string tag) {
  pm.Push<Hist1D>(
      Axis(25,0.0,100.0, "met", "p_{T}^{miss} [GeV]", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  //pm.Push<Hist1D>(
  //    Axis(25,-3.1416,3.1416, "met_phi", "p_{T}^{miss} #phi", {}), 
  //    "met>50"&&selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(5,-0.5,4.5, "nbdfm", "N_{b medium}", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(5,-0.5,4.5, "njet", "N_{jet}", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(3,1.5,4.5, "nlep", "N_{lep}", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
}

//Adds plots of variables using in ggH BDT (8 total):
// photon IDMVA, resolution, drmin, drmax, costheta, cosTheta, phi, pT/m lly
void AddGGHBDTPlots(PlotMaker& pm, const NamedFunc& selection, 
                    const vector<shared_ptr<Process>>& procs, 
                    const vector<PlotOpt>& ops,
                    const NamedFunc& weight, const string tag) {
  pm.Push<Hist1D>(
      Axis(25,0.0,1.5, llphoton_rel_pt, "ll#gamma System p_{T}/m", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,0.0,1.0, "photon_idmva[0]", "Photon IDMVA", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,0.0,0.2, photon_relpterr, "Photon #sigma_{E}/E", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,0.0,4.0, "photon_drmin[0]", "Photon #Delta R_{min}(l)", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,0.0,5.0, "photon_drmax[0]", "Photon #Delta R_{max}(l)", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,-1.0,1.0, "llphoton_costheta[0]", "ll#gamma cos(#theta)", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,-1.0,1.0, llphoton_coscaptheta0, "ll#gamma cos(#Theta)", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,-3.1416,3.1416, "llphoton_psi[0]", "ll#gamma #phi", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
}

//Adds plots of dijet variables (15 total): jet pt eta phi m, deta, dphi,
// dr(ph j1), dr(ph j2), dphi(jj,lly), zep, balance
void AddDijetPlots(PlotMaker& pm, const NamedFunc& selection, 
                    const vector<shared_ptr<Process>>& procs, 
                    const vector<PlotOpt>& ops,
                    const NamedFunc& weight, const string tag) {
  string dijet_balance_str = "|#sum_{ll,#gamma,j1,j2} #vec{p}_{T}|"
      "/(#sum_{ll,#gamma,j1,j2} |#vec{p}_{T}|)";
  dijet_balance_str = "Dijet balance (vector p_{T} sum over scalar p_{T} sum)";

  pm.Push<Hist1D>(
      Axis(25,30.0,200.0, lead_jet_pt, "Lead jet p_{T} [GeV]", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,30.0,100.0, sublead_jet_pt, "Sublead jet p_{T} [GeV]", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,-5.0,5.0, lead_jet_eta, "Lead jet #eta", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,-5.0,5.0, sublead_jet_eta, "Sublead jet #eta", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,-3.1416,3.1416, lead_jet_phi, "Lead jet #phi", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,-3.1416,3.1416, sublead_jet_phi, "Sublead jet #phi", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,0.0,50.0, lead_jet_m, "Lead jet mass [GeV]", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,0.0,30.0, sublead_jet_m, "Sublead jet mass [GeV]", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,-3.1416,3.1416, "dijet_dphi", "Dijet #Delta #phi", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,0.0,8.0, "dijet_deta", "Dijet #Delta #eta", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,0.0,6.0, "photon_jet1_dr[0]", "#Delta R(#gamma, jet 1)", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,0.0,6.0, "photon_jet2_dr[0]", "#Delta R(#gamma, jet 2)", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,-3.1416,3.1416, "llphoton_dijet_dphi[0]", 
      "#Delta #phi(ll#gamma, jj)", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,0.0,5.0, "photon_zeppenfeld[0]", 
      "|#eta_{#gamma}-(#eta_{j1}+#eta_{j2})/2|", {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
  pm.Push<Hist1D>(
      Axis(25,0.0,1.0, "llphoton_dijet_balance[0]", dijet_balance_str, {}), 
      selection, procs, ops).Weight(weight).Tag(tag);
}

float get_w_jet(float year, int type, int njet, bool ph_isjet) {
  //DY+jet
  if (type >= 6000 && type < 7000 && ph_isjet) {
    //run 2
    if (year < 2020) {
      if (njet==0) return 0.9111605383159621;
      else if (njet==1) return 1.1376003438448439; 
      else if (njet==2) return 1.713219935063525; 
      else return 3.1788706672250915;
    }
    //run 3
    else {
      if (njet==0) return 0.9728015266693347; 
      else if (njet==1) return 1.0002059092739153; 
      else if (njet==2) return 1.3489124392588017; 
      else return 2.9050241586119054;
    }
  }
  //DY+PU
  if (type >= 6000 && type < 7000 && ph_isjet) {
    //run 2
    if (year < 2020) {
      if (njet==0) return 1.0258646318505702;
      else if (njet==1) return 0.9269326556305092;
      else if (njet==2) return 0.7428198678341376;
      else return 0.5975578113064747;
    }
    //run 3
    else {
      if (njet==0) return 1.0026662278276175;
      else if (njet==1) return 0.9723535601861952;
      else if (njet==2) return 0.9755613569795534;
      else return 0.767196839119967;
    }
  }
  //ZG
  else if (type >= 17000 && type < 18000) {
    //run 2
    if (year < 2020) {
      if (njet==0) return 1.002995001426364;
      else if (njet==1) return 0.9112824642345808; 
      else if (njet==2) return 1.0720911090823475; 
      else return 1.7016023204619064;
    }
    //run 3
    else {
      if (njet==0) return 1.0265093217058043; 
      else if (njet==1) return 0.885799943406867;
      else if (njet==2) return 1.0100105581168564;
      else return 1.5602249486008501;
    }
  }
  return 1.0;
}

//Gets fake photon weight
float get_w_fakephoton(int run, bool this_photon_isjet, float ph_pt, 
                       float ph_abseta, float ph_idmva, float ph_res) {
  static const std::vector<float> fakeph_ptbins = {15.0,20.0,30.0,500.0};
  static const std::vector<float> fakeph_absetabins = {0.0,0.8,1.5,2.0,2.5};
  static const std::vector<float> fakeph_idbins = {0.14,0.5,0.8,0.9,1.0};
  static const std::vector<float> fakeph_resbins = {0.0,0.025,0.05,1.0};

  const std::vector<float> fakeph_jet_sfs_run2 = {
    0.87784536,0.86090972,0.83359961,1.00905351,
    1.20656524,1.15353438,1.28262424,1.15435585,
    0.52378337,0.87839728,0.63154299,1.13529129};
  const std::vector<float> fakeph_pu_sfs_run2 = {
    1.38054607,1.29143011,1.22872195,1.26953528,
    1.11197443,1.34174906,1.19494895,1.25048567,
    2.59609818,1.86582932,3.16422829,1.41518427};
  const std::vector<float> fakeph_jet_sfs_run3 = {
    0.73009776,1.02800363,0.76673460,1.29045474,
    1.03924426,1.12586378,0.90955539,1.10089699,
    0.74353182,0.34798856,0.79232260,0.69352779};
  const std::vector<float> fakeph_pu_sfs_run3 = {
    1.26284615,1.17956303,1.28810042,1.39744805,
    1.45429152,1.37924462,1.47740148,1.47429645,
    2.11270871,3.20168632,1.77614373,1.49915009};

  ph_idmva += 0.0;
  ph_res += 0.0;
  //static const std::vector<float> fakeph_idres_sfs_run2 = {
  //  0.86703542,1.13279172,1.03335977,
  //  0.99237755,1.18276864,1.01707066,
  //  1.18351245,1.28974280,1.14538561,
  //  1.37581143,1.39472730,0.97234205};
  //static const std::vector<float> fakeph_idres_sfs_run3 = {
  //  1.08025324,1.16182780,1.42111267,
  //  1.17597234,1.24471569,1.33403825,
  //  1.09437428,1.51553394,2.26460969,
  //  1.14807467,1.31767408,2.92946390};

  const std::vector<float> *fakeph_sfs;
  //const std::vector<float> *fakeph_idres_sfs;
  if (run == 2 && this_photon_isjet) {
    fakeph_sfs = &fakeph_jet_sfs_run2;
    //fakeph_idres_sfs = &fakeph_idres_sfs_run2;
  }
  else if (run == 2 && !this_photon_isjet) {
    fakeph_sfs = &fakeph_pu_sfs_run2;
    //fakeph_idres_sfs = &fakeph_idres_sfs_run2;
  }
  else if (run == 3 && this_photon_isjet) {
    fakeph_sfs = &fakeph_jet_sfs_run3;
    //fakeph_idres_sfs = &fakeph_idres_sfs_run3;
  }
  else if (run == 3 && !this_photon_isjet) {
    fakeph_sfs = &fakeph_pu_sfs_run3;
    //fakeph_idres_sfs = &fakeph_idres_sfs_run3;
  }
  else return 1.0;

  float sf_pteta = 1.0;
  float sf_idres = 1.0;
  for (unsigned ipt = 0; ipt < (fakeph_ptbins.size()-1); ipt++) {
    for (unsigned ieta = 0; ieta < (fakeph_absetabins.size()-1); ieta++) {
      if (ph_pt >= fakeph_ptbins[ipt] && ph_pt < fakeph_ptbins[ipt+1] 
          && ph_abseta >= fakeph_absetabins[ieta] 
          && ph_abseta < fakeph_absetabins[ieta+1]) {
        sf_pteta = fakeph_sfs->at(ipt*(fakeph_absetabins.size()-1)+ieta);
        break;
      }
    }
  }
  //for (unsigned iid = 0; iid < (fakeph_idbins.size()-1); iid++) {
  //  for (unsigned ires = 0; ires < (fakeph_resbins.size()-1); ires++) {
  //    if (ph_idmva >= fakeph_idbins[iid] && ph_idmva < fakeph_idbins[iid+1] 
  //        && ph_res >= fakeph_resbins[ires] 
  //        && ph_res < fakeph_resbins[ires+1]) {
  //      sf_idres = fakeph_idres_sfs->at(iid*(fakeph_resbins.size()-1)+ires);
  //      break;
  //    }
  //  }
  //}
  return sf_pteta*sf_idres;
}

//apply weights for photons with pT<20 in pinnacles_v0
float get_w_photon_lowpt(vector<bool> &photon_sig, vector<int> &photon_pflavor, 
                         vector<float> &photon_pt, vector<float> &photon_eta,
                         float year) {
  if (year>2022.9) return 1.0;
  float w = 1.0;
  for (unsigned iph = 0; iph < photon_pflavor.size(); iph++) {
    if (photon_pflavor[iph]==1 && photon_pt[iph]>15 && photon_pt[iph]<20) {
      float eta = photon_eta[iph];
      bool sig = photon_sig[iph];
      if (year > 2015.9 && year < 2016.1) {
        if (sig && eta > -2.5 && eta < -2.0)      w *= 1.0637626399875102;
        else if (sig && eta > -2.0 && eta < -1.5) w *= 0.980941745169525;
        else if (sig && eta > -1.5 && eta < -0.8) w *= 0.9713806263485965; 
        else if (sig && eta > -0.8 && eta < 0.0)  w *= 1.016705800510692;
        else if (sig && eta > 0.0 && eta < 0.8)   w *= 1.013255390999418;
        else if (sig && eta > 0.8 && eta < 1.5)   w *= 1.004959318390466; 
        else if (sig && eta > 1.5 && eta < 2.0)   w *= 0.9028498618621597;
        else if (sig && eta > 2.0 && eta < 2.5)   w *= 1.0078660115088238;
        else if (!sig && eta > -2.5 && eta < -2.0) w *= 0.6458359529885689;
        else if (!sig && eta > -2.0 && eta < -1.5) w *= 1.0620333583712274;
        else if (!sig && eta > -1.5 && eta < -0.8) w *= 1.1284296993607925;
        else if (!sig && eta > -0.8 && eta < 0.0)  w *= 0.8791299133785585;
        else if (!sig && eta > 0.0 && eta < 0.8)   w *= 0.9078126534514417;
        else if (!sig && eta > 0.8 && eta < 1.5)   w *= 0.9790640106166821;
        else if (!sig && eta > 1.5 && eta < 2.0)   w *= 1.3181105865876797;
        else if (!sig && eta > 2.0 && eta < 2.5)   w *= 0.9540851991993182;
      }
      else if (year > 2016.4 && year < 2016.6) {
        if (sig && eta > -2.5 && eta < -2.0)      w *= 0.9285812754096276;
        else if (sig && eta > -2.0 && eta < -1.5) w *= 1.0339806085618999;
        else if (sig && eta > -1.5 && eta < -0.8) w *= 0.9610075707366368; 
        else if (sig && eta > -0.8 && eta < 0.0)  w *= 0.9631877412261108;
        else if (sig && eta > 0.0 && eta < 0.8)   w *= 1.0046065422234682;
        else if (sig && eta > 0.8 && eta < 1.5)   w *= 1.0217023475344034;
        else if (sig && eta > 1.5 && eta < 2.0)   w *= 0.9611024969289972;
        else if (sig && eta > 2.0 && eta < 2.5)   w *= 0.9499203568038671;
        else if (!sig && eta > -2.5 && eta < -2.0) w *= 1.3690933970619903;
        else if (!sig && eta > -2.0 && eta < -1.5) w *= 0.8979370999676594;
        else if (!sig && eta > -1.5 && eta < -0.8) w *= 1.1627555236944658;
        else if (!sig && eta > -0.8 && eta < 0.0)  w *= 1.248132746431878;
        else if (!sig && eta > 0.0 && eta < 0.8)   w *= 0.9690851995285417;
        else if (!sig && eta > 0.8 && eta < 1.5)   w *= 0.9137491261701439;
        else if (!sig && eta > 1.5 && eta < 2.0)   w *= 1.1156632360862162;
        else if (!sig && eta > 2.0 && eta < 2.5)   w *= 1.2645235937513208;
      }
      else if (year > 2016.9 && year < 2017.1) {
        if (sig && eta > -2.5 && eta < -2.0)      w *= 1.0038600498776258;
        else if (sig && eta > -2.0 && eta < -1.5) w *= 0.973441200861423;
        else if (sig && eta > -1.5 && eta < -0.8) w *= 0.9566560634147776; 
        else if (sig && eta > -0.8 && eta < 0.0)  w *= 0.9949802364843384;
        else if (sig && eta > 0.0 && eta < 0.8)   w *= 0.9822339389814418;
        else if (sig && eta > 0.8 && eta < 1.5)   w *= 0.9735398543153528;
        else if (sig && eta > 1.5 && eta < 2.0)   w *= 0.9589739931516128;
        else if (sig && eta > 2.0 && eta < 2.5)   w *= 1.0588456176736205;
        else if (!sig && eta > -2.5 && eta < -2.0) w *= 0.9850895146655023;
        else if (!sig && eta > -2.0 && eta < -1.5) w *= 1.0711514790710508;
        else if (!sig && eta > -1.5 && eta < -0.8) w *= 1.1481185013060329;
        else if (!sig && eta > -0.8 && eta < 0.0)  w *= 1.028288938158254;
        else if (!sig && eta > 0.0 && eta < 0.8)   w *= 1.0985109885634243;
        else if (!sig && eta > 0.8 && eta < 1.5)   w *= 1.0891240442031525;
        else if (!sig && eta > 1.5 && eta < 2.0)   w *= 1.1097718747709093;
        else if (!sig && eta > 2.0 && eta < 2.5)   w *= 0.7665019031151598;
      }
      else if (year > 2017.9 && year < 2018.1) {
        if (sig && eta > -2.5 && eta < -2.0)      w *= 1.0026715628389615;
        else if (sig && eta > -2.0 && eta < -1.5) w *= 0.9589399238464549;
        else if (sig && eta > -1.5 && eta < -0.8) w *= 0.9599955091347154; 
        else if (sig && eta > -0.8 && eta < 0.0)  w *= 0.9771583854366053;
        else if (sig && eta > 0.0 && eta < 0.8)   w *= 0.9912104271076769;
        else if (sig && eta > 0.8 && eta < 1.5)   w *= 0.9725093784891765;
        else if (sig && eta > 1.5 && eta < 2.0)   w *= 0.9239176798268427;
        else if (sig && eta > 2.0 && eta < 2.5)   w *= 0.965953120237465;
        else if (!sig && eta > -2.5 && eta < -2.0) w *= 0.9897845251239432;
        else if (!sig && eta > -2.0 && eta < -1.5) w *= 1.105854764771766;
        else if (!sig && eta > -1.5 && eta < -0.8) w *= 1.1355540272325784;
        else if (!sig && eta > -0.8 && eta < 0.0)  w *= 1.1233277647359423;
        else if (!sig && eta > 0.0 && eta < 0.8)   w *= 1.0459648424357058;
        else if (!sig && eta > 0.8 && eta < 1.5)   w *= 1.0901855371506401;
        else if (!sig && eta > 1.5 && eta < 2.0)   w *= 1.1947301894556555;
        else if (!sig && eta > 2.0 && eta < 2.5)   w *= 1.1346617949968554;
      }
      else if (year > 2021.9 && year < 2022.1) {
        if (sig && eta > -2.5 && eta < -2.0)      w *= 0.9952717636788198;
        else if (sig && eta > -2.0 && eta < -1.5) w *= 0.9969438529681379;
        else if (sig && eta > -1.5 && eta < -0.8) w *= 0.9101314138722953; 
        else if (sig && eta > -0.8 && eta < 0.0)  w *= 0.9696422862626193;
        else if (sig && eta > 0.0 && eta < 0.8)   w *= 0.9857683103368693;
        else if (sig && eta > 0.8 && eta < 1.5)   w *= 1.0480699104733633;
        else if (sig && eta > 1.5 && eta < 2.0)   w *= 1.0242583178501214;
        else if (sig && eta > 2.0 && eta < 2.5)   w *= 1.0899023852844136;
        else if (!sig && eta > -2.5 && eta < -2.0) w *= 1.0144283882263259;
        else if (!sig && eta > -2.0 && eta < -1.5) w *= 1.0106494907269516;
        else if (!sig && eta > -1.5 && eta < -0.8) w *= 1.2020556700335518;
        else if (!sig && eta > -0.8 && eta < 0.0)  w *= 1.0949431857039922;
        else if (!sig && eta > 0.0 && eta < 0.8)   w *= 1.0471284971827695;
        else if (!sig && eta > 0.8 && eta < 1.5)   w *= 0.8866928706108905;
        else if (!sig && eta > 1.5 && eta < 2.0)   w *= 0.916163665320336;
        else if (!sig && eta > 2.0 && eta < 2.5)   w *= 0.6775232430131973;
      }
      else if (year > 2022.4 && year < 2022.6) {
        if (sig && eta > -2.5 && eta < -2.0)      w *= 1.0052678054133821;
        else if (sig && eta > -2.0 && eta < -1.5) w *= 0.962797814171453;
        else if (sig && eta > -1.5 && eta < -0.8) w *= 0.9263783918441049; 
        else if (sig && eta > -0.8 && eta < 0.0)  w *= 0.9374179718000482;
        else if (sig && eta > 0.0 && eta < 0.8)   w *= 0.9664231011166211;
        else if (sig && eta > 0.8 && eta < 1.5)   w *= 0.961355060779199;
        else if (sig && eta > 1.5 && eta < 2.0)   w *= 0.9799981025513707;
        else if (sig && eta > 2.0 && eta < 2.5)   w *= 0.9599848186660465;
        else if (!sig && eta > -2.5 && eta < -2.0) w *= 0.9825607111336123;
        else if (!sig && eta > -2.0 && eta < -1.5) w *= 1.1309203201239058;
        else if (!sig && eta > -1.5 && eta < -0.8) w *= 1.1664761752729338;
        else if (!sig && eta > -0.8 && eta < 0.0)  w *= 1.1990325160412474;
        else if (!sig && eta > 0.0 && eta < 0.8)   w *= 1.1100673514134627;
        else if (!sig && eta > 0.8 && eta < 1.5)   w *= 1.0907780198621746;
        else if (!sig && eta > 1.5 && eta < 2.0)   w *= 1.073268543888653;
        else if (!sig && eta > 2.0 && eta < 2.5)   w *= 1.1424960034725913;
      }
    }
  }
  return w;
}

//Gets llphoton pt weight
float get_w_llph_pt(int run, int type, bool this_photon_isjet, 
                    float llphoton_pt) {
  static const std::vector<float> llph_ptbins = {
      0.0,5.0,10.0,15.0,20.0,30.0,40.0,60.0,80.0,100.0,500.0};
  static const std::vector<float> llph_zg_sfs_run2 = {
      0.80095164,0.97383617,1.07866292,1.10863927,1.05103401,1.01643762,
      1.00813146,1.02230660,1.15133346,0.99455812};
  static const std::vector<float> llph_jt_sfs_run2 = {
      0.43717314,0.94298997,0.98101131,0.89560594,1.04299537,1.05185543,
      1.04641453,1.21183438,2.03966944,2.19504543};
  static const std::vector<float> llph_pu_sfs_run2 = {
      1.12644807,0.88215348,0.96358980,0.93909889,1.02689214,1.00167445,
      1.11180094,1.06075877,0.50111905,0.61013672};
  static const std::vector<float> llph_zg_sfs_run3 = {
      0.79753245,0.99954447,1.08457808,1.12089606,1.11421871,0.95435079,
      0.96804314,1.02642743,0.99670823,1.01414788};
  static const std::vector<float> llph_jt_sfs_run3 = {
      0.38981804,0.81881434,0.82474122,1.14330376,1.11759692,1.07523328,
      1.24399480,1.12016585,1.19354538,1.11502050};
  static const std::vector<float> llph_pu_sfs_run3 = {
      1.21388660,0.96262377,0.95149128,0.90487912,0.95656613,1.08625187,
      1.01390092,1.04426156,1.32007871,1.27830503};

  const std::vector<float> *llph_sfs;
  bool is_dyjet = (type >= 6000 && type < 7000 && this_photon_isjet);
  bool is_dypu = (type >= 6000 && type < 7000 && !this_photon_isjet);
  bool is_zg = (type >= 17000 && type < 18000);
  if (run == 2 && is_dyjet) llph_sfs = &llph_jt_sfs_run2;
  else if (run == 2 && is_dypu) llph_sfs = &llph_pu_sfs_run2;
  else if (run == 2 && is_zg) llph_sfs = &llph_zg_sfs_run2;
  else if (run == 3 && is_dyjet) llph_sfs = &llph_jt_sfs_run3;
  else if (run == 3 && is_dypu) llph_sfs = &llph_pu_sfs_run3;
  else if (run == 3 && is_zg) llph_sfs = &llph_zg_sfs_run3;
  else return 1.0;
  for (unsigned ipt = 0; ipt < (llph_ptbins.size()-1); ipt++) {
    if (llphoton_pt >= llph_ptbins[ipt] && llphoton_pt < llph_ptbins[ipt+1]) {
      return llph_sfs->at(ipt);
      break;
    }
  }
  return 1.0;
}

//Gets llphotonjj pt weight
//TODO change name to balance
float get_w_llpjj_pt(int run, int type, float llpjj_pt) {
  static const std::vector<float> llpjj_ptbins = {
      0.0,0.04,0.08,0.12,0.16,0.2,0.25,0.3,0.35,0.4,0.5,10.0};
  const std::vector<float> llpjj_zg_sfs_run2 = {
      0.70453953,0.73207184,0.93890318,0.88896998,0.99276156,1.32036498,
      1.02410519,1.20595101,1.89196258,1.29782599,2.42313670};
  const std::vector<float> llpjj_dy_sfs_run2 = {
      0.79156433,1.08118413,0.79126029,1.01151693,0.99202259,0.77389878,
      1.39061548,1.40040932,0.89391995,1.45405375,0.54732485};
  const std::vector<float> llpjj_zg_sfs_run3 = {
      0.70713096,0.74551468,0.85878921,0.96392069,0.91411744,1.09116120,
      1.47887854,1.79511457,1.57404356,1.63071287,1.71266180};
  const std::vector<float> llpjj_dy_sfs_run3 = {
      0.51945187,0.95900398,0.89510479,0.80162416,1.33578125,1.02527533,
      0.95155497,0.72837226,1.17599449,1.33340708,1.59637861};

  const std::vector<float> *llpjj_sfs;
  bool is_dy = (type >= 6000 && type < 7000);
  bool is_zg = (type >= 17000 && type < 18000);
  if (run == 2 && is_dy) llpjj_sfs = &llpjj_dy_sfs_run2;
  else if (run == 2 && is_zg) llpjj_sfs = &llpjj_zg_sfs_run2;
  else if (run == 3 && is_dy) llpjj_sfs = &llpjj_dy_sfs_run3;
  else if (run == 3 && is_zg) llpjj_sfs = &llpjj_zg_sfs_run3;
  else return 1.0;
  for (unsigned ipt = 0; ipt < (llpjj_ptbins.size()-1); ipt++) {
    if (llpjj_pt >= llpjj_ptbins[ipt] && llpjj_pt < llpjj_ptbins[ipt+1]) {
      return llpjj_sfs->at(ipt);
      break;
    }
  }
  return 1.0;
}

//evaluates Chebyshev polynomial (first kind)
float eval_cheby(vector<float> coefs, float x) {
  float result = 0.0;
  if (coefs.size()>0)
    result += coefs[0];
  if (coefs.size()>1)
    result += coefs[1]*x;
  if (coefs.size()>2)
    result += coefs[2]*(2.0*x*x-1.0);
  if (coefs.size()>3)
    result += coefs[3]*(4.0*x*x*x-3.0*x);
  if (coefs.size()>4)
    result += coefs[4]*(8.0*x*x*x*x-8.0*x*x+1.0);
  return result;
}

int main() {

  //---------------------------------------------------------------------------
  //                                       settings
  //---------------------------------------------------------------------------

  string production = "pinnaclesv0";
  string years = "Run3";

  //---------------------------------------------------------------------------
  //                                    initialization
  //---------------------------------------------------------------------------

  gErrorIgnoreLevel = 6000;

  //set lumi string
  string lumi_string = "138";
  set<string> year_set;
  set<string> data_dirs;
  if (years == "2016APV") {
    lumi_string = "20";
    year_set = {"2016APV"};
    data_dirs = {"/net/cms11/cms11r0/pico/NanoAODv9UCSB2/"
                 "htozgamma_pinnacles_v0/2016APV/"};
  }
  else if (years == "2016") {
    lumi_string = "17";
    year_set = {"2016"};
    data_dirs = {"/net/cms11/cms11r0/pico/NanoAODv9/"
                 "htozgamma_pinnacles_v0/2016/"};
  }
  else if (years == "2017") {
    lumi_string = "41";
    year_set = {"2017"};
    data_dirs = {"/net/cms11/cms11r0/pico/NanoAODv9/"
                 "htozgamma_pinnacles_v0/2017/"};
  }
  else if (years == "2018") {
    lumi_string = "60";
    year_set = {"2018"};
    data_dirs = {"/net/cms11/cms11r0/pico/NanoAODv9/"
                 "htozgamma_pinnacles_v0/2018/"};
  }
  else if (years == "2022") {
    lumi_string = "8";
    year_set = {"2022"};
    data_dirs = {"/net/cms11/cms11r0/pico/NanoAODv12/"
                 "htozgamma_pinnacles_v0/2022/"};
  }
  else if (years == "2022EE") {
    lumi_string = "27";
    year_set = {"2022EE"};
    data_dirs = {"/net/cms11/cms11r0/pico/NanoAODv12/"
                 "htozgamma_pinnacles_v0/2022EE/"};
  }
  else if (years == "2023") {
    lumi_string = "18";
    year_set = {"2023"};
    data_dirs = {"/net/cms11/cms11r0/pico/NanoAODv12/"
                 "htozgamma_pinnacles_v0/2023/"};
  }
  else if (years == "2023BPix") {
    lumi_string = "10";
    year_set = {"2023BPix"};
    data_dirs = {"/net/cms11/cms11r0/pico/NanoAODv12/"
                 "htozgamma_pinnacles_v0/2023BPix/"};
  }
  else if (years == "Run2") {
    lumi_string = "138";
    year_set = {"2016APV","2016","2017","2018"};
    data_dirs = {"/net/cms11/cms11r0/pico/NanoAODv9UCSB2/"
                 "htozgamma_pinnacles_v1/2016APV/",
                 "/net/cms11/cms11r0/pico/NanoAODv9/"
                 "htozgamma_pinnacles_v0/2016/",
                 "/net/cms11/cms11r0/pico/NanoAODv9/"
                 "htozgamma_pinnacles_v0/2017/",
                 "/net/cms11/cms11r0/pico/NanoAODv9/"
                 "htozgamma_pinnacles_v0/2018/"};
  }
  else if (years == "Run3") {
    lumi_string = "62";
    year_set = {"2022","2022EE","2023","2023BPix"};
    data_dirs = {"/net/cms11/cms11r0/pico/NanoAODv12/"
                 "htozgamma_pinnacles_v0/2022/",
                 "/net/cms11/cms11r0/pico/NanoAODv12/"
                 "htozgamma_pinnacles_v0/2022EE/",
                 "/net/cms11/cms11r0/pico/NanoAODv12/"
                 "htozgamma_pinnacles_v0/2023/",
                 "/net/cms11/cms11r0/pico/NanoAODv12/"
                 "htozgamma_pinnacles_v0/2023BPix/"};
  }
  else {
    throw invalid_argument("Unknown year");
  }

  //load processes
  vector<shared_ptr<Process>> procs;
  vector<shared_ptr<Process>> procs_mc;
  vector<shared_ptr<Process>> procs_data;
  //string region = "llcr";
  string region = "All";
  if (production == "pinnaclesv0") {
    procs = ZgSampleLoader().SetMacro("YEARS",year_set)
        .SetMacro("DATA_PRODUCTION_YEARS",data_dirs)
        .LoadSamples("txt/samples_zgamma.txt",region);
    procs_mc = ZgSampleLoader().SetMacro("YEARS",year_set)
        .SetMacro("DATA_PRODUCTION_YEARS",data_dirs)
        .LoadSamples("txt/samples_zgamma.txt",region,
        {Process::Type::background});
    procs_data = ZgSampleLoader().SetMacro("YEARS",year_set)
        .SetMacro("DATA_PRODUCTION_YEARS",data_dirs)
        .LoadSamples("txt/samples_zgamma.txt",region,
        {Process::Type::data});
    SetProcessesBackground(procs_data);
  }
  else if (production == "lassenv0") {
    set<string> production_names = {
        "/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/",
        "/net/cms11/cms11r0/pico/NanoAODv11p9/htozgamma_lassen_v0/"};
    procs = ZgSampleLoader().SetMacro("YEARS",year_set)
        .SetMacro("PRODUCITON_NAME",production_names)
        .LoadSamples("txt/samples_zgamma.txt",region);
    procs_mc = ZgSampleLoader().SetMacro("YEARS",year_set)
        .SetMacro("PRODUCITON_NAME",production_names)
        .LoadSamples("txt/samples_zgamma.txt",region,
        {Process::Type::background});
    procs_data = ZgSampleLoader().SetMacro("YEARS",year_set)
        .SetMacro("PRODUCITON_NAME",production_names)
        .LoadSamples("txt/samples_zgamma.txt",region,{Process::Type::data});
    SetProcessesBackground(procs_data);
  }
  else if (production == "kingscanyonv1") {
    procs = ZgSampleLoader().SetMacro("YEARS",year_set)
        .LoadSamples("txt/samples_zgamma_kingscanyonv1.txt",region);
    procs_mc = ZgSampleLoader().SetMacro("YEARS",year_set)
        .LoadSamples("txt/samples_zgamma_kingscanyonv1.txt",region,
        {Process::Type::background});
    procs_data = ZgSampleLoader().SetMacro("YEARS",year_set)
        .LoadSamples("txt/samples_zgamma_kingscanyonv1.txt",region,
        {Process::Type::data});
    SetProcessesBackground(procs_data);
  }
  else if (production == "kingscanyonv0") {
    procs = ZgSampleLoader().SetMacro("YEARS",year_set)
        .LoadSamples("txt/samples_zgamma_kingscanyonv0.txt",region);
    procs_mc = ZgSampleLoader().SetMacro("YEARS",year_set)
        .LoadSamples("txt/samples_zgamma_kingscanyonv0.txt",region,
        {Process::Type::background});
    procs_data = ZgSampleLoader().SetMacro("YEARS",year_set)
        .LoadSamples("txt/samples_zgamma_kingscanyonv0.txt",region,
        {Process::Type::data});
    SetProcessesBackground(procs_data);
  }
  else {
    throw invalid_argument("unknown production");
  }

  //load plot opts
  vector<PlotOpt> ops = {PlotOpt("txt/plot_styles.txt","LinLumiData")
                         .RatioMinimum(0.7).RatioMaximum(1.3)};
  vector<PlotOpt> ops_mllg = {PlotOpt("txt/plot_styles.txt","LinLumiData")
                              .RatioMinimum(0.7).RatioMaximum(1.3)
                              .CanvasSize(1200,600)};
  vector<PlotOpt> ops_root = {PlotOpt("txt/plot_styles.txt",
                                      "LinLumiDataRoot")};
  vector<PlotOpt> ops_log = {PlotOpt("txt/plot_styles.txt","LogLumiDataRoot")};
  vector<PlotOpt> ops2d = {PlotOpt("txt/plot_styles.txt","LinLumi2DRoot")};
  //for debugging
  vector<PlotOpt> ops_mc = {PlotOpt("txt/plot_styles.txt","LinLumi")}; 

  //---------------------------------------------------------------------------
  //                                   NamedFuncs
  //---------------------------------------------------------------------------

  NamedFunc zee_ctrlregion = NamedFunc("nel>=2&&nphoton==0&&pass"
      "&&ll_m[0]>80&&ll_m[0]<100").Name("Ztoee_CR");
  NamedFunc zmm_ctrlregion = NamedFunc("nmu>=2&&nphoton==0&&pass"
      "&&ll_m[0]>80&&ll_m[0]<100").Name("Ztomumu_CR");
  NamedFunc zeey_ctrlregion = NamedFunc("nel>=2&&nphoton==1&&pass"
      "&&llphoton_m[0]>80&&llphoton_m[0]<100").Name("Ztoeegamma_CR");
  NamedFunc zmmy_ctrlregion = NamedFunc("nmu>=2&&nphoton==1&&pass"
      "&&llphoton_m[0]>80&&llphoton_m[0]<100").Name("Ztomumugamma_CR");
  NamedFunc zlly_ctrlregion = NamedFunc("nlep>=2&&nphoton==1&&pass"
      "&&llphoton_m[0]>80&&llphoton_m[0]<100").Name("Ztoeegamma_CR");
  NamedFunc zeeyjj_ctrlregion = NamedFunc(zeey_ctrlregion&&"njet>=2")
                                .Name("Ztoeegamma_jj_CR");
  NamedFunc zmmyjj_ctrlregion = NamedFunc(zmmy_ctrlregion&&"njet>=2")
                                .Name("Ztomumugamma_jj_CR");
  NamedFunc zllyjj_ctrlregion = NamedFunc(zlly_ctrlregion&&"njet>=2")
                                .Name("Ztollgamma_jj_CR");
  NamedFunc sideband = NamedFunc(zg_baseline
      &&"(llphoton_m[0]<120||llphoton_m[0]>130)").Name("sideband");

  //in case MET looks weird
  NamedFunc pass_ra2b_filters = NamedFunc("pass_muon_jet"
      "&&pass_low_neutral_jet&&pass_htratio_dphi_tight&&pass_ecalnoisejet")
      .Name("RA2b_filters");


  //returns ll mass calculated with corrected muon pts
  const NamedFunc ll_corrected_m("ll_corrected_m",
      [](const Baby &b) -> NamedFunc::ScalarType{
    //only recalculate for Muons
    if (b.ll_lepid()->at(0)==11) return b.ll_m()->at(0);
    TLorentzVector m1, m2;
    m1.SetPtEtaPhiM(b.mu_corrected_pt()->at(b.ll_i1()->at(0)),
                    b.mu_eta()->at(b.ll_i1()->at(0)),
                    b.mu_phi()->at(b.ll_i1()->at(0)),
                    0.106);
    m2.SetPtEtaPhiM(b.mu_corrected_pt()->at(b.ll_i2()->at(0)),
                    b.mu_eta()->at(b.ll_i2()->at(0)),
                    b.mu_phi()->at(b.ll_i2()->at(0)),
                    0.106);
    return (m1+m2).M();
  });

  //returns llgamma mass calculated with corrected muon pts
  const NamedFunc llphoton_corrected_m("llphoton_corrected_m",
      [](const Baby &b) -> NamedFunc::ScalarType{
    //only recalculate for Muons
    if (b.ll_lepid()->at(0)==11) return b.llphoton_m()->at(0);
    TLorentzVector m1, m2, ph;
    m1.SetPtEtaPhiM(b.mu_corrected_pt()->at(b.ll_i1()->at(0)),
                    b.mu_eta()->at(b.ll_i1()->at(0)),
                    b.mu_phi()->at(b.ll_i1()->at(0)),
                    0.106);
    m2.SetPtEtaPhiM(b.mu_corrected_pt()->at(b.ll_i2()->at(0)),
                    b.mu_eta()->at(b.ll_i2()->at(0)),
                    b.mu_phi()->at(b.ll_i2()->at(0)),
                    0.106);
    ph.SetPtEtaPhiM(b.photon_pt()->at(0),
                    b.photon_eta()->at(0),
                    b.photon_phi()->at(0),
                    0.0);
    return (m1+m2+ph).M();
  });

  //Z pt reweighting for DY amcatnlo sample
  const NamedFunc w_z_pt("w_z_pt",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleTypeString().Contains("-")) 
      return 1.; //data
    float z_pt = 0.0;
    if (b.type() >= 6000 && b.type() < 7000 && b.nll()>0)
      z_pt = b.ll_pt()->at(0);
    else if (b.type() >= 17000 && b.type() < 18000 && b.nllphoton()>0)
      z_pt = b.llphoton_pt()->at(0);
    else
      return 1;
    if (z_pt<2.0)                      {return 0.751271;}
    else if (z_pt>2.0 && z_pt<4.0)     {return 0.793257;}
    else if (z_pt>4.0 && z_pt<6.0)     {return 0.891502;}
    else if (z_pt>6.0 && z_pt<8.0)     {return 0.997759;}
    else if (z_pt>8.0 && z_pt<10.0)    {return 1.06628;}
    else if (z_pt>10.0 && z_pt<12.0)   {return 1.09597;}
    else if (z_pt>12.0 && z_pt<14.0)   {return 1.10345;}
    else if (z_pt>14.0 && z_pt<16.0)   {return 1.09802;}
    else if (z_pt>16.0 && z_pt<18.0)   {return 1.08284;}
    else if (z_pt>18.0 && z_pt<20.0)   {return 1.0655;}
    else if (z_pt>20.0 && z_pt<25.0)   {return 1.03228;}
    else if (z_pt>25.0 && z_pt<30.0)   {return 0.987943;}
    else if (z_pt>30.0 && z_pt<35.0)   {return 0.956658;}
    else if (z_pt>35.0 && z_pt<40.0)   {return 0.944815;}
    else if (z_pt>40.0 && z_pt<50.0)   {return 0.954598;}
    else if (z_pt>50.0 && z_pt<60.0)   {return 0.971563;}
    else if (z_pt>60.0 && z_pt<80.0)   {return 0.975438;}
    else if (z_pt>80.0 && z_pt<100.0)  {return 0.974402;}
    else if (z_pt>100.0 && z_pt<200.0) {return 0.966944;}
    return 0.947066;
  });

  const vector<float> zpt_zgpt_bins = {0.0,4.0,6.0,8.0,10.0,12.0,16.0,20.0,25.0,
      30.0,35.0,40.0,60.0,80.0,100.0,150.0,9999.0};

  const vector<float> zpt_zg_sfs = {0.7602224389896548, 0.8594182659156125, 
      0.9427405172937634, 1.0078515058173994, 1.027402291435319, 
      1.1024283988487924, 1.1321218631094403, 1.075153620886966, 
      1.129844225261892, 1.0452747081887952, 1.074880559496453, 
      1.0290776535654291, 1.0505512517367457, 1.041691435347305, 
      0.9693604102862317, 0.847823729890565};
  const vector<float> zpt_dyjet_sfs = {0.6056258849612535, 0.8723458323590593, 
      1.041518303031137, 0.9559362174022588, 0.9033533240952537, 
      1.1194133551143086, 0.9705102439928782, 1.1073018406475452, 
      1.024652623300174, 0.9853077563387657, 1.0820597934763732, 
      1.014267176881946, 0.8650985095902827, 1.006883269665257, 
      1.1189048719028651, 1.3281832181175972};
  const vector<float> zpt_dypu_sfs = {0.9174545597728025, 0.75042633700729, 
      0.9979109395726279, 0.8974956610295334, 1.1310018200826795, 
      0.9302281426954941, 0.9679594310605467, 1.040383206905858, 
      1.0052032757888913, 1.0340676886040894, 0.9214644257580706, 
      1.056340889095911, 1.0423250530462163, 0.948309513739706, 
      0.7361310785401234, 0.5876774774828196};

  //Z+photon pt reweighting for DY and ZG amcatnlo samples
  const NamedFunc w_zph_pt("w_zph_pt", [zpt_zgpt_bins, zpt_zg_sfs,
      zpt_dyjet_sfs, zpt_dypu_sfs](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleTypeString().Contains("-")) 
      return 1.; //data
    float w_z = 1.0;
    float llph_pt = b.llphoton_pt()->at(0);
    bool b_photon_isjet = static_cast<bool>(photon_isjet.GetScalar(b));
    float w_jet = get_w_jet(2018, b.type(), b.njet(), b_photon_isjet);
    if (b.type() >= 6000 && b.type() < 7000) {
      //DYJets
      for (unsigned izpt = 0; izpt < (zpt_zgpt_bins.size()-1); izpt++) {
        if (zpt_zgpt_bins[izpt] <= llph_pt 
            && llph_pt < zpt_zgpt_bins[izpt+1]) {
          if (b_photon_isjet)
            w_z = zpt_dyjet_sfs[izpt];
          else
            w_z = zpt_dypu_sfs[izpt];
          break;
        }
      }
    }
    else if (b.type() >= 17000 && b.type() < 18000) {
      //ZG
      for (unsigned izpt = 0; izpt < (zpt_zgpt_bins.size()-1); izpt++) {
        if (zpt_zgpt_bins[izpt] <= llph_pt 
            && llph_pt < zpt_zgpt_bins[izpt+1]) {
          w_z = zpt_zg_sfs[izpt];
          break;
        }
      }
    }
    return w_z*w_jet;
  });

  const vector<float> cheby_coefs_zg = {0.913740696152836, 
      -0.24408929264188717, 0.13299000925940987, -0.017620931559798977, 
      0.0006806255043007578};
  const vector<float> cheby_coefs_dyjet = {0.336742109797091, 
      0.37744911284451144, 0.011178097387588092, -0.010384841937448942, 
      0.0006970229671639716};
  const vector<float> cheby_coefs_dypu = {1.4712130417026292, 
      -0.8298994909259443, 0.23902461673380465, -0.02525488379476952, 
      0.0008301694904925589};

  //Z+photon pt and njet reweighting for DY and ZG amcatnlo samples
  const NamedFunc w_isr("w_isr", [cheby_coefs_zg, cheby_coefs_dyjet, 
      cheby_coefs_dypu](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleTypeString().Contains("-")) 
      return 1.; //data
    bool this_photon_isjet = static_cast<bool>(photon_isjet.GetScalar(b));
    float w_jet = get_w_jet(2018, b.type(), b.njet(), this_photon_isjet);
    float log_ptzg_eff = log(fmin(b.llphoton_pt()->at(0), 500));
    //DY+Jet
    if (b.type() >= 6000 && b.type() < 7000 && this_photon_isjet) {
      return w_jet*eval_cheby(cheby_coefs_dyjet, log_ptzg_eff);
    }
    //DY+PU
    if (b.type() >= 6000 && b.type() < 7000 && !this_photon_isjet) {
      return w_jet*eval_cheby(cheby_coefs_dypu, log_ptzg_eff);
    }
    //ZG
    else if (b.type() >= 17000 && b.type() < 18000) {
      return w_jet*eval_cheby(cheby_coefs_zg, log_ptzg_eff);
    }
    return 1.0;
  });

  const vector<float> w_jeteta_leadjet_r2 = {1.1274049329327025, 
      -0.19597764474249468, 0.0653894519044776, -0.006472005991385048};
  const vector<float> w_jeteta_subljet_r2 = {1.183628992449538, 
      -0.2285281251867077, 0.062497885895323166, -0.005879335867725894};
  const vector<float> w_jeteta_leadjet_r3 = {1.1167523451318244, 
      -0.18307119012762088, 0.061012049687068075, -0.0060435642695586995};
  const vector<float> w_jeteta_subljet_r3 = {1.1705324056169006, 
      -0.23539659257931417, 0.07952325410568736, -0.00816763499268695};

  //Jet eta weighting
  const NamedFunc w_jeteta("w_jeteta", [w_jeteta_leadjet_r2, 
      w_jeteta_subljet_r2, w_jeteta_leadjet_r3, w_jeteta_subljet_r3]
      (const Baby &b) -> NamedFunc::ScalarType{
    if (b.njet() < 2)
      return 1.;
    if (b.SampleTypeString().Contains("-")) 
      return 1.; //data
    bool is_run3 = false;
    if (b.SampleTypeString().Contains("202")) 
      is_run3 = true;
    if (!is_run3) {
      return eval_cheby(w_jeteta_leadjet_r2, fabs(lead_jet_eta.GetScalar(b)))*
             eval_cheby(w_jeteta_subljet_r2, 
                        fabs(sublead_jet_eta.GetScalar(b)));
    }
    return eval_cheby(w_jeteta_leadjet_r3, fabs(lead_jet_eta.GetScalar(b)))*
           eval_cheby(w_jeteta_subljet_r3, fabs(sublead_jet_eta.GetScalar(b)));
  });

  const NamedFunc w_photon_lowpt("w_photon_lowpt",
      [](const Baby &b) -> NamedFunc::ScalarType{
    float year = 2023.0;
    if (b.SampleTypeString().Contains("-")) 
      return 1.; //data
    if (b.SampleTypeString() == "2016APV") 
      year = 2016.0;
    else if (b.SampleTypeString() == "2016") 
      year = 2016.5;
    else if (b.SampleTypeString() == "2017") 
      year = 2017.0;
    else if (b.SampleTypeString() == "2018") 
      year = 2018.0;
    else if (b.SampleTypeString() == "2022") 
      year = 2022.0;
    else if (b.SampleTypeString() == "2022EE") 
      year = 2022.5;
    return get_w_photon_lowpt(*b.photon_sig(), *b.photon_pflavor(), 
                              *b.photon_pt(), *b.photon_eta(), year);
  });

  const rw_mmp dnn_photon_weighter;

  const NamedFunc w_photon_shape("w_photon_shape",
      [dnn_photon_weighter](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleTypeString().Contains("-")) 
      return 1.; //data
    if (b.photon_pflavor()->at(0) != 1)
      return 1.;
    vector<float> dnn_input = {b.photon_pt()->at(0), 
        fabs(b.photon_eta()->at(0)), b.photon_idmva()->at(0),
        static_cast<float>(b.photon_energyErr()->at(0)/(b.photon_pt()->at(0)
         *TMath::CosH(b.photon_eta()->at(0))))};
    float dnn_output = dnn_photon_weighter.evaluate(dnn_input);
    return (dnn_output/(1.0-dnn_output));
  });

  //NamedFunc wrapper to get fake photon weight
  const NamedFunc w_fake("w_fake",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleTypeString().Contains("-")) 
      return 1.0; //data
    if (b.type() >= 6000 && b.type() < 7000) { //DYJets
      int run = 2;
      if (b.SampleTypeString().Contains("202")) {
        run = 3;
      }
      //DYJets
      //if (b.photon_pflavor()->at(0)==1) return 1.0;
      bool b_photon_isjet = static_cast<bool>(photon_isjet.GetScalar(b));
      float ph_pt = b.photon_pt()->at(0);
      float ph_abseta = fabs(b.photon_eta()->at(0));
      float ph_idmva = b.photon_idmva()->at(0);
      float ph_res = b.photon_energyErr()->at(0)/(b.photon_pt()->at(0));
      return get_w_fakephoton(run, b_photon_isjet, ph_pt, ph_abseta, ph_idmva,
                              ph_res);
    }
    return 1.0;
  });

  //NamedFunc wrapper to get llphoton pt weight
  const NamedFunc w_llph_pt("w_llph_pt",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleTypeString().Contains("-")) 
      return 1.0; //data
    int run = 2;
    if (b.SampleTypeString().Contains("202"))
      run = 3;
    if (!((b.type() >= 6000 && b.type() < 7000) 
          || (b.type() > 17000 && b.type() < 18000))) return 1.0;
    bool b_photon_isjet = static_cast<bool>(photon_isjet.GetScalar(b));
    return get_w_llph_pt(run, b.type(), b_photon_isjet, 
                         b.llphoton_pt()->at(0));
  });

  //NamedFunc wrapper to get llphotonjj pt weight
  const NamedFunc w_llpjj_pt("w_llpjj_pt",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleTypeString().Contains("-")) 
      return 1.0; //data
    int run = 2;
    if (b.SampleTypeString().Contains("202"))
      run = 3;
    if (!((b.type() >= 6000 && b.type() < 7000) 
          || (b.type() > 17000 && b.type() < 18000))) return 1.0;
    //float pt_x = b.llphoton_pt()->at(0)*cos(b.llphoton_phi()->at(0));
    //float pt_y = b.llphoton_pt()->at(0)*sin(b.llphoton_phi()->at(0));
    //float b_lead_jet_pt = -999;
    //float b_subl_jet_pt = -999;
    //float b_lead_jet_phi = 0.0;
    //float b_subl_jet_phi = 0.0;
    //for (unsigned ijet = 0; ijet<b.jet_pt()->size(); ijet++) {
    //  if (b.jet_isgood()->at(ijet)) {
    //    if (b.jet_pt()->at(ijet) > b_lead_jet_pt) {
    //      b_subl_jet_pt = b_lead_jet_pt;
    //      b_subl_jet_phi = b_lead_jet_phi;
    //      b_lead_jet_pt = b.jet_pt()->at(ijet);
    //      b_lead_jet_phi = b.jet_phi()->at(ijet);
    //    }
    //    else if (b.jet_pt()->at(ijet) > b_subl_jet_pt) {
    //      b_subl_jet_pt = b.jet_pt()->at(ijet);
    //      b_subl_jet_phi = b.jet_phi()->at(ijet);
    //    }
    //  }
    //}
    //pt_x += b_lead_jet_pt*cos(b_lead_jet_phi);
    //pt_y += b_lead_jet_pt*sin(b_lead_jet_phi);
    //pt_x += b_subl_jet_pt*cos(b_subl_jet_phi);
    //pt_y += b_subl_jet_pt*sin(b_subl_jet_phi);
    //float llpjj_pt =  hypot(pt_x, pt_y);
    if (b.njet() < 2) return 1.0;
    return get_w_llpjj_pt(run, b.type(), b.llphoton_dijet_balance()->at(0));
  });

  const NamedFunc w_lumi_years("w_lumi*w_years",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type()<1000) return 1.0;
    return b.w_lumi()*w_years.GetScalar(b);
  });

  const NamedFunc w_lumi_lep_years("w_lumi*w_lep*w_years",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type()<1000) return 1.0;
    return b.w_lumi()*b.w_lep()*w_years.GetScalar(b);
  });

  const NamedFunc w_lumi_lep_trig_years("w_lumi*w_lep*w_trig*w_years",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type()<1000) return 1.0;
    return b.w_lumi()*b.w_lep()*b.w_trig()*w_years.GetScalar(b);
  });

  const NamedFunc w_lumi_lep_trig_pu_years("w_lumi*w_lep*w_trig*w_pu*w_years",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type()<1000) return 1.0;
    return b.w_lumi()*b.w_lep()*b.w_trig()*b.w_pu()*w_years.GetScalar(b);
  });

  const NamedFunc w_lumi_lep_trig_pu_prefire_years(
      "w_lumi*w_lep*w_trig*w_pu*w_prefire*w_years",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type()<1000) return 1.0;
    return b.w_lumi()*b.w_lep()*b.w_trig()*b.w_pu()*b.w_prefire()
           *w_years.GetScalar(b);
  });

  const NamedFunc weight_syslep_up("weight_syslep_up",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type()<1000) return 1.0;
    return b.weight()/b.w_lep()*b.sys_lep()->at(0)*w_years.GetScalar(b);
  });

  const NamedFunc weight_syslep_dn("weight_syslep_dn",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type()<1000) return 1.0;
    return b.weight()/b.w_lep()*b.sys_lep()->at(1)*w_years.GetScalar(b);
  });

  const NamedFunc weight_systrig_up("weight_systrig_up",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type()<1000) return 1.0;
    return b.weight()/b.w_trig()*b.sys_trig()->at(0)*w_years.GetScalar(b);
  });

  const NamedFunc weight_systrig_dn("weight_systrig_dn",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type()<1000) return 1.0;
    return b.weight()/b.w_trig()*b.sys_trig()->at(1)*w_years.GetScalar(b);
  });

  const NamedFunc weight_sysprefire_up("weight_sysprefire_up",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type()<1000) return 1.0;
    return b.weight()/b.w_prefire()*b.sys_prefire()->at(0)
           *w_years.GetScalar(b);
  });

  const NamedFunc weight_sysprefire_dn("weight_sysprefire_dn",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type()<1000) return 1.0;
    return b.weight()/b.w_prefire()*b.sys_prefire()->at(1)
           *w_years.GetScalar(b);
  });
  
  //---------------------------------------------------------------------------
  //                                   plots and tables
  //---------------------------------------------------------------------------
  

  string tag = ("zgvalidate_"+production+"_"+years);

  PlotMaker pm;
  pm.multithreaded_ = true;
  pm.min_print_ = true;

  bool derive_zreweight = false;
  if (derive_zreweight) {
    vector<double> ll_reweight_bins = {0.0,2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,
       18.0,20.0,25.0,30.0,35.0,40.0,50.0,60.0,80.0,100.0,200.0,300.0};
    pm.Push<Hist1D>(
        Axis(ll_reweight_bins, "ll_pt[0]", "Z candidate p_{T} [GeV]", {}), 
        zee_ctrlregion, procs, ops_root).Weight("weight"*w_years).Tag(tag);
  }

  //do plots from Z->ll(y) control region to validate physics object modelling
  bool do_validate = false;
  //for (unsigned iwgt = 0; iwgt < 7; iwgt++) {
  if (do_validate) {
    unsigned iwgt = 0;

    //set weight
    NamedFunc weight = "weight"*w_years;
    //NamedFunc weight = "weight"*w_years*w_z_pt;
    if (iwgt==1) weight = w_lumi_years;
    if (iwgt==2) weight = w_lumi_lep_years;
    if (iwgt==3) weight = w_lumi_lep_trig_years;
    if (iwgt==4) weight = "weight"*w_years*w_z_pt;
    //if (iwgt==4) weight = w_lumi_lep_trig_pu_prefire_years;
    //if (iwgt==6) weight = w_lumi_lep_trig_pu_years;
    //if (iwgt==7) weight = weight_sysprefire_up*w_years*w_z_pt;
    //if (iwgt==8) weight = weight_sysprefire_dn*w_years*w_z_pt;
    //if (iwgt==7) weight = weight_syslep_up*w_years*w_z_pt;
    //if (iwgt==8) weight = weight_syslep_dn*w_years*w_z_pt;
    //if (iwgt==9) weight = weight_systrig_up*w_years*w_z_pt;
    //if (iwgt==10) weight = weight_systrig_dn*w_years*w_z_pt;

    //10 plots for Z->ee
    AddElectronPlots(pm, zee_ctrlregion, procs, ops, weight, tag);
    AddDilepPlots(pm, zee_ctrlregion, procs, ops, weight, tag);
    //10 plots for Z->MuMu
    AddMuonPlots(pm, zmm_ctrlregion, procs, ops, weight, tag);
    AddDilepPlots(pm, zmm_ctrlregion, procs, ops, weight, tag);
    //13+12=25 plots for Z->eey
    AddElectronPlots(pm, zeey_ctrlregion, procs, ops, weight, tag);
    AddPhotonPlots(pm, zeey_ctrlregion, procs, ops, weight, tag);
    AddDilepphotonPlots(pm, zeey_ctrlregion, procs, ops, weight, tag, false);
    AddMultiplicityPlots(pm, zeey_ctrlregion, procs, ops, weight, tag);
    AddGGHBDTPlots(pm, zeey_ctrlregion, procs, ops, weight, tag);
    //13+12=25 plots for Z->mumuy
    AddMuonPlots(pm, zmmy_ctrlregion, procs, ops, weight, tag);
    AddPhotonPlots(pm, zmmy_ctrlregion, procs, ops, weight, tag);
    AddDilepphotonPlots(pm, zmmy_ctrlregion, procs, ops, weight, tag, false);
    AddMultiplicityPlots(pm, zmmy_ctrlregion, procs, ops, weight, tag);
    AddGGHBDTPlots(pm, zmmy_ctrlregion, procs, ops, weight, tag);
    //15 plots for Z->lly + dijet
    AddDijetPlots(pm, zllyjj_ctrlregion, procs, ops, weight, tag);
  }

  //do plots from regular analysis sideband
  bool do_sideband = true;
  if (do_sideband) {
    vector<NamedFunc> weights = {"weight"*w_years*w_photon_shape*w_photon_lowpt
                                 *w_fake*w_llph_pt*w_llpjj_pt};
    //vector<NamedFunc> weights = {"weight"*w_years*w_photon_shape*w_photon_lowpt
    //                             *w_fake*w_llph_pt*w_llpjj_pt,
    //                             "weight"*w_years*w_photon_shape*w_photon_lowpt
    //                             };
    for (NamedFunc& weight : weights) {
      AddElectronPlots(pm, sideband&&"ll_lepid[0]==11", procs, ops, weight, 
                       tag);
      AddMuonPlots(pm, sideband&&"ll_lepid[0]==13", procs, ops, weight, tag);
      AddPhotonPlots(pm, sideband, procs, ops, weight, tag);
      AddDilepPlots(pm, sideband, procs, ops, weight, tag);
      AddDilepphotonPlots(pm, sideband, procs, ops, weight, tag, true);
      AddMultiplicityPlots(pm, sideband, procs, ops, weight, tag);
      AddGGHBDTPlots(pm, sideband, procs, ops, weight, tag);
      AddDijetPlots(pm, sideband&&"njet>=2", procs, ops, weight, tag);
    }
  }

  pm.SetLuminosityTag(lumi_string).MakePlots(1.0,"validation");

  return 0;
}
