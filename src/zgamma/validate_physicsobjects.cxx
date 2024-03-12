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

#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/plot_maker.hpp"
#include "core/process.hpp"
#include "zgamma/zg_functions.hpp"
#include "zgamma/zg_utilities.hpp"

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
using ZgFunctions::w_years;
using ZgUtilities::GetBackgroundProcesses;
using ZgUtilities::GetDataProcesses;
using ZgUtilities::SetProcessesBackground;
using ZgUtilities::ZgSampleLoader;

int main() {

  //------------------------------------------------------------------------------------
  //                                       settings
  //------------------------------------------------------------------------------------

  string production = "kingscanyonv1";
  string years = "2016APV";
  set<string> year_set;

  //------------------------------------------------------------------------------------
  //                                    initialization
  //------------------------------------------------------------------------------------

  gErrorIgnoreLevel = 6000;

  //set lumi string
  string lumi_string = "138";
  if (years == "2016APV") {
    lumi_string = "20";
    year_set = {"2016APV"};
  }
  else if (years == "2016") {
    lumi_string = "17";
    year_set = {"2016"};
  }
  else if (years == "2017") {
    lumi_string = "41";
    year_set = {"2017"};
  }
  else if (years == "2018") {
    lumi_string = "60";
    year_set = {"2018"};
  }
  else if (years == "2022") {
    lumi_string = "8";
    year_set = {"2022"};
  }
  else if (years == "2022EE") {
    lumi_string = "27";
    year_set = {"2022EE"};
  }
  else if (years == "2023") {
    lumi_string = "18";
    year_set = {"2023"};
  }
  else if (years == "2023BPix") {
    lumi_string = "10";
    year_set = {"2023BPix"};
  }
  else if (years == "Run2") {
    lumi_string = "138";
    year_set = {"2016APV","2016","2017","2018"};
  }
  else if (years == "Run3") {
    lumi_string = "62";
    year_set = {"2022","2022EE","2023","2023BPix"};
  }
  else {
    throw invalid_argument("Unknown year");
  }

  //load processes
  vector<shared_ptr<Process>> procs;
  vector<shared_ptr<Process>> procs_mc;
  vector<shared_ptr<Process>> procs_data;
  if (production == "kingscanyonv1") {
    procs = ZgSampleLoader().SetMacro("YEARS",year_set)
        .LoadSamples("txt/samples_zgamma.txt","llcr");
    procs_mc = GetBackgroundProcesses(
        ZgSampleLoader().SetMacro("YEARS",year_set)
        .LoadSamples("txt/samples_zgamma.txt","llcr"));
    procs_data = GetDataProcesses(
        ZgSampleLoader().SetMacro("YEARS",year_set)
        .LoadSamples("txt/samples_zgamma.txt","llcr"));
    SetProcessesBackground(procs_data);
  }
  else if (production == "kingscanyonv0") {
    procs = ZgSampleLoader().SetMacro("YEARS",year_set)
        .LoadSamples("txt/samples_zgamma_kingscanyonv0.txt","llcr");
    procs_mc = GetBackgroundProcesses(
        ZgSampleLoader().SetMacro("YEARS",year_set)
        .LoadSamples("txt/samples_zgamma_kingscanyonv0.txt","llcr"));
    procs_data = GetDataProcesses(
        ZgSampleLoader().SetMacro("YEARS",year_set)
        .LoadSamples("txt/samples_zgamma_kingscanyonv0.txt","llcr"));
    SetProcessesBackground(procs_data);
  }
  else {
    throw invalid_argument("unknown production");
  }

  //load plot opts
  vector<PlotOpt> ops = {PlotOpt("txt/plot_styles.txt","LinLumiDataRoot")};
  vector<PlotOpt> ops_log = {PlotOpt("txt/plot_styles.txt","LogLumiDataRoot")};
  vector<PlotOpt> ops2d = {PlotOpt("txt/plot_styles.txt","LinLumi2DRoot")};

  //------------------------------------------------------------------------------------
  //                                   NamedFuncs
  //------------------------------------------------------------------------------------

  NamedFunc zee_ctrlregion = NamedFunc("nel==2&&nphoton==0&&ll_m[0]>81&&ll_m[0]<101").Name("Ztoee_CR");
  NamedFunc zmm_ctrlregion = NamedFunc("nmu==2&&nphoton==0&&ll_m[0]>81&&ll_m[0]<101").Name("Ztomumu_CR");
  NamedFunc zeey_ctrlregion = NamedFunc("nel==2&&nphoton==1&&llphoton_m[0]>81&&llphoton_m[0]<101").Name("Ztoeegamma_CR");
  NamedFunc zmmy_ctrlregion = NamedFunc("nmu==2&&nphoton==1&&llphoton_m[0]>81&&llphoton_m[0]<101").Name("Ztomumugamma_CR");
  NamedFunc pass_ra2b_filters = NamedFunc(
      "pass_muon_jet&&pass_low_neutral_jet&&pass_htratio_dphi_tight&&pass_ecalnoisejet")
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
    if (b.ll_pt()->at(0)<4.0) return 0.797367794342;
    else if (b.ll_pt()->at(0)>4.0 && b.ll_pt()->at(0)<8.0) return 0.963481849777;
    else if (b.ll_pt()->at(0)>8.0 && b.ll_pt()->at(0)<12.0) return 1.10947026956;
    else if (b.ll_pt()->at(0)>12.0 && b.ll_pt()->at(0)<16.0) return 1.13303095085;
    else if (b.ll_pt()->at(0)>16.0 && b.ll_pt()->at(0)<20.0) return 1.10507933272;
    else if (b.ll_pt()->at(0)>20.0 && b.ll_pt()->at(0)<24.0) return 1.06546778102;
    else if (b.ll_pt()->at(0)>24.0 && b.ll_pt()->at(0)<28.0) return 1.02690512253;
    else if (b.ll_pt()->at(0)>28.0 && b.ll_pt()->at(0)<32.0) return 0.998144202223;
    else if (b.ll_pt()->at(0)>32.0 && b.ll_pt()->at(0)<36.0) return 0.980250142177;
    else if (b.ll_pt()->at(0)>36.0 && b.ll_pt()->at(0)<40.0) return 0.973450539142;
    else if (b.ll_pt()->at(0)>40.0 && b.ll_pt()->at(0)<44.0) return 0.97617079356;
    else if (b.ll_pt()->at(0)>44.0 && b.ll_pt()->at(0)<48.0) return 0.981375335328;
    else if (b.ll_pt()->at(0)>48.0 && b.ll_pt()->at(0)<52.0) return 0.993206582579;
    else if (b.ll_pt()->at(0)>52.0 && b.ll_pt()->at(0)<56.0) return 1.00126345214;
    else if (b.ll_pt()->at(0)>56.0 && b.ll_pt()->at(0)<60.0) return 0.992131806641;
    else if (b.ll_pt()->at(0)>60.0 && b.ll_pt()->at(0)<64.0) return 1.00412729589;
    else if (b.ll_pt()->at(0)>64.0 && b.ll_pt()->at(0)<68.0) return 0.996317799849;
    else if (b.ll_pt()->at(0)>68.0 && b.ll_pt()->at(0)<72.0) return 1.00097024039;
    else if (b.ll_pt()->at(0)>72.0 && b.ll_pt()->at(0)<76.0) return 0.991976342728;
    else if (b.ll_pt()->at(0)>76.0 && b.ll_pt()->at(0)<80.0) return 1.00189295984;
    else if (b.ll_pt()->at(0)>80.0 && b.ll_pt()->at(0)<84.0) return 0.999540995194;
    else if (b.ll_pt()->at(0)>84.0 && b.ll_pt()->at(0)<88.0) return 1.00989788201;
    else if (b.ll_pt()->at(0)>88.0 && b.ll_pt()->at(0)<92.0) return 0.998043283075;
    else if (b.ll_pt()->at(0)>92.0 && b.ll_pt()->at(0)<96.0) return 0.991865704292;
    return 0.98350435542;
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

  const NamedFunc w_lumi_lep_trig_pu_prefire_years("w_lumi*w_lep*w_trig*w_pu*w_prefire*w_years",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type()<1000) return 1.0;
    return b.w_lumi()*b.w_lep()*b.w_trig()*b.w_pu()*b.w_prefire()*w_years.GetScalar(b);
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
    return b.weight()/b.w_prefire()*b.sys_prefire()->at(0)*w_years.GetScalar(b);
  });

  const NamedFunc weight_sysprefire_dn("weight_sysprefire_dn",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type()<1000) return 1.0;
    return b.weight()/b.w_prefire()*b.sys_prefire()->at(1)*w_years.GetScalar(b);
  });
  
  //------------------------------------------------------------------------------------
  //                                   plots and tables
  //------------------------------------------------------------------------------------

  PlotMaker pm;
  pm.multithreaded_ = true;
  pm.min_print_ = true;

  for (unsigned iwgt = 0; iwgt < 9; iwgt++) {

    //set weight
    NamedFunc weight = "weight"*w_years;
    if (iwgt==1) weight = w_lumi_years;
    if (iwgt==2) weight = w_lumi_lep_years;
    if (iwgt==3) weight = w_lumi_lep_trig_years;
    if (iwgt==4) weight = w_lumi_lep_trig_pu_years;
    if (iwgt==5) weight = w_lumi_lep_trig_pu_prefire_years;
    if (iwgt==6) weight = "weight"*w_years*w_z_pt;
    if (iwgt==7) weight = weight_sysprefire_up*w_years*w_z_pt;
    if (iwgt==8) weight = weight_sysprefire_dn*w_years*w_z_pt;
    //if (iwgt==7) weight = weight_syslep_up*w_years*w_z_pt;
    //if (iwgt==8) weight = weight_syslep_dn*w_years*w_z_pt;
    //if (iwgt==9) weight = weight_systrig_up*w_years*w_z_pt;
    //if (iwgt==10) weight = weight_systrig_dn*w_years*w_z_pt;

    //Z->ee electron plots
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "el_pt[0]", "Lead electron p_{T} [GeV]", {}), 
        zee_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "el_pt[1]", "Sublead electron p_{T} [GeV]", {}), 
        zee_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "el_pt[0]", "Lead electron p_{T} [GeV]", {}), 
        zee_ctrlregion&&"trig_double_el", procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "el_pt[1]", "Sublead electron p_{T} [GeV]", {}), 
        zee_ctrlregion&&"trig_double_el", procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "el_pt[0]", "Lead electron p_{T} [GeV]", {}), 
        zee_ctrlregion&&"trig_single_el", procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "el_pt[1]", "Sublead electron p_{T} [GeV]", {}), 
        zee_ctrlregion&&"trig_single_el", procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-2.5,2.5, "el_eta[0]", "Lead electron #eta", {}), 
        zee_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-2.5,2.5, "el_eta[1]", "Sublead electron #eta", {}), 
        zee_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "el_phi[0]", "Lead electron #phi", {}), 
        zee_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "el_phi[1]", "Sublead electron #phi", {}), 
        zee_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,0.5, "el_reliso[0]", "Lead electron I_{rel}", {}), 
        zee_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,0.5, "el_reliso[1]", "Sublead electron I_{rel}", {}), 
        zee_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-1.0,1.0, "el_dz[0]", "Lead electron d_{z} [cm]", {}), 
        zee_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-1.0,1.0, "el_dz[1]", "Sublead electron d_{z} [cm]", {}), 
        zee_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-0.5,0.5, "el_dxy[0]", "Lead electron d_{xy} [cm]", {}), 
        zee_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-0.5,0.5, "el_dxy[1]", "Sublead electron d_{xy} [cm]", {}), 
        zee_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist2D>(Axis(25,-2.5,2.5, "el_eta[0]", "Lead electron #eta", {}), 
        Axis(25,-3.1416,3.1416, "el_phi[0]", "Lead electron #phi", {}), 
        zee_ctrlregion, procs_mc, ops2d).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years+"_mc");
    pm.Push<Hist2D>(Axis(25,-2.5,2.5, "el_eta[1]", "Sublead electron #eta", {}), 
        Axis(25,-3.1416,3.1416, "el_phi[1]", "Sublead electron #phi", {}), 
        zee_ctrlregion, procs_mc, ops2d).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years+"_mc");
    pm.Push<Hist2D>(Axis(25,-2.5,2.5, "el_eta[0]", "Lead electron #eta", {}), 
        Axis(25,-3.1416,3.1416, "el_phi[0]", "Lead electron #phi", {}), 
        zee_ctrlregion, procs_data, ops2d).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years+"_data");
    pm.Push<Hist2D>(Axis(25,-2.5,2.5, "el_eta[1]", "Sublead electron #eta", {}), 
        Axis(25,-3.1416,3.1416, "el_phi[1]", "Sublead electron #phi", {}), 
        zee_ctrlregion, procs_data, ops2d).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years+"_data");

    //Z->ee Z plots
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "ll_pt[0]", "Z candidate p_{T} [GeV]", {}), 
        zee_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-5.0,5.0, "ll_eta[0]", "Z candidate #eta", {}), 
        zee_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "ll_phi[0]", "Z candidate #phi", {}), 
        zee_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,71,111, "ll_m[0]", "Z candidate m [GeV]", {}), 
        "nel==2&&nphoton==0", procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);

    //Z->mumu muon plots
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "mu_pt[0]", "Lead muon p_{T} [GeV]", {}), 
        zmm_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "mu_pt[1]", "Sublead muon p_{T} [GeV]", {}), 
        zmm_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "mu_corrected_pt[0]", "Lead muon corrected p_{T} [GeV]", {}), 
        zmm_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "mu_corrected_pt[1]", "Sublead muon corrected p_{T} [GeV]", {}), 
        zmm_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "mu_corrected_pt[0]", "Lead muon corrected p_{T} [GeV]", {}), 
        zmm_ctrlregion&&"trig_single_mu", procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "mu_corrected_pt[1]", "Sublead muon corrected p_{T} [GeV]", {}), 
        zmm_ctrlregion&&"trig_single_mu", procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "mu_corrected_pt[0]", "Lead muon corrected p_{T} [GeV]", {}), 
        zmm_ctrlregion&&"trig_double_mu", procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "mu_corrected_pt[1]", "Sublead muon corrected p_{T} [GeV]", {}), 
        zmm_ctrlregion&&"trig_double_mu", procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-2.5,2.5, "mu_eta[0]", "Lead muon #eta", {}), 
        zmm_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-2.5,2.5, "mu_eta[1]", "Sublead muon #eta", {}), 
        zmm_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "mu_phi[0]", "Lead muon #phi", {}), 
        zmm_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "mu_phi[1]", "Sublead muon #phi", {}), 
        zmm_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,0.5, "mu_reliso[0]", "Lead muon I_{rel}", {}), 
        zmm_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,0.5, "mu_reliso[1]", "Sublead muon I_{rel}", {}), 
        zmm_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-1.0,1.0, "mu_dz[0]", "Lead muon d_{z} [cm]", {}), 
        zmm_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-1.0,1.0, "mu_dz[1]", "Sublead muon d_{z} [cm]", {}), 
        zmm_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-0.5,0.5, "mu_dxy[0]", "Lead muon d_{xy} [cm]", {}), 
        zmm_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-0.5,0.5, "mu_dxy[1]", "Sublead muon d_{xy} [cm]", {}), 
        zmm_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist2D>(Axis(25,-2.5,2.5, "mu_eta[0]", "Lead muon #eta", {}), 
        Axis(25,-3.1416,3.1416, "mu_phi[0]", "Lead muon #phi", {}), 
        zmm_ctrlregion, procs_mc, ops2d).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years+"_mc");
    pm.Push<Hist2D>(Axis(25,-2.5,2.5, "mu_eta[1]", "Sublead muon #eta", {}), 
        Axis(25,-3.1416,3.1416, "mu_phi[1]", "Sublead muon #phi", {}), 
        zmm_ctrlregion, procs_mc, ops2d).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years+"_mc");
    pm.Push<Hist2D>(Axis(25,-2.5,2.5, "mu_eta[0]", "Lead muon #eta", {}), 
        Axis(25,-3.1416,3.1416, "mu_phi[0]", "Lead muon #phi", {}), 
        zmm_ctrlregion, procs_data, ops2d).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years+"_data");
    pm.Push<Hist2D>(Axis(25,-2.5,2.5, "mu_eta[1]", "Sublead muon #eta", {}), 
        Axis(25,-3.1416,3.1416, "mu_phi[1]", "Sublead muon #phi", {}), 
        zmm_ctrlregion, procs_data, ops2d).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years+"_data");

    //Z->mumu Z plots
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "ll_pt[0]", "Z candidate p_{T} [GeV]", {}), 
        zmm_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-5.0,5.0, "ll_eta[0]", "Z candidate #eta", {}), 
        zmm_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "ll_phi[0]", "Z candidate #phi", {}), 
        zmm_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,71,111, "ll_m[0]", "Z candidate m [GeV]", {}), 
        "nmu==2&&nphoton==0", procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,71,111, ll_corrected_m, "Z candidate corrected m [GeV]", {}), 
        "nmu==2&&nphoton==0", procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);

    //Z->eey electron plots
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "el_pt[0]", "Lead electron p_{T} [GeV]", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "el_pt[1]", "Sublead electron p_{T} [GeV]", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-2.5,2.5, "el_eta[0]", "Lead electron #eta", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-2.5,2.5, "el_eta[1]", "Sublead electron #eta", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "el_phi[0]", "Lead electron #phi", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "el_phi[1]", "Sublead electron #phi", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,0.5, "el_reliso[0]", "Lead electron I_{rel}", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,0.5, "el_reliso[1]", "Sublead electron I_{rel}", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-1.0,1.0, "el_dz[0]", "Lead electron d_{z} [cm]", {}), 
        zeey_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-1.0,1.0, "el_dz[1]", "Sublead electron d_{z} [cm]", {}), 
        zeey_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-0.5,0.5, "el_dxy[0]", "Lead electron d_{xy} [cm]", {}), 
        zeey_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-0.5,0.5, "el_dxy[1]", "Sublead electron d_{xy} [cm]", {}), 
        zeey_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist2D>(Axis(25,-2.5,2.5, "el_eta[0]", "Lead electron #eta", {}), 
        Axis(25,-3.1416,3.1416, "el_phi[0]", "Lead electron #phi", {}), 
        zeey_ctrlregion, procs_mc, ops2d).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years+"_mc");
    pm.Push<Hist2D>(Axis(25,-2.5,2.5, "el_eta[1]", "Sublead electron #eta", {}), 
        Axis(25,-3.1416,3.1416, "el_phi[1]", "Sublead electron #phi", {}), 
        zeey_ctrlregion, procs_mc, ops2d).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years+"_mc");
    pm.Push<Hist2D>(Axis(25,-2.5,2.5, "el_eta[0]", "Lead electron #eta", {}), 
        Axis(25,-3.1416,3.1416, "el_phi[0]", "Lead electron #phi", {}), 
        zeey_ctrlregion, procs_data, ops2d).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years+"_data");
    pm.Push<Hist2D>(Axis(25,-2.5,2.5, "el_eta[1]", "Sublead electron #eta", {}), 
        Axis(25,-3.1416,3.1416, "el_phi[1]", "Sublead electron #phi", {}), 
        zeey_ctrlregion, procs_data, ops2d).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years+"_data");

    //Z->eey photon plots
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "photon_pt[0]", "Photon p_{T} [GeV]", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-2.5,2.5, "photon_eta[0]", "Photon #eta", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "photon_phi[0]", "Photon #phi", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-1.0,1.0, "photon_idmva[0]", "Photon IDMVA", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist2D>(Axis(25,-2.5,2.5, "photon_eta[0]", "Photon #eta", {}), 
        Axis(25,-3.1416,3.1416, "photon_phi[0]", "Photon #phi", {}), 
        zeey_ctrlregion, procs_mc, ops2d).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years+"_mc");
    pm.Push<Hist2D>(Axis(25,-2.5,2.5, "photon_eta[0]", "Photon #eta", {}), 
        Axis(25,-3.1416,3.1416, "photon_phi[0]", "Photon #phi", {}), 
        zeey_ctrlregion, procs_data, ops2d).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years+"_mc");

    //Z->eey jet/met/Z plots
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "jet_pt", "Jet p_{T} [GeV]", {}), 
        zeey_ctrlregion&&"jet_isgood", procs, ops).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-5.0,5.0, "jet_eta", "Jet #eta", {}), 
        zeey_ctrlregion&&"jet_isgood", procs, ops).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "jet_phi", "Jet #phi", {}), 
        zeey_ctrlregion&&"njet>0", procs, ops).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist2D>(Axis(25,-5.0,5.0, "jet_eta", "Jet #eta", {}), 
        Axis(25,-3.1416,3.1416, "jet_phi", "Jet #phi", {}), 
        zeey_ctrlregion&&"jet_isgood", procs_mc, ops2d).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years+"_mc");
    pm.Push<Hist2D>(Axis(25,-5.0,5.0, "jet_eta", "Jet #eta", {}), 
        Axis(25,-3.1416,3.1416, "jet_phi", "Jet #phi", {}), 
        zeey_ctrlregion&&"jet_isgood", procs_data, ops2d).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years+"_mc");
    pm.Push<Hist1D>(Axis(5,-0.5,4.5, "njet", "N_{jet}", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(5,-0.5,4.5, "nbdfl", "N_{b loose}", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(5,-0.5,4.5, "nbdfm", "N_{b medium}", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(5,-0.5,4.5, "nbdft", "N_{b tight}", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "met", "p_{T}^{miss} [GeV]", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "met", "p_{T}^{miss} [GeV]", {}), 
        zeey_ctrlregion&&"pass", procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "met", "p_{T}^{miss} [GeV]", {}), 
        zeey_ctrlregion&&"pass"&&pass_ra2b_filters, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "met_phi", "p_{T}^{miss} #phi", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "met_phi", "p_{T}^{miss} #phi", {}), 
        zeey_ctrlregion&&"pass", procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "met_phi", "p_{T}^{miss} #phi", {}), 
        zeey_ctrlregion&&"pass"&&pass_ra2b_filters, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "llphoton_pt[0]", "Z candidate p_{T} [GeV]", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-5.0,5.0, "llphoton_eta[0]", "Z candidate #eta", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "llphoton_phi[0]", "Z candidate #phi", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,71,111, "llphoton_m[0]", "Z candidate m [GeV]", {}), 
        "nel==2&&nphoton==1", procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);

    //Z->mumuy muon plots
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "mu_pt[0]", "Lead muon p_{T} [GeV]", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "mu_pt[1]", "Sublead muon p_{T} [GeV]", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "mu_corrected_pt[0]", "Lead muon corrected p_{T} [GeV]", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "mu_corrected_pt[1]", "Sublead muon corrected p_{T} [GeV]", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-2.5,2.5, "mu_eta[0]", "Lead muon #eta", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-2.5,2.5, "mu_eta[1]", "Sublead muon #eta", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "mu_phi[0]", "Lead muon #phi", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "mu_phi[1]", "Sublead muon #phi", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,0.5, "mu_reliso[0]", "Lead muon I_{rel}", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,0.5, "mu_reliso[1]", "Sublead muon I_{rel}", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-1.0,1.0, "mu_dz[0]", "Lead muon d_{z} [cm]", {}), 
        zmmy_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-1.0,1.0, "mu_dz[1]", "Sublead muon d_{z} [cm]", {}), 
        zmmy_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-0.5,0.5, "mu_dxy[0]", "Lead muon d_{xy} [cm]", {}), 
        zmmy_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-0.5,0.5, "mu_dxy[1]", "Sublead muon d_{xy} [cm]", {}), 
        zmmy_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist2D>(Axis(25,-2.5,2.5, "mu_eta[0]", "Lead muon #eta", {}), 
        Axis(25,-3.1416,3.1416, "mu_phi[0]", "Lead muon #phi", {}), 
        zmmy_ctrlregion, procs_mc, ops2d).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years+"_mc");
    pm.Push<Hist2D>(Axis(25,-2.5,2.5, "mu_eta[1]", "Sublead muon #eta", {}), 
        Axis(25,-3.1416,3.1416, "mu_phi[1]", "Sublead muon #phi", {}), 
        zmmy_ctrlregion, procs_mc, ops2d).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years+"_mc");
    pm.Push<Hist2D>(Axis(25,-2.5,2.5, "mu_eta[0]", "Lead muon #eta", {}), 
        Axis(25,-3.1416,3.1416, "mu_phi[0]", "Lead muon #phi", {}), 
        zmmy_ctrlregion, procs_data, ops2d).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years+"_data");
    pm.Push<Hist2D>(Axis(25,-2.5,2.5, "mu_eta[1]", "Sublead muon #eta", {}), 
        Axis(25,-3.1416,3.1416, "mu_phi[1]", "Sublead muon #phi", {}), 
        zmmy_ctrlregion, procs_data, ops2d).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years+"_data");

    //Z->mumuy photon plots
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "photon_pt[0]", "Photon p_{T} [GeV]", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-2.5,2.5, "photon_eta[0]", "Photon #eta", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "photon_phi[0]", "Photon #phi", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-1.0,1.0, "photon_idmva[0]", "Photon IDMVA", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist2D>(Axis(25,-2.5,2.5, "photon_eta[0]", "Photon #eta", {}), 
        Axis(25,-3.1416,3.1416, "photon_phi[0]", "Photon #phi", {}), 
        zmmy_ctrlregion, procs_mc, ops2d).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years+"_mc");
    pm.Push<Hist2D>(Axis(25,-2.5,2.5, "photon_eta[0]", "Photon #eta", {}), 
        Axis(25,-3.1416,3.1416, "photon_phi[0]", "Photon #phi", {}), 
        zmmy_ctrlregion, procs_data, ops2d).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years+"_mc");

    //Z->mumuy jet/met/Z plots
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "jet_pt", "Jet p_{T} [GeV]", {}), 
        zmmy_ctrlregion&&"jet_isgood", procs, ops).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-5.0,5.0, "jet_eta", "Jet #eta", {}), 
        zmmy_ctrlregion&&"jet_isgood", procs, ops).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "jet_phi", "Jet #phi", {}), 
        zmmy_ctrlregion&&"njet>0", procs, ops).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist2D>(Axis(25,-5.0,5.0, "jet_eta", "Jet #eta", {}), 
        Axis(25,-3.1416,3.1416, "jet_phi", "Jet #phi", {}), 
        zmmy_ctrlregion&&"jet_isgood", procs_mc, ops2d).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years+"_mc");
    pm.Push<Hist2D>(Axis(25,-5.0,5.0, "jet_eta", "Jet #eta", {}), 
        Axis(25,-3.1416,3.1416, "jet_phi", "Jet #phi", {}), 
        zmmy_ctrlregion&&"jet_isgood", procs_data, ops2d).Weight(weight)
        .Tag("zgvalidate_"+production+"_"+years+"_mc");
    pm.Push<Hist1D>(Axis(5,-0.5,4.5, "njet", "N_{jet}", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(5,-0.5,4.5, "nbdfl", "N_{b loose}", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(5,-0.5,4.5, "nbdfm", "N_{b medium}", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(5,-0.5,4.5, "nbdft", "N_{b tight}", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "met", "p_{T}^{miss} [GeV]", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "met", "p_{T}^{miss} [GeV]", {}), 
        zmmy_ctrlregion&&"pass", procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "met", "p_{T}^{miss} [GeV]", {}), 
        zmmy_ctrlregion&&"pass"&&pass_ra2b_filters, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "met_phi", "p_{T}^{miss} #phi", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "met_phi", "p_{T}^{miss} #phi", {}), 
        zmmy_ctrlregion&&"pass", procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "met_phi", "p_{T}^{miss} #phi", {}), 
        zmmy_ctrlregion&&"pass"&&pass_ra2b_filters, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "llphoton_pt[0]", "Z candidate p_{T} [GeV]", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-5.0,5.0, "llphoton_eta[0]", "Z candidate #eta", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "llphoton_phi[0]", "Z candidate #phi", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,71,111, "llphoton_m[0]", "Z candidate m [GeV]", {}), 
        "nmu==2&&nphoton==1", procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,71,111, llphoton_corrected_m, "Z candidate corrected m [GeV]", {}), 
        "nmu==2&&nphoton==1", procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
  }

  pm.SetLuminosityTag(lumi_string).MakePlots(1.0,"validation");

  return 0;
}
