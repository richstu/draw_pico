/**
 * script to make plots related to problematic photons in APV dataset
 */

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "TError.h"
#include "TColor.h"
#include "Math/Vector4D.h"

#include "core/baby.hpp"
#include "core/baby_nano.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/named_func_utilities.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/table_row.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/utilities.hpp"
#include "zgamma/nano_functions.hpp"
#include "zgamma/zg_functions.hpp"

using std::shared_ptr;
using std::string;
using std::vector;

using NamedFuncUtilities::FilterNamedFunc;
using NamedFuncUtilities::MapNamedFunc;
using NamedFuncUtilities::ReduceNamedFunc;
using NamedFuncUtilities::MultiReduceNamedFunc;
using NamedFuncUtilities::reduce_max;
using NamedFuncUtilities::reduce_maxfirst;
using NamedFuncUtilities::reduce_sublead;
using NamedFuncUtilities::reduce_subleadfirst;
using NamedFuncUtilities::reduce_sum;

using ZgFunctions::w_years;

using NanoFunctions::Lead_SignalMuon_pt;
using NanoFunctions::Muon_sig;
using NanoFunctions::nSignalMuon;


double dabs(double value) {
  //TODO add ascii art of Usain Bolt
  if (value < 0.0) return -1.0*value;
  return value;
}

int main() {

  //--------------------------------------------------------------------------
  //                                    initialization
  //--------------------------------------------------------------------------

  //setup
  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors_zgamma_official.txt", "default");
  
  const string year = "2016";
  string lumi_tag = "20";
  if (year == "2016")
    lumi_tag = "17";

  NanoFunctions::GoldenJsonLoader golden;

  vector<shared_ptr<Process>> procs;
  procs.push_back(Process::MakeShared<Baby_nano>("Data",
                       Process::Type::data, 1,
                       {"/net/cms11/cms11r0/pico/NanoAODv9/nano/"+year+"/data/*SingleMuon*"},
                       golden.pass_json()));
  procs.push_back(Process::MakeShared<Baby_nano>("DY+fake photon",
                       Process::Type::background, colors("dyjets"),
                       {"/net/cms11/cms11r0/pico/NanoAODv9/nano/"+year+"/mc/*DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX*"},
                       "1"));
  procs.back()->SetLineColor(colors("dyjets"));
  procs.push_back(Process::MakeShared<Baby_nano>("ll#gamma",
                       Process::Type::background, colors("zgtollg"),
                       {"/net/cms11/cms11r0/pico/NanoAODv9/nano/"+year+"/mc/*ZGToLLG_01J_5f_lowMLL*"},
                       "1"));
  procs.back()->SetLineColor(colors("zgtollg"));

  vector<shared_ptr<Process>> procs_data;
  procs_data.push_back(Process::MakeShared<Baby_nano>("Data",
                       Process::Type::data, 1,
                       {"/net/cms11/cms11r0/pico/NanoAODv9/nano/"+year+"/data/*SingleMuon*"},
                       golden.pass_json()));

  std::vector<PlotOpt> ops = {PlotOpt("txt/plot_styles.txt","LinLumiData")}; 
  std::vector<PlotOpt> ops_data = {PlotOpt("txt/plot_styles.txt","LinLumi")}; 

  //------------------------------------------------------------------------------------
  //                                   NamedFuncs
  //------------------------------------------------------------------------------------
  
  //Minimum delta r between photon and a signal muon
  const NamedFunc Photon_drmin("Photon_drmin",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    std::vector<double> Photon_drmin_;
    for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
      float drmin = 999;
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

  //this differs from the version in nano_functions because the pT threshold 
  //is higher and electrons are not considered
  const NamedFunc Photon_sig("Photon_sig",[Photon_drmin](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> Photon_drmin_ = Photon_drmin.GetVector(b);
    std::vector<double> Photon_sig_;
    for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
      double sig = 1;
      if (b.Photon_pt()->at(iph) < 20) sig = 0;
      else if (!(b.Photon_isScEtaEB()->at(iph) || b.Photon_isScEtaEE()->at(iph))) sig = 0;
      else if (!(b.Photon_mvaID_WP80()->at(iph))) sig = 0;
      else if (!b.Photon_electronVeto()->at(iph)) sig = 0;
      else if (Photon_drmin_[iph]<0.3) sig = 0;
      Photon_sig_.push_back(sig);
    }
    return Photon_sig_;
  });

  const NamedFunc Photon_sig_nopt("Photon_sig_nopt",[Photon_drmin](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> Photon_drmin_ = Photon_drmin.GetVector(b);
    std::vector<double> Photon_sig_;
    for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
      double sig = 1;
      if (!(b.Photon_isScEtaEB()->at(iph) || b.Photon_isScEtaEE()->at(iph))) sig = 0;
      else if (!(b.Photon_mvaID_WP80()->at(iph))) sig = 0;
      else if (!b.Photon_electronVeto()->at(iph)) sig = 0;
      else if (Photon_drmin_[iph]<0.3) sig = 0;
      Photon_sig_.push_back(sig);
    }
    return Photon_sig_;
  });

  //number of signal photons
  const NamedFunc nSignalPhoton = ReduceNamedFunc(
      Photon_sig, reduce_sum).Name("nSignalPhoton");

  //number of any pt signal photons
  const NamedFunc nNoptSignalPhoton = ReduceNamedFunc(
      Photon_sig_nopt, reduce_sum).Name("nNoptSignalPhoton");

  //signal photon parton flavor (1= photon, 11=electron, 0=other)
  const NamedFunc SignalPhoton_genPartFlav = FilterNamedFunc(
      "Photon_genPartFlav",Photon_sig).Name("SignalPhoton_genPartFlav");

  //signal photon pt
  const NamedFunc SignalPhoton_pt = FilterNamedFunc(
      "Photon_pt",Photon_sig).Name("SignalPhoton_pt");

  //signal photon eta
  const NamedFunc SignalPhoton_eta = FilterNamedFunc(
      "Photon_eta",Photon_sig).Name("SignalPhoton_eta");

  //signal photon phi
  const NamedFunc SignalPhoton_phi = FilterNamedFunc(
      "Photon_phi",Photon_sig).Name("SignalPhoton_phi");

  //signal photon r9
  const NamedFunc SignalPhoton_r9 = FilterNamedFunc(
      "Photon_r9",Photon_sig).Name("SignalPhoton_r9");

  //signal photon sieie
  const NamedFunc SignalPhoton_sieie = FilterNamedFunc(
      "Photon_sieie",Photon_sig).Name("SignalPhoton_sieie");

  //signal photon pfRelIso03_all
  const NamedFunc SignalPhoton_pfRelIso03_all = FilterNamedFunc(
      "Photon_pfRelIso03_all",Photon_sig).Name("SignalPhoton_pfRelIso03_all");

  //signal photon hoe
  const NamedFunc SignalPhoton_hoe = FilterNamedFunc(
      "Photon_hoe",Photon_sig).Name("SignalPhoton_hoe");

  //signal photon mvaID
  const NamedFunc SignalPhoton_mvaID = FilterNamedFunc(
      "Photon_mvaID",Photon_sig).Name("SignalPhoton_mvaID");

  //lead signal photon pt
  const NamedFunc Lead_SignalPhoton_pt = ReduceNamedFunc(
      SignalPhoton_pt, reduce_max).Name("Lead_SignalPhoton_pt");

  //lead signal photon eta
  const NamedFunc Lead_SignalPhoton_eta = MultiReduceNamedFunc(
      {SignalPhoton_pt,SignalPhoton_eta}, reduce_maxfirst)
      .Name("Lead_SignalPhoton_eta");

  //lead signal photon phi
  const NamedFunc Lead_SignalPhoton_phi = MultiReduceNamedFunc(
      {SignalPhoton_pt,SignalPhoton_phi}, reduce_maxfirst)
      .Name("Lead_SignalPhoton_phi");

  //lead signal photon r9
  const NamedFunc Lead_SignalPhoton_r9 = MultiReduceNamedFunc(
      {SignalPhoton_pt,SignalPhoton_r9}, reduce_maxfirst)
      .Name("Lead_SignalPhoton_r9");

  //lead signal photon sieie
  const NamedFunc Lead_SignalPhoton_sieie = MultiReduceNamedFunc(
      {SignalPhoton_pt,SignalPhoton_sieie}, reduce_maxfirst)
      .Name("Lead_SignalPhoton_sieie");

  //lead signal photon pfRelIso03_all
  const NamedFunc Lead_SignalPhoton_pfRelIso03_all = MultiReduceNamedFunc(
      {SignalPhoton_pt,SignalPhoton_pfRelIso03_all}, reduce_maxfirst)
      .Name("Lead_SignalPhoton_pfRelIso03_all");

  //lead signal photon hoe
  const NamedFunc Lead_SignalPhoton_hoe = MultiReduceNamedFunc(
      {SignalPhoton_pt,SignalPhoton_hoe}, reduce_maxfirst)
      .Name("Lead_SignalPhoton_hoe");

  //lead signal photon mvaID
  const NamedFunc Lead_SignalPhoton_mvaID = MultiReduceNamedFunc(
      {SignalPhoton_pt,SignalPhoton_mvaID}, reduce_maxfirst)
      .Name("Lead_SignalPhoton_mvaID");

  //lead signal photon abseta
  const NamedFunc Lead_SignalPhoton_abseta = MapNamedFunc(
      Lead_SignalPhoton_eta,dabs).Name("Lead_SignalPhoton_abseta");

  //this differs from the version in nano_functions 
  const NamedFunc Remove_overlap("Remove_overlap",[SignalPhoton_genPartFlav](const Baby &b) -> NamedFunc::ScalarType{
    //skip data
    if (b.SampleTypeString().Contains("-")) return 1;
    //for MC, just use existence of true photon
    std::vector<double> photon_pflavor = SignalPhoton_genPartFlav.GetVector(b);
    bool is_real_photon = false;
    if (photon_pflavor.size()>0)
      if (photon_pflavor[0] == 1)
        is_real_photon = true;
    if ((b.FirstFileName().find("DYJets") != std::string::npos) && is_real_photon) return 0;
    return 1;
  });

  //returns a weight that will scale events appropriately for 1 fb-1
  const NamedFunc w_lumi("w_lumi",[year](const Baby &b) -> NamedFunc::ScalarType{
    //skip data
    if (b.SampleTypeString().Contains("-")) return 1;
    //for MC, base on sample name
    float xs = 0.0; //pb
    float genEventSumw = 0.0; //from NanoAOD files
    if (b.FirstFileName().find("DYJetsToLL") != std::string::npos) {
      xs = 6077.22;
      if (year == "2016APV")
        genEventSumw = 1.54570798156e+12;
      else 
        genEventSumw = 1.22093462796e+12; //2016
    }
    if (b.FirstFileName().find("ZGToLLG") != std::string::npos) {
      xs = 172.4;
      if (year == "2016APV")
        genEventSumw = 13003821894.7;
      else
        genEventSumw = 13107255164.3; //2016
    }
    return b.genWeight()/genEventSumw*(xs*1000.0);
  });

  //returns photon scale factors (WP80 * CSEV)
  const NamedFunc w_photon("w_photon",[Lead_SignalPhoton_pt,Lead_SignalPhoton_eta,Lead_SignalPhoton_r9,year](const Baby &b) -> NamedFunc::ScalarType{
    float photon_pt = Lead_SignalPhoton_pt.GetScalar(b);
    float photon_eta = Lead_SignalPhoton_eta.GetScalar(b);
    float photon_r9 = Lead_SignalPhoton_r9.GetScalar(b);
    //skip data
    if (b.SampleTypeString().Contains("-")) return 1.0;
    //2016APV
    if (year == "2016APV") {
      float csev_sf = 1.0;
      if (fabs(photon_eta) < 1.5 && photon_r9 < 0.94) csev_sf = 0.9933314919471741;
      if (fabs(photon_eta) < 1.5 && photon_r9 > 0.94) csev_sf = 0.9969298243522644;
      if (fabs(photon_eta) > 1.5 && photon_r9 < 0.94) csev_sf = 0.9939557313919067;
      if (fabs(photon_eta) > 1.5 && photon_r9 > 0.94) csev_sf = 0.9926892518997192;
      if (photon_pt > 20 && photon_pt < 35  && photon_eta > -2.5 && photon_eta < -2.0)     return csev_sf*1.0024752616882324;
      if (photon_pt > 35 && photon_pt < 50  && photon_eta > -2.5 && photon_eta < -2.0)     return csev_sf*1.0034924745559692;
      if (photon_pt > 50 && photon_pt < 80  && photon_eta > -2.5 && photon_eta < -2.0)     return csev_sf*1.0057936906814575;
      if (photon_pt > 80 && photon_pt < 120 && photon_eta > -2.5 && photon_eta < -2.0)     return csev_sf*1.0311778783798218;
      if (photon_pt > 120                   && photon_eta > -2.5 && photon_eta < -2.0)     return csev_sf*0.9780092835426331;
      if (photon_pt > 20 && photon_pt < 35  && photon_eta > -2.0 && photon_eta < -1.566)   return csev_sf*0.9536231756210327;
      if (photon_pt > 35 && photon_pt < 50  && photon_eta > -2.0 && photon_eta < -1.566)   return csev_sf*0.9685863852500916;
      if (photon_pt > 50 && photon_pt < 80  && photon_eta > -2.0 && photon_eta < -1.566)   return csev_sf*0.9739243984222412;
      if (photon_pt > 80 && photon_pt < 120 && photon_eta > -2.0 && photon_eta < -1.566)   return csev_sf*0.9781209826469421;
      if (photon_pt > 120                   && photon_eta > -2.0 && photon_eta < -1.566)   return csev_sf*1.0259740352630615;
      if (photon_pt > 20 && photon_pt < 35  && photon_eta > -1.566 && photon_eta < -1.444) return csev_sf*1.0212417840957642;
      if (photon_pt > 35 && photon_pt < 50  && photon_eta > -1.566 && photon_eta < -1.444) return csev_sf*1.0039734840393066;
      if (photon_pt > 50 && photon_pt < 80  && photon_eta > -1.566 && photon_eta < -1.444) return csev_sf*0.9881734848022461;
      if (photon_pt > 80 && photon_pt < 120 && photon_eta > -1.566 && photon_eta < -1.444) return csev_sf*0.9838709831237793;
      if (photon_pt > 120                   && photon_eta > -1.566 && photon_eta < -1.444) return csev_sf*0.9812080264091492;
      if (photon_pt > 20 && photon_pt < 35  && photon_eta > -1.444 && photon_eta < -0.8)   return csev_sf*0.9790209531784058;
      if (photon_pt > 35 && photon_pt < 50  && photon_eta > -1.444 && photon_eta < -0.8)   return csev_sf*0.9789081811904907;
      if (photon_pt > 50 && photon_pt < 80  && photon_eta > -1.444 && photon_eta < -0.8)   return csev_sf*0.9793438911437988;
      if (photon_pt > 80 && photon_pt < 120 && photon_eta > -1.444 && photon_eta < -0.8)   return csev_sf*1.0096617937088013;
      if (photon_pt > 120                   && photon_eta > -1.444 && photon_eta < -0.8)   return csev_sf*0.9868578314781189;
      if (photon_pt > 20 && photon_pt < 35  && photon_eta > -0.8 && photon_eta < 0.0)      return csev_sf*0.9732262492179871;
      if (photon_pt > 35 && photon_pt < 50  && photon_eta > -0.8 && photon_eta < 0.0)      return csev_sf*0.9724220633506775;
      if (photon_pt > 50 && photon_pt < 80  && photon_eta > -0.8 && photon_eta < 0.0)      return csev_sf*0.9681227803230286;
      if (photon_pt > 80 && photon_pt < 120 && photon_eta > -0.8 && photon_eta < 0.0)      return csev_sf*0.9953051805496216;
      if (photon_pt > 120                   && photon_eta > -0.8 && photon_eta < 0.0)      return csev_sf*0.9964747428894043;
      if (photon_pt > 20 && photon_pt < 35  && photon_eta > 0.0 && photon_eta < 0.8)       return csev_sf*0.9677852392196655;
      if (photon_pt > 35 && photon_pt < 50  && photon_eta > 0.0 && photon_eta < 0.8)       return csev_sf*0.96875;
      if (photon_pt > 50 && photon_pt < 80  && photon_eta > 0.0 && photon_eta < 0.8)       return csev_sf*0.9656398296356201;
      if (photon_pt > 80 && photon_pt < 120 && photon_eta > 0.0 && photon_eta < 0.8)       return csev_sf*0.9976218938827515;
      if (photon_pt > 120                   && photon_eta > 0.0 && photon_eta < 0.8)       return csev_sf*0.9988207817077637;
      if (photon_pt > 20 && photon_pt < 35  && photon_eta > 0.8 && photon_eta < 1.444)     return csev_sf*0.9729729890823364;
      if (photon_pt > 35 && photon_pt < 50  && photon_eta > 0.8 && photon_eta < 1.444)     return csev_sf*0.9761605858802795;
      if (photon_pt > 50 && photon_pt < 80  && photon_eta > 0.8 && photon_eta < 1.444)     return csev_sf*0.9767156839370728;
      if (photon_pt > 80 && photon_pt < 120 && photon_eta > 0.8 && photon_eta < 1.444)     return csev_sf*1.0;
      if (photon_pt > 120                   && photon_eta > 0.8 && photon_eta < 1.444)     return csev_sf*1.0;
      if (photon_pt > 20 && photon_pt < 35  && photon_eta > 1.444 && photon_eta < 1.566)   return csev_sf*0.949999988079071;
      if (photon_pt > 35 && photon_pt < 50  && photon_eta > 1.444 && photon_eta < 1.566)   return csev_sf*0.9959999918937683;
      if (photon_pt > 50 && photon_pt < 80  && photon_eta > 1.444 && photon_eta < 1.566)   return csev_sf*0.9920739531517029;
      if (photon_pt > 80 && photon_pt < 120 && photon_eta > 1.444 && photon_eta < 1.566)   return csev_sf*0.9945651888847351;
      if (photon_pt > 120                   && photon_eta > 1.444 && photon_eta < 1.566)   return csev_sf*1.0809327363967896;
      if (photon_pt > 20 && photon_pt < 35  && photon_eta > 1.566 && photon_eta < 2.0)     return csev_sf*0.9486803412437439;
      if (photon_pt > 35 && photon_pt < 50  && photon_eta > 1.566 && photon_eta < 2.0)     return csev_sf*0.9668874144554138;
      if (photon_pt > 50 && photon_pt < 80  && photon_eta > 1.566 && photon_eta < 2.0)     return csev_sf*0.9722222089767456;
      if (photon_pt > 80 && photon_pt < 120 && photon_eta > 1.566 && photon_eta < 2.0)     return csev_sf*1.0092226266860962;
      if (photon_pt > 120                   && photon_eta > 1.566 && photon_eta < 2.0)     return csev_sf*0.9935064911842346;
      if (photon_pt > 20 && photon_pt < 35  && photon_eta > 2.0 && photon_eta < 2.5)       return csev_sf*1.003731369972229;
      if (photon_pt > 35 && photon_pt < 50  && photon_eta > 2.0 && photon_eta < 2.5)       return csev_sf*1.0023282766342163;
      if (photon_pt > 50 && photon_pt < 80  && photon_eta > 2.0 && photon_eta < 2.5)       return csev_sf*1.0115875005722046;
      if (photon_pt > 80 && photon_pt < 120 && photon_eta > 2.0 && photon_eta < 2.5)       return csev_sf*1.0520833730697632;
      if (photon_pt > 120                   && photon_eta > 2.0 && photon_eta < 2.5)       return csev_sf*1.0467835664749146;
      return 1.0;
    }
    else if (year == "2016") {
      float csev_sf = 1.0;
      if (fabs(photon_eta) < 1.5 && photon_r9 < 0.94) csev_sf = 0.9962645769119263;
      if (fabs(photon_eta) < 1.5 && photon_r9 > 0.94) csev_sf = 0.9962645769119263;
      if (fabs(photon_eta) > 1.5 && photon_r9 < 0.94) csev_sf = 0.9833487868309021;
      if (fabs(photon_eta) > 1.5 && photon_r9 > 0.94) csev_sf = 0.9887195229530334;
      if (photon_pt > 20 && photon_pt < 35  && photon_eta > -2.5 && photon_eta < -2.0)     return csev_sf*0.986369252204895;
      if (photon_pt > 35 && photon_pt < 50  && photon_eta > -2.5 && photon_eta < -2.0)     return csev_sf*0.9906976819038391;
      if (photon_pt > 50 && photon_pt < 80  && photon_eta > -2.5 && photon_eta < -2.0)     return csev_sf*0.9907727837562561;
      if (photon_pt > 80 && photon_pt < 120 && photon_eta > -2.5 && photon_eta < -2.0)     return csev_sf*1.0069364309310913;
      if (photon_pt > 120                   && photon_eta > -2.5 && photon_eta < -2.0)     return csev_sf*1.0261067152023315;
      if (photon_pt > 20 && photon_pt < 35  && photon_eta > -2.0 && photon_eta < -1.566)   return csev_sf*0.9362319111824036;
      if (photon_pt > 35 && photon_pt < 50  && photon_eta > -2.0 && photon_eta < -1.566)   return csev_sf*0.9565789699554443;
      if (photon_pt > 50 && photon_pt < 80  && photon_eta > -2.0 && photon_eta < -1.566)   return csev_sf*0.9634941220283508;
      if (photon_pt > 80 && photon_pt < 120 && photon_eta > -2.0 && photon_eta < -1.566)   return csev_sf*0.9935897588729858;
      if (photon_pt > 120                   && photon_eta > -2.0 && photon_eta < -1.566)   return csev_sf*1.015424132347107;
      if (photon_pt > 20 && photon_pt < 35  && photon_eta > -1.566 && photon_eta < -1.444) return csev_sf*0.9869494438171387;
      if (photon_pt > 35 && photon_pt < 50  && photon_eta > -1.566 && photon_eta < -1.444) return csev_sf*0.9893758296966553;
      if (photon_pt > 50 && photon_pt < 80  && photon_eta > -1.566 && photon_eta < -1.444) return csev_sf*0.9828495979309082;
      if (photon_pt > 80 && photon_pt < 120 && photon_eta > -1.566 && photon_eta < -1.444) return csev_sf*1.0928667783737183;
      if (photon_pt > 120                   && photon_eta > -1.566 && photon_eta < -1.444) return csev_sf*0.9521276354789734;
      if (photon_pt > 20 && photon_pt < 35  && photon_eta > -1.444 && photon_eta < -0.8)   return csev_sf*0.9609484076499939;
      if (photon_pt > 35 && photon_pt < 50  && photon_eta > -1.444 && photon_eta < -0.8)   return csev_sf*0.9678217768669128;
      if (photon_pt > 50 && photon_pt < 80  && photon_eta > -1.444 && photon_eta < -0.8)   return csev_sf*0.9672330021858215;
      if (photon_pt > 80 && photon_pt < 120 && photon_eta > -1.444 && photon_eta < -0.8)   return csev_sf*1.010830283164978;
      if (photon_pt > 120                   && photon_eta > -1.444 && photon_eta < -0.8)   return csev_sf*0.9740259647369385;
      if (photon_pt > 20 && photon_pt < 35  && photon_eta > -0.8 && photon_eta < 0.0)      return csev_sf*0.9588313698768616;
      if (photon_pt > 35 && photon_pt < 50  && photon_eta > -0.8 && photon_eta < 0.0)      return csev_sf*0.9642431735992432;
      if (photon_pt > 50 && photon_pt < 80  && photon_eta > -0.8 && photon_eta < 0.0)      return csev_sf*0.9623529314994812;
      if (photon_pt > 80 && photon_pt < 120 && photon_eta > -0.8 && photon_eta < 0.0)      return csev_sf*0.9917840361595154;
      if (photon_pt > 120                   && photon_eta > -0.8 && photon_eta < 0.0)      return csev_sf*0.9859976768493652;
      if (photon_pt > 20 && photon_pt < 35  && photon_eta > 0.0 && photon_eta < 0.8)       return csev_sf*0.97428959608078;
      if (photon_pt > 35 && photon_pt < 50  && photon_eta > 0.0 && photon_eta < 0.8)       return csev_sf*0.9733333587646484;
      if (photon_pt > 50 && photon_pt < 80  && photon_eta > 0.0 && photon_eta < 0.8)       return csev_sf*0.970202624797821;
      if (photon_pt > 80 && photon_pt < 120 && photon_eta > 0.0 && photon_eta < 0.8)       return csev_sf*1.004756212234497;
      if (photon_pt > 120                   && photon_eta > 0.0 && photon_eta < 0.8)       return csev_sf*0.9904988408088684;
      if (photon_pt > 20 && photon_pt < 35  && photon_eta > 0.8 && photon_eta < 1.444)     return csev_sf*0.9615384340286255;
      if (photon_pt > 35 && photon_pt < 50  && photon_eta > 0.8 && photon_eta < 1.444)     return csev_sf*0.9723618030548096;
      if (photon_pt > 50 && photon_pt < 80  && photon_eta > 0.8 && photon_eta < 1.444)     return csev_sf*0.977886974811554;
      if (photon_pt > 80 && photon_pt < 120 && photon_eta > 0.8 && photon_eta < 1.444)     return csev_sf*0.9939467310905457;
      if (photon_pt > 120                   && photon_eta > 0.8 && photon_eta < 1.444)     return csev_sf*1.0084540843963623;
      if (photon_pt > 20 && photon_pt < 35  && photon_eta > 1.444 && photon_eta < 1.566)   return csev_sf*0.9784052968025208;
      if (photon_pt > 35 && photon_pt < 50  && photon_eta > 1.444 && photon_eta < 1.566)   return csev_sf*0.9906542301177979;
      if (photon_pt > 50 && photon_pt < 80  && photon_eta > 1.444 && photon_eta < 1.566)   return csev_sf*0.9947159886360168;
      if (photon_pt > 80 && photon_pt < 120 && photon_eta > 1.444 && photon_eta < 1.566)   return csev_sf*0.9525745511054993;
      if (photon_pt > 120                   && photon_eta > 1.444 && photon_eta < 1.566)   return csev_sf*0.926630437374115;
      if (photon_pt > 20 && photon_pt < 35  && photon_eta > 1.566 && photon_eta < 2.0)     return csev_sf*0.9335302710533142;
      if (photon_pt > 35 && photon_pt < 50  && photon_eta > 1.566 && photon_eta < 2.0)     return csev_sf*0.9639519453048706;
      if (photon_pt > 50 && photon_pt < 80  && photon_eta > 1.566 && photon_eta < 2.0)     return csev_sf*0.9669749140739441;
      if (photon_pt > 80 && photon_pt < 120 && photon_eta > 1.566 && photon_eta < 2.0)     return csev_sf*0.9934123754501343;
      if (photon_pt > 120                   && photon_eta > 1.566 && photon_eta < 2.0)     return csev_sf*0.9935317039489746;
      if (photon_pt > 20 && photon_pt < 35  && photon_eta > 2.0 && photon_eta < 2.5)       return csev_sf*0.9924905896186829;
      if (photon_pt > 35 && photon_pt < 50  && photon_eta > 2.0 && photon_eta < 2.5)       return csev_sf*0.9964953064918518;
      if (photon_pt > 50 && photon_pt < 80  && photon_eta > 2.0 && photon_eta < 2.5)       return csev_sf*1.0104408264160156;
      if (photon_pt > 80 && photon_pt < 120 && photon_eta > 2.0 && photon_eta < 2.5)       return csev_sf*1.0093023777008057;
      if (photon_pt > 120                   && photon_eta > 2.0 && photon_eta < 2.5)       return csev_sf*1.0228049755096436;
      return 1.0;
    }
    return 1.0;
  });

  //Z (mumuph) candidate mass
  const NamedFunc ZCandidate_mass("ZCandidate_mass",[Photon_sig](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> muon_sig = Muon_sig.GetVector(b);
    std::vector<double> photon_sig = Photon_sig.GetVector(b);
    ROOT::Math::PtEtaPhiMVector mu1, mu2, ph;
    bool first_muon = true;
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      if (muon_sig[imu]) {
        if (first_muon) {
          mu1.SetCoordinates(b.Muon_pt()->at(imu),b.Muon_eta()->at(imu),b.Muon_phi()->at(imu),0.106);
          first_muon = false;
        }
        else {
          mu2.SetCoordinates(b.Muon_pt()->at(imu),b.Muon_eta()->at(imu),b.Muon_phi()->at(imu),0.106);
          break;
        }
      }
    }
    for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
      if (photon_sig[iph]) {
        ph.SetCoordinates(b.Photon_pt()->at(iph),b.Photon_eta()->at(iph),b.Photon_phi()->at(iph),0.0);
        break;
      }
    }
    return (mu1+mu2+ph).M();
  });
  
  //------------------------------------------------------------------------------------
  //                                   plots and tables
  //------------------------------------------------------------------------------------

  PlotMaker pm;
  pm.multithreaded_ = true;
  //pm.max_entries_ = 1000; //for debugging

  ////inclusive
  //pm.Push<Hist1D>(Axis(30,80.0,100.0, ZCandidate_mass, "m_{#mu#mu#gamma} [GeV]", {}), 
  //    Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //    ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25, 
  //    procs, ops).Weight(w_lumi*w_years).Tag("zgapv_"+year);
  //pm.Push<Hist1D>(Axis(30,20.0,80.0, Lead_SignalPhoton_pt, "Photon p_{T} [GeV]", {}), 
  //    Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //    ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25, 
  //    procs, ops).Weight(w_lumi*w_years).Tag("zgapv_"+year);
  //pm.Push<Hist1D>(Axis(30,-2.5,2.5, Lead_SignalPhoton_eta, "Photon #eta", {}), 
  //    Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //    ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25, 
  //    procs, ops).Weight(w_lumi*w_years).Tag("zgapv_"+year);
  //pm.Push<Hist1D>(Axis(30,-3.16,3.16, Lead_SignalPhoton_phi, "Photon #phi", {}), 
  //    Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //    ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25, 
  //    procs, ops).Weight(w_lumi*w_years).Tag("zgapv_"+year);
  //pm.Push<Hist1D>(Axis(30,0.0,1.0, Lead_SignalPhoton_r9, "Photon R_{9}", {}), 
  //    Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //    ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25, 
  //    procs, ops).Weight(w_lumi*w_years).Tag("zgapv_"+year);
  //pm.Push<Hist1D>(Axis(30,0.0,0.07, Lead_SignalPhoton_sieie, "Photon #sigma_{i#eta i#eta}", {}), 
  //    Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //    ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25, 
  //    procs, ops).Weight(w_lumi*w_years).Tag("zgapv_"+year);
  //pm.Push<Hist1D>(Axis(30,0.0,1.0, Lead_SignalPhoton_hoe, "Photon H/E", {}), 
  //    Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //    ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25, 
  //    procs, ops).Weight(w_lumi*w_years).Tag("zgapv_"+year);
  //pm.Push<Hist1D>(Axis(30,0.0,1.0, Lead_SignalPhoton_pfRelIso03_all, "Photon I_{rel 0.3}", {}), 
  //    Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //    ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25, 
  //    procs, ops).Weight(w_lumi*w_years).Tag("zgapv_"+year);
  //pm.Push<Hist1D>(Axis(30,0.0,1.0, Lead_SignalPhoton_mvaID, "Photon Fall17v2 MVA ID score", {}), 
  //    Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //    ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25, 
  //    procs, ops).Weight(w_lumi*w_years).Tag("zgapv_"+year);

  ////pT-binned
    vector<double> pt_bins = {20.0, 35.0, 50.0, 80.0};
    vector<string> pt_bin_names = {"20<photon_pt<35","35<photon_pt<50","50<photon_pt<80"};

  //for (unsigned ipt = 0; ipt < pt_bins.size()-1; ipt++) {
  //  NamedFunc photon_pt_cut = NamedFunc(Lead_SignalPhoton_pt>pt_bins[ipt]&&
  //      Lead_SignalPhoton_pt<pt_bins[ipt+1]).Name(pt_bin_names[ipt]);
  //  pm.Push<Hist1D>(Axis(30,-2.5,2.5, Lead_SignalPhoton_eta, "Photon #eta", {}), 
  //      Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //      ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25&&
  //      photon_pt_cut, procs, ops).Weight(w_lumi*w_years).Tag("zgapv_"+year);
  //  //comment these out if uninteresting
  //  pm.Push<Hist1D>(Axis(30,80.0,100.0, ZCandidate_mass, "m_{#mu#mu#gamma} [GeV]", {}), 
  //      Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //      ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25&&
  //      photon_pt_cut, procs, ops).Weight(w_lumi*w_years).Tag("zgapv_"+year);
  //  pm.Push<Hist1D>(Axis(30,-3.16,3.16, Lead_SignalPhoton_phi, "Photon #phi", {}), 
  //      Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //      ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25&&
  //      photon_pt_cut, procs, ops).Weight(w_lumi*w_years).Tag("zgapv_"+year);
  //  pm.Push<Hist1D>(Axis(30,0.0,1.0, Lead_SignalPhoton_r9, "Photon R_{9}", {}), 
  //      Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //      ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25&&
  //      photon_pt_cut, procs, ops).Weight(w_lumi*w_years).Tag("zgapv_"+year);
  //  pm.Push<Hist1D>(Axis(30,0.0,0.07, Lead_SignalPhoton_sieie, "Photon #sigma_{i#eta i#eta}", {}), 
  //      Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //      ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25&& 
  //      photon_pt_cut, procs, ops).Weight(w_lumi*w_years).Tag("zgapv_"+year);
  //  pm.Push<Hist1D>(Axis(30,0.0,1.0, Lead_SignalPhoton_hoe, "Photon H/E", {}), 
  //      Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //      ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25&& 
  //      photon_pt_cut, procs, ops).Weight(w_lumi*w_years).Tag("zgapv_"+year);
  //  pm.Push<Hist1D>(Axis(30,0.0,1.0, Lead_SignalPhoton_pfRelIso03_all, "Photon I_{rel 0.3}", {}), 
  //      Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //      ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25&& 
  //      photon_pt_cut, procs, ops).Weight(w_lumi*w_years).Tag("zgapv_"+year);
  //  pm.Push<Hist1D>(Axis(30,0.0,1.0, Lead_SignalPhoton_mvaID, "Photon Fall17v2 MVA ID score", {}), 
  //      Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //      ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25&& 
  //      photon_pt_cut, procs, ops).Weight(w_lumi*w_years).Tag("zgapv_"+year);
  //}

  ////eta-binned
    vector<double> abseta_bins = {0.0, 0.8, 1.5, 2.0, 2.5};
    vector<string> abseta_bin_names = {"0<photon_abseta<0.8","0.8<photon_abseta<1.444","1.556<photon_abseta<2.0","2.0<photon_abseta<2.5"};

  //for (unsigned ieta = 0; ieta < abseta_bins.size()-1; ieta++) {
  //  NamedFunc photon_eta_cut = NamedFunc(Lead_SignalPhoton_abseta>abseta_bins[ieta]&&
  //      Lead_SignalPhoton_abseta<abseta_bins[ieta+1]).Name(abseta_bin_names[ieta]);
  //  pm.Push<Hist1D>(Axis(30,20.0,80.0, Lead_SignalPhoton_pt, "Photon p_{T} [GeV]", {}), 
  //      Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //      ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25&&
  //      photon_eta_cut, procs, ops).Weight(w_lumi*w_years).Tag("zgapv_"+year);
  //  pm.Push<Hist1D>(Axis(30,80.0,100.0, ZCandidate_mass, "m_{#mu#mu#gamma} [GeV]", {}), 
  //      Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //      ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25&&
  //      photon_eta_cut, procs, ops).Weight(w_lumi*w_years).Tag("zgapv_"+year);
  //  pm.Push<Hist1D>(Axis(30,-3.16,3.16, Lead_SignalPhoton_phi, "Photon #phi", {}), 
  //      Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //      ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25&&
  //      photon_eta_cut, procs, ops).Weight(w_lumi*w_years).Tag("zgapv_"+year);
  //  pm.Push<Hist1D>(Axis(30,0.0,1.0, Lead_SignalPhoton_r9, "Photon R_{9}", {}), 
  //      Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //      ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25&&
  //      photon_eta_cut, procs, ops).Weight(w_lumi*w_years).Tag("zgapv_"+year);
  //  pm.Push<Hist1D>(Axis(30,0.0,0.07, Lead_SignalPhoton_sieie, "Photon #sigma_{i#eta i#eta}", {}), 
  //      Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //      ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25&& 
  //      photon_eta_cut, procs, ops).Weight(w_lumi*w_years).Tag("zgapv_"+year);
  //  pm.Push<Hist1D>(Axis(30,0.0,1.0, Lead_SignalPhoton_hoe, "Photon H/E", {}), 
  //      Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //      ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25&& 
  //      photon_eta_cut, procs, ops).Weight(w_lumi*w_years).Tag("zgapv_"+year);
  //  pm.Push<Hist1D>(Axis(30,0.0,1.0, Lead_SignalPhoton_pfRelIso03_all, "Photon I_{rel 0.3}", {}), 
  //      Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //      ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25&& 
  //      photon_eta_cut, procs, ops).Weight(w_lumi*w_years).Tag("zgapv_"+year);
  //  pm.Push<Hist1D>(Axis(30,0.0,1.0, Lead_SignalPhoton_mvaID, "Photon Fall17v2 MVA ID score", {}), 
  //      Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //      ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25&& 
  //      photon_eta_cut, procs, ops).Weight(w_lumi*w_years).Tag("zgapv_"+year);
  //}

  //with weights
  //pm.Push<Hist1D>(Axis(30,80.0,100.0, ZCandidate_mass, "m_{#mu#mu#gamma} [GeV]", {}), 
  //    Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //    ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25, 
  //    procs, ops).Weight(w_lumi*w_years*w_photon*"L1PreFiringWeight_Nom").Tag("zgapv_"+year);
  //pm.Push<Hist1D>(Axis(30,20.0,80.0, Lead_SignalPhoton_pt, "Photon p_{T} [GeV]", {}), 
  //    Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //    ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25, 
  //    procs, ops).Weight(w_lumi*w_years*w_photon*"L1PreFiringWeight_Nom").Tag("zgapv_"+year);
  //pm.Push<Hist1D>(Axis(30,-2.5,2.5, Lead_SignalPhoton_eta, "Photon #eta", {}), 
  //    Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //    ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25, 
  //    procs, ops).Weight(w_lumi*w_years*w_photon*"L1PreFiringWeight_Nom").Tag("zgapv_"+year);
  //pm.Push<Hist1D>(Axis(30,-3.16,3.16, Lead_SignalPhoton_phi, "Photon #phi", {}), 
  //    Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //    ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25, 
  //    procs, ops).Weight(w_lumi*w_years*w_photon*"L1PreFiringWeight_Nom").Tag("zgapv_"+year);
  //pm.Push<Hist1D>(Axis(30,0.0,1.0, Lead_SignalPhoton_r9, "Photon R_{9}", {}), 
  //    Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //    ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25, 
  //    procs, ops).Weight(w_lumi*w_years*w_photon*"L1PreFiringWeight_Nom").Tag("zgapv_"+year);
  //pm.Push<Hist1D>(Axis(30,0.0,0.07, Lead_SignalPhoton_sieie, "Photon #sigma_{i#eta i#eta}", {}), 
  //    Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //    ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25, 
  //    procs, ops).Weight(w_lumi*w_years*w_photon*"L1PreFiringWeight_Nom").Tag("zgapv_"+year);
  //pm.Push<Hist1D>(Axis(30,0.0,1.0, Lead_SignalPhoton_hoe, "Photon H/E", {}), 
  //    Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //    ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25, 
  //    procs, ops).Weight(w_lumi*w_years*w_photon*"L1PreFiringWeight_Nom").Tag("zgapv_"+year);
  //pm.Push<Hist1D>(Axis(30,0.0,1.0, Lead_SignalPhoton_pfRelIso03_all, "Photon I_{rel 0.3}", {}), 
  //    Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //    ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25, 
  //    procs, ops).Weight(w_lumi*w_years*w_photon*"L1PreFiringWeight_Nom").Tag("zgapv_"+year);
  //pm.Push<Hist1D>(Axis(30,0.0,1.0, Lead_SignalPhoton_mvaID, "Photon Fall17v2 MVA ID score", {}), 
  //    Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //    ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25, 
  //    procs, ops).Weight(w_lumi*w_years*w_photon*"L1PreFiringWeight_Nom").Tag("zgapv_"+year);
  //for (unsigned ipt = 0; ipt < pt_bins.size()-1; ipt++) {
  //  NamedFunc photon_pt_cut = NamedFunc(Lead_SignalPhoton_pt>pt_bins[ipt]&&
  //      Lead_SignalPhoton_pt<pt_bins[ipt+1]).Name(pt_bin_names[ipt]);
  //  pm.Push<Hist1D>(Axis(30,-2.5,2.5, Lead_SignalPhoton_eta, "Photon #eta", {}), 
  //      Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //      ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25&&
  //      photon_pt_cut, procs, ops).Weight(w_lumi*w_years*w_photon*"L1PreFiringWeight_Nom").Tag("zgapv_"+year);
  //}
  //for (unsigned ieta = 0; ieta < abseta_bins.size()-1; ieta++) {
  //  NamedFunc photon_eta_cut = NamedFunc(Lead_SignalPhoton_abseta>abseta_bins[ieta]&&
  //      Lead_SignalPhoton_abseta<abseta_bins[ieta+1]).Name(abseta_bin_names[ieta]);
  //  pm.Push<Hist1D>(Axis(30,20.0,80.0, Lead_SignalPhoton_pt, "Photon p_{T} [GeV]", {}), 
  //      Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //      ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25&&
  //      photon_eta_cut, procs, ops).Weight(w_lumi*w_years*w_photon*"L1PreFiringWeight_Nom").Tag("zgapv_"+year);
  //  pm.Push<Hist1D>(Axis(30,80.0,100.0, ZCandidate_mass, "m_{#mu#mu#gamma} [GeV]", {}), 
  //      Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //      ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25&&
  //      photon_eta_cut, procs, ops).Weight(w_lumi*w_years*w_photon*"L1PreFiringWeight_Nom").Tag("zgapv_"+year);
  //  pm.Push<Hist1D>(Axis(30,-3.16,3.16, Lead_SignalPhoton_phi, "Photon #phi", {}), 
  //      Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //      ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25&&
  //      photon_eta_cut, procs, ops).Weight(w_lumi*w_years*w_photon*"L1PreFiringWeight_Nom").Tag("zgapv_"+year);
  //  pm.Push<Hist1D>(Axis(30,0.0,1.0, Lead_SignalPhoton_r9, "Photon R_{9}", {}), 
  //      Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //      ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25&&
  //      photon_eta_cut, procs, ops).Weight(w_lumi*w_years*w_photon*"L1PreFiringWeight_Nom").Tag("zgapv_"+year);
  //  pm.Push<Hist1D>(Axis(30,0.0,0.07, Lead_SignalPhoton_sieie, "Photon #sigma_{i#eta i#eta}", {}), 
  //      Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //      ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25&& 
  //      photon_eta_cut, procs, ops).Weight(w_lumi*w_years*w_photon*"L1PreFiringWeight_Nom").Tag("zgapv_"+year);
  //  pm.Push<Hist1D>(Axis(30,0.0,1.0, Lead_SignalPhoton_hoe, "Photon H/E", {}), 
  //      Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //      ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25&& 
  //      photon_eta_cut, procs, ops).Weight(w_lumi*w_years*w_photon*"L1PreFiringWeight_Nom").Tag("zgapv_"+year);
  //  pm.Push<Hist1D>(Axis(30,0.0,1.0, Lead_SignalPhoton_pfRelIso03_all, "Photon I_{rel 0.3}", {}), 
  //      Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //      ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25&& 
  //      photon_eta_cut, procs, ops).Weight(w_lumi*w_years*w_photon*"L1PreFiringWeight_Nom").Tag("zgapv_"+year);
  //  pm.Push<Hist1D>(Axis(30,0.0,1.0, Lead_SignalPhoton_mvaID, "Photon Fall17v2 MVA ID score", {}), 
  //      Remove_overlap&&nSignalMuon==2&&nSignalPhoton==1&&ZCandidate_mass>80.0&&
  //      ZCandidate_mass<100.0&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25&& 
  //      photon_eta_cut, procs, ops).Weight(w_lumi*w_years*w_photon*"L1PreFiringWeight_Nom").Tag("zgapv_"+year);
  //}

  //ecorr plots
  for (unsigned ieta = 0; ieta < abseta_bins.size()-1; ieta++) {
    NamedFunc photon_eta_cut(abseta_bin_names[ieta],[abseta_bins,ieta](const Baby &b) -> NamedFunc::VectorType{
      std::vector<double> Photon_sig_;
      for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
        double sig = 0;
        if (fabs(b.Photon_eta()->at(iph))>abseta_bins[ieta] && fabs(b.Photon_eta()->at(iph))<abseta_bins[ieta+1]) sig = 1;
        Photon_sig_.push_back(sig);
      }
      return Photon_sig_;
    });
    pm.Push<Hist1D>(Axis(24,0.0,1.2, "Photon_eCorr", "Photon p_{T}^{corrected}/p_{T}^{original}", {}), 
        nSignalMuon==2&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25&&photon_eta_cut,
        procs_data, ops_data).Weight("1").Tag("zgapv_"+year);
    pm.Push<Hist1D>(Axis(24,0.0,1.2, "Photon_eCorr", "Photon p_{T}^{corrected}/p_{T}^{original}", {}), 
        nSignalMuon==2&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25&&photon_eta_cut&&"Photon_r9>0.98",
        procs_data, ops_data).Weight("1").Tag("zgapv_"+year);
    pm.Push<Hist1D>(Axis(24,0.0,1.2, "Photon_eCorr", "Photon p_{T}^{corrected}/p_{T}^{original}", {}), 
        nSignalMuon==2&&"HLT_IsoMu24"&&Lead_SignalMuon_pt>25&&photon_eta_cut&&"Photon_r9<0.98",
        procs_data, ops_data).Weight("1").Tag("zgapv_"+year);
  }

  pm.min_print_ = true;
  pm.SetLuminosityTag(lumi_tag).MakePlots(1.0);

  return 0;
}
