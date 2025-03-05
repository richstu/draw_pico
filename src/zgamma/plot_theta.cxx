#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <bitset>
#include <unistd.h>
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
#include "core/utilities.hpp"
using namespace std;
using namespace PlotOptTypes;


bool checkBit(int i, int n) {
  return((i%static_cast<int>(pow(2,n+1)))/static_cast<int>(pow(2,n)));
}

int main() {
  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");
  Process::Type back =  Process::Type::background;
  Process::Type sig =  Process::Type::signal;
//   Process::Type data =  Process::Type::data;
  string bfolder("/net/cms29/cms29r0/pico/NanoAODv5/zgamma_channelIslandsv3/2016/");
  string mc_path( bfolder+"mc/merged_zgmc_llg/");
  string sig_path(bfolder+"HToZG/merged_zgmc_llg/");
  NamedFunc el_trigs("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL");
  NamedFunc mu_trigs("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ");
  NamedFunc trigs(el_trigs || mu_trigs);
  auto proc_smzg  = Process::MakeShared<Baby_pico>("SM Z#gamma",       back, 
                       TColor::GetColor("#16bac5"),{mc_path+"*ZGToLLG*"}, trigs);
  auto proc_ewkzg = Process::MakeShared<Baby_pico>("EWK Z#gamma",      back, 
                       TColor::GetColor("#39a9df"),{mc_path+"*LLAJJ*"},   trigs);
  auto proc_dy    = Process::MakeShared<Baby_pico>("DY",               back, 
                       TColor::GetColor("#ffb400"),{mc_path+"*DYJets*"},  trigs && "stitch_dy");
  auto proc_ttg   = Process::MakeShared<Baby_pico>("ttbar",            back, 
                       TColor::GetColor("#ED702D"),{mc_path+"*TT_Tune*"}, trigs);
  auto proc_hzg   = Process::MakeShared<Baby_pico>("HToZ#gamma(x1000)", sig, 
                       TColor::GetColor("#ff0000"),{sig_path+"*.root"},   trigs);
  proc_smzg->SetLineWidth (0);
  proc_ttg->SetLineWidth  (0);
  proc_dy->SetLineWidth   (0);
  proc_ewkzg->SetLineWidth(0);
  proc_hzg->SetLineWidth(3);
  vector<shared_ptr<Process>> procs = {proc_dy, proc_ttg, proc_smzg, proc_ewkzg, proc_hzg};
  vector<shared_ptr<Process>> theory_procs = {proc_smzg, proc_hzg};
  NamedFunc costh_my("cos(#theta_{MY})",[](const Baby &b) -> NamedFunc::ScalarType{
    TLorentzVector lminus, lplus;
    if(b.ll_lepid()->at(0) == 11){
      if(b.el_charge()->at(b.ll_i1()->at(0)) < 0) {
        lminus.SetPtEtaPhiM(b.el_pt() ->at(b.ll_i1()->at(0)),
                            b.el_eta()->at(b.ll_i1()->at(0)),
                            b.el_phi()->at(b.ll_i1()->at(0)),
                            0.00511);
        lplus.SetPtEtaPhiM( b.el_pt() ->at(b.ll_i2()->at(0)),
                            b.el_eta()->at(b.ll_i2()->at(0)),
                            b.el_phi()->at(b.ll_i2()->at(0)),
                            0.00511);
      }
      else {
        lplus.SetPtEtaPhiM( b.el_pt() ->at(b.ll_i1()->at(0)),
                            b.el_eta()->at(b.ll_i1()->at(0)),
                            b.el_phi()->at(b.ll_i1()->at(0)),
                            0.00511);
        lminus.SetPtEtaPhiM(b.el_pt() ->at(b.ll_i2()->at(0)),
                            b.el_eta()->at(b.ll_i2()->at(0)),
                            b.el_phi()->at(b.ll_i2()->at(0)),
                            0.00511);
      }
    }
    else if(b.ll_lepid()->at(0) == 13){
      if(b.mu_charge()->at(b.ll_i1()->at(0)) < 0) {
        lminus.SetPtEtaPhiM(b.mu_pt() ->at(b.ll_i1()->at(0)),
                            b.mu_eta()->at(b.ll_i1()->at(0)),
                            b.mu_phi()->at(b.ll_i1()->at(0)),
                            0.105);
        lplus.SetPtEtaPhiM( b.mu_pt() ->at(b.ll_i2()->at(0)),
                            b.mu_eta()->at(b.ll_i2()->at(0)),
                            b.mu_phi()->at(b.ll_i2()->at(0)),
                            0.105);
      }
      else {
        lplus.SetPtEtaPhiM( b.mu_pt() ->at(b.ll_i1()->at(0)),
                            b.mu_eta()->at(b.ll_i1()->at(0)),
                            b.mu_phi()->at(b.ll_i1()->at(0)),
                            0.105);
        lminus.SetPtEtaPhiM(b.mu_pt() ->at(b.ll_i2()->at(0)),
                            b.mu_eta()->at(b.ll_i2()->at(0)),
                            b.mu_phi()->at(b.ll_i2()->at(0)),
                            0.105);
      }
    }
    TLorentzVector photon;
    photon.SetPtEtaPhiM(b.photon_pt()->at(0),b.photon_eta()->at(0),b.photon_phi()->at(0),0);
    TLorentzVector llg = lminus + lplus + photon;
    TLorentzVector ll = lminus + lplus;
    ll.Boost(-1*llg.BoostVector());
    lminus.Boost(-1*llg.BoostVector());
    lminus.Boost(-1*ll.BoostVector());
    return abs(cos(ll.Angle(lminus.Vect())));
  }); 
  NamedFunc costh_pt1("cos(#theta_{pt1})",[](const Baby &b) -> NamedFunc::ScalarType{
    TLorentzVector l1, l2;
    if(b.ll_lepid()->at(0) == 11){
      l1.SetPtEtaPhiM(b.el_pt() ->at(b.ll_i1()->at(0)),
                      b.el_eta()->at(b.ll_i1()->at(0)),
                      b.el_phi()->at(b.ll_i1()->at(0)),
                      0.00511);
      l2.SetPtEtaPhiM(b.el_pt() ->at(b.ll_i2()->at(0)),
                      b.el_eta()->at(b.ll_i2()->at(0)),
                      b.el_phi()->at(b.ll_i2()->at(0)),
                      0.00511);
    }
    else if(b.ll_lepid()->at(0) == 13){
      l1.SetPtEtaPhiM(b.mu_pt() ->at(b.ll_i1()->at(0)),
                      b.mu_eta()->at(b.ll_i1()->at(0)),
                      b.mu_phi()->at(b.ll_i1()->at(0)),
                      0.105);
      l2.SetPtEtaPhiM(b.mu_pt() ->at(b.ll_i2()->at(0)),
                      b.mu_eta()->at(b.ll_i2()->at(0)),
                      b.mu_phi()->at(b.ll_i2()->at(0)),
                      0.105);
    }
    TLorentzVector photon;
    photon.SetPtEtaPhiM(b.photon_pt()->at(0),b.photon_eta()->at(0),b.photon_phi()->at(0),0);
    TLorentzVector llg = l1 + l2 + photon;
    TLorentzVector ll  = l1 + l2;
    ll.Boost(-1*llg.BoostVector());
    l1.Boost(-1*llg.BoostVector());
    l1.Boost(-1*ll.BoostVector());
    return cos(ll.Angle(l1.Vect()));
//     TLorentzVector ll = l1 + l2;
//     photon.Boost(-ll.BoostVector());
//     TVector3 l(l1.Vect()), p(photon.Vect());
//     double costj = -1*l*p/(l.Mag()*p.Mag());
//     return costj;
  }); 
  NamedFunc abs_costh_pt1("abs(cos(#theta_{pt1}))",[](const Baby &b) -> NamedFunc::ScalarType{
    TLorentzVector l1, l2;
    if(b.ll_lepid()->at(0) == 11){
      l1.SetPtEtaPhiM(b.el_pt() ->at(b.ll_i1()->at(0)),
                      b.el_eta()->at(b.ll_i1()->at(0)),
                      b.el_phi()->at(b.ll_i1()->at(0)),
                      0.00511);
      l2.SetPtEtaPhiM(b.el_pt() ->at(b.ll_i2()->at(0)),
                      b.el_eta()->at(b.ll_i2()->at(0)),
                      b.el_phi()->at(b.ll_i2()->at(0)),
                      0.00511);
    }
    else if(b.ll_lepid()->at(0) == 13){
      l1.SetPtEtaPhiM(b.mu_pt() ->at(b.ll_i1()->at(0)),
                      b.mu_eta()->at(b.ll_i1()->at(0)),
                      b.mu_phi()->at(b.ll_i1()->at(0)),
                      0.105);
      l2.SetPtEtaPhiM(b.mu_pt() ->at(b.ll_i2()->at(0)),
                      b.mu_eta()->at(b.ll_i2()->at(0)),
                      b.mu_phi()->at(b.ll_i2()->at(0)),
                      0.105);
    }
    TLorentzVector photon;
    photon.SetPtEtaPhiM(b.photon_pt()->at(0),b.photon_eta()->at(0),b.photon_phi()->at(0),0);
    TLorentzVector llg = l1 + l2 + photon;
    TLorentzVector ll  = l1 + l2;
    ll.Boost(-1*llg.BoostVector());
    l1.Boost(-1*llg.BoostVector());
    l1.Boost(-1*ll.BoostVector());
    return abs(cos(ll.Angle(l1.Vect())));
  }); 
  NamedFunc costh_pt2("cos(#theta_{pt2})",[](const Baby &b) -> NamedFunc::ScalarType{
    TLorentzVector l1, l2;
    if(b.ll_lepid()->at(0) == 11){
      l1.SetPtEtaPhiM(b.el_pt() ->at(b.ll_i2()->at(0)),
                      b.el_eta()->at(b.ll_i2()->at(0)),
                      b.el_phi()->at(b.ll_i2()->at(0)),
                      0.00511);
      l2.SetPtEtaPhiM(b.el_pt() ->at(b.ll_i1()->at(0)),
                      b.el_eta()->at(b.ll_i1()->at(0)),
                      b.el_phi()->at(b.ll_i1()->at(0)),
                      0.00511);
    }
    else if(b.ll_lepid()->at(0) == 13){
      l1.SetPtEtaPhiM(b.mu_pt() ->at(b.ll_i2()->at(0)),
                      b.mu_eta()->at(b.ll_i2()->at(0)),
                      b.mu_phi()->at(b.ll_i2()->at(0)),
                      0.105);
      l2.SetPtEtaPhiM(b.mu_pt() ->at(b.ll_i1()->at(0)),
                      b.mu_eta()->at(b.ll_i1()->at(0)),
                      b.mu_phi()->at(b.ll_i1()->at(0)),
                      0.105);
    }
    TLorentzVector photon;
    photon.SetPtEtaPhiM(b.photon_pt()->at(0),b.photon_eta()->at(0),b.photon_phi()->at(0),0);
    TLorentzVector llg = l1 + l2 + photon;
    TLorentzVector ll  = l1 + l2;
    ll.Boost(-1*llg.BoostVector());
    l1.Boost(-1*llg.BoostVector());
    l1.Boost(-1*ll.BoostVector());
    return cos(ll.Angle(l1.Vect()));
//     TLorentzVector ll = l1 + l2;
//     photon.Boost(-ll.BoostVector());
//     TVector3 l(l1.Vect()), p(photon.Vect());
//     double costj = -1*l*p/(l.Mag()*p.Mag());
//     return costj;
  }); 
  NamedFunc costh_pico("cos(#theta_{pico})",[](const Baby &b) -> NamedFunc::ScalarType{
    return abs(b.llphoton_costhj()->at(0));
  });
  NamedFunc llphoton_cuts("llphoton_m[0]+ll_m[0]>=185 && llphoton_m[0] > 100 && llphoton_m[0] < 180 && photon_pt[0]/llphoton_m[0] >= 15./110 && photon_drmin[0] > 0.4");
  NamedFunc baseline("nphoton > 0" && "ll_m[0] > 50");
  vector<NamedFunc> lep = {"ll_lepid[0] == 11 && el_pt[ll_i1[0]] > 15 && el_pt[ll_i2[0]] > 15",
                           "ll_lepid[0] == 13 && mu_pt[ll_i1[0]] > 15 && mu_pt[ll_i2[0]] > 15"};
  vector<NamedFunc> cut = {"llphoton_m[0] + ll_m[0] >= 185",
                           "llphoton_m[0] > 100 && llphoton_m[0] < 180",
                           "photon_pt[0]/llphoton_m[0] >= 15./110",
                           "photon_drmin[0] > 0.4"};
  PlotOpt log_lumi("txt/plot_styles.txt","CMSPaper");
  log_lumi.Title(TitleType::info)
          .YAxis(YAxisType::log)
          .Stack(StackType::shapes)
          .Overflow(OverflowType::none)
          .YTitleOffset(1.)
          .AutoYAxis(false)
          .UseCMYK(false)
          .LeftMargin(0.17)
          .LegendColumns(1)
          .CanvasWidth(1077)
          .FileExtensions({"pdf"});
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt lin_stack = lin_lumi().Stack(StackType::signal_overlay);
  PlotOpt log_stack = log_lumi().Stack(StackType::signal_overlay);
  PlotOpt tlin_stack = lin_stack().CanvasWidth(700);
  PlotOpt tlin_lumi = lin_lumi().CanvasWidth(700);
  vector<PlotOpt> ops = {lin_stack,log_lumi};
  vector<PlotOpt> tops = {tlin_stack,tlin_lumi};
  NamedFunc wgt("w_lumi",[](const Baby &b) -> NamedFunc::ScalarType{ 
    double weight = b.w_lumi();
    if(b.type() >= 200000 && b.type() <= 205000)
      return 1000*weight;
    else if(b.type() == 6200)
      return weight/1.664;
    return        weight;
  });
  PlotMaker pm;
  for(int i(0); i < 2; i++) {
    NamedFunc selection = lep.at(i) && baseline && llphoton_cuts;
//     pm.Push<Hist1D>(Axis(20,0,1, costh_my,  "|cos(#theta_{MY})|"  ,{}),     
//                     selection, procs, ops).Weight(wgt);
    pm.Push<Hist1D>(Axis(20,0,1, costh_pt1,  "cos(#theta_{highest p_{T}})" ,{}),    
                    selection, procs, ops).Weight(wgt).Tag("high");
//     pm.Push<Hist1D>(Axis(20,0,1, costh_pt1,  "cos(#theta_{highest p_{T}})" ,{}),    
//                     selection, procs, ops).Weight(wgt).Tag("AN-high");
//     pm.Push<Hist1D>(Axis(20,0,1, abs_costh_pt1,  "|cos(#theta_{highest p_{T}})|" ,{}),    
//                     selection, procs, ops).Weight(wgt).Tag("AN-abs-high");
//     pm.Push<Hist1D>(Axis(40,-1,1, costh_pt2,  "cos(#theta_{lowest p_{T}})" ,{}),    
//                     selection, procs, ops).Weight(wgt).Tag("low");
//     pm.Push<Hist1D>(Axis(40,-1,1, "llphoton_costhj[0]","cos(#theta_{pico})",{}), 
//                     selection, procs, ops).Weight(wgt);
//     pm.Push<Hist1D>(Axis(40,-1,1, "llphoton_costhj[0]","cos(#theta_{pico})",{}), 
//                     lep.at(i) && baseline && "photon_drmin[0] > 0.4 && llphoton_m[0] > 100", theory_procs, tops).Weight(wgt).Tag("theory");
  }
  pm.min_print_ = true;
  pm.MakePlots(35.9);
}

