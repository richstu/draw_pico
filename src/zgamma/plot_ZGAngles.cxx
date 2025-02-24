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
#include "TMVA/Reader.h"
#include "TMVA/Configurable.h"
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
#include "zgamma/zg_utilities.hpp"
using namespace std;
using namespace PlotOptTypes;
using namespace ZgUtilities;


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
  auto proc_dy_nos= Process::MakeShared<Baby_pico>("DY",               back, 
                       TColor::GetColor("#ffb400"),{mc_path+"*DYJets*"},  trigs);
  auto proc_ttg   = Process::MakeShared<Baby_pico>("ttbar",            back, 
                       TColor::GetColor("#ED702D"),{mc_path+"*TT_Tune*"}, trigs);
  auto proc_hzg   = Process::MakeShared<Baby_pico>("HToZ#gamma(x1000)", sig,
                       kRed     ,{sig_path+"*.root"},   trigs);
  auto proc_hzg_2d = Process::MakeShared<Baby_pico>("HToZ#gamma(x1000)", back,
                       kRed     ,{sig_path+"*.root"},   trigs);
//   auto proc_hzg_vbf = Process::MakeShared<Baby_pico>("HToZ#gamma VBF(x1000)", sig, 
//                        kMagenta                   ,{sig_path+"*VBF*.root"},   trigs);
  auto proc_hzg_gg  = Process::MakeShared<Baby_pico>("HToZ#gamma ggF(x1000)", sig, 
                       TColor::GetColor("#ff0000"),{sig_path+"*GluGlu*.root"},   trigs);
  proc_smzg->SetLineWidth (1); proc_ttg->SetLineWidth  (1);
  proc_dy->SetLineWidth   (1); proc_ewkzg->SetLineWidth(1);
  proc_hzg->SetLineWidth(3);
  proc_dy_nos->SetLineWidth(1);
  proc_hzg_gg->SetLineWidth(3);
  vector<shared_ptr<Process>> procs = {proc_dy, proc_smzg, proc_hzg};
  vector<shared_ptr<Process>> theory_procs = {proc_dy, proc_smzg, proc_hzg_gg};
  vector<shared_ptr<Process>> sig_procs = {proc_hzg_2d};
//   vector<shared_ptr<Process>> theory_procs = {proc_hzg_vbf};
  NamedFunc llphoton_pz("llphoton_pz",[](const Baby &b) -> NamedFunc::ScalarType{
    return abs(AssignH(b).Pz());
  });
  NamedFunc llphotonjj_pz("llphotonjj_pz",[](const Baby &b) -> NamedFunc::ScalarType{
    TLorentzVector llphoton, dijet, llgjj;
    llphoton = AssignH(b);
    dijet.SetPtEtaPhiM(b.dijet_pt() ,
                       b.dijet_eta(),
                       b.dijet_phi(),
                       b.dijet_m()  );
    llgjj = llphoton + dijet;
    return abs(llgjj.Pz());
  });
  NamedFunc llphotonjj_pt("llphotonjj_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    TLorentzVector llphoton, dijet, llgjj;
    llphoton = AssignH(b);
    dijet.SetPtEtaPhiM(b.dijet_pt() ,
                       b.dijet_eta(),
                       b.dijet_phi(),
                       b.dijet_m()  );
    llgjj = llphoton + dijet;
    return abs(llgjj.Pt());
  });
  NamedFunc phi("phi",[](const Baby &b) -> NamedFunc::ScalarType{
    return Getphi(b);
  });
  NamedFunc phi_gen("phi_gen",[](const Baby &b) -> NamedFunc::ScalarType{
    return Getphi(b,true);
  });
  NamedFunc costheta("costheta",[](const Baby &b) -> NamedFunc::ScalarType{
    return -1*cos_theta(b);
  });
  NamedFunc costheta_gen("costheta_gen",[](const Baby &b) -> NamedFunc::ScalarType{
    return -1*cos_theta(b,true);
  });
  NamedFunc cosT("coscapTheta",[](const Baby &b) -> NamedFunc::ScalarType{
    return cos_Theta(b);
  });
  NamedFunc cosT_gen("coscapTheta_gen",[](const Baby &b) -> NamedFunc::ScalarType{
    return cos_Theta(b,true);
  });
  NamedFunc pTt("pTt",[](const Baby &b) -> NamedFunc::ScalarType{
    TVector3 g = AssignGamma(b).Vect();
    TVector3 h = AssignH(b).Vect();
    TVector3 z = AssignZ(b).Vect();
    g.SetZ(0); h.SetZ(0); z.SetZ(0);
//     if(z.Angle(g) < 3.14159/2)
//       return h.Mag();
    return h.Cross((z-g).Unit()).Mag();
  });
  NamedFunc pTt2("pTt2",[](const Baby &b) -> NamedFunc::ScalarType{
    TVector3 g = AssignGamma(b).Vect();
    TVector3 h = AssignH(b).Vect();
    TVector3 z = AssignZ(b).Vect();
    g.SetZ(0); h.SetZ(0); z.SetZ(0);
    return h.Cross(z-g).Mag()/h.Mag();
  });
  PlotOpt log_lumi("txt/plot_styles.txt","CMSPaper");
  log_lumi.Title(TitleType::info)
          .YAxis(YAxisType::log)
          .Stack(StackType::shapes)
//           .Overflow(OverflowType::none)
          .YTitleOffset(1.)
          .AutoYAxis(false)
          .UseCMYK(false)
          .LeftMargin(0.17)
          .LegendColumns(1)
//           .CanvasWidth(1077)
          .FileExtensions({"pdf"});
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt lin_stack = lin_lumi().Stack(StackType::signal_overlay);
  PlotOpt log_stack = log_lumi().Stack(StackType::signal_overlay);
  vector<PlotOpt> ops = {lin_stack, lin_lumi};
  PlotOpt lin_lumi_2d("txt/plot_styles.txt", "Scatter");
  vector<PlotOpt> ops2d = {lin_lumi_2d().Title(TitleType::info).YAxis(YAxisType::linear).Overflow(OverflowType::overflow)};
  vector<NamedFunc> loose_lep = {"ll_lepid[0] == 11 && el_pt[ll_i1[0]] > 15 && el_pt[ll_i2[0]] > 15",
                                 "ll_lepid[0] == 13 && mu_pt[ll_i1[0]] > 15 && mu_pt[ll_i2[0]] > 15"};
  NamedFunc llphoton_cuts("llphoton_m[0]+ll_m[0]>185 && llphoton_m[0] > 100 && llphoton_m[0] < 180 && photon_pt[0]/llphoton_m[0] >= 15./110 && photon_drmin[0] > 0.4");
  NamedFunc baseline("nphoton > 0 && ll_m[0] > 50");
  vector<NamedFunc> lep = {"ll_lepid[0] == 11 && el_pt[ll_i1[0]] > 25 && el_pt[ll_i2[0]] > 15",
                           "ll_lepid[0] == 13 && mu_pt[ll_i1[0]] > 20 && mu_pt[ll_i2[0]] > 10"};
  NamedFunc wgt("w_lumi",[](const Baby &b) -> NamedFunc::ScalarType{ 
    double weight = b.w_lumi();
    if(b.type() >= 200000 && b.type() <= 205000)
      return 1000*weight;
//     else if(b.type() == 6200)
//       return weight/1.664;
    return        weight;
  });
  PlotMaker pm;
  string tag;
  for(int i(0); i < 2; i++) {
    for(int j(0); j < 2; j++) {
      NamedFunc cut = loose_lep.at(i) && baseline;
      if(j == 1) cut = lep.at(i) && baseline && llphoton_cuts;
      tag = "Angles";
      if(i == 0) tag += "_el";
      else       tag += "_mu";
      if(j == 0) tag += "_loose";
      else       tag += "_baseline";
      pm.Push<Hist1D>(Axis(20,-1.0,1.0, cosT,        "cos #Theta"      ,{}), cut, theory_procs, ops).Weight(wgt).Tag(tag);
      pm.Push<Hist1D>(Axis(20,-1.0,1.0, costheta,    "cos #theta"      ,{}), cut, theory_procs, ops).Weight(wgt).Tag(tag);
      pm.Push<Hist1D>(Axis(32,-3.2,3.2, phi,         "#phi"            ,{}), cut, theory_procs, ops).Weight(wgt).Tag(tag);
      tag += "_gen";
      pm.Push<Hist1D>(Axis(20,-1.0,1.0, cosT_gen,    "cos #Theta_{gen}",{}), cut, theory_procs, ops).Weight(wgt).Tag(tag);
      pm.Push<Hist1D>(Axis(20,-1.0,1.0, costheta_gen,"cos #theta_{gen}",{}), cut, theory_procs, ops).Weight(wgt).Tag(tag);
      pm.Push<Hist1D>(Axis(32,-3.2,3.2, phi_gen     ,"#phi_{gen}"      ,{}), cut, theory_procs, ops).Weight(wgt).Tag(tag);
      pm.Push<Hist2D>(Axis(20,-1.0,1.0, costheta_gen,"cos #theta_{gen}",{}),
                      Axis(32,-3.2,3.2, phi_gen     ,"#phi_{gen}"      ,{}), cut, sig_procs   , ops2d).Weight(wgt).Tag(tag);
      pm.Push<Hist2D>(Axis(20,-1.0,1.0, costheta    ,"cos #theta"      ,{}),
                      Axis(32,-3.2,3.2, phi         ,"#phi"            ,{}), cut, sig_procs   , ops2d).Weight(wgt).Tag(tag);
//       pm.Push<Hist1D>(Axis(50,0,100, pTt,       "p_{Tt} [GeV]",{}), cut, procs, ops).Weight(wgt);
//       pm.Push<Hist1D>(Axis(50,0,100, pTt2,     "p_{Tt}2 [GeV]",{}), cut, procs, ops).Weight(wgt);
    }
  }
  pm.min_print_ = true;
  pm.MakePlots(35.9);
}

