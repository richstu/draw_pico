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
#include "core/named_func_utilities.hpp"
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
  auto proc_ttg   = Process::MakeShared<Baby_pico>("ttbar",            back, 
                       TColor::GetColor("#ED702D"),{mc_path+"*TT_Tune*"}, trigs);
  auto proc_hzg   = Process::MakeShared<Baby_pico>("HToZ#gamma(x100)", sig, 
                       TColor::GetColor("#ff0000"),{sig_path+"*.root"},   trigs);
  proc_smzg->SetLineWidth (1);
  proc_ttg->SetLineWidth  (1);
  proc_dy->SetLineWidth   (1);
  proc_ewkzg->SetLineWidth(1);
  proc_hzg->SetLineWidth(3);
  NamedFunc ctheta("costheta",   [](const Baby &b) -> NamedFunc::ScalarType{ return cos_theta(b); });
  NamedFunc cTheta("coscapTheta",[](const Baby &b) -> NamedFunc::ScalarType{ return cos_Theta(b); });
  NamedFunc phi(   "phi",        [](const Baby &b) -> NamedFunc::ScalarType{ return Getphi(b); });
  NamedFunc hpt_hm("h_ptdh_m",   [](const Baby &b) -> NamedFunc::ScalarType{ return b.llphoton_pt()->at(0)/b.llphoton_m()->at(0); });
  NamedFunc gpt_hm("g_ptdh_m",   [](const Baby &b) -> NamedFunc::ScalarType{ return b.photon_pt()->at(0)/b.llphoton_m()->at(0); });
  NamedFunc mll_mllg("mll_mllg", [](const Baby &b) -> NamedFunc::ScalarType{ return b.ll_m()->at(0)+b.llphoton_m()->at(0); });
  NamedFunc p_pte( "photon_pte", [](const Baby &b) -> NamedFunc::ScalarType{ return b.photon_pterr()->at(0)/b.photon_pt()->at(0); });
  NamedFunc p_eb(  "photon_EB",  [](const Baby &b) -> NamedFunc::ScalarType{ return abs(b.photon_eta()->at(0)) < 1.4442; });
  NamedFunc p_ee(  "photon_EE",  [](const Baby &b) -> NamedFunc::ScalarType{ return abs(b.photon_eta()->at(0)) > 1.566; });
  NamedFunc photon_drmax("photon_drmax",[](const Baby &b) -> NamedFunc::ScalarType{ 
    TVector3 photon = AssignGamma(b).Vect();
    TVector3 l1     = AssignL1(b).Vect();
    TVector3 l2     = AssignL2(b).Vect();
    return max(photon.DeltaR(l1),photon.DeltaR(l2));
  });
  NamedFunc sysbal("sysbal",[](const Baby &b) -> NamedFunc::ScalarType{ 
    TVector3 z = AssignZ(b).Vect();
    TVector3 g = AssignGamma(b).Vect();
    TVector3 j1, j2;
    int i1(-1), i2(-1);
    for(size_t ij(0); ij < b.jet_pt()->size(); ij++) 
      if(b.jet_isgood()->at(ij) && i1 == -1) i1 = ij;
      else if(b.jet_isgood()->at(ij) && i2 == -1) i2 = ij;
    j1.SetPtEtaPhi(b.jet_pt()->at(i1),b.jet_eta()->at(i1),b.jet_phi()->at(i1));
    j2.SetPtEtaPhi(b.jet_pt()->at(i2),b.jet_eta()->at(i2),b.jet_phi()->at(i2));
//     z.SetZ(0); g.SetZ(0); j1.SetZ(0); j2.SetZ(0);
    return (z+g+j1+j2).Pt()/(z.Pt()+g.Pt()+j1.Pt()+j2.Pt());
  });
  NamedFunc jjllg_dphi("jjllg_dphi",[](const Baby &b) -> NamedFunc::ScalarType{ 
    TVector3 h = AssignH(b).Vect();
    TVector3 jj;
    TVector3 j1, j2;
    int i1(-1), i2(-1);
    for(size_t ij(0); ij < b.jet_pt()->size(); ij++) 
      if(b.jet_isgood()->at(ij) && i1 == -1) i1 = ij;
      else if(b.jet_isgood()->at(ij) && i2 == -1) i2 = ij;
    j1.SetPtEtaPhi(b.jet_pt()->at(i1),b.jet_eta()->at(i1),b.jet_phi()->at(i1));
    j2.SetPtEtaPhi(b.jet_pt()->at(i2),b.jet_eta()->at(i2),b.jet_phi()->at(i2));
    jj = j1 + j2;
    return h.DeltaPhi(jj);
  });
  NamedFunc photon_jdr("photon_jdr",[](const Baby &b) -> NamedFunc::ScalarType{ 
    TVector3 g = AssignGamma(b).Vect();
    TVector3 j1, j2;
    int i1(-1), i2(-1);
    for(size_t ij(0); ij < b.jet_pt()->size(); ij++) 
      if(b.jet_isgood()->at(ij) && i1 == -1) i1 = ij;
      else if(b.jet_isgood()->at(ij) && i2 == -1) i2 = ij;
    j1.SetPtEtaPhi(b.jet_pt()->at(i1),b.jet_eta()->at(i1),b.jet_phi()->at(i1));
    j2.SetPtEtaPhi(b.jet_pt()->at(i2),b.jet_eta()->at(i2),b.jet_phi()->at(i2));
    return min(g.DeltaR(j1),g.DeltaR(j2));
  });
  NamedFunc photon_zep("photon_zep",[](const Baby &b) -> NamedFunc::ScalarType{ 
    int i1(-1), i2(-1);
    for(size_t ij(0); ij < b.jet_pt()->size(); ij++) 
      if(b.jet_isgood()->at(ij) && i1 == -1) i1 = ij;
      else if(b.jet_isgood()->at(ij) && i2 == -1) i2 = ij;
    return abs(b.photon_eta()->at(0) - (b.jet_eta()->at(i1)+b.jet_eta()->at(i2))/2);
  });
  NamedFunc pTt2("pTt2",[](const Baby &b) -> NamedFunc::ScalarType{
    TVector3 g = AssignGamma(b).Vect();
    TVector3 h = AssignH(b).Vect();
    TVector3 z = AssignZ(b).Vect();
    g.SetZ(0); h.SetZ(0); z.SetZ(0);
    return h.Cross(z-g).Mag()/h.Mag();
  });
  NamedFunc jet1_pt("jet1_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    int i1(-1);
    double pt(0);
    for(size_t ij(0); ij < b.jet_pt()->size(); ij++) 
      if(b.jet_isgood()->at(ij) && i1 == -1) {
        i1 = ij;
        pt = b.jet_pt()->at(ij);
      }
    return pt;
  });
  NamedFunc jet2_pt("jet2_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    int i1(-1), i2(-1);
    double pt(0);
    for(size_t ij(0); ij < b.jet_pt()->size(); ij++) 
      if(b.jet_isgood()->at(ij) && i1 == -1)  i1 = ij;
      else if(b.jet_isgood()->at(ij) && i2 == -1) {
        i2 = ij;
        pt = b.jet_pt()->at(ij);
      }
    return pt;
  });
  vector<shared_ptr<Process>> procs = {proc_dy, proc_ttg, proc_smzg,proc_hzg};
  NamedFunc llphoton_cuts("llphoton_m[0]+ll_m[0]>=185 && llphoton_m[0] > 100 && llphoton_m[0] < 180 && photon_pt[0]/llphoton_m[0] >= 15./110 && photon_drmin[0] > 0.4");
//   NamedFunc photon("photon_selection",[](const Baby &b) -> NamedFunc::ScalarType{ 
//    return 1;
//   });

  NamedFunc baseline("nphoton > 0" && "ll_m[0] > 50");
  vector<NamedFunc> lep = {"ll_lepid[0] == 11 && el_pt[ll_i1[0]] > 25 && el_pt[ll_i2[0]] > 15",
                           "ll_lepid[0] == 13 && mu_pt[ll_i1[0]] > 20 && mu_pt[ll_i2[0]] > 10"};
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
          .CanvasWidth(900)
          .FileExtensions({"pdf"});
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt lin_stack = lin_lumi().Stack(StackType::signal_overlay);
  PlotOpt log_stack = log_lumi().Stack(StackType::signal_overlay);
  vector<PlotOpt> ops = {lin_stack};
  NamedFunc wgt("w_lumi",[](const Baby &b) -> NamedFunc::ScalarType{ 
    double weight = b.w_lumi();
    if(b.type() >= 200000 && b.type() <= 205000)
      return 100*weight;
    return       weight;
  });
  PlotMaker pm;
  for(int i(0); i < 2; i++) {
    NamedFunc selection = lep.at(i) && baseline && llphoton_cuts;
//     if(i == 0) {
//       pm.Push<Hist1D>(Axis(15,0,150,"el_pt[ll_i1[0]]",           "Leading e p_{T}",            {}), selection, procs, ops).Weight(wgt);
//       pm.Push<Hist1D>(Axis(15,0,150,"el_pt[ll_i2[0]]",           "Trailng e p_{T}",            {}), selection, procs, ops).Weight(wgt);
//       pm.Push<Hist1D>(Axis(12,-3,3, "el_eta[ll_i1[0]]",          "Leading e #eta",             {}), selection, procs, ops).Weight(wgt);
//     }
//     if(i == 1) {
//       pm.Push<Hist1D>(Axis(15,0,150,"mu_pt[ll_i1[0]]",           "Leading mu p_{T}",           {}), selection, procs, ops).Weight(wgt);
//       pm.Push<Hist1D>(Axis(15,0,150,"mu_pt[ll_i2[0]]",           "Trailng mu p_{T}",           {}), selection, procs, ops).Weight(wgt);
//     }
//     pm.Push<Hist1D>(Axis(20,0,100, "photon_pt[0]", "Leading #gamma p_{T}",       {}), selection, procs, ops).Weight(wgt);
//     pm.Push<Hist1D>(Axis(20,0,1, abs(costh),                          "|cos(#theta_{MY})|",              {}), selection, procs, ops).Weight(wgt);
    pm.Push<Hist1D>(Axis(20,100,180, "llphoton_m[0]", "m_{ll#gamma}",{}), selection, procs, ops).Weight(wgt);
    // Kinematic MVA training variables
    pm.Push<Hist1D>(Axis(20,0,1,      cTheta,                         "cos(#Theta)",                {}), selection         , procs, ops).Weight(wgt).Tag("kin-1");
    pm.Push<Hist1D>(Axis(20,-1,1,     ctheta,                         "cos(#theta)",                {}), selection         , procs, ops).Weight(wgt).Tag("kin-2");
    pm.Push<Hist1D>(Axis(30,-3,3,     phi,                            "#phi",                       {}), selection         , procs, ops).Weight(wgt).Tag("kin-3");
    pm.Push<Hist1D>(Axis(40,0,1,      "llphoton_pt[0]/llphoton_m[0]", "p_{T,ll#gamma}/m_{ll#gamma}",{}), selection         , procs, ops).Weight(wgt).Tag("kin-4");
    if(i == 0){
    pm.Push<Hist1D>(Axis(20,-2.5,2.5, "el_eta[ll_i1[0]]",             "Leading e #eta",             {}), selection         , procs, ops).Weight(wgt).Tag("kin-5");
    pm.Push<Hist1D>(Axis(20,-2.5,2.5, "el_eta[ll_i2[0]]",             "Trailing e #eta",            {}), selection         , procs, ops).Weight(wgt).Tag("kin-6");
    }
    else if(i == 1){
    pm.Push<Hist1D>(Axis(20,-2.5,2.5, "mu_eta[ll_i1[0]]",             "Leading mu #eta",            {}), selection         , procs, ops).Weight(wgt).Tag("kin-5");
    pm.Push<Hist1D>(Axis(20,-2.5,2.5, "mu_eta[ll_i2[0]]",             "Trailing mu #eta",           {}), selection         , procs, ops).Weight(wgt).Tag("kin-6");
    }
    pm.Push<Hist1D>(Axis(20,-2.5,2.5, "photon_eta[0]",                "Photon #eta",                {}), selection         , procs, ops).Weight(wgt).Tag("kin-7");
    pm.Push<Hist1D>(Axis(35,-0.4,1,   "photon_idmva[0]",              "Photon IDMVA (EB)",          {}), selection && p_eb , procs, ops).Weight(wgt).Tag("kin-8");
    pm.Push<Hist1D>(Axis(20,-0.6,1,   "photon_idmva[0]",              "Photon IDMVA (EE)",          {}), selection && p_ee , procs, ops).Weight(wgt).Tag("kin-9");
    pm.Push<Hist1D>(Axis(20,0.4,3.4,  "photon_drmin[0]",              "Min #DeltaR(#gamma,l)",      {}), selection         , procs, ops).Weight(wgt).Tag("kin-10");
    pm.Push<Hist1D>(Axis(20,1.4,4.4,   photon_drmax,                  "Max #DeltaR(#gamma,l)",      {}), selection         , procs, ops).Weight(wgt).Tag("kin-11");
    pm.Push<Hist1D>(Axis(20,0.01,0.11,"photon_pterr[0]/photon_pt[0]", "#sigma_{#gamma}",            {}), selection         , procs, ops).Weight(wgt).Tag("kin-12");
    // Dijet MVA training variables
    selection = "njet >= 2" && lep.at(i) && baseline && llphoton_cuts;
    pm.Push<Hist1D>(Axis(12,30,210 ,jet1_pt     ,"Leading jet p_{T} [GeV]"   ,{}), selection, procs, ops).Weight(wgt).Tag("dijet-1");
    pm.Push<Hist1D>(Axis(12,30,210 ,jet2_pt     ,"Subleading jet p_{T} [GeV]",{}), selection, procs, ops).Weight(wgt).Tag("dijet-2");
    pm.Push<Hist1D>(Axis(10, 0,  1 , sysbal     ,"System balance"            ,{}), selection, procs, ops).Weight(wgt).Tag("dijet-3");
    pm.Push<Hist1D>(Axis(12,0,9    ,"dijet_deta","#Delta#eta(j_{1},j_{2})"   ,{}), selection, procs, ops).Weight(wgt).Tag("dijet-4");
    pm.Push<Hist1D>(Axis(10,0,3    ,"dijet_dphi","#Delta#phi(j_{1},j_{2})"   ,{}), selection, procs, ops).Weight(wgt).Tag("dijet-5");
    pm.Push<Hist1D>(Axis(10,0,3    ,jjllg_dphi  ,"#Delta#phi(Z#gamma,jj)"    ,{}), selection, procs, ops).Weight(wgt).Tag("dijet-6");
    pm.Push<Hist1D>(Axis(15,0.4,4.4,photon_jdr  ,"#DeltaR(#gamma,j)"         ,{}), selection, procs, ops).Weight(wgt).Tag("dijet-7");
    pm.Push<Hist1D>(Axis(20,0,140  ,pTt2        ,"p_{Tt} [GeV]"              ,{}), selection, procs, ops).Weight(wgt).Tag("dijet-8");
    pm.Push<Hist1D>(Axis( 8,0,6    ,photon_zep  ,"#gamma zeppenfeld"         ,{}), selection, procs, ops).Weight(wgt).Tag("dijet-9");
//     pm.Push<Hist1D>(Axis(20,0,1, "photon_r9[0]",   "Photon R9",                  {}), selection, procs, ops).Weight(wgt);
//     pm.Push<Hist1D>(Axis(4,-0.5,3.5,"nllphoton",                 "N_{ll#gamma}",               {}), selection, procs, ops).Weight(wgt);
//     pm.Push<Hist1D>(Axis(4,-0.5,3.5,"nphoton",                   "N_{#gamma}",                 {}), selection, procs, ops).Weight(wgt);
//     pm.Push<Hist1D>(Axis(4,-0.5,3.5,"nll",                       "N_{ll}",                     {}), selection, procs, ops).Weight(wgt);
//     pm.Push<Hist1D>(Axis(4,-0.5,3.5,"nel",                       "N_{e}",                      {}), selection, procs, ops).Weight(wgt);
//     pm.Push<Hist1D>(Axis(4,-0.5,3.5,"nmu",                       "N_{#mu}",                    {}), selection, procs, ops).Weight(wgt);
//     selection = lep.at(i) && baseline;
// //     pm.Push<Hist1D>(Axis(40,0,1,   "photon_idmva[0]", "Leading #gamma p_{T}",       {}), selection, procs, ops).Weight(wgt);
//     pm.Push<Hist1D>(Axis(25,50,300, "llphoton_m[0]+ll_m[0]",     "m_{ll#gamma}+m_{ll} [GeV]",{185}),     selection && NamedFuncUtilities::Nminus1(cut,0), procs, ops).Weight(wgt).Tag("NM1");
//     pm.Push<Hist1D>(Axis(20,50,250, "llphoton_m[0]",             "m_{ll#gamma} [GeV]",       {100,180}), selection && NamedFuncUtilities::Nminus1(cut,1), procs, ops).Weight(wgt).Tag("NM1");
//     pm.Push<Hist1D>(Axis(40,0,2,    "photon_pt[0]/llphoton_m[0]","p_{T,#gamma}/m_{ll#gamma}",{0.136}),   selection && NamedFuncUtilities::Nminus1(cut,2), procs, ops).Weight(wgt).Tag("NM1");
//     pm.Push<Hist1D>(Axis(20,0.05,2.05,"photon_drmin[0]",           "Min #DeltaR(#gamma,l)[0]",    {0.4}),  selection && NamedFuncUtilities::Nminus1(cut,3), procs, ops).Weight(wgt).Tag("NM1");
//     pm.Push<Hist1D>(Axis(40,0.05,4.05,"photon_drmin[0]",           "Min #DeltaR(#gamma,l)[0]",    {0.4}),  selection && NamedFuncUtilities::Nminus1(cut,3), procs, ops).Weight(wgt).Tag("NM1-bin");
  }
  pm.min_print_ = true;
  pm.MakePlots(35.9);
}

