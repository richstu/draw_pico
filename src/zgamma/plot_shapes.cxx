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
                       kRed                      ,{sig_path+"*.root"},   trigs);
  auto proc_hzg_vbf = Process::MakeShared<Baby_pico>("HToZ#gamma VBF(x1000)", sig, 
                       kMagenta                   ,{sig_path+"*VBF*.root"},   trigs);
  auto proc_hzg_gg  = Process::MakeShared<Baby_pico>("HToZ#gamma ggF(x1000)", sig, 
                       TColor::GetColor("#ff0000"),{sig_path+"*GluGlu*.root"},   trigs);
  proc_smzg->SetLineWidth (1); proc_ttg->SetLineWidth  (1);
  proc_dy->SetLineWidth   (1); proc_ewkzg->SetLineWidth(1);
  proc_hzg->SetLineWidth(4);
  proc_dy_nos->SetLineWidth(1);
  proc_hzg_gg->SetLineWidth(4);
  proc_hzg_vbf->SetLineWidth(4);
  vector<shared_ptr<Process>> procs = {proc_dy, proc_ttg, proc_smzg, proc_ewkzg, proc_hzg};
  vector<shared_ptr<Process>> theory_procs = {proc_smzg, proc_dy, proc_ewkzg, proc_hzg};
  PlotOpt log_lumi("txt/plot_styles.txt","CMSPaper");
  log_lumi.Title(TitleType::info)
          .YAxis(YAxisType::log)
          .Stack(StackType::shapes)
          .YTitleOffset(1.)
          .AutoYAxis(false)
          .UseCMYK(false)
          .LeftMargin(0.17)
          .LegendColumns(1)
          .FileExtensions({"pdf"});
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  vector<PlotOpt> ops = {lin_lumi};
  NamedFunc llphoton_cuts("llphoton_m[0]+ll_m[0]>185 && llphoton_m[0] > 100 && llphoton_m[0] < 180 && photon_pt[0]/llphoton_m[0] >= 15./110 && photon_drmin[0] > 0.4");
  vector<NamedFunc> sig_sel = {"llphoton_m[0]+ll_m[0]>185",
                               "llphoton_m[0] > 100 && llphoton_m[0] < 180",
                               "photon_pt[0]/llphoton_m[0] >= 15./110",
                               "photon_drmin[0] > 0.4"};
  NamedFunc baseline("nphoton > 0 && ll_m[0] > 50");
  vector<NamedFunc> lep = {"ll_lepid[0] == 11 && el_pt[ll_i1[0]] > 25 && el_pt[ll_i2[0]] > 15",
                           "ll_lepid[0] == 13 && mu_pt[ll_i1[0]] > 20 && mu_pt[ll_i2[0]] > 10"};
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
    pm.Push<Hist1D>(Axis(50,50,150, "ll_m[0]", "m_{ll} [GeV]",{80,100}),
                    selection, procs, ops).Weight(wgt);
    pm.Push<Hist1D>(Axis(60,100,400, "llphoton_m[0]+ll_m[0]", "m_{ll}+m_{ll#gamma} [GeV]",{185}),
                    lep.at(i) && baseline && NamedFuncUtilities::Nminus1(sig_sel,0), procs, ops).Weight(wgt);
    pm.Push<Hist1D>(Axis(44,80,300, "llphoton_m[0]", "m_{ll#gamma} [GeV]",{100,180}),
                    lep.at(i) && baseline && NamedFuncUtilities::Nminus1(sig_sel,1), procs, ops).Weight(wgt);
    pm.Push<Hist1D>(Axis(25,0,1, "photon_pt[0]/llphoton_m[0]", "p_{T,#gamma}/m_{ll#gamma}",{0.136}),
                    lep.at(i) && baseline && NamedFuncUtilities::Nminus1(sig_sel,2), procs, ops).Weight(wgt);
    pm.Push<Hist1D>(Axis(20,0,2, "photon_drmin[0]", "#DeltaR_{min}(l,#gamma)",{0.4}),
                    lep.at(i) && baseline && NamedFuncUtilities::Nminus1(sig_sel,3), procs, ops).Weight(wgt);
  }
  pm.min_print_ = true;
  pm.MakePlots(35.9);
}

