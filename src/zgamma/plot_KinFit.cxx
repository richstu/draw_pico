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
#include "zgamma/zg_utilities.hpp"
using namespace std;
using namespace PlotOptTypes;
using namespace ZgUtilities;

int main() {
  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");
  Process::Type sig  = Process::Type::signal;
  string bfolder("/net/cms17/cms17r0/pico/");
  string ul_path(bfolder+"NanoAODv9/zgamma_signal/2016/signal/merged_zgmc_llg/");

  NamedFunc kinFit("kinFit", [](const Baby &b) -> NamedFunc::ScalarType{
      double mZ_fit = KinRefit(b);
      return mZ_fit;
    });
  NamedFunc kinFitH("kinFitH", [](const Baby &b) -> NamedFunc::ScalarType{
      std::vector<TLorentzVector> reFit = RefitP4(b);
      TLorentzVector H = reFit[0] + reFit[1] + AssignGamma(b);
      return H.M();
    });
  NamedFunc wgt("weight",[](const Baby &b) -> NamedFunc::ScalarType{
      double weight = b.w_lumi();
      return weight*100;
    });

  NamedFunc el_trigs("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
  NamedFunc mu_trigs("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ");
  NamedFunc trigs(el_trigs || mu_trigs);
  NamedFunc mass_cuts("ll_m[llphoton_ill[0]] > 50 && "
  		      "llphoton_m[0]+ll_m[llphoton_ill[0]]>=185 && "
  		      "llphoton_m[0] > 100 && llphoton_m[0] < 180 && "
  		      "photon_pt[llphoton_iph[0]]/llphoton_m[0] >= 15./110");
  NamedFunc baseline("nllphoton > 0");
  vector<NamedFunc> lep = {"ll_lepid[llphoton_ill[0]] == 11 && "
  			   "el_pt[ll_i1[llphoton_ill[0]]] > 25 && "
  			   "el_pt[ll_i2[llphoton_ill[0]]] > 15 && "
  			   "el_sig[ll_i1[llphoton_ill[0]]] && "
  			   "el_sig[ll_i2[llphoton_ill[0]]]",
                           "ll_lepid[llphoton_ill[0]] == 13 && "
  			   "mu_pt[ll_i1[llphoton_ill[0]]] > 20 && "
  			   "mu_pt[ll_i2[llphoton_ill[0]]] > 10 && "
  			   "mu_sig[ll_i1[llphoton_ill[0]]] && "
  			   "mu_sig[ll_i2[llphoton_ill[0]]]"};
  NamedFunc pho("photon_pt[llphoton_iph[0]] > 15 && "
  		"photon_drmin[llphoton_iph[0]] > 0.4 && "
  		"photon_sig[llphoton_iph[0]]");

  PlotOpt log_lumi("txt/plot_styles.txt","CMSPaper");
  log_lumi.Title(TitleType::info)
          .YAxis(YAxisType::log)
          .Stack(StackType::shapes)
          .Overflow(OverflowType::none)
          .YTitleOffset(1.75)
          .AutoYAxis(false)
          .UseCMYK(false)
          .LeftMargin(0.17)
          .LegendColumns(1)
          .CanvasWidth(800)
          .CanvasHeight(800)
          .FileExtensions({"pdf"});

  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt lin_stack = lin_lumi().Stack(StackType::signal_overlay);
  PlotOpt log_stack = log_lumi().Stack(StackType::signal_overlay);
  vector<PlotOpt> ops = {lin_stack};

  auto proc_Red     = Process::MakeShared<Baby_pico>("HToZ#gamma (UL)", sig, 
						     TColor::GetColor("#ff0000"), {ul_path+"*.root"}, trigs);

  auto proc_Blue    = Process::MakeShared<Baby_pico>("HToZ#gamma (UL)", sig,
						     TColor::GetColor("#0000ff"), {ul_path+"*.root"}, trigs);

  proc_Red->SetLineWidth(3);
  proc_Blue->SetLineWidth(3);

  vector<shared_ptr<Process>> procs_Red = {proc_Red};
  vector<shared_ptr<Process>> procs_Blue = {proc_Blue};

  PlotMaker pm;
  string lepName[2] = {"el/", "mu/"};
  for(int i(0); i < 2; i++) {
    NamedFunc selection = baseline && lep.at(i) && pho && mass_cuts;
    pm.Push<Hist1D>(Axis(80, 100, 180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), selection, procs_Red, ops).Weight(wgt).Tag(lepName[i]+"llphoton_m");
    pm.Push<Hist1D>(Axis(100, 50, 150, "ll_m[llphoton_ill[0]]", "m_{ll} [GeV]", {}), selection, procs_Red, ops).Weight(wgt).Tag(lepName[i]+"ll_m");
    pm.Push<Hist1D>(Axis(80, 100, 180, kinFitH, "m_{ll#gamma} [GeV]", {}), selection, procs_Blue, ops).Weight(wgt).Tag(lepName[i]+"llphoton_m_Refit");
    pm.Push<Hist1D>(Axis(100, 50, 150, kinFit, "m_{ll} [GeV]", {}), selection, procs_Blue, ops).Weight(wgt).Tag(lepName[i]+"ll_m_Refit");
  }
  pm.min_print_ = true;
  pm.multithreaded_ = false;
  pm.MakePlots(35.9);
}

