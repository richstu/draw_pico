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
  Process::Type data =  Process::Type::data;
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
  auto proc_mc_back = Process::MakeShared<Baby_pico>("Pseudodata", data, 
                       kBlack, {mc_path+"*ZGToLLG*", mc_path+"*LLAJJ*",mc_path+"*DYJets*"},
                       trigs && "(type != 6200 || stitch_dy)");
  auto proc_hzg   = Process::MakeShared<Baby_pico>("HToZ#gamma(x20)", sig, 
                       kRed     ,{sig_path+"*.root"},   trigs);
  auto proc_hzg_vbf = Process::MakeShared<Baby_pico>("HToZ#gamma VBF(x1000)", sig, 
                       kMagenta                   ,{sig_path+"*VBF*.root"},   trigs);
  auto proc_hzg_gg  = Process::MakeShared<Baby_pico>("HToZ#gamma ggF(x1000)", sig, 
                       TColor::GetColor("#ff0000"),{sig_path+"*GluGlu*.root"},   trigs);
  proc_smzg->SetLineWidth (1); proc_dy->SetLineWidth   (1); proc_ewkzg->SetLineWidth(1);
  proc_hzg->SetLineWidth(3); proc_hzg_gg->SetLineWidth(3); proc_hzg_vbf->SetLineWidth(3);
  proc_mc_back->SetMarkerSize(1); vector<shared_ptr<Process>> procs = {proc_dy, proc_smzg, proc_hzg};
  vector<shared_ptr<Process>> cat_procs = {proc_mc_back, proc_hzg};
  vector<shared_ptr<Process>> theory_procs = {proc_dy, proc_smzg, proc_hzg_vbf, proc_hzg_gg};
  PlotOpt lin_lumi("txt/plot_styles.txt","CMSPaper");
  lin_lumi.Title(TitleType::info)
          .Stack(StackType::signal_overlay)
          .Overflow(OverflowType::none)
          .Bottom(BottomType::diff)
          .UseCMYK(false)
          .LegendColumns(1)
          .ShowBackgroundError(false)
          .FileExtensions({"pdf","root"});
  PlotOpt lin_stack = lin_lumi().Stack(StackType::signal_overlay);
//   vector<PlotOpt> ops = {lin_stack, lin_lumi};
  vector<PlotOpt> ops = {lin_lumi};
  NamedFunc pTt("pTt",[](const Baby &b) -> NamedFunc::ScalarType{
    TVector3 g = AssignGamma(b).Vect();
    TVector3 h = AssignH(b).Vect();
    TVector3 z = AssignZ(b).Vect();
    g.SetZ(0); h.SetZ(0); z.SetZ(0);
    return h.Cross((z-g).Unit()).Mag();
  });
  NamedFunc pTt2("pTt2",[](const Baby &b) -> NamedFunc::ScalarType{
    TVector3 g = AssignGamma(b).Vect();
    TVector3 h = AssignH(b).Vect();
    TVector3 z = AssignZ(b).Vect();
    g.SetZ(0); h.SetZ(0); z.SetZ(0);
    return h.Cross(z-g).Mag()/h.Mag();
  });
  NamedFunc el_check("el_check",[](const Baby &b) -> NamedFunc::ScalarType{
    if(abs(b.el_eta()->at(b.ll_i1()->at(0)) - b.el_eta()->at(b.ll_i2()->at(0))) < 0.075 &&
       abs(b.el_phi()->at(b.ll_i1()->at(0)) - b.el_phi()->at(b.ll_i2()->at(0))) < 0.125)
       return false;
    return true;
  });
  NamedFunc relpT("p_{T}^{#gamma}/m_{Z#gamma}",[](const Baby &b) -> NamedFunc::ScalarType{
    if(b.nphoton() == 1) return b.photon_pt()->at(0)/b.llphoton_m()->at(0);
    double relp(-1);
    if(b.photon_pt()->at(1) > b.photon_pt()->at(0)) 
      if(b.photon_drmin()->at(1) > 0.3) 
        for(int i = 0; i < b.nllphoton(); i++)
          if(b.llphoton_iph()->at(i) == 1) {
            relp = b.photon_pt()->at(1)/b.llphoton_m()->at(i);
            break;
          }
    if(relp < 0) relp = b.photon_pt()->at(0)/b.llphoton_m()->at(0);
    return relp;
  });
  NamedFunc llp_m("llp_m",[](const Baby &b) -> NamedFunc::ScalarType{
    if(b.nphoton() == 1) return b.llphoton_m()->at(0);
    double m(-1);
    if(b.photon_pt()->at(1) > b.photon_pt()->at(0)) 
      if(b.photon_drmin()->at(1) > 0.3) 
        for(size_t i = 0; i < b.llphoton_pt()->size(); i++)
          if(b.llphoton_iph()->at(i) == 1) {
            m = b.llphoton_m()->at(i);
            break;
          }
    if(m < 0) m = b.llphoton_m()->at(0);
    return m;
  });
  NamedFunc baseline("baseline",[](const Baby &b) -> NamedFunc::ScalarType{
    bool base(b.nphoton() > 0 && b.photon_id()->at(0) && b.photon_elveto()->at(0));
    vector<bool> bline = {b.photon_pt()->at(0)/b.llphoton_m()->at(0) > 0.12,
                          b.ll_m()->at(0) > 81.2 && b.ll_m()->at(0) < 101.2,
                          b.photon_drmin()->at(0) > 0.3};
//                           b.llphoton_m()->at(0) > 105 && b.llphoton_m()->at(0) < 160};
    for(size_t i = 0; i < bline.size(); i++)
      base = base && bline.at(i);
    return base;
  });
  NamedFunc els("electron channel",[](const Baby &b) -> NamedFunc::ScalarType{
    bool e =  b.ll_lepid()->at(0) == 11 && b.el_pt()->at(b.ll_i1()->at(0)) > 26 && b.el_pt()->at(b.ll_i2()->at(0)) > 17;
    return e;
  });
  NamedFunc mus("muon channel",[](const Baby &b) -> NamedFunc::ScalarType{
    bool mu = b.ll_lepid()->at(0) == 13 && b.mu_pt()->at(b.ll_i1()->at(0)) > 26 && b.mu_pt()->at(b.ll_i2()->at(0)) > 8;
    return mu;
  });
  NamedFunc vbf("VBF cut",[](const Baby &b) -> NamedFunc::ScalarType{
    return b.vbf_mva() > -0.99;
  });
  NamedFunc leps("e-or-#mu",[](const Baby &b) -> NamedFunc::ScalarType{
    bool e =  b.ll_lepid()->at(0) == 11 && b.el_pt()->at(b.ll_i1()->at(0)) > 26 && b.el_pt()->at(b.ll_i2()->at(0)) > 17;
    bool mu = b.ll_lepid()->at(0) == 13 && b.mu_pt()->at(b.ll_i1()->at(0)) > 26 && b.mu_pt()->at(b.ll_i2()->at(0)) > 8;
    return e || mu;
  });
  vector<NamedFunc> cat = { baseline && vbf && leps,
                            baseline && relpT > 0.4 && leps && !vbf               ,
                            baseline &&   pTt > 40  &&  els && !vbf && relpT < 0.4,
                            baseline &&   pTt < 40  &&  els && !vbf && relpT < 0.4,
                            baseline &&   pTt > 40  &&  mus && !vbf && relpT < 0.4,
                            baseline &&   pTt < 40  &&  mus && !vbf && relpT < 0.4};
  NamedFunc wgt("w_lumi",[](const Baby &b) -> NamedFunc::ScalarType{ 
    double weight = b.w_lumi();
    if(b.type() >= 200000 && b.type() <= 205000)
      return weight;
    else if(b.type() == 6200)
      return 0.823*weight/1.664;
    return   0.823*weight;
  });
  PlotMaker pm;
  vector<NamedFunc> lep = {els, mus};
//   for(int i(0); i < 2; i++) {
//       NamedFunc cut = lep.at(i) && baseline;
//       pm.Push<Hist1D>(Axis(50,0,100, pTt2, "Alternative p_{Tt} [GeV]",{40}), cut, procs, ops).Weight(wgt);
//       pm.Push<Hist1D>(Axis(40,0,1, "photon_pt[0]/llphoton_m[0]", "p_{T}^{#gamma}/m_{Z#gamma}",{0.12,0.4}), cut, procs, ops).Weight(wgt);
//       pm.Push<Hist1D>(Axis(40,0,1, relpT, "modified p_{T}^{#gamma}/m_{Z#gamma}",{0.12,0.4}), cut, procs, ops).Weight(wgt);
//   }
//   pm.Push<Hist1D>(Axis(20,100,180, "llphoton_m[0]", "m_{Z#gamma} [GeV]", {110,155}), cat.at(0), cat_procs, ops).Weight(wgt).Tag("cat1");
//   pm.Push<Hist1D>(Axis(40,100,180, "llphoton_m[0]", "m_{Z#gamma} [GeV]", {105,155}), cat.at(1), cat_procs, ops).Weight(wgt).Tag("cat2");
//   pm.Push<Hist1D>(Axis(20,100,180, "llphoton_m[0]", "m_{Z#gamma} [GeV]", {115,145}), cat.at(2), cat_procs, ops).Weight(wgt).Tag("cat3");
//   pm.Push<Hist1D>(Axis(40,100,180, "llphoton_m[0]", "m_{Z#gamma} [GeV]", {125.3}), cat.at(3), cat_procs, ops).Weight(wgt).Tag("cat4");
//   pm.Push<Hist1D>(Axis(20,100,180, "llphoton_m[0]", "m_{Z#gamma} [GeV]", {115,160}), cat.at(4), cat_procs, ops).Weight(wgt).Tag("cat5");
//   pm.Push<Hist1D>(Axis(40,100,180, "llphoton_m[0]", "m_{Z#gamma} [GeV]", {125.3}), cat.at(5), cat_procs, ops).Weight(wgt).Tag("cat6");
  pm.Push<Hist1D>(Axis(20,100,180, "llphoton_m[0]", "m_{Z#gamma} [GeV]", {}), cat.at(0), cat_procs, ops).Weight(wgt).Tag("cat1");
  pm.Push<Hist1D>(Axis(40,100,180, "llphoton_m[0]", "m_{Z#gamma} [GeV]", {}), cat.at(1), cat_procs, ops).Weight(wgt).Tag("cat2");
  pm.Push<Hist1D>(Axis(20,100,180, "llphoton_m[0]", "m_{Z#gamma} [GeV]", {}), cat.at(2), cat_procs, ops).Weight(wgt).Tag("cat3");
  pm.Push<Hist1D>(Axis(40,100,180, "llphoton_m[0]", "m_{Z#gamma} [GeV]", {}), cat.at(3), cat_procs, ops).Weight(wgt).Tag("cat4");
  pm.Push<Hist1D>(Axis(20,100,180, "llphoton_m[0]", "m_{Z#gamma} [GeV]", {}), cat.at(4), cat_procs, ops).Weight(wgt).Tag("cat5");
  pm.Push<Hist1D>(Axis(40,100,180, "llphoton_m[0]", "m_{Z#gamma} [GeV]", {}), cat.at(5), cat_procs, ops).Weight(wgt).Tag("cat6");
//   pm.Push<Hist1D>(Axis(40,100,180, "llphoton_m[0]", "m_{Z#gamma} [GeV]", {115,170}), baseline && !vbf && relpT < 0.4, cat_procs, ops).Weight(wgt).Tag("catuntagged");
//   pm.Push<Hist1D>(Axis(55,115,170, "llphoton_m[0]", "m_{Z#gamma} [GeV]", {}), cat.at(3), cat_procs, ops).Weight(wgt).Tag("cms-cat4");
//   pm.Push<Hist1D>(Axis(55,115,170, "llphoton_m[0]", "m_{Z#gamma} [GeV]", {}), cat.at(5), cat_procs, ops).Weight(wgt).Tag("cms-cat6");
  pm.min_print_ = true;
  pm.MakePlots(36);
}

