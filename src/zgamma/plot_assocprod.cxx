#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <bitset>
#include <unistd.h>
#include <getopt.h>
#include "TError.h"
#include "TColor.h"
#include "TH1.h"
#include "TH2.h"
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

int main() {

  //setup
  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");
  Process::Type back =  Process::Type::background;
  //Process::Type data =  Process::Type::data;
  Process::Type sig =  Process::Type::signal;

  string bfolder("/net/cms17/cms17r0/pico/");
  string mc_path( bfolder+"NanoAODv9/htozgamma_deathvalley_v2/2017/mc/skim_llg/");
  string sig_path(bfolder+"NanoAODv2/zgamma_signal/2017/signal/unskimmed/");

  //NamedFunc el_trigs = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL||HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ||HLT_Ele35_WPTight_Gsf||HLT_Ele32_WPTight_Gsf_L1DoubleEG";
  //NamedFunc mu_trigs = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8||HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8||HLT_IsoMu27||HLT_IsoMu24";
  NamedFunc el_trigs = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL||HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ";
  NamedFunc mu_trigs = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ";
  //NamedFunc el_trigs("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
  //NamedFunc mu_trigs("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8");
  NamedFunc trigs(el_trigs || mu_trigs);

  //NamedFunc to remove events with no ll or llphoton objects, which can happen in prodcutions from prasanna's earlier branch with ll sign req.
  NamedFunc sanitize("sanitize",[](const Baby &b) -> NamedFunc::ScalarType{ 
    unsigned ll_size = b.ll_m()->size();
    unsigned llg_size = b.llphoton_m()->size();
    if (b.nlep()>=2 && b.nphoton()>=1 && (ll_size == 0 || llg_size == 0)) return 0;
    return 1;
  });

  //note: need to modify for ZH
  NamedFunc z_to_ll("z_to_ll",[](const Baby &b) -> NamedFunc::ScalarType{ 
    for (unsigned imc = 0; imc < b.mc_pt()->size(); imc++) {
      if ((abs(b.mc_id()->at(imc))==11 || abs(b.mc_id()->at(imc))==13)&&b.mc_mom()->at(imc)==23) {
        return 1;
      }
    }
    return 0;
  });

  auto proc_smzg        = Process::MakeShared<Baby_pico>("Z+#gamma",       back, 
                             TColor::GetColor("#16bac5"),{mc_path+"*ZGToLLG*",mc_path+"*EWKZ2Jets*"}, trigs);
  auto proc_dy          = Process::MakeShared<Baby_pico>("Z+Fake Photon",               back, 
                             TColor::GetColor("#ffb400"),{mc_path+"*DYJets*"},  trigs && "stitch_dy");
  //poor man's stitch
  auto proc_tt          = Process::MakeShared<Baby_pico>("tt",               back, 
                             TColor::GetColor("#ff6600"),{mc_path+"*TTTo2L2Nu*"},  trigs && "photon_pflavor[0]!=1");
  auto proc_ttg         = Process::MakeShared<Baby_pico>("t/tt+#gamma",               back, 
                             TColor::GetColor("#ffa366"),{mc_path+"*TGJets*"},  trigs);
  auto proc_ttz         = Process::MakeShared<Baby_pico>("ttZ",               back, 
                             TColor::GetColor("#ff0080"),{mc_path+"*ttZJets*"},  trigs);
  auto proc_ttw         = Process::MakeShared<Baby_pico>("ttW",               back, 
                             TColor::GetColor("#993d00"),{mc_path+"*ttWJets_Tune*"},  trigs);
  //poor man's stitch for now
  auto proc_vv          = Process::MakeShared<Baby_pico>("Multiboson+Fake Photon",          back, 
                             TColor::GetColor("#5500ff"),{mc_path+"*WW_Tune*",mc_path+"*WZ_Tune*",mc_path+"*ZZ_Tune*",
                             mc_path+"*WWW_4F_Tune*",mc_path+"*WZZ*",mc_path+"*ZZZ*"},  trigs && "photon_pflavor[0]!=1");
  auto proc_wwz          = Process::MakeShared<Baby_pico>("Multiboson+Fake Photon (WWZ)",          back, 
                             TColor::GetColor("#5500ff"),{mc_path+"*WWZ_4F_Tune*"},  trigs && "photon_pflavor[0]!=1");
  auto proc_vvg         = Process::MakeShared<Baby_pico>("Multiboson+#gamma",               back, 
                             TColor::GetColor("#9966ff"),{mc_path+"*WWG*",mc_path+"*WZG*",mc_path+"*ZZG*"},  trigs );
  auto proc_higgs       = Process::MakeShared<Baby_pico>("Other Higgs",               back, 
                             TColor::GetColor("#29a329"),{mc_path+"*GluGluH*"},  trigs);

  auto proc_hzg_tth_x1000     = Process::MakeShared<Baby_pico>("ttH H#rightarrow Z#gamma (x1000)", sig, 
                           kRed     ,{sig_path+"*ttHToZG*"}, sanitize&&trigs&&z_to_ll);  
  
  auto proc_hzg_tth     = Process::MakeShared<Baby_pico>("ttH H#rightarrow Z#gamma (x50)", sig, 
                           kRed     ,{sig_path+"*ttHToZG*"}, sanitize&&trigs&&z_to_ll);  
  auto proc_hzg_other   = Process::MakeShared<Baby_pico>("Other H#rightarrow Z#gamma (x50)", sig, 
                           kOrange  ,{sig_path+"*GluGluHToZG*",sig_path+"*VBFHToZG*",sig_path+"*WplusH*",sig_path+"*WminusH*",sig_path+"*ZH*"},   sanitize&&trigs&&z_to_ll);

  auto proc_hzg_tth_full = Process::MakeShared<Baby_pico>("ttH H#rightarrow Z#gamma (x50)", sig, 
                           kRed     ,{sig_path+"*ttHToZG*"}, sanitize&&z_to_ll);  

  //proc_smzg->SetLineWidth(1); proc_dy->SetLineWidth(1); 
  //proc_mc_back->SetMarkerSize(1); 
  proc_hzg_tth->SetLineWidth(3); proc_hzg_other->SetLineWidth(3); 
  vector<shared_ptr<Process>> procs = {proc_dy, proc_smzg, proc_tt, proc_ttg, proc_ttz, proc_ttw, proc_vv, proc_wwz, proc_vvg, proc_higgs, proc_hzg_tth, proc_hzg_other};
  vector<shared_ptr<Process>> procs_mult = {proc_dy, proc_smzg, proc_tt, proc_ttg, proc_ttz, proc_ttw, proc_vv, proc_vvg, proc_higgs, proc_hzg_tth_x1000};
  vector<shared_ptr<Process>> procs_sig_full = {proc_hzg_tth_full};

  PlotOpt lin_lumi("txt/plot_styles.txt","CMSPaper");
  lin_lumi.Title(TitleType::info)
          .Stack(StackType::signal_overlay)
          .Overflow(OverflowType::none)
          .Bottom(BottomType::off)
          .UseCMYK(false)
          .LegendColumns(1)
          .ShowBackgroundError(false)
          .FileExtensions({"pdf"});
  PlotOpt lin_sorb = lin_lumi;
  lin_sorb.Bottom(BottomType::sorb).FileExtensions({"pdf","root"});
  PlotOpt lin_sorb_upper = lin_lumi;
  lin_sorb_upper.Bottom(BottomType::sorb_cut_upper).FileExtensions({"pdf","root"});
  PlotOpt lin_shapes = lin_lumi;
  lin_shapes.Stack(StackType::shapes);
  PlotOpt log_lumi = lin_lumi;
  log_lumi.YAxis(YAxisType::log);
  vector<PlotOpt> ops_shapes = {lin_shapes};
  vector<PlotOpt> ops = {lin_lumi};
  vector<PlotOpt> ops_log = {log_lumi};
  vector<PlotOpt> ops_sorb = {lin_sorb};
  vector<PlotOpt> ops_sorb_upper = {lin_sorb_upper};
  PlotOpt twodim_lin_lumi("txt/plot_styles.txt", "Scatter");
  twodim_lin_lumi.Title(TitleType::info)
                 .YAxis(YAxisType::log)
                 .Overflow(OverflowType::overflow)
                 .LogMinimum(0.001);
  vector<PlotOpt> twodim_ops = {twodim_lin_lumi};

  //NamedFuncs
  //NamedFunc pTt("pTt",[](const Baby &b) -> NamedFunc::ScalarType{
  //  TVector3 g = AssignGamma(b).Vect();
  //  TVector3 h = AssignH(b).Vect();
  //  TVector3 z = AssignZ(b).Vect();
  //  g.SetZ(0); h.SetZ(0); z.SetZ(0);
  //  return h.Cross((z-g).Unit()).Mag();
  //});
  //NamedFunc pTt2("pTt2",[](const Baby &b) -> NamedFunc::ScalarType{
  //  TVector3 g = AssignGamma(b).Vect();
  //  TVector3 h = AssignH(b).Vect();
  //  TVector3 z = AssignZ(b).Vect();
  //  g.SetZ(0); h.SetZ(0); z.SetZ(0);
  //  return h.Cross(z-g).Mag()/h.Mag();
  //});
  //NamedFunc el_check("el_check",[](const Baby &b) -> NamedFunc::ScalarType{
  //  if(abs(b.el_eta()->at(b.ll_i1()->at(0)) - b.el_eta()->at(b.ll_i2()->at(0))) < 0.075 &&
  //     abs(b.el_phi()->at(b.ll_i1()->at(0)) - b.el_phi()->at(b.ll_i2()->at(0))) < 0.125)
  //     return false;
  //  return true;
  //});
  //NamedFunc relpT("p_{T}^{#gamma}/m_{Z#gamma}",[](const Baby &b) -> NamedFunc::ScalarType{
  //  if(b.nphoton() == 1) return b.photon_pt()->at(0)/b.llphoton_m()->at(0);
  //  double relp(-1);
  //  if(b.photon_pt()->at(1) > b.photon_pt()->at(0)) 
  //    if(b.photon_drmin()->at(1) > 0.3) 
  //      for(int i = 0; i < b.nllphoton(); i++)
  //        if(b.llphoton_iph()->at(i) == 1) {
  //          relp = b.photon_pt()->at(1)/b.llphoton_m()->at(i);
  //          break;
  //        }
  //  if(relp < 0) relp = b.photon_pt()->at(0)/b.llphoton_m()->at(0);
  //  return relp;
  //});
  //NamedFunc llp_m("llp_m",[](const Baby &b) -> NamedFunc::ScalarType{
  //  if(b.nphoton() == 1) return b.llphoton_m()->at(0);
  //  double m(-1);
  //  if(b.photon_pt()->at(1) > b.photon_pt()->at(0)) 
  //    if(b.photon_drmin()->at(1) > 0.3) 
  //      for(size_t i = 0; i < b.llphoton_pt()->size(); i++)
  //        if(b.llphoton_iph()->at(i) == 1) {
  //          m = b.llphoton_m()->at(i);
  //          break;
  //        }
  //  if(m < 0) m = b.llphoton_m()->at(0);
  //  return m;
  //});
  //NamedFunc els("electron channel",[](const Baby &b) -> NamedFunc::ScalarType{
  //  bool e =  b.ll_lepid()->at(0) == 11 && b.el_pt()->at(b.ll_i1()->at(0)) > 26 && b.el_pt()->at(b.ll_i2()->at(0)) > 17;
  //  return e;
  //});
  //NamedFunc mus("muon channel",[](const Baby &b) -> NamedFunc::ScalarType{
  //  bool mu = b.ll_lepid()->at(0) == 13 && b.mu_pt()->at(b.ll_i1()->at(0)) > 26 && b.mu_pt()->at(b.ll_i2()->at(0)) > 8;
  //  return mu;
  //});
  //NamedFunc vbf("VBF cut",[](const Baby &b) -> NamedFunc::ScalarType{
  //  return b.vbf_mva() > -0.99;
  //});
  //NamedFunc leps("e-or-#mu",[](const Baby &b) -> NamedFunc::ScalarType{
  //  bool e =  b.ll_lepid()->at(0) == 11 && b.el_pt()->at(b.ll_i1()->at(0)) > 26 && b.el_pt()->at(b.ll_i2()->at(0)) > 17;
  //  bool mu = b.ll_lepid()->at(0) == 13 && b.mu_pt()->at(b.ll_i1()->at(0)) > 26 && b.mu_pt()->at(b.ll_i2()->at(0)) > 8;
  //  return e || mu;
  //});
  //NamedFunc wgt("weight",[](const Baby &b) -> NamedFunc::ScalarType{ 
  //  double weight = b.w_lumi();
  //  if(b.type() >= 200000 && b.type() <= 205000)
  //    return weight;
  //  else if(b.type() == 6200)
  //    return 0.823*weight/1.664;
  //  return   0.823*weight;
  //});
  
  NamedFunc mc_higgs_pt("mc_higgs_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
      if (b.mc_id()->at(imc)==25) {
        return b.mc_pt()->at(imc);
      }
    }
    return -999;
  });

  NamedFunc mc_photon_pt("mc_photon_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
      if (b.mc_id()->at(imc)==22 && b.mc_mom()->at(imc)==25) {
        return b.mc_pt()->at(imc);
      }
    }
    return -999;
  });

  NamedFunc mc_deltar_zg("mc_deltar_zg",[](const Baby &b) -> NamedFunc::ScalarType{
    TVector3 ph;
    TVector3 Z;
    for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
      if (b.mc_id()->at(imc)==22 && b.mc_mom()->at(imc)==25) {
        ph.SetPtEtaPhi(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc));
      }
      if (b.mc_id()->at(imc)==23 && b.mc_mom()->at(imc)==25) {
        Z.SetPtEtaPhi(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc));
      }
    }
    return ph.DeltaR(Z);
  });

  NamedFunc mc_llg_in_acceptance("mc_llg_in_acceptance",[](const Baby &b) -> NamedFunc::ScalarType{
    for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
      if (b.mc_id()->at(imc)==22 && b.mc_mom()->at(imc)==25) {
        if (!(b.mc_pt()->at(imc)>15 && abs(b.mc_eta()->at(imc))<2.5))
          return 0;
      }
      if ((abs(b.mc_id()->at(imc))==11 || abs(b.mc_id()->at(imc))==13)&&b.mc_mom()->at(imc)==23) {
        int mom_idx = b.mc_momidx()->at(imc);
        if (mom_idx != -1) {
          if (b.mc_mom()->at(mom_idx)==25) {
            float eta_cut = 2.5;
            float pt_cut = 15;
            if (abs(b.mc_id()->at(imc))==13) {
              eta_cut = 2.4;
              pt_cut = 10;
            }
            if (!(b.mc_pt()->at(imc)>pt_cut && abs(b.mc_eta()->at(imc))<eta_cut))
              return 0;
          }
        }
      }
    }
    return 1;
  });

  NamedFunc mc_llg_reco("mc_llg_reco",[](const Baby &b) -> NamedFunc::ScalarType{
    TVector3 l1, l2, ph;
    TVector3 reco_obj;
    int z_decay_pdgid = 11;
    int lep_idx = 0;
    bool found_ph(false), found_l1(false), found_l2(false);
    for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
      if (b.mc_id()->at(imc)==22 && b.mc_mom()->at(imc)==25) {
        ph.SetPtEtaPhi(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc));
      }
      if ((abs(b.mc_id()->at(imc))==11 || abs(b.mc_id()->at(imc))==13)&&b.mc_mom()->at(imc)==23) {
        int mom_idx = b.mc_momidx()->at(imc);
        if (mom_idx != -1) {
          if (b.mc_mom()->at(mom_idx)==25) {
            z_decay_pdgid = abs(b.mc_id()->at(imc));
            if (lep_idx==0) {
              l1.SetPtEtaPhi(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc));
              lep_idx = 1;
            }
            else {
              l2.SetPtEtaPhi(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc));
            }
          }
        }
      }
    }
    for (unsigned iph = 0; iph < b.photon_pt()->size(); iph++) {
      reco_obj.SetPtEtaPhi(b.photon_pt()->at(iph),b.photon_eta()->at(iph),b.photon_phi()->at(iph));
      if (reco_obj.DeltaR(ph)<0.2) {
        found_ph = true;
      }
    }
    if (z_decay_pdgid==13) {
      for (unsigned imu = 0; imu < b.mu_pt()->size(); imu++) {
        reco_obj.SetPtEtaPhi(b.mu_pt()->at(imu),b.mu_eta()->at(imu),b.mu_phi()->at(imu));
        if (reco_obj.DeltaR(l1)<0.2) {
          found_l1 = true;
        }
        else if (reco_obj.DeltaR(l2)<0.2) {
          found_l2 = true;
        }
      }
    }
    else {
      for (unsigned iel = 0; iel < b.el_pt()->size(); iel++) {
        reco_obj.SetPtEtaPhi(b.el_pt()->at(iel),b.el_eta()->at(iel),b.el_phi()->at(iel));
        if (reco_obj.DeltaR(l1)<0.2) {
          found_l1 = true;
        }
        else if (reco_obj.DeltaR(l2)<0.2) {
          found_l2 = true;
        }
      }
    }
    if (found_ph&&found_l1&&found_l2) return 1;
    return 0;
  });

  NamedFunc sigphoton_mindr("sigphoton_mindr",[](const Baby &b) -> NamedFunc::ScalarType{
    for (unsigned iph = 0; iph < b.photon_pt()->size(); iph++) {
      if (b.photon_pt()->at(iph)<15) continue;
      if (abs(b.photon_eta()->at(iph))>2.5) continue;
      if (abs(b.photon_eta()->at(iph))<1.4442 && b.photon_idmva()->at(iph) < -0.4) continue;
      if (abs(b.photon_eta()->at(iph))>1.4442 && abs(b.photon_eta()->at(iph))<1.566) continue;
      if (abs(b.photon_eta()->at(iph))>1.566 && b.photon_idmva()->at(iph) < -0.58) continue;
      if (!b.photon_elveto()->at(iph)) continue;
      //use first signal photon
      TVector3 ph, lep;
      ph.SetPtEtaPhi(b.photon_pt()->at(iph),b.photon_eta()->at(iph),b.photon_phi()->at(iph));
      float mindr = 999;
      for (unsigned iel = 0; iel < b.el_pt()->size(); iel++) {
        if (b.el_sig()->at(iel)) {
          lep.SetPtEtaPhi(b.el_pt()->at(iel),b.el_eta()->at(iel),b.el_phi()->at(iel));
          if (ph.DeltaR(lep)<mindr) mindr = ph.DeltaR(lep);
        }
      }
      for (unsigned imu = 0; imu < b.mu_pt()->size(); imu++) {
        if (b.mu_sig()->at(imu)) {
          lep.SetPtEtaPhi(b.mu_pt()->at(imu),b.mu_eta()->at(imu),b.mu_phi()->at(imu));
          if (ph.DeltaR(lep)<mindr) mindr = ph.DeltaR(lep);
        }
      }
      return mindr;
    }
    return 999;
  });

  NamedFunc coslowertheta("costheta",[](const Baby &b) -> NamedFunc::ScalarType{
    return cos_theta(b);
  });

  NamedFunc coscaptheta("coscaptheta",[](const Baby &b) -> NamedFunc::ScalarType{
    return cos_Theta(b);
  });

  NamedFunc llgphi("phi",[](const Baby &b) -> NamedFunc::ScalarType{
    return Getphi(b);
  });

  NamedFunc wgt("weight",[](const Baby &b) -> NamedFunc::ScalarType{
    double w_lumi = b.w_lumi();
    //fix mis-weighted assoc. prod
    if(b.type() >= 200200 && b.type() <= 200500)
      w_lumi /= 0.100974;
    //fix mis-weighted ttW and WWZ
    if (b.FirstFileName().find("ttW") != std::string::npos)
      w_lumi = 610.5/27662138;
    if (b.FirstFileName().find("WWZ") != std::string::npos)
      w_lumi = 354.0/10116400;
    return w_lumi;
  });

  //multiplicative factor for signal
  NamedFunc sigx50("sigx50",[](const Baby &b) -> NamedFunc::ScalarType{
    double w_lumi = 1.0;
    if(b.type() >= 200000 && b.type() <= 200500)
      w_lumi *= 50.0;
    return w_lumi;
  });

  NamedFunc sigx500("sigx500",[](const Baby &b) -> NamedFunc::ScalarType{
    double w_lumi = 1.0;
    if(b.type() >= 200000 && b.type() <= 200500)
      w_lumi *= 500.0;
    return w_lumi;
  });

  NamedFunc higgs_pt("higgs_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    if(b.type() >= 200000 && b.type() <= 205000)
      return b.w_lumi();
    return b.w_lumi();
  });

  NamedFunc decorr_photon_pt("decorr_photon_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.photon_pt()->size()==0) return 0;
    if (b.llphoton_m()->size()==0) return 0;
    //return (b.photon_pt()->at(0))-0.202*(b.llphoton_m()->at(0));
    return (b.photon_pt()->at(0))-0.207*(b.llphoton_m()->at(0));
  });

  NamedFunc photon_pt_mass("photon_pt_mass",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.photon_pt()->size()==0) return 0;
    if (b.llphoton_m()->size()==0) return 0;
    //return (b.photon_pt()->at(0))-0.202*(b.llphoton_m()->at(0));
    return (b.photon_pt()->at(0))/(b.llphoton_m()->at(0));
  });

  NamedFunc truth_mllg("truth_mllg",[](const Baby &b) -> NamedFunc::ScalarType{
    //calculate truth mllg for dy (pre-stitch) events with ISR
    int lep_num = 0;
    float max_ph_pt = 0;
    TLorentzVector lep1, lep2, ph;
    for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
      if (abs(b.mc_id()->at(imc))==11) {
        if (lep_num==0)
          lep1.SetPtEtaPhiM(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc),0.000511);
        else
          lep2.SetPtEtaPhiM(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc),0.000511);
        lep_num++;
      }
      if (abs(b.mc_id()->at(imc))==13) {
        if (lep_num==0)
          lep1.SetPtEtaPhiM(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc),0.106);
        else
          lep2.SetPtEtaPhiM(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc),0.106);
        lep_num++;
      }
      if (b.mc_id()->at(imc)==22) {
        if (b.mc_pt()->at(imc) > max_ph_pt) {
          max_ph_pt = b.mc_pt()->at(imc);
          ph.SetPtEtaPhiM(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc),0.0);
        }
      }
    }
    if (lep_num<2 || max_ph_pt <= 0) return -999;
    return (lep1+lep2+ph).M();
  });

  NamedFunc truth_photon_pt("truth_photon_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    float max_ph_pt = 0;
    for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
      if (b.mc_id()->at(imc)==22) {
        if (b.mc_pt()->at(imc) > max_ph_pt) {
          max_ph_pt = b.mc_pt()->at(imc);
        }
      }
    }
    return max_ph_pt;
  });

  NamedFunc truth_photon_abseta("truth_photon_abseta",[](const Baby &b) -> NamedFunc::ScalarType{
    float max_ph_pt = 0;
    float ph_abseta = 0;
    for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
      if (b.mc_id()->at(imc)==22) {
        if (b.mc_pt()->at(imc) > max_ph_pt) {
          max_ph_pt = b.mc_pt()->at(imc);
          ph_abseta = abs(b.mc_eta()->at(imc)); 
        }
      }
    }
    return ph_abseta;
  });

  NamedFunc stitch_zg("stitch_zg",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type()==6200) return b.stitch_dy();
    return 1;
  });

  NamedFunc truth_lepton_minpt("truth_lepton_minpt",[](const Baby &b) -> NamedFunc::ScalarType{
    float min_lep_pt = 999;
    for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
      if (abs(b.mc_id()->at(imc))==11||abs(b.mc_id()->at(imc))==13) {
        if (b.mc_pt()->at(imc) < min_lep_pt) {
          min_lep_pt = b.mc_pt()->at(imc);
        }
      }
    }
    return min_lep_pt;
  });

  NamedFunc truth_lepton_eta_acc("truth_lepton_eta_acc",[](const Baby &b) -> NamedFunc::ScalarType{
    for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
      if (abs(b.mc_id()->at(imc))==11||abs(b.mc_id()->at(imc))==13) {
        if (abs(b.mc_eta()->at(imc))>2.4)
          return 0.0;
      }
    }
    return 1.0;
  });

  NamedFunc mc_photon_drmin("mc_photon_drmin",[](const Baby &b) -> NamedFunc::ScalarType{
    float ph_eta = -999;
    float ph_phi = 0;
    float min_dr = 999;
    float dr = 999;
    bool found_ph = false;
    for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
      if (b.mc_id()->at(imc)==22 && b.mc_mom()->at(imc)==25) {
        ph_eta = b.mc_eta()->at(imc);
        ph_phi = b.mc_phi()->at(imc);
        found_ph = true;
      }
    }
    if (!found_ph) return 999;
    for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
      if (abs(b.mc_id()->at(imc))==11||abs(b.mc_id()->at(imc))==13) {
        dr = deltaR(ph_eta, ph_phi, b.mc_eta()->at(imc), b.mc_phi()->at(imc));
        if (dr < min_dr)
          min_dr = dr;
      }
    }
    return min_dr;
  });

  NamedFunc photon_sig_nodr("photon_sig_nodr",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> photon_sig_nodr_;
    for (unsigned iph = 0; iph < b.photon_pt()->size(); iph++) {
      float abseta = abs(b.photon_eta()->at(iph));
      float idmva = b.photon_idmva()->at(iph);
      float pt = b.photon_pt()->at(iph);
      if (pt>15 && ((abseta<1.4442&&idmva>-0.4)||(abseta>1.566&&abseta<2.5&&idmva>-0.58)))
        photon_sig_nodr_.push_back(1);
      else
        photon_sig_nodr_.push_back(0);
    }
    return photon_sig_nodr_;
  });

  NamedFunc nphoton_nodr("nphoton_nodr",[photon_sig_nodr](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> photon_sig_nodr_ = photon_sig_nodr.GetVector(b);
    int nph = 0;
    for (unsigned iph = 0; iph < b.photon_pt()->size(); iph++) {
      if (photon_sig_nodr_[iph]>0)
        nph++;
    }
    return nph;
  });

  NamedFunc true_photon_recod("true_photon_recod",[photon_sig_nodr](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> photon_sig_nodr_ = photon_sig_nodr.GetVector(b);
    for (unsigned iph = 0; iph < b.photon_pt()->size(); iph++) {
      if (photon_sig_nodr_[iph]>0) {
        if (b.photon_pflavor()->at(iph)==1)
          return 1;
      }
    }
    return 0;
  });

  NamedFunc leps_recod("leps_recod",[](const Baby &b) -> NamedFunc::ScalarType{
    int reco_num = 0;
    float true_lep_eta = 999;
    float true_lep_phi = 0;
    float true_antilep_eta = 999;
    float true_antilep_phi = 0;
    for (unsigned imc = 0; imc < b.mc_pt()->size(); imc++) {
      int pdgid = b.mc_id()->at(imc);
      if ((pdgid == 11 || pdgid == 13) && b.mc_mom()->at(imc)==23) {
        true_lep_eta = b.mc_eta()->at(imc);
        true_lep_phi= b.mc_phi()->at(imc);
      }
      if ((pdgid == -11 || pdgid == -13) && b.mc_mom()->at(imc)==23) {
        true_antilep_eta = b.mc_eta()->at(imc);
        true_antilep_phi= b.mc_phi()->at(imc);
      }
    }
    for (unsigned iel = 0; iel < b.el_pt()->size(); iel++) {
      if (b.el_sig()->at(iel)) {
        float dr_lep = deltaR(b.el_eta()->at(iel), b.el_phi()->at(iel), true_lep_eta, true_lep_phi);
        float dr_antilep = deltaR(b.el_eta()->at(iel), b.el_phi()->at(iel), true_antilep_eta, true_antilep_phi);
        if (dr_lep < 0.3 || dr_antilep < 0.3)
          reco_num++;
      }
    }
    for (unsigned imu = 0; imu < b.mu_pt()->size(); imu++) {
      if (b.mu_sig()->at(imu)) {
        float dr_lep = deltaR(b.mu_eta()->at(imu), b.mu_phi()->at(imu), true_lep_eta, true_lep_phi);
        float dr_antilep = deltaR(b.mu_eta()->at(imu), b.mu_phi()->at(imu), true_antilep_eta, true_antilep_phi);
        if (dr_lep < 0.3 || dr_antilep < 0.3)
          reco_num++;
      }
    }
    return reco_num;
  });

  NamedFunc min_electron_pt("min_electron_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    float min_electron_pt_ = 999;
    for (unsigned iel = 0; iel < b.el_pt()->size(); iel++) {
      if (b.el_sig()->at(iel)) {
        if (b.el_pt()->at(iel)<min_electron_pt_)
          min_electron_pt_ = b.el_pt()->at(iel);
      }
    }
    return min_electron_pt_;
  });

  NamedFunc min_muon_pt("min_muon_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    float min_muon_pt_ = 999;
    for (unsigned imu = 0; imu < b.mu_pt()->size(); imu++) {
      if (b.mu_sig()->at(imu)) {
        if (b.mu_pt()->at(imu)<min_muon_pt_)
          min_muon_pt_ = b.mu_pt()->at(imu);
      }
    }
    return min_muon_pt_;
  });

  NamedFunc assoc_electron_pt("assoc_electron_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    float min_electron_pt_ = 999;
    std::vector<int> skip_indices;
    if (b.ll_m()->size()>0) {
      if (b.ll_lepid()->at(0)==11) {
        skip_indices.push_back(b.ll_i1()->at(0));
        skip_indices.push_back(b.ll_i2()->at(0));
      }
    }
    if (skip_indices.size()==0) {
      skip_indices.push_back(-1);
      skip_indices.push_back(-1);
    }
    for (unsigned iel = 0; iel < b.el_pt()->size(); iel++) {
      if (b.el_sig()->at(iel) && static_cast<int>(iel) != skip_indices[0] && static_cast<int>(iel) != skip_indices[1]) {
        if (b.el_pt()->at(iel)<min_electron_pt_)
          min_electron_pt_ = b.el_pt()->at(iel);
      }
    }
    return min_electron_pt_;
  });

  NamedFunc assoc_muon_pt("assoc_muon_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    float min_muon_pt_ = 999;
    std::vector<int> skip_indices;
    if (b.ll_m()->size()>0) {
      if (b.ll_lepid()->at(0)==13) {
        skip_indices.push_back(b.ll_i1()->at(0));
        skip_indices.push_back(b.ll_i2()->at(0));
      }
    }
    if (skip_indices.size()==0) {
      skip_indices.push_back(-1);
      skip_indices.push_back(-1);
    }
    for (unsigned imu = 0; imu < b.mu_pt()->size(); imu++) {
      if (b.mu_sig()->at(imu) && static_cast<int>(imu) != skip_indices[0] && static_cast<int>(imu) != skip_indices[1]) {
        if (b.mu_pt()->at(imu)<min_muon_pt_)
          min_muon_pt_ = b.mu_pt()->at(imu);
      }
    }
    return min_muon_pt_;
  });

  NamedFunc photon_res = "photon_pterr[0]/photon_pt[0]";

  NamedFunc baseline_nolep = "(ll_m[0]>50) && (photon_pt[0]/llphoton_m[0]>=15.0/110.0) && ((llphoton_m[0]+ll_m[0])>=185) && (photon_drmin[0]>0.4)";
  NamedFunc baseline_el = "(ll_lepid[0]==11) && (el_pt[ll_i1[0]]>25) && (el_pt[ll_i2[0]]>15)" && baseline_nolep;
  NamedFunc baseline_mu = "(ll_lepid[0]==13) && (mu_pt[ll_i1[0]]>20) && (mu_pt[ll_i2[0]]>10)" && baseline_nolep;
  baseline_el.Name("electron_baseline");
  baseline_mu.Name("muon_baseline");
  NamedFunc baseline = (baseline_el||baseline_mu);
  baseline.Name("baseline");
  std::vector<NamedFunc> baseline_lep = {baseline_el, baseline_mu, baseline_el||baseline_mu};

  NamedFunc higgs_window = "llphoton_m[0]>122&&llphoton_m[0]<128";
  NamedFunc photon_pt = "photon_pt[0]";
  std::vector<float> photon_pt_bins = {15,20,25,30,40,50,200};

  //note skim_llg contains nlep>1 && nphoton>0 (delta r in nphoton?)
  //baseline include mZ>50 leadleppt>25(e)20(mu) ptgam/mllg>0.14 mll+mllg>185
  PlotMaker pm;

  bool plot_sig_accept = false;
  bool plot_other_objs = true;

  if (plot_sig_accept) {
    pm.Push<Hist1D>(Axis(40,0,4.0, mc_photon_drmin, "Truth Min #Delta R_{#gamma, l}", {}), 
        "1", procs_sig_full, ops).Weight(wgt*sigx50).Tag("zgtth");
    pm.Push<Table>("zgtth__reco", vector<TableRow>{
      TableRow("Inclusive", 
          "1",0,0,wgt),
      TableRow("1 Leptons Reconstructed", 
          leps_recod>=1,0,0,wgt),
      TableRow("2 Leptons Reconstructed", 
          leps_recod>=2,0,0,wgt),
      TableRow("Photon reconstructed (no dr)", 
          leps_recod>=2&&true_photon_recod,0,0,wgt),
      TableRow("Photon reconstructed (dr)", 
          leps_recod>=2&&"nphoton>0"&&true_photon_recod,0,0,wgt),
      TableRow("Photon reconstructed in slot 0 (no dr)", 
          leps_recod>=2&&nphoton_nodr>0.0&&"photon_pflavor[0]==1",0,0,wgt),
      TableRow("Photon reconstructed in slot 0 (dr)", 
          leps_recod>=2&&"nphoton>0&&photon_pflavor[0]==1",0,0,wgt),
      TableRow("Triggers", 
          leps_recod>=2&&"nphoton>0"&&true_photon_recod&&trigs,0,0,wgt),
      TableRow("b-Jet", 
          leps_recod>=2&&"nphoton>0"&&true_photon_recod&&trigs&&"nbdfl>=1",0,0,wgt),
      TableRow("Lepton", 
          leps_recod>=2&&"nlep>=3&&nphoton>0"&&true_photon_recod&&trigs,0,0,wgt),
      TableRow("Lepton+b-Jet", 
          leps_recod>=2&&"nlep>=3&&nphoton>0"&&true_photon_recod&&trigs&&"nbdfl>=1",0,0,wgt),
    },procs_sig_full,false,true,false,false,false,true).Precision(3);
  }

  if (plot_other_objs) {

    pm.Push<Hist1D>(Axis(30,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        "nlep>=2&&nphoton>=1&&nbdfl>=1", procs, ops).Weight(wgt*sigx50).Tag("zgtth");
    pm.Push<Hist1D>(Axis(30,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        "nlep>=3&&nphoton>=1&&nbdfl>=1", procs, ops).Weight(wgt*sigx50).Tag("zgtth");
    //to determine b-tagging
    pm.Push<Table>("zgtth__btag", vector<TableRow>{
      TableRow("Inclusive", 
          "nlep>=2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130",0,0,wgt),
      TableRow("1 Loose", 
          "nlep>=2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfl>=1",0,0,wgt),
      TableRow("1 Medium", 
          "nlep>=2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfm>=1",0,0,wgt),
      TableRow("1 Tight", 
          "nlep>=2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdft>=1",0,0,wgt),
      TableRow("2 Loose", 
          "nlep>=2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfl>=2",0,0,wgt),
      TableRow("Loose+Medium", 
          "nlep>=2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfl>=2&&nbdfm>=1",0,0,wgt),
      TableRow("Loose+Tight", 
          "nlep>=2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfl>=2&&nbdft>=1",0,0,wgt),
      TableRow("2 Medium", 
          "nlep>=2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfm>=2",0,0,wgt),
      TableRow("Medium+Tight", 
          "nlep>=2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfm>=2&&nbdft>=1",0,0,wgt),
      TableRow("2 Tight", 
          "nlep>=2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdft>=2",0,0,wgt),

      TableRow("2 L", 
          "nlep==2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130",0,0,wgt),
      TableRow("1 Loose", 
          "nlep==2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfl>=1",0,0,wgt),
      TableRow("1 Medium", 
          "nlep==2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfm>=1",0,0,wgt),
      TableRow("1 Tight", 
          "nlep==2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdft>=1",0,0,wgt),
      TableRow("2 Loose", 
          "nlep==2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfl>=2",0,0,wgt),
      TableRow("Loose+Medium", 
          "nlep==2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfl>=2&&nbdfm>=1",0,0,wgt),
      TableRow("Loose+Tight", 
          "nlep==2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfl>=2&&nbdft>=1",0,0,wgt),
      TableRow("2 Medium", 
          "nlep==2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfm>=2",0,0,wgt),
      TableRow("Medium+Tight", 
          "nlep==2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfm>=2&&nbdft>=1",0,0,wgt),
      TableRow("2 Tight", 
          "nlep==2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdft>=2",0,0,wgt),

      TableRow("3 L", 
          "nlep>=3&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130",0,0,wgt),
      TableRow("1 Loose", 
          "nlep>=3&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfl>=1",0,0,wgt),
      TableRow("1 Medium", 
          "nlep>=3&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfm>=1",0,0,wgt),
      TableRow("1 Tight", 
          "nlep>=3&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdft>=1",0,0,wgt),
      TableRow("2 Loose", 
          "nlep>=3&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfl>=2",0,0,wgt),
      TableRow("Loose+Medium", 
          "nlep>=3&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfl>=2&&nbdfm>=1",0,0,wgt),
      TableRow("Loose+Tight", 
          "nlep>=3&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfl>=2&&nbdft>=1",0,0,wgt),
      TableRow("2 Medium", 
          "nlep>=3&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfm>=2",0,0,wgt),
      TableRow("Medium+Tight", 
          "nlep>=3&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfm>=2&&nbdft>=1",0,0,wgt),
      TableRow("2 Tight", 
          "nlep>=3&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdft>=2",0,0,wgt),
    },procs,false,true,false,false,false,true).Precision(3);

    pm.Push<Hist1D>(Axis(30,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        "nlep==2&&nphoton>=1&&nbdfl>=2&&nbdfm>=1", procs, ops).Weight(wgt*sigx50).Tag("zgtth");
    pm.Push<Hist1D>(Axis(30,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        "nlep>=3&&nphoton>=1&&nbdfm>=1", procs, ops).Weight(wgt*sigx50).Tag("zgtth");
    pm.Push<Hist1D>(Axis(30,-1.0,1.0, "photon_idmva[0]", "Photon IDMVA", {}), 
        "nlep==2&&nphoton>=1&&nbdfl>=2&&nbdfm>=1", procs, ops).Weight(wgt*sigx50).Tag("zgtth");
    pm.Push<Hist1D>(Axis(30,-1.0,1.0, "photon_idmva[0]", "Photon IDMVA", {}), 
        "nlep>=3&&nphoton>=1&&nbdfm>=1", procs, ops).Weight(wgt*sigx50).Tag("zgtth");
    pm.Push<Hist1D>(Axis(50,0.0,150.0, "ll_m[0]", "m_{ll} [GeV]", {}), 
        "nlep==2&&nphoton>=1&&nbdfl>=2&&nbdfm>=1", procs, ops).Weight(wgt*sigx50).Tag("zgtth");
    pm.Push<Hist1D>(Axis(50,0.0,150.0, "ll_m[0]", "m_{ll} [GeV]", {}), 
        "nlep>=3&&nphoton>=1&&nbdfm>=1", procs, ops).Weight(wgt*sigx50).Tag("zgtth");

    pm.Push<Hist1D>(Axis(30,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        "nlep==2&&nphoton>=1&&nbdfl>=2&&nbdfm>=1&&photon_idmva[0]>0.5", procs, ops).Weight(wgt*sigx50).Tag("zgtth");
    pm.Push<Hist1D>(Axis(30,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        "nlep>=3&&nphoton>=1&&nbdfm>=1&&photon_idmva[0]>0.5", procs, ops).Weight(wgt*sigx50).Tag("zgtth");
    pm.Push<Hist1D>(Axis(30,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        "nlep==2&&nphoton>=1&&nbdfl>=2&&nbdfm>=1&&photon_idmva[0]>0.5&&ll_m[0]>80&&ll_m[0]<100", procs, ops).Weight(wgt*sigx50).Tag("zgtth");
    pm.Push<Hist1D>(Axis(30,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        "nlep>=3&&nphoton>=1&&nbdfm>=1&&photon_idmva[0]>0.5&&ll_m[0]>80&&ll_m[0]<100", procs, ops).Weight(wgt*sigx50).Tag("zgtth");

    pm.Push<Hist1D>(Axis(10,0.0,50.0, min_electron_pt, "min electron p_{T} [GeV]", {}), 
        "nlep>=3&&nel>=1&&nphoton>=1&&nbdfm>=1&&photon_idmva[0]>0.5&&ll_m[0]>80&&ll_m[0]<100", procs, ops).Weight(wgt*sigx50).Tag("zgtth");
    pm.Push<Hist1D>(Axis(10,0.0,50.0, min_muon_pt, "min muon p_{T} [GeV]", {}), 
        "nlep>=3&&nmu>=1&&nphoton>=1&&nbdfm>=1&&photon_idmva[0]>0.5&&ll_m[0]>80&&ll_m[0]<100", procs, ops).Weight(wgt*sigx50).Tag("zgtth");
    pm.Push<Hist1D>(Axis(10,0.0,50.0, assoc_electron_pt, "associated electron p_{T} [GeV]", {}), 
        "nlep>=3&&nphoton>=1&&nbdfm>=1&&photon_idmva[0]>0.5&&ll_m[0]>80&&ll_m[0]<100", procs, ops).Weight(wgt*sigx50).Tag("zgtth");
    pm.Push<Hist1D>(Axis(10,0.0,50.0, assoc_muon_pt, "associated muon p_{T} [GeV]", {}), 
        "nlep>=3&&nphoton>=1&&nbdfm>=1&&photon_idmva[0]>0.5&&ll_m[0]>80&&ll_m[0]<100", procs, ops).Weight(wgt*sigx50).Tag("zgtth");

    pm.Push<Hist1D>(Axis(30,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        ("nlep>=4"||("nlep>=3"&&min_electron_pt>15&&min_muon_pt>10))&&"nphoton>=1&&nbdfm>=1&&photon_idmva[0]>0.5&&ll_m[0]>80&&ll_m[0]<100", procs, ops).Weight(wgt*sigx50).Tag("zgtth");

    pm.Push<Hist1D>(Axis(30,0.0,200.0, "met", "p_{T}^{miss} [GeV]", {}), 
        "nlep==2&&nphoton>=1&&nbdfl>=2&&nbdfm>=1&&photon_idmva[0]>0.5&&ll_m[0]>80&&ll_m[0]<100", procs, ops).Weight(wgt*sigx50).Tag("zgtth");
    pm.Push<Hist1D>(Axis(30,0.0,200.0, "met", "p_{T}^{miss} [GeV]", {}), 
        "nlep>=3&&nphoton>=1&&nbdfm>=1&&photon_idmva[0]>0.5&&ll_m[0]>80&&ll_m[0]<100", procs, ops).Weight(wgt*sigx50).Tag("zgtth");
    pm.Push<Hist1D>(Axis(30,0.0,200.0, "mht", "H_{T}^{miss} [GeV]", {}), 
        "nlep==2&&nphoton>=1&&nbdfl>=2&&nbdfm>=1&&photon_idmva[0]>0.5&&ll_m[0]>80&&ll_m[0]<100", procs, ops).Weight(wgt*sigx50).Tag("zgtth");
    pm.Push<Hist1D>(Axis(30,0.0,200.0, "mht", "H_{T}^{miss} [GeV]", {}), 
        "nlep>=3&&nphoton>=1&&nbdfm>=1&&photon_idmva[0]>0.5&&ll_m[0]>80&&ll_m[0]<100", procs, ops).Weight(wgt*sigx50).Tag("zgtth");
    pm.Push<Hist1D>(Axis(30,0.0,400.0, "mht", "H_{T} [GeV]", {}), 
        "nlep==2&&nphoton>=1&&nbdfl>=2&&nbdfm>=1&&photon_idmva[0]>0.5&&ll_m[0]>80&&ll_m[0]<100", procs, ops).Weight(wgt*sigx50).Tag("zgtth");
    pm.Push<Hist1D>(Axis(30,0.0,400.0, "mht", "H_{T} [GeV]", {}), 
        "nlep>=3&&nphoton>=1&&nbdfm>=1&&photon_idmva[0]>0.5&&ll_m[0]>80&&ll_m[0]<100", procs, ops).Weight(wgt*sigx50).Tag("zgtth");
    pm.Push<Hist1D>(Axis(30,0.0,400.0, "njet", "N_{j}", {}), 
        "nlep==2&&nphoton>=1&&nbdfl>=2&&nbdfm>=1&&photon_idmva[0]>0.5&&ll_m[0]>80&&ll_m[0]<100", procs, ops).Weight(wgt*sigx50).Tag("zgtth");
    pm.Push<Hist1D>(Axis(30,0.0,400.0, "njet", "N_{j}", {}), 
        "nlep>=3&&nphoton>=1&&nbdfm>=1&&photon_idmva[0]>0.5&&ll_m[0]>80&&ll_m[0]<100", procs, ops).Weight(wgt*sigx50).Tag("zgtth");
  }

  pm.multithreaded_ = true;
  pm.min_print_ = true;
  pm.MakePlots(340.0);

  return 0;
}
