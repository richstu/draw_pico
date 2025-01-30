#include <vector>

#include "core/baby.hpp"
#include "core/named_func.hpp"
#include "core/named_func_utilities.hpp"
#include "zgamma/zg_functions.hpp"
#include "zgamma/zg_utilities.hpp"

using NamedFuncUtilities::ReduceNamedFunc;
using NamedFuncUtilities::MultiReduceNamedFunc;
using NamedFuncUtilities::FilterNamedFunc;
using NamedFuncUtilities::reduce_sum;
using NamedFuncUtilities::reduce_max;
using NamedFuncUtilities::reduce_maxfirst;
using NamedFuncUtilities::reduce_sublead;
using NamedFuncUtilities::reduce_subleadfirst;

namespace ZgFunctions {

  //isolated dielectron triggers for run 2
  const NamedFunc HLT_pass_dielectron("dielectron triggers",[](const Baby &b) -> NamedFunc::ScalarType{
    return b.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL()||b.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ();
  });

  //isolated dimuon triggers for run 2
  const NamedFunc HLT_pass_dimuon("dimuon triggers",[](const Baby &b) -> NamedFunc::ScalarType{
    if (abs(b.SampleType())==2016)
      return b.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ()||b.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ();
    return b.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8()||b.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8();
  });

  //isolated dilepton triggers for run 2
  const NamedFunc HLT_pass_dilepton = HLT_pass_dielectron||HLT_pass_dimuon;

  //isolated single electron triggers for run 2
  const NamedFunc HLT_pass_singleelectron("single electron triggers",[](const Baby &b) -> NamedFunc::ScalarType{
    if (abs(b.SampleType())==2016)
      return b.HLT_Ele27_WPTight_Gsf();
    if (abs(b.SampleType())==2017)
      return b.HLT_Ele35_WPTight_Gsf()||b.HLT_Ele32_WPTight_Gsf_L1DoubleEG();
    return b.HLT_Ele32_WPTight_Gsf();
  });

  //isolated single muon triggers for run 2
  const NamedFunc HLT_pass_singlemuon("single muon triggers",[](const Baby &b) -> NamedFunc::ScalarType{
    if (abs(b.SampleType())==2016)
      return b.HLT_IsoMu24();
    if (abs(b.SampleType())==2017)
      return b.HLT_IsoMu27()||b.HLT_IsoMu24();
    return b.HLT_IsoMu24();
  });

  //isolated single lepton triggers for run 2
  const NamedFunc HLT_pass_singlelepton = HLT_pass_singleelectron||HLT_pass_singlemuon;

  //year integrated lumi weights
  const NamedFunc w_years("w_years", [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleTypeString().Contains("-")) 
      return 1.; //data
    if (b.SampleTypeString()=="2016APV")
      return 19.51; //from brilcalc 
    if (b.SampleTypeString()=="2016")
      return 16.80; //from brilcalc
    else if (b.SampleTypeString()=="2017")
      return 41.48;
    else if (b.SampleTypeString()=="2018")
      return 59.83;
    else if (b.SampleTypeString()=="2022")
      return 8.17;
    else if (b.SampleTypeString()=="2022EE")
      return 27.01;
    else if (b.SampleTypeString()=="2023")
      return 17.68;
    //else if (b.SampleTypeString()=="2023BPix")
    return 9.53;
  });

  //year integrated lumi weights
  //const NamedFunc w_years_noapv("w_years", [](const Baby &b) -> NamedFunc::ScalarType{
  //  if (b.SampleType()<0) return 1.; //data
  //  if (b.SampleType()==2016)
  //    return 36.32264; 
  //  else if (b.SampleType()==2017)
  //    return 41.52756;
  //  return 59.67377;
  //});

  //Run 3 weight-up
  const NamedFunc w_run3("w_run3", [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleType()<0) return 1.; //data
    return 2.472;
  });

  //leading muon pt
  const NamedFunc lead_mu_pt = ReduceNamedFunc(FilterNamedFunc("mu_pt","mu_sig"),
      reduce_max).Name("lead_mu_pt");

  //subleading muon pt
  const NamedFunc sublead_mu_pt = ReduceNamedFunc(FilterNamedFunc("mu_pt","mu_sig"),
      reduce_sublead).Name("sublead_mu_pt");

  //leading electron pt
  const NamedFunc lead_el_pt = ReduceNamedFunc(FilterNamedFunc("el_pt","el_sig"),
      reduce_max).Name("lead_el_pt");

  //subleading electron pt
  const NamedFunc sublead_el_pt = ReduceNamedFunc(FilterNamedFunc("el_pt","el_sig"),
      reduce_sublead).Name("sublead_el_pt");

  //signal jet pt
  const NamedFunc sig_jet_pt = FilterNamedFunc("jet_pt","jet_isgood").Name("sig_jet_pt");

  //signal jet eta
  const NamedFunc sig_jet_eta = FilterNamedFunc("jet_eta","jet_isgood").Name("sig_jet_eta");

  //signal jet phi
  const NamedFunc sig_jet_phi = FilterNamedFunc("jet_phi","jet_isgood").Name("sig_jet_phi");

  //signal jet m
  const NamedFunc sig_jet_m = FilterNamedFunc("jet_m","jet_isgood").Name("sig_jet_m");

  //leading jet pt
  const NamedFunc lead_jet_pt = ReduceNamedFunc(sig_jet_pt,reduce_max).Name("lead_jet_pt");

  //subleading jet pt
  const NamedFunc sublead_jet_pt = ReduceNamedFunc(sig_jet_pt,reduce_sublead).Name("sublead_jet_pt");

  //lead jet eta
  const NamedFunc lead_jet_eta = MultiReduceNamedFunc(
      {sig_jet_pt,sig_jet_eta},reduce_maxfirst).Name("lead_jet_eta");

  //sublead jet eta
  const NamedFunc sublead_jet_eta = MultiReduceNamedFunc(
      {sig_jet_pt,sig_jet_eta},reduce_subleadfirst).Name("sublead_jet_eta");

  //lead jet phi
  const NamedFunc lead_jet_phi = MultiReduceNamedFunc(
      {sig_jet_pt,sig_jet_phi},reduce_maxfirst).Name("lead_jet_phi");

  //sublead jet phi
  const NamedFunc sublead_jet_phi = MultiReduceNamedFunc(
      {sig_jet_pt,sig_jet_phi},reduce_subleadfirst).Name("sublead_jet_phi");

  //lead jet m
  const NamedFunc lead_jet_m = MultiReduceNamedFunc(
      {sig_jet_pt,sig_jet_m},reduce_maxfirst).Name("lead_jet_m");

  //sublead jet m
  const NamedFunc sublead_jet_m = MultiReduceNamedFunc(
      {sig_jet_pt,sig_jet_m},reduce_subleadfirst).Name("sublead_jet_m");

  //trigger-paired lepton pt plateau cuts
  const NamedFunc trig_plateau_cuts("trig_plateau_cuts", [](const Baby &b) -> NamedFunc::ScalarType{
    double el_cut = 35;
    double mu_cut = 25;
    if (b.SampleTypeString()=="2016APV" || b.SampleTypeString()=="-2016APV" ||
        b.SampleTypeString()=="2016" || b.SampleTypeString()=="-2016")
      el_cut = 30;
    else if (b.SampleTypeString()=="2017" || b.SampleTypeString()=="-2017")
      mu_cut = 28;
    if (b.nel()>=1)
      if (b.trig_single_el() && b.el_pt()->at(0)>el_cut) return true;
    if (b.nel()>=2)
      if (b.trig_double_el() && b.el_pt()->at(0)>25 && b.el_pt()->at(1)>15) return true;
    if (b.nmu()>=1)
      if (b.trig_single_mu() && b.mu_pt()->at(0)>mu_cut) return true;
    if (b.nmu()>=2)
      if (b.trig_double_mu() && b.mu_pt()->at(0)>20 && b.mu_pt()->at(1)>10) return true;
    return false;
  });

  //common NamedFuncs for run 2 baseline selection
  NamedFunc zg_baseline_nolep = "nlep>=2 && nphoton>=1 && (ll_m[0]>50) && ((photon_pt[0]/llphoton_m[0])>=15.0/110.0) && ((llphoton_m[0]+ll_m[0])>=185)";
  NamedFunc zg_el_cuts = lead_el_pt>25&&sublead_el_pt>15;
  NamedFunc zg_mu_cuts = lead_mu_pt>20&&sublead_mu_pt>10;
  const NamedFunc zg_baseline_el_run2 = NamedFunc(zg_el_cuts && zg_baseline_nolep).Name("electron_baseline_run2");
  const NamedFunc zg_baseline_mu_run2 = NamedFunc(zg_mu_cuts && zg_baseline_nolep).Name("muon_baseline_run2");
  const NamedFunc zg_baseline_run2 = NamedFunc((zg_el_cuts || zg_mu_cuts) && zg_baseline_nolep).Name("baseline_run2");

  //new working baseline selection
  const NamedFunc zg_baseline = NamedFunc("(nel>=2||nmu>=2) && nphoton>=1 && ((photon_pt[0]/llphoton_m[0])>=15.0/110.0) && ll_m[0]>80 && ll_m[0]<100 && ((llphoton_m[0]+ll_m[0])>=185) && llphoton_m[0]>100 && llphoton_m[0]<180"&&trig_plateau_cuts).Name("baseline");

  //master stitch variable, updated for KingsCanyon_v0, later to be updated for v1
  const NamedFunc stitch("stitch",[](const Baby &b) -> NamedFunc::ScalarType{
    //remove ZGToLLG from DYJets
    if(b.type() >= 6000 && b.type() < 7000)
      return b.use_event();
    //remove DYJets from ZGToLLG
    if (b.type() >= 17000 && b.type() < 18000)
      return b.use_event();
    //remove ttG from TTJets - poor man for now
    if (b.type() >= 1000 && b.type() < 2000)
      //return b.stitch_photon(); currently bugged since photon doesn't exempt itself
      if (b.photon_pflavor()->size()>0)
        if (b.photon_pflavor()->at(0)==1)
          return 0;
    //various bugs in VVG samples, so don't worry about those
    //no WJets or EWKZ-ZG2J here either
    return 1.0;
  });

  //master stitch variable for deathvalley_v3
  const NamedFunc stitch_deathvalley("stitch_deathvalley",[](const Baby &b) -> NamedFunc::ScalarType{
    //remove ZGToLLG from DYJets
    if(b.type() >= 6000 && b.type() < 7000)
      return b.stitch_dy();
    //remove DYJets from ZGToLLG
    if (b.type() >= 17000 && b.type() < 18000)
      return !b.stitch_dy();
    //remove ttG from TTJets
    if (b.type() >= 1000 && b.type() < 2000)
      //return b.stitch_photon(); currently bugged since photon doesn't exempt itself
      if (b.photon_pflavor()->size()>0)
        if (b.photon_pflavor()->at(0)==1)
          return 0;
    //remove WWG from WW
    //if (b.type() >= 14000 && b.type() < 15000)
    //  return b.stitch_dy();
    //  //if (b.photon_pflavor()->size()>0)
    //  //  if (b.photon_pflavor()->at(0)==1)
    //  //    return 0;
    //remove WZG from WZ - poor man's stitch
    //if (b.type() >= 15000 && b.type() < 16000)
    //  return b.stitch_dy();
    //  //if (b.photon_pflavor()->size()>0)
    //  //  if (b.photon_pflavor()->at(0)==1)
    //  //    return 0;
    //don't use ZZG because only leptonic decays
    //remove ZZG from ZZ - poor man's stitch
    //if (b.type() >= 16000 && b.type() < 17000)
    //  if (b.photon_pflavor()->size()>0)
    //    if (b.photon_pflavor()->at(0)==1)
    //      return 0;
    return 1.0;
  });

  //drmax of lead photon
  const NamedFunc photon_drmax("photon_drmax",[](const Baby &b) -> NamedFunc::ScalarType{
    return ZgUtilities::pdrmax(b);
  });

  //relative pt uncertainty of lead photon for kingscanyon_v0 productions and earlier
  const NamedFunc photon_relpterr_deathvalley("photon_relpterr",[](const Baby &b) -> NamedFunc::ScalarType{
    return b.photon_energyErr()->at(0)/(b.photon_pt()->at(0)*TMath::CosH(b.photon_eta()->at(0)));
  });

  //relative pt uncertainty of lead photon
  const NamedFunc photon_relpterr("photon_relpterr",[](const Baby &b) -> NamedFunc::ScalarType{
    return b.photon_energyErr()->at(0)/(b.photon_pt()->at(0)*TMath::CosH(b.photon_eta()->at(0)));
  });

  //lead lepton eta (=lep_eta[0], but this isn't saved in slims =( )
  //only works for 2 lepton events
  const NamedFunc lead_lepton_eta("lead_lepton_eta",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.ll_lepid()->size() < 1) return 0;
    if (b.ll_lepid()->at(0) == 11) return b.el_eta()->at(b.ll_i1()->at(0));
    return b.mu_eta()->at(b.ll_i1()->at(0));
  });

  //sublead lepton eta (=lep_eta[0], but this isn't saved in slims =( )
  //only works for 2 lepton events
  const NamedFunc sublead_lepton_eta("sublead_lepton_eta",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.ll_lepid()->size() < 1) return 0;
    if (b.ll_lepid()->at(0) == 11) return b.el_eta()->at(b.ll_i2()->at(0));
    return b.mu_eta()->at(b.ll_i2()->at(0));
  });

  //pT/m of Higgs candidate
  const NamedFunc llphoton_rel_pt = NamedFunc("llphoton_pt[0]/llphoton_m[0]").Name("llphoton_rel_pt");

  //OR of triggers used in H->Zgamma analysis
  const NamedFunc trig = NamedFunc("trig_single_el||trig_single_mu||trig_double_el||trig_double_mu")
                         .Name("trig");

  //modified ptt variable
  const NamedFunc llphoton_pttmod("llphoton_pttmod",[](const Baby &b) -> NamedFunc::ScalarType{
    TVector3 zgamma, z, gamma;
    gamma.SetPtEtaPhi(b.photon_pt()->at(0),b.photon_eta()->at(0),b.photon_phi()->at(0)); 
    z.SetPtEtaPhi(b.ll_pt()->at(0),b.ll_eta()->at(0),b.ll_phi()->at(0)); 
    zgamma = z+gamma;
    return (zgamma.Cross(z-gamma)).Mag()/zgamma.Mag();
  });

  //maximum lepton mini isolation
  const NamedFunc max_lep_miniso("max_lep_miniso",[](const Baby &b) -> NamedFunc::ScalarType{
    double max_miniso = 0;
    for (unsigned iel = 0; iel < b.el_pt()->size(); iel++) {
      if (b.el_sig()->at(iel))
        if (b.el_miniso()->at(iel) > max_miniso)
          max_miniso = b.el_miniso()->at(iel);
    }
    for (unsigned imu = 0; imu < b.mu_pt()->size(); imu++) {
      if (b.mu_sig()->at(imu))
        if (b.mu_miniso()->at(imu) > max_miniso)
          max_miniso = b.mu_miniso()->at(imu);
    }
    return max_miniso;
  });

  //vector of whether GenParticles are the first copy
  const NamedFunc mc_isFirstCopy("mc_isFirstCopy",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> isFirstCopy;
    for (unsigned imc = 0; imc < b.mc_pt()->size(); imc++) {
      if ((b.mc_statusflag()->at(imc) & 0x1000) != 0)
        isFirstCopy.push_back(1.0);
      else
        isFirstCopy.push_back(0.0);
    }
    return isFirstCopy;
  });

  //vector mother PDGID of mother of GenParticle
  const NamedFunc mc_mommom("mc_mommom",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> mommom;
    for (unsigned imc = 0; imc < b.mc_pt()->size(); imc++) {
      int prev_idx = static_cast<int>(imc);
      int next_idx = b.mc_momidx()->at(imc);
      char level = 0;
      bool found_mommom = false;
      while (!found_mommom) {
        if (next_idx == -1) {
          mommom.push_back(-1);
          found_mommom = true;
        }
        else if (b.mc_id()->at(next_idx)==b.mc_id()->at(prev_idx)) {
          prev_idx = next_idx;
          next_idx = b.mc_momidx()->at(next_idx);
        }
        else if (level==0) {
          prev_idx = next_idx;
          next_idx = b.mc_momidx()->at(next_idx);
          level++;
        }
        else { //level = 1
          mommom.push_back(static_cast<double>(b.mc_id()->at(next_idx)));
          found_mommom = true;
        }
      }
    } //loop over mc particles
    return mommom;
  });

}







