#include "core/baby.hpp"
#include "core/named_func.hpp"
#include "zgamma/zg_functions.hpp"
#include "zgamma/zg_utilities.hpp"

//Includes for plots
#include "core/hist1d.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/process.hpp"

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
      return 17.61;
    //else if (b.SampleTypeString()=="2023BPix")
    return 9.53;
  });

 //This function will scale MC to 1 fb^-1 
 const NamedFunc wgt_1invfb("wgt_1invfb",[](const Baby &b) -> NamedFunc::ScalarType{ 
    if(b.SampleTypeString().Contains("-")) {
      return 1;
    }
    double w_year = w_years.GetScalar(b);
    double total_lum = 137.61;
    if(b.SampleTypeString().Contains("202")){total_lum = 62.32;}
    return b.weight()*w_year/total_lum;
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

  //common NamedFuncs for run 2 baseline selection - deprecated now?
  NamedFunc zg_baseline_nolep = "nlep>=2 && nphoton>=1 && (ll_m[0]>50) && ((photon_pt[0]/llphoton_m[0])>=15.0/110.0) && ((llphoton_m[0]+ll_m[0])>=185) && (photon_drmin[0]>0.4)";
  NamedFunc zg_el_cuts = "(ll_lepid[0]==11) && (el_pt[ll_i1[0]]>25) && (el_pt[ll_i2[0]]>15)";
  NamedFunc zg_mu_cuts = "(ll_lepid[0]==13) && (mu_pt[ll_i1[0]]>20) && (mu_pt[ll_i2[0]]>10)";
  const NamedFunc zg_baseline_el = NamedFunc(zg_el_cuts && zg_baseline_nolep).Name("electron_baseline");
  const NamedFunc zg_baseline_mu = NamedFunc(zg_mu_cuts && zg_baseline_nolep).Name("muon_baseline");
  const NamedFunc zg_baseline = NamedFunc((zg_el_cuts || zg_mu_cuts) && zg_baseline_nolep).Name("baseline");

  //master stitch variable
  const NamedFunc stitch("stitch",[](const Baby &b) -> NamedFunc::ScalarType{
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
    if (b.type() >= 14000 && b.type() < 15000)
      return b.stitch_dy();
      //if (b.photon_pflavor()->size()>0)
      //  if (b.photon_pflavor()->at(0)==1)
      //    return 0;
    //remove WZG from WZ - poor man's stitch
    if (b.type() >= 15000 && b.type() < 16000)
      return b.stitch_dy();
      //if (b.photon_pflavor()->size()>0)
      //  if (b.photon_pflavor()->at(0)==1)
      //    return 0;
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


  //Turns a vector<NamedFunc> into one usable for a cutflow table
  std::vector<NamedFunc> progressive_cuts(std::vector<NamedFunc> vector_NamedFunc){
    for(unsigned int idx = 1; idx<vector_NamedFunc.size(); idx++){ vector_NamedFunc[idx] = vector_NamedFunc[idx-1] && vector_NamedFunc[idx];}
    return vector_NamedFunc;
  }

  //This function adds all selections and reverses the selection at reverse
  NamedFunc Nreverse1(std::vector<NamedFunc> vector_NamedFunc, unsigned int reverse){
    NamedFunc return_NamedFunc = "1";
    for(unsigned int idx = 0; idx<vector_NamedFunc.size(); idx++){
      if(idx==reverse){return_NamedFunc = return_NamedFunc && !(vector_NamedFunc[idx]); continue;}  
      return_NamedFunc = return_NamedFunc && vector_NamedFunc[idx];
    }
    return return_NamedFunc;
  }

  //Returns a NamedFunc with all but one selection (marked by skip)
  NamedFunc Nminus1(std::vector<NamedFunc> vector_NamedFunc, unsigned int skip){
    NamedFunc return_NamedFunc = "1";
    for(unsigned int idx = 0; idx<vector_NamedFunc.size(); idx++){
      if(idx==skip){continue;}  
      return_NamedFunc = return_NamedFunc && vector_NamedFunc[idx];
    }
    return return_NamedFunc;
  }

  //Returns a NamedFunc without all selections in the vector skip
  NamedFunc Nminusk(std::vector<NamedFunc> vector_NamedFunc, std::vector<unsigned int> skip){
    unsigned int idx_s = 0;
    NamedFunc return_NamedFunc = "1";
    for(unsigned int idx = 0; idx<vector_NamedFunc.size(); idx++){
      if(idx_s < skip.size() && idx==skip[idx_s]){idx_s++; continue;}  
      return_NamedFunc = return_NamedFunc && vector_NamedFunc[idx];
    }
    return return_NamedFunc;
  }


  //pT/m of Higgs candidate
  const NamedFunc llphoton_rel_pt = NamedFunc("llphoton_pt[0]/llphoton_m[0]").Name("llphoton_rel_pt");

  //This version of llgamma mass is used when nphoton can be 0. This defaults to llphoton_m[0] otherwise
  const NamedFunc mlly("mlly",[](const Baby &b) -> NamedFunc::ScalarType{ return ZgUtilities::Findlly(b); });
  const NamedFunc pTy_mlly("pTy_mlly",[](const Baby &b) -> NamedFunc::ScalarType{ return (b.photon_pt()->at(0))/ZgUtilities::Findlly(b); });
  const NamedFunc mll_mlly("mll_mlly",[](const Baby &b) -> NamedFunc::ScalarType{ return (b.ll_m() ->at(0)) + ZgUtilities::Findlly(b); });
  const NamedFunc cosTheta("cosTheta",[](const Baby &b) -> NamedFunc::ScalarType{ return ZgUtilities::cos_Theta(b); });
  const NamedFunc costheta("costheta",[](const Baby &b) -> NamedFunc::ScalarType{ return ZgUtilities::cos_theta(b); });
  const NamedFunc Phi("Phi",[](const Baby &b) -> NamedFunc::ScalarType{ return ZgUtilities::Getphi(b); });

  //Below two functions are those that check if the triggers/trigger pT cuts are passed
  const NamedFunc pass_trigs = "trig_double_el||trig_double_mu||trig_single_el||trig_single_mu";
  const NamedFunc pass_singlelepton_trigs_and_pt = "(trig_single_mu && nmu > 0 && mu_pt[0] > 30) || (trig_single_el && nel > 0 && el_pt[0] > 37) ";
  const NamedFunc pass_dilepton_trigs_and_pt = "(trig_double_mu && nmu > 1 && mu_pt[0] > 20 && mu_pt[1] > 10) || (trig_double_el && nel > 1 && el_pt[0] > 25 && el_pt[1] > 15)";
  const NamedFunc pass_trigs_and_pt = pass_singlelepton_trigs_and_pt || pass_dilepton_trigs_and_pt;
  
  //This is a vector that gives the individual selections in the baseline selection
  const std::vector<NamedFunc> vector_run2_baseline      = {"nll>0", pass_trigs_and_pt, "ll_m[0] > 50", "nphoton > 0",  pTy_mlly > 15.0/110, mll_mlly > 185, mlly > 100 && mlly < 180};
  const std::vector<NamedFunc> vector_tightened_baseline = {"nll>0", pass_trigs_and_pt, "ll_m[0] > 50", "nphoton > 0 && photon_id80[0]",  pTy_mlly > 15.0/110, "ll_m[0] > 80 && ll_m[0] < 100",
                                                       mlly > 100 && mlly < 180};

  //Below is the baseline used by HIG-19-014
  const NamedFunc run2_baseline     = "nllphoton>0" && pass_trigs_and_pt && "ll_m[0] > 50 && photon_pt[0]/llphoton_m[0] > 15.0/110 && ll_m[0] + llphoton_m[0] > 185 && llphoton_m[0] > 100 && llphoton_m[0] < 180";
//  const NamedFunc run2_baseline      = "zg_cutBitMap[0]"; 

  //Below is the tightened baseline with the photon wp80 and the 80 GeV < mll < 100 GeV selection
  const NamedFunc tightened_baseline= "nllphoton>0" && pass_trigs_and_pt && "ll_m[0] > 50 && photon_id80[0] && photon_pt[0]/llphoton_m[0] > 15.0/110 && ll_m[0] > 80 && ll_m[0] < 100 && llphoton_m[0] > 100 && llphoton_m[0] < 180";

  //These are the control region plots made for the ggF category
  //Note mllg is not llphoton_m[0] because often control regions include photons that don't pass the photon object selections
  void ggF_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels, int ll_flavor){
    //-- Kinematic IDMVA Variables --//
    pm.Push<Hist1D>(Axis(40,-1,1,     cosTheta, "cos(#Theta)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_ggF_4_cTHETA_" + labels);
    pm.Push<Hist1D>(Axis(40,-1,1,     costheta, "cos(#theta)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_ggF_6_ctheta_" + labels);
    pm.Push<Hist1D>(Axis(40,-3.2,3.2, Phi,      "#phi",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_ggF_Z11_phi_" + labels);
    pm.Push<Hist1D>(Axis(80,60,120,   "ll_m[0]","m_{ll} [GeV]",  {80,100}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_ggF_Z12_mll" + labels);
    pm.Push<Hist1D>(Axis(40,100,180,  mlly     ,"m_{ll#gamma} [GeV]",  {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_ggF_Z13_mlly" + labels);

    //Creating the correct eta plots
    if(ll_flavor!=13){  
      pm.Push<Hist1D>(Axis(40,5,65,     "el_pt[ll_i1[0]]",  "p_{T}(e_{1})", {}), selection && "ll_lepid[0]==11", processes, ops).Weight(wgt).Tag("ShortName:CR_ggF_8_pt_el1_" + labels);
      pm.Push<Hist1D>(Axis(40,5,65,     "el_pt[ll_i2[0]]",  "p_{T}(e_{2})", {}), selection && "ll_lepid[0]==11", processes, ops).Weight(wgt).Tag("ShortName:CR_ggF_8_pt_el2_" + labels);
      pm.Push<Hist1D>(Axis(40,-2.6,2.6, "el_eta[ll_i1[0]]", "#eta(e_{1})",  {}), selection && "ll_lepid[0]==11", processes, ops).Weight(wgt).Tag("ShortName:CR_ggF_9_eta_el1_" + labels);
      pm.Push<Hist1D>(Axis(40,-2.6,2.6, "el_eta[ll_i2[0]]", "#eta(e_{2})",  {}), selection && "ll_lepid[0]==11", processes, ops).Weight(wgt).Tag("ShortName:CR_ggF_9_eta_el2_" + labels);

    }
    if(ll_flavor!=11){
      pm.Push<Hist1D>(Axis(40,0,75,     "mu_pt[ll_i1[0]]",  "p_{T}(#mu_{1})", {}), selection && "ll_lepid[0]==13", processes, ops).Weight(wgt).Tag("ShortName:CR_ggF_8_pt_mu1_" + labels);
      pm.Push<Hist1D>(Axis(40,0,75,     "mu_pt[ll_i2[0]]",  "p_{T}(#mu_{2})", {}), selection && "ll_lepid[0]==13", processes, ops).Weight(wgt).Tag("ShortName:CR_ggF_8_pt_mu2_" + labels);
      pm.Push<Hist1D>(Axis(40,-2.6,2.6, "mu_eta[ll_i1[0]]", "#eta(#mu_{1})",  {}), selection && "ll_lepid[0]==13", processes, ops).Weight(wgt).Tag("ShortName:CR_ggF_9_eta_mu1_" + labels);
      pm.Push<Hist1D>(Axis(40,-2.6,2.6, "mu_eta[ll_i2[0]]", "#eta(#mu_{2})",  {}), selection && "ll_lepid[0]==13", processes, ops).Weight(wgt).Tag("ShortName:CR_ggF_9_eta_mu2_" + labels);
    }
  
    pm.Push<Hist1D>(Axis(45,0.0,1,    "photon_idmva[0]","#gamma - IDMVA",                         {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_ggF_1_IDMVA_" + labels);
    pm.Push<Hist1D>(Axis(35,0,3.5,    "photon_drmin[0]","min( #Delta R(#gamma,l) )",              {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_ggF_2_minDR_" + labels);
    pm.Push<Hist1D>(Axis(20,0.01,0.25,"photon_pterr[0]/photon_pt[0]", "#sigma(#gamma)/E(#gamma)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_ggF_5_relunc_gamma_" + labels);
    pm.Push<Hist1D>(Axis(40,10,60,    "photon_pt[0]",   "p_{T}(#gamma)",                          {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_ggF_5b_pT_gamma_" + labels);
    pm.Push<Hist1D>(Axis(40,-2.6,2.6, "photon_eta[0]",  "#eta(#gamma)",                           {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_ggF_Z10_eta_gamma_" + labels);

    pm.Push<Hist1D>(Axis(40,0.05,0.25, pTy_mlly,     "p_{T}(ll#gamma)/m_{ll#gamma}", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_ggF_3_ptlly_mlly_" + labels);
    pm.Push<Hist1D>(Axis(45,0.5,5,     photon_drmax, "max( #Delta R(#gamma,l) )",    {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_ggF_7_maxDR_" + labels);

  }


  //This function again adds control region plots for the VBF category based on variables used in Run 2. Some overlap with ggF region but removing variables specific to the "kinematic BDT"
  void VBF_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> processes, std::vector<PlotOpt> ops, NamedFunc wgt, std::string labels, int ll_flavor){
    pm.Push<Hist1D>(Axis(12,400,800, "dijet_m",    "m_{jj} [GeV]",            {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_1_dijet_m" + labels);
    pm.Push<Hist1D>(Axis(12,0,9,     "dijet_deta", "#Delta#eta(j_{1},j_{2})", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_2_dijet_deta" + labels);
    pm.Push<Hist1D>(Axis(10,0,3,     "dijet_dphi", "#Delta#phi(j_{1},j_{2})", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_3_dijet_dphi" + labels);

    pm.Push<Hist1D>(Axis(30, 0, 3,   "llphoton_dijet_dphi[0]" ,"#Delta#phi(Z#gamma,jj)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_dijet_4_jjlly_dphi" + labels);
    pm.Push<Hist1D>(Axis(30, 0, 1.0, "llphoton_dijet_balance[0]" ,"System balance",         {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_dijet_5_sysbal" + labels);
    pm.Push<Hist1D>(Axis(20, 0, 140, "llphoton_pTt2[0]"       ,"p_{Tt} [GeV]",           {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_dijet_6_pTt2" + labels);

    pm.Push<Hist1D>(Axis( 8,   0,  6, "photon_zeppenfeld[0]", "#gamma zeppenfeld",      {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_dijet_7_zep" + labels);
    pm.Push<Hist1D>(Axis(15, 0.4,4.4, "photon_jet_mindr[0]",  "Min. #DeltaR(#gamma,j)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_dijet_8_yjdr" + labels);

    pm.Push<Hist1D>(Axis(90,30,300,      "jet_pt[0]",  "p_{T}(j_{1}) [GeV]",{}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_9_j1_pt" + labels);
    pm.Push<Hist1D>(Axis(90,30,300,      "jet_pt[1]",  "p_{T}(j_{2}) [GeV]",{}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_9_j2_pt" + labels);
    pm.Push<Hist1D>(Axis(80,-4.7,4.7,    "jet_eta[0]", "#eta(j_{1})", {}),       selection, processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_Z10_j1_eta_" + labels);
    pm.Push<Hist1D>(Axis(80,-4.7,4.7,    "jet_eta[1]", "#eta(j_{2})", {}),       selection, processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_Z10_j2_eta_" + labels);
    pm.Push<Hist1D>(Axis(63, -3.15, 3.15, "jet_phi[0]", "#phi(j_{1})", {}),       selection, processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_Z11_j1_phi_" + labels);
    pm.Push<Hist1D>(Axis(63, -3.15, 3.15, "jet_phi[1]", "#phi(j_{2})", {}),       selection, processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_Z11_j2_phi_" + labels);

    if(ll_flavor!=13){  
      pm.Push<Hist1D>(Axis(40,5,65,     "el_pt[ll_i1[0]]",  "p_{T}(e_{1})", {}), selection && "ll_lepid[0]==11", processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_12_pt_el1_" + labels);
      pm.Push<Hist1D>(Axis(40,5,65,     "el_pt[ll_i2[0]]",  "p_{T}(e_{2})", {}), selection && "ll_lepid[0]==11", processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_12_pt_el2_" + labels);
      pm.Push<Hist1D>(Axis(40,-2.6,2.6, "el_eta[ll_i1[0]]", "#eta(e_{1})",  {}), selection && "ll_lepid[0]==11", processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_13_eta_el1_" + labels);
      pm.Push<Hist1D>(Axis(40,-2.6,2.6, "el_eta[ll_i2[0]]", "#eta(e_{2})",  {}), selection && "ll_lepid[0]==11", processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_13_eta_el2_" + labels);

    }
    if(ll_flavor!=11){
      pm.Push<Hist1D>(Axis(40,0,75,     "mu_pt[ll_i1[0]]",  "p_{T}(#mu_{1})", {}), selection && "ll_lepid[0]==13", processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_12_pt_mu1_" + labels);
      pm.Push<Hist1D>(Axis(40,0,75,     "mu_pt[ll_i2[0]]",  "p_{T}(#mu_{2})", {}), selection && "ll_lepid[0]==13", processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_12_pt_mu2_" + labels);
      pm.Push<Hist1D>(Axis(40,-2.6,2.6, "mu_eta[ll_i1[0]]", "#eta(#mu_{1})",  {}), selection && "ll_lepid[0]==13", processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_13_eta_mu1_" + labels);
      pm.Push<Hist1D>(Axis(40,-2.6,2.6, "mu_eta[ll_i2[0]]", "#eta(#mu_{2})",  {}), selection && "ll_lepid[0]==13", processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_13_eta_mu2_" + labels);
    }

    pm.Push<Hist1D>(Axis(40,0.05,0.25, pTy_mlly,     "p_{T}(ll#gamma)/m_{ll#gamma}", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_Z14_ptlly_mlly_" + labels);
    pm.Push<Hist1D>(Axis(45,0.5,5,     photon_drmax, "max( #Delta R(#gamma,l) )",    {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_Z15_maxDR_" + labels);
    pm.Push<Hist1D>(Axis(80,60,120,    "ll_m[0]",    "m_{ll} [GeV]",           {80,100}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_Z16_mll" + labels);
    pm.Push<Hist1D>(Axis(40,100,180,   mlly,         "m_{ll#gamma} [GeV]",           {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_Z17_mlly" + labels);

    pm.Push<Hist1D>(Axis(45,0.0,1,    "photon_idmva[0]","#gamma - IDMVA",                         {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_Zph1_IDMVA_" + labels);
    pm.Push<Hist1D>(Axis(35,0,3.5,    "photon_drmin[0]","min( #Delta R(#gamma,l) )",              {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_Zph2_minDR_" + labels);
    pm.Push<Hist1D>(Axis(20,0.01,0.25,"photon_pterr[0]/photon_pt[0]", "#sigma(#gamma)/E(#gamma)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_Zph3_relunc_gamma_" + labels);
    pm.Push<Hist1D>(Axis(40,10,60,    "photon_pt[0]",   "p_{T}(#gamma)",                          {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_Zph4_pT_gamma_" + labels);
    pm.Push<Hist1D>(Axis(40,-2.6,2.6, "photon_eta[0]",  "#eta(#gamma)",                           {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_VBF_Zph5_eta_gamma_" + labels);
  }

  //This function, as of 11/03/23 is a work in progress. No control region plots were made in HIG-19-014
  void Lep_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> processes, std::vector<PlotOpt> ops, NamedFunc wgt, std::string labels, int ll_flavor){
    pm.Push<Hist1D>(Axis(80,60,120,    "ll_m[0]",    "m_{ll} [GeV]",           {80,100}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_LEP_1_mll" + labels);
    pm.Push<Hist1D>(Axis(40,100,180,   mlly,         "m_{ll#gamma} [GeV]",           {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_LEP_2_mlly" + labels);
    pm.Push<Hist1D>(Axis(40,0.05,0.25, pTy_mlly,     "p_{T}(ll#gamma)/m_{ll#gamma}", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_LEP_3_ptlly_mlly_" + labels);
    pm.Push<Hist1D>(Axis(45,0.5,5,     photon_drmax, "max( #Delta R(#gamma,l) )",    {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_LEP_4_maxDR_" + labels);

    pm.Push<Hist1D>(Axis(45,0.0,1,    "photon_idmva[0]","#gamma - IDMVA",                         {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_LEP_5_IDMVA_" + labels);
    pm.Push<Hist1D>(Axis(35,0,3.5,    "photon_drmin[0]","min( #Delta R(#gamma,l) )",              {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_LEP_6_minDR_" + labels);
    pm.Push<Hist1D>(Axis(20,0.01,0.25,"photon_pterr[0]/photon_pt[0]", "#sigma(#gamma)/E(#gamma)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_LEP_7_relunc_gamma_" + labels);
    pm.Push<Hist1D>(Axis(40,10,60,    "photon_pt[0]",   "p_{T}(#gamma)",                          {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_LEP_8_pT_gamma_" + labels);
    pm.Push<Hist1D>(Axis(40,-2.6,2.6, "photon_eta[0]",  "#eta(#gamma)",                           {}), selection, processes, ops).Weight(wgt).Tag("ShortName:CR_LEP_9_eta_gamma_" + labels);

    if(ll_flavor!=13){  
      pm.Push<Hist1D>(Axis(40,5,65,     "el_pt[ll_i1[0]]",  "p_{T}(e_{1})", {}), selection && "ll_lepid[0]==11", processes, ops).Weight(wgt).Tag("ShortName:CR_LEP_10_pt_el1_" + labels);
      pm.Push<Hist1D>(Axis(40,5,65,     "el_pt[ll_i2[0]]",  "p_{T}(e_{2})", {}), selection && "ll_lepid[0]==11", processes, ops).Weight(wgt).Tag("ShortName:CR_LEP_10_pt_el2_" + labels);
      pm.Push<Hist1D>(Axis(40,-2.6,2.6, "el_eta[ll_i1[0]]", "#eta(e_{1})",  {}), selection && "ll_lepid[0]==11", processes, ops).Weight(wgt).Tag("ShortName:CR_LEP_11_eta_el1_" + labels);
      pm.Push<Hist1D>(Axis(40,-2.6,2.6, "el_eta[ll_i2[0]]", "#eta(e_{2})",  {}), selection && "ll_lepid[0]==11", processes, ops).Weight(wgt).Tag("ShortName:CR_LEP_11_eta_el2_" + labels);

    }
    if(ll_flavor!=11){
      pm.Push<Hist1D>(Axis(40,0,75,     "mu_pt[ll_i1[0]]",  "p_{T}(#mu_{1})", {}), selection && "ll_lepid[0]==13", processes, ops).Weight(wgt).Tag("ShortName:CR_LEP_10_pt_mu1_" + labels);
      pm.Push<Hist1D>(Axis(40,0,75,     "mu_pt[ll_i2[0]]",  "p_{T}(#mu_{2})", {}), selection && "ll_lepid[0]==13", processes, ops).Weight(wgt).Tag("ShortName:CR_LEP_10_pt_mu2_" + labels);
      pm.Push<Hist1D>(Axis(40,-2.6,2.6, "mu_eta[ll_i1[0]]", "#eta(#mu_{1})",  {}), selection && "ll_lepid[0]==13", processes, ops).Weight(wgt).Tag("ShortName:CR_LEP_11_eta_mu1_" + labels);
      pm.Push<Hist1D>(Axis(40,-2.6,2.6, "mu_eta[ll_i2[0]]", "#eta(#mu_{2})",  {}), selection && "ll_lepid[0]==13", processes, ops).Weight(wgt).Tag("ShortName:CR_LEP_11_eta_mu2_" + labels);
    }


  }

}







