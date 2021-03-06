// Lists the MC cross sections

#include <iostream>
#include "core/cross_sections.hpp"

using namespace std;

namespace xsec{

  float crossSection(const TString &file, bool is2016){
    float xsec(999999.), Htobb(0.5824);

    if (is2016) {
        if(file.Contains("SMS-T1tttt_mGluino-1200_mLSP-800_Tune")) xsec = 0.0985;
        if(file.Contains("SMS-T1tttt_mGluino-2000_mLSP-100_Tune")) xsec = 0.00101;
        
        if(file.Contains("RPV") && file.Contains("1000"))  xsec = 0.325388;
        if(file.Contains("RPV") && file.Contains("1100"))  xsec = 0.163491;
        if(file.Contains("RPV") && file.Contains("1200"))  xsec = 0.0856418;
        if(file.Contains("RPV") && file.Contains("1300"))  xsec = 0.0460525;
        if(file.Contains("RPV") && file.Contains("1400"))  xsec = 0.0252977;


        //  Cross-section taken from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO
        // Alternative option: https://twiki.cern.ch/twiki/bin/view/Sandbox/FullNNLOcrossSections#Top_cross_section_for_13_TeV
        if(file.Contains("TTJets_Tune") || file.Contains("TT_"))  xsec = 815.96;
        if(file.Contains("TTJets_HT")){//LO cross sections with k-factor of 1.625 already applied
          if(file.Contains("2500toInf")) xsec = 0.0023234211;
          if(file.Contains("1200to2500")) xsec = 0.194972521;
          if(file.Contains("800to1200")) xsec = 1.07722318;
          if(file.Contains("600to800")) xsec = 2.61537118;

        }
        // The efficiency of the mtt>1000 cut is taken from sigma(mtt>1000)/sigma(inclusive) from mcm
        double mtt_1000_eff=(11.16/670.3);
        if(file.Contains("TTJets_Mtt-1000toInf")) xsec = 815.96*mtt_1000_eff;

        if(file.Contains("TTJets_DiLept") || file.Contains("TTTo2L2Nu")) xsec = 85.66; // (3*0.108)^2*815.96
        if(file.Contains("TTJets_SingleLept")) xsec = 178.7; //(1- ((1-3*0.108)^2+(3*0.108)^2))*815.96*0.5 per half
        if(file.Contains("TTToSemiLeptonic")) xsec = 357.4;
        
        if(file.Contains("TTJets_DiLept_genMET-150")) xsec = 0.06392*85.66; // filter_eff*(3*0.108)^2*815.96
        if(file.Contains("TTJets_SingleLept") && file.Contains("genMET-150")) xsec = 0.05217*178.7; //filter_eff*(1- ((1-3*0.108)^2+(3*0.108)^2))*815.96*0.5 per half

        // cross sections from mcm
        if(file.Contains("TTG")) xsec = 3.697;                
        if(file.Contains("TTTT_Tune")) xsec = 0.009103;
        // mcm cross section times the same kfactors used for leptonic samples
        if(file.Contains("WJetsToQQ_HT-600ToInf")) xsec = 95.14*1.21;
        if(file.Contains("ZJetsToQQ_HT600toInf")) xsec = 5.67*1.23;
        
        // From https://cms-pdmv.cern.ch/mcm
        // k-factors from https://mangano.web.cern.ch/mangano/public/MECCA/samples_50ns_miniaod.txt
        // k-factors are ratio of https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
        // NLO/NNLO cross-sections to that of an inclusive sample in mcm at lower order (LO/NLO)

        if(file.Contains("WJetsToLNu_Tune")) xsec=61526.7; //NNLO from Lesya's summary table 

        //cross-section per slice changed due to change in genHT definition
        if(file.Contains("WJetsToLNu_HT-70To100"))  xsec = 1372.*1.21; 
        if(file.Contains("WJetsToLNu_HT-100To200"))  xsec = 1347.*1.21; 
        if(file.Contains("WJetsToLNu_HT-200To400"))  xsec = 360.*1.21;
        if(file.Contains("WJetsToLNu_HT-400To600"))  xsec = 48.98*1.21;
        if(file.Contains("WJetsToLNu_HT-600ToInf"))  xsec = 18.77*1.21;
        if(file.Contains("WJetsToLNu_HT-600To800"))  xsec = 12.05*1.21;
        if(file.Contains("WJetsToLNu_HT-800To1200"))  xsec = 5.501*1.21;
        if(file.Contains("WJetsToLNu_HT-1200To2500"))  xsec = 1.329*1.21;
        if(file.Contains("WJetsToLNu_HT-2500ToInf"))  xsec = 0.03216*1.21;

        //updated 02-2019 with XSDB
        if(file.Contains("QCD_HT100to200_Tune")) xsec = 28060000;
        if(file.Contains("QCD_HT200to300_Tune"))   xsec = 1710000;
        if(file.Contains("QCD_HT300to500_Tune"))   xsec = 347500;
        if(file.Contains("QCD_HT500to700_Tune"))   xsec = 32060;
        if(file.Contains("QCD_HT700to1000_Tune"))  xsec = 6829;
        if(file.Contains("QCD_HT1000to1500_Tune")) xsec = 1207;
        if(file.Contains("QCD_HT1500to2000_Tune")) xsec = 120.0;
        if(file.Contains("QCD_HT2000toInf_Tune"))  xsec = 25.25;

        // Cross sections from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopRefXsec
        // multiplied by BF(W->mu,e,tau) = 0.324
        if (file.Contains("ST_s-channel_4f_leptonDecays"))     xsec = 3.34;
        if (file.Contains("ST_t-channel_antitop_4f_inclusiveDecays")) xsec = 80.95;
        if (file.Contains("ST_t-channel_top_4f_inclusiveDecays")) xsec = 136.02;
        if (file.Contains("ST_tW_antitop_5f_NoFullyHadronicDecays"))     xsec = 35.85*0.543;
        if (file.Contains("ST_tW_top_5f_NoFullyHadronicDecays"))     xsec = 35.85*0.543; 

        if(file.Contains("DYJetsToLL_M-10to50_Tune")) xsec = 18610*1.23;
        if(file.Contains("DYJetsToLL_M-50_Tune"))     xsec = 4895*1.23;

        if(file.Contains("DYJetsToLL_M-50_HT-70to100"))    xsec = 175.3*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-100to200"))    xsec = 139.4*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-200to400"))    xsec = 42.75*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-400to600"))    xsec = 5.497*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-600to800"))    xsec = 1.363*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-800to1200"))    xsec = 0.6759*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-1200to2500"))    xsec = 0.116*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-2500toInf"))    xsec = 0.002592*1.23;
        if(file.Contains("DYJetsToLL_M-50_HT-600toInf"))    xsec = 2.21*1.23;

        if(file.Contains("ZJetsToNuNu_HT-100To200"))  xsec = 280.35*1.27;
        if(file.Contains("ZJetsToNuNu_HT-200To400"))  xsec = 77.67*1.27;
        if(file.Contains("ZJetsToNuNu_HT-400To600"))  xsec = 10.73*1.27;
        if(file.Contains("ZJetsToNuNu_HT-600To800"))  xsec = 2.536*1.27;
        if(file.Contains("ZJetsToNuNu_HT-800To1200"))  xsec = 1.161*1.27;
        if(file.Contains("ZJetsToNuNu_HT-1200To2500"))  xsec = 0.2824*1.27;
        if(file.Contains("ZJetsToNuNu_HT-2500ToInf"))  xsec = 0.006459*1.27;
        if(file.Contains("ZJetsToNuNu_HT-600ToInf"))  xsec = 3.986*1.27;

        if(file.Contains("TTZToQQ"))                xsec = 0.5297;
        if(file.Contains("TTZToLLNuNu_M-10"))       xsec = 0.2529;
        if(file.Contains("TTWJetsToQQ"))            xsec = 0.4062;
        if(file.Contains("TTWJetsToLNu"))           xsec = 0.2043;
       
        //https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Diboson
        if(file.Contains("WWTo2L2Nu"))   xsec = 12.178; //NNLO
        if(file.Contains("WWToLNuQQ"))   xsec = 49.997; //NNLO
        if(file.Contains("ttHTobb_M125"))   xsec = 0.2934;

        if(file.Contains("WZTo1L3Nu"))   xsec = 3.05;
        if(file.Contains("WZTo1L1Nu2Q"))   xsec = 10.96;
        if(file.Contains("WZTo2L2Q"))   xsec = 5.595;
        if(file.Contains("WZTo3LNu"))   xsec = 4.42965;
        if(file.Contains("VVTo2L2Nu"))   xsec = 11.95;
        if(file.Contains("ZZ_Tune"))   xsec = 16.523;

        // Calculated at 13 TeV in
        // https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt1314TeV
        // Higgs branching ratios from
        // https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR
        if(file.Contains("ZH_HToBB_ZToLL_M-125"))      xsec = 0.883*Htobb*0.033658;
        if(file.Contains("ZH_HToBB_ZToNuNu_M-125"))    xsec = 0.883*Htobb*0.2;
        if(file.Contains("WH_HToBB_WToLNu_M-125"))     xsec = 1.373*Htobb*(0.1071+0.1063+0.1138);
        if(file.Contains("ZH_HToBB_ZToNuNu_M125"))    xsec = 0.883*Htobb*0.2;
        if(file.Contains("WH_HToBB_WToLNu_M125"))     xsec = 1.373*Htobb*(0.1071+0.1063+0.1138);
    } else {
        if(file.Contains("SMS-T1tttt_mGluino-1200_mLSP-800_Tune")) xsec = 0.0985;
        if(file.Contains("SMS-T1tttt_mGluino-2000_mLSP-100_Tune")) xsec = 0.00101;

        //  Cross-section taken from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO
        // Alternative option: https://twiki.cern.ch/twiki/bin/view/Sandbox/FullNNLOcrossSections#Top_cross_section_for_13_TeV
        if(file.Contains("TTJets_Tune") || file.Contains("TT_"))  xsec = 815.96;

        if(file.Contains("TTJets_DiLept") || file.Contains("TTTo2L2Nu")) xsec = 85.66; // (3*0.108)^2*815.96
        if(file.Contains("TTJets_SingleLept")) xsec = 178.7; //(1- ((1-3*0.108)^2+(3*0.108)^2))*815.96*0.5 per half
        
        if(file.Contains("TTJets_DiLept_genMET-150")) xsec = 0.0676543*85.66; // filter_eff*(3*0.108)^2*815.96
        if(file.Contains("TTJets_SingleLept") && file.Contains("genMET-150")) xsec = 0.0568246*178.7; //filter_eff*(1- ((1-3*0.108)^2+(3*0.108)^2))*815.96*0.5 per half

        if(file.Contains("TTJets_DiLept_genMET-80")) xsec = 0.412882*85.66; // filter_eff*(3*0.108)^2*815.96
        if(file.Contains("TTJets_SingleLept") && file.Contains("genMET-80")) xsec = 0.293137*178.7; //filter_eff*(1- ((1-3*0.108)^2+(3*0.108)^2))*815.96*0.5 per half

        // from cross XSDB
        if(file.Contains("TTG")) xsec = 4.078;                
        if(file.Contains("TTTT_Tune")) xsec = 0.008213;
        
        // From https://cms-pdmv.cern.ch/mcm
        // k-factors from https://mangano.web.cern.ch/mangano/public/MECCA/samples_50ns_miniaod.txt
        // k-factors are ratio of https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
        // NLO/NNLO cross-sections to that of an inclusive sample in mcm at lower order (LO/NLO)

        if(file.Contains("WJetsToLNu_Tune")) xsec=20508.9*3; //NNLO from Lesya's summary table 

        //cross-section per slice based on inclusive sample, roughly 10% higher than 2016, less in extreme tail
        if(file.Contains("WJetsToLNu_HT-100To200"))  xsec = 0.0262096*20508.9*3;
        if(file.Contains("WJetsToLNu_HT-200To400"))  xsec = 0.00772818*20508.9*3;
        if(file.Contains("WJetsToLNu_HT-400To600"))  xsec = 0.00109366*20508.9*3;
        if(file.Contains("WJetsToLNu_HT-600To800"))  xsec = 0.000272388*20508.9*3;
        if(file.Contains("WJetsToLNu_HT-800To1200"))  xsec = 0.000122233*20508.9*3;
        if(file.Contains("WJetsToLNu_HT-1200To2500"))  xsec = 2.71060e-5*20508.9*3;
        if(file.Contains("WJetsToLNu_HT-2500ToInf"))  xsec = 3.94174e-07*20508.9*3;

        if(file.Contains("QCD_HT100to200_Tune")) xsec = 23700000;
        if(file.Contains("QCD_HT200to300_Tune"))   xsec = 1547000;
        if(file.Contains("QCD_HT300to500_Tune"))   xsec = 322600;
        if(file.Contains("QCD_HT500to700_Tune"))   xsec = 29980;
        if(file.Contains("QCD_HT700to1000_Tune"))  xsec = 6334;
        if(file.Contains("QCD_HT1000to1500_Tune")) xsec = 1088;
        if(file.Contains("QCD_HT1500to2000_Tune")) xsec = 99.11;
        if(file.Contains("QCD_HT2000toInf_Tune"))  xsec = 20.23;

        // Cross sections from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopRefXsec
        // multiplied by BF(W->mu,e,tau) = 0.324
        if (file.Contains("ST_s-channel_4f_leptonDecays"))     xsec = 3.34;
        if (file.Contains("ST_t-channel_antitop_4f_inclusiveDecays") ||
            file.Contains("ST_t-channel_antitop_4f_InclusiveDecays")) xsec = 80.95;
        if (file.Contains("ST_t-channel_top_4f_inclusiveDecays") || 
            file.Contains("ST_t-channel_top_4f_InclusiveDecays")) xsec = 136.02;
        if (file.Contains("ST_tW_antitop_5f_NoFullyHadronicDecays"))     xsec = 35.85*0.543;
        if (file.Contains("ST_tW_top_5f_NoFullyHadronicDecays"))     xsec = 35.85*0.543; 

        if(file.Contains("DYJetsToLL_M-50_Tune"))     xsec = 2075.14*3;

        if(file.Contains("DYJetsToLL_M-50_HT-100to200"))    xsec = 0.0302083*2075.14*3;
        if(file.Contains("DYJetsToLL_M-50_HT-200to400"))    xsec = 0.00907651*2075.14*3;
        if(file.Contains("DYJetsToLL_M-50_HT-400to600"))    xsec = 0.00129238*2075.14*3;
        if(file.Contains("DYJetsToLL_M-50_HT-600to800"))    xsec = 0.000316039*2075.14*3;
        if(file.Contains("DYJetsToLL_M-50_HT-800to1200"))    xsec = 0.000137432*2075.14*3;
        if(file.Contains("DYJetsToLL_M-50_HT-1200to2500"))    xsec = 3.09368e-05*2075.14*3;
        if(file.Contains("DYJetsToLL_M-50_HT-2500toInf"))    xsec = 4.39860e-07*2075.14*3;

        // k-factor from DYJets 1.165
        if(file.Contains("ZJetsToNuNu_HT-100To200"))  xsec = 302.8*1.165;
        if(file.Contains("ZJetsToNuNu_HT-200To400"))  xsec = 92.59*1.165;
        if(file.Contains("ZJetsToNuNu_HT-400To600"))  xsec = 13.18*1.165;
        if(file.Contains("ZJetsToNuNu_HT-600To800"))  xsec = 3.257*1.165;
        if(file.Contains("ZJetsToNuNu_HT-800To1200"))  xsec = 1.49*1.165;
        if(file.Contains("ZJetsToNuNu_HT-1200To2500"))  xsec = 0.3419*1.165;
        if(file.Contains("ZJetsToNuNu_HT-2500ToInf"))  xsec = 0.005146*1.165;

        if(file.Contains("TTZToQQ"))                xsec = 0.5104;
        if(file.Contains("TTZToLLNuNu_M-10"))       xsec = 0.2432;
        if(file.Contains("TTWJetsToQQ"))            xsec = 0.4316;
        if(file.Contains("TTWJetsToLNu"))           xsec = 0.2149;
       
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Diboson
        if(file.Contains("WWTo2L2Nu"))   xsec = 12.178; //NNLO
        if(file.Contains("WWToLNuQQ"))   xsec = 49.997; //NNLO
        if(file.Contains("ttHTobb_M125"))   xsec = 0.2934;

        if(file.Contains("WZTo1L3Nu"))   xsec = 3.05;
        if(file.Contains("WZTo1L1Nu2Q"))   xsec = 10.73;
        if(file.Contains("WZTo2L2Q"))   xsec = 5.606;
        if(file.Contains("WZTo3LNu"))   xsec = 4.42965;
        if(file.Contains("VVTo2L2Nu"))   xsec = 11.95;
        if(file.Contains("ZZ_Tune"))   xsec = 16.523; //from twiki NLO

        // Calculated at 13 TeV in
        // https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt1314TeV
        // Higgs branching ratios from
        // https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR
        if(file.Contains("ZH_HToBB_ZToLL_M-125"))      xsec = 0.883*Htobb*0.033658;
        if(file.Contains("ZH_HToBB_ZToNuNu_M-125"))    xsec = 0.883*Htobb*0.2;
        if(file.Contains("WH_HToBB_WToLNu_M-125"))     xsec = 1.373*Htobb*(0.1071+0.1063+0.1138);
        if(file.Contains("ZH_HToBB_ZToNuNu_M125"))    xsec = 0.883*Htobb*0.2;
        if(file.Contains("WH_HToBB_WToLNu_M125"))     xsec = 1.373*Htobb*(0.1071+0.1063+0.1138);

    }
    if(xsec<=0) std::cout<<"BABYMAKER: Cross section not found for "<<file<<std::endl;

    return xsec;
  }

  float fractionNegWeights(const TString &file){
    float fneg(0.);

    // ttH, ttZ, ttW, ttGamma
    if(file.Contains("ttHJetTobb_M125_13TeV_amcatnloFXFX"))     fneg = 0.3515;
    if(file.Contains("TTZToQQ"))                                fneg = 0.2657;
    if(file.Contains("TTZToLLNuNu_M-10"))                       fneg = 0.2676;
    if(file.Contains("TTWJetsToQQ"))                            fneg = 0.2412;
    if(file.Contains("TTWJetsToLNu"))                           fneg = 0.2433;
    if(file.Contains("TTG"))                                    fneg = 0.342; // from MCM

    if(file.Contains("TTTT_TuneCUETP8M1_13TeV-amcatnlo"))       fneg = 0.41; // from MCM
    if(file.Contains("VVTo2L2Nu_13TeV_amcatnloFXFX"))       fneg = 0.20; // from MCM
    if(file.Contains("TTJets_Mtt-1000toInf"))                   fneg = 0.376996;

    // Single top
    if (file.Contains("ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8")) fneg = 0.1884;

    return fneg;
  }


  void signalCrossSection(int glu_mass, double &xsec, double &xsec_unc){
    if      (glu_mass == 595  ) { xsec = 0.       ; xsec_unc = 0.       ; return; } // we shouldn't have these points
    else if (glu_mass == 600. ) { xsec = 0.113E+02; xsec_unc = 8.62 /100; return; }
    else if (glu_mass == 605. ) { xsec = 0.108E+02; xsec_unc = 8.66 /100; return; }
    else if (glu_mass == 610. ) { xsec = 0.102E+02; xsec_unc = 8.7 /100; return; }
    else if (glu_mass == 615. ) { xsec = 0.974E+01; xsec_unc = 8.74 /100; return; }
    else if (glu_mass == 620. ) { xsec = 0.926E+01; xsec_unc = 8.78 /100; return; }
    else if (glu_mass == 625. ) { xsec = 0.881E+01; xsec_unc = 8.82 /100; return; }
    else if (glu_mass == 630. ) { xsec = 0.839E+01; xsec_unc = 8.86 /100; return; }
    else if (glu_mass == 635. ) { xsec = 0.799E+01; xsec_unc = 8.9 /100; return; }
    else if (glu_mass == 640. ) { xsec = 0.761E+01; xsec_unc = 8.94 /100; return; }
    else if (glu_mass == 645. ) { xsec = 0.725E+01; xsec_unc = 8.98 /100; return; }
    else if (glu_mass == 650. ) { xsec = 0.690E+01; xsec_unc = 9.02 /100; return; }
    else if (glu_mass == 655. ) { xsec = 0.658E+01; xsec_unc = 9.06 /100; return; }
    else if (glu_mass == 660. ) { xsec = 0.627E+01; xsec_unc = 9.1 /100; return; }
    else if (glu_mass == 665. ) { xsec = 0.598E+01; xsec_unc = 9.15 /100; return; }
    else if (glu_mass == 670. ) { xsec = 0.571E+01; xsec_unc = 9.19 /100; return; }
    else if (glu_mass == 675. ) { xsec = 0.544E+01; xsec_unc = 9.23 /100; return; }
    else if (glu_mass == 680. ) { xsec = 0.520E+01; xsec_unc = 9.27 /100; return; }
    else if (glu_mass == 685. ) { xsec = 0.496E+01; xsec_unc = 9.31 /100; return; }
    else if (glu_mass == 690. ) { xsec = 0.474E+01; xsec_unc = 9.35 /100; return; }
    else if (glu_mass == 695. ) { xsec = 0.452E+01; xsec_unc = 9.39 /100; return; }
    else if (glu_mass == 700. ) { xsec = 0.432E+01; xsec_unc = 9.43 /100; return; }
    else if (glu_mass == 705. ) { xsec = 0.413E+01; xsec_unc = 9.46 /100; return; }
    else if (glu_mass == 710. ) { xsec = 0.395E+01; xsec_unc = 9.5 /100; return; }
    else if (glu_mass == 715. ) { xsec = 0.377E+01; xsec_unc = 9.54 /100; return; }
    else if (glu_mass == 720. ) { xsec = 0.361E+01; xsec_unc = 9.58 /100; return; }
    else if (glu_mass == 725. ) { xsec = 0.345E+01; xsec_unc = 9.61 /100; return; }
    else if (glu_mass == 730. ) { xsec = 0.330E+01; xsec_unc = 9.65 /100; return; }
    else if (glu_mass == 735. ) { xsec = 0.316E+01; xsec_unc = 9.69 /100; return; }
    else if (glu_mass == 740. ) { xsec = 0.302E+01; xsec_unc = 9.72 /100; return; }
    else if (glu_mass == 745. ) { xsec = 0.289E+01; xsec_unc = 9.76 /100; return; }
    else if (glu_mass == 750. ) { xsec = 0.277E+01; xsec_unc = 9.8 /100; return; }
    else if (glu_mass == 755. ) { xsec = 0.265E+01; xsec_unc = 9.83 /100; return; }
    else if (glu_mass == 760. ) { xsec = 0.254E+01; xsec_unc = 9.87 /100; return; }
    else if (glu_mass == 765. ) { xsec = 0.243E+01; xsec_unc = 9.91 /100; return; }
    else if (glu_mass == 770. ) { xsec = 0.233E+01; xsec_unc = 9.94 /100; return; }
    else if (glu_mass == 775. ) { xsec = 0.223E+01; xsec_unc = 9.98 /100; return; }
    else if (glu_mass == 780. ) { xsec = 0.214E+01; xsec_unc = 10.01 /100; return; }
    else if (glu_mass == 785. ) { xsec = 0.205E+01; xsec_unc = 10.05 /100; return; }
    else if (glu_mass == 790. ) { xsec = 0.197E+01; xsec_unc = 10.09 /100; return; }
    else if (glu_mass == 795. ) { xsec = 0.188E+01; xsec_unc = 10.12 /100; return; }
    else if (glu_mass == 800. ) { xsec = 0.181E+01; xsec_unc = 10.16 /100; return; }
    else if (glu_mass == 805. ) { xsec = 0.173E+01; xsec_unc = 10.2 /100; return; }
    else if (glu_mass == 810. ) { xsec = 0.166E+01; xsec_unc = 10.23 /100; return; }
    else if (glu_mass == 815. ) { xsec = 0.160E+01; xsec_unc = 10.27 /100; return; }
    else if (glu_mass == 820. ) { xsec = 0.153E+01; xsec_unc = 10.31 /100; return; }
    else if (glu_mass == 825. ) { xsec = 0.147E+01; xsec_unc = 10.34 /100; return; }
    else if (glu_mass == 830. ) { xsec = 0.141E+01; xsec_unc = 10.38 /100; return; }
    else if (glu_mass == 835. ) { xsec = 0.136E+01; xsec_unc = 10.42 /100; return; }
    else if (glu_mass == 840. ) { xsec = 0.130E+01; xsec_unc = 10.45 /100; return; }
    else if (glu_mass == 845. ) { xsec = 0.125E+01; xsec_unc = 10.49 /100; return; }
    else if (glu_mass == 850. ) { xsec = 0.120E+01; xsec_unc = 10.53 /100; return; }
    else if (glu_mass == 855. ) { xsec = 0.115E+01; xsec_unc = 10.57 /100; return; }
    else if (glu_mass == 860. ) { xsec = 0.111E+01; xsec_unc = 10.6 /100; return; }
    else if (glu_mass == 865. ) { xsec = 0.107E+01; xsec_unc = 10.64 /100; return; }
    else if (glu_mass == 870. ) { xsec = 0.103E+01; xsec_unc = 10.68 /100; return; }
    else if (glu_mass == 875. ) { xsec = 0.986E+00; xsec_unc = 10.71 /100; return; }
    else if (glu_mass == 880. ) { xsec = 0.948E+00; xsec_unc = 10.75 /100; return; }
    else if (glu_mass == 885. ) { xsec = 0.912E+00; xsec_unc = 10.79 /100; return; }
    else if (glu_mass == 890. ) { xsec = 0.877E+00; xsec_unc = 10.82 /100; return; }
    else if (glu_mass == 895. ) { xsec = 0.844E+00; xsec_unc = 10.86 /100; return; }
    else if (glu_mass == 900. ) { xsec = 0.812E+00; xsec_unc = 10.89 /100; return; }
    else if (glu_mass == 905. ) { xsec = 0.781E+00; xsec_unc = 10.93 /100; return; }
    else if (glu_mass == 910. ) { xsec = 0.752E+00; xsec_unc = 10.97 /100; return; }
    else if (glu_mass == 915. ) { xsec = 0.723E+00; xsec_unc = 11.0 /100; return; }
    else if (glu_mass == 920. ) { xsec = 0.696E+00; xsec_unc = 11.04 /100; return; }
    else if (glu_mass == 925. ) { xsec = 0.670E+00; xsec_unc = 11.07 /100; return; }
    else if (glu_mass == 930. ) { xsec = 0.646E+00; xsec_unc = 11.11 /100; return; }
    else if (glu_mass == 935. ) { xsec = 0.622E+00; xsec_unc = 11.14 /100; return; }
    else if (glu_mass == 940. ) { xsec = 0.599E+00; xsec_unc = 11.18 /100; return; }
    else if (glu_mass == 945. ) { xsec = 0.577E+00; xsec_unc = 11.21 /100; return; }
    else if (glu_mass == 950. ) { xsec = 0.556E+00; xsec_unc = 11.25 /100; return; }
    else if (glu_mass == 955. ) { xsec = 0.535E+00; xsec_unc = 11.28 /100; return; }
    else if (glu_mass == 960. ) { xsec = 0.516E+00; xsec_unc = 11.32 /100; return; }
    else if (glu_mass == 965. ) { xsec = 0.497E+00; xsec_unc = 11.35 /100; return; }
    else if (glu_mass == 970. ) { xsec = 0.479E+00; xsec_unc = 11.39 /100; return; }
    else if (glu_mass == 975. ) { xsec = 0.462E+00; xsec_unc = 11.42 /100; return; }
    else if (glu_mass == 980. ) { xsec = 0.445E+00; xsec_unc = 11.46 /100; return; }
    else if (glu_mass == 985. ) { xsec = 0.430E+00; xsec_unc = 11.49 /100; return; }
    else if (glu_mass == 990. ) { xsec = 0.414E+00; xsec_unc = 11.53 /100; return; }
    else if (glu_mass == 995. ) { xsec = 0.399E+00; xsec_unc = 11.56 /100; return; }
    else if (glu_mass == 1000 ) { xsec = 0.385E+00; xsec_unc = 11.6 /100; return; }
    else if (glu_mass == 1005 ) { xsec = 0.372E+00; xsec_unc = 11.63 /100; return; }
    else if (glu_mass == 1010 ) { xsec = 0.359E+00; xsec_unc = 11.67 /100; return; }
    else if (glu_mass == 1015 ) { xsec = 0.346E+00; xsec_unc = 11.7 /100; return; }
    else if (glu_mass == 1020 ) { xsec = 0.334E+00; xsec_unc = 11.74 /100; return; }
    else if (glu_mass == 1025 ) { xsec = 0.322E+00; xsec_unc = 11.78 /100; return; }
    else if (glu_mass == 1030 ) { xsec = 0.311E+00; xsec_unc = 11.81 /100; return; }
    else if (glu_mass == 1035 ) { xsec = 0.300E+00; xsec_unc = 11.85 /100; return; }
    else if (glu_mass == 1040 ) { xsec = 0.290E+00; xsec_unc = 11.88 /100; return; }
    else if (glu_mass == 1045 ) { xsec = 0.280E+00; xsec_unc = 11.92 /100; return; }
    else if (glu_mass == 1050 ) { xsec = 0.270E+00; xsec_unc = 11.95 /100; return; }
    else if (glu_mass == 1055 ) { xsec = 0.261E+00; xsec_unc = 11.99 /100; return; }
    else if (glu_mass == 1060 ) { xsec = 0.252E+00; xsec_unc = 12.02 /100; return; }
    else if (glu_mass == 1065 ) { xsec = 0.243E+00; xsec_unc = 12.06 /100; return; }
    else if (glu_mass == 1070 ) { xsec = 0.235E+00; xsec_unc = 12.09 /100; return; }
    else if (glu_mass == 1075 ) { xsec = 0.227E+00; xsec_unc = 12.13 /100; return; }
    else if (glu_mass == 1080 ) { xsec = 0.219E+00; xsec_unc = 12.17 /100; return; }
    else if (glu_mass == 1085 ) { xsec = 0.212E+00; xsec_unc = 12.2 /100; return; }
    else if (glu_mass == 1090 ) { xsec = 0.205E+00; xsec_unc = 12.24 /100; return; }
    else if (glu_mass == 1095 ) { xsec = 0.198E+00; xsec_unc = 12.27 /100; return; }
    else if (glu_mass == 1100 ) { xsec = 0.191E+00; xsec_unc = 12.31 /100; return; }
    else if (glu_mass == 1105 ) { xsec = 0.185E+00; xsec_unc = 12.34 /100; return; }
    else if (glu_mass == 1110 ) { xsec = 0.179E+00; xsec_unc = 12.38 /100; return; }
    else if (glu_mass == 1115 ) { xsec = 0.173E+00; xsec_unc = 12.42 /100; return; }
    else if (glu_mass == 1120 ) { xsec = 0.167E+00; xsec_unc = 12.45 /100; return; }
    else if (glu_mass == 1125 ) { xsec = 0.162E+00; xsec_unc = 12.49 /100; return; }
    else if (glu_mass == 1130 ) { xsec = 0.156E+00; xsec_unc = 12.53 /100; return; }
    else if (glu_mass == 1135 ) { xsec = 0.151E+00; xsec_unc = 12.56 /100; return; }
    else if (glu_mass == 1140 ) { xsec = 0.146E+00; xsec_unc = 12.6 /100; return; }
    else if (glu_mass == 1145 ) { xsec = 0.141E+00; xsec_unc = 12.64 /100; return; }
    else if (glu_mass == 1150 ) { xsec = 0.137E+00; xsec_unc = 12.67 /100; return; }
    else if (glu_mass == 1155 ) { xsec = 0.132E+00; xsec_unc = 12.71 /100; return; }
    else if (glu_mass == 1160 ) { xsec = 0.128E+00; xsec_unc = 12.74 /100; return; }
    else if (glu_mass == 1165 ) { xsec = 0.124E+00; xsec_unc = 12.78 /100; return; }
    else if (glu_mass == 1170 ) { xsec = 0.120E+00; xsec_unc = 12.82 /100; return; }
    else if (glu_mass == 1175 ) { xsec = 0.116E+00; xsec_unc = 12.85 /100; return; }
    else if (glu_mass == 1180 ) { xsec = 0.112E+00; xsec_unc = 12.89 /100; return; }
    else if (glu_mass == 1185 ) { xsec = 0.109E+00; xsec_unc = 12.92 /100; return; }
    else if (glu_mass == 1190 ) { xsec = 0.105E+00; xsec_unc = 12.96 /100; return; }
    else if (glu_mass == 1195 ) { xsec = 0.102E+00; xsec_unc = 13.0 /100; return; }
    else if (glu_mass == 1200 ) { xsec = 0.985E-01; xsec_unc = 13.03 /100; return; }
    else if (glu_mass == 1205 ) { xsec = 0.953E-01; xsec_unc = 13.07 /100; return; }
    else if (glu_mass == 1210 ) { xsec = 0.923E-01; xsec_unc = 13.1 /100; return; }
    else if (glu_mass == 1215 ) { xsec = 0.894E-01; xsec_unc = 13.14 /100; return; }
    else if (glu_mass == 1220 ) { xsec = 0.866E-01; xsec_unc = 13.17 /100; return; }
    else if (glu_mass == 1225 ) { xsec = 0.838E-01; xsec_unc = 13.21 /100; return; }
    else if (glu_mass == 1230 ) { xsec = 0.812E-01; xsec_unc = 13.24 /100; return; }
    else if (glu_mass == 1235 ) { xsec = 0.786E-01; xsec_unc = 13.27 /100; return; }
    else if (glu_mass == 1240 ) { xsec = 0.762E-01; xsec_unc = 13.31 /100; return; }
    else if (glu_mass == 1245 ) { xsec = 0.738E-01; xsec_unc = 13.34 /100; return; }
    else if (glu_mass == 1250 ) { xsec = 0.715E-01; xsec_unc = 13.38 /100; return; }
    else if (glu_mass == 1255 ) { xsec = 0.692E-01; xsec_unc = 13.41 /100; return; }
    else if (glu_mass == 1260 ) { xsec = 0.671E-01; xsec_unc = 13.45 /100; return; }
    else if (glu_mass == 1265 ) { xsec = 0.650E-01; xsec_unc = 13.48 /100; return; }
    else if (glu_mass == 1270 ) { xsec = 0.630E-01; xsec_unc = 13.51 /100; return; }
    else if (glu_mass == 1275 ) { xsec = 0.610E-01; xsec_unc = 13.55 /100; return; }
    else if (glu_mass == 1280 ) { xsec = 0.591E-01; xsec_unc = 13.58 /100; return; }
    else if (glu_mass == 1285 ) { xsec = 0.573E-01; xsec_unc = 13.62 /100; return; }
    else if (glu_mass == 1290 ) { xsec = 0.556E-01; xsec_unc = 13.65 /100; return; }
    else if (glu_mass == 1295 ) { xsec = 0.539E-01; xsec_unc = 13.69 /100; return; }
    else if (glu_mass == 1300 ) { xsec = 0.522E-01; xsec_unc = 13.72 /100; return; }
    else if (glu_mass == 1305 ) { xsec = 0.506E-01; xsec_unc = 13.76 /100; return; }
    else if (glu_mass == 1310 ) { xsec = 0.491E-01; xsec_unc = 13.79 /100; return; }
    else if (glu_mass == 1315 ) { xsec = 0.476E-01; xsec_unc = 13.83 /100; return; }
    else if (glu_mass == 1320 ) { xsec = 0.461E-01; xsec_unc = 13.86 /100; return; }
    else if (glu_mass == 1325 ) { xsec = 0.447E-01; xsec_unc = 13.9 /100; return; }
    else if (glu_mass == 1330 ) { xsec = 0.434E-01; xsec_unc = 13.94 /100; return; }
    else if (glu_mass == 1335 ) { xsec = 0.421E-01; xsec_unc = 13.97 /100; return; }
    else if (glu_mass == 1340 ) { xsec = 0.408E-01; xsec_unc = 14.01 /100; return; }
    else if (glu_mass == 1345 ) { xsec = 0.396E-01; xsec_unc = 14.04 /100; return; }
    else if (glu_mass == 1350 ) { xsec = 0.384E-01; xsec_unc = 14.08 /100; return; }
    else if (glu_mass == 1355 ) { xsec = 0.372E-01; xsec_unc = 14.11 /100; return; }
    else if (glu_mass == 1360 ) { xsec = 0.361E-01; xsec_unc = 14.15 /100; return; }
    else if (glu_mass == 1365 ) { xsec = 0.350E-01; xsec_unc = 14.19 /100; return; }
    else if (glu_mass == 1370 ) { xsec = 0.340E-01; xsec_unc = 14.22 /100; return; }
    else if (glu_mass == 1375 ) { xsec = 0.330E-01; xsec_unc = 14.26 /100; return; }
    else if (glu_mass == 1380 ) { xsec = 0.320E-01; xsec_unc = 14.3 /100; return; }
    else if (glu_mass == 1385 ) { xsec = 0.310E-01; xsec_unc = 14.33 /100; return; }
    else if (glu_mass == 1390 ) { xsec = 0.301E-01; xsec_unc = 14.37 /100; return; }
    else if (glu_mass == 1395 ) { xsec = 0.292E-01; xsec_unc = 14.4 /100; return; }
    else if (glu_mass == 1400 ) { xsec = 0.284E-01; xsec_unc = 14.44 /100; return; }
    else if (glu_mass == 1405 ) { xsec = 0.275E-01; xsec_unc = 14.48 /100; return; }
    else if (glu_mass == 1410 ) { xsec = 0.267E-01; xsec_unc = 14.51 /100; return; }
    else if (glu_mass == 1415 ) { xsec = 0.259E-01; xsec_unc = 14.55 /100; return; }
    else if (glu_mass == 1420 ) { xsec = 0.252E-01; xsec_unc = 14.59 /100; return; }
    else if (glu_mass == 1425 ) { xsec = 0.244E-01; xsec_unc = 14.63 /100; return; }
    else if (glu_mass == 1430 ) { xsec = 0.237E-01; xsec_unc = 14.66 /100; return; }
    else if (glu_mass == 1435 ) { xsec = 0.230E-01; xsec_unc = 14.7 /100; return; }
    else if (glu_mass == 1440 ) { xsec = 0.224E-01; xsec_unc = 14.74 /100; return; }
    else if (glu_mass == 1445 ) { xsec = 0.217E-01; xsec_unc = 14.77 /100; return; }
    else if (glu_mass == 1450 ) { xsec = 0.211E-01; xsec_unc = 14.81 /100; return; }
    else if (glu_mass == 1455 ) { xsec = 0.205E-01; xsec_unc = 14.85 /100; return; }
    else if (glu_mass == 1460 ) { xsec = 0.199E-01; xsec_unc = 14.88 /100; return; }
    else if (glu_mass == 1465 ) { xsec = 0.193E-01; xsec_unc = 14.92 /100; return; }
    else if (glu_mass == 1470 ) { xsec = 0.187E-01; xsec_unc = 14.96 /100; return; }
    else if (glu_mass == 1475 ) { xsec = 0.182E-01; xsec_unc = 15.0 /100; return; }
    else if (glu_mass == 1480 ) { xsec = 0.177E-01; xsec_unc = 15.03 /100; return; }
    else if (glu_mass == 1485 ) { xsec = 0.172E-01; xsec_unc = 15.07 /100; return; }
    else if (glu_mass == 1490 ) { xsec = 0.167E-01; xsec_unc = 15.11 /100; return; }
    else if (glu_mass == 1495 ) { xsec = 0.162E-01; xsec_unc = 15.15 /100; return; }
    else if (glu_mass == 1500 ) { xsec = 0.157E-01; xsec_unc = 15.18 /100; return; }
    else if (glu_mass == 1505 ) { xsec = 0.153E-01; xsec_unc = 15.22 /100; return; }
    else if (glu_mass == 1510 ) { xsec = 0.148E-01; xsec_unc = 15.26 /100; return; }
    else if (glu_mass == 1515 ) { xsec = 0.144E-01; xsec_unc = 15.3 /100; return; }
    else if (glu_mass == 1520 ) { xsec = 0.140E-01; xsec_unc = 15.33 /100; return; }
    else if (glu_mass == 1525 ) { xsec = 0.136E-01; xsec_unc = 15.37 /100; return; }
    else if (glu_mass == 1530 ) { xsec = 0.132E-01; xsec_unc = 15.41 /100; return; }
    else if (glu_mass == 1535 ) { xsec = 0.128E-01; xsec_unc = 15.45 /100; return; }
    else if (glu_mass == 1540 ) { xsec = 0.125E-01; xsec_unc = 15.48 /100; return; }
    else if (glu_mass == 1545 ) { xsec = 0.121E-01; xsec_unc = 15.52 /100; return; }
    else if (glu_mass == 1550 ) { xsec = 0.118E-01; xsec_unc = 15.56 /100; return; }
    else if (glu_mass == 1555 ) { xsec = 0.115E-01; xsec_unc = 15.6 /100; return; }
    else if (glu_mass == 1560 ) { xsec = 0.111E-01; xsec_unc = 15.64 /100; return; }
    else if (glu_mass == 1565 ) { xsec = 0.108E-01; xsec_unc = 15.67 /100; return; }
    else if (glu_mass == 1570 ) { xsec = 0.105E-01; xsec_unc = 15.71 /100; return; }
    else if (glu_mass == 1575 ) { xsec = 0.102E-01; xsec_unc = 15.75 /100; return; }
    else if (glu_mass == 1580 ) { xsec = 0.993E-02; xsec_unc = 15.79 /100; return; }
    else if (glu_mass == 1585 ) { xsec = 0.966E-02; xsec_unc = 15.83 /100; return; }
    else if (glu_mass == 1590 ) { xsec = 0.939E-02; xsec_unc = 15.87 /100; return; }
    else if (glu_mass == 1595 ) { xsec = 0.912E-02; xsec_unc = 15.9 /100; return; }
    else if (glu_mass == 1600 ) { xsec = 0.887E-02; xsec_unc = 15.94 /100; return; }
    else if (glu_mass == 1605 ) { xsec = 0.862E-02; xsec_unc = 15.98 /100; return; }
    else if (glu_mass == 1610 ) { xsec = 0.838E-02; xsec_unc = 16.02 /100; return; }
    else if (glu_mass == 1615 ) { xsec = 0.815E-02; xsec_unc = 16.06 /100; return; }
    else if (glu_mass == 1620 ) { xsec = 0.792E-02; xsec_unc = 16.1 /100; return; }
    else if (glu_mass == 1625 ) { xsec = 0.770E-02; xsec_unc = 16.13 /100; return; }
    else if (glu_mass == 1630 ) { xsec = 0.749E-02; xsec_unc = 16.17 /100; return; }
    else if (glu_mass == 1635 ) { xsec = 0.728E-02; xsec_unc = 16.21 /100; return; }
    else if (glu_mass == 1640 ) { xsec = 0.708E-02; xsec_unc = 16.25 /100; return; }
    else if (glu_mass == 1645 ) { xsec = 0.689E-02; xsec_unc = 16.29 /100; return; }
    else if (glu_mass == 1650 ) { xsec = 0.670E-02; xsec_unc = 16.33 /100; return; }
    else if (glu_mass == 1655 ) { xsec = 0.651E-02; xsec_unc = 16.37 /100; return; }
    else if (glu_mass == 1660 ) { xsec = 0.633E-02; xsec_unc = 16.41 /100; return; }
    else if (glu_mass == 1665 ) { xsec = 0.616E-02; xsec_unc = 16.45 /100; return; }
    else if (glu_mass == 1670 ) { xsec = 0.599E-02; xsec_unc = 16.49 /100; return; }
    else if (glu_mass == 1675 ) { xsec = 0.583E-02; xsec_unc = 16.53 /100; return; }
    else if (glu_mass == 1680 ) { xsec = 0.567E-02; xsec_unc = 16.56 /100; return; }
    else if (glu_mass == 1685 ) { xsec = 0.551E-02; xsec_unc = 16.6 /100; return; }
    else if (glu_mass == 1690 ) { xsec = 0.536E-02; xsec_unc = 16.64 /100; return; }
    else if (glu_mass == 1695 ) { xsec = 0.521E-02; xsec_unc = 16.68 /100; return; }
    else if (glu_mass == 1700 ) { xsec = 0.507E-02; xsec_unc = 16.72 /100; return; }
    else if (glu_mass == 1705 ) { xsec = 0.493E-02; xsec_unc = 16.76 /100; return; }
    else if (glu_mass == 1710 ) { xsec = 0.480E-02; xsec_unc = 16.81 /100; return; }
    else if (glu_mass == 1715 ) { xsec = 0.467E-02; xsec_unc = 16.85 /100; return; }
    else if (glu_mass == 1720 ) { xsec = 0.454E-02; xsec_unc = 16.89 /100; return; }
    else if (glu_mass == 1725 ) { xsec = 0.442E-02; xsec_unc = 16.93 /100; return; }
    else if (glu_mass == 1730 ) { xsec = 0.430E-02; xsec_unc = 16.97 /100; return; }
    else if (glu_mass == 1735 ) { xsec = 0.418E-02; xsec_unc = 17.01 /100; return; }
    else if (glu_mass == 1740 ) { xsec = 0.407E-02; xsec_unc = 17.05 /100; return; }
    else if (glu_mass == 1745 ) { xsec = 0.396E-02; xsec_unc = 17.09 /100; return; }
    else if (glu_mass == 1750 ) { xsec = 0.385E-02; xsec_unc = 17.13 /100; return; }
    else if (glu_mass == 1755 ) { xsec = 0.375E-02; xsec_unc = 17.18 /100; return; }
    else if (glu_mass == 1760 ) { xsec = 0.365E-02; xsec_unc = 17.22 /100; return; }
    else if (glu_mass == 1765 ) { xsec = 0.355E-02; xsec_unc = 17.26 /100; return; }
    else if (glu_mass == 1770 ) { xsec = 0.345E-02; xsec_unc = 17.3 /100; return; }
    else if (glu_mass == 1775 ) { xsec = 0.336E-02; xsec_unc = 17.34 /100; return; }
    else if (glu_mass == 1780 ) { xsec = 0.327E-02; xsec_unc = 17.39 /100; return; }
    else if (glu_mass == 1785 ) { xsec = 0.318E-02; xsec_unc = 17.43 /100; return; }
    else if (glu_mass == 1790 ) { xsec = 0.310E-02; xsec_unc = 17.47 /100; return; }
    else if (glu_mass == 1795 ) { xsec = 0.301E-02; xsec_unc = 17.51 /100; return; }
    else if (glu_mass == 1800 ) { xsec = 0.293E-02; xsec_unc = 17.56 /100; return; }
    else if (glu_mass == 1805 ) { xsec = 0.286E-02; xsec_unc = 17.6 /100; return; }
    else if (glu_mass == 1810 ) { xsec = 0.278E-02; xsec_unc = 17.64 /100; return; }
    else if (glu_mass == 1815 ) { xsec = 0.271E-02; xsec_unc = 17.69 /100; return; }
    else if (glu_mass == 1820 ) { xsec = 0.263E-02; xsec_unc = 17.73 /100; return; }
    else if (glu_mass == 1825 ) { xsec = 0.256E-02; xsec_unc = 17.77 /100; return; }
    else if (glu_mass == 1830 ) { xsec = 0.249E-02; xsec_unc = 17.82 /100; return; }
    else if (glu_mass == 1835 ) { xsec = 0.243E-02; xsec_unc = 17.86 /100; return; }
    else if (glu_mass == 1840 ) { xsec = 0.236E-02; xsec_unc = 17.9 /100; return; }
    else if (glu_mass == 1845 ) { xsec = 0.230E-02; xsec_unc = 17.95 /100; return; }
    else if (glu_mass == 1850 ) { xsec = 0.224E-02; xsec_unc = 17.99 /100; return; }
    else if (glu_mass == 1855 ) { xsec = 0.218E-02; xsec_unc = 18.04 /100; return; }
    else if (glu_mass == 1860 ) { xsec = 0.212E-02; xsec_unc = 18.08 /100; return; }
    else if (glu_mass == 1865 ) { xsec = 0.207E-02; xsec_unc = 18.13 /100; return; }
    else if (glu_mass == 1870 ) { xsec = 0.201E-02; xsec_unc = 18.17 /100; return; }
    else if (glu_mass == 1875 ) { xsec = 0.196E-02; xsec_unc = 18.22 /100; return; }
    else if (glu_mass == 1880 ) { xsec = 0.191E-02; xsec_unc = 18.26 /100; return; }
    else if (glu_mass == 1885 ) { xsec = 0.186E-02; xsec_unc = 18.31 /100; return; }
    else if (glu_mass == 1890 ) { xsec = 0.181E-02; xsec_unc = 18.35 /100; return; }
    else if (glu_mass == 1895 ) { xsec = 0.176E-02; xsec_unc = 18.4 /100; return; }
    else if (glu_mass == 1900 ) { xsec = 0.171E-02; xsec_unc = 18.45 /100; return; }
    else if (glu_mass == 1905 ) { xsec = 0.167E-02; xsec_unc = 18.49 /100; return; }
    else if (glu_mass == 1910 ) { xsec = 0.163E-02; xsec_unc = 18.54 /100; return; }
    else if (glu_mass == 1915 ) { xsec = 0.158E-02; xsec_unc = 18.59 /100; return; }
    else if (glu_mass == 1920 ) { xsec = 0.154E-02; xsec_unc = 18.63 /100; return; }
    else if (glu_mass == 1925 ) { xsec = 0.150E-02; xsec_unc = 18.68 /100; return; }
    else if (glu_mass == 1930 ) { xsec = 0.146E-02; xsec_unc = 18.73 /100; return; }
    else if (glu_mass == 1935 ) { xsec = 0.142E-02; xsec_unc = 18.78 /100; return; }
    else if (glu_mass == 1940 ) { xsec = 0.139E-02; xsec_unc = 18.82 /100; return; }
    else if (glu_mass == 1945 ) { xsec = 0.135E-02; xsec_unc = 18.87 /100; return; }
    else if (glu_mass == 1950 ) { xsec = 0.131E-02; xsec_unc = 18.92 /100; return; }
    else if (glu_mass == 1955 ) { xsec = 0.128E-02; xsec_unc = 18.97 /100; return; }
    else if (glu_mass == 1960 ) { xsec = 0.125E-02; xsec_unc = 19.02 /100; return; }
    else if (glu_mass == 1965 ) { xsec = 0.121E-02; xsec_unc = 19.07 /100; return; }
    else if (glu_mass == 1970 ) { xsec = 0.118E-02; xsec_unc = 19.12 /100; return; }
    else if (glu_mass == 1975 ) { xsec = 0.115E-02; xsec_unc = 19.17 /100; return; }
    else if (glu_mass == 1980 ) { xsec = 0.112E-02; xsec_unc = 19.22 /100; return; }
    else if (glu_mass == 1985 ) { xsec = 0.109E-02; xsec_unc = 19.27 /100; return; }
    else if (glu_mass == 1990 ) { xsec = 0.106E-02; xsec_unc = 19.32 /100; return; }
    else if (glu_mass == 1995 ) { xsec = 0.104E-02; xsec_unc = 19.37 /100; return; }
    else if (glu_mass == 2000 ) { xsec = 0.101E-02; xsec_unc = 19.42 /100; return; }
    else if (glu_mass == 2005 ) { xsec = 0.983E-03; xsec_unc = 19.48 /100; return; }
    else if (glu_mass == 2010 ) { xsec = 0.957E-03; xsec_unc = 19.53 /100; return; }
    else if (glu_mass == 2015 ) { xsec = 0.933E-03; xsec_unc = 19.58 /100; return; }
    else if (glu_mass == 2020 ) { xsec = 0.908E-03; xsec_unc = 19.64 /100; return; }
    else if (glu_mass == 2025 ) { xsec = 0.885E-03; xsec_unc = 19.69 /100; return; }
    else if (glu_mass == 2030 ) { xsec = 0.862E-03; xsec_unc = 19.74 /100; return; }
    else if (glu_mass == 2035 ) { xsec = 0.840E-03; xsec_unc = 19.8 /100; return; }
    else if (glu_mass == 2040 ) { xsec = 0.818E-03; xsec_unc = 19.85 /100; return; }
    else if (glu_mass == 2045 ) { xsec = 0.797E-03; xsec_unc = 19.91 /100; return; }
    else if (glu_mass == 2050 ) { xsec = 0.776E-03; xsec_unc = 19.96 /100; return; }
    else if (glu_mass == 2055 ) { xsec = 0.756E-03; xsec_unc = 20.02 /100; return; }
    else if (glu_mass == 2060 ) { xsec = 0.737E-03; xsec_unc = 20.07 /100; return; }
    else if (glu_mass == 2065 ) { xsec = 0.718E-03; xsec_unc = 20.13 /100; return; }
    else if (glu_mass == 2070 ) { xsec = 0.699E-03; xsec_unc = 20.19 /100; return; }
    else if (glu_mass == 2075 ) { xsec = 0.681E-03; xsec_unc = 20.25 /100; return; }
    else if (glu_mass == 2080 ) { xsec = 0.664E-03; xsec_unc = 20.3 /100; return; }
    else if (glu_mass == 2085 ) { xsec = 0.647E-03; xsec_unc = 20.36 /100; return; }
    else if (glu_mass == 2090 ) { xsec = 0.630E-03; xsec_unc = 20.42 /100; return; }
    else if (glu_mass == 2095 ) { xsec = 0.614E-03; xsec_unc = 20.48 /100; return; }
    else if (glu_mass == 2100 ) { xsec = 0.598E-03; xsec_unc = 20.54 /100; return; }
    else if (glu_mass == 2105 ) { xsec = 0.583E-03; xsec_unc = 20.6 /100; return; }
    else if (glu_mass == 2110 ) { xsec = 0.568E-03; xsec_unc = 20.66 /100; return; }
    else if (glu_mass == 2115 ) { xsec = 0.553E-03; xsec_unc = 20.72 /100; return; }
    else if (glu_mass == 2120 ) { xsec = 0.539E-03; xsec_unc = 20.78 /100; return; }
    else if (glu_mass == 2125 ) { xsec = 0.525E-03; xsec_unc = 20.84 /100; return; }
    else if (glu_mass == 2130 ) { xsec = 0.512E-03; xsec_unc = 20.9 /100; return; }
    else if (glu_mass == 2135 ) { xsec = 0.499E-03; xsec_unc = 20.97 /100; return; }
    else if (glu_mass == 2140 ) { xsec = 0.486E-03; xsec_unc = 21.03 /100; return; }
    else if (glu_mass == 2145 ) { xsec = 0.473E-03; xsec_unc = 21.09 /100; return; }
    else if (glu_mass == 2150 ) { xsec = 0.461E-03; xsec_unc = 21.16 /100; return; }
    else if (glu_mass == 2155 ) { xsec = 0.449E-03; xsec_unc = 21.22 /100; return; }
    else if (glu_mass == 2160 ) { xsec = 0.438E-03; xsec_unc = 21.29 /100; return; }
    else if (glu_mass == 2165 ) { xsec = 0.427E-03; xsec_unc = 21.35 /100; return; }
    else if (glu_mass == 2170 ) { xsec = 0.416E-03; xsec_unc = 21.42 /100; return; }
    else if (glu_mass == 2175 ) { xsec = 0.405E-03; xsec_unc = 21.48 /100; return; }
    else if (glu_mass == 2180 ) { xsec = 0.395E-03; xsec_unc = 21.55 /100; return; }
    else if (glu_mass == 2185 ) { xsec = 0.385E-03; xsec_unc = 21.62 /100; return; }
    else if (glu_mass == 2190 ) { xsec = 0.375E-03; xsec_unc = 21.69 /100; return; }
    else if (glu_mass == 2195 ) { xsec = 0.365E-03; xsec_unc = 21.76 /100; return; }
    else if (glu_mass == 2200 ) { xsec = 0.356E-03; xsec_unc = 21.83 /100; return; }
    else if (glu_mass == 2205 ) { xsec = 0.347E-03; xsec_unc = 21.9 /100; return; }
    else if (glu_mass == 2210 ) { xsec = 0.338E-03; xsec_unc = 21.97 /100; return; }
    else if (glu_mass == 2215 ) { xsec = 0.330E-03; xsec_unc = 22.04 /100; return; }
    else if (glu_mass == 2220 ) { xsec = 0.321E-03; xsec_unc = 22.11 /100; return; }
    else if (glu_mass == 2225 ) { xsec = 0.313E-03; xsec_unc = 22.18 /100; return; }
    else if (glu_mass == 2230 ) { xsec = 0.305E-03; xsec_unc = 22.25 /100; return; }
    else if (glu_mass == 2235 ) { xsec = 0.297E-03; xsec_unc = 22.33 /100; return; }
    else if (glu_mass == 2240 ) { xsec = 0.290E-03; xsec_unc = 22.4 /100; return; }
    else if (glu_mass == 2245 ) { xsec = 0.283E-03; xsec_unc = 22.47 /100; return; }
    else if (glu_mass == 2250 ) { xsec = 0.275E-03; xsec_unc = 22.55 /100; return; }
    else if (glu_mass == 2255 ) { xsec = 0.268E-03; xsec_unc = 22.63 /100; return; }
    else if (glu_mass == 2260 ) { xsec = 0.262E-03; xsec_unc = 22.7 /100; return; }
    else if (glu_mass == 2265 ) { xsec = 0.255E-03; xsec_unc = 22.78 /100; return; }
    else if (glu_mass == 2270 ) { xsec = 0.248E-03; xsec_unc = 22.86 /100; return; }
    else if (glu_mass == 2275 ) { xsec = 0.242E-03; xsec_unc = 22.94 /100; return; }
    else if (glu_mass == 2280 ) { xsec = 0.236E-03; xsec_unc = 23.02 /100; return; }
    else if (glu_mass == 2285 ) { xsec = 0.230E-03; xsec_unc = 23.1 /100; return; }
    else if (glu_mass == 2290 ) { xsec = 0.224E-03; xsec_unc = 23.18 /100; return; }
    else if (glu_mass == 2295 ) { xsec = 0.219E-03; xsec_unc = 23.26 /100; return; }
    else if (glu_mass == 2300 ) { xsec = 0.213E-03; xsec_unc = 23.34 /100; return; }
    else if (glu_mass == 2305 ) { xsec = 0.208E-03; xsec_unc = 23.43 /100; return; }
    else if (glu_mass == 2310 ) { xsec = 0.202E-03; xsec_unc = 23.51 /100; return; }
    else if (glu_mass == 2315 ) { xsec = 0.197E-03; xsec_unc = 23.6 /100; return; }
    else if (glu_mass == 2320 ) { xsec = 0.192E-03; xsec_unc = 23.68 /100; return; }
    else if (glu_mass == 2325 ) { xsec = 0.187E-03; xsec_unc = 23.77 /100; return; }
    else if (glu_mass == 2330 ) { xsec = 0.183E-03; xsec_unc = 23.86 /100; return; }
    else if (glu_mass == 2335 ) { xsec = 0.178E-03; xsec_unc = 23.95 /100; return; }
    else if (glu_mass == 2340 ) { xsec = 0.174E-03; xsec_unc = 24.04 /100; return; }
    else if (glu_mass == 2345 ) { xsec = 0.169E-03; xsec_unc = 24.13 /100; return; }
    else if (glu_mass == 2350 ) { xsec = 0.165E-03; xsec_unc = 24.22 /100; return; }
    else if (glu_mass == 2355 ) { xsec = 0.161E-03; xsec_unc = 24.31 /100; return; }
    else if (glu_mass == 2360 ) { xsec = 0.157E-03; xsec_unc = 24.41 /100; return; }
    else if (glu_mass == 2365 ) { xsec = 0.153E-03; xsec_unc = 24.5 /100; return; }
    else if (glu_mass == 2370 ) { xsec = 0.149E-03; xsec_unc = 24.6 /100; return; }
    else if (glu_mass == 2375 ) { xsec = 0.145E-03; xsec_unc = 24.7 /100; return; }
    else if (glu_mass == 2380 ) { xsec = 0.142E-03; xsec_unc = 24.79 /100; return; }
    else if (glu_mass == 2385 ) { xsec = 0.138E-03; xsec_unc = 24.89 /100; return; }
    else if (glu_mass == 2390 ) { xsec = 0.134E-03; xsec_unc = 24.99 /100; return; }
    else if (glu_mass == 2395 ) { xsec = 0.131E-03; xsec_unc = 25.09 /100; return; }
    else if (glu_mass == 2400 ) { xsec = 0.128E-03; xsec_unc = 25.19 /100; return; }
    else if (glu_mass == 2405 ) { xsec = 0.125E-03; xsec_unc = 25.3 /100; return; }
    else if (glu_mass == 2410 ) { xsec = 0.121E-03; xsec_unc = 25.4 /100; return; }
    else if (glu_mass == 2415 ) { xsec = 0.118E-03; xsec_unc = 25.5 /100; return; }
    else if (glu_mass == 2420 ) { xsec = 0.115E-03; xsec_unc = 25.61 /100; return; }
    else if (glu_mass == 2425 ) { xsec = 0.113E-03; xsec_unc = 25.71 /100; return; }
    else if (glu_mass == 2430 ) { xsec = 0.110E-03; xsec_unc = 25.82 /100; return; }
    else if (glu_mass == 2435 ) { xsec = 0.107E-03; xsec_unc = 25.93 /100; return; }
    else if (glu_mass == 2440 ) { xsec = 0.104E-03; xsec_unc = 26.04 /100; return; }
    else if (glu_mass == 2445 ) { xsec = 0.102E-03; xsec_unc = 26.15 /100; return; }
    else if (glu_mass == 2450 ) { xsec = 0.991E-04; xsec_unc = 26.26 /100; return; }
    else if (glu_mass == 2455 ) { xsec = 0.966E-04; xsec_unc = 26.37 /100; return; }
    else if (glu_mass == 2460 ) { xsec = 0.941E-04; xsec_unc = 26.49 /100; return; }
    else if (glu_mass == 2465 ) { xsec = 0.918E-04; xsec_unc = 26.6 /100; return; }
    else if (glu_mass == 2470 ) { xsec = 0.895E-04; xsec_unc = 26.72 /100; return; }
    else if (glu_mass == 2475 ) { xsec = 0.872E-04; xsec_unc = 26.84 /100; return; }
    else if (glu_mass == 2480 ) { xsec = 0.850E-04; xsec_unc = 26.96 /100; return; }
    else if (glu_mass == 2485 ) { xsec = 0.829E-04; xsec_unc = 27.08 /100; return; }
    else if (glu_mass == 2490 ) { xsec = 0.808E-04; xsec_unc = 27.2 /100; return; }
    else if (glu_mass == 2495 ) { xsec = 0.788E-04; xsec_unc = 27.33 /100; return; }
    else if (glu_mass == 2500 ) { xsec = 0.768E-04; xsec_unc = 27.45 /100; return; }
    else if (glu_mass == 2505 ) { xsec = 0.749E-04; xsec_unc = 27.58 /100; return; }
    else if (glu_mass == 2510 ) { xsec = 0.730E-04; xsec_unc = 27.71 /100; return; }
    else if (glu_mass == 2515 ) { xsec = 0.712E-04; xsec_unc = 27.84 /100; return; }
    else if (glu_mass == 2520 ) { xsec = 0.694E-04; xsec_unc = 27.97 /100; return; }
    else if (glu_mass == 2525 ) { xsec = 0.677E-04; xsec_unc = 28.11 /100; return; }
    else if (glu_mass == 2530 ) { xsec = 0.660E-04; xsec_unc = 28.24 /100; return; }
    else if (glu_mass == 2535 ) { xsec = 0.643E-04; xsec_unc = 28.38 /100; return; }
    else if (glu_mass == 2540 ) { xsec = 0.627E-04; xsec_unc = 28.52 /100; return; }
    else if (glu_mass == 2545 ) { xsec = 0.611E-04; xsec_unc = 28.66 /100; return; }
    else if (glu_mass == 2550 ) { xsec = 0.596E-04; xsec_unc = 28.8 /100; return; }
    else if (glu_mass == 2555 ) { xsec = 0.581E-04; xsec_unc = 28.94 /100; return; }
    else if (glu_mass == 2560 ) { xsec = 0.566E-04; xsec_unc = 29.09 /100; return; }
    else if (glu_mass == 2565 ) { xsec = 0.552E-04; xsec_unc = 29.23 /100; return; }
    else if (glu_mass == 2570 ) { xsec = 0.538E-04; xsec_unc = 29.38 /100; return; }
    else if (glu_mass == 2575 ) { xsec = 0.525E-04; xsec_unc = 29.53 /100; return; }
    else if (glu_mass == 2580 ) { xsec = 0.512E-04; xsec_unc = 29.68 /100; return; }
    else if (glu_mass == 2585 ) { xsec = 0.499E-04; xsec_unc = 29.84 /100; return; }
    else if (glu_mass == 2590 ) { xsec = 0.486E-04; xsec_unc = 29.99 /100; return; }
    else if (glu_mass == 2595 ) { xsec = 0.474E-04; xsec_unc = 30.15 /100; return; }
    else if (glu_mass == 2600 ) { xsec = 0.462E-04; xsec_unc = 30.31 /100; return; }
    else if (glu_mass == 2605 ) { xsec = 0.451E-04; xsec_unc = 30.47 /100; return; }
    else if (glu_mass == 2610 ) { xsec = 0.439E-04; xsec_unc = 30.63 /100; return; }
    else if (glu_mass == 2615 ) { xsec = 0.428E-04; xsec_unc = 30.8 /100; return; }
    else if (glu_mass == 2620 ) { xsec = 0.418E-04; xsec_unc = 30.97 /100; return; }
    else if (glu_mass == 2625 ) { xsec = 0.407E-04; xsec_unc = 31.13 /100; return; }
    else if (glu_mass == 2630 ) { xsec = 0.397E-04; xsec_unc = 31.3 /100; return; }
    else if (glu_mass == 2635 ) { xsec = 0.387E-04; xsec_unc = 31.48 /100; return; }
    else if (glu_mass == 2640 ) { xsec = 0.377E-04; xsec_unc = 31.65 /100; return; }
    else if (glu_mass == 2645 ) { xsec = 0.368E-04; xsec_unc = 31.83 /100; return; }
    else if (glu_mass == 2650 ) { xsec = 0.359E-04; xsec_unc = 32.01 /100; return; }
    else if (glu_mass == 2655 ) { xsec = 0.350E-04; xsec_unc = 32.19 /100; return; }
    else if (glu_mass == 2660 ) { xsec = 0.341E-04; xsec_unc = 32.37 /100; return; }
    else if (glu_mass == 2665 ) { xsec = 0.332E-04; xsec_unc = 32.56 /100; return; }
    else if (glu_mass == 2670 ) { xsec = 0.324E-04; xsec_unc = 32.74 /100; return; }
    else if (glu_mass == 2675 ) { xsec = 0.316E-04; xsec_unc = 32.93 /100; return; }
    else if (glu_mass == 2680 ) { xsec = 0.308E-04; xsec_unc = 33.12 /100; return; }
    else if (glu_mass == 2685 ) { xsec = 0.300E-04; xsec_unc = 33.32 /100; return; }
    else if (glu_mass == 2690 ) { xsec = 0.293E-04; xsec_unc = 33.52 /100; return; }
    else if (glu_mass == 2695 ) { xsec = 0.285E-04; xsec_unc = 33.71 /100; return; }
    else if (glu_mass == 2700 ) { xsec = 0.278E-04; xsec_unc = 33.92 /100; return; }
    else if (glu_mass == 2705 ) { xsec = 0.271E-04; xsec_unc = 34.13 /100; return; }
    else if (glu_mass == 2710 ) { xsec = 0.265E-04; xsec_unc = 34.34 /100; return; }
    else if (glu_mass == 2715 ) { xsec = 0.258E-04; xsec_unc = 34.56 /100; return; }
    else if (glu_mass == 2720 ) { xsec = 0.251E-04; xsec_unc = 34.77 /100; return; }
    else if (glu_mass == 2725 ) { xsec = 0.245E-04; xsec_unc = 34.99 /100; return; }
    else if (glu_mass == 2730 ) { xsec = 0.239E-04; xsec_unc = 35.22 /100; return; }
    else if (glu_mass == 2735 ) { xsec = 0.233E-04; xsec_unc = 35.44 /100; return; }
    else if (glu_mass == 2740 ) { xsec = 0.227E-04; xsec_unc = 35.67 /100; return; }
    else if (glu_mass == 2745 ) { xsec = 0.221E-04; xsec_unc = 35.89 /100; return; }
    else if (glu_mass == 2750 ) { xsec = 0.216E-04; xsec_unc = 36.12 /100; return; }
    else if (glu_mass == 2755 ) { xsec = 0.211E-04; xsec_unc = 36.35 /100; return; }
    else if (glu_mass == 2760 ) { xsec = 0.205E-04; xsec_unc = 36.59 /100; return; }
    else if (glu_mass == 2765 ) { xsec = 0.200E-04; xsec_unc = 36.82 /100; return; }
    else if (glu_mass == 2770 ) { xsec = 0.195E-04; xsec_unc = 37.05 /100; return; }
    else if (glu_mass == 2775 ) { xsec = 0.190E-04; xsec_unc = 37.29 /100; return; }
    else if (glu_mass == 2780 ) { xsec = 0.185E-04; xsec_unc = 37.53 /100; return; }
    else if (glu_mass == 2785 ) { xsec = 0.181E-04; xsec_unc = 37.76 /100; return; }
    else if (glu_mass == 2790 ) { xsec = 0.176E-04; xsec_unc = 38.0 /100; return; }
    else if (glu_mass == 2795 ) { xsec = 0.172E-04; xsec_unc = 38.24 /100; return; }
    else if (glu_mass == 2800 ) { xsec = 0.168E-04; xsec_unc = 38.48 /100; return; }
    else if (glu_mass == 2805 ) { xsec = 0.163E-04; xsec_unc = 38.72 /100; return; }
    else if (glu_mass == 2810 ) { xsec = 0.159E-04; xsec_unc = 38.95 /100; return; }
    else if (glu_mass == 2815 ) { xsec = 0.155E-04; xsec_unc = 39.19 /100; return; }
    else if (glu_mass == 2820 ) { xsec = 0.151E-04; xsec_unc = 39.43 /100; return; }
    else if (glu_mass == 2825 ) { xsec = 0.148E-04; xsec_unc = 39.66 /100; return; }
    else if (glu_mass == 2830 ) { xsec = 0.144E-04; xsec_unc = 39.9 /100; return; }
    else if (glu_mass == 2835 ) { xsec = 0.140E-04; xsec_unc = 40.14 /100; return; }
    else if (glu_mass == 2840 ) { xsec = 0.137E-04; xsec_unc = 40.38 /100; return; }
    else if (glu_mass == 2845 ) { xsec = 0.133E-04; xsec_unc = 40.62 /100; return; }
    else if (glu_mass == 2850 ) { xsec = 0.130E-04; xsec_unc = 40.86 /100; return; }
    else if (glu_mass == 2855 ) { xsec = 0.127E-04; xsec_unc = 41.1 /100; return; }
    else if (glu_mass == 2860 ) { xsec = 0.124E-04; xsec_unc = 41.34 /100; return; }
    else if (glu_mass == 2865 ) { xsec = 0.121E-04; xsec_unc = 41.59 /100; return; }
    else if (glu_mass == 2870 ) { xsec = 0.118E-04; xsec_unc = 41.83 /100; return; }
    else if (glu_mass == 2875 ) { xsec = 0.115E-04; xsec_unc = 42.07 /100; return; }
    else if (glu_mass == 2880 ) { xsec = 0.112E-04; xsec_unc = 42.32 /100; return; }
    else if (glu_mass == 2885 ) { xsec = 0.109E-04; xsec_unc = 42.56 /100; return; }
    else if (glu_mass == 2890 ) { xsec = 0.106E-04; xsec_unc = 42.81 /100; return; }
    else if (glu_mass == 2895 ) { xsec = 0.104E-04; xsec_unc = 43.05 /100; return; }
    else if (glu_mass == 2900 ) { xsec = 0.101E-04; xsec_unc = 43.3 /100; return; }
    else if (glu_mass == 2905 ) { xsec = 0.986E-05; xsec_unc = 43.55 /100; return; }
    else if (glu_mass == 2910 ) { xsec = 0.961E-05; xsec_unc = 43.79 /100; return; }
    else if (glu_mass == 2915 ) { xsec = 0.937E-05; xsec_unc = 44.04 /100; return; }
    else if (glu_mass == 2920 ) { xsec = 0.914E-05; xsec_unc = 44.29 /100; return; }
    else if (glu_mass == 2925 ) { xsec = 0.891E-05; xsec_unc = 44.54 /100; return; }
    else if (glu_mass == 2930 ) { xsec = 0.869E-05; xsec_unc = 44.79 /100; return; }
    else if (glu_mass == 2935 ) { xsec = 0.848E-05; xsec_unc = 45.04 /100; return; }
    else if (glu_mass == 2940 ) { xsec = 0.827E-05; xsec_unc = 45.29 /100; return; }
    else if (glu_mass == 2945 ) { xsec = 0.806E-05; xsec_unc = 45.54 /100; return; }
    else if (glu_mass == 2950 ) { xsec = 0.786E-05; xsec_unc = 45.8 /100; return; }
    else if (glu_mass == 2955 ) { xsec = 0.767E-05; xsec_unc = 46.05 /100; return; }
    else if (glu_mass == 2960 ) { xsec = 0.748E-05; xsec_unc = 46.3 /100; return; }
    else if (glu_mass == 2965 ) { xsec = 0.729E-05; xsec_unc = 46.56 /100; return; }
    else if (glu_mass == 2970 ) { xsec = 0.711E-05; xsec_unc = 46.81 /100; return; }
    else if (glu_mass == 2975 ) { xsec = 0.694E-05; xsec_unc = 47.07 /100; return; }
    else if (glu_mass == 2980 ) { xsec = 0.677E-05; xsec_unc = 47.32 /100; return; }
    else if (glu_mass == 2985 ) { xsec = 0.660E-05; xsec_unc = 47.58 /100; return; }
    else if (glu_mass == 2990 ) { xsec = 0.644E-05; xsec_unc = 47.84 /100; return; }
    else if (glu_mass == 2995 ) { xsec = 0.628E-05; xsec_unc = 48.09 /100; return; }
    else if (glu_mass == 3000 ) { xsec = 0.612E-05; xsec_unc = 48.35 /100; return; }
    else {xsec = 0.; xsec_unc = 0.;} 
  }

  //void higgsino2DCrossSection(int hig_mass, double &xsec, double &xsec_unc) {
  //  if(hig_mass ==127) { xsec = .5824*.5824*1.44725; xsec_unc = 0.0395277; return;}
  //  else if(hig_mass ==150) { xsec = .5824*.5824*0.71514; xsec_unc = 0.0421496; return;}
  //  else if(hig_mass ==175) { xsec = .5824*.5824*0.419059; xsec_unc = 0.0453279; return;}
  //  else if(hig_mass ==200) { xsec = .5824*.5824*0.244213; xsec_unc = 0.047925; return;}
  //  else if(hig_mass ==225) { xsec = .5824*.5824*0.156286; xsec_unc = 0.0502876; return;}
  //  else if(hig_mass ==250) { xsec = .5824*.5824*0.104252; xsec_unc = 0.0526169; return;}
  //  else if(hig_mass ==275) { xsec = .5824*.5824*0.0719125; xsec_unc = 0.0549666; return;}
  //  else if(hig_mass ==300) { xsec = .5824*.5824*0.0509994; xsec_unc = 0.0572762; return;}
  //  else if(hig_mass ==325) { xsec = .5824*.5824*0.0369715; xsec_unc = 0.0590317; return;}
  //  else if(hig_mass ==350) { xsec = .5824*.5824*0.0273286; xsec_unc = 0.0607766; return;}
  //  else if(hig_mass ==375) { xsec = .5824*.5824*0.0205429; xsec_unc = 0.0625031; return;}
  //  else if(hig_mass ==400) { xsec = .5824*.5824*0.0156691; xsec_unc = 0.0642085; return;}
  //  else if(hig_mass ==425) { xsec = .5824*.5824*0.0120965; xsec_unc = 0.0657801; return;}
  //  else if(hig_mass ==450) { xsec = .5824*.5824*0.00944017; xsec_unc = 0.0674544; return;}
  //  else if(hig_mass ==475) { xsec = .5824*.5824*0.00743587; xsec_unc = 0.0686033; return;}
  //  else if(hig_mass ==500) { xsec = .5824*.5824*0.00590757; xsec_unc = 0.0699909; return;}
  //  else if(hig_mass ==525) { xsec = .5824*.5824*0.00469101; xsec_unc = 0.0713704; return;}
  //  else if(hig_mass ==550) { xsec = .5824*.5824*0.0038167; xsec_unc = 0.0722834; return;}
  //  else if(hig_mass ==575) { xsec = .5824*.5824*0.003073; xsec_unc = 0.0739957; return;}
  //  else if(hig_mass ==600) { xsec = .5824*.5824*0.00253015; xsec_unc = 0.0754291; return;}
  //  else if(hig_mass ==625) { xsec = .5824*.5824*0.00206136; xsec_unc = 0.0763466; return;}
  //  else if(hig_mass ==650) { xsec = .5824*.5824*0.00171418; xsec_unc = 0.0775695; return;}
  //  else if(hig_mass ==675) { xsec = .5824*.5824*0.00140934; xsec_unc = 0.0783375; return;}
  //  else if(hig_mass ==700) { xsec = .5824*.5824*0.00118113; xsec_unc = 0.0796388; return;}
  //  else if(hig_mass ==725) { xsec = .5824*.5824*0.000979349; xsec_unc = 0.0809883; return;}
  //  else if(hig_mass ==750) { xsec = .5824*.5824*0.000826366; xsec_unc = 0.081879; return;}
  //  else if(hig_mass ==775) { xsec = .5824*.5824*0.000690208; xsec_unc = 0.0842049; return;}
  //  else if(hig_mass ==800) { xsec = .5824*.5824*0.000586211; xsec_unc = 0.0862527; return;}
  //  else if(hig_mass ==825) { xsec = .5824*.5824*0.00049277; xsec_unc = 0.0864444; return;}
  //  else if(hig_mass ==850) { xsec = .5824*.5824*0.000420556; xsec_unc = 0.085742; return;}
  //  else if(hig_mass ==875) { xsec = .5824*.5824*0.000358734; xsec_unc = 0.0889174; return;}
  //  else if(hig_mass ==900) { xsec = .5824*.5824*0.000305935; xsec_unc = 0.0912439; return;}
  //  else if(hig_mass ==925) { xsec = .5824*.5824*0.000260948; xsec_unc = 0.091372; return;}
  //  else if(hig_mass ==950) { xsec = .5824*.5824*0.00022285; xsec_unc = 0.0919538; return;}
  //  else if(hig_mass ==975) { xsec = .5824*.5824*0.000189681; xsec_unc = 0.0938108; return;}
  //  else if(hig_mass ==1000) { xsec = .5824*.5824*0.00016428; xsec_unc = 0.0954285; return;}
  //  else if(hig_mass ==1025) { xsec = .5824*.5824*0.000142206; xsec_unc = 0.0957231; return;}
  //  else if(hig_mass ==1050) { xsec = .5824*.5824*0.000120971; xsec_unc = 0.0968997; return;}
  //  else if(hig_mass ==1075) { xsec = .5824*.5824*0.000105301; xsec_unc = 0.0979041; return;}
  //  else if(hig_mass ==1100) { xsec = .5824*.5824*9.12469e-05; xsec_unc = 0.0964142; return;}
  //  else if(hig_mass ==1125) { xsec = .5824*.5824*7.9765e-05; xsec_unc = 0.099902; return;}
  //  else if(hig_mass ==1150) { xsec = .5824*.5824*6.78234e-05; xsec_unc = 0.101061; return;}
  //  else if(hig_mass ==1175) { xsec = .5824*.5824*5.9016e-05; xsec_unc = 0.102051; return;}
  //  else if(hig_mass ==1200) { xsec = .5824*.5824*5.16263e-05; xsec_unc = 0.102499; return;}
  //  else if(hig_mass ==1225) { xsec = .5824*.5824*4.5147e-05; xsec_unc = 0.10403; return;}
  //  else if(hig_mass ==1250) { xsec = .5824*.5824*3.88343e-05; xsec_unc = 0.105206; return;}
  //  else if(hig_mass ==1275) { xsec = .5824*.5824*3.41304e-05; xsec_unc = 0.10619; return;}
  //  else if(hig_mass ==1300) { xsec = .5824*.5824*2.99353e-05; xsec_unc = 0.10783; return;}
  //  else if(hig_mass ==1325) { xsec = .5824*.5824*2.63637e-05; xsec_unc = 0.108024; return;}
  //  else if(hig_mass ==1350) { xsec = .5824*.5824*2.26779e-05; xsec_unc = 0.109016; return;}
  //  else if(hig_mass ==1375) { xsec = .5824*.5824*1.99318e-05; xsec_unc = 0.109822; return;}
  //  else if(hig_mass ==1400) { xsec = .5824*.5824*1.75031e-05; xsec_unc = 0.111631; return;}
  //  else if(hig_mass ==1425) { xsec = .5824*.5824*1.53974e-05; xsec_unc = 0.111417; return;}
  //  else if(hig_mass ==1450) { xsec = .5824*.5824*1.3245e-05; xsec_unc = 0.112313; return;}
  //  else if(hig_mass ==1475) { xsec = .5824*.5824*1.16416e-05; xsec_unc = 0.113058; return;}
  //  else{ xsec = 0; xsec_unc = 0;}
  //}

  // Unit: xsec=pb, xsec_unc=%
  void higgsino2DCrossSection(int hig_mass, double &xsec, double &xsec_unc) {
    if(hig_mass ==127) { xsec = .5824*.5824*1.44725; xsec_unc = 0.0395277; return;}
    else if(hig_mass ==150) { xsec = .5824*.5824*0.71514; xsec_unc = 0.0421496; return;}
    else if(hig_mass ==175) { xsec = .5824*.5824*0.419059; xsec_unc = 0.0453279; return;}
    else if(hig_mass ==200) { xsec = .5824*.5824*0.244213; xsec_unc = 0.047925; return;}
    else if(hig_mass ==225) { xsec = .5824*.5824*0.156286; xsec_unc = 0.0502876; return;}
    else if(hig_mass ==250) { xsec = .5824*.5824*0.104252; xsec_unc = 0.0526169; return;}
    else if(hig_mass ==275) { xsec = .5824*.5824*0.0719125; xsec_unc = 0.0549666; return;}
    else if(hig_mass ==300) { xsec = .5824*.5824*0.0509994; xsec_unc = 0.0572762; return;}
    else if(hig_mass ==325) { xsec = .5824*.5824*0.0369715; xsec_unc = 0.0590317; return;}
    else if(hig_mass ==350) { xsec = .5824*.5824*0.0273286; xsec_unc = 0.0607766; return;}
    else if(hig_mass ==375) { xsec = .5824*.5824*0.0205429; xsec_unc = 0.0625031; return;}
    else if(hig_mass ==400) { xsec = .5824*.5824*0.0156691; xsec_unc = 0.0642085; return;}
    else if(hig_mass ==425) { xsec = .5824*.5824*0.0120965; xsec_unc = 0.0657801; return;}
    else if(hig_mass ==450) { xsec = .5824*.5824*0.00944017; xsec_unc = 0.0674544; return;}
    else if(hig_mass ==475) { xsec = .5824*.5824*0.00743587; xsec_unc = 0.0686033; return;}
    else if(hig_mass ==500) { xsec = .5824*.5824*0.00590757; xsec_unc = 0.0699909; return;}
    else if(hig_mass ==525) { xsec = .5824*.5824*0.00473235; xsec_unc = 0.0713166; return;}
    else if(hig_mass ==550) { xsec = .5824*.5824*0.0038167; xsec_unc = 0.0722834; return;}
    else if(hig_mass ==575) { xsec = .5824*.5824*0.00309847; xsec_unc = 0.0739435; return;}
    else if(hig_mass ==600) { xsec = .5824*.5824*0.00253015; xsec_unc = 0.0754291; return;}
    else if(hig_mass ==625) { xsec = .5824*.5824*0.00207755; xsec_unc = 0.0763142; return;}
    else if(hig_mass ==650) { xsec = .5824*.5824*0.00171418; xsec_unc = 0.0775695; return;}
    else if(hig_mass ==675) { xsec = .5824*.5824*0.0014199; xsec_unc = 0.0782907; return;}
    else if(hig_mass ==700) { xsec = .5824*.5824*0.00118113; xsec_unc = 0.0796388; return;}
    else if(hig_mass ==725) { xsec = .5824*.5824*0.00098639; xsec_unc = 0.0809291; return;}
    else if(hig_mass ==750) { xsec = .5824*.5824*0.000826366; xsec_unc = 0.081879; return;}
    else if(hig_mass ==775) { xsec = .5824*.5824*0.000694985; xsec_unc = 0.0841352; return;}
    else if(hig_mass ==800) { xsec = .5824*.5824*0.000586211; xsec_unc = 0.0862527; return;}
    else if(hig_mass ==825) { xsec = .5824*.5824*0.000495914; xsec_unc = 0.0863944; return;}
    else if(hig_mass ==850) { xsec = .5824*.5824*0.000420556; xsec_unc = 0.085742; return;}
    else if(hig_mass ==875) { xsec = .5824*.5824*0.000361029; xsec_unc = 0.0888678; return;}
    else if(hig_mass ==900) { xsec = .5824*.5824*0.000305935; xsec_unc = 0.0912439; return;}
    else if(hig_mass ==925) { xsec = .5824*.5824*0.000262621; xsec_unc = 0.0913227; return;}
    else if(hig_mass ==950) { xsec = .5824*.5824*0.00022285; xsec_unc = 0.0919538; return;}
    else if(hig_mass ==975) { xsec = .5824*.5824*0.0001909; xsec_unc = 0.0937616; return;}
    else if(hig_mass ==1000) { xsec = .5824*.5824*0.00016428; xsec_unc = 0.0954285; return;}
    else if(hig_mass ==1025) { xsec = .5824*.5824*0.00014139; xsec_unc = 0.095765; return;}
    else if(hig_mass ==1050) { xsec = .5824*.5824*0.000121865; xsec_unc = 0.0967595; return;}
    else if(hig_mass ==1075) { xsec = .5824*.5824*0.000105913; xsec_unc = 0.0978621; return;}
    else if(hig_mass ==1100) { xsec = .5824*.5824*9.12469e-05; xsec_unc = 0.0964142; return;}
    else if(hig_mass ==1125) { xsec = .5824*.5824*7.93058e-05; xsec_unc = 0.0999433; return;}
    else if(hig_mass ==1150) { xsec = .5824*.5824*6.84561e-05; xsec_unc = 0.103594; return;}
    else if(hig_mass ==1175) { xsec = .5824*.5824*5.93602e-05; xsec_unc = 0.10201; return;}
    else if(hig_mass ==1200) { xsec = .5824*.5824*5.16263e-05; xsec_unc = 0.102499; return;}
    else if(hig_mass ==1225) { xsec = .5824*.5824*4.4906e-05; xsec_unc = 0.104071; return;}
    else if(hig_mass ==1250) { xsec = .5824*.5824*3.91587e-05; xsec_unc = 0.104736; return;}
    else if(hig_mass ==1275) { xsec = .5824*.5824*3.43135e-05; xsec_unc = 0.10615; return;}
    else if(hig_mass ==1300) { xsec = .5824*.5824*2.99353e-05; xsec_unc = 0.10783; return;}
    else if(hig_mass ==1325) { xsec = .5824*.5824*2.62223e-05; xsec_unc = 0.108061; return;}
    else if(hig_mass ==1350) { xsec = .5824*.5824*2.28072e-05; xsec_unc = 0.109427; return;}
    else if(hig_mass ==1375) { xsec = .5824*.5824*2.00393e-05; xsec_unc = 0.109789; return;}
    else if(hig_mass ==1400) { xsec = .5824*.5824*1.75031e-05; xsec_unc = 0.111631; return;}
    else if(hig_mass ==1425) { xsec = .5824*.5824*1.53144e-05; xsec_unc = 0.11145; return;}
    else if(hig_mass ==1450) { xsec = .5824*.5824*1.34572e-05; xsec_unc = 0.11084; return;}
    else if(hig_mass ==1475) { xsec = .5824*.5824*1.17047e-05; xsec_unc = 0.113027; return;}
    else{ xsec = 0; xsec_unc = 0;}
  }

  // Unit: xsec=pb, xsec_unc=%
  void higgsinoCrossSection(int hig_mass, double &xsec, double &xsec_unc) {
    if(hig_mass ==127) { xsec = .5824*.5824*7.6022; xsec_unc = 0.0393921; return;}
    else if(hig_mass ==150) { xsec = .5824*.5824*3.83231; xsec_unc = 0.0413612; return;}
    else if(hig_mass ==175) { xsec = .5824*.5824*2.26794; xsec_unc = 0.044299; return;}
    else if(hig_mass ==200) { xsec = .5824*.5824*1.33562; xsec_unc = 0.0474362; return;}
    else if(hig_mass ==225) { xsec = .5824*.5824*0.860597; xsec_unc = 0.0504217; return;}
    else if(hig_mass ==250) { xsec = .5824*.5824*0.577314; xsec_unc = 0.0532731; return;}
    else if(hig_mass ==275) { xsec = .5824*.5824*0.400107; xsec_unc = 0.0560232; return;}
    else if(hig_mass ==300) { xsec = .5824*.5824*0.284855; xsec_unc = 0.0586867; return;}
    else if(hig_mass ==325) { xsec = .5824*.5824*0.20736; xsec_unc = 0.0613554; return;}
    else if(hig_mass ==350) { xsec = .5824*.5824*0.153841; xsec_unc = 0.0640598; return;}
    else if(hig_mass ==375) { xsec = .5824*.5824*0.116006; xsec_unc = 0.066892; return;}
    else if(hig_mass ==400) { xsec = .5824*.5824*0.0887325; xsec_unc = 0.0697517; return;}
    else if(hig_mass ==425) { xsec = .5824*.5824*0.0686963; xsec_unc = 0.0723531; return;}
    else if(hig_mass ==450) { xsec = .5824*.5824*0.0537702; xsec_unc = 0.0748325; return;}
    else if(hig_mass ==475) { xsec = .5824*.5824*0.0424699; xsec_unc = 0.0775146; return;}
    else if(hig_mass ==500) { xsec = .5824*.5824*0.0338387; xsec_unc = 0.0802572; return;}
    else if(hig_mass ==525) { xsec = .5824*.5824*0.0271867; xsec_unc = 0.0825803; return;}
    else if(hig_mass ==550) { xsec = .5824*.5824*0.0219868; xsec_unc = 0.0849278; return;}
    else if(hig_mass ==575) { xsec = .5824*.5824*0.0179062; xsec_unc = 0.087561; return;}
    else if(hig_mass ==600) { xsec = .5824*.5824*0.0146677; xsec_unc = 0.0900693; return;}
    else if(hig_mass ==625) { xsec = .5824*.5824*0.012062; xsec_unc = 0.091959; return;}
    else if(hig_mass ==650) { xsec = .5824*.5824*0.00996406; xsec_unc = 0.094065; return;}
    else if(hig_mass ==675) { xsec = .5824*.5824*0.00828246; xsec_unc = 0.0957436; return;}
    else if(hig_mass ==700) { xsec = .5824*.5824*0.00689981; xsec_unc = 0.0982894; return;}
    else if(hig_mass ==725) { xsec = .5824*.5824*0.00578355; xsec_unc = 0.0999915; return;}
    else if(hig_mass ==750) { xsec = .5824*.5824*0.0048731; xsec_unc = 0.101211; return;}
    else if(hig_mass ==775) { xsec = .5824*.5824*0.00409781; xsec_unc = 0.104646; return;}
    else if(hig_mass ==800) { xsec = .5824*.5824*0.00346143; xsec_unc = 0.107618; return;}
    else if(hig_mass ==825) { xsec = .5824*.5824*0.0029337; xsec_unc = 0.108353; return;}
    else if(hig_mass ==850) { xsec = .5824*.5824*0.0024923; xsec_unc = 0.110016; return;}
    else if(hig_mass ==875) { xsec = .5824*.5824*0.00213679; xsec_unc = 0.112636; return;}
    else if(hig_mass ==900) { xsec = .5824*.5824*0.00180616; xsec_unc = 0.1134; return;}
    else if(hig_mass ==925) { xsec = .5824*.5824*0.00155453; xsec_unc = 0.116949; return;}
    else if(hig_mass ==950) { xsec = .5824*.5824*0.00132692; xsec_unc = 0.117027; return;}
    else if(hig_mass ==975) { xsec = .5824*.5824*0.00112975; xsec_unc = 0.121244; return;}
    else if(hig_mass ==1000) { xsec = .5824*.5824*0.000968853; xsec_unc = 0.126209; return;}
    else if(hig_mass ==1025) { xsec = .5824*.5824*0.000840602; xsec_unc = 0.121654; return;}
    else if(hig_mass ==1050) { xsec = .5824*.5824*0.000731306; xsec_unc = 0.118502; return;}
    else if(hig_mass ==1075) { xsec = .5824*.5824*0.000627083; xsec_unc = 0.127723; return;}
    else if(hig_mass ==1100) { xsec = .5824*.5824*0.000538005; xsec_unc = 0.134099; return;}
    else if(hig_mass ==1125) { xsec = .5824*.5824*0.00046747; xsec_unc = 0.133755; return;}
    else if(hig_mass ==1150) { xsec = .5824*.5824*0.000405108; xsec_unc = 0.120607; return;}
    else if(hig_mass ==1175) { xsec = .5824*.5824*0.000348261; xsec_unc = 0.139744; return;}
    else if(hig_mass ==1200) { xsec = .5824*.5824*0.000299347; xsec_unc = 0.162604; return;}
    else if(hig_mass ==1225) { xsec = .5824*.5824*0.000265935; xsec_unc = 0.137575; return;}
    else if(hig_mass ==1250) { xsec = .5824*.5824*0.000240471; xsec_unc = 0.119271; return;}
    else if(hig_mass ==1275) { xsec = .5824*.5824*0.000190411; xsec_unc = 0.138061; return;}
    else if(hig_mass ==1300) { xsec = .5824*.5824*0.000160765; xsec_unc = 0.122224; return;}
    else if(hig_mass ==1325) { xsec = .5824*.5824*0.000136272; xsec_unc = 0.138533; return;}
    else if(hig_mass ==1350) { xsec = .5824*.5824*0.000111174; xsec_unc = 0.177681; return;}
    else if(hig_mass ==1375) { xsec = .5824*.5824*9.74728e-05; xsec_unc = 0.138992; return;}
    else if(hig_mass ==1400) { xsec = .5824*.5824*7.80263e-05; xsec_unc = 0.118718; return;}
    else if(hig_mass ==1425) { xsec = .5824*.5824*6.96843e-05; xsec_unc = 0.139439; return;}
    else if(hig_mass ==1450) { xsec = .5824*.5824*6.96962e-05; xsec_unc = 0.198887; return;}
    else if(hig_mass ==1475) { xsec = .5824*.5824*4.98006e-05; xsec_unc = 0.139874; return;}
    else{ xsec = 0; xsec_unc = 0;}
  }
  
  void gluinoCrossSection(int glu_mass, double &xsec, double &xsec_unc) {
    //xsec = 0.5825*0.5825;
    //xsec_unc = 0.5825*0.5825;
    //temp fix for fullSIM, will make better later
    xsec = 1.0;
    xsec_unc = 1.0;
    if (glu_mass == 500) { xsec *= 33.8; xsec_unc *= 0.07869822485207102; return;}
    else if (glu_mass == 505) { xsec *= 31.9; xsec_unc *= 0.07899686520376176; return;}
    else if (glu_mass == 510) { xsec *= 30.1; xsec_unc *= 0.07906976744186046; return;}
    else if (glu_mass == 515) { xsec *= 28.4; xsec_unc *= 0.07957746478873239; return;}
    else if (glu_mass == 520) { xsec *= 26.8; xsec_unc *= 0.07985074626865672; return;}
    else if (glu_mass == 525) { xsec *= 25.3; xsec_unc *= 0.08023715415019762; return;}
    else if (glu_mass == 530) { xsec *= 24.0; xsec_unc *= 0.08083333333333333; return;}
    else if (glu_mass == 535) { xsec *= 22.7; xsec_unc *= 0.08105726872246696; return;}
    else if (glu_mass == 540) { xsec *= 21.4; xsec_unc *= 0.08130841121495327; return;}
    else if (glu_mass == 545) { xsec *= 20.3; xsec_unc *= 0.08177339901477831; return;}
    else if (glu_mass == 550) { xsec *= 19.2; xsec_unc *= 0.08229166666666668; return;}
    else if (glu_mass == 555) { xsec *= 18.2; xsec_unc *= 0.08241758241758242; return;}
    else if (glu_mass == 560) { xsec *= 17.2; xsec_unc *= 0.08313953488372093; return;}
    else if (glu_mass == 565) { xsec *= 16.3; xsec_unc *= 0.0834355828220859; return;}
    else if (glu_mass == 570) { xsec *= 15.5; xsec_unc *= 0.08387096774193549; return;}
    else if (glu_mass == 575) { xsec *= 14.7; xsec_unc *= 0.08435374149659865; return;}
    else if (glu_mass == 580) { xsec *= 13.9; xsec_unc *= 0.0841726618705036; return;}
    else if (glu_mass == 585) { xsec *= 13.2; xsec_unc *= 0.08484848484848487; return;}
    else if (glu_mass == 590) { xsec *= 12.6; xsec_unc *= 0.08492063492063492; return;}
    else if (glu_mass == 595) { xsec *= 11.9; xsec_unc *= 0.08571428571428572; return;}
    else if (glu_mass == 600) { xsec *= 11.3; xsec_unc *= 0.08619469026548672; return;}
    else if (glu_mass == 605) { xsec *= 10.8; xsec_unc *= 0.08657407407407407; return;}
    else if (glu_mass == 610) { xsec *= 10.2; xsec_unc *= 0.0869607843137255; return;}
    else if (glu_mass == 615) { xsec *= 9.74; xsec_unc *= 0.08737166324435318; return;}
    else if (glu_mass == 620) { xsec *= 9.26; xsec_unc *= 0.08779697624190064; return;}
    else if (glu_mass == 625) { xsec *= 8.81; xsec_unc *= 0.08819523269012486; return;}
    else if (glu_mass == 630) { xsec *= 8.39; xsec_unc *= 0.08855780691299164; return;}
    else if (glu_mass == 635) { xsec *= 7.99; xsec_unc *= 0.08898623279098873; return;}
    else if (glu_mass == 640) { xsec *= 7.61; xsec_unc *= 0.08935611038107753; return;}
    else if (glu_mass == 645) { xsec *= 7.25; xsec_unc *= 0.08979310344827586; return;}
    else if (glu_mass == 650) { xsec *= 6.9; xsec_unc *= 0.09014492753623188; return;}
    else if (glu_mass == 655) { xsec *= 6.58; xsec_unc *= 0.09057750759878419; return;}
    else if (glu_mass == 660) { xsec *= 6.27; xsec_unc *= 0.09106858054226476; return;}
    else if (glu_mass == 665) { xsec *= 5.98; xsec_unc *= 0.09147157190635452; return;}
    else if (glu_mass == 670) { xsec *= 5.71; xsec_unc *= 0.09194395796847636; return;}
    else if (glu_mass == 675) { xsec *= 5.44; xsec_unc *= 0.09227941176470587; return;}
    else if (glu_mass == 680) { xsec *= 5.2; xsec_unc *= 0.09269230769230768; return;}
    else if (glu_mass == 685) { xsec *= 4.96; xsec_unc *= 0.09314516129032259; return;}
    else if (glu_mass == 690) { xsec *= 4.74; xsec_unc *= 0.09345991561181434; return;}
    else if (glu_mass == 695) { xsec *= 4.52; xsec_unc *= 0.09380530973451329; return;}
    else if (glu_mass == 700) { xsec *= 4.32; xsec_unc *= 0.09421296296296296; return;}
    else if (glu_mass == 705) { xsec *= 4.13; xsec_unc *= 0.09467312348668282; return;}
    else if (glu_mass == 710) { xsec *= 3.95; xsec_unc *= 0.09493670886075949; return;}
    else if (glu_mass == 715) { xsec *= 3.77; xsec_unc *= 0.09549071618037135; return;}
    else if (glu_mass == 720) { xsec *= 3.61; xsec_unc *= 0.09584487534626038; return;}
    else if (glu_mass == 725) { xsec *= 3.45; xsec_unc *= 0.09623188405797102; return;}
    else if (glu_mass == 730) { xsec *= 3.3; xsec_unc *= 0.09636363636363637; return;}
    else if (glu_mass == 735) { xsec *= 3.16; xsec_unc *= 0.09683544303797467; return;}
    else if (glu_mass == 740) { xsec *= 3.02; xsec_unc *= 0.09735099337748344; return;}
    else if (glu_mass == 745) { xsec *= 2.89; xsec_unc *= 0.09757785467128026; return;}
    else if (glu_mass == 750) { xsec *= 2.77; xsec_unc *= 0.09783393501805054; return;}
    else if (glu_mass == 755) { xsec *= 2.65; xsec_unc *= 0.09811320754716982; return;}
    else if (glu_mass == 760) { xsec *= 2.54; xsec_unc *= 0.09881889763779528; return;}
    else if (glu_mass == 765) { xsec *= 2.43; xsec_unc *= 0.09917695473251027; return;}
    else if (glu_mass == 770) { xsec *= 2.33; xsec_unc *= 0.09957081545064378; return;}
    else if (glu_mass == 775) { xsec *= 2.23; xsec_unc *= 0.1; return;}
    else if (glu_mass == 780) { xsec *= 2.14; xsec_unc *= 0.09999999999999999; return;}
    else if (glu_mass == 785) { xsec *= 2.05; xsec_unc *= 0.10048780487804879; return;}
    else if (glu_mass == 790) { xsec *= 1.97; xsec_unc *= 0.10101522842639594; return;}
    else if (glu_mass == 795) { xsec *= 1.88; xsec_unc *= 0.10106382978723405; return;}
    else if (glu_mass == 800) { xsec *= 1.81; xsec_unc *= 0.1016574585635359; return;}
    else if (glu_mass == 805) { xsec *= 1.73; xsec_unc *= 0.10173410404624277; return;}
    else if (glu_mass == 810) { xsec *= 1.66; xsec_unc *= 0.10240963855421688; return;}
    else if (glu_mass == 815) { xsec *= 1.6; xsec_unc *= 0.1025; return;}
    else if (glu_mass == 820) { xsec *= 1.53; xsec_unc *= 0.10326797385620914; return;}
    else if (glu_mass == 825) { xsec *= 1.47; xsec_unc *= 0.10340136054421768; return;}
    else if (glu_mass == 830) { xsec *= 1.41; xsec_unc *= 0.10354609929078014; return;}
    else if (glu_mass == 835) { xsec *= 1.36; xsec_unc *= 0.10441176470588234; return;}
    else if (glu_mass == 840) { xsec *= 1.3; xsec_unc *= 0.10461538461538462; return;}
    else if (glu_mass == 845) { xsec *= 1.25; xsec_unc *= 0.1048; return;}
    else if (glu_mass == 850) { xsec *= 1.2; xsec_unc *= 0.10500000000000001; return;}
    else if (glu_mass == 855) { xsec *= 1.15; xsec_unc *= 0.10608695652173913; return;}
    else if (glu_mass == 860) { xsec *= 1.11; xsec_unc *= 0.10630630630630629; return;}
    else if (glu_mass == 865) { xsec *= 1.07; xsec_unc *= 0.10654205607476636; return;}
    else if (glu_mass == 870) { xsec *= 1.03; xsec_unc *= 0.10679611650485436; return;}
    else if (glu_mass == 875) { xsec *= 0.986; xsec_unc *= 0.1075050709939148; return;}
    else if (glu_mass == 880) { xsec *= 0.948; xsec_unc *= 0.10759493670886076; return;}
    else if (glu_mass == 885) { xsec *= 0.912; xsec_unc *= 0.10789473684210527; return;}
    else if (glu_mass == 890) { xsec *= 0.877; xsec_unc *= 0.10820980615735462; return;}
    else if (glu_mass == 895) { xsec *= 0.844; xsec_unc *= 0.10864928909952608; return;}
    else if (glu_mass == 900) { xsec *= 0.812; xsec_unc *= 0.10886699507389162; return;}
    else if (glu_mass == 905) { xsec *= 0.781; xsec_unc *= 0.10934699103713189; return;}
    else if (glu_mass == 910) { xsec *= 0.752; xsec_unc *= 0.10970744680851065; return;}
    else if (glu_mass == 915) { xsec *= 0.723; xsec_unc *= 0.1099585062240664; return;}
    else if (glu_mass == 920) { xsec *= 0.696; xsec_unc *= 0.1103448275862069; return;}
    else if (glu_mass == 925) { xsec *= 0.67; xsec_unc *= 0.11074626865671641; return;}
    else if (glu_mass == 930) { xsec *= 0.646; xsec_unc *= 0.11114551083591331; return;}
    else if (glu_mass == 935) { xsec *= 0.622; xsec_unc *= 0.11141479099678457; return;}
    else if (glu_mass == 940) { xsec *= 0.599; xsec_unc *= 0.11185308848080135; return;}
    else if (glu_mass == 945) { xsec *= 0.577; xsec_unc *= 0.1121317157712305; return;}
    else if (glu_mass == 950) { xsec *= 0.556; xsec_unc *= 0.11258992805755395; return;}
    else if (glu_mass == 955) { xsec *= 0.535; xsec_unc *= 0.11271028037383177; return;}
    else if (glu_mass == 960) { xsec *= 0.516; xsec_unc *= 0.11317829457364341; return;}
    else if (glu_mass == 965) { xsec *= 0.497; xsec_unc *= 0.11348088531187123; return;}
    else if (glu_mass == 970) { xsec *= 0.479; xsec_unc *= 0.1139874739039666; return;}
    else if (glu_mass == 975) { xsec *= 0.462; xsec_unc *= 0.11428571428571428; return;}
    else if (glu_mass == 980) { xsec *= 0.445; xsec_unc *= 0.1146067415730337; return;}
    else if (glu_mass == 985) { xsec *= 0.43; xsec_unc *= 0.11488372093023255; return;}
    else if (glu_mass == 990) { xsec *= 0.414; xsec_unc *= 0.11521739130434783; return;}
    else if (glu_mass == 995) { xsec *= 0.399; xsec_unc *= 0.11553884711779448; return;}
    else if (glu_mass == 1000) { xsec *= 0.385; xsec_unc *= 0.1161038961038961; return;}
    else if (glu_mass == 1005) { xsec *= 0.372; xsec_unc *= 0.11639784946236559; return;}
    else if (glu_mass == 1010) { xsec *= 0.359; xsec_unc *= 0.11671309192200557; return;}
    else if (glu_mass == 1015) { xsec *= 0.346; xsec_unc *= 0.1170520231213873; return;}
    else if (glu_mass == 1020) { xsec *= 0.334; xsec_unc *= 0.11736526946107784; return;}
    else if (glu_mass == 1025) { xsec *= 0.322; xsec_unc *= 0.11770186335403728; return;}
    else if (glu_mass == 1030) { xsec *= 0.311; xsec_unc *= 0.11800643086816721; return;}
    else if (glu_mass == 1035) { xsec *= 0.3; xsec_unc *= 0.11866666666666667; return;}
    else if (glu_mass == 1040) { xsec *= 0.29; xsec_unc *= 0.11896551724137933; return;}
    else if (glu_mass == 1045) { xsec *= 0.28; xsec_unc *= 0.11928571428571427; return;}
    else if (glu_mass == 1050) { xsec *= 0.27; xsec_unc *= 0.11962962962962963; return;}
    else if (glu_mass == 1055) { xsec *= 0.261; xsec_unc *= 0.11992337164750957; return;}
    else if (glu_mass == 1060) { xsec *= 0.252; xsec_unc *= 0.12023809523809524; return;}
    else if (glu_mass == 1065) { xsec *= 0.243; xsec_unc *= 0.1205761316872428; return;}
    else if (glu_mass == 1070) { xsec *= 0.235; xsec_unc *= 0.12085106382978725; return;}
    else if (glu_mass == 1075) { xsec *= 0.227; xsec_unc *= 0.1211453744493392; return;}
    else if (glu_mass == 1080) { xsec *= 0.219; xsec_unc *= 0.1219178082191781; return;}
    else if (glu_mass == 1085) { xsec *= 0.212; xsec_unc *= 0.12216981132075472; return;}
    else if (glu_mass == 1090) { xsec *= 0.205; xsec_unc *= 0.12243902439024391; return;}
    else if (glu_mass == 1095) { xsec *= 0.198; xsec_unc *= 0.12272727272727271; return;}
    else if (glu_mass == 1100) { xsec *= 0.191; xsec_unc *= 0.12303664921465969; return;}
    else if (glu_mass == 1105) { xsec *= 0.185; xsec_unc *= 0.12324324324324325; return;}
    else if (glu_mass == 1110) { xsec *= 0.179; xsec_unc *= 0.12402234636871509; return;}
    else if (glu_mass == 1115) { xsec *= 0.173; xsec_unc *= 0.12427745664739884; return;}
    else if (glu_mass == 1120) { xsec *= 0.167; xsec_unc *= 0.1245508982035928; return;}
    else if (glu_mass == 1125) { xsec *= 0.162; xsec_unc *= 0.12469135802469135; return;}
    else if (glu_mass == 1130) { xsec *= 0.156; xsec_unc *= 0.125; return;}
    else if (glu_mass == 1135) { xsec *= 0.151; xsec_unc *= 0.12582781456953643; return;}
    else if (glu_mass == 1140) { xsec *= 0.146; xsec_unc *= 0.12602739726027398; return;}
    else if (glu_mass == 1145) { xsec *= 0.141; xsec_unc *= 0.12624113475177307; return;}
    else if (glu_mass == 1150) { xsec *= 0.137; xsec_unc *= 0.12700729927007298; return;}
    else if (glu_mass == 1155) { xsec *= 0.132; xsec_unc *= 0.12727272727272726; return;}
    else if (glu_mass == 1160) { xsec *= 0.128; xsec_unc *= 0.12734374999999998; return;}
    else if (glu_mass == 1165) { xsec *= 0.124; xsec_unc *= 0.1274193548387097; return;}
    else if (glu_mass == 1170) { xsec *= 0.12; xsec_unc *= 0.12833333333333335; return;}
    else if (glu_mass == 1175) { xsec *= 0.116; xsec_unc *= 0.12844827586206897; return;}
    else if (glu_mass == 1180) { xsec *= 0.112; xsec_unc *= 0.12857142857142856; return;}
    else if (glu_mass == 1185) { xsec *= 0.109; xsec_unc *= 0.12935779816513762; return;}
    else if (glu_mass == 1190) { xsec *= 0.105; xsec_unc *= 0.1295238095238095; return;}
    else if (glu_mass == 1195) { xsec *= 0.102; xsec_unc *= 0.1303921568627451; return;}
    else if (glu_mass == 1200) { xsec *= 0.0985; xsec_unc *= 0.1299492385786802; return;}
    else if (glu_mass == 1205) { xsec *= 0.0953; xsec_unc *= 0.1311647429171039; return;}
    else if (glu_mass == 1210) { xsec *= 0.0923; xsec_unc *= 0.13109425785482123; return;}
    else if (glu_mass == 1215) { xsec *= 0.0894; xsec_unc *= 0.13087248322147652; return;}
    else if (glu_mass == 1220) { xsec *= 0.0866; xsec_unc *= 0.13163972286374134; return;}
    else if (glu_mass == 1225) { xsec *= 0.0838; xsec_unc *= 0.1324582338902148; return;}
    else if (glu_mass == 1230) { xsec *= 0.0812; xsec_unc *= 0.1330049261083744; return;}
    else if (glu_mass == 1235) { xsec *= 0.0786; xsec_unc *= 0.13231552162849872; return;}
    else if (glu_mass == 1240) { xsec *= 0.0762; xsec_unc *= 0.13254593175853016; return;}
    else if (glu_mass == 1245) { xsec *= 0.0738; xsec_unc *= 0.13333333333333333; return;}
    else if (glu_mass == 1250) { xsec *= 0.0715; xsec_unc *= 0.13384615384615386; return;}
    else if (glu_mass == 1255) { xsec *= 0.0692; xsec_unc *= 0.13410404624277458; return;}
    else if (glu_mass == 1260) { xsec *= 0.0671; xsec_unc *= 0.13442622950819672; return;}
    else if (glu_mass == 1265) { xsec *= 0.065; xsec_unc *= 0.13476923076923078; return;}
    else if (glu_mass == 1270) { xsec *= 0.063; xsec_unc *= 0.13507936507936508; return;}
    else if (glu_mass == 1275) { xsec *= 0.061; xsec_unc *= 0.13557377049180327; return;}
    else if (glu_mass == 1280) { xsec *= 0.0591; xsec_unc *= 0.1358714043993232; return;}
    else if (glu_mass == 1285) { xsec *= 0.0573; xsec_unc *= 0.13612565445026178; return;}
    else if (glu_mass == 1290) { xsec *= 0.0556; xsec_unc *= 0.1365107913669065; return;}
    else if (glu_mass == 1295) { xsec *= 0.0539; xsec_unc *= 0.13692022263450834; return;}
    else if (glu_mass == 1300) { xsec *= 0.0522; xsec_unc *= 0.13716475095785438; return;}
    else if (glu_mass == 1305) { xsec *= 0.0506; xsec_unc *= 0.13754940711462452; return;}
    else if (glu_mass == 1310) { xsec *= 0.0491; xsec_unc *= 0.13788187372708757; return;}
    else if (glu_mass == 1315) { xsec *= 0.0476; xsec_unc *= 0.13823529411764704; return;}
    else if (glu_mass == 1320) { xsec *= 0.0461; xsec_unc *= 0.1386117136659436; return;}
    else if (glu_mass == 1325) { xsec *= 0.0447; xsec_unc *= 0.13892617449664432; return;}
    else if (glu_mass == 1330) { xsec *= 0.0434; xsec_unc *= 0.13940092165898615; return;}
    else if (glu_mass == 1335) { xsec *= 0.0421; xsec_unc *= 0.13966745843230405; return;}
    else if (glu_mass == 1340) { xsec *= 0.0408; xsec_unc *= 0.14019607843137255; return;}
    else if (glu_mass == 1345) { xsec *= 0.0396; xsec_unc *= 0.1404040404040404; return;}
    else if (glu_mass == 1350) { xsec *= 0.0384; xsec_unc *= 0.14088541666666668; return;}
    else if (glu_mass == 1355) { xsec *= 0.0372; xsec_unc *= 0.14112903225806453; return;}
    else if (glu_mass == 1360) { xsec *= 0.0361; xsec_unc *= 0.14155124653739612; return;}
    else if (glu_mass == 1365) { xsec *= 0.035; xsec_unc *= 0.142; return;}
    else if (glu_mass == 1370) { xsec *= 0.034; xsec_unc *= 0.14205882352941177; return;}
    else if (glu_mass == 1375) { xsec *= 0.033; xsec_unc *= 0.1427272727272727; return;}
    else if (glu_mass == 1380) { xsec *= 0.032; xsec_unc *= 0.143125; return;}
    else if (glu_mass == 1385) { xsec *= 0.031; xsec_unc *= 0.14322580645161292; return;}
    else if (glu_mass == 1390) { xsec *= 0.0301; xsec_unc *= 0.14385382059800664; return;}
    else if (glu_mass == 1395) { xsec *= 0.0292; xsec_unc *= 0.14383561643835616; return;}
    else if (glu_mass == 1400) { xsec *= 0.0284; xsec_unc *= 0.1443661971830986; return;}
    else if (glu_mass == 1405) { xsec *= 0.0275; xsec_unc *= 0.14472727272727273; return;}
    else if (glu_mass == 1410) { xsec *= 0.0267; xsec_unc *= 0.1449438202247191; return;}
    else if (glu_mass == 1415) { xsec *= 0.0259; xsec_unc *= 0.14555984555984555; return;}
    else if (glu_mass == 1420) { xsec *= 0.0252; xsec_unc *= 0.14603174603174604; return;}
    else if (glu_mass == 1425) { xsec *= 0.0244; xsec_unc *= 0.14631147540983605; return;}
    else if (glu_mass == 1430) { xsec *= 0.0237; xsec_unc *= 0.14641350210970464; return;}
    else if (glu_mass == 1435) { xsec *= 0.023; xsec_unc *= 0.14695652173913046; return;}
    else if (glu_mass == 1440) { xsec *= 0.0224; xsec_unc *= 0.14732142857142858; return;}
    else if (glu_mass == 1445) { xsec *= 0.0217; xsec_unc *= 0.147926267281106; return;}
    else if (glu_mass == 1450) { xsec *= 0.0211; xsec_unc *= 0.14786729857819905; return;}
    else if (glu_mass == 1455) { xsec *= 0.0205; xsec_unc *= 0.14829268292682926; return;}
    else if (glu_mass == 1460) { xsec *= 0.0199; xsec_unc *= 0.1487437185929648; return;}
    else if (glu_mass == 1465) { xsec *= 0.0193; xsec_unc *= 0.14922279792746115; return;}
    else if (glu_mass == 1470) { xsec *= 0.0187; xsec_unc *= 0.1497326203208556; return;}
    else if (glu_mass == 1475) { xsec *= 0.0182; xsec_unc *= 0.15; return;}
    else if (glu_mass == 1480) { xsec *= 0.0177; xsec_unc *= 0.15028248587570622; return;}
    else if (glu_mass == 1485) { xsec *= 0.0172; xsec_unc *= 0.1505813953488372; return;}
    else if (glu_mass == 1490) { xsec *= 0.0167; xsec_unc *= 0.1508982035928144; return;}
    else if (glu_mass == 1495) { xsec *= 0.0162; xsec_unc *= 0.15123456790123457; return;}
    else if (glu_mass == 1500) { xsec *= 0.0157; xsec_unc *= 0.1515923566878981; return;}
    else if (glu_mass == 1505) { xsec *= 0.0153; xsec_unc *= 0.15228758169934642; return;}
    else if (glu_mass == 1510) { xsec *= 0.0148; xsec_unc *= 0.1527027027027027; return;}
    else if (glu_mass == 1515) { xsec *= 0.0144; xsec_unc *= 0.1527777777777778; return;}
    else if (glu_mass == 1520) { xsec *= 0.014; xsec_unc *= 0.15357142857142858; return;}
    else if (glu_mass == 1525) { xsec *= 0.0136; xsec_unc *= 0.1536764705882353; return;}
    else if (glu_mass == 1530) { xsec *= 0.0132; xsec_unc *= 0.1537878787878788; return;}
    else if (glu_mass == 1535) { xsec *= 0.0128; xsec_unc *= 0.1546875; return;}
    else if (glu_mass == 1540) { xsec *= 0.0125; xsec_unc *= 0.1552; return;}
    else if (glu_mass == 1545) { xsec *= 0.0121; xsec_unc *= 0.15537190082644628; return;}
    else if (glu_mass == 1550) { xsec *= 0.0118; xsec_unc *= 0.15593220338983052; return;}
    else if (glu_mass == 1555) { xsec *= 0.0115; xsec_unc *= 0.15565217391304348; return;}
    else if (glu_mass == 1560) { xsec *= 0.0111; xsec_unc *= 0.15675675675675674; return;}
    else if (glu_mass == 1565) { xsec *= 0.0108; xsec_unc *= 0.15648148148148147; return;}
    else if (glu_mass == 1570) { xsec *= 0.0105; xsec_unc *= 0.15714285714285714; return;}
    else if (glu_mass == 1575) { xsec *= 0.0102; xsec_unc *= 0.15784313725490196; return;}
    else if (glu_mass == 1580) { xsec *= 0.00993; xsec_unc *= 0.1581067472306143; return;}
    else if (glu_mass == 1585) { xsec *= 0.00966; xsec_unc *= 0.15838509316770186; return;}
    else if (glu_mass == 1590) { xsec *= 0.00939; xsec_unc *= 0.15867944621938232; return;}
    else if (glu_mass == 1595) { xsec *= 0.00912; xsec_unc *= 0.15899122807017543; return;}
    else if (glu_mass == 1600) { xsec *= 0.00887; xsec_unc *= 0.15896279594137544; return;}
    else if (glu_mass == 1605) { xsec *= 0.00862; xsec_unc *= 0.16009280742459397; return;}
    else if (glu_mass == 1610) { xsec *= 0.00838; xsec_unc *= 0.15990453460620524; return;}
    else if (glu_mass == 1615) { xsec *= 0.00815; xsec_unc *= 0.16073619631901842; return;}
    else if (glu_mass == 1620) { xsec *= 0.00792; xsec_unc *= 0.16161616161616163; return;}
    else if (glu_mass == 1625) { xsec *= 0.0077; xsec_unc *= 0.16103896103896104; return;}
    else if (glu_mass == 1630) { xsec *= 0.00749; xsec_unc *= 0.16154873164218958; return;}
    else if (glu_mass == 1635) { xsec *= 0.00728; xsec_unc *= 0.1620879120879121; return;}
    else if (glu_mass == 1640) { xsec *= 0.00708; xsec_unc *= 0.16242937853107345; return;}
    else if (glu_mass == 1645) { xsec *= 0.00689; xsec_unc *= 0.16255442670537007; return;}
    else if (glu_mass == 1650) { xsec *= 0.0067; xsec_unc *= 0.1626865671641791; return;}
    else if (glu_mass == 1655) { xsec *= 0.00651; xsec_unc *= 0.16436251920122888; return;}
    else if (glu_mass == 1660) { xsec *= 0.00633; xsec_unc *= 0.16429699842022116; return;}
    else if (glu_mass == 1665) { xsec *= 0.00616; xsec_unc *= 0.16396103896103897; return;}
    else if (glu_mass == 1670) { xsec *= 0.00599; xsec_unc *= 0.1649415692821369; return;}
    else if (glu_mass == 1675) { xsec *= 0.00583; xsec_unc *= 0.1653516295025729; return;}
    else if (glu_mass == 1680) { xsec *= 0.00567; xsec_unc *= 0.16560846560846562; return;}
    else if (glu_mass == 1685) { xsec *= 0.00551; xsec_unc *= 0.16606170598911071; return;}
    else if (glu_mass == 1690) { xsec *= 0.00536; xsec_unc *= 0.1664179104477612; return;}
    else if (glu_mass == 1695) { xsec *= 0.00521; xsec_unc *= 0.16679462571976966; return;}
    else if (glu_mass == 1700) { xsec *= 0.00507; xsec_unc *= 0.16725838264299803; return;}
    else if (glu_mass == 1705) { xsec *= 0.00493; xsec_unc *= 0.16754563894523325; return;}
    else if (glu_mass == 1710) { xsec *= 0.0048; xsec_unc *= 0.16812500000000002; return;}
    else if (glu_mass == 1715) { xsec *= 0.00467; xsec_unc *= 0.16852248394004285; return;}
    else if (glu_mass == 1720) { xsec *= 0.00454; xsec_unc *= 0.16894273127753304; return;}
    else if (glu_mass == 1725) { xsec *= 0.00442; xsec_unc *= 0.1692307692307692; return;}
    else if (glu_mass == 1730) { xsec *= 0.0043; xsec_unc *= 0.1697674418604651; return;}
    else if (glu_mass == 1735) { xsec *= 0.00418; xsec_unc *= 0.17009569377990433; return;}
    else if (glu_mass == 1740) { xsec *= 0.00407; xsec_unc *= 0.1705159705159705; return;}
    else if (glu_mass == 1745) { xsec *= 0.00396; xsec_unc *= 0.17095959595959595; return;}
    else if (glu_mass == 1750) { xsec *= 0.00385; xsec_unc *= 0.17142857142857143; return;}
    else if (glu_mass == 1755) { xsec *= 0.00375; xsec_unc *= 0.17173333333333335; return;}
    else if (glu_mass == 1760) { xsec *= 0.00365; xsec_unc *= 0.17232876712328768; return;}
    else if (glu_mass == 1765) { xsec *= 0.00355; xsec_unc *= 0.17267605633802818; return;}
    else if (glu_mass == 1770) { xsec *= 0.00345; xsec_unc *= 0.17304347826086958; return;}
    else if (glu_mass == 1775) { xsec *= 0.00336; xsec_unc *= 0.17351190476190476; return;}
    else if (glu_mass == 1780) { xsec *= 0.00327; xsec_unc *= 0.17400611620795106; return;}
    else if (glu_mass == 1785) { xsec *= 0.00318; xsec_unc *= 0.1742138364779874; return;}
    else if (glu_mass == 1790) { xsec *= 0.0031; xsec_unc *= 0.17483870967741935; return;}
    else if (glu_mass == 1795) { xsec *= 0.00301; xsec_unc *= 0.1750830564784053; return;}
    else if (glu_mass == 1800) { xsec *= 0.00293; xsec_unc *= 0.17576791808873723; return;}
    else if (glu_mass == 1805) { xsec *= 0.00286; xsec_unc *= 0.17587412587412585; return;}
    else if (glu_mass == 1810) { xsec *= 0.00278; xsec_unc *= 0.17625899280575538; return;}
    else if (glu_mass == 1815) { xsec *= 0.00271; xsec_unc *= 0.17675276752767527; return;}
    else if (glu_mass == 1820) { xsec *= 0.00263; xsec_unc *= 0.17718631178707225; return;}
    else if (glu_mass == 1825) { xsec *= 0.00256; xsec_unc *= 0.17773437499999997; return;}
    else if (glu_mass == 1830) { xsec *= 0.00249; xsec_unc *= 0.1783132530120482; return;}
    else if (glu_mass == 1835) { xsec *= 0.00243; xsec_unc *= 0.1786008230452675; return;}
    else if (glu_mass == 1840) { xsec *= 0.00236; xsec_unc *= 0.17881355932203388; return;}
    else if (glu_mass == 1845) { xsec *= 0.0023; xsec_unc *= 0.17956521739130435; return;}
    else if (glu_mass == 1850) { xsec *= 0.00224; xsec_unc *= 0.17991071428571428; return;}
    else if (glu_mass == 1855) { xsec *= 0.00218; xsec_unc *= 0.18027522935779816; return;}
    else if (glu_mass == 1860) { xsec *= 0.00212; xsec_unc *= 0.18066037735849055; return;}
    else if (glu_mass == 1865) { xsec *= 0.00207; xsec_unc *= 0.1811594202898551; return;}
    else if (glu_mass == 1870) { xsec *= 0.00201; xsec_unc *= 0.181592039800995; return;}
    else if (glu_mass == 1875) { xsec *= 0.00196; xsec_unc *= 0.18214285714285716; return;}
    else if (glu_mass == 1880) { xsec *= 0.00191; xsec_unc *= 0.18272251308900525; return;}
    else if (glu_mass == 1885) { xsec *= 0.00186; xsec_unc *= 0.18333333333333332; return;}
    else if (glu_mass == 1890) { xsec *= 0.00181; xsec_unc *= 0.18342541436464088; return;}
    else if (glu_mass == 1895) { xsec *= 0.00176; xsec_unc *= 0.18409090909090908; return;}
    else if (glu_mass == 1900) { xsec *= 0.00171; xsec_unc *= 0.1842105263157895; return;}
    else if (glu_mass == 1905) { xsec *= 0.00167; xsec_unc *= 0.18502994011976046; return;}
    else if (glu_mass == 1910) { xsec *= 0.00163; xsec_unc *= 0.18527607361963191; return;}
    else if (glu_mass == 1915) { xsec *= 0.00158; xsec_unc *= 0.1860759493670886; return;}
    else if (glu_mass == 1920) { xsec *= 0.00154; xsec_unc *= 0.18636363636363637; return;}
    else if (glu_mass == 1925) { xsec *= 0.0015; xsec_unc *= 0.18666666666666665; return;}
    else if (glu_mass == 1930) { xsec *= 0.00146; xsec_unc *= 0.18698630136986305; return;}
    else if (glu_mass == 1935) { xsec *= 0.00142; xsec_unc *= 0.18802816901408448; return;}
    else if (glu_mass == 1940) { xsec *= 0.00139; xsec_unc *= 0.18848920863309354; return;}
    else if (glu_mass == 1945) { xsec *= 0.00135; xsec_unc *= 0.18888888888888888; return;}
    else if (glu_mass == 1950) { xsec *= 0.00131; xsec_unc *= 0.18931297709923667; return;}
    else if (glu_mass == 1955) { xsec *= 0.00128; xsec_unc *= 0.18984374999999998; return;}
    else if (glu_mass == 1960) { xsec *= 0.00125; xsec_unc *= 0.1904; return;}
    else if (glu_mass == 1965) { xsec *= 0.00121; xsec_unc *= 0.19090909090909092; return;}
    else if (glu_mass == 1970) { xsec *= 0.00118; xsec_unc *= 0.19152542372881354; return;}
    else if (glu_mass == 1975) { xsec *= 0.00115; xsec_unc *= 0.19130434782608696; return;}
    else if (glu_mass == 1980) { xsec *= 0.00112; xsec_unc *= 0.19196428571428573; return;}
    else if (glu_mass == 1985) { xsec *= 0.00109; xsec_unc *= 0.1926605504587156; return;}
    else if (glu_mass == 1990) { xsec *= 0.00106; xsec_unc *= 0.19339622641509435; return;}
    else if (glu_mass == 1995) { xsec *= 0.00104; xsec_unc *= 0.1932692307692308; return;}
    else if (glu_mass == 2000) { xsec *= 0.00101; xsec_unc *= 0.19405940594059404; return;}
    else if (glu_mass == 2005) { xsec *= 0.000983; xsec_unc *= 0.19430315361139372; return;}
    else if (glu_mass == 2010) { xsec *= 0.000957; xsec_unc *= 0.19540229885057472; return;}
    else if (glu_mass == 2015) { xsec *= 0.000933; xsec_unc *= 0.19614147909967847; return;}
    else if (glu_mass == 2020) { xsec *= 0.000908; xsec_unc *= 0.1960352422907489; return;}
    else if (glu_mass == 2025) { xsec *= 0.000885; xsec_unc *= 0.1966101694915254; return;}
    else if (glu_mass == 2030) { xsec *= 0.000862; xsec_unc *= 0.19721577726218098; return;}
    else if (glu_mass == 2035) { xsec *= 0.00084; xsec_unc *= 0.1976190476190476; return;}
    else if (glu_mass == 2040) { xsec *= 0.000818; xsec_unc *= 0.1980440097799511; return;}
    else if (glu_mass == 2045) { xsec *= 0.000797; xsec_unc *= 0.19949811794228356; return;}
    else if (glu_mass == 2050) { xsec *= 0.000776; xsec_unc *= 0.19974226804123713; return;}
    else if (glu_mass == 2055) { xsec *= 0.000756; xsec_unc *= 0.19973544973544974; return;}
    else if (glu_mass == 2060) { xsec *= 0.000737; xsec_unc *= 0.20081411126187243; return;}
    else if (glu_mass == 2065) { xsec *= 0.000718; xsec_unc *= 0.201949860724234; return;}
    else if (glu_mass == 2070) { xsec *= 0.000699; xsec_unc *= 0.20171673819742492; return;}
    else if (glu_mass == 2075) { xsec *= 0.000681; xsec_unc *= 0.2026431718061674; return;}
    else if (glu_mass == 2080) { xsec *= 0.000664; xsec_unc *= 0.2033132530120482; return;}
    else if (glu_mass == 2085) { xsec *= 0.000647; xsec_unc *= 0.20401854714064915; return;}
    else if (glu_mass == 2090) { xsec *= 0.00063; xsec_unc *= 0.20476190476190473; return;}
    else if (glu_mass == 2095) { xsec *= 0.000614; xsec_unc *= 0.20521172638436483; return;}
    else if (glu_mass == 2100) { xsec *= 0.000598; xsec_unc *= 0.205685618729097; return;}
    else if (glu_mass == 2105) { xsec *= 0.000583; xsec_unc *= 0.2058319039451115; return;}
    else if (glu_mass == 2110) { xsec *= 0.000568; xsec_unc *= 0.20598591549295772; return;}
    else if (glu_mass == 2115) { xsec *= 0.000553; xsec_unc *= 0.20795660036166366; return;}
    else if (glu_mass == 2120) { xsec *= 0.000539; xsec_unc *= 0.2077922077922078; return;}
    else if (glu_mass == 2125) { xsec *= 0.000525; xsec_unc *= 0.20761904761904765; return;}
    else if (glu_mass == 2130) { xsec *= 0.000512; xsec_unc *= 0.208984375; return;}
    else if (glu_mass == 2135) { xsec *= 0.000499; xsec_unc *= 0.21042084168336675; return;}
    else if (glu_mass == 2140) { xsec *= 0.000486; xsec_unc *= 0.20987654320987653; return;}
    else if (glu_mass == 2145) { xsec *= 0.000473; xsec_unc *= 0.2109936575052854; return;}
    else if (glu_mass == 2150) { xsec *= 0.000461; xsec_unc *= 0.21149674620390455; return;}
    else if (glu_mass == 2155) { xsec *= 0.000449; xsec_unc *= 0.21224944320712694; return;}
    else if (glu_mass == 2160) { xsec *= 0.000438; xsec_unc *= 0.213013698630137; return;}
    else if (glu_mass == 2165) { xsec *= 0.000427; xsec_unc *= 0.21358313817330207; return;}
    else if (glu_mass == 2170) { xsec *= 0.000416; xsec_unc *= 0.21418269230769232; return;}
    else if (glu_mass == 2175) { xsec *= 0.000405; xsec_unc *= 0.21481481481481482; return;}
    else if (glu_mass == 2180) { xsec *= 0.000395; xsec_unc *= 0.21544303797468353; return;}
    else if (glu_mass == 2185) { xsec *= 0.000385; xsec_unc *= 0.21610389610389613; return;}
    else if (glu_mass == 2190) { xsec *= 0.000375; xsec_unc *= 0.2168; return;}
    else if (glu_mass == 2195) { xsec *= 0.000365; xsec_unc *= 0.21753424657534248; return;}
    else if (glu_mass == 2200) { xsec *= 0.000356; xsec_unc *= 0.21825842696629216; return;}
    else if (glu_mass == 2205) { xsec *= 0.000347; xsec_unc *= 0.21902017291066286; return;}
    else if (glu_mass == 2210) { xsec *= 0.000338; xsec_unc *= 0.21982248520710063; return;}
    else if (glu_mass == 2215) { xsec *= 0.00033; xsec_unc *= 0.22030303030303033; return;}
    else if (glu_mass == 2220) { xsec *= 0.000321; xsec_unc *= 0.22118380062305298; return;}
    else if (glu_mass == 2225) { xsec *= 0.000313; xsec_unc *= 0.22172523961661342; return;}
    else if (glu_mass == 2230) { xsec *= 0.000305; xsec_unc *= 0.22262295081967212; return;}
    else if (glu_mass == 2235) { xsec *= 0.000297; xsec_unc *= 0.22323232323232323; return;}
    else if (glu_mass == 2240) { xsec *= 0.00029; xsec_unc *= 0.22413793103448273; return;}
    else if (glu_mass == 2245) { xsec *= 0.000283; xsec_unc *= 0.2247349823321555; return;}
    else if (glu_mass == 2250) { xsec *= 0.000275; xsec_unc *= 0.22545454545454546; return;}
    else if (glu_mass == 2255) { xsec *= 0.000268; xsec_unc *= 0.22611940298507463; return;}
    else if (glu_mass == 2260) { xsec *= 0.000262; xsec_unc *= 0.22709923664122136; return;}
    else if (glu_mass == 2265) { xsec *= 0.000255; xsec_unc *= 0.22784313725490196; return;}
    else if (glu_mass == 2270) { xsec *= 0.000248; xsec_unc *= 0.22862903225806452; return;}
    else if (glu_mass == 2275) { xsec *= 0.000242; xsec_unc *= 0.22933884297520662; return;}
    else if (glu_mass == 2280) { xsec *= 0.000236; xsec_unc *= 0.23008474576271187; return;}
    else if (glu_mass == 2285) { xsec *= 0.00023; xsec_unc *= 0.2308695652173913; return;}
    else if (glu_mass == 2290) { xsec *= 0.000224; xsec_unc *= 0.23169642857142858; return;}
    else if (glu_mass == 2295) { xsec *= 0.000219; xsec_unc *= 0.23242009132420088; return;}
    else if (glu_mass == 2300) { xsec *= 0.000213; xsec_unc *= 0.23333333333333334; return;}
    else if (glu_mass == 2305) { xsec *= 0.000208; xsec_unc *= 0.23413461538461539; return;}
    else if (glu_mass == 2310) { xsec *= 0.000202; xsec_unc *= 0.23514851485148516; return;}
    else if (glu_mass == 2315) { xsec *= 0.000197; xsec_unc *= 0.23604060913705585; return;}
    else if (glu_mass == 2320) { xsec *= 0.000192; xsec_unc *= 0.23697916666666666; return;}
    else if (glu_mass == 2325) { xsec *= 0.000187; xsec_unc *= 0.23743315508021393; return;}
    else if (glu_mass == 2330) { xsec *= 0.000183; xsec_unc *= 0.23879781420765026; return;}
    else if (glu_mass == 2335) { xsec *= 0.000178; xsec_unc *= 0.2393258426966292; return;}
    else if (glu_mass == 2340) { xsec *= 0.000174; xsec_unc *= 0.24022988505747125; return;}
    else if (glu_mass == 2345) { xsec *= 0.000169; xsec_unc *= 0.2414201183431953; return;}
    else if (glu_mass == 2350) { xsec *= 0.000165; xsec_unc *= 0.24242424242424246; return;}
    else if (glu_mass == 2355) { xsec *= 0.000161; xsec_unc *= 0.24285714285714285; return;}
    else if (glu_mass == 2360) { xsec *= 0.000157; xsec_unc *= 0.2439490445859873; return;}
    else if (glu_mass == 2365) { xsec *= 0.000153; xsec_unc *= 0.24509803921568624; return;}
    else if (glu_mass == 2370) { xsec *= 0.000149; xsec_unc *= 0.24630872483221478; return;}
    else if (glu_mass == 2375) { xsec *= 0.000145; xsec_unc *= 0.24689655172413794; return;}
    else if (glu_mass == 2380) { xsec *= 0.000142; xsec_unc *= 0.24788732394366197; return;}
    else if (glu_mass == 2385) { xsec *= 0.000138; xsec_unc *= 0.2485507246376812; return;}
    else if (glu_mass == 2390) { xsec *= 0.000134; xsec_unc *= 0.25; return;}
    else if (glu_mass == 2395) { xsec *= 0.000131; xsec_unc *= 0.2511450381679389; return;}
    else if (glu_mass == 2400) { xsec *= 0.000128; xsec_unc *= 0.25156249999999997; return;}
    else if (glu_mass == 2405) { xsec *= 0.000125; xsec_unc *= 0.2528; return;}
    else if (glu_mass == 2410) { xsec *= 0.000121; xsec_unc *= 0.2537190082644628; return;}
    else if (glu_mass == 2415) { xsec *= 0.000118; xsec_unc *= 0.2550847457627119; return;}
    else if (glu_mass == 2420) { xsec *= 0.000115; xsec_unc *= 0.25652173913043474; return;}
    else if (glu_mass == 2425) { xsec *= 0.000113; xsec_unc *= 0.2575221238938053; return;}
    else if (glu_mass == 2430) { xsec *= 0.00011; xsec_unc *= 0.2581818181818182; return;}
    else if (glu_mass == 2435) { xsec *= 0.000107; xsec_unc *= 0.2588785046728972; return;}
    else if (glu_mass == 2440) { xsec *= 0.000104; xsec_unc *= 0.2605769230769231; return;}
    else if (glu_mass == 2445) { xsec *= 0.000102; xsec_unc *= 0.26176470588235295; return;}
    else if (glu_mass == 2450) { xsec *= 9.91e-05; xsec_unc *= 0.2623612512613522; return;}
    else if (glu_mass == 2455) { xsec *= 9.66e-05; xsec_unc *= 0.2639751552795031; return;}
    else if (glu_mass == 2460) { xsec *= 9.41e-05; xsec_unc *= 0.26461211477151964; return;}
    else if (glu_mass == 2465) { xsec *= 9.18e-05; xsec_unc *= 0.2657952069716776; return;}
    else if (glu_mass == 2470) { xsec *= 8.95e-05; xsec_unc *= 0.26703910614525145; return;}
    else if (glu_mass == 2475) { xsec *= 8.72e-05; xsec_unc *= 0.268348623853211; return;}
    else if (glu_mass == 2480) { xsec *= 8.5e-05; xsec_unc *= 0.26941176470588235; return;}
    else if (glu_mass == 2485) { xsec *= 8.29e-05; xsec_unc *= 0.27020506634499397; return;}
    else if (glu_mass == 2490) { xsec *= 8.08e-05; xsec_unc *= 0.2722772277227723; return;}
    else if (glu_mass == 2495) { xsec *= 7.88e-05; xsec_unc *= 0.27284263959390864; return;}
    else if (glu_mass == 2500) { xsec *= 7.68e-05; xsec_unc *= 0.27473958333333337; return;}
    else if (glu_mass == 2505) { xsec *= 7.49e-05; xsec_unc *= 0.27636849132176233; return;}
    else if (glu_mass == 2510) { xsec *= 7.3e-05; xsec_unc *= 0.27671232876712326; return;}
    else if (glu_mass == 2515) { xsec *= 7.12e-05; xsec_unc *= 0.27808988764044945; return;}
    else if (glu_mass == 2520) { xsec *= 6.94e-05; xsec_unc *= 0.2795389048991354; return;}
    else if (glu_mass == 2525) { xsec *= 6.77e-05; xsec_unc *= 0.28064992614475626; return;}
    else if (glu_mass == 2530) { xsec *= 6.6e-05; xsec_unc *= 0.2818181818181818; return;}
    else if (glu_mass == 2535) { xsec *= 6.43e-05; xsec_unc *= 0.28304821150855364; return;}
    else if (glu_mass == 2540) { xsec *= 6.27e-05; xsec_unc *= 0.28548644338118023; return;}
    else if (glu_mass == 2545) { xsec *= 6.11e-05; xsec_unc *= 0.2864157119476268; return;}
    else if (glu_mass == 2550) { xsec *= 5.96e-05; xsec_unc *= 0.28859060402684567; return;}
    else if (glu_mass == 2555) { xsec *= 5.81e-05; xsec_unc *= 0.28915662650602403; return;}
    else if (glu_mass == 2560) { xsec *= 5.66e-05; xsec_unc *= 0.2915194346289753; return;}
    else if (glu_mass == 2565) { xsec *= 5.52e-05; xsec_unc *= 0.29166666666666663; return;}
    else if (glu_mass == 2570) { xsec *= 5.38e-05; xsec_unc *= 0.2936802973977695; return;}
    else if (glu_mass == 2575) { xsec *= 5.25e-05; xsec_unc *= 0.29523809523809524; return;}
    else if (glu_mass == 2580) { xsec *= 5.12e-05; xsec_unc *= 0.296875; return;}
    else if (glu_mass == 2585) { xsec *= 4.99e-05; xsec_unc *= 0.2985971943887776; return;}
    else if (glu_mass == 2590) { xsec *= 4.86e-05; xsec_unc *= 0.3004115226337449; return;}
    else if (glu_mass == 2595) { xsec *= 4.74e-05; xsec_unc *= 0.30168776371308015; return;}
    else if (glu_mass == 2600) { xsec *= 4.62e-05; xsec_unc *= 0.30303030303030304; return;}
    else if (glu_mass == 2605) { xsec *= 4.51e-05; xsec_unc *= 0.30376940133037694; return;}
    else if (glu_mass == 2610) { xsec *= 4.39e-05; xsec_unc *= 0.3052391799544419; return;}
    else if (glu_mass == 2615) { xsec *= 4.28e-05; xsec_unc *= 0.30841121495327106; return;}
    else if (glu_mass == 2620) { xsec *= 4.18e-05; xsec_unc *= 0.30861244019138756; return;}
    else if (glu_mass == 2625) { xsec *= 4.07e-05; xsec_unc *= 0.31203931203931207; return;}
    else if (glu_mass == 2630) { xsec *= 3.97e-05; xsec_unc *= 0.31234256926952136; return;}
    else if (glu_mass == 2635) { xsec *= 3.87e-05; xsec_unc *= 0.3152454780361757; return;}
    else if (glu_mass == 2640) { xsec *= 3.77e-05; xsec_unc *= 0.3156498673740053; return;}
    else if (glu_mass == 2645) { xsec *= 3.68e-05; xsec_unc *= 0.3179347826086956; return;}
    else if (glu_mass == 2650) { xsec *= 3.59e-05; xsec_unc *= 0.3203342618384401; return;}
    else if (glu_mass == 2655) { xsec *= 3.5e-05; xsec_unc *= 0.3228571428571429; return;}
    else if (glu_mass == 2660) { xsec *= 3.41e-05; xsec_unc *= 0.3225806451612903; return;}
    else if (glu_mass == 2665) { xsec *= 3.32e-05; xsec_unc *= 0.3253012048192771; return;}
    else if (glu_mass == 2670) { xsec *= 3.24e-05; xsec_unc *= 0.3271604938271605; return;}
    else if (glu_mass == 2675) { xsec *= 3.16e-05; xsec_unc *= 0.3291139240506329; return;}
    else if (glu_mass == 2680) { xsec *= 3.08e-05; xsec_unc *= 0.33116883116883117; return;}
    else if (glu_mass == 2685) { xsec *= 3e-05; xsec_unc *= 0.33333333333333337; return;}
    else if (glu_mass == 2690) { xsec *= 2.93e-05; xsec_unc *= 0.3351535836177474; return;}
    else if (glu_mass == 2695) { xsec *= 2.85e-05; xsec_unc *= 0.3371929824561403; return;}
    else if (glu_mass == 2700) { xsec *= 2.78e-05; xsec_unc *= 0.33920863309352517; return;}
    else if (glu_mass == 2705) { xsec *= 2.71e-05; xsec_unc *= 0.3413284132841328; return;}
    else if (glu_mass == 2710) { xsec *= 2.65e-05; xsec_unc *= 0.3433962264150943; return;}
    else if (glu_mass == 2715) { xsec *= 2.58e-05; xsec_unc *= 0.3457364341085271; return;}
    else if (glu_mass == 2720) { xsec *= 2.51e-05; xsec_unc *= 0.347808764940239; return;}
    else if (glu_mass == 2725) { xsec *= 2.45e-05; xsec_unc *= 0.3497959183673469; return;}
    else if (glu_mass == 2730) { xsec *= 2.39e-05; xsec_unc *= 0.3523012552301255; return;}
    else if (glu_mass == 2735) { xsec *= 2.33e-05; xsec_unc *= 0.35450643776824037; return;}
    else if (glu_mass == 2740) { xsec *= 2.27e-05; xsec_unc *= 0.3568281938325991; return;}
    else if (glu_mass == 2745) { xsec *= 2.21e-05; xsec_unc *= 0.35882352941176476; return;}
    else if (glu_mass == 2750) { xsec *= 2.16e-05; xsec_unc *= 0.3611111111111111; return;}
    else if (glu_mass == 2755) { xsec *= 2.11e-05; xsec_unc *= 0.36350710900473926; return;}
    else if (glu_mass == 2760) { xsec *= 2.05e-05; xsec_unc *= 0.36585365853658536; return;}
    else if (glu_mass == 2765) { xsec *= 2e-05; xsec_unc *= 0.36799999999999994; return;}
    else if (glu_mass == 2770) { xsec *= 1.95e-05; xsec_unc *= 0.37025641025641026; return;}
    else if (glu_mass == 2775) { xsec *= 1.9e-05; xsec_unc *= 0.3731578947368421; return;}
    else if (glu_mass == 2780) { xsec *= 1.85e-05; xsec_unc *= 0.37513513513513513; return;}
    else if (glu_mass == 2785) { xsec *= 1.81e-05; xsec_unc *= 0.37734806629834255; return;}
    else if (glu_mass == 2790) { xsec *= 1.76e-05; xsec_unc *= 0.3801136363636364; return;}
    else if (glu_mass == 2795) { xsec *= 1.72e-05; xsec_unc *= 0.38255813953488366; return;}
    else if (glu_mass == 2800) { xsec *= 1.68e-05; xsec_unc *= 0.38452380952380955; return;}
    else if (glu_mass == 2805) { xsec *= 1.63e-05; xsec_unc *= 0.3871165644171779; return;}
    else if (glu_mass == 2810) { xsec *= 1.59e-05; xsec_unc *= 0.3893081761006289; return;}
    else if (glu_mass == 2815) { xsec *= 1.55e-05; xsec_unc *= 0.39161290322580644; return;}
    else if (glu_mass == 2820) { xsec *= 1.51e-05; xsec_unc *= 0.39403973509933776; return;}
    else if (glu_mass == 2825) { xsec *= 1.48e-05; xsec_unc *= 0.39662162162162157; return;}
    else if (glu_mass == 2830) { xsec *= 1.44e-05; xsec_unc *= 0.3993055555555556; return;}
    else if (glu_mass == 2835) { xsec *= 1.4e-05; xsec_unc *= 0.40142857142857147; return;}
    else if (glu_mass == 2840) { xsec *= 1.37e-05; xsec_unc *= 0.4036496350364964; return;}
    else if (glu_mass == 2845) { xsec *= 1.33e-05; xsec_unc *= 0.406015037593985; return;}
    else if (glu_mass == 2850) { xsec *= 1.3e-05; xsec_unc *= 0.4084615384615385; return;}
    else if (glu_mass == 2855) { xsec *= 1.27e-05; xsec_unc *= 0.4110236220472441; return;}
    else if (glu_mass == 2860) { xsec *= 1.24e-05; xsec_unc *= 0.4137096774193548; return;}
    else if (glu_mass == 2865) { xsec *= 1.21e-05; xsec_unc *= 0.415702479338843; return;}
    else if (glu_mass == 2870) { xsec *= 1.18e-05; xsec_unc *= 0.4186440677966102; return;}
    else if (glu_mass == 2875) { xsec *= 1.15e-05; xsec_unc *= 0.42086956521739133; return;}
    else if (glu_mass == 2880) { xsec *= 1.12e-05; xsec_unc *= 0.42321428571428577; return;}
    else if (glu_mass == 2885) { xsec *= 1.09e-05; xsec_unc *= 0.4256880733944953; return;}
    else if (glu_mass == 2890) { xsec *= 1.06e-05; xsec_unc *= 0.4283018867924528; return;}
    else if (glu_mass == 2895) { xsec *= 1.04e-05; xsec_unc *= 0.4307692307692308; return;}
    else if (glu_mass == 2900) { xsec *= 1.01e-05; xsec_unc *= 0.43267326732673267; return;}
    else if (glu_mass == 2905) { xsec *= 9.86e-06; xsec_unc *= 0.4350912778904665; return;}
    else if (glu_mass == 2910) { xsec *= 9.61e-06; xsec_unc *= 0.43808532778355885; return;}
    else if (glu_mass == 2915) { xsec *= 9.37e-06; xsec_unc *= 0.4407684098185699; return;}
    else if (glu_mass == 2920) { xsec *= 9.14e-06; xsec_unc *= 0.4431072210065645; return;}
    else if (glu_mass == 2925) { xsec *= 8.91e-06; xsec_unc *= 0.44556677890011226; return;}
    else if (glu_mass == 2930) { xsec *= 8.69e-06; xsec_unc *= 0.44764096662830843; return;}
    else if (glu_mass == 2935) { xsec *= 8.48e-06; xsec_unc *= 0.45047169811320753; return;}
    else if (glu_mass == 2940) { xsec *= 8.27e-06; xsec_unc *= 0.45344619105199513; return;}
    else if (glu_mass == 2945) { xsec *= 8.06e-06; xsec_unc *= 0.4553349875930521; return;}
    else if (glu_mass == 2950) { xsec *= 7.86e-06; xsec_unc *= 0.45801526717557256; return;}
    else if (glu_mass == 2955) { xsec *= 7.67e-06; xsec_unc *= 0.4602346805736637; return;}
    else if (glu_mass == 2960) { xsec *= 7.48e-06; xsec_unc *= 0.46256684491978606; return;}
    else if (glu_mass == 2965) { xsec *= 7.29e-06; xsec_unc *= 0.4650205761316873; return;}
    else if (glu_mass == 2970) { xsec *= 7.11e-06; xsec_unc *= 0.46835443037974683; return;}
    else if (glu_mass == 2975) { xsec *= 6.94e-06; xsec_unc *= 0.47118155619596547; return;}
    else if (glu_mass == 2980) { xsec *= 6.77e-06; xsec_unc *= 0.4726735598227474; return;}
    else if (glu_mass == 2985) { xsec *= 6.6e-06; xsec_unc *= 0.4757575757575757; return;}
    else if (glu_mass == 2990) { xsec *= 6.44e-06; xsec_unc *= 0.4782608695652174; return;}
    else if (glu_mass == 2995) { xsec *= 6.28e-06; xsec_unc *= 0.48089171974522293; return;}
    else if (glu_mass == 3000) { xsec *= 6.12e-06; xsec_unc *= 0.48366013071895425; return;}
    else { xsec = 0; xsec_unc = 0; }
  }

  void stopCrossSection(int stop_mass, double &xsec, double &xsec_unc){
    
    if(stop_mass == 100) { xsec = 1521.11; xsec_unc = 0.154038; return; }
    else if(stop_mass == 105) { xsec = 1233.18; xsec_unc = 0.154059; return; }
    else if(stop_mass == 110) { xsec = 1013.76; xsec_unc = 0.154088; return; }
    else if(stop_mass == 115) { xsec = 832.656; xsec_unc = 0.151503; return; }
    else if(stop_mass == 120) { xsec = 689.799; xsec_unc = 0.15044; return; }
    else if(stop_mass == 125) { xsec = 574.981; xsec_unc = 0.149895; return; }
    else if(stop_mass == 130) { xsec = 481.397; xsec_unc = 0.148906; return; }
    else if(stop_mass == 135) { xsec = 405.159; xsec_unc = 0.148952; return; }
    else if(stop_mass == 140) { xsec = 342.865; xsec_unc = 0.149119; return; }
    else if(stop_mass == 145) { xsec = 291.752; xsec_unc = 0.148022; return; }
    else if(stop_mass == 150) { xsec = 249.409; xsec_unc = 0.147477; return; }
    else if(stop_mass == 155) { xsec = 214.221; xsec_unc = 0.145928; return; }
    else if(stop_mass == 160) { xsec = 184.623; xsec_unc = 0.145821; return; }
    else if(stop_mass == 165) { xsec = 159.614; xsec_unc = 0.147859; return; }
    else if(stop_mass == 170) { xsec = 139.252; xsec_unc = 0.14547; return; }
    else if(stop_mass == 175) { xsec = 121.416; xsec_unc = 0.146341; return; }
    else if(stop_mass == 180) { xsec = 106.194; xsec_unc = 0.142033; return; }
    else if(stop_mass == 185) { xsec = 93.3347; xsec_unc = 0.144893; return; }
    else if(stop_mass == 190) { xsec = 82.2541; xsec_unc = 0.144677; return; }
    else if(stop_mass == 195) { xsec = 72.7397; xsec_unc = 0.144452; return; }
    else if(stop_mass == 200) { xsec = 64.5085; xsec_unc = 0.144098; return; }
    else if(stop_mass == 205) { xsec = 57.2279; xsec_unc = 0.144191; return; }
    else if(stop_mass == 210) { xsec = 50.9226; xsec_unc = 0.142457; return; }
    else if(stop_mass == 215) { xsec = 45.3761; xsec_unc = 0.14344; return; }
    else if(stop_mass == 220) { xsec = 40.5941; xsec_unc = 0.142634; return; }
    else if(stop_mass == 225) { xsec = 36.3818; xsec_unc = 0.142189; return; }
    else if(stop_mass == 230) { xsec = 32.6679; xsec_unc = 0.141592; return; }
    else if(stop_mass == 235) { xsec = 29.3155; xsec_unc = 0.142233; return; }
    else if(stop_mass == 240) { xsec = 26.4761; xsec_unc = 0.141723; return; }
    else if(stop_mass == 245) { xsec = 23.8853; xsec_unc = 0.139482; return; }
    else if(stop_mass == 250) { xsec = 21.5949; xsec_unc = 0.140595; return; }
    else if(stop_mass == 255) { xsec = 19.5614; xsec_unc = 0.138755; return; }
    else if(stop_mass == 260) { xsec = 17.6836; xsec_unc = 0.139505; return; }
    else if(stop_mass == 265) { xsec = 16.112; xsec_unc = 0.139531; return; }
    else if(stop_mass == 270) { xsec = 14.6459; xsec_unc = 0.139278; return; }
    else if(stop_mass == 275) { xsec = 13.3231; xsec_unc = 0.142549; return; }
    else if(stop_mass == 280) { xsec = 12.1575; xsec_unc = 0.141584; return; }
    else if(stop_mass == 285) { xsec = 11.0925; xsec_unc = 0.140904; return; }
    else if(stop_mass == 290) { xsec = 10.1363; xsec_unc = 0.138967; return; }
    else if(stop_mass == 295) { xsec = 9.29002; xsec_unc = 0.139107; return; }
    else if(stop_mass == 300) { xsec = 8.51615; xsec_unc = 0.139223; return; }
    else if(stop_mass == 305) { xsec = 7.81428; xsec_unc = 0.138996; return; }
    else if(stop_mass == 310) { xsec = 7.17876; xsec_unc = 0.139357; return; }
    else if(stop_mass == 315) { xsec = 6.60266; xsec_unc = 0.139256; return; }
    else if(stop_mass == 320) { xsec = 6.08444; xsec_unc = 0.137957; return; }
    else if(stop_mass == 325) { xsec = 5.60471; xsec_unc = 0.138144; return; }
    else if(stop_mass == 330) { xsec = 5.17188; xsec_unc = 0.136954; return; }
    else if(stop_mass == 335) { xsec = 4.77871; xsec_unc = 0.137554; return; }
    else if(stop_mass == 340) { xsec = 4.41629; xsec_unc = 0.137945; return; }
    else if(stop_mass == 345) { xsec = 4.08881; xsec_unc = 0.137075; return; }
    else if(stop_mass == 350) { xsec = 3.78661; xsec_unc = 0.136877; return; }
    else if(stop_mass == 355) { xsec = 3.50911; xsec_unc = 0.138089; return; }
    else if(stop_mass == 360) { xsec = 3.25619; xsec_unc = 0.138002; return; }
    else if(stop_mass == 365) { xsec = 3.02472; xsec_unc = 0.137093; return; }
    else if(stop_mass == 370) { xsec = 2.8077; xsec_unc = 0.138064; return; }
    else if(stop_mass == 375) { xsec = 2.61162; xsec_unc = 0.138477; return; }
    else if(stop_mass == 380) { xsec = 2.43031; xsec_unc = 0.136999; return; }
    else if(stop_mass == 385) { xsec = 2.26365; xsec_unc = 0.13728; return; }
    else if(stop_mass == 390) { xsec = 2.10786; xsec_unc = 0.13732; return; }
    else if(stop_mass == 395) { xsec = 1.9665; xsec_unc = 0.134737; return; }
    else if(stop_mass == 400) { xsec = 1.83537; xsec_unc = 0.136985; return; }
    else if(stop_mass == 405) { xsec = 1.70927; xsec_unc = 0.137114; return; }
    else if(stop_mass == 410) { xsec = 1.60378; xsec_unc = 0.135468; return; }
    else if(stop_mass == 415) { xsec = 1.49798; xsec_unc = 0.134453; return; }
    else if(stop_mass == 420) { xsec = 1.39688; xsec_unc = 0.136719; return; }
    else if(stop_mass == 425) { xsec = 1.31169; xsec_unc = 0.135013; return; }
    else if(stop_mass == 430) { xsec = 1.22589; xsec_unc = 0.133237; return; }
    else if(stop_mass == 435) { xsec = 1.14553; xsec_unc = 0.135478; return; }
    else if(stop_mass == 440) { xsec = 1.07484; xsec_unc = 0.137238; return; }
    else if(stop_mass == 445) { xsec = 1.01019; xsec_unc = 0.134187; return; }
    else if(stop_mass == 450) { xsec = 0.948333; xsec_unc = 0.134559; return; }
    else if(stop_mass == 455) { xsec = 0.890847; xsec_unc = 0.134587; return; }
    else if(stop_mass == 460) { xsec = 0.836762; xsec_unc = 0.134468; return; }
    else if(stop_mass == 465) { xsec = 0.787221; xsec_unc = 0.134149; return; }
    else if(stop_mass == 470) { xsec = 0.740549; xsec_unc = 0.134127; return; }
    else if(stop_mass == 475) { xsec = 0.697075; xsec_unc = 0.133926; return; }
    else if(stop_mass == 480) { xsec = 0.655954; xsec_unc = 0.134392; return; }
    else if(stop_mass == 485) { xsec = 0.618562; xsec_unc = 0.133705; return; }
    else if(stop_mass == 490) { xsec = 0.582467; xsec_unc = 0.133914; return; }
    else if(stop_mass == 495) { xsec = 0.549524; xsec_unc = 0.133691; return; }
    else if(stop_mass == 500) { xsec = 0.51848; xsec_unc = 0.133797; return; }
    else if(stop_mass == 505) { xsec = 0.489324; xsec_unc = 0.133608; return; }
    else if(stop_mass == 510) { xsec = 0.462439; xsec_unc = 0.133046; return; }
    else if(stop_mass == 515) { xsec = 0.436832; xsec_unc = 0.133703; return; }
    else if(stop_mass == 520) { xsec = 0.412828; xsec_unc = 0.13272; return; }
    else if(stop_mass == 525) { xsec = 0.390303; xsec_unc = 0.133443; return; }
    else if(stop_mass == 530) { xsec = 0.368755; xsec_unc = 0.133769; return; }
    else if(stop_mass == 535) { xsec = 0.348705; xsec_unc = 0.132706; return; }
    else if(stop_mass == 540) { xsec = 0.330157; xsec_unc = 0.132981; return; }
    else if(stop_mass == 545) { xsec = 0.312672; xsec_unc = 0.13277; return; }
    else if(stop_mass == 550) { xsec = 0.296128; xsec_unc = 0.132687; return; }
    else if(stop_mass == 555) { xsec = 0.280734; xsec_unc = 0.132363; return; }
    else if(stop_mass == 560) { xsec = 0.266138; xsec_unc = 0.13193; return; }
    else if(stop_mass == 565) { xsec = 0.251557; xsec_unc = 0.131731; return; }
    else if(stop_mass == 570) { xsec = 0.238537; xsec_unc = 0.133409; return; }
    else if(stop_mass == 575) { xsec = 0.226118; xsec_unc = 0.132741; return; }
    else if(stop_mass == 580) { xsec = 0.214557; xsec_unc = 0.131697; return; }
    else if(stop_mass == 585) { xsec = 0.203566; xsec_unc = 0.133257; return; }
    else if(stop_mass == 590) { xsec = 0.193079; xsec_unc = 0.132037; return; }
    else if(stop_mass == 595) { xsec = 0.183604; xsec_unc = 0.130973; return; }
    else if(stop_mass == 600) { xsec = 0.174599; xsec_unc = 0.132074; return; }
    else if(stop_mass == 605) { xsec = 0.166131; xsec_unc = 0.130154; return; }
    else if(stop_mass == 610) { xsec = 0.158242; xsec_unc = 0.13142; return; }
    else if(stop_mass == 615) { xsec = 0.150275; xsec_unc = 0.13285; return; }
    else if(stop_mass == 620) { xsec = 0.142787; xsec_unc = 0.130642; return; }
    else if(stop_mass == 625) { xsec = 0.136372; xsec_unc = 0.127962; return; }
    else if(stop_mass == 630) { xsec = 0.129886; xsec_unc = 0.132957; return; }
    else if(stop_mass == 635) { xsec = 0.123402; xsec_unc = 0.13016; return; }
    else if(stop_mass == 640) { xsec = 0.11795; xsec_unc = 0.127132; return; }
    else if(stop_mass == 645) { xsec = 0.112008; xsec_unc = 0.12808; return; }
    else if(stop_mass == 650) { xsec = 0.107045; xsec_unc = 0.129232; return; }
    else if(stop_mass == 655) { xsec = 0.102081; xsec_unc = 0.130012; return; }
    else if(stop_mass == 660) { xsec = 0.09725; xsec_unc = 0.129038; return; }
    else if(stop_mass == 665) { xsec = 0.0927515; xsec_unc = 0.129548; return; }
    else if(stop_mass == 670) { xsec = 0.0885084; xsec_unc = 0.130218; return; }
    else if(stop_mass == 675) { xsec = 0.0844877; xsec_unc = 0.130703; return; }
    else if(stop_mass == 680) { xsec = 0.0806192; xsec_unc = 0.131131; return; }
    else if(stop_mass == 685) { xsec = 0.0769099; xsec_unc = 0.131517; return; }
    else if(stop_mass == 690) { xsec = 0.0734901; xsec_unc = 0.132344; return; }
    else if(stop_mass == 695) { xsec = 0.0701805; xsec_unc = 0.132716; return; }
    else if(stop_mass == 700) { xsec = 0.0670476; xsec_unc = 0.133429; return; }
    else if(stop_mass == 705) { xsec = 0.0641426; xsec_unc = 0.13363; return; }
    else if(stop_mass == 710) { xsec = 0.0612942; xsec_unc = 0.133941; return; }
    else if(stop_mass == 715) { xsec = 0.0585678; xsec_unc = 0.134663; return; }
    else if(stop_mass == 720) { xsec = 0.0560753; xsec_unc = 0.134984; return; }
    else if(stop_mass == 725) { xsec = 0.0536438; xsec_unc = 0.135804; return; }
    else if(stop_mass == 730) { xsec = 0.0513219; xsec_unc = 0.135682; return; }
    else if(stop_mass == 735) { xsec = 0.0491001; xsec_unc = 0.136268; return; }
    else if(stop_mass == 740) { xsec = 0.0470801; xsec_unc = 0.136895; return; }
    else if(stop_mass == 745) { xsec = 0.045061; xsec_unc = 0.136816; return; }
    else if(stop_mass == 750) { xsec = 0.0431418; xsec_unc = 0.137455; return; }
    else if(stop_mass == 755) { xsec = 0.0413447; xsec_unc = 0.137833; return; }
    else if(stop_mass == 760) { xsec = 0.0396264; xsec_unc = 0.138518; return; }
    else if(stop_mass == 765) { xsec = 0.0379036; xsec_unc = 0.138537; return; }
    else if(stop_mass == 770) { xsec = 0.0363856; xsec_unc = 0.139334; return; }
    else if(stop_mass == 775) { xsec = 0.0348796; xsec_unc = 0.139597; return; }
    else if(stop_mass == 780) { xsec = 0.0334669; xsec_unc = 0.140267; return; }
    else if(stop_mass == 785) { xsec = 0.0320548; xsec_unc = 0.140406; return; }
    else if(stop_mass == 790) { xsec = 0.0307373; xsec_unc = 0.14115; return; }
    else if(stop_mass == 795) { xsec = 0.0295348; xsec_unc = 0.141397; return; }
    else if(stop_mass == 800) { xsec = 0.0283338; xsec_unc = 0.14171; return; }
    else if(stop_mass == 805) { xsec = 0.0272206; xsec_unc = 0.14241; return; }
    else if(stop_mass == 810) { xsec = 0.0261233; xsec_unc = 0.142891; return; }
    else if(stop_mass == 815) { xsec = 0.0251107; xsec_unc = 0.143632; return; }
    else if(stop_mass == 820) { xsec = 0.0241099; xsec_unc = 0.143805; return; }
    else if(stop_mass == 825) { xsec = 0.0230866; xsec_unc = 0.144428; return; }
    else if(stop_mass == 830) { xsec = 0.0221834; xsec_unc = 0.144791; return; }
    else if(stop_mass == 835) { xsec = 0.0213766; xsec_unc = 0.145511; return; }
    else if(stop_mass == 840) { xsec = 0.0204715; xsec_unc = 0.146131; return; }
    else if(stop_mass == 845) { xsec = 0.0197653; xsec_unc = 0.146602; return; }
    else if(stop_mass == 850) { xsec = 0.0189612; xsec_unc = 0.14702; return; }
    else if(stop_mass == 855) { xsec = 0.0182516; xsec_unc = 0.147648; return; }
    else if(stop_mass == 860) { xsec = 0.0175509; xsec_unc = 0.147944; return; }
    else if(stop_mass == 865) { xsec = 0.0168336; xsec_unc = 0.148528; return; }
    else if(stop_mass == 870) { xsec = 0.0162314; xsec_unc = 0.148772; return; }
    else if(stop_mass == 875) { xsec = 0.015625; xsec_unc = 0.149567; return; }
    else if(stop_mass == 880) { xsec = 0.0150143; xsec_unc = 0.150389; return; }
    else if(stop_mass == 885) { xsec = 0.0144112; xsec_unc = 0.150614; return; }
    else if(stop_mass == 890) { xsec = 0.0138979; xsec_unc = 0.151; return; }
    else if(stop_mass == 895) { xsec = 0.0133962; xsec_unc = 0.151325; return; }
    else if(stop_mass == 900) { xsec = 0.0128895; xsec_unc = 0.152026; return; }
    else if(stop_mass == 905) { xsec = 0.0123843; xsec_unc = 0.152968; return; }
    else if(stop_mass == 910) { xsec = 0.0119837; xsec_unc = 0.153089; return; }
    else if(stop_mass == 915) { xsec = 0.0114713; xsec_unc = 0.153678; return; }
    else if(stop_mass == 920) { xsec = 0.0110688; xsec_unc = 0.154082; return; }
    else if(stop_mass == 925) { xsec = 0.0106631; xsec_unc = 0.154806; return; }
    else if(stop_mass == 930) { xsec = 0.0102629; xsec_unc = 0.155313; return; }
    else if(stop_mass == 935) { xsec = 0.0098874; xsec_unc = 0.156066; return; }
    else if(stop_mass == 940) { xsec = 0.00952142; xsec_unc = 0.156055; return; }
    else if(stop_mass == 945) { xsec = 0.00916636; xsec_unc = 0.156849; return; }
    else if(stop_mass == 950) { xsec = 0.00883465; xsec_unc = 0.157177; return; }
    else if(stop_mass == 955) { xsec = 0.00851073; xsec_unc = 0.158094; return; }
    else if(stop_mass == 960) { xsec = 0.00820884; xsec_unc = 0.15844; return; }
    else if(stop_mass == 965) { xsec = 0.00791403; xsec_unc = 0.159216; return; }
    else if(stop_mass == 970) { xsec = 0.00763112; xsec_unc = 0.159742; return; }
    else if(stop_mass == 975) { xsec = 0.00735655; xsec_unc = 0.160548; return; }
    else if(stop_mass == 980) { xsec = 0.00710317; xsec_unc = 0.160626; return; }
    else if(stop_mass == 985) { xsec = 0.00684867; xsec_unc = 0.16144; return; }
    else if(stop_mass == 990) { xsec = 0.00660695; xsec_unc = 0.161813; return; }
    else if(stop_mass == 995) { xsec = 0.00637546; xsec_unc = 0.162158; return; }
    else if(stop_mass == 1000) { xsec = 0.00615134; xsec_unc = 0.162953; return; }
    else if(stop_mass == 1005) { xsec = 0.00593765; xsec_unc = 0.163716; return; }
    else if(stop_mass == 1010) { xsec = 0.00572452; xsec_unc = 0.163857; return; }
    else if(stop_mass == 1015) { xsec = 0.00553094; xsec_unc = 0.164628; return; }
    else if(stop_mass == 1020) { xsec = 0.00533968; xsec_unc = 0.164963; return; }
    else if(stop_mass == 1025) { xsec = 0.00514619; xsec_unc = 0.165762; return; }
    else if(stop_mass == 1030) { xsec = 0.00497235; xsec_unc = 0.165838; return; }
    else if(stop_mass == 1035) { xsec = 0.00479906; xsec_unc = 0.166646; return; }
    else if(stop_mass == 1040) { xsec = 0.00463806; xsec_unc = 0.166947; return; }
    else if(stop_mass == 1045) { xsec = 0.00447537; xsec_unc = 0.167071; return; }
    else if(stop_mass == 1050) { xsec = 0.00432261; xsec_unc = 0.167859; return; }
    else if(stop_mass == 1055) { xsec = 0.00417983; xsec_unc = 0.168637; return; }
    else if(stop_mass == 1060) { xsec = 0.00403886; xsec_unc = 0.168981; return; }
    else if(stop_mass == 1065) { xsec = 0.0038962; xsec_unc = 0.169794; return; }
    else if(stop_mass == 1070) { xsec = 0.00376343; xsec_unc = 0.169764; return; }
    else if(stop_mass == 1075) { xsec = 0.00364174; xsec_unc = 0.170634; return; }
    else if(stop_mass == 1080) { xsec = 0.00352093; xsec_unc = 0.170908; return; }
    else if(stop_mass == 1085) { xsec = 0.00339813; xsec_unc = 0.171929; return; }
    else if(stop_mass == 1090) { xsec = 0.00328695; xsec_unc = 0.172274; return; }
    else if(stop_mass == 1095) { xsec = 0.00317628; xsec_unc = 0.172617; return; }
    else if(stop_mass == 1100) { xsec = 0.00307413; xsec_unc = 0.173377; return; }
    else if(stop_mass == 1105) { xsec = 0.00297377; xsec_unc = 0.173822; return; }
    else if(stop_mass == 1110) { xsec = 0.00287148; xsec_unc = 0.174725; return; }
    else if(stop_mass == 1115) { xsec = 0.00278078; xsec_unc = 0.175091; return; }
    else if(stop_mass == 1120) { xsec = 0.00268873; xsec_unc = 0.175883; return; }
    else if(stop_mass == 1125) { xsec = 0.00260821; xsec_unc = 0.176126; return; }
    else if(stop_mass == 1130) { xsec = 0.00251529; xsec_unc = 0.176836; return; }
    else if(stop_mass == 1135) { xsec = 0.00243484; xsec_unc = 0.177128; return; }
    else if(stop_mass == 1140) { xsec = 0.00236295; xsec_unc = 0.177977; return; }
    else if(stop_mass == 1145) { xsec = 0.00228192; xsec_unc = 0.178507; return; }
    else if(stop_mass == 1150) { xsec = 0.00221047; xsec_unc = 0.179259; return; }
    else if(stop_mass == 1155) { xsec = 0.00213907; xsec_unc = 0.180255; return; }
    else if(stop_mass == 1160) { xsec = 0.00206845; xsec_unc = 0.180518; return; }
    else if(stop_mass == 1165) { xsec = 0.0020063; xsec_unc = 0.180954; return; }
    else if(stop_mass == 1170) { xsec = 0.00194569; xsec_unc = 0.181194; return; }
    else if(stop_mass == 1175) { xsec = 0.0018741; xsec_unc = 0.182145; return; }
    else if(stop_mass == 1180) { xsec = 0.00182266; xsec_unc = 0.183074; return; }
    else if(stop_mass == 1185) { xsec = 0.00176211; xsec_unc = 0.183375; return; }
    else if(stop_mass == 1190) { xsec = 0.00170006; xsec_unc = 0.184075; return; }
    else if(stop_mass == 1195) { xsec = 0.00164968; xsec_unc = 0.184438; return; }
    else if(stop_mass == 1200) { xsec = 0.00159844; xsec_unc = 0.185209; return; }
    else if(stop_mass == 1205) { xsec = 0.0015472; xsec_unc = 0.185977; return; }
    else if(stop_mass == 1210) { xsec = 0.00149657; xsec_unc = 0.186485; return; }
    else if(stop_mass == 1215) { xsec = 0.00145544; xsec_unc = 0.187347; return; }
    else if(stop_mass == 1220) { xsec = 0.00140288; xsec_unc = 0.188774; return; }
    else if(stop_mass == 1225) { xsec = 0.00136155; xsec_unc = 0.18989; return; }
    else if(stop_mass == 1230) { xsec = 0.00131271; xsec_unc = 0.188763; return; }
    else if(stop_mass == 1235) { xsec = 0.0012717; xsec_unc = 0.189588; return; }
    else if(stop_mass == 1240) { xsec = 0.00123066; xsec_unc = 0.19049; return; }
    else if(stop_mass == 1245) { xsec = 0.00119994; xsec_unc = 0.191442; return; }
    else if(stop_mass == 1250) { xsec = 0.0011583; xsec_unc = 0.193006; return; }
    else if(stop_mass == 1255) { xsec = 0.00112694; xsec_unc = 0.194441; return; }
    else if(stop_mass == 1260) { xsec = 0.00108716; xsec_unc = 0.194141; return; }
    else if(stop_mass == 1265) { xsec = 0.00105517; xsec_unc = 0.196361; return; }
    else if(stop_mass == 1270) { xsec = 0.00102241; xsec_unc = 0.196297; return; }
    else if(stop_mass == 1275) { xsec = 0.000991293; xsec_unc = 0.19762; return; }
    else if(stop_mass == 1280) { xsec = 0.000961012; xsec_unc = 0.197926; return; }
    else if(stop_mass == 1285) { xsec = 0.000932394; xsec_unc = 0.198682; return; }
    else if(stop_mass == 1290) { xsec = 0.000903404; xsec_unc = 0.199924; return; }
    else if(stop_mass == 1295) { xsec = 0.000876957; xsec_unc = 0.200777; return; }
    else if(stop_mass == 1300) { xsec = 0.000850345; xsec_unc = 0.201604; return; }
    else if(stop_mass == 1305) { xsec = 0.00082443; xsec_unc = 0.202883; return; }
    else if(stop_mass == 1310) { xsec = 0.00079983; xsec_unc = 0.20373; return; }
    else if(stop_mass == 1315) { xsec = 0.000775222; xsec_unc = 0.204622; return; }
    else if(stop_mass == 1320) { xsec = 0.000751372; xsec_unc = 0.205919; return; }
    else if(stop_mass == 1325) { xsec = 0.000728912; xsec_unc = 0.206884; return; }
    else if(stop_mass == 1330) { xsec = 0.000706867; xsec_unc = 0.207763; return; }
    else if(stop_mass == 1335) { xsec = 0.000685372; xsec_unc = 0.208587; return; }
    else if(stop_mass == 1340) { xsec = 0.000664649; xsec_unc = 0.209879; return; }
    else if(stop_mass == 1345) { xsec = 0.000644804; xsec_unc = 0.211487; return; }
    else if(stop_mass == 1350) { xsec = 0.000625155; xsec_unc = 0.212761; return; }
    else if(stop_mass == 1355) { xsec = 0.000606802; xsec_unc = 0.213529; return; }
    else if(stop_mass == 1360) { xsec = 0.000588512; xsec_unc = 0.214428; return; }
    else if(stop_mass == 1365) { xsec = 0.000570506; xsec_unc = 0.216584; return; }
    else if(stop_mass == 1370) { xsec = 0.000553379; xsec_unc = 0.216036; return; }
    else if(stop_mass == 1375) { xsec = 0.000536646; xsec_unc = 0.21775; return; }
    else if(stop_mass == 1380) { xsec = 0.000521404; xsec_unc = 0.218383; return; }
    else if(stop_mass == 1385) { xsec = 0.000505008; xsec_unc = 0.219675; return; }
    else if(stop_mass == 1390) { xsec = 0.000490353; xsec_unc = 0.221444; return; }
    else if(stop_mass == 1395) { xsec = 0.000476164; xsec_unc = 0.222016; return; }
    else if(stop_mass == 1400) { xsec = 0.000461944; xsec_unc = 0.222704; return; }
    else if(stop_mass == 1405) { xsec = 0.000448172; xsec_unc = 0.224911; return; }
    else if(stop_mass == 1410) { xsec = 0.000435082; xsec_unc = 0.225606; return; }
    else if(stop_mass == 1415) { xsec = 0.000422967; xsec_unc = 0.226095; return; }
    else if(stop_mass == 1420) { xsec = 0.000410381; xsec_unc = 0.22797; return; }
    else if(stop_mass == 1425) { xsec = 0.000398106; xsec_unc = 0.228949; return; }
    else if(stop_mass == 1430) { xsec = 0.000386792; xsec_unc = 0.231319; return; }
    else if(stop_mass == 1435) { xsec = 0.000375724; xsec_unc = 0.231724; return; }
    else if(stop_mass == 1440) { xsec = 0.000364616; xsec_unc = 0.232234; return; }
    else if(stop_mass == 1445) { xsec = 0.000353965; xsec_unc = 0.234637; return; }
    else if(stop_mass == 1450) { xsec = 0.000343923; xsec_unc = 0.234948; return; }
    else if(stop_mass == 1455) { xsec = 0.000333885; xsec_unc = 0.235468; return; }
    else if(stop_mass == 1460) { xsec = 0.000324344; xsec_unc = 0.23771; return; }
    else if(stop_mass == 1465) { xsec = 0.0003153; xsec_unc = 0.238004; return; }
    else if(stop_mass == 1470) { xsec = 0.00030583; xsec_unc = 0.240064; return; }
    else if(stop_mass == 1475) { xsec = 0.000296811; xsec_unc = 0.240314; return; }
    else if(stop_mass == 1480) { xsec = 0.000288149; xsec_unc = 0.239248; return; }
    else if(stop_mass == 1485) { xsec = 0.000279711; xsec_unc = 0.241257; return; }
    else if(stop_mass == 1490) { xsec = 0.000271724; xsec_unc = 0.241274; return; }
    else if(stop_mass == 1495) { xsec = 0.000264275; xsec_unc = 0.243545; return; }
    else if(stop_mass == 1500) { xsec = 0.000256248; xsec_unc = 0.24372; return; }
    else if(stop_mass == 1505) { xsec = 0.000248853; xsec_unc = 0.245827; return; }
    else if(stop_mass == 1510) { xsec = 0.000241844; xsec_unc = 0.246187; return; }
    else if(stop_mass == 1515) { xsec = 0.000234438; xsec_unc = 0.248442; return; }
    else if(stop_mass == 1520) { xsec = 0.000227374; xsec_unc = 0.248909; return; }
    else if(stop_mass == 1525) { xsec = 0.000221045; xsec_unc = 0.250895; return; }
    else if(stop_mass == 1530) { xsec = 0.000214431; xsec_unc = 0.248728; return; }
    else if(stop_mass == 1535) { xsec = 0.000208092; xsec_unc = 0.251043; return; }
    else if(stop_mass == 1540) { xsec = 0.000201748; xsec_unc = 0.253207; return; }
    else if(stop_mass == 1545) { xsec = 0.000196399; xsec_unc = 0.255641; return; }
    else if(stop_mass == 1550) { xsec = 0.000190474; xsec_unc = 0.255213; return; }
    else if(stop_mass == 1555) { xsec = 0.000185188; xsec_unc = 0.257329; return; }
    else if(stop_mass == 1560) { xsec = 0.000179263; xsec_unc = 0.256931; return; }
    else if(stop_mass == 1565) { xsec = 0.000174021; xsec_unc = 0.259111; return; }
    else if(stop_mass == 1570) { xsec = 0.000169176; xsec_unc = 0.258106; return; }
    else if(stop_mass == 1575) { xsec = 0.000163861; xsec_unc = 0.260597; return; }
    else if(stop_mass == 1580) { xsec = 0.000159583; xsec_unc = 0.262958; return; }
    else if(stop_mass == 1585) { xsec = 0.000154719; xsec_unc = 0.26195; return; }
    else if(stop_mass == 1590) { xsec = 0.000150506; xsec_unc = 0.264111; return; }
    else if(stop_mass == 1595) { xsec = 0.000145626; xsec_unc = 0.263077; return; }
    else if(stop_mass == 1600) { xsec = 0.000141382; xsec_unc = 0.265291; return; }
    else if(stop_mass == 1605) { xsec = 0.000137131; xsec_unc = 0.267424; return; }
    else if(stop_mass == 1610) { xsec = 0.000132187; xsec_unc = 0.26668; return; }
    else if(stop_mass == 1615) { xsec = 0.000127929; xsec_unc = 0.269117; return; }
    else if(stop_mass == 1620) { xsec = 0.000124086; xsec_unc = 0.267738; return; }
    else if(stop_mass == 1625) { xsec = 0.00011982; xsec_unc = 0.270483; return; }
    else if(stop_mass == 1630) { xsec = 0.000116042; xsec_unc = 0.268071; return; }
    else if(stop_mass == 1635) { xsec = 0.000112767; xsec_unc = 0.27127; return; }
    else if(stop_mass == 1640) { xsec = 0.000108936; xsec_unc = 0.269351; return; }
    else if(stop_mass == 1645) { xsec = 0.000105746; xsec_unc = 0.271783; return; }
    else if(stop_mass == 1650) { xsec = 0.000102693; xsec_unc = 0.27292; return; }
    else if(stop_mass == 1655) { xsec = 0.000100112; xsec_unc = 0.274445; return; }
    else if(stop_mass == 1660) { xsec = 9.75763e-05; xsec_unc = 0.275431; return; }
    else if(stop_mass == 1665) { xsec = 9.52062e-05; xsec_unc = 0.276946; return; }
    else if(stop_mass == 1670) { xsec = 9.29857e-05; xsec_unc = 0.277869; return; }
    else if(stop_mass == 1675) { xsec = 9.08285e-05; xsec_unc = 0.279347; return; }
    else if(stop_mass == 1680) { xsec = 8.87433e-05; xsec_unc = 0.281539; return; }
    else if(stop_mass == 1685) { xsec = 8.66618e-05; xsec_unc = 0.283509; return; }
    else if(stop_mass == 1690) { xsec = 8.46535e-05; xsec_unc = 0.284432; return; }
    else if(stop_mass == 1695) { xsec = 8.27102e-05; xsec_unc = 0.28591; return; }
    else if(stop_mass == 1700) { xsec = 8.07774e-05; xsec_unc = 0.287497; return; }
    else if(stop_mass == 1705) { xsec = 7.8666e-05; xsec_unc = 0.288194; return; }
    else if(stop_mass == 1710) { xsec = 7.6572e-05; xsec_unc = 0.290265; return; }
    else if(stop_mass == 1715) { xsec = 7.45994e-05; xsec_unc = 0.291193; return; }
    else if(stop_mass == 1720) { xsec = 7.25199e-05; xsec_unc = 0.293013; return; }
    else if(stop_mass == 1725) { xsec = 7.05189e-05; xsec_unc = 0.293697; return; }
    else if(stop_mass == 1730) { xsec = 6.85712e-05; xsec_unc = 0.294972; return; }
    else if(stop_mass == 1735) { xsec = 6.67296e-05; xsec_unc = 0.296167; return; }
    else if(stop_mass == 1740) { xsec = 6.49184e-05; xsec_unc = 0.297686; return; }
    else if(stop_mass == 1745) { xsec = 6.30949e-05; xsec_unc = 0.298524; return; }
    else if(stop_mass == 1750) { xsec = 6.13637e-05; xsec_unc = 0.299789; return; }
    else if(stop_mass == 1755) { xsec = 5.97301e-05; xsec_unc = 0.300928; return; }
    else if(stop_mass == 1760) { xsec = 5.80751e-05; xsec_unc = 0.302585; return; }
    else if(stop_mass == 1765) { xsec = 5.65479e-05; xsec_unc = 0.30366; return; }
    else if(stop_mass == 1770) { xsec = 5.49998e-05; xsec_unc = 0.305241; return; }
    else if(stop_mass == 1775) { xsec = 5.35686e-05; xsec_unc = 0.306718; return; }
    else if(stop_mass == 1780) { xsec = 5.20828e-05; xsec_unc = 0.306799; return; }
    else if(stop_mass == 1785) { xsec = 5.07079e-05; xsec_unc = 0.309201; return; }
    else if(stop_mass == 1790) { xsec = 4.93948e-05; xsec_unc = 0.310043; return; }
    else if(stop_mass == 1795) { xsec = 4.80635e-05; xsec_unc = 0.31138; return; }
    else if(stop_mass == 1800) { xsec = 4.67492e-05; xsec_unc = 0.312291; return; }
    else if(stop_mass == 1805) { xsec = 4.55055e-05; xsec_unc = 0.314321; return; }
    else if(stop_mass == 1810) { xsec = 4.42835e-05; xsec_unc = 0.315499; return; }
    else if(stop_mass == 1815) { xsec = 4.30744e-05; xsec_unc = 0.316302; return; }
    else if(stop_mass == 1820) { xsec = 4.19954e-05; xsec_unc = 0.317151; return; }
    else if(stop_mass == 1825) { xsec = 4.08527e-05; xsec_unc = 0.319048; return; }
    else if(stop_mass == 1830) { xsec = 3.97561e-05; xsec_unc = 0.319718; return; }
    else if(stop_mass == 1835) { xsec = 3.87041e-05; xsec_unc = 0.322028; return; }
    else if(stop_mass == 1840) { xsec = 3.76008e-05; xsec_unc = 0.32268; return; }
    else if(stop_mass == 1845) { xsec = 3.66914e-05; xsec_unc = 0.324529; return; }
    else if(stop_mass == 1850) { xsec = 3.56995e-05; xsec_unc = 0.325039; return; }
    else if(stop_mass == 1855) { xsec = 3.47689e-05; xsec_unc = 0.326767; return; }
    else if(stop_mass == 1860) { xsec = 3.38528e-05; xsec_unc = 0.328878; return; }
    else if(stop_mass == 1865) { xsec = 3.29644e-05; xsec_unc = 0.328975; return; }
    else if(stop_mass == 1870) { xsec = 3.20679e-05; xsec_unc = 0.329608; return; }
    else if(stop_mass == 1875) { xsec = 3.12583e-05; xsec_unc = 0.331541; return; }
    else if(stop_mass == 1880) { xsec = 3.04342e-05; xsec_unc = 0.333117; return; }
    else if(stop_mass == 1885) { xsec = 2.96516e-05; xsec_unc = 0.332866; return; }
    else if(stop_mass == 1890) { xsec = 2.88952e-05; xsec_unc = 0.336279; return; }
    else if(stop_mass == 1895) { xsec = 2.81145e-05; xsec_unc = 0.336845; return; }
    else if(stop_mass == 1900) { xsec = 2.73974e-05; xsec_unc = 0.338247; return; }
    else if(stop_mass == 1905) { xsec = 2.66796e-05; xsec_unc = 0.339708; return; }
    else if(stop_mass == 1910) { xsec = 2.59941e-05; xsec_unc = 0.339526; return; }
    else if(stop_mass == 1915) { xsec = 2.52784e-05; xsec_unc = 0.341137; return; }
    else if(stop_mass == 1920) { xsec = 2.46598e-05; xsec_unc = 0.342714; return; }
    else if(stop_mass == 1925) { xsec = 2.39932e-05; xsec_unc = 0.342328; return; }
    else if(stop_mass == 1930) { xsec = 2.33737e-05; xsec_unc = 0.34394; return; }
    else if(stop_mass == 1935) { xsec = 2.27623e-05; xsec_unc = 0.345138; return; }
    else if(stop_mass == 1940) { xsec = 2.21454e-05; xsec_unc = 0.346933; return; }
    else if(stop_mass == 1945) { xsec = 2.15924e-05; xsec_unc = 0.350815; return; }
    else if(stop_mass == 1950) { xsec = 2.10232e-05; xsec_unc = 0.349444; return; }
    else if(stop_mass == 1955) { xsec = 2.05211e-05; xsec_unc = 0.350155; return; }
    else if(stop_mass == 1960) { xsec = 1.98996e-05; xsec_unc = 0.352135; return; }
    else if(stop_mass == 1965) { xsec = 1.9408e-05; xsec_unc = 0.353328; return; }
    else if(stop_mass == 1970) { xsec = 1.88974e-05; xsec_unc = 0.354643; return; }
    else if(stop_mass == 1975) { xsec = 1.84612e-05; xsec_unc = 0.357904; return; }
    else if(stop_mass == 1980) { xsec = 1.79562e-05; xsec_unc = 0.358898; return; }
    else if(stop_mass == 1985) { xsec = 1.75673e-05; xsec_unc = 0.35989; return; }
    else if(stop_mass == 1990) { xsec = 1.70612e-05; xsec_unc = 0.360953; return; }
    else if(stop_mass == 1995) { xsec = 1.66228e-05; xsec_unc = 0.364709; return; }
    else if(stop_mass == 2000) { xsec = 1.62355e-05; xsec_unc = 0.365277; return; }
    else { xsec = 0.; xsec_unc = 0.; }
  }
}
