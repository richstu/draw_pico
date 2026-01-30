// Joseph's code implementing on-the-fly energy error scale and smear for run 3

#include "core/correction.hpp"
#include "core/named_func.hpp"
#include "zgamma/scalesmear.hpp"

#include "TMath.h"

#include <cmath>
#include <memory>
#include <string>
#include <vector>

using std::string;
using std::vector;

std::unique_ptr<correction::CorrectionSet> cs_scale_syst_2022;
std::unique_ptr<correction::CorrectionSet> cs_scale_syst_2022EE;
std::unique_ptr<correction::CorrectionSet> cs_scale_syst_2023;
std::unique_ptr<correction::CorrectionSet> cs_scale_syst_2023BPix;

correction::CompoundCorrection::Ref map_scale_2022;
correction::Correction::Ref map_smearing_2022;
correction::CompoundCorrection::Ref map_scale_2022EE;
correction::Correction::Ref map_smearing_2022EE;
correction::CompoundCorrection::Ref map_scale_2023;
correction::Correction::Ref map_smearing_2023;
correction::CompoundCorrection::Ref map_scale_2023BPix;
correction::Correction::Ref map_smearing_2023BPix;

void cs_setter(){
cs_scale_syst_2022 = correction::CorrectionSet::from_file(
    "txt/corrections/2022/photonSS_EtDependent.json");
map_scale_2022 = cs_scale_syst_2022->compound().at(
    "EGMScale_Compound_Pho_2022preEE");
map_smearing_2022 = cs_scale_syst_2022->at(
    "EGMSmearAndSyst_PhoPTsplit_2022preEE");

cs_scale_syst_2022EE = correction::CorrectionSet::from_file(
    "txt/corrections/2022EE/photonSS_EtDependent.json");
map_scale_2022EE = cs_scale_syst_2022EE->compound().at(
    "EGMScale_Compound_Pho_2022postEE");
map_smearing_2022EE = cs_scale_syst_2022EE->at(
    "EGMSmearAndSyst_PhoPTsplit_2022postEE");

cs_scale_syst_2023 = correction::CorrectionSet::from_file(
    "txt/corrections/2023/photonSS_EtDependent.json");
map_scale_2023 = cs_scale_syst_2023->compound().at(
    "EGMScale_Compound_Pho_2023preBPIX");
map_smearing_2023 = cs_scale_syst_2023->at(
    "EGMSmearAndSyst_PhoPTsplit_2023preBPIX");

cs_scale_syst_2023BPix = correction::CorrectionSet::from_file(
    "txt/corrections/2023BPix/photonSS_EtDependent.json");
map_scale_2023BPix = cs_scale_syst_2023BPix->compound().at(
    "EGMScale_Compound_Pho_2023postBPIX");
map_smearing_2023BPix = cs_scale_syst_2023BPix->at(
    "EGMSmearAndSyst_PhoPTsplit_2023postBPIX");
}

NamedFunc photon_relpterr_corrected("photon_pterr_corrected",[](const Baby &b)->NamedFunc::ScalarType{
  float smear, scale, smearing;
  float energy_raw = b.photon_pt_raw()->at(0)*cosh(b.photon_eta()->at(0));
  float energyErr_corrected = -999.;
  float energyErr = b.photon_energyErr()->at(0);
  if(b.type()<1000 && b.type()>-1000){
    if (b.SampleTypeString().Contains("2022")) smear = map_smearing_2022->evaluate({"smear",b.photon_pt()->at(0),b.photon_r9()->at(0),fabs(b.photon_eta()->at(0))});
    else if (b.SampleTypeString().Contains("2022EE")) smear = map_smearing_2022EE->evaluate({"smear",b.photon_pt()->at(0),b.photon_r9()->at(0),fabs(b.photon_eta()->at(0))});
    else if (b.SampleTypeString().Contains("2023")) smear = map_smearing_2023->evaluate({"smear",b.photon_pt()->at(0),b.photon_r9()->at(0),fabs(b.photon_eta()->at(0))});
    else if (b.SampleTypeString().Contains("2023BPix")) smear = map_smearing_2023BPix->evaluate({"smear",b.photon_pt()->at(0),b.photon_r9()->at(0),fabs(b.photon_eta()->at(0))});
    else smear = 0;
    scale = b.photon_pt()->at(0)/b.photon_pt_raw()->at(0);
    energyErr_corrected = sqrt((energyErr*energyErr) + (energy_raw*energy_raw*smear*smear))*scale;
  } else if (b.type()>=1000 || b.type()<=-1000){
    if (b.SampleTypeString().Contains("2022")) smear = map_smearing_2022->evaluate({"smear",b.photon_pt_raw()->at(0),b.photon_r9()->at(0),fabs(b.photon_eta()->at(0))});
    else if (b.SampleTypeString().Contains("2022EE")) smear = map_smearing_2022EE->evaluate({"smear",b.photon_pt_raw()->at(0),b.photon_r9()->at(0),fabs(b.photon_eta()->at(0))});
    else if (b.SampleTypeString().Contains("2023")) smear = map_smearing_2023->evaluate({"smear",b.photon_pt_raw()->at(0),b.photon_r9()->at(0),fabs(b.photon_eta()->at(0))});
    else if (b.SampleTypeString().Contains("2023BPix")) smear = map_smearing_2023BPix->evaluate({"smear",b.photon_pt_raw()->at(0),b.photon_r9()->at(0),fabs(b.photon_eta()->at(0))});
    else smear = 0;
    smearing = b.photon_pt()->at(0)/b.photon_pt_raw()->at(0);//To avoid having to regenerate random numbers
    energyErr_corrected = sqrt((energyErr*energyErr) + (energy_raw*energy_raw*smear*smear))*smearing;
  }
  return energyErr_corrected;
});

//get photon rel pt error with variation
NamedFunc assign_variation_lead_photon_relpterr_corrected(
    const NamedFunc &ph_sig, const NamedFunc &ph_pt, const string &name) {
  return NamedFunc(("sys_lead_photon_relpterr_"+name).c_str(),
      [&ph_sig, &ph_pt](const Baby &b) -> NamedFunc::ScalarType{
    vector<double> photon_pt = ph_pt.GetVector(b);
    vector<double> photon_sig = ph_sig.GetVector(b);
    double lead_photon_eta(0.0), lead_photon_pt(-999.0);
    double lead_photon_energyErr(0.0), lead_photon_pt_raw(-999.0);
    double lead_photon_r9(0.0);
    for (unsigned iph = 0; iph < photon_sig.size(); iph++) {
      if (photon_sig[iph]) {
        if (lead_photon_pt < photon_pt[iph]) {
          lead_photon_pt = photon_pt[iph];
          lead_photon_r9 = b.photon_r9()->at(iph);
          lead_photon_eta = b.photon_eta()->at(iph);
          lead_photon_pt_raw = b.photon_pt_raw()->at(iph);
          lead_photon_energyErr = b.photon_energyErr()->at(iph);
        }
      }
    }
    double smear(0.0), scale(1.0), smearing(1.0);
    double energy_raw = lead_photon_pt_raw*TMath::CosH(lead_photon_eta);
    double energyErr_corrected = lead_photon_energyErr;
    if (b.SampleType() < 2022 && b.SampleType() > -2022)
      return lead_photon_energyErr/(lead_photon_pt
                                    *TMath::CosH(lead_photon_eta));
    if(b.type()<1000 && b.type()>-1000){
      if (b.SampleTypeString().Contains("2022")) 
        smear = map_smearing_2022->evaluate({"smear",lead_photon_pt,
            lead_photon_r9,fabs(lead_photon_eta)});
      else if (b.SampleTypeString().Contains("2022EE")) 
        smear = map_smearing_2022EE->evaluate({"smear",lead_photon_pt,
            lead_photon_r9,fabs(lead_photon_eta)});
      else if (b.SampleTypeString().Contains("2023")) 
        smear = map_smearing_2023->evaluate({"smear",lead_photon_pt,
            lead_photon_r9,fabs(lead_photon_eta)});
      else if (b.SampleTypeString().Contains("2023BPix")) 
        smear = map_smearing_2023BPix->evaluate({"smear",lead_photon_pt,
            lead_photon_r9,fabs(lead_photon_eta)});
      scale = lead_photon_pt/lead_photon_pt_raw;
      energyErr_corrected = sqrt((lead_photon_energyErr*lead_photon_energyErr) 
          + (energy_raw*energy_raw*smear*smear))*scale;
    } else if (b.type()>=1000 || b.type()<=-1000){
      if (b.SampleTypeString().Contains("2022")) 
        smear = map_smearing_2022->evaluate({"smear",lead_photon_pt_raw,
            lead_photon_r9,fabs(lead_photon_eta)});
      else if (b.SampleTypeString().Contains("2022EE")) 
        smear = map_smearing_2022EE->evaluate({"smear",lead_photon_pt_raw,
            lead_photon_r9,fabs(lead_photon_eta)});
      else if (b.SampleTypeString().Contains("2023")) 
        smear = map_smearing_2023->evaluate({"smear",lead_photon_pt_raw,
            lead_photon_r9,fabs(lead_photon_eta)});
      else if (b.SampleTypeString().Contains("2023BPix")) 
        smear = map_smearing_2023BPix->evaluate({"smear",lead_photon_pt_raw,
            lead_photon_r9,fabs(lead_photon_eta)});
      // To avoid having to regenerate random numbers
      smearing = lead_photon_pt/lead_photon_pt_raw;
      energyErr_corrected = sqrt((lead_photon_energyErr*lead_photon_energyErr) 
          + (energy_raw*energy_raw*smear*smear))*smearing;
    }
    return energyErr_corrected/(lead_photon_pt*TMath::CosH(lead_photon_eta));
  }).EnableCaching(true);
}
