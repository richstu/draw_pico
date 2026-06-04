#include "zgamma/eventweighter.hpp"

#include "core/baby.hpp"
#include "core/correction.hpp"
#include "core/named_func.hpp"
#include "core/named_func_utilities.hpp"
#include "core/utilities.hpp"

#include <cmath>
#include <iostream>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

using std::cout;
using std::endl;
using std::lock_guard;
using std::make_unique;
using std::max;
using std::mutex;
using std::string;
using std::shared_ptr;
using std::unique_ptr;
using std::unordered_map;
using std::vector;

using NamedFuncUtilities::ReduceNamedFuncCached;

double reduce_index0(vector<double> inputs) {
  return inputs[0];
}

double reduce_index1(vector<double> inputs) {
  return inputs[1];
}

double reduce_index2(vector<double> inputs) {
  return inputs[2];
}

ElectronEventWeighter::ElectronEventWeighter(std::string year) {
  post_bpix_                = false;
  if (year=="2016APV") {
    in_file_electron_         = "txt/data/zgamma/2016APV/hzg_elid_2016APV_scalefactors.json";
    in_file_electron_reco_    = "txt/data/zgamma/2016APV/electron_recoSF2016preVFP.json";
    in_file_electron_iso0p10_ = "txt/data/zgamma/2016APV/hzg_eliso0p1_2016APV_efficiencies.json";
    in_file_electron_iso0p15_ = "txt/data/zgamma/2016APV/hzg_eliso0p15_2016APV_efficiencies.json";
  } else if (year=="2016") {
    in_file_electron_         = "txt/data/zgamma/2016/hzg_elid_2016_scalefactors.json";
    in_file_electron_reco_    = "txt/data/zgamma/2016/electron_recoSF2016postVFP.json";
    in_file_electron_iso0p10_ = "txt/data/zgamma/2016/hzg_eliso0p1_2016_efficiencies.json";
    in_file_electron_iso0p15_ = "txt/data/zgamma/2016/hzg_eliso0p15_2016_efficiencies.json";
  } else if (year=="2017") {
    in_file_electron_         = "txt/data/zgamma/2017/hzg_elid_2017_scalefactors.json";
    in_file_electron_reco_    = "txt/data/zgamma/2017/electron_recoSF2017.json";
    in_file_electron_iso0p10_ = "txt/data/zgamma/2017/hzg_eliso0p1_2017_efficiencies.json";
    in_file_electron_iso0p15_ = "txt/data/zgamma/2017/hzg_eliso0p15_2017_efficiencies.json";
  } else if (year=="2018") {
    in_file_electron_         = "txt/data/zgamma/2018/hzg_elid_2018_scalefactors.json";
    in_file_electron_reco_    = "txt/data/zgamma/2018/electron_recoSF2018.json";
    in_file_electron_iso0p10_ = "txt/data/zgamma/2018/hzg_eliso0p1_2018_efficiencies.json";
    in_file_electron_iso0p15_ = "txt/data/zgamma/2018/hzg_eliso0p15_2018_efficiencies.json";
  } else if (year=="2022"){
    in_file_electron_         = "txt/data/zgamma/2022/hzg_elid_2022_scalefactors.json";
    in_file_electron_reco_    = "txt/data/zgamma/2022/electron_recoSF2022.json";
    in_file_electron_iso0p10_ = "txt/data/zgamma/2022/hzg_eliso0p1_2022_efficiencies.json";
    in_file_electron_iso0p15_ = "txt/data/zgamma/2022/hzg_eliso0p15_2022_efficiencies.json";
  } else if (year=="2022EE"){
    in_file_electron_         = "txt/data/zgamma/2022EE/hzg_elid_2022EE_scalefactors.json";
    in_file_electron_reco_    = "txt/data/zgamma/2022EE/electron_recoSF2022EE.json";
    in_file_electron_iso0p10_ = "txt/data/zgamma/2022EE/hzg_eliso0p1_2022EE_efficiencies.json";
    in_file_electron_iso0p15_ = "txt/data/zgamma/2022EE/hzg_eliso0p15_2022EE_efficiencies.json";
  } else if (year=="2023"){
    in_file_electron_         = "txt/data/zgamma/2023/hzg_elid_2023_scalefactors.json";
    in_file_electron_reco_    = "txt/data/zgamma/2023/electron_recoSF2023.json";
    in_file_electron_iso0p10_ = "txt/data/zgamma/2023/hzg_eliso0p1_2023_efficiencies.json";
    in_file_electron_iso0p15_ = "txt/data/zgamma/2023/hzg_eliso0p15_2023_efficiencies.json";
  } else if (year=="2023BPix"){
    in_file_electron_         = "txt/data/zgamma/2023BPix/hzg_elid_2023BPix_scalefactors.json";
    in_file_electron_reco_    = "txt/data/zgamma/2023BPix/electron_recoSF2023BPix.json";
    in_file_electron_iso0p10_ = "txt/data/zgamma/2023BPix/hzg_eliso0p1_2023BPix_efficiencies.json";
    in_file_electron_iso0p15_ = "txt/data/zgamma/2023BPix/hzg_eliso0p15_2023BPix_efficiencies.json";
    cs_electron_bpixhole_     = correction::CorrectionSet::from_file(
        "txt/data/zgamma/2023BPix/hzg_elid_2023BPixHole_scalefactors.json");
    cs_el_hole_iso0p10_       = correction::CorrectionSet::from_file(
        "txt/data/zgamma/2023BPix/hzg_eliso0p1_2023BPixHole_efficiencies.json");
    cs_el_hole_iso0p15_       = correction::CorrectionSet::from_file(
        "txt/data/zgamma/2023BPix/hzg_eliso0p15_2023BPixHole_efficiencies.json");
    post_bpix_                     = true;
  } else {
    cout<<"Year has not been implemented in event_weighter"<<endl;
  }
  cs_electron_              = correction::CorrectionSet::from_file(in_file_electron_);
  cs_electron_reco_         = correction::CorrectionSet::from_file(in_file_electron_reco_);
  cs_el_iso0p10_            = correction::CorrectionSet::from_file(in_file_electron_iso0p10_);
  cs_el_iso0p15_            = correction::CorrectionSet::from_file(in_file_electron_iso0p15_);
  map_electron_id_pass_     = cs_electron_->at("sf_pass");
  map_electron_id_pass_unc_ = cs_electron_->at("unc_pass");
  map_electron_id_fail_     = cs_electron_->at("sf_fail");
  map_electron_id_fail_unc_ = cs_electron_->at("unc_fail");
  map_electron_reco_pass_     = cs_electron_reco_->at("sf_pass");
  map_electron_reco_pass_unc_ = cs_electron_reco_->at("unc_pass");
  map_electron_reco_fail_     = cs_electron_reco_->at("sf_fail");
  map_electron_reco_fail_unc_ = cs_electron_reco_->at("unc_fail");
  if (year == "2023BPix") {
    map_electron_hole_id_pass_     = cs_electron_bpixhole_->at("sf_pass");
    map_electron_hole_id_pass_unc_ = cs_electron_bpixhole_->at("unc_pass");
    map_electron_hole_id_fail_     = cs_electron_bpixhole_->at("sf_fail");
    map_electron_hole_id_fail_unc_ = cs_electron_bpixhole_->at("unc_fail");
  }
}

vector<double> ElectronEventWeighter::ElectronSF(const Baby &b) {
  double sf_tot = 1.0;
  double sf_tot_up = 1.0;
  double sf_tot_dn = 1.0;
  //SF logic: for each true electron, weight by prob(MC)/prob(data) which is
  //prob(reco and pass id) or prob(fail)=1-prob(reco and pass id)
  for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
    if (abs(b.mc_id()->at(imc))==11 && 
        ((b.mc_statusflag()->at(imc) & 0x2000)!=0) &&
        ((b.mc_statusflag()->at(imc) & 0x1) != 0)) {
      //is electron and last copy and prompt
      if ((b.mc_pt()->at(imc)<7.0f) 
          || (fabs(b.mc_eta()->at(imc))>2.5f)) continue;
      bool pass_id = false;
      float reco_pt = -999;
      float reco_eta = -999;
      float reco_phi = -999;
      float min_dr = 999;
      for(unsigned iel = 0; iel < b.el_sig()->size(); ++iel){
        if (b.el_sig()->at(iel)) {
          float dr = deltaR(b.mc_eta()->at(imc),b.mc_phi()->at(imc),
                            b.el_eta()->at(iel),b.el_phi()->at(iel));
          if (dr < 0.4f && dr < min_dr) {
            reco_pt = b.el_pt()->at(iel);
            reco_eta = b.el_eta()->at(iel);
            reco_phi = b.el_phi()->at(iel);
            min_dr = dr;
            pass_id = true;
          }
        }
      }
      if (reco_pt < 0) {
        reco_pt = b.mc_pt()->at(imc);
        reco_eta = b.mc_eta()->at(imc);
        reco_phi = b.mc_phi()->at(imc);
      }
      float sf_reco = 1.0;
      float unc_reco = 1.0;
      float sf = 1.0;
      float unc = 1.0;
      float sf_up = 1.0;
      float sf_dn = 1.0;
      bool in_bpix_region = (reco_eta > -1.566 && reco_eta < 0.0
                             && reco_phi > -1.2 && reco_phi < -0.8);
      if (pass_id) {
        if (post_bpix_ && in_bpix_region) {
          sf = map_electron_hole_id_pass_->evaluate({reco_pt,reco_eta});
          unc = map_electron_hole_id_pass_unc_->evaluate({reco_pt,reco_eta});
          sf_reco = map_electron_reco_pass_->evaluate({reco_pt,reco_eta});
          unc_reco = map_electron_reco_pass_unc_->evaluate({reco_pt,reco_eta});
        }
        else {
          sf = map_electron_id_pass_->evaluate({reco_pt,reco_eta});
          unc = map_electron_id_pass_unc_->evaluate({reco_pt,reco_eta});
          sf_reco = map_electron_reco_pass_->evaluate({reco_pt,reco_eta});
          unc_reco = map_electron_reco_pass_unc_->evaluate({reco_pt,reco_eta});
        }
      }
      else {
        if (post_bpix_ && in_bpix_region) {
          sf = map_electron_hole_id_pass_->evaluate({reco_pt,reco_eta});
          unc = -1.0*map_electron_hole_id_pass_unc_->evaluate({reco_pt,
                                                               reco_eta});
          sf_reco = map_electron_reco_fail_->evaluate({reco_pt,reco_eta});
          unc_reco = -1.0*map_electron_reco_fail_unc_->evaluate({reco_pt,
                                                                 reco_eta});
        }
        else {
          sf = map_electron_id_fail_->evaluate({reco_pt,reco_eta});
          unc = -1.0*map_electron_id_fail_unc_->evaluate({reco_pt,reco_eta});
        }
      }
      sf_up = (sf+unc)*(sf_reco+unc_reco);
      sf_dn = (sf-unc)*(sf_reco-unc_reco);
      if (isinf(sf) || isnan(sf)) sf = 1.0;
      if (isinf(sf_reco) || isnan(sf_reco)) sf_reco = 1.0;
      if (isinf(sf_up) || isnan(sf_up)) sf_up = 1.0;
      if (isinf(sf_dn) || isnan(sf_dn)) sf_dn = 1.0;
      sf_tot *= sf*sf_reco;
      sf_tot_up *= sf_up;
      sf_tot_dn *= fmax(sf_dn,0.0);
    }
  }
  //return is {nominal, up, dn}
  return {sf_tot, sf_tot_up, sf_tot_dn};
}

vector<double> ElectronEventWeighter::ElectronMinisoSF(const Baby &b) {
  double sf_tot = 1.0;
  double sf_tot_up = 1.0;
  double sf_tot_dn = 1.0;
  for (unsigned iel = 0; iel < b.el_pt()->size(); iel++) {
    if (!b.el_sig()->at(iel)) continue;
    float pt = b.el_pt()->at(iel);
    float eta = b.el_eta()->at(iel);
    float phi = b.el_phi()->at(iel);
    float miniso = b.el_miniso()->at(iel);
    float sf = 1.0;
    float sf_up = 1.0;
    float sf_dn = 1.0;
    bool in_bpix_region = (eta > -1.566 && eta < 0.0
                           && phi > -1.2 && phi < -0.8);
    float data_eff_0p10 = 1.0;
    float simu_eff_0p10 = 1.0;
    float data_eff_0p15 = 1.0;
    float simu_eff_0p15 = 1.0;
    float data_unc_0p10 = 1.0;
    float simu_unc_0p10 = 1.0;
    float data_unc_0p15 = 1.0;
    float simu_unc_0p15 = 1.0;
    if (post_bpix_ && in_bpix_region) {
      data_eff_0p10 = cs_el_hole_iso0p10_->at("effdata")->evaluate({pt,eta});
      simu_eff_0p10 = cs_el_hole_iso0p10_->at("effmc")->evaluate({pt,eta});
      data_eff_0p15 = cs_el_hole_iso0p15_->at("effdata")->evaluate({pt,eta});
      simu_eff_0p15 = cs_el_hole_iso0p15_->at("effmc")->evaluate({pt,eta});
      data_unc_0p10 = cs_el_hole_iso0p10_->at("systdata")->evaluate({pt,eta});
      simu_unc_0p10 = cs_el_hole_iso0p10_->at("systmc")->evaluate({pt,eta});
      data_unc_0p15 = cs_el_hole_iso0p15_->at("systdata")->evaluate({pt,eta});
      simu_unc_0p15 = cs_el_hole_iso0p15_->at("systmc")->evaluate({pt,eta});
    }
    else {
      data_eff_0p10 = cs_el_iso0p10_->at("effdata")->evaluate({pt,eta});
      simu_eff_0p10 = cs_el_iso0p10_->at("effmc")->evaluate({pt,eta});
      data_eff_0p15 = cs_el_iso0p15_->at("effdata")->evaluate({pt,eta});
      simu_eff_0p15 = cs_el_iso0p15_->at("effmc")->evaluate({pt,eta});
      data_unc_0p10 = cs_el_iso0p10_->at("systdata")->evaluate({pt,eta});
      simu_unc_0p10 = cs_el_iso0p10_->at("systmc")->evaluate({pt,eta});
      data_unc_0p15 = cs_el_iso0p15_->at("systdata")->evaluate({pt,eta});
      simu_unc_0p15 = cs_el_iso0p15_->at("systmc")->evaluate({pt,eta});
    }
    float data_eff = 1.0;
    float data_eff_up = 1.0;
    float data_eff_dn = 1.0;
    float simu_eff = 1.0;
    float simu_eff_up = 1.0;
    float simu_eff_dn = 1.0;
    if (miniso < 0.1) {
      data_eff = data_eff_0p10;
      simu_eff = simu_eff_0p10;
      data_eff_up = (data_eff_0p10+data_unc_0p10);
      data_eff_dn = (data_eff_0p10-data_unc_0p10);
      simu_eff_up = (simu_eff_0p10+simu_unc_0p10);
      simu_eff_dn = (simu_eff_0p10-simu_unc_0p10);
    }
    else if (miniso < 0.15) {
      data_eff = (data_eff_0p15-data_eff_0p10);
      data_eff_up = ((data_eff_0p15+data_unc_0p15)
                    -(data_eff_0p10+data_unc_0p10));
      data_eff_dn = ((data_eff_0p15-data_unc_0p15)
                    -(data_eff_0p10-data_unc_0p10));
      simu_eff = (simu_eff_0p15-simu_eff_0p10);
      simu_eff_up = ((simu_eff_0p15+simu_unc_0p15)
                    -(simu_eff_0p10+simu_unc_0p10));
      simu_eff_dn = ((simu_eff_0p15-simu_unc_0p15)
                    -(simu_eff_0p10-simu_unc_0p10));
    }
    else {
      data_eff = (1.0-data_eff_0p15);
      data_eff_up = (1.0-(data_eff_0p15+data_unc_0p15));
      data_eff_dn = (1.0-(data_eff_0p15-data_unc_0p15));
      simu_eff = (1.0-simu_eff_0p15);
      simu_eff_up = (1.0-(simu_eff_0p15+simu_unc_0p15));
      simu_eff_dn = (1.0-(simu_eff_0p15-simu_unc_0p15));
    }
    data_eff = max(data_eff, 0.0f);
    data_eff_up = max(data_eff_up, 0.0f);
    data_eff_dn = max(data_eff_dn, 0.0f);
    simu_eff = max(simu_eff, 0.0f);
    simu_eff_up = max(simu_eff_up, 0.0f);
    simu_eff_dn = max(simu_eff_dn, 0.0f);
    sf = data_eff/simu_eff;
    sf_up = data_eff_up/simu_eff_dn;
    sf_dn = data_eff_dn/simu_eff_up;
    if (isinf(sf) || isnan(sf)) sf = 1.0;
    if (isinf(sf_up) || isnan(sf_up)) sf_up = 1.0;
    if (isinf(sf_dn) || isnan(sf_dn)) sf_dn = 1.0;
    sf_tot *= sf;
    sf_tot_up *= sf_up;
    sf_tot_dn *= sf_dn;
  }
  return {sf_tot, sf_tot_up, sf_tot_dn};
}

vector<double> ElectronEventWeighter::ElectronCombinedSF(const Baby &b) {
  vector<double> w_elid = ElectronSF(b);
  vector<double> w_eliso = ElectronMinisoSF(b);
  return {w_elid[0]*w_eliso[0]/b.w_el(), 
          w_elid[1]*w_eliso[1]/w_elid[0]/w_eliso[0], 
          w_elid[2]*w_eliso[2]/w_elid[0]/w_eliso[0]};
}

const NamedFunc sys_el = NamedFunc("sys_el",[](const Baby &b) 
    -> NamedFunc::VectorType{
  static unordered_map<string, unique_ptr<ElectronEventWeighter>> weighters;
  static mutex weighter_mutex;
  vector<double> ret;
  //if (b.SampleTypeString().Contains("-")) 
  if (b.SampleType() < 0.0) 
    return {1.0f, 1.0f, 1.0f};
  {
    lock_guard<mutex> lock(weighter_mutex);
    string year_string(b.SampleTypeString().Data());
    if (weighters.count(year_string) == 0) {
      weighters[year_string] = make_unique<ElectronEventWeighter>(year_string);
    }
    return weighters[year_string]->ElectronCombinedSF(b);
  }
}).EnableCaching(true);

const NamedFunc w_el = ReduceNamedFuncCached(sys_el, 
    reduce_index0).Name("w_el");

const NamedFunc sys_el_up = ReduceNamedFuncCached(sys_el, 
    reduce_index1).Name("sys_el_up");

const NamedFunc sys_el_dn = ReduceNamedFuncCached(sys_el, 
    reduce_index2).Name("sys_el_dn");

float get_w_fakephoton(string run, int ph_isjet, float ph_pt, float ph_eta) {
  static unique_ptr<correction::CorrectionSet> cs_fakephoton_;
  static shared_ptr<const correction::Correction> map_fakephoton_;
  static mutex fakephoton_mutex;
  float ret = 0.0;
  {
    lock_guard<mutex> lock(fakephoton_mutex);
    if (cs_fakephoton_ == nullptr) {
      cs_fakephoton_ = correction::CorrectionSet::from_file("txt/data/zgamma/fakephoton.json");
      map_fakephoton_ = cs_fakephoton_->at("fakephoton_corrections");
    }
    ret = map_fakephoton_->evaluate({run, ph_isjet, ph_pt, fabs(ph_eta)});
  }
  return ret;
}
