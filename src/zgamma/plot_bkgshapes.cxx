//This script generates some tables for run3 H->ZGamma trigger studies
//
//Arguments
// --single_thread (-s) to run single thread for debugging
// --unblind (-u) to unblind (not done for AN plots)
// --year (-y) yearname to select data year (2016, 2017, 2018, run2)
// --tag (-t) to add a tag to produced plots
// --string_options (-o) to specify what to plot
//
//possible string options: 
//

#include "core/test.hpp"

#include <algorithm>
#include <bitset>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>
#include <getopt.h>

#include "Math/Vector4D.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TError.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TPad.h"
#include "TPie.h"
#include "TVector2.h"

#include "core/baby.hpp"
#include "core/cross_sections.hpp"
#include "core/efficiency_plot.hpp"
#include "core/event_scan.hpp"
#include "core/functions.hpp"
#include "core/gamma_params.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/named_func.hpp"
#include "core/palette.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/process.hpp"
#include "core/table.hpp"
#include "core/utilities.hpp"

#include "zgamma/zg_functions.hpp"
#include "zgamma/zg_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;

//manual bring this in for now
namespace script_utilities { 
  struct ArgStruct {
    bool single_thread;
    std::string year_string;
    std::string string_options;
    std::string tag;
    bool unblind;
  };

  ArgStruct get_options(int argc, char* argv[], std::string default_string_options) {
    //
    //initialize to defaults
    ArgStruct options;
    options.single_thread = false;
    options.year_string = "run2";
    options.string_options = default_string_options;
    options.tag = "";
    options.unblind = false;

    //parse arguments
    while(true){
      static struct option long_options[] = {
        {"single_thread", no_argument, 0, 's'},
        {"year", required_argument, 0, 'y'},
        {"unblind", no_argument, 0, 'u'},
        {"tag", required_argument, 0, 't'},
        {"string_options", required_argument, 0, 'o'},
        {0, 0, 0, 0}
      };

      char opt = -1;
      int option_index;
      opt = getopt_long(argc, argv, "sy:ut:o:", long_options, &option_index);
      if (opt == -1) break;

      std::string optname;
      switch(opt){
      case 's':
        options.single_thread = true;
        break;
      case 'u':
        options.unblind = true;
        break;
      case 'y':
        options.year_string = optarg;
        break;
      case 't':
        options.tag = optarg;
        break;
      case 'o':
        options.string_options = optarg;
        break;
      case 0:
        //handle cases with no short argument form
        std::cout << "Bad option! Found option name " << optname << std::endl;
        break;
      default:
        //std::cout << "Bad option! getopt_long returned character code 0" << opt << std::endl;
        printf("Bad option! getopt_long returned character code 0%o\n", opt);
        break;
      }
    }
    return options;
  }

  std::vector<PlotOpt> plot_lin(bool unblind) {
    PlotOpt lin_norm("txt/plot_styles.txt", "CMSPaper");
    lin_norm.Title(PlotOptTypes::TitleType::info)   
        .Bottom(PlotOptTypes::BottomType::off)
        .YAxis(PlotOptTypes::YAxisType::linear)
        .Stack(PlotOptTypes::StackType::signal_overlay)
        .Overflow(PlotOptTypes::OverflowType::none)
        .LegendColumns(3);
    if (unblind) {
      lin_norm.Stack(PlotOptTypes::StackType::data_norm);
      lin_norm.Bottom(PlotOptTypes::BottomType::ratio);
    }
    return {lin_norm};
  }

  std::vector<PlotOpt> plot_log(bool unblind) {
    PlotOpt log_norm("txt/plot_styles.txt", "CMSPaper");
    log_norm.Title(PlotOptTypes::TitleType::info)   
        .Bottom(PlotOptTypes::BottomType::off)
        .YAxis(PlotOptTypes::YAxisType::log)
        .LogMinimum(.2)
        .Stack(PlotOptTypes::StackType::signal_overlay)
        .Overflow(PlotOptTypes::OverflowType::none)
        .LegendColumns(3);
    if (unblind) {
      log_norm.Stack(PlotOptTypes::StackType::data_norm);
      log_norm.Bottom(PlotOptTypes::BottomType::ratio);
    }
    return {log_norm};
  }

  std::vector<PlotOpt> plot_shapes() {
    PlotOpt shapes_norm("txt/plot_styles.txt", "CMSPaper");
    shapes_norm.Title(PlotOptTypes::TitleType::info)   
        .Bottom(PlotOptTypes::BottomType::off)
        .YAxis(PlotOptTypes::YAxisType::linear)
        .Stack(PlotOptTypes::StackType::shapes)
        .LegendColumns(3);
    return {shapes_norm};
  }

  std::vector<PlotOpt> plot_log_shapes() {
    PlotOpt shapes_norm("txt/plot_styles.txt", "CMSPaper");
    shapes_norm.Title(PlotOptTypes::TitleType::info)   
        .Bottom(PlotOptTypes::BottomType::off)
        .YAxis(PlotOptTypes::YAxisType::log)
        .LogMinimum(.2)
        .Stack(PlotOptTypes::StackType::shapes)
        .LegendColumns(3);
    return {shapes_norm};
  }
}

//count number of times a value appears in a vector
template <class vec_type>
unsigned vector_count(std::vector<vec_type> vec, vec_type value) {
  unsigned counter = 0;
  for (unsigned i = 0; i < vec.size(); i++) {
    if (vec[i]==value)
      counter++;
  }
  return counter;
}

//Get mass from pdg id - only works for SM particles
double GetMass(int pdgid) {
  switch (pdgid) {
    case 1: //fall through
    case -1:
      return 0.00467;
      break;
    case 2: //fall through
    case -2:
      return 0.00216; 
      break;
    case 3: //fall through
    case -3:
      return 0.0934;
      break;
    case 4: //fall through
    case -4:
      return 1.27;
      break;
    case 5: //fall through
    case -5:
      return 4.18;
      break;
    case 6: //fall through
    case -6:
      return 172.69;
      break;
    case 11: //fall through
    case -11:
      return 0.000511;
      break;
    case 13: //fall through
    case -13:
      return 0.106;
      break;
    case 15: //fall through
    case -15:
      return 1.777;
      break;
    case 23:
      return 91.188;
      break;
    case 24: //fall through
    case -24:
      return 80.377;
      break;
    case 25:
      return 125.25;
      break;
    default:
      return 0;
  }
  return 0;
}

//get pdgid of mother ignoring FSR
int DistinctMomPdgId(std::vector<int>* GenPart_pdgId, std::vector<int>* GenPart_genPartIdxMother, int idx) {
  int self_pdgid = GenPart_pdgId->at(idx);
  if (GenPart_genPartIdxMother->at(idx) != -1) {
    idx = GenPart_genPartIdxMother->at(idx);
  }
  bool ancestors_end = false;
  while (!ancestors_end) {
    if (GenPart_pdgId->at(idx) != self_pdgid) {
      return GenPart_pdgId->at(idx);
    }
    else if (GenPart_genPartIdxMother->at(idx) != -1) {
      idx = GenPart_genPartIdxMother->at(idx);
    }
    else {
      ancestors_end = true;
    }
  }
  return 0;
}

//get idx of mother ignoring FSR
int DistinctMomIdx(std::vector<int>* GenPart_pdgId, std::vector<int>* GenPart_genPartIdxMother, int idx) {
  int self_pdgid = GenPart_pdgId->at(idx);
  bool ancestors_end = false;
  while (!ancestors_end) {
    if (GenPart_pdgId->at(idx) != self_pdgid) {
      return idx;
    }
    else if (GenPart_genPartIdxMother->at(idx) != -1) {
      idx = GenPart_genPartIdxMother->at(idx);
    }
    else {
      ancestors_end = true;
    }
  }
  return 0;
}

//------------------------------------------------------------------------------------
//                      named funcs to be moved to a common location
//------------------------------------------------------------------------------------

//Returns a scalar NamedFunc that is the sum of the entries of vector_named_func
NamedFunc SumNamedFunc(NamedFunc vector_named_func) {
  return NamedFunc("SumNamedFunc("+vector_named_func.Name()+")",[vector_named_func](const Baby &b) -> NamedFunc::ScalarType{
    double sum = 0;
    std::vector<double> vector_named_func_ = vector_named_func.GetVector(b);
    for (double named_func_entry : vector_named_func_) {
      sum += named_func_entry;
    }
    return sum;
  });
}

////Returns a scalar NamedFunc that is element index of vector_named_func or -999 if no index index
//NamedFunc AtNamedFunc(NamedFunc vector_named_func, unsigned index) {
//  return NamedFunc("AtNamedFunc("+vector_named_func.Name()+","+std::to_string(index)+")",[vector_named_func](const Baby &b) -> NamedFunc::ScalarType{
//    std::vector<double> vector_named_func_ = vector_named_func.GetVector(b);
//    if (index < 0 || index >= vector_named_func_.size()) {
//      return -999;
//    }
//    return vector_named_func_[index];
//  });
//}

//Returns a vector named func that is vector_named_func filtered with filter_named_func
NamedFunc FilterNamedFunc(NamedFunc vector_named_func, NamedFunc filter_named_func) {
  return NamedFunc("FilterNamedFunc("+vector_named_func.Name()+","+filter_named_func.Name()+")",[vector_named_func,filter_named_func](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> mapped_named_func;
    std::vector<double> vector_named_func_ = vector_named_func.GetVector(b);
    std::vector<double> filter_named_func_ = filter_named_func.GetVector(b);
    for (unsigned i = 0; i < vector_named_func_.size(); i++) {
      if (filter_named_func_[i]) {
        mapped_named_func.push_back(vector_named_func_[i]);
      }
    }
    return mapped_named_func;
  });
}

std::function<double(std::vector<double>)> FnMax = [](std::vector<double> a) {
  double ret = -999;
  for (unsigned i = 0; i < a.size(); i++) {
    if (a[i]>ret)
      ret = a[i];
  } 
  return ret;
};

//Returns a vector named func that is vector_named_func with map_function applied to it
NamedFunc MapNamedFunc(NamedFunc vector_named_func, std::function<double(double)> map_function) {
  return NamedFunc("MapNamedFunc("+vector_named_func.Name()+")",[vector_named_func,map_function](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> mapped_named_func;
    std::vector<double> vector_named_func_ = vector_named_func.GetVector(b);
    for (unsigned i = 0; i < vector_named_func_.size(); i++) {
      mapped_named_func.push_back(map_function(vector_named_func_[i]));
    }
    return mapped_named_func;
  });
}

//Returns a scalar named func that is vector_named_func with reduce_function applied to it
NamedFunc ReduceNamedFunc(NamedFunc vector_named_func, std::function<double(std::vector<double>)> reduce_function) {
  return NamedFunc("ReduceNamedFunc("+vector_named_func.Name()+")",[vector_named_func,reduce_function](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> vector_named_func_ = vector_named_func.GetVector(b);
    return reduce_function(vector_named_func_);
  });
}

//Returns a scalar NamedFunc that is element index of vector_named_func or -999 if no index index
//only considering indices for which mask_named_func is true
NamedFunc AtNamedFunc(NamedFunc vector_named_func, NamedFunc mask_named_func, unsigned index) {
  return NamedFunc("AtNamedFunc("+vector_named_func.Name()+","+mask_named_func.Name()+","+std::to_string(index)+")",[vector_named_func, mask_named_func, index](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> vector_named_func_ = vector_named_func.GetVector(b);
    std::vector<double> mask_named_func_ = mask_named_func.GetVector(b);
    unsigned sig_index = 0;
    for (unsigned i = 0; i < vector_named_func_.size(); i++) {
      if (mask_named_func_[i]) {
        if (sig_index==index) {
          return vector_named_func_[i];
        }
        sig_index++;
      }
    }
    return -999;
  });
}

//returns a vector NamedFunc that is the list of indices from list2 matched by deltaR to list1
//can exclude elements of list2 using NamedFunc list2_mask
NamedFunc DeltaRMatchNamedFunc(NamedFunc list1_eta, NamedFunc list1_phi, NamedFunc list2_eta, NamedFunc list2_phi, NamedFunc list2_mask) {
  return NamedFunc("DeltaRMatchNamedFunc("+list1_eta.Name().substr(0,list1_eta.Name().size()-4)+")",[list1_eta, list1_phi, list2_eta, list2_phi, list2_mask](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> list2_indices;
    std::vector<double> list1_eta_ = list1_eta.GetVector(b);
    std::vector<double> list1_phi_ = list1_phi.GetVector(b);
    std::vector<double> list2_eta_ = list2_eta.GetVector(b);
    std::vector<double> list2_phi_ = list2_phi.GetVector(b);
    std::vector<double> list2_mask_ = list2_mask.GetVector(b);
    for (unsigned i1 = 0; i1 < list1_eta_.size(); i1++) {
      double index = -1;
      float min_dr = 999;
      for (unsigned i2 = 0; i2 < list2_eta_.size(); i2++) {
        if (list2_mask_[i2]) {
          float dr = deltaR(list1_eta_[i1],list1_phi_[i1],list2_eta_[i2],list2_phi_[i2]);
          if (dr < min_dr) {
            min_dr = dr;
            index = static_cast<double>(i2);
          }
        }
      }
      list2_indices.push_back(index);
    }
    return list2_indices;
  });
}

////Returns a scalar NamedFunc that is element index of vector_named_func or -999 if no index index
////only considering indices for which mask_named_func is true
//NamedFunc SortNamedFunc(NamedFunc vector_named_func, NamedFunc sort_criteria) {
//  return NamedFunc("SortNamedFunc("+vector_named_func.Name()+","+mask_named_func.Name()+")",[vector_named_func, sort_criteria](const Baby &b) -> NamedFunc::VectorType{
//    std::vector<double> vector_named_func_ = vector_named_func.GetVector(b);
//    std::vector<double> sort_criteria_ = sort_criteria.GetVector(b);
//    std::vector<std::pair<double,double>> sorting_vector;
//    std::vector<double> return_vector;
//    for (unsigned i = 0; i < vector_named_func_.size(); i++) {
//      sorting_vector.push_back(std::pair<double,double>(vector_named_func_, sort_criteria_));
//    }
//    std::sort(sorting_vector.begin(), sorting_vector.end(), [](std::pair<double,double> a, std::pair<double,double> b){a[1]>b[1]});
//    for (unsigned i = 0; i < vector_named_func_.size(); i++) {
//      return_vector.push_back(sorting_vector[i][0]);
//    }
//    return -999;
//  });
//}

int main(int argc, char *argv[]){
  //------------------------------------------------------------------------------------
  //                                    constants
  //------------------------------------------------------------------------------------
  //can only go as low as Electron/Muon_sig
  //const std::vector<double> el_pt_bins = {7.0, 10.0, 12.0, 15.0, 20.0, 25.0, 27.0, 30.0, 35.0, 40.0, 50.0, 80.0, 120.0, 200.0};
  const std::vector<double> el_pt_bins = {7.0, 12.0, 15.0, 20.0, 25.0, 27.0, 30.0, 35.0, 40.0, 50.0, 80.0, 120.0, 200.0};
  const std::vector<double> el_pt_bins_ptcut = {7.0, 15.0, 20.0, 25.0, 27.0, 30.0, 35.0, 40.0, 50.0, 80.0, 120.0, 200.0};
  const std::vector<double> el_abseta_bins = {0.0, 0.8, 1.4442, 1.566, 2.55};
  const std::vector<double> mu_pt_bins = {5.0, 7.0, 10.0, 15.0, 20.0, 23.0, 25.0, 30.0, 40.0, 50.0, 60.0, 120.0, 300.0};
  const std::vector<double> mu_pt_bins_ptcut = {5.0, 17.0, 20.0, 23.0, 25.0, 30.0, 40.0, 50.0, 60.0, 120.0, 300.0}; //for pt cut efficiencies, can only go down to 10 GeV due to cut on muon TrigObj
  const std::vector<double> mu_abseta_bins = {0.0, 0.9, 1.2, 2.1, 2.4};

  //------------------------------------------------------------------------------------
  //                                    namedfuncs needed by samples
  //------------------------------------------------------------------------------------
  /**
   * New (truth-based) Overlap removal namedfunc
   * 1 for all samples other than DY
   * 0 for DY iff there is a photon that would be duplicated in the ZGToLLG sample
   */
  const NamedFunc stitch_dy("stitch_dy",[](const Baby &b) -> NamedFunc::ScalarType{
    //extra factor for weighting is xs/nevts_effective where effective events are calculated including negative weights
    //0.007519850838359999 GluGluH->ZG->llG xs (pb)
    //967054
    if (b.FirstFileName().find("DYJetsToLL") != std::string::npos) {
      for (int imc = 0; imc < b.nGenPart(); imc++) {
        if (b.GenPart_pdgId()->at(imc)==22 && b.GenPart_pt()->at(imc)>15 && abs(b.GenPart_eta()->at(imc))<2.6) {
          bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
          if (mc_statusFlags[0]||mc_statusFlags[8]) {
            bool found_other = false;
            for (int imc2 = 0; imc2 < b.nGenPart(); imc2++) {
              bitset<15> mc2_statusFlags(b.GenPart_statusFlags()->at(imc2));
              if (b.GenPart_pt()->at(imc2)>5.0 && mc2_statusFlags[8] && b.GenPart_pdgId()->at(imc2) != 22) {
                if (deltaR(b.GenPart_eta()->at(imc), b.GenPart_phi()->at(imc), b.GenPart_eta()->at(imc2), b.GenPart_phi()->at(imc2))<0.05) {
                  found_other = true;
                  break;
                } //other particle delta R too small
              } //other particle is hardprocess, hard, and not photon
            } // /GenPart 2 loop
            if (!found_other) {
              return 0;
            }
          } // /prompt or fromhardprocess
        } // /if photon in pt-eta range of ZGToLLG
      } // /GenPart loop
    } // /if DYJetsToLL
    return 1;
  });

  const NamedFunc GenPhoton_isFsr("GenPhoton_isFsr",[](const Baby &b) -> NamedFunc::ScalarType{
    double max_ph_pt = 0;
    double is_fsr = 0;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      if (b.GenPart_pdgId()->at(imc)==22) {
        bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
        if (mc_statusFlags[0]||mc_statusFlags[8]) {
          if (b.GenPart_pt()->at(imc) > max_ph_pt) {
            max_ph_pt = b.GenPart_pt()->at(imc);
            int mom_id = abs(DistinctMomPdgId(b.GenPart_pdgId(),b.GenPart_genPartIdxMother(),imc));
            int mom_idx = imc;
            while (mom_id == 11 || mom_id == 13) {
              mom_idx = DistinctMomIdx(b.GenPart_pdgId(),b.GenPart_genPartIdxMother(),mom_idx);
              mom_id = abs(DistinctMomPdgId(b.GenPart_pdgId(),b.GenPart_genPartIdxMother(),mom_idx));
            }
            if (mom_id==23) {
              is_fsr = 1;
            }
            else {
              is_fsr = 0;
            }
          }
        }
      }
    }
    return is_fsr;
  });

  const NamedFunc Photon_sig("Photon_sig",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> Photon_sig_;
    for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
      double sig = 1;
      if (b.Photon_pt()->at(iph) < 15) sig = 0;
      if (abs(b.Photon_eta()->at(iph)) > 2.5) sig = 0;
      if (!(b.Photon_isScEtaEB()->at(iph) || b.Photon_isScEtaEE()->at(iph))) sig = 0;
      // if (b.Photon_pfRelIso03_all()->at(iph)==PhotonRelIsoCut) continue; // no isolation cut in 2016...?
      if (!(b.Photon_mvaID()->at(iph) > 0.2 && 
          b.Photon_electronVeto()->at(iph))) sig = 0;
      Photon_sig_.push_back(sig);
    }
    return Photon_sig_;
  });

  //1 iff photon is within deltar 0.4 of prompt photon with pt>10 GeV
  const NamedFunc Photon_truthMatch("Photon_truthMatch",[Photon_sig](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> ph_sig_ = Photon_sig.GetVector(b);
    double max_ph_pt = -1;
    double ph_eta = 0;
    double ph_phi = 0;
    for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
      if (ph_sig_[iph]) {
        if (b.Photon_pt()->at(iph)>max_ph_pt) {
          max_ph_pt = b.Photon_pt()->at(iph);
          ph_eta = b.Photon_eta()->at(iph);
          ph_phi = b.Photon_phi()->at(iph);
        }
      }
    }
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      if (b.GenPart_pdgId()->at(imc)==22) {
        bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
        if (mc_statusFlags[0]||mc_statusFlags[8]) {
          if (b.GenPart_pt()->at(imc)>10) {
            double this_deltar = deltaR(b.GenPart_eta()->at(imc), b.GenPart_phi()->at(imc), ph_eta, ph_phi);
            if (this_deltar < 0.4) {
              return 1;
            }
          }
        }
      }
    }
    return 0;
  });

  //1 iff photon is within deltar 0.4 of jet with pt>10 GeV
  const NamedFunc Photon_isJet("Photon_isJet",[Photon_sig](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> ph_sig_ = Photon_sig.GetVector(b);
    double max_ph_pt = -1;
    double ph_eta = 0;
    double ph_phi = 0;
    for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
      if (ph_sig_[iph]) {
        if (b.Photon_pt()->at(iph)>max_ph_pt) {
          max_ph_pt = b.Photon_pt()->at(iph);
          ph_eta = b.Photon_eta()->at(iph);
          ph_phi = b.Photon_phi()->at(iph);
        }
      }
    }
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      if ((abs(b.GenPart_pdgId()->at(imc))>=1 && abs(b.GenPart_pdgId()->at(imc))<=5) || b.GenPart_pdgId()->at(imc)==21 ) {
        bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
        if (b.GenPart_pt()->at(imc)>10) {
          double this_deltar = deltaR(b.GenPart_eta()->at(imc), b.GenPart_phi()->at(imc), ph_eta, ph_phi);
          double pt_ratio = max_ph_pt/b.GenPart_pt()->at(imc);
          if (this_deltar < 0.4 && pt_ratio>0.5 && pt_ratio<1.5) {
            return 1;
          }
        }
      }
    }
    return 0;
  });

  //1 iff photon is within deltar 0.4 of jet with pt>15 GeV
  //pico version
  const NamedFunc photon_isjet("photon_isjet",[Photon_sig](const Baby &b) -> NamedFunc::ScalarType{
    if (b.photon_pflavor()->at(0)==1 || b.photon_pflavor()->at(0)==11) return 0;
    double ph_eta = b.photon_eta()->at(0);
    double ph_phi = b.photon_phi()->at(0);
    for (unsigned imc = 0; imc < b.mc_pt()->size(); imc++) {
      if ((abs(b.mc_id()->at(imc))>=1 && abs(b.mc_id()->at(imc))<=5) || b.mc_id()->at(imc)==21 ) {
        if (b.mc_pt()->at(imc)>15) {
          double this_deltar = deltaR(b.mc_eta()->at(imc), b.mc_phi()->at(imc), ph_eta, ph_phi);
          if (this_deltar < 0.4) {
            return 1;
          }
        }
      }
    }
    return 0;
  });

  //true if the higgs candidate has a photon and an electron that are the same
  const NamedFunc overlapping_phel("overlapping_phel",[Photon_sig](const Baby &b) -> NamedFunc::ScalarType{
    if (b.ll_lepid()->at(0)!=11) return 0;
    if (b.photon_elidx()->at(0)==b.ll_i1()->at(0) || 
        b.photon_elidx()->at(0)==b.ll_i2()->at(0)) return 1;
    return 0;
  });

  //------------------------------------------------------------------------------------
  //                                    initialization
  //------------------------------------------------------------------------------------
  gErrorIgnoreLevel = 6000;
  time_t begtime, endtime;
  time(&begtime);
  Palette colors("txt/colors.txt","default");
  script_utilities::ArgStruct options = script_utilities::get_options(
      argc, argv, "");
  std::vector<PlotOpt> plt_lin = script_utilities::plot_lin(false);
  std::vector<PlotOpt> plt_lin_over = script_utilities::plot_lin(false);
  plt_lin_over[0] = plt_lin_over[0].Overflow(OverflowType::both);
  std::vector<PlotOpt> plt_log = script_utilities::plot_log(false);
  std::vector<PlotOpt> plt_shapes = script_utilities::plot_shapes();
  std::vector<PlotOpt> plt_log_shapes = script_utilities::plot_log_shapes();
  PlotOpt style2D("txt/plot_styles.txt", "Scatter");
  std::vector<PlotOpt> twodim_plotopts = {style2D().Title(TitleType::info)
      .YAxis(YAxisType::linear).Overflow(OverflowType::overflow).LogMinimum(0.001)};
  std::vector<PlotOpt> twodim_log_plotopts = {style2D().Title(TitleType::info)
      .YAxis(YAxisType::log).Overflow(OverflowType::overflow).LogMinimum(0.001)};
  std::vector<PlotOpt> plt_lin_logx = {PlotOpt("txt/plot_styles.txt","CMSPaper")
      .Title(TitleType::info).Bottom(PlotOptTypes::BottomType::off)
      .XAxis(PlotOptTypes::YAxisType::log)
      .Overflow(PlotOptTypes::OverflowType::none).LegendColumns(3)};

  //set<int> years;
  //HigUtilities::parseYears(options.year_string, years);
  //int lumi_precision = 0;
  //string total_luminosity_string = HigUtilities::getLuminosityString(options.year_string, lumi_precision);
  std::string total_luminosity_string = "41.5";

  std::string base_folder = "/net/cms17/cms17r0/pico/NanoAODv7/nano/2017/signal/";
  std::string base_folder_ul = "/net/cms17/cms17r0/pico/NanoAODv2/nano/2017/";
  //std::string base_folder_ul = "/net/cms17/cms17r0/pico/NanoAODv9/nano/2017/";
  std::string base_folder_v7data = "/net/cms17/cms17r0/pico/NanoAODv7/nano/2017/data/";
  std::string base_folder_v7 = "/net/cms17/cms17r0/pico/NanoAODv7/nano/";

  std::vector<std::shared_ptr<Process>> procs;
  procs.push_back(Process::MakeShared<Baby_nano>("HToZG", 
      Process::Type::signal, kBlack,
      {base_folder+"*GluGluHToZG_ZToLL*"},"1"));

  //2017 triggers
  NamedFunc hlt_el_trigger = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL||HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ";
  NamedFunc hlt_mu_trigger = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8||HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8";
  NamedFunc hlt_single_el_trigger = "HLT_Ele35_WPTight_Gsf||HLT_Ele32_WPTight_Gsf_L1DoubleEG";
  NamedFunc hlt_single_mu_trigger = "HLT_IsoMu27||HLT_IsoMu24";
  //
  NamedFunc hlt_el_trigger_plus = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL||HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ||HLT_Ele35_WPTight_Gsf||HLT_Ele32_WPTight_Gsf_L1DoubleEG";
  NamedFunc hlt_mu_trigger_plus = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8||HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8||HLT_IsoMu27||HLT_IsoMu24";
  NamedFunc hlt_mu_trigger_plusplus = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8||HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8||HLT_IsoMu27||HLT_IsoMu24||HLT_Mu17_Photon30_IsoCaloId";
  NamedFunc l1_el_trigger = "L1_SingleEG24||L1_SingleEG34er2p1||L1_SingleIsoEG24er2p1||L1_SingleIsoEG24||L1_DoubleEG_18_17||L1_DoubleEG_22_10||L1_DoubleEG_LooseIso23_10||L1_DoubleEG_LooseIso24_10";
  NamedFunc l1_mu_trigger = "L1_DoubleMu_12_5||L1_DoubleMu_15_7_SQ||L1_DoubleMu_15_7_SQ_Mass_Min4";
  NamedFunc l1_mu_trigger_singlemu = "L1_SingleMu22||L1_SingleMu25";
  NamedFunc hlt_double_el_noiso_trigger = "HLT_DoubleEle25_CaloIdL_MW";
  NamedFunc hlt_double_mu_noiso_trigger = "HLT_Mu37_TkMu27";
  NamedFunc hlt_single_el_noiso_trigger = "HLT_Ele115_CaloIdVT_GsfTrkIdT";
  NamedFunc hlt_single_mu_noiso_trigger = "HLT_Mu50||HLT_Mu55";
  NamedFunc hlt_single_el_prescale_trigger = "HLT_Ele25_WPTight_Gsf||HLT_Ele27_WPTight_Gsf||HLT_Ele28_WPTight_Gsf||HLT_Ele32_WPTight_Gsf||HLT_Ele20_WPLoose_Gsf||HLT_Ele45_WPLoose_Gsf||"
                                             "HLT_Ele25_eta2p1_WPTight_Gsf||HLT_Ele27_eta2p1_WPTight_Gsf||HLT_Ele20_eta2p1_WPLoose_Gsf||HLT_Ele25_eta2p1_WPLoose_Gsf||HLT_Ele27_eta2p1_WPLoose_Gsf";
  NamedFunc hlt_single_mu_prescale_trigger = "HLT_IsoMu20||HLT_IsoMu22||HLT_IsoMu22_eta2p1||HLT_IsoMu24_eta2p1||HLT_IsoTkMu20||HLT_IsoTkMu22||HLT_IsoTkMu24||HLT_Mu45_eta2p1||HLT_TkMu50";
  NamedFunc hlt_mu_ph_trigger = "HLT_Mu17_Photon30_IsoCaloId";
  NamedFunc hlt_doubleph_trigger = "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55||HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90";
  NamedFunc l1_lep_trigger = l1_el_trigger||l1_mu_trigger;
  NamedFunc hlt_lep_trigger = hlt_el_trigger||hlt_mu_trigger;
  NamedFunc hlt_lep_trigger_plus = hlt_el_trigger_plus||hlt_mu_trigger_plus;
  NamedFunc l1_mu_eg_trigger = "L1_Mu5_EG23||L1_Mu5_LooseIsoEG20";

  //2018 triggers
  NamedFunc hlt_singleel_trigger_2018 = "HLT_Ele32_WPTight_Gsf";
  NamedFunc hlt_singlemu_trigger_2018 = "HLT_IsoMu27";
  NamedFunc hlt_doubleel_trigger_2018 = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL||HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ";
  NamedFunc hlt_doublemu_trigger_2018 = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8";
  NamedFunc hlt_doubleph_trigger_2018 = "HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto||HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90";
  NamedFunc l1_eg_trigger_2018 = "L1_SingleEG26er2p5||L1_SingleIsoEG24er2p1||L1_SingleIsoEG26er2p5||L1_DoubleEG_22_10_er2p5";
  NamedFunc l1_doublemu_trigger_2018 = "L1_DoubleMu_12_5||L1_DoubleMu_15_7";
  NamedFunc l1_mu_eg_trigger_2018 = "L1_Mu5_EG23er2p5||L1_Mu5_LooseIsoEG20er2p5";
  //l1 single mu, hlt mu+ph same

  std::vector<std::shared_ptr<Process>> full_procs;
  full_procs.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+jets", 
      Process::Type::background, TColor::GetColor("#ffb400"),
      {base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__0*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__1*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__2*"},stitch_dy));
  full_procs.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+#gamma", 
      Process::Type::background, TColor::GetColor("#16bac5"),
      {base_folder_ul+"mc/ZGToLLG_01J_5f_TuneCP5*"},"1"));
  full_procs.push_back(Process::MakeShared<Baby_nano>("gg#rightarrow H#rightarrow ZG (x100)", 
      Process::Type::signal, TColor::GetColor("#ff0000"),
      {base_folder_ul+"signal/GluGluHToZG*"},"1"));

  std::vector<std::shared_ptr<Process>> full_procs_noscale;
  full_procs_noscale.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+jets", 
      Process::Type::background, TColor::GetColor("#ffb400"),
      {base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__0*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__1*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__2*"},stitch_dy));
  full_procs_noscale.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+#gamma", 
      Process::Type::background, TColor::GetColor("#16bac5"),
      {base_folder_ul+"mc/ZGToLLG_01J_5f_TuneCP5*"},"1"));
  full_procs_noscale.push_back(Process::MakeShared<Baby_nano>("gg#rightarrow H#rightarrow ZG", 
      Process::Type::signal, TColor::GetColor("#ff0000"),
      {base_folder_ul+"signal/GluGluHToZG*"},"1"));

  std::vector<std::shared_ptr<Process>> procs_dy;
  procs_dy.push_back(Process::MakeShared<Baby_nano>("Z/#gamma* (DYJets) (FSR)", 
      Process::Type::background, TColor::GetColor("#e6cdfe"),
      {base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__0*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__1*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__2*"},stitch_dy&&GenPhoton_isFsr));
  procs_dy.push_back(Process::MakeShared<Baby_nano>("Z/#gamma* (ZGToLLG) (FSR)", 
      Process::Type::background, TColor::GetColor("#cc99fe"),
      {base_folder_ul+"mc/ZGToLLG_01J_5f_TuneCP5*"},stitch_dy&&GenPhoton_isFsr));
  procs_dy.push_back(Process::MakeShared<Baby_nano>("Z/#gamma* (DYJets) (ISR)", 
      Process::Type::background, TColor::GetColor("#30dce8"),
      {base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__0*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__1*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__2*"},stitch_dy&&!GenPhoton_isFsr));
  procs_dy.push_back(Process::MakeShared<Baby_nano>("Z/#gamma* (ZGToLLG) (ISR)", 
      Process::Type::background, TColor::GetColor("#16bac5"),
      {base_folder_ul+"mc/ZGToLLG_01J_5f_TuneCP5*"},stitch_dy&&!GenPhoton_isFsr));

  std::vector<std::shared_ptr<Process>> procs_dy_incl;
  procs_dy_incl.push_back(Process::MakeShared<Baby_nano>("Z/#gamma* (DYJets)", 
      Process::Type::background, TColor::GetColor("#16bac5"),
      {base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__0*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__1*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__2*"},stitch_dy));
  procs_dy_incl.push_back(Process::MakeShared<Baby_nano>("Z/#gamma* (ZGToLLG)", 
      Process::Type::background, TColor::GetColor("#16bac5"),
      {base_folder_ul+"mc/ZGToLLG_01J_5f_TuneCP5*"},stitch_dy));

  std::vector<std::shared_ptr<Process>> procs_dy_nostitch;
  procs_dy_nostitch.push_back(Process::MakeShared<Baby_nano>("Z/#gamma* (DYJets,no stitch) (FSR)", 
      Process::Type::background, TColor::GetColor("#cc99fe"),
      {base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__0*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__1*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__2*"},GenPhoton_isFsr));
  procs_dy_nostitch.push_back(Process::MakeShared<Baby_nano>("Z/#gamma* (DYJets,no stitch) (ISR)", 
      Process::Type::background, TColor::GetColor("#16bac5"),
      {base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__0*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__1*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__2*"},!GenPhoton_isFsr));

  std::vector<std::shared_ptr<Process>> procs_photon_all;
  procs_photon_all.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+Fake Photon (Jet)", 
      Process::Type::background, TColor::GetColor("#ffb400"),
      {base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__0*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__1*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__2*",
      base_folder_ul+"mc/ZGToLLG_01J_5f_TuneCP5*"},stitch_dy&&!Photon_truthMatch&&Photon_isJet));
  procs_photon_all.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+Fake Photon (Other)", 
      Process::Type::background, TColor::GetColor("#ffe956"),
      {base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__0*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__1*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__2*",
      base_folder_ul+"mc/ZGToLLG_01J_5f_TuneCP5*"},stitch_dy&&!Photon_truthMatch&&!Photon_isJet));
  procs_photon_all.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+#gamma (DYJets) (FSR)", 
      Process::Type::background, TColor::GetColor("#e6cdfe"),
      {base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__0*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__1*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__2*"},stitch_dy&&Photon_truthMatch&&GenPhoton_isFsr));
  procs_photon_all.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+#gamma (DYJets) (ISR)", 
      Process::Type::background, TColor::GetColor("#30dce8"),
      {base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__0*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__1*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__2*"},stitch_dy&&Photon_truthMatch&&!GenPhoton_isFsr));
  procs_photon_all.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+#gamma (ZGToLLG) (FSR)", 
      Process::Type::background, TColor::GetColor("#cc99fe"),
      {base_folder_ul+"mc/ZGToLLG_01J_5f_TuneCP5*"},stitch_dy&&Photon_truthMatch&&GenPhoton_isFsr));
  procs_photon_all.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+#gamma (ZGToLLG) (ISR)", 
      Process::Type::background, TColor::GetColor("#16bac5"),
      {base_folder_ul+"mc/ZGToLLG_01J_5f_TuneCP5*"},stitch_dy&&Photon_truthMatch&&!GenPhoton_isFsr));

  std::vector<std::shared_ptr<Process>> procs_photon_jet;
  procs_photon_jet.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+Fake Photon (Jet)", 
      Process::Type::background, TColor::GetColor("#ffb400"),
      {base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__0*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__1*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__2*",
      base_folder_ul+"mc/ZGToLLG_01J_5f_TuneCP5*"},stitch_dy&&!Photon_truthMatch&&Photon_isJet));

  std::vector<std::shared_ptr<Process>> procs_photon_other;
  procs_photon_other.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+Fake Photon (Other)", 
      Process::Type::background, TColor::GetColor("#ffe956"),
      {base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__0*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__1*",
      base_folder_ul+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__2*",
      base_folder_ul+"mc/ZGToLLG_01J_5f_TuneCP5*"},stitch_dy&&!Photon_truthMatch&&!Photon_isJet));

  std::string prod_folder("/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v0p5/");
  std::set<std::string> years = {"2018"};
  NamedFunc trig_and_stitch = (ZgFunctions::HLT_pass_dilepton||ZgFunctions::HLT_pass_singlelepton)
                              &&"use_event";
  Palette mc_colors("txt/colors_zgamma.txt","default");
  auto proc_zg =   Process::MakeShared<Baby_pico>("ZGToLLG", Process::Type::background, 
                   mc_colors("zgtollg"), attach_folder(prod_folder,years,"mc/merged_zgmc_llg",
                   {"*ZGToLLG_01J_5f_Tune*"}), 
                   trig_and_stitch);
  auto proc_dyal = Process::MakeShared<Baby_pico>("DYJets (real photon)", Process::Type::background, 
                   mc_colors("dyjets"), attach_folder(prod_folder,years,"mc/merged_zgmc_llg",
                   {"*DYJets*amcatnlo*"}), 
                   trig_and_stitch);
  //auto proc_dyph = Process::MakeShared<Baby_pico>("DYJets (real photon)", Process::Type::background, 
  //                 kYellow-9, attach_folder(prod_folder,years,"mc/merged_zgmc_llg",
  //                 {"*DYJets*amcatnlo*"}), 
  //                 trig_and_stitch&&"photon_pflavor[0]==1");
  //auto proc_dyel = Process::MakeShared<Baby_pico>("DYJets (photon is electron)", 
  //                 Process::Type::background, 
  //                 kViolet-1, attach_folder(prod_folder,years,"mc/merged_zgmc_llg",
  //                 {"*DYJets*amcatnlo*"}), 
  //                 trig_and_stitch&&"photon_pflavor[0]==11");
  //auto proc_dyjt = Process::MakeShared<Baby_pico>("DYJets (photon is jet)", 
  //                 Process::Type::background, 
  //                 mc_colors("dyjets"), attach_folder(prod_folder,years,"mc/merged_zgmc_llg",
  //                 {"*DYJets*amcatnlo*"}), 
  //                 trig_and_stitch&&photon_isjet);
  //auto proc_dyot = Process::MakeShared<Baby_pico>("DYJets (other fake photon)", 
  //                 Process::Type::background, 
  //                 kGreen+2, attach_folder(prod_folder,years,"mc/merged_zgmc_llg",
  //                 {"*DYJets*amcatnlo*"}), 
  //                 trig_and_stitch&&"photon_pflavor[0]==0"&&!photon_isjet);
  auto proc_hzg  = Process::MakeShared<Baby_pico>("H#rightarrow Z#gamma", 
                   Process::Type::signal, 
                   kRed, attach_folder(prod_folder,years,"mc/merged_zgmc_llg",
                   {"*GluGluHToZG*M-125*"}), 
                   trig_and_stitch);
  //std::vector<std::shared_ptr<Process>> procs_dysplit = {proc_dyel, proc_dyph, proc_dyjt, 
  //                                                       proc_dyot, proc_zg, proc_hzg};
  std::vector<std::shared_ptr<Process>> procs_dysplit = {proc_dyal, proc_zg, proc_hzg};

  //------------------------------------------------------------------------------------
  //                      named funcs to be moved to a common location
  //------------------------------------------------------------------------------------

  const NamedFunc Electron_sig("Electron_sig",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> Electron_sig_;
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      bool is_sig = true;
      if (b.Electron_pt()->at(iel) < 7) is_sig = false; //was 15
      float etasc = b.Electron_deltaEtaSC()->at(iel) + b.Electron_eta()->at(iel);
      if (abs(etasc) > 2.5) is_sig = false;
      if (abs(b.Electron_dz()->at(iel))>1.0) is_sig = false;
      if (abs(b.Electron_dxy()->at(iel))>0.5) is_sig = false; 
      double wp[2][3] = {{-0.145237, -0.0315746, -0.032173},
                         { 0.604775,  0.628743,   0.896462}};
      int ipt(1), ieta(2);
      if(b.Electron_pt()->at(iel)>10) ipt = 0;
      if(fabs(etasc) < 0.8) ieta = 0;
      else if(fabs(etasc) < 1.479) ieta = 1;
      double mva = b.Electron_mvaFall17V2Iso()->at(iel);
      if (mva <= wp[ipt][ieta])
        is_sig = false;
      if (is_sig)
        Electron_sig_.push_back(1.0);
      else
        Electron_sig_.push_back(0.0);
    }
    return Electron_sig_;
  });

  const NamedFunc Muon_sig("Muon_sig",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> Muon_sig_;
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      bool is_sig = true;
      if (b.Muon_pt()->at(imu) < 5) is_sig = false; //was 10
      if (abs(b.Muon_eta()->at(imu)) > 2.4) is_sig = false;
      if (abs(b.Muon_dz()->at(imu))>1.0)  is_sig = false;
      if (abs(b.Muon_dxy()->at(imu))>0.5) is_sig = false; 
      if (!((b.Muon_looseId()->at(imu) || (b.Muon_pt()->at(imu) > 200 && b.Muon_highPtId()->at(imu))) && 
          b.Muon_pfRelIso03_all()->at(imu) < 0.35 &&
          b.Muon_sip3d()->at(imu) < 4)) is_sig = false;
      if (is_sig)
        Muon_sig_.push_back(1.0);
      else
        Muon_sig_.push_back(0.0);
    }
    return Muon_sig_;
  });

  const NamedFunc NMuon_sig = SumNamedFunc(Muon_sig); //pt>5
  const NamedFunc NElectron_sig = SumNamedFunc(Electron_sig); //pt>7
  const NamedFunc NLepton_sig = NMuon_sig+NElectron_sig;
  const NamedFunc NPhoton_sig = SumNamedFunc(Photon_sig); //pt>15

  const NamedFunc NMuon_trig = SumNamedFunc("TrigObj_id==13");
  const NamedFunc NElectron_trig = SumNamedFunc("TrigObj_id==11");
  const NamedFunc NPhoton_trig = SumNamedFunc("TrigObj_id==22");

  //const NamedFunc Lead_Electron_pt = AtNamedFunc("Electron_pt", Electron_sig, 0);
  //const NamedFunc Sublead_Electron_pt = AtNamedFunc("Electron_pt", Electron_sig, 1);
  //const NamedFunc Lead_Muon_pt = AtNamedFunc("Muon_pt", Electron_sig, 0);
  //const NamedFunc Sublead_Muon_pt = AtNamedFunc("Muon_pt", Electron_sig, 1);
  //const NamedFunc Lead_Photon_pt = AtNamedFunc("Photon_pt", Photon_sig, 0);
  
  //const NamedFunc Lead_Electron_pt = FilterNamedFunc("Electron_pt", Electron_sig)[static_cast<int>(0)];
  //const NamedFunc Sublead_Electron_pt = FilterNamedFunc("Electron_pt", Electron_sig)[1];
  //const NamedFunc Lead_Muon_pt = FilterNamedFunc("Muon_pt", Muon_sig)[static_cast<int>(0)];
  //const NamedFunc Sublead_Muon_pt = FilterNamedFunc("Muon_pt", Muon_sig)[1];
  //const NamedFunc Lead_Photon_pt = FilterNamedFunc("Photon_pt", Photon_sig)[static_cast<int>(0)];

  //TODO: fix these
  const NamedFunc Lead_Electron_pt("Lead_Electron_pt",[Electron_sig](const Baby &b) -> NamedFunc::ScalarType{
    double el_pt = 0;
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_sig_[iel]) {
        if (b.Electron_pt()->at(iel) > el_pt)
          el_pt = b.Electron_pt()->at(iel);
      }
    }
    return el_pt;
  });

  const NamedFunc Sublead_Electron_pt("Sublead_Electron_pt",[Electron_sig](const Baby &b) -> NamedFunc::ScalarType{
    double sublead_el_pt = 0;
    double lead_el_pt = 0;
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_sig_[iel]) {
        if (b.Electron_pt()->at(iel) > lead_el_pt) {
          sublead_el_pt = lead_el_pt;
          lead_el_pt = b.Electron_pt()->at(iel);
        }
        else if (b.Electron_pt()->at(iel) > sublead_el_pt)
          sublead_el_pt = b.Electron_pt()->at(iel);
      }
    }
    return sublead_el_pt;
  });

  const NamedFunc Lead_Muon_pt("Lead_Muon_pt",[Muon_sig](const Baby &b) -> NamedFunc::ScalarType{
    double mu_pt = 0;
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      if (Muon_sig_[imu]) {
        if (b.Muon_pt()->at(imu) > mu_pt)
          mu_pt = b.Muon_pt()->at(imu);
      }
    }
    return mu_pt;
  });

  const NamedFunc Sublead_Muon_pt("Sublead_Muon_pt",[Muon_sig](const Baby &b) -> NamedFunc::ScalarType{
    double lead_mu_pt = 0;
    double sublead_mu_pt = 0;
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      if (Muon_sig_[imu]) {
        if (b.Muon_pt()->at(imu) > lead_mu_pt) {
          sublead_mu_pt = lead_mu_pt;
          lead_mu_pt = b.Muon_pt()->at(imu);
        }
        else if (b.Muon_pt()->at(imu) > sublead_mu_pt)
          sublead_mu_pt = b.Muon_pt()->at(imu);
      }
    }
    return sublead_mu_pt;
  });

  const NamedFunc Lead_Electron_abseta("Lead_Electron_abseta",[Electron_sig](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_sig_[iel]) {
        return abs(b.Electron_eta()->at(iel));
      }
    }
    return -999;
  });

  const NamedFunc Sublead_Electron_abseta("Sublead_Electron_abseta",[Electron_sig](const Baby &b) -> NamedFunc::ScalarType{
    int sig_idx = 0;
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_sig_[iel]) {
        if (sig_idx==1) {
          return abs(b.Electron_eta()->at(iel));
        }
        sig_idx++;
      }
    }
    return -999;
  });

  const NamedFunc Lead_Muon_abseta("Lead_Muon_abseta",[Muon_sig](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      if (Muon_sig_[imu]) {
        return abs(b.Muon_eta()->at(imu));
      }
    }
    return -999;
  });

  const NamedFunc Sublead_Muon_abseta("Sublead_Muon_abseta",[Muon_sig](const Baby &b) -> NamedFunc::ScalarType{
    int sig_idx = 0;
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      if (Muon_sig_[imu]) {
        if (sig_idx==1)
          return abs(b.Muon_eta()->at(imu));
        sig_idx++;
      }
    }
    return -999;
  });

  const NamedFunc Lead_Photon_pt("Lead_Photon_pt",[Photon_sig](const Baby &b) -> NamedFunc::ScalarType{
    double ph_pt = 0;
    std::vector<double> ph_sig_ = Photon_sig.GetVector(b);
    for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
      if (ph_sig_[iph]) {
        if (b.Photon_pt()->at(iph) > ph_pt)
          ph_pt = b.Photon_pt()->at(iph);
      }
    }
    return ph_pt;
  });

  //------------------------------------------------------------------------------------
  //                                     named funcs
  //------------------------------------------------------------------------------------
  
  const NamedFunc weight("weight",[](const Baby &b) -> NamedFunc::ScalarType{
    double h_zg_br = 0.001533;
    double z_ll_br = 0.100974;
    //extra factor for weighting is xs/nevts_effective where effective events are calculated including negative weights
    //0.007519850838359999 GluGluH->ZG->llG xs (pb)
    //967054
    //double xs = 7.7760403e-09; //GluGluHToZG
    double xs = 0.00751985;
    double nents = 400000.0;
    if (b.FirstFileName().find("VBFHToZG") != std::string::npos) {
      xs = h_zg_br*z_ll_br*3.782;
      nents = 200000;
    }
    if (b.FirstFileName().find("WminusH_HToZG") != std::string::npos) {
      xs = h_zg_br*0.527;
      nents = 299276;
    }
    if (b.FirstFileName().find("WplusH_HToZG") != std::string::npos) {
      xs = h_zg_br*0.831;
      nents = 299978;
    }
    if (b.FirstFileName().find("ZH_HToZG") != std::string::npos) {
      xs = h_zg_br*0.8839;
      nents = 297389;
    }
    if (b.FirstFileName().find("ttHToZG") != std::string::npos) {
      xs = h_zg_br*0.5071;
      nents = 200000;
    }
    if (b.FirstFileName().find("ZGToLLG") != std::string::npos) {
      xs = 117.864;
      //nents = 2.8095231e+09;
      //nents = 1240092;
      //nents = 9315642;
      nents = 29885702;
    }
    if (b.FirstFileName().find("DYJetsToLL") != std::string::npos) {
      xs = 6077.22;
      //nents = 3.2011619e+11;
      //nents = 3535863;
      //nents = 7680238;
      nents = 19014526;
      if (b.SampleType()==2016) nents = 3865459; 
    }
    //float w_lumi = 1.0;
    float w_year = 1.0;
    //if (b.Generator_weight()<0) w_lumi = -1.0;
    if (b.SampleType()==2016) w_year = 36.32264; 
    else if (b.SampleType()==2017) w_year = 41.52756;
    else if (b.SampleType()==2018) w_year = 59.67377;
    return w_year*xs*1000.0/nents;
  });

  const NamedFunc PicoType("PicoType",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.FirstFileName().find("ZGToLLG") != std::string::npos) {
      return 17100;
    }
    if (b.FirstFileName().find("DYJetsToLL") != std::string::npos) {
      return 6000;
    }
    return -1;
  });

  const NamedFunc Signal_x20("Signal_x20",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.FirstFileName().find("HToZG") != std::string::npos) {
      return 20.0;
    }
    return 1.0;
  });

  const NamedFunc z_decay_pdgid("z_decay_pdgid",[](const Baby &b) -> NamedFunc::ScalarType{
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid != 21 && abs_pdgid != 22 && abs_pdgid != 23) { //skip radiated particles and Z
        int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
        if (mom_idx != -1) {
          if ((b.GenPart_pdgId()->at(mom_idx))==23) {
            int z_mom_idx = b.GenPart_genPartIdxMother()->at(mom_idx);
            while (true) {
              if (z_mom_idx != -1) {
                int z_mom_id = b.GenPart_pdgId()->at(z_mom_idx);
                if (z_mom_id == 23) {
                  z_mom_idx = b.GenPart_genPartIdxMother()->at(z_mom_idx);
                }
                else {
                  if (z_mom_id == 25) { //Z is from H
                    return abs_pdgid;
                  }
                  break;
                }
              }
              else {
                break;
              }
            } //looking for Z parent
          }
        }
      }
    }
    return 0;
  });

  const NamedFunc Lead_Photon_mvaID("Lead_Photon_mvaID",[Photon_sig](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> Photon_sig_ = Photon_sig.GetVector(b);
    for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
      if (Photon_sig_[iph]) {
        return b.Photon_mvaID()->at(iph);
      }
    }
    return 0;
  });

  const NamedFunc LeadMuon_absEta("LeadMuon_absEta",[Muon_sig](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      if (Muon_sig_[imu]) {
        return fabs(b.Muon_eta()->at(imu));
      }
    }
    return 0;
  });

  const NamedFunc TrigObj_filterBit2("TrigObj_filterBit2",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> filterbit2;
    for (unsigned itrig = 0; itrig < b.TrigObj_pt()->size(); itrig++) {
      filterbit2.push_back(static_cast<double>((b.TrigObj_filterBits()->at(itrig) & 0x02) >> 1));
    }
    return filterbit2;
  });

  const NamedFunc Jet_sig("Jet_sig",[Electron_sig, Muon_sig, Photon_sig](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> Jet_sig_;
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    std::vector<double> ph_sig_ = Photon_sig.GetVector(b);
    for (unsigned ijet = 0; ijet < b.Jet_pt()->size(); ijet++) {
      double sig = 1;
      if (b.Jet_pt()->at(ijet) < 30) sig = 0;
      if (abs(b.Jet_eta()->at(ijet)) > 2.4) sig = 0;
      for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
        if (Electron_sig_[iel] > 0.5) {
          if (deltaR(b.Electron_eta()->at(iel), b.Electron_phi()->at(iel), b.Jet_eta()->at(ijet), b.Jet_phi()->at(ijet)) < 0.4)
            sig = 0;
        }
      }
      for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
        if (Muon_sig_[imu] > 0.5) {
          if (deltaR(b.Muon_eta()->at(imu), b.Muon_phi()->at(imu), b.Jet_eta()->at(ijet), b.Jet_phi()->at(ijet)) < 0.4)
            sig = 0;
        }
      }
      for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
        if (ph_sig_[iph] > 0.5) {
          if (deltaR(b.Photon_eta()->at(iph), b.Photon_phi()->at(iph), b.Jet_eta()->at(ijet), b.Jet_phi()->at(ijet)) < 0.4)
            sig = 0;
        }
      }
      Jet_sig_.push_back(sig);
    }
    return Jet_sig_;
  });

  const NamedFunc nJet_sig("nJet_sig",[Jet_sig](const Baby &b) -> NamedFunc::ScalarType{
    int njet_sig = 0;
    std::vector<double> jet_sig_ = Jet_sig.GetVector(b);
    for (double ijet_sig_ : jet_sig_) {
      if (ijet_sig_ > 0.5)
        njet_sig++;
    }
    return njet_sig;
  });

  const NamedFunc Electron_IsoPhotonPt("Electron_IsoPhotonPt",[Photon_sig](const Baby &b) -> NamedFunc::VectorType{
    //-1 - no HLT object, 0 - fails all, 1 - pass CaloIdL_TrackIdL_IsoVL, 2 - pass WPTight
    std::vector<double> el_isophpt;
    std::vector<double> ph_sig_ = Photon_sig.GetVector(b);
    double isophpt = 0;
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      isophpt = 0;
      for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
        if (ph_sig_[iph]>0.5) {
          if (deltaR(b.Electron_eta()->at(iel), b.Electron_phi()->at(iel), b.Photon_eta()->at(iph), b.Photon_phi()->at(iph)) < 0.3) {
            isophpt += b.Photon_pt()->at(iph);
          }
        }
      }
      el_isophpt.push_back(isophpt);
    }
    return el_isophpt;
  });

  const NamedFunc Electron_IsoElectronPt("Electron_IsoElectronPt",[Electron_sig](const Baby &b) -> NamedFunc::VectorType{
    //-1 - no HLT object, 0 - fails all, 1 - pass CaloIdL_TrackIdL_IsoVL, 2 - pass WPTight
    std::vector<double> el_isoelpt;
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    double isoelpt = 0;
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      isoelpt = 0;
      for (unsigned iel2 = 0; iel2 < b.Electron_pt()->size(); iel2++) {
        if (iel == iel2) continue;
        if (Electron_sig_[iel2]>0.5) {
          if (deltaR(b.Electron_eta()->at(iel), b.Electron_phi()->at(iel), b.Electron_eta()->at(iel2), b.Electron_phi()->at(iel2)) < 0.3) {
            isoelpt += b.Electron_pt()->at(iel2);
          }
        }
      }
      el_isoelpt.push_back(isoelpt);
    }
    return el_isoelpt;
  });

  const NamedFunc ZCand_mass("ZCand_mass",[Electron_sig, Muon_sig](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    double m_ll = -999;
    for (unsigned iel = 1; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_sig_[iel]) {
        for (unsigned iel2 = 0; iel2 < iel; iel2++) {
          if (Electron_sig_[iel2]) {
            ROOT::Math::PtEtaPhiMVector p1(b.Electron_pt()->at(iel),b.Electron_eta()->at(iel),b.Electron_phi()->at(iel),0.000511);
            ROOT::Math::PtEtaPhiMVector p2(b.Electron_pt()->at(iel2),b.Electron_eta()->at(iel2),b.Electron_phi()->at(iel2),0.000511);
            double this_mll = (p1+p2).M();
            if (fabs(this_mll-91.2)<fabs(m_ll-91.2))
              m_ll = this_mll;
          }
        }
      }
    }
    for (unsigned imu = 1; imu < b.Muon_pt()->size(); imu++) {
      if (Muon_sig_[imu]) {
        for (unsigned imu2 = 0; imu2 < imu; imu2++) {
          if (Muon_sig_[imu2]) {
            ROOT::Math::PtEtaPhiMVector p1(b.Muon_pt()->at(imu),b.Muon_eta()->at(imu),b.Muon_phi()->at(imu),0.106);
            ROOT::Math::PtEtaPhiMVector p2(b.Muon_pt()->at(imu2),b.Muon_eta()->at(imu2),b.Muon_phi()->at(imu2),0.106);
            double this_mll = (p1+p2).M();
            if (fabs(this_mll-91.2)<fabs(m_ll-91.2))
              m_ll = this_mll;
          }
        }
      }
    }
    return m_ll;
  });

  const NamedFunc HiggsCand_mass("HiggsCand_mass",[Electron_sig, Muon_sig, Photon_sig](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    std::vector<double> ph_sig_ = Photon_sig.GetVector(b);
    ROOT::Math::PtEtaPhiMVector ll_p;
    double m_ll = -999;
    double m_llg = -999;
    for (unsigned iel = 1; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_sig_[iel]) {
        for (unsigned iel2 = 0; iel2 < iel; iel2++) {
          if (Electron_sig_[iel2]) {
            ROOT::Math::PtEtaPhiMVector p1(b.Electron_pt()->at(iel),b.Electron_eta()->at(iel),b.Electron_phi()->at(iel),0.000511);
            ROOT::Math::PtEtaPhiMVector p2(b.Electron_pt()->at(iel2),b.Electron_eta()->at(iel2),b.Electron_phi()->at(iel2),0.000511);
            double this_mll = (p1+p2).M();
            if (fabs(this_mll-91.2)<fabs(m_ll-91.2)) {
              ll_p = p1+p2;
              m_ll = this_mll;
            }
          }
        }
      }
    }
    for (unsigned imu = 1; imu < b.Muon_pt()->size(); imu++) {
      if (Muon_sig_[imu]) {
        for (unsigned imu2 = 0; imu2 < imu; imu2++) {
          if (Muon_sig_[imu2]) {
            ROOT::Math::PtEtaPhiMVector p1(b.Muon_pt()->at(imu),b.Muon_eta()->at(imu),b.Muon_phi()->at(imu),0.106);
            ROOT::Math::PtEtaPhiMVector p2(b.Muon_pt()->at(imu2),b.Muon_eta()->at(imu2),b.Muon_phi()->at(imu2),0.106);
            double this_mll = (p1+p2).M();
            if (fabs(this_mll-91.2)<fabs(m_ll-91.2)) {
              ll_p = p1+p2; 
              m_ll = this_mll;
            }
          }
        }
      }
    }
    if (m_ll > 0) {
      for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
        if (ph_sig_[iph]) {
          ROOT::Math::PtEtaPhiMVector pph(b.Photon_pt()->at(iph),b.Photon_eta()->at(iph),b.Photon_phi()->at(iph),0.0);
          double this_mass = (pph+ll_p).M();
          if (fabs(this_mass-125.3)<fabs(m_llg-125.3))
            m_llg = this_mass;
        }
      }
    }
    return m_llg;
  });

  const NamedFunc HiggsCand_ptOverMass("HiggsCand_ptOverMass",[Electron_sig, Muon_sig, Photon_sig](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    std::vector<double> ph_sig_ = Photon_sig.GetVector(b);
    ROOT::Math::PtEtaPhiMVector ll_p;
    double m_ll = -999;
    double pt_llg = 0;
    double m_llg = -999;
    for (unsigned iel = 1; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_sig_[iel]) {
        for (unsigned iel2 = 0; iel2 < iel; iel2++) {
          if (Electron_sig_[iel2]) {
            ROOT::Math::PtEtaPhiMVector p1(b.Electron_pt()->at(iel),b.Electron_eta()->at(iel),b.Electron_phi()->at(iel),0.000511);
            ROOT::Math::PtEtaPhiMVector p2(b.Electron_pt()->at(iel2),b.Electron_eta()->at(iel2),b.Electron_phi()->at(iel2),0.000511);
            double this_mll = (p1+p2).M();
            if (fabs(this_mll-91.2)<fabs(m_ll-91.2)) {
              ll_p = p1+p2;
              m_ll = this_mll;
            }
          }
        }
      }
    }
    for (unsigned imu = 1; imu < b.Muon_pt()->size(); imu++) {
      if (Muon_sig_[imu]) {
        for (unsigned imu2 = 0; imu2 < imu; imu2++) {
          if (Muon_sig_[imu2]) {
            ROOT::Math::PtEtaPhiMVector p1(b.Muon_pt()->at(imu),b.Muon_eta()->at(imu),b.Muon_phi()->at(imu),0.106);
            ROOT::Math::PtEtaPhiMVector p2(b.Muon_pt()->at(imu2),b.Muon_eta()->at(imu2),b.Muon_phi()->at(imu2),0.106);
            double this_mll = (p1+p2).M();
            if (fabs(this_mll-91.2)<fabs(m_ll-91.2)) {
              ll_p = p1+p2; 
              m_ll = this_mll;
            }
          }
        }
      }
    }
    if (m_ll > 0) {
      for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
        if (ph_sig_[iph]) {
          ROOT::Math::PtEtaPhiMVector pph(b.Photon_pt()->at(iph),b.Photon_eta()->at(iph),b.Photon_phi()->at(iph),0.0);
          double this_mass = (pph+ll_p).M();
          if (fabs(this_mass-125.3)<fabs(m_llg-125.3)) {
            pt_llg = (pph+ll_p).Pt();
            m_llg = this_mass;
          }
        }
      }
    }
    return pt_llg/m_llg;
  });

  const NamedFunc delta_r_cut("delta_r_cut",[Electron_sig, Muon_sig, Photon_sig](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    std::vector<double> ph_sig_ = Photon_sig.GetVector(b);
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_sig_[iel]) {
        for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
          if (ph_sig_[iph]) {
            double this_deltar = deltaR(b.Electron_eta()->at(iel), b.Electron_phi()->at(iel), b.Photon_eta()->at(iph), b.Photon_phi()->at(iph));
            if (this_deltar < 0.4)
              return 0;
          }
        }
      }
    }
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      if (Muon_sig_[imu]) {
        for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
          if (ph_sig_[iph]) {
            double this_deltar = deltaR(b.Muon_eta()->at(imu), b.Muon_phi()->at(imu), b.Photon_eta()->at(iph), b.Photon_phi()->at(iph));
            if (this_deltar < 0.4)
              return 0;
          }
        }
      }
    }
    return 1;
  });

  const NamedFunc Photon_minDr("Photon_minDr",[Electron_sig, Muon_sig, Photon_sig](const Baby &b) -> NamedFunc::ScalarType{
    std::vector<double> Electron_sig_ = Electron_sig.GetVector(b);
    std::vector<double> Muon_sig_ = Muon_sig.GetVector(b);
    std::vector<double> ph_sig_ = Photon_sig.GetVector(b);
    double max_ph_pt = -1;
    double ph_eta = 0;
    double ph_phi = 0;
    double min_dr = 999;
    for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
      if (ph_sig_[iph]) {
        if (b.Photon_pt()->at(iph)>max_ph_pt) {
          max_ph_pt = b.Photon_pt()->at(iph);
          ph_eta = b.Photon_eta()->at(iph);
          ph_phi = b.Photon_phi()->at(iph);
        }
      }
    }
    if (max_ph_pt < 0) return -999;
    for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
      if (Electron_sig_[iel]) {
        double this_deltar = deltaR(b.Electron_eta()->at(iel), b.Electron_phi()->at(iel), ph_eta, ph_phi);
        if (this_deltar < min_dr)
          min_dr = this_deltar;
      }
    }
    for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
      if (Muon_sig_[imu]) {
        double this_deltar = deltaR(b.Muon_eta()->at(imu), b.Muon_phi()->at(imu), ph_eta, ph_phi);
        if (this_deltar < 0.4)
          min_dr = this_deltar;
      }
    }
    return min_dr;
  });


  //GenPhoton - properties of highest pt gen photon

  const NamedFunc GenPhoton_pt("GenPhoton_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    double max_ph_pt = 0;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      if (b.GenPart_pdgId()->at(imc)==22) {
        bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
        if (mc_statusFlags[0]||mc_statusFlags[8]) {
          if (b.GenPart_pt()->at(imc) > max_ph_pt)
            max_ph_pt = b.GenPart_pt()->at(imc);
        }
      }
    }
    return max_ph_pt;
  });

  const NamedFunc GenPhoton_abseta("GenPhoton_abseta",[](const Baby &b) -> NamedFunc::ScalarType{
    double max_ph_pt = 0;
    double ph_abseta = 0;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      if (b.GenPart_pdgId()->at(imc)==22) {
        bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
        if (mc_statusFlags[0]||mc_statusFlags[8]) {
          if (b.GenPart_pt()->at(imc) > max_ph_pt) {
            max_ph_pt = b.GenPart_pt()->at(imc);
            ph_abseta = abs(b.GenPart_eta()->at(imc));
          }
        }
      }
    }
    return ph_abseta;
  });

  const NamedFunc GenPhoton_minDr("GenPhoton_minDr",[](const Baby &b) -> NamedFunc::ScalarType{
    double max_ph_pt = 0;
    double ph_mindr = 999;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      if (b.GenPart_pdgId()->at(imc)==22) {
        bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
        if (mc_statusFlags[0]||mc_statusFlags[8]) {
          if (b.GenPart_pt()->at(imc) > max_ph_pt) {
            max_ph_pt = b.GenPart_pt()->at(imc);
            ph_mindr = 999;
            for (unsigned imc2 = 0; imc2 < b.GenPart_pdgId()->size(); imc2++) {
              if (abs(b.GenPart_pdgId()->at(imc2))==11 || abs(b.GenPart_pdgId()->at(imc2))==13) {
                bitset<15> mc2_statusFlags(b.GenPart_statusFlags()->at(imc2));
                if (mc2_statusFlags[13]&&(mc2_statusFlags[0]||mc2_statusFlags[8])) {
                  if (b.GenPart_pt()->at(imc2)>5) {
                    double dr = deltaR(b.GenPart_eta()->at(imc2),b.GenPart_phi()->at(imc2),b.GenPart_eta()->at(imc),b.GenPart_phi()->at(imc));
                    if (dr < ph_mindr)
                      ph_mindr = dr;
                  }
                }
              }
            }
          }
        }
      }
    }
    return ph_mindr;
  });

  //GenLepton - properties of generator level leptons

  const NamedFunc GenLepton_minpt("GenLepton_minpt",[](const Baby &b) -> NamedFunc::ScalarType{
    double max_lep_pt = -1;
    double max_antilep_pt = -1;
    for (int imc = 0; imc < b.nGenPart(); imc++) {
      if (b.GenPart_pdgId()->at(imc)==11 || b.GenPart_pdgId()->at(imc)==13) {
        bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
        if (mc_statusFlags[13]&&(mc_statusFlags[0]||mc_statusFlags[8])) {
          if (b.GenPart_pt()->at(imc)>max_lep_pt) {
            max_lep_pt = b.GenPart_pt()->at(imc);
          }
        }
      }
      if (b.GenPart_pdgId()->at(imc)==-11 || b.GenPart_pdgId()->at(imc)==-13) {
        bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
        if (mc_statusFlags[13]&&(mc_statusFlags[0]||mc_statusFlags[8])) {
          if (b.GenPart_pt()->at(imc)>max_antilep_pt) {
            max_antilep_pt = b.GenPart_pt()->at(imc);
          }
        }
      }
    }
    if (max_lep_pt > max_antilep_pt)
      return max_antilep_pt;
    return max_lep_pt;
  });

  const NamedFunc GenLepton_maxpt("GenLepton_maxpt",[](const Baby &b) -> NamedFunc::ScalarType{
    double max_lep_pt = -1;
    double max_antilep_pt = -1;
    for (int imc = 0; imc < b.nGenPart(); imc++) {
      if (b.GenPart_pdgId()->at(imc)==11 || b.GenPart_pdgId()->at(imc)==13) {
        bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
        if (mc_statusFlags[13]&&(mc_statusFlags[0]||mc_statusFlags[8])) {
          if (b.GenPart_pt()->at(imc)>max_lep_pt) {
            max_lep_pt = b.GenPart_pt()->at(imc);
          }
        }
      }
      if (b.GenPart_pdgId()->at(imc)==-11 || b.GenPart_pdgId()->at(imc)==-13) {
        bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
        if (mc_statusFlags[13]&&(mc_statusFlags[0]||mc_statusFlags[8])) {
          if (b.GenPart_pt()->at(imc)>max_antilep_pt) {
            max_antilep_pt = b.GenPart_pt()->at(imc);
          }
        }
      }
    }
    if (max_lep_pt > max_antilep_pt)
      return max_lep_pt;
    return max_antilep_pt;
  });

  //pdgid of highest pt lepton
  const NamedFunc GenLepton_pdgId("GenLepton_pdgId",[](const Baby &b) -> NamedFunc::ScalarType{
    double max_lep_pt = -1;
    int pdgid = 0;
    for (int imc = 0; imc < b.nGenPart(); imc++) {
      if (abs(b.GenPart_pdgId()->at(imc))==11 || abs(b.GenPart_pdgId()->at(imc))==13) {
        bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
        if (mc_statusFlags[13]&&(mc_statusFlags[0]||mc_statusFlags[8])) {
          if (b.GenPart_pt()->at(imc)>max_lep_pt) {
            max_lep_pt = b.GenPart_pt()->at(imc);
            pdgid = abs(b.GenPart_pdgId()->at(imc));
          }
        }
      }
    }
    return pdgid;
  });

  const NamedFunc GenLepton_etaacc("GenLepton_etaacc",[](const Baby &b) -> NamedFunc::ScalarType{
    float max_lep_pt = 0;
    float max_antilep_pt = 0;
    bool lep_acc = true;
    bool antilep_acc = true;
    for (int imc = 0; imc < b.nGenPart(); imc++) {
      int pdgid = b.GenPart_pdgId()->at(imc);
      if (pdgid==11 || pdgid==13) {
        bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
        if (mc_statusFlags[13]&&(mc_statusFlags[0]||mc_statusFlags[8])) {
          if (b.GenPart_pt()->at(imc)>max_lep_pt) {
            max_lep_pt = b.GenPart_pt()->at(imc);
            float abseta = abs(b.GenPart_eta()->at(imc));
            if (pdgid==11 && abseta>2.5) lep_acc = false;
            else if (pdgid==13 && abseta>2.4) lep_acc = false;
            else lep_acc = true;
          }
        }
      }
      if (pdgid==-11 || pdgid==-13) {
        bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
        if (mc_statusFlags[13]&&(mc_statusFlags[0]||mc_statusFlags[8])) {
          if (b.GenPart_pt()->at(imc)>max_antilep_pt) {
            max_antilep_pt = b.GenPart_pt()->at(imc);
            float abseta = abs(b.GenPart_eta()->at(imc));
            if (pdgid==-11 && abseta>2.5) antilep_acc = false;
            else if (pdgid==-13 && abseta>2.4) antilep_acc = false;
            else antilep_acc = true;
          }
        }
      }
    }
    return (antilep_acc && lep_acc);
  });

  //GenZCand - properties of gen ll system

  const NamedFunc GenZCand_mass("GenZCand_mass",[](const Baby &b) -> NamedFunc::ScalarType{
    double max_lep_pt = 0, max_antilep_pt = 0;
    TLorentzVector l1, l2;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      if (b.GenPart_pdgId()->at(imc)==11||b.GenPart_pdgId()->at(imc)==13) {
        bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
        if (mc_statusFlags[13]&&(mc_statusFlags[0]||mc_statusFlags[8])) {
          if (b.GenPart_pt()->at(imc) > max_lep_pt) {
            double m = GetMass(b.GenPart_pdgId()->at(imc));
            max_lep_pt = b.GenPart_pt()->at(imc);
            l1.SetPtEtaPhiM(b.GenPart_pt()->at(imc),b.GenPart_eta()->at(imc),b.GenPart_phi()->at(imc),m);
          }
        }
      }
      else if (b.GenPart_pdgId()->at(imc)==-11||b.GenPart_pdgId()->at(imc)==-13) {
        bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
        if (mc_statusFlags[13]&&(mc_statusFlags[0]||mc_statusFlags[8])) {
          if (b.GenPart_pt()->at(imc) > max_antilep_pt) {
            double m = GetMass(b.GenPart_pdgId()->at(imc));
            max_antilep_pt = b.GenPart_pt()->at(imc);
            l2.SetPtEtaPhiM(b.GenPart_pt()->at(imc),b.GenPart_eta()->at(imc),b.GenPart_phi()->at(imc),m);
          }
        }
      }
    }
    return (l1+l2).M();
  });

  const NamedFunc GenZCand_decayPdgId("GenZCand_decayPdgId",[](const Baby &b) -> NamedFunc::ScalarType{
    //use highest pt lepton
    double max_lep_pt = 0;
    int pdgid = 0;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
      if (abs_pdgid==11 || abs_pdgid==13 || abs_pdgid==15) {
        if (b.GenPart_pt()->at(imc) > max_lep_pt) {
          max_lep_pt = b.GenPart_pt()->at(imc);
          pdgid = abs_pdgid;
        }
      }
    }
    return static_cast<double>(pdgid);
  });

  //GenHiggsCand - properties of highest pt gen llg system

  const NamedFunc GenHiggsCand_mass("GenHiggsCand_mass",[](const Baby &b) -> NamedFunc::ScalarType{
    double max_ph_pt = 0, max_lep_pt = 0, max_antilep_pt = 0;
    int lep_flav = -1, antilep_flav = -2;
    TLorentzVector ph, l1, l2;
    for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
      if (b.GenPart_pdgId()->at(imc)==22) {
        bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
        if (mc_statusFlags[0]||mc_statusFlags[8]) {
          if (b.GenPart_pt()->at(imc) > max_ph_pt) {
            max_ph_pt = b.GenPart_pt()->at(imc);
            ph.SetPtEtaPhiM(b.GenPart_pt()->at(imc),b.GenPart_eta()->at(imc),b.GenPart_phi()->at(imc),0.0);
          }
        }
      }
      else if (b.GenPart_pdgId()->at(imc)==11||b.GenPart_pdgId()->at(imc)==13) {
        bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
        if (mc_statusFlags[13]&&(mc_statusFlags[0]||mc_statusFlags[8])) {
          if (b.GenPart_pt()->at(imc) > max_lep_pt) {
            double m = GetMass(b.GenPart_pdgId()->at(imc));
            max_lep_pt = b.GenPart_pt()->at(imc);
            lep_flav = abs(b.GenPart_pdgId()->at(imc));
            l1.SetPtEtaPhiM(b.GenPart_pt()->at(imc),b.GenPart_eta()->at(imc),b.GenPart_phi()->at(imc),m);
          }
        }
      }
      else if (b.GenPart_pdgId()->at(imc)==-11||b.GenPart_pdgId()->at(imc)==-13) {
        bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
        if (mc_statusFlags[13]&&(mc_statusFlags[0]||mc_statusFlags[8])) {
          if (b.GenPart_pt()->at(imc) > max_antilep_pt) {
            double m = GetMass(b.GenPart_pdgId()->at(imc));
            max_antilep_pt = b.GenPart_pt()->at(imc);
            antilep_flav = abs(b.GenPart_pdgId()->at(imc));
            l2.SetPtEtaPhiM(b.GenPart_pt()->at(imc),b.GenPart_eta()->at(imc),b.GenPart_phi()->at(imc),m);
          }
        }
      }
    }
    if (lep_flav != antilep_flav) return -1;
    return (ph+l1+l2).M();
  });

  const NamedFunc MinDeltaPt("MinDeltaPt",[](const Baby &b) -> NamedFunc::ScalarType{
    double ph_pt = b.photon_pt()->at(0);
    double min_dpt = 999.0;
    for (unsigned iel = 0; iel < b.el_pt()->size(); iel++) {
      if (!b.el_sig()->at(iel)) continue;
      double dpt = fabs(ph_pt-b.el_pt()->at(iel));
      if (dpt < min_dpt)
        min_dpt = dpt;
    }
    return min_dpt;
  });

  //plot_cutflow
  
  NamedFunc ObjectPreselection = (NElectron_sig>=2||NMuon_sig>=2)&&(Lead_Electron_pt>25||Lead_Muon_pt>20)&&(Sublead_Electron_pt>15||Sublead_Muon_pt>10)&&NPhoton_sig>=1&&(Lead_Photon_pt>15);
  ObjectPreselection.Name("Reco ee#gamma or #mu#mu#gamma");

  //NamedFunc baseline_selection = ((ZCand_mass>50)&&((Lead_Electron_pt>25)||(Lead_Muon_pt>20))&&(Lead_Photon_pt/HiggsCand_mass>0.14)&&delta_r_cut&&((HiggsCand_mass+ZCand_mass)>185));
  //NamedFunc baseline_selection_nopt = ((ZCand_mass>50)&&(Lead_Photon_pt/HiggsCand_mass>0.14)&&delta_r_cut&&((HiggsCand_mass+ZCand_mass)>185));
  //NamedFunc baseline_selection_nodr = ((ZCand_mass>50)&&((Lead_Electron_pt>25)||(Lead_Muon_pt>20))&&(Lead_Photon_pt/HiggsCand_mass>0.14)&&((HiggsCand_mass+ZCand_mass)>185));
  //NamedFunc midleadpt = ((Lead_Electron_pt>25)||(Lead_Muon_pt>20));
  //NamedFunc highleadpt = ((Lead_Electron_pt>38)||(Lead_Muon_pt>29));
  //NamedFunc highpt = ((Lead_Electron_pt>38&&Sublead_Electron_pt>38)||(Lead_Muon_pt>29&&Sublead_Muon_pt>29));
  //NamedFunc dy_selection = ((ZCand_mass>81)&&(ZCand_mass<101))&&(NElectron_sig>=2||NMuon_sig>=2);
  
  NamedFunc sublead_lep_ptcut = (GenLepton_pdgId==11&&GenLepton_minpt>15)||(GenLepton_pdgId==13&&GenLepton_minpt>10);
  sublead_lep_ptcut.Name("Sublead lepton p_{T} cut");
  NamedFunc lead_lep_ptcut = (GenLepton_pdgId==11&&GenLepton_maxpt>25)||(GenLepton_pdgId==13&&GenLepton_maxpt>20);
  lead_lep_ptcut.Name("Lead lepton p_{T} cut");

  const NamedFunc tighter_baseline = NamedFunc(ZgFunctions::zg_baseline&&"photon_id80[0]&&ll_m[0]>80&&ll_m[0]<100").Name("new_baseline");

  //------------------------------------------------------------------------------------
  //                                     make plots and pie charts
  //------------------------------------------------------------------------------------
  
  PlotMaker pm;

  bool plot_stitchcheck = false;
  bool plot_ptbins = false;
  bool plot_cutflow = false;
  bool plot_dysplit = true;

  if (plot_stitchcheck) {
    pm.Push<Hist1D>(Axis(40, 80, 160, GenHiggsCand_mass, "Truth m_{ll#gamma} [GeV]", {}), 
        "1", procs_dy_nostitch, plt_lin).Weight(weight).Tag("zgfit_nostitch");
    pm.Push<Hist1D>(Axis(30, 100, 160, GenHiggsCand_mass, "Truth m_{ll#gamma} [GeV]", {}), 
        "1", procs_dy_nostitch, plt_lin).Weight(weight).Tag("zgfit_nostitch_compressedaxis");
    pm.Push<Hist1D>(Axis(40, 0, 100, GenPhoton_pt, "Truth photon p_{T} [GeV]", {}), 
        "1", procs_dy_nostitch, plt_lin).Weight(weight).Tag("zgfit_nostitch");
    pm.Push<Hist1D>(Axis(40, 5, 100, GenPhoton_pt, "Truth photon p_{T} [GeV]", {}), 
        "1", procs_dy_nostitch, plt_lin).Weight(weight).Tag("zgfit_nostitch_compressedaxis");
    pm.Push<Hist1D>(Axis(40, 80, 160, GenHiggsCand_mass, "Truth m_{ll#gamma} [GeV]", {}), 
        "1", procs_dy, plt_lin).Weight(weight).Tag("zgfit");
    pm.Push<Hist1D>(Axis(30, 100, 160, GenHiggsCand_mass, "Truth m_{ll#gamma} [GeV]", {}), 
        "1", procs_dy, plt_lin).Weight(weight).Tag("zgfit_compressedaxis");
    pm.Push<Hist1D>(Axis(40, 0, 100, GenPhoton_pt, "Truth photon p_{T} [GeV]", {}), 
        "1", procs_dy, plt_lin).Weight(weight).Tag("zgfit");
    pm.Push<Hist1D>(Axis(40, 5, 100, GenPhoton_pt, "Truth photon p_{T} [GeV]", {}), 
        "1", procs_dy, plt_lin).Weight(weight).Tag("zgfit_compressedaxis");
  }

  if (plot_ptbins) {
    for (unsigned iaxis = 0; iaxis < 2; iaxis++) {
      int nbins = 40;
      int nbins_fake = 20;
      double lower_bound = 80;
      std::string tag = "zgfit";
      if (iaxis==1) {
        nbins = 30;
        nbins_fake = 15;
        lower_bound = 100;
        tag = "zgfit_compressedaxis";
      }
      std::vector<float> phpt_bins = {0,5,10,15,20,25,35,999};
      for (unsigned iphpt = 0; iphpt < phpt_bins.size()-1; iphpt++) {
        pm.Push<Hist1D>(Axis(nbins, lower_bound, 160, GenHiggsCand_mass, "Truth m_{ll#gamma} [GeV]", {}), 
            GenPhoton_pt>phpt_bins[iphpt]&&GenPhoton_pt<phpt_bins[iphpt+1], procs_dy, plt_lin).Weight(weight).Tag(tag);
        pm.Push<Hist1D>(Axis(nbins, lower_bound, 160, GenHiggsCand_mass, "Truth m_{ll#gamma} [GeV]", {}), 
            GenPhoton_pt>phpt_bins[iphpt], procs_dy, plt_lin).Weight(weight).Tag(tag);
        pm.Push<Hist1D>(Axis(nbins, lower_bound, 160, GenHiggsCand_mass, "Truth m_{ll#gamma} [GeV]", {}), 
            GenPhoton_pt>phpt_bins[iphpt]&&GenPhoton_pt<phpt_bins[iphpt+1], procs_dy_nostitch, plt_lin).Weight(weight).Tag(tag+"_nostitch");
        pm.Push<Hist1D>(Axis(nbins, lower_bound, 160, GenHiggsCand_mass, "Truth m_{ll#gamma} [GeV]", {}), 
            GenPhoton_pt>phpt_bins[iphpt], procs_dy_nostitch, plt_lin).Weight(weight).Tag(tag+"_nostitch");
        //fake photons
        pm.Push<Hist1D>(Axis(nbins_fake, lower_bound, 160, HiggsCand_mass, "Reco m_{ll#gamma} [GeV]", {}), 
            ObjectPreselection&&Photon_minDr>0.4&&Lead_Photon_pt>phpt_bins[iphpt], procs_photon_jet, plt_lin).Weight(weight).Tag(tag+"_photonjet");
        pm.Push<Hist1D>(Axis(nbins_fake, lower_bound, 160, HiggsCand_mass, "Reco m_{ll#gamma} [GeV]", {}), 
            ObjectPreselection&&Photon_minDr>0.4&&Lead_Photon_pt>phpt_bins[iphpt]&&Lead_Photon_pt<phpt_bins[iphpt+1], procs_photon_jet, plt_lin).Weight(weight).Tag(tag+"_photonjet");
        pm.Push<Hist1D>(Axis(nbins_fake, lower_bound, 160, HiggsCand_mass, "Reco m_{ll#gamma} [GeV]", {}), 
            ObjectPreselection&&Photon_minDr>0.4&&Lead_Photon_pt>phpt_bins[iphpt], procs_photon_other, plt_lin).Weight(weight).Tag(tag+"_photonother");
        pm.Push<Hist1D>(Axis(nbins_fake, lower_bound, 160, HiggsCand_mass, "Reco m_{ll#gamma} [GeV]", {}), 
            ObjectPreselection&&Photon_minDr>0.4&&Lead_Photon_pt>phpt_bins[iphpt]&&Lead_Photon_pt<phpt_bins[iphpt+1], procs_photon_other, plt_lin).Weight(weight).Tag(tag+"_photonother");
        pm.Push<Hist1D>(Axis(nbins_fake, lower_bound, 160, HiggsCand_mass, "Reco m_{ll#gamma} [GeV]", {}), 
            ObjectPreselection&&Photon_minDr>0.4&&Lead_Photon_pt>phpt_bins[iphpt], procs_photon_all, plt_lin).Weight(weight).Tag(tag+"_photonall");
        pm.Push<Hist1D>(Axis(nbins_fake, lower_bound, 160, HiggsCand_mass, "Reco m_{ll#gamma} [GeV]", {}), 
            ObjectPreselection&&Photon_minDr>0.4&&Lead_Photon_pt>phpt_bins[iphpt]&&Lead_Photon_pt<phpt_bins[iphpt+1], procs_photon_all, plt_lin).Weight(weight).Tag(tag+"_photonall");
      }
    }
  }

  if (plot_cutflow) {
    //various plots for debugging weirdness
    //turned out to be parent-less leptons
    //pm.Push<Hist1D>(Axis(40, 0, 100, GenLepton_minpt, "Truth Min Lepton p_{T} [GeV]", {}), 
    //    GenPhoton_pt>15&&GenPhoton_abseta<2.5, procs_dy, plt_lin).Weight(weight).Tag("zgfit");
    //pm.Push<Hist1D>(Axis(40, -20, 100, GenLepton_minpt, "Truth Min Lepton p_{T} [GeV]", {}), 
    //    GenPhoton_pt>15&&GenPhoton_abseta<2.5, procs_dy, plt_lin).Weight(weight).Tag("zgfit_extendedaxis");
    //pm.Push<Hist1D>(Axis(40, 0, 100, GenLepton_minpt, "Truth Min Lepton p_{T} [GeV]", {}), 
    //    GenPhoton_pt>15&&GenPhoton_abseta<2.5, procs_dy_nostitch, plt_lin).Weight(weight).Tag("zgfit_nostitch");
    //pm.Push<Hist2D>(Axis(40, 0, 100, GenLepton_minpt, "Truth Min Lepton p_{T} [GeV]", {}), 
    //    Axis(30, 100, 160, GenHiggsCand_mass, "Truth m_{ll#gamma} [GeV]", {}), 
    //    GenPhoton_pt>15&&GenPhoton_abseta<2.5&&GenPhoton_isFsr, procs_dy, twodim_log_plotopts).Weight(weight).Tag("zgfit");
    //pm.Push<Hist2D>(Axis(40, 0, 100, GenLepton_minpt, "Truth Min Lepton p_{T} [GeV]", {}), 
    //    Axis(30, 100, 160, GenHiggsCand_mass, "Truth m_{ll#gamma} [GeV]", {}), 
    //    GenPhoton_pt>15&&GenPhoton_abseta<2.5&&!GenPhoton_isFsr, procs_dy, twodim_log_plotopts).Weight(weight).Tag("zgfit");
    //pm.Push<Hist2D>(Axis(40, 0, 100, GenLepton_minpt, "Truth Min Lepton p_{T} [GeV]", {}), 
    //    Axis(30, 100, 160, GenHiggsCand_mass, "Truth m_{ll#gamma} [GeV]", {}), 
    //    GenPhoton_pt>15&&GenPhoton_abseta<2.5, procs_dy, twodim_log_plotopts).Weight(weight).Tag("zgfit");
    ////Z vs H mass
    //pm.Push<Hist2D>(Axis(30, 60, 120, GenZCand_mass, "Truth m_{ll} [GeV]", {}), 
    //    Axis(30, 100, 160, GenHiggsCand_mass, "Truth m_{ll#gamma} [GeV]", {}), 
    //    PicoType==6000&&!GenPhoton_isFsr, procs_dy, twodim_plotopts).Weight(weight).Tag("zgfit");
    //pm.Push<Hist2D>(Axis(30, 60, 120, GenZCand_mass, "Truth m_{ll} [GeV]", {}), 
    //    Axis(30, 100, 160, GenHiggsCand_mass, "Truth m_{ll#gamma} [GeV]", {}), 
    //    PicoType!=6000&&!GenPhoton_isFsr, procs_dy, twodim_plotopts).Weight(weight).Tag("zgfit");
    //pm.Push<Hist2D>(Axis(30, 60, 120, GenZCand_mass, "Truth m_{ll} [GeV]", {}), 
    //    Axis(30, 100, 160, GenHiggsCand_mass, "Truth m_{ll#gamma} [GeV]", {}), 
    //    PicoType==6000&&GenPhoton_isFsr, procs_dy, twodim_plotopts).Weight(weight).Tag("zgfit");
    //pm.Push<Hist2D>(Axis(30, 60, 120, GenZCand_mass, "Truth m_{ll} [GeV]", {}), 
    //    Axis(30, 100, 160, GenHiggsCand_mass, "Truth m_{ll#gamma} [GeV]", {}), 
    //    PicoType!=6000&&GenPhoton_isFsr, procs_dy, twodim_plotopts).Weight(weight).Tag("zgfit");
    //pm.Push<Hist2D>(Axis(30, 60, 120, GenZCand_mass, "Truth m_{ll} [GeV]", {}), 
    //    Axis(30, 100, 160, GenHiggsCand_mass, "Truth m_{ll#gamma} [GeV]", {}), 
    //    "1", procs_dy, twodim_plotopts).Weight(weight).Tag("zgfit");
    //pm.Push<Hist2D>(Axis(30, 60, 120, GenZCand_mass, "Truth m_{ll} [GeV]", {}), 
    //    Axis(30, 100, 160, GenHiggsCand_mass, "Truth m_{ll#gamma} [GeV]", {}), 
    //    PicoType==6000, procs_dy, twodim_plotopts).Weight(weight).Tag("zgfit");
    //pm.Push<Hist2D>(Axis(30, 60, 120, GenZCand_mass, "Truth m_{ll} [GeV]", {}), 
    //    Axis(30, 100, 160, GenHiggsCand_mass, "Truth m_{ll#gamma} [GeV]", {}), 
    //    PicoType!=6000, procs_dy, twodim_plotopts).Weight(weight).Tag("zgfit");
    //pm.Push<Hist2D>(Axis(30, 60, 120, GenZCand_mass, "Truth m_{ll} [GeV]", {}), 
    //    Axis(30, 100, 160, GenHiggsCand_mass, "Truth m_{ll#gamma} [GeV]", {}), 
    //    GenPhoton_isFsr, procs_dy, twodim_plotopts).Weight(weight).Tag("zgfit");
    //pm.Push<Hist2D>(Axis(30, 60, 120, GenZCand_mass, "Truth m_{ll} [GeV]", {}), 
    //    Axis(30, 100, 160, GenHiggsCand_mass, "Truth m_{ll#gamma} [GeV]", {}), 
    //    !GenPhoton_isFsr, procs_dy, twodim_plotopts).Weight(weight).Tag("zgfit");
    for (unsigned iaxis = 0; iaxis < 2; iaxis++) {
      int nbins = 40;
      int nbins_fake = 20;
      double lower_bound = 80;
      std::string base_tag = "zgfit";
      if (iaxis==1) {
        nbins = 30;
        nbins_fake = 15;
        lower_bound = 100;
        base_tag = "zgfit_compressedaxis";
      }

      for (unsigned lep_idx = 0; lep_idx < 3; lep_idx++) {
        NamedFunc lep_select = "1";
        NamedFunc reco_lep_select = "1";
        std::string tag = base_tag;
        if (lep_idx == 0) {
          lep_select = GenZCand_decayPdgId==11;
          reco_lep_select = NElectron_sig>=2;
          tag += "_el";
        }
        else if (lep_idx == 1) {
          lep_select = GenZCand_decayPdgId==13;
          reco_lep_select = NMuon_sig>=2;
          tag += "_mu";
        }

        pm.Push<Hist1D>(Axis(nbins, lower_bound, 160, GenHiggsCand_mass, "Truth m_{ll#gamma} [GeV]", {}), 
            lep_select, procs_dy, plt_lin).Weight(weight).Tag("FixName:"+tag+"__genmllg");
        pm.Push<Hist1D>(Axis(nbins, lower_bound, 160, GenHiggsCand_mass, "Truth m_{ll#gamma} [GeV]", {}), 
            lep_select&&GenPhoton_pt>15, procs_dy, plt_lin).Weight(weight).Tag("FixName:"+tag+"__genmllg_ptph");
        pm.Push<Hist1D>(Axis(nbins, lower_bound, 160, GenHiggsCand_mass, "Truth m_{ll#gamma} [GeV]", {}), 
            lep_select&&GenPhoton_pt>15&&GenPhoton_abseta<2.5, procs_dy, plt_lin).Weight(weight).Tag("FixName:"+tag+"__genmllg_ptph_etaph");
        pm.Push<Hist1D>(Axis(nbins, lower_bound, 160, GenHiggsCand_mass, "Truth m_{ll#gamma} [GeV]", {}), 
            lep_select&&GenPhoton_pt>15&&GenPhoton_abseta<2.5&&sublead_lep_ptcut, procs_dy, plt_lin).Weight(weight).Tag("FixName:"+tag+"__genmllg_ptph_etaph_slleppt");
        pm.Push<Hist1D>(Axis(nbins, lower_bound, 160, GenHiggsCand_mass, "Truth m_{ll#gamma} [GeV]", {}), 
            lep_select&&GenPhoton_pt>15&&GenPhoton_abseta<2.5&&sublead_lep_ptcut&&lead_lep_ptcut, procs_dy, plt_lin).Weight(weight).Tag("FixName:"+tag+"__genmllg_ptph_etaph_slleppt_lleppt");
        pm.Push<Hist1D>(Axis(nbins, lower_bound, 160, GenHiggsCand_mass, "Truth m_{ll#gamma} [GeV]", {}), 
            lep_select&&GenPhoton_pt>15&&GenPhoton_abseta<2.5&&sublead_lep_ptcut&&lead_lep_ptcut&&GenLepton_etaacc, procs_dy, plt_lin).Weight(weight).Tag("FixName:"+tag+"__genmllg_ptph_etaph_slleppt_lleppt_lepeta");
        pm.Push<Hist1D>(Axis(nbins, lower_bound, 160, GenHiggsCand_mass, "Truth m_{ll#gamma} [GeV]", {}), 
            lep_select&&GenPhoton_pt>15&&GenPhoton_abseta<2.5&&sublead_lep_ptcut&&lead_lep_ptcut&&GenLepton_etaacc&&GenPhoton_minDr>0.4, procs_dy, plt_lin).Weight(weight).Tag("FixName:"+tag+"__genmllg_ptph_etaph_slleppt_lleppt_lepeta_mindr");

        pm.Push<Hist1D>(Axis(nbins, lower_bound, 160, HiggsCand_mass, "Reco m_{ll#gamma} [GeV]", {}), 
            lep_select&&ObjectPreselection&&Photon_minDr>0.4&&Photon_truthMatch, procs_dy, plt_lin).Weight(weight).Tag("FixName:"+tag+"__recomllg_mindr");
        pm.Push<Hist1D>(Axis(nbins, lower_bound, 160, HiggsCand_mass, "Reco m_{ll#gamma} [GeV]", {}), 
            lep_select&&ObjectPreselection&&Photon_minDr>0.4&&Photon_truthMatch&&ZCand_mass>50, procs_dy, plt_lin).Weight(weight).Tag("FixName:"+tag+"__recomllg_mindr_mz");

        pm.Push<Hist1D>(Axis(nbins, lower_bound, 160, HiggsCand_mass, "Reco m_{ll#gamma} [GeV]", {}), 
            lep_select&&ObjectPreselection&&Photon_minDr>0.4&&Photon_truthMatch&&ZCand_mass>50&&(HiggsCand_mass+ZCand_mass>185), procs_dy, plt_lin).Weight(weight).Tag("FixName:"+tag+"__recomllg_mindr_mz_mzmh");
        pm.Push<Hist1D>(Axis(nbins, lower_bound, 160, HiggsCand_mass, "Reco m_{ll#gamma} [GeV]", {}), 
            lep_select&&ObjectPreselection&&Photon_minDr>0.4&&Photon_truthMatch&&ZCand_mass>50&&(HiggsCand_mass+ZCand_mass>185)&&(Lead_Photon_pt/HiggsCand_mass>0.14), procs_dy, plt_lin).Weight(weight).Tag("FixName:"+tag+"__recomllg_mindr_mz_mzmh_phpt");
        pm.Push<Hist1D>(Axis(nbins, lower_bound, 160, HiggsCand_mass, "Reco m_{ll#gamma} [GeV]", {}), 
            lep_select&&ObjectPreselection&&Photon_minDr>0.4&&Photon_truthMatch&&ZCand_mass>50&&(HiggsCand_mass+ZCand_mass>185)&&(Lead_Photon_pt/HiggsCand_mass>0.14)&&Lead_Photon_mvaID>0.5, procs_dy, plt_lin).Weight(weight).Tag("FixName:"+tag+"__recomllg_mindr_mz_mzmh_phpt_phid");
        //fake photons (jets)
        pm.Push<Hist1D>(Axis(nbins_fake, lower_bound, 160, HiggsCand_mass, "Reco m_{ll#gamma} [GeV]", {}), 
            reco_lep_select&&ObjectPreselection&&Photon_minDr>0.4, procs_photon_jet, plt_lin).Weight(weight).Tag("FixName:"+tag+"_photonjet__recomllg_mindr");
        pm.Push<Hist1D>(Axis(nbins_fake, lower_bound, 160, HiggsCand_mass, "Reco m_{ll#gamma} [GeV]", {}), 
            reco_lep_select&&ObjectPreselection&&Photon_minDr>0.4&&ZCand_mass>50, procs_photon_jet, plt_lin).Weight(weight).Tag("FixName:"+tag+"_photonjet__recomllg_mindr_mz");
        pm.Push<Hist1D>(Axis(nbins_fake, lower_bound, 160, HiggsCand_mass, "Reco m_{ll#gamma} [GeV]", {}), 
            reco_lep_select&&ObjectPreselection&&Photon_minDr>0.4&&ZCand_mass>50&&(HiggsCand_mass+ZCand_mass>185), procs_photon_jet, plt_lin).Weight(weight).Tag("FixName:"+tag+"_photonjet__recomllg_mindr_mz_mzmh");
        pm.Push<Hist1D>(Axis(nbins_fake, lower_bound, 160, HiggsCand_mass, "Reco m_{ll#gamma} [GeV]", {}), 
            reco_lep_select&&ObjectPreselection&&Photon_minDr>0.4&&ZCand_mass>50&&(HiggsCand_mass+ZCand_mass>185)&&(Lead_Photon_pt/HiggsCand_mass>0.14), procs_photon_jet, plt_lin).Weight(weight).Tag("FixName:"+tag+"_photonjet__recomllg_mindr_mz_mzmh_phpt");
        pm.Push<Hist1D>(Axis(nbins_fake, lower_bound, 160, HiggsCand_mass, "Reco m_{ll#gamma} [GeV]", {}), 
            reco_lep_select&&ObjectPreselection&&Photon_minDr>0.4&&ZCand_mass>50&&(HiggsCand_mass+ZCand_mass>185)&&(Lead_Photon_pt/HiggsCand_mass>0.14)&&Lead_Photon_mvaID>0.5, procs_photon_jet, plt_lin).Weight(weight).Tag("FixName:"+tag+"_photonjet__recomllg_mindr_mz_mzmh_phpt_phid");
        //fake photons (other)
        pm.Push<Hist1D>(Axis(nbins_fake, lower_bound, 160, HiggsCand_mass, "Reco m_{ll#gamma} [GeV]", {}), 
            reco_lep_select&&ObjectPreselection&&Photon_minDr>0.4, procs_photon_other, plt_lin).Weight(weight).Tag("FixName:"+tag+"_photonother__recomllg_mindr");
        pm.Push<Hist1D>(Axis(nbins_fake, lower_bound, 160, HiggsCand_mass, "Reco m_{ll#gamma} [GeV]", {}), 
            reco_lep_select&&ObjectPreselection&&Photon_minDr>0.4&&ZCand_mass>50, procs_photon_other, plt_lin).Weight(weight).Tag("FixName:"+tag+"_photonother__recomllg_mindr_mz");
        pm.Push<Hist1D>(Axis(nbins_fake, lower_bound, 160, HiggsCand_mass, "Reco m_{ll#gamma} [GeV]", {}), 
            reco_lep_select&&ObjectPreselection&&Photon_minDr>0.4&&ZCand_mass>50&&(HiggsCand_mass+ZCand_mass>185), procs_photon_other, plt_lin).Weight(weight).Tag("FixName:"+tag+"_photonother__recomllg_mindr_mz_mzmh");
        pm.Push<Hist1D>(Axis(nbins_fake, lower_bound, 160, HiggsCand_mass, "Reco m_{ll#gamma} [GeV]", {}), 
            reco_lep_select&&ObjectPreselection&&Photon_minDr>0.4&&ZCand_mass>50&&(HiggsCand_mass+ZCand_mass>185)&&(Lead_Photon_pt/HiggsCand_mass>0.14), procs_photon_other, plt_lin).Weight(weight).Tag("FixName:"+tag+"_photonother__recomllg_mindr_mz_mzmh_phpt");
        pm.Push<Hist1D>(Axis(nbins_fake, lower_bound, 160, HiggsCand_mass, "Reco m_{ll#gamma} [GeV]", {}), 
            reco_lep_select&&ObjectPreselection&&Photon_minDr>0.4&&ZCand_mass>50&&(HiggsCand_mass+ZCand_mass>185)&&(Lead_Photon_pt/HiggsCand_mass>0.14)&&Lead_Photon_mvaID>0.5, procs_photon_other, plt_lin).Weight(weight).Tag("FixName:"+tag+"_photonother__recomllg_mindr_mz_mzmh_phpt_phid");
        //all processes
        pm.Push<Hist1D>(Axis(nbins_fake, lower_bound, 160, HiggsCand_mass, "Reco m_{ll#gamma} [GeV]", {}), 
            reco_lep_select&&ObjectPreselection&&Photon_minDr>0.4, procs_photon_all, plt_lin).Weight(weight).Tag("FixName:"+tag+"_photonall__recomllg_mindr");
        pm.Push<Hist1D>(Axis(nbins_fake, lower_bound, 160, HiggsCand_mass, "Reco m_{ll#gamma} [GeV]", {}), 
            reco_lep_select&&ObjectPreselection&&Photon_minDr>0.4&&ZCand_mass>50, procs_photon_all, plt_lin).Weight(weight).Tag("FixName:"+tag+"_photonall__recomllg_mindr_mz");
        pm.Push<Hist1D>(Axis(nbins_fake, lower_bound, 160, HiggsCand_mass, "Reco m_{ll#gamma} [GeV]", {}), 
            reco_lep_select&&ObjectPreselection&&Photon_minDr>0.4&&ZCand_mass>50&&(HiggsCand_mass+ZCand_mass>185), procs_photon_all, plt_lin).Weight(weight).Tag("FixName:"+tag+"_photonall__recomllg_mindr_mz_mzmh");
        pm.Push<Hist1D>(Axis(nbins_fake, lower_bound, 160, HiggsCand_mass, "Reco m_{ll#gamma} [GeV]", {}), 
            reco_lep_select&&ObjectPreselection&&Photon_minDr>0.4&&ZCand_mass>50&&(HiggsCand_mass+ZCand_mass>185)&&(Lead_Photon_pt/HiggsCand_mass>0.14), procs_photon_all, plt_lin).Weight(weight).Tag("FixName:"+tag+"_photonall__recomllg_mindr_mz_mzmh_phpt");
        pm.Push<Hist1D>(Axis(nbins_fake, lower_bound, 160, HiggsCand_mass, "Reco m_{ll#gamma} [GeV]", {}), 
            reco_lep_select&&ObjectPreselection&&Photon_minDr>0.4&&ZCand_mass>50&&(HiggsCand_mass+ZCand_mass>185)&&(Lead_Photon_pt/HiggsCand_mass>0.14)&&Lead_Photon_mvaID>0.5, procs_photon_all, plt_lin).Weight(weight).Tag("FixName:"+tag+"_photonall__recomllg_mindr_mz_mzmh_phpt_phid");
      }
    }
  }

  if (plot_dysplit) {
    //pm.Push<Hist1D>(Axis(60.0,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
    //    tighter_baseline, procs_dysplit, plt_lin)
    //    .Weight("w_lumi"*ZgFunctions::w_years).Tag("zgshapes");
    //pm.Push<Hist1D>(Axis(60.0,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
    //    tighter_baseline&&"photon_elidx!=-1", procs_dysplit, plt_lin)
    //    .Weight("w_lumi"*ZgFunctions::w_years).Tag("zgshapes");
    //pm.Push<Hist1D>(Axis(60.0,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
    //    tighter_baseline&&overlapping_phel, procs_dysplit, plt_lin)
    //    .Weight("w_lumi"*ZgFunctions::w_years).Tag("zgshapes");
    //pm.Push<Hist1D>(Axis(60.0,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
    //    tighter_baseline&&!overlapping_phel, procs_dysplit, plt_lin)
    //    .Weight("w_lumi"*ZgFunctions::w_years).Tag("zgshapes");
    pm.Push<Hist1D>(Axis(40,100.0,180.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&overlapping_phel&&"photon_pflavor[0]==11", procs_dysplit, plt_lin)
        .Weight("w_lumi"*ZgFunctions::w_years).Tag("zgshapes");
    pm.Push<Hist1D>(Axis(40,100.0,180.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&!overlapping_phel&&"photon_pflavor[0]==11", procs_dysplit, plt_lin)
        .Weight("w_lumi"*ZgFunctions::w_years).Tag("zgshapes");
    pm.Push<Hist1D>(Axis(25,0.4,3.1, "photon_drmin[0]", "Min #Delta R(l, #gamma)", {}), 
        tighter_baseline&&overlapping_phel&&"photon_pflavor[0]==11", procs_dysplit, plt_lin)
        .Weight("w_lumi"*ZgFunctions::w_years).Tag("zgshapes");
    pm.Push<Hist1D>(Axis(25,0.4,3.1, "photon_drmin[0]", "Min #Delta R(l, #gamma)", {}), 
        tighter_baseline&&!overlapping_phel&&"photon_pflavor[0]==11", procs_dysplit, plt_lin)
        .Weight("w_lumi"*ZgFunctions::w_years).Tag("zgshapes");
    pm.Push<Hist1D>(Axis(25,0.0,50.0, MinDeltaPt, "Min #Delta p_{T}(l, #gamma) [GeV]", {}), 
        tighter_baseline&&overlapping_phel&&"photon_pflavor[0]==11", procs_dysplit, plt_lin)
        .Weight("w_lumi"*ZgFunctions::w_years).Tag("zgshapes");
    pm.Push<Hist1D>(Axis(25,0.0,50.0, MinDeltaPt, "Min #Delta p_{T}(l,#gamma) [GeV]", {}), 
        tighter_baseline&&!overlapping_phel&&"photon_pflavor[0]==11", procs_dysplit, plt_lin)
        .Weight("w_lumi"*ZgFunctions::w_years).Tag("zgshapes");
    pm.Push<Hist2D>(Axis(20,0.4,3.1, "photon_drmin[0]", "Min #Delta R(l, #gamma)", {}), 
        Axis(20,100.0,180.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&overlapping_phel&&"photon_pflavor[0]==11", procs_dysplit, twodim_log_plotopts)
        .Weight("w_lumi"*ZgFunctions::w_years).Tag("zgshapes");
    pm.Push<Hist2D>(Axis(20,0.4,3.1, "photon_drmin[0]", "Min #Delta R(l, #gamma)", {}), 
        Axis(20,100.0,180.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&!overlapping_phel&&"photon_pflavor[0]==11", procs_dysplit, twodim_log_plotopts)
        .Weight("w_lumi"*ZgFunctions::w_years).Tag("zgshapes");
    pm.Push<Hist2D>(Axis(25,0.0,50.0, MinDeltaPt, "Min #Delta p_{T}(l,#gamma) [GeV]", {}), 
        Axis(20,100.0,180.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&overlapping_phel&&"photon_pflavor[0]==11", procs_dysplit, twodim_log_plotopts)
        .Weight("w_lumi"*ZgFunctions::w_years).Tag("zgshapes");
    pm.Push<Hist2D>(Axis(25,0.0,50.0, MinDeltaPt, "Min #Delta p_{T}(l,#gamma) [GeV]", {}), 
        Axis(20,100.0,180.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&!overlapping_phel&&"photon_pflavor[0]==11", procs_dysplit, twodim_log_plotopts)
        .Weight("w_lumi"*ZgFunctions::w_years).Tag("zgshapes");
  }

  pm.multithreaded_ = !options.single_thread;
  pm.min_print_ = true;
  //pm.MakePlots(340.0);
  pm.MakePlots(1.0);

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
