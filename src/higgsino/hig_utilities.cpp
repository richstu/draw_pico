#include <algorithm>
#include <stdlib.h>
#include <regex>
#include "higgsino/hig_utilities.hpp"
#include "core/hist1d.hpp"
#include "core/utilities.hpp"
#include "core/palette.hpp"
#include "core/plot_opt.hpp"

namespace HigUtilities {
  using std::string;
  using std::to_string;
  using std::map;
  using std::vector;
  using std::set;
  using std::pair;
  using std::cout;
  using std::endl;
  using std::shared_ptr;
  using std::max;
  using std::make_pair;


  //HistInformation::HistInformation(Axis axis, NamedFunc cut, PlotOpt plot_opt) :
  //    axis_(axis), cut_(cut), plot_opt_(plot_opt) {
  //  figure_index = -1;
  //}

  int stringToVectorString(std::string const& inString, std::vector<std::string>& outputVector, std::string const & delimiter)
  {
      if( outputVector.size() != 0) {
      std::cout<<"[Error] StringFunctions::divideString() => outputVector size is not 0."<<std::endl;
      return 1;
    }
  
    if(delimiter.size()==0){
      std::cout<<"[Error] StringFunctions::divideString() => delimiter size is 0."<<std::endl;
      return 1;
    }
  
    size_t start = 0;
    size_t end = 0;
    while((end=inString.find(delimiter,start)) != string::npos) {
      if(start!=end) {
        //std::cout<<"found: "<<removeSpaces(inString.substr(start, end-start))<<std::endl;
        outputVector.push_back(removeSpaces(inString.substr(start, end-start)));
      }
      start = end + delimiter.length();
    }
    if(start!=inString.length()) {
      //std::cout<<"found: "<<removeSpaces(inString.substr(start))<<std::endl;
      outputVector.push_back(removeSpaces(inString.substr(start)));
    }
    return 0;
  }
  bool is_in_string_options(std::string const & string_options, std::string const & option) {
    std::vector<std::string> result;
    stringToVectorString(string_options, result, ",");
    return result.end() != std::find(result.begin(), result.end(), option);
  }
  
  int vectorStringToString(std::vector<std::string> const & inVector, std::string &outString, std::string const & delimiter){
    if (inVector.size()==0) {
      outString = "";
    } else {
      outString = "";
      for(unsigned iVector=0; iVector<inVector.size()-1; iVector++){
        outString += inVector[iVector] + delimiter;
      }
      outString += inVector[inVector.size()-1];
    }
    return 0;
  }
  
  std::string removeSpaces(std::string inString){
    inString.erase(std::remove(inString.begin(),inString.end(),' '),inString.end());
    return inString;
  }

  std::string nom2sys_string(std::string nom_string, std::string sys_idx) {
    std::vector<std::string> target_strings = {"met","ht","njet","nbl","nbm","nbt",
        "low_dphi_met","hig_cand_am[0]","hig_cand_dm[0]","hig_cand_drmax[0]"};
    std::vector<std::string> replacement_strings = {"sys_met["+sys_idx+"]",
        "sys_ht["+sys_idx+"]","sys_njet["+sys_idx+"]","sys_nbl["+sys_idx+"]",
        "sys_nbm["+sys_idx+"]","sys_nbt["+sys_idx+"]",
        "sys_low_dphi_met["+sys_idx+"]","sys_hig_cand_am["+sys_idx+"]",
        "sys_hig_cand_dm["+sys_idx+"]","sys_hig_cand_drmax["+sys_idx+"]"};
    std::string sys_string = nom_string;
    //std::cout << "DEBUG: sys_string: " << sys_string << std::endl;
    for (unsigned repl_idx = 0; repl_idx<target_strings.size(); repl_idx++) {
      //std::cout << "DEBUG: target_string: " << target_strings[repl_idx] << std::endl;
      size_t start_pos = sys_string.find(target_strings[repl_idx]); 
      while (start_pos != std::string::npos) {
        //int dummy;
        //std::cout << "DEBUG: start_pos: " << start_pos << std::endl;
        //std::cin >> dummy;
        //check for problematic substrings: m[ht], weig[ht], [met]_calo, low_dphi_[met], [met]_tru
        if (start_pos != 0) {
          if (target_strings[repl_idx] == "ht" && sys_string.at(start_pos-1)=='m') {
            start_pos = sys_string.find(target_strings[repl_idx], start_pos+target_strings[repl_idx].size()); 
            continue;
          }
          if (target_strings[repl_idx] == "ht" && sys_string.at(start_pos-1)=='g') {
            start_pos = sys_string.find(target_strings[repl_idx], start_pos+target_strings[repl_idx].size()); 
            continue;
          }
          if (target_strings[repl_idx] == "met" && sys_string.at(start_pos-1)=='_') {
            start_pos = sys_string.find(target_strings[repl_idx], start_pos+target_strings[repl_idx].size()); 
            continue;
          }
          if (target_strings[repl_idx] == "met" && sys_string.at(start_pos-1)=='_') {
            start_pos = sys_string.find(target_strings[repl_idx], start_pos+target_strings[repl_idx].size()); 
            continue;
          }
        }
        if (start_pos+target_strings[repl_idx].size() < sys_string.size()) {
          if (target_strings[repl_idx] == "met" && sys_string.at(start_pos+target_strings[repl_idx].size())=='_') {
            start_pos = sys_string.find(target_strings[repl_idx], start_pos+target_strings[repl_idx].size()); 
            continue;
          }
        }
        sys_string.replace(start_pos, target_strings[repl_idx].length(),
            replacement_strings[repl_idx]);
        start_pos = sys_string.find(target_strings[repl_idx], start_pos+replacement_strings[repl_idx].size()); 
      }
    }
    return sys_string;
  }

  std::vector<std::pair<std::string, std::string>> nom2sys_bins(
      std::vector<std::pair<std::string, std::string>> sample_bins, std::string sys_idx) {
    std::vector<std::pair<std::string, std::string>> sys_bins;
    for (unsigned bin_idx = 0; bin_idx < sample_bins.size(); bin_idx++) {
      sys_bins.push_back(std::pair<std::string, std::string>(sample_bins[bin_idx].first, nom2sys_string(sample_bins[bin_idx].second, sys_idx)));
    }
    return sys_bins;
  }


  // based on pass_run2 with removal of non-existing branches
  const NamedFunc pass_2016("pass_2016", [](const Baby &b) -> NamedFunc::ScalarType{
    bool pass_ =  b.pass_muon_jet() && (b.met()/b.met_calo()<5);
    if (b.type()<1000 && b.type()>0) { // Data
      pass_ = pass_ && b.pass_jets() && b.pass_goodv() && b.pass_hbhe() && 
              b.pass_hbheiso() && b.pass_ecaldeadcell() && b.pass_badpfmu() && b.pass_eebadsc();
    } else {
      if (b.type()>=100e3) { // FastSim
        //pass_ = pass_ && b.pass_fsjets() && b.pass_goodv() && b.pass_hbhe() && 
        //        b.pass_hbheiso() && b.pass_ecaldeadcell() && b.pass_badpfmu();
        pass_ = pass_ && b.pass_jets() && b.pass_goodv() && b.pass_hbhe() && 
                b.pass_hbheiso() && b.pass_ecaldeadcell() && b.pass_badpfmu();
      } else { //FullSim
        pass_ = pass_ && b.pass_jets() && b.pass_goodv() && b.pass_hbhe() && 
                b.pass_hbheiso() && b.pass_ecaldeadcell() && b.pass_badpfmu();
      }
    }
    if (pass_) return 1.;
    else return 0.;
  });

  // based on r136 of https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
  const NamedFunc pass_run2("pass_run2", [](const Baby &b) -> NamedFunc::ScalarType{
    bool pass_ = b.pass_muon_jet() && (b.met()/b.met_calo()<5);
    if (b.SampleType()<0) { // Data
      pass_ = pass_ && b.pass_jets() && b.pass_goodv() && b.pass_cschalo_tight() && b.pass_hbhe() && 
              b.pass_hbheiso() && b.pass_ecaldeadcell() && b.pass_badpfmu() && b.pass_eebadsc();
      if (b.SampleType()!=-2016) pass_ = pass_ && b.pass_badcalib();
    } else {
      if (b.type()>=100e3) { // FastSim
        //pass_ = pass_ && b.pass_fsjets() && b.pass_goodv() && b.pass_hbhe() && 
        //        b.pass_hbheiso() && b.pass_ecaldeadcell() && b.pass_badpfmu();
        pass_ = pass_ && b.pass_jets() && b.pass_goodv() && b.pass_hbhe() && 
                b.pass_hbheiso() && b.pass_ecaldeadcell() && b.pass_badpfmu();
        if (b.SampleType()!=2016) pass_ = pass_ && b.pass_badcalib();
      } else { //FullSim
        pass_ = pass_ && b.pass_jets() && b.pass_goodv() && b.pass_cschalo_tight() && b.pass_hbhe() && 
                b.pass_hbheiso() && b.pass_ecaldeadcell() && b.pass_badpfmu();
        if (b.SampleType()!=2016) pass_ = pass_ && b.pass_badcalib();
      }
    }
    if (pass_) return 1.;
    else return 0.;
  });
  
  const NamedFunc weight_2016("weight_2016", [](const Baby &b) -> NamedFunc::ScalarType{
    // Data
    if (b.type()<1000 && b.type()>0) return 1.;
    else return b.weight()/b.w_btag()*b.w_bhig();
  });

  //const NamedFunc w_years("w_years", [](const Baby &b) -> NamedFunc::ScalarType{
  //  if (b.SampleType()<0) return 1.;

  //  double wgt = 1;
  //  if (b.type()==106000) {
  //    return 137.;
  //  }
  //  if (b.SampleType()==2016){
  //    return wgt*35.9;
  //  } else if (b.SampleType()==2017){
  //    return wgt*41.5;
  //  } else {
  //    return wgt*59.6;
  //  }
  //});
  
  const NamedFunc pass_nhig_cand("pass_nhig_cand", [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.hig_cand_drmax()->size() == 0) return 0;
    else return 1;
  });

  const NamedFunc pass_nhig_df_cand("pass_nhig_df_cand", [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.hig_df_cand_drmax()->size() == 0) return 0;
    else return 1;
  });

  const NamedFunc w_CNToN1N2("w_CNToN1N2", [](const Baby &b) -> NamedFunc::ScalarType{
    if(b.type() != 106000) return 1;
    int mprod_rounded = ((b.mprod()+12)/25)*25;
    if(mprod_rounded ==125) return 1.44725/7.6022;
    else if(mprod_rounded ==150) return 0.71514/3.83231;
    else if(mprod_rounded ==175) return 0.419059/2.26794;
    else if(mprod_rounded ==200) return 0.244213/1.33562;
    else if(mprod_rounded ==225) return 0.156286/0.860597;
    else if(mprod_rounded ==250) return 0.104252/0.577314;
    else if(mprod_rounded ==275) return 0.0719125/0.400107;
    else if(mprod_rounded ==300) return 0.0509994/0.284855;
    else if(mprod_rounded ==325) return 0.0369715/0.20736;
    else if(mprod_rounded ==350) return 0.0273286/0.153841;
    else if(mprod_rounded ==375) return 0.0205429/0.116006;
    else if(mprod_rounded ==400) return 0.0156691/0.0887325;
    else if(mprod_rounded ==425) return 0.0120965/0.0686963;
    else if(mprod_rounded ==450) return 0.00944017/0.0537702;
    else if(mprod_rounded ==475) return 0.00743587/0.0424699;
    else if(mprod_rounded ==500) return 0.00590757/0.0338387;
    else if(mprod_rounded ==525) return 0.00473235/0.0271867;
    else if(mprod_rounded ==550) return 0.0038167/0.0219868;
    else if(mprod_rounded ==575) return 0.00309847/0.0179062;
    else if(mprod_rounded ==600) return 0.00253015/0.0146677;
    else if(mprod_rounded ==625) return 0.00207755/0.012062;
    else if(mprod_rounded ==650) return 0.00171418/0.00996406;
    else if(mprod_rounded ==675) return 0.0014199/0.00828246;
    else if(mprod_rounded ==700) return 0.00118113/0.00689981;
    else if(mprod_rounded ==725) return 0.00098639/0.00578355;
    else if(mprod_rounded ==750) return 0.000826366/0.0048731;
    else if(mprod_rounded ==775) return 0.000694985/0.00409781;
    else if(mprod_rounded ==800) return 0.000586211/0.00346143;
    else if(mprod_rounded ==825) return 0.000495914/0.0029337;
    else if(mprod_rounded ==850) return 0.000420556/0.0024923;
    else if(mprod_rounded ==875) return 0.000361029/0.00213679;
    else if(mprod_rounded ==900) return 0.000305935/0.00180616;
    else if(mprod_rounded ==925) return 0.000262621/0.00155453;
    else if(mprod_rounded ==950) return 0.00022285/0.00132692;
    else if(mprod_rounded ==975) return 0.0001909/0.00112975;
    else if(mprod_rounded ==1000) return 0.00016428/0.000968853;
    else if(mprod_rounded ==1025) return 0.00014139/0.000840602;
    else if(mprod_rounded ==1050) return 0.000121865/0.000731306;
    else if(mprod_rounded ==1075) return 0.000105913/0.000627083;
    else if(mprod_rounded ==1100) return 9.12469e-05/0.000538005;
    else if(mprod_rounded ==1125) return 7.93058e-05/0.00046747;
    else if(mprod_rounded ==1150) return 6.84561e-05/0.000405108;
    else if(mprod_rounded ==1175) return 5.93602e-05/0.000348261;
    else if(mprod_rounded ==1200) return 5.16263e-05/0.000299347;
    else if(mprod_rounded ==1225) return 4.4906e-05/0.000265935;
    else if(mprod_rounded ==1250) return 3.91587e-05/0.000240471;
    else if(mprod_rounded ==1275) return 3.43135e-05/0.000190411;
    else if(mprod_rounded ==1300) return 2.99353e-05/0.000160765;
    else if(mprod_rounded ==1325) return 2.62223e-05/0.000136272;
    else if(mprod_rounded ==1350) return 2.28072e-05/0.000111174;
    else if(mprod_rounded ==1375) return 2.00393e-05/9.74728e-05;
    else if(mprod_rounded ==1400) return 1.75031e-05/7.80263e-05;
    else if(mprod_rounded ==1425) return 1.53144e-05/6.96843e-05;
    else if(mprod_rounded ==1450) return 1.34572e-05/6.96962e-05;
    else if(mprod_rounded ==1475) return 1.17047e-05/4.98006e-05;
    else return 0;
  });

  //const NamedFunc w_CNToN1N2("w_CNToN1N2", [](const Baby &b) -> NamedFunc::ScalarType{
  //  if(b.type() != 106000) return 1;
  //  if(b.mprod() ==127) return 1.44725/7.6022;
  //  if(b.mprod() ==150) return 0.71514/3.83231;
  //  if(b.mprod() ==175) return 0.419059/2.26794;
  //  if(b.mprod() ==200) return 0.244213/1.33562;
  //  if(b.mprod() ==225) return 0.156286/0.860597;
  //  if(b.mprod() ==250) return 0.104252/0.577314;
  //  if(b.mprod() ==275) return 0.0719125/0.400107;
  //  if(b.mprod() ==300) return 0.0509994/0.284855;
  //  if(b.mprod() ==325) return 0.0369715/0.20736;
  //  if(b.mprod() ==350) return 0.0273286/0.153841;
  //  if(b.mprod() ==375) return 0.0205429/0.116006;
  //  if(b.mprod() ==400) return 0.0156691/0.0887325;
  //  if(b.mprod() ==425) return 0.0120965/0.0686963;
  //  if(b.mprod() ==450) return 0.00944017/0.0537702;
  //  if(b.mprod() ==475) return 0.00743587/0.0424699;
  //  if(b.mprod() ==500) return 0.00590757/0.0338387;
  //  if(b.mprod() ==526) return 0.00469101/0.0269524;
  //  if(b.mprod() ==550) return 0.0038167/0.0219868;
  //  if(b.mprod() ==576) return 0.003073/0.0177611;
  //  if(b.mprod() ==600) return 0.00253015/0.0146677;
  //  if(b.mprod() ==626) return 0.00206136/0.0119691;
  //  if(b.mprod() ==650) return 0.00171418/0.00996406;
  //  if(b.mprod() ==676) return 0.00140934/0.00822165;
  //  if(b.mprod() ==700) return 0.00118113/0.00689981;
  //  if(b.mprod() ==726) return 0.000979349/0.00574289;
  //  if(b.mprod() ==750) return 0.000826366/0.0048731;
  //  if(b.mprod() ==776) return 0.000690208/0.00407012;
  //  if(b.mprod() ==800) return 0.000586211/0.00346143;
  //  if(b.mprod() ==826) return 0.00049277/0.00291514;
  //  if(b.mprod() ==850) return 0.000420556/0.0024923;
  //  if(b.mprod() ==876) return 0.000358734/0.00212322;
  //  if(b.mprod() ==900) return 0.000305935/0.00180616;
  //  if(b.mprod() ==926) return 0.000260948/0.00154462;
  //  if(b.mprod() ==950) return 0.00022285/0.00132692;
  //  if(b.mprod() ==976) return 0.000189681/0.00112253;
  //  if(b.mprod() ==1000) return 0.00016428/0.000968853;
  //  if(b.mprod() ==1024) return 0.000142206/0.000845522;
  //  if(b.mprod() ==1052) return 0.000120971/0.000717628;
  //  if(b.mprod() ==1076) return 0.000105301/0.000623403;
  //  if(b.mprod() ==1100) return 9.12469e-05/0.000538005;
  //  if(b.mprod() ==1124) return 7.9765e-05/0.00047022;
  //  if(b.mprod() ==1152) return 6.78234e-05/0.000398787;
  //  if(b.mprod() ==1176) return 5.9016e-05/0.000346209;
  //  if(b.mprod() ==1200) return 5.16263e-05/0.000299347;
  //  if(b.mprod() ==1224) return 4.5147e-05/0.000267704;
  //  if(b.mprod() ==1252) return 3.88343e-05/0.000222061;
  //  if(b.mprod() ==1276) return 3.41304e-05/0.00018915;
  //  if(b.mprod() ==1300) return 2.99353e-05/0.000160765;
  //  if(b.mprod() ==1324) return 2.63637e-05/0.000137188;
  //  if(b.mprod() ==1352) return 2.26779e-05/0.000113724;
  //  if(b.mprod() ==1376) return 1.99318e-05/9.68213e-05;
  //  if(b.mprod() ==1400) return 1.75031e-05/7.80263e-05;
  //  if(b.mprod() ==1424) return 1.53974e-05/7.0157e-05;
  //  if(b.mprod() ==1452) return 1.3245e-05/5.81247e-05;
  //  if(b.mprod() ==1476) return 1.16416e-05/4.94646e-05;
  //  else return 0;
  //});


  TString nom2genmet(TString ibin){
    ibin.ReplaceAll("met", "met_tru");
    //fix unintended replacement...
    ibin.ReplaceAll("sys_met_tru[0]", "met_tru");
    ibin.ReplaceAll("sys_met_tru[1]", "met_tru");
    ibin.ReplaceAll("sys_met_tru[2]", "met_tru");
    ibin.ReplaceAll("sys_met_tru[3]", "met_tru");
    ibin.ReplaceAll("met_tru/met_tru_calo", "met/met_calo");
    ibin.ReplaceAll("low_dphi_met_tru", "low_dphi_met");
    ibin.ReplaceAll("sys_met_tru[0]/met_tru_calo", "sys_met[0]/met_calo");
    ibin.ReplaceAll("sys_met_tru[1]/met_tru_calo", "sys_met[1]/met_calo");
    ibin.ReplaceAll("sys_met_tru[2]/met_tru_calo", "sys_met[2]/met_calo");
    ibin.ReplaceAll("sys_met_tru[3]/met_tru_calo", "sys_met[3]/met_calo");
    ibin.ReplaceAll("low_dphi_met_tru", "low_dphi_met");
    return ibin;
  }
  
  TString nom2sys_bin(TString ibin, size_t shift_index){
    ibin.ReplaceAll("met", "sys_met["+to_string(shift_index)+"]");
    ibin.ReplaceAll("mt", "sys_mt["+to_string(shift_index)+"]");
    ibin.ReplaceAll("st", "sys_st["+to_string(shift_index)+"]");
    ibin.ReplaceAll("mj14", "sys_mj14["+to_string(shift_index)+"]");
    ibin.ReplaceAll("njets", "sys_njets["+to_string(shift_index)+"]");
    ibin.ReplaceAll("nbdm", "sys_nbdm["+to_string(shift_index)+"]");
    return ibin;
  }

  string getBaseFolder(string in_base_folder)
  {
    string baseFolder = "";
    string hostName = execute("echo $HOSTNAME");
    if(Contains(hostName, "cms") || Contains(hostName, "compute-")) baseFolder = in_base_folder;
    return baseFolder;
  }

  void parseMassPoints(string mass_points_string, vector<pair<string, string> > & mass_points)
  {
    if (mass_points_string!="") {
      size_t found = mass_points_string.find(",");
      size_t start = 0;
      string _tmp = "";
      while (found != string::npos) {
        _tmp = mass_points_string.substr(start, found-start);
        mass_points.push_back(make_pair(_tmp.substr(0,_tmp.find("_")), _tmp.substr(_tmp.find("_")+1)));
        cout<<"Adding mass point: mgluino = "<<mass_points.back().first<<" mlsp = "<<mass_points.back().second<<endl;
        start = found+1;
        found = mass_points_string.find(",",start);
      } 
      _tmp = mass_points_string.substr(start, found-start);
      mass_points.push_back(make_pair(_tmp.substr(0,_tmp.find("_")), _tmp.substr(_tmp.find("_")+1)));
      cout<<"Adding mass point: mgluino = "<<mass_points.back().first<<" mlsp = "<<mass_points.back().second<<endl;
    } else {
      mass_points.push_back(make_pair("127","1"));
      for (int ichi(150); ichi < 1301; ichi +=25) {
        mass_points.push_back(make_pair(to_string(ichi),"1"));
      }
    }
  }
  
  void parseYears(string years_string, set<int> & years)
  {
    if (years_string == "run2") years = {2016, 2017, 2018};
    else if (years_string!="")
    {
      size_t found = years_string.find(",");
      size_t start = 0;
      string _tmp = "";
      while (found != string::npos) {
        _tmp = years_string.substr(start, found-start);
        years.insert(atoi(_tmp.c_str()));
        start = found+1;
        found = years_string.find(",",start);
      } 
      _tmp = years_string.substr(start, found-start);
      years.insert(atoi(_tmp.c_str()));
    }
  }

  void findMassPoints(std::string signal_folder, std::vector<std::pair<std::string, std::string> > & mass_points) {
    set<string> signal_files = Glob(signal_folder+"/*.root");
    string mLSP, mChi;
    for (string signal_file : signal_files) {
      filenameToMassPoint(signal_file, mChi, mLSP);
      mass_points.push_back({mChi, mLSP});
    }
  }

  void filenameToMassPoint(std::string filename, string & mChi, string & mLSP) {
    size_t start_mChi = filename.find("mChi")+5;
    size_t end_mChi = filename.find("_", start_mChi);
    string mChi_string = filename.substr(start_mChi, (end_mChi-start_mChi));
    size_t start_mLSP = filename.find("mLSP")+5;
    size_t end_mLSP = filename.find("_", start_mLSP);
    if (end_mLSP == std::string::npos) end_mLSP = filename.find(".", start_mLSP);
    string mLSP_string = filename.substr(start_mLSP, (end_mLSP-start_mLSP));
    mChi = mChi_string;
    mLSP = mLSP_string;
  }

  void filenameToMassPoint(std::string filename, string & mGluino, string & mChi, string & mLSP) {
    size_t start_mGluino = filename.find("mGluino")+8;
    size_t end_mGluino = filename.find("_", start_mGluino);
    string mGluino_string = filename.substr(start_mGluino, (end_mGluino-start_mGluino));
    size_t start_mChi = filename.find("mChi")+5;
    size_t end_mChi = filename.find("_", start_mChi);
    string mChi_string = filename.substr(start_mChi, (end_mChi-start_mChi));
    size_t start_mLSP = filename.find("mLSP")+5;
    size_t end_mLSP = filename.find("_", start_mLSP);
    if (end_mLSP == std::string::npos) end_mLSP = filename.find(".", start_mLSP);
    string mLSP_string = filename.substr(start_mLSP, (end_mLSP-start_mLSP));
    mGluino = mGluino_string;
    mChi = mChi_string;
    mLSP = mLSP_string;
  }

  std::string setProcessName(std::string const & model, int const & mGluino, int const & mLSP)
  {
    return model+"_"+to_string(mGluino)+"_"+to_string(mLSP);
  }

  std::string setProcessNameLong(std::string const & model, int const & mGluino, int const & mLSP)
  {
    // Tune is added for utilites::parseMasses
    if (model == "T5HH") return model+"_mGluino-"+to_string(mGluino)+"_mLSP-"+to_string(mLSP)+"_Tune";
    return model+"_mChi-"+to_string(mGluino)+"_mLSP-"+to_string(mLSP)+"_Tune";
  }

  void getInfoFromProcessName(std::string const & processName, std::string & model, int & mGluino, int & mLSP)
  {
    vector<string> names;
    stringToVectorString(processName, names, "_");
    model = names[0];
    mGluino = atoi(names[1].c_str());
    mLSP = atoi(names[2].c_str());
  }

  void setABCDBins(map<string, string> xBins, map<string, string> yBins, map<string, vector<pair<string, string> > > dimensionBins, vector<pair<string, string> > & sampleBins)
  {
    // Combine dimensions to one list of dimensions.
    combineDimensionBins(dimensionBins);
  
    auto dimension = dimensionBins.begin();
    for (unsigned iBin = 0; iBin < dimension->second.size(); ++iBin)
    {
      for (auto yBin: yBins)
      {
        for (auto xBin : xBins)
        {
          string label = "x"+xBin.first+"_y"+yBin.first+"_"+(dimension->second)[iBin].first;
          string cut = xBin.second+"&&"+yBin.second+"&&"+(dimension->second)[iBin].second;
          //cout<<"bins: "<<label<<" "<<cut<<endl;
          sampleBins.push_back({label, cut});
        }
      }
    }
  }

  // 1. RSR, RCR, BSR, BCR
  // 2. BSR, BCR, RSR, RCR
  // 3. RSR, BSR, RCR, BCR
  // 4. BSR, RSR, BCR, RCR
  void setABCDBinsPriority(map<string, string> xBins, map<string, string> yBins, map<string, vector<pair<string, string> > > dimensionBins, vector<pair<string, string> > & sampleBins, int priority)
  {
    // Combine dimensions to one list of dimensions.
    combineDimensionBins(dimensionBins);
  
    auto dimension = dimensionBins.begin();
    for (unsigned iBin = 0; iBin < dimension->second.size(); ++iBin)
    {
      for (auto yBin: yBins)
      {
        for (auto xBin : xBins)
        {
          string label = "x"+xBin.first+"_y"+yBin.first+"_"+(dimension->second)[iBin].first;
          // Check if label is signal region or sideband region
          bool isSignalRegion = false;
          if (label.find("xsig")!=string::npos && label.find("ysig")!=string::npos) isSignalRegion = true;
          // Apply cut according to SR,CR and priority
          string cut = xBin.second+"&&"+yBin.second+"&&"+(dimension->second)[iBin].second;
          if (priority==1) cut = cut;
          else if (priority==2) cut = cut+"&&!boostSignalRegion&&!boostControlRegion";
          else if (priority==3) {
            if (isSignalRegion) cut = cut;
            else cut = cut+"&&!boostSignalRegion";
          } else if (priority==4) {
            if (isSignalRegion) cut = cut+"&&!boostSignalRegion";
            else cut = cut+"&&!boostSignalRegion&&!boostControlRegion";
          }
          //cout<<"bins: "<<label<<" "<<cut<<endl;
          sampleBins.push_back({label, cut});
        }
      }
    }
  }

  void combineDimensionBins(map<string, vector<pair<string, string> > > & dimensionBins)
  {
    unsigned nDimensions = dimensionBins.size();
    if (nDimensions == 1) return;
    auto dimension1 = dimensionBins.begin();
    auto dimension2 = dimension1++;
    // Combine dimensions
    dimensionBins[dimension1->first+"_"+dimension2->first];
    for (unsigned iDimension1 = 0; iDimension1 < dimension1->second.size(); ++iDimension1)
    {
      for (unsigned iDimension2 = 0; iDimension2 < dimension2->second.size(); ++iDimension2)
      {
        dimensionBins[dimension1->first+"_"+dimension2->first].push_back({dimension1->second[iDimension1].first+"_"+dimension2->second[iDimension2].first,dimension1->second[iDimension1].second+"&&"+dimension2->second[iDimension2].second});
        //cout<<dimension1->first+"_"+dimension2->first<<":"<<dimension1->second[iDimension1].first+"_"+dimension2->second[iDimension2].first<<" "<<dimension1->second[iDimension1].second+"&&"+dimension2->second[iDimension2].second<<endl;
      }
    }
    dimensionBins.erase(dimension1);
    dimensionBins.erase(dimension2);
    return combineDimensionBins(dimensionBins);
  }

  void setMcProcesses(set<int> years, map<string, string> samplePaths, NamedFunc filters, map<string, vector<shared_ptr<Process> > > & sampleProcesses)
  {
    // mcTag['ttx'] = {*.root,*.root}
    map<string, set<string> > mcTags;
  
    // mcTags["t#bar{t}+X"]     = set<string>({"*TTJets_*Lept*"});
    mcTags["t#bar{t}+X"]     = set<string>({"*TTJets_*Lept*", "*_TTZ*.root", "*_TTW*.root",
                                       "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root"});
    mcTags["V+jets"]   = set<string>({"*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root"});
    mcTags["QCD"]     = set<string>({"*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                                   "*_QCD_HT700to1000_*", "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                                   "*_QCD_HT2000toInf_*"});
    mcTags["Other"]   = set<string>({"*_WH*.root", "*_ZH_HToBB*.root", "*_ST_*.root",
                                       "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});
  
    map<string, int> mcColors;
    Palette colors("txt/colors.txt", "default");
    mcColors["t#bar{t}+X"] = colors("tt_1l");
    mcColors["V+jets"] = kOrange+1;
    mcColors["QCD"] = colors("other");
    mcColors["Other"] = kGreen+1;
  
    // mcPaths['ttx'] = {folder/*.root,folder/*.root}
    map<string, set<string> > mcPaths;
    for (auto mcTag : mcTags)
    {
      for (auto file : mcTag.second)
      {
        for (auto year : years)
        {
          mcPaths[mcTag.first].insert(samplePaths["mc_"+to_string(year)]+file);
          cout<<"Adding "<<samplePaths["mc_"+to_string(year)]+file<<endl;
        }
      }
    }
  
    sampleProcesses["mc"] = vector<shared_ptr<Process> > ();
    for (auto mcPath : mcPaths)
    {
      sampleProcesses["mc"].push_back(Process::MakeShared<Baby_pico>(mcPath.first, Process::Type::background, mcColors[mcPath.first], mcPath.second, filters));
    }
  }
  
  void setDataProcesses(set<int> years, map<string, string> samplePaths, NamedFunc filters, map<string, vector<shared_ptr<Process> > > & sampleProcesses)
  {
    // dataPaths['data'] = {folder/*.root,folder/*.root}
    map<string, set<string> > dataPaths;
    for (auto year : years)
    {
      dataPaths["Data"].insert(samplePaths["data_"+to_string(year)]+"*.root");
      cout<<"Adding "<<samplePaths["data_"+to_string(year)]+"*.root"<<endl;
    }
  
    sampleProcesses["data"] = vector<shared_ptr<Process> > ();
    for (auto dataPath : dataPaths)
    {
      sampleProcesses["data"].push_back(Process::MakeShared<Baby_pico>(dataPath.first, Process::Type::data, 1, dataPath.second, filters));
    }
  }
  
  void setSignalProcesses(vector<pair<string, string> > massPoints, set<int> years, map<string, string> samplePaths, NamedFunc filters, map<string, vector<shared_ptr<Process> > > & sampleProcesses)
  {
    // mcPaths['ttx'] = {folder/*.root,folder/*.root}
    map<string, set<string> > signalPaths;
    for (auto massPoint : massPoints)
    {
     for (auto year : years)
     {
       signalPaths[HigUtilities::setProcessName("TChiHH",atoi(massPoint.first.c_str()),atoi(massPoint.second.c_str()))].insert(samplePaths["signal_"+to_string(year)]+"*mChi-"+massPoint.first+"_mLSP-"+massPoint.second+"_*.root");
       cout<<"Adding "<<samplePaths["signal_"+to_string(year)]+"*mChi-"+massPoint.first+"_mLSP-"+massPoint.second+"_*.root"<<endl;
     }
    }
  
    sampleProcesses["signal"] = vector<shared_ptr<Process> > ();
    for (auto signalPath : signalPaths)
    {
      sampleProcesses["signal"].push_back(Process::MakeShared<Baby_pico>(signalPath.first, Process::Type::background, kBlack, signalPath.second, filters));
    }
  }

  void setSignalProcessesT5HH(vector<pair<string, string> > massPoints, set<int> years, map<string, string> samplePaths, NamedFunc filters, map<string, vector<shared_ptr<Process> > > & sampleProcesses)
  {
    map<string, set<string> > signalPaths;
    for (auto massPoint : massPoints)
    {
     for (auto year : years)
     {
       signalPaths[HigUtilities::setProcessName("T5HH",atoi(massPoint.first.c_str()),atoi(massPoint.second.c_str()))].insert(samplePaths["signal_"+to_string(year)]+"*mGluino-"+massPoint.first+"_mChi-"+to_string(stoi(massPoint.first)-50)+"_mLSP-"+massPoint.second+"_*.root");
       cout << "Adding " << samplePaths["signal_"+to_string(year)]+"*mGluino-"+massPoint.first+"_mChi-"+to_string(stoi(massPoint.first)-50)+"_mLSP-"+massPoint.second+"_*.root" << endl;
     }
    }
  
    sampleProcesses["signal"] = vector<shared_ptr<Process> > ();
    for (auto signalPath : signalPaths)
    {
      sampleProcesses["signal"].push_back(Process::MakeShared<Baby_pico>(signalPath.first, Process::Type::background, kBlack, signalPath.second, filters));
    }
  }

  void addBinCuts(vector<pair<string, string> > sampleBins, string baseline, NamedFunc weight, string tag, RowInformation & cutRows)
  {
    for (auto binCut : sampleBins)
    {
      cutRows.labels.push_back(tag +"_" + binCut.first);
      cutRows.tableRows.push_back(TableRow("", baseline&&binCut.second,0,0,weight));
      //cout<<cutRows.labels.back()<<" "<<cutRows.tableRows.back().cut_.Name()<<endl;
    }
  }
  
  void addBinCuts(vector<pair<string, string> > sampleBins, string baseline, NamedFunc weight, string tag, TString (*replaceFunc)(TString), RowInformation & cutRows)
  {
    for (auto binCut : sampleBins)
    {
      cutRows.labels.push_back(tag +"_" + binCut.first);
      cutRows.tableRows.push_back(TableRow("", replaceFunc(baseline+"&&"+binCut.second),0,0,weight));
      //cout<<cutRows.labels.back()<<" "<<cutRows.tableRows.back().cut_.Name()<<endl;
    }
  }
  
  void addBinCuts(vector<pair<string, string> > sampleBins, string baseline, NamedFunc weight, string tag, TString (*replaceFunc)(TString, size_t) ,RowInformation & cutRows)
  {
    for (auto binCut : sampleBins)
    {
      cutRows.labels.push_back(tag +"_up_" + binCut.first);
      cutRows.tableRows.push_back(TableRow("", replaceFunc(baseline+"&&"+binCut.second,1),0,0,weight));
      cutRows.labels.push_back(tag +"_down_" + binCut.first);
      cutRows.tableRows.push_back(TableRow("", replaceFunc(baseline+"&&"+binCut.second,2),0,0,weight));
      //cout<<cutRows.labels.back()<<" "<<cutRows.tableRows.back().cut_.Name()<<endl;
    }
  }
  
  void makePlots(map<string, RowInformation > & cutTable, map<string, vector<shared_ptr<Process> > > & sampleProcesses, float luminosity, PlotMaker & pm, bool verbose)
  {
    int tableIndex = 0;
    for (auto cutRow : cutTable)
    {
      // type = mc, data, signal
      string const & type = cutRow.first;
      pm.Push<Table>(type, cutTable[type].tableRows, sampleProcesses[type], true, false);
      cutTable[type].tableIndex= tableIndex;
      tableIndex++;
      if (verbose) for (auto tableRow : cutTable[type].tableRows) cout<<tableRow.cut_.Name()<<endl;
    }
    pm.multithreaded_ = true;
    pm.min_print_ = true; 
    pm.MakePlots(luminosity);  
  }

  void makePlots(map<string, RowInformation > & cutTable, map<string, HistInformation> & histInfo, map<string, vector<shared_ptr<Process> > > & sampleProcesses, float luminosity, PlotMaker & pm, bool verbose)
  {
    int tableIndex = 0;
    for (auto cutRow : cutTable)
    {
      // type = mc, data, signal
      string const & type = cutRow.first;
      pm.Push<Table>(type, cutTable[type].tableRows, sampleProcesses[type], true, false);
      cutTable[type].tableIndex= tableIndex;
      tableIndex++;
      if (verbose) for (auto tableRow : cutTable[type].tableRows) cout<<tableRow.cut_.Name()<<endl;
    }
    for (auto single_hist_info : histInfo)
    {
      // type = mc, data, signal
      string const & type = single_hist_info.first;
      std::vector<PlotOpt> plot_opts = {*(histInfo[type].plot_opt_)};
      pm.Push<Hist1D>(*(histInfo[type].axis_), *(histInfo[type].cut_), 
          sampleProcesses[type], plot_opts).Weight(*(histInfo[type].weight_)).DrawPlot(false);
      histInfo[type].figure_index = tableIndex;
      tableIndex++;
      //if (verbose) for (auto tableRow : cutTable[type].tableRows) cout<<tableRow.cut_.Name()<<endl;
    }
    pm.multithreaded_ = true;
    pm.min_print_ = true; 
    pm.MakePlots(luminosity);  
  }

  void addToMapYields(string const & label, GammaParams & yield, TableRow & yieldMeta, map<string, pair<GammaParams, TableRow> > & mYields)
  {
      if (mYields.find(label) != mYields.end())
      {
        cout<<"[Warning] addYieldData(): "<<label<<" already exists in mYields. Will not add data."<<endl;
      }
      pair<GammaParams, TableRow> yieldData = {yield, yieldMeta};
      mYields.insert(pair<string, pair<GammaParams, TableRow> > (label, yieldData));
  }
  
  void fillDataYields(PlotMaker & pm, RowInformation & dataRow, map<string, pair<GammaParams, TableRow> > & mYields, bool verbose)
  {
    Table * yieldTable;
    yieldTable = static_cast<Table*>(pm.Figures()[dataRow.tableIndex].get());
    vector<GammaParams> yields = yieldTable->DataYield();
    for (unsigned ipar(0); ipar< yields.size(); ipar++) 
    {
      string const & label = dataRow.labels[ipar];
      addToMapYields(label, yields[ipar], dataRow.tableRows[ipar], mYields);
      if (verbose) cout<<"[data] "<<label<<": "<<mYields.at(label).first.Yield()<<endl;
      //cout<<label<<": "<<mYields.at(label).first.Yield()<<" "<<mYields.at(label).second.cut_.Name()<<" "<<mYields.at(label).second.weight_.Name()<<endl;
    }
  }
  
  void fillMcYields(PlotMaker & pm, float luminosity, RowInformation & mcRow, map<string, pair<GammaParams, TableRow> > & mYields, bool verbose)
  {
    Table * yieldTable;
    yieldTable = static_cast<Table*>(pm.Figures()[mcRow.tableIndex].get());
    vector<GammaParams> yields = yieldTable->BackgroundYield(luminosity);
    for (unsigned ipar(0); ipar< yields.size(); ipar++) 
    {
      string const & label = mcRow.labels[ipar];
      addToMapYields(label, yields[ipar], mcRow.tableRows[ipar], mYields);
      if (verbose) cout<<"[mc] "<<label<<": "<<mYields.at(label).first.Yield()<<endl;
      //cout<<label<<": "<<mYields.at(label).first.Yield()<<" "<<mYields.at(label).second.cut_.Name()<<" "<<mYields.at(label).second.weight_.Name()<<endl;
    }
  }
  
  void fillSignalYieldsProcesses(PlotMaker & pm, float luminosity, vector<shared_ptr<Process> > & signalProcesses, RowInformation & signalRow, map<string, pair<GammaParams, TableRow> > & mYields, bool verbose)
  {
    Table * yieldTable;
    yieldTable = static_cast<Table*>(pm.Figures()[signalRow.tableIndex].get());
    for (auto & process : signalProcesses)
    {
      vector<GammaParams> yields = yieldTable->Yield(process.get(), luminosity);
      for (unsigned ipar(0); ipar< yields.size(); ipar++) 
      {
        string const & label = process->name_ + "_" +signalRow.labels[ipar];
        addToMapYields(label, yields[ipar], signalRow.tableRows[ipar], mYields);
        if (verbose) {
          cout<<"[sig] "<<label<<": "<<mYields.at(label).first.Yield()<<endl;
        }
        //cout<<label<<": "<<mYields.at(label).first.Yield()<<" "<<mYields.at(label).second.cut_.Name()<<" "<<mYields.at(label).second.weight_.Name()<<endl;
      }
    }
  }
  
  void fillAverageGenMetYields(vector<shared_ptr<Process> > & processes, vector<pair<string, string> > sampleBins, string const & signalTag, string const & signalGenMetTag, string const & signalAverageGenMetTag, map<string, pair<GammaParams, TableRow> > & mYields, bool verbose)
  {
    for (auto & process : processes)
    {
      string const & processName = process->name_;
      for(auto & bin : sampleBins)
      {
        string yieldLabel = processName + "_" + signalTag + "_" + bin.first;
        string genMetYieldLabel = processName + "_" + signalGenMetTag + "_" + bin.first;
        //cout<<yieldLabel<<" "<<genMetYieldLabel<<endl;
        float averageValue = 0.5*(mYields.at(yieldLabel).first.Yield() + mYields.at(genMetYieldLabel).first.Yield());
        //float averageError = hypot(0.5*mYields.at(yieldLabel).first.Uncertainty(),0.5*mYields.at(genMetYieldLabel).first.Uncertainty());
        float averageError = max(mYields.at(yieldLabel).first.Uncertainty(),mYields.at(genMetYieldLabel).first.Uncertainty());
        GammaParams average;
        average.SetYieldAndUncertainty(averageValue,averageError);
        TableRow blank("0.5("+yieldLabel+"+"+genMetYieldLabel+")");
        addToMapYields(processName+"_"+signalAverageGenMetTag+"_"+bin.first, average, blank, mYields);
        if (verbose) cout<<processName+"_"+signalAverageGenMetTag+"_"+bin.first<<" "<<averageValue<<" "<<averageError<<endl;
      }
    }
  }

  string to_pico_names(string in_cuts)
  {
    map<string, string> convert_map;
    convert_map["nbdt"] = "nbt";
    convert_map["nbdm"] = "nbm";
    convert_map["nbdl"] = "nbl";
    convert_map["higd_dm"] = "hig_cand_dm[0]";
    convert_map["higd_drmax"] = "hig_cand_drmax[0]";
    convert_map["higd_am"] = "hig_cand_am[0]";
    convert_map["nvleps"] = "nvlep";
    convert_map["njets"] = "njet";
    convert_map["nels"] = "nel";
    convert_map["nmus"] = "nmu";
    convert_map["mgluino"] = "mprod";
    convert_map["w_btag_deep"] = "w_btag";
    convert_map["w_bhig_deep"] = "w_bhig";

    for (auto convert : convert_map)
    {
      in_cuts = std::regex_replace(in_cuts, std::regex(convert.first), convert.second);
    }

    return in_cuts;
  }

  float signal_lepton_pt(std::vector<float>* const lep_pt, std::vector<bool>* const lep_sig) {
    float r_signal_lepton_pt = 0;
    for (unsigned int lep_idx = 0; lep_idx < lep_sig->size(); lep_idx++) {
      if (lep_sig->at(lep_idx)) {
        r_signal_lepton_pt = lep_pt->at(lep_idx);
        break;
      }
    }
    return r_signal_lepton_pt;
  }

  string getLuminosityString(string const & year_string, int const & decimals) {
    set<int> years;
    HigUtilities::parseYears(year_string, years);
    float total_luminosity = 0;
    for (auto const & year : years) {
      if (year == 2016) total_luminosity += 36.22;
      if (year == 2017) total_luminosity += 41.49;
      if (year == 2018) total_luminosity += 59.65;
    }
    string total_luminosity_string = RoundNumber(total_luminosity, decimals, 1).Data();
    return total_luminosity_string;
  }

}
