#ifndef H_PATHS
#define H_PATHS

#include <cstddef>

#include <string>

#include "core/named_func.hpp"
#include <map>
#include "core/process.hpp"

using namespace std;

namespace bbgg_paths{
 
 extern string bfoldersigHH;
 extern string bfoldersigZH;
 extern string bfolderbkg;
 extern string bfolderbkg_smh;
 extern string bfolderbkg_gjet2;
 extern string bfolderdata;

 struct ProcInfo {
      string basepath;
      string proclegend;
      string proccolor;
      string special_proccuts; 
      Process::Type type;
    };

 extern map<string, ProcInfo> SampleInfo;
 

}

namespace paths_4b{
  extern string sig_path;
  extern string bkg_met150;
  extern string bkg_1lep;
  extern string bkg_2lep;
  extern string data_met150;
  extern string data_1lep;
  extern string data_2lep;

  struct ProcInfo {
       string basepath;
       string proclegend;
       string proccolor;
       string special_proccuts; 
       Process::Type type;
     };
  extern map<string, ProcInfo> signal_info;
  extern map<vector<string>, ProcInfo> sample_info_met150;
  extern map<vector<string>, ProcInfo> sample_info_1lep;
  extern map<vector<string>, ProcInfo> sample_info_2lep;
  extern map<vector<string>, ProcInfo> sample_info_fitting;
  extern map<string, ProcInfo> signal_info_all;

}

#endif
