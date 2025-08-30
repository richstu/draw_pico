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

#endif
