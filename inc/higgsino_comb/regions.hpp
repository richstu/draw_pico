#ifndef H_REGIONS
#define H_REGIONS

#include <cstddef>

#include <string>

#include "core/named_func.hpp"

using namespace std;

namespace regions_4b{

 extern const vector<pair<string, NamedFunc>> cuts_4b_res;
 extern const vector<pair<string, NamedFunc>> cuts_4b_CR1;
 extern const vector<pair<string, NamedFunc>> cuts_4b_CR2;

 extern const string res_baseline;
 extern const string res_nm_higam;
 extern const string res_nm_met_higam;

 extern const string CR1_el;
 extern const string CR1_mu;
 extern const string CR2_el;
 extern const string CR2_mu;
 extern const string CR3;

}

#endif

