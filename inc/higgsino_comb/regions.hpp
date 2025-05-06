#ifndef H_REGIONS
#define H_REGIONS

#include <cstddef>

#include <string>

#include "core/named_func.hpp"
#include "higgsino_comb/vardef.hpp"

using namespace std;

namespace regions_4b{
 
 extern const NamedFunc sig_decay_4b;
 extern const NamedFunc dphi_res;

 extern vector<pair<string, NamedFunc>> cuts_4b_res;
 extern vector<pair<string, NamedFunc>> cuts_4b_CR1;
 extern vector<pair<string, NamedFunc>> cuts_4b_CR2;
 extern vector<pair<string, NamedFunc>> cuts_4b_CR3;
 extern vector<pair<string, NamedFunc>> cuts_4b_CR1b;

 extern const NamedFunc res_baseline;
 extern const NamedFunc res_nm_higam;
 extern const NamedFunc res_nm_met_higam;
 extern const NamedFunc res_test;

 extern const NamedFunc boo_skimcuts;
 extern const NamedFunc boo_baseline;

 extern const NamedFunc CR1_el;
 extern const NamedFunc CR1_mu;
 extern const NamedFunc CR2_el;
 extern const NamedFunc CR2_mu;
 extern const NamedFunc CR3;
 extern const NamedFunc CR1b_el;
 extern const NamedFunc CR1b_mu;

}

#endif

