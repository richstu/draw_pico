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
 extern const NamedFunc res_nm_met;
 extern const NamedFunc res_nm_njet;
 extern const NamedFunc res_nm_nb;
 extern const NamedFunc res_nm_fakemet;
 extern const NamedFunc res_nm_higam;
 extern const NamedFunc res_nm_higdm;
 extern const NamedFunc res_nm_higdrmax;
 extern const NamedFunc res_nm_met_higam;
 extern const NamedFunc res_test;

 extern const NamedFunc boo_skimcuts;
 extern const NamedFunc boo_baseline;
 extern const NamedFunc boo_nm_met;
 extern const NamedFunc boo_nm_ht;
 extern const NamedFunc boo_nm_ak8pt;
 extern const NamedFunc boo_nm_msd;

 extern const NamedFunc CR1_el;
 extern const NamedFunc CR1_mu;
 extern const NamedFunc CR2_el;
 extern const NamedFunc CR2_mu;
 extern const NamedFunc CR3;
 extern const NamedFunc CR1b_el;
 extern const NamedFunc CR1b_mu;

}

#endif

