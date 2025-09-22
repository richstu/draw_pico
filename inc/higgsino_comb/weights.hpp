#ifndef H_WEIGHTS
#define H_WEIGHTS

#include <cstddef>

#include <string>

#include "core/named_func.hpp"
#include "higgsino_comb/weights.hpp"

using namespace std;

namespace weights{
  extern const NamedFunc w_run2_datavmc;
  extern const NamedFunc w_lumi_hb;
  extern const NamedFunc w_lumi_hb_mettrig;
  extern const NamedFunc w_lumi_gmsb;
  extern const NamedFunc w_lumi_gmsb_mettrig;
}

namespace triggers{
  extern const NamedFunc single_ele_trig;
  extern const NamedFunc single_muon_trig;
  extern const NamedFunc ptmiss_htmiss_trig;
  extern const NamedFunc diphoton_trig;
}

namespace filters{
  extern const NamedFunc pass_hemveto;
  extern const NamedFunc pass_filters;
}

#endif
