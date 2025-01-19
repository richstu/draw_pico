#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <bitset>
#include <unistd.h>
#include <getopt.h>
#include "TError.h"
#include "TColor.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "core/baby.hpp"
#include "core/named_func.hpp"
#include "core/table_row.hpp"
#include "core/utilities.hpp"
#include "core/plot_opt.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"

#include "core/process.hpp"
#include "core/plot_maker.hpp"
#include "core/palette.hpp"

#include "zgamma/zg_utilities.hpp"
#include "zgamma/zg_functions.hpp"
#include "zgamma/categorization_utilities.hpp"

using namespace ZgUtilities;
using namespace CatUtilities;
using namespace ZgFunctions;

namespace CRUtilities{

bool select_era(int argc, char *argv[]);

extern const NamedFunc hasNanoPhoton;
extern const NamedFunc hasNll;
extern const NamedFunc Nllgamma_g0_minsel;

extern const NamedFunc higgs_window;
extern const NamedFunc cr_sideband;
extern const NamedFunc cr_fake_photon;
extern const NamedFunc cr_Zfsr_photon;
extern const NamedFunc cr_soft_photon;

extern const std::vector<NamedFunc>   controlregions_vec;  
extern const std::vector<std::string> controlregions_str_vec;
extern const std::vector<NamedFunc>   controlregions_refit_vec;  

//This code defines additional control regions split by the selections listed below
extern const std::vector<NamedFunc> map_pair_selections;
extern const std::unordered_map<int,int>    category_specific_control_regions;
extern const std::unordered_map<int,int>    category_specific_control_regions_whbs;
extern const std::unordered_map<int,std::string> category_specific_control_regions_names; 

std::vector<PlotOpt> an_plotting_options(std::string select_options="");


void constructCutflowTable(std::vector<TableRow> & tablerows, NamedFunc weight, int electron_or_muon, bool isMC = false, NamedFunc current_cut="1");

//std::vector<shared_ptr<Process>> control_region_procs(std::string year_select, int category=-1);
//extern const std::tuple<std::vector<NamedFunc>,std::vector<std::string>> create_control_region_selections_and_labels(std::string label, bool plot_all, int category=-1);

}
