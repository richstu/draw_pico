#include <algorithm>
#include <stdlib.h>
#include <regex>
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/functions.hpp"
#include "core/hist1d.hpp"
#include "core/utilities.hpp"
#include "core/palette.hpp"
#include "core/plot_opt.hpp"
#include "TVector2.h"
#include "TMath.h"
#include "Math/Vector4D.h"

#include "core/utilities.hpp"
#include "core/config_parser.hpp"
#include "higgsino_comb/vardef.hpp"

using namespace std;

namespace vardef{

 const NamedFunc xsec_ratio("xsec_ratio", [](const Baby &b) -> NamedFunc::ScalarType{ // ratio of (2d xsec)/(GMSB xsec)
   float hig_mass = 0;
   float r;
   for (size_t i = 0; i < b.mc_pt()->size(); i++){
     if (b.mc_id()->at(i) == 1000023) { hig_mass = b.mc_mass()->at(i);}
   }
   if      (hig_mass >= 125 && hig_mass <= 129)   { r = 1.44725/7.6022; }
   else if (hig_mass >= 148 && hig_mass <= 152)   { r = 0.71514/3.83231; }
   else if (hig_mass >= 173 && hig_mass <= 177)   { r = 0.419059/2.26794; }
   else if (hig_mass >= 198 && hig_mass <= 202)   { r = 0.244213/1.33562; }
   else if (hig_mass >= 223 && hig_mass <= 227)   { r = 0.156286/0.860597; }
   else if (hig_mass >= 248 && hig_mass <= 252)   { r = 0.104252/0.577314; }
   else if (hig_mass >= 273 && hig_mass <= 277)   { r = 0.0719125/0.400107; }
   else if (hig_mass >= 298 && hig_mass <= 302)   { r = 0.0509994/0.284855; }
   else if (hig_mass >= 323 && hig_mass <= 327)   { r = 0.0369715/0.20736; }
   else if (hig_mass >= 348 && hig_mass <= 352)   { r = 0.273286/0.153841; }
   else if (hig_mass >= 373 && hig_mass <= 377)   { r = 0.0205429/0.116006; }
   else if (hig_mass >= 398 && hig_mass <= 402)   { r = 0.0156691/0.0887325; }
   else if (hig_mass >= 423 && hig_mass <= 427)   { r = 0.0120965/0.0686963; }
   else if (hig_mass >= 448 && hig_mass <= 452)   { r = 0.00944017/0.0537702; }
   else if (hig_mass >= 473 && hig_mass <= 477)   { r = 0.00743587/0.0424699; }
   else if (hig_mass >= 498 && hig_mass <= 502)   { r = 0.00590757/0.0338387; }
   else if (hig_mass >= 523 && hig_mass <= 527)   { r = 0.00473235/0.0271867; }
   else if (hig_mass >= 548 && hig_mass <= 552)   { r = 0.0038167/0.0219868; }
   else if (hig_mass >= 573 && hig_mass <= 577)   { r = 0.00309847/0.0179062; }
   else if (hig_mass >= 598 && hig_mass <= 602)   { r = 0.00253015/0.0146677; }
   else if (hig_mass >= 623 && hig_mass <= 627)   { r = 0.00207755/0.012062; }
   else if (hig_mass >= 648 && hig_mass <= 652)   { r = 0.00171418/0.00997406; }
   else if (hig_mass >= 673 && hig_mass <= 677)   { r = 0.0014199/0.00828246; }
   else if (hig_mass >= 698 && hig_mass <= 702)   { r = 0.00118113/0.00689981; }
   else if (hig_mass >= 723 && hig_mass <= 727)   { r = 0.00098639/0.0057833; }
   else if (hig_mass >= 748 && hig_mass <= 752)   { r = 0.000826366/0.0048731; }
   else if (hig_mass >= 773 && hig_mass <= 777)   { r = 0.000694985/0.00409781; }
   else if (hig_mass >= 798 && hig_mass <= 802)   { r = 0.000586211/0.00346143; }
   else if (hig_mass >= 823 && hig_mass <= 827)   { r = 0.000495914/0.0029337; }
   else if (hig_mass >= 848 && hig_mass <= 852)   { r = 0.000420556/0.0024923; }
   else if (hig_mass >= 873 && hig_mass <= 877)   { r = 0.000361029/0.00213679; }
   else if (hig_mass >= 898 && hig_mass <= 902)   { r = 0.000305935/0.00180616; }
   else if (hig_mass >= 923 && hig_mass <= 927)   { r = 0.000262621/0.00155453; }
   else if (hig_mass >= 948 && hig_mass <= 952)   { r = 0.00022285/0.00132692; }
   else if (hig_mass >= 973 && hig_mass <= 977)   { r = 0.0001909/0.00112975; }
   else if (hig_mass >= 998 && hig_mass <= 1002)  { r = 0.00016428/0.000968853; }
   else if (hig_mass >= 1023 && hig_mass <= 1027) { r = 0.00014139/0.000840602; }
   else if (hig_mass >= 1048 && hig_mass <= 1052) { r = 0.000121865/0.00731306; }
   else if (hig_mass >= 1073 && hig_mass <= 1077) { r = 0.000105913/0.000627083; }
   else if (hig_mass >= 1098 && hig_mass <= 1102) { r = 0.0000912469/0.000538005; }
   else if (hig_mass >= 1123 && hig_mass <= 1127) { r = 0.0000793058/0.00046747; }
   else if (hig_mass >= 1148 && hig_mass <= 1152) { r = 0.0000684561/0.000405108; }
   else if (hig_mass >= 1173 && hig_mass <= 1177) { r = 0.0000593602/0.000348261; }
   else if (hig_mass >= 1198 && hig_mass <= 1202) { r = 0.0000516263/0.000299347; }
   else if (hig_mass >= 1223 && hig_mass <= 1227) { r = 0.000044906/0.000265935; }
   else if (hig_mass >= 1248 && hig_mass <= 1252) { r = 0.0000391587/0.00024071; }
   else if (hig_mass >= 1273 && hig_mass <= 1277) { r = 0.0000343135/0.000190411; }
   else if (hig_mass >= 1298 && hig_mass <= 1302) { r = 0.0000299353/0.000160765; }
   else if (hig_mass >= 1323 && hig_mass <= 1327) { r = 0.0000262223/0.000136272; }
   else if (hig_mass >= 1348 && hig_mass <= 1352) { r = 0.0000228072/0.000111174; }
   else if (hig_mass >= 1373 && hig_mass <= 1377) { r = 0.0000200393/0.0000974728; }
   else if (hig_mass >= 1398 && hig_mass <= 1402) { r = 0.0000175031/0.0000780263; }
   else if (hig_mass >= 1423 && hig_mass <= 1427) { r = 0.0000153144/0.0000696843; }
   else if (hig_mass >= 1448 && hig_mass <= 1452) { r = 0.0000134572/0.0000696962; }
   else if (hig_mass >= 1473 && hig_mass <= 1477) { r = 0.0000117047/0.0000498006; }
   else {r = 0;}
   return r;
 });

 const NamedFunc jets_hh_pt("jets_hh_pt", [](const Baby &b) -> NamedFunc::VectorType{
   vector<double> pt;
   pt.push_back(b.jet_pt()->at(b.jet_h1_indices()->at(0)));
   pt.push_back(b.jet_pt()->at(b.jet_h1_indices()->at(1)));
   pt.push_back(b.jet_pt()->at(b.jet_h2_indices()->at(0)));
   pt.push_back(b.jet_pt()->at(b.jet_h2_indices()->at(1)));
   return pt;
 });

 const NamedFunc jets_hh_eta("jets_hh_eta", [](const Baby &b) -> NamedFunc::VectorType{
   vector<double> eta;
   eta.push_back(b.jet_eta()->at(b.jet_h1_indices()->at(0)));
   eta.push_back(b.jet_eta()->at(b.jet_h1_indices()->at(1)));
   eta.push_back(b.jet_eta()->at(b.jet_h2_indices()->at(0)));
   eta.push_back(b.jet_eta()->at(b.jet_h2_indices()->at(1)));
   return eta;
 });

 const NamedFunc jets_hh_phi("jets_hh_phi", [](const Baby &b) -> NamedFunc::VectorType{
   vector<double> phi;
   phi.push_back(b.jet_phi()->at(b.jet_h1_indices()->at(0)));
   phi.push_back(b.jet_phi()->at(b.jet_h1_indices()->at(1)));
   phi.push_back(b.jet_phi()->at(b.jet_h2_indices()->at(0)));
   phi.push_back(b.jet_phi()->at(b.jet_h2_indices()->at(1)));
   return phi;
 });

 const NamedFunc deepflav_hh_bscore("deepflav_bscore", [](const Baby &b) -> NamedFunc::VectorType{
   vector<double> deepflav;
   deepflav.push_back(b.jet_deepflav()->at(b.jet_h1_indices()->at(0)));
   deepflav.push_back(b.jet_deepflav()->at(b.jet_h1_indices()->at(1)));
   deepflav.push_back(b.jet_deepflav()->at(b.jet_h2_indices()->at(0)));
   deepflav.push_back(b.jet_deepflav()->at(b.jet_h2_indices()->at(1)));
   return deepflav;
 });

 const NamedFunc num_b("num_b", [](const Baby &b) -> NamedFunc::ScalarType{
   int nb;
   if (b.nbdft()==2 && b.nbdfm()==2) { nb = 2; }
   else if (b.nbdft()>=2 && b.nbdfm()==3 && b.nbdfl()==3) { nb = 3; }
   else if (b.nbdft()>=2 && b.nbdfm()>=3 && b.nbdfl()>=4) { nb = 4; }
   else if (b.nbdfm()==1) { nb = 1;}
   else { nb = 0; }
   return nb;
 });

 const NamedFunc dphi1("dphi1", [](const Baby &b) -> NamedFunc::ScalarType{
   return fabs(b.jet_phi()->at(b.jet_ordered_pt_indices()->at(0)) - b.met_phi());
 });

 const NamedFunc dphi2("dphi2", [](const Baby &b) -> NamedFunc::ScalarType{
   return fabs(b.jet_phi()->at(b.jet_ordered_pt_indices()->at(1)) - b.met_phi());
 });

 const NamedFunc dphi3("dphi3", [](const Baby &b) -> NamedFunc::ScalarType{
   return fabs(b.jet_phi()->at(b.jet_ordered_pt_indices()->at(2)) - b.met_phi());
 });

 const NamedFunc dphi4("dphi4", [](const Baby &b) -> NamedFunc::ScalarType{
   return fabs(b.jet_phi()->at(b.jet_ordered_pt_indices()->at(3)) - b.met_phi());
 });

 // push all dphi's into a vector, vector size depends on number of jets
 const NamedFunc dphi_vec("dphi_vec", [](const Baby &b) -> NamedFunc::VectorType{
   vector<double> dphi;
   for (int i=0; i < min(b.njet(), 4); i++){
     dphi.push_back( fabs(b.jet_phi()->at(b.jet_ordered_pt_indices()->at(i)) - b.met_phi()) );
   }
   return dphi;
 });
  
}

