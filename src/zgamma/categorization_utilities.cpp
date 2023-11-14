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
#include "core/hist1d.hpp"
#include "core/utilities.hpp"
#include "zgamma/zg_utilities.hpp"
#include "zgamma/categorization_utilities.hpp"

using namespace ZgUtilities;


namespace CatUtilities {

  //This function checks if a deep flavor value passes the loose, medium, or tight working point.
  //The lmt is 0=loose, 1=medium, 2=tight
  bool btag_DF_pass(TString year,float deep_flav,int lmt=1){
    year.ReplaceAll("-","");
    std::map<TString, std::vector<float>> btag_df_wpts{
      {"2016", std::vector<float>({0.0614, 0.3093, 0.7221})},
      {"2016APV", std::vector<float>({0.0614, 0.3093, 0.7221})},
      {"2017", std::vector<float>({0.0521, 0.3033, 0.7489})},
      {"2018", std::vector<float>({0.0494, 0.2770, 0.7264})}
    };
  
    if(deep_flav < btag_df_wpts[year][lmt]){return false;}
    return true;
  }

  //This returns the first jet that passes the bjet_medium criteria
  TLorentzVector bjet_cand(const Baby &b,int lmt = 1){
    TLorentzVector b_jet;b_jet.SetPtEtaPhiM(0,0,0,0);
    for(unsigned int idx = 0; idx < b.jet_pt() -> size(); idx++){
      if(!(b.jet_isgood() -> at(idx))){continue;}
      if(!btag_DF_pass( std::to_string(b.SampleType()) , b.jet_deepflav()->at(idx),lmt ) ){continue;}
      b_jet.SetPtEtaPhiM(b.jet_pt() -> at(idx), b.jet_eta() -> at(idx),b.jet_phi() -> at(idx),b.jet_m() -> at(idx));
      break;
    }
    
    return b_jet;
  }

  //Returns the first two b-jets (highest pT) added together
  TLorentzVector bb_cand(const Baby &b, int lmt=1){
    int count = 0;
    TLorentzVector bb,b_jet;bb.SetPtEtaPhiM(0,0,0,0);


    for(unsigned int idx = 0; idx < b.jet_pt() -> size(); idx++){
      if(!(b.jet_isgood() -> at(idx))){continue;}
      if(!btag_DF_pass( b.SampleTypeString() , b.jet_deepflav()->at(idx),lmt ) ){continue;}
      b_jet.SetPtEtaPhiM(b.jet_pt() -> at(idx), b.jet_eta() -> at(idx),b.jet_phi() -> at(idx),b.jet_m() -> at(idx));
      bb += b_jet;
      count++;
      if(count>1){break;}
    }

    return bb;
  }

  //Assigns the nth jet (nth using the jet_num variable)
  TLorentzVector AssignJet(const Baby &b,int jet_num=1){
     TLorentzVector jet; jet.SetPtEtaPhiM(0,0,0,0);
     int count = 0;
     for(unsigned int i = 0; i < b.jet_pt() -> size(); i++){
       if(!(b.jet_isgood() -> at(i))){continue;}
       if(count < jet_num-1){count++; continue;}
       jet.SetPtEtaPhiM(b.jet_pt() -> at(i),b.jet_eta() -> at(i) ,b.jet_phi() -> at(i), b.jet_m() -> at(i) );
       break;
     }
     return  jet;
  } 

  //Returns the lepton of flavor flav and index idx
  TLorentzVector AssignLep(const Baby &b, int flav, int idx){
    TLorentzVector vec; vec.SetPtEtaPhiM(0,0,0,0);
    if(idx==-1){return vec;}
    if(flav==11){vec.SetPtEtaPhiM(b.el_pt() ->at(idx),b.el_eta() ->at(idx),b.el_phi() ->at(idx), 0.000511);}
    if(flav==13){vec.SetPtEtaPhiM(b.mu_pt() ->at(idx),b.mu_eta() ->at(idx),b.mu_phi() ->at(idx), 0.1057);}
    return vec;
  }

  //This returns the flavor and idces of the leptons that make up the ll not used for the lly candidate
  std::vector<int> Z2_idces(const Baby &b){
    int idx_keep = Z2_idx(b);
    if(idx_keep==-1){ return {-1}; }
    return {b.ll_lepid() -> at(idx_keep), b.ll_i1() -> at(idx_keep),b.ll_i2() -> at(idx_keep)};
  }

  //Returns TLorentzVector of second Z candidate
  TLorentzVector AssignZ2(const Baby &b){
    int idx_z2 = Z2_idx(b);
    TLorentzVector Z2(0,0,0,0); 
    if(idx_z2==-1){ return Z2; }
    Z2.SetPtEtaPhiM(b.ll_pt() ->at(idx_z2), b.ll_eta() ->at(idx_z2), b.ll_phi() ->at(idx_z2), b.ll_m() ->at(idx_z2));
    return Z2;
  }

  //Finds the first electron ignoring the idces indices 
  int el_firstIdx(const Baby &b, std::vector<int> idces={-1,-1}){
    for(int idx = 0; idx < b.nel(); idx++){
      for(auto jdx : idces){ if(idces[jdx]==idx){continue;} }
      if( !(b.el_sig() -> at(idx)) ){continue;}
      return idx;
    }
    return -1;
  }

  //Finds the first muon ignoring the idces indices 
  int mu_firstIdx(const Baby &b, std::vector<int> idces={-1,-1}){
    for(int idx = 0; idx < b.nmu(); idx++){
      for(auto jdx : idces){ if(idces[jdx]==idx){continue;} }
      if( !(b.mu_sig() -> at(idx)) ){continue;}
      return idx;
    }
    return -1;
  }


  //This returns the index of the third lepton (if it is an electron)
  int el3_idx(const Baby &b){
    std::vector<int> l12 = {b.ll_i1() -> at(0),b.ll_i2() -> at(0)};
    if(b.ll_lepid() -> at(0)==13 && b.nel()>0){        return el_firstIdx(b);}
    else if(b.ll_lepid() -> at(0)==11 && b.nel()>2){ return el_firstIdx(b,l12);} 
    return -1;
  }

  //This returns the index of the third lepton (if it is a muon)
  int mu3_idx(const Baby &b){
    std::vector<int> l12 = {b.ll_i1() -> at(0),b.ll_i2() -> at(0)};
    if(b.ll_lepid() -> at(0)==13 && b.nmu()>2){ return mu_firstIdx(b,l12);} 
    else if(b.ll_lepid() -> at(0)==11 && b.nmu() > 0){ return mu_firstIdx(b);}
    return -1;
  }


  //Assigns the third lepton ensuring the leptons from the Z are not selected
  TLorentzVector AssignL3(const Baby &b){
    int l3_idx = -1;
    TLorentzVector l3; l3.SetPtEtaPhiM(0,0,0,0);
    std::vector<int> l12 = {b.ll_i1() -> at(0),b.ll_i2() -> at(0)};

    if(b.nlep() > 3){
      l12 = Z2_idces(b);
      l3 = AssignLep(b,l12[0],l12[1]);
      return l3;
    }

    //This handles the case when there are only 3 leptons and returns third lepton
    //The priority is given to muons
    if(b.ll_lepid() -> at(0)==13){
      if(b.nmu()>2){ l3_idx=mu_firstIdx(b,l12);l3=AssignLep(b,13,l3_idx);} 
      else if(b.nel()>0){l3_idx=el_firstIdx(b);l3=AssignLep(b,11,l3_idx);}
    } else if(b.ll_lepid() -> at(0)==11){
      if(b.nmu()>0){l3_idx=mu_firstIdx(b);l3=AssignLep(b,13,l3_idx);}
      else if(b.nel()>2){ l3_idx=el_firstIdx(b,l12);l3=AssignLep(b,11,l3_idx);} 
    }

    return l3;
  }


  //This returns the index and flavor of the third lepton
  std::vector<int> L3_Idx(const Baby &b){
    std::vector<int> l12 = {b.ll_i1() -> at(0),b.ll_i2() -> at(0)};
    if(b.nlep() > 3){
      l12 = Z2_idces(b);
      return {l12[0],l12[1]};
    }

    if(b.ll_lepid() -> at(0)==13){
      if(b.nmu()>2){ return {13,mu_firstIdx(b,l12)}; } 
      else if(b.nel()>0){ return {11,el_firstIdx(b)}; } 
    } else if(b.ll_lepid() -> at(0)==11){
      if(b.nmu()>0){ return{13,mu_firstIdx(b)}; }
      else if(b.nel()>2){ return {11,el_firstIdx(b,l12)}; }
    }
    return {0,0};
  }

  //Assigns the fourth lepton ensuring the leptons from the Z or AssignL3 function are not selected
  TLorentzVector AssignL4(const Baby &b){
    TLorentzVector l4;
    std::vector<int> l34 = Z2_idces(b);
    l4 = AssignLep(b,l34[0],l34[2]);
    return l4;
  }

  //returns the index and flavor of the 4th lepton
  std::vector<int> L4_Idx(const Baby &b){
    TLorentzVector l4; l4.SetPtEtaPhiM(0,0,0,0);
    std::vector<int> l34 = Z2_idces(b);
    return {l34[0],l34[2]};
  }

  //This takes the previous functions and returns the flavor of the third lepton
  int L3_flav_i(const Baby &b){
    return L3_Idx(b)[0];
  }

  //This takes the previous functions and returns the flavor of the 4th lepton
  int L4_flav_i(const Baby &b){
    return L4_Idx(b)[0];
  }

  //This checks if theh third and fourth lepton match flavors (mostly) deprecated 
  bool Z2_flav_b(const Baby &b){
    return L4_flav_i(b) == L3_flav_i(b);
  }

  //Returns the difference in phi between to values
  float get_dphi(float phi1, float phi2){
     const double PI = 3.1415;
     double dphi = fmod(fabs(phi2-phi1), 2.*PI);
     dphi = dphi>PI ? 2.*PI-dphi : dphi;

     return dphi;
   }

  //Returns Delta Phi(H,MET)
  float dphi_Hmet_d(const Baby &b){
     TLorentzVector H; H = AssignH(b);
     return get_dphi(H.Phi(),b.met_phi());
   }

  //Returns the maximum deepflavor score of all jets in the event
  double max_deepflav_d(const Baby &b){
    float mx = -1;
    for(auto dfv : *b.jet_deepflav() ){
      mx = std::max(mx,dfv); 
    }
    return mx;
  }

  //This function returns the mass of the "best" top candidate
  //First the two most likely b-jets are found by finding max and 2nd max b-jet score
  //Then the best W candidate is found via dijet mass
  //The two b-jets are compared to see which one has the most compatible mass with the top
  double m_top(const Baby &b){
    std::vector <int> b_idx(2);
    std::vector <double> b_scores = {-1,-1};

    TLorentzVector temp1,temp2,W_cand,b_max,b_2max;
    double compDouble = 0;
    for(unsigned int idx_k = 0; idx_k < b.jet_pt() -> size(); idx_k++){
      if( !(b.jet_isgood()->at(idx_k)) ){continue;}
      compDouble = b.jet_deepflav()->at(idx_k);

      if(compDouble > b_scores[1] ){

        if( compDouble > b_scores[0]){
          b_idx[1] = b_idx[0]; b_scores[1] = b_scores[0];
          b_idx[0] = idx_k; b_scores[0] = compDouble;
          continue;
        }
 
        b_idx[1] = idx_k; b_scores[1] = compDouble;
       
      } 
    }

    compDouble = 1000;
    for(unsigned int idx_k = 0; idx_k < b.jet_pt() -> size(); idx_k++){
      if( !(b.jet_isgood()->at(idx_k))  || (static_cast<int>(idx_k)==b_idx[1]) || (static_cast<int>(idx_k)==b_idx[0]) ){continue;}
      temp1.SetPtEtaPhiM( b.jet_pt() -> at(idx_k), b.jet_eta() -> at(idx_k), b.jet_phi() -> at(idx_k), b.jet_m() -> at(idx_k) );

      for(unsigned int idx_j = idx_k+1; idx_j < b.jet_pt() -> size(); idx_j++){
        if( !(b.jet_isgood()->at(idx_j)) || (static_cast<int>(idx_j)==b_idx[1]) || (static_cast<int>(idx_j)==b_idx[0]) ){continue;}
        temp2.SetPtEtaPhiM(b.jet_pt() ->at(idx_j), b.jet_eta() -> at(idx_j), b.jet_phi()->at(idx_j), b.jet_m() -> at(idx_j));
        if( fabs((temp1+temp2).M() - 80.4) <  compDouble){
          W_cand = temp1+temp2;
          compDouble = W_cand.M();           
        }   
      }      
    }
    
    b_max.SetPtEtaPhiM( b.jet_pt() -> at(b_idx[0]), b.jet_eta() -> at(b_idx[0]), b.jet_phi() -> at(b_idx[0]), b.jet_m() -> at(b_idx[0]) );
    b_2max.SetPtEtaPhiM( b.jet_pt() -> at(b_idx[1]), b.jet_eta() -> at(b_idx[1]), b.jet_phi() -> at(b_idx[1]), b.jet_m() -> at(b_idx[1]) );
    double mtop1 = (b_max + W_cand).M();
    double mtop2 = (b_2max + W_cand).M();

    return fabs( mtop1 - 172.7) < fabs( mtop2 - 172.7) ? mtop1 : mtop2;
  }


  //sums jet pT individually
  double jet_pt_SUM(const Baby &b,int nlim=3){
    if(b.njet() < nlim){
      return -100;
    }
    int count = 0;
    double pt_sum = 0; 
    for(unsigned int i = 0; i < b.jet_pt() -> size(); i++){
       if(count>nlim){break;}
       if(!(b.jet_isgood() -> at(i))){continue;}
       pt_sum += b.jet_pt() -> at(i);
       count++;
    }
    return pt_sum;
  }

  //Returns mass of jets added to one TLorentzVector
  double m_jet_SUM(const Baby &b,int nlim=3){
     if(b.njet() < nlim){
      return -1;
     }

     int count = 0;
     TLorentzVector j_sum,j_it;
     j_sum.SetPtEtaPhiM(0,0,0,0);
     for(unsigned int i = 0; i < b.jet_pt() -> size(); i++){
       if(count>nlim){break;}
       if(!(b.jet_isgood() -> at(i))){continue;}
       j_it.SetPtEtaPhiM(b.jet_pt() -> at(i),b.jet_eta() -> at(i) ,b.jet_phi() -> at(i), b.jet_m() -> at(i) );
       j_sum += j_it;
       count++;
     }
     return  j_sum.M();
  }

  //Returns pT sum over m_nj
  double ptj_m4j_d(const Baby &b){
    return jet_pt_SUM(b,4)/m_jet_SUM(b,4);
  }


  //Loop trough ll pairs and finds first pair with both leptons not matching the idces of the leptons used for the lly candidate
  //This should be the ll closest to the mass of the Z separate from our first ll pair
  int Z2_idx(const Baby &b){
    int idx_keep = -1;
    for(unsigned int idx_ll = 0; idx_ll < b.ll_lepid() -> size(); idx_ll++){
    if(b.ll_i1() ->at(idx_ll) == b.ll_i1() ->at(0) || b.ll_i2() ->at(idx_ll) == b.ll_i2() -> at(0)){continue;}
      idx_keep = idx_ll; break;
    }
    return idx_keep;
  }

  //This function returns the pT of the second Z
  double Z2_pt_d(const Baby &b){ return b.ll_pt()->at(Z2_idx(b)); }
  double Z2_m_d(const Baby &b){  return b.ll_m()->at(Z2_idx(b)); }

  //Checks the Delta R separation between the second Z cand. and the photon used in the lly cand. 
  double Z2y_dR_d(const Baby &b){
    TLorentzVector Z2,y;
    y = AssignGamma(b);
    Z2 = AssignZ2(b); 
    return y.DeltaR(Z2);
  }

  //This returns Delta R between Z candidates
  double ZZ_dR_d(const Baby &b){
    TLorentzVector ll1,ll2;
    ll1 = AssignZ(b);
    ll2 = AssignZ2(b);
    return ll1.DeltaR(ll2);
  }

  //Returns the pt of the third lepton 
  double l3_pt_d(const Baby &b){
    TLorentzVector l3;
    l3 = AssignL3(b);
    return l3.Pt();
  }

  //Returns the pT of the fourth lepton
  double l4_pt_d(const Baby &b){
    TLorentzVector l4;
    l4 = AssignL4(b);
    return l4.Pt();
  }

  //Returns the mll for the third lepton and the leptons used as the ll candidate
  double l3_ZCheck_d(const Baby &b){
    TLorentzVector l1,l2,l3;
    l1 = AssignL1(b);
    l2 = AssignL2(b);
    l3 = AssignL3(b);
    return std::min( abs((l1+l3).M()-91.2) , abs((l2+l3).M() - 91.2) );
  }

  //Returns Delta R separation between the Z and the third lepton
  double l3_Z_DR_d(const Baby &b){
    TLorentzVector Z, l3;
    Z  = AssignZ(b); 
    l3 = AssignL3(b);
    return Z.DeltaR(l3);
  }

  //Returns Delta R between the fourth lepton and the Z candidate
  double l4_Z_DR_d(const Baby &b){
    TLorentzVector Z, l4;
    Z  = AssignZ(b); 
    l4 = AssignL4(b);
    return Z.DeltaR(l4);
  }

  //Returns the el_idmva score for the third lepton if it is an electron 
  //Note: for muons this returns 2. 
  double l3_idmva_d(const Baby &b){
    std::vector<int> l3_i = L3_Idx(b);
    if(l3_i[0] == 13 || l3_i[1]==-1){ return 2; }
    return b.el_idmva() -> at(l3_i[1]);
  }

  //Returns the el_idmva score for the fourth lepton if it is an electron 
  //Note: for muons this returns 2. 
  double l4_idmva_d(const Baby &b){
    std::vector<int> l4_i = L4_Idx(b);
    if(l4_i[0] == 13 || l4_i[1]==-1){ return 2; }
    return b.el_idmva() -> at(l4_i[1]);
  }

  //This function returns mini iso. for the third lepton candidate
  double l3_Imini_d(const Baby &b){
    std::vector<int> idx_l3 = L3_Idx(b);
    if(idx_l3[1]==-1){return 100;}
    if(idx_l3[0]==13){ return b.mu_miniso() -> at(idx_l3[1]); }
    if(idx_l3[0]==11){ return b.el_miniso() -> at(idx_l3[1]); }
    return 100;
  }

  //This function returns the mini iso. for the 4th lepton candidate
  double l4_Imini_d(const Baby &b){
    std::vector<int> idx_l4 = L4_Idx(b);
    if(idx_l4[1]==-1){return 100;}
    if(idx_l4[0]==13){ return b.mu_miniso() -> at(idx_l4[1]); }
    if(idx_l4[0]==11){ return b.el_miniso() -> at(idx_l4[1]); }
    return 100;
  }

  double lly_ptom_d(const Baby &b){ 
    TLorentzVector H = AssignH(b);
    return H.Pt()/(H.M());
  }


  const NamedFunc mbb("mbb",[](const Baby &b) -> NamedFunc::ScalarType{return bb_cand(b,0).M();});
  const NamedFunc mtop("mtop",[](const Baby &b) -> NamedFunc::ScalarType{return m_top(b);});
  const NamedFunc lly_pT("NF_lly_pT",[](const Baby &b) -> NamedFunc::ScalarType{return b.llphoton_pt()->at(0);});
  const NamedFunc lly_dR("lly_dR",[](const Baby &b) -> NamedFunc::ScalarType{return (AssignZ(b).DeltaR(AssignGamma(b)));});
  const NamedFunc pTlly_mlly("pTlly_mlly",[](const Baby &b) -> NamedFunc::ScalarType{return lly_ptom_d(b);});

  const NamedFunc l3_pt("l3_pt",[](const Baby &b) -> NamedFunc::ScalarType{return l3_pt_d(b);});
  const NamedFunc l3_idmva("l3_idmva",[](const Baby &b) -> NamedFunc::ScalarType{return l3_idmva_d(b);});
  const NamedFunc l3_mini("l3_mini",[](const Baby &b) -> NamedFunc::ScalarType{return l3_Imini_d(b);});
  const NamedFunc l3_flav("l3_flav",[](const Baby &b) -> NamedFunc::ScalarType{return L3_flav_i(b);});
  const NamedFunc l4_pt("l4_pt",[](const Baby &b) -> NamedFunc::ScalarType{return l4_pt_d(b);});
  const NamedFunc l4_idmva("l4_idmva",[](const Baby &b) -> NamedFunc::ScalarType{return l4_idmva_d(b);});
  const NamedFunc l4_mini("l4_mini",[](const Baby &b) -> NamedFunc::ScalarType{return l4_Imini_d(b);});
  const NamedFunc l4_flav("l4_flav",[](const Baby &b) -> NamedFunc::ScalarType{return L4_flav_i(b);});

  const NamedFunc l3_ZCheck("l3_ZCheck",[](const Baby &b) -> NamedFunc::ScalarType{return l3_ZCheck_d(b);});
  const NamedFunc l3_Z_DR("l3_Z_DR",[](const Baby &b) -> NamedFunc::ScalarType{return l3_Z_DR_d(b);});

  const NamedFunc ptj_m4j("ptj_m4j",[](const Baby &b) -> NamedFunc::ScalarType{return ptj_m4j_d(b);});
  const NamedFunc max_deepflav("max_deepflav",[](const Baby &b) -> NamedFunc::ScalarType{return max_deepflav_d(b);});
  const NamedFunc Z2_pt("Z2_pt",[](const Baby &b) -> NamedFunc::ScalarType{return Z2_pt_d(b);});
  const NamedFunc Z2_m("Z2_m",[](const Baby &b) -> NamedFunc::ScalarType{return Z2_m_d(b);});
  const NamedFunc Z2_flav("Z2_flav",[](const Baby &b) -> NamedFunc::ScalarType{return Z2_flav_b(b);});
  const NamedFunc Z2y_dR("Z2y_dR",[](const Baby &b) -> NamedFunc::ScalarType{return Z2y_dR_d(b);});
  const NamedFunc ZZ_dR("ZZ_dR",[](const Baby &b) -> NamedFunc::ScalarType{return ZZ_dR_d(b);});
  const NamedFunc dphi_Hmet("dphi_Hmet",[](const Baby &b) -> NamedFunc::ScalarType{return dphi_Hmet_d(b);});
 

  //This function is the method Run 2 used for categorization
  int R2_Cats(const Baby &b){ 
    if(b.nlep() > 2){ return 1;}
    if(b.njet() > 1){ return 2;} 
    return 3;
  }

  //NamedFuncs returning whether the event belongs to the respective Run 2 category
  const NamedFunc run2_lep("run2_lep", [](const Baby &b) -> NamedFunc::ScalarType{ return R2_Cats(b)==1; });
  const NamedFunc run2_VBF("run2_VBF", [](const Baby &b) -> NamedFunc::ScalarType{ return R2_Cats(b)==2; });
  const NamedFunc run2_ggF("run2_ggF", [](const Baby &b) -> NamedFunc::ScalarType{ return R2_Cats(b)==3; });

  //Vector containing the Run 2 categories
  const std::vector<NamedFunc> run2_categories = {run2_ggF, run2_VBF, run2_lep};

  //This function is the object selections used for the categories
  int cat_obj_sels(const Baby &b){
    if( b.njet()>0 && b.nlep()>3  && b.nbdfm()>0){  return 0; }
    if( b.njet()>2 && b.nlep()==3 && b.nbdfm()>1){  return 1; }
    if( b.njet()>2 && b.nlep()==3 && b.nbdfm()==1){ return 2; }
    if( b.nlep()>3  ){ return 3; }
    if( b.nlep()==3 ){ return 4; }

    if( b.njet()>4 && b.nbdfm()>1 ){  return 5; }
    if( b.njet()>4 && b.nbdfm()==1 ){ return 6; }
    if( b.njet()>1 && b.nbdfl()==2){  return 7; }
    if( b.met() > 90){  return 8; }
    if( b.njet()>1){ return 9;}
    return 10;
  }

  //These are functions for categories from the object selection above
  const NamedFunc cat_ttH_lep_4l("cat_ttH_lep_4l",[](const Baby &b) -> NamedFunc::ScalarType{      return cat_obj_sels(b)==0; });
  const NamedFunc cat_ttH_lep_3l2b("cat_ttH_lep_3l2b",[](const Baby &b) -> NamedFunc::ScalarType{  return cat_obj_sels(b)==1; });
  const NamedFunc cat_ttH_lep_3l1b("cat_ttH_lep_3l1b",[](const Baby &b) -> NamedFunc::ScalarType{  return cat_obj_sels(b)==2; });
  const NamedFunc cat_VH_4l("cat_VH_4l",[](const Baby &b) -> NamedFunc::ScalarType{                return cat_obj_sels(b)==3; });
  const NamedFunc cat_VH_3l("cat_VH_3l",[](const Baby &b) -> NamedFunc::ScalarType{                return cat_obj_sels(b)==4; });
  const NamedFunc cat_ttH_had_2b("cat_ttH_had_2b",[](const Baby &b) -> NamedFunc::ScalarType{      return cat_obj_sels(b)==5; });
  const NamedFunc cat_ttH_had_1b("cat_ttH_had_1b",[](const Baby &b) -> NamedFunc::ScalarType{      return cat_obj_sels(b)==6; });
  const NamedFunc cat_ZH_had("cat_ZH_had",[](const Baby &b) -> NamedFunc::ScalarType{              return cat_obj_sels(b)==7; });
  const NamedFunc cat_ZH_met("cat_ZH_met",[](const Baby &b) -> NamedFunc::ScalarType{              return cat_obj_sels(b)==8; });
  const NamedFunc cat_VBF("cat_VBF",[](const Baby &b) -> NamedFunc::ScalarType{                    return cat_obj_sels(b)==9; });
  const NamedFunc cat_ggF("cat_ggF",[](const Baby &b) -> NamedFunc::ScalarType{                    return cat_obj_sels(b)==10; });

  //Vector returning all the categories used for Run 3
  const std::vector<NamedFunc> run3_category_vector = {(cat_ttH_lep_4l || cat_ttH_lep_3l2b || cat_ttH_lep_3l1b) && "ht > 160" && l3_mini < 0.1, 
                                                       cat_VH_4l && "ll_m[0] > 85 && ll_m[0] < 95" &&  Z2_flav  && l3_mini < 0.1 && l4_mini < 0.1 && Z2_pt > 50, 
                                                       cat_VH_3l && "met > 50" && l3_pt > 10 && dphi_Hmet > 2.0,
                                                       ((cat_ttH_had_2b && max_deepflav > 0.75) || cat_ttH_had_1b) && "ht > 250",
                                                       cat_ZH_had && "photon_pt[0] > 20" && max_deepflav > 0.8 && "photon_drmin[0] < 2.0",
                                                       cat_ZH_met && dphi_Hmet > 2.5 && pTlly_mlly > 0.6 && "nbdfl==0",
                                                       cat_VBF && "dijet_m > 500",
                                                       cat_VBF && "dijet_m < 500",
                                                       cat_ggF};
  const std::vector<std::string> run3_category_labels = {"_cat_ttH_lep","_cat_VH_4l","_cat_VH_3l","_cat_ttH_had","_cat_VH_2bl","_cat_VH_MET","_cat_VBF_hmjj","_cat_VBF_lmjj","_cat_ggF"};                          


  //Some repeat functions not used outside of categorization utilities  
  const NamedFunc pTy_mlly("pTy_mlly",[](const Baby &b) -> NamedFunc::ScalarType{ return (b.photon_pt()->at(0))/ZgUtilities::Findlly(b); });
  const NamedFunc mlly("mlly",[](const Baby &b) -> NamedFunc::ScalarType{   return AssignH(b).M(); });
  const NamedFunc pT_lly("pT_lly",[](const Baby &b) -> NamedFunc::ScalarType{   return AssignH(b).Pt(); });
  const NamedFunc eta_lly("eta_lly",[](const Baby &b) -> NamedFunc::ScalarType{ return AssignH(b).Eta(); });
  const NamedFunc phi_lly("phi_lly",[](const Baby &b) -> NamedFunc::ScalarType{ return AssignH(b).Phi(); });

  //Photon plots added to each category's plotting function
  void sample_photon_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    pm.Push<Hist1D>(Axis(45,0.0,1,    "photon_idmva[0]",              "#gamma - IDMVA",           {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_phIDMVA__");
    pm.Push<Hist1D>(Axis(35,0,3.5,    "photon_drmin[0]",              "min( #Delta R(#gamma,l) )",{}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_minDR_");
    pm.Push<Hist1D>(Axis(20,0.01,0.25,"photon_pterr[0]/photon_pt[0]", "#sigma(#gamma)/E(#gamma)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_relunc_");
    pm.Push<Hist1D>(Axis(40,10,60,    "photon_pt[0]",                 "p_{T}(#gamma)",            {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_pT_");
    pm.Push<Hist1D>(Axis(40,-2.6,2.6, "photon_eta[0]",                "#eta(#gamma)",             {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_eta_");
    pm.Push<Hist1D>(Axis(40,-2.6,2.6, "photon_phi[0]",                "#phi(#gamma)",             {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_phi_");
  }


  //ll and llphoton plots added to each category's photon plots
  void sample_llphoton_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    pm.Push<Hist1D>(Axis(40,100,180, mlly,      "m_{ll#gamma} [GeV]", {122,128}), selection, processes ,ops).Weight(wgt).Tag("ShortName:" + labels + "_lly_m_");
    pm.Push<Hist1D>(Axis(50,0,250,   pT_lly,    "p_{T}(ll#gamma)", {}),           selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_lly_pt_");
    pm.Push<Hist1D>(Axis(52,-2.6,2.6,   eta_lly,   "#eta(ll#gamma)", {}),            selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_lly_eta_");
    pm.Push<Hist1D>(Axis(63,-3.15,3.15,   phi_lly,   "#phi(ll#gamma)", {}),            selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_lly_phi_");
    pm.Push<Hist1D>(Axis(40,0,4.5,   pTlly_mlly,"p_{T}(ll#gamma)/m_{ll#gamma}", {1}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_lly_ptom_");
    pm.Push<Hist1D>(Axis(40,0,3.0,   pTy_mlly,  "p_{T}(#gamma)/m_{ll#gamma}", {}),    selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_y_ptom_");
    pm.Push<Hist1D>(Axis(45,0,4.5,   lly_dR,    "#Delta R(ll,#gamma)", {}),           selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ll_y_dR_");

    pm.Push<Hist1D>(Axis(50,0,100, "ll_pt[0]",  "p_{T}(ll)",{}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ll_pt_");
    pm.Push<Hist1D>(Axis(52,-2.6,2.6, "ll_eta[0]", "#eta(ll)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ll_eta_");
    pm.Push<Hist1D>(Axis(63,-3.15,3.15, "ll_phi[0]", "#phi(ll)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ll_phi_");
    pm.Push<Hist1D>(Axis(50,0,100, "ll_m[0]",   "m_{ll}",   {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ll_m_");

  }


  //Sample ttH leptonic category plots
  void ttH_lep_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    sample_photon_plots(pm, selection, processes, ops, wgt, "run3cat_ttH_lep");
    sample_llphoton_plots(pm, selection, processes, ops, wgt, "run3cat_ttH_lep");

    pm.Push<Hist1D>(Axis(80,0,500, "ht",    "H_{T}", {120}),  selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_ttH_lep_HT_" + labels); 
    pm.Push<Hist1D>(Axis(75,0,200, "met",   "MET", {}),       selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_ttH_lep_MET" + labels);
    pm.Push<Hist1D>(Axis(10,0,10,  "njet",  "N_{jet}", {3}),  selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_ttH_lep_Njet_" + labels);
    pm.Push<Hist1D>(Axis(10,0,10,  "nlep",  "N_{lep}", {}),   selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_ttH_lep_Nlep_" + labels);
    pm.Push<Hist1D>(Axis(10,0,6 ,  "nbdfm", "N_{b,med}", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_ttH_3l2b_Nbdfm" + labels);

    pm.Push<Hist1D>(Axis(50,0,1,   l3_mini, "I_{mini} - l3 [GeV]", {}), selection,                processes, ops).Weight(wgt).Tag("ShortName:run3cat_ttH_lep_Imini" + labels);
    pm.Push<Hist1D>(Axis(55,-0.1,1,l3_idmva, "l_{3} - e IDMVA ", {}),   selection && l3_flav==11, processes, ops).Weight(wgt).Tag("ShortName:run3cat_ttH_lep_e3idmva" + labels); 
    pm.Push<Hist1D>(Axis(35,5,75,  l3_pt, "p_{T}(l_{3})", {}),          selection && l3_flav==11, processes, ops).Weight(wgt).Tag("ShortName:run3cat_ttH_lep_e3_pt" + labels);
    pm.Push<Hist1D>(Axis(35,5,75,  l3_pt, "p_{T}(l_{3})", {}),          selection && l3_flav==13, processes, ops).Weight(wgt).Tag("ShortName:run3cat_ttH_lep_mu3_pt" + labels);
  }

  //Sample VH 4l category plots
  void VH_4l_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    sample_photon_plots(pm, selection, processes, ops, wgt, "run3cat_VH_4l");
    sample_llphoton_plots(pm, selection, processes, ops, wgt, "run3cat_VH_4l");

    pm.Push<Hist1D>(Axis(50,0,1,    l3_mini,  "I_{mini} - l3 [GeV]", {}),     selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_4l_Imini_3l" + labels);
    pm.Push<Hist1D>(Axis(35,5,75,   l3_pt,    "p_{T}(e_{3})", {}),            selection && l3_flav==11, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_4l_e3_pt" + labels);
    pm.Push<Hist1D>(Axis(35,5,75,   l3_pt,    "p_{T}(#mu_{3})", {}),          selection && l3_flav==13, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_4l_mu3_pt" + labels);
    pm.Push<Hist1D>(Axis(55,-0.1,1, l3_idmva, "l_{3} - e IDMVA ", {}),        selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_4l_e3idmva" + labels);

    pm.Push<Hist1D>(Axis(50,0,1,    l4_mini, "I_{mini} - l4 [GeV]", {}),      selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_4l_Imini_4l" + labels);
    pm.Push<Hist1D>(Axis(35,5,75,   l4_pt,   "p_{T}(e_{4})", {}),             selection && l4_flav==11, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_4l_e4_pt" + labels);
    pm.Push<Hist1D>(Axis(35,5,75,   l4_pt,   "p_{T}(#mu_{4})", {}),           selection && l4_flav==13, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_4l_mu4_pt" + labels);
    pm.Push<Hist1D>(Axis(55,-0.1,1, l4_idmva, "l_{4} - e IDMVA ", {}),        selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_4l_e4idmva" + labels); 

    pm.Push<Hist1D>(Axis(80,0,400,  Z2_pt,  "p_{T}(ll) - Z2", {50}),          selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_4l_ptll_Z2" + labels);
    pm.Push<Hist1D>(Axis(70,50,120, Z2_m,   "m_{l_3 l_4} - Z2", {80,100}),    selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_4l_mll_Z2" + labels); 
    pm.Push<Hist1D>(Axis(50,0.0,5,  Z2y_dR, "#Delta R(Z_{2},#gamma)", {0.2}), selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_4l_Z2y_dr" + labels);
    pm.Push<Hist1D>(Axis(50,0.0,5,  ZZ_dR,  "#Delta R(Z_{1},Z_{2})", {}),     selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_4l_ZZ_dr" + labels);
  }

  //Sample VH 3l category plots
  void VH_3l_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    sample_photon_plots(pm, selection, processes, ops, wgt, "run3cat_VH_3l");
    sample_llphoton_plots(pm, selection, processes, ops, wgt, "run3cat_VH_3l");

    pm.Push<Hist1D>(Axis(30,0,60,   l3_ZCheck, "min(m_{l_il_j} - m_{Z}) [GeV]", {45}), selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_3l_ZCheck" + labels);
    pm.Push<Hist1D>(Axis(40,0,4.0,  l3_Z_DR,   "#Delta R(Z,l_{3})", {45}),             selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_3l_Zl3_DR" + labels);
    pm.Push<Hist1D>(Axis(50,0,1,    l3_mini,   "I_{mini} - l3 [GeV]", {}),             selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_3l_lep_l3_Imini" + labels);
    pm.Push<Hist1D>(Axis(55,-0.1,1, l3_idmva,  "l_{3} - e IDMVA ", {}),                selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_3l_e3idmva" + labels);
    pm.Push<Hist1D>(Axis(35,5,75,   l3_pt,     "p_{T}(e_{3})", {}),                    selection && l3_flav==11, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_3l_e3_pt" + labels);
    pm.Push<Hist1D>(Axis(35,5,75,   l3_pt,     "p_{T}(#mu_{3})", {}),                  selection && l3_flav==13, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_3l_mu3_pt" + labels);

    pm.Push<Hist1D>(Axis(45,0,4.5,  lly_dR,    "#Delta R(ll,#gamma)", {}),             selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_3l_ll_y_dR" + labels);

    pm.Push<Hist1D>(Axis(75,0,200,  "met",     "MET", {}),                             selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_3l_MET" + labels);
    pm.Push<Hist1D>(Axis(40,0,3.2,  dphi_Hmet, "#Delta #phi(MET,H)", {2.0}),           selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_3l_dphi_Hmet" + labels);
  }

  //Sample ttH had category plots
  void ttH_had_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    sample_photon_plots(pm, selection, processes, ops, wgt, "run3cat_ttH_had");
    sample_llphoton_plots(pm, selection, processes, ops, wgt, "run3cat_ttH_had");

    pm.Push<Hist1D>(Axis(75,0,200,   "met", "MET", {}),                                 selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_ttH_had_2b_MET" + labels);
    pm.Push<Hist1D>(Axis(80,0,500,   "ht", "H_{T}", {150}),                             selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_ttH_had_2b_HT_" + labels);

    pm.Push<Hist1D>(Axis(10,0,10,    "njet", "N_{jet}", {4}),                           selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_ttH_had_2b_Njet_" + labels);
    pm.Push<Hist1D>(Axis(50,100,300, mtop, "cand. m_{t} [GeV]", {130,220}),          selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_ttH_had_2b_mtop" + labels);
    pm.Push<Hist1D>(Axis(50,0,1,     ptj_m4j, "#Sigma(p_{T}(j_{i}))/m_{4j} [GeV]", {0.6}), selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_ttH_had_2b_pt_m4j" + labels);
    pm.Push<Hist1D>(Axis(40,0,1,     max_deepflav, "Deep Flavor - max", {}),         selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_ttH_had_2b_deepflav_BASE_MAX" + labels);
  }

  //Sample VH 2bl category plots
  void VH_2bl_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    sample_photon_plots(pm, selection, processes, ops, wgt, "run3cat_VH_2bl");
    sample_llphoton_plots(pm, selection, processes, ops, wgt, "run3cat_VH_2bl");

    pm.Push<Hist1D>(Axis(10,0,10,  "njet",        "N_{jet}", {2}),            selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_2bl_Njet_" + labels);
    pm.Push<Hist1D>(Axis(10,0,10,  "nbdfm",       "N_{bdfm}", {2}),           selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_2bl_Nbdfm_" + labels); 

    pm.Push<Hist1D>(Axis(40,0,1,    max_deepflav, "Max. Deep Flavor", {0.4}), selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_2bl_deepflav_BASE_MAX" + labels);
    pm.Push<Hist1D>(Axis(50,50,150, mbb,          "m_{bb}", {}),              selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_2bl_mbb" + labels);
    pm.Push<Hist1D>(Axis(75,0,200, "met",         "p_{T}^{miss}", {}),        selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_2bl_MET" + labels);
    pm.Push<Hist1D>(Axis(80,0,500, "ht",          "H_{T}", {}),               selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_2bl_HT_" + labels); 
  }

  //Sample VH MET category plots
  void VH_MET_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    sample_photon_plots(pm, selection, processes, ops, wgt, "run3cat_VH_MET");
    sample_llphoton_plots(pm, selection, processes, ops, wgt, "run3cat_VH_MET");

    pm.Push<Hist1D>(Axis(6,0,6,     "njet",    "N_{jet}", {}),               selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_MET_njet_" + labels);
    pm.Push<Hist1D>(Axis(6,0,6,     "nbdfm",   "N_{b,med}", {}),             selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_MET_nbdfm_" + labels);
    pm.Push<Hist1D>(Axis(6,0,6,     "nbdfl",   "N_{b,loose}", {}),           selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_MET_nbdfl_" + labels);

    pm.Push<Hist1D>(Axis(55,90,200, "met",     "p_{T}^{miss}", {}),          selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_MET_MET_" + labels);
    pm.Push<Hist1D>(Axis(40,0,3.2,  dphi_Hmet, "#Delta #phi(MET,H)", {2.5}), selection, processes, ops).Weight(wgt).Tag("ShortName:run3cat_VH_MET_dphi_Hmet" + labels);
  }

}


