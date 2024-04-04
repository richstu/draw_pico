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

  std::vector<int> Z_avg_idces(const Baby &b){
    int Z1 = -1; int Z2 = -1;
    double min_avg = 10000;
    TLorentzVector ll1, ll2;
    for(unsigned int idx_ll1 = 0; idx_ll1 < b.ll_lepid() ->size(); idx_ll1++){
      for(unsigned int idx_ll2 = idx_ll1+1; idx_ll2 < b.ll_lepid() ->size(); idx_ll2++){
        if(b.ll_lepid()->at(idx_ll1) == b.ll_lepid()->at(idx_ll2)){
          if(b.ll_i1() ->at(idx_ll1) == b.ll_i1() ->at(idx_ll2)){continue;}
          if(b.ll_i2() ->at(idx_ll1) == b.ll_i2() ->at(idx_ll2)){continue;}
          if(b.ll_i2() ->at(idx_ll1) == b.ll_i1() ->at(idx_ll2)){continue;}
          if(b.ll_i1() ->at(idx_ll1) == b.ll_i2() ->at(idx_ll2)){continue;}
        }

        ll1.SetPtEtaPhiM(b.ll_pt() ->at(idx_ll1), b.ll_eta() ->at(idx_ll1), b.ll_phi() ->at(idx_ll1), b.ll_m() ->at(idx_ll1));
        ll2.SetPtEtaPhiM(b.ll_pt() ->at(idx_ll2), b.ll_eta() ->at(idx_ll2), b.ll_phi() ->at(idx_ll2), b.ll_m() ->at(idx_ll2));

        if( fabs((ll1.M() + ll2.M())/2.0 - 91.2) < min_avg){
          Z1 = idx_ll1; Z2 = idx_ll2; min_avg = 10000;
        }
      }
    }

    TLorentzVector gam = AssignGamma(b);
    if(Z1==-1 || Z2==-1){return {-1,-1};}
    ll1.SetPtEtaPhiM(b.ll_pt() ->at(Z1), b.ll_eta() ->at(Z1), b.ll_phi() ->at(Z1), b.ll_m() ->at(Z1));
    ll2.SetPtEtaPhiM(b.ll_pt() ->at(Z2), b.ll_eta() ->at(Z2), b.ll_phi() ->at(Z2), b.ll_m() ->at(Z2));
    std::vector<int> ll_vec = (fabs(ll1.M() - 91.2) < fabs(ll2.M() - 91.2)) ? std::vector<int>({Z1,Z2}) : std::vector<int>({Z2,Z1});

    //(ll1.DeltaR(gam) < ll2.DeltaR(gam)) ? std::vector<int>({Z1,Z2}) : std::vector<int>({Z2,Z1});
    return ll_vec;
  }

  TLorentzVector Z1_avg(const Baby &b){ 
    int ll_Z1 = Z_avg_idces(b)[0];
    TLorentzVector Z1(0,0,0,0); 
    if(ll_Z1==-1){return Z1;}
    Z1.SetPtEtaPhiM(b.ll_pt() ->at(ll_Z1), b.ll_eta() ->at(ll_Z1), b.ll_phi() ->at(ll_Z1), b.ll_m() ->at(ll_Z1));
    return Z1;
  }

  TLorentzVector Z2_avg(const Baby &b){ 
    int ll_Z2 = Z_avg_idces(b)[0];
    TLorentzVector Z2(0,0,0,0);
    if(ll_Z2==-1){return Z2;}
    Z2.SetPtEtaPhiM(b.ll_pt() ->at(ll_Z2), b.ll_eta() ->at(ll_Z2), b.ll_phi() ->at(ll_Z2), b.ll_m() ->at(ll_Z2));
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
      if(l12[0]==-1){return l3;}
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
      if(l12[0]==-1){return {-1,-1};}
      return {l12[0],l12[1]};
    }

    if(b.ll_lepid() -> at(0)==13){
      if(b.nmu()>2){ return {13,mu_firstIdx(b,l12)}; } 
      else if(b.nel()>0){ return {11,el_firstIdx(b)}; } 
    } else if(b.ll_lepid() -> at(0)==11){
      if(b.nmu()>0){ return{13,mu_firstIdx(b)}; }
      else if(b.nel()>2){ return {11,el_firstIdx(b,l12)}; }
    }
    return {-1,-1};
  }

  //Assigns the fourth lepton ensuring the leptons from the Z or AssignL3 function are not selected
  TLorentzVector AssignL4(const Baby &b){
    TLorentzVector l4; l4.SetPtEtaPhiM(0,0,0,0);
    std::vector<int> l34 = Z2_idces(b);
    if(l34[0]==-1){return l4;}
    l4 = AssignLep(b,l34[0],l34[2]);
    return l4;
  }

  //returns the index and flavor of the 4th lepton
  std::vector<int> L4_Idx(const Baby &b){
    TLorentzVector l4; l4.SetPtEtaPhiM(0,0,0,0);
    std::vector<int> l34 = Z2_idces(b);
    if(l34[0]==-1){return {-1,-1};}
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
      if((static_cast<int>(idx_k)==b_idx[1]) || (static_cast<int>(idx_k)==b_idx[0]) ){continue;}
      if( !(b.jet_isgood()->at(idx_k))){continue;}
      temp1.SetPtEtaPhiM( b.jet_pt() -> at(idx_k), b.jet_eta() -> at(idx_k), b.jet_phi() -> at(idx_k), b.jet_m() -> at(idx_k) );

      for(unsigned int idx_j = idx_k+1; idx_j < b.jet_pt() -> size(); idx_j++){
        if((static_cast<int>(idx_k)==b_idx[1]) || (static_cast<int>(idx_k)==b_idx[0]) ){continue;}
        if( !(b.jet_isgood()->at(idx_j)) ){continue;}
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

  //This function returns the mass of the "best" top candidate where there are 3 leptons
  //First the most likely b-jets are found by finding max and 2nd max b-jet score
  //Then the best W candidate is found via dijet mass
  //The two b-jets are compared to see which one has the most compatible mass with the top
  double m_top_l3(const Baby &b){
    int b_idx = -1;
    double b_sc_max = -1.0;

    TLorentzVector temp1,temp2,W_cand,b_max;
    double compDouble = 0;
    for(unsigned int idx_k = 0; idx_k < b.jet_pt() -> size(); idx_k++){
      if( !(b.jet_isgood()->at(idx_k)) ){continue;}
      compDouble = b.jet_deepflav()->at(idx_k);

      if( compDouble > b_sc_max){
        b_sc_max = compDouble; b_idx = idx_k;
      }
    }

    compDouble = 1000;
    for(unsigned int idx_k = 0; idx_k < b.jet_pt() -> size(); idx_k++){
      if( !(b.jet_isgood()->at(idx_k))  || (static_cast<int>(idx_k)==b_idx) ){continue;}
      temp1.SetPtEtaPhiM( b.jet_pt() -> at(idx_k), b.jet_eta() -> at(idx_k), b.jet_phi() -> at(idx_k), b.jet_m() -> at(idx_k) );

      for(unsigned int idx_j = idx_k+1; idx_j < b.jet_pt() -> size(); idx_j++){
        if( !(b.jet_isgood()->at(idx_j)) || (static_cast<int>(idx_j)==b_idx) ){continue;}
        temp2.SetPtEtaPhiM(b.jet_pt() ->at(idx_j), b.jet_eta() -> at(idx_j), b.jet_phi()->at(idx_j), b.jet_m() -> at(idx_j));
        if( fabs((temp1+temp2).M() - 80.4) <  compDouble){
          W_cand = temp1+temp2;
          compDouble = W_cand.M();           
        }   
      }      
    }
    
    b_max.SetPtEtaPhiM( b.jet_pt() -> at(b_idx), b.jet_eta() -> at(b_idx), b.jet_phi() -> at(b_idx), b.jet_m() -> at(b_idx) );
    return (b_max + W_cand).M();
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

  double ptj_m5j_d(const Baby &b){
    return jet_pt_SUM(b,5)/m_jet_SUM(b,5);
  }


  //Loop trough ll pairs and finds first pair with both leptons not matching the idces of the leptons used for the lly candidate
  //This should be the ll closest to the mass of the Z separate from our first ll pair
  int Z2_idx(const Baby &b){
    int idx_keep = -1;
    for(unsigned int idx_ll = 0; idx_ll < b.ll_lepid() -> size(); idx_ll++){
    if((b.ll_lepid() -> at(idx_ll) && b.ll_lepid() -> at(0)) &&
       (b.ll_i1() ->at(idx_ll) == b.ll_i1() ->at(0) || b.ll_i2() ->at(idx_ll) == b.ll_i2() -> at(0)) ){continue;}
      idx_keep = idx_ll; break;
    }
    return idx_keep;
  }

int N_SF_leppairs(const Baby &b){
  //int two = 2;
  return b.nel()/2 + b.nmu()/2;
}

int N_OSSF_leppairs(const Baby &b){
  int os_lep_pairs = 0;
  std::unordered_map<int,bool> el_idces = {};
  std::unordered_map<int,bool> mu_idces = {};

  for(unsigned int idx_ll=0; idx_ll < b.ll_pt() ->size(); idx_ll++){
    if(b.ll_lepid() ->at(idx_ll)==11){
      if(el_idces.find(b.ll_i1()->at(idx_ll))!=el_idces.end()){continue;}
      if(el_idces.find(b.ll_i2()->at(idx_ll))!=el_idces.end()){continue;}
      if( !(b.el_sig()->at(b.ll_i1()->at(idx_ll))) || !(b.el_sig()->at(b.ll_i2()->at(idx_ll))) ){continue;}
      if( b.el_charge()->at(b.ll_i1()->at(idx_ll)) + b.el_charge()->at(b.ll_i2()->at(idx_ll)) != 0){continue;}
      os_lep_pairs++;
      el_idces.insert({b.ll_i1()->at(idx_ll),true});
      el_idces.insert({b.ll_i2()->at(idx_ll),true});
    }

    if(b.ll_lepid() ->at(idx_ll)==13){
      if(mu_idces.find(b.ll_i1()->at(idx_ll))!=mu_idces.end()){continue;}
      if(mu_idces.find(b.ll_i2()->at(idx_ll))!=mu_idces.end()){continue;}
      if( !(b.mu_sig()->at(b.ll_i1()->at(idx_ll))) || !(b.mu_sig()->at(b.ll_i2()->at(idx_ll))) ){continue;}
      if( b.mu_charge()->at(b.ll_i1()->at(idx_ll)) + b.mu_charge()->at(b.ll_i2()->at(idx_ll)) != 0){continue;}
      os_lep_pairs++;
      mu_idces.insert({b.ll_i1()->at(idx_ll),true});
      mu_idces.insert({b.ll_i2()->at(idx_ll),true});
    }
  }

  return os_lep_pairs;
}


  //This function returns the pT of the second Z
  double Z2_pt_d(const Baby &b){ int idx_2 = Z2_idx(b); return idx_2 != -1 ? b.ll_pt()->at(idx_2) : -1; }
  double Z2_m_d(const Baby &b){  int idx_2 = Z2_idx(b); return idx_2 != -1 ? b.ll_m()->at(idx_2)  : -1; }

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
  double max_Imini_d(const Baby &b){
    double comp = -1.0;
    for(int idx_i = 0; idx_i < b.nel(); idx_i++){
      if(!(b.el_sig() -> at(idx_i))){ continue; }
      comp = b.el_miniso() -> at(idx_i) < comp ? comp : b.el_miniso() -> at(idx_i);
    }
    for(int idx_i = 0; idx_i < b.nmu(); idx_i++){
      if(!(b.mu_sig() -> at(idx_i))){ continue; }
      comp = b.mu_miniso() -> at(idx_i) < comp ? comp : b.mu_miniso() -> at(idx_i);
    }

    return comp;
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

  double pTt_corr_d(const Baby &b){ 
    TVector3 H = AssignH(b).Vect();
    TLorentzVector Z = AssignZ(b);
    TLorentzVector gamma = AssignGamma(b);
    TVector3 ZmG = (Z-gamma).Vect();
    H.SetZ(0); ZmG.SetZ(0); 
    return H.Cross(ZmG.Unit()).Mag();
  }


  double cross_HZmG(const Baby &b){ 
    TVector3 H = AssignH(b).Vect();
    TLorentzVector Z = AssignZ(b);
    TLorentzVector gamma = AssignGamma(b);
    TVector3 ZmG = (Z-gamma).Vect();
    
    return (H.Cross(ZmG)).Mag()/(H.Mag());
  }

  double cross_HZmG_cms(const Baby &b){ 
    TVector3 H = AssignH(b).Vect();
    TLorentzVector Z = AssignZ(b);
    TLorentzVector gamma = AssignGamma(b);
    TVector3 ZmG = (Z-gamma).Vect();

    return (H.Cross(ZmG)).Pz()/(H.Mag());
  }



  double cross_HZmG_norm(const Baby &b){ 
    TVector3 H = AssignH(b).Vect();
    TLorentzVector Z = AssignZ(b);
    TLorentzVector gamma = AssignGamma(b);
    TVector3 ZmG = (Z-gamma).Vect();

    return (H.Cross(ZmG)).Mag()/(H.Mag()*ZmG.Mag());
  }

  double cross_HZmG_ZmGnorm(const Baby &b){ 
    TVector3 H = AssignH(b).Vect();
    TLorentzVector Z = AssignZ(b);
    TLorentzVector gamma = AssignGamma(b);
    TVector3 ZmG = (Z-gamma).Vect();

    return (H.Cross(ZmG)).Mag()/(ZmG.Mag());
  }

  double dot_HZmG(const Baby &b){ 
    TVector3 H = AssignH(b).Vect();
    TLorentzVector Z = AssignZ(b);
    TLorentzVector gamma = AssignGamma(b);
    TVector3 ZmG = (Z-gamma).Vect();

    return (H.Dot(ZmG))/(H.Mag());
  }

  int Nlep_truth(const Baby &b){
    int nlep = 0;
    for(unsigned int idx_mc = 0; idx_mc < b.mc_id() -> size(); idx_mc++){
      if(abs(b.mc_id() -> at(idx_mc)) == 11 || abs(b.mc_id() -> at(idx_mc)) == 13){ nlep++; } 
    }
    return nlep;
  }

  bool checkBit(int i, int n) {
    return((i%static_cast<int>(pow(2,n+1)))/static_cast<int>(pow(2,n)));
  }


  int Nb_truth(const Baby &b){
    int nb = 0; int bitmap;

    for(unsigned int idx_mc = 0; idx_mc < b.mc_id() -> size(); idx_mc++){
      if(abs(b.mc_id() -> at(idx_mc)) != 5){ continue;} 
      bitmap = b.mc_statusflag() ->at(idx_mc);
      if( checkBit(bitmap,12)){ continue; } 
      nb++;
    }

    return nb;
  }

  int Nel_truth(const Baby &b){
    int nb = 0; int bitmap;

    for(unsigned int idx_mc = 0; idx_mc < b.mc_id() -> size(); idx_mc++){
      if(abs(b.mc_id() -> at(idx_mc)) != 11){ continue;} 
      bitmap = b.mc_statusflag() ->at(idx_mc);
      if( checkBit(bitmap,12)){ continue; } 
      nb++;
    }

    return nb;
  }

  const NamedFunc N_SF_lpairs("N_SF_lpairs",[](const Baby &b) -> NamedFunc::ScalarType{return N_SF_leppairs(b);});
  const NamedFunc N_OSSF_lpairs("N_OSSF_lpairs",[](const Baby &b) -> NamedFunc::ScalarType{return N_OSSF_leppairs(b);});
  const NamedFunc Ntrub("Ntrub",[](const Baby &b) -> NamedFunc::ScalarType{return Nb_truth(b);});
  const NamedFunc Ntruel_NF("Ntruel_NF",[](const Baby &b) -> NamedFunc::ScalarType{return Nel_truth(b);});

  const NamedFunc mbb("mbb",[](const Baby &b) -> NamedFunc::ScalarType{return bb_cand(b,0).M();});
  const NamedFunc dr_bb_H("dr_bb_H",[](const Baby &b) -> NamedFunc::ScalarType{return bb_cand(b,0).DeltaR(AssignH(b));});
  const NamedFunc mtop("mtop",[](const Baby &b) -> NamedFunc::ScalarType{return m_top(b);});
  const NamedFunc mtop_l3("mtop_l3",[](const Baby &b) -> NamedFunc::ScalarType{return m_top_l3(b);});

  const NamedFunc pTt_corr("pTt_corr",[](const Baby &b) -> NamedFunc::ScalarType{return pTt_corr_d(b);});
  const NamedFunc HxZmG("HxZmG",[](const Baby &b) -> NamedFunc::ScalarType{return cross_HZmG(b);});
  const NamedFunc HxZmG_norm("HxZmG_norm",[](const Baby &b) -> NamedFunc::ScalarType{return cross_HZmG_norm(b);});
  const NamedFunc HxZmG_ZmGnorm("HxZmG_ZmGnorm",[](const Baby &b) -> NamedFunc::ScalarType{return cross_HZmG_ZmGnorm(b);});
  const NamedFunc HxZmG_cms("HxZmG_cms",[](const Baby &b) -> NamedFunc::ScalarType{return cross_HZmG_cms(b);});

  const NamedFunc HdZmG("HdZmG",[](const Baby &b) -> NamedFunc::ScalarType{return dot_HZmG(b);});
  const NamedFunc lly_pT("NF_lly_pT",[](const Baby &b) -> NamedFunc::ScalarType{return b.llphoton_pt()->at(0);});
  const NamedFunc lly_dR("lly_dR",[](const Baby &b) -> NamedFunc::ScalarType{return (AssignZ(b).DeltaR(AssignGamma(b)));});
  const NamedFunc H_Z_dR("H_Z_dR",[](const Baby &b) -> NamedFunc::ScalarType{return (AssignZ(b).DeltaR(AssignH(b)));});
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

  const NamedFunc max_Imini("max_Imini",[](const Baby &b) -> NamedFunc::ScalarType{return max_Imini_d(b);});
  const NamedFunc ptj_m4j("ptj_m4j",[](const Baby &b) -> NamedFunc::ScalarType{return ptj_m4j_d(b);});
  const NamedFunc ptj_m5j("ptj_m5j",[](const Baby &b) -> NamedFunc::ScalarType{return ptj_m5j_d(b);});
  const NamedFunc max_deepflav("max_deepflav",[](const Baby &b) -> NamedFunc::ScalarType{return max_deepflav_d(b);});
  const NamedFunc Z2_pt("Z2_pt",[](const Baby &b) -> NamedFunc::ScalarType{return Z2_pt_d(b);});
  const NamedFunc Z2_m("Z2_m",[](const Baby &b) -> NamedFunc::ScalarType{return Z2_m_d(b);});
  const NamedFunc Z2_flav("Z2_flav",[](const Baby &b) -> NamedFunc::ScalarType{return Z2_flav_b(b);});
  const NamedFunc Z2y_dR("Z2y_dR",[](const Baby &b) -> NamedFunc::ScalarType{return Z2y_dR_d(b);});
  const NamedFunc ZZ_dR("ZZ_dR",[](const Baby &b) -> NamedFunc::ScalarType{return ZZ_dR_d(b);});
  const NamedFunc dphi_Hmet("dphi_Hmet",[](const Baby &b) -> NamedFunc::ScalarType{return dphi_Hmet_d(b);});

  const NamedFunc Z1_m_avg("Z1_m_avg",[](const Baby &b) -> NamedFunc::ScalarType{    return Z1_avg(b).M();});
  const NamedFunc Z1_pt_avg("Z1_pt_avg",[](const Baby &b) -> NamedFunc::ScalarType{  return Z1_avg(b).Pt();});
  const NamedFunc Z1_eta_avg("Z1_eta_avg",[](const Baby &b) -> NamedFunc::ScalarType{return Z1_avg(b).Eta();});
  const NamedFunc Z1_phi_avg("Z1_phi_avg",[](const Baby &b) -> NamedFunc::ScalarType{return Z1_avg(b).Phi();});
  const NamedFunc Z2_m_avg("Z2_m_avg",[](const Baby &b) -> NamedFunc::ScalarType{    return Z2_avg(b).M();});
  const NamedFunc Z2_pt_avg("Z2_pt_avg",[](const Baby &b) -> NamedFunc::ScalarType{  return Z2_avg(b).Pt();});
  const NamedFunc Z2_eta_avg("Z2_eta_avg",[](const Baby &b) -> NamedFunc::ScalarType{return Z2_avg(b).Eta();});
  const NamedFunc Z2_phi_avg("Z2_phi_avg",[](const Baby &b) -> NamedFunc::ScalarType{return Z2_avg(b).Phi();});

  const NamedFunc Zg_m_avg("Zg_m_avg",[](const Baby &b) -> NamedFunc::ScalarType{    return (Z1_avg(b) + AssignGamma(b)).M();});
  const NamedFunc Zg_pt_avg("Zg_pt_avg",[](const Baby &b) -> NamedFunc::ScalarType{  return (Z1_avg(b) + AssignGamma(b)).Pt();});
  const NamedFunc Zg_eta_avg("Zg_eta_avg",[](const Baby &b) -> NamedFunc::ScalarType{return (Z1_avg(b) + AssignGamma(b)).Eta();});
  const NamedFunc Zg_phi_avg("Zg_phi_avg",[](const Baby &b) -> NamedFunc::ScalarType{return (Z1_avg(b) + AssignGamma(b)).Phi();});

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
  const std::vector<NamedFunc>   run2_categories = {run2_ggF, run2_VBF, run2_lep};
  const std::vector<std::string> run2_cat_labs   = {"_run2_ggF", "_run2_VBF", "_run2_lep"};

  const NamedFunc cat_ttH_lep = "(nlep==3 && njet>=3 && nbdfm>=1) || (nlep>=4 && nbdfm>=1)";
  const NamedFunc cat_ttH_had = "(nlep==2 && njet>=5 && nbdfm>=1)";
  const NamedFunc cat_ZH_met  = "(nlep==2 && njet<=1 && met>90)"; 
  const NamedFunc cat_WH_3l   = "(nlep>=3 && nbdfm==0)";
  const NamedFunc cat_VBF     = "(nlep==2 && njet>=2 && nbdfm==0)";
  const NamedFunc cat_ggF     = "(nlep==2 && njet<=1 && met<90)";

  const NamedFunc ZH_4l_eptcuts  = l3_pt > 20 && l4_pt > 20 && l3_flav==11 && l4_flav==11;
  const NamedFunc ZH_4l_muptcuts = l3_pt > 15 && l4_pt > 15 && l3_flav==13 && l4_flav==13;

  const NamedFunc baseline_ttH_lep =  max_Imini < 0.1 && "ll_m[0] > 85 && ll_m[0] < 95";
  const NamedFunc baseline_ttH_had = "ll_m[0] > 85 && ll_m[0] < 95";
  const NamedFunc baseline_ZH_met  = pTlly_mlly > 0.4 && HxZmG < 60 && "ll_m[0] > 85 && ll_m[0] < 95"; 
  const NamedFunc baseline_WH_3l   = max_Imini < 0.15 && pTlly_mlly > 0.3 && "ll_m[0] > 85 && ll_m[0] < 95" && "met > 30";
  const NamedFunc baseline_VBF     = "dijet_m > 600";

  const NamedFunc baseline_nomll_ttH_lep =  max_Imini < 0.1;
  const NamedFunc baseline_nomll_ttH_had = "1";
  const NamedFunc baseline_nomll_ZH_met  = pTlly_mlly > 0.4 && HxZmG < 60;
  const NamedFunc baseline_nomll_WH_3l   = max_Imini < 0.15 && pTlly_mlly > 0.3 && "met > 30";
  const NamedFunc baseline_nomll_VBF     = "dijet_m > 600";

  const std::vector<NamedFunc> run3_category_vector       = {cat_ttH_lep,cat_ttH_had,cat_ZH_met,cat_WH_3l,cat_VBF,cat_ggF};
  const std::vector<NamedFunc> run3_catwsel_vector        = {cat_ttH_lep && baseline_ttH_lep, cat_ttH_had && baseline_ttH_had,
                                                             cat_ZH_met && baseline_ZH_met, cat_WH_3l && baseline_WH_3l, cat_VBF, cat_ggF};
  const std::vector<NamedFunc> run3_catwsel_nomll_vector  = {cat_ttH_lep && baseline_nomll_ttH_lep, cat_ttH_had && baseline_nomll_ttH_had,
                                                             cat_ZH_met && baseline_nomll_ZH_met, cat_WH_3l && baseline_nomll_WH_3l, cat_VBF, cat_ggF}; 

  const std::vector<std::string> run3_category_labels = {"cat_ttH_lep","cat_ttH_had","cat_ZH_MET","cat_WH_3l","cat_VBF","cat_ggF"}; 


  //Some repeat functions not used outside of categorization utilities  
  const NamedFunc pTy_mlly("pTy_mlly",[](const Baby &b) -> NamedFunc::ScalarType{ return (b.photon_pt()->at(0))/ZgUtilities::Findlly(b); });
  const NamedFunc mlly("mlly",[](const Baby &b) -> NamedFunc::ScalarType{   return AssignH(b).M(); });
  const NamedFunc pT_lly("pT_lly",[](const Baby &b) -> NamedFunc::ScalarType{   return AssignH(b).Pt(); });
  const NamedFunc eta_lly("eta_lly",[](const Baby &b) -> NamedFunc::ScalarType{ return AssignH(b).Eta(); });
  const NamedFunc phi_lly("phi_lly",[](const Baby &b) -> NamedFunc::ScalarType{ return AssignH(b).Phi(); });
  const NamedFunc sigy_pTy("sigy_pTy",[](const Baby &b) -> NamedFunc::ScalarType{ return (b.photon_energyErr() -> at(0))/(b.photon_pt()->at(0)); });

  //Photon plots added to each category's plotting function
  void sample_photon_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    pm.Push<Hist1D>(Axis(45,0.0,1,      "photon_idmva[0]", "#gamma - IDMVA",            {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_IDMVA_");
    pm.Push<Hist1D>(Axis(35,0,3.5,      "photon_drmin[0]", "min( #Delta R(#gamma,l) )", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_minDR_");
    pm.Push<Hist1D>(Axis(35,0,3.5,      "photon_drmax[0]", "max( #Delta R(#gamma,l) )", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_maxDR_");
    pm.Push<Hist1D>(Axis(40,10,60,      "photon_pt[0]",    "p_{T}(#gamma) [GeV]",       {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_pT_");
    pm.Push<Hist1D>(Axis(40,-2.6,2.6,   "photon_eta[0]",   "#eta(#gamma)",              {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_eta_");
    pm.Push<Hist1D>(Axis(63,-3.15,3.15, "photon_phi[0]",   "#phi(#gamma)",              {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_phi_");
    pm.Push<Hist1D>(Axis(20,0.01,0.25,  "photon_energyErr[0]/photon_pt[0]", "#sigma(#gamma)/E(#gamma)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_relunc_");
  }


  //ll and llphoton plots added to each category's photon plots
  void sample_llphoton_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    pm.Push<Hist1D>(Axis(40,100,180,    mlly,          "m_{ll#gamma} [GeV]",     {122,128}), selection, processes ,ops).Weight(wgt).Tag("ShortName:" + labels + "_lly_m_");
    pm.Push<Hist1D>(Axis(50,0,250,      pT_lly,        "p_{T}(ll#gamma) [GeV]",         {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_lly_pt_");
    pm.Push<Hist1D>(Axis(52,-2.6,2.6,   eta_lly,       "#eta(ll#gamma)",                {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_lly_eta_");
    pm.Push<Hist1D>(Axis(63,-3.15,3.15, phi_lly,       "#phi(ll#gamma)",                {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_lly_phi_");
    pm.Push<Hist1D>(Axis(40,0,2.5,      pTlly_mlly,    "p_{T}(ll#gamma)/m_{ll#gamma}", {1}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_lly_ptom_");
    pm.Push<Hist1D>(Axis(40,0,1.2,      pTy_mlly,      "p_{T}(#gamma)/m_{ll#gamma}",    {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_y_ptom_");
    pm.Push<Hist1D>(Axis(45,0,4.5,      lly_dR,        "#Delta R(ll,#gamma)",           {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ll_y_dR_");
    //pm.Push<Hist1D>(Axis(50,0,100,      "llphoton_pTt","p^{t}_{T}",                     {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_pTt_");

    pm.Push<Hist1D>(Axis(50,0,100,      "ll_pt[0]",  "p_{T}(ll) [GeV]",{}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ll_pt_");
    pm.Push<Hist1D>(Axis(52,-2.6,2.6,   "ll_eta[0]", "#eta(ll)",       {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ll_eta_");
    pm.Push<Hist1D>(Axis(63,-3.15,3.15, "ll_phi[0]", "#phi(ll)",       {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ll_phi_");
    pm.Push<Hist1D>(Axis(50,80,100,     "ll_m[0]",   "m_{ll} [GeV]",   {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ll_m_");

  }

  //Basic lepton kinematics plots
  void sample_lepton_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    pm.Push<Hist1D>(Axis(40,5,65,     "el_pt[ll_i1[0]]",  "p_{T}(e_{1})",   {}), selection && "ll_lepid[0]==11", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_10_pt_el1_");
    pm.Push<Hist1D>(Axis(40,5,65,     "el_pt[ll_i2[0]]",  "p_{T}(e_{2})",   {}), selection && "ll_lepid[0]==11", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_10_pt_el2_");
    pm.Push<Hist1D>(Axis(40,-2.6,2.6, "el_eta[ll_i1[0]]", "#eta(e_{1})",    {}), selection && "ll_lepid[0]==11", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_11_eta_el1_");
    pm.Push<Hist1D>(Axis(40,-2.6,2.6, "el_eta[ll_i2[0]]", "#eta(e_{2})",    {}), selection && "ll_lepid[0]==11", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_11_eta_el2_");
    pm.Push<Hist1D>(Axis(40,0,75,     "mu_pt[ll_i1[0]]",  "p_{T}(#mu_{1})", {}), selection && "ll_lepid[0]==13", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_10_pt_mu1_");
    pm.Push<Hist1D>(Axis(40,0,75,     "mu_pt[ll_i2[0]]",  "p_{T}(#mu_{2})", {}), selection && "ll_lepid[0]==13", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_10_pt_mu2_");
    pm.Push<Hist1D>(Axis(40,-2.6,2.6, "mu_eta[ll_i1[0]]", "#eta(#mu_{1})",  {}), selection && "ll_lepid[0]==13", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_11_eta_mu1_");
    pm.Push<Hist1D>(Axis(40,-2.6,2.6, "mu_eta[ll_i2[0]]", "#eta(#mu_{2})",  {}), selection && "ll_lepid[0]==13", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_11_eta_mu2_");
  }

  //Sample ttH leptonic category plots
  void ttH_lep_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    sample_photon_plots(  pm, selection, processes, ops, wgt, labels + "_ttH_lep_");
    sample_lepton_plots(  pm, selection, processes, ops, wgt, labels + "_ttH_lep_");
    sample_llphoton_plots(pm, selection, processes, ops, wgt, labels + "_ttH_lep_");

    pm.Push<Hist1D>(Axis(80,0,500, "ht",    "H_{T} [GeV]", {120}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ttH_lep_HT_"); 
    pm.Push<Hist1D>(Axis(75,0,200, "met",   "MET [GeV]",      {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ttH_lep_MET_");
    pm.Push<Hist1D>(Axis(10,0,10,  "njet",  "N_{jet}",       {3}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ttH_lep_Njet_");
    pm.Push<Hist1D>(Axis(6,0,6,    "nlep",  "N_{lep}",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ttH_lep_Nlep_");
    pm.Push<Hist1D>(Axis(6,0,6 ,   "nbdfm", "N_{b,med}",      {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ttH_lep_Nbdfm");
    pm.Push<Hist1D>(Axis(6,0,6 ,   "ntruel","N_{e,truth}",    {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ttH_lep_Ntruel");
    pm.Push<Hist1D>(Axis(6,0,6 ,   "ntrumu","N_{#mu,truth}",  {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ttH_lep_Ntrumu");
    pm.Push<Hist1D>(Axis(6,0,6 ,   Ntrub,   "N_{b,truth}",    {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ttH_lep_Ntrub");
    pm.Push<Hist1D>(Axis(6,0,6 ,   Ntruel_NF,"N_{e,truth}",   {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ttH_lep_Ntruel_NF");

    pm.Push<Hist1D>(Axis(50,100,300, mtop_l3,   "cand. m_{t} [GeV]", {130,220}), selection && "nlep==3",   processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ttH_lep_mtop");
    pm.Push<Hist1D>(Axis(50,0,1,     l3_mini,   "I_{mini} - l3",            {}), selection,                processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ttH_lep_l3_Imini");
    pm.Push<Hist1D>(Axis(55,-0.1,1,  l3_idmva,  "l_{3} - e IDMVA ",         {}), selection && l3_flav==11, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ttH_lep_e3idmva"); 
    pm.Push<Hist1D>(Axis(35,5,75,    l3_pt,     "p_{T}(l_{3}) [GeV]",       {}), selection && l3_flav==11, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ttH_lep_e3_pt");
    pm.Push<Hist1D>(Axis(35,5,75,    l3_pt,     "p_{T}(l_{3}) [GeV]",       {}), selection && l3_flav==13, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ttH_lep_mu3_pt");
    pm.Push<Hist1D>(Axis(45,0,4.5,   H_Z_dR,    "#Delta R(ll,#gamma)",      {}), selection,                processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ttH_lep_H_Z_dR_");
    pm.Push<Hist1D>(Axis(50,0,1,     max_Imini, "max. I_{mini}",            {}), selection,                processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ttH_lep_max_Imini");

  }

  //Sample ttH had category plots
  void ttH_had_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    sample_photon_plots(  pm, selection, processes, ops, wgt, labels + "_ttH_had_");
    sample_lepton_plots(  pm, selection, processes, ops, wgt, labels + "_ttH_had_");
    sample_llphoton_plots(pm, selection, processes, ops, wgt, labels + "_ttH_had_");

    pm.Push<Hist1D>(Axis(10,0,10,    "njet",       "N_{jet}",                             {4}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ttH_had_Njet_");
    pm.Push<Hist1D>(Axis(75,0,200,   "met",        "MET",                                  {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ttH_had_MET");
    pm.Push<Hist1D>(Axis(80,150,500, "ht",         "H_{T}",                             {150}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ttH_had_HT_");
    pm.Push<Hist1D>(Axis(50,100,300, mtop,         "cand. m_{t} [GeV]",             {130,220}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ttH_had_mtop");
    pm.Push<Hist1D>(Axis(50,0,1,     ptj_m4j,      "#Sigma(p_{T}(j_{i}))/m_{4j} [GeV]", {0.6}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ttH_had_ptj_m4j");
    pm.Push<Hist1D>(Axis(50,0,1,     ptj_m5j,      "#Sigma(p_{T}(j_{i}))/m_{5j} [GeV]", {0.6}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ttH_had_ptj_m5j");
    pm.Push<Hist1D>(Axis(40,0,1,     max_deepflav, "Deep Flavor - max",                    {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ttH_had_deepflav_BASE_MAX");
    pm.Push<Hist1D>(Axis(50,0,1,     max_Imini,    "max. I_{mini}",                        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ttH_had_max_Imini");

  }


  //Sample VH 4l category plots
  void ZH_4l_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    sample_photon_plots(  pm, selection, processes, ops, wgt, labels + "_ZH_4l_");
    sample_lepton_plots(  pm, selection, processes, ops, wgt, labels + "_ZH_4l_");
    sample_llphoton_plots(pm, selection, processes, ops, wgt, labels + "_ZH_4l_");

    pm.Push<Hist1D>(Axis(50,0,1,    l3_mini,  "I_{mini} - l3", {}),           selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_4l_Imini_3l");
    pm.Push<Hist1D>(Axis(35,5,75,   l3_pt,    "p_{T}(e_{3}) [GeV]", {}),      selection && l3_flav==11, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_4l_e3_pt");
    pm.Push<Hist1D>(Axis(35,5,75,   l3_pt,    "p_{T}(#mu_{3}) [GeV]", {}),    selection && l3_flav==13, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_4l_mu3_pt");
    pm.Push<Hist1D>(Axis(55,-0.1,1, l3_idmva, "l_{3} - e IDMVA ", {}),        selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_4l_e3idmva");

    pm.Push<Hist1D>(Axis(50,0,1,    l4_mini, "I_{mini} - l4", {}),            selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_4l_Imini_4l");
    pm.Push<Hist1D>(Axis(35,5,75,   l4_pt,   "p_{T}(e_{4}) [GeV]", {}),       selection && l4_flav==11, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_4l_e4_pt");
    pm.Push<Hist1D>(Axis(35,5,75,   l4_pt,   "p_{T}(#mu_{4}) [GeV]", {}),     selection && l4_flav==13, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_4l_mu4_pt");
    pm.Push<Hist1D>(Axis(55,-0.1,1, l4_idmva, "l_{4} - e IDMVA ", {}),        selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_4l_e4idmva"); 

    pm.Push<Hist1D>(Axis(80,0,400,  Z2_pt,  "p_{T}(ll) - Z_{2} [GeV]", {50}),           selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_4l_ptll_Z2");
    pm.Push<Hist1D>(Axis(35,50,120, Z2_m,   "m_{l_{3} l_{4}} - Z_{2} [GeV]", {80,100}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_4l_mll_Z2"); 
    pm.Push<Hist1D>(Axis(50,0.0,5,  Z2y_dR, "#Delta R(Z_{2},#gamma) [GeV]", {0.2}),     selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_4l_Z2y_dr");
    pm.Push<Hist1D>(Axis(50,0.0,5,  ZZ_dR,  "#Delta R(Z_{1},Z_{2}) [GeV]", {}),         selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_4l_ZZ_dr");
    pm.Push<Hist1D>(Axis(50,0,1,     max_Imini, "max. I_{mini}", {}),         selection,                processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_4l_max_Imini");

    pm.Push<Hist1D>(Axis(35,50,120,   Z1_m_avg,   "m_{Z_{1}} [GeV]", {80,100}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_4l_mll_Z1_AVG"); 
    pm.Push<Hist1D>(Axis(80,0,400,    Z1_pt_avg,  "p_{T}(Z_{1}) [GeV]", {50}),  selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_4l_ptll_Z1_AVG");
    pm.Push<Hist1D>(Axis(80,-4,4,     Z1_eta_avg,  "#eta(Z_{1}) [GeV]", {}),    selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_4l_etall_Z1_AVG");
    pm.Push<Hist1D>(Axis(6.4,-3.2,3.2,Z1_phi_avg,  "#phi(Z_{1}) [GeV]", {}),    selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_4l_phill_Z1_AVG");

    pm.Push<Hist1D>(Axis(35,50,120,   Z2_m_avg,   "m_{Z_{2}} [GeV]", {80,100}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_4l_mll_Z2_AVG"); 
    pm.Push<Hist1D>(Axis(80,0,400,    Z2_pt_avg,  "p_{T}(Z_{2}) [GeV]", {50}),  selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_4l_ptll_Z2_AVG");
    pm.Push<Hist1D>(Axis(80,-4,4,     Z2_eta_avg,  "#eta(Z_{2}) [GeV]", {}),    selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_4l_etall_Z2_AVG");
    pm.Push<Hist1D>(Axis(6.4,-3.2,3.2,Z2_phi_avg,  "#phi(Z_{2}) [GeV]", {}),    selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_4l_phill_Z2_AVG");

    pm.Push<Hist1D>(Axis(35,100,180,  Zg_m_avg,   "m_{Z#gamma} [GeV]", {80,100}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_4l_mll_Zg_AVG"); 
    pm.Push<Hist1D>(Axis(80,0,400,    Zg_pt_avg,  "p_{T}(Z#gamma) [GeV]", {50}),  selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_4l_ptll_Zg_AVG");
    pm.Push<Hist1D>(Axis(80,-4,4,     Zg_eta_avg,  "#eta(Z#gamma) [GeV]", {}),    selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_4l_etall_Zg_AVG");
    pm.Push<Hist1D>(Axis(6.4,-3.2,3.2,Zg_phi_avg,  "#phi(Z#gamma) [GeV]", {}),    selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_4l_phill_Zg_AVG");


  }

  //Sample WH 3l category plots
  void WH_3l_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    sample_photon_plots(pm,   selection, processes, ops, wgt, labels + "_WH_3l_");
    sample_lepton_plots(pm,   selection, processes, ops, wgt, labels + "_WH_3l_");
    sample_llphoton_plots(pm, selection, processes, ops, wgt, labels + "_WH_3l_");

    pm.Push<Hist1D>(Axis(75,0,200,  "met",     "MET",                   {}), selection,                processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _WH_3l_MET");
    pm.Push<Hist1D>(Axis(40,0,3.2,  dphi_Hmet, "#Delta #phi(MET,H)", {2.0}), selection,                processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _WH_3l_dphi_Hmet");
    pm.Push<Hist1D>(Axis(40,0,4.0,  l3_Z_DR,   "#Delta R(Z,l_{3})",   {45}), selection,                processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _WH_3l_Zl3_DR");
    pm.Push<Hist1D>(Axis(50,0,1,    l3_mini,   "I_{mini} - l3 [GeV]",   {}), selection,                processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _WH_3l_lep_l3_Imini");
    pm.Push<Hist1D>(Axis(55,-0.1,1, l3_idmva,  "l_{3} - e IDMVA ",      {}), selection,                processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _WH_3l_e3idmva");
    pm.Push<Hist1D>(Axis(50,0,1,    max_Imini, "max. I_{mini}",         {}), selection,                processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _WH_3l_max_Imini");
    pm.Push<Hist1D>(Axis(35,5,75,   l3_pt,     "p_{T}(e_{3})",          {}), selection && l3_flav==11, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _WH_3l_e3_pt");
    pm.Push<Hist1D>(Axis(35,5,75,   l3_pt,     "p_{T}(#mu_{3})",        {}), selection && l3_flav==13, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _WH_3l_mu3_pt");

  }

  //Sample ZH 2bl category plots
  void ZH_2bl_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    sample_photon_plots(  pm, selection, processes, ops, wgt, labels + "_ZH_2bl");
    sample_lepton_plots(  pm, selection, processes, ops, wgt, labels + "_ZH_2bl");
    sample_llphoton_plots(pm, selection, processes, ops, wgt, labels + "_ZH_2bl");

    pm.Push<Hist1D>(Axis(10,0,10,  "njet",        "N_{jet}", {2}),            selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_2bl_Njet_");
    pm.Push<Hist1D>(Axis(10,0,10,  "nbdfm",       "N_{bdfm}", {2}),           selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_2bl_Nbdfm_"); 

    pm.Push<Hist1D>(Axis(40,0,4.0,  dr_bb_H,      "#Delta R(Z_{bb},H_{ll#gamma})", {45}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_2bl_Zbb_H_DR");
    pm.Push<Hist1D>(Axis(40,0,1,    max_deepflav, "Max. Deep Flavor", {0.4}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_2bl_deepflav_BASE_MAX");
    pm.Push<Hist1D>(Axis(50,50,150, mbb,          "m_{bb}", {}),              selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_2bl_mbb");
    pm.Push<Hist1D>(Axis(75,0,200, "met",         "p_{T}^{miss}", {}),        selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_2bl_MET");
    pm.Push<Hist1D>(Axis(80,0,500, "ht",          "H_{T}", {}),               selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_2bl_HT_");
    pm.Push<Hist1D>(Axis(50,0,1,    max_Imini, "max. I_{mini}", {}),          selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_2bl_max_Imini");

  }

  //Sample ZH MET category plots
  void ZH_MET_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    sample_photon_plots(  pm, selection, processes, ops, wgt, labels + "_ZH_MET_");
    sample_lepton_plots(  pm, selection, processes, ops, wgt, labels + "_ZH_MET_");
    sample_llphoton_plots(pm, selection, processes, ops, wgt, labels + "_ZH_MET_");

    pm.Push<Hist1D>(Axis(6,0,6,     "njet",    "N_{jet}", {}),               selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_MET_njet_");
    pm.Push<Hist1D>(Axis(6,0,6,     "nbdfm",   "N_{b,med}", {}),             selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_MET_nbdfm_");
    pm.Push<Hist1D>(Axis(6,0,6,     "nbdfl",   "N_{b,loose}", {}),           selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_MET_nbdfl_");

    pm.Push<Hist1D>(Axis(55,90,200, "met",     "p_{T}^{miss}", {}),          selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_MET_MET_");
    pm.Push<Hist1D>(Axis(40,0,3.2,  dphi_Hmet, "#Delta #phi(MET,H)", {2.5}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_MET_dphi_Hmet");
    pm.Push<Hist1D>(Axis(50,0,1,    max_Imini, "max. I_{mini}", {}),         selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + " _ZH_MET_max_Imini");

  }

  //These are the control region plots made for the ggF category
  //Note mllg is not llphoton_m[0] because often control regions include photons that don't pass the photon object selections
  void ggF_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    sample_photon_plots(  pm, selection, processes, ops, wgt, labels + "_ggF_");
    sample_lepton_plots(  pm, selection, processes, ops, wgt, labels + "_ggF_");
    sample_llphoton_plots(pm, selection, processes, ops, wgt, labels + "_ggF_");

    //--kinematic angles--//
    pm.Push<Hist1D>(Axis(40,-1,1,      "llphoton_cosTheta", "cos(#Theta)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_cTHETA_");
    pm.Push<Hist1D>(Axis(40,-1,1,      "llphoton_costheta", "cos(#theta)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ctheta_");
    pm.Push<Hist1D>(Axis(40,-3.2,3.2,  "llphoton_psi",      "#phi",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_phi_");
  }

  //This function again adds control region plots for the VBF category based on variables used in Run 2. Some overlap with ggF region but removing variables specific to the "kinematic BDT"
  void VBF_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    sample_photon_plots(  pm, selection, processes, ops, wgt, labels + "_VBF_");
    sample_lepton_plots(  pm, selection, processes, ops, wgt, labels + "_VBF_");
    sample_llphoton_plots(pm, selection, processes, ops, wgt, labels + "_VBF_");

    pm.Push<Hist1D>(Axis(80,0,1200,  "dijet_m",    "m_{jj} [GeV]",            {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_dijet_m");
    pm.Push<Hist1D>(Axis(12,0,9,     "dijet_deta", "#Delta#eta(j_{1},j_{2})", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_deta");
    pm.Push<Hist1D>(Axis(10,0,3,     "dijet_dphi", "#Delta#phi(j_{1},j_{2})", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_dphi");

    pm.Push<Hist1D>(Axis(30, 0, 3,   "llphoton_dijet_dphi[0]" ,"#Delta#phi(Z#gamma,jj)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_jjlly_dphi");
    pm.Push<Hist1D>(Axis(30, 0, 1.0, "llphoton_dijet_balance[0]" ,"System balance",         {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_sysbal");
    //pm.Push<Hist1D>(Axis(20, 0, 140, "llphoton_pTt[0]"       ,"p^{t}_{T} [GeV]",           {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_pTt2");

    pm.Push<Hist1D>(Axis( 8,   0,  6, "photon_zeppenfeld[0]", "#gamma zeppenfeld",      {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_zep");
    pm.Push<Hist1D>(Axis(15, 0.4,4.4, "photon_jet_mindr[0]",  "Min. #DeltaR(#gamma,j)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_yjdr");

    pm.Push<Hist1D>(Axis(90,    30,  300, "jet_pt[0]",  "p_{T}(j_{1}) [GeV]",{}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_j1_pt");
    pm.Push<Hist1D>(Axis(90,    30,  300, "jet_pt[1]",  "p_{T}(j_{2}) [GeV]",{}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_j2_pt");
    pm.Push<Hist1D>(Axis(80,  -4.7,  4.7, "jet_eta[0]", "#eta(j_{1})",       {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_j1_eta_");
    pm.Push<Hist1D>(Axis(80,  -4.7,  4.7, "jet_eta[1]", "#eta(j_{2})",       {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_j2_eta_");
    pm.Push<Hist1D>(Axis(63, -3.15, 3.15, "jet_phi[0]", "#phi(j_{1})",       {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_j1_phi_");
    pm.Push<Hist1D>(Axis(63, -3.15, 3.15, "jet_phi[1]", "#phi(j_{2})",       {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_j2_phi_");
  }

  //This function, as of 11/03/23 is a work in progress. No control region plots were made in HIG-19-014
  void Lep_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    sample_photon_plots(  pm, selection, processes, ops, wgt, labels + "_lep_");
    sample_lepton_plots(  pm, selection, processes, ops, wgt, labels + "_lep_");
    sample_llphoton_plots(pm, selection, processes, ops, wgt, labels + "_lep_");
  }


}


