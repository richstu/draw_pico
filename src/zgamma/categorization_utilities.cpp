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
      {"2018", std::vector<float>({0.0494, 0.2770, 0.7264})},
      {"2022", std::vector<float>({0.047,  0.245,  0.6734})},
      {"2022EE", std::vector<float>({0.0499, 0.2605, 0.6915})},  
      {"2023", std::vector<float>({0.0358, 0.1917, 0.6172})},
      {"2023BPix", std::vector<float>({0.0359, 0.1919, 0.6133})}

    };
  
    if(deep_flav < btag_df_wpts[year][lmt]){return false;}
    return true;
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

  //This returns the flavor and idces of the leptons that make up the ll not used for the lly candidate
  std::vector<int> Z2_idces(const Baby &b){
    int idx_keep = Z2_idx(b);
    if(idx_keep==-1){ return {-1}; }
    return {b.ll_lepid() -> at(idx_keep), b.ll_i1() -> at(idx_keep),b.ll_i2() -> at(idx_keep)};
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

  //This function returns the mass of the "best" top candidate
  //First the two most likely b-jets are found by finding max and 2nd max b-jet score
  //Then the best W candidate is found via dijet mass
  //The two b-jets are compared to see which one has the most compatible mass with the top
  double m_top_r3(const Baby &b){
    std::vector <int> b_idx(2);
    std::vector <double> b_scores = {-1,-1};

    TLorentzVector temp1,temp2,W_cand,b_max,b_2max;
    double compDouble = 0;
    for(unsigned int idx_k = 0; idx_k < b.jet_pt() -> size(); idx_k++){
      if( !(b.jet_isgood()->at(idx_k)) ){continue;}
      compDouble = b.jet_btagpnetb()->at(idx_k);

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



  double llgamma_dijet_dphi_d(const Baby &b){
    //commented out now but will be available in later picos
    if(b.nllphoton() > 0){
      return b.llphoton_dijet_dphi() -> at(0);
    }
    TLorentzVector v1, v2; v1.SetPtEtaPhiM(b.photon_pt() ->at(0), b.photon_eta() ->at(0), b.photon_phi() ->at(0), 0);
    v2.SetPtEtaPhiM(b.ll_pt() ->at(0), b.ll_eta() ->at(0), b.ll_phi() ->at(0), b.ll_m() -> at(0));
    v1=v1+v2;
    v2.SetPtEtaPhiM(b.dijet_pt(), b.dijet_eta(), b.dijet_phi(), b.dijet_m());
     
    return v1.DeltaPhi(v2);
  }


  double llgamma_dijet_deltaR_d(const Baby &b){
    //commented out now but will be available in later picos
    if(b.nllphoton() > 0){
      return b.llphoton_dijet_dr() -> at(0);
    }
    TLorentzVector v1, v2; v1.SetPtEtaPhiM(b.photon_pt() ->at(0), b.photon_eta() ->at(0), b.photon_phi() ->at(0), 0);
    v2.SetPtEtaPhiM(b.ll_pt() ->at(0), b.ll_eta() ->at(0), b.ll_phi() ->at(0), b.ll_m() -> at(0));
    v1=v1+v2;
    v2.SetPtEtaPhiM(b.dijet_pt(), b.dijet_eta(), b.dijet_phi(), b.dijet_m());
     
    return v1.DeltaR(v2);
  }


  double llgamma_dijet_balance_d(const Baby &b){
    //commented out now but will be available in later picos
    if(b.nllphoton() > 0){
      return b.llphoton_dijet_balance() -> at(0);
    }
    return -1;
  }

  double photon_zeppenfeld_d(const Baby &b){
    //commented out now but will be available in later picos
    if(b.nllphoton() > 0){
      return b.photon_zeppenfeld() -> at(0);
    }
    return -1;
  }

  double photon_jet_mindr_d(const Baby &b){
    //commented out now but will be available in later picos
    if(b.nllphoton() > 0){
      return b.photon_jet_mindr() -> at(0);
    }
    double min_dr = 1000;
    double dr = 1000;
    TLorentzVector ph, jet; 
    ph.SetPtEtaPhiM(b.photon_pt() ->at(0), b.photon_eta() ->at(0), b.photon_phi() ->at(0), 0);
    for(size_t idx_jet = 0;idx_jet < b.jet_pt() -> size(); idx_jet++){
      if(!(b.jet_isgood()->at(idx_jet))){continue;}
      jet.SetPtEtaPhiM(b.jet_pt() ->at(idx_jet), b.jet_eta() ->at(idx_jet), b.jet_phi() ->at(idx_jet), b.jet_m()->at(idx_jet));
      dr = ph.DeltaR(jet);
      min_dr = dr < min_dr ? dr : min_dr;
    }
     
    return min_dr;  

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
      if(l12[0]!=-1){ //return l3;}
        l3 = AssignLep(b,l12[0],l12[1]);
        return l3;
      }
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
      if(l12[0]!=-1){
        return {l12[0],l12[1]};
      }
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

  //Returns the maximum deepflavor score of all jets in the event
  double max_particlenet_d(const Baby &b){
    float mx = -1;
    for(auto dfv : *b.jet_btagpnetb() ){
      mx = std::max(mx,dfv); 
    }
    return mx;
  }

  double l3_pt_d(const Baby &b){
    TLorentzVector l3;
    l3 = AssignL3(b);
    return l3.Pt();
  }

  double l3_eta_d(const Baby &b){
    TLorentzVector l3;
    l3 = AssignL3(b);
    return l3.Eta();
  }

  double l3_phi_d(const Baby &b){
    TLorentzVector l3;
    l3 = AssignL3(b);
    return l3.Phi();
  }

  double l3_idmva_d(const Baby &b){
    std::vector<int> l3_i = L3_Idx(b);
    if(l3_i[0] == 13 || l3_i[1]==-1){ return 2; }
    return b.el_idmva() -> at(l3_i[1]);
  }

  //Returns Delta R separation between the Z and the third lepton
  double l3_Z_DR_d(const Baby &b){
    TLorentzVector Z, l3;
    Z  = AssignZ(b); 
    l3 = AssignL3(b);
    return Z.DeltaR(l3);
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

  //Returns the pT of the fourth lepton
  double l4_pt_d(const Baby &b){
    TLorentzVector l4;
    l4 = AssignL4(b);
    return l4.Pt();
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


  double lly_ptom_d(const Baby &b){ 
    TLorentzVector H = AssignH(b);
    return H.Pt()/(H.M());
  }

  double lly_refit_ptom_d(const Baby &b){ 
    TLorentzVector H; H.SetPtEtaPhiM(b.llphoton_refit_pt(), b.llphoton_refit_eta(), b.llphoton_refit_phi(), b.llphoton_refit_m());
    return H.Pt()/(H.M());
  }


  size_t idx_j1(const Baby &b){
    for(size_t idx_jet = 0; idx_jet < b.jet_isgood() -> size(); idx_jet++){
      if(b.jet_isgood() ->at(idx_jet)){ return idx_jet; }
    }
    return -1;
  }

  size_t idx_j2(const Baby &b){
    bool first = false;
    for(size_t idx_jet = 0; idx_jet < b.jet_isgood() -> size(); idx_jet++){
      if(!first && b.jet_isgood() ->at(idx_jet)){first=true; continue;}
      if(first && b.jet_isgood() ->at(idx_jet)){ return idx_jet; }
    }
    return -1;
  }

  size_t idx_jN(const Baby &b, int num_jet){
    int count_jet = 0;
    for(size_t idx_jet = 0; idx_jet < b.jet_isgood() -> size(); idx_jet++){
      if( b.jet_isgood() ->at(idx_jet)){ count_jet++; }
      if(count_jet >= num_jet){ return idx_jet; }
    }
    return -1;
  }

  const NamedFunc j1_pt("j1_pt",[](const Baby &b) -> NamedFunc::ScalarType{return b.jet_pt() -> at(idx_j1(b));});
  const NamedFunc j2_pt("j2_pt",[](const Baby &b) -> NamedFunc::ScalarType{return b.jet_pt() -> at(idx_j2(b));});
  const NamedFunc j3_pt("j3_pt",[](const Baby &b) -> NamedFunc::ScalarType{return b.jet_pt() -> at(idx_jN(b,3));});
  const NamedFunc j4_pt("j4_pt",[](const Baby &b) -> NamedFunc::ScalarType{return b.jet_pt() -> at(idx_jN(b,4));});
  const NamedFunc j5_pt("j5_pt",[](const Baby &b) -> NamedFunc::ScalarType{return b.jet_pt() -> at(idx_jN(b,5));});

  const NamedFunc j1_eta("j1_eta",[](const Baby &b) -> NamedFunc::ScalarType{return b.jet_eta() -> at(idx_j1(b));});
  const NamedFunc j2_eta("j2_eta",[](const Baby &b) -> NamedFunc::ScalarType{return b.jet_eta() -> at(idx_j2(b));});
  const NamedFunc j3_eta("j3_eta",[](const Baby &b) -> NamedFunc::ScalarType{return b.jet_eta() -> at(idx_jN(b,3));});
  const NamedFunc j4_eta("j4_eta",[](const Baby &b) -> NamedFunc::ScalarType{return b.jet_eta() -> at(idx_jN(b,4));});
  const NamedFunc j5_eta("j5_eta",[](const Baby &b) -> NamedFunc::ScalarType{return b.jet_eta() -> at(idx_jN(b,5));});

  const NamedFunc j1_phi("j1_phi",[](const Baby &b) -> NamedFunc::ScalarType{return b.jet_phi() -> at(idx_j1(b));});
  const NamedFunc j2_phi("j2_phi",[](const Baby &b) -> NamedFunc::ScalarType{return b.jet_phi() -> at(idx_j2(b));});
  const NamedFunc j3_phi("j3_phi",[](const Baby &b) -> NamedFunc::ScalarType{return b.jet_phi() -> at(idx_jN(b,3));});
  const NamedFunc j4_phi("j4_phi",[](const Baby &b) -> NamedFunc::ScalarType{return b.jet_phi() -> at(idx_jN(b,4));});
  const NamedFunc j5_phi("j5_phi",[](const Baby &b) -> NamedFunc::ScalarType{return b.jet_phi() -> at(idx_jN(b,5));});

  const NamedFunc j1_m("j1_m",[](const Baby &b) -> NamedFunc::ScalarType{return b.jet_m() -> at(idx_j1(b));});
  const NamedFunc j2_m("j2_m",[](const Baby &b) -> NamedFunc::ScalarType{return b.jet_m() -> at(idx_j2(b));});
  const NamedFunc j3_m("j3_m",[](const Baby &b) -> NamedFunc::ScalarType{return b.jet_m() -> at(idx_jN(b,3));});
  const NamedFunc j4_m("j4_m",[](const Baby &b) -> NamedFunc::ScalarType{return b.jet_m() -> at(idx_jN(b,4));});
  const NamedFunc j5_m("j5_m",[](const Baby &b) -> NamedFunc::ScalarType{return b.jet_m() -> at(idx_jN(b,5));});

  const NamedFunc lly_pT("NF_lly_pT",[](const Baby &b) -> NamedFunc::ScalarType{return b.llphoton_pt()->at(0);});
  const NamedFunc lly_dR("lly_dR",[](const Baby &b) -> NamedFunc::ScalarType{return (AssignZ(b).DeltaR(AssignGamma(b)));});
  const NamedFunc pTlly_mlly("pTlly_mlly",[](const Baby &b) -> NamedFunc::ScalarType{return lly_ptom_d(b);});
  const NamedFunc pTlly_mlly_refit("pTlly_mlly_refit",[](const Baby &b) -> NamedFunc::ScalarType{return lly_refit_ptom_d(b);});

  const NamedFunc l3_pt("l3_pt",[](const Baby &b) -> NamedFunc::ScalarType{return l3_pt_d(b);});
  const NamedFunc l3_eta("l3_eta",[](const Baby &b) -> NamedFunc::ScalarType{return l3_eta_d(b);});
  const NamedFunc l3_phi("l3_phi",[](const Baby &b) -> NamedFunc::ScalarType{return l3_phi_d(b);});
  const NamedFunc l3_idmva("l3_idmva",[](const Baby &b) -> NamedFunc::ScalarType{return l3_idmva_d(b);});
  const NamedFunc l3_mini("l3_mini",[](const Baby &b) -> NamedFunc::ScalarType{return l3_Imini_d(b);});
  const NamedFunc l3_flav("l3_flav",[](const Baby &b) -> NamedFunc::ScalarType{return L3_flav_i(b);});
  const NamedFunc l3_Z_DR("l3_Z_DR",[](const Baby &b) -> NamedFunc::ScalarType{return l3_Z_DR_d(b);});
  const NamedFunc H_Z_dR("H_Z_dR",[](const Baby &b) -> NamedFunc::ScalarType{return (AssignZ(b).DeltaR(AssignH(b)));});

  const NamedFunc l4_pt("l4_pt",[](const Baby &b) -> NamedFunc::ScalarType{return l4_pt_d(b);});
  const NamedFunc l4_idmva("l4_idmva",[](const Baby &b) -> NamedFunc::ScalarType{return l4_idmva_d(b);});
  const NamedFunc l4_mini("l4_mini",[](const Baby &b) -> NamedFunc::ScalarType{return l4_Imini_d(b);});
  const NamedFunc l4_flav("l4_flav",[](const Baby &b) -> NamedFunc::ScalarType{return L4_flav_i(b);});

  const NamedFunc max_Imini("max_Imini",[](const Baby &b) -> NamedFunc::ScalarType{return max_Imini_d(b);});
  const NamedFunc ptj_m4j("ptj_m4j",[](const Baby &b) -> NamedFunc::ScalarType{return ptj_m4j_d(b);});
  const NamedFunc ptj_m5j("ptj_m5j",[](const Baby &b) -> NamedFunc::ScalarType{return ptj_m5j_d(b);});
  const NamedFunc max_deepflav("max_deepflav",[](const Baby &b) -> NamedFunc::ScalarType{return max_deepflav_d(b);});
  const NamedFunc max_particlenet("max_particlenet",[](const Baby &b) -> NamedFunc::ScalarType{return max_particlenet_d(b);});
  const NamedFunc mtop("mtop",[](const Baby &b) -> NamedFunc::ScalarType{return m_top(b);});

  const NamedFunc mtop_r3("mtop_r3",[](const Baby &b) -> NamedFunc::ScalarType{return m_top_r3(b);});
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
  const std::vector<NamedFunc>   run2_categories = {run2_ggF, run2_VBF, run2_lep};
  const std::vector<std::string> run2_cat_labs   = {"_run2_ggF", "_run2_VBF", "_run2_lep"};

  const NamedFunc cat_ttH_lep = "(nlep==3 && njet>=3 && nbdfm>=1) || (nlep>=4 && nbdfm>=1)";
  const NamedFunc cat_ttH_had = "(nlep==2 && njet>=5 && nbdfm>=1)";
  const NamedFunc cat_ZH_met  = "(nlep==2 && njet<=1 && met>90)"; 
  const NamedFunc cat_WH_3l   = "(nlep>=3 && nbdfm==0)";
  const NamedFunc cat_VBF     = "(nlep==2 && njet>=2 && nbdfm==0)";
  const NamedFunc cat_ggF     = "(nlep==2 && njet<=1 && met<90)";

  const NamedFunc baseline_ttH_lep =  max_Imini < 0.1 && "ll_m[0] > 85 && ll_m[0] < 95";
  const NamedFunc baseline_ttH_had = "ll_m[0] > 85 && ll_m[0] < 95";
  const NamedFunc baseline_ZH_met  = pTlly_mlly > 0.4 && "ll_m[0] > 85 && ll_m[0] < 95"; // && llgamma_pTt < 60 ; 
  const NamedFunc baseline_WH_3l   = max_Imini < 0.15 && pTlly_mlly > 0.3 && "ll_m[0] > 85 && ll_m[0] < 95" && "met > 30";
  const NamedFunc baseline_VBF     = "dijet_m > 600";

  const NamedFunc baseline_refit_ttH_lep =  max_Imini < 0.1 && "ll_refit_m > 85 && ll_refit_m < 95";
  const NamedFunc baseline_refit_ttH_had = "ll_refit_m > 85 && ll_refit_m < 95";
  const NamedFunc baseline_refit_ZH_met  = pTlly_mlly_refit > 0.4 && "ll_refit_m > 85 && ll_refit_m < 95";
  const NamedFunc baseline_refit_WH_3l   = max_Imini < 0.15 && pTlly_mlly_refit > 0.3 && "ll_refit_m > 85 && ll_refit_m < 95" && "met > 30";
  const NamedFunc baseline_refit_VBF     = "dijet_m > 600";

  const NamedFunc baseline_nomll_ttH_lep =  max_Imini < 0.1;
  const NamedFunc baseline_nomll_ttH_had = "1";
  const NamedFunc baseline_nomll_ZH_met  = pTlly_mlly > 0.4; 
  const NamedFunc baseline_nomll_WH_3l   = max_Imini < 0.15 && pTlly_mlly > 0.3 && "met > 30";
  const NamedFunc baseline_nomll_VBF     = "dijet_m > 600";


  const std::vector<std::vector<NamedFunc>> categories_selections_vector = {{max_Imini < 0.1, "ll_m[0] > 85 && ll_m[0] < 95"}, 
                                                                            {"ll_m[0] > 85 && ll_m[0] < 95"},
                                                                            {pTlly_mlly > 0.4, "ll_m[0] > 85 && ll_m[0] < 95"}, 
                                                                            {max_Imini < 0.15, pTlly_mlly > 0.3, "ll_m[0] > 85 && ll_m[0] < 95", "met > 30"},
                                                                            {"dijet_m > 600"}};

  const std::vector<std::vector<NamedFunc>> categories_refit_selections_vector = {{max_Imini < 0.1, "ll_refit_m > 85 && ll_refit_m < 95"}, 
                                                                                  {"ll_refit_m > 85 && ll_refit_m < 95"},
                                                                                  {pTlly_mlly_refit > 0.4, "ll_refit_m > 85 && ll_refit_m < 95"}, 
                                                                                  {max_Imini < 0.15, pTlly_mlly_refit > 0.3, "ll_refit_m > 85 && ll_refit_m < 95", "met > 30"},
                                                                                  {"dijet_m > 600"}};
  
  const std::vector<NamedFunc> run3_category_vector       = {cat_ttH_lep, cat_ttH_had, cat_ZH_met, cat_WH_3l, cat_VBF, cat_ggF};
  const std::vector<NamedFunc> run3_catwsel_vector        = {cat_ttH_lep && baseline_ttH_lep, cat_ttH_had && baseline_ttH_had,
                                                             cat_ZH_met && baseline_ZH_met, cat_WH_3l && baseline_WH_3l, cat_VBF, cat_ggF};

  const std::vector<NamedFunc> run3_refit_catwsel_vector  = {cat_ttH_lep && baseline_refit_ttH_lep, cat_ttH_had && baseline_refit_ttH_had,
                                                             cat_ZH_met && baseline_refit_ZH_met, cat_WH_3l && baseline_refit_WH_3l, cat_VBF, cat_ggF};

  const std::vector<NamedFunc> run3_catwsel_nomll_vector  = {cat_ttH_lep && baseline_nomll_ttH_lep, cat_ttH_had && baseline_nomll_ttH_had,
                                                             cat_ZH_met && baseline_nomll_ZH_met, cat_WH_3l && baseline_nomll_WH_3l, cat_VBF, cat_ggF}; 

  const std::vector<std::string> run3_category_labels = {"cat_ttH_lep","cat_ttH_had","cat_ZH_MET","cat_WH_3l","cat_VBF","cat_ggF"}; 

  const std::vector<NamedFunc> run3_catbs_vector = {baseline_ttH_lep, baseline_ttH_had, baseline_ZH_met, baseline_WH_3l, "1", "1"};
  const std::vector<NamedFunc> run3_catbs_nomll_vector  = {baseline_nomll_ttH_lep, baseline_nomll_ttH_had, baseline_nomll_ZH_met, baseline_nomll_WH_3l, "1", "1"}; 


  //Some repeat functions not used outside of categorization utilities  
  const NamedFunc pTy_mlly("pTy_mlly",[](const Baby &b) -> NamedFunc::ScalarType{ return (b.photon_pt()->at(0))/ZgUtilities::Findlly(b); });
  const NamedFunc mlly("mlly",[](const Baby &b) -> NamedFunc::ScalarType{         return AssignH(b).M(); });
  const NamedFunc pT_lly("pT_lly",[](const Baby &b) -> NamedFunc::ScalarType{     return AssignH(b).Pt(); });
  const NamedFunc eta_lly("eta_lly",[](const Baby &b) -> NamedFunc::ScalarType{   return AssignH(b).Eta(); });
  const NamedFunc rap_lly("rap_lly",[](const Baby &b) -> NamedFunc::ScalarType{   return AssignH(b).Rapidity(); });
  const NamedFunc phi_lly("phi_lly",[](const Baby &b) -> NamedFunc::ScalarType{   return AssignH(b).Phi(); });
  const NamedFunc sigy_pTy("sigy_pTy",[](const Baby &b) -> NamedFunc::ScalarType{ return (b.photon_energyErr() -> at(0))/(b.photon_pt()->at(0)); });

  const NamedFunc rap_ll("rap_ll",[](const Baby &b) -> NamedFunc::ScalarType{   return AssignZ(b).Rapidity(); });


  const NamedFunc dR_lly_dijet("dR_lly_dijet",[](const Baby &b) -> NamedFunc::ScalarType{ return llgamma_dijet_deltaR_d(b); });
  const NamedFunc dphi_lly_dijet("dphi_lly_dijet",[](const Baby &b) -> NamedFunc::ScalarType{ return fabs(llgamma_dijet_dphi_d(b)); });
  const NamedFunc diphi_dijet("diphi_dijet",[](const Baby &b) -> NamedFunc::ScalarType{ return fabs(b.dijet_dphi()); });
  const NamedFunc balance_lly_dijet("balance_lly_dijet",[](const Baby &b) -> NamedFunc::ScalarType{ return llgamma_dijet_balance_d(b); });
  const NamedFunc mindR_y_jet("mindR_y_jet",[](const Baby &b) -> NamedFunc::ScalarType{ return photon_jet_mindr_d(b); });
  const NamedFunc zeppenfeld_y_jet("zeppenfeld_y_jet",[](const Baby &b) -> NamedFunc::ScalarType{ return photon_zeppenfeld_d(b); });

  const NamedFunc l1_eta("l1_eta",[](const Baby &b) -> NamedFunc::ScalarType{
      if (b.nll()>0){
	if (b.ll_lepid()->at(0) == 11)
	  return b.el_eta()->at(b.ll_i1()->at(0));
	else if (b.ll_lepid()->at(0) == 13)
	  return b.mu_eta()->at(b.ll_i1()->at(0));
      }
      return 0;
    });

  const NamedFunc l2_eta("l2_eta",[](const Baby &b) -> NamedFunc::ScalarType{
      if (b.nll()>0){
	if (b.ll_lepid()->at(0) == 11)
	  return b.el_eta()->at(b.ll_i2()->at(0));
	else if (b.ll_lepid()->at(0) == 13)
	  return b.mu_eta()->at(b.ll_i2()->at(0));
      }
      return 0;
    });


  const NamedFunc rap_lly_refit("rap_lly_refit",[](const Baby &b) -> NamedFunc::ScalarType{
    TLorentzVector h; h.SetPtEtaPhiM(b.llphoton_refit_pt(), b.llphoton_refit_eta(), b.llphoton_refit_phi(), b.llphoton_refit_m());
    return h.Rapidity(); 
    });


  const NamedFunc lly_costheta("lly_costheta",[](const Baby &b) -> NamedFunc::ScalarType{
    if(b.nllphoton()>0){
      return b.llphoton_costheta() -> at(0);
    }
    return ZgUtilities::cos_theta(b);
  });

  const NamedFunc lly_cosTheta("lly_cosTheta",[](const Baby &b) -> NamedFunc::ScalarType{
    if(b.nllphoton()>0){
      return b.llphoton_cosTheta() -> at(0);
    }
    return ZgUtilities::cos_Theta(b);
  });

  const NamedFunc lly_angphi("lly_angphi",[](const Baby &b) -> NamedFunc::ScalarType{
    if(b.nllphoton()>0){
      return b.llphoton_psi() -> at(0);
    }
    return ZgUtilities::Getphi(b);
  });

  //Photon plots added to each category's plotting function
  void sample_photon_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    pm.Push<Hist1D>(Axis(45,0.0,1,      "photon_idmva[0]", "#gamma - IDMVA",            {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_IDMVA_");
    pm.Push<Hist1D>(Axis(35,0,3.5,      "photon_drmin[0]", "min( #Delta R(#gamma,l) )", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_minDR_");
    pm.Push<Hist1D>(Axis(60,0,6.0,      "photon_drmax[0]", "max( #Delta R(#gamma,l) )", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_maxDR_");
    pm.Push<Hist1D>(Axis(40,10,60,      "photon_pt[0]",    "p_{T}(#gamma) [GeV]",       {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_pT_");
    pm.Push<Hist1D>(Axis(40,-2.6,2.6,   "photon_eta[0]",   "#eta(#gamma)",              {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_eta_");
    pm.Push<Hist1D>(Axis(63,-3.15,3.15, "photon_phi[0]",   "#phi(#gamma)",              {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_phi_");
    pm.Push<Hist1D>(Axis(40,0.01,0.25,  "photon_energyErr[0]/photon_pt[0]", "#sigma(#gamma)/E(#gamma)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_relunc_");
  }

  //Photon plots added to each category's plotting function
  void sample_photon_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, int nbins, std::string labels){
    pm.Push<Hist1D>(Axis(nbins,0.0,1,      "photon_idmva[0]", "#gamma - IDMVA",            {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_IDMVA_");
    pm.Push<Hist1D>(Axis(nbins,0,3.5,      "photon_drmin[0]", "min( #Delta R(#gamma,l) )", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_minDR_");
    pm.Push<Hist1D>(Axis(nbins,0,6.0,      "photon_drmax[0]", "max( #Delta R(#gamma,l) )", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_maxDR_");
    pm.Push<Hist1D>(Axis(nbins,10,60,      "photon_pt[0]",    "p_{T}(#gamma) [GeV]",       {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_pT_");
    pm.Push<Hist1D>(Axis(nbins,-2.6,2.6,   "photon_eta[0]",   "#eta(#gamma)",              {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_eta_");
    pm.Push<Hist1D>(Axis(nbins,-3.15,3.15, "photon_phi[0]",   "#phi(#gamma)",              {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_phi_");
    pm.Push<Hist1D>(Axis(nbins,0.01,0.25,  "photon_energyErr[0]/photon_pt[0]", "#sigma(#gamma)/E(#gamma)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_relunc_");
  }

  //ll and llphoton plots added to each category's photon plots
  void sample_llphoton_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    //checks for a specific case for control regions where we want 50 < mlly < 100
    if(Contains(labels,"Zfsrpeak")){
      pm.Push<Hist1D>(Axis(50,50,100,   mlly,        "m_{ll#gamma} [GeV]",     {120,130}), selection, processes ,ops).Weight(wgt).Tag("ShortName:" + labels + "_lly_m_");
    }else{
      pm.Push<Hist1D>(Axis(40,100,180,  mlly,        "m_{ll#gamma} [GeV]",     {120,130}), selection, processes ,ops).Weight(wgt).Tag("ShortName:" + labels + "_lly_m_");
    }
    pm.Push<Hist1D>(Axis(50,0,250,      pT_lly,      "p_{T}(ll#gamma) [GeV]",         {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_lly_pt_");
    pm.Push<Hist1D>(Axis(52,-5,5,       eta_lly,     "#eta(ll#gamma)",                {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_lly_eta_");
    pm.Push<Hist1D>(Axis(52,-3,3,       rap_lly,     "y(ll#gamma)",                   {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_lly_rap_");
    pm.Push<Hist1D>(Axis(63,-3.15,3.15, phi_lly,     "#phi(ll#gamma)",                {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_lly_phi_");
    pm.Push<Hist1D>(Axis(40,0,2.5,      pTlly_mlly,  "p_{T}(ll#gamma)/m_{ll#gamma}", {1}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_lly_ptom_");
    pm.Push<Hist1D>(Axis(40,0,1.2,      pTy_mlly,    "p_{T}(#gamma)/m_{ll#gamma}",    {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_y_ptom_");
    pm.Push<Hist1D>(Axis(45,0,4.5,      lly_dR,      "#Delta R(ll,#gamma)",           {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ll_y_dR_");

    pm.Push<Hist1D>(Axis(50,0,150,      "ll_pt[0]",  "p_{T}(ll) [GeV]", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ll_pt_");
    pm.Push<Hist1D>(Axis(50,-5, 5,      "ll_eta[0]", "#eta(ll)",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ll_eta_");
    pm.Push<Hist1D>(Axis(52,-3,3,       rap_ll,     "y(ll)",            {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ll_rap_");
    pm.Push<Hist1D>(Axis(63,-3.15,3.15, "ll_phi[0]", "#phi(ll)",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ll_phi_");
    pm.Push<Hist1D>(Axis(50,80,100,     "ll_m[0]",   "m_{ll} [GeV]",    {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ll_m_");
  }

  void sample_llphoton_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, int nbins, std::string labels){
    //checks for a specific case for control regions where we want 50 < mlly < 100
    if(Contains(labels,"Zfsrpeak")){
      pm.Push<Hist1D>(Axis(nbins, 50, 100,   mlly,        "m_{ll#gamma} [GeV]", {120,130}), selection, processes ,ops).Weight(wgt).Tag("ShortName:" + labels + "_lly_m_");
      pm.Push<Hist1D>(Axis(nbins, 50,  80,   "ll_m[0]",   "m_{ll} [GeV]",       {}),        selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ll_m_");
    }else{
      int nbins_mlly = nbins;
      if(nbins < 16){ nbins_mlly = 8; 
      } else if(nbins < 20){ nbins_mlly=16; 
      }else if(nbins < 40){ nbins_mlly=20; }

      pm.Push<Hist1D>(Axis(nbins_mlly, 100, 180,  mlly,      "m_{ll#gamma} [GeV]", {120,130}), selection, processes ,ops).Weight(wgt).Tag("ShortName:" + labels + "_lly_m_");
      pm.Push<Hist1D>(Axis(nbins,       80, 100,  "ll_m[0]", "m_{ll} [GeV]",       {}),        selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ll_m_");
    }
    pm.Push<Hist1D>(Axis(nbins,0,250,      pT_lly,      "p_{T}(ll#gamma) [GeV]",         {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_lly_pt_");
    pm.Push<Hist1D>(Axis(nbins,-5,5,       eta_lly,     "#eta(ll#gamma)",                {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_lly_eta_");
    pm.Push<Hist1D>(Axis(nbins,-3,3,       rap_lly,     "y(ll#gamma)",                   {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_lly_rap_");
    pm.Push<Hist1D>(Axis(nbins,-3.15,3.15, phi_lly,     "#phi(ll#gamma)",                {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_lly_phi_");
    pm.Push<Hist1D>(Axis(nbins,0,2.5,      pTlly_mlly,  "p_{T}(ll#gamma)/m_{ll#gamma}", {1}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_lly_ptom_");
    pm.Push<Hist1D>(Axis(nbins,0,1.2,      pTy_mlly,    "p_{T}(#gamma)/m_{ll#gamma}",    {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_y_ptom_");
    pm.Push<Hist1D>(Axis(nbins,0,4.5,      lly_dR,      "#Delta R(ll,#gamma)",           {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ll_y_dR_");

    pm.Push<Hist1D>(Axis(nbins,0,150,      "ll_pt[0]",  "p_{T}(ll) [GeV]", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ll_pt_");
    pm.Push<Hist1D>(Axis(nbins,-5, 5,      "ll_eta[0]", "#eta(ll)",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ll_eta_");
    pm.Push<Hist1D>(Axis(nbins,-3,3,       rap_ll,     "y(ll)",            {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ll_rap_");
    pm.Push<Hist1D>(Axis(nbins,-3.15,3.15, "ll_phi[0]", "#phi(ll)",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ll_phi_");

  }


  void sample_kinrefit_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    //checks for a specific case for control regions where we want 50 < mlly < 100
    pm.Push<Hist1D>(Axis(80,100,180,   "llphoton_refit_m",                   "m_{ll#gamma} [GeV]",     {120,130}), selection, processes ,ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_lly_m_");
    pm.Push<Hist1D>(Axis(50,0,250,     "llphoton_refit_pt",                  "p_{T}(ll#gamma) [GeV]",         {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_lly_pt_");
    pm.Push<Hist1D>(Axis(52,-5,5,      "llphoton_refit_eta",                 "#eta(ll#gamma)",                {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_lly_eta_");
    pm.Push<Hist1D>(Axis(52,-2.5,2.5,  rap_lly_refit,                        "y(ll#gamma)",                   {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_lly_eta_");
    pm.Push<Hist1D>(Axis(63,-3.15,3.15,"llphoton_refit_phi",                 "#phi(ll#gamma)",                {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_lly_phi_");
    pm.Push<Hist1D>(Axis(40,0,2.5,     "llphoton_refit_pt/llphoton_refit_m", "p_{T}(ll#gamma)/m_{ll#gamma}", {1}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_lly_ptom_");
    pm.Push<Hist1D>(Axis(40,0,1.2,     "photon_pt[0]/llphoton_refit_m",      "p_{T}(#gamma)/m_{ll#gamma}",    {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_y_ptom_");
    pm.Push<Hist1D>(Axis(45,0,4.5,     "llphoton_refit_dr",                  "#Delta R(ll,#gamma)",           {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_ll_y_dR_");
    pm.Push<Hist1D>(Axis(45,-3.15,3.15,"llphoton_refit_dphi",                "#Delta R(ll,#gamma)",           {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_ll_y_dphi_");

    pm.Push<Hist1D>(Axis(40,-1,1,      "llphoton_refit_cosTheta", "cos(#Theta)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_cTHETA_");
    pm.Push<Hist1D>(Axis(40,-1,1,      "llphoton_refit_costheta", "cos(#theta)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_ctheta_");
    pm.Push<Hist1D>(Axis(40,-3.2,3.2,  "llphoton_refit_psi",      "#phi",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_phi_");
    pm.Push<Hist1D>(Axis(40,0,140,     "llphoton_refit_pTt",      "p^{t}_{T}",   {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_pTt_");

    pm.Push<Hist1D>(Axis(50,0,150,      "ll_refit_pt",  "p_{T}(ll) [GeV]",{}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_ll_pt_");
    pm.Push<Hist1D>(Axis(50,-5, 5,      "ll_refit_eta", "#eta(ll)",       {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_ll_eta_");
    pm.Push<Hist1D>(Axis(63,-3.15,3.15, "ll_refit_phi", "#phi(ll)",       {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_ll_phi_");
    pm.Push<Hist1D>(Axis(50,80,100,     "ll_refit_m",   "m_{ll} [GeV]",   {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_ll_m_");

  }

  void sample_kinrefit_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, int nbins, std::string labels){
    //checks for a specific case for control regions where we want 50 < mlly < 100
    pm.Push<Hist1D>(Axis(nbins,100,180,   "llphoton_refit_m",                   "m_{ll#gamma} [GeV]",     {120,130}), selection, processes ,ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_lly_m_");
    pm.Push<Hist1D>(Axis(nbins,0,250,     "llphoton_refit_pt",                  "p_{T}(ll#gamma) [GeV]",         {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_lly_pt_");
    pm.Push<Hist1D>(Axis(nbins,-5,5,      "llphoton_refit_eta",                 "#eta(ll#gamma)",                {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_lly_eta_");
    pm.Push<Hist1D>(Axis(nbins,-2.5,2.5,  rap_lly_refit,                        "y(ll#gamma)",                   {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_lly_eta_");
    pm.Push<Hist1D>(Axis(nbins,-3.15,3.15,"llphoton_refit_phi",                 "#phi(ll#gamma)",                {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_lly_phi_");
    pm.Push<Hist1D>(Axis(nbins,0,2.5,     "llphoton_refit_pt/llphoton_refit_m", "p_{T}(ll#gamma)/m_{ll#gamma}", {1}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_lly_ptom_");
    pm.Push<Hist1D>(Axis(nbins,0,1.2,     "photon_pt[0]/llphoton_refit_m",      "p_{T}(#gamma)/m_{ll#gamma}",    {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_y_ptom_");
    pm.Push<Hist1D>(Axis(nbins,0,4.5,     "llphoton_refit_dr",                  "#Delta R(ll,#gamma)",           {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_ll_y_dR_");
    pm.Push<Hist1D>(Axis(nbins,-3.15,3.15,"llphoton_refit_dphi",                "#Delta R(ll,#gamma)",           {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_ll_y_dphi_");

    pm.Push<Hist1D>(Axis(nbins,-1,1,      "llphoton_refit_cosTheta", "cos(#Theta)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_cTHETA_");
    pm.Push<Hist1D>(Axis(nbins,-1,1,      "llphoton_refit_costheta", "cos(#theta)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_ctheta_");
    pm.Push<Hist1D>(Axis(nbins,-3.2,3.2,  "llphoton_refit_psi",      "#phi",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_phi_");
    pm.Push<Hist1D>(Axis(nbins,0,140,     "llphoton_refit_pTt",      "p^{t}_{T}",   {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_pTt_");

    pm.Push<Hist1D>(Axis(nbins,0,150,      "ll_refit_pt",  "p_{T}(ll) [GeV]",{}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_ll_pt_");
    pm.Push<Hist1D>(Axis(nbins,-5, 5,      "ll_refit_eta", "#eta(ll)",       {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_ll_eta_");
    pm.Push<Hist1D>(Axis(nbins,-3.15,3.15, "ll_refit_phi", "#phi(ll)",       {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_ll_phi_");
    pm.Push<Hist1D>(Axis(nbins,80,100,     "ll_refit_m",   "m_{ll} [GeV]",   {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_refit_ll_m_");

  }

  //Basic lepton kinematics plots
  void sample_lepton_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    pm.Push<Hist1D>(Axis(50,     5,  100, "el_pt[ll_i1[0]]",  "p_{T}(e_{1})",   {}), selection && "ll_lepid[0]==11", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_pt_el1_");
    pm.Push<Hist1D>(Axis(40,  -2.6,  2.6, "el_eta[ll_i1[0]]", "#eta(e_{1})",    {}), selection && "ll_lepid[0]==11", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_eta_el1_");
    pm.Push<Hist1D>(Axis(62, -3.15, 3.15, "el_phi[ll_i1[0]]", "#phi(e_{1})",    {}), selection && "ll_lepid[0]==11", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_phi_el1_");

    pm.Push<Hist1D>(Axis(40,     5,   65, "el_pt[ll_i2[0]]",  "p_{T}(e_{2})",   {}), selection && "ll_lepid[0]==11", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_pt_el2_");
    pm.Push<Hist1D>(Axis(40,  -2.6,  2.6, "el_eta[ll_i2[0]]", "#eta(e_{2})",    {}), selection && "ll_lepid[0]==11", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_eta_el2_");
    pm.Push<Hist1D>(Axis(62, -3.15, 3.15, "el_phi[ll_i2[0]]", "#phi(e_{2})",    {}), selection && "ll_lepid[0]==11", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_phi_el2_");

    pm.Push<Hist1D>(Axis(65,     0,  110, "mu_pt[ll_i1[0]]",  "p_{T}(#mu_{1})", {}), selection && "ll_lepid[0]==13", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_pt_mu1_");
    pm.Push<Hist1D>(Axis(40,  -2.6,  2.6, "mu_eta[ll_i1[0]]", "#eta(#mu_{1})",  {}), selection && "ll_lepid[0]==13", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_eta_mu1_");
    pm.Push<Hist1D>(Axis(62, -3.15, 3.15, "mu_phi[ll_i1[0]]", "#phi(#mu_{1})",  {}), selection && "ll_lepid[0]==13", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_phi_mu1_");

    pm.Push<Hist1D>(Axis(40,     0,   75, "mu_pt[ll_i2[0]]",  "p_{T}(#mu_{2})", {}), selection && "ll_lepid[0]==13", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_pt_mu2_");
    pm.Push<Hist1D>(Axis(40,  -2.6,  2.6, "mu_eta[ll_i2[0]]", "#eta(#mu_{2})",  {}), selection && "ll_lepid[0]==13", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_eta_mu2_");
    pm.Push<Hist1D>(Axis(62, -3.15, 3.15, "mu_phi[ll_i2[0]]", "#phi(#mu_{2})",  {}), selection && "ll_lepid[0]==13", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_phi_mu2_");
  }

  void sample_lepton_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, int nbins, std::string labels){
    pm.Push<Hist1D>(Axis(nbins,     5,  100, "el_pt[ll_i1[0]]",  "p_{T}(e_{1})",   {}), selection && "ll_lepid[0]==11", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_pt_el1_");
    pm.Push<Hist1D>(Axis(nbins,  -2.6,  2.6, "el_eta[ll_i1[0]]", "#eta(e_{1})",    {}), selection && "ll_lepid[0]==11", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_eta_el1_");
    pm.Push<Hist1D>(Axis(nbins, -3.15, 3.15, "el_phi[ll_i1[0]]", "#phi(e_{1})",    {}), selection && "ll_lepid[0]==11", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_phi_el1_");

    pm.Push<Hist1D>(Axis(nbins,     5,   65, "el_pt[ll_i2[0]]",  "p_{T}(e_{2})",   {}), selection && "ll_lepid[0]==11", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_pt_el2_");
    pm.Push<Hist1D>(Axis(nbins,  -2.6,  2.6, "el_eta[ll_i2[0]]", "#eta(e_{2})",    {}), selection && "ll_lepid[0]==11", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_eta_el2_");
    pm.Push<Hist1D>(Axis(nbins, -3.15, 3.15, "el_phi[ll_i2[0]]", "#phi(e_{2})",    {}), selection && "ll_lepid[0]==11", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_phi_el2_");

    pm.Push<Hist1D>(Axis(nbins,     0,  110, "mu_pt[ll_i1[0]]",  "p_{T}(#mu_{1})", {}), selection && "ll_lepid[0]==13", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_pt_mu1_");
    pm.Push<Hist1D>(Axis(nbins,  -2.6,  2.6, "mu_eta[ll_i1[0]]", "#eta(#mu_{1})",  {}), selection && "ll_lepid[0]==13", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_eta_mu1_");
    pm.Push<Hist1D>(Axis(nbins, -3.15, 3.15, "mu_phi[ll_i1[0]]", "#phi(#mu_{1})",  {}), selection && "ll_lepid[0]==13", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_phi_mu1_");

    pm.Push<Hist1D>(Axis(nbins,     0,   75, "mu_pt[ll_i2[0]]",  "p_{T}(#mu_{2})", {}), selection && "ll_lepid[0]==13", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_pt_mu2_");
    pm.Push<Hist1D>(Axis(nbins,  -2.6,  2.6, "mu_eta[ll_i2[0]]", "#eta(#mu_{2})",  {}), selection && "ll_lepid[0]==13", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_eta_mu2_");
    pm.Push<Hist1D>(Axis(nbins, -3.15, 3.15, "mu_phi[ll_i2[0]]", "#phi(#mu_{2})",  {}), selection && "ll_lepid[0]==13", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_phi_mu2_");
  }


  //Sample ttH leptonic category plots
  void ttH_lep_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels, bool run3){
    sample_photon_plots(  pm, selection, processes, ops, wgt, labels + "_ttH_lep_");
    sample_lepton_plots(  pm, selection, processes, ops, wgt, labels + "_ttH_lep_");
    sample_llphoton_plots(pm, selection, processes, ops, wgt, labels + "_ttH_lep_");
    //sample_kinrefit_plots(pm, selection, processes, ops, wgt, labels + "_ttH_lep_");
    //mvp_objects(          pm, selection, processes, ops, wgt, labels + "_ttH_lep_");

    pm.Push<Hist1D>(Axis(80,0,500, "ht",    "H_{T} [GeV]", {120}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_HT_"); 
    pm.Push<Hist1D>(Axis(75,0,200, "met",   "MET [GeV]",      {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_MET_");
    pm.Push<Hist1D>(Axis(10,0,10,  "njet",  "N_{jet}",       {3}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_Njet_");
    pm.Push<Hist1D>(Axis(6,0,6,    "nlep",  "N_{lep}",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_Nlep_");
    pm.Push<Hist1D>(Axis(6,0,6 ,   "nbdfm", "N_{b,med}",      {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_Nbdfm");

    if(run3==true){
      pm.Push<Hist1D>(Axis(40,0,1,   max_particlenet, "Particle Net - max", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_pnet_BASE_MAX");
    } else {
      pm.Push<Hist1D>(Axis(40,0,1,   max_deepflav, "Deep Flavor - max", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_deepflav_BASE_MAX");
    }

    pm.Push<Hist1D>(Axis(50,0,1,     l3_mini,   "I_{mini} - l3",            {}), selection,                processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_l3_Imini");
    pm.Push<Hist1D>(Axis(55,-0.1,1,  l3_idmva,  "l_{3} - e IDMVA ",         {}), selection && l3_flav==11, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_e3idmva"); 

    pm.Push<Hist1D>(Axis(35,     5,   75, l3_pt,  "p_{T}(e_{3}) [GeV]", {}), selection && l3_flav==11, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_e3_pt");
    pm.Push<Hist1D>(Axis(40,  -2.6,  2.6, l3_eta, "#eta(e_{3})",        {}), selection && l3_flav==11, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_e3_eta");
    pm.Push<Hist1D>(Axis(62, -3.15, 3.15, l3_phi, "#phi(e_{3})",        {}), selection && l3_flav==11, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_e3_phi");

    pm.Push<Hist1D>(Axis(35,     5,   75, l3_pt,  "p_{T}(#mu_{3}) [GeV]", {}), selection && l3_flav==13, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_mu3_pt");
    pm.Push<Hist1D>(Axis(40,  -2.6,  2.6, l3_eta, "#eta(#mu_{3})",        {}), selection && l3_flav==13, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_mu3_eta");
    pm.Push<Hist1D>(Axis(62, -3.15, 3.15, l3_phi, "#phi(#mu_{3})",        {}), selection && l3_flav==13, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_mu3_phi");

    pm.Push<Hist1D>(Axis(35,     5,   75, l3_pt,  "p_{T}(l_{3}) [GeV]", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_l3_pt");
    pm.Push<Hist1D>(Axis(40,  -2.6,  2.6, l3_eta, "#eta(l_{3})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_l3_eta");
    pm.Push<Hist1D>(Axis(62, -3.15, 3.15, l3_phi, "#phi(l_{3})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_l3_phi");

    pm.Push<Hist1D>(Axis(60,    30,  150, j1_pt,  "p_{T}(j_{1}) [GeV]", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_j1_pt");
    pm.Push<Hist1D>(Axis(50, -5.0,  5.0,  j1_eta, "#eta(j_{1})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_j1_eta_");
    pm.Push<Hist1D>(Axis(63, -3.15, 3.15, j1_phi, "#phi(j_{1})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_j1_phi_");
    pm.Push<Hist1D>(Axis(40,     0,   40, j1_m,   "m(j_{1}) [GeV]",     {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_j1_m");

    pm.Push<Hist1D>(Axis(60,    30,  150, j2_pt,  "p_{T}(j_{2}) [GeV]", {}), selection && "nlep==3", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_j2_pt");
    pm.Push<Hist1D>(Axis(50,  -5.0,  5.0, j2_eta, "#eta(j_{2})",        {}), selection && "nlep==3", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_j2_eta_");
    pm.Push<Hist1D>(Axis(63, -3.15, 3.15, j2_phi, "#phi(j_{2})",        {}), selection && "nlep==3", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_j2_phi_");
    pm.Push<Hist1D>(Axis(40,     0,   40, j2_m,   "m(j_{2}) [GeV]",     {}), selection && "nlep==3", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_j2_m");

    pm.Push<Hist1D>(Axis(60,    30,  150, j3_pt,  "p_{T}(j_{3}) [GeV]", {}), selection && "nlep==3", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_j3_pt");
    pm.Push<Hist1D>(Axis(50,  -5.0,  5.0, j3_eta, "#eta(j_{3})",        {}), selection && "nlep==3", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_j3_eta_");
    pm.Push<Hist1D>(Axis(63, -3.15, 3.15, j3_phi, "#phi(j_{3})",        {}), selection && "nlep==3", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_j3_phi_");
    pm.Push<Hist1D>(Axis(40,     0,   40, j3_m,   "m(j_{3}) [GeV]",     {}), selection && "nlep==3", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_j3_m");

    pm.Push<Hist1D>(Axis(45,0,4.5,   H_Z_dR,    "#Delta R(ll,#gamma)",      {}), selection,                processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_H_Z_dR_");
    pm.Push<Hist1D>(Axis(50,0,1,     max_Imini, "max. I_{mini}",            {}), selection,                processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_max_Imini");

    //pm.Push<Hist1D>(Axis(6,0,6 ,   "ntruel","N_{e,truth}",    {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_Ntruel");
    //pm.Push<Hist1D>(Axis(6,0,6 ,   "ntrumu","N_{#mu,truth}",  {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_Ntrumu");
    //pm.Push<Hist1D>(Axis(6,0,6 ,   Ntrub,   "N_{b,truth}",    {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_Ntrub");
    //pm.Push<Hist1D>(Axis(6,0,6 ,   Ntruel_NF,"N_{e,truth}",   {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_Ntruel_NF");
  }

  void ttH_lep_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels, int nbins, bool run3){
    sample_photon_plots(  pm, selection, processes, ops, wgt, nbins, labels + "_ttH_lep_");
    sample_lepton_plots(  pm, selection, processes, ops, wgt, nbins, labels + "_ttH_lep_");
    sample_llphoton_plots(pm, selection, processes, ops, wgt, nbins, labels + "_ttH_lep_");
    //sample_kinrefit_plots(pm, selection, processes, ops, wgt, nbins, labels + "_ttH_lep_");
    //mvp_objects(          pm, selection, processes, ops, wgt, labels + "_ttH_lep_");

    pm.Push<Hist1D>(Axis(nbins,0,500, "ht",    "H_{T} [GeV]", {120}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_HT_"); 
    pm.Push<Hist1D>(Axis(nbins,0,200, "met",   "MET [GeV]",      {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_MET_");
    pm.Push<Hist1D>(Axis(10,0,10,  "njet",  "N_{jet}",       {3}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_Njet_");
    pm.Push<Hist1D>(Axis(6,0,6,    "nlep",  "N_{lep}",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_Nlep_");
    pm.Push<Hist1D>(Axis(6,0,6 ,   "nbdfm", "N_{b,med}",      {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_Nbdfm");

    if(run3==true){
      pm.Push<Hist1D>(Axis(nbins,0,1,   max_particlenet, "Particle Net - max", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_pnet_BASE_MAX");
    } else {
      pm.Push<Hist1D>(Axis(nbins,0,1,   max_deepflav, "Deep Flavor - max", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_deepflav_BASE_MAX");
    }

    pm.Push<Hist1D>(Axis(nbins,0,1,     l3_mini,   "I_{mini} - l3",            {}), selection,                processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_l3_Imini");
    pm.Push<Hist1D>(Axis(nbins,-0.1,1,  l3_idmva,  "l_{3} - e IDMVA ",         {}), selection && l3_flav==11, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_e3idmva"); 

    pm.Push<Hist1D>(Axis(nbins,     5,   75, l3_pt,  "p_{T}(e_{3}) [GeV]", {}), selection && l3_flav==11, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_e3_pt");
    pm.Push<Hist1D>(Axis(nbins,  -2.6,  2.6, l3_eta, "#eta(e_{3})",        {}), selection && l3_flav==11, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_e3_eta");
    pm.Push<Hist1D>(Axis(nbins, -3.15, 3.15, l3_phi, "#phi(e_{3})",        {}), selection && l3_flav==11, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_e3_phi");

    pm.Push<Hist1D>(Axis(nbins,    5,   75, l3_pt,  "p_{T}(#mu_{3}) [GeV]", {}), selection && l3_flav==13, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_mu3_pt");
    pm.Push<Hist1D>(Axis(nbins, -2.6,  2.6, l3_eta, "#eta(#mu_{3})",        {}), selection && l3_flav==13, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_mu3_eta");
    pm.Push<Hist1D>(Axis(nbins,-3.15, 3.15, l3_phi, "#phi(#mu_{3})",        {}), selection && l3_flav==13, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_mu3_phi");

    pm.Push<Hist1D>(Axis(nbins,     5,   75, l3_pt,  "p_{T}(l_{3}) [GeV]", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_l3_pt");
    pm.Push<Hist1D>(Axis(nbins,  -2.6,  2.6, l3_eta, "#eta(l_{3})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_l3_eta");
    pm.Push<Hist1D>(Axis(nbins, -3.15, 3.15, l3_phi, "#phi(l_{3})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_l3_phi");

    pm.Push<Hist1D>(Axis(nbins,    30,  150, j1_pt,  "p_{T}(j_{1}) [GeV]", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_j1_pt");
    pm.Push<Hist1D>(Axis(nbins, -5.0,  5.0,  j1_eta, "#eta(j_{1})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_j1_eta_");
    pm.Push<Hist1D>(Axis(nbins, -3.15, 3.15, j1_phi, "#phi(j_{1})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_j1_phi_");
    pm.Push<Hist1D>(Axis(nbins,     0,   40, j1_m,   "m(j_{1}) [GeV]",     {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_j1_m");

    pm.Push<Hist1D>(Axis(nbins,    30,  150, j2_pt,  "p_{T}(j_{2}) [GeV]", {}), selection && "nlep==3", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_j2_pt");
    pm.Push<Hist1D>(Axis(nbins,  -5.0,  5.0, j2_eta, "#eta(j_{2})",        {}), selection && "nlep==3", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_j2_eta_");
    pm.Push<Hist1D>(Axis(nbins, -3.15, 3.15, j2_phi, "#phi(j_{2})",        {}), selection && "nlep==3", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_j2_phi_");
    pm.Push<Hist1D>(Axis(nbins,     0,   40, j2_m,   "m(j_{2}) [GeV]",     {}), selection && "nlep==3", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_j2_m");

    pm.Push<Hist1D>(Axis(nbins,    30,  150, j3_pt,  "p_{T}(j_{3}) [GeV]", {}), selection && "nlep==3", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_j3_pt");
    pm.Push<Hist1D>(Axis(nbins,  -5.0,  5.0, j3_eta, "#eta(j_{3})",        {}), selection && "nlep==3", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_j3_eta_");
    pm.Push<Hist1D>(Axis(nbins, -3.15, 3.15, j3_phi, "#phi(j_{3})",        {}), selection && "nlep==3", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_j3_phi_");
    pm.Push<Hist1D>(Axis(nbins,     0,   40, j3_m,   "m(j_{3}) [GeV]",     {}), selection && "nlep==3", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_j3_m");

    pm.Push<Hist1D>(Axis(nbins,0,4.5,   H_Z_dR,    "#Delta R(ll,#gamma)",      {}), selection,                processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_H_Z_dR_");
    pm.Push<Hist1D>(Axis(nbins,0,1,     max_Imini, "max. I_{mini}",            {}), selection,                processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_max_Imini");
  }

  //Sample ttH had category plots
  void ttH_had_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels, bool run3){
    sample_photon_plots(  pm, selection, processes, ops, wgt, labels + "_ttH_had_");
    sample_lepton_plots(  pm, selection, processes, ops, wgt, labels + "_ttH_had_");
    sample_llphoton_plots(pm, selection, processes, ops, wgt, labels + "_ttH_had_");
    //sample_kinrefit_plots(pm, selection, processes, ops, wgt, labels + "_ttH_had_");
    //mvp_objects(          pm, selection, processes, ops, wgt, labels + "_ttH_had_");

    pm.Push<Hist1D>(Axis(60,    30,  250, j1_pt,  "p_{T}(j_{1}) [GeV]", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j1_pt");
    pm.Push<Hist1D>(Axis(100, -5.0,  5.0, j1_eta, "#eta(j_{1})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j1_eta_");
    pm.Push<Hist1D>(Axis(63, -3.15, 3.15, j1_phi, "#phi(j_{1})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j1_phi_");
    pm.Push<Hist1D>(Axis(40,     0,   40, j1_m,   "m(j_{1}) [GeV]",     {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j1_m");

    pm.Push<Hist1D>(Axis(60,    30,  230, j2_pt,  "p_{T}(j_{2}) [GeV]", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j2_pt");
    pm.Push<Hist1D>(Axis(100, -5.0,  5.0, j2_eta, "#eta(j_{2})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j2_eta_");
    pm.Push<Hist1D>(Axis(63, -3.15, 3.15, j2_phi, "#phi(j_{2})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j2_phi_");
    pm.Push<Hist1D>(Axis(40,     0,   40, j2_m,   "m(j_{2}) [GeV]",     {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j2_m");

    pm.Push<Hist1D>(Axis(60,    30,  200, j3_pt,  "p_{T}(j_{3}) [GeV]", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j3_pt");
    pm.Push<Hist1D>(Axis(100, -5.0,  5.0, j3_eta, "#eta(j_{3})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j3_eta_");
    pm.Push<Hist1D>(Axis(63, -3.15, 3.15, j3_phi, "#phi(j_{3})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j3_phi_");
    pm.Push<Hist1D>(Axis(40,     0,   40, j3_m,   "m(j_{3}) [GeV]",     {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j3_m");

    pm.Push<Hist1D>(Axis(60,    30,  120, j4_pt,  "p_{T}(j_{4}) [GeV]", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j4_pt");
    pm.Push<Hist1D>(Axis(100, -5.0,  5.0, j4_eta, "#eta(j_{4})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j4_eta_");
    pm.Push<Hist1D>(Axis(63, -3.15, 3.15, j4_phi, "#phi(j_{4})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j4_phi_");
    pm.Push<Hist1D>(Axis(40,     0,   40, j4_m,   "m(j_{4}) [GeV]",     {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j4_m");

    pm.Push<Hist1D>(Axis(60,    30,  80, j5_pt,  "p_{T}(j_{5}) [GeV]", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j5_pt");
    pm.Push<Hist1D>(Axis(100, -5.0,  5.0, j5_eta, "#eta(j_{5})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j5_eta_");
    pm.Push<Hist1D>(Axis(63, -3.15, 3.15, j5_phi, "#phi(j_{5})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j5_phi_");
    pm.Push<Hist1D>(Axis(40,     0,   40, j5_m,   "m(j_{5}) [GeV]",     {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j5_m");

    pm.Push<Hist1D>(Axis(10,0,10,    "njet",       "N_{jet}",                             {4}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_Njet_");
    pm.Push<Hist1D>(Axis(75,0,200,   "met",        "MET",                                  {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_MET");
    pm.Push<Hist1D>(Axis(80,150,1000, "ht",        "H_{T}",                             {150}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_HT_");

    pm.Push<Hist1D>(Axis(50,0,1,     ptj_m4j,      "#Sigma(p_{T}(j_{i}))/m_{4j} [GeV]", {0.6}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_ptj_m4j");
    pm.Push<Hist1D>(Axis(50,0,1,     ptj_m5j,      "#Sigma(p_{T}(j_{i}))/m_{5j} [GeV]", {0.6}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_ptj_m5j");
    pm.Push<Hist1D>(Axis(50,0,1,     max_Imini,    "max. I_{mini}",                        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_max_Imini");

    if(run3==true){
      pm.Push<Hist1D>(Axis(40,  0,   1, max_particlenet, "Particle Net - max", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_pnet_BASE_MAX");
      pm.Push<Hist1D>(Axis(75, 50, 350, mtop_r3,                        "m_{t}", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_mtop_");
    } else {
      pm.Push<Hist1D>(Axis(40,  0,   1, max_deepflav, "Deep Flavor - max", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_lep_deepflav_BASE_MAX");
      pm.Push<Hist1D>(Axis(75, 50, 350, mtop,         "m_{t}",             {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_mtop_");
    }
  }

  void ttH_had_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels, int nbins, bool run3){
    sample_photon_plots(  pm, selection, processes, ops, wgt, nbins, labels + "_ttH_had_");
    sample_lepton_plots(  pm, selection, processes, ops, wgt, nbins, labels + "_ttH_had_");
    sample_llphoton_plots(pm, selection, processes, ops, wgt, nbins, labels + "_ttH_had_");
    //sample_kinrefit_plots(pm, selection, processes, ops, wgt, nbins, labels + "_ttH_had_");
    //mvp_objects(          pm, selection, processes, ops, wgt, labels + "_ttH_had_");

    pm.Push<Hist1D>(Axis(nbins,    30,  250, j1_pt,  "p_{T}(j_{1}) [GeV]", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j1_pt");
    pm.Push<Hist1D>(Axis(nbins, -5.0,  5.0, j1_eta, "#eta(j_{1})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j1_eta_");
    pm.Push<Hist1D>(Axis(nbins, -3.15, 3.15, j1_phi, "#phi(j_{1})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j1_phi_");
    pm.Push<Hist1D>(Axis(nbins,     0,   40, j1_m,   "m(j_{1}) [GeV]",     {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j1_m");

    pm.Push<Hist1D>(Axis(nbins,    30,  230, j2_pt,  "p_{T}(j_{2}) [GeV]", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j2_pt");
    pm.Push<Hist1D>(Axis(nbins, -5.0,  5.0, j2_eta, "#eta(j_{2})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j2_eta_");
    pm.Push<Hist1D>(Axis(nbins, -3.15, 3.15, j2_phi, "#phi(j_{2})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j2_phi_");
    pm.Push<Hist1D>(Axis(nbins,     0,   40, j2_m,   "m(j_{2}) [GeV]",     {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j2_m");

    pm.Push<Hist1D>(Axis(nbins,    30,  200, j3_pt,  "p_{T}(j_{3}) [GeV]", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j3_pt");
    pm.Push<Hist1D>(Axis(nbins, -5.0,  5.0, j3_eta, "#eta(j_{3})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j3_eta_");
    pm.Push<Hist1D>(Axis(nbins, -3.15, 3.15, j3_phi, "#phi(j_{3})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j3_phi_");
    pm.Push<Hist1D>(Axis(nbins,     0,   40, j3_m,   "m(j_{3}) [GeV]",     {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j3_m");

    pm.Push<Hist1D>(Axis(nbins,    30,  120, j4_pt,  "p_{T}(j_{4}) [GeV]", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j4_pt");
    pm.Push<Hist1D>(Axis(nbins, -5.0,   5.0, j4_eta, "#eta(j_{4})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j4_eta_");
    pm.Push<Hist1D>(Axis(nbins, -3.15, 3.15, j4_phi, "#phi(j_{4})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j4_phi_");
    pm.Push<Hist1D>(Axis(nbins,     0,   40, j4_m,   "m(j_{4}) [GeV]",     {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j4_m");

    pm.Push<Hist1D>(Axis(nbins,    30,   80, j5_pt,  "p_{T}(j_{5}) [GeV]", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j5_pt");
    pm.Push<Hist1D>(Axis(nbins, -5.0,   5.0, j5_eta, "#eta(j_{5})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j5_eta_");
    pm.Push<Hist1D>(Axis(nbins, -3.15, 3.15, j5_phi, "#phi(j_{5})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j5_phi_");
    pm.Push<Hist1D>(Axis(nbins,     0,   40, j5_m,   "m(j_{5}) [GeV]",     {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_j5_m");

    pm.Push<Hist1D>(Axis(nbins,0,10,    "njet",       "N_{jet}",                             {4}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_Njet_");
    pm.Push<Hist1D>(Axis(nbins,0,200,   "met",        "MET",                                  {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_MET");
    pm.Push<Hist1D>(Axis(nbins,150,1000, "ht",        "H_{T}",                             {150}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_HT_");

    pm.Push<Hist1D>(Axis(nbins,0,1,     ptj_m4j,      "#Sigma(p_{T}(j_{i}))/m_{4j} [GeV]", {0.6}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_ptj_m4j");
    pm.Push<Hist1D>(Axis(nbins,0,1,     ptj_m5j,      "#Sigma(p_{T}(j_{i}))/m_{5j} [GeV]", {0.6}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_ptj_m5j");
    pm.Push<Hist1D>(Axis(nbins,0,1,     max_Imini,    "max. I_{mini}",                        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_max_Imini");

    if(run3==true){
      pm.Push<Hist1D>(Axis(nbins,  0,   1, max_particlenet, "Particle Net - max", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_pnet_BASE_MAX");
      pm.Push<Hist1D>(Axis(nbins, 50, 350, mtop_r3,                        "m_{t}", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_mtop_");
    } else {
      pm.Push<Hist1D>(Axis(nbins,  0,   1, max_deepflav, "Deep Flavor - max", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_deepflav_BASE_MAX");
      pm.Push<Hist1D>(Axis(nbins, 50, 350, mtop,         "m_{t}",             {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ttH_had_mtop_");
    }
  }




  //Sample WH 3l category plots
  void WH_3l_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    sample_photon_plots(pm,   selection, processes, ops, wgt, labels + "_WH_3l_");
    sample_lepton_plots(pm,   selection, processes, ops, wgt, labels + "_WH_3l_");
    sample_llphoton_plots(pm, selection, processes, ops, wgt, labels + "_WH_3l_");
    //sample_kinrefit_plots(pm, selection, processes, ops, wgt, labels + "_WH_3l_");
    //mvp_objects(          pm, selection, processes, ops, wgt, labels + "_WH_3l_");

    pm.Push<Hist1D>(Axis(75,0,200,  "met",     "MET",                   {}), selection,                processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_MET");
    pm.Push<Hist1D>(Axis(40,0,3.2,  dphi_Hmet, "#Delta #phi(MET,H)", {2.0}), selection,                processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_dphi_Hmet");
    pm.Push<Hist1D>(Axis(40,0,4.0,  l3_Z_DR,   "#Delta R(Z,l_{3})",   {45}), selection,                processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_Zl3_DR");
    pm.Push<Hist1D>(Axis(50,0,1,    l3_mini,   "I_{mini} - l3 [GeV]",   {}), selection,                processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_lep_l3_Imini");
    pm.Push<Hist1D>(Axis(55,-0.1,1, l3_idmva,  "l_{3} - e IDMVA ",      {}), selection,                processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_e3idmva");
    pm.Push<Hist1D>(Axis(50,0,1,    max_Imini, "max. I_{mini}",         {}), selection,                processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_max_Imini");

    pm.Push<Hist1D>(Axis(30,     5,   75,  l3_pt,  "p_{T}(e_{3}) [GeV]", {}), selection && l3_flav==11, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_e3_pt");
    pm.Push<Hist1D>(Axis(5.2,  -2.6,  2.6, l3_eta, "#eta(e_{3})",        {}), selection && l3_flav==11, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_e3_eta");
    pm.Push<Hist1D>(Axis(6.3, -3.15, 3.15, l3_phi, "#phi(e_{3})",        {}), selection && l3_flav==11, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_e3_phi");

    pm.Push<Hist1D>(Axis(30,    5,   75,  l3_pt,  "p_{T}(#mu_{3}) [GeV]", {}), selection && l3_flav==13, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_mu3_pt");
    pm.Push<Hist1D>(Axis(5.2, -2.6,  2.6, l3_eta, "#eta(#mu_{3})",        {}), selection && l3_flav==13, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_mu3_eta");
    pm.Push<Hist1D>(Axis(6.3,-3.15, 3.15, l3_phi, "#phi(#mu_{3})",        {}), selection && l3_flav==13, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_mu3_phi");

    pm.Push<Hist1D>(Axis(30,     5,   75,  l3_pt,  "p_{T}(l_{3}) [GeV]", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_l3_pt");
    pm.Push<Hist1D>(Axis(5.2,  -2.6,  2.6, l3_eta, "#eta(l_{3})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_l3_eta");
    pm.Push<Hist1D>(Axis(6.3, -3.15, 3.15, l3_phi, "#phi(l_{3})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_l3_phi");

    pm.Push<Hist1D>(Axis(50,    30,  150, j1_pt,  "p_{T}(j_{1}) [GeV]", {}), selection && "njet>0", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_j1_pt");
    pm.Push<Hist1D>(Axis(50,  -5.0,  5.0, j1_eta, "#eta(j_{1})",        {}), selection && "njet>0", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_j1_eta_");
    pm.Push<Hist1D>(Axis(63, -3.15, 3.15, j1_phi, "#phi(j_{1})",        {}), selection && "njet>0", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_j1_phi_");
    pm.Push<Hist1D>(Axis(50,     0,   40, j1_m,   "m(j_{1}) [GeV]",     {}), selection && "njet>0", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_j1_m");

    pm.Push<Hist1D>(Axis(50,    30,  150, j2_pt,  "p_{T}(j_{2}) [GeV]", {}), selection && "njet>0", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_j2_pt");
    pm.Push<Hist1D>(Axis(50,  -5.0,  5.0, j2_eta, "#eta(j_{2})",        {}), selection && "njet>0", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_j2_eta_");
    pm.Push<Hist1D>(Axis(63, -3.15, 3.15, j2_phi, "#phi(j_{2})",        {}), selection && "njet>0", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_j2_phi_");
    pm.Push<Hist1D>(Axis(50,     0,   40, j2_m,   "m(j_{2}) [GeV]",     {}), selection && "njet>0", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_j2_m");
  }

  void WH_3l_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels, int nbins){
    sample_photon_plots(pm,   selection, processes, ops, wgt, nbins, labels + "_WH_3l_");
    sample_lepton_plots(pm,   selection, processes, ops, wgt, nbins, labels + "_WH_3l_");
    sample_llphoton_plots(pm, selection, processes, ops, wgt, nbins, labels + "_WH_3l_");
    //sample_kinrefit_plots(pm, selection, processes, ops, wgt, nbins, labels + "_WH_3l_");
    //mvp_objects(          pm, selection, processes, ops, wgt, labels + "_WH_3l_");

    pm.Push<Hist1D>(Axis(nbins,0,200,  "met",     "MET",                   {}), selection,                processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_MET");
    pm.Push<Hist1D>(Axis(nbins,0,3.2,  dphi_Hmet, "#Delta #phi(MET,H)", {2.0}), selection,                processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_dphi_Hmet");
    pm.Push<Hist1D>(Axis(nbins,0,4.0,  l3_Z_DR,   "#Delta R(Z,l_{3})",   {45}), selection,                processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_Zl3_DR");
    pm.Push<Hist1D>(Axis(nbins,0,1,    l3_mini,   "I_{mini} - l3 [GeV]",   {}), selection,                processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_lep_l3_Imini");
    pm.Push<Hist1D>(Axis(nbins,-0.1,1, l3_idmva,  "l_{3} - e IDMVA ",      {}), selection,                processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_e3idmva");
    pm.Push<Hist1D>(Axis(nbins,0,1,    max_Imini, "max. I_{mini}",         {}), selection,                processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_max_Imini");

    pm.Push<Hist1D>(Axis(nbins,     5,   75, l3_pt,  "p_{T}(e_{3}) [GeV]", {}), selection && l3_flav==11, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_e3_pt");
    pm.Push<Hist1D>(Axis(nbins,  -2.6,  2.6, l3_eta, "#eta(e_{3})",        {}), selection && l3_flav==11, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_e3_eta");
    pm.Push<Hist1D>(Axis(nbins, -3.15, 3.15, l3_phi, "#phi(e_{3})",        {}), selection && l3_flav==11, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_e3_phi");

    pm.Push<Hist1D>(Axis(nbins,    5,   75, l3_pt,  "p_{T}(#mu_{3}) [GeV]", {}), selection && l3_flav==13, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_mu3_pt");
    pm.Push<Hist1D>(Axis(nbins, -2.6,  2.6, l3_eta, "#eta(#mu_{3})",        {}), selection && l3_flav==13, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_mu3_eta");
    pm.Push<Hist1D>(Axis(nbins,-3.15, 3.15, l3_phi, "#phi(#mu_{3})",        {}), selection && l3_flav==13, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_mu3_phi");

    pm.Push<Hist1D>(Axis(nbins,     5,   75, l3_pt,  "p_{T}(l_{3}) [GeV]", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_l3_pt");
    pm.Push<Hist1D>(Axis(nbins,  -2.6,  2.6, l3_eta, "#eta(l_{3})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_l3_eta");
    pm.Push<Hist1D>(Axis(nbins, -3.15, 3.15, l3_phi, "#phi(l_{3})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_l3_phi");

    pm.Push<Hist1D>(Axis(nbins,    30,  150, j1_pt,  "p_{T}(j_{1}) [GeV]", {}), selection && "njet>0", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_j1_pt");
    pm.Push<Hist1D>(Axis(nbins,  -5.0,  5.0, j1_eta, "#eta(j_{1})",        {}), selection && "njet>0", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_j1_eta_");
    pm.Push<Hist1D>(Axis(nbins, -3.15, 3.15, j1_phi, "#phi(j_{1})",        {}), selection && "njet>0", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_j1_phi_");
    pm.Push<Hist1D>(Axis(nbins,     0,   40, j1_m,   "m(j_{1}) [GeV]",     {}), selection && "njet>0", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_j1_m");

    pm.Push<Hist1D>(Axis(nbins,    30,  150, j2_pt,  "p_{T}(j_{2}) [GeV]", {}), selection && "njet>1", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_j2_pt");
    pm.Push<Hist1D>(Axis(nbins,  -5.0,  5.0, j2_eta, "#eta(j_{2})",        {}), selection && "njet>1", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_j2_eta_");
    pm.Push<Hist1D>(Axis(nbins, -3.15, 3.15, j2_phi, "#phi(j_{2})",        {}), selection && "njet>1", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_j2_phi_");
    pm.Push<Hist1D>(Axis(nbins,     0,   40, j2_m,   "m(j_{2}) [GeV]",     {}), selection && "njet>1", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_WH_3l_j2_m");
  }

  //Sample ZH MET category plots
  void ZH_MET_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    sample_photon_plots(  pm, selection, processes, ops, wgt, labels + "_ZH_MET_");
    sample_lepton_plots(  pm, selection, processes, ops, wgt, labels + "_ZH_MET_");
    sample_llphoton_plots(pm, selection, processes, ops, wgt, labels + "_ZH_MET_");
    //sample_kinrefit_plots(pm, selection, processes, ops, wgt, labels + "_ZH_MET_");
    //mvp_objects(          pm, selection, processes, ops, wgt, labels + "_ZH_MET_");

    pm.Push<Hist1D>(Axis(6, 0, 6, "njet",    "N_{jet}",     {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ZH_MET_njet_");
    pm.Push<Hist1D>(Axis(6, 0, 6, "nbdfm",   "N_{b,med}",   {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ZH_MET_nbdfm_");
    pm.Push<Hist1D>(Axis(6, 0, 6, "nbdfl",   "N_{b,loose}", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ZH_MET_nbdfl_");

    pm.Push<Hist1D>(Axis(55, 90, 200, "met",     "p_{T}^{miss}", {}),          selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ZH_MET_MET_");
    pm.Push<Hist1D>(Axis(40,  0, 3.2, dphi_Hmet, "#Delta #phi(MET,H)", {2.5}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ZH_MET_dphi_Hmet");
    pm.Push<Hist1D>(Axis(50,  0,   1, max_Imini, "max. I_{mini}", {}),         selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ZH_MET_max_Imini");

    pm.Push<Hist1D>(Axis(50,    30,  150, j1_pt,  "p_{T}(j_{1}) [GeV]", {}), selection && "njet>0", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ZH_MET_j1_pt");
    pm.Push<Hist1D>(Axis(50,  -5.0,  5.0, j1_eta, "#eta(j_{1})",        {}), selection && "njet>0", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ZH_MET_j1_eta_");
    pm.Push<Hist1D>(Axis(63, -3.15, 3.15, j1_phi, "#phi(j_{1})",        {}), selection && "njet>0", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ZH_MET_j1_phi_");
    pm.Push<Hist1D>(Axis(50,     0,   40, j1_m,   "m(j_{1}) [GeV]",     {}), selection && "njet>0", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ZH_MET_j1_m");
  }

  void ZH_MET_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels, int nbins){
    sample_photon_plots(  pm, selection, processes, ops, wgt, nbins, labels + "_ZH_MET_");
    sample_lepton_plots(  pm, selection, processes, ops, wgt, nbins, labels + "_ZH_MET_");
    sample_llphoton_plots(pm, selection, processes, ops, wgt, nbins, labels + "_ZH_MET_");
    //sample_kinrefit_plots(pm, selection, processes, ops, wgt, nbins, labels + "_ZH_MET_");
    //mvp_objects(          pm, selection, processes, ops, wgt, labels + "_ZH_MET_");

    pm.Push<Hist1D>(Axis(6, 0, 6, "njet",    "N_{jet}",     {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ZH_MET_njet_");
    pm.Push<Hist1D>(Axis(6, 0, 6, "nbdfm",   "N_{b,med}",   {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ZH_MET_nbdfm_");
    pm.Push<Hist1D>(Axis(6, 0, 6, "nbdfl",   "N_{b,loose}", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ZH_MET_nbdfl_");

    pm.Push<Hist1D>(Axis(nbins, 90, 200, "met",     "p_{T}^{miss}", {}),          selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ZH_MET_MET_");
    pm.Push<Hist1D>(Axis(nbins,  0, 3.2, dphi_Hmet, "#Delta #phi(MET,H)", {2.5}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ZH_MET_dphi_Hmet");
    pm.Push<Hist1D>(Axis(nbins,  0,   1, max_Imini, "max. I_{mini}", {}),         selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ZH_MET_max_Imini");

    pm.Push<Hist1D>(Axis(nbins,    30,  150, j1_pt,  "p_{T}(j_{1}) [GeV]", {}), selection && "njet>0", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ZH_MET_j1_pt");
    pm.Push<Hist1D>(Axis(nbins,  -5.0,  5.0, j1_eta, "#eta(j_{1})",        {}), selection && "njet>0", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ZH_MET_j1_eta_");
    pm.Push<Hist1D>(Axis(nbins, -3.15, 3.15, j1_phi, "#phi(j_{1})",        {}), selection && "njet>0", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ZH_MET_j1_phi_");
    pm.Push<Hist1D>(Axis(nbins,     0,   40, j1_m,   "m(j_{1}) [GeV]",     {}), selection && "njet>0", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ZH_MET_j1_m");

  }

  //These are the control region plots made for the ggF category
  //Note mllg is not llphoton_m[0] because often control regions include photons that don't pass the photon object selections
  void ggF_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    sample_photon_plots(  pm, selection, processes, ops, wgt, labels + "_ggF_");
    sample_lepton_plots(  pm, selection, processes, ops, wgt, labels + "_ggF_");
    sample_llphoton_plots(pm, selection, processes, ops, wgt, labels + "_ggF_");
    //sample_kinrefit_plots(pm, selection, processes, ops, wgt, labels + "_ggF_");
    //mvp_objects(          pm, selection, processes, ops, wgt, labels + "_ggF_");

    //--kinematic angles--//
    pm.Push<Hist1D>(Axis(40,   -1  ,1,  lly_cosTheta, "cos(#Theta)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ggF_cTHETA_");
    pm.Push<Hist1D>(Axis(40,   -1,  1,  lly_costheta, "cos(#theta)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ggF_ctheta_");
    pm.Push<Hist1D>(Axis(40, -3.2,3.2,  lly_angphi,   "#phi",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ggF_phi_");

    pm.Push<Hist1D>(Axis(45,     0,   90, "met",  "p_{T}^{miss} [GeV]", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ggF_met");
    pm.Push<Hist1D>(Axis(50,    30,  150, j1_pt,  "p_{T}(j_{1}) [GeV]", {}), selection && "njet>0", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ggF_j1_pt");
    pm.Push<Hist1D>(Axis(50,  -5.0,  5.0, j1_eta, "#eta(j_{1})",        {}), selection && "njet>0", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ggF_j1_eta_");
    pm.Push<Hist1D>(Axis(50, -3.15, 3.15, j1_phi, "#phi(j_{1})",        {}), selection && "njet>0", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ggF_j1_phi_");
    pm.Push<Hist1D>(Axis(50,     0,   40, j1_m,   "m(j_{1}) [GeV]",     {}), selection && "njet>0", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ggF_j1_m");

  }

  void ggF_input_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    pm.Push<Hist1D>(Axis(45,0.0,1,      "photon_idmva[0]", "#gamma - IDMVA",            {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_IDMVA_");
    pm.Push<Hist1D>(Axis(35,0,3.5,      "photon_drmin[0]", "min( #Delta R(#gamma,l) )", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_minDR_");
    pm.Push<Hist1D>(Axis(60,0,6.0,      "photon_drmax[0]", "max( #Delta R(#gamma,l) )", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_maxDR_");
    pm.Push<Hist1D>(Axis(40,-2.6,2.6,   "photon_eta[0]",   "#eta(#gamma)",              {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_eta_");
    pm.Push<Hist1D>(Axis(40,0.01,0.25,  "photon_energyErr[0]/photon_pt[0]", "#sigma(#gamma)/E(#gamma)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_gamma_relunc_");
    pm.Push<Hist1D>(Axis(40,   -1  ,1,  lly_cosTheta, "cos(#Theta)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ggF_cTHETA_");
    pm.Push<Hist1D>(Axis(40,   -1,  1,  lly_costheta, "cos(#theta)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ggF_ctheta_");
    pm.Push<Hist1D>(Axis(40, -3.2,3.2,  lly_angphi,   "#phi",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ggF_phi_");
    pm.Push<Hist1D>(Axis(40,100,180,  mlly,        "m_{ll#gamma} [GeV]",     {120,130}), selection, processes ,ops).Weight(wgt).Tag("ShortName:" + labels + "_lly_m_");

    pm.Push<Hist1D>(Axis(40,0,2.5,      pTlly_mlly,  "p_{T}(ll#gamma)/m_{ll#gamma}", {1}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_lly_ptom_");

    pm.Push<Hist1D>(Axis(40,  -2.6,  2.6, l1_eta, "#eta(l_{1})",    {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_eta_l1_");
    pm.Push<Hist1D>(Axis(40,  -2.6,  2.6, l2_eta, "#eta(l_{2})",    {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_eta_l2_");
  }


  void VBF_input_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    pm.Push<Hist1D>(Axis(45,0.0,1,      "photon_idmva[0]", "#gamma - IDMVA",            {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_gamma_IDMVA_");
    pm.Push<Hist1D>(Axis(35,0,3.5,      "photon_drmin[0]", "min( #Delta R(#gamma,l) )", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_gamma_minDR_");
    pm.Push<Hist1D>(Axis(60,0,5.0,      "photon_drmax[0]", "max( #Delta R(#gamma,l) )", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_gamma_maxDR_");
    pm.Push<Hist1D>(Axis(40,-2.6,2.6,   "photon_eta[0]",   "#eta(#gamma)",              {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_gamma_eta_");
    pm.Push<Hist1D>(Axis(40,0.01,0.15,  "photon_energyErr[0]/photon_pt[0]", "#sigma(#gamma)/E(#gamma)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_gamma_relunc_");
    pm.Push<Hist1D>(Axis(40,   -1  ,1,  lly_cosTheta, "cos(#Theta)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_cTHETA_");
    pm.Push<Hist1D>(Axis(40,   -1,  1,  lly_costheta, "cos(#theta)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_ctheta_");
    pm.Push<Hist1D>(Axis(40, -3.2,3.2,  lly_angphi,   "#phi",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_phi_");
    pm.Push<Hist1D>(Axis(40,100,180,  mlly,        "m_{ll#gamma} [GeV]",     {120,130}), selection, processes ,ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_lly_m_");

    pm.Push<Hist1D>(Axis(40,0,2.5,      pTlly_mlly,  "p_{T}(ll#gamma)/m_{ll#gamma}", {1}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_lly_ptom_");

    pm.Push<Hist1D>(Axis(40,  -2.6,  2.6, l1_eta, "#eta(l_{1})",    {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_eta_l1_");
    pm.Push<Hist1D>(Axis(40,  -2.6,  2.6, l2_eta, "#eta(l_{2})",    {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_eta_l2_");
    pm.Push<Hist1D>(Axis(5, 0, 5, "njet", "N_{jet}", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_njet_");

    pm.Push<Hist1D>(Axis(40,    0, 1500, "dijet_m",    "m_{jj} [GeV]",            {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_dijet_m");
    //    pm.Push<Hist1D>(Axis(nbins,    0,  300, "dijet_pt",   "p_{T}(jj) [GeV]",         {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_dijet_pt");
    pm.Push<Hist1D>(Axis(40,    0,    9, "dijet_deta", "#Delta#eta(j_{1},j_{2})", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_deta");
    pm.Push<Hist1D>(Axis(40, 0,  3.15, diphi_dijet, "#Delta#phi(j_{1},j_{2}))", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_dphi");

    pm.Push<Hist1D>(Axis(40,0,3.15, dphi_lly_dijet,    "#Delta#phi(Z#gamma,jj))", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_jjlly_dphi");
    //    pm.Push<Hist1D>(Axis(40,    0, 4.5, dR_lly_dijet,      "#Delta R(Z#gamma,jj)",   {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_jjlly_dphi");
    pm.Push<Hist1D>(Axis(40,    0, 1.0, balance_lly_dijet, "System balance",         {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_sysbal");
    pm.Push<Hist1D>(Axis(40,    0,   6, zeppenfeld_y_jet,  "#gamma zeppenfeld",      {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_zep");
    //    pm.Push<Hist1D>(Axis(40,  0.4, 4.4, mindR_y_jet,       "Min. #DeltaR(#gamma,j)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_yjdr");

    pm.Push<Hist1D>(Axis(40,    30, 150, j1_pt,  "p_{T}(j_{1}) [GeV]", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_j1_pt");
    //    pm.Push<Hist1D>(Axis(nbins, -5.0,  5.0, j1_eta, "#eta(j_{1})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_j1_eta_");
    //    pm.Push<Hist1D>(Axis(nbins, -3.15, 3.15, j1_phi, "#phi(j_{1})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_j1_phi_");
    pm.Push<Hist1D>(Axis(40,0.4,6,      "photon_jet1_dr[0]", "#DeltaR(#gamma,j1)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_photon_j1_dr");
    pm.Push<Hist1D>(Axis(40,0.4,6,      "photon_jet2_dr[0]", "#DeltaR(#gamma,j2)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_photon_j2_dr");

    pm.Push<Hist1D>(Axis(40,   30,  150, j2_pt,  "p_{T}(j_{2}) [GeV]", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_j2_pt");
    //    pm.Push<Hist1D>(Axis(nbins, -5.0,  5.0, j2_eta, "#eta(j_{2})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_j2_eta_");
    //    pm.Push<Hist1D>(Axis(nbins, -3.15, 3.15, j2_phi, "#phi(j_{2})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_j2_phi_");

  }


  void ggF_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels, int nbins){
    sample_photon_plots(  pm, selection, processes, ops, wgt, nbins, labels + "_ggF_");
    sample_lepton_plots(  pm, selection, processes, ops, wgt, nbins, labels + "_ggF_");
    sample_llphoton_plots(pm, selection, processes, ops, wgt, nbins, labels + "_ggF_");
    //sample_kinrefit_plots(pm, selection, processes, ops, wgt, nbins, labels + "_ggF_");
    //mvp_objects(          pm, selection, processes, ops, wgt, labels + "_ggF_");

    //--kinematic angles--//
    pm.Push<Hist1D>(Axis(nbins,   -1  ,1,  lly_cosTheta, "cos(#Theta)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ggF_cTHETA_");
    pm.Push<Hist1D>(Axis(nbins,   -1,  1,  lly_costheta, "cos(#theta)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ggF_ctheta_");
    pm.Push<Hist1D>(Axis(nbins, -3.2,3.2,  lly_angphi,   "#phi",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ggF_phi_");

    pm.Push<Hist1D>(Axis(nbins,     0,   90, "met",  "p_{T}^{miss} [GeV]", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ggF_met");
    pm.Push<Hist1D>(Axis(nbins,    30,  150, j1_pt,  "p_{T}(j_{1}) [GeV]", {}), selection && "njet>0", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ggF_j1_pt");
    pm.Push<Hist1D>(Axis(nbins,  -5.0,  5.0, j1_eta, "#eta(j_{1})",        {}), selection && "njet>0", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ggF_j1_eta_");
    pm.Push<Hist1D>(Axis(nbins, -3.15, 3.15, j1_phi, "#phi(j_{1})",        {}), selection && "njet>0", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ggF_j1_phi_");
    pm.Push<Hist1D>(Axis(nbins,     0,   40, j1_m,   "m(j_{1}) [GeV]",     {}), selection && "njet>0", processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_ggF_j1_m");

  }

  //This function again adds control region plots for the VBF category based on variables used in Run 2. Some overlap with ggF region but removing variables specific to the "kinematic BDT"
  void VBF_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    sample_photon_plots(  pm, selection, processes, ops, wgt, labels + "_VBF_");
    sample_lepton_plots(  pm, selection, processes, ops, wgt, labels + "_VBF_");
    sample_llphoton_plots(pm, selection, processes, ops, wgt, labels + "_VBF_");
    //sample_kinrefit_plots(pm, selection, processes, ops, wgt, labels + "_VBF_");
    //mvp_objects(          pm, selection, processes, ops, wgt, labels + "_VBF_");

    pm.Push<Hist1D>(Axis(6, 0, 6, "njet", "N_{jet}", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_njet_");

    pm.Push<Hist1D>(Axis(40,    0, 1000, "dijet_m",    "m_{jj} [GeV]",            {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_dijet_m");
    pm.Push<Hist1D>(Axis(40,    0,  300, "dijet_pt",   "p_{T}(jj) [GeV]",         {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_dijet_pt");
    pm.Push<Hist1D>(Axis(90,    0,    9, "dijet_deta", "#Delta#eta(j_{1},j_{2})", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_deta");
    pm.Push<Hist1D>(Axis(64, -3.2,  3.2, "dijet_dphi", "#Delta#phi(j_{1},j_{2})", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_dphi");

    pm.Push<Hist1D>(Axis(30,-3.15,3.15, dphi_lly_dijet,    "#Delta#phi(Z#gamma,jj)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_jjlly_dphi");
    pm.Push<Hist1D>(Axis(30,    0, 4.5, dR_lly_dijet,      "#Delta R(Z#gamma,jj)",   {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_jjlly_dphi");
    pm.Push<Hist1D>(Axis(30,    0, 1.0, balance_lly_dijet, "System balance",         {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_sysbal");
    pm.Push<Hist1D>(Axis( 8,    0,   6, zeppenfeld_y_jet,  "#gamma zeppenfeld",      {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_zep");
    pm.Push<Hist1D>(Axis(15,  0.4, 4.4, mindR_y_jet,       "Min. #DeltaR(#gamma,j)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_yjdr");

    pm.Push<Hist1D>(Axis(60,    30,  250, j1_pt,  "p_{T}(j_{1}) [GeV]", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_j1_pt");
    pm.Push<Hist1D>(Axis(100, -5.0,  5.0, j1_eta, "#eta(j_{1})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_j1_eta_");
    pm.Push<Hist1D>(Axis(63, -3.15, 3.15, j1_phi, "#phi(j_{1})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_j1_phi_");

    pm.Push<Hist1D>(Axis(60,    30,  200, j2_pt,  "p_{T}(j_{2}) [GeV]", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_j2_pt");
    pm.Push<Hist1D>(Axis(100, -5.0,  5.0, j2_eta, "#eta(j_{2})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_j2_eta_");
    pm.Push<Hist1D>(Axis(63, -3.15, 3.15, j2_phi, "#phi(j_{2})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_j2_phi_");

    //pm.Push<Hist1D>(Axis(70,   0, 140, llgamma_pTt,        "p^{t}_{T} [GeV]",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_pTt");
  }

  void VBF_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels, int nbins){
    sample_photon_plots(  pm, selection, processes, ops, wgt, nbins, labels + "_VBF_");
    sample_lepton_plots(  pm, selection, processes, ops, wgt, nbins, labels + "_VBF_");
    sample_llphoton_plots(pm, selection, processes, ops, wgt, nbins, labels + "_VBF_");
    //sample_kinrefit_plots(pm, selection, processes, ops, wgt, nbins, labels + "_VBF_");
    //mvp_objects(          pm, selection, processes, ops, wgt, labels + "_VBF_");

    pm.Push<Hist1D>(Axis(6, 0, 6, "njet", "N_{jet}", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_njet_");

    pm.Push<Hist1D>(Axis(nbins,    0, 1500, "dijet_m",    "m_{jj} [GeV]",            {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_dijet_m");
    pm.Push<Hist1D>(Axis(nbins,    0,  300, "dijet_pt",   "p_{T}(jj) [GeV]",         {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_dijet_pt");
    pm.Push<Hist1D>(Axis(nbins,    0,    9, "dijet_deta", "#Delta#eta(j_{1},j_{2})", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_deta");
    pm.Push<Hist1D>(Axis(nbins, -3.2,  3.2, "dijet_dphi", "#Delta#phi(j_{1},j_{2})", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_dphi");

    pm.Push<Hist1D>(Axis(nbins,-3.15,3.15, dphi_lly_dijet,    "#Delta#phi(Z#gamma,jj)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_jjlly_dphi");
    pm.Push<Hist1D>(Axis(nbins,    0, 4.5, dR_lly_dijet,      "#Delta R(Z#gamma,jj)",   {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_jjlly_dphi");
    pm.Push<Hist1D>(Axis(nbins,    0, 1.0, balance_lly_dijet, "System balance",         {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_sysbal");
    pm.Push<Hist1D>(Axis(nbins,    0,   6, zeppenfeld_y_jet,  "#gamma zeppenfeld",      {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_zep");
    pm.Push<Hist1D>(Axis(nbins,  0.4, 4.4, mindR_y_jet,       "Min. #DeltaR(#gamma,j)", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_yjdr");

    pm.Push<Hist1D>(Axis(nbins,    30, 150, j1_pt,  "p_{T}(j_{1}) [GeV]", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_j1_pt");
    pm.Push<Hist1D>(Axis(nbins, -5.0,  5.0, j1_eta, "#eta(j_{1})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_j1_eta_");
    pm.Push<Hist1D>(Axis(nbins, -3.15, 3.15, j1_phi, "#phi(j_{1})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_j1_phi_");

    pm.Push<Hist1D>(Axis(nbins,   30,  150, j2_pt,  "p_{T}(j_{2}) [GeV]", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_j2_pt");
    pm.Push<Hist1D>(Axis(nbins, -5.0,  5.0, j2_eta, "#eta(j_{2})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_j2_eta_");
    pm.Push<Hist1D>(Axis(nbins, -3.15, 3.15, j2_phi, "#phi(j_{2})",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_VBF_j2_phi_");
  }


  //This function, as of 11/03/23 is a work in progress. No control region plots were made in HIG-19-014
  void Lep_controlregion_plots(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    sample_photon_plots(  pm, selection, processes, ops, wgt, labels + "_lep_");
    sample_lepton_plots(  pm, selection, processes, ops, wgt, labels + "_lep_");
    sample_llphoton_plots(pm, selection, processes, ops, wgt, labels + "_lep_");
  }

  void mvp_procs_and_corrections(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    pm.Push<Hist1D>(Axis(40, 0.8,  1.2, "w_lep",     "w_{lep}",     {}), selection, processes, ops).Weight(wgt).Tag("ShortName:MVP_procs_and_corrections_" + labels + "_w_lep");
    pm.Push<Hist1D>(Axis(40, 0.8,  1.2, "w_fs_lep",  "w_{fs lep}",  {}), selection, processes, ops).Weight(wgt).Tag("ShortName:MVP_procs_and_corrections_" + labels + "_w_fs_lep");
    pm.Push<Hist1D>(Axis(40, 0.8,  1.2, "w_bhig_df", "w_{bhig df}", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:MVP_procs_and_corrections_" + labels + "_w_bhig_df");
    pm.Push<Hist1D>(Axis(20, 0.8,  1.2, "w_isr",     "w_{isr}",     {}), selection, processes, ops).Weight(wgt).Tag("ShortName:MVP_procs_and_corrections_" + labels + "_w_isr");
    pm.Push<Hist1D>(Axis(40, 0.6,  1.4, "w_pu",      "w_{pu}",      {}), selection, processes, ops).Weight(wgt).Tag("ShortName:MVP_procs_and_corrections_" + labels + "_w_pu");
    pm.Push<Hist1D>(Axis(40, 0.8,  1.4, "w_prefire", "w_{prefire}", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:MVP_procs_and_corrections_" + labels + "_w_prefire");
    pm.Push<Hist1D>(Axis(40, 0.8,  1.2, "w_photon",  "w_{photon}",  {}), selection, processes, ops).Weight(wgt).Tag("ShortName:MVP_procs_and_corrections_" + labels + "_w_photon");
  }

  void mvp_objects(PlotMaker &pm, NamedFunc selection, std::vector<std::shared_ptr<Process>> &processes, std::vector<PlotOpt> &ops, NamedFunc wgt, std::string labels){
    pm.Push<Hist1D>(Axis(5, 0, 5, "nel",       "N_{e}",        {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_nel");
    pm.Push<Hist1D>(Axis(5, 0, 5, "nmu",       "N_{#mu}",      {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_nmu");
    pm.Push<Hist1D>(Axis(5, 0, 5, "nll",       "N_{ll}",       {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_nll");
    pm.Push<Hist1D>(Axis(7, 0, 7, "njet",      "N_{jet}",      {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_njet");
    pm.Push<Hist1D>(Axis(4, 0, 4, "nbdfl",     "N_{b,loose}",  {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_nbdfl");
    pm.Push<Hist1D>(Axis(4, 0, 4, "nbdfm",     "N_{b,medium}", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_nbdfm");
    pm.Push<Hist1D>(Axis(4, 0, 4, "nbdft",     "N_{b,tight}",  {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_nbdft");
    pm.Push<Hist1D>(Axis(2, 0, 2, "nphoton",   "N_{#gamma}",   {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_nphoton");
    pm.Push<Hist1D>(Axis(2, 0, 2, "nllphoton", "N_{ll#gamma}", {}), selection, processes, ops).Weight(wgt).Tag("ShortName:" + labels + "_nllphoton");
  }


  std::vector<std::vector<float>> plotting_bins(bool plotting_run3){
    std::vector<std::vector<float>> vec_nbins_run2 = {
    {5,   5, 10},
    {12, 12, 20},
    {12, 12, 20},
    {12, 12, 20},
    {50, 50, 75},
    {100,100,100}};

    std::vector<std::vector<float>> vec_nbins_run3 = {
    {5,   5,   5},
    {8,   8,  12},
    {8,   8,  12},
    {8,   8,  12},
    {50, 50,  75},
    {100,100,100}};

    if(!plotting_run3){ return vec_nbins_run2; 
    } else { return vec_nbins_run3; }
  }



//--------------------------------Skim Selector-----------------------------//
/*unordered_map<int,string> run2_skims = {{0, "Run2SkimttHlep"}, {1, "Run2SkimttHhad"}, {2, "Run2SkimZHMET"}, {3, "Run2SkimWH3l"}, {4, "Run2SkimVBF"},{5, "Run2SkimggF"}};
unordered_map<int,string> run3_skims = {{0, "Run3SkimttHlep"}, {1, "Run3SkimttHhad"}, {2, "Run3SkimZHMET"}, {3, "Run3SkimWH3l"}, {4, "Run3SkimVBF"},{5, "Run3SkimggF"}};

procs = ZgUtilities::procs_with_sig_scale("txt/samples_zgamma_skims.txt", run3_skims[category],10);
vector<shared_ptr<Process>> select_proper_skim_files(bool plotting_run3, string year, int category){

}
ZgUtilities::procs_sig_with_scale()
*/


}


