#include "zgamma/zg_utilities.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <regex>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>

#include "RooAbsPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooFFTConvPdf.h"
#include "RooGenericPdf.h"
#include "RooRealVar.h"

#include "core/baby.hpp"
#include "core/fastforest.hpp"
#include "core/named_func.hpp"
#include "core/named_func_utilities.hpp"
#include "core/mva_wrapper.hpp"
#include "core/palette.hpp"
#include "core/process.hpp"
#include "core/sample_loader.hpp"
#include "core/utilities.hpp"
#include "zgamma/zg_functions.hpp"
#include "zgamma/KinZfitter.hpp"

namespace ZgUtilities {
  using std::vector;
  using std::string;
  using std::map;
  using std::to_string;
  using std::cout;
  using std::endl;
  using std::shared_ptr;
  using std::function;
  using fastforest::FastForest;
  using NamedFuncUtilities::MapNamedFunc;
  // Returns negative lepton 4-momentum
  TLorentzVector AssignL1(const Baby &b, bool gen) {
    TLorentzVector l1;
    if(gen) {
      for(size_t i = 0; i < b.mc_id()->size(); i++)
        if(b.mc_id()->at(i) == 11 || 
           b.mc_id()->at(i) == 13
           //|| b.mc_id()->at(i) == 15
           ) {
            if(b.mc_mom()->at(i) == 23) {
              l1.SetPtEtaPhiM(b.mc_pt()->at(i),
                              b.mc_eta()->at(i),
                              b.mc_phi()->at(i),
                              b.mc_mass()->at(i));
              break;
            }
          }
    }
    else{
      int il(-1);
      if(b.ll_lepid()->at(0) == 11){
        if(b.el_charge()->at(b.ll_i1()->at(0)) < 0) il = b.ll_i1()->at(0);
        else                                        il = b.ll_i2()->at(0);
        l1.SetPtEtaPhiM(b.el_pt() ->at(il), 
                        b.el_eta()->at(il),
                        b.el_phi()->at(il), 0.00511);
      }
      else if(b.ll_lepid()->at(0) == 13){
        if(b.mu_charge()->at(b.ll_i1()->at(0)) < 0) il = b.ll_i1()->at(0);
        else                                        il = b.ll_i2()->at(0);
        l1.SetPtEtaPhiM(b.mu_pt() ->at(il),
                        b.mu_eta()->at(il),
                        b.mu_phi()->at(il), 0.105);
      }
    }
    return l1;
  }

/*  TLorentzVector AssignJJ(const Baby &b, bool gen){
    TLorentzVector jj;
    double  jj_m = 100000;
    int r_idx = -1;

    for(size_t i; i < b.dijet_m()->size(); i++){
      if( fabs( jj_m - 91.2 ) < fabs( (b.dijet_m() -> at(i)) - 91.2)){
        jj_m = b.dijet_m()->at(i);
        r_idx = i;
      }

    }
   
   jj.SetPtEtaPhiM(b.dijet_pt() ->at(r_idx),
                   b.dijet_eta()->at(r_idx),
                   b.dijet_phi()->at(r_idx),
   return jj;    
  }
*/
  // Returns positive lepton 4-momentum
  TLorentzVector AssignL2(const Baby &b, bool gen) {
    TLorentzVector l2;
    if(gen) {
      for(size_t i = 0; i < b.mc_id()->size(); i++)
        if(b.mc_id()->at(i) == -11 || 
           b.mc_id()->at(i) == -13
           //|| b.mc_id()->at(i) == -15
           ) {
            if(b.mc_mom()->at(i) == 23) {
              //cout<<"AssignL2 Z index: "<<b.mc_momidx()->at(i)<<endl;
              //cout<<"AssignL2 Z mom index: "<<b.mc_momidx()->at(b.mc_momidx()->at(i))<<endl;
              //cout<<"AssignL2 Z mom pid: "<<b.mc_id()->at(b.mc_momidx()->at(b.mc_momidx()->at(i)))<<endl;
              //cout<<"AssignL2 Z mom pt: "<<b.mc_pt()->at(b.mc_momidx()->at(b.mc_momidx()->at(i)))<<endl;
              l2.SetPtEtaPhiM(b.mc_pt()->at(i),
                              b.mc_eta()->at(i),
                              b.mc_phi()->at(i),
                              b.mc_mass()->at(i));
              break;
            }
          }
    }
    else{
      int il(-1);
      if(b.ll_lepid()->at(0) == 11){
        if(b.el_charge()->at(b.ll_i1()->at(0)) < 0) il = b.ll_i2()->at(0);
        else                                        il = b.ll_i1()->at(0);
        l2.SetPtEtaPhiM(b.el_pt() ->at(il), 
                        b.el_eta()->at(il),
                        b.el_phi()->at(il), 0.00511);
      }
      else if(b.ll_lepid()->at(0) == 13){
        if(b.mu_charge()->at(b.ll_i1()->at(0)) < 0) il = b.ll_i2()->at(0);
        else                                        il = b.ll_i1()->at(0);
        l2.SetPtEtaPhiM(b.mu_pt() ->at(il),
                        b.mu_eta()->at(il),
                        b.mu_phi()->at(il), 0.105);
      }
    }
    return l2;
  }

  double AssignL1Error(const Baby &b) {
    //std::cout << " CRASH LEP" << std::endl;
    double l1Err = 0.0;
    int il(-1);
    int flav = b.ll_lepid() ->at(0);

    TLorentzVector l1;
    if(flav == 11){
      if(b.el_charge()->at(b.ll_i1()->at(0)) < 0){ il = b.ll_i1()->at(0);}// lead = true; }
      else                                          { il = b.ll_i2()->at(0);}// lead = false; }

      l1.SetPtEtaPhiM(b.el_pt() ->at(il), b.el_eta()->at(il), b.el_phi()->at(il), 0.000511);
      l1Err = b.el_energyErr()->at(il) * l1.Pt() / l1.P();

    }
    else if(flav == 13){
      if(b.mu_charge()->at(b.ll_i1()->at(0)) < 0){ il = b.ll_i1()->at(0);}//lead=true;}
      else                                          { il = b.ll_i2()->at(0);}//lead=false;}

      l1Err = b.mu_ptErr()->at(il);
    }
    return l1Err;
  }

  // Returns positive lepton 4-momentum error
  double AssignL2Error(const Baby &b) {
    //std::cout << " CRASH SLEP" << std::endl;
    double l2Err = 0.0;
    int il(-1);
    int flav = b.ll_lepid() ->at(0);
    
    TLorentzVector l2;
    if(flav == 11){
      if(b.el_charge()->at(b.ll_i1()->at(0)) < 0) {il = b.ll_i2()->at(0);}//lead=false;}
      else                                           {il = b.ll_i1()->at(0);}//lead=true;}
      l2.SetPtEtaPhiM(b.el_pt() ->at(il), b.el_eta()->at(il), b.el_phi()->at(il), 0.000511);
      l2Err = b.el_energyErr()->at(il) * l2.Pt() / l2.P();
    }
    else if(flav == 13){
      if(b.mu_charge()->at(b.ll_i1()->at(0)) < 0) {il = b.ll_i2()->at(0);}//lead=false;}
      else                                           {il = b.ll_i1()->at(0);}//lead=true;}

      l2Err = b.mu_ptErr()->at(il);
    }
    return l2Err;
  }


  // Returns Z 4-momentum
  TLorentzVector AssignZ(const Baby &b, bool gen) {
    TLorentzVector ll;
    if(gen) 
        ll = AssignL1(b,gen) + AssignL2(b,gen);
    else
      ll.SetPtEtaPhiM(b.ll_pt() ->at(0),
                      b.ll_eta()->at(0),
                      b.ll_phi()->at(0),
                      b.ll_m()  ->at(0));
    return ll;
  }

  // Returns photon 4-momentum
  TLorentzVector AssignGamma(const Baby &b, bool gen) {
    TLorentzVector gamma;
    bool FoundGamma(false);
    if(gen) {
      for(size_t i = 0; i < b.mc_id()->size(); i++)
        if(b.mc_id()->at(i) == 22 && b.mc_pt()->at(i) > 5) {
          if (!FoundGamma) gamma.SetPtEtaPhiM(b.mc_pt()->at(i),
                                              b.mc_eta()->at(i),
                                              b.mc_phi()->at(i),
                                              b.mc_mass()->at(i));
          FoundGamma = true;
          if(b.mc_mom()->at(i) == 25) {
            //cout<<"Higgs pt:"<<b.mc_pt()->at(b.mc_momidx()->at(i))<<endl;
            //cout<<"AssignGamma H index: "<<b.mc_momidx()->at(i)<<endl;
            //cout<<"AssignGamma H pt: "<<b.mc_pt()->at(b.mc_momidx()->at(i))<<endl;
            gamma.SetPtEtaPhiM(b.mc_pt()->at(i),
                               b.mc_eta()->at(i),
                               b.mc_phi()->at(i),
                               b.mc_mass()->at(i));
            break;
          }
      }
    } else {
      gamma.SetPtEtaPhiM(b.photon_pt() ->at(0),
                         b.photon_eta()->at(0),
                         b.photon_phi()->at(0), 0);
    }
    return gamma;
  }

  TLorentzVector FindGamma(const Baby &b, bool gen) {
    TLorentzVector gamma;

    if(gen){
    			gamma.SetPtEtaPhiM(b.photon_pt()->at(0), b.photon_eta()->at(0), b.photon_phi()->at(0), 0);
    }

    if(!gen){
      for(unsigned int i = 0; i < b.photon_pt()->size(); i++) {
	  	  if(b.photon_pt() -> at(i) > 10 && fabs(b.photon_eta() -> at(i) < 2.6)){ 
         			gamma.SetPtEtaPhiM(b.photon_pt()->at(i), b.photon_eta()->at(i), b.photon_phi()->at(i), 0);
              break;

        }
	    }
      //gamma.SetPtEtaPhiM(0,0,0,0); -- See if this helps seg faults
    }
    return gamma;

  }
     
  double Findlly(const Baby &b) { 
    TLorentzVector photon,ll;
    if(b.nllphoton() > 0){
      return (b.llphoton_m() -> at(0));
    } else {
      photon = FindGamma(b);
      ll = AssignZ(b);
    }
      return (ll + photon).M();
  }



  // Returns Higgs 4-momentum
  TLorentzVector AssignH(const Baby &b, bool gen) {
    TLorentzVector h, y;
    if(gen) {
      bool FoundH(false);
      for(size_t i = 0; i < b.mc_id()->size(); i++) {
        if(b.mc_id()->at(i) == 25) {
          //cout<<"AssignH index:"<<i<<endl;
          //cout<<"AssignH pt:"<<b.mc_pt()->at(i)<<endl;
          FoundH = true;
          h.SetPtEtaPhiM(b.mc_pt()->at(i),
                         b.mc_eta()->at(i),
                         b.mc_phi()->at(i),
                         b.mc_mass()->at(i));
        }
      }
      if(!FoundH) h = AssignZ(b,gen) + AssignGamma(b,gen);
    }
    else {  
	 TLorentzVector z = AssignZ(b);
	 TLorentzVector photon;

   if(b.llphoton_m() -> size() > 0){
  	 //for(unsigned int i = 0; i < b.llphoton_pt()->size(); i++) {
	
        h.SetPtEtaPhiM(b.llphoton_pt()->at(0),
                       b.llphoton_eta()->at(0),
                       b.llphoton_phi()->at(0),
                       b.llphoton_m()->at(0));
      //}
   } else if ( (b.photon_pt() -> size() > 0) && (b.ll_m() -> size() > 0) && (b.llphoton_m() -> size() == 0) ){
        h = AssignZ(b,gen) + AssignGamma(b,gen);

   } else {
        h.SetPtEtaPhiM(0,0,0,0);
   }
  
  }
   return h;  
  }

  //
  // Variables used for defining kinematic angles presented in https://arxiv.org/pdf/1108.2274.pdf
  //

  // Returns 4-momentum of q1 (quark from gluon-gluon fusion)
  //  Defined in Equation 4
  TLorentzVector AssignQ1(const Baby &b, bool gen) {
    TLorentzVector h = AssignH(b,gen);
    TVector3 htran = h.BoostVector();
    htran.SetZ(0);
    h.Boost(-1*htran);
    TLorentzVector k1;
    double pz, E;
    pz = h.Pz() + h.E();
    E  = h.E()  + h.Pz();
    k1.SetPxPyPzE(0,0,pz/2,E/2);
    k1.Boost(htran);
    return k1;
  }


  // Returns 4-momentum of q2 (quark from gluon-gluon fusion)
  //  Defined in Equation 5
  TLorentzVector AssignQ2(const Baby &b, bool gen) {
    TLorentzVector k2;
    TLorentzVector h = AssignH(b,gen);
    TVector3 htran = h.BoostVector();
    htran.SetZ(0);
    h.Boost(-1*htran);
    double pz, E;
    pz = h.Pz() - h.E();
    E  = h.E()  - h.Pz();
    k2.SetPxPyPzE(0,0,pz/2,E/2);
    k2.Boost(htran);
    return k2;
  }
  
  // Returns magnitude of Z candidate 3-momentum 
  //  Defined in Equation 7
  double lambdaZ(const Baby &b, bool gen) {
    TLorentzVector P = AssignH(b, gen);
    TLorentzVector l1 = AssignL1(b,gen);
    TLorentzVector l2 = AssignL2(b,gen);
    TLorentzVector Z = AssignZ(b, gen);
    double M = P.M(), mll = Z.M();
    return sqrt(pow(P.Dot(Z)/M,2)-pow(mll,2));
  }

  // Cosine of angle between lepton 1 and parent Z in Higgs frame 
  //  Defined in Equation 13
  double cos_theta(const Baby &b, bool gen) {
    TLorentzVector P =  AssignH(b,gen);
    TLorentzVector l1 = AssignL1(b,gen);
    TLorentzVector l2 = AssignL2(b,gen);
    double M = P.M();
    double lZ = lambdaZ(b,gen);
    double ctheta = P.Dot(l1-l2)/(M*lZ);
    if(ctheta > 1) ctheta = 0.999;
    if(ctheta <-1) ctheta = -0.999;
    return ctheta;
  }

  // Cosine of angle between incoming quarks and outgoing Zs in higgs frame 
  //  Defined in Equation 8
  double cos_Theta(const Baby &b, bool gen) {
    TLorentzVector H  = AssignH(b,gen);
    TLorentzVector Z  = AssignZ(b,gen);
    TLorentzVector q1 = AssignQ1(b,gen);
    TLorentzVector q2 = AssignQ2(b,gen);
    //cout<<"a H: "<<H.Pt()<<" Z: "<<Z.Pt()<<endl;
    double M = H.M();
    double lZ = lambdaZ(b,gen);
    double cosTheta = Z.Dot(q1-q2)/(M*lZ);
    //cout<<"a q1: "<<q1.Pt()<<" q2: "<<q2.Pt()<<" M: "<<M<<" lZ:"<<lZ<<" cosT: "<<cosTheta<<endl;
    if(abs(cosTheta) > 1.01) cout << "ERROR: cTheta = " << cosTheta <<  endl;
    return cosTheta;
  }

  // Angle of the Z decay plane from the z-axis (defined in Equation 1) in the higgs frame
  //  Defined in Equation 21+22
  double Getphi(const Baby &b, bool gen) {
    TVector3 l1 = AssignL1(b, gen).Vect();
    TVector3 l2 = AssignL2(b, gen).Vect();
    TVector3 q1 = AssignQ1(b, gen).Vect();
    TVector3 Z  = AssignZ( b, gen).Vect();
    double cosphi, sinphi;
    cosphi = -1*l1.Cross(l2).Dot(q1.Cross(Z))/l1.Cross(l2).Mag()/q1.Cross(Z).Mag();
    sinphi = -1*l1.Cross(l2).Dot(q1)/l1.Cross(l2).Mag()/q1.Mag();
    double phi(0);
    if(abs(cosphi) > 1.01) cout << "ERROR: cphi = " << cosphi <<  endl;
    if(cosphi > 1) cosphi = 1;
    if(cosphi < -1) cosphi = -1;
    if(sinphi < 0) phi = -1*acos(cosphi);
    else           phi = acos(cosphi);
    return phi;
  }

  double pdrmax(const Baby &b){ 
    //TVector3 photon = AssignGamma(b).Vect();
    //TVector3 l1     = AssignL1(b).Vect();
    //TVector3 l2     = AssignL2(b).Vect();
    //return max(photon.DeltaR(l1),photon.DeltaR(l2));
    float l1eta = 0, l1phi = 0, l2eta = 0, l2phi = 0;
    if( b.ll_lepid() -> at(0) == 11){
      l1eta = b.el_eta() -> at( b.ll_i1() -> at(0) );
      l2eta = b.el_eta() -> at( b.ll_i2() -> at(0) );
      l1phi = b.el_phi() -> at( b.ll_i1() -> at(0) );
      l2phi = b.el_phi() -> at( b.ll_i2() -> at(0) );
    } else if( b.ll_lepid() -> at(0) == 13 ) {
      l1eta = b.mu_eta() -> at( b.ll_i1() -> at(0) );
      l2eta = b.mu_eta() -> at( b.ll_i2() -> at(0) );
      l1phi = b.mu_phi() -> at( b.ll_i1() -> at(0) );
      l2phi = b.mu_phi() -> at( b.ll_i2() -> at(0) );
    }
    float dr1 = deltaR(b.photon_eta()->at(0),b.photon_phi()->at(0),l1eta,l1phi);
    float dr2 = deltaR(b.photon_eta()->at(0),b.photon_phi()->at(0),l2eta,l2phi);
    return std::max(dr1, dr2);
  }

  // Returns 4-momentum of q1 (quark from gluon-gluon fusion)
  //  Defined in Equation 4
  TLorentzVector get_q1(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma) {
    TLorentzVector llg_p4 = lep_minus + lep_plus + gamma;
    TVector3 htran = llg_p4.BoostVector();
    htran.SetZ(0);
    llg_p4.Boost(-1*htran);
    double pz, E;
    pz = llg_p4.Pz() + llg_p4.E();
    E  = llg_p4.E()  + llg_p4.Pz();
    TLorentzVector k1;
    k1.SetPxPyPzE(0,0,pz/2,E/2);
    k1.Boost(htran);
    return k1;
  }
  
  // Returns 4-momentum of q2 (quark from gluon-gluon fusion)
  //  Defined in Equation 5
  TLorentzVector get_q2(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma) {
    TLorentzVector llg_p4 = lep_minus + lep_plus + gamma;
    TVector3 htran = llg_p4.BoostVector();
    htran.SetZ(0);
    llg_p4.Boost(-1*htran);
    double pz, E;
    pz = llg_p4.Pz() - llg_p4.E();
    E  = llg_p4.E()  - llg_p4.Pz();
    TLorentzVector k2;
    k2.SetPxPyPzE(0,0,pz/2,E/2);
    k2.Boost(htran);
    return k2;
  }
  
  // Returns magnitude of Z candidate 3-momentum 
  //  Defined in Equation 7
  double get_lambdaZ(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma) {
    TLorentzVector ll_p4 = lep_minus + lep_plus;
    TLorentzVector llg_p4 = ll_p4 + gamma;
    return sqrt(pow(llg_p4.Dot(ll_p4)/llg_p4.M(),2)-pow(ll_p4.M(),2));
  }

  // Cosine of angle between incoming quarks and outgoing Zs in higgs frame 
  //  Defined in Equation 8
  double get_cosTheta(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma) {
    TLorentzVector ll_p4 = lep_minus + lep_plus;
    TLorentzVector llg_p4 = ll_p4 + gamma;
    TLorentzVector q1_p4 = get_q1(lep_minus, lep_plus, gamma);
    TLorentzVector q2_p4 = get_q2(lep_minus, lep_plus, gamma);
    double lambdaZ = get_lambdaZ(lep_minus, lep_plus, gamma);
    double cosTheta = ll_p4.Dot(q2_p4-q1_p4)/(llg_p4.M()*lambdaZ);
    if(abs(cosTheta) > 1.01) std::cout << "ERROR: cTheta = " << cosTheta <<  std::endl;
    return cosTheta;
  }
  
  // Cosine of angle between lepton 1 and parent Z in Higgs frame 
  //  Defined in Equation 13
  double get_costheta(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma) {
    TLorentzVector llg_p4 = lep_minus + lep_plus + gamma;
    double lambdaZ = get_lambdaZ(lep_minus, lep_plus, gamma);
    double ctheta = llg_p4.Dot(lep_plus-lep_minus)/(llg_p4.M()*lambdaZ);
    if(ctheta > 1) ctheta = 0.999;
    if(ctheta <-1) ctheta = -0.999;
    return ctheta;
  }
  
  // Angle of the Z decay plane from the z-axis (defined in Equation 1) in the higgs frame
  //  Defined in Equation 21+22
  double get_phi(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma) {
    TLorentzVector llg_p4 = lep_minus + lep_plus + gamma;
    TLorentzVector q1_p4 = get_q1(lep_minus, lep_plus, gamma);
    TLorentzVector l1_p4 = lep_minus;
    TLorentzVector l2_p4 = lep_plus;
    TLorentzVector ll_p4 = lep_minus + lep_plus;
  
    // Boost l1, l2, q1, ll to llg frame
    TVector3 llgBoost = llg_p4.BoostVector();
    l1_p4.Boost(-1*llgBoost);
    l2_p4.Boost(-1*llgBoost);
    q1_p4.Boost(-1*llgBoost);
    ll_p4.Boost(-1*llgBoost);
  
    TVector3 l1_p3 = l1_p4.Vect();
    TVector3 l2_p3 = l2_p4.Vect();
    TVector3 q1_p3 = q1_p4.Vect();
    TVector3 ll_p3  = ll_p4.Vect();
  
    double sinTheta = sqrt(1-pow(get_cosTheta(lep_minus, lep_plus, gamma),2));
    double cosphi, sinphi;
    cosphi = -1*l1_p3.Cross(l2_p3).Dot(q1_p3.Cross(ll_p3))/l1_p3.Cross(l2_p3).Mag()/q1_p3.Cross(ll_p3).Mag();
    sinphi = -1*l1_p3.Cross(l2_p3).Dot(q1_p3)/l1_p3.Cross(l2_p3).Mag()/q1_p3.Mag()/sinTheta;
    double phi(0);
    if(abs(cosphi) > 1.01) std::cout << "ERROR: cphi = " << cosphi <<  std::endl;
    if(cosphi > 1) cosphi = 1;
    if(cosphi < -1) cosphi = -1;
    if(sinphi < 0) phi = -1*acos(cosphi);
    else           phi = acos(cosphi);
    if (phi < 0) phi += 2*TMath::Pi();
    return phi;
  }

  //returns working version of kinematic BDT
  std::shared_ptr<MVAWrapper> KinematicBdt() {
    std::shared_ptr<MVAWrapper> kin_bdt_reader = std::make_shared<MVAWrapper>("kinematic_bdt");
    kin_bdt_reader->SetVariable("photon_mva","photon_idmva[0]");
    kin_bdt_reader->SetVariable("min_dR","photon_drmin[0]");
    kin_bdt_reader->SetVariable("max_dR","photon_drmax[0]");
    kin_bdt_reader->SetVariable("pt_mass","llphoton_pt[0]/llphoton_m[0]");
    kin_bdt_reader->SetVariable("cosTheta","llphoton_cosTheta[0]");
    kin_bdt_reader->SetVariable("costheta","llphoton_costheta[0]");
    kin_bdt_reader->SetVariable("phi","llphoton_psi[0]");
    kin_bdt_reader->SetVariable("photon_res",ZgFunctions::photon_relpterr);
    kin_bdt_reader->SetVariable("photon_rapidity","photon_eta[0]");
    kin_bdt_reader->SetVariable("l1_rapidity",ZgFunctions::lead_lepton_eta);
    kin_bdt_reader->SetVariable("l2_rapidity",ZgFunctions::sublead_lepton_eta);
    kin_bdt_reader->BookMVA("/homes/oshiro/public_weights/shuffled_phidcomp_post_phidcomp_post_BDT.weights.xml");
    return kin_bdt_reader;
  }

  //returns NamedFunc that selects low BDT score category "ggF/untagged 4"
  NamedFunc category_ggh4_old(std::shared_ptr<MVAWrapper> kinematic_bdt) {
    NamedFunc kinematic_bdt_score = kinematic_bdt->GetDiscriminant();
    return (kinematic_bdt_score>-0.1&&kinematic_bdt_score<0.04);
  }

  //returns NamedFunc that selects medium BDT score category "ggF/untagged 3"
  NamedFunc category_ggh3_old(std::shared_ptr<MVAWrapper> kinematic_bdt) {
    NamedFunc kinematic_bdt_score = kinematic_bdt->GetDiscriminant();
    return (kinematic_bdt_score>0.04&&kinematic_bdt_score<0.16);
  }

  //returns NamedFunc that selects high BDT score category "ggF/untagged 2"
  NamedFunc category_ggh2_old(std::shared_ptr<MVAWrapper> kinematic_bdt) {
    NamedFunc kinematic_bdt_score = kinematic_bdt->GetDiscriminant();
    return (kinematic_bdt_score>0.16&&kinematic_bdt_score<0.26);
  }

  //returns NamedFunc that selects very high BDT score category "ggF/untagged 1"
  NamedFunc category_ggh1_old(std::shared_ptr<MVAWrapper> kinematic_bdt) {
    NamedFunc kinematic_bdt_score = kinematic_bdt->GetDiscriminant();
    return (kinematic_bdt_score>0.26);
  }

  const NamedFunc dijet_absdphi = MapNamedFunc("dijet_dphi", fabsf);
  const NamedFunc llphoton_dijet_absdphi = MapNamedFunc(
      "llphoton_dijet_dphi[0]", fabsf);

  //returns working version of dijet BDT
  vector<shared_ptr<MVAWrapper>> VbfBdts() {
    vector<shared_ptr<MVAWrapper>> vbf_bdt_readers = {
      std::make_shared<MVAWrapper>("vbf_bdt0"),
      std::make_shared<MVAWrapper>("vbf_bdt1"),
      std::make_shared<MVAWrapper>("vbf_bdt2"),
      std::make_shared<MVAWrapper>("vbf_bdt3")};
    for (std::shared_ptr<MVAWrapper> bdt_reader : vbf_bdt_readers) {
      bdt_reader->SetVariable("dijet_deta","dijet_deta");
      bdt_reader->SetVariable("dijet_dphi",dijet_absdphi);
      bdt_reader->SetVariable("j1_pt",ZgFunctions::lead_jet_pt);
      bdt_reader->SetVariable("j2_pt",ZgFunctions::sublead_jet_pt);
      bdt_reader->SetVariable("llphoton_dijet_balance",
                              "llphoton_dijet_balance[0]");
      bdt_reader->SetVariable("llphoton_dijet_dphi",llphoton_dijet_absdphi);
      bdt_reader->SetVariable("phi","llphoton_psi[0]");
      bdt_reader->SetVariable("min_dR","photon_drmin[0]");
      bdt_reader->SetVariable("max_dR","photon_drmax[0]");
      bdt_reader->SetVariable("costheta","llphoton_costheta[0]");
      bdt_reader->SetVariable("cosTheta","llphoton_cosTheta[0]");
      bdt_reader->SetVariable("pt_mass","llphoton_pt[0]/llphoton_m[0]");
      bdt_reader->SetVariable("l1_rapidity",ZgFunctions::lead_lepton_eta);
      bdt_reader->SetVariable("l2_rapidity",ZgFunctions::sublead_lepton_eta);
      bdt_reader->SetVariable("photon_rapidity","photon_eta[0]");
      bdt_reader->SetVariable("photon_mva","photon_idmva[0]");
      bdt_reader->SetVariable("photon_res",ZgFunctions::photon_relpterr);
      bdt_reader->SetVariable("photon_zeppenfeld","photon_zeppenfeld[0]");
      bdt_reader->SetVariable("photon_jet1_dr","photon_jet1_dr[0]");
      bdt_reader->SetVariable("photon_jet2_dr","photon_jet2_dr[0]");
    }
    vbf_bdt_readers[0]->BookMVA("/net/cms37/data1/rui/Training/dataset_2JClassic_pinnacles_run2p3/weights/TMVAClassification_BDT_0.weights.xml");
    vbf_bdt_readers[1]->BookMVA("/net/cms37/data1/rui/Training/dataset_2JClassic_pinnacles_run2p3/weights/TMVAClassification_BDT_1.weights.xml");
    vbf_bdt_readers[2]->BookMVA("/net/cms37/data1/rui/Training/dataset_2JClassic_pinnacles_run2p3/weights/TMVAClassification_BDT_2.weights.xml");
    vbf_bdt_readers[3]->BookMVA("/net/cms37/data1/rui/Training/dataset_2JClassic_pinnacles_run2p3/weights/TMVAClassification_BDT_3.weights.xml");
    return vbf_bdt_readers;
  }

  //returns NamedFunc that returns VBF score
  NamedFunc vbf_bdt_score(vector<shared_ptr<MVAWrapper>> vbf_bdts) {
    vector<NamedFunc> bdt_scores;
    for (shared_ptr<MVAWrapper> vbf_bdt : vbf_bdts) {
      bdt_scores.push_back(vbf_bdt->GetDiscriminant());
    }
    return NamedFunc("vbf_bdtscore",[bdt_scores](const Baby &b) 
        -> NamedFunc::ScalarType{
      if (b.njet()<2) return -1.0;
      int bdt_index = (((b.event()%314159)+1)%4);
      return bdt_scores[bdt_index].GetScalar(b);
    });
  }

  //returns NamedFunc that selects very high BDT score VBF category
  NamedFunc category_vbf1(vector<shared_ptr<MVAWrapper>> vbf_bdts) {
    vector<NamedFunc> bdt_scores;
    for (shared_ptr<MVAWrapper> vbf_bdt : vbf_bdts) {
      bdt_scores.push_back(vbf_bdt->GetDiscriminant());
    }
    return NamedFunc("vbf1",[bdt_scores](const Baby &b) 
        -> NamedFunc::ScalarType{
      if (b.njet()<2) return -1.0;
      int bdt_index = (((b.event()%314159)+1)%4);
      float score = bdt_scores[bdt_index].GetScalar(b);
      if (score > 0.489) return 1.0;
      return 0.0;
    });
  }

  //returns NamedFunc that selects high BDT score VBF category
  NamedFunc category_vbf2(vector<shared_ptr<MVAWrapper>> vbf_bdts) {
    vector<NamedFunc> bdt_scores;
    for (shared_ptr<MVAWrapper> vbf_bdt : vbf_bdts) {
      bdt_scores.push_back(vbf_bdt->GetDiscriminant());
    }
    return NamedFunc("vbf2",[bdt_scores](const Baby &b) 
        -> NamedFunc::ScalarType{
      if (b.njet()<2) return -1.0;
      int bdt_index = (((b.event()%314159)+1)%4);
      float score = bdt_scores[bdt_index].GetScalar(b);
      if (score > 0.286 && score < 0.489) return 1.0;
      return 0.0;
    });
  }

  //returns NamedFunc that selects medium BDT score VBF category
  NamedFunc category_vbf3(vector<shared_ptr<MVAWrapper>> vbf_bdts) {
    vector<NamedFunc> bdt_scores;
    for (shared_ptr<MVAWrapper> vbf_bdt : vbf_bdts) {
      bdt_scores.push_back(vbf_bdt->GetDiscriminant());
    }
    return NamedFunc("vbf3",[bdt_scores](const Baby &b) 
        -> NamedFunc::ScalarType{
      if (b.njet()<2) return -1.0;
      int bdt_index = (((b.event()%314159)+1)%4);
      float score = bdt_scores[bdt_index].GetScalar(b);
      if (score > 0.083 && score < 0.286) return 1.0;
      return 0.0;
    });
  }

  //returns NamedFunc that selects low BDT score VBF category
  NamedFunc category_vbf4(vector<shared_ptr<MVAWrapper>> vbf_bdts) {
    vector<NamedFunc> bdt_scores;
    for (shared_ptr<MVAWrapper> vbf_bdt : vbf_bdts) {
      bdt_scores.push_back(vbf_bdt->GetDiscriminant());
    }
    return NamedFunc("vbf4",[bdt_scores](const Baby &b) 
        -> NamedFunc::ScalarType{
      if (b.njet()<2) return -1.0;
      int bdt_index = (((b.event()%314159)+1)%4);
      float score = bdt_scores[bdt_index].GetScalar(b);
      if (score < 0.083) return 1.0;
      return 0.0;
    });
  }

  //returns XGBoost BDTs
  const vector<FastForest> XGBoostBDTs() {
    vector<string> xgb_features{"f0","f1","f2","f3","f4","f5","f6","f7","f8",
                                "f9","f10"};
    //in order vars are photon_mva, photon_res, max_dR, min_dR, pt_mass, 
    //cosTheta, costheta, phi, l1_rapidity, l2_rapidity, and photon_rapidity
    cout << "Loading XGBoost BDTs from "
         << "/homes/oshiro/public_weights/xgb_bdt_*.txt" << endl;
    cout << "Evaluator is NOT thread safe, ensure PlotMaker is using 1 thread"
         << endl;
    return {fastforest::load_txt(
        "/homes/oshiro/public_weights/xgb_bdt_0.txt", xgb_features),
        fastforest::load_txt(
        "/homes/oshiro/public_weights/xgb_bdt_1.txt", xgb_features),
        fastforest::load_txt(
        "/homes/oshiro/public_weights/xgb_bdt_2.txt", xgb_features),
        fastforest::load_txt(
        "/homes/oshiro/public_weights/xgb_bdt_3.txt", xgb_features)};
  }

  //hacky machinery to cache XGBoost score since it is slow
  bool compare_float_vectors_delta(vector<float> a, vector<float> b, 
                                   float delta=1.0e-5) {
    //all the BDT inputs are O(1) numbers
    if (a.size() != b.size()) return false;
    for (unsigned i = 0; i < a.size(); i++) {
      if (fabs(a[i]-b[i])>delta) return false;
    }
    return true;
  }

  //gets XGBoost score
  float GetXGBoostBDTScore(const vector<FastForest> &xgb_bdts, const Baby &b) {
    static vector<vector<float>> xgb_bdt_input_cache(4, vector<float>(11, 0));
    static vector<float> xgb_bdt_output_cache(4, 0.0);
    int bdt_idx = (b.event()%314159)%4;
    vector<float> bdt_inputs = {
        b.photon_idmva()->at(0),
        static_cast<float>(ZgFunctions::photon_relpterr.GetScalar(b)),
        b.photon_drmax()->at(0),
        b.photon_drmin()->at(0),
        b.llphoton_pt()->at(0)/b.llphoton_m()->at(0),
        b.llphoton_cosTheta()->at(0),
        b.llphoton_costheta()->at(0),
        b.llphoton_phi()->at(0),
        static_cast<float>(ZgFunctions::lead_lepton_eta.GetScalar(b)),
        static_cast<float>(ZgFunctions::sublead_lepton_eta.GetScalar(b)),
        b.photon_eta()->at(0)};
    if (!compare_float_vectors_delta(bdt_inputs, 
                                     xgb_bdt_input_cache[bdt_idx])) {
      for (unsigned iin = 0; iin < 11; iin++) {
        xgb_bdt_input_cache[bdt_idx][iin] = bdt_inputs[iin];
      }
      xgb_bdt_output_cache[bdt_idx] = xgb_bdts[bdt_idx](bdt_inputs.data())+0.5;
    }
    return xgb_bdt_output_cache[bdt_idx];
  }

  //Returns NamedFunc that returns XGBoost score
  NamedFunc XGBoostBDTScore(const vector<FastForest> &xgb_bdts) {
    return NamedFunc("xgb_bdtscore",[xgb_bdts](const Baby &b) 
        -> NamedFunc::ScalarType{
      return GetXGBoostBDTScore(xgb_bdts, b);
    });
  }

  //returns NamedFunc that selects low BDT score category "ggF 4"
  NamedFunc category_ggf4(const vector<FastForest> &xgb_bdts) {
    return NamedFunc("ggf4",[xgb_bdts](const Baby &b) 
        -> NamedFunc::ScalarType{
      float score = GetXGBoostBDTScore(xgb_bdts, b);
      if (score<0.47) return 1.0;
      return 0.0;
    });
  }

  //returns NamedFunc that selects medium BDT score category "ggF 3"
  NamedFunc category_ggf3(const vector<FastForest> &xgb_bdts) {
    return NamedFunc("ggf3",[xgb_bdts](const Baby &b) 
        -> NamedFunc::ScalarType{
      float score = GetXGBoostBDTScore(xgb_bdts, b);
      if (score>0.47&&score<0.64) return 1.0;
      return 0.0;
    });
  }

  //returns NamedFunc that selects high BDT score category "ggF 2"
  NamedFunc category_ggf2(const vector<FastForest> &xgb_bdts) {
    return NamedFunc("ggf2",[xgb_bdts](const Baby &b) 
        -> NamedFunc::ScalarType{
      float score = GetXGBoostBDTScore(xgb_bdts, b);
      if (score>0.64&&score<0.81) return 1.0;
      return 0.0;
    });
  }

  //returns NamedFunc that selects very high BDT score category "ggF 1"
  NamedFunc category_ggf1(const vector<FastForest> &xgb_bdts) {
    return NamedFunc("ggf1",[xgb_bdts](const Baby &b) 
        -> NamedFunc::ScalarType{
      float score = GetXGBoostBDTScore(xgb_bdts, b);
      if (score>0.81) return 1.0;
      return 0.0;
    });
  }

  //returns a sample loader that has the H->Zy colors pre-sets and NamedFuncs loaded
  SampleLoader ZgSampleLoader() {
    SampleLoader zg_sample_loader;
    zg_sample_loader.LoadNamedFunc("HLT_pass_dilepton",ZgFunctions::HLT_pass_dilepton);
    zg_sample_loader.LoadNamedFunc("stitch",ZgFunctions::stitch);
    zg_sample_loader.LoadNamedFunc("HLT_pass_dilepton&&stitch",ZgFunctions::HLT_pass_dilepton&&ZgFunctions::stitch);
    zg_sample_loader.LoadNamedFunc("(HLT_pass_dilepton||HLT_pass_singlelepton)&&stitch",
        (ZgFunctions::HLT_pass_dilepton||ZgFunctions::HLT_pass_singlelepton)&&ZgFunctions::stitch);
    zg_sample_loader.LoadNamedFunc("hzg_elchannel",ZgFunctions::Ztoee);
    zg_sample_loader.LoadNamedFunc("hzg_muchannel",ZgFunctions::ZtoMuMu);
    zg_sample_loader.LoadNamedFunc("use_event&&(trig_single_el||trig_single_mu"
        "||trig_double_el||trig_double_mu)&&photon_isjet",
        ("use_event&&(trig_single_el||trig_single_mu||trig_double_el||"
        "trig_double_mu)")&&ZgFunctions::photon_isjet);
    zg_sample_loader.LoadNamedFunc("use_event&&(trig_single_el||trig_single_mu"
        "||trig_double_el||trig_double_mu)&&!photon_isjet",
        ("use_event&&(trig_single_el||trig_single_mu||trig_double_el||"
        "trig_double_mu)")&&!ZgFunctions::photon_isjet);
    //zg_sample_loader.LoadNamedFunc("use_event&&trig&&photon_isjet",
    //    "use_event"&&ZgFunction::trig&&ZgFunctions::photon_isjet);
    //zg_sample_loader.LoadNamedFunc("use_event&&trig&&photon_isother",
    //    "use_event"&&ZgFunction::trig&&ZgFunctions::photon_isother);
    //zg_sample_loader.LoadNamedFunc("use_event&&trig&&photon_isisr",
    //    "use_event"&&ZgFunction::trig&&ZgFunctions::photon_isisr);
    //zg_sample_loader.LoadNamedFunc("use_event&&trig&&photon_isfsr",
    //    "use_event"&&ZgFunction::trig&&ZgFunctions::photon_isfsr);
    zg_sample_loader.LoadPalette("txt/colors_zgamma_official.txt","default");
    return zg_sample_loader;
  }
  
  void rename_signal(std::vector<std::shared_ptr<Process>> &procs, int factor) {
    for (std::shared_ptr<Process> &proc : procs) {
      if (proc->type_ == Process::Type::signal && factor != 1) {
        proc->name_ = proc->name_ + " (x" + std::to_string(factor) + ")";
      }
    }
  }

  std::vector<std::shared_ptr<Process>> procs_with_sig_scale(const std::string file_name, const std::string loader_name, int sig_factor){
    std::vector<std::shared_ptr<Process>> procs_loaded = ZgSampleLoader().LoadSamples(file_name, loader_name);
    rename_signal(procs_loaded,sig_factor);
    return procs_loaded;
  }


  const std::map<std::string, std::vector<float>> btag_df_wpts{
      {"2016APV", std::vector<float>({0.0508, 0.2598, 0.6502})},
      {"2016", std::vector<float>({0.0480, 0.2489, 0.6377})},
      {"2017", std::vector<float>({0.0532, 0.3040, 0.7476})},
      {"2018", std::vector<float>({0.0490, 0.2783, 0.7100})},
      {"2022", std::vector<float>({0.0583, 0.3086, 0.7183})},
      {"2022EE", std::vector<float>({0.0614, 0.3196, 0.73})},
      {"2023", std::vector<float>({0.0479, 0.2431, 0.6553})},
      {"2023BPix", std::vector<float>({0.048, 0.2435, 0.6563})}
      };

  //returns WP. 1 is loose, 2 is medium, 3 is tight
  float get_btag_wp_deepjet(const std::string& year, int wp) {
    if (wp<1 || wp>3)
      return 0;
    std::string year_clean(year);
    if (year_clean[0]=='-') {
      year_clean = year_clean.substr(1, string::npos);
    }
    if (btag_df_wpts.count(year_clean)>0)
      return btag_df_wpts.at(year_clean).at(wp-1);
    return 0;
  }

  double dR(double eta1, double eta2,double phi1,double phi2){
    return ((eta1-eta2)*(eta1-eta2) + (phi1-phi2)*(phi1-phi2));
  }

  std::map<unsigned int, TLorentzVector> fsrphoton_ret(const Baby &b){
    TLorentzVector fsr,fsr2, l1, l2;
    std::map<unsigned int, TLorentzVector> return_map;
    fsr2.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
    fsr.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
    if(b.ll_lepid() -> at(0) == 11){
      //return_map[0] = fsr; return_map[1] = fsr2;
      return return_map;
    }
    vector<int> fsridces= {-1,-1};
    for(int idx_fsr=0; idx_fsr < (b.nfsrphoton()); idx_fsr++){//was +999 for before lassen_v0
      if(dR(b.fsrphoton_eta() ->at(idx_fsr),b.photon_eta() -> at(0), b.fsrphoton_phi() ->at(idx_fsr),b.photon_phi() -> at(0)) < 0.2 ){continue;}

      if( dR(b.fsrphoton_eta() ->at(idx_fsr),b.mu_eta() -> at(b.ll_i1()->at(0)), b.fsrphoton_phi() ->at(idx_fsr),b.mu_phi() -> at(b.ll_i1()->at(0)) < 0.4 ) ){
        if(fsridces[0]>-1 && b.fsrphoton_droveret2() -> at(idx_fsr) > b.fsrphoton_droveret2() -> at(fsridces[0])){continue;}
        fsridces[0] = (idx_fsr);
        continue;
      }
      if( dR(b.fsrphoton_eta() ->at(idx_fsr),b.mu_eta() -> at(b.ll_i2()->at(0)), b.fsrphoton_phi() ->at(idx_fsr),b.mu_phi() -> at(b.ll_i2()->at(0)) < 0.4 ) ){
        if(fsridces[1]>-1 && b.fsrphoton_droveret2() -> at(idx_fsr) > b.fsrphoton_droveret2() -> at(fsridces[1])){continue;}
        fsridces[1] = (idx_fsr);
        continue;
      }

    }

    if(fsridces[0]>-1){
      fsr.SetPtEtaPhiM( b.fsrphoton_pt() -> at(fsridces[0]),b.fsrphoton_eta() -> at(fsridces[0]),b.fsrphoton_phi() -> at(fsridces[0]),0 );
    }
    if(fsridces[1]>-1){
      fsr2.SetPtEtaPhiM( b.fsrphoton_pt() -> at(fsridces[1]),b.fsrphoton_eta() -> at(fsridces[1]),b.fsrphoton_phi() -> at(fsridces[1]),0 );
    }
    return_map[0] = fsr; return_map[1] = fsr2;
    return return_map;
  }

  double KinRefit(const Baby &b) {
    KinZfitter *kinZfitter;
    kinZfitter = new KinZfitter();

    std::map<unsigned int, TLorentzVector> selectedLeptons;
    selectedLeptons[0] = AssignL1(b);
    selectedLeptons[1] = AssignL2(b);
    std::map<unsigned int, double> errorLeptons;
    errorLeptons[0] = AssignL1Error(b);
    errorLeptons[1] = AssignL2Error(b);
    std::map<unsigned int, TLorentzVector> selectedFsrMap = fsrphoton_ret(b);

    kinZfitter->Setup(selectedLeptons, selectedFsrMap, errorLeptons);
    kinZfitter->KinRefitZ1();
    double massZ1REFIT = kinZfitter->GetRefitMZ1();
    delete kinZfitter;
    return massZ1REFIT;
  }

  double KinRefit(const Baby &b, TString txtFile) {
    KinZfitter *kinZfitter;
    kinZfitter = new KinZfitter(txtFile);

    std::map<unsigned int, TLorentzVector> selectedLeptons;
    selectedLeptons[0] = AssignL1(b);
    selectedLeptons[1] = AssignL2(b);
    std::map<unsigned int, double> errorLeptons;
    errorLeptons[0] = AssignL1Error(b);
    errorLeptons[1] = AssignL2Error(b);
    std::map<unsigned int, TLorentzVector> selectedFsrMap = fsrphoton_ret(b);

    kinZfitter->Setup(selectedLeptons, selectedFsrMap, errorLeptons);
    kinZfitter->KinRefitZ1();

    double massZ1REFIT = kinZfitter->GetRefitMZ1();
    delete kinZfitter;
    return massZ1REFIT;
  }

  std::vector<TLorentzVector> RefitP4(const Baby &b, TString txtFile) {
    KinZfitter *kinZfitter;
    kinZfitter = new KinZfitter(txtFile);

    std::map<unsigned int, TLorentzVector> selectedLeptons;
    selectedLeptons[0] = AssignL1(b);
    selectedLeptons[1] = AssignL2(b);
    std::map<unsigned int, double> errorLeptons;
    errorLeptons[0] = AssignL1Error(b);
    errorLeptons[1] = AssignL2Error(b);
    std::map<unsigned int, TLorentzVector> selectedFsrMap = fsrphoton_ret(b);

    kinZfitter->Setup(selectedLeptons, selectedFsrMap, errorLeptons);
    kinZfitter->KinRefitZ1();
    std::vector<TLorentzVector> reFit = kinZfitter->GetRefitP4s();
    delete kinZfitter;
    return reFit;
  }

  double difference_check(const Baby &b, TString txtFile){
    KinZfitter *kinZfitter;
    kinZfitter = new KinZfitter(txtFile);

    std::map<unsigned int, TLorentzVector> selectedLeptons;
    selectedLeptons[0] = AssignL1(b);
    selectedLeptons[1] = AssignL2(b);
    std::map<unsigned int, double> errorLeptons;
    errorLeptons[0] = AssignL1Error(b);
    errorLeptons[1] = AssignL2Error(b);
    std::map<unsigned int, TLorentzVector> selectedFsrMap = fsrphoton_ret(b);

    kinZfitter->Setup(selectedLeptons, selectedFsrMap, errorLeptons);
    kinZfitter->KinRefitZ1();
    std::vector<TLorentzVector> reFit = kinZfitter->GetRefitP4s();
    delete kinZfitter;
    return ((reFit[0]+reFit[1]).M() - (selectedLeptons[0]+selectedLeptons[1]).M());
  }

  double difference_check_lly(const Baby &b, TString txtFile){
    KinZfitter *kinZfitter;
    kinZfitter = new KinZfitter(txtFile);

    std::map<unsigned int, TLorentzVector> selectedLeptons;
    selectedLeptons[0] = AssignL1(b);
    selectedLeptons[1] = AssignL2(b);
    std::map<unsigned int, double> errorLeptons;
    errorLeptons[0] = AssignL1Error(b);
    errorLeptons[1] = AssignL2Error(b);
    std::map<unsigned int, TLorentzVector> selectedFsrMap = fsrphoton_ret(b);


    kinZfitter->Setup(selectedLeptons, selectedFsrMap, errorLeptons);
    kinZfitter->KinRefitZ1();
    std::vector<TLorentzVector> reFit = kinZfitter->GetRefitP4s();
    TLorentzVector ph = AssignGamma(b);
    delete kinZfitter;
    return ((reFit[0]+reFit[1]+ph).M() - (selectedLeptons[0]+selectedLeptons[1]+ph).M());
  }

  //This function is used to assign the corrected muon pt for negative lepton
  TLorentzVector AssignCorrL1(const Baby &b) {
    TLorentzVector l1;
    int il(-1);
    if(b.ll_lepid()->at(0) == 11){
        if(b.el_charge()->at(b.ll_i1()->at(0)) < 0) il = b.ll_i1()->at(0);
        else                                        il = b.ll_i2()->at(0);
        l1.SetPtEtaPhiM(b.el_pt() ->at(il),
                        b.el_eta()->at(il),
                        b.el_phi()->at(il), 0.00511);
    }
    else if(b.ll_lepid()->at(0) == 13){
        if(b.mu_charge()->at(b.ll_i1()->at(0)) < 0) il = b.ll_i1()->at(0);
        else                                        il = b.ll_i2()->at(0);
        l1.SetPtEtaPhiM(b.mu_corrected_pt() ->at(il),
                        b.mu_eta()->at(il),
                        b.mu_phi()->at(il), 0.105);
    }
    return l1;
  }


  // Returns positive corrected lepton 4-momentum
  TLorentzVector AssignCorrL2(const Baby &b) {
    TLorentzVector l2;
    int il(-1);
    if(b.ll_lepid()->at(0) == 11){
      if(b.el_charge()->at(b.ll_i1()->at(0)) < 0) il = b.ll_i2()->at(0);
      else                                        il = b.ll_i1()->at(0);
        l2.SetPtEtaPhiM(b.el_pt() ->at(il),
                        b.el_eta()->at(il),
                        b.el_phi()->at(il), 0.00511);
    }
    else if(b.ll_lepid()->at(0) == 13){
      if(b.mu_charge()->at(b.ll_i1()->at(0)) < 0) il = b.ll_i2()->at(0);
      else                                        il = b.ll_i1()->at(0);
        l2.SetPtEtaPhiM(b.mu_corrected_pt() ->at(il),
                        b.mu_eta()->at(il),
                        b.mu_phi()->at(il), 0.105);
    }
    return l2;
  }

  double AssignCorrL1Error(const Baby &b) {
    double l1Err = 0.0;
    int il(-1);
    int flav = b.ll_lepid() ->at(0);

    TLorentzVector l1;
    if(flav == 11){
      if(b.el_charge()->at(b.ll_i1()->at(0)) < 0){ il = b.ll_i1()->at(0);}// lead = true; }
      else                                          { il = b.ll_i2()->at(0);}// lead = false; }

      l1.SetPtEtaPhiM(b.el_pt() ->at(il), b.el_eta()->at(il), b.el_phi()->at(il), 0.000511);
      l1Err = b.el_energyErr()->at(il) * l1.Pt() / l1.P();

    }
    else if(flav == 13){
      if(b.mu_charge()->at(b.ll_i1()->at(0)) < 0){ il = b.ll_i1()->at(0);}//lead=true;}
      else                                          { il = b.ll_i2()->at(0);}//lead=false;}

      l1Err = b.mu_corrected_ptErr()->at(il);
    }
    return l1Err;
  }

  // Returns positive lepton 4-momentum error
  double AssignCorrL2Error(const Baby &b) {
    double l2Err = 0.0;
    int il(-1);
    int flav = b.ll_lepid() ->at(0);

    TLorentzVector l2;
    if(flav == 11){
      if(b.el_charge()->at(b.ll_i1()->at(0)) < 0) {il = b.ll_i2()->at(0);}//lead=false;}
      else                                           {il = b.ll_i1()->at(0);}//lead=true;}
      l2.SetPtEtaPhiM(b.el_pt() ->at(il), b.el_eta()->at(il), b.el_phi()->at(il), 0.000511);
      l2Err = b.el_energyErr()->at(il) * l2.Pt() / l2.P();
    }
    else if(flav == 13){
      if(b.mu_charge()->at(b.ll_i1()->at(0)) < 0) {il = b.ll_i2()->at(0);}//lead=false;}
      else                                           {il = b.ll_i1()->at(0);}//lead=true;}

      //eta_mu = b.mu_eta() -> at(il);
      l2Err = b.mu_corrected_ptErr()->at(il);
    }
    return l2Err;
  }

  double KinRefitCorrected(const Baby &b, TString txtFile) {
    KinZfitter *kinZfitter;
    kinZfitter = new KinZfitter(txtFile);

    std::map<unsigned int, TLorentzVector> selectedLeptons;

    selectedLeptons[0] = AssignCorrL1(b);
    selectedLeptons[1] = AssignCorrL2(b);
    std::map<unsigned int, double> errorLeptons;
    errorLeptons[0] = AssignCorrL1Error(b);
    errorLeptons[1] = AssignCorrL2Error(b);
    std::map<unsigned int, TLorentzVector> selectedFsrMap = fsrphoton_ret(b);

    kinZfitter->Setup(selectedLeptons, selectedFsrMap, errorLeptons);
    kinZfitter->KinRefitZ1();
    double massZ1REFIT = kinZfitter->GetRefitMZ1();
    delete kinZfitter;
    return massZ1REFIT;
  }

  std::vector<TLorentzVector> RefitP4Corr(const Baby &b, TString txtFile) {
    KinZfitter *kinZfitter;
    kinZfitter = new KinZfitter(txtFile);

    std::map<unsigned int, TLorentzVector> selectedLeptons;
    selectedLeptons[0] = AssignCorrL1(b);
    selectedLeptons[1] = AssignCorrL2(b);
    std::map<unsigned int, double> errorLeptons;
    errorLeptons[0] = AssignCorrL1Error(b);
    errorLeptons[1] = AssignCorrL2Error(b);
    std::map<unsigned int, TLorentzVector> selectedFsrMap = fsrphoton_ret(b);

    kinZfitter->Setup(selectedLeptons, selectedFsrMap, errorLeptons);
    kinZfitter->KinRefitZ1();
    std::vector<TLorentzVector> reFit = kinZfitter->GetRefitP4s();
    delete kinZfitter;
    return reFit;
  }

  //sets all processes to background, useful for making colz 2D plots of data
  void SetProcessesBackground(std::vector<std::shared_ptr<Process>> &processes) {
    for (std::shared_ptr<Process> process : processes) {
      process->type_ = Process::Type::background;
    }
  }

  //returns lepton pt invariant mass with custom refit
  vector<double> get_lep_pt_custom_refit(const Baby &b, 
      shared_ptr<KinZfitter> kinZfitter, NamedFunc el_pt_var, 
      NamedFunc mu_pt_var) {

    std::map<unsigned int, TLorentzVector> selectedLeptons;
    selectedLeptons[0] = TLorentzVector();
    selectedLeptons[1] = TLorentzVector();
    int l1_idx = b.ll_i1()->at(0);
    int l2_idx = b.ll_i2()->at(0);
    if (b.ll_lepid()->at(0) == 11) {
      vector<double> el_pt = el_pt_var.GetVector(b);
      selectedLeptons[0].SetPtEtaPhiM(el_pt[l1_idx], b.el_eta()->at(l1_idx),
                                      b.el_phi()->at(l1_idx), 0.000511);
      selectedLeptons[1].SetPtEtaPhiM(el_pt[l2_idx], b.el_eta()->at(l2_idx),
                                      b.el_phi()->at(l2_idx), 0.000511);
    }
    else {
      vector<double> mu_pt = mu_pt_var.GetVector(b);
      selectedLeptons[0].SetPtEtaPhiM(mu_pt[l1_idx], b.mu_eta()->at(l1_idx),
                                      b.mu_phi()->at(l1_idx), 0.105);
      selectedLeptons[1].SetPtEtaPhiM(mu_pt[l2_idx], b.mu_eta()->at(l2_idx),
                                      b.mu_phi()->at(l2_idx), 0.105);
    }
    map<unsigned int, double> errorLeptons;
    errorLeptons[0] = AssignL1Error(b);
    errorLeptons[1] = AssignL2Error(b);
    map<unsigned int, TLorentzVector> selectedFsrMap = fsrphoton_ret(b);

    kinZfitter->Setup(selectedLeptons, selectedFsrMap, errorLeptons);
    kinZfitter->KinRefitZ1();
    vector<TLorentzVector> reFit = kinZfitter->GetRefitP4s();
    vector<double> refit_pt = {reFit[0].Pt(), reFit[1].Pt()};
    return refit_pt;
  }

}

