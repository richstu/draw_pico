#include "zgamma/zg_utilities.hpp"

#include <algorithm>
#include <stdlib.h>
#include <regex>
#include <vector>
#include "core/baby.hpp"
#include "core/palette.hpp"
#include "core/sample_loader.hpp"
#include "core/utilities.hpp"
#include "zgamma/zg_functions.hpp"
#include "zgamma/KinZfitter.h"

namespace ZgUtilities {
  using std::vector;
  using std::string;
  using std::to_string;
  using std::cout;
  using std::endl;
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
    std::cout << " CRASH LEP" << std::endl;
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
    std::cout << " CRASH SLEP" << std::endl;
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

  SampleLoader ZgSampleLoader() {
    SampleLoader zg_sample_loader;
    zg_sample_loader.LoadNamedFunc("HLT_pass_dilepton",ZgFunctions::HLT_pass_dilepton);
    zg_sample_loader.LoadNamedFunc("stitch",ZgFunctions::stitch);
    zg_sample_loader.LoadNamedFunc("HLT_pass_dilepton&&stitch",ZgFunctions::HLT_pass_dilepton&&ZgFunctions::stitch);
    zg_sample_loader.LoadNamedFunc("(HLT_pass_dilepton||HLT_pass_singlelepton)&&stitch",
        (ZgFunctions::HLT_pass_dilepton||ZgFunctions::HLT_pass_singlelepton)&&ZgFunctions::stitch);
    zg_sample_loader.LoadPalette("txt/colors_zgamma.txt","default");
    return zg_sample_loader;
  }
  
  void rename_signal(std::vector<std::shared_ptr<Process>> &procs, int factor) {
    for (std::shared_ptr<Process> &proc : procs) {
      if (proc->type_ == Process::Type::signal) {
        proc->name_ = proc->name_ + " (x" + std::to_string(factor) + ")";
      }
    }
  }

  std::vector<std::shared_ptr<Process>> procs_with_sig_scale(const std::string file_name, const std::string loader_name, int sig_factor){
    std::vector<std::shared_ptr<Process>> procs_loaded = ZgSampleLoader().LoadSamples(file_name, loader_name);
    rename_signal(procs_loaded,sig_factor);
    return procs_loaded;
  }


  //0 is untagged, 1 is loose, 2 is medium, 3 is tight
  int get_btag_wp_deepjet(int year, float discriminator_value) {
    if (abs(year)==2016) {
      if (discriminator_value > 0.7221) return 3;
      else if (discriminator_value > 0.3093) return 2;
      else if (discriminator_value > 0.0614) return 1;
      return 0;
    }
    else if (abs(year)==2017) {
      if (discriminator_value > 0.7489) return 3;
      else if (discriminator_value > 0.3033) return 2;
      else if (discriminator_value > 0.0521) return 1;
      return 0;
    }
    else if (abs(year)==2018) {
      if (discriminator_value > 0.7264) return 3;
      else if (discriminator_value > 0.2770) return 2;
      else if (discriminator_value > 0.0494) return 1;
      return 0;
    }
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
    bool isData = false;
    kinZfitter = new KinZfitter(isData);

    std::map<unsigned int, TLorentzVector> selectedLeptons;
    selectedLeptons[0] = AssignL1(b);
    selectedLeptons[1] = AssignL2(b);
    std::map<unsigned int, double> errorLeptons;
    errorLeptons[0] = AssignL1Error(b);
    errorLeptons[1] = AssignL2Error(b);
    std::map<unsigned int, TLorentzVector> selectedFsrMap = fsrphoton_ret(b);

    kinZfitter->Setup(selectedLeptons, selectedFsrMap, errorLeptons,b.ll_lepid()->at(0),b.SampleTypeString() );
    kinZfitter->KinRefitZ1();
    double massZ1REFIT = kinZfitter->GetRefitMZ1();
    delete kinZfitter;
    return massZ1REFIT;
  }

  double KinRefit(const Baby &b, TString txtFile) {
    KinZfitter *kinZfitter;
    bool isData = false;
    kinZfitter = new KinZfitter(isData,txtFile);

    std::map<unsigned int, TLorentzVector> selectedLeptons;
    selectedLeptons[0] = AssignL1(b);
    selectedLeptons[1] = AssignL2(b);
    std::map<unsigned int, double> errorLeptons;
    errorLeptons[0] = AssignL1Error(b);
    errorLeptons[1] = AssignL2Error(b);
    std::map<unsigned int, TLorentzVector> selectedFsrMap = fsrphoton_ret(b);

    kinZfitter->Setup(selectedLeptons, selectedFsrMap, errorLeptons,b.ll_lepid()->at(0),b.SampleTypeString()  );
    kinZfitter->KinRefitZ1();

    double massZ1REFIT = kinZfitter->GetRefitMZ1();
    delete kinZfitter;
    return massZ1REFIT;
  }

  std::vector<TLorentzVector> RefitP4(const Baby &b, TString txtFile) {
    KinZfitter *kinZfitter;
    bool isData = false;
    kinZfitter = new KinZfitter(isData,txtFile);

    std::map<unsigned int, TLorentzVector> selectedLeptons;
    selectedLeptons[0] = AssignL1(b);
    selectedLeptons[1] = AssignL2(b);
    std::map<unsigned int, double> errorLeptons;
    errorLeptons[0] = AssignL1Error(b);
    errorLeptons[1] = AssignL2Error(b);
    std::map<unsigned int, TLorentzVector> selectedFsrMap = fsrphoton_ret(b);

    kinZfitter->Setup(selectedLeptons, selectedFsrMap, errorLeptons,b.ll_lepid()->at(0) ,b.SampleTypeString() );
    kinZfitter->KinRefitZ1();
    std::vector<TLorentzVector> reFit = kinZfitter->GetRefitP4s();
    delete kinZfitter;
    return reFit;
  }

  double difference_check(const Baby &b, TString txtFile){
    KinZfitter *kinZfitter;
    bool isData = false;
    kinZfitter = new KinZfitter(isData,txtFile);

    std::map<unsigned int, TLorentzVector> selectedLeptons;
    selectedLeptons[0] = AssignL1(b);
    selectedLeptons[1] = AssignL2(b);
    std::map<unsigned int, double> errorLeptons;
    errorLeptons[0] = AssignL1Error(b);
    errorLeptons[1] = AssignL2Error(b);
    std::map<unsigned int, TLorentzVector> selectedFsrMap = fsrphoton_ret(b);

    kinZfitter->Setup(selectedLeptons, selectedFsrMap, errorLeptons,b.ll_lepid()->at(0) ,b.SampleTypeString() );
    kinZfitter->KinRefitZ1();
    std::vector<TLorentzVector> reFit = kinZfitter->GetRefitP4s();
    delete kinZfitter;
    return ((reFit[0]+reFit[1]).M() - (selectedLeptons[0]+selectedLeptons[1]).M());
  }

  double difference_check_lly(const Baby &b, TString txtFile){
    KinZfitter *kinZfitter;
    bool isData = false;
    kinZfitter = new KinZfitter(isData,txtFile);

    std::map<unsigned int, TLorentzVector> selectedLeptons;
    selectedLeptons[0] = AssignL1(b);
    selectedLeptons[1] = AssignL2(b);
    std::map<unsigned int, double> errorLeptons;
    errorLeptons[0] = AssignL1Error(b);
    errorLeptons[1] = AssignL2Error(b);
    std::map<unsigned int, TLorentzVector> selectedFsrMap = fsrphoton_ret(b);


    kinZfitter->Setup(selectedLeptons, selectedFsrMap, errorLeptons,b.ll_lepid()->at(0) ,b.SampleTypeString() );
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
    bool isData = false;
    kinZfitter = new KinZfitter(isData,txtFile);

    std::map<unsigned int, TLorentzVector> selectedLeptons;

    selectedLeptons[0] = AssignCorrL1(b);
    selectedLeptons[1] = AssignCorrL2(b);
    std::map<unsigned int, double> errorLeptons;
    errorLeptons[0] = AssignCorrL1Error(b);
    errorLeptons[1] = AssignCorrL2Error(b);
    std::map<unsigned int, TLorentzVector> selectedFsrMap = fsrphoton_ret(b);

    kinZfitter->Setup(selectedLeptons, selectedFsrMap, errorLeptons,b.ll_lepid()->at(0),b.SampleTypeString() );
    kinZfitter->KinRefitZ1();
    double massZ1REFIT = kinZfitter->GetRefitMZ1();
    delete kinZfitter;
    return massZ1REFIT;
  }

  std::vector<TLorentzVector> RefitP4Corr(const Baby &b, TString txtFile) {
    KinZfitter *kinZfitter;
    bool isData = false;
    kinZfitter = new KinZfitter(isData,txtFile);

    std::map<unsigned int, TLorentzVector> selectedLeptons;
    selectedLeptons[0] = AssignCorrL1(b);
    selectedLeptons[1] = AssignCorrL2(b);
    std::map<unsigned int, double> errorLeptons;
    errorLeptons[0] = AssignCorrL1Error(b);
    errorLeptons[1] = AssignCorrL2Error(b);
    std::map<unsigned int, TLorentzVector> selectedFsrMap = fsrphoton_ret(b);

    kinZfitter->Setup(selectedLeptons, selectedFsrMap, errorLeptons,b.ll_lepid()->at(0) ,b.SampleTypeString() );
    kinZfitter->KinRefitZ1();
    std::vector<TLorentzVector> reFit = kinZfitter->GetRefitP4s();
    delete kinZfitter;
    return reFit;
  }

}

