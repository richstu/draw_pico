#include <vector>

#include "zgamma/apply_zg_trigeffs.hpp"


std::vector<float> eff_elhltpt12(float pt, float abseta) {
  float errup=0.0, errdown=0.0, eff=1.0;
  if (abseta > 0 && abseta < 0.8) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 10 && pt < 15) {eff = 1; errup = 0; errdown = 0.108691;}
    else if (pt > 15 && pt < 20){eff = 1; errup = 0; errdown = 0.0100098;} 
    else if (pt > 20 && pt < 25) {eff = 1; errup = 0; errdown = 0.00443704;}
    else if (pt > 25 && pt < 27) {eff = 1; errup = 0; errdown = 0.00354779;}
    else if (pt > 27 && pt < 30) {eff = 1; errup = 0; errdown = 7.46625e-05;}
    else if (pt > 30 && pt < 40) {eff = 1; errup = 0; errdown = 5.3645e-06;}
    else if (pt > 40 && pt < 50) {eff = 1; errup = 0; errdown = 3.63544e-06;}
    else if (pt > 50 && pt < 80) {eff = 1; errup = 0; errdown = 1.6325e-05;}
    else if (pt > 80 && pt < 120) {eff = 1; errup = 0; errdown = 0.000196966;}
    else if (pt > 120 && pt < 200) {eff = 1; errup = 0; errdown = 0.00056788;}
    else if (pt > 200 && pt < 9999) {eff = 1; errup = 0; errdown = 0.00274402;}
  } else if (abseta > 0.8 && abseta < 1.4442) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 1; errup = 0; errdown = 0.841345;}
    else if (pt > 10 && pt < 15) {eff = 1; errup = 0; errdown = 0.184992;}
    else if (pt > 15 && pt < 20) {eff = 1; errup = 0; errdown = 0.0252456;}
    else if (pt > 20 && pt < 25) {eff = 1; errup = 0; errdown = 0.00694931;}
    else if (pt > 25 && pt < 27) {eff = 1; errup = 0; errdown = 0.00573667;}
    else if (pt > 27 && pt < 30) {eff = 1; errup = 0; errdown = 0.000151825;}
    else if (pt > 30 && pt < 40) {eff = 1; errup = 0; errdown = 7.85275e-06;}
    else if (pt > 40 && pt < 50) {eff = 1; errup = 0; errdown = 5.29723e-06;}
    else if (pt > 50 && pt < 80) {eff = 1; errup = 0; errdown = 2.37129e-05;}
    else if (pt > 80 && pt < 120) {eff = 1; errup = 0; errdown = 0.000284199;}
    else if (pt > 120 && pt < 200) {eff = 1; errup = 0; errdown = 0.000788811;}
    else if (pt > 200 && pt < 9999) {eff = 1; errup = 0; errdown = 0.00413786;}
  } else if (abseta > 1.4442 && abseta < 1.566) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 10 && pt < 15) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 15 && pt < 20) {eff = 1; errup = 0; errdown = 0.15411;}
    else if (pt > 20 && pt < 25) {eff = 1; errup = 0; errdown = 0.00949362;}
    else if (pt > 25 && pt < 27) {eff = 1; errup = 0; errdown = 0.00601796;}
    else if (pt > 27 && pt < 30) {eff = 1; errup = 0; errdown = 0.00125076;}
    else if (pt > 30 && pt < 40) {eff = 1; errup = 0; errdown = 8.0461e-05;}
    else if (pt > 40 && pt < 50) {eff = 1; errup = 0; errdown = 5.94263e-05;}
    else if (pt > 50 && pt < 80) {eff = 1; errup = 0; errdown = 0.000246063;}
    else if (pt > 80 && pt < 120) {eff = 1; errup = 0; errdown = 0.00293661;}
    else if (pt > 120 && pt < 200) {eff = 1; errup = 0; errdown = 0.00708299;}
    else if (pt > 200 && pt < 9999) {eff = 1; errup = 0; errdown = 0.0419109;}
  } else if (abseta > 1.566 && abseta < 2.5) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 1; errup = 0; errdown = 0.458642;}
    else if (pt > 10 && pt < 15) {eff = 1; errup = 0; errdown = 0.184992;}
    else if (pt > 15 && pt < 20) {eff = 1; errup = 0; errdown = 0.0252456;}
    else if (pt > 20 && pt < 25) {eff = 1; errup = 0; errdown = 0.00868727;}
    else if (pt > 25 && pt < 27) {eff = 1; errup = 0; errdown = 0.00657692;}
    else if (pt > 27 && pt < 30) {eff = 1; errup = 0; errdown = 0.000225424;}
    else if (pt > 30 && pt < 40) {eff = 1; errup = 0; errdown = 9.4774e-06;}
    else if (pt > 40 && pt < 50) {eff = 1; errup = 0; errdown = 6.04032e-06;}
    else if (pt > 50 && pt < 80) {eff = 1; errup = 0; errdown = 2.66823e-05;}
    else if (pt > 80 && pt < 120) {eff = 1; errup = 0; errdown = 0.00036234;}
    else if (pt > 120 && pt < 200) {eff = 1; errup = 0; errdown = 0.00113509;}
    else if (pt > 200 && pt < 9999) {eff = 1; errup = 0; errdown = 0.00667227;}
  }
  return std::vector<float>({eff, errup, errdown});
}

std::vector<float> eff_elhltpt23(float pt, float abseta) {
  float errup=0.0, errdown=0.0, eff=1.0;
  if (abseta > 0 && abseta < 0.8) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 10 && pt < 15) {eff = 0.125; errup = 0.141688; errdown = 0.0803106;}
    else if (pt > 15 && pt < 20) {eff = 0.0163934; errup = 0.0156891; errdown = 0.00890888;}
    else if (pt > 20 && pt < 25) {eff = 0.333333; errup = 0.0247703; errdown = 0.0239245;}
    else if (pt > 25 && pt < 27) {eff = 0.986486; errup = 0.00497045; errdown = 0.00720284;}
    else if (pt > 27 && pt < 30) {eff = 0.999959; errup = 3.35502e-05; errdown = 9.32544e-05;}
    else if (pt > 30 && pt < 40) {eff = 1; errup = 0; errdown = 5.3645e-06;}
    else if (pt > 40 && pt < 50) {eff = 1; errup = 0; errdown = 3.63544e-06;}
    else if (pt > 50 && pt < 80) {eff = 1; errup = 0; errdown = 1.6325e-05;}
    else if (pt > 80 && pt < 120) {eff = 1; errup = 0; errdown = 0.000196966;}
    else if (pt > 120 && pt < 200) {eff = 1; errup = 0; errdown = 0.00056788;}
    else if (pt > 200 && pt < 9999) {eff = 1; errup = 0; errdown = 0.00274402;}
  } else if (abseta > 0.8 && abseta < 1.4442) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 0; errup = 0.841345; errdown = 0;}
    else if (pt > 10 && pt < 15) {eff = 0; errup = 0.184992; errdown = 0;}
    else if (pt > 15 && pt < 20) {eff = 0.0277778; errup = 0.0354652; errdown = 0.0179217;}
    else if (pt > 20 && pt < 25) {eff = 0.378788; errup = 0.0321745; errdown = 0.0312003;}
    else if (pt > 25 && pt < 27) {eff = 0.98125; errup = 0.00741272; errdown = 0.0110312;}
    else if (pt > 27 && pt < 30) {eff = 0.999918; errup = 6.82266e-05; errdown = 0.000189626;}
    else if (pt > 30 && pt < 40) {eff = 1; errup = 0; errdown = 7.85275e-06;}
    else if (pt > 40 && pt < 50) {eff = 1; errup = 0; errdown = 5.29723e-06;}
    else if (pt > 50 && pt < 80) {eff = 1; errup = 0; errdown = 2.37129e-05;}
    else if (pt > 80 && pt < 120) {eff = 1; errup = 0; errdown = 0.000284199;}
    else if (pt > 120 && pt < 200) {eff = 1; errup = 0; errdown = 0.000788811;}
    else if (pt > 200 && pt < 9999) {eff = 1; errup = 0; errdown = 0.00413786;}
  } else if (abseta > 1.4442 && abseta < 1.566) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 10 && pt < 15) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 15 && pt < 20) {eff = 0.818182; errup = 0.116508; errdown = 0.191402;}
    else if (pt > 20 && pt < 25) {eff = 0.963731; errup = 0.0132804; errdown = 0.0189895;}
    else if (pt > 25 && pt < 27) {eff = 0.990164; errup = 0.00534846; errdown = 0.00947512;}
    else if (pt > 27 && pt < 30) {eff = 0.99864; errup = 0.00087814; errdown = 0.00179042;}
    else if (pt > 30 && pt < 40) {eff = 1; errup = 0; errdown = 8.0461e-05;}
    else if (pt > 40 && pt < 50) {eff = 1; errup = 0; errdown = 5.94263e-05;}
    else if (pt > 50 && pt < 80) {eff = 1; errup = 0; errdown = 0.000246063;}
    else if (pt > 80 && pt < 120) {eff = 1; errup = 0; errdown = 0.00293661;}
    else if (pt > 120 && pt < 200) {eff = 1; errup = 0; errdown = 0.00708299;}
    else if (pt > 200 && pt < 9999) {eff = 1; errup = 0; errdown = 0.0419109;}
  } else if (abseta > 1.566 && abseta < 2.5) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 0.333333; errup = 0.414535; errdown = 0.277375;}
    else if (pt > 10 && pt < 15) {eff = 0.111111; errup = 0.211558; errdown = 0.0920993;}
    else if (pt > 15 && pt < 20) {eff = 0.0416667; errup = 0.0388735; errdown = 0.0225916;}
    else if (pt > 20 && pt < 25) {eff = 0.341232; errup = 0.0357313; errdown = 0.034122;}
    else if (pt > 25 && pt < 27) {eff = 0.960573; errup = 0.0115904; errdown = 0.0153911;}
    else if (pt > 27 && pt < 30) {eff = 0.999633; errup = 0.000199933; errdown = 0.000357229;}
    else if (pt > 30 && pt < 40) {eff = 1; errup = 0; errdown = 9.4774e-06;}
    else if (pt > 40 && pt < 50) {eff = 1; errup = 0; errdown = 6.04032e-06;}
    else if (pt > 50 && pt < 80) {eff = 1; errup = 0; errdown = 2.66823e-05;}
    else if (pt > 80 && pt < 120) {eff = 1; errup = 0; errdown = 0.00036234;}
    else if (pt > 120 && pt < 200) {eff = 1; errup = 0; errdown = 0.00113509;}
    else if (pt > 200 && pt < 9999) {eff = 1; errup = 0; errdown = 0.00667227;}
  }
  return std::vector<float>({eff, errup, errdown});
}

std::vector<float> eff_dielleg12(float pt, float abseta) {
  float errup=0.0, errdown=0.0, eff=1.0;
  if (abseta > 0 && abseta < 0.8) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 0.0129345; errup = 0.00269325; errdown = 0.0022645;}
    else if (pt > 10 && pt < 15) {eff = 0.457011; errup = 0.00484239; errdown = 0.00483436;}
    else if (pt > 15 && pt < 20) {eff = 0.783823; errup = 0.00261984; errdown = 0.00264286;}
    else if (pt > 20 && pt < 25) {eff = 0.816224; errup = 0.00177837; errdown = 0.00179172;}
    else if (pt > 25 && pt < 27) {eff = 0.832487; errup = 0.00223836; errdown = 0.00226231;}
    else if (pt > 27 && pt < 30) {eff = 0.838551; errup = 0.00156632; errdown = 0.00157861;}
    else if (pt > 30 && pt < 40) {eff = 0.873193; errup = 0.00054017; errdown = 0.00054214;}
    else if (pt > 40 && pt < 50) {eff = 0.889888; errup = 0.000611085; errdown = 0.000614062;}
    else if (pt > 50 && pt < 80) {eff = 0.904008; errup = 0.00259569; errdown = 0.00265902;}
    else if (pt > 80 && pt < 120) {eff = 0.903786; errup = 0.00739647; errdown = 0.00791876;}
    else if (pt > 120 && pt < 200) {eff = 0.927954; errup = 0.013974; errdown = 0.016692;}
    else if (pt > 200 && pt < 9999) {eff = 1; errup = 0; errdown = 0.0802771;}
  } else if (abseta > 0.8 && abseta < 1.4442) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 0.014419; errup = 0.00289578; errdown = 0.00244838;}
    else if (pt > 10 && pt < 15) {eff = 0.354597; errup = 0.00566697; errdown = 0.00562694;}
    else if (pt > 15 && pt < 20) {eff = 0.751839; errup = 0.00339701; errdown = 0.00342816;}
    else if (pt > 20 && pt < 25) {eff = 0.812305; errup = 0.00221919; errdown = 0.0022394;}
    else if (pt > 25 && pt < 27) {eff = 0.827735; errup = 0.00278399; errdown = 0.00281973;}
    else if (pt > 27 && pt < 30) {eff = 0.8437; errup = 0.00188438; errdown = 0.00190294;}
    else if (pt > 30 && pt < 40) {eff = 0.875227; errup = 0.000642558; errdown = 0.0006454;}
    else if (pt > 40 && pt < 50) {eff = 0.897445; errup = 0.000725873; errdown = 0.000730434;}
    else if (pt > 50 && pt < 80) {eff = 0.903035; errup = 0.00313583; errdown = 0.00322737;}
    else if (pt > 80 && pt < 120) {eff = 0.929825; errup = 0.0078062; errdown = 0.00864677;}
    else if (pt > 120 && pt < 200) {eff = 0.958525; errup = 0.0134374; errdown = 0.0183674;}
    else if (pt > 200 && pt < 9999) {eff = 1; errup = 0; errdown = 0.102638;}
  } else if (abseta > 1.4442 && abseta < 1.566) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 0.0106007; errup = 0.0102039; errdown = 0.00576385;}
    else if (pt > 10 && pt < 15) {eff = 0.320755; errup = 0.0163027; errdown = 0.0158909;}
    else if (pt > 15 && pt < 20) {eff = 0.677504; errup = 0.0108681; errdown = 0.0110588;}
    else if (pt > 20 && pt < 25) {eff = 0.783815; errup = 0.00709654; errdown = 0.00726578;}
    else if (pt > 25 && pt < 27) {eff = 0.812594; errup = 0.0088474; errdown = 0.00917132;}
    else if (pt > 27 && pt < 30) {eff = 0.81372; errup = 0.00628643; errdown = 0.00645094;}
    else if (pt > 30 && pt < 40) {eff = 0.852343; errup = 0.00213861; errdown = 0.00216431;}
    else if (pt > 40 && pt < 50) {eff = 0.880169; errup = 0.00247644; errdown = 0.00252093;}
    else if (pt > 50 && pt < 80) {eff = 0.884615; errup = 0.0108525; errdown = 0.0117679;}
    else if (pt > 80 && pt < 120) {eff = 0.902913; errup = 0.029471; errdown = 0.0385366;}
    else if (pt > 120 && pt < 200) {eff = 0.904762; errup = 0.0612701; errdown = 0.112063;}
    else if (pt > 200 && pt < 9999) {eff = 1; errup = 0; errdown = 1;}
  } else if (abseta > 1.566 && abseta < 2.5) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 0.0109307; errup = 0.00216322; errdown = 0.00183179;}
    else if (pt > 10 && pt < 15) {eff = 0.268893; errup = 0.00447329; errdown = 0.00442724;}
    else if (pt > 15 && pt < 20) {eff = 0.710259; errup = 0.00323862; errdown = 0.00326003;}
    else if (pt > 20 && pt < 25) {eff = 0.817144; errup = 0.00211865; errdown = 0.00213774;}
    else if (pt > 25 && pt < 27) {eff = 0.850592; errup = 0.00259424; errdown = 0.00263153;}
    else if (pt > 27 && pt < 30) {eff = 0.874341; errup = 0.00173578; errdown = 0.00175639;}
    else if (pt > 30 && pt < 40) {eff = 0.908584; errup = 0.000597416; errdown = 0.000600935;}
    else if (pt > 40 && pt < 50) {eff = 0.931163; errup = 0.000624897; errdown = 0.000630169;}
    else if (pt > 50 && pt < 80) {eff = 0.939379; errup = 0.00268275; errdown = 0.00279587;}
    else if (pt > 80 && pt < 120) {eff = 0.941463; errup = 0.00823565; errdown = 0.00939174;}
    else if (pt > 120 && pt < 200) {eff = 0.95935; errup = 0.0174478; errdown = 0.0265636;}
    else if (pt > 200 && pt < 9999) {eff = 1; errup = 0; errdown = 0.308024;}
  }
  return std::vector<float>({eff, errup, errdown});
}

std::vector<float> eff_muhltpt8(float pt, float abseta) {
  float errup=0.0, errdown=0.0, eff=1.0;
  if (abseta > 0 && abseta < 0.9) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 10 && pt < 15) {eff = 1; errup = 0; errdown = 1;}
    else if (pt > 15 && pt < 20) {eff = 1; errup = 0; errdown = 0.010465;}
    else if (pt > 20 && pt < 25) {eff = 1; errup = 0; errdown = 4.42842e-05;}
    else if (pt > 25 && pt < 30) {eff = 1; errup = 0; errdown = 1.08576e-05;}
    else if (pt > 30 && pt < 40) {eff = 1; errup = 0; errdown = 2.40087e-06;}
    else if (pt > 40 && pt < 50) {eff = 1; errup = 0; errdown = 1.99445e-06;}
    else if (pt > 50 && pt < 60) {eff = 1; errup = 0; errdown = 8.73674e-06;}
    else if (pt > 60 && pt < 120) {eff = 1; errup = 0; errdown = 1.70205e-05;}
    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.000264366;}
  } else if (abseta > 0.9 && abseta < 1.2) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 10 && pt < 15) {eff = 1; errup = 0; errdown = 1;}
    else if (pt > 15 && pt < 20) {eff = 1; errup = 0; errdown = 0.0141701;}
    else if (pt > 20 && pt < 25) {eff = 1; errup = 0; errdown = 0.000157394;}
    else if (pt > 25 && pt < 30) {eff = 1; errup = 0; errdown = 4.22748e-05;}
    else if (pt > 30 && pt < 40) {eff = 1; errup = 0; errdown = 8.90816e-06;}
    else if (pt > 40 && pt < 50) {eff = 1; errup = 0; errdown = 6.5278e-06;}
    else if (pt > 50 && pt < 60) {eff = 1; errup = 0; errdown = 2.84834e-05;}
    else if (pt > 60 && pt < 120) {eff = 1; errup = 0; errdown = 5.62833e-05;}
    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.000844149;}
  } else if (abseta > 1.2 && abseta < 2.1) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 1; errup = 0; errdown = 0.841345;}
    else if (pt > 10 && pt < 15) {eff = 1; errup = 0; errdown = 0.368878;}
    else if (pt > 15 && pt < 20) {eff = 1; errup = 0; errdown = 0.00352739;}
    else if (pt > 20 && pt < 25) {eff = 1; errup = 0; errdown = 6.39334e-05;}
    else if (pt > 25 && pt < 30) {eff = 1; errup = 0; errdown = 1.72152e-05;}
    else if (pt > 30 && pt < 40) {eff = 1; errup = 0; errdown = 4.12349e-06;}
    else if (pt > 40 && pt < 50) {eff = 1; errup = 0; errdown = 2.97083e-06;}
    else if (pt > 50 && pt < 60) {eff = 1; errup = 0; errdown = 1.29624e-05;}
    else if (pt > 60 && pt < 120) {eff = 1; errup = 0; errdown = 2.62452e-05;}
    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.000452793;}
  } else if (abseta > 2.1 && abseta < 2.4) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 1; errup = 0; errdown = 1;}
    else if (pt > 10 && pt < 15) {eff = 1; errup = 0; errdown = 0.264229;}
    else if (pt > 15 && pt < 20) {eff = 1; errup = 0; errdown = 0.0127035;}
    else if (pt > 20 && pt < 25) {eff = 1; errup = 0; errdown = 0.000406683;}
    else if (pt > 25 && pt < 30) {eff = 1; errup = 0; errdown = 8.13139e-05;}
    else if (pt > 30 && pt < 40) {eff = 1; errup = 0; errdown = 1.92119e-05;}
    else if (pt > 40 && pt < 50) {eff = 1; errup = 0; errdown = 1.7265e-05;}
    else if (pt > 50 && pt < 60) {eff = 1; errup = 0; errdown = 7.4936e-05;}
    else if (pt > 60 && pt < 120) {eff = 1; errup = 0; errdown = 0.000163662;}
    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.0035963;}
  }
  return std::vector<float>({eff, errup, errdown});
}

std::vector<float> eff_muhltpt17(float pt, float abseta) {
  float errup=0.0, errdown=0.0, eff=1.0;
  if (abseta > 0 && abseta < 0.9) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 10 && pt < 15) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 15 && pt < 20) {eff = 1; errup = 0; errdown = 0.010465;}
    else if (pt > 20 && pt < 25) {eff = 1; errup = 0; errdown = 4.42842e-05;}
    else if (pt > 25 && pt < 30) {eff = 1; errup = 0; errdown = 1.08576e-05;}
    else if (pt > 30 && pt < 40) {eff = 1; errup = 0; errdown = 2.40087e-06;}
    else if (pt > 40 && pt < 50) {eff = 1; errup = 0; errdown = 1.99445e-06;}
    else if (pt > 50 && pt < 60) {eff = 1; errup = 0; errdown = 8.73674e-06;}
    else if (pt > 60 && pt < 120) {eff = 1; errup = 0; errdown = 1.70205e-05;}
    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.000264366;}
  } else if (abseta > 0.9 && abseta < 1.2) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 10 && pt < 15) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 15 && pt < 20) {eff = 1; errup = 0; errdown = 0.0141701;}
    else if (pt > 20 && pt < 25) {eff = 1; errup = 0; errdown = 0.000157394;}
    else if (pt > 25 && pt < 30) {eff = 1; errup = 0; errdown = 4.22748e-05;}
    else if (pt > 30 && pt < 40) {eff = 1; errup = 0; errdown = 8.90816e-06;}
    else if (pt > 40 && pt < 50) {eff = 1; errup = 0; errdown = 6.5278e-06;}
    else if (pt > 50 && pt < 60) {eff = 1; errup = 0; errdown = 2.84834e-05;}
    else if (pt > 60 && pt < 120) {eff = 1; errup = 0; errdown = 5.62833e-05;}
    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.000844149;}
  } else if (abseta > 1.2 && abseta < 2.1) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 1; errup = 0; errdown = 0.841345;}
    else if (pt > 10 && pt < 15) {eff = 1; errup = 0; errdown = 0.368878;}
    else if (pt > 15 && pt < 20) {eff = 1; errup = 0; errdown = 0.00352739;}
    else if (pt > 20 && pt < 25) {eff = 1; errup = 0; errdown = 6.39334e-05;}
    else if (pt > 25 && pt < 30) {eff = 1; errup = 0; errdown = 1.72152e-05;}
    else if (pt > 30 && pt < 40) {eff = 1; errup = 0; errdown = 4.12349e-06;}
    else if (pt > 40 && pt < 50) {eff = 1; errup = 0; errdown = 2.97083e-06;}
    else if (pt > 50 && pt < 60) {eff = 1; errup = 0; errdown = 1.29624e-05;}
    else if (pt > 60 && pt < 120) {eff = 1; errup = 0; errdown = 2.62452e-05;}
    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.000452793;}
  } else if (abseta > 2.1 && abseta < 2.4) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 0; errup = 1; errdown = 0;}
    else if (pt > 10 && pt < 15) {eff = 1; errup = 0; errdown = 0.264229;}
    else if (pt > 15 && pt < 20) {eff = 1; errup = 0; errdown = 0.0127035;}
    else if (pt > 20 && pt < 25) {eff = 1; errup = 0; errdown = 0.000406683;}
    else if (pt > 25 && pt < 30) {eff = 1; errup = 0; errdown = 8.13139e-05;}
    else if (pt > 30 && pt < 40) {eff = 1; errup = 0; errdown = 1.92119e-05;}
    else if (pt > 40 && pt < 50) {eff = 1; errup = 0; errdown = 1.7265e-05;}
    else if (pt > 50 && pt < 60) {eff = 1; errup = 0; errdown = 7.4936e-05;}
    else if (pt > 60 && pt < 120) {eff = 1; errup = 0; errdown = 0.000163662;}
    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.0035963;}
  }
  return std::vector<float>({eff, errup, errdown});
}

std::vector<float> eff_dimuleg8(float pt, float abseta) {
  float errup=0.0, errdown=0.0, eff=1.0;
  if (abseta > 0 && abseta < 0.9) {
    if (pt > 5 && pt < 7) {eff = 0.00385356; errup = 0.00505989; errdown = 0.00248866;}
    else if (pt > 7 && pt < 10) {eff = 0.697792; errup = 0.00915289; errdown = 0.00930941;}
    else if (pt > 10 && pt < 15) {eff = 0.886263; errup = 0.00249071; errdown = 0.00253858;}
    else if (pt > 15 && pt < 20) {eff = 0.906334; errup = 0.00138409; errdown = 0.00140252;}
    else if (pt > 20 && pt < 25) {eff = 0.912274; errup = 0.000898098; errdown = 0.000906438;}
    else if (pt > 25 && pt < 30) {eff = 0.912271; errup = 0.000641987; errdown = 0.000646245;}
    else if (pt > 30 && pt < 40) {eff = 0.911749; errup = 0.000348264; errdown = 0.000349507;}
    else if (pt > 40 && pt < 50) {eff = 0.910437; errup = 0.000441892; errdown = 0.000443861;}
    else if (pt > 50 && pt < 60) {eff = 0.91305; errup = 0.00251375; errdown = 0.0025802;}
    else if (pt > 60 && pt < 120) {eff = 0.900718; errup = 0.00333859; errdown = 0.00343965;}
    else if (pt > 120 && pt < 9999) {eff = 0.867391; errup = 0.0160412; errdown = 0.0177462;}
  } else if (abseta > 0.9 && abseta < 1.2) {
    if (pt > 5 && pt < 7) {eff = 0.00197628; errup = 0.00452972; errdown = 0.00163493;}
    else if (pt > 7 && pt < 10) {eff = 0.697398; errup = 0.0110908; errdown = 0.0113198;}
    else if (pt > 10 && pt < 15) {eff = 0.901621; errup = 0.00356399; errdown = 0.00368046;}
    else if (pt > 15 && pt < 20) {eff = 0.911847; errup = 0.00234571; errdown = 0.00240264;}
    else if (pt > 20 && pt < 25) {eff = 0.922572; errup = 0.00159069; errdown = 0.00162087;}
    else if (pt > 25 && pt < 30) {eff = 0.923591; errup = 0.00117866; errdown = 0.00119543;}
    else if (pt > 30 && pt < 40) {eff = 0.926994; errup = 0.000601835; errdown = 0.00060642;}
    else if (pt > 40 && pt < 50) {eff = 0.928959; errup = 0.000711608; errdown = 0.000718217;}
    else if (pt > 50 && pt < 60) {eff = 0.929361; errup = 0.00422294; errdown = 0.00446185;}
    else if (pt > 60 && pt < 120) {eff = 0.921064; errup = 0.00552283; errdown = 0.00588584;}
    else if (pt > 120 && pt < 9999) {eff = 0.917293; errup = 0.0240424; errdown = 0.0312689;}
  } else if (abseta > 1.2 && abseta < 2.1) {
    if (pt > 5 && pt < 7) {eff = 0.00661813; errup = 0.00182762; errdown = 0.00146445;}
    else if (pt > 7 && pt < 10) {eff = 0.623311; errup = 0.00586251; errdown = 0.00589839;}
    else if (pt > 10 && pt < 15) {eff = 0.891151; errup = 0.00210883; errdown = 0.00214492;}
    else if (pt > 15 && pt < 20) {eff = 0.913368; errup = 0.0013848; errdown = 0.00140496;}
    else if (pt > 20 && pt < 25) {eff = 0.91741; errup = 0.000997062; errdown = 0.00100806;}
    else if (pt > 25 && pt < 30) {eff = 0.91911; errup = 0.000758351; errdown = 0.000764857;}
    else if (pt > 30 && pt < 40) {eff = 0.923037; errup = 0.000416828; errdown = 0.000418901;}
    else if (pt > 40 && pt < 50) {eff = 0.925504; errup = 0.000474513; errdown = 0.000477299;}
    else if (pt > 50 && pt < 60) {eff = 0.925551; errup = 0.00288916; errdown = 0.00299383;}
    else if (pt > 60 && pt < 120) {eff = 0.926176; errup = 0.00375231; errdown = 0.00393136;}
    else if (pt > 120 && pt < 9999) {eff = 0.918919; errup = 0.017104; errdown = 0.0206986;}
  } else if (abseta > 2.1 && abseta < 2.4) {
    if (pt > 5 && pt < 7) {eff = 0.00426621; errup = 0.00287575; errdown = 0.00184154;}
    else if (pt > 7 && pt < 10) {eff = 0.565879; errup = 0.0100901; errdown = 0.0101439;}
    else if (pt > 10 && pt < 15) {eff = 0.868991; errup = 0.00405957; errdown = 0.00416734;}
    else if (pt > 15 && pt < 20) {eff = 0.884012; errup = 0.00290405; errdown = 0.00296771;}
    else if (pt > 20 && pt < 25) {eff = 0.894054; errup = 0.00210315; errdown = 0.00214019;}
    else if (pt > 25 && pt < 30) {eff = 0.89672; errup = 0.00160715; errdown = 0.00162939;}
    else if (pt > 30 && pt < 40) {eff = 0.899588; errup = 0.000962316; errdown = 0.000970536;}
    else if (pt > 40 && pt < 50) {eff = 0.901388; errup = 0.00126356; errdown = 0.00127804;}
    else if (pt > 50 && pt < 60) {eff = 0.907641; errup = 0.00751323; errdown = 0.00807885;}
    else if (pt > 60 && pt < 120) {eff = 0.901299; errup = 0.0108462; errdown = 0.0119507;}
    else if (pt > 120 && pt < 9999) {eff = 0.941176; errup = 0.048713; errdown = 0.12258;}
  }
  return std::vector<float>({eff, errup, errdown});
}

std::vector<float> eff_isoel3235(float pt, float abseta) {
  float errup=0.0, errdown=0.0, eff=1.0;
  if (abseta > 0 && abseta < 0.8) {
    if (pt > 5 && pt < 7) {eff = 0.00140449; errup = 0.00322219; errdown = 0.00116189;}
    else if (pt > 7 && pt < 10) {eff = 0.00103734; errup = 0.00238132; errdown = 0.000858155;}
    else if (pt > 10 && pt < 15) {eff = 0.00156986; errup = 0.00360059; errdown = 0.0012987;}
    else if (pt > 15 && pt < 20) {eff = 0; errup = 0.00357535; errdown = 0;}
    else if (pt > 20 && pt < 25) {eff = 0; errup = 0.00984914; errdown = 0;}
    else if (pt > 25 && pt < 27) {eff = 0; errup = 0.00713791; errdown = 0;}
    else if (pt > 27 && pt < 30) {eff = 0.519118; errup = 0.0198511; errdown = 0.0199094;}
    else if (pt > 30 && pt < 40) {eff = 0.753176; errup = 0.0187794; errdown = 0.0197427;}
    else if (pt > 40 && pt < 50) {eff = 0.79952; errup = 0.0140911; errdown = 0.0148408;}
    else if (pt > 50 && pt < 80) {eff = 0.821168; errup = 0.0166547; errdown = 0.0178908;}
    else if (pt > 80 && pt < 120) {eff = 0.881057; errup = 0.0153914; errdown = 0.0171894;}
    else if (pt > 120 && pt < 200) {eff = 0.845638; errup = 0.0303193; errdown = 0.0354625;}
    else if (pt > 200 && pt < 9999) {eff = 0; errup = 0; errdown = 0;}
  } else if (abseta > 0.8 && abseta < 1.4442) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 0.00385215; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 0.00167785; errup = 0.0038476; errdown = 0.00138804;}
    else if (pt > 10 && pt < 15) {eff = 0; errup = 0.00418489; errdown = 0;}
    else if (pt > 15 && pt < 20) {eff = 0; errup = 0.00459198; errdown = 0;}
    else if (pt > 20 && pt < 25) {eff = 0; errup = 0.0170586; errdown = 0;}
    else if (pt > 25 && pt < 27) {eff = 0; errup = 0.0107085; errdown = 0;}
    else if (pt > 27 && pt < 30) {eff = 0.397614; errup = 0.023005; errdown = 0.0225799;}
    else if (pt > 30 && pt < 40) {eff = 0.750716; errup = 0.0237977; errdown = 0.0253206;}
    else if (pt > 40 && pt < 50) {eff = 0.781991; errup = 0.0167243; errdown = 0.0176572;}
    else if (pt > 50 && pt < 80) {eff = 0.80826; errup = 0.0218627; errdown = 0.0238032;}
    else if (pt > 80 && pt < 120) {eff = 0.863333; errup = 0.0201738; errdown = 0.0227942;}
    else if (pt > 120 && pt < 200) {eff = 0.780488; errup = 0.0476882; errdown = 0.0553797;}
    else if (pt > 200 && pt < 9999) {eff = 0.519118; errup = 0.0153914; errdown = 0.0171894;}
  } else if (abseta > 1.4442 && abseta < 1.566) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 0.0233265; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 0; errup = 0.0184243; errdown = 0;}
    else if (pt > 10 && pt < 15) {eff = 0; errup = 0.0230346; errdown = 0;}
    else if (pt > 15 && pt < 20) {eff = 0; errup = 0.0368748; errdown = 0;}
    else if (pt > 20 && pt < 25) {eff = 0; errup = 0.0923495; errdown = 0;}
    else if (pt > 25 && pt < 27) {eff = 0.0416667; errup = 0.0893855; errdown = 0.0344944;}
    else if (pt > 27 && pt < 30) {eff = 0.381818; errup = 0.0762757; errdown = 0.0714409;}
    else if (pt > 30 && pt < 40) {eff = 0.552632; errup = 0.0908274; errdown = 0.0939923;}
    else if (pt > 40 && pt < 50) {eff = 0.707692; errup = 0.0600147; errdown = 0.0672026;}
    else if (pt > 50 && pt < 80) {eff = 0.777778; errup = 0.0729951; errdown = 0.0911652;}
    else if (pt > 80 && pt < 120) {eff = 0.681818; errup = 0.109079; errdown = 0.12875;}
    else if (pt > 120 && pt < 200) {eff = 0.833333; errup = 0.138285; errdown = 0.28735;}
    else if (pt > 200 && pt < 9999) {eff = 0.519118; errup = 0.0153914; errdown = 0.0171894;}
  } else if (abseta > 1.566 && abseta < 2.5) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 0.00315827; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 0.00145138; errup = 0.0033295; errdown = 0.00120068;}
    else if (pt > 10 && pt < 15) {eff = 0; errup = 0.00378094; errdown = 0;}
    else if (pt > 15 && pt < 20) {eff = 0; errup = 0.00507278; errdown = 0;}
    else if (pt > 20 && pt < 25) {eff = 0; errup = 0.0133482; errdown = 0;}
    else if (pt > 25 && pt < 27) {eff = 0.00537634; errup = 0.0122537; errdown = 0.00444799;}
    else if (pt > 27 && pt < 30) {eff = 0.282222; errup = 0.0228143; errdown = 0.0217971;}
    else if (pt > 30 && pt < 40) {eff = 0.656627; errup = 0.0270269; errdown = 0.028023;}
    else if (pt > 40 && pt < 50) {eff = 0.737705; errup = 0.0192155; errdown = 0.0201225;}
    else if (pt > 50 && pt < 80) {eff = 0.771127; errup = 0.0256314; errdown = 0.0276711;}
    else if (pt > 80 && pt < 120) {eff = 0.795276; errup = 0.0259964; errdown = 0.028496;}
    else if (pt > 120 && pt < 200) {eff = 0.853659; errup = 0.0564709; errdown = 0.0771695;}
    else if (pt > 200 && pt < 9999) {eff = 0.519118; errup = 0.0153914; errdown = 0.0171894;}
  }
  return std::vector<float>({eff, errup, errdown});
}

std::vector<float> eff_isomu2427(float pt, float abseta) {
  float errup=0.0, errdown=0.0, eff=1.0;
  if (abseta > 0 && abseta < 0.9) {
    if (pt > 5 && pt < 7) {eff = 0.0016; errup = 0.00366954; errdown = 0.00132363;}
    else if (pt > 7 && pt < 10) {eff = 0.00414938; errup = 0.00401971; errdown = 0.00225741;}
    else if (pt > 10 && pt < 15) {eff = 0.00178253; errup = 0.00234618; errdown = 0.00115127;}
    else if (pt > 15 && pt < 20) {eff = 0.00117371; errup = 0.00269375; errdown = 0.000970967;}
    else if (pt > 20 && pt < 25) {eff = 0.103107; errup = 0.0127297; errdown = 0.0115402;}
    else if (pt > 25 && pt < 30) {eff = 0.781513; errup = 0.0172733; errdown = 0.0182655;}
    else if (pt > 30 && pt < 40) {eff = 0.830366; errup = 0.0116332; errdown = 0.0122773;}
    else if (pt > 40 && pt < 50) {eff = 0.861438; errup = 0.0126449; errdown = 0.0136385;}
    else if (pt > 50 && pt < 60) {eff = 0.866559; errup = 0.0138113; errdown = 0.0150583;}
    else if (pt > 60 && pt < 120) {eff = 0.857143; errup = 0.0091918; errdown = 0.0096932;}
    else if (pt > 120 && pt < 9999) {eff = 0.808929; errup = 0.0169152; errdown = 0.0180766;}
  } else if (abseta > 0.9 && abseta < 1.2) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 0.0115844; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 0; errup = 0.0082586; errdown = 0;}
    else if (pt > 10 && pt < 15) {eff = 0; errup = 0.0055299; errdown = 0;}
    else if (pt > 15 && pt < 20) {eff = 0.00416667; errup = 0.00951566; errdown = 0.00344712;}
    else if (pt > 20 && pt < 25) {eff = 0.0809524; errup = 0.0234566; errdown = 0.0189743;}
    else if (pt > 25 && pt < 30) {eff = 0.698925; errup = 0.0350786; errdown = 0.0373806;}
    else if (pt > 30 && pt < 40) {eff = 0.772727; errup = 0.0254723; errdown = 0.0275096;}
    else if (pt > 40 && pt < 50) {eff = 0.816964; errup = 0.026494; errdown = 0.0295599;}
    else if (pt > 50 && pt < 60) {eff = 0.819672; errup = 0.0291864; errdown = 0.0330027;}
    else if (pt > 60 && pt < 120) {eff = 0.851685; errup = 0.0171206; errdown = 0.0188043;}
    else if (pt > 120 && pt < 9999) {eff = 0.845588; errup = 0.0317581; errdown = 0.0374164;}
  } else if (abseta > 1.2 && abseta < 2.1) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 0.00350057; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 0; errup = 0.00291799; errdown = 0;}
    else if (pt > 10 && pt < 15) {eff = 0.00235018; errup = 0.00309128; errdown = 0.00151785;}
    else if (pt > 15 && pt < 20) {eff = 0.0030722; errup = 0.00403759; errdown = 0.00198411;}
    else if (pt > 20 && pt < 25) {eff = 0.100182; errup = 0.0145043; errdown = 0.0129459;}
    else if (pt > 25 && pt < 30) {eff = 0.651584; errup = 0.0234085; errdown = 0.0241276;}
    else if (pt > 30 && pt < 40) {eff = 0.76738; errup = 0.0157356; errdown = 0.0164808;}
    else if (pt > 40 && pt < 50) {eff = 0.790927; errup = 0.0184278; errdown = 0.0196367;}
    else if (pt > 50 && pt < 60) {eff = 0.826316; errup = 0.0198178; errdown = 0.0216486;}
    else if (pt > 60 && pt < 120) {eff = 0.823654; errup = 0.0125525; errdown = 0.0132647;}
    else if (pt > 120 && pt < 9999) {eff = 0.830325; errup = 0.0230458; errdown = 0.025616;}
  } else if (abseta > 2.1 && abseta < 2.4) {
    if (pt > 5 && pt < 7) {eff = 0.0125; errup = 0.0162484; errdown = 0.00806978;}
    else if (pt > 7 && pt < 10) {eff = 0; errup = 0.0107085; errdown = 0;}
    else if (pt > 10 && pt < 15) {eff = 0; errup = 0.00818516; errdown = 0;}
    else if (pt > 15 && pt < 20) {eff = 0; errup = 0.0105853; errdown = 0;}
    else if (pt > 20 && pt < 25) {eff = 0.0538462; errup = 0.0277989; errdown = 0.0196464;}
    else if (pt > 25 && pt < 30) {eff = 0.447368; errup = 0.0512155; errdown = 0.0502109;}
    else if (pt > 30 && pt < 40) {eff = 0.673077; errup = 0.0394325; errdown = 0.041831;}
    else if (pt > 40 && pt < 50) {eff = 0.734694; errup = 0.0468096; errdown = 0.0521118;}
    else if (pt > 50 && pt < 60) {eff = 0.696203; errup = 0.0549362; errdown = 0.0604646;}
    else if (pt > 60 && pt < 120) {eff = 0.725581; errup = 0.0315616; errdown = 0.0338138;}
    else if (pt > 120 && pt < 9999) {eff = 0.625; errup = 0.0842352; errdown = 0.0913836;}
  }
  return std::vector<float>({eff, errup, errdown});
}

float eff_hig19014_ele12(float pt, float abseta) {
  //2017 values
  float eff=0.0;
  if (abseta > 0 && abseta < 0.8) {
    if (pt > 5 && pt < 7) eff = 0.0;
    else if (pt > 7 && pt < 12) eff = 0.0;
    else if (pt > 12 && pt < 15) eff = 0.0;
    else if (pt > 15 && pt < 17) eff = 0.775; 
    else if (pt > 17 && pt < 20) eff = 0.795;
    else if (pt > 20 && pt < 25) eff = 0.805;
    else if (pt > 25 && pt < 27) eff = 0.820;
    else if (pt > 27 && pt < 30) eff = 0.845;
    else if (pt > 30 && pt < 40) eff = 0.880;
    else if (pt > 40 && pt < 50) eff = 0.910;
    else if (pt > 50 && pt < 80) eff = 0.920;
    else if (pt > 80 && pt < 120) eff = 0.945;
    else if (pt > 120 && pt < 200) eff = 0.950;
    else if (pt > 200 && pt < 9999) eff = 0.955;
  } else if (abseta > 0.8 && abseta < 1.444) {
    if (pt > 5 && pt < 7) eff = 0.0;
    else if (pt > 7 && pt < 12) eff = 0.0;
    else if (pt > 12 && pt < 15) eff = 0.0;
    else if (pt > 15 && pt < 17) eff = 0.695; 
    else if (pt > 17 && pt < 20) eff = 0.755;
    else if (pt > 20 && pt < 25) eff = 0.795;
    else if (pt > 25 && pt < 27) eff = 0.820;
    else if (pt > 27 && pt < 30) eff = 0.840;
    else if (pt > 30 && pt < 40) eff = 0.870;
    else if (pt > 40 && pt < 50) eff = 0.905;
    else if (pt > 50 && pt < 80) eff = 0.915;
    else if (pt > 80 && pt < 120) eff = 0.945;
    else if (pt > 120 && pt < 200) eff = 0.945;
    else if (pt > 200 && pt < 9999) eff = 0.955;
  } else if (abseta > 1.444 && abseta < 1.566) {
    if (pt > 5 && pt < 7) eff = 0.0;
    else if (pt > 7 && pt < 12) eff = 0.0;
    else if (pt > 12 && pt < 15) eff = 0.0;
    else if (pt > 15 && pt < 17) eff = 0.500; 
    else if (pt > 17 && pt < 20) eff = 0.550;
    else if (pt > 20 && pt < 25) eff = 0.655;
    else if (pt > 25 && pt < 27) eff = 0.710;
    else if (pt > 27 && pt < 30) eff = 0.750;
    else if (pt > 30 && pt < 40) eff = 0.895;
    else if (pt > 40 && pt < 50) eff = 0.850;
    else if (pt > 50 && pt < 80) eff = 0.805;
    else if (pt > 80 && pt < 120) eff = 0.850;
    else if (pt > 120 && pt < 200) eff = 0.895;
    else if (pt > 200 && pt < 9999) eff = 0.910;
  } else if (abseta > 1.566 && abseta < 2.500) {
    if (pt > 5 && pt < 7) eff = 0.0;
    else if (pt > 7 && pt < 12) eff = 0.0;
    else if (pt > 12 && pt < 15) eff = 0.340;
    else if (pt > 15 && pt < 17) eff = 0.575; 
    else if (pt > 17 && pt < 20) eff = 0.715;
    else if (pt > 20 && pt < 25) eff = 0.805;
    else if (pt > 25 && pt < 27) eff = 0.850;
    else if (pt > 27 && pt < 30) eff = 0.855;
    else if (pt > 30 && pt < 40) eff = 0.900;
    else if (pt > 40 && pt < 50) eff = 0.940;
    else if (pt > 50 && pt < 80) eff = 0.950;
    else if (pt > 80 && pt < 120) eff = 0.960;
    else if (pt > 120 && pt < 200) eff = 0.960;
    else if (pt > 200 && pt < 9999) eff = 0.980;
  }
  return eff;
}


float eff_hig19014_ele23(float pt, float abseta) {
  //2017 values
  float eff=1.0;
  if (abseta > 0 && abseta < 0.8) {
    if (pt > 5 && pt < 10) eff = 0.0;
    else if (pt > 10 && pt < 15) eff = 0.0;
    else if (pt > 15 && pt < 17) eff = 0.0; 
    else if (pt > 17 && pt < 20) eff = 0.0;
    else if (pt > 20 && pt < 25) eff = 0.240;
    else if (pt > 25 && pt < 27) eff = 0.805;
    else if (pt > 27 && pt < 30) eff = 0.840;
    else if (pt > 30 && pt < 40) eff = 0.880;
    else if (pt > 40 && pt < 50) eff = 0.910;
    else if (pt > 50 && pt < 80) eff = 0.920;
    else if (pt > 80 && pt < 120) eff = 0.940;
    else if (pt > 120 && pt < 200) eff = 0.945;
    else if (pt > 200 && pt < 9999) eff = 0.955;
  } else if (abseta > 0.8 && abseta < 1.444) {
    if (pt > 5 && pt < 10) eff = 0.0;
    else if (pt > 10 && pt < 15) eff = 0.0;
    else if (pt > 17 && pt < 20) eff = 0.0;
    else if (pt > 20 && pt < 25) eff = 0.180;
    else if (pt > 25 && pt < 27) eff = 0.760;
    else if (pt > 27 && pt < 30) eff = 0.840;
    else if (pt > 30 && pt < 40) eff = 0.860;
    else if (pt > 40 && pt < 50) eff = 0.905;
    else if (pt > 50 && pt < 80) eff = 0.910;
    else if (pt > 80 && pt < 120) eff = 0.940;
    else if (pt > 120 && pt < 200) eff = 0.945;
    else if (pt > 200 && pt < 9999) eff = 0.955;
  } else if (abseta > 1.444 && abseta < 1.566) {
    if (pt > 5 && pt < 10) eff = 0.0;
    else if (pt > 10 && pt < 15) eff = 0.0;
    else if (pt > 17 && pt < 20) eff = 0.0;
    else if (pt > 20 && pt < 25) eff = 0.255;
    else if (pt > 25 && pt < 27) eff = 0.645;
    else if (pt > 27 && pt < 30) eff = 0.700;
    else if (pt > 30 && pt < 40) eff = 0.795;
    else if (pt > 40 && pt < 50) eff = 0.840;
    else if (pt > 50 && pt < 80) eff = 0.805;
    else if (pt > 80 && pt < 120) eff = 0.825;
    else if (pt > 120 && pt < 200) eff = 0.900;
    else if (pt > 200 && pt < 9999) eff = 0.910;
  } else if (abseta > 1.566 && abseta < 2.500) {
    if (pt > 5 && pt < 10) eff = 0.0;
    else if (pt > 10 && pt < 15) eff = 0.0;
    else if (pt > 17 && pt < 20) eff = 0.005;
    else if (pt > 20 && pt < 25) eff = 0.150;
    else if (pt > 25 && pt < 27) eff = 0.795;
    else if (pt > 27 && pt < 30) eff = 0.845;
    else if (pt > 30 && pt < 40) eff = 0.900;
    else if (pt > 40 && pt < 50) eff = 0.945;
    else if (pt > 50 && pt < 80) eff = 0.950;
    else if (pt > 80 && pt < 120) eff = 0.965;
    else if (pt > 120 && pt < 200) eff = 0.965;
    else if (pt > 200 && pt < 9999) eff = 0.980;
  }
  return eff;
}

float eff_hig19014_mu8(float pt, float abseta) {
  //2017 values
  float eff=1.0;
  if (abseta > 0 && abseta < 0.9) {
    if (pt > 5 && pt < 8) eff = 0.025;
    else if (pt > 8 && pt < 10) eff = 0.920;
    else if (pt > 10 && pt < 15) eff = 0.925;
    else if (pt > 15 && pt < 20) eff = 0.925; 
    else if (pt > 20 && pt < 25) eff = 0.925;
    else if (pt > 25 && pt < 30) eff = 0.930;
    else if (pt > 30 && pt < 40) eff = 0.930;
    else if (pt > 40 && pt < 50) eff = 0.930;
    else if (pt > 50 && pt < 60) eff = 0.925;
    else if (pt > 60 && pt < 120) eff = 0.920;
    else if (pt > 120 && pt < 9999) eff = 0.905;
  } else if (abseta > 0.9 && abseta < 1.2) {
    if (pt > 5 && pt < 8) eff = 0.025;
    else if (pt > 8 && pt < 10) eff = 0.920;
    else if (pt > 10 && pt < 15) eff = 0.925;
    else if (pt > 15 && pt < 20) eff = 0.930; 
    else if (pt > 20 && pt < 25) eff = 0.935;
    else if (pt > 25 && pt < 30) eff = 0.935;
    else if (pt > 30 && pt < 40) eff = 0.935;
    else if (pt > 40 && pt < 50) eff = 0.940;
    else if (pt > 50 && pt < 60) eff = 0.940;
    else if (pt > 60 && pt < 120) eff = 0.940;
    else if (pt > 120 && pt < 9999) eff = 0.920;
  } else if (abseta > 1.2 && abseta < 2.1) {
    if (pt > 5 && pt < 8) eff = 0.020;
    else if (pt > 8 && pt < 10) eff = 0.900;
    else if (pt > 10 && pt < 15) eff = 0.925;
    else if (pt > 15 && pt < 20) eff = 0.930; 
    else if (pt > 20 && pt < 25) eff = 0.935;
    else if (pt > 25 && pt < 30) eff = 0.935;
    else if (pt > 30 && pt < 40) eff = 0.935;
    else if (pt > 40 && pt < 50) eff = 0.940;
    else if (pt > 50 && pt < 60) eff = 0.940;
    else if (pt > 60 && pt < 120) eff = 0.945;
    else if (pt > 120 && pt < 9999) eff = 0.945;
  } else if (abseta > 2.1 && abseta < 2.4) {
    if (pt > 5 && pt < 8) eff = 0.040;
    else if (pt > 8 && pt < 10) eff = 0.840;
    else if (pt > 10 && pt < 15) eff = 0.890;
    else if (pt > 15 && pt < 20) eff = 0.905; 
    else if (pt > 20 && pt < 25) eff = 0.920;
    else if (pt > 25 && pt < 30) eff = 0.920;
    else if (pt > 30 && pt < 40) eff = 0.925;
    else if (pt > 40 && pt < 50) eff = 0.925;
    else if (pt > 50 && pt < 60) eff = 0.925;
    else if (pt > 60 && pt < 120) eff = 0.920;
    else if (pt > 120 && pt < 9999) eff = 0.905;
  }
  return eff;
}

float eff_hig19014_mu17(float pt, float abseta) {
  //2017 values
  float eff=1.0;
  if (abseta > 0 && abseta < 0.9) {
    if (pt > 5 && pt < 10) eff = 0.0;
    else if (pt > 10 && pt < 15) eff = 0.0;
    else if (pt > 15 && pt < 17) eff = 0.015;
    else if (pt > 17 && pt < 20) eff = 0.905; 
    else if (pt > 20 && pt < 25) eff = 0.925;
    else if (pt > 25 && pt < 30) eff = 0.925;
    else if (pt > 30 && pt < 40) eff = 0.925;
    else if (pt > 40 && pt < 50) eff = 0.925;
    else if (pt > 50 && pt < 60) eff = 0.920;
    else if (pt > 60 && pt < 120) eff = 0.920;
    else if (pt > 120 && pt < 9999) eff = 0.905;
  } else if (abseta > 0.9 && abseta < 1.2) {
    if (pt > 5 && pt < 10) eff = 0.0;
    else if (pt > 10 && pt < 15) eff = 0.0;
    else if (pt > 15 && pt < 17) eff = 0.040;
    else if (pt > 17 && pt < 20) eff = 0.905; 
    else if (pt > 20 && pt < 25) eff = 0.920;
    else if (pt > 25 && pt < 30) eff = 0.940;
    else if (pt > 30 && pt < 40) eff = 0.940;
    else if (pt > 40 && pt < 50) eff = 0.940;
    else if (pt > 50 && pt < 60) eff = 0.940;
    else if (pt > 60 && pt < 120) eff = 0.940;
    else if (pt > 120 && pt < 9999) eff = 0.920;
  } else if (abseta > 1.2 && abseta < 2.1) {
    if (pt > 5 && pt < 10) eff = 0.0;
    else if (pt > 10 && pt < 15) eff = 0.0;
    else if (pt > 15 && pt < 17) eff = 0.050;
    else if (pt > 17 && pt < 20) eff = 0.865; 
    else if (pt > 20 && pt < 25) eff = 0.915;
    else if (pt > 25 && pt < 30) eff = 0.920;
    else if (pt > 30 && pt < 40) eff = 0.925;
    else if (pt > 40 && pt < 50) eff = 0.940;
    else if (pt > 50 && pt < 60) eff = 0.940;
    else if (pt > 60 && pt < 120) eff = 0.940;
    else if (pt > 120 && pt < 9999) eff = 0.940;
  } else if (abseta > 2.1 && abseta < 2.4) {
    if (pt > 5 && pt < 10) eff = 0.0;
    else if (pt > 10 && pt < 15) eff = 0.0;
    else if (pt > 15 && pt < 17) eff = 0.190;
    else if (pt > 17 && pt < 20) eff = 0.710; 
    else if (pt > 20 && pt < 25) eff = 0.830;
    else if (pt > 25 && pt < 30) eff = 0.850;
    else if (pt > 30 && pt < 40) eff = 0.870;
    else if (pt > 40 && pt < 50) eff = 0.885;
    else if (pt > 50 && pt < 60) eff = 0.895;
    else if (pt > 60 && pt < 120) eff = 0.900;
    else if (pt > 120 && pt < 9999) eff = 0.885;
  }
  return eff;
}
