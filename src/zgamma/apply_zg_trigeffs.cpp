#include <vector>

#include "zgamma/apply_zg_trigeffs.hpp"

//std::vector<float> eff_elhltpt23(float pt, float abseta) {
//  float errup=0.0, errdown=0.0, eff=1.0;
//  if (abseta > 0 && abseta < 0.8) {
//    if (pt > 7 && pt < 10) {eff = 0; errup = 0.000744473; errdown = 0;}
//    else if (pt > 10 && pt < 12) {eff = 0; errup = 0.000605416; errdown = 0;}
//    else if (pt > 12 && pt < 15) {eff = 0.000128783; errup = 0.000296077; errdown = 0.000106535;}
//    else if (pt > 15 && pt < 20) {eff = 0.000120739; errup = 0.000117432; errdown = 6.57096e-05;}
//    else if (pt > 20 && pt < 25) {eff = 0.00321233; errup = 0.000280772; errdown = 0.000259065;}
//    else if (pt > 25 && pt < 27) {eff = 0.0192974; errup = 0.000857901; errdown = 0.000822798;}
//    else if (pt > 27 && pt < 30) {eff = 0.457064; errup = 0.00212608; errdown = 0.00212452;}
//    else if (pt > 30 && pt < 35) {eff = 0.801033; errup = 0.00104061; errdown = 0.0010447;}
//    else if (pt > 35 && pt < 40) {eff = 0.906172; errup = 0.000573499; errdown = 0.000576648;}
//    else if (pt > 40 && pt < 50) {eff = 0.919065; errup = 0.000362995; errdown = 0.000364482;}
//    else if (pt > 50 && pt < 80) {eff = 0.926712; errup = 0.000740108; errdown = 0.000747018;}
//    else if (pt > 80 && pt < 120) {eff = 0.932497; errup = 0.00248923; errdown = 0.00257564;}
//    else if (pt > 120 && pt < 9999) {eff = 0.951389; errup = 0.00366803; errdown = 0.00393955;}
//  } else if (abseta > 0.8 && abseta < 1.4442) {
//    if (pt > 7 && pt < 10) {eff = 0; errup = 0.000780782; errdown = 0;}
//    else if (pt > 10 && pt < 12) {eff = 0; errup = 0.000836478; errdown = 0;}
//    else if (pt > 12 && pt < 15) {eff = 0.000194099; errup = 0.000446194; errdown = 0.000160569;}
//    else if (pt > 15 && pt < 20) {eff = 0.000245263; errup = 0.000193884; errdown = 0.000117376;}
//    else if (pt > 20 && pt < 25) {eff = 0.00369537; errup = 0.000376959; errdown = 0.00034357;}
//    else if (pt > 25 && pt < 27) {eff = 0.0186315; errup = 0.0010476; errdown = 0.000994197;}
//    else if (pt > 27 && pt < 30) {eff = 0.336783; errup = 0.00246486; errdown = 0.00245606;}
//    else if (pt > 30 && pt < 35) {eff = 0.772262; errup = 0.00132153; errdown = 0.00132694;}
//    else if (pt > 35 && pt < 40) {eff = 0.911008; errup = 0.000670697; errdown = 0.00067527;}
//    else if (pt > 40 && pt < 50) {eff = 0.928774; errup = 0.000414421; errdown = 0.000416653;}
//    else if (pt > 50 && pt < 80) {eff = 0.934931; errup = 0.000845724; errdown = 0.000856004;}
//    else if (pt > 80 && pt < 120) {eff = 0.945397; errup = 0.0027228; errdown = 0.00285347;}
//    else if (pt > 120 && pt < 9999) {eff = 0.96668; errup = 0.00362313; errdown = 0.00402363;}
//  } else if (abseta > 1.4442 && abseta < 1.566) {
//    if (pt > 7 && pt < 10) {eff = 0; errup = 0.00648426; errdown = 0;}
//    else if (pt > 10 && pt < 12) {eff = 0; errup = 0.00708299; errdown = 0;}
//    else if (pt > 12 && pt < 15) {eff = 0; errup = 0.00286353; errdown = 0;}
//    else if (pt > 15 && pt < 20) {eff = 0.00419507; errup = 0.00206229; errdown = 0.00145042;}
//    else if (pt > 20 && pt < 25) {eff = 0.0555069; errup = 0.00419511; errdown = 0.00392574;}
//    else if (pt > 25 && pt < 27) {eff = 0.155746; errup = 0.00857083; errdown = 0.00821316;}
//    else if (pt > 27 && pt < 30) {eff = 0.389774; errup = 0.0079737; errdown = 0.00791614;}
//    else if (pt > 30 && pt < 35) {eff = 0.737773; errup = 0.00428998; errdown = 0.00433519;}
//    else if (pt > 35 && pt < 40) {eff = 0.856978; errup = 0.0025992; errdown = 0.00263873;}
//    else if (pt > 40 && pt < 50) {eff = 0.896748; errup = 0.00162101; errdown = 0.00164365;}
//    else if (pt > 50 && pt < 80) {eff = 0.882448; errup = 0.00347724; errdown = 0.00356721;}
//    else if (pt > 80 && pt < 120) {eff = 0.909605; errup = 0.0108685; errdown = 0.0121005;}
//    else if (pt > 120 && pt < 9999) {eff = 0.946996; errup = 0.0133303; errdown = 0.0169006;}
//  } else if (abseta > 1.566 && abseta < 2.5) {
//    if (pt > 7 && pt < 10) {eff = 0; errup = 0.000574795; errdown = 0;}
//    else if (pt > 10 && pt < 12) {eff = 0; errup = 0.000574615; errdown = 0;}
//    else if (pt > 12 && pt < 15) {eff = 0; errup = 0.000265013; errdown = 0;}
//    else if (pt > 15 && pt < 20) {eff = 0.000101041; errup = 0.00013325; errdown = 6.52627e-05;}
//    else if (pt > 20 && pt < 25) {eff = 0.00215234; errup = 0.000284377; errdown = 0.000252853;}
//    else if (pt > 25 && pt < 27) {eff = 0.0141624; errup = 0.000910277; errdown = 0.000857548;}
//    else if (pt > 27 && pt < 30) {eff = 0.224702; errup = 0.0022028; errdown = 0.00218766;}
//    else if (pt > 30 && pt < 35) {eff = 0.640574; errup = 0.00157228; errdown = 0.0015753;}
//    else if (pt > 35 && pt < 40) {eff = 0.889546; errup = 0.000803088; errdown = 0.000808213;}
//    else if (pt > 40 && pt < 50) {eff = 0.936953; errup = 0.000424953; errdown = 0.000427632;}
//    else if (pt > 50 && pt < 80) {eff = 0.949704; errup = 0.000809845; errdown = 0.000822279;}
//    else if (pt > 80 && pt < 120) {eff = 0.952719; errup = 0.00290805; errdown = 0.0030827;}
//    else if (pt > 120 && pt < 9999) {eff = 0.954786; errup = 0.00504891; errdown = 0.00561401;}
//  }
//  return std::vector<float>({eff, errup, errdown});
//}
//
//std::vector<float> eff_dielleg12(float pt, float abseta) {
//  float errup=0.0, errdown=0.0, eff=1.0;
//  if (abseta > 0 && abseta < 0.8) {
//    if (pt > 7 && pt < 10) {eff = 0.012945; errup = 0.00269541; errdown = 0.00226632;}
//    else if (pt > 10 && pt < 12) {eff = 0.0253289; errup = 0.00318164; errdown = 0.00285136;}
//    else if (pt > 12 && pt < 15) {eff = 0.626014; errup = 0.0055386; errdown = 0.00557143;}
//    else if (pt > 15 && pt < 20) {eff = 0.783797; errup = 0.00262011; errdown = 0.00264313;}
//    else if (pt > 20 && pt < 25) {eff = 0.816201; errup = 0.00177857; errdown = 0.00179192;}
//    else if (pt > 25 && pt < 27) {eff = 0.832487; errup = 0.00223836; errdown = 0.00226231;}
//    else if (pt > 27 && pt < 30) {eff = 0.838534; errup = 0.00156647; errdown = 0.00157877;}
//    else if (pt > 30 && pt < 35) {eff = 0.860641; errup = 0.000903503; errdown = 0.00090842;}
//    else if (pt > 35 && pt < 40) {eff = 0.881148; errup = 0.000671396; errdown = 0.000674683;}
//    else if (pt > 40 && pt < 50) {eff = 0.889903; errup = 0.000611071; errdown = 0.000614048;}
//    else if (pt > 50 && pt < 80) {eff = 0.904008; errup = 0.00259569; errdown = 0.00265902;}
//    else if (pt > 80 && pt < 120) {eff = 0.903786; errup = 0.00739647; errdown = 0.00791876;}
//    else if (pt > 120 && pt < 9999) {eff = 0.922636; errup = 0.0144127; errdown = 0.0170728;}
//  } else if (abseta > 0.8 && abseta < 1.4442) {
//    if (pt > 7 && pt < 10) {eff = 0.0144251; errup = 0.002897; errdown = 0.00244941;}
//    else if (pt > 10 && pt < 12) {eff = 0.0259091; errup = 0.00384786; errdown = 0.00338801;}
//    else if (pt > 12 && pt < 15) {eff = 0.494953; errup = 0.00706297; errdown = 0.00706098;}
//    else if (pt > 15 && pt < 20) {eff = 0.751793; errup = 0.00339753; errdown = 0.00342869;}
//    else if (pt > 20 && pt < 25) {eff = 0.812339; errup = 0.00221921; errdown = 0.00223943;}
//    else if (pt > 25 && pt < 27) {eff = 0.827707; errup = 0.0027844; errdown = 0.00282014;}
//    else if (pt > 27 && pt < 30) {eff = 0.843696; errup = 0.00188443; errdown = 0.00190299;}
//    else if (pt > 30 && pt < 35) {eff = 0.862414; errup = 0.00108516; errdown = 0.00109237;}
//    else if (pt > 35 && pt < 40) {eff = 0.883127; errup = 0.000794395; errdown = 0.00079909;}
//    else if (pt > 40 && pt < 50) {eff = 0.897446; errup = 0.000725904; errdown = 0.000730465;}
//    else if (pt > 50 && pt < 80) {eff = 0.903024; errup = 0.00313616; errdown = 0.00322771;}
//    else if (pt > 80 && pt < 120) {eff = 0.928967; errup = 0.00784714; errdown = 0.00868498;}
//    else if (pt > 120 && pt < 9999) {eff = 0.958525; errup = 0.0134374; errdown = 0.0183674;}
//  } else if (abseta > 1.4442 && abseta < 1.566) {
//    if (pt > 7 && pt < 10) {eff = 0.0106007; errup = 0.0102039; errdown = 0.00576385;}
//    else if (pt > 10 && pt < 12) {eff = 0.138996; errup = 0.0249223; errdown = 0.0218915;}
//    else if (pt > 12 && pt < 15) {eff = 0.394081; errup = 0.0202201; errdown = 0.0198771;}
//    else if (pt > 15 && pt < 20) {eff = 0.677504; errup = 0.0108681; errdown = 0.0110588;}
//    else if (pt > 20 && pt < 25) {eff = 0.783815; errup = 0.00709654; errdown = 0.00726578;}
//    else if (pt > 25 && pt < 27) {eff = 0.813004; errup = 0.00884207; errdown = 0.00916659;}
//    else if (pt > 27 && pt < 30) {eff = 0.81372; errup = 0.00628643; errdown = 0.00645094;}
//    else if (pt > 30 && pt < 35) {eff = 0.839259; errup = 0.00357671; errdown = 0.00364136;}
//    else if (pt > 35 && pt < 40) {eff = 0.860535; errup = 0.00266377; errdown = 0.00270662;}
//    else if (pt > 40 && pt < 50) {eff = 0.880162; errup = 0.00247658; errdown = 0.00252106;}
//    else if (pt > 50 && pt < 80) {eff = 0.884615; errup = 0.0108525; errdown = 0.0117679;}
//    else if (pt > 80 && pt < 120) {eff = 0.902913; errup = 0.029471; errdown = 0.0385366;}
//    else if (pt > 120 && pt < 9999) {eff = 0.904762; errup = 0.0612701; errdown = 0.112063;}
//  } else if (abseta > 1.566 && abseta < 2.5) {
//    if (pt > 7 && pt < 10) {eff = 0.0109307; errup = 0.00216322; errdown = 0.00183179;}
//    else if (pt > 10 && pt < 12) {eff = 0.0184202; errup = 0.00269498; errdown = 0.00237454;}
//    else if (pt > 12 && pt < 15) {eff = 0.384394; errup = 0.00592518; errdown = 0.00589149;}
//    else if (pt > 15 && pt < 20) {eff = 0.710316; errup = 0.00323867; errdown = 0.00326009;}
//    else if (pt > 20 && pt < 25) {eff = 0.817141; errup = 0.00211885; errdown = 0.00213794;}
//    else if (pt > 25 && pt < 27) {eff = 0.850584; errup = 0.00259437; errdown = 0.00263166;}
//    else if (pt > 27 && pt < 30) {eff = 0.874331; errup = 0.00173591; errdown = 0.00175652;}
//    else if (pt > 30 && pt < 35) {eff = 0.895512; errup = 0.00100114; errdown = 0.00100964;}
//    else if (pt > 35 && pt < 40) {eff = 0.91735; errup = 0.000737991; errdown = 0.000744006;}
//    else if (pt > 40 && pt < 50) {eff = 0.931165; errup = 0.000624937; errdown = 0.00063021;}
//    else if (pt > 50 && pt < 80) {eff = 0.939379; errup = 0.00268275; errdown = 0.00279587;}
//    else if (pt > 80 && pt < 120) {eff = 0.941463; errup = 0.00823565; errdown = 0.00939174;}
//    else if (pt > 120 && pt < 9999) {eff = 0.95935; errup = 0.0174478; errdown = 0.0265636;}
//  }
//  return std::vector<float>({eff, errup, errdown});
//}
//
//std::vector<float> eff_muhltpt17(float pt, float abseta) {
//  float errup=0.0, errdown=0.0, eff=1.0;
//  if (abseta > 0 && abseta < 0.9) {
//    if (pt > 5 && pt < 7) {eff = 0; errup = 0.00354096; errdown = 0;}
//    else if (pt > 7 && pt < 10) {eff = 0; errup = 0.00071332; errdown = 0;}
//    else if (pt > 10 && pt < 15) {eff = 0; errup = 0.000112691; errdown = 0;}
//    else if (pt > 15 && pt < 20) {eff = 0.00389236; errup = 0.000318329; errdown = 0.000295153;}
//    else if (pt > 20 && pt < 23) {eff = 0.301194; errup = 0.00205709; errdown = 0.00204917;}
//    else if (pt > 23 && pt < 25) {eff = 0.537094; errup = 0.00225954; errdown = 0.00226105;}
//    else if (pt > 25 && pt < 30) {eff = 0.82903; errup = 0.000833374; errdown = 0.000836601;}
//    else if (pt > 30 && pt < 40) {eff = 0.866481; errup = 0.00036174; errdown = 0.000362569;}
//    else if (pt > 40 && pt < 50) {eff = 0.887467; errup = 0.000309959; errdown = 0.000310705;}
//    else if (pt > 50 && pt < 60) {eff = 0.902724; errup = 0.000609708; errdown = 0.000613125;}
//    else if (pt > 60 && pt < 120) {eff = 0.901655; errup = 0.000855842; errdown = 0.000862497;}
//    else if (pt > 120 && pt < 9999) {eff = 0.888141; errup = 0.00356491; errdown = 0.00366523;}
//  } else if (abseta > 0.9 && abseta < 1.2) {
//    if (pt > 5 && pt < 7) {eff = 0; errup = 0.00363177; errdown = 0;}
//    else if (pt > 7 && pt < 10) {eff = 0; errup = 0.00104076; errdown = 0;}
//    else if (pt > 10 && pt < 15) {eff = 0; errup = 0.000261698; errdown = 0;}
//    else if (pt > 15 && pt < 20) {eff = 0.00871875; errup = 0.000836734; errdown = 0.000766833;}
//    else if (pt > 20 && pt < 23) {eff = 0.295943; errup = 0.00381382; errdown = 0.00378579;}
//    else if (pt > 23 && pt < 25) {eff = 0.538557; errup = 0.00430175; errdown = 0.00430746;}
//    else if (pt > 25 && pt < 30) {eff = 0.811972; errup = 0.00169079; errdown = 0.00170249;}
//    else if (pt > 30 && pt < 40) {eff = 0.862518; errup = 0.000704098; errdown = 0.000707133;}
//    else if (pt > 40 && pt < 50) {eff = 0.887386; errup = 0.000561027; errdown = 0.000563471;}
//    else if (pt > 50 && pt < 60) {eff = 0.90209; errup = 0.00110518; errdown = 0.00111634;}
//    else if (pt > 60 && pt < 120) {eff = 0.897924; errup = 0.00158096; errdown = 0.00160277;}
//    else if (pt > 120 && pt < 9999) {eff = 0.887084; errup = 0.00641928; errdown = 0.00674366;}
//  } else if (abseta > 1.2 && abseta < 2.1) {
//    if (pt > 5 && pt < 7) {eff = 0; errup = 0.000609222; errdown = 0;}
//    else if (pt > 7 && pt < 10) {eff = 0; errup = 0.000264594; errdown = 0;}
//    else if (pt > 10 && pt < 15) {eff = 0.00013693; errup = 0.000133178; errdown = 7.45212e-05;}
//    else if (pt > 15 && pt < 20) {eff = 0.0124952; errup = 0.000570298; errdown = 0.000546212;}
//    else if (pt > 20 && pt < 23) {eff = 0.251228; errup = 0.00219043; errdown = 0.00217788;}
//    else if (pt > 23 && pt < 25) {eff = 0.513877; errup = 0.00262943; errdown = 0.00263019;}
//    else if (pt > 25 && pt < 30) {eff = 0.780436; errup = 0.00112003; errdown = 0.00112414;}
//    else if (pt > 30 && pt < 40) {eff = 0.830997; errup = 0.000511629; errdown = 0.000512864;}
//    else if (pt > 40 && pt < 50) {eff = 0.858819; errup = 0.000410034; errdown = 0.00041103;}
//    else if (pt > 50 && pt < 60) {eff = 0.874947; errup = 0.00081753; errdown = 0.000822119;}
//    else if (pt > 60 && pt < 120) {eff = 0.879603; errup = 0.0011485; errdown = 0.00115798;}
//    else if (pt > 120 && pt < 9999) {eff = 0.876618; errup = 0.00485486; errdown = 0.00502094;}
//  } else if (abseta > 2.1 && abseta < 2.4) {
//    if (pt > 5 && pt < 7) {eff = 0; errup = 0.0015696; errdown = 0;}
//    else if (pt > 7 && pt < 10) {eff = 0; errup = 0.000737612; errdown = 0;}
//    else if (pt > 10 && pt < 15) {eff = 0.000717669; errup = 0.000485207; errdown = 0.000309955;}
//    else if (pt > 15 && pt < 20) {eff = 0.0117705; errup = 0.00105803; errdown = 0.000974872;}
//    else if (pt > 20 && pt < 23) {eff = 0.107686; errup = 0.00299285; errdown = 0.00292249;}
//    else if (pt > 23 && pt < 25) {eff = 0.324534; errup = 0.00470441; errdown = 0.00466961;}
//    else if (pt > 25 && pt < 30) {eff = 0.593711; errup = 0.00252587; errdown = 0.00253081;}
//    else if (pt > 30 && pt < 40) {eff = 0.705915; errup = 0.00123867; errdown = 0.00124172;}
//    else if (pt > 40 && pt < 50) {eff = 0.751504; errup = 0.00114816; errdown = 0.00115171;}
//    else if (pt > 50 && pt < 60) {eff = 0.769389; errup = 0.00235506; errdown = 0.00237191;}
//    else if (pt > 60 && pt < 120) {eff = 0.781308; errup = 0.00344268; errdown = 0.00348176;}
//    else if (pt > 120 && pt < 9999) {eff = 0.795699; errup = 0.0160855; errdown = 0.0170372;}
//  }
//  return std::vector<float>({eff, errup, errdown});
//}
//
//std::vector<float> eff_dimuleg8(float pt, float abseta) {
//  float errup=0.0, errdown=0.0, eff=1.0;
//  if (abseta > 0 && abseta < 0.9) {
//    if (pt > 5 && pt < 7) {eff = 0.00385356; errup = 0.00505989; errdown = 0.00248866;}
//    else if (pt > 7 && pt < 10) {eff = 0.698062; errup = 0.00915229; errdown = 0.00930908;}
//    else if (pt > 10 && pt < 15) {eff = 0.886263; errup = 0.00249071; errdown = 0.00253858;}
//    else if (pt > 15 && pt < 20) {eff = 0.906415; errup = 0.00138361; errdown = 0.00140205;}
//    else if (pt > 20 && pt < 23) {eff = 0.912981; errup = 0.00125743; errdown = 0.00127396;}
//    else if (pt > 23 && pt < 25) {eff = 0.911771; errup = 0.00128278; errdown = 0.00129971;}
//    else if (pt > 25 && pt < 30) {eff = 0.912375; errup = 0.000641682; errdown = 0.000645941;}
//    else if (pt > 30 && pt < 40) {eff = 0.911811; errup = 0.000348166; errdown = 0.000349409;}
//    else if (pt > 40 && pt < 50) {eff = 0.910492; errup = 0.000441784; errdown = 0.000443754;}
//    else if (pt > 50 && pt < 60) {eff = 0.912906; errup = 0.00251545; errdown = 0.00258187;}
//    else if (pt > 60 && pt < 120) {eff = 0.900718; errup = 0.00333859; errdown = 0.00343965;}
//    else if (pt > 120 && pt < 9999) {eff = 0.85439; errup = 0.0165766; errdown = 0.0181918;}
//  } else if (abseta > 0.9 && abseta < 1.2) {
//    if (pt > 5 && pt < 7) {eff = 0.00197628; errup = 0.00452972; errdown = 0.00163493;}
//    else if (pt > 7 && pt < 10) {eff = 0.697398; errup = 0.0110908; errdown = 0.0113198;}
//    else if (pt > 10 && pt < 15) {eff = 0.901621; errup = 0.00356399; errdown = 0.00368046;}
//    else if (pt > 15 && pt < 20) {eff = 0.911909; errup = 0.00234504; errdown = 0.00240199;}
//    else if (pt > 20 && pt < 23) {eff = 0.922192; errup = 0.002215; errdown = 0.00227337;}
//    else if (pt > 23 && pt < 25) {eff = 0.923252; errup = 0.00228533; errdown = 0.00234844;}
//    else if (pt > 25 && pt < 30) {eff = 0.923627; errup = 0.00117842; errdown = 0.0011952;}
//    else if (pt > 30 && pt < 40) {eff = 0.927034; errup = 0.000601698; errdown = 0.000606283;}
//    else if (pt > 40 && pt < 50) {eff = 0.929016; errup = 0.000711366; errdown = 0.000717976;}
//    else if (pt > 50 && pt < 60) {eff = 0.929361; errup = 0.00422294; errdown = 0.00446185;}
//    else if (pt > 60 && pt < 120) {eff = 0.921064; errup = 0.00552283; errdown = 0.00588584;}
//    else if (pt > 120 && pt < 9999) {eff = 0.903704; errup = 0.0256738; errdown = 0.0324938;}
//  } else if (abseta > 1.2 && abseta < 2.1) {
//    if (pt > 5 && pt < 7) {eff = 0.00662032; errup = 0.00182822; errdown = 0.00146494;}
//    else if (pt > 7 && pt < 10) {eff = 0.623401; errup = 0.00586265; errdown = 0.00589856;}
//    else if (pt > 10 && pt < 15) {eff = 0.891232; errup = 0.00210824; errdown = 0.00214434;}
//    else if (pt > 15 && pt < 20) {eff = 0.913543; errup = 0.00138369; errdown = 0.00140385;}
//    else if (pt > 20 && pt < 23) {eff = 0.917482; errup = 0.00137978; errdown = 0.00140091;}
//    else if (pt > 23 && pt < 25) {eff = 0.917532; errup = 0.00144223; errdown = 0.00146534;}
//    else if (pt > 25 && pt < 30) {eff = 0.919152; errup = 0.000758192; errdown = 0.000764699;}
//    else if (pt > 30 && pt < 40) {eff = 0.923075; errup = 0.000416741; errdown = 0.000418815;}
//    else if (pt > 40 && pt < 50) {eff = 0.925559; errup = 0.000474368; errdown = 0.000477155;}
//    else if (pt > 50 && pt < 60) {eff = 0.925551; errup = 0.00288916; errdown = 0.00299383;}
//    else if (pt > 60 && pt < 120) {eff = 0.925987; errup = 0.00375637; errdown = 0.0039353;}
//    else if (pt > 120 && pt < 9999) {eff = 0.902622; errup = 0.0183565; errdown = 0.0216761;}
//  } else if (abseta > 2.1 && abseta < 2.4) {
//    if (pt > 5 && pt < 7) {eff = 0.00426621; errup = 0.00287575; errdown = 0.00184154;}
//    else if (pt > 7 && pt < 10) {eff = 0.566333; errup = 0.0100928; errdown = 0.010147;}
//    else if (pt > 10 && pt < 15) {eff = 0.869241; errup = 0.00405685; errdown = 0.00416472;}
//    else if (pt > 15 && pt < 20) {eff = 0.884084; errup = 0.00290338; errdown = 0.00296706;}
//    else if (pt > 20 && pt < 23) {eff = 0.892403; errup = 0.00292143; errdown = 0.00299181;}
//    else if (pt > 23 && pt < 25) {eff = 0.896152; errup = 0.00303193; errdown = 0.00311095;}
//    else if (pt > 25 && pt < 30) {eff = 0.896845; errup = 0.0016064; errdown = 0.00162865;}
//    else if (pt > 30 && pt < 40) {eff = 0.899698; errup = 0.000961904; errdown = 0.000970127;}
//    else if (pt > 40 && pt < 50) {eff = 0.901419; errup = 0.00126342; errdown = 0.00127791;}
//    else if (pt > 50 && pt < 60) {eff = 0.907641; errup = 0.00751323; errdown = 0.00807885;}
//    else if (pt > 60 && pt < 120) {eff = 0.901299; errup = 0.0108462; errdown = 0.0119507;}
//    else if (pt > 120 && pt < 9999) {eff = 0.941176; errup = 0.048713; errdown = 0.12258;}
//  }
//  return std::vector<float>({eff, errup, errdown});
//}
//
//std::vector<float> eff_isoel3235(float pt, float abseta) {
//  float errup=0.0, errdown=0.0, eff=1.0;
//  if (abseta > 0 && abseta < 0.8) {
//    if (pt > 7 && pt < 10) {eff = 0.00140449; errup = 0.00322219; errdown = 0.00116189;}
//    else if (pt > 10 && pt < 12) {eff = 0; errup = 0.00436343; errdown = 0;}
//    else if (pt > 12 && pt < 15) {eff = 0.00184162; errup = 0.004222; errdown = 0.00152352;}
//    else if (pt > 15 && pt < 20) {eff = 0.00156986; errup = 0.00360059; errdown = 0.0012987;}
//    else if (pt > 20 && pt < 25) {eff = 0; errup = 0.00357535; errdown = 0;}
//    else if (pt > 25 && pt < 27) {eff = 0; errup = 0.00984914; errdown = 0;}
//    else if (pt > 27 && pt < 30) {eff = 0; errup = 0.00713791; errdown = 0;}
//    else if (pt > 30 && pt < 35) {eff = 0.312676; errup = 0.0265349; errdown = 0.0254211;}
//    else if (pt > 35 && pt < 40) {eff = 0.744615; errup = 0.0248944; errdown = 0.0264924;}
//    else if (pt > 40 && pt < 50) {eff = 0.753176; errup = 0.0187794; errdown = 0.0197427;}
//    else if (pt > 50 && pt < 80) {eff = 0.79952; errup = 0.0140911; errdown = 0.0148408;}
//    else if (pt > 80 && pt < 120) {eff = 0.821168; errup = 0.0166547; errdown = 0.0178908;}
//    else if (pt > 120 && pt < 9999) {eff = 0.881057; errup = 0.0153914; errdown = 0.0171894;}
//  } else if (abseta > 0.8 && abseta < 1.4442) {
//    if (pt > 7 && pt < 10) {eff = 0; errup = 0.00385215; errdown = 0;}
//    else if (pt > 10 && pt < 12) {eff = 0.00381679; errup = 0.00872168; errdown = 0.00315765;}
//    else if (pt > 12 && pt < 15) {eff = 0; errup = 0.00549688; errdown = 0;}
//    else if (pt > 15 && pt < 20) {eff = 0; errup = 0.00418489; errdown = 0;}
//    else if (pt > 20 && pt < 25) {eff = 0; errup = 0.00459198; errdown = 0;}
//    else if (pt > 25 && pt < 27) {eff = 0; errup = 0.0170586; errdown = 0;}
//    else if (pt > 27 && pt < 30) {eff = 0; errup = 0.0107085; errdown = 0;}
//    else if (pt > 30 && pt < 35) {eff = 0.197026; errup = 0.027285; errdown = 0.0248654;}
//    else if (pt > 35 && pt < 40) {eff = 0.628205; errup = 0.033073; errdown = 0.0342397;}
//    else if (pt > 40 && pt < 50) {eff = 0.750716; errup = 0.0237977; errdown = 0.0253206;}
//    else if (pt > 50 && pt < 80) {eff = 0.781991; errup = 0.0167243; errdown = 0.0176572;}
//    else if (pt > 80 && pt < 120) {eff = 0.80826; errup = 0.0218627; errdown = 0.0238032;}
//    else if (pt > 120 && pt < 9999) {eff = 0.863333; errup = 0.0201738; errdown = 0.0227942;}
//  } else if (abseta > 1.4442 && abseta < 1.566) {
//    if (pt > 7 && pt < 10) {eff = 0; errup = 0.0233265; errdown = 0;}
//    else if (pt > 10 && pt < 12) {eff = 0; errup = 0.0409781; errdown = 0;}
//    else if (pt > 12 && pt < 15) {eff = 0; errup = 0.0329191; errdown = 0;}
//    else if (pt > 15 && pt < 20) {eff = 0; errup = 0.0230346; errdown = 0;}
//    else if (pt > 20 && pt < 25) {eff = 0; errup = 0.0368748; errdown = 0;}
//    else if (pt > 25 && pt < 27) {eff = 0; errup = 0.0923495; errdown = 0;}
//    else if (pt > 27 && pt < 30) {eff = 0.0416667; errup = 0.0893855; errdown = 0.0344944;}
//    else if (pt > 30 && pt < 35) {eff = 0.25; errup = 0.108468; errdown = 0.0872543;}
//    else if (pt > 35 && pt < 40) {eff = 0.518519; errup = 0.1114; errdown = 0.112994;}
//    else if (pt > 40 && pt < 50) {eff = 0.552632; errup = 0.0908274; errdown = 0.0939923;}
//    else if (pt > 50 && pt < 80) {eff = 0.707692; errup = 0.0600147; errdown = 0.0672026;}
//    else if (pt > 80 && pt < 120) {eff = 0.777778; errup = 0.0729951; errdown = 0.0911652;}
//    else if (pt > 120 && pt < 9999) {eff = 0.681818; errup = 0.109079; errdown = 0.12875;}
//  } else if (abseta > 1.566 && abseta < 2.5) {
//    if (pt > 7 && pt < 10) {eff = 0; errup = 0.00315827; errdown = 0;}
//    else if (pt > 10 && pt < 12) {eff = 0.0031348; errup = 0.00717134; errdown = 0.00259339;}
//    else if (pt > 12 && pt < 15) {eff = 0; errup = 0.00496338; errdown = 0;}
//    else if (pt > 15 && pt < 20) {eff = 0; errup = 0.00378094; errdown = 0;}
//    else if (pt > 20 && pt < 25) {eff = 0; errup = 0.00507278; errdown = 0;}
//    else if (pt > 25 && pt < 27) {eff = 0; errup = 0.0133482; errdown = 0;}
//    else if (pt > 27 && pt < 30) {eff = 0.00537634; errup = 0.0122537; errdown = 0.00444799;}
//    else if (pt > 30 && pt < 35) {eff = 0.0983607; errup = 0.0229311; errdown = 0.0192944;}
//    else if (pt > 35 && pt < 40) {eff = 0.5; errup = 0.0371673; errdown = 0.0371673;}
//    else if (pt > 40 && pt < 50) {eff = 0.656627; errup = 0.0270269; errdown = 0.028023;}
//    else if (pt > 50 && pt < 80) {eff = 0.737705; errup = 0.0192155; errdown = 0.0201225;}
//    else if (pt > 80 && pt < 120) {eff = 0.771127; errup = 0.0256314; errdown = 0.0276711;}
//    else if (pt > 120 && pt < 9999) {eff = 0.795276; errup = 0.0259964; errdown = 0.028496;}
//  }
//  return std::vector<float>({eff, errup, errdown});
//}
//
//std::vector<float> eff_isomu2427(float pt, float abseta) {
//  float errup=0.0, errdown=0.0, eff=1.0;
//  if (abseta > 0 && abseta < 0.9) {
//    if (pt > 5 && pt < 7) {eff = 0.0016; errup = 0.00366954; errdown = 0.00132363;}
//    else if (pt > 7 && pt < 10) {eff = 0.00414938; errup = 0.00401971; errdown = 0.00225741;}
//    else if (pt > 10 && pt < 15) {eff = 0.00178253; errup = 0.00234618; errdown = 0.00115127;}
//    else if (pt > 15 && pt < 20) {eff = 0.00117371; errup = 0.00269375; errdown = 0.000970967;}
//    else if (pt > 20 && pt < 23) {eff = 0.0021645; errup = 0.00495958; errdown = 0.00179065;}
//    else if (pt > 23 && pt < 25) {eff = 0.292683; errup = 0.0318821; errdown = 0.0300831;}
//    else if (pt > 25 && pt < 30) {eff = 0.781513; errup = 0.0172733; errdown = 0.0182655;}
//    else if (pt > 30 && pt < 40) {eff = 0.830366; errup = 0.0116332; errdown = 0.0122773;}
//    else if (pt > 40 && pt < 50) {eff = 0.861438; errup = 0.0126449; errdown = 0.0136385;}
//    else if (pt > 50 && pt < 60) {eff = 0.866559; errup = 0.0138113; errdown = 0.0150583;}
//    else if (pt > 60 && pt < 120) {eff = 0.857143; errup = 0.0091918; errdown = 0.0096932;}
//    else if (pt > 120 && pt < 9999) {eff = 0.808929; errup = 0.0169152; errdown = 0.0180766;}
//  } else if (abseta > 0.9 && abseta < 1.2) {
//    if (pt > 5 && pt < 7) {eff = 0; errup = 0.0115844; errdown = 0;}
//    else if (pt > 7 && pt < 10) {eff = 0; errup = 0.0082586; errdown = 0;}
//    else if (pt > 10 && pt < 15) {eff = 0; errup = 0.0055299; errdown = 0;}
//    else if (pt > 15 && pt < 20) {eff = 0.00416667; errup = 0.00951566; errdown = 0.00344712;}
//    else if (pt > 20 && pt < 23) {eff = 0.00714286; errup = 0.0162324; errdown = 0.00590966;}
//    else if (pt > 23 && pt < 25) {eff = 0.228571; errup = 0.0613516; errdown = 0.0525767;}
//    else if (pt > 25 && pt < 30) {eff = 0.698925; errup = 0.0350786; errdown = 0.0373806;}
//    else if (pt > 30 && pt < 40) {eff = 0.772727; errup = 0.0254723; errdown = 0.0275096;}
//    else if (pt > 40 && pt < 50) {eff = 0.816964; errup = 0.026494; errdown = 0.0295599;}
//    else if (pt > 50 && pt < 60) {eff = 0.819672; errup = 0.0291864; errdown = 0.0330027;}
//    else if (pt > 60 && pt < 120) {eff = 0.851685; errup = 0.0171206; errdown = 0.0188043;}
//    else if (pt > 120 && pt < 9999) {eff = 0.845588; errup = 0.0317581; errdown = 0.0374164;}
//  } else if (abseta > 1.2 && abseta < 2.1) {
//    if (pt > 5 && pt < 7) {eff = 0; errup = 0.00350057; errdown = 0;}
//    else if (pt > 7 && pt < 10) {eff = 0; errup = 0.00291799; errdown = 0;}
//    else if (pt > 10 && pt < 15) {eff = 0.00235018; errup = 0.00309128; errdown = 0.00151785;}
//    else if (pt > 15 && pt < 20) {eff = 0.0030722; errup = 0.00403759; errdown = 0.00198411;}
//    else if (pt > 20 && pt < 23) {eff = 0.00310559; errup = 0.00710486; errdown = 0.00256923;}
//    else if (pt > 23 && pt < 25) {eff = 0.237885; errup = 0.0316479; errdown = 0.0291654;}
//    else if (pt > 25 && pt < 30) {eff = 0.651584; errup = 0.0234085; errdown = 0.0241276;}
//    else if (pt > 30 && pt < 40) {eff = 0.76738; errup = 0.0157356; errdown = 0.0164808;}
//    else if (pt > 40 && pt < 50) {eff = 0.790927; errup = 0.0184278; errdown = 0.0196367;}
//    else if (pt > 50 && pt < 60) {eff = 0.826316; errup = 0.0198178; errdown = 0.0216486;}
//    else if (pt > 60 && pt < 120) {eff = 0.823654; errup = 0.0125525; errdown = 0.0132647;}
//    else if (pt > 120 && pt < 9999) {eff = 0.830325; errup = 0.0230458; errdown = 0.025616;}
//  } else if (abseta > 2.1 && abseta < 2.4) {
//    if (pt > 5 && pt < 7) {eff = 0.0125; errup = 0.0162484; errdown = 0.00806978;}
//    else if (pt > 7 && pt < 10) {eff = 0; errup = 0.0107085; errdown = 0;}
//    else if (pt > 10 && pt < 15) {eff = 0; errup = 0.00818516; errdown = 0;}
//    else if (pt > 15 && pt < 20) {eff = 0; errup = 0.0105853; errdown = 0;}
//    else if (pt > 20 && pt < 23) {eff = 0; errup = 0.0219368; errdown = 0;}
//    else if (pt > 23 && pt < 25) {eff = 0.148936; errup = 0.0709778; errdown = 0.0532465;}
//    else if (pt > 25 && pt < 30) {eff = 0.447368; errup = 0.0512155; errdown = 0.0502109;}
//    else if (pt > 30 && pt < 40) {eff = 0.673077; errup = 0.0394325; errdown = 0.041831;}
//    else if (pt > 40 && pt < 50) {eff = 0.734694; errup = 0.0468096; errdown = 0.0521118;}
//    else if (pt > 50 && pt < 60) {eff = 0.696203; errup = 0.0549362; errdown = 0.0604646;}
//    else if (pt > 60 && pt < 120) {eff = 0.725581; errup = 0.0315616; errdown = 0.0338138;}
//    else if (pt > 120 && pt < 9999) {eff = 0.625; errup = 0.0842352; errdown = 0.0913836;}
//  }
//  return std::vector<float>({eff, errup, errdown});
//}

std::vector<float> eff_elhltpt23(float pt, float abseta) {
  float errup=0.0, errdown=0.0, eff=1.0;
  if (abseta > 0 && abseta < 0.8) {
    if (pt > 7 && pt < 15) {eff = 0.0304; errup = 0.016; errdown = 0.011;}
    else if (pt > 15 && pt < 20) {eff = 0.0304348; errup = 0.0160084; errdown = 0.011157;}
    else if (pt > 20 && pt < 25) {eff = 0.341004; errup = 0.0230534; errdown = 0.022357;}
    else if (pt > 25 && pt < 27) {eff = 0.993421; errup = 0.00283884; errdown = 0.00442608;}
    else if (pt > 27 && pt < 30) {eff = 0.999921; errup = 5.10658e-05; errdown = 0.000104266;}
    else if (pt > 30 && pt < 35) {eff = 1; errup = 0; errdown = 1.56117e-05;}
    else if (pt > 35 && pt < 40) {eff = 1; errup = 0; errdown = 7.84081e-06;}
    else if (pt > 40 && pt < 50) {eff = 1; errup = 0; errdown = 3.52132e-06;}
    else if (pt > 50 && pt < 80) {eff = 1; errup = 0; errdown = 1.57584e-05;}
    else if (pt > 80 && pt < 120) {eff = 1; errup = 0; errdown = 0.000191635;}
    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.000554705;}
  } else if (abseta > 0.8 && abseta < 1.4442) {
    if (pt > 7 && pt < 15) {eff = 0.0303; errup = 0.0285; errdown = 0.016;}
    else if (pt > 15 && pt < 20) {eff = 0.030303; errup = 0.028599; errdown = 0.0164473;}
    else if (pt > 20 && pt < 25) {eff = 0.405316; errup = 0.0302401; errdown = 0.0295754;}
    else if (pt > 25 && pt < 27) {eff = 0.989362; errup = 0.0045876; errdown = 0.00713259;}
    else if (pt > 27 && pt < 30) {eff = 0.999771; errup = 0.000124517; errdown = 0.000222508;}
    else if (pt > 30 && pt < 35) {eff = 1; errup = 0; errdown = 2.31433e-05;}
    else if (pt > 35 && pt < 40) {eff = 1; errup = 0; errdown = 1.10004e-05;}
    else if (pt > 40 && pt < 50) {eff = 1; errup = 0; errdown = 5.13872e-06;}
    else if (pt > 50 && pt < 80) {eff = 1; errup = 0; errdown = 2.3365e-05;}
    else if (pt > 80 && pt < 120) {eff = 1; errup = 0; errdown = 0.000274947;}
    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.00077259;}
  } else if (abseta > 1.4442 && abseta < 1.566) {
    if (pt > 7 && pt < 15) {eff = 0.084; errup = 0.103; errdown = 0.084;}
    else if (pt > 15 && pt < 20) {eff = 0.75; errup = 0.103407; errdown = 0.133748;}
    else if (pt > 20 && pt < 25) {eff = 0.956835; errup = 0.0121489; errdown = 0.0159118;}
    else if (pt > 25 && pt < 27) {eff = 0.99569; errup = 0.0027836; errdown = 0.00565666;}
    else if (pt > 27 && pt < 30) {eff = 1; errup = 0; errdown = 0.000910982;}
    else if (pt > 30 && pt < 35) {eff = 1; errup = 0; errdown = 0.000207792;}
    else if (pt > 35 && pt < 40) {eff = 1; errup = 0; errdown = 0.000109215;}
    else if (pt > 40 && pt < 50) {eff = 1; errup = 0; errdown = 6.21192e-05;}
    else if (pt > 50 && pt < 80) {eff = 1; errup = 0; errdown = 0.000294332;}
    else if (pt > 80 && pt < 120) {eff = 1; errup = 0; errdown = 0.0030535;}
    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.00692318;}
  } else if (abseta > 1.566 && abseta < 2.55) {
    if (pt > 7 && pt < 15) {eff = 0.084; errup = 0.038; errdown = 0.0286;}
    else if (pt > 15 && pt < 20) {eff = 0.0842105; errup = 0.0389062; errdown = 0.0286115;}
    else if (pt > 20 && pt < 25) {eff = 0.347107; errup = 0.0332729; errdown = 0.0319275;}
    else if (pt > 25 && pt < 27) {eff = 0.967005; errup = 0.0089534; errdown = 0.0116542;}
    else if (pt > 27 && pt < 30) {eff = 0.999659; errup = 0.000185318; errdown = 0.000331124;}
    else if (pt > 30 && pt < 35) {eff = 1; errup = 0; errdown = 2.98393e-05;}
    else if (pt > 35 && pt < 40) {eff = 1; errup = 0; errdown = 1.31921e-05;}
    else if (pt > 40 && pt < 50) {eff = 1; errup = 0; errdown = 5.9143e-06;}
    else if (pt > 50 && pt < 80) {eff = 1; errup = 0; errdown = 2.76862e-05;}
    else if (pt > 80 && pt < 120) {eff = 1; errup = 0; errdown = 0.000358181;}
    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.00112882;}
  }
  return std::vector<float>({eff, errup, errdown});
}

std::vector<float> eff_dielleg12(float pt, float abseta) {
  float errup=0.0, errdown=0.0, eff=1.0;
  if (abseta > 0 && abseta < 0.8) {
    if (pt > 7 && pt < 12) {eff = 0.0235973; errup = 0.00218208; errdown = 0.00200828;}
    else if (pt > 12 && pt < 15) {eff = 0.618448; errup = 0.00550047; errdown = 0.00553067;}
    else if (pt > 15 && pt < 20) {eff = 0.777711; errup = 0.00263657; errdown = 0.00265892;}
    else if (pt > 20 && pt < 25) {eff = 0.815401; errup = 0.00177492; errdown = 0.00178814;}
    else if (pt > 25 && pt < 27) {eff = 0.82854; errup = 0.00224629; errdown = 0.00226968;}
    else if (pt > 27 && pt < 30) {eff = 0.837537; errup = 0.00156586; errdown = 0.00157805;}
    else if (pt > 30 && pt < 35) {eff = 0.858253; errup = 0.00090979; errdown = 0.000914673;}
    else if (pt > 35 && pt < 40) {eff = 0.880238; errup = 0.000674997; errdown = 0.000678289;}
    else if (pt > 40 && pt < 50) {eff = 0.888684; errup = 0.000613081; errdown = 0.00061604;}
    else if (pt > 50 && pt < 80) {eff = 0.895795; errup = 0.00267898; errdown = 0.00274036;}
    else if (pt > 80 && pt < 120) {eff = 0.894515; errup = 0.00759585; errdown = 0.0080903;}
    else if (pt > 120 && pt < 9999) {eff = 0.923295; errup = 0.0142932; errdown = 0.0169347;}
  } else if (abseta > 0.8 && abseta < 1.4442) {
    if (pt > 7 && pt < 12) {eff = 0.0181859; errup = 0.0022153; errdown = 0.00198969;}
    else if (pt > 12 && pt < 15) {eff = 0.488943; errup = 0.00693829; errdown = 0.00693409;}
    else if (pt > 15 && pt < 20) {eff = 0.755107; errup = 0.00334278; errdown = 0.00337362;}
    else if (pt > 20 && pt < 25) {eff = 0.810641; errup = 0.00221806; errdown = 0.002238;}
    else if (pt > 25 && pt < 27) {eff = 0.827926; errup = 0.0027799; errdown = 0.00281559;}
    else if (pt > 27 && pt < 30) {eff = 0.841489; errup = 0.00188037; errdown = 0.00189852;}
    else if (pt > 30 && pt < 35) {eff = 0.860308; errup = 0.00108443; errdown = 0.0010915;}
    else if (pt > 35 && pt < 40) {eff = 0.882528; errup = 0.00078721; errdown = 0.000791792;}
    else if (pt > 40 && pt < 50) {eff = 0.895206; errup = 0.000723806; errdown = 0.00072823;}
    else if (pt > 50 && pt < 80) {eff = 0.897602; errup = 0.00319002; errdown = 0.00327899;}
    else if (pt > 80 && pt < 120) {eff = 0.923826; errup = 0.00794359; errdown = 0.00873614;}
    else if (pt > 120 && pt < 9999) {eff = 0.949772; errup = 0.0147252; errdown = 0.019454;}
  } else if (abseta > 1.4442 && abseta < 1.566) {
    if (pt > 7 && pt < 12) {eff = 0.0851789; errup = 0.0131347; errdown = 0.0116183;}
    else if (pt > 12 && pt < 15) {eff = 0.4757; errup = 0.0199239; errdown = 0.0198496;}
    else if (pt > 15 && pt < 20) {eff = 0.738636; errup = 0.0103565; errdown = 0.0106215;}
    else if (pt > 20 && pt < 25) {eff = 0.805836; errup = 0.00677726; errdown = 0.00695776;}
    else if (pt > 25 && pt < 27) {eff = 0.83835; errup = 0.00818349; errdown = 0.00852172;}
    else if (pt > 27 && pt < 30) {eff = 0.835424; errup = 0.00585415; errdown = 0.00602265;}
    else if (pt > 30 && pt < 35) {eff = 0.861461; errup = 0.00328862; errdown = 0.00335454;}
    else if (pt > 35 && pt < 40) {eff = 0.879919; errup = 0.00247947; errdown = 0.00252396;}
    else if (pt > 40 && pt < 50) {eff = 0.899791; errup = 0.00251065; errdown = 0.00256703;}
    else if (pt > 50 && pt < 80) {eff = 0.888095; errup = 0.0109847; errdown = 0.0119587;}
    else if (pt > 80 && pt < 120) {eff = 0.883495; errup = 0.0321372; errdown = 0.0406717;}
    else if (pt > 120 && pt < 9999) {eff = 0.869565; errup = 0.0701228; errdown = 0.110814;}
  } else if (abseta > 1.566 && abseta < 2.55) {
    if (pt > 7 && pt < 12) {eff = 0.0152856; errup = 0.00155663; errdown = 0.00142051;}
    else if (pt > 12 && pt < 15) {eff = 0.356028; errup = 0.00552024; errdown = 0.00548267;}
    else if (pt > 15 && pt < 20) {eff = 0.6996; errup = 0.00319689; errdown = 0.00321628;}
    else if (pt > 20 && pt < 25) {eff = 0.80927; errup = 0.0021147; errdown = 0.00213265;}
    else if (pt > 25 && pt < 27) {eff = 0.852071; errup = 0.00254126; errdown = 0.00257749;}
    else if (pt > 27 && pt < 30) {eff = 0.872859; errup = 0.0017212; errdown = 0.00174118;}
    else if (pt > 30 && pt < 35) {eff = 0.894648; errup = 0.00099533; errdown = 0.00100365;}
    else if (pt > 35 && pt < 40) {eff = 0.915473; errup = 0.000738998; errdown = 0.000744881;}
    else if (pt > 40 && pt < 50) {eff = 0.928823; errup = 0.000641385; errdown = 0.000646741;}
    else if (pt > 50 && pt < 80) {eff = 0.935343; errup = 0.00279817; errdown = 0.0029129;}
    else if (pt > 80 && pt < 120) {eff = 0.94492; errup = 0.00801411; errdown = 0.00918603;}
    else if (pt > 120 && pt < 9999) {eff = 0.946429; errup = 0.021046; errdown = 0.0306273;}
  }
  return std::vector<float>({eff, errup, errdown});
}

std::vector<float> eff_muhltpt17(float pt, float abseta) {
  float errup=0.0, errdown=0.0, eff=1.0;
  if (abseta > 0 && abseta < 0.9) {
    if (pt > 5 && pt < 17) {eff = 1; errup = 0; errdown = 0.0101758;}
    else if (pt > 17 && pt < 20) {eff = 1; errup = 0; errdown = 0.000121089;}
    else if (pt > 20 && pt < 23) {eff = 1; errup = 0; errdown = 6.96752e-05;}
    else if (pt > 23 && pt < 25) {eff = 1; errup = 0; errdown = 1.08459e-05;}
    else if (pt > 25 && pt < 30) {eff = 1; errup = 0; errdown = 2.40426e-06;}
    else if (pt > 30 && pt < 40) {eff = 1; errup = 0; errdown = 1.99943e-06;}
    else if (pt > 40 && pt < 50) {eff = 1; errup = 0; errdown = 8.65264e-06;}
    else if (pt > 50 && pt < 60) {eff = 1; errup = 0; errdown = 1.68671e-05;}
    else if (pt > 60 && pt < 120) {eff = 1; errup = 0; errdown = 0.000264404;}
    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.0;}
  } else if (abseta > 0.9 && abseta < 1.2) {
    if (pt > 5 && pt < 17) {eff = 1; errup = 0; errdown = 0.841345;}
    else if (pt > 17 && pt < 20) {eff = 1; errup = 0; errdown = 0.0139553;}
    else if (pt > 20 && pt < 23) {eff = 1; errup = 0; errdown = 0.000421005;}
    else if (pt > 23 && pt < 25) {eff = 1; errup = 0; errdown = 0.000251612;}
    else if (pt > 25 && pt < 30) {eff = 1; errup = 0; errdown = 4.24483e-05;}
    else if (pt > 30 && pt < 40) {eff = 1; errup = 0; errdown = 8.92435e-06;}
    else if (pt > 40 && pt < 50) {eff = 1; errup = 0; errdown = 6.53689e-06;}
    else if (pt > 50 && pt < 60) {eff = 1; errup = 0; errdown = 2.82452e-05;}
    else if (pt > 60 && pt < 120) {eff = 1; errup = 0; errdown = 5.59141e-05;}
    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.000849211;}
  } else if (abseta > 1.2 && abseta < 2.1) {
    if (pt > 5 && pt < 17) {eff = 1; errup = 0; errdown = 0.132046;}
    else if (pt > 17 && pt < 20) {eff = 1; errup = 0; errdown = 0.00365339;}
    else if (pt > 20 && pt < 23) {eff = 1; errup = 0; errdown = 0.000182552;}
    else if (pt > 23 && pt < 25) {eff = 1; errup = 0; errdown = 9.76414e-05;}
    else if (pt > 25 && pt < 30) {eff = 1; errup = 0; errdown = 1.71325e-05;}
    else if (pt > 30 && pt < 40) {eff = 1; errup = 0; errdown = 4.10622e-06;}
    else if (pt > 40 && pt < 50) {eff = 1; errup = 0; errdown = 2.94918e-06;}
    else if (pt > 50 && pt < 60) {eff = 1; errup = 0; errdown = 1.27312e-05;}
    else if (pt > 60 && pt < 120) {eff = 1; errup = 0; errdown = 2.58244e-05;}
    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.000449367;}
  } else if (abseta > 2.1 && abseta < 2.4) {
    if (pt > 5 && pt < 17) {eff = 1; errup = 0; errdown = 0.123222;}
    else if (pt > 17 && pt < 20) {eff = 1; errup = 0; errdown = 0.0133482;}
    else if (pt > 20 && pt < 23) {eff = 1; errup = 0; errdown = 0.00146471;}
    else if (pt > 23 && pt < 25) {eff = 1; errup = 0; errdown = 0.000549408;}
    else if (pt > 25 && pt < 30) {eff = 1; errup = 0; errdown = 8.03732e-05;}
    else if (pt > 30 && pt < 40) {eff = 1; errup = 0; errdown = 1.89466e-05;}
    else if (pt > 40 && pt < 50) {eff = 1; errup = 0; errdown = 1.69806e-05;}
    else if (pt > 50 && pt < 60) {eff = 1; errup = 0; errdown = 7.24957e-05;}
    else if (pt > 60 && pt < 120) {eff = 1; errup = 0; errdown = 0.000158929;}
    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.00343525;}
  }
  return std::vector<float>({eff, errup, errdown});
}

std::vector<float> eff_dimuleg8(float pt, float abseta) {
  float errup=0.0, errdown=0.0, eff=1.0;
  if (abseta > 0 && abseta < 0.9) {
    if (pt > 5 && pt < 7) {eff = 0.00394477; errup = 0.00517911; errdown = 0.00254755;}
    else if (pt > 7 && pt < 10) {eff = 0.694476; errup = 0.00910381; errdown = 0.00925509;}
    else if (pt > 10 && pt < 15) {eff = 0.886696; errup = 0.00247647; errdown = 0.00252401;}
    else if (pt > 15 && pt < 20) {eff = 0.906604; errup = 0.00138046; errdown = 0.00139886;}
    else if (pt > 20 && pt < 23) {eff = 0.912717; errup = 0.00125707; errdown = 0.00127352;}
    else if (pt > 23 && pt < 25) {eff = 0.911736; errup = 0.00128089; errdown = 0.00129776;}
    else if (pt > 25 && pt < 30) {eff = 0.911967; errup = 0.000642548; errdown = 0.000646796;}
    else if (pt > 30 && pt < 40) {eff = 0.911438; errup = 0.000349094; errdown = 0.000350338;}
    else if (pt > 40 && pt < 50) {eff = 0.910278; errup = 0.000442899; errdown = 0.000444873;}
    else if (pt > 50 && pt < 60) {eff = 0.912147; errup = 0.00253293; errdown = 0.00259961;}
    else if (pt > 60 && pt < 120) {eff = 0.900545; errup = 0.00334203; errdown = 0.00344308;}
    else if (pt > 120 && pt < 9999) {eff = 0.858079; errup = 0.0165571; errdown = 0.0182232;}
  } else if (abseta > 0.9 && abseta < 1.2) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 0.00355464; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 0.690919; errup = 0.0109721; errdown = 0.0111862;}
    else if (pt > 10 && pt < 15) {eff = 0.900716; errup = 0.00355662; errdown = 0.00367138;}
    else if (pt > 15 && pt < 20) {eff = 0.911025; errup = 0.00235294; errdown = 0.00240963;}
    else if (pt > 20 && pt < 23) {eff = 0.920621; errup = 0.0022307; errdown = 0.00228859;}
    else if (pt > 23 && pt < 25) {eff = 0.922566; errup = 0.00230053; errdown = 0.00236386;}
    else if (pt > 25 && pt < 30) {eff = 0.923307; errup = 0.00118194; errdown = 0.00119874;}
    else if (pt > 30 && pt < 40) {eff = 0.926376; errup = 0.000604257; errdown = 0.000608836;}
    else if (pt > 40 && pt < 50) {eff = 0.928886; errup = 0.000712551; errdown = 0.00071917;}
    else if (pt > 50 && pt < 60) {eff = 0.928262; errup = 0.004262; errdown = 0.00450123;}
    else if (pt > 60 && pt < 120) {eff = 0.921757; errup = 0.00551998; errdown = 0.00588627;}
    else if (pt > 120 && pt < 9999) {eff = 0.894737; errup = 0.026977; errdown = 0.0337182;}
  } else if (abseta > 1.2 && abseta < 2.1) {
    if (pt > 5 && pt < 7) {eff = 0.00514966; errup = 0.00163045; errdown = 0.00127199;}
    else if (pt > 7 && pt < 10) {eff = 0.625476; errup = 0.00580069; errdown = 0.00583653;}
    else if (pt > 10 && pt < 15) {eff = 0.889697; errup = 0.00210416; errdown = 0.00213953;}
    else if (pt > 15 && pt < 20) {eff = 0.913046; errup = 0.00138107; errdown = 0.00140103;}
    else if (pt > 20 && pt < 23) {eff = 0.91745; errup = 0.00137465; errdown = 0.00139561;}
    else if (pt > 23 && pt < 25) {eff = 0.916608; errup = 0.00144643; errdown = 0.00146938;}
    else if (pt > 25 && pt < 30) {eff = 0.918414; errup = 0.000758885; errdown = 0.000765339;}
    else if (pt > 30 && pt < 40) {eff = 0.922473; errup = 0.000417114; errdown = 0.000419174;}
    else if (pt > 40 && pt < 50) {eff = 0.925005; errup = 0.000473895; errdown = 0.000476654;}
    else if (pt > 50 && pt < 60) {eff = 0.926583; errup = 0.00286453; errdown = 0.00296902;}
    else if (pt > 60 && pt < 120) {eff = 0.923575; errup = 0.00379112; errdown = 0.00396696;}
    else if (pt > 120 && pt < 9999) {eff = 0.914815; errup = 0.0171477; errdown = 0.0205477;}
  } else if (abseta > 2.1 && abseta < 2.4) {
    if (pt > 5 && pt < 7) {eff = 0.00437828; errup = 0.00295102; errdown = 0.00188988;}
    else if (pt > 7 && pt < 10) {eff = 0.57137; errup = 0.0101514; errdown = 0.0102106;}
    else if (pt > 10 && pt < 15) {eff = 0.867599; errup = 0.00406624; errdown = 0.00417297;}
    else if (pt > 15 && pt < 20) {eff = 0.883745; errup = 0.00288995; errdown = 0.00295282;}
    else if (pt > 20 && pt < 23) {eff = 0.893085; errup = 0.00289799; errdown = 0.00296776;}
    else if (pt > 23 && pt < 25) {eff = 0.897214; errup = 0.00300117; errdown = 0.00307952;}
    else if (pt > 25 && pt < 30) {eff = 0.895741; errup = 0.00160517; errdown = 0.00162712;}
    else if (pt > 30 && pt < 40) {eff = 0.900261; errup = 0.000953104; errdown = 0.000961229;}
    else if (pt > 40 && pt < 50) {eff = 0.901632; errup = 0.00125176; errdown = 0.00126601;}
    else if (pt > 50 && pt < 60) {eff = 0.906128; errup = 0.00749666; errdown = 0.00804915;}
    else if (pt > 60 && pt < 120) {eff = 0.905685; errup = 0.0105964; errdown = 0.0117083;}
    else if (pt > 120 && pt < 9999) {eff = 0.9; errup = 0.0643201; errdown = 0.116971;}
  }
  return std::vector<float>({eff, errup, errdown});
}

std::vector<float> eff_isoel3235(float pt, float abseta) {
  float errup=0.0, errdown=0.0, eff=1.0;
  if (abseta > 0 && abseta < 0.8) {
    if (pt > 7 && pt < 12) {eff = 0.000882613; errup = 0.00202664; errdown = 0.00073015;}
    else if (pt > 12 && pt < 15) {eff = 0.00184162; errup = 0.004222; errdown = 0.00152352;}
    else if (pt > 15 && pt < 20) {eff = 0.00156986; errup = 0.00360059; errdown = 0.0012987;}
    else if (pt > 20 && pt < 25) {eff = 0; errup = 0.00357535; errdown = 0;}
    else if (pt > 25 && pt < 27) {eff = 0; errup = 0.00984914; errdown = 0;}
    else if (pt > 27 && pt < 30) {eff = 0; errup = 0.00713791; errdown = 0;}
    else if (pt > 30 && pt < 35) {eff = 0.312676; errup = 0.0265349; errdown = 0.0254211;}
    else if (pt > 35 && pt < 40) {eff = 0.744615; errup = 0.0248944; errdown = 0.0264924;}
    else if (pt > 40 && pt < 50) {eff = 0.753176; errup = 0.0187794; errdown = 0.0197427;}
    else if (pt > 50 && pt < 80) {eff = 0.79952; errup = 0.0140911; errdown = 0.0148408;}
    else if (pt > 80 && pt < 120) {eff = 0.821168; errup = 0.0166547; errdown = 0.0178908;}
    else if (pt > 120 && pt < 9999) {eff = 0.881057; errup = 0.0153914; errdown = 0.0171894;}
  } else if (abseta > 0.8 && abseta < 1.4442) {
    if (pt > 7 && pt < 12) {eff = 0.00135318; errup = 0.00310473; errdown = 0.00111944;}
    else if (pt > 12 && pt < 15) {eff = 0; errup = 0.00549688; errdown = 0;}
    else if (pt > 15 && pt < 20) {eff = 0; errup = 0.00418489; errdown = 0;}
    else if (pt > 20 && pt < 25) {eff = 0; errup = 0.00459198; errdown = 0;}
    else if (pt > 25 && pt < 27) {eff = 0; errup = 0.0170586; errdown = 0;}
    else if (pt > 27 && pt < 30) {eff = 0; errup = 0.0107085; errdown = 0;}
    else if (pt > 30 && pt < 35) {eff = 0.197026; errup = 0.027285; errdown = 0.0248654;}
    else if (pt > 35 && pt < 40) {eff = 0.628205; errup = 0.033073; errdown = 0.0342397;}
    else if (pt > 40 && pt < 50) {eff = 0.750716; errup = 0.0237977; errdown = 0.0253206;}
    else if (pt > 50 && pt < 80) {eff = 0.781991; errup = 0.0167243; errdown = 0.0176572;}
    else if (pt > 80 && pt < 120) {eff = 0.80826; errup = 0.0218627; errdown = 0.0238032;}
    else if (pt > 120 && pt < 9999) {eff = 0.863333; errup = 0.0201738; errdown = 0.0227942;}
  } else if (abseta > 1.4442 && abseta < 1.566) {
    if (pt > 7 && pt < 12) {eff = 0; errup = 0.0149771; errdown = 0;}
    else if (pt > 12 && pt < 15) {eff = 0; errup = 0.0329191; errdown = 0;}
    else if (pt > 15 && pt < 20) {eff = 0; errup = 0.0230346; errdown = 0;}
    else if (pt > 20 && pt < 25) {eff = 0; errup = 0.0368748; errdown = 0;}
    else if (pt > 25 && pt < 27) {eff = 0; errup = 0.0923495; errdown = 0;}
    else if (pt > 27 && pt < 30) {eff = 0.0416667; errup = 0.0893855; errdown = 0.0344944;}
    else if (pt > 30 && pt < 35) {eff = 0.25; errup = 0.108468; errdown = 0.0872543;}
    else if (pt > 35 && pt < 40) {eff = 0.518519; errup = 0.1114; errdown = 0.112994;}
    else if (pt > 40 && pt < 50) {eff = 0.552632; errup = 0.0908274; errdown = 0.0939923;}
    else if (pt > 50 && pt < 80) {eff = 0.707692; errup = 0.0600147; errdown = 0.0672026;}
    else if (pt > 80 && pt < 120) {eff = 0.777778; errup = 0.0729951; errdown = 0.0911652;}
    else if (pt > 120 && pt < 9999) {eff = 0.681818; errup = 0.109079; errdown = 0.12875;}
  } else if (abseta > 1.566 && abseta < 2.55) {
    if (pt > 7 && pt < 12) {eff = 0.00110988; errup = 0.00254752; errdown = 0.000918161;}
    else if (pt > 12 && pt < 15) {eff = 0; errup = 0.00496338; errdown = 0;}
    else if (pt > 15 && pt < 20) {eff = 0; errup = 0.00378094; errdown = 0;}
    else if (pt > 20 && pt < 25) {eff = 0; errup = 0.00507278; errdown = 0;}
    else if (pt > 25 && pt < 27) {eff = 0; errup = 0.0133482; errdown = 0;}
    else if (pt > 27 && pt < 30) {eff = 0.00537634; errup = 0.0122537; errdown = 0.00444799;}
    else if (pt > 30 && pt < 35) {eff = 0.0983607; errup = 0.0229311; errdown = 0.0192944;}
    else if (pt > 35 && pt < 40) {eff = 0.5; errup = 0.0371673; errdown = 0.0371673;}
    else if (pt > 40 && pt < 50) {eff = 0.656627; errup = 0.0270269; errdown = 0.028023;}
    else if (pt > 50 && pt < 80) {eff = 0.737705; errup = 0.0192155; errdown = 0.0201225;}
    else if (pt > 80 && pt < 120) {eff = 0.771127; errup = 0.0256314; errdown = 0.0276711;}
    else if (pt > 120 && pt < 9999) {eff = 0.795276; errup = 0.0259964; errdown = 0.028496;}
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
    else if (pt > 20 && pt < 23) {eff = 0.0021645; errup = 0.00495958; errdown = 0.00179065;}
    else if (pt > 23 && pt < 25) {eff = 0.292683; errup = 0.0318821; errdown = 0.0300831;}
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
    else if (pt > 20 && pt < 23) {eff = 0.00714286; errup = 0.0162324; errdown = 0.00590966;}
    else if (pt > 23 && pt < 25) {eff = 0.228571; errup = 0.0613516; errdown = 0.0525767;}
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
    else if (pt > 20 && pt < 23) {eff = 0.00310559; errup = 0.00710486; errdown = 0.00256923;}
    else if (pt > 23 && pt < 25) {eff = 0.237885; errup = 0.0316479; errdown = 0.0291654;}
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
    else if (pt > 20 && pt < 23) {eff = 0; errup = 0.0219368; errdown = 0;}
    else if (pt > 23 && pt < 25) {eff = 0.148936; errup = 0.0709778; errdown = 0.0532465;}
    else if (pt > 25 && pt < 30) {eff = 0.447368; errup = 0.0512155; errdown = 0.0502109;}
    else if (pt > 30 && pt < 40) {eff = 0.673077; errup = 0.0394325; errdown = 0.041831;}
    else if (pt > 40 && pt < 50) {eff = 0.734694; errup = 0.0468096; errdown = 0.0521118;}
    else if (pt > 50 && pt < 60) {eff = 0.696203; errup = 0.0549362; errdown = 0.0604646;}
    else if (pt > 60 && pt < 120) {eff = 0.725581; errup = 0.0315616; errdown = 0.0338138;}
    else if (pt > 120 && pt < 9999) {eff = 0.625; errup = 0.0842352; errdown = 0.0913836;}
  }
  return std::vector<float>({eff, errup, errdown});
}

//std::vector<float> effsig_elhltpt23(float pt, float abseta) {
//  float errup=0.0, errdown=0.0, eff=1.0;
//  if (abseta > 0 && abseta < 0.8) {
//    if (pt > 7 && pt < 15) {eff = 0.00724638; errup = 0.00167264; errdown = 0.00138259;}
//    else if (pt > 15 && pt < 20) {eff = 0.0025385; errup = 0.000837766; errdown = 0.000647606;}
//    else if (pt > 20 && pt < 25) {eff = 0.3921; errup = 0.00494285; errdown = 0.00492105;}
//    else if (pt > 25 && pt < 27) {eff = 0.984879; errup = 0.00164644; errdown = 0.00183288;}
//    else if (pt > 27 && pt < 30) {eff = 0.991322; errup = 0.000941662; errdown = 0.00104874;}
//    else if (pt > 30 && pt < 35) {eff = 0.997816; errup = 0.000327728; errdown = 0.000380883;}
//    else if (pt > 35 && pt < 40) {eff = 0.998246; errup = 0.000286852; errdown = 0.000337983;}
//    else if (pt > 40 && pt < 50) {eff = 0.998956; errup = 0.000166363; errdown = 0.000195216;}
//    else if (pt > 50 && pt < 80) {eff = 0.999566; errup = 9.62023e-05; errdown = 0.00012031;}
//    else if (pt > 80 && pt < 120) {eff = 0.999583; errup = 0.000226693; errdown = 0.000405025;}
//    else if (pt > 120 && pt < 9999) {eff = 0.999545; errup = 0.000376195; errdown = 0.00104493;}
//  } else if (abseta > 0.8 && abseta < 1.4442) {
//    if (pt > 7 && pt < 15) {eff = 0.0122699; errup = 0.00238707; errdown = 0.00202674;}
//    else if (pt > 15 && pt < 20) {eff = 0.00745763; errup = 0.00153038; errdown = 0.0012883;}
//    else if (pt > 20 && pt < 25) {eff = 0.347758; errup = 0.0057846; errdown = 0.00574057;}
//    else if (pt > 25 && pt < 27) {eff = 0.97301; errup = 0.0026918; errdown = 0.00296566;}
//    else if (pt > 27 && pt < 30) {eff = 0.991879; errup = 0.00110865; errdown = 0.0012698;}
//    else if (pt > 30 && pt < 35) {eff = 0.996342; errup = 0.000530981; errdown = 0.0006139;}
//    else if (pt > 35 && pt < 40) {eff = 0.998287; errup = 0.000354346; errdown = 0.000436322;}
//    else if (pt > 40 && pt < 50) {eff = 0.998902; errup = 0.000213879; errdown = 0.000260163;}
//    else if (pt > 50 && pt < 80) {eff = 0.99941; errup = 0.000141605; errdown = 0.000180471;}
//    else if (pt > 80 && pt < 120) {eff = 0.999772; errup = 0.000188741; errdown = 0.000524449;}
//    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.001413;}
//  } else if (abseta > 1.4442 && abseta < 1.566) {
//    if (pt > 7 && pt < 15) {eff = 0.0268097; errup = 0.0112204; errdown = 0.00828156;}
//    else if (pt > 15 && pt < 20) {eff = 0.0780669; errup = 0.013347; errdown = 0.0116542;}
//    else if (pt > 20 && pt < 25) {eff = 0.51677; errup = 0.0182006; errdown = 0.0182437;}
//    else if (pt > 25 && pt < 27) {eff = 0.954545; errup = 0.0107703; errdown = 0.0134984;}
//    else if (pt > 27 && pt < 30) {eff = 0.984263; errup = 0.00465375; errdown = 0.00624858;}
//    else if (pt > 30 && pt < 35) {eff = 0.997638; errup = 0.00128533; errdown = 0.00229243;}
//    else if (pt > 35 && pt < 40) {eff = 0.999214; errup = 0.00065036; errdown = 0.00180546;}
//    else if (pt > 40 && pt < 50) {eff = 1.0; errup = 0; errdown = 0.000863146;}
//    else if (pt > 50 && pt < 80) {eff = 0.999204; errup = 0.000514446; errdown = 0.00104955;}
//    else if (pt > 80 && pt < 120) {eff = 1.0; errup = 0; errdown = 0.00501748;}
//    else if (pt > 120 && pt < 9999) {eff = 1.0; errup = 0; errdown = 0.0216785;}
//  } else if (abseta > 1.566 && abseta < 2.55) {
//    if (pt > 7 && pt < 15) {eff = 0.0122093; errup = 0.00217379; errdown = 0.00186818;}
//    else if (pt > 15 && pt < 20) {eff = 0.00805659; errup = 0.00145798; errdown = 0.00124947;}
//    else if (pt > 20 && pt < 25) {eff = 0.429984; errup = 0.00583826; errdown = 0.00581907;}
//    else if (pt > 25 && pt < 27) {eff = 0.982967; errup = 0.00214326; errdown = 0.00242561;}
//    else if (pt > 27 && pt < 30) {eff = 0.994318; errup = 0.000953913; errdown = 0.00112803;}
//    else if (pt > 30 && pt < 35) {eff = 0.997884; errup = 0.00042869; errdown = 0.000525489;}
//    else if (pt > 35 && pt < 40) {eff = 0.999305; errup = 0.000240402; errdown = 0.000342461;}
//    else if (pt > 40 && pt < 50) {eff = 0.99979; errup = 0.000100636; errdown = 0.000166237;}
//    else if (pt > 50 && pt < 80) {eff = 0.999822; errup = 8.50874e-05; errdown = 0.000140557;}
//    else if (pt > 80 && pt < 120) {eff = 1; errup = 0; errdown = 0.000623036;}
//    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.00268793;}
//  }
//  return std::vector<float>({eff, errup, errdown});
//}
//
//std::vector<float> effsig_dielleg12(float pt, float abseta) {
//  float errup=0.0, errdown=0.0, eff=1.0;
//  if (abseta > 0 && abseta < 0.8) {
//    if (pt > 7 && pt < 12) {eff = 0.0630279; errup = 0.00682027; errdown = 0.00622281;}
//    else if (pt > 12 && pt < 15) {eff = 0.697987; errup = 0.0110214; errdown = 0.0112484;}
//    else if (pt > 15 && pt < 20) {eff = 0.810188; errup = 0.00538091; errdown = 0.00549824;}
//    else if (pt > 20 && pt < 25) {eff = 0.846785; errup = 0.0037829; errdown = 0.00385984;}
//    else if (pt > 25 && pt < 27) {eff = 0.859051; errup = 0.00501509; errdown = 0.00516566;}
//    else if (pt > 27 && pt < 30) {eff = 0.857296; errup = 0.00383211; errdown = 0.00391848;}
//    else if (pt > 30 && pt < 35) {eff = 0.882569; errup = 0.00250314; errdown = 0.0025497;}
//    else if (pt > 35 && pt < 40) {eff = 0.900163; errup = 0.0024278; errdown = 0.00248072;}
//    else if (pt > 40 && pt < 50) {eff = 0.921381; errup = 0.00226005; errdown = 0.00232013;}
//    else if (pt > 50 && pt < 80) {eff = 0.902297; errup = 0.00570022; errdown = 0.00600292;}
//    else if (pt > 80 && pt < 120) {eff = 0.8753; errup = 0.0164066; errdown = 0.0183373;}
//    else if (pt > 120 && pt < 9999) {eff = 0.787234; errup = 0.0439033; errdown = 0.0507369;}
//  } else if (abseta > 0.8 && abseta < 1.4442) {
//    if (pt > 7 && pt < 12) {eff = 0.0513001; errup = 0.00653704; errdown = 0.00586873;}
//    else if (pt > 12 && pt < 15) {eff = 0.553753; errup = 0.0132221; errdown = 0.0132967;}
//    else if (pt > 15 && pt < 20) {eff = 0.762376; errup = 0.00693025; errdown = 0.00706954;}
//    else if (pt > 20 && pt < 25) {eff = 0.817531; errup = 0.00495402; errdown = 0.005059;}
//    else if (pt > 25 && pt < 27) {eff = 0.832073; errup = 0.00668484; errdown = 0.00689902;}
//    else if (pt > 27 && pt < 30) {eff = 0.843811; errup = 0.00485738; errdown = 0.00498137;}
//    else if (pt > 30 && pt < 35) {eff = 0.855324; errup = 0.00341711; errdown = 0.00348457;}
//    else if (pt > 35 && pt < 40) {eff = 0.886279; errup = 0.00319141; errdown = 0.00327017;}
//    else if (pt > 40 && pt < 50) {eff = 0.903391; errup = 0.00309467; errdown = 0.00318419;}
//    else if (pt > 50 && pt < 80) {eff = 0.88322; errup = 0.00770456; errdown = 0.00815502;}
//    else if (pt > 80 && pt < 120) {eff = 0.937799; errup = 0.0167442; errdown = 0.021519;}
//    else if (pt > 120 && pt < 9999) {eff = 0.857143; errup = 0.0602801; errdown = 0.085095;}
//  } else if (abseta > 1.4442 && abseta < 1.566) {
//    if (pt > 7 && pt < 12) {eff = 0.0572519; errup = 0.018201; errdown = 0.014381;}
//    else if (pt > 12 && pt < 15) {eff = 0.408451; errup = 0.0312131; errdown = 0.030531;}
//    else if (pt > 15 && pt < 20) {eff = 0.604096; errup = 0.0208523; errdown = 0.0212222;}
//    else if (pt > 20 && pt < 25) {eff = 0.69759; errup = 0.016289; errdown = 0.0167827;}
//    else if (pt > 25 && pt < 27) {eff = 0.70235; errup = 0.0240816; errdown = 0.0251955;}
//    else if (pt > 27 && pt < 30) {eff = 0.74024; errup = 0.0173547; errdown = 0.0181073;}
//    else if (pt > 30 && pt < 35) {eff = 0.776224; errup = 0.0125052; errdown = 0.0130049;}
//    else if (pt > 35 && pt < 40) {eff = 0.811434; errup = 0.0125608; errdown = 0.0132105;}
//    else if (pt > 40 && pt < 50) {eff = 0.830857; errup = 0.0128456; errdown = 0.0136354;}
//    else if (pt > 50 && pt < 80) {eff = 0.802548; errup = 0.0327593; errdown = 0.0369825;}
//    else if (pt > 80 && pt < 120) {eff = 0.833333; errup = 0.106875; errdown = 0.17901;}
//    else if (pt > 120 && pt < 9999) {eff = 0.83333; errup = 0.106875; errdown = 0.17901;}
//  } else if (abseta > 1.566 && abseta < 2.55) {
//    if (pt > 7 && pt < 12) {eff = 0.0683761; errup = 0.00595204; errdown = 0.00552454;}
//    else if (pt > 12 && pt < 15) {eff = 0.560363; errup = 0.0117027; errdown = 0.0117686;}
//    else if (pt > 15 && pt < 20) {eff = 0.77093; errup = 0.0062948; errdown = 0.0064166;}
//    else if (pt > 20 && pt < 25) {eff = 0.825067; errup = 0.00479366; errdown = 0.0048977;}
//    else if (pt > 25 && pt < 27) {eff = 0.856566; errup = 0.00647652; errdown = 0.00672287;}
//    else if (pt > 27 && pt < 30) {eff = 0.875331; errup = 0.00473853; errdown = 0.00489473;}
//    else if (pt > 30 && pt < 35) {eff = 0.896078; errup = 0.00333336; errdown = 0.0034289;}
//    else if (pt > 35 && pt < 40) {eff = 0.924764; errup = 0.0030313; errdown = 0.00314526;}
//    else if (pt > 40 && pt < 50) {eff = 0.933022; errup = 0.00297978; errdown = 0.00310506;}
//    else if (pt > 50 && pt < 80) {eff = 0.908954; errup = 0.00794582; errdown = 0.00859033;}
//    else if (pt > 80 && pt < 120) {eff = 0.891473; errup = 0.0277863; errdown = 0.0346796;}
//    else if (pt > 120 && pt < 9999) {eff = 0.888889; errup = 0.0714316; errdown = 0.128174;}
//  }
//  return std::vector<float>({eff, errup, errdown});
//}
//
//std::vector<float> effsig_muhltpt17(float pt, float abseta) {
//  float errup=0.0, errdown=0.0, eff=1.0;
//  if (abseta > 0 && abseta < 0.9) {
//    if (pt > 5 && pt < 17) {eff = 0.00763889; errup = 0.00132795; errdown = 0.00114404;}
//    else if (pt > 17 && pt < 20) {eff = 0.983089; errup = 0.0020323; errdown = 0.00228727;}
//    else if (pt > 20 && pt < 23) {eff = 0.999648; errup = 0.000227269; errdown = 0.000463896;}
//    else if (pt > 23 && pt < 25) {eff = 0.999589; errup = 0.000265146; errdown = 0.000541175;}
//    else if (pt > 25 && pt < 30) {eff = 1; errup = 0; errdown = 0.000111395;}
//    else if (pt > 30 && pt < 40) {eff = 0.999958; errup = 2.73846e-05; errdown = 5.59161e-05;}
//    else if (pt > 40 && pt < 50) {eff = 1; errup = 0; errdown = 4.22157e-05;}
//    else if (pt > 50 && pt < 60) {eff = 1; errup = 0; errdown = 6.07921e-05;}
//    else if (pt > 60 && pt < 120) {eff = 1; errup = 0; errdown = 6.23886e-05;}
//    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.000731458;}
//  } else if (abseta > 0.9 && abseta < 1.2) {
//    if (pt > 5 && pt < 17) {eff = 0.0212404; errup = 0.00340596; errdown = 0.00297018;}
//    else if (pt > 17 && pt < 20) {eff = 0.970082; errup = 0.00466012; errdown = 0.00541952;}
//    else if (pt > 20 && pt < 23) {eff = 0.991648; errup = 0.00212726; errdown = 0.00274564;}
//    else if (pt > 23 && pt < 25) {eff = 0.995935; errup = 0.00161127; errdown = 0.00242006;}
//    else if (pt > 25 && pt < 30) {eff = 0.997731; errup = 0.000673359; errdown = 0.000909737;}
//    else if (pt > 30 && pt < 40) {eff = 0.998883; errup = 0.000285027; errdown = 0.000368927;}
//    else if (pt > 40 && pt < 50) {eff = 0.999026; errup = 0.000277145; errdown = 0.000369866;}
//    else if (pt > 50 && pt < 60) {eff = 0.999546; errup = 0.000217033; errdown = 0.000358443;}
//    else if (pt > 60 && pt < 120) {eff = 0.999879; errup = 0.000100394; errdown = 0.000279013;}
//    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.00269581;}
//  } else if (abseta > 1.2 && abseta < 2.1) {
//    if (pt > 5 && pt < 17) {eff = 0.0159456; errup = 0.00156442; errdown = 0.00143212;}
//    else if (pt > 17 && pt < 20) {eff = 0.956369; errup = 0.00337066; errdown = 0.00362795;}
//    else if (pt > 20 && pt < 23) {eff = 0.997066; errup = 0.000802513; errdown = 0.00105804;}
//    else if (pt > 23 && pt < 25) {eff = 0.998004; errup = 0.000735807; errdown = 0.00107342;}
//    else if (pt > 25 && pt < 30) {eff = 0.998601; errup = 0.000357047; errdown = 0.000462095;}
//    else if (pt > 30 && pt < 40) {eff = 0.999424; errup = 0.000142549; errdown = 0.000183038;}
//    else if (pt > 40 && pt < 50) {eff = 0.999518; errup = 0.000137049; errdown = 0.00018294;}
//    else if (pt > 50 && pt < 60) {eff = 0.999824; errup = 9.59215e-05; errdown = 0.000171416;}
//    else if (pt > 60 && pt < 120) {eff = 0.999809; errup = 0.000103695; errdown = 0.000185306;}
//    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.00168604;}
//  } else if (abseta > 2.1 && abseta < 2.4) {
//    if (pt > 5 && pt < 17) {eff = 0.0236686; errup = 0.00355498; errdown = 0.00312516;}
//    else if (pt > 17 && pt < 20) {eff = 0.920956; errup = 0.00823274; errdown = 0.00904963;}
//    else if (pt > 20 && pt < 23) {eff = 0.997517; errup = 0.00135128; errdown = 0.0024098;}
//    else if (pt > 23 && pt < 25) {eff = 0.996466; errup = 0.00192249; errdown = 0.00342523;}
//    else if (pt > 25 && pt < 30) {eff = 0.997329; errup = 0.000984409; errdown = 0.00143553;}
//    else if (pt > 30 && pt < 40) {eff = 0.998447; errup = 0.000507668; errdown = 0.000708343;}
//    else if (pt > 40 && pt < 50) {eff = 0.999123; errup = 0.000419489; errdown = 0.000692586;}
//    else if (pt > 50 && pt < 60) {eff = 0.998981; errup = 0.000554726; errdown = 0.000990572;}
//    else if (pt > 60 && pt < 120) {eff = 0.999576; errup = 0.000350977; errdown = 0.000974933;}
//    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.0165973;}
//  }
//  return std::vector<float>({eff, errup, errdown});
//}
//
//std::vector<float> effsig_dimuleg8(float pt, float abseta) {
//  float errup=0.0, errdown=0.0, eff=1.0;
//  if (abseta > 0 && abseta < 0.9) {
//    if (pt > 5 && pt < 7) {eff = 0.032; errup = 0.015403; errdown = 0.0109986;}
//    else if (pt > 7 && pt < 10) {eff = 0.64775; errup = 0.021786; errdown = 0.0223903;}
//    else if (pt > 10 && pt < 15) {eff = 0.928474; errup = 0.00552597; errdown = 0.00593245;}
//    else if (pt > 15 && pt < 20) {eff = 0.94332; errup = 0.0032194; errdown = 0.00339546;}
//    else if (pt > 20 && pt < 23) {eff = 0.94717; errup = 0.0031607; errdown = 0.00334391;}
//    else if (pt > 23 && pt < 25) {eff = 0.944058; errup = 0.00351148; errdown = 0.00372444;}
//    else if (pt > 25 && pt < 30) {eff = 0.949019; errup = 0.00184759; errdown = 0.00191194;}
//    else if (pt > 30 && pt < 40) {eff = 0.951852; errup = 0.00112251; errdown = 0.00114761;}
//    else if (pt > 40 && pt < 50) {eff = 0.948571; errup = 0.0016979; errdown = 0.00175167;}
//    else if (pt > 50 && pt < 60) {eff = 0.935252; errup = 0.0054123; errdown = 0.00584829;}
//    else if (pt > 60 && pt < 120) {eff = 0.94639; errup = 0.0060448; errdown = 0.00672059;}
//    else if (pt > 120 && pt < 9999) {eff = 0.915888; errup = 0.0269782; errdown = 0.0360474;}
//  } else if (abseta > 0.9 && abseta < 1.2) {
//    if (pt > 5 && pt < 7) {eff = 0.0215054; errup = 0.0276606; errdown = 0.0138784;}
//    else if (pt > 7 && pt < 10) {eff = 0.654275; errup = 0.0301908; errdown = 0.0314085;}
//    else if (pt > 10 && pt < 15) {eff = 0.913846; errup = 0.009052; errdown = 0.00994886;}
//    else if (pt > 15 && pt < 20) {eff = 0.930911; errup = 0.00621956; errdown = 0.00675739;}
//    else if (pt > 20 && pt < 23) {eff = 0.937035; errup = 0.00634943; errdown = 0.0069727;}
//    else if (pt > 23 && pt < 25) {eff = 0.954254; errup = 0.00593592; errdown = 0.00671419;}
//    else if (pt > 25 && pt < 30) {eff = 0.956844; errup = 0.00325468; errdown = 0.00349715;}
//    else if (pt > 30 && pt < 40) {eff = 0.963099; errup = 0.00187782; errdown = 0.00197186;}
//    else if (pt > 40 && pt < 50) {eff = 0.964016; errup = 0.00266773; errdown = 0.00286475;}
//    else if (pt > 50 && pt < 60) {eff = 0.954939; errup = 0.00865013; errdown = 0.0103837;}
//    else if (pt > 60 && pt < 120) {eff = 0.962667; errup = 0.00975722; errdown = 0.0125564;}
//    else if (pt > 120 && pt < 9999) {eff = 0.958333; errup = 0.0344944; errdown = 0.0893855;}
//  } else if (abseta > 1.2 && abseta < 2.1) {
//    if (pt > 5 && pt < 7) {eff = 0.025; errup = 0.00986135; errdown = 0.0073761;}
//    else if (pt > 7 && pt < 10) {eff = 0.626601; errup = 0.015539; errdown = 0.0157964;}
//    else if (pt > 10 && pt < 15) {eff = 0.921727; errup = 0.00495568; errdown = 0.00524992;}
//    else if (pt > 15 && pt < 20) {eff = 0.929107; errup = 0.00376887; errdown = 0.00395791;}
//    else if (pt > 20 && pt < 23) {eff = 0.934867; errup = 0.00412268; errdown = 0.00437184;}
//    else if (pt > 23 && pt < 25) {eff = 0.942826; errup = 0.00441778; errdown = 0.0047492;}
//    else if (pt > 25 && pt < 30) {eff = 0.938366; errup = 0.00261198; errdown = 0.00271722;}
//    else if (pt > 30 && pt < 40) {eff = 0.940006; errup = 0.0016673; errdown = 0.00171118;}
//    else if (pt > 40 && pt < 50) {eff = 0.943494; errup = 0.00235609; errdown = 0.00245008;}
//    else if (pt > 50 && pt < 60) {eff = 0.937558; errup = 0.00742272; errdown = 0.00828933;}
//    else if (pt > 60 && pt < 120) {eff = 0.96248; errup = 0.00767382; errdown = 0.00934163;}
//    else if (pt > 120 && pt < 9999) {eff = 0.974359; errup = 0.0212212; errdown = 0.0565052;}
//  } else if (abseta > 2.1 && abseta < 2.4) {
//    if (pt > 5 && pt < 7) {eff = 0.00571429; errup = 0.0130167; errdown = 0.00472761;}
//    else if (pt > 7 && pt < 10) {eff = 0.616246; errup = 0.0267568; errdown = 0.0274423;}
//    else if (pt > 10 && pt < 15) {eff = 0.918592; errup = 0.00913419; errdown = 0.0101105;}
//    else if (pt > 15 && pt < 20) {eff = 0.929094; errup = 0.00697686; errdown = 0.00763698;}
//    else if (pt > 20 && pt < 23) {eff = 0.939919; errup = 0.00761905; errdown = 0.00857467;}
//    else if (pt > 23 && pt < 25) {eff = 0.933234; errup = 0.0096708; errdown = 0.0110544;}
//    else if (pt > 25 && pt < 30) {eff = 0.939333; errup = 0.00538936; errdown = 0.00585438;}
//    else if (pt > 30 && pt < 40) {eff = 0.934152; errup = 0.00402317; errdown = 0.00425742;}
//    else if (pt > 40 && pt < 50) {eff = 0.94913; errup = 0.00570444; errdown = 0.00634088;}
//    else if (pt > 50 && pt < 60) {eff = 0.946429; errup = 0.0173078; errdown = 0.0235078;}
//    else if (pt > 60 && pt < 120) {eff = 0.935897; errup = 0.0274106; errdown = 0.041038;}
//    else if (pt > 120 && pt < 9999) {eff = 0.935897; errup = 0.0274106; errdown = 0.041038;}
//  }
//  return std::vector<float>({eff, errup, errdown});
//}
//
//std::vector<float> effsig_isoel3235(float pt, float abseta) {
//  float errup=0.0, errdown=0.0, eff=1.0;
//  if (abseta > 0 && abseta < 0.8) {
//    if (pt > 7 && pt < 12) {eff = 0.0125294; errup = 0.00394757; errdown = 0.00308833;}
//    else if (pt > 12 && pt < 15) {eff = 0.016129; errup = 0.00442683; errdown = 0.00355869;}
//    else if (pt > 15 && pt < 20) {eff = 0.0133279; errup = 0.00272527; errdown = 0.00229774;}
//    else if (pt > 20 && pt < 25) {eff = 0.00995223; errup = 0.00240009; errdown = 0.00197073;}
//    else if (pt > 25 && pt < 27) {eff = 0.0143488; errup = 0.00513425; errdown = 0.00391309;}
//    else if (pt > 27 && pt < 30) {eff = 0.00997698; errup = 0.00358069; errdown = 0.00272399;}
//    else if (pt > 30 && pt < 35) {eff = 0.37397; errup = 0.011681; errdown = 0.0115392;}
//    else if (pt > 35 && pt < 40) {eff = 0.798653; errup = 0.0111114; errdown = 0.0115735;}
//    else if (pt > 40 && pt < 50) {eff = 0.818634; errup = 0.0103951; errdown = 0.0108639;}
//    else if (pt > 50 && pt < 80) {eff = 0.845096; errup = 0.0110805; errdown = 0.0117392;}
//    else if (pt > 80 && pt < 120) {eff = 0.871658; errup = 0.024908; errdown = 0.0293026;}
//    else if (pt > 120 && pt < 9999) {eff = 0.773585; errup = 0.0604081; errdown = 0.0722734;}
//  } else if (abseta > 0.8 && abseta < 1.4442) {
//    if (pt > 7 && pt < 12) {eff = 0.0108696; errup = 0.00728622; errdown = 0.00468716;}
//    else if (pt > 12 && pt < 15) {eff = 0.00817996; errup = 0.00642045; errdown = 0.0039108;}
//    else if (pt > 15 && pt < 20) {eff = 0.0150435; errup = 0.00426169; errdown = 0.00340524;}
//    else if (pt > 20 && pt < 25) {eff = 0.0149834; errup = 0.00344192; errdown = 0.00285151;}
//    else if (pt > 25 && pt < 27) {eff = 0.0280612; errup = 0.00721808; errdown = 0.00588481;}
//    else if (pt > 27 && pt < 30) {eff = 0.0257689; errup = 0.00542272; errdown = 0.00456287;}
//    else if (pt > 30 && pt < 35) {eff = 0.382271; errup = 0.0117785; errdown = 0.0116449;}
//    else if (pt > 35 && pt < 40) {eff = 0.760475; errup = 0.0108156; errdown = 0.0111508;}
//    else if (pt > 40 && pt < 50) {eff = 0.795504; errup = 0.0084703; errdown = 0.00873234;}
//    else if (pt > 50 && pt < 80) {eff = 0.81698; errup = 0.00880576; errdown = 0.00913736;}
//    else if (pt > 80 && pt < 120) {eff = 0.831884; errup = 0.0205354; errdown = 0.0225942;}
//    else if (pt > 120 && pt < 9999) {eff = 0.808219; errup = 0.0478339; errdown = 0.057445;}
//  } else if (abseta > 1.4442 && abseta < 1.566) {
//    if (pt > 7 && pt < 12) {eff = 0; errup = 0.0323409; errdown = 0;}
//    else if (pt > 12 && pt < 15) {eff = 0.0172414; errup = 0.0385306; errdown = 0.0142673;}
//    else if (pt > 15 && pt < 20) {eff = 0; errup = 0.0132521; errdown = 0;}
//    else if (pt > 20 && pt < 25) {eff = 0.030303; errup = 0.0176611; errdown = 0.0119554;}
//    else if (pt > 25 && pt < 27) {eff = 0.045045; errup = 0.0293233; errdown = 0.0193206;}
//    else if (pt > 27 && pt < 30) {eff = 0.0927152; errup = 0.0299748; errdown = 0.0238499;}
//    else if (pt > 30 && pt < 35) {eff = 0.41637; errup = 0.0314438; errdown = 0.030814;}
//    else if (pt > 35 && pt < 40) {eff = 0.648889; errup = 0.0332612; errdown = 0.0346731;}
//    else if (pt > 40 && pt < 50) {eff = 0.670455; errup = 0.0259267; errdown = 0.0269483;}
//    else if (pt > 50 && pt < 80) {eff = 0.722543; errup = 0.024802; errdown = 0.0261632;}
//    else if (pt > 80 && pt < 120) {eff = 0.6; errup = 0.0635978; errdown = 0.0667734;}
//    else if (pt > 120 && pt < 9999) {eff = 0.894737; errup = 0.0676897; errdown = 0.122322;}
//  } else if (abseta > 1.566 && abseta < 2.55) {
//    if (pt > 7 && pt < 12) {eff = 0.012987; errup = 0.0168719; errdown = 0.00838402;}
//    else if (pt > 12 && pt < 15) {eff = 0.017341; errup = 0.0165803; errdown = 0.00942304;}
//    else if (pt > 15 && pt < 20) {eff = 0.0190641; errup = 0.00755187; errdown = 0.00563308;}
//    else if (pt > 20 && pt < 25) {eff = 0.0278716; errup = 0.00564841; errdown = 0.0047808;}
//    else if (pt > 25 && pt < 27) {eff = 0.0239617; errup = 0.0077947; errdown = 0.00607655;}
//    else if (pt > 27 && pt < 30) {eff = 0.033123; errup = 0.00582353; errdown = 0.00502962;}
//    else if (pt > 30 && pt < 35) {eff = 0.463689; errup = 0.0100338; errdown = 0.0100051;}
//    else if (pt > 35 && pt < 40) {eff = 0.745782; errup = 0.00852334; errdown = 0.00871171;}
//    else if (pt > 40 && pt < 50) {eff = 0.782019; errup = 0.00600705; errdown = 0.00612675;}
//    else if (pt > 50 && pt < 80) {eff = 0.831632; errup = 0.00503684; errdown = 0.00515776;}
//    else if (pt > 80 && pt < 120) {eff = 0.820896; errup = 0.0131773; errdown = 0.0139462;}
//    else if (pt > 120 && pt < 9999) {eff = 0.819383; errup = 0.0261713; errdown = 0.0292196;}
//  }
//  return std::vector<float>({eff, errup, errdown});
//}
//
//std::vector<float> effsig_isomu2427(float pt, float abseta) {
//  float errup=0.0, errdown=0.0, eff=1.0;
//  if (abseta > 0 && abseta < 0.9) {
//    if (pt > 5 && pt < 7) {eff = 0; errup = 0.00833336; errdown = 0;}
//    else if (pt > 7 && pt < 10) {eff = 0; errup = 0.00260428; errdown = 0;}
//    else if (pt > 10 && pt < 15) {eff = 0.000465116; errup = 0.00106873; errdown = 0.000384769;}
//    else if (pt > 15 && pt < 20) {eff = 0.000301841; errup = 0.000693746; errdown = 0.000249698;}
//    else if (pt > 20 && pt < 23) {eff = 0.000466418; errup = 0.00107172; errdown = 0.000385846;}
//    else if (pt > 23 && pt < 25) {eff = 0.389011; errup = 0.0136406; errdown = 0.0134735;}
//    else if (pt > 25 && pt < 30) {eff = 0.882835; errup = 0.00574128; errdown = 0.00598905;}
//    else if (pt > 30 && pt < 40) {eff = 0.905059; errup = 0.0043194; errdown = 0.00449815;}
//    else if (pt > 40 && pt < 50) {eff = 0.921097; errup = 0.00556485; errdown = 0.00593366;}
//    else if (pt > 50 && pt < 60) {eff = 0.934378; errup = 0.00779056; errdown = 0.00869475;}
//    else if (pt > 60 && pt < 120) {eff = 0.923459; errup = 0.00843587; errdown = 0.00932718;}
//    else if (pt > 120 && pt < 9999) {eff = 0.954023; errup = 0.0218756; errdown = 0.0348613;}
//  } else if (abseta > 0.9 && abseta < 1.2) {
//    if (pt > 5 && pt < 7) {eff = 0.0188679; errup = 0.0420516; errdown = 0.0156137;}
//    else if (pt > 7 && pt < 10) {eff = 0; errup = 0.0121182; errdown = 0;}
//    else if (pt > 10 && pt < 15) {eff = 0.00240964; errup = 0.00551903; errdown = 0.00199345;}
//    else if (pt > 15 && pt < 20) {eff = 0.00251256; errup = 0.00330425; errdown = 0.00162272;}
//    else if (pt > 20 && pt < 23) {eff = 0.00499168; errup = 0.00483166; errdown = 0.00271545;}
//    else if (pt > 23 && pt < 25) {eff = 0.364606; errup = 0.0235687; errdown = 0.0229645;}
//    else if (pt > 25 && pt < 30) {eff = 0.842059; errup = 0.0109887; errdown = 0.0116202;}
//    else if (pt > 30 && pt < 40) {eff = 0.903894; errup = 0.00658346; errdown = 0.00699646;}
//    else if (pt > 40 && pt < 50) {eff = 0.937324; errup = 0.00646185; errdown = 0.00711129;}
//    else if (pt > 50 && pt < 60) {eff = 0.953523; errup = 0.00817033; errdown = 0.00965411;}
//    else if (pt > 60 && pt < 120) {eff = 0.9374; errup = 0.00975599; errdown = 0.0112749;}
//    else if (pt > 120 && pt < 9999) {eff = 0.977273; errup = 0.0188087; errdown = 0.0503283;}
//  } else if (abseta > 1.2 && abseta < 2.1) {
//    if (pt > 5 && pt < 7) {eff = 0; errup = 0.0189946; errdown = 0;}
//    else if (pt > 7 && pt < 10) {eff = 0.00574713; errup = 0.0130908; errdown = 0.00475478;}
//    else if (pt > 10 && pt < 15) {eff = 0.00583658; errup = 0.00564474; errdown = 0.00317483;}
//    else if (pt > 15 && pt < 20) {eff = 0.00480307; errup = 0.00323617; errdown = 0.00207311;}
//    else if (pt > 20 && pt < 23) {eff = 0.00849057; errup = 0.00385356; errdown = 0.0027716;}
//    else if (pt > 23 && pt < 25) {eff = 0.431843; errup = 0.0169699; errdown = 0.0168164;}
//    else if (pt > 25 && pt < 30) {eff = 0.795406; errup = 0.00765207; errdown = 0.00786566;}
//    else if (pt > 30 && pt < 40) {eff = 0.829117; errup = 0.00460275; errdown = 0.0047017;}
//    else if (pt > 40 && pt < 50) {eff = 0.862206; errup = 0.00455469; errdown = 0.00468231;}
//    else if (pt > 50 && pt < 60) {eff = 0.878103; errup = 0.00549537; errdown = 0.00571166;}
//    else if (pt > 60 && pt < 120) {eff = 0.892613; errup = 0.00548426; errdown = 0.00573486;}
//    else if (pt > 120 && pt < 9999) {eff = 0.886719; errup = 0.0200973; errdown = 0.0234085;}
//  } else if (abseta > 2.1 && abseta < 2.4) {
//    if (pt > 5 && pt < 7) {eff = 0; errup = 0.15411; errdown = 0;}
//    else if (pt > 7 && pt < 10) {eff = 0; errup = 0.0923495; errdown = 0;}
//    else if (pt > 10 && pt < 15) {eff = 0; errup = 0.0485397; errdown = 0;}
//    else if (pt > 15 && pt < 20) {eff = 0.0178571; errup = 0.0230669; errdown = 0.0115258;}
//    else if (pt > 20 && pt < 23) {eff = 0; errup = 0.0148562; errdown = 0;}
//    else if (pt > 23 && pt < 25) {eff = 0.422819; errup = 0.0442283; errdown = 0.0431113;}
//    else if (pt > 25 && pt < 30) {eff = 0.661404; errup = 0.0203857; errdown = 0.0209765;}
//    else if (pt > 30 && pt < 40) {eff = 0.773057; errup = 0.00964611; errdown = 0.00993667;}
//    else if (pt > 40 && pt < 50) {eff = 0.817049; errup = 0.00871766; errdown = 0.0090428;}
//    else if (pt > 50 && pt < 60) {eff = 0.818491; errup = 0.0101637; errdown = 0.0106112;}
//    else if (pt > 60 && pt < 120) {eff = 0.82084; errup = 0.010623; errdown = 0.011121;}
//    else if (pt > 120 && pt < 9999) {eff = 0.824324; errup = 0.0457325; errdown = 0.0557456;}
//  }
//  return std::vector<float>({eff, errup, errdown});
//}

std::vector<float> effsig_elhltpt23(float pt, float abseta) {
  float errup=0.0, errdown=0.0, eff=1.0;
  if (abseta > 0 && abseta < 0.8) {
    if (pt > 7 && pt < 15) {eff = 0.00951777; errup = 0.00312646; errdown = 0.00242343;}
    else if (pt > 15 && pt < 20) {eff = 0.00164948; errup = 0.0013023; errdown = 0.000789256;}
    else if (pt > 20 && pt < 25) {eff = 0.46317; errup = 0.00782595; errdown = 0.00780813;}
    else if (pt > 25 && pt < 27) {eff = 0.976713; errup = 0.00319091; errdown = 0.00364777;}
    else if (pt > 27 && pt < 30) {eff = 0.994443; errup = 0.00117365; errdown = 0.00145015;}
    else if (pt > 30 && pt < 35) {eff = 0.997128; errup = 0.000594034; errdown = 0.000731193;}
    else if (pt > 35 && pt < 40) {eff = 0.999644; errup = 0.00019358; errdown = 0.000345881;}
    else if (pt > 40 && pt < 50) {eff = 0.999467; errup = 0.000184502; errdown = 0.000262852;}
    else if (pt > 50 && pt < 80) {eff = 0.999347; errup = 0.000185925; errdown = 0.000248163;}
    else if (pt > 80 && pt < 120) {eff = 1; errup = 0; errdown = 0.000653326;}
    else if (pt > 120 && pt < 9999) {eff = 0.998937; errup = 0.000879131; errdown = 0.00243942;}
  } else if (abseta > 0.8 && abseta < 1.4442) {
    if (pt > 7 && pt < 15) {eff = 0.0105178; errup = 0.00377339; errdown = 0.00287124;}
    else if (pt > 15 && pt < 20) {eff = 0.00493151; errup = 0.00224407; errdown = 0.00161109;}
    else if (pt > 20 && pt < 25) {eff = 0.462442; errup = 0.0095972; errdown = 0.00956996;}
    else if (pt > 25 && pt < 27) {eff = 0.966392; errup = 0.00472495; errdown = 0.00541101;}
    else if (pt > 27 && pt < 30) {eff = 0.989374; errup = 0.0019891; errdown = 0.00239581;}
    else if (pt > 30 && pt < 35) {eff = 0.998071; errup = 0.0005994; errdown = 0.000821967;}
    else if (pt > 35 && pt < 40) {eff = 0.999089; errup = 0.000393618; errdown = 0.000616093;}
    else if (pt > 40 && pt < 50) {eff = 0.999584; errup = 0.000199027; errdown = 0.000328714;}
    else if (pt > 50 && pt < 80) {eff = 0.999485; errup = 0.000204162; errdown = 0.000307324;}
    else if (pt > 80 && pt < 120) {eff = 1; errup = 0; errdown = 0.00101944;}
    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.00373492;}
  } else if (abseta > 1.4442 && abseta < 1.566) {
    if (pt > 7 && pt < 15) {eff = 0.00628931; errup = 0.0143129; errdown = 0.0052034;}
    else if (pt > 15 && pt < 20) {eff = 0.092511; errup = 0.0234313; errdown = 0.0194408;}
    else if (pt > 20 && pt < 25) {eff = 0.576471; errup = 0.0279879; errdown = 0.0284615;}
    else if (pt > 25 && pt < 27) {eff = 0.92638; errup = 0.0205516; errdown = 0.0265437;}
    else if (pt > 27 && pt < 30) {eff = 0.988889; errup = 0.00604109; errdown = 0.0106898;}
    else if (pt > 30 && pt < 35) {eff = 0.997972; errup = 0.00167805; errdown = 0.00464876;}
    else if (pt > 35 && pt < 40) {eff = 0.997988; errup = 0.00166454; errdown = 0.00461147;}
    else if (pt > 40 && pt < 50) {eff = 1; errup = 0; errdown = 0.00213843;}
    else if (pt > 50 && pt < 80) {eff = 1; errup = 0; errdown = 0.00185228;}
    else if (pt > 80 && pt < 120) {eff = 1; errup = 0; errdown = 0.0132521;}
    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.0376284;}
  } else if (abseta > 1.566 && abseta < 2.55) {
    if (pt > 7 && pt < 15) {eff = 0.0118694; errup = 0.00374129; errdown = 0.00292621;}
    else if (pt > 15 && pt < 20) {eff = 0.0294944; errup = 0.00413256; errdown = 0.00366392;}
    else if (pt > 20 && pt < 25) {eff = 0.335416; errup = 0.00877989; errdown = 0.00866959;}
    else if (pt > 25 && pt < 27) {eff = 0.867746; errup = 0.00887882; errdown = 0.00939394;}
    else if (pt > 27 && pt < 30) {eff = 0.960295; errup = 0.00395923; errdown = 0.00435551;}
    else if (pt > 30 && pt < 35) {eff = 0.988755; errup = 0.00153333; errdown = 0.00175507;}
    else if (pt > 35 && pt < 40) {eff = 0.995232; errup = 0.00100729; errdown = 0.00124491;}
    else if (pt > 40 && pt < 50) {eff = 0.998564; errup = 0.000426166; errdown = 0.000575989;}
    else if (pt > 50 && pt < 80) {eff = 0.999889; errup = 9.17737e-05; errdown = 0.000255059;}
    else if (pt > 80 && pt < 120) {eff = 1; errup = 0; errdown = 0.00159822;}
    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.0073664;}
  }
  return std::vector<float>({eff, errup, errdown});
}

std::vector<float> effsig_dielleg12(float pt, float abseta) {
  float errup=0.0, errdown=0.0, eff=1.0;
  if (abseta > 0 && abseta < 0.8) {
    if (pt > 7 && pt < 12) {eff = 0.0711744; errup = 0.0125667; errdown = 0.0109171;}
    else if (pt > 12 && pt < 15) {eff = 0.683812; errup = 0.0185395; errdown = 0.0191184;}
    else if (pt > 15 && pt < 20) {eff = 0.846727; errup = 0.0083855; errdown = 0.00876625;}
    else if (pt > 20 && pt < 25) {eff = 0.874692; errup = 0.00584446; errdown = 0.00608126;}
    else if (pt > 25 && pt < 27) {eff = 0.877791; errup = 0.00800226; errdown = 0.00846247;}
    else if (pt > 27 && pt < 30) {eff = 0.875758; errup = 0.00609084; errdown = 0.00635086;}
    else if (pt > 30 && pt < 35) {eff = 0.900166; errup = 0.00369246; errdown = 0.00381543;}
    else if (pt > 35 && pt < 40) {eff = 0.919858; errup = 0.00324965; errdown = 0.0033718;}
    else if (pt > 40 && pt < 50) {eff = 0.93026; errup = 0.00227767; errdown = 0.00234737;}
    else if (pt > 50 && pt < 80) {eff = 0.933462; errup = 0.00207147; errdown = 0.00213211;}
    else if (pt > 80 && pt < 120) {eff = 0.92729; errup = 0.00574432; errdown = 0.00617602;}
    else if (pt > 120 && pt < 9999) {eff = 0.938776; errup = 0.00993616; errdown = 0.0115554;}
  } else if (abseta > 0.8 && abseta < 1.4442) {
    if (pt > 7 && pt < 12) {eff = 0.0834915; errup = 0.0138542; errdown = 0.0121509;}
    else if (pt > 12 && pt < 15) {eff = 0.588454; errup = 0.0219745; errdown = 0.022318;}
    else if (pt > 15 && pt < 20) {eff = 0.784727; errup = 0.011231; errdown = 0.0116586;}
    else if (pt > 20 && pt < 25) {eff = 0.830604; errup = 0.00818804; errdown = 0.00850629;}
    else if (pt > 25 && pt < 27) {eff = 0.849192; errup = 0.0108422; errdown = 0.011495;}
    else if (pt > 27 && pt < 30) {eff = 0.854414; errup = 0.00808214; errdown = 0.00845975;}
    else if (pt > 30 && pt < 35) {eff = 0.876201; errup = 0.00494905; errdown = 0.00512098;}
    else if (pt > 35 && pt < 40) {eff = 0.896662; errup = 0.00445872; errdown = 0.00463142;}
    else if (pt > 40 && pt < 50) {eff = 0.910021; errup = 0.00316758; errdown = 0.00326934;}
    else if (pt > 50 && pt < 80) {eff = 0.921457; errup = 0.00276793; errdown = 0.00285836;}
    else if (pt > 80 && pt < 120) {eff = 0.926991; errup = 0.00710404; errdown = 0.00776644;}
    else if (pt > 120 && pt < 9999) {eff = 0.927492; errup = 0.0143507; errdown = 0.0172023;}
  } else if (abseta > 1.4442 && abseta < 1.566) {
    if (pt > 7 && pt < 12) {eff = 0.048; errup = 0.0275697; errdown = 0.0188766;}
    else if (pt > 12 && pt < 15) {eff = 0.467391; errup = 0.057499; errdown = 0.0567216;}
    else if (pt > 15 && pt < 20) {eff = 0.563063; errup = 0.0351552; errdown = 0.03576;}
    else if (pt > 20 && pt < 25) {eff = 0.687898; errup = 0.0270649; errdown = 0.0283321;}
    else if (pt > 25 && pt < 27) {eff = 0.662338; errup = 0.0400803; errdown = 0.0423587;}
    else if (pt > 27 && pt < 30) {eff = 0.780876; errup = 0.0268558; errdown = 0.0292587;}
    else if (pt > 30 && pt < 35) {eff = 0.769392; errup = 0.019723; errdown = 0.0209124;}
    else if (pt > 35 && pt < 40) {eff = 0.802905; errup = 0.0184795; errdown = 0.0198067;}
    else if (pt > 40 && pt < 50) {eff = 0.810875; errup = 0.0136661; errdown = 0.0144328;}
    else if (pt > 50 && pt < 80) {eff = 0.851685; errup = 0.0120587; errdown = 0.0128858;}
    else if (pt > 80 && pt < 120) {eff = 0.833333; errup = 0.0359043; errdown = 0.0424496;}
    else if (pt > 120 && pt < 9999) {eff = 0.878788; errup = 0.0570897; errdown = 0.0855134;}
  } else if (abseta > 1.566 && abseta < 2.55) {
    if (pt > 7 && pt < 12) {eff = 0.0615385; errup = 0.00985507; errdown = 0.00864789;}
    else if (pt > 12 && pt < 15) {eff = 0.484195; errup = 0.0196703; errdown = 0.0196232;}
    else if (pt > 15 && pt < 20) {eff = 0.770238; errup = 0.0103933; errdown = 0.0107242;}
    else if (pt > 20 && pt < 25) {eff = 0.84328; errup = 0.00782547; errdown = 0.00814735;}
    else if (pt > 25 && pt < 27) {eff = 0.889524; errup = 0.0097611; errdown = 0.0105391;}
    else if (pt > 27 && pt < 30) {eff = 0.896901; errup = 0.00747555; errdown = 0.00796744;}
    else if (pt > 30 && pt < 35) {eff = 0.906859; errup = 0.00464458; errdown = 0.00485615;}
    else if (pt > 35 && pt < 40) {eff = 0.918054; errup = 0.00440738; errdown = 0.00462771;}
    else if (pt > 40 && pt < 50) {eff = 0.932145; errup = 0.00316855; errdown = 0.00330833;}
    else if (pt > 50 && pt < 80) {eff = 0.952956; errup = 0.00251428; errdown = 0.00264505;}
    else if (pt > 80 && pt < 120) {eff = 0.942824; errup = 0.00796593; errdown = 0.00907471;}
    else if (pt > 120 && pt < 9999) {eff = 0.926136; errup = 0.0198191; errdown = 0.025339;}
  }
  return std::vector<float>({eff, errup, errdown});
}

std::vector<float> effsig_muhltpt17(float pt, float abseta) {
  float errup=0.0, errdown=0.0, eff=1.0;
  if (abseta > 0 && abseta < 0.9) {
    if (pt > 5 && pt < 17) {eff = 0.00527937; errup = 0.00199877; errdown = 0.00150053;}
    else if (pt > 17 && pt < 20) {eff = 0.982598; errup = 0.00325012; errdown = 0.00390704;}
    else if (pt > 20 && pt < 23) {eff = 0.999118; errup = 0.000569814; errdown = 0.0011624;}
    else if (pt > 23 && pt < 25) {eff = 1; errup = 0; errdown = 0.000973097;}
    else if (pt > 25 && pt < 30) {eff = 0.999848; errup = 0.000125455; errdown = 0.000348643;}
    else if (pt > 30 && pt < 40) {eff = 1; errup = 0; errdown = 9.83561e-05;}
    else if (pt > 40 && pt < 50) {eff = 1; errup = 0; errdown = 0.000105612;}
    else if (pt > 50 && pt < 60) {eff = 1; errup = 0; errdown = 0.00015088;}
    else if (pt > 60 && pt < 120) {eff = 1; errup = 0; errdown = 0.00015908;}
    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.00186731;}
  } else if (abseta > 0.9 && abseta < 1.2) {
    if (pt > 5 && pt < 17) {eff = 0.0111111; errup = 0.00442637; errdown = 0.00328959;}
    else if (pt > 17 && pt < 20) {eff = 0.974955; errup = 0.00656785; errdown = 0.0084952;}
    else if (pt > 20 && pt < 23) {eff = 0.995726; errup = 0.00232491; errdown = 0.00413945;}
    else if (pt > 23 && pt < 25) {eff = 0.998201; errup = 0.0014879; errdown = 0.00412358;}
    else if (pt > 25 && pt < 30) {eff = 1; errup = 0; errdown = 0.000991438;}
    else if (pt > 30 && pt < 40) {eff = 0.998884; errup = 0.000442774; errdown = 0.000666253;}
    else if (pt > 40 && pt < 50) {eff = 1; errup = 0; errdown = 0.000381565;}
    else if (pt > 50 && pt < 60) {eff = 0.99885; errup = 0.000550493; errdown = 0.000908684;}
    else if (pt > 60 && pt < 120) {eff = 1; errup = 0; errdown = 0.000574078;}
    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.00682058;}
  } else if (abseta > 1.2 && abseta < 2.1) {
    if (pt > 5 && pt < 17) {eff = 0.0119974; errup = 0.00229779; errdown = 0.00195515;}
    else if (pt > 17 && pt < 20) {eff = 0.955326; errup = 0.00543073; errdown = 0.00609611;}
    else if (pt > 20 && pt < 23) {eff = 0.996339; errup = 0.00134918; errdown = 0.00196635;}
    else if (pt > 23 && pt < 25) {eff = 0.998588; errup = 0.000912247; errdown = 0.00185985;}
    else if (pt > 25 && pt < 30) {eff = 1; errup = 0; errdown = 0.000425973;}
    else if (pt > 30 && pt < 40) {eff = 0.999385; errup = 0.000226856; errdown = 0.000331209;}
    else if (pt > 40 && pt < 50) {eff = 0.9999; errup = 8.27744e-05; errdown = 0.000230053;}
    else if (pt > 50 && pt < 60) {eff = 1; errup = 0; errdown = 0.000260623;}
    else if (pt > 60 && pt < 120) {eff = 0.999678; errup = 0.000207751; errdown = 0.000424071;}
    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.00433262;}
  } else if (abseta > 2.1 && abseta < 2.4) {
    if (pt > 5 && pt < 17) {eff = 0.0225641; errup = 0.00582471; errdown = 0.00474025;}
    else if (pt > 17 && pt < 20) {eff = 0.926209; errup = 0.0132816; errdown = 0.0156563;}
    else if (pt > 20 && pt < 23) {eff = 0.993711; errup = 0.00342095; errdown = 0.00607987;}
    else if (pt > 23 && pt < 25) {eff = 1; errup = 0; errdown = 0.00559715;}
    else if (pt > 25 && pt < 30) {eff = 0.997899; errup = 0.00135683; errdown = 0.00276412;}
    else if (pt > 30 && pt < 40) {eff = 0.999085; errup = 0.000591198; errdown = 0.00120598;}
    else if (pt > 40 && pt < 50) {eff = 1; errup = 0; errdown = 0.00103317;}
    else if (pt > 50 && pt < 60) {eff = 0.999107; errup = 0.000738625; errdown = 0.00205013;}
    else if (pt > 60 && pt < 120) {eff = 1; errup = 0; errdown = 0.00196917;}
    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.042887;}
  }
  return std::vector<float>({eff, errup, errdown});
}

std::vector<float> effsig_dimuleg8(float pt, float abseta) {
  float errup=0.0, errdown=0.0, eff=1.0;
  if (abseta > 0 && abseta < 0.9) {
    if (pt > 5 && pt < 7) {eff = 0.0491803; errup = 0.045533; errdown = 0.026647;}
    else if (pt > 7 && pt < 10) {eff = 0.693694; errup = 0.0321845; errdown = 0.0340512;}
    else if (pt > 10 && pt < 15) {eff = 0.937788; errup = 0.00823984; errdown = 0.00931878;}
    else if (pt > 15 && pt < 20) {eff = 0.942036; errup = 0.00526642; errdown = 0.00573331;}
    else if (pt > 20 && pt < 23) {eff = 0.950636; errup = 0.0049023; errdown = 0.00538473;}
    else if (pt > 23 && pt < 25) {eff = 0.93836; errup = 0.00590895; errdown = 0.00646014;}
    else if (pt > 25 && pt < 30) {eff = 0.948921; errup = 0.00295975; errdown = 0.00312614;}
    else if (pt > 30 && pt < 40) {eff = 0.946839; errup = 0.00172322; errdown = 0.00177667;}
    else if (pt > 40 && pt < 50) {eff = 0.946058; errup = 0.00180795; errdown = 0.0018659;}
    else if (pt > 50 && pt < 60) {eff = 0.947856; errup = 0.00214365; errdown = 0.0022284;}
    else if (pt > 60 && pt < 120) {eff = 0.941776; errup = 0.00238021; errdown = 0.00247307;}
    else if (pt > 120 && pt < 9999) {eff = 0.879455; errup = 0.0106474; errdown = 0.0114819;}
  } else if (abseta > 0.9 && abseta < 1.2) {
    if (pt > 5 && pt < 7) {eff = 0.0588235; errup = 0.0723694; errdown = 0.0379027;}
    else if (pt > 7 && pt < 10) {eff = 0.691057; errup = 0.0438473; errdown = 0.0472391;}
    else if (pt > 10 && pt < 15) {eff = 0.93765; errup = 0.0118987; errdown = 0.0142071;}
    else if (pt > 15 && pt < 20) {eff = 0.946299; errup = 0.00862063; errdown = 0.0100245;}
    else if (pt > 20 && pt < 23) {eff = 0.954003; errup = 0.0086629; errdown = 0.0103606;}
    else if (pt > 23 && pt < 25) {eff = 0.95122; errup = 0.010162; errdown = 0.0123788;}
    else if (pt > 25 && pt < 30) {eff = 0.958819; errup = 0.00521766; errdown = 0.00588878;}
    else if (pt > 30 && pt < 40) {eff = 0.960847; errup = 0.00284251; errdown = 0.00304713;}
    else if (pt > 40 && pt < 50) {eff = 0.964251; errup = 0.00289055; errdown = 0.0031242;}
    else if (pt > 50 && pt < 60) {eff = 0.965795; errup = 0.00333384; errdown = 0.00366194;}
    else if (pt > 60 && pt < 120) {eff = 0.959878; errup = 0.00384465; errdown = 0.00421355;}
    else if (pt > 120 && pt < 9999) {eff = 0.893281; errup = 0.0196734; errdown = 0.0230906;}
  } else if (abseta > 1.2 && abseta < 2.1) {
    if (pt > 5 && pt < 7) {eff = 0.0231214; errup = 0.0179039; errdown = 0.0110334;}
    else if (pt > 7 && pt < 10) {eff = 0.631455; errup = 0.0241941; errdown = 0.0248413;}
    else if (pt > 10 && pt < 15) {eff = 0.918964; errup = 0.00793836; errdown = 0.00867513;}
    else if (pt > 15 && pt < 20) {eff = 0.916047; errup = 0.00642864; errdown = 0.00688902;}
    else if (pt > 20 && pt < 23) {eff = 0.932577; errup = 0.00632511; errdown = 0.00689737;}
    else if (pt > 23 && pt < 25) {eff = 0.940818; errup = 0.00699255; errdown = 0.00780746;}
    else if (pt > 25 && pt < 30) {eff = 0.933562; errup = 0.00422985; errdown = 0.00448657;}
    else if (pt > 30 && pt < 40) {eff = 0.929039; errup = 0.00256421; errdown = 0.00265102;}
    else if (pt > 40 && pt < 50) {eff = 0.932984; errup = 0.00267392; errdown = 0.00277455;}
    else if (pt > 50 && pt < 60) {eff = 0.94372; errup = 0.00296832; errdown = 0.00311886;}
    else if (pt > 60 && pt < 120) {eff = 0.937291; errup = 0.00340937; errdown = 0.00358625;}
    else if (pt > 120 && pt < 9999) {eff = 0.876214; errup = 0.0164531; errdown = 0.0184132;}
  } else if (abseta > 2.1 && abseta < 2.4) {
    if (pt > 5 && pt < 7) {eff = 0.0166667; errup = 0.0372819; errdown = 0.0137916;}
    else if (pt > 7 && pt < 10) {eff = 0.538462; errup = 0.0470512; errdown = 0.0476918;}
    else if (pt > 10 && pt < 15) {eff = 0.910526; errup = 0.0147857; errdown = 0.0171368;}
    else if (pt > 15 && pt < 20) {eff = 0.902394; errup = 0.0128639; errdown = 0.0144515;}
    else if (pt > 20 && pt < 23) {eff = 0.918367; errup = 0.0139469; errdown = 0.0162748;}
    else if (pt > 23 && pt < 25) {eff = 0.892193; errup = 0.019165; errdown = 0.0223603;}
    else if (pt > 25 && pt < 30) {eff = 0.909091; errup = 0.0104467; errdown = 0.011575;}
    else if (pt > 30 && pt < 40) {eff = 0.903766; errup = 0.00678699; errdown = 0.00722559;}
    else if (pt > 40 && pt < 50) {eff = 0.912961; errup = 0.0071497; errdown = 0.00769771;}
    else if (pt > 50 && pt < 60) {eff = 0.911379; errup = 0.00947284; errdown = 0.0104248;}
    else if (pt > 60 && pt < 120) {eff = 0.894309; errup = 0.0114285; errdown = 0.0125602;}
    else if (pt > 120 && pt < 9999) {eff = 0.85; errup = 0.0578397; errdown = 0.0788355;}
  }
  return std::vector<float>({eff, errup, errdown});
}

std::vector<float> effsig_isoel3235(float pt, float abseta) {
  float errup=0.0, errdown=0.0, eff=1.0;
  if (abseta > 0 && abseta < 0.8) {
    if (pt > 7 && pt < 12) {eff = 0.00214592; errup = 0.00491716; errdown = 0.00177528;}
    else if (pt > 12 && pt < 15) {eff = 0.0118812; errup = 0.00702898; errdown = 0.0047029;}
    else if (pt > 15 && pt < 20) {eff = 0.00938478; errup = 0.00425662; errdown = 0.00306288;}
    else if (pt > 20 && pt < 25) {eff = 0.0122324; errup = 0.00460886; errdown = 0.00347059;}
    else if (pt > 25 && pt < 27) {eff = 0.00769231; errup = 0.00742581; errdown = 0.00418357;}
    else if (pt > 27 && pt < 30) {eff = 0.0308008; errup = 0.00997287; errdown = 0.00779584;}
    else if (pt > 30 && pt < 35) {eff = 0.477987; errup = 0.0206124; errdown = 0.0205405;}
    else if (pt > 35 && pt < 40) {eff = 0.784708; errup = 0.0188249; errdown = 0.0200317;}
    else if (pt > 40 && pt < 50) {eff = 0.807626; errup = 0.0167052; errdown = 0.0178266;}
    else if (pt > 50 && pt < 80) {eff = 0.843267; errup = 0.0173689; errdown = 0.0189802;}
    else if (pt > 80 && pt < 120) {eff = 0.875; errup = 0.0376679; errdown = 0.0485347;}
    else if (pt > 120 && pt < 9999) {eff = 1; errup = 0; errdown = 0.142229;}
  } else if (abseta > 0.8 && abseta < 1.4442) {
    if (pt > 7 && pt < 12) {eff = 0.00518135; errup = 0.0118131; errdown = 0.00428665;}
    else if (pt > 12 && pt < 15) {eff = 0.0045045; errup = 0.0102815; errdown = 0.00372664;}
    else if (pt > 15 && pt < 20) {eff = 0.0100806; errup = 0.00676189; errdown = 0.00434749;}
    else if (pt > 20 && pt < 25) {eff = 0.00539084; errup = 0.0042419; errdown = 0.00257824;}
    else if (pt > 25 && pt < 27) {eff = 0.0217391; errup = 0.011513; errdown = 0.00798301;}
    else if (pt > 27 && pt < 30) {eff = 0.0457666; errup = 0.012319; errdown = 0.0100053;}
    else if (pt > 30 && pt < 35) {eff = 0.469188; errup = 0.0194074; errdown = 0.019318;}
    else if (pt > 35 && pt < 40) {eff = 0.757329; errup = 0.0176647; errdown = 0.0185415;}
    else if (pt > 40 && pt < 50) {eff = 0.764576; errup = 0.0143134; errdown = 0.0149179;}
    else if (pt > 50 && pt < 80) {eff = 0.824934; errup = 0.0140459; errdown = 0.0149484;}
    else if (pt > 80 && pt < 120) {eff = 0.892308; errup = 0.0275795; errdown = 0.0344342;}
    else if (pt > 120 && pt < 9999) {eff = 0.911765; errup = 0.0476324; errdown = 0.0784416;}
  } else if (abseta > 1.4442 && abseta < 1.566) {
    if (pt > 7 && pt < 12) {eff = 0; errup = 0.0472931; errdown = 0;}
    else if (pt > 12 && pt < 15) {eff = 0; errup = 0.108691; errdown = 0;}
    else if (pt > 15 && pt < 20) {eff = 0.0192308; errup = 0.0428344; errdown = 0.0159141;}
    else if (pt > 20 && pt < 25) {eff = 0.0140845; errup = 0.0316413; errdown = 0.0116543;}
    else if (pt > 25 && pt < 27) {eff = 0.147059; errup = 0.0872484; errdown = 0.0620071;}
    else if (pt > 27 && pt < 30) {eff = 0.296875; errup = 0.0679783; errdown = 0.060837;}
    else if (pt > 30 && pt < 35) {eff = 0.418605; errup = 0.0596803; errdown = 0.0575964;}
    else if (pt > 35 && pt < 40) {eff = 0.627907; errup = 0.0559374; errdown = 0.0592187;}
    else if (pt > 40 && pt < 50) {eff = 0.692913; errup = 0.0430384; errdown = 0.0463518;}
    else if (pt > 50 && pt < 80) {eff = 0.775194; errup = 0.0381514; errdown = 0.0428444;}
    else if (pt > 80 && pt < 120) {eff = 0.689655; errup = 0.0935426; errdown = 0.108876;}
    else if (pt > 120 && pt < 9999) {eff = 0.916667; errup = 0.0690403; errdown = 0.16652;}
  } else if (abseta > 1.566 && abseta < 2.55) {
    if (pt > 7 && pt < 12) {eff = 0.0298507; errup = 0.0380185; errdown = 0.0192575;}
    else if (pt > 12 && pt < 15) {eff = 0.0441176; errup = 0.0410576; errdown = 0.0239151;}
    else if (pt > 15 && pt < 20) {eff = 0.020202; errup = 0.015685; errdown = 0.00964386;}
    else if (pt > 20 && pt < 25) {eff = 0.0224719; errup = 0.0094347; errdown = 0.00694883;}
    else if (pt > 25 && pt < 27) {eff = 0.0114068; errup = 0.0109711; errdown = 0.00620172;}
    else if (pt > 27 && pt < 30) {eff = 0.0561122; errup = 0.0122922; errdown = 0.0103423;}
    else if (pt > 30 && pt < 35) {eff = 0.320554; errup = 0.0159521; errdown = 0.0155568;}
    else if (pt > 35 && pt < 40) {eff = 0.633367; errup = 0.0155852; errdown = 0.0158602;}
    else if (pt > 40 && pt < 50) {eff = 0.732684; errup = 0.0104338; errdown = 0.0106922;}
    else if (pt > 50 && pt < 80) {eff = 0.791193; errup = 0.00893809; errdown = 0.00922114;}
    else if (pt > 80 && pt < 120) {eff = 0.880682; errup = 0.0175236; errdown = 0.0198604;}
    else if (pt > 120 && pt < 9999) {eff = 0.853333; errup = 0.0418939; errdown = 0.0527554;}
  }
  return std::vector<float>({eff, errup, errdown});
}

std::vector<float> effsig_isomu2427(float pt, float abseta) {
  float errup=0.0, errdown=0.0, eff=1.0;
  if (abseta > 0 && abseta < 0.9) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 0.0189946; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 0; errup = 0.00669653; errdown = 0;}
    else if (pt > 10 && pt < 15) {eff = 0; errup = 0.00199694; errdown = 0;}
    else if (pt > 15 && pt < 20) {eff = 0; errup = 0.00132646; errdown = 0;}
    else if (pt > 20 && pt < 23) {eff = 0; errup = 0.00211631; errdown = 0;}
    else if (pt > 23 && pt < 25) {eff = 0.39576; errup = 0.0216108; errdown = 0.0212271;}
    else if (pt > 25 && pt < 30) {eff = 0.872209; errup = 0.00934876; errdown = 0.00994538;}
    else if (pt > 30 && pt < 40) {eff = 0.907205; errup = 0.00682098; errdown = 0.00728329;}
    else if (pt > 40 && pt < 50) {eff = 0.93441; errup = 0.00790607; errdown = 0.00883853;}
    else if (pt > 50 && pt < 60) {eff = 0.940594; errup = 0.0118099; errdown = 0.0142172;}
    else if (pt > 60 && pt < 120) {eff = 0.931122; errup = 0.0128701; errdown = 0.0152854;}
    else if (pt > 120 && pt < 9999) {eff = 0.939394; errup = 0.0390483; errdown = 0.074402;}
  } else if (abseta > 0.9 && abseta < 1.2) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 0.0683597; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 0; errup = 0.0279261; errdown = 0;}
    else if (pt > 10 && pt < 15) {eff = 0; errup = 0.0105248; errdown = 0;}
    else if (pt > 15 && pt < 20) {eff = 0.00361011; errup = 0.0082522; errdown = 0.00298664;}
    else if (pt > 20 && pt < 23) {eff = 0; errup = 0.00872845; errdown = 0;}
    else if (pt > 23 && pt < 25) {eff = 0.420732; errup = 0.0419836; errdown = 0.0409445;}
    else if (pt > 25 && pt < 30) {eff = 0.854; errup = 0.0160321; errdown = 0.017536;}
    else if (pt > 30 && pt < 40) {eff = 0.893643; errup = 0.0108821; errdown = 0.0118984;}
    else if (pt > 40 && pt < 50) {eff = 0.934783; errup = 0.0105685; errdown = 0.0122779;}
    else if (pt > 50 && pt < 60) {eff = 0.942966; errup = 0.0143272; errdown = 0.0181346;}
    else if (pt > 60 && pt < 120) {eff = 0.938462; errup = 0.014953; errdown = 0.0187519;}
    else if (pt > 120 && pt < 9999) {eff = 0.96; errup = 0.0331137; errdown = 0.0860512;}
  } else if (abseta > 1.2 && abseta < 2.1) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 0.0527078; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 0; errup = 0.0204732; errdown = 0;}
    else if (pt > 10 && pt < 15) {eff = 0; errup = 0.0082586; errdown = 0;}
    else if (pt > 15 && pt < 20) {eff = 0.00221729; errup = 0.0050801; errdown = 0.00183432;}
    else if (pt > 20 && pt < 23) {eff = 0; errup = 0.00449116; errdown = 0;}
    else if (pt > 23 && pt < 25) {eff = 0.374286; errup = 0.027635; errdown = 0.0268784;}
    else if (pt > 25 && pt < 30) {eff = 0.782112; errup = 0.0129264; errdown = 0.0134829;}
    else if (pt > 30 && pt < 40) {eff = 0.843468; errup = 0.00709467; errdown = 0.00735935;}
    else if (pt > 40 && pt < 50) {eff = 0.839983; errup = 0.00769909; errdown = 0.00800225;}
    else if (pt > 50 && pt < 60) {eff = 0.853472; errup = 0.00941047; errdown = 0.00991949;}
    else if (pt > 60 && pt < 120) {eff = 0.861341; errup = 0.00957946; errdown = 0.0101454;}
    else if (pt > 120 && pt < 9999) {eff = 0.883721; errup = 0.0351207; errdown = 0.0454633;}
  } else if (abseta > 2.1 && abseta < 2.4) {
    if (pt > 5 && pt < 7) {eff = 0; errup = 0.264229; errdown = 0;}
    else if (pt > 7 && pt < 10) {eff = 0; errup = 0.264229; errdown = 0;}
    else if (pt > 10 && pt < 15) {eff = 0; errup = 0.0923495; errdown = 0;}
    else if (pt > 15 && pt < 20) {eff = 0; errup = 0.0461088; errdown = 0;}
    else if (pt > 20 && pt < 23) {eff = 0.0142857; errup = 0.0320826; errdown = 0.0118208;}
    else if (pt > 23 && pt < 25) {eff = 0.245283; errup = 0.0734189; errdown = 0.0624139;}
    else if (pt > 25 && pt < 30) {eff = 0.666667; errup = 0.0325597; errdown = 0.0341201;}
    else if (pt > 30 && pt < 40) {eff = 0.746309; errup = 0.0162592; errdown = 0.0169477;}
    else if (pt > 40 && pt < 50) {eff = 0.790875; errup = 0.0147207; errdown = 0.0154897;}
    else if (pt > 50 && pt < 60) {eff = 0.768212; errup = 0.0175209; errdown = 0.0184509;}
    else if (pt > 60 && pt < 120) {eff = 0.786972; errup = 0.017518; errdown = 0.0185791;}
    else if (pt > 120 && pt < 9999) {eff = 0.857143; errup = 0.0551652; errdown = 0.0755709;}
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
  } else if (abseta > 1.566 && abseta < 2.550) {
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
  } else if (abseta > 1.566 && abseta < 2.550) {
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
