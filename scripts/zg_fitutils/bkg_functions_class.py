import ROOT

#### Bernstein
class Bern2Class:
    def __init__(self, x, gauss_mu, cat="", p0=10, p_init=0.3, bond=20, sigma_init=7., step_init=105.):
        self.p0 = ROOT.RooRealVar("b2p0_" + cat, "b2p1_" + cat, p0)
        self.p1 = ROOT.RooRealVar("b2p1_" + cat, "b2p1_" + cat, p_init,-bond, bond)
        self.p2 = ROOT.RooRealVar("b2p2_" + cat, "b2p2_" + cat, p_init,-bond, bond)
        self.sigma = ROOT.RooRealVar("sigma_bern2_" + cat,"sigma_bern2_" + cat,sigma_init,  0.1, 15.)
        self.stepval = ROOT.RooRealVar("stepval_bern2_" + cat, "stepval_bern2_" + cat, step_init, 95., 115.)
        self.pdf = ROOT.RooGaussStepBernstein("bern2_" +cat + "_model", "Bernstein2 (X) gauss " + cat, x,  gauss_mu, self.sigma, self.stepval, ROOT.RooArgList(self.p0,self.p1,self.p2))
 
class Bern3Class:
    def __init__(self, x, gauss_mu, cat="", p0=10, p_init=0.3, bond=20, sigma_init=7., step_init=105.):
        self.p0 = ROOT.RooRealVar("b3p0_" + cat, "b3p1_" + cat, p0)
        self.p1 = ROOT.RooRealVar("b3p1_" + cat, "b3p1_" + cat, p_init,-bond, bond)
        self.p2 = ROOT.RooRealVar("b3p2_" + cat, "b3p2_" + cat, p_init,-bond, bond)
        self.p3 = ROOT.RooRealVar("b3p3_" + cat, "b3p3_" + cat, p_init,-bond, bond)
        self.sigma = ROOT.RooRealVar("sigma_bern2_" + cat,"sigma_bern2_" + cat,sigma_init,  0.1, 15.)
        self.stepval = ROOT.RooRealVar("stepval_bern2_" + cat, "stepval_bern2_" + cat, step_init, 95., 115.)
        self.pdf = ROOT.RooGaussStepBernstein("bern2_" +cat + "_model", "Bernstein3 (X) gauss " + cat, x,  gauss_mu, self.sigma, self.stepval, ROOT.RooArgList(self.p0,self.p1,self.p2,self.p3))

class Bern4Class:
    def __init__(self, x, gauss_mu, cat="", p0=10, p_init=0.3, bond=20, sigma_init=7., step_init=105.):
        self.p0 = ROOT.RooRealVar("b4p0_" + cat, "b4p1_" + cat, p0)
        self.p1 = ROOT.RooRealVar("b4p1_" + cat, "b4p1_" + cat, p_init,-bond, bond)
        self.p2 = ROOT.RooRealVar("b4p2_" + cat, "b4p2_" + cat, p_init,-bond, bond)
        self.p3 = ROOT.RooRealVar("b4p3_" + cat, "b4p3_" + cat, p_init,-bond, bond)
        self.p4 = ROOT.RooRealVar("b4p4_" + cat, "b4p4_" + cat, p_init,-bond, bond)
        self.sigma = ROOT.RooRealVar("sigma_bern2_" + cat,"sigma_bern2_" + cat,sigma_init,  0.1, 15.)
        self.stepval = ROOT.RooRealVar("stepval_bern2_" + cat, "stepval_bern2_" + cat, step_init, 95., 115.)
        self.pdf = ROOT.RooGaussStepBernstein("bern2_" +cat + "_model", "Bernstein4 (X) gauss " + cat, x,  gauss_mu, self.sigma, self.stepval, ROOT.RooArgList(self.p0,self.p1,self.p2, self.p3, self.p4))

class Bern5Class:
    def __init__(self, x, gauss_mu, cat="", p0=10, p_init=0.3, bond=20, sigma_init=7., step_init=105.):
        self.p0 = ROOT.RooRealVar("b5p0_" + cat, "b5p1_" + cat, p0)
        self.p1 = ROOT.RooRealVar("b5p1_" + cat, "b5p1_" + cat, p_init,-bond, bond)
        self.p2 = ROOT.RooRealVar("b5p2_" + cat, "b5p2_" + cat, p_init,-bond, bond)
        self.p3 = ROOT.RooRealVar("b4p3_" + cat, "b4p3_" + cat, p_init,-bond, bond)
        self.p4 = ROOT.RooRealVar("b4p4_" + cat, "b4p4_" + cat, p_init,-bond, bond)
        self.p5 = ROOT.RooRealVar("b4p5_" + cat, "b4p5_" + cat, p_init,-bond, bond)
        self.sigma = ROOT.RooRealVar("sigma_bern2_" + cat,"sigma_bern2_" + cat,sigma_init,  0.1, 15.)
        self.stepval = ROOT.RooRealVar("stepval_bern2_" + cat, "stepval_bern2_" + cat, step_init, 95., 115.)
        self.pdf = ROOT.RooGaussStepBernstein("bern2_" +cat + "_model", "Bernstein5 (X) gauss " + cat, x,  gauss_mu, self.sigma, self.stepval, ROOT.RooArgList(self.p0,self.p1,self.p2, self.p3, self.p4, self.p5))


#### Pow1 series
class Pow1Class:
    def __init__(self, x, gauss_mu, cat="", sigma_init = 7., step_init = 105., p_init = -7.):
        self.t = ROOT.RooRealVar("pow1t_" + cat, "t pow1" + cat, step_init, 95., 115.)
        self.p = ROOT.RooRealVar("pow1p_" + cat, "p1 pow1" + cat, p_init, -11., -1.)
        self.sigma = ROOT.RooRealVar("sigma_pow1_" + cat,"sigma_pow1_"+cat, sigma_init,  1., 15.)
        self.gauss = ROOT.RooGaussian("gaussxpow1_"+cat, "gaussian PDF pow1 " + cat, x, gauss_mu, self.sigma)
        self.step = ROOT.RooGenericPdf("step_pow1_" + cat, "step_pow1_" + cat,\
                                        "( ((@0-@1)*153.85 <0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85/4) ) )*(@0^@2)", ROOT.RooArgList(x,self.t,self.p))
        x.setBins(20000, "cache")
        self.pdf = ROOT.RooFFTConvPdf("pow1_" + cat + "_model", "step pow1 (X) gauss " + cat, x, self.step, self.gauss)
        self.pdf.setBufferFraction(0.5)

#### Exp series
class Exp1Class:
    def __init__(self, x, gauss_mu, cat="", sigma_init = 7., step_init = 105., p_init = -0.02):
        self.sigma = ROOT.RooRealVar("sigma_exp1_" + cat,"sigma_exp1" + cat       ,sigma_init,  1., 15.)
        self.t = ROOT.RooRealVar("exp_t1_" + cat, "t exp1 " + cat, step_init, 90., 115.)
        self.p1 = ROOT.RooRealVar("exp_p1_" + cat, "p1 exp1 " +cat, p_init, -2., 0)
        self.step = ROOT.RooGenericPdf("step_exp1_" + cat, "step_exp1 " + cat, "( ((@0-@1)*153.85<0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*(TMath::Exp(@0*@2))", ROOT.RooArgList(x,self.t,self.p1))
        self.gauss = ROOT.RooGaussian("gaussxexp1_" + cat, "gaussian PDF exp1 " + cat, x, gauss_mu, self.sigma)

        x.setBins(20000, "cache")
        self.pdf = ROOT.RooFFTConvPdf("exp1_"+cat + "_model", "step exp1 (X) gauss " + cat, x, self.step, self.gauss)
        self.pdf.setBufferFraction(0.5)


class Exp2Class:
    def __init__(self, x, gauss_mu, cat="", sigma_init = 7., step_init = 105., p1_init = -0.02, p2_init = -0.02, f_init = 0.5):
        self.sigma = ROOT.RooRealVar("sigma_exp2_" + cat,"sigma_exp2" + cat       ,sigma_init,  1., 15.)
        self.t = ROOT.RooRealVar("exp2_t1_" + cat, "t exp2 " + cat, step_init, 90., 115.)
        self.p1 = ROOT.RooRealVar("exp2_p1_" + cat, "p1 exp2 " +cat, p1_init, -2., 0)
        self.p2 = ROOT.RooRealVar("exp2_p2_" + cat, "p2 exp2 " +cat, p2_init, -3., 0)
        self.f = ROOT.RooRealVar("exp2_f_" + cat, "f exp2 "+ cat, f_init, 0., 1.)
        self.step = ROOT.RooGenericPdf("step_exp2_" + cat, "step_exp2_" + cat, "( ((@0-@1)*153.85<0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*(@4 * TMath::Exp(@0*@2) + (1-@4)*TMath::Exp(@0*@3))", \
                                      ROOT.RooArgList(x,self.t,self.p1,self.p2,self.f))
        self.gauss = ROOT.RooGaussian("gaussxexp2_" + cat, "gaussian PDF exp2 " + cat, x, gauss_mu, self.sigma)

        x.setBins(20000, "cache")
        self.pdf = ROOT.RooFFTConvPdf("exp2_"+cat + "_model", "step exp2 (X) gauss " + cat, x, self.step, self.gauss)
        self.pdf.setBufferFraction(0.5)

#### Laurent series
class Lau1Class:
    def __init__(self, x, gauss_mu, cat="", sigma_init = 7., step_init = 105., p1 = -7, p2 = -6, f_init = 0.5):
        self.sigma = ROOT.RooRealVar("sigma_lau1_" + cat,"sigma_lau1_" + cat       ,sigma_init,  1., 15.)
        self.t = ROOT.RooRealVar("lau1_t_" + cat, "t lau1 " + cat, step_init, 90., 115.)
        self.f = ROOT.RooRealVar("lau1_f_" + cat, "f lau1 " + cat, f_init, 0., 1.)
        self.p1 = ROOT.RooRealVar("lau1_p1_" + cat, "p1 lau1 " + cat, p1)
        self.p2 = ROOT.RooRealVar("lau1_p2_" + cat, "p2 lau1 " + cat, p2)
        self.step = ROOT.RooGenericPdf("step_lau1_" + cat, "step_lau1_" + cat, "( ((@0-@1)*153.85<0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*(@2*@0^(@3) + (1-@2)*@0^(@4))", \
                                        ROOT.RooArgList(x,self.t,self.f, self.p1, self.p2))
        self.gauss = ROOT.RooGaussian("gaussxlau1_" + cat, "gaussian PDF lau1 " + cat, x, gauss_mu, self.sigma)

        x.setBins(20000, "cache")
        self.pdf = ROOT.RooFFTConvPdf("lau1_" +cat+ "_model", "step lau1 (X) gauss " + cat, x, self.step, self.gauss)
        self.pdf.setBufferFraction(0.5)

class Lau2Class:
    def __init__(self, x, gauss_mu, cat="", sigma_init = 7., step_init = 105., p1 = -8, p2 = -7, p3 = -6, f_init = 0.05):
        self.sigma = ROOT.RooRealVar("sigma_lau2_" + cat,"sigma_lau2_" + cat       ,sigma_init,  1., 15.)
        self.t = ROOT.RooRealVar("lau2_t_" + cat, "t lau2 " + cat, step_init, 90., 115.)
        self.f1 = ROOT.RooRealVar("lau2_f1_" + cat, "f1 lau2 " + cat, f_init, 0., 1.)
        self.f2 = ROOT.RooRealVar("lau2_f2_" + cat, "f2 lau2 " + cat, f_init, 0., 1.)
        self.p1 = ROOT.RooRealVar("lau2_p1_" + cat, "p1 lau2 " + cat, p1)
        self.p2 = ROOT.RooRealVar("lau2_p2_" + cat, "p2 lau2 " + cat, p2)
        self.p3 = ROOT.RooRealVar("lau2_p3_" + cat, "p3 lau2 " + cat, p3)
        self.step = ROOT.RooGenericPdf("step_lau2_" + cat, "step_lau2_" + cat, "( ((@0-@1)*153.85<0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*(@2*@0^(@4) + @3*@0^(@5) + (1-@2-@3)*@0^(@6))", \
                                    ROOT.RooArgList(x,self.t, self.f1, self.f2, self.p1, self.p2, self.p3))
        self.gauss = ROOT.RooGaussian("gaussxlau2_" + cat, "gaussian PDF lau2 " + cat, x, gauss_mu, self.sigma)

        x.setBins(20000, "cache")
        self.pdf = ROOT.RooFFTConvPdf("lau2_" +cat+ "_model", "step lau2 (X) gauss " + cat, x, self.step, self.gauss)
        self.pdf.setBufferFraction(0.5)
    
class Lau3Class:
    def __init__(self, x, gauss_mu, cat="", sigma_init = 7., step_init = 105., p1 = -8, p2 = -7, p3 = -6, p4 = -5, f_init = 0.05):
        self.sigma = ROOT.RooRealVar("sigma_lau3_" + cat,"sigma_lau3_" + cat       ,sigma_init,  1., 15.)
        self.t = ROOT.RooRealVar("lau3_t_" + cat, "t lau3 " + cat, step_init, 90., 115.)
        self.f1 = ROOT.RooRealVar("lau3_f1_" + cat, "f1 lau3 " + cat, f_init, 0., 1.)
        self.f2 = ROOT.RooRealVar("lau3_f2_" + cat, "f2 lau3 " + cat, f_init, 0., 1.)
        self.f3 = ROOT.RooRealVar("lau3_f3_" + cat, "f3 lau3 " + cat, f_init, 0., 1.)
        self.p1 = ROOT.RooRealVar("lau3_p1_" + cat, "p1 lau3 " + cat, p1)
        self.p2 = ROOT.RooRealVar("lau3_p2_" + cat, "p2 lau3 " + cat, p2)
        self.p3 = ROOT.RooRealVar("lau3_p3_" + cat, "p3 lau3 " + cat, p3)
        self.p4 = ROOT.RooRealVar("lau3_p4_" + cat, "p4 lau3 " + cat, p4)
        self.step = ROOT.RooGenericPdf("step_lau3_" + cat, "step_lau3_" + cat, "( ((@0-@1)*153.85<0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*(@2*@0^(@5) + @3*@0^(@6) + @4*@0^(@7) + (1-@2-@3-@4)*@0^(@8))", \
                                      ROOT.RooArgList(x,self.t, self.f1, self.f2, self.f3, self.p1, self.p2, self.p3, self.p4))
        self.gauss = ROOT.RooGaussian("gaussxlau3_" + cat, "gaussian PDF lau3 " + cat, x, gauss_mu, self.sigma)

        x.setBins(20000, "cache")
        self.pdf = ROOT.RooFFTConvPdf("lau3_" +cat+ "_model", "step lau3 (X) gauss " + cat, x, self.step, self.gauss)
        self.pdf.setBufferFraction(0.5)
