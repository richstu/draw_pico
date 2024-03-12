import ROOT
from Xc_Minimizer import *
#### For some reason, in order to evaluate RooAbsPdf in the RooChi2Var or RooNLLVar, all the relavent RooRealVar's must be in the same scope

#### Bernstein functions
def Bern2Minization(x, gauss_mu, histogram, method = "Chi2", error = "Poisson", cat = "", p0 = 10., p_init = 0.3, bond = 50., sigma_init = 7., step_init = 105.,\
                     printLevel = -1, eps = 1, offSet = False, strategy = 0):
    b2p0 = ROOT.RooRealVar("b2p0_" + cat, "b2p1_" + cat, p0)
    b2p1 = ROOT.RooRealVar("b2p1_" + cat, "b2p1_" + cat, p_init,-bond,bond)
    b2p2 = ROOT.RooRealVar("b2p2_" + cat, "b2p2_" + cat, p_init,-bond,bond)
    sigma_bern2 = ROOT.RooRealVar("sigma_bern2_" + cat,"sigma_bern2_" + cat       ,sigma_init,  0.1, 15.)
    stepval_bern2 = ROOT.RooRealVar("stepval_bern2_" + cat, "stepval_bern2_" + cat, step_init, 95., 115.)
    bern2_model = ROOT.RooGaussStepBernstein("bern2_" +cat + "_model", "Bernstein2 (X) gauss " + cat, x, gauss_mu, sigma_bern2, stepval_bern2, ROOT.RooArgList(b2p0,b2p1,b2p2))
    if method == "Chi2":
        if error == "Poisson": chi2 = ROOT.RooChi2Var("chi2_bern2_" + cat, "chi2 bern2 " + cat, bern2_model,  histogram, ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
        elif error == "SumW2": chi2 = ROOT.RooChi2Var("chi2_bern2_" + cat, "chi2 bern2 " + cat, bern2_model,  histogram, ROOT.RooFit.DataError(ROOT.RooAbsData.SumwW2))
        Minimizer_Chi2(chi2, -1, 100, False, strategy)
        return [Minimizer_Chi2(chi2, printLevel, eps, offSet, strategy), chi2.getVal()]
    elif method == "NLL":
        nll = ROOT.RooNLLVar("bern2_NLL_" + cat, "bern2_NLL_" + cat, bern2_model, histogram)
        Minimizer_NLL(chi2, -1, 100, False, strategy)
        return [Minimizer_NLL(chi2, printLevel, eps, offSet, strategy), nll.getVal()]


def Bern3Minization(x, gauss_mu, histogram, method = "Chi2", error = "Poisson", cat = "", p0 = 10., p_init = 0.3, bond = 50., sigma_init = 7., step_init = 105.,\
                     printLevel = -1, eps = 1, offSet = False, strategy = 0):
    b3p0 = ROOT.RooRealVar("b3p0_" + cat, "b3p1_" + cat, p0)
    b3p1 = ROOT.RooRealVar("b3p1_" + cat, "b3p1_" + cat, p_init,-bond,bond)
    b3p2 = ROOT.RooRealVar("b3p2_" + cat, "b3p2_" + cat, p_init,-bond,bond)
    b3p3 = ROOT.RooRealVar("b3p3_" + cat, "b3p3_" + cat, p_init,-bond,bond)
    sigma_bern3 = ROOT.RooRealVar("sigma_bern3_" + cat,"sigma_bern3_" + cat       ,sigma_init,  0.1, 15.)
    stepval_bern3 = ROOT.RooRealVar("stepval_bern3_" + cat, "stepval_bern3_" + cat, step_init, 95., 115.)
    bern3_model = ROOT.RooGaussStepBernstein("bern3_" +cat + "_model", "Bernstein3 (X) gauss " + cat, x, gauss_mu, sigma_bern3, stepval_bern3, ROOT.RooArgList(b3p0,b3p1,b3p2,b3p3))
    if method == "Chi2":
        if error == "Poisson": chi2 = ROOT.RooChi2Var("chi2_bern3_" + cat, "chi2 bern3 " + cat, bern3_model,  histogram, ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
        elif error == "SumW2": chi2 = ROOT.RooChi2Var("chi2_bern3_" + cat, "chi2 bern3 " + cat, bern3_model,  histogram, ROOT.RooFit.DataError(ROOT.RooAbsData.SumwW2))
        Minimizer_Chi2(chi2, -1, 100, False, strategy)
        return [Minimizer_Chi2(chi2, printLevel, eps, offSet, strategy), chi2.getVal()]
    elif method == "NLL":
        nll = ROOT.RooNLLVar("bern3_NLL_" + cat, "bern3_NLL_" + cat, bern3_model, histogram)
        Minimizer_NLL(chi2, -1, 100, False, strategy)
        return [Minimizer_NLL(chi2, printLevel, eps, offSet, strategy), nll.getVal()]

def Bern4Minization(x, gauss_mu, histogram, method = "Chi2", error = "Poisson", cat = "", p0 = 10., p_init = 0.3, bond = 50., sigma_init = 7., step_init = 105.,\
                     printLevel = -1, eps = 1, offSet = False, strategy = 0):
    b4p0 = ROOT.RooRealVar("b4p0_" + cat, "b4p1_" + cat, p0)
    b4p1 = ROOT.RooRealVar("b4p1_" + cat, "b4p1_" + cat, p_init,-bond,bond)
    b4p2 = ROOT.RooRealVar("b4p2_" + cat, "b4p2_" + cat, p_init,-bond,bond)
    b4p3 = ROOT.RooRealVar("b4p3_" + cat, "b4p3_" + cat, p_init,-bond,bond)
    b4p4 = ROOT.RooRealVar("b4p4_" + cat, "b4p4_" + cat, p_init,-bond,bond)
    sigma_bern4 = ROOT.RooRealVar("sigma_bern4_" + cat,"sigma_bern4_" + cat       ,sigma_init,  0.1, 15.)
    stepval_bern4 = ROOT.RooRealVar("stepval_bern4_" + cat, "stepval_bern4_" + cat, step_init, 95., 115.)
    bern4_model = ROOT.RooGaussStepBernstein("bern4_" +cat + "_model", "Bernstein4 (X) gauss " + cat, x, gauss_mu, sigma_bern4, stepval_bern4, ROOT.RooArgList(b4p0,b4p1,b4p2,b4p3, b4p4))
    if method == "Chi2":
        if error == "Poisson": chi2 = ROOT.RooChi2Var("chi2_bern4_" + cat, "chi2 bern4 " + cat, bern4_model,  histogram, ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
        elif error == "SumW2": chi2 = ROOT.RooChi2Var("chi2_bern4_" + cat, "chi2 bern4 " + cat, bern4_model,  histogram, ROOT.RooFit.DataError(ROOT.RooAbsData.SumwW2))
        Minimizer_Chi2(chi2, -1, 100, False, strategy)
        return [Minimizer_Chi2(chi2, printLevel, eps, offSet, strategy), chi2.getVal()]
    elif method == "NLL":
        nll = ROOT.RooNLLVar("bern4_NLL_" + cat, "bern4_NLL_" + cat, bern4_model, histogram)
        Minimizer_NLL(chi2, -1, 100, False, strategy)
        return [Minimizer_NLL(chi2, printLevel, eps, offSet, strategy), nll.getVal()]

def Bern5Minization(x, gauss_mu, histogram, method = "Chi2", error = "Poisson", cat = "", p0 = 10., p_init = 0.3, bond = 50., sigma_init = 7., step_init = 105.,\
                     printLevel = -1, eps = 1, offSet = False, strategy = 0):
    b5p0 = ROOT.RooRealVar("b5p0_" + cat, "b5p1_" + cat, p0)
    b5p1 = ROOT.RooRealVar("b5p1_" + cat, "b5p1_" + cat, p_init,-bond,bond)
    b5p2 = ROOT.RooRealVar("b5p2_" + cat, "b5p2_" + cat, p_init,-bond,bond)
    b5p3 = ROOT.RooRealVar("b5p3_" + cat, "b5p3_" + cat, p_init,-bond,bond)
    b5p4 = ROOT.RooRealVar("b5p4_" + cat, "b5p4_" + cat, p_init,-bond,bond)
    b5p5 = ROOT.RooRealVar("b5p5_" + cat, "b5p5_" + cat, p_init,-bond,bond)
    sigma_bern5 = ROOT.RooRealVar("sigma_bern5_" + cat,"sigma_bern5_" + cat       ,sigma_init,  0.1, 15.)
    stepval_bern5 = ROOT.RooRealVar("stepval_bern5_" + cat, "stepval_bern5_" + cat, step_init, 95., 115.)
    bern5_model = ROOT.RooGaussStepBernstein("bern5_" +cat + "_model", "Bernstein5 (X) gauss " + cat, x, gauss_mu, sigma_bern5, stepval_bern5, ROOT.RooArgList(b5p0,b5p1,b5p2,b5p3, b5p4, b5p5))
    if method == "Chi2":
        if error == "Poisson": chi2 = ROOT.RooChi2Var("chi2_bern5_" + cat, "chi2 bern5 " + cat, bern5_model,  histogram, ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
        elif error == "SumW2": chi2 = ROOT.RooChi2Var("chi2_bern5_" + cat, "chi2 bern5 " + cat, bern5_model,  histogram, ROOT.RooFit.DataError(ROOT.RooAbsData.SumwW2))
        Minimizer_Chi2(chi2, -1, 100, False, strategy)
        return [Minimizer_Chi2(chi2, printLevel, eps, offSet, strategy), chi2.getVal()]
    elif method == "NLL":
        nll = ROOT.RooNLLVar("bern5_NLL_" + cat, "bern5_NLL_" + cat, bern5_model, histogram)
        Minimizer_NLL(chi2, -1, 100, False, strategy)
        return [Minimizer_NLL(chi2, printLevel, eps, offSet, strategy), nll.getVal()]

#### Pow series
def Pow1Minization(x, gauss_mu, histogram, method = "Chi2", error = "Poisson", cat = "", sigma_init = 7., step_init = 105., p_init = -7.,\
                    printLevel = -1, eps = 1, offSet = False, strategy = 0):
    pow1t = ROOT.RooRealVar("pow1t_" + cat, "t pow1" + cat, step_init, 95., 115.)
    pow1p = ROOT.RooRealVar("pow1p_" + cat, "p1 pow1" + cat, p_init, -11., -1.)
    sigma_pow1 = ROOT.RooRealVar("sigma_pow1_" + cat,"sigma_pow1_"+cat, sigma_init,  1., 15.)
    gauss_pow1 = ROOT.RooGaussian("gaussxpow1_"+cat, "gaussian PDF pow1 " + cat, x, gauss_mu, sigma_pow1)
    step_pow1 = ROOT.RooGenericPdf("step_pow1_" + cat, "step_pow1_" + cat, "( ((@0-@1)*153.85 <0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85/4) ) )*(@0^@2)", ROOT.RooArgList(x,pow1t,pow1p))
    
    x.setBins(20000, "cache")
    pow1_model = ROOT.RooFFTConvPdf("pow1_" + cat + "_model", "step pow1 (X) gauss " + cat, x, step_pow1, gauss_pow1)
    pow1_model.setBufferFraction(0.5)
    if method == "Chi2":
        if error == "Poisson": chi2 = ROOT.RooChi2Var("chi2_pow1_" + cat, "chi2 pow1 " + cat, pow1_model,  histogram, ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
        elif error == "SumW2": chi2 = ROOT.RooChi2Var("chi2_pow1_" + cat, "chi2 pow1 " + cat, pow1_model,  histogram, ROOT.RooFit.DataError(ROOT.RooAbsData.SumwW2))
        Minimizer_Chi2(chi2, -1, 100, False, strategy)
        return [Minimizer_Chi2(chi2, printLevel, eps, offSet, strategy), chi2.getVal()]
    elif method == "NLL":
        nll = ROOT.RooNLLVar("pow1_NLL_" + cat, "pow1_NLL_" + cat, pow1_model, histogram)
        Minimizer_NLL(chi2, -1, 100, False, strategy)
        return [Minimizer_NLL(chi2, printLevel, eps, offSet, strategy), nll.getVal()]

#### Exp series 
def Exp1Minization(x, mu_gauss, histogram, method = "Chi2", error = "Poisson", cat = "", sigma_init = 7., step_init = 105., p_init = -0.02,\
                    printLevel = -1, eps = 1, offSet = False, strategy = 0):
    sigma_exp1 = ROOT.RooRealVar("sigma_exp1_" + cat,"sigma_exp1" + cat       ,sigma_init,  1., 15.)
    exp1_t = ROOT.RooRealVar("exp_t1_" + cat, "t exp1 " + cat, step_init, 90., 115.)
    exp1_p1 = ROOT.RooRealVar("exp_p1_" + cat, "p1 exp1 " +cat, p_init, -2., 0)
    step_exp1 = ROOT.RooGenericPdf("step_exp1_" + cat, "step_exp1 " + cat, "( ((@0-@1)*153.85<0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*(TMath::Exp(@0*@2))", ROOT.RooArgList(x,exp1_t,exp1_p1))
    gauss_exp1 = ROOT.RooGaussian("gaussxexp1_" + cat, "gaussian PDF exp1 " + cat, x, mu_gauss, sigma_exp1)

    x.setBins(20000, "cache")
    exp1_model = ROOT.RooFFTConvPdf("exp1_"+cat + "_model", "step exp1 (X) gauss " + cat, x, step_exp1, gauss_exp1)
    exp1_model.setBufferFraction(0.5)
    if method == "Chi2":
        if error == "Poisson": chi2 = ROOT.RooChi2Var("chi2_exp1_" + cat, "chi2 exp1 " + cat, exp1_model,  histogram, ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
        elif error == "SumW2": chi2 = ROOT.RooChi2Var("chi2_exp1_" + cat, "chi2 exp1 " + cat, exp1_model,  histogram, ROOT.RooFit.DataError(ROOT.RooAbsData.SumwW2))
        Minimizer_Chi2(chi2, -1, 100, False, strategy)
        return [Minimizer_Chi2(chi2, printLevel, eps, offSet, strategy), chi2.getVal()]
    elif method == "NLL":
        nll = ROOT.RooNLLVar("exp1_NLL_" + cat, "exp1_NLL_" + cat, exp1_model, histogram)
        Minimizer_NLL(chi2, -1, 100, False, strategy)
        return [Minimizer_NLL(chi2, printLevel, eps, offSet, strategy), nll.getVal()]

def Exp2Minization(x, mu_gauss, histogram, method = "Chi2", error = "Poisson", cat = "", sigma_init = 7., step_init = 105., p1_init = -0.02, p2_init = -0.02, f_init = 0.5,\
                    printLevel = -1, eps = 1, offSet = False, strategy = 0):
    sigma_exp2 = ROOT.RooRealVar("sigma_exp2_" + cat,"sigma_exp2" + cat       ,sigma_init,  1., 15.)
    exp2_t = ROOT.RooRealVar("exp_t1_" + cat, "t exp2 " + cat, step_init, 90., 115.)
    exp2_p1 = ROOT.RooRealVar("exp_p1_" + cat, "p1 exp2 " +cat, p1_init, -2., 0)
    exp2_p2 = ROOT.RooRealVar("exp_p2_" + cat, "p2 exp2 " +cat, p2_init, -3., 0)
    exp2_f = ROOT.RooRealVar("exp2_f_" + cat, "f exp2 "+ cat, f_init, 0., 1.)
    step_exp2 = ROOT.RooGenericPdf("step_exp2_" + cat, "step_exp2_" + cat, "( ((@0-@1)*153.85<0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*(@4 * TMath::Exp(@0*@2) + (1-@4)*TMath::Exp(@0*@3))", \
                                      ROOT.RooArgList(x,exp2_t,exp2_p1,exp2_p2,exp2_f))
    gauss_exp2 = ROOT.RooGaussian("gaussxexp2_" + cat, "gaussian PDF exp2 " + cat, x, mu_gauss, sigma_exp2)

    x.setBins(20000, "cache")
    exp2_model = ROOT.RooFFTConvPdf("exp2_" +cat + "_model", "step exp2 (X) gauss " + cat, x, step_exp2, gauss_exp2)
    exp2_model.setBufferFraction(0.5)
    if method == "Chi2":
        if error == "Poisson": chi2 = ROOT.RooChi2Var("chi2_exp2_" + cat, "chi2 exp2 " + cat, exp2_model,  histogram, ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
        elif error == "SumW2": chi2 = ROOT.RooChi2Var("chi2_exp2_" + cat, "chi2 exp2 " + cat, exp2_model,  histogram, ROOT.RooFit.DataError(ROOT.RooAbsData.SumwW2))
        Minimizer_Chi2(chi2, -1, 100, False, strategy)
        return [Minimizer_Chi2(chi2, printLevel, eps, offSet, strategy), chi2.getVal()]
    elif method == "NLL":
        nll = ROOT.RooNLLVar("exp2_NLL_" + cat, "exp2_NLL_" + cat, exp2_model, histogram)
        Minimizer_NLL(chi2, -1, 100, False, strategy)
        return [Minimizer_NLL(chi2, printLevel, eps, offSet, strategy), nll.getVal()]


#### Laurent series
def Lau1Minization(x, mu_gauss, histogram, method = "Chi2", error = "Poisson", cat = "", sigma_init = 7., step_init = 105., p1 = -7, p2 = -6, f_init = 0.5,\
                    printLevel = -1, eps = 1, offSet = False, strategy = 0):
    sigma_lau1_ = ROOT.RooRealVar("sigma_lau1_" + cat,"sigma_lau1_" + cat       ,sigma_init,  1., 15.)
    lau1_t_ = ROOT.RooRealVar("lau1_t_" + cat, "t lau1 " + cat, step_init, 90., 115.)
    lau1_f_ = ROOT.RooRealVar("lau1_f_" + cat, "f lau1 " + cat, f_init, 0., 1.)
    lau1_p1_ = ROOT.RooRealVar("lau1_p1_" + cat, "p1 lau1 " + cat, p1)
    lau1_p2_ = ROOT.RooRealVar("lau1_p2_" + cat, "p2 lau1 " + cat, p2)
    step_lau1_ = ROOT.RooGenericPdf("step_lau1_" + cat, "step_lau1_" + cat, "( ((@0-@1)*153.85<0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*(@2*@0^(@3) + (1-@2)*@0^(@4))", \
                                    ROOT.RooArgList(x,lau1_t_,lau1_f_, lau1_p1_, lau1_p2_))
    gauss_lau1_ = ROOT.RooGaussian("gaussxlau1_" + cat, "gaussian PDF lau1 " + cat, x, mu_gauss, sigma_lau1_)

    x.setBins(20000, "cache")
    lau1_model = ROOT.RooFFTConvPdf("lau1_" +cat+ "_model", "step lau1 (X) gauss " + cat, x, step_lau1_, gauss_lau1_)
    lau1_model.setBufferFraction(0.5)
    if method == "Chi2":
        if error == "Poisson": chi2 = ROOT.RooChi2Var("chi2_lau1_" + cat, "chi2 lau1 " + cat, lau1_model,  histogram, ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
        elif error == "SumW2": chi2 = ROOT.RooChi2Var("chi2_lau1_" + cat, "chi2 lau1 " + cat, lau1_model,  histogram, ROOT.RooFit.DataError(ROOT.RooAbsData.SumwW2))
        Minimizer_Chi2(chi2, -1, 100, False, strategy)
        return [Minimizer_Chi2(chi2, printLevel, eps, offSet, strategy), chi2.getVal()]
    elif method == "NLL":
        nll = ROOT.RooNLLVar("lau1_NLL_" + cat, "lau1_NLL_" + cat, lau1_model, histogram)
        Minimizer_NLL(chi2, -1, 100, False, strategy)
        return [Minimizer_NLL(chi2, printLevel, eps, offSet, strategy), nll.getVal()]


def Lau2Minization(x, mu_gauss, histogram, method = "Chi2", error = "Poisson", cat = "", sigma_init = 7., step_init = 105., p1 = -8, p2 = -7, p3 = -6, f_init = 0.05,\
                    printLevel = -1, eps = 1, offSet = False, strategy = 0):
    sigma_lau2_ = ROOT.RooRealVar("sigma_lau2_" + cat,"sigma_lau2_" + cat       ,sigma_init,  1., 15.)
    lau2_t_ = ROOT.RooRealVar("lau2_t_" + cat, "t lau2 " + cat, step_init, 90., 115.)
    lau2_f1_ = ROOT.RooRealVar("lau2_f1_" + cat, "f1 lau2 " + cat, f_init, 0., 1.)
    lau2_f2_ = ROOT.RooRealVar("lau2_f2_" + cat, "f2 lau2 " + cat, f_init, 0., 1.)
    lau2_p1_ = ROOT.RooRealVar("lau2_p1_" + cat, "p1 lau2 " + cat, p1)
    lau2_p2_ = ROOT.RooRealVar("lau2_p2_" + cat, "p2 lau2 " + cat, p2)
    lau2_p3_ = ROOT.RooRealVar("lau2_p3_" + cat, "p3 lau2 " + cat, p3)
    step_lau2_ = ROOT.RooGenericPdf("step_lau2_" + cat, "step_lau2_" + cat, "( ((@0-@1)*153.85<0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*(@2*@0^(@4) + @3*@0^(@5) + (1-@2-@3)*@0^(@6))", \
                                    ROOT.RooArgList(x,lau2_t_,lau2_f1_, lau2_f2_, lau2_p1_, lau2_p2_, lau2_p3_))
    gauss_lau2_ = ROOT.RooGaussian("gaussxlau2_" + cat, "gaussian PDF lau2 " + cat, x, mu_gauss, sigma_lau2_)

    x.setBins(20000, "cache")
    lau2_model = ROOT.RooFFTConvPdf("lau2_" +cat+ "_model", "step lau2 (X) gauss " + cat, x, step_lau2_, gauss_lau2_)
    lau2_model.setBufferFraction(0.5)
    if method == "Chi2":
        if error == "Poisson": chi2 = ROOT.RooChi2Var("chi2_lau2_" + cat, "chi2 lau2 " + cat, lau2_model,  histogram, ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
        elif error == "SumW2": chi2 = ROOT.RooChi2Var("chi2_lau2_" + cat, "chi2 lau2 " + cat, lau2_model,  histogram, ROOT.RooFit.DataError(ROOT.RooAbsData.SumwW2))
        Minimizer_Chi2(chi2, -1, 100, False, strategy)
        return [Minimizer_Chi2(chi2, printLevel, eps, offSet, strategy), chi2.getVal()]
    elif method == "NLL":
        nll = ROOT.RooNLLVar("lau2_NLL_" + cat, "lau2_NLL_" + cat, lau2_model, histogram)
        Minimizer_NLL(chi2, -1, 100, False, strategy)
        return [Minimizer_NLL(chi2, printLevel, eps, offSet, strategy), nll.getVal()]

def Lau3Minization(x, mu_gauss, histogram, method = "Chi2", error = "Poisson", cat = "", sigma_init = 7., step_init = 105., p1 = -8, p2 = -7, p3 = -6, p4 = -5, f_init = 0.05,\
                    printLevel = -1, eps = 1, offSet = False, strategy = 0):
    sigma_lau3_ = ROOT.RooRealVar("sigma_lau3_" + cat,"sigma_lau3_" + cat       ,sigma_init,  1., 15.)
    lau3_t_ = ROOT.RooRealVar("lau3_t_" + cat, "t lau3 " + cat, step_init, 90., 115.)
    lau3_f1_ = ROOT.RooRealVar("lau3_f1_" + cat, "f1 lau3 " + cat, f_init, 0., 1.)
    lau3_f2_ = ROOT.RooRealVar("lau3_f2_" + cat, "f2 lau3 " + cat, f_init, 0., 1.)
    lau3_f3_ = ROOT.RooRealVar("lau3_f3_" + cat, "f3 lau3 " + cat, f_init, 0., 1.)
    lau3_p1_ = ROOT.RooRealVar("lau3_p1_" + cat, "p1 lau3 " + cat, p1)
    lau3_p2_ = ROOT.RooRealVar("lau3_p2_" + cat, "p2 lau3 " + cat, p2)
    lau3_p3_ = ROOT.RooRealVar("lau3_p3_" + cat, "p3 lau3 " + cat, p3)
    lau3_p4_ = ROOT.RooRealVar("lau3_p4_" + cat, "p4 lau3 " + cat, p4)
    step_lau3_ = ROOT.RooGenericPdf("step_lau3_" + cat, "step_lau3_" + cat, "( ((@0-@1)*153.85<0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*(@2*@0^(@5) + @3*@0^(@6) + @4*@0^(@7) + (1-@2-@3-@4)*@0^(@8))", \
                                      ROOT.RooArgList(x,lau3_t_,lau3_f1_, lau3_f2_, lau3_f3_, lau3_p1_, lau3_p2_, lau3_p3_, lau3_p4_))
    gauss_lau3_ = ROOT.RooGaussian("gaussxlau3_" + cat, "gaussian PDF lau3 " + cat, x, mu_gauss, sigma_lau3_)

    x.setBins(20000, "cache")
    lau3_model = ROOT.RooFFTConvPdf("lau3_" +cat+ "_model", "step lau3 (X) gauss " + cat, x, step_lau3_, gauss_lau3_)
    lau3_model.setBufferFraction(0.5)
    if method == "Chi2":
        if error == "Poisson": chi2 = ROOT.RooChi2Var("chi2_lau3_" + cat, "chi2 lau3 " + cat, lau3_model,  histogram, ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
        elif error == "SumW2": chi2 = ROOT.RooChi2Var("chi2_lau3_" + cat, "chi2 lau3 " + cat, lau3_model,  histogram, ROOT.RooFit.DataError(ROOT.RooAbsData.SumwW2))
        Minimizer_Chi2(chi2, -1, 100, False, strategy)
        return [Minimizer_Chi2(chi2, printLevel, eps, offSet, strategy), chi2.getVal()]
    elif method == "NLL":
        nll = ROOT.RooNLLVar("lau3_NLL_" + cat, "lau3_NLL_" + cat, lau3_model, histogram)
        Minimizer_NLL(chi2, -1, 100, False, strategy)
        return [Minimizer_NLL(chi2, printLevel, eps, offSet, strategy), nll.getVal()]
