#/usr/bin/env python
'''Script for generating plots of corrections used in H->Zgamma analysis
Note, this depends on JSON files in nano2pico as well as a python3 version 
of correction lib being accesible
'''

from argparse import ArgumentParser
from array import array
from correctionlib import CorrectionSet
from os import path
import ROOT

def get_palette_hig21001():
  return [ROOT.TColor.GetColor('#de5a6a'), ROOT.TColor.GetColor('#ffcc66'), 
          ROOT.TColor.GetColor('#c0e684'), ROOT.TColor.GetColor('#64c0e8'), 
          ROOT.TColor.GetColor('#9999cc'), ROOT.TColor.GetColor('#ffccff')]

def get_palette_lines(nbins):
  if (nbins <= 6):
    return [ROOT.TColor.GetColor('#de394d'), ROOT.TColor.GetColor('#ffbd38'), 
            ROOT.TColor.GetColor('#aee82a'), ROOT.TColor.GetColor('#28aee8'), 
            ROOT.TColor.GetColor('#446bcc'), ROOT.TColor.GetColor('#7346cc')]
  if (nbins <= 8):
    return [ROOT.TColor.GetColor('#de394d'), ROOT.TColor.GetColor('#ff9c38'), 
            ROOT.TColor.GetColor('#f6c24b'), ROOT.TColor.GetColor('#aee82a'), 
            ROOT.TColor.GetColor('#28aee8'), ROOT.TColor.GetColor('#446bcc'), 
            ROOT.TColor.GetColor('#7346cc'), ROOT.TColor.GetColor('#ee66ac')]
  return [ROOT.TColor.GetColor('#de394d'), ROOT.TColor.GetColor('#ff9c38'), 
          ROOT.TColor.GetColor('#f6c24b'), ROOT.TColor.GetColor('#aee82a'), 
          ROOT.TColor.GetColor('#1b9e48'), ROOT.TColor.GetColor('#28aee8'), 
          ROOT.TColor.GetColor('#446bcc'), ROOT.TColor.GetColor('#50378f'),
          ROOT.TColor.GetColor('#7346cc'), ROOT.TColor.GetColor('#ee66ac')]

if __name__=='__main__':

  #parse arguments and get set up
  argument_parser = ArgumentParser(prog='plot_corrections.py',
      description='Generates plots of corrections used by H->Zgamma analysis.')
  argument_parser.add_argument('-n','--nano2pico_dir',default='../nano2pico/')
  args = argument_parser.parse_args()
  if not path.isdir(args.nano2pico_dir):
    raise Exception('nano2pico path not found, please specify with -n option')
  ROOT.gStyle.SetOptStat(0)

  #settings for each plot to be made

  filenames = []
  primary_var_binnings = []
  secondary_var_binnings = []
  primary_var_names = []
  secondary_var_names = []
  correction_names = []
  efficiency_names = []
  lumi_values = []
  plot_names = []
  ordering = []
  log_x = []

  #2016APV ele12
  filenames.append(args.nano2pico_dir+'/data/zgamma/2016preVFP_UL/trigeff_ele12.json')
  secondary_var_binnings.append([0,0.8,1.4442,1.566,2.5])
  primary_var_binnings.append([7,11,12,13,14,16,20,25,30,40,50,100,200,500])
  primary_var_names.append('Electron p_{T} [GeV]')
  secondary_var_names.append('Electron |#eta|')
  correction_names.append('trigeff_ele12')
  efficiency_names.append('Ele12_CaloIdL_TrackIdL_IsoVL')
  lumi_values.append('20')
  plot_names.append('trigeff_ele12_2016APV.pdf')
  ordering.append(0)
  log_x.append(True)

  #2016APV ele23
  filenames.append(args.nano2pico_dir+'/data/zgamma/2016preVFP_UL/trigeff_ele23.json')
  secondary_var_binnings.append([0,0.8,1.4442,1.566,2.5])
  primary_var_binnings.append([7,21,22,23,24,25,27,30,35,40,50,100,200,500])
  primary_var_names.append('Electron p_{T} [GeV]')
  secondary_var_names.append('Electron |#eta|')
  correction_names.append('trigeff_ele23')
  efficiency_names.append('Ele23_CaloIdL_TrackIdL_IsoVL')
  lumi_values.append('20')
  plot_names.append('trigeff_ele23_2016APV.pdf')
  ordering.append(0)
  log_x.append(True)

  #2016APV ele27
  filenames.append(args.nano2pico_dir+'/data/zgamma/2016preVFP_UL/trigeff_ele27.json')
  secondary_var_binnings.append([0,0.8,1.4442,1.566,2.5])
  primary_var_binnings.append([7,25,26,27,28,29,31,35,40,50,100,200,500])
  primary_var_names.append('Electron p_{T} [GeV]')
  secondary_var_names.append('Electron |#eta|')
  correction_names.append('trigeff_ele27')
  efficiency_names.append('Ele27_WPTight_Gsf')
  lumi_values.append('20')
  plot_names.append('trigeff_ele27_2016APV.pdf')
  ordering.append(0)
  log_x.append(True)

  #2016APV mu8
  filenames.append(args.nano2pico_dir+'/data/zgamma/2016preVFP_UL/trigeff_mu8.json')
  secondary_var_binnings.append([0,0.9,1.2,2.1,2.4])
  primary_var_binnings.append([5,7.75,8,8.1,8.25,8.5,10,15,40,100,200,500])
  primary_var_names.append('Muon p_{T} [GeV]')
  secondary_var_names.append('Muon |#eta|')
  correction_names.append('trigeff_mu8')
  efficiency_names.append('(Tk)Mu8_TrkIsoVVL')
  lumi_values.append('20')
  plot_names.append('trigeff_mu8_2016APV.pdf')
  ordering.append(0)
  log_x.append(True)

  #2016APV mu17
  filenames.append(args.nano2pico_dir+'/data/zgamma/2016preVFP_UL/trigeff_mu17.json')
  secondary_var_binnings.append([0,0.9,1.2,2.1,2.4])
  primary_var_binnings.append([5,16.75,17,17.1,17.25,18,20,25,40,100,200,500])
  primary_var_names.append('Muon p_{T} [GeV]')
  secondary_var_names.append('Muon |#eta|')
  correction_names.append('trigeff_mu17')
  efficiency_names.append('Mu17_TrkIsoVVL')
  lumi_values.append('20')
  plot_names.append('trigeff_mu17_2016APV.pdf')
  ordering.append(0)
  log_x.append(True)

  #2016APV mu24
  filenames.append(args.nano2pico_dir+'/data/zgamma/2016preVFP_UL/trigeff_mu24.json')
  secondary_var_binnings.append([0,0.9,1.2,2.1,2.4])
  primary_var_binnings.append([5,23,23.75,24,24.25,24.5,26,30,40,60,120,200,500])
  primary_var_names.append('Muon p_{T} [GeV]')
  secondary_var_names.append('Muon |#eta|')
  correction_names.append('trigeff_mu24')
  efficiency_names.append('Iso(Tk)Mu24')
  lumi_values.append('20')
  plot_names.append('trigeff_mu24_2016APV.pdf')
  ordering.append(0)
  log_x.append(True)

  #2016 ele12
  filenames.append(args.nano2pico_dir+'/data/zgamma/2016postVFP_UL/trigeff_ele12.json')
  secondary_var_binnings.append([0,0.8,1.4442,1.566,2.5])
  primary_var_binnings.append([7,11,12,13,14,16,20,25,30,40,50,100,200,500])
  primary_var_names.append('Electron p_{T} [GeV]')
  secondary_var_names.append('Electron |#eta|')
  correction_names.append('trigeff_ele12')
  efficiency_names.append('Ele12_CaloIdL_TrackIdL_IsoVL')
  lumi_values.append('17')
  plot_names.append('trigeff_ele12_2016.pdf')
  ordering.append(0)
  log_x.append(True)

  #2016 ele23
  filenames.append(args.nano2pico_dir+'/data/zgamma/2016postVFP_UL/trigeff_ele23.json')
  secondary_var_binnings.append([0,0.8,1.4442,1.566,2.5])
  primary_var_binnings.append([7,21,22,23,24,25,27,30,35,40,50,100,200,500])
  primary_var_names.append('Electron p_{T} [GeV]')
  secondary_var_names.append('Electron |#eta|')
  correction_names.append('trigeff_ele23')
  efficiency_names.append('Ele23_CaloIdL_TrackIdL_IsoVL')
  lumi_values.append('17')
  plot_names.append('trigeff_ele23_2016.pdf')
  ordering.append(0)
  log_x.append(True)

  #2016 ele27
  filenames.append(args.nano2pico_dir+'/data/zgamma/2016postVFP_UL/trigeff_ele27.json')
  secondary_var_binnings.append([0,0.8,1.4442,1.566,2.5])
  primary_var_binnings.append([7,25,26,27,28,29,31,35,40,50,100,200,500])
  primary_var_names.append('Electron p_{T} [GeV]')
  secondary_var_names.append('Electron |#eta|')
  correction_names.append('trigeff_ele27')
  efficiency_names.append('Ele27_WPTight_Gsf')
  lumi_values.append('17')
  plot_names.append('trigeff_ele27_2016.pdf')
  ordering.append(0)
  log_x.append(True)

  #2016 mu8
  filenames.append(args.nano2pico_dir+'/data/zgamma/2016postVFP_UL/trigeff_mu8.json')
  secondary_var_binnings.append([0,0.9,1.2,2.1,2.4])
  primary_var_binnings.append([5,7.75,8,8.1,8.25,8.5,10,15,40,100,200,500])
  primary_var_names.append('Muon p_{T} [GeV]')
  secondary_var_names.append('Muon |#eta|')
  correction_names.append('trigeff_mu8')
  efficiency_names.append('(Tk)Mu8_TrkIsoVVL')
  lumi_values.append('17')
  plot_names.append('trigeff_mu8_2016.pdf')
  ordering.append(0)
  log_x.append(True)

  #2016 mu17
  filenames.append(args.nano2pico_dir+'/data/zgamma/2016postVFP_UL/trigeff_mu17.json')
  secondary_var_binnings.append([0,0.9,1.2,2.1,2.4])
  primary_var_binnings.append([5,16.75,17,17.1,17.25,18,20,25,40,100,200,500])
  primary_var_names.append('Muon p_{T} [GeV]')
  secondary_var_names.append('Muon |#eta|')
  correction_names.append('trigeff_mu17')
  efficiency_names.append('Mu17_TrkIsoVVL')
  lumi_values.append('17')
  plot_names.append('trigeff_mu17_2016.pdf')
  ordering.append(0)
  log_x.append(True)

  #2016 mu24
  filenames.append(args.nano2pico_dir+'/data/zgamma/2016postVFP_UL/trigeff_mu24.json')
  secondary_var_binnings.append([0,0.9,1.2,2.1,2.4])
  primary_var_binnings.append([5,23,23.75,24,24.25,24.5,26,30,40,60,120,200,500])
  primary_var_names.append('Muon p_{T} [GeV]')
  secondary_var_names.append('Muon |#eta|')
  correction_names.append('trigeff_mu24')
  efficiency_names.append('Iso(Tk)Mu24')
  lumi_values.append('17')
  plot_names.append('trigeff_mu24_2016.pdf')
  ordering.append(0)
  log_x.append(True)

  #2017 ele12
  filenames.append(args.nano2pico_dir+'/data/zgamma/2017_UL/trigeff_ele12.json')
  secondary_var_binnings.append([0,0.8,1.4442,1.566,2.5])
  primary_var_binnings.append([7,11,12,13,14,16,20,25,30,40,50,100,200,500])
  primary_var_names.append('Electron p_{T} [GeV]')
  secondary_var_names.append('Electron |#eta|')
  correction_names.append('trigeff_ele12')
  efficiency_names.append('Ele12_CaloIdL_TrackIdL_IsoVL')
  lumi_values.append('41')
  plot_names.append('trigeff_ele12_2017.pdf')
  ordering.append(0)
  log_x.append(True)

  #2017 ele23
  filenames.append(args.nano2pico_dir+'/data/zgamma/2017_UL/trigeff_ele23.json')
  secondary_var_binnings.append([0,0.8,1.4442,1.566,2.5])
  primary_var_binnings.append([7,21,22,23,24,25,27,30,35,40,50,100,200,500])
  primary_var_names.append('Electron p_{T} [GeV]')
  secondary_var_names.append('Electron |#eta|')
  correction_names.append('trigeff_ele23')
  efficiency_names.append('Ele23_CaloIdL_TrackIdL_IsoVL')
  lumi_values.append('41')
  plot_names.append('trigeff_ele23_2017.pdf')
  ordering.append(0)
  log_x.append(True)

  #2017 ele32
  filenames.append(args.nano2pico_dir+'/data/zgamma/2017_UL/trigeff_ele32.json')
  secondary_var_binnings.append([0,0.8,1.4442,1.566,2.5])
  primary_var_binnings.append([7,31,32,33,34,35,38,45,80,120,200,500])
  primary_var_names.append('Electron p_{T} [GeV]')
  secondary_var_names.append('Electron |#eta|')
  correction_names.append('trigeff_ele32')
  efficiency_names.append('Ele32_WPTight_Gsf_L1DoubleEG')
  lumi_values.append('41')
  plot_names.append('trigeff_ele32_2017.pdf')
  ordering.append(0)
  log_x.append(True)

  #2017 mu8
  filenames.append(args.nano2pico_dir+'/data/zgamma/2017_UL/trigeff_mu8.json')
  secondary_var_binnings.append([0,0.9,1.2,2.1,2.4])
  primary_var_binnings.append([5,7.75,8,8.1,8.25,8.5,10,15,40,100,200,500])
  primary_var_names.append('Muon p_{T} [GeV]')
  secondary_var_names.append('Muon |#eta|')
  correction_names.append('trigeff_mu8')
  efficiency_names.append('Mu8_TrkIsoVVL')
  lumi_values.append('41')
  plot_names.append('trigeff_mu8_2017.pdf')
  ordering.append(0)
  log_x.append(True)

  #2017 mu17
  filenames.append(args.nano2pico_dir+'/data/zgamma/2017_UL/trigeff_mu17.json')
  secondary_var_binnings.append([0,0.9,1.2,2.1,2.4])
  primary_var_binnings.append([5,16.75,17,17.1,17.25,18,20,25,40,100,200,500])
  primary_var_names.append('Muon p_{T} [GeV]')
  secondary_var_names.append('Muon |#eta|')
  correction_names.append('trigeff_mu17')
  efficiency_names.append('Mu17_TrkIsoVVL')
  lumi_values.append('41')
  plot_names.append('trigeff_mu17_2017.pdf')
  ordering.append(0)
  log_x.append(True)

  #2017 mu27
  filenames.append(args.nano2pico_dir+'/data/zgamma/2017_UL/trigeff_mu27.json')
  secondary_var_binnings.append([0,0.9,1.2,2.1,2.4])
  primary_var_binnings.append([5,26.75,27,27.25,27.5,29,32,40,60,120,200,500])
  primary_var_names.append('Muon p_{T} [GeV]')
  secondary_var_names.append('Muon |#eta|')
  correction_names.append('trigeff_mu27')
  efficiency_names.append('IsoMu27')
  lumi_values.append('41')
  plot_names.append('trigeff_mu27_2017.pdf')
  ordering.append(0)
  log_x.append(True)

  #2018 ele12
  filenames.append(args.nano2pico_dir+'/data/zgamma/2018_UL/trigeff_ele12.json')
  secondary_var_binnings.append([0,0.8,1.4442,1.566,2.5])
  primary_var_binnings.append([7,11,12,13,14,16,20,25,30,40,50,100,200,500])
  primary_var_names.append('Electron p_{T} [GeV]')
  secondary_var_names.append('Electron |#eta|')
  correction_names.append('trigeff_ele12')
  efficiency_names.append('Ele12_CaloIdL_TrackIdL_IsoVL')
  lumi_values.append('60')
  plot_names.append('trigeff_ele12_2018.pdf')
  ordering.append(0)
  log_x.append(True)

  #2018 ele23
  filenames.append(args.nano2pico_dir+'/data/zgamma/2018_UL/trigeff_ele23.json')
  secondary_var_binnings.append([0,0.8,1.4442,1.566,2.5])
  primary_var_binnings.append([7,21,22,23,24,25,27,30,35,40,50,100,200,500])
  primary_var_names.append('Electron p_{T} [GeV]')
  secondary_var_names.append('Electron |#eta|')
  correction_names.append('trigeff_ele23')
  efficiency_names.append('Ele23_CaloIdL_TrackIdL_IsoVL')
  lumi_values.append('60')
  plot_names.append('trigeff_ele23_2018.pdf')
  ordering.append(0)
  log_x.append(True)

  #2018 ele32
  filenames.append(args.nano2pico_dir+'/data/zgamma/2018_UL/trigeff_ele32.json')
  secondary_var_binnings.append([0,0.8,1.4442,1.566,2.5])
  primary_var_binnings.append([7,31,32,33,34,35,38,45,80,120,200,500])
  primary_var_names.append('Electron p_{T} [GeV]')
  secondary_var_names.append('Electron |#eta|')
  correction_names.append('trigeff_ele32')
  efficiency_names.append('Ele32_WPTight_Gsf')
  lumi_values.append('60')
  plot_names.append('trigeff_ele32_2018.pdf')
  ordering.append(0)
  log_x.append(True)

  #2018 mu8
  filenames.append(args.nano2pico_dir+'/data/zgamma/2018_UL/trigeff_mu8.json')
  secondary_var_binnings.append([0,0.9,1.2,2.1,2.4])
  primary_var_binnings.append([5,7.75,8,8.1,8.25,8.5,10,15,40,100,200,500])
  primary_var_names.append('Muon p_{T} [GeV]')
  secondary_var_names.append('Muon |#eta|')
  correction_names.append('trigeff_mu8')
  efficiency_names.append('Mu8_TrkIsoVVL')
  lumi_values.append('60')
  plot_names.append('trigeff_mu8_2018.pdf')
  ordering.append(0)
  log_x.append(True)

  #2018 mu17
  filenames.append(args.nano2pico_dir+'/data/zgamma/2018_UL/trigeff_mu17.json')
  secondary_var_binnings.append([0,0.9,1.2,2.1,2.4])
  primary_var_binnings.append([5,16.75,17,17.1,17.25,18,20,25,40,100,200,500])
  primary_var_names.append('Muon p_{T} [GeV]')
  secondary_var_names.append('Muon |#eta|')
  correction_names.append('trigeff_mu17')
  efficiency_names.append('Mu17_TrkIsoVVL')
  lumi_values.append('60')
  plot_names.append('trigeff_mu17_2018.pdf')
  ordering.append(0)
  log_x.append(True)

  #2018 mu24
  filenames.append(args.nano2pico_dir+'/data/zgamma/2018_UL/trigeff_mu24.json')
  secondary_var_binnings.append([0,0.9,1.2,2.1,2.4])
  primary_var_binnings.append([5,23,23.75,24,24.25,24.5,26,30,40,60,120,200,500])
  primary_var_names.append('Muon p_{T} [GeV]')
  secondary_var_names.append('Muon |#eta|')
  correction_names.append('trigeff_mu24')
  efficiency_names.append('IsoMu24')
  lumi_values.append('60')
  plot_names.append('trigeff_mu24_2018.pdf')
  ordering.append(0)
  log_x.append(True)

  #2016APV WPL
  filenames.append(args.nano2pico_dir+'/data/zgamma/2016preVFP_UL/electron_WPL.json')
  secondary_var_binnings.append([-2.5,-2.0,-1.566,-1.4442,-0.8,0.0,0.8,1.442,1.566,2.0,2.5])
  primary_var_binnings.append([7,10,20,35,50,100,200,500])
  primary_var_names.append('Electron p_{T} [GeV]')
  secondary_var_names.append('Electron #eta')
  correction_names.append('ElectronWPL')
  efficiency_names.append('Fall17v2 Iso ID WPL')
  lumi_values.append('20')
  plot_names.append('electron_WPL_2016APV.pdf')
  ordering.append(0)
  log_x.append(True)

  #2016 WPL
  filenames.append(args.nano2pico_dir+'/data/zgamma/2016postVFP_UL/electron_WPL.json')
  secondary_var_binnings.append([-2.5,-2.0,-1.566,-1.4442,-0.8,0.0,0.8,1.442,1.566,2.0,2.5])
  primary_var_binnings.append([7,10,20,35,50,100,200,500])
  primary_var_names.append('Electron p_{T} [GeV]')
  secondary_var_names.append('Electron #eta')
  correction_names.append('ElectronWPL')
  efficiency_names.append('Fall17v2 Iso ID WPL')
  lumi_values.append('17')
  plot_names.append('electron_WPL_2016.pdf')
  ordering.append(0)
  log_x.append(True)

  #2017 WPL
  filenames.append(args.nano2pico_dir+'/data/zgamma/2017_UL/electron_WPL.json')
  secondary_var_binnings.append([-2.5,-2.0,-1.566,-1.4442,-0.8,0.0,0.8,1.442,1.566,2.0,2.5])
  primary_var_binnings.append([7,10,20,35,50,100,200,500])
  primary_var_names.append('Electron p_{T} [GeV]')
  secondary_var_names.append('Electron #eta')
  correction_names.append('ElectronWPL')
  efficiency_names.append('Fall17v2 Iso ID WPL')
  lumi_values.append('41')
  plot_names.append('electron_WPL_2017.pdf')
  ordering.append(0)
  log_x.append(True)

  #2018 WPL
  filenames.append(args.nano2pico_dir+'/data/zgamma/2018_UL/electron_WPL.json')
  secondary_var_binnings.append([-2.5,-2.0,-1.566,-1.4442,-0.8,0.0,0.8,1.442,1.566,2.0,2.5])
  primary_var_binnings.append([7,10,20,35,50,100,200,500])
  primary_var_names.append('Electron p_{T} [GeV]')
  secondary_var_names.append('Electron #eta')
  correction_names.append('ElectronWPL')
  efficiency_names.append('Fall17v2 Iso ID WPL')
  lumi_values.append('60')
  plot_names.append('electron_WPL_2018.pdf')
  ordering.append(0)
  log_x.append(True)

  #2016APV WPL eta
  filenames.append(args.nano2pico_dir+'/data/zgamma/2016preVFP_UL/electron_WPL.json')
  primary_var_binnings.append([-2.5,-2.0,-1.566,-1.4442,-0.8,0.0,0.8,1.442,1.566,2.0,2.5])
  secondary_var_binnings.append([7,10,20,35,50,100,200,500])
  secondary_var_names.append('Electron p_{T} [GeV]')
  primary_var_names.append('Electron #eta')
  correction_names.append('ElectronWPL')
  efficiency_names.append('Fall17v2 Iso ID WPL')
  lumi_values.append('20')
  plot_names.append('electron_WPL_eta_2016APV.pdf')
  ordering.append(1)
  log_x.append(False)

  #2016 WPL eta
  filenames.append(args.nano2pico_dir+'/data/zgamma/2016postVFP_UL/electron_WPL.json')
  primary_var_binnings.append([-2.5,-2.0,-1.566,-1.4442,-0.8,0.0,0.8,1.442,1.566,2.0,2.5])
  secondary_var_binnings.append([7,10,20,35,50,100,200,500])
  secondary_var_names.append('Electron p_{T} [GeV]')
  primary_var_names.append('Electron #eta')
  correction_names.append('ElectronWPL')
  efficiency_names.append('Fall17v2 Iso ID WPL')
  lumi_values.append('17')
  plot_names.append('electron_WPL_eta_2016.pdf')
  ordering.append(1)
  log_x.append(False)

  #2017 WPL eta
  filenames.append(args.nano2pico_dir+'/data/zgamma/2017_UL/electron_WPL.json')
  primary_var_binnings.append([-2.5,-2.0,-1.566,-1.4442,-0.8,0.0,0.8,1.442,1.566,2.0,2.5])
  secondary_var_binnings.append([7,10,20,35,50,100,200,500])
  secondary_var_names.append('Electron p_{T} [GeV]')
  primary_var_names.append('Electron #eta')
  correction_names.append('ElectronWPL')
  efficiency_names.append('Fall17v2 Iso ID WPL')
  lumi_values.append('41')
  plot_names.append('electron_WPL_eta_2017.pdf')
  ordering.append(1)
  log_x.append(False)

  #2018 WPL eta
  filenames.append(args.nano2pico_dir+'/data/zgamma/2018_UL/electron_WPL.json')
  primary_var_binnings.append([-2.5,-2.0,-1.566,-1.4442,-0.8,0.0,0.8,1.442,1.566,2.0,2.5])
  secondary_var_binnings.append([7,10,20,35,50,100,200,500])
  secondary_var_names.append('Electron p_{T} [GeV]')
  primary_var_names.append('Electron #eta')
  correction_names.append('ElectronWPL')
  efficiency_names.append('Fall17v2 Iso ID WPL')
  lumi_values.append('60')
  plot_names.append('electron_WPL_eta_2018.pdf')
  ordering.append(1)
  log_x.append(False)

  for iplot in range(len(filenames)):

    #load file and set up binning
    correction_set = CorrectionSet.from_file(filenames[iplot])
    corr_name = correction_names[iplot]
    var1_bins = primary_var_binnings[iplot]
    var2_bins = secondary_var_binnings[iplot]
    var1_bin_centers = [(var1_bins[ivar1]+var1_bins[ivar1+1])/2.0 for ivar1 in range(len(var1_bins)-1)]
    var2_bin_centers = [(var2_bins[ivar2]+var2_bins[ivar2+1])/2.0 for ivar2 in range(len(var2_bins)-1)]

    #fill x and y values from JSON
    x_vals = array('d')
    x_sup = array('d')
    x_sdn = array('d')
    y_vals = [array('d') for ivar2 in range(len(var2_bin_centers))]
    y_sup = [array('d') for ivar2 in range(len(var2_bin_centers))]
    y_sdn = [array('d') for ivar2 in range(len(var2_bin_centers))]
    ratio_vals = [array('d') for ivar2 in range(len(var2_bin_centers))]
    ratio_sup = [array('d') for ivar2 in range(len(var2_bin_centers))]
    ratio_sdn = [array('d') for ivar2 in range(len(var2_bin_centers))]
    for ivar1 in range(len(var1_bin_centers)):
      var1_bin_center = var1_bin_centers[ivar1]
      x_vals.append(var1_bin_center)
      x_sdn.append(var1_bin_center-var1_bins[ivar1])
      x_sup.append(var1_bins[ivar1+1]-var1_bin_center)
      for ivar2 in range(len(var2_bin_centers)):
        var2_bin_center = var2_bin_centers[ivar2]
        data_eff = 0.0
        data_unc = 0.0
        simu_eff = 0.0
        simu_unc = 0.0
        if ordering[iplot]==0:
          data_eff = correction_set[corr_name].evaluate('effdata',
              var2_bin_center,var1_bin_center)
          data_unc = correction_set[corr_name].evaluate('systdata',
              var2_bin_center,var1_bin_center)
          simu_eff = correction_set[corr_name].evaluate('effmc',
              var2_bin_center,var1_bin_center)
          simu_unc = correction_set[corr_name].evaluate('systmc',
              var2_bin_center,var1_bin_center)
        elif ordering[iplot]==1:
          data_eff = correction_set[corr_name].evaluate('effdata',
              var1_bin_center,var2_bin_center)
          data_unc = correction_set[corr_name].evaluate('systdata',
              var1_bin_center,var2_bin_center)
          simu_eff = correction_set[corr_name].evaluate('effmc',
              var1_bin_center,var2_bin_center)
          simu_unc = correction_set[corr_name].evaluate('systmc',
              var1_bin_center,var2_bin_center)
        data_sup = min(data_unc,1.0-data_eff)
        data_sdn = min(data_unc,data_eff)
        simu_sup = min(simu_unc,1.0-simu_eff)
        simu_sdn = min(simu_unc,simu_eff)
        y_vals[ivar2].append(data_eff)
        y_sup[ivar2].append(data_sup)
        y_sdn[ivar2].append(data_sdn)
        ratio = 0.0
        if (simu_eff != 0):
          ratio = data_eff/simu_eff
        ratio_up = ratio
        ratio_dn = ratio
        if (simu_eff-simu_sdn != 0) and (simu_eff+simu_sup != 0):
          ratio_up = (data_eff+data_sup)/(simu_eff-simu_sdn)
          ratio_dn = (data_eff-data_sdn)/(simu_eff+simu_sup)
        ratio_vals[ivar2].append(ratio)
        ratio_sup[ivar2].append(ratio_up-ratio)
        ratio_sdn[ivar2].append(ratio-ratio_dn)

    #create and style graphs
    upper_graphs = []
    lower_graphs = []
    colors = get_palette_lines(len(var2_bin_centers))
    for ivar2 in range(len(var2_bin_centers)):
      upper_graphs.append(ROOT.TGraphAsymmErrors(len(x_vals),x_vals,
        y_vals[ivar2],x_sdn,x_sup,y_sdn[ivar2],y_sup[ivar2]))
      upper_graphs[-1].SetLineWidth(3)
      upper_graphs[-1].SetLineColor(colors[ivar2])
      upper_graphs[-1].SetMarkerColor(colors[ivar2])
      upper_graphs[-1].SetMarkerStyle(ROOT.kFullCircle)

      lower_graphs.append(ROOT.TGraphAsymmErrors(len(x_vals),x_vals,
        ratio_vals[ivar2],x_sdn,x_sup,ratio_sdn[ivar2],ratio_sup[ivar2]))
      lower_graphs[-1].SetLineWidth(3)
      lower_graphs[-1].SetLineColor(colors[ivar2])
      lower_graphs[-1].SetMarkerColor(colors[ivar2])
      lower_graphs[-1].SetMarkerStyle(ROOT.kFullCircle)

    #do plotting
    dummy_hist_upper = ROOT.TH1D('','',100,var1_bins[0],var1_bins[-1])
    dummy_hist_upper.SetMinimum(0.0)
    dummy_hist_upper.SetMaximum(1.2)
    dummy_hist_upper.SetLabelSize(0,'x')
    dummy_hist_upper.SetLabelSize(0.028,'y')
    dummy_hist_upper.SetTitleSize(0.025,'y')
    dummy_hist_upper.GetYaxis().SetTitle(efficiency_names[iplot]+' Efficiency (Data)')
    dummy_hist_upper.GetYaxis().SetNdivisions(606)
    dummy_hist_upper.GetXaxis().SetNdivisions(606)

    dummy_hist_lower = ROOT.TH1D('','',100,var1_bins[0],var1_bins[-1])
    dummy_hist_lower.SetMinimum(0.25)
    dummy_hist_lower.SetMaximum(1.75)
    dummy_hist_lower.SetLabelSize(0.028,'x')
    dummy_hist_lower.SetLabelSize(0.028,'y')
    dummy_hist_lower.SetTitleSize(0.025,'y')
    dummy_hist_lower.GetYaxis().SetTitle('Efficiency Data/MC')
    dummy_hist_lower.GetYaxis().SetNdivisions(502)
    dummy_hist_lower.GetYaxis().SetLimits(0.3,1.6)
    dummy_hist_lower.GetXaxis().SetTitle(primary_var_names[iplot])
    dummy_hist_lower.GetXaxis().SetNdivisions(606)

    can = ROOT.TCanvas('c','c',600,600)
    top_pad = ROOT.TPad('top_pad','',0.0,0.0,1.0,1.0)
    top_pad.SetTicks(1,1)
    top_pad.SetMargin(0.14,0.06,0.31,0.05)
    top_pad.SetFillStyle(4000)
    top_pad.SetLogx(log_x[iplot])

    bot_pad = ROOT.TPad('bot_pad','',0.0,0.0,1.0,1.0)
    bot_pad.SetTicks(1,1)
    bot_pad.SetMargin(0.14,0.06,0.1,0.71)
    bot_pad.SetFillStyle(4000)
    bot_pad.SetLogx(log_x[iplot])

    top_pad.Draw()
    top_pad.cd()
    dummy_hist_upper.Draw()
    leg = ROOT.TLegend(0.6,0.35,0.9,0.55)
    leg.SetEntrySeparation(0)
    for ivar2 in range(len(var2_bin_centers)):
      upper_graphs[ivar2].Draw('same P')
      var_nounit = secondary_var_names[iplot]
      var_unit = ''
      if var_nounit.find('[') != -1:
        var_unit = ' '+var_nounit[var_nounit.find('[')+1:var_nounit.find(']')]
        var_nounit = var_nounit[:var_nounit.find('[')-1]
      eta_text = (str(var2_bins[ivar2])+' < '+var_nounit
          +' < '+str(var2_bins[ivar2+1])+var_unit)
      leg.AddEntry(upper_graphs[ivar2], eta_text, 'LP')
    leg.SetBorderSize(0)
    leg.Draw('same')
    label = ROOT.TLatex()
    label.SetTextSize(0.035)
    label.SetNDC(ROOT.kTRUE)
    label.SetTextAlign(11)
    label.DrawLatex(0.15,0.96,'#font[62]{CMS} #scale[0.8]{#font[52]{Preliminary}}')
    label.SetTextAlign(31)
    label.SetTextSize(0.03)
    label.DrawLatex(0.93,0.96,'#font[42]{'+lumi_values[iplot]+' fb^{-1} (13 TeV)}')
    top_pad.Modified()

    can.cd()
    bot_pad.Draw('same')
    bot_pad.cd()
    dummy_hist_lower.Draw()
    for ivar2 in range(len(var2_bin_centers)):
      lower_graphs[ivar2].Draw('same P')
    bot_pad.Modified()
    can.Draw()
    can.SaveAs('plots/zgcorr__'+plot_names[iplot])

