#!/bin/env python
import glob
import ROOT
import os, sys
import subprocess
import json
import argparse

  #log_file.write(+"\n")
def find_nanoaod_file(nanoaod_folder, run, event, log_file):
  print("[Info] Finding ucsb nanoaod file")
  log_file.write("[Info] Finding ucsb nanoaod file"+"\n")
  nanoaod_filenames = glob.glob(nanoaod_folder+"/*"+dataset+"*.root")
  print("Finding "+nanoaod_folder+"/*"+dataset+"*.root")
  found_nanoaod = []
  for nanoaod_filename in nanoaod_filenames:
    nanoaod = ROOT.TFile(nanoaod_filename)
    tree_nanoaod = nanoaod.Get("Events")
    #print(nanoaod_filename)
    #print(tree_nanoaod.GetEntries("run=="+str(run)+"&&event=="+str(event)))
    if (tree_nanoaod.GetEntries("run=="+str(run)+"&&event=="+str(event)) != 0):
      print("  Found: "+nanoaod_filename)
      log_file.write("  Found: "+nanoaod_filename+"\n")
      found_nanoaod.append(nanoaod_filename)
  return found_nanoaod

def find_pico(run, event, dataset, pico_folder, log_file):
  # Find which pico file has event
  print("[Info] Finding ucsb pico file")
  log_file.write("[Info] Finding ucsb pico file"+"\n")
  pico_filenames = glob.glob(pico_folder+"/*"+dataset+"*.root")
  print("Finding run:"+str(run)+" event: "+str(event)+" from "+pico_folder+"/*"+dataset+"*.root")
  log_file.write("Finding run:"+str(run)+" event: "+str(event)+" from "+pico_folder+"/*"+dataset+"*.root"+"\n")
  found_pico = []
  for pico_filename in pico_filenames:
    pico = ROOT.TFile(pico_filename)
    tree_pico = pico.Get("tree")
    #print(pico_filename)
    #print(tree_pico.GetEntries("run=="+str(run)+"&&event=="+str(event)))
    if (tree_pico.GetEntries("run=="+str(run)+"&&event=="+str(event)) != 0):
      print("  "+pico_filename)
      log_file.write("  "+pico_filename+"\n")
      found_pico.append(pico_filename)
  return found_pico

def get_miniaod(run, event, dataset, nanoaod_folder, pico_folder, pico_prefix, output_folder, log_file):
  # Find which pico file has event
  #print("[Info] Finding ucsb pico file")
  #pico_filenames = glob.glob(pico_folder+"/*"+dataset+"*.root")
  #found_pico = []
  #for pico_filename in pico_filenames:
  #  pico = ROOT.TFile(pico_filename)
  #  tree_pico = pico.Get("tree")
  #  #print(pico_filename)
  #  #print(tree_pico.GetEntries("run=="+str(run)+"&&event=="+str(event)))
  #  if (tree_pico.GetEntries("run=="+str(run)+"&&event=="+str(event)) != 0):
  #    print("  "+pico_filename)
  #    found_pico.append(pico_filename)

  if pico_folder!="none":
    found_pico = find_pico(run,event,dataset, pico_folder, log_file)
    #found_pico = ["/cms29r0/pico/picov5/nano/2016/mc/TTJets_DiLept_genMET-150_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISummer16picov5__PUMoriond17_Nano1June2019_102X_mcRun2_asymptotic_v7-v1__40000__2956186E-4503-6A43-8C69-8F1F10C19A06.root"]
    if len(found_pico) == 0:
      print("[Error] Could not find event")
      log_file.write("[Error] Could not find event\n")
      return False
    if len(found_pico) != 1:
      print("[Error] Found multiple files with event")
      log_file.write("[Error] Found multiple files with event\n")
      return False
    # Find nano
    found_nanoaod = nanoaod_folder + found_pico[0].replace(pico_prefix,"").replace(pico_folder,"")
    if not os.path.exists(found_nanoaod):
      print('[Error] Could not find nanoaod: '+found_nanoaod)
      log_file.write('[Error] Could not find nanoaod: '+found_nanoaod+'\n')
      return False
  else:
    # Find nano
    found_nanoaods = find_nanoaod_file(nanoaod_folder, run, event, log_file)
    if len(found_nanoaods) == 0:
      print("[Error] Could not find event")
      log_file.write("[Error] Could not find event\n")
      return False
    if len(found_nanoaods) != 1:
      print("[Error] Found multiple files with event")
      log_file.write("[Error] Found multiple files with event\n")
      return False
    found_nanoaod = found_nanoaods[0]
  
  # Find lumi for event
  print("[Info] Find lumiblock")
  log_file.write("[Info] Find lumiblock"+'\n')
  nanoaod = ROOT.TFile(found_nanoaod)
  tree_nanoaod = nanoaod.Get("Events")
  tree_nanoaod.Draw("run:luminosityBlock:event","event=="+str(event),"goff")
  lumi = int(tree_nanoaod.GetV2()[0])
  run = int(tree_nanoaod.GetV1()[0])
  print("  "+str(lumi))

  # Find parent miniAOD files
  print("[Info] Finding parent miniAOD files")
  log_file.write("[Info] Finding parent miniAOD files"+'\n')
  nanoaod_ucsbname = found_nanoaod
  nanoaod_ucsbname_split = nanoaod_ucsbname.split("/")[-1].split("__")
  #modified for data
  nanoaod_dasname = "/store/data/"+nanoaod_ucsbname_split[1]+"/"+nanoaod_ucsbname_split[0]+"/NANOAOD/"+nanoaod_ucsbname_split[2]+"/"+nanoaod_ucsbname_split[3]+"/"+nanoaod_ucsbname_split[4]
  if not ('Run' in nanoaod_ucsbname_split[1]):
    nanoaod_dasname = "/store/mc/"+nanoaod_ucsbname_split[1]+"/"+nanoaod_ucsbname_split[0]+"/NANOAODSIM/"+nanoaod_ucsbname_split[2]+"/"+nanoaod_ucsbname_split[3]+"/"+nanoaod_ucsbname_split[4]
  #print(nanoaod_dasname)
  command = 'dasgoclient -query="parent file='+nanoaod_dasname+'"'
  print('Running: '+command)
  miniaod_dasnames = subprocess.check_output(command, shell=True).split()
  print("[Info] Finding miniAOD with lumi")
  log_file.write("[Info] Finding miniAOD with lumi"+'\n')
  # Find miniAOD with lumi
  found_miniaods = []
  for miniaod_dasname in miniaod_dasnames:
    #don't know if this can be done in MC, but it must be done in data to avoid wrong runs(?): skip if wrong run
    #run_mini_string = str(run)[0:3]+'/'+str(run)[3:6]
    #if not (run_mini_string in miniaod_dasname):
    #  continue
    command = 'dasgoclient -query="lumi file='+miniaod_dasname+'" -json'
    print('Running: '+command)
    result = subprocess.check_output(command, shell=True)
    result_dict = json.loads(result)
    #store listo of tuples [(lumiblock,run),(lumiblock,run),...]
    lumis = [(result_dict[i]['lumi'][0]['number'],result_dict[i]['lumi'][0]['run_number']) for i in range(len(result_dict))]
    print(lumis)
    if (lumi,run) in lumis:
      print("  "+miniaod_dasname)
      log_file.write("  "+miniaod_dasname+'\n')
      found_miniaods.append(miniaod_dasname)

  if len(found_miniaods) == 0:
    print("[Error] Could not find event")
    log_file.write("[Error] Could not find event"+'\n')
    return False
  if len(found_miniaods) != 1:
    print("[Error] Found multiple miniAODs with same lumi")
    print(found_miniaods)
    log_file.write("[Error] Found multiple miniAODs with same lumi"+'\n')
    log_file.write(found_miniaods+'\n')
    return False

  print('')

  #get miniAOD -mo
  command = 'xrdcp root://cms-xrd-global.cern.ch//'+found_miniaods[0]+' '+output_folder+'/'
  print('running '+command)
  os.system(command)
  
  #don't know what any of this is, but commenting it out -mo
  ##found_miniaods = ["/store/mc/RunIISummer16MiniAODv3/TTJets_DiLept_genMET-150_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/80000/7C362B40-6BF0-E811-938D-AC1F6B0DE228.root"]
  #output_filename = dataset.replace("*","")+"_run"+str(run)+"_lumi"+str(lumi)+"_event"+str(event)+".root"
  #print("[Info] Running below commands to get event")
  #log_file.write("[Info] Running below commands to get event"+'\n')
  #edm_command = "edmCopyPickMerge inputFiles="+found_miniaods[0]+" eventsToProcess="+str(run)+":"+str(event)+" outputFile="+output_filename
  #print('')
  ##print('  ssh uaf-8;voms-proxy-init -voms cms -valid 168:0')
  #print('  ssh uaf-8 "source setCMSEnv &&'+edm_command+'" && scp uaf-8:'+output_filename+" "+output_filename+" && ssh uaf-8 rm "+output_filename)
  #log_file.write('  ssh uaf-8 "source setCMSEnv &&'+edm_command+'" && scp uaf-8:'+output_filename+" "+output_filename+" && ssh uaf-8 rm "+output_filename+'\n')
  #os.system('ssh uaf-8 "source setCMSEnv &&'+edm_command+'" && scp uaf-8:'+output_filename+" "+output_folder+" && ssh uaf-8 rm "+output_filename)
  #print('')
  #print("[Info] Running below commands to print decay tree")
  #log_file.write("[Info] Running below commands to print decay tree"+'\n')
  #print('  cmsRun scripts/cmssw_treedraw.py inputFiles=file:'+output_folder+"/"+output_filename+'>'+output_folder+'/decay_tree_'+output_filename.replace('.root','')+'.txt')
  #log_file.write('  cmsRun scripts/cmssw_treedraw.py inputFiles=file:'+output_folder+"/"+output_filename+'>'+output_folder+'/decay_tree_'+output_filename.replace('.root','')+'.txt'+'\n')
  #os.system('cmsRun scripts/cmssw_treedraw.py inputFiles=file:'+output_folder+"/"+output_filename+'>'+output_folder+'/decay_tree_'+output_filename.replace('.root','')+'.txt')
  return True


if __name__ == "__main__":
  # Example: ./scripts/jb_scripts/get_miniAOD.py -e "1:160306" -p none -o trash "SingleLeptFromT_TuneCUETP8M1*v1__100*"
  # Example: ./scripts/get_miniAOD.py -el eventlist_ra2b -o ra2b_events "TTJets_SingleLeptFromT_Tune"

  parser = argparse.ArgumentParser(description='''\
Script that downloads MiniAOD files that contain a specific run and event number.
Also uses scripts/cmssw_treedraw.py to write the decay chain in a text file.
Uses uaf-8 to download the file.
Use either -r and -e or -el.
[Example]
./get_miniAOD.py -el eventlist_ra2b -o ra2b_events "TTJets_SingleLeptFromT_Tune
''', formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('dataset', help='Datset name the event is in.')
  parser.add_argument('-e', '--run_event', help='Run and event number. Should written like "1:2,1:4,1:10", where first number is run')
  parser.add_argument('-el', '--eventlist', help='An alternative method to specify run and event number using a file. The file content should be like below.\n run event\n 1 2\n 1 4\n1 10')
  parser.add_argument('-n', '--nanoaod_folder', default='/net/cms17/cms17r0/pico/NanoAODv9/nano/2018/data/', help='Folder that has nanoaod files. Used to identify the lumiblock number')
  parser.add_argument('-p', '--pico_folder', default='/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath_v3/2018/data/skim_met150/', help="Use 'none' if searching only with nanoaod")
  parser.add_argument('-pp', '--pico_prefix', default='pico_met150_', help='The prefix that is applied to the pico files. Used to find the dataset name.')
  parser.add_argument('-o', '--output_folder', default='./miniaod', help='Folder that will containt the output miniAOD.')
  args = parser.parse_args()
  
  log_filename = 'get_miniAOD.log'
  log_file = open(args.output_folder+"/"+log_filename, 'w')
  #run = 1
  #event = 421491
  #dataset = "TTJets_DiLept_genMET-150"
  #nanoaod_folder = "/cms29r0/pico/NanoAODv5/nano/2016/mc"
  dataset = args.dataset
  nanoaod_folder = args.nanoaod_folder
  pico_folder = args.pico_folder
  pico_prefix = args.pico_prefix
  output_folder = args.output_folder

  run_events = []
  if args.eventlist:
    os.system("cp "+args.eventlist+" "+args.output_folder+"/"+os.path.basename(args.eventlist))
    with open(args.eventlist) as file_eventlist:
      line_split = file_eventlist.readline().split()
      run_index = line_split.index('run')
      event_index = line_split.index('event')
      for line in file_eventlist:
        line_split = line.split()
        run = line_split[run_index]
        event = line_split[event_index]
        if run == 'run': continue
        run_events.append(run+":"+event)
  else:
    #run_events = json.loads(args.run_event)
    for line in args.run_event.split(","):
      line = line.rstrip()
      run_events.append(line)


  print("Run below commands first")
  print("  voms-proxy-init -voms cms -valid 168:0")
  print("  ssh uaf-8")
  print("  voms-proxy-init -voms cms -valid 168:0")
  raw_input("Press Enter")
  for run_event in run_events:
    run = run_event.split(':')[0]
    event = run_event.split(':')[1]
    get_miniaod(run, event, dataset, nanoaod_folder, pico_folder, pico_prefix, output_folder, log_file)
