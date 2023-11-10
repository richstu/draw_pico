#!/bin/env python3
import os
import subprocess
import sys

# envDict: { key: value }
def findEnviornment(scriptname, envDict):
  if not os.path.isfile(scriptname):
    print ("[Error] Can't find script:"+scriptname)

  command = ['env', '-i', 'bash', '-c', 'source '+scriptname+' && env']
  proc = subprocess.Popen(command, stdout = subprocess.PIPE, shell=True)
  for line in proc.stdout:
    if '{' in line.decode('utf-8'): continue
    if '}' in line.decode('utf-8'): continue
    t_array= line.decode('utf-8').split("=",1)
    key=t_array[0]
    value=t_array[1].rstrip('\n')
    envDict[key] = value
  proc.communicate()

def returnEnviornment(scriptname):
  envDict = {}
  findEnviornment(scriptname, envDict)
  return envDict

def addRootEnv(_env):
  _env.Append (CCFLAGS = '-isystem `root-config --incdir`' )
  #_env.Append (CCFLAGS = '-g' ) #debug symbols cause big executables, use only when gdb is needed
  _env.Append (CCFLAGS = '`root-config --cflags`' )
  _env.Append (LINKFLAGS = '`root-config --glibs`') 
  _env.Append (LINKFLAGS = '`root-config --ldflags`')
  _env.Append (LINKFLAGS = ['-lRooFit', '-lGenVector', '-lRooStats', '-lRooFitCore', '-lMathMore', '-lTMVA'])
  #_env.Append (LINKFLAGS = ['-lRooFit', '-lGenVector', '-lRooStats', '-lRooFitCore', '-lTMVA']) #MathMore not in ROOT 6.18??

def addWarningEnv(_env):
  _env.Append (CCFLAGS = ['-pedantic', 
                          '-Wall', '-Wextra', '-Wshadow', '-Woverloaded-virtual', '-Wold-style-cast', 
                          '-Wcast-align', '-Wcast-qual', '-Wdisabled-optimization', 
                          '-Wformat=2', '-Wformat-nonliteral', '-Wformat-security', 
                          '-Wformat-y2k', '-Winit-self', '-Winvalid-pch', 
                          '-Wmissing-format-attribute', '-Wmissing-include-dirs', '-Wmissing-noreturn', 
                          '-Wpacked', '-Wpointer-arith', '-Wredundant-decls', '-Wstack-protector', 
                          '-Wswitch-default', '-Wundef', '-Wvariadic-macros', 
                          '-Wwrite-strings', '-Wctor-dtor-privacy', '-Wnon-virtual-dtor', '-Wsign-promo', '-Wsign-compare', 
                          '-Wunused', '-Werror', '-Wlong-long','-Wswitch-enum', '-Wunreachable-code', 
                          #'-Wunsafe-loop-optimizations', '-Wfloat-equal', '-Wsign-conversion', # Makes too many errors
                         ])

def addExternalEnv(_env):
  _env.Append (CCFLAGS = '-isystem external_inc' )

def addBasicEnv(_env):
  _env.Append (CCFLAGS = '-O2')

def addKernelEnv(_env):
  _env['kernel'] = getKernel()

def getKernel():
  return subprocess.check_output("uname -r | cut -d '-' -f1", shell=True, universal_newlines=True).rstrip()

def addCombineEnv(_env):
  kernel = getKernel()
  if (kernel=='2.6.32'): #slc6
    _env.Append (CCFLAGS = '-I/net/cms29/cms29r0/pico/cc7/CMSSW_10_2_11_patch1/src/HiggsAnalysis/CombinedLimit/interface' )
    _env.Append (LINKFLAGS = '-L/net/cms29/cms29r0/pico/CMSSW_10_2_11_patch1/lib/slc6_amd64_gcc700 -lHiggsAnalysisCombinedLimit') 
  else: #assume cc7
    _env.Append (CCFLAGS = '-I/net/cms29/cms29r0/pico/CMSSW_10_2_11_patch1/src/HiggsAnalysis/CombinedLimit/interface' )
    _env.Append (LINKFLAGS = '-L/net/cms29/cms29r0/pico/cc7/CMSSW_10_2_11_patch1/lib/slc7_amd64_gcc700 -lHiggsAnalysisCombinedLimit') 

SConsignFile('kernel/'+getKernel()+'/sconsign.dblite')

# [Requried] Export SET_ENV_PATH to script that sets root and scons environment
analysisEnv = Environment(ENV = returnEnviornment(os.environ['SET_ENV_PATH']))

addBasicEnv(analysisEnv)
addKernelEnv(analysisEnv)
addRootEnv(analysisEnv)
addExternalEnv(analysisEnv)
addWarningEnv(analysisEnv)
addCombineEnv(analysisEnv)

exportEnv = analysisEnv
SConscript('SConscript', variant_dir='build/'+analysisEnv['kernel'], duplicate=0, exports="exportEnv")
