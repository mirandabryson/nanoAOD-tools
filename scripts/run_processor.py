#!/usr/bin/env python

import sys
import re

from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor   import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel       import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop       import Module

from PhysicsTools.NanoAODTools.postprocessing.modules.tW_scattering.ObjectSelection import *
from PhysicsTools.NanoAODTools.postprocessing.modules.tW_scattering.GenAnalyzer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.tW_scattering.lumiWeightProducer import *

from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2       import *


files       = sys.argv[1].split(',')
lumiWeight  = float(sys.argv[2])
isData      = int(sys.argv[3]) == 1
year        = int(sys.argv[4])
era         = sys.argv[5]
isFastSim   = int(sys.argv[6]) == 1

if files[0].startswith('/store/'):
    print "Had to add a prefix"
    files = [ 'root://cmsxrootd.fnal.gov/' + f for f in files ]

#json support to be added

print "Sumweight:", lumiWeight
print "Data:", isData
print "Year:", year
print "FastSim", isFastSim 
print "Files:", files

jetmet = createJMECorrector(isMC=(not isData), dataYear=year, runPeriod=era, jesUncert="Total", jetType = "AK4PFchs", applySmearing = True, isFastSim = isFastSim )

modules = [\
    lumiWeightProd(lumiWeight, isData),
    jetmet()
    ]

if not isData:
    isW = False
    isWExt = False
    for f in files:
        print f
        if re.search("W[1-4]Jets", f): isW=True
        if re.search("NuPt", f): isWExt=True
        print isW, isWExt
    modules += [genAnalyzer(isW, isWExt)]

modules += [\
    selector(year, isData),
    ]

if isData:
    if year==2018:
        jsonInput='PhysicsTools/NanoAODTools/python/postprocessing/modules/tW_scattering/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'
    else:
        jsonInput=None # FIXME
else:
    jsonInput=None

# apply PV requirement	
cut  = 'PV_ndof>4 && sqrt(PV_x*PV_x+PV_y*PV_y)<=2 && abs(PV_z)<=24'
# loose skim	
#cut += '&& MET_pt>200'
cut += '&& Sum$(Jet_pt>30&&abs(Jet_eta<2.4))>=2'

p = PostProcessor('./', files, cut=cut, modules=modules,fwkJobReport=True, prefetch=True,\
#    branchsel='PhysicsTools/NanoAODTools/python/postprocessing/modules/tW_scattering/keep_and_drop_in.txt',\
    outputbranchsel='PhysicsTools/NanoAODTools/python/postprocessing/modules/tW_scattering/keep_and_drop.txt',
    jsonInput=jsonInput )

p.run()
