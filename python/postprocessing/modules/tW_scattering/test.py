from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor   import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel       import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop       import Module
from PhysicsTools.NanoAODTools.postprocessing.modules.tW_scattering.ObjectSelection import *
from PhysicsTools.NanoAODTools.postprocessing.modules.tW_scattering.GenAnalyzer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.tW_scattering.lumiWeightProducer import *

from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2       import *

jetmet = createJMECorrector(isMC=True, dataYear=2018, jesUncert="Total", jetType = "AK4PFchs", applySmearing = True, isFastSim = False )

#json support to be added
modules = [\
    lumiWeightProd(1),
    jetmet(),
    genAnalyzer(),
    selector2018(),
    ]
# apply PV requirement
cut  = 'PV_ndof>4 && sqrt(PV_x*PV_x+PV_y*PV_y)<=2 && abs(PV_z)<=24'

# loose skim
cut += '&& (Sum$(Electron_pt>30&&abs(Electron_eta)<2.4&&Electron_miniPFRelIso_all<0.1&&Electron_cutBased>=3)+Sum$(Muon_pt>25&&abs(Muon_eta)<2.4&&Muon_mediumId>0&&Muon_miniPFRelIso_all<0.1))>0'
cut += '&& ( (Sum$(Jet_pt>25&&abs(Jet_eta)<2.4)>=4) || (Sum$(Jet_pt>25&&abs(Jet_eta)<2.4)>=2 && (Sum$(Electron_pt>10&&abs(Electron_eta)<2.4)+Sum$(Muon_pt>10&&abs(Muon_eta)<2.4&&Muon_mediumId>0))>=3) )'

p = PostProcessor('./', ['/hadoop/cms/store/user/dspitzba/tW_scattering/tW_scattering/nanoAOD/tW_scattering_nanoAOD_500.root'], cut=cut, modules=modules,\
#    branchsel='keep_and_drop_in.txt',\
    outputbranchsel='keep_and_drop.txt' )

p.run()
