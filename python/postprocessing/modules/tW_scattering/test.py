from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor   import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel       import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop       import Module
from PhysicsTools.NanoAODTools.postprocessing.modules.tW_scattering.ObjectSelection import *
from PhysicsTools.NanoAODTools.postprocessing.modules.tW_scattering.GenAnalyzer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.tW_scattering.lumiWeightProducer import *

#json support to be added
modules = [\
    lumiWeightProd(1),
    genAnalyzer(),
    selector2018(),
    ]
# apply PV requirement
cut  = 'PV_ndof>4 && sqrt(PV_x*PV_x+PV_y*PV_y)<=2 && abs(PV_z)<=24'

# loose skim
cut += '&& (nElectron+nMuon)>0 && Sum$(Jet_pt>25&&abs(Jet_eta)<2.4)>=4'

p = PostProcessor('./', ['/hadoop/cms/store/user/dspitzba/tW_scattering/tW_scattering/nanoAOD/tW_scattering_nanoAOD_500.root'], cut=cut, modules=modules,\
    branchsel='keep_and_drop_in.txt',\
    outputbranchsel='keep_and_drop.txt' )

p.run()
