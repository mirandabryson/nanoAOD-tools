import ROOT
import os
import numpy as np
import pandas as pd
import math
import glob
import itertools
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

def hasBit(value,bit):
  """Check if i'th bit is set to 1, i.e. binary of 2^(i-1),
  from the right to the left, starting from position i=0."""
  # https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/mc102X_doc.html#GenPart
  # Gen status flags, stored bitwise, are:
  #    0: isPrompt,                          8: fromHardProcess,
  #    1: isDecayedLeptonHadron,             9: isHardProcessTauDecayProduct,
  #    2: isTauDecayProduct,                10: isDirectHardProcessTauDecayProduct,
  #    3: isPromptTauDecayProduct,          11: fromHardProcessBeforeFSR,
  #    4: isDirectTauDecayProduct,          12: isFirstCopy,
  #    5: isDirectPromptTauDecayProduct,    13: isLastCopy,
  #    6: isDirectHadronDecayProduct,       14: isLastCopyBeforeFSR
  #    7: isHardProcess,
  ###return bin(value)[-bit-1]=='1'
  ###return format(value,'b').zfill(bit+1)[-bit-1]=='1'
  return (value & (1 << bit))>0

class GenAnalyzer(Module):

    def __init__(self):
        pass

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        ## Define a first minimum set of objects needed for the analysis
        #FIXME objects should have cross-reference to full collection

        # New collection of Muons
        self.out.branch("GenL_pt", "F",         lenVar="nGenL")
        self.out.branch("GenL_eta", "F",        lenVar="nGenL")
        self.out.branch("GenL_phi", "F",        lenVar="nGenL")
        self.out.branch("GenL_pdgId", "I",      lenVar="nGenL")
        self.out.branch("GenL_fromTop", "I",    lenVar="nGenL")
        self.out.branch("GenL_fromTau", "I",    lenVar="nGenL")
        self.out.branch("GenL_fromZ", "I",    lenVar="nGenL")
        self.out.branch("GenL_fromW", "I",    lenVar="nGenL")

        self.out.branch("Spectator_pt", "F",         lenVar="nSpectator")
        self.out.branch("Spectator_eta", "F",        lenVar="nSpectator")
        self.out.branch("Spectator_phi", "F",        lenVar="nSpectator")
        self.out.branch("Spectator_pdgId", "I",      lenVar="nSpectator")

        self.out.branch("Scatter_pt", "F",         lenVar="nScatter")
        self.out.branch("Scatter_eta", "F",        lenVar="nScatter")
        self.out.branch("Scatter_phi", "F",        lenVar="nScatter")
        self.out.branch("Scatter_pdgId", "I",      lenVar="nScatter")

        self.out.branch("W_pt", "F",         lenVar="nW")
        self.out.branch("W_eta", "F",        lenVar="nW")
        self.out.branch("W_phi", "F",        lenVar="nW")
        self.out.branch("W_pdgId", "I",    lenVar="nW")
        self.out.branch("W_fromTop", "I",    lenVar="nW")

        self.out.branch("Top_pt", "F",         lenVar="nTop")
        self.out.branch("Top_eta", "F",        lenVar="nTop")
        self.out.branch("Top_phi", "F",        lenVar="nTop")
        self.out.branch("Top_pdgId", "I",        lenVar="nTop")

        # Counter for good b-tags
        self.out.branch("nLepFromTop",     "I")
        self.out.branch("nLepFromTau",    "I")
        self.out.branch("nLepFromW",    "I")
        self.out.branch("nLepFromZ",    "I")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def deltaPhi(self, phi1, phi2):
        dphi = phi2-phi1
        if  dphi > math.pi:
            dphi -= 2.0*math.pi
        if dphi <= -math.pi:
            dphi += 2.0*math.pi
        return abs(dphi)

    def deltaR2(self, l1, l2):
        return self.deltaPhi(l1.phi, l2.phi)**2 + (l1.eta - l2.eta)**2

    def deltaR(self, l1, l2):
        return math.sqrt(self.deltaR2(l1,l2))

    def hasAncestor(self, p, ancestorPdg, genParts):
        motherIdx = p.genPartIdxMother
        while motherIdx>0:
            if (abs(genParts[motherIdx].pdgId) == ancestorPdg): return True
            motherIdx = genParts[motherIdx].genPartIdxMother
        return False
        

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        spectator = []

        GenParts    = Collection(event, "GenPart")

        tops = [ p for p in GenParts if (abs(p.pdgId)==6 and hasBit(p.statusFlags,13) ) ] # last copy ts
        #print len(tops)

        Ws = [ p for p in GenParts if (abs(p.pdgId)==24 and hasBit(p.statusFlags,13) ) ] # last copy Ws

        leptons = [ p for p in GenParts if ((abs(p.pdgId)==11 or abs(p.pdgId)==13) and hasBit(p.statusFlags,13) and (hasBit(p.statusFlags,0) or hasBit(p.statusFlags, 2)) ) ]

        scatter = [ p for p in GenParts if (p.genPartIdxMother==0 or p.genPartIdxMother==1) ] # get the intial particles from 2->N scattering
        
        # spectator: light quark from index 0/1 in gen collection
        spectator = [ p for p in scatter if ( abs(p.pdgId)<6 or abs(p.pdgId)==21 ) ]

        for W in Ws:
            #fromTop = self.hasAncestor(W, 6, GenParts)
            W.fromTop = ( 1 if self.hasAncestor(W, 6, GenParts) else 0 )

        for lep in leptons:
            lep.fromTop = ( 1 if self.hasAncestor(lep, 6, GenParts) else 0 )
            lep.fromTau = ( 1 if self.hasAncestor(lep, 15, GenParts) else 0 )
            lep.fromZ = ( 1 if self.hasAncestor(lep, 23, GenParts) else 0 )
            lep.fromW = ( 1 if self.hasAncestor(lep, 24, GenParts) else 0 )

        #print len(tops), len(Ws), len(scatter), len(spectator), len(leptons)

        ## spectator for our signal results in forward jet :)
        #if len(spectator)>0:
        #    print spectator[0].pt, spectator[0].eta

        ## make pandas dataframe out of list
        #leptons_pd = pd.DataFrame(leptons)

        self.out.fillBranch("nLepFromTop", sum( [ l.fromTop for l in leptons ] ) )
        self.out.fillBranch("nLepFromTau", sum( [ l.fromTau for l in leptons ] ) )
        self.out.fillBranch("nLepFromW",   sum( [ l.fromW for l in leptons ] ) )
        self.out.fillBranch("nLepFromZ",   sum( [ l.fromZ for l in leptons ] ) )

        self.out.fillBranch("nGenL",          len(leptons) )
        if len(leptons)>0:
            self.out.fillBranch("GenL_pt",        [ l.pt for l in leptons ])
            self.out.fillBranch("GenL_eta",       [ l.eta for l in leptons ])
            self.out.fillBranch("GenL_phi",       [ l.phi for l in leptons ])
            self.out.fillBranch("GenL_pdgId",     [ l.pdgId for l in leptons ])
            self.out.fillBranch("GenL_fromTop",   [ l.fromTop for l in leptons ])
            self.out.fillBranch("GenL_fromTau",   [ l.fromTau for l in leptons ])
            self.out.fillBranch("GenL_fromW",     [ l.fromW for l in leptons ])
            self.out.fillBranch("GenL_fromZ",     [ l.fromZ for l in leptons ])


        self.out.fillBranch("nScatter",          len(scatter) )
        if len(scatter)>0:
            self.out.fillBranch("Scatter_pt",        [ p.pt    for p in scatter ])
            self.out.fillBranch("Scatter_eta",       [ p.eta   for p in scatter ])
            self.out.fillBranch("Scatter_phi",       [ p.phi   for p in scatter ])
            self.out.fillBranch("Scatter_pdgId",     [ p.pdgId for p in scatter ])

        self.out.fillBranch("nSpectator",          len(spectator) )
        if len(spectator)>0:
            self.out.fillBranch("Spectator_pt",        [ p.pt    for p in spectator ])
            self.out.fillBranch("Spectator_eta",       [ p.eta   for p in spectator ])
            self.out.fillBranch("Spectator_phi",       [ p.phi   for p in spectator ])
            self.out.fillBranch("Spectator_pdgId",     [ p.pdgId for p in spectator ])

        self.out.fillBranch("nW",          len(Ws) )
        if len(Ws)>0:
            self.out.fillBranch("W_pt",        [ p.pt    for p in Ws ])
            self.out.fillBranch("W_eta",       [ p.eta   for p in Ws ])
            self.out.fillBranch("W_phi",       [ p.phi   for p in Ws ])
            self.out.fillBranch("W_pdgId",     [ p.pdgId for p in Ws ])
            self.out.fillBranch("W_fromTop",   [ p.fromTop for p in Ws ])

        self.out.fillBranch("nTop",          len(tops) )
        if len(tops)>0:
            self.out.fillBranch("Top_pt",        [ p.pt    for p in tops ])
            self.out.fillBranch("Top_eta",       [ p.eta   for p in tops ])
            self.out.fillBranch("Top_phi",       [ p.phi   for p in tops ])
            self.out.fillBranch("Top_pdgId",     [ p.pdgId for p in tops ])

        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

genAnalyzer = lambda : GenAnalyzer( )
