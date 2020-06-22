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
from PhysicsTools.NanoAODTools.postprocessing.framework.mt2Calculator import mt2Calculator


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

        self.out.branch("H_pt", "F",         lenVar="nH")
        self.out.branch("H_eta", "F",        lenVar="nH")
        self.out.branch("H_phi", "F",        lenVar="nH")
        self.out.branch("H_pdgId", "I",    lenVar="nH")

        self.out.branch("b_pt", "F",         lenVar="nb")
        self.out.branch("b_eta", "F",        lenVar="nb")
        self.out.branch("b_phi", "F",        lenVar="nb")
        self.out.branch("b_pdgId", "I",    lenVar="nb")
        self.out.branch("b_fromH", "I",    lenVar="nb")

        self.out.branch("j_pt", "F",         lenVar="nj")
        self.out.branch("j_eta", "F",        lenVar="nj")
        self.out.branch("j_phi", "F",        lenVar="nj")
        self.out.branch("j_pdgId", "I",    lenVar="nj")
        self.out.branch("j_fromW", "I",    lenVar="nj")

        self.out.branch("Top_pt", "F",         lenVar="nTop")
        self.out.branch("Top_eta", "F",        lenVar="nTop")
        self.out.branch("Top_phi", "F",        lenVar="nTop")
        self.out.branch("Top_pdgId", "I",        lenVar="nTop")

        # Counter for good b-tags
        self.out.branch("nLepFromTop",      "I")
        self.out.branch("nLepFromTau",      "I")
        self.out.branch("nLepFromW",        "I")
        self.out.branch("nLepFromZ",        "I")

        #miranda time
        self.out.branch("DeltaRbb", "F",   lenVar="nH")
        self.out.branch("MCT", "F",        lenVar="nH")
        self.out.branch("MT2_WH", "F",     lenVar="nH")
        self.out.branch("MT2_bbjj", "F",     lenVar="nH")
        self.out.branch("MT2_bjj_b", "F",     lenVar="nH")
        self.out.branch("Mbb", "F",         lenVar="nH")


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
    
    def MCT2(self, b1, b2):
      return 2*b1.pt*b2.pt*(1+math.cos(self.deltaPhi(b1.phi,b2.phi)))

    def MCT(self, b1, b2):
      return math.sqrt(self.MCT2(b1,b2))

    def Mbb(self, b1, b2):
      bjet1 = ROOT.TLorentzVector()
      bjet1.SetPtEtaPhiM(b1.pt, b1.eta, b1.phi, 0)
      bjet2 = ROOT.TLorentzVector()
      bjet2.SetPtEtaPhiM(b2.pt, b2.eta, b2.phi, 0)
      return (bjet1 + bjet2).M()

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

#        print 'event %d *dr phil voice* shut the hell up bitch go kill yourself'%( event.event)

        #MET
        met_pt  = event.MET_pt
        met_phi = event.MET_phi
        
        #Particles
        spectator = []

        GenParts    = Collection(event, "GenPart")

        tops = [ p for p in GenParts if (abs(p.pdgId)==6 and hasBit(p.statusFlags,13) ) ] # last copy ts
        #print len(tops)

        Ws = [ p for p in GenParts if (abs(p.pdgId)==24 and hasBit(p.statusFlags,13) ) ] # last copy Ws

        js = [ p for p in GenParts if ((abs(p.pdgId) == 1 or abs(p.pdgId) == 2 or abs(p.pdgId) == 3 or abs(p.pdgId) == 4 or abs(p.pdgId) == 5 or abs(p.pdgId) == 6) and hasBit(p.statusFlags, 7))] 

        Hs = [ p for p in GenParts if (abs(p.pdgId)==25 and hasBit(p.statusFlags,13) )] # last copy Hs

        bs = [ p for p in GenParts if (abs(p.pdgId) == 5 and hasBit(p.statusFlags, 7))]  #hard scatter bs

        leptons = [ p for p in GenParts if ((abs(p.pdgId)==11 or abs(p.pdgId)==13) and hasBit(p.statusFlags,13) and (hasBit(p.statusFlags,0) or hasBit(p.statusFlags, 2)) ) ]

        scatter = [ p for p in GenParts if (p.genPartIdxMother==0 or p.genPartIdxMother==1) ] # get the intial particles from 2->N scattering
        
        # spectator: light quark from index 0/1 in gen collection
        spectator = [ p for p in scatter if ( abs(p.pdgId)<6 or abs(p.pdgId)==21 ) ]

        for W in Ws:
            #fromTop = self.hasAncestor(W, 6, GenParts)
          W.fromTop = ( 1 if self.hasAncestor(W, 6, GenParts) else 0 )

        for j in js:
          j.fromW = (1 if self.hasAncestor(j, 24, GenParts) else 0 )

        for b in bs:
          b.fromH = ( 1 if self.hasAncestor(b, 25, GenParts) else 0 )



        deltaRbb = []
        Mct = []
        mbb = []

        for b,s in itertools.combinations(bs, 2):
          if (b.fromH == 1 and s.fromH == 1):
            dr = self.deltaR(b, s)
            mct = self.MCT(b,s)
            Mbb = self.Mbb(b,s)
          else:
            continue
          deltaRbb.append(dr)
          Mct.append(mct)
          mbb.append(Mbb)

        for lep in leptons:
            lep.fromTop = ( 1 if self.hasAncestor(lep, 6, GenParts) else 0 )
            lep.fromTau = ( 1 if self.hasAncestor(lep, 15, GenParts) else 0 )
            lep.fromZ = ( 1 if self.hasAncestor(lep, 23, GenParts) else 0 )
            lep.fromW = ( 1 if self.hasAncestor(lep, 24, GenParts) else 0 )



        #MT2WH
            
        mt2WH = []

        mt2Calculator.setMet(met_pt, met_phi)
        
        for H in Hs:
          mt2Calculator.setJet1(H.pt, H.eta, H.phi)
        for W in Ws:
          mt2Calculator.setJet2(W.pt, W.eta, W.phi)
        mt2WH.append(mt2Calculator.mt2jj())
        
        #MT2bbjj

        mt2bbjj =[]
        mt2bjjb = []

        for b,s in itertools.combinations(bs, 2):
          if (b.fromH == 1 and s.fromH == 1):
            mt2Calculator.setBJet1(b.pt, b.eta, b.phi)
            mt2Calculator.setBJet2(s.pt, s.eta, s.phi)
        for j,s in itertools.combinations(js, 2):
          if (j.fromW == 1 and s.fromW == 1):
            mt2Calculator.setJet1(j.pt, j.eta, j.phi)
            mt2Calculator.setJet2(s.pt, s.eta, s.phi)
        mt2bbjj.append(mt2Calculator.mt2bbjj())
        mt2bjjb.append(mt2Calculator.mt2bjjb())

        

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

        self.out.fillBranch("DeltaRbb",     deltaRbb)
        self.out.fillBranch("MCT",          Mct)
        self.out.fillBranch("MT2_WH",       mt2WH)
        self.out.fillBranch("MT2_bbjj",       mt2bbjj)
        self.out.fillBranch("MT2_bjj_b",       mt2bjjb)
        self.out.fillBranch("Mbb",           mbb)

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

        self.out.fillBranch("nH",          len(Hs) )
        if len(Hs)>0:
            self.out.fillBranch("H_pt",        [ p.pt    for p in Hs ])
            self.out.fillBranch("H_eta",       [ p.eta   for p in Hs ])
            self.out.fillBranch("H_phi",       [ p.phi   for p in Hs ])
            self.out.fillBranch("H_pdgId",     [ p.pdgId for p in Hs ])

        self.out.fillBranch("nb",          len(bs) )
        if len(bs)>0:
            self.out.fillBranch("b_pt",        [ p.pt    for p in bs ])
            self.out.fillBranch("b_eta",       [ p.eta   for p in bs ])
            self.out.fillBranch("b_phi",       [ p.phi   for p in bs ])
            self.out.fillBranch("b_pdgId",     [ p.pdgId for p in bs ])
            self.out.fillBranch("b_fromH",     [ p.fromH for p in bs ])

        self.out.fillBranch("nj",          len(js) )
        if len(js)>0:
            self.out.fillBranch("j_pt",        [ p.pt    for p in js ])
            self.out.fillBranch("j_eta",       [ p.eta   for p in js ])
            self.out.fillBranch("j_phi",       [ p.phi   for p in js ])
            self.out.fillBranch("j_pdgId",     [ p.pdgId for p in js ])
            self.out.fillBranch("j_fromW",     [ p.fromW for p in js ])

        self.out.fillBranch("nTop",          len(tops) )
        if len(tops)>0:
            self.out.fillBranch("Top_pt",        [ p.pt    for p in tops ])
            self.out.fillBranch("Top_eta",       [ p.eta   for p in tops ])
            self.out.fillBranch("Top_phi",       [ p.phi   for p in tops ])
            self.out.fillBranch("Top_pdgId",     [ p.pdgId for p in tops ])


        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

genAnalyzer = lambda : GenAnalyzer( )
