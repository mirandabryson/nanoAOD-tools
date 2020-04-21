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



class PhysicsObjects(Module):

    def __init__(self, year=2018):
        self.year = year
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
        self.out.branch("Lepton_pt", "F", lenVar="nLepton")
        self.out.branch("Lepton_eta", "F", lenVar="nLepton")
        self.out.branch("Lepton_phi", "F", lenVar="nLepton")
        self.out.branch("Lepton_pdgId", "I", lenVar="nLepton")
        self.out.branch("Lepton_miniIso", "F", lenVar="nLepton")
        self.out.branch("Lepton_muIndex", "I", lenVar="nLepton")
        self.out.branch("Lepton_elIndex", "I", lenVar="nLepton")

        ## New collection of Jets. #FIXME overlap removal with GoodLeptons
        #self.out.branch("GoodJet_pt", "F", lenVar="nGoodJet")
        #self.out.branch("GoodJet_eta", "F", lenVar="nGoodJet")
        #self.out.branch("GoodJet_phi", "F", lenVar="nGoodJet")
        #self.out.branch("GoodJet_btag", "F", lenVar="nGoodJet")

        ##
        self.out.branch("Muon_isVeto",      "F", lenVar="nMuon")
        self.out.branch("Muon_isTight",     "F", lenVar="nMuon")
        self.out.branch("Electron_isVeto",  "F", lenVar="nElectron")
        self.out.branch("Electron_isTight", "F", lenVar="nElectron")
        self.out.branch("Jet_isGoodJet",    "F", lenVar="nJet")
        self.out.branch("Jet_isGoodBJet",   "F", lenVar="nJet")
        self.out.branch("Jet_crossClean",   "F", lenVar="nJet")


        # Counter for good b-tags
        self.out.branch("nLepton",      "I")
        self.out.branch("nGoodJet",     "I")
        self.out.branch("nGoodBTag",    "I")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def isGoodJet(self, jet):
        return (jet.pt > 30 and abs(jet.eta)<2.4 and jet.jetId>1)

    def isGoodBJet(self, jet):
        if self.year == 2018:
            threshold = 0.4184
        return (self.isGoodJet(jet) and jet.btagDeepB > threshold)

    def isVetoMuon(self, muon):
        return (muon.looseId and muon.pt>5 and abs(muon.eta)<2.4 and muon.miniPFRelIso_all < 0.2 and abs(muon.dxy)<0.1 and abs(muon.dz)<0.5)

    def isVetoElectron(self, electron):
        return (electron.cutBased>0 and electron.miniPFRelIso_all < 0.2)

    def isTightMuon(self, muon):
        return (muon.pt > 25 and muon.mediumId and abs(muon.eta)<2.4 and muon.miniPFRelIso_all < 0.1)

    def isTightElectron(self, electron):
        return (electron.pt > 30 and electron.cutBased >= 3 and abs(electron.eta) < 2.4 and electron.miniPFRelIso_all < 0.1)# and electron.sip3d < 4.0 and abs(electron.dxy) < 0.05 and abs(electron.dz) < 0.1)

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

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        muons       = Collection(event, "Muon")
        electrons   = Collection(event, "Electron")
        jets        = Collection(event, "Jet")

        # tight lepton collection, will be sorted by pt
        leptons     = []

        isTightMuon = []
        isVetoMuon = []
        for i,mu in enumerate(muons):
            mu.isTight  = self.isTightMuon(mu)
            mu.isVeto   = self.isVetoMuon(mu)
            isTightMuon.append(mu.isTight)
            isVetoMuon.append(mu.isVeto)
            if self.isTightMuon(mu):
                leptons.append({'pt':mu.pt, 'eta':mu.eta, 'phi':mu.phi, 'pdgId':mu.pdgId, 'miniIso':mu.miniPFRelIso_all, 'muIndex':i, 'elIndex':-1})


        isTightElectron = []
        isVetoElectron = []
        for i,el in enumerate(electrons):
            el.isTight  = self.isTightElectron(el)
            el.isVeto   = self.isVetoElectron(el)
            isTightElectron.append(el.isTight)
            isVetoElectron.append(el.isVeto)
            if self.isTightElectron(el):
                leptons.append({'pt':el.pt, 'eta':el.eta, 'phi':el.phi, 'pdgId':el.pdgId, 'miniIso':el.miniPFRelIso_all, 'muIndex':-1, 'elIndex':i})

        

        cleanMaskV  = []
        isGoodJet   = []
        isGoodBJet  = []

        for j in jets:

            j.cleanMask = 1
            for coll in [electrons, muons]:
                for p in coll:
                    if p.isVeto:
                        if self.deltaR(j, p) < 0.4:
                            j.cleanMask = 0

            isGoodJet.append(1 if (self.isGoodJet(j) and j.cleanMask) else 0)
            isGoodBJet.append(1 if (self.isGoodBJet(j) and j.cleanMask) else 0)
            
            cleanMaskV.append(j.cleanMask)

        self.out.fillBranch("Muon_isTight",      isTightMuon)
        self.out.fillBranch("Muon_isVeto",      isVetoMuon)
        self.out.fillBranch("Electron_isTight",  isTightElectron)
        self.out.fillBranch("Electron_isVeto",  isVetoElectron)
        self.out.fillBranch("Jet_crossClean",   cleanMaskV)
        self.out.fillBranch("Jet_isGoodJet",    isGoodJet)
        self.out.fillBranch("Jet_isGoodBJet",   isGoodBJet)

        # make pandas dataframe out of list
        leptons_pd = pd.DataFrame(leptons)

        self.out.fillBranch("nLepton",          len(leptons_pd) )
        if len(leptons_pd)>0:
            self.out.fillBranch("Lepton_pt",        leptons_pd.sort_values(by='pt', ascending=False)['pt'].tolist() )
            self.out.fillBranch("Lepton_eta",       leptons_pd.sort_values(by='pt', ascending=False)['eta'].tolist() )
            self.out.fillBranch("Lepton_phi",       leptons_pd.sort_values(by='pt', ascending=False)['phi'].tolist() )
            self.out.fillBranch("Lepton_pdgId",     leptons_pd.sort_values(by='pt', ascending=False)['pdgId'].tolist() )
            self.out.fillBranch("Lepton_miniIso",   leptons_pd.sort_values(by='pt', ascending=False)['miniIso'].tolist() )
            self.out.fillBranch("Lepton_muIndex",   leptons_pd.sort_values(by='pt', ascending=False)['muIndex'].tolist() )
            self.out.fillBranch("Lepton_elIndex",   leptons_pd.sort_values(by='pt', ascending=False)['elIndex'].tolist() )

        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

selector2018 = lambda : PhysicsObjects( year=2018 )
