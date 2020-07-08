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
        

        self.out.branch("Lepton_pt", "F", lenVar="nLepton")
        self.out.branch("Lepton_eta", "F", lenVar="nLepton")
        self.out.branch("Lepton_phi", "F", lenVar="nLepton")
        self.out.branch("Lepton_pdgId", "I", lenVar="nLepton")
        self.out.branch("Lepton_miniIso", "F", lenVar="nLepton")
        self.out.branch("Lepton_muIndex", "I", lenVar="nLepton")
        self.out.branch("Lepton_elIndex", "I", lenVar="nLepton")

        self.out.branch("Muon_isVeto",      "I", lenVar="nMuon")
        self.out.branch("Muon_isTight",     "I", lenVar="nMuon")
        self.out.branch("Electron_isVeto",  "I", lenVar="nElectron")
        self.out.branch("Electron_isTight", "I", lenVar="nElectron")
        
        self.out.branch("Jet_isGoodJet",    "I", lenVar="nJet")
        self.out.branch("Jet_isGoodBJet",   "I", lenVar="nJet")
        self.out.branch("Jet_leptonClean",  "I", lenVar="nJet")

        self.out.branch("FatJet_isGoodJet", "I", lenVar="nFatJet")
        self.out.branch("FatJet_fromW", "I", lenVar="nFatJet")
        self.out.branch("FatJet_fromH", "I", lenVar="nFatJet")
        self.out.branch("FatJet_leptonClean",  "I", lenVar="nFatJet")

        #MIRANDA'S ADD-ONS
        self.out.branch("Tau_isVeto",       "F", lenVar="nTau")
        self.out.branch("nVetoTau",         "I")
        self.out.branch("IsoTrack_isVeto",  "F", lenVar="nIsoTrack")
        self.out.branch("nVetoIsoTrack",    "I")

        # Counter for good b-tags
        self.out.branch("nLepton",      "I")
        self.out.branch("nGoodJet",     "I")
        self.out.branch("nGoodBTag",    "I")


        self.out.branch("GoodJet_pt"    ,"F",  lenVar="nGoodJet")
        self.out.branch("GoodJet_eta"   ,"F",  lenVar="nGoodJet")
        self.out.branch("GoodJet_phi"   ,"F",  lenVar="nGoodJet")


        self.out.branch("BJet_pt"    ,"F",  lenVar="nBJet")
        self.out.branch("BJet_eta"   ,"F",  lenVar="nBJet")
        self.out.branch("BJet_phi"   ,"F",  lenVar="nBJet")

        self.out.branch("GoodFatJet_pt"    ,"F",  lenVar="nGoodFatJet")
        self.out.branch("GoodFatJet_eta"   ,"F",  lenVar="nGoodFatJet")
        self.out.branch("GoodFatJet_phi"   ,"F",  lenVar="nGoodFatJet")

        self.out.branch("FatJetfromW_pt"    ,"F",  lenVar="nFatJetfromW")
        self.out.branch("FatJetfromW_eta"   ,"F",  lenVar="nFatJetfromW")
        self.out.branch("FatJetfromW_phi"   ,"F",  lenVar="nFatJetfromW")

        self.out.branch("FatJetfromH_pt"    ,"F",  lenVar="nFatJetfromH")
        self.out.branch("FatJetfromH_eta"   ,"F",  lenVar="nFatJetfromH")
        self.out.branch("FatJetfromH_phi"   ,"F",  lenVar="nFatJetfromH")

        self.out.branch("bb_pt"               ,"F",  lenVar="nbb")
        self.out.branch("bb_eta"               ,"F",  lenVar="nbb")
        self.out.branch("bb_phi"               ,"F",  lenVar="nbb")
        self.out.branch("bb_mass"               ,"F",  lenVar="nbb")
        self.out.branch("MCT"               ,"F",  lenVar="nbb")

        self.out.branch("jj_pt"               ,"F",  lenVar="njj")
        self.out.branch("jj_eta"               ,"F",  lenVar="njj")
        self.out.branch("jj_phi"               ,"F",  lenVar="njj")
        self.out.branch("jj_mass"               ,"F",  lenVar="njj")
        self.out.branch("MCTjj"               ,"F",  lenVar="njj")




    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def isGoodJet(self, jet):
        return (jet.pt > 30 and abs(jet.eta)<2.4 and jet.jetId>0)

    def isGoodBJet(self, jet):
        if self.year == 2016:
            threshold = 0.6321
        if self.year == 2017:
            threshold = 0.4941
        if self.year == 2018:
            threshold = 0.4184
        return (self.isGoodJet(jet) and jet.btagDeepB > threshold)

    def isGoodFatJet(self, fatjet):
        return (fatjet.pt > 200 and fatjet.jetId>0)

    def isFatJetfromW(self, fatjet):
        return(self.isGoodFatJet(fatjet) and fatjet.deepTagMD_WvsQCD > 0.9)

    def isFatJetfromH(self, fatjet):
        if self.year == 2016:
            threshold = 0.8945
        if self.year == 2017:
            threshold = 0.8695
        if self.year == 2018:
            threshold = 0.8365
        return(self.isGoodFatJet(fatjet) and fatjet.deepTagMD_HbbvsQCD > threshold)

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

    #MIRANDA'S ADD ONS

    def isVetoTau(self, tau):
        return(tau.pt > 20 and abs(tau.eta) < 2.3 and tau.idDecayMode and tau.idMVAoldDM2017v1 == 8)
    
    def isVetoIsoTrack(self, isotrack):
        return(isotrack.pt > 10 and abs(isotrack.eta) < 2.4 and (isotrack.miniPFRelIso_all < (0.1*isotrack.pt) or isotrack.miniPFRelIso_all < 6))

    def invMass(self, o1, o2):
        v1 = ROOT.TLorentzVector()
        v2 = ROOT.TLorentzVector()
        v1.SetPtEtaPhi(o1['pt'], o1['eta'], o1['phi'], 0)
        v2.SetPtEtaPhi(o2['pt'], o2['eta'], o2['phi'], 0)
        return (v1+v2).M()

    def MCT2(self, b1, b2):
      return 2*b1.pt*b2.pt*(1+math.cos(self.deltaPhi(b1.phi,b2.phi)))

    def MCT(self, b1, b2):
      return math.sqrt(self.MCT2(b1,b2))

    def Mjj(self, b1, b2):
      bjet1 = ROOT.TLorentzVector()
      bjet1.SetPtEtaPhiM(b1.pt, b1.eta, b1.phi, 0)
      bjet2 = ROOT.TLorentzVector()
      bjet2.SetPtEtaPhiM(b2.pt, b2.eta, b2.phi, 0)
      return (bjet1 + bjet2).M()
    
    def Ptjj(self, b1, b2):
      bjet1 = ROOT.TLorentzVector()
      bjet1.SetPtEtaPhiM(b1.pt, b1.eta, b1.phi, 0)
      bjet2 = ROOT.TLorentzVector()
      bjet2.SetPtEtaPhiM(b2.pt, b2.eta, b2.phi, 0)
      return (bjet1 + bjet2).Pt()
    
    def Etajj(self, b1, b2):
      bjet1 = ROOT.TLorentzVector()
      bjet1.SetPtEtaPhiM(b1.pt, b1.eta, b1.phi, 0)
      bjet2 = ROOT.TLorentzVector()
      bjet2.SetPtEtaPhiM(b2.pt, b2.eta, b2.phi, 0)
      return (bjet1 + bjet2).Eta()
    
    def Phijj(self, b1, b2):
      bjet1 = ROOT.TLorentzVector()
      bjet1.SetPtEtaPhiM(b1.pt, b1.eta, b1.phi, 0)
      bjet2 = ROOT.TLorentzVector()
      bjet2.SetPtEtaPhiM(b2.pt, b2.eta, b2.phi, 0)
      return (bjet1 + bjet2).Phi()
    

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

#        print 'event %d *dr phil voice* shut the hell up bitch go kill yourself'%( event.event)

        muons       = Collection(event, "Muon")
        electrons   = Collection(event, "Electron")
        taus        = Collection(event, "Tau")
        jets        = Collection(event, "Jet")
        isotracks   = Collection(event, "IsoTrack")
        fatjets     = Collection(event, "FatJet")

        # MET
        met_pt  = event.MET_pt
        met_phi = event.MET_phi

        #create collection for taus, get quality cuts, similar to jets, if its already a lepton, not a veto tau anymore

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


        #TAU TIME
        isVetoTau   = []
        

        for t in taus:
            t.cleanmask = 1
            for coll in [electrons, muons]:
                for p in coll:
                    if self.deltaR(t, p) < 0.4:
                        t.cleanmask = 0
            isVetoTau.append(1 if (self.isVetoTau(t) and t.cleanmask) else 0)

        #TRACK TIME
            
        isVetoIsoTrack  = []
 
        for t in isotracks:
            t.cleanmask = 1
            for coll in [electrons, muons]:
                for p in coll:
                    if self.deltaR(t, p) > 0.4:
                        t.cleanmask = 0
            isVetoIsoTrack.append(1 if (self.isVetoIsoTrack(t)and t.cleanmask) else 0)


        cleanMaskV   = []
        isGoodJet    = []
        isGoodBJet   = []

        goodjets    = []
        bjets       = []

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

            # Fill the other collections
            if j.cleanMask:

                if self.isGoodJet(j):
                    goodjets.append({'pt':j.pt, 'eta':j.eta, 'phi':j.phi})
                    
                    if self.isGoodBJet(j):
                        bjets.append({'pt':j.pt, 'eta':j.eta, 'phi':j.phi})

        Mct = []
        bb = []
            
        ibb = 0
        for b,j in itertools.combinations(jets, 2):
            if self.isGoodBJet(b) and self.isGoodBJet(j):
                ibb += 1
                mct = self.MCT(b,j)
                Mbb = self.Mjj(b,j)
                ptbb = self.Ptjj(b,j)
                etabb = self.Etajj(b,j)
                phibb = self.Phijj(b,j)
                Mct.append(mct)
                bb.append({'pt': ptbb, 'eta': etabb, 'phi': phibb, 'mass': Mbb})
        
        Mctjj = []
        jj = []
            
        ijj = 0
        for b,j in itertools.combinations(jets, 2):
            if self.isGoodJet(b) and self.isGoodJet(j):
                ijj += 1
                mct = self.MCT(b,j)
                Mjj = self.Mjj(b,j)
                ptjj = self.Ptjj(b,j)
                etajj = self.Etajj(b,j)
                phijj = self.Phijj(b,j)
                Mctjj.append(mct)
                jj.append({'pt': ptjj, 'eta': etajj, 'phi': phijj, 'mass': Mjj})
        
        cleanMaskW  = []
        isGoodFatJet   = []
        isFatJetfromW  = []
        isFatJetfromH = []
        
        goodfatjets = []
        fatjetsfromW = []
        fatjetsfromH = []
        
        for f in fatjets:

            f.cleanMask = 1
            for coll in [electrons, muons]:
                for p in coll:
                    if p.isVeto:
                        if self.deltaR(f, p) < 0.4:
                            f.cleanMask = 0

            isGoodFatJet.append(1 if (self.isGoodFatJet(f) and f.cleanMask) else 0)
            isFatJetfromW.append(1 if (self.isFatJetfromW(f) and f.cleanMask) else 0)
            isFatJetfromH.append(1 if (self.isFatJetfromH(f) and f.cleanMask) else 0)
            
            cleanMaskW.append(f.cleanMask)

            # Fill the other collections
            if f.cleanMask:

                if self.isGoodFatJet(f):
                    goodfatjets.append({'pt':f.pt, 'eta':f.eta, 'phi':f.phi})
                    
                    if self.isFatJetfromW(f):
                        fatjetsfromW.append({'pt':f.pt, 'eta':f.eta, 'phi':f.phi})
                
                    if self.isFatJetfromH(f):
                        fatjetsfromH.append({'pt':f.pt, 'eta':f.eta, 'phi':f.phi})
                       

        # make sure the jets are properly sorted. they _should_ be sorted, but this can change once we reapply JECs if necessary
        bjets       = sorted(bjets, key = lambda i: i['pt'], reverse=True)
        goodjets    = sorted(goodjets, key = lambda i: i['pt'], reverse=True)
        goodfatjets = sorted(goodfatjets, key = lambda i: i['pt'], reverse=True)
        fatjetsfromW = sorted(fatjetsfromW, key = lambda i: i['pt'], reverse=True)
        fatjetsfromH = sorted(fatjetsfromH, key = lambda i: i['pt'], reverse=True)



        self.out.fillBranch("Muon_isTight",         isTightMuon)
        self.out.fillBranch("Muon_isVeto",          isVetoMuon)
        self.out.fillBranch("Electron_isTight",     isTightElectron)
        self.out.fillBranch("Electron_isVeto",      isVetoElectron)
        self.out.fillBranch("Jet_leptonClean",      cleanMaskV)
        self.out.fillBranch("Jet_isGoodJet",        isGoodJet)
        self.out.fillBranch("Jet_isGoodBJet",       isGoodBJet)
        self.out.fillBranch("FatJet_leptonClean",   cleanMaskW)
        self.out.fillBranch("FatJet_isGoodJet",     isGoodFatJet)
        self.out.fillBranch("FatJet_fromW",         isFatJetfromW)
        self.out.fillBranch("FatJet_fromH",         isFatJetfromH)
        self.out.fillBranch("nGoodBTag",            sum(isGoodBJet))
        self.out.fillBranch("nGoodJet",             sum(isGoodJet))

        self.out.fillBranch("Tau_isVeto",           isVetoTau)
        self.out.fillBranch("nVetoTau",             sum(isVetoTau))
        self.out.fillBranch("IsoTrack_isVeto",      isVetoIsoTrack)
        self.out.fillBranch("nVetoIsoTrack",        sum(isVetoIsoTrack))


        goodjets_pd = pd.DataFrame(goodjets)
        self.out.fillBranch("nGoodJet",          len(goodjets_pd) )
        if len(goodjets_pd)>0:
            self.out.fillBranch("GoodJet_pt",        goodjets_pd.sort_values(by='pt', ascending=False)['pt'].tolist() )
            self.out.fillBranch("GoodJet_eta",       goodjets_pd.sort_values(by='pt', ascending=False)['eta'].tolist() )
            self.out.fillBranch("GoodJet_phi",       goodjets_pd.sort_values(by='pt', ascending=False)['phi'].tolist() )
            
        
        bjets_pd = pd.DataFrame(bjets)
        self.out.fillBranch("nBJet",          len(bjets_pd) )
        if len(bjets_pd)>0:
            self.out.fillBranch("BJet_pt",        bjets_pd.sort_values(by='pt', ascending=False)['pt'].tolist() )
            self.out.fillBranch("BJet_eta",       bjets_pd.sort_values(by='pt', ascending=False)['eta'].tolist() )
            self.out.fillBranch("BJet_phi",       bjets_pd.sort_values(by='pt', ascending=False)['phi'].tolist() )

        bb_pd = pd.DataFrame(bb)
        self.out.fillBranch("nbb",          ibb )
        if len(bb_pd)>0:
            self.out.fillBranch("bb_pt",        bb_pd.sort_values(by='pt', ascending=False)['pt'].tolist() )
            self.out.fillBranch("bb_eta",       bb_pd.sort_values(by='pt', ascending=False)['eta'].tolist() )
            self.out.fillBranch("bb_phi",       bb_pd.sort_values(by='pt', ascending=False)['phi'].tolist() )
            self.out.fillBranch("bb_mass",       bb_pd.sort_values(by='pt', ascending=False)['mass'].tolist() )
            self.out.fillBranch("MCT",          Mct)

        jj_pd = pd.DataFrame(jj)
        self.out.fillBranch("njj",          ijj )
        if len(jj_pd)>0:
            self.out.fillBranch("jj_pt",        jj_pd.sort_values(by='pt', ascending=False)['pt'].tolist() )
            self.out.fillBranch("jj_eta",       jj_pd.sort_values(by='pt', ascending=False)['eta'].tolist() )
            self.out.fillBranch("jj_phi",       jj_pd.sort_values(by='pt', ascending=False)['phi'].tolist() )
            self.out.fillBranch("jj_mass",       jj_pd.sort_values(by='pt', ascending=False)['mass'].tolist() )
            self.out.fillBranch("MCTjj",          Mctjj)

        goodfatjets_pd = pd.DataFrame(goodfatjets)
        self.out.fillBranch("nGoodFatJet",          len(goodfatjets_pd) )
        if len(goodfatjets_pd)>0:
            self.out.fillBranch("GoodFatJet_pt",        goodfatjets_pd.sort_values(by='pt', ascending=False)['pt'].tolist() )
            self.out.fillBranch("GoodFatJet_eta",       goodfatjets_pd.sort_values(by='pt', ascending=False)['eta'].tolist())
            self.out.fillBranch("GoodFatJet_phi",       goodfatjets_pd.sort_values(by='pt', ascending=False)['phi'].tolist())

        fatjetsfromW_pd = pd.DataFrame(fatjetsfromW)
        self.out.fillBranch("nFatJetfromW",          len(fatjetsfromW_pd) )
        if len(fatjetsfromW_pd)>0:
            self.out.fillBranch("FatJetfromW_pt",       fatjetsfromW_pd.sort_values(by='pt', ascending=False)['pt'].tolist())
            self.out.fillBranch("FatJetfromW_eta",      fatjetsfromW_pd.sort_values(by='pt', ascending=False)['eta'].tolist() )
            self.out.fillBranch("FatJetfromW_phi",      fatjetsfromW_pd.sort_values(by='pt', ascending=False)['phi'].tolist() )

        fatjetsfromH_pd = pd.DataFrame(fatjetsfromH)
        self.out.fillBranch("nFatJetfromH",          len(fatjetsfromH_pd) )
        if len(fatjetsfromH_pd)>0:
            self.out.fillBranch("FatJetfromH_pt",       fatjetsfromH_pd.sort_values(by='pt', ascending=False)['pt'].tolist() )
            self.out.fillBranch("FatJetfromH_eta",      fatjetsfromH_pd.sort_values(by='pt', ascending=False)['eta'].tolist() )
            self.out.fillBranch("FatJetfromH_phi",      fatjetsfromH_pd.sort_values(by='pt', ascending=False)['phi'].tolist() )


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

selector2016 = lambda : PhysicsObjects( year=2016 )
selector2017 = lambda : PhysicsObjects( year=2017 )
selector2018 = lambda : PhysicsObjects( year=2018 )

