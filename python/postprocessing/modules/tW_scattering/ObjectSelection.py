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

class Object():
    def __init__(self):
        pass

    @classmethod
    def fromDict(cls, d):
        obj = cls()
        for key in d.keys():
            setattr(obj, key, d[key])
        return obj

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
        self.out.branch("Lepton_mass", "F", lenVar="nLepton")
        self.out.branch("Lepton_pdgId", "I", lenVar="nLepton")
        self.out.branch("Lepton_miniIso", "F", lenVar="nLepton")
        self.out.branch("Lepton_muIndex", "I", lenVar="nLepton")
        self.out.branch("Lepton_elIndex", "I", lenVar="nLepton")

        self.out.branch("RecoW_pt", "F", lenVar="nRecoW")
        self.out.branch("RecoW_eta", "F", lenVar="nRecoW")
        self.out.branch("RecoW_phi", "F", lenVar="nRecoW")
        self.out.branch("RecoW_mass", "F", lenVar="nRecoW")
        self.out.branch("RecoW_genPt", "F", lenVar="nRecoW")
        self.out.branch("RecoW_genEta", "F", lenVar="nRecoW")
        self.out.branch("RecoW_genPhi", "F", lenVar="nRecoW")
        self.out.branch("RecoW_genMass", "F", lenVar="nRecoW")
        self.out.branch("RecoW_qglSum", "F", lenVar="nRecoW")
        self.out.branch("RecoW_qglProd", "F", lenVar="nRecoW")


        ## New collection of Jets. #FIXME overlap removal with GoodLeptons
        #self.out.branch("GoodJet_pt", "F", lenVar="nGoodJet")
        #self.out.branch("GoodJet_eta", "F", lenVar="nGoodJet")
        #self.out.branch("GoodJet_phi", "F", lenVar="nGoodJet")
        #self.out.branch("GoodJet_btag", "F", lenVar="nGoodJet")

        ##
        self.out.branch("Muon_isVeto",      "I", lenVar="nMuon")
        self.out.branch("Muon_isTight",     "I", lenVar="nMuon")
        self.out.branch("Electron_isVeto",  "I", lenVar="nElectron")
        self.out.branch("Electron_isTight", "I", lenVar="nElectron")
        self.out.branch("Jet_isGoodJet",    "I", lenVar="nJet")
        self.out.branch("Jet_isGoodJetAll", "I", lenVar="nJet")
        self.out.branch("Jet_isGoodBJet",   "I", lenVar="nJet")
        self.out.branch("Jet_leptonClean",  "I", lenVar="nJet")
        self.out.branch("Jet_WIdx",  "I", lenVar="nJet")


        # Counter for good b-tags
        self.out.branch("nLepton",      "I")
        self.out.branch("nVetoLepton",  "I")
        self.out.branch("nGoodJet",     "I")
        self.out.branch("nGoodBTag",    "I")

        self.out.branch("isSingleLep",  "I")
        self.out.branch("isDiLep",      "I")
        self.out.branch("isTriLep",     "I")
        self.out.branch("isSS",         "I")

        self.out.branch("diWness",      "F")
        self.out.branch("MT",      "F")
        self.out.branch("MT_puppi",      "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def isGoodJet(self, jet):
        return (jet.pt > 25 and abs(jet.eta)<2.4 and jet.jetId>1)

    def isFwdJet(self, jet):
        return (((jet.pt > 25 and abs(jet.eta)<2.7) or (jet.pt>60 and abs(jet.eta)>=2.7 and abs(jet.eta)<3.0) or (jet.pt>25 and abs(jet.eta)>=3.0 and abs(jet.eta)<5.0))  and jet.jetId>1)

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

    def invMass(self, o1, o2):
        v1 = ROOT.TLorentzVector()
        v2 = ROOT.TLorentzVector()
        v1.SetPtEtaPhiM(o1['pt'], o1['eta'], o1['phi'], o1['mass'])
        v2.SetPtEtaPhiM(o2['pt'], o2['eta'], o2['phi'], o2['mass'])
        return (v1+v2).M(), (v1+v2).Pt(), (v1+v2).Eta(), (v1+v2).Phi()

    def getW(self, jet1, jet2):
        mass, pt, eta, phi = self.invMass(jet1,jet2)
        qgl_sum = jet1['qgl'] + jet2['qgl']
        qgl_prod = jet1['qgl'] * jet2['qgl']
        return {'chi2': abs(mass-80.)/20., 'mass':mass, 'pt':pt, 'eta':eta, 'phi':phi, 'qgl_sum':qgl_sum, 'qgl_prod': qgl_prod}

    def getWcandidates(self, jets):
        uniqueCombs = [(0,1,2,3), (1,2,3,0), (0,2,1,3)]
        minChi2 = 9999
        jetIndices = range(len(jets))
        sets = [ comb for comb in itertools.combinations(jetIndices, 4) ]
        W_cand = []
        #print "Number of jets:", len(jets)
        #print "All combinations of the 4 jets", sets
        for s in sets:
            # selection of 4 jets
            for combs in uniqueCombs:
                # combination of the 4 jets
                indices = [ s[x] for x in combs ]
                #print "One of the combinations", indices
                W = [self.getW(jets[indices[0]], jets[indices[1]]), self.getW(jets[indices[2]], jets[indices[3]])]
                chi2 = W[0]['chi2'] + W[1]['chi2']
                W_cand.append({'chi2': chi2, 'indices':indices, 'W':W})

        W_cand = sorted(W_cand, key=lambda x: x['chi2'])

        return W_cand
            
    def getRealWs(self, jets, genWs):
        combs = [ comb for comb in itertools.combinations(jets, 2) ]
        Wcands = []
        for j1, j2 in combs:
            Wcand = self.getW(j1, j2)
            if j1['WIdx']==j2['WIdx'] and j1['WIdx']>-1:
                Wcand['genPt'] = genWs[j1['WIdx']].pt
                Wcand['genEta'] = genWs[j1['WIdx']].eta
                Wcand['genPhi'] = genWs[j1['WIdx']].phi
                Wcand['genMass'] = genWs[j1['WIdx']].mass
            else:
                Wcand['genPt']      = -9999
                Wcand['genEta']     = -9999
                Wcand['genPhi']     = -9999
                Wcand['genMass']    = -9999
            Wcands.append(Wcand)

        return Wcands
                


    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        muons       = Collection(event, "Muon")
        electrons   = Collection(event, "Electron")
        jets        = Collection(event, "Jet")
        genjets     = Collection(event, "GenJet")
        genW        = Collection(event, "W")
        
        # MET
        met_pt  = event.MET_pt
        met_phi = event.MET_phi

        # Puppi MET
        puppi_met_pt = event.PuppiMET_pt
        puppi_met_phi = event.PuppiMET_phi

        # tight lepton collection, will be sorted by pt
        leptons     = []
        vleptons    = []

        isTightMuon = []
        isVetoMuon = []
        for i,mu in enumerate(muons):
            mu.isTight  = self.isTightMuon(mu)
            mu.isVeto   = self.isVetoMuon(mu)
            isTightMuon.append(mu.isTight)
            isVetoMuon.append(mu.isVeto)
            if self.isTightMuon(mu):
                leptons.append({'pt':mu.pt, 'eta':mu.eta, 'phi':mu.phi, 'pdgId':mu.pdgId, 'miniIso':mu.miniPFRelIso_all, 'muIndex':i, 'elIndex':-1, 'mass':mu.mass})
            if self.isVetoMuon(mu):
                vleptons.append({'pt':mu.pt, 'eta':mu.eta, 'phi':mu.phi, 'pdgId':mu.pdgId, 'miniIso':mu.miniPFRelIso_all, 'muIndex':i, 'elIndex':-1, 'mass':mu.mass})


        isTightElectron = []
        isVetoElectron = []
        for i,el in enumerate(electrons):
            el.isTight  = self.isTightElectron(el)
            el.isVeto   = self.isVetoElectron(el)
            isTightElectron.append(el.isTight)
            isVetoElectron.append(el.isVeto)
            if self.isTightElectron(el):
                leptons.append({'pt':el.pt, 'eta':el.eta, 'phi':el.phi, 'pdgId':el.pdgId, 'miniIso':el.miniPFRelIso_all, 'muIndex':-1, 'elIndex':i, 'mass':el.mass})
            if self.isVetoElectron(el):
                vleptons.append({'pt':el.pt, 'eta':el.eta, 'phi':el.phi, 'pdgId':el.pdgId, 'miniIso':el.miniPFRelIso_all, 'muIndex':-1, 'elIndex':i, 'mass':el.mass})

        leptons = sorted(leptons, key = lambda i: i['pt'], reverse=True)
        
        isSingleLep = len(leptons)==1
        isDiLep     = len(leptons)==2
        isTriLep    = len(leptons)==3

        isSS = False
        if len(leptons)>1:
            isSS = leptons[0]['pdgId']*leptons[1]['pdgId']>0

        MT = -1
        MT_puppi = -1
        if len(leptons)>0:
            MT = math.sqrt(2*leptons[0]['pt']*met_pt*(1-math.cos(leptons[0]['phi']-met_phi)))
            MT_puppi = math.sqrt(2*leptons[0]['pt']*puppi_met_pt*(1-math.cos(leptons[0]['phi']-puppi_met_phi)))

        cleanMaskV  = []
        isGoodJet   = []
        isGoodBJet  = []
        isGoodJetAll = []

        fwdjets     = []
        jets_out    = []
        bjets       = []
        nonbjets    = []
        alljets = []
        WIdx = []

        for j in jets:

            j.cleanMask = 1
            for coll in [electrons, muons]:
                for p in coll:
                    if p.isVeto:
                        if self.deltaR(j, p) < 0.4:
                            j.cleanMask = 0

            isGoodJet.append(1 if (self.isGoodJet(j) and j.cleanMask) else 0)
            isGoodJetAll.append(1 if (self.isFwdJet(j) and j.cleanMask) else 0)
            isGoodBJet.append(1 if (self.isGoodBJet(j) and j.cleanMask) else 0)
            
            cleanMaskV.append(j.cleanMask)

            # now get some info from the GenJets
            j.WIdx = -1
            if j.genJetIdx>-1 and j.genJetIdx<len(genjets): # apparently lowest pt genjets are not stored...
                j.WIdx = genjets[j.genJetIdx].WIdx
            WIdx.append(j.WIdx)

            alljets.append({'pt':j.pt, 'eta':j.eta, 'phi':j.phi, 'qgl':j.qgl, 'WIdx':j.WIdx})

            # Fill the other collections
            if j.cleanMask:

                if self.isGoodJet(j):
                    jets_out.append({'pt':j.pt_nom, 'eta':j.eta, 'phi':j.phi, 'qgl':j.qgl, 'WIdx':j.WIdx, 'mass':j.mass_nom})
                    
                    if self.isGoodBJet(j):
                        bjets.append({'pt':j.pt, 'eta':j.eta, 'phi':j.phi})
                    else:
                        nonbjets.append({'pt':j.pt_nom, 'eta':j.eta, 'phi':j.phi, 'qgl':j.qgl, 'WIdx':j.WIdx, 'mass':j.mass_nom})
                
                if self.isFwdJet(j):
                    fwdjets.append({'pt':j.pt, 'eta':j.eta, 'phi':j.phi})

            #for idx, gjet in enumerate(genjets):
            #    if 
                       

        # make sure the jets are properly sorted. they _should_ be sorted, but this can change once we reapply JECs if necessary
        bjets       = sorted(bjets, key = lambda i: i['pt'], reverse=True)
        nonbjets    = sorted(nonbjets, key = lambda i: i['pt'], reverse=True)
        fwdjets     = sorted(fwdjets, key = lambda i: i['pt'], reverse=True) # all jets, including forward ones
        jets_out    = sorted(jets_out, key = lambda i: i['pt'], reverse=True)

        # calculate MT, Mlb, Mlb_max, M_jj_b1, M_jj_b2 etc
        Mlb_00 = -99 # leading lepton, first b-jet
        Mlb_01 = -99 # leading lepton, second b-jet
        Mlb_10 = -99 # trailing lepton, first b-jet
        Mlb_11 = -99 # trailing lepton, second b-jet

        #if len(jets)>0:
        #    if len(leptons)>0:
        #        Mlb_00 = 
        #
        #if len(leptons)>1

        ## W candidates
        # get any combination of 4 non b-jets
        # this is *very* inefficient. Need to think of a better way to reconstruct two Ws
        recoWs = []
        if len(jets_out)>3:
            #W_cands = self.getWcandidates(nonbjets)
            recoWs = self.getRealWs(jets_out, genW)
            

        ### This doesn't really work
        diWness = 9999.
        #recoWs = []
        #perfectMatch = False
        #if len(W_cands)>0:
        #    for W_cand in W_cands:
        #        genMatch1 = False
        #        genMatch2 = False
        #        for W in genW:
        #            if not genMatch1:
        #                genMatch1 = self.deltaR(Object.fromDict(W_cand['W'][0]), W)<0.4
        #            if not genMatch2:
        #                genMatch2 = self.deltaR(Object.fromDict(W_cand['W'][1]), W)<0.4
        #        if genMatch1 and genMatch2:
        #            perfectMatch = W_cand
        # #           print "Nice, found the right pair"
        #            break
        #    
        #    for recoW in W_cands[0]['W']:
        #        genMatch = False
        #        for W in genW:
        #            genMatch = self.deltaR(Object.fromDict(recoW), W)<0.4
        #            if genMatch: break
        #        recoWs.append({'pt':recoW['pt'], 'eta':recoW['eta'], 'phi':recoW['phi'], 'mass':recoW['mass'], 'genMatch':genMatch})
        #        diWness = W_cands[0]['chi2']


        self.out.fillBranch("Muon_isTight",     isTightMuon)
        self.out.fillBranch("Muon_isVeto",      isVetoMuon)
        self.out.fillBranch("Electron_isTight", isTightElectron)
        self.out.fillBranch("Electron_isVeto",  isVetoElectron)
        self.out.fillBranch("Jet_leptonClean",   cleanMaskV)
        self.out.fillBranch("Jet_isGoodJet",    isGoodJet)
        self.out.fillBranch("Jet_isGoodJetAll", isGoodJetAll)
        self.out.fillBranch("Jet_isGoodBJet",   isGoodBJet)
        self.out.fillBranch("Jet_WIdx",   WIdx)
        self.out.fillBranch("nGoodBTag",        sum(isGoodBJet))
        self.out.fillBranch("nGoodJet",         sum(isGoodJet))

        self.out.fillBranch("isSingleLep",      isSingleLep)
        self.out.fillBranch("isDiLep",          isDiLep)
        self.out.fillBranch("isTriLep",         isTriLep)
        self.out.fillBranch("isSS",             isSS)
        self.out.fillBranch("diWness",          diWness)
        self.out.fillBranch("MT",               MT)
        self.out.fillBranch("MT_puppi",         MT_puppi)

        # make pandas dataframe out of list
        leptons_pd = pd.DataFrame(leptons)

        self.out.fillBranch("nVetoLepton",          len(vleptons) )
        self.out.fillBranch("nLepton",          len(leptons_pd) )
        if len(leptons_pd)>0:
            self.out.fillBranch("Lepton_pt",        leptons_pd.sort_values(by='pt', ascending=False)['pt'].tolist() )
            self.out.fillBranch("Lepton_eta",       leptons_pd.sort_values(by='pt', ascending=False)['eta'].tolist() )
            self.out.fillBranch("Lepton_phi",       leptons_pd.sort_values(by='pt', ascending=False)['phi'].tolist() )
            self.out.fillBranch("Lepton_mass",      leptons_pd.sort_values(by='pt', ascending=False)['mass'].tolist() )
            self.out.fillBranch("Lepton_pdgId",     leptons_pd.sort_values(by='pt', ascending=False)['pdgId'].tolist() )
            self.out.fillBranch("Lepton_miniIso",   leptons_pd.sort_values(by='pt', ascending=False)['miniIso'].tolist() )
            self.out.fillBranch("Lepton_muIndex",   leptons_pd.sort_values(by='pt', ascending=False)['muIndex'].tolist() )
            self.out.fillBranch("Lepton_elIndex",   leptons_pd.sort_values(by='pt', ascending=False)['elIndex'].tolist() )

        recoWs_pd = pd.DataFrame(recoWs)

        self.out.fillBranch("nRecoW",          len(recoWs_pd) )
        if len(recoWs_pd)>0:
            self.out.fillBranch("RecoW_pt",        recoWs_pd['pt'].tolist() )
            self.out.fillBranch("RecoW_eta",       recoWs_pd['eta'].tolist() )
            self.out.fillBranch("RecoW_phi",       recoWs_pd['phi'].tolist() )
            self.out.fillBranch("RecoW_mass",      recoWs_pd['mass'].tolist() )
            self.out.fillBranch("RecoW_genPt",     recoWs_pd['genPt'].tolist() )
            self.out.fillBranch("RecoW_genEta",     recoWs_pd['genEta'].tolist() )
            self.out.fillBranch("RecoW_genPhi",     recoWs_pd['genPhi'].tolist() )
            self.out.fillBranch("RecoW_genMass",     recoWs_pd['genMass'].tolist() )
            self.out.fillBranch("RecoW_qglSum",     recoWs_pd['qgl_sum'].tolist() )
            self.out.fillBranch("RecoW_qglProd",     recoWs_pd['qgl_prod'].tolist() )


        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

selector2018 = lambda : PhysicsObjects( year=2018 )
