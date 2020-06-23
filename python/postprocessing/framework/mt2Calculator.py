import ROOT, array
from math import pi, sqrt, cos, sin

#wrapper class for MT2 variables
class MT2Calculator:
    def __init__(self):
        ROOT.gROOT.ProcessLine(".L /home/users/mbryson/CMSSW_10_2_9/src/NanoAODTools/python/postprocessing/framework/mt2_bisect.cpp")
        self.mt2 = ROOT.mt2()
        self.reset()
        self.leptonMass = 0.
        self.jetMass = 0.
        self.bjetMass = 0.
        self.mt2Mass_ll   = 0.   #probe mass for the daughter system
        self.mt2Mass_jj   = 0.
        self.mt2Mass_bb   = 80.4 #probe mass for the daughter system
        self.mt2Mass_blbl = 0.   #probe mass for the daughter system
        self.mt2Mass_bbjj = 0.
        self.mt2Mass_bjjb = 0.
    def reset(self):
        self.met=None
        self.lepton1=None
        self.lepton2=None
        self.bjet1=None
        self.bjet2=None

#Setters
    def setMet(self, pt, phi):
        self.met = ROOT.TVector2(pt*cos(phi), pt*sin(phi))
    def setJet1(self, pt1, eta1, phi1):
        self.jet1 = ROOT.TLorentzVector()
        self.jet1.SetPtEtaPhiM(pt1, eta1, phi1, self.jetMass)
    def setJet2(self, pt2, eta2, phi2):
        self.jet2 = ROOT.TLorentzVector()
        self.jet2.SetPtEtaPhiM(pt2, eta2, phi2, self.jetMass)
    def setJets(self, pt1, eta1, phi1, pt2, eta2, phi2):
        self.setJet1(pt1,eta1,phi1)
        self.setJet2(pt2,eta2,phi2)
    def setLepton1(self, pt1, eta1, phi1):
        self.lepton1 = ROOT.TLorentzVector()
        self.lepton1.SetPtEtaPhiM(pt1, eta1, phi1, self.leptonMass)
    def setLepton2(self, pt2, eta2, phi2):
        self.lepton2 = ROOT.TLorentzVector()
        self.lepton2.SetPtEtaPhiM(pt2, eta2, phi2, self.leptonMass)
    def setLeptons(self, pt1, eta1, phi1, pt2, eta2, phi2):
        self.setLepton1(pt1,eta1,phi1)
        self.setLepton2(pt2,eta2,phi2)
    def setBJet1(self, pt1, eta1, phi1):
        self.bjet1 = ROOT.TLorentzVector()
        self.bjet1.SetPtEtaPhiM(pt1, eta1, phi1, self.bjetMass)
    def setBJet2(self, pt2, eta2, phi2):
        self.bjet2 = ROOT.TLorentzVector()
        self.bjet2.SetPtEtaPhiM(pt2, eta2, phi2, self.bjetMass)
    def setBJets(self, pt1, eta1, phi1, pt2, eta2, phi2):
        self.setBJet1(pt1,eta1,phi1)
        self.setBJet2(pt2,eta2,phi2)

#Traditional MT2
    def mt2ll(self):
        assert self.met and self.lepton1 and self.lepton2, "Incomplete specification, need met/lepton1/lepton2"
        pmiss  = array.array('d',[  0., self.met.Px(), self.met.Py()] )
        l1     = array.array('d',[  0., self.lepton1.Px(), self.lepton1.Py()] )
        l2     = array.array('d',[  0., self.lepton2.Px(), self.lepton2.Py()] )
        self.mt2.set_mn(self.mt2Mass_ll)
        self.mt2.set_momenta(l1, l2, pmiss)
        return self.mt2.get_mt2()

#MT2jj
    def mt2jj(self):
        assert self.met and self.jet1 and self.jet2, "Incomplete specification, need met/jet1/jet2"
        pmiss  = array.array('d',[  0., self.met.Px(), self.met.Py()] )
        jet1     = array.array('d',[  0., self.jet1.Px(), self.jet1.Py()] )
        jet2     = array.array('d',[  0., self.jet2.Px(), self.jet2.Py()] )
        self.mt2.set_mn(self.mt2Mass_jj)
        self.mt2.set_momenta(jet1, jet2, pmiss)
        return self.mt2.get_mt2()

#MT2bb (treating leptons invisibly, endpoint at top mass)
    def mt2bb(self):
        assert self.met and self.lepton1 and self.lepton2 and self.bjet1 and self.bjet2, "Incomplete specification, need met/lepton1/lepton2/bjet1/bjet2"
        pmiss_vec = self.met+ROOT.TVector2(self.lepton1.Px()+self.lepton2.Px(), self.lepton1.Py()+self.lepton2.Py())
        pmiss  = array.array('d',[  0., pmiss_vec.Px(), pmiss_vec.Py()] )
        b1     = array.array('d',[  0., self.bjet1.Px(), self.bjet1.Py()] )
        b2     = array.array('d',[  0., self.bjet2.Px(), self.bjet2.Py()] )
        self.mt2.set_mn(self.mt2Mass_bb)
        self.mt2.set_momenta(b1, b2, pmiss)
        return self.mt2.get_mt2()
#MT2blbl (Brians variant)
    def mt2blbl(self, strategy="minMaxMass"):
        assert self.met and self.lepton1 and self.lepton2 and self.bjet1 and self.bjet2, "Incomplete specification, need met/lepton1/lepton2/bjet1/bjet2"

        #select lepton/bjet pairing by minimizing maximum mass
        if strategy=="minMaxMass":
            max1 = max([(self.lepton1 + self.bjet1).M(), (self.lepton2 + self.bjet2).M()])
            max2 = max([(self.lepton1 + self.bjet2).M(), (self.lepton2 + self.bjet1).M()])
            if max1<max2: #Choose pairing with smaller invariant mass
                bl1,bl2 = self.lepton1+self.bjet1, self.lepton2+self.bjet2
            else:
                bl1,bl2 = self.lepton1+self.bjet2, self.lepton2+self.bjet1
        else:
            assert False, "only minMaxMass implemented"

        pmiss  = array.array('d',[  0., self.met.Px(), self.met.Py()] )
        bl1     = array.array('d',[  0., bl1.Px(), bl1.Py()] )
        bl2     = array.array('d',[  0., bl2.Px(), bl2.Py()] )
        self.mt2.set_mn(self.mt2Mass_blbl)
        self.mt2.set_momenta(bl1, bl2, pmiss)
        return self.mt2.get_mt2()

#MT2bbjj
    def mt2bbjj(self, strategy="minMaxMass"):
        assert self.met and self.jet1 and self.jet2 and self.bjet1 and self.bjet2, "Incomplete specification, need met/lepton1/lepton2/bjet1/bjet2"

        #select lepton/bjet pairing by minimizing maximum mass
        if strategy=="minMaxMass":
            max1 = max([(self.bjet1 + self.bjet2).M(), (self.jet1 + self.jet2).M()])
            max2 = max([(self.jet1 + self.bjet2).M(), (self.jet2 + self.bjet1).M()])
            max3 = max([(self.bjet1 + self.jet1).M(), (self.bjet2 + self.jet2).M()])
            if (max1<max2 and max1<max3): #Choose pairing with smallest invariant mass
                b1b2,j1j2 = self.bjet1+self.bjet2, self.jet2+self.bjet2
            elif (max2<max1 and max2<max3):
                b1b2,j1j2 = self.jet1+self.bjet2, self.jet2+self.bjet1
            elif (max3<max1 and max3<max2):
                b1b2,j1j2 = self.bjet1+self.jet1, self.bjet2+self.jet2
        else:
            assert False, "only minMaxMass implemented"

        pmiss  = array.array('d',[  0., self.met.Px(), self.met.Py()] )
        b1b2     = array.array('d',[  0., b1b2.Px(), b1b2.Py()] )
        j1j2     = array.array('d',[  0., j1j2.Px(), j1j2.Py()] )
        self.mt2.set_mn(self.mt2Mass_bbjj)
        self.mt2.set_momenta(b1b2, j1j2, pmiss)
        return self.mt2.get_mt2()

    def mt2bjjb(self, strategy="minMaxMass"):
        assert self.met and self.jet1 and self.jet2 and self.bjet1 and self.bjet2, "Incomplete specification, need met/jet1/jet2/bjet1/bjet2"

        #select bjetjetjet/bjet pairing by minimizing maximum mass
        if strategy=="minMaxMass":
            max1 = max([(self.bjet1 + self.jet1 + self.jet2).M(), (self.bjet2).M()])
            max2 = max([(self.jet1 + self.jet2 + self.bjet2).M(), (self.bjet1).M()])
            if max1<max2: #Choose pairing with smaller invariant mass
                b1jj,b2jj = self.jet1+self.jet2+self.bjet1, self.bjet2
            else:
                b1jj,b2jj = self.jet1+self.jet2+self.bjet2, self.bjet1
        else:
            assert False, "only minMaxMass implemented"
            
        pmiss  = array.array('d',[  0., self.met.Px(), self.met.Py()] )
        b1jj     = array.array('d',[  0., b1jj.Px(), b1jj.Py()] )
        b2jj     = array.array('d',[  0., b2jj.Px(), b2jj.Py()] )
        self.mt2.set_mn(self.mt2Mass_bjjb)
        self.mt2.set_momenta(b1jj, b2jj, pmiss)
        return self.mt2.get_mt2()

mt2Calculator = MT2Calculator()
