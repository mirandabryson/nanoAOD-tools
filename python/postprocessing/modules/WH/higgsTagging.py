import ROOT
import os
import numpy as np
import pandas as pd
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class higgsTagging(Module):
    def __init__(self, year=2016, isData=False):
        self.year = year
        self.isData = isData
        if year == 2016:
            self.htagWP = 0.8945
        elif year == 2017:
            self.htagWP = 0.8695
        elif year == 2018:
            self.htagWP = 0.8365
        else:
            print "Don't know year %s"%year


    def beginJob(self):
        self.SFs = pd.DataFrame.from_csv(os.path.expandvars('$CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/modules/WH/deepak8v2_bbvslight.csv'), header=0, sep=' ', parse_dates=False, index_col=None)

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("hasNano", "I")
        self.out.branch("nHiggs", "I")
        if not self.isData:
            self.out.branch("jup_nHiggs", "I")
            self.out.branch("jdown_nHiggs", "I")
        self.out.branch("w_higgsSF", "F")
        self.out.branch("w_higgsSFUp", "F")
        self.out.branch("w_higgsSFDown", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def getSF(self, pt, sigma=0):
        SF = 1
        if sigma == 0:
            SF = float(self.SFs[(self.SFs['year']==self.year) & (self.SFs['wp']=='mp') & (self.SFs['ptmin']<pt) & (self.SFs['ptmax']>=pt)]['sf'])
        if sigma == -1:
            SF = float(self.SFs[(self.SFs['year']==self.year) & (self.SFs['wp']=='mp') & (self.SFs['ptmin']<pt) & (self.SFs['ptmax']>=pt)]['sf']) - float(self.SFs[(self.SFs['year']==self.year) & (self.SFs['wp']=='mp') & (self.SFs['ptmin']<pt) & (self.SFs['ptmax']>=pt)]['sflow'])
        if sigma == 1:
            SF = float(self.SFs[(self.SFs['year']==self.year) & (self.SFs['wp']=='mp') & (self.SFs['ptmin']<pt) & (self.SFs['ptmax']>=pt)]['sf']) + float(self.SFs[(self.SFs['year']==self.year) & (self.SFs['wp']=='mp') & (self.SFs['ptmin']<pt) & (self.SFs['ptmax']>=pt)]['sfhigh'])
        return SF


    def analyze(self, event):

        """process event, return True (go to next module) or False (fail, go to next event)"""
        jets = Collection(event, "FatJet")
        ptmin = 250

        nHiggs = 0
        jup_nHiggs = 0
        jdown_nHiggs = 0

        w_higgsSF = 1
        w_higgsSFUp = 1
        w_higgsSFDown = 1

        for j in jets:
            if j.pt_nom > ptmin and j.deepTagMD_HbbvsQCD > self.htagWP:
                nHiggs += 1
                if not self.isData:
                    w_higgsSF *= self.getSF(j.pt_nom, sigma=0)
                    w_higgsSFUp *= self.getSF(j.pt_nom, sigma=1)
                    w_higgsSFDown *= self.getSF(j.pt_nom, sigma=-1)
            if not self.isData:
                if j.pt_jesTotalUp > ptmin and j.deepTagMD_HbbvsQCD > self.htagWP:
                    jup_nHiggs += 1
                if j.pt_jesTotalDown > ptmin and j.deepTagMD_HbbvsQCD > self.htagWP:
                    jdown_nHiggs += 1
            
        self.out.fillBranch("nHiggs", nHiggs)
        if not self.isData:
            self.out.fillBranch("jup_nHiggs", jup_nHiggs)
            self.out.fillBranch("jdown_nHiggs", jdown_nHiggs)

        self.out.fillBranch("w_higgsSF", w_higgsSF)
        self.out.fillBranch("w_higgsSFUp", w_higgsSFUp)
        self.out.fillBranch("w_higgsSFDown", w_higgsSFDown)

        self.out.fillBranch("hasNano", 1)

        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

higgsSF16 = lambda : higgsTagging( 2016 )
higgsSF17 = lambda : higgsTagging( 2017 )
higgsSF18 = lambda : higgsTagging( 2018 )

