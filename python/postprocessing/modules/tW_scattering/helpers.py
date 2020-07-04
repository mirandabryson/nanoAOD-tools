#!/usr/bin/env python                                                                              
import os, sys
import subprocess

from PhysicsTools.NanoAODTools.postprocessing.modules.tW_scattering.ObjectSelection import *

def chooseselector(datayear):
    if datayear == 2016:
        print "2016 selected as dataset year"
        selector = lambda : PhysicsObjects( year=2016 )
    elif datayear == 2017:
        print "2017 selected as dataset year"
        selector = lambda : PhysicsObjects( year=2017 )
    elif datayear == 2018:
        print "2018 selected as dataset year"
        selector = lambda : PhysicsObjects( year=2018 )
    else:
        selector = lambda : PhysicsObjects( year=2018 )
    
    return selector

