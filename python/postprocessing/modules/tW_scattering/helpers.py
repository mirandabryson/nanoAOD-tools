#!/usr/bin/env python                                                                              
import os, sys
import subprocess

from PhysicsTools.NanoAODTools.postprocessing.modules.tW_scattering.ObjectSelection import *

def chooseselector(datayear):
    if datayear == 2016:
        selector = lambda : PhysicsObjects( year=2016 )
    elif datayear == 2017:
        selector = lambda : PhysicsObjects( year=2017 )
    elif datayear == 2018:
        selector = lambda : PhysicsObjects( year=2018 )
    else:
        selector = lambda : PhysicsObjects( year=2018 )
    
    return selector

