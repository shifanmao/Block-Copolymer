#!/usr/bin/python
#for use of polymem simulations
#adapted from Pete
from copy import copy
import subprocess
import math
import csv
import os

SUBMIT = open('condorsubmit','w')

import csv
with open('epsvec') as f:
    reader = csv.reader(f, delimiter=",")
    epsvec = list(reader)

EPS1=0;
EPSF=102;
G=5;
for testnum in range(EPS1,EPSF):
    EPS=float(epsvec[testnum][0])
    EPS=EPS/G

    print " TESTNUM # ", testnum

    ## Make directory to move things over in
    foldername = 'sdata-%s' %(testnum)
    SUBMIT.write('cp *.m %s/mfiles/. \n' %foldername)
    SUBMIT.write('condor_submit %s/condor* \n\n' %foldername)

    if not os.access('%s' %foldername,os.F_OK):
        os.mkdir('%s' %foldername)
    os.chdir('%s' %foldername)

    fname='condor.scalc-%s' %(testnum)
    CONDOR = open(fname,'w')
    CONDOR.write("Universe = vanilla\n")
    CONDOR.write("Getenv = True \n")
    CONDOR.write("Executable = /usr/local/bin/matlab \n")
    CONDOR.write("Initialdir = /tower12/home/shifan/wlc-phase-03-11-15/%s/mfiles/ \n" %foldername)
    CONDOR.write("input = scalc.m\n")
    CONDOR.write("output = out.out\n")
    CONDOR.write("error = out.err\n")
    CONDOR.write("Log = out.log\n")
    CONDOR.write("Arguments = -nodisplay \n")
    CONDOR.write("+NeverSuspend = True\n")
    CONDOR.write("+NeverPreempt = True\n")
    CONDOR.write("+TimeBeforePreempt = 72\n")
    CONDOR.write("Queue")
    CONDOR.close()
    
    ## Next, write out the param files
    if not os.access('mfiles', os.F_OK):
        os.mkdir('mfiles')
        os.mkdir('Sdata')
    os.chdir('mfiles')

    #fname='param.loop%s' %runname
    fname = 'scalc.m'
    PARAM = open(fname,'w')
    PARAM.write('%------ calculate structure factor with rotational averaes ------ \n')
    PARAM.write('randphase_func(%s)' %(EPS))

    os.chdir('../../')
