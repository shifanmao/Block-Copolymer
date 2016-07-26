#!/usr/bin/python
#for use of polymem simulations
#adapted from Pete
from copy import copy
import subprocess
import math
import csv
import os

lk=50
SUBMIT = open('condorsubmit','w')

testnumv = [6]
cpnumv   = [1]
vfacv    = [2]

simdir = 'polymem-exp-05-31-16/rand-pt-exp-05-31-16'
calcdir = 'polymem-exp-05-31-16/scalcbatch-05-31-16'

for i in range(0,len(testnumv)):
    testnum = testnumv[i]
    cpnum = cpnumv[i]
    vfac = vfacv[i]
    
    print " TESTNUM # ", testnum, " WITH CPNUM # ", cpnum
    
    ## Make directory to move things over in
    foldername = 'sdata-%s-%s' %(testnum,cpnum)
    SUBMIT.write('cp /tower12/home/shifan/%s-%s-%s/chilist %s/Sdata/. \n' %(simdir,testnum,cpnum,foldername))
    SUBMIT.write('cp *.m %s/mfiles/. \n' %foldername)
    SUBMIT.write('condor_submit %s/condor* \n\n' %foldername)
    
    if not os.access('%s' %foldername,os.F_OK):
        os.mkdir('%s' %foldername)
        os.chdir('%s' %foldername)
        
        fname='condor.scalc-%s-%s' %(testnum,cpnum)
        CONDOR = open(fname,'w')
        CONDOR.write("Universe = vanilla\n")
        CONDOR.write("Getenv = True \n")
        CONDOR.write("Executable = /usr/local/bin/matlab \n")
        CONDOR.write("Initialdir = /tower12/home/shifan/%s/%s/mfiles/ \n" %(calcdir,foldername))
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
    PARAM.write('scalcfun(%s,%s,%s,%s)' %(testnum,cpnum,lk,vfac))

    os.chdir('../../')
