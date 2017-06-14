#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 13:13:14 2017

@author: dongliu
"""

import numpy
import pylab
import os
import sys
import matplotlib.pyplot as plt


def readexpvalues(expdat,nstep):
    array=numpy.genfromtxt(expdat, delimiter='\t')#pay attention here, every time the type of the data changes, we need to change this line too.
    x=[]    
    for i in range(0,nstep):
        x.append(array[i,1::2].astype('float'))
    return x

def readp(expdat):
    array=numpy.genfromtxt(expdat, delimiter='\t')#pay attention here, every time the type of the data changes, we need to change this line too.
    return array

def readplist(expdat,nstep):
    array=numpy.genfromtxt(expdat, delimiter='\t')#pay attention here, every time the type of the data changes, we need to change this line too.
    x=[]    
    for i in range(0,nstep):
        x.append(array[i,1::1].astype('float'))
    return x    
    
if __name__=="__main__":
    

    nstep=33

    plotrange=100
    
    plotpath="/Users/dongliu/ClionProjects/HFPx2D/cmake-build-debug/"
    #plotpath="/Users/dongliu/Documents/coderesults/Partially filled not convergent"
    
    figp=plt.figure('l_c')
    resultlc=plotpath+"cracklength.txt"
    xlc=readp(resultlc)
    plt.plot(range(0,nstep),xlc,'o')
    plt.show()
    
    figp=plt.figure('iteration')
    resultiteration=plotpath+"iteration.txt"
    xitera=readp(resultiteration)
    plt.plot(range(0,nstep),xitera,'o')
    plt.show()
    
    
    figl=plt.figure('l')
    results=plotpath+"outputcn1.txt"
    x=readexpvalues(results,nstep)
    for i in range(0,nstep):
        plt.plot(x[i])
    plt.show()

    figvl=plt.figure('volume_list')
    resultvl=plotpath+"outputvlist.txt"
    xvl=readplist(resultvl,nstep)
    for i in range(0,nstep):
        plt.plot(xvl[i])  
    plt.show()
    
    
    figel=plt.figure('elas_residual_list')
    resultel=plotpath+"outputelas.txt"
    xel=readexpvalues(resultel,nstep)
    for i in range(0,nstep):
        plt.plot(xel[i])
    plt.show()
#
#    figp=plt.figure('p')
#    resultp=plotpath+"outputpressurecn1.txt"
#    x1=readp(resultp)
#    plt.plot(range(0,nstep+1),x1,'o')
##    plt.yscale('log')
##    plt.xscale('log')
#    plt.show()

    figplist=plt.figure('p_list')
    resultplist=plotpath+"outputplist.txt"
    xplist=readplist(resultplist,nstep+1)
    for i in range(0,nstep+1):
        plt.plot(xplist[i])
    plt.xlim([0,plotrange])
    plt.show()

    
    figcoh=plt.figure('coh')
    resultcoh=plotpath+"outputlcohcn1.txt"
    x2=readp(resultcoh)
    plt.plot(range(0,nstep),x2)
#    plt.yscale('log')
#    plt.xscale('log')
    plt.show()
    
#    figenergy=plt.figure('energy_g J integration over cohesive zone')
#    resultenergy=plotpath+"outputenergyg.txt"
#    xenergy=readp(resultenergy)
#    plt.plot(range(0,nstep),xenergy,'o')
##    plt.yscale('log')
##    plt.xscale('log')
#    plt.show()
#    
#    figenergyp=plt.figure('energy_p G*t from work done by pressure')
#    resultenergyp=plotpath+"outputenergyp.txt"
#    xenergyp=readp(resultenergyp)
#    plt.plot(range(0,nstep),xenergyp,'o')
##    plt.yscale('log')
##    plt.xscale('log')
#    plt.show()
#    
#    figenergyjint=plt.figure('Energy balance')
#    resultenergyjint=plotpath+"outputenergyj.txt"
#    xenergyjint=readp(resultenergyjint)
#    plt.plot(range(0,nstep+1),xenergyjint)
##    plt.yscale('log')
##    plt.xscale('log')
#    plt.show()
##    
#    figenergyf=plt.figure('energy_f G from K calculated by pressure')
#    resultenergyf=plotpath+"outputenergyf.txt"
#    xenergyf=readp(resultenergyf)
#    plt.plot(range(0,nstep),xenergyf)
##    plt.yscale('log')
##    plt.xscale('log')
#    plt.show()
#    
#    figenergycoh=plt.figure('energy_coh G from K calculated by cohesive force')
#    resultenergycoh=plotpath+"outputenergycoh.txt"
#    xenergycoh=readp(resultenergycoh)
#    plt.plot(range(0,nstep),xenergycoh)
##    plt.yscale('log')
##    plt.xscale('log')
#    plt.show()
    
#    figcohf=plt.figure('cohf')
#    resultcohf=plotpath+"outputcohfcn1.txt"
#    xf=readexpvalues(resultcohf,nstep)
##    plt.plot(xf[1])
#    for j in range(0,nstep):
#        plt.plot(xf[j])
#    plt.show()
    
    figstress=plt.figure('stress')
    resultstress=plotpath+"outputstresscn1.txt"
    xstress=readexpvalues(resultstress,nstep)
#    plt.plot(xf[1])
    for m in range(0,nstep):
        plt.plot(xstress[m])
    plt.show()
    
    
    figerrm=plt.figure('errormatrix')
    resulterrm=plotpath+"outputviscoitera.txt"
    errmx=readplist(resulterrm,nstep)
    for i in range(0,nstep):
        plt.plot(errmx[i])
    plt.show()
    
    
#    
#    
#    figv=plt.figure('volume change')
#    resultvol=plotpath+"outputvol.txt"
#    xvol=readp(resultvol)
#    plt.plot(range(0,nstep),xvol,'o')
##    plt.yscale('log')
##    plt.xscale('log')
#    plt.show()