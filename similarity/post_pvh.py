#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 18 07:57:31 2022
post analysis script for pv+h2o system
sample sum file
#  id    H           O           Mg           Si    chi
#    sol l  w
#    1  2   3     4   5   6     7   8  9      10  11  12     13
0    0 191 25    301 391 496    90 89 181    108 102 150    0.0094


@author: jiedeng
"""
import matplotlib.pyplot as plt

def show_atoms_num(ch):
    fig,ax = plt.subplots(1,3,figsize=(15,4),sharey=True)
    ax[0].set_title('solid')
    ax[0].plot(ch[:,0],ch[:,7],label='Mg = {0:.1f}'.format(ch[:,7].mean()))
    ax[0].plot(ch[:,0],ch[:,10],label='Si = {0:.1f}'.format(ch[:,10].mean()))
    ax[0].plot(ch[:,0],ch[:,4],label='O = {0:.1f}'.format(ch[:,4].mean()))
    ax[0].plot(ch[:,0],ch[:,1],label='H = {0:.1f}'.format(ch[:,1].mean()))
    ax[0].set_xlabel('step');ax[0].set_ylabel('# of atoms')
    ax[0].legend()
    ax[1].set_title('liquid')
    ax[1].plot(ch[:,0],ch[:,8],label='Mg = {0:.1f}'.format(ch[:,8].mean()))
    ax[1].plot(ch[:,0],ch[:,11],label='Si =  {0:.1f}'.format(ch[:,11].mean()))
    ax[1].plot(ch[:,0],ch[:,5],label='O = {0:.1f}'.format(ch[:,5].mean()))
    ax[1].plot(ch[:,0],ch[:,2],label='H = {0:.1f}'.format(ch[:,2].mean()))
    ax[1].legend()
    ax[1].set_xlabel('step');
    ax[2].set_title('interface')
    ax[2].plot(ch[:,0],ch[:,9],label='Mg= {0:.1f}'.format(ch[:,9].mean()))
    ax[2].plot(ch[:,0],ch[:,12],label='Si= {0:.1f}'.format(ch[:,12].mean()))
    ax[2].plot(ch[:,0],ch[:,6],label='O= {0:.1f}'.format(ch[:,6].mean()))
    ax[2].plot(ch[:,0],ch[:,3],label='H = {0:.1f}'.format(ch[:,3].mean()))
    ax[2].legend()
    ax[2].set_xlabel('step');
    fig.savefig('atoms.pdf',bbox_inches='tight')

    

def show_water_content(ch):
    xs=ch[:,1]/2/(ch[:,10]+ch[:,1]/2)  # water molal fraction mol of water/(mole of H/2 + mol of Si)
    xl=ch[:,2]/2/(ch[:,11]+ch[:,2]/2)
    xw=ch[:,3]/2/(ch[:,12]+ch[:,3]/2)
    fig,ax = plt.subplots(1,1,figsize=(5,4))
    ax.plot(ch[:,0],xs,label='s = {0:.5f}'.format(xs.mean()))
    ax.plot(ch[:,0],xl,label='l = {0:.5f}'.format(xl.mean()))
    ax.plot(ch[:,0],xw,label='4w = {0:.5f}'.format(xw.mean()))
    ax.set_xlabel('step');ax.set_ylabel('water content' )
    ax.legend()
    fig.savefig('water.pdf',bbox_inches='tight')
    return xs, xl, xw

def show(ch,beg=0,end=-1):
    ch = ch[beg:end]
    print('++'*20)
    print(' '*10 + 'Do select frames based on chi, and nonzero_sl' + ' '*10)
    print('++'*20)
    print(' '*10 + 'Step 1: # of atoms  in each phase vs. step' + ' '*10)
    show_atoms_num(ch)
    
    # print(' '*10 + 'Step 2: z length, volume, and density of each phase' + ' '*10)
    # vs, vl, vw, rhos, rhol, rhow = show_length_scale(ch)
    
    print(' '*10 + 'Step 3: water content in each phase (H/2)/(Si+H/2)' + ' '*10)
    xs, xl, xw = show_water_content(ch)
    
    # print(' '*10 + 'Step 4: water partitioning vs. step' + ' '*10)
    # w2step = show_water_step(ch,xs/xl)
    # print('++'*20)
    
    print(' '*10 + 'Step 5: Print out statistics' + ' '*10)
    
    # print('solid: 2Mg+4Si+1H = {0}, 2O = {1}'.format(ch[:,16].mean()*2 + ch[:,19].mean()*4 + ch[:,7].mean(), 2*ch[:,22].mean()))
    # print('liquid: 2Mg+4Si+1H = {0}, 2O = {1}'.format(ch[:,17].mean()*2 + ch[:,20].mean()*4 + ch[:,8].mean(), 2*ch[:,23].mean()))
    # print('interface: 2Mg+4Si+1H = {0}, 2O = {1}'.format(ch[:,18].mean()*2 + ch[:,21].mean()*4 + ch[:,9].mean(), 2*ch[:,24].mean()))

    # print("solid (A^3):", vs.mean())
    # print("liquid (A^3):", vl.mean())
    # print("interface (A^3):", vw.mean())

    print("solid xs:, water content(ppm)", xs.mean(), xs.mean()*18/(xs.mean()*18+100)*1e6 ) # add block average with errorbar
    print("liquid xs:, liquid content(ppm)", xl.mean(), xl.mean()*18/(xl.mean()*18+100)*1e6) # add block average with errorbar
    # print("D", xs.mean()/xl.mean(), w2step[-1])


