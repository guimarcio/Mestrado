#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import numpy as np
from scipy.signal import correlate as cc
import os
os.system('cms command')

def length_ajust(s1,s2,x1,x2):
    ns1 = len(s1)
    ns2 = len(s2)
    nx1 = len(x1)
    nx2 = len(x2)
    l = [ns1, ns2, nx1, nx2]
    lmax = max(l)
    s1 = np.append(s1, np.zeros(lmax-ns1))
    s2 = np.append(s2, np.zeros(lmax-ns2))
    x1 = np.append(x1, np.zeros(lmax-nx1))
    x2 = np.append(x2, np.zeros(lmax-nx2))
    return s1,s2,x1,x2

# s1 com y1 e s2 com y2 devem estar com correlação cruzada com kmax=0
def metricas(s1,s2,y1,y2,noise=False,n1=0,n2=0):
    # ajuste de comprimento   
    s1, s2, y1, y2 = length_ajust(s1,s2,y1,y2)
    
    # s_target
    Ps1y1 = np.dot(y1,s1)*s1/np.dot(s1,s1)
    Ps2y2 = np.dot(y2,s2)*s2/np.dot(s2,s2)
    
    # e_interf
    Rss = np.zeros([2,2])
    vecy1 = np.zeros([2,1])
    vecy2 = np.zeros([2,1])
    
    Rss[0,0] = np.dot(s1,s1)
    Rss[0,1] = np.dot(s1,s2)
    Rss[1,0] = Rss[0,1]
    Rss[1,1] = np.dot(s2,s2)
    Rinv = np.linalg.inv(Rss)
    
    vecy1[0] = np.dot(y1,s1)
    vecy1[1] = np.dot(y1,s2)
    vecy2[0] = np.dot(y2,s1)
    vecy2[1] = np.dot(y2,s2)
    
    c1 = Rinv@vecy1
    c2 = Rinv@vecy2
    
    Psy1 = c1[0]*s1 + c1[1]*s2
    Psy2 = c2[0]*s1 + c2[1]*s2
    
    e_inty1 = Psy1 - Ps1y1
    e_inty2 = Psy2 - Ps2y2
    
    # e_noise
    if noise == True:
        Psny1 = Psy1 + np.dot(y1,n1)*n1/np.dot(n1,n1) + np.dot(y1,n2)*n2/np.dot(n2,n2)
        Psny2 = Psy2 + np.dot(y2,n1)*n1/np.dot(n1,n1) + np.dot(y2,n2)*n2/np.dot(n2,n2)
        e_noise1 = Psny1 - Psy1
        e_noise2 = Psny2 - Psy2
    else:
        e_noise1 = np.zeros(y1.shape[0])
        e_noise2 = np.zeros(y2.shape[0])
    
    # e_artif
    if noise == True:
        e_artf1 = y1 - Psny1
        e_artf2 = y2 - Psny2
    else:
        e_artf1 = y1 - Psy1
        e_artf2 = y2 - Psy2
    
    # metricas
    SDRy1 = 10 * np.log10(((np.sum(Ps1y1**2)) / (np.sum((e_inty1+e_noise1+e_artf1)**2)))+1e-10)
    SDRy2 = 10 * np.log10(((np.sum(Ps2y2**2)) / (np.sum((e_inty2+e_noise2+e_artf2)**2)))+1e-10)
    
    SIRy1 = 10 * np.log10(((np.sum(Ps1y1**2)) / (np.sum(e_inty1**2)))+1e-10)
    SIRy2 = 10 * np.log10(((np.sum(Ps2y2**2)) / (np.sum(e_inty2**2)))+1e-10)
    
    SARy1 = 10 * np.log10(((np.sum((Ps1y1+e_inty1+e_noise1)**2)) / ((np.sum(e_artf1**2))) + 1e-10))
    SARy2 = 10 * np.log10(((np.sum((Ps2y2+e_inty2+e_noise2)**2)) / ((np.sum(e_artf1**2))) + 1e-10))
        
    if noise == True:
        SNRy1 = 10 * np.log10(((np.sum((Ps1y1+e_inty1)**2)) / (np.sum(e_noise1**2)))+1e-10)
        SNRy2 = 10 * np.log10(((np.sum((Ps2y2+e_inty2)**2)) / (np.sum(e_noise2**2)))+1e-10)

        
        return SDRy1,SDRy2,SIRy1,SIRy2,SNRy1,SNRy2,SARy1,SARy2
    else:
        return SDRy1,SDRy2,SIRy1,SIRy2,SARy1,SARy2
    
def normalizar(x):
    return x / max(abs(x))

def sisdr(ref,x):
    nr = np.linalg.norm(ref)**2
    num = np.linalg.norm(x@ref*ref / nr)
    den = np.linalg.norm(x@ref*ref / nr - x) + 1e-10
    return 10 * np.log10(num / den)

def time_align(ref,x):
    Nr = ref.shape[0]
    Nx = x.shape[0]
    if Nr > Nx: 
        raise Exception('Sinal de referência maior que sinal em análise.')
    pad = Nx - Nr
    ref_pad = np.append(ref, np.zeros(pad))
    crs = cc(x, ref_pad)
    kadj = x.shape[0] - 1
    max_idx = np.where(abs(crs) == max(abs(crs)))[0][0] - kadj
    k = np.arange(0,len(crs),1) - kadj

    if max_idx < 0:
        x = np.append(np.zeros(abs(max_idx)), x)
    
    if crs[abs(max_idx-kadj)] != max(crs):
        corr = -1
    else:
        corr = 1
    return ref[0:-max_idx], corr*x[max_idx:len(ref)]

def metricasSI(s1,s2,y1,y2):
    # ajuste de comprimento   
    s1, s2, y1, y2 = length_ajust(s1,s2,y1,y2)
    
    # s_target
    Ps1y1 = np.dot(y1,s1)*s1/np.dot(s1,s1)
    Ps2y2 = np.dot(y2,s2)*s2/np.dot(s2,s2)
    
    # e_interf
    Rss = np.zeros([2,2])
    vecy1 = np.zeros([2,1])
    vecy2 = np.zeros([2,1])
    
    Rss[0,0] = np.dot(s1,s1)
    Rss[0,1] = np.dot(s1,s2)
    Rss[1,0] = Rss[0,1]
    Rss[1,1] = np.dot(s2,s2)
    Rinv = np.linalg.inv(Rss)
    
    vecy1[0] = np.dot(y1,s1)
    vecy1[1] = np.dot(y1,s2)
    vecy2[0] = np.dot(y2,s1)
    vecy2[1] = np.dot(y2,s2)
    
    c1 = Rinv@vecy1
    c2 = Rinv@vecy2
    
    Psy1 = c1[0]*s1 + c1[1]*s2
    Psy2 = c2[0]*s1 + c2[1]*s2
    
    e_inty1 = Psy1 - Ps1y1
    e_inty2 = Psy2 - Ps2y2
    
    nr = np.linalg.norm(s1)**2
    num1 = np.linalg.norm(y1@s1*s1 / nr)
    den1 = np.linalg.norm(y1@s1*s1 / nr - y1) + 1e-10
    SISDR1 = 10 * np.log10(num1 / den1)
    
    nr = np.linalg.norm(s2)**2
    num2 = np.linalg.norm(y2@s2*s2 / nr)
    den2 = np.linalg.norm(y2@s2*s2 / nr - y2) + 1e-10
    SISDR2 = 10 * np.log10(num2 / den2)
    
    SISIR1 = 10 * np.log10(((np.sum((s1-y1)**2)) / (np.sum(e_inty1**2)))+1e-10)
    SISIR2 = 10 * np.log10(((np.sum((s2-y2)**2)) / (np.sum(e_inty2**2)))+1e-10)
    
    SISAR1 = -10 * np.log10( (10**(-SISDR1/10) - 10**(-SISIR1/10)) )
    SISAR2 = -10 * np.log10( (10**(-SISDR2/10) - 10**(-SISIR2/10)) )
    
    return SISDR1, SISDR2, SISIR1, SISIR2, SISAR1, SISAR2

def PESQ(fs, ori, est):
    pesq = os.popen('pesq +'+str(fs)+" "+ori+" "+est).read()[-6:-1]
    try:
        pesq = float(pesq)
        return pesq
    except:
        return -5
       