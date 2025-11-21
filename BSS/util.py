#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import numpy as np
import scipy.signal as ss
from fractions import Fraction

# Ajuste de comprimento do sinal
def length_match(sig1, sig2):
    comp = sig1.shape[0] - sig2.shape[0]
    if comp >= 0:
        sig2 = np.r_[sig2, np.zeros(abs(comp))]
    else:
        sig1 = np.r_[sig1, np.zeros(abs(comp))]
    return sig1, sig2

# Normalização do sinal
def norm_sig(sig):
    sig = sig/np.max(np.abs(sig))
    return sig

# Mudança de valor da frequencia de amostragem do sinal
def changeFs(sig, old_Fs, new_Fs):
    frac = Fraction(new_Fs, old_Fs)
    dws = frac.denominator
    ups = frac.numerator
    return new_Fs, ss.resample_poly(sig,ups,dws)

def corr(y1,y2):
    m1 = y1.mean()
    m2 = y2.mean()
    m12 = (y1*y2).mean()
    sd1 = y1.std()
    sd2 = y2.std()
    return (m12 - m1*m2) / (sd1 * sd2)

# Solução do problema de permutação via correlação dos envelopes
def perm_solver(W,Y):
    nFreq = Y.shape[0]
    nMic = Y.shape[1]
    Y = np.abs(Y)
    P = np.zeros([nFreq, nMic, nMic])
    P[0,...] += np.eye(nMic)
    perm = np.array([[0,1],[1,0]])
    for f in range(1,nFreq):
        y1 = Y[f-1,0,:]
        y2 = Y[f-1,1,:]
        y1a = Y[f,0,:]
        y2a = Y[f,1,:]
        c11 = corr(y1,y1a)
        c21 = corr(y2,y1a)
        c12 = corr(y1,y2a)
        c22 = corr(y2,y2a)
        if (c11 + c22) < (c12 + c21):
            P[f,...] += perm 
            Y[f,...] = perm @ Y[f,...]
        else:
            P[f,...] += np.eye(nMic)
    return P @ W # troca as linhas de W

# Solução do problema de escala via princípio da mínima distorção (N saídas)
# Simplificação do problema de projection back com imagens das fontes de cada mic (NxM saídas)
def scale_solver(W):
    A = np.linalg.inv(W)
    Adgn = np.diagonal(A,axis1=1,axis2=2)
    Adgn = Adgn[...,np.newaxis] * np.eye(2)
    return Adgn @ W # princípio da mínima distorção

# PCA real ou complexo
def pca(X):
    I = np.zeros([X.shape[0],2,2], dtype='complex') + np.eye(2)
    
    # Centralização
    mX = np.mean(X, axis=2)
    mX = mX[..., np.newaxis]
    X = X - mX
    
    # Covariância
    covX = X @ np.conjugate(X.transpose(0,2,1)) / X.shape[2]
    
    # Tranformação
    eig = np.linalg.eig(covX)
    D = eig[0][..., np.newaxis] * I
    E = eig[1]
    V = np.linalg.inv(np.sqrt(D)) @ np.conjugate(E.transpose(0,2,1)) # transformação whitening V = D^(-1/2).Eh
    Z = V @ X # dados transformados pelo PCA
    return V, Z

class util():
    def __init__(self, X):
        self.X = X # misturas
        self.nTime = X.shape[1] # comprimento
        self.nSig = 2 # numero de sinais
        
    def stft(self, f_s=1, wind='hann', tam_jan=1024, sobrep_jan=512, pad=True, psd=False):
        self.tfmix = [] # tensor de misturas em T-F
        self.f = [] # Frequencias
        self.tfr = [] # Quadros
        
        #STFT para cada mistura
        for i in range(self.nSig):
            f, tfr, Zxx = (ss.stft(self.X[i,:], fs=f_s, window=wind, nperseg=tam_jan, 
                                 noverlap=sobrep_jan, padded=pad))
            if psd == True:
                Zxx = np.abs(Zxx)**2
            self.f.append(f)
            self.tfr.append(tfr)
            self.tfmix.append(Zxx)    
        self.f = np.asarray(self.f)
        self.t = np.asarray(self.tfr)
        self.tfmix = np.asarray(self.tfmix)
    
    def istft(self, Xedt, f_s=1, wind='hann', tam_jan=1024, sobrep_jan=512):
        self.Xrec = []
        self.t = []
        for i in range(self.nSig):
            t, sig = ss.istft(Xedt[i], fs=f_s, window=wind, nperseg=tam_jan, noverlap=sobrep_jan)
            sig = sig/max(abs(sig))
            self.Xrec.append(sig)
            self.t.append(t)
        self.Xrec = np.asarray(self.Xrec)

