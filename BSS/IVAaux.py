#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import numpy as np
import util

class IVAaux:
    def __init__(self, X, nIter=100, eps=1e-10, scale=True):
        self.X = X.transpose(1,0,2)
        self.XX = (X[:,np.newaxis,...] * X[np.newaxis,...].conj()).transpose(2,0,1,3) 
        self.nMic = X.shape[0]
        self.nSrc = self.nMic
        self.nFreq = X.shape[1]
        self.nFrames = X.shape[2]
        self.nIter = nIter
        self.eps = eps
        self.scale = scale
        self.W = np.tile(np.eye(self.nMic, dtype='complex128'), (self.nFreq,1,1))
        self.e = np.eye(self.nMic)
        self.Jhist = []
        self.Whist = []
        
    def separation(self):
        self.Y = self.W @ self.X
    
    def calc_covs(self):
        self.Gr = np.sum(np.abs(self.Y)**2,axis=0)**0.5 # r = sqrt(sumI||y_j,n||²)
        self.gy = np.ones_like(self.Y) / (2*self.Gr + self.eps)
        aux = self.gy[:,:,np.newaxis,np.newaxis,:] * self.XX[:,np.newaxis,...]
        self.V = np.mean(aux,axis=-1).transpose(1,0,2,3) # shape (nSrc,nFreq,nMic,nMic)
        
    def calc_cost(self):
        logdet = self.nFrames * np.sum(np.log(np.abs(np.linalg.det(self.W)) + self.eps))
        self.J = np.sum(self.Gr) - 2 * logdet
        self.Jhist = np.append(self.Jhist,self.J)
    
    def updatew1(self,src):
        wcol = np.linalg.inv(self.W @ self.V[src,...]) @ (self.e[src].T).reshape([2,1])
        wcol_h = wcol.conj().transpose(0,2,1)
        wcol = wcol / ( (wcol_h @ self.V[src,...] @ wcol) ** 0.5 )
        self.W[:,np.newaxis,src,:] = wcol.conj().transpose(0,2,1)
        
    def scaling(self):
        self.W = util.scale_solver(self.W)
        
    def execute(self):
        for it in range(self.nIter):
            self.separation()
            self.calc_covs()
            self.calc_cost()
            self.Whist.append(self.W)
            for nsrc in range(self.nSrc):
                self.updatew1(nsrc)
                
        if self.scale == True:
            self.scaling()

        # separação após a resolver ambiguidade de escala
        self.separation()
        self.Y = self.Y.transpose(1,0,2)


