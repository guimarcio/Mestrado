#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import numpy as np
import util

class fdicaMLnat:
    def __init__(self, X, nIter = 100, eta = 0.1, eps = 1e-10, scale = True, perm = True, alfa=0.1):
        self.X = np.moveaxis(X,[0,1,2], [1,0,2]) # shape (nFreq,nMic,nFrames)
        self.nIter = nIter
        self.eta = eta
        self.eps = eps
        self.alfa = alfa
        self.scale = scale
        self.perm = perm
        self.nMic = self.X.shape[1]
        self.nFreq = self.X.shape[0]
        self.nFrames = self.X.shape[2]
        self.W = np.zeros([self.nFreq, self.nMic, self.nMic], dtype='complex') + np.eye(self.nMic)
        self.I = np.zeros([self.nFreq, self.nMic, self.nMic], dtype='complex') + np.eye(self.nMic)
        self.Whist = []
        
    def contrasts(self):
        # G(y) = -ln( exp(- sqrt(|yn|² + alfa))/b ) : b=1 ; alfa=0.1
        self.Gy = np.sqrt( np.abs(self.Y)**2 + self.alfa )
        # g(y) = G'(y)
        self.gy = self.Y / (2 * self.Gy + self.eps)
    
    def cost_func(self):
        # Função de custo por Max. Likel.
        self.J = self.nFrames*(np.sum(np.mean(self.Gy,axis=2),axis=1) - 
                               2*np.log(np.abs(np.linalg.det(self.W))) + self.eps)
        self.J = self.J[..., np.newaxis] # shape(nFreq,1)
        
    def separation(self):
        self.Y = self.W @ self.X
    
    def updateW(self):
        Yh = np.conjugate(np.moveaxis(self.Y,[0,1,2],[0,2,1])) # hermitiana do vetor (conjugado transposto)
        self.W = self.W - self.eta*((self.gy @ Yh) / self.nFrames - self.I)@self.W
    
    def permutation(self):
        self.W = util.perm_solver(self.W, self.Y)
    
    def scaling(self):
        self.W = util.scale_solver(self.W)
    
    def execute(self):
        # Separação
        self.separation()
        self.contrasts()
        self.cost_func()
        self.Whist.append(self.W)
        self.Jhist = self.J
        for it in range(1,self.nIter):
            self.updateW()
            self.separation()
            self.contrasts()
            self.cost_func()
            self.Whist.append(self.W)
            self.Jhist = np.hstack((self.Jhist, self.J))
            
        # Resolvendo permutação e escala
        if self.perm == True:
            self.permutation()
        
        if self.scale == True:
            self.scaling()
        
        # Separação depois do pós processamento
        self.separation()
        self.Y = np.moveaxis(self.Y, [0,1,2], [1,0,2])

