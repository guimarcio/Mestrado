#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import numpy as np
import util

class IVAng:
    def __init__(self, X, nIter=100, eta=0.1, alfa=0, eps=1e-10, scale=True):
        self.X = X.transpose(1,0,2)
        self.nIter = nIter
        self.eta = eta
        self.alfa = alfa
        self.eps = eps
        self.scale = scale
        self.nMic = X.shape[0]
        self.nFreq = X.shape[1]
        self.nFrames = X.shape[2]
        self.nSrc = self.nMic
        self.W = np.tile(np.eye(2,dtype='complex128'), (self.nFreq,1,1))
        self.I = np.tile(np.eye(2), (self.nFreq,1,1))
        self.Jhist = []
        self.Whist = []
        
    def separation(self):
        self.Y = self.W @ self.X
    
    def cost_func(self):
        # norma do vetor (ao longo das frequências) ou seja, distribuição SSL
        self.Gy = np.linalg.norm(self.Y, axis=0)
        # derivada da função contraste
        self.gy = self.Y / (2 * self.Gy[np.newaxis,...] + self.eps)
        # função de custo
        self.J = np.sum(self.Gy) - 2 * self.nFrames * sum(np.log(abs(np.linalg.det(self.W) + self.eps)))
        # histórico
        self.Jhist = np.append(self.Jhist, self.J)
    
    def updateW(self):
        update = ( self.gy @ self.Y.conj().transpose(0,2,1) / self.nFrames - self.I ) @ self.W
        self.W = self.W - self.eta * update
        
    def scaling(self):
        self.W = util.scale_solver(self.W)
        
    def execute(self):
        # Separação
        self.separation()
        self.cost_func()
        self.Whist.append(self.W)
        for it in range(1,self.nIter):
            self.updateW()
            self.separation()
            self.cost_func()
            self.Whist.append(self.W)
    
        if self.scale == True:
            self.scaling()
    
        # separação após a resolver ambiguidade de escala
        self.separation()
        self.Y = self.Y.transpose(1,0,2)