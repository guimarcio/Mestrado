#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import util

class fdicaMLfast:
    def __init__(self, X, nIter=100, eps=1e-10, scale=True, perm=True, alfa=0.1):
        self.V, self.Z = util.pca(X.transpose(1,0,2)) # X.shape = (nFreq,nMic,nFrames)
        self.X = X.transpose(1,0,2)
        self.nIter = nIter
        self.eps = eps
        self.alfa = alfa
        self.scale = scale
        self.perm = perm
        self.nMic = self.Z.shape[1]
        self.nFreq = self.Z.shape[0]
        self.nFrames = self.Z.shape[2]
        self.U = np.zeros([self.nFreq, self.nMic, self.nMic], dtype='complex') + np.eye(self.nMic)
        self.Rot=np.zeros([self.nFreq,self.nMic,self.nMic], dtype='complex')+np.array([[0,1],[-1,0]])
        Uvec = self.U[..., np.newaxis]
        self.u1 = Uvec[:,0,...] # N = 2 e M = 2 u1.shape = (nfreq,2,1)
        self.Uhist = []
        # U = [u1,...,un]^H
        
    def contrasts(self):
        # G(y) = -ln( exp(- sqrt(|yn|² + alfa))/b ) : b=1 ; alfa=0.1
        self.Gy = np.sqrt( np.abs(self.y1)**2 + self.alfa )
        # g(y) = G'(y)
        self.gy = np.conjugate(self.y1) / (2 * self.Gy + self.eps)
        # g'(y) = G''(y)
        self.gy2 = 1 / (2*self.Gy + self.eps) * ( 1 - 0.5*(self.Gy**2 - self.alfa) / (self.Gy**2 + self.eps) )
        
    def cost_func(self):
        # Função de custo por Max. Likel.
        self.J = self.nFrames*(np.mean(self.Gy,axis=2))
        self.J = self.J[..., np.newaxis] # shape(nFreq,1)
        
    def separation1(self):
        self.y1 = np.conjugate(self.u1).transpose(0,2,1) @ self.Z
        
    def separation(self):
        self.Y = self.W @ self.X
                               
    def updateu1(self):
        e1 = np.mean(self.gy * self.Z, axis=2, keepdims=True)
        e2 = np.mean(self.gy2,axis=2,keepdims=True) * self.u1
        self.u1 = -e1 + e2        
        self.u1 = self.u1 / np.linalg.norm(self.u1,axis=1, keepdims=True)
        
    def permutation(self):
        self.W = util.perm_solver(self.W, self.Y)
    
    def scaling(self):
        self.W = util.scale_solver(self.W)
    
    def save_U(self):
        self.Uhist.append(self.u1)
    
    def projection(self):
        u2 = self.Rot @ self.u1 # cálculo da rotação ortogonal 2D de u1
        self.U = np.conjugate(np.append(self.u1,u2,axis=2).transpose(0,2,1))
        self.W = self.U @ self.V
        self.separation()

    def execution(self):
        # Separação
        self.separation1()
        self.contrasts()
        self.cost_func()
        self.save_U()
        self.Jhist = self.J
        for it in range(1,self.nIter):
            self.updateu1()
            self.separation1()
            self.contrasts()
            self.cost_func()
            self.save_U()
            self.Jhist = np.hstack((self.Jhist, self.J))
        
        # Projeções dos parâmetros e sinais
        self.projection()
                               
        # Resolvendo permutação e escala
        if self.perm == True:
            self.permutation()
        
        if self.scale == True:
            self.scaling()
        
        # Separação depois do pós processamento
        self.separation()
        self.Y = np.moveaxis(self.Y, [0,1,2], [1,0,2])

