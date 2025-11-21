#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import numpy as np
import util

class IVAfast:
    def __init__(self, X, nIter=100, eps=1e-10, scale=True):
        self.X = X.transpose(1,0,2)
        self.V, self.Z = util.pca(self.X) 
        self.nIter = nIter
        self.eps = eps
        self.scale = scale
        self.nMic = self.X.shape[1]
        self.nFreq = self.X.shape[0]
        self.nFrames = self.X.shape[2]
        self.nSrc = self.nMic
        self.U = np.tile(np.eye(self.nMic, dtype='complex128'), (self.nFreq,1,1))
        self.Jhist = []
        self.Uhist = []
        
    def separation(self):
        self.Y = self.W @ self.X
    
    def separation1(self):
        self.Y = self.U @ self.Z
    
    def ortogonalizar(self):
        # A operação sqrt(inv(W @ Wh)) deve ser feita por decomposição em autovalores e autovetores
        Uh = np.conjugate(self.U.transpose(0,2,1))
        D, E = np.linalg.eig(self.U@Uh)
        Eh = np.conjugate(E.transpose(0,2,1))
        Dsqrt = np.sqrt(np.linalg.inv(D[:,np.newaxis] * np.tile(np.eye(2),(self.nFreq,1,1))))
        self.U = ((E@Dsqrt)@Eh)@self.U
         
    def projection(self):
        self.W = self.U @ self.V
    
    def cost_func(self):
        self.Yabs = np.abs(self.Y)**2
        # norma do vetor (ao longo das frequências) ou seja, distribuição SSL
        self.Gy = (np.sum(self.Yabs, axis=0)) ** 0.5
        # derivada da função contraste
        self.gy = 1 / (2 * self.Gy[np.newaxis,...] + self.eps)
        # derivada segunda da função contraste
        term1 = (np.sum(self.Yabs,axis=0)**-1.5)[np.newaxis,...]
        self.gy2 = -0.25 * term1
        # função de custo
        self.J = np.sum(self.Gy)
        # histórico
        self.Jhist = np.append(self.Jhist, self.J)

    def update1(self, lin):
        ucol = self.U[:,lin,np.newaxis,:].conj().transpose(0,2,1)
        e1 = np.mean(self.gy[:,lin,np.newaxis,:] + self.Yabs[:,lin,np.newaxis,:] * self.gy2[:,lin,np.newaxis,:], axis=2, keepdims=True) * ucol
        e2 = np.mean(self.Y[:,lin,np.newaxis,:].conj() * self.gy[:,lin,np.newaxis,:] * self.Z, axis=2, keepdims=True)
        ucol = e1 - e2
        norm = np.linalg.norm(ucol,axis=1,keepdims=True)
        ucol = ucol / (norm + self.eps)
        self.U[:,lin,np.newaxis,:] = ucol.conj().transpose(0,2,1)

    def scaling(self):
        self.W = util.scale_solver(self.W)
    
    def execute(self):
        # Separação
        self.separation1()
        self.cost_func()
        self.Uhist.append(self.U)
        for it in range(1,self.nIter):
            self.update1(0)
            self.update1(1)
            self.ortogonalizar()
            self.separation1()
            self.cost_func()
            self.Uhist.append(self.U)
            
        # projeção do PCA
        self.projection()
        
        # resolve ambbiguidade de escala
        if self.scale == True:
            self.scaling()
        
        # Separação final
        self.separation()
        self.Y = self.Y.transpose(1,0,2)
        
        
        
        