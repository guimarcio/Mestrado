#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import numpy as np
import util

class fdicaMLaux:
    def __init__(self, X, nIter = 100, eps = 1e-10, scale = True, perm = True, alfa=0.1):
        self.X = X # shape (nMic,nFreq,nFrames)
        self.nIter = nIter
        self.eps = eps
        self.alfa = alfa
        self.scale = scale
        self.perm = perm
        self.nMic = self.X.shape[0]
        self.nSrc = self.nMic
        self.nFreq = self.X.shape[1]
        self.nFrames = self.X.shape[2]
        self.W = np.zeros([self.nFreq, self.nMic, self.nMic], dtype='complex128') + np.eye(self.nMic)
        self.Whist = []
        
    def cost_func(self):
        logdet = np.log(np.abs(np.linalg.det(self.W)) + self.eps) 
        self.G = np.sqrt( np.abs(self.Y)**2 )  # (n_sources, n_bins, n_frames)
        self.gy = 2 * np.ones_like(self.Y)
        self.J = np.sum(np.mean(self.G, axis=2), axis=1) - 2 * logdet
        self.J = self.J[...,np.newaxis]
    
    def calc_covs(self):
        XX = self.X[:, np.newaxis, :, :] * self.X[np.newaxis, :, :, :].conj()
        XX = XX.transpose(2, 0, 1, 3)  # (n_bins, n_channels, n_channels, n_frames)
        Y_abs = self.G
        denom = (2 * Y_abs + self.eps)
        varphi = self.gy / denom  # (n_sources, n_bins, n_frames) G'r/2r
        # varphi = varphi.transpose(1, 0, 2)  # (n_bins, n_sources, n_frames)
        GXX = varphi[:, :, np.newaxis, np.newaxis, :] * XX[:, np.newaxis, :, :, :]
        self.Vn = np.mean(GXX, axis=-1)  # (n_bins, n_sources, n_channels, n_channels)
        
    def separation(self):
        self.Y = self.W @ self.X.transpose(1,0,2)
    
    def updateW(self):
        U = self.Vn

        E = np.eye(self.nSrc, self.nMic)  # (n_sources, n_channels)
        E = np.tile(E, reps=(self.nFreq, 1, 1))  # (n_bins, n_sources, n_channels)

        for src_idx in range(self.nSrc):
            w_n_Hermite = self.W[:, src_idx, :]  # (n_bins, n_channels)
            U_n = U[:, src_idx, :, :]
            e_n = E[:, src_idx, :]  # (n_bins, n_n_channels)

            WU = self.W @ U_n
            w_n = np.linalg.solve(WU, e_n)  # (n_bins, n_channels)
            wUw = w_n[:, np.newaxis, :].conj() @ U_n @ w_n[:, :, np.newaxis]
            wUw = np.real(wUw[..., 0])
            wUw = np.maximum(wUw, 0)
            denom = np.sqrt(wUw)
            denom = denom + self.eps
            w_n_Hermite = w_n.conj() / denom
            self.W[:, src_idx, :] = w_n_Hermite
    
    def permutation(self):
        self.W = util.perm_solver(self.W, self.Y)
    
    def scaling(self):
        self.W = util.scale_solver(self.W)
    
    def execution(self):
        # Separação
        self.separation()
        self.cost_func()
        self.Jhist = self.J
        self.Whist.append(self.W)
        self.calc_covs()
        for it in range(1,self.nIter):
            self.updateW()
            self.separation()
            self.cost_func()
            self.Jhist = np.hstack((self.Jhist, self.J))
            self.Whist.append(self.W)
            self.calc_covs()
            
        # Resolvendo permutação e escala
        if self.perm == True:
            self.permutation()
        
        if self.scale == True:
            self.scaling()
        
        # Separação depois do pós processamento
        self.separation()
        self.Y = np.moveaxis(self.Y, [0,1,2], [1,0,2])

