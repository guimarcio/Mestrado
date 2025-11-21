#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import numpy as np

class ICA_negen():
    def __init__(self, X, gama=1):
        self.X = X # misturas 2 x nTime
        self.nTime = self.X.shape[1] # sinais misturados em cada linha, cada coluna uma amostra
        self.gama = gama
        
        self.w = np.array([1,0]).reshape(2,1) # vetor de separação
        self.w2 = np.array([0,1]).reshape(2,1) # vetor de separação 2 (ortogonal a w)
        self.whist = [] # histórico do vetor de separação
        self.negen = 0
        self.negen_hist = [] # histórico da negentropia
        
    def preprocessing(self):
        # centralização em 0 para cada mistura
        for j in range(2):
            self.X[j,:] = self.X[j,:] - np.mean(self.X[j,:])
        
        # covariância e transformação whitening
        self.XX = (1/self.nTime)*np.einsum('it,jt->ij',self.X, self.X) # covariância das misturas
        self.eig = np.linalg.eig(self.XX) # decomposição em autovals e autovets
        self.D = np.einsum('i,ij->ij', self.eig[0], np.eye(2)) # Matriz diag. dos autovals
        self.E = np.einsum('ij->ji', self.eig[1]) # Matriz dos autovets transposta
        self.V = (np.linalg.inv(np.sqrt(self.D)))@self.E # transformação whitening
        self.Z = self.V@self.X # dados transformados pelo PCA (descorr. e centralizados em 0)
        
    def execICA(self, nIter):
        for o in range(nIter):
            self.whist.append(self.w.copy())
            self.calcNegen()
            self.negen_hist.append(self.negen)
            self.updatew()
            '''
            if o%10 == 0:
                print('iteração:', o)
            '''
        self.s1ch = self.w.T@self.Z
        self.s2ch = self.w2.T@self.Z
    
    def calcNegen(self):
        # projeção dos dados no vetor w (como é um vetor normalizado e deseja-se o 
        # valor da norma de cada dado, então a fórmula completa não é necessária)
        self.proj = self.w.T@self.Z
        
        # Atualização do gama
        self.gama = (np.mean(np.log10(np.cosh(self.proj))) - 0.162) - self.gama
        #print(self.gama)
        
        # Aproximação da negentropia por momentos não polinomiais log10(cosh(y)) (não cresce tão rápido)
        # A constante é E[G(y)] quando y tem pdf normal \mu=0 e \sigma = 1
        self.negen = (np.mean(np.log10(np.cosh(self.proj))) - 0.162)**2 
        self.negen_hist.append(self.negen)
        
    def updatew(self):
        self.w = self.w + self.gama * self.Z@np.tanh(self.w.T@self.Z).T/self.nTime
        self.w = self.w/(np.linalg.norm(self.w)) # normalização do vetor
        self.w2 = (self.w.T@np.array([[0,-1],[1,0]])).reshape([2,1]) # w2 ortogonal a w