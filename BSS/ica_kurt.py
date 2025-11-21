#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import numpy as np

class ICA_kurt():
    def __init__(self, X):
        self.X = X # misturas 2 x nTime
        self.nTime = self.X.shape[1] # sinais misturados em cada linha, cada coluna uma amostra
        
        self.w = np.array([1,0]).reshape(2,1) # vetor de separação
        self.w2 = np.array([0,1]).reshape(2,1) # vetor de separação 2 (ortogonal a w)
        self.whist = [] # histórico do vetor de separação
        self.kurtose = 0
        self.kurt_hist = [] # histórico da kurtose
        
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
        
    def execICA(self, nIter, eta):
        for o in range(nIter):
            self.whist.append(self.w.copy())
            self.calcKurt()
            self.kurt_hist.append(self.kurtose)
            self.updatew(eta)
            '''
            if o%10 == 0:
                print('iteração:', o)
            '''
        self.s1ch = self.w.T@self.Z
        self.s2ch = self.w2.T@self.Z
    
    def calcKurt(self):
        self.proj = self.w.T@self.Z # projeção dos dados no vetor w (como é um vetor normalizado e deseja-se o 
                                    # valor da norma de cada dado, então a fórmula completa não é necessária)
            
        self.kurtose = np.mean(self.proj**4) - 3 # curtose no espaço descorrelacionado e normalizado
    
    def updatew(self, eta):
        self.w = self.w + eta*np.sign(self.kurtose)*(self.Z@((self.w.T@self.Z).T)**3)/self.nTime
        self.w = self.w/(np.linalg.norm(self.w)) # normalização do vetor
        self.w2 = (self.w.T@np.array([[0,-1],[1,0]])).reshape([2,1]) # w2 ortogonal a w