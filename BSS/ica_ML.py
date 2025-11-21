#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import numpy as np

class ICA_ML():
    def __init__(self, X):
        self.X = X # observações T x 2 misturas 
        self.nTime = self.X.shape[1] # sinais misturados em cada linha, cada coluna uma amostra
        self.B = np.eye(2) # matriz de separação
        self.Bhist = [] # histórico do vetor de separação
        self.delBhist = [] # histórico do termo de atualização
        self.Jhist = []
        
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
        #print('Z:'+str(self.Z.shape[0])+'x'+str(self.Z.shape[1]))
        
    def cost(self):
        Sabs = np.abs(self.S)
        detB = self.nTime * np.log(abs(np.linalg.det(self.B)) + 1e-10)
        J = np.sum(Sabs) - detB
        self.Jhist.append(J)      
    
    def execICA(self, nIter, eta):
        for o in range(nIter):
            self.S = self.B@self.Z
            self.cost()
            self.Bhist.append(self.B.copy())
            self.updateB(eta)
            self.delBhist.append(self.delB.copy())
            '''
            if o%10 == 0:
                print('iteração:', o)
            '''
            self.B[0] = self.B[0]/np.linalg.norm(self.B[0])
            self.B[1] = self.B[1]/np.linalg.norm(self.B[1]) 
        self.S = self.B@self.Z
        self.s1ch = self.S[0,:]
        self.s2ch = self.S[1,:]
               
    def updateB(self,eta):
        self.delB = np.linalg.inv(self.B.T) + (1/self.nTime) * -2*np.tanh(self.B@self.Z)@self.Z.T
        self.B = self.B + eta*self.delB 