
import numpy as np
import sympy as sp

class Adinkra:

    def __init__(obj,matr):
        #The constructor for the Adinkra class.  Takes an adjacency matrix as input.
        obj.ADJMatrix = matr
            #this might not need to be saved, as all the info is held in the LMatrices array
        obj.open      = len(matr)
        obj.closed    = len(matr[0])
        obj.nodes     = obj.open+obj.closed
        obj.colors    = abs(np.max(matr)) 
            #This won't work if all of the last color is negative (all of the connections of the highest-indexed color are dashed). 
            #Does this ever happen?  I can't find a method to take the absolute value of a matrix.  It could be done manually if necessary.
        obj.LMatrices = np.zeros((obj.colors, obj.open, obj.closed)) 
        obj.hand      = obj.colors%2
        
        for i in range(0, obj.colors):
            obj.LMatrices[i] = Adinkra.LColBlock(obj, i+1)
        
    def LColBlock(a,color):
        #Calculates the adjacency matrices for each of the color-induced subgraphs of the adinkra.
        h    = a.open;
        w    = a.closed;
        lcol = np.zeros((h, w));
        
        for i in range(0, h):
            for j in range(0, w):
                
                if (a.ADJMatrix[i,j]==color):
                    lcol[i,j] = 1
                elif (a.ADJMatrix[i,j]==color*-1):
                    lcol[i,j] = -1
                    
        return lcol
    
    def liftMatrix(a, m, w):
            #m is the coefficient, w is a list of word parameters.
        d = a.closed
        M = sp.Matrix(np.identity(d, dtype=int))
        for i in w:
            bin = i
            for j in range(d-1, -1, -1):
                if (bin >= (2**j)):
                    M[j,j] *= m
                    bin -= 2**j
        return M
    
    def LTilde(a, color, MLb, MLf):
            #Here, color is the actual number in the adjacency matrix.  We access LMatrices[color-1] because the matrix indeces are off by one.
            #MLb and MLf are the already-calculated node-lifting operators.
        L = a.LMatrices[color-1]
        L = MLb*L*MLf
        return L
    
    def RTilde(a, color, MRf, MRb):
        R = np.transpose(a.LMatrices[color-1])
        R = MRf*R*MRb
        return R
    
    def BanchoffL(a, MLb, MLf, MRf, MRb):
        I=a.colors
        BL = a.LTilde(I, MLb, MLf)
        for i in range(I-1,0,-1):
            if (i%2==1):
                BL = BL*(a.RTilde(i, MRf, MRb))
            else:
                BL = BL*(a.LTilde(i, MLb, MLf))
        return BL
    
    def BanchoffR(a, MLb, MLf, MRf, MRb):
        I=a.colors
        BR = a.RTilde(I, MRf, MRb)
        for i in range(I-1,0,-1):
            if (i%2==1):
                BR = BR*(a.LTilde(i, MLb, MLf))
            else:
                BR = BR*(a.RTilde(i, MRf, MRb))
        return BR
    
    def hymns(a,wb,wf):
        #initialize constants
        mb   = sp.Symbol('mb')
        rhoB = sp.Symbol('rhoB')
        mub  = mb/rhoB
        muf  = sp.Symbol('muf')
        rhoF = sp.Symbol('rhoF')
        mf   = muf/rhoF
        
        #calculate lifting matrices
        MLb = a.liftMatrix(mb,    wb)
        MLf = a.liftMatrix(1/mf,  wf)
        MRf = a.liftMatrix(muf,   wf)
        MRb = a.liftMatrix(1/mub, wb)
        
        #calculate Banchoff matrices
        BL = a.BanchoffL(MLb, MLf, MRf, MRb)
        BR = a.BanchoffR(MLb, MLf, MRf, MRb)
        
        #get eigenvalues
        if (a.hand==0):
            return BL.eigenvals(), BR.eigenvals() #if adinkra has even # of colors
        elif (a.hand==1):
            return (BL*BR).eigenvals(), (BR*BL).eigenvals() #if adinkra has odd # of colors
        
matr = np.genfromtxt('4x8adk-fig8.csv', dtype='int', delimiter=',')
a=Adinkra(matr)
wb=[254,128] #11111110 and 10000000 -- 8th boson is lifted twice, bosons 2-7 are lifted once.
wf=[240]     #11110000 -- fermions 4-8 are lifted once.
print(a.hymns(wb,wf))
    
