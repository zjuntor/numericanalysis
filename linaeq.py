# -*- coding: utf-8 -*-

#Copyright (c) 2013 Juntao ZH<zjuntor@gmail.com>

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in
#all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#THE SOFTWARE.


#This program, as <numerical analysis> homework, is used to solve linear equation group, using Gauss-Seidel method.




#using scipy & numpy lib to solve matrix operation


from scipy import linalg
import numpy as np
 



#@parm : a,b are the matrix and vector in linear equation group ax=b
#@parm : tolerace using the norm of the difference between two successive iterations as tolerace.
#@return : return the root vector of the linear equation group, or null for this method being not convergent. 

def gaussSiedel(a,b,tolerace):
    size = len(b)
   
    #to judge it is convergent or not using Gauss-Seidel method.
    L = np.asmatrix(np.array([([0.0]*size)]*size))
    U = np.asmatrix(np.array([([0.0]*size)]*size))
    D = np.asmatrix(np.array([([0.0]*size)]*size))
    for i in range(0,size):
        for j in range(0,size):
            if i == j:
                D[i,j] = a[i,j]
            if i < j:
                U[i,j] = a[i,j]
            if i > j:
                L[i,j] = a[i,j]
    G = -linalg.inv(D+L).dot(U)
    eigs,eigvectors = linalg.eig(G)
    
    r = max(map(abs,eigs))
    if r >= 1:
        print 'is not convergent using Gauss-Seidel method, Try Jacobi please.'    
        return null

    # begin to iterate to solve the equation group.
    x = np.array([0.0]*size)
    xpre = np.array([0.0]*size)
    
    norm = 10000000
    while norm > tolerace:
        for i in range(0,size):
            sum1 = 0.0
            
            for j in range(0,size):
                sum1 = sum1 + a[i,j]*x[j]
            
            xpre[i] = x[i]
            x[i] = x[i] + (b[i] - sum1)/a[i,i]
            
        norm = max(map(abs,x-xpre))
    return x
