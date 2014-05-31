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


#This program, as <numerical analysis> homework, is used to integrate 1-d function by being given lower and upper boundarys and tolerance, using Trapezia, Simpton, Cotes and Romberg method. 




import math


#This Method will be used as a template Decorator, not for client user directly.
#@parm : fun is one n-th divided integration method function of Trapezia or Simpson
#@return : the method to integrate by give f- function to be integrate, a,b- lower and upper boundaris, then t- the given tolerance
def wapper(fun):
    def f(f,a,b,t):
        i = 1
        pre = fun(f,a,b,i)
        cur = fun(f,a,b,i+1)
        while abs(cur-pre)>t :
            i = i+1
            pre = cur
            cur = fun(f,a,b,i+1)     
        return cur
    return f




#This method using Trapezia method is finally converse to trapezia(f,a,b,t) by Decorator @wapper, the arguments f,a,b,t are described in comments of wapper above 
#param f : function to be integrate
#param a,b : lower and upper boundaris
#param n : n-th divided integration
#return : the integration result
@wapper
def trapezia(f,a,b,n):
    n = int(math.pow(2,n))
    h = (b-a)*1.0/n
    res = 0
    for i in range(0,n+1):
        res = res + 2*f(a+i*h)
    res = (res - f(a) - f(b))*h/2
    return res 

#This method using Simpson method is finally converse to simpson(f,a,b,t) by Decorator @wapper, the arguments f,a,b,t are described in comments of wapper above 
#param f : function to be integrate
#param a,b : lower and upper boundaris
#param n : n-th divided integration
#return : the integration result
@wapper
def simpson(f,a,b,n):
    h = (b-a)*1.0/(2*n)
    res = 0
    for i in range(1,n+1):
        res = res + 4 * f(a+(2*i-1)*h) + 2 * f(a+2*i*h)
    res = (res +  f(a) -f(b))*h/3
    return res    


#This method is using Romberg method
#param f : function to be integrate
#param a,b : lower and upper boundaris
#param t : the given tolerance
#return : the integration result
def romberg(f,a,b,t):
    res = {'trapezia':[],'simpson':[],'cotes':[],'romberg':[]}
    res['trapezia'].append((f(a)+f(b))*(b-a)/2.0)
    i = 1
    while(True):
        tmp = 0
        tmp2 = (b-a)/math.pow(2,i)
    
        for j in range(0,int(math.pow(2,i-1))):
            tmp = tmp + f(a + (2*j+1)*tmp2)
    
        res['trapezia'].append(res['trapezia'][i-1]/2.0 + tmp*tmp2)
        res['simpson'].append(res['trapezia'][i]*4/3.0 - res['trapezia'][i-1]/3)
        if i>1:
            res['cotes'].append(res['simpson'][i-1]*16/15.0 - res['simpson'][i-2]/15.0)
        if i>2:
            res['romberg'].append(res['cotes'][i-2]*64/63.0 - res['cotes'][i-3]/63.0)
        if i>3:
            if abs(res['romberg'][-1]-res['romberg'][-2])<t:    
                return res['romberg'][-1]
        i = i + 1
