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



#This program, as <Numerical Analysisis> homework, used for finding root of unitary equation.






#@parm f : the function object in equation f(x) = 0.
#@parm l : the low boundary of the section the root in.
#@parm h : the high boundary of the section the root in.
#@tolerace : the torerace between the root had found and the exact one.
#@return : the root of the equation.

def bisection(f,l,h,tolerace):
    mid = (l + h) / 2.0;
    while abs(f(mid)) > tolerace:
        if f(mid)*f(l) < 0:
            h = mid
        else:
            l = mid
        mid = (l + h) / 2.0
    return mid



#@parm f : the function object in equation f(x) = 0.
#@parm f1 : the first-order deivative of f.
#@parm x : the initial value of the root.
#@tolerace : the torerace between the root had found and the exact one.
#@return : the root of the equation.

def newtown(f,f1,x,tolerace):
    while abs(f(x)) > tolerace:
        x = x - f(x) / f1(x)
    return x



#@parm f : the function object in equation f(x) = 0.
#@parm a : the first value in iteration series
#@parm b : the second value in iteration series.
#@tolerace : the torerace between the root had found and the exact one.
#@return : the root of the equation.

def secant(f,a,b,tolerace):
    while abs(f(b)) > tolerace:
        tmp = b - (b - a) * f(b) / (f(b) - f(a))
        a = b
        b = tmp
    return b
