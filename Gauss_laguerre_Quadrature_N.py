# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 23:19:30 2021

@author: DELL
"""
import scipy.integrate as integral
import numpy as np
import numpy.polynomial.laguerre as lag
import math
import matplotlib.pyplot as plt
#Recurrence for laguerre
def L(n,x):
    if n==0:
        return 1
    elif n==1:
        return 1-x
    else:
        return ((2*n-1-x)*L(n-1,x)-(n-1)*L(n-2,x))/(n)
    
#question: f(x)=(2+x^2)^(-1) ; hence
        #  g(x)=exp(x)*(2+x^2)^(-1)
def f(x):
    f=1/(2+x**2)
    return f

def func(x):
    g=np.exp(x)*f(x)
    return g

# Recurrence for L'(n,x):x*L'(n,x)=(x-n-1)*L(n,x)+(n+1)*L(n+1,x)
dl_dx=lambda n,x:((x-n-1)*L(n,x)+(n+1)*L(n+1,x))/(x)


n=int(input("enter n"))

"""finding roots of nth laguerre polynomial and storing it in t"""
Q=[]
for i in range (n):
    Q.append(0)
Q.append(1)
t=lag.lagroots(Q)
print("x found by quadrature rule",t)

#finding the solution
sol=0
for i in range(n):
    #weight
    w=-1/(n*dl_dx(n,t[i])*L(n-1,t[i]))
    #print(w)
    sol1=w*func(t[i])
    sol=sol+sol1   



     
main_sol=sol
print("the value of integral of the function with n=",n,"is\n",main_sol)
man_sol=integral.quad(f, 0, math.inf)
man_sol=round(man_sol[0],10)
print("manual sol",man_sol)

y = lambda x:1/(2+x**2)
x = np.linspace(0,1000,1)
plt.plot(x,f(x))