# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 19:03:13 2021

@author: Abhay Singh Rawat
19/17026
"""
from scipy.integrate import odeint
import matplotlib.pyplot as plt 


def dy(x,y):
    dy=x+y**2 #diff eq
    return dy


x1=0
x2=1

#stepsize
h=float(input('enter stepsize'))

n=(x2-x1)/h
n=int(n)

x=x1
y0=1
X=[]
X.append(x)
for i in range(n):
    x=x+h
    X.append(x)


y=odeint(dy,y0,X)
print('value of y at',X[20],'\t',y[20])
plt.plot(X,y,'r')
plt.plot(X[2],y[2],'*')
plt.ylim(1,2)
print(y)