# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 20:40:33 2021

@author: Abhay Singh Rawat
Roll No: 19/17026
Fourier series 

1) Square Wave
"""
print("Square Wave")
import scipy.integrate as integrate
import numpy as np
import matplotlib.pyplot as plt
def f(x):   #function for the series
    f=4
    return f

#xi=float(input("Enter the initial limit"))   #limits of the function
#xf=float(input("Enter the final limit"))
xi=-np.pi
xf=np.pi
n=int(input("Enter the no. of terms you want to be used in the interval"))


A=np.zeros(n)   #an array for the values of a(n) 
B=np.zeros(n)     #an array for the values of b(n)
a0=integrate.quad(f, xi,xf)    #a(0)
a0=a0[0]/(xf-xi)
for i in range(1,n+1):
    fa1=lambda x:np.sin(2*i*np.pi*x/(xf-xi))*f(x)   #part to be integrated for a(n)
    fa2=lambda x:np.sin(2*i*np.pi*x/(xf-xi))*0      #since function has different values in diff intervals
    a1=integrate.quad(fa1,0,xf)                     #integration wrt limits for an
    a2=integrate.quad(fa2,xi,0)
    fb1=lambda x:np.cos(2*i*np.pi*x/(xf-xi))*f(x)   #part to be integrated for b(n)
    fb2=lambda x:np.cos(2*i*np.pi*x/(xf-xi))*0      #since function has different values in diff intervals
    b1=integrate.quad(fb1, 0,xf)                    #integration wrt limits for bn
    b2=integrate.quad(fb2,-xf,0)
    
    A[i-1]=(2/(xf-xi))*(a1[0]+a2[0])  #storing all a(n)
    A[i-1]=round(A[i-1],10) #rouding to 10 digits
    B[i-1]=(2/(xf-xi))*(b1[0]+b2[0]) #storing all b(n)
    B[i-1]=round(B[i-1],10)
#print('a',A,B)

n1=int(input('No. of divisions in the range of x,n1 '))    #no. of divisions in x we want
h=(xf-(-xf))/n1
x=np.arange(-xf,xf,h)   #storing n1 in x
F=np.zeros(n1) #a matrix for the final values of f(x)

for j in range(n1):
    sum1=0
    sum2=0
    for i in range(n):
        sum1=sum1+A[i]*np.sin(2*(i+1)*np.pi*x[j]/(xf-xi))   #sum of all A for a given x
        sum2=sum2+B[i]*np.cos(2*(i+1)*np.pi*x[j]/(xf-xi))   #sum of all B for a given x
    F[j]=sum1+sum2+a0/2                                     #f(x)=a0/2+sum(Ai)+sum(Bi)
#plotting
plt.subplot(1,2,1)
plt.plot(x,F)
plt.title("square wave(0,4)")
plt.grid()



"""Full Wave Rectifier"""
print("*************************\nFull Wave Rectifier")

n=int(input("Enter the no. of terms you want to be used in the interval"))
f2=lambda x:np.sin(x)
fl=lambda x:-np.sin(x)
A1=np.zeros(n)   #an array for the values of a(n) 
B1=np.zeros(n)     #an array for the values of b(n)
a01=(integrate.quad(f2, 0,xf))[0]+(integrate.quad(fl, xi,0))[0]    #a(0)
a01=a01*2/(xf-xi)
for i in range(1,n+1):
    a1=integrate.quad(lambda x:np.sin(2*i*np.pi*x/(xf-xi))*np.sin(x),0,xf)  #integration wrt limits
    a2=integrate.quad(lambda x:np.sin(2*i*np.pi*x/(xf-xi))*(-np.sin(x)),xi,0)
    b1=integrate.quad(lambda x:np.cos(2*i*np.pi*x/(xf-xi))*np.sin(x), 0,xf)
    b2=integrate.quad(lambda x:np.cos(2*i*np.pi*x/(xf-xi))*(-np.sin(x)),-xf,0)
    
    A1[i-1]=(2/(xf-xi))*(a1[0]+a2[0])  #a(n)
    A1[i-1]=round(A1[i-1],10)
    B1[i-1]=(2/(xf-xi))*(b1[0]+b2[0]) #b(n)
    B1[i-1]=round(B1[i-1],10)
#print('a',A,B1)

n1=int(input('No. of divisions in the range of x,n1 '))
h=(xf-(-xf))/n1
x=np.arange(-xf,xf,h)
F=np.zeros(n1)

for j in range(n1):
    sum1=0
    sum2=0
    for i in range(n):
        sum1=sum1+A1[i]*np.sin(2*(i+1)*np.pi*x[j]/(xf-xi))
        sum2=sum2+B1[i]*np.cos(2*(i+1)*np.pi*x[j]/(xf-xi))
    F[j]=sum1+sum2+a01/2
plt.subplot(1,2,2)
plt.plot(x,F)
plt.title("Full Wave rectifier")
plt.grid()
#CrankNicolesonplt.savefig("Fourier",dpi=1000)