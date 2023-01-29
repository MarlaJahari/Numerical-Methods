# -*- coding: utf-8 -*-
"""
Created on Sun Mar  7 19:55:25 2021

@author: Abhay Singh Rawat
Roll No.-19/17026
topic-Gauss Hermite Quadrature
"""
import matplotlib.pyplot as plt
import numpy as np
import numpy.polynomial.hermite as her


s=float(input('enter sigma'))
#defining a function for fcatorial
def fact(n):
    f=1
    for i in range(1,n+1):
        f=f*i
    return f

#Recurrence for hermite
def H(n,x):
    if n==0:
        return 1
    elif n==1:
        return 2*x
    else:
        return (2*x*H(n-1,x)-2*(n-1)*H(n-2,x))

#question: f(x)=exp(-((x-2)^2)/2s^2)*(x+3)
#taking y=(x-2)/(sqrt(2*s)), we get the following function:

#the required function will have y as roots of nth hermite sol
def g(y):
    g=(2*y*(s**2)+(5*s*np.sqrt(2)))/(np.sqrt(2*(np.pi)*s**2))
    return g

    


n=int(input("enter n"))

"""finding roots of nth hermite polynomial and storing it in t"""
Q=[]
for i in range (n):
    Q.append(0)
Q.append(1)
t=her.hermroots(Q)

print("y found by quadrature rule",t)
#finding the solution
sol=0
for i in range(0,n):
    #weight
    w=(2**(n-1))*fact(n)*np.sqrt(np.pi)/((n*H(n-1,t[i]))**2)
    print(w)
    sol1=w*g(t[i])
    #print(sol1)
    sol=sol+sol1   


main_sol=sol
print("the value of integral of the function with n=",n,"is\n",main_sol)

#inbuilt funcm hermgauss(n) returns weights&roots for n sample points    
def inbuilt(f):
   [p,q] = her.hermgauss(n)    #[xi,wi for Hn(x)]
   return sum(q*f(p))      #returning the quadrature sum
print("Inbuilt hermgauss()-",inbuilt(g))

#Plotting
f=lambda x,s:(np.exp(-((x-2)**2)/2*s**2))*(x+3)/(np.sqrt(2*np.pi*s**2))
x=np.linspace(-10,10,100000)

plt.plot(x,f(x,0.01))
plt.xlim(-5,5)
plt.grid()