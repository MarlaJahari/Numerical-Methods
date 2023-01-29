# -*- coding: utf-8 -*-
"""
"""

import matplotlib.pyplot as plt            
import numpy as np            


#Solving second order BVP using finite difference method

print("GENERAL EXPRESSION")
print("\ny''=p(x)*y'+q(x)*y+F(x)")
p=int(input('Enter p(x)'))
q=int(input('Enter q(x)'))
n=int(input('Enter the no: of steps:'))

#Boundary conditions

a=int(input('What is your a:'))
b=int(input('What is your b:'))
h=(b-a)/n
y1=int(input('Value of y at a:'))
y2=int(input('Value of y at b '))

A=-2*q*h**2-4         
B=h*p+2
C=-h*p+2
x=a                            

X=np.zeros(n+1)
Y=np.zeros(n)
Z=np.zeros(n-1)

for i in range(n-1):
    Z[i]=C
    
for i in range(n):
    Y[i]=A
    
for i in range(n+1):
    X[i]=B

M=np.diag(X,0)+np.diag(Y,1)+np.diag(Z,2)
print(M)

M=np.delete(M,(n-1,n),axis=0) 
M=np.delete(M,(0,n),axis=1)
  
print('Tri-diagonal matrix consisting of the coefficients of y\n',M)

#Implementing Thomas Method

A=[]
B=[]
C=[]
for i in range (n-2):                    #Upper and lower diagonals
    a=M[i+1][i]
    A.append(a)
    c=M[i][i+1]
    C.append(c)
for i in range(n-1):
    b=M[i][i]
    B.append(b)                          #Central diagonal

          
print(A,B,C)                             #Matrices

#Finding values of y at different x

f=x*(x-4)-(h*p+2)*y1 #given function
X=[x]
D=[f]

for i in range(n-1):
    x=x+h
    if i==n-3:
        f=x*(x-4)-2*y2+h*p*y2
        D.append(f)
    else:
        f=x*(x-4)
        D.append(2*f*h**2)
    X.append(x)
    
D=np.delete(D,(0,n))                     #Deleting unnecessary values

print('\n........\n')
print('Constant matrix',D)  
W=np.linalg.solve(M,D)

D[0]=D[0]/B[0]
C[0]=C[0]/B[0]
for i in range(1,n-2):
     C[i]=C[i]/(B[i]-A[i-1]*C[i-1])
for i in range(1,n-1):
     D[i]=(D[i]-A[i-1]*D[i-1])/(B[i]-A[i-1]*C[i-1])
     
y=D[n-2]
Yin=[]
Y=[y1]
Yin.append(y)
for i in range(0,n-2):
    y=D[n-3-i]-y*C[n-3-i]
    Yin.append(y)

for i in range(n-1):
    y=Yin[n-i-2]
    Y.append(y)

Y.append(y2)    

print('\n........\n')
print('Values of x at which y is found',X)
print('\n........\n')
print('Discrete values of y found using Thomas algorithm')
print(Y)
print('\n........\n')

print('Solution using in-built function:',W)
plt.plot(X,Y,label="Solution by thomas algorithm:")

for i in range(5):
    a=(n//4)*i
    plt.plot(X[a],Y[a],'*')
