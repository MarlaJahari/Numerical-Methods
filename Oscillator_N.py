# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 12:01:20 2021

@author: DELL
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 11:13:46 2021

@author: DELL
"""
import sympy as sp
import numpy as np
from sympy.physics.qho_1d import E_n,psi_n
from sympy.abc import y

omega=1
m=1
hb= sp.Symbol("hb")
x=hb/2


n=1

eva = E_n(n,omega)
wf=psi_n(n,y,m,omega)

print(eva)
print(wf)

E= eva.subs({hb:1})
print(E/hbar)
