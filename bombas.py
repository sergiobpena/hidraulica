#!/usr/bin/env python
# coding: utf-8

# In[83]:
# def add(self,bombas):
#     minima_h_existente=None
#     if len(self.bombas)>0:
#         for b in self.bombas:
#             if minima_h_existente:
#                 if np.amin(b.h) > minima_h_existente:
#                     minima_h_existente=np.amin(b.h)
#     for b in bombas:
#         if minima_h_existente:
#             if np.amin(b.h) < minima_h_existente:
#                 minima_h_existente=np.amin(b.h)
#         else:
#             minima_h_existente=np.amin(b.h)

import numpy as np 
from sympy import symbols,solve , nsolve
from mpmath import mp
from sympy.plotting import plot
#ebara 4n410
# q_curva_ebara_1 = np.array([1.19,1.64,2.03])
# h_curva_ebara_1 = np.array([185,172,157])
#grundfos SP5A-60
ccgrundfos=[
    (0,383.8)
    (0.508,338.5)
    (1.271,259.9),
    (1.406,235.04),
    (1.503,215),
    (1.552,203.5),
    (1.593,193.7),
    (1.642,181),
    (1.692,167.6),
    (1.755,149.7),
    (1.876,112)
]
q_curva_ebara_1 = np.array([1,1.46,1.68])
h_curva_ebara_1 = np.array([295.14,237.85,173.61])

q_colector_1=np.array([2,1.6,0])
h_colector_1=np.array([31.277,27.4,19.64])

q_colector_2=np.array([2,1.6,0])
h_colector_2=np.array([39.9,36.2,28.9])

q_impulsion=np.array([4,3.2,0])
h_impulsion=np.array([169.157,157.756,135.067])

x, y = symbols('x y')
def axuste_cuadrado(array_q,array_h):
    x, y = symbols('x y')
    coeficientes=np.polyfit(array_q,array_h,2,full=True)
    axuste=coeficientes[0][0]*x**2 + coeficientes[0][1]*x + coeficientes[0][2] - y
    return axuste

ce1=axuste_cuadrado(q_curva_ebara_1,h_curva_ebara_1)


cc1=axuste_cuadrado(q_colector_1,h_colector_1)
cc2=axuste_cuadrado(q_colector_2,h_colector_2)
ci=axuste_cuadrado(q_impulsion,h_impulsion)

q_funcionamento=np.array([1.19,1.64,2.03])


def obter_reducida(curva_b,curva_t,qq):
    q=[]
    h=[]
    for n in qq:
        q.append(n)
        s_c_b=solve(curva_b.subs(x,n),y)[0]
        s_c_t=solve(curva_t.subs(x,n),y)[0]
        xx=s_c_b-s_c_t
        h.append(xx)
    a_q=np.array(q,dtype=np.float)
    a_h=np.array(h,dtype=np.float)
    ac=axuste_cuadrado(a_q,a_h)
    return ac
#     return axuste_cuadrado(a_q,a_h)
    
c_bomba_1=obter_reducida(ce1,cc1,q_funcionamento)
c_bomba_2=obter_reducida(ce1,cc2,q_funcionamento)
# t_1=[0.2,0.3,0.4]
# a_1=np.array(t_1)
# print(a_1)
h_combinacion=np.array([130,140,150,160])
qcomb=[]
hcomb=[]
for n in h_combinacion:
    sc1=solve(c_bomba_1.subs(y,n),x)[0]
    sc2=solve(c_bomba_2.subs(y,n),x)[0]
    ss=sc1+sc2
    hcomb.append(n)
    qcomb.append(ss)
mp.dps = 15
c_combinada=axuste_cuadrado(np.array(qcomb,dtype=np.float),np.array(hcomb,dtype=np.float))
conxunto=nsolve((c_combinada,ci),(x,y),(1,6))
print('Funcionamento Instalacion : Q : {} \t H: {}'.format(conxunto[0],conxunto[1]))
f_bomba_1=nsolve((c_bomba_1,ci),(x,y),(1,4))
q_f_b1=solve(c_bomba_1.subs(y,conxunto[1]),x)
h_f_b1=solve(ce1.subs(x,q_f_b1[1]),y)
print('Funcionamento Pozo q : Q : {} \t H: {}'.format(q_f_b1[1],h_f_b1[0]))
q_f_b2=solve(c_bomba_2.subs(y,conxunto[1]),x)
h_f_b2=solve(ce1.subs(x,q_f_b2[1]),y)
print('Funcionamento Pozo q : Q : {} \t H: {}'.format(q_f_b2[1],h_f_b2[0]))
# print('Funcionamento Bomba 1 reduciida : {}'.format(f_bomba_1))
f_bomba_2=nsolve((c_bomba_2,ci),(x,y),(1,4))
# print('Funcionamento Bomba 1 reduciida : {}'.format(f_bomba_2))
gbo = plot(ce1.subs(y,0), show=False,xlim=(1,4),ylim=(0,225),line_color='green')
gim= plot(ci.subs(y,0),show=False,line_color='red')
gbc=plot(c_combinada.subs(y,0),show=False,line_color='yellow')
gbr1=plot(c_bomba_1.subs(y,0),show=False,line_color='brown')
gbr2=plot(c_bomba_2.subs(y,0),show=False,line_color='orange')
gbo.append(gim[0])
gbo.append(gbr1[0])
gbo.append(gbr2[0])
gbo.append(gbc[0])
gbo.show()


# In[ ]:





# In[ ]:





# In[ ]:




