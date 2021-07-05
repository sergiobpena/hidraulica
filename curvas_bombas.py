import math
from sympy import symbols, log, solve, nsolve
from sympy.plotting import plot as syplt
from scipy import constants
from mpmath import mp
import copy
import numpy as np
# import matplotlib.pyplot as plt

mp.dps = 15
# q, h = symbols('q h')


class Configuracion:
    K = 0.02 * 10 ** -3
    g = constants.g
    visc = 1.01 * 10 ** -6


class Tramo_abstracto:
    def __init__(self):
        self.curva = 0
        self.dhloc = 0
        self.dhlt = 0
        self.dhlu = 0
        self.xeometrica = 0
        self.velocidade = 0
        self.lonxitude = 0

    def get_h(self, caudal):
        h = self.get_curva().subs(Tramo.q, caudal)
        return h

    def get_q(self, altura):
        resultante = self.get_curva() - altura
        caudal = nsolve(resultante, Tramo.q, 0)
        return caudal

    def get_plot(self):
        graf = syplt(self.get_curva(), 0, 6 / 1000, show=False)
        return graf

    def eval_q_np_array(self, np_array):
        ls = []
        for p in np_array:
            ls.append(self.get_h(p))
        return np.array(ls)


class Tramo(Tramo_abstracto):
    q = symbols('q')
    h = symbols('h')
    f = symbols('f')

    def __init__(self, diametro, lonxitude, xeometrica, c_locais, u_caudal='l', u_diam='mm'):
        super().__init__()
        self._parse_diam(diametro, u_diam)
        self.lonxitude = lonxitude
        self.xeometrica = xeometrica
        self.c_locais = c_locais
        self._area()
        self._velocidade()

    def u_caudal(self, u_caudal):
        if u_caudal == 'l':
            self.ud_caudal = 0.001

    def _parse_diam(self, diametro, u_diam):
        if u_diam == 'mm':
            self.diametro = diametro * 0.001
        else:
            self.diametro = diametro

    def _area(self):
        self.area = math.pi * (self.diametro / 2) ** 2

    def _velocidade(self):
        self.velocidade = Tramo.q / self.area

    def reynolds(self, q):
        reynolds = (self.velocidade * self.diametro) / Configuracion.visc
        s = reynolds.subs(Tramo.q, q)
        return s

    def f_c_w(self, q):
        ps = Configuracion.K / (3.71 * self.diametro)
        rei=self.reynolds(q)
        print(ps)
        print(rei)
        ss = 2.51 / (rei * Tramo.f ** 0.5)
        print(ss)
        f_s = (0.25 / ((log(ps + ss)/math.log(10)) ** 2))-Tramo.f
        r=f_s.subs(Tramo.q,q)
        res=nsolve(r,Tramo.f,10**-6)
        return res

    def get_dhlu(self, q):
        vel = self._velocidade()
        vel = vel.subs(Tramo.q, q)
        s = self.get_f_aprox(q) * ((vel ** 2) / (2 * Configuracion.g * self.diametro))
        return s

    def get_f_aprox(self, q):
        reynolds = self.reynolds(q)
        # print('Reynolds : {}'.format(reynolds))
        ss = 5.74 / (reynolds ** 0.9)
        ps = Configuracion.K / (3.71 * self.diametro)
        suma = ss + ps
        divisor=log(suma)/math.log(10)
        # divisor = log(suma)
        x = 0.25 / (divisor ** 2)
        r = x.subs(Tramo.q, q)
        #         self.j=self.f*((self.velocidade**2)/(2*Configuracion.g*self.diametro))
        return r

    def get_dhlt(self, q):
        dhlt = self.get_f_aprox(q) * (
        ((self.velocidade ** 2) / (2 * Configuracion.g * self.diametro))) * self.lonxitude
        res=dhlt.subs(Tramo.q,q)
        return res

    def get_dhloc(self, q):
        dhloc = self.c_locais * (self.velocidade ** 2) / (2 * Configuracion.g)
        res = dhloc.subs(Tramo.q, q)
        return res

    #     def get_curva(self):
    #         self.curva=self.get_dhlt()+self.get_dhloc()+self.xeometrica
    #         return self.curva
    def get_dh(self, caudal):
        h = self.get_dhlt(caudal) + self.get_dhloc(caudal)
        return h

    def get_h(self, caudal):
        h = self.get_dhlt(caudal) + self.get_dhloc(caudal) + self.xeometrica
        return h


#     def get_q(self,altura):
#         resultante=self.get_curva()-altura
#         caudal=nsolve(resultante,Tramo.q,0)
#         return caudal
#     def get_plot(self):
#         graf=syplt(self.get_curva(),0,6/1000,show=False)
#         return graf


class Tramo_composto(Tramo_abstracto):
    def __init__(self, tramos):
        self.tramos = []
        self.add(tramos)

    def add(self, tramos):
        for t in tramos:
            self.tramos.append(t)

    def get_h(self, caudal):
        h = 0
        for t in self.tramos:
            h += t.get_h(caudal)
        return h

    def eval_q_np_array(self, caudais):
        alturas = []
        for c in caudais:
            a = 0
            for t in self.tramos:
                a += t.get_h(c)
            alturas.append(a)
        altos = np.array(alturas)
        return altos

if __name__ == '__main__':

    conexion_pozo1 = Tramo(diametro=55.40, lonxitude=728.61, xeometrica=19.364 + (22 - 3.5), c_locais=7.9, u_diam='mm')
    conexion_pozo2 = Tramo(diametro=55.40, lonxitude=293.401, xeometrica=28.933, c_locais=5.4, u_diam='mm')

    imp_pn10 = Tramo(diametro=79.20, lonxitude=1096.10, xeometrica=50.719 + 3.5, c_locais=4.7, u_diam='mm')
    imp_pn16 = Tramo(diametro=73.60, lonxitude=850, xeometrica=33.162, c_locais=2.4, u_diam='mm')
    imp_pn25 = Tramo(diametro=65.40, lonxitude=590, xeometrica=47.685, c_locais=6, u_diam='mm')
    caudais=np.arange(0,4,0.25)
    caudais_m=caudais/1000
    # print(caudais)
    # print(imp_pn10.get_dh(0.0032))
    tc=Tramo_composto([imp_pn10,imp_pn16,imp_pn25])
    print(imp_pn10.f_c_w(0.0032))
    print(imp_pn10.get_f_aprox(0.0032))
    # print(imp_pn10.area)