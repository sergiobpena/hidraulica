import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from curvas_bombas import Tramo,Tramo_composto
ccgrundfos=[
    (0,383.8),
    (0.508,338.5),
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
class Bomba:
    def __init__(self,curva):
        self.q=None
        self.h=None
        self.lee_puntos(curva)
    def lee_puntos(self,puntos):
        q=[]
        h=[]
        for f in puntos:
            q.append(f[0])
            h.append(f[1])
        self.q=np.array(q)
        self.h=np.array(h)
        self.curva=interpolate.UnivariateSpline(self.q,self.h)
    def get_h(self,q):
        h_new = interpolate.splev(q, self.curva, der=0)
        return h_new
    def get_q(self,h):
        aux=self.h - h
        splaux=interpolate.UnivariateSpline(self.q,aux,s=0)
        sol=splaux.roots()
        if len(sol)>0 and len(sol)<2 :
            return sol[0]
    def reducida(self,tramo):
        self.h_reducida=self.h - tramo.eval_q_np_array(self.q/1000)
        self.curva_reducida=interpolate.UnivariateSpline(self.q,self.h_reducida,s=0)
        return self.curva_reducida
    def get_h_reducida(self,q):
        h_new=interpolate.splev(q,self.curva_reducida,der=0)
        return h_new
    def get_q_reducida(self,h):
        aux=self.h_reducida-h
        splaux=interpolate.UnivariateSpline(self.q,aux,s=0)
        sol=splaux.roots()
        return sol[0]
    def get_x_interpolada(self,spline,y):
        pass
class BombaParalela:
    def __init__(self,bombas):
        self.h=None
        self.q=None
        self.curva=None
        self.bombas=[]
        self.add(bombas)
    def add(self,bombas):
        for b in bombas:
            self.bombas.add(bombas)
        self.bomba_baixa=None
        self.bomba_alta=None
        h_max=0
        h0=0
        h_min=0
        for b in self.bombas:
            if h_max == 0 :
                h_max=b.h_reducida[0]
                self.bomba_alta=b
            else:
                if h_max < b.h[0]:
                    h_max = b.h_reducida[0]
                    self.bomba_alta=b
            if h_min == 0 :
                h_min=np.amin(b.h_reducida)
                self.bomba_baixa=b
            else:
                if h_min > np.amin(b.h_reducida):
                    h_min=np.amin(b.h_reducida)
                    self.bomba_baixa=b
        h0=self.bomba_baixa.h_reducida[0]


        eval_h=np.arange(h_min,h0,1.0)
        res_q1 = []
        # curva_baixa = self.bomba_baixa.curva_reducida
        # curva_alta = self.bomba_alta.curva_reducida

        for i in eval_h:
            qh=0
            for b in self.bombas:
                qh+=b.get_get_q_reducida(i)
            res_q1.append(qh)
        self.curva_combinada=interpolate.UnivariateSpline(np.array(res_q1),eval_h,s=0)
    def get_y(o,a):
        pass







if __name__ == '__main__':
    b=Bomba(ccgrundfos)
    bomba_2=Bomba(ccgrundfos)

    tcr_p=interpolate.UnivariateSpline(b.q, b.h, s=0)


    conexion_pozo1 = Tramo(diametro=55.40, lonxitude=728.61, xeometrica=19.364 + (22 - 3.5), c_locais=7.9, u_diam='mm')
    red=b.reducida(conexion_pozo1)
    # pozo 1
    ca=np.arange(0,2,0.1)
    h_red=red(ca)
    # plt.plot(ca,h_red)
    # plt.plot(b.q,b.h)
    # plt.show()
    #pozo 2
    conexion_pozo2 = Tramo(diametro=55.40, lonxitude=293.401, xeometrica=28.933, c_locais=5.4, u_diam='mm')
    red2=bomba_2.reducida(conexion_pozo2)
    h_red2=red2(ca)

    #impulsion

    imp_pn10 = Tramo(diametro=79.20, lonxitude=1096.10, xeometrica=50.719 + 3.5, c_locais=4.7, u_diam='mm')
    imp_pn16 = Tramo(diametro=73.60, lonxitude=850, xeometrica=33.162, c_locais=2.4, u_diam='mm')
    imp_pn25 = Tramo(diametro=65.40, lonxitude=590, xeometrica=47.685, c_locais=6, u_diam='mm')
    tc=Tramo_composto([imp_pn10,imp_pn16,imp_pn25])
    ca_tot=np.arange(0,4,0.1)
    h_comp=tc.eval_q_np_array(ca_tot/1000)

    ############### curva composta b1 b2
    h_min=0
    x = np.linspace(0, 2, 100)

    # Note that even in the OO-style, we use `.pyplot.figure` to create the figure.
    fig, ax = plt.subplots()  # Create a figure and an axes.
    ax.plot(ca_tot, h_comp, label='Impulsion')  # Plot some data on the axes.
    ax.plot(ca, h_red, label='Pozo 1')  # Plot more data on the axes...
    ax.plot(ca, h_red2, label='Pozo 2')  # ... and some more.
    ax.set_xlabel('Caudal')  # Add an x-label to the axes.
    ax.set_ylabel('Altura')  # Add a y-label to the axes.
    ax.set_title("Bombeos")  # Add a title to the axes.
    ax.legend()  # Add a legend.
    plt.grid(True)
    plt.minorticks_on()
    plt.show()
    print(fig)

