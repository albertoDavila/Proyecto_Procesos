# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.

from tkinter import Label, Button, Entry, Frame, Tk, messagebox, filedialog
import tkinter as tk
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,NavigationToolbar2Tk)
# from libFracMix import FracMix
from tkinter.ttk import Combobox
import pandas as pd
import numpy as np


METODOS=['Elija una opcion','Minimos cuadrados', 'Minimo cuadrado ponderado']
FORMAS=['Elija una opción','Lote', 'Recursivo']
PROYECTOS=['Elija una opción','Simulación Minimos cuadrados', 'Planta Radioactiva']


class MainFrame(Frame):
    def __init__(self, master=None):
        super().__init__(master, width=600, height=600)
        self.master = master
        self.grid()
        #self.pack()
        self.pantalla_inicio()

        self.Proyecto_Elegido='Nada'

    def Seleccion_proyecto(self):
        if (self.cmbProyectos.get() == PROYECTOS[1]):
            self.create_widgets()
        elif(self.cmbProyectos.get() == PROYECTOS[2]):
            self.Planta_widgets()

    def ElegirOpciones(self):
        if (self.cmbMetodos.get() == METODOS[1]) & (self.cmbFormas.get()== FORMAS[1]):
            self.batch_widgets()
        elif (self.cmbMetodos.get() == METODOS[1]) & (self.cmbFormas.get()== FORMAS[2]):
            self.recursivo_widgets()
        elif (self.cmbMetodos.get() == METODOS[2]) & (self.cmbFormas.get()== FORMAS[1]):
            self.batch_ponderado_widgets()
        elif (self.cmbMetodos.get() == METODOS[2]) & (self.cmbFormas.get() == FORMAS[2]):
            self.recursivo_ponderado_widgets()


            #messagebox.showinfo(title="Entro la funcion", message=" Valio madre")


    def LeerArchivo(self):
        file = filedialog.askopenfilename(title="Seleccionar archivo txt", filetypes=(("Archivos txt", "*.txt"), ("Todos los archivos", "*.*")))
        self.f = pd.read_csv(file, sep ="\t", header=None)

    def plot(self):
        if(self.Proyecto_Elegido=='Minimos Cuadrados'):
            fig = Figure(figsize=(5, 4), dpi=100)
            ax = fig.add_subplot(111)
            canvas_wid = FigureCanvasTkAgg(fig, master=root)
            canvas_wid.get_tk_widget().grid(row=0,column=1)
            canvas_wid.draw()
            self.redraw(canvas_wid,ax)
        elif(self.Proyecto_Elegido=='Planta'):
            fig = Figure(figsize=(5, 8), dpi=100)
            ax = fig.add_subplot(311)
            ax2 = fig.add_subplot(312)
            ax3 = fig.add_subplot(313)
            canvas_wid = FigureCanvasTkAgg(fig, master=root)
            canvas_wid.get_tk_widget().grid(row=0, column=1)
            canvas_wid.draw()
            self.redraw_radioactivo(canvas_wid, ax, ax2, ax3)

    def redraw_radioactivo(self, canvas,ax,ax2, ax3):
        ax.clear()
        #ax.set_xlabel('Tiempo')
        ax.set_ylabel('Valor')
        ax.plot(self.xenon_list)
        ax.legend(['Xenon'])

        ax2.clear()
        ax2.set_ylabel('Valor')
        ax2.plot(self.lodine_list)
        ax2.legend(['Iodine'])

        ax3.clear()
        ax3.set_xlabel('Tiempo')
        ax3.set_ylabel('Valor')
        ax3.plot(self.Num_neutrons_graph_list)
        ax3.legend(['# Neutrons'])
        canvas.draw()



    def redraw(self, canvas,ax):
        ax.clear()
        ax.set_xlabel('Iteraciones')
        ax.set_ylabel('Coeficientes')
        x = []
        for i in range(self.iter + 1):
            x.append(i)
        for i in range(len(self.coeficientes_a) + len(self.coeficientes_b)):
            ax.plot(x, self.theta[i, :self.iter + 1], label= self.entradas[i])
            ax.legend(loc='upper right', frameon=False)
        canvas.draw()



    def pausar(self):
        self.pausa=True
        #if (self.Proyecto_Elegido=='Planta'):
         #   self.Primera=False

    def continuar(self):
        self.pausa = False
        if (self.cmbMetodos.get() == METODOS[1]) & (self.cmbFormas.get() == FORMAS[2]):
            self.minimos_cuadrados_recursivo()
        elif (self.cmbMetodos.get() == METODOS[2]) & (self.cmbFormas.get() == FORMAS[2]):
            self.minimos_cuadrados_recursivo_ponderado()

    def iteracion_Recursivo(self):
        if(self.pausa==False):

            self.iter = self.iter + 1
            self.chi = []
            self.N = self.N + 1
            for i, value in enumerate(self.coef):
                if i < len(self.coeficientes_a):
                    self.chi.append(self.YN[self.N - (value + 1)])

                else:
                    self.chi.append(self.UN[self.N - (value + 1)])

            self.chi = np.matrix(self.chi).T
            self.beta = self.chi.T @ self.P @ self.chi
            self.P = self.P - (self.P @ self.chi * self.chi.T @ self.P) / (1 + self.beta)
            self.t_max = self.t_max - self.P @ self.chi @ (self.chi.T @ self.t_max - self.YN[self.N - 1])
            self.theta[:, self.iter] = self.t_max.ravel()

            self.error = np.asarray(self.YN[self.N-1] - self.chi.T@self.theta[:,self.iter])
            self.error_list.append(self.error)


            if abs(self.error)< abs(next(iter(self.best_error.values()))):
                self.best_error={self.N_estatica +self.iter: self.error}






            self.txtRes.delete(0, 'end')
            self.txtRes.insert(0, self.t_max.ravel())

            self.txtIter.delete(0, 'end')
            self.txtIter.insert(0, self.iter + self.N_estatica )
            self.plot()

            self.txtError.delete(0, 'end')
            self.txtError.insert(0, self.error)

            self.txtBest_itera.delete(0, 'end')
            self.txtBest_itera.insert(0, next(iter(self.best_error.keys())))

            self.txtBest_error.delete(0, 'end')
            self.txtBest_error.insert(0, next(iter(self.best_error.values())))

            root.after(1000, self.iteracion_Recursivo)



    def minimos_cuadrados_recursivo(self):
        self.error_list=[]
        self.best_error={}
        self.pausa = False
        self.iter=0
        self.N = int(self.txtN.get())
        self.N_estatica = self.N
        self.entradas_a = (self.txtn.get().split(" "))
        self.coeficientes_a = list(map(lambda x: x[1], self.entradas_a))
        self.coeficientes_a = list(map(int, self.coeficientes_a))

        self.entradas_b = (self.txtm.get().split(" "))
        self.entradas=self.entradas_a+ self.entradas_b
        self.coeficientes_b = list(map(lambda x: x[1], self.entradas_b))
        self.coeficientes_b = list(map(int, self.coeficientes_b))

        self.psi=np.zeros((self.N,len(self.coeficientes_a)+len(self.coeficientes_b)))

        self.YN = self.f[0]
        self.UN = self.f[1]

        for count, value in enumerate(self.coeficientes_a):
            self.psi[value:, count] = self.YN[0:self.N - value]

        for count, value in enumerate(self.coeficientes_b):
            self.psi[value :, count + len(self.coeficientes_a)] = self.UN[0:self.N - value]

        self.P = np.linalg.inv(self.psi.T @ self.psi)
        self.theta = np.zeros((len(self.coeficientes_a) + len(self.coeficientes_b), 1000))
        self.t = np.matrix((self.P @ self.psi.T @ self.YN[0:self.N]))

        self.t_max = self.t.reshape((self.t.size, 1))
        self.theta[:, self.iter] = self.t_max.ravel()
        #self.coeficientes_a, self.coeficientes_b = zip(*sorted(zip(self.coeficientes_a, self.coeficientes_b)))
        self.coef = self.coeficientes_a + self.coeficientes_b

        self.error = self.YN[self.N - 1] - self.psi[self.N - 1] @ self.theta[:,self.iter]
        self.best_error={self.N_estatica:self.error}
        self.error_list.append(self.error)


        self.txtRes.delete(0, 'end')
        self.txtRes.insert(0, self.t_max.ravel())

        self.txtIter.delete(0, 'end')
        self.txtIter.insert(0, self.iter)

        self.txtError.delete(0, 'end')
        self.txtError.insert(0, self.error)

        self.txtBest_itera.delete(0, 'end')
        self.txtBest_itera.insert(0, next(iter(self.best_error.keys())))

        self.txtBest_error.delete(0, 'end')
        self.txtBest_error.insert(0, next(iter(self.best_error.values())))




        self.btnIteracion = Button(self, text="Pausar", command=self.pausar)
        self.btnIteracion.place(x=250, y=310)

        if (self.pausa == False):
            self.iteracion_Recursivo()





    def minimos_cuadrados_batch(self):
        self.N = int(self.txtN.get())
        self.entradas_a = (self.txtn.get().split(" "))
        self.coeficientes_a = list(map(lambda x: x[1], self.entradas_a))
        self.coeficientes_a = list(map(int, self.coeficientes_a))

        self.entradas_b = (self.txtm.get().split(" "))
        self.coeficientes_b = list(map(lambda x: x[1], self.entradas_b))
        self.coeficientes_b = list(map(int, self.coeficientes_b))

        self.phi=np.zeros((self.N,len(self.coeficientes_a)+len(self.coeficientes_b)))

        self.YN = self.f[0]
        self.UN = self.f[1]

        for count, value in enumerate(self.coeficientes_a):
            self.phi[value:, count] = self.YN[0:self.N - value]

        for count, value in enumerate(self.coeficientes_b):
            self.phi[value :, count + len(self.coeficientes_a)] = self.UN[0:self.N - value]

        self.theta=np.linalg.inv(self.phi.T@self.phi)@self.phi.T@self.YN[:self.N]
        self.error = self.YN[self.N-1] - self.phi[self.N-1] @ self.theta




        #self.plot()

        # messagebox.showinfo(title="Entro la funcion", message=str(theta))
        # result = int(self.txtRes.get())

        self.txtRes.delete(0, 'end')
        self.txtRes.insert(0, self.theta)

        self.txtError.delete(0, 'end')
        self.txtError.insert(0, self.error)



    def minimos_cuadrados_batch_ponderado(self):
        self.N = int(self.txtN.get())
        self.gamma= float(self.txtGamma.get())
        self.wf=float(self.txtwf.get())
        self.entradas_a = (self.txtn.get().split(" "))
        self.coeficientes_a = list(map(lambda x: x[1], self.entradas_a))
        self.coeficientes_a = list(map(int, self.coeficientes_a))

        self.entradas_b = (self.txtm.get().split(" "))
        self.coeficientes_b = list(map(lambda x: x[1], self.entradas_b))
        self.coeficientes_b = list(map(int, self.coeficientes_b))

        self.phi=np.zeros((self.N,len(self.coeficientes_a)+len(self.coeficientes_b)))

        self.YN = self.f[0]
        self.UN = self.f[1]
        self.WN=np.zeros((self.N,self.N), float)

        for i in range(self.N):
            self.WN[i, i] = self.wf * (self.gamma ** (self.N - i - 1))

        for count, value in enumerate(self.coeficientes_a):
            self.phi[value:, count] = self.YN[0:self.N - value]

        for count, value in enumerate(self.coeficientes_b):
            self.phi[value :, count + len(self.coeficientes_a)] = self.UN[0:self.N - value]

        self.theta=np.linalg.inv(self.phi.T@self.WN@self.phi)@self.phi.T@self.WN@self.YN[:self.N]

        self.error = self.YN[self.N-1] - self.phi[self.N-1] @ self.theta




        # messagebox.showinfo(title="Entro la funcion", message=str(theta))
        # result = int(self.txtRes.get())

        self.txtRes.delete(0, 'end')
        self.txtRes.insert(0, self.theta)

        self.txtError.delete(0, 'end')
        self.txtError.insert(0, self.error)

    def batch_ponderado_widgets(self):
        for child in self.winfo_children():
            child.destroy()
        plt.close()

        self.create_widgets()

        self.btnArchivo = Button(self, text="Subir archivo", command = self.LeerArchivo)
        self.btnArchivo.place(x=150, y=180)
        Label(self, text="Cargar Archivo:").place(x=30, y=180)

        Label(self, text="Cantidad de datos (N):").place(x=30, y=220)
        self.txtN = Entry(self, width=15)
        self.txtN.place(x=370, y=220)

        Label(self, text="Valor del ratio de atenuación (Ɣ):").place(x=30, y=250)
        self.txtGamma = Entry(self, width=15)
        self.txtGamma.place(x=370, y=250)

        Label(self, text="Valor del factor de ponderación (a):").place(x=30, y=280)
        self.txtwf = Entry(self, width=15)
        self.txtwf.place(x=370, y=280)

        Label(self, text="Variables a's (empezar por a1 y separar con espacios)").place(x=30, y=310)
        self.txtn = Entry(self, width=15)
        self.txtn.place(x=370, y=310)

        Label(self, text="Variables b's (empezar por b0 y separar con espacios)").place(x=30, y=340)
        self.txtm = Entry(self, width=15)
        self.txtm.place(x=370, y=340)

        self.btnEmpezar = Button(self, text="Empezar", command=self.minimos_cuadrados_batch_ponderado)
        self.btnEmpezar.place(x=150, y=370)

        Label(self, text="Resultado").place(x=30, y=400)
        self.txtRes = Entry(self, width=50)
        self.txtRes.place(x=100, y=400)

        Label(self, text="Error").place(x=30, y=430)
        self.txtError = Entry(self, width=20)
        self.txtError.place(x=100, y=430)




    def recursivo_widgets(self):
        for child in self.winfo_children():
            child.destroy()

        plt.close()

        self.create_widgets()
        self.cmbFormas.set(FORMAS[2])
        self.cmbMetodos.set(METODOS[1])

        self.btnArchivo = Button(self, text="Subir archivo", command = self.LeerArchivo)
        self.btnArchivo.place(x=150, y=180)
        Label(self, text="Cargar Archivo:").place(x=30, y=180)

        Label(self, text="Cantidad de datos inciales (N):").place(x=30, y=220)
        self.txtN = Entry(self, width=15)
        self.txtN.place(x=370, y=220)

        Label(self, text="Variables a's (empezar por a1 y separar con espacios)").place(x=30, y=250)
        self.txtn = Entry(self, width=15)
        self.txtn.place(x=370, y=250)

        Label(self, text="Variables b's (empezar por b0 y separar con espacios)").place(x=30, y=280)
        self.txtm = Entry(self, width=15)
        self.txtm.place(x=370, y=280)

        self.btnEmpezar = Button(self, text="Empezar", command=self.minimos_cuadrados_recursivo)
        self.btnEmpezar.place(x=150, y=310)

        Label(self, text="Iteración actual").place(x=30, y=340)
        self.txtIter = Entry(self, width=5)
        self.txtIter.place(x=130, y=340)

        Label(self, text="Resultado").place(x=30, y=370)
        self.txtRes = Entry(self, width=50)
        self.txtRes.place(x=100, y=370)

        Label(self, text="Error").place(x=30, y=400)
        self.txtError = Entry(self, width=20)
        self.txtError.place(x=100, y=400)

        Label(self, text="Mejor resultado encontrado: ").place(x=30, y=450)
        Label(self, text="Iteracion: ").place(x=30, y=470)
        self.txtBest_itera = Entry(self, width=10)
        self.txtBest_itera.place(x=100, y=470)

        Label(self, text="Error: ").place(x=230, y=470)
        self.txtBest_error = Entry(self, width=15)
        self.txtBest_error.place(x=300, y=470)

    def batch_widgets(self):
        for child in self.winfo_children():
            child.destroy()
        plt.close()

        self.create_widgets()

        self.btnArchivo = Button(self, text="Subir archivo", command = self.LeerArchivo)
        self.btnArchivo.place(x=150, y=180)
        Label(self, text="Cargar Archivo:").place(x=30, y=180)

        Label(self, text="Cantidad de datos (N):").place(x=30, y=220)
        self.txtN = Entry(self, width=15)
        self.txtN.place(x=370, y=220)

        Label(self, text="Variables a's (empezar por a1 y separar con espacios)").place(x=30, y=250)
        self.txtn = Entry(self, width=15)
        self.txtn.place(x=370, y=250)

        Label(self, text="Variables b's (empezar por b0 y separar con espacios)").place(x=30, y=280)
        self.txtm = Entry(self, width=15)
        self.txtm.place(x=370, y=280)

        self.btnEmpezar = Button(self, text="Empezar", command=self.minimos_cuadrados_batch)
        self.btnEmpezar.place(x=150, y=310)

        Label(self, text="Resultado").place(x=30, y=340)
        self.txtRes = Entry(self, width=50)
        self.txtRes.place(x=100, y=340)

        Label(self, text="Error").place(x=30, y=370)
        self.txtError = Entry(self, width=20)
        self.txtError.place(x=100, y=370)

    def minimos_cuadrados_recursivo_ponderado(self):
        self.error_list = []
        self.best_error = {}
        self.pausa = False
        self.iter = 0
        self.N = int(self.txtN.get())
        self.N_estatica = self.N
        self.gamma = float(self.txtGamma.get())
        self.wf = float(self.txtwf.get())
        self.entradas_a = (self.txtn.get().split(" "))
        self.coeficientes_a = list(map(lambda x: x[1], self.entradas_a))
        self.coeficientes_a = list(map(int, self.coeficientes_a))

        self.entradas_b = (self.txtm.get().split(" "))
        self.entradas = self.entradas_a + self.entradas_b
        self.coeficientes_b = list(map(lambda x: x[1], self.entradas_b))
        self.coeficientes_b = list(map(int, self.coeficientes_b))

        self.psi = np.zeros((self.N, len(self.coeficientes_a) + len(self.coeficientes_b)))

        self.YN = self.f[0]
        self.UN = self.f[1]
        self.WN = np.zeros((self.N, self.N), float)

        for i in range(self.N):
            self.WN[i, i] = self.wf * (self.gamma ** (self.N - i - 1))

        for count, value in enumerate(self.coeficientes_a):
            self.psi[value:, count] = self.YN[0:self.N - value]

        for count, value in enumerate(self.coeficientes_b):
            self.psi[value:, count + len(self.coeficientes_a)] = self.UN[0:self.N - value]

        self.P = np.linalg.inv(self.psi.T @ self.WN @ self.psi)
        self.theta = np.zeros((len(self.coeficientes_a) + len(self.coeficientes_b), 1000))

        self.t = np.matrix((self.P @ self.psi.T @ self.WN @ self.YN[0:self.N]))
        self.t_max = self.t.reshape((self.t.size, 1))
        self.theta[:, self.iter] = self.t_max.ravel()
        # self.coeficientes_a, self.coeficientes_b = zip(*sorted(zip(self.coeficientes_a, self.coeficientes_b)))
        self.coef = self.coeficientes_a + self.coeficientes_b

        self.error = self.YN[self.N - 1] - self.psi[self.N - 1] @ self.theta[:,self.iter]
        self.best_error = {self.N_estatica: self.error}


        self.txtRes.delete(0, 'end')
        self.txtRes.insert(0, self.t_max.ravel())

        self.txtIter.delete(0, 'end')
        self.txtIter.insert(0, self.iter)

        self.txtError.delete(0, 'end')
        self.txtError.insert(0, self.error)

        self.txtBest_itera.delete(0, 'end')
        self.txtBest_itera.insert(0, next(iter(self.best_error.keys())))

        self.txtBest_error.delete(0, 'end')
        self.txtBest_error.insert(0, next(iter(self.best_error.values())))

        self.btnIteracion = Button(self, text="Pausar", command=self.pausar)
        self.btnIteracion.place(x=250, y=370)


        if (self.pausa == False):
            self.iteracion_Recursivo_ponderado()

        # messagebox.showinfo(title="Entro la funcion", message=str(theta))
        # result = int(self.txtRes.get())

    def iteracion_Recursivo_ponderado(self):
        if (self.pausa == False):
            self.iter = self.iter + 1
            self.chi = []
            self.N = self.N + 1
            for i, value in enumerate(self.coef):
                if i < len(self.coeficientes_a):
                    self.chi.append(self.YN[self.N - (value + 1)])

                else:
                    self.chi.append(self.UN[self.N - (value + 1)])

            self.chi = np.matrix(self.chi).T
            self.C = (1 / self.wf) + self.chi.T @ (self.P / self.gamma) @ self.chi

            self.L = (self.P @ self.chi) / (self.gamma * self.C)

            self.P = (self.P / self.gamma) - ((self.L@ self.chi.T@self.P) / self.gamma)

            self.t_max = self.t_max - (self.P @ self.chi) @ ((self.chi.T @ self.t_max) - self.YN[self.N - 1])
            self.theta[:, self.iter] = self.t_max.ravel()

            self.error = np.asarray(self.YN[self.iter] - self.chi.T@self.theta[:,self.iter])

            if abs(self.error)< abs(next(iter(self.best_error.values()))):
                self.best_error={self.N_estatica +self.iter: self.error}


            self.txtRes.delete(0, 'end')
            self.txtRes.insert(0, self.t_max.ravel())



            self.txtIter.delete(0, 'end')
            self.txtIter.insert(0, self.iter + self.N_estatica)
            self.plot()

            self.txtError.delete(0, 'end')
            self.txtError.insert(0, self.error)

            self.txtBest_itera.delete(0, 'end')
            self.txtBest_itera.insert(0, next(iter(self.best_error.keys())))

            self.txtBest_error.delete(0, 'end')
            self.txtBest_error.insert(0, next(iter(self.best_error.values())))

            root.after(1000, self.iteracion_Recursivo_ponderado)

    def recursivo_ponderado_widgets(self):
        for child in self.winfo_children():
            child.destroy()
        plt.close()
        self.create_widgets()

        self.cmbFormas.set(FORMAS[2])
        self.cmbMetodos.set(METODOS[2])
        self.btnArchivo = Button(self, text="Subir archivo", command=self.LeerArchivo)
        self.btnArchivo.place(x=150, y=180)
        Label(self, text="Cargar Archivo:").place(x=30, y=180)
        Label(self, text="Cantidad de datos (N):").place(x=30, y=220)
        self.txtN = Entry(self, width=15)
        self.txtN.place(x=370, y=220)
        Label(self, text="Valor del ratio de atenuación (Ɣ):").place(x=30, y=250)
        self.txtGamma = Entry(self, width=15)
        self.txtGamma.place(x=370, y=250)

        Label(self, text="Valor del factor de ponderación (a):").place(x=30, y=280)
        self.txtwf = Entry(self, width=15)
        self.txtwf.place(x=370, y=280)

        Label(self, text="Variables a's (empezar por a1 y separar con espacios)").place(x=30, y=310)
        self.txtn = Entry(self, width=15)
        self.txtn.place(x=370, y=310)

        Label(self, text="Variables b's (empezar por b0 y separar con espacios)").place(x=30, y=340)
        self.txtm = Entry(self, width=15)
        self.txtm.place(x=370, y=340)

        self.btnEmpezar = Button(self, text="Empezar", command=self.minimos_cuadrados_recursivo_ponderado)
        self.btnEmpezar.place(x=150, y=370)

        Label(self, text="Iteración actual").place(x=30, y=400)
        self.txtIter = Entry(self, width=5)
        self.txtIter.place(x=130, y=400)

        Label(self, text="Resultado").place(x=30, y=430)
        self.txtRes = Entry(self, width=50)
        self.txtRes.place(x=100, y=430)

        Label(self, text="Error").place(x=30, y=460)
        self.txtError = Entry(self, width=20)
        self.txtError.place(x=100, y=460)

        Label(self, text="Mejor resultado encontrado: ").place(x=30, y=490)
        Label(self, text="Iteracion: ").place(x=30, y=510)
        self.txtBest_itera = Entry(self, width=10)
        self.txtBest_itera.place(x=100, y=510)

        Label(self, text="Error: ").place(x=230, y=510)
        self.txtBest_error = Entry(self, width=15)
        self.txtBest_error.place(x=300, y=510)




    def create_widgets(self):

        Label(self, text="Metodo:").place(x=30, y=90)
        Label(self, text="Forma:").place(x=30, y=120)
        self.Proyecto_Elegido='Minimos Cuadrados'


        self.btnCalcular = Button(self, text="Configuracion", command=self.ElegirOpciones)
        self.btnCalcular.place(x=100, y=150)



        self.cmbMetodos = Combobox(self, width="20", values=METODOS, state="readonly")
        self.cmbMetodos.place(x=100, y=90)
        self.cmbMetodos.current(0)

        self.cmbFormas = Combobox(self, width="20", values=FORMAS, state="readonly")
        self.cmbFormas.place(x=100, y=120)
        self.cmbFormas.current(0)






    def pantalla_inicio(self):
        Label(self, text="Proyecto:").place(x=30, y=90)

        self.btnCalcular = Button(self, text="Siguiente", command=self.Seleccion_proyecto)
        self.btnCalcular.place(x=100, y=150)

        self.cmbProyectos = Combobox(self, width="20", values=PROYECTOS, state="readonly")
        self.cmbProyectos.place(x=100, y=90)
        self.cmbProyectos.current(0)

    def Planta_widgets(self):
        for child in self.winfo_children():
            child.destroy()
        plt.close()

        self.pantalla_inicio()
        self.Proyecto_Elegido = 'Planta'


        Label(self, text="Valor de Referencia :").place(x=30, y=180)
        self.txtRef = Entry(self, width=15)
        self.txtRef.place(x=300, y=180)

        Label(self, text="Neutrons").place(x=30, y=220)
        self.txtNeu = Entry(self, width=30)
        self.txtNeu.place(x=100, y=220)

        Label(self, text="Iodine").place(x=30, y=250)
        self.txtLod = Entry(self, width=30)
        self.txtLod.place(x=100, y=250)

        Label(self, text="Xenon").place(x=30, y=280)
        self.txtXen = Entry(self, width=30)
        self.txtXen.place(x=100, y=280)

        Label(self, text="Moderation").place(x=30, y=310)
        self.txtMod = Entry(self, width=30)
        self.txtMod.place(x=100, y=310)


        self.btnEmpezar = Button(self, text="Empezar", command=self.Planta_simulacion)
        self.btnEmpezar.place(x=150, y=340)

        self.btnIteracion = Button(self, text="Pausar", command=self.pausar)
        self.btnIteracion.place(x=250, y=340)

        self.btnEmpezar = Button(self, text="Continuar", command=self.Continuar)
        self.btnEmpezar.place(x=350, y=340)

        self.lodine = 0
        self.xenon = 0
        self.Num_neutrons = 1
        self.var_arr = 0
        self.kp = .01
        self.ki = .001
        self.moderation = 0
        self.var_loop_list=[]
        self.lodine_list=[]
        self.xenon_list=[]
        self.Num_neutrons_list=[]
        self.Num_neutrons_graph_list=[]
        self.Num_neutrons_list.append(self.Num_neutrons)
        self.Num_neutrons_graph_list.append(self.Num_neutrons)
        self.pausa = False

    def Continuar(self):

        self.pausa=False
        self.Planta_simulacion()


    def Planta_simulacion(self):


        if (self.pausa == False):


            self.lodine_list.append(self.lodine)
            self.xenon_list.append(self.xenon)

            self.var_loop_list.append(self.var_arr)

            self.lodine = np.float64(self.lodine_list[-1])
            self.xenon = np.float64(self.xenon_list[-1])
            self.Num_neutrons = np.float64(self.Num_neutrons_list[-1])
            self.var_arr=np.float64(self.var_loop_list[-1])

            self.ref = np.float64(self.txtRef.get())
            self.a = (self.ref - self.Num_neutrons)
            self.b = self.var_arr + self.a

            self.var_arr=self.b

            self.c = (self.a * self.kp) + (self.b * self.ki)



            if ((self.c > 0) & (self.c <= .60)):
                self.moderation = self.c
            elif (self.c > .60):
                self.c = .60
                self.moderation = self.c
            elif (self.c < 0):
                self.c = 0
                self.moderation = self.c




            self.d = (self.Num_neutrons + 1) * self.c
            self.e = ((1- self.c) * (self.Num_neutrons + 1)) - (self.Num_neutrons + 1)
            self.f = self.d * 3
            self.h = self.e + self.f

            self.Num_neutrons = self.h
            self.Num_neutrons_graph_list.append(self.Num_neutrons)

            if(self.Num_neutrons>= 1000):
                messagebox.showinfo(title="Advertencia de seguridad",
                                    message='El reactor se apagó por rebasar el limite de referencia')
                self.Planta_widgets()



            self.g = (self.lodine * (.995)) + (self.f * .11)
            self.m = self.lodine * (1-.995)
            self.i = (self.m * .5) + (self.xenon * .999)


            self.n = (1 / self.i) + (1 / self.h)

            self.j=-(1/self.n)




            self.Num_neutrons= self.h + self.j
            self.Num_neutrons_list.append(self.Num_neutrons)
            self.lodine = self.g

            self.xenon = self.j + self.i

            self.plot()

            self.txtNeu.delete(0, 'end')
            self.txtNeu.insert(0, self.Num_neutrons_graph_list[-1])

            self.txtLod.delete(0, 'end')
            self.txtLod.insert(0, self.lodine_list[-1])

            self.txtXen.delete(0, 'end')
            self.txtXen.insert(0, self.xenon_list[-1])

            self.txtMod.delete(0, 'end')
            self.txtMod.insert(0, self.moderation)
            #print('referencia utilizada ', self.ref)


            #print('Diferencia entre # Neutros', self.Num_neutrons_list[-2] - self.Num_neutrons_list[-1])

            if(self.ref>=1000):
                messagebox.showinfo(title="Advertencia de seguridad", message='No se puede reducir la referenica a este valor')
                self.Planta_widgets()


            root.after(1000, self.Planta_simulacion)












root = Tk()
root.title('Proyecto Procesos de Control')
root.geometry("800x700")
app = MainFrame(master=root)

app.mainloop()


