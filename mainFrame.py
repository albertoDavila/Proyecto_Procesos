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


METODOS=['Elija una opcion','Minimos cuadrados', 'otros']
FORMAS=['Elija una opción','Lote', 'Recursivo']



class MainFrame(Frame):
    def __init__(self, master=None):
        super().__init__(master, width=500, height=500)
        self.master = master
        self.grid()
        #self.pack()
        self.create_widgets()




    def ElegirOpciones(self):
        if (self.cmbMetodos.get() == METODOS[1]) & (self.cmbFormas.get()== FORMAS[1]):
            self.batch_widgets()
        elif (self.cmbMetodos.get() == METODOS[1]) & (self.cmbFormas.get()== FORMAS[2]):
            self.recursivo_widgets()
            #messagebox.showinfo(title="Entro la funcion", message=" Valio madre")


    def LeerArchivo(self):
        file = filedialog.askopenfilename(title="Seleccionar archivo txt", filetypes=(("Archivos txt", "*.txt"), ("Todos los archivos", "*.*")))
        self.f = pd.read_csv(file, sep ="\t", header=None)

    def plot(self):

        fig=Figure(figsize=(5,4), dpi=100)

        ax = fig.add_subplot(111)
        canvas_wid = FigureCanvasTkAgg(fig, master=root)
        canvas_wid.get_tk_widget().grid(row=0,column=1)
        canvas_wid.draw()
        self.redraw(canvas_wid,ax)





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
            print(self.theta[i, :self.iter + 1])
        canvas.draw()

    def pausar(self):
        self.pausa=True

    def continuar(self):
        self.pausa = False
        self.iteracion_Recursivo()

    def iteracion_Recursivo(self):
        if(self.pausa==False):
            self.iter = self.iter + 1
            self.chi = []
            self.N = self.N + 1
            for i, value in enumerate(self.coef):
                if i < len(self.coeficientes_a):
                    self.chi.append(self.YN[self.N - (value + 1)])

                else:
                    self.chi.append(self.UN[self.N - (value + 2)])

            self.chi = np.matrix(self.chi).T
            self.beta = self.chi.T @ self.P @ self.chi
            self.P = self.P - (self.P @ self.chi * self.chi.T @ self.P) / (1 + self.beta)
            self.t_max = self.t_max - self.P @ self.chi @ (self.chi.T @ self.t_max - self.YN[self.N - 1])
            self.theta[:, self.iter] = self.t_max.ravel()

            #self.error = self.YN[self.N - 1] - self.chi[self.N - 1] @ self.theta
            self.txtRes.delete(0, 'end')
            self.txtRes.insert(0, self.t_max.ravel())

            self.txtIter.delete(0, 'end')
            self.txtIter.insert(0, self.iter)
            self.plot()

            self.txtError.delete(0, 'end')
            self.txtError.insert(0, self.error)

            root.after(1000, self.iteracion_Recursivo)



    def minimos_cuadrados_recursivo(self):
        self.pausa = False
        self.iter=0
        self.N = int(self.txtN.get())
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
        self.theta = np.zeros((len(self.coeficientes_a) + len(self.coeficientes_b), 50))
        self.t = np.matrix((self.P @ self.psi.T @ self.YN[0:self.N]))

        self.t_max = self.t.reshape((self.t.size, 1))
        self.theta[:, self.iter] = self.t_max.ravel()
        #self.coeficientes_a, self.coeficientes_b = zip(*sorted(zip(self.coeficientes_a, self.coeficientes_b)))
        self.coef = self.coeficientes_a + self.coeficientes_b

        self.error = self.YN[self.N - 1] - self.psi[self.N - 1] @ self.theta

        self.txtRes.delete(0, 'end')
        self.txtRes.insert(0, self.t_max.ravel())

        self.txtIter.delete(0, 'end')
        self.txtIter.insert(0, self.iter)

        self.txtError.delete(0, 'end')
        self.txtError.insert(0, self.error)



        self.btnIteracion = Button(self, text="Continuar", command=self.continuar)
        self.btnIteracion.place(x=100, y=430)
        self.btnIteracion = Button(self, text="Pausar", command=self.pausar)
        self.btnIteracion.place(x=250, y=430)

        if (self.pausa == False):
            self.iteracion_Recursivo()




        # messagebox.showinfo(title="Entro la funcion", message=str(theta))
        # result = int(self.txtRes.get())


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

        # messagebox.showinfo(title="Entro la funcion", message=str(theta))
        # result = int(self.txtRes.get())

        self.txtRes.delete(0, 'end')
        self.txtRes.insert(0, self.theta)

        self.txtError.delete(0, 'end')
        self.txtError.insert(0, self.error)

    def recursivo_widgets(self):
        for child in self.winfo_children():
            child.destroy()

        self.create_widgets()

        self.btnArchivo = Button(self, text="Subir archivo", command = self.LeerArchivo)
        self.btnArchivo.place(x=150, y=180)
        Label(self, text="Cargar Archivo:").place(x=30, y=180)

        Label(self, text="Cantidad de datos inciales (N):").place(x=30, y=220)
        self.txtN = Entry(self, width=15)
        self.txtN.place(x=300, y=220)

        Label(self, text="Escribe las variables a's (empezando por a1)").place(x=30, y=250)
        self.txtn = Entry(self, width=15)
        self.txtn.place(x=300, y=250)

        Label(self, text="Escribe las variables b's (empezando por b0)").place(x=30, y=280)
        self.txtm = Entry(self, width=15)
        self.txtm.place(x=300, y=280)

        self.btnEmpezar = Button(self, text="Empezar", command=self.minimos_cuadrados_recursivo)
        self.btnEmpezar.place(x=150, y=310)

        Label(self, text="Iteración actual").place(x=30, y=340)
        self.txtIter = Entry(self, width=5)
        self.txtIter.place(x=130, y=340)

        Label(self, text="Resultado").place(x=30, y=370)
        self.txtRes = Entry(self, width=40)
        self.txtRes.place(x=100, y=370)

        Label(self, text="Error").place(x=30, y=400)
        self.txtError = Entry(self, width=20)
        self.txtError.place(x=100, y=400)

    def batch_widgets(self):
        for child in self.winfo_children():
            child.destroy()

        self.create_widgets()

        self.btnArchivo = Button(self, text="Subir archivo", command = self.LeerArchivo)
        self.btnArchivo.place(x=150, y=180)
        Label(self, text="Cargar Archivo:").place(x=30, y=180)

        Label(self, text="Cantidad de datos (N):").place(x=30, y=220)
        self.txtN = Entry(self, width=15)
        self.txtN.place(x=300, y=220)

        Label(self, text="Escribe las variables as (empezando por a1)").place(x=30, y=250)
        self.txtn = Entry(self, width=15)
        self.txtn.place(x=300, y=250)

        Label(self, text="Escribe las variables as (empezando por b0)").place(x=30, y=280)
        self.txtm = Entry(self, width=15)
        self.txtm.place(x=300, y=280)

        self.btnEmpezar = Button(self, text="Empezar", command=self.minimos_cuadrados_batch)
        self.btnEmpezar.place(x=150, y=310)

        Label(self, text="Resultado").place(x=30, y=340)
        self.txtRes = Entry(self, width=30)
        self.txtRes.place(x=100, y=340)

        Label(self, text="Error").place(x=30, y=370)
        self.txtError = Entry(self, width=20)
        self.txtError.place(x=100, y=370)


    def create_widgets(self):

        Label(self, text="Metodo:").place(x=30, y=90)
        Label(self, text="Forma:").place(x=30, y=120)


        self.btnCalcular = Button(self, text="Configuracion", command=self.ElegirOpciones)
        self.btnCalcular.place(x=100, y=150)


        self.cmbMetodos = Combobox(self, width="20", values=METODOS, state="readonly")
        self.cmbMetodos.place(x=100, y=90)
        self.cmbMetodos.current(0)

        self.cmbFormas = Combobox(self, width="20", values=FORMAS, state="readonly")
        self.cmbFormas.place(x=100, y=120)
        self.cmbFormas.current(0)




root = Tk()
root.title('Proyecto Procesos de Control')
root.geometry("800x700")
app = MainFrame(master=root)

app.mainloop()


