from tkinter import *
import matplotlib.pyplot as plt
import math
import numpy as np
import pandas as pd
import os


root = Tk()

# window size
root.geometry("1280x425+50+50") # window size and position
root.configure(bg="#5aa683")
root.title("Calculation of properties of Composite Laminates") # Title Name
root.resizable()

# Creating Variable value
zeroth_var = StringVar()
first_var = StringVar()
Second_var = StringVar()
Third_var = StringVar()
Foth_var = StringVar()
Fift_var = StringVar()
sixeth_var = StringVar()
seventh_var = StringVar()
eigth_var = StringVar()
nineth_var = StringVar()
tenth_var = StringVar()

# ----------------------For Dummy ----------------------#

# Creating Variable for Dummy value
value1 = StringVar()
value2 = StringVar()
value3 = StringVar()
value4 = StringVar()
value5 = StringVar()
value6 = StringVar()
value7 = StringVar()
value8 = StringVar()
value9 = StringVar()
value10 = StringVar()
value11 = StringVar()


# Funtion to get All Dummy value:
def getDummyvalue():
    try:
        print(float(value1.get()))
        print(float(value2.get()))
        print(float(value3.get()))
        print(float(value4.get()))
        print(float(value5.get()))
        print(float(value6.get()))
        print(float(value7.get()))
        print(float(value8.get()))
        print(float(value9.get()))
        print(float(value10.get()))
        print(float(value11.get()))
    except Exception as e:
        print(e)

#     Set All Dummy Value To Zero
    value1.set(0.0) # You Can Set any Value that you want, In my case i will set to 0.0, but you can leave it blank
    value2.set(0.0)
    value3.set(0.0)
    value4.set(0.0)
    value5.set(0.0)
    value6.set(0.0)
    value7.set(0.0)
    value8.set(0.0)
    value9.set(0.0)
    value10.set(0.0)
    value11.set(0.0)
# -------------------------END------------------------------#

result = IntVar()
result.set(2021.0)
# first_var.set(0.00)

# Creating Global Variable

Ef1 = 0.0
Ef2 = 0.0
nu_f = 0.0
Gf = 0.0
Em = 0.0
nu_m = 0.0
Vf = 0.0
# Gm = (0.5 * Em) / (1 + nu_m)


def out():
    print(Ef1, Ef2, nu_f, Gf, Em, nu_m, Vf, No_of_Layers, Thickness, Ply_Orientation,  Gm)

# -------------------------Your Code ----------------------------------------#
# ----------------------------Start------------------------------------------#
def run():
    class color:
       PURPLE = '\033[95m'
       CYAN = '\033[96m'
       DARKCYAN = '\033[36m'
       BLUE = '\033[94m'
       GREEN = '\033[92m'
       YELLOW = '\033[93m'
       RED = '\033[91m'
       BOLD = '\033[1m'
       UNDERLINE = '\033[4m'
       END = '\033[0m'
    
    def Effective_properties(Ef1, Ef2, Em, Vf):
        try:
            E1 = Vf * Ef1 + (1 - Vf) * Em
            E2 = (Ef2 * Em) / ((Ef2 * (1 - Vf)) + (Em * Vf))
            nu12 = Vf * nu_f + (1 - Vf) * nu_m
            G12 = (Gf * Gm) / ((Gm * Vf) + (Gf * (1 - Vf)))
        except Exception as error:
            print(error)
            pass
        print("Effective Properties of composite are:\nLongitudinal Stiffness of the composite,\n E1 =", E1)
        print("Transverse Stiffness of the composite, \n E2 =", E2)
        print("Poisson's ratio of the composite, \nnu12 =", nu12)
        print("Shear Modulus of the composite,\n G12 =", G12)
        return (E1, E2, nu12, G12)


    E1, E2, nu12, G12 = Effective_properties(Ef1, Ef2, Em, Vf)

    def compliance_matrix(E1, E2, nu12, G12):
        S11 = 1 / E1
        S12 = -nu12 / E1
        S16 = 0
        S21 = -nu12 / E1
        S22 = 1 / E2
        S26 = 0
        S61 = 0
        S62 = 0
        S66 = 1 / G12

        S = np.array([[S11, S12, S16], [S21, S22, S26], [S61, S62, S66]])

        return S

    def stiffness_matrix(E1, E2, nu12, G12):
        s_Matrix = compliance_matrix(E1, E2, nu12, G12)
        Q = np.linalg.inv(s_Matrix)
        print("The Compliance Matrix is=\n", s_Matrix)
        print("The Stiffness Matrix is=\n", Q)
        return Q

    S_Matrix = compliance_matrix(E1, E2, nu12, G12)
    Q_Matrix = stiffness_matrix(E1, E2, nu12, G12)

    def transformation_compliance(S, theta):
        c = math.cos(theta)
        s = math.sin(theta)

        S11 = S[0][0]
        S12 = S[0][1]
        S13 = S[0][2]
        S21 = S[1][0]
        S22 = S[1][1]
        S26 = S[0][2]
        S31 = S[2][0]
        S32 = S[2][1]
        S66 = S[2][2]

        S_xx = pow(c, 4) * S11 + pow(c, 2) * pow(s, 2) * (2 * S12 + S66) + pow(s, 4) * S22
        S_xy = (pow(c, 4) + pow(s, 4)) * S12 + pow(c, 2) * pow(s, 2) * (S11 + S22 - S66)
        S_yy = pow(s, 4) * S11 + pow(c, 2) * pow(s, 2) * (2 * S12 + S66) + pow(c, 4) * S22
        S_xs = 2 * pow(c, 3) * pow(s, 1) * (S11 - S12) + 2 * pow(c, 1) * pow(s, 3) * (S12 - S22) - pow(c, 1) * pow(s, 1) * (
                    pow(c, 2) - pow(s, 2)) * S66
        S_ys = 2 * pow(c, 1) * pow(s, 3) * (S11 - S12) + 2 * pow(c, 3) * pow(s, 1) * (S12 - S22) + pow(c, 1) * pow(s, 1) * (
                    pow(c, 2) - pow(s, 2)) * S66
        S_ss = 4 * pow(c, 2) * pow(s, 2) * (S11 - S12) - 4 * pow(c, 2) * pow(s, 2) * (S12 - S22) + pow(
            (pow(c, 2) - pow(s, 2)), 2) * S66

        S_yx = S_xy
        S_sx = S_xs
        S_sy = S_ys

        S_transformed = np.array([[S_xx, S_xy, S_xs], [S_yx, S_yy, S_ys], [S_sx, S_sy, S_ss]])

        print("Transformed Compliance Matrix is a general cordinate system is = \n", S_transformed)

        return S_transformed

    def transformation_stiffness(Q, theta):
        c = math.cos(theta)
        s = math.sin(theta)

        Q_11 = Q[0][0]
        Q_12 = Q[0][1]
        Q_16 = Q[0][2]
        Q_21 = Q[1][0]
        Q_22 = Q[1][1]
        Q_26 = Q[1][2]
        Q_61 = Q[2][0]
        Q_62 = Q[2][1]
        Q_66 = Q[2][2]

        Q_xx = pow(c, 4) * Q_11 + 2 * pow(c, 2) * pow(s, 2) * (Q_12 + 2 * Q_66) + pow(s, 4) * Q_22
        Q_xy = pow(c, 2) * pow(s, 2) * (Q_11 + Q_22 - 4 * Q_66) + Q_12 * (pow(c, 4) + pow(s, 4))
        Q_yy = pow(s, 4) * Q_11 + 2 * pow(c, 2) * pow(s, 2) * (Q_12 + 2 * Q_66) + pow(c, 4) * Q_22
        Q_xs = pow(c, 3) * s * (Q_11 - Q_12) + c * pow(s, 3) * (Q_12 - Q_22) - 2 * c * s * (pow(c, 2) - pow(s, 2)) * Q_66
        Q_ys = c * pow(s, 3) * (Q_11 - Q_12) + pow(c, 3) * s * (Q_12 - Q_22) + 2 * c * s * (pow(c, 2) - pow(s, 2)) * Q_66
        Q_ss = pow(c, 2) * pow(s, 2) * (Q_11 + Q_22 - 2 * Q_12) + pow((pow(c, 2) - pow(s, 2)), 2) * Q_66

        Q_yx = Q_xy
        Q_sx = Q_xs
        Q_sy = Q_ys

        Q_transformed = np.array([[Q_xx, Q_xy, Q_xs], [Q_yx, Q_yy, Q_ys], [Q_sx, Q_sy, Q_ss]])

        print("The rotated stiffness matrix in the X-Y coordinate system is = \n", Q_transformed)
        return Q_transformed

    S_xy = transformation_compliance(S_Matrix, 0)
    Q_xy = transformation_stiffness(Q_Matrix, 0)


    def Engineering_Constants_transformation(E1, E2, G12, nu12):
        theta = np.linspace(0, math.pi / 2, num=100)
        theta_degrees = np.linspace(0, 90, num=100)

        minus_eetaxy_x_matrix = []
        Gxy_by_G12_matrix = []
        ExbyE2_matrix = []
        nuxy_matrix = []
        E_x_matrix = []
        E_y_matrix = []
        G_xy_matrix = []

        for angle in theta:
            c = math.cos(angle)
            s = math.sin(angle)
            Ex = pow((c ** 4 / E1 + c ** 2 * s ** 2 * (1 / G12 - 2 * nu12 / E1) + s ** 4 / E2), -1)
            nuxy = Ex * ((nu12 * (s ** 4 + s ** 4)) / E1 - c ** 2 * s ** 2 * (1 / E1 + 1 / E2 - 1 / G12))
            Ey = pow((s ** 4 / E1 + c ** 2 * s ** 2 * (1 / G12 - 2 * nu12 / E1) + c ** 4 / E2), -1)
            eetaxy_x = Ex * (
                        c ** 3 * s * (2 / E1 + 2 * nu12 / E1 - 1 / G12) - c * s ** 3 * (2 / E2 + 2 * nu12 / E1 - 1 / G12))
            eetaxy_y = Ey * (
                        c * s ** 3 * (2 / E1 + 2 * nu12 / E1 - 1 / G12) - c ** 3 * s * (2 / E2 + 2 * nu12 / E1 - 1 / G12))
            Gxy = pow((4 * c ** 2 * s ** 2 * (1 / E1 + 1 / E2 + 2 * nu12 / E1) + (c ** 2 - s ** 2) ** 2 / G12), -1)

            E_x_matrix.append(Ex)
            E_y_matrix.append(Ey)
            G_xy_matrix.append(Gxy)

        plt.plot(theta_degrees, E_x_matrix)
        plt.plot(theta_degrees, E_y_matrix)
        plt.plot(theta_degrees, G_xy_matrix)
        plt.show()
        # ---------------NKS------------------#
        zeroth_var.set(0.0)
        first_var.set(0.0)
        Second_var.set(0.0)
        Third_var.set(0.0)
        Foth_var.set(0.0)
        Fift_var.set(0.0)
        sixeth_var.set(0.0)
        seventh_var.set(0.0)
        eigth_var.set(0.0)
        nineth_var.set(0.0)
        tenth_var.set(0)
        # ------------------------------------#


    Engineering_Constants_transformation(E1, E2, G12, nu12)


    def ABD_Matrix(No_of_Layers, totalThickness, Ply_Orientation):
        global z
        z = []

        N = No_of_Layers
        t = Thickness

        tk = t / N

        z0 = -N / 2 * tk
        z.append(z0)

        for i in range(0, N):
            presentCoordinate = z[-1] + tk
            z.append(presentCoordinate)

        Aij = np.array([[0] * 3] * 3)
        Bij = np.array([[0] * 3] * 3)
        Dij = np.array([[0] * 3] * 3)
        global Q_temp
        Q_temp = []

        for k in range(1, N + 1):
            Qij = transformation_stiffness(Q_Matrix, Ply_Orientation[k - 1])
            Aij = Aij + np.array(Qij) * (z[k] - z[k - 1])
            Bij = Bij + 0.5 * np.array(Qij) * (pow((z[k]), 2) - pow(z[k - 1], 2))
            Dij = Dij + (1 / 3) * np.array(Qij) * (pow((z[k]), 3) - pow(z[k - 1], 3))
            Q_temp.append(Qij)

        print("A matrix is:", Aij)
        print("B matrix is:", Bij)
        print("D matrix is:", Dij)

        return Aij, Bij, Dij

    No_of_Layers = 8
    Thickness = 1
    Ply_Orientation = [0, 0.767, -0.767, 1.57, 1.57, -0.767, -0.767, 0]

    A, B, D = ABD_Matrix(No_of_Layers, Thickness, Ply_Orientation)


# ----------------------------End-------------------------------------------#

def add():
    global Ef1, Ef2, nu_f, Gf,  Em, nu_m, Vf, No_of_Layers, Thickness, Ply_Orientation
    global Gm
    try:
        Ef1 = (float(zeroth_var.get()))
        Ef2 = (float(first_var.get()))
        nu_f = (float(Second_var.get()))
        Gf = (float(Third_var.get()))
        Em = (float(Foth_var.get()))
        nu_m = (float(Fift_var.get()))
        Vf = (float(sixeth_var.get()))
        # Gm = (float(seventh_var.get()))
        No_of_Layers = (float(eigth_var.get()))
        Thickness = (float(nineth_var.get()))
        temp = tenth_var.get()
        temp2 = (temp.split(","))
        Ply_Orientation =[]
        for nks in temp2:
            Ply_Orientation.append(float(nks))
        print(Ply_Orientation)
        Gm = (0.5 * Em) / (1 + nu_m)
        # result.set("")
    except:
        print("Error Found")
        pass
    # Rset to Zero all
    run()
    zeroth_var.set(0.0)
    first_var.set(0.0)
    Second_var.set(0.0)
    Third_var.set(0.0)
    Foth_var.set(0.0)
    Fift_var.set(0.0)
    sixeth_var.set(0.0)
    seventh_var.set(0.0)
    eigth_var.set(0.0)
    nineth_var.set(0.0)
    tenth_var.set(0)
    out()
    pass


def clear():
    zeroth_var.set(0.0)
    first_var.set(0.0)
    Second_var.set(0.0)
    Third_var.set(0.0)
    Foth_var.set(0.0)
    Fift_var.set(0.0)
    sixeth_var.set(0.0)
    seventh_var.set(0.0)
    eigth_var.set(0.0)
    nineth_var.set(0.0)
    tenth_var.set(0)
# Creating Entry Section

# Creating Text and Entry
font = "Roboto Mono"
label = Label(root,text="Longitudinal Modulus of Fiber, Ef1 (Gpa):" ,borderwidth=1,font=(font,15,'bold'),bg="#5a5663",fg="black", justify=LEFT, width=30)
label.grid(row=0,column=0,  ipadx=20, pady=2, padx=2)

label2 = Label(root,text="Transverse Modulus of Fiber, Ef2 (Gpa)" ,borderwidth=1,font=(font,15,'bold'),bg="#5a5663",fg="black", justify=LEFT, width=30)
label2.grid(row=1,column=0, ipadx=20, pady=2, padx=2)

label3 = Label(root,text="Poisson’s Ratio Fiber, nu_f",borderwidth=1, font=(font,15,'bold'),bg="#5a5663",fg="black", justify=LEFT, width=30)
label3.grid(row=2,column=0, sticky = W, ipadx=20, pady=2, padx=2)

label = Label(root,text="Shear Modulus of Fiber, Gf (Gpa)" ,borderwidth=1,font=(font,15,'bold'),bg="#5a5663",fg="black", justify=LEFT, width=30)
label.grid(row=3,column=0, sticky = W, ipadx=20, pady=2, padx=2)

label = Label(root,text="Modulus of Matrix, Em (Gpa)" ,borderwidth=1,font=(font,15,'bold'),bg="#5a5663",fg="black", justify=LEFT, width=30)
label.grid(row=4,column=0, sticky = W, ipadx=20, pady=2, padx=2)

label = Label(root,text="Poisson’s Ratio Matrix, nu_m" ,borderwidth=1,font=(font,15,'bold'),bg="#5a5663",fg="black", justify=LEFT, width=30)
label.grid(row=5,column=0, sticky = W, ipadx=20, pady=2, padx=2)

label = Label(root,text="Volume fraction of Fiber, Vf" ,borderwidth=1,font=(font,15,'bold'),bg="#5a5663",fg="black", justify=LEFT, width=30)
label.grid(row=6,column=0, sticky = W, ipadx=20, pady=2, padx=2)

label = Label(root,text="Number of Plies" ,borderwidth=1,font=(font,15,'bold'),bg="#5a5663",fg="black", justify=LEFT, width=30)
label.grid(row=7,column=0, sticky = W, ipadx=20, pady=2, padx=2)

label = Label(root,text="Thickness of Composite (mm)",borderwidth=1,font=(font,15,'bold'),bg="#5a5663",fg="black", justify=LEFT, width=30)
label.grid(row=8,column=0, sticky = W, ipadx=20, pady=2, padx=2)

label = Label(root,text="Orientation of Plies, put angles , seperated", borderwidth=1,font=(font,15,'bold'),bg="#5a5663",fg="black", justify=LEFT, width=30)
label.grid(row=9,column=0, sticky = W, ipadx=20, pady=2, padx=2)

button1= Button(root,borderwidth=1, text="Calculate", font=(font, 15, 'bold'), command = add,bg="#001254",fg="green", width=30)
button1.grid(row=10,column=0, pady=5, padx=2)

button1= Button(root,borderwidth=1, text="Clear Data", font=(font, 15, 'bold'), command = clear,bg="#001254",fg="green", width=30)
button1.grid(row=11,column=0, pady=5, padx=2)

# ------------------------------------------------ #


entry1=Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE,bg='white', textvariable=zeroth_var)
entry1.grid(row=0, column=1, padx=2)

entry1=Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE,bg='white', textvariable=first_var)
entry1.grid(row=1, column=1, padx=2)

entry1=Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE,bg='white', textvariable=Second_var)
entry1.grid(row=2, column=1, padx=2 )

entry1=Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE,bg='white', textvariable=Third_var)
entry1.grid(row=3, column=1, padx=2)

entry1=Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE,bg='white', textvariable=Foth_var)
entry1.grid(row=4, column=1, padx=2)

entry1=Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE,bg='white', textvariable=Fift_var)
entry1.grid(row=5, column=1, padx=2)

entry1=Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE,bg='white', textvariable=sixeth_var)
entry1.grid(row=6, column=1)

entry1=Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE,bg='white', textvariable=eigth_var)
entry1.grid(row=7, column=1, padx=2)

entry1=Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE,bg='white', textvariable=nineth_var)
entry1.grid(row=8, column=1, padx=2)

entry1=Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE,bg='white', textvariable=tenth_var)
entry1.grid(row=9, column=1, padx=2)

# entry1=Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE,bg='white', textvariable=result)
# entry1.grid(row=10, column=1, padx=2)

# For Another purpose

label = Label(root,text="Nx, (N-mm):" ,borderwidth=1,font=(font,15,'bold'),bg="#5a5663",fg="black", justify=LEFT, width=30)
label.grid(row=0,column=2, sticky = W, ipadx=20, pady=2, padx=2)

label = Label(root,text="Ny, (N-mm):" ,borderwidth=1,font=(font,15,'bold'),bg="#5a5663",fg="black", justify=LEFT, width=30)
label.grid(row=1,column=2, sticky = W, ipadx=20, pady=2, padx=2)

label = Label(root,text="Nxy, (N-mm):",borderwidth=1, font=(font,15,'bold'),bg="#5a5663",fg="black", justify=LEFT, width=30)
label.grid(row=2,column=2, sticky = W, ipadx=20, pady=2, padx=2)

label = Label(root,text="Mx, (N):" ,borderwidth=1,font=(font,15,'bold'),bg="#5a5663",fg="black", justify=LEFT, width=30)
label.grid(row=3,column=2, sticky = W, ipadx=20, pady=2, padx=2)

label = Label(root,text="My, (N):" ,borderwidth=1,font=(font,15,'bold'),bg="#5a5663",fg="black", justify=LEFT, width=30)
label.grid(row=4,column=2, sticky = W, ipadx=20, pady=2, padx=2)

label = Label(root,text="Mxy, (N):" ,borderwidth=1,font=(font,15,'bold'),bg="#5a5663",fg="black", justify=LEFT, width=30)
label.grid(row=5,column=2, sticky = W, ipadx=20, pady=2, padx=2)

label = Label(root,text="Longitudinal Tensile Strength, (MPa):" ,borderwidth=1,font=(font,15,'bold'),bg="#5a5663",fg="black", justify=LEFT, width=30)
label.grid(row=6,column=2, sticky = W, ipadx=20, pady=2, padx=2)

label = Label(root,text="Transverse Tensile Strength, (MPa):" ,borderwidth=1,font=(font,15,'bold'),bg="#5a5663",fg="black", justify=LEFT, width=30)
label.grid(row=7,column=2, sticky = W, ipadx=20, pady=2, padx=2)

label = Label(root,text="Longitudinal Compressive Strength, (MPa):",borderwidth=1,font=(font,15,'bold'),bg="#5a5663",fg="black", justify=LEFT, width=30)
label.grid(row=8,column=2, sticky = W, ipadx=20, pady=2, padx=2)

label = Label(root,text="Transverse Compressive Strength, (MPa):", borderwidth=1,font=(font,15,'bold'),bg="#5a5663",fg="black", justify=LEFT, width=30)
label.grid(row=9,column=2, sticky = W, ipadx=20, pady=2, padx=2)

label = Label(root,text="InterLaminar Shear Strength, (MPa)", borderwidth=1,font=(font,15,'bold'),bg="#5a5663",fg="black", justify=LEFT, width=30)
label.grid(row=9,column=2, sticky = W, ipadx=20, pady=2, padx=2)

# Note: put the Function name only, in command which you want to call. ie.command=function_name
#button1= Button(root,borderwidth=1, text="Submit", font=(font, 15, 'bold'), command =getDummyvalue,bg="#001254",fg="green", width=30)
#button1.grid(row=10,column=2, pady=5, padx=2)

#button1= Button(root,borderwidth=1, text="Clear Data", font=(font, 15, 'bold'), command =clear,bg="#001254",fg="green", width=30)
#button1.grid(row=11,column=2, pady=5, padx=2)

# ------------------------------------------------ #

entry1=Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE,bg='white', textvariable=value1)
entry1.grid(row=0, column=3, padx=2)

entry1=Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE,bg='white', textvariable=value2)
entry1.grid(row=1, column=3, padx=2)

entry1=Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE,bg='white', textvariable=value3)
entry1.grid(row=2, column=3, padx=2)

entry1=Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE,bg='white', textvariable=value4)
entry1.grid(row=3, column=3, padx=2)

entry1=Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE,bg='white', textvariable=value5)
entry1.grid(row=4, column=3, padx=2)

entry1=Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE,bg='white', textvariable=value6)
entry1.grid(row=5, column=3, padx=2)

entry1=Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE,bg='white', textvariable=value7)
entry1.grid(row=6, column=3, padx=2)

entry1=Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE,bg='white', textvariable=value8)
entry1.grid(row=7, column=3, padx=2)

entry1=Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE,bg='white', textvariable=value9)
entry1.grid(row=8, column=3, padx=2)

entry1=Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE,bg='white', textvariable=value10)
entry1.grid(row=9, column=3, padx=2)

entry1=Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE,bg='white', textvariable=value11)
entry1.grid(row=10, column=3, padx=2)

root.mainloop()


# -------------------------for save result in File-------------------------#
# G54 = 526
# file = open("result.txt", 'a')
# for i in range(10):
#     file.write("Shear Modulus of the composite, G12 =" + str(G54))
#     file.write("\n")
#
# file.write("________________________________________________________\n")
# file.close()

# --------------------------------------------------#
