from tkinter import *
import matplotlib.pyplot as plt
import math
import numpy as np
import pandas as pd
import os

root = Tk()

# window size
root.geometry("1280x425+50+50")  # window size and position
root.configure(bg="#5aa683")
root.title("Calculation of properties of Composite Laminates")  # Title Name
root.maxsize(1280, 425)
root.minsize(1280, 425)
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

Nx = 100.0
Ny = 0.0
Nxy = 0.0
Mx = 0.0
My = 0.0
Mxy = 0.0

Xt = 1500.0
Xc = 27.0
Yt = 570.0
Yc = 570.0
S12 = 34.0


# Gm = (0.5 * Em) / (1 + nu_m)


def out():
    print(Ef1, Ef2, nu_f, Gf, Em, nu_m, Vf, No_of_Layers, Thickness, Ply_Orientation, Gm)
    print(Nx, Ny, Nxy, Mx, My, Mxy, Xt, Xc, Yt, Yc, S12)


# -------------------------Your Code ----------------------------------------#
# ----------------------------Start------------------------------------------#
def run():
    f = open("result.txt", 'w')
    f.write("")
    f.close()
    file = open("result.txt", 'a')
    file.write("\n")

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

    # **********************************************Calculation of Effective Properties*****************************************************

    def Effective_properties(Ef1, Ef2, Em, Vf):  # Defined Function to calculate effective properties
        E1 = Vf * Ef1 + (1 - Vf) * Em  # Calculation of modulus in longitudinal direction
        E2 = (Ef2 * Em) / ((Ef2 * (1 - Vf)) + (Em * Vf))  # Calculation of modulus in Transverse direction
        nu12 = Vf * nu_f + (1 - Vf) * nu_m  # Calculation of Poisson's ratio of the laminate
        G12 = (Gf * Gm) / ((Gm * Vf) + (Gf * (1 - Vf)))  # Calculation of Shear Modulus of the laminate

        file.write("--------------------------------------------------------------------------------------------\n")
        file.write('The Effective Properties of the lamina are:\n')
        file.write("--------------------------------------------------------------------------------------------\n")
        file.write("Longitudinal Stiffness of the composite, E1 (Gpa) =" + str(E1))
        file.write("\n")
        file.write("Transverse Stiffness of the lamina, E2 (Gpa) =" + str(E2))
        file.write("\n")
        file.write("Poisson's ratio of the lamina, nu12 =" + str(nu12))
        file.write("\n")
        file.write("Shear Modulus of the lamina, G12 (Gpa) ="+ str(G12))
        file.write("\n")
        file.write("--------------------------------------------------------------------------------------------\n")
        return (E1, E2, nu12, G12)

    E1, E2, nu12, G12 = Effective_properties(Ef1, Ef2, Em, Vf)

    # **********************************************Calculation of Compliance Matrix*****************************************************

    def compliance_matrix(E1, E2, nu12, G12):  # Defined Function to calculate compliance matrix
        S11 = 1 / E1  # Defined formulation to calculate elements of compliance matrix
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

    S_Matrix = compliance_matrix(E1, E2, nu12, G12)  # Called function of Compliance matrix to define it golbally

    # **********************************************Calculation of Stiffness Matrix*****************************************************

    def stiffness_matrix(E1, E2, nu12, G12):  # Defined Function to calculate stiffness matrix
        s_Matrix = compliance_matrix(E1, E2, nu12, G12)  # Called function of compliance matrix
        Q = np.linalg.inv(s_Matrix)  # Defined formulation to calculate stiffness matrix
        # file.write("The Compliance Matrix is =\n\n",S_Matrix)
        # file.write("\n")
        # file.write("The Stiffness Matrix is =\n\n",Q)
        # file.write("\n")
        return Q

    Q_Matrix = stiffness_matrix(E1, E2, nu12, G12)  # Called function of Stiffness matrix to define it golbally

    # *************************************Calculation of Transformed Compliance Matrix****************************************************

    def transformation_compliance(S, theta):  # Defined Function to calculate transformation of compliance matrix
        theta_radian = theta * np.pi / 180
        c = np.cos(theta_radian)
        s = np.sin(theta_radian)

        S11 = S[0][0]  # Assigned values of elements of compliance matrix to variables for further calculations
        S12 = S[0][1]
        S13 = S[0][2]
        S21 = S[1][0]
        S22 = S[1][1]
        S26 = S[0][2]
        S31 = S[2][0]
        S32 = S[2][1]
        S66 = S[2][2]

        S_xx = c ** 4 * S11 + c ** 2 * s ** 2 * (
                2 * S12 + S66) + s ** 4 * S22  # Defined formulation to calculate elements of tranformed matrix
        S_xy = (c ** 4 + s ** 4) * S12 + c ** 2 * s ** 2 * (S11 + S22 - S66)
        S_yy = s ** 2 * S11 + c ** 2 * s ** 2 * (2 * S12 + S66) + c ** 4 * S22
        S_xs = 2 * c ** 3 * s * (S11 - S12) + 2 * c * s ** 3 * (S12 - S22) - c * s * (c ** 2 - s ** 2) * S66
        S_ys = 2 * c * s ** 3 * (S11 - S12) + 2 * c ** 3 * s * (S12 - S22) + c * s * (c ** 2 - s ** 2) * S66
        S_ss = 4 * c ** 2 * s ** 2 * (S11 - S12) - 4 * c ** 2 * s ** 2 * (S12 - S22) + ((c ** 2 - s ** 2) ** 2) * S66

        S_yx = S_xy
        S_sx = S_xs
        S_sy = S_ys

        S_transformed = np.array([[S_xx, S_xy, S_xs], [S_yx, S_yy, S_ys], [S_sx, S_sy, S_ss]]) / 10 ** 3

        # file.write("Transformed Compliance Matrix is a general cordinate system is = \n", S_transformed)
        # file.write("\n")
        return S_transformed

    # *************************************Calculation of Transformed Stiffness Matrix****************************************************

    def transformation_stiffness(Q, theta):  # Defined Function to calculate transformation of stiffness matrix
        theta_radian = theta * np.pi / 180
        c = np.cos(theta_radian)
        s = np.sin(theta_radian)

        Q_11 = Q[0][0]  # Assigned values of elements of stiffness matrix to variables for further calculations
        Q_12 = Q[0][1]
        Q_16 = Q[0][2]
        Q_21 = Q[1][0]
        Q_22 = Q[1][1]
        Q_26 = Q[1][2]
        Q_61 = Q[2][0]
        Q_62 = Q[2][1]
        Q_66 = Q[2][2]

        Q_xx = c ** 4 * Q_11 + 2 * c ** 2 * s ** 2 * (
                Q_12 + 2 * Q_66) + s ** 4 * Q_22  # Defined formulation to calculate elements of tranformed matrix
        Q_xy = c ** 2 * s ** 2 * (Q_11 + Q_22 - 4 * Q_66) + Q_12 * (c ** 4 + s ** 4)
        Q_yy = s ** 4 * Q_11 + 2 * c ** 2 * s ** 2 * (Q_12 + 2 * Q_66) + c ** 4 * Q_22
        Q_xs = c ** 3 * s * (Q_11 - Q_12) + c * s ** 3 * (Q_12 - Q_22) - 2 * c * s * (c ** 2 - s ** 2) * Q_66
        Q_ys = c * s ** 3 * (Q_11 - Q_12) + c ** 3 * s * (Q_12 - Q_22) + 2 * c * s * (c ** 2 - s ** 2) * Q_66
        Q_ss = c ** 2 * s ** 2 * (Q_11 + Q_22 - 2 * Q_12) + pow((c ** 2 - s ** 2), 2) * Q_66

        Q_yx = Q_xy
        Q_sx = Q_xs
        Q_sy = Q_ys

        Q_transformed = np.array([[Q_xx, Q_xy, Q_xs], [Q_yx, Q_yy, Q_ys], [Q_sx, Q_sy, Q_ss]]) * 10 ** 3
        # Unit of Q has been converted into MPa for convenience

        # file.write("The Transformed stiffness matrix in the X-Y coordinate system is = \n", Q_transformed)
        # file.write("\n")
        return Q_transformed

    S_xy = transformation_compliance(S_Matrix, 0)  # Called function of Compliance matrix to define it golbally
    Q_xy = transformation_stiffness(Q_Matrix, 0)  # Called function of Compliance matrix to define it golbally

    # *************************************Calculation of Transformed Engineering constants and plot graphs*******************************

    def Engineering_Constants_Transformation(E1, E2, G12,
                                             nu12):  # Defined function to calculate transformation of engineering constants

        theta = np.linspace(0, math.pi / 2, num=200)  # Defined list to set variation of theta
        theta_degrees = np.linspace(0, 90, num=200)  # Defined list to calculate variation of theta

        Q_xs_matrix = []  # Defined list to calculate variation of Q16
        E_x_matrix = []  # Defined list to calculate variation of transformation of E1
        E_y_matrix = []  # Defined list to calculate variation of transformation of E2
        G_xy_matrix = []  # Defined list to calculate variation of transformation of G12

        for angle in theta:
            c = math.cos(angle)
            s = math.sin(angle)
            Q11 = E1 / (1 - nu12 ** 2 * (E2 / E1))  # Defined formulation to calculate elements of stiffness matrix
            Q22 = E2 / (1 - nu12 ** 2 * (E2 / E1))
            Q12 = (nu12 * E2) / (1 - nu12 ** 2 * (E2 / E1))
            Q66 = G12
            Ex = (c ** 4 / E1 + c ** 2 * s ** 2 * (1 / G12 - 2 * nu12 / E1) + s ** 4 / E2) ** (
                -1)  # Defined formulation to calculate engineering contants
            nuxy = Ex * ((nu12 * (s ** 4 + s ** 4)) / E1 - c ** 2 * s ** 2 * (1 / E1 + 1 / E2 - 1 / G12))
            Ey = (s ** 4 / E1 + c ** 2 * s ** 2 * (1 / G12 - 2 * nu12 / E1) + c ** 4 / E2) ** (-1)
            eetaxy_x = Ex * (
                    c ** 3 * s * (2 / E1 + 2 * nu12 / E1 - 1 / G12) - c * s ** 3 * (2 / E2 + 2 * nu12 / E1 - 1 / G12))
            eetaxy_y = Ey * (
                    c * s ** 3 * (2 / E1 + 2 * nu12 / E1 - 1 / G12) - c ** 3 * s * (2 / E2 + 2 * nu12 / E1 - 1 / G12))
            Gxy = pow((4 * c ** 2 * s ** 2 * (1 / E1 + 1 / E2 + 2 * nu12 / E1) + (c ** 2 - s ** 2) ** 2 / G12), -1)
            Q_xs = c ** 3 * s ** 1 * (Q11 - Q12) + c ** 1 * s ** 3 * (Q12 - Q22) - 2 * c ** 1 * s ** 1 * (
                    c ** 2 - s ** 2) * Q66

            E_x_matrix.append(Ex)
            E_y_matrix.append(Ey)
            G_xy_matrix.append(Gxy)
            Q_xs_matrix.append(Q_xs)

        plot1 = plt.figure(1)  # Creation of plot of Ex, Ey, Gxy in variation with theta
        plt.grid()
        plt.plot(theta_degrees, E_x_matrix)
        plt.plot(theta_degrees, E_y_matrix)
        plt.plot(theta_degrees, G_xy_matrix)
        plt.title("Variation of different Modulus with \u03F4\n", fontsize=16)
        plt.xlabel("\u03F4 ", fontsize=12)
        plt.ylabel("Ex,Ey,Gxy ", fontsize=12)
        plt.legend(["Ex", "Ey", "Gxy"], prop={"size": 12}, loc=1, ncol=3)

        plot2 = plt.figure(2)  # Creation of plot of Q16 in variation with theta
        plt.grid()
        plt.plot(theta_degrees, Q_xs_matrix)
        plt.title("Variation of Q16 with \u03F4\n", fontsize=16)
        plt.xlabel("\u03F4 ", fontsize=12)
        plt.ylabel("Q16 ", fontsize=12)

        plt.show()

    Engineering_Constants_Transformation(E1, E2, G12, nu12)  # Called function to plot the Graph

    # ******************************************Calculation of ABD matrix***************************************************************

    def ABD_Matrix(No_of_Layers, Thickness, Ply_Orientation):  # Defined function to calculate ABD Matrix
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

        A = np.zeros((3, 3))  # Unit of A GN-mm
        B = np.zeros((3, 3))  # Unit of B GN
        D = np.zeros((3, 3))  # Unit of C GN/mm

        global Q_temp
        Q_temp = []

        for j in range(1, N + 1):
            # file.write(f"For Ply number {k}:\n")
            Q = transformation_stiffness(Q_Matrix,
                                         Ply_Orientation[j - 1])  # Defined formulation to calculate A, B, D Matrix
            A = A + np.array(Q) * (z[j] - z[j - 1])
            B = B + (1 / 2) * np.array(Q) * (((z[j]) ** 2) - (z[j - 1]) ** 2)
            D = D + (1 / 3) * np.array(Q) * (((z[j]) ** 3) - (z[j - 1]) ** 3)

            Q_temp.append(Q)

        # file.write("\n A Matrix is = \n", A)
        # file.write("\n")
        # file.write("\n B Matrix is = \n", B)
        # file.write("\n")
        # file.write("\n D Matrix is = \n", D)
        # file.write("\n")

        return A, B, D

    # **************************************Values for of ABD matrix Calculations****************************************************

    No_of_Layers = 8
    Thickness = 1
    Ply_Orientation = [0, 22, -22, 90, 90, -22, 22, 0]

    A, B, D = ABD_Matrix(No_of_Layers, Thickness, Ply_Orientation)

    #  ___________________New Program_____________________

    Row_1 = np.concatenate([A, B])
    Row_2 = np.concatenate([B, D])

    ABD_Matrix = np.concatenate((Row_1, Row_2), axis=1)  # Code for the formation of ABD complete matrix
    ABD_Matrix_1 = np.around(ABD_Matrix, decimals=5)

    # file.write("\nABD Matrix is: = \n",ABD_Matrix_1)

    # **************************************External Load definitions***************************************************

    # Nx = 100
    # Ny = 0
    # Nxy = 0
    # Mx = 0
    # My = 0
    # Mxy = 0

    # **************************************Calculations of NM matrix and Epsilon_0 Matrix***************************************************

    N_bar = np.array([[Nx], [Ny], [Nxy]])  # Formation of N matrix
    M_bar = np.array([[Mx], [My], [Mxy]])  # Formation of M matrix

    Load_Matrix = np.concatenate((N_bar, M_bar))
    # file.write("Load matrix = \n",Load_Matrix)
    Epsilon_K = np.matmul(np.linalg.inv(ABD_Matrix), Load_Matrix)
    abd_matrix = np.linalg.inv(ABD_Matrix)
    abd_matrix_1 = np.around(abd_matrix, decimals=3)
    Epsilon_0 = np.vsplit(Epsilon_K, 2)[0]  # Calculations of Epsilon0 matrix
    Kappa = np.split(Epsilon_K, 2)[1]  # Calculations of Kappa matrix

    file.write("\n--------------------------------------------------------------------------------------------\n")
    file.write('______________________ABD Matrix is:________________________\n')
    file.write("\n" + str(ABD_Matrix_1))
    file.write("\n For further calculations convenience the unit of A here is GPa-mm, B is GPa and D is GPa/mm. \n")
    file.write("\n--------------------------------------------------------------------------------------------\n")
    file.write(f"\n \n \n")
    # file.write("\n Epsilon K = \n",Epsilon_K)
    # file.write("\n abd_matrix is = \n",abd_matrix_1)
    # file.write("\n Epsilon_0 matrix is = \n",Epsilon_0)
    # file.write("\n Kappa matrix is = \n",Kappa)
    # file.write("\n No of layers = \n",No_of_Layers)

    # *********************************Calculations of Epsilon, Sigma k and Sigma k 1-2 direction**********************************************

    Epsilon_store = []
    Sigma_K_store = []
    T_sigma_store = []
    for i in range(0, No_of_Layers):
        z_average = (z[i] + z[i + 1]) / 2
        Epsilon = np.array(Epsilon_0) + z_average * np.array(Kappa)
        Epsilon_store.append(Epsilon)
        # file.write(f"\n Epsilon Matrix {i+1}\n", Epsilon)

    Sigma_Fiber_CS_store = []
    theta_radian_store = []
    for i in range(0, No_of_Layers):
        theta = Ply_Orientation[i]
        theta_radian = theta * np.pi / 180
        theta_radian_store.append(theta_radian)
        c = np.cos(theta_radian)
        s = np.sin(theta_radian)
        T_sigma = np.array(
            [[c ** 2, s ** 2, 2 * c * s], [s ** 2, c ** 2, -2 * c * s], [-c * s, c * s, c ** 2 - s ** 2]])
        T_sigma_store.append(T_sigma)
        Sigma_K = np.dot(Q_temp[i], Epsilon)
        Sigma_K_store.append(Sigma_K)
        Sigma_fiber_direction = np.dot(T_sigma, Sigma_K)
        Sigma_Fiber_CS_store.append(Sigma_fiber_direction)
        # file.write(f"\n T Sigma in fiber direction co-ordinate system for Layer {i+1}\n", T_sigma)
        # file.write(f"\n Sigma in fiber direction co-ordinate system for Layer {i+1}\n", Sigma_fiber_direction)
    # file.write(Sigma_K_store)
    # file.write(theta_radian_store)

    # **********************************************Definition of Strength Parameter**********************************************

    # Xt = 1500
    # Xc = 27
    # Yt = 570
    # Yc = 570
    # S12 = 34

    # ********************************Definition of Global list for Hashin criteria implementation********************************

    R_Tensile_M = []
    R_Compressive_M = []
    R_Tensile_F = []
    R_Compressive_F = []
    Temp_List = []
    R_Final = []

    # ***********************************Hashin criteria implementation & Margin of safety calculations************************************

    def Hashin_Criteria(Xt, Xc, Yt, Yc, S12):
        for i in range(0, No_of_Layers):
            S23 = S12
            Sigma_33 = 0
            Sigma_Fiber_CS = Sigma_Fiber_CS_store[i]
            Sigma_11 = Sigma_Fiber_CS[0]
            Sigma_22 = Sigma_Fiber_CS[1]
            Sigma_66 = Sigma_Fiber_CS[2]
            Condition_Matrix = Sigma_22 + Sigma_33
            Condition_Fiber = Sigma_11

            # file.write(Condition_Matrix)
            # file.write(Condition_Fiber)

            if (Condition_Matrix >= 0):
                # file.write("--------------------------------------------------------------------------------------------\n")
                # file.write("Checking Tensile failure mode in Matrix : \n")

                Failure_Value = (Sigma_22 / Yt ** 2) ** 2 + (Sigma_66 / S12) ** 2
                # file.write("Failure value is : \n", Failure_Value)

                if (Failure_Value >= 1):

                    # file.write(f"Matrix of Ply {i+1} failed in Tensile Mode\n")
                    # file.write(color.RED + color.BOLD + 'PLEASE CHECK' + color.END)

                    c1 = (Sigma_22 / Yt ** 2) ** 2 + (Sigma_66 / S12) ** 2
                    c2 = 0
                    c3 = -1
                    R = np.roots([c1, c2, c3])
                    if (R[0] > R[1]):
                        R_Tensile_M.append(R[0])
                    else:
                        R_Tensile_M.append(R[1])

                    # R_Tensile_M.append(R)
                    R_Compressive_M.append(100000000000)

                else:
                    c1 = (Sigma_22 / Yt ** 2) ** 2 + (Sigma_66 / S12) ** 2
                    c2 = 0
                    c3 = -1
                    R = np.roots([c1, c2, c3])

                    if (R[0] > R[1]):
                        R_Tensile_M.append(R[0])
                    else:
                        R_Tensile_M.append(R[1])

                    # R_Tensile_M.append(R)
                    R_Compressive_M.append(100000000000)

                    # file.write(f"Margin of safety for ply {i+1} is \n", R)

                    # file.write(f"Matrix of Ply {i+1} is safe from Tensile Mode Failure\n")


            elif (Condition_Matrix < 0):
                # file.write("Checking Compressive failure mode in Matrix : \n")

                Failure_Value = ((Yc / (2 * S23)) ** 2 - 1) * (Sigma_22 / Yt) + (Sigma_22 / (4 * S23)) ** 2 + (
                        Sigma_66 / S12) ** 2
                # file.write("Failure value is : \n", Failure_Value)

                if (Failure_Value >= 1):

                    # file.write(f"Matrix of Ply {i+1} + \033[1m + failed in Compressive Mode\n")
                    # file.write(color.RED + color.BOLD + 'PLEASE CHECK' + color.END)

                    c1 = (Sigma_22 / (4 * S23)) ** 2 + (Sigma_66 / S12) ** 2
                    c2 = ((Yc / (2 * S23)) ** 2 - 1) * (Sigma_22 / Yt)
                    c3 = -1
                    R = np.roots([c1, c2, c3])

                    if (R[0] > R[1]):
                        R_Compressive_M.append(R[0])
                    else:
                        R_Compressive_M.append(R[1])

                    # R_Compressive_M.append(R)
                    R_Tensile_M.append(100000000000)

                else:
                    c1 = (Sigma_22 / (4 * S23)) ** 2 + (Sigma_66 / S12) ** 2
                    c2 = ((Yc / (2 * S23)) ** 2 - 1) * (Sigma_22 / Yt)
                    c3 = -1
                    R = np.roots([c1, c2, c3])

                    if (R[0] > R[1]):
                        R_Compressive_M.append(R[0])
                    else:
                        R_Compressive_M.append(R[1])

                    # R_Compressive_M.append(R)
                    R_Tensile_M.append(100000000000)

                    # file.write(f"Margin of safety for ply {i+1} is \n", R)
                    # file.write(f"Matrix of Ply {i+1} is safe from Compressive Mode Failure\n")
                    # file.write("--------------------------------------------------------------------------------------------\n")

            if (Condition_Fiber >= 0):

                # file.write("Checking Tensile failure mode in Fiber : \n")

                Failure_Value = (Sigma_11 / Xt) ** 2 + (Sigma_66 / S12) ** 2

                # file.write("Failure value is : \n", Failure_Value)

                if (Failure_Value >= 1):

                    # file.write(f"Fibers of Ply {i+1} failed in Tensile Mode\n")
                    # file.write(color.RED + color.BOLD + 'PLEASE CHECK' + color.END)

                    c1 = (Sigma_11 / Xt) ** 2 + (Sigma_66 / S12) ** 2
                    c2 = 0
                    c3 = -1
                    R = np.roots([c1, c2, c3])

                    if (R[0] > R[1]):
                        R_Tensile_F.append(R[0])
                    else:
                        R_Tensile_F.append(R[1])

                    # R_Tensile_F.append(R)
                    R_Compressive_F.append(100000000000)

                else:
                    c1 = (Sigma_11 / Xt) ** 2 + (Sigma_66 / S12) ** 2
                    c2 = 0
                    c3 = -1
                    R = np.roots([c1, c2, c3])

                    if (R[0] > R[1]):
                        R_Tensile_F.append(R[0])
                    else:
                        R_Tensile_F.append(R[1])

                    # R_Tensile_F.append(R)
                    R_Compressive_F.append(100000000000)

                    # file.write(f"Margin of safety for ply {i+1} is \n", R)
                    # file.write(f"Fibers of Ply {i+1} is safe from Tensile Mode Failure\n")

            elif (Condition_Fiber < 0):
                # file.write("Checking Compressive failure mode in Fiber : \n")

                Failure_Value = (Sigma_11 / Xc) ** 2

                # file.write("Failure value is : \n", Failure_Value)

                if (Failure_Value >= 1):

                    # file.write(f"Fibers of Ply {i+1} failed in Compressive Mode\n")
                    # file.write(color.RED + color.BOLD + 'PLEASE CHECK' + color.END)

                    c1 = (Sigma_11 / Xc) ** 2
                    c2 = 0
                    c3 = -1
                    R = np.roots([c1, c2, c3])

                    if (R[0] > R[1]):
                        R_Compressive_F.append(R[0])
                    else:
                        R_Compressive_F.append(R[1])

                    # R_Compressive_F.append(R)
                    R_Tensile_F.append(100000000000)

                else:
                    c1 = (Sigma_11 / Xc) ** 2
                    c2 = 0
                    c3 = -1
                    Rf = np.roots([c1, c2, c3])

                    if (R[0] > R[1]):
                        R_Compressive_F.append(R[0])
                    else:
                        R_Compressive_F.append(R[1])

                    # R_Compressive_F.append(R)
                    R_Tensile_F.append(100000000000)

                    # file.write(f"Margin of safety for ply {i+1} is \n", R)
                    # file.write(f"Fibers of Ply {i+1} is safe from Compressive Mode Failure\n")

        for i in range(0, No_of_Layers):
            Temp_List = [R_Compressive_F[i], R_Tensile_F[i], R_Compressive_M[i], R_Tensile_M[i]]
            Temp_List.sort()
            # file.write(Temp_List)
            R_Final.append(Temp_List[0])

        # file.write(R_Final)
        file.write('Calculation of Margin of safety using Hashins Criteria\n')

        # for i in range(0, No_of_Layers):
        # file.write(f"Margin of Safety for Ply {i+1} is: ", R_Final[i])
        # file.write("-------------------------------------------------------------------\n")

        file.write("---------------------------------------------------------------------\n")
        file.write("Laminate no.        ||  Ply Orientation    ||  R_intact \n")
        file.write("+------------++------------------++-----------------------------------++\n")
        for i in range(0, No_of_Layers):
            file.write("||    "+str(i + 1).ljust(5)+"         ||        " + str(Ply_Orientation[i]).ljust(5)+ "        || " + str(R_Final[i]).ljust(5))
            file.write("\n")
        file.write("+---------------------------------++-------------------------------------++\n")
        file.write(f"Hence, for 0 degree ply the margin of safety is : " + str(R_Final[0])+"\n")
        file.write(f"Hence, for 90 degree ply the margin of safety is :" + str(R_Final[3])+"\n")
        file.write("----------------------------------------------------------------------\n")
        file.write('Since the Margin of Safety for 90 degree laminate is less than 1, hence it has failed')
        file.write("\n---------------------------------------------------------------------\n")
        file.close()

    Hashin_Criteria(Xt, Xc, Yt, Yc, S12)
    # file = open("result.txt", 'a')
    # file.write("\n")
    file.close()


# ----------------------------End-------------------------------------------#

def add():
    global Ef1, Ef2, nu_f, Gf, Em, nu_m, Vf, No_of_Layers, Thickness, Ply_Orientation
    global Gm, Nx, Ny, Nxy, Mx, My, Mxy, Xt, Xc, Yt, Yc, S12
    try:
        Ef1 = (float(zeroth_var.get()))
        Ef2 = (float(first_var.get()))
        nu_f = (float(Second_var.get()))
        Gf = (float(Third_var.get()))
        Em = (float(Foth_var.get()))
        nu_m = (float(Fift_var.get()))
        Vf = (float(sixeth_var.get()))
        No_of_Layers = (float(eigth_var.get()))
        Thickness = (float(nineth_var.get()))
        temp = tenth_var.get()
        temp2 = (temp.split(","))
        Ply_Orientation = []
        for nks in temp2:
            Ply_Orientation.append(float(nks))
        # file.write(Ply_Orientation)
        Nx = float(value1.get())
        Ny = float(value2.get())
        Nxy = float(value3.get())
        Mx = float(value4.get())
        My = float(value5.get())
        Mxy = float(value6.get())
        Xt = float(value7.get())
        Xc = float(value8.get())
        Yt = float(value9.get())
        Yc = float(value10.get())
        S12 = float(value11.get())

        Gm = (0.5 * Em) / (1 + nu_m)
        # result.set("")
    except:
        print("Error Found To get Value Line: 727")
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

    value1.set(0.0)
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

    value1.set(0.0)
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


# Creating Entry Section

# Creating Text and Entry
font = "Roboto Mono"
label = Label(root, text="Longitudinal Modulus of Fiber, Ef1 (Gpa):", borderwidth=1, font=(font, 15, 'bold'),
              bg="#5a5663", fg="black", justify=LEFT, width=30)
label.grid(row=0, column=0, ipadx=20, pady=2, padx=2)

label2 = Label(root, text="Transverse Modulus of Fiber, Ef2 (Gpa)", borderwidth=1, font=(font, 15, 'bold'),
               bg="#5a5663", fg="black", justify=LEFT, width=30)
label2.grid(row=1, column=0, ipadx=20, pady=2, padx=2)

label3 = Label(root, text="Poisson’s Ratio Fiber, nu_f", borderwidth=1, font=(font, 15, 'bold'), bg="#5a5663",
               fg="black", justify=LEFT, width=30)
label3.grid(row=2, column=0, sticky=W, ipadx=20, pady=2, padx=2)

label = Label(root, text="Shear Modulus of Fiber, Gf (Gpa)", borderwidth=1, font=(font, 15, 'bold'), bg="#5a5663",
              fg="black", justify=LEFT, width=30)
label.grid(row=3, column=0, sticky=W, ipadx=20, pady=2, padx=2)

label = Label(root, text="Modulus of Matrix, Em (Gpa)", borderwidth=1, font=(font, 15, 'bold'), bg="#5a5663",
              fg="black", justify=LEFT, width=30)
label.grid(row=4, column=0, sticky=W, ipadx=20, pady=2, padx=2)

label = Label(root, text="Poisson’s Ratio Matrix, nu_m", borderwidth=1, font=(font, 15, 'bold'), bg="#5a5663",
              fg="black", justify=LEFT, width=30)
label.grid(row=5, column=0, sticky=W, ipadx=20, pady=2, padx=2)

label = Label(root, text="Volume fraction of Fiber, Vf", borderwidth=1, font=(font, 15, 'bold'), bg="#5a5663",
              fg="black", justify=LEFT, width=30)
label.grid(row=6, column=0, sticky=W, ipadx=20, pady=2, padx=2)

label = Label(root, text="Number of Plies", borderwidth=1, font=(font, 15, 'bold'), bg="#5a5663", fg="black",
              justify=LEFT, width=30)
label.grid(row=7, column=0, sticky=W, ipadx=20, pady=2, padx=2)

label = Label(root, text="Thickness of Composite (mm)", borderwidth=1, font=(font, 15, 'bold'), bg="#5a5663",
              fg="black", justify=LEFT, width=30)
label.grid(row=8, column=0, sticky=W, ipadx=20, pady=2, padx=2)

label = Label(root, text="Orientation of Plies, put angles , seperated", borderwidth=1, font=(font, 15, 'bold'),
              bg="#5a5663", fg="black", justify=LEFT, width=30)
label.grid(row=9, column=0, sticky=W, ipadx=20, pady=2, padx=2)

button1 = Button(root, borderwidth=1, text="Calculate", font=(font, 15, 'bold'), command=add, bg="#001254", fg="green",
                 width=30)
button1.grid(row=10, column=0, pady=5, padx=2)

button1 = Button(root, borderwidth=1, text="Clear Data", font=(font, 15, 'bold'), command=clear, bg="#001254",
                 fg="green", width=30)
button1.grid(row=11, column=0, pady=5, padx=2)

# ------------------------------------------------ #


entry1 = Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE, bg='white', textvariable=zeroth_var)
entry1.grid(row=0, column=1, padx=2)

entry1 = Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE, bg='white', textvariable=first_var)
entry1.grid(row=1, column=1, padx=2)

entry1 = Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE, bg='white', textvariable=Second_var)
entry1.grid(row=2, column=1, padx=2)

entry1 = Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE, bg='white', textvariable=Third_var)
entry1.grid(row=3, column=1, padx=2)

entry1 = Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE, bg='white', textvariable=Foth_var)
entry1.grid(row=4, column=1, padx=2)

entry1 = Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE, bg='white', textvariable=Fift_var)
entry1.grid(row=5, column=1, padx=2)

entry1 = Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE, bg='white', textvariable=sixeth_var)
entry1.grid(row=6, column=1)

entry1 = Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE, bg='white', textvariable=eigth_var)
entry1.grid(row=7, column=1, padx=2)

entry1 = Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE, bg='white', textvariable=nineth_var)
entry1.grid(row=8, column=1, padx=2)

entry1 = Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE, bg='white', textvariable=tenth_var)
entry1.grid(row=9, column=1, padx=2)

# entry1=Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE,bg='white', textvariable=result)
# entry1.grid(row=10, column=1, padx=2)

# For Another purpose

label = Label(root, text="Nx, (N-mm):", borderwidth=1, font=(font, 15, 'bold'), bg="#5a5663", fg="black", justify=LEFT,
              width=30)
label.grid(row=0, column=2, sticky=W, ipadx=20, pady=2, padx=2)

label = Label(root, text="Ny, (N-mm):", borderwidth=1, font=(font, 15, 'bold'), bg="#5a5663", fg="black", justify=LEFT,
              width=30)
label.grid(row=1, column=2, sticky=W, ipadx=20, pady=2, padx=2)

label = Label(root, text="Nxy, (N-mm):", borderwidth=1, font=(font, 15, 'bold'), bg="#5a5663", fg="black", justify=LEFT,
              width=30)
label.grid(row=2, column=2, sticky=W, ipadx=20, pady=2, padx=2)

label = Label(root, text="Mx, (N):", borderwidth=1, font=(font, 15, 'bold'), bg="#5a5663", fg="black", justify=LEFT,
              width=30)
label.grid(row=3, column=2, sticky=W, ipadx=20, pady=2, padx=2)

label = Label(root, text="My, (N):", borderwidth=1, font=(font, 15, 'bold'), bg="#5a5663", fg="black", justify=LEFT,
              width=30)
label.grid(row=4, column=2, sticky=W, ipadx=20, pady=2, padx=2)

label = Label(root, text="Mxy, (N):", borderwidth=1, font=(font, 15, 'bold'), bg="#5a5663", fg="black", justify=LEFT,
              width=30)
label.grid(row=5, column=2, sticky=W, ipadx=20, pady=2, padx=2)

label = Label(root, text="Longitudinal Tensile Strength, (MPa):", borderwidth=1, font=(font, 15, 'bold'), bg="#5a5663",
              fg="black", justify=LEFT, width=30)
label.grid(row=6, column=2, sticky=W, ipadx=20, pady=2, padx=2)

label = Label(root, text="Transverse Tensile Strength, (MPa):", borderwidth=1, font=(font, 15, 'bold'), bg="#5a5663",
              fg="black", justify=LEFT, width=30)
label.grid(row=7, column=2, sticky=W, ipadx=20, pady=2, padx=2)

label = Label(root, text="Longitudinal Compressive Strength, (MPa):", borderwidth=1, font=(font, 15, 'bold'),
              bg="#5a5663", fg="black", justify=LEFT, width=30)
label.grid(row=8, column=2, sticky=W, ipadx=20, pady=2, padx=2)

label = Label(root, text="Transverse Compressive Strength, (MPa):", borderwidth=1, font=(font, 15, 'bold'),
              bg="#5a5663", fg="black", justify=LEFT, width=30)
label.grid(row=9, column=2, sticky=W, ipadx=20, pady=2, padx=2)

label = Label(root, text="InterLaminar Shear Strength, (MPa)", borderwidth=1, font=(font, 15, 'bold'), bg="#5a5663",
              fg="black", justify=LEFT, width=30)
label.grid(row=10, column=2, sticky=W, ipadx=20, pady=2, padx=2)

# Note: put the Function name only, in command which you want to call. ie.command=function_name
# button1= Button(root,borderwidth=1, text="Submit", font=(font, 15, 'bold'), command =getDummyvalue,bg="#001254",fg="green", width=30)
# button1.grid(row=10,column=2, pady=5, padx=2)

# button1= Button(root,borderwidth=1, text="Clear Data", font=(font, 15, 'bold'), command =clear,bg="#001254",fg="green", width=30)
# button1.grid(row=11,column=2, pady=5, padx=2)

# ------------------------------------------------ #

entry1 = Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE, bg='white', textvariable=value1)
entry1.grid(row=0, column=3, padx=2)

entry1 = Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE, bg='white', textvariable=value2)
entry1.grid(row=1, column=3, padx=2)

entry1 = Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE, bg='white', textvariable=value3)
entry1.grid(row=2, column=3, padx=2)

entry1 = Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE, bg='white', textvariable=value4)
entry1.grid(row=3, column=3, padx=2)

entry1 = Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE, bg='white', textvariable=value5)
entry1.grid(row=4, column=3, padx=2)

entry1 = Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE, bg='white', textvariable=value6)
entry1.grid(row=5, column=3, padx=2)

entry1 = Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE, bg='white', textvariable=value7)
entry1.grid(row=6, column=3, padx=2)

entry1 = Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE, bg='white', textvariable=value8)
entry1.grid(row=7, column=3, padx=2)

entry1 = Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE, bg='white', textvariable=value9)
entry1.grid(row=8, column=3, padx=2)

entry1 = Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE, bg='white', textvariable=value10)
entry1.grid(row=9, column=3, padx=2)

entry1 = Entry(root, borderwidth=1, font=(font, 15, 'bold'), relief=GROOVE, bg='white', textvariable=value11)
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
