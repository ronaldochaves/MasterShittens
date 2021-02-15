import os
import numpy as np
import pandas as pd
import tkinter as tk
import matplotlib.pyplot as plt
import matplotlib.patches as ptc
import matplotlib.lines as lines
import matplotlib.image as img
from scipy.interpolate import interp1d
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from tkinter import ttk


class point:
    T = 0.0
    p = 0.0
    T0 = 0.0
    p0 = 0.0
    Ts = 0.0
    rho = 0.0
    A = 0.0
    U = 0.0
    V = 0.0
    W = 0.0
    Va = 0.0
    Vu = 0.0
    Wa = 0.0
    Wu = 0.0
    alpha = 0.0
    beta = 0.0
    Mach = 0.0
    Machrel = 0.0
    mu = 0.0


class geometry:
    A = 0.0
    h = 0.0
    rm = 0.0
    rt = 0.0
    rh = 0.0
    rr = 0.0


class blades:
    h = 0.0
    c = 0.0
    s = 0.0
    w = 0.0
    rm = 0.0
    rt = 0.0
    rh = 0.0
    rr = 0.0
    n = int(0)
    hc = 0.0
    sc = 0.0
    hw = 0.0
    tc = 0.0
    kh = 0.0
    tes = 0.0
    teo = 0.0
    in_angle = 0.0
    out_angle = 0.0


class losses:
    Ypn = 0.0
    Ysn = 0.0
    Ykn = 0.0
    Yten = 0.0
    Ytn = 0.0
    lbdn = 0.0
    effn = 0.0
    Ypr = 0.0
    Ysr = 0.0
    Ykr = 0.0
    Yter = 0.0
    Ytr = 0.0
    lbdr = 0.0
    effr = 0.0


def graph_smith_chart():
    nrows = data_smith[0].size
    index_inflection = []
    for k in range(0, 7):
        inflection = np.max(data_smith[2*k+1])
        for i in range(0, nrows-1):
            if (data_smith[2*k+1][i] == inflection):
                index_inflection.append(i)
                break
        for j in range(nrows-1, i, -1):
            for jj in range(nrows-1, i+nrows-j, -1):
                if (data_smith[2*k+1][jj] > data_smith[2*k+1][jj-1]):
                    data_smith[2*k][jj], data_smith[2*k][jj -
                                                         1] = data_smith[2*k][jj-1], data_smith[2*k][jj]
                    data_smith[2*k+1][jj], data_smith[2*k+1][jj -
                                                             1] = data_smith[2*k+1][jj-1], data_smith[2*k+1][jj]

    fig, axs = plt.subplots(figsize=(6, 4))
    for i in range(0, 7):
        axs.plot(data_smith[2*i], data_smith[2*i+1], 'k-')
        axs.text(data_smith[2*i][index_inflection[i]], data_smith[2*i+1]
                 [index_inflection[i]]+0.05, "%d%%" % (94-i), color='black')
    axs.set(xlabel=r"$\phi$ = V$_a$ / U",
            ylabel=r"$\psi$ = 2c$_p$$\Delta$T$_0$ / U$^2$")
    axs.spines['top'].set_visible(False)
    axs.spines['right'].set_visible(False)
    plt.xlim(0.3, 1.2)
    plt.ylim(1.5, 5.5)
    return fig


def info_blade():
    figure3, axs = plt.subplots(figsize=(5.5, 0.5))
    axs.text(0.0, 0.0, r"$\dfrac{dh_0}{dr} \: = \: V_a \: \dfrac{dV_a}{dr} \: + \: V_w \: \dfrac{dV_w}{dr} \: + \: \dfrac{V_w^2}{r}$",
             horizontalalignment='center', color='black', fontsize=14)
    plt.xlim(-1.0, 1.0)
    plt.ylim(-0.5, 1.5)
    plt.axis('off')
    canvas3 = FigureCanvasTkAgg(figure3, master=info_blade_Gui)
    canvas3.draw()
    canvas3.get_tk_widget().pack()
    plt.close(figure3)

    figure4, axs = plt.subplots(figsize=(5.5, 0.5))
    axs.text(-0.5, 0.0, r"$V_w r \: = \: constant$",
             horizontalalignment='center', color='black', fontsize=14)
    axs.text(+0.5, 0.0, r"$V_a \: = \: constant$",
             horizontalalignment='center', color='black', fontsize=14)
    plt.xlim(-1.0, 1.0)
    plt.axis('off')
    canvas4 = FigureCanvasTkAgg(figure4, master=info_blade_Gui)
    canvas4.draw()
    canvas4.get_tk_widget().pack()
    plt.close(figure4)

    figure5, axs = plt.subplots(figsize=(5.5, 0.5))
    axs.text(-0.5, 0.0, r"$\alpha_2 \: = \: constant$",
             horizontalalignment='center', color='black', fontsize=14)
    axs.text(+0.5, 0.0, r"$V_{w_2}\:r^{\:sin^2\:\alpha_2} \: = \: constant$",
             horizontalalignment='center', color='black', fontsize=14)
    plt.xlim(-1.0, 1.0)
    plt.axis('off')
    canvas5 = FigureCanvasTkAgg(figure5, master=info_blade_Gui)
    canvas5.draw()
    canvas5.get_tk_widget().pack()
    plt.close(figure5)


def info_stress():
    figure7, axs = plt.subplots(figsize=(5, 1.2))
    axs.text(0.0, +0.4, r"$(\sigma_{ct})_{max} \: = \: \dfrac{2}{3} \: \pi \: N^2 \: \rho_b \: (A_2+A_3)$",
             horizontalalignment='center', color='black', fontsize=14)
    axs.text(0.0, -0.4, r"$(\sigma_{gb})_{max} \: = \: \dfrac{\dot m \: (V_{w_2}+V_{w_3}) \: h_R}{2 \: z \: n_R \: c_R^3}$",
             horizontalalignment='center', color='black', fontsize=14)
    plt.ylim(-0.6, 0.6)
    plt.xlim(-1.0, 1.0)
    plt.axis('off')
    canvas7 = FigureCanvasTkAgg(figure7, master=info_stress_Gui)
    canvas7.draw()
    canvas7.get_tk_widget().pack()
    plt.close(figure7)

    figure8 = graph_zBn()
    canvas8 = FigureCanvasTkAgg(figure8, master=info_stress_Gui)
    canvas8.draw()
    canvas8.get_tk_widget().pack()
    plt.close(figure8)


def graph_zBn():
    global data_B, data_n
    fig, axs1 = plt.subplots(figsize=(6, 4))
    axs1.plot(data_B[0], data_B[1], 'b-')
    axs2 = axs1.twinx()
    axs2.plot(data_n[0], data_n[1], 'g-')
    axs1.set_xlim(40, 120)
    axs1.set_ylim(200, 1000)
    axs2.set_ylim(1.0, 2.0)
    axs1.set_xticklabels('{:.0f}°'.format(i) for i in axs1.get_xticks())
    axs1.set_xlabel("Blade camber angle\n" +
                    r"( $\beta$$_2$$_r$ + $\beta$$_3$$_r$ )", labelpad=10)
    axs1.set_ylabel("Coefficient B", labelpad=0)
    axs2.set_ylabel("Coefficient n", labelpad=10)
    axs1.spines['top'].set_visible(False)
    axs2.spines['top'].set_visible(False)
    B_leg = lines.Line2D([], [], color='blue', linestyle='-')
    n_leg = lines.Line2D([], [], color='green', linestyle='-')
    axs1.legend((B_leg, n_leg), ("B", "n"), loc='lower left')
    axs1.text(
        80, 900, r"z = $\dfrac{I_{xx}}{y_{max}}$ = $\dfrac{1}{B}$ $\left(10\:\dfrac{t}{c}\right)^n$", fontsize=14)
    fig.tight_layout()
    return fig


def ita_logo_plot():
    fig, axs = plt.subplots(figsize=(6, 2))
    axs.imshow(img.imread("ita_logo.png"))
    plt.axis('off')
    return fig


def read_input():  # lê os inputs
    global m, ws, T01, p01, PR0
    m = 20
    ws = 166460
    T01 = 1100
    p01 = 4e5
    PR0 = 1.873


def read_hypth():  # lê as hipóteses
    global N, U, alpha1, alpha3, lbdN, phi, phi2, phi3
    global read_hypth_flag, iteration_flag
    read_hypth_flag = True
    N = 250 # Hz
    U = 340
    alpha1 = 0
    alpha3 = 10
    if (iteration_flag == int(0)):
        lbdN = 0.05
    else:
        text10.set("%.4f" % lbdN)
    if (axial_velocity_flag.get() == 1):
        phi = 0.8
    else:
        phi2 = 0.7
        phi3 = 1.0
    velocity_diagrams_mean()


def velocity_diagrams_mean():  # define cálculos do triangulo de velocidades baseado no input de variação da vel axial
    global DT0, psi, dor, lbdN
    DT0 = ws/cp
    psi = (2.0*cp*DT0)/(U**2.0)

    P3m.U = U
    if (axial_velocity_flag.get() == 1):
        P3m.Va = U*phi
    else:
        P3m.Va = U*phi3
    P3m.alpha = alpha3
    P3m.Vu = P3m.Va*np.tan(np.radians(P3m.alpha))
    P3m.V = np.sqrt(P3m.Va**2.0 + P3m.Vu**2.0)
    P3m.Wa = P3m.Va
    P3m.Wu = P3m.Vu + U
    P3m.W = np.sqrt(P3m.Wa**2.0 + P3m.Wu**2.0)
    P3m.beta = np.degrees(np.arctan(P3m.Wu/P3m.Wa))

    P2m.U = U
    if (axial_velocity_flag.get() == 1):
        P2m.Va = U*phi
    else:
        P2m.Va = U*phi2
    P2m.Vu = U*psi/2.0 - P3m.Vu
    P2m.V = np.sqrt(P2m.Va**2.0 + P2m.Vu**2.0)
    P2m.alpha = np.degrees(np.arctan(P2m.Vu/P2m.Va))
    P2m.Wa = P2m.Va
    P2m.Wu = P2m.Vu - U
    P2m.W = np.sqrt(P2m.Wa**2.0 + P2m.Wu**2.0)
    P2m.beta = np.degrees(np.arctan(P2m.Wu/P2m.Wa))

    P1m.V = P3m.V
    P1m.alpha = alpha1
    P1m.Va = P1m.V*np.cos(np.radians(P1m.alpha))
    P1m.Vu = P1m.V*np.sin(np.radians(P1m.alpha))

    P1m.T0 = T01
    P1m.T = P1m.T0 - (P1m.V**2.0)/(2.0*cp)
    P1m.p0 = p01
    P1m.p = P1m.p0*(P1m.T/P1m.T0)**(gamma/(gamma-1.0))
    P1m.rho = (P1m.p*10**5.0)/(R*P1m.T)
    G1.A = m/(P1m.rho*P1m.Va)

    P2m.T0 = T01
    P2m.T = P2m.T0 - (P2m.V**2.0)/(2.0*cp)
    P2m.Ts = P2m.T - (lbdN*P2m.V**2.0)/(2.0*cp)
    prn = (T01/P2m.Ts)**(gamma/(gamma-1.0))
    prc = ((gamma+1.0)/2.0)**(gamma/(gamma-1.0))
    P2m.p = p01/prn
    P2m.p0 = P2m.p*(P2m.T0/P2m.T)**(gamma/(gamma-1.0))
    P2m.rho = (P2m.p*10**5.0)/(R*P2m.T)
    G2.A = m/(P2m.rho*P2m.Va)

    P3m.T0 = T01 - DT0
    P3m.T = P3m.T0 - (P3m.V**2.0)/(2.0*cp)
    P3m.p0 = p01/PR0
    P3m.p = P3m.p0*(P3m.T/P3m.T0)**(gamma/(gamma-1.0))
    P3m.Ts = P2m.T*(P3m.p/P2m.p)**((gamma-1.0)/gamma)
    P3m.rho = (P3m.p*10**5.0)/(R*P3m.T)
    G3.A = m/(P3m.rho*P3m.Va)

    dor = (P2m.T - P3m.T)/(P1m.T - P3m.T)

    G1.rm, G1.h, G1.rt, G1.rh, G1.rr = calculate_geometry(G1, U, N)
    G2.rm, G2.h, G2.rt, G2.rh, G2.rr = calculate_geometry(G2, U, N)
    G3.rm, G3.h, G3.rt, G3.rh, G3.rr = calculate_geometry(G3, U, N)

    P1m.Mach = P1m.V/np.sqrt(gamma*R*P1m.T)
    P2m.Mach = P2m.V/np.sqrt(gamma*R*P2m.T)
    P3m.Mach = P3m.V/np.sqrt(gamma*R*P3m.T)

    P2m.Machrel = P2m.W/np.sqrt(gamma*R*P2m.T)
    P3m.Machrel = P3m.W/np.sqrt(gamma*R*P3m.T)

    global visc_air_temp_curve
    visc_air_temp_curve = interp1d(
        data_visc_air_temp[0], data_visc_air_temp[1], kind='linear')

    P1m.mu = visc_air_temp_curve(P1m.T - 273.15) * 10 ** (-6.0)
    P2m.mu = visc_air_temp_curve(P2m.T - 273.15) * 10 ** (-6.0)
    P3m.mu = visc_air_temp_curve(P3m.T - 273.15) * 10 ** (-6.0)

    global column2_frame, veloc_frame, noteb_veloc
    global mean_frame, mean_frame_flag, blade_frame_flag
    if (mean_frame_flag == False):
        column2_frame = tk.Frame(design_window)
        column2_frame.grid(row=0, column=1, sticky='n')
        veloc_frame = tk.LabelFrame(
            column2_frame, text=" VELOCITY DIAGRAMS ", font=framelab_font)
        veloc_frame.grid(row=0, column=0, padx=10, pady=10)
        noteb_veloc = ttk.Notebook(veloc_frame)
        noteb_veloc.grid(row=0, column=0, sticky='n')
        mean_frame = tk.Frame(noteb_veloc)
        noteb_veloc.add(mean_frame, text='Mean')
        infobtn0 = tk.Button(veloc_frame, text="Learn More",
                             command=info_veloc, font=standard_font, width=10)
        infobtn0.grid(row=1, column=0, sticky='w'+'e', padx=5, pady=5)

    global iteration_flag
    if (iteration_flag == 0):
        figure1 = velocity_thermoprop_table(P1m, P2m, P3m)
        canvas1 = FigureCanvasTkAgg(figure1, master=mean_frame)
        canvas1.draw()
        canvas1.get_tk_widget().grid(row=0, column=0)
        plt.close(figure1)

        figure2 = velocity_diagrams_mean_plot(P1m, P2m, P3m)
        canvas2 = FigureCanvasTkAgg(figure2, master=mean_frame)
        canvas2.draw()
        canvas2.get_tk_widget().grid(row=1, column=0)
        plt.close(figure2)

    mean_frame_flag = True

    if (blade_frame_flag == False):
        blade_design()
    elif (root_tip_frame_flag == True):
        velocity_diagrams_root_tip()


def calculate_geometry(G, U_, N_):  # função para cálculo de geometria
    G.rm = U_/(2.0*np.pi*N_)
    G.h = G.A/(2.0*np.pi*G.rm)
    G.rt = G.rm + G.h/2.0
    G.rh = G.rm - G.h/2.0
    G.rr = G.rt/G.rh
    return G.rm, G.h, G.rt, G.rh, G.rr


# apresenta as informações termodinâmicas
def velocity_thermoprop_table(P1, P2, P3):
    fig, table_mean = plt.subplots(figsize=(5, 4.5))
    label_cols = ["Variables", "Units", "Station 1", "Station 2", "Station 3"]
    table_data = [["T$_0$", "K", "%.1f" % P1.T0, "%.1f" % P2.T0, "%.1f" % P3.T0],
                  ["T", "K", "%.1f" % P1.T, "%.1f" % P2.T, "%.1f" % P3.T],
                  ["p$_0$", "bar", "%.3f" % P1.p0, "%.3f" %
                      P2.p0, "%.3f" % P3.p0],
                  ["p", "bar", "%.3f" % P1.p, "%.3f" % P2.p, "%.3f" % P3.p],
                  [r"$\rho$", "kg/m$^3$", "%.3f" %
                      P1.rho, "%.3f" % P2.rho, "%.3f" % P3.rho],
                  ["V", "m/s", "%.1f" % P1.V, "%.1f" % P2.V, "%.1f" % P3.V],
                  ["V$_a$", "m/s", "%.1f" % P1.Va, "%.1f" %
                      P2.Va, "%.1f" % P3.Va],
                  ["V$_w$", "m/s", "%.1f" % P1.Vu, "%.1f" %
                      P2.Vu, "%.1f" % P3.Vu],
                  ["W", "m/s", "-", "%.1f" % P2.W, "%.1f" % P3.W],
                  ["W$_a$", "m/s", "-", "%.1f" % P2.Wa, "%.1f" % P3.Wa],
                  ["W$_w$", "m/s", "-", "%.1f" % P2.Wu, "%.1f" % P3.Wu],
                  ["U", "m/s", "-", "%.1f" % P2.U, "%.1f" % P3.U],
                  [r"$\alpha$", "°", "%.2f" % P1.alpha, "%.2f" %
                      P2.alpha, "%.2f" % P3.alpha],
                  [r"$\beta$", "°", "-", "%.2f" % P2.beta, "%.2f" % P3.beta],
                  ["Ma", "-", "%.3f" % P1.Mach, "%.3f" %
                      P2.Mach, "%.3f" % P3.Mach],
                  ["Ma$_r$$_e$$_l$", "-", "-", "%.3f" % P2.Machrel, "%.3f" % P3.Machrel]]
    table1 = table_mean.table(cellText=table_data, colLabels=label_cols, edges='horizontal',
                              cellLoc='center', rowLoc='center', colLoc='center', loc='center')
    table1.auto_set_font_size(False)
    table1.set_fontsize(10)
    table1.scale(0.9, 1.5)
    plt.axis('off')
    return fig


def velocity_diagrams_mean_plot(P1, P2, P3):  # figura triangulos de vel
    fig, velocity_diagrams = plt.subplots(figsize=(5, 3))

    origin_x = [0.0, P2.U, 0.0, P3.U, 0.0]
    origin_y = [0.0, 0.0, 0.0, 0.0, 0.0]
    direction_x = [P2.U, P2.Wu, P2.Vu, -P3.Wu, -P3.Vu]
    direction_y = [0.0, -P2.Wa, -P2.Va, -P3.Wa, -P3.Va]
    velocity_diagrams.quiver(origin_x, origin_y, direction_x, direction_y, color=[
                             'black', 'red', 'blue', 'red', 'blue'], angles='xy', scale_units='xy', scale=1)

    velocity_diagrams.text(P2.U/2, 10, '$U$', color='black')
    velocity_diagrams.text(P2.U+P2.Wu/3+10, -P2.Wa/3-10, '$W_2$', color='red')
    velocity_diagrams.text(2*P2.Vu/3-30, -2*P2.Va/3-10, '$V_2$', color='blue')
    velocity_diagrams.text(P3.U-2*P3.Wu/3+10, -2*P3.Wa /
                           3-10, '$W_3$', color='red')
    velocity_diagrams.text(-P3.Vu/3-30, -P3.Va/3-10, '$V_3$', color='blue')

    if (P2.Vu > P2.U):
        horizontal_line1 = lines.Line2D(
            [P2.U, P2.Vu], [0, 0], color='gray', linestyle='--')
        velocity_diagrams.add_line(horizontal_line1)
    vertical_line1 = lines.Line2D(
        [P2.Vu, P2.Vu], [-P2.Va, 0], color='gray', linestyle='--')
    velocity_diagrams.add_line(vertical_line1)
    if (P3.Vu > 0.0):
        horizontal_line2 = lines.Line2D(
            [-P3.Vu, 0], [0, 0], color='gray', linestyle='--')
        velocity_diagrams.add_line(horizontal_line2)
    vertical_line2 = lines.Line2D(
        [-P3.Vu, -P3.Vu], [-P3.Va, 0], color='gray', linestyle='--')
    velocity_diagrams.add_line(vertical_line2)

    if (P2.beta > 0.0):
        tt1, tt2 = 0.0, P2.beta
    else:
        tt1, tt2 = P2.beta, 0.0
    beta2_arc = ptc.Arc((P2.Vu, -P2.Va), 3*P2.Va/5, 3*P2.Va/5,
                        angle=90, theta1=tt1, theta2=tt2, color='red')
    velocity_diagrams.add_patch(beta2_arc)
    if (P2.alpha > 0.0):
        tt1, tt2 = 0.0, P2.alpha
    else:
        tt1, tt2 = P2.alpha, 0.0
    alpha2_arc = ptc.Arc((P2.Vu, -P2.Va), 2*P2.Va/5, 2*P2.Va/5,
                         angle=90, theta1=tt1, theta2=tt2, color='blue')
    velocity_diagrams.add_patch(alpha2_arc)
    if (P3.beta > 0.0):
        tt1, tt2 = -P3.beta, 0.0
    else:
        tt1, tt2 = 0.0, -P3.beta
    beta3_arc = ptc.Arc((-P3.Vu, -P3.Va), 2*P3.Va/5, 2 *
                        P3.Va/5, angle=90, theta1=tt1, theta2=tt2, color='red')
    velocity_diagrams.add_patch(beta3_arc)
    if (P3.alpha > 0.0):
        tt1, tt2 = -P3.alpha, 0.0
    else:
        tt1, tt2 = 0.0, -P3.alpha
    alpha3_arc = ptc.Arc((-P3.Vu, -P3.Va), 3*P3.Va/5, 3 *
                         P3.Va/5, angle=90, theta1=tt1, theta2=tt2, color='blue')
    velocity_diagrams.add_patch(alpha3_arc)

    velocity_diagrams.text(P2.Vu-(7*P2.Va/20)*np.sin(np.radians(P2.beta/2+6)), -P2.Va+(
        7*P2.Va/20)*np.cos(np.radians(P2.beta/2+6)), r'$\beta_2$', color='red')
    velocity_diagrams.text(P2.Vu-(3*P2.Va/10)*np.sin(np.radians((P2.alpha+P2.beta)/2+6)), -P2.Va+(
        3*P2.Va/10)*np.cos(np.radians((P2.alpha+P2.beta)/2+6)), r'$\alpha_2$', color='blue')
    velocity_diagrams.text(-P3.Vu+(5*P3.Va/20)*np.sin(np.radians((P3.alpha+P3.beta)/2-6)), -P3.Va+(
        5*P3.Va/20)*np.cos(np.radians((P3.alpha+P3.beta)/2-6)), r'$\beta_3$', color='red')
    velocity_diagrams.text(-P3.Vu+(7*P3.Va/20)*np.sin(np.radians(P3.alpha/2-6)), -P3.Va+(
        7*P3.Va/20)*np.cos(np.radians(P3.alpha/2-6)), r'$\alpha_3$', color='blue')

    if (axial_velocity_flag.get() == 1):
        plt.title(label=r"$\phi$ = %.3f ; $\psi$ = %.3f ; $\Lambda$ = %.3f" % (
            phi, psi, dor), loc='center', fontsize=12)
    else:
        plt.title(label=r"$\phi_2$ = %.3f ; $\phi_3$ = %.3f ; $\psi$ = %.3f ; $\Lambda$ = %.3f" % (
            phi2, phi3, psi, dor), loc='center', fontsize=12)
    plt.axis('equal')
    plt.xlim(np.min([-P3.Vu, 0.0])-25, np.max([P2.Vu, P2.U])+25)
    plt.ylim(-np.max([P2.Va, P3.Va])-25, 50)
    plt.axis('off')
    return fig


def blade_design():  # define parâmetros da subjanela de blade design
    global blade_frame
    blade_frame = tk.LabelFrame(
        column1_frame, text=" BLADE DESIGN METHOD ", font=framelab_font)
    blade_frame.grid(row=2, column=0, padx=10, pady=10)
    blade_frame.grid_columnconfigure(0, minsize=150)
    blade_frame.grid_columnconfigure(1, minsize=10)
    blade_frame.grid_columnconfigure(2, minsize=50)
    blade_frame.grid_columnconfigure(3, minsize=80)
    blade_frame.grid_columnconfigure(4, minsize=60)
    global blade_design_flag, text14
    blade_design_flag = tk.IntVar(value=1)
    rdbpar1 = tk.Radiobutton(
        blade_frame, text="Free Vortex Design", variable=blade_design_flag, value=1)
    rdbpar1.grid(row=0, column=0, columnspan=5, sticky='w')
    rdbpar2 = tk.Radiobutton(
        blade_frame, text="Constant Nozzle Angle Design", variable=blade_design_flag, value=2)
    rdbpar2.grid(row=1, column=0, columnspan=5, sticky='w')
    rdbpar3 = tk.Radiobutton(
        blade_frame, text="First Power Design", variable=blade_design_flag, value=3)
    rdbpar3.grid(row=2, column=0, columnspan=5, sticky='w')

    nxtbtn3 = tk.Button(blade_frame, text="Next",
                        command=velocity_diagrams_root_tip, font=standard_font, width=10)
    nxtbtn3.grid(row=3, column=0, columnspan=5, sticky='e', padx=5, pady=5)

    infobtn3 = tk.Button(blade_frame, text="Info",
                         command=info_blade, font=standard_font, width=10)
    infobtn3.grid(row=3, column=0, columnspan=5, sticky='w', padx=5, pady=5)

    global blade_frame_flag
    blade_frame_flag = True


def velocity_diagrams_root_tip():
    global P1h, P2h, P3h
    P1h = point()
    P2h = point()
    P3h = point()

    global P1t, P2t, P3t
    P1t = point()
    P2t = point()
    P3t = point()

    P1h, P2h, P3h = calculate_tip_hub(
        G2.rh, G3.rh, G2.rm, G3.rm, P1m, P2m, P3m)
    P1t, P2t, P3t = calculate_tip_hub(
        G2.rt, G3.rt, G2.rm, G3.rm, P1m, P2m, P3m)

    global column2_frame, choke_frame, root_frame, tip_frame, root_tip_frame_flag
    if (root_tip_frame_flag == False):
        choke_frame = tk.LabelFrame(
            column2_frame, text=" CHOKE CONDITION ", font=framelab_font)
        choke_frame.grid(row=1, column=0, padx=10, pady=10)
        choke_frame.grid_columnconfigure(0, minsize=362)
        root_frame = tk.Frame(noteb_veloc)
        noteb_veloc.add(root_frame, text='Root')
        tip_frame = tk.Frame(noteb_veloc)
        noteb_veloc.add(tip_frame, text='Tip')

    global iteration_flag
    if (iteration_flag == 0):
        figure3 = velocity_thermoprop_table(P1h, P2h, P3h)
        canvas3 = FigureCanvasTkAgg(figure3, master=root_frame)
        canvas3.draw()
        canvas3.get_tk_widget().grid(row=0, column=0)
        plt.close(figure3)

        figure4 = velocity_diagrams_root_tip_plot(P1h, P2h, P3h)
        canvas4 = FigureCanvasTkAgg(figure4, master=root_frame)
        canvas4.draw()
        canvas4.get_tk_widget().grid(row=1, column=0)
        plt.close(figure4)

        figure5 = velocity_thermoprop_table(P1t, P2t, P3t)
        canvas5 = FigureCanvasTkAgg(figure5, master=tip_frame)
        canvas5.draw()
        canvas5.get_tk_widget().grid(row=0, column=0)
        plt.close(figure5)

        figure6 = velocity_diagrams_root_tip_plot(P1t, P2t, P3t)
        canvas6 = FigureCanvasTkAgg(figure6, master=tip_frame)
        canvas6.draw()
        canvas6.get_tk_widget().grid(row=1, column=0)
        plt.close(figure6)

        Mach_numbers = [P1h.Mach, P1m.Mach, P1t.Mach,
                        P2h.Mach, P2m.Mach, P3m.Mach,
                        P2h.Machrel, P2m.Machrel, P2t.Machrel,
                        P3h.Machrel, P3m.Machrel, P3t.Machrel]
        figure0 = graph_Mach(Mach_numbers)
        canvas0 = FigureCanvasTkAgg(figure0, master=choke_frame)
        canvas0.draw()
        canvas0.get_tk_widget().grid(row=0, column=0)
        plt.close(figure0)

    root_tip_frame_flag = True

    pitch_chord_nblades_p1()


def graph_Mach(Mach_numbers):
    fig, axs = plt.subplots(figsize=(5, 2.8))
    locations = [0.2, 0.4, 0.6, 0.9, 1.1, 1.3, 1.6, 1.8, 2.0, 2.3, 2.5, 2.7]
    colors = ['xkcd:green', 'xkcd:red', 'xkcd:blue',
              'xkcd:green', 'xkcd:red', 'xkcd:blue',
              'xkcd:green', 'xkcd:red', 'xkcd:blue',
              'xkcd:green', 'xkcd:red', 'xkcd:blue']
    labels = ["", "Station 1\n(Absolute)", "", "", "Station 1\n(Absolute)",
              "", "", "Station 2\n(Relative)", "", "", "Station 3\n(Relative)", "", ]
    axs.bar(locations, Mach_numbers, color=colors,
            width=0.2, tick_label=labels)
    axs.set(ylabel="Mach Number")
    root_leg = lines.Line2D([], [], color='xkcd:green', linestyle='-')
    mean_leg = lines.Line2D([], [], color='xkcd:red', linestyle='-')
    tip_leg = lines.Line2D([], [], color='xkcd:blue', linestyle='-')
    axs.legend((root_leg, mean_leg, tip_leg),
               ("Root", "Mean", "Tip"), loc='upper left')
    axs.spines['top'].set_visible(False)
    axs.spines['right'].set_visible(False)
    axs.set_ylim(0, 1.2)
    axs.set_xlim(0, 2.9)
    fig.tight_layout()
    return fig


def calculate_tip_hub(r2, r3, r2m, r3m, P1m_, P2m_, P3m_):
    P1 = point()
    P2 = point()
    P3 = point()

    P1 = P1m_

    if (blade_design_flag.get() == 1):
        P2.alpha = np.degrees(
            np.arctan((r2m/r2)*np.tan(np.radians(P2m_.alpha))))
        P2.beta = np.degrees(
            np.arctan((r2m/r2)*np.tan(np.radians(P2m_.alpha))-(r2/r2m)*(P2m_.U/P2m_.Va)))
        P3.alpha = np.degrees(
            np.arctan((r3m/r3)*np.tan(np.radians(P3m_.alpha))))
        P3.beta = np.degrees(
            np.arctan((r3m/r3)*np.tan(np.radians(P3m_.alpha))+(r3/r3m)*(P3m_.U/P3m_.Va)))

        P2.U = (r2/r2m)*P2m_.U
        P2.Va = P2m_.Va
        P2.Vu = P2.Va*np.tan(np.radians(P2.alpha))
        P2.V = np.sqrt(P2.Va**2.0 + P2.Vu**2.0)
        P2.Wa = P2.Va
        P2.Wu = P2.Wa*np.tan(np.radians(P2.beta))
        P2.W = np.sqrt(P2.Wa**2.0 + P2.Wu**2.0)

        P3.U = (r3/r3m)*P3m_.U
        P3.Va = P3m_.Va
        P3.Vu = P3.Va*np.tan(np.radians(P3.alpha))
        P3.V = np.sqrt(P3.Va**2.0 + P3.Vu**2.0)
        P3.Wa = P3.Va
        P3.Wu = P3.Wa*np.tan(np.radians(P3.beta))
        P3.W = np.sqrt(P3.Wa**2.0 + P3.Wu**2.0)

        P2.T0 = P1.T0
        P2.T = P2.T0 - (P2.V**2.0)/(2.0*cp)
        P2.Ts = P2.T - (lbdN*P2.V**2.0)/(2.0*cp)
        prn = (P1.T0/P2.Ts)**(gamma/(gamma-1.0))
        P2.p = P1.p0/prn
        P2.p0 = P2.p*(P2.T0/P2.T)**(gamma/(gamma-1.0))
        P2.rho = (P2.p*10**5.0)/(R*P2.T)

        P3.T0 = P2.T0 - DT0
        P3.T = P3.T0 - (P3.V**2.0)/(2.0*cp)
        P3.p0 = P1.p0/PR0
        P3.p = P3.p0*(P3.T/P3.T0)**(gamma/(gamma-1.0))
        P3.Ts = P2.T*(P3.p/P2.p)**((gamma-1.0)/gamma)
        P3.rho = (P3.p*10**5.0)/(R*P3.T)
    elif (blade_design_flag.get() == 2):
        P2.U = P2m_.U*(r2/r2m)
        P2.alpha = P2m_.alpha
        P2.Va = P2m_.Va*(r2m/r2)**((np.sin(np.radians(P2.alpha)))**2.0)
        P2.Vu = P2m_.Vu*(r2m/r2)**((np.sin(np.radians(P2.alpha)))**2.0)
        P2.V = np.sqrt(P2.Va**2.0 + P2.Vu**2.0)
        P2.Wa = P2.Va
        P2.Wu = P2.Vu - P2.U
        P2.W = np.sqrt(P2.Wa**2.0 + P2.Wu**2.0)
        P2.beta = np.degrees(np.arctan(P2.Wu/P2.Wa))

        P2.T0 = P1.T0
        P2.T = P2.T0 - (P2.V**2.0)/(2.0*cp)
        P2.Ts = P2.T - (lbdN*P2.V**2.0)/(2.0*cp)
        prn = (P1.T0/P2.Ts)**(gamma/(gamma-1.0))
        P2.p = P1.p0/prn
        P2.p0 = P2.p*(P2.T0/P2.T)**(gamma/(gamma-1.0))
        P2.rho = (P2.p*10**5.0)/(R*P2.T)

        P3.U = P3m_.U*(r3/r3m)
        Vanew = P3m_.Va
        P3.T0 = P2.T0 - DT0
        P3.p0 = P1.p0/PR0
        P3.Vu = ((P2m_.U*P2m_.Vu + P3m_.U*P3m_.Vu) - P2.U*P2.Vu)/P3.U
        while (np.abs(P3.Va - Vanew) > 0.001):
            P3.Va = Vanew
            P3.V = np.sqrt(P3.Va**2.0 + P3.Vu**2.0)
            P3.T = P3.T0 - (P3.V**2.0)/(2.0*cp)
            P3.p = P3.p0*(P3.T/P3.T0)**(gamma/(gamma-1.0))
            P3.Ts = P2.T*(P3.p/P2.p)**((gamma-1.0)/gamma)
            P3.rho = (P3.p*10**5.0)/(R*P3.T)
            vanew = (P2.rho*P2.Va*r2)/(P3.rho*r3)
        P3.alpha = np.degrees(np.arctan(P3.Vu/P3.Va))
        P3.Wa = P3.Va
        P3.Wu = P3.Vu + P2.U
        P3.W = np.sqrt(P3.Wa**2.0 + P3.Wu**2.0)
        P3.beta = np.degrees(np.arctan(P3.Wu/P3.Wa))

    P2.Mach = P2.V/np.sqrt(gamma*R*P2.T)
    P3.Mach = P3.V/np.sqrt(gamma*R*P3.T)

    P2.Machrel = P2.W/np.sqrt(gamma*R*P2.T)
    P3.Machrel = P3.W/np.sqrt(gamma*R*P3.T)

    global visc_air_temp_curve
    P1.mu = visc_air_temp_curve(P1.T-273.15)*10**(-6.0)
    P2.mu = visc_air_temp_curve(P2.T-273.15)*10**(-6.0)
    P3.mu = visc_air_temp_curve(P3.T-273.15)*10**(-6.0)

    return P1, P2, P3


def velocity_diagrams_root_tip_plot(P1, P2, P3):
    fig, velocity_diagrams = plt.subplots(figsize=(5, 3))

    offset_x = P3.U + 25
    origin_x = [offset_x, offset_x + P2.U, offset_x, 0.0, P3.U, 0.0]
    origin_y = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    direction_x = [P2.U, P2.Wu, P2.Vu, P3.U, -P3.Wu, -P3.Vu]
    direction_y = [0.0, -P2.Wa, -P2.Va, 0.0, -P3.Wa, -P3.Va]
    velocity_diagrams.quiver(origin_x, origin_y, direction_x, direction_y, color=[
                             'black', 'red', 'blue', 'black', 'red', 'blue'], angles='xy', scale_units='xy', scale=1)

    velocity_diagrams.text(offset_x + P2.U/2-20, 10, '$U_2$', color='black')
    velocity_diagrams.text(offset_x + P2.U+P2.Wu/3+20, -
                           P2.Wa/3-20, '$W_2$', color='red')
    velocity_diagrams.text(offset_x + 2*P2.Vu/3-40, -
                           2*P2.Va/3-20, '$V_2$', color='blue')
    velocity_diagrams.text(P3.U/2-20, 10, '$U_3$', color='black')
    velocity_diagrams.text(P3.U-2*P3.Wu/3+20, -2*P3.Wa /
                           3-20, '$W_3$', color='red')
    velocity_diagrams.text(-P3.Vu/3-40, -P3.Va/3-20, '$V_3$', color='blue')

    if (P2.beta > 0.0):
        horizontal_line1 = lines.Line2D(
            [offset_x + P2.U, offset_x + P2.Vu], [0, 0], color='gray', linestyle='--')
        velocity_diagrams.add_line(horizontal_line1)
    vertical_line1 = lines.Line2D(
        [offset_x + P2.Vu, offset_x + P2.Vu], [-P2.Va, 0], color='gray', linestyle='--')
    velocity_diagrams.add_line(vertical_line1)
    if (P3.alpha > 0.0):
        horizontal_line2 = lines.Line2D(
            [-P3.Vu, 0], [0, 0], color='gray', linestyle='--')
        velocity_diagrams.add_line(horizontal_line2)
    vertical_line2 = lines.Line2D(
        [-P3.Vu, -P3.Vu], [-P3.Va, 0], color='gray', linestyle='--')
    velocity_diagrams.add_line(vertical_line2)

    if (P2.beta > 0.0):
        tt1, tt2 = 0.0, P2.beta
    else:
        tt1, tt2 = P2.beta, 0.0
    beta2_arc = ptc.Arc((offset_x + P2.Vu, -P2.Va), 3*P2.Va/5,
                        3*P2.Va/5, angle=90, theta1=tt1, theta2=tt2, color='red')
    velocity_diagrams.add_patch(beta2_arc)
    if (P2.alpha > 0.0):
        tt1, tt2 = 0.0, P2.alpha
    else:
        tt1, tt2 = P2.alpha, 0.0
    alpha2_arc = ptc.Arc((offset_x + P2.Vu, -P2.Va), 2*P2.Va/5,
                         2*P2.Va/5, angle=90, theta1=tt1, theta2=tt2, color='blue')
    velocity_diagrams.add_patch(alpha2_arc)
    if (P3.beta > 0.0):
        tt1, tt2 = -P3.beta, 0.0
    else:
        tt1, tt2 = 0.0, -P3.beta
    beta3_arc = ptc.Arc((-P3.Vu, -P3.Va), 2*P3.Va/5, 2 *
                        P3.Va/5, angle=90, theta1=tt1, theta2=tt2, color='red')
    velocity_diagrams.add_patch(beta3_arc)
    if (P3.alpha > 0.0):
        tt1, tt2 = -P3.alpha, 0.0
    else:
        tt1, tt2 = 0.0, -P3.alpha
    alpha3_arc = ptc.Arc((-P3.Vu, -P3.Va), 3*P3.Va/5, 3 *
                         P3.Va/5, angle=90, theta1=tt1, theta2=tt2, color='blue')
    velocity_diagrams.add_patch(alpha3_arc)

    velocity_diagrams.text(offset_x + P2.Vu-(7*P2.Va/20)*np.sin(np.radians(P2.beta/2+9)), -
                           P2.Va+(7*P2.Va/20)*np.cos(np.radians(P2.beta/2+9)), r'$\beta_2$', color='red')
    velocity_diagrams.text(offset_x + P2.Vu-(7*P2.Va/20)*np.sin(np.radians((P2.alpha+P2.beta)/2+9)), -
                           P2.Va+(3*P2.Va/10)*np.cos(np.radians((P2.alpha+P2.beta)/2+9)), r'$\alpha_2$', color='blue')
    velocity_diagrams.text(-P3.Vu+(5*P3.Va/20)*np.sin(np.radians((P3.alpha+P3.beta)/2-9)), -P3.Va+(
        5*P3.Va/20)*np.cos(np.radians((P3.alpha+P3.beta)/2-9)), r'$\beta_3$', color='red')
    velocity_diagrams.text(-P3.Vu+(7*P3.Va/20)*np.sin(np.radians(P3.alpha/2-9)), -P3.Va+(
        7*P3.Va/20)*np.cos(np.radians(P3.alpha/2-9)), r'$\alpha_3$', color='blue')

    dor_ = (P2.T - P3.T)/(P1.T - P3.T)

    plt.title(label=r"$\Lambda$ = %.3f" % dor_, loc='center', fontsize=12)
    plt.axis('equal')
    plt.xlim(np.min([-P3.Vu, 0.0])-25, offset_x+np.max([P2.Vu, P2.U])+25)
    plt.ylim(-np.max([P2.Va, P3.Va])-25, 50)
    plt.axis('off')
    return fig


def pitch_chord_nblades_p1():  # calculos pitch/chord height/chord
    global nozzle, rotor
    nozzle = blades()
    rotor = blades()

    nozzle.h, nozzle.rm, nozzle.rt, nozzle.rh, nozzle.rr = calculate_blades(
        nozzle, G1, G2)
    rotor.h, rotor.rm, rotor.rt, rotor.rh, rotor.rr = calculate_blades(
        rotor, G2, G3)

    nozzle.sc = calculate_optsc(P1m.alpha, P2m.alpha)
    rotor.sc = calculate_optsc(P2m.beta, P3m.beta)

    figure7 = graph_optsc()
    canvas7 = FigureCanvasTkAgg(figure7, master=pcnbl_frame)
    canvas7.draw()
    canvas7.get_tk_widget().grid(row=0, column=0, columnspan=4)
    plt.close(figure7)

    pitch_chord_nblades_p2()


def pitch_chord_nblades_p2():
    nozzle.hc = 3.0
    rotor.hc = 3.0
    rotor.hw = 3.0
    rotor.sbw = 0.25

    nozzle.c = nozzle.h/nozzle.hc
    nozzle.s = nozzle.sc*nozzle.c
    nozzle.n = 2.0*np.pi*nozzle.rm/nozzle.s
    if (np.abs(2.0*np.floor(nozzle.n/2.0) - nozzle.n) < np.abs(2.0*np.ceil(nozzle.n/2.0) - nozzle.n)):
        nozzle.n = int(2.0*np.floor(nozzle.n/2.0))
    else:
        nozzle.n = int(2.0*np.ceil(nozzle.n/2.0))
    nozzle.s = 2.0*np.pi*nozzle.rm/nozzle.n
    nozzle.c = nozzle.s/nozzle.sc
    nozzle.hc = nozzle.h/nozzle.c

    rotor.c = rotor.h/rotor.hc
    rotor.s = rotor.sc*rotor.c
    rotor.n = 2.0*np.pi*rotor.rm/rotor.s
    rotor.n = nearest_prime(rotor.n)
    rotor.s = 2.0*np.pi*rotor.rm/rotor.n
    rotor.c = rotor.s/rotor.sc
    rotor.hc = rotor.h/rotor.c

    rotor.w = rotor.h/rotor.hw
    sbnr = rotor.sbw*rotor.w
    flare_angle = 2.0 * \
        np.degrees(np.arctan(((G3.h - G2.h)/2.0)/(sbnr + rotor.w)))
    nozzle.w = ((G2.h - G1.h)/2.0) / \
        np.tan(np.radians(flare_angle/2.0)) - sbnr
    nozzle.hw = nozzle.h/nozzle.w

    nozzle.in_angle = P1m.alpha
    nozzle.out_angle = P2m.alpha
    rotor.in_angle = P2m.beta
    rotor.out_angle = P3m.beta

    global sktch_frame, column4_frame, sktch_frame_flag
    if (sktch_frame_flag == False):
        column4_frame = tk.Frame(design_window)
        column4_frame.grid(row=0, column=3, sticky='n')
        sktch_frame = tk.LabelFrame(
            column4_frame, text=" STAGE SKETCH ", font=framelab_font)
        sktch_frame.grid(row=0, column=0, padx=10, pady=10)

    figure8 = sketch_stage(
        nozzle.h, rotor.h, nozzle.w, rotor.w, sbnr, flare_angle)
    canvas8 = FigureCanvasTkAgg(figure8, master=sktch_frame)
    canvas8.draw()
    canvas8.get_tk_widget().grid(row=7, column=0, columnspan=4)
    plt.close(figure8)

    figure9, blade_parameters = plt.subplots(figsize=(5, 3))
    label_cols = ["Variables", "Units", "Nozzle", "Rotor"]
    table_data = [["r$_m$", "m", "%.4f" % nozzle.rm, "%.4f" % rotor.rm],
                  ["h", "m", "%.4f" % nozzle.h, "%.4f" % rotor.h],
                  ["w", "m", "%.4f" % nozzle.w, "%.4f" % rotor.w],
                  ["s$_N$$_R$", "m", "%.4f" % sbnr, "%.4f" % sbnr],
                  [r"$\gamma$", "°", "%.2f" %
                      flare_angle, "%.2f" % flare_angle],
                  ["c", "m", "%.4f" % nozzle.c, "%.4f" % rotor.c],
                  ["s", "m", "%.4f" % nozzle.s, "%.4f" % rotor.s],
                  ["(h/c)", "-", "%.3f" %
                   nozzle.hc, "%.3f" % rotor.hc],
                  ["(s/c)", "-", "%.3f" %
                   nozzle.sc, "%.3f" % rotor.sc],
                  ["n", "-", "%d" % nozzle.n, "%d" % rotor.n]]
    table1 = blade_parameters.table(cellText=table_data, colLabels=label_cols, edges='horizontal',
                                    cellLoc='center', rowLoc='center', colLoc='center', loc='center')
    table1.auto_set_font_size(False)
    table1.set_fontsize(10)
    table1.scale(0.7, 1.5)
    plt.axis('off')
    canvas9 = FigureCanvasTkAgg(figure9, master=sktch_frame)
    canvas9.draw()
    canvas9.get_tk_widget().grid(row=8, column=0, columnspan=4)
    plt.close(figure9)

    global blpar_frame_flag
    if (blpar_frame_flag == False):
        shroud_flag = 1
    elif (stress_frame_flag == True):
        stress_analysis()


def calculate_blades(blade, Gin, Gout):
    blade.h = (Gin.h + Gout.h)/2.0
    blade.rm = (Gin.rm + Gout.rm)/2.0
    blade.rt = (Gin.rt + Gout.rt)/2.0
    blade.rh = (Gin.rh + Gout.rh)/2.0
    blade.rr = blade.rt/blade.rh
    return blade.h, blade.rm, blade.rt, blade.rh, blade.rr


def calculate_optsc(inlet, outlet):
    in_lower = int(np.floor(inlet/10.0))
    in_upper = int(np.floor(inlet/10.0)) + 1
    interp1 = np.poly1d(np.polyfit(
        data_sc[2*in_lower], data_sc[2*in_lower+1], 3))
    interp2 = np.poly1d(np.polyfit(
        data_sc[2*in_upper], data_sc[2*in_upper+1], 3))
    optimum_sc = interp1(outlet) + (interp2(outlet)-interp1(outlet)) * \
        (inlet-in_lower*10)/(in_upper*10-in_lower*10)
    return optimum_sc


def graph_optsc():
    fig, axs = plt.subplots(figsize=(6, 4))
    axs.text(37.5, 1.0, "Relative inlet gas\n" +
             r"angle ($\alpha_1$ or $\beta_2$)", color='black')
    for i in range(0, 7):
        axs.plot(data_sc[2*i], data_sc[2*i+1], 'k-')
        axs.text(data_sc[2*i][0]-2.5, data_sc[2*i+1]
                 [0]-0.01, "%d°" % (i*10), color='black')
    axs.set(xlabel=r"Relative efflux gas angle ($\alpha_2$ or $\beta_3$)",
            ylabel="'Optimum' s/c")
    axs.spines['top'].set_visible(False)
    axs.spines['right'].set_visible(False)
    plt.xlim(35, 75)
    plt.ylim(0.5, 1.1)
    global P2m, P3m
    axs.scatter(P2m.alpha, nozzle.sc, c='b', marker='o', alpha=0.5)
    nozzle_v = lines.Line2D([P2m.alpha, P2m.alpha], [
                            0.5, nozzle.sc], color='blue', linestyle='--', alpha=0.5)
    axs.add_line(nozzle_v)
    nozzle_h = lines.Line2D([35, P2m.alpha], [
                            nozzle.sc, nozzle.sc], color='blue', linestyle='--', alpha=0.5)
    axs.add_line(nozzle_h)
    axs.scatter(P3m.beta, rotor.sc, c='g', marker='o', alpha=0.5)
    rotor_v = lines.Line2D([P3m.beta, P3m.beta], [
                           0.5, rotor.sc], color='green', linestyle='--', alpha=0.5)
    axs.add_line(rotor_v)
    rotor_h = lines.Line2D(
        [35, P3m.beta], [rotor.sc, rotor.sc], color='green', linestyle='--', alpha=0.5)
    axs.add_line(rotor_h)

    axs.set_xticklabels('{:.0f}°'.format(i) for i in axs.get_xticks())
    nozzle_leg = lines.Line2D([], [], color='blue',
                              linestyle='--', alpha=0.5, marker='o')
    rotor_leg = lines.Line2D([], [], color='green',
                             linestyle='--', alpha=0.5, marker='o')
    axs.legend((nozzle_leg, rotor_leg), ("(s/c)$_N$ = %.3f" %
                                         nozzle.sc, "(s/c)$_R$ = %.3f" % rotor.sc), loc='upper right')
    return fig


def nearest_prime(value):
    primes = []
    for num in range(2, 200):
        flag = False
        for prime in primes:
            if (num % prime == 0):
                flag = True
                break
        if (flag == False):
            primes.append(num)
    index = np.argmin(np.abs(np.array(primes) - value))
    return primes[index]


def sketch_stage(hN, hR, wN, wR, sNR, flr):
    fig, sketch = plt.subplots(figsize=(5, 3.5))

    xpos1 = 0.0
    ypos1 = hN/2.0 - (sNR + wN/2.0)*np.tan(np.radians(flr/2.0))
    xpos2 = -ypos1/np.tan(np.radians(flr/2.0))
    ypos2 = 0.0
    auxiliary_line1 = lines.Line2D(
        [xpos1, xpos2], [ypos1, ypos2], color='gray', linestyle='--')
    sketch.add_line(auxiliary_line1)
    auxiliary_line2 = lines.Line2D(
        [xpos1, xpos2], [-ypos1, -ypos2], color='gray', linestyle='--')
    sketch.add_line(auxiliary_line2)
    xpos3 = 3.0*sNR + wN + wR
    ypos3 = 0.0
    auxiliary_line3 = lines.Line2D(
        [xpos2, xpos3], [ypos2, ypos3], color='gray', linestyle='-.')
    sketch.add_line(auxiliary_line3)
    flare_arc = ptc.Arc((xpos2, ypos2), np.sqrt(xpos2**2.0 + ypos1**2.0)/2.0, np.sqrt(
        xpos2**2.0 + ypos1**2.0)/2.0, angle=0, theta1=-flr/2.0, theta2=flr/2.0, color='black')
    sketch.add_patch(flare_arc)
    sketch.text(-5.0*np.sqrt(xpos2**2.0 + ypos1**2.0) /
                8.0, 0.0, r"$\gamma$", color='black')

    xpos1 = 0.0
    ypos1 = hN/2.0 - (sNR + wN/2.0)*np.tan(np.radians(flr/2.0))
    xpos2 = 3.0*sNR + wN + wR
    ypos2 = hR/2.0 + (sNR + wR/2.0)*np.tan(np.radians(flr/2.0))
    wall_line1 = lines.Line2D(
        [xpos1, xpos2], [ypos1, ypos2], color='black', linestyle='-')
    sketch.add_line(wall_line1)
    wall_line2 = lines.Line2D(
        [xpos1, xpos2], [-ypos1, -ypos2], color='black', linestyle='-')
    sketch.add_line(wall_line2)

    xpos1 = sNR
    ypos1 = hN/2.0 - (wN/2.0)*np.tan(np.radians(flr/2.0))
    nozzle_line1 = lines.Line2D(
        [xpos1, xpos1], [ypos1, -ypos1], color='black', linestyle='-')
    sketch.add_line(nozzle_line1)
    xpos2 = sNR + wN
    ypos2 = hN/2.0 + (wN/2.0)*np.tan(np.radians(flr/2.0))
    nozzle_line2 = lines.Line2D(
        [xpos2, xpos2], [ypos2, -ypos2], color='black', linestyle='-')
    sketch.add_line(nozzle_line2)
    xpos1 = sNR + wN/2.0
    ypos1 = hN/6.0
    sketch.text(xpos1-0.002, ypos1, "N", color='black')

    xpos1 = 2.0*sNR + wN
    ypos1 = hR/2.0 - (wR/2.0)*np.tan(np.radians(flr/2.0))
    rotor_line1 = lines.Line2D(
        [xpos1, xpos1], [ypos1-0.02*hR, -ypos1], color='black', linestyle='-')
    sketch.add_line(rotor_line1)
    xpos2 = 2.0*sNR + wN + wR
    ypos2 = hR/2.0 + (wR/2.0)*np.tan(np.radians(flr/2.0))
    rotor_line2 = lines.Line2D(
        [xpos2, xpos2], [ypos2-0.02*hR, -ypos2], color='black', linestyle='-')
    sketch.add_line(rotor_line2)
    rotor_line3 = lines.Line2D(
        [xpos1, xpos2], [ypos1-0.02*hR, ypos2-0.02*hR], color='black', linestyle='-')
    sketch.add_line(rotor_line3)
    xpos1 = 2.0*sNR + wN + wR/2.0
    ypos1 = hN/6.0
    sketch.text(xpos1-0.002, ypos1, "R", color='black')

    xpos1 = sNR/2.0
    ypos1 = G1.h/2.0
    stations_line1 = lines.Line2D([xpos1, xpos1], [
                                  ypos1+0.2*G1.h, -(ypos1 + + 0.2*G1.h)], color='blue', linestyle='--')
    sketch.add_line(stations_line1)
    sketch.text(xpos1-0.002, ypos1+0.3*G1.h, "1", color='blue')
    xpos1 = 3.0*sNR/2.0 + wN
    ypos1 = G2.h/2.0
    stations_line2 = lines.Line2D(
        [xpos1, xpos1], [ypos1+0.2*G1.h, -(ypos1+0.2*G1.h)], color='blue', linestyle='--')
    sketch.add_line(stations_line2)
    sketch.text(xpos1-0.002, ypos1+0.3*G1.h, "2", color='blue')
    xpos1 = 5.0*sNR/2.0 + wN + wR
    ypos1 = G3.h/2.0
    stations_line3 = lines.Line2D(
        [xpos1, xpos1], [ypos1+0.2*G1.h, -(ypos1+0.2*G1.h)], color='blue', linestyle='--')
    sketch.add_line(stations_line3)
    sketch.text(xpos1-0.002, ypos1+0.3*G1.h, "3", color='blue')

    cline1 = lines.Line2D([1.5*sNR + wN/2.0, 6.0*sNR + wN + wR],
                          [hN/2.0, hN/2.0], color='gray', linestyle='-')
    sketch.add_line(cline1)
    cline2 = lines.Line2D([1.5*sNR + wN/2.0, 6.0*sNR + wN + wR],
                          [-hN/2.0, -hN/2.0], color='gray', linestyle='-')
    sketch.add_line(cline2)
    cline3 = lines.Line2D([2.5*sNR + wN + wR/2.0, 9.0*sNR + wN + wR],
                          [hR/2.0-0.02*hR, hR/2.0-0.02*hR], color='gray', linestyle='-')
    sketch.add_line(cline3)
    cline4 = lines.Line2D([2.5*sNR + wN + wR/2.0, 9.0*sNR + wN + wR],
                          [-hR/2.0, -hR/2.0], color='gray', linestyle='-')
    sketch.add_line(cline4)
    origin_x = [5.0*sNR + wN + wR, 5.0*sNR + wN +
                wR, 8.0*sNR + wN + wR, 8.0*sNR + wN + wR]
    origin_y = [0.0, 0.0, 0.0, 0.0]
    direction_x = [0.0, 0.0, 0.0, 0.0]
    direction_y = [hN/2.0, -hN/2.0, hR/2.0-0.02*hR, -hR/2.0]
    sketch.quiver(origin_x, origin_y, direction_x, direction_y, color=[
                  'gray', 'gray', 'gray', 'gray'], angles='xy', scale_units='xy', scale=1)
    sketch.text(5.5*sNR + wN + wR, 0.0, "h$_N$", color='gray')
    sketch.text(8.5*sNR + wN + wR, 0.0, "h$_R$", color='gray')

    ypos2 = hR/2.0 + (wR/2.0)*np.tan(np.radians(flr/2.0))
    xpos1 = sNR
    ypos1 = hN/2.0 - (wN/2.0)*np.tan(np.radians(flr/2.0))
    cline5 = lines.Line2D(
        [xpos1, xpos1], [-ypos1-0.5*sNR, -ypos2-4.0*sNR], color='gray', linestyle='-')
    sketch.add_line(cline5)
    xpos1 = sNR + wN
    ypos1 = hN/2.0 + (wN/2.0)*np.tan(np.radians(flr/2.0))
    cline6 = lines.Line2D(
        [xpos1, xpos1], [-ypos1-0.5*sNR, -ypos2-4.0*sNR], color='gray', linestyle='-')
    sketch.add_line(cline6)
    xpos1 = 2.0*sNR + wN
    ypos1 = hR/2.0 - (wR/2.0)*np.tan(np.radians(flr/2.0))
    cline7 = lines.Line2D(
        [xpos1, xpos1], [-ypos1-0.5*sNR, -ypos2-4.0*sNR], color='gray', linestyle='-')
    sketch.add_line(cline7)
    xpos1 = 2.0*sNR + wN + wR
    ypos1 = hR/2.0 + (wR/2.0)*np.tan(np.radians(flr/2.0))
    cline8 = lines.Line2D(
        [xpos1, xpos1], [-ypos1-0.5*sNR, -ypos2-4.0*sNR], color='gray', linestyle='-')
    sketch.add_line(cline8)
    cline9 = lines.Line2D([sNR + wN, 2.0*sNR + wN], [-ypos1 -
                                                     3.5*sNR, -ypos2-3.5*sNR], color='gray', linestyle='-')
    sketch.add_line(cline9)
    origin_x = [sNR, sNR + wN, 2.0*sNR + wN, 2.0*sNR + wN + wR]
    origin_y = [-ypos2-3.5*sNR, -ypos2-3.5*sNR, -ypos2-3.5*sNR, -ypos2-3.5*sNR]
    direction_x = [wN, -wN, wR, -wR]
    direction_y = [0.0, 0.0, 0.0, 0.0]
    sketch.quiver(origin_x, origin_y, direction_x, direction_y, color=[
                  'gray', 'gray', 'gray', 'gray'], angles='xy', scale_units='xy', scale=1)
    sketch.text(sNR + wN/2.0 - 0.005, -ypos2-5.0*sNR, "w$_N$", color='gray')
    sketch.text(1.5*sNR + wN - 0.005, -ypos2 -
                5.0*sNR, "s$_N$$_R$", color='gray')
    sketch.text(2.0*sNR + wN + wR/2.0 - 0.005, -
                ypos2-5.0*sNR, "w$_R$", color='gray')

    plt.axis('equal')
    plt.axis('off')
    fig.tight_layout()
    return fig


def stress_analysis():
    if ((not text19.get()) | (not text20.get()) | (not text21.get()) | (not text22.get()) | (not text23.get()) | (not text24.get()) | (not text25.get())):
        labbpr121.grid(row=14, column=0, columnspan=5, sticky='w'+'e')
    else:
        labbpr121.grid_forget()
        global stress_frame, stress_frame_flag
        if (stress_frame_flag == False):
            stress_frame = tk.LabelFrame(
                column3_frame, text=" STRESS ANALYSIS ", font=framelab_font)
            stress_frame.grid(row=1, column=0, padx=10, pady=10)

        nozzle.tc = 0.20
        rotor.tc = 0.20
        rotor.kh = 0.02
        nozzle.tes = 0.02
        rotor.tes = 0.02
        nozzle.teo = 0.04
        rotor.teo = 0.04

        global sigma_ct, sigma_gb, data_B, data_n
        rho_bl = 8000  # (kg/m³)
        sigma_ct = (4.0*np.pi*N**2.0*rho_bl*(G2.A+G3.A)/6.0)/10**6.0
        camber = P2h.beta + P3h.beta
        B_curve = interp1d(data_B[0], data_B[1], kind='linear')
        B = B_curve(camber)
        n_curve = interp1d(data_n[0], data_n[1], kind='linear')
        n = n_curve(camber)
        z = (10.0*rotor.tc)**n/B
        sigma_gb = (m*(P2m.Vu + P3m.Vu)*rotor.h /
                    (2.0*rotor.n*z*rotor.c**3.0))/10**6.0

        figure10 = graph_stress()
        canvas10 = FigureCanvasTkAgg(figure10, master=stress_frame)
        canvas10.draw()
        canvas10.get_tk_widget().grid(row=0, column=0)
        plt.close(figure10)


def graph_stress():
    fig, axs = plt.subplots(figsize=(6, 5))
    axs.text(20, 92.5, "Turbine inlet temperature T$_0$$_1$", color='black')
    for i in range(0, 3):
        axs.plot(data_stress[2*i], data_stress[2*i+1], 'k-')
    axs.text(np.amax(data_stress[0])-30, np.amin(data_stress[1]
                                                 )+4, "1100 K", rotation=-27, color='black')
    axs.text(np.amax(data_stress[2])-30, np.amin(data_stress[3]
                                                 )+5, "1200 K", rotation=-35, color='black')
    axs.text(np.amax(data_stress[4])-20, np.amin(data_stress[5]
                                                 )+6, "1300 K", rotation=-55, color='black')
    axs.set(xlabel=r"Permissible $\sigma$$_c$$_t$ [MN/m$^2$]",
            ylabel=r"Permissible $\sigma$$_g$$_b$ [MN/m$^2$]")
    axs.set_title((r"$\sigma$$_c$$_t$ = %.3f MN/m$^2$ ; $\sigma$$_g$$_b$ = %.3f MN/m$^2$" + "\n") %
                  (sigma_ct, sigma_gb), fontsize=12)
    axs.spines['top'].set_visible(False)
    axs.spines['right'].set_visible(False)
    plt.xlim(0, 250)
    plt.ylim(0, 100)

    axs.scatter(sigma_ct, sigma_gb, color='red', marker='o', alpha=0.5)
    vertical_line = lines.Line2D([sigma_ct, sigma_ct], [
                                 0, sigma_gb], color='red', linestyle='--', alpha=0.5)
    axs.add_line(vertical_line)
    horizontal_line = lines.Line2D(
        [0, sigma_ct], [sigma_gb, sigma_gb], color='red', linestyle='--', alpha=0.5)
    axs.add_line(horizontal_line)
    axs.text(60, 10, "0.2% creep strain\n" + "in 10,000 h for\n" +
             "uncooled blade", horizontalalignment='center', color='black')
    return fig


def design_point_losses():
    global design_point_AM, design_point_DC, design_point_KO
    design_point_AM = losses()
    design_point_DC = losses()
    design_point_KO = losses()

    design_point_AM, design_point_DC, design_point_KO = loss_models(
        P1m, P2m, P3m, P1h.Mach, P2h.Machrel)

    global dloss_frame, noteb_dloss, losses_frame_flag
    global AM_frame, DC_frame, KO_frame
    if (losses_frame_flag == False):
        dloss_frame = tk.LabelFrame(
            losses_window, text=" DESIGN POINT LOSSES ", font=framelab_font)
        dloss_frame.grid(row=0, column=0, padx=10, pady=10, sticky='n')
        noteb_dloss = ttk.Notebook(dloss_frame)
        noteb_dloss.pack(side='top')
        AM_frame = tk.Frame(noteb_dloss)
        noteb_dloss.add(AM_frame, text='Ainley-Mathieson')
        AM_btn = tk.Button(AM_frame, text="Choose Model to Iterate",
                           command=AM_iteration, font=("times new roman", 12))
        AM_btn.grid(row=1, column=0, sticky='w'+'e', padx=5, pady=5)
        DC_frame = tk.Frame(noteb_dloss)
        noteb_dloss.add(DC_frame, text='Dunham-Came')
        DC_btn = tk.Button(DC_frame, text="Choose Model to Iterate",
                           command=DC_iteration, font=("times new roman", 12))
        DC_btn.grid(row=1, column=0, sticky='w'+'e', padx=5, pady=5)
        KO_frame = tk.Frame(noteb_dloss)
        noteb_dloss.add(KO_frame, text='Kacker-Okapuu')
        KO_btn = tk.Button(KO_frame, text="Choose Model to Iterate",
                           command=KO_iteration, font=("times new roman", 12))
        KO_btn.grid(row=1, column=0, sticky='w'+'e', padx=5, pady=5)

    global iteration_flag
    if (iteration_flag == int(0)):
        figure11 = graph_design_point_losses(design_point_AM, 1)
        canvas11 = FigureCanvasTkAgg(figure11, master=AM_frame)
        canvas11.draw()
        canvas11.get_tk_widget().grid(row=0, column=0)
        plt.close(figure11)

        figure12 = graph_design_point_losses(design_point_DC, 2)
        canvas12 = FigureCanvasTkAgg(figure12, master=DC_frame)
        canvas12.draw()
        canvas12.get_tk_widget().grid(row=0, column=0)
        plt.close(figure12)

        figure13 = graph_design_point_losses(design_point_KO, 3)
        canvas13 = FigureCanvasTkAgg(figure13, master=KO_frame)
        canvas13.draw()
        canvas13.get_tk_widget().grid(row=0, column=0)
        plt.close(figure13)
    elif (iteration_flag == int(1)):
        AM_iteration()
    elif (iteration_flag == int(2)):
        DC_iteration()
    elif (iteration_flag == int(3)):
        KO_iteration()

    losses_frame_flag = True


def AM_iteration():
    global iteration_flag, design_point_AM, lbdN
    iteration_flag = int(1)
    if (np.abs(lbdN - design_point_AM.lbdn) > 0.0001):
        lbdN = design_point_AM.lbdn
    else:
        iteration_flag = int(0)
    read_input()


def DC_iteration():
    global iteration_flag, design_point_DC, lbdN
    iteration_flag = int(2)
    if (np.abs(lbdN - design_point_DC.lbdn) > 0.0001):
        lbdN = design_point_DC.lbdn
    else:
        iteration_flag = int(0)
    read_input()


def KO_iteration():
    global iteration_flag, design_point_KO, lbdN
    iteration_flag = int(3)
    if (np.abs(lbdN - design_point_KO.lbdn) > 0.0001):
        lbdN = design_point_KO.lbdn
    else:
        iteration_flag = int(0)
    read_input()


def graph_design_point_losses(loss, option):

    # option: (1 - AM), (2 - DC), (3 - KO)

    fig, axs = plt.subplots(2, 2, figsize=(10, 12))
    if (option == 1):
        axs[0, 0].set_title("Ainley-Mathieson Model", fontsize=16)
    elif (option == 2):
        axs[0, 0].set_title("Dunham-Came Model", fontsize=16)
    else:
        axs[0, 0].set_title("Kacker-Okapuu Model", fontsize=16)

    total = loss.Ytn + loss.Ytr
    fracs = [loss.Ytn/total, loss.Ytr/total]
    color = ['blue', 'red']
    axs[0, 0].pie(fracs, colors=color, startangle=90,
                  autopct=None, wedgeprops=dict(width=0.4))
    axs[0, 0].legend(["Nozzle = %.2f%%" % (
        fracs[0]*100), " Rotor  = %.2f%%" % (fracs[1]*100)], loc='center', fontsize=12)
    axs[0, 0].axis('equal')

    table_data = [[r"Y$_N$ =", "%.4f" % loss.Ytn],
                  [r"$\lambda_N$ =", "%.4f" % loss.lbdn],
                  [r"Y$_R$ =", "%.4f" % loss.Ytr],
                  [r"$\lambda_R$ =", "%.4f" % loss.lbdr],
                  [r"$\eta_s$ =", "%.4f" % loss.eff]]
    table1 = axs[0, 1].table(cellText=table_data, colLabels=None, edges='open',
                             cellLoc='center', rowLoc='center', colLoc='center', loc='center')
    table1.auto_set_font_size(False)
    table1.set_fontsize(16)
    table1.scale(0.35, 2.0)
    axs[0, 1].axis('off')

    if (option == 3):
        total = loss.Ypn + loss.Ysn + loss.Yten
        fracs = [loss.Ypn/total, loss.Ysn/total, loss.Yten/total]
        color = ['gray', 'green', 'yellow']
    else:
        total = loss.Ypn + loss.Ysn
        fracs = [loss.Ypn/total, loss.Ysn/total]
        color = ['gray', 'green']
    axs[1, 0].pie(fracs, colors=color, startangle=90,
                  autopct=None, wedgeprops=dict(width=0.4))
    axs[1, 0].set_title("Nozzle", fontsize=16)
    if (option == 3):
        axs[1, 0].legend(['Y$_p$ = %.2f%%' % (fracs[0]*100), 'Y$_s$ = %.2f%%' %
                          (fracs[1]*100), 'Y$_t$$_e$ = %.2f%%' % (fracs[2]*100)], loc='center', fontsize=12)
    else:
        axs[1, 0].legend(['Y$_p$ = %.2f%%' % (
            fracs[0]*100), 'Y$_s$ = %.2f%%' % (fracs[1]*100)], loc='center', fontsize=12)
    axs[1, 0].axis('equal')

    if (option == 3):
        total = loss.Ypr + loss.Ysr + loss.Ykr + loss.Yter
        fracs = [loss.Ypr/total, loss.Ysr/total,
                 loss.Ykr/total, loss.Yter/total]
        color = ['gray', 'green', 'orange', 'yellow']
    else:
        total = loss.Ypr + loss.Ysr + loss.Ykr
        fracs = [loss.Ypr/total, loss.Ysr/total, loss.Ykr/total]
        color = ['gray', 'green', 'orange']
    axs[1, 1].pie(fracs, colors=color, startangle=90,
                  autopct=None, wedgeprops=dict(width=0.4))
    axs[1, 1].set_title("Rotor", fontsize=16)
    if (option == 3):
        axs[1, 1].legend(['Y$_p$ = %.2f%%' % (fracs[0]*100), 'Y$_s$ = %.2f%%' % (fracs[1]*100), 'Y$_k$ = %.2f%%' %
                          (fracs[2]*100), 'Y$_t$$_e$ = %.2f%%' % (fracs[3]*100)], loc='center', fontsize=12)
    else:
        axs[1, 1].legend(['Y$_p$ = %.2f%%' % (fracs[0]*100), 'Y$_s$ = %2.2f%%' %
                          (fracs[1]*100), 'Y$_k$ = %.2f%%' % (fracs[2]*100)], loc='center', fontsize=12)
    axs[1, 1].axis('equal')

    fig.tight_layout(pad=3.0)
    return fig


def tournier_genk(sigma_, beta_out):

    if (beta_out <= 63.2):
        A = -5.58*10**(-5.0)*beta_out**2.0 + 1.03*10**(-2.0)*beta_out - 0.275
        B = 1.553*10**(-5.0)*beta_out**2.0 - 2.32 * \
            10**(-3.0)*beta_out + 8.02*10**(-2.0)
        C = -8.54*10**(-3.0)*beta_out + 0.238
        D = 4.83*10**(-5.0)*beta_out**2.0 + 2.83 * \
            10**(-4.0)*beta_out - 2.93*10**(-2.0)
    else:
        A = 2.44*10**(-4.0)*beta_out**2.0 - 4.33*10**(-2.0)*beta_out + 1.92
        B = -4.02*10**(-5.0)*beta_out**2.0 + 6.94*10**(-3.0)*beta_out - 0.282
        C = -2.23*10**(-4.0)*beta_out**2.0 + 4.96*10**(-2.0)*beta_out - 2.548
        D = 5.39*10**(-5.0)*beta_out**2.0 - 1.57*10**(-2.0)*beta_out + 0.958
    Yp0_axial = A + B/sigma_ + sigma_*(C + D*sigma_)

    if (beta_out <= 60.0):
        sigma_min = -8.63*10**(-6.0)*beta_out**3.0 + 9.68*10**(-4.0) * \
            beta_out**2.0 - 3.76*10**(-2.0)*beta_out + 1.272
        if (sigma_ >= sigma_min):
            n = 1.524*10**(-4.0)*beta_out**2.0 - 0.031*beta_out + 2.992
            A = -2.91*10**(-3.0)*beta_out + 0.30260
        else:
            n = 3.271*10**(-3.0)*beta_out**2.0 - 0.3010*beta_out + 9.023
            A = 2.701*10**(-3.0)*beta_out**2.0 - 0.2456*beta_out + 5.909
    else:
        sigma_min = -5.14*10**(-4.0)*beta_out**2.0 + \
            5.48*10**(-2.0)*beta_out - 0.798
        if (sigma_ >= sigma_min):
            n = 1.524*10**(-4.0)*beta_out**2.0 - 0.031*beta_out + 2.992
            A = 5.407*10**(-3.0)*beta_out - 0.19642
        else:
            n = 1.174*10**(-2.0)*beta_out**2.0 - 1.5731*beta_out + 54.85
            A = 9.240*10**(-3.0)*beta_out**2.0 - 1.2067*beta_out + 40.04
    Yp0_impulse = 0.280*(1.0 - sigma_min) + A*np.abs(sigma_ - sigma_min)**n

    return Yp0_axial, Yp0_impulse


def graph_tournier_genk():
    fig1, axs1 = plt.subplots()
    fig2, axs2 = plt.subplots()
    for beta_out in ([40, 50, 60, 65, 70, 75, 80]):
        Yp_0_vector, Yp_1_vector = [], []
        sc_vector = np.arange(0.30, 1.21, 0.05)
        for sigma_ in sc_vector:
            Yp_0, Yp_1 = tournier_genk(sigma_, beta_out)
            Yp_0_vector.append(Yp_0)
            Yp_1_vector.append(Yp_1)
        axs1.plot(sc_vector, Yp_0_vector, 'k-')
        axs1.text(1.225, Yp_0_vector[len(Yp_0_vector)-1],
                  "%d°" % beta_out, color='black')
        if (beta_out <= 70):
            axs2.plot(sc_vector, Yp_1_vector, 'k-')
            axs2.text(1.225, Yp_1_vector[len(
                Yp_1_vector)-1], "%d°" % beta_out, color='black')
    axs1.text(0.80, 0.09, "Relative outlet gas\n" +
              r"angle ($\alpha_2$ or $\beta_3$)", color='black')
    axs1.set(xlabel="pitch/chord (s/c)",
             ylabel=r"Profile loss coefficient ($Y_p$)")
    axs1.set_xlim([0.3, 1.3])
    axs1.set_ylim([0.0, 0.10])
    axs1.spines['top'].set_visible(False)
    axs1.spines['right'].set_visible(False)
    axs2.text(0.80, 0.22, "Relative outlet gas\n" +
              r"angle ($\alpha_2$ or $\beta_3$)", color='black')
    axs2.set(xlabel="pitch/chord (s/c)",
             ylabel=r"Profile loss coefficient ($Y_p$)")
    axs2.set_xlim([0.3, 1.3])
    axs2.set_ylim([0.06, 0.24])
    axs2.spines['top'].set_visible(False)
    axs2.spines['right'].set_visible(False)
    return fig1, fig2


def ainley_mathieson(blade, Gin, Gout, flow_in_angle, flow_out_angle):

    #####                        PROFILE LOSSES                         #####

    Yp0_axial, Yp0_impulse = tournier_genk(blade.sc, flow_out_angle)
    Yp0 = (Yp0_axial + ((blade.in_angle/flow_out_angle)**2.0)*(Yp0_impulse -
                                                               Yp0_axial))*((blade.tc/0.2)**(blade.in_angle/flow_out_angle))

    ## CORREÇÕES das Fig. 7b, 7c, 7a e 6 ##
    ##                Xi                 ##
    Xi = 1.0
    Yp = Yp0*Xi

    #####               SECONDARY AND TIP CLEARANCE LOSSES               #####
    file_lambd = 'csv_data/lambd.csv'
    data_lambd = pd.read_csv(file_lambd, header=None, skiprows=[0, 1])
    lambd_curve = interp1d(data_lambd[0], data_lambd[1], kind='linear')
    lambd_parameter = ((Gout.A*np.cos(np.radians(flow_out_angle))) /
                       (Gin.A*np.cos(np.radians(flow_in_angle))))**2.0/(1.0 + 1/blade.rr)
    lambd = lambd_curve(lambd_parameter)

    if (int(shroud_flag.get()) == 1):
        B = 0.50
    else:
        B = 0.25

    flow_mean_angle = np.arctan(
        (np.tan(np.radians(flow_out_angle)) - np.tan(np.radians(flow_in_angle)))/2.0)
    Z = ((2.0*((np.tan(np.radians(flow_in_angle)) + np.tan(np.radians(flow_out_angle)))*np.cos(flow_mean_angle)))
         ** 2.0) * ((np.cos(np.radians(flow_out_angle)))**2.0/(np.cos(flow_mean_angle))**3.0)

    Ys = lambd*Z
    Yk = B*blade.kh*Z

    #####                          TOTAL LOSSES                          #####
    #####                 TRAILING EDGE THICKNESS EFFECT                 #####

    data_Xte = pd.read_csv('Xte.csv', header=None, skiprows=[0, 1])
    Xte_curve = interp1d(data_Xte[0], data_Xte[1], kind='linear')
    if (blade.tes < 0.12):
        Xte = Xte_curve(blade.tes)
    else:
        Xte = 1.7

    Yt = (Yp + Ys + Yk)*Xte

    return Yp0, Z, Xte, Yp, Ys, Yk, Yt


def dunham_came(blade, Re, Ma, Yp_AM, Z, Xte, flow_out_angle):

    #####                        PROFILE LOSSES                         #####

    if (Ma > 1.0):
        Yp = Yp_AM*(1 + 60*(Ma - 1)**2.0)
    else:
        Yp = Yp_AM

    Yp = Yp*(Re/200000)**(-0.2)

    #####               SECONDARY AND TIP CLEARANCE LOSSES               #####

    if (int(shroud_flag.get()) == 1):
        B = 0.47
    else:
        B = 0.37

    Ys0 = 0.0334*(1/blade.hc)*(np.cos(np.radians(flow_out_angle)) /
                               np.cos(np.radians(blade.in_angle)))*Z

    Ys = Ys0*(Re/200000)**(-0.2)
    Yk = B*(1/blade.hc)*((blade.kh*blade.hc)**(0.78))*Z

    #####                          TOTAL LOSSES                          #####
    #####                 TRAILING EDGE THICKNESS EFFECT                 #####

    Yt = (Yp + Ys + Yk)*Xte

    return Ys0, Yp, Ys, Yk, Yt


def kacker_okapuu(blade, Re, Ma_in, Ma_out, Mah_in, p_in, p_out, Yp0, Ys0, Yk_DC, flow_out_angle):

    #####                        PROFILE LOSSES                         #####

    if (Re <= 200000):
        XRe = (Re/200000)**(-0.4)
    elif (Re > 200000 and Re < 2000000):
        XRe = 1
    else:
        XRe = (Re/2000000)**(-0.2)

    if (Ma_out > 0.2):
        K1 = 1.0 - 1.25*(Ma_out - 0.2)
    else:
        K1 = 1.0

    K2 = (Ma_in/Ma_out)**2.0
    Kp = 1.0 - K2*(1.0 - K1)

    Yshock = 0.75*((Mah_in - 0.4)**(1.75))*(1/blade.rr)*(p_in/p_out)*(1 - (1 + ((gamma - 1)/2.0)*Ma_in**2.0)
                                                                      ** (gamma/(gamma - 1)))/(1 - (1 + ((gamma - 1)/2.0)*Ma_out**2.0)**(gamma/(gamma - 1)))
    Yp = 0.914*XRe*((2.0/3.0)*Kp*Yp0 + Yshock)

    #####               SECONDARY AND TIP CLEARANCE LOSSES               #####

    K3 = (1/blade.hw)**2.0
    Ks = 1.0 - K3*(1.0 - Kp)

    if (blade.hc <= 2.0):
        Xar = 1.0 - 0.25*np.sqrt(2.0 - blade.hc)
    else:
        Xar = 1.0

    Ys = 1.2*Ys0*Xar*Ks
    Yk = Yk_DC

    #####                 TRAILING EDGE THICKNESS EFFECT                 #####

    TEaxial = 0.59563*blade.teo**2.0 + 0.12264*blade.teo - 2.2796*10**(-3.0)
    TEimpulse = 0.31066*blade.teo**2.0 + 0.065617*blade.teo - 1.4318*10**(-3.0)
    TE_coefficient = TEaxial + (TEimpulse - TEaxial) * \
        (blade.in_angle/flow_out_angle)**2.0

    Yte = ((1 - ((gamma - 1)/2.0)*(Ma_out**2.0)*((1/(1 - TE_coefficient)) - 1))**((-gamma) /
                                                                                  (gamma - 1)) - 1)/(1 - (1 + ((gamma - 1)/2.0)*(Ma_out**2.0))**((-gamma)/(gamma - 1)))

    #####                          TOTAL LOSSES                          #####

    Yt = Yp + Ys + Yk + Yte

    return Yp, Ys, Yk, Yte, Yt


def loss_models(P1, P2, P3, Mach_1h, Machrel_2h):

    Ren = P2.rho*P2.V*nozzle.c/P2.mu
    Rer = P3.rho*P3.W*rotor.c/P3.mu

    losses_AM = losses()
    losses_DC = losses()
    losses_KO = losses()

    #####                        AINLEY-MATHIESON                        #####

    Yp0n, Zn, Xten, Ypn_AM, Ysn_AM, Ykn_AM, Ytn_AM = ainley_mathieson(
        nozzle, G1, G2, P1.alpha, P2.alpha)
    Yp0r, Zr, Xter, Ypr_AM, Ysr_AM, Ykr_AM, Ytr_AM = ainley_mathieson(
        rotor, G2, G3, P2.beta, P3.beta)

    lbdn_AM = Ytn_AM*P2.Ts/P2.T0
    lbdr_AM = Ytr_AM*P3.Ts/(P3.T + (P3.W**2.0)/(2.0*cp))

    eff_AM = 1/(1 + ((lbdr_AM*((P3.W**2.0)/(2.0*cp)) + (P3.T/P2.T)
                      * lbdn_AM*((P2.V**2.0)/(2.0*cp)))/(P1.T0 - P3.T0)))
    eff_AM = 1 - (1 - eff_AM)*(((Ren + Rer)/2.0)/200000)**(-0.2)

    losses_AM.Ypn, losses_AM.Ysn, losses_AM.Ykn, losses_AM.Ytn = Ypn_AM, Ysn_AM, Ykn_AM, Ytn_AM
    losses_AM.Ypr, losses_AM.Ysr, losses_AM.Ykr, losses_AM.Ytr = Ypr_AM, Ysr_AM, Ykr_AM, Ytr_AM
    losses_AM.lbdn, losses_AM.lbdr, losses_AM.eff = lbdn_AM, lbdr_AM, eff_AM

    #####                          DUNHAM-CAME                           #####

    Ys0n, Ypn_DC, Ysn_DC, Ykn_DC, Ytn_DC = dunham_came(
        nozzle, Ren, P2.Mach, Ypn_AM, Zn, Xten, P2.alpha)
    Ys0r, Ypr_DC, Ysr_DC, Ykr_DC, Ytr_DC = dunham_came(
        rotor, Rer, P3.Machrel, Ypr_AM, Zr, Xter, P3.beta)

    lbdn_DC = Ytn_DC*P2.Ts/P2.T0
    lbdr_DC = Ytr_DC*P3.Ts/(P3.T + (P3.W**2.0)/(2.0*cp))

    eff_DC = 1/(1 + ((lbdr_DC*((P3.W**2.0)/(2.0*cp)) + (P3.T/P2.T)
                      * lbdn_DC*((P2.V**2.0)/(2.0*cp)))/(P1.T0 - P3.T0)))

    losses_DC.Ypn, losses_DC.Ysn, losses_DC.Ykn, losses_DC.Ytn = Ypn_DC, Ysn_DC, Ykn_DC, Ytn_DC
    losses_DC.Ypr, losses_DC.Ysr, losses_DC.Ykr, losses_DC.Ytr = Ypr_DC, Ysr_DC, Ykr_DC, Ytr_DC
    losses_DC.lbdn, losses_DC.lbdr, losses_DC.eff = lbdn_DC, lbdr_DC, eff_DC

    #####                         KACKER-OKAPUU                          #####

    # Mach_hub = P1.Mach*(1.0 + 1.8*((1.0 - 1/nozzle.rr)**2.2)) # If nozzle, constant = 1.8
    Ypn_KO, Ysn_KO, Ykn_KO, Yten_KO, Ytn_KO = kacker_okapuu(
        nozzle, Ren, P1.Mach, P2.Mach, Mach_1h, P1.p, P2.p, Yp0n, Ys0n, Ykn_DC, P2.alpha)

    # Mach_hub = P2.Machrel*(1.0 + 5.2*((1.0 - 1/rotor.rr)**2.2)) # If rotor, constante = 5.2
    Ypr_KO, Ysr_KO, Ykr_KO, Yter_KO, Ytr_KO = kacker_okapuu(
        rotor, Rer, P2.Machrel, P3.Machrel, Machrel_2h, P2.p, P3.p, Yp0r, Ys0r, Ykr_DC, P3.beta)

    lbdn_KO = Ytn_KO*P2.Ts/P2.T0
    lbdr_KO = Ytr_KO*P3.Ts/(P3.T + (P3.W**2.0)/(2.0*cp))

    eff_KO = 1/(1 + ((lbdr_KO*((P3.W**2.0)/(2.0*cp)) + (P3.T/P2.T)
                      * lbdn_KO*((P2.V**2.0)/(2.0*cp)))/(P1.T0 - P3.T0)))

    losses_KO.Ypn, losses_KO.Ysn, losses_KO.Ykn, losses_KO.Yten, losses_KO.Ytn = Ypn_KO, Ysn_KO, Ykn_KO, Yten_KO, Ytn_KO
    losses_KO.Ypr, losses_KO.Ysr, losses_KO.Ykr, losses_KO.Yter, losses_KO.Ytr = Ypr_KO, Ysr_KO, Ykr_KO, Yter_KO, Ytr_KO
    losses_KO.lbdn, losses_KO.lbdr, losses_KO.eff = lbdn_KO, lbdr_KO, eff_KO

    return losses_AM, losses_DC, losses_KO


global data_sc, data_stress, data_B, data_n, data_smith
file_optimum_sc = 'csv_data/optimum_sc.csv'
file_perm_stress = 'csv_data/permissible_stress.csv'
file_B_coeff = 'csv_data/gas_bending_B_coefficient.csv'
file_n_coeff = 'csv_data/gas_bending_n_coefficient.csv'
file_visc = 'csv_data/visc_air_temp.csv'
file_smith = 'csv_data/smith_chart.csv'

data_sc = pd.read_csv(file_optimum_sc, header=None, skiprows=[0, 1])
data_stress = pd.read_csv(file_perm_stress, header=None, skiprows=[0, 1])
data_B = pd.read_csv(file_B_coeff, header=None, skiprows=[0, 1])
data_n = pd.read_csv(file_B_coeff, header=None, skiprows=[0, 1])
data_visc_air_temp = pd.read_csv(file_visc, header=None, skiprows=[0, 1])
data_smith = pd.read_csv(file_smith, header=None, skiprows=[0, 1])

global cp, gamma, R
cp = 1148  # (J/kg.K)
gamma = 4.0/3.0
R = 287  # (J/kg.K)

global P1m, P2m, P3m
P1m = point()
P2m = point()
P3m = point()

global G1, G2, G3
G1 = geometry()
G2 = geometry()
G3 = geometry()

hypth_frame_flag = False
read_hypth_flag = False
mean_frame_flag = False
blade_frame_flag = False
root_tip_frame_flag = False
pcnbl_frame_flag = False
sktch_frame_flag = False
blpar_frame_flag = False
stress_frame_flag = False
losses_frame_flag = False
iteration_flag = int(0)

standard_font = ("times new roman", 10)
framelab_font = ("times new roman", 12, "bold", "italic")
