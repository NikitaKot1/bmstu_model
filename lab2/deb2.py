from math import log, exp, pi
from prettytable import PrettyTable
from matplotlib import pyplot

import pandas as pd
import matplotlib.pyplot as plt

R = 0.35
l = 12
L_k = 187 * (10**(-6))
C_k = 268 * (10**(-6))
R_k = 0.25
R_k2 = - 0.35
U_co = 1400
I_o = 3
T_w = 2000
STEP = 1e-6


def I_analyt():
    pass


def f(t, I, U):
    return (U - (R_k2 + R) * I) / L_k


def f2(t, I, U, I_arr, T0_arr, m_arr, T_arr, sigma_arr):
    return (U - (R_k + find_R(I, I_arr, T0_arr, m_arr, T_arr, sigma_arr)) * I) / L_k


def phi(t, I):
    return -(I / C_k)


def runge4_withR(I0, U0, h, I_arr, T0_arr, m_arr, T_arr, sigma_arr, t0=0, t_max=0.01):
    I_n = I0
    U_n = U0
    t_n = t0

    t_res = []
    I_res = []
    U_res = []
    R0 = find_R(I0, I_arr, T0_arr, m_arr, T_arr, sigma_arr)
    T0, m = find_T0_m(I0, I_arr, T0_arr, m_arr)
    R_res = []
    T0_res = []

    while t_n < t_max:
        k1 = h * f2(t_n, I_n, U_n, I_arr, T0_arr, m_arr, T_arr, sigma_arr)
        q1 = h * phi(t_n, I_n)
        k2 = h * f2(t_n + h / 2, I_n + k1 / 2, U_n + q1 / 2, I_arr, T0_arr, m_arr, T_arr, sigma_arr)
        q2 = h * phi(t_n + h / 2, I_n + k1 / 2)
        k3 = h * f2(t_n + h / 2, I_n + k2 / 2, U_n + q2 / 2, I_arr, T0_arr, m_arr, T_arr, sigma_arr)
        q3 = h * phi(t_n + h / 2, I_n + k2 / 2)
        k4 = h * f2(t_n + h, I_n + k3, U_n + q3, I_arr, T0_arr, m_arr, T_arr, sigma_arr)
        q4 = h * phi(t_n + h / 2, I_n + k3)

        t_n = t_n + h
        I_n = I_n + (k1 + 2*k2 + 2*k3 + k4) / 6
        U_n = U_n + (q1 + 2*q2 + 2*q3 + q4) / 6

        R_p = find_R(I_n, I_arr, T0_arr, m_arr, T_arr, sigma_arr)
        T0, m = find_T0_m(I_n, I_arr, T0_arr, m_arr)

        t_res.append(t_n)
        I_res.append(I_n)
        U_res.append(U_n)
        R_res.append(R_p)
        T0_res.append(T0)

    return t_res, I_res, U_res, R_res, T0_res


def runge4(I0, U0, h, t0=0, t_max=0.01):
    I_n = I0
    U_n = U0
    t_n = t0

    t_res = [t0]
    I_res = [I0]
    U_res = [U0]

    while t_n < t_max:
        k1 = h * f(t_n, I_n, U_n)
        q1 = h * phi(t_n, I_n)
        k2 = h * f(t_n + h / 2, I_n + k1 / 2, U_n + q1 / 2)
        q2 = h * phi(t_n + h / 2, I_n + k1 / 2)
        k3 = h * f(t_n + h / 2, I_n + k2 / 2, U_n + q2 / 2)
        q3 = h * phi(t_n + h / 2, I_n + k2 / 2)
        k4 = h * f(t_n + h, I_n + k3, U_n + q3)
        q4 = h * phi(t_n + h / 2, I_n + k3)

        t_n = t_n + h
        I_n = I_n + (k1 + 2*k2 + 2*k3 + k4) / 6
        U_n = U_n + (q1 + 2*q2 + 2*q3 + q4) / 6

        t_res.append(t_n)
        I_res.append(I_n)
        U_res.append(U_n)

    return t_res, I_res, U_res


def runge2_withR(I0, U0, h, I_arr, T0_arr, m_arr, T_arr, sigma_arr, t0=0, t_max=0.01, beta=1/2):
    I_n = I0
    U_n = U0
    t_n = t0

    t_res = [t0]
    I_res = [I0]
    U_res = [U0]
    R0 = find_R(I0, I_arr, T0_arr, m_arr, T_arr, sigma_arr)
    T0, m = find_T0_m(I0, I_arr, T0_arr, m_arr)
    R_res = [R0]
    T0_res = [T0]

    while t_n < t_max:
        k1 = h * f2(t_n, I_n, U_n, I_arr, T0_arr, m_arr, T_arr, sigma_arr)
        q1 = h * phi(t_n, I_n)
        k2 = h * f2(t_n + h / (2*beta), I_n + k1 / (2*beta), U_n + q1 / (2*beta),
                    I_arr, T0_arr, m_arr, T_arr, sigma_arr)
        q2 = h * phi(t_n + h / (2*beta), I_n + k1 / (2*beta))

        t_n = t_n + h
        I_n = I_n + (1 - beta) * k1 + beta * k2
        U_n = U_n + (1 - beta) * q1 + beta * q2

        R_p = find_R(I_n, I_arr, T0_arr, m_arr, T_arr, sigma_arr)
        T0, m = find_T0_m(I_n, I_arr, T0_arr, m_arr)

        t_res.append(t_n)
        I_res.append(I_n)
        U_res.append(U_n)
        R_res.append(R_p)
        T0_res.append(T0)

    return t_res, I_res, U_res, R_res, T0_res


def runge2(I0, U0, h, t0=0, t_max=0.01, beta=1/2):
    I_n = I0
    U_n = U0
    t_n = t0

    t_res = [t0]
    I_res = [I0]
    U_res = [U0]


    while t_n < t_max:
        k1 = h * f(t_n, I_n, U_n)
        q1 = h * phi(t_n, I_n)
        k2 = h * f(t_n + h / (2*beta), I_n + k1 / (2*beta), U_n + q1 / (2*beta))
        q2 = h * phi(t_n + h / (2*beta), I_n + k1 / (2*beta))

        t_n = t_n + h
        I_n = I_n + (1 - beta) * k1 + beta * k2
        U_n = U_n + (1 - beta) * q1 + beta * q2

        t_res.append(t_n)
        I_res.append(I_n)
        U_res.append(U_n)

    return t_res, I_res, U_res


def euler_withR(I0, U0, h, I_arr, T0_arr, m_arr, T_arr, sigma_arr, t0=0, t_max=0.01):
    I_n = I0
    U_n = U0
    t_n = t0

    t_res = [t0]
    I_res = [I0]
    U_res = [U0]
    R0 = find_R(I0, I_arr, T0_arr, m_arr, T_arr, sigma_arr)
    T0, m = find_T0_m(I0, I_arr, T0_arr, m_arr)
    R_res = [R0]
    T0_res = [T0]

    while t_n < t_max:
        k1 = h * f2(t_n, I_n, U_n, I_arr, T0_arr, m_arr, T_arr, sigma_arr)
        q1 = h * phi(t_n, I_n)

        t_n = t_n + h
        I_n = I_n + k1
        U_n = U_n + q1

        R_p = find_R(I_n, I_arr, T0_arr, m_arr, T_arr, sigma_arr)
        T0, m = find_T0_m(I_n, I_arr, T0_arr, m_arr)
        t_res.append(t_n)
        I_res.append(I_n)
        U_res.append(U_n)
        R_res.append(R_p)
        T0_res.append(T0)

    return t_res, I_res, U_res, R_res, T0_res


def euler(I0, U0, h, t0=0, t_max=0.01):
    I_n = I0
    U_n = U0
    t_n = t0

    t_res = [t0]
    I_res = [I0]
    U_res = [U0]

    while t_n < t_max:
        k1 = h * f(t_n, I_n, U_n)
        q1 = h * phi(t_n, I_n)

        t_n = t_n + h
        I_n = I_n + k1
        U_n = U_n + q1

        t_res.append(t_n)
        I_res.append(I_n)
        U_res.append(U_n)

    return t_res, I_res, U_res


# функция T
def T_func(T0, z, m):
    return T0 + (T_w - T0) * z**m

# функция R
def R_func(S):
    return l/(2*pi*R**2 * S)
    #return (2*pi*R**2 * S)



def find_T0_m(I, I_arr, T0_arr, m_arr):
    n = len(I_arr)
    j = 0

    if I < I_arr[0]:
        m = m_arr[0]
        T0 = T0_arr[0]
        return T0, m
    # elif I > I_arr[n - 1]:
    #     m = m_arr[n - 1]
    #     T0 = T0_arr[n - 1]
    #     return T0, m


    while True:
        if I_arr[j] > I or j == n - 2:
            break
        j += 1
    j -= 1
    #while j < n - 1 and I_arr[j] > I:
        #j += 1

    if j < n - 1:
        dx = I_arr[j+1] - I_arr[j]
        di = I - I_arr[j]
        T0 = T0_arr[j] + ((T0_arr[j + 1] - T0_arr[j]) * di / dx)
        m = m_arr[j] + ((m_arr[j + 1] - m_arr[j]) * di / dx)
        #print(j, I, I_arr[j+1], I_arr[j])
    else:
        dx = I_arr[n-1] - I_arr[n-2]
        di = I - I_arr[n - 1]
        T0 = T0_arr[n - 2] + ((T0_arr[n - 1] - T0_arr[n - 2]) * di / dx)
        m = m_arr[n - 1]

    if m < 0:
        print(I_arr[-1])
        #print(m, I, fl)
    return T0, m


def find_sigma(T, T_arr, sigma_arr):
    n = len(T_arr)
    j = 0
    if T < T_arr[0]:
        sigma = sigma_arr[0]
        return sigma

    elif T > T_arr[n - 1]:
        sigma = sigma_arr[n - 1]
        return sigma

    while True:
        if T_arr[j] > T or j == n - 2:
            break
        j += 1
    j -= 1

    #while j < n - 1 and T_arr[j] > T:
        #j += 1
    if j < n - 1:
        dx = T_arr[j+1] - T_arr[j]
        di = T - T_arr[j]
        sigma = sigma_arr[j] + ((sigma_arr[j + 1] - sigma_arr[j]) * di / dx)
        #print(T, T_arr[j+1], T_arr[j])
    else:
        dx = T_arr[n - 1] - T_arr[n - 2]
        di = T - T_arr[n - 1]
        sigma = sigma_arr[n - 2] + ((sigma_arr[n - 1] - sigma_arr[n - 2]) * di / dx)

    return sigma


def find_R(I, I_arr, T0_arr, m_arr, T_arr, sigma_arr):
    T0, m = find_T0_m(I, I_arr, T0_arr, m_arr)
    h = 1 / 100
    zarr = []
    z = 0
    z_max = 1
    while z < z_max + h:
        zarr.append(z)
        z = z + h
    S = integral(zarr, m, T_arr, sigma_arr, T0)
    R = R_func(S)
    return R

# метод трапеций
def integral(zarr, m, T_arr, sigma_arr, T0):
    l = len(zarr)
    s = 0
    t2 = T = T_func(T0, zarr[0], m)
    s2 = find_sigma(t2, T_arr, sigma_arr)
    for i in range(l - 1):
        t1 = t2
        t2 = T_func(T0, zarr[i+1], m)
        s1 = s2
        s2 = find_sigma(t2, T_arr, sigma_arr)
        s += ((s2 + s1) / 2) * (zarr[i+1] - zarr[i]) * ((zarr[i+1] + zarr[i]) / 2)
    return s


# линейная с вырвнивающими коэффициентами lnx, lny
def interpolation(x_arr, y_arr, h):
    x_arr = [log(x) for x in x_arr]
    y_arr = [log(y) for y in y_arr]
    res_x = []
    res_y = []
    for i in range(len(x_arr) - 1):
        dx = x_arr[i + 1] - x_arr[i]
        dy = y_arr[i + 1] - y_arr[i]
        k = dy / dx
        x = x_arr[i]
        y = y_arr[i]
        #k *= y / x
        while (x + h) < x_arr[i + 1]:
            res_x.append(x)
            res_y.append(y)
            x += h
            y += k * h
        if i == len(x_arr) - 1:
            res_x.append(x)
            res_y.append(y)
    res_x = [exp(x) for x in res_x]
    res_y = [exp(y) for y in res_y]
    return res_x, res_y


if __name__ == "__main__":

    STEP1 = 0.0001
    MAX1 = 0.01
    STEP2 = 1e-6
    MAX2 = 780*1e-6
    tb = PrettyTable()
    #ch = 17
    I = [0.5, 1, 5, 10, 50, 200, 400, 800, 1200]
    T0 = [6730, 6790, 7150, 7270, 8010, 9185, 10010, 11140, 12010]
    m = [0.50, 0.55, 1.7, 3, 11, 32, 40, 41, 39]
    T = [4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000]
    sigma = [0.031, 0.27, 2.05, 6.06, 12.0, 19.9, 29.6, 41.1, 54.1, 67.7, 81.5]

    print("Task №1")
    print("     1) Calculate by Euler")
    print("     2) Calculate by Runge2")
    print("     3) Calculate by Runge4")
    print("Task №2")
    print("     4) Calculate 'I(t) RK4 Rp+Rk=0'")
    print("Task №3")
    print("     5) Calculate 'I(t) with Rk=200'")

    ch = int(input("Choose numbers (1 - 5): "))

    if ch in range(1, 4):

        fig, axes = plt.subplots(nrows=2, ncols=2, sharex=False, sharey=False, figsize=(10, 10))
        #ax1, ax2, ax3, ax4, ax5, ax6 = axes.flatten()
        ax1, ax2, ax3, ax4 = axes.flatten()
        colors = ['tab:orange', 'tab:red', 'tab:blue', 'tab:green', 'tab:purple', 'tab:yellow']

    iT, isigma = interpolation(T, sigma, 0.0001)
    if ch == 1:
        x_res, y_res, z_res, r_res, t_res = euler_withR(I_o, U_co, STEP2, I, T0, m, T, sigma, 0, MAX2)

        ax1.set_title("I(t) RK1")
        ax1.plot(x_res, y_res, colors[0])
        ax2.set_title("U(t) RK1")
        ax2.plot(x_res, z_res, colors[1])
        x_res2 = [x_res[i] for i in range(4, len(x_res))]
        r_res2 = [r_res[i] for i in range(4, len(x_res))]
        ax3.set_title("Rp(t) RK1")
        ax3.plot(x_res2, r_res2, colors[2])
        ir_res = [y_res[i] * r_res[i] for i in range(len(y_res))]
        # ax4.set_title("I(t) * Rp(t) RK1")
        # ax4.plot(x_res, ir_res, colors[3])
        ax4.set_title("T0(t) RK1")
        ax4.plot(x_res, t_res, colors[4])

    elif ch == 2:
        x_res, y_res, z_res, r_res, t_res = runge2_withR(I_o, U_co, STEP2, I, T0, m, T, sigma, 0, MAX2)

        ax1.set_title("I(t) RK2")
        ax1.plot(x_res, y_res, colors[0])
        ax2.set_title("U(t) RK2")
        ax2.plot(x_res, z_res, colors[1])
        x_res2 = [x_res[i] for i in range(4, len(x_res))]
        r_res2 = [r_res[i] for i in range(4, len(x_res))]
        ax3.set_title("Rp(t) RK2")
        ax3.plot(x_res2, r_res2, colors[2])
        ir_res = [y_res[i] * r_res[i] for i in range(len(y_res))]
        # ax4.set_title("I(t) * Rp(t) RK2")
        # ax4.plot(x_res, ir_res, colors[3])
        ax4.set_title("T0(t) RK2")
        ax4.plot(x_res, t_res, colors[4])

    elif ch == 3:
        x_res, y_res, z_res, r_res, t_res = runge4_withR(I_o, U_co, STEP2, I, T0, m, T, sigma, 0, MAX2)
        
        ax1.set_title("I(t) RK4")
        ax1.plot(x_res, y_res, colors[0])
        ax2.set_title("U(t) RK4")
        ax2.plot(x_res, z_res, colors[1])
        x_res2 = [x_res[i] for i in range(3, len(x_res))]
        r_res2 = [r_res[i] for i in range(3, len(x_res))]
        ax3.set_title("Rp(t) RK4")
        ax3.plot(x_res2, r_res2, colors[2])
        ir_res = [y_res[i] * r_res[i] for i in range(len(y_res))]
        # ax4.set_title("I(t) * Rp(t) RK4")
        # ax4.plot(x_res, ir_res, colors[3])
        ax4.set_title("T0(t) RK4")
        ax4.plot(x_res, t_res, colors[4])

    elif ch == 4:
        x_res, y_res, z_res = runge4(I_o, U_co, STEP1 / 100, 0, 1)
        pyplot.plot(x_res, y_res)
        pyplot.title("I(t) RK4 Rp+Rk=0")

    elif ch == 5:
        R_k = 200
        x_res, y_res, z_res, r_res, t_res = runge4_withR(I_o, U_co, 1e-7, I, T0, m, T, sigma, 0, 20*1e-6)
        pyplot.plot(x_res, y_res)
        pyplot.title("I(t) with Rk=200 RK4")

    plt.show()