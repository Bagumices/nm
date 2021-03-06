import numpy as np


def f1(x, y1, y2):
    return np.sin(alpha * y1 * y1) + x + y2


def f2(x, y1, y2):
    return x + y1 - alpha * y2 * y2 + 1


n = 10
alpha = 3

t0 = 0
t1 = 1
t = np.zeros(n + 1)
dt = (t1 - t0) / n
for i in range(n + 1):
    t[i] = t0 + i * dt

U1 = np.zeros(n + 1)
U2 = np.zeros(n + 1)
U3 = np.zeros(n + 1)
U4 = np.zeros(n + 1)

U1[0] = 1
U2[0] = 0.5
U3[0] = 1
U4[0] = 0.5


for i in range(n):
    k1_1 = f1(t[i], U1[i], U2[i])
    k1_2 = f2(t[i], U1[i], U2[i])

    k2_1 = f1(t[i] + dt / 2, U1[i] + 0.5 * dt * k1_1, U2[i] + 0.5 * dt * k1_2)
    k2_2 = f2(t[i] + dt / 2, U1[i] + 0.5 * dt * k1_1, U2[i] + 0.5 * dt * k1_2)

    k3_1 = f1(t[i] + dt / 2, U1[i] + dt * k2_1 / 2, U2[i] + dt * k2_2 / 2)
    k3_2 = f2(t[i] + dt / 2, U1[i] + dt * k2_1 / 2, U2[i] + dt * k2_2 / 2)

    k4_1 = f1(t[i] + dt, U1[i] + dt * k3_1, U2[i] + dt * k3_2)
    k4_2 = f2(t[i] + dt, U1[i] + dt * k3_1, U2[i] + dt * k3_2)

    U1[i + 1] = U1[i] + dt * (k1_1 + 2 * k2_1 + 2 * k3_1 + k4_1) / 6
    U2[i + 1] = U2[i] + dt * (k1_2 + 2 * k2_2 + 2 * k3_2 + k4_2) / 6
    U3[i + 1] = U1[i] + dt * f1(t[i], U3[i], U4[i])
    U4[i + 1] = U2[i] + dt * f2(t[i], U3[i], U4[i])
    U3[i + 1] = U1[i] + 0.5 * dt * (f1(t[i], U3[i], U4[i]) + f1(t[i + 1], U3[i + 1], U4[i + 1]))
    U4[i + 1] = U2[i] + 0.5 * dt * (f2(t[i], U3[i], U4[i]) + f2(t[i + 1], U3[i + 1], U4[i + 1]))


def pretty_print(t, u_1_result, u_2_result, method_name):
    print(method_name)

    print('t\t u_1\t u_2')
    for i in range(len(t)):
        print("%.2f %.4f %.4f" % (t[i], u_1_result[i], u_2_result[i]))
    print()


pretty_print(t, U1, U2, "Runge-Kutta order 4")
pretty_print(t, U3, U4, "Hoine")
