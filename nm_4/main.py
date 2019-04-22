from typing import List, Tuple
from math import sin, pow


def generate_grid(start: int, end: int, h: float) -> List[float]:
    result = [start]
    i = 1
    while start + i * h <= end:
        result.append(start + i * h)
        i += 1
    return result


def f1(t: float, u_1: float, u_2: float) -> float:
    return sin(3 * pow(u_1, 2)) + t + u_2


def f2(t: float, u_1: float, u_2: float):
    return t + u_1 - 3 * pow(u_2, 2) + 1


def euler_method(
        grid: List[float],
        h: float,
        u_1_0: float,
        u_2_0: float) -> Tuple[List[float], List[float]]:
    u_1_result: List[float] = [u_1_0]
    u_2_result: List[float] = [u_2_0]
    grid_len = len(grid)

    for i in range(grid_len):
        u_1 = u_1_result[i-1]
        u_2 = u_2_result[i-1]

        u_1_next = u_1 + h * u_1
        u_2_next = u_2 + h * u_2

        u_1_result.append(u_1_next)
        u_2_result.append(u_2_next)

    return u_1_result, u_2_result


def runge_kutta_method(grid: List[float], u_1_0: float, u_2_0: float, h: float) -> Tuple[List[float], List[float]]:
    tau = h / 2
    u_1_result = [u_1_0]
    u_2_result = [u_2_0]

    for i in range(len(grid)):
        # 2 order
        # k_1_u_1 = u_1_result[i]
        # k_2_u_1 = du_dt_1(grid[i] + tau, u_1_result[i] + k_1_u_1, u_2_result[i] + k_1_u_1)
        #
        # k_1_u_2 = u_2_result[i]
        # k_2_u_2 = du_dt_2(grid[i] + tau, u_1_result[i] + k_1_u_2, u_2_result[i] + k_1_u_2)
        #
        # u_1_next = u_1_result[i] + 0.5 * tau * (k_1_u_1 + k_2_u_1)
        # u_2_next = u_2_result[i] + 0.5 * tau * (k_1_u_2 + k_2_u_2)
        #
        # u_1_result.append(u_1_next)
        # u_2_result.append(u_2_next)
        # 4 order
        k_1_u_1 = u_1_result[i]
        k_2_u_1 = f1(grid[i] + tau / 4, u_1_result[i] + (tau * k_1_u_1) / 4, u_2_result[i] + (tau * k_1_u_1) / 4)
        k_3_u_1 = f1(grid[i] + tau / 2, u_1_result[i] + (tau + k_2_u_1) / 2, u_2_result[i] + (tau * k_2_u_1) / 2)
        k_4_u_1 = f1(grid[i] + tau, u_1_result[i] + tau * k_1_u_1 - 2 * tau * k_2_u_1 + 2 * tau * k_3_u_1,
                     u_2_result[i] + tau * k_1_u_1 - 2 * tau * k_2_u_1 + 2 * tau * k_3_u_1)

        k_1_u_2 = u_2_result[i]
        k_2_u_2 = f2(grid[i] + tau / 4, u_1_result[i] + (tau * k_1_u_2) / 4, u_2_result[i] + (tau * k_1_u_2) / 4)
        k_3_u_2 = f2(grid[i] + tau / 2, u_1_result[i] + (tau + k_2_u_2) / 2, u_2_result[i] + (tau * k_2_u_2) / 2)
        k_4_u_2 = f2(grid[i] + tau, u_1_result[i] + tau * k_1_u_2 - 2 * tau * k_2_u_2 + 2 * tau * k_3_u_2,
                     u_2_result[i] + tau * k_1_u_2 - 2 * tau * k_2_u_2 + 2 * tau * k_3_u_2)

        u_1_next = tau/6 * (k_1_u_1 + 4*k_3_u_1 + k_4_u_1) + u_1_result[i]
        u_2_next = tau/6 * (k_1_u_2 + 4*k_3_u_2 + k_4_u_2) + u_2_result[i]

        u_1_result.append(u_1_next)
        u_2_result.append(u_2_next)

    return u_1_result, u_2_result


def choine_method(grid: List[float], u_1_0: float, u_2_0: float, h: float) -> Tuple[List[float], List[float]]:
    tau = h/2
    u_1_result = [u_1_0]
    u_2_result = [u_2_0]

    for i in range(len(grid)):
        u_1_i = f1(grid[i], u_1_result[i], u_2_result[i])
        u_2_i = f2(grid[i], u_1_result[i], u_2_result[i])

        u_1_next = u_1_result[i] + 0.5 * tau * u_1_result[i] + 0.5 * tau * u_1_i + 0.5 * pow(tau, 2) * u_1_i
        u_2_next = u_2_result[i] + 0.5 * tau * u_2_result[i] + 0.5 * tau * u_2_i + 0.5 * pow(tau, 2) * u_2_i

        u_1_result.append(u_1_next)
        u_2_result.append(u_2_next)

    return u_1_result, u_2_result


def pretty_print(grid: List[float], u_1_result: List[float], u_2_result: List[float], method_name: str):
    print(method_name)

    print('t\t u_1\t u_2')
    for i in range(len(grid)):
        print("%.2f %.4f %.4f" % (grid[i], u_1_result[i], u_2_result[i]))
    print()


def main():
    h = 0.1
    start = 0
    end = 1
    u_1_0 = 1
    u_2_0 = 0.5

    grid = generate_grid(start=start, end=end, h=h)

    u_1_result, u_2_result = euler_method(grid=grid, h=h, u_1_0=u_1_0, u_2_0=u_2_0)
    pretty_print(grid, u_1_result, u_2_result, "Euler method")

    u_1_result, u_2_result = runge_kutta_method(grid=grid, u_1_0=u_1_0, u_2_0=u_2_0, h=h)
    pretty_print(grid, u_1_result, u_2_result, "Runge-Kutta method")

    u_1_result, u_2_result = choine_method(grid=grid, u_1_0=u_1_0, u_2_0=u_2_0, h=h)
    pretty_print(grid, u_1_result, u_2_result, "Choine method")


if __name__ == '__main__':
    main()
