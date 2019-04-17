from math import pow
from typing import List, Tuple

from pandas import DataFrame


def generate_grid(h: float, tau: float, n: int) -> List[Tuple[float, float]]:
    return [(i * h, j * tau) for i in range(n+1) for j in range(n+1)]


def alpha(x: float) -> float:
    return 1/(1+x)


def d2_alpha_dx2(x: float) -> float:
    return 2/pow(1+x, 3)


def beta(x: float) -> float:
    return -0.7/pow(1+x, 2)


def y_0_j(t_j: float) -> float:
    return 1/(0.7*t_j + 1)


def y_1_j(t_j: float) -> float:
    return 1/(0.7*t_j + 2)


def y_i_0(x_i: float) -> float:
    return alpha(x_i)


def f(x: float, t: float) -> float:
    return 2 * (pow(0.7, 2) - 1)/pow(x + 0.7*t + 1, 3)


def y_i_1(x_i: float, tau: float) -> float:
    return alpha(x_i) + tau * beta(x_i) + pow(tau, 2)/2 * (d2_alpha_dx2(x_i) + f(x_i, 0))


def grid_method(
        x_sections: List[float],
        y_sections: List[float],
        s: float,
        n: int,
        tau: float):
    def calculate_border_condition() -> List[List[float]]:
        for i, y in enumerate(y_sections):
            matrix[i][0] = y_0_j(y)
            matrix[i][n] = y_1_j(y)

        return matrix

    def calculate_first_two_line() -> List[List[float]]:
        for i in range(1, n):
            matrix[0][i] = y_i_0(x_sections[i])
            matrix[1][i] = y_i_1(x_sections[i], tau)

        return matrix

    matrix: List[List[float]] = [[0 for i in range(n+1)] for i in range(n+1)]
    matrix = calculate_border_condition()
    matrix = calculate_first_two_line()

    for i in range(1, n):
        for j in range(1, n):
            matrix[i+1][j] = s * matrix[i][j+1] + 2 * (1 - s) * matrix[i][j] + s * matrix[i][j-1] - matrix[i-1][j] + \
                             pow(tau, 2) * f(x_sections[j], y_sections[i])

    print(DataFrame(matrix))


def main():
    n = 10
    l = 1
    t = 1
    h = l/n
    tau = t/n
    s = pow(tau, 2)/pow(h, 2)

    # grid = generate_grid(h, tau, n)
    # section_method(grid, tau, h, s, n)
    x_sections = [i * h for i in range(n+1)]
    y_sections = [i * tau for i in range(n+1)]
    grid_method(x_sections, y_sections, s, n, tau)


if __name__ == '__main__':
    main()
