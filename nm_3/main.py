from typing import List, Tuple

from pandas import DataFrame


def generate_t_grid(m: int, t: float) -> List[float]:
    tau = t / (m-1)
    return [j * tau for j in range(m)]


def generate_x_grid(n: int, l: float) -> List[float]:
    h = l / (n-1)
    return [i * h for i in range(n)]


def section_method(n: int, m: int, x_grid: List[float], t_grid: List[float], tau: float, h: float) -> List[List[float]]:
    def calculate_border_conditions() -> List[List[float]]:
        # u(0, t) = 0, u(1, t) = t
        for i in range(m):
            matrix[i][n-1] = t_grid[i]
        return matrix

    def f_i(x_i: float, t_j: float, y_t_previous: float) -> float:
        return tau * (pow(x_i, 2) - 2 * t_j) + y_t_previous

    def calculate_driving_factors(
            current_t_layer: int,
            previous_y_layer: List[float]) -> Tuple[List[float], List[float]]:
        a_i = tau/pow(h, 2)
        b_i = a_i
        c_i = 1 + (2 * a_i)

        psi_result = [0.0 for i in range(n)]
        eta_result = [0.0 for i in range(n)]

        psi_next = 0
        eta_next = t_grid[current_t_layer]

        for i in reversed(range(1, n-1)):
            psi_result[i] = a_i / (c_i - b_i * psi_next)
            eta_result[i] = (eta_next * b_i + f_i(x_grid[i], t_grid[current_t_layer], previous_y_layer[i])) \
                            / (c_i - b_i * psi_next)

            psi_next = psi_result[i]
            eta_next = eta_result[i]

        return psi_result, eta_result

    def calculate_y() -> List[List[float]]:
        for j in range(1, m):
            psi_list, eta_list = calculate_driving_factors(j, matrix[j-1])
            for i in range(0, n-2):
                matrix[j][i+1] = psi_list[i+1] * matrix[j][i] + eta_list[i+1]

        return matrix

    matrix: List[List[float]] = [[0.0 for i in range(n)] for j in range(m)]
    matrix = calculate_border_conditions()
    matrix = calculate_y()

    return matrix


def print_matrix(matrix: List[List[float]], x_grid: List[float], t_grid: List[float], n: int, m: int):
    for j in reversed(range(m)):
        print("t(%.2f)" % t_grid[j], end=' ')
        for i in range(n):
            print("%.4f" % matrix[j][i], end='    ')
        print()

    print(end='        ')
    for i in range(n):
        print("x(%.2f)" % x_grid[i], end='   ')


def main():
    n = 11
    m = 6
    l = 1
    tau = 0.05 / m
    h = l / n
    t = 0.05

    x_grid = generate_x_grid(n, l)
    t_grid = generate_t_grid(m, t)
    matrix = section_method(n, m, x_grid, t_grid, tau, h)
    print_matrix(matrix, x_grid, t_grid, n, m)


if __name__ == '__main__':
    main()
