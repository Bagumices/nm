from math import sin, pi, cos, pow
from typing import List, Callable
from pprint import pprint

from numpy.linalg import solve
from pandas import DataFrame


def get_intervals(start: float, end: float, n: int) -> List[List[float]]:
    step = (end - start) / n
    return [[start + (step * i), start + (step * (i+1))] for i in range(0, n)]


def f(x: float, variant_number: int) -> float:
    return (variant_number + x) / 3.5


def p(x: float, variant_number: int) -> float:
    return (4-0.1*x)/(x * x + variant_number / 16)


def q(x: float, variant_number: int) -> float:
    return (x+5)/(x * x + 0.9 * variant_number)


def phi_0(x: float, a: float, b: float) -> float:
    return 2.142857 - 3.285714 * cos((pi/2.0)*((x - a) / (b-a)) + (pi/2.0))


def phi_k(x: float, a: float, b: float, k: int) -> float:
    return cos(-pi * k * ((x - a) / (b - a)) + (pi / 2.0))


def phi_0_diff(x: float, a: float, b: float) -> float:
    return 3.285714 * sin((-pi / 2.0)*((x - a) / (b - a)) + (pi / 2.0))*((-pi/2.0) * (1.0/(b - a)))


def phi_k_diff(x: float, a: float, b: float, k: int):
    return sin((-pi * k * (x - a)) / (b - a) + (pi / 2.0)) * ((pi * k) / (b - a))


def trapezoid_method(
        integrand_function: Callable,
        a: float,
        b: float,
        n: int,
        step: float) -> float:

    x_i = [a + step * i for i in range(0, n)]
    middle_result = 0

    for i in range(1, n):
        middle_result += integrand_function(x_i[i])

    return (step / 2) * (integrand_function(a) + 2 * middle_result + integrand_function(b))


def d_zero(a: float, b: float, variant_number: int, step: float, n: int) -> float:
    def integrand_function(x: float) -> float:
        return p(x, variant_number) * pow(phi_0_diff(x, a, b), 2) + q(x, variant_number) * pow(phi_0(x, a, b), 2) \
               + 2 * f(x, variant_number) * phi_0(x, a, b)

    return trapezoid_method(integrand_function, a, b, n, step)


def b_k(a: float, b: float, variant_number: int, step: float, n: int, k: int) -> float:
    def integrand_function(x: float) -> float:
        return p(x, variant_number) * phi_0_diff(x, a, b) * phi_k_diff(x, a, b, k) + q(x, variant_number) \
               * phi_0(x, a, b) * phi_k(x, a, b, k) + f(x, variant_number) * phi_k(x, a, b, k)

    return trapezoid_method(integrand_function, a, b, n, step)


def a_ik(a: float, b: float, i: int, k: int, variant_number: int, n: int, step: float) -> float:
    def integrand_function(x: float) -> float:
        return p(x, variant_number) * phi_k_diff(x, a, b, k) * phi_k_diff(x, a, b, i) + q(x, variant_number) \
               * phi_k(x, a, b, k) * phi_k(x, a, b, i)

    return trapezoid_method(integrand_function, a, b, n, step)


def ritz_method(a: float, b: float, variant_number: int, step: float, n: int) -> List[float]:
    a_matrix: List[List[float]] = []
    b_matrix: List[float] = []

    # fill A_matrix
    for i in range(n):
        current_row = []

        for k in range(n):
            current_result = a_ik(a, b, i, k, variant_number, n, step)
            # print("with params a={}, b={}, i={}, k={} step={} ==> result = {}".format(a, b, i, k, step, current_result))
            current_row.append(current_result)

        a_matrix.append(current_row)

    # fill B_matrix
    for k in range(n):
        b_matrix.append(b_k(a, b, variant_number, step, n, k))

    print("A matrix")
    print(DataFrame(a_matrix))
    print()
    print("B matrix")
    print(DataFrame(b_matrix))

    # solve linear system
    result = solve(a_matrix, b_matrix)
    print()
    print("C coefficients")
    print(DataFrame(result))
    return result


def main():
    eps = 10e-6
    variant_number = 2
    n = 4
    a = 3/5 - variant_number/13
    b = 2 - variant_number/13
    step = (b - a) / 5
    nu_1 = 15/(variant_number+3)
    nu_2 = -(6*variant_number)/21
    ritz_method(a, b, variant_number, step, n)


if __name__ == '__main__':
    main()
