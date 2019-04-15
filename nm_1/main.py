from math import sin, pi, cos, pow
from typing import List, Callable, Tuple

from numpy.linalg import solve
from pandas import DataFrame


def get_intervals(start: float, end: float, n: int) -> List[float]:
    step = (end - start) / (n-1)
    return [start + (step * i) for i in range(0, n)]


def f(x: float, variant_number: int) -> float:
    return (variant_number + x) / 3.5


def p(x: float, variant_number: int) -> float:
    return (4-0.1*x)/(x * x + variant_number / 16)


def q(x: float, variant_number: int) -> float:
    return (x+5)/(x * x + 0.9 * variant_number)


def phi_0(x: float, a: float, b: float, nu_1: float, nu_2: float) -> float:
    return nu_1 - (nu_2-nu_1)*cos((pi/2.0)*((x - a) / (b-a)) + (pi/2.0))


def phi_k(x: float, a: float, b: float, k: int) -> float:
    return cos(-pi * k * ((x - a) / (b - a)) + (pi / 2.0))


def phi_0_diff(x: float, a: float, b: float, nu_1: float, nu_2: float) -> float:
    return (nu_2-nu_1) * sin((-pi / 2.0)*((x - a) / (b - a)) + (pi / 2.0))*((-pi/2.0) * (1.0/(b - a)))


def phi_k_diff(x: float, a: float, b: float, k: int):
    return sin((-pi * k * (x - a)) / (b - a) + (pi / 2.0)) * ((pi * k) / (b - a))


def trapezoid_method(
        integrand_function: Callable,
        a: float,
        b: float,
        integral_split_count: int,
        step: float) -> float:

    x_i = [a + step * i for i in range(0, integral_split_count)]
    middle_result = 0

    for i in range(1, integral_split_count):
        middle_result += integrand_function(x_i[i])

    return (step / 2) * (integrand_function(a) + 2 * middle_result + integrand_function(b))


def d_zero(
        a: float,
        b: float,
        variant_number: int,
        step: float,
        integral_split_count: int,
        nu_1: float,
        nu_2: float) -> float:
    def integrand_function(x: float) -> float:
        return p(x, variant_number) * pow(phi_0_diff(x, a, b, nu_1, nu_2), 2) + q(x, variant_number) \
               * pow(phi_0(x, a, b, nu_1, nu_2), 2) + 2 * f(x, variant_number) * phi_0(x, a, b, nu_1, nu_2)

    return trapezoid_method(integrand_function, a, b, integral_split_count, step)


def b_k(
        a: float,
        b: float,
        variant_number: int,
        step: float,
        integral_split_count: int,
        k: int,
        nu_1: float,
        nu_2: float) -> float:
    def integrand_function(x: float) -> float:
        return p(x, variant_number) * phi_0_diff(x, a, b, nu_1, nu_2) * phi_k_diff(x, a, b, k) + q(x, variant_number) \
               * phi_0(x, a, b, nu_1, nu_2) * phi_k(x, a, b, k) + f(x, variant_number) * phi_k(x, a, b, k)

    return trapezoid_method(integrand_function, a, b, integral_split_count, step)


def a_ik(a: float, b: float, i: int, k: int, variant_number: int, integral_split_count: int, step: float) -> float:
    def integrand_function(x: float) -> float:
        return p(x, variant_number) * phi_k_diff(x, a, b, k) * phi_k_diff(x, a, b, i) + q(x, variant_number) \
               * phi_k(x, a, b, k) * phi_k(x, a, b, i)

    return trapezoid_method(integrand_function, a, b, integral_split_count, step)


def ritz_method(
        a: float,
        b: float,
        variant_number: int,
        step: float,
        integral_split_count: int,
        system_function_count: int,
        nu_1: float,
        nu_2: float) -> Tuple[List[List[float]], List[float], List[float]]:

    a_matrix: List[List[float]] = []
    b_matrix: List[float] = []

    for i in range(1, system_function_count+1):
        current_row = []

        for k in range(1, system_function_count+1):
            current_row.append(a_ik(a, b, i, k, variant_number, integral_split_count, step))
        a_matrix.append(current_row)

    for k in range(1, system_function_count+1):
        b_matrix.append(b_k(a, b, variant_number, step, integral_split_count, k, nu_1, nu_2))

    return a_matrix, b_matrix, solve(a_matrix, b_matrix)


def u(x: float, coefficients: List[float], a: float, b: float, nu_1: float, nu_2: float) -> float:
    result = 0
    for i, c in enumerate(coefficients, start=1):
        result += c * phi_k(x, a, b, i)
    return phi_0(x, a, b, nu_1, nu_2) + result


def calculate_u_in_interval(
        coefficients: List[float],
        a: float,
        b: float,
        section_split: int,
        nu_1: float,
        nu_2: float) -> Tuple[List[float], List[float]]:
    intervals = get_intervals(a, b, section_split)
    result = []

    for point in intervals:
        result.append(u(point, coefficients, a, b, nu_1, nu_2))

    return intervals, result


def main():
    variant_number = 2
    system_function_count = 4
    section_split = 5
    integral_split_count = 5

    a = 3/5 - variant_number/13
    b = 2 - variant_number/13
    step = (b - a) / integral_split_count

    nu_1 = 15/(variant_number+3)
    nu_2 = -(6*variant_number)/21

    a_matrix, b_matrix, result = ritz_method(a, b, variant_number, step, integral_split_count, system_function_count,
                                             nu_1, nu_2)

    # condition number = 44.2619

    print("A matrix")
    print(DataFrame(a_matrix))

    print("\nB matrix")
    print(DataFrame(b_matrix))

    print("\nC coefficients")
    print(DataFrame(result))

    print("\nnu_1 = {}, nu_2 = {}, a = {}, b = {}\n".format(nu_1, nu_2, a, b))

    interval, result = calculate_u_in_interval(result, a, b, section_split, nu_1, nu_2)
    for i in range(len(interval)):
        print("u({}) = {}".format(interval[i], result[i]))


if __name__ == '__main__':
    main()
