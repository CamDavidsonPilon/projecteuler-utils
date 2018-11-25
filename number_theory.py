from fractions import gcd
from collections import defaultdict
import math
from itertools import count
import numpy as np
from prime_numbers import coprime, all_prime_divisors, primesfrom2to
from utils import infinite_product, PHI, is_int, fast_2matrix_expon_mod_m


def pythagorean_triples():
    """
    returns (a,b,c) s.t. a**2 + b**2 == c**2

    """
    def _gen_k():
        for k in count(start=1):
            yield k

    def _gen_m_n():
        for m in count(start=1):
            if m % 2 == 0:
                gen_n = xrange(1, m, 2)
            else:
                gen_n = xrange(1, m)
            for n in gen_n:
                if (m - n) % 2 == 1 and coprime(m, n):
                    yield (m, n)

    for (k, (m, n)) in infinite_product(_gen_k(), _gen_m_n()):
        a = k*(m**2 - n**2)
        b = 2*k*m*n
        c = k*(m**2 + n**2)
        assert a**2 + b**2 == c**2  # yess so glad I wrote this assert - self documenting
        yield a, b, c


def primitive_pythagorean_triples():

    def _gen_m_n():
        for m in count(start=1):
            if m % 2 == 0:
                gen_n = xrange(1, m, 2)
            else:
                gen_n = xrange(1, m)
            for n in gen_n:
                if (m - n) % 2 == 1 and coprime(m, n):
                    yield (m, n)

    for (m, n) in _gen_m_n():
        a = (m**2 - n**2)
        b = 2*m*n
        c = (m**2 + n**2)
        assert a**2 + b**2 == c**2
        yield a, b, c


def totient_function(n, prime_cache=None):
    if n % 2 == 0:
        return 2 * totient_function(n / 2) if (n / 2) % 2 == 0 else totient_function(n / 2)

    numerator = 1
    denominator = 1
    for p in all_prime_divisors(n, prime_cache):
        numerator *= p - 1
        denominator *= p
    return n * numerator / denominator


class TotientDict(dict):
    def __init__(self, max_value, *args):
        dict.__init__(self, args)
        self.max_value = max_value

    def __getitem__(self, key):
        if key > self.max_value:
            raise KeyError

        mod4 = key % 4
        if key == 1:
            val = 1
        elif mod4 == 2:
            val = self[key/2]
        elif mod4 == 0:
            val = 2 * self[key/2]
        else:
            val = dict.__getitem__(self, key)
        return val


def totient_lookup(up_to):
    PRIMES = primesfrom2to(up_to)
    totients = TotientDict(up_to)
    totients[1] = 1

    for p in PRIMES:
        if p == 2:
            continue  # TotientDict handles this
        for pn in xrange(p, up_to, 2*p):
            totients[pn] = totients.get(pn, pn)/p * (p - 1)

    return totients


def lcm(*args):
    def _lcm(a, b):
        return a * b // gcd(a, b)
    return reduce(_lcm, args)


def mobius_lookup(up_to):
    sqrt_up_to = int(math.sqrt(up_to))
    mu = defaultdict(lambda: 1)
    mu[1] = 1
    for i in xrange(2, sqrt_up_to+1):
        if mu[i] == 1:
            for j in xrange(i, up_to, i):
                mu[j] *= -i
            for j in xrange(i**2, up_to, i**2):
                mu[j] = 0

    for i in xrange(2, up_to):
        if mu[i] == i:
            mu[i] = 1
        elif mu[i] == -i:
            mu[i] = -1
        elif mu[i] < 0:
            mu[i] = 1
        elif mu[i] > 0:
            mu[i] = -1

    return mu


def xgcd(a, b):
    """Extended GCD:
    Returns (gcd, x, y) where gcd is the greatest common divisor of a and b
    with the sign of b if b is nonzero, and with the sign of a if b is 0.
    The numbers x,y are such that gcd = ax+by."""
    prevx, x = 1, 0
    prevy, y = 0, 1
    while b:
        q, r = divmod(a, b)
        x, prevx = prevx - q*x, x
        y, prevy = prevy - q*y, y
        a, b = b, r
    gcd_ = a
    return gcd_, prevx, prevy


def modular_multiplicate_inverse(n, p):
    # solves n*x = 1 mod(p)
    # ex: 3*x = 1 mod 5 return x = 2
    assert gcd(n, p) == 1, 'inputs must be coprime or no solution exists.'
    sol = xgcd(n, p)[1]
    if sol < 0:
        return p + sol
    else:
        return sol


def chinese_remainder_solver(input):
    """
    Finds the unique solution to
    x = a1 mod(m1)
    x = a2 mod(m2)
    ...
    x = an mod(mn)

    where m1,m2,.. are pairwise coprime

    input is a list of the form [(a1, m1), (a2, m2), ...]
    returns x, lcm(m1,m2,...)
    """
    def binary_chinese_remainder_solver((a1, m1), (a2, m2)):
        (_gcd, n1, n2) = xgcd(m1, m2)
        assert _gcd == 1, "m1 and m2 should be coprime (gcd == 1)"
        return (a1*n2*m2 + a2*m1*n1, m1*m2)

    sol, lcm = reduce(binary_chinese_remainder_solver, input)
    return sol % lcm, lcm


def linear_congruence_solver(a, b, m):
    """
    solves ax = b mod m
    """
    def solutions(sol_mod_m, m):
        while True:
            yield sol_mod_m
            sol_mod_m += m

    g = gcd(a, m)
    if g == 1:
        inverse_a = modular_multiplicate_inverse(a, m)
        return solutions(b * inverse_a % m, m)
    elif b % g == 0:
        return linear_congruence_solver(a/g, b/g, m/g)
    else:
        return iter([])


class Fibonacci():

    def __init__(self):
        self._cache = {}
        self._n = 0
        self._fib_generator = self.fib_generator()

    def fib(self, k):
        if k in self._cache:
            return self._cache[k]

        for fib in self._fib_generator:
            self._n += 1
            self._cache[self._n] = fib

            if self._n == k:
                break

        return fib

    @staticmethod
    def fib_generator():
        yield 1
        yield 1

        fib_1 = fib_2 = 1

        while True:
            fib = fib_1 + fib_2
            yield fib

            fib_2 = fib_1
            fib_1 = fib

    @staticmethod
    def fib_pair_generator():
        yield (1, 1)

        fib_1 = fib_2 = 1

        while True:
            fib = fib_1 + fib_2
            yield (fib_1, fib)

            fib_2 = fib_1
            fib_1 = fib

    def index(self, n):
        v = np.log(n * np.sqrt(5) + 0.5)/np.log(PHI)
        # for large values the above becomes unstable
        if abs(v - np.round(v)) < 1e-8:
            return int(np.round(v))
        else:
            return int(np.floor(v))

    def find_largest_fib_below_n(self, n):
        return self.fib(self.index(n))

    def zeckendorf(self, n):
        if n == 0:
            return []
        else:
            largest_fib_below_n = self.find_largest_fib_below_n(n)
            return [largest_fib_below_n] + self.zeckendorf(n - largest_fib_below_n)

    def zeckendorf_digit(self, n):
        base = ['0'] * (self.index(n) - 1)
        zeckendorf_fibs = self.zeckendorf(n)
        for fib in zeckendorf_fibs:
            base[self.index(n) - self.index(fib)] = '1'
        return ''.join(base)

    def zeckendorf_digit_to_decimal(self, z):
        running_sum = 0
        for i, char in enumerate(reversed(z), start=1):
            if char == '1':
                running_sum += self.fib(i+1)
        return running_sum

    def fib_mod_m(self, k, mod):
        """
        Can compute arbitrarily large Fib numbers, mod m, using
        fast matrix multiplication. Backed by a cache too.
        """
        FIBMATRIX = ((1, 1), (1, 0))
        return fast_2matrix_expon_mod_m(FIBMATRIX, k, mod)[0][1]


def linear_diophantine_solver(a, b, c):
    """
    solves a*x + b*y = c for (x,y)
    If a single solution exists, then an infinite number of solutions exists, indexed
    by an integer k. This functions returns a _function_ that accepts k
    """

    class NoSolution(Exception):
        pass

    if c % gcd(a, b) != 0:
        raise NoSolution()

    # find single solution
    gcd_, x, y = xgcd(a, b)
    x = x * abs(c)
    y = y * abs(c)

    u, v = a / gcd_, b / gcd_
    return lambda k: (x + k * v, y - k * u)


def diophantine_count(a, n):
    # from https://math.stackexchange.com/questions/30638/count-the-number-of-positive-solutions-for-a-linear-diophantine-equation
    """Computes the number of nonnegative solutions (x) of the linear
    Diophantine equation
        a[0] * x[0] + ... a[N-1] * x[N-1] = n

    Theory: For natural numbers a[0], a[2], ..., a[N - 1], n, and j,
    let p(a, n, j) be the number of nonnegative solutions.

    Then one has:
        p(a, m, j) = sum p(a[1:], m - k * a[0], j - 1), where the sum is taken
        over 0 <= k <= floor(m // a[0])

    Examples
    --------
    >>> diophantine_count([3, 2, 1, 1], 47)
    3572
    >>> diophantine_count([3, 2, 1, 1], 40)
    2282
    """

    def p(a, m, j):
        if j == 0:
            return int(m == 0)
        else:
            return sum([p(a[1:], m - k * a[0], j - 1)
                        for k in xrange(1 + m // a[0])])

    return p(a, n, len(a))
