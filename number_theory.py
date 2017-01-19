from fractions import gcd
from collections import defaultdict
import math
from itertools import count
from prime_numbers import coprime, all_prime_divisors, primesfrom2to
from utils import infinite_product


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
    return a, prevx, prevy


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