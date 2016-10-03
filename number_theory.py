from collections import defaultdict
import math
from itertools import count
from prime_numbers import coprime, all_prime_divisors
from utils import infinite_product

def pythagorean_triples():

    def gen_k():
        for _ in count(start=1):
            yield _

    def gen_m_n():
        for m in count(start=1):
            for n in xrange(1, m):
                if coprime(m, n) and (m - n) % 2 == 1:
                    yield (m, n)

    for (k, (m, n)) in infinite_product(gen_k(), gen_m_n()):
        a = k*(m**2 - n**2)
        b = 2*k*m*n
        c = k*(m**2 + n**2)
        assert a**2 + b**2 == c**2
        yield a, b, c


def primitive_pythagorean_triples():

    def gen_m_n():
        for m in count(start=1):
            for n in xrange(1, m):
                if coprime(m, n) and (m - n) % 2 == 1:
                    yield (m, n)

    for (m, n) in gen_m_n():
        a = (m**2 - n**2)
        b = 2*m*n
        c = (m**2 + n**2)
        assert a**2 + b**2 == c**2
        yield a, b, c


def totient_function(n, prime_cache=None):
    product = n
    for p in all_prime_divisors(n, prime_cache):
        product *= (1 - 1./p)
    return int(round(product))


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