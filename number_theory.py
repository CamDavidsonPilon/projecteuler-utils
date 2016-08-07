from prime_numbers import coprime, all_prime_divisors
from itertools import count
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


def totient_function(n):
    product = n
    for p in all_prime_divisors(n):
        product *= (1 - 1./p)
    return int(round(product))
