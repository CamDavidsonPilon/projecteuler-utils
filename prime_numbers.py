# prime_numbers.py
import fractions
from math import sqrt
from operator import mul
import numpy as np
from collections import defaultdict
from primality_tests import is_probably_prime
from utils import product



def eratosthenes(limit):
    # naive and slow, use primesfrom2to
    numbers_up_to_limit, still_valid = range(2, limit+1), [True]*(limit-1)
    assert len(numbers_up_to_limit) == len(still_valid)

    for i, n in enumerate(numbers_up_to_limit):
        is_prime = still_valid[i]
        if is_prime:
            for j in xrange(i+n, limit-1, n):
                still_valid[j] = False
    return [numbers_up_to_limit[i] for i in xrange(limit-1) if still_valid[i]]


def primesfrom2to(n):
    # http://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n-in-python/3035188#3035188
    """ Input n>=1 Returns a array of primes, 2 <= p < n """
    if n == 1:
        return []
    elif n == 2:
        return []
    elif n == 3:
        return [2]
    elif n == 4:
        return [2, 3]
    elif n == 5:
        return [2, 3]
    sieve = np.ones(n/3 + (n % 6 == 2), dtype=np.bool)
    sieve[0] = False
    for i in xrange(int(n**0.5)/3+1):
        if sieve[i]:
            k = 3 * i + 1 | 1
            sieve[      ((k*k)/3)      ::2*k] = False
            sieve[(k*k+4*k-2*k*(i&1))/3::2*k] = False
    return map(int, np.r_[2, 3, ((3*np.nonzero(sieve)[0]+1) | 1)])


def primal_decomposition(n, prime_cache=None):

    original_n = n
    output = []

    if prime_cache is None:
        primes = primesfrom2to(int(n/2) + 2)
    else:
        primes = prime_cache

    for p in primes:
        # if using a prime cache, I can stop after int(n/2) assuming the cache is sorted
        if p > int(original_n/2)+1:
            break

        count = 0
        while n % p == 0:
            count += 1
            n = n / p

        if count > 0:
            output.append((p, count))

        if n == 1:
            break

    if output == []:
        # n is a prime then
        output = [(n, 1)]
    return output


def pair_divisors(n):
    """
    This can be done faster. I don't _need_ to check all the numbers below sqrt(n) in some cases
    like 2^60 - 1 = 3^2 * 5^2 * 7 * 11 * 13 * 31 * 41 * 61 * 151 * 331 * 1321

    """
    primes_that_factor_into_n = []
    max_power_of_prime = 1

    yield (1, n)
    for i in xrange(2, int(sqrt(n))+1):
        if n % i == 0:
            yield (i, n/i)

def all_divisors(n):
    for d1, d2 in pair_divisors(n):
        yield d1
        if d1 != d2:
            yield d2


def all_prime_divisors(n, prime_cache=None):
    if prime_cache is None:
        prime_cache = primesfrom2to(n + 1)
    for p in prime_cache:
        if n % p == 0:
            n = n / p
            yield p
        if n == 1:
            raise StopIteration()


def radical(n, prime_cache=None):
    if n == 1:
        return 1
    return reduce(mul, all_prime_divisors(n, prime_cache))


def divisor_0(n, prime_cache=None):
    decomp = primal_decomposition(n, prime_cache)
    return reduce(lambda x, y: x*y, map(lambda (p, a): a+1, decomp), 1)


def divisor_k(n, k, prime_cache=None):
    assert k > 0
    if n == 1:
        return 1
    if n == 0:
        return 0
    product = 1
    decomp = primal_decomposition(n, prime_cache)
    for p, a in decomp:
        product *= (p**(k*a + k) - 1) / (p**k-1)
    return product


def divisor_k_lookup(up_to, k):
    """
    Creates a cache for looking up divisor_k.
    Uses an n**2 algorithm to create this.
    """
    div = defaultdict(lambda: 1)
    div[1] = 1

    for i in xrange(2, up_to):
        for j in xrange(i, up_to, i):
            div[j] += i**k

    return div


def coprime(m, n):
    return fractions.gcd(m, n) == 1


def divisors(factors):
    """
    Generates all divisors, unordered, from the prime factorization.
    From https://numericalrecipes.wordpress.com/2009/05/13/divisors/
    """
    ps = sorted(set(factors))
    omega = len(ps)

    def rec_gen(n=0):
        if n == omega:
            yield 1
        else:
            pows = [1]
            for j in xrange(factors.count(ps[n])):
                pows += [pows[-1] * ps[n]]
            for q in rec_gen(n + 1):
                for p in pows:
                    yield p * q

    for p in rec_gen():
        yield p
