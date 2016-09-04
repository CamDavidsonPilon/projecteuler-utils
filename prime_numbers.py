# prime_numbers.py
import fractions
from math import sqrt
import numpy as np
from collections import defaultdict


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
    return np.r_[2, 3, ((3*np.nonzero(sieve)[0]+1) | 1)]


def primal_decomposition(n, prime_cache=None):

    original_n = n
    output = []

    if prime_cache is None:
        primes = primesfrom2to(int(n/2) + 2)
    else:
        primes = prime_cache

    for p in primes:
        # if using a prime cache, I can stop after int(n/2)
        #if p > int(original_n/2)+1:
        #    break
        # nope, can't assume they are sorted

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
    yield (1, n)
    for i in xrange(2, int(sqrt(n))+1):
        if n % i == 0:
            yield (i, n/i)


def all_divisors(n):
    for d1, d2 in pair_divisors(n):
        yield d1
        if d1 != d2:
            yield d2


def all_prime_divisors(n):
    primes = set(primesfrom2to(n + 1)) # this is shitty. I need to really only check if n is a prime, not primes between n/2 and n.
    for d in all_divisors(n):
        if d in primes:
            yield d


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
    Uses an efficient algorithm to create this. 
    """
    div = defaultdict(lambda: 0)
    for i in xrange(1, up_to):
        for j in xrange(i, up_to, i):
            div[j] += i**k

    return div


def totient_function(n):
    product = n
    for p in all_prime_divisors(n):
        product *= (1 - 1./p)
    return int(round(product))


def coprime(m, n):
    return fractions.gcd(m, n) == 1

