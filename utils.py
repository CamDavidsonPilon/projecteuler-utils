import collections
import functools
from itertools import chain, combinations, imap, count
import math
from operator import mul
from fractions import gcd
import random
import numpy as np

MILLION = 1000000
BILLION = 1000000000
PHI = 0.5 + 0.5 * np.sqrt(5)


def product(list):
    return functools.reduce(mul, list, 1)


def is_square(x):
    return is_int(math.sqrt(x))


def is_int(x):
    if math.isnan(x):
        return False
    return int(x) == x


def digital_root(x):
    return int(x - 9*int((x-1)/9))


def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))


class memoized(object):
    '''Decorator. Caches a function's return value each time it is called.
    If called later with the same arguments, the cached value is returned
    (not reevaluated).
    '''
    def __init__(self, func):
        self.func = func
        self.cache = {}

    def __call__(self, *args):
        if not isinstance(args, collections.Hashable):
            # uncacheable. a list, for instance.
            # better to not cache than blow up.
            return self.func(*args)
        if args in self.cache:
            return self.cache[args]
        else:
            value = self.func(*args)
            self.cache[args] = value
            return value

    def __repr__(self):
        '''Return the function's docstring.'''
        return self.func.__doc__

    def __get__(self, obj, objtype):
        '''Support instance methods.'''
        return functools.partial(self.__call__, obj)


def promote_digits(digits):
    """
    see problem 20, 16, 104
    """
    for i in xrange(digits.shape[0]):
        d = digits[i]
        if d >= 10 and digits.shape[0] > i + 1:
            remainder = d % 10
            divisor = d / 10
            digits[i] = remainder
            digits[i+1] += divisor
    return digits


def digits(n):
    return map(int, str(n))


def reverse_number(n):
    reversed_n = 0
    power = int(np.log10(n))
    while power >= 0:
        n, r = divmod(n, 10)
        reversed_n += r * 10 ** power
        power -= 1
    return reversed_n


def to_digit(array):
    return sum(10**(len(array) - i - 1)*int(array[i]) for i in range(len(array)))


def nCr(n, r):
    if r > n:
        return 0
    if r < 0 or n < 0:
        return 0

    f = math.factorial
    return f(n) / f(r) / f(n-r)


def infinite_product(iterx, itery):
    intermediate_x, intermediate_y = [], []
    intercept = 0

    while True:
        try:
            intermediate_y.append(itery.next())
            intermediate_x.append(iterx.next())

            for x in range(intercept+1):
                y = intercept - x
                yield (intermediate_x[x], intermediate_y[y])

            intercept += 1
        except StopIteration:
            break


def random_choice_from_dist(a):
    r = random.random()
    a = np.array(a)
    a = 1.0*a/a.sum()
    running_sum = 0
    i = -1
    while r > running_sum:
        i += 1
        running_sum += a[i]
    return i


def decreasing_elements(array, M, max_length):
    """
    Replaces code like

        M = 20

        for a1 in xrange(M, 0, -1):
            for a2 in xrange(a1, 0, -1):
                for a3 in xrange(a2, 0, -1):
                    for a4 in xrange(a3, 0, -1):
                        yield (a1,a2,a3,a4)

    with

        decreasing_elements([], M, 4)

    """
    if len(array) == max_length:
        yield array
    else:
        for x in range(M+1):
            for _array in decreasing_elements(array + [x], x, max_length):
                yield _array


def iflatmap(func, iterable):
    return chain.from_iterable((func(x) for x in iterable))


def continued_fraction_expansion(x, tol=10e-8):
    # this is breaking due to precision errors. If determine sqrt, use
    # sqrt_continued_fraction_expansion
    while abs(x - int(x)) > tol:
        yield int(x)
        x = 1./(x - int(x))


def sqrt_continued_fraction_expansion(n):
    """
    determines the continued fraction expansion of sqrt(n),
    that is, return ai s.t.

    a1 + 1/(a2 + 1/(a3 + 1/(....))) == sqrt(n)
    """
    m = 0.
    d = 1.
    a0 = a = int(math.sqrt(n))
    yield a
    while True:
        m = d*a - m
        d = (n - m**2) / d
        a = int((a0 + m)/d)
        yield a


def continued_fraction_convergents(x):
    """
    https://en.wikipedia.org/wiki/Continued_fraction#Infinite_continued_fractions_and_convergents

    returns the values
    a1, a1 + 1/a2, a1 + 1/(a2 + (1/a3)), ... ==  x

    """
    hn_1, hn_2 = 1, 0
    kn_1, kn_2 = 0, 1

    for an in continued_fraction_expansion(x):
        hn = an*hn_1 + hn_2
        kn = an*kn_1 + kn_2
        yield (hn, kn)

        hn_1, hn_2 = hn, hn_1
        kn_1, kn_2 = kn, kn_1


def sqrt_continued_fraction_convergents(n):
    """
    https://en.wikipedia.org/wiki/Continued_fraction#Infinite_continued_fractions_and_convergents

    returns the values
    a1, a1 + 1/a2, a1 + 1/(a2 + (1/a3)), ... ==  sqrt(n)

    """
    hn_1, hn_2 = 1, 0
    kn_1, kn_2 = 0, 1

    for an in sqrt_continued_fraction_expansion(n):
        hn = an*hn_1 + hn_2
        kn = an*kn_1 + kn_2
        yield (hn, kn)

        hn_1, hn_2 = hn, hn_1
        kn_1, kn_2 = kn, kn_1


def sqrt_continued_fraction_to_decimal_expansion(n):
    hn_1, hn_2 = 1, 0
    kn_1, kn_2 = 0, 1

    for an in sqrt_continued_fraction_expansion(n):
        hn = an*hn_1 + hn_2
        kn = an*kn_1 + kn_2
        if kn_1 > 0 and hn/kn == hn_1/kn_1:
            d = hn/kn
            yield d

            hn -= d*kn
            hn *= 10

            hn_1 -= d*kn_1
            hn_1 *= 10

        hn_1, hn_2 = hn, hn_1
        kn_1, kn_2 = kn, kn_1


def polygonal_iterator(d):

    if d == 3:
        f = lambda n: n*(n+1)/2
    elif d == 4:
        f = lambda n: n**2
    elif d == 5:
        f = lambda n: n*(3*n-1)/2
    elif d == 6:
        f = lambda n: n*(2*n-1)
    elif d == 7:
        f = lambda n: n*(5*n-3)/2
    elif d == 8:
        f = lambda n: n*(3*n-2)
    else:
        raise ValueError("no function for %d" % d)
    return imap(f, count(start=1))


# http://stackoverflow.com/questions/16344284/how-to-generate-a-list-of-palindrome-numbers-within-a-given-range
def palindrome_number_generator():
    yield 0
    lower = 1
    while True:
        higher = lower*10
        for i in xrange(lower, higher):
            s = str(i)
            yield int(s+s[-2::-1])
        for i in xrange(lower, higher):
            s = str(i)
            yield int(s+s[::-1])
        lower = higher


def palindromes(lower, upper):
    all_palindrome_numbers = palindrome_number_generator()
    for p in all_palindrome_numbers:
        if p >= lower:
            break
    palindrome_list = [p]
    for p in all_palindrome_numbers:
        # Because we use the same generator object,
        # p continues where the previous loop halted.
        if p >= upper:
            break
        palindrome_list.append(p)
    return palindrome_list


def is_terminating(num, denom):
    if gcd(num, denom) != 1:
        denom = denom/gcd(num, denom)

    while denom % 5 == 0:
        denom = denom / 5
    while denom % 2 == 0:
        denom = denom / 2
    return denom == 1


def fast_expon_mod_m(x, n, m):
    if n == 1:
        return x
    elif n == 0:
        return 1
    elif n % 2 == 0:
        return fast_expon_mod_m(x * x % m, n / 2, m)
    elif n % 2 == 1:
        return x * fast_expon_mod_m(x * x % m, (n - 1) / 2, m) % m


def fast_matrix_expon_mod_m(M, n, m):
    """
    Compute M**n mod m for square numpy matrix/vector M

    """
    if n == 0:
        return np.eye(M.shape[0], dtype=int)
    if n % 2 == 1:
        return (M.dot(fast_matrix_expon_mod_m(M, n-1, m))) % m
    else:
        D = fast_matrix_expon_mod_m(M, n/2, m)
        return (D.dot(D)) % m


@memoized
def fast_2matrix_expon_mod_m(M, n, m):
    """
    Compute M**n mod m for square (2x2) list of lists

    """
    def matrix_mul(
                   ((a11, a12), (a21, a22)),
                   ((b11, b12), (b21, b22)),
                  ):
        return (
               (a11*b11 + a12*b21, a11*b12 + a12*b22),
               (a21*b11 + a22*b21, a21*b12 + a22*b22)
               )

    def matrix_mod(((a11, a12), (a21, a22)), m):
        return ((a11 % m, a12 % m), (a21 % m, a22 % m))

    if n == 0:
        return ((1, 0), (0, 1))
    if n % 2 == 1:
        return matrix_mod(matrix_mul(M, fast_2matrix_expon_mod_m(M, n-1, m)), m)
    else:
        D = fast_2matrix_expon_mod_m(M, n/2, m)
        return matrix_mod(matrix_mul(D, D), m)
