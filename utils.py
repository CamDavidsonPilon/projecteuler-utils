import collections
import functools
import numpy as np
from itertools import chain, combinations
import math

MILLION = 1000000
BILLION = 1000000000
PHI = 0.5 + 0.5 * np.sqrt(5)


def is_int(x):
    return int(x) == x


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


def nCr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)


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

