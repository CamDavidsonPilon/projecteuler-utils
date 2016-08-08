# miller_rabin primality test
import random


def signed_mod(a, n):
    return a % n - n

def _find_ds(n):
    power = 0
    while True:
        if (n-1) % 2**power == 0 and ((n-1) / 2**power) % 2 == 1:
            return (n-1) / 2**power, power
        power +=1


def miller_rabin(n, checks=1):
    d,s = _find_ds(n)
    truth = True
    for i in range(checks):
        a = random.randint(1, min(n,20)-1)
        test_1 = a ** d % n != 1
        test_2 = all(signed_mod(a ** (d*2**r), n) != -1 for r in range(0,s))
        truth *= not (test_1 and test_2)
    return bool(truth)


def is_probably_prime(n, checks=5):
    return miller_rabin(n, 5)