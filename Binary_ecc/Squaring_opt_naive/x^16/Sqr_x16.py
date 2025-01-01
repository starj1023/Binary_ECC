import math
from math import log2

from projectq import MainEngine
from projectq.ops import H, CNOT, Measure, Toffoli, X, All, T, Tdag, S, Swap
from projectq.backends import CircuitDrawer, ResourceCounter, ClassicalSimulator
from projectq.meta import Loop, Compute, Uncompute, Control


def Squraing_opt_naive(eng):

    n = 16

    x1 = eng.allocate_qureg(n)
    x2 = eng.allocate_qureg(n)
    x3 = eng.allocate_qureg(n)

    if (resource_check != 1):
        Round_constant_XOR(eng, x1, 0xabab, n)

    # Compiler-friendly implementation (Comment or uncomment)
    x2 = Squaring_2(eng, x1, x2, n)

    # Naive implementation (Comment or uncomment)
    # x3 = Squaring_2_slow(eng, x1, x3, n)

    if (resource_check != 1):
        print_state(eng, x2, n) # Comment or uncomment
        # print_state(eng, x3, n) # Comment or uncomment


def Squaring_2(eng, input, output, n):

    matrix = [1,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,
            1,1,1,1,0,0,0,0,1,0,0,1,0,0,0,0,
            1,0,1,1,1,1,0,1,1,0,0,0,0,0,0,0,
            0,1,1,0,1,1,0,0,0,0,0,1,0,0,0,0,
            0,1,1,0,0,0,1,0,1,0,1,0,0,0,1,0,
            1,0,1,0,1,1,0,0,0,0,1,1,0,0,0,0,
            0,0,1,1,1,0,1,1,1,0,0,0,0,0,0,0,
            0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,
            1,1,1,0,0,1,0,0,0,1,0,0,0,1,0,0,
            0,1,0,1,1,0,0,0,0,1,1,0,0,0,0,0,
            0,1,1,1,0,1,1,1,0,0,0,0,0,0,0,0,
            0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,
            1,1,0,0,1,0,0,0,1,0,0,0,1,0,0,0,
            1,0,1,1,0,0,0,0,1,1,0,0,0,0,0,0,
            1,1,1,0,1,1,1,0,0,0,0,0,0,0,0,0,
            0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0]

    for i in range(n):
        for j in range(n):
            if (matrix[n * j + ((n - 1) - j - i) % n] == 1):
                CNOT | (input[(i + j) % n], output[j])

    return output

def Squaring_2_slow(eng, input, output, n):

    matrix = [1,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,
            1,1,1,1,0,0,0,0,1,0,0,1,0,0,0,0,
            1,0,1,1,1,1,0,1,1,0,0,0,0,0,0,0,
            0,1,1,0,1,1,0,0,0,0,0,1,0,0,0,0,
            0,1,1,0,0,0,1,0,1,0,1,0,0,0,1,0,
            1,0,1,0,1,1,0,0,0,0,1,1,0,0,0,0,
            0,0,1,1,1,0,1,1,1,0,0,0,0,0,0,0,
            0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,
            1,1,1,0,0,1,0,0,0,1,0,0,0,1,0,0,
            0,1,0,1,1,0,0,0,0,1,1,0,0,0,0,0,
            0,1,1,1,0,1,1,1,0,0,0,0,0,0,0,0,
            0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,
            1,1,0,0,1,0,0,0,1,0,0,0,1,0,0,0,
            1,0,1,1,0,0,0,0,1,1,0,0,0,0,0,0,
            1,1,1,0,1,1,1,0,0,0,0,0,0,0,0,0,
            0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0]

    for i in range(n):
        for j in range(n):
            if (matrix[n * (i + 1) - 1 - j] == 1):
                CNOT | (input[j], output[i])

    return output

def Squaring_plus_input(eng, input, output, n):
    matrix = [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
              1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,
              1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0,
              1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0,
              0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0,
              0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0,
              0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,
              0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0,
              1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0,
              0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
              0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
              0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
              0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    for i in range(n):
        for j in range(n):
            if (matrix[n * j + ((n - 1) - j - i) % n] == 1):
                CNOT | (input[(i + j) % n], output[j])

    return output

def Round_constant_XOR(eng, k, rc, bit):
    for i in range(bit):
        if (rc >> i & 1):
            X | k[i]

def print_state(eng, b, n):
    All(Measure) | b
    print('Result : ', end='')
    for i in range(n):
        print(int(b[n-1-i]), end='')
    print()

global resource_check
global AND_check
global NCT

NCT = 1
resource_check = 0
classic = ClassicalSimulator()
eng = MainEngine(classic)
Squraing_opt_naive(eng)
eng.flush()

resource_check = 1
NCT = 1
AND_check = 0 # If AND = 1, the number of qubits must be manually counted.
Resource = ResourceCounter()
eng = MainEngine(Resource)
Squraing_opt_naive(eng)
print('\n')
print(Resource)