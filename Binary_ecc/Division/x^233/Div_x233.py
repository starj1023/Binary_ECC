import math
from math import log2

from projectq import MainEngine
from projectq.ops import H, CNOT, Measure, Toffoli, X, All, T, Tdag, S, Swap
from projectq.backends import CircuitDrawer, ResourceCounter, ClassicalSimulator
from projectq.meta import Loop, Compute, Uncompute, Control
from Matrix import Matrix, Matrix_2_4, Matrix_2_2, Matrix_2_8, Matrix_2_16, Matrix_2_40, Matrix_2_32, Matrix_2_64, Matrix_2_104

def Inversion(eng):

    n = 233

    f = eng.allocate_qureg(n)
    g = eng.allocate_qureg(n)
    h = eng.allocate_qureg(n)
    ancilla = eng.allocate_qureg(12180)

    if (resource_check != 1):
        Round_constant_XOR(eng, f, 0xf, n)
        Round_constant_XOR(eng, g, 0xf, n)

    f_inverse, ancilla, free_anc1, mul = Inverison_Itoh_Tsujii_based(eng, f, ancilla)  # f_inverse = f^-1

    count = 0
    fg, count, ancilla = recursive_karatsuba(eng, f_inverse, g, n, count, ancilla)  # fg = f^-1 * g
    fg = Reduction(eng, fg)

    CNOT233(eng, fg, h)

    if (resource_check != 1):
        print_state(eng, h, n)  # Should be 1.

def CNOT233(eng, a, b):
    for i in range(233):
        CNOT | (a[i], b[i])

def Inverison_Itoh_Tsujii_based(eng, a, ancilla) :

    n = 233

    a2 = eng.allocate_qureg(n)
    A_2_2 = eng.allocate_qureg(n)

    a2 = Squaring(eng, a, a2, n) # a^2

    count = 0
    A, count, ancilla = recursive_karatsuba(eng, a, a2, n, count, ancilla) # a * a^2
    A = Reduction(eng, A)

    A_2_2 = Squaring_2_2(eng, A, A_2_2, n)
    a2 = Squaring(eng, a, a2, n) # Reverse

    count = 0
    B, count, ancilla = recursive_karatsuba(eng, A_2_2, A, n, count, ancilla)  # A * A^4
    B = Reduction(eng, B)

    #
    A_2_2 = Squaring_2_2(eng, A, A_2_2, n) # reverse
    B_2_4 = Squaring_2_4(eng, B, a2, n) # recycle

    count = 0
    C, count, ancilla = recursive_karatsuba(eng, B, B_2_4, n, count, ancilla)  # B * B_2_4
    C = Reduction(eng, C)

    B_2_4 = Squaring_2_4(eng, B, a2, n)  # reverse
    C_2_8 = Squaring_2_8(eng, C, A_2_2, n) # recycle

    count = 0
    D, count, ancilla = recursive_karatsuba(eng, C, C_2_8, n, count, ancilla)  # B * B_2_4
    D = Reduction(eng, D)

    C_2_8 = Squaring_2_8(eng, C, A_2_2, n)  # reverse
    D_2_16 = Squaring_2_16(eng, D, a2, n) # recycle

    count = 0
    E, count, ancilla = recursive_karatsuba(eng, D, D_2_16, n, count, ancilla)  # B * B_2_4 
    E = Reduction(eng, E)

    D_2_16 = Squaring_2_16(eng, D, a2, n)  # reverse
    E_2_32 = Squaring_2_32(eng, E, A_2_2, n) # recycle

    count = 0
    F, count, ancilla = recursive_karatsuba(eng, E, E_2_32, n, count, ancilla)  # B * B_2_4
    F = Reduction(eng, F)

    E_2_32 = Squaring_2_32(eng, E, A_2_2, n)  # reverse
    F_2_64 = Squaring_2_64(eng, F, a2, n) # recycle
    E_2_8 = Squaring_2_8(eng, E, A_2_2, n)  # recycle

    count = 0
    G, count, ancilla = recursive_karatsuba(eng, F, F_2_64, n, count, ancilla)  # B * B_2_4
    G = Reduction(eng, G)

    F_2_64 = Squaring_2_64(eng, F, a2, n)  # reverse
    F_2_40 = Squaring_2_40(eng, F, a2, n)  # recycle


    count = 0
    mul1, count, ancilla = recursive_karatsuba(eng, C, E_2_8, n, count, ancilla)  # B * B_2_4
    mul1 = Reduction(eng, mul1)

    E_2_8 = Squaring_2_8(eng, E, A_2_2, n)  # reverse
    G_2_104 = Squaring_2_104(eng, G, A_2_2, n)  # recycle

    count = 0
    mul2, count, ancilla = recursive_karatsuba(eng, mul1, F_2_40, n, count, ancilla)  # B * B_2_4
    mul2 = Reduction(eng, mul2)

    F_2_40 = Squaring_2_40(eng, F, a2, n)  # reverse

    count = 0
    mul3, count, ancilla = recursive_karatsuba(eng, mul2, G_2_104, n, count, ancilla)  # B * B_2_4
    mul3 = Reduction(eng, mul3)

    Result = Squaring(eng, mul3, a2, n) # recycle

    return Result, ancilla, A_2_2, mul3

def copy(eng, a, b, n):
    for i in range(n):
        CNOT | (a[i], b[i])

def Reduction(eng, result): # x^16 + x^5 + x^3 + x + 1
    n = 233

    for i in range(n - 1):
        if (i < n):
            CNOT | (result[i + n], result[i])

    for i in range(n - 1):
        if (i + 74 < n):
            CNOT | (result[i + n], result[i + 74])

    Modular_small(eng, result[392:], result, 73)

    return result[0:2 * n - 1]

def Modular_small(eng, input, result, size):

    for i in range(size):
        CNOT | (input[i], result[0 + i])
    for i in range(size):
        CNOT | (input[i], result[74 + i])

def Squaring_2_8(eng, input, output, n):

    for i in range(n):
        for j in range(n):
            if (Matrix_2_8[n * j + ((n - 1) - j - i) % n] == 1):
                CNOT | (input[(i + j) % n], output[j])

    return output

def Squaring_2_4(eng, input, output, n):

    for i in range(n):
        for j in range(n):
            if (Matrix_2_4[n * j + ((n - 1) - j - i) % n] == 1):
                CNOT | (input[(i + j) % n], output[j])

    return output

def Squaring_2_16(eng, input, output, n):

    for i in range(n):
        for j in range(n):
            if (Matrix_2_16[n * j + ((n - 1) - j - i) % n] == 1):
                CNOT | (input[(i + j) % n], output[j])

    return output

def Squaring_2_64(eng, input, output, n):

    for i in range(n):
        for j in range(n):
            if (Matrix_2_64[n * j + ((n - 1) - j - i) % n] == 1):
                CNOT | (input[(i + j) % n], output[j])

    return output

def Squaring_2_104(eng, input, output, n):

    for i in range(n):
        for j in range(n):
            if (Matrix_2_104[n * j + ((n - 1) - j - i) % n] == 1):
                CNOT | (input[(i + j) % n], output[j])

    return output

def Squaring_2_32(eng, input, output, n):

    for i in range(n):
        for j in range(n):
            if (Matrix_2_32[n * j + ((n - 1) - j - i) % n] == 1):
                CNOT | (input[(i + j) % n], output[j])

    return output

def Squaring_2_40(eng, input, output, n):

    for i in range(n):
        for j in range(n):
            if (Matrix_2_40[n * j + ((n - 1) - j - i) % n] == 1):
                CNOT | (input[(i + j) % n], output[j])

    return output

def Squaring_2_2(eng, input, output, n):

    for i in range(n):
        for j in range(n):
            if (Matrix_2_2[n * j + ((n - 1) - j - i) % n] == 1):
                CNOT | (input[(i + j) % n], output[j])

    return output

def Squaring(eng, input, output, n):

    for i in range(n):
        for j in range(n):
            if (Matrix[n * j + ((n - 1) - j - i) % n] == 1):
                CNOT | (input[(i + j) % n], output[j])

    return output

def Squaring_plus_input(eng, input, output, n):

    Matrix_input = Matrix

    for i in range(n):
        Matrix_input[((n - 1) -i) + n * i] = Matrix_input[((n - 1) - i) + n * i] ^ 1

    for i in range(n):
        for j in range(n):
            if (Matrix_input[n * j + ((n - 1) - j - i) % n] == 1):
                CNOT | (input[(i + j) % n], output[j])

    return output

def recursive_karatsuba(eng, a, b, n, count, ancilla):

    if(n==1):
        c = eng.allocate_qubit()
        Toffoli_gate(eng, a, b, c)

        return c, count, ancilla

    c_len = 3**math.log(n, 2) 
    r_low = n//2    

    if(n%2!=0):
        r_low = r_low +1 

    r_a = []
    r_b = []

    # Provide rooms and prepare operands
    r_a = ancilla[count:count + r_low]

    count = count + r_low

    r_b = ancilla[count:count + r_low]

    count = count + r_low

    with Compute(eng):
        for i in range(r_low):
            CNOT | (a[i], r_a[i])
        for i in range(n//2):
            CNOT | (a[r_low + i], r_a[i])
        for i in range(r_low):
            CNOT | (b[i], r_b[i])
        for i in range(n//2):
            CNOT | (b[r_low + i], r_b[i])

    # upper-part setting
    if(r_low == 1):
        c = eng.allocate_qureg(3)
        Toffoli_gate(eng, a[0], b[0], c[0])
        Toffoli_gate(eng, a[1], b[1], c[2])
        CNOT | (c[0], c[1])
        CNOT | (c[2], c[1])
        Toffoli_gate(eng, r_a, r_b, c[1])

        Uncompute(eng)
        return c, count, ancilla

    c_a = []
    c_b = []
    c_r = []

    c_a, count, ancilla = recursive_karatsuba(eng, a[0:r_low], b[0:r_low], r_low, count, ancilla)
    c_b, count, ancilla = recursive_karatsuba(eng, a[r_low:n], b[r_low:n], n//2, count, ancilla)
    c_r, count, ancilla = recursive_karatsuba(eng, r_a[0:r_low], r_b[0:r_low], r_low, count, ancilla) 

    Uncompute(eng)

    result = []
    result = combine(eng, c_a, c_b, c_r, n)

    return result, count, ancilla

def combine(eng, a, b, r, n):
    if (n % 2 != 0):
        for i in range(n):
            CNOT | (a[i], r[i])
        for i in range(n - 2):
            CNOT | (b[i], r[i])

        for i in range(n // 2):
            CNOT | (a[n // 2 + 1 + i], r[i])
        for i in range(n // 2):
            CNOT | (b[i], r[n // 2 + 1 + i])

        out = []
        for i in range(n // 2 + 1):  
            out.append(a[i])
        for i in range(n):
            out.append(r[i])
        for i in range((2 * n - 1) - n // 2 - 1 - n):
            out.append(b[n // 2 + i])

        return out

    half_n = int(n/2) 
    for i in range(n-1):
        CNOT | (a[i], r[i])
        CNOT | (b[i], r[i])
    for i in range(half_n-1):
        CNOT | (a[half_n+i], r[i])
        CNOT | (b[i], r[half_n+i])

    result = []
    for i in range(half_n):
        result.append(a[i])
    for i in range(n-1):
        result.append(r[i])
    for i in range(half_n):
        result.append(b[half_n-1+i])

    return result

def room(eng, length):

    room = eng.allocate_qureg(length)

    return room

def copy(eng, a, b, length):
    for i in range(length):
        CNOT | (a[i], b[i])

def Toffoli_gate(eng, a, b, c):

    if(NCT):
        Toffoli | (a, b, c)
    else:
        if (resource_check):
            if(AND_check):
                ancilla = eng.allocate_qubit()
                H | c
                CNOT | (b, ancilla)
                CNOT | (c, a)
                CNOT | (c, b)
                CNOT | (a, ancilla)
                Tdag | a
                Tdag | b
                T | c
                T | ancilla
                CNOT | (a, ancilla)
                CNOT | (c, b)
                CNOT | (c, a)
                CNOT | (b, ancilla)
                H | c
                S | c

            else:
                Tdag | a
                Tdag | b
                H | c
                CNOT | (c, a)
                T | a
                CNOT | (b, c)
                CNOT | (b, a)
                T | c
                Tdag | a
                CNOT | (b, c)
                CNOT | (c, a)
                T | a
                Tdag | c
                CNOT | (b, a)
                H | c


def Round_constant_XOR(eng, k, rc, bit):
    for i in range(bit):
        if (rc >> i & 1):
            X | k[i]

def print_state(eng, b, n):
    All(Measure) | b
    print('Result : ', end='')
    for i in range(n):
        print(int(b[n-1-i]), end='')
    print('\n')

global resource_check
global AND_check
global NCT

NCT = 1
resource_check = 0
classic = ClassicalSimulator()
eng = MainEngine(classic)
Inversion(eng)
eng.flush()

resource_check = 1
NCT = 1
AND_check = 0 # If AND = 1, the number of qubits must be manually counted.
Resource = ResourceCounter()
eng = MainEngine(Resource)
Inversion(eng)
print('\n')
print(Resource)