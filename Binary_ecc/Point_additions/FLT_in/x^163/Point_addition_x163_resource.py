import math
from math import log2

from projectq import MainEngine
from projectq.ops import H, CNOT, Measure, Toffoli, X, All, T, Tdag, S, Swap
from projectq.backends import CircuitDrawer, ResourceCounter, ClassicalSimulator
from projectq.meta import Loop, Compute, Uncompute, Control
from Matrix import Matrix, Matrix_2_4, Matrix_2_2, Matrix_2_8, Matrix_2_16, Matrix_2_34, Matrix_2_32, Matrix_2_64

def Point_addition(eng):

    a = 3 #0b11
    b = 2
    n = 163

    q = eng.allocate_qubit()
    x1 = eng.allocate_qureg(n)
    y1 = eng.allocate_qureg(n)
    lambda_0 = eng.allocate_qureg(n)
    ancilla = eng.allocate_qureg(8448)

    x2 = 0xff00
    y2 = 0x00ff

    if (resource_check != 1):
        print('## Point Addition Result ##')
        X | q  # Comment/Uncomment

    if (resource_check != 1):
        Round_constant_XOR(eng, x1, 0xffff, n)
        Round_constant_XOR(eng, y1, 0xffff, n)

    dummy0 = eng.allocate_qureg(n)
    dummy1 = eng.allocate_qureg(n)

    Round_constant_XOR(eng, x1, x2, n)

    i_count = 0
    for i in range(n):
        if (y2 >> i & 1):
            i_count = i_count + 1

    with Compute(eng):
        copylist = copy_parallel(eng, q, ancilla, i_count)

    i_count = 0
    for i in range(n):
        if (y2 >> i & 1):
            with Control(eng, copylist[i_count]):
                X | y1[i]
            i_count = i_count + 1

    Uncompute(eng)

    anc0, ancilla, free_anc1, mul = Inverison_Itoh_Tsujii_based(eng, x1, ancilla)  # need reverse

    count = 0
    anc1, count, ancilla = recursive_karatsuba(eng, anc0, y1, n, count, ancilla)  # need reverse
    anc1 = Reduction(eng, anc1)
    CNOT163(eng, anc1, lambda_0)

    # reverse
    count = 0
    dummy_result, count, ancilla = recursive_karatsuba(eng, dummy0, dummy1, n, count, ancilla)
    dummy_result = Reduction(eng, dummy_result)
    dummy0, ancilla, free_anc1, mul = Inverison_Itoh_Tsujii_based(eng, dummy0, ancilla)

    count = 0
    anc1, count, ancilla = recursive_karatsuba(eng, x1, lambda_0, n, count, ancilla)  # need reverse
    anc1 = Reduction(eng, anc1)
    CNOT163(eng, anc1, y1)

    # reverse
    count = 0
    dummy_result, count, ancilla = recursive_karatsuba(eng, dummy0, dummy1, n, count, ancilla)
    dummy_result = Reduction(eng, dummy_result)

    with Compute(eng):
        Squaring_plus_input(eng, lambda_0, y1, n)
        Round_constant_XOR(eng, y1, a, n)
        Round_constant_XOR(eng, y1, x2, n)

    with Compute(eng):
        copylist1 = copy_parallel(eng, q, ancilla, n)

    for i in range(n):
        Toffoli_gate(eng, copylist1[i], y1[i], x1[i])

    Uncompute(eng)
    Uncompute(eng)

    count = 0
    anc2, count, ancilla = recursive_karatsuba(eng, x1, lambda_0, n, count, ancilla)  # need reverse
    anc2 = Reduction(eng, anc2)

    CNOT163(eng, anc2, y1)

    # reverse
    count = 0
    dummy_result, count, ancilla = recursive_karatsuba(eng, dummy0, dummy1, n, count, ancilla)
    dummy_result = Reduction(eng, dummy_result)

    anc0, ancilla, free_anc1, mul = Inverison_Itoh_Tsujii_based(eng, x1, ancilla)  # need reverse

    count = 0
    anc1, count, ancilla = recursive_karatsuba(eng, y1, anc0, n, count, ancilla)  # need reverse
    anc1 = Reduction(eng, anc1)

    CNOT163(eng, anc1, lambda_0)

    # reverse
    count = 0
    dummy_result, count, ancilla = recursive_karatsuba(eng, dummy0, dummy1, n, count, ancilla)
    dummy_result = Reduction(eng, dummy_result)
    dummy0, ancilla, free_anc1, mul = Inverison_Itoh_Tsujii_based(eng, dummy0, ancilla)

    i_count = 0
    for i in range(n):
        if (y2 >> i & 1):
            i_count = i_count + 1

    with Compute(eng):
        copylist2 = copy_parallel(eng, q, ancilla, n)

    i_count = 0
    for i in range(n):
        if (y2 >> i & 1):
            with Control(eng, copylist2[i_count]):
                X | y1[i]
            i_count = i_count + 1

    Round_constant_XOR(eng, x1, x2, n)

    for i in range(n):
        Toffoli_gate(eng, copylist2[i], x1[i], y1[i])

    Uncompute(eng)

def copy_parallel(eng, value, ancillas, n):

    divide = int(log2(n))
    last = n - 2 ** divide

    copy_list = []
    copy_list.append(value)
    for i in range(n - 1):
        copy_list.append(ancillas[i])

    for i in range(divide):
        for j in range(2 ** i):
            if (i == 0):
                CNOT | (copy_list[0], copy_list[1])
            else:
                CNOT | (copy_list[j], copy_list[2 ** i + j])

    for i in range(last):
        CNOT | (copy_list[i], copy_list[2 ** divide + i])

    return copy_list

def Copy(eng, a, b, n):
    for i in range(n):
        CNOT | (a, b[i])

def CSWAP(eng, a, b, c):
    CNOT | (c, b)
    Toffoli_gate(eng, a, b, c)
    CNOT | (c, b)

def CNOT163(eng, a, b):
    for i in range(163):
        CNOT | (a[i], b[i])

def Inverison_Itoh_Tsujii_based(eng, a, ancilla) :

    n = 163
    ancilla = eng.allocate_qureg(8448)

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
    E_2_2 = Squaring_2_2(eng, E, A_2_2, n)  # recycle
    F_2_64 = Squaring_2_64(eng, F, a2, n) # recycle

    count = 0
    G, count, ancilla = recursive_karatsuba(eng, F, F_2_64, n, count, ancilla)  # B * B_2_4
    G = Reduction(eng, G)

    F_2_64 = Squaring_2_64(eng, F, a2, n)  # reverse
    G_2_34 = Squaring_2_34(eng, G, a2, n) # recycle

    count = 0
    mul1, count, ancilla = recursive_karatsuba(eng, A, E_2_2, n, count, ancilla)  # B * B_2_4
    mul1 = Reduction(eng, mul1)

    E_2_2 = Squaring_2_2(eng, E, A_2_2, n)  # reverse

    count = 0
    mul2, count, ancilla = recursive_karatsuba(eng, mul1, G_2_34, n, count, ancilla)  # B * B_2_4
    mul2 = Reduction(eng, mul2)

    Result = Squaring(eng, mul2, A_2_2, n)

    return Result, ancilla, a2, mul2

def copy(eng, a, b, n):
    for i in range(n):
        CNOT | (a[i], b[i])

def Reduction(eng, result): # x^163 + x^7 + x^6 + x^3 + 1
    n = 163

    for i in range(n-1):
        if(i<163):
            CNOT | (result[i+163], result[i])

    for i in range(n - 1):
        if (i+3 < 163):
            CNOT | (result[i + 163], result[i+3])

    for i in range(n - 1):
        if (i+6 < 163):
            CNOT | (result[i + 163], result[i+6])

    for i in range(n - 1):
        if (i+7 < 163):
            CNOT | (result[i + 163], result[i+7])

    Modular_small(eng, result[323:], result, 2)
    Modular_small(eng, result[320:], result, 5)
    Modular_small(eng, result[319:], result, 6)

    return result[0:2 * n - 1]

def Modular_small(eng, input, result, size):

    for i in range(size):
        CNOT | (input[i], result[0 + i])
    for i in range(size):
        CNOT | (input[i], result[3 + i])
    for i in range(size):
        CNOT | (input[i], result[6 + i])
    for i in range(size):
        CNOT | (input[i], result[7 + i])

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

def Squaring_2_32(eng, input, output, n):

    for i in range(n):
        for j in range(n):
            if (Matrix_2_32[n * j + ((n - 1) - j - i) % n] == 1):
                CNOT | (input[(i + j) % n], output[j])

    return output

def Squaring_2_34(eng, input, output, n):

    for i in range(n):
        for j in range(n):
            if (Matrix_2_34[n * j + ((n - 1) - j - i) % n] == 1):
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
        Matrix_input[((n - 1) - i) + n * i] = Matrix_input[((n - 1) - i) + n * i] ^ 1

    for i in range(n):
        for j in range(n):
            if (Matrix_input[n * j + ((n - 1) - j - i) % n] == 1):
                CNOT | (input[(i + j) % n], output[j])

    for i in range(n):
        Matrix_input[((n - 1) - i) + n * i] = Matrix_input[((n - 1) - i) + n * i] ^ 1

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

resource_check = 1
NCT = 0
AND_check = 0 # If AND = 1, the number of qubits must be manually counted.
Resource = ResourceCounter()
eng = MainEngine(Resource)
Point_addition(eng)
print('\n')
print(Resource)