#! /usr/bin/env python3
# -*- coding: utf-8 -*-

from basic_SPN import *
import linear_cryptanalysis_lib as lc_lib
from math import fabs, ceil

# modify accordingly
def do_sbox(number):
    return sbox[number]

# modify accordingly
def do_inv_sbox(number):
    return sbox_inv[number]

# modify accordingly
def do_pbox(state):
    state_temp = 0
    for bitIdx in range(0,16):
        if(state & (1 << bitIdx)):
            state_temp |= (1 << pbox[bitIdx])
    state = state_temp
    return state

def main():

    NUM_P_C_PAIRS = 5000
    SBOX_BITS  = 4
    NUM_SBOXES = 4
    NUM_ROUNDS = 4
    MIN_BIAS = 0.008
    MAX_BLOCKS_TO_BF = 3

    lc_lib.initialize(NUM_P_C_PAIRS,
                      SBOX_BITS,
                      NUM_SBOXES,
                      NUM_ROUNDS,
                      MIN_BIAS,
                      MAX_BLOCKS_TO_BF,
                      do_sbox,
                      do_inv_sbox,
                      do_pbox)


    print('analizing cipher...')
    # there is no need to do this each time
    linear_aproximations = lc_lib.analize_cipher()

    print('\nbest 10 linear aproximations:')
    # just for demonstration
    for i in range(10):
        print(linear_aproximations[i])

    print('\nthe linear aproximation with the best bias will be used:')
    # you may choose anyone you like
    linear_aproximation = linear_aproximations[0]

    # just for demonstration
    end_sboxs = ', '.join(list(map(str, linear_aproximation[2])))
    print('ε: {:f}\nstart: sbox n°{:d}\nend sboxes:{}'.format(linear_aproximation[0], linear_aproximation[1][0], end_sboxs))
    Nl = ceil( pow(pow(linear_aproximation[0], -1), 2))
    print('\nneeded plaintext/ciphertext pairs (according to Matsui\'s paper):\nNl ≈ 1/(ε^2): {:d}'.format(ceil(Nl)))

    # this will be different with another cipher
    key = keyGeneration()
    # k is the last round key
    k = key[-3:]
    k = int(k, 16)

    # the 'encrypt' function might be different for you
    p_c_pairs = []
    for pt in range(NUM_P_C_PAIRS):
        p_c_pairs.append( [pt, encrypt(pt, key)] )

    print('\nbreaking cipher...\n')
    # obtain the biases given the p/c pairs and the linear aproximation
    biases = lc_lib.get_biases(p_c_pairs, linear_aproximation)

    # get the key with the highest bias
    maxResult, maxIdx = 0, 0
    for rIdx, result in enumerate(biases):
        if result > maxResult:
            maxResult = result
            maxIdx    = rIdx

    # maxIdx won't be equal to k if the final sboxes aren't consecutive!
    # in this example, they are
    if maxIdx == k:
        print('Success!')
        print('obtained key bits: {:d}'.format(maxIdx))
    else:
        print('Failure')


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        pass
