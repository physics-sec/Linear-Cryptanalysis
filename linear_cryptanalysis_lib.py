#! /usr/bin/env python3
# -*- coding: utf-8 -*-

from math import fabs, ceil
import multiprocessing
import concurrent.futures

initialized = False

# i know this is ugly
def initialize(num_p_c_pairs, sbox_bits , num_sboxes, num_rounds, min_bias, max_blocks_to_bf, do_sbox_param, do_inv_sbox_param, do_pbox_param):
    global NUM_P_C_PAIRS, SBOX_BITS , NUM_SBOXES, NUM_ROUNDS, MIN_BIAS, MAX_BLOCKS_TO_BF, do_sbox, do_inv_sbox, do_pbox, initialized
    NUM_P_C_PAIRS = num_p_c_pairs
    SBOX_BITS = sbox_bits
    NUM_SBOXES = num_sboxes
    NUM_ROUNDS = num_rounds
    MIN_BIAS = min_bias
    MAX_BLOCKS_TO_BF = max_blocks_to_bf
    do_sbox = do_sbox_param
    do_inv_sbox = do_inv_sbox_param
    do_pbox = do_pbox_param
    initialized = True

def apply_mask(value, mask):
    #retrieve the parity of mask/value
    interValue = value & mask
    total = 0
    while interValue > 0:
        temp = interValue & 1
        interValue = interValue >> 1
        if temp == 1:
            total = total ^ 1
    return total

# create the bias table for the sbox
# (if there are multiple sboxs, then this won't work)
def create_bias_table():
    if not initialized: exit('initialize the library first!')

    tablesize = 1 << SBOX_BITS
    bias_table = []
   
    table = []
    for x in range(1, tablesize):
        for y in range(1, tablesize):
            matches = 0
            for num in range(tablesize):
                # calculate the parity of the number before going in the sbox
                in_mask  = apply_mask(num, x)
                # calculate the parity of the number after going out the sbox
                # do_sbox is supposed to perform the substitution, make sure is well defined!
                out_mask = apply_mask(do_sbox(num), y)
                # if the parity is the same in both cases, add a match
                if in_mask == out_mask:
                    matches += 1

            # calculate the bias
            bias = (matches / tablesize) - 1/2
            # if the bias is greater than 0, save the 'x' and 'y' combination
            if bias > 0:
                table.append( [x, y, bias] )

    # return the table of biases that are greater than 0
    return table

# from an sbox and the "output" of a bias y,
# calculate which sboxs will be reached and in which bits
def get_destination(num_sbox, y):
    # pass 'y' through the permutation
    offset = (NUM_SBOXES - (num_sbox-1) - 1) * SBOX_BITS
    Y = y << offset
    # do_pbox is supposed to transpose the state, make sure is well defined!
    permuted = do_pbox(Y)

    sboxes_reached = {}
    # sboxes go from 1 to NUM_SBOXES from left to right
    # bits go from 1 to SBOX_BITS from left to right
    for sbox in range(1, NUM_SBOXES + 1):
        for bit in range(SBOX_BITS):
            bits_offset = ((NUM_SBOXES - (sbox-1) - 1) * SBOX_BITS) + bit
            # if 'sbox' has a 1 in the position 'bit' then take note of that
            if permuted & (1 << bits_offset) != 0:
                if sbox not in sboxes_reached:
                    sboxes_reached[sbox] = []
                sboxes_reached[sbox].append(SBOX_BITS - bit)
    # return which sboxes where reached and in which bit
    return sboxes_reached

# convert a list of bits to an integer
def bits_to_num(inputbits):
    Y_input = 0
    for input_pos in inputbits:
        Y_input |= 1 << (SBOX_BITS - input_pos)
    return Y_input

# convert an integer to a list of it's bits
def num_to_bits(num):
    bits = []
    for index in range(SBOX_BITS):
        if (1 << index) & num > 0:
            bits.append( SBOX_BITS - index )
    return bits

# this function eliminates the linear approximations that have
# a bias below the MIN_BIAS threshold and then sorts the result
def sort_linear_aproximations(linear_aproximations):
    if not initialized: exit('initialize the library first!')

    sorted_linear_aproximations = []
    for linear_aproximation in linear_aproximations:
        biases = linear_aproximation['biases']

        # calculate the resulting bias following the Piling-Up Lemma
        resulting_bias = 1
        for x, y, bias in biases:
            resulting_bias *= bias
        resulting_bias *= 1 << (len(biases) - 1)

        # construct the element of the list resulting list
        x, y, _ = biases[0]
        _, num_sbox = linear_aproximation['start']
        entry = [resulting_bias, [num_sbox, num_to_bits(x)], linear_aproximation['state']]
        # keep the entry only if has a bias grater than MIN_BIAS
        if resulting_bias > MIN_BIAS:
            sorted_linear_aproximations.append( entry )

    # sort and return the result
    sorted_linear_aproximations = sorted(sorted_linear_aproximations, key=lambda elem: fabs(elem[0]), reverse=True)
    return sorted_linear_aproximations

# calculate all the possible linear aproximations given a bias table
def get_linear_aproximations(bias_table, current_states=None, depth=1):
    if not initialized: exit('initialize the library first!')

    # run for NUM_ROUNDS - 1 times
    if depth == NUM_ROUNDS:
        # delete elements that involve more than MAX_BLOCKS_TO_BF final sboxes
        current_states = [elem for elem in current_states if len(elem['state']) <= MAX_BLOCKS_TO_BF]
        # return the linear aproximations that reach to no more than MAX_BLOCKS_TO_BF sboxes
        return current_states

    # at the beginnig, only one sbox can be chosen
    # (this could be done differently)
    if depth == 1:
        # for each bias and each sbox, calculate which sboxes are reached (in the lower layer)
        # this will be the next step's new initial state
        current_states = []
        for x, y, bias in bias_table:

            for num_sbox in range(1, NUM_SBOXES + 1):

                sboxes_reached = get_destination(num_sbox, y)

                entry = {}
                entry['start']  = [depth, num_sbox]
                entry['biases'] = [[x, y, bias]]
                entry['state']  = sboxes_reached

                current_states.append( entry )
        # call the function recursevely with the new current state and new depth
        return get_linear_aproximations(bias_table, current_states, depth + 1)

    else:
        # for each set of possible states it will do the following:
        #   for each sbox that we last reached,
        #   it will calculate all possible moves according to the bias table.
        #   then it will calculate all possible the combinations of choices
        # this set of combinations, will be our next 'current_states'
        # lastly, it will call itself recursevely
        next_states = []
        for current_state in current_states:

            curr_pos = current_state['state']

            # calculate all possible moves from 'curr_sbox'
            possible_step_per_sbox = {}
            for curr_sbox in curr_pos:

                inputs  = curr_pos[curr_sbox]
                Y_input = bits_to_num(inputs)

                possible_steps = []

                # only use the biases which input matches the current sbox
                possible_biases = [ elem for elem in bias_table if elem[0] == Y_input ]
                for x, y, bias in possible_biases:

                    sboxes_reached = get_destination(curr_sbox, y)

                    step = {'to': sboxes_reached, 'path': [x, y, bias]}

                    possible_steps.append(step)

                possible_step_per_sbox[curr_sbox] = possible_steps


            # combine all the possible choises of each sbox in all possible ways
            # for example, if there are 2 sboxes and each has 4 possible moves
            # then calculate all 16 (4x4) possible combinations.

            possible_steps_combinations = []

            # initialize the possible_steps_combinations array
            # with all the first sbox's possible steps
            first_sbox = list(possible_step_per_sbox.keys())[0]
            for possible_step in possible_step_per_sbox[first_sbox]:
                step = possible_step
                step['from'] = [depth, first_sbox]
                possible_steps_combinations.append( [step] )

            # for each 'curr_sbox' that is not 'first_sbox'...
            for curr_sbox in possible_step_per_sbox:
                if curr_sbox == first_sbox:
                    continue

                # save the combinations calcualted so far
                combinations_so_far = possible_steps_combinations
                possible_steps_combinations = []

                # add each curr_sbox's possible step to the possible_steps_combinations array
                for possible_step in possible_step_per_sbox[curr_sbox]:

                    for step_taken in combinations_so_far:

                        possible_step['from'] = [depth, curr_sbox]
                        add_step = step_taken.copy()
                        add_step.append(possible_step)
                        possible_steps_combinations.append( add_step )

            # now, for each combination, check to which sboxes we reached and what are their inputs
            # this will be the next state
            for possible_step in possible_steps_combinations:

                # save the first sbox and the previous biases
                entry = {}
                entry['start'] = current_state['start']
                entry['biases'] = current_state['biases'].copy()
                entry['state'] = {}

                # add the new biases
                for elem in possible_step:
                    entry['biases'].append( elem['path'] )

                    # add the final sboxes and their inputs
                    for destination in elem['to']:
                        if destination not in entry['state']:
                            entry['state'][destination] = []

                        new_bits = elem['to'][destination]
                        entry['state'][destination] += new_bits

                # update the next_states
                next_states.append( entry )


        return get_linear_aproximations(bias_table, next_states, depth + 1)

def analize_cipher():
    if not initialized: exit('initialize the library first!')

    # analize the sbox and create the bias table
    table = create_bias_table()
    table_sorted = sorted(table, key=lambda elem: fabs(elem[2]), reverse=True)

    # take the best 1000 results (so that the following algorithm finishes quickly)
    # TODO: this could be done in a better way
    max_size = 1000
    table_len = len(table_sorted)
    if table_len > max_size:
        print('\n[*] reducing bias table size from {:d} to {:d}\n'.format(table_len, max_size))
        table_sorted = table_sorted[:max_size]

    # calculate all possible linear aproximations (with a bias > 0)
    linear_aproximations = get_linear_aproximations(table_sorted)
    # sort the list from the best approximations to the worst
    linear_aproximations_sorted = sort_linear_aproximations(linear_aproximations)
    # return the sorted list of approximations
    return linear_aproximations_sorted

def bit(num, n):
    # get the nth bit of num
    return (num >> (SBOX_BITS - n)) & 1

def get_xor(plaintext, ciphertext, key, linear_aproximation):
    bias, p_data, c_data = linear_aproximation

    # get the plaintext block
    p_block_num, p_bits = p_data
    pt = plaintext >> ((NUM_SBOXES - p_block_num) * SBOX_BITS)
    pt = pt & ((1 << SBOX_BITS) - 1)

    # calculate the plaintext's part of the xor
    xor_pt = 0
    for b in p_bits:
        xor_pt = xor_pt ^ bit(pt, b)

    # for each final sbox, get the according ciphertext block
    xor_u = 0
    i = len(c_data) - 1
    for c_block_num in c_data:
        c_bits = c_data[c_block_num]

        # get the ciphertext block
        ct = ciphertext >> ((NUM_SBOXES - c_block_num) * SBOX_BITS)
        ct = ct & ((1 << SBOX_BITS) - 1)

        # get the key block that corresponds with the sbox
        k = key >> (i * SBOX_BITS)
        k = k & ((1 << SBOX_BITS) - 1)

        # xor the key and the ciphertext to get v (the sbox output)
        v = ct ^ k

        # get the sbox input
        # do_inv_sbox is supposed to calculate the inverse of the substitution, make sure is well defined!
        u = do_inv_sbox(v)

        # calculate the input of the sbox's part of the xor
        for b in c_bits:
            xor_u = xor_u ^ bit(u, b)

        i -= 1

    # return the result of the full xor
    return xor_pt ^ xor_u

def get_biases_for_key_space(keystart, keyend, p_c_pairs, linear_aproximation):

    try:
        hits = [0] * (keyend - keystart)
    except OverflowError:
        exit('the amount of key bits to brute force is too large.')

    # get the result of the aproximation for each possible key
    for key in range(keystart, keyend):
        for plaintext, ciphertext in p_c_pairs:
            xor = get_xor(plaintext, ciphertext, key, linear_aproximation)
            if xor == 0:
                hits[keystart - key] += 1

    result = {'start': keystart, 'end': keyend, 'hits': hits}
    return result

def get_biases(p_c_pairs, linear_aproximation):
    if not initialized: exit('initialize the library first!')

    # calculate how many key bits must be brute forced
    key_bits = len(linear_aproximation[2]) * SBOX_BITS
    try:
        # get the key's maximum size
        key_max  = 1 << key_bits
    except MemoryError:
        exit('the amount of key bits to brute force is too large.')

    num_cores = multiprocessing.cpu_count()

    sub_key_space = key_max // num_cores

    bias_lists = []

    # run in num_cores threads
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_cores) as executor:

        future_list = []
        # divide the key space into num_cores parts
        for core in range(num_cores):
            start = sub_key_space * core
            end   = start + sub_key_space

            future = executor.submit(get_biases_for_key_space, start, end, p_c_pairs, linear_aproximation)
            future_list.append(future)

        # get the result for each thread
        for future in concurrent.futures.as_completed(future_list):
            bias_lists.append(future.result())

    try:
        hits = [0] * key_max
    except OverflowError:
        exit('the amount of key bits to brute force is too large.')

    # join all the results
    for result in bias_lists:
        start = result['start']
        end   = result['end']
        array = result['hits']
        for hit in range(start, end):
            hits[hit] = array[start - hit]

    # calculate the bias for each key
    bias = [ fabs(num_hits - float(NUM_P_C_PAIRS/2)) / float(NUM_P_C_PAIRS) for num_hits in hits ]

    return bias


if __name__ == "__main__":
    print('import this in from script')
