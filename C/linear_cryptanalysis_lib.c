#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <sys/sysinfo.h>
#include "linear_cryptanalysis_lib.h"

uint32_t ARR_SIZE = 3000;

// retrieve the parity of mask/value
uint8_t apply_mask(uint64_t value, uint64_t mask)
{
    uint64_t interValue = value & mask;
    uint8_t  total = 0;
    uint64_t temp;
    while (interValue > 0)
    {
        temp = interValue & 1;
        interValue = interValue >> 1;
        if (temp == 1)
        {
            total = total ^ 1;
        }
    }
    return total;
}

// create the bias table for the sbox
// (if there are multiple sboxs, then this won't work)
uint64_t create_bias_table(struct sbox_aprox table[], uint64_t tablesize){
    uint64_t x, y, num;
    uint64_t index = 0;
    for (x = 1; x < tablesize; x++)
    {
        for (y = 1; y < tablesize; y++)
        {
            uint64_t matches = 0;
            for (num = 0; num < tablesize; num++)
            {
                // calculate the parity of the number before going in the sbox
                uint64_t in_mask  = apply_mask(num, x);
                // calculate the parity of the number after going out the sbox
                // do_sbox is supposed to perform the substitution, make sure is well defined!
                uint64_t out_mask = apply_mask(do_sbox(num), y);
                // if the parity is the same in both cases, add a match
                if (in_mask == out_mask)
                {
                    matches += 1;
                }
            }
            // calculate the bias
            double bias = ((double)matches / (double)tablesize) - (double)0.5;
            // if the bias is greater than 0, save the 'x' and 'y' combination
            if (bias > 0)
            {
                struct sbox_aprox elem;
                elem.x = x;
                elem.y = y;
                elem.bias = bias;
                table[index++] = elem;
            }
        }
    }
    return index;
}

// from an sbox and the "output" of a bias y,
// calculate which sboxs will be reached and in which bits
void get_destination(uint8_t sboxes_reached[NUM_SBOXES][SBOX_BITS], uint64_t pos_sbox, uint64_t y)
{
    // pass 'y' through the permutation
    uint64_t offset = (NUM_SBOXES - pos_sbox - 1) * SBOX_BITS;
    uint64_t Y = y << offset;
    // do_pbox is supposed to transpose the state, make sure is well defined!
    uint64_t permuted = do_pbox(Y);

    // sboxes go from 1 to NUM_SBOXES from left to right
    // bits go from 1 to SBOX_BITS from left to right
    for (uint8_t sbox = 1; sbox <= NUM_SBOXES; sbox++)
    {
        for (uint16_t bit = 0; bit < SBOX_BITS; bit++)
        {
            uint64_t bits_offset = ((NUM_SBOXES - (sbox-1) - 1) * SBOX_BITS) + bit;
            // if 'sbox' has a 1 in the position 'bit' then take note of that
            if ((permuted & (1 << bits_offset)) != 0)
            {
                sboxes_reached[sbox-1][SBOX_BITS - bit-1] = 1;
            }
            else
            {
                sboxes_reached[sbox-1][SBOX_BITS - bit-1] = 0;
            }
        }
    }
}

// convert a list of bits to an integer
uint64_t bits_to_num(uint8_t inputbits[])
{
    uint64_t Y_input = 0;
    for (uint16_t i = 0; i < SBOX_BITS; i++)
    {
        if (inputbits[i] == 1)
        {
            Y_input |= 1 << (SBOX_BITS - i-1);    
        }
    }
    return Y_input;
}

// convert an integer to a list of it's bits
void num_to_bits(unsigned num, uint64_t bits[])
{
    for (uint16_t index = 0; index < SBOX_BITS; index++)
    {
        if (((1 << index) & num) != 0)
        {
            bits[SBOX_BITS - index-1] = 1;
        }
        else
        {
            bits[SBOX_BITS - index-1] = 0;
        }
    }
}

int cmpfunc(const void * a, const void * b)
{
    double rest = ((struct sbox_aprox*)a)->bias - ((struct sbox_aprox*)b)->bias;
    if (rest == 0)
    {
        return 0;
    }
    if (rest > 0)
    {
        return -1;
    }
    else
    {
        return 1;
    }
}

void sort_sbox_laprox(struct sbox_aprox linear_aproximations[], uint64_t tableSize)
{
    qsort(linear_aproximations, tableSize, sizeof(struct sbox_aprox), cmpfunc);
}

struct state** resize(struct state* states[])
{
    ARR_SIZE *= 2;
    states = realloc(states, ARR_SIZE * sizeof(struct state*));
    if (states == NULL)
    {
        printf("Error! memory not allocated.");
        exit(-1);
    }
    memset(states + (int)(ARR_SIZE/2), 0, (int)(ARR_SIZE/2));
    return states;
}

// calculate all the possible linear aproximations given a bias table
struct state** get_linear_aproximations(struct sbox_aprox bias_table[], uint64_t tableSize, struct state* current_states[], uint64_t depth)
{
    // run for NUM_ROUNDS - 1 times
    if (depth == NUM_ROUNDS) return current_states;

    // at the beginnig, only one sbox can be chosen
    // (this could be done differently)
    if (depth == 1)
    {
        current_states = calloc(ARR_SIZE, sizeof(struct state*));
        if (current_states == NULL)
        {
            printf("Error! memory not allocated.");
            exit(-1);
        }

        // for each bias and each sbox, calculate which sboxes are reached (in the lower layer)
        // this will be the next step's new initial state
        uint64_t index = 0;
        for (uint64_t i = 0; i < tableSize; i++)
        {
            struct sbox_aprox s_aprox = bias_table[i];
            // if the current bias is less than MIN_BIAS, discard it
            if (s_aprox.bias < MIN_BIAS) continue;

            for (uint8_t pos_sbox = 0; pos_sbox < NUM_SBOXES; pos_sbox++)
            {
                struct state* first_state = calloc(1, sizeof(struct state));
                if (first_state == NULL)
                {
                    printf("Error! memory not allocated.");
                    exit(-1);
                }
                first_state->first_sbox = pos_sbox;
                first_state->first_x = s_aprox.x;
                first_state->biasesMult = s_aprox.bias;
                first_state->cantBiases = 1;
                get_destination(first_state->position, pos_sbox, s_aprox.y);

                current_states[index++] = first_state;
                if (index == ARR_SIZE) current_states = resize(current_states);
            }
        }
        // call the function recursevely with the new current state and new depth
        return get_linear_aproximations(bias_table, tableSize, current_states, 2);
    }
    else
    {
        // for each set of possible states it will do the following:
        //   for each sbox that we last reached,
        //   it will calculate all possible moves according to the bias table.
        //   then it will calculate all possible the combinations of choices
        // this set of combinations, will be our next 'current_states'
        // lastly, it will call itself recursevely
        uint64_t len_states = 0;
        while (current_states[len_states++] != NULL);
        len_states--;

        struct state** next_states = calloc(ARR_SIZE, sizeof(struct state*));
        uint64_t next_state_pos = 0;
        if (next_states == NULL)
        {
            printf("Error! memory not allocated.");
            exit(-1);
        }

        // calculate all possible moves from each 'curr_state'
        for (uint64_t current_state_pos = 0; current_state_pos < len_states; current_state_pos++)
        {
            struct state curr_state = *current_states[current_state_pos];

            struct step possible_step_per_sbox[NUM_SBOXES][tableSize];
            uint64_t num_possible_step_per_sbox[NUM_SBOXES];
            uint64_t num_combinations = 1;
            uint64_t cant_start_sboxes = 0;
            uint8_t  start_sboxes[NUM_SBOXES];

            for (uint8_t sbox_pos = 0; sbox_pos < NUM_SBOXES; sbox_pos++)
            {
                uint8_t* inputs = curr_state.position[sbox_pos];
                uint64_t Y_input = bits_to_num(inputs); // is uint64_t always enough? uint64_t may be better?
                if (Y_input == 0)
                {
                    // this sbox has no inputs
                    num_possible_step_per_sbox[sbox_pos] = 0;
                    continue;
                }

                uint64_t count = 0;
                for (uint64_t j = 0; j < tableSize; j++)
                {
                    struct sbox_aprox lin_aprox = bias_table[j];
                    // only use linear approximations which input matches the current sbox
                    if (lin_aprox.x != Y_input) continue;

                    struct step possible_step;
                    possible_step.from = sbox_pos;
                    possible_step.path = lin_aprox;
                    get_destination(possible_step.to, sbox_pos, lin_aprox.y);

                    possible_step_per_sbox[sbox_pos][count++] = possible_step;
                }
                num_possible_step_per_sbox[sbox_pos] = count;
                if (count > 0)
                {
                    num_combinations *= count;
                    start_sboxes[cant_start_sboxes++] = sbox_pos;
                }                
            }

            if(cant_start_sboxes == 0)
            {
                printf("there are no possible paths! Too few linear aproximations?");
                exit(1);
            }

            // combine all the possible choises of each sbox in all possible ways
            // for example, if there are 2 sboxes and each has 4 possible moves
            // then calculate all 16 (4x4) possible combinations.

            struct step possible_steps_combinations[num_combinations][NUM_SBOXES];

            for (uint64_t comb_num = 0; comb_num < num_combinations; comb_num++)
            {
                possible_steps_combinations[comb_num][start_sboxes[0]] = possible_step_per_sbox[start_sboxes[0]][comb_num % num_possible_step_per_sbox[start_sboxes[0]]];
                for (uint64_t sbox = 1; sbox < cant_start_sboxes; sbox++)
                {
                    uint64_t real_sbox = start_sboxes[sbox];
                    uint64_t mod = 1;
                    for (uint64_t prev_sbox = 0; prev_sbox < sbox; prev_sbox++)
                    {
                        mod *= num_possible_step_per_sbox[start_sboxes[prev_sbox]];
                    }
                    uint64_t index = (comb_num / mod) % num_possible_step_per_sbox[real_sbox];
                    possible_steps_combinations[comb_num][real_sbox] = possible_step_per_sbox[real_sbox][index];
                }
            }

            // now, for each combination, check to which sboxes we reached and what are their inputs
            // this will be the next state

            for (uint64_t comb = 0; comb < num_combinations; comb++)
            {
                struct step* combination = possible_steps_combinations[comb];
                struct state* new_state = calloc(1, sizeof(struct state));
                if (new_state == NULL)
                {
                    printf("Error! memory not allocated.");
                    exit(-1);
                }

                new_state->first_sbox = curr_state.first_sbox;
                new_state->first_x = curr_state.first_x;
                new_state->biasesMult = curr_state.biasesMult;
                new_state->cantBiases = curr_state.cantBiases;

                uint64_t new_biases = 0;
                for (uint8_t sbox = 0; sbox < NUM_SBOXES; sbox++)
                {
                    if (num_possible_step_per_sbox[sbox] == 0) continue;

                    struct step step = combination[sbox];
                    new_state->biasesMult *= step.path.bias;
                    new_state->cantBiases += 1;

                    for (uint8_t sb = 0; sb < NUM_SBOXES; sb++)
                    {
                        for (uint16_t bit = 0; bit < SBOX_BITS; bit++)
                        {
                            if (step.to[sb][bit] == 1)
                            {
                                new_state->position[sb][bit] = 1;
                            }
                        }
                    }
                }

                // if the current bias is less than MIN_BIAS, discard the state
                double biasTotal = new_state->biasesMult * (1 << (new_state->cantBiases - 1));
                if (biasTotal < MIN_BIAS)
                {
                    free(new_state);
                    continue;
                }

                next_states[next_state_pos++] = new_state;
                if (next_state_pos == ARR_SIZE) next_states = resize(next_states);
            }
            free(current_states[current_state_pos]);
        }
        free(current_states);
        return get_linear_aproximations(bias_table, tableSize, next_states, depth+1);
    }
}

int cmpLinearAprox(const void * a, const void * b)
{   
    struct state s1 = *(*(struct state**)a);
    struct state s2 = *(*(struct state**)b);

    double biasTotal1 = s1.biasesMult * (1 << (s1.cantBiases - 1));
    double biasTotal2 = s2.biasesMult * (1 << (s2.cantBiases - 1));

    double rest = biasTotal1 - biasTotal2;
    if (rest == 0)
    {
        return 0;
    }
    if (rest > 0)
    {
        return -1;
    }
    else
    {
        return 1;
    }
}

void sort_linear_aproximations(struct state* linear_aproximations[])
{
    uint64_t size = 0;
    while (1)
    {
        if (linear_aproximations[size] == 0)
        {
            break;
        }
        if (size == ARR_SIZE)
        {
            break;
        }
        size++;
    }
    qsort(linear_aproximations, (int)size, sizeof(struct state*), cmpLinearAprox);
}

struct state** analize_cipher(void)
{
    uint64_t tablesize = 1 << SBOX_BITS;
    struct sbox_aprox table [tablesize*tablesize];
    uint64_t tableSize = create_bias_table(table, tablesize);

    sort_sbox_laprox(table, tableSize);
    uint16_t maxsize = 100;
    if (tableSize > maxsize)
    {
        printf("\n[*] reducing bias table size from %d to %d\n", tableSize, maxsize);
        tableSize = maxsize;
    }

    struct state** linear_aproximations =  get_linear_aproximations(table, tableSize, NULL, 1);
    sort_linear_aproximations(linear_aproximations);
    return linear_aproximations;
}

// get the bitth bit of num, from left to right
uint8_t getBit(uint64_t num, uint64_t bit)
{
    return (num >> (SBOX_BITS - bit - 1)) & 1;
}

// get the xor acoording to the plaintext, the ciphertext and linear aproximation
// uint64_t is not always enough, a char array should be used
uint64_t get_xor(uint64_t plaintext, uint64_t ciphertext, uint64_t key, uint64_t cantBlocks, struct state linear_aproximation)
{
    uint64_t pt, ct, k, v, u;
    // get the plaintext block
    pt = plaintext >> ((NUM_SBOXES - linear_aproximation.first_sbox - 1) * SBOX_BITS);
    pt = pt & ((1 << SBOX_BITS) - 1);
    uint64_t xor_pt = 0;
    for (uint16_t bit = 0; bit < SBOX_BITS; bit++)
    {
        uint16_t i = SBOX_BITS - bit - 1;
        if ((1 << i) & linear_aproximation.first_x)
        {
            xor_pt = xor_pt ^ getBit(pt, bit);
        }
    }

    uint64_t keyblock = cantBlocks - 1;
    uint64_t xor_u = 0;
    // for each final sbox, get the according ciphertext block
    for (uint8_t sbox = 0; sbox < NUM_SBOXES; sbox++)
    {
        uint64_t sbox_input = bits_to_num(linear_aproximation.position[sbox]);
        if (sbox_input == 0) continue;

        // get the ciphertext block
        ct = ciphertext >> ((NUM_SBOXES - sbox - 1) * SBOX_BITS);
        ct = ct & ((1 << SBOX_BITS) - 1);

        // get the key block that corresponds with the sbox
        k = key >> (keyblock * SBOX_BITS);
        k = k & ((1 << SBOX_BITS) - 1);
        keyblock--;

        // xor the key and the ciphertext to get v (the sbox output)
        v = ct ^ k;

        // get the sbox input
        // do_inv_sbox is supposed to calculate the inverse of the substitution, make sure is well defined!
        u = do_inv_sbox(v);

        // calculate the input of the sbox's part of the xor
        for (uint16_t bit = 0; bit < SBOX_BITS; bit++)
        {
            uint16_t i = SBOX_BITS - bit - 1;
            if ((1 << i) & sbox_input)
            {
                xor_u = xor_u ^ getBit(u, bit);
            }
        }
    }
    // return the xor between the plaintext and the ciphertext parts}
    return xor_pt ^ xor_u;
}

void* get_biases_for_key_space(void* args)
{
    // get the arguments struct
    struct threadParam params = *((struct threadParam*)args);

    double* hits = calloc(params.keyend - params.keystart, sizeof(double));
    if (hits == NULL)
    {
        printf("Error! memory not allocated.");
        exit(-1);
    }

    // for each possible key
    for (uint64_t key = params.keystart; key < params.keyend; key++)
    {
        for (uint64_t i = 0; i < NUM_P_C_PAIRS; i++)
        {
            // check if the linear aproximation checks out for each pair of plaintext/ciphertext
            uint64_t xor = get_xor(params.plaintexts[i], params.ciphertexts[i], key, params.cantBlocks, params.linear_aproximation);
            if (xor == 0)
            {
                hits[key - params.keystart]++;
            }
        }
    }

    // calculate the resulting bias following the Piling-Up Lemma
    for (uint64_t hit = params.keystart; hit < params.keyend; hit++)
    {
        hits[hit - params.keystart] = fabs(hits[hit - params.keystart] - (double)(NUM_P_C_PAIRS/(double)2)) / (double)(NUM_P_C_PAIRS);
    }

    // set the result struct
    struct partialResult* result = calloc(3, sizeof(struct partialResult));
    if (result == NULL)
    {
        printf("Error! memory not allocated.");
        exit(-1);
    }
    result->keystart = params.keystart;
    result->ketend   = params.keyend;
    result->hits     = hits;
    return (void*)result;
}

double* get_biases(uint64_t plaintexts[], uint64_t ciphertexts[], struct state linear_aproximation)
{
    // calculate how many key blocks must be brute forced
    uint64_t cantBlocks = 0;
    for (uint8_t sbox = 0; sbox < NUM_SBOXES; sbox++)
    {
        uint64_t num = bits_to_num(linear_aproximation.position[sbox]);
        if (num > 0) cantBlocks++;
    }

    // calculate how many key bits must be brute forced
    uint32_t key_bits = cantBlocks * SBOX_BITS;
    uint64_t key_max = 1 << key_bits;

    // run in num_cores threads
    uint8_t num_cores = get_nprocs_conf();
    num_cores *= 2;
    printf("using %d threads\n\n", num_cores);
    uint64_t sub_key_space = key_max / num_cores;
    pthread_t t_ids[num_cores];
    struct threadParam param[num_cores];
    struct partialResult* results[num_cores];

    for (uint8_t core = 0; core < num_cores; core++)
    {
        param[core].keystart = sub_key_space * core;
        param[core].keyend = param[core].keystart + sub_key_space;
        param[core].cantBlocks = cantBlocks;
        param[core].plaintexts = plaintexts;
        param[core].ciphertexts = ciphertexts;
        param[core].linear_aproximation = linear_aproximation;
        // each thread will handle a segment of the key space
        pthread_create(&t_ids[core], NULL, get_biases_for_key_space, (void*)&param[core]);
    }
    for (uint8_t core = 0; core < num_cores; core++)
    {
        pthread_join(t_ids[core], (void*)&results[core]);
    }

    double* hits = calloc(key_max, sizeof(double));
    if (hits == NULL)
    {
        printf("Error! memory not allocated.");
        exit(-1);
    }

    // get and combine the result for each thread
    for (uint8_t core = 0; core < num_cores; core++)
    {
        struct partialResult result = *results[core];
        uint64_t start = result.keystart;
        uint64_t end   = result.ketend;
        for (uint64_t i = 0; i < (end - start); i++)
        {
            hits[start + i] = result.hits[i];
        }
        free(result.hits);
        free(results[core]);
    }

    // return the array of biases
    return hits;
}

void printState(struct state state)
{
    printf("first_sbox:%d\n", state.first_sbox);
    printf("first_x:%d\n", state.first_x);

    double biasTotal = state.biasesMult * (1 << (state.cantBiases - 1));
    printf("bias:%f\n", biasTotal);

    for (uint8_t sbox = 0; sbox < NUM_SBOXES; sbox++)
    {
        uint64_t input = bits_to_num(state.position[sbox]);
        if (input > 0)
        {
            printf("sbox:%d -> %d\n", sbox, input);
        }
    }
    printf("\n");
}

void freeMem(struct state** linear_aproximations)
{
    for (uint64_t i = 0; i < ARR_SIZE; i++)
    {
        if (linear_aproximations[i] != 0)
        {
            free(linear_aproximations[i]);
        }
        else
        {
            break;
        }
    }
    free(linear_aproximations);
}
