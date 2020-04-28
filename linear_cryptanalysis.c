#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "linear_cryptanalysis.h"

int ARR_SIZE = 10000;
int sbox[] = {0xE, 0x4, 0xD, 0x1, 0x2, 0xF, 0xB, 0x8, 0x3, 0xA, 0x6, 0xC, 0x5, 0x9, 0x0, 0x7};
int pbox[] = {0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15};

// modify accordingly
int do_sbox(int number)
{
    return sbox[number];
}

// modify accordingly
int do_pbox(int state)
{
    int state_temp = 0;
    for (int bitIdx = 0; bitIdx < 16; bitIdx++)
    {
        if (state & (1 << bitIdx))
        {
            state_temp |= (1 << pbox[bitIdx]);
        }
    }
    state = state_temp;
    return state;
}

// retrieve the parity of mask/value
int apply_mask(int value, int mask)
{
    int interValue = value & mask;
    int total = 0;
    int temp;
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
int create_bias_table(struct sbox_aprox table[], int tablesize){
    int x, y, num;
    int index = 0;
    for (x = 1; x < tablesize; x++)
    {
        for (y = 1; y < tablesize; y++)
        {
            int matches = 0;
            for (num = 0; num < tablesize; num++)
            {
                // calculate the parity of the number before going in the sbox
                int in_mask  = apply_mask(num, x);
                // calculate the parity of the number after going out the sbox
                // do_sbox is supposed to perform the substitution, make sure is well defined!
                int out_mask = apply_mask(do_sbox(num), y);
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

void printTable(struct sbox_aprox table[], int tableSize)
{
    printf("table:\n");
    for (int i = 0; i < tableSize; i++)
    {
        printf("x:%d y:%d bias:%f\n", table[i].x, table[i].y, table[i].bias);
    }
}

// from an sbox and the "output" of a bias y,
// calculate which sboxs will be reached and in which bits
void get_destination(int sboxes_reached[NUM_SBOXES][SBOX_BITS], int pos_sbox, int y)
{
    // pass 'y' through the permutation
    int offset = (NUM_SBOXES - pos_sbox - 1) * SBOX_BITS;
    int Y = y << offset;
    // do_pbox is supposed to transpose the state, make sure is well defined!
    int permuted = do_pbox(Y);

    //int sboxes_reached[NUM_SBOXES][SBOX_BITS];
    // sboxes go from 1 to NUM_SBOXES from left to right
    // bits go from 1 to SBOX_BITS from left to right
    for (int sbox = 1; sbox <= NUM_SBOXES; sbox++)
    {
        for (int bit = 0; bit < SBOX_BITS; bit++)
        {
            int bits_offset = ((NUM_SBOXES - (sbox-1) - 1) * SBOX_BITS) + bit;
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
int bits_to_num(int inputbits[])
{
    int Y_input = 0;
    for (int i = 0; i < SBOX_BITS; i++)
    {
        if (inputbits[i] == 1)
        {
            Y_input |= 1 << (SBOX_BITS - i-1);    
        }
    }
    return Y_input;
}

// convert an integer to a list of it's bits
void num_to_bits(unsigned num, int bits[])
{
    for (int index = 0; index < SBOX_BITS; index++)
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

void sort_sbox_laprox(struct sbox_aprox linear_aproximations[], int tableSize)
{
    qsort(linear_aproximations, tableSize, sizeof(struct sbox_aprox), cmpfunc);
}

struct state** get_linear_aproximations(struct sbox_aprox bias_table[], int tableSize, struct state* current_states[], int depth)
{
    if (depth == NUM_ROUNDS) return current_states;

    if (depth == 1)
    {
        current_states = calloc(ARR_SIZE, sizeof(struct state*));
        if (current_states == NULL)
        {
            printf("Error! memory not allocated.");
            exit(-1);
        }
        int index = 0;
        for (int i = 0; i < tableSize; i++)
        {
            struct sbox_aprox s_aprox = bias_table[i];
            for (int pos_sbox = 0; pos_sbox < NUM_SBOXES; pos_sbox++)
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
                if (index == ARR_SIZE)
                {
                    current_states = realloc(current_states, 2 * ARR_SIZE * sizeof(struct state*));
                    ARR_SIZE *= 2;
                    if (current_states == NULL)
                    {
                        printf("Error! memory not allocated.");
                        exit(-1);
                    }
                }
            }
        }
        return get_linear_aproximations(bias_table, tableSize, current_states, 2);
    }
    else
    {
        int len_states = 0;
        while (current_states[len_states++] != NULL);
        len_states--;

        struct state** next_states = calloc(ARR_SIZE, sizeof(struct state*));
        if (next_states == NULL)
        {
            printf("Error! memory not allocated.");
            exit(-1);
        }
        int next_state_pos = 0;

        for (int current_state_pos = 0; current_state_pos < len_states; current_state_pos++)
        {
            struct state curr_state = *current_states[current_state_pos];

            struct step possible_step_per_sbox[NUM_SBOXES][tableSize];
            int num_possible_step_per_sbox[NUM_SBOXES];
            int num_combinations = 1;
            int cant_start_sboxes = 0;

            for (int sbox_pos = 0; sbox_pos < NUM_SBOXES; sbox_pos++)
            {
                int* inputs = curr_state.position[sbox_pos];
                int Y_input = bits_to_num(inputs); // is int always enough? unsigned long may be better?
                if (Y_input == 0)
                {
                    // this sbox has no inputs
                    num_possible_step_per_sbox[sbox_pos] = 0;
                    continue;
                }

                int count = 0;
                for (int j = 0; j < tableSize; j++)
                {
                    struct sbox_aprox lin_aprox = bias_table[j];
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
                    cant_start_sboxes++;
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

            int start_sboxes[cant_start_sboxes];
            int index = 0;
            for (int i = 0; i < NUM_SBOXES; i++)
            {
                if (num_possible_step_per_sbox[i] > 0)
                {
                    start_sboxes[index++] = i;
                }
            }

            struct step possible_steps_combinations[num_combinations][NUM_SBOXES];

            for (int comb_num = 0; comb_num < num_combinations; comb_num++)
            {
                possible_steps_combinations[comb_num][start_sboxes[0]] = possible_step_per_sbox[start_sboxes[0]][comb_num % num_possible_step_per_sbox[start_sboxes[0]]];
                for (int sbox = 1; sbox < cant_start_sboxes; sbox++)
                {
                    int real_sbox = start_sboxes[sbox];
                    int mod = 1;
                    for (int prev_sbox = 0; prev_sbox < sbox; prev_sbox++)
                    {
                        mod *= num_possible_step_per_sbox[start_sboxes[prev_sbox]];
                    }
                    int index = (comb_num / mod) % num_possible_step_per_sbox[real_sbox];
                    possible_steps_combinations[comb_num][real_sbox] = possible_step_per_sbox[real_sbox][index];
                }
            }

            // now, for each combination, check to which sboxes we reached and what are their inputs
            // this will be the next state

            for (int comb = 0; comb < num_combinations; comb++)
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

                int new_biases = 0;
                for (int sbox = 0; sbox < NUM_SBOXES; sbox++)
                {
                    if (num_possible_step_per_sbox[sbox] == 0) continue;

                    struct step step = combination[sbox];
                    new_state->biasesMult *= step.path.bias;
                    new_state->cantBiases += 1;

                    for (int sb = 0; sb < NUM_SBOXES; sb++)
                    {
                        for (int bit = 0; bit < SBOX_BITS; bit++)
                        {
                            if (step.to[sb][bit] == 1)
                            {
                                new_state->position[sb][bit] = 1;
                            }
                        }
                    }
                }

                double biasTotal = new_state->biasesMult * (1 << (new_state->cantBiases - 1));
                if (biasTotal >= MIN_BIAS)
                {
                    next_states[next_state_pos++] = new_state;
                }

                //next_states[next_state_pos++] = new_state;

                if (next_state_pos == ARR_SIZE)
                {
                    ARR_SIZE *= 2;
                    next_states = realloc(next_states, ARR_SIZE * sizeof(struct state*));
                    if (next_states == NULL)
                    {
                        printf("Error! memory not allocated.");
                        exit(-1);
                    }
                    memset(next_states + (int)(ARR_SIZE/2), 0, (int)(ARR_SIZE/2));
                }
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
    int size = 0;
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
    qsort(linear_aproximations, size, sizeof(struct state*), cmpLinearAprox);
}

struct state** analize_cipher(void)
{
    int tablesize = 1 << SBOX_BITS;
    struct sbox_aprox table [tablesize*tablesize];
    int tableSize = create_bias_table(table, tablesize);

    sort_sbox_laprox(table, tableSize);
    int maxsize = 100;
    if (tableSize > maxsize)
    {
        printf("\n[*] reducing bias table size from %d to %d\n", tableSize, maxsize);
        tableSize = maxsize;
    }

    struct state** linear_aproximations =  get_linear_aproximations(table, tableSize, NULL, 1);
    sort_linear_aproximations(linear_aproximations);
    return linear_aproximations;
}

// get the nth bit of num
int bit(int num, int n)
{
    return (num >> (SBOX_BITS - n)) & 1;
}

void get_xor(void)
{
    return;
}

void get_biases_for_key_space(void)
{
    return;
}

void get_biases(void)
{
    return;
}

void printState(struct state state)
{
    printf("first_sbox:%d\n", state.first_sbox);
    printf("first_x:%d\n", state.first_x);

    double biasTotal = state.biasesMult * (1 << (state.cantBiases - 1));
    printf("bias:%f\n", biasTotal);

    for (int sbox = 0; sbox < NUM_SBOXES; sbox++)
    {
        int input = bits_to_num(state.position[sbox]);
        if (input > 0)
        {
            printf("sbox:%d -> %d\n", sbox, input);
        }
    }
    printf("\n");
}

int main()
{
    struct state** linear_aproximations = analize_cipher();
    for (int i = 0; i < 30; i++)
    {
        printState(*linear_aproximations[i]);
    }
    return 0;
}
