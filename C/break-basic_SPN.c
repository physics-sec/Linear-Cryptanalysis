#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include "break-basic_SPN.h"
#include "basic_SPN.h"
#include "linear_cryptanalysis_lib.h"

int pbox[];
int sbox[];
int sbox_inv[];

// modify accordingly
int do_sbox(int number)
{
    return sbox[number];
}

// modify accordingly
int do_inv_sbox(int number)
{
    return sbox_inv[number];
}

// modify accordingly
int do_pbox(int state)
{
    int newSstate = 0;
    for (int bitIdx = 0; bitIdx < 16; bitIdx++)
    {
        if (state & (1 << bitIdx))
        {
            newSstate |= (1 << pbox[bitIdx]);
        }
    }
    return newSstate;
}

int main()
{
    struct timeval start, end, seed;
    printf("analizing cipher...\n\n");
    gettimeofday(&start, NULL);
    struct state** linear_aproximations = analize_cipher();
    gettimeofday(&end, NULL);
    printf("done. took %lu us\n", (end.tv_sec - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec);

    if (linear_aproximations[0] == 0)
    {
        printf("no linear aproximations could be found!\n");
        exit(1);
    }

    printf("\nthe linear aproximation with the best bias will be used\n");
    struct state linear_aproximation = *linear_aproximations[0];
    printState(linear_aproximation);

    // key generation
    gettimeofday(&seed, NULL);
    srand(seed.tv_usec);
    unsigned char key[] = {rand(), rand(), rand(), rand(), rand(), rand(), rand(), rand(), rand(), rand()};

    // calculate the keybits that should be recovered
    int lastRoundKey = (key[8] << 8) | key[9];
    int keySearch = 0;
    for (int sbox = 0; sbox < NUM_SBOXES; sbox++)
    {
        int input = bits_to_num(linear_aproximation.position[sbox]);
        if (input == 0) continue;
        int keyblock;
        keyblock = lastRoundKey >> ((NUM_SBOXES - sbox - 1) * SBOX_BITS);
        keyblock = keyblock & ((1 << SBOX_BITS) - 1);
        keyblock = keyblock << (((NUM_SBOXES - sbox - 1) * SBOX_BITS));
        keySearch |= keyblock;
    }

    // generate plaintext/ciphertext pairs
    unsigned long ciphertexts[NUM_P_C_PAIRS];
    unsigned long plaintexts[NUM_P_C_PAIRS];
    for (int i = 0; i < NUM_P_C_PAIRS; i++)
    {
        ciphertexts[i] = encrypt(i, key);
        plaintexts[i]  = i;
    }

    printf("\nbreaking cipher...\n\n");
    // obtain the biases given the p/c pairs and the linear aproximation
    gettimeofday(&start, NULL);
    double* biases = get_biases(plaintexts, ciphertexts, linear_aproximation);
    gettimeofday(&end, NULL);
    printf("done. took %lu us\n", (end.tv_sec - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec);

    // get the key with the most hits
    double maxResult = 0;
    int maxIdx = 0;

    for (int i = 0; i < 4096; i ++)
    {
        if (biases[i] > maxResult)
        {
            maxResult = biases[i];
            maxIdx = i;
        }
    }

    if (maxIdx == keySearch)
    {
        printf("key found: %d\n", maxIdx);
    }
    else
    {
        printf("Fail!\n");
    }

    freeMem(linear_aproximations);
    free(biases);
    return 0;
}
