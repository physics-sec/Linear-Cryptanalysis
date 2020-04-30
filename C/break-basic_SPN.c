#include <stdio.h>
#include <stdint.h>
#include <sys/time.h>
#include <stdlib.h>
#include "basic_SPN.h"
#include "linear_cryptanalysis_lib.h"

// modify accordingly
uint16_t do_sbox(uint16_t number)
{
    return sbox[number];
}

// modify accordingly
uint16_t do_inv_sbox(uint16_t number)
{
    return sbox_inv[number];
}

// modify accordingly
uint16_t do_pbox(uint16_t state)
{
    uint16_t newSstate = 0;
    for (uint8_t bitIdx = 0; bitIdx < 16; bitIdx++)
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
    printf("done. took %lu μs\n", (end.tv_sec - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec);

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
    uint8_t key[] = {rand(), rand(), rand(), rand(), rand(), rand(), rand(), rand(), rand(), rand()};

    // calculate the keybits that should be recovered
    uint16_t lastRoundKey = (key[8] << 8) | key[9];
    uint16_t keySearch = 0;
    uint8_t numKeyBlocks = 0;
    for (uint8_t sbox = 0; sbox < NUM_SBOXES; sbox++)
    {
        uint64_t input = bits_to_num(linear_aproximation.position[sbox]);
        if (input == 0) continue;
        numKeyBlocks++;
        uint16_t keyblock;
        keyblock = lastRoundKey >> ((NUM_SBOXES - sbox - 1) * SBOX_BITS);
        keyblock = keyblock & ((1 << SBOX_BITS) - 1);
        keyblock = keyblock << (((NUM_SBOXES - sbox - 1) * SBOX_BITS));
        keySearch |= keyblock;
    }

    // generate plaintext/ciphertext pairs
    uint64_t ciphertexts[NUM_P_C_PAIRS];
    uint64_t plaintexts[NUM_P_C_PAIRS];
    for (uint32_t i = 0; i < NUM_P_C_PAIRS; i++)
    {
        ciphertexts[i] = encrypt(i, key);
        plaintexts[i]  = i;
    }

    printf("\nbreaking cipher...\n\n");
    // obtain the biases given the p/c pairs and the linear aproximation
    gettimeofday(&start, NULL);
    double* biases = get_biases(plaintexts, ciphertexts, linear_aproximation);
    gettimeofday(&end, NULL);
    printf("done. took %lu μs\n", (end.tv_sec - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec);

    // get the key with the most hits
    double maxResult = 0;
    uint64_t maxIdx = 0;

    uint64_t maxKey = 1 << (numKeyBlocks * SBOX_BITS);
    for (int i = 0; i < maxKey; i ++)
    {
        if (biases[i] > maxResult)
        {
            maxResult = biases[i];
            maxIdx = i;
        }
    }

    if (maxIdx == keySearch)
    {
        printf("\nkey found: %d\n", maxIdx);
    }
    else
    {
        printf("\nFail!\n");
    }

    freeMem(linear_aproximations);
    free(biases);
    return 0;
}
