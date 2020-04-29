#include <stdio.h>
#include <time.h>
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
    printf("analizing cipher...\n");
    struct state** linear_aproximations = analize_cipher();
    if (linear_aproximations[0] == 0)
    {
        printf("no linear aproximations could be found!\n");
        exit(1);
    }
    printf("\nbest linear aproximations:\n");
    for (int i = 0; i < 10; i++)
    {
        if (linear_aproximations[i] != 0)
        {
            printState(*linear_aproximations[i]);
        }
    }

    printf("\nthe linear aproximation with the best bias will be used\n");
    struct state linear_aproximation = *linear_aproximations[0];
    printState(linear_aproximation);

    // key generation
    srand(time(NULL));
    unsigned char key[] = {rand(), rand(), rand(), rand(), rand(), rand(), rand(), rand(), rand(), rand()};

    int lastRoundKey = (key[8] << 8) | key[9];
    int keySearch = lastRoundKey & 0x0fff;

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
    double* biases = get_biases(plaintexts, ciphertexts, linear_aproximation);

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
