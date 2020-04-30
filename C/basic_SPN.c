#include <stdint.h>
#include <stdio.h>
#include "basic_SPN.h"

uint8_t pbox[]     = {0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15};
uint8_t sbox[]     = {0xE, 0x4, 0xD, 0x1, 0x2, 0xF, 0xB, 0x8, 0x3, 0xA, 0x6, 0xC, 0x5, 0x9, 0x0, 0x7};
uint8_t sbox_inv[] = {0xE, 3, 4, 8, 1, 0xC, 0xA, 0xF, 7, 0xD, 9, 6, 0xB, 2,0 , 5};

uint16_t _apply_sbox(uint16_t state, uint8_t* type_sbox)
{
    uint16_t newState, temp;
    newState = 0;    
    for (uint8_t idx = 0; idx < 4; idx++)
    {
        temp = (state & (0xf << (idx * 4))) >> (idx * 4);
        temp = type_sbox[temp];
        newState |= temp << (idx * 4);
    }
    return newState;
}

uint16_t apply_sbox(uint16_t state)
{
    return _apply_sbox(state, sbox);
}

uint16_t apply_sbox_inv(uint16_t state)
{
    return _apply_sbox(state, sbox_inv);
}

uint16_t encrypt(uint16_t pt, uint8_t k[])
{
    uint16_t state = pt;
    uint16_t subkey;

    for(uint8_t round = 0; round < 3; round++)
    {
        subkey = (k[round*2] << 8) | k[1 + round*2];

        // XOR state with round key (3, subkeys 1,..,4)
        state = state ^ subkey;

        // Break state into nibbles, perform sbox on each nibble, write to state (1)
        state = apply_sbox(state);

        // Permute the state bitwise (2)
        uint16_t state_temp = 0;
        for (uint8_t bitIdx = 0; bitIdx < 16; bitIdx++)
        {
            if(state & (1 << bitIdx))
            {
                state_temp |= (1 << pbox[bitIdx]);
            }
        }
        state = state_temp;
    }
    // Final round of SPN cipher (k4, sbox, s5)
    subkey = (k[6] << 8) | k[7];

    state = state ^ subkey;

    state = apply_sbox(state);

    subkey = (k[8] << 8) | k[9];

    // Final subkey (key 5) mixing
    state = state ^ subkey;

    return state;
}

uint16_t decrypt(uint16_t ct, uint8_t k[])
{
    uint16_t state = ct;
    uint16_t subkey;

    subkey = (k[8] << 8) | k[9];
    // Undo final round key
    state = state ^ subkey;

    // Apply inverse s-box
    state = apply_sbox_inv(state);

    // Undo first 3 rounds of simple SPN cipher
    for (uint8_t round = 3; round > 0; round--)
    {
        subkey = (k[round*2] << 8) | k[1 + round*2];

        // XOR state with round key (3, subkeys 4,..,0)
        state = state ^ subkey;

        // Un-permute the state bitwise (2)
        uint16_t state_temp = 0;
        for (uint8_t bitIdx = 0; bitIdx < 16; bitIdx++)
        {
            if(state & (1 << bitIdx))
            {
                state_temp |= (1 << pbox[bitIdx]);
            }
        }
        state = state_temp;

        // Apply inverse s-box
        state = apply_sbox_inv(state);
    }
    subkey = (k[0] << 8) | k[1];

    state = state ^ subkey;

    return state;
}
