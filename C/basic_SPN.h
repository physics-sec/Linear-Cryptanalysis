#include <stdint.h>
#include <stdio.h>

#ifndef basic_SPN
#define basic_SPN

uint16_t _apply_sbox(uint16_t state, int* type_sbox);

uint16_t apply_sbox(uint16_t state);

uint16_t apply_sbox_inv(uint16_t state);

uint16_t encrypt(uint16_t pt, unsigned char k[]);

uint16_t decrypt(uint16_t ct, unsigned char k[]);

#endif