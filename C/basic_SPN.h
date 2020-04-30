#ifndef basic_SPN
#define basic_SPN

uint8_t pbox[];
uint8_t sbox[];
uint8_t sbox_inv[];

uint16_t _apply_sbox(uint16_t state, uint8_t* type_sbox);

uint16_t apply_sbox(uint16_t state);

uint16_t apply_sbox_inv(uint16_t state);

uint16_t encrypt(uint16_t pt, uint8_t k[]);

uint16_t decrypt(uint16_t ct, uint8_t k[]);

#endif