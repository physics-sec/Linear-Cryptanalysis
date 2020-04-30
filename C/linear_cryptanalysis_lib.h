
#ifndef LINEARCRYPTLIB
#define LINEARCRYPTLIB

#define NUM_P_C_PAIRS 10000
#define SBOX_BITS 4
#define NUM_SBOXES 4
#define NUM_ROUNDS 4
#define MIN_BIAS 0.008
#define MAX_BLOCKS_TO_BF 3

struct sbox_aprox {
  uint8_t x, y;
  double bias;
};

struct state {
  uint8_t first_sbox;
  uint8_t first_x;
  double biasesMult;
  uint8_t cantBiases;
  uint8_t position[NUM_SBOXES][SBOX_BITS];
};

struct step {
  uint8_t from;
  uint8_t to[NUM_SBOXES][SBOX_BITS];
  struct sbox_aprox path;
};

struct partialResult {
  uint64_t keystart;
  uint64_t ketend;
  double* hits;;
};

struct threadParam {
  uint64_t keystart;
  uint64_t keyend;
  uint64_t cantBlocks;
  uint64_t* plaintexts;
  uint64_t* ciphertexts;
  struct state linear_aproximation;
};

struct state** analize_cipher(void);
void freeMem(struct state** linear_aproximations);
void printState(struct state state);
double* get_biases(uint64_t plaintexts[], uint64_t ciphertexts[], struct state linear_aproximation);
uint64_t bits_to_num(uint8_t inputbits[]);

#endif