
#define NUM_P_C_PAIRS 10000
#define SBOX_BITS 4
#define NUM_SBOXES 4
#define NUM_ROUNDS 4
#define MIN_BIAS 0.008
#define MAX_BLOCKS_TO_BF 3

struct sbox_aprox {
  unsigned char x, y;
  double bias;
};

struct state {
  unsigned char first_sbox;
  unsigned char first_x;
  double biasesMult;
  int cantBiases;
  int position[NUM_SBOXES][SBOX_BITS];
};

struct step {
  unsigned char from;
  int to[NUM_SBOXES][SBOX_BITS];
  struct sbox_aprox path;
};
