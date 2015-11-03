// BORROWED FROM THE SOFTWARE JELLYFISH

#ifndef JELLYFISH_LIB
#define JELLYFISH_LIB

#define R -1
#define I -2
#define O -3
#define MER_A 0
#define C 1
#define G 2
#define MER_T 3
const char dna_codes[256] = {
  O, O, O, O, O, O, O, O, O, O, I, O, O, O, O, O,
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
  O, O, O, O, O, O, O, O, O, O, O, O, O, R, O, O,
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
  O, MER_A, R, C, R, O, O, G, R, O, O, R, O, R, R, O,
  O, O, R, R, MER_T, O, R, R, R, R, O, O, O, O, O, O,
  O, MER_A, R, C, R, O, O, G, R, O, O, R, O, R, R, O,
  O, O, R, R, MER_T, O, R, R, R, R, O, O, O, O, O, O,
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O
};


static const size_t CODE_A = 0;
static const size_t CODE_C = 0;
static const size_t CODE_G = 0;
static const size_t CODE_T = 0;
// Non DNA codes have the MSB on
static const size_t CODE_RESET   = (size_t)-1;
static const size_t CODE_IGNORE  = (size_t)-2;
static const size_t CODE_COMMENT = (size_t)-3;
static const size_t CODE_NOT_DNA = ((size_t)1) << (sizeof(size_t) - 1);

static uint64_t mer_string_to_binary(const char *in, size_t klen) {
  uint64_t res = 0;
  for(size_t i = 0; i < klen; i++) {
    const size_t c = dna_codes[(size_t)*in++];
    if(c & CODE_NOT_DNA)
      return 0;
    res = (res << 2) | c;
  }
  return res;
}


static void mer_binary_to_string(uint64_t mer, size_t klen, char *out) {
  static const char table[4] = { 'A', 'C', 'G', 'T' };

  for(unsigned int i = 0 ; i < klen; i++) {
    out[klen-1-i] = table[mer & (uint64_t)0x3];
    mer >>= 2;
  }
  out[klen] = '\0';
}


static string mer_binary_to_string(uint64_t mer, size_t klen) {
  static const char table[4] = { 'A', 'C', 'G', 'T' };

  string out(klen, '-');

  for(unsigned int i = 0 ; i < klen; i++) {
    out[klen-1-i] = table[mer & (uint64_t)0x3];
    mer >>= 2;
  }
  return out;
}

#endif
