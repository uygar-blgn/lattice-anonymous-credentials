#include <stdio.h>
#include "sep.h"
#include "randombytes.h"
#include "random.h"

#define NTESTS 1
#define NSUBTESTS 3

static int sep_test(void)
{
  int rval = 1;
  sep_sk_t sk;
  sep_pk_t pk;
  sep_sig_t sig;
  uint8_t state[STATE_BYTES], msg[PARAM_M*PARAM_N/8];
  randombytes(state, STATE_BYTES);

  printf("\nsep_test\n");

  sep_keys_init(&pk, &sk);
  sep_sig_init(&sig);

  for (int i = 0; i < NSUBTESTS; i++)
  {
    sep_keygen(&pk, &sk);
    for (int j = 0; j < NSUBTESTS; j++)
    {
      randombytes(msg, PARAM_M*PARAM_N/8);
      sep_sign(&sig, state, &sk, &pk, msg);
      if (!sep_verify(&sig, msg, &pk))
      {
        printf("FAIL: sep_verify returned 0 for a valid signature.\n");
        rval = 0;
        goto sep_test_cleanup;
      }
      msg[0] ^= 1;
      if (sep_verify(&sig, msg, &pk))
      {
        printf("FAIL: sep_verify returned 1 for a wrong message.\n");
        rval = 0;
        goto sep_test_cleanup;
      }
      printf(".");
      fflush(stdout);
    }
  }

sep_test_cleanup:
  sep_sig_clear(&sig);
  sep_keys_clear(&pk, &sk);
  return rval;
}

int main(void) {
  arith_setup();
  random_init();
  printf("sep sign/verify test\n");
  int pass = 1;
  for (int i = 0; i < NTESTS; i++)
  {
    pass &= sep_test();
    if (!pass) {
      printf("\nFAILED\n");
      break;
    }
  }
  if (pass) {
    printf("\npassed\n");
  }
  arith_teardown();
  return pass ? 0 : 1;
}
