#include <stdio.h>
#include "arith.h"
#include "sep.h"
#include "randombytes.h"
#include "random.h"

#define NTESTS 1
#define NSUBTESTS 5

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
        printf("sep_verify returned zero for a valid signature.\n");
        rval = 0;
        goto sep_test_cleanup;
      }
      msg[0] ^= 1;
      if (sep_verify(&sig, msg, &pk))
      {
        printf("sep_verify returned non-zero for a signature on the wrong message.\n");
        rval = 0;
        goto sep_test_cleanup;
      }
      printf(":");
      fflush(stdout);
    }
  }

  sep_test_cleanup:
  sep_sig_clear(&sig);
  sep_keys_clear(&pk, &sk);
  return rval;
}

int main(void)
{
  int rc = 0;

  arith_setup();
  random_init();

  printf("Testing SEP with TSampler (Phases 1-5)\n");
  printf("=======================================\n");

  if (!sep_test())
  {
    printf("\nFAIL: sep_test\n");
    rc = 1;
  }
  else
  {
    printf("\nPASS: sep_test\n");
  }

  arith_teardown();
  return rc;
}
