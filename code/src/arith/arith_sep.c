/* Lightweight arith setup/teardown and cross-type helpers for sep-only builds.
 * Initializes only the q, real, and z rings needed by sep.c / covariance.c / sampling.c.
 * Does not touch the not-yet-updated qiss/qshow layers. */

#include "arith.h"
#include "covariance.h"

void poly_z_mat_d_d_from_poly_q_mat_d_d(poly_z_mat_d_d res, const poly_q_mat_d_d arg) {
  size_t i,j,k;
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_D; j++) {
      for (k = 0; k < PARAM_N; k++) {
        poly_z_set_coeff_si(res->rows[i]->entries[j], k,
          poly_q_get_coeff_centered(arg->rows[i]->entries[j], k));
      }
    }
  }
}

void poly_real_from_poly_q(poly_real res, const poly_q arg) {
  size_t i;
  for (i = 0; i < PARAM_N; i++) {
    poly_real_set_si(res, i, poly_q_get_coeff_centered(arg, i));
  }
}

void poly_q_from_poly_real(poly_q res, const poly_real arg) {
  size_t i;
  for (i = 0; i < PARAM_N; i++) {
    poly_q_set_coeff(res, i, poly_real_get_coeff_rounded(arg, i));
  }
}

void poly_q_samplefz(poly_q res, const poly_real f, const poly_real c) {
  poly_real tmp;
  poly_real_init(tmp);
  poly_real_samplefz(tmp, f, c);
  poly_q_from_poly_real(res, tmp);
  poly_real_clear(tmp);
}

void poly_real_sub_poly_real_poly_q(poly_real res, const poly_q lhs, const poly_real rhs) {
  poly_real tmp;
  poly_real_init(tmp);
  poly_real_from_poly_q(tmp, lhs);
  poly_real_sub(res, tmp, rhs);
  poly_real_clear(tmp);
}

void arith_setup(void) {
  arith_q_setup();
  arith_real_setup();
  arith_z_setup();
  fft_precomp_setup();

  poly_q_vec_d_setup();
  poly_q_vec_m_setup();
  poly_q_vec_k_setup();
  poly_q_mat_d_m_setup();
  poly_q_mat_d_k_setup();
  poly_q_mat_d_d_setup();

  poly_z_vec_d_setup();
  poly_z_mat_d_d_setup();
}

void arith_teardown(void) {
  arith_q_teardown();
  arith_real_teardown();
  arith_z_teardown();
  fft_precomp_teardown();

  poly_q_vec_d_teardown();
  poly_q_vec_m_teardown();
  poly_q_vec_k_teardown();
  poly_q_mat_d_m_teardown();
  poly_q_mat_d_k_teardown();
  poly_q_mat_d_d_teardown();

  poly_z_vec_d_teardown();
  poly_z_mat_d_d_teardown();
}
