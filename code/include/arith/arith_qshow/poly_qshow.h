#ifndef POLY_QSHOW_H
#define POLY_QSHOW_H

#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz.h>
#include <gmp.h>

#include "params.h"

#ifdef __SIZEOF_INT128__
__extension__ typedef __int128 int128;
__extension__ typedef unsigned __int128 uint128;
#endif

// 9 bytes per coefficient: ceil(65/8) = 9 (qπ₂ requires 65 bits)
#define COEFFQSHOW_PACKEDBYTES 9
#define POLYQSHOW_PACKEDBYTES (PARAM_N_SHOW * COEFFQSHOW_PACKEDBYTES)

// Backing polynomial type: arbitrary-precision modulus via fmpz_mod_poly
typedef fmpz_mod_poly_t poly_qshow;

// Coefficient type: int128 to hold values in [0, qπ₂-1] where qπ₂ ≈ 2^65
typedef int128 coeff_qshow;

// Accessor for the shared fmpz_mod context (initialized in arith_qshow_setup)
const fmpz_mod_ctx_t *arith_qshow_ctx(void);

void arith_qshow_setup(void);
void arith_qshow_teardown(void);

void poly_qshow_init(poly_qshow res);
void poly_qshow_clear(poly_qshow arg);
void poly_qshow_zero(poly_qshow res);
void poly_qshow_set(poly_qshow res, const poly_qshow arg);
coeff_qshow poly_qshow_get_coeff(const poly_qshow arg, size_t n);
coeff_qshow poly_qshow_get_coeff_centered(const poly_qshow arg, size_t n);
void poly_qshow_set_coeff(poly_qshow arg, size_t n, coeff_qshow c);
void poly_qshow_neg(poly_qshow res, const poly_qshow arg);
void poly_qshow_add(poly_qshow res, const poly_qshow lhs, const poly_qshow rhs);
void poly_qshow_neg(poly_qshow res, const poly_qshow arg);
void poly_qshow_sub(poly_qshow res, const poly_qshow lhs, const poly_qshow rhs);
void poly_qshow_mul(poly_qshow res, const poly_qshow lhs, const poly_qshow rhs);
void poly_qshow_mul_x(poly_qshow res, const poly_qshow arg);
void poly_qshow_mul_scalar(poly_qshow out, const poly_qshow lhs, const coeff_qshow rhs);
void poly_qshow_muladd_constant(poly_qshow arg, const coeff_qshow c0_lhs, const coeff_qshow c0_rhs);
void poly_qshow_shift_left(poly_qshow res, const poly_qshow arg, size_t n);
void poly_qshow_conjugate(poly_qshow out, const poly_qshow arg);
int poly_qshow_equal(const poly_qshow lhs, const poly_qshow rhs);
void poly_qshow_dump(const poly_qshow arg);
uint128 poly_qshow_sq_norm2(const poly_qshow arg);
void poly_qshow_pack(uint8_t buf[POLYQSHOW_PACKEDBYTES], const poly_qshow arg);
void coeff_qshow_pack(uint8_t buf[COEFFQSHOW_PACKEDBYTES], const coeff_qshow arg);

#endif /* POLY_QSHOW_H */
