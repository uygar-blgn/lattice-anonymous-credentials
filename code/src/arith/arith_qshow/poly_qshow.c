#include "arith_qshow.h"
#include "macros.h"

// Shared fmpz_mod context for qπ₂ = 2^65 - 163
static fmpz_mod_ctx_t  QSHOW_CTX;
// f = x^PARAM_N_SHOW + 1
static fmpz_mod_poly_t POLY_F;
// Precomputed reverse-inverse of f for fast mulmod
static fmpz_mod_poly_t POLY_F_REV_INV;

/*************************************************
* Name:        arith_qshow_ctx
*
* Description: Return a pointer to the shared fmpz_mod context for qπ₂.
*              Used by other poly_qshow_* functions that need the context
*              externally (e.g., vec/mat wrappers).
**************************************************/
const fmpz_mod_ctx_t *arith_qshow_ctx(void) {
  return (const fmpz_mod_ctx_t *)&QSHOW_CTX;
}

/*************************************************
* Name:        fmpz_mod_poly_invert  [static]
*
* Description: Compute res = arg^{-1} mod (mod, QSHOW_CTX)
*              using extended GCD.
**************************************************/
static void fmpz_mod_poly_invert(fmpz_mod_poly_t res,
                                  const fmpz_mod_poly_t arg,
                                  const fmpz_mod_poly_t mod) {
  fmpz_mod_poly_t G, S;
  fmpz_mod_poly_init(G, QSHOW_CTX);
  fmpz_mod_poly_init(S, QSHOW_CTX);
  // xgcd(G, S, T, A, B) computes G = S·A + T·B.
  // Pass A=mod, B=arg so G = S·mod + res·arg → res = arg^{-1} mod mod.
  fmpz_mod_poly_xgcd(G, S, res, mod, arg, QSHOW_CTX);
  ASSERT_DEBUG(fmpz_mod_poly_is_one(G, QSHOW_CTX), "GCD != 1");
  fmpz_mod_poly_mulmod(G, res, arg, mod, QSHOW_CTX);
  ASSERT_DEBUG(fmpz_mod_poly_is_one(G, QSHOW_CTX), "arg * res % mod != 1");
  fmpz_mod_poly_clear(G, QSHOW_CTX);
  fmpz_mod_poly_clear(S, QSHOW_CTX);
}

/*************************************************
* Name:        arith_qshow_setup
*
* Description: Initialize and setup the backend for arithmetic
*              modulo qπ₂ = 2^65 - 163.  Must be called once before
*              any other poly_qshow_* function is used.
**************************************************/
void arith_qshow_setup(void) {
  // Initialize context with qπ₂ from its decimal string representation
  fmpz_t q;
  fmpz_init(q);
  fmpz_set_str(q, PARAM_Q_SHOW_STR, 10);
  fmpz_mod_ctx_init(QSHOW_CTX, q);
  fmpz_clear(q);

  // Build f = x^PARAM_N_SHOW + 1
  fmpz_mod_poly_init(POLY_F, QSHOW_CTX);
  fmpz_mod_poly_set_coeff_ui(POLY_F, 0,            1, QSHOW_CTX);
  fmpz_mod_poly_set_coeff_ui(POLY_F, PARAM_N_SHOW, 1, QSHOW_CTX);

  // Build f_len = x^(PARAM_N_SHOW+1) (used as modulus for xgcd)
  fmpz_mod_poly_t f_len;
  fmpz_mod_poly_init(f_len, QSHOW_CTX);
  fmpz_mod_poly_set_coeff_ui(f_len, PARAM_N_SHOW + 1, 1, QSHOW_CTX);

  // Build f_rev = reverse(f, PARAM_N_SHOW+1)
  fmpz_mod_poly_t f_rev;
  fmpz_mod_poly_init(f_rev, QSHOW_CTX);
  fmpz_mod_poly_reverse(f_rev, POLY_F, PARAM_N_SHOW + 1, QSHOW_CTX);

  // Compute POLY_F_REV_INV = f_rev^{-1} mod f_len
  fmpz_mod_poly_init(POLY_F_REV_INV, QSHOW_CTX);
  fmpz_mod_poly_invert(POLY_F_REV_INV, f_rev, f_len);

  fmpz_mod_poly_clear(f_len, QSHOW_CTX);
  fmpz_mod_poly_clear(f_rev, QSHOW_CTX);
}

/*************************************************
* Name:        arith_qshow_teardown
*
* Description: Clean up and release resources for the qπ₂ arithmetic backend.
*              Must be called once at the end.
**************************************************/
void arith_qshow_teardown(void) {
  fmpz_mod_poly_clear(POLY_F,         QSHOW_CTX);
  fmpz_mod_poly_clear(POLY_F_REV_INV, QSHOW_CTX);
  fmpz_mod_ctx_clear(QSHOW_CTX);
}

/*************************************************
* Name:        poly_qshow_init
*
* Description: Initialize polynomial and set it to zero.
*              Must be called before any operations on the polynomial.
*
* Arguments:   - poly_qshow res: polynomial to be initialized
**************************************************/
void poly_qshow_init(poly_qshow res) {
  fmpz_mod_poly_init(res, QSHOW_CTX);
}

/*************************************************
* Name:        poly_qshow_clear
*
* Description: Clear a polynomial and release all associated memory.
*
* Arguments:   - poly_qshow arg: polynomial to be cleared
**************************************************/
void poly_qshow_clear(poly_qshow arg) {
  fmpz_mod_poly_clear(arg, QSHOW_CTX);
}

/*************************************************
* Name:        poly_qshow_zero
*
* Description: Set an initialized polynomial to zero.
*
* Arguments:   - poly_qshow res: polynomial to be zeroed (initialized)
**************************************************/
void poly_qshow_zero(poly_qshow res) {
  fmpz_mod_poly_zero(res, QSHOW_CTX);
}

/*************************************************
* Name:        poly_qshow_set
*
* Description: Set a polynomial equal to another polynomial.
*
* Arguments:   - poly_qshow res: output polynomial (initialized)
*              - const poly_qshow arg: input polynomial
**************************************************/
void poly_qshow_set(poly_qshow res, const poly_qshow arg) {
  fmpz_mod_poly_set(res, arg, QSHOW_CTX);
}

/*************************************************
* Name:        int128_from_fmpz  [static]
*
* Description: Convert a non-negative fmpz value in [0, qπ₂-1] to int128.
*              Uses GMP multi-limb extraction for values ≥ 2^63.
**************************************************/
static coeff_qshow int128_from_fmpz(const fmpz_t c) {
  if (fmpz_fits_si(c)) {
    return (coeff_qshow)fmpz_get_si(c);
  }
  // Value is in [2^63, qπ₂-1] ⊂ [2^63, 2^65-1]: needs two 64-bit limbs
  mpz_t z;
  mpz_init(z);
  fmpz_get_mpz(z, c);
  mp_limb_t lo = mpz_getlimbn(z, 0);
  mp_limb_t hi = (mpz_size(z) > 1) ? mpz_getlimbn(z, 1) : 0;
  coeff_qshow res = (coeff_qshow)(((uint128)hi << 64) | (uint128)lo);
  mpz_clear(z);
  return res;
}

/*************************************************
* Name:        fmpz_from_int128  [static]
*
* Description: Set fmpz f to the value represented by a coeff_qshow (int128).
*              Handles both positive and negative values.
*              The caller must then reduce mod qπ₂ if needed.
**************************************************/
static void fmpz_from_int128(fmpz_t f, coeff_qshow c) {
  if (c >= 0) {
    uint128 uc = (uint128)c;
    ulong hi = (ulong)(uc >> 64);
    ulong lo = (ulong)(uc & (uint128)UINT64_MAX);
    if (hi == 0) {
      fmpz_set_ui(f, lo);
    } else {
      fmpz_set_ui(f, hi);
      fmpz_mul_2exp(f, f, 64);
      fmpz_add_ui(f, f, lo);
    }
  } else {
    // Negate, convert positive part, then negate the fmpz
    int128 pos = -c;
    uint128 uc = (uint128)pos;
    ulong hi = (ulong)(uc >> 64);
    ulong lo = (ulong)(uc & (uint128)UINT64_MAX);
    if (hi == 0) {
      fmpz_set_ui(f, lo);
    } else {
      fmpz_set_ui(f, hi);
      fmpz_mul_2exp(f, f, 64);
      fmpz_add_ui(f, f, lo);
    }
    fmpz_neg(f, f);
  }
}

/*************************************************
* Name:        poly_qshow_get_coeff
*
* Description: Get coefficient of x^n of a polynomial in [0, qπ₂-1].
*              condition: 0 <= n < PARAM_N_SHOW
*
* Arguments:   - const poly_qshow arg: polynomial to be read
*              - size_t n: degree of the coefficient
*
* Returns:     coeff_qshow (int128) value in [0, qπ₂-1]
**************************************************/
coeff_qshow poly_qshow_get_coeff(const poly_qshow arg, size_t n) {
  ASSERT_DEBUG(n < PARAM_N_SHOW,
    "Illegal argument: cannot get coefficient of poly at given position.");
  fmpz_t c;
  fmpz_init(c);
  fmpz_mod_poly_get_coeff_fmpz(c, arg, (slong)n, QSHOW_CTX);
  coeff_qshow res = int128_from_fmpz(c);
  fmpz_clear(c);
  return res;
}

/*************************************************
* Name:        poly_qshow_get_coeff_centered
*
* Description: Get coefficient of x^n in centered representation [-qπ₂/2, qπ₂/2].
*              condition: 0 <= n < PARAM_N_SHOW
*
* Arguments:   - const poly_qshow arg: polynomial to be read
*              - size_t n: degree of the coefficient
*
* Returns:     coeff_qshow in centered range
**************************************************/
coeff_qshow poly_qshow_get_coeff_centered(const poly_qshow arg, size_t n) {
  static const uint128 Q  = ((uint128)PARAM_Q_SHOW_HIGH64 << 64)
                          |  (uint128)PARAM_Q_SHOW_LOW64;
  static const uint128 QH = Q / 2;
  coeff_qshow c = poly_qshow_get_coeff(arg, n);  // c in [0, Q-1]
  uint128 uc = (uint128)c;
  return (coeff_qshow)(uc > QH ? (int128)(uc - Q) : (int128)uc);
}

/*************************************************
* Name:        poly_qshow_set_coeff
*
* Description: Set coefficient of x^n of a polynomial.
*              condition: 0 <= n < PARAM_N_SHOW
*              The coefficient is reduced mod qπ₂ into [0, qπ₂-1].
*
* Arguments:   - poly_qshow arg: polynomial (initialized)
*              - size_t n: degree of the coefficient
*              - coeff_qshow c: new coefficient value (may be negative or ≥ qπ₂)
**************************************************/
void poly_qshow_set_coeff(poly_qshow arg, size_t n, coeff_qshow c) {
  ASSERT_DEBUG(n < PARAM_N_SHOW,
    "Illegal argument: cannot set coefficient of poly at given position.");
  static const uint128 Q = ((uint128)PARAM_Q_SHOW_HIGH64 << 64)
                         |  (uint128)PARAM_Q_SHOW_LOW64;
  // Reduce c into [0, Q-1] using unsigned 128-bit arithmetic
  // Cast signed int128 to unsigned, add Q to handle negatives, reduce mod Q
  uint128 uc = ((uint128)((int128)c % (int128)Q) + Q) % Q;

  fmpz_t f;
  fmpz_init(f);
  ulong hi = (ulong)(uc >> 64);
  ulong lo = (ulong)(uc & (uint128)UINT64_MAX);
  if (hi == 0) {
    fmpz_set_ui(f, lo);
  } else {
    fmpz_set_ui(f, hi);
    fmpz_mul_2exp(f, f, 64);
    fmpz_add_ui(f, f, lo);
  }
  fmpz_mod_poly_set_coeff_fmpz(arg, (slong)n, f, QSHOW_CTX);
  fmpz_clear(f);
}

/*************************************************
* Name:        poly_qshow_neg
*
* Description: Negate a polynomial coefficient-wise mod qπ₂.
*
* Arguments:   - poly_qshow res: output polynomial (initialized)
*              - const poly_qshow arg: polynomial to be negated
**************************************************/
void poly_qshow_neg(poly_qshow res, const poly_qshow arg) {
  fmpz_mod_poly_neg(res, arg, QSHOW_CTX);
}

/*************************************************
* Name:        poly_qshow_add
*
* Description: Add two polynomials mod qπ₂.
*
* Arguments:   - poly_qshow res: output polynomial (initialized)
*              - const poly_qshow lhs: first summand
*              - const poly_qshow rhs: second summand
**************************************************/
void poly_qshow_add(poly_qshow res, const poly_qshow lhs, const poly_qshow rhs) {
  fmpz_mod_poly_add(res, lhs, rhs, QSHOW_CTX);
}

/*************************************************
* Name:        poly_qshow_sub
*
* Description: Subtract two polynomials mod qπ₂.
*
* Arguments:   - poly_qshow res: output polynomial (initialized)
*              - const poly_qshow lhs: minuend
*              - const poly_qshow rhs: subtrahend
**************************************************/
void poly_qshow_sub(poly_qshow res, const poly_qshow lhs, const poly_qshow rhs) {
  fmpz_mod_poly_sub(res, lhs, rhs, QSHOW_CTX);
}

/*************************************************
* Name:        poly_qshow_mul
*
* Description: Multiply two polynomials reduced modulo x^PARAM_N_SHOW + 1
*              and mod qπ₂.
*
* Arguments:   - poly_qshow res: output polynomial (initialized)
*              - const poly_qshow lhs: first factor
*              - const poly_qshow rhs: second factor
**************************************************/
void poly_qshow_mul(poly_qshow res, const poly_qshow lhs, const poly_qshow rhs) {
  ASSERT_DEBUG(fmpz_mod_poly_degree(lhs, QSHOW_CTX) < PARAM_N_SHOW,
    "Argument to `poly_qshow_mul` must already be reduced for FLINT.");
  ASSERT_DEBUG(fmpz_mod_poly_degree(rhs, QSHOW_CTX) < PARAM_N_SHOW,
    "Argument to `poly_qshow_mul` must already be reduced for FLINT.");
  fmpz_mod_poly_mulmod_preinv(res, lhs, rhs, POLY_F, POLY_F_REV_INV, QSHOW_CTX);
}

/*************************************************
* Name:        poly_qshow_mul_x
*
* Description: Multiply a polynomial by x, reduced modulo x^PARAM_N_SHOW + 1
*              and mod qπ₂.
*
* Arguments:   - poly_qshow res: output polynomial (initialized)
*              - const poly_qshow arg: input polynomial
**************************************************/
void poly_qshow_mul_x(poly_qshow res, const poly_qshow arg) {
  // Shift left by 1: degree may reach PARAM_N_SHOW
  fmpz_mod_poly_shift_left(res, arg, 1, QSHOW_CTX);
  // Wrap: coefficient of x^PARAM_N_SHOW contributes -coeff to x^0 (since x^N = -1)
  fmpz_t top;
  fmpz_init(top);
  fmpz_mod_poly_get_coeff_fmpz(top, res, PARAM_N_SHOW, QSHOW_CTX);
  if (!fmpz_is_zero(top)) {
    // coeff[0] -= top  (i.e., coeff[0] += Q - top, then reduce)
    fmpz_t c0;
    fmpz_init(c0);
    fmpz_mod_poly_get_coeff_fmpz(c0, res, 0, QSHOW_CTX);
    fmpz_sub(c0, c0, top);
    fmpz_mod(c0, c0, fmpz_mod_ctx_modulus(QSHOW_CTX));
    fmpz_mod_poly_set_coeff_fmpz(res, 0, c0, QSHOW_CTX);
    fmpz_clear(c0);
    // Clear the degree-N coefficient
    fmpz_zero(top);
    fmpz_mod_poly_set_coeff_fmpz(res, PARAM_N_SHOW, top, QSHOW_CTX);
  }
  fmpz_clear(top);
}

/*************************************************
* Name:        poly_qshow_mul_scalar
*
* Description: Multiply a polynomial by a scalar (coeff_qshow = int128).
*
* Arguments:   - poly_qshow out: output polynomial (initialized)
*              - const poly_qshow lhs: input polynomial
*              - coeff_qshow rhs: scalar multiplier
**************************************************/
void poly_qshow_mul_scalar(poly_qshow out, const poly_qshow lhs,
                            const coeff_qshow rhs) {
  fmpz_t s;
  fmpz_init(s);
  fmpz_from_int128(s, rhs);
  fmpz_mod(s, s, fmpz_mod_ctx_modulus(QSHOW_CTX));
  fmpz_mod_poly_scalar_mul_fmpz(out, lhs, s, QSHOW_CTX);
  fmpz_clear(s);
}

/*************************************************
* Name:        poly_qshow_muladd_constant
*
* Description: Increment constant coefficient (x^0) of a polynomial
*              by the product of two scalars: arg[0] += c0_lhs * c0_rhs mod qπ₂.
*
* Arguments:   - poly_qshow arg: polynomial to be incremented (initialized)
*              - const coeff_qshow c0_lhs: first scalar factor
*              - const coeff_qshow c0_rhs: second scalar factor
**************************************************/
void poly_qshow_muladd_constant(poly_qshow arg,
                                 const coeff_qshow c0_lhs,
                                 const coeff_qshow c0_rhs) {
  // Product in 128-bit signed arithmetic, then reduce mod qπ₂
  int128 tmp = (int128)c0_lhs * (int128)c0_rhs;

  // Convert product to fmpz and reduce
  fmpz_t f, coeff;
  fmpz_init(f);
  fmpz_init(coeff);
  fmpz_from_int128(f, (coeff_qshow)tmp);
  fmpz_mod(f, f, fmpz_mod_ctx_modulus(QSHOW_CTX));

  // Add to existing constant coefficient
  fmpz_mod_poly_get_coeff_fmpz(coeff, arg, 0, QSHOW_CTX);
  fmpz_add(coeff, coeff, f);
  fmpz_mod(coeff, coeff, fmpz_mod_ctx_modulus(QSHOW_CTX));
  fmpz_mod_poly_set_coeff_fmpz(arg, 0, coeff, QSHOW_CTX);

  fmpz_clear(f);
  fmpz_clear(coeff);
}

/*************************************************
* Name:        poly_qshow_shift_left
*
* Description: Shift polynomial coefficients left by n places (multiply by x^n),
*              NOT reduced mod x^PARAM_N_SHOW + 1.
*
* Arguments:   - poly_qshow res: output polynomial (initialized)
*              - const poly_qshow arg: input polynomial
*              - size_t n: shift amount
**************************************************/
void poly_qshow_shift_left(poly_qshow res, const poly_qshow arg, size_t n) {
  fmpz_mod_poly_shift_left(res, arg, (slong)n, QSHOW_CTX);
}

/*************************************************
* Name:        poly_qshow_conjugate
*
* Description: Compute the conjugate of a polynomial (evaluation at x^{-1}).
*              In the ring Z[x]/(x^N+1), x^{-1} = -x^{N-1}, so the conjugate
*              of sum a_i x^i is a_0 - sum_{i=1}^{N-1} a_i x^{N-i}.
*
* Arguments:   - poly_qshow out: output polynomial (initialized)
*              - const poly_qshow arg: input polynomial
**************************************************/
void poly_qshow_conjugate(poly_qshow out, const poly_qshow arg) {
  // coeff 0 is unchanged; coeff N/2 negates; rest reverse-negate
  coeff_qshow c = poly_qshow_get_coeff(arg, 0);
  coeff_qshow d = poly_qshow_get_coeff(arg, PARAM_N_SHOW / 2);
  poly_qshow_set_coeff(out, 0, c);
  poly_qshow_set_coeff(out, PARAM_N_SHOW / 2, -(coeff_qshow)d);
  for (size_t i = 1; i < PARAM_N_SHOW / 2; ++i) {
    c = poly_qshow_get_coeff(arg, PARAM_N_SHOW - i);
    d = poly_qshow_get_coeff(arg, i);
    poly_qshow_set_coeff(out,                i, -(coeff_qshow)c);
    poly_qshow_set_coeff(out, PARAM_N_SHOW - i, -(coeff_qshow)d);
  }
}

/*************************************************
* Name:        poly_qshow_equal
*
* Description: Test equality of two polynomials.
*
* Arguments:   - const poly_qshow lhs: first polynomial
*              - const poly_qshow rhs: second polynomial
*
* Returns 1 if equal, 0 otherwise.
**************************************************/
int poly_qshow_equal(const poly_qshow lhs, const poly_qshow rhs) {
  return fmpz_mod_poly_equal(lhs, rhs, QSHOW_CTX);
}

/*************************************************
* Name:        poly_qshow_dump
*
* Description: Print a polynomial.
*
* Arguments:   - const poly_qshow arg: polynomial to be printed
**************************************************/
void poly_qshow_dump(const poly_qshow arg) {
  fmpz_mod_poly_print(arg, QSHOW_CTX);
}

/*************************************************
* Name:        poly_qshow_sq_norm2
*
* Description: Compute the squared ℓ₂ norm of a polynomial using centered
*              representation.
*
* Arguments:   - const poly_qshow arg: input polynomial
*
* Returns uint128 with the squared ℓ₂ norm.
**************************************************/
uint128 poly_qshow_sq_norm2(const poly_qshow arg) {
  uint128 sq_norm2 = 0;
  for (size_t i = 0; i < PARAM_N_SHOW; i++) {
    coeff_qshow c = poly_qshow_get_coeff_centered(arg, i);
    // c is Gaussian-distributed; products fit well within uint128 in practice
    sq_norm2 += (uint128)((int128)c * (int128)c);
  }
  return sq_norm2;
}

/*************************************************
* Name:        poly_qshow_pack
*
* Description: Pack a polynomial mod qπ₂ into a byte array.
*              Uses 9 bytes per coefficient (65-bit value).
*
* Arguments:   - uint8_t buf[POLYQSHOW_PACKEDBYTES]: output buffer
*              - const poly_qshow arg: input polynomial
**************************************************/
void poly_qshow_pack(uint8_t buf[POLYQSHOW_PACKEDBYTES], const poly_qshow arg) {
  for (size_t i = 0; i < PARAM_N_SHOW; i++) {
    // Get coefficient as uint128 (value in [0, qπ₂-1] < 2^65)
    uint128 x = (uint128)poly_qshow_get_coeff(arg, i);
    for (int b = 0; b < COEFFQSHOW_PACKEDBYTES; b++) {
      buf[COEFFQSHOW_PACKEDBYTES * i + b] = (uint8_t)(x & 0xFF);
      x >>= 8;
    }
  }
}

/*************************************************
* Name:        coeff_qshow_pack
*
* Description: Pack a single coefficient mod qπ₂ into a byte array.
*              Uses 9 bytes (65-bit value).  Negative centered values are
*              mapped to their positive representative in [0, qπ₂-1] first.
*
* Arguments:   - uint8_t buf[COEFFQSHOW_PACKEDBYTES]: output buffer
*              - const coeff_qshow arg: input coefficient
**************************************************/
void coeff_qshow_pack(uint8_t buf[COEFFQSHOW_PACKEDBYTES],
                       const coeff_qshow arg) {
  static const uint128 Q = ((uint128)PARAM_Q_SHOW_HIGH64 << 64)
                         |  (uint128)PARAM_Q_SHOW_LOW64;
  // Map arg into [0, Q-1]
  uint128 uc = ((uint128)((int128)arg % (int128)Q) + Q) % Q;
  for (int b = 0; b < COEFFQSHOW_PACKEDBYTES; b++) {
    buf[b] = (uint8_t)(uc & 0xFF);
    uc >>= 8;
  }
}
