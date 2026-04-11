#ifndef PARAMS_H
#define PARAMS_H

/*************************************************
 * Domain separators for XOF expansion
 **************************************************/
// sep
#define DOMAIN_SEPARATOR_A 0
#define DOMAIN_SEPARATOR_R 1
#define DOMAIN_SEPARATOR_A3 2
#define DOMAIN_SEPARATOR_U 3
#define DOMAIN_SEPARATOR_D 4
// osig
#define DOMAIN_SEPARATOR_DS 5
#define DOMAIN_SEPARATOR_S 6
#define DOMAIN_SEPARATOR_RAND 7
#define DOMAIN_SEPARATOR_A1_ISS 8
#define DOMAIN_SEPARATOR_A2_ISS 9
#define DOMAIN_SEPARATOR_BYG_ISS 10
#define DOMAIN_SEPARATOR_B_ISS 11
#define DOMAIN_SEPARATOR_CHAL1_ISS 12
#define DOMAIN_SEPARATOR_CHAL2_ISS 13
#define DOMAIN_SEPARATOR_CHAL3_ISS 14
#define DOMAIN_SEPARATOR_CHAL4_ISS 15
#define DOMAIN_SEPARATOR_RAND_S2_ISS 16
#define DOMAIN_SEPARATOR_RAND_G_ISS 17
// show
#define DOMAIN_SEPARATOR_A1_SHOW 18
#define DOMAIN_SEPARATOR_A2_SHOW 19
#define DOMAIN_SEPARATOR_BYG_SHOW 20
#define DOMAIN_SEPARATOR_B_SHOW 21
#define DOMAIN_SEPARATOR_CHAL1_SHOW 22
#define DOMAIN_SEPARATOR_CHAL2_SHOW 23
#define DOMAIN_SEPARATOR_CHAL3_SHOW 24
#define DOMAIN_SEPARATOR_CHAL4_SHOW 25
#define DOMAIN_SEPARATOR_RAND_S2_SHOW 26
#define DOMAIN_SEPARATOR_RAND_G_SHOW 27

/*************************************************
 * Signature parameters
 **************************************************/
// Ring degree for the signature
#define PARAM_N 256
// Modulus for the signature (q = q_L * q_H = 2401 * 169)
#define PARAM_Q 405769L
// Low modulus (q_L = b_L^k_L = 49^2)
#define PARAM_QL 2401L
// High modulus (q_H = b_H^k_H = 13^2)
#define PARAM_QH 169L
// Modulus bit-length upper bound for uniform sampling
#define PARAM_Q_BITLEN 19
// Module rank for the signature
#define PARAM_D 4
// G_H gadget dimension
#define PARAM_KH 2
// G_L gadget dimension
#define PARAM_KL 2
// G_H gadget base
#define PARAM_BH 13
// G_L gadget base
#define PARAM_BL 49
// Dimension of the message vector (without usk)
#define PARAM_M 10
// Number of iterations for the spectral norm estimation
#define PARAM_IT_SPEC_NORM 6
// Hamming weight of the tags
#define PARAM_W 5
// Bound on the square spectral norm of R
#define PARAM_R_MAX_SQ_SPECTRAL_NORM 3898.89415332639373446000
// Difference s_1^2 - s_L^2 (for perturbation sampler p1_p2[i>=1])
#define PARAM_S1SQ_SLSQ 22816850.56626987829804420471
// Squared Gaussian parameter s_3^2 (for perturbation sampler p3)
#define PARAM_S3SQ 22816838.86761695519089698792
// Gaussian parameter s_4 for v_3
#define PARAM_S4 76.49922142570575545051
// Gaussian parameter s_L/b_L for Z-sampling of zL
#define PARAM_SL_DIV_BL 3.38718927843394679300
// Gaussian parameter s_H/b_H for Z-sampling of zH
#define PARAM_SH_DIV_BH 3.38718927843394679300
// Gaussian width for p_2 (sqrt(s_2^2 - s_L^2))
#define PARAM_SQRT_S2SQ_SLSQ 3.42032935794184789202
// Gaussian width for p_4 (sqrt(s_4^2 - s_H^2))
#define PARAM_SQRT_S4SQ_SHSQ 62.55545719277940008851
// Negated ratio -s_H^2*s_4^2/(s_4^2 - s_H^2) (for covariance)
#define PARAM_NEG_SHSQS4SQ_DIV_S4SQ_SHSQ -2899.67458311796281122952
// Negated ratio -s_H^2/(s_4^2 - s_H^2) (for covariance)
#define PARAM_NEG_SHSQ_DIV_S4SQ_SHSQ -0.49549038516081478356
// Squared verification bound on v_{1,1}
#define PARAM_B11SQ 18149347086L
// Squared verification bound on v_{1,2}
#define PARAM_B12SQ 4652284329L
// Squared verification bound on v_2
#define PARAM_B2SQ 2238819L
// Squared verification bound on v_3
#define PARAM_B3SQ 367252L

// Legacy compatibility constants for OSIG/SHOW layers (to be removed in Phases 6-7)
// These were the old single-gadget SEP parameters: gadget dimension K=5, base B=14
#define PARAM_K 5
#define PARAM_B 14

// Length of the public and secret seeds
#define SEED_BYTES 32
// Length of the public seed for CRS expansion
#define CRS_SEED_BYTES 32
// Length of the state
#define STATE_BYTES 64

/*************************************************
 * [ISSUANCE] Zero-Knowledge proof parameters
 **************************************************/
// Ring degree for the issuance proof
#define PARAM_N_ISS 64
// Ring degree gap between the issuance proof and the signature (subring embedding)
#define PARAM_K_ISS 4
// Modulus for the issuance proof
#define PARAM_Q_ISS 223205310001L
// Modulus bit-length upper bound for uniform sampling
#define PARAM_Q_ISS_BITLEN 38
// Modulus factor for the issuance proof
#define PARAM_Q1_ISS 524201L
// Modulus bit-length upper bound for uniform sampling mod q_1
#define PARAM_Q1_ISS_BITLEN 20
// Second modulus factor for the issuance proof
#define PARAM_Q2_ISS PARAM_Q
// Inverse of q_1 modulo q_2
#define PARAM_Q1_INVMOD_Q2_ISS 343579L
// Inverse of q_2 modulo q_1
#define PARAM_Q2_INVMOD_Q1_ISS 101223L
// Module rank for the issuance proof
#define PARAM_D_ISS 16
// Witness dimension
#define PARAM_M1_ISS 105
// Scaled witness dimension (m_1 / k_hat)
#define PARAM_M1_K_ISS 26
// ABDLOP commitment randomness dimension
#define PARAM_M2_ISS 56
// Soundness amplification dimension
#define PARAM_L_ISS 2
// Dimension for Approximate Range Proof
#define PARAM_ARP_ISS 256
// Rank for Approximate Range Proof (256 / n)
#define PARAM_ARP_DIV_N_ISS 4
// 256 / n + l
#define PARAM_ARP_DIV_N_L_ISS 6
// Gaussian mask width for cs_1
#define PARAM_S1_ISS 369050.89730269293067976832
// Squared Gaussian mask width for cs_1
#define PARAM_S1SQ_ISS 136198564799.92280578613281250000
// Gaussian mask width for cs_2
#define PARAM_S2_ISS 275602.77920886297943070531
// Squared Gaussian mask width for cs_2
#define PARAM_S2SQ_ISS 75956891907.64927673339843750000
// Gaussian mask width for Rs_1 (ARP)
#define PARAM_S3_ISS 72848.10643310110026504844
// Squared Gaussian mask width for Rs_1 (ARP)
#define PARAM_S3SQ_ISS 5306846610.88842582702636718750
// Rejection sampling rate for y_1
#define PARAM_REJ1_ISS 2
// Rejection sampling rate for y_2
#define PARAM_REJ2_ISS 2
// Rejection sampling rate for y_3
#define PARAM_REJ3_ISS 2
// Squared verification bound for z_1
#define PARAM_B1SQ_ISS 180657566352359UL
// Squared verification bound for z_2
#define PARAM_B2SQ_ISS 60411097594469UL
// Squared verification bound for z_3
#define PARAM_B3SQ_ISS 584702794673UL
// Infinity norm of challenges
#define PARAM_RHO_ISS 8
// Manhattan-like norm of challenges
#define PARAM_ETA_ISS 93

/*************************************************
 * [SHOW] Zero-Knowledge proof parameters
 **************************************************/
// Ring degree for the show proof
#define PARAM_N_SHOW 64
// Ring degree gap between the issuance proof and the signature (subring embedding)
#define PARAM_K_SHOW 4
// Modulus for the show proof: qπ₂ = 2^65 - 163 (exceeds 64-bit range; stored in two halves)
// Low 64 bits of qπ₂ = 2^64 - 163
#define PARAM_Q_SHOW_LOW64 18446744073709551453UL
// High 64 bits of qπ₂ (= 1, since qπ₂ = 1*2^64 + (2^64-163))
#define PARAM_Q_SHOW_HIGH64 1UL
// Decimal string for fmpz_set_str in arith_qshow_setup
#define PARAM_Q_SHOW_STR "36893488147419103069"
// Modulus bit-length upper bound for uniform sampling (65 bits)
#define PARAM_Q_SHOW_BITLEN 65
// Modulus factor for the show proof: q₁ = floor(qπ₂ / q) = floor(36893488147419103069 / 405769)
#define PARAM_Q1_SHOW 90922392167511LL
// Modulus bit-length upper bound for uniform sampling mod q_1 (47 bits)
#define PARAM_Q1_SHOW_BITLEN 47
// Second modulus factor for the show proof
#define PARAM_Q2_SHOW PARAM_Q
// Inverse of q_1 modulo q_2: q₁⁻¹ mod q
#define PARAM_Q1_INVMOD_Q2_SHOW 363290L
// Inverse of q_2 modulo q_1: q⁻¹ mod q₁
#define PARAM_Q2_INVMOD_Q1_SHOW 9518450884330LL
// Module rank for the show proof
#define PARAM_D_SHOW 20
// Witness dimension
#define PARAM_M1_SHOW 149
// ABDLOP commitment randomness dimension
#define PARAM_M2_SHOW 72
// Soundness amplification dimension
#define PARAM_L_SHOW 1
// Dimension for Approximate Range Proof
#define PARAM_ARP_SHOW 256
// Rank for Approximate Range Proof (256 / n)
#define PARAM_ARP_DIV_N_SHOW 4
// 256 / n + l
#define PARAM_ARP_DIV_N_L_SHOW 5
// Gaussian mask width σ₁ for cs_1 (Table C.3)
#define PARAM_S1_SHOW 29907260.15800000190734863281
// Squared Gaussian mask width σ₁² for cs_1
#define PARAM_S1SQ_SHOW 894443992211544.87500000000000000000
// Gaussian mask width σ₂ for cs_2 (Table C.3)
#define PARAM_S2_SHOW 13440.06600000000022737367
// Squared Gaussian mask width σ₂² for cs_2
#define PARAM_S2SQ_SHOW 180635374.08404800295829772949
// Gaussian mask width σ₃ for Rs_1 ARP (Table C.3)
#define PARAM_S3_SHOW 1417550244.79400014877319335938
// Squared Gaussian mask width σ₃² for Rs_1 ARP
#define PARAM_S3SQ_SHOW 2009448521065028608.00000000000000000000
// Rejection sampling rate for y_1
#define PARAM_REJ1_SHOW 2
// Rejection sampling rate for y_2
#define PARAM_REJ2_SHOW 2
// Rejection sampling rate for y_3
#define PARAM_REJ3_SHOW 2
// Compression parameter γ (mask commitment), Table C.3
#define PARAM_GAMMA_SHOW 1129442642L
// Compression parameter D (witness commitment), Table C.3
// Note: distinct from PARAM_D_SHOW (module rank = 20)
#define PARAM_D_COM_SHOW 22
// Squared verification bound for z_1 (= REJ1² · 2 · m1 · n · σ₁²; t = 128·ln2, chi-sq tail)
// B1SQ = 10333582971951253504 (fits in uint64, HIGH64 = 0)
#define PARAM_B1SQ_SHOW_LOW64 10333582971951253504UL
// Squared verification bound for z_1 (high 64 bits) — zero since B1SQ < 2^64
#define PARAM_B1SQ_SHOW_HIGH64 0UL
// Squared verification bound for z_2 (= REJ2² · 2 · m2 · n · σ₂²)
// B2SQ = 1095417987986
#define PARAM_B2SQ_SHOW 1095417987986UL
// Squared verification bound for z_3 (= REJ3² · 2 · ARP · σ₃²)
// B3SQ = 1476668782078250450944 > 2^64; split into HIGH64=80, LOW64=929256181486321664
#define PARAM_B3SQ_SHOW_LOW64 929256181486321664UL
#define PARAM_B3SQ_SHOW_HIGH64 80UL
// Infinity norm of challenges
#define PARAM_RHO_SHOW 8
// Manhattan-like norm of challenges
#define PARAM_ETA_SHOW 93
// Proof bound B'_{1,1}² = (B_{1,1} + sqrt(nd))² = (134719.51 + 32)² (Show protocol, §A.3)
#define PARAM_B11_PRIME_SQ ((int64_t)18157969448LL)
// Proof bound B'_{1,2}² = (B_{1,2} + sqrt(nd))² = (68207.66 + 32)² (Show protocol, §A.3)
#define PARAM_B12_PRIME_SQ ((int64_t)4656651197LL)

/*************************************************
 * Testing
 **************************************************/

// The modulus factor
#define Q_1_MOD 47
// The modular inverse of Q_1_MOD mod Q_MIN_MOD (used for CRT reconstruction)
#define Q_1_MOD_INV 2
// The smallest modulus factor
#define Q_MIN_MOD 31
// The modular inverse of Q_MIN_MOD mod Q_1 (used for CRT reconstruction)
#define Q_MIN_MOD_INV 44
// The modulus for the proof system
#define Q_HAT_MOD (Q_1_MOD * Q_MIN_MOD)
// The degree of the ring for proofs
#define N_HAT_RING 4

#endif /* PARAMS_H */
