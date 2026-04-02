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
// Modulus for the signature (= q_L * q_H = 2401 * 169)
#define PARAM_Q 405769L
// Low modulus factor q_L = 7^4
#define PARAM_QL 2401L
// High modulus factor q_H = 13^2
#define PARAM_QH 169L
// Modulus bit-length upper bound for uniform sampling
#define PARAM_Q_BITLEN 19
// Module rank for the signature
#define PARAM_D 4
// High gadget dimension (k_H)
#define PARAM_KH 2
// Low gadget dimension (k_L)
#define PARAM_KL 2
// High gadget base (b_H)
#define PARAM_BH 13
// Low gadget base (b_L)
#define PARAM_BL 49
// Dimension of the message vector (without usk)
#define PARAM_M 10
// Number of iterations for the spectral norm estimation
#define PARAM_IT_SPEC_NORM 5
// Hamming weight of the tags
#define PARAM_W 5
// Bound on the square spectral norm of R
#define PARAM_R_MAX_SQ_SPECTRAL_NORM 3898.89415332639373446000
// Squared difference s_1^2 - s_L^2 (perturbation sampler)
#define PARAM_S1SQ_SLSQ 22816850.56626987829804420471
// Squared Gaussian width s_3^2 (for v_12)
#define PARAM_S3SQ 22816838.86761695519089698792
// Gaussian width s_4 (for v_2 and v_3)
#define PARAM_S4 76.49922142570575545051
// Ratio s_L / b_L
#define PARAM_SL_DIV_BL 3.38718927843394723709
// Ratio s_H / b_H
#define PARAM_SH_DIV_BH 3.38718927843394679300
// sqrt(s_2^2 - s_L^2)
#define PARAM_SQRT_S2SQ_SLSQ 3.42032935794184789202
// sqrt(s_4^2 - s_H^2)
#define PARAM_SQRT_S4SQ_SHSQ 62.55545719277940008851
// -s_H^2 * s_4^2 / (s_4^2 - s_H^2)
#define PARAM_NEG_SHSQS4SQ_DIV_S4SQ_SHSQ -2899.67458311796281122952
// -s_H^2 / (s_4^2 - s_H^2)
#define PARAM_NEG_SHSQ_DIV_S4SQ_SHSQ -0.49549038516081478356
// Squared verification bound on v_11
#define PARAM_B11SQ 18149347086L
// Squared verification bound on v_12
#define PARAM_B12SQ 4652284329L
// Squared verification bound on v_2
#define PARAM_B2SQ 2238819L
// Squared verification bound on v_3
#define PARAM_B3SQ 367252L

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
// Modulus for the issuance proof (single prime)
#define PARAM_Q_ISS 562949953421189L
// Modulus bit-length upper bound for uniform sampling
#define PARAM_Q_ISS_BITLEN 49
// Module rank for the issuance proof
#define PARAM_D_ISS 20
// Witness dimension
#define PARAM_M1_ISS 105
// Scaled witness dimension (m_1 / k_hat)
#define PARAM_M1_K_ISS 26
// ABDLOP commitment randomness dimension
#define PARAM_M2_ISS 58
// Soundness amplification dimension
#define PARAM_L_ISS 2
// Dimension for Approximate Range Proof
#define PARAM_ARP_ISS 256
// Rank for Approximate Range Proof (256 / n)
#define PARAM_ARP_DIV_N_ISS 4
// 256 / n + l
#define PARAM_ARP_DIV_N_L_ISS 6
// Gaussian mask width for cs_1
#define PARAM_S1_ISS 16230.42121646810664969962
// Squared Gaussian mask width for cs_1
#define PARAM_S1SQ_ISS 263426572.86397805809974670410
// Gaussian mask width for cs_2
#define PARAM_S2_ISS 12062.82807641120962216519
// Squared Gaussian mask width for cs_2
#define PARAM_S2SQ_ISS 145511821.20105457305908203125
// Gaussian mask width for Rs_1 (ARP)
#define PARAM_S3_ISS 3522387.00453344406560063362
// Squared Gaussian mask width for Rs_1 (ARP)
#define PARAM_S3SQ_ISS 12407210209706.08984375000000000000
// Rejection sampling rate for y_1
#define PARAM_REJ1_ISS 2
// Rejection sampling rate for y_2
#define PARAM_REJ2_ISS 2
// Rejection sampling rate for y_3
#define PARAM_REJ3_ISS 2
// Squared verification bound for z_1
#define PARAM_B1SQ_ISS 1409647146440UL
// Squared verification bound for z_2
#define PARAM_B2SQ_ISS 462921986980UL
// Infinity norm of challenges
#define PARAM_RHO_ISS 8
// Manhattan-like norm of challenges
#define PARAM_ETA_ISS 93
// Squared verification bound for z_3 (ARP integers); ≈ 2 * PARAM_ARP_ISS * PARAM_S3SQ_ISS
#define PARAM_B3SQ_ISS 6352491635051564UL

/*************************************************
* [SHOW] Zero-Knowledge proof parameters
**************************************************/
// Ring degree for the show proof
#define PARAM_N_SHOW 64
// Ring degree gap between the show proof and the signature (subring embedding)
#define PARAM_K_SHOW 4
// Modulus for the show proof (single prime)
#define PARAM_Q_SHOW 288230376151711717L
// Modulus bit-length upper bound for uniform sampling
#define PARAM_Q_SHOW_BITLEN 58
// Module rank for the show proof
#define PARAM_D_SHOW 23
// Witness dimension
#define PARAM_M1_SHOW 149
// ABDLOP commitment randomness dimension
#define PARAM_M2_SHOW 74
// Soundness amplification dimension
#define PARAM_L_SHOW 2
// Dimension for Approximate Range Proof
#define PARAM_ARP_SHOW 256
// Rank for Approximate Range Proof (256 / n)
#define PARAM_ARP_DIV_N_SHOW 4
// 256 / n + l
#define PARAM_ARP_DIV_N_L_SHOW 6
// Gaussian mask width for cs_1
#define PARAM_S1_SHOW 29907260.15818377584218978882
// Squared Gaussian mask width for cs_1
#define PARAM_S1SQ_SHOW 894444210169286.62500000000000000000
// Gaussian mask width for cs_2
#define PARAM_S2_SHOW 13625.45460733394611452240
// Squared Gaussian mask width for cs_2
#define PARAM_S2SQ_SHOW 185653013.25651785731315612793
// Gaussian mask width for Rs_1 (ARP)
#define PARAM_S3_SHOW 48203799.59694476425647735596
// Squared Gaussian mask width for Rs_1 (ARP)
#define PARAM_S3SQ_SHOW 2323606295582412.00000000000000000000
// Rejection sampling rate for y_1
#define PARAM_REJ1_SHOW 2
// Rejection sampling rate for y_2
#define PARAM_REJ2_SHOW 2
// Rejection sampling rate for y_3
#define PARAM_REJ3_SHOW 2
// Squared verification bound for z_1
#define PARAM_B1SQ_SHOW 6559758562648003584UL
// Squared verification bound for z_2
#define PARAM_B2SQ_SHOW 729389531316UL
// Infinity norm of challenges
#define PARAM_RHO_SHOW 8
// Manhattan-like norm of challenges
#define PARAM_ETA_SHOW 93

/*************************************************
* Legacy (to be removed in Phase 5 with poly_q_vec_k / poly_q_mat_d_k)
**************************************************/
// Old single-gadget dimension – still referenced by poly_q_vec_k.h / poly_q_mat_d_k.h
// and the not-yet-updated osig / show layers.
#define PARAM_K 5
// Old Gaussian width s_2 for poly_q_vec_k_sample_gaussian_s2 (unused, kept for compilation)
#define PARAM_S2 PARAM_S4

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
