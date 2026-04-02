# -*- coding: utf-8 -*-
"""
Parameter computation for integrating the tag-friendly sampler
into the lattice-anonymous-credentials scheme.

This script is self-contained — it does NOT require SageMath or the lattice
estimator.  The estimator is only needed for security level validation; all
cryptographic parameters (Gaussian widths, verification bounds, witness
dimensions) can be computed with standard Python math.  Stub functions
replace the estimator calls so the parameter classes run without Sage.

Run from any directory:
    python compute_params_tag_friendly.py
"""
import sys, os
from math import sqrt, exp, log, log2, floor, ceil, pi, e
from scipy.special import comb
from sympy import isprime, prevprime, nextprime, divisors

# ---------------------------------------------------------------------------
# Stubs for the lattice estimator (estimate_SIS_LWE functions + ND class)
# These are only used for security level printing; the crypto parameters
# themselves don't depend on them.
# ---------------------------------------------------------------------------
class _ND:
    @staticmethod
    def CenteredBinomial(k): return None
    @staticmethod
    def UniformMod(q):       return None
    @staticmethod
    def Uniform(lo, hi):     return None

ND = _ND()

def estimate_LWE(n, q, Xs, Xe, m, cost_model, rough=False):
    return (0, 0.0, 0.0, 0.0)

def estimate_SIS(n, m, q, beta, cost_model, betainf=None):
    return (0, 0.0, 0, 0.0, 0.0)

def estimate_ISIS(n, m, q, beta, cost_model):
    return (0, 0.0, 0, 0.0, 0.0)

# ---------------------------------------------------------------------------
# Paste the relevant class bodies from parameters_anonymous_credentials.py
# (with the security computation sections left intact — they just call our
#  stubs and produce zeros, which is fine).
# ---------------------------------------------------------------------------

COST_MODEL  = 'realistic_sieving'
ROUGH       = False
QUANTUM     = False
LOG2_EPS    = -40
SPECTRAL_SLACK = 6

def c_star(dim:int, sec:float):
    """Bisection: find c>0 s.t. sec + dim*(log2(c*sqrt(2pi))+(1/2-pi*c^2)*log2(e)) = 0.
    f is maximised at c_max = 1/sqrt(2pi) with f(c_max) = sec > 0.
    We return the right root (the tailcut rate), just above c_max.
    """
    f = lambda c: sec + dim * (log2(c * sqrt(2*pi)) + (1/2 - pi*c**2)*log2(e))
    a = 1.0 / sqrt(2*pi)          # c_max: f(a) = sec > 0  always
    b = a + 1.0                   # start to the right and expand until f < 0
    while f(b) > 0:
        b *= 2
    for _ in range(200):
        m = (a + b) / 2
        if f(m) > 0: a = m
        else:         b = m
    return (a + b) / 2

# ---- SIG_Parameters (verbatim from parameters_anonymous_credentials.py) ---

class SIG_Parameters:
    def __init__(self, target_bitsec:int, n:int, d:int, b_H:int, k_H:int,
                 b_L:int, k_L:int, m:int, no_guessing=False):
        self.n   = n
        self.d   = d
        self.m_s = 2 * self.d
        self.m   = 2 * self.d + m
        self.Q   = 2 ** 32

        assert isprime(b_H) and (b_H % 8 == 5) and b_H >= 2
        self.b_H = b_H;  self.k_H = k_H;  self.q_H = self.b_H ** self.k_H
        self.b_L = b_L;  self.k_L = k_L;  self.q_L = self.b_L ** self.k_L
        self.q   = self.q_L * self.q_H

        w = 1
        while comb(self.n, w) < self.Q:
            w += 1
        self.w = w

        self.smoothing_ndkLkH = sqrt((log(2*self.n*self.d*(self.k_L+self.k_H)) - LOG2_EPS*log(2))/pi)
        self.smoothing_ndkL   = sqrt((log(2*self.n*self.d*self.k_L) - LOG2_EPS*log(2))/pi)
        self.smoothing_ndkH   = sqrt((log(2*self.n*self.d*self.k_H) - LOG2_EPS*log(2))/pi)

        B_R      = 3/4 * (sqrt(self.d*self.n) + sqrt(self.k_H*self.d*self.n) + SPECTRAL_SLACK)
        norm_Si_m = sqrt(self.d*self.n/2) * sqrt(self.m*self.n)
        self.square_bound_sk = B_R**2

        self.s_L = self.smoothing_ndkL   * self.b_L
        self.s_H = self.smoothing_ndkH   * self.b_H
        self.s_1 = sqrt(self.s_L**4/(self.s_L**2 - self.smoothing_ndkLkH**2)
                        + 3*self.s_H**4/(self.s_H**2 - self.smoothing_ndkLkH**2) * B_R**2)
        self.s_2 = self.s_L**2 / sqrt(self.s_L**2 - self.smoothing_ndkLkH**2)
        self.s_3 = self.s_H**2 / sqrt(self.s_H**2 - self.smoothing_ndkLkH**2) * sqrt(3) * B_R
        self.s_4 = self.s_H**2 / sqrt(self.s_H**2 - self.smoothing_ndkLkH**2) * sqrt(3)
        self.s_11 = sqrt(
            self.s_L**4/(self.s_L**2 - self.smoothing_ndkLkH**2)
            * (self.b_L**(2*self.k_L) - 1)/(self.b_L**2 - 1)
            + 3*self.s_H**4/(self.s_H**2 - self.smoothing_ndkLkH**2) * B_R**2)

        self.alpha_11 = self.s_11 / (norm_Si_m + sqrt(self.d*self.n))
        self.alpha_12 = self.s_3  / (norm_Si_m + sqrt(self.d*self.n))
        self.M_11 = exp(pi / self.alpha_11**2)
        self.M_12 = exp(pi / self.alpha_12**2)

        # Security estimates (stubs return 0, that's fine)
        Xs = ND.CenteredBinomial(1); Xe = ND.CenteredBinomial(1)
        res = estimate_LWE(n=self.n*self.d, q=self.q, Xs=Xs, Xe=Xe,
                           m=self.n*self.d, cost_model=COST_MODEL, rough=ROUGH)
        self.mlwe_S_blocksize = floor(res[0]); self.mlwe_S_rhf = res[1]
        self.mlwe_S_coresvp_c = res[2];        self.mlwe_S_coresvp_q = res[3]

        Xs = ND.UniformMod(self.q); Xe = ND.Uniform(0, 1)
        res = estimate_LWE(n=self.n*self.d, q=self.q, Xs=Xs, Xe=Xe,
                           m=2*self.n*self.d, cost_model=COST_MODEL, rough=ROUGH)
        self.mlwe_U_blocksize = floor(res[0]); self.mlwe_U_rhf = res[1]
        self.mlwe_U_coresvp_c = res[2];        self.mlwe_U_coresvp_q = res[3]

        tail_cut_prob_log = 20
        self.B_11_s = floor(c_star(self.d*self.n, tail_cut_prob_log)**2
                            * self.s_11**2 * (self.d*self.n))
        self.B_11   = sqrt(self.B_11_s)
        self.B_12_s = floor(c_star(self.d*self.n, tail_cut_prob_log)**2
                            * self.s_3**2  * (self.d*self.n))
        self.B_12   = sqrt(self.B_12_s)
        self.B_2_s  = floor(c_star(self.k_H*self.d*self.n, tail_cut_prob_log)**2
                            * self.s_4**2  * (self.k_H*self.d*self.n))
        self.B_2    = sqrt(self.B_2_s)
        self.B_3_s  = floor(c_star(self.n, tail_cut_prob_log)**2
                            * self.s_4**2  * self.n)
        self.B_3    = sqrt(self.B_3_s)

        self.B_11_prime_s = floor((self.B_11 + sqrt(self.d*self.n))**2)
        self.B_11_prime   = sqrt(self.B_11_prime_s)
        self.B_12_prime_s = floor((self.B_12 + sqrt(self.d*self.n))**2)
        self.B_12_prime   = sqrt(self.B_12_prime_s)

        # M-ISIS and M-SIS stubs
        self.misis_beta = sqrt(floor(2*self.n*self.d))
        res = estimate_ISIS(n=self.n*self.d, m=2*self.n*self.d, q=self.q,
                            beta=self.misis_beta, cost_model=COST_MODEL)
        self.misis_subdim=res[0]; self.misis_rhf=res[1]; self.misis_bkz=res[2]
        self.misis_coresvp_c=res[3]; self.misis_coresvp_q=res[4]

        self.msis_beta_I = sqrt(floor(
            (sqrt(self.B_11_prime_s + self.B_12_prime_s) + sqrt(self.d*self.n)*self.B_2)**2
            + self.B_3_s + self.n + 1))
        res = estimate_SIS(n=self.n*self.d, m=self.n*(2*self.d+2+self.m),
                           q=self.q, beta=self.msis_beta_I, cost_model=COST_MODEL)
        self.msis_subdim_I=res[0]; self.msis_rhf_I=res[1]; self.msis_bkz_I=res[2]
        self.msis_coresvp_c_I=res[3]; self.msis_coresvp_q_I=res[4]

        self.msis_beta_II = sqrt(floor(
            (2*sqrt(self.B_11_prime_s + self.B_12_prime_s)
             + sqrt(self.d*self.n)*sqrt(4*self.B_2_s + self.n))**2
            + 4*self.B_3_s))
        res = estimate_SIS(n=self.n*self.d, m=self.n*(2*self.d+1),
                           q=self.q, beta=self.msis_beta_II, cost_model=COST_MODEL)
        self.msis_subdim_II=res[0]; self.msis_rhf_II=res[1]; self.msis_bkz_II=res[2]
        self.msis_coresvp_c_II=res[3]; self.msis_coresvp_q_II=res[4]

        self.tag_space_size   = comb(self.n, self.w)
        self.forgery_advantage = 0.0
        self.security = 0.0

        self.pk_bitsize  = self.d*self.k_H*self.d*self.n*ceil(log2(self.q)) + 256
        self.sk_bitsize  = 2*self.d*k_H*self.d*self.n*ceil(log2(3))
        self.tag_bitsize = self.n
        self.v_12_bitsize  = ceil(self.n*self.d*(1/2 + log2(self.s_3)))
        self.v_23_bitsize  = ceil((self.k_H*self.d + 1)*self.n*(1/2 + log2(self.s_4)))
        self.v_2_bitsize   = ceil(self.k_H*self.d*self.n*(1/2 + log2(self.s_4)))
        self.v_3_bitsize   = ceil(self.n*(1/2 + log2(self.s_4)))
        self.sig_bitsize   = self.tag_bitsize + self.v_12_bitsize + self.v_23_bitsize


# ---- Issue_ZKP_Parameters (verbatim) -------------------------------------

class Issue_ZKP_Parameters:
    def __init__(self, target_bitsec:int, n:int, d:int, m_2:int, q_start:int,
                 n_attr:int, gamma:int, D:int, sig,
                 bimodal=True, compression=True, garbage=True):
        self.sec    = target_bitsec
        self.n      = n
        self.d      = d
        self.n_attr = n_attr
        self.kappa  = 2

        if compression:
            found_q_gamma = False
            q = q_start + 1
            while not found_q_gamma:
                q = prevprime(q)
                while q % 8 != 5:
                    q = prevprime(q)
                divs = divisors(q - 1)
                for div in divs:
                    if (gamma <= div <= 5*gamma/4) and (div % 2 == 0):
                        self.gamma = div; found_q_gamma = True; break
            self.D = D
        else:
            self.gamma = 0; self.D = 0
            q = prevprime(q_start + 1)
            while (q % (4*self.kappa) != 2*self.kappa + 1) or (q < (2*sqrt(self.kappa))**self.kappa):
                q = prevprime(q)
        self.q     = q
        self.q_min = self.q

        if garbage:
            l = ceil(self.sec / log2(self.q_min)); self.l = ceil(l/2)
        else:
            self.l = ceil(self.sec / log2(self.q_min))

        self.rho = ceil(1/2 * (2**(2*(self.sec+1)/self.n) - 1))
        self.eta = {64:93, 128:42, 256:37, 512:57, 1024:84}[self.n]
        self.challenge_space_size = (2*self.rho + 1)**(self.n // 2) / 2
        self.k = sig.n // self.n

        self.m_1  = 2 * sig.d * self.k
        self.m_1 += (sig.m - self.n_attr) * self.k
        if bimodal: self.m_1 += 1

        self.m_2    = m_2
        self.xi_s2  = 1
        self.bound_witness = sqrt(self.n * self.m_1)
        self.B_j = 1/sig.q * (sig.q/2 * sqrt(sig.n*sig.d) +
                               sqrt(sig.n*sig.d) +
                               sig.q/2 * sig.n*sig.d*sqrt(sig.n*sig.d) +
                               sig.q/2 * sig.n*sig.m*sqrt(sig.n*sig.d))

        self.M_1 = 2; self.M_2 = 2; self.M_3 = 2
        if bimodal:
            self.eps_1=0; self.eps_2=0; self.eps_3=0
            self.alpha_1 = sqrt(pi/log(self.M_1))
            self.alpha_2 = sqrt(pi/log(self.M_2))
            self.alpha_3 = sqrt(pi/log(self.M_3))
        else:
            self.eps_1=2**-130; self.eps_2=2**-130; self.eps_3=2**-130
            self.alpha_1 = sqrt(pi)/log(self.M_1)*(sqrt(log(1/self.eps_1)+log(self.M_1))+sqrt(log(1/self.eps_1)))
            self.alpha_2 = sqrt(pi)/log(self.M_2)*(sqrt(log(1/self.eps_2)+log(self.M_2))+sqrt(log(1/self.eps_2)))
            self.alpha_3 = sqrt(pi)/log(self.M_3)*(sqrt(log(1/self.eps_3)+log(self.M_3))+sqrt(log(1/self.eps_3)))

        self.s_1 = self.alpha_1 * self.eta * self.bound_witness
        self.s_2 = self.alpha_2 * self.eta * self.xi_s2 * sqrt(self.n * self.m_2)
        self.s_3 = self.alpha_3 * sqrt(337) * sqrt(self.bound_witness**2 + self.B_j**2)

        self.B_256_s = floor(c_star(256, self.sec+3)**2 * self.s_3**2 * 256)
        self.B_256   = sqrt(self.B_256_s)
        self.B_arp   = self.B_256 * 2/sqrt(26)
        self.cond_1_bound = 41 * self.n * self.m_1 * self.B_arp
        self.cond_2_bound = self.B_arp**2 + self.B_arp * sqrt(self.n*self.m_1)
        self.cond_3_bound = self.B_arp**2 - 1
        self.cond_4_bound = (sig.q/2 + 1 + sig.q/2*sig.n*sig.d +
                             sig.q/2*sig.n*sig.m + sig.q*self.B_arp)

        assert (self.q > self.cond_1_bound and self.q > self.cond_2_bound and
                self.q > self.cond_3_bound and self.q > self.cond_4_bound), \
               (f"ISS ZKP modulus too small (q={self.q}).\n"
                f"  cond_1={self.cond_1_bound:.3e}\n  cond_2={self.cond_2_bound:.3e}\n"
                f"  cond_3={self.cond_3_bound:.3e}\n  cond_4={self.cond_4_bound:.3e}")

        self.B_1_s = 4 * floor(c_star(self.m_1*self.n, self.sec+3)**2 * self.s_1**2 * (self.m_1*self.n))
        if self.D != 0 and self.gamma != 0:
            self.B_2_s = floor((2*sqrt(floor(c_star(self.m_2*self.n, self.sec+3)**2 * self.s_2**2 * (self.m_2*self.n)))
                                + (2**self.D*self.eta + self.gamma)*sqrt(self.n*self.d))**2)
        else:
            self.B_2_s = 4 * floor(c_star(self.m_2*self.n, self.sec+3)**2 * self.s_2**2 * (self.m_2*self.n))

        # Security stubs
        Xs = ND.CenteredBinomial(self.xi_s2); Xe = ND.CenteredBinomial(self.xi_s2)
        res = estimate_LWE(n=self.n*(self.m_2-(self.d+floor(256/self.n)+self.l+1)),
                           q=self.q, Xs=Xs, Xe=Xe, m=self.n*self.m_2,
                           cost_model=COST_MODEL, rough=ROUGH)
        self.mlwe_blocksize=floor(res[0]); self.mlwe_rhf=res[1]
        self.mlwe_coresvp_c=res[2];        self.mlwe_coresvp_q=res[3]

        self.msis_beta = 4 * self.eta * sqrt(self.B_1_s + self.B_2_s)
        res = estimate_SIS(n=self.n*self.d, m=self.n*(self.m_1+self.m_2),
                           q=self.q, beta=self.msis_beta, cost_model=COST_MODEL)
        self.msis_subdim=res[0]; self.msis_rhf=res[1]; self.msis_bkz=res[2]
        self.msis_coresvp_c=res[3]; self.msis_coresvp_q=res[4]

        self.soundness_error = (self.q_min**(-self.n/self.kappa) +
                                2/self.challenge_space_size + 0.0)
        if garbage: self.soundness_error += self.q_min**(-2*self.l)
        else:       self.soundness_error += self.q_min**(-self.l)


# ---- Show_ZKP_Parameters (verbatim) --------------------------------------

class Show_ZKP_Parameters:
    def __init__(self, target_bitsec:int, n:int, d:int, m_2:int, q_start:int,
                 n_attr:int, gamma:int, D:int, sig,
                 bimodal=True, compression=True, garbage=True):
        self.sec    = target_bitsec
        self.n      = n
        self.d      = d
        self.n_attr = n_attr
        self.kappa  = 2

        if compression:
            found_q_gamma = False
            q = q_start + 1
            while not found_q_gamma:
                q = prevprime(q)
                while q % 8 != 5:
                    q = prevprime(q)
                divs = divisors(q - 1)
                for div in divs:
                    if (gamma <= div <= 5*gamma/4) and (div % 2 == 0):
                        self.gamma = div; found_q_gamma = True; break
            self.D = D
        else:
            self.gamma = 0; self.D = 0
            q = prevprime(q_start + 1)
            while (q % (4*self.kappa) != 2*self.kappa + 1) or (q < (2*sqrt(self.kappa))**self.kappa):
                q = prevprime(q)
        self.q     = q
        self.q_min = q

        if garbage:
            l = ceil(self.sec / log2(self.q_min)); self.l = ceil(l/2)
        else:
            self.l = ceil(self.sec / log2(self.q_min))

        self.rho = ceil(1/2 * (2**(2*(self.sec+1)/self.n) - 1))
        self.eta = {64:93, 128:42, 256:37, 512:57, 1024:84}[self.n]
        self.challenge_space_size = (2*self.rho + 1)**(self.n // 2) / 2
        self.k = sig.n // self.n

        self.m_1  = sig.d * self.k + 1          # v_11 + a_11
        self.m_1 += sig.d * self.k + 1          # v_12 + a_12
        self.m_1 += sig.d * sig.k_H * self.k + 1  # v_2 + a_2
        self.m_1 += self.k + 1                  # v_3 + a_3
        self.m_1 += self.k                      # tag
        self.m_1 += (sig.m - self.n_attr) * self.k  # msg
        if bimodal: self.m_1 += 1

        self.m_2   = m_2
        self.xi_s2 = 1
        self.bound_witness = sqrt(sig.B_11_prime_s + sig.B_12_prime_s +
                                  sig.B_2_s + sig.B_3_s + sig.w +
                                  sig.n*(sig.m - self.n_attr) + (1 if bimodal else 0))
        self.B_j = 1/sig.q * (sig.q/2 * sqrt(sig.n*sig.d) +
                               sig.q/2 * sig.n*sig.m*sqrt(sig.n*sig.d) +
                               sig.B_11_prime +
                               sig.q/2 * sqrt(sig.n*sig.d) * sig.B_12_prime +
                               sig.q/2 * sqrt(sig.n*sig.d*sig.k_H) * sig.B_2 +
                               sig.q/2 * sqrt(sig.n) * sig.B_3 +
                               sig.q_L*(sig.q_H-1)/(sig.b_H-1)*sig.w*sig.B_2*sqrt(sig.n*sig.d))

        self.M_1 = 2; self.M_2 = 2; self.M_3 = 2
        if bimodal:
            self.eps_1=0; self.eps_2=0; self.eps_3=0
            self.alpha_1 = sqrt(pi/log(self.M_1))
            self.alpha_2 = sqrt(pi/log(self.M_2))
            self.alpha_3 = sqrt(pi/log(self.M_3))
        else:
            self.eps_1=2**-130; self.eps_2=2**-130; self.eps_3=2**-130
            self.alpha_1 = sqrt(pi)/log(self.M_1)*(sqrt(log(1/self.eps_1)+log(self.M_1))+sqrt(log(1/self.eps_1)))
            self.alpha_2 = sqrt(pi)/log(self.M_2)*(sqrt(log(1/self.eps_2)+log(self.M_2))+sqrt(log(1/self.eps_2)))
            self.alpha_3 = sqrt(pi)/log(self.M_3)*(sqrt(log(1/self.eps_3)+log(self.M_3))+sqrt(log(1/self.eps_3)))

        self.s_1 = self.alpha_1 * self.eta * self.bound_witness
        self.s_2 = self.alpha_2 * self.eta * self.xi_s2 * sqrt(self.n * self.m_2)
        self.s_3 = self.alpha_3 * sqrt(337) * sqrt(self.bound_witness**2 + self.B_j**2)

        self.B_256_s = floor(c_star(256, self.sec+3)**2 * self.s_3**2 * 256)
        self.B_256   = sqrt(self.B_256_s)
        self.B_arp   = self.B_256 * 2/sqrt(26)
        self.cond_1_bound = max(sig.B_11_prime_s, sig.B_12_prime_s, sig.B_2_s, sig.B_3_s, sig.w)
        self.cond_2_bound = 41 * self.n * self.m_1 * self.B_arp
        self.cond_3_bound = self.B_arp**2 - min(sig.B_11_prime_s, sig.B_12_prime_s,
                                                 sig.B_2_s, sig.B_3_s, sig.w)
        self.cond_4_bound = sig.w + sqrt(sig.w*self.n*self.k)
        self.cond_5_bound = self.B_arp**2 + self.B_arp*sqrt(self.n*self.k*(sig.m - self.n_attr))
        self.cond_6_bound = self.B_arp**2 - 1
        self.cond_7_bound = (sig.q/2 + sig.q/2*sig.n*sig.m +
                             sig.B_11_prime +
                             sig.q/2*sqrt(sig.n*sig.d)*sig.B_12_prime +
                             sig.q/2*sqrt(sig.n*sig.d*sig.k_H)*sig.B_2 +
                             sig.q/2*sqrt(sig.n)*sig.B_3 +
                             sig.q_L*(sig.q_H-1)/(sig.b_H-1)*sig.w*sig.B_2 +
                             sig.q*self.B_arp)

        assert all(self.q > c for c in [self.cond_1_bound, self.cond_2_bound,
                                         self.cond_3_bound, self.cond_4_bound,
                                         self.cond_5_bound, self.cond_6_bound,
                                         self.cond_7_bound]), \
               (f"SHOW ZKP modulus too small (q={self.q}).\n"
                f"  cond_1={self.cond_1_bound:.3e}\n  cond_2={self.cond_2_bound:.3e}\n"
                f"  cond_3={self.cond_3_bound:.3e}\n  cond_4={self.cond_4_bound:.3e}\n"
                f"  cond_5={self.cond_5_bound:.3e}\n  cond_6={self.cond_6_bound:.3e}\n"
                f"  cond_7={self.cond_7_bound:.3e}")

        self.B_1_s = 4 * floor(c_star(self.m_1*self.n, self.sec+3)**2 * self.s_1**2 * (self.m_1*self.n))
        if self.D != 0 and self.gamma != 0:
            self.B_2_s = floor((2*sqrt(floor(c_star(self.m_2*self.n, self.sec+3)**2 * self.s_2**2 * (self.m_2*self.n)))
                                + (2**self.D*self.eta + self.gamma)*sqrt(self.n*self.d))**2)
        else:
            self.B_2_s = 4 * floor(c_star(self.m_2*self.n, self.sec+3)**2 * self.s_2**2 * (self.m_2*self.n))

        # Security stubs
        Xs = ND.CenteredBinomial(self.xi_s2); Xe = ND.CenteredBinomial(self.xi_s2)
        res = estimate_LWE(n=self.n*(self.m_2-(self.d+floor(256/self.n)+self.l+1)),
                           q=self.q, Xs=Xs, Xe=Xe, m=self.n*self.m_2,
                           cost_model=COST_MODEL, rough=ROUGH)
        self.mlwe_blocksize=floor(res[0]); self.mlwe_rhf=res[1]
        self.mlwe_coresvp_c=res[2];        self.mlwe_coresvp_q=res[3]

        self.msis_beta = 4 * self.eta * sqrt(self.B_1_s + self.B_2_s)
        res = estimate_SIS(n=self.n*self.d, m=self.n*(self.m_1+self.m_2),
                           q=self.q, beta=self.msis_beta, cost_model=COST_MODEL)
        self.msis_subdim=res[0]; self.msis_rhf=res[1]; self.msis_bkz=res[2]
        self.msis_coresvp_c=res[3]; self.msis_coresvp_q=res[4]

        self.soundness_error = (self.q_min**(-self.n/2) +
                                2/self.challenge_space_size + 0.0)
        if garbage: self.soundness_error += self.q_min**(-2*self.l)
        else:       self.soundness_error += self.q_min**(-self.l)


# ===========================================================================
# Instantiate and print
# ===========================================================================

print("Computing signature parameters...")
sig = SIG_Parameters(
    target_bitsec = 128,
    n    = 256,
    d    = 4,
    b_H  = 13,
    k_H  = 2,
    b_L  = 49,
    k_L  = 2,
    m    = 10,
)

print(f"\n{'='*70}")
print("SIGNATURE PARAMETERS")
print(f"{'='*70}")
print(f"  q       = {sig.q}")
print(f"  q_L     = {sig.q_L}")
print(f"  q_H     = {sig.q_H}")
print(f"  b_L     = {sig.b_L},  k_L = {sig.k_L}")
print(f"  b_H     = {sig.b_H},  k_H = {sig.k_H}")
print(f"  w       = {sig.w}")
print(f"  sq. spectral norm bound = {sig.square_bound_sk:.20f}")
print(f"  s_L     = {sig.s_L:.20f}")
print(f"  s_H     = {sig.s_H:.20f}")
print(f"  s_1     = {sig.s_1:.20f}  (perturbation width 1)")
print(f"  s_2     = {sig.s_2:.20f}  (perturbation width 2)")
print(f"  s_3     = {sig.s_3:.20f}  (width for v12)")
print(f"  s_4     = {sig.s_4:.20f}  (width for v2 and v3)")
print(f"  s_11    = {sig.s_11:.20f} (final width for v11)")
print(f"  s_L/b_L = {sig.s_L/sig.b_L:.20f}")
print(f"  s_H/b_H = {sig.s_H/sig.b_H:.20f}")
print(f"  s_1^2 - s_L^2          = {sig.s_1**2 - sig.s_L**2:.20f}")
print(f"  s_3^2                  = {sig.s_3**2:.20f}")
print(f"  sqrt(s_2^2 - s_L^2)   = {sqrt(sig.s_2**2 - sig.s_L**2):.20f}")
print(f"  sqrt(s_4^2 - s_H^2)   = {sqrt(sig.s_4**2 - sig.s_H**2):.20f}")
print(f"  -s_H^2*s_4^2/(s_4^2-s_H^2) = {-sig.s_H**2*sig.s_4**2/(sig.s_4**2-sig.s_H**2):.20f}")
print(f"  -s_H^2/(s_4^2-s_H^2)       = {-sig.s_H**2/(sig.s_4**2-sig.s_H**2):.20f}")
print(f"  B_11_s  = {sig.B_11_s}")
print(f"  B_12_s  = {sig.B_12_s}")
print(f"  B_2_s   = {sig.B_2_s}")
print(f"  B_3_s   = {sig.B_3_s}")
print(f"  B_11_prime_s = {sig.B_11_prime_s}")
print(f"  B_12_prime_s = {sig.B_12_prime_s}")

print("\nComputing ISS parameters (compression=False, bimodal=True, garbage=True)...")
# cond_2 = B_arp^2 ~ (1.45e7)^2 ~ 2.1e14.  Need q > 2.1e14.
# Use q_start = 2^49 ≈ 5.6e14 for safety.
issue = Issue_ZKP_Parameters(
    target_bitsec = 128,
    n       = 64,
    d       = 20,
    m_2     = 58,
    q_start = ceil(2**49),
    n_attr  = 0,
    gamma   = 0,
    D       = 0,
    sig     = sig,
    bimodal     = True,
    compression = False,
    garbage     = True,
)

print(f"\n{'='*70}")
print("ISSUANCE ZKP PARAMETERS")
print(f"{'='*70}")
print(f"  n    = {issue.n},  k = {issue.k},  d = {issue.d}")
print(f"  q    = {issue.q}  ({ceil(log2(issue.q))} bits)")
print(f"  l    = {issue.l},  rho = {issue.rho},  eta = {issue.eta}")
print(f"  m_1  = {issue.m_1},  m_2 = {issue.m_2}")
print(f"  ARP_DIV_N = {256//issue.n},  ARP_DIV_N_L = {256//issue.n + issue.l}")
print(f"  bound_witness = {issue.bound_witness:.6f}")
print(f"  B_j           = {issue.B_j:.6e}")
print(f"  B_arp         = {issue.B_arp:.6e}")
print(f"  cond_4        = {issue.cond_4_bound:.6e}  (must be < q)")
print(f"  s_1  = {issue.s_1:.20f}")
print(f"  s_2  = {issue.s_2:.20f}")
print(f"  s_3  = {issue.s_3:.20f}")
print(f"  B_1_s = {issue.B_1_s}")
print(f"  B_2_s = {issue.B_2_s}")

print("\nComputing SHOW parameters (compression=False, bimodal=True, garbage=True)...")
# B_arp ~ 2e8, so dominant condition is B_arp^2 ~ 4e16 ~ 2^55.
# Also cond_7 ~ 8e13. Use q_start = 2^58 ≈ 2.9e17, keeping q within 64 bits.
show = Show_ZKP_Parameters(
    target_bitsec = 128,
    n       = 64,
    d       = 23,
    m_2     = 74,
    q_start = ceil(2**58),
    n_attr  = 0,
    gamma   = 0,
    D       = 0,
    sig     = sig,
    bimodal     = True,
    compression = False,
    garbage     = True,
)

print(f"\n{'='*70}")
print("SHOW ZKP PARAMETERS")
print(f"{'='*70}")
print(f"  n    = {show.n},  k = {show.k},  d = {show.d}")
print(f"  q    = {show.q}  ({ceil(log2(show.q))} bits)")
print(f"  l    = {show.l},  rho = {show.rho},  eta = {show.eta}")
print(f"  m_1  = {show.m_1},  m_2 = {show.m_2}")
print(f"  ARP_DIV_N = {256//show.n},  ARP_DIV_N_L = {256//show.n + show.l}")
print(f"  bound_witness = {show.bound_witness:.6f}")
print(f"  B_j           = {show.B_j:.6e}")
print(f"  B_arp         = {show.B_arp:.6e}")
print(f"  cond_7        = {show.cond_7_bound:.6e}  (must be < q)")
print(f"  s_1  = {show.s_1:.20f}")
print(f"  s_2  = {show.s_2:.20f}")
print(f"  s_3  = {show.s_3:.20f}")
print(f"  B_1_s = {show.B_1_s}")
print(f"  B_2_s = {show.B_2_s}")

# ===========================================================================
# Final summary: generate the new params.h block
# ===========================================================================
print(f"\n\n{'#'*70}")
print("# NEW params.h CONSTANTS")
print(f"{'#'*70}\n")

print("/* Signature parameters */")
print(f"#define PARAM_Q          {sig.q}L")
print(f"#define PARAM_QL         {sig.q_L}L")
print(f"#define PARAM_QH         {sig.q_H}L")
print(f"#define PARAM_Q_BITLEN   {ceil(log2(sig.q))}")
print(f"#define PARAM_D          {sig.d}")
print(f"#define PARAM_KH         {sig.k_H}")
print(f"#define PARAM_KL         {sig.k_L}")
print(f"#define PARAM_BH         {sig.b_H}")
print(f"#define PARAM_BL         {sig.b_L}")
print(f"#define PARAM_M          10")
print(f"#define PARAM_W          {sig.w}")
print(f"#define PARAM_R_MAX_SQ_SPECTRAL_NORM {sig.square_bound_sk:.20f}")
print(f"#define PARAM_S1SQ_SLSQ  {sig.s_1**2 - sig.s_L**2:.20f}")
print(f"#define PARAM_S3SQ       {sig.s_3**2:.20f}")
print(f"#define PARAM_S4         {sig.s_4:.20f}")
print(f"#define PARAM_SL_DIV_BL  {sig.s_L/sig.b_L:.20f}")
print(f"#define PARAM_SH_DIV_BH  {sig.s_H/sig.b_H:.20f}")
print(f"#define PARAM_SQRT_S2SQ_SLSQ {sqrt(sig.s_2**2 - sig.s_L**2):.20f}")
print(f"#define PARAM_SQRT_S4SQ_SHSQ {sqrt(sig.s_4**2 - sig.s_H**2):.20f}")
print(f"#define PARAM_NEG_SHSQS4SQ_DIV_S4SQ_SHSQ {-sig.s_H**2*sig.s_4**2/(sig.s_4**2-sig.s_H**2):.20f}")
print(f"#define PARAM_NEG_SHSQ_DIV_S4SQ_SHSQ {-sig.s_H**2/(sig.s_4**2-sig.s_H**2):.20f}")
print(f"#define PARAM_B11SQ      {sig.B_11_s}L")
print(f"#define PARAM_B12SQ      {sig.B_12_s}L")
print(f"#define PARAM_B2SQ       {sig.B_2_s}L")
print(f"#define PARAM_B3SQ       {sig.B_3_s}L")
print()
print("/* ISS proof parameters */")
print(f"#define PARAM_N_ISS              {issue.n}")
print(f"#define PARAM_K_ISS              {issue.k}")
print(f"#define PARAM_Q_ISS              {issue.q}L")
print(f"#define PARAM_Q_ISS_BITLEN       {ceil(log2(issue.q))}")
print(f"#define PARAM_D_ISS              {issue.d}")
print(f"#define PARAM_M1_ISS             {issue.m_1}")
print(f"#define PARAM_M1_K_ISS           {issue.m_1 // issue.k}")
print(f"#define PARAM_M2_ISS             {issue.m_2}")
print(f"#define PARAM_L_ISS              {issue.l}")
print(f"#define PARAM_ARP_ISS            256")
print(f"#define PARAM_ARP_DIV_N_ISS      {256 // issue.n}")
print(f"#define PARAM_ARP_DIV_N_L_ISS    {256 // issue.n + issue.l}")
print(f"#define PARAM_S1_ISS             {issue.s_1:.20f}")
print(f"#define PARAM_S1SQ_ISS           {issue.s_1**2:.20f}")
print(f"#define PARAM_S2_ISS             {issue.s_2:.20f}")
print(f"#define PARAM_S2SQ_ISS           {issue.s_2**2:.20f}")
print(f"#define PARAM_S3_ISS             {issue.s_3:.20f}")
print(f"#define PARAM_S3SQ_ISS           {issue.s_3**2:.20f}")
print(f"#define PARAM_REJ1_ISS           {int(issue.M_1)}")
print(f"#define PARAM_REJ2_ISS           {int(issue.M_2)}")
print(f"#define PARAM_REJ3_ISS           {int(issue.M_3)}")
print(f"#define PARAM_B1SQ_ISS           {issue.B_1_s}UL")
print(f"#define PARAM_B2SQ_ISS           {issue.B_2_s}UL")
print(f"#define PARAM_RHO_ISS            {issue.rho}")
print(f"#define PARAM_ETA_ISS            {issue.eta}")
print()
print("/* SHOW proof parameters */")
print(f"#define PARAM_N_SHOW             {show.n}")
print(f"#define PARAM_K_SHOW             {show.k}")
print(f"#define PARAM_Q_SHOW             {show.q}L")
print(f"#define PARAM_Q_SHOW_BITLEN      {ceil(log2(show.q))}")
print(f"#define PARAM_D_SHOW             {show.d}")
print(f"#define PARAM_M1_SHOW            {show.m_1}")
print(f"#define PARAM_M2_SHOW            {show.m_2}")
print(f"#define PARAM_L_SHOW             {show.l}")
print(f"#define PARAM_ARP_SHOW           256")
print(f"#define PARAM_ARP_DIV_N_SHOW     {256 // show.n}")
print(f"#define PARAM_ARP_DIV_N_L_SHOW   {256 // show.n + show.l}")
print(f"#define PARAM_S1_SHOW            {show.s_1:.20f}")
print(f"#define PARAM_S1SQ_SHOW          {show.s_1**2:.20f}")
print(f"#define PARAM_S2_SHOW            {show.s_2:.20f}")
print(f"#define PARAM_S2SQ_SHOW          {show.s_2**2:.20f}")
print(f"#define PARAM_S3_SHOW            {show.s_3:.20f}")
print(f"#define PARAM_S3SQ_SHOW          {show.s_3**2:.20f}")
print(f"#define PARAM_REJ1_SHOW          {int(show.M_1)}")
print(f"#define PARAM_REJ2_SHOW          {int(show.M_2)}")
print(f"#define PARAM_REJ3_SHOW          {int(show.M_3)}")
print(f"#define PARAM_B1SQ_SHOW          {show.B_1_s}UL")
print(f"#define PARAM_B2SQ_SHOW          {show.B_2_s}UL")
print(f"#define PARAM_RHO_SHOW           {show.rho}")
print(f"#define PARAM_ETA_SHOW           {show.eta}")
