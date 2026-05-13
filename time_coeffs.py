import numpy as np
from scipy.linalg import expm, eigh
from params import hbar, B0, e, L, m_star, w_c, AU_TO_MEV, AU_TO_NM
from hamiltonian import build_H
from wave_function import QHO_basis

"""
Ve = ne * (hbar*B0/m_star), for ne=0.5 -> 1st flat band
Ve -> Ve*cos(w_ac*t)
w_ac: ac frequency
time_operator: exp(i*H*dt)
"""

w_ac = 1 * w_c
T = 20 * np.pi / w_ac    ## consider 10 full period
nstep = 1000
dt = T / nstep
nmax = 201

ky = -0.15 * AU_TO_NM            ## tunable
Ve = 0.5 * (hbar*B0/m_star)      ## tunable
level = 1                        ## tunable: 0 for ground state

def build_H_ac(ky, Ve, t):
    Ve_ac = Ve * np.cos(w_ac * t)      ## tunable: sin or cos
    return build_H(ky, Ve_ac, nmax)

def time_evolute(ky, Ve, level):
    """
    evolution of ground state coefficient
    return:
        1. c_snapshots(coefficient of "level", shape: (nstep, nmax))
        2. t_array(shape: nstep)
    """
    H = build_H_ac(ky, Ve, 0)
    _, evecs0 = eigh(H)
    coeffs = evecs0[:, level].astype(complex)     ## initial state coefficient

    c_snapshots = []
    t_array = []

    for step in range(nstep):
        t = step * dt
        H_t = build_H_ac(ky, Ve, t + dt/2)        ## mid of each dt
        U = expm (-1j * H_t * dt / hbar)
        coeffs = U @ coeffs
        c_snapshots.append(coeffs.copy())
        t_array.append(t)

    return np.array(c_snapshots), np.array(t_array)

if __name__ == '__main__':
    import time
    print(f"Evolving ky = -0.15 for 1/4 period ({nstep} steps)")
    
    ## runtime
    t0 = time.time()
    c_snaps, t_arr = time_evolute(ky, Ve, level)
    print(f"Done in {time.time()-t0:.1f}s. Shape: {c_snaps.shape}")

    ## check norm (expect ~1.0 throughout)
    norm = np.array([np.linalg.norm(coeffs) for coeffs in c_snaps])
    print(f"Norm: min={norm.min():.6f}, max={norm.max():.6f}")


def compute_prob_evolute(ky, Ve, level):
    """
    Returns:
        prob  : shape ()
        t_arr : shape (nstep)
        x     : x grid in atomic units
    """
    x = np.linspace(-2**0.5*L, 2**0.5*L, 401)        ## x in au
    c_snaps, t_arr = time_evolute(ky, Ve, level)
    all_basis = QHO_basis(x, nmax)
    ## shape: c_snaps=(nstep, nmax), all_basis=(nmax, lan(x)), psi_t=(nstep, len(x))
    psi_t = c_snaps @ all_basis
    prob = np.abs(psi_t)**2

    return prob, t_arr, x

if __name__ == '__main__':
    ## quick check: prob should stay normalized over time
    prob, t_arr, x = compute_prob_evolute(ky, Ve, level)
    dx = x[1] - x[0]
    norms = prob.sum(axis=1) * dx
    print(f"Spatial norm: min={norms.min():.4f}, max={norms.max():.4f}")  ## expect ~1.0