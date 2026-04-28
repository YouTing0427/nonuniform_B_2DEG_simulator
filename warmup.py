import numpy as np
from scipy.linalg import expm
import matplotlib.pyplot as plt

dt = 0.5
f = 1
times, p0s, p1s = [], [], []
psi = np.array([1,0], dtype = complex)

def Hamiltonian(t):
    return np.array([
        [1, 1j * np.cos(f * t)],
        [-1j * np.cos(f * t), 2]
    ], dtype = complex)

def evolve(psi, H, dt):
    U = expm(-1j * H * dt)
    return U @ psi

print("f=", f, "tgrid=", dt, "\nt | psi0 | psi1 | norm")

for tgrid in range(10):
    t = tgrid * dt

    psi = evolve(psi, Hamiltonian(t), dt)
    times.append(t)
    p0s.append(abs(psi[0])**2)
    p1s.append(abs(psi[1])**2)

for tgrid in range(100):
    t = tgrid * dt
    psi = evolve(psi, Hamiltonian(t), dt)
    p0 = abs(psi[0])**2
    p1 = abs(psi[1])**2
    norm = np.linalg.norm(psi)

    print(f"{t:.2f} | {p0:.4f} | {p1:.4f} | {norm:.6f}")


plt.plot(times, p0s, label = '|psi0|^2')
plt.plot(times, p1s, label = '|psi1|^2')
plt.xlabel('time')
plt.ylabel('probability')
plt.legend()
plt.tight_layout()
plt.savefig('warmup3.png', dpi = 150)
plt.show()