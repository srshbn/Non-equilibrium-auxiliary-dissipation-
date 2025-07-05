import numpy as np
import matplotlib.pyplot as plt
from qutip import *

# Parameters
N_system = 4
N_aux = 2
N_total = N_system + N_aux
J_aux = 0.7
h = 1
gamma = 0.5
tlist = np.linspace(0, 10, 500)
num_cycles = 5
aux_indices = [N_system, N_system + 1]

# Pauli operators
si = qeye(2)
sx = sigmax()
sz = sigmaz()

# ⬇️ Function: build the Hamiltonian
def construct_hamiltonian(J, J_aux, h, N_system, N_aux):
    N_total = N_system + N_aux
    H = 0
    for i in range(N_system - 1):
        H += -J * tensor([si]*i + [sx, sx] + [si]*(N_total - i - 2))
    for i in range(N_system):
        H += -h * tensor([si]*i + [sz] + [si]*(N_total - i - 1))
    H += -J_aux * tensor([sx, sx] + [si]*(N_total - 2))  # edge 0 <-> aux
    H += -J_aux * tensor([si]*(N_system - 1) + [sx, sx] + [si]*(N_aux - 1))  # edge 3 <-> aux
    return H

# ⬇️ Function: reset an auxiliary qubit to |0>
def reset_qubit(rho, index, N_total):
    P0 = basis(2, 0) * basis(2, 0).dag()
    ops = [si for _ in range(N_total)]
    ops[index] = P0
    reset_op = tensor(ops)
    return (reset_op * rho * reset_op.dag()).unit()

# ⬇️ Function: construct collapse ops
def construct_c_ops(gamma, N_total):
    return [np.sqrt(gamma) * tensor([si]*i + [destroy(2)] + [si]*(N_total - i - 1)) for i in range(N_total)]

# ⬇️ Function: apply 1 dissipative cycle
def dissipative_cycle(rho, H, c_ops, tlist, aux_indices, N_total):
    result = mesolve(H, rho, tlist, c_ops)
    rho_final = result.states[-1]
    for i in aux_indices:
        rho_final = reset_qubit(rho_final, i, N_total)
    return rho_final

# ⬇️ Function: evolve over multiple cycles
def simulate_dissipative_evolution(H, tlist, c_ops, N_total, aux_indices, num_cycles=5):
    rho = tensor([basis(2, 0)] * N_total)  # all qubits in |0⟩
    rho = rho * rho.dag()
    for _ in range(num_cycles):
        rho = dissipative_cycle(rho, H, c_ops, tlist, aux_indices, N_total)
    return rho

# ⬇️ Scan over J/h
ratios = np.linspace(0.1, 2.0, 30)
fidelities = []

for r in ratios:
    J = r * h
    H = construct_hamiltonian(J, J_aux, h, N_system, N_aux)
    c_ops = construct_c_ops(gamma, N_system + N_aux)
    rho = simulate_dissipative_evolution(H, tlist, c_ops, N_system + N_aux, aux_indices, num_cycles)
    _, ground_state = H.groundstate()
    fidelities.append(fidelity(rho, ground_state))

# ⬇️ Plot the result
plt.figure(figsize=(8, 5))
plt.plot(ratios, fidelities, marker='o')
plt.xlabel("J/h")
plt.ylabel("Fidelity with Ground State")
plt.title("Fidelity vs J/h")
plt.grid(True)
plt.tight_layout()


# ⬇️ Sweep over J, h, and gamma
J_vals = np.linspace(0.2, 2.0, 2)
h_vals = np.linspace(0.2, 2.0, 2)
gamma_vals = np.linspace(0.1, 1.0, 2)

fidelity_grid = np.zeros((len(h_vals), len(J_vals)))

for i, h in enumerate(h_vals):
    for j, J in enumerate(J_vals):
        best_fid = 0
        for gamma in gamma_vals:
            H = construct_hamiltonian(J, J_aux, h, N_system, N_aux)
            c_ops = construct_c_ops(gamma, N_total)
            rho = simulate_dissipative_evolution(H, tlist, c_ops, N_total, aux_indices, num_cycles)
            _, ground = H.groundstate()
            f = fidelity(rho, ground)
            best_fid = max(best_fid, f)
        fidelity_grid[i, j] = best_fid

# ⬇️ Plot 2D heatmap of max fidelity over gamma
plt.figure(figsize=(8, 6))
plt.imshow(
    fidelity_grid,
    origin='lower',
    aspect='auto',
    extent=[J_vals[0], J_vals[-1], h_vals[0], h_vals[-1]],
    cmap='viridis'
)
plt.colorbar(label='Max Fidelity (over γ)')
plt.xlabel("J")
plt.ylabel("h")
plt.title("Max Fidelity vs J and h (optimized over γ)")
plt.tight_layout()

plt.show()
