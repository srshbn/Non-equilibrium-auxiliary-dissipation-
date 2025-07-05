# Non-equilibrium-auxiliary-dissipation

# Simulating Dissipative Quantum Dynamics with Auxiliary Qubits
This script simulates a chain of qubits interacting with auxiliary qubits under dissipation, using QuTiP. It models open quantum system evolution and studies how fidelity with the ground state varies with system parameters.

# Description
Constructs a spin Hamiltonian for a 1D chain with auxiliary edge couplings.
Simulates dissipation using collapse operators and periodic qubit resets.
Runs multiple dissipative cycles to evolve toward a steady state.
Sweeps over J/h ratios and computes fidelity with the ground state.
Performs 2D parameter sweep over J, h, and γ (decay rate), and visualizes fidelity heatmaps.
O# utput
A line plot showing fidelity vs J/h.
A 2D heatmap showing maximum fidelity for combinations of J and h, optimized over γ.
# Requirements
Install the required Python packages:

pip install numpy matplotlib qutip
How to Run
Simply run the script:

python dissipative_qutip_simulation.py
The script will:

Build the system and auxiliary qubit model
Simulate open-system evolution
Plot fidelity results
# Customize
Change N_system and N_aux to adjust system size.
Tune num_cycles, gamma, and other physical parameters.
Modify the Hamiltonian to simulate different qubit arrangements or interactions.
This script is useful for prototyping dissipative protocols, variational ansatz prep, and hybrid quantum dynamics models.
