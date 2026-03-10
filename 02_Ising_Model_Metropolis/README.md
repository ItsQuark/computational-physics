# 2D Ising Model: Phase Transitions and Metropolis Algorithm

This repository contains an implementation of the 2D Ising Model, a mathematical model of ferromagnetism in statistical mechanics. The project explores the emergence of macroscopic magnetic order from microscopic spin interactions using Monte Carlo simulations.

## Physical Model

The model consists of a square lattice of size $N \times N$ where each site $(i, j)$ contains a discrete variable $s_{i,j} \in \{+1, -1\}$, representing the magnetic dipole moment of an atomic "spin". 

The energy of a given configuration $S$ is defined by the Hamiltonian:

$$E(S) = -J \sum_{\langle i,j \rangle} s_i s_j$$

where $\langle i,j \rangle$ indicates a sum over nearest neighbors. For this simulation, we assume ferromagnetic coupling ($J > 1$) and periodic boundary conditions to simulate an infinite lattice:
- $s_{0,j} = s_{N,j}$
- $s_{N+1,j} = s_{1,j}$

At a given temperature $T$, the probability of finding the system in configuration $S$ follows the Boltzmann distribution:
$P(S) \propto e^{-\beta E(S)}$, with $\beta = 1/k_B T$.

## Numerical Method: Monte Carlo Metropolis

Due to the exponential growth of the state space ($2^{N^2}$), analytical computation of the partition function $Z$ is intractable for large $N$. This project utilizes the **Metropolis-Hastings Algorithm** to generate a Markov Chain of configurations that sample the equilibrium distribution.

### The Metropolis Algorithm
1. **Initialization:** Start with a random (high $T$) or ordered (low $T$) spin configuration.
2. **Selection:** Select a random spin $s_{n,m}$ in the lattice.
3. **Energy Calculation:** Compute the energy change $\Delta E$ if the spin were to be flipped:
   $$\Delta E = 2s_{n,m} \sum_{neighbors} s_{neighbor}$$
4. **Acceptance Criterion:** Generate a uniform random number $\eta \in [0, 1]$. Accept the flip if:
   $$\eta < \min(1, e^{-\beta \Delta E})$$
5. **Iteration:** Repeat the process. 

A "Monte Carlo Step" (MCS) is defined as $N^2$ trial flips, ensuring that, on average, every spin has had the opportunity to change state once.

## Detailed Balance and Master Equation
The algorithm is constructed to satisfy the **detailed balance condition**:
$$P(S) T(S \to S') = P(S') T(S' \to S)$$
where $T$ is the transition probability. This ensures that the Boltzmann distribution is the stationary solution of the Master Equation governing the Markov process, allowing the system to relax toward thermodynamic equilibrium.

## Technical Specifications
- **Language:** Fortran 90/95
- **Visualization:** Gnuplot (Real-time monitoring of spin domains)
- **Key Parameters:** $N \times N$ lattice size, temperature $T \in [0, 5]$, and $10^6$ Monte Carlo steps for thermalization and averaging.
