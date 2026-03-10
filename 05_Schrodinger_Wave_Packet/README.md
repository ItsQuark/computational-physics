# Time-Dependent Schrödinger Equation: Cayley's Form and Crank-Nicolson Method

This project focuses on the numerical resolution of the one-dimensional Time-Dependent Schrödinger Equation (TDSE) for a quantum particle. The implementation ensures the conservation of the wave function's norm by using a unitary evolution operator based on Cayley's approximation.

## Theoretical Background

The state of a particle is described by the complex wave function $\psi(x, t)$. In dimensionless units ($\hbar = 1$, $m = 1/2$), the TDSE is expressed as:

$$i \frac{\partial \psi(x, t)}{\partial t} = \left[ -\frac{\partial^2}{\partial x^2} + V(x) \right] \psi(x, t) = H\psi(x, t)$$

The formal solution involves the time-evolution operator:
$$\psi(x, t) = e^{-i(t-t_0)H}\psi(x, t_0)$$

A critical requirement for any numerical scheme in quantum mechanics is to maintain the normalization $\int |\psi|^2 dx = 1$, which corresponds to the conservation of total probability.

## Numerical Method

### Discretization and Unitarity
The space is discretized into a lattice $x_j = j\Delta x$ and time into steps $t_n = n\Delta t$. To avoid the non-unitary nature of standard Taylor expansions, we employ **Cayley's approximation** (equivalent to the Crank-Nicolson scheme):

$$e^{-i\Delta t H} \approx \frac{1 - i\Delta t H / 2}{1 + i\Delta t H / 2}$$

This operator is exactly unitary and accurate up to $O(\Delta t^2)$.

### The Tridiagonal Algorithm
The evolution step $\psi^{n+1} = \frac{1 - i\Delta t H / 2}{1 + i\Delta t H / 2} \psi^n$ is reformulated as a linear system:
$$\left[ 1 + \frac{i\Delta t H}{2} \right] q^n = 2\psi^n$$
where $\psi^{n+1} = q^n - \psi^n$.

Since the discretized Hamiltonian $H_D$ only involves nearest neighbors, the resulting system is **tridiagonal**. It is solved efficiently using a specialized recursion (Thomas Algorithm variant) to invert the matrix at each time step.

## Implementation Details

### Initial Conditions
Two types of initial states are simulated:
1. **Hamiltonian Eigenstates:** Stationary states of a particle in a 1D box.
2. **Gaussian Wave Packets:** A plane wave modulated by a Gaussian envelope, allowing for the study of dispersion and packet propagation.

### Observables and Validation
The simulation tracks and validates several physical quantities:
- **Norm Conservation:** Verification that $\int P(x, t) dx = 1$ remains constant over time.
- **Expectation Values:** Calculation of $\langle x(t) \rangle$ and $\langle p(t) \rangle$.
- **Heisenberg Uncertainty Principle:** Monitoring the product $\Delta x(t) \Delta p(t)$.

## Technical Specifications
- **Language:** Fortran 90/95 (utilizing `double complex` precision).
- **Spatial Resolution:** $S = 1000$ points.
- **Visualization:** Gnuplot animations of $\text{Re}[\psi]$, $\text{Im}[\psi]$, and $|\psi|^2$.
