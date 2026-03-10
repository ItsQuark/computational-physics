# Restricted Three-Body Problem: Earth-Moon-Spacecraft Trajectory

This project implements a numerical solution for the restricted three-body problem, specifically simulating the trajectory of a spacecraft from Earth to the Moon. The model assumes a simplified system where the spacecraft's mass is negligible compared to the Earth and the Moon, thus not perturbing their motion.

## Physical Model

The system is defined by three bodies: the Earth (MT), the Moon (ML), and the spacecraft (m). The Earth is considered fixed at the origin, while the Moon follows a circular orbit with constant angular velocity $\omega$. 

The dynamics are derived using the Hamiltonian formalism in polar coordinates $(r, \phi)$, where $r$ is the distance from the Earth. The distance from the spacecraft to the Moon ($r_L$) is given by:

$$r_L(r, \phi, t) = \sqrt{r^2 + d_{TL}^2 - 2rd_{TL}\cos(\phi - \omega t)}$$

The Hamiltonian of the system is:

$$H = \frac{p_r^2}{2m} + \frac{p_\phi^2}{2mr^2} - G\frac{m M_T}{r} - G\frac{m M_L}{r_L(r, \phi, t)}$$

The resulting equations of motion form a set of four coupled non-linear ordinary differential equations (ODEs) for $\dot{r}$, $\dot{\phi}$, $\dot{p_r}$, and $\dot{p_\phi}$.

## Numerical Implementation

### Runge-Kutta 4th Order (RK4)
The system of ODEs is solved using a classic 4th-order Runge-Kutta algorithm. This method provides a balance between computational efficiency and numerical precision, with a local truncation error of $O(h^5)$.

### Adaptive Step Size
To handle the disparate time scales of the mission (short time scales near planetary bodies and long scales during translunar injection), an adaptive step size control is implemented. The algorithm monitors the relative error $\epsilon$ by comparing a full step $h$ with two half-steps $h/2$:

$$\epsilon \approx \frac{16}{15} |y(t+h; h/2) - y(t+h; h)|$$

The step size is dynamically adjusted to keep the error within a user-defined tolerance $\epsilon_{max}$.

## Conservation Laws and Validation
The Jacobi integral, $H_0 = H - \omega p_\phi$, is a constant of motion in the rotating frame. This conservation law is used as the primary metric to validate the accuracy of the numerical integration and to monitor the accumulated error of the adaptive RK4 algorithm.

## Technical Requirements
- **Language:** Fortran 90/95
- **Visualization:** Gnuplot
- **Physical Constants:** Extracted from Planetary Data System (NASA)
