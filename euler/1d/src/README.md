# Solves 1D Euler Equations for the sod shock tube

Author: Prashant Goswami
supervisor: Prof. Ashish Bhole

#1. Governing equations

The solver computes solutions to the one-dimensional compressible Euler
equations in conservative form:

```
dU/dt + dF(U)/dx = 0
```

where the conserved variable vector and flux function are

```
U = (rho, rho*u, rho*E)^T

F(U) = ( rho*u ,  rho*u^2 + p ,  u*(rho*E + p) )^T
```

with "rho" the density, "u" the velocity, "p" the pressure, and "rho*E" the
total energy density. The system is closed with the ideal gas equation of
state,

```
p = (gamma - 1) * ( rho*E - 0.5*rho*u^2 ),      gamma = 1.4
```

which also gives the local speed of sound, c = sqrt(gamma*p/rho).

Working in the conserved variables (rho, rho*u, rho*E) rather than the
primitive variables (rho, u, p) ensures that mass, momentum, and energy are
transported exactly according to their physical flux balance at each grid
cell - a property required for the numerical solution to correctly capture
shocks and contact discontinuities with the correct propagation speeds and
jump conditions.

#2. Spatial discretization

The domain x ∈ [0,1] is discretized on a uniform grid of Np points using
PETSc's DMDA (Distributed Array) with "dof = 3" (one block of 3 conserved
variables per grid point) and DM_BOUNDARY_GHOSTED, which allocates ghost
storage at both domain edges without imposing any built-in periodicity or
extrapolation.

The spatial derivative of the flux, dF/dx, is computed by a selectable
finite-difference operator, applied independently to each of the three flux
components:

The derivative operator itself is independent of the equation being solved —
it accepts any array and returns its spatial derivative, which allows the
same schemes, with identical coefficients, to serve both the Euler
solver and the underlying linear advection solver.

#3. Boundary conditions

The Sod shock tube is posed as an effectively unbounded domain: sufficiently
far from the initial discontinuity, the gas remains at its original,
undisturbed state for the duration of the simulation. This is enforced via a
"Dirichlet boundary condition": the ghost cells at the left and right
boundaries are fixed, at every timestep, to the corresponding conservative
form of the prescribed left and right states, "(rho_L, u_L, p_L)" and
"(rho_R, u_R, p_R)".

#4. Time integration and stability

Time integration uses PETSc's TS (Time Stepper) framework with an
explicit SSP Runge-Kutta scheme (TSSSP). The timestep is set from the CFL
condition using the maximum characteristic wave speed of the system,

```
dt = cfl * dx / max( |u| + c )
```

evaluated from the prescribed left/right states, since the Sod problem's
initial data consists of only two constant regions and this bound is
sufficient for the parameters used here.

#5. Shock capturing: artificial dissipation

Central and compact finite-difference schemes assume local smoothness of the
solution; at a genuine discontinuity (a shock or contact) this assumption is
violated, producing spurious oscillatory overshoot (the Gibbs phenomenon)
that is amplified, rather than damped, by these dissipation-free schemes
leading to unbounded growth and numerical blow-up if left uncorrected.

To stabilize the solution near discontinuities, a second-order artificial
dissipation term is added to the right-hand side of each conserved variable's
evolution equation:

```
res_k = -dF_k/dx + nu_art * d^2(U_k)/dx^2,     nu_art = nu_coeff * wave_speed_max * dx
```

where d^2/dx^2 is a standard 3-point discrete Laplacian and nu_coeff is a
runtime-configurable coefficient. Scaling the coefficient by the local grid
spacing keeps the magnitude of the added dissipation matched to grid
resolution under refinement. This is a standard artificial-viscosity
approach in the spirit of von Neumann & Richtmyer (1950); the solver
currently uses a fixed, uniform coefficient rather than a
discontinuity-adaptive sensor (e.g. Jameson, Schmidt & Turkel, 1981).

#6. Verification

The implementation is verified against the classical Sod shock tube
benchmark: domain [0,1], discontinuity at x = 0.5, left state
(rho,u,p) = (1, 0, 1), right state (rho,u,p) = (0.125, 0, 0.1),
gamma = 1.4. An exact solution (Toro's Riemann solver: Newton-Raphson
solution for the star-region pressure and velocity, followed by
similarity-variable sampling of the resulting rarefaction fan, contact
discontinuity, and shock) is used as the reference solution for direct
comparison, in both conservative and primitive variables.

The numerical solution correctly reproduces the three-wave structure of the
Sod problem — a left-propagating rarefaction fan, a contact discontinuity,
and a right-propagating shock — with density, velocity, and pressure
profiles matching the exact solution closely across all three regions.

#7. Known limitations and ongoing work

- The artificial dissipation coefficient is uniform across the domain and
  must be manually tuned per grid resolution and scheme order; a
  discontinuity-adaptive dissipation scheme is under consideration as a
  more robust alternative.
- The CFL timestep is currently computed only once from the initial condition and used throughout the computation,
  rather than adaptively from the evolving solution's local wave speeds.
