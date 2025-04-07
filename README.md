# Numerical Methods Solver Suite

A curated MATLAB collection of numerical analysis solvers implementing root-finding, integration, systems of nonlinear equations, ODE solvers, and fluid dynamics simulations. Developed as part of advanced coursework in numerical methods.

##  Features

### Root-Finding Techniques
- **Bisection**, **Secant**, and **Newton-Raphson** methods for solving non-linear equations (e.g., spherical tank depth, projectile motion).

### Numerical Integration
- **Trapezoidal**, **Simpson’s**, and **Gauss-Legendre Quadrature** to estimate physical quantities like magnetic field from Biot-Savart law.

### Nonlinear Systems
- Newton-Raphson method to solve a system of biochemical equations involving substrate and biomass concentration.

### ODE Solvers
- Simulation of spring-mass-damper system using **Euler** and **Runge-Kutta (RK4)** methods.
- Analysis of projectile trajectory under gravity and air drag.

### Finite Difference Method (FDM)
- Couette flow profile analysis using second-order FDM and validation with exact solution.

---

##  File Overview

| File             | Description                                                         |
|------------------|---------------------------------------------------------------------|
| `Exercise1.m`    | Root finding for a spherical tank volume                           |
| `Exercise2.m`    | Biot-Savart magnetic field calculation using 3 integration methods  |
| `Exercise3.m`    | Solves nonlinear biochemical system with Newton-Raphson            |
| `Exercise4.m`    | Solves spring-mass-damper with Euler and RK4 methods               |
| `Exercise5.m`    | Couette flow velocity profile using Finite Difference Method        |
| `Exercise6.m`    | Peak altitude computation for projectile under drag using Newton-Raphson |
| `test.m`         | Duplicate of `Exercise4` for quick testing                          |
| `FinalExamReport.pdf` | Includes explanations, methodology, and results for all exercises |

---

##  How to Run

1. Open MATLAB.
2. Place all files in the same folder.
3. Run any script by typing:

```matlab
Exercise1
```

Replace `1` with `2` to `6` for other exercises.

---

##  Author

**Yasmine Elsisi**  

