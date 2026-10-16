# 2D Viscous-Inviscid Interaction (VII) Aerodynamic Solver

![MATLAB](https://img.shields.io/badge/MATLAB-R2023a-orange.svg)
![Aerodynamics](https://img.shields.io/badge/Domain-Aerodynamics-blue.svg)

## Overview
This repository contains the MATLAB source code for a 2D aerodynamic solver developed as part of a Bachelor's Thesis in Mechanical Engineering. 

The code evaluates the aerodynamic performance of 2D airfoils by relying on a **weakly coupled Viscous-Inviscid Interaction (VII)** approach. It circumvents the computational cost of full Navier-Stokes (CFD) simulations by coupling a potential flow solver with an integral boundary layer solver, providing fast and reliable estimations of pressure distributions, lift, drag, and flow separation for pre-stall conditions.

> **Full Thesis:** For a deep dive into the mathematical formulation and theoretical background, you can [read the full Thesis PDF here (Italian)](docs/Tesi_VII_Solver_Bunino.pdf).
> *Note: The source code in this repository has been translated to English and slightly refactored for an international audience. As a result, specific variable names and script structures may differ slightly from the original Italian code snippets mentioned in the thesis.* 

## Scope & Assumptions
Before using the solver, please note the following physical and mathematical assumptions underlying the code:
- **Incompressible Flow:** The fluid is assumed to have constant density (Mach number < 0.3).
- **Steady-State Flow:** The equations (e.g., von Kármán integral equation) are formulated for stationary conditions ($\partial/\partial t = 0$).
- **Reynolds Number Regime:** The integral boundary layer methods are suitable for moderate to high Reynolds numbers (e.g., $Re \sim 10^5 - 10^7$) typical of standard airfoil operations, where the thin boundary layer approximation holds.
- **Attached Flow / Pre-Stall:** The solver reliably predicts performance in the linear lift region. It is not designed to model massive flow separation or deep post-stall behaviors, as it lacks a dead-air wake model.

## Key Features
- **Inviscid Potential Flow Solver:** A low-order surface-panel method utilizing **constant-strength distributed sources and doublets** (unlike higher-order linear/quadratic methods). It is closed by the Kutta condition to ensure physically accurate circulation and lift generation.
- **Laminar Boundary Layer:** Solved using **Thwaites' Integral Method**, allowing for the calculation of the momentum thickness ($\theta$), displacement thickness ($\delta^*$), and wall shear stress ($\tau_w$).
- **Turbulent Boundary Layer & Transition:** Extended via the **Buri parameter** to handle turbulent regimes, with flow separation predicted using **Alber's criterion**.
- **VII Coupling Algorithm:** A weak coupling strategy where the inviscid flow provides the external velocity distribution, and the boundary layer solver feeds back the displacement thickness ($\delta^\*$) to virtually "thicken" the airfoil geometry. **When flow separation is predicted, $\delta^*$ is kept constant** (meaning the wake is not explicitly modeled). Subsequently, a **trailing edge "pinch"** (artificial closure of the displacement thickness) is enforced to prevent numerical divergence and ensure solver convergence. Due to this simplified post-separation treatment, the code is **not suitable for bluff bodies**, as demonstrated by the circular cylinder test case (you can try it yourself ;) ).
- **Advanced Visualization:** High-quality plotting capabilities, including 3D boundary layer extrusion and validation against the Blasius exact analytical solution. The boundary layer displacement thickness is artificially scaled up for better visual appreciation.

## Configuration & Input Parameters
The solver can be easily configured by modifying the initial parameters directly within the main script. Users can adjust the following key variables to set up their custom simulation:
- **Flow Conditions:** Reynolds number, fluid density, dynamic/kinematic viscosity, and angle of attack ($\alpha$).
- **Geometry & Discretization:** Chord length, number of panels for the inviscid solver, and the specific airfoil geometry file.
- **Study Case:** Switch between standard airfoil analysis (`'AIRFOIL'`) or flat plate validation (`'BLASIUS'`).
- **Visualization Scaling:** Because the physical boundary layer is extremely thin, its thickness ($\delta^*$) is scaled up in the plots for better visual clarity. You can customize these multipliers by changing the `thickness_plot_scale` variable for 3D plots and the `scale_2d` variable for 2D plots.

### Available Airfoil Coordinates
To run the solver on a specific geometry, you must manually update the coordinate file variable in the main script with the exact name of the `.dat` file. The repository comes pre-loaded with a variety of standard and high-performance airfoils:

| Filename | Geometry / Airfoil Type |
| :--- | :--- |
| `ag16_coords.dat` | AG16 |
| `ah83159_coords.dat` | AH 83-159 |
| `b737d_coords.dat` | Boeing 737 (b737d) |
| `circle_coords.dat` | Circular cylinder |
| `e421_coords.dat` | Eppler 421 |
| `n0010_coords.dat` | NACA 0010 (Symmetric) |
| `n2412_coords.dat` | NACA 2412 |
| `n23012_coords.dat` | NACA 23012 |
| `n4418_coords.dat` | NACA 4418 |
| `n64a410_coords.dat` | NACA 64A410 |
| `raf89_coords.dat` | RAF 89 |


## Validation & Results
The model has been extensively validated against canonical theoretical cases (e.g., Blasius flat plate) and `foil.tools` data for standard airfoils (e.g., NACA 0010, NACA 23012).

> **Figure 1:** *Validation of the Thwaites integral method against the Blasius exact analytical solution for a flat plate.*
> <img width="800" alt="BL" src="https://github.com/user-attachments/assets/34aad194-b57f-46d9-8e36-f0eacce93099" />

> **Figure 2:** *Aerodynamic polars (Lift and Drag coefficients) for the symmetric NACA 0010 airfoil, compared with foil.tools data.*
> <img width="800" alt="N0010" src="https://github.com/user-attachments/assets/8a039882-a9f6-40bb-a35e-bfe8520e7e37" />

> **Figure 3:** *Aerodynamic polars for the cambered NACA 23012 airfoil, showing good agreement in the linear pre-stall region.*
> <img width="800" alt="N23012" src="https://github.com/user-attachments/assets/a2d42baf-b016-4d2d-a493-640497ef9c5c" />

### Solver Limitations (Post-Stall & Massive Separation)
To transparently highlight the physical boundaries of the weakly coupled VII approach, the solver was also tested on high-camber geometries (e.g., the **Eppler 421** airfoil) and at high angles of attack. Under these conditions, massive flow separation occurs. Because the current model does not include a dead-air wake formulation, the weak coupling breaks down, leading to an overprediction of lift and an inability to capture deep post-stall aerodynamics.

> **Figure 4:** *Example of solver limitations on the Eppler 421 airfoil. The divergence from reference data highlights the breakdown of the weak viscous-inviscid coupling during massive separation.*
> <img width="2574" height="1664" alt="E421" src="https://github.com/user-attachments/assets/54272c57-aa6d-4ada-ac77-7ced6a25d7b2" />

> **Figure 5:** *Comparison of the Pressure Coefficient ($C_p$) distribution.*
> <img width="2528" height="1682" alt="Cp" src="https://github.com/user-attachments/assets/c984a5ac-3661-4a65-a8b0-1690eaa202fc" />

> **Figure 6a:** *3D visualization of the boundary layer displacement thickness for a b737d (Boeing) airfoil. Note that the boundary layer thickness has been artificially scaled up for better visual appreciation.*
> <img width="3008" height="1408" alt="3D" src="https://github.com/user-attachments/assets/f6742d55-5ca4-4093-8802-aa1c37a736a5" />

> **Figure 6b:** *2D visualization of the same scaled boundary layer displacement thickness over the airfoil.*
> <img width="3136" height="1344" alt="2D" src="https://github.com/user-attachments/assets/7979104c-d319-4968-b160-55c9f5aae35e" />


## Getting Started

### Prerequisites
- MATLAB (The code has been written using standard MATLAB libraries. No specific toolboxes are strictly required, but verify compatibility with your version).

