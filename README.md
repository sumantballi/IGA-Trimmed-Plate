# IGA-Trimmed-Plate
NURBS-based IGA plate with trimmed hole and dynamic loading

# NURBS-based IGA plate with trimmed hole and dynamic loading

Minimal isogeometric analysis (IGA) workflow for a thin plate with a trimmed region,
implemented in Python. The project focuses on thin-walled behaviour, validation against
classical beam theory, and a first step towards trimmed IGA solids under dynamic loading.

## Introduction

This project implements a small isogeometric analysis (IGA) workflow for a thin, rectangular plate. The plate mid-surface is represented by a NURBS-based IGA model, which is used to:

validate the formulation on a cantilever benchmark against Euler–Bernoulli beam theory, and

study the static and dynamic response of a plate with a central trimmed (cut-out) region under in-plane loading.

The focus is on thin-walled behaviour, trimmed geometry handling, and explicit dynamic response, using only open-source tools (Python, NumPy, Matplotlib)

The workflow includes:

- A rectangular thin plate discretised with a NURBS surface (isogeometric plane-stress).
- A cantilever benchmark where IGA mid-height deflections are compared with
  Euler–Bernoulli beam theory.
- A plate with a central trimmed region (circular cut-out approximated by excluded
  integration) under in-plane tension.
- An explicit dynamic response of the trimmed plate under a ramped edge load.


## What to look for

- **Cantilever benchmark**:
  - IGA mid-height deflection closely matches Euler–Bernoulli beam theory.
  - Shows that the NURBS-based assembly and BCs are correct.

- **Trimmed plate (static)**:
  - Displacement field for a plate with a circular soft/trimmed region under tension.
  - Larger displacements and “necking” near the cut-out.

- **Trimmed plate (dynamic)**:
  - Time history of edge displacement under a ramped load.
  - Undamped / damped behaviour depending on load ramp and Rayleigh damping settings.


## Scope

- Isogeometric **plane-stress** plate with NURBS shape functions.
- Approximate **trimmed region** via skipped integration points inside a circular hole.
- Linear elastic material, small-strain kinematics.
- Explicit dynamic response with optional Rayleigh damping.

Not included:

- Full 3D solids or shell-specific kinematics.
- Nonlinear material laws (plasticity, damage).
- Contact or multi-part assemblies.
- Industrial crash models – this is an IGA “pilot” component only.

## Results (quick look)

- **Cantilever**: IGA mid-height deflection vs Euler–Bernoulli:
  - Tip deflections differ only slightly (mesh/discretisation effects).
- **Static trimmed plate**: displacement magnitude field with a clear soft region
  around the hole.
- **Dynamic trimmed plate**: displacement–time curve at the mid-right node, showing
  either ringing (fast ramp, no damping) or a smoother approach (slower ramp / damping).

See `figs/` for the plots.

## Key methods & equations

- **NURBS basis**: standard B-spline basis functions with weights, assembled into a
  2D NURBS surface for the plate mid-surface.

- **Plane-stress stiffness**:
  - The element stiffness is assembled from
    \[
      \mathbf{K}_e = \int_{\Omega_e} \mathbf{B}^T \mathbf{C} \mathbf{B}\, t\,\mathrm{d}\Omega
    \]
    with:
    - `B` – strain–displacement matrix built from NURBS derivatives (mapped to x,y),
    - `C` – plane-stress elasticity matrix,
    - `t` – plate thickness.

- **Consistent mass + lumped mass**:
  - A consistent mass matrix is first assembled from NURBS shape functions,
    then row-sum lumping is applied for explicit dynamics.

- **Trimmed region (approximate)**:
  - Gauss points whose physical coordinates satisfy `r < r_hole` are skipped in the
    stiffness and mass integration, effectively softening/removing the material
    inside the circular cut-out.

- **Explicit dynamics**:
  - Uses a central-difference-type explicit update for
    \[
      \mathbf{M} \ddot{\mathbf{u}} + \mathbf{K} \mathbf{u} = \mathbf{f}(t),
    \]
    with optional Rayleigh damping
    \(\mathbf{C} = \alpha \mathbf{M} + \beta \mathbf{K}\).

## Future work

The current implementation is intentionally minimal and leaves several natural extensions:

1. **True trimmed IGA integration**
   - Replace the “skip Gauss points” approximation by exact integration over the trimmed domain (e.g. subdivision, boundary-fitted quadrature, or element-splitting techniques).
   - Compare the approximate and exact trimming strategies in terms of accuracy and robustness.

2. **Plate and shell kinematics**
   - Extend the formulation from plane-stress to plate/shell kinematics (Kirchhoff–Love or Mindlin–Reissner) to capture bending more realistically in thin-walled components.
   - Validate against analytical plate solutions and classical shell elements.

3. **Nonlinear material and failure**
   - Introduce elastoplastic or damage models to study localized failure around the cut-out.
   - Compare stress/strain and failure patterns against standard FEM shell models for the same geometry.

4. **More realistic dynamic loading**
   - Use more representative load histories (impact pulses, crash-like decelerations) and study mode participation and energy distribution.
   - Systematically investigate the influence of damping, ramp time, and mesh refinement on the dynamic response.

5. **Component-level comparison with industrial FEM**
   - Build an equivalent shell or solid model (e.g. in LS-DYNA/Abaqus) and compare stiffness, natural frequencies, and local stress concentrations with the IGA model.
   - Explore simple coupling strategies between IGA patches and standard FE meshes.

6. **Automation and workflow tooling**
   - Wrap the current scripts into a small Python package or command-line tool that reads case definitions from JSON/YAML.
   - Add automated generation of plots and summary tables to make the workflow reproducible and easy to reuse.
## References

1. Du, X., Zhao, G., Zhang, R., Wang, W., & Yang, J. (2022).  
   **“Numerical implementation for isogeometric analysis of thin-walled structures based on a Bézier extraction framework: nligaStruct.”**  
   *Thin-Walled Structures, 180, 109844.*  
  

2. Wang, Y., Jin, L., Yang, H., Hao, P., Ji, Y., & Wang, B. (2023).  
   **“Isogeometric-based mapping modeling and buckling analysis of stiffened panels with trimmed surfaces.”**  
   *Thin-Walled Structures, 186, 110676.*  
    

3. Chasapi, M., Antolín, P., & Buffa, A. (2024).  
   **“Fast parametric analysis of trimmed multi-patch isogeometric Kirchhoff–Love shells using a local reduced basis method.”**  
   *Engineering with Computers, 40, 3623–3650.*  
  

4. Yokoyama, Y., Nagasaka, K., Masuda, I., et al. (2024).  
   **“Isogeometric Analysis in Structural Deformation and Automobile Crash.”** *Preprints.org.*  
   

