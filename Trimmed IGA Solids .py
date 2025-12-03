#!/usr/bin/env python
# coding: utf-8

# In[4]:


"""
IGA thin plate project: cantilever benchmark + trimmed hole + dynamics
======================================================================

Part A – Cantilever benchmark (no hole)
---------------------------------------
- NURBS-based isogeometric plate model (plane-stress).
- Thin rectangular plate used as a cantilever beam.
- Left edge clamped, right edge loaded in y.
- Compare IGA vertical deflection along the mid-height vs Euler–Bernoulli
  beam theory.

Part B – Trimmed plate with hole + dynamic loading
--------------------------------------------------
- Same IGA plate but with a circular "trimmed" hole (integration skipped inside).
- Static tension on right edge, left edge clamped.
- Explicit central-difference dynamic simulation with ramped load.
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import leggauss
import matplotlib.tri as mtri


# ----------------------------------------------------------------------
# B-spline / NURBS utilities
# ----------------------------------------------------------------------

def find_span(n, p, u, U):
    """Find knot span index for parameter u (B-spline)."""
    if u == U[n+1]:
        return n
    low = p
    high = n + 1
    mid = (low + high) // 2
    while u < U[mid] or u >= U[mid+1]:
        if u < U[mid]:
            high = mid
        else:
            low = mid
        mid = (low + high) // 2
    return mid


def basis_funs(i, u, p, U):
    """B-spline basis functions N_i,p(u) for span i."""
    N = np.zeros(p+1)
    left = np.zeros(p+1)
    right = np.zeros(p+1)
    N[0] = 1.0
    for j in range(1, p+1):
        left[j] = u - U[i+1-j]
        right[j] = U[i+j] - u
        saved = 0.0
        for r in range(j):
            temp = N[r] / (right[r+1] + left[j-r])
            N[r] = saved + right[r+1] * temp
            saved = left[j-r] * temp
        N[j] = saved
    return N


def ders_basis_funs(i, u, p, n, U):
    """
    Derivatives of B-spline basis up to order n.
    Returns ders[k,j] = d^k N_{i-p+j}(u) / du^k
    """
    ndu = np.zeros((p+1, p+1))
    ndu[0, 0] = 1.0
    left = np.zeros(p+1)
    right = np.zeros(p+1)

    for j in range(1, p+1):
        left[j] = u - U[i+1-j]
        right[j] = U[i+j] - u
        saved = 0.0
        for r in range(j):
            ndu[j, r] = right[r+1] + left[j-r]
            temp = ndu[r, j-1] / ndu[j, r]
            ndu[r, j] = saved + right[r+1] * temp
            saved = left[j-r] * temp
        ndu[j, j] = saved

    ders = np.zeros((n+1, p+1))
    for j in range(p+1):
        ders[0, j] = ndu[j, p]

    a = np.zeros((2, p+1))
    for r in range(p+1):
        s1 = 0
        s2 = 1
        a[0, 0] = 1.0
        for k in range(1, n+1):
            d = 0.0
            rk = r - k
            pk = p - k
            if r >= k:
                a[s2, 0] = a[s1, 0] / ndu[pk+1, rk]
                d = a[s2, 0] * ndu[rk, pk]
            j1 = 1 if rk >= -1 else -rk
            j2 = k-1 if (r-1) <= pk else p-r
            for j in range(j1, j2+1):
                a[s2, j] = (a[s1, j] - a[s1, j-1]) / ndu[pk+1, rk+j]
                d += a[s2, j] * ndu[rk+j, pk]
            if r <= pk:
                a[s2, k] = -a[s1, k-1] / ndu[pk+1, r]
                d += a[s2, k] * ndu[r, pk]
            ders[k, r] = d
            s1, s2 = s2, s1

    # factorial factors
    r = p
    for k in range(1, n+1):
        for j in range(p+1):
            ders[k, j] *= r
        r *= (p-k)
    return ders


class NURBSSurface:
    """
    Simple NURBS surface in 2D (x,y) with given degree p,q, knots U,V,
    control points (n_u x n_v x 2) and optional weights (n_u x n_v).
    """

    def __init__(self, p, q, U, V, control_points, weights=None):
        self.p = p
        self.q = q
        self.U = np.array(U, float)
        self.V = np.array(V, float)
        self.ctrl = np.array(control_points, float)  # (n_u, n_v, 2)
        self.n_u = self.ctrl.shape[0]
        self.n_v = self.ctrl.shape[1]
        if weights is None:
            self.w = np.ones((self.n_u, self.n_v))
        else:
            self.w = np.array(weights, float)

    def span_u(self, u):
        return find_span(self.n_u - 1, self.p, u, self.U)

    def span_v(self, v):
        return find_span(self.n_v - 1, self.q, v, self.V)

    def basis_and_derivs(self, u, v):
        """
        Evaluate NURBS basis and derivatives at (u,v):
            R, dR_dxi, dR_deta, J1, x_phys, spans
        """
        p, q = self.p, self.q
        U, V = self.U, self.V

        uspan = self.span_u(u)
        vspan = self.span_v(v)

        Nu = basis_funs(uspan, u, p, U)
        Nv = basis_funs(vspan, v, q, V)
        dNu = ders_basis_funs(uspan, u, p, 1, U)[1]
        dNv = ders_basis_funs(vspan, v, q, 1, V)[1]

        nen = (p+1)*(q+1)
        R = np.zeros(nen)
        dR_dxi = np.zeros(nen)
        dR_deta = np.zeros(nen)

        w_sum = 0.0
        dw_dxi = 0.0
        dw_deta = 0.0
        pts = []

        idx = 0
        for j in range(q+1):
            v_idx = vspan - q + j
            for i in range(p+1):
                u_idx = uspan - p + i
                w_ = self.w[u_idx, v_idx]
                Nuv = Nu[i] * Nv[j]
                dN_dxi = dNu[i] * Nv[j]
                dN_deta = Nu[i] * dNv[j]

                R[idx] = Nuv * w_
                dR_dxi[idx] = dN_dxi * w_
                dR_deta[idx] = dN_deta * w_
                w_sum += R[idx]
                dw_dxi += dR_dxi[idx]
                dw_deta += dR_deta[idx]
                pts.append(self.ctrl[u_idx, v_idx])
                idx += 1

        pts = np.array(pts)  # (nen,2)

        # Rationalize
        R /= w_sum
        dR_dxi = (dR_dxi * w_sum - R * dw_dxi) / (w_sum**2)
        dR_deta = (dR_deta * w_sum - R * dw_deta) / (w_sum**2)

        # Map to physical space
        x_phys = np.zeros(2)
        dx_dxi = np.zeros(2)
        dx_deta = np.zeros(2)
        for a in range(nen):
            x_phys += R[a] * pts[a]
            dx_dxi += dR_dxi[a] * pts[a]
            dx_deta += dR_deta[a] * pts[a]

        J1 = np.column_stack((dx_dxi, dx_deta))  # 2x2
        return R, dR_dxi, dR_deta, J1, x_phys, (uspan, vspan)


# ----------------------------------------------------------------------
# Mesh & assembly utilities
# ----------------------------------------------------------------------

def make_open_uniform_knots(n_ctrl, degree):
    """Open uniform knot vector for given number of control points and degree."""
    n_el = n_ctrl - degree
    kv = [0.0]*(degree+1)
    for i in range(1, n_el):
        kv.append(i/float(n_el))
    kv += [1.0]*(degree+1)
    return np.array(kv, float)


def build_knot_spans(U, p):
    """Return list of (span_index, u_start, u_end) for non-zero knot intervals."""
    spans = []
    m = len(U) - 1
    for i in range(p, m-p):
        if U[i] != U[i+1]:
            spans.append((i, U[i], U[i+1]))
    return spans


def assemble_iga_plane_stress(surface, thickness, E, nu, rho,
                              hole_center=None, hole_radius=0.0,
                              ngauss=3):
    """
    Assemble global stiffness K and consistent mass M for 2D plane-stress IGA
    on a NURBS surface, with optional circular 'hole' trimming.
    """
    p, q = surface.p, surface.q
    U, V = surface.U, surface.V
    spans_u = build_knot_spans(U, p)
    spans_v = build_knot_spans(V, q)

    ncp = surface.n_u * surface.n_v
    ndof = 2 * ncp
    K = np.zeros((ndof, ndof))
    M = np.zeros((ndof, ndof))

    # Plane-stress material matrix
    C = E / (1 - nu**2) * np.array([
        [1,   nu,          0],
        [nu,  1,           0],
        [0,   0,  (1 - nu)/2]
    ])

    # Gauss points in [-1,1]
    gp, gw = leggauss(ngauss)

    for iu_span, u1, u2 in spans_u:
        for iv_span, v1, v2 in spans_v:
            # Local connectivity (control point indices)
            loc_ctrl = []
            for j in range(q+1):
                v_idx = iv_span - q + j
                for i in range(p+1):
                    u_idx = iu_span - p + i
                    loc_ctrl.append((u_idx, v_idx))
            loc_ctrl = np.array(loc_ctrl, int)
            nen = loc_ctrl.shape[0]

            Ke = np.zeros((2*nen, 2*nen))
            Me = np.zeros((2*nen, 2*nen))

            for a in range(ngauss):
                xi = gp[a]; wi = gw[a]
                u = 0.5 * ((u2 - u1)*xi + (u2 + u1))
                du_dxi = 0.5 * (u2 - u1)
                for b in range(ngauss):
                    eta = gp[b]; wj = gw[b]
                    v = 0.5 * ((v2 - v1)*eta + (v2 + v1))
                    dv_deta = 0.5 * (v2 - v1)
                    J2 = du_dxi * dv_deta  # param jacobian

                    R, dR_dxi, dR_deta, J1, x_phys, _ = surface.basis_and_derivs(u, v)
                    detJ1 = np.linalg.det(J1)
                    if detJ1 <= 0:
                        raise ValueError("Negative or zero Jacobian in physical mapping")

                    # Trimming: skip points inside the circular hole
                    if hole_center is not None and hole_radius > 0.0:
                        dx = x_phys[0] - hole_center[0]
                        dy = x_phys[1] - hole_center[1]
                        if dx*dx + dy*dy < hole_radius**2:
                            continue

                    invJ1 = np.linalg.inv(J1)
                    # Map to physical derivatives dR/dx, dR/dy
                    dRdx = invJ1[0, 0]*dR_dxi + invJ1[0, 1]*dR_deta
                    dRdy = invJ1[1, 0]*dR_dxi + invJ1[1, 1]*dR_deta

                    B = np.zeros((3, 2*nen))
                    for a_loc in range(nen):
                        B[0, 2*a_loc]     = dRdx[a_loc]   # exx
                        B[1, 2*a_loc+1]   = dRdy[a_loc]   # eyy
                        B[2, 2*a_loc]     = dRdy[a_loc]   # exy
                        B[2, 2*a_loc+1]   = dRdx[a_loc]

                    weight = wi * wj * detJ1 * J2 * thickness

                    Ke += B.T @ C @ B * weight

                    # Consistent mass matrix
                    N = R
                    for i_loc in range(nen):
                        for j_loc in range(nen):
                            m_ij = rho * thickness * N[i_loc] * N[j_loc] * wi * wj * detJ1 * J2
                            Me[2*i_loc,     2*j_loc]     += m_ij
                            Me[2*i_loc + 1, 2*j_loc + 1] += m_ij

            # Assemble element into global K, M
            for a_loc, (iu, iv) in enumerate(loc_ctrl):
                A = 2 * (iv * surface.n_u + iu)
                for b_loc, (ju, jv) in enumerate(loc_ctrl):
                    Bgi = 2 * (jv * surface.n_u + ju)
                    K[A:A+2, Bgi:Bgi+2] += Ke[2*a_loc:2*a_loc+2,
                                              2*b_loc:2*b_loc+2]
                    M[A:A+2, Bgi:Bgi+2] += Me[2*a_loc:2*a_loc+2,
                                              2*b_loc:2*b_loc+2]

    return K, M


# ----------------------------------------------------------------------
# BCs, static & dynamic solvers
# ----------------------------------------------------------------------

def apply_bc_and_solve_static(K, F, fixed_dofs):
    """Solve K u = F with homogeneous Dirichlet BCs at fixed_dofs."""
    fixed_dofs = np.array(fixed_dofs, int)
    all_dofs = np.arange(K.shape[0])
    free = np.setdiff1d(all_dofs, fixed_dofs)

    Kff = K[np.ix_(free, free)]
    Ff = F[free]

    u = np.zeros(K.shape[0])
    u[free] = np.linalg.solve(Kff, Ff)
    return u


def lump_mass(M):
    """Row-sum lumping of consistent mass matrix."""
    m_lumped = np.sum(M, axis=1)
    return np.diag(m_lumped)


def central_difference(K, M_lumped, F_time, fixed_dofs, dt, n_steps):
    """
    Explicit central-difference time integration:
        M u_ddot + K u = F(t)
    with lumped mass (diagonal M_lumped).
    """
    ndof = K.shape[0]
    fixed_dofs = np.array(fixed_dofs, int)
    free = np.setdiff1d(np.arange(ndof), fixed_dofs)

    Kf = K[np.ix_(free, free)]
    Mf_diag = np.diag(M_lumped)[free]
    invM = 1.0 / Mf_diag

    u = np.zeros(ndof)
    u_old = np.zeros(ndof)
    U_hist = []

    for n in range(n_steps):
        t = n*dt
        F = F_time(t)
        Ff = F[free]
        Rf = Kf @ u[free]   # internal forces
        a_n = invM * (Ff - Rf)

        # central difference update
        u_new_free = 2*u[free] - u_old[free] + dt*dt*a_n
        u_old[free] = u[free]
        u[free] = u_new_free

        U_hist.append(u.copy())

    return np.array(U_hist)


# ----------------------------------------------------------------------
# PART A: Cantilever benchmark (IGA vs Euler–Bernoulli)
# ----------------------------------------------------------------------

def run_cantilever_benchmark():
    """
    IGA cantilever plate acting as a beam.
    Compares IGA mid-height vertical deflection vs Euler–Bernoulli theory.
    Geometry & material follow your QUAD8 cantilever setup.
    """
    print("\n=== PART A: IGA cantilever benchmark ===")

    groupnumber = 27
    L = 4.0 + 0.05*groupnumber   # Length (same as your MATLAB script)
    b = 0.5                      # Height
    t = 0.01                     # Thickness
    E = 200e9
    nu = 0.25
    rho = 7800.0
    F_total = 1.0                # Total vertical load at free end

    # Second moment of area (same convention as your MATLAB code)
    I = t * b**3 / 12.0

    # NURBS discretization
    p = q = 2
    n_u, n_v = 17, 5      # control points in each direction (fairly coarse)
    U = make_open_uniform_knots(n_u, p)
    V = make_open_uniform_knots(n_v, q)

    ctrl = np.zeros((n_u, n_v, 2))
    for i in range(n_u):
        x = L * i/(n_u-1)
        for j in range(n_v):
            y = b * j/(n_v-1)
            ctrl[i, j, 0] = x
            ctrl[i, j, 1] = y

    surface = NURBSSurface(p, q, U, V, ctrl)

    # Assemble K (no hole)
    K, M = assemble_iga_plane_stress(surface, t, E, nu, rho,
                                     hole_center=None, hole_radius=0.0,
                                     ngauss=3)
    ndof = K.shape[0]

    # Load: distribute F_total vertically on right edge
    F = np.zeros(ndof)
    right_nodes = []
    for j in range(n_v):
        iu = n_u - 1
        iv = j
        gnode = iv * n_u + iu
        right_nodes.append(gnode)
    Fy_per = -F_total / len(right_nodes)   # negative = downward
    for gnode in right_nodes:
        F[2*gnode + 1] += Fy_per  # y-direction DOF

    # BC: clamp left edge (u = v = 0)
    fixed_dofs = []
    for j in range(n_v):
        iu = 0
        iv = j
        gnode = iv * n_u + iu
        fixed_dofs += [2*gnode, 2*gnode+1]

    # Solve static problem
    u_static = apply_bc_and_solve_static(K, F, fixed_dofs)

    # Extract mid-height deflection along beam (y-direction)
    iv_mid = n_v // 2
    x_list = []
    v_iga = []
    for iu in range(n_u):
        gnode = iv_mid * n_u + iu
        x_list.append(ctrl[iu, iv_mid, 0])
        v_iga.append(u_static[2*gnode + 1])
    x_arr = np.array(x_list)
    v_iga = np.array(v_iga)

    # Euler–Bernoulli analytical deflection for cantilever with end load F
    # v(x) = -F x^2 (3L - x) / (6 E I)
    v_eb = -F_total * x_arr**2 * (3*L - x_arr) / (6.0 * E * I)

    # Print tip deflections
    print(f"Tip deflection (IGA)    : {v_iga[-1]:.6e} m")
    print(f"Tip deflection (EB beam): {v_eb[-1]:.6e} m")

    # Plot comparison
    plt.figure()
    plt.plot(x_arr, v_eb, '-o', label='Euler–Bernoulli')
    plt.plot(x_arr, v_iga, '-s', label='IGA (mid-height)')
    plt.xlabel("x [m]")
    plt.ylabel("Vertical deflection v [m]")
    plt.title("Cantilever beam: IGA vs Euler–Bernoulli")
    plt.grid(True)
    plt.legend()


# ----------------------------------------------------------------------
# PART B: Trimmed plate with hole + dynamic loading
# ----------------------------------------------------------------------

def run_trimmed_plate_and_dynamic():
    """
    IGA plate with circular trimmed hole:
      - static tension on right edge, left edge clamped,
      - explicit dynamic ramped load.
    """
    print("\n=== PART B: Trimmed plate + dynamic loading ===")

    # Geometry
    Lx, Ly = 1.0, 0.5
    thickness = 0.01

    # Material
    E = 210e9
    nu = 0.3
    rho = 7800.0

    # NURBS discretization
    p = q = 2
    n_u, n_v = 12, 8
    U = make_open_uniform_knots(n_u, p)
    V = make_open_uniform_knots(n_v, q)

    ctrl = np.zeros((n_u, n_v, 2))
    for i in range(n_u):
        x = Lx * i/(n_u-1)
        for j in range(n_v):
            y = Ly * j/(n_v-1)
            ctrl[i, j, 0] = x
            ctrl[i, j, 1] = y

    surface = NURBSSurface(p, q, U, V, ctrl)

    hole_center = (Lx/2.0, Ly/2.0)
    hole_radius = 0.1

    # Assemble K, M
    print("Assembling IGA stiffness and mass matrices...")
    K, M = assemble_iga_plane_stress(surface, thickness, E, nu, rho,
                                     hole_center=hole_center,
                                     hole_radius=hole_radius,
                                     ngauss=2)
    ndof = K.shape[0]

    # Static load case: pull on right edge in x-direction
    F = np.zeros(ndof)
    total_F = 1e4  # total force in x
    right_nodes = []
    for j in range(n_v):
        iu = n_u - 1
        iv = j
        gnode = iv * n_u + iu
        right_nodes.append(gnode)
    fx_per = total_F / len(right_nodes)
    for gnode in right_nodes:
        F[2*gnode] += fx_per  # x-direction

    # Clamp left edge
    fixed_dofs = []
    for j in range(n_v):
        iu = 0
        iv = j
        gnode = iv * n_u + iu
        fixed_dofs += [2*gnode, 2*gnode+1]

    # Static solve
    print("Solving static problem (plate with hole)...")
    u_static = apply_bc_and_solve_static(K, F, fixed_dofs)
    umag = np.sqrt(u_static[0::2]**2 + u_static[1::2]**2)
    print("Max static displacement magnitude:", umag.max())

    # Plot static displacement magnitude on control net
    fig, ax = plt.subplots()
    X = ctrl[:, :, 0]
    Y = ctrl[:, :, 1]
    Ux = u_static[0::2].reshape(n_v, n_u).T
    Uy = u_static[1::2].reshape(n_v, n_u).T
    Umag = np.sqrt(Ux**2 + Uy**2)

    # Mask points inside the hole so they don't show up in the plot
    R2 = (X - hole_center[0])**2 + (Y - hole_center[1])**2
    mask = R2 < hole_radius**2
    Umag_masked = np.ma.masked_where(mask, Umag)

    tri = mtri.Triangulation(X.flatten(), Y.flatten())
    h = ax.tripcolor(tri, Umag_masked.flatten(), shading='gouraud')
    
    plt.colorbar(h, ax=ax, label='|u| [m]')
    ax.set_title("Static displacement magnitude (IGA, trimmed hole)")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")

    # Dynamic analysis
    print("Setting up dynamic problem...")
    M_l = lump_mass(M)

    fixed = np.array(fixed_dofs, int)
    free = np.setdiff1d(np.arange(ndof), fixed)
    K_diag = np.diag(K)[free]
    M_diag = np.diag(M_l)[free]
    omega_sq = K_diag / M_diag
    omega_max = np.sqrt(omega_sq.max())
    dt_crit = 2.0 / omega_max
    dt = 0.2 * dt_crit
    print(f"Estimated critical dt ~ {dt_crit:.3e}, using dt = {dt:.3e}")

    # Time-dependent load: ramp to total_F over t_ramp, then hold
    t_ramp = 0.02
    def F_time(t):
        factor = min(t/t_ramp, 1.0)
        Fv = np.zeros(ndof)
        for gnode in right_nodes:
            Fv[2*gnode] += fx_per * factor
        return Fv

    n_steps = 500
    print("Running explicit dynamic simulation...")
    U_hist = central_difference(K, M_l, F_time, fixed_dofs, dt, n_steps)

    # Track x-displacement at mid-right node
    mid_j = n_v // 2
    gnode_mid_right = mid_j * n_u + (n_u-1)
    dof_track = 2*gnode_mid_right
    u_track = U_hist[:, dof_track]
    t_arr = np.arange(n_steps) * dt

    fig2, ax2 = plt.subplots()
    ax2.plot(t_arr, u_track)
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel("u_x at mid-right node [m]")
    ax2.set_title("Dynamic response (IGA plate with trimmed hole)")
    ax2.grid(True)


# ----------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------

if __name__ == "__main__":
    run_cantilever_benchmark()
    run_trimmed_plate_and_dynamic()
    plt.tight_layout()
    plt.show()


# In[ ]:




