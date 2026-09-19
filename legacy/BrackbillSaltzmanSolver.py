import numpy as np

class BrackbillSaltzmanSolver:
    def __init__(self, m=100, n=100):
        self.m = m
        self.n = n

    def relax_grid(self, R_init, Z_init, W_monitor, inner_iters=300, omega_sor=0.2, adapt_weight=1.0):
        """
        Performs Gauss-Seidel point relaxation on the coordinate fields
        to minimize smoothness, orthogonality, and volume variation.
        """
        # Create working copies of the physical coordinate grids
        R = np.copy(R_init)
        Z = np.copy(Z_init)

        # Grid spacings in computational domain (assumed unit square delta = 1)
        dxi = 1.0
        deta = 1.0

        # Safety epsilon to prevent divide-by-zero in calm fluid regions
        epsilon = 1e-8

        # ========================================================
        # NEW: Smooth the monitor function to remove numerical "cliffs"
        # ========================================================
        W_smooth = np.copy(W_monitor)
        for _ in range(4): # 4 passes of gentle mathematical blurring
            W_smooth[1:-1, 1:-1] = 0.5 * W_smooth[1:-1, 1:-1] + 0.125 * (
                W_smooth[2:, 1:-1] + W_smooth[:-2, 1:-1] +
                W_smooth[1:-1, 2:] + W_smooth[1:-1, :-2]
            )

        # Update boundaries of the smoothed array
        W_smooth[0, :] = W_smooth[1, :]
        W_smooth[-1, :] = W_smooth[-2, :]
        W_smooth[:, 0] = W_smooth[:, 1]
        W_smooth[:, -1] = W_smooth[:, -2]

        # 1. Pre-calculate the monitor function gradients using W_SMOOTH
        P_forcing = np.zeros_like(W_smooth)
        Q_forcing = np.zeros_like(W_smooth)

        # Central differences for internal nodes of the weight map
        for i in range(1, self.m - 1):
            for j in range(1, self.n - 1):
                dW_dxi = (W_smooth[i+1, j] - W_smooth[i-1, j]) / (2.0 * dxi)
                dW_deta = (W_smooth[i, j+1] - W_smooth[i, j-1]) / (2.0 * deta)

                # B&S localized pull parameters with safety epsilon
                P_forcing[i, j] = dW_dxi / (W_smooth[i, j] + epsilon)
                Q_forcing[i, j] = dW_deta / (W_smooth[i, j] + epsilon)

        print(f"Max P Forcing: {np.max(np.abs(P_forcing)):.4f}")
        print(f"Max Q Forcing: {np.max(np.abs(Q_forcing)):.4f}")

        # 2. Main Gauss-Seidel Relaxation Loop
        for sweep in range(inner_iters):
            # Loop over all interior nodes (leaving boundaries anchored)
            for i in range(1, self.m - 1):
                for j in range(1, self.n - 1):

                    # Compute localized 1st-order transformation metrics
                    r_xi = (R[i+1, j] - R[i-1, j]) / 2.0
                    r_eta = (R[i, j+1] - R[i, j-1]) / 2.0
                    z_xi = (Z[i+1, j] - Z[i-1, j]) / 2.0
                    z_eta = (Z[i, j+1] - Z[i, j-1]) / 2.0

                    # Compute the local Grid Jacobian (Cell Volume tracking)
                    J = r_xi * z_eta - r_eta * z_xi
                    if abs(J) < 1e-9:
                        J = 1.0 if J >= 0 else -1.0 # Protect against zero-volume cell collapse

                    # Compute B&S elliptic metric coefficients
                    alpha = r_eta**2 + z_eta**2
                    beta  = r_xi * r_eta + z_xi * z_eta
                    gamma = r_xi**2 + z_xi**2

                    # Compute mixed cross-derivatives (for the non-orthogonality penalty)
                    r_xi_eta = (R[i+1, j+1] - R[i+1, j-1] - R[i-1, j+1] + R[i-1, j-1]) / 4.0
                    z_xi_eta = (Z[i+1, j+1] - Z[i+1, j-1] - Z[i-1, j+1] + Z[i-1, j-1]) / 4.0

                    # 3. Formulate the explicit updates for R and Z nodes
                    # Added 'adapt_weight' to scale down the extreme right-hand source terms
                    # 1. Restore J**2 to keep the math dimensionally stable
                    r_rhs = 2.0 * beta * r_xi_eta - adapt_weight * abs(J) * (P_forcing[i, j] * r_xi + Q_forcing[i, j] * r_eta)
                    z_rhs = 2.0 * beta * z_xi_eta - adapt_weight * abs(J) * (P_forcing[i, j] * z_xi + Q_forcing[i, j] * z_eta)

                    # 2. Calculate the raw target coordinates
                    r_target_raw = (alpha * (R[i+1, j] + R[i-1, j]) + gamma * (R[i, j+1] + R[i, j-1]) - r_rhs) / (2.0 * (alpha + gamma))
                    z_target_raw = (alpha * (Z[i+1, j] + Z[i-1, j]) + gamma * (Z[i, j+1] + Z[i, j-1]) - z_rhs) / (2.0 * (alpha + gamma))

                    # 3. Calculate the raw requested physical shift
                    r_shift = r_target_raw - R[i, j]
                    z_shift = z_target_raw - Z[i, j]

                    # 4. THE PROPORTIONAL LIMITER (Strict 10%)
                    max_r_step = 0.10 * abs(R[i+1, j] - R[i-1, j])
                    max_z_step = 0.10 * abs(Z[i, j+1] - Z[i, j-1])

                    r_shift_clipped = np.clip(r_shift, -max_r_step, max_r_step)
                    z_shift_clipped = np.clip(z_shift, -max_z_step, max_z_step)

                    # 5. Apply the safe shift
                    R[i, j] = R[i, j] + omega_sor * r_shift_clipped
                    Z[i, j] = Z[i, j] + omega_sor * z_shift_clipped

            # 4. Enforce Boundary Sliding Rules after each sweep
            R, Z = self._apply_boundary_conditions(R, Z)

        return R, Z

    def _apply_boundary_conditions(self, R, Z):
        """
        Enforces physical anchors while letting lines slide along boundaries
        orthogonally, preventing boundary node overlap.
        """
        # Left Axis (Symmetry): R locked at 0, Z perfectly mirrors interior (orthogonal)
        R[0, :] = 0.0
        Z[0, :] = Z[1, :]

        # Right Wall (Cold Mold): R locked at 1, Z perfectly mirrors interior
        R[-1, :] = 1.0
        Z[-1, :] = Z[-2, :]

        # Bottom Wall (Baseplate): Z locked at 0, R perfectly mirrors interior
        Z[:, 0] = 0.0
        R[:, 0] = R[:, 1]

        # Top Boundary (Melt Pool Surface): Z locked flat, R perfectly mirrors interior
        Z[:, -1] = 1.5
        R[:, -1] = R[:, -2]

        # ==========================================
        # NEW: ABSOLUTE CORNER LOCKS (The Safety Pins)
        # ==========================================
        # This guarantees the physical domain never warps or shrinks
        R[0, 0], Z[0, 0]   = 0.0, 0.0  # Bottom-Left
        R[-1, 0], Z[-1, 0] = 1.0, 0.0  # Bottom-Right
        R[0, -1], Z[0, -1] = 0.0, 1.5  # Top-Left
        R[-1, -1], Z[-1, -1] = 1.0, 1.5  # Top-Right

        return R, Z