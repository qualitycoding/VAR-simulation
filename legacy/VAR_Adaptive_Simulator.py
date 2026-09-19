import numpy as np
from scipy.interpolate import griddata
import os
import matplotlib.pyplot as plt

class VAR_Adaptive_Simulator:
    def __init__(self, m=100, n=100):
        self.m, self.n = m, n
        self.grid_solver = BrackbillSaltzmanSolver(m, n)
        # Initialize fields
        self.R = np.zeros((m, n)); self.Z = np.zeros((m, n))
        self.psi = np.zeros((m, n)); self.omega = np.zeros((m, n))
        self.theta = np.zeros((m, n)); self.F = np.zeros((m, n))
        self.U_comp = np.zeros((m, n)); self.W_comp = np.zeros((m, n))

    def calculate_grid_metrics(self):
        """
        Calculates the covariant and contravariant metric tensors required 
        to map the physical domain (R, Z) to the computational domain (xi, eta).
        Updates the grid Jacobian (J) and smoothing coefficients (alpha, beta, gamma).
        """
        R_xi = (self.R[2:, 1:-1] - self.R[:-2, 1:-1]) / 2.0
        R_eta = (self.R[1:-1, 2:] - self.R[1:-1, :-2]) / 2.0
        Z_xi = (self.Z[2:, 1:-1] - self.Z[:-2, 1:-1]) / 2.0
        Z_eta = (self.Z[1:-1, 2:] - self.Z[1:-1, :-2]) / 2.0

        self.alpha = R_eta**2 + Z_eta**2
        self.beta = R_xi * R_eta + Z_xi * Z_eta
        self.gamma = R_xi**2 + Z_xi**2
        self.J = R_xi * Z_eta - R_eta * Z_xi


    def solve_temperature(self, r_0=0.5, Pr=0.14, Nu_base=40.0, N_rad=0, epsilon=0.01):
        """
        Solves the energy equation using upwind differencing for convection 
        and central differencing for diffusion. Applies a Gaussian heat flux 
        at the top boundary to simulate the electrical arc, and convective 
        cooling boundary conditions at the crucible walls.
        """
        theta_old = self.theta.copy()
        theta_temp = self.theta.copy()

        a = self.alpha
        b = self.beta
        g = self.gamma
        J_int = self.J

        psi_xi  = (self.psi[2:, 1:-1] - self.psi[:-2, 1:-1]) / 2.0
        psi_eta = (self.psi[1:-1, 2:] - self.psi[1:-1, :-2]) / 2.0
        R_xi  = (self.R[2:, 1:-1] - self.R[:-2, 1:-1]) / 2.0
        R_eta = (self.R[1:-1, 2:] - self.R[1:-1, :-2]) / 2.0
        Z_xi  = (self.Z[2:, 1:-1] - self.Z[:-2, 1:-1]) / 2.0
        Z_eta = (self.Z[1:-1, 2:] - self.Z[1:-1, :-2]) / 2.0

        psi_Z = (psi_eta * R_xi - psi_xi * R_eta) / J_int
        psi_R = (psi_xi * Z_eta - psi_eta * Z_xi) / J_int

        R_interior = self.R[1:-1, 1:-1] + 1e-10
        U = psi_Z / R_interior
        W = -psi_R / R_interior

        theta_R_upwind = np.where(U > 0,
                                  self.theta[1:-1, 1:-1] - self.theta[:-2, 1:-1],
                                  self.theta[2:, 1:-1] - self.theta[1:-1, 1:-1])

        theta_Z_upwind = np.where(W > 0,
                                  self.theta[1:-1, 1:-1] - self.theta[1:-1, :-2],
                                  self.theta[1:-1, 2:] - self.theta[1:-1, 1:-1])

        convection = Pr * (U * theta_R_upwind + W * theta_Z_upwind)

        theta_cross = (self.theta[:-2, :-2] - self.theta[2:, :-2] -
                       self.theta[:-2, 2:] + self.theta[2:, 2:]) / 4.0

        theta_new_interior = (a * (self.theta[:-2, 1:-1] + self.theta[2:, 1:-1]) -
                              2.0 * b * theta_cross +
                              g * (self.theta[1:-1, :-2] + self.theta[1:-1, 2:]) -
                              convection * J_int**2) / (2.0 * (a + g))

        # 1. Stage the interior
        theta_temp[1:-1, 1:-1] = theta_new_interior
        # ==========================================
        # 2. CORNER-SAFE BOUNDARY CONDITIONS
        # ==========================================
        dZ_top = self.Z[:, -1] - self.Z[:, -2]
        R_top = self.R[:, -1]

        # Apply the stable Arc Flux
        gaussian_flux = Nu_base * np.exp(-(R_top / r_0)**2)

        # A. Top Surface (Flux). Applies to origin [0, -1] but stops before the right wall [:-1]
        theta_temp[:-1, -1] = (4.0 * theta_temp[:-1, -2] - theta_temp[:-1, -3] + 2.0 * dZ_top[:-1] * gaussian_flux[:-1]) / 3.0

        # B. Axis of Symmetry. Explicitly excludes corners [1:-1]
        theta_temp[0, 1:-1] = (4.0 * theta_temp[1, 1:-1] - theta_temp[2, 1:-1]) / 3.0

        # C. Crucible Walls: CONVECTIVE COOLING
        # Lowering Bi from 1.5 to 0.2 allows the ingot to retain more heat globally,
        # pushing the 1.0 melt contour deeper into the interior computational nodes.
        Bi = 0.2
        dR_right = self.R[-1, :] - self.R[-2, :]
        dZ_bottom = self.Z[:, 1] - self.Z[:, 0]

        # Right Wall cooling
        theta_temp[-1, :] = (4.0 * theta_temp[-2, :] - theta_temp[-3, :]) / (3.0 + 2.0 * dR_right * Bi)

        # Bottom Wall cooling
        theta_temp[:, 0] = (4.0 * theta_temp[:, 1] - theta_temp[:, 2]) / (3.0 + 2.0 * dZ_bottom * Bi)

        # 3. Global Under-relaxation
        self.relax_theta = 0.005
        self.theta = (1.0 - self.relax_theta) * theta_old + self.relax_theta * theta_temp

        return np.max(np.abs(self.theta - theta_old))


    def solve_magnetic_induction(self, Rm=10.0):
        """
        Solves the magnetic induction equation to determine the magnetic field 
        distribution (F) within the domain. Driven by the magnetic Reynolds 
        number (Rm) and advected by the fluid velocities.
        """
        F_old = self.F.copy()
        F_temp = self.F.copy() # Temporary array

        a = self.alpha
        g = self.gamma
        b = self.beta
        j2 = self.J**2

        u = self.U_comp[1:-1, 1:-1]
        w = self.W_comp[1:-1, 1:-1]

        f_plus_xi   = self.F[2:, 1:-1]
        f_minus_xi  = self.F[:-2, 1:-1]
        f_plus_eta  = self.F[1:-1, 2:]
        f_minus_eta = self.F[1:-1, :-2]

        f_cross = (self.F[2:, 2:] - self.F[:-2, 2:] -
                   self.F[2:, :-2] + self.F[:-2, :-2]) / 4.0

        advection = (u * (f_plus_xi - f_minus_xi) + w * (f_plus_eta - f_minus_eta)) / 2.0
        diffusion = (a * (f_plus_xi + f_minus_xi) + g * (f_plus_eta + f_minus_eta) - 2.0 * b * f_cross)

        rhs = diffusion - (Rm * advection * j2)

        # Stage Interior
        F_temp[1:-1, 1:-1] = rhs / (2.0 * (a + g))

        # Stage Boundaries (Assuming zero outer walls, symmetric axis)
        F_temp[0, :] = (4.0 * F_temp[1, :] - F_temp[2, :]) / 3.0
        F_temp[-1, :] = 0.0
        F_temp[:, 0] = 0.0
        F_temp[:, -1] = 0.0 # Modify this top boundary if your arc injects F here!

        # Under-relax the entire field to squash oscillations
        relax_F = 0.01
        self.F = (1.0 - relax_F) * F_old + relax_F * F_temp

        return np.max(np.abs(self.F - F_old))


    def solve_stream_function(self):
        """
        Solves the Poisson equation linking vorticity to the stream function.
        Applies a mathematical damping factor (liquid_fraction) to ensure 
        fluid velocity is strictly zero within the solid regions of the domain.
        """
        psi_old = self.psi.copy()

        a = self.alpha
        g = self.gamma
        b = self.beta
        j2 = self.J**2

        p_plus_xi   = self.psi[2:, 1:-1]
        p_minus_xi  = self.psi[:-2, 1:-1]
        p_plus_eta  = self.psi[1:-1, 2:]
        p_minus_eta = self.psi[1:-1, :-2]

        p_cross = (self.psi[2:, 2:] - self.psi[:-2, 2:] -
                   self.psi[2:, :-2] + self.psi[:-2, :-2]) / 4.0

        laplacian = (a * (p_plus_xi + p_minus_xi) + g * (p_plus_eta + p_minus_eta) - 2.0 * b * p_cross)

        # Calculate the raw stream function
        rhs = laplacian + (self.omega[1:-1, 1:-1] * j2)
        psi_new = rhs / (2.0 * (a + g))

        # ==========================================
        # NEW: Smooth Phase Change (Mushy Zone)
        # ==========================================
        # Ramps from 0.0 (fully solid at T=0.95) to 1.0 (fully liquid at T=1.05)
        liquid_fraction = np.clip((self.theta[1:-1, 1:-1] - 0.85) / 0.15, 0.0, 1.0)

        # Smoothly damp the flow in the solid/mushy region
        psi_new = psi_new * liquid_fraction

        self.psi[1:-1, 1:-1] = psi_new

        # Re-enforce Physical Boundaries
        self.psi[0, :] = 0.0
        self.psi[-1, :] = 0.0
        self.psi[:, 0] = 0.0
        self.psi[:, -1] = 0.0

        return np.max(np.abs(self.psi - psi_old))


    def solve_vorticity(self, Gr=2.9e8, Ma=5000.0):
        """
        Solves the fluid vorticity transport equation. 
        Couples thermal buoyancy (Grashof number), Marangoni surface shear 
        stress at the top boundary, and electromagnetic Lorentz forces. 
        Strictly confines flow generation to the liquid and mushy zones.
        """
        omega_old = self.omega.copy()

        a = self.alpha
        g = self.gamma
        b = self.beta
        j2 = self.J**2

        u = self.U_comp[1:-1, 1:-1]
        w = self.W_comp[1:-1, 1:-1]

        w_plus_xi   = self.omega[2:, 1:-1]
        w_minus_xi  = self.omega[:-2, 1:-1]
        w_plus_eta  = self.omega[1:-1, 2:]
        w_minus_eta = self.omega[1:-1, :-2]

        w_cross = (self.omega[2:, 2:] - self.omega[:-2, 2:] -
                   self.omega[2:, :-2] + self.omega[:-2, :-2]) / 4.0

        convection = (u * (w_plus_xi - w_minus_xi) + w * (w_plus_eta - w_minus_eta)) / 2.0

        # ==========================================
        # NEW: Smooth Phase Change (Mushy Zone)
        # ==========================================
        liquid_fraction = np.clip((self.theta[1:-1, 1:-1] - 0.85) / 0.15, 0.0, 1.0)

        # Only allow buoyancy to generate force in the liquid/mushy zone
        buoyancy = Gr * (self.theta[2:, 1:-1] - self.theta[:-2, 1:-1]) / 2.0
        buoyancy = buoyancy * liquid_fraction


        # ==========================================
        # NEW: LORENTZ FORCE (Electromagnetic Coupling)
        # ==========================================
        # 1. Calculate spatial derivatives of F^2 in computational space
        F2_xi  = (self.F[2:, 1:-1]**2 - self.F[:-2, 1:-1]**2) / 2.0
        F2_eta = (self.F[1:-1, 2:]**2 - self.F[1:-1, :-2]**2) / 2.0

        R_xi   = (self.R[2:, 1:-1] - self.R[:-2, 1:-1]) / 2.0
        R_eta  = (self.R[1:-1, 2:] - self.R[1:-1, :-2]) / 2.0

        # 2. Chain rule to find the true vertical gradient: d(F^2)/dZ
        dF2_dZ = (R_xi * F2_eta - R_eta * F2_xi) / self.J

        R_interior = self.R[1:-1, 1:-1] + 1e-10

        # 3. S_lor is the Electromagnetic Interaction Parameter
        # (Check your thesis for the exact Stuart number, but 50.0 is a safe start)
        S_lor = 50.0
        lorentz = -S_lor * (dF2_dZ / (R_interior**2))

        # 4. Multiply by the liquid fraction so magnetic fields don't push solid metal!
        lorentz = lorentz * liquid_fraction



        diffusion = (a * (w_plus_xi + w_minus_xi) + g * (w_plus_eta + w_minus_eta) - 2.0 * b * w_cross)

        # Calculate the raw vorticity
        rhs = diffusion - (convection * j2) - (buoyancy * j2) + (lorentz * j2)
        omega_new = rhs / (2.0 * (a + g))




        # Smoothly damp the vorticity in the solid region
        omega_new = omega_new * liquid_fraction

        self.omega[1:-1, 1:-1] = omega_new

        # ==========================================
        # NEW: MARANGONI BOUNDARY CONDITION (The Engine)
        # ==========================================
        # Vorticity at the top surface is driven by the radial temperature gradient.
        # Ma is the Marangoni number. (1000.0 to 5000.0 is a typical starting range)
        Ma = 5000.0

        # Calculate dTheta / dR along the top surface, excluding corners
        dR_top = self.R[2:, -1] - self.R[:-2, -1]
        dTheta_top = self.theta[2:, -1] - self.theta[:-2, -1]

        # Apply the Marangoni shear stress to the top boundary
        self.omega[1:-1, -1] = -Ma * (dTheta_top / dR_top)

        # Ensure other boundaries remain solid/no-slip
        self.omega[:, 0] = 0.0
        self.omega[0, :] = 0.0
        self.omega[-1, :] = 0.0


        return np.max(np.abs(self.omega - omega_old))



    def update_velocities(self):
        """Derives physical velocities from the stream function."""
        dpsi_dz = (self.psi[3:-1, 2:-2] - self.psi[1:-3, 2:-2]) / 2.0
        dpsi_dr = (self.psi[2:-2, 3:-1] - self.psi[2:-2, 1:-3]) / 2.0

        R_vals = self.R[2:-2, 2:-2] + 1e-10 # Prevent division by zero
        self.U_comp[2:-2, 2:-2] = dpsi_dz / R_vals
        self.W_comp[2:-2, 2:-2] = -dpsi_dr / R_vals

    def sync_boundaries(self):
        """
        Enforces zero-gradient (Neumann) boundary conditions using
        2nd-order one-sided finite differences.
        """
        fields = [self.theta, self.F, self.omega, self.U_comp, self.W_comp]

        for field in fields:
            # Left / Right boundaries (axis and crucible wall)
            # field[:, 1] = (4 * field[:, 2] - field[:, 3]) / 3.0
            field[:, 1] = (4.0 * field[:, 2] - field[:, 3]) / 3.0
            field[:, -2] = (4.0 * field[:, -3] - field[:, -4]) / 3.0

            # Top / Bottom boundaries (melt surface and solidus pool base)
            field[1, :] = (4.0 * field[2, :] - field[3, :]) / 3.0
            field[-2, :] = (4.0 * field[-3, :] - field[-4, :]) / 3.0

        # Note: The stream function (self.psi) does NOT get this update,
        # because it uses a Dirichlet condition (psi = 0) at all walls!

    def step(self):
        W = self.compute_flow_monitor()
        R_new, Z_new = self.grid_solver.relax_grid(self.R, self.Z, W, inner_iters=300, omega_sor=0.2, adapt_weight=25.0)
        self.remap_physics(R_new, Z_new)

    def remap_physics(self, R_new, Z_new):
        """
        Interpolates the active physical fields (temperature, vorticity, etc.) 
        from the old grid coordinates to the newly adapted Brackbill-Saltzman 
        grid coordinates using linear interpolation with a nearest-neighbor fallback.
        """
        old_points = np.column_stack((self.R.flatten(), self.Z.flatten()))
        new_points = np.column_stack((R_new.flatten(), Z_new.flatten()))

        fields = [self.psi, self.omega, self.theta, self.F]
        new_fields = []

        for field in fields:
            # 1. Primary linear interpolation
            interp_field = griddata(old_points, field.flatten(), new_points, method='linear')

            # 2. Nearest-neighbor fallback for floating-point boundary NaNs
            nan_mask = np.isnan(interp_field)
            if np.any(nan_mask):
                nearest_field = griddata(old_points, field.flatten(), new_points, method='nearest')
                interp_field[nan_mask] = nearest_field[nan_mask]

            new_fields.append(interp_field.reshape(self.m, self.n))

        self.psi, self.omega, self.theta, self.F = new_fields
        self.R, self.Z = R_new, Z_new

    def save_checkpoint(self, filename):
        full_path = f"/content/drive/MyDrive/VAR_Simulations/{filename}"
        np.savez_compressed(full_path, R=self.R, Z=self.Z, theta=self.theta, omega=self.omega, psi=self.psi, F=self.F)
        print(f"--> Saved to {filename}")

    def load_checkpoint(self, filename):
        full_path = f"/content/drive/MyDrive/VAR_Simulations/{filename}"
        data = np.load(full_path)
        self.R, self.Z = data['R'], data['Z']
        self.theta, self.omega = data['theta'], data['omega']
        self.psi, self.F = data['psi'], data['F']
        print(f"--> Loaded {filename}")

    #def compute_flow_monitor(self):
    #    return 1.0 + 50.0 * (np.abs(self.omega) / (np.max(np.abs(self.omega)) + 1e-9))

    def compute_flow_monitor(self):
        # Normalize vorticity
        norm_omega = np.abs(self.omega) / (np.max(np.abs(self.omega)) + 1e-9)

        # Normalize temperature gradients (to track the melt pool boundary)
        dTheta_R = np.gradient(self.theta, axis=0)
        dTheta_Z = np.gradient(self.theta, axis=1)
        grad_theta = np.sqrt(dTheta_R**2 + dTheta_Z**2)
        norm_grad_theta = grad_theta / (np.max(grad_theta) + 1e-9)

        # Blend them together! (50x gravity for flow, 50x gravity for the melt boundary)
        return 1.0 + 50.0 * norm_omega + 50.0 * norm_grad_theta



def plot_simulation_state(sim, step):
    """Visualizes the grid structure and thermal profile."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Plot Grid
    ax1.scatter(sim.R, sim.Z, s=0.5, c='k', alpha=0.3)
    ax1.set_title(f"Step {step}: Adaptive Grid Distribution")

    # Plot Temperature
    contour = ax2.contourf(sim.R, sim.Z, sim.theta, levels=50, cmap='inferno')
    fig.colorbar(contour, ax=ax2)
    ax2.set_title(f"Step {step}: Thermal Profile")

    plt.show()

def plot_simulation_state(sim, step):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # Grid structure
    ax1.scatter(sim.R, sim.Z, s=0.5, c='k', alpha=0.2)
    ax1.set_title(f"Step {step}: Adaptive Grid Distribution")

    # Physics field (Temperature)
    contour = ax2.contourf(sim.R, sim.Z, sim.theta, levels=60, cmap='magma')
    fig.colorbar(contour, ax=ax2, label='Temperature')
    ax2.set_title(f"Step {step}: Thermal Profile")

if __name__ == "__main__":
    # 1. Initialize and Load
    sim = VAR_Adaptive_Simulator(m=100, n=100)
    #sim.load_checkpoint("var_sim_step_8_fully_converged.npz")
    #sim.load_checkpoint("var_sim_checkpoint_step_99900.npz")
    sim.load_checkpoint("var_sim_FINAL_step_149999.npz")

    # 2. Failsafe Initialization
    # (Ensure velocity fields exist in case they weren't in the checkpoint)
    if not hasattr(sim, 'U_comp'):
        sim.U_comp = np.zeros((100, 100))
        sim.W_comp = np.zeros((100, 100))


    # ==========================================
    # 2. THE HOT START PATCH (Bypass the Glass Floor)
    # ==========================================
    r_pool = 0.5   # Width of our starting liquid bowl
    z_top = 1.5    # The Z-coordinate of your top surface

    # Calculate the geometric distance from the top-center arc origin
    dist_from_arc = np.sqrt(sim.R**2 + (sim.Z - z_top)**2)

    # Force the interior of this bowl to be fully melted (T = 1.1)
    sim.theta[dist_from_arc < r_pool] = np.maximum(sim.theta[dist_from_arc < r_pool], 1.1)

    print("Hot Start applied! Liquid pool initialized.")



    # Calculate baseline metrics before any solvers run
    sim.calculate_grid_metrics()

    # Diagnostic prints
    print("\n--- DEBUG DIAGNOSTICS ---")
    print(f"Loaded Grid Shape (R): {sim.R.shape}")
    print(f"Loaded Grid Shape (Z): {sim.Z.shape}")
    print(f"Min J right after loading: {np.min(sim.J):.6f}")
    print(f"Max J right after loading: {np.max(sim.J):.6f}")
    print("-------------------------\n")

    if np.min(sim.J) < 0:
        print(f"\n[ALIGNMENT] Loaded grid has negative Jacobian (Min J: {np.min(sim.J):.4f}).")
        print("[ALIGNMENT] Axis mismatch detected. Transposing arrays to match adaptive coordinate system...")

        # Transpose all spatial and physical fields
        sim.R = sim.R.T
        sim.Z = sim.Z.T
        sim.theta = sim.theta.T
        sim.omega = sim.omega.T
        sim.psi = sim.psi.T
        sim.F = sim.F.T

        # Recalculate grid metrics with aligned data
        sim.calculate_grid_metrics()
        print(f"[ALIGNMENT] Transpose successful. New Min J: {np.min(sim.J):.6f}\n")

    # 3. Parameters
    max_total_steps = 100000
    adapt_freq = 200
    save_freq = 50
    flow_tol = 1e-5
    # 3.5. Update the Arc Power! (Override whatever was used previously)
    # We must melt the metal to wake up the fluid flow.
    #current_Nu_base = 40.0

    # 4. Resume the Main Loop
    start_step = 140000
    max_total_steps = 150000

    # ==========================================
    # STEADY-STATE DETECTOR SETTINGS
    # ==========================================
    steady_tolerance = 1e-4       # The residual threshold for "equilibrium"
    required_steady_steps = 100   # How many consecutive steps it must stay below the tolerance
    steady_counter = 0            # Tracks consecutive stable steps

    print(f"Starting simulation. Steady-state target: {steady_tolerance} for {required_steady_steps} steps.")
    print("Starting High-Fidelity VAR Simulation Loop...")

    for step in range(start_step, max_total_steps):
        # Define how long the physics should stabilize before adapting again
        # 1000 steps gives the residual plenty of time to decay smoothly
        adapt_epoch_interval = 1000

        # ==========================================
        # A. GRID ADAPTATION (Multi-Epoch Refinement)
        # ==========================================
        # We also trigger on step 10 just to get the initial cluster formed quickly
        if step == 10 or (step > 0 and step % adapt_epoch_interval == 0):
            epoch_num = step // adapt_epoch_interval if step >= adapt_epoch_interval else 0
            print("\n" + "="*60)
            print(f"[GRID EPOCH {epoch_num}] Output stabilized. Pausing physics.")
            print(f"[GRID EPOCH {epoch_num}] Re-calculating flow monitor and adapting grid...")

            # 1. Compute the new monitor based on the stabilized physics
            W_monitor = sim.compute_flow_monitor()

            # 2. Relax the grid using our smoothed, proportional-limited solver
            # (Assuming you initialized your BrackbillSaltzmanSolver as 'solver' earlier)
            new_R, new_Z = sim.grid_solver.relax_grid(sim.R, sim.Z, W_monitor)

            # 3. Safely map the physics over using our bulletproof remap function
            sim.remap_physics(new_R, new_Z)

            # 4. Update the metric derivatives for the PDE solvers
            sim.calculate_grid_metrics()

            print(f"[GRID EPOCH {epoch_num}] Mesh successfully refined. Resuming solvers...")
            print("="*60 + "\n")

        # B. PHYSICS SOLVERS (Coupled Order)
        res_F = sim.solve_magnetic_induction()
        res_T = sim.solve_temperature()
        # Print the maximum temperature of the INTERIOR nodes (ignoring the surface)
        #print(f"Max Interior Temp: {np.max(sim.theta[1:-1, 1:-1]):.4f}")
        res_W = sim.solve_vorticity()
        res_P = sim.solve_stream_function()

        # ==========================================

        # C. KINEMATIC & BOUNDARY UPDATES
        sim.update_velocities()
        sim.sync_boundaries()

        # D. MONITORING
        #max_res = max(res_F, res_T, res_W, res_P)
        # ==========================================
        # STEADY-STATE CHECK
        # ==========================================
        # We only care about the fluid and thermal fields for the physical pool shape
        max_res = max(res_T, res_W, res_P)

        if max_res < steady_tolerance:
            steady_counter += 1
        else:
            steady_counter = 0  # Reset if it bounces back up!

        # Check if we hit the target
        if steady_counter >= required_steady_steps:
            print(f"\n============================================================")
            print(f"[STEADY STATE REACHED] Simulation converged at Step {step}!")
            print(f"All primary physics residuals held below {steady_tolerance} for {required_steady_steps} steps.")
            print(f"============================================================\n")

            # Force a final save and break the loop
            np.savez(f"var_sim_FINAL_steady_state_step_{step}.npz",
                 theta=sim.theta, psi=sim.psi, omega=sim.omega, F=sim.F,
                 R=sim.R, Z=sim.Z)
            break
        # ==========================================

        # FAILSAFE: Catch Blowups instantly
        if np.isnan(max_res) or np.min(sim.J) < 0:
            print(f"\nCRITICAL FAILURE at Step {step}: Simulation blew up or grid folded (Min J: {np.min(sim.J):.4f}).")
            break

        if step % 10 == 0:
            print(f"Step {step:03d} | Max Res: {max_res:.2e} | Min J: {np.min(sim.J):.2e} "
                      f"(F:{res_F:.6f}, Th:{res_T:.6f}, W:{res_W:.6f}, Psi:{res_P:.6f})")

            # Plotting and memory management
            try:
                plot_simulation_state(sim, step)
                import matplotlib.pyplot as plt
                plt.close('all') # This fixes the memory leak!
            except NameError:
                pass

        # E. PERSISTENCE
        if step > 0 and step % save_freq == 0:
            sim.save_checkpoint(f"var_sim_checkpoint_step_{step}.npz")

        # F. CONVERGENCE
        if max_res < flow_tol:
            print(f"Convergence reached at step {step}!")
            break

    # ==========================================
    # POST-LOOP WRAP-UP (OUTSIDE THE FOR LOOP)
    # ==========================================
    print("Simulation loop concluded. Saving final state...")

    # 1. Guarantee a final save
    sim.save_checkpoint(f"var_sim_FINAL_step_{step}.npz")

    # 2. Guarantee a final plot
    try:
        plot_simulation_state(sim, step)
        plt.savefig(f"/content/drive/MyDrive/VAR_Simulations/snapshot_FINAL.png")
        print("Final plot generated and saved.")
    except NameError:
        print("Note: plot_simulation_state function not defined. Skipping final plot.")

    # 3. Grid Health Check
    if np.min(sim.J) <= 1e-5:
        print("WARNING: Grid collapse detected (J <= 0). The mesh has folded.")
        print("Recommendation: Decrease your flow monitor weight (W) or increase grid smoothing (omega_sor) in BrackbillSaltzmanSolver.")

    # ==========================================
    # FINAL CONCLUSION & PLOTTING
    # ==========================================
    print("\nSimulation Loop Terminated. Generating Final High-Fidelity Plots...")

    import matplotlib.pyplot as plt

    # Create a wide 3-panel figure to show the full picture
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(22, 6))

    # --- Panel 1: The Final Adapted Mesh ---
    # Plotting the actual grid lines to show the density around the vortex

    # Clear the axis first
    ax1.clear()

    # Plot the horizontal grid lines
    for i in range(sim.R.shape[0]):
        ax1.plot(sim.R[i, :], sim.Z[i, :], color='black', linewidth=0.5, alpha=0.3)

    # Plot the vertical grid lines
    for j in range(sim.R.shape[1]):
        ax1.plot(sim.R[:, j], sim.Z[:, j], color='black', linewidth=0.5, alpha=0.3)

    ax1.set_title("Final Adapted Grid")

    #for i in range(sim.m):
    #    ax1.plot(sim.R[i, :], sim.Z[i, :], 'k-', alpha=0.1, linewidth=0.5)
    #for j in range(sim.n):
    #    ax1.plot(sim.R[:, j], sim.Z[:, j], 'k-', alpha=0.1, linewidth=0.5)
    #ax1.set_title("Final Adapted Grid Distribution")
    #ax1.set_xlabel("Radius (R)")
    #ax1.set_ylabel("Depth (Z)")
    #ax1.set_xlim([0, 1.0])
    #ax1.set_ylim([0, 1.5])

    # --- Panel 2: The Final Thermal Profile ---
    t_plot = ax2.contourf(sim.R, sim.Z, sim.theta, levels=100, cmap='inferno')
    fig.colorbar(t_plot, ax=ax2, label='Non-Dimensional Temperature')
    ax2.set_title("Final Resolved Thermal Field")
    ax2.set_xlabel("Radius (R)")

    # --- Panel 3: Flow Streamlines (Vortex Flow) ---
    # Using 80+ contour levels guarantees we capture the tightly packed recirculation
    psi_plot = ax3.contour(sim.R, sim.Z, sim.psi, levels=80, cmap='viridis', linewidths=1.2)
    fig.colorbar(psi_plot, ax=ax3, label='Stream Function')
    ax3.set_title("Final Fluid Streamlines")
    ax3.set_xlabel("Radius (R)")

    plt.tight_layout()
    plt.savefig("var_sim_FINAL_CONVERGED_RESULTS.png", dpi=300, bbox_inches='tight')
    plt.show()

    print("--> Final plots saved to var_sim_FINAL_CONVERGED_RESULTS.png")