import tkinter as tk
from tkinter import ttk, messagebox
import math

class BridgeFeasibilityCalc:
    def __init__(self, root):
        self.root = root
        self.root.title("HV Cable Feasibility & Force Calculator (Detailed)")
        self.root.geometry("900x900")#what does this do exactly?
        
        # --- TAB CONTROL ---
        tab_control = ttk.Notebook(root)
        self.tab_bridge = ttk.Frame(tab_control)
        self.tab_sc = ttk.Frame(tab_control)
        
        tab_control.add(self.tab_bridge, text="1. Bridge Feasibility (CIGRE 669)")
        tab_control.add(self.tab_sc, text="2. Short Circuit Forces (IEC 61914)")
        tab_control.pack(expand=1, fill="both")
        
        self.setup_bridge_tab()
        self.setup_sc_tab()

    def setup_bridge_tab(self):
        """
        Sets up the GUI for the Thermal Expansion and Geometric Feasibility Check.
        References: CIGRE Technical Brochure 669 (Mechanical forces in large cross section cable systems).
        """
        frame = ttk.LabelFrame(self.tab_bridge, text="Thermal & Geometric Parameters (Ref: CIGRE TB 669)", padding=15)
        frame.pack(fill="x", padx=10, pady=10)

        # --- ROW 0: Conductor Area ---
        # Used for Axial Thrust Calculation: F = E * A * alpha * dT
        ttk.Label(frame, text="Conductor Area [mm²]:").grid(row=0, column=0, sticky="w")
        self.area_var = tk.DoubleVar(value=2500)
        ttk.Entry(frame, textvariable=self.area_var).grid(row=0, column=1)

        # --- ROW 1: Cable Diameter ---
        # Used for:
        # 1. Geometric Snake Width (Space check)
        # 2. Minimum bending radius for Radial Force check
        ttk.Label(frame, text="Cable Outer Diameter (De) [mm]:").grid(row=1, column=0, sticky="w")
        self.dia_var = tk.DoubleVar(value=130)
        self.dia_entry = ttk.Entry(frame, textvariable=self.dia_var)
        self.dia_entry.grid(row=1, column=1)
        # Automatically update the spacing in Tab 2 when diameter changes
        self.dia_entry.bind("<FocusOut>", self.update_default_spacing)

        # --- ROW 2: Young's Modulus (E) ---
        # Reference: CIGRE 669, Page 100, Table 7 "Recommended values for Effective Modulus"
        # Typical values: 20,000 - 35,000 N/mm² for large XLPE cables.
        ttk.Label(frame, text="Young's Modulus (E) [N/mm²]:").grid(row=2, column=0, sticky="w")
        self.E_var = tk.DoubleVar(value=30000)
        ttk.Entry(frame, textvariable=self.E_var).grid(row=2, column=1)
        ttk.Label(frame, text="(Ref: CIGRE 669 Table 7)").grid(row=2, column=2, sticky="w", padx=5)

        # --- ROW 3: Thermal Expansion (Alpha) ---
        # Reference: CIGRE 669, Page 100, Table 6 "Coefficients of thermal expansion"
        # Copper: 18e-6 1/K (Standard Engineering Value for Cable Systems)
        # Aluminum: 23e-6 1/K
        ttk.Label(frame, text="Expansion Coeff (α) [x10⁻⁶ 1/K]:").grid(row=3, column=0, sticky="w")
        self.alpha_input_var = tk.DoubleVar(value=20)#Could also be 18 if you dont want to be conservative
        ttk.Entry(frame, textvariable=self.alpha_input_var).grid(row=3, column=1)
        ttk.Label(frame, text="(Ref: CIGRE 669 Table 6)").grid(row=3, column=2, sticky="w", padx=5)

        # --- ROW 4: Temperature Delta ---
        # Definition: T_max_operation - T_installation (e.g., 90°C - 10°C = 80K)
        # Reference: CIGRE 669 Page 13 (Formula definitions)
        ttk.Label(frame, text="Temp. Difference (ΔT) [K]:").grid(row=4, column=0, sticky="w")
        self.dt_var = tk.DoubleVar(value=70)#Installation at 20C° and max conducter temp at 90C° differance is 20K°
        ttk.Entry(frame, textvariable=self.dt_var).grid(row=4, column=1)

        # --- ROW 5: Span Length ---
        # Distance between cleats. Used for "Snaking Amplitude" calculation.
        ttk.Label(frame, text="Span / Cleat Distance (L) [mm]:").grid(row=5, column=0, sticky="w")
        self.span_var = tk.DoubleVar(value=5000)#every 5m there is a cleate
        ttk.Entry(frame, textvariable=self.span_var).grid(row=5, column=1)

        # --- ROW 6: Duct Width ---
        # User constraint for geometric check.
        ttk.Label(frame, text="Available Duct Width [mm]:").grid(row=6, column=0, sticky="w")
        self.width_constraint_var = tk.DoubleVar(value=600)
        ttk.Entry(frame, textvariable=self.width_constraint_var).grid(row=6, column=1)

        # --- ROW 7: Bending Radius (NEW) ---
        # Used for Radial Force Check at the 90 degree corner
        ttk.Label(frame, text="Bend Radius at Corner (R) [mm]:").grid(row=7, column=0, sticky="w")
        self.bend_radius_var = tk.DoubleVar(value=2000) # Default placeholder
        ttk.Entry(frame, textvariable=self.bend_radius_var).grid(row=7, column=1)
        ttk.Label(frame, text="(Standard Min: 15 * OD)").grid(row=7, column=2, sticky="w", padx=5)

        # --- ROW 8: Number of Cleats at Bend (NEW) ---
        ttk.Label(frame, text="Number of Cleats at 90° Bend:").grid(row=8, column=0, sticky="w")
        self.num_cleats_bend_var = tk.IntVar(value=5) # Default 5
        ttk.Entry(frame, textvariable=self.num_cleats_bend_var).grid(row=8, column=1)
        ttk.Label(frame, text="(To distribute radial load)").grid(row=8, column=2, sticky="w", padx=5)

        # --- BUTTON ---
        btn = ttk.Button(frame, text="CALCULATE BRIDGE FEASIBILITY", command=self.calc_bridge)
        btn.grid(row=9, column=0, columnspan=3, pady=15)

        # --- RESULTS AREA ---
        self.bridge_res_frame = ttk.LabelFrame(self.tab_bridge, text="Analysis Results", padding=15)
        self.bridge_res_frame.pack(fill="both", expand=True, padx=10, pady=10)
        
        self.lbl_rigid = ttk.Label(self.bridge_res_frame, text="Result A (Rigid)...", font=("Courier", 10), justify="left")
        self.lbl_rigid.pack(anchor="w", pady=5)
        
        ttk.Separator(self.bridge_res_frame, orient='horizontal').pack(fill='x', pady=10)
        
        self.lbl_flex = ttk.Label(self.bridge_res_frame, text="Result B (Flexible)...", font=("Courier", 10), justify="left")
        self.lbl_flex.pack(anchor="w", pady=5)

    def setup_sc_tab(self):
        """
        Sets up the GUI for Short Circuit Force Calculation.
        References: IEC 61914 Annex B and DNV Training Material.
        """
        frame = ttk.LabelFrame(self.tab_sc, text="Short Circuit Parameters (Ref: IEC 61914)", padding=15)
        frame.pack(fill="x", padx=10, pady=10)
        
        # Configuration Type
        ttk.Label(frame, text="Configuration:").grid(row=0, column=0, sticky="w")
        self.sc_type_var = tk.StringVar(value="Flat")
        type_combo = ttk.Combobox(frame, textvariable=self.sc_type_var, state="readonly")
        type_combo['values'] = ("Trefoil", "Flat")
        type_combo.grid(row=0, column=1, sticky="w", pady=5)
        
        # Peak Current (ip)
        # Note: Formulas use ip (Peak), not Ik'' (RMS).
        # Reference: IEC 61914 Annex B formulas use 'ip'.
        ttk.Label(frame, text="Peak Current (ip) [kA]:").grid(row=1, column=0, sticky="w")
        self.ip_var = tk.DoubleVar(value=160)
        ttk.Entry(frame, textvariable=self.ip_var).grid(row=1, column=1, pady=5)
        
        # Spacing (S)
        # Reference: IEC 61914 Annex B defines S as "centre-to-centre distance".
        ttk.Label(frame, text="Conductor Center-to-Center Spacing (S) [mm]:").grid(row=2, column=0, sticky="w")
        self.sc_spacing_var = tk.DoubleVar(value=130) # Default to diameter
        ttk.Entry(frame, textvariable=self.sc_spacing_var).grid(row=2, column=1, pady=5)
        
        ttk.Label(frame, text="(Distance between the centers of the cables)", font=("Arial", 8, "italic")).grid(row=3, column=1, sticky="w")

        # Cleat Distance Info
        ttk.Label(frame, text="Cleat Distance (L): Uses 'Span' from Tab 1").grid(row=4, column=0, columnspan=2, pady=15, sticky="w")

        btn = ttk.Button(frame, text="CALCULATE SHORT CIRCUIT FORCES", command=self.calc_sc)
        btn.grid(row=5, column=0, columnspan=2, pady=10)
        
        # Result Label Area
        self.sc_res_label = ttk.Label(self.tab_sc, text="Results will appear here...", font=("Courier", 10), justify="left")
        self.sc_res_label.pack(pady=20, padx=20, anchor="w")

    def update_default_spacing(self, event):
        """Updates default spacing and radius when diameter changes in Tab 1."""
        try:
            d = self.dia_var.get()
            self.sc_spacing_var.set(d)           # Update Short Circuit Spacing
            self.bend_radius_var.set(15 * d)     # Update Bend Radius to Min Standard
        except:
            pass

    def calc_bridge(self):
        """
        Performs the Thermal and Geometric Feasibility Calculation.
        
        Physics & Equations:
        1. Axial Thrust (Rigid): CIGRE 669 Page 13, Eq 3: F = E * A * alpha * dT
        2. Radial Force (Bend): Basic Mechanics: F_radial = F_axial / Radius
        3. Flexible Snake: CIGRE 669 Section 4.2 (Approximation of snake width)
        """
        try:
            # --- INPUTS ---
            area = self.area_var.get()          # mm²
            dt = self.dt_var.get()              # K
            od = self.dia_var.get()             # mm
            span = self.span_var.get()          # mm
            avail_width = self.width_constraint_var.get() # mm
            
            E_eff = self.E_var.get()            # N/mm²
            alpha = self.alpha_input_var.get() * 1e-6 # 1/K

            # --- A. RIGID INSTALLATION CALCULATION ---
            # Source: CIGRE TB 669, Equation 3, Page 13
            # F_th = E_eff * Ac * alpha * dT
            F_rigid_N = E_eff * area * alpha * dt
            F_rigid_kN = F_rigid_N / 1000.0
            
            # Total Thrust (3 Phases)
            # Assumption: All 3 phases are anchored to the same abutment.
            F_total_tons = (F_rigid_kN * 3) / 9.81 

            # Radial Force at 90° Bend
            # Source: Mechanics of curved beams under axial load.
            # Approximation: Distributed Radial Force = Axial Force / Bend Radius
            # --- 1. Get Bend Radius Input ---
            bend_radius_mm = self.bend_radius_var.get()
            
            # --- 2. Check Min Radius Warning (15 * OD) ---
            min_radius = 15 * od
            if bend_radius_mm < min_radius:
                messagebox.showwarning(
                    "Bending Radius Warning", 
                    f"Selected Radius ({bend_radius_mm:.0f} mm) is less than the standard minimum "
                    f"of 15x Cable Diameter ({min_radius:.0f} mm).\n\n"
                    "Risk: Damage to cable insulation or sheath."
                )

            # --- 3. Force per Cleat (User Input) ---
            num_cleats_bend = self.num_cleats_bend_var.get()
            
            # Conversion to meters
            bend_radius_m = bend_radius_mm / 1000.0

            # Validation
            if num_cleats_bend <= 0:
                messagebox.showerror("Input Error", "Number of cleats at bend must be at least 1.")
                return

            # A. Calculate Distributed Radial Force (q) in kN/m
            # This is the "pressure" trying to straighten the cable
            q_radial_knm = F_rigid_kN / bend_radius_m

            # B. Calculate Arc Length of the 90 degree bend
            arc_length_m = (2 * math.pi * bend_radius_m) / 4.0

            # C. Calculate Pull-Out Force per Cleat (Point Load)
            # Total Force on Arc = q * Length. Divide by count.
            F_cleat_bend_kN = (q_radial_knm * arc_length_m) / num_cleats_bend

            rigid_text = "--- OPTION A: RIGID INSTALLATION (Straight) ---\n"
            rigid_text += f"1. AXIAL THRUST (Push from Bridge):\n"
            rigid_text += f"   Force per Cable: {F_rigid_kN:.1f} kN\n"
            rigid_text += f"   TOTAL LOAD on Abutment: {F_total_tons:.1f} Tons\n\n"
            rigid_text += f"2. 90° BEND 'KILLER' LOAD:\n"
            rigid_text += f"   (Assuming R={bend_radius_m:.1f}m and {num_cleats_bend} cleats on corner)\n"
            rigid_text += f"   > Distributed Radial Force (q): {q_radial_knm:.1f} kN/m\n"
            rigid_text += f"   > PULL-OUT FORCE PER CLEAT: {F_cleat_bend_kN:.1f} kN\n"
            rigid_text += "   (NOTE: This is continuous tension, not short-circuit shear!)"
            
            if F_total_tons > 30:
                rigid_text += "VERDICT: IMPOSSIBLE. Structure likely cannot hold >30 tons shear."
            else:
                rigid_text += "VERDICT: High load, check structural limits."
            
            self.lbl_rigid.config(text=rigid_text)

            # --- B. FLEXIBLE INSTALLATION CALCULATION ---
            # Source: CIGRE TB 669, Section 4.2 "Flexible Systems"
            # Step 1: Calculate Free Expansion (dL) over the span
            dL = alpha * span * dt 
            
            # Step 2: Calculate Snake Geometry
            # Heuristic approximation for sine-wave snake amplitude to absorb dL.
            # Initial Amplitude (A0) usually ~ 1.2 x OD
            A_initial = 1.2 * od 
            # Final Amplitude (A_final)
            A_final = math.sqrt(A_initial**2 + (0.375 * span * dL))
            
            # Step 3: Required Width
            # Width = Effective Cable Width (Trefoil ~2*OD) + Total Swing (2 * Amplitude)
            req_width = (2 * od) + (2 * A_final)
            
            flex_text = "--- OPTION B: FLEXIBLE INSTALLATION (Snaking) ---\n"
            flex_text += f"Geometric Calculation (Ref: CIGRE 669 Sec 4.2):\n"
            flex_text += f"  > Expansion over {span}mm span: {dL:.2f} mm\n"
            flex_text += f"  > Required Snake Amplitude: {A_final:.0f} mm\n"
            flex_text += f"  > MINIMUM DUCT WIDTH REQUIRED: {req_width:.0f} mm\n"
            
            if req_width > avail_width:
                flex_text += f"VERDICT: IMPOSSIBLE. Needs {req_width:.0f} mm, only {avail_width:.0f} available."
            else:
                flex_text += "VERDICT: Fits geometrically. Requires Swivel Cleats."

            self.lbl_flex.config(text=flex_text)

        except ValueError:
            messagebox.showerror("Error", "Check inputs in Tab 1.")

    def calc_sc(self):
        """
        Performs the Short Circuit Force Calculation.
        
        Physics & Equations:
        Reference: IEC 61914 Edition 2.0 2015-11, Annex B "Formulae for calculation of forces"
        
        Variables:
        - ip: Peak short-circuit current [kA]
        - S: Centre-to-centre spacing [m]
        - Ft: Maximum force per unit length [N/m]
        
        Equations implemented:
        1. 2-Phase Fault (Phase-to-Phase): Eq B.4 -> Factor 0.20
        2. 3-Phase Flat (Outer): Eq B.5           -> Factor 0.16
        3. 3-Phase Flat (Middle): Eq B.6          -> Factor 0.17
        4. 3-Phase Trefoil: Eq B.7                -> Factor 0.17
        """
        try:
            # Inputs
            ip = self.ip_var.get()          # kA
            od = self.dia_var.get()         # mm
            span = self.span_var.get()      # mm (Cleat spacing)
            config = self.sc_type_var.get() # "Flat" or "Trefoil"
            
            S_mm = self.sc_spacing_var.get()
            
            # Sanity Check
            if S_mm < od:
                messagebox.showwarning("Physics Warning", f"Spacing ({S_mm}mm) is less than Cable Diameter ({od}mm)!")

            S_m = S_mm / 1000.0   # Convert to meters
            span_m = span / 1000.0 # Convert to meters
            
            # --- CALCULATION LOGIC (IEC 61914 Annex B) ---
            
            # 1. Two-Phase Fault (Equation B.4)
            # "For a phase-to-phase short circuit... the maximum force..."
            # Formula: F = 0.2 * ip^2 / S
            f_2ph_nm = 0.2 * (ip**2) / S_m
            f_2ph_knm = f_2ph_nm / 1000.0
            
            # 2. Three-Phase Flat - Outer Conductor (Equation B.5)
            # "For a three-phase short circuit... on the outer conductors..."
            # Formula: F = 0.16 * ip^2 / S
            f_3ph_flat_outer_nm = 0.16 * (ip**2) / S_m
            f_3ph_flat_outer_knm = f_3ph_flat_outer_nm / 1000.0
            
            # 3. Three-Phase Flat - Middle Conductor (Equation B.6)
            # "For a three-phase short circuit... on the centre conductor..."
            # Formula: F = 0.17 * ip^2 / S
            f_3ph_flat_middle_nm = 0.17 * (ip**2) / S_m
            f_3ph_flat_middle_knm = f_3ph_flat_middle_nm / 1000.0
            
            # 4. Three-Phase Trefoil (Equation B.7)
            # "For a three-phase short circuit... in trefoil formation..."
            # Formula: F = 0.17 * ip^2 / S
            f_3ph_trefoil_nm = 0.17 * (ip**2) / S_m
            f_3ph_trefoil_knm = f_3ph_trefoil_nm / 1000.0
            
            # --- RESULT GENERATION ---
            res = f"--- RESULTS FOR {config.upper()} CONFIGURATION ---\n"
            res += f"Spacing S = {S_mm} mm | Peak Current ip = {ip} kA\n\n"
            
            # Always show 2-Phase (Worst Case usually)
            res += "1. 2-PHASE FAULT (Phase-to-Phase)\n"
            res += f"   Ref: IEC 61914 Eq B.4 (Factor 0.2)\n"
            res += f"   Force per Meter: {f_2ph_knm:.2f} kN/m\n"
            res += f"   Force on Cleat (L={span}mm): {(f_2ph_knm * span_m):.2f} kN\n\n"
            
            if config == "Flat":
                res += "2. 3-PHASE FAULT (Flat Formation)\n"
                res += "   A) Outer Conductor:\n"
                res += f"      Ref: IEC 61914 Eq B.5 (Factor 0.16)\n"
                res += f"      Force per Meter: {f_3ph_flat_outer_knm:.2f} kN/m\n"
                res += f"      Force on Cleat: {(f_3ph_flat_outer_knm * span_m):.2f} kN\n"
                res += "   B) Middle Conductor (Max 3-Ph):\n"
                res += f"      Ref: IEC 61914 Eq B.6 (Factor 0.17)\n"
                res += f"      Force per Meter: {f_3ph_flat_middle_knm:.2f} kN/m\n"
                res += f"      Force on Cleat: {(f_3ph_flat_middle_knm * span_m):.2f} kN\n"
            
            elif config == "Trefoil":
                res += "2. 3-PHASE FAULT (Trefoil Formation)\n"
                res += f"   Ref: IEC 61914 Eq B.7 (Factor 0.17)\n"
                res += f"   Force per Meter: {f_3ph_trefoil_knm:.2f} kN/m\n"
                res += f"   Force on Cleat: {(f_3ph_trefoil_knm * span_m):.2f} kN\n"

            res += "\n" + "="*40 + "\n"
            
            # Determine Max Design Force
            if config == "Flat":
                max_force_nm = max(f_2ph_nm, f_3ph_flat_middle_nm)
                design_case = "2-Phase (Eq B.4)" if f_2ph_nm > f_3ph_flat_middle_nm else "3-Phase Middle (Eq B.6)"
            else:
                max_force_nm = max(f_2ph_nm, f_3ph_trefoil_nm)
                design_case = "2-Phase (Eq B.4)" if f_2ph_nm > f_3ph_trefoil_nm else "3-Phase Trefoil (Eq B.7)"
            
            max_force_cleat = (max_force_nm / 1000.0) * span_m
            
            res += f"DESIGN LOAD (WORST CASE): {max_force_cleat:.2f} kN\n"
            res += f"Governing Case: {design_case}"
            
            if max_force_cleat > 25:
                 res += "\n\nWARNING: Load > 25 kN. Standard M10/M12 bolts may shear."

            # Update Label & Popup
            self.sc_res_label.config(text=res)
            
            # Popup for clearer visibility
            messagebox.showinfo(f"Detailed Results ({config})", res)
            
        except ValueError:
            messagebox.showerror("Error", "Check inputs in Tab 2.")

if __name__ == "__main__":
    root = tk.Tk()
    app = BridgeFeasibilityCalc(root)
    root.mainloop()
