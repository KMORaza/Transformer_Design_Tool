from PyQt5.QtWidgets import QDialog, QWidget, QHBoxLayout, QVBoxLayout, QFormLayout, QLabel, QLineEdit, QPushButton, QTextEdit, QComboBox, QSizePolicy
from PyQt5.QtGui import QFont
from PyQt5.QtCore import Qt
from transformer import Transformer
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import math, traceback

class VibrationAndAcousticNoiseModelWindow(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Vibration and Acoustic Noise Model")
        self.setGeometry(200, 200, 1000, 600)
        self.setFont(QFont("Consolas", 10))
        
        # Apply dark theme
        self.setStyleSheet("""
            QDialog {
                background-color: #2E2E2E;
            }
            QWidget#controlPanel, QWidget#outputPanel {
                background: qlineargradient(x1:0, y1:0, x2:2, y2:1, stop:0 #4A4A4A, stop:1 #3A3A3A);
                border: 2px solid #1E1e1E;
                border-radius: 5px;
                padding: 10px;
                box-shadow: 5px 5px 10px #1E1E1E;
            }
            QLabel {
                color: #D3D3D3;
                font: 10pt "Consolas";
            }
            QLineEdit, QComboBox, QTextEdit {
                background-color: #3C3C3C;
                color: #FFFFFF;
                border: 1px solid #5A5A5A;
                border-radius: 3px;
                padding: 3px;
                font: 10pt "Consolas";
            }
            QLineEdit, QComboBox {
                min-width: 150px;
            }
            QPushButton {
                background: qlineargradient(x1:0, y1:0, x2:1, y2:1, stop:0 #6A6A6A, stop:1 #5A5A5A);
                color: #FFFFFF;
                border: 1px solid #4A4A4A;
                border-radius: 5px;
                padding: 5px;
                font: 10pt "Consolas";
            }
            QPushButton:hover {
                background: qlineargradient(x1:0, y1:0, x2:1, y2:1, stop:0 #7A7A7A, stop:1 #6A6A6A);
            }
            QPushButton:pressed {
                background: qlineargradient(x1:0, y1:0, x2:1, y2:1, stop:0 #5A5A5A, stop:1 #4A4A4A);
            }
        """)
        
        # Main layout
        main_widget = QWidget()
        self.setLayout(QHBoxLayout())
        self.layout().addWidget(main_widget)
        main_layout = QHBoxLayout()
        main_widget.setLayout(main_layout)
        
        # Control panel
        control_panel = QWidget(objectName="controlPanel")
        control_layout = QVBoxLayout()
        control_panel.setLayout(control_layout)
        control_panel.setMinimumWidth(300)
        control_panel.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)
        
        # Input form
        form_layout = QFormLayout()
        self.transformer_type = QComboBox()
        self.transformer_type.addItems(["Power", "Distribution", "Isolation", "Current", "Potential", "High-Frequency"])
        self.transformer_type.setFont(QFont("Consolas", 10))
        self.transformer_type.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        
        self.core_type = QComboBox()
        self.core_type.addItems(["Core-Type", "Shell-Type"])
        self.core_type.setFont(QFont("Consolas", 10))
        self.core_type.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        
        self.core_material = QComboBox()
        self.core_material.addItems(["Silicon Steel", "Ferrite"])
        self.core_material.setFont(QFont("Consolas", 10))
        self.core_material.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        
        self.winding_config = QComboBox()
        self.winding_config.addItems(["Concentric", "Interleaved"])
        self.winding_config.setFont(QFont("Consolas", 10))
        self.winding_config.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        
        self.voltage_in = QLineEdit("230")
        self.voltage_in.setFont(QFont("Consolas", 10))
        self.voltage_in.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.voltage_out = QLineEdit("12")
        self.voltage_out.setFont(QFont("Consolas", 10))
        self.voltage_out.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.frequency = QLineEdit("50")
        self.frequency.setFont(QFont("Consolas", 10))
        self.frequency.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.power = QLineEdit("1000")
        self.power.setFont(QFont("Consolas", 10))
        self.power.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.efficiency = QLineEdit("0.95")
        self.efficiency.setFont(QFont("Consolas", 10))
        self.efficiency.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.core_youngs_modulus = QLineEdit("200e9")
        self.core_youngs_modulus.setFont(QFont("Consolas", 10))
        self.core_youngs_modulus.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.winding_youngs_modulus = QLineEdit("110e9")
        self.winding_youngs_modulus.setFont(QFont("Consolas", 10))
        self.winding_youngs_modulus.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.magnetostriction_coeff = QLineEdit("4e-6")
        self.magnetostriction_coeff.setFont(QFont("Consolas", 10))
        self.magnetostriction_coeff.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.damping_ratio = QLineEdit("0.02")
        self.damping_ratio.setFont(QFont("Consolas", 10))
        self.damping_ratio.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.ambient_temp = QLineEdit("20")
        self.ambient_temp.setFont(QFont("Consolas", 10))
        self.ambient_temp.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.operating_temp = QLineEdit("80")
        self.operating_temp.setFont(QFont("Consolas", 10))
        self.operating_temp.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        
        form_layout.addRow(QLabel("Transformer Type:"), self.transformer_type)
        form_layout.addRow(QLabel("Core Configuration:"), self.core_type)
        form_layout.addRow(QLabel("Core Material:"), self.core_material)
        form_layout.addRow(QLabel("Winding Configuration:"), self.winding_config)
        form_layout.addRow(QLabel("Primary Voltage (V):"), self.voltage_in)
        form_layout.addRow(QLabel("Secondary Voltage (V):"), self.voltage_out)
        form_layout.addRow(QLabel("Frequency (Hz):"), self.frequency)
        form_layout.addRow(QLabel("Power (VA):"), self.power)
        form_layout.addRow(QLabel("Efficiency (0-1):"), self.efficiency)
        form_layout.addRow(QLabel("Core Young’s Modulus (Pa):"), self.core_youngs_modulus)
        form_layout.addRow(QLabel("Winding Young’s Modulus (Pa):"), self.winding_youngs_modulus)
        form_layout.addRow(QLabel("Magnetostriction Coeff:"), self.magnetostriction_coeff)
        form_layout.addRow(QLabel("Damping Ratio (0-1):"), self.damping_ratio)
        form_layout.addRow(QLabel("Ambient Temp (°C):"), self.ambient_temp)
        form_layout.addRow(QLabel("Operating Temp (°C):"), self.operating_temp)
        
        # Simulate button
        self.simulate_btn = QPushButton("Simulate")
        self.simulate_btn.setFont(QFont("Consolas", 10))
        self.simulate_btn.clicked.connect(self.run_simulation)
        self.simulate_btn.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        control_layout.addLayout(form_layout)
        control_layout.addWidget(self.simulate_btn)
        
        # Plot panel
        self.plot_widget = VibrationPlot()
        self.plot_widget.setStyleSheet("""
            background-color: #3C3C3C;
            border: 2px solid #1E1E1E;
            border-radius: 5px;
            padding: 10px;
        """)
        self.plot_widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        
        # Output panel
        output_panel = QWidget(objectName="outputPanel")
        output_layout = QVBoxLayout()
        output_panel.setLayout(output_layout)
        output_panel.setMinimumWidth(300)
        output_panel.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)
        
        output_label = QLabel("Vibration and Acoustic Noise Results")
        output_label.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.output_text = QTextEdit()
        self.output_text.setFont(QFont("Consolas", 10))
        self.output_text.setReadOnly(True)
        self.output_text.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        output_layout.addWidget(output_label)
        output_layout.addWidget(self.output_text)
        
        # Add panels to main layout
        main_layout.addWidget(control_panel, stretch=2)
        main_layout.addWidget(self.plot_widget, stretch=5)
        main_layout.addWidget(output_panel, stretch=2)
        
        self.transformer = Transformer()
        self.vibration_results = None

    def run_simulation(self):
        try:
            # Get inputs
            transformer_type = self.transformer_type.currentText()
            core_type = self.core_type.currentText()
            core_material = self.core_material.currentText()
            winding_config = self.winding_config.currentText()
            Vp = float(self.voltage_in.text())
            Vs = float(self.voltage_out.text())
            freq = float(self.frequency.text())
            power = float(self.power.text())
            eff = float(self.efficiency.text())
            core_E = float(self.core_youngs_modulus.text())
            winding_E = float(self.winding_youngs_modulus.text())
            lambda_s = float(self.magnetostriction_coeff.text())
            zeta = float(self.damping_ratio.text())
            ambient_temp = float(self.ambient_temp.text())
            operating_temp = float(self.operating_temp.text())
            
            # Validate inputs
            if Vp <= 0 or Vs <= 0 or freq <= 0 or power <= 0 or eff <= 0 or eff > 1:
                self.output_text.setText("Error: All inputs must be positive, and efficiency must be between 0 and 1.")
                return
            if core_E <= 0 or winding_E <= 0:
                self.output_text.setText("Error: Young’s modulus must be positive.")
                return
            if lambda_s <= 0:
                self.output_text.setText("Error: Magnetostriction coefficient must be positive.")
                return
            if zeta < 0 or zeta > 1:
                self.output_text.setText("Error: Damping ratio must be between 0 and 1.")
                return
            if ambient_temp < -50 or operating_temp < -50 or operating_temp > 200:
                self.output_text.setText("Error: Temperatures must be between -50°C and 200°C.")
                return
            
            # Get transformer parameters
            results = self.transformer.calculate(
                transformer_type, core_type, Vp, [Vs], freq, [power], eff,
                core_material=core_material, winding_config=winding_config,
                ambient_temp=ambient_temp, operating_temp=operating_temp
            )
            core_area = results['core_area']
            B = results.get('max_flux_density', 1.5)
            I_primary = power / (Vp * eff)
            B_leakage = results.get('leakage_flux_density', 0.01)
            
            print(f"Simulation Inputs: core_area={core_area:.4f} m², B={B:.2f} T, I_primary={I_primary:.2f} A, B_leakage={B_leakage:.4f} T")
            
            # Run vibration simulation
            self.vibration_results = self._run_vibration_simulation(
                core_area=core_area,
                B=B,
                I_primary=I_primary,
                B_leakage=B_leakage,
                core_E=core_E,
                winding_E=winding_E,
                lambda_s=lambda_s,
                zeta=zeta,
                freq=freq
            )
            
            # Update output
            max_displacement = np.max(np.sqrt(self.vibration_results['displacement_x']**2 + self.vibration_results['displacement_y']**2)) * 1e6
            max_velocity = np.max(np.sqrt(self.vibration_results['velocity_x']**2 + self.vibration_results['velocity_y']**2))
            spl = self.vibration_results['spl']
            
            output = f"Transformer Type: {transformer_type}\n"
            output += f"Core Type: {core_type}\n"
            output += f"Core Material: {core_material}\n"
            output += f"Winding Configuration: {winding_config}\n"
            output += f"Max Vibration Displacement: {max_displacement:.2f} μm\n"
            output += f"Max Vibration Velocity: {max_velocity:.4f} m/s\n"
            output += f"Sound Pressure Level (1 m): {spl:.1f} dB\n"
            output += f"Dominant Frequencies: {2*freq:.0f}, {4*freq:.0f} Hz\n"
            self.output_text.setText(output)
            
            # Plot vibration displacement
            displacement_magnitude = np.sqrt(self.vibration_results['displacement_x']**2 + self.vibration_results['displacement_y']**2)
            print(f"Displacement Magnitude: min={np.min(displacement_magnitude):.2e} m, max={np.max(displacement_magnitude):.2e} m")
            self.plot_widget.plot_vibration(
                self.vibration_results['x'],
                self.vibration_results['y'],
                displacement_magnitude
            )
            
        except ValueError:
            self.output_text.setText("Error: Invalid input. Please enter numeric values.")
        except Exception as e:
            self.output_text.setText(f"Error: {str(e)}\n{traceback.format_exc()}")

    def _run_vibration_simulation(self, core_area, B, I_primary, B_leakage, core_E, winding_E, lambda_s, zeta, freq):
        # 2D domain
        core_side = math.sqrt(core_area)
        winding_thickness = 0.01
        domain_width = core_side + 2 * winding_thickness
        domain_height = core_side
        
        # Mesh (20x10 grid)
        nx, ny = 20, 10
        x = np.linspace(0, domain_width, nx)
        y = np.linspace(0, domain_height, ny)
        X, Y = np.meshgrid(x, y)
        nodes = np.vstack((X.ravel(), Y.ravel())).T
        num_nodes = len(nodes)
        
        # Triangular elements
        triangles = []
        for i in range(ny - 1):
            for j in range(nx - 1):
                n1 = i * nx + j
                n2 = n1 + 1
                n3 = n1 + nx
                n4 = n3 + 1
                triangles.append([n1, n2, n3])
                triangles.append([n2, n4, n3])
        triangles = np.array(triangles)
        
        # Material properties
        E = np.ones(num_nodes) * core_E
        density = np.ones(num_nodes) * 7650
        poisson = np.ones(num_nodes) * 0.3
        core_region = (nodes[:, 0] >= winding_thickness) & (nodes[:, 0] <= core_side + winding_thickness) & (nodes[:, 1] >= 0) & (nodes[:, 1] <= core_side)
        winding_region = ~core_region
        E[winding_region] = winding_E
        density[winding_region] = 8960
        poisson[winding_region] = 0.34
        
        # Forces
        mu_0 = 4 * math.pi * 1e-7
        B_s = 2.0
        F = np.zeros((num_nodes, 2))
        epsilon_m = lambda_s * (B / B_s)**2
        stress_m = E[core_region] * epsilon_m
        F_em = (B**2 * core_area) / (2 * mu_0 * 0.1)
        F[core_region, 0] += stress_m / 0.1
        F[core_region, 1] += F_em / core_area
        F_lorentz = I_primary * B_leakage / 0.1
        F[winding_region, 0] += F_lorentz
        print(f"Forces: Core stress_m={np.mean(stress_m):.2e} Pa, F_em={F_em:.2e} N/m³, Winding F_lorentz={F_lorentz:.2e} N/m³")
        
        # FEM: [M]ü + [C]u̇ + [K]u = F
        M = sp.lil_matrix((2*num_nodes, 2*num_nodes))
        C = sp.lil_matrix((2*num_nodes, 2*num_nodes))
        K = sp.lil_matrix((2*num_nodes, 2*num_nodes))
        F_vec = np.zeros(2*num_nodes)
        
        for tri in triangles:
            n1, n2, n3 = tri
            xe = nodes[[n1, n2, n3], 0]
            ye = nodes[[n1, n2, n3], 1]
            area = 0.5 * abs(xe[0]*(ye[1] - ye[2]) + xe[1]*(ye[2] - ye[0]) + xe[2]*(ye[0] - ye[1]))
            b = np.array([ye[1] - ye[2], ye[2] - ye[0], ye[0] - ye[1]]) / (2 * area)
            c = np.array([xe[2] - xe[1], xe[0] - xe[2], xe[1] - xe[0]]) / (2 * area)
            E_elem = np.mean(E[[n1, n2, n3]])
            nu_elem = np.mean(poisson[[n1, n2, n3]])
            rho_elem = np.mean(density[[n1, n2, n3]])
            D = (E_elem / (1 - nu_elem**2)) * np.array([
                [1, nu_elem, 0],
                [nu_elem, 1, 0],
                [0, 0, (1 - nu_elem)/2]
            ])
            B_elem = np.zeros((3, 6))
            for i in range(3):
                B_elem[0, 2*i] = b[i]
                B_elem[1, 2*i + 1] = c[i]
                B_elem[2, 2*i] = c[i]
                B_elem[2, 2*i + 1] = b[i]
            Ke = area * 0.1 * (B_elem.T @ D @ B_elem)
            Me = (rho_elem * area * 0.1 / 3) * np.eye(6)
            for i in range(3):
                for j in range(3):
                    idx_i = [2*(n1+i), 2*(n1+i)+1]
                    idx_j = [2*(n1+j), 2*(n1+j)+1]
                    K[np.ix_(idx_i, idx_j)] += Ke[2*i:2*i+2, 2*j:2*j+2]
                    M[np.ix_(idx_i, idx_j)] += Me[2*i:2*i+2, 2*j:2*j+2]
            F_elem = np.mean(F[[n1, n2, n3]], axis=0) * (area * 0.1 / 3)
            for i in range(3):
                F_vec[2*(n1+i)] += F_elem[0]
                F_vec[2*(n1+i)+1] += F_elem[1]
        
        # Boundary conditions: Fix nodes at y=0
        fixed_nodes = np.where(nodes[:, 1] == 0)[0]
        print(f"Applying boundary conditions to {len(fixed_nodes)} nodes at y=0")
        for node in fixed_nodes:
            K[2*node, 2*node] += 1e12
            K[2*node+1, 2*node+1] += 1e12
            F_vec[2*node] = 0
            F_vec[2*node+1] = 0
        
        # Damping: Rayleigh damping
        omega = 2 * math.pi * (2 * freq)
        C = 2 * zeta * omega * M
        
        # Harmonic solution
        A = -omega**2 * M + 1j * omega * C + K
        A = A.tocsr()
        cond = np.linalg.cond(A.toarray()) if num_nodes < 50 else "N/A (large matrix)"
        print(f"Matrix condition number: {cond}")
        u, istop, itn = spla.lsqr(A, F_vec)[:3]
        print(f"LSQR solver: istop={istop}, iterations={itn}")
        u_x = u[::2].real
        u_y = u[1::2].real
        v_x = (omega * u[::2].imag)
        v_y = (omega * u[1::2].imag)
        if not np.all(np.isfinite(u)):
            print("Warning: Non-finite displacements detected")
            u_x = np.zeros_like(u_x)
            u_y = np.zeros_like(u_y)
            v_x = np.zeros_like(v_x)
            v_y = np.zeros_like(v_y)
        print(f"Displacements: u_x min={np.min(u_x):.2e} m, max={np.max(u_x):.2e} m, u_y min={np.min(u_y):.2e} m, max={np.max(u_y):.2e} m")
        
        # Acoustic noise
        rho_air = 1.225
        c_air = 343
        S = 2 * (domain_width + domain_height) * 0.1
        v_rms = np.sqrt(np.mean(v_x**2 + v_y**2))
        P_acoustic = rho_air * c_air * S * v_rms**2
        r = 1.0
        p_rms = np.sqrt((P_acoustic * rho_air * c_air) / (4 * math.pi * r**2))
        p_0 = 20e-6
        spl = 20 * np.log10(p_rms / p_0) if p_rms > 0 else 0
        
        return {
            'x': X,
            'y': Y,
            'displacement_x': u_x.reshape(X.shape),
            'displacement_y': u_y.reshape(X.shape),
            'velocity_x': v_x.reshape(X.shape),
            'velocity_y': v_y.reshape(X.shape),
            'spl': spl
        }

class VibrationPlot(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.canvas.updateGeometry()
        self.layout = QVBoxLayout()
        self.layout.addWidget(self.canvas)
        self.setLayout(self.layout)
        
        self.ax = self.figure.add_subplot(111)
        self.ax.set_title("Vibration Displacement Magnitude (μm)")
        self.ax.set_xlabel("x (m)")
        self.ax.set_ylabel("y (m)")
        
    def plot_vibration(self, x, y, displacement):
        try:
            print(f"Plotting: displacement min={np.min(displacement):.2e} m, max={np.max(displacement):.2e} m")
            if np.all(displacement == 0) or np.any(np.isnan(displacement)):
                print("Warning: Displacement contains all zeros or NaNs")
                displacement = np.ones_like(displacement) * 1e-9  # Fallback for visualization
            
            self.ax.clear()
            displacement_um = displacement * 1e6
            levels = np.linspace(max(np.min(displacement_um), 0), max(np.max(displacement_um), 1e-6), 50)
            contour = self.ax.contourf(x, y, displacement_um, levels=levels, cmap='viridis')
            cbar = self.figure.colorbar(contour, ax=self.ax, label='Displacement (μm)')
            
            self.ax.set_title("Vibration Displacement Magnitude (μm)")
            self.ax.set_xlabel("x (m)")
            self.ax.set_ylabel("y (m)")
            
            self.figure.tight_layout()
            self.canvas.draw()
            self.canvas.draw_idle()
            self.canvas.flush_events()
            print("Plot updated successfully")
            
        except Exception as e:
            print(f"Plot error: {str(e)}\n{traceback.format_exc()}")