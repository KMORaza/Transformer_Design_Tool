from PyQt5.QtWidgets import QDialog, QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QLabel, QLineEdit, QPushButton, QTextEdit, QComboBox, QSizePolicy, QCheckBox
from PyQt5.QtGui import QFont
from PyQt5.QtCore import Qt
from transformer import Transformer
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import math
import traceback

class ThermalSimulationFEMWindow(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Thermal Simulation (FEM)")
        self.setGeometry(200, 200, 1000, 600)
        self.setFont(QFont("Consolas", 10))
        
        # Apply dark theme stylesheet
        self.setStyleSheet("""
            QDialog {
                background-color: #2E2E2E;
            }
            QWidget#controlPanel, QWidget#outputPanel {
                background: qlineargradient(x1:0, y1:0, x2:1, y2:1, stop:0 #4A4A4A, stop:1 #3A3A3A);
                border: 2px solid #1E1E1E;
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
            QCheckBox {
                color: #D3D3D3;
                font: 10pt "Consolas";
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
        self.ambient_temp = QLineEdit("20")
        self.ambient_temp.setFont(QFont("Consolas", 10))
        self.ambient_temp.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.operating_temp = QLineEdit("80")
        self.operating_temp.setFont(QFont("Consolas", 10))
        self.operating_temp.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.core_thermal_k = QLineEdit("30")
        self.core_thermal_k.setFont(QFont("Consolas", 10))
        self.core_thermal_k.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.winding_thermal_k = QLineEdit("400")
        self.winding_thermal_k.setFont(QFont("Consolas", 10))
        self.winding_thermal_k.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.emissivity = QLineEdit("0.8")
        self.emissivity.setFont(QFont("Consolas", 10))
        self.emissivity.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.convection_coeff = QLineEdit("10")
        self.convection_coeff.setFont(QFont("Consolas", 10))
        self.convection_coeff.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.airflow_velocity = QLineEdit("0.5")
        self.airflow_velocity.setFont(QFont("Consolas", 10))
        self.airflow_velocity.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        
        form_layout.addRow(QLabel("Transformer Type:"), self.transformer_type)
        form_layout.addRow(QLabel("Core Configuration:"), self.core_type)
        form_layout.addRow(QLabel("Core Material:"), self.core_material)
        form_layout.addRow(QLabel("Winding Configuration:"), self.winding_config)
        form_layout.addRow(QLabel("Primary Voltage (V):"), self.voltage_in)
        form_layout.addRow(QLabel("Secondary Voltage (V):"), self.voltage_out)
        form_layout.addRow(QLabel("Frequency (Hz):"), self.frequency)
        form_layout.addRow(QLabel("Power (VA):"), self.power)
        form_layout.addRow(QLabel("Efficiency (0-1):"), self.efficiency)
        form_layout.addRow(QLabel("Ambient Temp (°C):"), self.ambient_temp)
        form_layout.addRow(QLabel("Operating Temp (°C):"), self.operating_temp)
        form_layout.addRow(QLabel("Core Thermal k (W/m·K):"), self.core_thermal_k)
        form_layout.addRow(QLabel("Winding Thermal k (W/m·K):"), self.winding_thermal_k)
        form_layout.addRow(QLabel("Emissivity (0-1):"), self.emissivity)
        form_layout.addRow(QLabel("Convection Coeff (W/m²·K):"), self.convection_coeff)
        form_layout.addRow(QLabel("Airflow Velocity (m/s):"), self.airflow_velocity)
        
        # Show heat flux checkbox
        self.show_heat_flux = QCheckBox("Show Heat Flux Vectors")
        self.show_heat_flux.setFont(QFont("Consolas", 10))
        self.show_heat_flux.setChecked(False)
        
        # Simulate button
        self.simulate_btn = QPushButton("Simulate")
        self.simulate_btn.setFont(QFont("Consolas", 10))
        self.simulate_btn.clicked.connect(self.run_simulation)
        self.simulate_btn.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        
        control_layout.addLayout(form_layout)
        control_layout.addWidget(self.show_heat_flux)
        control_layout.addWidget(self.simulate_btn)
        
        # Plot panel
        self.plot_widget = ThermalPlot()
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
        
        output_label = QLabel("Thermal Simulation Results")
        output_label.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.output_text = QTextEdit()
        self.output_text.setFont(QFont("Consolas", 10))
        self.output_text.setReadOnly(True)
        self.output_text.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        output_layout.addWidget(output_label)
        output_layout.addWidget(self.output_text)
        
        # Add panels to main layout with stretch factors
        main_layout.addWidget(control_panel, stretch=2)
        main_layout.addWidget(self.plot_widget, stretch=5)
        main_layout.addWidget(output_panel, stretch=2)
        
        self.transformer = Transformer()

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
            ambient_temp = float(self.ambient_temp.text())
            operating_temp = float(self.operating_temp.text())
            core_thermal_k = float(self.core_thermal_k.text())
            winding_thermal_k = float(self.winding_thermal_k.text())
            emissivity = float(self.emissivity.text())
            convection_coeff = float(self.convection_coeff.text())
            airflow_velocity = float(self.airflow_velocity.text())
            show_heat_flux = self.show_heat_flux.isChecked()
            
            # Validate inputs
            if Vp <= 0 or Vs <= 0 or freq <= 0 or power <= 0 or eff <= 0 or eff > 1:
                self.output_text.setText("Error: All inputs must be positive, and efficiency must be between 0 and 1.")
                return
            if ambient_temp < -50 or operating_temp < -50 or operating_temp > 200:
                self.output_text.setText("Error: Temperatures must be between -50°C and 200°C.")
                return
            if core_thermal_k <= 0 or winding_thermal_k <= 0:
                self.output_text.setText("Error: Thermal conductivities must be positive.")
                return
            if emissivity < 0 or emissivity > 1:
                self.output_text.setText("Error: Emissivity must be between 0 and 1.")
                return
            if convection_coeff <= 0 or airflow_velocity <= 0:
                self.output_text.setText("Error: Convection coefficient and airflow velocity must be positive.")
                return
            
            # Get transformer losses
            results = self.transformer.calculate(
                transformer_type, core_type, Vp, [Vs], freq, [power], eff,
                core_material=core_material, winding_config=winding_config,
                ambient_temp=ambient_temp, operating_temp=operating_temp
            )
            core_loss = results['core_loss']
            copper_loss = results['total_copper_loss']
            core_area = results['core_area']
            
            # Run FEM thermal simulation
            thermal_results = self._run_fem_simulation(
                core_loss=core_loss,
                copper_loss=copper_loss,
                core_area=core_area,
                core_thermal_k=core_thermal_k,
                winding_thermal_k=winding_thermal_k,
                emissivity=emissivity,
                convection_coeff=convection_coeff,
                ambient_temp=ambient_temp + 273.15,  # Convert to Kelvin
                airflow_velocity=airflow_velocity,
                operating_temp=operating_temp
            )
            
            # Update output
            max_temp_c = thermal_results['max_temperature'] - 273.15
            avg_temp_c = thermal_results['avg_temperature'] - 273.15
            conduction_loss = thermal_results['conduction_loss']
            convection_loss = thermal_results['convection_loss']
            radiation_loss = thermal_results['radiation_loss']
            hotspot_x, hotspot_y = thermal_results['hotspot_position']
            
            output = f"Transformer Type: {transformer_type}\n"
            output += f"Core Type: {core_type}\n"
            output += f"Core Material: {core_material}\n"
            output += f"Winding Configuration: {winding_config}\n"
            output += f"Max Temperature (Hotspot): {max_temp_c:.2f} °C\n"
            output += f"Average Temperature: {avg_temp_c:.2f} °C\n"
            output += f"Hotspot Location: ({hotspot_x:.3f}, {hotspot_y:.3f}) m\n"
            output += f"Conduction Heat Loss: {conduction_loss:.2f} W\n"
            output += f"Convection Heat Loss: {convection_loss:.2f} W\n"
            output += f"Radiation Heat Loss: {radiation_loss:.4f} W\n"
            self.output_text.setText(output)
            
            # Update plot
            self.plot_widget.update_plot(
                thermal_results['x'],
                thermal_results['y'],
                thermal_results['temperature'],
                thermal_results['hotspot_position'],
                thermal_results['heat_flux'],
                show_heat_flux
            )
            
        except ValueError:
            self.output_text.setText("Error: Invalid input. Please enter numeric values.")
        except Exception as e:
            self.output_text.setText(f"Error: {str(e)}\n{traceback.format_exc()}")

    def _run_fem_simulation(self, core_loss, copper_loss, core_area, core_thermal_k, winding_thermal_k, emissivity, convection_coeff, ambient_temp, airflow_velocity, operating_temp):
        # Simplified 2D domain: core (square) with windings (annulus)
        core_side = math.sqrt(core_area)
        winding_thickness = 0.01  # Assume 10 mm winding thickness
        domain_width = core_side + 2 * winding_thickness
        domain_height = core_side
        
        # Generate mesh (20x10 grid for simplicity)
        nx, ny = 20, 10
        x = np.linspace(0, domain_width, nx)
        y = np.linspace(0, domain_height, ny)
        X, Y = np.meshgrid(x, y)
        nodes = np.vstack((X.ravel(), Y.ravel())).T
        num_nodes = len(nodes)
        
        # Generate triangular elements
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
        
        # Assign material properties and heat sources
        thermal_k = np.ones(num_nodes) * core_thermal_k
        heat_source = np.zeros(num_nodes)
        core_region = (nodes[:, 0] >= winding_thickness) & (nodes[:, 0] <= core_side + winding_thickness) & (nodes[:, 1] >= 0) & (nodes[:, 1] <= core_side)
        winding_region = ~core_region
        thermal_k[winding_region] = winding_thermal_k
        
        # Distribute heat sources (W/m³)
        core_volume = core_area * 0.1  # Assume 0.1 m depth
        winding_volume = (domain_width * domain_height - core_area) * 0.1
        heat_source[core_region] = core_loss / core_volume
        heat_source[winding_region] = copper_loss / winding_volume
        
        # Assemble FEM system: [K]*{T} = {F}
        K = sp.lil_matrix((num_nodes, num_nodes))
        F = np.zeros(num_nodes)
        
        for tri in triangles:
            # Element nodes
            n1, n2, n3 = tri
            xe = nodes[[n1, n2, n3], 0]
            ye = nodes[[n1, n2, n3], 1]
            
            # Element area
            area = 0.5 * abs(xe[0]*(ye[1] - ye[2]) + xe[1]*(ye[2] - ye[0]) + xe[2]*(ye[0] - ye[1]))
            
            # Shape function derivatives
            b = np.array([ye[1] - ye[2], ye[2] - ye[0], ye[0] - ye[1]]) / (2 * area)
            c = np.array([xe[2] - xe[1], xe[0] - xe[2], xe[1] - xe[0]]) / (2 * area)
            
            # Element conductivity (average of nodes)
            k_elem = np.mean(thermal_k[[n1, n2, n3]])
            
            # Element stiffness matrix
            Ke = k_elem * area * (np.outer(b, b) + np.outer(c, c))
            for i in range(3):
                for j in range(3):
                    K[n1+i, n1+j] += Ke[i, j]
            
            # Element heat source
            q_elem = np.mean(heat_source[[n1, n2, n3]])
            Fe = q_elem * area / 3 * np.ones(3)
            F[[n1, n2, n3]] += Fe
        
        # Boundary conditions (convection and radiation)
        sigma = 5.67e-8  # Stefan-Boltzmann constant
        boundary_nodes = np.where(
            (nodes[:, 0] < 1e-6) | (nodes[:, 0] > domain_width - 1e-6) |
            (nodes[:, 1] < 1e-6) | (nodes[:, 1] > domain_height - 1e-6)
        )[0]
        
        # Approximate surface length per boundary node
        perimeter = 2 * (domain_width + domain_height)
        num_boundary_nodes = len(boundary_nodes)
        length_per_node = perimeter / num_boundary_nodes if num_boundary_nodes > 0 else 1e-6
        
        # Linearize radiation: h_rad ≈ 4*ε*σ*T_ref^3
        T_ref = (ambient_temp + operating_temp + 273.15) / 2
        h_rad = 4 * emissivity * sigma * T_ref**3
        h_total = convection_coeff + h_rad
        
        for n in boundary_nodes:
            K[n, n] += h_total * length_per_node
            F[n] += h_total * ambient_temp * length_per_node
        
        # Solve system
        K = K.tocsr()
        T = spla.spsolve(K, F)
        
        # Calculate heat fluxes
        grad_x = np.zeros(num_nodes)
        grad_y = np.zeros(num_nodes)
        for tri in triangles:
            n1, n2, n3 = tri
            xe = nodes[[n1, n2, n3], 0]
            ye = nodes[[n1, n2, n3], 1]
            Te = T[[n1, n2, n3]]
            area = 0.5 * abs(xe[0]*(ye[1] - ye[2]) + xe[1]*(ye[2] - ye[0]) + xe[2]*(ye[0] - ye[1]))
            b = np.array([ye[1] - ye[2], ye[2] - ye[0], ye[0] - ye[1]]) / (2 * area)
            c = np.array([xe[2] - xe[1], xe[0] - xe[2], xe[1] - xe[0]]) / (2 * area)
            grad_x_elem = np.dot(b, Te)
            grad_y_elem = np.dot(c, Te)
            grad_x[[n1, n2, n3]] += grad_x_elem / 3
            grad_y[[n1, n2, n3]] += grad_y_elem / 3
        heat_flux = -thermal_k[:, np.newaxis] * np.vstack((grad_x, grad_y)).T
        
        # Calculate heat losses
        conduction_loss = np.sum(heat_source) * (domain_width * domain_height * 0.1)
        convection_loss = np.sum([convection_coeff * length_per_node * (T[n] - ambient_temp) for n in boundary_nodes])
        radiation_loss = np.sum([emissivity * sigma * length_per_node * (T[n]**4 - ambient_temp**4) for n in boundary_nodes])
        
        # Find hotspot
        max_temp_idx = np.argmax(T)
        hotspot_position = nodes[max_temp_idx]
        
        return {
            'x': X,
            'y': Y,
            'temperature': T.reshape(X.shape),
            'max_temperature': np.max(T),
            'avg_temperature': np.mean(T),
            'hotspot_position': hotspot_position,
            'heat_flux': heat_flux.reshape(X.shape + (2,)),
            'conduction_loss': conduction_loss,
            'convection_loss': convection_loss,
            'radiation_loss': radiation_loss
        }

class ThermalPlot(QWidget):
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
        self.ax.set_title("Temperature Distribution (°C)")
        self.ax.set_xlabel("x (m)")
        self.ax.set_ylabel("y (m)")
        
    def update_plot(self, x, y, temperature, hotspot_position, heat_flux, show_heat_flux):
        try:
            self.ax.clear()
            T_celsius = temperature - 273.15
            contour = self.ax.contourf(x, y, T_celsius, levels=50, cmap='hot')
            self.figure.colorbar(contour, ax=self.ax, label='Temperature (°C)')
            
            # Mark hotspot
            self.ax.plot(hotspot_position[0], hotspot_position[1], 'r*', markersize=15, label=f'Hotspot: {np.max(T_celsius):.1f}°C')
            
            # Plot heat flux vectors
            if show_heat_flux:
                U, V = heat_flux[:, :, 0], heat_flux[:, :, 1]
                self.ax.quiver(x, y, U, V, color='blue', scale=1e3)
            
            self.ax.set_title("Temperature Distribution (°C)")
            self.ax.set_xlabel("x (m)")
            self.ax.set_ylabel("y (m)")
            self.ax.legend()
            self.figure.tight_layout()
            self.canvas.draw()
            self.canvas.flush_events()
        except Exception as e:
            print(f"Plot error: {str(e)}\n{traceback.format_exc()}")