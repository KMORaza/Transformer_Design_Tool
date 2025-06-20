from PyQt5.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QLabel, QLineEdit, QPushButton, QTextEdit, QComboBox, QSizePolicy, QMessageBox
from PyQt5.QtGui import QFont
from PyQt5.QtCore import Qt
from transformer import Transformer
from plot import TransformerPlot
from MultiWinding import MultiWindingWindow
from MagneticFluxDensity import MagneticFluxDensityWindow
from CoreSaturationAnalyzer import CoreSaturationAnalyzerWindow
from LeakageAndStrayFieldEst import LeakageAndStrayFieldEstWindow
from WindingConfigOptimizer import WindingConfigOptimizerWindow
from CoreLossCopperLossEst import CoreLossCopperLossEstWindow
from TempCorrectedResistanceModeling import TempCorrectedResistanceModelingWindow
from ThermalSimulationFEM import ThermalSimulationFEMWindow
from VibrationAndAcousticNoiseModel import VibrationAndAcousticNoiseModelWindow
from HFRA import HFRAWindow
import traceback

class TransformerDesignApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Transformer Design Tool")
        self.setGeometry(100, 100, 1200, 600)
        self.setFont(QFont("Consolas", 10))
        
        # Apply dark theme stylesheet
        self.setStyleSheet("""
            QMainWindow {
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
        """)
        
        # Main widget and layout
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        main_layout = QHBoxLayout()
        main_widget.setLayout(main_layout)
        
        # Control panel
        control_panel = QWidget(objectName="controlPanel")
        control_layout = QFormLayout()
        control_panel.setLayout(control_layout)
        control_panel.setMinimumWidth(300)
        control_panel.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)
        
        # Transformer type selection
        self.transformer_type = QComboBox()
        self.transformer_type.addItems(["Power", "Distribution", "Isolation", "Current", "Potential", "High-Frequency"])
        self.transformer_type.setFont(QFont("Consolas", 10))
        self.transformer_type.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        
        # Core configuration selection
        self.core_type = QComboBox()
        self.core_type.addItems(["Core-Type", "Shell-Type"])
        self.core_type.setFont(QFont("Consolas", 10))
        self.core_type.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        
        # Input fields
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
        
        control_layout.addRow(QLabel("Transformer Type:"), self.transformer_type)
        control_layout.addRow(QLabel("Core Configuration:"), self.core_type)
        control_layout.addRow(QLabel("Primary Voltage (V):"), self.voltage_in)
        control_layout.addRow(QLabel("Secondary Voltage (V):"), self.voltage_out)
        control_layout.addRow(QLabel("Frequency (Hz):"), self.frequency)
        control_layout.addRow(QLabel("Power (VA):"), self.power)
        control_layout.addRow(QLabel("Efficiency (0-1):"), self.efficiency)
        
        # Calculate button
        self.calculate_btn = QPushButton("Calculate")
        self.calculate_btn.setFont(QFont("Consolas", 10))
        self.calculate_btn.clicked.connect(self.calculate_design)
        self.calculate_btn.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        control_layout.addRow(self.calculate_btn)
        
        # Multi-Winding button
        self.multi_winding_btn = QPushButton("Multi-Winding")
        self.multi_winding_btn.setFont(QFont("Consolas", 10))
        self.multi_winding_btn.clicked.connect(self.open_multi_winding_window)
        self.multi_winding_btn.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        control_layout.addRow(self.multi_winding_btn)
        
        # Magnetic Flux Density button
        self.flux_density_btn = QPushButton("Magnetic Flux Density")
        self.flux_density_btn.setFont(QFont("Consolas", 10))
        self.flux_density_btn.clicked.connect(self.open_flux_density_window)
        self.flux_density_btn.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        control_layout.addRow(self.flux_density_btn)
        
        # Core Saturation Analysis button
        self.core_saturation_btn = QPushButton("Core Saturation Analysis")
        self.core_saturation_btn.setFont(QFont("Consolas", 10))
        self.core_saturation_btn.clicked.connect(self.open_core_saturation_window)
        self.core_saturation_btn.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        control_layout.addRow(self.core_saturation_btn)
        
        # Leakage and Stray Field button
        self.leakage_stray_btn = QPushButton("Leakage and Stray Field")
        self.leakage_stray_btn.setFont(QFont("Consolas", 10))
        self.leakage_stray_btn.clicked.connect(self.open_leakage_stray_window)
        self.leakage_stray_btn.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        control_layout.addRow(self.leakage_stray_btn)
        
        # Winding Optimization button
        self.winding_optimization_btn = QPushButton("Winding Optimization")
        self.winding_optimization_btn.setFont(QFont("Consolas", 10))
        self.winding_optimization_btn.clicked.connect(self.open_winding_optimization_window)
        self.winding_optimization_btn.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        control_layout.addRow(self.winding_optimization_btn)
        
        # Core and Copper Loss button
        self.loss_estimation_btn = QPushButton("Core and Copper Loss")
        self.loss_estimation_btn.setFont(QFont("Consolas", 10))
        self.loss_estimation_btn.clicked.connect(self.open_loss_estimation_window)
        self.loss_estimation_btn.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        control_layout.addRow(self.loss_estimation_btn)
        
        # Temperature-Corrected Resistance button
        self.temp_resistance_btn = QPushButton("Temperature-Corrected Resistance")
        self.temp_resistance_btn.setFont(QFont("Consolas", 10))
        self.temp_resistance_btn.clicked.connect(self.open_temp_resistance_window)
        self.temp_resistance_btn.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        control_layout.addRow(self.temp_resistance_btn)
        
        # Thermal Simulation (FEM) button
        self.thermal_simulation_btn = QPushButton("Thermal Simulation (FEM)")
        self.thermal_simulation_btn.setFont(QFont("Consolas", 10))
        self.thermal_simulation_btn.clicked.connect(self.open_thermal_simulation_window)
        self.thermal_simulation_btn.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        control_layout.addRow(self.thermal_simulation_btn)
        
        # Vibration and Acoustic Noise Model button
        self.vibration_noise_btn = QPushButton("Vibration and Acoustic Noise")
        self.vibration_noise_btn.setFont(QFont("Consolas", 10))
        self.vibration_noise_btn.clicked.connect(self.open_vibration_noise_window)
        self.vibration_noise_btn.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        control_layout.addRow(self.vibration_noise_btn)
        
        # Harmonic and Frequency Response button
        self.hfra_btn = QPushButton("Harmonic and Frequency Response")
        self.hfra_btn.setFont(QFont("Consolas", 10))
        self.hfra_btn.clicked.connect(self.open_hfra_window)
        self.hfra_btn.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        control_layout.addRow(self.hfra_btn)
        
        # Output panel
        output_panel = QWidget(objectName="outputPanel")
        output_layout = QVBoxLayout()
        output_panel.setLayout(output_layout)
        output_panel.setMinimumWidth(300)
        output_panel.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)
        
        output_label = QLabel("Design Results")
        output_label.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.output_text = QTextEdit()
        self.output_text.setFont(QFont("Consolas", 10))
        self.output_text.setReadOnly(True)
        self.output_text.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        output_layout.addWidget(output_label)
        output_layout.addWidget(self.output_text)
        
        # Plot panel
        self.plot_widget = TransformerPlot()
        self.plot_widget.setStyleSheet("""
            background-color: #3C3C3C;
            border: 2px solid #1E1E1E;
            border-radius: 5px;
            padding: 10px;
        """)
        self.plot_widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        
        # Add panels to main layout with stretch factors
        main_layout.addWidget(control_panel, stretch=2)
        main_layout.addWidget(self.plot_widget, stretch=5)
        main_layout.addWidget(output_panel, stretch=2)
        
        self.transformer = Transformer()

    def calculate_design(self):
        try:
            transformer_type = self.transformer_type.currentText()
            core_type = self.core_type.currentText()
            Vp = float(self.voltage_in.text())
            Vs = float(self.voltage_out.text())
            freq = float(self.frequency.text())
            power = float(self.power.text())
            eff = float(self.efficiency.text())
            
            # Validate inputs
            if Vp <= 0 or Vs <= 0 or freq <= 0 or power <= 0 or eff <= 0 or eff > 1:
                self.output_text.setText("Error: All inputs must be positive, and efficiency must be between 0 and 1.")
                return
            
            # Perform calculations
            results = self.transformer.calculate(transformer_type, core_type, Vp, [Vs], freq, [power], eff)
            
            # Update output
            output = f"Transformer Type: {transformer_type}\n"
            output += f"Core Type: {core_type}\n"
            output += f"Core Area: {results['core_area']:.4f} m²\n"
            output += f"Turns Ratio: {results['turns_ratios'][0]:.2f}\n"
            output += f"Primary Turns: {results['primary_turns']:.0f}\n"
            output += f"Secondary Turns: {results['secondary_turns'][0]:.0f}\n"
            output += f"Primary Wire Gauge: {results['primary_wire_gauge']} AWG\n"
            output += f"Secondary Wire Gauge: {results['secondary_wire_gauges'][0]} AWG\n"
            output += f"Interleaving Pattern: {results['interleaving_pattern']}\n"
            output += f"Primary Layers: {results['primary_layers']}\n"
            output += f"Secondary Layers: {results['secondary_layers'][0]}\n"
            output += f"Core Loss: {results['core_loss']:.2f} W\n"
            output += f"Efficiency: {results['efficiency']*100:.2f}%\n"
            output += f"Peak Flux Density: {results['max_flux_density']:.2f} T\n"
            self.output_text.setText(output)
            
            # Update plot with secondary flux density and time points
            self.plot_widget.update_plot(results['flux_densities'][1], results['time_points'], 1.0 / freq)
            
        except ValueError:
            self.output_text.setText("Error: Invalid input. Please enter numeric values.")
        except Exception as e:
            self.output_text.setText(f"Error: {str(e)}\n{traceback.format_exc()}")

    def open_multi_winding_window(self):
        try:
            self.multi_winding_window = MultiWindingWindow(self)
            self.multi_winding_window.show()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to open Multi-Winding window: {str(e)}\n{traceback.format_exc()}")

    def open_flux_density_window(self):
        try:
            self.flux_density_window = MagneticFluxDensityWindow(self)
            self.flux_density_window.show()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to open Magnetic Flux Density window: {str(e)}\n{traceback.format_exc()}")

    def open_core_saturation_window(self):
        try:
            self.core_saturation_window = CoreSaturationAnalyzerWindow(self)
            self.core_saturation_window.show()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to open Core Saturation Analysis window: {str(e)}\n{traceback.format_exc()}")

    def open_leakage_stray_window(self):
        try:
            self.leakage_stray_window = LeakageAndStrayFieldEstWindow(self)
            self.leakage_stray_window.show()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to open Leakage and Stray Field window: {str(e)}\n{traceback.format_exc()}")

    def open_winding_optimization_window(self):
        try:
            self.winding_optimization_window = WindingConfigOptimizerWindow(self)
            self.winding_optimization_window.show()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to open Winding Optimization window: {str(e)}\n{traceback.format_exc()}")

    def open_loss_estimation_window(self):
        try:
            self.loss_estimation_window = CoreLossCopperLossEstWindow(self)
            self.loss_estimation_window.show()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to open Core and Copper Loss window: {str(e)}\n{traceback.format_exc()}")

    def open_temp_resistance_window(self):
        try:
            self.temp_resistance_window = TempCorrectedResistanceModelingWindow(self)
            self.temp_resistance_window.show()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to open Temperature-Corrected Resistance window: {str(e)}\n{traceback.format_exc()}")

    def open_thermal_simulation_window(self):
        try:
            self.thermal_simulation_window = ThermalSimulationFEMWindow(self)
            self.thermal_simulation_window.show()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to open Thermal Simulation (FEM) window: {str(e)}\n{traceback.format_exc()}")

    def open_vibration_noise_window(self):
        try:
            self.vibration_noise_window = VibrationAndAcousticNoiseModelWindow(self)
            self.vibration_noise_window.show()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to open Vibration and Acoustic Noise window: {str(e)}\n{traceback.format_exc()}")

    def open_hfra_window(self):
        try:
            self.hfra_window = HFRAWindow(self)
            self.hfra_window.show()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to open Harmonic and Frequency Response window: {str(e)}\n{traceback.format_exc()}")