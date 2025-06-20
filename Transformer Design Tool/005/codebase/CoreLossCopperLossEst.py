from PyQt5.QtWidgets import QDialog, QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QLabel, QLineEdit, QPushButton, QTextEdit, QComboBox, QSizePolicy
from PyQt5.QtGui import QFont
from PyQt5.QtCore import Qt
from transformer import Transformer
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import math
import traceback

class CoreLossCopperLossEstWindow(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Core & Copper Loss Estimation")
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
        
        form_layout.addRow(QLabel("Transformer Type:"), self.transformer_type)
        form_layout.addRow(QLabel("Core Configuration:"), self.core_type)
        form_layout.addRow(QLabel("Core Material:"), self.core_material)
        form_layout.addRow(QLabel("Winding Configuration:"), self.winding_config)
        form_layout.addRow(QLabel("Primary Voltage (V):"), self.voltage_in)
        form_layout.addRow(QLabel("Secondary Voltage (V):"), self.voltage_out)
        form_layout.addRow(QLabel("Frequency (Hz):"), self.frequency)
        form_layout.addRow(QLabel("Power (VA):"), self.power)
        form_layout.addRow(QLabel("Efficiency (0-1):"), self.efficiency)
        
        # Calculate button
        self.calculate_btn = QPushButton("Calculate")
        self.calculate_btn.setFont(QFont("Consolas", 10))
        self.calculate_btn.clicked.connect(self.calculate_losses)
        self.calculate_btn.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        control_layout.addLayout(form_layout)
        control_layout.addWidget(self.calculate_btn)
        
        # Plot panel
        self.plot_widget = LossPlot()
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
        
        output_label = QLabel("Loss Estimation Results")
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

    def calculate_losses(self):
        try:
            transformer_type = self.transformer_type.currentText()
            core_type = self.core_type.currentText()
            core_material = self.core_material.currentText()
            winding_config = self.winding_config.currentText()
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
            results = self.transformer.calculate(
                transformer_type, core_type, Vp, [Vs], freq, [power], eff,
                core_material=core_material, winding_config=winding_config
            )
            
            # Update output
            output = f"Transformer Type: {transformer_type}\n"
            output += f"Core Type: {core_type}\n"
            output += f"Core Material: {core_material}\n"
            output += f"Winding Configuration: {winding_config}\n"
            output += f"Total Core Loss: {results['core_loss']:.2f} W\n"
            output += f"Hysteresis Loss: {results['hysteresis_loss']:.2f} W\n"
            output += f"Eddy Current Loss: {results['eddy_current_loss']:.2f} W\n"
            output += f"Total Copper Loss: {results['total_copper_loss']:.2f} W\n"
            output += f"Primary Copper Loss: {results['copper_loss_primary']:.2f} W\n"
            output += f"Secondary Copper Loss: {results['copper_loss_secondary']:.2f} W\n"
            output += f"Overall Efficiency: {results['efficiency']*100:.2f}%\n"
            self.output_text.setText(output)
            
            # Update plot
            self.plot_widget.update_plot(
                results['hysteresis_loss'], results['eddy_current_loss'],
                results['copper_loss_primary'], results['copper_loss_secondary']
            )
            
        except ValueError:
            self.output_text.setText("Error: Invalid input. Please enter numeric values.")
        except Exception as e:
            self.output_text.setText(f"Error: {str(e)}\n{traceback.format_exc()}")

class LossPlot(QWidget):
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
        self.ax.set_title("Loss Breakdown")
        self.ax.set_ylabel("Power Loss (W)")
        self.ax.grid(True)
        
    def update_plot(self, hysteresis_loss, eddy_current_loss, copper_loss_primary, copper_loss_secondary):
        try:
            self.ax.clear()
            categories = ['Hysteresis', 'Eddy Current', 'Primary Copper', 'Secondary Copper']
            losses = [hysteresis_loss, eddy_current_loss, copper_loss_primary, copper_loss_secondary]
            self.ax.bar(categories, losses, color=['#FF6666', '#66B2FF', '#FFCC66', '#66FF99'])
            self.ax.set_title("Loss Breakdown")
            self.ax.set_ylabel("Power Loss (W)")
            self.ax.grid(True)
            for i, v in enumerate(losses):
                self.ax.text(i, v + max(losses) * 0.02, f'{v:.2f}', ha='center')
            self.figure.tight_layout()
            self.canvas.draw()
            self.canvas.flush_events()
        except Exception as e:
            print(f"Plot error: {str(e)}\n{traceback.format_exc()}")