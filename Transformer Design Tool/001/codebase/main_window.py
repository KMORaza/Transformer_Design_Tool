from PyQt5.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QLabel, QLineEdit, QPushButton, QTextEdit, QComboBox, QSizePolicy, QMessageBox
from PyQt5.QtGui import QFont
from PyQt5.QtCore import Qt
from transformer import Transformer
from plot import TransformerPlot
from MultiWinding import MultiWindingWindow
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
            output += f"Efficiency: {results['efficiency']*100:.2f}%"
            self.output_text.setText(output)
            
            # Update plot
            self.plot_widget.update_plot(results['flux_densities'][0])
            
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