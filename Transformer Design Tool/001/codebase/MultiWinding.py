from PyQt5.QtWidgets import QDialog, QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QLabel, QLineEdit, QPushButton, QTextEdit, QComboBox, QSizePolicy, QTableWidget, QTableWidgetItem
from PyQt5.QtGui import QFont
from PyQt5.QtCore import Qt
from transformer import Transformer
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import math
import traceback

class MultiWindingWindow(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Multi-Winding Transformer Design")
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
            QLineEdit, QComboBox, QTextEdit, QTableWidget {
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
        
        self.voltage_in = QLineEdit("230")
        self.voltage_in.setFont(QFont("Consolas", 10))
        self.voltage_in.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.frequency = QLineEdit("50")
        self.frequency.setFont(QFont("Consolas", 10))
        self.frequency.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.total_power = QLineEdit("1000")
        self.total_power.setFont(QFont("Consolas", 10))
        self.total_power.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.efficiency = QLineEdit("0.95")
        self.efficiency.setFont(QFont("Consolas", 10))
        self.efficiency.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        
        form_layout.addRow(QLabel("Transformer Type:"), self.transformer_type)
        form_layout.addRow(QLabel("Core Configuration:"), self.core_type)
        form_layout.addRow(QLabel("Primary Voltage (V):"), self.voltage_in)
        form_layout.addRow(QLabel("Frequency (Hz):"), self.frequency)
        form_layout.addRow(QLabel("Total Power (VA):"), self.total_power)
        form_layout.addRow(QLabel("Efficiency (0-1):"), self.efficiency)
        
        # Secondary windings table
        self.winding_table = QTableWidget()
        self.winding_table.setRowCount(2)
        self.winding_table.setColumnCount(2)
        self.winding_table.setHorizontalHeaderLabels(["Secondary Voltage (V)", "Power (VA)"])
        self.winding_table.setFont(QFont("Consolas", 10))
        self.winding_table.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.winding_table.horizontalHeader().setStretchLastSection(True)
        for row in range(2):
            self.winding_table.setItem(row, 0, QTableWidgetItem("12"))
            self.winding_table.setItem(row, 1, QTableWidgetItem("500"))
        
        # Add/Remove winding buttons
        winding_btn_layout = QHBoxLayout()
        self.add_winding_btn = QPushButton("Add Winding")
        self.add_winding_btn.setFont(QFont("Consolas", 10))
        self.add_winding_btn.clicked.connect(self.add_winding)
        self.add_winding_btn.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.remove_winding_btn = QPushButton("Remove Winding")
        self.remove_winding_btn.setFont(QFont("Consolas", 10))
        self.remove_winding_btn.clicked.connect(self.remove_winding)
        self.remove_winding_btn.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        winding_btn_layout.addWidget(self.add_winding_btn)
        winding_btn_layout.addWidget(self.remove_winding_btn)
        
        # Calculate button
        self.calculate_btn = QPushButton("Calculate")
        self.calculate_btn.setFont(QFont("Consolas", 10))
        self.calculate_btn.clicked.connect(self.calculate_design)
        self.calculate_btn.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        
        control_layout.addLayout(form_layout)
        control_layout.addWidget(QLabel("Secondary Windings:"))
        control_layout.addWidget(self.winding_table)
        control_layout.addLayout(winding_btn_layout)
        control_layout.addWidget(self.calculate_btn)
        
        # Plot panel
        self.plot_widget = MultiWindingPlot()
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
        
        output_label = QLabel("Design Results")
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

    def add_winding(self):
        try:
            row_count = self.winding_table.rowCount()
            self.winding_table.insertRow(row_count)
            self.winding_table.setItem(row_count, 0, QTableWidgetItem("12"))
            self.winding_table.setItem(row_count, 1, QTableWidgetItem("500"))
        except Exception as e:
            self.output_text.setText(f"Error adding winding: {str(e)}\n{traceback.format_exc()}")

    def remove_winding(self):
        try:
            if self.winding_table.rowCount() > 1:
                self.winding_table.removeRow(self.winding_table.rowCount() - 1)
        except Exception as e:
            self.output_text.setText(f"Error removing winding: {str(e)}\n{traceback.format_exc()}")

    def calculate_design(self):
        try:
            transformer_type = self.transformer_type.currentText()
            core_type = self.core_type.currentText()
            Vp = float(self.voltage_in.text())
            freq = float(self.frequency.text())
            total_power = float(self.total_power.text())
            eff = float(self.efficiency.text())
            
            # Validate secondary windings data
            Vs_list = []
            power_list = []
            for row in range(self.winding_table.rowCount()):
                voltage_item = self.winding_table.item(row, 0)
                power_item = self.winding_table.item(row, 1)
                if voltage_item is None or power_item is None or not voltage_item.text() or not power_item.text():
                    self.output_text.setText("Error: All secondary winding fields must be filled.")
                    return
                Vs = float(voltage_item.text())
                power = float(power_item.text())
                if Vs <= 0 or power <= 0:
                    self.output_text.setText("Error: Secondary voltages and powers must be positive.")
                    return
                Vs_list.append(Vs)
                power_list.append(power)
            
            # Validate total power
            if sum(power_list) > total_power:
                self.output_text.setText("Error: Sum of secondary powers exceeds total power.")
                return
            
            # Perform calculations
            results = self.transformer.calculate(transformer_type, core_type, Vp, Vs_list, freq, power_list, eff)
            
            # Update output
            output = f"Transformer Type: {transformer_type}\n"
            output += f"Core Type: {core_type}\n"
            output += f"Core Area: {results['core_area']:.4f} m²\n"
            output += f"Primary Turns: {results['primary_turns']:.0f}\n"
            output += f"Primary Wire Gauge: {results['primary_wire_gauge']} AWG\n"
            output += f"Primary Layers: {results['primary_layers']}\n"
            output += f"Interleaving Pattern: {results['interleaving_pattern']}\n"
            output += f"Core Loss: {results['core_loss']:.2f} W\n"
            output += f"Efficiency: {results['efficiency']*100:.2f}%\n"
            output += "\nSecondary Windings:\n"
            for i, (Vs, tr, turns, gauge, layers, mutual_ind) in enumerate(zip(Vs_list, results['turns_ratios'], results['secondary_turns'], results['secondary_wire_gauges'], results['secondary_layers'], results['mutual_inductances'])):
                output += f"Winding {i+1}:\n"
                output += f"  Voltage: {Vs:.2f} V\n"
                output += f"  Turns Ratio: {tr:.2f}\n"
                output += f"  Turns: {turns:.0f}\n"
                output += f"  Wire Gauge: {gauge} AWG\n"
                output += f"  Layers: {layers}\n"
                output += f"  Mutual Inductance (to primary): {mutual_ind:.6f} H\n"
            self.output_text.setText(output)
            
            # Update plot
            self.plot_widget.update_plot(results['flux_densities'], len(Vs_list))
            
        except Exception as e:
            self.output_text.setText(f"Error: {str(e)}\n{traceback.format_exc()}")

class MultiWindingPlot(QWidget):
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
        self.ax.set_title("Magnetic Flux Density (Multiple Windings)")
        self.ax.set_xlabel("Time (ms)")
        self.ax.set_ylabel("Flux Density (T)")
        self.ax.grid(True)
        
    def update_plot(self, flux_densities, num_windings):
        try:
            self.ax.clear()
            for i in range(num_windings):
                self.ax.plot(flux_densities[i], label=f"Winding {i+1}")
            self.ax.set_title("Magnetic Flux Density (Multiple Windings)")
            self.ax.set_xlabel("Time (ms)")
            self.ax.set_ylabel("Flux Density (T)")
            self.ax.grid(True)
            self.ax.legend()
            self.canvas.draw()
        except Exception as e:
            print(f"Plot error: {str(e)}\n{traceback.format_exc()}")