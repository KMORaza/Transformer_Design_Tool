from PyQt5.QtWidgets import QWidget, QVBoxLayout, QSizePolicy
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import traceback

class TransformerPlot(QWidget):
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
        self.ax.set_title("Magnetic Flux Density")
        self.ax.set_xlabel("Time (s)")
        self.ax.set_ylabel("Flux Density (T)")
        self.ax.grid(True)
        
    def update_plot(self, flux_density, time_points, time_range):
        try:
            self.ax.clear()
            self.ax.plot(time_points, flux_density, label="Secondary")
            self.ax.set_xlim(0, time_range)
            self.ax.set_ylim(-max(abs(min(flux_density)), abs(max(flux_density))) * 1.1,
                            max(abs(min(flux_density)), abs(max(flux_density))) * 1.1)
            self.ax.set_title("Magnetic Flux Density")
            self.ax.set_xlabel("Time (s)")
            self.ax.set_ylabel("Flux Density (T)")
            self.ax.grid(True)
            self.ax.legend()
            self.figure.tight_layout()  # Adjust layout to prevent label cutoff
            self.canvas.draw()
            self.canvas.flush_events()  # Ensure immediate redraw
        except Exception as e:
            print(f"Plot error: {str(e)}\n{traceback.format_exc()}")