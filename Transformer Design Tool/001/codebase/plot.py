from PyQt5.QtWidgets import QWidget, QVBoxLayout, QSizePolicy
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

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
        self.ax.set_xlabel("Time (ms)")
        self.ax.set_ylabel("Flux Density (T)")
        self.ax.grid(True)
        
    def update_plot(self, flux_density):
        self.ax.clear()
        self.ax.plot(flux_density)
        self.ax.set_title("Magnetic Flux Density")
        self.ax.set_xlabel("Time (ms)")
        self.ax.set_ylabel("Flux Density (T)")
        self.ax.grid(True)
        self.canvas.draw()