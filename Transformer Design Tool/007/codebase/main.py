import sys
from PyQt5.QtWidgets import QApplication
from main_window import TransformerDesignApp

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = TransformerDesignApp()
    window.show()
    sys.exit(app.exec_())