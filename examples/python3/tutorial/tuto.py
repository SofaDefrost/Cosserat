import sys
from PyQt5.QtWidgets import QApplication, QWidget

app = QApplication(sys.argv)

root = QWidget()

root.resize(250, 250)
root.setWindowTitle("Hello world!")
root.show()

sys.exit(app.exec_())
